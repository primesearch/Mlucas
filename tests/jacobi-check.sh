#!/usr/bin/env bash
# End-to-end tests for the LL Jacobi residue check (help.txt sections [4], [11], [12], [14b]).
#
# Usage: tests/jacobi-check.sh <Mlucas binary> [<Mlucas binary built with -DMLUCAS_FAULT_INJECT>] [<mlucas.cfg>]
#
# Runs small LL tests in a scratch directory and asserts on the .stat file, the savefiles and the
# results.txt JSON line. Exponents: M86243 and M216091 (known primes, whole LL in seconds) and 86249 /
# 216103 (prime exponents with composite Mersenne numbers, so a results line is emitted). The
# fault-injection tests need the second binary; they are skipped with a notice if it is not given.
#
# What is asserted, and why, is spelled out next to each test. Every expectation is one the check's
# design makes: a check that passes at every checkpoint of a clean run; a corrupted savefile that is
# refused on read; an injected fault caught at the next check with a rollback to the previous checkpoint
# and a correct final result; an injected fault that passes its own checkpoint but is caught at the next
# and recovered through the .J1 generation; and - deliberately - an injected fault that the check MISSES
# and a run that then ends wrong. The last one documents the check's known blind spot rather than hiding
# it: the Jacobi symbol is invariant under the LL recurrence from one iteration after a corruption on,
# so a corruption that passes the first check after it is never caught by a later one.
set -u -o pipefail

MLUCAS=${1:?usage: $0 <Mlucas> [<Mlucas-fault-inject>] [<mlucas.cfg>]}
MLUCAS_FI=${2:-}
CFG=${3:-}
MLUCAS=$(readlink -f "$MLUCAS"); [[ -n $MLUCAS_FI ]] && MLUCAS_FI=$(readlink -f "$MLUCAS_FI")
[[ -n $CFG ]] && CFG=$(readlink -f "$CFG")
CPU=${JACOBI_TEST_CPU:-0:3}
WORK=${JACOBI_TEST_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/jacobi-check.XXXXXX")}
mkdir -p "$WORK"
PASS=0; FAIL=0; SKIP=0

ok()   { echo "  ok   - $*"; PASS=$((PASS+1)); }
bad()  { echo "  FAIL - $*"; FAIL=$((FAIL+1)); }
skip() { echo "  skip - $*"; SKIP=$((SKIP+1)); }
expect_grep()   { local f=$1 pat=$2 what=$3; if grep -q -- "$pat" "$f"; then ok "$what"; else bad "$what (pattern not found: $pat)"; fi; }
expect_nogrep() { local f=$1 pat=$2 what=$3; if grep -q -- "$pat" "$f"; then bad "$what (unexpected: $(grep -m1 -- "$pat" "$f"))"; else ok "$what"; fi; }
expect_count()  { local f=$1 pat=$2 min=$3 what=$4; local n; n=$(grep -c -- "$pat" "$f"); if (( n >= min )); then ok "$what ($n)"; else bad "$what: $n < $min"; fi; }

# Fresh run directory with the binary, cfg (if any) and the mlucas.ini options the tests rely on.
setup() {	# setup <dir> <exponent> <CheckInterval> [<binary>]
	local d=$WORK/$1 p=$2 ci=$3 bin=${4:-$MLUCAS}
	rm -rf "$d"; mkdir -p "$d"
	ln -s "$bin" "$d/Mlucas"
	[[ -n $CFG && -f $CFG ]] && cp "$CFG" "$d/mlucas.cfg"
	printf 'Test=%s\n' "$p" > "$d/worktodo.txt"
	printf 'CheckInterval = %s\nJacobiCheckHours = 0\n' "$ci" > "$d/mlucas.ini"
	echo "$d"
}
run() {	# run <dir> [env...] [-- <extra Mlucas args>] : run Mlucas to completion (or until it exits) in <dir>, log to run.log
	local d=$1; shift; local envs=() extra=()
	while (( $# )); do if [[ $1 == -- ]]; then shift; extra=("$@"); break; fi; envs+=("$1"); shift; done
	( cd "$d" && env "${envs[@]}" timeout 900 ./Mlucas -cpu "$CPU" "${extra[@]}" > run.log 2>&1 ); echo $? > "$d/exit"
}
# Run in the background and stop it with SIGINT once the .stat file shows the given iteration.
run_until_iter_then_interrupt() {	# <dir> <exponent> <iteration>
	local d=$1 p=$2 it=$3 pid i
	( cd "$d" && exec ./Mlucas -cpu "$CPU" > run.log 2>&1 ) & pid=$!
	for ((i = 0; i < 600; ++i)); do
		sleep 0.5
		grep -q "Iter# = $it " "$d/p$p.stat" 2>/dev/null && break
		kill -0 "$pid" 2>/dev/null || break
	done
	kill -INT "$pid" 2>/dev/null; wait "$pid" 2>/dev/null; echo $? > "$d/exit"
}
res64_of() { grep -o '"res64":"[0-9A-F]*"' "$1/results.txt" 2>/dev/null | tail -1 | cut -d'"' -f4; }

echo "== Jacobi residue check: end-to-end tests (work dir $WORK)"
echo "   Mlucas: $MLUCAS"; echo "   fault-inject build: ${MLUCAS_FI:-(none - fault-injection tests skipped)}"

# ---------------------------------------------------------------------------------------------
# T1 and T2 use exponents whose *default* FFT length is 4K, not 2K. An AVX-512 build cannot run 2K at
# all: the teensy-FFT guard in mers_mod_square.c wants complex-length/radix_final >= 16*RE_IM_STRIDE,
# which is 128 there against 64 for AVX/AVX2, and no 2K radix set reaches it. On such a host every 2K
# radix set is rejected, no mlucas.cfg entry for 2K can exist, the remedial timing self-test at 2K
# fails the same way, and the run aborts before doing any work. Measured on a Zen 4: 'Mlucas -s tt'
# passes 4 of 26 cases there and writes cfg entries for 4K and 8K only.
#
# Forcing the length with -fft does not work around it: for a Mersenne production run the 9/8 rule in
# ernstMain() discards any forced length more than one size above the default, so '-fft 4' on a 2K
# exponent is silently reverted to 2K (8*4 > 9*2). The exponent has to be one whose default is already
# usable. M86243 is the Mersenne prime in the 4K band (66742 < p <= 88438) and 86249 is a prime whose
# Mersenne number is composite, so T2 still gets its results line.
echo "-- T1: clean LL of a Mersenne prime (M86243), check at every checkpoint"
d=$(setup t1 86243 1000); run "$d"
S=$d/p86243.stat
expect_grep  "$S" "M86243 is a known MERSENNE PRIME" "correct verdict"
expect_count "$S" "Jacobi check passed" 87 "a pass logged at every checkpoint plus the final residue"
expect_nogrep "$S" "FAILED" "no failure on a clean run"
expect_grep  "$S" "At iteration 86241, shift = [0-9]*: Jacobi check passed" "final residue checked before the verdict"
if [[ -f $d/p86243.J && -f $d/p86243.J1 ]]; then ok ".J and .J1 written"; else bad ".J/.J1 missing"; fi
if [[ -f $d/q86243 ]]; then bad "q86243 left behind after completion"; else ok "q86243 removed at completion"; fi

# ---------------------------------------------------------------------------------------------
echo "-- T2: clean LL of a composite Mersenne number (M86249): results line carries the Jacobi count"
d=$(setup t2 86249 1000); run "$d"
S=$d/p86249.stat
expect_grep "$S" "M86249 is not prime" "correct verdict"
expect_count "$S" "Jacobi check passed" 87 "a pass logged at every checkpoint plus the final residue"
if [[ -f $d/results.txt ]]; then ok "results.txt written"; else bad "results.txt missing"; fi
if grep -qE '"error-code":"[0-9A-F]{6}0[0-9A-F]"' "$d/results.txt"; then ok "clean Jacobi nibble in error-code"; else bad "Jacobi nibble != 0 on a clean run"; fi
expect_grep "$d/results.txt" '"errors":{"Roundoff":0, "jacobi":0}' "jacobi count 0 in the errors object"
CLEAN_86249=$(res64_of "$d"); if [[ ${#CLEAN_86249} -eq 16 ]]; then ok "clean Res64 recorded: $CLEAN_86249"; else bad "no Res64 in results.txt"; fi

# ---------------------------------------------------------------------------------------------
echo "-- T3: interrupt and resume (M216091): the loaded residue is Jacobi-checked before use"
d=$(setup t3 216091 10000)
run_until_iter_then_interrupt "$d" 216091 30000
S=$d/p216091.stat
expect_grep "$S" "Received SIGINT signal: writing savefile" "interrupt wrote a savefile"
if grep -q "Received SIGINT signal: writing savefile at Iter = \([0-9]*\)" "$S"; then
	it=$(grep -o "writing savefile at Iter = [0-9]*" "$S" | grep -o "[0-9]*$")
	# The interrupt-driven checkpoint must not run the check (shutdown must stay fast). Order-aware: when the signal
	# lands between intervals the interrupt is reported at the iteration of the regular checkpoint just written, whose
	# own check legitimately ran *before* the signal - so look only at what follows the interrupt message:
	if sed -n '/Received SIGINT signal/,$p' "$S" | grep -q "Jacobi check passed"; then bad "a Jacobi check ran on the interrupt checkpoint"; else ok "no Jacobi check on the interrupt checkpoint"; fi
	run "$d"
	# ...and the resume must check the residue it loads:
	expect_grep "$S" "Restart file p216091 (iteration $it) passed the Jacobi check" "restart-read Jacobi check of the saved residue"
	expect_grep "$S" "Restarting M216091 at iteration = $it" "resumed at the interrupted iteration"
	expect_grep "$S" "M216091 is a known MERSENNE PRIME" "correct verdict after resume"
else
	bad "could not interrupt the run in time"
fi

# ---------------------------------------------------------------------------------------------
echo "-- T4: corrupted savefile is refused on read and the chain falls through (M216091)"
d=$(setup t4 216091 10000)
run_until_iter_then_interrupt "$d" 216091 30000
S=$d/p216091.stat
if [[ -f $d/p216091 && -f $d/q216091 && -f $d/p216091.J ]]; then
	# Damage p and q identically in the residue body, which breaks their checksum triplet; .J is
	# intact and must be the file the run resumes from. Seek to the middle of the file rather than a
	# fixed byte: the header length is not the same for every build, and a fixed small offset can land
	# in it, where the damage does not break the residue checksum and the file is accepted.
	for f in p216091 q216091; do
		sz=$(wc -c < "$d/$f")
		printf '\xff' | dd of="$d/$f" bs=1 seek=$((sz / 2)) conv=notrunc status=none
	done
	run "$d"
	expect_grep "$S" "read_ppm1_savefiles Failed on savefile p216091" "damaged primary rejected"
	expect_grep "$S" "read_ppm1_savefiles Failed on savefile q216091" "damaged secondary rejected"
	expect_grep "$S" "Restart file p216091.J (iteration [0-9]*) passed the Jacobi check" "resumed from the last Jacobi-passed checkpoint"
	expect_grep "$S" "M216091 is a known MERSENNE PRIME" "correct verdict after falling through the chain"
else
	bad "expected p/q/.J savefiles after the interrupt"
fi

# ---------------------------------------------------------------------------------------------
if [[ -z $MLUCAS_FI ]]; then
	skip "T5: fault-injection tests need a -DMLUCAS_FAULT_INJECT build as the 2nd argument"
else
	echo "-- T5: fault injection at iteration 50000 (verdict classified at run time)"
	# Which digit produces which verdict is not a fixed property. The injector perturbs one digit of the
	# residue as stored, and where that digit lands depends on the FFT length and radix set the host
	# chose - which vary with the build mode and with whatever the timing self-test wrote to mlucas.cfg.
	# Pinning a digit per verdict was reproducible only on the machine the classification was done on:
	# digit 13 gives "caught at its own checkpoint" on one build and "caught at the next one" on another,
	# from a byte-identical residue at the injected iteration, and the fixed mapping failed on AVX-512,
	# AVX2 and ASIMD alike.
	#
	# So classify at run time and assert what is true of every digit:
	#   - if any check fires, the rollback chain must recover and the final residue must match the clean
	#     run, and the error-code Jacobi nibble must equal the number of failures;
	#   - if none fires, the run must end with a WRONG residue. That is the check's documented blind
	#     spot - the Jacobi symbol is invariant under the LL recurrence from one iteration after a
	#     corruption on - and asserting it keeps the limit visible rather than hidden.
	# One global assertion keeps the suite honest: at least one digit must be caught, so a check that
	# never fires cannot pass this test.
	d=$(setup t5clean 216103 10000); run "$d" -- -shift 0
	CLEAN=$(res64_of "$d"); if [[ ${#CLEAN} -eq 16 ]]; then ok "clean M216103 run Res64 $CLEAN"; else bad "clean run produced no results line"; fi

	ncaught=0; CAUGHT_AT_OWN=''
	for w in 0 3 13 21; do
		d=$(setup "t5w$w" 216103 10000 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_ITER=50000 "MLUCAS_FAULT_WORD=$w" -- -shift 0
		S=$d/p216103.stat
		if [[ ! -f $S ]]; then bad "digit $w: no .stat file produced"; continue; fi
		expect_grep "$S" "FAULT INJECTION: added 1.0 to residue digit $w at iteration 50000" "digit $w: injection fired"
		nfail=$(grep -c "Jacobi check at iteration [0-9]* FAILED" "$S")
		got=$(res64_of "$d")
		if (( nfail > 0 )); then
			ncaught=$((ncaught + 1))
			# A digit caught at its own checkpoint is the one T6 needs, so remember the first such.
			if [[ -z $CAUGHT_AT_OWN ]] && grep -q "Jacobi check at iteration 50000 FAILED" "$S"; then CAUGHT_AT_OWN=$w; fi
			if [[ $got == "$CLEAN" ]]; then ok "digit $w: caught ($nfail check failure(s)), chain recovered the clean residue"; else bad "digit $w: caught but final Res64 $got != clean $CLEAN"; fi
			expect_grep "$S" "Restarting from" "digit $w: a rollback was announced"
			nib=$(printf '%X' $(( nfail > 15 ? 15 : nfail )))
			if grep -qE "\"error-code\":\"[0-9A-F]{6}${nib}[0-9A-F]\"" "$d/results.txt"; then ok "digit $w: Jacobi nibble $nib matches the $nfail failure(s)"; else bad "digit $w: Jacobi nibble != $nib: $(grep -o '"error-code":"[0-9A-F]*"' "$d/results.txt")"; fi
		else
			if [[ $got != "$CLEAN" ]]; then ok "digit $w: not caught - run ends with a wrong residue $got, the check's documented blind spot"; else bad "digit $w: not caught, yet the final residue is clean?!"; fi
			if grep -qE '"error-code":"[0-9A-F]{6}0[0-9A-F]"' "$d/results.txt"; then ok "digit $w: Jacobi nibble 0, nothing was detected"; else bad "digit $w: nibble nonzero with no failure: $(grep -o '"error-code":"[0-9A-F]*"' "$d/results.txt")"; fi
		fi
	done
	if (( ncaught > 0 )); then ok "the check caught $ncaught of the 4 injected faults"; else bad "no injected fault was caught at all - the check is inert"; fi
fi

# ---------------------------------------------------------------------------------------------
if [[ -n $MLUCAS_FI ]]; then
	echo "-- T6: a reproducible fault: the rollback chain is walked to the end and the run aborts rather than loop"
	# MLUCAS_FAULT_REPEAT re-fires the injection at every visit of iteration 50000, so every retry fails there
	# again: current savefile -> .J -> .J1 -> scratch -> fifth failure aborts with a hardware warning.
	# The digit must be one that fails the check at its own checkpoint, which T5 established above for
	# this host; hard-coding one made this test pass only where that classification happened to hold.
	if [[ -z ${CAUGHT_AT_OWN:-} ]]; then
		skip "T6: no digit was caught at its own checkpoint on this build, so the retry chain cannot be driven"
	else
	d=$(setup t6 216103 10000 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_ITER=50000 "MLUCAS_FAULT_WORD=$CAUGHT_AT_OWN" MLUCAS_FAULT_REPEAT=1 -- -shift 0
	S=$d/p216103.stat
	expect_count "$S" "FAULT INJECTION: added 1.0 to residue digit $CAUGHT_AT_OWN at iteration 50000" 5 "injection re-fired on every retry"
	expect_count "$S" "Jacobi check at iteration 50000 FAILED" 4 "four failures each followed by a rollback"
	expect_grep "$S" "Restarting from the current savefile" "1st rollback: current savefile"
	expect_grep "$S" "Restarting from the last Jacobi-passed savefile" "2nd rollback: .J"
	expect_grep "$S" "Restarting from the previous Jacobi-passed savefile" "3rd rollback: .J1"
	expect_grep "$S" "Restarting from scratch" "4th rollback: scratch"
	expect_grep "$S" "failed 5 times in a row, the last after restarting from scratch" "5th failure aborts with the hardware warning"
	if [[ $(cat "$d/exit") != 0 ]]; then ok "run stopped (exit $(cat "$d/exit"))"; else bad "run did not stop"; fi
	if [[ -f $d/p216103 && -f $d/p216103.J ]]; then ok "savefiles left in place"; else bad "savefiles missing after the abort"; fi
	fi

	echo "-- T7: PRP run whose run flag is cleared between intervals (a signal during the checkpoint write) must stop, not run on frozen"
	d=$WORK/t7; rm -rf "$d"; mkdir -p "$d"; ln -s "$MLUCAS_FI" "$d/Mlucas"; [[ -n $CFG && -f $CFG ]] && cp "$CFG" "$d/mlucas.cfg"
	printf 'PRP=1,2,216091,-1,75,0,3,1\n' > "$d/worktodo.txt"; printf 'CheckInterval = 10000\n' > "$d/mlucas.ini"
	run "$d" MLUCAS_FAULT_STOP_AT=50000; S=$d/p216091.stat
	expect_grep "$S" "FAULT INJECTION: run flag cleared between intervals at iteration 50000" "flag cleared between intervals"
	expect_grep "$S" "Received SIGINT signal: writing savefile at Iter = 50000 and exiting" "treated as an interrupt at the last completed iteration"
	expect_nogrep "$S" "Iter# = 60000" "no further intervals were 'completed' after the stop"
	expect_nogrep "$S" "MaxErr = 0.000000000" "no frozen-residue checkpoints"
	if [[ $(cat "$d/exit") == 0 ]]; then ok "clean exit after the savefile write"; else bad "exit $(cat "$d/exit")"; fi
	run "$d"
	expect_grep "$S" "Restarting M216091 at iteration = 50000" "resumed at the interrupted iteration"
	expect_grep "$S" "Gerbicz check passed" "the run's Gerbicz check passed after the resume"
	expect_grep "$S" "M216091 is a known MERSENNE PRIME" "correct PRP verdict"
else
	skip "T5/T6/T7 need the -DMLUCAS_FAULT_INJECT build"
fi

# ---------------------------------------------------------------------------------------------
echo "== $PASS passed, $FAIL failed, $SKIP skipped (work dir $WORK)"
(( FAIL == 0 ))
