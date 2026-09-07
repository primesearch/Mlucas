#!/usr/bin/env bash
# End-to-end tests for the LL Jacobi residue check (help.txt sections [4], [11], [12], [14b]).
#
# Usage: tests/jacobi-check.sh <Mlucas binary> [<Mlucas binary built with -DMLUCAS_FAULT_INJECT>] [<mlucas.cfg>]
#
# Runs small LL tests in a scratch directory and asserts on the .stat file, the savefiles and the
# results.txt JSON line. Exponents: M44497 and M216091 (known primes, whole LL in seconds) and 44501 /
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
	for i in $(seq 1 600); do
		sleep 0.5
		grep -q "Iter# = $it " "$d/p$p.stat" 2>/dev/null && break
		kill -0 $pid 2>/dev/null || break
	done
	kill -INT $pid 2>/dev/null; wait $pid 2>/dev/null; echo $? > "$d/exit"
}
res64_of() { grep -o '"res64":"[0-9A-F]*"' "$1/results.txt" 2>/dev/null | tail -1 | cut -d'"' -f4; }

echo "== Jacobi residue check: end-to-end tests (work dir $WORK)"
echo "   Mlucas: $MLUCAS"; echo "   fault-inject build: ${MLUCAS_FI:-(none - fault-injection tests skipped)}"

# ---------------------------------------------------------------------------------------------
echo "-- T1: clean LL of a Mersenne prime (M44497), check at every checkpoint"
d=$(setup t1 44497 1000); run "$d"
S=$d/p44497.stat
expect_grep  "$S" "M44497 is a known MERSENNE PRIME" "correct verdict"
expect_count "$S" "Jacobi check passed" 45 "a pass logged at every checkpoint plus the final residue"
expect_nogrep "$S" "FAILED" "no failure on a clean run"
expect_grep  "$S" "At iteration 44495, shift = [0-9]*: Jacobi check passed" "final residue checked before the verdict"
[[ -f $d/p44497.J && -f $d/p44497.J1 ]] && ok ".J and .J1 written" || bad ".J/.J1 missing"
[[ -f $d/q44497 ]] && bad "q44497 left behind after completion" || ok "q44497 removed at completion"

# ---------------------------------------------------------------------------------------------
echo "-- T2: clean LL of a composite Mersenne number (M44501): results line carries the Jacobi count"
d=$(setup t2 44501 1000); run "$d"
S=$d/p44501.stat
expect_grep "$S" "M44501 is not prime" "correct verdict"
expect_count "$S" "Jacobi check passed" 45 "a pass logged at every checkpoint plus the final residue"
[[ -f $d/results.txt ]] && ok "results.txt written" || bad "results.txt missing"
if grep -qE '"error-code":"[0-9A-F]{6}0[0-9A-F]"' "$d/results.txt"; then ok "clean Jacobi nibble in error-code"; else bad "Jacobi nibble != 0 on a clean run"; fi
expect_grep "$d/results.txt" '"errors":{"Roundoff":0, "jacobi":0}' "jacobi count 0 in the errors object"
CLEAN_44501=$(res64_of "$d"); [[ ${#CLEAN_44501} -eq 16 ]] && ok "clean Res64 recorded: $CLEAN_44501" || bad "no Res64 in results.txt"

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
	# Damage p and q identically in the residue body (byte 100 lives inside the residue), which
	# breaks their checksum triplet; .J is intact and must be the file the run resumes from.
	for f in p216091 q216091; do printf '\xff' | dd of="$d/$f" bs=1 seek=100 conv=notrunc status=none; done
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
	echo "-- T5: fault injection at iteration 50000 (fixed, deterministic cases)"
	# Each (exponent, iteration, digit, residue shift) gives a deterministic verdict, so the cases are fixed
	# rather than searched for. The runs use -shift 0: the injector perturbs one digit of the residue *as
	# stored*, so with a random shift the perturbed bit position - and hence the verdict - would change from
	# run to run. The injection lands exactly at a checked checkpoint, so the check sees two independent
	# values: the symbol at the corrupted iteration itself (the check at 50000) and the one shared by every
	# later check (from 60000 on). One case per class, all on M216103 (digits classified with -shift 0):
	#   caught      - digit 3: fails at 50000; rollback to the previous checkpoint (40000), correct result,
	#                 1 error counted;
	#   caught_late - digit 13: passes at 50000 (so the corrupt residue is saved to p, q and .J) but fails
	#                 at 60000 and on every retry: the chain must walk p -> .J -> .J1 (40000, clean),
	#                 3 errors counted;
	#   missed      - digit 0: passes both, so no check ever fires, and the run ends with a WRONG residue.
	#                 Asserted on purpose: it is the check's known limit, and this documents it rather than
	#                 hiding it.
	d=$(setup t5clean 216103 10000); run "$d" -- -shift 0
	CLEAN=$(res64_of "$d"); [[ ${#CLEAN} -eq 16 ]] && ok "clean M216103 run Res64 $CLEAN" || bad "clean run produced no results line"

	d=$(setup t5caught 216103 10000 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_ITER=50000 MLUCAS_FAULT_WORD=3 -- -shift 0
	S=$d/p216103.stat
	expect_grep "$S" "FAULT INJECTION: added 1.0 to residue digit 3 at iteration 50000" "caught: injection fired"
	expect_grep "$S" "Jacobi check at iteration 50000 FAILED" "caught: failure at the corrupted checkpoint"
	expect_grep "$S" "Restarting from the current savefile" "caught: rollback to the current savefile announced"
	expect_grep "$S" "Restart file p216103 (iteration 40000) passed the Jacobi check" "caught: rolled back to the previous (clean) checkpoint"
	[[ $(res64_of "$d") == "$CLEAN" ]] && ok "caught: final Res64 matches the clean run" || bad "caught: final Res64 $(res64_of "$d") != clean $CLEAN"
	# error-code bits 4-7 = 7th hex digit; other nibbles may legitimately be nonzero (e.g. a startup roundoff warning)
	if grep -qE '"error-code":"[0-9A-F]{6}1[0-9A-F]"' "$d/results.txt"; then ok "caught: one Jacobi error in the error-code nibble (bits 4-7)"; else bad "caught: Jacobi nibble != 1: $(grep -o '"error-code":"[0-9A-F]*"' "$d/results.txt")"; fi
	expect_grep "$d/results.txt" '"jacobi":1' "caught: jacobi count 1 in the errors object"

	d=$(setup t5late 216103 10000 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_ITER=50000 MLUCAS_FAULT_WORD=13 -- -shift 0
	S=$d/p216103.stat
	expect_grep "$S" "FAULT INJECTION: added 1.0 to residue digit 13 at iteration 50000" "caught_late: injection fired"
	expect_nogrep "$S" "Jacobi check at iteration 50000 FAILED" "caught_late: passes the check at its own checkpoint"
	expect_count "$S" "Jacobi check at iteration 60000 FAILED" 3 "caught_late: the same failure recurs on each retry from a poisoned file"
	expect_grep "$S" "Restart file p216103 (iteration 50000) passed the Jacobi check" "caught_late: p (poisoned) passes its on-read check, as the symbol invariance predicts"
	expect_grep "$S" "Restart file p216103.J (iteration 50000) passed the Jacobi check" "caught_late: .J (poisoned) likewise"
	expect_grep "$S" "Restart file p216103.J1 (iteration 40000) passed the Jacobi check" "caught_late: chain reaches .J1, the clean checkpoint"
	[[ $(res64_of "$d") == "$CLEAN" ]] && ok "caught_late: final Res64 matches the clean run" || bad "caught_late: final Res64 $(res64_of "$d") != clean $CLEAN"
	if grep -qE '"error-code":"[0-9A-F]{6}3[0-9A-F]"' "$d/results.txt"; then ok "caught_late: three Jacobi errors in the error-code nibble"; else bad "caught_late: Jacobi nibble != 3: $(grep -o '"error-code":"[0-9A-F]*"' "$d/results.txt")"; fi

	d=$(setup t5missed 216103 10000 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_ITER=50000 MLUCAS_FAULT_WORD=0 -- -shift 0
	S=$d/p216103.stat
	expect_grep "$S" "FAULT INJECTION: added 1.0 to residue digit 0 at iteration 50000" "missed: injection fired"
	expect_nogrep "$S" "FAILED" "missed: no check ever fires (symbol invariant under the recurrence)"
	[[ $(res64_of "$d") != "$CLEAN" ]] && ok "missed: run ends with a WRONG residue $(res64_of "$d") - this is the documented limit of the check" || bad "missed: final residue equals the clean one?!"
	if grep -qE '"error-code":"[0-9A-F]{6}0[0-9A-F]"' "$d/results.txt"; then ok "missed: Jacobi nibble 0, as expected - nothing was detected"; else bad "missed: Jacobi nibble != 0"; fi
fi

# ---------------------------------------------------------------------------------------------
if [[ -n $MLUCAS_FI ]]; then
	echo "-- T6: a reproducible fault: the rollback chain is walked to the end and the run aborts rather than loop"
	# MLUCAS_FAULT_REPEAT re-fires the injection at every visit of iteration 50000, so every retry fails there
	# again: current savefile -> .J -> .J1 -> scratch -> fifth failure aborts with a hardware warning.
	d=$(setup t6 216103 10000 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_ITER=50000 MLUCAS_FAULT_WORD=3 MLUCAS_FAULT_REPEAT=1 -- -shift 0
	S=$d/p216103.stat
	expect_count "$S" "FAULT INJECTION: added 1.0 to residue digit 3 at iteration 50000" 5 "injection re-fired on every retry"
	expect_count "$S" "Jacobi check at iteration 50000 FAILED" 4 "four failures each followed by a rollback"
	expect_grep "$S" "Restarting from the current savefile" "1st rollback: current savefile"
	expect_grep "$S" "Restarting from the last Jacobi-passed savefile" "2nd rollback: .J"
	expect_grep "$S" "Restarting from the previous Jacobi-passed savefile" "3rd rollback: .J1"
	expect_grep "$S" "Restarting from scratch" "4th rollback: scratch"
	expect_grep "$S" "failed 5 times in a row, the last after restarting from scratch" "5th failure aborts with the hardware warning"
	[[ $(cat "$d/exit") != 0 ]] && ok "run stopped (exit $(cat "$d/exit"))" || bad "run did not stop"
	[[ -f $d/p216103 && -f $d/p216103.J ]] && ok "savefiles left in place" || bad "savefiles missing after the abort"

	echo "-- T7: PRP run whose run flag is cleared between intervals (a signal during the checkpoint write) must stop, not run on frozen"
	d=$WORK/t7; rm -rf "$d"; mkdir -p "$d"; ln -s "$MLUCAS_FI" "$d/Mlucas"; [[ -n $CFG && -f $CFG ]] && cp "$CFG" "$d/mlucas.cfg"
	printf 'PRP=1,2,216091,-1,75,0,3,1\n' > "$d/worktodo.txt"; printf 'CheckInterval = 10000\n' > "$d/mlucas.ini"
	run "$d" MLUCAS_FAULT_STOP_AT=50000; S=$d/p216091.stat
	expect_grep "$S" "FAULT INJECTION: run flag cleared between intervals at iteration 50000" "flag cleared between intervals"
	expect_grep "$S" "Received SIGINT signal: writing savefile at Iter = 50000 and exiting" "treated as an interrupt at the last completed iteration"
	expect_nogrep "$S" "Iter# = 60000" "no further intervals were 'completed' after the stop"
	expect_nogrep "$S" "MaxErr = 0.000000000" "no frozen-residue checkpoints"
	[[ $(cat "$d/exit") == 0 ]] && ok "clean exit after the savefile write" || bad "exit $(cat "$d/exit")"
	run "$d"
	expect_grep "$S" "Restarting M216091 at iteration = 50000" "resumed at the interrupted iteration"
	expect_grep "$S" "Gerbicz check passed" "the run's Gerbicz check passed after the resume"
	expect_grep "$S" "M216091 is a known MERSENNE PRIME" "correct PRP verdict"
else
	skip "T6/T7 need the -DMLUCAS_FAULT_INJECT build"
fi

# ---------------------------------------------------------------------------------------------
echo "== $PASS passed, $FAIL failed, $SKIP skipped (work dir $WORK)"
(( FAIL == 0 ))
