#!/usr/bin/env bash
# End-to-end tests for the p-1 stage 1 Gerbicz check and the p-1 Jacobi integrity checks (help.txt section [9]).
#
# Usage: tests/jacobi-pm1.sh <Mlucas binary> [<Mlucas binary built with -DMLUCAS_FAULT_INJECT>] [<mlucas.cfg>]
#
# Exponent 44497 with B1 = B2 = 1500000: stage 1 is 2165373 iterations (~30 s at a 2K FFT), so with the default
# check interval of 10^6 there are Gerbicz checks at 1000000 and 2000000 plus the end-of-run check at 2165000.
# Stage 2 is not run (B2 = B1).
#
# What is asserted, and why: a clean run passes every check and the final residue passes its Jacobi check; a
# run interrupted mid-epoch resumes with the stored check-product (the epoch continues, and the next check
# still passes); a savefile whose residue is damaged is refused by the Jacobi read check and the run falls back
# to the secondary; and an injected fault at a checkpoint is caught by the next Gerbicz check with certainty
# (unlike the LL Jacobi check, this one is exact), the run rolls back to the last .G checkpoint and finishes
# with the same stage 1 residue as the clean run.
set -u -o pipefail

MLUCAS=${1:?usage: $0 <Mlucas> [<Mlucas-fault-inject>] [<mlucas.cfg>]}
MLUCAS_FI=${2:-}
CFG=${3:-}
MLUCAS=$(readlink -f "$MLUCAS"); [[ -n $MLUCAS_FI ]] && MLUCAS_FI=$(readlink -f "$MLUCAS_FI")
[[ -n $CFG ]] && CFG=$(readlink -f "$CFG")
CPU=${JACOBI_TEST_CPU:-0:3}
WORK=${JACOBI_TEST_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/jacobi-pm1.XXXXXX")}
mkdir -p "$WORK"
PASS=0; FAIL=0; SKIP=0
P=44497; B1=1500000; LAST=2165373; LASTCHK=2165000

ok()   { echo "  ok   - $*"; PASS=$((PASS+1)); }
bad()  { echo "  FAIL - $*"; FAIL=$((FAIL+1)); }
skip() { echo "  skip - $*"; SKIP=$((SKIP+1)); }
expect_grep()   { local f=$1 pat=$2 what=$3; if grep -q -- "$pat" "$f"; then ok "$what"; else bad "$what (pattern not found: $pat)"; fi; }
expect_nogrep() { local f=$1 pat=$2 what=$3; if grep -q -- "$pat" "$f"; then bad "$what (unexpected: $(grep -m1 -- "$pat" "$f"))"; else ok "$what"; fi; }
expect_count()  { local f=$1 pat=$2 min=$3 what=$4; local n; n=$(grep -c -- "$pat" "$f"); if (( n >= min )); then ok "$what ($n)"; else bad "$what: $n < $min"; fi; }

setup() {	# setup <dir> [<binary>]
	local d=$WORK/$1 bin=${2:-$MLUCAS}
	rm -rf "$d"; mkdir -p "$d"
	ln -s "$bin" "$d/Mlucas"
	[[ -n $CFG && -f $CFG ]] && cp "$CFG" "$d/mlucas.cfg"
	printf 'Pminus1=1,2,%s,-1,%s,%s\n' "$P" "$B1" "$B1" > "$d/worktodo.txt"
	printf 'CheckInterval = 10000\n' > "$d/mlucas.ini"
	echo "$d"
}
run() {	# run <dir> [env...]
	local d=$1; shift
	( cd "$d" && env "$@" timeout 1800 ./Mlucas -cpu "$CPU" > run.log 2>&1 ); echo $? > "$d/exit"
}
run_until_iter_then_interrupt() {	# <dir> <iteration>
	local d=$1 it=$2 pid i
	( cd "$d" && exec ./Mlucas -cpu "$CPU" > run.log 2>&1 ) & pid=$!
	for i in $(seq 1 1200); do
		sleep 0.5
		grep -q "S1 bit = $it " "$d/p$P.stat" 2>/dev/null && break
		kill -0 $pid 2>/dev/null || break
	done
	kill -INT $pid 2>/dev/null; wait $pid 2>/dev/null; echo $? > "$d/exit"
}
final_res64() { grep "S1 bit = $LAST " "$1/p$P.stat" | grep -o "Res64: [0-9A-F]*" | tail -1 | cut -d' ' -f2; }

echo "== p-1 stage 1 Gerbicz / Jacobi checks: end-to-end tests (work dir $WORK)"

# ---------------------------------------------------------------------------------------------
echo "-- P1: clean p-1 stage 1 (M$P, B1 = $B1)"
d=$(setup p1); run "$d"; S=$d/p$P.stat
expect_grep  "$S" "At iteration 1000000, shift = 0: Gerbicz check passed" "check at 10^6 passed"
expect_grep  "$S" "At iteration 2000000, shift = 0: Gerbicz check passed" "check at 2*10^6 passed"
expect_grep  "$S" "At iteration $LASTCHK, shift = 0: Gerbicz check passed" "end-of-run check at the last block boundary passed"
expect_nogrep "$S" "Gerbicz check iteration" "no failures"
expect_grep  "$S" "Stage 1 final residue passed the Jacobi check" "final-residue Jacobi (conversion-path) check passed"
expect_grep  "$S" "GCD" "stage 1 GCD ran"
[[ -f $d/p$P.G ]] && ok ".G (last-good-check) savefile written" || bad ".G missing"
expect_grep "$d/results.txt" '"errors":{"Roundoff":[0-9]*, "gerbicz":0, "jacobi":0}' "results line carries zero Gerbicz and Jacobi counts"
CLEAN=$(final_res64 "$d"); [[ ${#CLEAN} -eq 16 ]] && ok "clean final stage 1 Res64 $CLEAN" || bad "no final stage 1 Res64 in the .stat file"

# ---------------------------------------------------------------------------------------------
echo "-- P2: interrupt mid-epoch and resume: the stored check-product carries the epoch across the restart"
d=$(setup p2); run_until_iter_then_interrupt "$d" 1300000; S=$d/p$P.stat
if grep -q "Received SIGINT signal: writing savefile at Iter = " "$S"; then
	run "$d"
	expect_grep  "$S" "Restart file p$P (stage 1 iteration [0-9]*) passed the Jacobi check" "loaded residue passed the Jacobi read check"
	expect_nogrep "$S" "carries no Gerbicz check-product" "the epoch continued from the stored product (no new epoch)"
	expect_grep  "$S" "At iteration 2000000, shift = 0: Gerbicz check passed" "the check spanning the restart passed"
	[[ $(final_res64 "$d") == "$CLEAN" ]] && ok "final stage 1 Res64 matches the clean run" || bad "final Res64 $(final_res64 "$d") != clean $CLEAN"
else
	bad "could not interrupt the run in time"
fi

# ---------------------------------------------------------------------------------------------
echo "-- P3: damaged primary savefile is refused by the Jacobi read check; the secondary is used"
d=$(setup p3); run_until_iter_then_interrupt "$d" 1300000; S=$d/p$P.stat
if [[ -f $d/p$P && -f $d/q$P ]]; then
	# Flip a residue byte in p only and re-derive nothing: the S-H triplet catches this first, so also patch the
	# triplet-consistent case would need the corruptor tool; here the point is the fall-through to q.
	printf '\xff' | dd of="$d/p$P" bs=1 seek=100 conv=notrunc status=none
	run "$d"
	expect_grep "$S" "read_ppm1_savefiles Failed on savefile p$P" "damaged primary rejected"
	expect_grep "$S" "Restart file q$P (stage 1 iteration [0-9]*) passed the Jacobi check" "secondary passed the Jacobi read check"
	[[ $(final_res64 "$d") == "$CLEAN" ]] && ok "final stage 1 Res64 matches the clean run" || bad "final Res64 $(final_res64 "$d") != clean $CLEAN"
else
	bad "expected p/q savefiles after the interrupt"
fi

# ---------------------------------------------------------------------------------------------
if [[ -z $MLUCAS_FI ]]; then
	skip "P4: fault-injection test needs a -DMLUCAS_FAULT_INJECT build as the 2nd argument"
else
	echo "-- P4: fault injected at iteration 1500000: the Gerbicz check at 2000000 must catch it (exact, not a coin flip)"
	d=$(setup p4 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_ITER=1500000 MLUCAS_FAULT_WORD=0; S=$d/p$P.stat
	expect_grep "$S" "FAULT INJECTION: added 1.0 to residue digit 0 at iteration 1500000" "injection fired"
	expect_grep "$S" "Gerbicz check iteration 2000000 failed! Restarting from last-good-Gerbicz-check data" "caught at the next check"
	expect_grep "$S" "Restart file p$P.G (stage 1 iteration 1000000) passed the Jacobi check" "rolled back to the .G checkpoint at 10^6"
	expect_grep "$S" "At iteration 2000000, shift = 0: Gerbicz check passed" "the retried block passed"
	[[ $(final_res64 "$d") == "$CLEAN" ]] && ok "final stage 1 Res64 matches the clean run" || bad "final Res64 $(final_res64 "$d") != clean $CLEAN"
	expect_grep "$d/results.txt" '"gerbicz":1' "one Gerbicz error counted in the results line"
	# bits 20-23 = 3rd hex digit of the 8-digit code (digits are bits 31-28, 27-24, 23-20, ...)
	if grep -qE '"error-code":"[0-9A-F]{2}1[0-9A-F]{5}"' "$d/results.txt"; then ok "Gerbicz nibble (bits 20-23) = 1 in the error-code"; else bad "Gerbicz nibble != 1: $(grep -o '"error-code":"[0-9A-F]*"' "$d/results.txt")"; fi
fi

# ---------------------------------------------------------------------------------------------
# Stage 2: B1 = 20000 (fast stage 1), B2 = 2000000. The exact input checks - ladder recomputation, static-table
# checksums - and the accumulator sanity tests run at every stage 2 checkpoint here (JacobiCheckHours = 0).
setup2() {	# setup2 <dir> [<binary>]
	local d=$WORK/$1 bin=${2:-$MLUCAS}
	rm -rf "$d"; mkdir -p "$d"
	ln -s "$bin" "$d/Mlucas"
	[[ -n $CFG && -f $CFG ]] && cp "$CFG" "$d/mlucas.cfg"
	printf 'Pminus1=1,2,%s,-1,20000,2000000\n' "$P" > "$d/worktodo.txt"
	printf 'CheckInterval = 10000\nJacobiCheckHours = 0\n' > "$d/mlucas.ini"
	echo "$d"
}
s2_final() { grep -o '"B2":[0-9]*, "error-code":"[0-9A-F]*"' "$1/results.txt" 2>/dev/null | tail -1; grep "S2 at q" "$1/p$P.stat" | tail -1 | grep -o "Res64: [0-9A-F]*"; }

echo "-- P5: clean stage 2 (B1 = 20000, B2 = 2000000): ladder, table and accumulator checks pass at every checkpoint"
d=$(setup2 p5); run "$d"; S=$d/p$P.stat
expect_count "$S" "Stage 2 ladder check passed at q = " 5 "ladder recomputation passed at the stage 2 checkpoints"
expect_nogrep "$S" "ERROR: M$P stage 2" "no stage 2 check failure"
expect_grep "$d/results.txt" '"B2":2000000' "stage 2 completed and reported"
CLEAN2=$(s2_final "$d" | tail -1); [[ -n $CLEAN2 ]] && ok "clean stage 2 final accumulator $CLEAN2" || bad "no stage 2 Res64 in the .stat file"

if [[ -n $MLUCAS_FI ]]; then
	echo "-- P6: stage 2 ladder corrupted at q ~ 1000000: caught exactly by the recomputation; restart finishes clean"
	d=$(setup2 p6 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_S2=ladder MLUCAS_FAULT_S2_Q=1000000; S=$d/p$P.stat
	expect_grep "$S" "FAULT INJECTION: stage 2 ladder perturbed" "injection fired"
	expect_grep "$S" "stage 2 ladder check FAILED" "ladder mismatch detected"
	[[ $(cat "$d/exit") != 0 ]] && ok "run stopped (exit $(cat "$d/exit")) with the .s2 checkpoint intact" || bad "run did not stop"
	[[ -f $d/p$P.s2 ]] && ok ".s2 checkpoint present" || bad ".s2 checkpoint missing"
	run "$d"	# restart: rebuilds ladder and tables from the stage 1 residue, resumes the accumulator from .s2
	# (A restarted stage 2 rebuilds its prime-pairing window from the checkpoint q, so its final accumulator is a
	# legitimately different product from an uninterrupted run's - the same is true of a plain interrupt/resume -
	# hence recovery is asserted through completion and clean checks, not through accumulator equality.)
	[[ $(grep -c "stage 2 ladder check FAILED" "$S") -eq 1 ]] && ok "no repeat failure on restart" || bad "ladder failure recurred after restart"
	expect_grep "$d/results.txt" '"B2":2000000' "restarted stage 2 completed"
	expect_count "$S" "Stage 2 ladder check passed at q = " 3 "ladder checks pass after the restart"

	echo "-- P7: stage 2 table entry corrupted: caught by the checksum; restart finishes clean"
	d=$(setup2 p7 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_S2=table MLUCAS_FAULT_S2_Q=1000000; S=$d/p$P.stat
	expect_grep "$S" "FAULT INJECTION: stage 2 table perturbed" "injection fired"
	expect_grep "$S" "stage 2 table entry 1 of [0-9]* fails its checksum" "checksum mismatch detected on the corrupted entry"
	run "$d"
	[[ $(grep -c "fails its checksum" "$S") -eq 1 ]] && ok "no repeat failure on restart" || bad "checksum failure recurred after restart"
	expect_grep "$d/results.txt" '"B2":2000000' "restarted stage 2 completed"

	echo "-- P8: stage 2 accumulator corrupted: NOT caught (no invariant exists) - the run ends with a different result"
	d=$(setup2 p8 "$MLUCAS_FI"); run "$d" MLUCAS_FAULT_S2=accum MLUCAS_FAULT_S2_Q=1000000; S=$d/p$P.stat
	expect_grep "$S" "FAULT INJECTION: stage 2 accum perturbed" "injection fired"
	expect_nogrep "$S" "ERROR: M$P stage 2" "nothing catches it - this is the documented limit of the stage 2 checks"
	[[ $(s2_final "$d" | tail -1) != "$CLEAN2" ]] && ok "run ends with a WRONG accumulator $(s2_final "$d" | tail -1) - documented, not hidden" || bad "accumulator equals the clean one?!"
else
	skip "P6-P8: stage 2 fault-injection tests need a -DMLUCAS_FAULT_INJECT build"
fi

echo "== $PASS passed, $FAIL failed, $SKIP skipped (work dir $WORK)"
(( FAIL == 0 ))
