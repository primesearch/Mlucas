#!/usr/bin/env bash
#
# Build Mlucas for one SIMD mode under AddressSanitizer and assert that it computes the
# *right answers* - not merely that it emits no sanitizer report.
#
# Usage:  .github/scripts/asan-simd-check.sh <avx2|avx512> [runner...]
#
#   <mode>       SIMD build mode, passed straight through to makemake.sh.
#   [runner...]  Optional launcher prefixed to the binary, e.g. "$HOME/sde/sde64 -skx --".
#                Needed for avx512 because GitHub's Azure runner SKUs are inconsistent
#                about AVX-512 support; Intel SDE removes that dependency entirely.
#
# Reads CC / CFLAGS / LDFLAGS from the environment (makemake.sh's generated Makefile uses
# "CFLAGS ?=", so an exported CFLAGS wins outright). Do NOT export CPPFLAGS: makemake.sh
# sets "CPPFLAGS ?= -D_GNU_SOURCE ...", and without -D_GNU_SOURCE threadpool.c loses
# CPU_ZERO/CPU_SET/sched_setaffinity on current glibc.
#
# WHY A RESIDUE ASSERTION, AND NOT A "grep for ERROR: AddressSanitizer" JOB
# ------------------------------------------------------------------------
# The defect class this job exists to catch is the undeclared inline-asm clobber: an asm
# block writes a GPR, opmask or high vector register that it does not name in its clobber
# list, so the compiler keeps a live value there across the block. AddressSanitizer cannot
# see that at all - no invalid memory access ever occurs, the arithmetic is simply wrong.
# Two such defects were live on Mlucas main simultaneously and produced zero sanitizer
# reports while silently returning wrong residues for most AVX-512 FFT lengths a real user
# would pick. A report-grepping job passes such a tree with flying colours.
#
# Consequently:
#   * every case is checked against a pinned Res64, and a mismatch fails the job;
#   * the FFT radices actually selected are checked too, so a future change to
#     get_fft_radices() cannot silently move a case off the code path it was chosen for;
#   * the exponent is pinned with -m, so the reference residues cannot drift with the
#     self-test tables;
#   * a missing Res64 line is a failure. Exit status alone is NOT sufficient: Mlucas
#     prints "ERROR ERROR...Halting test" and then still exits 0.
#
# DO NOT ADD -fsanitize-recover=address. It is not a neutral convenience flag here: adding
# it and nothing else makes a *broken* tree produce correct results, on both clang 18 and
# clang 22 - i.e. it hides precisely the defect class this job exists for.
#
set -euo pipefail

MODE=${1:?usage: $0 <avx2|avx512> [runner...]}
shift
RUNNER=("$@")

# Per-case timeout. SIGKILL, not the default SIGTERM: a process running under Intel SDE
# (Pin/ptrace) can end up parked where SIGTERM is never delivered and a plain `timeout` never
# fires - one such run sat wedged for 13 hours. This is not hypothetical for this job in
# particular: on a *broken* tree the SIGSEGV case sometimes wedges Pin instead of dying, which
# is precisely when the timeout has to work.
#
# 600 s is ~15x the slowest passing case measured (128K under SDE with ASan), and 5 x 600 s
# still fits inside the job's timeout-minutes, so a fully wedged run fails the job rather than
# burning the whole runner allocation.
CASE_TIMEOUT=${CASE_TIMEOUT:-600}

# fftlen(Kdbl) | radset | exponent | expected radices | expected Res64 at 100 iters, shift 0
#
# Why these five cells, and no others - each is seconds to a couple of minutes, even under SDE:
#
#   32/0   32 32 16     Trailing radix 16. Control: passes even on a broken tree, so a red
#                       rs1/rs2 beside a green rs0 localises the fault to one routine.
#   32/1   32 16 32     Trailing radix 32 -> selects radix32_wrapper_square. THE discriminator
#   32/2   16 32 32     for the AVX-512 clobber defect in radix32_wrapper_square_gcc64.h. The
#                       discriminator is the *trailing* radix, not the leading one, so a matrix
#                       built from default radix sets (or one FFT length per leading radix) can
#                       miss this entirely - hence the explicit -radset. Both fail on unfixed
#                       main under clang/avx512 (roundoff error at iteration 1 / SIGSEGV).
#   32/3   16 8 8 16    Four-pass, trailing 16. Second control, different pass structure.
#   128/0  128 16 32    The AVX-512 opmask/high-vector carry-macro clobbers in carry_gcc64.h.
#                       gcc-only in practice - clang never reached it in any configuration -
#                       so this cell is why the gcc leg of the matrix cannot be dropped as
#                       redundant. Fails on unfixed main under gcc/avx512 with a wild wtA
#                       pointer.
#
# All four 32K radsets share one exponent and therefore one expected residue, which is what
# makes them a clean same-binary A/B rather than five unrelated numbers.
CASES=(
	"32|0|673469|32 32 16|1A4EF8A0D172FBAA"
	"32|1|673469|32 16 32|1A4EF8A0D172FBAA"
	"32|2|673469|16 32 32|1A4EF8A0D172FBAA"
	"32|3|673469|16 8 8 16|1A4EF8A0D172FBAA"
	"128|0|2614999|128 16 32|040918890E98F8DA"
)

OBJ=obj_$MODE
rm -rf -- "$OBJ"

echo "::group::Build $MODE under AddressSanitizer"
echo "CC=${CC:-gcc}"
echo "CFLAGS=${CFLAGS:-<makemake default>}"
bash -e -o pipefail -- makemake.sh "$MODE"
echo "::endgroup::"

# A silently-unsanitised binary must not be able to pass as a sanitised one.
nsan=$(nm "$OBJ/Mlucas" | grep -c __asan_ || true)
echo "$OBJ/Mlucas: $nsan __asan_ symbols"
if [[ $nsan -eq 0 ]]; then
	echo "FAIL: $OBJ/Mlucas contains no __asan_ symbols - it was not instrumented." >&2
	exit 1
fi

rc=0
for case in "${CASES[@]}"; do
	IFS='|' read -r len radset expo want_radices want_res64 <<<"$case"

	# Inside obj_*/, which .gitignore already covers, so a local run leaves no stray files.
	log="$OBJ/asan_${MODE}_${len}_${radset}.log"
	echo "::group::$MODE  M$expo  fftlen ${len}K  radset $radset  (expect radices $want_radices)"
	crc=0
	ASAN_OPTIONS=detect_leaks=0:abort_on_error=0 \
		timeout -s KILL "$CASE_TIMEOUT" \
		"${RUNNER[@]}" "./$OBJ/Mlucas" -m "$expo" -fftlen "$len" -radset "$radset" \
		-shift 0 -iters 100 -cpu 0 >"$log" 2>&1 || crc=$?
	cat -- "$log"
	echo "::endgroup::"

	got_res64=$(awk '{ for (i = 1; i <= NF; i++) if ($i == "Res64:") { r = $(i+1); sub(/\..*/, "", r); print r; exit } }' "$log")
	got_radices=$(awk -F'radices' '/Using complex FFT radices/ { n = split($2, a, " "); s = a[1]; for (i = 2; i <= n; i++) s = s " " a[i]; print s; exit }' "$log")

	if [[ -z $got_radices ]]; then
		echo "FAIL $MODE fftlen ${len}K radset $radset: died before printing its radices (exit $crc)." >&2
		grep -E 'Roundoff warning|SDE ERROR|unaligned memory reference|ERROR: AddressSanitizer' -- "$log" >&2 || true
		rc=1
		continue
	fi
	if [[ $got_radices != "$want_radices" ]]; then
		echo "FAIL $MODE fftlen ${len}K radset $radset: radices [$got_radices] != expected [$want_radices]." >&2
		echo "     This case no longer exercises the code path it was chosen for; re-pin it." >&2
		rc=1
		continue
	fi
	if [[ -z $got_res64 ]]; then
		echo "FAIL $MODE fftlen ${len}K radset $radset: no Res64 line (exit $crc)." >&2
		grep -E 'Roundoff warning|SDE ERROR|unaligned memory reference|ERROR: AddressSanitizer' -- "$log" >&2 || true
		rc=1
		continue
	fi
	if [[ $got_res64 != "$want_res64" ]]; then
		echo "FAIL $MODE fftlen ${len}K radset $radset: Res64 $got_res64 != expected $want_res64." >&2
		rc=1
		continue
	fi
	# Secondary, and genuinely secondary: a correct residue with a sanitizer report is still
	# a failure. Kept last so its absence can never be mistaken for a pass.
	if grep -q 'ERROR: AddressSanitizer' -- "$log"; then
		echo "FAIL $MODE fftlen ${len}K radset $radset: correct Res64 but AddressSanitizer reported:" >&2
		grep -A2 'ERROR: AddressSanitizer' -- "$log" >&2
		rc=1
		continue
	fi
	echo "PASS $MODE fftlen ${len}K radset $radset  radices [$got_radices]  Res64 $got_res64"
done

exit "$rc"
