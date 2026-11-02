#!/bin/bash

# Shell script for checking that tests started at an explicit residue shift (-shift) give the right
# residue, e.g. bash ../check-shifted-seed.sh ./Mlucas -cpu 0

################################################################################
#                                                                              #
#   (C) 2026 by Mark Rose.                                                     #
#                                                                              #
#  This program is free software; you can redistribute it and/or modify it     #
#  under the terms of the GNU General Public License as published by the       #
#  Free Software Foundation; either version 2 of the License, or (at your      #
#  option) any later version.                                                  #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but WITHOUT #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   #
#  more details.                                                               #
#                                                                              #
#  You should have received a copy of the GNU General Public License along     #
#  with this program; see the file GPL.txt.  If not, you may view one at       #
#  http://www.fsf.org/licenses/licenses.html, or obtain one by writing to the  #
#  Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA     #
#  02111-1307, USA.                                                            #
#                                                                              #
################################################################################

# Usage: bash check-shifted-seed.sh [Mlucas binary or directory containing it] [Mlucas arguments...]
#
# The shifted initial seed (4 for LL, the PRP base, 3 for Pepin) is stored in the one residue word
# that holds bit -shift. Which word that is, and where it sits in the SIMD-permuted data array,
# depends on the FFT radices and the instruction set, so a window of shifts about one word wide
# misbehaves on a build with a bug in the seeding code, and the window moves from build to build.
# The sweeps below are therefore ranges that cover the window on every layout, and include the
# shifts that put the seed into the imaginary slot of the complex datum for Fermat-mod.
# Each run takes about a tenth of a second, nearly all of it Mlucas start-up, so the runs are spread over
# all the CPUs (-cpu is not needed: Mlucas runs single-threaded when told nothing). If the arguments pin
# the threads with -cpu, -core or -nthread, the runs happen one at a time instead, as they would all land
# on the same cores. Set JOBS to override the number of simultaneous runs.
#
# Every run is 100 iterations with a fixed radix set (-radset), and the Res64 is compared with the
# pre-computed value that Mlucas itself checks for -iters 100. Each run happens in a fresh temporary
# directory, which is removed afterwards.
# Set MLUCAS_LAUNCHER to run the binary under an emulator, e.g. MLUCAS_LAUNCHER='sde64 -skx --'.
# Exits non-zero, naming each case and shift, if any run aborts or gives a different Res64.

# Mlucas
MLUCAS=./Mlucas

# Number of iterations (the expected residues below are for 100)
ITERS=100

if [[ $# -gt 0 && $1 != -* ]]; then
	MLUCAS=$1
	shift
fi
if [[ -d $MLUCAS ]]; then
	MLUCAS=$MLUCAS/Mlucas
fi
if [[ ! -x $MLUCAS ]]; then
	echo "Error: Mlucas binary '$MLUCAS' not found or not executable." >&2
	exit 2
fi
# The runs happen in other directories, so make the path absolute
MLUCAS=$(cd -- "$(dirname -- "$MLUCAS")" && pwd)/$(basename -- "$MLUCAS")

# Mlucas arguments
ARGS=(
	"$@"
	# Add desired -cpu or -core settings here, or as following arguments, e.g. bash ../check-shifted-seed.sh ./Mlucas -cpu 0
)

# Optional launcher (emulator) command
LAUNCHER=()
if [[ -n ${MLUCAS_LAUNCHER-} ]]; then
	read -r -a LAUNCHER <<<"$MLUCAS_LAUNCHER"
fi

if [[ -z ${JOBS-} ]]; then
	JOBS=$(nproc 2>/dev/null || echo 1)
	for arg in "${ARGS[@]}"; do
		if [[ $arg == -cpu || $arg == -core || $arg == -nthread ]]; then
			JOBS=1
		fi
	done
fi

TMPDIR_ROOT=$(mktemp -d)
cleanup() {
	local pids
	mapfile -t pids < <(jobs -p)
	if [[ ${#pids[@]} -gt 0 ]]; then
		kill -- "${pids[@]}" 2>/dev/null || true
		wait 2>/dev/null || true
	fi
	rm -rf -- "$TMPDIR_ROOT"
}
trap cleanup EXIT

RUNS=0
declare -A SKIPPED
TOLERATE_BASELINE_FAILURE=

# run_one INDEX NAME RESULT_LABEL EXPECTED_RES64 SHIFT ARGS...
# Writes a message to "$TMPDIR_ROOT/fail.INDEX" if the run is wrong.
# RESULT_LABEL is the text before "Res64:" on the result line: empty for LL and Fermat, "PRP-3 " for PRP.
run_one() {
	local index=$1 name=$2 label=$3 expected=$4 shift_count=$5
	shift 5
	local dir output status=0 line found problems
	dir=$(mktemp -d -p "$TMPDIR_ROOT")
	output=$(cd -- "$dir" && "${LAUNCHER[@]}" "$MLUCAS" "$@" -iters "$ITERS" "${ARGS[@]}" -shift "$shift_count" 2>&1) || status=$?
	rm -rf -- "$dir"
	# grep exits 1 when it finds nothing, which is the normal case. Under 'bash -e -o pipefail', as CI runs this script,
	# that would end the function here, before the residue is checked, so a clean run has to be allowed to match nothing:
	problems=$(grep -i 'assert\|error\|warn' <<<"$output" | head -n 2 | tr '\n' ' ') || true
	if [[ $status -ne 0 ]]; then
		echo "$name, -shift $shift_count: Mlucas exited with status $status: $problems" >"$TMPDIR_ROOT/fail.$index"
	elif [[ -n $problems ]]; then
		echo "$name, -shift $shift_count: $problems" >"$TMPDIR_ROOT/fail.$index"
	else
		found=
		while IFS= read -r line; do
			if [[ $line == "${label}Res64: "* ]]; then
				found=${line#"${label}Res64: "}
				found=${found%%.*}
				break
			fi
		done <<<"$output"
		if [[ $found != "$expected" ]]; then
			echo "$name, -shift $shift_count: Res64 ${found:-not printed}, expected $expected" >"$TMPDIR_ROOT/fail.$index"
		fi
	fi
}

# check_case NAME RESULT_LABEL EXPECTED_RES64 SHIFT_START SHIFT_END ARGS...
# Runs Mlucas with -shift from SHIFT_START to SHIFT_END inclusive.
check_case() {
	local name=$1 label=$2 expected=$3 start=$4 end=$5
	shift 5
	local shift_count
	if [[ -n ${SKIPPED[$name]-} ]]; then
		return
	fi
	if [[ $start -eq 0 && -n $TOLERATE_BASELINE_FAILURE ]]; then
		# Without a shift, first: if that is already wrong, the radix set is broken for another reason
		run_one base "$name" "$label" "$expected" 0 "$@"
		if [[ -e $TMPDIR_ROOT/fail.base ]]; then
			SKIPPED[$name]=$(<"$TMPDIR_ROOT/fail.base")
			rm -f -- "$TMPDIR_ROOT/fail.base"
			echo "SKIPPED: ${SKIPPED[$name]} (it fails without a shift too, so this is not a shifted-seed problem)" >&2
			if [[ -n ${GITHUB_ACTIONS-} ]]; then
				echo "::warning::check-shifted-seed.sh skipped ${SKIPPED[$name]}"
			fi
			return
		fi
	fi
	for ((shift_count = start; shift_count <= end; ++shift_count)); do
		while [[ $(jobs -r | wc -l) -ge $JOBS ]]; do
			wait -n || true	# A finished job's own status is not an error here: failures are recorded in files
		done
		RUNS=$((RUNS + 1))
		run_one "$(printf '%05d' "$RUNS")" "$name" "$label" "$expected" "$shift_count" "$@" &
	done
}

start=$SECONDS

# Mersenne M86243 at 4K, radices 16,8,16: the shifted seed lands in physical index 1 of the data array
# at shifts 43-63 with SSE2/AVX2/AVX-512 and at 22-42 with scalar double. 86243 / 2 = 43121.
for s in 0:80 43121:43121 86242:86242; do
	check_case "M86243 LL 4K" "" CA7BFFCCF8FBC07C "${s%:*}" "${s#*:}" -m 86243 -fft 4 -radset 0
	check_case "M86243 PRP-3 4K" "PRP-3 " 7AD1A0FCFC20DD4F "${s%:*}" "${s#*:}" -m 86243 -fft 4 -radset 0 -prp
done

# Fermat F16 at 4K (power-of-2 length, no DWT weights): shifts 19-31 with AVX2 and the twist applied twice.
# 65536 / 2 = 32768, from where the seed is in the imaginary slot of the complex datum: 32787-32799 with AVX2.
for s in 0:64 32768:32832 65535:65535; do
	check_case "F16 4K" "" AAE76C15C2B37465 "${s%:*}" "${s#*:}" -f 16 -fft 4 -radset 0
done

# Fermat F18 at 14K (non-power-of-2 length, with DWT weights), radices 56,8,16: shifts 19-36 and 131091-131108
# with AVX2. 262144 / 2 = 131072.
# AVX-512 builds get a fatal roundoff error at iteration 14 on this radix set at every shift, including 0
# (radix56_main_carry_loop.h, fixed by PR #319), and no other non-power-of-2 length is usable for F18 with
# AVX-512 (radix 28 and 60 have no AVX-512 Fermat-mod routine). Such a build says so and skips this case
# rather than reporting a failure that has nothing to do with the shift. Delete this line once #319 is in.
TOLERATE_BASELINE_FAILURE=1
for s in 0:64 131072:131136 262143:262143; do
	check_case "F18 14K" "" 7C6B681485EB86DB "${s%:*}" "${s#*:}" -f 18 -fft 14 -radset 0
done
TOLERATE_BASELINE_FAILURE=

wait

shopt -s nullglob
failed=("$TMPDIR_ROOT"/fail.*)
if [[ ${#failed[@]} -ne 0 ]]; then
	for file in "${failed[@]}"; do
		echo "FAIL: $(<"$file")" >&2
	done
	echo "${#failed[@]} of $RUNS shifted-seed runs FAILED." >&2
	exit 1
fi
echo "All $RUNS shifted-seed runs gave the expected residue in $((SECONDS - start)) s${SKIPPED[*]:+ (${#SKIPPED[@]} case skipped, see above)}."
