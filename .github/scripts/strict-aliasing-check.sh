#!/usr/bin/env bash
# Compile every translation unit the way makemake.sh builds it - each x86 SIMD mode, and the Mfactor
# sources for each word size - with gcc's -Wstrict-aliasing=2, and fail if gcc reports any
# type-punned access that "will" or "might" break the strict-aliasing rules.
#
# Level 2, not the default: the -Wstrict-aliasing=3 that -Wall turns on reports none of the
# violations this tree has had, including the 14 "will break" sites fixed in the carry loops. Level 2
# finds them, at the cost of false positives, which are worked around in the source rather than
# exempted here, so that the check stays all-or-nothing.
#
# The warning is a heuristic and cannot see every violation - the twopmodq96_q4 miscompile was a
# real strict-aliasing bug that gcc does not warn about at any level - so a clean run here is a floor,
# not a proof. makemake.sh's -fno-strict-aliasing covers the rest.
#
# Usage: .github/scripts/strict-aliasing-check.sh [mode...]   (default: all of the modes below)
# Needs gcc: clang accepts -Wstrict-aliasing but does not implement these checks.
set -u -o pipefail

CC=${CC:-gcc}
cd "$(dirname -- "$0")/../.." || exit 2

declare -A MODE_FLAGS=(
	[nosimd]=''
	[sse2]='-DUSE_SSE2 -msse2'
	[avx]='-DUSE_AVX -mavx'
	[avx2]='-DUSE_AVX2 -mavx2 -mfma'
	[avx512]='-DUSE_AVX512 -mavx512f -mavx512cd -mfma'
)
if (( $# )); then MODES=("$@"); else MODES=(nosimd sse2 avx avx2 avx512); fi

# The object lists, taken from makemake.sh so the two cannot drift apart:
read -ra OBJS <<<"$(sed -n 's/^OBJS=//p' makemake.sh)"
read -ra OBJS_MFAC <<<"$(sed -n 's/^OBJS_MFAC=//p' makemake.sh)"
(( ${#OBJS[@]} && ${#OBJS_MFAC[@]} )) || { echo "Could not read the object lists from makemake.sh" >&2; exit 2; }

BASE=(-std=gnu99 -O2 -D_GNU_SOURCE -DUSE_THREADS -Wstrict-aliasing=2 -fdiagnostics-color=never)
LOG=$(mktemp)
trap 'rm -f "$LOG" "$LOG".*' EXIT

units=()
for mode in "${MODES[@]}"; do
	[[ -n ${MODE_FLAGS[$mode]+x} ]] || { echo "Unknown mode '$mode'" >&2; exit 2; }
	# Mlucas: every object, as makemake.sh's .c.o rule compiles it.
	for o in "${OBJS[@]}"; do
		units+=("$mode|${MODE_FLAGS[$mode]}|src/${o%.o}.c")
	done
	# Mfactor: factor.c for each word size, with makemake.sh's TRYQ choice, and the other Mfactor
	# objects for the multiword sizes (which get -DFACTOR_STANDALONE and the word macro there).
	for w in '' -DP2WORD -DP3WORD -DP4WORD -DNWORD; do
		case $w in -DP4WORD|-DNWORD) tryq= ;; *) tryq=-DTRYQ=4 ;; esac
		units+=("$mode${w:+ $w}|${MODE_FLAGS[$mode]} -DFACTOR_STANDALONE $w $tryq|src/factor.c")
		if [[ -n $w ]]; then
			for o in "${OBJS_MFAC[@]}"; do
				[[ $o == factor.o ]] && continue
				units+=("$mode $w|${MODE_FLAGS[$mode]} $w -DFACTOR_STANDALONE|src/${o%.o}.c")
			done
		fi
	done
done

compile() {	# compile <tag|flags|source>: print any strict-aliasing warnings, prefixed with the tag
	local tag=${1%%|*} rest=${1#*|} flags
	local src=${rest##*|}
	read -ra flags <<<"${rest%|*}"
	"$CC" "${BASE[@]}" "${flags[@]}" -c "$src" -o /dev/null 2>&1 \
		| grep -E 'break strict-aliasing rules' | sed "s|^|[$tag] |"
	return 0
}

# A simple job pool: at most $(nproc) compiles at a time, each writing its own output file.
NJOBS=$(nproc)
for i in "${!units[@]}"; do
	while (( $(jobs -rp | wc -l) >= NJOBS )); do wait -n; done
	compile "${units[$i]}" > "$LOG.$i" &
done
wait
cat "$LOG".* > "$LOG" 2>/dev/null

echo "Compiled ${#units[@]} translation units with $("$CC" --version | head -n1)."
if [[ -s $LOG ]]; then
	n=$(sed -E 's/^\[[^]]*\] //' "$LOG" | sort -u | wc -l)
	echo "::error::$n distinct strict-aliasing warning(s) at -Wstrict-aliasing=2:"
	sort "$LOG" | uniq
	exit 1
fi
echo "No strict-aliasing warnings at -Wstrict-aliasing=2."
