#!/usr/bin/env bash
#
# Build every SIMD mode this compiler and architecture support, for both the
# Mlucas and Mfactor targets, and run --help only on the modes this host can
# actually execute.
#
# Two separate questions, deliberately kept apart:
#
#   Can the COMPILER emit it?  -> decides whether we build.  Probed, because
#       gcc >= 16 and clang >= 19 have dropped -march=knl, while older ones
#       still accept it.
#   Can the HOST run it?       -> decides whether we execute the binary.  A
#       runner may lack AVX-512 (or have it masked off by a hypervisor), and
#       running such a binary is an immediate SIGILL that says nothing about
#       the code.
#
# Every build failure is a failure.  There is no '|| true' here: masking is
# what let 'P2WORD may not be used together with USE_FMADD' sit unnoticed on
# every FMA-capable host, and what let the K1OM job report success while
# producing no Mfactor binary at all.
#
set -u -o pipefail

cd "$(dirname "$0")/../.." || exit 1

CC=${CC:-gcc}
MFAC_WORDS=('' 1word 2word 3word 4word nword)

case "$(uname -m)" in
    x86_64|amd64) CAND=(nosimd sse2 avx avx2 avx512 avx512_knl) ;;
    aarch64|arm64) CAND=(nosimd asimd) ;;
    *)            CAND=(nosimd) ;;
esac

mode_cflags() {
    case $1 in
        nosimd)     echo "" ;;
        sse2)       echo "-DUSE_SSE2 -msse2" ;;
        avx)        echo "-DUSE_AVX -mavx" ;;
        avx2)       echo "-DUSE_AVX2 -mavx2 -mfma" ;;
        avx512)     echo "-DUSE_AVX512 -mavx512f -mavx512cd -mavx512dq -mavx512bw -mavx512vl -mfma" ;;
        avx512_knl) echo "-DUSE_AVX512 -march=knl -mavx512f -mavx512cd -mavx512er -mfma" ;;
        asimd)      echo "-DUSE_ARM_V8_SIMD" ;;
    esac
}

# Can the compiler emit this mode at all?
compiler_supports() {
    # shellcheck disable=SC2046
    $CC $(mode_cflags "$1") -x c -c /dev/null -o /dev/null 2>/dev/null
}

# Can this host execute a binary built for this mode?
native=$($CC -dM -E -x c /dev/null 2>/dev/null)
host_supports() {
    case $1 in
        nosimd)     return 0 ;;
        sse2)       grep -q '__SSE2__'    <<<"$native" ;;
        avx)        grep -q '__AVX__'     <<<"$native" ;;
        avx2)       grep -q '__AVX2__'    <<<"$native" ;;
        avx512*)    grep -q '__AVX512F__' <<<"$native" ;;
        asimd)      grep -q '__ARM_NEON'  <<<"$native" ;;
        *)          return 1 ;;
    esac
}

MODES=()
SKIPPED=()
for m in "${CAND[@]}"; do
    if compiler_supports "$m"; then MODES+=("$m"); else SKIPPED+=("$m"); fi
done

echo "compiler : $($CC --version | head -1)"
echo "arch     : $(uname -m)"
echo "building : ${MODES[*]}"
[[ ${#SKIPPED[@]} -gt 0 ]] && echo "unsupported by this compiler, not built: ${SKIPPED[*]}"
echo

fail=0
results=()

record() { results+=("$1"); [[ $1 == FAIL* ]] && fail=1; return 0; }

for mode in "${MODES[@]}"; do
    echo "=== $mode ==="

    # makemake.sh refuses to build over an existing target ("Error: ... already
    # exists"), and 'make clean' leaves the binary behind -- so a second local run
    # would fail on every mode. Start each build from an empty object dir.
    rm -rf "obj${mode:+_$mode}"
    if bash -e -o pipefail -- makemake.sh "$mode"; then
        record "PASS  Mlucas   $mode"
    else
        record "FAIL  Mlucas   $mode"
    fi
    dir="obj${mode:+_$mode}"
    [[ -d $dir ]] && (cd "$dir" && make clean >/dev/null 2>&1)

    for word in "${MFAC_WORDS[@]}"; do
        # An empty word must be omitted, not passed as an empty argument:
        # makemake.sh counts positional parameters and rejects a blank one.
        rm -rf "obj_mfac${word:+_$word}${mode:+_$mode}"
        args=(mfac)
        [[ -n $word ]] && args+=("$word")
        args+=("$mode")
        if bash -e -o pipefail -- makemake.sh "${args[@]}"; then
            record "PASS  Mfactor $mode ${word:-base}"
        else
            record "FAIL  Mfactor $mode ${word:-base}"
            continue
        fi
        dir="obj_mfac${word:+_$word}${mode:+_$mode}"
        bin="$dir/Mfactor${word:+_$word}"
        if [[ -x $bin ]] && host_supports "$mode"; then
            "./$bin" -h >/dev/null || record "FAIL  run -h  $mode ${word:-base}"
        fi
        [[ -d $dir ]] && (cd "$dir" && make clean >/dev/null 2>&1)
    done
done

echo
echo "==================== build matrix ===================="
printf '%s\n' "${results[@]}"
echo "======================================================"
printf '%s built, %s failed\n' \
    "$(printf '%s\n' "${results[@]}" | grep -c '^PASS')" \
    "$(printf '%s\n' "${results[@]}" | grep -c '^FAIL')"

exit "$fail"
