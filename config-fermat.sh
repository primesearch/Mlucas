#!/bin/bash
# Shell script for generating fermat.cfg; Mlucas output saved to config-fermat.log
ITERS=10000     # Number of iterations; use 100, 1000, or 10000 to match pre-computed values
MIN=14          # Minimum Fermat number (14 or greater)
MAX=29          # Maximum Fermat number (33 or less)
ARGS=(
    "$@"
    # Add desired -cpu or -core settings here, or as following arguments, e.g. ./config-fermat.sh -cpu 0:3
)
FFTS=([1]=14 [2]=15 [4]=16 [7]=17 [8]=17)   # First, tiny FFT lengths for F14 to F17;
for ((n = 0; n < 16; ++n)); do              # Then, from small up to egregious FFTs for F18 to F33.
    m=$((1 << n))                           # The largest FFT reached is 512M, if MAX is set to 33.
    f=$((18 + n))                           # Note that large FFTs require considerable runtime at 10000 iterations.
    for k in 15 16; do
        if [[ $k -eq 15 && $n -lt 11 ]]; then
            FFTS[14 * m]=$f                 # k = 7 multiples (7K, 14K, ...) become unworkable after F28 (14M).
        fi
        FFTS[k * m]=$f                      # k = 15, 16 should both be supported up to at least F32.
        if [[ $k -eq 15 && $n -gt 5 ]]; then
            FFTS[63 * m >> 2]=$f            # k = 63 is supported for F24 and above (1008K).
        fi
    done
done
for fft in "${!FFTS[@]}"; do
    f=${FFTS[fft]}
    if [[ -n $MIN && $f -lt $MIN ]]; then
        continue
    elif [[ -n $MAX && $f -gt $MAX ]]; then
        break
    fi
    printf '\n\tTesting F%s (%s),\tFFT length: %sK\n\n' "$f" $((1 << f)) "$fft"
    args=("${ARGS[@]}")
    if [[ $f -le 17 || $f -ge 32 ]]; then
        args+=(-shift 0)
    fi
    time ./Mlucas -f "$f" -fft "$fft" -iters $ITERS "${args[@]}" 2>&1 | tee -a config-fermat.log | grep -i 'error\|warn\|assert\|writing'
done
