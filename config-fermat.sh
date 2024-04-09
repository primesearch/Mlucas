#!/bin/bash

# Shell script for generating fermat.cfg; Mlucas output saved to config-fermat.log

################################################################################
#                                                                              #
#   (C) 2024 by Catherine Cowie and Teal Dulcet.                               #
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

# Mlucas
MLUCAS=./Mlucas

# Number of iterations
# use 100, 1000, or 10000 to match pre-computed values
ITERS=100

# Minimum Fermat number (14 or greater)
MIN=14

# Maximum Fermat number (33 or less)
MAX=29

# Mlucas arguments
ARGS=(
	"$@"
	# Add desired -cpu or -core settings here, or as following arguments, e.g. ../config-fermat.sh -cpu 0:3

)

# First, tiny FFT lengths for F14 to F17;
FFTS=([1]=14 [2]=15 [4]=16 [7]=17 [8]=17)
# Then, from small up to egregious FFTs for F18 to F33.
# The largest FFT reached is 512M, if MAX is set to 33.
# Note that large FFTs require considerable runtime at 10000 iterations.
for ((n = 0; n < 16; ++n)); do
	m=$((1 << n))
	f=$((18 + n))
	for k in 15 16; do
		if [[ $k -eq 15 && $n -lt 11 ]]; then
			# k = 7 multiples (7K, 14K, ...) become unworkable after F28 (14M).
			FFTS[14 * m]=$f
		fi
		# k = 15, 16 should both be supported up to at least F32.
		FFTS[k * m]=$f
		if [[ $k -eq 15 && $n -gt 5 ]]; then
			# k = 63 is supported for F24 and above (1008K).
			FFTS[63 * m >> 2]=$f
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
	printf '\n\tTesting F%s (2^%s + 1),\tFFT length: %sK\n\n' "$f" $((1 << f)) "$fft"
	args=("${ARGS[@]}")
	if [[ $f -le 17 || $f -ge 32 ]]; then
		args+=(-shift 0)
	fi
	time $MLUCAS -f "$f" -fft "$fft" -iters $ITERS "${args[@]}" 2>&1 | tee -a config-fermat.log | grep -i 'error\|warn\|assert\|writing\|pmax_rec\|fft radices'
done
