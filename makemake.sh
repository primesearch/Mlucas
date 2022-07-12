#!/bin/bash

# EWM: This is a makefile-only cutdown of Teal Dulcet's much more extensive Mlucas install/build/tune
# script, available at https://raw.github.com/tdulcet/Distributed-Computing-Scripts/master/mlucas.sh ;
# he does not explicitly use the GPL boilerplate as below, but assures me his version is GPL-covered.

################################################################################
#                                                                              #
#   (C) 2021 by Ernst W. Mayer.                                                #
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

# Exit if any of the sommands fail:
set -e

# for mode in avx512 avx2 avx sse2; do
	# if grep -iq "$mode" /proc/cpuinfo; then
		# echo -e "The CPU supports the ${mode^^} SIMD build mode.\n"
		# ARGS+=( "-DUSE_${mode^^}" )
		# break
	# fi
# done

# $0 contains script-name, but $@ starts with first ensuing cmd-line arg, if it exists:
echo "Total number of arguments : $#"
for ((i=1; i<=$#; ++i)); do
	echo "Arg[$i] = ${!i}"
done
if [[ $# -gt 1 ]]; then
	echo "Usage: $0 [SIMD build mode]" >&2
	echo "Only 1 optional argument supported, it must be one of the supported SIMD-arithmetic types:" >&2
	echo -e "\t[x86_64: avx512_skylake avx512_knl avx2 avx sse2]; [Armv8: asimd]; or 'nosimd' for scalar-double build.\n" >&2
	exit 1
fi

DIR=obj
EXE=Mlucas
ARGS=()

if echo "$OSTYPE" | grep -iq 'darwin'; then
	echo -e "MacOS detected for build host.\n"
	CPU_THREADS=$(sysctl -n hw.ncpu)
else # echo "$OSTYPE" | grep -iq 'linux'
	echo -e "Assuming OS = Linux for build host.\n"
	CPU_THREADS=$(nproc --all)
fi

# Thx to tdulcet for streamlined case-based syntax here, but ugh - non-matching ')', really?:
if [[ $# -eq 1 ]]; then

	if [ "$1" = 'avx512_skylake' ]; then
		echo "Building for avx512_skylake SIMD in directory obj_$1; the executable will be named Mlucas_$1"
		ARGS+=( "-DUSE_AVX512" -march=skylake-avx512 )
	elif [ "$1" = 'avx512_knl' ]; then
		echo "Building for avx2 SIMD in directory obj_$1; the executable will be named Mlucas_$1"
		ARGS+=( "-DUSE_AVX512" -march=knl )
	elif [ "$1" = 'avx2' ]; then
		echo "Building for avx2 SIMD in directory obj_$1; the executable will be named Mlucas_$1"
		ARGS+=( "-DUSE_AVX2" -mavx2 )
	elif [ "$1" = 'avx' ]; then
		echo "Building for avx SIMD in directory obj_$1; the executable will be named Mlucas_$1"
		ARGS+=( "-DUSE_AVX" -mavx )
	elif [ "$1" = 'sse2' ]; then
		echo "Building for sse2 SIMD in directory obj_$1; the executable will be named Mlucas_$1"
		ARGS+=( "-DUSE_SSE2" )
	elif [ "$1" = 'asimd' ]; then
		echo "Building for avx2 SIMD in directory obj_$1; the executable will be named Mlucas_$1"
		ARGS+=( "-DUSE_ARM_V8_SIMD" )
	elif [ "$1" = 'nosimd' ]; then
		echo "Building in scalar-double (no-SIMD) mode in directory obj_$1; the executable will be named Mlucas_$1"
		# This one's a no-op
	else
		echo "Unrecognized SIMD-build flag ... aborting."
		exit 1
	fi

	DIR+="_$1"
	EXE+="_$1"

elif uname -a | grep -iq 'Mac'; then

	echo -e "MacOS detected.\n"
	if sysctl -a | grep machdep.cpu.features | grep -iq 'avx512'; then
		echo -e "The CPU supports the AVX512 SIMD build mode.\n"
		ARGS+=( "-DUSE_AVX512" -march=native )
	elif sysctl -a | grep machdep.cpu.features | grep -iq 'avx2'; then
		echo -e "The CPU supports the AVX2 SIMD build mode.\n"
		ARGS+=( "-DUSE_AVX2" -march=native -mavx2 )
	elif sysctl -a | grep machdep.cpu.features | grep -iq 'avx'; then
		echo -e "The CPU supports the AVX SIMD build mode.\n"
		ARGS+=( "-DUSE_AVX" -march=native -mavx )
	elif sysctl -a | grep machdep.cpu.features | grep -iq 'sse2'; then
		echo -e "The CPU supports the SSE2 SIMD build mode.\n"
		ARGS+=( "-DUSE_SSE2" -march=native )
	elif sysctl -a | grep machdep.cpu.features | grep -iq 'asimd'; then
		echo -e "The CPU supports the ASIMD build mode.\n"
		ARGS+=( "-DUSE_ARM_V8_SIMD" -march=native )
	else
		echo -e "The CPU supports no Mlucas-recognized ASIMD build mode ... building in scalar-double mode.\n"
		ARGS+=( -march=native )
	fi

else

	echo -e "Assuming OS = Linux.\n"
	if grep -iq 'avx512' /proc/cpuinfo; then
		echo -e "The CPU supports the AVX512 SIMD build mode.\n"
		ARGS+=( "-DUSE_AVX512" -march=native )
	elif grep -iq 'avx2' /proc/cpuinfo; then
		echo -e "The CPU supports the AVX2 SIMD build mode.\n"
		ARGS+=( "-DUSE_AVX2" -march=native -mavx2 )
	elif grep -iq 'avx' /proc/cpuinfo; then
		echo -e "The CPU supports the AVX SIMD build mode.\n"
		ARGS+=( "-DUSE_AVX" -march=native -mavx )
	elif grep -iq 'sse2' /proc/cpuinfo; then
		echo -e "The CPU supports the SSE2 SIMD build mode.\n"
		ARGS+=( "-DUSE_SSE2" -march=native )
	elif grep -iq 'asimd' /proc/cpuinfo; then
		echo -e "The CPU supports the ASIMD build mode.\n"
		ARGS+=( "-DUSE_ARM_V8_SIMD" -march=native )
	else
		echo -e "The CPU supports no Mlucas-recognized ASIMD build mode ... building in scalar-double mode.\n"
		ARGS+=( -march=native )
	fi

fi

# -p prevents "File exists" warning if obj-dir already exists:
mkdir -p "$DIR"
cd "$DIR"

# Clang-under-MacOS linker barfs if one tries to explicitly invoke standard libs - h/t tdulcet for the
# conditional-inline syntax. Some OSes put the GMP headers in /usr/local/include, so -I that path in the
# compile command. If said path does not exist, make silently ignores it.
# Re. the -g flag to include the debugging symbols, they bloat executable size but if someone's Mlucas
# crashes/segfaults, one can rerun with GDB (gdb -ex=r ./Mlucas) to see the filename, line number and
# stack trace of the issue. If one wishes, one can run 'strip -g Mlucas' to remove the debugging symbols:
cat << EOF > Makefile
CC?=gcc
OBJS=\$(patsubst ../src/%.c, %.o, \$(wildcard ../src/*.c))

$EXE: \$(OBJS)
	\$(CC) -Wall -g -o \$@ \$(OBJS) $(echo "$OSTYPE" | grep -iq 'darwin' || echo "-lm -lpthread -lrt") -lgmp
%.o: ../src/%.c
	\$(CC) -Wall -g -c -I/usr/local/include -O3 ${ARGS[@]} -DUSE_THREADS \$<
clean:
	rm -f *.o
EOF

echo -e "\nBuilding Mlucas"
printf "%s CPU cores detected ... parallel-building using that number of make threads.\n" "$CPU_THREADS"
if ! make -j "$CPU_THREADS" > build.log 2>&1; then
	echo -e "There were build errors - see build.log for details.\n"
	exit 1
fi
