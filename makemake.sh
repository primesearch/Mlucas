#!/bin/bash

# EWM: This is a makefile-only cutdown of Teal Dulcet's much more extensive Mlucas install/build/tune
# script, available at https://raw.github.com/tdulcet/Distributed-Computing-Scripts/master/mlucas.sh ;
# he does not explicitly use the GPL boilerplate as below, but assures me his version is GPL-covered.

################################################################################
#                                                                              #
#   (C) 2021 by Ernst W. Mayer and Teal Dulcet.                                #
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

# Exit if any of the commands fail:
set -e

shopt -s nocasematch

# for mode in avx512 avx2 avx sse2; do
	# if grep -iq "$mode" /proc/cpuinfo; then
		# echo -e "The CPU supports the ${mode^^} SIMD build mode.\n"
		# ARGS+=( "-DUSE_${mode^^}" )
		# break
	# fi
# done

DIR=obj
Mlucas=Mlucas
Mfactor=Mfactor
TARGET=$Mlucas
ARGS=(-DUSE_THREADS) # Optional compile args
WORDS=''
# Optional link args
LARG=()

MODES=()
GMP=1
HWLOC=0

if echo "$OSTYPE" | grep -iq 'darwin'; then
	echo -e "MacOS detected for build host.\n"
	CPU_THREADS=$(sysctl -n hw.ncpu)
else # echo "$OSTYPE" | grep -iq 'linux'
	echo -e "Assuming OS = Linux for build host.\n"
	CPU_THREADS=$(nproc --all)
fi

if ! command -v make >/dev/null; then
	echo "Error: This script requires Make" >&2
	echo "On Ubuntu and Debian run: 'sudo apt-get update' and 'sudo apt-get install build-essential -y'" >&2
	exit 1
fi
if [[ -n $CC ]]; then
	if ! command -v "$CC" >/dev/null; then
		echo "Error: $CC is not installed." >&2
		exit 1
	fi
elif ! command -v gcc >/dev/null; then
	echo "Error: This script requires the GNU C compiler" >&2
	echo "On Ubuntu and Debian run: 'sudo apt-get update' and 'sudo apt-get install build-essential -y'" >&2
	exit 1
fi

# $0 contains script-name, but $@ starts with first ensuing cmd-line arg, if it exists:
echo "Total number of input parameters = $#"

# v21: Keep the cross-platform-build arch-specifying command-line flag, but now also need to
# support several added ones for 3rd-party-library usage. This needs to be in arbitrary argument
# order fashion, [details snipped]
arglist=("$@") # Local array into which we copy cmd-line args in order to be able to manipulate them
for i in "${!arglist[@]}"; do
	echo "Arg[$i] = ${arglist[i]}"
done
# Now loop over the optional args and execute the above-described preprocessing step:
for arg in "$@"; do

	case ${arg} in
		'no_gmp')
			GMP=0
			;;
		'use_hwloc')
			HWLOC=1
			;;
		'avx512_skylake' | 'avx512_knl' | 'k1om' | 'avx2' | 'avx' | 'sse2' | 'asimd' | 'nosimd')
			MODES+=("$arg")
			;;
		'mfac')
			TARGET=$Mfactor
			;;
		'1word' | '2word' | '3word' | '4word' | 'nword')
			WORDS=$arg
			;;
		*)
			echo "Usage: $0 [SIMD build mode]" >&2
			echo "Optional arguments must be 'no_gmp', 'use_hwloc' or one and only one of the supported SIMD-arithmetic types:" >&2
			echo -e "\t[x86_64: avx512_skylake avx512_knl k1om avx2 avx sse2]; [Armv8: asimd]; or 'nosimd' for scalar-double build.\n" >&2
			exit 1
			;;
	esac

done

if ((GMP)); then
	LARG+=(-lgmp)
else
	echo "Building sans Gnu-MP ... this means no GCDs will be taken in p-1 work."
	ARGS+=(-DINCLUDE_GMP=0)
fi

if ((HWLOC)); then
	echo "Building with HWLOC hardware-topology support."
	ARGS+=(-DINCLUDE_HWLOC=1)
	LARG+=(-lhwloc)
fi

if [[ $TARGET == "$Mfactor" ]]; then
	DIR+=_mfac
	# trap "rm $PWD/src/factor.c" EXIT
	# cp -vf src/factor.c{.txt,}
fi

if [[ -n $WORDS ]]; then
	if [[ $TARGET == "$Mfactor" ]]; then
		arg=$WORDS
		if [[ ${arg} == 'nword' ]]; then
			WORDS=-DNWORD
		else
			WORDS=-DP"${arg::1}"WORD
		fi
		Mfactor+="_$arg"
		TARGET=$Mfactor
	else
		echo "Error: The argument '$WORDS' requires 'mfac'." >&2
		exit 1
	fi
fi

# First if/elif clause handles cross-platform builds and non-default values for "Use GMP?" and "Use HWLOC?":
# o "Use GMP" = TRUE is default in link step, 'no_gmp' overrides;
# o "Use HWLOC" = FALSE is default, 'use_hwloc' overrides.
# Thx to tdulcet for offering a streamlined case-based syntax here, but ugh - non-matching ')', really?:
if [[ ${#MODES[*]} -gt 1 ]]; then
	echo -e "Only one arch-specifying optional argument is allowed ... aborting." >&2
	exit 1
fi

if [[ ${#MODES[*]} -eq 1 ]]; then

	arg=${MODES[0]}

	case ${arg} in
		'avx512_skylake')
			echo "Building for avx512_skylake SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_AVX512 -march=skylake-avx512)
			;;
		'avx512_knl')
			echo "Building for avx512_knl SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_AVX512 -march=knl)
			;;
		'k1om')
			echo "Building for 1st-gen Xeon Phi 512-bit SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_IMCI512)
			;;
		'avx2')
			echo "Building for avx2 SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_AVX2 -mavx2)
			;;
		'avx')
			echo "Building for avx SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_AVX -mavx)
			;;
		'sse2')
			echo "Building for sse2 SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_SSE2)
			;;
		'asimd')
			echo "Building for asimd SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_ARM_V8_SIMD)
			;;
		'nosimd')
			echo "Building in scalar-double (no-SIMD) mode in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			# This one's a no-op
			;;
		*)
			echo "Unrecognized SIMD-build flag ... aborting." >&2
			exit 1
			;;
	esac

	DIR+="_$arg"

elif echo "$OSTYPE" | grep -iq 'darwin'; then

	# MacOS:
	if sysctl -a | grep machdep.cpu.features | grep -iq 'avx512'; then
		echo -e "The CPU supports the AVX512 SIMD build mode.\n"
		ARGS+=(-DUSE_AVX512 -march=native)
	elif sysctl -a | grep machdep.cpu.features | grep -iq 'avx2'; then
		echo -e "The CPU supports the AVX2 SIMD build mode.\n"
		ARGS+=(-DUSE_AVX2 -march=native -mavx2)
	elif sysctl -a | grep machdep.cpu.features | grep -iq 'avx'; then
		echo -e "The CPU supports the AVX SIMD build mode.\n"
		ARGS+=(-DUSE_AVX -march=native -mavx)
	elif sysctl -a | grep machdep.cpu.features | grep -iq 'sse2'; then
		echo -e "The CPU supports the SSE2 SIMD build mode.\n"
		# On my Core2Duo Mac, 'native' gives "error: bad value for -march= switch":
		ARGS+=(-DUSE_SSE2 -march=core2)
	elif sysctl -a | grep machdep.cpu.features | grep -iq 'asimd'; then
		echo -e "The CPU supports the ASIMD build mode.\n"
		ARGS+=(-DUSE_ARM_V8_SIMD -march=native)
	else
		echo -e "The CPU supports no Mlucas-recognized ASIMD build mode ... building in scalar-double mode.\n"
		ARGS+=(-march=native)
	fi

else

	# Linux:
	if grep -iq 'avx512' /proc/cpuinfo; then
		echo -e "The CPU supports the AVX512 SIMD build mode.\n"
		ARGS+=(-DUSE_AVX512 -march=native)
	elif grep -iq 'avx2' /proc/cpuinfo; then
		echo -e "The CPU supports the AVX2 SIMD build mode.\n"
		ARGS+=(-DUSE_AVX2 -march=native -mavx2)
	elif grep -iq 'avx' /proc/cpuinfo; then
		echo -e "The CPU supports the AVX SIMD build mode.\n"
		ARGS+=(-DUSE_AVX -march=native -mavx)
	elif grep -iq 'sse2' /proc/cpuinfo; then
		echo -e "The CPU supports the SSE2 SIMD build mode.\n"
		ARGS+=(-DUSE_SSE2 -march=native)
	elif grep -iq 'asimd' /proc/cpuinfo; then
		echo -e "The CPU supports the ASIMD build mode.\n"
		ARGS+=(-DUSE_ARM_V8_SIMD -march=native)
	else
		echo -e "The CPU supports no Mlucas-recognized ASIMD build mode ... building in scalar-double mode.\n"
		ARGS+=(-march=native)
	fi

fi

if [[ -d $DIR ]]; then
	echo "Warning: '$DIR' already exists"
fi

# -p prevents "File exists" warning if obj-dir already exists:
mkdir -p "$DIR"
cd "$DIR"

if [[ -x $TARGET ]]; then
	echo "Error: '$DIR/$TARGET' already exists." >&2
	exit 1
fi

# Clang-under-MacOS linker barfs if one tries to explicitly invoke standard libs - h/t tdulcet for the
# conditional-inline syntax. Some OSes put the GMP headers in /usr/local/include, so -I that path in the
# compile command. If said path does not exist, make silently ignores it.
# Re. the -g flag to include the debugging symbols, they bloat executable size but if someone's Mlucas
# crashes/segfaults, one can rerun with GDB (gdb -ex=r ./Mlucas) to see the filename, line number and
# stack trace of the issue. If one wishes, one can run 'strip -g Mlucas' to remove the debugging symbols:
cat <<EOF >Makefile
CC?=gcc
CFLAGS=-fdiagnostics-color -Wall -g -O3 # -flto=auto
CPPFLAGS=-I/usr/local/include
LDLIBS=$(echo "$OSTYPE" | grep -iq 'darwin' || echo "-lm -lpthread -lrt") ${LARG[@]}

OBJS=br.o dft_macro.o fermat_mod_square.o fgt_m61.o get_cpuid.o get_fft_radices.o get_fp_rnd_const.o get_preferred_fft_radix.o getRealTime.o imul_macro.o mers_mod_square.o mi64.o Mlucas.o pairFFT_mul.o pair_square.o pm1.o qfloat.o radix1008_ditN_cy_dif1.o radix1024_ditN_cy_dif1.o radix104_ditN_cy_dif1.o radix10_ditN_cy_dif1.o radix112_ditN_cy_dif1.o radix11_ditN_cy_dif1.o radix120_ditN_cy_dif1.o radix128_ditN_cy_dif1.o radix12_ditN_cy_dif1.o radix13_ditN_cy_dif1.o radix144_ditN_cy_dif1.o radix14_ditN_cy_dif1.o radix15_ditN_cy_dif1.o radix160_ditN_cy_dif1.o radix16_dif_dit_pass.o radix16_ditN_cy_dif1.o radix16_dyadic_square.o radix16_pairFFT_mul.o radix16_wrapper_ini.o radix16_wrapper_square.o radix176_ditN_cy_dif1.o radix17_ditN_cy_dif1.o radix18_ditN_cy_dif1.o radix192_ditN_cy_dif1.o radix208_ditN_cy_dif1.o radix20_ditN_cy_dif1.o radix224_ditN_cy_dif1.o radix22_ditN_cy_dif1.o radix240_ditN_cy_dif1.o radix24_ditN_cy_dif1.o radix256_ditN_cy_dif1.o radix26_ditN_cy_dif1.o radix288_ditN_cy_dif1.o radix28_ditN_cy_dif1.o radix30_ditN_cy_dif1.o radix31_ditN_cy_dif1.o radix320_ditN_cy_dif1.o radix32_dif_dit_pass.o radix32_ditN_cy_dif1.o radix32_dyadic_square.o radix32_wrapper_ini.o radix32_wrapper_square.o radix352_ditN_cy_dif1.o radix36_ditN_cy_dif1.o radix384_ditN_cy_dif1.o radix4032_ditN_cy_dif1.o radix40_ditN_cy_dif1.o radix44_ditN_cy_dif1.o radix48_ditN_cy_dif1.o radix512_ditN_cy_dif1.o radix52_ditN_cy_dif1.o radix56_ditN_cy_dif1.o radix5_ditN_cy_dif1.o radix60_ditN_cy_dif1.o radix63_ditN_cy_dif1.o radix64_ditN_cy_dif1.o radix6_ditN_cy_dif1.o radix72_ditN_cy_dif1.o radix768_ditN_cy_dif1.o radix7_ditN_cy_dif1.o radix80_ditN_cy_dif1.o radix88_ditN_cy_dif1.o radix8_dif_dit_pass.o radix8_ditN_cy_dif1.o radix960_ditN_cy_dif1.o radix96_ditN_cy_dif1.o radix992_ditN_cy_dif1.o radix9_ditN_cy_dif1.o rng_isaac.o threadpool.o twopmodq100.o twopmodq128_96.o twopmodq128.o twopmodq160.o twopmodq192.o twopmodq256.o twopmodq64_test.o twopmodq80.o twopmodq96.o twopmodq.o types.o util.o
OBJS_MFAC=getRealTime.o get_cpuid.o get_fft_radices.o get_fp_rnd_const.o imul_macro.o mi64.o qfloat.o rng_isaac.o twopmodq100.o twopmodq128_96.o twopmodq128.o twopmodq160.o twopmodq192.o twopmodq256.o twopmodq64_test.o twopmodq80.o twopmodq96.o twopmodq.o types.o util.o threadpool.o factor.o

$Mlucas: \$(OBJS)
	\$(CC) \$(LDFLAGS) \$(CFLAGS) -o \$@ \$^ \$(LDLIBS)
$Mfactor: \$(OBJS_MFAC)
	\$(CC) \$(LDFLAGS) \$(CFLAGS) -o \$@ \$^ \$(LDLIBS)
factor.o: ../src/factor.c
	\$(CC) \$(CFLAGS) \$(CPPFLAGS) -c ${ARGS[@]} -DFACTOR_STANDALONE $WORDS -DTRYQ=4 \$<
%.o: ../src/%.c
	\$(CC) \$(CFLAGS) \$(CPPFLAGS) -c ${ARGS[@]} \$<
clean:
	rm -f \$(OBJS) \$(OBJS_MFAC)

.phony: clean
EOF

# if [[ -e build.log ]]; then
	# cp -vf --backup=t build.log{,}
# fi

echo -e "Building $TARGET"
printf "%'d CPU cores detected ... parallel-building using that number of make threads.\n" "$CPU_THREADS"
if ! time make -j "$CPU_THREADS" "$TARGET" >build.log 2>&1; then
	echo -e "\n*** There were build errors - see '${DIR}/build.log' for details. ***\n" >&2
	grep -A 2 'error:' build.log || tail build.log
	# exit 1
fi

echo -e "\nWarnings:\n"
grep 'warning:' build.log | awk '{ print $NF }' | sort | uniq -c | sort -nr

echo -e "\nErrors:\n"
grep -A 2 'error:' build.log || echo "None"

echo
