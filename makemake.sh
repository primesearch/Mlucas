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

DIR=obj
Mlucas=Mlucas
Mfactor=Mfactor
TARGET=$Mlucas
ARGS=(-DUSE_THREADS) # Optional compile args
WORDS=''
# Optional link args
LD_ARGS=()
# Optional Make args
MAKE_ARGS=()

MODES=()
GMP=1
HWLOC=0

case $OSTYPE in
	darwin*)
		echo -e "MacOS detected for build host.\n"
		CPU_THREADS=$(sysctl -n hw.ncpu)
		;;
	msys | cygwin)
		echo -e "Windows detected for build host.\n"
		CPU_THREADS=$NUMBER_OF_PROCESSORS
		;;
	linux* | *)
		echo -e "Assuming OS = Linux for build host.\n"
		CPU_THREADS=$(nproc --all)
		;;
esac

MAKE=make
if ! command -v $MAKE >/dev/null && command -v mingw32-make >/dev/null; then
	MAKE=mingw32-make
fi
if ! command -v $MAKE >/dev/null; then
	echo "Error: This script requires Make" >&2
	echo "On Ubuntu and Debian run: 'sudo apt update' and 'sudo apt install -y build-essential'" >&2
	exit 1
fi
if [[ -n $CC ]]; then
	if ! command -v "$CC" >/dev/null; then
		echo "Error: $CC is not installed." >&2
		exit 1
	fi
elif ! command -v gcc >/dev/null; then
	echo "Error: This script requires the GNU C compiler" >&2
	echo "On Ubuntu and Debian run: 'sudo apt update' and 'sudo apt install -y build-essential'" >&2
	exit 1
fi

# Returns success iff $CC (default gcc) accepts the given flag(s) for a full compile-and-link of a
# trivial program - used below to auto-detect toolchain-version-dependent flag/feature support instead
# of hardcoding version-number cutoffs (which drift out of date and vary by distro/backport):
try_flag() {
	printf 'int main(void){return 0;}\n' | "${CC:-gcc}" "$@" -x c -o /dev/null - >/dev/null 2>&1
}

# Returns success iff $CC's assembler accepts the AVX-512 constructs Mlucas's inline-asm kernels use:
# the "extended" register names (zmm16-31, xmm16-31, k0-k7), and the AVX-512ER reciprocal instruction
# vrcp28pd (mi64_modmul53_batch, compiled and called under plain USE_AVX512). Some older Clang releases
# reject one or both even when otherwise AVX-512-aware - e.g. clang 3.8 rejects the extended register
# names, clang 5.0 assembles those but rejects vrcp28pd - so probe both, exactly as the CI does (#73):
try_avx512_asm() {
	"${CC:-gcc}" -mavx512f -x c -o /dev/null - >/dev/null 2>&1 <<'EOF'
int main(void)
{
	__asm__ __volatile__(
		"vpxord %%zmm31,%%zmm31,%%zmm31\n\t"
		"vmovdqa64 %%xmm16,%%xmm17\n\t"
		"kmovw %%k1,%%eax\n\t"
		"vrcp28pd %%zmm8,%%zmm0"
		::: "zmm31","xmm16","xmm17","k1","eax","zmm0","zmm8"
	);
	return 0;
}
EOF
}

# Returns success iff $CC can compile two separate translation units with the given -flto[=...] variant
# and link them together. A single-file compile-and-link (like try_flag) is too weak a test here: LTO
# breakage on some older or misconfigured toolchains (missing/mismatched ar/nm/ranlib plugin support,
# old binutils) only shows up once the linker actually has to combine LTO object files from more than
# one translation unit - which is exactly what building Mlucas's ~90 source files does:
try_lto() {
	local cc=${CC:-gcc} flag=${1:--flto} tmpdir
	tmpdir=$(mktemp -d) || return 1
	trap 'rm -rf "$tmpdir"' RETURN
	printf 'int mm_lto_probe_helper(void){return 0;}\n' >"$tmpdir/a.c"
	printf 'int mm_lto_probe_helper(void);\nint main(void){return mm_lto_probe_helper();}\n' >"$tmpdir/b.c"
	(
		cd "$tmpdir" && \
		"$cc" "$flag" -c a.c -o a.o && \
		"$cc" "$flag" -c b.c -o b.o && \
		"$cc" "$flag" a.o b.o -o out
	) >/dev/null 2>&1
}

# GNU Make's -O (synchronize parallel-job output) flag needs Make >= 4.0 - probe for the flag itself
# rather than assuming by version number (which drifts, and varies by distro/backport):
MAKE_SUPPORTS_dashO=0
"$MAKE" --help 2>/dev/null | grep -wq -- '-O' && MAKE_SUPPORTS_dashO=1

if [[ ! $OSTYPE == darwin* ]]; then
	LD_ARGS+=(-lm -lpthread)
	if [[ $OSTYPE != msys && $OSTYPE != cygwin ]]; then
		LD_ARGS+=(-lrt)
	fi
fi
((MAKE_SUPPORTS_dashO)) && MAKE_ARGS+=(-O)
MAKE_ARGS+=(-j "$CPU_THREADS")

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
		no_gmp)
			GMP=0
			;;
		use_hwloc)
			HWLOC=1
			;;
		avx512_skylake | avx512_knl | avx512 | k1om | avx2 | avx | sse2 | asimd | nosimd)
			MODES+=("$arg")
			;;
		mfac)
			TARGET=$Mfactor
			;;
		[1-4n]word)
			WORDS=$arg
			;;
		*)
			echo "Usage: $0 [SIMD build mode]" >&2
			echo "Optional arguments must be 'no_gmp', 'use_hwloc' or one and only one of the supported SIMD-arithmetic types:" >&2
			echo -e "\t[x86_64: avx512 k1om avx2 avx sse2]; [Armv8: asimd]; or 'nosimd' for scalar-double build.\n" >&2
			exit 1
			;;
	esac

done

if ((GMP)); then
	LD_ARGS+=(-lgmp)
else
	echo "Building sans Gnu-MP ... this means no GCDs will be taken in p-1 work."
	ARGS+=(-DINCLUDE_GMP=0)
fi

if ((HWLOC)); then
	echo "Building with HWLOC hardware-topology support."
	ARGS+=(-DINCLUDE_HWLOC=1)
	LD_ARGS+=(-lhwloc)
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

if [[ $OSTYPE == msys || $OSTYPE == cygwin ]]; then
	Mlucas+=.exe
	Mfactor+=.exe
	TARGET+=.exe
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
		avx512_skylake)
			echo "Building for avx512_skylake SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			echo "Warning: The 'avx512_skylake' option is deprecated, use 'avx512' instead."
			ARGS+=(-DUSE_AVX512 -march=skylake-avx512 -mavx512f -mavx512cd -mavx512dq -mavx512bw -mavx512vl -mfma)
			;;
		avx512_knl)
			echo "Building for avx512_knl SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			echo "Warning: The 'avx512_knl' option is deprecated, use 'avx512' instead."
			ARGS+=(-DUSE_AVX512 -march=knl -mavx512f -mavx512cd -mavx512er -mfma)
			;;
		avx512)
			echo "Building for AVX512 SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_AVX512 -mavx512f -mavx512cd -mavx512dq -mavx512bw -mavx512vl -mfma)
			;;
		k1om)
			echo "Building for 1st-gen Xeon Phi 512-bit SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_IMCI512)
			;;
		avx2)
			echo "Building for AVX2 SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_AVX2 -mavx2 -mfma)
			;;
		avx)
			echo "Building for AVX SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_AVX -mavx)
			;;
		sse2)
			echo "Building for SSE2 SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_SSE2 -msse2)
			;;
		asimd)
			echo "Building for ASIMD SIMD in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			ARGS+=(-DUSE_ARM_V8_SIMD)
			;;
		nosimd)
			echo "Building in scalar-double (no-SIMD) mode in directory '${DIR}_${arg}'; the executable will be named '${TARGET}'"
			# This one's a no-op
			;;
		*)
			echo "Unrecognized SIMD-build flag ... aborting." >&2
			exit 1
			;;
	esac

	# The AVX-512 kernels use "extended" register names (zmm16-31/xmm16-31/k0-7) that some older
	# Clang releases reject even when otherwise AVX-512-aware; probe rather than let the user hit a
	# wall of asm errors deep in the build:
	if [[ $arg == avx512* ]] && ! try_avx512_asm; then
		echo "Error: ${CC:-gcc}'s assembler does not support the AVX-512 extended register names (zmm16-31/xmm16-31/k0-7) needed for this build mode ... aborting. Try a newer compiler, or build with 'avx2' instead." >&2
		exit 1
	fi

	DIR+="_$arg"

elif [[ $OSTYPE == darwin* ]]; then

	# MacOS:
	if (($(sysctl -n hw.optional.avx512f))) && try_avx512_asm; then
		echo -e "The CPU supports the AVX512 SIMD build mode.\n"
		ARGS+=(-DUSE_AVX512 -march=native -mavx512f -mavx512cd -mavx512dq -mavx512bw -mavx512vl -mfma)
	elif (($(sysctl -n hw.optional.avx512f))) || (($(sysctl -n hw.optional.avx2_0))); then
		if (($(sysctl -n hw.optional.avx512f))); then
			echo "Warning: CPU supports AVX-512 but ${CC:-gcc}'s assembler rejects the extended register names needed ... falling back to AVX2." >&2
		fi
		echo -e "The CPU supports the AVX2 SIMD build mode.\n"
		ARGS+=(-DUSE_AVX2 -march=native -mavx2 -mfma)
	elif (($(sysctl -n hw.optional.avx1_0))); then
		echo -e "The CPU supports the AVX SIMD build mode.\n"
		ARGS+=(-DUSE_AVX -march=native -mavx)
	elif (($(sysctl -n hw.optional.sse2))); then
		echo -e "The CPU supports the SSE2 SIMD build mode.\n"
		# On my Core2Duo Mac, 'native' gives "error: bad value for -march= switch":
		ARGS+=(-DUSE_SSE2 -march=core2 -msse2)
	elif (($(sysctl -n hw.optional.neon))); then
		echo -e "The CPU supports the ASIMD build mode.\n"
		ARGS+=(-DUSE_ARM_V8_SIMD)
		if try_flag -mcpu=native; then
			ARGS+=(-mcpu=native)
		elif try_flag -march=native; then
			ARGS+=(-march=native)
		fi
		# else: no arch flag - aarch64 has NEON/ASIMD in its baseline ISA, and ancient clang (e.g. 3.8) supports
		# neither -mcpu=native nor -march=native, so building without either still yields a working ASIMD binary
	else
		echo -e "The CPU supports no Mlucas-recognized SIMD build mode ... building in scalar-double mode.\n"
		echo "Warning: If this is a 64-bit x86 or ARM system, this likely means there is a bug in this script. Please report!"
		ARGS+=(-march=native)
	fi

elif [[ $OSTYPE == linux* ]]; then

	# Linux:
	if grep -iq 'avx512' /proc/cpuinfo && try_avx512_asm; then
		echo -e "The CPU supports the AVX512 SIMD build mode.\n"
		ARGS+=(-DUSE_AVX512 -march=native -mavx512f -mavx512cd -mavx512dq -mavx512bw -mavx512vl -mfma)
	elif grep -iq 'avx512\|avx2' /proc/cpuinfo; then
		if grep -iq 'avx512' /proc/cpuinfo; then
			echo "Warning: CPU supports AVX-512 but ${CC:-gcc}'s assembler rejects the extended register names needed ... falling back to AVX2." >&2
		fi
		echo -e "The CPU supports the AVX2 SIMD build mode.\n"
		ARGS+=(-DUSE_AVX2 -march=native -mavx2 -mfma)
	elif grep -iq 'avx' /proc/cpuinfo; then
		echo -e "The CPU supports the AVX SIMD build mode.\n"
		ARGS+=(-DUSE_AVX -march=native -mavx)
	elif grep -iq 'sse2' /proc/cpuinfo; then
		echo -e "The CPU supports the SSE2 SIMD build mode.\n"
		ARGS+=(-DUSE_SSE2 -march=native -msse2)
	elif grep -iq 'asimd' /proc/cpuinfo && [[ $HOSTTYPE == aarch64 ]]; then
		echo -e "The CPU supports the ASIMD build mode.\n"
		ARGS+=(-DUSE_ARM_V8_SIMD)
		if try_flag -mcpu=native; then
			ARGS+=(-mcpu=native)
		elif try_flag -march=native; then
			ARGS+=(-march=native)
		fi
		# else: no arch flag - aarch64 has NEON/ASIMD in its baseline ISA, and ancient clang (e.g. 3.8) supports
		# neither -mcpu=native nor -march=native, so building without either still yields a working ASIMD binary
	else
		echo -e "The CPU supports no Mlucas-recognized SIMD build mode ... building in scalar-double mode.\n"
		echo "Warning: If this is a 64-bit x86 or ARM system, this likely means there is a bug in this script. Please report!"
		ARGS+=(-march=native)
	fi

else

	# Fallback path for hosts without /proc/cpuinfo or sysctl (notably Windows/MSYS2/Cygwin): compile a tiny
	# probe that just reports the CPU's highest SIMD level as a keyword, then map that to build flags in the
	# shell below - reusing the same try_avx512_asm / try_flag probes as the Linux and Darwin branches so the
	# AVX-512 extended-register-name check and the -mcpu/-march fallback apply here too (see #60, #67).
	# Adapted from: https://stackoverflow.com/a/28939692
	tmpdir=$(mktemp -d)
	trap 'rm -rf "$tmpdir"' EXIT
	cat <<'EOF' >"$tmpdir/simd.c"
#include <stdio.h>
int main()
{
// defined(__amd64) || defined(__amd64__) || defined(_M_AMD64) || defined(_M_EMT64) || defined(__x86_64) || defined(__x86_64__)
#ifdef __x86_64__
	#ifdef __AVX512F__
		puts("avx512");
	#elif defined __AVX2__
		puts("avx2");
	#elif defined __AVX__
		puts("avx");
	#elif defined __SSE2__
		puts("sse2");
	#else
		puts("none_x86");
	#endif
#elif defined(__aarch64__)
	#ifdef __ARM_NEON
		puts("asimd");
	#else
		puts("none_arm");
	#endif
#else
	puts("none");
#endif
	return 0;
}
EOF

	args=()
	case $HOSTTYPE in
		aarch64 | arm*)
			args+=(-mcpu=native)
			;;
		x86_64 | *)
			args+=(-march=native)
			;;
	esac
	"${CC:-gcc}" -Wall -g -O3 "${args[@]}" -o "$tmpdir/simd" "$tmpdir/simd.c"
	if ! output=$("$tmpdir/simd"); then
		echo "Error: Unable to detect the SIMD build mode" >&2
		exit 1
	fi

	case $output in
		avx512)
			if try_avx512_asm; then
				echo -e "The CPU supports the AVX512 SIMD build mode.\n"
				ARGS+=(-DUSE_AVX512 -march=native -mavx512f -mavx512cd -mavx512dq -mavx512bw -mavx512vl -mfma)
			else
				echo "Warning: CPU supports AVX-512 but ${CC:-gcc}'s assembler rejects the extended register names needed ... falling back to AVX2." >&2
				echo -e "The CPU supports the AVX2 SIMD build mode.\n"
				ARGS+=(-DUSE_AVX2 -march=native -mavx2 -mfma)
			fi
			;;
		avx2)
			echo -e "The CPU supports the AVX2 SIMD build mode.\n"
			ARGS+=(-DUSE_AVX2 -march=native -mavx2 -mfma)
			;;
		avx)
			echo -e "The CPU supports the AVX SIMD build mode.\n"
			ARGS+=(-DUSE_AVX -march=native -mavx)
			;;
		sse2)
			echo -e "The CPU supports the SSE2 SIMD build mode.\n"
			ARGS+=(-DUSE_SSE2 -march=native -msse2)
			;;
		asimd)
			echo -e "The CPU supports the ASIMD build mode.\n"
			ARGS+=(-DUSE_ARM_V8_SIMD)
			if try_flag -mcpu=native; then
				ARGS+=(-mcpu=native)
			elif try_flag -march=native; then
				ARGS+=(-march=native)
			fi
			# else: no arch flag - aarch64 has NEON/ASIMD in its baseline ISA, and ancient clang (e.g. 3.8) supports
			# neither -mcpu=native nor -march=native, so building without either still yields a working ASIMD binary
			;;
		none_arm)
			echo -e "The CPU supports no Mlucas-recognized SIMD build mode ... building in scalar-double mode.\n"
			echo "Warning: This likely means there is a bug in this script. Please report!" >&2
			if try_flag -mcpu=native; then
				ARGS+=(-mcpu=native)
			elif try_flag -march=native; then
				ARGS+=(-march=native)
			fi
			;;
		*)
			echo -e "The CPU supports no Mlucas-recognized SIMD build mode ... building in scalar-double mode.\n"
			echo "Warning: If this is a 64-bit x86 or ARM system, this likely means there is a bug in this script. Please report!" >&2
			ARGS+=(-march=native)
			;;
	esac
fi

if [[ -d $DIR ]]; then
	echo "Warning: Directory '$DIR' already exists"
fi

# -p prevents "File exists" warning if obj-dir already exists:
mkdir -p "$DIR"
cd "$DIR"

if [[ -x $TARGET ]]; then
	echo "Error: '$DIR/$TARGET' already exists." >&2
	exit 1
fi

# -fdiagnostics-color needs GCC >= 4.9 (or a recent-enough Clang); -flto is broken/absent on some older
# or misconfigured toolchains (notably some Clang-on-old-glibc and MSYS2-Clang combos) - probe for both
# instead of assuming. CI jobs that need a different CFLAGS entirely (sanitizer builds) should export a
# CFLAGS environment variable before invoking this script - the generated Makefile's "CFLAGS ?=" already
# defers to a pre-set environment CFLAGS instead of the computed value below. Prefer -flto=auto (parallel
# LTO codegen, see #56) over plain -flto when supported:
CFLAGS_PROBED=(-Wall -g -O3)
try_flag -fdiagnostics-color && CFLAGS_PROBED=(-fdiagnostics-color "${CFLAGS_PROBED[@]}")
if try_lto -flto=auto; then
	CFLAGS_PROBED+=(-flto=auto)
elif try_lto -flto; then
	CFLAGS_PROBED+=(-flto)
else
	echo "Warning: ${CC:-gcc} does not support (or reliably link with) -flto ... building without LTO." >&2
fi

# Clang-under-MacOS linker barfs if one tries to explicitly invoke standard libs - h/t tdulcet for the
# conditional-inline syntax. Some OSes put the GMP headers in /usr/local/include, so -I that path in the
# compile command. If said path does not exist, make silently ignores it.
# Re. the -g flag to include the debugging symbols, they bloat executable size but if someone's Mlucas
# crashes/segfaults, one can rerun with GDB (gdb -ex=r ./Mlucas) to see the filename, line number and
# stack trace of the issue. If one wishes, one can run 'strip -g Mlucas' to remove the debugging symbols:
cat <<EOF >Makefile
CC ?= gcc
CFLAGS ?= ${CFLAGS_PROBED[*]}
CPPFLAGS ?= -D_GNU_SOURCE -I/usr/local/include -I/opt/homebrew/include
LDFLAGS ?= -L/opt/homebrew/lib
LDLIBS ?= ${LD_ARGS[@]} # -static

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

echo -e "Building $TARGET"
printf "%'d CPU cores detected ... parallel-building using that number of make threads.\n" "$CPU_THREADS"
if ! time $MAKE "${MAKE_ARGS[@]}" "$TARGET" >build.log 2>&1; then
	echo -e "\n*** There were build errors - see '${DIR}/build.log' for details. ***\n" >&2
	grep -A 2 '[Ee]rror:' build.log || tail build.log
	exit 1
fi

echo -e "\nWarnings:\n"
grep 'warning:' build.log | grep '#warning' | sed 's/^.*warning://p' | sort -u
echo -e "\nWarning counts:\n"
grep 'warning:' build.log | awk '{ print $NF }' | sort | uniq -c | sort -nr || echo "None"

echo -e "\nErrors:\n"
grep -A 2 '[Ee]rror:' build.log || echo "None"

echo
