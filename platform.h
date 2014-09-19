/*******************************************************************************
*                                                                              *
*   (C) 1997-2013 by Ernst W. Mayer.                                           *
*                                                                              *
*  This program is free software; you can redistribute it and/or modify it     *
*  under the terms of the GNU General Public License as published by the       *
*  Free Software Foundation; either version 2 of the License, or (at your      *
*  option) any later version.                                                  *
*                                                                              *
*  This program is distributed in the hope that it will be useful, but WITHOUT *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   *
*  more details.                                                               *
*                                                                              *
*  You should have received a copy of the GNU General Public License along     *
*  with this program; see the file GPL.txt.  If not, you may view one at       *
*  http://www.fsf.org/licenses/licenses.html, or obtain one by writing to the  *
*  Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA     *
*  02111-1307, USA.                                                            *
*                                                                              *
*******************************************************************************/

/****************************************************************************
 * Principal platform (CPU/compiler/OS) identification header file.
 ****************************************************************************/
#ifndef platform_h_included
#define platform_h_included

#undef	PLATFORM_DEBUG
#undef	OS_DEBUG
#undef	OS_BITS_DEBUG
#undef	CPU_DEBUG
#undef	CMPLR_DEBUG

// Make PLATFORM_DEBUG def'd the default for CUDA builds, since those use 2 distinct passes which it is useful to clearly separate:
#ifdef __CUDACC__
	#define PLATFORM_DEBUG
#endif

// Define PLATFORM_DEBUG at compile time to enable these internal platfrom-related diagnostic #warning prints:
#ifdef PLATFORM_DEBUG
	/* Set = 1 to print brief OS summary at compile time: */
	#define	OS_DEBUG	1
	/* Set = 1 to print brief OS-bitness summary at compile time: */
	#define	OS_BITS_DEBUG	1
	/* Set = 1 to print brief OS summary at compile time: */
	#define	CPU_DEBUG	1
	/* Set = 1 to print brief OS summary at compile time: */
	#define	CMPLR_DEBUG	1
#endif

/* Platform-dependent #defines summarizing key performance-related features: */
#undef	FP_MANTISSA_BITS_DOUBLE	/* Number of significand bits of a register double.
									Typically = 53 for IEEE-compliant RISC,
									            64 for x86-style FPU.
								*/
#undef	MUL64_SUCKS	/* 64-bit integer MUL available but crappy? (E.g. Sparc) */
#undef	MULH64_FAST	/* Fast hardware support for upper half (or both halves simultaneously)
					of 64x64=>128-bit unsigned integer product?
					This implies that MUL64x32 aliases to MUL_LOHI, i.e. that the
					2nd argument to a MUL64x32 call need not be restricted to 32 bits. */

#undef	HARDWARE_FMADD	/* Indicates whether the hardware in question supports a fused multiply-add operation. */

#undef	USE_FMADD	/* Only set if HARDWARE_FMADD is set *and* user has #def'd FMADD_YES.
					If set, user is expected to provide a suitable associated
					set of macros in the file float_intrin.h */

/* Locally defined OS types: */
#undef	OS_NAME

#undef	OS_TYPE
#undef	OS_TYPE_UNKNOWN
#undef	OS_TYPE_WINDOWS
#undef	OS_TYPE_LINUX
#undef	OS_TYPE_MACOSX
#undef	OS_TYPE_DECOSF		/* DEC OSF (later HP TruUnix) */
#undef	OS_TYPE_DECVMS		/* DEC VMS (originally for VAX) */
#undef	OS_TYPE_HPUX
#undef	OS_TYPE_SUN
#undef	OS_TYPE_AIX
#undef	OS_TYPE_IRIX

#undef	OS_VERSION

/* Some of these are currently identified via the compiler (below) rather than directly: */
#if(defined(WINDOWS) || defined(_WINDOWS) || defined(WIN32) || defined(_WIN32))
	#define	OS_TYPE
	#define	OS_TYPE_WINDOWS
#elif(defined(linux) || defined(__linux__) || defined(__linux))
	#define	OS_TYPE
	#define	OS_TYPE_LINUX
#elif(defined(__APPLE__))
	#define	OS_TYPE
	#define	OS_TYPE_MACOSX
#elif(defined(_AIX))
	#define	OS_TYPE
	#define	OS_TYPE_AIX
#elif(defined(__osf__))
	#define OS_TYPE
	#define	OS_TYPE_DECOSF
#elif(defined(__VMS))
	#define	OS_TYPE
	#define	OS_TYPE_DECVMS
#elif(defined(sun))	/* 20 Nov 2006: removed " && defined(sparc)" from here, since Solaris now runs on other platforms, e.g. AMD64 */
	#define	OS_TYPE
	#define	OS_TYPE_SUN
#else
	#error Unable to determine OS type!
#endif

// Posix-compliant OS?
#undef OS_POSIX_COMPLIANT
#ifndef OS_TYPE_WINDOWS
	#define OS_POSIX_COMPLIANT
#endif

#undef	OS_BITS

#undef MULTITHREAD	/* Master #define controlling multithreaded mode -- enable by invoking USE_THREADS
					and the appropriate compiler options (e.g. -openmp; compiler-dependent) at build time. */

/*
Locally defined CPU types - undef these here to prevent namespace collisions.
Each of the above CPUs typically comes in multiple flavors: try to restrict
our list here to the salient ones, that actually need to be differentiated
based on different key capabilities. Default CPU subtypes are as indicated.
*/
#undef	CPU_TYPE
#undef	CPU_NAME
#undef	CPU_SUBTYPE_NAME

#undef	CPU_IS_UNKNOWN
	#undef	CPU_SUBTYPE_UNKNOWN

#undef	CPU_IS_ALFA
	#undef	CPU_SUBTYPE_EV4
	#undef	CPU_SUBTYPE_EV5
	#undef	CPU_SUBTYPE_EV6

#undef	CPU_IS_MIPS
	#undef	CPU_SUBTYPE_PRE_R10K
	#undef	CPU_SUBTYPE_R10K
	#undef	CPU_SUBTYPE_R12K
	#undef	CPU_SUBTYPE_R15K

#undef	CPU_IS_SPARC
	#undef	CPU_SUBTYPE_ULTRA1
	#undef	CPU_SUBTYPE_ULTRA2
	#undef	CPU_SUBTYPE_ULTRA3

#undef	CPU_IS_PPC
	#undef	CPU_SUBTYPE_PPC32
	#undef	CPU_SUBTYPE_PPC64	/* This is the main differentiator we need for now,
								i.e. 64-bit or not? Separate check of AltiVec-or-not
								differentiates between G3 and G4 (both of which are PPC32). */
#undef	CPU_IS_X86
	#undef	CPU_SUBTYPE_CLASSIC
	#undef	CPU_SUBTYPE_PPRO
	#undef	CPU_SUBTYPE_P2
	#undef	CPU_SUBTYPE_P3
	#undef	CPU_SUBTYPE_P4

#undef	CPU_IS_X86_64
	#undef	CPU_SUBTYPE_AMD32
	#undef	CPU_SUBTYPE_AMD64
/* DO WE NEED ANY OF THESE?
	#undef	CPU_SUBTYPE_K7
	#undef	CPU_SUBTYPE_K8
	#undef	CPU_SUBTYPE_OPTERON
*/

#undef	CPU_IS_IA64	// Itanium
	#undef	CPU_SUBTYPE_IT1
	#undef	CPU_SUBTYPE_IT2

#undef	CPU_IS_POWER
	#undef	CPU_SUBTYPE_POWER1
	#undef	CPU_SUBTYPE_POWER2
	#undef	CPU_SUBTYPE_POWER3
	#undef	CPU_SUBTYPE_POWER4
	#undef	CPU_SUBTYPE_POWER5

#undef	CPU_IS_HPPA
	#undef	CPU_SUBTYPE_PA1
	#undef	CPU_SUBTYPE_PA2

// SIMD-functionality-related flags: We currently only care about the x86 family SSE2 / AVX / AVX2 instruction sets.
// We define these in a "grandfathered-in" fashion, i.e. each higher-functionality flag (= higher up in our descending-ordering)
// will automatically turn on the lower ones, e.g. AVX2 enables both AVX and SSE2. This is because we code things such that
// SIMD-enabled #ifdefs wrap variable defs shared by all these specific SIMD ISAs, but e.g. for AVX/AVX2 the vec_dbl type
// encodes a 4-double struct, whereas for SSE2, vec_dbl is a 2-double struct. In places where we really need to invoke
// separate code sections based on these #defs, we check them in inverse order, i.e. highest first. That way if e.g.
// AVX2, AVX and SSE2 have different versions of a common-named arithmetic macro, the fact that all > 1 of these #defines
// are enabled yields no conflict, since the highest-level one wins, as it were:

// We don't allow the user to set more than one of these defines:
#ifdef USE_AVX2
	#if defined(USE_AVX) || defined(USE_SSE2)
		#error Only one of USE_SSE2 / USE_AVX / USE_AVX2 may be defined at compile time!
	#endif

	#define USE_AVX
	#define	USE_SSE2

#elif defined(USE_AVX)
	#if defined(USE_AVX2) || defined(USE_SSE2)
		#error Only one of USE_SSE2 / USE_AVX / USE_AVX2 may be defined at compile time!
	#endif

	#define	USE_SSE2

#elif defined(USE_SSE2)
	#if defined(USE_AVX2) || defined(USE_AVX)
		#error Only one of USE_SSE2 / USE_AVX / USE_AVX2 may be defined at compile time!
	#endif
#endif

// GPU-specific flagging is done analogously to SIMD: REALLY_GPU implies USE_GPU, but not v.v.:
#ifdef REALLY_GPU
	#define	USE_GPU
#endif


/* Locally defined Compiler types and versions: */
#undef	COMPILER_NAME

#undef	COMPILER_TYPE
#undef	COMPILER_TYPE_UNKNOWN

#undef	COMPILER_VERSION

#undef	COMPILER_TYPE_APPLEC	/* Apple C compiler */
#undef	COMPILER_TYPE_DECC		/* DEC C compiler */
#undef	COMPILER_TYPE_GCC		/* Gnu C compiler */
#undef	COMPILER_TYPE_HPC		/* HP C compiler */
#undef	COMPILER_TYPE_ICC		/* Intel compiler */
#undef	COMPILER_TYPE_MSVC		/* Microsoft Visual C++ (later .NET) compiler */
#undef	COMPILER_TYPE_MWERKS	/* Metrowerks Codewarrior */
#undef	COMPILER_TYPE_NVCC		/* nVidia C++ (for GPU) compiler */
#undef	COMPILER_TYPE_SUNC		/* SunPro C compiler */
#undef	COMPILER_TYPE_XLC		/* IBM XL C compiler */
/* Miscellaneous notes on the key-but-lesser-known features of the compilers:
	- MWERKS and MSVC/.NET share the same inline ASM syntax, as do GCC and ICC.
	-
*/

// For GPU builds via NVCC, must allow for the 2-compile-pass paradigm:
//	Pass 1 [host code] should pass through to same predefines-handling as a non-GPU build on the host arch would;
//	Pass 2 [device code] is handled here.
// We can differentiate between the 2 passes by using that __CUDACC__ is def'd for both, but __CUDA_ARCH__ only for Pass 2:
#if defined(__CUDACC__) && !defined(__CUDA_ARCH__)

	#warning CUDA: Host-code compile pass

#endif
#ifdef __CUDA_ARCH__

	#warning CUDA: Device-code compile pass

	#define CPU_TYPE
	#define CPU_IS_NVIDIA

	#define COMPILER_TYPE
	#define COMPILER_TYPE_NVCC

	#define	USE_RINT	// Ensures that float/doule round-to-nearest proceeds via the efficient cuda-stdlib rint() function

  // If GCC version predefines __SIZEOF_POINTER__, that is most reliable (the LONG_MAX-based test below failed under mingw64 because that defined LONG_MAX and LONG_LONG_MAX differently in 64-bit mode:
  #ifdef __SIZEOF_POINTER__

	#if __SIZEOF_POINTER__ == 4
		#define OS_BITS 32
	#elif __SIZEOF_POINTER__ == 8
		#define OS_BITS 64
	#else
		#error __SIZEOF_POINTER__ defined but returns unrecognized value! Use gcc -dM -E < /dev/null | grep __SIZEOF_POINTER__ to examine.
	#endif

  // Otherwise see if can use the value of __LONG_MAX__ in limits.h to quickly and portably determine whether it's 32-bit or 64-it OS:
  #else
	// Syntax here is GCC and SunStudio/MSVC, respectively:
	#if !defined(__LONG_MAX__) && !defined(LONG_MAX)
		#include <limits.h>
	#endif

	#if !defined(__LONG_MAX__) &&  defined(LONG_MAX)
		#define __LONG_MAX__  LONG_MAX
	#endif

	#ifdef __LONG_MAX__
		#if __LONG_MAX__ == 2147483647L
			#define OS_BITS 32
		#elif __LONG_MAX__ == 9223372036854775807L
			#define OS_BITS 64
		#else
			#error  __LONG_MAX__ defined but value unrecognized!
		#endif
	#else
		#error platform.h: failed to properly set OS_BITS for NVCC build!
	#endif

  #endif

/* 32-bit X86, Gnu C or Metrowerks CodeWarrior C compiler: */
#elif(defined(__i386) || defined(__i386__))

	#ifndef OS_BITS
		#define OS_BITS 32
	#endif

	#define CPU_TYPE
	#define CPU_IS_X86
	#define	INTEGER_MUL_32

	#ifndef COMPILER_TYPE	// NVCC may already have been defined
		/* Sun C compiler - needs '-xarch=amd64' compiler flag in order for __amd64 to be def'd: */
		#if(defined(__SUNPRO_C))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_SUNC
		/* Gnu C compiler */
		#elif(defined(__GNUC__) || defined(__GNUG__))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_GCC
		/* Metrowerks CodeWarrior C compiler */
		#elif(defined(__MWERKS__))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_MWERKS
		/* Intel C compiler: */
		#elif(defined(__INTEL_COMPILER))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_ICC

			#include <xmmintrin.h>	/* Principal header file for Streaming SIMD Extensions intrinsics - we need it for the _MM_HINT predefines used in prefetch.h */
		/* Unknown: */
		#else
			#error __i386__ defined for unexpected compiler type!

			#define COMPILER_TYPE
			#define COMPILER_TYPE_UNKNOWN
		#endif
	#endif

/* AMD64 ISA - includes "Intel 64" (formerly known as "EMT64", formerly formerly known as "IA32e", yada, yada): */
#elif(defined(__amd64) || defined(__amd64__) || defined(_M_AMD64) || defined(_M_EMT64) || defined(__x86_64) || defined(__x86_64__))

	#ifndef OS_BITS
		#define OS_BITS 64
	#endif

	#define CPU_TYPE
	#define CPU_IS_X86_64

	#define MULH64_FAST
	#define	HARDWARE_FMADD

	#ifndef COMPILER_TYPE	// NVCC may already have been defined
		/* Sun C compiler - needs '-xarch=amd64' compiler flag in order for __amd64 to be def'd: */
		#if(defined(__SUNPRO_C))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_SUNC
		// Gnu C compiler: Note this gets triggered as desired for llvm/clang, but for that we use non-GCC flags
		// to set COMPILER_NAME to reflect that while the compiler is gcc-compatible, it is not gcc:
		#elif(defined(__GNUC__) || defined(__GNUG__))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_GCC
		/* Metrowerks CodeWarrior C compiler */
		#elif(defined(__MWERKS__))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_MWERKS
		/* Intel C compiler: */
		#elif(defined(__INTEL_COMPILER))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_ICC
		/* MS Visual C/Studio/.NET/Whatever: */
		#elif(defined(_MSC_VER))
			#define COMPILER_TYPE
			#define COMPILER_TYPE_MSVC
		/* Unknown: */
		#else
			#define COMPILER_TYPE
			#define COMPILER_TYPE_UNKNOWN
		#endif
	#endif

/* 32-bit X86, Intel C or MSVC/.NET compiler */
#elif(defined(_M_IX86) || defined(_MSC_VER))

	#ifndef OS_BITS
		#define OS_BITS 32
	#endif

	#define CPU_TYPE
	#define CPU_IS_X86
	#define	INTEGER_MUL_32

	/* Intel C compiler: */
	#if(defined(__INTEL_COMPILER))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_ICC

		#include <xmmintrin.h>	/* Principal header file for Streaming SIMD Extensions intrinsics - we need it for the _MM_HINT predefines used in prefetch.h */
		#include <emmintrin.h>
		#ifndef _INCLUDED_EMM
			#error Failed to include <emmintrin.h> header file!
		#endif

		/*
		Various header files containing ICC intrinsics for IA32:
		we include individual ones of these as needed by the various
		functionality modules (e.g. prefetch.h and sse2.h):
		#include <dvec.h>					* SSE 2 intrinsics for Class Libraries												*
		#include <fvec.h>					* SSE intrinsics for Class Libraries												*
		#include <ia32intrin.h>				* Miscellaneous IA-32 intrinsics.													*
		#include <ivec.h>					* MMX™ instructions intrinsics for Class Libraries									*
		#include <mathf.h>					* Principal header file for legacy Intel Math Library								*
		#include <mathimf.h>				* Principal header file for current Intel Math Library								*
		#include <mmintrin.h>				* Intrinsics for MMX instructions													*
		#include <omp.h>					* Principal header file OpenMP*														*
		#include <pgouser.h>				* For use in the instrumentation compilation phase of profile-guided optimizations	*
		#include <pmmintrin.h>				* Principal header file for SSE3 intrinsics											*
		#include <sse2mmx.h>				* Emulated versions of 128-bit SSE2 intrinsics (for P3 and Itanium)					*
		#include <stdarg.h>					* Replacement header for standard stdarg.h											*
		#include <varargs.h>				* Replacement header for standard varargs.h											*
		*/

	/* MS Visual C/Studio/.NET/Whatever: */
	#elif(defined(_MSC_VER))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_MSVC
	/* Unknown: */
	#else
		#define COMPILER_TYPE
		#define COMPILER_TYPE_UNKNOWN
	#endif
	/* TODO: p3/mmx stuff here */

	/* Important Pentium CPU subtypes 600 == p3, 700 == p4, etc: */
	#if(_M_IX86 >= 700)
		#define	CPU_SUBTYPE
		#define CPU_SUBTYPE_P4
		#define	CPU_SUBTYPE_NAME "P4"
	#endif

/* Itanium: */
#elif(defined(__ia64) || defined(__ia64__))

	#ifndef OS_BITS
		#define OS_BITS 64
	#endif

	#define CPU_TYPE
	#define CPU_IS_IA64

	#define MULH64_FAST
	#define	HARDWARE_FMADD

	/* HP C or C++ compiler for HPUX: */
	#if(defined(__HP_cc) || defined(__HP_aCC))

		#define	OS_TYPE
		#define	OS_TYPE_HPUX
		#define COMPILER_TYPE
		#define COMPILER_TYPE_HPC

		#include <fenv.h>
		#include <machine/sys/inline.h>

	/* IA-64 (Itanium), Intel C compiler (ICC v7 or above).
	We need to put this before the IA64/GCC case because, bizarrely,
	ICC also #defines __GNUC__ and __GNUG__ on the Itanium even
	though gcc-style inline ASM is not supported on that platform
	(which is why we need to include the ICC intrinsics header files).
	*/
	#elif(defined(__INTEL_COMPILER))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_ICC

		#include <ia64intrin.h>	/* Needed for access to e.g. _m64_xmahu intrinsic */
		#include <xmmintrin.h>	/* Principal header file for Streaming SIMD Extensions intrinsics - we need it for the _MM_HINT predefines used in prefetch.h */

		/*
		Various header files containing ICC intrinsics for IA64:
		we include individual ones of these as needed by the various
		functionality modules (e.g. prefetch.h and sse2.h):
		#include <emmintrin.h>				* Principal header file for SSE2 intrinsics										*
		#include <fvec.h>					* SSE intrinsics for Class Libraries											*
		#include <ia64intrin.h>				* Miscellaneous IA-64 intrinsics.												*
		#include <ia64regs.h>				* Header file for registers on Itanium-based systems							*
		#include <ivec.h>					* MMX™ instructions intrinsics for Class Libraries								*
		#include <mathimf.h>				* Principal header file for current Intel Math Library							*
		#include <mmclass.h>				* Principal header files for MMX™ instructions with Class Libraries				*
		#include <mmintrin.h>				* Intrinsics for MMX instructions												*
		#include <omp.h>					* Principal header file OpenMP*													*
		#include <sse2mmx.h>				* Emulated versions of 128-bit SSE2 intrinsics (for P3 and Itanium)					*
		*/

	/* IA-64 (Itanium), Gnu C or C++ compiler */
	#elif(defined(__GNUC__) || defined(__GNUG__))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_GCC
	#else
		#define COMPILER_TYPE
		#define COMPILER_TYPE_UNKNOWN
	#endif

/* HP PA-RISC UNDER HPUX */
#elif(defined(_PA_RISC1_0) || defined(_PA_RISC1_1) || defined(_PA_RISC2_0) || defined(hppa) || defined(__hppa__))

	#define	OS_TYPE
	#define	OS_TYPE_HPUX
	#define CPU_TYPE
	#define CPU_IS_HPPA
	#define COMPILER_TYPE
	#define COMPILER_TYPE_HPC

/* IBM Power: Note that I found that only __xlc__ was properly defined on the power5 I used for my tests.
We Deliberately put this #elif ahead of the PowerPC one because the Power/AIX compiler also defines some PowerPC flags: */
#elif(defined(_POWER) || defined(_ARCH_PWR))

	#ifndef	OS_TYPE_AIX
		#error ARCH = POWER defined but not OS_TYPE_AIX!
	#endif

	#define CPU_TYPE
	#define CPU_IS_POWER

	#define	HARDWARE_FMADD

	/* IBM XLC compiler */
	#if(defined(xlc) || defined(__xlc__) || defined(__xlc) || defined(XLC) || defined(__XLC__) || defined(__XLC) || defined(__xlC__))

		#define COMPILER_TYPE
		#define COMPILER_TYPE_XLC

	#elif(defined(__GNUC__) || defined(__GNUG__))

		#define COMPILER_TYPE
		#define COMPILER_TYPE_GCC

	#else
		#error unknown compiler for Power.
	#endif

/* PowerPC: */
/* TODO:***need way of distinguishing e.g. G3,4,5***/
/*
gcc for 32-bit PowerPC defines __ppc__, while gcc in 64-bit mode
defines __ppc64__ (and *not* __ppc__). IBM's xlc documentation
states that the compiler defines __64BIT__ if compiling in 64-bit mode.
*/
#elif(defined(__ppc__) || defined(__powerpc__) || defined(__PPC__) || defined(__powerc) || defined(__ppc64__))

	#define CPU_TYPE
	#define CPU_IS_PPC

	#define	HARDWARE_FMADD

/* TODO:***do we need separate ifdefs for gcc and applec?***/
	/* Gnu C compiler for OS-X: */
	#if(defined(__GNUC__) || defined(__GNUG__))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_GCC

		/* Get the hardware intrinsics file, usually at
		/usr/include/gcc/darwin/default/ppc_intrinsics.h .
		The XLC compiler includes the same intrinsics automatically.
		Under GCC and CW we need to define them ourselves. */
		#include <ppc_intrinsics.h>

	/* IBM XLC compiler */
	#elif(defined(xlc) || defined(__xlc__) || defined(__xlc) || defined(XLC) || defined(__XLC__) || defined(__XLC))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_XLC

	/* Metrowerks CodeWarrior C compiler */
	#elif(defined(__MWERKS__))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_MWERKS
		#include "ppc_intrinsics.h"

	/* Unknown: */
	#else
		#define COMPILER_TYPE
		#define COMPILER_TYPE_UNKNOWN
	#endif

	/* 32-bit vs. 64-bit: */
	#if(defined(__powerpc64__) || defined(__ppc64) || defined(__PPC64) || defined(__ppc64__) || defined(__PPC64__) || defined(__64BIT__))
		#ifndef OS_BITS
			#define OS_BITS 64
		#endif
		#define MULH64_FAST	/* "Fast" on the G5 is relative to 32x32-bit... */
		#define	CPU_SUBTYPE
		#define CPU_SUBTYPE_PPC64
		#define	CPU_SUBTYPE_NAME "64-bit"
	#else
		#ifndef OS_BITS
			#define OS_BITS 32
		#endif
		#define	CPU_SUBTYPE
		#define CPU_SUBTYPE_PPC32
		#define	CPU_SUBTYPE_NAME "32-bit"
		#define	INTEGER_MUL_32
	#endif

/* Alpha: */
#elif(defined(__alpha) || defined(__alpha__) || defined(__alpha) || defined(__ALPHA) || defined(__ALPHA__) || defined(__ALPHA))

	#ifndef OS_BITS
		#define OS_BITS 64
	#endif

	#define CPU_TYPE
	#define CPU_IS_ALFA

	#define MULH64_FAST

	/* Compaq (native) C compiler: */
	#if(defined(__DECC))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_DECC

		/* Alpha system under TruUNIX (formerly known as OSF) or VMS. */
		#if(defined(OS_TYPE_DECOSF) || defined(OS_TYPE_DECVMS))
			/* This assumes a link of form machine -> arch/alpha exists in /usr/include... */
			#include <machine/builtins.h>
			/* ...this doesn't. */
			#ifndef __BUILTINS_LOADED
				#include <arch/alpha/builtins.h>
			#endif
		/* Alpha system under Linux. Compile using ccc, not cc. */
		#elif(defined(OS_TYPE_LINUX))
			/* This assumes a link of form machine -> alpha exists in /usr/include... */
			#include <machine/builtins.h>
			/* ...this doesn't. */
			#ifndef __BUILTINS_LOADED
				#include <alpha/builtins.h>
			#endif
		#else
			#error unknown OS type for Alpha.
		#endif

		#include <c_asm.h> /* Needed for assembly code macros (e.g. prefetch)  */
		#ifndef __C_ASM_LOADED
			#error	Unable to locate or <c_asm.h> file from /usr/include !
		#endif

		/* If the #includes platform.h failed to capture __UMULH, can try to use this: */
		#ifndef __BUILTINS_LOADED
				#define  __UMULH(__x, __y) asm("umulh %a0, %a1, %v0;", (__x), (__y))
		#endif

	#elif(defined(__GNUC__) || defined(__GNUG__))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_GCC

		#define __UMULH(x, y)	\
		  ({ uint64 rslt;\
		  __asm__("umulh %1, %2, %0" : "=r" (rslt) : "r" (x), "r" (y)); rslt; })
	#else
		#error unknown compiler for Alpha.
		/*
		#define COMPILER_TYPE
		#define COMPILER_TYPE_UNKNOWN
		*/
	#endif

/* Ultrasparc, GNU C or C++ compiler */
#elif(defined(__sparc))

	#define CPU_TYPE
	#define CPU_IS_SPARC

	#define MUL64_SUCKS

	#if(defined(__SUNPRO_C))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_SUNC
	#elif(defined(__GNUC__) || defined(__GNUG__))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_GCC
	#else
		#define COMPILER_TYPE
		#define COMPILER_TYPE_UNKNOWN
	#endif

	/* Ultrasparc v8+ - we designate all v8+ as "sparc2," including ultra3: */
	#if(defined(__sparc_v9__))
		#define	CPU_SUBTYPE
		#define CPU_SUBTYPE_ULTRA2
		#define	CPU_SUBTYPE_NAME "Ultra2 or higher"
	#else
		#define	CPU_SUBTYPE
		#define CPU_VERSION_ULTRA1
		#define	CPU_SUBTYPE_NAME "Ultra1"
	#endif

#else
	/* Default is to assume crappy 64-bit integer MUL: */
	#define MUL64_SUCKS

	#define OS_TYPE
	#define OS_TYPE_UNKNOWN
	#define CPU_TYPE
	#define CPU_IS_UNKNOWN
	#define COMPILER_TYPE
	#define COMPILER_TYPE_UNKNOWN
#endif

/* Use fused floating mul/add? */
#ifdef FMADD_YES
	#ifdef HARDWARE_FMADD
		#define	USE_FMADD
	#else
		#error Fused mul/add not implemented/supported on this platform - please rebuild without FMADD_YES flag!
	#endif
#endif

/* Check that OS_TYPE, CPU_TYPE, COMPILER_TYPE are all defined: */
#ifndef OS_TYPE
	#error platform.h : OS_TYPE not defined!
#endif
#ifndef CPU_TYPE
	#error platform.h : CPU_TYPE not defined!
#endif
#ifndef COMPILER_TYPE
	#error platform.h : COMPILER_TYPE not defined!
#endif

/* Time to name names: */
#if(defined(OS_TYPE_WINDOWS))
  #if OS_DEBUG
	#warning	OS_NAME "Windows"
  #endif
	#define	OS_NAME "Windows"
#elif(defined(OS_TYPE_LINUX))
  #if OS_DEBUG
	#warning	OS_NAME "Linux"
  #endif
	#define	OS_NAME "Linux"
#elif(defined(OS_TYPE_MACOSX))
  #if OS_DEBUG
	#warning	OS_NAME "OS X"
  #endif
	#define	OS_NAME "OS X"
#elif(defined(OS_TYPE_DECOSF))/* DEC OSF (later HP TruUnix) */
  #if OS_DEBUG
	#warning	OS_NAME "DEC OSF / HP TruUnix"
  #endif
	#define	OS_NAME "DEC OSF / HP TruUnix"
#elif(defined(OS_TYPE_DECVMS))/* DEC VMS (originally for VAX) */
  #if OS_DEBUG
	#warning	OS_NAME "DEC VMS"
  #endif
	#define	OS_NAME "DEC VMS"
#elif(defined(OS_TYPE_HPUX))
  #if OS_DEBUG
	#warning	OS_NAME "HPUX"
  #endif
	#define	OS_NAME "HPUX"
#elif(defined(OS_TYPE_SUN))
  #if OS_DEBUG
	#warning	OS_NAME "SunOS / Solaris"
  #endif
	#define	OS_NAME "SunOS / Solaris"
#elif(defined(OS_TYPE_AIX))
  #if OS_DEBUG
	#warning	OS_NAME "AIX"
  #endif
	#define	OS_NAME "AIX"
#elif(defined(OS_TYPE_IRIX))
  #if OS_DEBUG
	#warning	OS_NAME "Irix"
  #endif
	#define	OS_NAME "Irix"
#elif(defined(OS_TYPE_UNKNOWN))
  #if OS_DEBUG
	#warning	OS_NAME "Unknown"
  #endif
	#define	OS_NAME "Unknown"
#endif

#ifndef OS_VERSION
	#define OS_VERSION "[Unknown]"
#endif

#if OS_BITS_DEBUG
	#ifndef OS_BITS
		#error OS_BITS not defined!
	#elif OS_BITS == 32
		#warning Compiling in 32-bit mode
	#elif OS_BITS == 64
		#warning Compiling in 64-bit mode
	#else
		#error OS_BITS defined but value not supported!
	#endif
#endif

#ifdef USE_SSE2

  #if defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

	#ifndef OS_BITS
		#error	Unable to determine if OS is 32-bit or 64-bit!
	#endif

	// For SSE2/GCC, double-check the value of OS_BITS against various int-type compiler intrinsics:
	#ifdef COMPILER_TYPE_GCC

	  #undef OS_BIT2

	  /* If GCC version predefines __SIZEOF_POINTER__, that is most reliable (the LONG_MAX-based test below failed under mingw64 because that defined LONG_MAX and LONG_LONG_MAX differently in 64-bit mode: */
	  #ifdef __SIZEOF_POINTER__

		#if __SIZEOF_POINTER__ == 4
			#define OS_BIT2 32
		#elif __SIZEOF_POINTER__ == 8
			#define OS_BIT2 64
		#else
			#error __SIZEOF_POINTER__ defined but returns unrecognized value! Use gcc -dM -E < /dev/null | grep __SIZEOF_POINTER__ to examine.
		#endif

	  /* Otherwise see if can use the value of __LONG_MAX__ in limits.h to quickly and portably determine whether it's 32-bit or 64-it OS: */
	  #else
		/* Syntax here is GCC and SunStudio/MSVC, respectively: */
		#if !defined(__LONG_MAX__) && !defined(LONG_MAX)
			#include <limits.h>
		#endif

		#if !defined(__LONG_MAX__) &&  defined(LONG_MAX)
			#define __LONG_MAX__  LONG_MAX
		#endif

		#ifdef __LONG_MAX__
			#if __LONG_MAX__ == 2147483647L
				#define OS_BIT2 32
			#elif __LONG_MAX__ == 9223372036854775807L
				#define OS_BIT2 64
			#else
				#error  __LONG_MAX__ defined but value unrecognized!
			#endif
		#else
			#error platform.h: failed to properly set OS_BIT2 for SSE2/GCC build!
		#endif

	  #endif

	  #if OS_BITS != OS_BIT2
		#error OS_BITS != OS_BIT2 for SSE2/GCC build!
	  #endif

	#endif

  #else

	#error SSE2 code not supported for this compiler!

  #endif

#endif

/* Multihreading support. */
extern int MAX_THREADS;
extern int NTHREADS;

// Needed to add the GCC flag here, otherwise MSVC tries to #include headers, etc even if USE_THREADS not defined:
#ifdef USE_THREADS

  #if defined(COMPILER_TYPE_GCC) && defined(OS_POSIX_COMPLIANT)

	// OpenMP requires USE_OMP to be def'd in addition to USE_THREADS:
	#ifdef USE_OMP

		// Found OpenMP header? The predefines here are for Linux, Sun Studio and Windows/MSVC, respectively:
		#include <omp.h>
		#if(defined(__OMP_H) || defined(_OPENMP) || defined(_OMPAPI))
			#define MULTITHREAD
			#define USE_OMP
		#else
			#error OpenMP header <omp.h> not detected - Platform may not provide OpenMP multithreading support.
		#endif

	#else	// Default is to use pthreads:

		/* Online docs [e.g. http://www.kernel.org/doc/man-pages/online/pages/man2/sched_setaffinity.2.html]
		tell us that in order to access macros for `cpu_set', we must #define _GNU_SOURCE before including <sched.h>.
		However, whether I define the above or not, on my distro (Fedora v16), I get these compile-time warnings
		indicating that the affinity stuff is not being included:

			"warning: implicit declaration of function ÔCPU_ZEROÕ", etc.

		Some sleuthing reveals that (at least in my distro) sched.h #includes <features.h>,
		and the latter header then defines __USE_GNU if _GNU_SOURCE is defined.
		In other words, #define _GNU_SOURCE simply causes __USE_GNU to also get defined.
		Or at least it should, but my above-quoted build diagnostics indicate otherwise.
		Even more bizarrely, when (in addition to defining just before including sched.h
		in my threading-related header file) I add -D_GNU_SOURCE to my compile command line,
		now all of a sudden the compiler sees both definitions and says

			platform.h:804:0: warning: "_GNU_SOURCE" redefined [enabled by default]
			<command-line>:0:0: note: this is the location of the previous definition

		...and the "implicit" warnings disppear. Anyway, add that one to the build-related mental "WTF?" bin
		and just #define __USE_GNU instead.
		*/
		#define __USE_GNU
		#include <sched.h>

		#include <sys/sysctl.h>
		#include <unistd.h>	// Needed for Posix sleep() command, among other things
		#ifdef OS_TYPE_LINUX

			// These additional Linux-only includes make sure __NR_gettid, used in our syscall-based get-thread-ID, is defined:
			#include <linux/unistd.h>
			#include <asm/unistd.h>

		#elif defined(OS_TYPE_MACOSX)

			// Mac OS X affinity API is via these:
			#include <mach/thread_policy.h>
			#include <mach/mach.h>
			// Gah - the osx mach/*h header tree sets its version of various CPU_IS_*** #defs,
			// incl. CPU_IS_SPARC, which then causes FP_MANTISSA_BITS_DOUBLE to get set to 53
			// rather the x86_64-required 64 ... that caused me to rename CPU_IS_* to CPU_IS_*.
		#endif

		#include <pthread.h>
		// Found pthread header?
		#if(defined(_PTHREAD_H))
			#define MULTITHREAD
			#define USE_PTHREAD
		#else
			#error Pthreads header <pthread.h> not detected - Platform may not provide multithreading support.
		#endif

	  #endif	// OpenMP or Pthread

	#else

		#error Multithreading currently only supported for Linux/GCC builds!

	#endif

#endif

#if defined(CPU_IS_ALFA)
  #if CPU_DEBUG
	#warning	CPU_NAME "Alpha"
  #endif
	#define	CPU_NAME "Alpha"
#elif defined(CPU_IS_MIPS)
  #if CPU_DEBUG
	#warning	CPU_NAME "Mips"
  #endif
	#define	CPU_NAME "Mips"
#elif defined(CPU_IS_SPARC)
  #if CPU_DEBUG
	#warning	CPU_NAME "Sparc"
  #endif
	#define	CPU_NAME "Sparc"
#elif defined(CPU_IS_PPC)
  #if CPU_DEBUG
	#warning	CPU_NAME "PowerPC"
  #endif
	#define	CPU_NAME "PowerPC"
#elif defined(CPU_IS_X86)
  #if CPU_DEBUG
	#warning	CPU_NAME "x86"
  #endif
	#define	CPU_NAME "x86"
	#define	FP_MANTISSA_BITS_DOUBLE	64
#elif defined(CPU_IS_X86_64)
  #if CPU_DEBUG
	#warning	CPU_NAME "x86_64"
  #endif
	#define	CPU_NAME "x86_64"
	#define	FP_MANTISSA_BITS_DOUBLE	64
#elif defined(CPU_IS_POWER)
  #if CPU_DEBUG
	#warning	CPU_NAME "Power"
  #endif
	#define	CPU_NAME "Power"
#elif defined(CPU_IS_HPPA)
  #if CPU_DEBUG
	#warning	CPU_NAME "HPPA"
  #endif
	#define	CPU_NAME "HPPA"
#elif defined(CPU_IS_IA64)
  #if CPU_DEBUG
	#warning	CPU_NAME "Itanium"
  #endif
	#define	CPU_NAME "Itanium"
	#define	FP_MANTISSA_BITS_DOUBLE	64
#elif defined(CPU_IS_NVIDIA)
  #if CPU_DEBUG
	#warning	CPU_NAME "nVidia"
  #endif
	#define	CPU_NAME "nVidia"
#elif defined(CPU_IS_UNKNOWN)
  #if CPU_DEBUG
	#warning	CPU_NAME "Unknown"
  #endif
	#define	CPU_NAME "Unknown"
#endif

#ifndef	FP_MANTISSA_BITS_DOUBLE
	#define	FP_MANTISSA_BITS_DOUBLE	53
#endif

#if(defined(COMPILER_TYPE_NVCC))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "nVidia C++ (for GPU)"
  #endif
	#define COMPILER_NAME "nVidia C++ (for GPU)"
#elif(defined(COMPILER_TYPE_GCC))
  #ifdef __clang__
   #if CMPLR_DEBUG
	#warning	COMPILER_NAME "Gnu-C-compatible [llvm/clang]"
   #endif
	#define COMPILER_NAME "Gnu-C-compatible [llvm/clang]"
	#ifdef	__clang_version__
		#define	COMPILER_VERSION	__clang_version__
	#endif
  #else
   #if CMPLR_DEBUG
	#warning	COMPILER_NAME "Gnu C [or other compatible]"
   #endif
	#define COMPILER_NAME "Gnu C [or other compatible]"
	#ifdef	__VERSION__
		#define	COMPILER_VERSION	__VERSION__
	#endif
  #endif
#elif(defined(COMPILER_TYPE_MWERKS))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "Metrowerks Codewarrior"
  #endif
	#define COMPILER_NAME "Metrowerks Codewarrior"
#elif(defined(COMPILER_TYPE_ICC))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "Intel C"
  #endif
	#define COMPILER_NAME "Intel C"
#elif(defined(COMPILER_TYPE_HPC))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "HP C"
  #endif
	#define COMPILER_NAME "HP C"
#elif(defined(COMPILER_TYPE_APPLEC))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "Apple C"
  #endif
	#define COMPILER_NAME "Apple C"
#elif(defined(COMPILER_TYPE_XLC))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "IBM XL-C"
  #endif
	#define COMPILER_NAME "IBM XL-C"
#elif(defined(COMPILER_TYPE_DECC))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "DECC/HP-C"
  #endif
	#define COMPILER_NAME "DECC/HP-C"
#elif(defined(COMPILER_TYPE_SUNC))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "SunPro C"
  #endif
	#define COMPILER_NAME "SunPro C"
#elif(defined(COMPILER_TYPE_MSVC))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "MSVC/.NET"
  #endif
	#define COMPILER_NAME "MSVC/.NET"

  #if(_MSC_VER < 1300)	/* Use of MSVC 7+ required, mainly so we don't get compiler errors about 'll' extensions to 64-bit const ints. */
	#error Your version of MSVC appears to be out of date - Use of MSVC/.NET v7.0 or later is ***required***
  #endif

	#ifndef _MSC_VER
		#error _MSC_VER not defined!
	#elif(_MSC_VER == 1020)
		#define COMPILER_VERSION "4.2"
	#elif(_MSC_VER <= 1100)
		#define COMPILER_VERSION "5"
	#elif(_MSC_VER <= 1200)
		#define COMPILER_VERSION "6"
	#elif(_MSC_VER <= 1300)
		#define COMPILER_VERSION "7"
	#elif(_MSC_VER <= 1400)
		#define COMPILER_VERSION "8"
	#else
		#define COMPILER_VERSION "Unrecognized"
	#endif

#elif(defined(COMPILER_TYPE_UNKNOWN))
  #if CMPLR_DEBUG
	#warning	COMPILER_NAME "Unknown Compiler"
  #endif
	#define	COMPILER_NAME "Unknown Compiler"
#endif

#ifndef COMPILER_VERSION
	#define COMPILER_VERSION "[Unknown]"
#endif

// SIMD code only available for 64-bit GCC build:
#if defined(USE_SSE2) && !defined(COMPILER_TYPE_GCC)
	#error  SIMD code only available for GCC (and compatible) build!
#endif

// SIMD+AVX code only available for 64-bit build:
#if defined(USE_AVX) && (OS_BITS == 32)
	#error  SIMD+AVX code only available for 64-bit build!
#endif

#endif	/* platform_h_included */

/* List of headers used by Mprime, by way of reference:

	From: 	George Woltman <woltman@alum.mit.edu>
	Subject: 	Re: thread affinity
	Date: 	November 5, 2012 8:10:10 PM PST
	To: 	E. Mayer <ewmayer@aol.com>

	[snip of some thread-affinity-related stuff]

	The last mprime built for the Mac included these libs:

	LIBS   = ../gwnum/release/gwnum.a -lm -lpthread -lcurl -lIOKit -framework CoreFoundation -lstdc++

	and in case CPU_SET etc are macros, I include these files at compile time:

	// Include files needed by all ports
	#include "prime.h"
	#include <ctype.h>
	#include <fcntl.h>
	#include <math.h>
	#include <memory.h>
	#include <signal.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <sys/stat.h>

	// Required Linux files
	#ifdef __linux__
	#include <dirent.h>
	#include <unistd.h>
	#include <linux/unistd.h>
	#include <asm/unistd.h>
	#define __USE_GNU
	#include <sched.h>
	#include <sys/resource.h>
	#include <sys/sysinfo.h>
	#include <sys/time.h>
	#include <sys/timeb.h>
	#endif

	// Required Mac OS X files
	#ifdef __APPLE__
	#include <dirent.h>
	#include <pthread.h>
	#include <sched.h>
	#include <unistd.h>
	#include <sys/types.h>
	#include <sys/resource.h>
	#include <sys/sysctl.h>
	#include <sys/time.h>
	#include <sys/timeb.h>
	#include <CoreFoundation/CoreFoundation.h>
	#include <IOKit/ps/IOPowerSources.h>
	#include <IOKit/ps/IOPSKeys.h>
	#endif

	// Required FreeBSD files
	#ifdef __FreeBSD__
	#include <dirent.h>
	#include <pthread.h>
	#include <sched.h>
	#include <unistd.h>
	#include <sys/param.h>
	#include <sys/cpuset.h>
	#include <sys/resource.h>
	#include <sys/sysctl.h>
	#include <sys/types.h>
	#include <sys/time.h>
	#include <sys/timeb.h>
	#endif

	Hope that helps,
	George
*/

