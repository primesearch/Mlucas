/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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

/* Only one of the following 3 should be set = 1 at any time.
   If > 1 is set, only the first onesuch will be respected. */
/* Set = 1 to print brief OS summary at compile time and exit: */
#undef	OS_DEBUG
#define	OS_DEBUG	0

/* Set = 1 to print brief OS summary at compile time and exit: */
#undef	CPU_DEBUG
#define	CPU_DEBUG	0

/* Set = 1 to print brief OS summary at compile time and exit: */
#undef	CMPLR_DEBUG
#define	CMPLR_DEBUG	0

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
#endif

/* See if can use the value of __LONG_MAX__ in limits.h to quickly and portably determine whether it's 32-bit or 64-it OS.
Currently this is needed only for the SSE2 GCC-style assembly code, so wrap it in an appropriate #define:
*/
//#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC)

	#undef	OS_BITS

	/* Syntax here is GCC and SunStudio/MSVC, respectively: */
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
		#if(sizeof(long) == 4)
			#define OS_BITS 32
		#elif(sizeof(long) == 8)
			#define OS_BITS 64
		#else
			#error  Value of sizeof(long) not recognized!
		#endif
	#endif

	#ifndef OS_BITS
		#error platform.h: failed to properly set OS_BITS!
	#endif
/*
	#if OS_BITS == 32
		#error 32-bit OS detected
	#elif OS_BITS == 64
		#error 64-bit OS detected
	#endif
*/
//#endif

/* Multihreading support. */
extern int MAX_THREADS;
extern int NTHREADS;

#undef MULTITHREAD	/* Master #define controlling multithreaded mode -- enable by invoking USE_THREADS
					and the appropriate compiler options (e.g. -openmp; compiler-dependent) at build time. */

#if defined(USE_THREADS)

	#include <omp.h>
	/* Found OpenMP header? The predefines here are for Linux, Sun Studio and Windows/MSVC, respectively: */
	#if(defined(__OMP_H) || defined(_OPENMP) || defined(_OMPAPI))
		#define MULTITHREAD
	#else
		#error OpenMP header <omp.h> not detected - Platform may not provide multithreading support.
	#endif

#endif

/*
Locally defined CPU types - undef these here to prevent namespace collisions.
Each of the above CPUs typically comes in multiple flavors: try to restrict
our list here to the salient ones, that actually need to be differentiated
based on different key capabilities. Default CPU subtypes are as indicated.

More than one of the subtypes may be def'd simultaneously, e.g. for CPU_TYPE_PPC we
might have CPU_SUBTYPE_PPC32 and CPU_HAS_ALTIVEC def'd simultaneously (for G4, not G3).
Mutual exclusivity of subtypes is implied via their being listed in one block, e.g.
(again using the PowerPC as an example) CPU_SUBTYPE_PPC32 and CPU_SUBTYPE_PPC64
are mutually incompatible.
*/
#undef	CPU_NAME
#undef	CPU_SUBTYPE_NAME
#undef	CPU_TYPE
#undef	CPU_SUBTYPE
#undef	CPU_TYPE_UNKNOWN
	#undef	CPU_SUBTYPE_UNKNOWN
#undef	CPU_TYPE_ALFA
	#undef	CPU_SUBTYPE_EV4
	#undef	CPU_SUBTYPE_EV5
	#undef	CPU_SUBTYPE_EV6
#undef	CPU_TYPE_MIPS
	#undef	CPU_SUBTYPE_PRE_R10K
	#undef	CPU_SUBTYPE_R10K
	#undef	CPU_SUBTYPE_R12K
	#undef	CPU_SUBTYPE_R15K
#undef	CPU_TYPE_SPARC
	#undef	CPU_SUBTYPE_ULTRA1
	#undef	CPU_SUBTYPE_ULTRA2
	#undef	CPU_SUBTYPE_ULTRA3
#undef	CPU_TYPE_PPC
	#undef	CPU_SUBTYPE_PPC32
	#undef	CPU_SUBTYPE_PPC64	/* This is the main differentiator we need for now,
								i.e. 64-bit or not? Separate check of AltiVec-or-not
								differentiates between G3 and G4 (both of which are PPC32). */
#undef	CPU_TYPE_IA32
	#undef	CPU_SUBTYPE_CLASSIC
	#undef	CPU_SUBTYPE_PPRO
	#undef	CPU_SUBTYPE_P2
	#undef	CPU_SUBTYPE_P3
	#undef	CPU_SUBTYPE_P4
#undef	CPU_TYPE_IA64
	#undef	CPU_SUBTYPE_IT1
	#undef	CPU_SUBTYPE_IT2
	#undef	CPU_TYPE_X86_64

#undef	CPU_TYPE_AMD32
#undef	CPU_TYPE_AMD64
/* DO WE NEED ANY OF THESE?
	#undef	CPU_SUBTYPE_K7
	#undef	CPU_SUBTYPE_K8
	#undef	CPU_SUBTYPE_OPTERON
*/
#undef	CPU_TYPE_POWER
#undef	CPU_TYPE_POWER
	#undef	CPU_SUBTYPE_POWER1
	#undef	CPU_SUBTYPE_POWER2
	#undef	CPU_SUBTYPE_POWER3
	#undef	CPU_SUBTYPE_POWER4
	#undef	CPU_SUBTYPE_POWER5
#undef	CPU_TYPE_HPPA
	#undef	CPU_SUBTYPE_PA1
	#undef	CPU_SUBTYPE_PA2

/* SIMD COPROCESSOR AND MULTIMEDIA EXTENSION FLAGS:
   AltiVec is shared among several CPU families (e.g. PowerPC and IBM Power)
   AS IS MMX (Pentium and IA64) so don't attach them to any specific CPU type:
*/
#undef	CPU_HAS_SSE1
#undef	CPU_HAS_SSE2
#undef	CPU_HAS_SSE3
#undef	CPU_HAS_ALTIVEC
	/* According to the Motorola AltiVec PIM, any AVec-enabled compiler predefines:
	#if(defined(__VEC__) && (__VEC__ == 10205))
	*/
	/* On the G5/Mac OS X, __VEC__ == 10206, and __ALTIVEC__ is defined -
	probably better to use the latter #define to identify AVec on PowerPC chips:
	*/
	#if(defined(__ALTIVEC__))
		#define	CPU_HAS_ALTIVEC		1
	#endif
#undef	CPU_HAS_MMX
#undef	CPU_HAS_CELL


/* Locally defined Compiler types and versions: */
#undef	COMPILER_NAME

#undef	COMPILER_TYPE
#undef	COMPILER_TYPE_UNKNOWN

#undef	COMPILER_VERSION

#undef	COMPILER_TYPE_GCC		/* Gnu C compiler */
#undef	COMPILER_TYPE_MWERKS	/* Metrowerks Codewarrior */
#undef	COMPILER_TYPE_ICC		/* Intel compiler */
#undef	COMPILER_TYPE_HPC		/* HP C compiler */
#undef	COMPILER_TYPE_APPLEC	/* Apple C compiler */
#undef	COMPILER_TYPE_XLC		/* IBM XL C compiler */
#undef	COMPILER_TYPE_DECC		/* DEC C compiler */
#undef	COMPILER_TYPE_SUNC		/* SunPro C compiler */
#undef	COMPILER_TYPE_MSVC		/* Microsoft Visual C++ (later .NET) compiler */
/* Miscellaneous notes on the key-but-lesser-known features of the compilers:
	- MWERKS and MSVC/.NET share the same inline ASM syntax, as do GCC and ICC.
	-
*/

/* 32-bit X86, Gnu C or Metrowerks CodeWarrior C compiler: */
#if(defined(__i386) || defined(__i386__))

	#ifndef OS_BITS
		#define OS_BITS 32
	#endif

	#define CPU_TYPE
	#define CPU_TYPE_IA32
	#define	INTEGER_MUL_32

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

/* AMD64 ISA - includes "Intel 64" (formerly known as "EMT64", formerly formerly known as "IA32e", yada, yada): */
#elif(defined(__amd64) || defined(__amd64__) || defined(_M_AMD64) || defined(_M_EMT64) || defined(__x86_64) || defined(__x86_64__))

	#ifndef OS_BITS
		#define OS_BITS 64
	#endif

	#define CPU_TYPE
	#define CPU_TYPE_AMD64

	#define MULH64_FAST
	#define	HARDWARE_FMADD

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
	/* MS Visual C/Studio/.NET/Whatever: */
	#elif(defined(_MSC_VER))
		#define COMPILER_TYPE
		#define COMPILER_TYPE_MSVC
	/* Unknown: */
	#else
		#define COMPILER_TYPE
		#define COMPILER_TYPE_UNKNOWN
	#endif

/* 32-bit X86, Intel C or MSVC/.NET compiler */
#elif(defined(_M_IX86) || defined(_MSC_VER))

	#ifndef OS_BITS
		#define OS_BITS 32
	#endif

	#define CPU_TYPE
	#define CPU_TYPE_IA32
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

		#if(defined (__linux__))
			#include <sched.h>		/* Needed to set CPU affinity in multithreaded mode */
		#endif
		/*
		Various header files containing ICC intrinsics for IA32:
		we include individual ones of these as needed by the various
		functionality modules (e.g. prefetch.h and sse2.h):
		#include <dvec.h>					* SSE 2 intrinsics for Class Libraries												*
		#include <fvec.h>					* SSE intrinsics for Class Libraries												*
		#include <ia32intrin.h>				* Miscellaneous IA-32 intrinsics.													*
		#include <ivec.h>					* MMXô instructions intrinsics for Class Libraries									*
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
	#define CPU_TYPE_IA64

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

		#if(defined (__linux__))
			#include <sched.h>		/* Needed to set CPU affinity in multithreaded mode */
		#endif
		/*
		Various header files containing ICC intrinsics for IA64:
		we include individual ones of these as needed by the various
		functionality modules (e.g. prefetch.h and sse2.h):
		#include <emmintrin.h>				* Principal header file for SSE2 intrinsics										*
		#include <fvec.h>					* SSE intrinsics for Class Libraries											*
		#include <ia64intrin.h>				* Miscellaneous IA-64 intrinsics.												*
		#include <ia64regs.h>				* Header file for registers on Itanium-based systems							*
		#include <ivec.h>					* MMXô instructions intrinsics for Class Libraries								*
		#include <mathimf.h>				* Principal header file for current Intel Math Library							*
		#include <mmclass.h>				* Principal header files for MMXô instructions with Class Libraries				*
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
	#define CPU_TYPE_HPPA
	#define COMPILER_TYPE
	#define COMPILER_TYPE_HPC

/* IBM Power: Note that I found that only __xlc__ was properly defined on the power5 I used for my tests.
We Deliberately put this #elif ahead of the PowerPC one because the Power/AIX compiler also defines some PowerPC flags: */
#elif(defined(_POWER) || defined(_ARCH_PWR))

	#ifndef	OS_TYPE_AIX
		#error ARCH = POWER defined but not OS_TYPE_AIX!
	#endif

	#define CPU_TYPE
	#define CPU_TYPE_POWER

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
	#define CPU_TYPE_PPC

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
	#define CPU_TYPE_ALFA

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
	#define CPU_TYPE_SPARC

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
	#define CPU_TYPE_UNKNOWN
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
	#error	OS_NAME "Windows"
  #else
	#define	OS_NAME "Windows"
  #endif
#elif(defined(OS_TYPE_LINUX))
  #if OS_DEBUG
	#error	OS_NAME "Linux"
  #else
	#define	OS_NAME "Linux"
  #endif
#elif(defined(OS_TYPE_MACOSX))
  #if OS_DEBUG
	#error	OS_NAME "OS X"
  #else
	#define	OS_NAME "OS X"
  #endif
#elif(defined(OS_TYPE_DECOSF))/* DEC OSF (later HP TruUnix) */
  #if OS_DEBUG
	#error	OS_NAME "DEC OSF / HP TruUnix"
  #else
	#define	OS_NAME "DEC OSF / HP TruUnix"
  #endif
#elif(defined(OS_TYPE_DECVMS))/* DEC VMS (originally for VAX) */
  #if OS_DEBUG
	#error	OS_NAME "DEC VMS"
  #else
	#define	OS_NAME "DEC VMS"
  #endif
#elif(defined(OS_TYPE_HPUX))
  #if OS_DEBUG
	#error	OS_NAME "HPUX"
  #else
	#define	OS_NAME "HPUX"
  #endif
#elif(defined(OS_TYPE_SUN))
  #if OS_DEBUG
	#error	OS_NAME "SunOS / Solaris"
  #else
	#define	OS_NAME "SunOS / Solaris"
  #endif
#elif(defined(OS_TYPE_AIX))
  #if OS_DEBUG
	#error	OS_NAME "AIX"
  #else
	#define	OS_NAME "AIX"
  #endif
#elif(defined(OS_TYPE_IRIX))
  #if OS_DEBUG
	#error	OS_NAME "Irix"
  #else
	#define	OS_NAME "Irix"
  #endif
#elif(defined(OS_TYPE_UNKNOWN))
  #if OS_DEBUG
	#error	OS_NAME "Unknown"
  #else
	#define	OS_NAME "Unknown"
  #endif
#endif

#ifndef OS_VERSION
	#define OS_VERSION "[Unknown]"
#endif

#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC)
	#ifndef OS_BITS
		#error	Unable to determine if OS is 32-bit or 64-bit!
	#endif
#endif

#if(defined(CPU_TYPE_ALFA))
  #if CPU_DEBUG
	#error	CPU_NAME "Alpha"
  #else
	#define	CPU_NAME "Alpha"
  #endif
#elif(defined(CPU_TYPE_MIPS))
  #if CPU_DEBUG
	#error	CPU_NAME "Mips"
  #else
	#define	CPU_NAME "Mips"
  #endif
#elif(defined(CPU_TYPE_SPARC))
  #if CPU_DEBUG
	#error	CPU_NAME "Sparc"
  #else
	#define	CPU_NAME "Sparc"
  #endif
#elif(defined(CPU_TYPE_PPC))
  #if CPU_DEBUG
	#error	CPU_NAME "PowerPC"
  #else
	#define	CPU_NAME "PowerPC"
  #endif
#elif(defined(CPU_TYPE_IA32))
  #if CPU_DEBUG
	#error	CPU_NAME "Pentium"
  #else
	#define	CPU_NAME "Pentium"
	#define	FP_MANTISSA_BITS_DOUBLE	64
  #endif
#elif(defined(CPU_TYPE_AMD))
  #if CPU_DEBUG
	#error	CPU_NAME "AMD"
  #else
	#define	CPU_NAME "AMD"
	#define	FP_MANTISSA_BITS_DOUBLE	64
  #endif
#elif(defined(CPU_TYPE_AMD64))
  #if CPU_DEBUG
	#error	CPU_NAME "AMD64"
  #else
	#define	CPU_NAME "AMD64"
	#define	FP_MANTISSA_BITS_DOUBLE	64
  #endif
#elif(defined(CPU_TYPE_X86_64))
  #if CPU_DEBUG
	#error	CPU_NAME "x86-64"
  #else
	#define	CPU_NAME "x86-64"
	#define	FP_MANTISSA_BITS_DOUBLE	64
  #endif
#elif(defined(CPU_TYPE_POWER))
  #if CPU_DEBUG
	#error	CPU_NAME "Power"
  #else
	#define	CPU_NAME "Power"
  #endif
#elif(defined(CPU_TYPE_HPPA))
  #if CPU_DEBUG
	#error	CPU_NAME "HPPA"
  #else
	#define	CPU_NAME "HPPA"
  #endif
#elif(defined(CPU_TYPE_IA64))
  #if CPU_DEBUG
	#error	CPU_NAME "Itanium"
  #else
	#define	CPU_NAME "Itanium"
	#define	FP_MANTISSA_BITS_DOUBLE	64
  #endif
#elif(defined(CPU_TYPE_UNKNOWN))
  #if CPU_DEBUG
	#error	CPU_NAME "Unknown"
  #else
	#define	CPU_NAME "Unknown"
  #endif
#endif

#ifndef	FP_MANTISSA_BITS_DOUBLE
	#define	FP_MANTISSA_BITS_DOUBLE	53
#endif

#if(defined(COMPILER_TYPE_GCC))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "Gnu C"
  #else
	#define COMPILER_NAME "Gnu C"
  #endif
	#ifdef	__VERSION__
		#define	COMPILER_VERSION	__VERSION__
	#endif
#elif(defined(COMPILER_TYPE_MWERKS))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "Metrowerks Codewarrior"
  #else
	#define COMPILER_NAME "Metrowerks Codewarrior"
  #endif
#elif(defined(COMPILER_TYPE_ICC))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "Intel C"
  #else
	#define COMPILER_NAME "Intel C"
  #endif
#elif(defined(COMPILER_TYPE_HPC))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "HP C"
  #else
	#define COMPILER_NAME "HP C"
  #endif
#elif(defined(COMPILER_TYPE_APPLEC))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "Apple C"
  #else
	#define COMPILER_NAME "Apple C"
  #endif
#elif(defined(COMPILER_TYPE_XLC))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "IBM XL-C"
  #else
	#define COMPILER_NAME "IBM XL-C"
  #endif
#elif(defined(COMPILER_TYPE_DECC))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "DECC/HP-C"
  #else
	#define COMPILER_NAME "DECC/HP-C"
  #endif
#elif(defined(COMPILER_TYPE_SUNC))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "SunPro C"
  #else
	#define COMPILER_NAME "SunPro C"
  #endif
#elif(defined(COMPILER_TYPE_MSVC))
  #if CMPLR_DEBUG
	#error	COMPILER_NAME "MSVC/.NET"
  #else
	#define COMPILER_NAME "MSVC/.NET"
  #endif

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
	#error	COMPILER_NAME "Unknown Compiler"
  #else
	#define	COMPILER_NAME "Unknown Compiler"
  #endif
#endif

#ifndef COMPILER_VERSION
	#define COMPILER_VERSION "[Unknown]"
#endif

#endif	/* platform_h_included */

