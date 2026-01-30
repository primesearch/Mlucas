/*******************************************************************************
*                                                                              *
*   (C) 1997-2021 by Ernst W. Mayer.                                           *
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

/* Implements some basic prefetching on various architectures: */

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef prefetch_h_included
#define prefetch_h_included

#include "platform.h"

#undef	PFETCH

/* Pentium3 family, Gnu C or C++ compiler */
#if(defined(CPU_IS_X86) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC) ) )

	# define PFETCH 1

	# define CACHE_LINE_INTS    8
	# define CACHE_LINE_DOUBLES 4
	# define CACHE_LINE_COMPLEX 2

	# define prefetch_p_doubles(_pd)	__asm__ volatile ("prefetchnta (%0)": :"r"(_pd+2));

/* Pentium3 family, Metrowerks CodeWarrior C compiler: */
#elif(defined(CPU_IS_X86) && defined(COMPILER_TYPE_MWERKS))

	# define PFETCH 1

	# define CACHE_LINE_INTS    8
	# define CACHE_LINE_DOUBLES 4
	# define CACHE_LINE_COMPLEX 2

	# define prefetch_p_doubles(_pd)\
	{\
		/* EWM - need the extra {} to get CW to inline these sans errors: */\
		{	__asm	mov edx, _pd	};	\
		{	__asm	prefetcht0 [edx]};	\
	}

/* Pentium3 family, Microsoft Visual C compiler - uses CodeWarrior-style prefetch macros, sans extra {}.
Note that MSVC doesn't like the argument to the macro to be an arithmetic expression, i.e. need _pd
to point to a simple pointer variable:
*/
#elif(defined(CPU_IS_X86) && defined(COMPILER_TYPE_MSVC))

  #ifdef MULTITHREAD
	# define PFETCH 0	/* "error C3011: inline assembly not allowed directly within a parallel region" - Gee, thanks, MSFT. */
  #else
	# define PFETCH 1	/* Explicit prefetch is consistently slower on my p4 using MSVC 2005 */
  #endif

	# define CACHE_LINE_INTS    8
	# define CACHE_LINE_DOUBLES 4
	# define CACHE_LINE_COMPLEX 2

  #ifdef PFETCH
	# define prefetch_p_doubles(_pd)\
	{\
		__asm	mov edx, _pd   	\
		__asm	prefetcht0 [edx]\
	}
  #else
	# define prefetch_p_doubles(_pd)	/* NO-OP! */ do { (void)(_pd); } while (0)
  #endif

/* AMD64/Solaris, Sun C compiler */
#elif(defined(CPU_IS_X86_64) && defined(CPU_SUBTYPE_AMD64) && defined(COMPILER_TYPE_SUNC))

    # define PFETCH 1

    # define CACHE_LINE_INTS    8
    # define CACHE_LINE_DOUBLES 4
    # define CACHE_LINE_COMPLEX 2

    # define prefetch_p_doubles(_pd)    __asm volatile ("prefetchw (%0)": :"r"(_pd+2));

/* X86 family, Intel C compiler (ICC v7 or above). We need to put this
before the COMPILER_TYPE_GCC case because, bizarrely, ICC also #defines __GNUC__ .
*/
#elif((defined(CPU_IS_X86) || defined(CPU_IS_IA64) || defined(CPU_IS_X86_64)) && defined(COMPILER_TYPE_ICC))

	# define CACHE_LINE_INTS    16
	# define CACHE_LINE_DOUBLES  8
	# define CACHE_LINE_COMPLEX  4

	#if(defined(CPU_IS_IA64))
		#define PFETCH 1
		#define prefetch_p_doubles(_pd)		__lfetch(0x01, _pd + 64)
	#else
	/*
		The _MM_HINT_* constants for use with _mm_prefetch are not defined under ecc v6,
		   but since the pre-version-7 compiler was crap anyway, it's moot.
		   The various options here are as follows:

			_MM_HINT_T1     0    Temporal locality, level 1
			_MM_HINT_NT1    1    No temporal locality, level 1
			_MM_HINT_NT2    2    No temporal locality, level 2
			_MM_HINT_NTA    3    No temporal locality, all levels

		In my timing tests of 4/05/2005 using the v8 compiler,
		default (i.e. no-explicit-prefetch), T1, NT1 and NT2 all gave roughly identical results
		(NT1 seems a tad faster overall than no-explicit, but we're talking the 1-2% range -
		anyway, that's why we use it as the preferred option here), but NTA proved to be
		consistently 1-5% slower than the above - NOT RECOMMENDED.
	*/
		#define PFETCH 1
		# define prefetch_p_doubles(_pd)	_mm_prefetch((char *)(_pd + 2), 1);
	#endif

/* AMD64 family, Gnu C or C++ compiler */
#elif(defined(CPU_IS_X86_64) && defined(CPU_SUBTYPE_AMD64) && defined(COMPILER_TYPE_GCC))

	# define PFETCH 1

	# define CACHE_LINE_INTS    8
	# define CACHE_LINE_DOUBLES 4
	# define CACHE_LINE_COMPLEX 2

	# define prefetch_p_doubles(_pd)	__asm__ volatile ("prefetchw (%0)": :"r"(_pd+2));

/* AMD64 family, MSVC */
#elif(defined(CPU_IS_X86_64) && defined(CPU_SUBTYPE_AMD64) && defined(COMPILER_TYPE_MSVC))

# define PFETCH 1

	# define CACHE_LINE_INTS    8
	# define CACHE_LINE_DOUBLES 4
	# define CACHE_LINE_COMPLEX 2

  #ifdef PFETCH
	# define prefetch_p_doubles(_pd)\
	{\
		__asm	mov edx, _pd   	\
		__asm	prefetchw [edx]\
	}
  #else
	# define prefetch_p_doubles(_pd)	/* NO-OP! */ do { (void)(_pd); } while (0)
  #endif

/* AMD Athlon/Duron family, Metrowerks CodeWarrior C compiler */
/* Since the x86 prefetch seems to work as well as any of the Athlon-specific ones, skip this for now. */
#elif(defined(CPU_IS_AMD) && defined(COMPILER_TYPE_MWERKS))

	#error	Prefetch for CPU_IS_AMD currently not supported!

	# define PFETCH 1

	# define CACHE_LINE_INTS    8
	# define CACHE_LINE_DOUBLES 4
	# define CACHE_LINE_COMPLEX 2

	# define prefetch_p_doubles(_pd)\
	{\
		{	__asm	mov edx, _pd	};	\
		{	__asm	prefetchW [edx]};	\
	}

/* IA-64 (Itanium), HP C or C++ compiler for HPUX: */
#elif(defined(CPU_IS_IA64) && defined(COMPILER_TYPE_HPC))

	# define CACHE_LINE_INTS    16
	# define CACHE_LINE_DOUBLES  8
	# define CACHE_LINE_COMPLEX  4

	#if 1
	/*
		From HP's Intel Itanium® processor family Performance Tuning Guide White Paper:

			// request a cache line to be loaded (4th one from where
			// we are =((128bytes/cacheline)*4 cachelines)/sizeof float)

			__lfetch(0x01, v+128);		<*** What does the '0x01' specify?

		...so for doubles, prefetching 4 lines ahead would mean __lfetch(0x01, address + 64).

		An additional important note in the above document (emphasis mine):

			"If you are developing code from scratch, or performing a significant rewrite of existing code, it makes
		sense to try to organize the data structures so that data that is accessed together is also grouped within
		the same cache line. This will make it much easier for the compiler to effectively pre-fetch. In addition, IT IS
		VERY IMPORTANT TO SEGREGATE INTEGER AND FLOATING-POINT VARIABLES SO THAT THEY DO NOT FALL IN THE SAME CACHE
		LINE. Remember, integer loads/stores are serviced out of the L1 cache, and floating-point loads/stores are
		serviced out of the L2 cache. If a floating-point store is done to a cache line that is resident in both L1 and
		L2, the line in L1 will be cast out. A subsequent reference to integer data in that cache line will cause the
		line to be reloaded, and an unnecessary stall will result."
	*/
		#define PFETCH 1
		#define prefetch_p_doubles(_pd)		__lfetch(0x01, _pd + 64)
	#else
	/*
		From HP's "inline assembly for Itanium-based HP-UX":

		Table 1-14 _Asm_lftype – legal lfetch lftype completers:

		typedef enum {
			_LFTYPE_NONE  = 0,	// Ignore faults
			_LFTYPE_FAULT = 1	// Raise faults
		} _Asm_lftype;

		Table 1-15 _Asm_lfhint – legal lfetch lfhint completers:

		typedef enum {
			_LFHINT_NONE = 0,	// Temporal locality, level 1
			_LFHINT_NT1 = 1,	// No temporal locality, level 1
			_LFHINT_NT2 = 2,	// No temporal locality, level 2
			_LFHINT_NTA = 3		// No temporal locality, all levels
		} _Asm_lfhint;

		Table 1-31 Memory Management Opcodes:
		...
		void _Asm_lfetch		( _Asm_lftype, _Asm_lfhint, void * r3);
		void _Asm_lfetch_excl	( _Asm_lftype, _Asm_lfhint, void * r3);
	*/
		#define PFETCH 1
		#define prefetch_p_doubles(_pd)	_Asm_lfetch_excl(_LFTYPE_NONE, _LFHINT_NT2 , (void *)(_pd    +  4))
	#endif

/* IA-64 (Itanium), Gnu C or C++ compiler. We need to put this
after the IA64/ICC case because, bizarrely, ICC also #defines __GNUC__ .
*/
#elif(defined(CPU_IS_IA64) && defined(COMPILER_TYPE_GCC))

	# define PFETCH 1

	# define CACHE_LINE_INTS    16
	# define CACHE_LINE_DOUBLES  8
	# define CACHE_LINE_COMPLEX  4

	# define prefetch_p_doubles(_pd)	__asm__ volatile ("lfetch.nt2 [%0]": :"r"(_pd+2));

/* PowerPC, Apple/Gnu C or IBM XLC compiler: */
#elif(defined(CPU_IS_PPC) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_XLC)))

	# define PFETCH 1

	# define CACHE_LINE_INTS    8
	# define CACHE_LINE_DOUBLES 4
	# define CACHE_LINE_COMPLEX 2

	/* Use << 3 here rather than *sizeof(double). */
	# define prefetch_p_doubles(_pd)	__asm__ volatile ("dcbt  0,%1": :"r"(_pd), "r"(CACHE_LINE_DOUBLES << 3));

/* PowerPC, Metrowerks CodeWarrior or MSVC compiler */
#elif(defined(CPU_IS_PPC) && (defined(COMPILER_TYPE_MWERKS) || defined(COMPILER_TYPE_MSVC)))

	# define PFETCH 1

	# define CACHE_LINE_INTS    8
	# define CACHE_LINE_DOUBLES 4
	# define CACHE_LINE_COMPLEX 2

	/* Use << 3 here rather than *sizeof(double). */
	# define prefetch_p_doubles(_pd)	__dcbt(_pd,CACHE_LINE_DOUBLES << 3);

/* Alpha system, DEC/Compaq C compiler */
#elif(defined(CPU_IS_ALFA) && defined(COMPILER_TYPE_DECC))

	# define PFETCH 1

	# define CACHE_LINE_INTS    16
	# define CACHE_LINE_DOUBLES  8
	# define CACHE_LINE_COMPLEX  4

	# define prefetch_block(_array,_k)	asm("fetch_m   0(%a0)",(_array + _k));
	# define prefetch_p_doubles(_pd)	asm("ldt %f31,32(%a0)",(_pd));

/* Alpha system, Gnu C or C++ compiler */
#elif(defined(CPU_IS_ALFA) && defined(COMPILER_TYPE_GCC))

	# define PFETCH 1

	# define CACHE_LINE_INTS    16
	# define CACHE_LINE_DOUBLES  8
	# define CACHE_LINE_COMPLEX  4

	# define prefetch_p_doubles(_pd)	__asm__ volatile ("ldt $f31,32(%0)": :"r"(_pd));

/* Ultrasparc pre-v9 (no prefetch instruction per se), Sunpro compiler */
/* *** NEED CORRECT ASM SYNTAX HERE - JSP SAYS NATIVE C ALLOWS ONLY VERY LIMITED ASM INSERTION ***
#elif(defined(sparc) && !defined(__sparc_v9__) && defined(__SUNPRO_C))

# define SPARCV8 1
# define CACHE_LINE_INTS    16
# define CACHE_LINE_DOUBLES  8
# define CACHE_LINE_COMPLEX  4
# define prefetch_data(_array,_k) asm ("ldd [%0   ],%%f30": :"r"(_array+_k));
# define prefetch_p_doubles(_pd)  asm ("ldd [%0+64],%%f30": :"r"(_pd      ));
*/

/* Ultrasparc v8 and above, Sunpro compiler: */
#elif(defined(CPU_IS_SPARC) && defined(COMPILER_TYPE_SUNC))

  #if(defined(__SUN_PREFETCH))

	# include <sun_prefetch.h>

	# define PFETCH 1

	# define CACHE_LINE_INTS    16
	# define CACHE_LINE_DOUBLES  8
	# define CACHE_LINE_COMPLEX  4

	# define prefetch_data(_array,_k)	/* NO-OP! */
	# define prefetch_p_doubles(_pd)	sparc_prefetch_write_many(_pd);	/* Alternatives here are {read|write}_{once|many} */

  #else

	# define CACHE_LINE_INTS	0
	# define CACHE_LINE_DOUBLES	0
	# define CACHE_LINE_COMPLEX	0

	# define prefetch_data(_array,_k)	/* NO-OP! Use #pragma prefetch _array[_k] here??? */
	# define prefetch_p_doubles(_pd)	/* NO-OP! */ do { (void)(_pd); } while (0)

  #endif

/* Ultrasparc pre-v8 (no prefetch instruction per se), GNU C or C++ compiler */
#elif(defined(CPU_IS_SPARC) && defined(CPU_SUBTYPE_SPARC1) && defined(COMPILER_TYPE_GCC))

	# define PFETCH 1

	# define CACHE_LINE_INTS    16
	# define CACHE_LINE_DOUBLES  8
	# define CACHE_LINE_COMPLEX  4

	# define prefetch_p_doubles(_pd)	__asm__ volatile ("ldd [%0+16],%%f30": :"r"(_pd      ):"%f30");

/* Ultrasparc v8+, GNU C or C++ compiler */
#elif(defined(CPU_IS_SPARC) && defined(CPU_SUBTYPE_SPARC2) && defined(COMPILER_TYPE_GCC))

	# define PFETCH 1

	# define CACHE_LINE_INTS    16
	# define CACHE_LINE_DOUBLES  8
	# define CACHE_LINE_COMPLEX  4

	# define prefetch_p_doubles(_pd)	__asm__ volatile ("prefetch [%0+64],3": :"g"(_pd));

/*
	HP PA-RISC UNDER HPUX:

	From http://docs.hp.com/hpux/onlinedocs/2212/A-03-37relnotes.html:

	Gather/Scatter Prefetch pragma
	A pragma is now supported to prefetch specified cache lines. The behavior of this pragma is similar to +Odataprefetch but the
	prefetch pragma can access specific elements in indexed arrays that are stored in cache.  In addition, any valid lvalue can be used
	as an argument, but the intent of the pragma is to support array processing.

	Syntax:

	#pragma prefetch

	There can be only one argument per pragma.  The compiler generates instructions to prefetch the cache lines starting from the address
	given in the argument. The array element values prefetched must be valid.  Reading outside the boundaries of an array results
	in undefined behavior at runtime.

	Example:
	The function below will prefetch ia and b, but not a[ia[i]] when compiled with +O2 +Odataprefetch +DA2.0 (or +DA2.0W).

	void testprefc2(int n, double *a, int *ia, double *b)
	{
		for(int i=0; i<n, i++)
		{
			b[i]=a[ia[i]];
		}
	}

	Recording this routine as

	#define USER_SPECIFIED 30
	void testprefc2(int n, double *a, int *ia, double *b)
	{
		int dist=(int)USER_SPECIFIED;
		int nend=MAX(0,n_dist);	// so as not to read past the end of ia
		for(i=0;i<nend;i++)	// original loop is for (i=0;i<n;i++)
		{
			#pragma prefetch ia[i+4*dist]
			#pragma prefetch a[ia[i+dist]]
			b[i]=a[ia[i]];
		}
		// finish up last part with no prefetching
		for (int i=nend;i<n;i++)
		{
		   b[i]=a[ia[i]];
		}
	}

	The two pragma statements allow a[ia[i]] to be prefetched. Note that the compiler continues to unroll the loops as in the original code.

	There can be problems using the prefetch pragma when the kernel cannot allocate large pages.  Without large pages, there can be
	performance lost to Translation Lookaside Faults (TLB).  The optimal page size varies with different applications but 4MB page size is
	a good average.

	TLB faults occur when a particular page address does not reside in the TLB buffer.  This buffer contains the mapping of the
	virtual addresses to the absolute addresses of the pages recently fetched in the cache.  A TLB fault happens when a reference to a
	particular virtual page address cannot be translated to an absolute address in the buffer.

	Even when all the TLB and prefetch features are working, you are still limited by the memory bandwidth of the system.  The top
	bandwidth may be reduced by failing to load all the memory slots in some PA-RISC systems.  The memory controller depends on having all slots loaded to get the best bank interleaving.
*/
#elif(defined(CPU_IS_HPPA))

	# define CACHE_LINE_INTS	0
	# define CACHE_LINE_DOUBLES	0
	# define CACHE_LINE_COMPLEX	0

	# define prefetch_data(_array,_k)	/* NO-OP! Use #pragma prefetch _array[_k] here??? */
	# define prefetch_p_doubles(_pd)	/* NO-OP! */ do { (void)(_pd); } while (0)

/* IBM Power: */
#elif(defined(CPU_IS_POWER))

	# define CACHE_LINE_INTS	0
	# define CACHE_LINE_DOUBLES	0
	# define CACHE_LINE_COMPLEX	0

	# define prefetch_data(_array,_k)	/* NO-OP! Use #pragma prefetch _array[_k] here??? */
	# define prefetch_p_doubles(_pd)	/* NO-OP! */ do { (void)(_pd); } while (0)

#else

//	#error No platform-specific prefetch block found in prefetch.h!
	# define PFETCH 1

	# define CACHE_LINE_INTS     0
	# define CACHE_LINE_DOUBLES  0
	# define CACHE_LINE_COMPLEX  0

	# define prefetch_data(_array,_k)	/* */
	# define prefetch_p_doubles(_pd)	do { (void)(_pd); } while (0)

#endif

#endif	/* prefetch_h_included */

