/*******************************************************************************
*                                                                              *
*   (C) 1997-2020 by Ernst W. Mayer.                                           *
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

/* NOTE: Default should simply include this sourcefile as part of an Mlucas build.
To build the sieve factoring code in standalone mode, see the compile instructions below!
*/
#ifndef FACTOR_STANDALONE

	#include "factor.h"
	#include "align.h"

#else

	#include "Mlucas.h"

	/* Define FFT-related globals (declared in Mdata.h) */
	uint32 N2,NRT,NRT_BITS,NRTM1;
	uint32 NRADICES, RADIX_VEC[10];	/* RADIX_VEC[] stores sequence of complex FFT radices used.	*/
  #ifdef MULTITHREAD
	uint64 CORE_SET[MAX_CORES>>6];	// Bitmap for user-controlled affinity setting, as specified via the -cpu flag
  #endif

  #ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
	const uint32 mask02 = 0xfffffff8,
		br8[8]    = {0,4,1,5,2,6,3,7},	// length-8 index-scramble array for mapping from scalar-complex to AVX (re,re,re,re,im,im,im,im)
		brinv8[8] = {0,2,4,6,1,3,5,7};	// length-8 index-unscramble array: br[brinv[i]] = brinv[br[i]] = i .
  #endif
  #ifdef USE_SSE2
	const uint32 mask01 = 0xfffffffc,
		br4[4]    = {0,2,1,3};	// length-4 index-scramble array for mapping from scalar-complex (re,im,re,im) to SSE2 (re,re,im,im)
								// For length-4 this is its own inverse.
  #endif

	FILE *fp, *fq;
	const char OFILE     [] = "results.txt";
	/* Restart file name:
	This is set at runtime, based on either command-line -file flag argument
	or in the default case (i.e. -file {} not used) the exponent being processed: t[exponent].
	*/
	char RESTARTFILE[STR_MAX_LEN] = "";
	/* These should all be set to a valid (nonzero) value at the time the appropriate test is begun */
	uint32 TEST_TYPE		= 0;
	uint32 MODULUS_TYPE		= 0;
	uint32 TRANSFORM_TYPE	= 0;
	int INTERACT;

#endif


#undef RTIME
#undef CTIME

#ifdef MULTITHREAD
	#define RTIME	/* In multithreaded mode, need to use real (wall-clock) time */
#else
	#define CTIME	/* In single-thread mode, prefer cycle-based time because of its finer granularity */
#endif

// Use x86_64 inline-asm?
#undef YES_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
#endif

#ifdef MULTITHREAD

  #ifndef USE_PTHREAD
	#error Only Pthreads supported!
  #endif

	#include "threadpool.h"

	struct fac_thread_data_t{
		uint64*count;	// Put function-return value for unthreaded here
		int tid;	// Thread ID
		uint32 pass;
		uint64 interval_lo;
		uint64 interval_hi;
		double fbits_in_2p;
		uint32 nclear;
		uint32 sieve_len;
		uint32 p_last_small;	//largest odd prime appearing in the product; that is; the (nclear)th odd prime.
		uint32 nprime;			// #sieving primes (counting from 3)
		uint32 MAX_SIEVING_PRIME;
	#ifdef USE_AVX512
		uint32 *psmall;
	#endif
		uint8 *pdiff;
		uint32*startval;
		uint64*k_to_try;
		uint64*factor_k;		// List of found factors for each p gets updated (we hope) within
		uint32*nfactor;			// Here the '*' is to denote a writeable scalar
		uint32 findex;
		const char *pstring;
		uint64*p;
		uint64*two_p;
		uint64*q;
		uint64*q2;
		uint64*u64_arr;
		uint32 lenP;
		uint32 lenQ;
		uint32*kdeep;
		uint32*ndeep;
		uint64 countmask;
		uint32 CMASKBITS;
		uint32 incr;
		uint64 kstart;
		uint64*bit_map;
		uint64*bit_map2;
		double *tdiff;
		int    MODULUS_TYPE;
		const char  *VERSION;
		const char  *OFILE;
	};
#endif

// GPU and legacy-CPU version uses [p,k] (mod 60) classes, leading to 16 passes;
// in a manycore setting we will want to compile-time -DTF_CLASSES=4620 to keep all the threads busy:
#if defined(USE_GPU) || !defined(TF_CLASSES) || (TF_CLASSES != 4620)
	#define TF_CLASSES	60
	#define TF_PASSES	16
	#define TF_CLSHIFT	6
#else	// 4620 = 60*7*11 classes, leading to 960 passes: for CPU-builds targeting manycore architectures
	#define TF_CLASSES	4620
	#define TF_PASSES	960
	#define TF_CLSHIFT	12
#endif

#ifdef	USE_FMADD	// Need to add 100-bit modpow routines before enabling this for build of this file
	#warning USE_FMADD set in factor.c ... Using 100-bit FMA-based modpow.
#endif

#define SPOT_CHECK	0	// Enable periodic Spot-check (PRP or composite) of factor candidates

// printf character buffers - when using to print args in a single printf, need a separate buffer for each arg:
char char_buf0[STR_MAX_LEN], char_buf1[STR_MAX_LEN], char_buf2[STR_MAX_LEN];
//#define FAC_DEBUG
#ifdef FAC_DEBUG
	#warning FAC_DEBUG enabled for this build.
	/* Set k_targ to some known-factor k to debug a specific missed-factor case: */
	uint64 k_targ = 2*439880504ull;	// Factor k: for M(p) q = 2.k.p + 1; for F(n) q = k*2^(n+2) + 1 .
	uint32 pass_targ = 0xffffffff;	/* Init to a known-invalid value; if user specifies
	 								a known (valid) test factor via k_targ, pass_targ will
	 								be set to the associated (legitimate) value between 0 and TF_PASSES-1;
	 								(pass_targ < TF_PASSES) can subsequently be used as a quick
	 								test for whether an appropriate test factor has been set. */
	#define DBG_SIEVE	1
  #if DBG_SIEVE
	uint32 *startval_incr;
	uint32 i64_targ, bit_targ;
  #endif
#endif

// Oct 2015: GCD-associated self-tests provides a fair bit of added coverage of the mi64 library, so always include:
#ifdef INCLUDE_PM1
	#include "gcd_lehmer.h"
#endif

// Tables of Fermat-base-2 pseudoprimes needed by the small-primes-diff/2 computation:
#include "f2psp_3_5.h"

/* Factoring-only globals: */
int restart;

#ifdef FACTOR_STANDALONE
	/* In standalone, need to locally define some Mdata.h/Mlucas.h-declared globals
	which would normally be defined in Mlucas.c. */

	// System-related globals:
	uint32 SYSTEM_RAM = 0, MAX_RAM_USE = 90;	// Total usable main memory size, and max. amount of that to use per instance, in MB.
												// Default value of the latter is 90 (as in, use up to 90% of available RAM).

	char ESTRING[STR_MAX_LEN];	/* Exponent in string form */
	char PSTRING[STR_MAX_LEN];	/* [exponent | index] of [Mersenne | Fermat ] Number being tested in string form, typically
								 estring concatenated with other descriptors, e.g. strcat("M",estring) for Mersenne case
								*/
	uint64 PMIN;	/* minimum #bits allowed for FFT-based mul */
	uint64 PMAX;	/* maximum #bits allowed depends on max. FFT length allowed
					  and will be determined at runtime, via call to given_N_get_maxP(). */
	char cbuf[STR_MAX_LEN*2],cstr[STR_MAX_LEN];
	char in_line[STR_MAX_LEN];
	/* Declare a blank STATFILE string to ease program logic: */
	char STATFILE[] = "";
	// Prefixes corr. to #defines in Mdata.h
	const char*NUM_PREFIX[MODULUS_TYPE_MAX+1] = {"Unknown modulus type!","M","MM","F"};
	// These externs declared in platform.h:
	int MAX_THREADS = 0;		// Max. allowable No. of threads.
	int NTHREADS = 1;			// Actual No. of threads. If multithreading disabled, set = 1. Thus init = 1 and set

#endif

/* Adjust the # of sieving primes to reflect that modmul cost
goes up very quickly with increasing # of 64-bit words in q.

Note that the total allocated memory for the sieve is roughly (17 bytes * NUM_SIEVING_PRIME).

There are 203280220 odd primes < 2^32, thus to use all 32-bit primes needs 3.2-3.3 GB of memory, hence
such a maximal-depth sieve run should just fit into the working memory of a machine having 4GB RAM.
A FAC_DEBUG-enabled build needs 4 more bytes per sieving prime, which together with other needed memory
allocated will likely need slightly more than 4GB of working memory with such a maximal-depth sieve.
*/
#ifndef NUM_SIEVING_PRIME
	#ifdef NWORD	/* N-word P, N presumed to be > 1: */
		#ifdef FAC_DEBUG
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	2000000
		#endif
	#elif(defined(P4WORD))	/* 4-word P: */
		#ifdef FAC_DEBUG
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	1000000
		#endif
	#elif(defined(P3WORD))	/* 3-word P: */
		#ifdef FAC_DEBUG
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	500000
		#endif
	#elif(defined(P2WORD))	/* 2-word P: */
		#ifdef FAC_DEBUG
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	200000
		#endif
	#else	/* 1-word P limit is set by #of bits p can have so 120^2*p < 2^64: */
		#if defined(FAC_DEBUG) || defined(USE_GPU)
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	100000
		#endif
	#endif
#endif

#if NUM_SIEVING_PRIME > 203280220
	#error NUM_SIEVING_PRIME exceeds 32-bit-primes-set limit of 203280220!
#endif

/*...Known Mersenne prime exponents. This array must be null-terminated.	*/
const uint32 knowns[] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941
	,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593
	,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,0x0};

/*** Code to small candidate factors q=2*k*p+1 of Mersenne numbers M(p) = 2^p - 1.

	Uses the Left-to-Right binary method of exponentiation
	to quickly test whether [M(p)+1] mod q = 1, which is equivalent to
	testing whether M(p) mod q = 0, i.e. whether q divides M(p).

	We want the sieving code to be cache-friendly, which means minimizing
	the size of the sieve-related bit arrays in working memory. One way to
	achieve this is to break a single sequential sieving loop that tests
	all candidate factors of form q = 2*k*p + 1 satisfying some smoothness
	property (i.e. q not divsible by small primes below some predetermined
	threshold - note that it is easily seen that direct nonfactorial compositeness
	testing of each q is *not* cost-effective, since checking whether e.g.
	2^(q-1) == 1 mod q is as or more expensive than checking whether q divides
	2^p-1, i.e. whether 2^p == 1 mod q) into a number of smaller sieving
	steps, each of which checks only q's having k == 1 modulo some product
	of small primes. By eliminating k's which lie in residue classes such
	that the resulting q cannot possibly be prime we increase the density of
	candidate q's in our smaller "sievelets", and at the same time require
	smaller bit arrays to be held in working memory.

	Here, we use (p,k) mod 60 correlations to reduce the number of possible k's
	by roughly three-fourths. Since p prime it immediately follows that p must == 1 or
	(odd prime > 5) (mod 60), yielding 16 residue classes:
		p == 1,7,11,13,17,19,23,29,31,37,41,43,47,53,59 (mod 60).
	For each of these, there is an analogous residue classes on k: Specifically,
	a factor candidate q=2*k*p+1 can only be prime if the value of (k mod 60)
	appears in the appropriate (p mod 60) row of the following table:

	p%60	Acceptable values of k%60:
	--	-----------------------------------------------
	 1	00,03,08,11,15,20,23,24,35,36,39,44,48,51,56,59
	 7	00,05,08,09,12,17,20,24,29,32,33,44,45,48,53,57
	11	00,01,04,09,13,16,21,24,25,28,33,36,40,45,48,49
	13	00,03,08,11,12,15,20,23,27,32,35,36,47,48,51,56
	17	00,03,04,07,12,15,19,24,27,28,39,40,43,48,52,55
	19	00,05,09,12,17,20,21,24,29,32,36,41,44,45,56,57
	23	00,01,12,13,16,21,25,28,33,36,37,40,45,48,52,57
	29	00,04,07,12,15,16,19,24,27,31,36,39,40,51,52,55
	31	00,05,08,09,20,21,24,29,33,36,41,44,45,48,53,56
	37	00,03,08,12,15,20,23,24,27,32,35,39,44,47,48,59
	41	00,03,04,15,16,19,24,28,31,36,39,40,43,48,51,55
	43	00,05,08,12,17,20,21,32,33,36,41,45,48,53,56,57
	47	00,04,09,12,13,24,25,28,33,37,40,45,48,49,52,57
	49	00,11,12,15,20,24,27,32,35,36,39,44,47,51,56,59
	53	00,03,07,12,15,16,27,28,31,36,40,43,48,51,52,55
	59	00,01,04,09,12,16,21,24,25,36,37,40,45,49,52,57

Where did these come from?
Let pm := p%60, and km := k%60. For p prime > 5, we know p cannot be divisible by 2, 3 or 5.
Working modulo the product of small primes 2^2.3.5 = 60, we require GCD(p%60, 60) = 1, which allows
pm = 1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59. Now, since we also know that prime-exponent
Mersenne factor candidates must have form

q = 2*k*p + 1 = 2*(i*60 + km)*(j*60 + pm) + 1 == 2*km*pm + 1 (modulo 120) .

For any value of pm, the only possible value of km are those for which
GCD(2*km*pm + 1, 120) = 1, i.e. (2*km*pm + 1) is not divisible by 3 or 5.

Let's try pm = 1 as an example:

pm := p%60 = 1:

km := k%60 :            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59
2*km*pm+1  :            1  3  5  7  9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81 83 85 87 89 91 93 95 97 99.01.03.05.07.09.11.13.15.17.19
GCD(2*km*pm+1, 120)   : 1  3  5  1  3  1  1 15  1  1  3  1  5  3  1  1  3  5  1  3  1  1  5  1  1  3  1  5  3  1  1  3  5  1  3  1  1  5  1  1  3  1  5  3  1  1  3  5  1  3  1  1 15  1  1  3  1  5  3  1
Acceptable? (* = yes) : *        *     *  *     *  *     *        *  *        *     *  *     *  *     *        *  *        *     *  *     *  *     *        *  *        *     *  *     *  *     *        *

Thus, 32 k-residue classes survive the GCD criterion.
We further use quadratic residues: From the form of N = 2^p - 1, we have that
2^p == 1 (mod N). Multiplying by 2, we see that 2^(p+1) == 2 (mod N). Since p is odd
(primality of p is not crucial here, just oddness), the LHS is an even power of 2,
hence a perfect square, which implies that 2 is a QR mod N, and thus also a QR mod
any factor q of N. That immediately implies that q == +-1 (mod 8), i.e. that any
factor q must be of the form 8*n +- 1. Thus our set of acceptable km's is cut in half:

km := k%60 :            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59
2*km*pm+1 % 8         :+1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1 +1 +3 -3 -1
Acceptable? (* = yes) : *        *              *        *           *              *        *  *                                *  *        *              *           *        *              *        *

Thus, the divisibility-by 3,5 test leaves 32 acceptable residue classes,
precisely half of which are eliminated by the QR criterion.
An analogous 32/16 pattern occurs for the other 15 possible values of p%60.

NOTE that since the smallest possible k value = 1, this corresponds to the zero bit of
our sieve, i.e. k%60 = 0 is aliased to k%60 = 60 - that's why the k%60 == 0 pass is
done last among the 16 k%60 passes.

Mar 2015: Added extension of above scheme to (mod 4620), i.e. 60*7*11, thus adding primes 7,11 to the earlier pair 3,5.
There are now 16*77*(1 - 1/7)*(1 - 1/11) = 16*60 = 960 possible residue classes for the p's, and analogously for the k's.
Unlike for (mod 60), use simple utility functions to manage these, rather than a precomputed giant 960 x 960 table.

                        *     *     *     *     *

   REVISION CONTROL:

   01/04/2003:	First C version. Limit = 2^65. ASM mul macros only for Alpha.

   06/--/2003:	Nonportable (Alpha and Itanium only) C version. Limit = 2^ 96.

   12/31/2003:	Nonportable (Alpha and Itanium only) C version.
				Limit = 2^128, p up to 2^114 via -DP2WORD compile-time #define.

   06/01/2005:	Nonportable (Alpha and Itanium only) C version.
				Limit = 2^192, p up to 2^128 via -DP3WORD compile-time #define.
				Can also handle case where p (actually, p*120) << 2^128 but q > 2^128
				via combination of -DP2WORD and -DQ3WORD compile-time #defines.
				This invokes a fast version of the 8-operand modular powering routine
				which assumes that only the lower 128 bits of each octet of q's differ.

	05/15/2006:	***First Public Release***
   				Portable C version with ASM integer-multiply macros for a decent variety of 64-bit
 				architectures (Alpha, Itanium, AMD64, PowerPC, IBM Power64).
 				Factor limits of 2^128/192/256 via PxWORD compile-time #define with x = 2,3,4.

	01/15/2006:	Arbitrary-length integer arithmetic (but optimal for "small" integer lengths,
                say on the order of ten 64-bit words) via enhancements to the mi64 library.
                Invoke by defining NWORD, e.g. via -DNWORD.

	2009-2019	Added support for x86 SSE2 floating-double SIMD hardware; improved inline ASM for x86 uint64 mode.

	2013		Added support for AVX to float-based (78-bit modmul) routines.

	2015		Added support for 4620 ( = 60*7*11) p-mod classes, yielding 960 passes. This saves ~25%
				(precisely: 17/60*100%) effort versus legacy mod-60, and lends itself to manycore/manythread
				CPU (not GPU - there we use a different threading model) implemeentation.

				Added support for nVidia CUDA, geared toward 32-bit-int and floating-double (esp. FMA-arithmetic)
				capability of Fermi+ families of GPUs. Also added AVX2-based modmul routines based on capability
				of FMA-double math to exactly compute 106-bit product of 53-bit signed inputs via the sequence

					[input doubles x,y]
					double hi = x*y;
					double lo = FMA(x,y, -hi);

				This still needs normalization to properly split the resultant significant bits across lo,hi,
				e.g. in a base-2^52 signed-int multiword context:

					hh  = DNINT(hi*TWO52FLINV);		// This part remains in hi...
					cy  = FMA(hh ,-TWO52FLOAT,hi);	// hi - hh*2^52 is 'backward carry' from hi into lo, needed for proper base-normalization
					lo += cy;
					hi  = hh;

				Thus the cost of each 100+ - bit product is 1 ADD, 2 MUL, 2 FMA, 1 DNINT. (We assume/hope the negations are free).

                        *     *     *     *     *

	COMPILING THE PROGRAM IN STANDALONE MODE: Below are the best compile options for various
	platforms on which the code has been run successfully. If your particular hardware/software
	is not represented, please let me know what kinds of results you get with your compiler.

	The compile sequences assume that the sourcefile is called factor.c
	and produce an executable named "factor".

  *** Alpha: ***

	NOTES: Alpha 21264 (ev6) and above should give by far the best per-cycle performance
	of all the Alphas, because it fully pipelines the MULL64 and MULH64 hardware integer
	multiply instructions, i.e. needs effectively just 2 cycles to get a full-length
	128-bit product of unsigned 64-bit integers.
	The 21164 (ev5) only partially pipelines these - it can start a new integer mul instruction
	every 8 cycles, so needs effectively 16 cycles for a 128-bit product.
	The 21064 (ev4) doesn't pipeline integer muls at all, so will need 30-40 cycles for such
	a product. Although integer muls are the rate-limiting operation on the slower machines,
	they are not on the 21264, so in practice the 21264 gives ~4-5x better per-cycle factoring
	performance than the 21164 and ~12-15x better than the 21064.

  * (1A) Digital/Compaq/HP C compiler (cc) for TruUnix V5.0+: if you're not sure what type of Alpha
	CPU you have, substitute 'host' for the architecture-specific (ev4/generic/ev6) flag below.

	21064/ev4/ev45:			cc -o factor -DFACTOR_STANDALONE -O5 -fast -arch ev4     -tune ev4     factor.c -lm
	21164/ev5/ev56/generic:	cc -o factor -DFACTOR_STANDALONE -O5 -fast -arch generic -tune generic factor.c -lm
	21164/ev6/ev67/ev68:	cc -o factor -DFACTOR_STANDALONE -O5 -fast -arch ev6     -tune ev6     factor.c -lm

  * (1B) Digital/Compaq/HP C compiler (ccc) for Linux: same as for Unix, but replace 'cc' with 'ccc'.

  * (1C) Digital/Compaq/HP C compiler (cc) for VMS (Since this is integer-dominated code I'm not sure
	if all these are necessary, but some of the floating-point ones apparently are, and the others certainly won't hurt):

	*** NEED TO FIGURE OUT HOW TO INCLUDE FACTOR_STANDALONE #DEFINE UNDER VMS! ***
	cc/decc/lis/mach/nodebug/float=ieee_float/ieee_mode=fast/arch=host/assume=accuracy_sensitive/opt=(inline=speed,level=5,tune=host,pipeline,unroll=1) factor.c

  * (1D) Gnu C for Alpha/Linux: I'm not sure what the best compile options are. Based on past
	experience, I strongly suggest using either the native Digital/Compaq/HP C compiler,
	or a prebuilt binary based on it, as gcc generally gives suboptimal performance on most platforms.


  *** IA64 ***

  * (1A) Intel C (icc) for IA64: please make sure to use v8.0 or later!

	icc -o factor -static -DFACTOR_STANDALONE -O3 factor.c


  * (2B) HP C for IA64/HPUX:

	cc -o factor -DFACTOR_STANDALONE +O3 factor.c -lm

	(Note that you'll see a bazillion ignorable warnings of the form
	 Warning 1001: "factor.c", line 6574 # Conversion from 'int' to
	 '__fpreg' performed using intermediate conversion to __float80. )


  *** Athlon-64 (a.k.a. Opteron) ***

	gcc -o factor -DFACTOR_STANDALONE -m64 -O3 factor.c
*/

/* Various compile-time #defines: (define via -D[option], unless noted otherwise) -
	Subordinate #defines are indicated via indentation (e.g. PIPELINE_MUL192 is only
	meaningful if P3WORD is in effect) :

	FACTOR_STANDALONE - build as a standalone program, rather than a callable
						factor() subroutine.

	FAC_DEBUG - define at compile time to enable debugging diagnostic prints
				and assertions of the form ASSERT.

	DBG_SIEVE - define = 1 to debug the sieve without actually spending time
				doing trial divides. NOTE: Requires FAC_DEBUG to be def'd.

	NOBRANCH - define to invoke branchless versions of key code segments.

	QUIT_WHEN_FACTOR_FOUND - Default is to finish all 16 passes to the specified depth,
	i.e. to find *all* factors below the factoring limit. To instead stop immediately
	if a factor is found, define this one.

***	USE_FLOAT - Only relevant if no PxWORD with x > 1 is defined,
	i.e. when dealing with 1-word p's (more specifically, p such that 120*p < 2^64.)
	Define this flag to use specialized floating-double-based code (cf. twopmodq80.c)
	to do the bulk of each modular powering. This limits factor candidates to < 2^78.
	If this flag is invoked, the user may not invoke the USE_65BIT or USE_95BIT flags.

	USE_FLOATING_MULH64 - currently unsupported

	USE_65BIT - Only relevant if no PxWORD with x > 1 is defined,
	i.e. when dealing with 1-word p's (more spedifically, p such that 120*p < 2^64.)
	Default is to use generic 96-bit factor routines for all q's > 64 bits.
	Define this flag to use specialized 65-bit code to handle q's in [2^64, 2^65).

***CURRENTLY UNSUPPORTED***	USE_95BIT - Only relevant if no PxWORD with x > 1 is defined,
	i.e. when dealing with 1-word p's (more specifically, p such that 120*p < 2^64.)
	Default is to use generic 96-bit factor routines for all q's > 64 bits.
	Define this flag to use specialized 95-bit code to handle q's in [2^65, 2^95).

	USE_128x96 - Only relevant if no PxWORD with x > 2 is defined, and only for q's in [2^64, 2^96].
	Invoke to replace calls to the full 128-bit modmul routines with ones to
	special 96-bit or 128/96-bit-hybrid hyroutines when the operands allow, specifically q < 2^96.
	There are 3 currently supported values:

		0 (or undef'd) - use  fully 128-bit routines for q's in [2^64, 2^96]

		1 - use strictly   96-bit routines for q's in [2^64, 2^96]

		2 - use hybrid 128_96-bit routines for q's in [2^64, 2^96]

	P1WORD - p small enough such that p*120 < 2^64, factor limit q < 2^96

	P2WORD - factor limit q < 2^128, i.e. q needs 2 full 64-bit words of storage.
	Also needed if p*120 is sufficiently close to 2^64 that the assumption of the
	96-bit modmul routines that the high 32 bits of q_j = 2*(k%60 + 60*j)*p + 1
	change only rarely with increasing j is no longer tenable.

	P3WORD - factor limit q < 2^192, i.e. q needs 3 full 64-bit words of storage.

		PIPELINE_MUL192 - when def'd != 0 while building twopmodq192.c, uses pipelined versions of 192-bit MUL macros.

	P4WORD - factor limit q < 2^256, i.e. q needs 4 full 64-bit words of storage.

		To-Do: Add support for PIPELINE_MUL256

	NWORD - Arbitrary-length p and q, only restriction is that (as for all other size ranges) kmax < 2^64 .
*/

/*********************************************************************************************************************************/

/* NOTE: Exponents > 64 bits *require* standalone-mode build: */

#ifndef FACTOR_STANDALONE

  #ifdef FAC_DEBUG
	#error FAC_DEBUG only permitted in standalone mode!
  #endif

int factor(char *pstring, double bmin, double bmax)
{
	ASSERT(HERE, 0, "TF currently not supported as part of Mlucas, only via standalone Mfactor build - please delete any .o files and retry USING 'makemake.sh mfac' from Mluas dir above /src.");
	return 1;
}

#else	// Standalone build - we leave some FACTOR_STANDALONE-wrapped stuff inside code below for future Mlucas integration

int main(int argc, char *argv[])
{
	static int first_entry = TRUE;

  #ifdef FACTOR_STANDALONE

	/*...file pointer	*/
	FILE *fp, *fq;
	char stFlag[STR_MAX_LEN];
	/* Allocate storage for any needed Globals declared in Mdata.h
	(in non-standalone mode these are instead defined in Mlucas.c): */
	int MODULUS_TYPE   = 0;
	char pstring[STR_MAX_LEN] = "";
	/*...program version with patch suffix... */
	const char VERSION[] = "3.0x";			/* <--- a suffix of x, y, or z indicates a beta code. */
	const char OFILE  [] = "results.txt";	/* ASCII logfile containing factors found and/or
											final factoring-run result ONLY for each assignment */
	double bmin = 0.0, bmax = 0.0;	/* store log2 of (min|max) factoring bound */
  #endif

  #ifdef MULTITHREAD

	static struct fac_thread_data_t *tdat = 0x0;
	int thr_id;
	// For Threadpool-based dispatch:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)PerPass_tfSieve, NULL, 0x0};

  #endif

	/* Make these large enough to hold the max. p,q supported by the software (currently,
	256 bits) then use as many of the subfields as needed, depending on the values of P,QWORDS.
	*/
	uint32	findex = 0, nbits_in_p = 0, nbits_in_q = 0, lenP = 0, lenQ = 0, bits_in_pq2 = 0;
	double fbits_in_2p = 0;
	uint64	*factor_ptmp = 0x0, *factor_k = 0x0;	/* Use this array to store the factor k's (assumed < 2^64) of any factors found. */
	uint64	*p = 0x0, *two_p = 0x0, *p2NC = 0x0, *q = 0x0, *q2 = 0x0, *k_to_try = 0x0;

  #ifndef TRYQ
	#define TRYQ	1
  #elif TRYQ < 1
	#undef TRYQ
	#define TRYQ	1
  #endif

  #ifdef USE_AVX512
	#ifndef USE_FLOAT
		#warning USE_AVX512 only meaningful USE_FLOAT also defined at compile time - setting this #define.
		#define USE_FLOAT
	#endif
	#define MAX_TRYQ	64
  #elif defined(USE_AVX)
	#define MAX_TRYQ	16
  #elif defined(USE_SSE2) || defined(X64_ASM)
	#define MAX_TRYQ	 8
  #else
	#define MAX_TRYQ	 4
  #endif

  #if (TRYQ > MAX_TRYQ) && !defined(USE_ARM_V8_SIMD)	// Ignore on ARMv8, since no TF support there anyway
	#error TRYQ exceeds MAX_TRYQ for this build mode!
  #endif

	int nargs,itmp;
	/* pdsum stores the sums needed for the base (%30 == 0) candidate of each length-30 interval;
	pdiff stores the diffs/2 of these absolute offsets, circularly shifted to the right one place,
	since to update we need diff[j] = curr_p(current) - curr_p(previous) = (pdsum[j] - pdsum[j-1])/2,
	i.e. pdiff_8[i] = (pdsum_8[i]-pdsum_8[i-1])/2, i computed using mod-8 arithmetic and pdsum_8[] using mod-30 arithmetic:
	*/
	uint32 itmp32, pdiff_8[8] = {2,1,2,1,2,3,1,3}, pdsum_8[8] = { 0, 2, 6, 8,12,18,20,26};
	uint64 itmp64;
	uint32 pmodNC,kmodNC,incr[TF_PASSES];
	uint64 two64modp;
	uint32 bit,curr_p,i,ihi,m,ncopies,qmod8,regs_todo;
	uint32 l,i64,nfactor,word;
	uint32 nprime = NUM_SIEVING_PRIME;	/* Don't make this const, in case we need to reduce it to satisfy MAX_SIEVING_PRIME < (q_min = 2*p+1) */
	uint32 f2psp_idx = 0;	/* Index to next-expected Fermat base-2 pseudoprime in the precomputed table */

	/* LEN is the number of 64-bit words in our sieving bitmap. Successive
	bits of the sieve bitmap correpond to successive factor k-values, with bit 0
	corresponding to the smallest possible factor, k = 1 (note the unit-offset!)
	In practice, this full-sized sieve is split into 16 smaller (i.e. hopefully
	cache-sized) "sievelets", each of which contains bits corresponding to
	successive k's from one of the 16 allowable (k mod TF_CLASSES) families for the given
	exponent p (more specifically, the given (p mod TF_CLASSES) value.) Both full-length
	sieve and 1/16-th length sievelets share the property that each time we run
	through the bits of onesuch we run through 64*LEN worth of k's.
	*/
	const uint32 nclear=6, len = 3*5*7*11*13*17;	// LEN = product of first NCLEAR odd primes [= 255255],
	const uint32 p_last_small = 17;					// p_last_small = largest prime appearing in the product [= 17].
	uint32 prime[] = {3,5,7,11,13,17};	// Also need the first [nclear] odd primes - since 'const' is treated as a read-only flag
								// on *variables*, use [] instead of [nclear] to avoid compiler 'variable-length array init' errors.
  #if TF_CLASSES == 60
	const uint32 bit_len = (len << TF_CLSHIFT)/TF_CLASSES; 	// 255255*64  /  60 = 272272: Number of bits in each of the  16 mod-  60 sievelets
  #else	// 4620 classes:
	const uint32 bit_len = (len << TF_CLSHIFT)/TF_CLASSES;	// 255255*64^2/4620 = 226304: Number of bits in each of the 960 mod-4620 sievelets
  #endif
	//   bits cleared of multiples of 3,5,7,11,13, 17 and q mod 8 = +-1 are here:
	uint64 *temp_late = 0x0;		/* Even though we know how large to make this, it's only needed
									for data inits, so we calloc it at runtime and free it later. */
	uint32 on_bits = 0;
	uint64 *bit_map, *bit_map2, *bit_atlas = 0x0;
	uint32 pass = 0xffffffff, passmin = 0, passnow = 0, passmax = TF_PASSES-1;
	uint64 count = 0,countmask,j,k,kmin = 0,kmax = 0,know = 0,kplus = 0;
	uint32 CMASKBITS;	// This is set at runtime based on the operand sizes, but treat as read-only subsequently.

	/* If restart file found, use these to store bmin/max, kmin/max, passmin/max
	data contained therein (e.g. for comparison with command-line args and determination
	of run status). We don't init them here since if a valid restart file is found,
	all these values should be read from it - that way if we get an "uninitialized"
	warning from the compiler, we know we've done something wrong in the code logic.
	*/
	int incomplete_run = FALSE;
	int curr_line;
	double bmin_file, bmax_file;
	uint64 tf_passes_file, kmin_file, know_file = 0, kmax_file = 0;
	uint32 passmin_file = 0, passnow_file = 0, passmax_file = 0;

	uint64 interval_lo,interval_now,interval_hi,ninterval;

// This stuff is for the small-primes sieve:
	uint32 max_diff;
  #ifdef USE_AVX512	// Use vector-int math and gather-load/scatter-store to accelerate the bit-clearing
	uint32 *psmall;
  #endif
	uint8 *pdiff;	/* Compact table storing the (difference/2) between adjacent odd primes.
							http://mathworld.wolfram.com/PrimeGaps.html shows the first >256-gap at ~387 million
							and the first >512-gap at ~300 billion, so using the half-of-even-gap trick makes
							the difference between an algo which is safe for all 32-bit primes and one with a
							limit of < 2^30. */
	uint32 *startval, *pinv, *kdeep = 0x0, ndeep = 0;
	uint32 MAX_SIEVING_PRIME = 0;
	uint64 *u64_arr = 0x0;	/* generic array allowing us to hook into the mi64 routines */

  #ifdef P1WORD
	double twop_float = 0,fqlo,fqhi;
	uint128 p128,q128,t128;	// Despite the naming, these are needed for nominal 1-word runs with moduli exceeding 64 bits
  #endif
  #ifdef P3WORD
	uint192 p192,q192,t192;
  #ifdef USE_FLOAT
	uint256 x256;	// Needed to hold result of twopmodq200_8WORD_DOUBLE
  #endif
  #endif
  #ifdef P4WORD
	uint256 p256,q256,t256;
  #endif

	/*...time-related stuff	*/
  #ifdef CTIME
	clock_t clock1, clock2;
  #else	// Multithreaded needs wall-clock, not CPU time:
	double clock1, clock2;	// Jun 2014: Switched to getRealTime() code
  #endif
	double td, tdiff;

  #if TEST_TRIALDIV
	double citer;
  #endif
	char *char_addr;

/* Set == 1 to test the trial-div stuff: */
  #define	TEST_TRIALDIV	0
  #if TEST_TRIALDIV
	#define MAX_ARRAY_DIM 10000
	uint32	vec_len = MAX_ARRAY_DIM;
	uint64*	xvec = (uint64 *)calloc(MAX_ARRAY_DIM, sizeof(uint64));
	uint32 tryq[8];
  #endif

  #ifdef USE_GPU
	cudaError_t cudaError = cudaGetLastError();	// Call this to reset error flag to 0
	if(cudaError != cudaSuccess)
	{
		printf("ERROR: cudaGetLastError() returned %d: %s\n", cudaError, cudaGetErrorString(cudaError));
		ASSERT(HERE, 0, "factor.c : GPU-side error detected!");
	}
  #endif

  #ifdef macintosh
	argc = ccommand(&argv);	/* Macintosh CW */
  #endif

/* Allocate factor_k array and align on 16-byte boundary: */
	factor_ptmp = ALLOC_UINT64(factor_ptmp, 24);
	factor_k = ALIGN_UINT64(factor_ptmp);	factor_ptmp = 0x0;
	ASSERT(HERE, ((uint64)factor_k & 0x3f) == 0, "factor_k not 64-byte aligned!");

/*...initialize logicals and factoring parameters...	*/
	restart = FALSE;

  #ifdef FACTOR_STANDALONE
	host_init();
  #endif

/***********************************************************************/
/******* In standalone mode, process any command-line arguments: *******/
/***********************************************************************/
/*
********** Mfactor command line options: **********

REQUIRED:
		* One (and ONLY one) of -m|mm|f, followed by a valid numerical exponent;
		* One (and ONLY one) of -bmax|kmax, unless it's a restart, i.e. a valid checkpoint file
			for the number in question exists. Iff -bmax|kmax specified, an optional lower-bound
			argument -bmin|kmin may also be specified, which must not exceed the upper bound.
		* If neither -bmax|kmax specified, it is assumed that a valid checkpoint file
			for the number in question exists. The data in this file will indicate either
			an as-yet-uncompleted factoring run for the number in question (in which case
			the run is resumed at the point at which it left off), or a completed run. In
			the latter instance, if a -kplus argument was specified on the command line,
			the k-bounds of the previous completed run are incremented and a new run with
			k-bounds [kmax_previous, kmax_previous + kplus] is begun. If -kplus is specified
			but the restart-file data indicate an as-yet-uncompleted run, a warning is issued,
			the -kplus argument ignored, and the incomplete run is resumed.

Others are optional and in some cases mutually exclusive:

	-h          Prints this help menu and exits

	-m [int]    Trial-factor the Mersenne number M(int) = 2^int - 1, with int < 2^MAX_BITS_P.

	-mm [int]   Trial-factor the double-Mersenne number M(M(int)) = 2^(2^int) - 1, with M(int) < 2^MAX_BITS_P.

	-f [int]    Trial-factor the Fermat number F(int) = 2^(2^int) + 1, with int <= %u.\n",MAX_BITS_P.
			NOTE:
				* Fermat number Trial-factoring not currently supported (future release.)

	-bmin [float] Log2(min factor to try), in floating form (>= 0, default = 0).
	-bmax [float] Log2(max factor to try), in floating form ( < 64*NWORDS).
			NOTES:
				* If -bmin/bmax used to set lower/upper bounds for factoring, -kmin/kmax disallowed.
				* bmin/bmax form of bounds-setting only allowed for single-word-p case, since
				  multiword p may cause float approximations to p, 2*p etc to overflow.

	-kmin  [int] Lowest  factor K value to be tried in each pass ( > 0, default = 1).
	-kmax  [int] Highest factor K value to be tried in each pass ( < 2^64).
			NOTE:
				* If -kmin/kmax used to set lower/upper bounds for factoring, -bmin/bmax disallowed.

	-kplus [int] Added   factor K value to be tried in each pass ( < 2^64), for an exponent
				for which one or more previous shallower factoring runs have already been done
				(specifically, a checkpoint file for a previous run exists.)
			NOTES:
				* If -bmin/bmax or -kmin/kmax used to set bounds for factoring, -kplus disallowed (and v.v.)
				* If -kmin|kmax from a previous run of the number in question found
				 in a checkpoint file, that old kmax serves as kmin for the new run
				 and (old kmax) + (kplus) serves as kmax for the new run.

	-passmin [int]  Maximum factoring pass for the run (0-TF_PASSES-1, default =  0).
	-passmax [int]  Maximum factoring pass for the run (0-TF_PASSES-1, default = 15).
			NOTE:
				* If passmin|max from a previous run of the number in question found
				 in a checkpoint file and those pass bounds conflict with the ones
				 given via the command line, an error message is printed and the
				 run aborted. This is done as a precaution against inadvertently
				 skipping a range of trial-factoring bounds in a multipart series
				 of partial factoring runs. In this event the user shhould carefully
				 compare the checkpoint file(s) for the number in question they have
				 saved from previous runs with their current command line and modify
				 one or the other so as to remove any pass-range conflicts.

	-q [int]    A known factor for the number (only used if FAC_DEBUG def'd).
*/
  #ifdef FACTOR_STANDALONE

	nargs = 1;
	if(!argv[nargs])
		goto MFACTOR_HELP;
	while(argv[nargs])
	{
		strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

		if(stFlag[0] != '-')
		{
			fprintf(stderr,"*** ERROR: Illegal command-line option %s\n", stFlag);
			fprintf(stderr,"*** All command-line options must be of form -{flag} [argument]\n\n");
			goto MFACTOR_HELP;
		}

		if(STREQ(stFlag, "-h"))
		{
			goto MFACTOR_HELP;
		}
		/* Type of number to be trial-factored: */
		else if(STREQ(stFlag, "-m"))	/* Mersenne */
		{
			strncpy(pstring, argv[nargs++], STR_MAX_LEN);
			MODULUS_TYPE = MODULUS_TYPE_MERSENNE;
		}
		else if(STREQ(stFlag, "-mm"))	/* Double-Mersenne */
		{
			strncpy(pstring, argv[nargs++], STR_MAX_LEN);
			MODULUS_TYPE = MODULUS_TYPE_MERSMERS;
		}
		else if(STREQ(stFlag, "-f"))	/* Fermat */
		{
			strncpy(pstring, argv[nargs++], STR_MAX_LEN);
			MODULUS_TYPE = MODULUS_TYPE_FERMAT;
		}

		/* Factor bounds, in log2(qmin/qmax) (floating double) form: */
		else if(STREQ(stFlag, "-bmin"))
		{
			if(kmin || kmax || kplus)
			{
				fprintf(stderr,"*** ERROR: If -kmin/kmax or -kplus used to set bounds for factoring, -bmin [and -bmax] disallowed.\n");
				goto MFACTOR_HELP;
			}
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			bmin = convert_base10_char_double(stFlag);
		  #ifdef FAC_DEBUG
			printf("bmin = %lf\n", bmin);
		  #endif
		}
		else if(STREQ(stFlag, "-bmax"))
		{
			if(kmin || kmax || kplus)
			{
				fprintf(stderr,"*** ERROR: If -kmin/kmax or -kplus used to set bounds for factoring, -bmax [and -bmin] disallowed.\n");
				goto MFACTOR_HELP;
			}
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			bmax = convert_base10_char_double(stFlag);
		  #ifdef FAC_DEBUG
			printf("bmax = %lf\n", bmax);
		  #endif
		}

		/* Factor bounds, in kmin/kmax (uint64) form: */
		else if(STREQ(stFlag, "-kmin"))
		{
			if(bmin || bmax || kplus)
			{
				fprintf(stderr,"*** ERROR: If -bmin/bmax or -kplus used to set bounds for factoring, -kmin/kmax disallowed.\n");
				goto MFACTOR_HELP;
			}
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			kmin = convert_base10_char_uint64(stFlag);
		}
		else if(STREQ(stFlag, "-kmax"))
		{
			if(bmin || bmax || kplus)
			{
				fprintf(stderr,"*** ERROR: If -bmin/bmax or -kplus used to set bounds for factoring, -kmin/kmax disallowed.\n");
				goto MFACTOR_HELP;
			}
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			kmax = convert_base10_char_uint64(stFlag);
		}

		else if(STREQ(stFlag, "-kplus"))
		{
			if(bmin || bmax || kmin || kmax)
			{
				fprintf(stderr,"*** ERROR: If -bmin/bmax or -kmin/kmax used to set bounds for factoring, -kplus disallowed.\n");
				goto MFACTOR_HELP;
			}
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			kplus = convert_base10_char_uint64(stFlag);
		}

		/* Pass bounds: */
		else if(STREQ(stFlag, "-passmin"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			passmin = (uint32)convert_base10_char_uint64(stFlag);
			ASSERT(HERE, passmin < TF_PASSES,"factor.c: passmin < TF_PASSES");
		}
		else if(STREQ(stFlag, "-passmax"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			passmax = (uint32)convert_base10_char_uint64(stFlag);
			ASSERT(HERE, passmax < TF_PASSES,"factor.c: passmax < TF_PASSES");
			ASSERT(HERE, passmax >= passmin       ,"factor.c: passmax >= passmin");
		}

		// Number of threads to use?
		else if(STREQ(stFlag, "-nthread"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
		#ifndef MULTITHREAD
			fprintf(stderr,"Multithreading not enabled; ignoring -nthread argument.\n");
			NTHREADS = 1;
		#else
			itmp = (uint32)convert_base10_char_uint64(stFlag);
			if(itmp > TF_PASSES) {
				fprintf(stderr, "factor.c: Specifed #nthreads [%u] > TF_PASSES [%u] ... using just %u threads.\n",itmp,TF_PASSES,TF_PASSES);
				NTHREADS = TF_PASSES;
			} else if(itmp < 1) {
				fprintf(stderr, "factor.c: Specifed #nthreads [%u] < minimum of 1 ... running 1-threaded.\n",itmp);
				NTHREADS = 1;
			} else {
				NTHREADS = itmp;
			}
		  #ifdef NWORD
			ASSERT(HERE, NTHREADS == 1, "Arbitrary-precision build currently only supports single-threaded runs!");
		  #endif
		#endif
		}
		else	// Come again?
		{
			fprintf(stderr,"*** ERROR: Unrecognized command-line option %s\n", stFlag);
			fprintf(stderr,"*** All command-line options must be of form -{flag} [argument]\n\n");
			goto MFACTOR_HELP;
		}
	}

  #else

	/* If non-standalone mode, make sure statfile name is non-empty: */
	ASSERT(HERE, STRNEQ(STATFILE, ""), "STATFILE string empty");
	fp = mlucas_fopen(STATFILE, "a");
	if(!fp) {
		fprintf(stderr,"ERROR: Unable to open statfile %s for writing.\n",STATFILE);
		ASSERT(HERE, 0,"0");
	} else {
		fclose(fp); fp = 0x0;
	}

  #endif	/* #ifdef FACTOR_STANDALONE */

	// One-time allocs and inits:
	if(first_entry)
	{
		first_entry = FALSE;
	#ifndef MULTITHREAD
		#warning Building factor.c in unthreaded (i.e. single-main-thread) mode.
		ASSERT(HERE, NTHREADS == 1, "NTHREADS must == 1 in single-threaded mode!");
		k_to_try = (uint64 *)calloc(TRYQ * NTHREADS, sizeof(uint64));
	#else
		MAX_THREADS = get_num_cores();
		ASSERT(HERE, MAX_THREADS > 0, "Illegal #Cores value stored in MAX_THREADS");
		ASSERT(HERE, MAX_THREADS <= MAX_CORES,"MAX_THREADS exceeds the MAX_CORES setting in Mdata.h .");

		if(!NTHREADS) {
			NTHREADS = 1;
			fprintf(stderr,"No CPU set or threadcount specified ... running single-threaded.\n");
			// Use the same affinity-setting code here as for the -cpu option, but simply for cores [0:NTHREADS-1]:
		} else if(NTHREADS > MAX_CORES) {
			sprintf(cbuf,"FATAL: NTHREADS = %d exceeds the MAX_CORES setting in Mdata.h = %d\n", NTHREADS, MAX_CORES);
			ASSERT(HERE, 0, cbuf);
		} else {	// In timing-test mode, allow #threads > #cores
			if(NTHREADS > MAX_THREADS) {
				fprintf(stderr,"WARN: NTHREADS = %d exceeds number of cores = %d\n", NTHREADS, MAX_THREADS);
			}
			fprintf(stderr,"NTHREADS = %d\n", NTHREADS);
		}
		sprintf(cbuf,"0:%d",NTHREADS-1);
		parseAffinityString(cbuf);
		k_to_try = (uint64 *)calloc(TRYQ * NTHREADS, sizeof(uint64));

		// Up to TF_PASSES work units (perhaps fewer if a restart) get done by a pool of NTHREADS threads.  Yypically have
		// NTHREADS <= TF_PASSES, i.e. pool threads get reassigned a fresh work unit as they complete their current one.
		// We do the || work in discrete (endpoint-sync'ed) 'waves' of NTHREADS each, thus need only that many thread-data allocs:
		if(tdat) {
			free((void *)tdat); tdat = 0x0;	// Not sure if we might ever have occasion to realloc here, but easy enough to set up for it
		}
		tdat = (struct fac_thread_data_t *)calloc(NTHREADS, sizeof(struct fac_thread_data_t));

		// Alloc threadpool of NTHREADS threads, which will concurrently/asynchronally
		// do TF_PASSES 'work units' (factoring passes for various (k mod TF_CLASSES) k-classes:
		main_work_units = 0;
		pool_work_units = NTHREADS;
		ASSERT(HERE, 0x0 != (tpool = threadpool_init(NTHREADS, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
		printf("Factor.c: Init threadpool of %d threads\n", NTHREADS);

		// Apr 2015: Init-calls to any inline-asm-using modpow functions:
		int thr_id = -1;
		twopmodq96_q4(0ull, 0ull, 0ull, 0ull, 0ull, NTHREADS, thr_id);
	  #ifdef USE_SSE2
		twopmodq78_3WORD_DOUBLE_q2 (0ull, 0ull,0ull, NTHREADS, thr_id);
		twopmodq78_3WORD_DOUBLE_q4 (0ull, 0ull,0ull,0ull,0ull, NTHREADS, thr_id);
		twopmodq78_3WORD_DOUBLE_q8 (0ull, k_to_try, NTHREADS, thr_id);
	  #endif
	  #ifdef USE_AVX
		twopmodq78_3WORD_DOUBLE_q16(0ull, k_to_try, NTHREADS, thr_id);
	  #endif
	  #ifdef USE_AVX512
		twopmodq78_3WORD_DOUBLE_q32(0ull, k_to_try, NTHREADS, thr_id);
		twopmodq78_3WORD_DOUBLE_q64(0ull, k_to_try, NTHREADS, thr_id);
	  #endif
	#endif
	}	// End (inits)

/* Do a quick series of self-tests: */
  #if 1//def FAC_DEBUG
	test_fac();
  #endif

// Oct 2015: GCD-associated self-tests provides a fair bit of added coverage of the mi64 library, so always include:
  #ifdef INCLUDE_PM1
	/* Simple self-tester for GCD routines in gcd_lehmer.c: */
	ASSERT(HERE, test_gcd() == 0, "Factor_init : GCD test failed.\n");
exit(0);
  #endif

	/* Make sure a valid exponent string has been given - if this is the only
	command-line parameter, will attempt to read the other needed run parameters
	from the corresponding checkpoint file:
	*/
	ASSERT(HERE, STRNEQ(pstring,""),"factor.c : pstring empty!");

	/* -bmin/bmax used to set bounds for factoring: */
	if(bmin || bmax) {
		ASSERT(HERE, (kmin==0 && kmax==0 && kplus==0),"(kmin==0 && kmax==0 && kplus==0)");

		if(bmin < 0) {
			fprintf(stderr,"ERROR: log2(min factor) must be >= 0. Offending entry = %lf.\n", bmin);		ASSERT(HERE, 0,"0");
		} else if(bmin >= MAX_BITS_Q) {
			fprintf(stderr,"ERROR: log2(min factor) exceeds allowable limit of %u. Offending entry = %lf.\n", MAX_BITS_Q, bmin);	ASSERT(HERE, 0,"0");
		}

		if(bmax <= 0) {
			fprintf(stderr,"ERROR: log2(max factor) must be > 0. Offending entry = %lf.\n", bmax);		ASSERT(HERE, 0,"0");
		} else if(bmax > MAX_BITS_Q) {
			fprintf(stderr,"ERROR: log2(max factor) exceeds allowable limit of %u. Offending entry = %lf.\n", MAX_BITS_Q, bmax);	ASSERT(HERE, 0,"0");
		}

		if(bmax < bmin) {
			fprintf(stderr,"ERROR: (bmax = %lf) < (bmin = %lf)!\n", bmax, bmin);	ASSERT(HERE, 0,"0");
		}
	}

	/* -kmin/kmax used to set bounds for factoring: */
	if(kmin || kmax) {
		ASSERT(HERE, kmax != 0 ,"factor.c: kmax not set!");
		ASSERT(HERE, (int64)kmax > 0, "kmax must be 63 bits or less!");
		ASSERT(HERE, (bmin==0 && bmax==0 && kplus==0),"(bmin==0 && bmax==0 && kplus==0)");

		if(kmax < kmin) {
			fprintf(stderr,"ERROR: (kmax = %s) < (kmin = %s)!\n", &char_buf0[convert_uint64_base10_char(char_buf0, kmax)], &char_buf1[convert_uint64_base10_char(char_buf1, kmin)]);
			ASSERT(HERE, 0,"0");
		}
	}

	ASSERT(HERE, bmax > 0.0 || kmax != 0 ,"factor.c: One of bmax or kmax must be set!");

	ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			  || (MODULUS_TYPE == MODULUS_TYPE_MERSMERS)
			  || (MODULUS_TYPE ==   MODULUS_TYPE_FERMAT)
				, "Unsupported modulus type!");

	lenQ = ((uint32)MAX_BITS_Q + 63)>>6;	// This sets upper bound on #words needed to store max. factor candidate

	// Convert power-of-2 exponent to unsigned int form and allocate the exponent-storage vector.
	// We use MAX_BITS_P (defined in Mdata.h) to set the allocated storage here, but use the user-set
	// exponent to set the number of words of that allocated storage which are actually used:
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		findex = convert_base10_char_uint64(pstring);
		nbits_in_p = findex+1;
		lenP = (nbits_in_p + 63)>>6;
		p     = (uint64 *)calloc( ((uint32)MAX_BITS_P + 63)>>6, sizeof(uint64));
		p[0] = 1;	mi64_shl(p,p,findex,lenP);	// p = (uint64)1 << findex;
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS)
	{
		findex = convert_base10_char_uint64(pstring);	// This var was really named as abbreviation of "Fermat index", but re-use for MMp
		nbits_in_p = findex;
		if(findex > 1000) {	// Large MMp need deeper sieving on each k passing the default sieve
			kdeep = (uint32 *)calloc( 1024, sizeof(uint32));
			ASSERT(HERE, kdeep != 0x0, "Calloc of kdeep[] failed!");
		}
		lenP = (nbits_in_p + 63)>>6;
		p     = (uint64 *)calloc( ((uint32)MAX_BITS_P + 63)>>6, sizeof(uint64));
		p[0] = 1;	mi64_shl(p,p,findex,lenP);
		mi64_sub_scalar(p,1,p,lenP);	// p = 2^findex - 1;
	#ifdef FAC_DEBUG
		printf("%s(%s) = M(p) with p = %s\n", NUM_PREFIX[MODULUS_TYPE], pstring, &char_buf0[convert_mi64_base10_char(char_buf0, p, lenP, 0)]);
	#endif
	} else {
		// Convert stringified exponent to mi64 form, using same #limbs as for factor candidates:
		p = convert_base10_char_mi64(pstring, &lenQ);	// This does the mem-alloc for us in this case
		lenP = mi64_getlen(p, lenQ); ASSERT(HERE, lenP > 0, "factor.c: Error converting pstring!");
		nbits_in_p = (lenP<<6) - mi64_leadz(p, lenP);
	}

	// Allocate the other modulus-dependent vectors:
	two_p   = (uint64 *)calloc(lenQ, sizeof(uint64));
	p2NC    = (uint64 *)calloc(lenQ, sizeof(uint64));
	q       = (uint64 *)calloc(lenQ * NTHREADS, sizeof(uint64));
	q2      = (uint64 *)calloc(lenQ * NTHREADS, sizeof(uint64));
	u64_arr = (uint64 *)calloc(lenQ * NTHREADS, sizeof(uint64));

	// Now use the just-allocated vector storage to compute how many words are really needed for qmax.
	// Since the sieving always proceeds in full passes through the bit-cleared sieve, the actual kmax used
	// may be up to (len*64)-1 larger than the user-specified kmax:
	if(kmax) {
		interval_hi = (uint64)ceil( (double)kmax / ((uint64)len << TF_CLSHIFT) );	// Copied from restart-file code below
		// Actual kmax used at runtime = interval_hi*(len << TF_CLSHIFT);
		u64_arr[lenP] = mi64_mul_scalar( p, 2*interval_hi*(len << TF_CLSHIFT), u64_arr, lenP);
		lenQ = lenP + (u64_arr[lenP] != 0);
	} else {
		lenQ = ( (uint32)(ceil(bmax)) + 63 ) >> 6;
	}

	// Mersenne numbers must have odd (check primality further on) exponents:
	if((MODULUS_TYPE != MODULUS_TYPE_FERMAT) && (p[0] & 1) == 0)
    {
		fprintf(stderr,"p must be odd! Offending p = %s\n", pstring); ASSERT(HERE, 0,"0");
	}

	/* For purposes of the bits-in-p limit, treat Fermat numbers as having 2^findex rather than 2^findex + 1 bits: */
	if((nbits_in_p - (MODULUS_TYPE == MODULUS_TYPE_FERMAT)) > MAX_BITS_P)
	{
		fprintf(stderr,"p too large - limit is %u bits. Offending p = %s\n", MAX_BITS_P, pstring);
		ASSERT(HERE, 0,"0");
	}
	// To track lg(q) = lg(2.k.p+1), use approximation q ~= 2.k.p, thus lg(q) ~= lg(2.p) + lg(k).
	fbits_in_2p = (double)mi64_extract_lead64(p, lenP, &itmp64) - 64;
//printf("fbits_in_2p = mi64_extract_lead64[= %10u] - 64 = %10.4f\n",mi64_extract_lead64(p, lenP, &itmp64),fbits_in_2p);
	fbits_in_2p += log((double)itmp64)*ILG2 + 1;	// Add 1 to lg(p) to get lg(2p)
//printf("fbits_in_2p += log((double)itmp64)*ILG2 [= %10.4f] = %10.4f\n",log((double)itmp64)*ILG2,fbits_in_2p);
  #if 0	// 11/2013: No clue what I was thinking here...
	// If 2p < 2^64 we left-justify the leading bits to make result lie in [2^63, 2^64), so result here must always be > 2^63:
	ASSERT(HERE, fbits_in_2p >= 63, "fbits_in_2p out of range!");
	fbits_in_2p += nbits_in_p - 64.0;	// lg(2.p) ... Cast 64 to double to avoid signed-int subtract of RHS terms.
  #endif
	// Do some quick sanity tests of exponent for the various kinds of moduli:
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, findex == mi64_trailz(p, lenP), "Internal Fermat-exponent bad power of 2!");
		mi64_shrl(p, q, findex, lenP,lenP);
		mi64_sub_scalar(q, 1ull, q, lenP);
		ASSERT(HERE, mi64_iszero(q, lenP), "Internal Fermat-exponent not a power of 2!");
	}
	else
	{
		// For M(M(p)), make sure the M(p) is actually prime:
		if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS)
		{
			for(i=0; knowns[i] != 0; i++)
			{
				if(findex == knowns[i]) { break; }
			}
			ASSERT(HERE, (knowns[i] != 0), "Double-Mersenne exponent not a known Mersenne prime!");

			// And now proceed to all-binary-ones test of vector-form M(p):
			mi64_add_scalar(p, 1ull, q, lenP);
			ASSERT(HERE, findex == mi64_trailz(q, lenP), "Internal M(M(p))-exponent bad power of 2!");
			mi64_shrl(q, q, findex, lenP,lenP);
			mi64_sub_scalar(q, 1ull, q, lenP);
			ASSERT(HERE, mi64_iszero(q, lenP), "Internal M(M(p))-exponent fails all-binary-ones check!");
		}
		// We can use a lookup table vs known M(p) for all cases, but if Mersenne or M(M(p)) with suitably small p,
		// add a base-2 Fermat PRP test, more as a self-test of the various modpow routines than anything else:
		if(lenP < 1000) {
			mi64_sub_scalar(p, 1ull, q, lenP);	/* q = p-1 */
			if(!mi64_twopmodq(q, lenP, 0, p, lenP, 0x0))
			{
				fprintf(stderr,"WARNING: p = %s is not prime ... proceeding anyway, on presumption user wants this.\n", pstring);
			//	ASSERT(HERE, 0,"0");	Dec 2019 ... allowing odd composite exponents can still be useful, e.g. ATH used to TF M(p^2) for known Mersenne primes
			}
		}
	}

	/* 2*p: Don't need to worry about overflow here since we've allocated two*p, p2NC, q, etc to be of lenQ, not lenP: */
	two_p[lenP] = mi64_add(p, p, two_p, lenP);	// Need to account for fact that 2p may have 1 more word than p
												// (I.e. use lenQ rather than lenP for multiword ops on two_p).
	/* 2*p*[number of classes]: */
	p2NC[lenP] = mi64_mul_scalar(p, (uint64)2*TF_CLASSES, p2NC, lenP);

  #ifdef FAC_DEBUG
	u64_arr[lenP] = mi64_mul_scalar(two_p,k_targ,u64_arr,lenP);	u64_arr[0] += 1;	// q = 2.k.p + 1
	printf("FAC_DEBUG: Doing targeted debug-TF of %s(%s) with target factor candidate q = %s\n",NUM_PREFIX[MODULUS_TYPE],pstring,&char_buf0[convert_mi64_base10_char(char_buf0, u64_arr, lenQ, 0)]);
	printf("two_p        = %s\n", &char_buf0[convert_mi64_base10_char(char_buf0, two_p, lenQ, 0)]);
	printf("TF_CLASSES = %u\n",(uint32)TF_CLASSES);
  #endif

	// p mod TF_CLASSES:
	if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS)
	{
		pmodNC = twopmmodq64(findex, (uint64)TF_CLASSES) - 1;	// For double-Mersenne factoring need M(p) mod #TF_CLASSES
		/*
		The above routine computes 2^p (mod 60) via Montgomery-mul-based powering. That requires an odd modulus,
		so the actual powering step computes 2^(p-2) (mod 15), multiplies the result by 4 to get 2^p (mod 60),
		and subtracts 1 to get M(p) (mod 60).

		Note that there is a shortcut to obtaining 2^(p-2) (mod 15), namely summing the hex digits of 2^(p-2)
		(mod 15). Since the hexadecimal base 16 == 1 (mod 15), this sum gives the desired result (mod 15).
		Since 2^(p-2) i binary is just a 1 followed by (p-2) binary zeros, in base-16 it is just a leading hex
		digit d = 2^((p-2)%4) followed by a string of (p-2)/4 hexadecimal zeros, the latter of which contribute
		0 to the hex-digit sum (mod 15). Thus

			2^(p-2) == 2^((p-2)%4) (mod 15), whence

			M(p) == 4*2^((p-2)%4) - 1 (mod 60).

		For example, for p = 521 we have 2^(p-2) == 2^3 == 8 (mod 15), whence M(p) == 31 (mod 60).

		For p = 607 we 2^(p-2) == 2^1 == 2 (mod 15), whence M(p) == 7 (mod 60).

		In fact for any odd-exponent M(p) there are only the two possibilities p == 1 or 3 (mod 4), for
		which 2^(p-2) == 2^3 or 2^1 (mod 15) and M(p) == 31 or 7 (mod 60), respectively.

		For these two M(p) (mod 60) values the respective sets of 16 eligible k (mod 60) values of possible factors
		q = 2.k.M(p) + 1 are

			p == 1 (mod 4), M(p) == 31 (mod 60) : k (mod 60) = (any of) 0, 5, 8, 9,20,21,24,29,33,36,41,44,45,48,53,56

			p == 3 (mod 4), M(p) ==  7 (mod 60) : k (mod 60) = (any of) 0, 5, 8, 9,12,17,20,24,29,32,33,44,45,48,53,57 .

		M(p) == 31 or 7 (mod 60) also implies that the following check is not needed, but include it for formal completeness,
		and in case someone else modifies this code for a purpose where the above call might in fact return 0 (mod 60):
		*/
		printf("p mod %u = %d\n", TF_CLASSES, pmodNC);
		itmp32 = mi64_div_y32(p, TF_CLASSES, 0x0, lenP);
		if(pmodNC != itmp32) {
			printf("p mod %u v2 = %d\n", TF_CLASSES, itmp32);
			fprintf(stderr,"Warning: Differing (p %% TF_CLASSES) values from Powering and direct-long-div! Proceeding using the 2nd result (%u).\n",itmp32);
			pmodNC = itmp32;
		}
	} else {
		pmodNC = mi64_div_y32(p, TF_CLASSES, 0x0, lenP);
	}

	// If user-set kmax, test factoring range vs internal limits
	if(kmax) {
		interval_hi = (uint64)ceil((double)kmax/((uint64)len << TF_CLSHIFT));	// Copied from restart-file code below
		u64_arr[lenP] = mi64_mul_scalar( p, 2*interval_hi*(len << TF_CLSHIFT), u64_arr, lenP);
		ASSERT(HERE, lenQ == lenP+(u64_arr[lenP] != 0), "");

		nbits_in_q = (lenQ<<6) - mi64_leadz(u64_arr, lenQ);

		if(nbits_in_q > MAX_BITS_Q)
		{
			fprintf(stderr,"qmax too large - limit is %u bits. Offending p, kmax = %s, %s\n", MAX_BITS_Q, pstring, &char_buf0[convert_uint64_base10_char(char_buf0, kmax)]);
			ASSERT(HERE, 0,"0");
		}
	}

	/* log2[nearest power of 2 to (nbits_in_p)*lenQ^2)] */
	bits_in_pq2 = nbits_in_p*lenQ*lenQ;
	bits_in_pq2 = 32 - leadz32(bits_in_pq2);
	CMASKBITS = (30 - (bits_in_pq2>>1));
	countmask = (1ull << CMASKBITS) - 1;

/*****************************************************/
/****************** RESTART STUFF: *******************/
/*****************************************************/

	/* Restart file for a given exponent is named 't{exponent}'.
	Since Fermat-number exponents are so much smaller than Mersenne-number ones,
	we assume there is no overlap, i.e. if pstring <= MAX_BITS_P, it's a
	Fermat-number factoring run, pstring > MAX_BITS_P is a Mersenne-number run.
	*/
	RESTARTFILE[0] = 't'; RESTARTFILE[1] = '\0'; strcat(RESTARTFILE, pstring);
	// Checkpointing only supported for single-threaded runs:
	if(NTHREADS > 1)
		fprintf(stderr,"WARN: Checkpointing only supported for single-threaded runs!\n");
	else
		fprintf(stderr,"INFO: Will write checkpoint data to savefile %s.\n",RESTARTFILE);

	fprintf(stderr,"INFO: Will write savefile %s every 2^%u = %llu factor candidates tried.\n",RESTARTFILE,CMASKBITS,countmask+1);

	/**** process restart-file and any command-line params: ****/
	// Note: return value of read_savefile is signed:
	itmp = read_savefile(RESTARTFILE, pstring, &bmin_file,&bmax_file, &kmin_file,&know_file,&kmax_file, &passmin_file,&passnow_file,&passmax_file, &count);
	if(itmp == -1) {
		sprintf(cbuf,"INFO: No factoring savefile %s found ... starting from scratch.\n",RESTARTFILE);
		fprintf(stderr,"%s",cbuf);
	#ifndef FACTOR_STANDALONE
		fq = mlucas_fopen(STATFILE,"a"); fprintf(fq,"%s",cbuf); fclose(fq); fq = 0x0;
	#endif
		// Init savefile with above read_savefile fields so ensuing checkpoint-writes only need to update the pass# and k:
//		ASSERT(HERE,0 == init_savefile(RESTARTFILE, pstring, bmin,bmax, kmin,know,kmax, passmin,passnow,passmax, count),"init_savefile failed!");
	} else {
		ASSERT(HERE,!itmp,"There were errors reading the savefile ... aborting");
		count = 0ull;	// Need to reset == 0 prior to sieving so kvector-fill code works properly

		/* If previous run is not yet complete, ignore any increased factor-bound-related
		command-line parameters and instead proceed to complete the previous run first:
		*/
		if((know_file < kmax_file) || (passnow_file < passmax_file)) {
			incomplete_run = TRUE;
			fprintf(stderr,"INFO: Previous run to kmax = %s not yet complete.\n"  , &char_buf0[convert_uint64_base10_char(char_buf0, kmax_file)]);
			fprintf(stderr,"Ignoring any increased factor-bound-related command-line parameters and proceeding to complete previous run.\n");
			bmin = bmin_file; bmax = bmax_file;
			passmin = passmin_file; passnow = passnow_file; passmax = passmax_file;
			kmin = kmin_file; know = know_file; kmax = kmax_file;
			kplus = 0;
		} else {
			/**** Previous run was completed - check that current params satisfy one (and only one)
			of the following sets of conditions:
				1) -bmin/bmax used to set bounds for factoring:
					In this case we expect any command-line bmin will be >= that in the restart file
					(in fact we expect bmin >= bmax_file, i.e. that the runs are nonoverlapping -
					if not we warn and set bmin = bmax_file), and that bmax > bmax_file.
			****/
			if(bmin || bmax) {
			#if(!defined(P1WORD))
			//	ASSERT(HERE, 0,"bmin/bmax form of bounds-setting only allowed for single-word-p case!");
			#endif
				ASSERT(HERE, (kmin==0 && kmax==0 && kplus==0),"(kmin==0 && kmax==0 && kplus==0) - please delete any restart files for this p and retry debug run.");

				if(bmin) {
					ASSERT(HERE, bmin >= bmin_file - 0.0000000001,"bmin >= bmin_file");
					if(bmin < bmax_file)
						fprintf(stderr,"WARNING: Specified bmin (%lf) smaller than previous-run bmax = %lf. Setting equal to avoid overlapping runs.\n", bmin, bmax_file);
				}
				bmin = bmax_file;
				/* We expect any command-line bmax will be > that in the restart file: */
				if(bmax)
					ASSERT(HERE, bmax > bmax_file - 0.0000000001,"bmax >= bmax_file");
			}

			/****
				2) -kmin/kmax used to set bounds for factoring:
					In this case we expect any command-line kmin will be >= that in the restart file
					(in fact we expect kmin >= kmax_file, i.e. that the runs are nonoverlapping -
					if not we warn and set kmin = kmax_file), and that kmax > kmax_file.
			****/
			if(kmin || kmax) {
				ASSERT(HERE, (bmin==0 && bmax==0 && kplus==0),"(bmin==0 && bmax==0 && kplus==0)");
				if(kmin) {
					ASSERT(HERE, kmin >= kmin_file,"kmin >= kmin_file");
					if(kmin < kmax_file)
						fprintf(stderr,"WARNING: Specified kmin (%s) smaller than previous-run kmax = %s. Setting equal to avoid overlapping runs.\n", &char_buf0[convert_uint64_base10_char(char_buf0, kmax)], &char_buf1[convert_uint64_base10_char(char_buf1, kmax_file)]);
				}
				kmin = kmax_file;
				/* We expect any command-line kmax will be > that in the restart file: */
				if(kmax)
					ASSERT(HERE, kmax > kmax_file,"kmax >= kmax_file");
			}

			/****
				3) -kplus used to increment an upper bound from a previous factoring run:
			****/
			if(kplus) {
				ASSERT(HERE, (bmin==0 && bmax==0 && kmin==0 && kmax==0),"(bmin==0 && bmax==0 && kmin==0 && kmax==0)");
				kmin = kmax_file;
				/* Ensure incremented value kmax fits into a 64-bit unsigned int: */
				ASSERT(HERE, (kmin + kplus) > kplus, "kmax_file + kplus exceeds 2^64!");
				kmax = kmin + kplus;
				kplus = 0;	/* If kplus != 0 detected further on, that indicates that no valid restart
							file was found for factoring-bounds incrementing. */
			}
		}
		/* Successfully processed restart file: */
		restart = TRUE;
	}

/************************ END(RESTART STUFF) *******************/

  #warning bmax/kmax-synchro needs re-do!
	/* If it's not a restart of an as-yet-uncompleted run, synchronize the factoring-bound params: */
	if(!incomplete_run)
	{
		/* Double-check factoring pass bounds: */
		if(passmin > (TF_PASSES-1) )
		{
			fprintf(stderr,"ERROR: passmin must be <= %u. Offending entry = %u.\n", TF_PASSES-1, passmin);
			ASSERT(HERE, 0,"0");
		}

		if(passmax < passmin)
		{
			fprintf(stderr,"ERROR: (passmax = %u) < (passmin = %u)!\n", passmax, passmin);
			ASSERT(HERE, 0,"0");
		}
		if(passmax > (TF_PASSES-1) )
		{
			fprintf(stderr,"ERROR: passmax must be <= %u. Offending entry = %u.\n", TF_PASSES-1, passmax);
			ASSERT(HERE, 0,"0");
		}

		/**** Process factor candidate bounds: ****/

		/* If any of bmin|kmin, bmax|kmax nonzero, calculate its counterpart: */
	#ifdef P1WORD
		/* Find FP approximation to 2*p - can't use this for multiword case, because double approximation tp 2*p may overflow: */
		twop_float = (double)two_p[0];
	#endif
		/* Compute kmax if not already set: */
		if(!kmax) {
			ASSERT(HERE, bmax <= (nbits_in_p+65), "Specified bmax implies kmax > 64-bit, which exceeds the program's limit ... aborting.");
			kmax = given_b_get_k(bmax, two_p, lenQ);
			ASSERT(HERE, kmax > 0, "Something went wrong with the computation of kmax ... possibly your bmax implies kmax > 64-bit?");
		}
		if(kmin || bmin) {
			if(kmin == 0ull) {	/* Lower Bound given in log2rithmic form */
				ASSERT(HERE, bmin <= bmax, "bmin >= bmax!");
				kmin = given_b_get_k(bmin, two_p, lenQ);
			} else {
				ASSERT(HERE, kmin <= kmax, "kmin >= kmax!");
			#ifdef P1WORD
				fqlo = kmin*twop_float + 1.0;
				bmin = log(fqlo)*ILG2;
			#endif
			}
		} else {
		#ifdef P1WORD
			fqlo = 1.0;
		#endif
		}
ASSERT(HERE,0 == init_savefile(RESTARTFILE, pstring, bmin,bmax, kmin,know,kmax, passmin,passnow,passmax, count),"init_savefile failed!");
//**** Do savefile-init here? ******
		if(kmax || bmax) {
			if(kmax == 0ull) {	/* Upper Bound given in log2rithmic form */
				kmax = given_b_get_k(bmax, two_p, lenQ);
 			} else {
			#ifdef P1WORD
				fqhi = kmax*twop_float + 1.0;
				bmax = log(fqhi)*ILG2;
			#endif
			}
		} else
			ASSERT(HERE, 0 ,"factor.c : One of bmax, kmax must be nonzero!");

		/**** At this point the paired elements bmin|kmin, bmax|kmax are in synchrony. ****/

		/* If kplus given on command line, a valid restart file should have been found
		and kmax incremented at this point, i.e. kplus should have been reset to zero:
		*/
		ASSERT(HERE, kplus == 0, "kplus must be zero here!");

		know = kmin;
		passnow = passmin;

	}	/* endif(!incomplete_run) */

/*****************************************************/
/****************** SIEVE STUFF: *********************/
/*****************************************************/

	ASSERT(HERE, NUM_SIEVING_PRIME > 0, "factor.c : NUM_SIEVING_PRIME > 0");

/*   allocate the arrays and initialize the array of sieving primes	*/
	temp_late = (uint64 *)calloc(len, sizeof(uint64));

  #if TF_CLASSES == 60
	i = len/TF_CLASSES + 1;	// len not divisible by TF_CLASSES, so add a pad-word
  #else
	i = (len*64)/TF_CLASSES + 1;	// 64*len divisible by TF_CLASSES, no need for padding
		//**** Oct 2016: AVX-512 vector-bit-clear needs a padding element, so add one. ****
  #endif
	bit_map = (uint64 *)calloc(i * NTHREADS, sizeof(uint64));
	bit_map2= (uint64 *)calloc(i * NTHREADS, sizeof(uint64));	// 2nd alloc to give each thread 1 bit-clearable copy of master bit_map
	if (bit_map == NULL) {
		fprintf(stderr,"Memory allocation failure for BITMAP array");
		ASSERT(HERE, 0,"0");
	}
	bit_atlas = (uint64 *)calloc(i * TF_PASSES, sizeof(uint64));
	if (bit_atlas == NULL) {
		fprintf(stderr,"Memory allocation failure for TEMPLATE array");
		ASSERT(HERE, 0,"0");
	}
printf("Allocated %u words in master template, %u in per-pass bit_map [%u x that in bit_atlas]\n",len,i,TF_PASSES);

  #ifdef USE_AVX512	// Use vector-int math and gather-load/scatter-store to accelerate the bit-clearing
	psmall = (uint32 *)calloc(NUM_SIEVING_PRIME * NTHREADS, sizeof(uint32));
	if (psmall == NULL) {
		fprintf(stderr,"Memory allocation failure for PSMALL array");
		ASSERT(HERE, 0,"0");
	}
  #endif

	pdiff = (uint8 *)calloc(NUM_SIEVING_PRIME * NTHREADS, sizeof(uint8));
	if (pdiff == NULL) {
		fprintf(stderr,"Memory allocation failure for pdiff array");
		ASSERT(HERE, 0,"0");
	}

	startval = (uint32 *)calloc(NUM_SIEVING_PRIME * NTHREADS, sizeof(uint32));
	if (startval == NULL) {
		fprintf(stderr,"Memory allocation failure for STARTVAL array");
		ASSERT(HERE, 0,"0");
	}

	pinv = (uint32 *)calloc(NUM_SIEVING_PRIME, sizeof(uint32));
	if (pinv == NULL) {
		fprintf(stderr,"Memory allocation failure for PINV array");
		ASSERT(HERE, 0,"0");
	}

  #if DBG_SIEVE
	startval_incr = (uint32 *)calloc(NUM_SIEVING_PRIME, sizeof(uint32));
	if (startval_incr == NULL) {
		fprintf(stderr,"Memory allocation failure for STARTVAL_INCR array");
		ASSERT(HERE, 0,"0");
	}
  #endif

		/* Check integrity (at least in the sense of monotonicity) for the precomputed pseudoprime table: */
		for(i = 1; i < 9366; ++i) {
			ASSERT(HERE, f2psp[i] > f2psp[i-1],"Misplaced pseudoprime!");
		}

		/* Test some near-2^32 known-prime cases: */
		curr_p = (uint32)-5;
		itmp32 = twopmodq32(curr_p-1, curr_p);
		ASSERT(HERE, itmp32 == 1,"twopmodq32: 2^32 - 5 test fails!");
		curr_p = (uint32)-17;
		itmp32 = twopmodq32(curr_p-1, curr_p);
		ASSERT(HERE, itmp32 == 1,"twopmodq32: 2^32 -17 test fails!");
		curr_p = (uint32)-35;	/* Start of the last length-30 curr_p%30 == 11 interval < 2^32; the 6th candidate in that interval, 2^32-17, is prime */
		itmp32 = twopmodq32_x8(curr_p, curr_p+ 2, curr_p+ 6, curr_p+ 8, curr_p+12, curr_p+18, curr_p+20, curr_p+26);
		ASSERT(HERE, itmp32 ==32,"twopmodq32_x8: 2^32 -35 test fails!");

		fprintf(stderr,"Generating difference table of first %u small primes\n", nprime);
		curr_p = 3;	/* Current prime stored in l. */
		max_diff = 0;

		f2psp_idx = 0;	/* Index to next-expected Fermat base-2 pseudoprime in the precomputed table */

		/* Init first few diffs between 3/5, 5/7, 7/11, so can start loop with curr_p = 11 == 1 (mod 10), as required by twopmodq32_x8(): */
		pdiff[0] = 0;	pdiff[1] = pdiff[2] = 1;
		ihi = curr_p = 11;
	#ifdef USE_AVX512	// Use vector-int math and gather-load/scatter-store to accelerate the bit-clearing
		psmall[0] = 3; psmall[1] = 5; psmall[2] = 7;
	#endif
		/* Process chunks of length 30, starting with curr_p == 11 (mod 30). Applying the obvious divide-by-3,5 mini-sieve,
		we have 8 candidates in each interval: curr_p + [ 0, 2, 6, 8,12,18,20,26].
		For example: curr_p = 11 gives the 8 candidates: 11,13,17,19,23,29,31,37.
		*/
		for(i = 3; i < nprime; curr_p += 30)
		{
			/* Make sure (curr_p + 29) < 2^32: */
			if(curr_p > 0xffffffe3) {
				fprintf(stderr,"curr_p overflows 32 bits!");
				nprime = i;
				break;
			}

			/* Max sieving prime must be < smallest candidate factor of M(p) */
		#ifdef P1WORD
			if((curr_p+29) > two_p[0]) {
				nprime = i;
				break;
			}
		#endif

			/* Do a quick Fermat base-2 compositeness test before invoking the more expensive mod operations: */
			itmp32 = twopmodq32_x8(curr_p, curr_p+ 2, curr_p+ 6, curr_p+ 8, curr_p+12, curr_p+18, curr_p+20, curr_p+26);
			for(j = 0; j < 8; ++j)
			{
				if((itmp32 >> j)&0x1)	// It's a PRP, so check against the table of known pseudoprimes and
				{						// (if it's not a PSP) init for the next gap
					ASSERT(HERE, curr_p <= f2psp[f2psp_idx],"Error in pseudoprime sieve");
					if((curr_p + pdsum_8[j]) == f2psp[f2psp_idx])	/* It's a base-2 pseudoprime */
					{
						++f2psp_idx;
						pdiff[i] += pdiff_8[j];
						continue;
					}
					else	/* It's prime - add final increment to current pdiff[i] and then increment i: */
					{
						ihi = (curr_p + pdsum_8[j]);
					#ifdef USE_AVX512	// Use vector-int math and gather-load/scatter-store to accelerate the bit-clearing
						psmall[i] = ihi;
					#endif
						pdiff[i] += pdiff_8[j];
						if(pdiff[i] > max_diff)
						{
							max_diff = pdiff[i];
						#if DBG_SIEVE
							printf("pdiff = %d at curr_p = %u\n", 2*max_diff,ihi);
						#endif
						}
						if(++i == nprime)
						{
							break;
						}
					}
				}
				else
				{
					pdiff[i] += pdiff_8[j];
				}
			}
		}
		MAX_SIEVING_PRIME = ihi;
	#ifdef MULTITHREAD
		uint8 *byteptr = pdiff;	// Each thread gets its own copy of the pdiff data:
		for(thr_id = 1; thr_id < NTHREADS; thr_id++) {
			byteptr += NUM_SIEVING_PRIME;
			memcpy(byteptr, pdiff, NUM_SIEVING_PRIME);
		}
	#endif

	#if 1//def FAC_DEBUG
		printf("Using first %u odd primes; max gap = %u\n",nprime,2*max_diff);
		printf("max sieving prime = %u\n",MAX_SIEVING_PRIME);
	#endif

  #if 0
	// Oct 2015: Play with Smarandache numbers ():
	i = 2000000;	ASSERT(HERE, i <= nprime, "prime limit exceeded in testSmarandache!");
	testSmarandache(100001,101000, pdiff, i);
	exit(0);
  #endif
  #if 0
	// Oct 2018: Play with "sieve survivors" stats: lim(n --> oo) prod_(p <= n)(1-1/p)/(1/ln(p^2))
	i = 1000000000;	ASSERT(HERE, i <= MAX_SIEVING_PRIME, "prime limit exceeded in testSieveProdAsymp!");
	struct qfloat qfprod = QHALF, qt;
	double prod = 0.5, log_psq = log((double)i*i);
	for(m = 0, curr_p = 3; m < nprime; m++) {
		curr_p += (pdiff[m] << 1);
		if(curr_p > i) break;
		prod *= (1-1./curr_p);
//	printf("p = %u: prod = %18.15f\n",curr_p,prod);
		qt = qfsub(QONE,qf_rational_quotient(1ull,(uint64)curr_p));
		qfprod = qfmul(qfprod,qt);
	}
	printf("Used primes <= %u: 1/ln(p^2) = %18.15f, prod_(p <= n)(1-1/p) = %18.15f, ratio = %18.15f, qfprod = %18.15f\n",i,log_psq,prod,log_psq*qfdbl(qfprod),qfdbl(qfprod));
	exit(0);
  #endif
/* Time the vector trialdiv stuff: */
  #if TEST_TRIALDIV
	for(i = 0; i < vec_len; i++)
	{
		xvec[i]  = rng_isaac_rand();
	}

	clock1 = clock();
	curr_p = 3;
	for(m = 0; m < nprime; m++)
	{
		curr_p += (pdiff[m] << 1);
		if(mi64_is_div_by_scalar32(xvec, curr_p, vec_len) == TRUE)
			printf("mi64_is_div_by_scalar32 test: %10u is a divisor\n", curr_p);
	}

	clock2 = clock();	/* Assume sieve setup time < 2^32 cycles - even if that is violated it's no big deal at this point. */
	tdiff = (double)(clock2 - clock1);
	citer = tdiff*2000000000.0/CLOCKS_PER_SEC;
	citer /= (double)vec_len*nprime;
	printf	("Elapsed Time =%s; cycles/iter = %10.2f\n",get_time_str(tdiff),citer);

	clock1 = clock();
	curr_p = 3;
	for(m = 0; m < nprime; m+=4)
	{
		curr_p += (pdiff[m  ] << 1); tryq[0] = curr_p;
		curr_p += (pdiff[m+1] << 1); tryq[1] = curr_p;
		curr_p += (pdiff[m+2] << 1); tryq[2] = curr_p;
		curr_p += (pdiff[m+3] << 1); tryq[3] = curr_p;
		j = mi64_is_div_by_scalar32_x4(xvec, tryq[0], tryq[1], tryq[2], tryq[3], vec_len);
		if(j != 0)
		{
			for(i = 0; i < 4; ++i)
			{
				if((j >> i)&1)
					printf("mi64_is_div_by_scalar32_x4 test: %10u is a divisor\n", tryq[i]);
			}
		}
	}

	clock2 = clock();	/* Assume sieve setup time < 2^32 cycles - even if that is violated it's no big deal at this point. */
	tdiff = (double)(clock2 - clock1);
	citer = tdiff*2000000000.0/CLOCKS_PER_SEC;
	citer /= (double)vec_len*nprime;
	printf	("Elapsed Time =%s; cycles/iter = %10.2f\n",get_time_str(tdiff),citer);

	clock1 = clock();
	curr_p = 3;
	for(m = 0; m < nprime; m += 8)
	{
		curr_p += (pdiff[m  ] << 1); tryq[0] = curr_p;
		curr_p += (pdiff[m+1] << 1); tryq[1] = curr_p;
		curr_p += (pdiff[m+2] << 1); tryq[2] = curr_p;
		curr_p += (pdiff[m+3] << 1); tryq[3] = curr_p;
		curr_p += (pdiff[m+4] << 1); tryq[4] = curr_p;
		curr_p += (pdiff[m+5] << 1); tryq[5] = curr_p;
		curr_p += (pdiff[m+6] << 1); tryq[6] = curr_p;
		curr_p += (pdiff[m+7] << 1); tryq[7] = curr_p;
		j = mi64_is_div_by_scalar32_x8(xvec, tryq[0], tryq[1], tryq[2], tryq[3], tryq[4], tryq[5], tryq[6], tryq[7], vec_len);
		if(j != 0)
		{
			for(i = 0; i < 8; ++i)
			{
				if((j >> i)&1)
					printf("mi64_is_div_by_scalar32_x8 test: %10u is a divisor\n", tryq[i]);
			}
		}
	}

	clock2 = clock();	/* Assume sieve setup time < 2^32 cycles - even if that is violated it's no big deal at this point. */
	tdiff = (double)(clock2 - clock1);
	citer = tdiff*2000000000.0/CLOCKS_PER_SEC;
	citer /= (double)vec_len*nprime;
	printf	("Elapsed Time =%s; cycles/iter = %10.2f\n",get_time_str(tdiff),citer);

	clock1 = clock();
	curr_p = 3;
	for(m = 0; m < nprime; m++)
	{
		curr_p += (pdiff[m] << 1);
		if(mi64_is_div_by_scalar64(xvec, (uint64)curr_p, vec_len) == TRUE)
			printf("mi64_is_div_by_scalar64 test: %10u is a divisor\n", curr_p);
	}
	clock2 = clock();	/* Assume sieve setup time < 2^32 cycles - even if that is violated it's no big deal at this point. */
	tdiff = (double)(clock2 - clock1);
	citer = tdiff*2000000000.0/CLOCKS_PER_SEC;
	citer /= (double)vec_len*nprime;
	printf	("Elapsed Time =%s; cycles/iter = %10.2f\n",get_time_str(tdiff),citer);

	free((void*)xvec);
  #endif

/*   for p < max prime in precomputed table, need to truncate the range of primes...	*/
	/*********** WHY WAS THIS HERE? *************
	curr_p = MAX_SIEVING_PRIME;
	for(;;)
	{
		if(p > curr_p) break;
		curr_p -= (pdiff[nprime--] << 1);
	#ifdef FAC_DEBUG
		ASSERT(HERE, curr_p == prime[nprime], "factor.c : curr_p == prime[nprime]");
	#endif
	}
	MAX_SIEVING_PRIME = curr_p;
	***********/

	/****************** KNOWN-TEST-FACTOR STUFF: *******************/
  #ifdef FAC_DEBUG

	if(k_targ) {
		// Could add check of whether associated q is prime, but assume user knows what he's soing,
		// and furthermore may be usefull to allow for composite products-of-smaller-prime-factors:
		kmodNC = k_targ%TF_CLASSES;

		printf("p mod %u = %d\n",TF_CLASSES, pmodNC);
		printf("k mod %u = %d\n",TF_CLASSES, kmodNC);

		/* ...and get the pass number on which the factor should be found.
		(Remember that the pass number returned by CHECK_PKMOD[60|4620] is unit-offset).
		If a known factor given, only process the given k/log2 range for that pass:
		*/
	#if TF_CLASSES == 60
		pass_targ = CHECK_PKMOD60  (p,lenP, k_targ, 0x0) - 1;
	#else	// 4620 classes:
		pass_targ = CHECK_PKMOD4620(p,lenP, k_targ, 0x0) - 1;
	#endif
		ASSERT(HERE, (pass_targ < TF_PASSES), "Candidate factor set via k_targ is not a possible factor for this exponent!");
		printf("Target pass for debug-factor = %u\n",pass_targ);
	}

  #endif

	itmp64 = (uint64)mi64_div_y32(p,TF_CLASSES,0x0,lenP);
//	printf("p %% 60 = %llu\n",itmp64);

  #if TF_CLASSES == 60
/*
	const int pmod_vec[] = { 1, 7,11,13,17,19,23,29,31,37,41,43,47,49,53,59, 2,4,8,16,32, 0x0};
	for(i = 0; pmod_vec[i] != 0; i++) {
		ASSERT(HERE, CHECK_PKMOD60(pmod_vec[i], k, incr) == 16, "CHECK_PKMOD60 returns something other than the expected #TF_PASSES = 16!\n");
	}
	exit(0);
Mersenne Mp: Acceptable km-values for the 16 possible pm (= p%60) values:
	pm =  1: 0, 3, 8,11,15,20,23,24,35,36,39,44,48,51,56,59
	pm =  7: 0, 5, 8, 9,12,17,20,24,29,32,33,44,45,48,53,57
	pm = 11: 0, 1, 4, 9,13,16,21,24,25,28,33,36,40,45,48,49
	pm = 13: 0, 3, 8,11,12,15,20,23,27,32,35,36,47,48,51,56
	pm = 17: 0, 3, 4, 7,12,15,19,24,27,28,39,40,43,48,52,55
	pm = 19: 0, 5, 9,12,17,20,21,24,29,32,36,41,44,45,56,57
	pm = 23: 0, 1,12,13,16,21,25,28,33,36,37,40,45,48,52,57
	pm = 29: 0, 4, 7,12,15,16,19,24,27,31,36,39,40,51,52,55
	pm = 31: 0, 5, 8, 9,20,21,24,29,33,36,41,44,45,48,53,56
	pm = 37: 0, 3, 8,12,15,20,23,24,27,32,35,39,44,47,48,59
	pm = 41: 0, 3, 4,15,16,19,24,28,31,36,39,40,43,48,51,55
	pm = 43: 0, 5, 8,12,17,20,21,32,33,36,41,45,48,53,56,57
	pm = 47: 0, 4, 9,12,13,24,25,28,33,37,40,45,48,49,52,57
	pm = 49: 0,11,12,15,20,24,27,32,35,36,39,44,47,51,56,59
	pm = 53: 0, 3, 7,12,15,16,27,28,31,36,40,43,48,51,52,55
	pm = 59: 0, 1, 4, 9,12,16,21,24,25,36,37,40,45,49,52,57
Fermat Fn (n > 0): 0,Acceptable km-values for the ? possible pm (= p%60) values:
	pm =  2: 0, 4,10,12,18,22,24,28,30,34,40,42,48,52,54,58
	pm =  4: 0, 2, 6,12,14,20,24,26,30,32,36,42,44,50,54,56
	pm =  8: 0, 6,10,12,16,18,22,28,30,36,40,42,46,48,52,58
	pm = 16: 0, 6, 8,14,18,20,24,26,30,36,38,44,48,50,54,56	<*** F36 factor has k = 20 ... why do I miss? ***
	pm = 32: 0, 4,10,12,18,22,24,28,30,34,40,42,48,52,54,58
*/
	i = CHECK_PKMOD60  (&itmp64,1, k, incr);
	ASSERT(HERE, i == TF_PASSES, "CHECK_PKMOD60 returns something other than the expected #TF_PASSES! Exponent not of the required form (odd prime or odd composite == any_of[1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59] (mod 60).\n");
/*
	printf("k mod 60 = [");
	for(i = 0, j = 0; i < 16; i++) {
		j += incr[i];
		printf("%3u",(uint32)j);
	}
	printf("]\n");
	exit(0);
*/
  #else	// 4620 classes:
	i = CHECK_PKMOD4620(&itmp64,1, k, incr);
	ASSERT(HERE, i == TF_PASSES, "CHECK_PKMOD4620 returns something other than the expected #TF_PASSES! Exponent not of the required form (odd prime or odd composite == any_of[960 possible values] (mod 4620).\n");
  #endif

	/* If it's a restart, interval_lo for the initial pass will be based
	on (know), rather than (kmin) - handle that just inside the pass-loop: */
	/* Using TF_CLASSES = 60 by way of example:
	sievelets have bit_len = (len<<6)/60 ... each bit represents a k-incr of +60, so each pass (a.k.a. 'interval')
	rep. a k-incr of len<<6. Thus to get interval corr. to given kmax, use kmax = interval*(len<<6); yielding
		interval = kmax/(len<< 6) = (kmax>> 6)/len .
	For TF_CLASSES = 4620, we have bit_len = (len<<12)/4620, thus kmax = interval*(len<<12), thus
		interval = kmax/(len<<12) = (kmax>>12)/len .
	*/
	interval_lo = kmin/((uint64)len << TF_CLSHIFT);
	interval_now= know/((uint64)len << TF_CLSHIFT);
	interval_hi = (uint64)ceil((double)kmax/((uint64)len << TF_CLSHIFT));

	/* Make sure we always do at least one full pass through the sieve
	(e.g. if kmin = kmax = 1, ceil(kmax << TF_CLSHIFT) gives 0, same as (kmin << TF_CLSHIFT): */
	if(interval_hi == interval_lo)
		interval_hi += 1;

	ninterval = interval_hi - interval_lo;

	/* The actual kmin/kmax values used in the run are exact multiples
	of len*64 - we always do at least one full pass through the sieve,
	i.e. we run through at least (len*64) worth of k's */
	kmin = interval_lo *(len << TF_CLSHIFT);
	know = interval_now*(len << TF_CLSHIFT);
	kmax = interval_hi *(len << TF_CLSHIFT);

	/* And now that we have the actual kmin/kmax, recalculate these: */
  #ifdef P1WORD
	fqlo = kmin*twop_float + 1.0;
	fqhi = kmax*twop_float + 1.0;
  #endif

	/* 11/14/05: Since we don't actually use bmin/bmax for anything other
	than setting sieving bounds (which then get modified via the above
	k-is-exact-multiple-of-sieve-length anyway), preserve any user-set
	values, since these are typically whole numbers, and look nicer
	in diagnostic and savefile printing:
	*/
	/*	bmin = log(fqlo)/log(2.0);*/
	/*	bmax = log(fqhi)/log(2.0);*/

  #ifdef FAC_DEBUG
	/* Make sure the range of k's for the run contains any target factor: */
	if(k_targ)
		ASSERT(HERE, (kmin <= k_targ) && (kmax >= k_targ),"k_targ not in [kmin, kmax]");
  #endif

  #ifdef FACTOR_STANDALONE
	fp = stderr;
  #else
	fp = mlucas_fopen(STATFILE,"a");
  #endif
	fq = mlucas_fopen(OFILE,"a");

  #ifdef P1WORD
	sprintf(char_buf0, "Searching in the interval k=[%s, %s], i.e. q=[%e, %e]\n", &char_buf1[convert_uint64_base10_char(char_buf1, kmin )], &char_buf2[convert_uint64_base10_char(char_buf2, kmax )],fqlo,fqhi);
  #else
	sprintf(char_buf0, "Searching in the interval k=[%s, %s]\n", &char_buf1[convert_uint64_base10_char(char_buf1, kmin )], &char_buf2[convert_uint64_base10_char(char_buf2, kmax )]);
  #endif
	fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);

	sprintf(char_buf0, "Each of %3u (p mod %u) passes will consist of %s intervals of length %u\n", passmax-passmin+1, TF_CLASSES, &char_buf1[convert_uint64_base10_char(char_buf1, ninterval)], bit_len);
	fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);

	if(passnow != passmin || know != kmin)
	{
		sprintf(char_buf0, "Resuming execution with pass %u and k = %s\n", passnow, &char_buf1[convert_uint64_base10_char(char_buf1, know )]);
		fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);
		sprintf(char_buf0, "#Q tried = %s\n", &char_buf1[convert_uint64_base10_char (char_buf1, count)] );
		fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);
	count = 0;	// Reset == 0 prior to sieving so kvector-fill code works properly
	}

  #ifndef FACTOR_STANDALONE
	fclose(fp);
  #endif
	fp = 0x0;
	fclose(fq); fq = 0x0;

	//...Init clock counter:
  #ifdef CTIME
	clock1 = clock();
  #else
	clock1 = getRealTime();
  #endif
	tdiff = 0.0;

	// quick way to set all the bits = 1:
	for(i = 0; i < len; i++) {
		temp_late[i] = ~(uint64)0;
	}

	// Now generate q = 2kp+1 for as many k as desired; trial divide only if q passes
    // a small-primes sieve and if q mod 8 = +1 or -1 [Mersenne] or k even [Fermat; this is because we retrofit
    // the known form of Fn factors q = j.2^(n+2)+1 into our q = 2.k.p+1 scheme, replacing p with 2^n]:

	if(MODULUS_TYPE ==  MODULUS_TYPE_FERMAT) {
		/* No-odd-k sieve is simplest, so do it first.
		Remember that k = 1 is in the zeroth bit here, i.e. we cycle through
		bits 0-3 which correspond to q = 2kp+1, 4kp+1, 6kp+1, 8kp+1, respectively;
		thus we keep only bits 1 (k = 2) and 3 (k = 4), by masking with 16 x 0b1010 = 16 x 0xA:
		*/
		temp_late[0] &= 0xAAAAAAAAAAAAAAAAull;
	}
	else
	{
		/* q mod 8 = +-1 sieve is simplest, so do it first. Note q mod 8 = +-1 is guaranteed
		for (p mod TF_CLASSES) sieve via choice of acceptable values of increment (incr), but doing
		q%8 = +-1 here is trivial and speeds search for multiples of small primes below.
		*/
		qmod8 = 1;
		for(i = 0; i < 4; i++) {
			/* Remember that k = 1 is in the zeroth bit here, i.e. we cycle through
			bits 0-3 which correspond to q = 2kp+1, 4kp+1, 6kp+1, 8kp+1, respectively.

			We don't need to worry about overflow-on-add of (two_p + qmod8),
			since an overflow won't affect the result, modulo 8.
			*/
			qmod8=(two_p[0] + qmod8) & 7;
			if(qmod8==3 || qmod8==5) {
				for(l = i; l < 64; l += 4) {
					temp_late[0] &= ~((uint64)1 << l); /* clear every fourth such bit from the first quadword...	*/
				}
			}
		}
	}

	// Small-primes sieve needs only three (q mod 8 = +-1)-sieved registers to begin with, make 2 copies of the result:
	temp_late[1]=temp_late[0];
	temp_late[2]=temp_late[0];

/*...next, find multiples of small primes.	*/

	regs_todo=1;

	curr_p = 3;
	for(m = 0; m < nclear; m++)
	{
		curr_p += (pdiff[m] << 1);

		two64modp = 0x8000000000000000ull%curr_p;
		two64modp = (two64modp + two64modp)%curr_p;

		regs_todo=regs_todo*curr_p;

		mi64_set_eq_scalar(q, 1ull, lenQ);

		for(k = 0; k < regs_todo; k++)	/* number of registers to run through while seaching for first multiple of (m)th prime...	*/
		{
			for(i = 0; i < 64; i++)			/* ...and run through the 64 bits of each register. */
			{
				mi64_add(q,two_p,q,lenQ);
				if(!mi64_div_y32(q,curr_p,0x0,lenQ)) 	// To-Do: Replace all these slow mi64_div calls with mi64_is_div_by_scalar32,
				{										// then consider vectorizing this step, and using SSE2 vector version! */
				#if DBG_SIEVE
				//	if(curr_p < 100) printf("0: Found a multiple of %u in bit %u of register %u\n", curr_p, i, (uint32)k);
				#endif
					for(l = (uint32)(k<<6)+i; l < (regs_todo<<6); l += curr_p)
					{
						i64=(l>>6);	/* l/64	*/
						bit=l&63;	/* l%64	*/
						temp_late[i64] &= ~((uint64)1 << bit);
					}
					goto KLOOP;
				}
			}
		}
		/* Should never reach this regular-loop-exit point: */
		fprintf(stderr,"ERROR: failed to find a multiple of prime %u\n", curr_p);
		ASSERT(HERE, 0,"0");

	KLOOP:
		/* Propagate copies of length (regs_todo) bit-cleared portion of sieve to remaining parts of sieve.
		Here is the output of the commented-out diagnostic print:

			Propagating 5 copies of first 3 template-array words.
			Propagating 7 copies of first 15 template-array words.
			Propagating 11 copies of first 105 template-array words.
			Propagating 13 copies of first 1155 template-array words.
			Propagating 17 copies of first 15015 template-array words.

		Thus we really only need a template length = [product of first (nclear-1) off primes], since
		that simply repeats n = [(nclear)th odd prime] times. We simply need to make sure our loop below
		which copies the needed bits from the template to the bit_atlas also makes use of this repetition.
		*/
		if(m < (nclear-1)) {
			ncopies = prime[m+1];
			l = regs_todo;
		//	printf("Propagating %u copies of first %u template-array words.\n", ncopies,l);
			for(i = 2; i <= ncopies; i++, l += regs_todo) {
				for(j = 0; j < regs_todo; j++) {
					temp_late[l+j] = temp_late[j];
				}
			}
		}
	}	/* endfor(m = 0; m < nclear; m++) */

	for(m = 0; m < len; m++) {
		on_bits += popcount64(temp_late[m]);
	}
	printf("%u ones bits of %u in master sieve template.\n", on_bits, len<<6);

  #ifdef FACTOR_STANDALONE
	 printf(   "TRYQ = %u, max sieving prime = %u\n",TRYQ,MAX_SIEVING_PRIME);
  #else
	ASSERT(HERE, fp == 0x0,"0");
	fp = mlucas_fopen(STATFILE,"a");
	fprintf(fp,"TRYQ = %u, max sieving prime = %u\n",TRYQ,MAX_SIEVING_PRIME);
	fclose(fp); fp = 0x0;
  #endif

	/* Init bitmap in atlas for each of the [TF_PASSES] k mod TF_CLASSES cases:
	We advance the current-bit-to-copy index [i] in increments based on the length-TF_PASSES incr[] array.
	Copy each such bit into the current-bit ([bit]th) bit of the current ([word*TF_PASSES + l]th) bit_atlas word,
	then increment l (which determines which of the TF_PASSES sievelets we are in).

	Here are the resulting stats based on the number of classes:

	TF_CLASSES = 60:
	Copy 16 of every 60 bits - every [3.75]th bit on average - from the master template to the compact bit_atlas.
	Since template has 64 x 255255 = 16336320 bits, and 16336320/60 = 272272 bits for each of 16 sievelets concatenated
	into our bit_atlas, which thus has ceil(272272/64) = 4255 words, with high word [4254] using just its low 16 bits.

	TF_CLASSES = 4620:
	Copy 960 of every 60*77 bits - every [4.8125]th bit on average - from the master template to the compact bit_atlas.
	Since template has 64^2 x 255255 = 1045524480 bits (= 130690560 bytes [~130MB]), and 1045524480/4620 = 226304 bits
	for each of 960 sievelets concatenated into our bit_atlas, which thus has 226304/64 = 3536 words, with high word
	[3535] using all its bits, as 226304 is divisible by 64.

	Now, we prefer not to have to allocate the full-length template array unless absolutely necessary, but for
	the 4620 case we can simply re-use the 60-class version 64 times.
	*/
	i = incr[0] - 1;
	j = k = 0;
	bit = word = 0;	// To shut up may-be-uninit warnings
	ncopies = 0;
	for(;;) /* K loops over 64-bit registers...	*/
	{
		l = 0;
L3:
	#if TF_CLASSES == 60
		if(k == 0 && ncopies == 1) break;
	#else	// 4620 classes:
		if(k == 0 && ncopies == 64) break;
	#endif
		for(;;) /* I loops over bits...	*/
		{
			word = (uint32)(j>>6);
			bit  = (uint32)(j&63);
			if((temp_late[k]>>i) & 1) {
		//	if(k<10)printf("Copying set bit at k = %u in word %u of temp_late to word %u*%u + %u, bit %u of bit_atlas\n",i+1,(uint32)k, word,TF_PASSES,l,bit);
				*(bit_atlas + (word * TF_PASSES) + l) |= ((uint64)1 << bit);
			}
			if(++l == TF_PASSES) {
				l = 0;
				j++;
			}
			i += incr[l];
			// I will not generally land on 64 here, so need a more general mod-needed check:
			if(i >= 64) {
				i &= 63;
				k++;
				if(k == len) {
					k = 0;
					ncopies++;
				}
				goto L3;
			}
		} /* end of I loop	*/
	}	/* end of K loop	*/
//printf("L3: template word %u [used %u copies] bit_atlas chart %u, word %u, bit %u\n",(uint32)k,ncopies,l,word,bit);	exit(0);
	// For 60|4620 classes expect to end at bit 15|63 of the last word of each of the TF_PASSES = 16|960 sievelets (a.k.a. charts in our atlas):
	ASSERT(HERE, (k == 0) && (l == 0), "bit_atlas init: Exit check 1 failed!");
  #if TF_CLASSES == 60
	ASSERT(HERE, (word == 4254) && (bit == 15), "bit_atlas init: Exit check 2 failed!");
  #else	// 4620 classes:
	ASSERT(HERE, (word == 3535) && (bit == 63), "bit_atlas init: Exit check 2 failed!");
  #endif

  #ifdef FAC_DEBUG
  #if TF_CLASSES == 60
	i = len/TF_CLASSES + 1;	// len not divisible by TF_CLASSES, so add a pad-word
	j = i*64 - 48;	// #bits
  #else
	i = (len*64)/TF_CLASSES;	// 64*len divisible by TF_CLASSES, no need for padding
	j = i*64;	// #bits
  #endif
	l = 0;	// accum popc
	for(m = 0; m < i; m++) {
		l += popcount64(bit_atlas[m]);
	}
//	printf("%u ones bits of %u [%6.2f%%] in bit_atlas set.\n",l,(uint32)j,100.*(float)l/j);	exit(0);
//	  60:	184349 ones bits of 272272 [ 67.71%] in bit_atlas set.
//	4620:	196610 ones bits of 226304 [ 86.88%] in bit_atlas set.
  #endif

/*...deallocate full-sized bit_atlas.	*/
	free((void *)temp_late); temp_late = 0x0;

/*...At this point, replace the relative with the absolute increments:	*/
	for(i = 1; i < TF_PASSES; i++) {	// Skip pass 0 here
		incr[i] = incr[i-1] + incr[i];
	}
  #if TF_CLASSES == 60
	i = 0;
	switch(pmodNC)
	{
		/*   p mod 12 = 1:	*/
		case  1:ASSERT(HERE, incr[i++]== 3&&incr[i++]== 8&&incr[i++]==11&&incr[i++]==15&&incr[i++]==20&&incr[i++]==23&&incr[i++]==24&&incr[i++]==35&&incr[i++]==36&&incr[i++]==39&&incr[i++]==44&&incr[i++]==48&&incr[i++]==51&&incr[i++]==56&&incr[i++]==59&&incr[i++]==60, "factor.c : case  1"); break;	/* k mod 5 .ne. 2	*/
		case 37:ASSERT(HERE, incr[i++]== 3&&incr[i++]== 8&&incr[i++]==12&&incr[i++]==15&&incr[i++]==20&&incr[i++]==23&&incr[i++]==24&&incr[i++]==27&&incr[i++]==32&&incr[i++]==35&&incr[i++]==39&&incr[i++]==44&&incr[i++]==47&&incr[i++]==48&&incr[i++]==59&&incr[i++]==60, "factor.c : case 37"); break;	/* k mod 5 .ne. 1	*/
		case 13:ASSERT(HERE, incr[i++]== 3&&incr[i++]== 8&&incr[i++]==11&&incr[i++]==12&&incr[i++]==15&&incr[i++]==20&&incr[i++]==23&&incr[i++]==27&&incr[i++]==32&&incr[i++]==35&&incr[i++]==36&&incr[i++]==47&&incr[i++]==48&&incr[i++]==51&&incr[i++]==56&&incr[i++]==60, "factor.c : case 13"); break;	/* k mod 5 .ne. 4	*/
		case 49:ASSERT(HERE, incr[i++]==11&&incr[i++]==12&&incr[i++]==15&&incr[i++]==20&&incr[i++]==24&&incr[i++]==27&&incr[i++]==32&&incr[i++]==35&&incr[i++]==36&&incr[i++]==39&&incr[i++]==44&&incr[i++]==47&&incr[i++]==51&&incr[i++]==56&&incr[i++]==59&&incr[i++]==60, "factor.c : case 49"); break;	/* k mod 5 .ne. 3	*/
		/*   p mod 12 == 7:	*/
		case 31:ASSERT(HERE, incr[i++]== 5&&incr[i++]== 8&&incr[i++]== 9&&incr[i++]==20&&incr[i++]==21&&incr[i++]==24&&incr[i++]==29&&incr[i++]==33&&incr[i++]==36&&incr[i++]==41&&incr[i++]==44&&incr[i++]==45&&incr[i++]==48&&incr[i++]==53&&incr[i++]==56&&incr[i++]==60, "factor.c : case 31"); break;	/* k mod 5 .ne. 2	*/
		case  7:ASSERT(HERE, incr[i++]== 5&&incr[i++]== 8&&incr[i++]== 9&&incr[i++]==12&&incr[i++]==17&&incr[i++]==20&&incr[i++]==24&&incr[i++]==29&&incr[i++]==32&&incr[i++]==33&&incr[i++]==44&&incr[i++]==45&&incr[i++]==48&&incr[i++]==53&&incr[i++]==57&&incr[i++]==60, "factor.c : case  7"); break;	/* k mod 5 .ne. 1	*/
		case 43:ASSERT(HERE, incr[i++]== 5&&incr[i++]== 8&&incr[i++]==12&&incr[i++]==17&&incr[i++]==20&&incr[i++]==21&&incr[i++]==32&&incr[i++]==33&&incr[i++]==36&&incr[i++]==41&&incr[i++]==45&&incr[i++]==48&&incr[i++]==53&&incr[i++]==56&&incr[i++]==57&&incr[i++]==60, "factor.c : case 43"); break;	/* k mod 5 .ne. 4	*/
		case 19:ASSERT(HERE, incr[i++]== 5&&incr[i++]== 9&&incr[i++]==12&&incr[i++]==17&&incr[i++]==20&&incr[i++]==21&&incr[i++]==24&&incr[i++]==29&&incr[i++]==32&&incr[i++]==36&&incr[i++]==41&&incr[i++]==44&&incr[i++]==45&&incr[i++]==56&&incr[i++]==57&&incr[i++]==60, "factor.c : case 19"); break;	/* k mod 5 .ne. 3	*/
		/*   p mod 12 == 5:	*/
		case 41:ASSERT(HERE, incr[i++]== 3&&incr[i++]== 4&&incr[i++]==15&&incr[i++]==16&&incr[i++]==19&&incr[i++]==24&&incr[i++]==28&&incr[i++]==31&&incr[i++]==36&&incr[i++]==39&&incr[i++]==40&&incr[i++]==43&&incr[i++]==48&&incr[i++]==51&&incr[i++]==55&&incr[i++]==60, "factor.c : case 41"); break;	/* k mod 5 .ne. 2	*/
		case 17:ASSERT(HERE, incr[i++]== 3&&incr[i++]== 4&&incr[i++]== 7&&incr[i++]==12&&incr[i++]==15&&incr[i++]==19&&incr[i++]==24&&incr[i++]==27&&incr[i++]==28&&incr[i++]==39&&incr[i++]==40&&incr[i++]==43&&incr[i++]==48&&incr[i++]==52&&incr[i++]==55&&incr[i++]==60, "factor.c : case 17"); break;	/* k mod 5 .ne. 1	*/
		case 53:ASSERT(HERE, incr[i++]== 3&&incr[i++]== 7&&incr[i++]==12&&incr[i++]==15&&incr[i++]==16&&incr[i++]==27&&incr[i++]==28&&incr[i++]==31&&incr[i++]==36&&incr[i++]==40&&incr[i++]==43&&incr[i++]==48&&incr[i++]==51&&incr[i++]==52&&incr[i++]==55&&incr[i++]==60, "factor.c : case 53"); break;	/* k mod 5 .ne. 4	*/
		case 29:ASSERT(HERE, incr[i++]== 4&&incr[i++]== 7&&incr[i++]==12&&incr[i++]==15&&incr[i++]==16&&incr[i++]==19&&incr[i++]==24&&incr[i++]==27&&incr[i++]==31&&incr[i++]==36&&incr[i++]==39&&incr[i++]==40&&incr[i++]==51&&incr[i++]==52&&incr[i++]==55&&incr[i++]==60, "factor.c : case 29"); break;	/* k mod 5 .ne. 3	*/
		/*   p mod 12 == 11:	*/
		case 11:ASSERT(HERE, incr[i++]== 1&&incr[i++]== 4&&incr[i++]== 9&&incr[i++]==13&&incr[i++]==16&&incr[i++]==21&&incr[i++]==24&&incr[i++]==25&&incr[i++]==28&&incr[i++]==33&&incr[i++]==36&&incr[i++]==40&&incr[i++]==45&&incr[i++]==48&&incr[i++]==49&&incr[i++]==60, "factor.c : case 11"); break;	/* k mod 5 .ne. 2	*/
		case 47:ASSERT(HERE, incr[i++]== 4&&incr[i++]== 9&&incr[i++]==12&&incr[i++]==13&&incr[i++]==24&&incr[i++]==25&&incr[i++]==28&&incr[i++]==33&&incr[i++]==37&&incr[i++]==40&&incr[i++]==45&&incr[i++]==48&&incr[i++]==49&&incr[i++]==52&&incr[i++]==57&&incr[i++]==60, "factor.c : case 47"); break;	/* k mod 5 .ne. 1	*/
		case 23:ASSERT(HERE, incr[i++]== 1&&incr[i++]==12&&incr[i++]==13&&incr[i++]==16&&incr[i++]==21&&incr[i++]==25&&incr[i++]==28&&incr[i++]==33&&incr[i++]==36&&incr[i++]==37&&incr[i++]==40&&incr[i++]==45&&incr[i++]==48&&incr[i++]==52&&incr[i++]==57&&incr[i++]==60, "factor.c : case 23"); break;	/* k mod 5 .ne. 4	*/
		case 59:ASSERT(HERE, incr[i++]== 1&&incr[i++]== 4&&incr[i++]== 9&&incr[i++]==12&&incr[i++]==16&&incr[i++]==21&&incr[i++]==24&&incr[i++]==25&&incr[i++]==36&&incr[i++]==37&&incr[i++]==40&&incr[i++]==45&&incr[i++]==49&&incr[i++]==52&&incr[i++]==57&&incr[i++]==60, "factor.c : case 59"); break;	/* k mod 5 .ne. 3	*/
		default:
			ASSERT(HERE, MODULUS_TYPE == MODULUS_TYPE_FERMAT,"Only Mersenne and fermat-number factoring supported!");
	}
  #endif

	clock2 = clock();	/* Assume sieve setup time < 2^32 cycles - even if that
							that is violated it's no big deal at this point. */
	/* Use td here, as tdiff is reserved for the total runtime from factoring start: */
	// Accumulate the cycle count in a floating double to avoid problems with integer overflow
	// of the clock() result, if clock_t happens to be 32-bit int on the host platform:
  #ifdef CTIME
	clock2 = clock();
	td = (double)(clock2 - clock1);
	clock1 = clock2;
  #else
	clock2 = getRealTime();
	td = clock2 - clock1;
  #endif

  #ifdef FACTOR_STANDALONE
	if(!restart)
		printf("Time to set up sieve =%s\n",get_time_str(td));
  #endif

/* Run through each of the 16 "sievelets" as many times as necessary, each time copying
the appropriate q mod 8 and small-prime bit-cleared bit_atlas into memory, clearing bits
corresponding to multiples of the larger tabulated primes, and trial-factoring any
candidate factors that survive sieving.	*/

	nfactor = 0;

  #ifdef FAC_DEBUG
	/* If a known factor given, only process the given k/log2 range for that pass: */
	if(pass_targ < TF_PASSES) {
		passmin = passnow = passmax = pass_targ;
		printf("Setting run parameters to execute only the debug-targeted pass %u.\n",pass_targ);
	}
  #endif

  #if 0//def USE_GPU *** Doing this here gives 'cudaGetLastError() returned 36: cannot set while device is active in this process' -
			// This is because of the start-of-run GPU-self-tests in util.c; thus moved upstream to immediately precede those. ***
	#error Wrong place for this!
	// Disable default spin-loop-wait-for-GPU:
	cudaSetDeviceFlags(cudaDeviceBlockingSync);
  #endif

  #ifdef MULTITHREAD

//	printf("start; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 10000000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	// Populate the work-unit-encoding data structs which will get done by our pool of threads.
	// Current assignment may be restart of a partially-completed run, in which case npass < TF_PASSES
  #if TF_CLASSES == 60
	m = len/TF_CLASSES + 1;	// len not divisible by TF_CLASSES, so add a pad-word
  #else
	m = (len*64)/TF_CLASSES;	// 64*len divisible by TF_CLASSES, no need for padding
  #endif
	uint32 npass = passmax - passnow + 1;
	fprintf(stderr, "INFO: %u passes to do; bit_map has %u 64-bit words.\n",npass,m);
	/*
	Break overall work into ceil(npass/NTHREADS) 'waves' of NTHREADS each, i.e. each wave using all
	threads in the pool, with thread-synchronization at the end of each thread's factoring interval.
	We need to sync up (rather than just allowing each thread to grab more work as soon as it's done
	with its current (k mod) interval) because we memalloc for only NTHREADS concurrent threads, thus
	e.g. using NTHREADS = 4 if we start with passes 0-3 and the pass 0 and 2 threads finish ahead of 1 and 3,
	those 2 faster threads will get assigned the pass 4 and 5 work, and since 5 == 1 (mod NTHREADS), that
	thread's memory footprint will overlap that of the not-yet-complete pass 1 thread's.
	*/
	uint32 wave, wave_max = ceil( (double)npass/NTHREADS );
	fprintf(stderr,"INFO: Doing %u threadpool-waves of %u pool threads each:\n",wave_max,NTHREADS);
	for(wave = 0; wave < wave_max; wave++)
	{
		for(thr_id = 0; thr_id < NTHREADS; thr_id++)	// Unique (within the context of the current wave) task ID
		{												// used for slotting thread-local accesses to shared data arrays
			i = wave*NTHREADS + thr_id;
			pass = passnow + i;
			// For partial-waves, easiest is to proceed as usual, 'init' the full NTHREADS pool tasks,
			// but make the extra ones no-ops. But also need to avoid reading nonexistent bitmap data:
			if(pass <= passmax) {
				// Load 'master copy' of sievelet for the current pass number:
				for(l = 0; l < m; l++) {
					bit_map[l + m * thr_id] = *(bit_atlas + (l * TF_PASSES) + pass);
				}

				/* Starting no.-of-times-through-sieve = kmin/(64*len) : */
				if(pass == passnow && (know > kmin)) {
					interval_lo = know/((uint64)len << TF_CLSHIFT);
					ASSERT(HERE, know == interval_lo *(len << TF_CLSHIFT),"know == interval_lo*(len << TF_CLSHIFT)");
				} else {
					interval_lo = kmin/((uint64)len << TF_CLSHIFT);
					ASSERT(HERE, kmin == interval_lo *(len << TF_CLSHIFT),"kmin == interval_lo*(len << TF_CLSHIFT)");
				}
			} else {
				interval_lo = interval_hi;	// This is what defines a 'no-op' pool task.
			}
			/* Set initial k for this pass to default value (= incr[pass]) + interval_lo*(64*len),
			(assume this could be as large as 64 bits), then use it to set initial q for this pass:
			*/
			ASSERT(HERE, (double)interval_lo*(len << TF_CLSHIFT) < TWO64FLOAT, "(double)interval_lo*len < TWO64FLOAT");
			k = (uint64)incr[pass] + interval_lo*(len << TF_CLSHIFT);
		//	fprintf(stderr," [*** Init pass %u data: k0 = %llu, word0 = %16llX\n",pass,k,bit_map[0]);
			struct fac_thread_data_t* targ = tdat + thr_id;
			targ->count = &count;
			targ->tid = thr_id;		// Within the per-thread TFing, only the pool-thread ID matters
			targ->pass = pass;
			targ->interval_lo = interval_lo;
			targ->interval_hi = interval_hi;
			targ->fbits_in_2p = fbits_in_2p;
		#ifdef USE_AVX512
			targ->psmall = psmall;
		#endif
			targ->nclear = nclear;
			targ->sieve_len = len;
			targ->p_last_small = p_last_small;
			targ->nprime = nprime;
			targ->MAX_SIEVING_PRIME = MAX_SIEVING_PRIME;
			targ->pdiff = pdiff + NUM_SIEVING_PRIME * thr_id;
			targ->startval = startval + NUM_SIEVING_PRIME * thr_id;
			targ->k_to_try = k_to_try + TRYQ              * thr_id;
			targ->factor_k = factor_k + TRYQ              * thr_id;
			targ->nfactor = &nfactor;
			targ->findex = findex;
			targ->pstring = pstring;
			targ->p     = p;
			targ->two_p = two_p;
			targ->q       = q       + lenQ * thr_id;
			targ->q2      = q2      + lenQ * thr_id;
			targ->u64_arr = u64_arr + lenQ * thr_id;
			targ->lenP = lenP;
			targ->lenQ = lenQ;
			targ->kdeep = kdeep;
			targ->ndeep = &ndeep;
			targ->countmask = countmask;
			targ->CMASKBITS = CMASKBITS;
			targ->incr = incr[pass];
			targ->kstart = k;
			targ->bit_map  = bit_map  + m * thr_id;
			targ->bit_map2 = bit_map2 + m * thr_id;
			targ->tdiff = &tdiff;	// In || mode update tdiff directly, but only from the 0-thread
			targ->MODULUS_TYPE = MODULUS_TYPE;
			targ->VERSION      = VERSION;
			targ->OFILE        = OFILE;

		}	// thr_id-loop
		// Use exit value of [pass] here: If doing only a subset of the full 'current wave' complement
		// of NTHREADS passes - this can only occur during the final wave - adjust pool_work_units accordingly:
		if(pass > passmax) {
			pool_work_units = NTHREADS - (pass - passmax);	// Subtract excess #passes from default pool_work_units value
			fprintf(stderr,"INFO: Final threadpool wave will use only %u of the %u pool threads.\n",pool_work_units,NTHREADS);
		}

		if(pool_work_units > 1)
			fprintf(stderr, "Passes %u - %u: ",pass-NTHREADS+1, pass-NTHREADS+pool_work_units);
		else
			fprintf(stderr, "Pass %u: ",pass);

		// For partial-waves, easiest is to proceed as usual, 'init' the full NTHREADS pool tasks,
		// but make the extra ones no-ops. Here that means adding the full complement of NTHREADS tasks to the pool:
		for(thr_id = 0; thr_id < NTHREADS; ++thr_id)
		{
			task_control.data = (void*)(&tdat[thr_id]);
			threadpool_add_task(tpool, &task_control, task_is_blocking);
		#if 0
			printf("adding pool task %d with pool ID [%d]\n",thr_id,((struct thread_init *)(&task_control)->data)->thread_num);
			struct fac_thread_data_t* targ = tdat + thr_id;
			printf("This task has: pass %u, interval_[lo,hi] = [%llu,%llu]\n",targ->pass,targ->interval_lo,targ->interval_hi);
			printf("; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
		#endif
		}

		while(tpool->free_tasks_queue.num_tasks != NTHREADS) {
			// Posix sleep() too granular here; use finer-resolution, declared in <time.h>; cf. http://linux.die.net/man/2/nanosleep
			ASSERT(HERE, 0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
		}
		fprintf(stderr,"\n");	// For pretty-printing, have the inline-pass-printing reflect || work, newlines reflect sync-points
	};	// wave-loop

  #else	// Single-threaded execution:

	for(pass = passnow; pass <= passmax; pass++)
	{
		// Load 'master copy' of sievelet for the current pass number:
	  #if TF_CLASSES == 60
		m = len/TF_CLASSES + 1;	// len not divisible by TF_CLASSES, so add a pad-word
	  #else
		m = (len*64)/TF_CLASSES;	// 64*len divisible by TF_CLASSES, no need for padding
	  #endif
		for(i = 0; i < m; i++) {
			bit_map[i] = *(bit_atlas + (i * TF_PASSES) + pass);
		}

	#if DBG_SIEVE
		/* If debugging sieve, make sure critical bit hasn't been cleared: */
		if( k_targ && (((bit_map[i64_targ] >> bit_targ) & 1) == 0) ) {
			fprintf(stderr,"Critical bit cleared in master bitmap!\n");
			ASSERT(HERE, 0,"0");
		}
	#endif

	#ifdef FACTOR_STANDALONE
		if(!restart) {
			printf("pass = %u",pass);	fflush(stdout);
		}
	#else
		ASSERT(HERE, fp == 0x0,"0");
		fp = mlucas_fopen(STATFILE,"a");
		fprintf(fp,"Starting Trial-factoring Pass %2u...\n",pass);
		fclose(fp); fp = 0x0;
	#endif

		/* Starting no.-of-times-through-sieve = kmin/(64*len) : */
		if(pass == passnow && (know > kmin)) {
			interval_lo = know/((uint64)len << TF_CLSHIFT);
			ASSERT(HERE, know == interval_lo*((uint64)len << TF_CLSHIFT),"know == interval_lo*((uint64)len << TF_CLSHIFT)");
		} else {
			interval_lo = kmin/((uint64)len << TF_CLSHIFT);
			ASSERT(HERE, kmin == interval_lo*((uint64)len << TF_CLSHIFT),"kmin == interval_lo*((uint64)len << TF_CLSHIFT)");
		}

		/* Set initial k for this pass to default value (= incr[pass]) + interval_lo*(64*len),
		(assume this could be as large as 64 bits), then use it to set initial q for this pass:
		*/
		ASSERT(HERE, (double)interval_lo*(len << TF_CLSHIFT) < TWO64FLOAT, "(double)interval_lo*len < TWO64FLOAT");
		k = (uint64)incr[pass] + interval_lo*(len << TF_CLSHIFT);

		i = nprime;	// Remember, MAX_SIEVING_PRIME is a *variable* and set at runtime, as opposed to the predef NUM_SIEVING_PRIME;
					// And for small exponents, the actual #sieving prime is in nprime, and may be < NUM_SIEVING_PRIME.
		count += PerPass_tfSieve(
			pstring,
			pass,
			interval_lo,  interval_hi,
			fbits_in_2p,
			nclear,
			len,	// #64-bit words in the full-length sieving bitmap setup by the (CPU-side) caller.
			p_last_small,	//largest odd prime appearing in the product; that is, the (nclear)th odd prime.
			i,	// #sieving primes (counting from 3)
			MAX_SIEVING_PRIME,
			pdiff,
			startval,	// This gets updated within
			k_to_try,	// Unused by GPU code; include to yield a (mostly) uniform API
			factor_k,	// List of found factors for each p gets updated (we hope) within
			&nfactor,	// Here the '*' is to denote a writeable scalar
			findex,
			p, two_p, lenP, lenQ, incr[pass],
		#ifndef USE_GPU
 			kdeep, &ndeep, countmask, CMASKBITS, q, q2, u64_arr,
		#endif
			k,
			bit_map, bit_map2,	// GPU version uses only 1st of these
			&tdiff,
			MODULUS_TYPE,
			VERSION,
			OFILE
		);

	#ifdef USE_GPU
		cudaError = cudaGetLastError();
		if(cudaError != cudaSuccess)
		{
			printf("ERROR: cudaGetLastError() returned %d: %s\n", cudaError, cudaGetErrorString(cudaError));
			ASSERT(HERE, 0, "factor.c : GPU-side error detected!");
		}
	#endif

	#ifdef FACTOR_STANDALONE
		if(!restart)
			printf("\n");
	#else
		fp = mlucas_fopen(STATFILE,"a");
		fprintf(fp,"Trial-factoring Pass %2u: time =%s\n", pass,get_time_str(td));
		fclose(fp); fp = 0x0;
	#endif

	/***********/
	#ifdef ONEPASS
		return 0;
	#endif
	/***********/
	}	/* end of pass loop	*/

  #endif	// MULTITHREAD ?

/*...all done.	*/
  #ifdef FACTOR_STANDALONE
	if(!restart)
	{
		printf(   "%s(%s) has %u factors in range k = [%llu, %llu], passes %u-%u\n",
	 	NUM_PREFIX[MODULUS_TYPE], pstring, nfactor, kmin, kmax, passmin, passmax);
		printf(   "Performed %s trial divides\n", &char_buf0[convert_uint64_base10_char(char_buf0, count)]);
		/* Since we're done accumulating cycle count, divide to get total time in seconds: */
		printf(   "Clocks =%s\n",get_time_str(tdiff));
	}
  #else
	ASSERT(HERE, fp == 0x0,"0");
	fp = mlucas_fopen(STATFILE,"a");
	fprintf(fp,"Performed %s trial divides\n", &char_buf0[convert_uint64_base10_char(char_buf0, count)]);
	/* Since we're done accumulating cycle count, divide to get total time in seconds: */
	fprintf(fp,"Clocks =%s\n",get_time_str(tdiff));
	fclose(fp); fp = 0x0;
  #endif

	fp = mlucas_fopen(   OFILE,"a");
  #ifdef P1WORD
	 fprintf(fp,"M(%s) has %u factors in range k = [%llu, %llu], passes %u-%u\n", pstring, nfactor, kmin, kmax, passmin, passmax);
  #else
	 fprintf(fp,"M(%s) has %u factors in range k = [%llu, %llu], passes %u-%u\n", pstring, nfactor, kmin, kmax, passmin, passmax);
  #endif
	fclose(fp); fp = 0x0;

  #ifdef FAC_DEBUG
	/* If a test factor was given, make sure we found at least one factor: */
	if(k_targ > 0)
	{
		ASSERT(HERE, nfactor > 0,"k_targ > 0 but failed to find at least one factor");
	}
  #endif

	// If in double-Mersenne deep-sieve mode, print sorted list of k-to-do:
	if((MODULUS_TYPE == MODULUS_TYPE_MERSMERS) && (findex > 10000)) {
		if(ndeep > 0) {
			qsort(kdeep, ndeep, sizeof(uint32), ncmp_uint32);
			printf("MM(%u): Do deep sieving for k = ",findex);
			for(i = 0; i < ndeep-1; i++) {
				printf("%u,",kdeep[i]);
			}
			printf("%u\n",kdeep[i]);
		}
	}

	// Free the allocated memory:
	free((void *)factor_ptmp);
	free((void *)p);
	free((void *)kdeep);
	free((void *)bit_map);
	free((void *)bit_map2);
	free((void *)bit_atlas);
	free((void *)pdiff);
	free((void *)startval);
	free((void *)pinv);
	free((void *)two_p);
	free((void *)p2NC);
	free((void *)q);
	free((void *)q2);
	free((void *)u64_arr);
  #ifdef MULTITHREAD
	free((void *)tdat); tdat = 0x0;
  #endif

	return(0);

	/* Only reachable from argc/argv section: */
  #ifdef FACTOR_STANDALONE
MFACTOR_HELP:
	printf(" Mfactor command line options ...\n");
	printf(" <CR>        Default mode: prompts for manual keyboard entry\n");
	printf("\n");
	printf(" -h          Prints this help file and exits\n");
	printf("\n");
	printf(" -m {num}    Trial-factor the Mersenne number M(num) = 2^num - 1.\n");
	printf("\n");
	printf(" -mm {num}   Trial-factor the double-Mersenne number M(M(num)) = 2^(2^num) - 1.\n");
	printf("\n");
	printf(" -f {num}    Trial-factor the Fermat number F(num) = 2^(2^num) + 1.\n");
	printf("\n");
	printf(" -file {string}    Name of checkpoint file (needed for restart-from-interrupt)\n");
	printf("\n");
  #ifdef P1WORD
	printf(" -bmin {num} Log2(minimum factor to try), in floating double form.\n");
	printf(" If > 10^9 its whole-number part is taken as the kmin value instead.\n");
	printf("\n");
	printf(" -bmax {num} Log2(maximum factor to try), in floating double form.\n");
	printf(" If > 10^9 its whole-number part is taken as the kmax value instead.\n");
	printf("\n");
  #endif
	printf(" -kmin {num}  Lowest factor K value to be tried in each pass ( > 0).\n");
	printf("\n");
	printf(" -kmax {num} Highest factor K value to be tried in each pass ( < 2^64).\n");
	printf("\n");
	printf(" -passmin {num}  Current factoring pass (0-%d).\n",TF_PASSES-1);
	printf("\n");
	printf(" -passmax {num}  Maximum pass for the run (0-%d).\n",TF_PASSES-1);
  #ifdef MULTITHREAD
	printf("\n");
	printf(" -nthread {num}  Number of threads to use (1-%u). Each pass gets done by\n\t\t a separate thread; if #passes > #threads, some threads will do multiple passes.\n",TF_PASSES);
  #endif
	/* If we reached here other than via explicit invocation of the help menu, assert: */
	if(!STREQ(stFlag, "-h"))
		ASSERT(HERE, 0,"Mfactor: Unrecognized command-line option!");
	return(0);
  #endif
}

#endif	/* #ifdef FACTOR_STANDALONE */

/******************/

#ifndef USE_GPU

  #ifndef MULTITHREAD

	uint64 PerPass_tfSieve(
		const char *pstring,
		const uint32 pass,
		const uint64 interval_lo, const uint64 interval_hi,
		const double fbits_in_2p,
		const uint32 nclear,
		const uint32 sieve_len,
		const uint32 p_last_small,	//largest odd prime appearing in the product; that is, the (nclear)th odd prime.
		const uint32 nprime,	// #sieving primes (counting from 3)
		const uint32 MAX_SIEVING_PRIME,
		const uint8 *pdiff,
			  uint32*startval,
			  uint64*k_to_try,
			  uint64*factor_k,	// List of found factors for each p gets updated (we hope) within
			  uint32*nfactor,	// Here the '*' is to denote a writeable scalar
		const uint32 findex,
		const uint64*p, const uint64*two_p, const uint32 lenP, const uint32 lenQ, const uint32 incr,
		uint32*kdeep, uint32*ndeep, const uint64 countmask, const uint32 CMASKBITS,
		uint64*q, uint64*q2, uint64*u64_arr,
		const uint64 kstart,
		const uint64*bit_map, uint64*bit_map2,
		double *tdiff,
		const int MODULUS_TYPE,
		const char*VERSION,
		const char*OFILE
	) {
		int    tid = 0;

  #else

	void*
	PerPass_tfSieve(void*thread_arg)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct fac_thread_data_t* targ = thread_arg;	// Ref'd as task->data in threadpool.c::worker_thr_routine() caller
		int    tid          = targ->tid;	// Thread ID (Use the pool-thread ID here rather than the task ID ... there are typically many more tasks than pool threads)
		uint32 pass         = targ->pass;	// Which factoring pass is the thread doing?
		uint64 interval_lo  = targ->interval_lo;
		uint64 interval_hi  = targ->interval_hi;
		double fbits_in_2p  = targ->fbits_in_2p;
	#ifdef USE_AVX512
		uint32 *psmall = targ->psmall;
	#endif
		uint32 nclear       = targ->nclear;
		uint32 sieve_len    = targ->sieve_len;
		uint32 p_last_small = targ->p_last_small;	//largest odd prime appearing in the product; that is = targ->; the (nclear)th odd prime.
		uint32 nprime       = targ->nprime;			// #sieving primes (counting from 3)
		uint32 MAX_SIEVING_PRIME = targ->MAX_SIEVING_PRIME;
		uint8 *pdiff        = targ->pdiff;
		uint32*startval     = targ->startval;
		uint64*k_to_try     = targ->k_to_try;
		uint64*factor_k     = targ->factor_k;		// List of found factors for each p gets updated (we hope) within
		uint32*nfactor      = targ->nfactor;		// Here the '*' is to denote a writeable scalar
		uint32 findex       = targ->findex;
		const char *pstring = targ->pstring;
		uint64*p            = targ->p;
		uint64*two_p        = targ->two_p;
		uint64*q            = targ->q;
		uint64*q2           = targ->q2;
		uint64*u64_arr      = targ->u64_arr;
		uint32 lenP         = targ->lenP;
		uint32 lenQ         = targ->lenQ;
		uint32*kdeep        = targ->kdeep;
		uint32*ndeep        = targ->ndeep;
		uint64 countmask    = targ->countmask;
		uint32 CMASKBITS    = targ->CMASKBITS;
		uint32 incr         = targ->incr;
		uint64 kstart       = targ->kstart;
		uint64*bit_map      = targ->bit_map ;
		uint64*bit_map2     = targ->bit_map2;
		double *tdiff       = targ->tdiff;
		int    MODULUS_TYPE = targ->MODULUS_TYPE;
		const char *VERSION = targ->VERSION;
		const char *OFILE   = targ->OFILE;
  #endif
		int found_pass = FALSE;
	#ifdef MULTITHREAD
		// Proper init (as opposed to no-init) key to avoiding deadlock here.
		// Started with 2 separate _checkpoint and _foundfactor mutexes here, but since both code sections
		// in question call some of the same mi64 functions, replaced with 'one mutex to rule them all' model:
		pthread_mutex_t mutex_mi64        = PTHREAD_MUTEX_INITIALIZER,
						mutex_updatecount = PTHREAD_MUTEX_INITIALIZER;	// No mi64 calls here.
	#endif
		FILE *fp = 0x0;
		char *char_addr;
	#if TF_CLASSES == 60
		const uint32 TRYQM1 = TRYQ-1, bit_len = (sieve_len << TF_CLSHIFT)/TF_CLASSES; 	// 255255*64  /  60 = 272272: Number of bits in each of the  16 mod-  60 sievelets
	#else	// 4620 classes:
		const uint32 TRYQM1 = TRYQ-1, bit_len = (sieve_len << TF_CLSHIFT)/TF_CLASSES;	// 255255*64^2/4620 = 226304: Number of bits in each of the 960 mod-4620 sievelets
	#endif
		int itmp;
		uint32 bit,bit_hi,curr_p,i,ihi,idx,j,l,m;
		uint64 count = 0ull, itmp64, k = 0ull, sweep, res;
		int32 q_index = -1;
		double fbits_in_k = 0, fbits_in_q = 0;
	#ifdef P1WORD
		double twop_float = 0,fqlo,fqhi;
		uint128 p128,q128,t128;	// Despite the naming, these are needed for nominal 1-word runs with moduli exceeding 64 bits
	#endif
	#ifdef P3WORD
		uint192 p192,q192,t192;
	  #ifdef USE_FLOAT
		uint256 x256;	// Needed to hold result of twopmodq200_8WORD_DOUBLE
	  #endif
	#endif
	#ifdef P4WORD
		uint256 p256,q256,t256;
	#endif
		char cbuf[STR_MAX_LEN*2], cbuf2[STR_MAX_LEN*2];
	#ifdef CTIME
		clock_t clock1, clock2;
	#else	// Multithreaded needs wall-clock, not CPU time:
		double clock1, clock2;	// Jun 2014: Switched to getRealTime() code
	#endif

		if(interval_lo == interval_hi) {
			printf("Thread %u immediate-return (no-op)\n",tid);
			return 0x0;
		}

	#if 0	/************** disable for now - need to sync with similar code in main() ***************/
	#error Need to sync this code with similar code in main()!
	// TF restart files are in HRF, not binary - checkpointing only supported for single-threaded runs:
	if(NTHREADS == 1)
	{
		fp = mlucas_fopen(RESTARTFILE, "r");
		if(fp) {
			// If file exists, it should have the proper first 2 lines:
			itmp = fscanf(fp,"%s\n",cstr);
			if(itmp <= 0 || !STREQ(cstr,pstring)) {
				sprintf(char_buf0,"Line 1 entry found in factoring savefile [%s] does not match exponent of run [%s].",cstr,pstring);
				ASSERT(HERE,0,char_buf0);
			}
			itmp = fscanf(fp,"%u\n",&i  );
			if(itmp <= 0 || i != TF_PASSES      ) {
				sprintf(char_buf0,"Line 1 entry found in factoring savefile [%d] does not match exponent of run [%d].",i,TF_PASSES);
				ASSERT(HERE,0,char_buf0);
			}
			// See if restart file has a pass/max-k-reached entry matching the current pass:
			while(fgets(cstr,STR_MAX_LEN,fp)) {
				if((char_addr = strstr(cstr,"Pass ")) != 0) {
					itmp = sscanf(char_addr,"%u",i);
					if(itmp <= 0) {
						fprintf(stderr,"ERROR: unable to read [Pass *: k] entry: offending line = [%s]\n",cstr); ASSERT(HERE, 0,"0");
					}
					if(i == pass) {	// Is the pass index the one we are updating? If yes, update the k-value
						ASSERT(HERE, !found_pass, "Multiple current-pass entry found in savefile!");
						found_pass = TRUE;
						// Read the max-k-reached value
						ASSERT(HERE,((char_addr = strstr(cstr,"Pass ")) != 0),"Expected : following pass number not found!");
						itmp = sscanf(char_addr,"%llu",k);
						ASSERT(HERE,itmp >= 0,"Unable to read max-k-reached value!");
						// Even if valid entry found, process rest of file to ensure no duplicate-pass-number entries
					}
				}
			}
			/* pstring*/
			++curr_line;
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				fprintf(stderr,"ERROR: unable to read Line %d (current exponent) of factoring restart file %s!\n", curr_line, RESTARTFILE);		ASSERT(HERE, 0,"0");
			}
			/* Strip the expected newline char from in_line: */
			char_addr = strstr(in_line, "\n");
			if(char_addr)
				*char_addr = '\0';
			/* Make sure restart-file and current-run pstring match: */
			if(STRNEQ(in_line, pstring)) {
				fprintf(stderr,"ERROR: current exponent %s != Line %d of factoring restart file %s!\n",pstring, curr_line, RESTARTFILE);		ASSERT(HERE, 0,"0");
			}

			/* bmin */
			++curr_line;
			fgets(cbuf, STR_MAX_LEN*2, fp);
			itmp = sscanf(cbuf, "%lf", &bmin_file);
			if(itmp != 1) {
				fprintf(stderr,"ERROR: unable to parse Line %d (bmin) of factoring restart file %s. Offending input = %s\n", curr_line, RESTARTFILE, cbuf);		ASSERT(HERE, 0,"0");
			}

			/* bmax */
			++curr_line;
			fgets(cbuf, STR_MAX_LEN*2, fp);
			itmp = sscanf(cbuf, "%lf", &bmax_file);
			if(itmp != 1) {
				fprintf(stderr,"ERROR: unable to parse Line %d (bmin) of factoring restart file %s. Offending input = %s\n", curr_line, RESTARTFILE, cbuf);		ASSERT(HERE, 0,"0");
			}

		/************************************
		LINE PAIRS 5/6 AND 7/8 ARE USED TO DETERMINE WHETHER A PREVIOUS
		FACTORING RUN OF THE SAME EXPONENT COMPLETED OR NOT: If know >= kmax
		and passnow = passmax then the previous run completed, in which case
		we allow a new run to a deeper bound, i.e. reset passnow = passmin
		and run passes passmin through passmax from bounds kmin to kmax.
		*************************************/
			/* KMin */
			++curr_line;
	GET_LINE4:
		/**** redo this ****/
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				fprintf(stderr,"ERROR: 'KMin' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			}
			char_addr = strstr(in_line, "KMin");
			/* Since the preceding fscanf call may leave us at the end of curr_line-1
			(rather than the beginning of curr_line), allow for a possible 2nd needed
			fgets call here: */
			if(!char_addr) {
				goto GET_LINE4;
			} else {
				char_addr = strstr(in_line, "=");
				if(!char_addr) {
					fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
				}
				char_addr++;
				kmin_file = convert_base10_char_uint64(char_addr);
			}

			/* KNow */
			++curr_line;
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				fprintf(stderr,"ERROR: unable to read Line %d (KNow) of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			}
			char_addr = strstr(in_line, "KNow");
			if(!char_addr) {
				fprintf(stderr,"ERROR: 'KNow' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			} else {
				char_addr = strstr(in_line, "=");
				if(!char_addr) {
					fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
				}
				char_addr++;
				know_file = convert_base10_char_uint64(char_addr);
			}

			/* KMax */
			++curr_line;
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				fprintf(stderr,"ERROR: unable to read Line %d (KMax) of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			}
			char_addr = strstr(in_line, "KMax");
			if(!char_addr) {
				fprintf(stderr,"ERROR: 'KMax' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			} else {
				char_addr = strstr(in_line, "=");
				if(!char_addr) {
					fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
				}
				char_addr++;
				kmax_file = convert_base10_char_uint64(char_addr);
			}

			/* PassMin */
			++curr_line;
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				fprintf(stderr,"ERROR: unable to read Line %d (PassMin) of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			}
			char_addr = strstr(in_line, "PassMin");
			if(!char_addr) {
				fprintf(stderr,"ERROR: 'PassMin' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			} else {
				char_addr = strstr(in_line, "=");
				if(!char_addr) {
					fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
				}
				char_addr++;
				passmin_file = (uint32)convert_base10_char_uint64(char_addr);
				ASSERT(HERE, passmin_file < TF_PASSES,"factor.c: passmin < TF_PASSES");
			}

			/* PassNow */
			++curr_line;
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				fprintf(stderr,"ERROR: unable to read Line %d (PassNow) of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			}
			char_addr = strstr(in_line, "PassNow");
			if(!char_addr) {
				fprintf(stderr,"ERROR: 'PassNow' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			} else {
				char_addr = strstr(in_line, "=");
				if(!char_addr) {
					fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
				}
				char_addr++;
				passnow_file = (uint32)convert_base10_char_uint64(char_addr);
				ASSERT(HERE, passnow_file < TF_PASSES,"factor.c: passnow < TF_PASSES");
				ASSERT(HERE, passnow_file >= passmin_file  ,"factor.c: passnow_file >= passmin_file");
			}

			/* PassMax */
			++curr_line;
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				fprintf(stderr,"ERROR: unable to read Line %d (PassMax) of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			}
			char_addr = strstr(in_line, "PassMax");
			if(!char_addr) {
				fprintf(stderr,"ERROR: 'PassMax' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			} else {
				char_addr = strstr(in_line, "=");
				if(!char_addr) {
					fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
				}
				char_addr++;
				passmax_file = (uint32)convert_base10_char_uint64(char_addr);
				ASSERT(HERE, passmax_file < TF_PASSES,"factor.c: passmax_file < TF_PASSES");
				ASSERT(HERE, passmax_file >= passnow_file  ,"factor.c: passmax_file >= passnow_file");
			}

			/* Number of q's tried: */
			++curr_line;
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				fprintf(stderr,"ERROR: unable to read Line %d (#Q tried) of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			}
			char_addr = strstr(in_line, "#Q tried");
			if(!char_addr) {
				fprintf(stderr,"ERROR: '#Q tried' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
			} else {
				char_addr = strstr(in_line, "=");
				if(!char_addr) {
					fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);	ASSERT(HERE, 0,"0");
				}
				char_addr++;
				count = convert_base10_char_uint64(char_addr);	// Need to reset == 0 prior to sieving so kvector-fill code works properly
			}
			fclose(fp); fp = 0x0;
		}	// endif(NTHREADS == 1)

		/**** process restart-file and any command-line params: ****/

		/* If previous run is not yet complete, ignore any increased factor-bound-related
		command-line parameters and instead proceed to complete the previous run first:
		*/
		if((know_file < kmax_file) || (passnow_file < passmax_file))
		{
			incomplete_run = TRUE;

			fprintf(stderr,"INFO: Previous run to kmax = %s not yet complete.\n"  , &char_buf0[convert_uint64_base10_char(char_buf0, kmax_file)]);
			fprintf(stderr,"Ignoring any increased factor-bound-related command-line parameters and proceeding to complete previous run.\n");

			bmin = bmin_file;
			bmax = bmax_file;

			passmin = passmin_file;
			passnow = passnow_file;
			passmax = passmax_file;

			kmin = kmin_file;
			know = know_file;
			kmax = kmax_file;
			kplus = 0;
		}
		else
		{
			/**** Previous run was completed - check that current params satisfy one (and only one)
			of the following sets of conditions:

				1) -bmin/bmax used to set bounds for factoring:
					In this case we expect any command-line bmin will be >= that in the restart file
					(in fact we expect bmin >= bmax_file, i.e. that the runs are nonoverlapping -
					if not we warn and set bmin = bmax_file), and that bmax > bmax_file.
			****/
			if(bmin || bmax)
			{
			#if(!defined(P1WORD))
			//	ASSERT(HERE, 0,"bmin/bmax form of bounds-setting only allowed for single-word-p case!");
			#endif
				ASSERT(HERE, (kmin==0 && kmax==0 && kplus==0),"(kmin==0 && kmax==0 && kplus==0) - please delete any restart files for this p and retry debug run.");

				if(bmin) {
					ASSERT(HERE, bmin >= bmin_file - 0.0000000001,"bmin >= bmin_file");
					if(bmin < bmax_file)
						fprintf(stderr,"WARNING: Specified bmin (%lf) smaller than previous-run bmax = %lf. Setting equal to avoid overlapping runs.\n", bmin, bmax_file);
				}
				bmin = bmax_file;

				/* We expect any command-line bmax will be > that in the restart file: */
				if(bmax)
					ASSERT(HERE, bmax > bmax_file - 0.0000000001,"bmax >= bmax_file");
			}

			/****
				2) -kmin/kmax used to set bounds for factoring:
					In this case we expect any command-line kmin will be >= that in the restart file
					(in fact we expect kmin >= kmax_file, i.e. that the runs are nonoverlapping -
					if not we warn and set kmin = kmax_file), and that kmax > kmax_file.
			****/
			if(kmin || kmax)
			{
				ASSERT(HERE, (bmin==0 && bmax==0 && kplus==0),"(bmin==0 && bmax==0 && kplus==0)");

				if(kmin) {
					ASSERT(HERE, kmin >= kmin_file,"kmin >= kmin_file");
					if(kmin < kmax_file)
						fprintf(stderr,"WARNING: Specified kmin (%s) smaller than previous-run kmax = %s. Setting equal to avoid overlapping runs.\n", &char_buf0[convert_uint64_base10_char(char_buf0, kmax)], &char_buf1[convert_uint64_base10_char(char_buf1, kmax_file)]);
				}
				kmin = kmax_file;

				/* We expect any command-line kmax will be > that in the restart file: */
				if(kmax)
					ASSERT(HERE, kmax > kmax_file,"kmax >= kmax_file");
			}

			/****
				3) -kplus used to increment an upper bound from a previous factoring run:
			****/
			if(kplus)
			{
				ASSERT(HERE, (bmin==0 && bmax==0 && kmin==0 && kmax==0),"(bmin==0 && bmax==0 && kmin==0 && kmax==0)");

				kmin = kmax_file;
				/* Ensure incremented value kmax fits into a 64-bit unsigned int: */
				ASSERT(HERE, (kmin + kplus) > kplus, "kmax_file + kplus exceeds 2^64!");
				kmax = kmin + kplus;
				kplus = 0;	/* If kplus != 0 detected further on, that indicates that no valid restart
							file was found for factoring-bounds incrementing. */
			}
		}
		/* Successfully processed restart file: */
		restart = TRUE;
	}

	// No restart file found, or no enty for this pass found in same, or max-k-reached-for-this-pass found therein <= kstart;
	// Take starting k value from kstart:
	if(!k)
		k = kstart;
	#endif	// #if 0
/************************ END(RESTART STUFF) *******************/

	#ifdef MULTITHREAD
	//	fprintf(stderr, "In PerPass_tfSieve task_id = %u, worker thread id %u\n", tid, ((struct thread_init *)targ)->thread_num);
	//	if(interval_hi > interval_lo) fprintf(stderr, "pass = %u",pass);	// Only print this diagnostic for non-empty tasks
	#endif
		// In || mode, only the 0-thread accumulates runtime, but do this for all threads to avoid uninit warnings:
	#ifdef CTIME
		clock1 = clock();
	#else
		clock1 = getRealTime();
	#endif

	  #ifdef FAC_DEBUG
		// compute qstart = 2.kstart.p + 1:
		ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,q,lenQ), "2.k.p overflows!");
		q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
		printf(" Initial q for this pass = %s.\n", &char_buf0[convert_mi64_base10_char(char_buf0, q, lenQ, 0)]);
	  #endif
//if(pass==4)
//	printf("\nPass %u: k0 = %u, word0 prior to deep-prime clearing = %16llX\n",pass,(uint32)kstart,bit_map[0]);
		// Compute startbit k (occurrence of first multiple of prime curr_p in first pass through the relevant sievelet:
		if((lenP == 1) && (p[0] <= MAX_SIEVING_PRIME))
			get_startval(MODULUS_TYPE, p[0], findex, two_p, lenQ, bit_len, interval_lo, incr, nclear, nprime, p_last_small, pdiff, startval);
		else
			get_startval(MODULUS_TYPE, 0ull, findex, two_p, lenQ, bit_len, interval_lo, incr, nclear, nprime, p_last_small, pdiff, startval);

		for(sweep = interval_lo; sweep < interval_hi; ++sweep)
		{
#ifdef MULTITHREAD
//if(tid == 0)
#endif
//	printf("sweep %llu: k0 = %llu, count %llu: k0-3 = %llu,%llu,%llu,%llu\n",sweep,kstart,count,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

			/* Accumulate the cycle count every so often to avoid problems with integer overflow
			of the clock() result, if clock_t happens to be a 32-bit int type on the host platform:
			*/
			if((tid == 0) && sweep && (sweep & 127) == 0) {	// Only the lowest-index thread does screen status-update prints
			#ifdef FACTOR_STANDALONE
				printf(".");	fflush(stdout);
			#endif
			}

		  #if TF_CLASSES == 60
			ihi = sieve_len/TF_CLASSES + 1;	// sieve_len not divisible by TF_CLASSES, so add a pad-word
		  #else
			ihi = (sieve_len*64)/TF_CLASSES;	// 64*sieve_len divisible by TF_CLASSES, no need for padding
		  #endif
			memcpy(bit_map2, bit_map, (ihi<<3));	// Load fresh copy of master sievelet

		#ifdef FAC_DEBUG	// If enable this, make sure to also uncomment "# survived" complement below!
			m = 0;	// accum popc
			for(i = 0; i < ihi; i++) {
				m += popcount64(bit_map2[i]);
			}
		  #ifdef MULTITHREAD
			if(tid == 0)
		  #endif
			if(sweep == 0) printf("%u ones bits of %u [%6.2f%%] in bit_map set.\n",m,bit_len,100.*(float)m/bit_len);
		#endif

			// To track lg(q) = lg(2.k.p+1), use approximation q ~= 2.k.p, thus lg(q) ~= lg(2.p) + lg(k).
			// At start of each pass through the k-based sieve, use 2.k.p with k = [starting k + sievebits]
			// to bound qmax from above, and compute lg(qmax) using the above logarithmic sum.
			fbits_in_k = log((double)k + TF_CLASSES*bit_len)*ILG2;	// Use k-value at end of upcoming pass thru sieve as upper-bound
			fbits_in_q = fbits_in_2p + fbits_in_k;
	//	if(fbits_in_q > 64)
	//	printf("sweep = %llu: fbits_in_q = fbits_in_2p [%10.4f] + fbits_in_k [%10.4f] = %10.4f\n",sweep,fbits_in_2p,fbits_in_k,fbits_in_q);

		/*********************************************/
		#if DBG_SIEVE
			if(k_targ)
			{
				/* See if k_targ falls in the range of k's for the current sieve interval: */
				k = (uint64)incr + sweep*(sieve_len<<6);	/* Starting k for the current interval: */

				/* If so, calculate the location of the critical bit
				and allow execution to proceed as normal to bit-clearing step:
				*/
				if((k <= k_targ) && (k_targ < (k+(sieve_len<<6))))
				{
					itmp64 = k_targ - k;
					ASSERT(HERE, itmp64%TF_CLASSES == 0,"(k_targ - k)%TF_CLASSES == 0");
					itmp64 /= TF_CLASSES;
					i64_targ = itmp64 >> 6;
					bit_targ = itmp64 & 63;
				}
				/* If debugging sieve & known factor will not occur in the current sieve interval,
				skip all of that bit-clearing nonsense and just increment the offset for each
				sieving prime to what it will be after another full pass through the sievelet:
				*/
				else
				{
					curr_p = p_last_small;
					for(m = nclear; m < nprime; m++)
					{
						curr_p += (pdiff[m] << 1);

						/* Special-handling code for p == curr_p case: */
						if(startval[m] == 0xffffffff)
						{
							continue;
						}

						/* Need to account for the fact that primes greater than the # of bits in the sieve
						may not hit *any* sieve bit on the current pass through the sieve:
						*/
						if(startval[m] >= bit_len) {
							startval[m] -= bit_len;
						} else {
							/* Compute new startvalue: */
							startval[m] += startval_incr[m] - curr_p;		/* Subtract off curr_p... */
							startval[m] += (-(startval_incr[m] >> 31)) & curr_p;	/* ...and restore it if result was < 0. */
						}
					}
					/* If reference factor not in current k-range, skip the TFing... */
					goto QUIT;
				}
			}
		#endif
		/*********************************************/

			/*   ...and clear the bits corresponding to the small primes.	*/

		#ifdef USE_AVX512	// Use vector-int math and gather-load/scatter-store to accelerate the bit-clearing
								// EWM: For pmax around the 'sweet spot', this 2-loop approach is barely faster than
								// above pure-C scalar-int code, though AVX-512 asm is a clear winner for large pmax.
			// Split our loop-over-primes into 2 parts, the 2nd of which handles primes > bit_len
			// [ = 272272 or 226304, resp., depending on whether TF_CLASSES = 60 or 4620].
			// We vectorize the 2nd loop, since each prime therein will hit at most one bit of the sievelet,
			// i.e. we require no while-loop, only an if(curr_p's startval < bit_len or not) conditional.
		// Loop #1:
			curr_p = p_last_small;
			for(m = nclear; m < nprime; m++)
			{
				curr_p += (pdiff[m] << 1);
				if(curr_p > bit_len && !((nprime - m)&63)) {	// 2nd clause is to make Loop #2 count a multiple of 64
					curr_p -= (pdiff[m] << 1);
					ASSERT(HERE, curr_p < p[0],"On Loop 1 exit: curr_p >= p!");
					break;
				}
				l = startval[m];
				// Need to account for the fact that primes greater than the # of bits in the sieve
				// may not hit *any* sieve bit on the current pass through the sieve:
				while(l < bit_len) {
					bit_clr32((uint32 *)bit_map2,l);	// 64-bit arithmetic offers no advamtage here
					l += curr_p;
				}
				startval[m] = l-bit_len;	// Save new startvalue
			}
		// Loop #2:
		//	printf("loop 2: m = %u, nprime = %u\n",m,nprime); exit(0);
		/*	for( ; m < nprime; m += 32) {	*/
			__asm__ volatile (\
			"	movl	%[__m] ,%%edx 	\n\t"\
			"	movq	%[__psmall]  ,%%rax 			\n\t	leaq (%%rax,%%rdx,4),%%rax	\n\t"/* &psmall[m] */\
			"	movq	%[__startval],%%rbx 			\n\t	leaq (%%rbx,%%rdx,4),%%rbx	\n\t"/* &startval[m] */\
			"	movq	%[__bit_map] ,%%rcx 			\n\t"\
			"	movl	$-2,%%edx						\n\t"\
			"	vpbroadcastd	%%edx,%%zmm31			\n\t"/* 0x111...110 x 16, shared across cols */
			"	movl	%[__bit_len] ,%%edx 			\n\t"\
			"	vpbroadcastd	%%edx,%%zmm30			\n\t"/* bit_len x 16, shared across cols */\
			"movl	%[__nprime] ,%%esi 	\n\t"/* ASM loop control structured as for(j = nprime-m; j != 0; j -= 32){...} */\
			"loop_bitclear_short:		\n\t"/* loop begin */\
				/* Each additional column to the right has zmm indices (except for data shared across cols, as noted) += 4, k-indices += 2: */\
				"	vmovups	(%%rbx),%%zmm0					\n\t	vmovups	0x40(%%rbx),%%zmm4				\n\t"/* next 16 startvals */\
				"	vpcmpd	$1,%%zmm30,%%zmm0,%%k1			\n\t	vpcmpd	$1,%%zmm30,%%zmm4,%%k3			\n\t"/* startval < bit_len ? If true, the current prime hits a bit in the sievelet. */\
				"	vmovdqa32 %%zmm0,%%zmm2%{%%k1%}			\n\t	vmovdqa32 %%zmm4,%%zmm6%{%%k3%}			\n\t"/* copy of 16 startvals, with elts which will not be touched set = 0. This is only to keep the resulting bitmap-dword-fetch index in range, though note the 0-word is 'live' data ... we will again use the same writemask to prevent the fetched bit_map[0] word from being modified. */\
				"	kmovw	%%k1,%%k2						\n\t	kmovw	%%k3,%%k4						\n\t"/* mask-reg zeroed by gather-load below, so save copy */\
				"	vpsrld	$5,%%zmm2,%%zmm2				\n\t	vpsrld	$5,%%zmm6,%%zmm6				\n\t"/* >>= 5 to convert startvals into dword indices */\
				"vpgatherdd (%%rcx,%%zmm2,4),%%zmm3%{%%k2%}	\n\t vpgatherdd (%%rcx,%%zmm6,4),%%zmm7%{%%k4%}	\n\t"/* gather-load the bitmap words */\
				"	vprolvd %%zmm0,%%zmm31,%%zmm1%{%%k1%}	\n\t	vprolvd %%zmm4,%%zmm31,%%zmm5%{%%k3%}	\n\t"/* 1 <<= bit (circular shift). Note no need to explicitly (mod 32) shift count */\
				"	vpandd	%%zmm1,%%zmm3 ,%%zmm3%{%%k1%}	\n\t	vpandd	%%zmm5,%%zmm7 ,%%zmm7%{%%k3%}	\n\t"/* AND with mask to clear the bits */\
				"	kmovw	%%k1,%%k2						\n\t	kmovw	%%k3,%%k4						\n\t"/* mask-reg zeroed by scatter-store below, so save copy */\
				"vpscatterdd %%zmm3,(%%rcx,%%zmm2,4)%{%%k2%}\n\t vpscatterdd %%zmm7,(%%rcx,%%zmm6,4)%{%%k4%}\n\t"/* scatter-store the bit-cleared bitmap words */\
				/* update the startvals - k1-mask determines which words get += curr_p; all words get -= bit_len */\
				"	vpaddd	(%%rax),%%zmm0,%%zmm0%{%%k1%}	\n\t vpaddd	0x40(%%rax),%%zmm4,%%zmm4%{%%k3%}	\n\t"/* startval{k1} += curr_p */\
				"	vpsubd	%%zmm30,%%zmm0,%%zmm0			\n\t	vpsubd	%%zmm30,%%zmm4,%%zmm4			\n\t"/* startval -= bit_len */\
				"	vmovups	%%zmm0,(%%rbx)					\n\t	vmovups	%%zmm4,0x40(%%rbx)				\n\t"/* write startvals */\
			"addq	$0x80,%%rax 	\n\t	addq	$0x80,%%rbx 	\n\t"\
			"subq	$32,%%rsi		\n\t"\
			"jnz loop_bitclear_short	\n\t"/* loop end; continue is via jump-back if rdi != 0 */\
			:	: [__psmall] "m" (psmall)	/* No outputs; All inputs from memory addresses here */\
				 ,[__startval] "m" (startval)	\
				 ,[__bit_map] "m" (bit_map2)	\
				 ,[__bit_len] "m" (bit_len)	\
				 ,[__m] "m" (m)	\
				 ,[__nprime] "nprime" (nprime-m)	\
				: "cc","memory","cl","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm30","xmm31"	/* Clobbered registers */\
			);
		/*	}	*/

		#else	/******** Non-SIMD (pre-AVX512) **********/

		  #ifdef USE_NCQ
			#warning Using 4-way bit-clear in PerPass_tfSieve.
			uint32 cq[4], ncq = 0;	// mnemonic: ncq = 'Number in to-be-Cleared Queue', cq[] stores said queue
		  #endif
			curr_p = p_last_small;
			for(m = nclear; m < nprime; m++)
			{
				curr_p += (pdiff[m] << 1);	//	if(current_prime == 107) printf("   prime %8d has offset = %8d\n", curr_p, startval[m]);	//	if(pass == 4 && startval[m] < 100) printf("1: Found a multiple of %u in bit %u\n", curr_p,startval[m]);
				l = startval[m];
				// Special-handling code for p == curr_p case:
			//	if(l == 0xffffffff) {	//	printf("hit p == curr_p case!\n");
			//		continue;
			//	}
				// Need to account for the fact that primes greater than the # of bits in the sieve
				// may not hit *any* sieve bit on the current pass through the sieve:
				if(l >= bit_len) {
					startval[m] -= bit_len & -(l != 0xffffffff);	// Try combining regular and (p == curr_p) cases
				} else {
					while(l < bit_len) {
					#ifdef USE_NCQ
						cq[ncq++] = l;
						if(ncq == 4) {
							bit_clr32_x4((uint32 *)bit_map2,cq[0],cq[1],cq[2],cq[3]);
							ncq = 0;
						}
					#else
						bit_clr32((uint32 *)bit_map2,l);
					#endif
					#if DBG_SIEVE
						if(k_targ && l == (i64_targ*64 + bit_targ)) {
							fprintf(stderr,"Critical bit being cleared by prime %u, with offset %u\n", curr_p, startval[m]);	ASSERT(HERE, 0,"0");
						}
					#endif
						l += curr_p;
					}
					/*...save new startvalue:	*/
				#if DBG_SIEVE
					ASSERT(HERE, (startval[m] + startval_incr[m]) < (curr_p + curr_p), "factor.c : (startval[m] + startval_incr[m]) < (curr_p + curr_p)");
					ASSERT(HERE, l-bit_len == (startval[m] + startval_incr[m])%curr_p, "factor.c : l-bit_len == (startval[m] + startval_incr[m])%curr_p");
				#endif
					startval[m] = l-bit_len;
				}
			}

		#endif	// USE_AVX512 ?

//	if(pass==4)printf("\nPass %u: word0 after deep-prime clearing = %16llX\n",pass,bit_map2[0]);

			// Now run through the bits of the current copy of the sieve, trial dividing if a bit = 1:
		  #if TF_CLASSES == 60
			ihi = sieve_len/TF_CLASSES + 1;	// sieve_len not divisible by TF_CLASSES, so add a pad-word
		  #else
			ihi = (sieve_len*64)/TF_CLASSES;	// 64*sieve_len divisible by TF_CLASSES, no need for padding
		  #endif
			ASSERT(HERE, ihi == ((bit_len+63)>>6), "Ihi value-check failed!");
		#ifdef FAC_DEBUG
			m = 0;	// accum popc
			for(i = 0; i < ihi; i++) {
				m += popcount64(bit_map2[i]);
			  #ifdef MULTITHREAD
			//	if(tid == 0)
			  #endif
			//	{ ui64_bitstr(bit_map2[i], cbuf);	printf("%4u: %s\n",i,cbuf); }
			}
		  #ifdef MULTITHREAD
			if(tid == 0)
		  #endif
			printf("%u [%6.2f%%] survived; count = %llu\n",m,100.*(float)m/bit_len,count);
		#endif

			bit_hi = 64;
			for(i = 0; i < ihi; i++)	/* K loops over 64-bit registers. Don't assume bit_len a multiple of 64.	*/
			{
				// Special casing for last sievelet word, which may be only partly full:
				if((bit_len - (i<<6)) < 64) {
					bit_hi = (bit_len - (i<<6));
				}
				for(bit = 0; bit < bit_hi; bit++)	// BIT loops over bits in each sievelet word
				{
				#ifdef FAC_DEBUG
					/* If a known factor is specified, here it is in the bitmap: */
					if(ABS((int64)(k-k_targ)) < 1000) printf("Trying k = %llu\n",k);
					if(k == k_targ) {
						printf("here it is: sweep = %s, bitmap word = %u, bit = %3u\n", &cbuf[convert_uint64_base10_char(cbuf, sweep)], i, bit);
						if((bit_map2[i] >> bit) & 1)
							printf("Trying k_targ = %llu...\n", k_targ);
						else
							ASSERT(HERE, 0,"0");
					}
				#endif

				// *** If current sieve bit=1, add q to trial-division queue: ***

					if((bit_map2[i] >> bit) & 1)
					{
						q_index = count++ & TRYQM1;	/* Post-increment count, so this will work. */

						/* Every so often (every [2^32 / (nearest power of 2 to lenP*lenQ^2)] q's seems a good interval)
						do some factor-candidate sanity checks.

						***NOTE: *** The (count & ...) here must be with a constant = 2^n-1, n >= 3.

						Due to the thread-unsafeness of the mi64 library used herein, in || mode, serialize execution here with a mutex.
						*/
						if((count & countmask) == 0)
						{
							fprintf(stderr,"[k = %llu]",k);
						#ifdef MULTITHREAD
							pthread_mutex_lock(&mutex_mi64);
						//	printf("Count = %u * 2^%u checkpoint: Thread %u locked mutex_mi64 ... ",(uint32)(count >> CMASKBITS),CMASKBITS,tid);
						#endif
							fp = mlucas_fopen(OFILE,"a");
							ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,q,lenQ), "2.k.p overflows!");
							q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
						#ifdef FAC_DEBUG
							sprintf(cbuf, " Count = %u * 2^%u: k = %llu, Current q = %s\n",
								(uint32)(count >> CMASKBITS),CMASKBITS,k,&cbuf[convert_mi64_base10_char(cbuf, q, lenQ, 0)]);
							fprintf(stderr, "%s", cbuf);
						#endif

							/* Do some sanity checks on the occasional factor candidate to ensure
							that the sieve is functioning properly and no overflows have occurred.
							*/

							/* In MMp mode, make sure results of fast and slow-debug modpow agree: */
							if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS) {
								res = (mi64_twopmodq_qmmp(findex, k, u64_arr) == 1);
								if(res != (mi64_twopmodq(p, lenP, k, q, lenQ, q2) == 1) || q2[0] != u64_arr[0]) {
									sprintf(cbuf, "ERROR: Spot-check k = %llu, Results of mi64_twopmodq_qmmp and mi64_twopmodq differ!\n", k);
									fprintf(fp,"%s", cbuf);
									ASSERT(HERE, 0, cbuf);
								}
							}

							/* Make sure that q == 1 (mod 2p) and that q is a base-2 PRP: */
							mi64_clear(u64_arr, lenQ);	// Use q2 for quotient [i.e. factor-candidate k] and u64_arr for remainder
							mi64_div(q,two_p,lenQ,lenQ,q2,u64_arr);
							if(mi64_getlen(q2, lenQ) != 1) {
								sprintf(cbuf, "ERROR: Count = %u * 2^%u: k = %llu, Current q = %s: k must be 64-bit!\n",
									(uint32)(count >> CMASKBITS),CMASKBITS,k,&cbuf[convert_mi64_base10_char(cbuf, q, lenQ, 0)]);
								fprintf(fp,"%s", cbuf);
								ASSERT(HERE, 0, cbuf);
							}
							if(!mi64_cmp_eq_scalar(u64_arr, 1ull, lenQ))
							{
								sprintf(cbuf, "ERROR: Count = %u * 2^%u: k = %llu, Current q = %s: q mod (2p) = %s != 1!\n",
									(uint32)(count >> CMASKBITS),CMASKBITS,k,&cbuf[convert_mi64_base10_char(cbuf, q, lenQ, 0)],
									&cbuf2[convert_mi64_base10_char(cbuf2, u64_arr, lenQ, 0)]);
								fprintf(fp,"%s", cbuf);
								ASSERT(HERE, 0, cbuf);
							}

							/* If q is composite [only check this in debug mode since it costs more than checking
							divisibility of 2^p-1 by q], make sure its smallest divisor
							is larger than the max sieving prime being used for the run: */
							mi64_set_eq    (q2, q, lenQ);
							mi64_sub_scalar(q2,1ull,q2,lenQ);	// Re-use q2 to store q-1
							if(mi64_twopmodq(q2, lenQ, 0, q, lenQ, 0x0) != 1) {
							#if SPOT_CHECK
								printf(" INFO: Spot-check q with k = %llu is composite\n",k);
							#endif
								l = 3;
								for(m = 0; m < nprime; m++) {
									l += (pdiff[m] << 1);
									// Is q % (current small sieving prime) == 0?
									// Cast-to-32-bit-array means doubling the length argument, but THAT IS DONE AUTOMATICALLY INSIDE THE FUNCTION
									if(mi64_is_div_by_scalar32((uint32 *)q, l, lenQ)) {
									#ifdef MULTITHREAD
									//	if(tid != 0) break;	// Can make thread-specific by fiddling the rhs of the !=
										printf("Thread %u, k = %llu: q = ",tid,k);
										if(lenQ > 1)printf("2^64 * %llu + ",q[1]);
										printf("%llu has a small divisor: %u\n",q[0], l);
										ASSERT(HERE, 0, "Abort...");
									#else
										sprintf(cbuf, "ERROR: Count = %u * 2^%u: k = %llu, Current q = %s has a small divisor: %u\n",
											(uint32)(count >> CMASKBITS),CMASKBITS,k,&cbuf[convert_mi64_base10_char(cbuf, q, lenQ, 0)],l);
										fprintf(fp,"%s", cbuf);
										ASSERT(HERE, 0, cbuf);
									#endif
									}
								}
							} else {
							#if SPOT_CHECK
								printf(" INFO: Spot-check q with k = %llu is base-2 PRP\n",k);
							#endif
							}
							fclose(fp); fp = 0x0;

						#ifdef MULTITHREAD
						//	printf("Thread %u unlocking mutex_mi64.\n",tid);
							pthread_mutex_unlock(&mutex_mi64);
						#endif
						}	/* endif((count & countmask) == 0) */

					/***************************************************************************************/
					#if(TRYQ == 0)	/* If testing speed of sieve alone, skip to incrementing of q. */
					/***************************************************************************************/

					#else

						k_to_try[q_index] = k;

						if(q_index == TRYQM1)
						{
							q_index = -1;	/* q_index = -1 indicates factor-candidate queue empty. (More precisely,
											it will be, after we test the current batch of candidates. */

					/***************************************************************************************/
						#if(TRYQ == 1)		/************** try 1 factor candidates at a time **************/
					/***************************************************************************************/

						  #ifdef NWORD

							if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS) {
								if(findex > 10000) {
									if(k < 1000) {
										printf("Do deep sieving for k = %u\n",(uint32)k);
									/****** Apr 2105: This all needs to be made thread-safe ******/
									ASSERT(HERE, 0, "This all needs to be made thread-safe!");
										kdeep[*ndeep++] = (uint32)k;
										ASSERT(HERE, *ndeep < 1024, "Increase allocation of kdeep[] array or use deeper sieving bound to reduce #candidate k's!");
									//	itmp64 = factor_qmmp_sieve64((uint32)findex, k, MAX_SIEVING_PRIME+2, 0x0001000000000000ull);
									//	if(itmp64) {
									//		printf("Q( k = %u ) has a small factor: %20llu\n",(uint32)k, itmp64);
									//	}
									}
									res = 0;
								} else {
									ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,q,lenQ), "2.k.p overflows!");
									q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
									res = (mi64_twopmodq_qmmp(findex, k, u64_arr) == 1);
								// Uncomment to debug by comparing the results of the slow and fast-MMp-optimized modmul routines
								/*
									if(res != (mi64_twopmodq(p, lenP, k, q, lenQ, q2) == 1) || q2[0] != u64_arr[0]) {
										ASSERT(HERE, 0, "bzzt!");
									}
								*/
								}
							} else {
								ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,q,lenQ), "2.k.p overflows!");
								q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
								res = mi64_twopmodq(p, lenP, k, q, lenQ, u64_arr);
							}

						  #elif(defined(P4WORD))

							ASSERT(HERE, 0ull == mi64_mul_scalar(two_p,k,(uint64*)&q256,lenQ), "2.k.p overflows!");
							q256.d0 += 1;	// No need to check for carry since 2.k.p even
							p256.d0 = p[0]; p256.d1 = p[1]; p256.d2 = p[2]; p256.d3 = p[3];
							t256 = twopmodq256(p256,q256);
							res = CMPEQ256(t256, ONE256);

						  #elif(defined(P3WORD))

						   #ifdef USE_FLOAT

							ASSERT(HERE, !p[2], "twopmodq200: p[2] nonzero!");
							x256 = twopmodq200_8WORD_qmmp(p,k);	res = (uint64)CMPEQ256(x256, ONE256);

						   #else

							ASSERT(HERE, 0ull == mi64_mul_scalar(two_p,k,(uint64*)&q192,lenQ), "2.k.p overflows!");
							q192.d0 += 1;	// No need to check for carry since 2.k.p even
							p192.d0 = p[0]; p192.d1 = p[1]; p192.d2 = p[2];
							t192 = twopmodq192(p192,q192);
							res = CMPEQ192(t192, ONE192);

						   #endif	/* USE_FLOAT */

						  #elif(defined(P2WORD))

						   #if USE_128x96 == 1
							/* Use strictly  96-bit routines: */
							if(p[1] == 0 && (q[1] >> 32) == 0)
								res = twopmodq96	(p[0],k);
							else
						   #elif USE_128x96 == 2
							/* Use hybrid 128_96-bit routines: */
							if(p[1] == 0 && (q[1] >> 32) == 0)
								res = twopmodq128_96(p[0],k);
							else
						   #endif
							/* Use fully 128-bit routines: */
							res = twopmodq128x2(p[0],k);

						  #else

							#ifdef USE_FMADD
								/* Use 50x50-bit FMADD-based modmul routines, if def'd: */
								res = twopmodq100_2WORD_DOUBLE(p[0],k);
							#elif(defined(USE_FLOAT))
								/* Otherwise use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE(p[0],k);
							#else
								/* Otherwise use pure-integer-based modmul: */
								if(fbits_in_q < 63) {
									itmp64 = k*(p[0]<<1) + 1;
									res = twopmodq63(p[0],itmp64);
									res = (res == 1);
								} else if(fbits_in_q < 64) {
									itmp64 = k*(p[0]<<1) + 1;
									res = twopmodq64(p[0],itmp64);
									res = (res == 1);
								}
							  #ifdef USE_65BIT
								else if(fbits_in_q < 65)
									res = twopmodq65(p[0],k);
							  #endif
								else
								{
									ASSERT(HERE, fbits_in_q < 96, "fbits_in_q exceeds allowable limit of 96!");
								  #if USE_128x96 == 1
									/* Use strictly  96-bit routines: */
									res = twopmodq96	(p[0],k);
								  #elif USE_128x96 == 2
									/* Use hybrid 128_96-bit routines: */
									res = twopmodq128_96(p[0],k);
								  #else
									/* Use fully 128-bit routines: */
								//	ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,(uint64*)&q128,lenQ), "2.k.p overflows!");
									res = twopmodq128x2(p,k);
								  #endif
								}
							#endif	/* #ifdef USE_FMADD */
						  #endif	/* #ifdef NWORD */

					/***************************************************************************************/
						#elif(TRYQ == 2)	/************** try 2 factor candidates at a time **************/
					/***************************************************************************************/

						  #ifdef P3WORD

							#ifndef USE_FLOAT
								#error	TRYQ = 2 / P3WORD only allowed if USE_FLOAT is defined!
							#endif	/* #ifdef USE_FMADD */

							ASSERT(HERE, !p[2], "twopmodq200: p[2] nonzero!");
							res  = twopmodq200_8WORD_qmmp_x2_sse2(p,k_to_try[0],k_to_try[1]);

						  #elif(defined(P1WORD))

							#ifdef USE_FMADD
								/* Use 50x50-bit FMADD-based modmul routines, if def'd: */
								res = twopmodq100_2WORD_DOUBLE_q2(p[0],k_to_try[0],k_to_try[1]);
							#elif(defined(USE_FLOAT))
								/* Otherwise use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE_q2(p[0],k_to_try[0],k_to_try[1], 0,tid);
							#else
								#error	TRYQ = 2 / P1WORD only allowed if USE_FLOAT or USE_FMADD is defined!
							#endif	/* #ifdef USE_FMADD */

						  #else

							#error TRYQ == 2 requires P1WORD or P3WORD to be defined!

						  #endif

					/***************************************************************************************/
						#elif(TRYQ == 4)	/************** try 4 factor candidates at a time **************/
					/***************************************************************************************/

						  #ifdef P3WORD

						//	ASSERT(HERE, !p[2], "twopmodq200: p[2] nonzero!");
							res = twopmodq192_q4(p,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

						  #elif(defined(P2WORD))

							/* Itanium seems to perform poorly with 96-bit variant: */
						   #if USE_128x96 == 1
							/* Use strictly  96-bit routines: */
							if(p[1] == 0 && fbits_in_q < 96)
								res = twopmodq96_q4     (p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3], 0,tid);
							else
						   #elif USE_128x96 == 2
							/* Use hybrid 128_96-bit routines: */
							if(p[1] == 0 && fbits_in_q < 96)
								res = twopmodq128_96_q4 (p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
							else
						   #endif
							/* Use fully 128-bit routines: */
								res = twopmodq128_q4    (p   ,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

						  #else	// Default single-word-p mode:

							#ifdef USE_FMADD
								/* Use 50x50-bit FMADD-based modmul routines, if def'd: */
								res = twopmodq100_2WORD_DOUBLE_q4(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

							#elif(defined(USE_FLOAT))
								/* Otherwise use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE_q4(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3], 0,tid);

							#else
								/* fbits_in_q is calculated for largest q of current
								set, so if that > 64, use 96-bit routines for all q's. */
							  #ifndef YES_ASM
								if(fbits_in_q < 63)
									res = twopmodq63_q4(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								else
							  #endif	// If x86_64, use faster inline-asm-ified 64-bit modpow for all q's < 2^64:
								if(fbits_in_q < 64)
									res = twopmodq64_q4(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

							  #ifdef USE_65BIT
								else if(fbits_in_q < 65)
									res = twopmodq65_q4(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

							  #endif
							  #ifdef USE_72BIT
								else if(fbits_in_q < 72)
									res = twopmodq72_q4(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
							  #endif
								else {
									ASSERT(HERE, fbits_in_q < 96, "fbits_in_q exceeds allowable limit of 96!");
								#if USE_128x96 == 1
									/* Use strictly  96-bit routines: */
									res = twopmodq96_q4		(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3], 0,tid);
								#elif USE_128x96 == 2
									/* Use hybrid 128_96-bit routines: */
									res = twopmodq128_96_q4	(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								#else
									/* Use fully 128-bit routines: */
									res = twopmodq128_q4	(p   ,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								#endif
								}
							#endif	/* #ifdef USE_FLOAT */
						  #endif	/* #ifdef NWORD */

					/***************************************************************************************/
						#elif(TRYQ == 8)	/************** try 8 factor candidates at a time **************/
					/***************************************************************************************/

						  #ifdef P3WORD
						/*
							if((q[2] >> 32) == 0)
								res = twopmodq160_q8(p,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
							else
						*/
								res = twopmodq192_q8(p,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);

						  #elif(defined(P2WORD))

							/* Itanium seems to perform poorly with 96-bit variant: */
						   #if USE_128x96 == 1
							/* Use strictly  96-bit routines: */
							if(p[1] == 0 && fbits_in_q < 96)
								res = twopmodq96_q8     (p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7], 0,tid);
							else
						   #elif USE_128x96 == 2
							/* Use hybrid 128_96-bit routines: */
							if(p[1] == 0 && fbits_in_q < 96)
								res = twopmodq128_96_q8 (p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
							else
						   #endif
							/* Use fully 128-bit routines: */
								res = twopmodq128_q8    (p   ,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);

						  #else

							/* fbits_in_q is calculated for largest q of current
							set, so if that > 64, use 65,78 or 96-bit routines for all q's. */
							#if defined(USE_FLOAT) && defined(USE_SSE2) && (OS_BITS == 64)
								/* Otherwise use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE_q8(p[0],k_to_try, 0,tid);
							#else
								if(fbits_in_q < 63)
									res = twopmodq63_q8(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								else if(fbits_in_q < 64)
									res = twopmodq64_q8(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
							  #ifdef USE_65BIT
								else if(fbits_in_q < 65)
									res = twopmodq65_q8(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
							  #endif
								else {
									ASSERT(HERE, fbits_in_q < 96, "fbits_in_q exceeds allowable limit of 96!");
								#if USE_128x96 == 1
									/* Use strictly  96-bit routines: */
									res = twopmodq96_q8		(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7], 0,tid);
								#elif USE_128x96 == 2
									/* Use hybrid 128_96-bit routines: */
									res = twopmodq128_96_q8	(p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								#else
									/* Use fully 128-bit routines: */
									p128.d0 = p[0]; p128.d1 = 0;
									res = twopmodq128_q8	(p   ,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								#endif
								}

							#endif	/* #ifdef USE_FLOAT */
						  #endif	/* #ifdef NWORD */

					/***************************************************************************************/
						#elif(TRYQ == 16)	/************** try 16 factor candidates at a time *************/
					/***************************************************************************************/

							#if(defined(NWORD) || defined(P4WORD) || defined(P3WORD) || defined(P2WORD))
								#error (TRYQ == 16) only supported for 64-bit/P1WORD/GCC/AVX builds!
							#elif defined(USE_FLOAT) && defined(USE_AVX) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
								/* Use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE_q16( p[0], k_to_try, 0,tid );
							#else
								#error (TRYQ == 16) only supported for 64-bit/P1WORD/GCC/AVX builds!
							#endif	/* #ifdef USE_FLOAT */

					/***************************************************************************************/
						#elif(TRYQ == 32)	/************** try 32 factor candidates at a time *************/
					/***************************************************************************************/

							#if(defined(NWORD) || defined(P4WORD) || defined(P3WORD) || defined(P2WORD))
								#error (TRYQ == 32) only supported for 64-bit/P1WORD/GCC/AVX512 builds!
							#elif defined(USE_FLOAT) && defined(USE_AVX512) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
								/* Use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE_q32( p[0], k_to_try, 0,tid );
							#else
								#error (TRYQ == 32) only supported for 64-bit/P1WORD/GCC/AVX512 builds!
							#endif	/* #ifdef USE_FLOAT */

					/***************************************************************************************/
						#elif(TRYQ == 64)	/************** try 64 factor candidates at a time *************/
					/***************************************************************************************/

							#if(defined(NWORD) || defined(P4WORD) || defined(P3WORD) || defined(P2WORD))
								#error (TRYQ == 64) only supported for 64-bit/P1WORD/GCC/AVX512 builds!
							#elif defined(USE_FLOAT) && defined(USE_AVX512) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
								/* Use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE_q64( p[0], k_to_try, 0,tid );
							#else
								#error (TRYQ == 64) only supported for 64-bit/P1WORD/GCC/AVX512 builds!
							#endif	/* #ifdef USE_FLOAT */

						#endif	/* endif(TRYQ == ...) */

							/* Print any factors that were found in the current batch: */
							for(l = 0; l < TRYQ; l++)
							{
								if((res >> l) & 1)	/* If Lth bit = 1, Lth candidate of the inputs is a factor */
								{
								#ifdef MULTITHREAD
									pthread_mutex_lock(&mutex_mi64);
								//	printf("Found Factor: Thread %u locked mutex_mi64 ... ",tid);
								#endif
									/* Recover the factor: */
									q[lenP] = mi64_mul_scalar( p, 2*k_to_try[l], q, lenP);
									q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
									if(mi64_twopmodq(p, lenP, k_to_try[l], q, lenQ, q2) != 1)
									{
										fprintf(stderr, "ERROR: k = %llu, post-check indicates this does not yield a factor.\n", k_to_try[l]);
									//	printf("Args sent to mi64_twopmodq:\n");
									//	printf("p = %s\n", &cbuf[convert_mi64_base10_char(cbuf, p, lenP, 0)]);
									//	printf("q = %s\n", &cbuf[convert_mi64_base10_char(cbuf, q, lenQ, 0)]);
									//	printf("res = %s\n", &cbuf[convert_mi64_base10_char(cbuf, q2, lenQ, 0)]);
									} else {
										/* Do a quick base-3 compositeness check (base-2 would be much faster due to
										our fast Montgomery arithmetic-based powering for that, but it's useless for
										weeding out composite Mersenne factors since those are all base-2 Fermat pseudoprimes).
										If it's composite we skip it, since we expect to recover the individual prime subfactors
										on subsequent passes (although this should only ever happen for small p and q > (2p+1)^2 :
										*/
										uint32 known_factor_div_check_done = 0;
									TEST_FAC_PRIM:
										if(mi64_pprimeF(q, 3ull, lenQ)) {
											factor_k[(*nfactor)++] = k_to_try[l];
											if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
												sprintf(cbuf,"\n\tFactor found: q = %s = 2^(%u+2)*%llu. This factor is a probable prime.\n",&cstr[convert_mi64_base10_char(cstr, q, lenQ, 0)],findex,k_to_try[l]/2);
											else
												sprintf(cbuf,"\n\tFactor found: q = %s = 2*p*k + 1 with k = %llu. This factor is a probable prime.\n",&cstr[convert_mi64_base10_char(cstr, q, lenQ, 0)],k_to_try[l]);
										#ifdef FAC_DEBUG
											if(TRYQM1 > 1)
												printf("factor was number %u of 0-%u in current batch.\n", l, TRYQM1);
										#endif
										} else {	// Composite factor; this should only occur in "single-word" (q < 2^96) mode:
											if(known_factor_div_check_done) {	// Already divided out all pvsly-found factors
												sprintf(cbuf,"\n\tComposite Factor found: q = %s; you will have to factor this one separately.\n",&cstr[convert_mi64_base10_char(cstr, q, lenQ, 0)]);
											} else {
												printf("\n\tComposite Factor found: q = %s; checking if any previously-found ones divide it...\n",&cstr[convert_mi64_base10_char(cstr, q, lenQ, 0)]);
												for(j = 0; j < *nfactor; j++) {
													q2[lenP] = mi64_mul_scalar( p, 2*factor_k[j], q2, lenP);
													ASSERT(HERE, lenP == 1 && q2[lenP] == 0ull, "Unexpected carryout in known-factor computation!");
													q2[0] += 1;	// q2 = 2.k.p + 1; No need to check for carry since 2.k.p even
													mi64_clear(u64_arr, lenQ);	// Use u64_arr for quotient; only care if remainder == 0 or not
													if(mi64_div(q,q2,lenQ,lenQ,u64_arr,0x0)) {
														/* in this case, need to update factor_k entry to reflect k of cofactor>
														Given factor q which is product of 2 factors f1 = 2.k1.p+1 and f2 = 2.k2.p+1,
														the first of which has been previously found, we have
														q = f1*f2 = (2.k1.p+1).(2.k2.p+1) = 4.k1.k2.p^2 + 2.(k1+k2).p + 1 = 2.k.p+1,
														so k = 2.k1.k2.p + (k1+k2) = k1 + k2.(2.k1.p + 1) = k1 + f1.k2 .
														Thus if have pvsly found f1 and now find the composite factor q = f1.f2,
														to get k2 from k and k1, use k2 = (k - k1)/f1: */
														factor_k[*nfactor-1] = (factor_k[*nfactor-1] - factor_k[j])/q2[0];
														if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
															sprintf(cbuf,"\n\tFactor divisible by previously-found factor 2^(%u+2)*%llu.\n",findex,factor_k[j]);
														else
															sprintf(cbuf,"\n\tFactor divisible by previously-found factor 2*p*k + 1 with k = %llu.\n",factor_k[j]);
													}
													mi64_set_eq(q, u64_arr, lenQ);
												}
												known_factor_div_check_done = 1;
												// If dividing out any previously-found factors leaves a nontrivial cofactor, send it back to above is-PRP check:
												if(!mi64_cmp_eq_scalar(q, 1ull, lenQ))
													goto TEST_FAC_PRIM;
											}
										}	/* endif(factor a probable prime?) */
									#ifdef FACTOR_STANDALONE
										fprintf(stderr,"%s", cbuf);
									#else
										fp = mlucas_fopen(STATFILE,"a");	ASSERT(HERE, fp != 0x0,"0");
										fprintf(fp,"%s", cbuf);
										fclose(fp); fp = 0x0;
									#endif
										fp = mlucas_fopen(   OFILE,"a");	ASSERT(HERE, fp != 0x0,"0");
										fprintf(fp,"%s", cbuf);
										fclose(fp); fp = 0x0;
									#ifdef QUIT_WHEN_FACTOR_FOUND
										return 0;
									#endif
									}	// end(L-loop)
								#ifdef MULTITHREAD
								//	printf("Thread %u unlocking mutex_mi64.",tid);
									pthread_mutex_unlock(&mutex_mi64);
								#endif
								}	/* endif((res >> l) & 1)		*/
							}	/* endfor(l = 0; l < TRYQ; l++)	*/
						}	/* endif(q_index == TRYQM1)		*/

					#endif	/* endif(TRYQ == 0) */

					}	/* endif((bit_map2[i] >> bit) & 1)	*/
					// Increment k, i.e. increment current q by 2*p*TF_CLASSES:
					k += TF_CLASSES;

				} /* end of BIT loop	*/
			}	/* end of K loop	*/

		#if 0	//(TRYQ > 1)	Aug 2022: This code no longer needed; use 'run -bmin 57 -bmax 64 -m 7962742673' to see why
			#error Aug 2022: This code no longer needed!
			/* Clean up any remaining in queue */
			if(q_index >= 0) {
			//	if(pass == 14 && k_to_try[0] > 16300000 && k_to_try[0] < 16340000)
			//		printf("Cleaning up remaining %u candidates in queue...\n",q_index+1);
				for(l = 0; l <= (uint32)q_index; l++) {
				#ifdef FAC_DEBUG
					itmp64 = mi64_mul_scalar(two_p,k_to_try[l],q,lenQ);
					// Should only happen benignly, for q just above a wordcount boundary due to padding at high end of current sieve interval
				//	if(itmp64)
				//		fprintf(stderr,"2.k.p overflows for k = %llu, result = %llu*2^64 + %llu\n",k_to_try[l],itmp64,q[0]);
					q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
				//	if(k_to_try[0] > 16300000 && k_to_try[0] < 16340000)printf("A: Trying k[%u] = %llu, q = %s\n",l,k_to_try[l],&cbuf[convert_mi64_base10_char(cbuf, q, lenQ, 0)]);
				#endif

				#ifdef P4WORD

					ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k_to_try[l],q,lenQ), "2.k.p overflows!");
					q[0] += 1;	// No need to check for carry since 2.k.p even
					t256 = twopmodq256(*(uint256*)p,*(uint256*)q);
					res = CMPEQ256(t256, ONE256);

				#elif(defined(P3WORD))

					ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k_to_try[l],q,lenQ), "2.k.p overflows!");
					q[0] += 1;	// No need to check for carry since 2.k.p even
					t192 = twopmodq192(*(uint192*)p,*(uint192*)q);
					res = CMPEQ192(t192, ONE192);

				#elif(defined(P2WORD))

					res = twopmodq128x2(p,k_to_try[l]);

				#else

					if(fbits_in_q > 64) {
						res = twopmodq128_96(p[0],k_to_try[l]);
					} else {
						q[0] = 2*p[0]*k_to_try[l] + 1;
						res = twopmodq64(p[0],q[0]);
						res = (res == 1);
					}

				#endif	/* endif(P1WORD) */

					if(res == 1)	/* If Lth bit = 1, Lth candidate of the inputs is a factor */
					{
					#ifdef MULTITHREAD
					//	pthread_mutex_lock(&mutex_mi64);
					#endif
						/* Check if it's a composite factor - if so, skip: */
						if(mi64_pprimeF(q, 3ull, lenQ))
						{
							sprintf(cbuf,"Factor: q = %s. Program: E%s\n",&cbuf[convert_mi64_base10_char(cbuf, q, lenQ, 0)],VERSION);
						#ifdef FACTOR_STANDALONE
							fprintf(stderr,"%s", cbuf);
						#else
							fp = mlucas_fopen(STATFILE,"a");	ASSERT(HERE, fp != 0x0,"0");
							fprintf(fp,"%s", cbuf);
							fclose(fp); fp = 0x0;
						#endif

							fp = mlucas_fopen(   OFILE,"a");	ASSERT(HERE, fp != 0x0,"0");
							fprintf(fp,"%s", cbuf);
							fclose(fp); fp = 0x0;

						#ifdef FAC_DEBUG
							printf("factor was number %u of 0-%u in current batch.\n", l, TRYQM1);
						#endif
							factor_k[(*nfactor)++] = k_to_try[l];

						#ifdef QUIT_WHEN_FACTOR_FOUND
							return 0;
						#endif
						}
					#ifdef MULTITHREAD
					//	pthread_mutex_unlock(&mutex_mi64);
					#endif
					}	/* endif(res == 1) */
				}	/* endfor(l = 0; l <= (uint32)q_index; l++) */
			}	/* endif(q_index >= 0) */
		#endif	/* end #if(TRYQ > 1) */

	/******************* MOVE CHKPT STUFF TO FACTOR-MAIN?? ***********************/
	// Checkpointing only supported for single-threaded runs:
	if(NTHREADS == 1) {
	#if !FAC_DEBUG
		// Every 1024th pass, write the checkpoint file, with format as described previously:
		if(((sweep + 1) %(1024/lenQ + 1)) == 0 || ((sweep + 1) == interval_hi)) {
			i = write_savefile(RESTARTFILE, pstring, pass, k, count);	// Only overwrite passnow, know and count fields of savefile
			ASSERT(HERE,!i,"There were errors writing the savefile ... aborting");
		}	/* Successfully wrote restart file. */
	#endif /* #if !FAC_DEBUG */
	}
#if 0
	#error Multithreaded checkpointing stuff needs debug & test!
	else {
	  #if (!FAC_DEBUG)
		char *char_addr, TMPFILE[STR_MAX_LEN] = "";
		strcpy(TMPFILE, RESTARTFILE);	strcat(TMPFILE, ".tmp");
		/*
		Every 1024th pass through the small-primes sieve, and also following the final pass
		through the sieve, write the checkpoint file, with format as described previously.

		Since we expect that multiple jobs and/or threads may be working on the same exponent,
		we make such checkpoint-updates atomic as follows:

		1. job/thread X acquires file lock and opens checkpoint file <filename> for reading, if it exists.
			No other job/thread may acquire file lock until X releases it;
		2. X also opens a 2nd, temporary, file <filename.tmp> for writing, beginning with line 1: exponent (= pstring)
			and line 2, Value of TF_PASSES in the build (16 or 960)
		3. If there was an existing savefile found in step [1], X copies its contents
			to the .tmp file in [2] line-by-line, only updating the single "Pass *: [max k reached]" entry
			corresponding to the current pass whose progress is being saved via checkpointing;
		4. X closes both files, renames <filename.tmp> to <filename>, thus overwriting the now-obsolete
			version of the latter;
		5. X releases the file lock and resumes processing factor candidates.
		*/
		if(((sweep + 1) %(1024/lenQ + 1)) == 0 || ((sweep + 1) == interval_hi))
		{
			/* TF restart files are in HRF, not binary: */
			fp = mlucas_fopen(RESTARTFILE, "r");
			fq = mlucas_fopen(    TMPFILE, "w");
			if(!fp) {
				fprintf(stderr,"INFO: factoring savefile %s not found - will create.\n",RESTARTFILE);
			} else {	// If file exists, it should have the proper first 2 lines:
				itmp = fscanf(fp,"%s\n",cstr); if(itmp <= 0 || !STREQ(cstr,pstring)) ASSERT(HERE,0,"Line 1 entry found in factoring savefile does not match exponent of run.");
				itmp = fscanf(fp,"%u\n",&i  ); if(itmp <= 0 || i != TF_PASSES      ) ASSERT(HERE,0,"Line 2 entry found in factoring savefile does not match TF_PASSES value of build.");
			}
			if(!fq) {
				fprintf(stderr,"INFO: Unable to open factoring savefile %s for reading and/or %s.tmp for writing...quitting.\n",RESTARTFILE,RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}

			curr_line = 0;
			/* pstring: */
			++curr_line;
			itmp = fprintf(fq,"%s\n",pstring);
			if(itmp <= 0) {
				fprintf(stderr,"ERROR: unable to write Line %d (current exponent) to %s.\n", curr_line, TMPFILE);
				ASSERT(HERE, 0,"0");
			}
			/* TF_PASSES: */
			++curr_line;
			itmp = fprintf(fq,"%u\n",TF_PASSES);
			if(itmp <= 0) {
				fprintf(stderr,"ERROR: unable to write Line %d (TF_PASSES of build) to %s!\n", curr_line, TMPFILE);
				ASSERT(HERE, 0,"0");
			}

			// Now copy any remaining entries in existing file, modifying only the one corr. to the current pass, if it exists:
			if(fp) {
				while(fgets(cstr,STR_MAX_LEN,fp)) {
					if((char_addr = strstr(cstr,"Pass ")) != 0) {
						itmp = sscanf(char_addr,"%u",i);
						if(itmp <= 0) {
							fprintf(stderr,"ERROR: unable to read [Pass *: k] entry: offending line = [%s]\n",cstr); ASSERT(HERE, 0,"0");
						}
						if(i == pass) {	// Is the pass index the one we are updating? If yes, update the k-value
							ASSERT(HERE, !found_pass, "Multiple current-pass entry found in savefile!");
							found_pass = TRUE;
							// Calculate the current k-value
							k = (uint64)incr + (sweep+1)*(sieve_len<<6);
							fprintf(fq,"Pass %u: %llu\n",pass,k);
						} else			// Otherwise just copy as-is
							fputs(cstr,fq);
					} else {	// Just copy as-is
						fputs(cstr,fq);
					}
				}
			}
			if(fp) { fclose(fp); fp = 0x0; }
			fclose(fq); fq = 0x0;
			if(rename(TMPFILE,RESTARTFILE)) {
				sprintf(cstr,"ERROR: unable to rename %s file ==> %s.\n",TMPFILE,RESTARTFILE);
				ASSERT(HERE,0,cstr);
			}
		}	// Successfully updated restart file.
	  #endif /* #if !FAC_DEBUG */
	}	// endif(NTHREADS > 1)
#endif	// #if 0|1

	#if DBG_SIEVE
		QUIT:
	#endif
		#ifdef MULTITHREAD
		  if(tid == 0) {	// In || mode, only the 0-thread accumulated runtime
		#endif
		#ifdef CTIME
			clock2 = clock();
		#else
			clock2 = getRealTime();
		#endif
			*tdiff += (double)(clock2 - clock1);
			clock1 = clock2;
//		exit(0);
		#ifdef MULTITHREAD
		  }
		#endif
			continue;
		} /* end of sweep loop	*/

	  #ifdef MULTITHREAD

		pthread_mutex_lock(&mutex_updatecount);
	//	printf("Thread %u locked mutex_updatecount ... Updating q-tried count: %llu + %llu = ",tid,*(targ->count),count);
		*(targ->count) += count;
	//	printf("%llu ... Thread %u done.\n",*(targ->count),tid);
		pthread_mutex_unlock(&mutex_updatecount);
		return 0x0;
	  #else
		return count;
	  #endif
	}

#endif	// USE_GPU ?

/******************/

/* For an exponent p and a factor index k, each either unmodded or (mod 60), does one of 2 things, depending
on the nullity (or not) of the input array pointer *incr:

[1] incr == 0x0: Checks validity of the input [p,k] (mod 60) combination, by returning the
factoring pass number (in unit-offset form, i.e. 1-16) on which the factor with the given k-mod value
should occur if it's one of the valid combinations of p%60 and k%60. If invalid, returns 0.

[2] incr != 0x0: Checks validity of the input p (mod 60) value, i.e. checks that p can possibly be prime
according to its (mod 60) residue. (We assume the unmodded p > 60 here.) If invalid, returns 0;
otherwise populates the arglist incr[] array with the 16 increments in k (mod 60) covering all the
valid residue classes for the given p (mod 60). These increments sum to 60, i.e. te final pass
will be the k == 0 (mod 60) one.
*/
uint32 CHECK_PKMOD60(uint64*p, uint32 lenP, uint64 k, uint32*incr)
{
	uint32 i,kcur, *iptr = 0x0;	// iptr will point to either the in-array (if one provided) or the following local array:
	uint32 iloc[16];
	uint32 pm = 0, km = k%60, FERMAT = 0;
	uint64 q = 0ull;
	pm = mi64_div_by_scalar64(p,60ull,lenP,0x0);
	 q = 2*km*pm + 1;
	// Note isPow2 works the same for unmodded and modded exponent, e.g. 2^[2,3,4,5,6,...]%60 = [4,8,16,32,4,...]:
	FERMAT = (pm > 1) && isPow2(pm);

	if((pm%3 == 0) || (pm%5 == 0))
		return 0;

	if(incr == 0x0) {	// If no incr-array, check the validity of the km := k (mod 60) value:
		iptr = iloc;
		if(FERMAT) {
			// Fermat: For a valid p-mod, only possible values of km in a factor q = 2.k.p+1 are those for which k even [as shown by Lucas]
			// and for which GCD(2*km*pm + 1, 2*60) = 1, i.e. (2*km*pm + 1) is not divisible by 3 or 5.
			if((km&1) == 0) {
				if((q%3 == 0) || (q%5 == 0))
					return 0;
			} else {
			//	printf("CHECK_PKMOD60: q mod 8 = %u ... invalid.\n",q&7);
				return 0;
			}
		} else {
			// Mersenne: For a valid p-mod, only possible values of km are those for which k == +-1 (mod 8) [by quadratic residuacity]
			// and for which GCD(2*km*pm + 1, 2*60) = 1, i.e. (2*km*pm + 1) is not divisible by 3 or 5.
			if(((q&7) == 1) || ((q&7) == 7)) {
				if((q%3 == 0) || (q%5 == 0))
					return 0;
			} else {
			//	printf("CHECK_PKMOD60: q mod 8 = %u ... invalid.\n",q&7);
				return 0;
			}
		}
	} else {
		iptr = incr;
	}

	// Populate the incr-array (either the arglist one, if provided, or the local one)
	// with the km-increments of valid factor candidates classes corresponding for the given pm-value:
//printf("pm = %u: Acceptable km-values = ",pm);
	i = 0;	// Index of current array slot
	kcur = 0;
	for(k = 1; k <= 60; k++) {	// Unit-offset pass values here!
		q = 2*k*pm + 1;
		if( FERMAT && ( (k&1) || (q%3 == 0) || (q%5 == 0) ) )
			continue;
		else if( ((q&7) == 3) || ((q&7) == 5) || (q%3 == 0) || (q%5 == 0) )
			continue;

		iptr[i++] = k - kcur;	kcur = k;	// iptr stores the *increments* between adjacent km-values
//printf("%u .. ",kcur);
		if(!incr && (k == km)) {	// In Check-pair-valid mode can exit as soon as we have the input-km's pass
			return i;
		}
	}
	ASSERT(HERE, i == 16, "Expect precisely 16 valid k (mod 60) classes!");
//printf("\n");
	return i;	// Nonzero return value indicates success
}

// Same as above, but (mod 4620) - i.e. small-primes 3,5,7,11 built into the sieve length - and 960 resulting passes:
uint32 CHECK_PKMOD4620(uint64*p, uint32 lenP, uint64 k, uint32*incr)
{
	uint32 i,kcur, *iptr = 0x0;	// iptr will point to either the in-array (if one provided) or the following local array:
	uint32 iloc[960];
	uint32 pm = 0, km = k%4620, FERMAT = 0;
	uint64 q = 0ull;
	pm = mi64_div_by_scalar64(p,4620ull,lenP,0x0);
	 q = 2*km*pm + 1;
	// For (mod 4620) we do not have the property that isPow2 works the same for unmodded and modded power-of-2
	// exponents that we do (mod 60), so must infer Fermat-ness based simply on evenness of exponent (modded or not).
	// or better, by call to the mi64 library:
	FERMAT = mi64_isPow2(p,lenP,&i);	// If p a power of 2, binary exponent returned in i

	if((pm%3 == 0) || (pm%5 == 0) || (pm%7 == 0) || (pm%11 == 0))
		return 0;

	if(incr == 0x0) {	// If no incr-array, check the validity of the km := k (mod 60) value:
		iptr = iloc;
		if(FERMAT) {
			// Fermat: For a valid p-mod, only possible values of km in a factor q = 2.k.p+1 are those for which k even [as shown by Lucas]
			// and for which GCD(2*km*pm + 1, 2*4620) = 1, i.e. (2*km*pm + 1) is not divisible by 3,5,7 or 11.
			if((km&1) == 0) {
				if((q%3 == 0) || (q%5 == 0) || (q%7 == 0) || (q%11 == 0))
					return 0;
			} else {
			//	printf("CHECK_PKMOD4620: q mod 8 = %u ... invalid.\n",q&7);
				return 0;
			}
		} else {
			// Mersenne: For a valid p-mod, the only possible value of km are those for which k == +-1 (mod 8) [by quadratic residuacity]
			// and for which GCD(2*km*pm + 1, 2*4620) = 1, i.e. (2*km*pm + 1) is not divisible by 3,5,7 or 11.
		//	printf("CHECK_PKMOD4620: pm,km = %u,%u: q = %llu [mod 8 = %u]\n",pm,km,q,(uint32)q&7);
			if(((q&7) == 1) || ((q&7) == 7)) {
				if((q%3 == 0) || (q%5 == 0) || (q%7 == 0) || (q%11 == 0))
					return 0;
			} else {
			//	printf("CHECK_PKMOD4620: q mod 8 = %u ... invalid.\n",q&7);
				return 0;
			}
		}
	} else {
		iptr = incr;
	}

	// Populate the incr-array (either the arglist one, if provided, or the local one)
	// with the km-increments of valid factor candidates classes corresponding for the given pm-value:
//printf("pm = %u: Acceptable km-values = ",pm);
	i = 0;	// Index of current array slot
	kcur = 0;
	for(k = 1; k <= 4620; k++) {	// Unit-offset pass values here!
		q = 2*k*pm + 1;
		if( FERMAT && ( (k&1) || (q%3 == 0) || (q%5 == 0) || (q%7 == 0) || (q%11 == 0) ) )
			continue;
		else if( ((q&7) == 3) || ((q&7) == 5) || (q%3 == 0) || (q%5 == 0) || (q%7 == 0) || (q%11 == 0) )
			continue;

		iptr[i++] = k - kcur;	kcur = k;	// iptr stores the *increments* between adjacent km-values
//printf("%u .. ",kcur);
		if(!incr && (k == km)) {	// In Check-pair-valid mode can exit as soon as we have the input-km's pass
			return i;
		}
	}
	ASSERT(HERE, i == 960, "Expect precisely 960 valid k (mod 4620) classes!");
	return i;	// Nonzero return value indicates success
}

// Computes 2*p (mod curr_p):
uint32 twop_mod_smallp(const int MODULUS_TYPE, const uint64*two_p, const uint32 findex, const uint32 len2P, const uint32 curr_p)
{
	uint32 r;
#ifdef P1WORD
	r = two_p[0] % curr_p;
#else
	// 26. Sep 2012: This step is horribly slow for larger MMp - Accelerate by using that 2p = 2*Mp for double Mersennes.
	// Ex: For p = 43112609, the binary-powering-based version is ~10000x faster than the long-div-based one:
	if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS) {
		r = twompmodq32(findex, curr_p);			// 2^-p (mod q)
		r = modinv32(r, curr_p)-1;	// 2^+p (mod q) = Mp (mod q)
		// modinv32 returns a signed result, do a conditional add to make >= 0
		r += (-((int32)r < 0)) & curr_p;
		r += r;
		if(r >= curr_p) { r -= curr_p; }
	//	ASSERT(HERE, r == mi64_div_y32(two_p, curr_p, 0x0, len2P), "Fast 2p (mod q) for MMp fails!");
	} else {
		r = mi64_div_y32(two_p, curr_p, 0x0, len2P);
	}
#endif
	return r;
}

void	get_startval(
	const int MODULUS_TYPE,
	const uint64 p,		// Only need LSW of this
	const uint32 findex,// Double-Mersenne and Fermat cases
	const uint64*two_p,	// Here, need the full multiword array (of which use just LSW if P!WORD def'd)
	const uint32 lenQ,	// Manyword case ... note we supply max #words needed to hold factor candidates here, i.e. two_p occupies no more than that
	const uint32 bit_len,
	const uint32 interval_lo, const uint32 incr,
	const uint32 nclear, const uint32 nprime, const uint32 p_last_small,
	const uint8 *pdiff,
	      uint32*startval
)
{
	uint32 i, m, curr_p;
	uint64 dstartval;
	/* startbit k (occurrence of first multiple of prime curr_p in first pass through the relevant sievelet) is defined by
		(offset[curr_ p] + k*prime) (mod TF_CLASSES) = incr(pass)-1, k = 0, ... ,TF_CLASSES-1 .
	*/
	curr_p = p_last_small;
	for(m = nclear; m < nprime; m++)
	{
		curr_p += (pdiff[m] << 1);	ASSERT(HERE, pprimeF(curr_p,2), "Alleged curr_p is Composite!");
		uint32 twop_mod_currp = twop_mod_smallp(MODULUS_TYPE, two_p, findex, lenQ, curr_p);	// This handles both the 1-word and multiword-exponent cases
		// Special-handling code for p == curr_p case - this is needed to prevent 0-input  assertion in the modinv computation below.
		// Dec 2019: Replaced (p == curr_p) with if() clause which also catches curr_p-divides-exponent for composite exponents:
		if(!twop_mod_currp) {
			startval[m] = 0xffffffff;
			continue;
		}
	/*
		Given a pass#, get inc = incr[pass] [3452 here] and seek an index x [604 here] such that

		1 + 2*p*( x*TF_CLASSES + inc ) is divisible by the current small-prime, curr_p.
		In other words, we need to find the smallest nonnegative integer x such that

			2*p*( x*TF_CLASSES + inc ) == -1 (mod curr_p),

		i.e. find x such that

			(2*p*TF_CLASSES)*x == (-1 - 2*p*inc) (mod curr_p).

		Letting A := 2*p*TF_CLASSES (mod curr_p) and B := (-1 - 2*p*inc) (mod curr_p), we need to solve A*x == B (mod curr_p)
		for nonnegative integer x.
		If we first find BI := inverse of the rhs constant B (mod curr_p), we multiply both sides by that to get
		(A*BI)*x == 1 (mod curr_p), and then we can simply use a second eGCD to find x = inverse of (A*BI) (mod curr_p):
	*/
		uint32 A = ((uint64)twop_mod_currp * (uint64)TF_CLASSES) % curr_p;
		uint32 B = ((uint64)twop_mod_currp * (uint64)incr + 1) % curr_p;
		if(B == 0) {	// Must guard against eGCD with identical args below
			i = 0;
		} else {
			// Debug NOTE: PARI has a poorly-findable-for-the-non-PARI-expert modinv functionality:
			// its eGCD is hidden behind the bezout(x,y) function; help for that ('? bezout') gives
			//	bezout(x,y): returns [u,v,d] such that d=gcd(x,y) and u*x+v*y=d.
			// Thus e.g. p = 933551; i=809470; bezout(i,p) ==> [-435932, 377991, 1], thus -435932 == 497619 is the inverse of i (mod p).
			uint32 BI = modinv32(curr_p - B, curr_p);	// For the above example (curr_p = 933551; (curr_p - B) = 809470) gives 4294531364, which is just the correct -435932 aliased to 2^32 - 435932
			uint32 ABI = ((uint64)A * (uint64)(curr_p - BI)) % curr_p;	// A*BI
			i = curr_p - modinv32(ABI, curr_p); i %= curr_p;
		}
		startval[m] = i;

	// 2nd Part only Needed for intervals not starting from the min-k for the pass in question:
	/*
		Calculate and store increment of offset for each sieving prime < bit_len,
		used to quickly find what offset will be after another pass through sievelet.
		For each sieving prime curr_p, we hop through the bit_len bits of the
		sievelet in strides of curr_p, starting at bit = startval[curr_p].
		Letting k := ceil(bit_len/curr_p) (i.e. how many times we hit curr_p
		on a single pass through the sieve, rounded up), the change in offset
		due to a single pass through the sieve is

			d(startval) = k*curr_p - bit_len , which we note can also be written as
						= curr_p - (bitlen % curr_p).

		Example: bit_len = 100 (unrealistic number, but that's not important here)
				 curr_p = 17

		Then:	k = ceil(bit_len/curr_p) = ceil(5.88...) = 6,
		and
			d(startval) = k*curr_p - bit_len		 = 6*17 - 100 = 102 - 100 = 2 , or
						= curr_p - (bitlen % curr_p) = 17 - (100%17) = 17- 15 = 2 .

		(For primes > bitlen, (bitlen % curr_p) = bitlen, so the mod is superfluous.)

		Thus the offset at the beginning of the next pass through the sieve is

			startval' = (startval + d(startval))%curr_p .

		If we want to accomodate arbitrarily large kmin values for the start
		of our sieving runs, we need to calculate how many passes through the
		sieve the given kmin value corresponds to - that is stored in the
		interval_lo:

			interval_lo = floor(kmin/64.0/len) ,

		and thus the offset at the beginning of the initial pass through the sieve is

			startval' = (startval + interval_lo*d(startval)))%curr_p ,

		where we'll probably want to do a mod-curr_p of interval_lo prior to
		the multiply by (k*curr_p - bit_len) to keep intermediates < (curr_p)^2,
		which means < 2^64 is we allow curr_p as large as 32 bits.
	*/
		if(interval_lo != 0) {
			/* bit_len is a uint32, so use i (also a 32-bit) in place of k (64-bit) here: */
			i = ceil(1.0*bit_len/curr_p);
			ASSERT(HERE, i*curr_p - bit_len == curr_p - (bit_len % curr_p), "i*curr_p - bit_len == curr_p - (bit_len % curr_p)");

			/* Now calculate dstartval for the actual current-pass kmin value,
			according to the number of times we'd need to run through the sieve
			(starting with k = 0) to get to kmin: */
			dstartval = (uint64)(i*curr_p - bit_len);
			dstartval = (interval_lo*dstartval) % curr_p;
			dstartval += startval[m];
			if(dstartval >= curr_p)
				startval[m] = dstartval - curr_p;
			else
				startval[m] = dstartval;

		#ifdef FAC_DEBUG
			ASSERT(HERE, startval     [m] < curr_p, "factor.c : startval     [m] < curr_p");
		  #if DBG_SIEVE
			startval_incr[m] = i*curr_p - bit_len;
			ASSERT(HERE, startval_incr[m] < curr_p, "factor.c : startval_incr[m] < curr_p");
		  #endif
		#endif
		}
	}	/* endfor(m = nclear; m < nprime; m++) */
}

uint64 given_b_get_k(double bits, const uint64 two_p[], uint32 len)
{
	int i,l;
	uint64 itmp64, k;
	double fqlo, twop_float;
#ifdef P1WORD
	/* Find FP approximation to 2*p - can't use this for multiword case, because double approximation tp 2*p may overflow: */
	twop_float = (double)two_p[0];
	fqlo = pow(2.0, bits);
	k = (uint64)(fqlo/twop_float);
#else
	// In the multiword case, need to compute a double-precision approximation to 2^bits/(2*p)
	// while avoiding possible overflow of a double-exponent field. (I.e. can't directly compute
	// pow(2.0, bits) because b may exceed __DBL_MAX_EXP__ = 1024).
	i = mi64_extract_lead64(two_p, len, &itmp64);	// i has bitlength of 2*p; itmp64 has leading 64 bits
	// Number of low-order bits we discarded in retaining just the leading 64 ... this needs to be signed in case 2p < 2^64
	l = i-64;
	k = (uint64)(pow(2.0, bits-l)/(double)itmp64);
//	convert_uint64_base2_char(cbuf, itmp64);
//	printf("2*p = %16llX has %u bits, lead64 = %s ==> k = %16llu.\n",itmp64,i,cbuf,k);
#endif
	return k;
}

/* The factoring checkpoint file is assumed to have the format:
	Line 1:		{string containing the current exponent stored in pstring}
	Line 2:		TF_PASSES; supporte values are 16 and 960

	Line 3:		bmin = {Log2(minimum factor to try), in floating double form}
				If > 10^9 its whole-number part is taken as the KMin value instead.
	Line 4:		bmax = {Log2(maximum factor to try), in floating double form}
				If > 10^9 its whole-number part is taken as the KMax value instead.

	Line 5:		KMin = {smallest factor K value to be tried in each pass}
	Line 6:		KNow = { largest factor K value tried so far during current pass}
	Line 7:		KMax = { largest factor K value to be tried in each pass}

	Line 8:		PassMin = {maximum pass for the run (typically TF_PASSES-1, but perhaps not, e.g. for a factoring assignment split over multiple CPUs.}
	Line 9:		PassNow = {current factoring pass}
	Line 10:	PassMax = {maximum pass for the run (typically TF_PASSES-1, but perhaps not, e.g. for a factoring assignment split over multiple CPUs.}

	Line 11:	Number of q's tried so far during the run

	Line 12+:	Any diagnostic info not needed for restarting from interrupt
				(mainly, in standalone mode can use this in place of STATFILE.)
*/
int read_savefile(const char*fname, const char*pstring, double*bmin, double*bmax,
uint64*kmin, uint64*know, uint64*kmax, uint32*passmin, uint32*passnow, uint32*passmax, uint64*count)
{
	int itmp;
	uint32 curr_line = 0, nerr = 0;
	uint64 tf_passes = 0;
	char *char_addr;
	/* TF restart files are in HRF, not binary: */
	fp = mlucas_fopen(fname, "r");
	if(!fp) {
		return -1;
	} else {
		sprintf(cbuf,"Factoring savefile %s found ... reading ...\n",fname);
		fprintf(stderr,"%s",cbuf);
	#ifndef FACTOR_STANDALONE
		fq = mlucas_fopen(STATFILE,"a"); fprintf(fq,"%s",cbuf); fclose(fq); fq = 0x0;
	#endif
		/* Line 1: pstring */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (current exponent) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Strip the expected newline char from in_line: */
		char_addr = strstr(in_line, "\n");
		if(char_addr)
			*char_addr = '\0';
		/* Make sure restart-file and current-run pstring match: */
		if(STRNEQ(in_line, pstring)) {
			++nerr; fprintf(stderr,"ERROR: current exponent %s != Line %d of factoring restart file %s!\n",pstring,curr_line,fname);
		}

		/* Line 2: TF_PASSES */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (TF_PASSES) of factoring restart file %s!\n",curr_line,fname);
		}
		char_addr = strstr(in_line, "tf_passes");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'tf_passes' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			tf_passes = convert_base10_char_uint64(char_addr);
			if(tf_passes != TF_PASSES) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s: TF_PASSES value [%llu] mismatches that of build [%u]!\n",curr_line,fname, tf_passes, (uint32)TF_PASSES);
			}
		}

		/* Line 3: bmin */
		++curr_line;
		fgets(in_line, STR_MAX_LEN, fp);
		char_addr = strstr(in_line, "bmin");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'bmin' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
		}
		itmp = sscanf(char_addr, "%lf",bmin);
		if(itmp != 1) {
			++nerr; fprintf(stderr,"ERROR: unable to parse Line %d (bmin) of factoring restart file %s. Offending input = %s\n",curr_line,fname, in_line);
		}

		/* Line 4: bmax */
		++curr_line;
		fgets(in_line, STR_MAX_LEN, fp);
		char_addr = strstr(in_line, "bmax");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'bmax' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
		}
		itmp = sscanf(char_addr, "%lf",bmax);
		if(itmp != 1) {
			++nerr; fprintf(stderr,"ERROR: unable to parse Line %d (bmax) of factoring restart file %s. Offending input = %s\n",curr_line,fname, in_line);
		}

	/************************************
	LINE PAIRS 5/6 AND 7/8 ARE USED TO DETERMINE WHETHER A PREVIOUS
	FACTORING RUN OF THE SAME EXPONENT COMPLETED OR NOT: If know >= kmax
	and passnow = passmax then the previous run completed, in which case
	we allow a new run to a deeper bound, i.e. reset passnow = passmin
	and run passes passmin through passmax from bounds kmin to kmax.
	*************************************/
		/* Line 5: kmin */
		++curr_line;
		fgets(in_line, STR_MAX_LEN, fp);
		char_addr = strstr(in_line, "kmin");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'kmin' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			*kmin = convert_base10_char_uint64(char_addr);
		}

		/* Line 6: know */
		++curr_line;
		fgets(in_line, STR_MAX_LEN, fp);
		char_addr = strstr(in_line, "know");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'know' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			*know = convert_base10_char_uint64(char_addr);
		}

		/* Line 7: kmax */
		++curr_line;
		fgets(in_line, STR_MAX_LEN, fp);
		char_addr = strstr(in_line, "kmax");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'kmax' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			*kmax = convert_base10_char_uint64(char_addr);
		}

		/* Line 8: passmin */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (PassMin) of factoring restart file %s!\n",curr_line,fname);
		}
		char_addr = strstr(in_line, "passmin");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'passmin' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			*passmin = (uint32)convert_base10_char_uint64(char_addr);
			if(*passmin >= TF_PASSES)	{ ++nerr; fprintf(stderr,"factor.c: Require passmin[%u] < TF_PASSES[%u]",*passmin,TF_PASSES); }
		}

		/* Line 9: passnow */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (PassMin) of factoring restart file %s!\n",curr_line,fname);
		}
		char_addr = strstr(in_line, "passnow");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'passnow' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			*passnow = (uint32)convert_base10_char_uint64(char_addr);
			if(*passnow >= TF_PASSES)	{ ++nerr; fprintf(stderr,"factor.c: Require passnow[%u] < TF_PASSES[%u]",*passnow,TF_PASSES); }
			if(*passnow <  *passmin  )	{ ++nerr; fprintf(stderr,"factor.c: Require passnow[%u] >= passmin[%u]",*passnow,*passmin); }
		}

		/* Line 10: passmax */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (PassMin) of factoring restart file %s!\n",curr_line,fname);
		}
		char_addr = strstr(in_line, "passmax");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'passmax' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			*passmax = (uint32)convert_base10_char_uint64(char_addr);
			if(*passmax >= TF_PASSES)	{ ++nerr; fprintf(stderr,"factor.c: Require passmax[%u] < TF_PASSES[%u]",*passmax,TF_PASSES); }
			if(*passmax <  *passnow  )	{ ++nerr; fprintf(stderr,"factor.c: Require passmax[%u] >= passnow[%u]",*passmax,*passnow); }
		}

		/* Line 11: Number of q's tried: */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (#Q tried) of factoring restart file %s!\n",curr_line,fname);
		}
		char_addr = strstr(in_line, "#Q tried");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: '#Q tried' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			*count = convert_base10_char_uint64(char_addr);	// Need to reset == 0 prior to sieving so kvector-fill code works properly
		}
		fclose(fp); fp = 0x0;
	}
	return (int)nerr;
}

// Init savefile with above read_savefile fields so ensuing checkpoint-writes only need to update the pass# and k:
int init_savefile(const char*fname, const char*pstring, double bmin, double bmax,
uint64 kmin, uint64 know, uint64 kmax, uint32 passmin, uint32 passnow, uint32 passmax, uint64 count)
{
	int itmp;
	uint32 curr_line = 0, nerr = 0;
	/* TF restart files are in HRF, not binary: */
	fp = mlucas_fopen(fname,"w");	// Open in write mode
	if(!fp) {
	#ifndef FACTOR_STANDALONE
		fp = mlucas_fopen(STATFILE,"a");
		fprintf(	fp,"INFO: Unable to open factoring savefile %s for writing...quitting.\n",fname);
		fclose(fp); fp = 0x0;
	#endif
		fprintf(stderr,"INFO: Unable to open factoring savefile %s for writing...quitting.\n",fname);
		return -1;
	} else {
		/* Line 1: pstring: */
		++curr_line; itmp = fprintf(fp,"%s\n",pstring);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (current exponent) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 2: TF_PASSES: */
		++curr_line; itmp = fprintf(fp,"tf_passes = %u\n",TF_PASSES);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (TF_PASSES of build) to %s!\n",curr_line,fname);
		}
		/* Line 3: bmin: */
		++curr_line; itmp = fprintf(fp, "bmin = %lf\n", bmin);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (bmin) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 4: bmax: */
		++curr_line; itmp = fprintf(fp, "bmax = %lf\n", bmax);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (bmax) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 5: kmin: */
		++curr_line; itmp = fprintf(fp,"kmin = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, kmin)]);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (kmin) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 6: know: */
		++curr_line; itmp = fprintf(fp,"know = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, know)]);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (know) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 7: kmax: */
		++curr_line; itmp = fprintf(fp,"kmax = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, kmax)]);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (kmax) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 8: passmin: */
		++curr_line; itmp = fprintf(fp,"passmin = %u\n", passmin);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (passmin) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 9: passnow: */
		++curr_line; itmp = fprintf(fp,"passnow = %u\n", passnow);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (passnow) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 10: passmax: */
		++curr_line; itmp = fprintf(fp,"passmax = %u\n", passmax);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (passmax) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Line 11: Number of q's tried: */
		++curr_line; itmp = fprintf(fp,"#Q tried = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, count)]);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (#Q tried) of factoring restart file %s!\n",curr_line,fname);
		}
		fclose(fp); fp = 0x0;
		return (int)nerr;
	}
}

// Only overwrite passnow, know and count fields of savefile:
int write_savefile(const char*fname, const char*pstring, uint32 passnow, uint64 know, uint64 count)
{
	 int itmp;
	uint32 curr_line = 0, nerr = 0, passnow_file;
	uint64 know_file, count_file;
	char *char_addr;
	/* TF restart files are in HRF, not binary: */
	fp = mlucas_fopen(fname,"r+");	// Open in update ("read plus") mode
	if(!fp) {
	#ifndef FACTOR_STANDALONE
		fp = mlucas_fopen(STATFILE,"a");
		fprintf(	fp,"INFO: Unable to open factoring savefile %s for writing...quitting.\n",fname);
		fclose(fp); fp = 0x0;
	#endif
		fprintf(stderr,"INFO: Unable to open factoring savefile %s for writing...quitting.\n",fname);
		return -1;
	} else {
		/* Line 1: pstring */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (current exponent) of factoring restart file %s!\n",curr_line,fname);
		}
		/* Strip the expected newline char from in_line: */
		char_addr = strstr(in_line, "\n");
		if(char_addr)
			*char_addr = '\0';
		/* Make sure restart-file and current-run pstring match: */
		if(STRNEQ(in_line, pstring)) {
			++nerr; fprintf(stderr,"ERROR: current exponent %s != Line %d of factoring restart file %s!\n",pstring,curr_line,fname);
		}

		/* Line 6: know */
		while(++curr_line < 6) {
			if(!fgets(in_line, STR_MAX_LEN, fp)) {
				++nerr; fprintf(stderr,"ERROR: unable to read Line %d of factoring restart file %s!\n",curr_line,fname);
			}
		}
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (know) of factoring restart file %s!\n",curr_line,fname);
		}
		char_addr = strstr(in_line, "know");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'know' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			know_file = convert_base10_char_uint64(char_addr);
		}
		itmp = fprintf(fp,"know = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, know)]);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (know) of factoring restart file %s!\n",curr_line,fname);
		}

		/* Line 7: kmax: */
		++curr_line; fgets(in_line, STR_MAX_LEN, fp);
		/* Line 8: passmin: */
		++curr_line; fgets(in_line, STR_MAX_LEN, fp);

		/* Line 9: passnow: */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (passnow) of factoring restart file %s!\n",curr_line,fname);
		}
		char_addr = strstr(in_line, "passnow");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: 'passnow' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			passnow_file = convert_base10_char_uint64(char_addr);
		}
		itmp = fprintf(fp,"passnow = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, passnow)]);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (passnow) of factoring restart file %s!\n",curr_line,fname);
		}

		// Check progress: compared to previous checkpoint, passnow should be same and know greater, or passnow should be greater:
		if(passnow == passnow_file && know > know_file) {
			/* No-op */
		} else if(passnow > passnow_file) {
			/* No-op */
		} else {
			++nerr; fprintf(stderr,"ERROR: In factoring restart file %s: compared to previous checkpoint, passnow[%u] should be same as file[%u] and know[%llu] greater than file[%llu], or passnow should be greater!\n",fname,passnow,passnow_file,know,know_file);
		}

		/* Line 10: passmax: */
		++curr_line; fgets(in_line, STR_MAX_LEN, fp);

		/* Line 11: Number of q's tried: */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			++nerr; fprintf(stderr,"ERROR: unable to read Line %d (#Q tried) of factoring restart file %s!\n",curr_line,fname);
		}
		char_addr = strstr(in_line, "#Q tried");
		if(!char_addr) {
			++nerr; fprintf(stderr,"ERROR: '#Q tried' not found in Line %d of factoring restart file %s!\n",curr_line,fname);
		} else {
			char_addr = strstr(in_line, "=");
			if(!char_addr) {
				++nerr; fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n",curr_line,fname);
			}
			char_addr++;
			count_file = convert_base10_char_uint64(char_addr);	// Need to reset == 0 prior to sieving so kvector-fill code works properly
		}
		++curr_line; itmp = fprintf(fp,"#Q tried = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, count)]);
		if(itmp <= 0) {
			++nerr; fprintf(stderr,"ERROR: unable to write Line %d (#Q tried) of factoring restart file %s!\n",curr_line,fname);
		}
		fclose(fp); fp = 0x0;
		return (int)nerr;
	}
}

/* This is actually an auxiliary source file, but give it a .h extension to allow wildcarded project builds of form 'gcc -c *.c' */
#include "factor_test.h"

#undef YES_ASM
