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

/* NOTE: Default should simply include this sourcefile as part of an Mlucas build.
To build the sieve factoring code in standalone mode, see the compile instructions below!
*/

#include "factor.h"
#include "align.h"

#define FAC_DEBUG	0
#if FAC_DEBUG
	#define DBG_SIEVE	0
#endif

#ifdef INCLUDE_PM1
	#include "gcd_lehmer.h"
#endif

#define ALLOW_PRP	0	/* Set nonzero to show Fermat base-2 pseudoprimes during small-prime-gen step
							(and to have the latter run 10-30x more slowly). The name here is perhaps
							a tad confusing: "allow" means "don't weed out via use of prestored tables."
						*/

#if(defined(USE_SSE2) && defined(COMPILER_TYPE_MSVC) && OS_BITS == 32)

	#define FERMAT_PSP2	0	/* Set nonzero to activate the special base-2 Fermat-pseudoprime functionality */
	#if FERMAT_PSP2
		#include "f2psp.h"
	#endif

#endif

/* Tables of Fermat-base-2 pseudoprimes needed by the small-primes sieve code: */
#include "f2psp_3_5.h"

/* Factoring-only globals: */
int restart;
uint64 checksum1 = 0;	/* Sum (mod 2^64) of all q's tried during a predetermined interval */
uint64 checksum2 = 0;	/* Sum (mod 2^64) of all (2^p mod q) result valus for a predetermined interval */

uint64 chk1,chk2,chk3,chk4;
#ifdef FACTOR_STANDALONE
	/* In standalone, need to locally define some Mdata.h/Mlucas.h-declared globals
	which would normally be defined in Mlucas.c. */
	char ESTRING[STR_MAX_LEN];	/* Exponent in string form */
	char PSTRING[STR_MAX_LEN];	/* [exponent | index] of [Mersenne | Fermat ] Number being tested in string form, typically
								 estring concatenated with other descriptors, e.g. strcat("M",estring) for Mersenne case
								*/
	uint32 PMIN;	/* minimum #bits allowed for FFT-based mul */
	uint32 PMAX;	/* maximum #bits allowed depends on max. FFT length allowed
					  and will be determined at runtime, via call to given_N_get_maxP(). */
	char cbuf[STR_MAX_LEN];
	char in_line[STR_MAX_LEN];
	/* Declare a blank STATFILE string to ease program logic: */
	char STATFILE[] = "";
	const char*NUM_PREFIX[MODULUS_TYPE_MAX+1] = {"Unknown modulus type!","M","MM","F"};	// Prefixes corr. to #defines in Mdata.h
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
	#if(defined(NWORD))	/* N-word P, N presumed to be > 1: */
		#if FAC_DEBUG
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	2000000
		#endif
	#elif(defined(P4WORD))	/* 4-word P: */
		#if FAC_DEBUG
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	2000000
		#endif
	#elif(defined(P3WORD))	/* 3-word P: */
		#if FAC_DEBUG
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	2000000
		#endif
	#elif(defined(P2WORD))	/* 2-word P: */
		#if FAC_DEBUG
			#define NUM_SIEVING_PRIME	10000
		#else
			#define NUM_SIEVING_PRIME	1000000
		#endif
	#else	/* 1-word P limit is set by #of bits p can have so 120^2*p < 2^64: */
		#if FAC_DEBUG
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
	,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,0x0};

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
	by roughly three-fourths. Specifically, a number of form q=2*k*p+1
	can only be prime if the value of k mod 60 appears in the appropriate
	p mod 60 row of the following table:

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
GCD(2*km*pm + 1, 120) = 1, i.e. 2*km*pm is not divisible by 3 or 5.

Let's try pm = 1 as an example:

pm := p%60 = 1:

km := k%60 :            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59
2*km*pm+1  :            1  3  5  7  9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81 83 85 87 89 91 93 95 97 99.01.03.05.07.09.11.13.15.17.19
GCD(2*km*pm+1, 120)   : 1  3  5  1  3  1  1 15  1  1  3  1  5  3  1  1  3  5  1  3  1  1  5  1  1  3  1  5  3  1  1  3  5  1  3  1  1  5  1  1  3  1  5  3  1  1  3  5  1  3  1  1 15  1  1  3  1  5  3  1
Acceptable? (* = yes) : *        *     *  *     *  *     *        *  *        *     *  *     *  *     *        *  *        *     *  *     *  *     *        *  *        *     *  *     *  *     *        *

We further use quadratic residues: From the form of N = 2^p - 1, we have that
2^p == 1 (mod N). Multiplying by 2, we see that 2^(p+1) == 2 (mod N). Since p is odd
(primality of p is not crucial here, just oddity), the LHS is an even power of 2,
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

	2009-2011	Added support for x86 SSE2 floating-double SIMD hardware; improved inline ASM for x86 uint64 mode.

	2012		Added support for nVidia CUDA, geared toward floating-double capability of Fermi family of GPUs.

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

	FAC_DEBUG - define = 1 in factor.h to enable debugging diagnostic prints
				and assertions of the form ASSERT.

	DBG_SIEVE - define = 1 to debug the sieve without actually spending time
				doing trial divides. NOTE: Requires FAC_DEBUG = 1.

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

		PIPELINE_MUL192 - when defined, uses pipelined versions of 192-bit MUL macros.

	P4WORD - factor limit q < 2^256, i.e. q needs 4 full 64-bit words of storage.

		PIPELINE_MUL256 - when defined, uses pipelined versions of 256-bit MUL macros.

	NWORD - Arbitrary-length p and q, only restriction is that (as for all other size ranges) kmax < 2^64 .
*/

/*********************************************************************************************************************************/

/* NOTE: Exponents > 64 bits *require* standalone-mode build: */

#ifdef FACTOR_STANDALONE

	int main(int argc, char *argv[])
	{

#else

  #if FAC_DEBUG
	#error FAC_DEBUG only permitted in standalone mode!
  #endif

	int factor(char *pstring, double bmin, double bmax)
	{

#endif

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
	/* Restart file name:
	This is set at runtime, based on either command-line -file flag argument
	or (if -file {} not used) the exponent being processed. */
	char RESTARTFILE[STR_MAX_LEN] = "";

	double bmin = 0.0, bmax = 0.0;	/* store log2 of (min|max) factoring bound */
#endif
	/* Make these large enough to hold the max. p,q supported by the software (currently,
	256 bits) then use as many of the subfields as needed, depending on the values of P,QWORDS.
	*/
	uint32	findex = 0, nbits_in_p = 0, nbits_in_q = 0, lenP = 0, lenQ = 0, bits_in_pq2 = 0;
	const double ILG2 = 1.0/log(2.0);
	double	fbits_in_2p = 0, fbits_in_k = 0, fbits_in_q = 0;
	uint64	*factor_ptmp = 0x0, *factor_k = 0x0;	/* Use this array to store the factor k's (assumed < 2^64) of any factors found. */
	uint64	*p = 0x0, *two_p = 0x0, *p120 = 0x0, *q = 0x0, *q2 = 0x0;

#if(!defined(TRYQ))
	#define TRYQ	1
#endif
#if(TRYQ < 1)
	#undef TRYQ
	#define TRYQ	1
#endif

	uint32 TRYQM1 = TRYQ-1;

// This is for the GPU-style vectorized sieve
#ifdef USE_GPU
	const uint32 kvsz = 0x1 << 20, kvm1 = kvsz-1;
	uint32 kvidx, kv_numq, pshift, jshift, zshift, start_index, lead7;
	uint64 *cpu_kvec = 0x0;	// Vector of factor k's needing testing to see if q = 2.p.k+1 divides 2^p-1; "cpu" prefix means host copy.
	uint64 *gpu_kvec = 0x0;	// GPU device copy of kvec.
	uint8 *cpu_rvec = 0x0;	// output vector (given vector set of input TF candidates q = 2.k.p+1, resulting 2^p mod q, in binary "is factor?" form)
	uint8 *gpu_rvec = 0x0;	// GPU device copy of rvec.
	int threadsPerBlock = 0, blocksPerGrid = 0;
#endif

	uint64 *k_to_try = 0x0;
	int32 q_index;

	int   nargs,itmp;
	/* pdsum stores the sums needed for the base (%30 == 0) candidate of each length-30 interval;
	pdiff stores the diffs/2 of these absolute offsets, circularly shifted to the right one place,
	since to update we need diff[j] = curr_p(current) - curr_p(previous) = (pdsum[j] - pdsum[j-1])/2,
	i.e. pdiff_8[i] = (pdsum_8[i]-pdsum_8[i-1])/2, i computed using mod-8 arithmetic and pdsum_8[] using mod-30 arithmetic:
	*/
	uint32 itmp32, pdiff_8[8] = {2,1,2,1,2,3,1,3}, pdsum_8[8] = { 0, 2, 6, 8,12,18,20,26};
	uint64 itmp64;
	uint32 pmod60,kmod60,qmod8;
	uint64 res,two64modp;
	uint32 bit,bit_hi,curr_p,i,ihi,idx,incr[16],m,ncopies,regs_todo;
	uint32 twop_mod_currp,currp_mod_60,lmod60;
	uint32 l,i64,nfactor,word;
	uint32 nprime = NUM_SIEVING_PRIME;	/* Don't make this const, in case we need to reduce it to satisfy MAX_SIEVING_PRIME < (q_min = 2*p+1) */
#if ALLOW_PRP
	uint32 isprime,num_pseudo,sqrt_p;
#else
	uint32 f2psp_idx = 0;	/* Index to next-expected Fermat base-2 pseudoprime in the precomputed table */
#endif
	/* LEN is the number of 64-bit words in our sieving bitmap. Successive
	bits of the sieve bitmap correpond to successive factor k-values, with bit 0
	corresponding to the smallest possible factor, k = 1 (note the unit-offset!)
	In practice, this full-sized sieve is split into 16 smaller (i.e. hopefully
	cache-sized) "sievelets", each of which contains bits corresponding to
	successive k's from one of the 16 allowable k mod 60 families for the given
	exponent p (more specifically, the given p mod 60 value.) Both full-length
	sieve and 1/16-th length sievelets share the property that each time we run
	through the bits of onesuch we run through 64*LEN worth of k's.
	*/
	const uint32 nclear=6, len=3*5*7*11*13*17, p_last_small = 17;	/* Make sure LEN = product of first NCLEAR odd primes,
									p_last_small = largest odd prime appearing in the product. */
	uint32 prime[] = {3,5,7,11,13,17};	/* Also need the first [nclear] odd primes - since 'const' is treated as a readonly flag
						on *variables*, use [] instead of [nclear] to avoid variable-length array init errors. */
	const uint32 bit_len=64*len/60; 	/* Number of bits/quadwords in each of the 16 mod-60 sievelets	*/
/*   bits cleared of multiples of 3,5,7,11,13, 17 and q mod 8 = +-1 are here:	*/
	uint64 *temp_late = 0x0;		/* Even though we know how large to make this, it's only needed
									for data inits, so we calloc it at runtime and free it later. */
/*	uint64 bit_map[len/60+1],bit_atlas[len/60+1][16];	* bit_maps corresponding to 16 separate K mod 60 cases are here.*/
	uint64 *bit_map, *bit_atlas = 0x0;
	uint32 pass,passmin = 0,passnow = 0,passmax = 15;
	uint64 count = 0,countmask,j,k,kmin = 0,kmax = 0,know = 0,kplus = 0,mask;
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
	uint64 kmin_file, know_file = 0, kmax_file = 0;
	uint32 passmin_file = 0, passnow_file = 0, passmax_file = 0;

	uint64 sweep,interval_lo,interval_now,interval_hi,ninterval;

/*   this stuff is for the small-primes sieve...	*/
	uint32 max_diff;
	uint64 dstartval;
	/*
	unsigned char pdiff[NUM_SIEVING_PRIME];	// Compact table storing the difference between adjacent odd primes.
	uint32 offset[NUM_SIEVING_PRIME], startval[NUM_SIEVING_PRIME], startval_incr[NUM_SIEVING_PRIME];
	*/
	unsigned char *pdiff;	/* Compact table storing the (difference/2) between adjacent odd primes.
							http://mathworld.wolfram.com/PrimeGaps.html shows the first >256-gap at ~387 million
							and the first >512-gap at ~300 billion, so using the half-of-even-gap trick makes
							the difference between an algo which is safe for all 32-bit primes and one with a
							limit of < 2^30. */
	uint32 *offset, *startval, *pinv, *kdeep = 0x0, ndeep = 0;
	uint32 MAX_SIEVING_PRIME = 0;
	uint64 *u64_arr = 0x0;	/* generic array allowing us to hook into the mi64 routines */

#ifdef P1WORD
	double twop_float = 0,fqlo,fqhi;
	uint128 p128,q128,t128;	// Despite the naming, these are needed for nominal 1-word runs with moduli exceeding 64 bits
#endif
#ifdef P3WORD
	uint192 p192,q192,t192;
  #if(defined(USE_FLOAT))
	uint256 x256;	// Needed to hold result of twopmodq200_8WORD_DOUBLE
  #endif
#endif
#ifdef P4WORD
	uint256 p256,q256,t256;
#endif
#if FAC_DEBUG
	uint32 on_bits = 0;
	/* Set k_targ to some known-factor k to debug a specific missed-factor case: */
	uint64 k_targ = 0;
	uint32 pass_targ = 0xffffffff;	/* Init to a known-invalid value; if user specifies
	 								a known test factor via k_targ, pass_targ will
	 								be set to the associated (legitimate) value between 0 and 15;
	 								(pass_targ < 16) can subsequently be used as a quick
	 								test for whether a known test factor has been set. */
#endif

#if DBG_SIEVE
	uint32 *startval_incr;
	uint32 i64_targ, bit_targ;
#endif

	/*...time-related stuff	*/
	clock_t clock1, clock2;
	double td, tdiff;
#if TEST_TRIALDIV
	double citer;
#endif
	/* printf character buffers - when using to print args in a single printf,
	need a separate buffer for each arg: */
	char char_buf0[STR_MAX_LEN], char_buf1[STR_MAX_LEN], char_buf2[STR_MAX_LEN];
	char *char_addr;

/* Set == 1 to test the trial-div stuff: */
#define	TEST_TRIALDIV	0
#if TEST_TRIALDIV
	#define MAX_ARRAY_DIM 10000
	uint32	vec_len = MAX_ARRAY_DIM;
	uint64*	xvec = (uint64 *)calloc(MAX_ARRAY_DIM, sizeof(uint64));
	uint32 tryq[8];
#endif

#ifdef macintosh
	argc = ccommand(&argv);	/* Macintosh CW */
#endif

/* Get thread count and use to allocate the row dimension of the [NTHREAD] x [TRYQ] factor-candidates array: */
	uint32 NTHREAD = 1;	/* TODO: Replace with dynamic thread-count init */

	k_to_try = ALLOC_UINT64(k_to_try, NTHREAD*TRYQ);
	k_to_try = ALIGN_UINT64(k_to_try);	/* Don't care about possible small memleak here */
	ASSERT(HERE, k_to_try != 0, "k_to_try alloc failes!");

#ifdef USE_GPU
	// Alloc host copy of kvec:
	cpu_kvec = ALLOC_UINT64(cpu_kvec, kvsz);	cpu_kvec = ALIGN_UINT64(cpu_kvec);
	ASSERT(HERE, (cpu_kvec != 0) && (((uint64)cpu_kvec & 7) == 0), "cpu_kvec alloc/align failed!");
	// GPU device copy; this malloc is bytewise:
	cudaMalloc(&gpu_kvec, (kvsz << 3));

	// Allocate output vector (resulting 2^p mod q, in binary "is factor?" form) in host memory:
	cpu_rvec = (uint8 *)malloc(kvsz);	// Until impl packed-bitmap scheme for device code return values, use byte array for return values
	cudaMalloc(&gpu_rvec, kvsz);
#endif

/* Allocate factor_k array and align on 16-byte boundary: */
	factor_ptmp = ALLOC_UINT64(factor_ptmp, 24);
	factor_k = ALIGN_UINT64(factor_ptmp);	factor_ptmp = 0x0;
	ASSERT(HERE, ((uint64)factor_k & 0x3f) == 0, "factor_k not 64-byte aligned!");

/*...initialize logicals and factoring parameters...	*/
	restart=FALSE;

#ifdef FACTOR_STANDALONE
	host_init();
#endif

#if INCLUDE_PM1
	/* Simple self-tester for GCD routines in gcd_lehmer.c: */
	fprintf(stderr, "INFO: testing GCD routines...\n");
	if(test_gcd() != 0)
	{
		sprintf(cbuf, "Factor_init : GCD test failed.\n");
		fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	} else {
		exit(0);
	}
#endif

/* Do a quick series of self-tests: */
#if 1//FAC_DEBUG
	test_fac();
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

	-passmin [int]  Maximum factoring pass for the run (0-15, default =  0).
	-passmax [int]  Maximum factoring pass for the run (0-15, default = 15).
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

	-q [int]    A known factor for the number (only used if FAC_DEBUG = 1).
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
/*** TODO: implement the needed modular powering routines ***/
ASSERT(HERE, 0,"Fermat-number trial factoring not currently supported.");

			strncpy(pstring, argv[nargs++], STR_MAX_LEN);
			MODULUS_TYPE = MODULUS_TYPE_FERMAT;
		}

		/* Factor bounds, in log2(qmin/qmax) (floating double) form: */
		else if(STREQ(stFlag, "-bmin"))
		{
		#if(defined(P1WORD))
			if(kmin || kmax || kplus)
			{
				fprintf(stderr,"*** ERROR: If -kmin/kmax or -kplus used to set bounds for factoring, -bmin [and -bmax] disallowed.\n");
				goto MFACTOR_HELP;
			}
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			bmin = convert_base10_char_double(stFlag);
		  #if FAC_DEBUG
			printf("bmin = %lf\n", bmin);
		  #endif
		#else
			ASSERT(HERE, 0,"bmin/bmax form of bounds-setting only allowed for single-word-p case!");
		#endif
		}
		else if(STREQ(stFlag, "-bmax"))
		{
		#if(defined(P1WORD))
			if(kmin || kmax || kplus)
			{
				fprintf(stderr,"*** ERROR: If -kmin/kmax or -kplus used to set bounds for factoring, -bmax [and -bmin] disallowed.\n");
				goto MFACTOR_HELP;
			}
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			bmax = convert_base10_char_double(stFlag);
		  #if FAC_DEBUG
			printf("bmax = %lf\n", bmax);
		  #endif
		#else
			ASSERT(HERE, 0,"bmin/bmax form of bounds-setting only allowed for single-word-p case!");
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
			ASSERT(HERE, passmin < FACTOR_PASS_MAX,"factor.c: passmin < FACTOR_PASS_MAX");
		}
		else if(STREQ(stFlag, "-passmax"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);
			passmax = (uint32)convert_base10_char_uint64(stFlag);
			ASSERT(HERE, passmax < FACTOR_PASS_MAX,"factor.c: passmax < FACTOR_PASS_MAX");
			ASSERT(HERE, passmax >= passmin       ,"factor.c: passmax >= passmin");
		}

		/* Come again? */
		else
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
	if(!fp)
	{
		fprintf(stderr,"ERROR: Unable to open statfile %s for writing.\n",STATFILE);
		ASSERT(HERE, 0,"0");
	}
	else
	{
		fclose(fp); fp = 0x0;
	}

#endif	/* #ifdef FACTOR_STANDALONE */

	/* Make sure a valid exponent string has been given - if this is the only
	command-line parameter, will attempt to read the other needed run parameters
	from the corresponding checkpoint file:
	*/
	ASSERT(HERE, STRNEQ(pstring,""),"factor.c : pstring empty!");

	/* -bmin/bmax used to set bounds for factoring: */
	if(bmin || bmax)
	{
		ASSERT(HERE, (kmin==0 && kmax==0 && kplus==0),"(kmin==0 && kmax==0 && kplus==0)");

		if(bmin < 0)
		{
			fprintf(stderr,"ERROR: log2(min factor) must be >= 0. Offending entry = %lf.\n", bmin);
			ASSERT(HERE, 0,"0");
		}
		else if(bmin >= MAX_BITS_Q)
		{
			fprintf(stderr,"ERROR: log2(min factor) exceeds allowable limit of %u. Offending entry = %lf.\n", MAX_BITS_Q, bmin);
			ASSERT(HERE, 0,"0");
		}

		if(bmax <= 0)
		{
			fprintf(stderr,"ERROR: log2(max factor) must be > 0. Offending entry = %lf.\n", bmax);
			ASSERT(HERE, 0,"0");
		}
		else if(bmax > MAX_BITS_Q)
		{
			fprintf(stderr,"ERROR: log2(max factor) exceeds allowable limit of %u. Offending entry = %lf.\n", MAX_BITS_Q, bmax);
			ASSERT(HERE, 0,"0");
		}

		if(bmax < bmin)
		{
			fprintf(stderr,"ERROR: (bmax = %lf) < (bmin = %lf)!\n", bmax, bmin);
			ASSERT(HERE, 0,"0");
		}
	}

	/* -kmin/kmax used to set bounds for factoring: */
	if(kmin || kmax)
	{
		ASSERT(HERE, kmax != 0 ,"factor.c: kmax not set!");
		ASSERT(HERE, (int64)kmax > 0, "kmax must be 63 bits or less!");
		ASSERT(HERE, (bmin==0 && bmax==0 && kplus==0),"(bmin==0 && bmax==0 && kplus==0)");

		if(kmax < kmin)
		{
			fprintf(stderr,"ERROR: (kmax = %s) < (kmin = %s)!\n", &char_buf0[convert_uint64_base10_char(char_buf0, kmax)], &char_buf1[convert_uint64_base10_char(char_buf1, kmin)]);
			ASSERT(HERE, 0,"0");
		}
	}

	ASSERT(HERE, bmax > 0.0 || kmax != 0 ,"factor.c: One of bmax or kmax must be set!");

/* Temporary - needed until I finish implementing Fermat-number capability: */
ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		  || (MODULUS_TYPE == MODULUS_TYPE_MERSMERS)
		/*|| (MODULUS_TYPE ==   MODULUS_TYPE_FERMAT)*/
			, "Unsupported modulus type!");

	// Convert power-of-2 exponent to unsigned int form and allocate the exponent-storage vector.
	// We use MAX_BITS_P (defined in Mdata.h) to set the allocated storage here, but use the user-set
	// exponent to set the number of words of that allocate storage which are actually used:
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
	#if FAC_DEBUG
		printf("%s(%s) = M(p) with p = %s\n", NUM_PREFIX[MODULUS_TYPE], pstring, &char_buf0[convert_mi64_base10_char(char_buf0, p, lenP, 0)]);
	#endif
	}
	else
	{
		/* Convert stringified exponent to mi64 form and check the length of p and of the max. factor: */
		p = convert_base10_char_mi64(pstring, &lenP);	// This does the mem-alloc for us in this case
		ASSERT(HERE, lenP > 0, "factor.c: Error converting pstring!");
		nbits_in_p = (lenP<<6) - mi64_leadz(p, lenP);
	}

	// Allocate the other modulus-dependent vectors:
	lenQ = ((uint32)MAX_BITS_Q + 63)>>6;	// This sets upper bound on #words needed to store max. factor candidate

	two_p = (uint64 *)calloc(lenQ, sizeof(uint64));
	p120  = (uint64 *)calloc(lenQ, sizeof(uint64));
	q     = (uint64 *)calloc(lenQ, sizeof(uint64));
	q2    = (uint64 *)calloc(lenQ, sizeof(uint64));
	u64_arr = (uint64 *)calloc(lenQ, sizeof(uint64));

	// Now use the just-allocated vector storage to compute how many words are really needed for qmax.
	// Since the sieving always proceeds in full passes through the bit-cleared sieve, the actual kmax used
	// may be up to (len*64)-1 larger than the user-specified kmax:
	if(kmax)
	{
		interval_hi = (uint64)ceil((double)(kmax>>6)/len);	// Copied from restart-file code below
		u64_arr[lenP] = mi64_mul_scalar( p, 2*interval_hi*(len<<6), u64_arr, lenP);	// Actual kmax used at runtime = interval_hi*(len<<6);
		lenQ = lenP + (u64_arr[lenP] != 0);
	}
	else
	{
		lenQ = ( (uint32)(ceil(bmax)) + 63 ) >> 6;
	}

	if((p[0] & 1) == 0)
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
		mi64_shrl(p, q, findex, lenP);
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
			mi64_shrl(q, q, findex, lenP);
			mi64_sub_scalar(q, 1ull, q, lenP);
			ASSERT(HERE, mi64_iszero(q, lenP), "Internal M(M(p))-exponent fails all-binary-ones check!");
		}
		// We can use a lookup table vs known M(p) for all cases, but if Mersenne or M(M(p)) with suitably small p,
		// add a base-2 Fermat PRP test, more as a self-test of the various modpow routines than anything else:
		if(lenP < 1000) {
			mi64_sub_scalar(p, 1ull, q, lenP);	/* q = p-1 */
			if(!mi64_twopmodq(q, lenP, 0, p, lenP, 0x0))
			{
				fprintf(stderr,"p is not prime! Offending p = %s\n", pstring);
				ASSERT(HERE, 0,"0");
			}
		}
	}

	/* 2*p: Don't need to worry about overflow here since we've allocated two*p, p120, q, etc to be of lenQ, not lenP: */
	two_p[lenP] = mi64_add(p, p, two_p, lenP);

	/* 120*p: */
	p120[lenP] = mi64_mul_scalar(p, 120ull, p120, lenP);

  #if FAC_DEBUG
	printf("two_p    = %s\n", &char_buf0[convert_mi64_base10_char(char_buf0, two_p, lenP, 0)]);
	printf("120*p    = %s\n", &char_buf0[convert_mi64_base10_char(char_buf0, p120 , lenP, 0)]);
  #endif

	/* p mod 60: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS)
	{
		pmod60 = twopmmodq64(findex, 60ull) - 1;		// For double-Mersenne factoring need M(p) mod 60
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
		printf("p mod 60 = %d\n", pmod60);
		if(pmod60 != mi64_div_y32(p, 60, 0x0, lenP)) {
			printf("p mod 60 B = %d\n", mi64_div_y32(p, 60, 0x0, lenP));
			ASSERT(HERE, 0, "Powering and direct-long-dive give differing 9p % 60) values!");
		}

	} else {
		pmod60 = mi64_div_y32(p, 60, 0x0, lenP);
	}

	// If user-set kmax, test factoring range vs internal limits
	if(kmax) {
		interval_hi = (uint64)ceil((double)(kmax>>6)/len);	// Copied from restart-file code below
		u64_arr[lenP] = mi64_mul_scalar( p, 2*interval_hi*(len<<6), u64_arr, lenP);	// Actual kmax used at runtime = interval_hi*(len<<6);
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
	CMASKBITS = (30-bits_in_pq2);
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

	/* TF restart files are in HRF, not binary: */
	fp = mlucas_fopen(RESTARTFILE, "r");
	if(!fp)
	{
									fprintf(stderr,"INFO: No factoring savefile %s found ... starting from scratch.\n",RESTARTFILE);
	#ifndef FACTOR_STANDALONE
		fq = mlucas_fopen(STATFILE,"a");	fprintf(	fq,"INFO: No factoring savefile %s found ... starting from scratch.\n",RESTARTFILE);	fclose(fq); fq = 0x0;
	#endif
	}
	else
	{
									fprintf(stderr,"Factoring savefile %s found ... reading ...\n",RESTARTFILE);
	#ifndef FACTOR_STANDALONE
		fq = mlucas_fopen(STATFILE,"a");	fprintf(	fq,"Factoring savefile %s found ... reading ...\n",RESTARTFILE);	fclose(fq); fq = 0x0;
	#endif
		/* The factoring checkpoint file is assumed to have the format:
			Line 1:		{string containing the current exponent stored in pstring}

			Line 2:		{Log2(minimum factor to try), in floating double form}
						If > 10^9 its whole-number part is taken as the KMin value instead.
			Line 3:		{Log2(maximum factor to try), in floating double form}
						If > 10^9 its whole-number part is taken as the KMax value instead.

			Line 4:		KMin = {smallest factor K value to be tried in each pass}
			Line 5:		KNow = { largest factor K value tried so far during current pass}
			Line 6:		KMax = { largest factor K value to be tried in each pass}

			Line 7:		PassMin = {maximum pass for the run (typically 15, but perhaps not, e.g. for a factoring assignment split over multiple CPUs.}
			Line 8:		PassNow = {current factoring pass}
			Line 9:		PassMax = {maximum pass for the run (typically 15, but perhaps not, e.g. for a factoring assignment split over multiple CPUs.}

			Line 10:	Number of q's tried so far during the run
			Line 11:	64-bit (sum of trial q)%2^64 checksum

			Line 12+:	Any diagnostic info not needed for restarting from interrupt
						(mainly, in standalone mode can use this in place of STATFILE.)
		*/
		curr_line = 0;

		/* pstring*/
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (current exponent) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		/* Strip the expected newline char from in_line: */
		char_addr = strstr(in_line, "\n");
		if(char_addr)
			*char_addr = '\0';
		/* Make sure restart-file and current-run pstring match: */
		if(STRNEQ(in_line, pstring))
		{
			fprintf(stderr,"ERROR: current exponent %s != Line %d of factoring restart file %s!\n",pstring, curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}

		/* bmin */
		++curr_line;
		fgets(cbuf, STR_MAX_LEN, fp);
		itmp = sscanf(cbuf, "%lf", &bmin_file);
		if(itmp != 1)
		{
			fprintf(stderr,"ERROR: unable to parse Line %d (bmin) of factoring restart file %s. Offending input = %s\n", curr_line, RESTARTFILE, cbuf);
			ASSERT(HERE, 0,"0");
		}

		/* bmax */
		++curr_line;
		fgets(cbuf, STR_MAX_LEN, fp);
		itmp = sscanf(cbuf, "%lf", &bmax_file);
		if(itmp != 1)
		{
			fprintf(stderr,"ERROR: unable to parse Line %d (bmin) of factoring restart file %s. Offending input = %s\n", curr_line, RESTARTFILE, cbuf);
			ASSERT(HERE, 0,"0");
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
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: 'KMin' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "KMin");
		/* Since the preceding fscanf call may leave us at the end of curr_line-1
		(rather than the beginning of curr_line), allow for a possible 2nd needed
		fgets call here: */
		if(!char_addr)
		{
			goto GET_LINE4;
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			kmin_file = convert_base10_char_uint64(char_addr);
		}

		/* KNow */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (KNow) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "KNow");
		if(!char_addr)
		{
			fprintf(stderr,"ERROR: 'KNow' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			know_file = convert_base10_char_uint64(char_addr);
		}

		/* KMax */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (KMax) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "KMax");
		if(!char_addr)
		{
			fprintf(stderr,"ERROR: 'KMax' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			kmax_file = convert_base10_char_uint64(char_addr);
		}

		/* PassMin */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (PassMin) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "PassMin");
		if(!char_addr)
		{
			fprintf(stderr,"ERROR: 'PassMin' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			passmin_file = (uint32)convert_base10_char_uint64(char_addr);
			ASSERT(HERE, passmin_file < FACTOR_PASS_MAX,"factor.c: passmin < FACTOR_PASS_MAX");
		}

		/* PassNow */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (PassNow) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "PassNow");
		if(!char_addr)
		{
			fprintf(stderr,"ERROR: 'PassNow' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			passnow_file = (uint32)convert_base10_char_uint64(char_addr);
			ASSERT(HERE, passnow_file < FACTOR_PASS_MAX,"factor.c: passnow < FACTOR_PASS_MAX");
			ASSERT(HERE, passnow_file >= passmin_file  ,"factor.c: passnow_file >= passmin_file");
		}

		/* PassMax */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (PassMax) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "PassMax");
		if(!char_addr)
		{
			fprintf(stderr,"ERROR: 'PassMax' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			passmax_file = (uint32)convert_base10_char_uint64(char_addr);
			ASSERT(HERE, passmax_file < FACTOR_PASS_MAX,"factor.c: passmax_file < FACTOR_PASS_MAX");
			ASSERT(HERE, passmax_file >= passnow_file  ,"factor.c: passmax_file >= passnow_file");
		}

		/* Number of q's tried: */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (#Q tried) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "#Q tried");
		if(!char_addr)
		{
			fprintf(stderr,"ERROR: '#Q tried' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			count = convert_base10_char_uint64(char_addr);	// Need to reset == 0 prior to sieving so kvector-fill code works properly
		}

		/* 64-bit (sum of trial q)%2^64 checksum: */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (Checksum1) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "Checksum1");
		if(!char_addr)
		{
			fprintf(stderr,"ERROR: 'Checksum1' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			itmp = sscanf(char_addr, "%16llX", &checksum1);
			if(itmp != 1)
			{
				fprintf(stderr,"ERROR: unable to parse hex entry of Line %d (Checksum1) of factoring restart file %s. Offending input = %s\n", curr_line, RESTARTFILE, in_line);
				ASSERT(HERE, 0,"0");
			}
		}

		/* 64-bit (sum of 2^p % q)%2^64 checksum: */
		++curr_line;
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			fprintf(stderr,"ERROR: unable to read Line %d (Checksum2) of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		char_addr = strstr(in_line, "Checksum2");
		if(!char_addr)
		{
			fprintf(stderr,"ERROR: 'Checksum2' not found in Line %d of factoring restart file %s!\n", curr_line, RESTARTFILE);
			ASSERT(HERE, 0,"0");
		}
		else
		{
			char_addr = strstr(in_line, "=");
			if(!char_addr)
			{
				fprintf(stderr,"ERROR: Line %d of factoring restart file %s lacks the required = sign!\n", curr_line, RESTARTFILE);
				ASSERT(HERE, 0,"0");
			}
			char_addr++;
			itmp = sscanf(char_addr, "%16llX", &checksum2);
			if(itmp != 1)
			{
				fprintf(stderr,"ERROR: unable to parse hex entry of Line %d (Checksum2) of factoring restart file %s. Offending input = %s\n", curr_line, RESTARTFILE, in_line);
				ASSERT(HERE, 0,"0");
			}
		}

		fclose(fp); fp = 0x0;

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
				ASSERT(HERE, 0,"bmin/bmax form of bounds-setting only allowed for single-word-p case!");
			#endif
				ASSERT(HERE, (kmin==0 && kmax==0 && kplus==0),"(kmin==0 && kmax==0 && kplus==0) - please delete any restart files for this p and retry debug run.");

				if(bmin)
				{
					ASSERT(HERE, bmin >= bmin_file - 0.0000000001,"bmin >= bmin_file");
					if(bmin < bmax_file)
					{
						fprintf(stderr,"WARNING: Specified bmin (%lf) smaller than previous-run bmax = %lf. Setting equal to avoid overlapping runs.\n", bmin, bmax_file);
					}
				}
				bmin = bmax_file;

				/* We expect any command-line bmax will be > that in the restart file: */
				if(bmax)
				{
					ASSERT(HERE, bmax > bmax_file - 0.0000000001,"bmax >= bmax_file");
				}
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

				if(kmin)
				{
					ASSERT(HERE, kmin >= kmin_file,"kmin >= kmin_file");
					if(kmin < kmax_file)
					{
						fprintf(stderr,"WARNING: Specified kmin (%s) smaller than previous-run kmax = %s. Setting equal to avoid overlapping runs.\n", &char_buf0[convert_uint64_base10_char(char_buf0, kmax)], &char_buf1[convert_uint64_base10_char(char_buf1, kmax_file)]);
					}
				}
				kmin = kmax_file;

				/* We expect any command-line kmax will be > that in the restart file: */
				if(kmax)
				{
					ASSERT(HERE, kmax > kmax_file,"kmax >= kmax_file");
				}
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
		restart=TRUE;
	}

/************************ END(RESTART STUFF) *******************/

	/* If it's not a restart of an as-yet-uncompleted run, synchronize the factoring-bound params: */
	if(!incomplete_run)
	{
		/* Double-check factoring pass bounds: */
		if(passmin > (FACTOR_PASS_MAX-1) )
		{
			fprintf(stderr,"ERROR: passmin must be <= %u. Offending entry = %u.\n", FACTOR_PASS_MAX-1, passmin);
			ASSERT(HERE, 0,"0");
		}

		if(passmax < passmin)
		{
			fprintf(stderr,"ERROR: (passmax = %u) < (passmin = %u)!\n", passmax, passmin);
			ASSERT(HERE, 0,"0");
		}
		if(passmax > (FACTOR_PASS_MAX-1) )
		{
			fprintf(stderr,"ERROR: passmax must be <= %u. Offending entry = %u.\n", FACTOR_PASS_MAX-1, passmax);
			ASSERT(HERE, 0,"0");
		}

		/**** Process factor candidate bounds: ****/

		/* If any of bmin|kmin, bmax|kmax nonzero, calculate its counterpart: */
	#if(defined(P1WORD))
		/* Find FP approximation to 2*p - can't use this for multiword case, because double approximation tp 2*p may overflow: */
		twop_float = (double)p[0]+p[0];
		/* Now that have twop_float, compute kmax if not already set: */
		if(!kmax)
		{
			fqlo = pow(2.0, bmax);
			kmax = (uint64)(fqlo/twop_float);
		}
		if(kmin || bmin)
		{
			if(kmin == 0)	/* Lower Bound given in log2rithmic form */
			{
				ASSERT(HERE, bmin <= bmax, "bmin >= bmax!");
				fqlo = pow(2.0, bmin);
				kmin = (uint64)(fqlo/twop_float);
			}
			else
			{
				ASSERT(HERE, kmin <= kmax, "kmin >= kmax!");
				fqlo = kmin*twop_float + 1.0;
				bmin = log(fqlo)*ILG2;
			}
		}
		else
		{
			fqlo = 1.0;
		}

		if(kmax || bmax)
		{
			if(kmax == 0)	/* Upper Bound given in log2rithmic form */
			{
				fqhi = pow(2.0, bmax);
				kmax = (uint64)(fqhi/twop_float);
			}
			else
			{
				fqhi = kmax*twop_float + 1.0;
				bmax = log(fqhi)*ILG2;
			}
		}
		else
		{
			ASSERT(HERE, 0 ,"factor.c : One of bmax, kmax must be nonzero!");
		}
	#endif

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

	bit_map = (uint64 *)calloc(len/60+1, sizeof(uint64));
	if (bit_map == NULL) {
		fprintf(stderr,"Memory allocation failure for BITMAP array");
		ASSERT(HERE, 0,"0");
	}

	bit_atlas = (uint64 *)calloc((len/60+1)*16, sizeof(uint64));
	if (bit_atlas == NULL) {
		fprintf(stderr,"Memory allocation failure for TEMPLATE array");
		ASSERT(HERE, 0,"0");
	}

	pdiff = (unsigned char *)calloc(NUM_SIEVING_PRIME, sizeof(unsigned char));
	if (pdiff == NULL) {
		fprintf(stderr,"Memory allocation failure for pdiff array");
		ASSERT(HERE, 0,"0");
	}

	offset = (uint32 *)calloc(NUM_SIEVING_PRIME, sizeof(uint32));
	if (offset == NULL) {
		fprintf(stderr,"Memory allocation failure for OFFSET array");
		ASSERT(HERE, 0,"0");
	}

	startval = (uint32 *)calloc(NUM_SIEVING_PRIME, sizeof(uint32));
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

	/*
	switch(MAX_SIEVING_PRIME){
		case 0:
			nprime=10; break;
		case 10000:
			nprime=1228; break;
		case 20000:
			nprime=1228+1033; break;
		case 30000:
			nprime=1228+1033+983; break;
		case 40000:
			nprime=1228+1033+983+958; break;
		case 50000:
			nprime=1228+1033+983+958+930; break;
		case 60000:
			nprime=1228+1033+983+958+930+924; break;
		case 70000:
			nprime=1228+1033+983+958+930+924+485; break;
		default:
			fprintf(stderr,"ERROR: MAX_SIEVING_PRIME must be a multiple of 10000, max = 70000.\n");
			ASSERT(HERE, 0,"0");
	}
	offset = (uint32 *)calloc(nprime, sizeof(uint32));
	*/

	/*
	nprime = NUM_SIEVING_PRIME;
	int PWORDS = NWORD;	// Number of 64-bit words needed to store p
	for(i = 0; i < NUM_SIEVING_PRIME[PWORDS-1]; i++)
	{
		pdiff[i] = 0;
	}
	*/

	#if 1
		/* Check integrity (at least in the sense of monotonicity) for the precomputed pseudoprime table: */
		for(i = 1; i < 9366; ++i)
		{
			ASSERT(HERE, f2psp[i] > f2psp[i-1],"Misplaced pseudoprime!");
		}
	#endif

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
	#if ALLOW_PRP
		fprintf(stderr,"Base-2 Fermat pseudoprimes = \n");
	#endif

		curr_p = 3;	/* Current prime stored in l. */
		max_diff = 0;

	#if ALLOW_PRP
		num_pseudo = 0;
	#else
		f2psp_idx = 0;	/* Index to next-expected Fermat base-2 pseudoprime in the precomputed table */
	#endif

	#if ALLOW_PRP
		/* pdiff[0] = 0 represents the diff between 2 and 3, so start loop at i=1: */
		ihi = curr_p = 3;
		for(i = 1; i < nprime; ++i)
		{
			curr_p += 2;
			++pdiff[i];
	#else
		/* Init first few diffs between 3/5, 5/7, 7/11, so can start loop with curr_p = 11 == 1 (mod 10), as required by twopmodq32_x8(): */
		pdiff[1] = 1;
		pdiff[2] = 1;
		ihi = curr_p = 11;
		/* Process chunks of length 30, starting with curr_p == 11 (mod 30). Applying the obvious divide-by-3,5 mini-sieve,
		we have 8 candidates in each interval: curr_p + [ 0, 2, 6, 8,12,18,20,26].
		For example: curr_p = 11 gives the 8 candidates: 11,13,17,19,23,29,31,37.
		*/
		for(i = 3; i < nprime; curr_p += 30)
		{
	#endif
			/* Make sure (curr_p + 29) < 2^32: */
			if(curr_p > 0xffffffe3)
			{
				fprintf(stderr,"curr_p overflows 32 bits!");
				nprime = i;
				break;
			}

			/* Max sieving prime must be < smallest candidate factor of M(p) */
		#if(defined(P1WORD))
			if((curr_p+29) > two_p[0])
			{
				nprime = i;
				break;
			}
		#endif

		#if !ALLOW_PRP

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
			continue;

		#elif ALLOW_PRP

			if(twopmodq32(curr_p-1, curr_p) == 0)
			{
				--i;
				continue;
			}

			if(curr_p < 341)
			{
				if(pdiff[i] > max_diff)
					max_diff = pdiff[i];
				continue;
			}

		  #if 0
			/* Modding w.r.to the small primes up to 43 eliminates roughly half the pseudoprimes: */
			if(curr_p%3 == 0 || curr_p%5 == 0 || curr_p%7 == 0 || curr_p%11== 0 || curr_p%13== 0 || curr_p%17== 0 || curr_p%19== 0 || curr_p%23== 0 || curr_p%29== 0 || curr_p%31== 0 || curr_p%37== 0 || curr_p%41== 0 || curr_p%43== 0)
			{
				++num_pseudo;
				fprintf(stderr,"%u[<47]...",curr_p);
				--i;
				continue;
			}
			l = 47;
			idx = 13;
		  #elif 1
			/* To speed the small-prime sieve init, test divisibility by smaller primes via 64-bit gcd rather than serial modulo.
				Use that a uint64 can store the product of odd primes 3 through 53 = 16294579238595022365.
			*/
			itmp64 = 16294579238595022365ull;
			j = (uint32)gcd64((uint64)curr_p, itmp64);
			if(j > 1)
			{
				++num_pseudo;
	//	fprintf(stderr,"%u...",curr_p);
				--i;
				continue;
			}
			l = 59;
			idx = 15;
		  #else
			l = 3;
			idx = 0;
		  #endif

			/* Trial-divide the next-prime candidate by all odds <= sqrt(curr_p).
			If none divide it, curr_p is prime - increment outer-loop parameters.
			If curr_p not prime, move on to next odd number; repeat until we find next odd prime. */
			sqrt_p = (uint32)(sqrt(1.0*curr_p) + 0.5);	/* Add a small fudge factor to the result of the sqrt to make sure RO errors don't
								cause the sqrt to be one smaller than it should be once we truncate to integer. */
			isprime = TRUE;

			while(l <= sqrt_p)
			{
				/* It's composite - go to next-larger odd number */
				if(curr_p%l == 0)
				{
					++num_pseudo;
		//	fprintf(stderr,"%u...",curr_p);
					isprime = FALSE;
					break;
				}
				l += (pdiff[++idx] << 1);
			}
			if(!isprime)
			{
				--i;
			}
			else
			{
				ihi = curr_p;
				if(pdiff[i] > max_diff)
				{
					max_diff = pdiff[i];
				#if 0
					printf("\npdiff = %d at curr_p = %u\n", 2*max_diff,ihi);
				#endif
				}
			}

		#elif 0	// A bit of test code:

			ASSERT(HERE, curr_p <= f2psp[f2psp_idx],"Error in pseudoprime sieve");
			if(curr_p == f2psp[f2psp_idx])
			{
				++f2psp_idx;
				--i;
				continue;
			}
			else
			{
				ihi = curr_p;
				if(pdiff[i] > max_diff)
					max_diff = pdiff[i];
			}

		#endif	// ALLOW_PRP

		}
		MAX_SIEVING_PRIME = ihi;

	#if 1//FAC_DEBUG
		printf("Using first %u odd primes; max gap = %u\n",nprime,2*max_diff);
		printf("max sieving prime = %u\n",MAX_SIEVING_PRIME);
	  #if ALLOW_PRP
		printf("number of base-2 Fermat pseudoprimes found = %u\n", num_pseudo);
	  #endif
	#endif

#if FERMAT_PSP2
	#include "f2psp.txt"
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
	#if FAC_DEBUG
		ASSERT(HERE, curr_p == prime[nprime], "factor.c : curr_p == prime[nprime]");
	#endif
	}
	MAX_SIEVING_PRIME = curr_p;
	***********/

	/****************** KNOWN-TEST-FACTOR STUFF: *******************/
#if FAC_DEBUG

	if(k_targ)
	{
		// Could add check of whether associated q is prime, but assume user knows what he's soing,
		// and furthermore may be useful to allow for composite products-of-smaller-prime-factors:
		kmod60 = k_targ%60;

		printf("p mod 60 = %d\n", pmod60);
		printf("k mod 60 = %d\n", kmod60);

		/* ...and get the pass number on which the factor should be found.
		(Remember that the pass number returned by CHECK_PKMOD60 is unit-offset).
		If a known factor given, only process the given k/log2 range for that pass:
		*/
		pass_targ = CHECK_PKMOD60(pmod60, kmod60) - 1;
	}
#endif	/* end #if(FAC_DEBUG) (q_targ processing) */
/* Uncomment and set k if want to test a given p,k pair for legality:
	k = 201;
	kmod60 = k%60;
	printf("pmod60, kmod60 = %u, %u\n",pmod60,kmod60);
	ASSERT(HERE, CHECK_PKMOD60(pmod60, kmod60), "K must be one of the admissible values (mod 60) for this p (mod 60)!");
	exit(0);
*/
	/* The 16 relevant p mod 60 cases. It's times like these that I really long for f90-style vector array initializations: */
	i = 0;
	switch(pmod60){
		case  1:incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1; break;
		case  7:incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3; break;
		case 11:incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11; break;
		case 13:incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4; break;
		case 17:incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5; break;
		case 19:incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3; break;
		case 23:incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3; break;
		case 29:incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5; break;
		case 31:incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4; break;
		case 37:incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1; break;
		case 41:incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5; break;
		case 43:incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3; break;
		case 47:incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3; break;
		case 49:incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1; break;
		case 53:incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5; break;
		case 59:incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3;incr[i++]= 4;incr[i++]= 5;incr[i++]= 3;incr[i++]= 1;incr[i++]=11;incr[i++]= 1;incr[i++]= 3;incr[i++]= 5;incr[i++]= 4;incr[i++]= 3;incr[i++]= 5;incr[i++]= 3; break;
		default:
			fprintf(stderr,"%u not an acceptable value for P mod 60. P likely nonprime.\n", pmod60); ASSERT(HERE, 0,"0");
	}

	/* If it's a restart, interval_lo for the initial pass will be based
	on (know), rather than (kmin) - handle that just inside the pass-loop: */
	interval_lo = (kmin>>6)/(uint64)len;
	interval_now= (know>>6)/(uint64)len;
	interval_hi = (uint64)ceil((double)(kmax>>6)/len);

	/* Make sure we always do at least one full pass through the sieve
	(e.g. if kmin = kmax = 1, ceil(kmax>>6) gives 0, same as (kmin>>6): */
	if(interval_hi == interval_lo)
		interval_hi += 1;

	ninterval = interval_hi - interval_lo;

	/* The actual kmin/kmax values used in the run are exact multiples
	of len*64 - we always do at least one full pass through the sieve,
	i.e. we run through at least (len*64) worth of k's */
	kmin = interval_lo *(len<<6);
	know = interval_now*(len<<6);
	kmax = interval_hi *(len<<6);

	/* And now that we have the actual kmin/kmax, recalculate these: */
  #if(defined(P1WORD))
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

#if FAC_DEBUG
	/* Make sure the range of k's for the run contains any target factor: */
	if(k_targ)
		ASSERT(HERE, (kmin <= k_targ) && (kmax >= k_targ),"k_targ not in [kmin, kmax]");
#endif

	ASSERT(HERE, fp == 0x0,"0");
#ifdef FACTOR_STANDALONE
	fp = stderr;
#else
	fp = mlucas_fopen(STATFILE,"a");
#endif
	fq = mlucas_fopen(OFILE,"a");

#ifdef P1WORD
	sprintf(char_buf0, "searching in the interval k=[%s, %s], i.e. q=[%e, %e]\n", &char_buf1[convert_uint64_base10_char(char_buf1, kmin )], &char_buf2[convert_uint64_base10_char(char_buf2, kmax )],fqlo,fqhi);
#else
	sprintf(char_buf0, "searching in the interval k=[%s, %s]\n", &char_buf1[convert_uint64_base10_char(char_buf1, kmin )], &char_buf2[convert_uint64_base10_char(char_buf2, kmax )]);
#endif
	fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);

	sprintf(char_buf0, "each of %2u (p mod 60) passes will consist of %s intervals of length %u\n", passmax-passmin+1, &char_buf1[convert_uint64_base10_char(char_buf1, ninterval)], bit_len);
	fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);

	if(passnow != passmin || know != kmin)
	{
		sprintf(char_buf0, "Resuming execution with pass %u and k = %s\n", passnow, &char_buf1[convert_uint64_base10_char(char_buf1, know )]);
		fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);
		sprintf(char_buf0, "#Q tried = %s\n", &char_buf1[convert_uint64_base10_char (char_buf1, count)] );
	count = 0;	// Reset == 0 prior to sieving so kvector-fill code works properly
		fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);
		sprintf(char_buf0, "Checksum1 = %s, Checksum2 = %s\n", &char_buf1[convert_uint64_base16_char (char_buf1, checksum1)], &char_buf2[convert_uint64_base16_char (char_buf2, checksum2)]);
		fprintf(fp, "%s", char_buf0);	fprintf(fq, "%s", char_buf0);
	}

#ifndef FACTOR_STANDALONE
	fclose(fp);
#endif
	fp = 0x0;
	fclose(fq); fq = 0x0;

/*...init clocks, etc....*/
	clock1 = clock();
	tdiff = 0.0;

/* quick way to set all the bits = 1	*/
	for(i = 0; i < len; i++)
	{
		temp_late[i] = ~(uint64)0;
	}

/*...now generate q = 2kp+1 for as many k as desired; trial divide only if q passes
     a small-primes sieve and if q mod 8 = +1 or -1...	*/

	/*  q mod 8 = +-1 sieve is simplest, so do it first. Note q mod 8 = +-1 is guaranteed
	!   for p mod 60 sieve via choice of acceptable values of increment (incr), but doing
	!   q%8 = +-1 here is trivial and speeds search for multiples of small primes below.
	*/
	qmod8 = 1;
	for(i = 0; i < 4; i++)
	{
		/* Remember that k = 1 is in the zeroth bit here, i.e. we cycle through
		bits 0-3 which correspond to q = 2kp+1, 4kp+1, 6kp+1, 8kp+1, respectively.

		We don't need to worry about overflow-on-add of (two_p + qmod8),
		since an overflow won't affect the result, modulo 8.
		*/
		qmod8=(two_p[0] + qmod8) & 7;
		if(qmod8==3 || qmod8==5)
		{
			for(l = i; l < 64; l += 4)
			{
				temp_late[0] &= ~((uint64)1 << l); /* clear every fourth such bit from the first quadword...	*/
			}
		}
	}
/* ...Since the small-primes sieve needs only three (q mod 8 = +-1)-sieved registers to begin with, make 2 copies of the result...	*/
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
					if(curr_p < 100)
						printf("0: Found a multiple of %u in bit %u of register %u\n", curr_p, i, (uint32)k);
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
		/*...Now propagate copies of length (regs_todo) bit-cleared portion of sieve
		to remaining parts of sieve.
		*/
		if(m < (nclear-1))
		{
			ncopies = prime[m+1];
			l = regs_todo;
			for(i = 2; i <= ncopies; i++, l += regs_todo)
			{
				for(j = 0; j < regs_todo; j++)
				{
					temp_late[l+j] = temp_late[j];
				}
			}
		}
	}	/* endfor(m = 0; m < nclear; m++) */

#if FAC_DEBUG
	for(m = 0; m < len; m++)
	{
		for(i = 0; i < 64; i++)
		{
			on_bits += (temp_late[m]>>i) & 1;
		}
	}
	printf("%u ones bits of %u in sieve.\n", on_bits, len<<6);
#endif

	/*...calculate offsets, needed for primes larger than prime(nclear). */
	curr_p = p_last_small;
	for(m = nclear; m < nprime; m++)
	{
//	if((m & 0xffff) == 0) { printf("Offsets: m = %u of %u primes; curr_p = %u\n",(uint32)m,(uint32)nprime, curr_p); }
		curr_p += (pdiff[m] << 1);

		two64modp = 0x8000000000000000ull%curr_p;
		/*two64modp = (two64modp + two64modp) - curr_p;  Replace this with the following, which eschews the expensive library mod call: */
		two64modp = (two64modp + two64modp) - curr_p;
		two64modp += (-(two64modp >> 63)) & curr_p;

		/*q=1;*/

		/* If NUM_SIEVING_PRIME is such that the largest sieveprime >= Mersenne exponent p,
		without special handling here the sieve-clearing routine winds up sending a zero argument
		to the extended GCD routine used to calculate sieve bit offsets, causing an error.
		The fix is simple - in such cases we assign offset[m] the special value of 0xffffffff,
		and trap for that in the sieve-clearing routine - if this value is encountered for an
		offset[] element, the sieve-clearing routine skips over the prime in question (which
		must necessarily equal the exponent being tested) as it processes the sieving primes -
		this is of course OK, since p can never divide a factor candidate q = 2*k*p+1 anyway.
		*/
	#ifdef P1WORD
		twop_mod_currp = two_p[0]%curr_p;
		if(twop_mod_currp == 0)
		{
			ASSERT(HERE, p[0] == curr_p ,"2*p == 0 mod curr_p but p != curr_p");
			offset[m] = 0xffffffff;
			continue;
		}
	#else
		// 26. Sep 2012: This step is horribly slow for larger MMp - Accelerate by using that 2p = 2*Mp for double Mersennes.
		// Ex: For p = 43112609, the binary-powering-based version is ~10000x faster than the long-div-based one:
		if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS) {
			twop_mod_currp = twompmodq32(findex, curr_p);			// 2^-p (mod q)
			twop_mod_currp = modinv32(twop_mod_currp, curr_p)-1;	// 2^+p (mod q) = Mp (mod q)
			// modinv32 returns a signed result, do a conditional add to make >= 0
			twop_mod_currp += (-((int32)twop_mod_currp < 0)) & curr_p;
			twop_mod_currp += twop_mod_currp;
			if(twop_mod_currp >= curr_p) { twop_mod_currp -= curr_p; }
		//	ASSERT(HERE, twop_mod_currp == mi64_div_y32(two_p, curr_p, 0x0, lenP), "Fast 2p (mod q) for MMp fails!");
		} else {
			twop_mod_currp = mi64_div_y32(two_p, curr_p, 0x0, lenP);
		}
		if(twop_mod_currp == 0)
		{
//	t256 = NIL256;	t256.d0 = curr_p;
			if(0)	//==================== To-do: This check only valid for p prime ====================
			{
				ASSERT(HERE, mi64_cmp_eq_scalar(p, curr_p, lenP),"2*p == 0 mod curr_p but p != curr_p");
			}
			offset[m] = 0xffffffff;
			continue;
		}
	#endif

	#if 0	/* Change this to 1 if you prefer to do things really slowly */
		uint32 q_mod_currp = 1;

		for(k = 0; k < len; k++)
		{
			for(i = 0; i < 64; i++)
			{
			/*
				q += two_p;
				q_mod_currp = q%curr_p;	// This system mod call is *really* slow, so replace it with a better sequence:
			*/
				q_mod_currp  =  q_mod_currp + twop_mod_currp - curr_p;	// Subtract off curr_p...
				q_mod_currp += (-(q_mod_currp >> 31)) & curr_p;		// ...and restore it if result was < 0.

				if(q_mod_currp == 0)
				{
				#if FAC_DEBUG
					/* Print diagnostic info for selected prime subranges: */
					if(curr_p < 300)
						printf("A: Found a multiple of %8u in bit %6u; modinv = %8u\n", curr_p, (k<<6)+i+1, -modinv32((uint32)twop_mod_currp, curr_p)-1);

					if((k<<6)+i != -modinv32((uint32)twop_mod_currp, curr_p)-1)
					{
						fprintf(stderr,"MISMATCH: %u != %u\n", (k<<6)+i, -modinv32((uint32)twop_mod_currp, curr_p)-1);
						ASSERT(HERE, 0,"0");
					}
				#endif

					/*...later, will move startvalue into upper 2 bytes of prime(m).	*/
					offset[m] = (k<<6)+i;
					goto OUTER;
				}
				/* If find a k such that 2*k*p == +1 (rather than -1) mod curr_p, i.e. q == 2 mod curr_p,
				then this is the mirror image (w.r.to the sieve) of the bit we want:
				*/
				else if(q_mod_currp == 2)
				{
				#if FAC_DEBUG
					/* Print diagnostic info for the smaller primes: */
					if(curr_p < 300)
						printf("B: Found a multiple of %8u in bit %6u; modinv = %8u\n", curr_p, (k<<6)+i+1, +modinv32((uint32)twop_mod_currp, curr_p)-1);

					if((k<<6)+i != +modinv32((uint32)twop_mod_currp, curr_p)-1)
					{
						fprintf(stderr,"MISMATCH: %u != %u\n", (k<<6)+i, +modinv32((uint32)twop_mod_currp, curr_p)-1);
						ASSERT(HERE, 0,"0");
					}
				#endif

					/* Since bit 0 of sieve corresponds to factor index k = 1, do the mirror-image
					index subtraction in unit-offset mode, i.e. if mirror image is in (unit-offset)
					bit B', then B = curr_p - B' = curr_p - [(k<<6)+i+1].
					Subtract one from this to get B in zero-offset form:
					*/
					offset[m] = curr_p - ((k<<6)+i+2);
					goto OUTER;
				}
			}
		}
		fprintf(stderr,"ERROR: failed to find a multiple of prime %u\n", curr_p);
		ASSERT(HERE, 0,"0");

		OUTER: continue;

	#else
		/*
		Q: can we find these sieve bit offsets in an a priori fashion, rather than via
		expensive trial and error?
		A: Indeed we can. For example, consider p = 2^89 + 89 = 618970019642690137449562201.
		A snippet of output from the trial-and-error small-prime disivor code:

		...
		Found a multiple of 224699 in bit  6 of register  254	<== k = (64*254 + 7) = 16263
		Found a multiple of 224711 in bit 51 of register 3156
		Found a multiple of 224717 in bit 51 of register 1454
		Found a multiple of 224729 in bit 10 of register 2536
		Found a multiple of 224737 in bit 33 of register 2270
		Found a multiple of 224743 in bit 26 of register 1669
		...

		Each bit of the initial sieve (let's number them from 0 to n-1) represents an
		increment of 2*p. Need to find a bit k such that q(k) := 2*p*k + 1 == 0 mod curr_p.
		If we precompute m = 2*p mod curr_p, then q(k) mod curr_p = (m*k + 1) mod curr_p.
		In our example, the prime curr_p = 224699 has 2*p == 198862 mod curr_p. We need to
		find the smallest integer k such that k*198862 == -1 mod 224699. This gives k=16263.
		So what we need is basically a fast way to find the mod-inverse of the
		2*p (mod curr_p), and then simply negate that - We can use the eGCD for this.

		Example PARI code to find inverse of 198862 mod 224699
		(negate result to get the desired sign, +16263):

		x = 224699; y = 198862;

		a = 1; b = 0; g = x; u = 0; v = 1; w = y; w

		q = (g-g%w)/w;
		d = a - q*u;
		e = b - q*v;
		f = g - q*w;
		a = u;
		b = v;
		g = w;
		u=d;u
		v=e;v
		w=f;w

		(repeat above 10-step sequence until w=0).

		The sequence of w's is 198862, 25837, 18003, 7834, 2335, 829, 677, 152, 69, 14, 13, 1, 0.
		Inverse (value of b on exit from the above loop) is -16263.

		Note that being able to quit as soon as we find either k or k' halves the time
		needed to set up the sieve.
		More importantly, it reveals one subtlety in the eGCD approach: If the inverse
		of (2*p) mod curr_p is negative, then the inverse = -k and we simply negate it
		to get k. If on the other hand the inverse is positive, inverse == k' (mod curr_p),
		and we get k = curr_p - k'.

		Lastly, in the actual implementation we add one to k to get the location of the
		sieve bit, since the zero bit of our sieve corresponds to k = 1, not k = 0.
		*/
		if((int)twop_mod_currp < 0)
		{
			fprintf(stderr,"twop_mod_currp = %u for curr_p = %u out of range!\n",twop_mod_currp, curr_p);
			ASSERT(HERE, 0,"0");
		}

		i = modinv32((uint32)twop_mod_currp, curr_p);
		ASSERT(HERE, i, "factor.c : i");

		/* i < 0 corresponds to the (q_mod_currp == 0) case in the above (#ifdef'ed out) section.
		   i > 0 corresponds to the (q_mod_currp == 2) case.
		*/
		if((int)i < 0)
		{
		#if DBG_SIEVE
			/* Print diagnostic info for selected prime subranges: */
		//	if(curr_p > 116965000 && curr_p < 116966000)
		//		printf("C: Found a multiple of %8u in bit %6u\n", curr_p, -i-1);
		#endif
			offset[m] =        - (i+1);
		}
		else
		{
		#if DBG_SIEVE
			/* Print diagnostic info for selected prime subranges: */
		//	if(curr_p > 116965000 && curr_p < 116966000)
		//		printf("D: Found a multiple of %8u in bit %6u\n", curr_p, curr_p-i-1);
		#endif

			offset[m] = curr_p - (i+1);
		}

	#endif	/* end(#if 0) */
	}	/* endfor(m = nclear; m < nprime; m++) */

#ifdef FACTOR_STANDALONE
	 printf(   "TRYQ = %u, max sieving prime = %u\n",TRYQ,MAX_SIEVING_PRIME);
#else
	ASSERT(HERE, fp == 0x0,"0");
	fp = mlucas_fopen(STATFILE,"a");
	fprintf(fp,"TRYQ = %u, max sieving prime = %u\n",TRYQ,MAX_SIEVING_PRIME);
	fclose(fp); fp = 0x0;
#endif

	/*...create bitmap in atlas for each of the sixteen K mod 60 cases...	*/
	for(i = 0; i <= len/60; i++)
	{
		for(j = 0; j < 16; j++)
		{
			/*bit_atlas[i][j]=0;*/
			*(bit_atlas+(i<<4)+j)=0;	/* calloc! */

		}
	}

	i=incr[0]-1;
	j=0;
	k=0;
	for(;;) /* K loops over 64-bit registers...	*/
	{
		l=0;
L3:		if(k==len)break;
		for(;;) /* I loops over bits...	*/
		{
			word=(uint32)(j>>6);
			bit =(uint32)(j&63);
/*
if(l == 10 && word==1024 && bit == 24)
{
if(!((temp_late[k]>>i) & 1))printf("critical bit = 0! Loc in temp_late: word %d, bit %3u\n", k, i);
}
*/
			if((temp_late[k]>>i) & 1)
			{
				/*bit_atlas[word][l] |= ((uint64)1 << bit);*/
				*(bit_atlas+(word<<4)+l) |= ((uint64)1 << bit);
			}
			l++;
			if(l>>4)
			{
				l &= 15;
				j++;
			}
			i=i+incr[l];
			if(i>>6)
			{
				i &= 63;
				k++;
				goto L3;
			}
		} /* end of I loop	*/
	}	/* end of K loop	*/

/*...deallocate full-sized bit_atlas.	*/
	free((void *)temp_late); temp_late = 0x0;

/*...At this point, replace the relative with the absolute increments:	*/
	for(i = 1; i < 16; i++)	/* Skip pass 0 here */
	{
		incr[i] = incr[i-1] + incr[i];
	}

	i = 0;
	switch(pmod60)
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
			fprintf(stderr,"%u not an acceptable value for P mod 60. P likely nonprime.\n",pmod60); ASSERT(HERE, 0,"0");
	}

	clock2 = clock();	/* Assume sieve setup time < 2^32 cycles - even if that
							that is violated it's no big deal at this point. */
	/* Use td here, as tdiff is reserved for the total runtime from factoring start: */
	td = (double)(clock2 - clock1);;
	clock1 = clock2;

#ifdef FACTOR_STANDALONE
	if(!restart)
		printf("Time to set up sieve =%s\n",get_time_str(td));
#endif
	td = 0;	/* Use td to accumulate time spent on each individual factoring pass */

/* Run through each of the 16 "sievelets" as many times as necessary, each time copying
the appropriate q mod 8 and small-prime bit-cleared bit_atlas into memory, clearing bits
corresponding to multiples of the larger tabulated primes, and trial-factoring any
candidate factors that survive sieving.	*/

#if(TRYQ > 0)
	q_index = -1;	/* Make this < 0 if factor queue is empty. */
#endif

	nfactor = 0;

#if FAC_DEBUG
	/* If a known factor given, only process the given k/log2 range for that pass: */
	if(pass_targ < 16)
		passmin = passnow = passmax = pass_targ;
#endif

	for(pass = passnow; pass <= passmax; pass++)
	{
	#ifdef FACTOR_STANDALONE
		if(!restart)
		{
			printf("pass = %d",pass);
			fflush(stdout);
		}
	#else
		ASSERT(HERE, fp == 0x0,"0");
		fp = mlucas_fopen(STATFILE,"a");
		fprintf(fp,"Starting Trial-factoring Pass %2u...\n",pass);
		fclose(fp); fp = 0x0;
	#endif

		/* Starting no.-of-times-through-sieve = kmin/(64*len) : */
		if(pass == passnow && (know > kmin))
		{
			interval_lo = (know>>6)/(uint64)len;
			ASSERT(HERE, know == interval_lo *(len<<6),"know == interval_lo *(len<<6)");
		}
		else
		{
			interval_lo = (kmin>>6)/(uint64)len;
			ASSERT(HERE, kmin == interval_lo *(len<<6),"kmin == interval_lo *(len<<6)");
		}

		/* Set initial k for this pass to default value (= incr[pass]) + interval_lo*(64*len),
		(assume this could be as large as 64 bits), then use it to set initial q for this pass:
		*/
		ASSERT(HERE, (double)interval_lo*(len<<6) < TWO64FLOAT, "(double)interval_lo*len < TWO64FLOAT");
		k = (uint64)incr[pass] + interval_lo*(len<<6);
	  #if FAC_DEBUG
		printf(" Initial k for this pass = %llu.\n", k);
	  #endif

		/* 2.k.p: */
		ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,q,lenQ), "2.k.p overflows!");
		q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
	  #if FAC_DEBUG
		printf(" Initial q for this pass = %s.\n", &char_buf0[convert_mi64_base10_char(char_buf0, q, lenQ, 0)]);
	  #endif

		/* startbit k (occurrence of first multiple of prime curr_p in first pass
		through the relevant sievelet) is defined by

			(offset[curr_ p] + k*prime) (mod 60) = incr(pass)-1, k=0,59 .
		*/
		curr_p = p_last_small;
		for(m = nclear; m < nprime; m++)
		{
			curr_p += (pdiff[m] << 1);
			currp_mod_60 = curr_p%60;
	/********** TO DO: COULD DO A SEPARATE FIND-MULTIPLE OF nTH PRIME STEP FOR EACH PASS HERE, OBVIATING NEED FOR OFFSET[] ARRAY*****/
			l=offset[m];

			/* Special-handling code for p == curr_p case: */
			if(l == 0xffffffff)
			{
				startval[m] = 0xffffffff;
				continue;
			}

			lmod60 = l%60;

			for(i = 0; i <= 60; i++)
			{
				/*lmod60 = l%60;	// This system mod call is *really* slow, so replace it with the better sequence below: */

				if(lmod60 == incr[pass]-1)
				{
					l=l/60;
					goto SET_START;
				}
				l += curr_p;

				lmod60 += currp_mod_60 - 60;		/* Subtract off 60... */
				lmod60 += (-(lmod60 >> 31)) & 60;	/* ...and restore it if result was < 0. */
			}
			fprintf(stderr,"ERROR: failed to find a multiple of prime %u\n", curr_p);
			ASSERT(HERE, 0,"0");

		#if FAC_DEBUG
			if(l >= len)
				fprintf(stderr,"WARNING: first multiple of prime %u lies in bit %u, outside sievelet length %u\n", curr_p, l, len);
		#endif

		SET_START: startval[m] = l;

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

			If we want to accommodate arbitrarily large kmin values for the start
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

			/* bit_len is a uint32, so use i (also a 32-bit) in place of k (64-bit) here: */
			i = ceil(1.0*bit_len/curr_p);
			ASSERT(HERE, i*curr_p - bit_len == curr_p - (bit_len % curr_p), "i*curr_p - bit_len == curr_p - (bit_len % curr_p)");

			/* Now calculate dstartval for the actual current-pass kmin value,
			according to the number of times we'd need to run through the sieve
			(starting with k = 0) to get to kmin: */
			dstartval = (uint64)(i*curr_p - bit_len);
			dstartval = (interval_lo*dstartval)%curr_p;
			dstartval += startval[m];
			if(dstartval >= curr_p)
				startval[m] = dstartval - curr_p;
			else
				startval[m] = dstartval;

		#if FAC_DEBUG
			ASSERT(HERE, startval     [m] < curr_p, "factor.c : startval     [m] < curr_p");
		  #if DBG_SIEVE
			startval_incr[m] = i*curr_p - bit_len;
			ASSERT(HERE, startval_incr[m] < curr_p, "factor.c : startval_incr[m] < curr_p");
		  #endif
		#endif

		}	/* endfor(m = nclear; m < nprime; m++) */

		for(sweep = interval_lo; sweep < interval_hi; ++sweep)
		{
		/*
			printf("sweep = %8d  kstart = %s, k_targ = %s\n", sweep, k, k_targ);
		*/
			// To track lg(q) = lg(2.k.p+1), use approximation q ~= 2.k.p, thus lg(q) ~= lg(2.p) + lg(k).
			// At start of each pass through the k-based sieve, use 2.k.p with k = [starting k + sievebits]
			// to bound qmax from above, and compute lg(qmax) using the above logarithmic sum.
			fbits_in_k = log((double)k + 60*bit_len)*ILG2;	// Use k-value at end of upcoming pass through sieve as upper-bound
			fbits_in_q = fbits_in_2p + fbits_in_k;
//printf("sweep = %d: fbits_in_q = fbits_in_2p [%10.4f] + fbits_in_k [%10.4f] = %10.4f\n",sweep,fbits_in_2p,fbits_in_k,fbits_in_q);
			/* Accumulate the cycle count every so often to avoid problems with integer overflow
			of the clock() result, if clock_t happens to be a 32-bit int type on the host platform:
			*/
			if((sweep & 127) == 0)
			{
				clock2 = clock();
				td += (double)(clock2 - clock1);
				clock1 = clock2;

			#ifdef FACTOR_STANDALONE
				if(!restart)
				{
					printf(".");
					fflush(stdout);
				}
			#endif
			}

		/*********************************************/
		#if DBG_SIEVE
			if(k_targ)
			{
				/* See if k_targ falls in the range of k's for the current sieve interval: */
				k = (uint64)incr[pass] + sweep*(len<<6);	/* Starting k for the current interval: */

				/* If so, calculate the location of the critical bit
				and allow execution to proceed as normal to bit-clearing step:
				*/
				if((k <= k_targ) && (k_targ < (k+(len<<6))))
				{
					itmp64 = k_targ - k;
					ASSERT(HERE, itmp64%60 == 0,"(k_targ - k)%60 == 0");
					itmp64 /= 60;
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
						if(startval[m] >= bit_len)
						{
							startval[m] -= bit_len;
						}
						else
						{
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

			/*   load the sievelet...	*/
			for(i = 0; i <= len/60; i++)
			{
				/*bit_map[i] = bit_atlas[i][pass];	// Use memcpy here? */
				bit_map[i] = *(bit_atlas+(i<<4)+pass);	/* Use memcpy here? */
			}

		#if DBG_SIEVE
			/* If debugging sieve, make sure critical bit hasn't been cleared: */
			if( k_targ && (((bit_map[i64_targ] >> bit_targ) & 1) == 0) )
			{
				fprintf(stderr,"Critical bit cleared in original bitmap!\n");
				ASSERT(HERE, 0,"0");
			}
		#endif
			/*   ...and clear the bits corresponding to the small primes.	*/
			curr_p = p_last_small;
			for(m = nclear; m < nprime; m++)
			{
				curr_p += (pdiff[m] << 1);
				/*
				if(current_prime == 107) printf("   prime %8d has offset = %8d\n", curr_p, startval[m]);
		`		*/

				l = startval[m];

				/* Special-handling code for p == curr_p case: */
				if(l == 0xffffffff)
				{
					continue;
				}

				/* Need to account for the fact that primes greater than the # of bits in the sieve
				may not hit *any* sieve bit on the current pass through the sieve:
				*/
				if(l >= bit_len)
				{
					startval[m] -= bit_len;
				}
				else
				{
					while(l < bit_len)
					{
						i64=l >> 6;	/* l/64	*/
						bit=l & 63;	/* l%64	*/
					#if DBG_SIEVE
						if(k_targ && (i64 == i64_targ && bit == bit_targ))
						{
							fprintf(stderr,"Critical bit being cleared by prime %u, with offset %u\n", curr_p, startval[m]);
							ASSERT(HERE, 0,"0");
						}
					#endif
						mask = ~((uint64)1 << bit);
						bit_map[i64] &= mask;

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

			/*...now run through the bits of the current copy of the sieve,
			trial dividing if a bit = 1.
			*/
			ihi = ((bit_len+63)>>6);
			bit_hi = 64;
			for(i = 0; i < ihi; i++)	/* K loops over 64-bit registers. Don't assume bit_len a multiple of 64.	*/
			{
				/* Special casing for last bit_atlas word, which may be only partly full */
				if((bit_len - (i<<6)) < 64)
					bit_hi = (bit_len - (i<<6));

				for(bit = 0; bit < bit_hi; bit++)		/* BIT loops over bits in each bit_atlas word...	*/
				{
				#if(FAC_DEBUG)
					/* If a known factor is specified, here it is in the bitmap: */
					if(k == k_targ)
					{
						printf("here it is: sweep = %s, bitmap word = %u, bit = %3u\n", &char_buf0[convert_uint64_base10_char(char_buf0, sweep)], i, bit);

						if((bit_map[i] >> bit) & 1)
							printf("Trying k_targ = %llu...\n", k_targ);
						else
							ASSERT(HERE, 0,"0");
					}
				#endif

				/*...If current sieve bit=1, add q to trial-division queue: */

					if((bit_map[i] >> bit) & 1)
					{

		  #ifdef USE_GPU	// Set nonzero to print each factor candidate tried

			kvidx = (count & kvm1);	// kvsz a power of 2, so %kvsz can be effected via &(kvsz-1)

			// If count a multiple of kvsz and (count > 0), this means we have "rolled over the odometer:
			// Test all the k's in the now-full k-vector, and start accumulating them anew:
			if(kvidx == 0)
			{
			GPU_FACTOR_BATCH:
				if(count == 0) {
					// Compute auxiliary TF data:
					ASSERT(HERE, (*p >> 32) == 0, "p must be < 2^32 for GPU TF!");
					pshift = (uint32)*p + 78;
					jshift = leadz32(pshift);
					/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
					lead7 = ((pshift<<jshift) >> 25);
					if(lead7 > 77) {
						lead7 >>= 1;
						start_index =  32-jshift-6;	/* Use only the leftmost 6 bits */
					} else {
						start_index =  32-jshift-7;
					}
					zshift = 77 - lead7;
					zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
					pshift = ~pshift;
				} else {
					if(kvidx != 0)
					{	// Adjust array upper bound if it's a final partially-full cleanup-batch of candidates
						kv_numq = kvidx;
					} else {			// Regular complete batch
						kv_numq = kvsz;
					}
				printf("count = %15llu * 2^20: Sending range k = [%20llu, %20llu] to GPU for testing...\n", (count >> 20), cpu_kvec[0],cpu_kvec[kv_numq-1]);

					//=========== to-do: Add factor-#words-wrapper here ===========

					// Since GPU code tests many candidates per batch, adopt the following return-value convention:
					//			|	= 0		No factors found in current batch
					//	retval	|	> 0		 1 factor in current batch, whose k-value is returned
					//			|	< 0		n = (-retval) (>1) factors in current batch, user should re-sieve interval in 'slow'
					//									CPU-based mode to find the factor k's
					//
					// Copy kvec from host to device memory and Invoke GPU-parallel kernel:
					threadsPerBlock = 1;	//==== Need to init from global set in util.c system-param-init ====
					blocksPerGrid = ((kv_numq/TRYQ) + threadsPerBlock - 1) / threadsPerBlock;

				fprintf(stderr,"cudaMemcpy into address 0x%16x...",gpu_kvec);
					cudaMemcpy(gpu_kvec, cpu_kvec, kvsz, cudaMemcpyHostToDevice);
				fprintf(stderr,"done.");
				//	for(int ki=0; ki<kvsz; ++ki) {
				//		ASSERT(HERE, gpu_kvec[ki] == cpu_kvec[ki], "cudaMemcpy cpu_kvec --> gpu_kvec failed!");
				//	}
				fprintf(stderr,"\nGPU_TF78<<< %d, %d >>>...",blocksPerGrid, threadsPerBlock);
					GPU_TF78<<<blocksPerGrid, threadsPerBlock>>>((uint32)*p, pshift, zshift, start_index, gpu_kvec, gpu_rvec, kv_numq);
				//	cudaThreadSynchronize();

					// Copy bytewise output vector from the device to host:
					cudaMemcpy(cpu_rvec, gpu_rvec, kvsz, cudaMemcpyDeviceToHost);
					for(m = 0; m < kv_numq; m++) {
						if(cpu_rvec[m] != 0)
						{
							factor_k[nfactor++] = cpu_kvec[m];
							sprintf(char_buf0,"\n\trvec[%d] = %d: Factor with k = %s. Program: E%s\n",m,cpu_rvec[m], &char_buf1[convert_uint64_base10_char(char_buf1, factor_k[nfactor-1])], VERSION);
							fprintf(stderr,"%s", char_buf0);
				exit(0);
						}
					}
					//=== cudaFree(gpu_kvec) ? ===

					// Adjust array upper bound if it's a final partially-full cleanup-batch of candidates
					if(kvidx != 0)
					{
						goto GPU_CLEANUP_DONE;
					}
					kvidx = 0;
				}
			}
			cpu_kvec[kvidx] = k;
			count++;
			goto INCR_K;
		  #endif

						q_index = count++ & TRYQM1;	/* Post-increment count, so this will work. */

						/* Every so often (every [2^32 / (nearest power of 2 to lenP*lenQ^2)] q's seems a good interval)
						do some factor-candidate sanity checks.

						***NOTE: *** The (count & ...) here must be with a constant = 2^n-1, n >= 3.
						*/
// 25. Sep 2012: Disable the checkpointing (and anything else which attempts to actually compute something mod-q) for large-MMp deep-candidate-sieve code-dev:
if(0)
//						if((count & countmask) == 0)
						{
							fp = mlucas_fopen(OFILE,"a");

							ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,q,lenQ), "2.k.p overflows!");
							q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
						#if(FAC_DEBUG)
							if(!restart) {
								sprintf(cbuf, " Count = %u * 2^%u: k = %llu, Current q = %s, checksum1 = %llu\n",
									(uint32)(count >> CMASKBITS),CMASKBITS,k,&char_buf1[convert_mi64_base10_char(char_buf1, q, lenQ, 0)],checksum1);
								fprintf(stderr, "%s", cbuf);
							}
						#endif

							/* Do some sanity checks on the occasional factor candidate to ensure
							that the sieve is functioning properly and no overflows have occurred.
							*/

							/* In MMp mode, make sure results of fast and slow-debug modpow agree: */
							if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS) {
printf("Exiting at first checkpoint.\n");
exit(0);
								res = (mi64_twopmodq_qmmp(findex, k, u64_arr) == 1);
								if(res != (mi64_twopmodq(p, lenP, k, q, lenQ, q2) == 1) || q2[0] != u64_arr[0]) {
									sprintf(cbuf, "ERROR: Spot-check k = %llu, Results of mi64_twopmodq_qmmp and mi64_twopmodq differ!\n", k);
									fprintf(fp,"%s", cbuf);
									if(!restart) {
										ASSERT(HERE, 0, cbuf);
									}
								}
							}

							/* Make sure that q == 1 (mod 2p) and that q is a base-2 PRP: */
							mi64_clear(u64_arr, lenQ);	// Use q2 for quotient [i.e. factor-candidate k] and u64_arr for remainder
							mi64_div(q,two_p,lenQ,lenQ,q2,u64_arr);
							if(mi64_getlen(q2, lenQ) != 1) {
								sprintf(cbuf, "ERROR: Count = %u * 2^%u: k = %llu, Current q = %s: k must be 64-bit!\n",
									(uint32)(count >> CMASKBITS),CMASKBITS,k,&char_buf1[convert_mi64_base10_char(char_buf1, q, lenQ, 0)]);
									fprintf(fp,"%s", cbuf);
									if(!restart) {
										ASSERT(HERE, 0, cbuf);
									}
							}
							if(!mi64_cmp_eq_scalar(u64_arr, 1ull, lenQ))
							{
								sprintf(cbuf, "ERROR: Count = %u * 2^%u: k = %llu, Current q = %s: q mod (2p) = %s != 1!\n",
									(uint32)(count >> CMASKBITS),CMASKBITS,k,&char_buf1[convert_mi64_base10_char(char_buf1, q, lenQ, 0)],
									&char_buf2[convert_mi64_base10_char(char_buf2, u64_arr, lenQ, 0)]);
									fprintf(fp,"%s", cbuf);
									if(!restart) {
										ASSERT(HERE, 0, cbuf);
									}
							}

							/* If q is composite [only check this in debug mode since it costs more than checking
							divisibility of 2^p-1 by q], make sure its smallest divisor
							is larger than the max sieving prime being used for the run: */
							mi64_set_eq    (q2, q, lenQ);
							mi64_sub_scalar(q2,1ull,q2,lenQ);	// Re-use q2 to store q-1
							if(mi64_twopmodq(q2, lenQ, 0, q, lenQ, 0x0) != 1)
							{
								if(!restart) {
									printf(" INFO: Spot-check q with k = %llu is composite\n",k);
								}

								l = 3;
								idx = 0;
								for(m = 0; m < nprime; m++)
								{
									/* Is q%(current small sieving prime) == 0? */
									if(mi64_div_y32(q, l, 0x0, lenQ) == 0)
									{
										sprintf(cbuf, "ERROR: Count = %u * 2^%u: k = %llu, Current q = %s has a small divisor: %u\n",
											(uint32)(count >> CMASKBITS),CMASKBITS,k,&char_buf1[convert_mi64_base10_char(char_buf1, q, lenQ, 0)],l);
										fprintf(fp,"%s", cbuf);
										if(!restart) {
											ASSERT(HERE, 0, cbuf);
										}
									}
									l += (pdiff[++idx] << 1);
								}
							} else {
								if(!restart) {
									printf(" INFO: Spot-check q with k = %llu is base-2 PRP\n",k);
								}
							}
							fclose(fp); fp = 0x0;

						}	/* endif((count & countmask) == 0) */

					#if(TRYQ == 0)	/* If testing speed of sieve alone, skip to incrementing of q. */

					#else

						k_to_try[q_index] = k;

						if(q_index == TRYQM1)
						{
							q_index = -1;	/* q_index = -1 indicates factor-candidate queue empty. (More precisely,
											it will be, after we test the current batch of candidates. */

						#if(TRYQ == 1)	/* try 1 factor candidate at a time */

						  #if(defined(NWORD))

							if(MODULUS_TYPE == MODULUS_TYPE_MERSMERS) {
								if(findex > 10000) {
									if(k < 1000) {
										kmod60 = k%60;
									//	printf("pmod60, kmod60 = %u, %u\n",pmod60,kmod60);
										ASSERT(HERE, CHECK_PKMOD60(pmod60, kmod60), "K must be one of the admissible values (mod 60) for this p (mod 60)!");
										printf("Do deep sieving for k = %u\n",(uint32)k);
										kdeep[ndeep++] = (uint32)k;
										ASSERT(HERE, ndeep < 1024, "Increase allocation of kdeep[] array or use deeper sieving bound to reduce #candidate k's!");
									//	itmp64 = factor_qmmp_sieve64((uint32)findex, k, MAX_SIEVING_PRIME+2, 0x0001000000000000ull);
									//	if(itmp64) {
									//		printf("Q( k = %u ) has a small factor: %20llu\n",(uint32)k, itmp64);
									//	}
									}
									res = 0;
								} else {
									ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,q,lenQ), "2.k.p overflows!");
									q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
									checksum1 += q[0];
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
								checksum1 += q[0];
								res = mi64_twopmodq(p, lenP, k, q, lenQ, u64_arr);
							}
							checksum2 += u64_arr[0];

						  #elif(defined(P4WORD))

							ASSERT(HERE, p[3] == 0 && 0 == mi64_mul_scalar(two_p,k,(uint64*)&q256,lenQ), "2.k.p overflows!");
							q256.d0 += 1;	// No need to check for carry since 2.k.p even
							p256.d0 = p[0]; p256.d1 = p[1]; p256.d1 = p[1]; p256.d2 = p[2]; p256.d3 = p[3];
							checksum1 += q256.d0;
							t256 = twopmodq256(p256,q256);
							checksum2 += t256.d0;
							res = CMPEQ256(t256, ONE256);

						  #elif(defined(P3WORD))

						   #if(defined(USE_FLOAT))

						    ASSERT(HERE, !p[2], "twopmodq200: p[2] nonzero!");
							x256 = twopmodq200_8WORD_qmmp(&chk1, &chk2, p,k);	res = (uint64)CMPEQ256(x256, ONE256);
							checksum1 += chk1;	checksum2 += chk2;

						   #if 0	/*** Set == 1 to enable check-vs-192-bit ***/
							/* Compare to strictly 192-bit-int routines: */
							q192.d1 = q192.d2 = 0ull;
							ASSERT(HERE, p[2] == 0 && 0 == mi64_mul_scalar(two_p,k,(uint64*)&q192,lenQ), "2.k.p overflows!");
							q192.d0 += 1;	// No need to check for carry since 2.k.p even
							p192.d0 = p[0]; p192.d1 = p[1]; p192.d1 = p[1]; p192.d2 = p[2];
							t192 = twopmodq192(p192,q192);
							// No need to cast x256 to uint192 for compare since compare macro simply examines 64-bit words 0,1,2:
							if(count > 0 && (!CMPEQ192(x256, t192) || x256.d3 != 0)) {
								fprintf(stderr, "200-bit-float result [%2d] vs 192-bit-int [%2] mismatch: count = %llu\n", count);
								ASSERT(HERE, 0, "abort");
							}
						   #endif	/* if(0) */

						   #else

							ASSERT(HERE, p[2] == 0 && 0 == mi64_mul_scalar(two_p,k,(uint64*)&q192,lenQ), "2.k.p overflows!");
							q192.d0 += 1;	// No need to check for carry since 2.k.p even
							p192.d0 = p[0]; p192.d1 = p[1]; p192.d1 = p[1]; p192.d2 = p[2];
							checksum1 += q192.d0;
							t192 = twopmodq192(p192,q192);
							checksum2 += t192.d0;
							res = CMPEQ192(t192, ONE192);

						   #endif	/* USE_FLOAT */

						  #elif(defined(P2WORD))

						   #if(  defined(USE_128x96) && USE_128x96 == 1)
							/* Use strictly  96-bit routines: */
							if(p[1] == 0 && (q[1] >> 32) == 0)
								res = twopmodq96	(&checksum1, &checksum2, p, k);
							else
						   #elif(defined(USE_128x96) && USE_128x96 == 2)
							/* Use hybrid 128_96-bit routines: */
							if(p[1] == 0 && (q[1] >> 32) == 0)
								res = twopmodq128_96(&checksum1, &checksum2, p, k);
							else
						   #endif
							/* Use fully 128-bit routines: */
							res = twopmodq128x2(&checksum1, &checksum2, p, k);

						  #else

							#ifdef USE_FMADD
								/* Use 50x50-bit FMADD-based modmul routines, if def'd: */
								res = twopmodq100_2WORD_DOUBLE(&checksum1, &checksum2, p,k);
							#elif(defined(USE_FLOAT))
								/* Otherwise use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE(&checksum1, &checksum2, p,k);
							#else
								/* Otherwise use pure-integer-based modmul: */
								if(fbits_in_q < 63)
								{
									itmp64 = k*p;	itmp64 += itmp64 + 1;
									checksum1 += itmp64;
									res = twopmodq63(p,itmp64);
									checksum2 += res;
									res = (res == 1);
								}
								else if(fbits_in_q < 64)
								{
									itmp64 = k*p;	itmp64 += itmp64 + 1;
									checksum1 += itmp64;
									res = twopmodq64(p,itmp64);
									checksum2 += res;
									res = (res == 1);
								}
							  #if(defined(USE_65BIT))
								else if(fbits_in_q < 65)
								{
									res = twopmodq65(&checksum1, &checksum2, p,k);
								}
							  #endif
								else
								{
									ASSERT(HERE, fbits_in_q < 96, "fbits_in_q exceeds allowable limit of 96!");
								  #if(  defined(USE_128x96) && USE_128x96 == 1)
									/* Use strictly  96-bit routines: */
									res = twopmodq96	(&checksum1, &checksum2, p,k);
								  #elif(defined(USE_128x96) && USE_128x96 == 2)
									/* Use hybrid 128_96-bit routines: */
									res = twopmodq128_96(&checksum1, &checksum2, p,k);
								  #else
									/* Use fully 128-bit routines: */
									ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k,(uint64*)&q128,lenQ), "2.k.p overflows!");
									q128.d0 += 1;	// No need to check for carry since 2.k.p even
									p128.d0 = p[0]; p128.d1 = p[1];
									checksum1 += q128.d0;
									t128 = twopmodq128x2(p128,q128);
									checksum2 += t128.d0;
									res = CMPEQ128(t128, ONE128);
								  #endif
								}
							#endif	/* #ifdef USE_FMADD */
						  #endif	/* #if(defined(P3WORD)) */

						#elif(TRYQ == 2)	/* try 2 factor candidates at a time */

						  #if(defined(P3WORD))

							#ifndef USE_FLOAT
								#error	TRYQ = 2 / P3WORD only allowed if USE_FLOAT is defined!
							#endif	/* #ifdef USE_FMADD */

						    ASSERT(HERE, !p[2], "twopmodq200: p[2] nonzero!");
							res  = twopmodq200_8WORD_qmmp_x2_sse2(&checksum1, &checksum2, p,k_to_try[0],k_to_try[1]);

						  #elif(defined(P1WORD))

							#ifdef USE_FMADD
								/* Use 50x50-bit FMADD-based modmul routines, if def'd: */
								res = twopmodq100_2WORD_DOUBLE_q2(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1]);
							#elif(defined(USE_FLOAT))
								/* Otherwise use 78-bit floating-double-based modmul: */
								res = twopmodq78_3WORD_DOUBLE_q2(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1]);
							#else
								#error	TRYQ = 2 / P1WORD only allowed if USE_FLOAT or USE_FMADD is defined!
							#endif	/* #ifdef USE_FMADD */

						  #else

							#error TRYQ == 2 requires P1WORD or P3WORD to be defined!

						  #endif

						#elif(TRYQ == 4)	/* try 4 factor candidates at a time */

						  #if(defined(P3WORD))

						//	ASSERT(HERE, !p[2], "twopmodq200: p[2] nonzero!");
							res = twopmodq192_q4(&checksum1, &checksum2, p,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

						  #elif(defined(P2WORD))

							/* Itanium seems to perform poorly with 96-bit variant: */
						   #if(  defined(USE_128x96) && USE_128x96 == 1)
							/* Use strictly  96-bit routines: */
							if(p[1] == 0 && fbits_in_q < 96)
							{
								res = twopmodq96_q4     (&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
							}
							else
						   #elif(defined(USE_128x96) && USE_128x96 == 2)
							/* Use hybrid 128_96-bit routines: */
							if(p[1] == 0 && fbits_in_q < 96)
								res = twopmodq128_96_q4 (&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
							else
						   #endif
							/* Use fully 128-bit routines: */
								res = twopmodq128_q4    (&checksum1, &checksum2, p   ,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

						  #else

							#ifdef USE_FMADD
								/* Use 50x50-bit FMADD-based modmul routines, if def'd: */
								res = twopmodq100_2WORD_DOUBLE_q4(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);

							#elif(defined(USE_FLOAT))
								/* Otherwise use 78-bit floating-double-based modmul: */
								// Quick check as to whether this routine is the one being used:
							  #if 0
								chk1 = chk2 = 0;
								ASSERT(HERE, 0, "Using twopmodq78_3WORD_DOUBLE_q4\n");
							  #endif

							  #if 1	/*** Set == 0 to enable check-vs-96-bit ***/
								res = twopmodq78_3WORD_DOUBLE_q4(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
							  #else
								res = twopmodq78_3WORD_DOUBLE_q4(&chk1,&chk2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								checksum1 += chk1;
								checksum2 += chk2;
								/* Compare to strictly 96-bit-int routines: */
								chk3 = chk4 = 0;
								res = twopmodq96_q4(&chk3, &chk4, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								if(count > 0 && (chk1!=chk3 || chk2!=chk4))
								{
									fprintf(stderr, "78-bit-float vs 96-bit-int mismatch: count = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, count)] );
								}
							  #endif	/* if(0) */
							#else
								/* fbits_in_q is calculated for largest q of current
								set, so if that > 64, use 96-bit routines for all q's. */
								if(fbits_in_q < 63)
								{
									res = twopmodq63_q4(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								}
								else if(fbits_in_q < 64)
								{
									res = twopmodq64_q4(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								}
							  #if(defined(USE_65BIT))
								else if(fbits_in_q < 65)
								{
									res = twopmodq65_q4(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								}
							  #endif
								else
								{
									ASSERT(HERE, fbits_in_q < 96, "fbits_in_q exceeds allowable limit of 96!");
								#if(defined(USE_128x96) && USE_128x96 == 1)
									/* Use strictly  96-bit routines: */
									res = twopmodq96_q4		(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								#elif(!defined(USE_128x96) || USE_128x96 == 2)
									/* Use hybrid 128_96-bit routines: */
									res = twopmodq128_96_q4	(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								#else
									/* Use fully 128-bit routines: */
									res = twopmodq128_q4	(&checksum1, &checksum2, p   ,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3]);
								#endif
								}
							#endif	/* #ifdef USE_FLOAT */
						  #endif	/* #if(defined(P3WORD)) */

						#elif(TRYQ == 8)	/* try 8 factor candidates at a time */

						  #if(defined(P3WORD))
						/*
							if((q[2] >> 32) == 0)
								res = twopmodq160_q8(&checksum1, &checksum2, p,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
							else
						*/
								res = twopmodq192_q8(&checksum1, &checksum2, p,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);

						  #elif(defined(P2WORD))

							/* Itanium seems to perform poorly with 96-bit variant: */
						   #if(  defined(USE_128x96) && USE_128x96 == 1)
							/* Use strictly  96-bit routines: */
							if(p[1] == 0 && fbits_in_q < 96)
								res = twopmodq96_q8     (&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
							else
						   #elif(defined(USE_128x96) && USE_128x96 == 2)
							/* Use hybrid 128_96-bit routines: */
							if(p[1] == 0 && fbits_in_q < 96)
								res = twopmodq128_96_q8 (&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
							else
						   #endif
							/* Use fully 128-bit routines: */
								res = twopmodq128_q8    (&checksum1, &checksum2, p   ,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);

						  #else

							/* fbits_in_q is calculated for largest q of current
							set, so if that > 64, use 65,78 or 96-bit routines for all q's. */
							#if defined(USE_FLOAT) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
								/* Otherwise use 78-bit floating-double-based modmul: */
								chk1 = chk2 = 0;
							  #if 0
								ASSERT(HERE, 0, "Using twopmodq78_3WORD_DOUBLE_q8\n");
							  #endif
								res = twopmodq78_3WORD_DOUBLE_q8(&chk1, &chk2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								checksum1 += chk1;
								checksum2 += chk2;
							#else
								if(fbits_in_q < 63)
								{
									res = twopmodq63_q8(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								}
								else if(fbits_in_q < 64)
								{
									res = twopmodq64_q8(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								}
							  #if(defined(USE_65BIT))
								else if(fbits_in_q < 65)
								{
									res = twopmodq65_q8(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								}
							  #endif
								else
								{
									ASSERT(HERE, fbits_in_q < 96, "fbits_in_q exceeds allowable limit of 96!");
								#if(defined(USE_128x96) && USE_128x96 == 1)
									/* Use strictly  96-bit routines: */
									res = twopmodq96_q8		(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								#elif(!defined(USE_128x96) || USE_128x96 == 2)
									/* Use hybrid 128_96-bit routines: */
									res = twopmodq128_96_q8	(&checksum1, &checksum2, p[0],k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								#else
									/* Use fully 128-bit routines: */
									p128.d0 = p[0]; p128.d1 = 0;
									res = twopmodq128_q8	(&checksum1, &checksum2, p   ,k_to_try[0],k_to_try[1],k_to_try[2],k_to_try[3],k_to_try[4],k_to_try[5],k_to_try[6],k_to_try[7]);
								#endif
								}

							#endif	/* #ifdef USE_FLOAT */
						  #endif	/* #if(defined(P3WORD)) */

						#elif(TRYQ == 16)	/* try 16 factor candidates at a time */

							#if(defined(NWORD) || defined(P4WORD) || defined(P3WORD) || defined(P2WORD))
								#error (TRYQ == 16) only supported for 64-bit/P1WORD/GCC/AVX builds!
							#elif defined(USE_FLOAT) && defined(USE_AVX) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
								/* Use 78-bit floating-double-based modmul: */
								chk1 = chk2 = 0;
								res = twopmodq78_3WORD_DOUBLE_q16(&chk1, &chk2, p[0]
									,k_to_try[  0],k_to_try[  1],k_to_try[  2],k_to_try[  3],k_to_try[  4],k_to_try[  5],k_to_try[  6],k_to_try[  7]
									,k_to_try[0x8],k_to_try[0x9],k_to_try[0xa],k_to_try[0xb],k_to_try[0xc],k_to_try[0xd],k_to_try[0xe],k_to_try[0xf]
								);
								checksum1 += chk1;
								checksum2 += chk2;
							#else
								#error (TRYQ == 16) only supported for 64-bit/P1WORD/GCC/AVX builds!
							#endif	/* #ifdef USE_FLOAT */

						#endif	/* endif(TRYQ == ...) */

							/* Print any factors that were found in the current batch: */
							for(l = 0; l < TRYQ; l++)
							{
								if((res >> l) & 1)	/* If Lth bit = 1, Lth candidate of the inputs is a factor */
								{
									/* Recover the factor: */
									q[lenP] = mi64_mul_scalar( p, 2*k_to_try[l], q, lenP);
									q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
								  #if FAC_DEBUG
									printf("Factor found: q = %s. Checking primality...\n", &char_buf0[convert_mi64_base10_char(char_buf0, q, lenQ, 0)]);
								  #endif

									/* Do a quick base-3 compositeness check (base-2 would be much faster due to
									our fast Montgomery arithmetic-based powering for that, but it's useless for
									weeding out composite Mersenne factors since those are all base-2 Fermat pseudoprimes).
									If it's composite we skip it, since we expect to recover the individual prime subfactors
									on subsequent passes (although this should only ever happen for small p and q > (2p+1)^2 :
									*/
									if(mi64_pprimeF(q, 3ull, lenQ))
									{
										factor_k[nfactor++] = k_to_try[l];
										sprintf(char_buf0,"\n\tFactor with k = %llu. Program: E%s\n", k_to_try[l], VERSION);

									#ifdef FACTOR_STANDALONE
									  if(!restart)
										fprintf(stderr,"%s", char_buf0);
									#else
										ASSERT(HERE, fp == 0x0,"0");
										fp = mlucas_fopen(STATFILE,"a");
										fprintf(fp,"%s", char_buf0);
										fclose(fp); fp = 0x0;
									#endif

										fp = mlucas_fopen(   OFILE,"a");
										fprintf(fp,"%s", char_buf0);
										fclose(fp); fp = 0x0;

									#if FAC_DEBUG
										if(TRYQM1 > 1)
											printf("factor was number %u of 0-%u in current batch.\n", l, TRYQM1);
									#endif

									#ifdef QUIT_WHEN_FACTOR_FOUND
										return 0;
									#endif
									}	/* endif(factor a probable prime?) */
								}	/* endif((res >> l) & 1)		*/
							}	/* endfor(l = 0; l < TRYQ; l++)	*/
						}	/* endif(q_index == TRYQM1)		*/

					#endif	/* endif(TRYQ == 0) */

					}	/* endif((bit_map[i] >> bit) & 1)	*/
		#ifdef USE_GPU
			INCR_K:
		#endif
					// Increment k, i.e. increment current q by p120:
					k += 60;

				} /* end of BIT loop	*/
			}	/* end of K loop	*/

		#if(TRYQ > 1)

			/* Clean up any remaining in queue */
			if(q_index >= 0)
			{
				for(l = 0; l <= (uint32)q_index; l++)
				{
				#if FAC_DEBUG
					ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k_to_try[l],q,lenQ), "2.k.p overflows!");
					q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
					printf("A: Trying q = %s\n", &char_buf0[convert_mi64_base10_char(char_buf0, q, lenQ, 0)]);
				#endif

				#if(defined(P4WORD))

					ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k_to_try[l],(uint64*)&q256,lenQ), "2.k.p overflows!");
					q256.d0 += 1;	// No need to check for carry since 2.k.p even
					p256.d0 = p[0]; p256.d1 = p[1]; p256.d1 = p[1]; p256.d2 = p[2]; p256.d3 = p[3];
					checksum1 += q256.d0;
					t256 = twopmodq256(p256,q256);
					checksum2 += t256.d0;
					res = CMPEQ256(t256, ONE256);

				#elif(defined(P3WORD))

					ASSERT(HERE, 0 == mi64_mul_scalar(two_p,k_to_try[l],(uint64*)&q192,lenQ), "2.k.p overflows!");
					q192.d0 += 1;	// No need to check for carry since 2.k.p even
					p192.d0 = p[0]; p192.d1 = p[1]; p192.d1 = p[1]; p192.d2 = p[2];
					checksum1 += q192.d0;
					t192 = twopmodq192(p192,q192);
					checksum2 += t192.d0;
					res = CMPEQ192(t192, ONE192);

				#elif(defined(P2WORD))

					res = twopmodq128x2(&checksum1, &checksum2, p,k_to_try[l]);

				#else

					if(fbits_in_q < 63)
					{
						itmp64 = p[0]*k_to_try[l];	itmp64 += itmp64 + 1;
						checksum1 += itmp64;
						res = twopmodq63(p[0],itmp64);
						checksum2 += res;
						res = (res == 1);
					}
					else if(fbits_in_q < 64)
					{
						itmp64 = p[0]*k_to_try[l];	itmp64 += itmp64 + 1;
						checksum1 += itmp64;
						res = twopmodq64(p[0],itmp64);
						checksum2 += res;
						res = (res == 1);
					}
					else
					{
						res = twopmodq128_96(&checksum1, &checksum2, p[0],k_to_try[l]);
					}

				#endif	/* endif(P1WORD) */

					if(res == 1)	/* If Lth bit = 1, Lth candidate of the inputs is a factor */
					{
						/* Check if it's a composite factor - if so, skip: */
						if(mi64_pprimeF(q, 3ull, lenQ))
						{
							sprintf(char_buf0,"Factor with k = %llu. Program: E%s\n", k, VERSION);
						#ifdef FACTOR_STANDALONE
						  if(!restart)
							fprintf(stderr,"%s", char_buf0);
						#else
							ASSERT(HERE, fp == 0x0,"0");
							fp = mlucas_fopen(STATFILE,"a");
							fprintf(fp,"%s", char_buf0);
							fclose(fp); fp = 0x0;
						#endif

							fp = mlucas_fopen(   OFILE,"a");
							fprintf(fp,"%s", char_buf0);
							fclose(fp); fp = 0x0;

						#if FAC_DEBUG
							printf("factor was number %u of 0-%u in current batch.\n", l, TRYQM1);
						#endif
							factor_k[nfactor++] = k;

						#ifdef QUIT_WHEN_FACTOR_FOUND
							return 0;
						#endif
						}
					}	/* endif(res == 1) */
				}	/* endfor(l = 0; l <= (uint32)q_index; l++) */
			}	/* endif(q_index >= 0) */
		#endif	/* end #if(TRYQ > 1) */

		#if(!FAC_DEBUG)
			/*
			Every 1024th pass, write the checkpoint file, with format as described previously.
			*/
			if(((sweep + 1) %(1024/lenQ + 1)) == 0 || ((sweep + 1) == interval_hi))
			{
				/* TF restart files are in HRF, not binary: */
				fp = mlucas_fopen(RESTARTFILE, "w");
				if(!fp)
				{
					fprintf(stderr,"INFO: Unable to open factoring savefile %s for writing...quitting.\n",RESTARTFILE);
				#ifndef FACTOR_STANDALONE
					fp = mlucas_fopen(STATFILE,"a");
					fprintf(	fp,"INFO: Unable to open factoring savefile %s for writing...quitting.\n",RESTARTFILE);
					fclose(fp); fp = 0x0;
				#endif
					ASSERT(HERE, 0,"0");
				}
				else
				{
					curr_line = 0;

					/* pstring: */
					++curr_line;
					itmp = fprintf(fp,"%s\n",pstring);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (current exponent) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* bmin: */
					++curr_line;
					itmp = fprintf(fp, "%lf\n", bmin);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (bmin) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* bmax: */
					++curr_line;
					itmp = fprintf(fp, "%lf\n", bmax);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (bmax) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* KMin: */
					++curr_line;
					itmp = fprintf(fp,"KMin = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, kmin)]);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (KMin) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* KNow: */
					++curr_line;
					/* Calculate the current k-value... */
					k = (uint64)incr[pass] + (sweep+1)*(len<<6);
					itmp = fprintf(fp,"KNow = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, k)]);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (KNow) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* KMax: */
					++curr_line;
					itmp = fprintf(fp,"KMax = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, kmax)]);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (KMax) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* PassMin: */
					++curr_line;
					itmp = fprintf(fp,"PassMin = %u\n", passmin);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (PassMin) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* PassNow: */
					++curr_line;
					itmp = fprintf(fp,"PassNow = %u\n", pass   );
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (PassNow) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* PassMax: */
					++curr_line;
					itmp = fprintf(fp,"PassMax = %u\n", passmax);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (PassMax) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* Number of q's tried: */
					++curr_line;
					itmp = fprintf(fp,"#Q tried = %s\n", &char_buf0[convert_uint64_base10_char (char_buf0, count)]);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (#Q tried) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* 64-bit (sum of trial q)%2^64 checksum: */
					++curr_line;
					itmp = fprintf(fp,"Checksum1 = %s\n", &char_buf0[convert_uint64_base16_char (char_buf0, checksum1)]);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (Checksum1) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					/* 64-bit (sum of 2^p % q)%2^64 checksum: */
					++curr_line;
					itmp = fprintf(fp,"Checksum2 = %s\n", &char_buf0[convert_uint64_base16_char (char_buf0, checksum2)]);
					if(itmp <= 0)
					{
						fprintf(stderr,"ERROR: unable to write Line %d (Checksum2) of factoring restart file %s!\n", curr_line, RESTARTFILE);
						ASSERT(HERE, 0,"0");
					}

					fclose(fp); fp = 0x0;
				}
			}	/* Successfully wrote restart file. */
		#endif /* #if !FAC_DEBUG */
	#if DBG_SIEVE
		QUIT:
	#endif
			continue;
		} /* end of sweep loop	*/

		clock2 = clock();
		td += (double)(clock2 - clock1);
		clock1 = clock2;
		tdiff += td;	/* Update total-time accumulator (tdiff)
						prior to dividing per-pass accumulator (td) */
	#ifdef FACTOR_STANDALONE
		if(!restart)
			printf("\n");
	#else

		fp = mlucas_fopen(STATFILE,"a");
		fprintf(fp,"Trial-factoring Pass %2u: time =%s\n", pass,get_time_str(td));
		fclose(fp); fp = 0x0;
	#endif
		td = 0;

	/***********/
	#ifdef ONEPASS
		return 0;
	#endif
	/***********/
	}	/* end of pass loop	*/

/* Any partial-batch factors left to be done in GPU mode? */
#ifdef USE_GPU	/* Set nonzero to print each factor candidate tried */
	kvidx = (count & kvm1);	// kvsz a power of 2, so %kvsz can be effected via &(kvsz-1)
	if(kvidx >= TRYQM1) {
		goto GPU_FACTOR_BATCH;
	}
	GPU_CLEANUP_DONE:
  #ifdef REALLY_GPU
	cudaThreadSynchronize();
  #endif
#endif

/*...all done.	*/
#ifdef FACTOR_STANDALONE
	if(!restart)
	{
	 printf(   "%s(%s) has %u factors in range k = [%llu, %llu], passes %u-%u\n",
	 	NUM_PREFIX[MODULUS_TYPE], pstring, nfactor, kmin, kmax, passmin, passmax);
	 printf(   "Performed %s trial divides\n", &char_buf0[convert_uint64_base10_char(char_buf0, count)]);
	 printf(   "Checksum1 = %s\n", &char_buf0[convert_uint64_base16_char (char_buf0, checksum1)]);
	 printf(   "Checksum2 = %s\n", &char_buf0[convert_uint64_base16_char (char_buf0, checksum2)]);

	/* Since we're done accumulating cycle count, divide to get total time in seconds: */
	 printf(   "Clocks =%s\n",get_time_str(tdiff));
	}
#else
	ASSERT(HERE, fp == 0x0,"0");
	fp = mlucas_fopen(STATFILE,"a");
	fprintf(fp,"Performed %s trial divides\n", &char_buf0[convert_uint64_base10_char(char_buf0, count)]);
	fprintf(fp,"Checksum1 = %s\n", &char_buf0[convert_uint64_base16_char (char_buf0, checksum1)]);
	fprintf(fp,"Checksum2 = %s\n", &char_buf0[convert_uint64_base16_char (char_buf0, checksum2)]);

	/* Since we're done accumulating cycle count, divide to get total time in seconds: */
	fprintf(fp,"Clocks =%s\n",get_time_str(tdiff));
	fclose(fp); fp = 0x0;
#endif

	fp = mlucas_fopen(   OFILE,"a");
  #if(defined(P1WORD))
	 fprintf(fp,"M(%s) has %u factors in range k = [%llu, %llu], passes %u-%u\n", pstring, nfactor, kmin, kmax, passmin, passmax);
  #else
	 fprintf(fp,"M(%s) has %u factors in range k = [%llu, %llu], passes %u-%u\n", pstring, nfactor, kmin, kmax, passmin, passmax);
  #endif
	fclose(fp); fp = 0x0;

#if FAC_DEBUG
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
  #if(defined(P1WORD))
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
	printf(" -passmin {num}  Current factoring pass (0-15).\n");
	printf("\n");
	printf(" -passmax {num}  Maximum pass for the run (0-15).\n");
	/* If we reached here other than via explicit invocation of the help menu, assert: */
	if(!STREQ(stFlag, "-h"))
		ASSERT(HERE, 0,"Mfactor: Unrecognized command-line option!");
	return(0);
#endif
}

/******************/

/* For an exponent p and a factor index k, both mod 60, returns the factoring
pass number (in unit-offset form, i.e. 1-16) on which the factor should occur
if it's one of the valid combinations of p%60 and k%60. If invalid, returns 0. */
uint32 CHECK_PKMOD60(uint32 p, uint32 k)
{
ASSERT(HERE, p < 60 && k < 60, "CHECK_PKMOD60: args must both be (mod 60)!");
if(p== 1){ if(k==0)return 16; else if(k== 3)return 1; else if(k== 8)return 2; else if(k==11)return 3; else if(k==15)return 4; else if(k==20)return 5; else if(k==23)return 6; else if(k==24)return 7; else if(k==35)return 8; else if(k==36)return 9; else if(k==39)return 10; else if(k==44)return 11; else if(k==48)return 12; else if(k==51)return 13; else if(k==56)return 14; else if(k==59)return 15; else return 0; };
if(p== 7){ if(k==0)return 16; else if(k== 5)return 1; else if(k== 8)return 2; else if(k== 9)return 3; else if(k==12)return 4; else if(k==17)return 5; else if(k==20)return 6; else if(k==24)return 7; else if(k==29)return 8; else if(k==32)return 9; else if(k==33)return 10; else if(k==44)return 11; else if(k==45)return 12; else if(k==48)return 13; else if(k==53)return 14; else if(k==57)return 15; else return 0; };
if(p==11){ if(k==0)return 16; else if(k== 1)return 1; else if(k== 4)return 2; else if(k== 9)return 3; else if(k==13)return 4; else if(k==16)return 5; else if(k==21)return 6; else if(k==24)return 7; else if(k==25)return 8; else if(k==28)return 9; else if(k==33)return 10; else if(k==36)return 11; else if(k==40)return 12; else if(k==45)return 13; else if(k==48)return 14; else if(k==49)return 15; else return 0; };
if(p==13){ if(k==0)return 16; else if(k== 3)return 1; else if(k== 8)return 2; else if(k==11)return 3; else if(k==12)return 4; else if(k==15)return 5; else if(k==20)return 6; else if(k==23)return 7; else if(k==27)return 8; else if(k==32)return 9; else if(k==35)return 10; else if(k==36)return 11; else if(k==47)return 12; else if(k==48)return 13; else if(k==51)return 14; else if(k==56)return 15; else return 0; };
if(p==17){ if(k==0)return 16; else if(k== 3)return 1; else if(k== 4)return 2; else if(k== 7)return 3; else if(k==12)return 4; else if(k==15)return 5; else if(k==19)return 6; else if(k==24)return 7; else if(k==27)return 8; else if(k==28)return 9; else if(k==39)return 10; else if(k==40)return 11; else if(k==43)return 12; else if(k==48)return 13; else if(k==52)return 14; else if(k==55)return 15; else return 0; };
if(p==19){ if(k==0)return 16; else if(k== 5)return 1; else if(k== 9)return 2; else if(k==12)return 3; else if(k==17)return 4; else if(k==20)return 5; else if(k==21)return 6; else if(k==24)return 7; else if(k==29)return 8; else if(k==32)return 9; else if(k==36)return 10; else if(k==41)return 11; else if(k==44)return 12; else if(k==45)return 13; else if(k==56)return 14; else if(k==57)return 15; else return 0; };
if(p==23){ if(k==0)return 16; else if(k== 1)return 1; else if(k==12)return 2; else if(k==13)return 3; else if(k==16)return 4; else if(k==21)return 5; else if(k==25)return 6; else if(k==28)return 7; else if(k==33)return 8; else if(k==36)return 9; else if(k==37)return 10; else if(k==40)return 11; else if(k==45)return 12; else if(k==48)return 13; else if(k==52)return 14; else if(k==57)return 15; else return 0; };
if(p==29){ if(k==0)return 16; else if(k== 4)return 1; else if(k== 7)return 2; else if(k==12)return 3; else if(k==15)return 4; else if(k==16)return 5; else if(k==19)return 6; else if(k==24)return 7; else if(k==27)return 8; else if(k==31)return 9; else if(k==36)return 10; else if(k==39)return 11; else if(k==40)return 12; else if(k==51)return 13; else if(k==52)return 14; else if(k==55)return 15; else return 0; };
if(p==31){ if(k==0)return 16; else if(k== 5)return 1; else if(k== 8)return 2; else if(k== 9)return 3; else if(k==20)return 4; else if(k==21)return 5; else if(k==24)return 6; else if(k==29)return 7; else if(k==33)return 8; else if(k==36)return 9; else if(k==41)return 10; else if(k==44)return 11; else if(k==45)return 12; else if(k==48)return 13; else if(k==53)return 14; else if(k==56)return 15; else return 0; };
if(p==37){ if(k==0)return 16; else if(k== 3)return 1; else if(k== 8)return 2; else if(k==12)return 3; else if(k==15)return 4; else if(k==20)return 5; else if(k==23)return 6; else if(k==24)return 7; else if(k==27)return 8; else if(k==32)return 9; else if(k==35)return 10; else if(k==39)return 11; else if(k==44)return 12; else if(k==47)return 13; else if(k==48)return 14; else if(k==59)return 15; else return 0; };
if(p==41){ if(k==0)return 16; else if(k== 3)return 1; else if(k== 4)return 2; else if(k==15)return 3; else if(k==16)return 4; else if(k==19)return 5; else if(k==24)return 6; else if(k==28)return 7; else if(k==31)return 8; else if(k==36)return 9; else if(k==39)return 10; else if(k==40)return 11; else if(k==43)return 12; else if(k==48)return 13; else if(k==51)return 14; else if(k==55)return 15; else return 0; };
if(p==43){ if(k==0)return 16; else if(k== 5)return 1; else if(k== 8)return 2; else if(k==12)return 3; else if(k==17)return 4; else if(k==20)return 5; else if(k==21)return 6; else if(k==32)return 7; else if(k==33)return 8; else if(k==36)return 9; else if(k==41)return 10; else if(k==45)return 11; else if(k==48)return 12; else if(k==53)return 13; else if(k==56)return 14; else if(k==57)return 15; else return 0; };
if(p==47){ if(k==0)return 16; else if(k== 4)return 1; else if(k== 9)return 2; else if(k==12)return 3; else if(k==13)return 4; else if(k==24)return 5; else if(k==25)return 6; else if(k==28)return 7; else if(k==33)return 8; else if(k==37)return 9; else if(k==40)return 10; else if(k==45)return 11; else if(k==48)return 12; else if(k==49)return 13; else if(k==52)return 14; else if(k==57)return 15; else return 0; };
if(p==49){ if(k==0)return 16; else if(k==11)return 1; else if(k==12)return 2; else if(k==15)return 3; else if(k==20)return 4; else if(k==24)return 5; else if(k==27)return 6; else if(k==32)return 7; else if(k==35)return 8; else if(k==36)return 9; else if(k==39)return 10; else if(k==44)return 11; else if(k==47)return 12; else if(k==51)return 13; else if(k==56)return 14; else if(k==59)return 15; else return 0; };
if(p==53){ if(k==0)return 16; else if(k== 3)return 1; else if(k== 7)return 2; else if(k==12)return 3; else if(k==15)return 4; else if(k==16)return 5; else if(k==27)return 6; else if(k==28)return 7; else if(k==31)return 8; else if(k==36)return 9; else if(k==40)return 10; else if(k==43)return 11; else if(k==48)return 12; else if(k==51)return 13; else if(k==52)return 14; else if(k==55)return 15; else return 0; };
if(p==59){ if(k==0)return 16; else if(k== 1)return 1; else if(k== 4)return 2; else if(k== 9)return 3; else if(k==12)return 4; else if(k==16)return 5; else if(k==21)return 6; else if(k==24)return 7; else if(k==25)return 8; else if(k==36)return 9; else if(k==37)return 10; else if(k==40)return 11; else if(k==45)return 12; else if(k==49)return 13; else if(k==52)return 14; else if(k==57)return 15; else return 0; };
return 0;
}

/* This is actually an auxiliary source file, but give it a .h extension to allow wildcarded project builds of form 'gcc -c *.c' */
#include "factor_test.h"

