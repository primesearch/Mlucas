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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef factor_h_included
#define factor_h_included

#include "util.h"

#ifdef __cplusplus
extern "C" {
#endif

/********
PLEASE REFER TO FACTOR.C FOR A DESCRIPTION OF THE APPLICABLE #DEFINES

	(except for FACTOR_PASS_MAX and MAX_BITS_P|Q, which are defined in Mdata.h,
	 and MAX_IMUL_WORDS and MUL_LOHI64_SUBROUTINE, defined in imul_macro.h).
********/

#undef P1WORD
#if(!defined(P2WORD) && !defined(P3WORD) && !defined(P4WORD) && !defined(NWORD))
	#define P1WORD
#endif

#ifdef P3WORD
  #ifndef	PIPELINE_MUL192	// User can build with -DPIPELINE_MUL192=0 to override the default
	#define PIPELINE_MUL192	1
  #endif
#endif

#if !defined(TRYQ) && (defined(NWORD) || defined(P4WORD) || defined(P3WORD))
	#undef TRYQ
	#define TRYQ	1
#endif

#ifdef USE_FLOAT
	/* This will need to be made platform-dependent at some point: */
  #ifndef TRYQ
   #ifdef USE_GPU
	#define TRYQ	1
   #else
	#define TRYQ	4
   #endif
  #elif(TRYQ != 1 && TRYQ != 2 && TRYQ != 4 && TRYQ != 8 && TRYQ != 16 && TRYQ != 32 && TRYQ != 64)
	#error USE_FLOAT option requires TRYQ = 1, 2, 4, 8, 16, 32 or 64
  #endif

	/* FP-based modmul currently only supported for one-word p's: */
	#ifdef P2WORD
		#error P2WORD may not be used together with USE_FLOAT!
	#endif
	#ifdef P3WORD
		#error P3WORD may not be used together with USE_FLOAT!
	#endif
	#ifdef P4WORD
		#error P4WORD may not be used together with USE_FLOAT!
	#endif

	#ifdef USE_65BIT
		#error USE_65BIT may not be used together with USE_FLOAT!
	#endif
	#ifdef USE_95BIT
		#error USE_65BIT may not be used together with USE_FLOAT!
	#endif
#elif defined(MUL_LOHI64_SUBROUTINE)
  #ifndef TRYQ
	#define TRYQ	4
  #elif(TRYQ != 1 && TRYQ != 4 && TRYQ != 4)
	#error MUL_LOHI64_SUBROUTINE option requires TRYQ = 1, 4 or 8
  #endif
#endif

#ifdef USE_FMADD
	/* This will need to be made platform-dependent at some point: */
  #ifndef TRYQ
	#define TRYQ	2
  #elif(TRYQ != 1 && TRYQ != 2 && TRYQ != 4)
	#error USE_FMADD option requires TRYQ = 1, 2 or 4
  #endif

	/* FMADD-based modmul currently only supported for one-word p's: */
	#ifdef P2WORD
		#error P2WORD may not be used together with USE_FMADD!
	#endif
	#ifdef P3WORD
		#error P3WORD may not be used together with USE_FMADD!
	#endif
	#ifdef P4WORD
		#error P4WORD may not be used together with USE_FMADD!
	#endif

	#ifdef USE_FLOAT
		#error USE_FLOAT may not be used together with USE_FMADD!
	#endif

	#ifdef USE_65BIT
		#error USE_65BIT may not be used together with USE_FMADD!
	#endif
	#ifdef USE_95BIT
		#error USE_65BIT may not be used together with USE_FMADD!
	#endif
#endif

// If using x86_64 inline-asm, use ASM-ized 96-bit modpow rather than 128-bit variants:
#if !defined(USE_128x96) && (defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define USE_128x96	1
#endif

/* Special-handling #define for 96-bit factor candidates: */
#ifdef USE_128x96
	/* Currently only values of 0,1,2 supported: */
	#if(USE_128x96 > 2)
		#error Unsupported value of USE_128x96!
	#endif

	/* If 2-word-or-greater factoring not enabled, make sure factors < 2^96: */
	#if((MAX_BITS_Q > 96) && defined(P1WORD))
		#error Factoring > 96 bits requires PxWORD with x > 2 to be defined!
	#endif
#endif

/* Key debug #define off value of EWM_DEBUG (set in masterdefs.h) and FACTOR_STANDALONE, set at compile time: */
#ifndef FAC_DEBUG
	#if EWM_DEBUG && defined(FACTOR_STANDALONE)
		#warning Setting FAC_DEBUG in factor.h
		#define FAC_DEBUG
	#endif
#endif

#ifndef TRYQ	/* # of candidates at a time to try. Must be of form 2^k, or zero to skip trial division step
				(e.g. to test sieve-related timings.) A value of 8 seems optimal on the Alpha 21264, which is
				reasonable - TRYQ < 8 may not allow all the integer MULs to complete by the time their results
				are needed, whereas TRYQ = 16 needs more registers than are available and at the same time
				doesn't pipeline significantly better than TRYQ = 8. */
	#if DBG_SIEVE
		#define TRYQ	0
	#elif defined(USE_GPU)
		#define TRYQ	1
	#elif defined(INTEGER_MUL_32) || !defined(USE_SSE2)
		#define TRYQ	4
	#else
		#define TRYQ	8
	#endif

#endif

/* Don't have 16-input versions of the twopmodq routines for q > 2^63, so only allow TRYQ up to 8. */
#ifndef TRYQ
	#error TRYQ undefined!
#endif
#if(TRYQ != 0 && TRYQ != 1 && TRYQ != 2 && TRYQ != 4 && TRYQ != 8 && TRYQ != 16 && TRYQ != 32 && TRYQ != 64)
	#error	Illegal value of TRYQ
#endif
#if((TRYQ == 2 || TRYQ > 8) && !defined(USE_FLOAT) && !defined(USE_FMADD))
	#error	TRYQ = 2 and TRYQ > 8 only allowed if USE_FLOAT is defined
#endif

/* Make sure the TRYQ = 4 fused macros are only used on the PPC32: */
#ifdef MOD_INI_Q4
	#ifndef CPU_SUBTYPE_PPC32
		#error TRYQ = 4 fused macros should only used on the PPC!
	#endif
#endif

/* Factoring-only globals: */
extern int restart;

#ifdef FACTOR_STANDALONE
	/* Declare a blank STATFILE string to ease program logic: */
	extern char STATFILE[];
#endif

/*******************************************************************************
   Function prototypes. The corresponding function definitions will either
   be in a {function name}.c file or (for cases where a .c file contains
   multiple function definitions) in the given .c file:
*******************************************************************************/

/* factor.c: */

#ifndef FACTOR_STANDALONE
	/* exponents > 64 bits require standalone-mode build: */
	int factor(char *pstring, double log2_min_factor, double log2_max_factor);
#endif

int		test_fac(void);
uint64	factor_qmmp_sieve64(uint32 p, uint64 k, uint64 imin, uint64 imax);

uint64	test_modsqr64    (uint64  x, uint64  q);
uint96	test_modsqr96    (uint96  x, uint96  q);
uint128	test_modsqr128_96(uint128 x, uint128 q);
uint128	test_modsqr128   (uint128 x, uint128 q);

// The single-trial-divisor versions of many of the modpow routines are in terms of exponent and *divisor*
// to allow them to serve both as TF modules and to be used more generally, e.g. for PRP testing.
// The multiple-trial-divisor routines are intended for TF use only, thus input their divisors
// q = 2.k.p+1 using the factor index k:
uint32	test_twopmodq64(uint32 imax);
uint64	twopmodq63    (uint64 p, uint64 q);	// q, not k!
uint64	twopmodq63_q4 (uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq63_q8 (uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);
uint64	twopmodq64    (uint64 p, uint64 q);	// q, not k!
uint64	twopmodq64_q4 (uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq64_q8 (uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);
uint64	twopmodq65    (uint64 p, uint64 k);
uint64	twopmodq65_q4 (uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq65_q8 (uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);

// Conventional positive-power version of twopmodq, returns true mod:
uint64	twopmmodq64    (uint64 p, uint64 q);	// q, not k!
// This variant of twopmodq64_q4 returns the 4 true mods, overwriting the inputs
void	twopmmodq64_q4(uint64 p, uint64 *i0, uint64 *i1, uint64 *i2, uint64 *i3, uint64 *qi0, uint64 *qi1, uint64 *qi2, uint64 *qi3);

uint64 twopmodq63_x8(uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint64 q4, uint64 q5, uint64 q6, uint64 q7);

uint64	twopmodq78_3WORD_DOUBLE   (uint64 p, uint64 k);
uint64	twopmodq78_3WORD_DOUBLE_q2(uint64 p, uint64 k0, uint64 k1, int init_sse2, int thr_id);
uint64	twopmodq78_3WORD_DOUBLE_q4(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, int init_sse2, int thr_id);
uint64	twopmodq78_3WORD_DOUBLE_q4_REF(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);

#if defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
uint64	twopmodq78_3WORD_DOUBLE_q8 (uint64 p, uint64 k[]
		, int init_sse2, int thr_id
	);
uint64	twopmodq78_3WORD_DOUBLE_q16(uint64 p, uint64 k[]
		, int init_sse2, int thr_id
	);
uint64	twopmodq78_3WORD_DOUBLE_q32(uint64 p, uint64 k[]
		, int init_sse2, int thr_id
	);
uint64	twopmodq78_3WORD_DOUBLE_q64(uint64 p, uint64 k[]
		, int init_sse2, int thr_id
	);

// Routines optimized for vector-FMA-based 100-plus-bit product computation (Intel AVX2 and beyond):
uint64	twopmodq100_2WORD_DOUBLE_q32(uint64 p, uint64 k[]
		, int init_sse2, int thr_id
	);
uint64	twopmodq100_2WORD_DOUBLE_q64(uint64 p, uint64 k[]
		, int init_sse2, int thr_id
	);
#endif

uint96	twopmodq96    (uint64 p, uint64 k);
// The'init_sse2' arg is ill-named here as it refers to 'init any local-mem needed for inline-asm routines' whether SSE/AVX used or
// not - here we use just x86_64 non-vector instructions - but use same name as in vector-radix* code to show analogous mem-handling:
uint64	twopmodq96_q4 (uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, int init_sse2, int thr_id);
uint64	twopmodq96_q8 (uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7, int init_sse2, int thr_id);
uint64	twopmodq96_q8_const_qhi(uint64 p, uint64 khi, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);

uint64	twopmodq128_96   (uint64 p, uint64 k);
uint64	twopmodq128_96_q4(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq128_96_q8(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);
uint64	twopmodq128_96_q8_const_qhi(uint64 p, uint64 khi, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);

uint128 twopmmodq128	(uint128 p, uint128 q);
uint128	twopmodq128		(uint128 p, uint128 q);	// q, not k!
uint64	twopmodq128x2	(uint64 *p, uint64 k);
uint64	twopmodq128x2B	(uint64 *p, uint128 q);	// q, not k!
uint64	twopmodq128_q4	(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq128_q8	(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);

uint64	twopmodq160   (uint64 *p, uint64 k);
uint64	twopmodq160_q4(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq160_q8(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);

uint192 twopmmodq192	(uint192 p, uint192 q);
uint192	twopmodq192		(uint192 p, uint192 q);	// q, not k!
uint64	twopmodq192_q4	(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq192_q4_qmmp(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq192_q8	(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);

uint256	twopmodq200_8WORD_DOUBLE(uint64 *p, uint64 k);
uint256	twopmodq200_8WORD_qmmp  (uint64 *p, uint64 k);

uint256 twopmmodq256	(uint256 p, uint256 q);
uint256	twopmodq256		(uint256 p, uint256 q);	// q, not k!
uint64	twopmodq256_q4	(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3);
uint64	twopmodq256_q8	(uint64 *p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7);

uint32	CHECK_PKMOD60  (uint64 *p, uint32 lenP, uint64 k, uint32*incr);
uint32	CHECK_PKMOD4620(uint64 *p, uint32 lenP, uint64 k, uint32*incr);

/******************************************/
/*            GPU-TF stuff:               */
/******************************************/
#ifdef FLOAT_PINV
	typedef float 	pinv_t;
#else
	typedef uint32	pinv_t;
#endif

	// Top-level routines for CPU-parallel and GPU-side sieving and testing of resulting factor candidates:
#ifdef MULTITHREAD

	void*				// Thread-arg pointer *must* be cast to void and specialized inside the function
	PerPass_tfSieve(void*thread_arg);

#else

  #ifdef USE_GPU	// GPU version currently only supports 32-bit p and 64-bit k
	__host__
  #endif
	uint64 PerPass_tfSieve(
		const char *pstring,
		const uint32 pass,
		const uint64 interval_lo, const uint64 interval_hi,
		const double fbits_in_2p,
		const uint32 nclear,
		const uint32 sieve_len,	// #64-bit words in the full-length sieving bitmap setup by the (CPU-side) caller.
								// This routine uses only a specific (len/60)-sized "sievelet" excerpt.
								// sieve_len = 64*(product of first NCLEAR odd primes)
		const uint32 p_last_small,	//largest odd prime appearing in the product; that is, the (nclear)th odd prime.
		const uint32 nprime,		// #sieving primes (counting from 3)
		const uint32 MAX_SIEVING_PRIME,
	#if defined(USE_AVX512) && !defined(USE_GPU)
		const uint32 *psmall,
	#endif
		const uint8 *pdiff,
		      uint32*startval,	// This gets updated within
			  uint64*k_to_try,	// Unused by GPU code; include to yield a (mostly) uniform API
			  uint64*factor_k,	// List of found factors for each p gets updated (we hope) within
			  uint32*nfactor,	// Here the '*' is to denote a writeable scalar
		const uint32 findex,
		const uint64*p, const uint64*two_p, const uint32 lenP, const uint32 lenQ, const uint32 incr,
	#ifndef USE_GPU
		uint32*kdeep, uint32*ndeep, const uint64 countmask, const uint32 CMASKBITS,
		uint64*q, uint64*q2, uint64*u64_arr,
	#endif
		const uint64 kstart,
		const uint64*bit_map, uint64*bit_map2,	// GPU version uses only 1st of these
		double *tdiff,
		const int MODULUS_TYPE,
		const char*VERSION,
		const char*OFILE
	);

#endif

#ifdef USE_GPU	// GPU-only sieving/TFing functions:

	__global__ void popcountBitmap      (uint32*bit_map, const uint32 nthr, const uint32 nwords);
	__global__ void clrPrimes_Lt64      (uint64*bit_map, const uint32 nthr, const pinv_t*pinv, const uint32*startval);
	__global__ void clrPrimes_Gt64_Lt128(uint32*bit_map, const uint32 nthr, const uint32 nwords, const uint32 np, const uint8*pdiff, const pinv_t*pinv, const uint32*startval);
	__global__ void clrPrimes_GtA_LtB   (uint32*bit_map, const uint32 nthr, const uint32 nwords, const uint32 np, const uint8*pdiff, const pinv_t*pinv, const uint32*startval, const uint32 curr_p);
	__global__ void clrPrimes_LuchaLibre(uint32*bit_map, const uint32 nthr, const uint32 bit_len, const uint32*prim, uint32*startval);

	// All-in-one version of the 3-function sequence below:
	__global__ void parseBitmap_InitKoff3in1(      uint32*bit_map, uint32*kvec,                        const uint32 nthr, const uint32 nwords);

	__global__ void parseBitmap_CountSetBits(const uint32*bit_map,                    uint32*popc_vec, const uint32 nthr, const uint32 nwords);
	__global__ void parseBitmap_PopcountSums(                      uint32*kvec,       uint32*popc_vec, const uint32 nthr);
	__global__ void parseBitmap_InitKoffsets(const uint32*bit_map, uint32*kvec, const uint32*popc_vec, const uint32 nthr, const uint32 nwords);

	__global__ void updateStartval(
		const uint32 bit_len,
		const uint32 nthr,
		const uint32 nclear,
		const uint32*p100,
		const uint8 *pdiff,
			  uint32*startval
	);

	__host__ uint32 GPU_tfBatch(	// TRYQ treated as 1 in this version
		const uint64*p, const uint64*two_p, const uint32 lenP, const uint32 lenQ,
		const uint64 k0,			// k-value corr. to the 0-bit in the bit_map
		const uint32 nclear,
		const uint32 nprime,		// #sieving primes (counting from 3)
		const uint32*gpu_prim,
		const uint32*gpu_p100,
		const uint8 *gpu_pdiff,
		const pinv_t*gpu_pinv,
		const double fbits_in_q,	// Bit-ness of largest candidate-q in each batch
		uint32*gpu_startval,
		uint32*gpu_bit_map,
		uint32*gpu_kvec, const uint32 kvsz,
		uint32*gpu_popc_vec,
		const uint32 bit_map_len,	// #32-bit words in GPU-bitmap (all bits assumed significant, i.e. no padding)
		uint64*kfac,
		int last_batch
	);

  #ifdef __CUDACC__

// Self-tests:
	__global__ void VecModpow64(const uint64*pvec, const uint64*pshft, const uint32*zshft, const uint32*stidx, const uint64*kvec, uint8*rvec, int N);
	__global__ void VecModpow78(const uint64*pvec, const uint64*pshft, const uint32*zshft, const uint32*stidx, const uint64*kvec, uint8*rvec, int N);
	__global__ void VecModpow96(const uint64*pvec, const uint64*pshft, const uint32*zshft, const uint32*stidx, const uint64*kvec, uint8*rvec, int N);

// Production TFing:

	// 64-bit, global functions (simple distributors-of-work-to-parallel-device-function-callees):
 __global__ void GPU_TF64_pop64(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64 k0, uint32*kvec, const uint32 N);
	__global__ void GPU_TF64   (const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	__global__ void GPU_TF64_q4(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	__global__ void GPU_TF64_q8(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	// 64-bit, device functions:
	__device__ uint32 twopmodq64_gpu(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);
	__device__ uint32 twopmodq64_q4_GPU(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k0, const uint64 k1, const uint64 k2, const uint64 k3,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);
	__device__ uint32 twopmodq64_q8_GPU(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k0, const uint64 k1, const uint64 k2, const uint64 k3, const uint64 k4, const uint64 k5, const uint64 k6, const uint64 k7,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);

	// 78-bit, global functions (simple distributors-of-work-to-parallel-device-function-callees):
 __global__ void GPU_TF78_pop64(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64 k0, uint32*kvec, const uint32 N);
	__global__ void GPU_TF78   (const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	__global__ void GPU_TF78_q4(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	__global__ void GPU_TF78_q8(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	// 78-bit, device functions:
	__device__ uint32 twopmodq78_gpu(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);
	__device__ uint32 twopmodq78_q4_GPU(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k0, const uint64 k1, const uint64 k2, const uint64 k3,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);
	__device__ uint32 twopmodq78_q8_GPU(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k0, const uint64 k1, const uint64 k2, const uint64 k3, const uint64 k4, const uint64 k5, const uint64 k6, const uint64 k7,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);

	// 96-bit, global functions (simple distributors-of-work-to-parallel-device-function-callees):
 __global__ void GPU_TF96_pop64(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64 k0, uint32*kvec, const uint32 N);
	__global__ void GPU_TF96   (const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	__global__ void GPU_TF96_q4(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	__global__ void GPU_TF96_q8(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N);
	// 96-bit, device functions:
	__device__ uint32 twopmodq96_gpu(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);
	__device__ uint32 twopmodq96_q4_GPU(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k0, const uint64 k1, const uint64 k2, const uint64 k3,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);
	__device__ uint32 twopmodq96_q8_GPU(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k0, const uint64 k1, const uint64 k2, const uint64 k3, const uint64 k4, const uint64 k5, const uint64 k6, const uint64 k7,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	);

  #endif	// __CUDACC__ ?

#endif	// USE_GPU ?

// Computes 2*p (mod curr_p):
uint32	twop_mod_smallp(const int MODULUS_TYPE, const uint64*two_p, const uint32 findex, const uint32 lenP, const uint32 curr_p);

void	get_startval(
	const int MODULUS_TYPE,
	const uint64 p,		// Only need LSW of this
	const uint32 findex,// Double-Mersenne and Fermat cases
	const uint64*two_p,	// Here, need the full multiword array (of which use just LSW if P!WORD def'd)
	const uint32 lenP,	// Manyword case
	const uint32 bit_len,
	const uint32 interval_lo, const uint32 incr,
	const uint32 nclear, const uint32 nprime, const uint32 p_last_small,
	const uint8 *pdiff,
	      uint32*startval
);

uint64 given_b_get_k(double bits, const uint64 two_p[], uint32 len);

int read_savefile(const char*fname, const char*pstring, double*bmin, double*bmax,
uint64*kmin, uint64*know, uint64*kmax, uint32*passmin, uint32*pass, uint32*passmax, uint64*count);

int init_savefile(const char*fname, const char*pstring, double bmin, double bmax,
uint64 kmin, uint64 know, uint64 kmax, uint32 passmin, uint32 passnow, uint32 passmax, uint64 count);

int write_savefile(const char*fname, const char*pstring, uint32 passnow, uint64 know, uint64 count);

#ifdef __cplusplus
}
#endif

#endif	/* factor_h_included */

