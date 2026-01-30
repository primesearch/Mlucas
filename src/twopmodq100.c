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

#ifdef USE_AVX512

	#include "factor.h"
	#include "twopmodq100.h"
	#include "align.h"
	
	/***********************************************************************************/
	/**** 100-BIT INPUTS, USING FLOATING-DOUBLE FMA-ARITHMETIC FOR MODMUL **************/
	/***********************************************************************************/

	uint64 twopmodq100_2WORD_DOUBLE_q32(uint64 p, uint64 k[], int init_sse2, int thr_id)
	{
//p = 16264097; k[0] = k[1] = k[2] = k[3] = k[4] = k[5] = k[6] = k[7] = 204476354235ull;	// Debug
		const char func[] = "twopmodq100_2WORD_DOUBLE_q32";
		 int32 j = 0;	/* This needs to be signed because of the LR binary exponentiation. */
		uint64 lead7, r = 0;
	#ifdef USE_AVX512_I
		uint64 tmp64;
	#endif
		static uint64 psave = 0, pshift;
		static uint32 start_index, zshift, first_entry = TRUE;
		const uint8* minv8_ptr = minv8;	// Ptr to Table of precomputed byte-inverses def'd in mi64.h
		static int max_threads = 1;	// Default local-array-init is for just a single thread ... caller can re-init for > 1 threads later, if desired.
	#ifdef USE_AVX512_I
	  #error AVX-512 IFMA instruction extensions not yet supported!
		const uint64 mask_lo52 = 0x000FFFFFFFFFFFFFull;
		static uint64 *sc_arr = 0x0, *sc_ptr;
		uint64 *itmp;
		vec_u64 *tmp;
	  #ifdef MULTITHREAD
		static uint64 *__r0;	// Base address for discrete per-thread local stores. *NOTE*: more convenient to use uint64* rather than vec_u64* here
	  #else
		static	// Following set of pointers only static-izable in single-thread mode
	  #endif
		uint64 *fq0[32]/* ,*fq1[32] */, *fqinv0[32]/* ,*fqinv1[32] */, *fx0[32],*fx1[32];
		for(j = 0; j < 32; j++) {
			ASSERT((k[j] >> 52) == 0ull, "Ks must be < 2^52!");
		}
	#else
		const double crnd = 3.0*0x4000000*0x2000000;	// Const used to emulate DNINT(x) and (when multiplied by BASE) 2^50 * DNINT(x*2^bpow2)
		const uint32 bpow2 = 48, pow0 = bpow2-32, pow1 = 32-pow0;	// Powers-of-2 used in double->uint32 conversion: say bpow2 = 52 -> need to concatenate
																	// high 20 (pow0) bits of fq0 and low 12 (pow1) bits of fq1 to get middle 32 bits of q.
		const double base = (double)(1ull<<bpow2), binv = 1.0/base, base0 = (double)(1u<<pow0), binv0 = 1.0/base0, base1 = (double)(1u<<pow1), binv1 = 1.0/base1;
		static double *sc_arr = 0x0, *sc_ptr;
		const double *dptr0,*dptr1,*dptr2,*dptr3,*dptr4,*dptr5,*dptr6,*dptr7;
		double *dptr,dtmp;
		vec_dbl *tmp,*tm1;
	  #ifdef MULTITHREAD
		static double *__r0;	// Base address for discrete per-thread local stores ...
	  #else
		static	// Following set of pointers only static-izable in single-thread mode
	  #endif
		double *fq0[32]/* ,*fq1[32] */, *fqinv0[32]/* ,*fqinv1[32] */, *fx0[32],*fx1[32], kdbl[32];
		// AVX-512 Foundation lacks the needed DQ extensions, so use HLL to convert kvec entries to double:
		for(j = 0; j < 32; j++) {
			ASSERT((k[j] >> 52) == 0ull, "Ks must be < 2^52!");
			kdbl[j] = (double)k[j];
		}
		ASSERT(base  == (double)(1ull<<bpow2) && base <= TWO48FLOAT && base *binv  == 1.0, "Current version only supports bpow2 <= 48; base*binv must == 1.0!");
		ASSERT(base0 == (double)(1ull<< pow0)                       && base0*binv0 == 1.0, "base0,binv0 check fails!");
		ASSERT(base1 == (double)(1ull<< pow1)                       && base1*binv1 == 1.0, "base1,binv1 check fails!");
	#endif
		uint32 mul_width = bpow2<<1;
		if(p != psave) {
		//	first_entry = FALSE;
			psave  = p;
			pshift = p + mul_width;
			j = leadz64(pshift);
			/* Extract leftmost 7 bits of pshift (if >= mul_width, use the leftmost 6) and subtract from mul_width-1: */
			lead7 = ((pshift<<j) >> 57);	// lead7 in [64,127]
			if(lead7 >= mul_width) {
				lead7 >>= 1;	// lead7 in [bpow2,63]
				start_index =  64-j-6;	/* Use only the leftmost 6 bits */
			} else {
				start_index =  64-j-7;
			}
			// lead7 in [bpow2, 2*bpow2):
			zshift = mul_width-1 - lead7;	// zshift in [0,bpow2-1]
			ASSERT(zshift < bpow2, "zshift out of expected range!");
			pshift = ~pshift;
		}

		/* In order to eschew complex thread-block-and-sync logic related to the local-store-init step in multithread mode,
		switch to a special-init-mode-call paradigm, in which the function is inited once (for as many threads as needed)
		prior to being executed:
		*/
		if(first_entry || init_sse2)	// Just check (init_sse2 != 0) here, to allow the *value* of init_sse2 to store #threads
		{
			first_entry = FALSE;
			if(init_sse2) {
			#ifndef MULTITHREAD
				max_threads = 1;
			#else
				max_threads = init_sse2;
			#endif
				fprintf(stderr, "%s: Setting up for as many as %d threads...\n",func,max_threads);
			#ifndef COMPILER_TYPE_GCC
				ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
			#endif
				ASSERT(max_threads >= NTHREADS, "Multithreading requires max_threads >= NTHREADS!");
				ASSERT(thr_id == -1, "Init-mode call must be outside of any multithreading!");
			}
			if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage (e.g. unthreaded call to the function)
				free((void *)sc_arr);	sc_arr=0x0;
			}
			// Alloc the local-memory block the #bytes multiplier has plenty of extra room built in, e.g. for debug-data-writes:
		#ifdef USE_AVX512_I

			sc_arr = ALLOC_UINT64(sc_arr, 0x140*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			sc_ptr = (uint64 *)ALIGN_VEC_U64(sc_arr);	// Force vec_u64-alignment
			ASSERT(((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		  #ifdef MULTITHREAD
			__r0  = sc_ptr;
		  #else
			/* Remember, these are POINTERS-TO-UINT64, so need an increment of 8 to span an AVX-512 register: */
			fq0   [0] = sc_ptr + 0x000;
			//fq1   [0] = sc_ptr + 0x020;
			// +0x40 - Insert another slot here for q.lo + base*q.hi approximant
			fqinv0[0] = sc_ptr + 0x060;
			//fqinv1[0] = sc_ptr + 0x080;
			fx0   [0] = sc_ptr + 0x0A0;
			fx1   [0] = sc_ptr + 0x0C0;
			// +0xE0,100 - Insert another 2 pairs of padding slots here for high-product-words (hi.lo,hi.hi) register spills
			// +0x120 - Insert another slot here for hi.lo + base*h.hi approximant
			for(j = 1, itmp = sc_ptr+1; j < 32; j++, itmp++) {
				fq0   [j] = itmp + 0x000;
				//fq1   [j] = itmp + 0x020;
				fqinv0[j] = itmp + 0x060;
				//fqinv1[j] = itmp + 0x080;
				fx0   [j] = itmp + 0x0A0;
				fx1   [j] = itmp + 0x0C0;
			}
			// +0xE0,100 - Insert another 2 pairs of padding slots here for high-product-words register spills
		  #endif

		/***************************************************/
		#else	// Default AVX-512 floating-point-FMA mode
		/***************************************************/

			sc_arr = ALLOC_DOUBLE(sc_arr, 0x140*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);	// Force vec_u64-alignment
			ASSERT(((uintptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		  #ifdef MULTITHREAD
			__r0  = sc_ptr;
		  #else
			/* Remember, these are POINTERS-TO-DOUBLE, so need an increment of 8 to span an AVX-512 register: */
			fq0   [0] = sc_ptr + 0x000;
			//fq1   [0] = sc_ptr + 0x020;
			// +0x40 - Insert another slot here for q.lo + base*q.hi approximant
			fqinv0[0] = sc_ptr + 0x060;
			//fqinv1[0] = sc_ptr + 0x080;
			fx0   [0] = sc_ptr + 0x0A0;
			fx1   [0] = sc_ptr + 0x0C0;
			// +0xE0,100 - Insert another 2 pairs of padding slots here for high-product-words (hi.lo,hi.hi) register spills
			// +0x120 - Insert another slot here for hi.lo + base*h.hi approximant
			for(j = 1, dptr = sc_ptr+1; j < 32; j++, dptr++) {
				fq0   [j] = dptr + 0x000;
				//fq1   [j] = dptr + 0x020;
				fqinv0[j] = dptr + 0x060;
				//fqinv1[j] = dptr + 0x080;
				fx0   [j] = dptr + 0x0A0;
				fx1   [j] = dptr + 0x0C0;
			}
		  #endif

		#endif
			if(init_sse2) return 0;
		}	/* end of inits */

		/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD

		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		sc_ptr = __r0 + thr_id*0x140;

	  #ifdef USE_AVX512_I
		// Cf. above init-block for memory layout
		/* Remember, these are POINTERS-TO-UINT64, so need an increment of 8 to span an AVX-512 register: */
		fq0   [0] = sc_ptr + 0x000;
		//fq1   [0] = sc_ptr + 0x020;
		fqinv0[0] = sc_ptr + 0x060;
		//fqinv1[0] = sc_ptr + 0x080;
		fx0   [0] = sc_ptr + 0x0A0;
		fx1   [0] = sc_ptr + 0x0C0;
		for(j = 1, itmp = sc_ptr+1; j < 32; j++, itmp++) {
			fq0   [j] = itmp + 0x000;
			//fq1   [j] = itmp + 0x020;
			fqinv0[j] = itmp + 0x060;
			//fqinv1[j] = itmp + 0x080;
			fx0   [j] = itmp + 0x0A0;
			fx1   [j] = itmp + 0x0C0;
		}

	  #else
		// Cf. above init-block for memory layout
		/* Remember, these are POINTERS-TO-DOUBLE, so need an increment of 8 to span an AVX-512 register: */
		fq0   [0] = sc_ptr + 0x000;
		//fq1   [0] = sc_ptr + 0x020;
		fqinv0[0] = sc_ptr + 0x060;
		//fqinv1[0] = sc_ptr + 0x080;
		fx0   [0] = sc_ptr + 0x0A0;
		fx1   [0] = sc_ptr + 0x0C0;
		for(j = 1, dptr = sc_ptr+1; j < 32; j++, dptr++) {
			fq0   [j] = dptr + 0x000;
			//fq1   [j] = dptr + 0x020;
			fqinv0[j] = dptr + 0x060;
			//fqinv1[j] = dptr + 0x080;
			fx0   [j] = dptr + 0x0A0;
			fx1   [j] = dptr + 0x0C0;
		}
		// +0xC0,E0 - Insert another 2 pairs of padding slots here for high-product-words register spills

	  #endif

	#endif

	#ifdef MUL_LOHI64_SUBROUTINE
		#error MUL_LOHI64_SUBROUTINE defined!
	#endif
		ASSERT((p >> 63) == 0, "twopmodq100_q32: p must be < 2^63!");

	#ifdef USE_AVX512_I

		tmp64 = p+p;	// uint64
		tmp = fq0[0];	// vec_u64*
		__asm__ volatile (\
			"movq	%[__fq0] ,%%rax				\n\t"\
			"movq	%[__p2]  ,%%rbx				\n\t"\
			"movq	%[__kvec],%%rcx				\n\t"\
			"movq	%[__mask26],%%rsi			\n\t"\
			"vpsrlq	$25,(%%rsi),%%zmm1			\n\t"/* zmm1 = 0x3ffffff >> 25 = 1 */\
			"vpxorq	%%zmm0,%%zmm0,%%zmm0		\n\t"/* zmm0 = 0 */\
			"vpbroadcastq	%%rbx,%%zmm2		\n\t"/* broadcast 2*p to all 8 qwords of zmm */\
			"vmovdqu64	(%%rcx),%%zmm3			\n\t"/* load the factor-candidate k's, 8 per zmm ... note kvec not nec. 16-byte aligned */\
		... 52-bit mull,mulh ...
			"vmovaps	%%zmm0,0x000(%%rax)		\n\t"/* store lo52 */\
			"vmovaps	%%zmm1,0x040(%%rax)		\n\t"/* store hi52 */\
			:					/* outputs: none */\
			: [__fq0]  "m" (tmp)	/* All inputs from memory addresses here */\
			 ,[__p2]   "m" (tmp64)\
			 ,[__kvec] "m" (k)\
			 ,[__mask26] "m" (mask_lo26)\
			: "cc","memory","rax","rbx","rcx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
		);

	#else
	/*** CF notes in twopmodq78_3WORD_DOUBLE_q32 routine on qinv-computation ... need to add high-32-bits accumulator to go beyond 96 bits here ***/
		dtmp = p+p;	dptr = &dtmp;	// dtmp = double; dptr = double*
		dptr0 = &base; dptr1 = &binv; dptr2 = &base0; dptr3 = &binv0; dptr4 = &base1; dptr5 = &binv1; dptr6 = &TWO32FLOAT; dptr7 = &TWO32FLINV;
		tmp = (vec_dbl*)fq0[0];	tm1 = (vec_dbl*)kdbl;	// Use vec_dbl* ptrs to hold address args ... asm-code cares not a whit about types
		__asm__ volatile (\
			"movq	%[__fq0] ,%%rax				\n\t"\
			"movq	%[__p2]  ,%%rbx				\n\t"\
			"movq	%[__kvec],%%rcx				\n\t"/* Each column uses 7 distinct vector registers, so increment all non-shared vec-reg indices by +7 for each additional rightward column: */\
			"vbroadcastsd	(%%rbx),%%zmm1		\n\t	vmovaps	%%zmm1,%%zmm8				\n\t	vmovaps	%%zmm1,%%zmm15				\n\t	vmovaps	%%zmm1,%%zmm22				\n\t"/* broadcast 2*p to all 8 qwords of zmm ... double-prec version is from mem-address */\
			"vmovupd	(%%rcx),%%zmm3			\n\t	vmovupd	0x40(%%rcx),%%zmm10			\n\t	vmovupd	0x80(%%rcx),%%zmm17			\n\t	vmovupd	0xc0(%%rcx),%%zmm24			\n\t"/* load the factor-candidate k's, 8 per zmm ... note kvec not nec. 16-byte aligned */\
			"movq	%[__base],%%rbx				\n\t	vbroadcastsd	(%%rbx),%%zmm30		\n\t"/* BASE ~= 50 bits */\
			"movq	%[__binv],%%rcx				\n\t	vbroadcastsd	(%%rcx),%%zmm31		\n\t"/* BINV */\
			"vmulpd		%%zmm1,%%zmm3,%%zmm2	\n\t	vmulpd		%%zmm8,%%zmm10,%%zmm9	\n\t	vmulpd		%%zmm15,%%zmm17,%%zmm16	\n\t	vmulpd		%%zmm22,%%zmm24,%%zmm23	\n\t"/* hi = x*y = 2p*k */\
			"vmovaps	%%zmm2,%%zmm0			\n\t	vmovaps	%%zmm9,%%zmm7				\n\t	vmovaps	%%zmm16,%%zmm14				\n\t	vmovaps	%%zmm23,%%zmm21				\n\t"/* cpy hi into lo-destination reg */\
		"vfmsub231pd	%%zmm1,%%zmm3,%%zmm0	\n\t vfmsub231pd	%%zmm8,%%zmm10,%%zmm7	\n\t vfmsub231pd	%%zmm15,%%zmm17,%%zmm14	\n\t vfmsub231pd	%%zmm22,%%zmm24,%%zmm21	\n\t"/* lo = fma(x,y,-hi; xmm1,3 FREE */\
			"vmovaps	%%zmm2,0x200(%%rax)		\n\t	vmovaps	%%zmm9,0x240(%%rax)			\n\t	vmovaps	%%zmm16,0x280(%%rax)		\n\t	vmovaps	%%zmm23,0x2c0(%%rax)		\n\t"/* write hi-approximant to full 2-word q into scratch slot set aside for it */\
			"movl	$1,%%esi					\n\t"\
			"vpbroadcastd	%%esi,%%zmm3		\n\t	vpbroadcastd	%%esi,%%zmm10		\n\t	vpbroadcastd	%%esi,%%zmm17		\n\t	vpbroadcastd	%%esi,%%zmm24		\n\t"/* Broadcast to all 8 int32 slots of zmm3, since KNL lacks AVX-512 VL extensions */\
			"vcvtdq2pd	%%ymm3,%%zmm3			\n\t	vcvtdq2pd	%%ymm10,%%zmm10			\n\t	vcvtdq2pd	%%ymm17,%%zmm17			\n\t	vcvtdq2pd	%%ymm24,%%zmm24			\n\t"/* zmm3 = (double)1, conversion from 4-way int32 in lower half of zmm3 */\
			"vaddpd		%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd		%%zmm10,%%zmm7,%%zmm7	\n\t	vaddpd		%%zmm17,%%zmm14,%%zmm14	\n\t	vaddpd		%%zmm24,%%zmm21,%%zmm21	\n\t"/* lo += 1 */\
		/* Normalize to properly split the resultant significant bits across lo,hi, use base-2^52 signed-int-stored-as-double: */\
			"vmovaps	%%zmm2,%%zmm1			\n\t	vmovaps	%%zmm9,%%zmm8				\n\t	vmovaps	%%zmm16,%%zmm15				\n\t	vmovaps	%%zmm23,%%zmm22				\n\t"/* cpy hi into hh-destination reg */\
			"vmulpd		%%zmm31,%%zmm1,%%zmm1	\n\t	vmulpd		%%zmm31,%%zmm8,%%zmm8	\n\t	vmulpd		%%zmm31,%%zmm15,%%zmm15	\n\t	vmulpd		%%zmm31,%%zmm22,%%zmm22	\n\t"/*            hi*BINV  */\
			"vrndscalepd	$1,%%zmm1,%%zmm1	\n\t	vrndscalepd	$1,%%zmm8,%%zmm8		\n\t	vrndscalepd	$1,%%zmm15,%%zmm15		\n\t	vrndscalepd	$1,%%zmm22,%%zmm22		\n\t"/* hh = FLOOR(hi*BINV) */\
			"vmovaps	%%zmm2,%%zmm3			\n\t	vmovaps	%%zmm9,%%zmm10				\n\t	vmovaps	%%zmm16,%%zmm17				\n\t	vmovaps	%%zmm23,%%zmm24				\n\t"/* cpy hi into cy-destination reg */\
		"vfnmadd231pd	%%zmm30,%%zmm1,%%zmm3	\n\t vfnmadd231pd	%%zmm30,%%zmm8,%%zmm10	\n\t vfnmadd231pd	%%zmm30,%%zmm15,%%zmm17	\n\t vfnmadd231pd	%%zmm30,%%zmm22,%%zmm24	\n\t"/* cy = fma(hh,-BASE,hi)	Backward carry from hi into lo */\
			"vaddpd		%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd		%%zmm10,%%zmm7,%%zmm7	\n\t	vaddpd		%%zmm17,%%zmm14,%%zmm14	\n\t	vaddpd		%%zmm24,%%zmm21,%%zmm21	\n\t"/* lo += cy */\
			"vmovaps	%%zmm1,0x100(%%rax)		\n\t	vmovaps	%%zmm8,0x140(%%rax)			\n\t	vmovaps	%%zmm15,0x180(%%rax)		\n\t	vmovaps	%%zmm22,0x1c0(%%rax)		\n\t"/* store hi half (hh) into fq1 */\
			"vmovaps	%%zmm0,0x000(%%rax)		\n\t	vmovaps	%%zmm7,0x040(%%rax)			\n\t	vmovaps	%%zmm14,0x080(%%rax)		\n\t	vmovaps	%%zmm21,0x0c0(%%rax)		\n\t"/* store lo half into fq0 */\
		/*** lo in zmm0, hi in zmm1. ***/\
		/* Assemble 32-bit pieces of q from the floating-double ones ... AVX-512F has just double->uint32 conversions, so emulate the needed double->quad conversion using floating-double arithmetic: */\
			"movq	%[__two32f],%%rbx			\n\t	vbroadcastsd	(%%rbx),%%zmm30		\n\t"/* 2^+32 */\
			"movq	%[__two32i],%%rcx			\n\t	vbroadcastsd	(%%rcx),%%zmm31		\n\t"/* 2^-32 */\
			"vmulpd		%%zmm31,%%zmm0,%%zmm2	\n\t	vmulpd	 	%%zmm31,%%zmm7,%%zmm9	\n\t	vmulpd		%%zmm31,%%zmm14,%%zmm16	\n\t	vmulpd		%%zmm31,%%zmm21,%%zmm23	\n\t"/*       lo*TWO32FLINV  */\
			"vrndscalepd	$1,%%zmm2,%%zmm2	\n\t	vrndscalepd	$1,%%zmm9,%%zmm9		\n\t	vrndscalepd	$1,%%zmm16,%%zmm16		\n\t	vrndscalepd	$1,%%zmm23,%%zmm23		\n\t"/* FLOOR(lo*TWO32FLINV) = upper pow0 bits of lo */\
		"vfnmadd231pd	%%zmm30,%%zmm2,%%zmm0	\n\t vfnmadd231pd	%%zmm30,%%zmm9,%%zmm7	\n\t vfnmadd231pd	%%zmm30,%%zmm16,%%zmm14	\n\tvfnmadd231pd	%%zmm30,%%zmm23,%%zmm21	\n\t"/* lo32, unsigned (>= 0) */\
			"movq	%[__base1],%%rbx			\n\t	vbroadcastsd	(%%rbx),%%zmm30		\n\t"/* 2^+b, b = #bits needed from low end of fq1 */\
			"movq	%[__binv1],%%rcx			\n\t	vbroadcastsd	(%%rcx),%%zmm31		\n\t"/* 2^-b */\
			"vmulpd	 	%%zmm31,%%zmm1,%%zmm3	\n\t	vmulpd	 	%%zmm31,%%zmm8,%%zmm10	\n\t	vmulpd		%%zmm31,%%zmm15,%%zmm17	\n\t	vmulpd		%%zmm31,%%zmm22,%%zmm24	\n\t"/*       lo*BINV1  */\
			"vrndscalepd	$1,%%zmm3,%%zmm3	\n\t	vrndscalepd	$1,%%zmm10,%%zmm10		\n\t	vrndscalepd	$1,%%zmm17,%%zmm17		\n\t	vrndscalepd	$1,%%zmm24,%%zmm24		\n\t"/* FLOOR(lo*BINV1) = upper 32 bits of hi */\
		"vfnmadd231pd	%%zmm30,%%zmm3,%%zmm1	\n\t vfnmadd231pd	%%zmm30,%%zmm10,%%zmm8	\n\t vfnmadd231pd	%%zmm30,%%zmm17,%%zmm15	\n\t vfnmadd231pd	%%zmm30,%%zmm24,%%zmm22	\n\t"/* lower pow1 bits of hi; (pow0+pow1) = 32 */\
			"vcvttpd2udq		  %%zmm0,%%ymm0	\n\t	vcvttpd2udq		  %%zmm7,%%ymm7		\n\t	vcvttpd2udq		  %%zmm14,%%ymm14	\n\t	vcvttpd2udq		  %%zmm21,%%ymm21	\n\t"/* lo 32 bits of q */\
			"vcvttpd2udq		  %%zmm1,%%ymm1	\n\t	vcvttpd2udq		  %%zmm8,%%ymm8		\n\t	vcvttpd2udq		  %%zmm15,%%ymm15	\n\t	vcvttpd2udq		  %%zmm22,%%ymm22	\n\t"/* lower pow1 bits of q.hi; (pow0+pow1) = 32 */\
			"vcvttpd2udq		  %%zmm2,%%ymm2	\n\t	vcvttpd2udq		  %%zmm9,%%ymm9		\n\t	vcvttpd2udq		  %%zmm16,%%ymm16	\n\t	vcvttpd2udq		  %%zmm23,%%ymm23	\n\t"/* upper pow0 bits of q.lo */\
			"vcvttpd2udq		  %%zmm3,%%ymm3	\n\t	vcvttpd2udq		  %%zmm10,%%ymm10	\n\t	vcvttpd2udq		  %%zmm17,%%ymm17	\n\t	vcvttpd2udq		  %%zmm24,%%ymm24	\n\t"/* hi 32 bits of q */\
	/*** Must use full-length zmm-regs here (even though upper 256 bits unused) due to AVX-512F full-width constraint (applying to ymm needs AVX-512VL extensions): ***/\
			"movl	%[__pow0],%%ebx				\n\t	vmovd	%%ebx,%%xmm31				\n\t"/* VPSLLD on 16-way int32 data ... inanely, there is no take-shift-count-from-GPR version of this, must stick count into low end of an xmm-reg */\
			"vpslld	%%xmm31,%%zmm1,%%zmm1		\n\t	vpslld	%%xmm31,%%zmm8,%%zmm8		\n\t	vpslld	%%xmm31,%%zmm15,%%zmm15		\n\t	vpslld	%%xmm31,%%zmm22,%%zmm22		\n\t"/* left-justify the low pow1 bits of q.hi... */\
			"vpaddd	%%zmm2,%%zmm1,%%zmm1		\n\t	vpaddd	%%zmm9,%%zmm8,%%zmm8		\n\t	vpaddd	%%zmm16,%%zmm15,%%zmm15		\n\t	vpaddd	%%zmm23,%%zmm22,%%zmm22		\n\t"/* ...and add to upper pow0 bits of q.lo to get q.md32... */\
			"vmovaps	%%zmm3,%%zmm2			\n\t	vmovaps	%%zmm10,%%zmm9				\n\t	vmovaps	%%zmm17,%%zmm16				\n\t	vmovaps	%%zmm24,%%zmm23				\n\t"/* ...and cpy hi 32 bits of q into m2,9,16,23. */\
		/* Zero-extend each uint32 in ymm0-2 to a uint64, leaving q.lo32,md32,hi32 in qword-form in zmm0,1,2: */\
			"vpmovzxdq	%%ymm0,%%zmm0			\n\t	vpmovzxdq	%%ymm7,%%zmm7			\n\t	vpmovzxdq	%%ymm14,%%zmm14			\n\t	vpmovzxdq	%%ymm21,%%zmm21			\n\t"
			"vpmovzxdq	%%ymm1,%%zmm1			\n\t	vpmovzxdq	%%ymm8,%%zmm8			\n\t	vpmovzxdq	%%ymm15,%%zmm15			\n\t	vpmovzxdq	%%ymm22,%%zmm22			\n\t"\
			"vpmovzxdq	%%ymm2,%%zmm2			\n\t	vpmovzxdq	%%ymm9,%%zmm9			\n\t	vpmovzxdq	%%ymm16,%%zmm16			\n\t	vpmovzxdq	%%ymm23,%%zmm23			\n\t"\
		/* Gather-load the needed byte-sized qinv initializers: */\
			"movq	%[__minv8],%%rbx			\n\t"\
			"movl	$-1,%%esi					\n\t"\
			"vpbroadcastd	%%esi,%%zmm3		\n\t	vpbroadcastd	%%esi,%%zmm10		\n\t	vpbroadcastd	%%esi,%%zmm17		\n\t	vpbroadcastd	%%esi,%%zmm24		\n\t"/* Mask-reg = 0x11...11; no longer use VPCMPEQD for this since in AVX-512
														that outputs one result bit per zmm-subfield-compare to a mask-reg */\
			"vpsrlq	$56,%%zmm3,%%zmm3			\n\t	vpsrlq	$56,%%zmm10,%%zmm10			\n\t	vpsrlq	$56,%%zmm17,%%zmm17			\n\t	vpsrlq	$56,%%zmm24,%%zmm24			\n\t"/* Each qword = 0x0000000000FF */\
			"vpandq	%%zmm3,%%zmm0,%%zmm4		\n\t	vpandq	%%zmm10,%%zmm7,%%zmm11		\n\t	vpandq	%%zmm17,%%zmm14,%%zmm18		\n\t	vpandq	%%zmm24,%%zmm21,%%zmm25		\n\t"/* AND with mask-reg leaves just low byte of each uint64 subfield */\
			"vpsrlq	$01,%%zmm4,%%zmm4			\n\t	vpsrlq	$01,%%zmm11,%%zmm11			\n\t	vpsrlq	$01,%%zmm18,%%zmm18			\n\t	vpsrlq	$01,%%zmm25,%%zmm25			\n\t"/* Lookup indices in minv byte-table */\
		/* Note VPGATHERDD may not use default-opmask k0: */
			"movl	$-1,%%esi					\n\t"\
			"kmovw	%%esi,%%k1					\n\t	kmovw	%%esi,%%k2					\n\t	kmovw	%%esi,%%k3					\n\t	kmovw	%%esi,%%k4					\n\t"/* Init opmask k1 (Only need the low byte) */\
	"vpgatherqq (%%rbx,%%zmm4),%%zmm5%{%%k1%}	\n\tvpgatherqq (%%rbx,%%zmm11),%%zmm12%{%%k2%}\n\tvpgatherqq (%%rbx,%%zmm18),%%zmm19%{%%k3%}\n\tvpgatherqq (%%rbx,%%zmm25),%%zmm26%{%%k4%}\n\t"/* minv8 is a byte-table ... instruction sets mask-reg = 0, so each col uses separate same-valued k-reg */\
			"vpandq	%%zmm3,%%zmm5,%%zmm3		\n\t	vpandq	%%zmm10,%%zmm12,%%zmm10		\n\t	vpandq	%%zmm17,%%zmm19,%%zmm17		\n\t	vpandq	%%zmm24,%%zmm26,%%zmm24		\n\t"/* Iterative qinv-computation doesn't care about the garbage upper 7 bytes of each qword
				resulting from the gather-load from the packed byte-array, but zero all-but-low-byte of each uint64 subfield to ease debug */\
		/* Newton iteration as described in comment above this asm: q0-2 in zmm0-2, qinv0-2 in zmm3-5; currently only have low byte of qinv: */\
			"movq	$2,%%rcx					\n\t"\
			"vpbroadcastq	%%rcx,%%zmm30		\n\t"/* vector-int64 2 */\
		/* 1. q-data are 32-bit (of which we need the low 16, but just use q32 and VPMULLD), qinv are 8-bit: */\
			"vpbroadcastd	%%esi,%%zmm31		\n\t"/* Mask-reg = 0x11...11 */\
			"vpsrlq	$48,%%zmm31,%%zmm31			\n\t"/* Each qword = 0x00000000FFFF */\
			"vpmulld	%%zmm0,%%zmm3,%%zmm4	\n\t	vpmulld	%%zmm7,%%zmm10,%%zmm11		\n\t	vpmulld	%%zmm14,%%zmm17,%%zmm18		\n\t	vpmulld	%%zmm21,%%zmm24,%%zmm25		\n\t"/* tmp16  = MULL16(q16,qinv8); Use MULL32 for simplicity */\
			"vpsubd		%%zmm4,%%zmm30,%%zmm4	\n\t	vpsubd		%%zmm11,%%zmm30,%%zmm11	\n\t	vpsubd		%%zmm18,%%zmm30,%%zmm18	\n\t	vpsubd		%%zmm25,%%zmm30,%%zmm25	\n\t"/* 2-tmp16 */\
			"vpmulld	%%zmm3,%%zmm4,%%zmm3	\n\t	vpmulld	%%zmm10,%%zmm11,%%zmm10		\n\t	vpmulld	%%zmm17,%%zmm18,%%zmm17		\n\t	vpmulld	%%zmm24,%%zmm25,%%zmm24		\n\t"/* MULL16(qinv8,(2-tmp16)); Use MULL32 for simplicity, then mask off hi48 bits to get... */\
			"vpandq		%%zmm3,%%zmm31,%%zmm3	\n\t	vpandq		%%zmm10,%%zmm31,%%zmm10	\n\t	vpandq		%%zmm17,%%zmm31,%%zmm17	\n\t	vpandq		%%zmm24,%%zmm31,%%zmm24	\n\t"/* qinv16 [Don't in fact need the masking-off of hi16 bits, but useful for debugging] */\
		/* 2. q-data are 32-bit, qinv are 16-bit: */\
			"vpmulld	%%zmm0,%%zmm3,%%zmm4	\n\t	vpmulld	%%zmm7,%%zmm10,%%zmm11		\n\t	vpmulld	%%zmm14,%%zmm17,%%zmm18		\n\t	vpmulld	%%zmm21,%%zmm24,%%zmm25		\n\t"/* tmp32  = MULL32(q32,qinv16); */\
			"vpsubd		%%zmm4,%%zmm30,%%zmm4	\n\t	vpsubd		%%zmm11,%%zmm30,%%zmm11	\n\t	vpsubd		%%zmm18,%%zmm30,%%zmm18	\n\t	vpsubd		%%zmm25,%%zmm30,%%zmm25	\n\t"/* 2-tmp32 */\
			"vpmulld	%%zmm3,%%zmm4,%%zmm3	\n\t	vpmulld	%%zmm10,%%zmm11,%%zmm10		\n\t	vpmulld	%%zmm17,%%zmm18,%%zmm17		\n\t	vpmulld	%%zmm24,%%zmm25,%%zmm24		\n\t"/* qinv.lo32 = MULL32(qinv16,(2-tmp32)); */\
		/* 3. q-data are 64-bit, qinv-data are 32-bit: */\
			"vpmuludq	%%zmm0,%%zmm3,%%zmm4	\n\t	vpmuludq	%%zmm7,%%zmm10,%%zmm11	\n\t	vpmuludq	%%zmm14,%%zmm17,%%zmm18	\n\t	vpmuludq	%%zmm21,%%zmm24,%%zmm25	\n\t"/* VPMULUDQ(q64.lo32,qinv32) */\
			"vpmulld	%%zmm1,%%zmm3,%%zmm5	\n\t	vpmulld	%%zmm8,%%zmm10,%%zmm12		\n\t	vpmulld	%%zmm15,%%zmm17,%%zmm19		\n\t	vpmulld	%%zmm22,%%zmm24,%%zmm26		\n\t"/* VPMULLD (q64.hi32,qinv32)<<32 */\
			"vpsllq	$32,%%zmm5,%%zmm5			\n\t	vpsllq	$32,%%zmm12,%%zmm12			\n\t	vpsllq	$32,%%zmm19,%%zmm19			\n\t	vpsllq	$32,%%zmm26,%%zmm26			\n\t"\
			"vpaddq		%%zmm4,%%zmm5,%%zmm4	\n\t	vpaddq		%%zmm11,%%zmm12,%%zmm11	\n\t	vpaddq		%%zmm18,%%zmm19,%%zmm18	\n\t	vpaddq		%%zmm25,%%zmm26,%%zmm25	\n\t"/* tmp64  = MULL64(q64,qinv32); */\
			/* */\
			"vpsubq		%%zmm4,%%zmm30,%%zmm4	\n\t	vpsubq		%%zmm11,%%zmm30,%%zmm11	\n\t	vpsubq		%%zmm18,%%zmm30,%%zmm18	\n\t	vpsubq		%%zmm25,%%zmm30,%%zmm25	\n\t"/* 2-tmp64 */\
			"vpmuludq	%%zmm4,%%zmm3,%%zmm5	\n\t	vpmuludq	%%zmm11,%%zmm10,%%zmm12	\n\t	vpmuludq	%%zmm18,%%zmm17,%%zmm19	\n\t	vpmuludq	%%zmm25,%%zmm24,%%zmm26	\n\t"/* VPMULUDQ((2-tmp64).lo32,qinv32) */\
			"vpsrlq	$32,%%zmm4,%%zmm4			\n\t	vpsrlq	$32,%%zmm11,%%zmm11			\n\t	vpsrlq	$32,%%zmm18,%%zmm18			\n\t	vpsrlq	$32,%%zmm25,%%zmm25			\n\t"/* (2-tmp64).hi32 */\
			"vpmulld	%%zmm4,%%zmm3,%%zmm4	\n\t	vpmulld	%%zmm11,%%zmm10,%%zmm11		\n\t	vpmulld	%%zmm18,%%zmm17,%%zmm18		\n\t	vpmulld	%%zmm25,%%zmm24,%%zmm25		\n\t"/* VPMULLD ((2-tmp64).hi32,qinv32)<<32 */\
			"vpsllq	$32,%%zmm4,%%zmm4			\n\t	vpsllq	$32,%%zmm11,%%zmm11			\n\t	vpsllq	$32,%%zmm18,%%zmm18			\n\t	vpsllq	$32,%%zmm25,%%zmm25			\n\t"\
			"vpaddq		%%zmm4,%%zmm5,%%zmm4	\n\t	vpaddq		%%zmm11,%%zmm12,%%zmm11	\n\t	vpaddq		%%zmm18,%%zmm19,%%zmm18	\n\t	vpaddq		%%zmm25,%%zmm26,%%zmm25	\n\t"/* qinv64 = MULL64(qinv32,(2-tmp64)); */\
			"vpsrlq	$32,%%zmm4,%%zmm4			\n\t	vpsrlq	$32,%%zmm11,%%zmm11			\n\t	vpsrlq	$32,%%zmm18,%%zmm18			\n\t	vpsrlq	$32,%%zmm25,%%zmm25			\n\t"/* qinv.md32 ['md' refers to middle 32-bits of the eventual 96-bit inverse; lo32 in zmm3] */\
		/* 4. 64 => 96-bits: q-data are 96-bit (in our integer-qinv-compute-then-truncate-to-78-bit and cvt-to-double model), qinv-data are 64-bit: */\
			"vpmuludq	%%zmm0,%%zmm3,%%zmm5	\n\t	vpmuludq	%%zmm7,%%zmm10,%%zmm12	\n\t	vpmuludq	%%zmm14,%%zmm17,%%zmm19	\n\t	vpmuludq	%%zmm21,%%zmm24,%%zmm26	\n\t"/* A. VPMULUDQ(q.lo32, qinv.lo32) */\
			"vpsrlq	$32,%%zmm5,%%zmm5			\n\t	vpsrlq	$32,%%zmm12,%%zmm12			\n\t	vpsrlq	$32,%%zmm19,%%zmm19			\n\t	vpsrlq	$32,%%zmm26,%%zmm26			\n\t"/* right-shift 64-bit result 32 bits to keep just the high halves; */\
			"vpmuludq	%%zmm1,%%zmm3,%%zmm6	\n\t	vpmuludq	%%zmm8,%%zmm10,%%zmm13	\n\t	vpmuludq	%%zmm15,%%zmm17,%%zmm20	\n\t	vpmuludq	%%zmm22,%%zmm24,%%zmm27	\n\t"/* B. VPMULUDQ(q.md32, qinv.lo32), VPADDQ to result of [A], */\
			"vpaddq		%%zmm6,%%zmm5,%%zmm5	\n\t	vpaddq		%%zmm13,%%zmm12,%%zmm12	\n\t	vpaddq		%%zmm20,%%zmm19,%%zmm19	\n\t	vpaddq		%%zmm27,%%zmm26,%%zmm26	\n\t"/* can ignore any carry because it would be into bit 96 of full-length product; */\
			"vpmuludq	%%zmm0,%%zmm4,%%zmm6	\n\t	vpmuludq	%%zmm7,%%zmm11,%%zmm13	\n\t	vpmuludq	%%zmm14,%%zmm18,%%zmm20	\n\t	vpmuludq	%%zmm21,%%zmm25,%%zmm27	\n\t"/* C. VPMULUDQ(q.lo32, qinv.md32), */\
			"vpaddq		%%zmm6,%%zmm5,%%zmm5	\n\t	vpaddq		%%zmm13,%%zmm12,%%zmm12	\n\t	vpaddq		%%zmm20,%%zmm19,%%zmm19	\n\t	vpaddq		%%zmm27,%%zmm26,%%zmm26	\n\t"/* VPADDQ to result of [B], again can ignore any carry, */\
			"vpsrlq	$32,%%zmm5,%%zmm5			\n\t	vpsrlq	$32,%%zmm12,%%zmm12			\n\t	vpsrlq	$32,%%zmm19,%%zmm19			\n\t	vpsrlq	$32,%%zmm26,%%zmm26			\n\t"/* right-shift 64-bit result 32 bits; */\
			"vpmulld	%%zmm1,%%zmm4,%%zmm6	\n\t	vpmulld	%%zmm8,%%zmm11,%%zmm13		\n\t	vpmulld	%%zmm15,%%zmm18,%%zmm20		\n\t	vpmulld	%%zmm22,%%zmm25,%%zmm27		\n\t"/* D. VPMULLD (q.md32, qinv.md32), */\
			"vpaddq		%%zmm6,%%zmm5,%%zmm5	\n\t	vpaddq		%%zmm13,%%zmm12,%%zmm12	\n\t	vpaddq		%%zmm20,%%zmm19,%%zmm19	\n\t	vpaddq		%%zmm27,%%zmm26,%%zmm26	\n\t"/* VPADDD to result of [C] automatically masks off high 32 bits of sum. */\
			"vpbroadcastd	%%esi,%%zmm31		\n\t"/* Mask-reg = 0x11...11 */\
			"vpsrlq	$32,%%zmm31,%%zmm31			\n\t"/* 0 all but bottom 32 bits of each qword */\
			"vpandq		%%zmm5,%%zmm31,%%zmm5	\n\t	vpandq		%%zmm12,%%zmm31,%%zmm12	\n\t	vpandq		%%zmm19,%%zmm31,%%zmm19	\n\t	vpandq		%%zmm26,%%zmm31,%%zmm26	\n\t"/* mp32  = MULH64(q.lo64 * qinv) % 2^32 */\
			"vpmulld	%%zmm2,%%zmm3,%%zmm6	\n\t	vpmulld	%%zmm9,%%zmm10,%%zmm13		\n\t	vpmulld	%%zmm16,%%zmm17,%%zmm20		\n\t	vpmulld	%%zmm23,%%zmm24,%%zmm27		\n\t"/* MULL32(q.hi32, qinv.lo32) */\
			"vpaddd		%%zmm5,%%zmm6,%%zmm5	\n\t	vpaddd		%%zmm12,%%zmm13,%%zmm12	\n\t	vpaddd		%%zmm19,%%zmm20,%%zmm19	\n\t	vpaddd		%%zmm26,%%zmm27,%%zmm26	\n\t"/* tmp32 += MULL32(q.hi32, qinv.lo32); */\
			"vpxorq		%%zmm6,%%zmm6,%%zmm6	\n\t	vpxorq		%%zmm13,%%zmm13,%%zmm13	\n\t	vpxorq		%%zmm20,%%zmm20,%%zmm20	\n\t	vpxorq		%%zmm27,%%zmm27,%%zmm27	\n\t"/* 0 */\
			"vpsubd		%%zmm3,%%zmm6,%%zmm6	\n\t	vpsubd		%%zmm10,%%zmm13,%%zmm13	\n\t	vpsubd		%%zmm17,%%zmm20,%%zmm20	\n\t	vpsubd		%%zmm24,%%zmm27,%%zmm27	\n\t"/* -qinv.lo32 */\
			"vpmulld	%%zmm6,%%zmm5,%%zmm5	\n\t	vpmulld	%%zmm13,%%zmm12,%%zmm12		\n\t	vpmulld	%%zmm20,%%zmm19,%%zmm19		\n\t	vpmulld	%%zmm27,%%zmm26,%%zmm26		\n\t"/* qinv.hi32 = MULL32(-qinv.lo32, tmp32); */\
		/* qinv.lo|md|hi32 in zmm3-5,10-12,17-19,24-26: */\
		/* Down-convert each uint64 field to a uint32, leaving qinv.lo32,md32,hi32 in dword-form in ymm3-5... */\
			"vpmovqd	%%zmm3,%%ymm3			\n\t	vpmovqd	%%zmm10,%%ymm10				\n\t	vpmovqd	%%zmm17,%%ymm17				\n\t	vpmovqd	%%zmm24,%%ymm24				\n\t"
			"vpmovqd	%%zmm4,%%ymm4			\n\t	vpmovqd	%%zmm11,%%ymm11				\n\t	vpmovqd	%%zmm18,%%ymm18				\n\t	vpmovqd	%%zmm25,%%ymm25				\n\t"\
			"vpmovqd	%%zmm5,%%ymm5			\n\t	vpmovqd	%%zmm12,%%ymm12				\n\t	vpmovqd	%%zmm19,%%ymm19				\n\t	vpmovqd	%%zmm26,%%ymm26				\n\t"\
		/* ...then up-convert packed uint32s to to packed doubles in zmm3-5: */\
			"vcvtudq2pd	%%ymm3,%%zmm3			\n\t	vcvtudq2pd	%%ymm10,%%zmm10			\n\t	vcvtudq2pd	%%ymm17,%%zmm17			\n\t	vcvtudq2pd	%%ymm24,%%zmm24			\n\t"/* (double)fqinv.lo32 */\
			"vcvtudq2pd	%%ymm4,%%zmm4			\n\t	vcvtudq2pd	%%ymm11,%%zmm11			\n\t	vcvtudq2pd	%%ymm18,%%zmm18			\n\t	vcvtudq2pd	%%ymm25,%%zmm25			\n\t"/* (double)fqinv.md32 */\
			"vcvtudq2pd	%%ymm5,%%zmm5			\n\t	vcvtudq2pd	%%ymm12,%%zmm12			\n\t	vcvtudq2pd	%%ymm19,%%zmm19			\n\t	vcvtudq2pd	%%ymm26,%%zmm26			\n\t"/* (double)fqinv.hi32 */\
		/* Reassemble the three 32-bit chunks into 2 bpow2-bit pieces and write to local data store: */\
			"movq	%[__base0],%%rbx			\n\t	vbroadcastsd	(%%rbx),%%zmm30		\n\t"/* 2^#bits needing to be moved into high end of 2-word fqinv0 */\
			"movq	%[__binv0],%%rcx			\n\t	vbroadcastsd	(%%rcx),%%zmm31		\n\t"\
			"vmulpd	 	%%zmm31,%%zmm4,%%zmm6	\n\t	vmulpd	 	%%zmm31,%%zmm11,%%zmm13	\n\t	vmulpd		%%zmm31,%%zmm18,%%zmm20	\n\t	vmulpd		%%zmm31,%%zmm25,%%zmm27	\n\t"/*       fqinv.md32*BINV0  */\
			"vrndscalepd	$1,%%zmm6,%%zmm6	\n\t	vrndscalepd	$1,%%zmm13,%%zmm13		\n\t	vrndscalepd	$1,%%zmm20,%%zmm20		\n\t	vrndscalepd	$1,%%zmm27,%%zmm27		\n\t"/* FLOOR(fqinv.md32*BINV0) = (upper pow1 bits of fqinv.md32), needing appending at low end of 2-word fqinv.hi */\
		"vfnmadd231pd	%%zmm30,%%zmm6,%%zmm4	\n\t vfnmadd231pd	%%zmm30,%%zmm13,%%zmm11	\n\t vfnmadd231pd	%%zmm30,%%zmm20,%%zmm18	\n\t vfnmadd231pd	%%zmm30,%%zmm27,%%zmm25	\n\t"/* lower pow0 bits of qinv.md32 */\
			"movq	%[__two32f],%%rbx			\n\t	vbroadcastsd	(%%rbx),%%zmm30		\n\t"/* 2^#bits needing to be moved into high end of 2-word fqinv0 */\
			"movq	%[__base1] ,%%rcx			\n\t	vbroadcastsd	(%%rcx),%%zmm31		\n\t"/* 2^pow1 */\
		"vfmadd231pd	%%zmm30,%%zmm4,%%zmm3	\n\t vfmadd231pd	%%zmm30,%%zmm11,%%zmm10	\n\t vfmadd231pd	%%zmm30,%%zmm18,%%zmm17	\n\t vfmadd231pd	%%zmm30,%%zmm25,%%zmm24	\n\t"/* fqinv.lo32 + (lower pow0 bits of fqinv.md32)<<32 */\
		"vfmadd231pd	%%zmm31,%%zmm5,%%zmm6	\n\t vfmadd231pd	%%zmm31,%%zmm12,%%zmm13	\n\t vfmadd231pd	%%zmm31,%%zmm19,%%zmm20	\n\t vfmadd231pd	%%zmm31,%%zmm26,%%zmm27	\n\t"/* fqinv.hi32<<pow1 + (upper pow1 bits of fqinv.md32) */\
		/* Store 2 words of floating-double qinv: */\
			"movq	%[__fqinv0] ,%%rax			\n\t"\
			"vmovaps	%%zmm3,0x000(%%rax)		\n\t	vmovaps	%%zmm10,0x040(%%rax)		\n\t	vmovaps	%%zmm17,0x080(%%rax)		\n\t	vmovaps	%%zmm24,0x0c0(%%rax)		\n\t"/* store fqinv.lo */\
			"vmovaps	%%zmm6,0x100(%%rax)		\n\t	vmovaps	%%zmm13,0x140(%%rax)		\n\t	vmovaps	%%zmm20,0x180(%%rax)		\n\t	vmovaps	%%zmm27,0x1c0(%%rax)		\n\t"/* store fqinv.hi */\
			:					/* outputs: none */\
			: [__fq0]  "m" (tmp)	/* All inputs from memory addresses here */\
			 ,[__fqinv0] "m" (fqinv0)	\
			 ,[__p2]   "m" (dptr)\
			 ,[__kvec] "m" (tm1)	/* vec-of-doubles copy of input kvec */\
			 ,[__minv8] "m" (minv8_ptr)	/* Ptr to Table of precomputed byte-inverses def'd in mi64.h */\
			 ,[__pow0] "m" (pow0)	/* uint32 (since VPBROADCASTD takes an actual int32, as opposed to a pointer to one) */\
			 ,[__base]   "m" (dptr0), [__binv]   "m" (dptr1)	/* double-ptrs */\
			 ,[__base0]  "m" (dptr2), [__binv0]  "m" (dptr3)	/* double-ptrs */\
			 ,[__base1]  "m" (dptr4), [__binv1]  "m" (dptr5)	/* double-ptrs */\
			 ,[__two32f] "m" (dptr6), [__two32i] "m" (dptr7)	/* double-ptrs */\
			: "cc","memory","rax","rbx","rcx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27", "xmm30","xmm31"	/* Clobbered registers */\
		);

	#endif
		// Init the modpow residue:
	ASSERT(zshift < 48, "zshift out of expected range!");
		dtmp = 1ull<<zshift;	VEC_DBL_INIT_8((vec_dbl*)fx0[0], dtmp);	VEC_DBL_INIT_8((vec_dbl*)fx0[8], dtmp);	VEC_DBL_INIT_8((vec_dbl*)fx0[16], dtmp);	VEC_DBL_INIT_8((vec_dbl*)fx0[24], dtmp);
		dtmp = 0ull;			VEC_DBL_INIT_8((vec_dbl*)fx1[0], dtmp);	VEC_DBL_INIT_8((vec_dbl*)fx1[8], dtmp);	VEC_DBL_INIT_8((vec_dbl*)fx1[16], dtmp);	VEC_DBL_INIT_8((vec_dbl*)fx1[24], dtmp);
/* Debug:
dtmp = 299706231258404;	VEC_DBL_INIT_8((vec_dbl*)fx0[0], dtmp);
dtmp = 5364;	VEC_DBL_INIT_8((vec_dbl*)fx1[0], dtmp);
*/
		// Modular powering loop, returns 2^-p (mod q):
		dptr0 = &base; dptr1 = &binv; dptr2 = &crnd;
		for(j = start_index-1; j >= 0; j--) {
			SSE2_twopmodq100_modmul_q32(fq0[0],fqinv0[0],fx0[0], dptr0,dptr1,dptr2, pshift,j);
		}
#if 0
====================================
bc-based mont_sqr:
a = 2^48; b = a^2
p = 16264097; k = 204476354235; q = 2*k*p+1; q0 = q%a; q1 = q/a
i = 184103730195591 + a*42018343201564
x = 2^33

j = 15:												match?
u = 78505709017419349267191180560; h = 28783919		yes
v = 75936893064530124722627314544 [mull(lo,qinv)]	yes
l = 6374942692452671604 [umulh(lo,q)]				...605 = 378894886445173 + a*22647; 1 too high! Low word should = 378894886445172

Assemble UMULH result in size-a pieces as code does:
l0 = 18075101070192; l1 = 269782038715967

t=q0*l0	 4957896702225535922795195664 = 156340089793808 + a* 17613987432076, discard lo half; code gives 17613987432076.555
u=q1*l0	           427096563187566768 =  99023517501616 + a*           1517, t/a + u     = 427114177174998844, code gives ?
v=q0*l1	73999668100078109760798910089 = 162270801356425 + a*262899633085844, t/a + u + v = 73999668100505223937973908933, code gives ?
																				/= a = 262899633087361.99088, code gives 262899633087362, incorrect round-upward!
w=q1*l1	          6374679792819584243 = 115995253357811 + a*          22647				^^^^^^^^ lg() = 47.9, i.e. only 5 bits for frac-part
****************************************** This causes e.g. 90 of 1030 of my 63-bit TF self-tests to fail ***************************************
Add lo halves of middle 2 subproducts to hi half of q0*l0:
17613987432076 + (99023517501616 + 162270801356425) = 278908306290117 = 0.99088*base, should give cy = 0 into upper 96 bits.
Add hi halves of middle 2 subproducts to lo half of q1*l1:
(1517+262899633085844) + 115995253357811 = 378894886445172, one less than the incorrect result given by the code.

Here is bitfield governing 2x in modpow: start_index = 18, so need low 18 bits, 111101001111111110:
j	bit	mod_sqr sequence																		xout			code gives same?
--	--	------------------------------------------------------------------------------------	---------------		----
24	1	
23	1	
22	1	
21	1	
20	1	
19	1	
18	1	
17	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	532072169072628280	
16	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	1510131481307217188	
15	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	552607653089827812	...810, 2 too low!
14	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	4464895331244842697	
13	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				3784880979688691560	
12	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	2183008738903427601	
11	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				1968888628584671906	
10	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				3670564915036980298	
9	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	2460527599989537806	
8	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	4480761456286949510	
7	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	1402897681754837451	
6	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	4165896150199927982	
5	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	2525760091890533381	
4	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	4723403449100646885	
3	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	137959594232784167	
2	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	1014048173418414995	
1	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	5164180314503158571	
0	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				3325623259484400796	
...and a final mod-doubling gives x = 1, as desired.
473!
====================================
Check if adding a fixed multiple (e.g. 3*base) of base to ensure y[1] > 0 produces montmul output matching that of the careful-carry version:
Ex:    a = 2^48, b = a^2 = 2^96,
	   q = 79065663617798299760500594651 = 34061411080155 + a*280897664658391
	qinv = 74272008985203238536365558867 = 72444849245267 + a*263867182273725
Init x = 96-bit random input:
	   x = 72839422333764741456945352976
====================================
y = x^2 = lo+b*hi; lo = 126972778979584+a*134886089954027 = 37967059028401327043499991296,
					hi = 66965852514392496320999605105 = 39607200400241+a*237910500240419 = -241867776310415+a*237910500240418
					*** my code gives hi = 321082177110897+a*237910500240418, i.e. lo-word has extra b, hi-word is -1. ***
	lo = mull(lo,qinv) = 30075412459459750629810877184 = 189201160057600+a*106849329240289
					*** my code gives -92273816653056+a*106849329240290, i.e. lo-word -= b, hi-word +=1. ***
compare opcount so far:
100-bit modmul: 10 mul, 11 fma, 23 add, 2 round
 78-bit modmul: 10 mul, 15 fma,  5 add, 7 round ... each round ~= 2 add, nets out to 4 more fma, 6 fewer add

	lo = umulh(lo,  q) = 30013727053912913255727430808 = 18786916039832+a*106630178656221, hi-lo = 36952125460479583065272174297 .
a=2^48; q0=34061411080155; q1=280897664658391; q =q0+a*q1
l0= -92273816653056; l1=106849329240290; l =l0+a*l1
Assemble in size-a pieces as code will do:
q0*l0	 -3142976400954592595841703680 = -154502197731072 + a*-11166095251818, discard lo half
q1*l0	-25919499606959979281206192896 = -260621830433536 + a*-92084560801310, this (full 96-bit) + q0*lo.lo = -25919499606959990447301444714, code gives -2.5919499606959991e+28
q0*l1	  3639438926892343434445444950 =  125199038928726 + a* 12929884458729, += -25919499606959979281206192896 = -22280060680067635846760747946, code gives -2.228006068006765e+28
																				/= a = -79154676342581, code gives -79154676342581.531, floor() = -79154676342582.
q1*l1	 30013727053912992410403773390 =   97941592382414 + a*106630178656221

Add lo halves of middle 2 subproducts to hi half of q0*l0:
-11166095251818 + (-260621830433536+125199038928726) = -146588886756628, gives cy = -1 into upper 96 bits.
Add hi halves of middle 2 subproducts to lo half of q1*l1:
(-92084560801310+12929884458729) + 97941592382414 = 18786916039833, += cy (-=1) ==> 18786916039832, this is < a so no cyin to high 48 bits.
Thus get result
lo =  18786916039832+a*106630178656221, as expected.
hi = 321082177110897+a*237910500240418, thus hi > lo, compute x = (hi-lo) = 302295261071065+a*131280321584197
bit = 1, so we double this (mod q):
2*x = 604590522142130+a*262560643168394, which is < q, so no modding needed.

q0 = 34061411080155; q1= 280897664658391
t0 =604590522142130;t1=262560643168394

opcount - umulh adds 5 mul, 3 fma, 5 add, 1 round:
100-bit modmul: 15 mul, 14 fma, 28 add, 3 round
 78-bit modmul: 15 mul, 20 fma,  6 add, 8 round ... each round ~= 2 add, nets out to 6 more fma, 12 fewer add
====================================
#endif
		/* Do a final mod-doubling and return. Code here is specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2 or -(q-1)/2:
		*/
		dptr = fq0[0];
		__asm__ volatile (\
		/* Assume x.lo,hi still in zmm1,zmm2; binv,base still in zmm30,zmm31: */\
			"movq	%[__fq0] ,%%rax	\n\t"\
			/* Mod-doubling code taken from SSE2_twopmodq78_modmul_q16, expects x.lo26,hi52 in zmm1,2, resp., so copy zmm0,4,8,12 -> zmm1,5,9,13 after computing x.hi52 (in zmm2,6,10,14): */\
			"vaddpd	%%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd	%%zmm5,%%zmm5,%%zmm5	\n\t	vaddpd	%%zmm9 ,%%zmm9 ,%%zmm9	\n\t	vaddpd	%%zmm13,%%zmm13,%%zmm13	\n\t"/* 2*x.lo */\
			"vaddpd	%%zmm2,%%zmm2,%%zmm2		\n\t	vaddpd	%%zmm6,%%zmm6,%%zmm6	\n\t	vaddpd	%%zmm10,%%zmm10,%%zmm10	\n\t	vaddpd	%%zmm14,%%zmm14,%%zmm14	\n\t"/* 2*x.hi */\
			/* And subtract q.lo,hi: */\
			"vsubpd	     (%%rax),%%zmm1,%%zmm1	\n\t vsubpd	0x040(%%rax),%%zmm5,%%zmm5	\n\t vsubpd	0x080(%%rax),%%zmm9 ,%%zmm9 \n\t vsubpd	0x0c0(%%rax),%%zmm13,%%zmm13\n\t"/* -= q.lo */\
			"vsubpd	0x100(%%rax),%%zmm2,%%zmm2	\n\t vsubpd	0x140(%%rax),%%zmm6,%%zmm6	\n\t vsubpd	0x180(%%rax),%%zmm10,%%zmm10\n\t vsubpd	0x1c0(%%rax),%%zmm14,%%zmm14\n\t"/* -= q.hi */\
			/* Final-residue computation needs lo-digit >= 0 normalization: */\
			"vmulpd	%%zmm30,%%zmm1,%%zmm1		\n\t	vmulpd	%%zmm30,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm30,%%zmm9,%%zmm9	\n\t	vmulpd	%%zmm30,%%zmm13,%%zmm13	\n\t"/* xlo *= binv */\
			"vrndscalepd	$1,%%zmm1,%%zmm0	\n\t	vrndscalepd	$1,%%zmm5,%%zmm4	\n\t	vrndscalepd	$1,%%zmm9,%%zmm8	\n\t	vrndscalepd	$1,%%zmm13,%%zmm12	\n\t"/* fcy */\
			"vsubpd	%%zmm0,%%zmm1,%%zmm1		\n\t	vsubpd	%%zmm4,%%zmm5,%%zmm5	\n\t	vsubpd	%%zmm8,%%zmm9,%%zmm9	\n\t	vsubpd	%%zmm12,%%zmm13,%%zmm13	\n\t"/* fx0 -= fcy */\
			"vaddpd	%%zmm0,%%zmm2,%%zmm2		\n\t	vaddpd	%%zmm4,%%zmm6,%%zmm6	\n\t	vaddpd	%%zmm8,%%zmm10,%%zmm10	\n\t	vaddpd	%%zmm12,%%zmm14,%%zmm14	\n\t"/* Add carry into xhi */\
			"vmulpd	%%zmm31,%%zmm1,%%zmm1		\n\t	vmulpd	%%zmm31,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm31,%%zmm9,%%zmm9	\n\t	vmulpd	%%zmm31,%%zmm13,%%zmm13	\n\t"/* fx0 *= base ... result must go into m0 */\
			/* (lo26 == 1 && hi52 == 0), save result in a bitmask: */
			"vmulpd	%%zmm31,%%zmm30,%%zmm3		\n\t"/* 1.0 ... no good way to generate 1.0 on-the-fly, so use base * binv */\
			"vpxorq	     %%zmm0,%%zmm0,%%zmm0	\n\t"/* 0.0 */\
			"vcmppd	$0,%%zmm1,%%zmm3,%%k0		\n\t	vcmppd	$0,%%zmm5,%%zmm3,%%k1	\n\t	vcmppd	$0,%%zmm9 ,%%zmm3,%%k2	\n\t	vcmppd	$0,%%zmm13,%%zmm3,%%k3	\n\t"/* low 26 bits == 1 ? In AVX512 version, mask-regs k0-3 replace dest-regs m1,5,9,13 */\
			"vcmppd	$0,%%zmm2,%%zmm0,%%k4		\n\t	vcmppd	$0,%%zmm6,%%zmm0,%%k5	\n\t	vcmppd	$0,%%zmm10,%%zmm0,%%k6	\n\t	vcmppd	$0,%%zmm14,%%zmm0,%%k7	\n\t"/* top 52 bits == 0 ? In AVX512 version, mask-regs k4-7 replace dest-regs m0,4,8,12 */\
			"kandw	%%k0,%%k4,%%k0				\n\t	kandw	%%k1,%%k5,%%k1			\n\t	kandw	%%k2,%%k6,%%k2			\n\t	kandw	%%k3,%%k7,%%k3			\n\t"/* and the 2 results together into four 8-bit bitmasks */\
			"kmovw	%%k0,%%eax					\n\t	kmovw	%%k1,%%ebx				\n\t	kmovw	%%k2,%%ecx				\n\t	kmovw	%%k3,%%edx				\n\t"\
			/* concatenate into 16-bit result and write: */\
			"shlq	$8 ,%%rbx		\n\t"\
			"shlq	$16,%%rcx		\n\t"\
			"shlq	$24,%%rdx		\n\t"\
			"addq	%%rbx,%%rax		\n\t"\
			"addq	%%rdx,%%rcx		\n\t"\
			"addq	%%rcx,%%rax		\n\t"\
			"movq	%%rax,%[__result]	\n\t"\
			:					/* outputs: none */\
			: [__fq0] "m" (dptr)	/* All inputs from memory addresses here */\
			 ,[__result] "m" (r)	\
			: "cc","memory","cl","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);

		return r;
	}

#endif	// USE_AVX512 ?

