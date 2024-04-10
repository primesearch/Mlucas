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

#include "factor.h"	/* Need this for include of platform-specific defs used in YES_ASM test below */

//*** This file contains specialized implementations of 64-bit modpow for 64-bit-int-capable GPU
//*** [e.g. nVidia CUDA with compute capability >= 2.0] and 64-bit x86:

#undef YES_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
#endif

#ifdef __CUDACC__

	// Simple GPU-ized version of twopmodq64.
	// Return 1 if q = 2.k.p+1 divides 2^p-1, 0 otherwise.
	__device__ uint32 twopmodq64_gpu(			// Sans noinline qualifier, compile time of GPU_TF64_pop64 using this function blows up.
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	)
	{
	#ifdef __CUDA_ARCH__	// Only allow device-code compilation of function body, otherwise get
							// no end of inline-asm-related errors, since nvcc doesn't speak that language
		 int32 j;	// Do 2* on k since p is 32-bit (e.g. (p<<1) = 0 for p = F31):
		uint64 q = p*(k<<1) + 1, qinv, x,y, qhalf = q>>1;	// q>>1 = (q-1)/2, since q odd.
		// q must be odd for Montgomery-style modmul to work:
		uint32 tmp32, q32 = (uint32)q;
		uint32 qinv32 = (q32 + q32 + q32) ^ (uint32)2;	// This gives 4-bit inverse
		// First 3 iterations (4->8, 8->16, 16->32 bits) can all use 32-bit arithmetic:
		for(j = 0; j < 3; j++)
		{
			tmp32 = q32*qinv32;
			qinv32 = qinv32*((uint32)2 - tmp32);
		}
		// Promote to 64-bit and do one more iteration to get a 64-bit inverse:
		qinv = qinv32;
		qinv = qinv*((uint64)2 - q*qinv);

		// Since zstart := 1<<zshift is a power of two < 2^128, use streamlined code sequence for the first iteration:
		j = start_index-1;
		// MULL64(zstart,qinv) simply amounts to a left-shift of the bits of qinv:
		x = qinv << zshift;
	#ifdef MUL_LOHI64_SUBROUTINE
		x = __MULH64(q,x);
	#else
		MULH64(q,x,x);
	#endif
		// hi = 0 in this instance, which simplifies things.
		y = q - x;

		if((pshift >> j) & 1)	// Faster than ANDing with a precomputed length-32 2^n array (n=0,...,31)
		{
		#ifdef NOBRANCH
			x = y + y;	x -= q & -(x >= q || x < y);
		#else
			x = y + y;	if(x >= q || x < y) x -= q;
		#endif
		} else {
			x = y;
		}

		for(j = start_index-2; j >= 0; j--)
		{
			// 3-MUL sequence to effect x = MONT_SQR64(x,q,qinv):
		#ifdef MUL_LOHI64_SUBROUTINE
			SQR_LOHI64(x,&x,&y);	// Lo half of x^2 overwrites x, hi half into y
			MULL64(qinv,x,x);
			x = MULH64(q,x);
		#else
			SQR_LOHI64(x, x, y);
			MULL64(qinv,x,x);
			MULH64(q,x,x);
		#endif
			x = y - x + ((-(int64)(y < x)) & q);	/* did we have a borrow from (y-x)? */

			if((pshift >> j) & 1)
			{
			#ifdef NOBRANCH
				x = x + x - (q & -(x > qhalf));
			#else
				if(x > qhalf)	// Combines overflow-on-add and need-to-subtract-q-from-sum checks
				{
					x = x + x;
					x -= q;
				}
				else
					x = x + x;
			#endif
			}
		}
		// Double and return. This is specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		x = x+x-q+FERMAT;	// In the case of interest, x = (q+1)/2 < 2^63, so x + x cannot overflow.
		return (x==1);
	#else	// ifndef __CUDA_ARCH__
		ASSERT(0, "Device code being called in host mode!");
		return 0;
	#endif
	}

  // Comment out for now - NVCC takes ungodly long to compile the MONT_SQR64_q4 sequence, and GPU seems to like 1-factor-at-a-time better anyway:
  #if 0
	// 4-operand version: Returns results in bits <0:3> of retval. A '1' bit indicates that the corresponding k-value yields a factor.
	__noinline__ __device__ uint32 twopmodq64_q4_GPU(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k0, const uint64 k1, const uint64 k2, const uint64 k3,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	)
	{
	#ifdef __CUDA_ARCH__	// Only allow device-code compilation of function body, otherwise get
							// no end of inline-asm-related errors, since nvcc doesn't speak that language
		uint8 r = 0;
		 int32 j;
		uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3
			, qinv0, qinv1, qinv2, qinv3
			, x0, x1, x2, x3
			, y0, y1, y2, y3;

		qinv0 = (q0+q0+q0) ^ (uint64)2;
		qinv1 = (q1+q1+q1) ^ (uint64)2;
		qinv2 = (q2+q2+q2) ^ (uint64)2;
		qinv3 = (q3+q3+q3) ^ (uint64)2;

		for(j = 0; j < 4; j++)
		{
			qinv0 = qinv0*((uint64)2 - q0*qinv0);
			qinv1 = qinv1*((uint64)2 - q1*qinv1);
			qinv2 = qinv2*((uint64)2 - q2*qinv2);
			qinv3 = qinv3*((uint64)2 - q3*qinv3);
		}

		// Since zstart := 1<<zshift is a power of two < 2^128, use streamlined code sequence for the first iteration:
		j = start_index-1;

		// MULL64(zstart,qinv) simply amounts to a left-shift of the bits of qinv:
		x0 = qinv0 << zshift;
		x1 = qinv1 << zshift;
		x2 = qinv2 << zshift;
		x3 = qinv3 << zshift;

	#ifdef MUL_LOHI64_SUBROUTINE
		x0 = __MULH64(q0,x0);
		x1 = __MULH64(q1,x1);
		x2 = __MULH64(q2,x2);
		x3 = __MULH64(q3,x3);
	#else
		MULH64(q0,x0,x0);
		MULH64(q1,x1,x1);
		MULH64(q2,x2,x2);
		MULH64(q3,x3,x3);
	#endif
		// hi = 0 in this instance, which simplifies things.
		y0 = q0 - x0;
		y1 = q1 - x1;
		y2 = q2 - x2;
		y3 = q3 - x3;

		if((pshift >> j) & 1)
		{
		#ifdef NOBRANCH
			x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
			x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
			x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
			x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
		#else
			x0 = y0 + y0;	if(x0 >= q0 || x0 < y0) x0 -= q0;
			x1 = y1 + y1;	if(x1 >= q1 || x1 < y1) x1 -= q1;
			x2 = y2 + y2;	if(x2 >= q2 || x2 < y2) x2 -= q2;
			x3 = y3 + y3;	if(x3 >= q3 || x3 < y3) x3 -= q3;
		#endif
		}
		else
		{
			x0 = y0;
			x1 = y1;
			x2 = y2;
			x3 = y3;
		}

		for(j = start_index-2; j >= 0; j--)
		{
//	MONT_SQR64_q4(x0,x1,x2,x3,q0,q1,q2,q3,qinv0,qinv1,qinv2,qinv3,y0,y1,y2,y3);
			SQR_LOHI64(x0, x0, y0);								\
			SQR_LOHI64(x1, x1, y1);								\
			SQR_LOHI64(x2, x2, y2);								\
			SQR_LOHI64(x3, x3, y3);								\
			MULL64(qinv0,x0,x0);								\
			MULL64(qinv1,x1,x1);								\
			MULL64(qinv2,x2,x2);								\
			MULL64(qinv3,x3,x3);								\
			MULH64(q0,x0,x0);									\
			MULH64(q1,x1,x1);									\
			MULH64(q2,x2,x2);									\
			MULH64(q3,x3,x3);									\
			/* did we have a borrow from (y-x)? */				\
			y0 = y0 - x0 + ((-(int64)(y0 < x0)) & q0);		\
			y1 = y1 - x1 + ((-(int64)(y1 < x1)) & q1);		\
			y2 = y2 - x2 + ((-(int64)(y2 < x2)) & q2);		\
			y3 = y3 - x3 + ((-(int64)(y3 < x3)) & q3);		\
			if((pshift >> j) & 1)
			{
			#ifdef NOBRANCH
				x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
				x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
				x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
				x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
			#else
				x0 = y0 + y0;	if(x0 >= q0 || x0 < y0) x0 -= q0;
				x1 = y1 + y1;	if(x1 >= q1 || x1 < y1) x1 -= q1;
				x2 = y2 + y2;	if(x2 >= q2 || x2 < y2) x2 -= q2;
				x3 = y3 + y3;	if(x3 >= q3 || x3 < y3) x3 -= q3;
			#endif
			}
			else
			{
				x0 = y0;
				x1 = y1;
				x2 = y2;
				x3 = y3;
			}
		}
		//...Double and return.	Specialize for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		r  = (x0+x0-q0+FERMAT == 1);
		r += (x1+x1-q1+FERMAT == 1);
		r += (x2+x2-q2+FERMAT == 1);
		r += (x3+x3-q3+FERMAT == 1);
		return r;
	#else	// ifndef __CUDA_ARCH__
		ASSERT(0, "Device code being called in host mode!");
		return 0;
	#endif
	}
  #endif 	// #if 0

	__global__ void VecModpow64(const uint64*pvec, const uint64*pshft, const uint32*zshft, const uint32*stidx, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			uint64 p           = pvec[i];
			uint64 pshift      = pshft[i];
			uint32 zshift      = zshft[i];
			uint32 start_index = stidx[i];
			uint64 k           = kvec[i];
			// Result returned in rvec[i]:
			rvec[i] = twopmodq64_gpu(p, pshift, 0, k, start_index, zshift, i);
		}
	}

	__global__ void GPU_TF64_pop64(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64 k0, uint32*kvec, const uint32 N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			uint64 k = k0 + kvec[i];
			// If this k yields a factor, save the k in the 'one-beyond' slot of the kvec array:
			if(twopmodq64_gpu(p, pshift, FERMAT, k, start_index, zshift, i)) {
				kvec[N  ] = (uint32)(k      );	// Use 'one-beyond' elt of kvec, assume at most one factor k per batch (of ~50000), thus no need for atomic update here.
				kvec[N+1] = (uint32)(k >> 32);	// Must break 64-bit k into 2 32-bit pieces since kvec (32-bit offsets) array is 32-bit
			}
		}
	}

	__global__ void GPU_TF64(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			uint64 k = kvec[i];
			// Result returned in rvec[i]:
			rvec[i] = twopmodq64_gpu(p, pshift, FERMAT, k, start_index, zshift, i);
		}
	}

	__global__ void GPU_TF64_q4(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			int i4 = i << 2;	// *= 4 here to yield base-ID for each k-quartet
			uint64 k0 = kvec[i4], k1 = kvec[i4+1], k2 = kvec[i4+2], k3 = kvec[i4+3];
			// Result returned in bits <0:3> of rvec[i]:
		#if 0	// Disable until fix slow-compile issue with the twopmodq78_q4_GPU source code
			rvec[i] = twopmodq64_q4_GPU(p, pshift, FERMAT, k0,k1,k2,k3, start_index, zshift, i);
		#else
			uint8 r0,r1,r2,r3;
			r0 = twopmodq64_gpu(p, pshift, FERMAT, k0, start_index, zshift, i);
			r1 = twopmodq64_gpu(p, pshift, FERMAT, k1, start_index, zshift, i);
			r2 = twopmodq64_gpu(p, pshift, FERMAT, k2, start_index, zshift, i);
			r3 = twopmodq64_gpu(p, pshift, FERMAT, k3, start_index, zshift, i);
			rvec[i] = r0 + (r1<<1) + (r2<<2) + (r3<<3);
		#endif
		}
	}
/*
	__global__ void GPU_TF64_q8(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			int i8 = i << 3;	// *= 8 here to yield base-ID for each k-octtet
			uint64 k0 = kvec[i8], k1 = kvec[i8+1], k2 = kvec[i8+2], k3 = kvec[i8+3], k4 = kvec[i8+4], k5 = kvec[i8+5], k6 = kvec[i8+6], k7 = kvec[i8+7];
			// Result returned in bits <0:7> (that is, all 8 bits) of rvec[i]:
			rvec[i] = twopmodq64_q8_GPU(p, pshift, FERMAT, k0,k1,k2,k3,k4,k5,k6,k7, start_index, zshift, i);
		}
	}
*/
#endif	// __CUDACC__

#ifdef YES_ASM
/*** 4-trial-factor version, optimized x86_64 asm for main powering loop ***/
uint64 twopmodq64_q4(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
//int dbg = ( (p == (1ull<<32)) && ( (k0 == 2958ull) || (k1 == 2958ull) || (k2 == 2958ull) || (k3 == 2958ull) ) );
//if(dbg) printf("Hit! k0-3 = %llu, %llu, %llu, %llu\n",k0, k1, k2, k3);
	int32 j;
	uint64 r = (p<<1), q0 = 1+r*k0, q1 = 1+r*k1, q2 = 1+r*k2, q3 = 1+r*k3
	, qinv0, qinv1, qinv2, qinv3
	, x0, x1, x2, x3
	, y0, y1, y2, y3
	, lo0, lo1, lo2, lo3;
	uint64 pshift;
	uint32 jshift, leadb, start_index, zshift;
	uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	pshift = p + 64;
	jshift = leadz64(pshift);
	// Extract leftmost 6 bits of pshift and subtract from 64:
	leadb = ((pshift<<jshift) >> (64-6));
	start_index = (64-6)-jshift;
	zshift = 63 - leadb;
	zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	pshift = ~pshift;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 , "even modulus!");

	// Compute 64-bit mod-inverses starting with 8-bits-good initializers:
	uint32 q32_0,q32_1,q32_2,q32_3, qi32_0,qi32_1,qi32_2,qi32_3;
	// minv8 = Table of precomputed byte-inverses def'd in mi64.h:
	q32_0  = q0; qi32_0 = minv8[(q0 & 0xff)>>1];
	q32_1  = q1; qi32_1 = minv8[(q1 & 0xff)>>1];
	q32_2  = q2; qi32_2 = minv8[(q2 & 0xff)>>1];
	q32_3  = q3; qi32_3 = minv8[(q3 & 0xff)>>1];
	// Iteration 1 yields 16-bits-good inverses:
	qi32_0 = qi32_0*(2u - q32_0*qi32_0);
	qi32_1 = qi32_1*(2u - q32_1*qi32_1);
	qi32_2 = qi32_2*(2u - q32_2*qi32_2);
	qi32_3 = qi32_3*(2u - q32_3*qi32_3);
	// Iteration 2 yields 32-bits-good inverses, which we store in uint64s to set up for final 64-bits-math iteration:
	qinv0 = qi32_0*(2u - q32_0*qi32_0);
	qinv1 = qi32_1*(2u - q32_1*qi32_1);
	qinv2 = qi32_2*(2u - q32_2*qi32_2);
	qinv3 = qi32_3*(2u - q32_3*qi32_3);
	// Iteration 3 yields 64-bits-good inverses:
	qinv0 = qinv0*(2ull - q0*qinv0);
	qinv1 = qinv1*(2ull - q1*qinv1);
	qinv2 = qinv2*(2ull - q2*qinv2);
	qinv3 = qinv3*(2ull - q3*qinv3);

	/* Since zstart is a power of two < 2^128, use a
	 streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL64(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	lo0 = qinv0 << zshift;
	lo1 = qinv1 << zshift;
	lo2 = qinv2 << zshift;
	lo3 = qinv3 << zshift;
	MULH64(q0,lo0,lo0);
	MULH64(q1,lo1,lo1);
	MULH64(q2,lo2,lo2);
	MULH64(q3,lo3,lo3);

	/* hi = 0 in this instance, which simplifies things. */
	y0 = q0 - lo0;
	y1 = q1 - lo1;
	y2 = q2 - lo2;
	y3 = q3 - lo3;

	if((pshift >> j) & (uint64)1)
	{
		x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
		x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
		x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
		x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
	}
	else
	{
		x0 = y0;
		x1 = y1;
		x2 = y2;
		x3 = y3;
	}
//if(dbg) printf("x0 = %llu\n",x0);

//printf("twopmodq64_q4 : x1 = %s\n", &str0[convert_uint64_base10_char(str0, x1)] );
//for(j = start_index-2; j >= 0; j--) {
	__asm__ volatile (\
	/* Load the q0|1|2 into rbx|rsi|rdi, keeping rcx free for loop index and rax:rdx for double-width MULs: */\
		"movq	%[__q0],%%rbx	\n\t"\
		"movq	%[__q1],%%rsi	\n\t"\
		"movq	%[__q2],%%rdi	\n\t"\
		/* Must load q3 as-needed due to lack of user register to store it */\
	/* Load the x's into r12-15: */\
		"movq	%[__x0],%%r12	\n\t"\
		"movq	%[__x1],%%r13	\n\t"\
		"movq	%[__x2],%%r14	\n\t"\
		"movq	%[__x3],%%r15	\n\t"\
	/* Pure-ASM loop control: for(j = start_index-2; j >= 0; j--) */\
		"movslq	%[__start_index], %%rcx		\n\t"\
		"subq $2,%%rcx						\n\t"\
		"test %%rcx, %%rcx					\n\t"\
		"jl LoopEnd4			\n\t"/* Skip if n < 0 */\
	"LoopBeg4:								\n\t"\
	/* SQR_LOHI_q4(x*, lo*, hi*): */\
		"movq	%%r12,%%rax	\n\t"/* x0-3 in r8-11. */\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r8 	\n\t"/* lo0 */\
		"movq	%%rdx,%%r12	\n\t"/* hi0 */\
		"movq	%%r13,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r9 	\n\t"/* lo1 */\
		"movq	%%rdx,%%r13	\n\t"/* hi1 */\
		"movq	%%r14,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r10	\n\t"/* lo2 */\
		"movq	%%rdx,%%r14	\n\t"/* hi2 */\
		"movq	%%r15,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r11	\n\t"/* lo3 */\
		"movq	%%rdx,%%r15	\n\t"/* hi3 */\
	/* MULL_q4((qinv*, lo*, lo*): */\
		"imulq	%[__qinv0],%%r8 	\n\t"\
		"imulq	%[__qinv1],%%r9 	\n\t"\
		"imulq	%[__qinv2],%%r10	\n\t"\
		"imulq	%[__qinv3],%%r11	\n\t"\
	/* UMULH_q4((q*, lo*, lo*): lo0-3 in r8-11: */\
		"movq	%%r8 ,%%rax	\n\t"/* lo0 */\
		"mulq	%%rbx		\n\t"/* q0 in rbx; q0*lo0 in rax:rdx */\
		"movq	%%rdx,%%r8 	\n\t"/* Discard low 64 bits [rax] */\
		"\n\t"\
		"movq	%%r9 ,%%rax	\n\t"/* lo1 */\
		"mulq	%%rsi		\n\t"/* q1 in rsi; q1*lo1 in rax:rdx */\
		"movq	%%rdx,%%r9 	\n\t"/* Discard low 64 bits [rax] */\
		"\n\t"\
		"movq	%%r10,%%rax	\n\t"/* lo2 */\
		"mulq	%%rdi		\n\t"/* q2 in rdi; q2*lo2 in rax:rdx */\
		"movq	%%rdx,%%r10	\n\t"/* Discard low 64 bits [rax] */\
		"\n\t"\
		"movq	%%r11,%%rax	\n\t"/* lo3 */\
		"mulq	%[__q3]		\n\t"\
		"movq	%%rdx,%%r11	\n\t"/* Discard low 64 bits [rax] */\
		"\n\t"\
	/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\
		"subq	%%r8 ,%%r12			\n\t"/* r12 = (h-l) */\
		"leaq (%%r12,%%rbx),%%r8 	\n\t"/* r8  = (h-l)+q */\
		"cmovcq %%r8 ,%%r12	\n\t"/* if carry = true (i.e. h > l), copy source (r8 = h-l+q) into dest (r12), else leave dest = h-l. */\
		"\n\t"\
		"subq	%%r9 ,%%r13			\n\t"\
		"leaq (%%r13,%%rsi),%%r9 	\n\t"\
		"cmovcq %%r9 ,%%r13	\n\t"\
		"\n\t"\
		"subq	%%r10,%%r14			\n\t"\
		"leaq (%%r14,%%rdi),%%r10	\n\t"\
		"cmovcq %%r10,%%r14	\n\t"\
		"\n\t"\
		"movq	%[__q3],%%rdx	\n\t"/* Re-use q3 several times below, so store in free register */\
		"subq	%%r11,%%r15			\n\t"\
		"leaq (%%r15,%%rdx),%%r11	\n\t"\
		"cmovcq %%r11,%%r15	\n\t"\
		"\n\t"\
	/* If current bit of pshift == 1, double each output modulo q: */\
		/* if((pshift >> j) & (uint64)1) { */\
		"movq	%[__pshift],%%rax	\n\t"/* Need to follow this with load-j-into-ecx if use HLL loop control in debug mode */\
		"shrq	%%cl,%%rax				\n\t"\
		"andq	$0x1,%%rax				\n\t"\
	"je	twopmodq64_q4_pshiftjmp			\n\t"\
		"\n\t"\
		"movq	%%r12,%%r8 	\n\t"/* r8  <- Copy of x */\
		"movq	%%r12,%%rax	\n\t"/* rax <- Copy of x */\
		"leaq (%%r12,%%r12),%%r12	\n\t"/* x+x */\
		"subq	%%rbx,%%r8	\n\t"/* r8  <- x-q */\
		"addq	%%rax,%%r8 	\n\t"/* r8 <- 2x-q */\
		"cmovcq %%r8 ,%%r12	\n\t"/* if carry (i.e. x+x needed modding), copy source (r8 = x+x-q) into dest (r12), else leave dest = x+x. */\
		"\n\t"\
		"movq	%%r13,%%r9 	\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"leaq (%%r13,%%r13),%%r13	\n\t"\
		"subq	%%rsi,%%r9	\n\t"\
		"addq	%%rax,%%r9 	\n\t"\
		"cmovcq %%r9 ,%%r13	\n\t"\
		"\n\t"\
		"movq	%%r14,%%r10	\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"leaq (%%r14,%%r14),%%r14	\n\t"\
		"subq	%%rdi,%%r10	\n\t"\
		"addq	%%rax,%%r10	\n\t"\
		"cmovcq %%r10,%%r14	\n\t"\
		"\n\t"\
		"movq	%%r15,%%r11	\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"leaq (%%r15,%%r15),%%r15	\n\t"\
		"subq	%[__q3],%%r11	\n\t"/* Weird...using rdx in place of __q3 here made timings 1-2 percent  worse. */\
		"addq	%%rax,%%r11	\n\t"\
		"cmovcq %%r11,%%r15	\n\t"\
		"\n\t"\
	"twopmodq64_q4_pshiftjmp:	\n\t"\
		/* } endif((pshift >> j) & (uint64)1) */\
		"subq	$1,%%rcx	\n\t"/* j-- */\
		"cmpq	$0,%%rcx	\n\t"/* compare j vs 0 */\
		"jge	LoopBeg4	\n\t"/* if (j >= 0), Loop */\
	"LoopEnd4:			\n\t"\
		"movq	%%r12,%[__x0]	\n\t"\
		"movq	%%r13,%[__x1]	\n\t"\
		"movq	%%r14,%[__x2]	\n\t"\
		"movq	%%r15,%[__x3]	\n\t"\
		:	/* outputs: none */\
		: [__q0] "m" (q0)	/* All inputs from memory addresses here */\
		 ,[__q1] "m" (q1)	\
		 ,[__q2] "m" (q2)	\
		 ,[__q3] "m" (q3)	\
		 ,[__qinv0] "m" (qinv0)	\
		 ,[__qinv1] "m" (qinv1)	\
		 ,[__qinv2] "m" (qinv2)	\
		 ,[__qinv3] "m" (qinv3)	\
		 ,[__x0] "m" (x0)	\
		 ,[__x1] "m" (x1)	\
		 ,[__x2] "m" (x2)	\
		 ,[__x3] "m" (x3)	\
		 ,[__pshift] "m" (pshift)	\
		 ,[__j] "m" (j)	/* Only need this if debug and explicit loop enabled, but try with/sans for each version of the asm, pivk faster one. */\
		 ,[__start_index] "m" (start_index)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"		/* Clobbered registers */\
		);
//printf("twopmodq64_q4 : x1 = %s\n", &str0[convert_uint64_base10_char(str0, x1)] );
//	}
//printf("twopmodq64_q4 : x1 = %s\n", &str0[convert_uint64_base10_char(str0, x1+x1-q1)] );
//exit(0);
	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
//if(dbg) printf("xout = %llu\n",x0+x0-q0+FERMAT);
	r = 0;
	if(x0+x0-q0+FERMAT == 1) r +=  1;
	if(x1+x1-q1+FERMAT == 1) r +=  2;
	if(x2+x2-q2+FERMAT == 1) r +=  4;
	if(x3+x3-q3+FERMAT == 1) r +=  8;
	return(r);
}


/*** 8-trial-factor version ***/
uint64 twopmodq64_q8(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
	 int32 j;
	uint64 r = (p<<1), q0 = 1+r*k0, q1 = 1+r*k1, q2 = 1+r*k2, q3 = 1+r*k3, q4 = 1+r*k4, q5 = 1+r*k5, q6 = 1+r*k6, q7 = 1+r*k7
		, qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, y0, y1, y2, y3, y4, y5, y6, y7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7;
	uint64 pshift;
	uint32 jshift, leadb, start_index, zshift;
	uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	pshift = p + 64;
	jshift = leadz64(pshift);
	// Extract leftmost 6 bits of pshift and subtract from 64:
	leadb = ((pshift<<jshift) >> (64-6));
	start_index = (64-6)-jshift;
	zshift = 63 - leadb;
	zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	pshift = ~pshift;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 && q4 & 1 && q5 & 1 && q6 & 1 && q7 & 1 , "even modulus!");

	// Compute 64-bit mod-inverses starting with 8-bits-good initializers:
	uint32 q32_0,q32_1,q32_2,q32_3,q32_4,q32_5,q32_6,q32_7, qi32_0,qi32_1,qi32_2,qi32_3,qi32_4,qi32_5,qi32_6,qi32_7;
	// minv8 = Table of precomputed byte-inverses def'd in mi64.h:
	q32_0  = q0; qi32_0 = minv8[(q0 & 0xff)>>1];
	q32_1  = q1; qi32_1 = minv8[(q1 & 0xff)>>1];
	q32_2  = q2; qi32_2 = minv8[(q2 & 0xff)>>1];
	q32_3  = q3; qi32_3 = minv8[(q3 & 0xff)>>1];
	q32_4  = q4; qi32_4 = minv8[(q4 & 0xff)>>1];
	q32_5  = q5; qi32_5 = minv8[(q5 & 0xff)>>1];
	q32_6  = q6; qi32_6 = minv8[(q6 & 0xff)>>1];
	q32_7  = q7; qi32_7 = minv8[(q7 & 0xff)>>1];
	// Iteration 1 yields 16-bits-good inverses:
	qi32_0 = qi32_0*(2u - q32_0*qi32_0);
	qi32_1 = qi32_1*(2u - q32_1*qi32_1);
	qi32_2 = qi32_2*(2u - q32_2*qi32_2);
	qi32_3 = qi32_3*(2u - q32_3*qi32_3);
	qi32_4 = qi32_4*(2u - q32_4*qi32_4);
	qi32_5 = qi32_5*(2u - q32_5*qi32_5);
	qi32_6 = qi32_6*(2u - q32_6*qi32_6);
	qi32_7 = qi32_7*(2u - q32_7*qi32_7);
	// Iteration 2 yields 32-bits-good inverses, which we store in uint64s to set up for final 64-bits-math iteration:
	qinv0 = qi32_0*(2u - q32_0*qi32_0);
	qinv1 = qi32_1*(2u - q32_1*qi32_1);
	qinv2 = qi32_2*(2u - q32_2*qi32_2);
	qinv3 = qi32_3*(2u - q32_3*qi32_3);
	qinv4 = qi32_4*(2u - q32_4*qi32_4);
	qinv5 = qi32_5*(2u - q32_5*qi32_5);
	qinv6 = qi32_6*(2u - q32_6*qi32_6);
	qinv7 = qi32_7*(2u - q32_7*qi32_7);
	// Iteration 3 yields 64-bits-good inverses:
	qinv0 = qinv0*(2ull - q0*qinv0);
	qinv1 = qinv1*(2ull - q1*qinv1);
	qinv2 = qinv2*(2ull - q2*qinv2);
	qinv3 = qinv3*(2ull - q3*qinv3);
	qinv4 = qinv4*(2ull - q4*qinv4);
	qinv5 = qinv5*(2ull - q5*qinv5);
	qinv6 = qinv6*(2ull - q6*qinv6);
	qinv7 = qinv7*(2ull - q7*qinv7);

	/* Since zstart is a power of two < 2^128, use a
	streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL64(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	lo0 = qinv0 << zshift;
	lo1 = qinv1 << zshift;
	lo2 = qinv2 << zshift;
	lo3 = qinv3 << zshift;
	lo4 = qinv4 << zshift;
	lo5 = qinv5 << zshift;
	lo6 = qinv6 << zshift;
	lo7 = qinv7 << zshift;

	MULH64(q0,lo0,lo0);
	MULH64(q1,lo1,lo1);
	MULH64(q2,lo2,lo2);
	MULH64(q3,lo3,lo3);
	MULH64(q4,lo4,lo4);
	MULH64(q5,lo5,lo5);
	MULH64(q6,lo6,lo6);
	MULH64(q7,lo7,lo7);

	/* hi = 0 in this instance, which simplifies things. */
	y0 = q0 - lo0;
	y1 = q1 - lo1;
	y2 = q2 - lo2;
	y3 = q3 - lo3;
	y4 = q4 - lo4;
	y5 = q5 - lo5;
	y6 = q6 - lo6;
	y7 = q7 - lo7;

	if((pshift >> j) & (uint64)1)
	{
		x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
		x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
		x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
		x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
		x4 = y4 + y4;	x4 -= q4 & -(x4 >= q4 || x4 < y4);
		x5 = y5 + y5;	x5 -= q5 & -(x5 >= q5 || x5 < y5);
		x6 = y6 + y6;	x6 -= q6 & -(x6 >= q6 || x6 < y6);
		x7 = y7 + y7;	x7 -= q7 & -(x7 >= q7 || x7 < y7);
	}
	else
	{
		x0 = y0;
		x1 = y1;
		x2 = y2;
		x3 = y3;
		x4 = y4;
		x5 = y5;
		x6 = y6;
		x7 = y7;
	}

	__asm__ volatile (\
	/* Load the x's into r8-15: */\
		"movq	%[__x0],%%r8 	\n\t"\
		"movq	%[__x1],%%r9 	\n\t"\
		"movq	%[__x2],%%r10	\n\t"\
		"movq	%[__x3],%%r11	\n\t"\
		"movq	%[__x4],%%r12	\n\t"\
		"movq	%[__x5],%%r13	\n\t"\
		"movq	%[__x6],%%r14	\n\t"\
		"movq	%[__x7],%%r15	\n\t"\
	/* Pure-ASM loop control: for(j = start_index-2; j >= 0; j--) */\
		"movslq	%[__start_index], %%rcx		\n\t"\
		"subq $2,%%rcx						\n\t"\
		"test %%rcx, %%rcx					\n\t"\
		"jl LoopEnd8		\n\t"/* Skip if n < 0 */\
	"LoopBeg8:								\n\t"\
	/* SQR_LOHI_q4(x*, lo*, hi*): */\
		"movq	%%r8 ,%%rax		\n\t"\
		"mulq	%%rax			\n\t"\
		"movq	%%rax,%%r8		\n\t"/* lo */\
		"movq	%%rdx,%[__x0]	\n\t"/* Move hi back into x */\
		"movq	%%r9 ,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r9  	\n\t"\
		"movq	%%rdx,%[__x1]	\n\t"\
		"movq	%%r10,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r10 	\n\t"\
		"movq	%%rdx,%[__x2]	\n\t"\
		"movq	%%r11,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r11 	\n\t"\
		"movq	%%rdx,%[__x3]	\n\t"\
		"movq	%%r12,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r12 	\n\t"\
		"movq	%%rdx,%[__x4]	\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r13 	\n\t"\
		"movq	%%rdx,%[__x5]	\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r14		\n\t"\
		"movq	%%rdx,%[__x6]	\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r15		\n\t"\
		"movq	%%rdx,%[__x7]	\n\t"\
		"\n\t"\
	/* MULL_q4((qinv*, lo*, lo*): */\
		"\n\t"\
		"imulq	%[__qinv0],%%r8 	\n\t"\
		"imulq	%[__qinv1],%%r9 	\n\t"\
		"imulq	%[__qinv2],%%r10	\n\t"\
		"imulq	%[__qinv3],%%r11	\n\t"\
		"imulq	%[__qinv4],%%r12	\n\t"\
		"imulq	%[__qinv5],%%r13	\n\t"\
		"imulq	%[__qinv6],%%r14	\n\t"\
		"imulq	%[__qinv7],%%r15	\n\t"\
		"\n\t"\
	/* UMULH_q4((q*, lo*, lo*): lo0-3 in r8-11: */\
		"\n\t"\
		"movq	%%r8 ,%%rax	\n\t"/* lo0 */\
		"mulq	%[__q0]		\n\t"/* q0*lo0 in rax:rdx */\
		"movq	%%rdx,%%r8 	\n\t"/* Discard low 64 bits [rax] */\
		"\n\t"\
		"movq	%%r9 ,%%rax	\n\t"\
		"mulq	%[__q1]		\n\t"\
		"movq	%%rdx,%%r9 	\n\t"\
		"\n\t"\
		"movq	%%r10,%%rax	\n\t"\
		"mulq	%[__q2]		\n\t"\
		"movq	%%rdx,%%r10	\n\t"\
		"\n\t"\
		"movq	%%r11,%%rax	\n\t"\
		"mulq	%[__q3]		\n\t"\
		"movq	%%rdx,%%r11	\n\t"\
		"\n\t"\
		"movq	%%r12,%%rax	\n\t"\
		"mulq	%[__q4]		\n\t"\
		"movq	%%rdx,%%r12	\n\t"\
		"\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"mulq	%[__q5]		\n\t"\
		"movq	%%rdx,%%r13	\n\t"\
		"\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"mulq	%[__q6]		\n\t"\
		"movq	%%rdx,%%r14	\n\t"\
		"\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"mulq	%[__q7]		\n\t"\
		"movq	%%rdx,%%r15	\n\t"\
		"\n\t"\
	/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\
		"movq	%[__x0],%%rax	\n\t" /* h */\
		"movq	%[__q0],%%rbx	\n\t" /* h */\
		"subq	%%r8 ,%%rax		\n\t" /* rax = (h-l), sets cf used by cmov */\
		"leaq (%%rax,%%rbx),%%r8 	\n\t" /* r8  = (h-l)+q, no cf set */\
		"cmovncq %%rax,%%r8 	\n\t" /* if carry = false (i.e. h >= l), copy source (rax = h-l) into dest (r8), else leave dest = (h-l+q. */\
		"\n\t"\
		"movq	%[__x1],%%rax	\n\t"\
		"movq	%[__q1],%%rbx	\n\t"\
		"subq	%%r9 ,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r9 	\n\t"\
		"cmovncq %%rax,%%r9 	\n\t"\
		"\n\t"\
		"movq	%[__x2],%%rax	\n\t"\
		"movq	%[__q2],%%rbx	\n\t"\
		"subq	%%r10,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r10	\n\t"\
		"cmovncq %%rax,%%r10	\n\t"\
		"\n\t"\
		"movq	%[__x3],%%rax	\n\t"\
		"movq	%[__q3],%%rbx	\n\t"\
		"subq	%%r11,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r11	\n\t"\
		"cmovncq %%rax,%%r11	\n\t"\
		"\n\t"\
		"movq	%[__x4],%%rax	\n\t"\
		"movq	%[__q4],%%rbx	\n\t"\
		"subq	%%r12,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r12	\n\t"\
		"cmovncq %%rax,%%r12	\n\t"\
		"\n\t"\
		"movq	%[__x5],%%rax	\n\t"\
		"movq	%[__q5],%%rbx	\n\t"\
		"subq	%%r13,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r13	\n\t"\
		"cmovncq %%rax,%%r13	\n\t"\
		"\n\t"\
		"movq	%[__x6],%%rax	\n\t"\
		"movq	%[__q6],%%rbx	\n\t"\
		"subq	%%r14,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r14	\n\t"\
		"cmovncq %%rax,%%r14	\n\t"\
		"\n\t"\
		"movq	%[__x7],%%rax	\n\t"\
		"movq	%[__q7],%%rbx	\n\t"\
		"subq	%%r15,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r15	\n\t"\
		"cmovncq %%rax,%%r15	\n\t"\
		"\n\t"\
	/* If current bit of pshift == 1, double each output modulo q: */\
		/* if((pshift >> j) & (uint64)1) { */\
		"movq	%[__pshift],%%rax	\n\t"/* Need to follow this with load-j-into-ecx if use HLL loop control in debug mode */\
		"shrq	%%cl,%%rax			\n\t"\
		"andq	$0x1,%%rax			\n\t"\
	"je	twopmodq64_q8_pshiftjmp		\n\t"\
		"\n\t"\
		"movq	%[__q0],%%rax  	\n\t"/* rax <- Copy of q */\
		"leaq (%%r8 ,%%r8 ),%%rbx	\n\t"/* rbx <- x+x, no CF set, leave r8 (x) intact */\
		"subq	%%r8 ,%%rax		\n\t"/* rax <- q-x, sets CF */\
		"subq	%%rax,%%r8 		\n\t"/* r8  <- x-(q-x) = x-q+x */\
		"cmovcq %%rbx,%%r8		\n\t" /* if borrow (i.e. x+x did not need modding), copy source (rax = x+x) into dest (r8 ), else leave dest = x+x-q. */\
		"\n\t"\
		"movq	%[__q1],%%rax  	\n\t"\
		"leaq (%%r9 ,%%r9 ),%%rbx	\n\t"\
		"subq	%%r9 ,%%rax		\n\t"\
		"subq	%%rax,%%r9 		\n\t"\
		"cmovcq %%rbx,%%r9 		\n\t"\
		"\n\t"\
		"movq	%[__q2],%%rax  	\n\t"\
		"leaq (%%r10,%%r10),%%rbx	\n\t"\
		"subq	%%r10,%%rax		\n\t"\
		"subq	%%rax,%%r10		\n\t"\
		"cmovcq %%rbx,%%r10 	\n\t"\
		"\n\t"\
		"movq	%[__q3],%%rax  	\n\t"\
		"leaq (%%r11,%%r11),%%rbx	\n\t"\
		"subq	%%r11,%%rax		\n\t"\
		"subq	%%rax,%%r11		\n\t"\
		"cmovcq %%rbx,%%r11		\n\t"\
		"\n\t"\
		"movq	%[__q4],%%rax  	\n\t"\
		"leaq (%%r12,%%r12),%%rbx	\n\t"\
		"subq	%%r12,%%rax		\n\t"\
		"subq	%%rax,%%r12		\n\t"\
		"cmovcq %%rbx,%%r12		\n\t"\
		"\n\t"\
		"movq	%[__q5],%%rax  	\n\t"\
		"leaq (%%r13,%%r13),%%rbx	\n\t"\
		"subq	%%r13,%%rax		\n\t"\
		"subq	%%rax,%%r13		\n\t"\
		"cmovcq %%rbx,%%r13		\n\t"\
		"\n\t"\
		"movq	%[__q6],%%rax  	\n\t"\
		"leaq (%%r14,%%r14),%%rbx	\n\t"\
		"subq	%%r14,%%rax		\n\t"\
		"subq	%%rax,%%r14		\n\t"\
		"cmovcq %%rbx,%%r14		\n\t"\
		"\n\t"\
		"movq	%[__q7],%%rax  	\n\t"\
		"leaq (%%r15,%%r15),%%rbx	\n\t"\
		"subq	%%r15,%%rax		\n\t"\
		"subq	%%rax,%%r15		\n\t"\
		"cmovcq %%rbx,%%r15		\n\t"\
		"\n\t"\
	"twopmodq64_q8_pshiftjmp:	\n\t"\
		/* } endif((pshift >> j) & (uint64)1) */\
		"subq	$1,%%rcx	\n\t"/* j-- */\
		"cmpq	$0,%%rcx	\n\t"/* compare j vs 0 */\
		"jge	LoopBeg8	\n\t"/* if (j >= 0), Loop */\
	"LoopEnd8:			\n\t"\
		"movq	%%r8 ,%[__x0]	\n\t"\
		"movq	%%r9 ,%[__x1]	\n\t"\
		"movq	%%r10,%[__x2]	\n\t"\
		"movq	%%r11,%[__x3]	\n\t"\
		"movq	%%r12,%[__x4]	\n\t"\
		"movq	%%r13,%[__x5]	\n\t"\
		"movq	%%r14,%[__x6]	\n\t"\
		"movq	%%r15,%[__x7]	\n\t"\
		:	/* outputs: none */\
		: [__q0] "m" (q0)	/* All inputs from memory addresses here */\
		 ,[__q1] "m" (q1)	\
		 ,[__q2] "m" (q2)	\
		 ,[__q3] "m" (q3)	\
		 ,[__q4] "m" (q4)	\
		 ,[__q5] "m" (q5)	\
		 ,[__q6] "m" (q6)	\
		 ,[__q7] "m" (q7)	\
		 ,[__qinv0] "m" (qinv0)	\
		 ,[__qinv1] "m" (qinv1)	\
		 ,[__qinv2] "m" (qinv2)	\
		 ,[__qinv3] "m" (qinv3)	\
		 ,[__qinv4] "m" (qinv4)	\
		 ,[__qinv5] "m" (qinv5)	\
		 ,[__qinv6] "m" (qinv6)	\
		 ,[__qinv7] "m" (qinv7)	\
		 ,[__x0] "m" (x0)	\
		 ,[__x1] "m" (x1)	\
		 ,[__x2] "m" (x2)	\
		 ,[__x3] "m" (x3)	\
		 ,[__x4] "m" (x4)	\
		 ,[__x5] "m" (x5)	\
		 ,[__x6] "m" (x6)	\
		 ,[__x7] "m" (x7)	\
		 ,[__pshift] "m" (pshift)	\
		 ,[__j] "m" (j)	/* Only need this if debug and explicit loop enabled, but try with/sans for each version of the asm, pivk faster one. */\
		 ,[__start_index] "m" (start_index)	\
		: "cc","memory","rax","rbx","rcx","rdx","r8","r9","r10","r11","r12","r13","r14","r15"		/* Clobbered registers */\
		);

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/

	r = 0;
	if(x0+x0-q0+FERMAT == 1) r +=  1;
	if(x1+x1-q1+FERMAT == 1) r +=  2;
	if(x2+x2-q2+FERMAT == 1) r +=  4;
	if(x3+x3-q3+FERMAT == 1) r +=  8;
	if(x4+x4-q4+FERMAT == 1) r += 16;
	if(x5+x5-q5+FERMAT == 1) r += 32;
	if(x6+x6-q6+FERMAT == 1) r += 64;
	if(x7+x7-q7+FERMAT == 1) r +=128;
	return(r);
}

#endif /* YES_ASM */

