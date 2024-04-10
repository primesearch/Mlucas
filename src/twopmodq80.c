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

#include "factor.h"
#include "twopmodq80.h"

// Set = 1 to use streamlined version of modmul sequence, leaving high 52 bits in one double whenever feasible
#define TRUNC_MUL78	1

//#define DBG_ASSERT ASSERT

#ifdef USE_SSE2
	#include "align.h"

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

#endif

#define FAC_DEBUG	0
#if FAC_DEBUG
	char char_buf[1024], str0[64], str1[64];
#endif


/***********************************************************************************/
/*** 78-BIT INPUTS, USING FLOATING-POINT-DOUBLE ARITHMETIC FOR MODMUL **************/
/***********************************************************************************/
/*
Function to find 2^(-p) mod q, where p and q are 64 and 78-bit unsigned integers,
respectively, with q stored in a uint96. Uses a floating-double-based Montgomery-style
modmul with a power-of-2 modulus = TWO26FLOAT^3 (i.e. our MODQ operation effects multiply
modulo TWO26FLOAT^3). Using balanced-digit representation in the floating-point form the
largest base we can use with our 3-word-input convolution algorithm is TWO26FLOAT = 2^26,
with the MSW of each input being non-balanced, i.e. strictly nonnegative.

The key 3-operation sequence here is as follows:

	SQR_LOHI78(x,lo,hi);	// Input   x has 78 bits; outputs lo & hi have 78 bits
	MULL78(lo,qinv,lo);		// Inputs lo & qinv, and output (overwrites lo) have 78 bits
	MULH78(q,lo,lo);		// Inputs  q &   lo, and output (overwrites lo) have 78 bits
*/

#ifdef __CUDACC__

	// Simple GPU-ized version of twopmodq78_3WORD_DOUBLE (sans the experimental pshift = p+96 - here use the standard p+78).
	// Return value: 1 if q = 2.k.p+1 divides 2^p-1, 0 otherwise.
	__device__ uint32 twopmodq78_gpu(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	)
	{
	#ifdef __CUDA_ARCH__	// Only allow device-code compilation of function body, otherwise get
							// no end of inline-asm-related errors, since nvcc doesn't speak that language
	#if FAC_DEBUG
		int dbg = (p == 16727479) && (k == 7946076362870052);
	#endif
		 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
		uint32 q32,qinv32,tmp32;
		uint64 hi64;
		uint96 q, qhalf, qinv, x, lo;
		uint192 prod192;
		// host variables - in this case the globals with these names - cannot be directly read in a device function:
		const uint96 ONE96  = {(uint64)1, (uint32)0};
		const double TWO26FLOAT = (double)0x04000000, TWO26FLINV = 1.0/TWO26FLOAT,
					TWO52FLOAT = TWO26FLOAT*TWO26FLOAT, TWO52FLINV = 1.0/TWO52FLOAT,
					TWO64FLOAT = (double)4.0*0x80000000*0x80000000;
		double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2,fqhalf, fx0,fx1,fx2, flo0,flo1,flo2,flohi52, fhi0;
	#if TRUNC_MUL78	// Streamlined version of modmul sequence, leaving high 52 bits in one double whenever feasible
		double fq_or_nil_lo26[2], fq_or_nil_hi52[2];
	#else
		double fhi1,fhi2;
	#endif
		int fidx;
	#if FAC_DEBUG
		if(dbg) {
			printf("twopmodq78_3WORD_DOUBLE with p = %u, k = %llu, tid = %u\n",p,k,i);
		}
	#endif
/*
if(k == 7946076362870052)printf("In twopmodq78_3WORD_DOUBLE with i = %u, p = %u, k = %llu\n",i,p,k);
*/
		q.d0 = p+p;
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(q.d0, k,&q.d0,&q.d1);	// Common HLL C-code IMUL macro shared by both device and host-code compile passes
	#else
		MUL_LOHI64(q.d0, k, q.d0, q.d1);	// Common HLL C-code IMUL macro shared by both device and host-code compile passes
	#endif
		q.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */

		/* Convert q to floating form: */
		CVT_UINT78_3WORD_DOUBLE(q, fq0,fq1,fq2);
	#if FAC_DEBUG
		if(dbg) {
			printf("fq = %20.15f, %20.15f, %20.15f\n",fq0,fq1,fq2);
		}
	#endif
	#if TRUNC_MUL78	// Streamlined version of modmul sequence, leaving high 52 bits in one double whenever feasible
		fq_or_nil_lo26[0] = 0.0; fq_or_nil_lo26[1] = fq0; fq_or_nil_hi52[0] = 0.0; fq_or_nil_hi52[1] = fq1 + TWO26FLOAT*fq2;
	#endif

		RSHIFT_FAST96(q, 1, qhalf);	/* = (q-1)/2, since q odd. */
		fqhalf = qhalf.d1*TWO64FLOAT + qhalf.d0;

		q32 = (uint32)q.d0;
		qinv32 = (q32 + q32 + q32) ^ (uint32)2;	// This gives 4-bit inverse
		// First 3 iterations (4->8, 8->16, 16->32 bits) can all use 32-bit arithmetic:
		for(j = 0; j < 3; j++)
		{
			tmp32 = q32*qinv32;
			qinv32 = qinv32*((uint32)2 - tmp32);
		}
		// Promote to 64-bit and do one more iteration to get a 64-bit inverse:
		qinv.d0 = qinv32;
		hi64 = q.d0*qinv.d0;
		qinv.d0 = qinv.d0*((uint64)2 - hi64);
		// 64->128-bit inverse, then truncate to the needed 78-bit:
		MULH64(q.d0, qinv.d0, hi64);
		qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + hi64);
		qinv.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
//	if(i == 0)printf("In twopmodq78_gpu with p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,k,zshift,start_index);
		/* Convert qinv to floating form: */
		CVT_UINT78_3WORD_DOUBLE(qinv, fqinv0,fqinv1,fqinv2);
	#if FAC_DEBUG
		if(dbg) {
			printf("fqinv = %20.15f, %20.15f, %20.15f\n",fqinv0,fqinv1,fqinv2);
		}
	#endif

		/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
		j = start_index-1;
		/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
		LSHIFT96(qinv, zshift, lo);
		lo.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */

		MUL_LOHI96_PROD192(q,lo,prod192);
		RSHIFT192(prod192,78,prod192);
		lo.d0 = prod192.d0;	lo.d1 = prod192.d1;
		/* hi = 0 in this instance, which simplifies things. */
		SUB96(q, lo, x);	/* Put the result in lo (rather than x), to ease overflow check below */

		if((pshift >> j) & 1)
		{
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x, qhalf)){ ADD96(x, x, x); SUB96(x, q, x); }else{ ADD96(x, x, x); }
		}
		/* Convert x to floating form: */
		CVT_UINT78_3WORD_DOUBLE(x, fx0,fx1,fx2);

		for(j = start_index-2; j >= 0; j--)
		{
		/********************************************************************************************************/
		#if TRUNC_MUL78	// Streamlined version of modmul sequence, leaving high 52 bits in one double whenever feasible
		/********************************************************************************************************/

			/* SQR_LOHI96(x,lo,hi); */
			SQR_LOHI78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED(fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fx1);
			/* MULL96(lo,qinv,lo); */
			MULL78_3WORD_DOUBLE(flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2);
			/* MULH78(q,lo,lo); */
//if(i==0)printf("MULH Ins: x0 = %f; x1 = %f; x2 = %f\n",fq0,fq1,fq2);
//if(i==0)printf("MULH Ins: y0 = %f; y1 = %f; y2 = %f\n",flo0,flo1,flo2);
			MULH78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flohi52);
//if(i==0)printf("MULH Out: z0 = %f; z12 = %f\n",flo0,flohi52);
#if 0
x0 = 15135097; x1 = 6567697; x2 = 59027372; x=x0+a*x1+b*x2
y0 = 22129433; y1 = -21457855; y2 = 15218570; y=y0+a*y1+b*y2
z0 = 6272576; z12 = 898312175313603; z=z0+a*z12	<*** z0 is +1 too large ***
#endif
			/* If h < l, then calculate q+(h-l) < q; otherwise calculate h-l. */
			fidx  = (fx1 < flohi52);
			fx1  -= flohi52;
			fhi0 -= flo0;
			fx1  +=        fq_or_nil_hi52[fidx];
			fx0   = fhi0 + fq_or_nil_lo26[fidx];
		/*
			// Normalize result and store in x-vector...use hi0 term as carry:
			fhi0  = DNINT(fx0*TWO26FLINV);
			fx0  -= fhi0*TWO26FLOAT;
			fx1  += fhi0;
		*/
			/* Branchless version of the double-and-subtract-q-if-(2x >= q) sequence:
				x = x + x - ((-(x > qhalf)) & q);
			In FP version replace integer and with array lookup.
			*/
		  #if 1	// On my GTX430 this with-branch sequence decidedly faster than the branchless variant below
			if((pshift >> j) & 1)
			{
				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				/* Fill FP approximation to x with as many SBs as possible prior to comparing */
				/* TRUE here means will need to subtract q from 2x: */
				fidx = (fqhalf < (fx1*TWO26FLOAT + fx0));
				fx1 *= 2;	fx0 += fx0;
				fx1 -= fq_or_nil_hi52[fidx];
				fx0 -= fq_or_nil_lo26[fidx];
			}
		  #else
				int mul2 = (pshift >> j) & 1;
				fidx = mul2 & (fqhalf < (fx1*TWO26FLOAT + fx0));
				fx0 *= mul2+1;
				fx1 *= mul2+1;
				fx1 -= fq_or_nil_hi52[fidx];
				fx0 -= fq_or_nil_lo26[fidx];
		  #endif
			/* Normalize the result - use currently-unused x2 coeff as carry: */
			/* Digit 0: */
			fx2  = DNINT(fx0*TWO26FLINV);
			fx0 -= fx2*TWO26FLOAT;
			/* Digit 1: */
			fx1 += fx2;
			fx2  = DNINT(fx1*TWO26FLINV);
			fx1 -= fx2*TWO26FLOAT;
			/* Digit 2 already in x2 term. */

		/********************************************************************************************************/
		#else	// Basic impl of 3-word-double-based 78-bit-integer modmul:
		/********************************************************************************************************/

			/*...x^2 mod q is returned in x. */
			/* SQR_LOHI96(x,lo,hi); */
			SQR_LOHI78_3WORD_DOUBLE(fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2);
			/* MULL96(lo,qinv,lo); */
			MULL78_3WORD_DOUBLE(flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2);
			/* MULH78(q,lo,lo); */
			MULH78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);

			/* If h < l, then calculate q+(h-l) < q; otherwise calculate h-l. */
			if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
			{
				SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
				ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
			}
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

			if((pshift >> j) & 1)
			{
				/* Convert fx to uint96 for purpose of this particular comparison: */
				CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x, qhalf))
				{
					ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
					SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				}
				NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);
			}

		#endif	/* #ifdef TRUNC_MUL78 */
		}
		/*...Double and return.	These are specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		*/
		CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);
		ADD96(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
		SUB96(x,q,x);
	#if FAC_DEBUG
		if(dbg) {
			printf("k = %llu: X_out = %u*2^64 + %llu\n", x.d1,x.d0);
		}
	#endif
		return (CMPEQ96(x, ONE96));
	#else	// ifndef __CUDA_ARCH__
		ASSERT(0, "Device code being called in host mode!");
		return 0;
	#endif
	}

	// 4-operand version: Returns results in bits <0:3> of retval. A '1' bit indicates that the corresponding k-value yields a factor.
	__device__ uint32 twopmodq78_q4_GPU(
		const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint64 k0, const uint64 k1, const uint64 k2, const uint64 k3,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	)
	{
	#ifdef __CUDA_ARCH__	// Only allow device-code compilation of function body, otherwise get
							// no end of inline-asm-related errors, since nvcc doesn't speak that language
		 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
		uint32 q32_0, qinv32_0, tmp32_0
			 , q32_1, qinv32_1, tmp32_1
			 , q32_2, qinv32_2, tmp32_2
			 , q32_3, qinv32_3, tmp32_3;
		uint64 tmp0, tmp1, tmp2, tmp3, r;
		uint96 q0, qinv0, qhalf0, x0, lo0
			 , q1, qinv1, qhalf1, x1, lo1
			 , q2, qinv2, qhalf2, x2, lo2
			 , q3, qinv3, qhalf3, x3, lo3;
		uint192 prod192;
		const uint96 ONE96  = {(uint64)1, (uint32)0};
		const double TWO26FLOAT = (double)0x04000000, TWO26FLINV = 1.0/TWO26FLOAT;
		const double TWO64FLOAT = (double)4.0*0x80000000*0x80000000;
		double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2,fqhalf0, fx0,fx1,fx2, flo0,flo1,flo2,flohi52, fhi0;
		double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2,fqhalf1, gx0,gx1,gx2, glo0,glo1,glo2,glohi52, ghi0;
		double hq0,hq1,hq2, hqinv0,hqinv1,hqinv2,fqhalf2, hx0,hx1,hx2, hlo0,hlo1,hlo2,hlohi52, hhi0;
		double iq0,iq1,iq2, iqinv0,iqinv1,iqinv2,fqhalf3, ix0,ix1,ix2, ilo0,ilo1,ilo2,ilohi52, ihi0;
		double fq_or_nil_lo26[2], fq_or_nil_hi52[2];
		double gq_or_nil_lo26[2], gq_or_nil_hi52[2];
		double hq_or_nil_lo26[2], hq_or_nil_hi52[2];
		double iq_or_nil_lo26[2], iq_or_nil_hi52[2];
		int fidx,gidx,hidx,iidx;

		q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;
	#ifdef MUL_LOHI64_SUBROUTINE
		// MUL_LOHI64 expects a 64-bit high-part pointer, in 32bit builds this buggers us if we try dumping hi-part directly into 32-bit q.d1
		MUL_LOHI64(q0.d0, k0,&q0.d0,&tmp0);	q0.d1 = tmp0;
		MUL_LOHI64(q1.d0, k1,&q1.d0,&tmp0);	q1.d1 = tmp0;
		MUL_LOHI64(q2.d0, k2,&q2.d0,&tmp0);	q2.d1 = tmp0;
		MUL_LOHI64(q3.d0, k3,&q3.d0,&tmp0);	q3.d1 = tmp0;
	#else
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
		MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
		MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);
	#endif
		q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
		q1.d0 += 1;
		q2.d0 += 1;
		q3.d0 += 1;

		q32_0 = (uint32)q0.d0;
		q32_1 = (uint32)q1.d0;
		q32_2 = (uint32)q2.d0;
		q32_3 = (uint32)q3.d0;

		/* Convert q to floating form: */
		CVT_UINT78_3WORD_DOUBLE(q0, fq0,fq1,fq2);
		CVT_UINT78_3WORD_DOUBLE(q1, gq0,gq1,gq2);
		CVT_UINT78_3WORD_DOUBLE(q2, hq0,hq1,hq2);
		CVT_UINT78_3WORD_DOUBLE(q3, iq0,iq1,iq2);

		fq_or_nil_lo26[0] = 0.0; fq_or_nil_lo26[1] = fq0; fq_or_nil_hi52[0] = 0.0; fq_or_nil_hi52[1] = fq1 + TWO26FLOAT*fq2;
		gq_or_nil_lo26[0] = 0.0; gq_or_nil_lo26[1] = gq0; gq_or_nil_hi52[0] = 0.0; gq_or_nil_hi52[1] = gq1 + TWO26FLOAT*gq2;
		hq_or_nil_lo26[0] = 0.0; hq_or_nil_lo26[1] = hq0; hq_or_nil_hi52[0] = 0.0; hq_or_nil_hi52[1] = hq1 + TWO26FLOAT*hq2;
		iq_or_nil_lo26[0] = 0.0; iq_or_nil_lo26[1] = iq0; iq_or_nil_hi52[0] = 0.0; iq_or_nil_hi52[1] = iq1 + TWO26FLOAT*iq2;

		RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
		RSHIFT_FAST96(q1, 1, qhalf1);
		RSHIFT_FAST96(q2, 1, qhalf2);
		RSHIFT_FAST96(q3, 1, qhalf3);

		fqhalf0 = qhalf0.d1*TWO64FLOAT + qhalf0.d0;
		fqhalf1 = qhalf1.d1*TWO64FLOAT + qhalf1.d0;
		fqhalf2 = qhalf2.d1*TWO64FLOAT + qhalf2.d0;
		fqhalf3 = qhalf3.d1*TWO64FLOAT + qhalf3.d0;

		// This gives 4-bit inverse:
		qinv32_0 = (q32_0 + q32_0 + q32_0) ^ (uint32)2;
		qinv32_1 = (q32_1 + q32_1 + q32_1) ^ (uint32)2;
		qinv32_2 = (q32_2 + q32_2 + q32_2) ^ (uint32)2;
		qinv32_3 = (q32_3 + q32_3 + q32_3) ^ (uint32)2;
		// First 3 iterations (4->8, 8->16, 16->32 bits) can all use 32-bit arithmetic:
		for(j = 0; j < 3; j++)
		{
			tmp32_0 = q32_0*qinv32_0;
			tmp32_1 = q32_1*qinv32_1;
			tmp32_2 = q32_2*qinv32_2;
			tmp32_3 = q32_3*qinv32_3;
			qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);
			qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);
			qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);
			qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);
		}
		// Promote to 64-bit and do one more iteration to get a 64-bit inverse:
		qinv0.d0 = (uint64)qinv32_0;
		qinv1.d0 = (uint64)qinv32_1;
		qinv2.d0 = (uint64)qinv32_2;
		qinv3.d0 = (uint64)qinv32_3;
		tmp0 = q0.d0*qinv0.d0;
		tmp1 = q1.d0*qinv1.d0;
		tmp2 = q2.d0*qinv2.d0;
		tmp3 = q3.d0*qinv3.d0;
		qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
		qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
		qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);
		qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);
		// 64->128-bit inverse, then truncate to the needed 78-bit:
	#ifdef MUL_LOHI64_SUBROUTINE
		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
		qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + __MULH64(q2.d0, qinv2.d0));
		qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + __MULH64(q3.d0, qinv3.d0));
	#else
		MULH64(q0.d0, qinv0.d0, tmp0);
		MULH64(q1.d0, qinv1.d0, tmp1);
		MULH64(q2.d0, qinv2.d0, tmp2);
		MULH64(q3.d0, qinv3.d0, tmp3);

		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);
		qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + tmp2);
		qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + tmp3);
	#endif
		/* 64 bits: */
		qinv0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
		qinv1.d1 &= 0x0000000000003fff;
		qinv2.d1 &= 0x0000000000003fff;
		qinv3.d1 &= 0x0000000000003fff;

		/* Convert qinv to floating form: */
		CVT_UINT78_3WORD_DOUBLE(qinv0, fqinv0,fqinv1,fqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv1, gqinv0,gqinv1,gqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv2, hqinv0,hqinv1,hqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv3, iqinv0,iqinv1,iqinv2);

		/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
		j = start_index-1;

		/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
		LSHIFT96(qinv0, zshift, lo0);	lo0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
		LSHIFT96(qinv1, zshift, lo1);	lo1.d1 &= 0x0000000000003fff;
		LSHIFT96(qinv2, zshift, lo2);	lo2.d1 &= 0x0000000000003fff;
		LSHIFT96(qinv3, zshift, lo3);	lo3.d1 &= 0x0000000000003fff;

		MUL_LOHI96_PROD192(q0,lo0,prod192);	RSHIFT192(prod192,78,prod192);	lo0.d0 = prod192.d0;	lo0.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q1,lo1,prod192);	RSHIFT192(prod192,78,prod192);	lo1.d0 = prod192.d0;	lo1.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q2,lo2,prod192);	RSHIFT192(prod192,78,prod192);	lo2.d0 = prod192.d0;	lo2.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q3,lo3,prod192);	RSHIFT192(prod192,78,prod192);	lo3.d0 = prod192.d0;	lo3.d1 = prod192.d1;

		/* hi = 0 in this instance, which simplifies things. */
		SUB96(q0, lo0, x0);	/* Put the result in lo (rather than x), to ease overflow check below */
		SUB96(q1, lo1, x1);
		SUB96(q2, lo2, x2);
		SUB96(q3, lo3, x3);

		if((pshift >> j) & (uint32)1)
		{
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
			if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
			if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
			if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
		}

		/* Convert x to floating form: */
		CVT_UINT78_3WORD_DOUBLE(x0, fx0,fx1,fx2);
		CVT_UINT78_3WORD_DOUBLE(x1, gx0,gx1,gx2);
		CVT_UINT78_3WORD_DOUBLE(x2, hx0,hx1,hx2);
		CVT_UINT78_3WORD_DOUBLE(x3, ix0,ix1,ix2);

		/*...x^2 mod q is returned in x. */
		/* All 3-word-double-form operands have components in the following size ranges:
			fword0,1 in [-2^25, +2^25]
			fword2   in [   -1, +2^26]
		*/
		for(j = start_index-2; j >= 0; j--)
		{
			/*...x^2 mod q is returned in x. */

		/*** Only provide TRUNC_MUL78-style streamlined version of modmul sequence in 4-q version ***/

			/* SQR_LOHI96(x,lo,hi); */
			SQR_LOHI78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(
				  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fx1
				, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,gx1
				, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hx1
				, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ix1
			);
			/* MULL96(lo,qinv,lo); */
			MULL78_3WORD_DOUBLE_q4(
				  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
				, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
				, hlo0,hlo1,hlo2, hqinv0,hqinv1,hqinv2, hlo0,hlo1,hlo2
				, ilo0,ilo1,ilo2, iqinv0,iqinv1,iqinv2, ilo0,ilo1,ilo2
			);
			/* MULH96(q,lo,lo); */
			MULH78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(
				  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flohi52
				, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glohi52
				, hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlohi52
				, iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilohi52
			);

			/* If h < l, then calculate q+(h-l) < q; otherwise calculate h-l. */
			fidx = (fx1 < flohi52);
			gidx = (gx1 < glohi52);
			hidx = (hx1 < hlohi52);
			iidx = (ix1 < ilohi52);

			fx1 -= flohi52;
			gx1 -= glohi52;
			hx1 -= hlohi52;
			ix1 -= ilohi52;

			fhi0 -= flo0;
			ghi0 -= glo0;
			hhi0 -= hlo0;
			ihi0 -= ilo0;

			fx1 += fq_or_nil_hi52[fidx];
			gx1 += gq_or_nil_hi52[gidx];
			hx1 += hq_or_nil_hi52[hidx];
			ix1 += iq_or_nil_hi52[iidx];

			fx0 = fhi0 + fq_or_nil_lo26[fidx];
			gx0 = ghi0 + gq_or_nil_lo26[gidx];
			hx0 = hhi0 + hq_or_nil_lo26[hidx];
			ix0 = ihi0 + iq_or_nil_lo26[iidx];
		/*
			// Normalize result and store in x-vector...use hi0 terms as carries:
			fhi0 = DNINT(fx0*TWO26FLINV);
			ghi0 = DNINT(gx0*TWO26FLINV);
			hhi0 = DNINT(hx0*TWO26FLINV);
			ihi0 = DNINT(ix0*TWO26FLINV);

			fx0 -= fhi0*TWO26FLOAT;
			gx0 -= ghi0*TWO26FLOAT;
			hx0 -= hhi0*TWO26FLOAT;
			ix0 -= ihi0*TWO26FLOAT;

			fx1 += fhi0;
			gx1 += ghi0;
			hx1 += hhi0;
			ix1 += ihi0;
		*/
			/* Branchless version of the double-and-subtract-q-if-(2x >= q) sequence:
				x = x + x - ((-(x > qhalf)) & q);
			In FP version replace integer and with array lookup.
			*/
			if((pshift >> j) & 1)
			{
				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				/* Fill FP approximation to x with as many SBs as possible prior to comparing */
				/* TRUE here means will need to subtract q from 2x: */
				fidx = (fqhalf0 < fx1*TWO26FLOAT + fx0);
				gidx = (fqhalf1 < gx1*TWO26FLOAT + gx0);
				hidx = (fqhalf2 < hx1*TWO26FLOAT + hx0);
				iidx = (fqhalf3 < ix1*TWO26FLOAT + ix0);

				fx1 *= 2;	fx0 += fx0;
				gx1 *= 2;	gx0 += gx0;
				hx1 *= 2;	hx0 += hx0;
				ix1 *= 2;	ix0 += ix0;

				fx1 -= fq_or_nil_hi52[fidx];
				gx1 -= gq_or_nil_hi52[gidx];
				hx1 -= hq_or_nil_hi52[hidx];
				ix1 -= iq_or_nil_hi52[iidx];

				fx0     -= fq_or_nil_lo26[fidx];
				gx0     -= gq_or_nil_lo26[gidx];
				hx0     -= hq_or_nil_lo26[hidx];
				ix0     -= iq_or_nil_lo26[iidx];
			}

			/* Normalize the result - use currently-unused x2 coeff as carry: */
			/* Digit 0: */
			fx2 = DNINT(fx0*TWO26FLINV);
			gx2 = DNINT(gx0*TWO26FLINV);
			hx2 = DNINT(hx0*TWO26FLINV);
			ix2 = DNINT(ix0*TWO26FLINV);

			fx0 -= fx2*TWO26FLOAT;
			gx0 -= gx2*TWO26FLOAT;
			hx0 -= hx2*TWO26FLOAT;
			ix0 -= ix2*TWO26FLOAT;

			/* Digit 1: */
			fx1 += fx2;
			gx1 += gx2;
			hx1 += hx2;
			ix1 += ix2;

			fx2 = DNINT(fx1*TWO26FLINV);
			gx2 = DNINT(gx1*TWO26FLINV);
			hx2 = DNINT(hx1*TWO26FLINV);
			ix2 = DNINT(ix1*TWO26FLINV);

			fx1 -= fx2*TWO26FLOAT;
			gx1 -= gx2*TWO26FLOAT;
			hx1 -= hx2*TWO26FLOAT;
			ix1 -= ix2*TWO26FLOAT;

			/* Digit 2 already in x2 term. */
		}

		/*...Double and return.	These are specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		*/
		CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
		CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
		CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
		CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);

		ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
		ADD96(x1,x1,x1);
		ADD96(x2,x2,x2);
		ADD96(x3,x3,x3);

		SUB96(x0,q0,x0);
		SUB96(x1,q1,x1);
		SUB96(x2,q2,x2);
		SUB96(x3,q3,x3);

		tmp0 = CMPEQ96(x0, ONE96);
		tmp1 = CMPEQ96(x1, ONE96);
		tmp2 = CMPEQ96(x2, ONE96);
		tmp3 = CMPEQ96(x3, ONE96);
		r = tmp0;
		r += tmp1 << 1;
		r += tmp2 << 2;
		r += tmp3 << 3;
		return r;
	#else	// ifndef __CUDA_ARCH__
		ASSERT(0, "Device code being called in host mode!");
		return 0;
	#endif
	}

	__global__ void VecModpow78(const uint64*pvec, const uint64*pshft, const uint32*zshft, const uint32*stidx, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			uint64 p           = pvec[i];
			uint64 pshift      = pshft[i];
			uint32 zshift      = zshft[i];
			uint32 start_index = stidx[i];
			uint64 k           = kvec[i];
			// Result returned in rvec[i]:
			rvec[i] = twopmodq78_gpu(p, pshift, 0, k, start_index, zshift, i);
		}
	}

	__global__ void GPU_TF78_pop64(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64 k0, uint32*kvec, const uint32 N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			uint64 k = k0 + kvec[i];
			// If this k yields a factor, save the k in the 'one-beyond' slot of the kvec array:
			if(twopmodq78_gpu(p, pshift, FERMAT, k, start_index, zshift, i)) {
				kvec[N] = k;	// Use 'one-beyond' elt of kvec, assume at most one factor k per batch (of ~50000), thus no need for atomic update here.
			}
		}
	}

	__global__ void GPU_TF78(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			uint64 k           = kvec[i];
			// Result returned in rvec[i]:
			rvec[i] = twopmodq78_gpu(p, pshift, FERMAT, k, start_index, zshift, i);
		}
	}

	__global__ void GPU_TF78_q4(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			int i4 = i << 2;	// *= 4 here to yield base-ID for each k-quartet
			uint64 k0 = kvec[i4], k1 = kvec[i4+1], k2 = kvec[i4+2], k3 = kvec[i4+3];
			// Result returned in bits <0:3> of rvec[i]:
		#if 1
			rvec[i] = twopmodq78_q4_GPU(p, pshift, FERMAT, k0,k1,k2,k3, start_index, zshift, i);
		#else
			uint8 r0,r1,r2,r3;
			r0 = twopmodq78_gpu(p, pshift, FERMAT, k0, start_index, zshift, i);
			r1 = twopmodq78_gpu(p, pshift, FERMAT, k1, start_index, zshift, i);
			r2 = twopmodq78_gpu(p, pshift, FERMAT, k2, start_index, zshift, i);
			r3 = twopmodq78_gpu(p, pshift, FERMAT, k3, start_index, zshift, i);
			rvec[i] = r0 + (r1<<1) + (r2<<2) + (r3<<3);
		#endif
		}
	}

/*
	__global__ void GPU_TF78_q8(const uint64 p, const uint64 pshift, const uint32 FERMAT, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		if(i < N) {	// Allow for partially-filled thread blocks
			int i8 = i << 3;	// *= 8 here to yield base-ID for each k-octtet
			uint64 k0 = kvec[i8], k1 = kvec[i8+1], k2 = kvec[i8+2], k3 = kvec[i8+3], k4 = kvec[i8+4], k5 = kvec[i8+5], k6 = kvec[i8+6], k7 = kvec[i8+7];
			// Result returned in bits <0:7> (that is, all 8 bits) of rvec[i]:
			rvec[i] = twopmodq78_q8_GPU(p, pshift, FERMAT, k0,k1,k2,k3,k4,k5,k6,k7, start_index, zshift, i);
		}
	}
*/
#endif	// __CUDACC__

/* 26 Dec 2014: Woohoo! "First light" for GPU-based Mfactor:

	INFO: No factoring savefile t7818977 found ... starting from scratch.
	Generating difference table of first 100000 small primes
	Using first 100000 odd primes; max gap = 114
	max sieving prime = 1299721
	searching in the interval k=[1596054739319040, 1596076891368960], i.e. q=[2.495903e+22, 2.495938e+22]
	each of 16 (p mod 60) passes will consist of 1356 intervals of length 272272
	TRYQ = 4, max sieving prime = 1299721
	Time to set up sieve = 00:00:00.060
	pass = 0...........
	pass = 1...........
	pass = 2...........
	pass = 3...........
	pass = 4...........
	pass = 5...........
	pass = 6...........
	pass = 7...........
	pass = 8...........
	pass = 9...........
	pass = 10...........
	pass = 11...........
	pass = 12...........
	pass = 13...........
	pass = 14...........
	pass = 15..k = 1596059422755180, result = 1
	GPU_TF78: Nonzero result[337938] 1 for p = 7818977, k = 1596059422755180

        rvec[337938] = 1: Factor with k = 1596059422755180. Program: E3.0x
*/

uint64 twopmodq78_3WORD_DOUBLE(uint64 p, uint64 k)
{
	const char func[] = "twopmodq78_3WORD_DOUBLE";
#if FAC_DEBUG
	int dbg = 0;//(p==2147483647 && k==20269004);
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead7, hi64;
	uint96 q, qhalf, qinv, x, lo;
	uint192 prod192;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	double fq0, fq1, fq2;
	double fqinv0, fqinv1, fqinv2;
	double fx0, fx1, fx2;
	double flo0,flo1,flo2,fhi0,fhi1,fhi2;
	/* Experimental: Vars needed for post-powering-loop x * 2^(96-78) modmul: To test, try changing code below to e.g. pshift = p+96i */
	double qdiv, qmul, kmul;
	uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

#if FAC_DEBUG
	if(dbg) {
		printf("%s with p = %llu, k = %llu\n",func,p,k);
	}
#endif
	ASSERT((p >> 63) == 0, "twopmodq78_q2 : p must be < 2^63!");
	q.d0 = p+p;
#ifdef MUL_LOHI64_SUBROUTINE
	// MUL_LOHI64 expects a 64-bit high-part pointer, in 32bit builds this buggers us if we try dumping hi-part directly into 32-bit q.d1
	MUL_LOHI64(q.d0, k,&q.d0,&hi64);	q.d1 = hi64;
#else
	MUL_LOHI64(q.d0, k, q.d0, q.d1);
#endif
	q.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	ASSERT((q.d1 >> 14) == 0, "twopmodq78 : (q.d1 >> 14) != 0");

	/* Convert q to floating form: */
	CVT_UINT78_3WORD_DOUBLE(q, fq0,fq1,fq2);
	RSHIFT_FAST96(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 78;	// Can use 78 or 96 here - cf. adjustment code imm. preceding return

	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 7 or 8 to account for the fact that we can do the powering for the
	!    leftmost 7 or 8 bits (depending on whether the leftmost 7 > 95 or not) via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming
		the leftmost 6/7 bits have already been processed via a shift.

		Since 7 bits with leftmost bit = 1 is guaranteed
		to be in [64,127], the shift count here is in [0, 63].
		That means that zstart < 2^64. Together with the fact that
		squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration.
		Every little bit counts (literally in this case :), right?
		*/
		j = leadz64(pshift);
		/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 77: */
		lead7 = ((pshift<<j) >> 57);
		if(lead7 > 77)
		{
			lead7 >>= 1;
			start_index =  64-j-6;	/* Use only the leftmost 6 bits */
		}
		else
			start_index =  64-j-7;

		zshift = 77 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	  #if FAC_DEBUG
		if(dbg) printf("pshift = 0x%llX\n",pshift);
	  #endif
	}

	/*
	!    Find modular inverse (mod 2^78) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv.d0 = (q.d0 + q.d0 + q.d0) ^ (uint64)2;	qinv.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*(2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
	defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
	*/
	for(j = 0; j < 4; j++)
	{
		hi64 = q.d0*qinv.d0;
		qinv.d0 = qinv.d0*((uint64)2 - hi64);
	}

	/* Now that have bottom 64 bits of qinv, do one more Newton iteration using full 96-bit operands.

	Since qinv.d1 = 0 and q.d0*qinv.d0 == 1 mod 2^64 here; can take advantage of this to speed the iteration:

		q*qinv	= (q.d1*2^64 + q.d0)*qinv.d0
				= q.d1*qinv.d0*2^64 + q.d0*qinv.d0, but this == 1 mod 2^64, i.e.
				= (q.d1*qinv.d0 + MULH64(q.d0,qinv.d0))*2^64 + 1 ,

	i.e. do a MULL64(q.d1,qinv.d0) + MULH64(q.d0,qinv.d0) to get bits 64-127 of q*qinv, set bits 0-63 = 1,
	simply negate bits 64-127 to get r = 2-q*qinv (mod 2^128), then use that (bits 0-63 of r) = 1 to
	simplify the MULL128(qinv,r) operation, i.e.

		qinv*r = qinv.d0*(r.d1*2^64 + 1) == (MULL64(qinv.d0,r.d1)*2^64 + qinv.d0) mod 2^128 .

	In the 96-bit case we then zero the upper 32 bits of the 128-bit result.
	This needs a total of just 3 MULs and 2 ALUs, compared to 8 MULs and 8 ALUs for the naive sequence.
	(May be able to squeeze a bit more out by using 32-bit MULLs rather than 64-bit MULLs in some places,
	but since the main action is in the j-loop below, don't bother with further speedups here.
	*/
	/* qinv has 128 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + __MULH64(q.d0, qinv.d0));
#else
	MULH64(q.d0, qinv.d0, hi64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + hi64);
#endif
	qinv.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */

//	printf("twopmodq78_3WORD_DOUBLE with p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", (uint32)p,(uint32)pshift,k,zshift,start_index);

	/* Convert qinv to floating form: */
/*	cvt_uint78_3word_double(qinv, &fqinv0,&fqinv1,&fqinv2);	*/
	CVT_UINT78_3WORD_DOUBLE(qinv, fqinv0,fqinv1,fqinv2);
#if FAC_DEBUG
	if(dbg) {
		printf("fq = %20.15f, %20.15f, %20.15f\n",fq0,fq1,fq2);
		printf("fqinv = %20.15f, %20.15f, %20.15f\n",fqinv0,fqinv1,fqinv2);
	}
#endif

	/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT96(qinv, zshift, lo);
	lo.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */

	MUL_LOHI96_PROD192(q,lo,prod192);
	RSHIFT192(prod192,78,prod192);
	lo.d0 = prod192.d0;	lo.d1 = prod192.d1;

	/* hi = 0 in this instance, which simplifies things. */
	SUB96(q, lo, x);	/* Put the result in lo (rather than x), to ease overflow check below */

	if((pshift >> j) & (uint64)1)
	{
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT96(x, qhalf)){ ADD96(x, x, x); SUB96(x, q, x); }else{ ADD96(x, x, x); }
	}

	/* Convert x to floating form: */
	CVT_UINT78_3WORD_DOUBLE(x, fx0,fx1,fx2);

	for(j = start_index-2; j >= 0; j--)
	{
	#if FAC_DEBUG
		if(dbg) {
			printf("J = %d, 2x = %d: fx = %20.15f, %20.15f, %20.15f\n",j,((pshift >> j) & (uint64)1) == 1,fx0,fx1,fx2);
		}
	#endif
		/*...x^2 mod q is returned in x. */
		/* SQR_LOHI96(x,lo,hi); */
		SQR_LOHI78_3WORD_DOUBLE(fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2);
		/* MULL96(lo,qinv,lo); */
		MULL78_3WORD_DOUBLE(flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2);
		/* MULH96(q,lo,lo); */
		MULH78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
		{
			SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
			ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
		}
		/* Normalize the result: */
		NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

		if((pshift >> j) & (uint64)1)
		{
			/* Convert fx to uint96 for purpose of this particular comparison: */
			CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x, qhalf))
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
			}
			/* Normalize the result: */
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);
		}
	}
#if FAC_DEBUG
	if(dbg) {
		printf("fx = %20.15f, %20.15f, %20.15f\n",fx0,fx1,fx2);
	}
#endif

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);

if(~pshift != p+78) {
	ASSERT(~pshift > (p+78), "twopmodq80 : Only support pshift >= true value!");
	ASSERT((~pshift - (p+78)) < 32, "twopmodq80 : Only support pshift-diff < 32!");
	qmul  = fx0 + fx1*TWO26FLOAT;
	qmul += fx2*TWO26FLOAT*TWO26FLOAT;
	// Extra power of 2 is because in this flow we do not do the final 2*x-q step in the 'else' below:
	kmul = (double)((uint64)1 << (~pshift - (p+77)));
	qmul *= kmul;
	// This is the multiplicative inverse of q, not to be confused with the Montgomery-inverse:
	qdiv  = fq0 + fq1*TWO26FLOAT;
	qdiv += fq2*TWO26FLOAT*TWO26FLOAT;
	// Compute approx. ratio:
	qmul /= qdiv;
	qmul = DNINT(qmul);
	// Use 96-bit variable to store q*qmul in preparation for exact-arithmetic subtract:
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64( x.d0, (uint64)kmul,&x.d0,&hi64);
#else
	MUL_LOHI64( x.d0, (uint64)kmul, x.d0, hi64);
#endif
	x.d1 = hi64 + kmul*x.d1;
	lo = q;
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(lo.d0, (uint64)qmul,&lo.d0,&hi64);
#else
	MUL_LOHI64(lo.d0, (uint64)qmul, lo.d0, hi64);
#endif
	lo.d1 = hi64 + qmul*lo.d1;
	lo.d0 -= FERMAT;
	SUB96(x,lo,x);
#if FAC_DEBUG
	if(dbg) {
		printf("X_out[A] = %u*2^64 + %llu\n", x.d1,x.d0);
	}
#endif
} else {
	ADD96(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
	q.d0 -= FERMAT;
	SUB96(x,q,x);
#if FAC_DEBUG
	if(dbg) {
		printf("X_out[B] = %u*2^64 + %llu\n", x.d1,x.d0);
	}
#endif
}
	return (CMPEQ96(x, ONE96));
}

/****************************************************************************************/
/********** x86-asm-using functions not buildable on k1om, replace with stubs: **********/
/****************************************************************************************/
#ifdef USE_IMCI512

	uint64 twopmodq78_3WORD_DOUBLE_q2(uint64 p, uint64 k0, uint64 k1, int init_sse2, int thr_id) {
		ASSERT(0,"twopmodq78_3WORD_DOUBLE_q2 cannot be called for k1om builds!");	return 0;
	}
	uint64 twopmodq78_3WORD_DOUBLE_q4(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, int init_sse2, int thr_id) {
		ASSERT(0,"twopmodq78_3WORD_DOUBLE_q4 cannot be called for k1om builds!");	return 0;
	}
	uint64 twopmodq78_3WORD_DOUBLE_q4_REF(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3) {
		ASSERT(0,"twopmodq78_3WORD_DOUBLE_q4_REF cannot be called for k1om builds!");	return 0;
	}
	uint64 twopmodq78_3WORD_DOUBLE_q8(uint64 p, uint64 k[], int init_sse2, int thr_id) {
		ASSERT(0,"twopmodq78_3WORD_DOUBLE_q8 cannot be called for k1om builds!");	return 0;
	}
	uint64 twopmodq78_3WORD_DOUBLE_q16(uint64 p, uint64 k[], int init_sse2, int thr_id) {
		ASSERT(0,"twopmodq78_3WORD_DOUBLE_q16 cannot be called for k1om builds!");	return 0;
	}
	uint64 twopmodq78_3WORD_DOUBLE_q32(uint64 p, uint64 k[], int init_sse2, int thr_id) {
		ASSERT(0,"twopmodq78_3WORD_DOUBLE_q32 cannot be called for k1om builds!");	return 0;
	}
	uint64 twopmodq78_3WORD_DOUBLE_q64(uint64 p, uint64 k[], int init_sse2, int thr_id) {
		ASSERT(0,"twopmodq78_3WORD_DOUBLE_q64 cannot be called for k1om builds!");	return 0;
	}

#else

	/*
	This routine uses an SSE2 version of the FLOOR function, which normally needs the rounding mode
	to be set to round-toward-minus-infinity, so we need to twiddle the MXCSR register -
	from http://softpixel.com/~cwright/programming/simd/sse.php :

	The MXCSR register is a 32-bit register containing flags for control and status information regarding SSE instructions.
	As of SSE3, only bits 0-15 have been defined.

	Mnemonic	Bit Location	Description				[EWM: Default value on MSVC/ia32 = 0x1FA0, so the bits marked with [x] are set:]
	--------	------------	---------------------	-----
		FZ		bit 15			Flush To Zero
		R+		bit<14:13> = 10	Round Positive
		R-		bit<14:13> = 01	Round Negative
		RZ		bit<14:13> = 11	Round To Zero
		RN		bit<14:13> = 00	Round To Nearest		[x]
		PM		bit 12			Precision Mask			[x]
		UM		bit 11			Underflow Mask			[x]
		OM		bit 10			Overflow Mask			[x]
		ZM		bit 9			Divide By Zero Mask		[x]
		DM		bit 8			Denormal Mask			[x]
		IM		bit 7			Invalid Operation Mask	[x]
		DAZ		bit 6			Denormals Are Zero
		PE		bit 5			Precision Flag			[x]
		UE		bit 4			Underflow Flag
		OE		bit 3			Overflow Flag
		ZE		bit 2			Divide By Zero Flag
		DE		bit 1			Denormal Flag
		IE		bit 0			Invalid Operation Flag

	So to enable round-toward-minus-infinity, we set bit 13.

	The problem with this is that we need to do it and undo it every pass through the squaring sequence:

		__asm	stmxcsr	mxcsr_ptr
		__asm	xor		mxcsr_ptr, 0x2000
		__asm	ldmxcsr	mxcsr_ptr
			__asm	movaps	xmm6,xmm2	// fcy = cpy of fprod2
			__asm cvtpd2dq	xmm6,xmm6	// Convert to 32-bit int, rounding toward -oo
			__asm cvtdq2pd	xmm6,xmm6	// ...and back to double.
		__asm	xor		mxcsr_ptr, 0x2000
		__asm	ldmxcsr	mxcsr_ptr

	So we instead emulate [round toward -oo](x) via DNINT(x-(0.5-epsilon)) and use the defsault round-to-nearest mode.
	With a fast floating-point emulation of the DNINT() function (our "DNINT" macro), this will generally run much faster
	than alternatives such as the above, or ones involving cast-to-integer and back.
	*/

	/*** 2-trial-factor version ***/
	uint64 twopmodq78_3WORD_DOUBLE_q2(uint64 p, uint64 k0, uint64 k1, int init_sse2, int thr_id)
	{
		const char func[] = "twopmodq78_3WORD_DOUBLE_q2";
	#if FAC_DEBUG
		int dbg = 0;//(p==2147483647 && (k0==20269004 || k1==20269004));
	#endif
		 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
		uint64 lead7, tmp0, tmp1, r;
		uint96 q0, q1, qinv0, qinv1, qhalf0, qhalf1, x0, x1, lo0, lo1;
		uint192 prod192;
		static uint64 psave = 0, pshift;
		static uint32 start_index, zshift, first_entry = TRUE;
		uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		const uint64 ihalf = 0x3FDfffffffffffffull;	/* Bitfield storing 0.5*(1-epsilon) in IEEE64 format */
		static int max_threads = 1;	// Default local-array-init is for just a single thread ... caller can re-init for > 1 threads later, if desired.
		static double *sc_arr = 0x0, *sc_ptr;
		vec_dbl *tmp;
		double dtmp;
	  #ifdef MULTITHREAD
		static double *__r0;	// Base address for discrete per-thread local stores ...
								// *NOTE*: more convenient to use double* rather than vec_dbl* here
	  #else
		static	// Following set of pointers only static-izable in single-thread mode
	  #endif
		double *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2
					, *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2
					, *half, *two26f, *two26i, *two13i, *sse2_rnd;
	#else

		double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2, fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2;
		double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2;
		// Note: In ||-init mode, *value* of init_sse2 to store #threads-to-init-for:
		if(init_sse2) {
			ASSERT(init_sse2 <= 1, "Multithreading currently only supported for SIMD builds!");
			return 0;	// In non-SIMD mode, ||-init call is a no-op
		}

	#endif

	#if FAC_DEBUG
		if(dbg) {
			printf("%s with p = %llu, k0 = %llu, k1 = %llu\n",func,p,k0,k1);
		}
	#endif

		if(p != psave)
		{
		//	first_entry = FALSE;
			psave  = p;
			pshift = p + 78;
		/*
		!    find number of leading zeros in p, use it to find the position of the leftmost
		!    ones bit, and subtract 7 or 8 to account for the fact that we can do the powering for the
		!    leftmost 7 or 8 bits (depending on whether the leftmost 7 > 95 or not) via a simple shift.
		*/
			/* Leftward bit at which to start the l-r binary powering, assuming the leftmost 6/7 bits have already been processed via a shift.

			Since 7 bits with leftmost bit = 1 is guaranteed to be in [64,127], the shift count here is in [0, 63].
			That means that zstart < 2^64. Together with the fact that squaring a power of two gives another power of two, we can
			simplify the modmul code sequence for the first iteration. Every little bit counts (literally in this case :), right?
			*/
			j = leadz64(pshift);
			/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 77: */
			lead7 = ((pshift<<j) >> 57);
			if(lead7 > 77)
			{
				lead7 >>= 1;
				start_index =  64-j-6;	/* Use only the leftmost 6 bits */
			}
			else
				start_index =  64-j-7;

			zshift = 77 - lead7;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

			pshift = ~pshift;
		  #if FAC_DEBUG
			if(dbg) printf("pshift = 0x%llX\n",pshift);
		  #endif
		}

	#ifdef USE_SSE2

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
			// Alloc the local-memory block:
			sc_arr = ALLOC_DOUBLE(sc_arr, 0x2c*max_threads + 4);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);	// Force vec_dbl-alignment
			ASSERT(((uintptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		#ifdef MULTITHREAD
			__r0  = sc_ptr;
			two13i  = sc_ptr + 0x18;
			two26f  = sc_ptr + 0x1a;
			two26i  = sc_ptr + 0x1c;
			sse2_rnd= sc_ptr + 0x1e;
			half    = sc_ptr + 0x20;
			dtmp = *(double*)&ihalf;
			for(j = 0; j < max_threads; ++j) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT_2((vec_dbl*)two13i  , TWO13FLINV);
				VEC_DBL_INIT_2((vec_dbl*)two26f  , TWO26FLOAT);
				VEC_DBL_INIT_2((vec_dbl*)two26i  , TWO26FLINV);
				VEC_DBL_INIT_2((vec_dbl*)sse2_rnd, crnd      );
				VEC_DBL_INIT_2((vec_dbl*)half    , dtmp      );
				// Move on to next thread's local store:
				two13i   += 0x2c;
				two26f   += 0x2c;
				two26i   += 0x2c;
				sse2_rnd += 0x2c;
				half     += 0x2c;
			}
		#else
			/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 2 to span an SSE register: */
			fq0    = sc_ptr + 0x00;		gq0    = sc_ptr + 0x01;
			fq1    = sc_ptr + 0x02;		gq1    = sc_ptr + 0x03;
			fq2    = sc_ptr + 0x04;		gq2    = sc_ptr + 0x05;
			fqhi52 = sc_ptr + 0x06;		gqhi52 = sc_ptr + 0x07;
			fqinv0 = sc_ptr + 0x08;		gqinv0 = sc_ptr + 0x09;
			fqinv1 = sc_ptr + 0x0a;		gqinv1 = sc_ptr + 0x0b;
			fqinv2 = sc_ptr + 0x0c;		gqinv2 = sc_ptr + 0x0d;
			fx0    = sc_ptr + 0x0e;		gx0    = sc_ptr + 0x0f;
			fx1    = sc_ptr + 0x10;		gx1    = sc_ptr + 0x11;
			fx2    = sc_ptr + 0x12;		gx2    = sc_ptr + 0x13;
			// +0x14,15,16,17 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
			two13i = sc_ptr + 0x18;
			two26f = sc_ptr + 0x1a;
			two26i = sc_ptr + 0x1c;
			sse2_rnd=sc_ptr + 0x1e;
			half   = sc_ptr + 0x20;
			/* Can premultiply each of the multiword-mul inputs by 1/sqrt(2^13) due to quadraticity: */
			*two13i++ = TWO13FLINV;		*two13i-- = TWO13FLINV;
			*two26f++ = TWO26FLOAT;		*two26f-- = TWO26FLOAT;
			*two26i++ = TWO26FLINV;		*two26i-- = TWO26FLINV;
			/* SSE2 math = 53-mantissa-bit IEEE double-float: */
			*sse2_rnd++ = crnd;		*sse2_rnd-- = crnd;
			/* We init "half" = 0.5-epsilon here, because emulating FLOOR(x) via DNINT(x-half) requires
			us to always round up if x is a whole number, and our DNINT emulation can round either way if fractional part = 0.5:
			*/
			dtmp = *(double*)&ihalf;
			*half++ = dtmp;			*half-- = dtmp;
		#endif
			if(init_sse2) return 0;
		}	/* end of inits */

		/* If multithreaded, set the local-store pointers needed for the current thread; */
	  #ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		sc_ptr = __r0 + thr_id*0x2c;
		/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 2 to span an SSE register: */
		fq0    = sc_ptr + 0x00;		gq0    = sc_ptr + 0x01;
		fq1    = sc_ptr + 0x02;		gq1    = sc_ptr + 0x03;
		fq2    = sc_ptr + 0x04;		gq2    = sc_ptr + 0x05;
		fqhi52 = sc_ptr + 0x06;		gqhi52 = sc_ptr + 0x07;
		fqinv0 = sc_ptr + 0x08;		gqinv0 = sc_ptr + 0x09;
		fqinv1 = sc_ptr + 0x0a;		gqinv1 = sc_ptr + 0x0b;
		fqinv2 = sc_ptr + 0x0c;		gqinv2 = sc_ptr + 0x0d;
		fx0    = sc_ptr + 0x0e;		gx0    = sc_ptr + 0x0f;
		fx1    = sc_ptr + 0x10;		gx1    = sc_ptr + 0x11;
		fx2    = sc_ptr + 0x12;		gx2    = sc_ptr + 0x13;
		// +0x14,15,16,17 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
		two13i = sc_ptr + 0x18;
		two26f = sc_ptr + 0x1a;
		two26i = sc_ptr + 0x1c;
		sse2_rnd=sc_ptr + 0x1e;
		half   = sc_ptr + 0x20;
	//	printf("Thr %d ONE96_PTR address = %llX; data.d0,d1 = %llu,%u\n",thr_id,(uint64)ONE96_PTR,ONE96_PTR->d0,ONE96_PTR->d1);
		tmp = (vec_dbl*)sse2_rnd; ASSERT((tmp->d0 == crnd) && (tmp->d1 == crnd), "Bad data at sse2_rnd address!");
	  #endif

	#endif

		ASSERT((p >> 63) == 0, "twopmodq78_q2 : p must be < 2^63!");
		q0.d0 = q1.d0 = p+p;
	#ifdef MUL_LOHI64_SUBROUTINE
		// MUL_LOHI64 expects a 64-bit high-part pointer, in 32bit builds this buggers us if we try dumping hi-part directly into 32-bit q.d1
		MUL_LOHI64(q0.d0, k0,&q0.d0,&tmp0);	q0.d1 = tmp0;
		MUL_LOHI64(q1.d0, k1,&q1.d0,&tmp0);	q1.d1 = tmp0;
	#else
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
	#endif
		q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
		q1.d0 += 1;
		ASSERT((q0.d1 >> 14) == 0, "twopmodq78_q2 : (q0.d1 >> 14) != 0");
		ASSERT((q1.d1 >> 14) == 0, "twopmodq78_q2 : (q1.d1 >> 14) != 0");

		/* Convert q to floating form: */
	#ifdef USE_SSE2

		CVT_UINT78_3WORD_DOUBLE(q0 ,*fq0,*fq1,*fq2);
		CVT_UINT78_3WORD_DOUBLE(q1 ,*gq0,*gq1,*gq2);

	  #ifdef USE_ARM_V8_SIMD
		// No TF support on ARMv8, just supply a stub macro:
		__asm__ volatile (\
			"ldr x0,%[__fq0]	\n\t"\
			:					// outputs: none
			: [__fq0] "m" (fq0)	// All inputs from memory addresses here
			: "cc","memory","x0"	/* Clobbered registers */\
		);

	  #elif OS_BITS == 32

		#error 32-bit OSes no longer supported for SIMD builds!

	  #else

		__asm__ volatile (\
			"movq	%[__fq0],%%rax	\n\t"\
			"movq	%[__two26f],%%rsi \n\t"\
			"movaps	    (%%rsi),%%xmm5	\n\t"\
			"movaps	0x10(%%rsi),%%xmm6	\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t"\
			"movaps	0x10(%%rax),%%xmm1	\n\t"\
			"movaps	0x20(%%rax),%%xmm2	\n\t"\
			"movaps	%%xmm2,%%xmm3		\n\t"\
			"mulpd	%%xmm5,%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm3		\n\t"\
			"mulpd	%%xmm6,%%xmm0		\n\t"\
			"mulpd	%%xmm6,%%xmm1		\n\t"\
			"mulpd	%%xmm6,%%xmm2		\n\t"\
			"movaps	%%xmm0,    (%%rax)	\n\t"\
			"movaps	%%xmm1,0x10(%%rax)	\n\t"\
			"movaps	%%xmm2,0x20(%%rax)	\n\t"\
			"movaps	%%xmm3,0x30(%%rax)	\n\t"\
			:					/* outputs: none */\
			: [__fq0] "m" (fq0)	/* All inputs from memory addresses here */\
			 ,[__two26f] "m" (two26f)\
			: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm5","xmm6"		/* Clobbered registers */\
		);

	  #endif

	  #if FAC_DEBUG
		if(dbg) {
			printf("fq = %20.15f, %20.15f, %20.15f\n", *fq0 * TWO26FLOAT, *fq1 * TWO26FLOAT, *fq2 * TWO26FLOAT);
		}
	  #endif

	#else
		CVT_UINT78_3WORD_DOUBLE(q0, fq0,fq1,fq2);
		CVT_UINT78_3WORD_DOUBLE(q1, gq0,gq1,gq2);
	  #if FAC_DEBUG
		if(dbg) {
			printf("fq = %20.15f, %20.15f, %20.15f\n",fq0,fq1,fq2);
		}
	  #endif
	#endif

		RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
		RSHIFT_FAST96(q1, 1, qhalf1);

		/*
		!    Find modular inverse (mod 2^78) of q in preparation for modular multiply.
		*/
		/* q must be odd for Montgomery-style modmul to work: */
		/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
		qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
		qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;

		/* Newton iteration involves repeated steps of form

			qinv = qinv*(2 - q*qinv);

		Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
		defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
		*/
		for(j = 0; j < 4; j++)
		{
			tmp0 = q0.d0*qinv0.d0;
			tmp1 = q1.d0*qinv1.d0;

			qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
			qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
		}

		/* Now that have bottom 64 bits of qinv, do one more Newton iteration using full 96-bit operands.
		qinv has 128 bits, but only the upper 64 get modified here. */
	#ifdef MUL_LOHI64_SUBROUTINE
		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
	#else
		MULH64(q0.d0, qinv0.d0, tmp0);
		MULH64(q1.d0, qinv1.d0, tmp1);

		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);
	#endif
		qinv0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
		qinv1.d1 &= 0x0000000000003fff;

		/* Convert qinv to floating form: */
	#ifdef USE_SSE2
		CVT_UINT78_3WORD_DOUBLE(qinv0 ,*fqinv0,*fqinv1,*fqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv1 ,*gqinv0,*gqinv1,*gqinv2);

	  #ifdef USE_ARM_V8_SIMD
		// No TF support on ARMv8, just supply a stub macro:
		__asm__ volatile (\
			"ldr x0,%[__fqinv0]	\n\t"\
			:					// outputs: none
			: [__fqinv0] "m" (fqinv0)	// All inputs from memory addresses here
			: "cc","memory","x0"	/* Clobbered registers */\
		);

	  #elif OS_BITS == 32

		#error 32-bit OSes no longer supported for SIMD builds!

	  #else

		__asm__ volatile (\
			"movq	%[__fqinv0],%%rax	\n\t"\
			"movq	%[__two26i],%%rsi	\n\t"\
			"movaps	    (%%rsi),%%xmm6	\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t"\
			"movaps	0x10(%%rax),%%xmm1	\n\t"\
			"movaps	0x20(%%rax),%%xmm2	\n\t"\
			"mulpd	%%xmm6,%%xmm0		\n\t"\
			"mulpd	%%xmm6,%%xmm1		\n\t"\
			"mulpd	%%xmm6,%%xmm2		\n\t"\
			"movaps	%%xmm0,    (%%rax)	\n\t"\
			"movaps	%%xmm1,0x10(%%rax)	\n\t"\
			"movaps	%%xmm2,0x20(%%rax)	\n\t"\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm6"	/* Clobbered registers */\
		);

	  #endif

	  #if FAC_DEBUG
		if(dbg) {
			printf("fqinv = %20.15f, %20.15f, %20.15f\n", *fqinv0 * TWO26FLOAT, *fqinv1 * TWO26FLOAT, *fqinv2 * TWO26FLOAT);
		}
	  #endif

	#else

		CVT_UINT78_3WORD_DOUBLE(qinv0, fqinv0,fqinv1,fqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv1, gqinv0,gqinv1,gqinv2);
	  #if FAC_DEBUG
		if(dbg) {
			printf("fqinv = %20.15f, %20.15f, %20.15f\n",fqinv0,fqinv1,fqinv2);
		}
	  #endif

	#endif

		/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
		j = start_index-1;

		/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
		LSHIFT96(qinv0, zshift, lo0);	lo0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
		LSHIFT96(qinv1, zshift, lo1);	lo1.d1 &= 0x0000000000003fff;

		MUL_LOHI96_PROD192(q0,lo0,prod192);	RSHIFT192(prod192,78,prod192);	lo0.d0 = prod192.d0;	lo0.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q1,lo1,prod192);	RSHIFT192(prod192,78,prod192);	lo1.d0 = prod192.d0;	lo1.d1 = prod192.d1;

		/* hi = 0 in this instance, which simplifies things. */
		SUB96(q0, lo0, x0);	/* Put the result in lo (rather than x), to ease overflow check below */
		SUB96(q1, lo1, x1);

		if((pshift >> j) & (uint32)1)
		{
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
			if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
		}

		/* Convert x to floating form: */
	#ifdef USE_SSE2
		CVT_UINT78_3WORD_DOUBLE(x0 ,*fx0,*fx1,*fx2);
		CVT_UINT78_3WORD_DOUBLE(x1 ,*gx0,*gx1,*gx2);
	#else
		CVT_UINT78_3WORD_DOUBLE(x0, fx0,fx1,fx2);
		CVT_UINT78_3WORD_DOUBLE(x1, gx0,gx1,gx2);
	#endif

		/*...x^2 mod q is returned in x. */
		/* All 3-word-double-form operands have components in the following size ranges:
			fword0,1 in [-2^25, +2^25]
			fword2   in [   -1, +2^26]
		*/
	#ifndef USE_SSE2

		for(j = start_index-2; j >= 0; j--)
		{
		#if FAC_DEBUG
			if(dbg) {
				printf("J = %d: fx = %20.15f, %20.15f, %20.15f\n",j,fx0,fx1,fx2);
			}
		#endif
			/*...x^2 mod q is returned in x. */
			/* SQR_LOHI96(x,lo,hi); */
			SQR_LOHI78_3WORD_DOUBLE_q2(
				  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2
				, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2
			);
			/* MULL96(lo,qinv,lo); */
			MULL78_3WORD_DOUBLE_q2(
				  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
				, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
			);
			/* MULH96(q,lo,lo); */
			MULH78_3WORD_DOUBLE_q2(
				  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2
				, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2
			);

			/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
			if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
			{
				SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
				ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
			}
			/* Normalize the result: */
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

			if(CMPLT78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2))
			{
				SUB78_3WORD_DOUBLE(gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2);
				ADD78_3WORD_DOUBLE(glo0,glo1,glo2, ghi0,ghi1,ghi2, gx0,gx1,gx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2, gx0,gx1,gx2);
			}
			/* Normalize the result: */
			NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

			if((pshift >> j) & (uint64)1)
			{
				/* Convert fx to uint96 for purpose of this particular comparison: */
				CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
				CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x0, qhalf0))
				{
					ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
					SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				}
				/* Normalize the result: */
				NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x1, qhalf1))
				{
					ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
					SUB78_3WORD_DOUBLE(gx0,gx1,gx2, gq0,gq1,gq2, gx0,gx1,gx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
				}
				/* Normalize the result: */
				NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);
			}
		}	/* for(j...) */

	#else	// USE_SSE2:

		/* Inner loop body needs 42 movaps, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */

	  #ifdef USE_ARM_V8_SIMD
		// No TF support on ARMv8, just supply a stub macro:
		__asm__ volatile (\
			"ldr x0,%[__fq0]	\n\t"\
			:					// outputs: none
			: [__fq0] "m" (fq0)	// All inputs from memory addresses here
			: "cc","memory","x0"	/* Clobbered registers */\
		);

	  #elif OS_BITS == 32

		#error 32-bit OSes no longer supported for SIMD builds!
		ASSERT((uint32)(~pshift) == 0, "p+78 must be 32-bit here for 32-bit ASM support!");

	  #else	// The 64-bit version of the macro is timing-suboptimal because I used it as a testbed:
				// This 2-TF-input/4-xxm-register version serves as the basis for an 8-input version
				// using all 16 xmm-regs available in 64-bit mode. Doing 4 independent 4-register-using
				// instruction streams in parallel should allo for better latency hiding than we get with
				// original 2-input/8-register implementation and its translation to 4-input/16-register,
				// which can be seen in the _q4 version of this routine.
		for(j = start_index-2; j >= 0; j--)
		{
		#if FAC_DEBUG
			if(dbg) {
				printf("J = %d, 2x = %d: fx = %20.15f, %20.15f, %20.15f\n",j,((pshift >> j) & (uint64)1) == 1,*fx0,*fx1,*fx2);
			}
		#endif
		__asm__ volatile (\
		/* SQR_LOHI78_3WORD_DOUBLE_q2(fx, flo,fhi). Inputs flo0,1,2 enter in xmm0,1,2: */\
			"movq	%[__fx0],%%rax		\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t"\
			"movaps	0x10(%%rax),%%xmm1	\n\t"\
			"movaps	0x20(%%rax),%%xmm2	\n\t"\
			"movq	%[__two26i],%%rbx	\n\t"\
			"movaps	-0x20(%%rbx),%%xmm3	\n\t"/* two13i */\
			"movaps	 0x10(%%rbx),%%xmm7	\n\t"/* xmm7 = RND ... only need to burn 5th reg here in SSE2 mode, since ROUNDPD only became available in SSE4.2 */\
			"/* fx0,1,2 assumed in xmm0,1,2 on loop entry */\n\t"\
			"mulpd	%%xmm3,%%xmm0		\n\t"/* scale [fx0 . 2^-13] */\
			"mulpd	%%xmm3,%%xmm1		\n\t"/* scale [fx1 . 2^-13] */\
			"mulpd	%%xmm3,%%xmm2		\n\t"/* scale [fx2 . 2^-13] */\
			"movaps	%%xmm0,%%xmm3		\n\t"/* cpy of fx0 */\
			"addpd	%%xmm3,%%xmm3		\n\t"/* 2.fx0 */\
			"movaps	%%xmm2,0x20(%%rax)	\n\t"/* Store fx2 to free up a register */\
			"/* Digit 0: */\n\t"\
			"mulpd	%%xmm0,%%xmm0		\n\t"/* [  fx0 *= fx0] / 2^26 */\
			"movaps	%%xmm0,%%xmm2		\n\t"/* Init: FCY = dnint(fx0*fx0) */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t"/* fprod0 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fprod0 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* FCY    *= two26i */\
			"movaps	%%xmm0,    (%%rax)	\n\t"/* Store fprod0 to free up a register */\
			"/* Digit 1: */\n\t"\
			"movaps	%%xmm1,%%xmm0		\n\t"/* cpy of fx1 */\
			"mulpd	%%xmm3,%%xmm0		\n\t"/* [fx1 *= 2.fx0] / 2^26 */\
			"addpd	%%xmm2,%%xmm0		\n\t"/* fprod1 = 2.fx0*fx1 + FCY */\
			"movaps	%%xmm0,%%xmm2		\n\t"/* FCY */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t"/* fprod1 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fprod1 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* FCY    *= two26i */\
			"movaps	%%xmm0,0x10(%%rax)	\n\t"/* Store fprod1 to free up a register */\
			"/* Digit 2: */\n\t"\
			"movaps	%%xmm1,%%xmm0		\n\t"/* cpy of fx1 */\
			"addpd	%%xmm0,%%xmm0		\n\t"/* 2.fx1 */\
			"mulpd	%%xmm1,%%xmm1		\n\t"/* [  fx1 *= fx1] / 2^26 */\
			"mulpd	0x20(%%rax),%%xmm3		\n\t"/* [2.fx0 *= fx2] / 2^26 */\
			"addpd	%%xmm2,%%xmm1		\n\t"/* fx1*fx1 += FCY */\
			"addpd	%%xmm3,%%xmm1		\n\t"/* fprod2 = 2.fx0*fx2 + fx1*fx1 + FCY */\
			"movaps	%%xmm1,%%xmm2		\n\t"/* fcy = cpy of fprod2 */\
			"subpd	0x20(%%rbx),%%xmm2	\n\t"/* fprod2 - 0.5 */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"/* fcy = floor(fprod2) */\
			"subpd	%%xmm2,%%xmm1		\n\t"/* fprod2 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm1	\n\t"/* fprod2 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* FCY    *= two26i */\
			"/* Digit 3: */\n\t"\
			"movaps	0x20(%%rax),%%xmm3	\n\t"/* Reload fx2 */\
			"mulpd	%%xmm3,%%xmm0		\n\t"/* [2.fx1 *= fx2] / 2^26 */\
			"addpd	%%xmm2,%%xmm0		\n\t"/* fprod3 = 2.fx1*fx2 + FCY */\
			"movaps	%%xmm0,%%xmm2		\n\t"/* FCY */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t"/* fprod3 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fprod3 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* FCY    *= two26i */\
			"/* Digit 4: */\n\t"\
			"mulpd	%%xmm3,%%xmm3		\n\t"/* [  fx2 *= fx2] / 2^26 */\
			"addpd	%%xmm2,%%xmm3		\n\t"/* fprod4 = fx2*fx2 + fcy */\
			"movaps	%%xmm3,%%xmm2		\n\t"/* FCY */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm2,%%xmm3		\n\t"/* fprod4 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm3	\n\t"/* fprod4 *= two26f */\
			"/* Digit 5 = the carry. CY in xmm2; fhi0,1 in xmm0,xmm3 */\n\t"\
			"mulpd	-0x10(%%rbx),%%xmm2	\n\t"/* fhi2 * two26f */\
			"addpd	%%xmm2,%%xmm3		\n\t"/* fhi, top 52 bits; xmm2 FREE */\
			"movaps	%%xmm0,0x30(%%rax)	\n\t"/* Store fhi0,1 = fprod3,4 to free up 2 more registers */\
			"movaps	%%xmm3,0x40(%%rax)	\n\t"/* Recall: fx0,1 in eax+0x[0,2]0, fx2 in xmm1. */\
		/* MULL78 section below needs flo0,1,2 in xmm0,1,2: */\
		/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */\
			"movq	%[__fqinv0],%%rdx	\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t"/* Reload flo0 */\
			"movaps	%%xmm1,%%xmm2		\n\t"/* Move flo2 into xmm2 */\
			"movaps	%%xmm0,%%xmm3		\n\t"/* xmm3 = cpy of x0 */\
			"mulpd	    (%%rdx),%%xmm0	\n\t"/* fprod0 = x0*y0 */\
			"movaps	%%xmm0,%%xmm1		\n\t"/* FCY */\
			"addpd	%%xmm7,%%xmm1		\n\t"\
			"subpd	%%xmm7,%%xmm1		\n\t"\
			"subpd	%%xmm1,%%xmm0		\n\t"/* fprod0 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fprod0 *= two26f */\
			"mulpd	     (%%rbx),%%xmm1	\n\t"/* fcy    *= two26i */\
			"movaps	%%xmm0,    (%%rax)	\n\t"/* Store fprod0 to free up a register */\
			"movaps	%%xmm3,%%xmm0		\n\t"/* xmm0 = cpy of x0 */\
			"/* Digit 1: */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm3	\n\t"/* x0*y1 */\
			"addpd	%%xmm1,%%xmm3		\n\t"/* x0*y1 + fcy; xmm1 FREE */\
			"movaps	0x10(%%rax),%%xmm1	\n\t"/* Reload flo1 */\
			"mulpd	    (%%rdx),%%xmm1	\n\t"/* x1*y0 */\
			"addpd	%%xmm3,%%xmm1		\n\t"/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */\
			"movaps	%%xmm1,%%xmm3		\n\t"/* FCY */\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm3,%%xmm1		\n\t"/* fprod1 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm1	\n\t"/* fprod0 *= two26f */\
			"mulpd	     (%%rbx),%%xmm3	\n\t"/* fcy    *= two26i */\
			"movaps	%%xmm1,0x20(%%rax)	\n\t"/* Store fprod1 [in unused flo2 slot] to free up a register */\
			"movaps	0x10(%%rax),%%xmm1	\n\t"/* Reload flo1 */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"mulpd	    (%%rdx),%%xmm2	\n\t"/* x2*y0 */\
			"mulpd	0x10(%%rdx),%%xmm1	\n\t"/* x1*y1 */\
			"mulpd	0x20(%%rdx),%%xmm0	\n\t"/* x0*y2 */\
			"addpd	%%xmm1,%%xmm2		\n\t"/* x1*y1 + x2*y0; xmm1 FREE */\
			"addpd	%%xmm0,%%xmm3		\n\t"/* x0*y2 + fcy; xmm0 FREE */\
			"addpd	%%xmm3,%%xmm2		\n\t"/* fprod2; xmm3 FREE */\
			"movaps	%%xmm2,%%xmm3		\n\t"/* fcy = cpy of fprod2 */\
			"subpd	0x20(%%rbx),%%xmm3	\n\t"/* fprod2 - 0.5 */\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm7,%%xmm3		\n\t"/* fcy = floor(fprod2) */\
			"subpd	%%xmm3,%%xmm2		\n\t"/* fprod2 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm2	\n\t"/* fprod2 *= two26f */\
			"movaps	    (%%rax),%%xmm0	\n\t"/* Reload fprod0 */\
			"movaps	0x20(%%rax),%%xmm1	\n\t"/* Reload fprod1 [in flo2 slot] */\
		/* MULH96(q,lo,lo) --> lo =(q*lo)/2^78: flo0,1,2 enter in xmm0,1,2. In the comments x = q, y = lo: */\
			"movq	%[__fq0],%%rdx		\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm3		\n\t"/* xmm3 = cpy of y0 */\
			"mulpd	    (%%rdx),%%xmm0	\n\t"/* fprod0 = y0*x0 */\
			"mulpd	    (%%rbx),%%xmm0	\n\t"/* CY    *= two26i */\
			"/* Digit 1: */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm3	\n\t"/* y0 *= x1 */\
			"addpd	%%xmm0,%%xmm3		\n\t"/* y0*x1 + CY; xmm0 FREE */\
			"mulpd	    (%%rdx),%%xmm1	\n\t"/* y1 *= x0 */\
			"addpd	%%xmm3,%%xmm1		\n\t"/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */\
			"mulpd	    (%%rbx),%%xmm1	\n\t"/* CY    *= two26i */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"movaps	%%xmm2,%%xmm3		\n\t"/* xmm3 = cpy of y2 */\
			"movaps	0x20(%%rax),%%xmm0	\n\t"/* Reload y1 [in flo2 slot] */\
			"mulpd	    (%%rdx),%%xmm2	\n\t"/* y2 *= x0 */\
			"mulpd	0x10(%%rdx),%%xmm0	\n\t"/* y1 *= x1 */\
			"addpd	%%xmm1,%%xmm2		\n\t"/* y2*x0 + CY; xmm1 FREE */\
			"movaps	    (%%rax),%%xmm1	\n\t"/* Reload y0 */\
			"mulpd	0x20(%%rdx),%%xmm1	\n\t"/* y0 *= x2 */\
			"addpd	%%xmm0,%%xmm1		\n\t"/* y1*x1 + y0*x2; xmm0 FREE */\
			"addpd	%%xmm1,%%xmm2		\n\t"/* fprod2; xmm1 FREE */\
			"subpd	0x20(%%rbx),%%xmm2	\n\t"/* fprod2 - 0.5 */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"/* fcy = floor(fprod2) */\
			"mulpd	    (%%rbx),%%xmm2	\n\t"/* CY    *= two26i */\
		/* Precompute partial products needed for upper half: */\
			"movaps	%%xmm3,%%xmm1		\n\t"/* xmm1 = cpy of y2 */\
			"movaps	0x20(%%rax),%%xmm0	\n\t"/* Reload y1 [in flo2 slot] */\
			"mulpd	0x20(%%rdx),%%xmm0	\n\t"/* y1 *= x2 */\
			"mulpd	0x10(%%rdx),%%xmm1	\n\t"/* y2 *= x1 */\
			"mulpd	0x20(%%rdx),%%xmm3	\n\t"/* y2 *= x2 */\
			"/* Digit 3: */\n\t"\
			"addpd	%%xmm2,%%xmm1		\n\t"/* y2*x1 += CY */\
			"addpd	%%xmm1,%%xmm0		\n\t"/* fprod3 = y1*x2 + y2*x1 + CY */\
			"movaps	%%xmm0,%%xmm2		\n\t"/* FCY */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t"/* fprod3 -= CY */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fprod3 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* CY    *= two26i */\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm2,%%xmm3		\n\t"/* CY = y2*x2 + CY */\
			"movaps	%%xmm3,%%xmm1		\n\t"/* fprod4 = cpy of CY */\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm3,%%xmm1		\n\t"/* fprod4 -= CY. Note we use the carry in the next section! */\
			"mulpd	-0x10(%%rbx),%%xmm1	\n\t"/* fprod4 *= two26f */\
		/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\
			"movaps	0x40(%%rax),%%xmm2	\n\t"/* fhi, top 52 bits */\
			"mulpd	-0x10(%%rbx),%%xmm3		\n\t"/* flo2 *= two26f */\
			"addpd	%%xmm1,%%xmm3		\n\t"/* flo, top 52 bits */\
			"movaps	0x30(%%rax),%%xmm1	\n\t"/* fhi, low 26 bits */\
			"cmppd	$0x1,%%xmm3,%%xmm2	\n\t"/* bitmask = (fhi < flo)? */\
			"movaps	%%xmm2,(%%rax)	\n\t"/* Store bitmask to free up a register */\
			"movaps	0x40(%%rax),%%xmm2	\n\t"/* fhi, top 52 bits */\
			"subpd	%%xmm0,%%xmm1		\n\t"/* (fhi - flo), low 26 bits, xmm0 FREE */\
			"subpd	%%xmm3,%%xmm2		\n\t"/* (fhi - flo), top 52 bits, xmm3 FREE */\
			"movaps	-0x10(%%rbx),%%xmm3	\n\t"/* 2^26 */\
			"mulpd	    (%%rdx),%%xmm3	\n\t"/* fq, low 26 bits */\
			"movaps	0x30(%%rdx),%%xmm0	\n\t"/* qhi52 :+ fq, top 52 bits */\
			"andpd	(%%rax),%%xmm0		\n\t"/* qhi52 & bitmask */\
			"addpd	%%xmm0,%%xmm2		\n\t"/* xhi = (h-l)hi + (qhi52 & bitmask) */\
			"movaps	%%xmm3,0x10(%%rax)	\n\t"/* Store qlo26 to free up a register */\
			"andpd	(%%rax),%%xmm3		\n\t"/* qlo26 &= bitmask */\
			"addpd	%%xmm3,%%xmm1		\n\t"/* xlo = (h-l)lo + (qlo26 & bitmask) */\
			"/* xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in 0x30(%%rdx),0x10(%%rax): */\n\t"\
		/* if((pshift >> j) & (uint64)1) { */\
			"movq	%[__pshift],%%rax	\n\t"\
			"movslq	%[__j],%%rcx		\n\t"\
			"shrq	%%cl,%%rax		\n\t"\
			"andq	$0x1,%%rax		\n\t"\
		"je twopmodq78_3wdq2_gcc64	\n\t"/* Sep 2016: with FAC_DEBUG enabled, gcc/clang both gave "error: invalid symbol redefinition" errors w.r.to this label - had to turn off optimization. */\
			"movq	%[__fx0],%%rax		\n\t"\
			"addpd	%%xmm2,%%xmm2		\n\t"/* top 52 bits */\
			"addpd	%%xmm1,%%xmm1		\n\t"/* low 26 bits */\
			"/* If x > q, subtract q: */\n\t"\
			"movaps	0x30(%%rdx),%%xmm0	\n\t"/* qhi52 */\
			"movaps	%%xmm0,%%xmm3		\n\t"/* cpy of qhi */\
			"cmppd	$0x2,%%xmm2,%%xmm3	\n\t"/* bitmask = (qhi <= xhi)? */\
			"andpd	%%xmm3,%%xmm0		\n\t"/* qhi52 & bitmask */\
			"andpd	0x10(%%rax),%%xmm3	\n\t"/* qlo26 & bitmask */\
			"subpd	%%xmm0,%%xmm2		\n\t"/* x % q, top 52 bits */\
			"subpd	%%xmm3,%%xmm1		\n\t"/* x % q, low 26 bits */\
		"twopmodq78_3wdq2_gcc64:	\n\t"\
		/* } */\
			"/* Normalize the result: */\n\t"\
			"movaps	     (%%rbx),%%xmm3	\n\t"/* two26i */\
			"mulpd	%%xmm3,%%xmm1		\n\t"/* xlo *= two26i */\
			"mulpd	%%xmm3,%%xmm2		\n\t"/* xhi *= two26i */\
			"movaps	%%xmm1,%%xmm0		\n\t"/* Init: fcy = cpy of ~xlo */\
			"addpd	%%xmm7,%%xmm1		\n\t"\
			"subpd	%%xmm7,%%xmm1		\n\t"\
			"subpd	%%xmm1,%%xmm0		\n\t"/* fx0 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fx0 *= two26f */\
			"mulpd	%%xmm3,%%xmm1		\n\t"/* fcy *= two26i */\
			"addpd	%%xmm1,%%xmm2		\n\t"/* Add carry into xhi */\
			"movaps	%%xmm2,%%xmm1		\n\t"/* Init: fcy = cpy of ~xhi */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"/* fx2 = fcy (no divide by 2^26) */\
			"subpd	%%xmm2,%%xmm1		\n\t"/* fx1 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm1	\n\t"/* fx1 *= two26f */\
			"movq	%[__fx0],%%rax		\n\t"\
			"movaps	%%xmm0,    (%%rax)	\n\t"\
			"movaps	%%xmm1,0x10(%%rax)	\n\t"\
			"movaps	%%xmm2,0x20(%%rax)	\n\t"\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			 ,[__fx0]	 "m" (fx0)		\
			 ,[__fq0]	 "m" (fq0)		\
			 ,[__pshift] "m" (pshift)	\
			 ,[__j]		 "m" (j)		\
			: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm7"	/* Clobbered registers */\
		);
		}	/* for(j...) */

	  #endif	/* OS_BITS */
		/*
		// Branchless version of the conditional doubling used in the preceding macro: //\
			"xorq	$0x1,%%rax		\n\t"\
			"movd	%%rax,%%xmm0		\n\t"\
			"pshufd	$0,%%xmm0,%%xmm0	\n\t"\
			"cvtdq2pd	%%xmm0,%%xmm0	\n\t"\
			"xorpd	%%xmm3,%%xmm3		\n\t"\
			"cmppd	$0,%%xmm3,%%xmm0	\n\t"\
			"movq	%[__fx0],%%rax		\n\t"\
		// Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in 0x30(%%rdx),0x10(%%rax): //\
			"movaps	%%xmm1,%%xmm3		\n\t"// cpy of xlo //\
			"andpd	%%xmm0,%%xmm3		\n\t"// xlo_ & bitmask //\
			"andpd	%%xmm2,%%xmm0		\n\t"// xhi_ & bitmask, overwrite bitmask with result //\
			"addpd	%%xmm3,%%xmm1		\n\t"// low 26 bits //\
			"addpd	%%xmm0,%%xmm2		\n\t"// top 52 bits //\
			"// If x > q, subtract q: //\n\t"\
			"movaps	0x30(%%rdx),%%xmm0	\n\t"// qhi52 //\
			"movaps	%%xmm0,%%xmm3		\n\t"// cpy of qhi //\
			"cmppd	$0x2,%%xmm2,%%xmm3	\n\t"// bitmask = (qhi <= xhi)? //\
			"andpd	%%xmm3,%%xmm0		\n\t"// qhi52 & bitmask //\
			"andpd	0x10(%%rax),%%xmm3	\n\t"// qlo26 & bitmask //\
			"subpd	%%xmm0,%%xmm2		\n\t"// x % q, top 52 bits //\
			"subpd	%%xmm3,%%xmm1		\n\t"// x % q, low 26 bits //\
		*/

	#endif	/* USE_SSE2 */

	#ifdef DEBUG_SSE2
	//	exit(0);
	#endif

		/*...Double and return.	These are specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		*/
	#ifdef USE_SSE2
		CVT78_3WORD_DOUBLE_UINT96(*fx0,*fx1,*fx2, x0);
		CVT78_3WORD_DOUBLE_UINT96(*gx0,*gx1,*gx2, x1);
	  #if FAC_DEBUG
		if(dbg) {
				printf("fx = %20.15f, %20.15f, %20.15f\n", *fx0, *fx1, *fx2);
				printf("gx = %20.15f, %20.15f, %20.15f\n", *gx0, *gx1, *gx2);
		}
	  #endif
	#else
		CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
		CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
	  #if FAC_DEBUG
		if(dbg) {
			printf("fq = %20.15f, %20.15f, %20.15f\n",fq0,fq1,fq2);
			printf("gq = %20.15f, %20.15f, %20.15f\n",gq0,gq1,gq2);
		}
	  #endif
	#endif
		ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
		ADD96(x1,x1,x1);
		q0.d0 -= FERMAT;
		q1.d0 -= FERMAT;
		SUB96(x0,q0,x0);
		SUB96(x1,q1,x1);

		tmp0 = CMPEQ96(x0, ONE96);
		tmp1 = CMPEQ96(x1, ONE96);
		r = tmp0;
		r += tmp1 << 1;
	#if FAC_DEBUG
		if(dbg) {
			printf("xout0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x0)]);
			printf("xout1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x1)]);
		}
	#endif
		return r;
	}

	/*** 4-trial-factor version ***/
	uint64 twopmodq78_3WORD_DOUBLE_q4(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3
			, int init_sse2, int thr_id
	)
	{
		const char func[] = "twopmodq78_3WORD_DOUBLE_q4";
		 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
		uint32 q32_0, qinv32_0, tmp32_0
			 , q32_1, qinv32_1, tmp32_1
			 , q32_2, qinv32_2, tmp32_2
			 , q32_3, qinv32_3, tmp32_3;
		uint64 lead7, tmp0, tmp1, tmp2, tmp3, r;
		uint96 q0, qinv0, qhalf0, x0, lo0
			 , q1, qinv1, qhalf1, x1, lo1
			 , q2, qinv2, qhalf2, x2, lo2
			 , q3, qinv3, qhalf3, x3, lo3;
		uint192 prod192;
		static uint64 psave = 0, pshift;
		static uint32 start_index, zshift, first_entry = TRUE;
		uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		const uint64 ihalf = 0x3FDfffffffffffffull;	/* Bitfield storing 0.5*(1-epsilon) in IEEE64 format */
		static int max_threads = 1;	// Default local-array-init is for just a single thread ... caller can re-init for > 1 threads later, if desired.
		static double *sc_arr = 0x0, *sc_ptr;
		vec_dbl *tmp;
		double dtmp;
	  #ifdef MULTITHREAD
		static double *__r0;	// Base address for discrete per-thread local stores ...
								// *NOTE*: more convenient to use double* rather than vec_dbl* here
	  #else
		static	// Following set of pointers only static-izable in single-thread mode
	  #endif
		double *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2, *flo0,*flo1,*flo2, *fhi0,*fhi1,*fhi2
			, *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2, *glo0,*glo1,*glo2, *ghi0,*ghi1,*ghi2
			, *hq0,*hq1,*hq2,*hqhi52, *hqinv0,*hqinv1,*hqinv2, *hx0,*hx1,*hx2, *hlo0,*hlo1,*hlo2, *hhi0,*hhi1,*hhi2
			, *iq0,*iq1,*iq2,*iqhi52, *iqinv0,*iqinv1,*iqinv2, *ix0,*ix1,*ix2, *ilo0,*ilo1,*ilo2, *ihi0,*ihi1,*ihi2
			, *half, *two26f, *two26i, *two13i, *sse2_rnd;

	#else

		double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2,fqhalf0, fx0,fx1,fx2, flo0,flo1,flo2,flohi52, fhi0;
		double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2,fqhalf1, gx0,gx1,gx2, glo0,glo1,glo2,glohi52, ghi0;
		double hq0,hq1,hq2, hqinv0,hqinv1,hqinv2,fqhalf2, hx0,hx1,hx2, hlo0,hlo1,hlo2,hlohi52, hhi0;
		double iq0,iq1,iq2, iqinv0,iqinv1,iqinv2,fqhalf3, ix0,ix1,ix2, ilo0,ilo1,ilo2,ilohi52, ihi0;
	  #if TRUNC_MUL78	// Streamlined version of modmul sequence, leaving high 52 bits in one double whenever feasible
		double fq_or_nil_lo26[2], fq_or_nil_hi52[2];
		double gq_or_nil_lo26[2], gq_or_nil_hi52[2];
		double hq_or_nil_lo26[2], hq_or_nil_hi52[2];
		double iq_or_nil_lo26[2], iq_or_nil_hi52[2];
	  #else
		double fhi1,fhi2;
		double ghi1,ghi2;
		double hhi1,hhi2;
		double ihi1,ihi2;
	  #endif
		int fidx,gidx,hidx,iidx;
		if(init_sse2) {
			ASSERT(init_sse2 <= 1, "Multithreading currently only supported for SIMD builds!");
			return 0;	// In non-SIMD mode, ||-init call is a no-op
		}

	#endif

		if(p != psave)
		{
		//	first_entry = FALSE;
			psave  = p;
			pshift = p + 78;
			j = leadz64(pshift);
			/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 77: */
			lead7 = ((pshift<<j) >> 57);
			if(lead7 > 77)
			{
				lead7 >>= 1;
				start_index =  64-j-6;	/* Use only the leftmost 6 bits */
			}
			else
				start_index =  64-j-7;

			zshift = 77 - lead7;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

			pshift = ~pshift;
		}

	#ifdef USE_SSE2

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
			// Alloc the local-memory block:
			sc_arr = ALLOC_DOUBLE(sc_arr, 0x50*max_threads + 4);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);	// Force vec_dbl-alignment
			ASSERT(((uintptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		#ifdef MULTITHREAD
			__r0  = sc_ptr;
			two13i = sc_ptr + 0x40;
			two26f = sc_ptr + 0x42;
			two26i = sc_ptr + 0x44;
			sse2_rnd=sc_ptr + 0x46;
			half   = sc_ptr + 0x48;
			dtmp = *(double*)&ihalf;
			for(j = 0; j < max_threads; ++j) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT_2((vec_dbl*)two13i  , TWO13FLINV);
				VEC_DBL_INIT_2((vec_dbl*)two26f  , TWO26FLOAT);
				VEC_DBL_INIT_2((vec_dbl*)two26i  , TWO26FLINV);
				VEC_DBL_INIT_2((vec_dbl*)sse2_rnd, crnd      );
				VEC_DBL_INIT_2((vec_dbl*)half    , dtmp      );
				// Move on to next thread's local store:
				two13i   += 0x50;
				two26f   += 0x50;
				two26i   += 0x50;
				sse2_rnd += 0x50;
				half     += 0x50;
			}
		#else
			/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 2 to span an SSE register: */
			fq0    = sc_ptr + 0x00;		gq0    = sc_ptr + 0x01;		hq0    = sc_ptr + 0x02;		iq0    = sc_ptr + 0x03;
			fq1    = sc_ptr + 0x04;		gq1    = sc_ptr + 0x05;		hq1    = sc_ptr + 0x06;		iq1    = sc_ptr + 0x07;
			fq2    = sc_ptr + 0x08;		gq2    = sc_ptr + 0x09;		hq2    = sc_ptr + 0x0a;		iq2    = sc_ptr + 0x0b;
			fqhi52 = sc_ptr + 0x0c;		gqhi52 = sc_ptr + 0x0d;		hqhi52 = sc_ptr + 0x0e;		iqhi52 = sc_ptr + 0x0f;
			fqinv0 = sc_ptr + 0x10;		gqinv0 = sc_ptr + 0x11;		hqinv0 = sc_ptr + 0x12;		iqinv0 = sc_ptr + 0x13;
			fqinv1 = sc_ptr + 0x14;		gqinv1 = sc_ptr + 0x15;		hqinv1 = sc_ptr + 0x16;		iqinv1 = sc_ptr + 0x17;
			fqinv2 = sc_ptr + 0x18;		gqinv2 = sc_ptr + 0x19;		hqinv2 = sc_ptr + 0x1a;		iqinv2 = sc_ptr + 0x1b;
			fx0    = sc_ptr + 0x1c;		gx0    = sc_ptr + 0x1d;		hx0    = sc_ptr + 0x1e;		ix0    = sc_ptr + 0x1f;
			fx1    = sc_ptr + 0x20;		gx1    = sc_ptr + 0x21;		hx1    = sc_ptr + 0x22;		ix1    = sc_ptr + 0x23;
			fx2    = sc_ptr + 0x24;		gx2    = sc_ptr + 0x25;		hx2    = sc_ptr + 0x26;		ix2    = sc_ptr + 0x27;
			flo0   = sc_ptr + 0x28;		glo0   = sc_ptr + 0x29;		hlo0   = sc_ptr + 0x2a;		ilo0   = sc_ptr + 0x2b;
			flo1   = sc_ptr + 0x2c;		glo1   = sc_ptr + 0x2d;		hlo1   = sc_ptr + 0x2e;		ilo1   = sc_ptr + 0x2f;
			flo2   = sc_ptr + 0x30;		glo2   = sc_ptr + 0x31;		hlo2   = sc_ptr + 0x32;		ilo2   = sc_ptr + 0x33;
			fhi0   = sc_ptr + 0x34;		ghi0   = sc_ptr + 0x35;		hhi0   = sc_ptr + 0x36;		ihi0   = sc_ptr + 0x37;
			fhi1   = sc_ptr + 0x38;		ghi1   = sc_ptr + 0x39;		hhi1   = sc_ptr + 0x3a;		ihi1   = sc_ptr + 0x3b;
			fhi2   = sc_ptr + 0x3c;		ghi2   = sc_ptr + 0x3d;		hhi2   = sc_ptr + 0x3e;		ihi2   = sc_ptr + 0x3f;
			two13i = sc_ptr + 0x40;
			two26f = sc_ptr + 0x42;
			two26i = sc_ptr + 0x44;
			sse2_rnd=sc_ptr + 0x46;
			half   = sc_ptr + 0x48;
			/* Can premultiply each of the multiword-mul inputs by 1/sqrt(2^13) due to quadraticity: */
			*two13i++ = TWO13FLINV;		*two13i-- = TWO13FLINV;
			*two26f++ = TWO26FLOAT;		*two26f-- = TWO26FLOAT;
			*two26i++ = TWO26FLINV;		*two26i-- = TWO26FLINV;
			/* SSE2 math = 53-mantissa-bit IEEE double-float: */
			*sse2_rnd++ = crnd;		*sse2_rnd-- = crnd;
			/* We init "half" = 0.5-epsilon here, because emulating FLOOR(x) via DNINT(x-half) requires
			us to always round up if x is a whole number, and our DNINT emulation can round either way if fractional part = 0.5:
			*/
			dtmp = *(double*)&ihalf;
			*half++ = dtmp;			*half-- = dtmp;
		#endif
			if(init_sse2) return 0;
		}	/* end of inits */

		/* If multithreaded, set the local-store pointers needed for the current thread; */
	  #ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		sc_ptr = __r0 + thr_id*0x50;
		/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 2 to span an SSE register: */
		fq0    = sc_ptr + 0x00;		gq0    = sc_ptr + 0x01;		hq0    = sc_ptr + 0x02;		iq0    = sc_ptr + 0x03;
		fq1    = sc_ptr + 0x04;		gq1    = sc_ptr + 0x05;		hq1    = sc_ptr + 0x06;		iq1    = sc_ptr + 0x07;
		fq2    = sc_ptr + 0x08;		gq2    = sc_ptr + 0x09;		hq2    = sc_ptr + 0x0a;		iq2    = sc_ptr + 0x0b;
		fqhi52 = sc_ptr + 0x0c;		gqhi52 = sc_ptr + 0x0d;		hqhi52 = sc_ptr + 0x0e;		iqhi52 = sc_ptr + 0x0f;
		fqinv0 = sc_ptr + 0x10;		gqinv0 = sc_ptr + 0x11;		hqinv0 = sc_ptr + 0x12;		iqinv0 = sc_ptr + 0x13;
		fqinv1 = sc_ptr + 0x14;		gqinv1 = sc_ptr + 0x15;		hqinv1 = sc_ptr + 0x16;		iqinv1 = sc_ptr + 0x17;
		fqinv2 = sc_ptr + 0x18;		gqinv2 = sc_ptr + 0x19;		hqinv2 = sc_ptr + 0x1a;		iqinv2 = sc_ptr + 0x1b;
		fx0    = sc_ptr + 0x1c;		gx0    = sc_ptr + 0x1d;		hx0    = sc_ptr + 0x1e;		ix0    = sc_ptr + 0x1f;
		fx1    = sc_ptr + 0x20;		gx1    = sc_ptr + 0x21;		hx1    = sc_ptr + 0x22;		ix1    = sc_ptr + 0x23;
		fx2    = sc_ptr + 0x24;		gx2    = sc_ptr + 0x25;		hx2    = sc_ptr + 0x26;		ix2    = sc_ptr + 0x27;
		flo0   = sc_ptr + 0x28;		glo0   = sc_ptr + 0x29;		hlo0   = sc_ptr + 0x2a;		ilo0   = sc_ptr + 0x2b;
		flo1   = sc_ptr + 0x2c;		glo1   = sc_ptr + 0x2d;		hlo1   = sc_ptr + 0x2e;		ilo1   = sc_ptr + 0x2f;
		flo2   = sc_ptr + 0x30;		glo2   = sc_ptr + 0x31;		hlo2   = sc_ptr + 0x32;		ilo2   = sc_ptr + 0x33;
		fhi0   = sc_ptr + 0x34;		ghi0   = sc_ptr + 0x35;		hhi0   = sc_ptr + 0x36;		ihi0   = sc_ptr + 0x37;
		fhi1   = sc_ptr + 0x38;		ghi1   = sc_ptr + 0x39;		hhi1   = sc_ptr + 0x3a;		ihi1   = sc_ptr + 0x3b;
		fhi2   = sc_ptr + 0x3c;		ghi2   = sc_ptr + 0x3d;		hhi2   = sc_ptr + 0x3e;		ihi2   = sc_ptr + 0x3f;
		two13i = sc_ptr + 0x40;
		two26f = sc_ptr + 0x42;
		two26i = sc_ptr + 0x44;
		sse2_rnd=sc_ptr + 0x46;
		half   = sc_ptr + 0x48;
		tmp = (vec_dbl*)sse2_rnd; ASSERT((tmp->d0 == crnd) && (tmp->d1 == crnd), "Bad data at sse2_rnd address!");
	  #endif

	#endif

		ASSERT((p >> 63) == 0, "twopmodq78_q4 : p must be < 2^63!");
		q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;
	#ifdef MUL_LOHI64_SUBROUTINE
		// MUL_LOHI64 expects a 64-bit high-part pointer, in 32bit builds this buggers us if we try dumping hi-part directly into 32-bit q.d1
		MUL_LOHI64(q0.d0, k0,&q0.d0,&tmp0);	q0.d1 = tmp0;
		MUL_LOHI64(q1.d0, k1,&q1.d0,&tmp0);	q1.d1 = tmp0;
		MUL_LOHI64(q2.d0, k2,&q2.d0,&tmp0);	q2.d1 = tmp0;
		MUL_LOHI64(q3.d0, k3,&q3.d0,&tmp0);	q3.d1 = tmp0;
	#else
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
		MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
		MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);
	#endif
		q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
		q1.d0 += 1;
		q2.d0 += 1;
		q3.d0 += 1;
		ASSERT((q0.d1 >> 14) == 0, "twopmodq78_q4 : (q0.d1 >> 14) != 0");
		ASSERT((q1.d1 >> 14) == 0, "twopmodq78_q4 : (q1.d1 >> 14) != 0");
		ASSERT((q2.d1 >> 14) == 0, "twopmodq78_q4 : (q2.d1 >> 14) != 0");
		ASSERT((q3.d1 >> 14) == 0, "twopmodq78_q4 : (q3.d1 >> 14) != 0");

		q32_0 = (uint32)q0.d0;
		q32_1 = (uint32)q1.d0;
		q32_2 = (uint32)q2.d0;
		q32_3 = (uint32)q3.d0;

		/* Convert q to floating form: */
	#ifdef USE_SSE2

		CVT_UINT78_3WORD_DOUBLE(q0 ,*fq0,*fq1,*fq2);
		CVT_UINT78_3WORD_DOUBLE(q1 ,*gq0,*gq1,*gq2);
		CVT_UINT78_3WORD_DOUBLE(q2 ,*hq0,*hq1,*hq2);
		CVT_UINT78_3WORD_DOUBLE(q3 ,*iq0,*iq1,*iq2);

	  #ifdef USE_ARM_V8_SIMD
		// No TF support on ARMv8, just supply a stub macro:
		__asm__ volatile (\
			"ldr x0,%[__fq0]	\n\t"\
			:					// outputs: none
			: [__fq0] "m" (fq0)	// All inputs from memory addresses here
			: "cc","memory","x0"	/* Clobbered registers */\
		);

	  #elif OS_BITS == 32

		#error 32-bit OSes no longer supported for SIMD builds!

	  #elif OS_BITS == 64

		__asm__ volatile (\
			"movq	%[__fq0],%%rax						\n\t"\
			"movq	%[__two26f],%%rsi					\n\t"\
			"movaps	    (%%rsi),%%xmm8	/* two26f */	\n\t"\
			"movaps	0x10(%%rsi),%%xmm9	/* two26i */	\n\t"\
			"movaps	    (%%rax),%%xmm0	/* fq0 */		\n\t	movaps	0x10(%%rax),%%xmm4	/* hq0 */		\n\t"\
			"movaps	0x20(%%rax),%%xmm1	/* fq1 */		\n\t	movaps	0x30(%%rax),%%xmm5	/* hq1 */		\n\t"\
			"movaps	0x40(%%rax),%%xmm2	/* fq2 */		\n\t	movaps	0x50(%%rax),%%xmm6	/* hq2 */		\n\t"\
			"movaps	%%xmm2,%%xmm3		/* cpy fg2 */	\n\t	movaps	%%xmm6,%%xmm7		/* cpy hi2 */	\n\t"\
			"mulpd	%%xmm8,%%xmm3						\n\t	mulpd	%%xmm8,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3	/* Hi 52 out bits */\n\t	addpd	%%xmm5,%%xmm7		/* Hi 52 */		\n\t"\
			"mulpd	%%xmm9,%%xmm0						\n\t	mulpd	%%xmm9,%%xmm4						\n\t"\
			"mulpd	%%xmm9,%%xmm1						\n\t	mulpd	%%xmm9,%%xmm5						\n\t"\
			"mulpd	%%xmm9,%%xmm2						\n\t	mulpd	%%xmm9,%%xmm6						\n\t"\
			"movaps	%%xmm0,    (%%rax)	/* fq0/2^26 */	\n\t	movaps	%%xmm4,0x10(%%rax)	/* hq0/2^26 */	\n\t"\
			"movaps	%%xmm1,0x20(%%rax)	/* fq1/2^26 */	\n\t	movaps	%%xmm5,0x30(%%rax)	/* hq1/2^26 */	\n\t"\
			"movaps	%%xmm2,0x40(%%rax)	/* fq2/2^26 */	\n\t	movaps	%%xmm6,0x50(%%rax)	/* hq2/2^26 */	\n\t"\
			"movaps	%%xmm3,0x60(%%rax)	/* fqhi52 */	\n\t	movaps	%%xmm7,0x70(%%rax)	/* hqhi52 */	\n\t"\
			:					/* outputs: none */\
			: [__fq0] "m" (fq0)	/* All inputs from memory addresses here */\
			 ,[__two26f] "m" (two26f)\
			: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
		);

	  #endif

	#else
		CVT_UINT78_3WORD_DOUBLE(q0, fq0,fq1,fq2);
		CVT_UINT78_3WORD_DOUBLE(q1, gq0,gq1,gq2);
		CVT_UINT78_3WORD_DOUBLE(q2, hq0,hq1,hq2);
		CVT_UINT78_3WORD_DOUBLE(q3, iq0,iq1,iq2);

		fq_or_nil_lo26[0] = 0.0; fq_or_nil_lo26[1] = fq0; fq_or_nil_hi52[0] = 0.0; fq_or_nil_hi52[1] = fq1 + TWO26FLOAT*fq2;
		gq_or_nil_lo26[0] = 0.0; gq_or_nil_lo26[1] = gq0; gq_or_nil_hi52[0] = 0.0; gq_or_nil_hi52[1] = gq1 + TWO26FLOAT*gq2;
		hq_or_nil_lo26[0] = 0.0; hq_or_nil_lo26[1] = hq0; hq_or_nil_hi52[0] = 0.0; hq_or_nil_hi52[1] = hq1 + TWO26FLOAT*hq2;
		iq_or_nil_lo26[0] = 0.0; iq_or_nil_lo26[1] = iq0; iq_or_nil_hi52[0] = 0.0; iq_or_nil_hi52[1] = iq1 + TWO26FLOAT*iq2;
	#endif

		RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
		RSHIFT_FAST96(q1, 1, qhalf1);
		RSHIFT_FAST96(q2, 1, qhalf2);
		RSHIFT_FAST96(q3, 1, qhalf3);

	#ifndef USE_SSE2
		fqhalf0 = qhalf0.d1*TWO64FLOAT + qhalf0.d0;
		fqhalf1 = qhalf1.d1*TWO64FLOAT + qhalf1.d0;
		fqhalf2 = qhalf2.d1*TWO64FLOAT + qhalf2.d0;
		fqhalf3 = qhalf3.d1*TWO64FLOAT + qhalf3.d0;
	#endif

		qinv32_0 = (q32_0 + q32_0 + q32_0) ^ (uint32)2;
		qinv32_1 = (q32_1 + q32_1 + q32_1) ^ (uint32)2;
		qinv32_2 = (q32_2 + q32_2 + q32_2) ^ (uint32)2;
		qinv32_3 = (q32_3 + q32_3 + q32_3) ^ (uint32)2;
		/* 4 bits: */
		tmp32_0 = q32_0*qinv32_0;
		tmp32_1 = q32_1*qinv32_1;
		tmp32_2 = q32_2*qinv32_2;
		tmp32_3 = q32_3*qinv32_3;
		qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);
		qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);
		qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);
		qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);
		/* 8 bits: */
		tmp32_0 = q32_0*qinv32_0;
		tmp32_1 = q32_1*qinv32_1;
		tmp32_2 = q32_2*qinv32_2;
		tmp32_3 = q32_3*qinv32_3;
		qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);
		qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);
		qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);
		qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);
		/* 16 bits: */
		tmp32_0 = q32_0*qinv32_0;
		tmp32_1 = q32_1*qinv32_1;
		tmp32_2 = q32_2*qinv32_2;
		tmp32_3 = q32_3*qinv32_3;
		qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);
		qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);
		qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);
		qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);
		/* 32 bits: */
		qinv0.d0 = (uint64)qinv32_0;
		qinv1.d0 = (uint64)qinv32_1;
		qinv2.d0 = (uint64)qinv32_2;
		qinv3.d0 = (uint64)qinv32_3;
		tmp0 = q0.d0*qinv0.d0;
		tmp1 = q1.d0*qinv1.d0;
		tmp2 = q2.d0*qinv2.d0;
		tmp3 = q3.d0*qinv3.d0;
		qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
		qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
		qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);
		qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);
		/* 64 bits: */
	#ifdef MUL_LOHI64_SUBROUTINE
		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
		qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + __MULH64(q2.d0, qinv2.d0));
		qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + __MULH64(q3.d0, qinv3.d0));
	#else
		MULH64(q0.d0, qinv0.d0, tmp0);
		MULH64(q1.d0, qinv1.d0, tmp1);
		MULH64(q2.d0, qinv2.d0, tmp2);
		MULH64(q3.d0, qinv3.d0, tmp3);

		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);
		qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + tmp2);
		qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + tmp3);
	#endif
		/* 64 bits: */
		qinv0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
		qinv1.d1 &= 0x0000000000003fff;
		qinv2.d1 &= 0x0000000000003fff;
		qinv3.d1 &= 0x0000000000003fff;

		/* Convert qinv to floating form: */
	#ifdef USE_SSE2

		CVT_UINT78_3WORD_DOUBLE(qinv0 ,*fqinv0,*fqinv1,*fqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv1 ,*gqinv0,*gqinv1,*gqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv2 ,*hqinv0,*hqinv1,*hqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv3 ,*iqinv0,*iqinv1,*iqinv2);

		#if defined(COMPILER_TYPE_MSVC) && (TRYQ == 4)

			#error twopmodq78_3WORD_DOUBLE_q4: 4-operand SSE2 modmul does not support MSVC inline assembler!

		#else	/* GCC-style inline ASM: */

		  #ifdef USE_ARM_V8_SIMD
			// No TF support on ARMv8, just supply a stub macro:
			__asm__ volatile (\
				"ldr x0,%[__fq0]	\n\t"\
				:					// outputs: none
				: [__fq0] "m" (fq0)	// All inputs from memory addresses here
				: "cc","memory","x0"	/* Clobbered registers */\
			);

		  #elif OS_BITS == 32

			#error 32-bit OSes no longer supported for SIMD builds!

		  #elif OS_BITS == 64

			__asm__ volatile (\
				"movq	%[__fqinv0],%%rax	\n\t"\
				"movq	%[__two26i],%%rsi	\n\t"\
				"movaps	(%%rsi),%%xmm6		\n\t"\
				"movaps	    (%%rax),%%xmm0	\n\t	movaps	0x10(%%rax),%%xmm3	\n\t"\
				"movaps	0x20(%%rax),%%xmm1	\n\t	movaps	0x30(%%rax),%%xmm4	\n\t"\
				"movaps	0x40(%%rax),%%xmm2	\n\t	movaps	0x50(%%rax),%%xmm5	\n\t"\
				"mulpd	%%xmm6,%%xmm0		\n\t	mulpd	%%xmm6,%%xmm3		\n\t"\
				"mulpd	%%xmm6,%%xmm1		\n\t	mulpd	%%xmm6,%%xmm4		\n\t"\
				"mulpd	%%xmm6,%%xmm2		\n\t	mulpd	%%xmm6,%%xmm5		\n\t"\
				"movaps	%%xmm0,    (%%rax)	\n\t	movaps	%%xmm3,0x10(%%rax)	\n\t"\
				"movaps	%%xmm1,0x20(%%rax)	\n\t	movaps	%%xmm4,0x30(%%rax)	\n\t"\
				"movaps	%%xmm2,0x40(%%rax)	\n\t	movaps	%%xmm5,0x50(%%rax)	\n\t"\
				:					/* outputs: none */\
				: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
				 ,[__two26i] "m" (two26i)	\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6"	/* Clobbered registers */\
			);

		  #endif

		#endif

	#else
		CVT_UINT78_3WORD_DOUBLE(qinv0, fqinv0,fqinv1,fqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv1, gqinv0,gqinv1,gqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv2, hqinv0,hqinv1,hqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv3, iqinv0,iqinv1,iqinv2);
	#endif

		/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
		j = start_index-1;

		/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
		LSHIFT96(qinv0, zshift, lo0);	lo0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
		LSHIFT96(qinv1, zshift, lo1);	lo1.d1 &= 0x0000000000003fff;
		LSHIFT96(qinv2, zshift, lo2);	lo2.d1 &= 0x0000000000003fff;
		LSHIFT96(qinv3, zshift, lo3);	lo3.d1 &= 0x0000000000003fff;

		MUL_LOHI96_PROD192(q0,lo0,prod192);	RSHIFT192(prod192,78,prod192);	lo0.d0 = prod192.d0;	lo0.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q1,lo1,prod192);	RSHIFT192(prod192,78,prod192);	lo1.d0 = prod192.d0;	lo1.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q2,lo2,prod192);	RSHIFT192(prod192,78,prod192);	lo2.d0 = prod192.d0;	lo2.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q3,lo3,prod192);	RSHIFT192(prod192,78,prod192);	lo3.d0 = prod192.d0;	lo3.d1 = prod192.d1;

		/* hi = 0 in this instance, which simplifies things. */
		SUB96(q0, lo0, x0);	/* Put the result in lo (rather than x), to ease overflow check below */
		SUB96(q1, lo1, x1);
		SUB96(q2, lo2, x2);
		SUB96(q3, lo3, x3);

		if((pshift >> j) & (uint32)1)
		{
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
			if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
			if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
			if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
		}

		/* Convert x to floating form: */
	#ifdef USE_SSE2

		CVT_UINT78_3WORD_DOUBLE(x0 ,*fx0,*fx1,*fx2);
		CVT_UINT78_3WORD_DOUBLE(x1 ,*gx0,*gx1,*gx2);
		CVT_UINT78_3WORD_DOUBLE(x2 ,*hx0,*hx1,*hx2);
		CVT_UINT78_3WORD_DOUBLE(x3 ,*ix0,*ix1,*ix2);

		#ifdef USE_ARM_V8_SIMD
			// No TF support on ARMv8
		#elif OS_BITS == 32
			// 32-bit loads these at start of each loop pass
		#elif OS_BITS == 64

		__asm__ volatile (\
			"movq	%[__fx0],%%rax	\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t	movaps	0x10(%%rax),%%xmm8	\n\t"\
			"movaps	0x20(%%rax),%%xmm2	\n\t	movaps	0x30(%%rax),%%xmm10	\n\t"\
			"movaps	0x40(%%rax),%%xmm4	\n\t	movaps	0x50(%%rax),%%xmm12	\n\t"\
			:				/* outputs: none */\
			: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
			: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
		);

		#endif

	#else
		CVT_UINT78_3WORD_DOUBLE(x0, fx0,fx1,fx2);
		CVT_UINT78_3WORD_DOUBLE(x1, gx0,gx1,gx2);
		CVT_UINT78_3WORD_DOUBLE(x2, hx0,hx1,hx2);
		CVT_UINT78_3WORD_DOUBLE(x3, ix0,ix1,ix2);
	#endif

		/*...x^2 mod q is returned in x. */
		/* All 3-word-double-form operands have components in the following size ranges:
			fword0,1 in [-2^25, +2^25]
			fword2   in [   -1, +2^26]
		*/
	#ifndef USE_SSE2

		for(j = start_index-2; j >= 0; j--)
		{
			/*...x^2 mod q is returned in x. */

		/********************************************************************************************************/
		#ifdef TRUNC_MUL78	// Streamlined version of modmul sequence, leaving high 52 bits in one double whenever feasible
		/********************************************************************************************************/

			/* SQR_LOHI96(x,lo,hi); */
			SQR_LOHI78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(
				  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fx1
				, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,gx1
				, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hx1
				, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ix1
			);
			/* MULL96(lo,qinv,lo); */
			MULL78_3WORD_DOUBLE_q4(
				  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
				, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
				, hlo0,hlo1,hlo2, hqinv0,hqinv1,hqinv2, hlo0,hlo1,hlo2
				, ilo0,ilo1,ilo2, iqinv0,iqinv1,iqinv2, ilo0,ilo1,ilo2
			);
			/* MULH96(q,lo,lo); */
			MULH78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(
				  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flohi52
				, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glohi52
				, hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlohi52
				, iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilohi52
			);

			/* If h < l, then calculate q+(h-l) < q; otherwise calculate h-l. */
			fidx = (fx1 < flohi52);
			gidx = (gx1 < glohi52);
			hidx = (hx1 < hlohi52);
			iidx = (ix1 < ilohi52);

			fx1 -= flohi52;
			gx1 -= glohi52;
			hx1 -= hlohi52;
			ix1 -= ilohi52;

			fhi0 -= flo0;
			ghi0 -= glo0;
			hhi0 -= hlo0;
			ihi0 -= ilo0;

			fx1 += fq_or_nil_hi52[fidx];
			gx1 += gq_or_nil_hi52[gidx];
			hx1 += hq_or_nil_hi52[hidx];
			ix1 += iq_or_nil_hi52[iidx];

			fx0 = fhi0 + fq_or_nil_lo26[fidx];
			gx0 = ghi0 + gq_or_nil_lo26[gidx];
			hx0 = hhi0 + hq_or_nil_lo26[hidx];
			ix0 = ihi0 + iq_or_nil_lo26[iidx];
			// Normalize result and store in x-vector...use hi0 terms as carries:
			fhi0 = DNINT(fx0*TWO26FLINV);
			ghi0 = DNINT(gx0*TWO26FLINV);
			hhi0 = DNINT(hx0*TWO26FLINV);
			ihi0 = DNINT(ix0*TWO26FLINV);

			fx0 -= fhi0*TWO26FLOAT;
			gx0 -= ghi0*TWO26FLOAT;
			hx0 -= hhi0*TWO26FLOAT;
			ix0 -= ihi0*TWO26FLOAT;

			fx1 += fhi0;
			gx1 += ghi0;
			hx1 += hhi0;
			ix1 += ihi0;

			/* Branchless version of the double-and-subtract-q-if-(2x >= q) sequence:
				x = x + x - ((-(x > qhalf)) & q);
			In FP version replace integer and with array lookup.
			*/
			if((pshift >> j) & (uint64)1)
			{
				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				/* Fill FP approximation to x with as many SBs as possible prior to comparing */
				/* TRUE here means will need to subtract q from 2x: */
				fidx = (fqhalf0 < fx1*TWO26FLOAT + fx0);
				gidx = (fqhalf1 < gx1*TWO26FLOAT + gx0);
				hidx = (fqhalf2 < hx1*TWO26FLOAT + hx0);
				iidx = (fqhalf3 < ix1*TWO26FLOAT + ix0);

				fx1 *= 2;	fx0 += fx0;
				gx1 *= 2;	gx0 += gx0;
				hx1 *= 2;	hx0 += hx0;
				ix1 *= 2;	ix0 += ix0;

				fx1 -= fq_or_nil_hi52[fidx];
				gx1 -= gq_or_nil_hi52[gidx];
				hx1 -= hq_or_nil_hi52[hidx];
				ix1 -= iq_or_nil_hi52[iidx];

				fx0     -= fq_or_nil_lo26[fidx];
				gx0     -= gq_or_nil_lo26[gidx];
				hx0     -= hq_or_nil_lo26[hidx];
				ix0     -= iq_or_nil_lo26[iidx];
			}

			/* Normalize the result - use currently-unused x2 coeff as carry: */
			/* Digit 0: */
			fx2 = DNINT(fx0*TWO26FLINV);
			gx2 = DNINT(gx0*TWO26FLINV);
			hx2 = DNINT(hx0*TWO26FLINV);
			ix2 = DNINT(ix0*TWO26FLINV);

			fx0 -= fx2*TWO26FLOAT;
			gx0 -= gx2*TWO26FLOAT;
			hx0 -= hx2*TWO26FLOAT;
			ix0 -= ix2*TWO26FLOAT;

			/* Digit 1: */
			fx1 += fx2;
			gx1 += gx2;
			hx1 += hx2;
			ix1 += ix2;

			fx2 = DNINT(fx1*TWO26FLINV);
			gx2 = DNINT(gx1*TWO26FLINV);
			hx2 = DNINT(hx1*TWO26FLINV);
			ix2 = DNINT(ix1*TWO26FLINV);

			fx1 -= fx2*TWO26FLOAT;
			gx1 -= gx2*TWO26FLOAT;
			hx1 -= hx2*TWO26FLOAT;
			ix1 -= ix2*TWO26FLOAT;

			/* Digit 2 already in x2 term. */

		/********************************************************************************************************/
		#else	// Basic impl of 3-word-double-based 78-bit-integer modmul:
		/********************************************************************************************************/

			/*...x^2 mod q is returned in x. */
			/* SQR_LOHI96(x,lo,hi); */
			SQR_LOHI78_3WORD_DOUBLE_q4(
				  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2
				, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2
				, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hhi1,hhi2
				, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ihi1,ihi2
			);
			/* MULL96(lo,qinv,lo); */
			MULL78_3WORD_DOUBLE_q4(
				  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
				, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
				, hlo0,hlo1,hlo2, hqinv0,hqinv1,hqinv2, hlo0,hlo1,hlo2
				, ilo0,ilo1,ilo2, iqinv0,iqinv1,iqinv2, ilo0,ilo1,ilo2
			);
			/* MULH78(q,lo,lo); */
			MULH78_3WORD_DOUBLE_q4(
				  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2
				, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2
				, hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlo1,hlo2
				, iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilo1,ilo2
			);

			/* If h < l, then calculate q+(h-l) < q; otherwise calculate h-l. */
			if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
			{
				SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
				ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
			}
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

			if(CMPLT78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2))
			{
				SUB78_3WORD_DOUBLE(gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2);
				ADD78_3WORD_DOUBLE(glo0,glo1,glo2, ghi0,ghi1,ghi2, gx0,gx1,gx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2, gx0,gx1,gx2);
			}
			NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

			if(CMPLT78_3WORD_DOUBLE(hhi0,hhi1,hhi2, hlo0,hlo1,hlo2))
			{
				SUB78_3WORD_DOUBLE(hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlo1,hlo2);
				ADD78_3WORD_DOUBLE(hlo0,hlo1,hlo2, hhi0,hhi1,hhi2, hx0,hx1,hx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(hhi0,hhi1,hhi2, hlo0,hlo1,hlo2, hx0,hx1,hx2);
			}
			NORMALIZE78_3WORD_DOUBLE(hx0,hx1,hx2);

			if(CMPLT78_3WORD_DOUBLE(ihi0,ihi1,ihi2, ilo0,ilo1,ilo2))
			{
				SUB78_3WORD_DOUBLE(iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilo1,ilo2);
				ADD78_3WORD_DOUBLE(ilo0,ilo1,ilo2, ihi0,ihi1,ihi2, ix0,ix1,ix2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(ihi0,ihi1,ihi2, ilo0,ilo1,ilo2, ix0,ix1,ix2);
			}
			NORMALIZE78_3WORD_DOUBLE(ix0,ix1,ix2);

			if((pshift >> j) & (uint64)1)
			{
				/* Convert fx to uint96 for purpose of this particular comparison: */
				CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
				CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
				CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
				CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x0, qhalf0))
				{
					ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
					SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				}
				NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x1, qhalf1))
				{
					ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
					SUB78_3WORD_DOUBLE(gx0,gx1,gx2, gq0,gq1,gq2, gx0,gx1,gx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
				}
				NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x2, qhalf2))
				{
					ADD78_3WORD_DOUBLE(hx0,hx1,hx2, hx0,hx1,hx2, hx0,hx1,hx2);
					SUB78_3WORD_DOUBLE(hx0,hx1,hx2, hq0,hq1,hq2, hx0,hx1,hx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(hx0,hx1,hx2, hx0,hx1,hx2, hx0,hx1,hx2);
				}
				NORMALIZE78_3WORD_DOUBLE(hx0,hx1,hx2);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x3, qhalf3))
				{
					ADD78_3WORD_DOUBLE(ix0,ix1,ix2, ix0,ix1,ix2, ix0,ix1,ix2);
					SUB78_3WORD_DOUBLE(ix0,ix1,ix2, iq0,iq1,iq2, ix0,ix1,ix2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(ix0,ix1,ix2, ix0,ix1,ix2, ix0,ix1,ix2);
				}
				NORMALIZE78_3WORD_DOUBLE(ix0,ix1,ix2);
			}

		#endif	/* #if(1) */
	//printf("j = %2d: x = %20.10f,%20.10f,%20.10f\n",j,fx0,fx1,fx2);
		}
	//exit(0);

	#elif defined(USE_SSE2)	/* Inner loop body needs 42 movaps, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */

		#ifdef USE_ARM_V8_SIMD
			// No TF support on ARMv8
		#elif OS_BITS == 32

			#error 32-bit OSes no longer supported for SIMD builds!

		#elif OS_BITS == 64	// 16 XMM registers, but cannot assumee a proper ROUNDPD instruction which we can use
							// for FLOOR, CEILING, DNINT without burning a register to hold any "magic constant",
							// since that was added only a few x86_64 chip generations in. (I did test it on
							// several systems supporting ROUNDPD, though, and using +- RND runs ~same speed here).

		  #define ASM_LOOP 0

		  #if ASM_LOOP	/* Pure-ASM loop control, and branchless mod-doubling segment: */
			__asm__ volatile (\
				"movslq	%[__start_index], %%rcx		\n\t"\
				"subq $2,%%rcx						\n\t"\
				"test %%rcx, %%rcx					\n\t"\
				"jl Loop4End		/* Skip if n < 0 */	\n\t"\
			"Loop4Beg:								\n\t"\
			/* SQR_LOHI78_3WORD_DOUBLE_q4(fx, flo,fhi). Inputs flo0,1,2 enter in xmm0,2,4 [lcol] and xmm8,10,12 [rcol]: */\
				"/* SQR_LOHI78_3WORD_DOUBLE_q4(): */\n\t"\
				"movq	%[__fx0],%%rax				\n\t"\
				"movq	%[__two26i],%%rbx			\n\t"\
				"movaps	-0x20(%%rbx),%%xmm6	/* two13i */	\n\t"\
			/* fx0,1,2 assumed in xmm0,2,4 on loop entry */\
				"mulpd	%%xmm6,%%xmm0			\n\t	mulpd	%%xmm6 ,%%xmm8			\n\t"\
				"mulpd	%%xmm6,%%xmm2			\n\t	mulpd	%%xmm6 ,%%xmm10			\n\t"\
				"mulpd	%%xmm6,%%xmm4			\n\t	mulpd	%%xmm6 ,%%xmm12			\n\t"\
				"movaps	0x10(%%rbx),%%xmm7		\n\t"/* xmm7 = rnd_const, shared between both columns */\
				"movaps	%%xmm0,%%xmm1			\n\t	movaps	%%xmm8 ,%%xmm9			\n\t"\
				"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"\
				"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9 ,%%xmm9			\n\t"\
				"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11			\n\t"\
				"movaps	%%xmm1,%%xmm5			\n\t	movaps	%%xmm9 ,%%xmm13			\n\t"\
				"mulpd	%%xmm0,%%xmm0			\n\t	mulpd	%%xmm8 ,%%xmm8			\n\t"\
				"mulpd	%%xmm2,%%xmm1			\n\t	mulpd	%%xmm10,%%xmm9			\n\t"\
				"mulpd	%%xmm2,%%xmm2			\n\t	mulpd	%%xmm10,%%xmm10			\n\t"\
				"mulpd	%%xmm4,%%xmm5			\n\t	mulpd	%%xmm12,%%xmm13			\n\t"\
				"mulpd	%%xmm4,%%xmm3			\n\t	mulpd	%%xmm12,%%xmm11			\n\t"\
				"mulpd	%%xmm4,%%xmm4			\n\t	mulpd	%%xmm12,%%xmm12			\n\t"\
				/* Move this part of Digit 2 computation here to free up xmm5,13: */\
				"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10			\n\t"\
				/* Digit 0: */\
				"movaps (%%rbx),%%xmm15			\n\t	movaps	-0x10(%%rbx),%%xmm5		\n\t"\
				"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
				"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
				"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5,%%xmm8		\n\t"\
				"mulpd	%%xmm15,%%xmm6		\n\t	mulpd	     %%xmm15,%%xmm14	\n\t"\
				/* Digit 1: */\
				"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
				"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm6,%%xmm1			\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
				"mulpd	%%xmm5,%%xmm1		\n\t	mulpd	%%xmm5,%%xmm9		\n\t"\
				"mulpd	     %%xmm15,%%xmm6		\n\t	mulpd	     %%xmm15,%%xmm14	\n\t"\
				/* Digit 2: */\
				"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
				"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t"\
				"subpd	0x20(%%rbx),%%xmm6		\n\t	subpd	0x20(%%rbx),%%xmm14		\n\t"\
				"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm6,%%xmm2			\n\t	subpd	%%xmm14,%%xmm10			\n\t"\
				"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t"\
				"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
				/* Digit 3: */\
				"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
				"movaps	%%xmm3,%%xmm6			\n\t	movaps	%%xmm11,%%xmm14			\n\t"\
				"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm6,%%xmm3			\n\t	subpd	%%xmm14,%%xmm11			\n\t"\
				"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t"\
				"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
				/* Digit 4: */\
				"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t"\
				"movaps	%%xmm4,%%xmm6			\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
				"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm6,%%xmm4			\n\t	subpd	%%xmm14,%%xmm12			\n\t"\
				"mulpd	%%xmm5,%%xmm4			\n\t	mulpd	%%xmm5 ,%%xmm12			\n\t"\
				/* Digit 5 = the carry. flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6 */\
				"mulpd	%%xmm5,%%xmm6			\n\t	mulpd	%%xmm5 ,%%xmm14			\n\t"\
				"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t"\
				"movaps	%%xmm3,0xc0(%%rax)		\n\t	movaps	%%xmm11,0xd0(%%rax)		\n\t"\
				"movaps	%%xmm4,0xe0(%%rax)		\n\t	movaps	%%xmm12,0xf0(%%rax)		\n\t"\
			/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */\
				"movq	%[__fqinv0],%%rdx		\n\t"\
				/* Digit 0: */\
				"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t"\
				"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t"\
				"mulpd	    (%%rdx),%%xmm0		\n\t	mulpd	0x10(%%rdx),%%xmm8		\n\t"\
				"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
				"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
				"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t"\
				"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
				/* Digit 1: */\
				"mulpd	0x20(%%rdx),%%xmm3		\n\t	mulpd	0x30(%%rdx),%%xmm11		\n\t"\
				"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
				"mulpd	    (%%rdx),%%xmm1		\n\t	mulpd	0x10(%%rdx),%%xmm9		\n\t"\
				"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
				"movaps	%%xmm1,%%xmm3			\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
				"addpd	%%xmm7,%%xmm3			\n\t	addpd	%%xmm7 ,%%xmm11			\n\t"\
				"subpd	%%xmm7,%%xmm3			\n\t	subpd	%%xmm7 ,%%xmm11			\n\t"\
				"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t"\
				"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t"\
				"mulpd	%%xmm15,%%xmm3			\n\t	mulpd	%%xmm15,%%xmm11			\n\t"\
				/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\
				"mulpd	    (%%rdx),%%xmm2		\n\t	mulpd	0x10(%%rdx),%%xmm10		\n\t"\
				"mulpd	0x20(%%rdx),%%xmm6		\n\t	mulpd	0x30(%%rdx),%%xmm14		\n\t"\
				"mulpd	0x40(%%rdx),%%xmm4		\n\t	mulpd	0x50(%%rdx),%%xmm12		\n\t"\
				"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
				"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11			\n\t"\
				"addpd	%%xmm3,%%xmm2			\n\t	addpd	%%xmm11,%%xmm10			\n\t"\
				"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"\
				"subpd	0x20(%%rbx),%%xmm3		\n\t	subpd	0x20(%%rbx),%%xmm11		\n\t"\
				"addpd	%%xmm7,%%xmm3			\n\t	addpd	%%xmm7 ,%%xmm11			\n\t"\
				"subpd	%%xmm7,%%xmm3			\n\t	subpd	%%xmm7 ,%%xmm11			\n\t"\
				"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t"\
				"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t"\
			/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */\
				"movq	%[__fq0],%%rdx		\n\t"\
				/* Digit 0: */\
				"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t"\
				"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t"\
				"mulpd	    (%%rdx),%%xmm0		\n\t	mulpd	0x10(%%rdx),%%xmm8		\n\t"\
				"addpd	%%xmm7,%%xmm0			\n\t	addpd	%%xmm7 ,%%xmm8			\n\t"\
				"subpd	%%xmm7,%%xmm0			\n\t	subpd	%%xmm7 ,%%xmm8			\n\t"\
				"mulpd	    %%xmm15,%%xmm0		\n\t	mulpd	    %%xmm15,%%xmm8		\n\t"\
				/* Digit 1: */\
				"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
				"mulpd	0x20(%%rdx),%%xmm3		\n\t	mulpd	0x30(%%rdx),%%xmm11		\n\t"\
				"addpd	%%xmm0,%%xmm3			\n\t	addpd	%%xmm8 ,%%xmm11			\n\t"\
				"mulpd	    (%%rdx),%%xmm1		\n\t	mulpd	0x10(%%rdx),%%xmm9		\n\t"\
				"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
				"addpd	%%xmm7,%%xmm1			\n\t	addpd	%%xmm7 ,%%xmm9			\n\t"\
				"subpd	%%xmm7,%%xmm1			\n\t	subpd	%%xmm7 ,%%xmm9			\n\t"\
				"mulpd	    %%xmm15,%%xmm1		\n\t	mulpd	    %%xmm15,%%xmm9		\n\t"\
				/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\
				"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t"\
				"movaps	%%xmm6,%%xmm0			\n\t	movaps	%%xmm14,%%xmm8			\n\t"\
				"mulpd	    (%%rdx),%%xmm2		\n\t	mulpd	0x10(%%rdx),%%xmm10		\n\t"\
				"mulpd	0x20(%%rdx),%%xmm6		\n\t	mulpd	0x30(%%rdx),%%xmm14		\n\t"\
				"mulpd	0x40(%%rdx),%%xmm4		\n\t	mulpd	0x50(%%rdx),%%xmm12		\n\t"\
				"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
				"addpd	%%xmm4,%%xmm1			\n\t	addpd	%%xmm12,%%xmm9			\n\t"\
				"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm9 ,%%xmm10			\n\t"\
				"subpd	0x20(%%rbx),%%xmm2		\n\t	subpd	0x20(%%rbx),%%xmm10		\n\t"\
				"addpd	%%xmm7,%%xmm2			\n\t	addpd	%%xmm7 ,%%xmm10			\n\t"\
				"subpd	%%xmm7,%%xmm2			\n\t	subpd	%%xmm7 ,%%xmm10			\n\t"\
				"mulpd	    %%xmm15,%%xmm2		\n\t	mulpd	    %%xmm15,%%xmm10		\n\t"\
				/* Precompute all the needed partial products: */\
				"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t"\
				"mulpd	0x40(%%rdx),%%xmm0		\n\t	mulpd	0x50(%%rdx),%%xmm8		\n\t"\
				"mulpd	0x20(%%rdx),%%xmm1		\n\t	mulpd	0x30(%%rdx),%%xmm9		\n\t"\
				"mulpd	0x40(%%rdx),%%xmm3		\n\t	mulpd	0x50(%%rdx),%%xmm11		\n\t"\
				/* Digit 3: */\
				"addpd	%%xmm2,%%xmm1			\n\t	addpd	%%xmm10,%%xmm9			\n\t"\
				"addpd	%%xmm1,%%xmm0			\n\t	addpd	%%xmm9 ,%%xmm8			\n\t"\
				"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
				"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
				"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
				"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t"\
				"mulpd	    %%xmm15,%%xmm6		\n\t	mulpd	    %%xmm15,%%xmm14		\n\t"\
				/* Digit 4: */\
				"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
				"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t"\
				"addpd	%%xmm7,%%xmm3			\n\t	addpd	%%xmm7 ,%%xmm11			\n\t"\
				"subpd	%%xmm7,%%xmm3			\n\t	subpd	%%xmm7 ,%%xmm11			\n\t"\
				"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t"\
				"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t"\
				/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\
				"movaps %%xmm5,%%xmm13	/* Need a copy of -0x10(%%rbx), a.k.a. %%xmm5 */\n\t"\
				"movaps	0xe0(%%rax),%%xmm2		\n\t	movaps	0xf0(%%rax),%%xmm10		\n\t"\
				"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t"\
				"movaps	0x60(%%rdx),%%xmm4		\n\t	movaps	0x70(%%rdx),%%xmm12		\n\t"\
				"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t"\
				"mulpd	    (%%rdx),%%xmm5		\n\t	mulpd	0x10(%%rdx),%%xmm13		\n\t"\
				"addpd	%%xmm1,%%xmm3			\n\t	addpd	%%xmm9 ,%%xmm11			\n\t"\
				"movaps	0xc0(%%rax),%%xmm1		\n\t	movaps	0xd0(%%rax),%%xmm9		\n\t"\
				"cmppd	$0x1,%%xmm3,%%xmm6		\n\t	cmppd	$0x1,%%xmm11,%%xmm14	\n\t"\
				"subpd	%%xmm0,%%xmm1			\n\t	subpd	%%xmm8 ,%%xmm9			\n\t"\
				"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t"\
				"movaps	%%xmm4,%%xmm0			\n\t	movaps	%%xmm12,%%xmm8			\n\t"\
				"andpd	%%xmm6,%%xmm4			\n\t	andpd	%%xmm14,%%xmm12			\n\t"\
				"addpd	%%xmm4,%%xmm2			\n\t	addpd	%%xmm12,%%xmm10			\n\t"\
				"andpd	%%xmm5,%%xmm6			\n\t	andpd	%%xmm13,%%xmm14			\n\t"\
				"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t"\
				/* qlo26 in xmm5, qhi52 in xmm0 */\
				/* if((pshift >> j) & (uint64)1) { */\
				"movq	%[__pshift],%%rax		\n\t"\
				"shrq	%%cl,%%rax				\n\t"\
				"andq	$0x1,%%rax				\n\t"\
				/* Branchless version of the conditional doubling: */\
				"xorq	$0x1,%%rax				\n\t"\
				"movd	%%rax,%%xmm6			\n\t"\
				"pshufd	$0,%%xmm6,%%xmm6		\n\t"\
				"cvtdq2pd	%%xmm6,%%xmm6		\n\t"\
				"xorpd	%%xmm3,%%xmm3			\n\t"\
				"cmppd	$0x0,%%xmm3,%%xmm6		\n\t"\
				/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm10,xmm9; qhi52,qlo26 in xmm8,xmm13 */\
				"movaps	%%xmm2,%%xmm4			\n\t	movaps	%%xmm10,%%xmm12			\n\t"\
				"movaps	%%xmm1,%%xmm3			\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
				"andpd	%%xmm6,%%xmm4			\n\t	andpd	%%xmm6 ,%%xmm12			\n\t"\
				"andpd	%%xmm6,%%xmm3			\n\t	andpd	%%xmm6 ,%%xmm11			\n\t"\
				"addpd	%%xmm4,%%xmm2			\n\t	addpd	%%xmm12,%%xmm10			\n\t"\
				"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
				/* If x > q, subtract q: */\
				"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
				"cmppd	$0x2,%%xmm2,%%xmm6		\n\t	cmppd	$0x2,%%xmm10,%%xmm14	\n\t"\
				"andpd	%%xmm6,%%xmm0			\n\t	andpd	%%xmm14,%%xmm8			\n\t"\
				"andpd	%%xmm6,%%xmm5			\n\t	andpd	%%xmm14,%%xmm13			\n\t"\
				"subpd	%%xmm0,%%xmm2			\n\t	subpd	%%xmm8 ,%%xmm10			\n\t"\
				"subpd	%%xmm5,%%xmm1			\n\t	subpd	%%xmm13,%%xmm9			\n\t"\
				/* } */\
				/* Normalize the result: */\
				"movaps	-0x10(%%rbx),%%xmm4		\n\t	movaps	-0x10(%%rbx),%%xmm12	\n\t"\
				"movaps	     %%xmm15,%%xmm3		\n\t	movaps	     %%xmm15,%%xmm11	\n\t"\
				"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm11,%%xmm9			\n\t"\
				"mulpd	%%xmm3,%%xmm2			\n\t	mulpd	%%xmm11,%%xmm10			\n\t"\
				"movaps	%%xmm1,%%xmm0			\n\t	movaps	%%xmm9 ,%%xmm8			\n\t"\
				"addpd	%%xmm7,%%xmm1			\n\t	addpd	%%xmm7 ,%%xmm9			\n\t"\
				"subpd	%%xmm7,%%xmm1			\n\t	subpd	%%xmm7 ,%%xmm9			\n\t"\
				"subpd	%%xmm1,%%xmm0			\n\t	subpd	%%xmm9 ,%%xmm8			\n\t"\
				"mulpd	%%xmm4,%%xmm0			\n\t	mulpd	%%xmm12,%%xmm8			\n\t"\
				"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm11,%%xmm9			\n\t"\
				"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm9 ,%%xmm10			\n\t"\
				"movaps	%%xmm2,%%xmm1			\n\t	movaps	%%xmm10,%%xmm9			\n\t"\
				"addpd	%%xmm7,%%xmm2			\n\t	addpd	%%xmm7 ,%%xmm10			\n\t"\
				"subpd	%%xmm7,%%xmm2			\n\t	subpd	%%xmm7 ,%%xmm10			\n\t"\
				"subpd	%%xmm2,%%xmm1			\n\t	subpd	%%xmm10,%%xmm9			\n\t"\
				"mulpd	%%xmm4,%%xmm1			\n\t	mulpd	%%xmm12,%%xmm9			\n\t"\
				"/* Move high 2 words of result into input registers expected by start of loop body: */\n\t"\
				"movaps	%%xmm2,%%xmm4 /* fx2 */	\n\t	movaps	%%xmm10,%%xmm12	/* hx2 */	\n\t"\
				"movaps	%%xmm1,%%xmm2 /* fx1 */	\n\t	movaps	%%xmm9 ,%%xmm10	/* hx1 */	\n\t"\
				"subq	$1,%%rcx	/* j-- */		\n\t"\
				"cmpq	$0,%%rcx	/* j > 0 ?	*/	\n\t"\
				"jge	Loop4Beg	/* if (j >= 0), loop */	\n\t"\
			"Loop4End:							\n\t"\
				:					/* outputs: none */\
				: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
				 ,[__two26i] "m" (two26i)	\
				 ,[__fx0]	 "m" (fx0)		\
				 ,[__fq0]	 "m" (fq0)		\
				 ,[__pshift] "m" (pshift)	\
				 ,[__start_index] "m" (start_index)	\
				: "cl","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
			);

		  #else	/* !ASM_LOOP: High-level loop construct, and explicit branch to do the modular doubling: */

			for(j = start_index-2; j >= 0; j--)
			{
				SSE2_twopmodq78_modmul_q4(fq0,pshift,j);
			}	/* for(j...) */

		  #endif	/* ASM_LOOP */

			__asm__ volatile (\
				"movq	%[__fx0],%%rax	\n\t"\
				"movaps	%%xmm0,    (%%rax)	\n\t	movaps	%%xmm8 ,0x10(%%rax)	\n\t"\
				"movaps	%%xmm2,0x20(%%rax)	\n\t	movaps	%%xmm10,0x30(%%rax)	\n\t"\
				"movaps	%%xmm4,0x40(%%rax)	\n\t	movaps	%%xmm12,0x50(%%rax)	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
			);

		#endif	/* OS_BITS = 32 / 64 ? */

	#endif	/* USE_SSE2 */

		/*...Double and return.	These are specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		*/
	#ifdef USE_SSE2
		CVT78_3WORD_DOUBLE_UINT96(*fx0,*fx1,*fx2, x0);
		CVT78_3WORD_DOUBLE_UINT96(*gx0,*gx1,*gx2, x1);
		CVT78_3WORD_DOUBLE_UINT96(*hx0,*hx1,*hx2, x2);
		CVT78_3WORD_DOUBLE_UINT96(*ix0,*ix1,*ix2, x3);
	#else
		CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
		CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
		CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
		CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);
	#endif
		ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
		ADD96(x1,x1,x1);
		ADD96(x2,x2,x2);
		ADD96(x3,x3,x3);
		q0.d0 -= FERMAT;
		q1.d0 -= FERMAT;
		q2.d0 -= FERMAT;
		q3.d0 -= FERMAT;
		SUB96(x0,q0,x0);
		SUB96(x1,q1,x1);
		SUB96(x2,q2,x2);
		SUB96(x3,q3,x3);

		tmp0 = CMPEQ96(x0, ONE96);
		tmp1 = CMPEQ96(x1, ONE96);
		tmp2 = CMPEQ96(x2, ONE96);
		tmp3 = CMPEQ96(x3, ONE96);
		r = tmp0;
		r += tmp1 << 1;
		r += tmp2 << 2;
		r += tmp3 << 3;
		return r;
	}

	/*** Older "reference" version ***/
	uint64 twopmodq78_3WORD_DOUBLE_q4_REF(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
	{
		 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
		uint64 lead7, tmp0, tmp1, tmp2, tmp3, r;
		uint96 q0, q1, q2, q3
			, qinv0, qhalf0, x0, lo0
			, qinv1, qhalf1, x1, lo1
			, qinv2, qhalf2, x2, lo2
			, qinv3, qhalf3, x3, lo3;
		uint192 prod192;
		static uint64 psave = 0, pshift;
		static uint32 start_index, zshift, first_entry = TRUE;
		double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2, fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2;
		double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2;
		double hq0,hq1,hq2, hqinv0,hqinv1,hqinv2, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hhi1,hhi2;
		double iq0,iq1,iq2, iqinv0,iqinv1,iqinv2, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ihi1,ihi2;
		uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

		ASSERT((p >> 63) == 0, "twopmodq78_q4 : p must be < 2^63!");
		q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;
	#ifdef MUL_LOHI64_SUBROUTINE
		// MUL_LOHI64 expects a 64-bit high-part pointer, in 32bit builds this buggers us if we try dumping hi-part directly into 32-bit q.d1
		MUL_LOHI64(q0.d0, k0,&q0.d0,&tmp0);	q0.d1 = tmp0;
		MUL_LOHI64(q1.d0, k1,&q1.d0,&tmp0);	q1.d1 = tmp0;
		MUL_LOHI64(q2.d0, k2,&q2.d0,&tmp0);	q2.d1 = tmp0;
		MUL_LOHI64(q3.d0, k3,&q3.d0,&tmp0);	q3.d1 = tmp0;
	#else
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
		MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
		MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);
	#endif
		q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
		q1.d0 += 1;
		q2.d0 += 1;
		q3.d0 += 1;
		ASSERT((q0.d1 >> 14) == 0, "twopmodq78_q4 : (q0.d1 >> 14) != 0");
		ASSERT((q1.d1 >> 14) == 0, "twopmodq78_q4 : (q1.d1 >> 14) != 0");
		ASSERT((q2.d1 >> 14) == 0, "twopmodq78_q4 : (q2.d1 >> 14) != 0");
		ASSERT((q3.d1 >> 14) == 0, "twopmodq78_q4 : (q3.d1 >> 14) != 0");

		/* Convert q to floating form: */
		CVT_UINT78_3WORD_DOUBLE(q0, fq0,fq1,fq2);
		CVT_UINT78_3WORD_DOUBLE(q1, gq0,gq1,gq2);
		CVT_UINT78_3WORD_DOUBLE(q2, hq0,hq1,hq2);
		CVT_UINT78_3WORD_DOUBLE(q3, iq0,iq1,iq2);

		RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
		RSHIFT_FAST96(q1, 1, qhalf1);
		RSHIFT_FAST96(q2, 1, qhalf2);
		RSHIFT_FAST96(q3, 1, qhalf3);

		if(first_entry || p != psave)
		{
			first_entry = FALSE;
			psave  = p;
			pshift = p + 78;
			j = leadz64(pshift);
			lead7 = ((pshift<<j) >> 57);
			if(lead7 > 77)
			{
				lead7 >>= 1;
				start_index =  64-j-6;	/* Use only the leftmost 6 bits */
			}
			else
				start_index =  64-j-7;

			zshift = 77 - lead7;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

			pshift = ~pshift;
		}

		qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
		qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;
		qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d1 = (uint64)0;
		qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d1 = (uint64)0;

		for(j = 0; j < 4; j++)
		{
			tmp0 = q0.d0*qinv0.d0;
			tmp1 = q1.d0*qinv1.d0;
			tmp2 = q2.d0*qinv2.d0;
			tmp3 = q3.d0*qinv3.d0;

			qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
			qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
			qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);
			qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);
		}

	#ifdef MUL_LOHI64_SUBROUTINE
		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
		qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + __MULH64(q2.d0, qinv2.d0));
		qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + __MULH64(q3.d0, qinv3.d0));
	#else
		MULH64(q0.d0, qinv0.d0, tmp0);
		MULH64(q1.d0, qinv1.d0, tmp1);
		MULH64(q2.d0, qinv2.d0, tmp2);
		MULH64(q3.d0, qinv3.d0, tmp3);

		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);
		qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + tmp2);
		qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + tmp3);
	#endif
		qinv0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
		qinv1.d1 &= 0x0000000000003fff;
		qinv2.d1 &= 0x0000000000003fff;
		qinv3.d1 &= 0x0000000000003fff;

		/* Convert qinv to floating form: */
		CVT_UINT78_3WORD_DOUBLE(qinv0, fqinv0,fqinv1,fqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv1, gqinv0,gqinv1,gqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv2, hqinv0,hqinv1,hqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv3, iqinv0,iqinv1,iqinv2);

		/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
		j = start_index-1;

		/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
		LSHIFT96(qinv0, zshift, lo0);	lo0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
		LSHIFT96(qinv1, zshift, lo1);	lo1.d1 &= 0x0000000000003fff;
		LSHIFT96(qinv2, zshift, lo2);	lo2.d1 &= 0x0000000000003fff;
		LSHIFT96(qinv3, zshift, lo3);	lo3.d1 &= 0x0000000000003fff;

		MUL_LOHI96_PROD192(q0,lo0,prod192);	RSHIFT192(prod192,78,prod192);	lo0.d0 = prod192.d0;	lo0.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q1,lo1,prod192);	RSHIFT192(prod192,78,prod192);	lo1.d0 = prod192.d0;	lo1.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q2,lo2,prod192);	RSHIFT192(prod192,78,prod192);	lo2.d0 = prod192.d0;	lo2.d1 = prod192.d1;
		MUL_LOHI96_PROD192(q3,lo3,prod192);	RSHIFT192(prod192,78,prod192);	lo3.d0 = prod192.d0;	lo3.d1 = prod192.d1;

		/* hi = 0 in this instance, which simplifies things. */
		SUB96(q0, lo0, x0);	/* Put the result in lo (rather than x), to ease overflow check below */
		SUB96(q1, lo1, x1);
		SUB96(q2, lo2, x2);
		SUB96(q3, lo3, x3);

		if((pshift >> j) & (uint64)1)
		{
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
			if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
			if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
			if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
		}

		/* Convert x to floating form: */
		CVT_UINT78_3WORD_DOUBLE(x0, fx0,fx1,fx2);
		CVT_UINT78_3WORD_DOUBLE(x1, gx0,gx1,gx2);
		CVT_UINT78_3WORD_DOUBLE(x2, hx0,hx1,hx2);
		CVT_UINT78_3WORD_DOUBLE(x3, ix0,ix1,ix2);

		for(j = start_index-2; j >= 0; j--)
		{
			/*...x^2 mod q is returned in x. */
			/* SQR_LOHI96(x,lo,hi); */
			SQR_LOHI78_3WORD_DOUBLE_q4(
				  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2
				, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2
				, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hhi1,hhi2
				, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ihi1,ihi2
			);
			/* MULL96(lo,qinv,lo); */
			MULL78_3WORD_DOUBLE_q4(
				  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
				, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
				, hlo0,hlo1,hlo2, hqinv0,hqinv1,hqinv2, hlo0,hlo1,hlo2
				, ilo0,ilo1,ilo2, iqinv0,iqinv1,iqinv2, ilo0,ilo1,ilo2
			);
			/* MULH96(q,lo,lo); */
			MULH78_3WORD_DOUBLE_q4(
				  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2
				, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2
				, hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlo1,hlo2
				, iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilo1,ilo2
			);
			/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
			if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
			{
				SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
				ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
			}
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

			if(CMPLT78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2))
			{
				SUB78_3WORD_DOUBLE(gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2);
				ADD78_3WORD_DOUBLE(glo0,glo1,glo2, ghi0,ghi1,ghi2, gx0,gx1,gx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2, gx0,gx1,gx2);
			}
			NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

			if(CMPLT78_3WORD_DOUBLE(hhi0,hhi1,hhi2, hlo0,hlo1,hlo2))
			{
				SUB78_3WORD_DOUBLE(hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlo1,hlo2);
				ADD78_3WORD_DOUBLE(hlo0,hlo1,hlo2, hhi0,hhi1,hhi2, hx0,hx1,hx2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(hhi0,hhi1,hhi2, hlo0,hlo1,hlo2, hx0,hx1,hx2);
			}
			NORMALIZE78_3WORD_DOUBLE(hx0,hx1,hx2);

			if(CMPLT78_3WORD_DOUBLE(ihi0,ihi1,ihi2, ilo0,ilo1,ilo2))
			{
				SUB78_3WORD_DOUBLE(iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilo1,ilo2);
				ADD78_3WORD_DOUBLE(ilo0,ilo1,ilo2, ihi0,ihi1,ihi2, ix0,ix1,ix2);
			}
			else
			{
				SUB78_3WORD_DOUBLE(ihi0,ihi1,ihi2, ilo0,ilo1,ilo2, ix0,ix1,ix2);
			}
			NORMALIZE78_3WORD_DOUBLE(ix0,ix1,ix2);

			if((pshift >> j) & (uint64)1)
			{
				/* Convert fx to uint96 for purpose of this particular comparison: */
				CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
				CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
				CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
				CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x0, qhalf0))
				{
					ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
					SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				}
				NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x1, qhalf1))
				{
					ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
					SUB78_3WORD_DOUBLE(gx0,gx1,gx2, gq0,gq1,gq2, gx0,gx1,gx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
				}
				NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x2, qhalf2))
				{
					ADD78_3WORD_DOUBLE(hx0,hx1,hx2, hx0,hx1,hx2, hx0,hx1,hx2);
					SUB78_3WORD_DOUBLE(hx0,hx1,hx2, hq0,hq1,hq2, hx0,hx1,hx2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(hx0,hx1,hx2, hx0,hx1,hx2, hx0,hx1,hx2);
				}
				NORMALIZE78_3WORD_DOUBLE(hx0,hx1,hx2);

				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x3, qhalf3))
				{
					ADD78_3WORD_DOUBLE(ix0,ix1,ix2, ix0,ix1,ix2, ix0,ix1,ix2);
					SUB78_3WORD_DOUBLE(ix0,ix1,ix2, iq0,iq1,iq2, ix0,ix1,ix2);
				}
				else
				{
					ADD78_3WORD_DOUBLE(ix0,ix1,ix2, ix0,ix1,ix2, ix0,ix1,ix2);
				}
				NORMALIZE78_3WORD_DOUBLE(ix0,ix1,ix2);

			}
		}

		/*...Double and return.	These are specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		*/
		CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
		CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
		CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
		CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);

		ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
		ADD96(x1,x1,x1);
		ADD96(x2,x2,x2);
		ADD96(x3,x3,x3);
		q0.d0 -= FERMAT;
		q1.d0 -= FERMAT;
		q2.d0 -= FERMAT;
		q3.d0 -= FERMAT;
		SUB96(x0,q0,x0);
		SUB96(x1,q1,x1);
		SUB96(x2,q2,x2);
		SUB96(x3,q3,x3);

		r = 0;
		if(CMPEQ96(x0, ONE96)) r +=  1;
		if(CMPEQ96(x1, ONE96)) r +=  2;
		if(CMPEQ96(x2, ONE96)) r +=  4;
		if(CMPEQ96(x3, ONE96)) r +=  8;
		return(r);
	}

	#ifdef USE_ARM_V8_SIMD

		// No TF support on ARMv8
		uint64 twopmodq78_3WORD_DOUBLE_q8(uint64 p, uint64 k[], int init_sse2, int thr_id)
		{
			ASSERT(0,"No TF support on ARMv8!");
		}

	#elif defined(X64_ASM) && defined(USE_SSE2)

	  #if 1	// Nov 2013: Now that have prototyped 2-TF-input/4-SSE-register modmul macro, use 8-input all-float in both SSE2 and AVX mode
			// Set = 0 here to enable hybrid (SSE2-mode only)

		// AVX/8-input analog of twopmodq78_3WORD_DOUBLE_q4 :
		uint64 twopmodq78_3WORD_DOUBLE_q8(uint64 p, uint64 k[], int init_sse2, int thr_id)
		{
			const char func[] = "twopmodq78_3WORD_DOUBLE_q8 [AVX]";
			 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
			uint32 q32_0, qinv32_0, tmp32_0
				 , q32_1, qinv32_1, tmp32_1
				 , q32_2, qinv32_2, tmp32_2
				 , q32_3, qinv32_3, tmp32_3
				 , q32_4, qinv32_4, tmp32_4
				 , q32_5, qinv32_5, tmp32_5
				 , q32_6, qinv32_6, tmp32_6
				 , q32_7, qinv32_7, tmp32_7;
			uint64 lead7, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r;
			uint96 q0, qinv0, qhalf0, x0, lo0
				 , q1, qinv1, qhalf1, x1, lo1
				 , q2, qinv2, qhalf2, x2, lo2
				 , q3, qinv3, qhalf3, x3, lo3
				 , q4, qinv4, qhalf4, x4, lo4
				 , q5, qinv5, qhalf5, x5, lo5
				 , q6, qinv6, qhalf6, x6, lo6
				 , q7, qinv7, qhalf7, x7, lo7;
			uint192 prod192;
			static uint64 psave = 0, pshift;
			static uint32 start_index, zshift, first_entry = TRUE;
			uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

			static int max_threads = 1;	// Default local-array-init is for just a single thread ... caller can re-init for > 1 threads later, if desired.
			static double *sc_arr = 0x0, *sc_ptr;
		  #ifdef MULTITHREAD
			static double *__r0;	// Base address for discrete per-thread local stores ...
									// *NOTE*: more convenient to use double* rather than vec_dbl* here
		  #else
			static	// Following set of pointers only static-izable in single-thread mode
		  #endif
			double *aq0,*aq1,*aq2,*aqhi52, *aqinv0,*aqinv1,*aqinv2, *ax0,*ax1,*ax2
				, *bq0,*bq1,*bq2,*bqhi52, *bqinv0,*bqinv1,*bqinv2, *bx0,*bx1,*bx2
				, *cq0,*cq1,*cq2,*cqhi52, *cqinv0,*cqinv1,*cqinv2, *cx0,*cx1,*cx2
				, *dq0,*dq1,*dq2,*dqhi52, *dqinv0,*dqinv1,*dqinv2, *dx0,*dx1,*dx2
				, *eq0,*eq1,*eq2,*eqhi52, *eqinv0,*eqinv1,*eqinv2, *ex0,*ex1,*ex2
				, *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2
				, *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2
				, *hq0,*hq1,*hq2,*hqhi52, *hqinv0,*hqinv1,*hqinv2, *hx0,*hx1,*hx2
				, *two26f, *two26i, *two13i;

			if(p != psave)
			{
			//	first_entry = FALSE;
				psave  = p;
				pshift = p + 78;
				j = leadz64(pshift);
				/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 77: */
				lead7 = ((pshift<<j) >> 57);
				if(lead7 > 77)
				{
					lead7 >>= 1;
					start_index =  64-j-6;	/* Use only the leftmost 6 bits */
				}
				else
					start_index =  64-j-7;

				zshift = 77 - lead7;
				zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

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
				// Alloc the local-memory block - SSE2 needs 6 fewer double-slots than AVX (since only need one copy each of two13i,two26f,two26i), but use same AVX-alloc for both:
				sc_arr = ALLOC_DOUBLE(sc_arr, 0x6c*max_threads + 4);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
				sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);	// Force vec_dbl-alignment
				ASSERT(((uintptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			#ifdef MULTITHREAD
				__r0  = sc_ptr;
			  #ifdef USE_AVX
				two13i = sc_ptr + 0x60;
				two26f = sc_ptr + 0x64;
				two26i = sc_ptr + 0x68;	// Equivalent of 27 AVX slots, alloc 28
			  #else	// SSE2:
				two13i = sc_ptr + 0x60;
				two26f = sc_ptr + 0x62;
				two26i = sc_ptr + 0x64;	// Equivalent of 51 SSE2 slots, alloc 56 = 28*2, same as AVX
			  #endif
				for(j = 0; j < max_threads; ++j) {
					// Need specific-length inits here since might be building overall in AVX-512 mode:
				  #ifdef USE_AVX
					VEC_DBL_INIT_4((vec_dbl*)two13i  , TWO13FLINV);
					VEC_DBL_INIT_4((vec_dbl*)two26f  , TWO26FLOAT);
					VEC_DBL_INIT_4((vec_dbl*)two26i  , TWO26FLINV);
				  #else	// SSE2:
					VEC_DBL_INIT_2((vec_dbl*)two13i  , TWO13FLINV);
					VEC_DBL_INIT_2((vec_dbl*)two26f  , TWO26FLOAT);
					VEC_DBL_INIT_2((vec_dbl*)two26i  , TWO26FLINV);
				  #endif
					// Move on to next thread's local store:
					two13i   += 0x6c;
					two26f   += 0x6c;
					two26i   += 0x6c;
				}
			#else
				/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 4 to span an AVX register: */
				aq0    = sc_ptr + 0x00;	bq0    = sc_ptr + 0x01;	cq0    = sc_ptr + 0x02;	dq0    = sc_ptr + 0x03;	eq0    = sc_ptr + 0x04;	fq0    = sc_ptr + 0x05;	gq0    = sc_ptr + 0x06;	hq0    = sc_ptr + 0x07;
				aq1    = sc_ptr + 0x08;	bq1    = sc_ptr + 0x09;	cq1    = sc_ptr + 0x0a;	dq1    = sc_ptr + 0x0b;	eq1    = sc_ptr + 0x0c;	fq1    = sc_ptr + 0x0d;	gq1    = sc_ptr + 0x0e;	hq1    = sc_ptr + 0x0f;
				aq2    = sc_ptr + 0x10;	bq2    = sc_ptr + 0x11;	cq2    = sc_ptr + 0x12;	dq2    = sc_ptr + 0x13;	eq2    = sc_ptr + 0x14;	fq2    = sc_ptr + 0x15;	gq2    = sc_ptr + 0x16;	hq2    = sc_ptr + 0x17;
				aqhi52 = sc_ptr + 0x18;	bqhi52 = sc_ptr + 0x19;	cqhi52 = sc_ptr + 0x1a;	dqhi52 = sc_ptr + 0x1b;	eqhi52 = sc_ptr + 0x1c;	fqhi52 = sc_ptr + 0x1d;	gqhi52 = sc_ptr + 0x1e;	hqhi52 = sc_ptr + 0x1f;
				aqinv0 = sc_ptr + 0x20;	bqinv0 = sc_ptr + 0x21;	cqinv0 = sc_ptr + 0x22;	dqinv0 = sc_ptr + 0x23;	eqinv0 = sc_ptr + 0x24;	fqinv0 = sc_ptr + 0x25;	gqinv0 = sc_ptr + 0x26;	hqinv0 = sc_ptr + 0x27;
				aqinv1 = sc_ptr + 0x28;	bqinv1 = sc_ptr + 0x29;	cqinv1 = sc_ptr + 0x2a;	dqinv1 = sc_ptr + 0x2b;	eqinv1 = sc_ptr + 0x2c;	fqinv1 = sc_ptr + 0x2d;	gqinv1 = sc_ptr + 0x2e;	hqinv1 = sc_ptr + 0x2f;
				aqinv2 = sc_ptr + 0x30;	bqinv2 = sc_ptr + 0x31;	cqinv2 = sc_ptr + 0x32;	dqinv2 = sc_ptr + 0x33;	eqinv2 = sc_ptr + 0x34;	fqinv2 = sc_ptr + 0x35;	gqinv2 = sc_ptr + 0x36;	hqinv2 = sc_ptr + 0x37;
				ax0    = sc_ptr + 0x38;	bx0    = sc_ptr + 0x39;	cx0    = sc_ptr + 0x3a;	dx0    = sc_ptr + 0x3b;	ex0    = sc_ptr + 0x3c;	fx0    = sc_ptr + 0x3d;	gx0    = sc_ptr + 0x3e;	hx0    = sc_ptr + 0x3f;
				ax1    = sc_ptr + 0x40;	bx1    = sc_ptr + 0x41;	cx1    = sc_ptr + 0x42;	dx1    = sc_ptr + 0x43;	ex1    = sc_ptr + 0x44;	fx1    = sc_ptr + 0x45;	gx1    = sc_ptr + 0x46;	hx1    = sc_ptr + 0x47;
				ax2    = sc_ptr + 0x48;	bx2    = sc_ptr + 0x49;	cx2    = sc_ptr + 0x4a;	dx2    = sc_ptr + 0x4b;	ex2    = sc_ptr + 0x4c;	fx2    = sc_ptr + 0x4d;	gx2    = sc_ptr + 0x4e;	hx2    = sc_ptr + 0x4f;
				// +0x50,58 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
			  #ifdef USE_AVX
				two13i = sc_ptr + 0x60;
				two26f = sc_ptr + 0x64;
				two26i = sc_ptr + 0x68;	// Equivalent of 27 AVX slots, alloc 28
			  #else	// SSE2:
				two13i = sc_ptr + 0x60;
				two26f = sc_ptr + 0x62;
				two26i = sc_ptr + 0x64;	// Equivalent of 51 SSE2 slots, alloc 56 = 28*2, same as AVX
			  #endif
				/* Can premultiply each of the multiword-mul inputs by 1/sqrt(2^13) due to quadraticity: */
				// Need specific-length inits here since might be building overall in AVX-512 mode:
			  #ifdef USE_AVX
				VEC_DBL_INIT_4((vec_dbl*)two13i  , TWO13FLINV);
				VEC_DBL_INIT_4((vec_dbl*)two26f  , TWO26FLOAT);
				VEC_DBL_INIT_4((vec_dbl*)two26i  , TWO26FLINV);
			  #else	// SSE2:
				VEC_DBL_INIT_2((vec_dbl*)two13i  , TWO13FLINV);
				VEC_DBL_INIT_2((vec_dbl*)two26f  , TWO26FLOAT);
				VEC_DBL_INIT_2((vec_dbl*)two26i  , TWO26FLINV);
			  #endif
			#endif
				if(init_sse2) return 0;
			}	/* end of inits */

			/* If multithreaded, set the local-store pointers needed for the current thread; */
		#ifdef MULTITHREAD
			ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
			sc_ptr = __r0 + thr_id*0x6c;
			/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 4 to span an AVX register: */
			aq0    = sc_ptr + 0x00;	bq0    = sc_ptr + 0x01;	cq0    = sc_ptr + 0x02;	dq0    = sc_ptr + 0x03;	eq0    = sc_ptr + 0x04;	fq0    = sc_ptr + 0x05;	gq0    = sc_ptr + 0x06;	hq0    = sc_ptr + 0x07;
			aq1    = sc_ptr + 0x08;	bq1    = sc_ptr + 0x09;	cq1    = sc_ptr + 0x0a;	dq1    = sc_ptr + 0x0b;	eq1    = sc_ptr + 0x0c;	fq1    = sc_ptr + 0x0d;	gq1    = sc_ptr + 0x0e;	hq1    = sc_ptr + 0x0f;
			aq2    = sc_ptr + 0x10;	bq2    = sc_ptr + 0x11;	cq2    = sc_ptr + 0x12;	dq2    = sc_ptr + 0x13;	eq2    = sc_ptr + 0x14;	fq2    = sc_ptr + 0x15;	gq2    = sc_ptr + 0x16;	hq2    = sc_ptr + 0x17;
			aqhi52 = sc_ptr + 0x18;	bqhi52 = sc_ptr + 0x19;	cqhi52 = sc_ptr + 0x1a;	dqhi52 = sc_ptr + 0x1b;	eqhi52 = sc_ptr + 0x1c;	fqhi52 = sc_ptr + 0x1d;	gqhi52 = sc_ptr + 0x1e;	hqhi52 = sc_ptr + 0x1f;
			aqinv0 = sc_ptr + 0x20;	bqinv0 = sc_ptr + 0x21;	cqinv0 = sc_ptr + 0x22;	dqinv0 = sc_ptr + 0x23;	eqinv0 = sc_ptr + 0x24;	fqinv0 = sc_ptr + 0x25;	gqinv0 = sc_ptr + 0x26;	hqinv0 = sc_ptr + 0x27;
			aqinv1 = sc_ptr + 0x28;	bqinv1 = sc_ptr + 0x29;	cqinv1 = sc_ptr + 0x2a;	dqinv1 = sc_ptr + 0x2b;	eqinv1 = sc_ptr + 0x2c;	fqinv1 = sc_ptr + 0x2d;	gqinv1 = sc_ptr + 0x2e;	hqinv1 = sc_ptr + 0x2f;
			aqinv2 = sc_ptr + 0x30;	bqinv2 = sc_ptr + 0x31;	cqinv2 = sc_ptr + 0x32;	dqinv2 = sc_ptr + 0x33;	eqinv2 = sc_ptr + 0x34;	fqinv2 = sc_ptr + 0x35;	gqinv2 = sc_ptr + 0x36;	hqinv2 = sc_ptr + 0x37;
			ax0    = sc_ptr + 0x38;	bx0    = sc_ptr + 0x39;	cx0    = sc_ptr + 0x3a;	dx0    = sc_ptr + 0x3b;	ex0    = sc_ptr + 0x3c;	fx0    = sc_ptr + 0x3d;	gx0    = sc_ptr + 0x3e;	hx0    = sc_ptr + 0x3f;
			ax1    = sc_ptr + 0x40;	bx1    = sc_ptr + 0x41;	cx1    = sc_ptr + 0x42;	dx1    = sc_ptr + 0x43;	ex1    = sc_ptr + 0x44;	fx1    = sc_ptr + 0x45;	gx1    = sc_ptr + 0x46;	hx1    = sc_ptr + 0x47;
			ax2    = sc_ptr + 0x48;	bx2    = sc_ptr + 0x49;	cx2    = sc_ptr + 0x4a;	dx2    = sc_ptr + 0x4b;	ex2    = sc_ptr + 0x4c;	fx2    = sc_ptr + 0x4d;	gx2    = sc_ptr + 0x4e;	hx2    = sc_ptr + 0x4f;
			// +0x50,58 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
		  #ifdef USE_AVX
			two13i = sc_ptr + 0x60;
			two26f = sc_ptr + 0x64;
			two26i = sc_ptr + 0x68;	// Equivalent of 27 AVX slots, alloc 28
		  #else	// SSE2:
			two13i = sc_ptr + 0x60;
			two26f = sc_ptr + 0x62;
			two26i = sc_ptr + 0x64;	// Equivalent of 51 SSE2 slots, alloc 56 = 28*2, same as AVX
		  #endif
		#endif

			ASSERT((p >> 63) == 0, "twopmodq78_q8 : p must be < 2^63!");
			q0.d0 = q1.d0 = q2.d0 = q3.d0 = q4.d0 = q5.d0 = q6.d0 = q7.d0 = p+p;
		#ifdef MUL_LOHI64_SUBROUTINE
			#error MUL_LOHI64_SUBROUTINE defined!
		#endif
			MUL_LOHI64(q0.d0, k[0], q0.d0, q0.d1);
			MUL_LOHI64(q1.d0, k[1], q1.d0, q1.d1);
			MUL_LOHI64(q2.d0, k[2], q2.d0, q2.d1);
			MUL_LOHI64(q3.d0, k[3], q3.d0, q3.d1);
			MUL_LOHI64(q4.d0, k[4], q4.d0, q4.d1);
			MUL_LOHI64(q5.d0, k[5], q5.d0, q5.d1);
			MUL_LOHI64(q6.d0, k[6], q6.d0, q6.d1);
			MUL_LOHI64(q7.d0, k[7], q7.d0, q7.d1);
			/* Since 2*p*k even, no need to check for overflow here */
			q0.d0 += 1;
			q1.d0 += 1;
			q2.d0 += 1;
			q3.d0 += 1;
			q4.d0 += 1;
			q5.d0 += 1;
			q6.d0 += 1;
			q7.d0 += 1;
			ASSERT((q0.d1 >> 14) == 0, "twopmodq78_q8 : (q0.d1 >> 14) != 0");
			ASSERT((q1.d1 >> 14) == 0, "twopmodq78_q8 : (q1.d1 >> 14) != 0");
			ASSERT((q2.d1 >> 14) == 0, "twopmodq78_q8 : (q2.d1 >> 14) != 0");
			ASSERT((q3.d1 >> 14) == 0, "twopmodq78_q8 : (q3.d1 >> 14) != 0");
			ASSERT((q4.d1 >> 14) == 0, "twopmodq78_q8 : (q4.d1 >> 14) != 0");
			ASSERT((q5.d1 >> 14) == 0, "twopmodq78_q8 : (q5.d1 >> 14) != 0");
			ASSERT((q6.d1 >> 14) == 0, "twopmodq78_q8 : (q6.d1 >> 14) != 0");
			ASSERT((q7.d1 >> 14) == 0, "twopmodq78_q8 : (q7.d1 >> 14) != 0");

			q32_0 = (uint32)q0.d0;
			q32_1 = (uint32)q1.d0;
			q32_2 = (uint32)q2.d0;
			q32_3 = (uint32)q3.d0;
			q32_4 = (uint32)q4.d0;
			q32_5 = (uint32)q5.d0;
			q32_6 = (uint32)q6.d0;
			q32_7 = (uint32)q7.d0;

			/* Convert q to floating form: */
			CVT_UINT78_3WORD_DOUBLE(q0 ,*aq0,*aq1,*aq2);
			CVT_UINT78_3WORD_DOUBLE(q1 ,*bq0,*bq1,*bq2);
			CVT_UINT78_3WORD_DOUBLE(q2 ,*cq0,*cq1,*cq2);
			CVT_UINT78_3WORD_DOUBLE(q3 ,*dq0,*dq1,*dq2);
			CVT_UINT78_3WORD_DOUBLE(q4 ,*eq0,*eq1,*eq2);
			CVT_UINT78_3WORD_DOUBLE(q5 ,*fq0,*fq1,*fq2);
			CVT_UINT78_3WORD_DOUBLE(q6 ,*gq0,*gq1,*gq2);
			CVT_UINT78_3WORD_DOUBLE(q7 ,*hq0,*hq1,*hq2);

		#ifdef USE_AVX

			__asm__ volatile (\
				"movq	%[__aq0],%%rax						\n\t"\
				"movq	%[__two26f],%%rsi					\n\t"\
				"vmovaps	    (%%rsi),%%ymm8	/* two26f */	\n\t"\
				"vmovaps	0x20(%%rsi),%%ymm9	/* two26i */	\n\t"\
				"vmovaps	    (%%rax),%%ymm0	/* aq0 */		\n\t	vmovaps	0x20(%%rax),%%ymm4	/* eq0 */		\n\t"\
				"vmovaps	0x40(%%rax),%%ymm1	/* aq1 */		\n\t	vmovaps	0x60(%%rax),%%ymm5	/* eq1 */		\n\t"\
				"vmovaps	0x80(%%rax),%%ymm2	/* aq2 */		\n\t	vmovaps	0xa0(%%rax),%%ymm6	/* eq2 */		\n\t"\
				"vmovaps	%%ymm2,%%ymm3		/* cpy aq2 */	\n\t	vmovaps	%%ymm6,%%ymm7		/* cpy eq2 */	\n\t"\
				"vmulpd	%%ymm8,%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm8,%%ymm7,%%ymm7				\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3/* Hi 52 out bits */\n\t	vaddpd	%%ymm5,%%ymm7,%%ymm7	/* Hi 52*/	\n\t"\
				"vmulpd	%%ymm9,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm9,%%ymm4,%%ymm4				\n\t"\
				"vmulpd	%%ymm9,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm9,%%ymm5,%%ymm5				\n\t"\
				"vmulpd	%%ymm9,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm9,%%ymm6,%%ymm6				\n\t"\
				"vmovaps	%%ymm0,    (%%rax)	/* aq0/2^26 */	\n\t	vmovaps	%%ymm4,0x20(%%rax)	/* eq0/2^26 */	\n\t"\
				"vmovaps	%%ymm1,0x40(%%rax)	/* aq1/2^26 */	\n\t	vmovaps	%%ymm5,0x60(%%rax)	/* eq1/2^26 */	\n\t"\
				"vmovaps	%%ymm2,0x80(%%rax)	/* aq2/2^26 */	\n\t	vmovaps	%%ymm6,0xa0(%%rax)	/* eq2/2^26 */	\n\t"\
				"vmovaps	%%ymm3,0xc0(%%rax)	/* aq[hi52] */	\n\t	vmovaps	%%ymm7,0xe0(%%rax)	/* eq[hi52] */	\n\t"\
				:					/* outputs: none */\
				: [__aq0] "m" (aq0)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
			);

		#else	// SSE2:

			__asm__ volatile (\
				"movq	%[__aq0],%%rax						\n\t"\
				"movq	%[__two26f],%%rsi					\n\t"\
				"movaps	    (%%rsi),%%xmm8	/* two26f */	\n\t"\
				"movaps	0x10(%%rsi),%%xmm9	/* two26i */	\n\t"\
				"movaps	    (%%rax),%%xmm0	/* aq0 */		\n\t	movaps	0x10(%%rax),%%xmm4	/* cq0 */		\n\t"\
				"movaps	0x40(%%rax),%%xmm1	/* aq1 */		\n\t	movaps	0x50(%%rax),%%xmm5	/* cq1 */		\n\t"\
				"movaps	0x80(%%rax),%%xmm2	/* aq2 */		\n\t	movaps	0x90(%%rax),%%xmm6	/* cq2 */		\n\t"\
				"movaps	%%xmm2,%%xmm3		/* cpy aq2 */	\n\t	movaps	%%xmm6,%%xmm7		/* cpy cq2 */	\n\t"\
				"mulpd	%%xmm8,%%xmm3						\n\t	mulpd	%%xmm8,%%xmm7						\n\t"\
				"addpd	%%xmm1,%%xmm3	/* Hi 52 out bits */\n\t	addpd	%%xmm5,%%xmm7		/* Hi 52 */		\n\t"\
				"mulpd	%%xmm9,%%xmm0						\n\t	mulpd	%%xmm9,%%xmm4						\n\t"\
				"mulpd	%%xmm9,%%xmm1						\n\t	mulpd	%%xmm9,%%xmm5						\n\t"\
				"mulpd	%%xmm9,%%xmm2						\n\t	mulpd	%%xmm9,%%xmm6						\n\t"\
				"movaps	%%xmm0,    (%%rax)	/* aq0/2^26 */	\n\t	movaps	%%xmm4,0x10(%%rax)	/* cq0/2^26 */	\n\t"\
				"movaps	%%xmm1,0x40(%%rax)	/* aq1/2^26 */	\n\t	movaps	%%xmm5,0x50(%%rax)	/* cq1/2^26 */	\n\t"\
				"movaps	%%xmm2,0x80(%%rax)	/* aq2/2^26 */	\n\t	movaps	%%xmm6,0x90(%%rax)	/* cq2/2^26 */	\n\t"\
				"movaps	%%xmm3,0xc0(%%rax)	/* aq[hi52] */	\n\t	movaps	%%xmm7,0xd0(%%rax)	/* cq[hi52] */	\n\t"\
				/* Now do 2nd set of doubles: */\
				"movaps	0x20(%%rax),%%xmm0	/* eq0 */		\n\t	movaps	0x30(%%rax),%%xmm4	/* gq0 */		\n\t"\
				"movaps	0x60(%%rax),%%xmm1	/* eq1 */		\n\t	movaps	0x70(%%rax),%%xmm5	/* gq1 */		\n\t"\
				"movaps	0xa0(%%rax),%%xmm2	/* eq2 */		\n\t	movaps	0xb0(%%rax),%%xmm6	/* gq2 */		\n\t"\
				"movaps	%%xmm2,%%xmm3		/* cpy eq2 */	\n\t	movaps	%%xmm6,%%xmm7		/* cpy gq2 */	\n\t"\
				"mulpd	%%xmm8,%%xmm3						\n\t	mulpd	%%xmm8,%%xmm7						\n\t"\
				"addpd	%%xmm1,%%xmm3	/* Hi 52 out bits */\n\t	addpd	%%xmm5,%%xmm7		/* Hi 52 */		\n\t"\
				"mulpd	%%xmm9,%%xmm0						\n\t	mulpd	%%xmm9,%%xmm4						\n\t"\
				"mulpd	%%xmm9,%%xmm1						\n\t	mulpd	%%xmm9,%%xmm5						\n\t"\
				"mulpd	%%xmm9,%%xmm2						\n\t	mulpd	%%xmm9,%%xmm6						\n\t"\
				"movaps	%%xmm0,0x20(%%rax)	/* eq0/2^26 */	\n\t	movaps	%%xmm4,0x30(%%rax)	/* gq0/2^26 */	\n\t"\
				"movaps	%%xmm1,0x60(%%rax)	/* eq1/2^26 */	\n\t	movaps	%%xmm5,0x70(%%rax)	/* gq1/2^26 */	\n\t"\
				"movaps	%%xmm2,0xa0(%%rax)	/* eq2/2^26 */	\n\t	movaps	%%xmm6,0xb0(%%rax)	/* gq2/2^26 */	\n\t"\
				"movaps	%%xmm3,0xe0(%%rax)	/* eq[hi52] */	\n\t	movaps	%%xmm7,0xf0(%%rax)	/* gq[hi52] */	\n\t"\
				:					/* outputs: none */\
				: [__aq0] "m" (aq0)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
			);

		#endif

			/* = (q-1)/2, since q odd. */
			RSHIFT_FAST96(q0, 1, qhalf0);
			RSHIFT_FAST96(q1, 1, qhalf1);
			RSHIFT_FAST96(q2, 1, qhalf2);
			RSHIFT_FAST96(q3, 1, qhalf3);
			RSHIFT_FAST96(q4, 1, qhalf4);
			RSHIFT_FAST96(q5, 1, qhalf5);
			RSHIFT_FAST96(q6, 1, qhalf6);
			RSHIFT_FAST96(q7, 1, qhalf7);

			qinv32_0 = (q32_0 + q32_0 + q32_0) ^ (uint32)2;					qinv32_4 = (q32_4 + q32_4 + q32_4) ^ (uint32)2;
			qinv32_1 = (q32_1 + q32_1 + q32_1) ^ (uint32)2;					qinv32_5 = (q32_5 + q32_5 + q32_5) ^ (uint32)2;
			qinv32_2 = (q32_2 + q32_2 + q32_2) ^ (uint32)2;					qinv32_6 = (q32_6 + q32_6 + q32_6) ^ (uint32)2;
			qinv32_3 = (q32_3 + q32_3 + q32_3) ^ (uint32)2;					qinv32_7 = (q32_7 + q32_7 + q32_7) ^ (uint32)2;
			/* 4 bits: */													/* 4 bits: */
			tmp32_0 = q32_0*qinv32_0;										tmp32_4 = q32_4*qinv32_4;
			tmp32_1 = q32_1*qinv32_1;										tmp32_5 = q32_5*qinv32_5;
			tmp32_2 = q32_2*qinv32_2;										tmp32_6 = q32_6*qinv32_6;
			tmp32_3 = q32_3*qinv32_3;										tmp32_7 = q32_7*qinv32_7;
			qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);						qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);
			qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);						qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);
			qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);						qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);
			qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);						qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);
			/* 8 bits: */													/* 8 bits: */
			tmp32_0 = q32_0*qinv32_0;										tmp32_4 = q32_4*qinv32_4;
			tmp32_1 = q32_1*qinv32_1;										tmp32_5 = q32_5*qinv32_5;
			tmp32_2 = q32_2*qinv32_2;										tmp32_6 = q32_6*qinv32_6;
			tmp32_3 = q32_3*qinv32_3;										tmp32_7 = q32_7*qinv32_7;
			qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);						qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);
			qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);						qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);
			qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);						qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);
			qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);						qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);
			/* 16 bits: */													/* 16 bits: */
			tmp32_0 = q32_0*qinv32_0;										tmp32_4 = q32_4*qinv32_4;
			tmp32_1 = q32_1*qinv32_1;										tmp32_5 = q32_5*qinv32_5;
			tmp32_2 = q32_2*qinv32_2;										tmp32_6 = q32_6*qinv32_6;
			tmp32_3 = q32_3*qinv32_3;										tmp32_7 = q32_7*qinv32_7;
			qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);						qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);
			qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);						qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);
			qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);						qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);
			qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);						qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);
			/* 32 bits: */													/* 32 bits: */
			qinv0.d0 = (uint64)qinv32_0;									qinv4.d0 = (uint64)qinv32_4;
			qinv1.d0 = (uint64)qinv32_1;									qinv5.d0 = (uint64)qinv32_5;
			qinv2.d0 = (uint64)qinv32_2;									qinv6.d0 = (uint64)qinv32_6;
			qinv3.d0 = (uint64)qinv32_3;									qinv7.d0 = (uint64)qinv32_7;
			tmp0 = q0.d0*qinv0.d0;											tmp4 = q4.d0*qinv4.d0;
			tmp1 = q1.d0*qinv1.d0;											tmp5 = q5.d0*qinv5.d0;
			tmp2 = q2.d0*qinv2.d0;											tmp6 = q6.d0*qinv6.d0;
			tmp3 = q3.d0*qinv3.d0;											tmp7 = q7.d0*qinv7.d0;
			qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);							qinv4.d0 = qinv4.d0*((uint64)2 - tmp4);
			qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);							qinv5.d0 = qinv5.d0*((uint64)2 - tmp5);
			qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);							qinv6.d0 = qinv6.d0*((uint64)2 - tmp6);
			qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);							qinv7.d0 = qinv7.d0*((uint64)2 - tmp7);
			/* 64 bits: */													/* 64 bits: */
			MULH64(q0.d0, qinv0.d0, tmp0);									MULH64(q4.d0, qinv4.d0, tmp4);
			MULH64(q1.d0, qinv1.d0, tmp1);									MULH64(q5.d0, qinv5.d0, tmp5);
			MULH64(q2.d0, qinv2.d0, tmp2);									MULH64(q6.d0, qinv6.d0, tmp6);
			MULH64(q3.d0, qinv3.d0, tmp3);									MULH64(q7.d0, qinv7.d0, tmp7);
			qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);					qinv4.d1 = -qinv4.d0*(q4.d1*qinv4.d0 + tmp4);
			qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);					qinv5.d1 = -qinv5.d0*(q5.d1*qinv5.d0 + tmp5);
			qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + tmp2);					qinv6.d1 = -qinv6.d0*(q6.d1*qinv6.d0 + tmp6);
			qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + tmp3);					qinv7.d1 = -qinv7.d0*(q7.d1*qinv7.d0 + tmp7);
			/* 64 bits - Only want the lower 14 bits here: */				/* 64 bits - Only want the lower 14 bits here: */
			qinv0.d1 &= 0x3fff;												qinv4.d1 &= 0x3fff;
			qinv1.d1 &= 0x3fff;												qinv5.d1 &= 0x3fff;
			qinv2.d1 &= 0x3fff;												qinv6.d1 &= 0x3fff;
			qinv3.d1 &= 0x3fff;												qinv7.d1 &= 0x3fff;
			/* Convert qinv to floating form: */							/* Convert qinv to floating form: */
			CVT_UINT78_3WORD_DOUBLE(qinv0 ,*aqinv0,*aqinv1,*aqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv1 ,*bqinv0,*bqinv1,*bqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv2 ,*cqinv0,*cqinv1,*cqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv3 ,*dqinv0,*dqinv1,*dqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv4 ,*eqinv0,*eqinv1,*eqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv5 ,*fqinv0,*fqinv1,*fqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv6 ,*gqinv0,*gqinv1,*gqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv7 ,*hqinv0,*hqinv1,*hqinv2);

		#ifdef USE_AVX

			__asm__ volatile (\
				"movq	%[__aqinv0],%%rax	\n\t"\
				"movq	%[__two26i],%%rsi	\n\t"\
				"vmovaps	(%%rsi),%%ymm6		\n\t"\
				"vmovaps	    (%%rax),%%ymm0	\n\t	vmovaps	0x20(%%rax),%%ymm3	\n\t"\
				"vmovaps	0x40(%%rax),%%ymm1	\n\t	vmovaps	0x60(%%rax),%%ymm4	\n\t"\
				"vmovaps	0x80(%%rax),%%ymm2	\n\t	vmovaps	0xa0(%%rax),%%ymm5	\n\t"\
				"vmulpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm6,%%ymm3,%%ymm3	\n\t"\
				"vmulpd	%%ymm6,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
				"vmovaps	%%ymm0,    (%%rax)	\n\t	vmovaps	%%ymm3,0x20(%%rax)	\n\t"\
				"vmovaps	%%ymm1,0x40(%%rax)	\n\t	vmovaps	%%ymm4,0x60(%%rax)	\n\t"\
				"vmovaps	%%ymm2,0x80(%%rax)	\n\t	vmovaps	%%ymm5,0xa0(%%rax)	\n\t"\
				:					/* outputs: none */\
				: [__aqinv0] "m" (aqinv0)	/* All inputs from memory addresses here */\
				 ,[__two26i] "m" (two26i)	\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6"	/* Clobbered registers */\
			);

		#else	// SSE2:

			__asm__ volatile (\
				"movq	%[__aqinv0],%%rax	\n\t"\
				"movq	%[__two26i],%%rsi	\n\t"\
				"movaps	(%%rsi),%%xmm12		\n\t"\
				"movaps	    (%%rax),%%xmm0	\n\t	movaps	0x10(%%rax),%%xmm3	\n\t	movaps	0x20(%%rax),%%xmm6	\n\t	movaps	0x30(%%rax),%%xmm9	\n\t"\
				"movaps	0x40(%%rax),%%xmm1	\n\t	movaps	0x50(%%rax),%%xmm4	\n\t	movaps	0x60(%%rax),%%xmm7	\n\t	movaps	0x70(%%rax),%%xmm10	\n\t"\
				"movaps	0x80(%%rax),%%xmm2	\n\t	movaps	0x90(%%rax),%%xmm5	\n\t	movaps	0xa0(%%rax),%%xmm8	\n\t	movaps	0xb0(%%rax),%%xmm11	\n\t"\
				"mulpd	%%xmm12,%%xmm0		\n\t	mulpd	%%xmm12,%%xmm3		\n\t	mulpd	%%xmm12,%%xmm6		\n\t	mulpd	%%xmm12,%%xmm9		\n\t"\
				"mulpd	%%xmm12,%%xmm1		\n\t	mulpd	%%xmm12,%%xmm4		\n\t	mulpd	%%xmm12,%%xmm7		\n\t	mulpd	%%xmm12,%%xmm10		\n\t"\
				"mulpd	%%xmm12,%%xmm2		\n\t	mulpd	%%xmm12,%%xmm5		\n\t	mulpd	%%xmm12,%%xmm8		\n\t	mulpd	%%xmm12,%%xmm11		\n\t"\
				"movaps	%%xmm0,    (%%rax)	\n\t	movaps	%%xmm3,0x10(%%rax)	\n\t	movaps	%%xmm6,0x20(%%rax)	\n\t	movaps	%%xmm9 ,0x30(%%rax)	\n\t"\
				"movaps	%%xmm1,0x40(%%rax)	\n\t	movaps	%%xmm4,0x50(%%rax)	\n\t	movaps	%%xmm7,0x60(%%rax)	\n\t	movaps	%%xmm10,0x70(%%rax)	\n\t"\
				"movaps	%%xmm2,0x80(%%rax)	\n\t	movaps	%%xmm5,0x90(%%rax)	\n\t	movaps	%%xmm8,0xa0(%%rax)	\n\t	movaps	%%xmm11,0xb0(%%rax)	\n\t"\
				:					/* outputs: none */\
				: [__aqinv0] "m" (aqinv0)	/* All inputs from memory addresses here */\
				 ,[__two26i] "m" (two26i)	\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12"	/* Clobbered registers */\
			);

		#endif

			/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
			j = start_index-1;

			/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
			LSHIFT96(qinv0, zshift, lo0);			LSHIFT96(qinv4, zshift, lo4);
			LSHIFT96(qinv1, zshift, lo1);			LSHIFT96(qinv5, zshift, lo5);
			LSHIFT96(qinv2, zshift, lo2);			LSHIFT96(qinv6, zshift, lo6);
			LSHIFT96(qinv3, zshift, lo3);			LSHIFT96(qinv7, zshift, lo7);
			/* Only want the lower 14 bits here */
			lo0.d1 &= 0x3fff;			lo4.d1 &= 0x3fff;
			lo1.d1 &= 0x3fff;			lo5.d1 &= 0x3fff;
			lo2.d1 &= 0x3fff;			lo6.d1 &= 0x3fff;
			lo3.d1 &= 0x3fff;			lo7.d1 &= 0x3fff;
			// These share the prod192 variable, so cannot be done in side-by-side fashion:
			// Inline a (... >> 78) specialization of the RSHIFT192 macro here, since have a constant shift count:
			MUL_LOHI96_PROD192(q0,lo0,prod192);	lo0.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo0.d1 = (prod192.d2 >> 14);
			MUL_LOHI96_PROD192(q1,lo1,prod192);	lo1.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo1.d1 = (prod192.d2 >> 14);
			MUL_LOHI96_PROD192(q2,lo2,prod192);	lo2.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo2.d1 = (prod192.d2 >> 14);
			MUL_LOHI96_PROD192(q3,lo3,prod192);	lo3.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo3.d1 = (prod192.d2 >> 14);
			MUL_LOHI96_PROD192(q4,lo4,prod192);	lo4.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo4.d1 = (prod192.d2 >> 14);
			MUL_LOHI96_PROD192(q5,lo5,prod192);	lo5.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo5.d1 = (prod192.d2 >> 14);
			MUL_LOHI96_PROD192(q6,lo6,prod192);	lo6.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo6.d1 = (prod192.d2 >> 14);
			MUL_LOHI96_PROD192(q7,lo7,prod192);	lo7.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo7.d1 = (prod192.d2 >> 14);

			/* hi = 0 in this instance, which simplifies things. */
			/* Put the result in lo (rather than x), to ease overflow check below */
			SUB96(q0, lo0, x0);						SUB96(q4, lo4, x4);
			SUB96(q1, lo1, x1);						SUB96(q5, lo5, x5);
			SUB96(q2, lo2, x2);						SUB96(q6, lo6, x6);
			SUB96(q3, lo3, x3);						SUB96(q7, lo7, x7);

			if((pshift >> j) & (uint32)1)
			{
				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
				if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
				if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
				if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
				if(CMPUGT96(x4, qhalf4)){ ADD96(x4, x4, x4); SUB96(x4, q4, x4); }else{ ADD96(x4, x4, x4); }
				if(CMPUGT96(x5, qhalf5)){ ADD96(x5, x5, x5); SUB96(x5, q5, x5); }else{ ADD96(x5, x5, x5); }
				if(CMPUGT96(x6, qhalf6)){ ADD96(x6, x6, x6); SUB96(x6, q6, x6); }else{ ADD96(x6, x6, x6); }
				if(CMPUGT96(x7, qhalf7)){ ADD96(x7, x7, x7); SUB96(x7, q7, x7); }else{ ADD96(x7, x7, x7); }
			}

			/* Convert x to floating form: */
			CVT_UINT78_3WORD_DOUBLE(x0 ,*ax0,*ax1,*ax2);
			CVT_UINT78_3WORD_DOUBLE(x1 ,*bx0,*bx1,*bx2);
			CVT_UINT78_3WORD_DOUBLE(x2 ,*cx0,*cx1,*cx2);
			CVT_UINT78_3WORD_DOUBLE(x3 ,*dx0,*dx1,*dx2);
			CVT_UINT78_3WORD_DOUBLE(x4 ,*ex0,*ex1,*ex2);
			CVT_UINT78_3WORD_DOUBLE(x5 ,*fx0,*fx1,*fx2);
			CVT_UINT78_3WORD_DOUBLE(x6 ,*gx0,*gx1,*gx2);
			CVT_UINT78_3WORD_DOUBLE(x7 ,*hx0,*hx1,*hx2);

		#ifdef USE_AVX

			__asm__ volatile (\
				"movq	%[__ax0],%%rax	\n\t"\
				"vmovaps	    (%%rax),%%ymm0	\n\t	vmovaps	0x20(%%rax),%%ymm8	\n\t"\
				"vmovaps	0x40(%%rax),%%ymm2	\n\t	vmovaps	0x60(%%rax),%%ymm10	\n\t"\
				"vmovaps	0x80(%%rax),%%ymm4	\n\t	vmovaps	0xa0(%%rax),%%ymm12	\n\t"\
				:				/* outputs: none */\
				: [__ax0] "m" (ax0)		/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
			);

		#else	// SSE2:

			__asm__ volatile (\
				"movq	%[__ax0],%%rax	\n\t"\
				"movaps	    (%%rax),%%xmm0	\n\t	movaps	0x10(%%rax),%%xmm4	\n\t	movaps	0x20(%%rax),%%xmm8	\n\t	movaps	0x30(%%rax),%%xmm12	\n\t"\
				"movaps	0x40(%%rax),%%xmm1	\n\t	movaps	0x50(%%rax),%%xmm5	\n\t	movaps	0x60(%%rax),%%xmm9	\n\t	movaps	0x70(%%rax),%%xmm13	\n\t"\
				"movaps	0x80(%%rax),%%xmm2	\n\t	movaps	0x90(%%rax),%%xmm6	\n\t	movaps	0xa0(%%rax),%%xmm10	\n\t	movaps	0xb0(%%rax),%%xmm14	\n\t"\
				:					/* outputs: none */\
				: [__ax0] "m" (ax0)	/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm12","xmm13","xmm14"	/* Clobbered registers */\
			);

		#endif

			/*...x^2 mod q is returned in x. */
			/* All 3-word-double-form operands have components in the following size ranges:
				fword0,1 in [-2^25, +2^25]
				fword2   in [   -1, +2^26]
			*/

		/* Inner loop body needs 42 movaps, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */
		/* !ASM_LOOP: High-level loop construct, and explicit branch to do the modular doubling: */

			for(j = start_index-2; j >= 0; j--)
			{
			#ifdef USE_AVX
				SSE2_twopmodq78_modmul_q8(aq0,pshift,j);
			#else	// SSE2:
				SSE2_twopmodq78_modmul_q8(aq0,aqinv0,ax0,two26i,pshift,j);
			#endif
			}	/* for(j...) */

		#ifdef USE_AVX

			__asm__ volatile (\
				"movq	%[__ax0],%%rax	\n\t"\
				"vmovaps	%%ymm0,    (%%rax)	\n\t	vmovaps	%%ymm8 ,0x20(%%rax)	\n\t"\
				"vmovaps	%%ymm2,0x40(%%rax)	\n\t	vmovaps	%%ymm10,0x60(%%rax)	\n\t"\
				"vmovaps	%%ymm4,0x80(%%rax)	\n\t	vmovaps	%%ymm12,0xa0(%%rax)	\n\t"\
				:				/* outputs: none */\
				: [__ax0] "m" (ax0)		/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
			);

		#else	// SSE2:

			__asm__ volatile (\
				"movq	%[__ax0],%%rax	\n\t"\
				"movaps	%%xmm0,    (%%rax)	\n\t	movaps	%%xmm4,0x10(%%rax)	\n\t	movaps	%%xmm8 ,0x20(%%rax)	\n\t	movaps	%%xmm12,0x30(%%rax)	\n\t"\
				"movaps	%%xmm1,0x40(%%rax)	\n\t	movaps	%%xmm5,0x50(%%rax)	\n\t	movaps	%%xmm9 ,0x60(%%rax)	\n\t	movaps	%%xmm13,0x70(%%rax)	\n\t"\
				"movaps	%%xmm2,0x80(%%rax)	\n\t	movaps	%%xmm6,0x90(%%rax)	\n\t	movaps	%%xmm10,0xa0(%%rax)	\n\t	movaps	%%xmm14,0xb0(%%rax)	\n\t"\
				:					/* outputs: none */\
				: [__ax0] "m" (ax0)	/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm12","xmm13","xmm14"	/* Clobbered registers */\
			);

		#endif

			/*...Double and return.	These are specialized for the case
			where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
			*/
			CVT78_3WORD_DOUBLE_UINT96(*ax0,*ax1,*ax2, x0);
			CVT78_3WORD_DOUBLE_UINT96(*bx0,*bx1,*bx2, x1);
			CVT78_3WORD_DOUBLE_UINT96(*cx0,*cx1,*cx2, x2);
			CVT78_3WORD_DOUBLE_UINT96(*dx0,*dx1,*dx2, x3);
			CVT78_3WORD_DOUBLE_UINT96(*ex0,*ex1,*ex2, x4);
			CVT78_3WORD_DOUBLE_UINT96(*fx0,*fx1,*fx2, x5);
			CVT78_3WORD_DOUBLE_UINT96(*gx0,*gx1,*gx2, x6);
			CVT78_3WORD_DOUBLE_UINT96(*hx0,*hx1,*hx2, x7);
			/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
			ADD96(x0,x0,x0);			ADD96(x4,x4,x4);
			ADD96(x1,x1,x1);			ADD96(x5,x5,x5);
			ADD96(x2,x2,x2);			ADD96(x6,x6,x6);
			ADD96(x3,x3,x3);			ADD96(x7,x7,x7);
			if(FERMAT) {
				q0.d0 -= 2;		q4.d0 -= 2;
				q1.d0 -= 2;		q5.d0 -= 2;
				q2.d0 -= 2;		q6.d0 -= 2;
				q3.d0 -= 2;		q7.d0 -= 2;
			}
			SUB96(x0,q0,x0);	SUB96(x4,q4,x4);
			SUB96(x1,q1,x1);	SUB96(x5,q5,x5);
			SUB96(x2,q2,x2);	SUB96(x6,q6,x6);
			SUB96(x3,q3,x3);	SUB96(x7,q7,x7);
			tmp0 = CMPEQ96(x0, ONE96);	tmp4 = CMPEQ96(x4, ONE96);
			tmp1 = CMPEQ96(x1, ONE96);	tmp5 = CMPEQ96(x5, ONE96);
			tmp2 = CMPEQ96(x2, ONE96);	tmp6 = CMPEQ96(x6, ONE96);
			tmp3 = CMPEQ96(x3, ONE96);	tmp7 = CMPEQ96(x7, ONE96);
			r  = tmp0;
			r += tmp1 << 1;
			r += tmp2 << 2;
			r += tmp3 << 3;
			r += tmp4 << 4;
			r += tmp5 << 5;
			r += tmp6 << 6;
			r += tmp7 << 7;
			return r;
		}

	  #elif defined(USE_SSE2)

		/************ 64-BIT SSE2 ONLY *** 8-trial-factor, hybrid int/float version **********************/
		/* The integer asm test-code is in twopmodq96.c, also wrapped in the same-named YES_ASM #define. */
		/* The floating-point code is the above 4-factor 64-bit SSE2 code, interleaved with the integer. */
		/*************************************************************************************************/
		uint64 twopmodq78_3WORD_DOUBLE_q8(uint64 p, uint64 k[], int init_sse2, int thr_id)
		{
		#error Needs ||ization before using!
			const char func[] = "twopmodq78_3WORD_DOUBLE_q8 [SSE2]";
			 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
			uint64 lead7, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r;
			static uint64 psave = 0, pshift;
			static uint32 start_index, zshift0,zshift1, first_entry = TRUE;
			/* Floating-point stuff: */
			uint32 q32_0, qinv32_0, tmp32_0
				 , q32_1, qinv32_1, tmp32_1
				 , q32_2, qinv32_2, tmp32_2
				 , q32_3, qinv32_3, tmp32_3
				 , q32_4, qinv32_4, tmp32_4
				 , q32_5, qinv32_5, tmp32_5
				 , q32_6, qinv32_6, tmp32_6
				 , q32_7, qinv32_7, tmp32_7;
			uint96 q0, qinv0, qhalf0, x0, lo0
				 , q1, qinv1, qhalf1, x1, lo1
				 , q2, qinv2, qhalf2, x2, lo2
				 , q3, qinv3, qhalf3, x3, lo3;
			uint192 prod192;
			/* Force 128-bit alignment needed for SSE work via the complex type, but also need double-pointer accesses to indivdual TF data: */
			static vec_dbl *sc_arr = 0x0;
			static double *sc_ptr;
			static double *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2, *flo0,*flo1,*flo2, *fhi0,*fhi1,*fhi2;
			static double *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2, *glo0,*glo1,*glo2, *ghi0,*ghi1,*ghi2;
			static double *hq0,*hq1,*hq2,*hqhi52, *hqinv0,*hqinv1,*hqinv2, *hx0,*hx1,*hx2, *hlo0,*hlo1,*hlo2, *hhi0,*hhi1,*hhi2;
			static double *iq0,*iq1,*iq2,*iqhi52, *iqinv0,*iqinv1,*iqinv2, *ix0,*ix1,*ix2, *ilo0,*ilo1,*ilo2, *ihi0,*ihi1,*ihi2;
			static double *half, *two26f, *two26i, *two13i, *sse2_rnd;
			static uint64 ihalf = 0x3FDfffffffffffffull;	/* Bitfield storing 0.5*(1-epsilon) in IEEE64 format */
			/* Integer stuff: */
			uint96 q4,q5,q6,q7;
			static uint64 *sm_ptr, *ptr64;
			static uint96 *ONE96_PTR
				,*qptr4,*qinv4,*qhalf4,*x4,*lo4,*hi4
				,*qptr5,*qinv5,*qhalf5,*x5,*lo5,*hi5
				,*qptr6,*qinv6,*qhalf6,*x6,*lo6,*hi6
				,*qptr7,*qinv7,*qhalf7,*x7,*lo7,*hi7;
			/* Vars needed for post-powering-loop x * 2^(96-78) modmul of the FP outputs resulting from use of common pshift = p+96 for both quartets of moduli */
			double qdiv0, qmul0, kmul0;
			double qdiv1, qmul1, kmul1;
			double qdiv2, qmul2, kmul2;
			double qdiv3, qmul3, kmul3;
			uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

			if(first_entry || p != psave)
			{
				first_entry = FALSE;
				psave  = p;
				pshift = p + 96;

			/* This pshift is correct for 96-bit modmul, but not for the 78-bit floating-point one. We need to use a common pshift
			for the pshift-bit-pattern-dependent conditional doubling inside the big ASM block to work the same for the integer and floating-point data,
			so the use of this 'wrong' pshift needs the floating-point loop outputs to be post-processed in order to yield the
			correct result, i.e. the one which would result from using pshift = p+78. Here is how that works:

			If we add k to the 'true' pshift, need to take the powering-loop result and multiply by 2^k to get the proper result, mod q.
			*/
				j = leadz64(pshift);
				/* Extract leftmost 7 bits of pshift (if > 77/95, use the leftmost 6) and subtract from 96: */
				lead7 = ((pshift<<j) >> 57);	/* In [64,127] */
				if(lead7 > 77) {
					lead7 >>= 1;	/* Guarantees that lead7 in [39,77] */
					start_index =  64-j-6;	/* Use only the leftmost 6 bits */
				} else {
					start_index =  64-j-7;
				}
				zshift0 = 77 - lead7;	zshift1 = 95 - lead7;	/* zshift0 in [0,38]; zshift1 in [18,56] */
				zshift0 <<= 1;			zshift1 <<= 1;			/* In [0,76]/[36,112]; Doubling the shift count here takes cares of the first SQR_LOHI */
				pshift = ~pshift;
				/* 40 16-byte slots for floats, 16 for ints: */
				sc_arr = ALLOC_VEC_DBL(sc_arr, 40+32);	ASSERT(sc_arr != 0x0, "ERROR: unable to allocate sc_arr!");
				sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);
				ASSERT(((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
				/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 2 to span an SSE register.
				The bytewise address offsets of the left-column pointers (relative to base address fq0) are in the right comment-column: */
																															/* Byte offset */
				fq0    = sc_ptr + 0x00;		gq0    = sc_ptr + 0x01;		hq0    = sc_ptr + 0x02;		iq0    = sc_ptr + 0x03;	/* 0x000 */
				fq1    = sc_ptr + 0x04;		gq1    = sc_ptr + 0x05;		hq1    = sc_ptr + 0x06;		iq1    = sc_ptr + 0x07;	/* 0x020 */
				fq2    = sc_ptr + 0x08;		gq2    = sc_ptr + 0x09;		hq2    = sc_ptr + 0x0a;		iq2    = sc_ptr + 0x0b;	/* 0x040 */
				fqhi52 = sc_ptr + 0x0c;		gqhi52 = sc_ptr + 0x0d;		hqhi52 = sc_ptr + 0x0e;		iqhi52 = sc_ptr + 0x0f;	/* 0x060 */
				fqinv0 = sc_ptr + 0x10;		gqinv0 = sc_ptr + 0x11;		hqinv0 = sc_ptr + 0x12;		iqinv0 = sc_ptr + 0x13;	/* 0x080 */
				fqinv1 = sc_ptr + 0x14;		gqinv1 = sc_ptr + 0x15;		hqinv1 = sc_ptr + 0x16;		iqinv1 = sc_ptr + 0x17;	/* 0x0a0 */
				fqinv2 = sc_ptr + 0x18;		gqinv2 = sc_ptr + 0x19;		hqinv2 = sc_ptr + 0x1a;		iqinv2 = sc_ptr + 0x1b;	/* 0x0c0 */
				fx0    = sc_ptr + 0x1c;		gx0    = sc_ptr + 0x1d;		hx0    = sc_ptr + 0x1e;		ix0    = sc_ptr + 0x1f;	/* 0x0e0 */
				fx1    = sc_ptr + 0x20;		gx1    = sc_ptr + 0x21;		hx1    = sc_ptr + 0x22;		ix1    = sc_ptr + 0x23;	/* 0x100 */
				fx2    = sc_ptr + 0x24;		gx2    = sc_ptr + 0x25;		hx2    = sc_ptr + 0x26;		ix2    = sc_ptr + 0x27;	/* 0x120 */
				flo0   = sc_ptr + 0x28;		glo0   = sc_ptr + 0x29;		hlo0   = sc_ptr + 0x2a;		ilo0   = sc_ptr + 0x2b;	/* 0x140 */
				flo1   = sc_ptr + 0x2c;		glo1   = sc_ptr + 0x2d;		hlo1   = sc_ptr + 0x2e;		ilo1   = sc_ptr + 0x2f;	/* 0x160 */
				flo2   = sc_ptr + 0x30;		glo2   = sc_ptr + 0x31;		hlo2   = sc_ptr + 0x32;		ilo2   = sc_ptr + 0x33;	/* 0x180 */
				fhi0   = sc_ptr + 0x34;		ghi0   = sc_ptr + 0x35;		hhi0   = sc_ptr + 0x36;		ihi0   = sc_ptr + 0x37;	/* 0x1a0 */
				fhi1   = sc_ptr + 0x38;		ghi1   = sc_ptr + 0x39;		hhi1   = sc_ptr + 0x3a;		ihi1   = sc_ptr + 0x3b;	/* 0x1c0 */
				fhi2   = sc_ptr + 0x3c;		ghi2   = sc_ptr + 0x3d;		hhi2   = sc_ptr + 0x3e;		ihi2   = sc_ptr + 0x3f;	/* 0x1e0 */
				two13i = sc_ptr + 0x40;	/* 0x200 */
				two26f = sc_ptr + 0x42;	/* 0x210 */
				two26i = sc_ptr + 0x44;	/* 0x220 */
				sse2_rnd=sc_ptr + 0x46;	/* 0x230 */
				half   = sc_ptr + 0x48;	/* 0x240 */
				/* Can premultiply each of the multiword-mul inputs by 1/sqrt(2^13) due to quadraticity: */
				*two13i++ = TWO13FLINV;		*two13i-- = TWO13FLINV;
				*two26f++ = TWO26FLOAT;		*two26f-- = TWO26FLOAT;
				*two26i++ = TWO26FLINV;		*two26i-- = TWO26FLINV;
				/* SSE2 math = 53-mantissa-bit IEEE double-float: */
				*sse2_rnd++ = 3.0*0x4000000*0x2000000;
				*sse2_rnd-- = 3.0*0x4000000*0x2000000;
				/* We init "half" = 0.5-epsilon here, because emulating FLOOR(x) via DNINT(x-half) requires
				us to always round up if x is a whole number, and our DNINT emulation can round either way if fractional part = 0.5:
				*/
				*half++ = *(double*)&ihalf;	*half-- = *(double*)&ihalf;

				/* Need both float and integer data to share same allocated chunk of memory, so can use a single base/offset scheme to manage both */
				sm_ptr = (uint64*)(sc_ptr + 0x50);	/* Contiguous offset w.r.to last float data above is 0x4a, but start ints at +0x50 for ease: */
				ASSERT((uint32)sm_ptr == ((uint32)sc_ptr +  0x280), "sm_ptr not offset as expected!");
				/* Remember, these are pointers-to-uint128, so need an increment of 2 to span a memory slot: */											/* Byte offsets: */
				qptr4  = (uint96*)(sm_ptr + 0x00);	qptr5  = (uint96*)(sm_ptr + 0x02);	qptr6  = (uint96*)(sm_ptr + 0x04);	qptr7  = (uint96*)(sm_ptr + 0x06);	/* 0x280 */
				qinv4  = (uint96*)(sm_ptr + 0x08);	qinv5  = (uint96*)(sm_ptr + 0x0a);	qinv6  = (uint96*)(sm_ptr + 0x0c);	qinv7  = (uint96*)(sm_ptr + 0x0e);	/* 0x2c0 */
				x4     = (uint96*)(sm_ptr + 0x10);	x5     = (uint96*)(sm_ptr + 0x12);	x6     = (uint96*)(sm_ptr + 0x14);	x7     = (uint96*)(sm_ptr + 0x16);	/* 0x300 */
				lo4    = (uint96*)(sm_ptr + 0x18);	lo5    = (uint96*)(sm_ptr + 0x1a);	lo6    = (uint96*)(sm_ptr + 0x1c);	lo7    = (uint96*)(sm_ptr + 0x1e);	/* 0x340 */
				qhalf4 = (uint96*)(sm_ptr + 0x20);	qhalf5 = (uint96*)(sm_ptr + 0x22);	qhalf6 = (uint96*)(sm_ptr + 0x24);	qhalf7 = (uint96*)(sm_ptr + 0x26);	/* 0x380 */
				hi4    = (uint96*)(sm_ptr + 0x28);	hi5    = (uint96*)(sm_ptr + 0x2a);	hi6    = (uint96*)(sm_ptr + 0x2c);	hi7    = (uint96*)(sm_ptr + 0x2e);	/* 0x3c0 */
				ONE96_PTR = (uint96*)(sm_ptr + 0x30);
				ptr64 = (uint64*)ONE96_PTR;	*ptr64++ = ONE96.d0;	*ptr64-- = ONE96.d1;
			}	/* first_entry */

			ASSERT((p >> 63) == 0, "p must be < 2^63!");
			q0.d0 = q1.d0 = q2.d0 = q3.d0 = q4.d0 = q5.d0 = q6.d0 = q7.d0 = p+p;
			MUL_LOHI64(q0.d0, k[0], q0.d0, q0.d1);
			MUL_LOHI64(q1.d0, k[1], q1.d0, q1.d1);
			MUL_LOHI64(q2.d0, k[2], q2.d0, q2.d1);
			MUL_LOHI64(q3.d0, k[3], q3.d0, q3.d1);
			MUL_LOHI64(q4.d0, k[4], q4.d0, q4.d1);
			MUL_LOHI64(q5.d0, k[5], q5.d0, q5.d1);
			MUL_LOHI64(q6.d0, k[6], q6.d0, q6.d1);
			MUL_LOHI64(q7.d0, k[7], q7.d0, q7.d1);

			q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
			q1.d0 += 1;
			q2.d0 += 1;
			q3.d0 += 1;
			q4.d0 += 1;
			q5.d0 += 1;
			q6.d0 += 1;
			q7.d0 += 1;
			ASSERT((q0.d1 >> 14) == 0, "twopmodq78_q8 : (q0.d1 >> 14) != 0");
			ASSERT((q1.d1 >> 14) == 0, "twopmodq78_q8 : (q1.d1 >> 14) != 0");
			ASSERT((q2.d1 >> 14) == 0, "twopmodq78_q8 : (q2.d1 >> 14) != 0");
			ASSERT((q3.d1 >> 14) == 0, "twopmodq78_q8 : (q3.d1 >> 14) != 0");
			ASSERT((q4.d1 >> 14) == 0, "twopmodq78_q8 : (q4.d1 >> 14) != 0");
			ASSERT((q5.d1 >> 14) == 0, "twopmodq78_q8 : (q5.d1 >> 14) != 0");
			ASSERT((q6.d1 >> 14) == 0, "twopmodq78_q8 : (q6.d1 >> 14) != 0");
			ASSERT((q7.d1 >> 14) == 0, "twopmodq78_q8 : (q7.d1 >> 14) != 0");

			/*****************************************************************************************************/
			/*** From here onward, q0-3 get processed via 78-bit float-based modmul, q4-7 via 96-bit pure-int: ***/
			/*****************************************************************************************************/
			q32_0 = (uint32)q0.d0;
			q32_1 = (uint32)q1.d0;
			q32_2 = (uint32)q2.d0;
			q32_3 = (uint32)q3.d0;

			ptr64 = (uint64*)qptr4;
			*ptr64++ = q4.d0;	*ptr64++ = q4.d1;		q32_4 = (uint32)q4.d0;
			*ptr64++ = q5.d0;	*ptr64++ = q5.d1;		q32_5 = (uint32)q5.d0;
			*ptr64++ = q6.d0;	*ptr64++ = q6.d1;		q32_6 = (uint32)q6.d0;
			*ptr64++ = q7.d0;	*ptr64++ = q7.d1;		q32_7 = (uint32)q7.d0;

			/* Convert q0-3 to floating form: */
			CVT_UINT78_3WORD_DOUBLE(q0 ,*fq0,*fq1,*fq2);
			CVT_UINT78_3WORD_DOUBLE(q1 ,*gq0,*gq1,*gq2);
			CVT_UINT78_3WORD_DOUBLE(q2 ,*hq0,*hq1,*hq2);
			CVT_UINT78_3WORD_DOUBLE(q3 ,*iq0,*iq1,*iq2);

			__asm__ volatile (\
				"movq	%[__fq0],%%rax						\n\t"\
				"movq	%[__two26f],%%rsi					\n\t"\
				"movaps	    (%%rsi),%%xmm8	/* two26f */	\n\t"\
				"movaps	0x10(%%rsi),%%xmm9	/* two26i */	\n\t"\
				"movaps	    (%%rax),%%xmm0	/* fq0 */		\n\t	movaps	0x10(%%rax),%%xmm4	/* hq0 */		\n\t"\
				"movaps	0x20(%%rax),%%xmm1	/* fq1 */		\n\t	movaps	0x30(%%rax),%%xmm5	/* hq1 */		\n\t"\
				"movaps	0x40(%%rax),%%xmm2	/* fq2 */		\n\t	movaps	0x50(%%rax),%%xmm6	/* hq2 */		\n\t"\
				"movaps	%%xmm2,%%xmm3		/* cpy fg2 */	\n\t	movaps	%%xmm6,%%xmm7		/* cpy hi2 */	\n\t"\
				"mulpd	%%xmm8,%%xmm3						\n\t	mulpd	%%xmm8,%%xmm7						\n\t"\
				"addpd	%%xmm1,%%xmm3	/* Hi 52 out bits */\n\t	addpd	%%xmm5,%%xmm7		/* Hi 52 */		\n\t"\
				"mulpd	%%xmm9,%%xmm0						\n\t	mulpd	%%xmm9,%%xmm4						\n\t"\
				"mulpd	%%xmm9,%%xmm1						\n\t	mulpd	%%xmm9,%%xmm5						\n\t"\
				"mulpd	%%xmm9,%%xmm2						\n\t	mulpd	%%xmm9,%%xmm6						\n\t"\
				"movaps	%%xmm0,    (%%rax)	/* fq0/2^26 */	\n\t	movaps	%%xmm4,0x10(%%rax)	/* hq0/2^26 */	\n\t"\
				"movaps	%%xmm1,0x20(%%rax)	/* fq1/2^26 */	\n\t	movaps	%%xmm5,0x30(%%rax)	/* hq1/2^26 */	\n\t"\
				"movaps	%%xmm2,0x40(%%rax)	/* fq2/2^26 */	\n\t	movaps	%%xmm6,0x50(%%rax)	/* hq2/2^26 */	\n\t"\
				"movaps	%%xmm3,0x60(%%rax)	/* fqhi52 */	\n\t	movaps	%%xmm7,0x70(%%rax)	/* hqhi52 */	\n\t"\
				:					/* outputs: none */\
				: [__fq0] "m" (fq0)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
			);

			RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
			RSHIFT_FAST96(q1, 1, qhalf1);
			RSHIFT_FAST96(q2, 1, qhalf2);
			RSHIFT_FAST96(q3, 1, qhalf3);
			/* The 2nd set of 4 needs ptr-version of the shift macro: */
			RSHIFT_FAST96_PTR(qptr4, 1, qhalf4);
			RSHIFT_FAST96_PTR(qptr5, 1, qhalf5);
			RSHIFT_FAST96_PTR(qptr6, 1, qhalf6);
			RSHIFT_FAST96_PTR(qptr7, 1, qhalf7);

			/* Initial seed has low 4 bits of mod-[power-of-2 modmul base] inverse correct: */
			qinv32_0 = (q32_0 + q32_0 + q32_0) ^ (uint32)2;
			qinv32_1 = (q32_1 + q32_1 + q32_1) ^ (uint32)2;
			qinv32_2 = (q32_2 + q32_2 + q32_2) ^ (uint32)2;
			qinv32_3 = (q32_3 + q32_3 + q32_3) ^ (uint32)2;
			qinv32_4 = (q32_4 + q32_4 + q32_4) ^ (uint32)2;
			qinv32_5 = (q32_5 + q32_5 + q32_5) ^ (uint32)2;
			qinv32_6 = (q32_6 + q32_6 + q32_6) ^ (uint32)2;
			qinv32_7 = (q32_7 + q32_7 + q32_7) ^ (uint32)2;

			/* Newton iteration involves repeated steps of form

				qinv = qinv*(2 - q*qinv);

			Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
			defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
			*/
			/* 8 bits: */
			tmp32_0 = q32_0*qinv32_0;
			tmp32_1 = q32_1*qinv32_1;
			tmp32_2 = q32_2*qinv32_2;
			tmp32_3 = q32_3*qinv32_3;
			tmp32_4 = q32_4*qinv32_4;
			tmp32_5 = q32_5*qinv32_5;
			tmp32_6 = q32_6*qinv32_6;
			tmp32_7 = q32_7*qinv32_7;
			qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);
			qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);
			qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);
			qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);
			qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);
			qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);
			qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);
			qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);
			/* 16 bits: */
			tmp32_0 = q32_0*qinv32_0;
			tmp32_1 = q32_1*qinv32_1;
			tmp32_2 = q32_2*qinv32_2;
			tmp32_3 = q32_3*qinv32_3;
			tmp32_4 = q32_4*qinv32_4;
			tmp32_5 = q32_5*qinv32_5;
			tmp32_6 = q32_6*qinv32_6;
			tmp32_7 = q32_7*qinv32_7;
			qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);
			qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);
			qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);
			qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);
			qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);
			qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);
			qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);
			qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);
			/* 32 bits: */
			tmp32_0 = q32_0*qinv32_0;
			tmp32_1 = q32_1*qinv32_1;
			tmp32_2 = q32_2*qinv32_2;
			tmp32_3 = q32_3*qinv32_3;
			tmp32_4 = q32_4*qinv32_4;
			tmp32_5 = q32_5*qinv32_5;
			tmp32_6 = q32_6*qinv32_6;
			tmp32_7 = q32_7*qinv32_7;
			qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);
			qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);
			qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);
			qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);
			qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);
			qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);
			qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);
			qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);
			/* 64 bits: */
			qinv0.d0 = (uint64)qinv32_0;
			qinv1.d0 = (uint64)qinv32_1;
			qinv2.d0 = (uint64)qinv32_2;
			qinv3.d0 = (uint64)qinv32_3;
			qinv4->d0 = (uint64)qinv32_4;
			qinv5->d0 = (uint64)qinv32_5;
			qinv6->d0 = (uint64)qinv32_6;
			qinv7->d0 = (uint64)qinv32_7;
			tmp0 = q0.d0*qinv0.d0;
			tmp1 = q1.d0*qinv1.d0;
			tmp2 = q2.d0*qinv2.d0;
			tmp3 = q3.d0*qinv3.d0;
			tmp4 = qptr4->d0 * qinv4->d0;
			tmp5 = qptr5->d0 * qinv5->d0;
			tmp6 = qptr6->d0 * qinv6->d0;
			tmp7 = qptr7->d0 * qinv7->d0;
			qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
			qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
			qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);
			qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);
			qinv4->d0 = qinv4->d0 * ((uint64)2 - tmp4);
			qinv5->d0 = qinv5->d0 * ((uint64)2 - tmp5);
			qinv6->d0 = qinv6->d0 * ((uint64)2 - tmp6);
			qinv7->d0 = qinv7->d0 * ((uint64)2 - tmp7);
			/* 128 bits - Since low 64 bits will not change, just compute the high 64: */
			MULH64(q0.d0, qinv0.d0, tmp0);
			MULH64(q1.d0, qinv1.d0, tmp1);
			MULH64(q2.d0, qinv2.d0, tmp2);
			MULH64(q3.d0, qinv3.d0, tmp3);
			MULH64(qptr4->d0, qinv4->d0, tmp4);
			MULH64(qptr5->d0, qinv5->d0, tmp5);
			MULH64(qptr6->d0, qinv6->d0, tmp6);
			MULH64(qptr7->d0, qinv7->d0, tmp7);

			qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);
			qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);
			qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + tmp2);
			qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + tmp3);
			qinv4->d1 = -qinv4->d0 * (qptr4->d1 * qinv4->d0 + tmp4);
			qinv5->d1 = -qinv5->d0 * (qptr5->d1 * qinv5->d0 + tmp5);
			qinv6->d1 = -qinv6->d0 * (qptr6->d1 * qinv6->d0 + tmp6);
			qinv7->d1 = -qinv7->d0 * (qptr7->d1 * qinv7->d0 + tmp7);

			/* 78/96 bits: */
			qinv0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
			qinv1.d1 &= 0x0000000000003fff;
			qinv2.d1 &= 0x0000000000003fff;
			qinv3.d1 &= 0x0000000000003fff;
			qinv4->d1 &= 0x00000000ffffffff;	/*  Only want the lower 32 bits here  */
			qinv5->d1 &= 0x00000000ffffffff;
			qinv6->d1 &= 0x00000000ffffffff;
			qinv7->d1 &= 0x00000000ffffffff;

			/* Convert qinv to floating form: */
			CVT_UINT78_3WORD_DOUBLE(qinv0 ,*fqinv0,*fqinv1,*fqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv1 ,*gqinv0,*gqinv1,*gqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv2 ,*hqinv0,*hqinv1,*hqinv2);
			CVT_UINT78_3WORD_DOUBLE(qinv3 ,*iqinv0,*iqinv1,*iqinv2);

				__asm__ volatile (\
					"movq	%[__fqinv0],%%rax	\n\t"\
					"movq	%[__two26i],%%rsi	\n\t"\
					"movaps	(%%rsi),%%xmm6		\n\t"\
					"movaps	    (%%rax),%%xmm0	\n\t	movaps	0x10(%%rax),%%xmm3	\n\t"\
					"movaps	0x20(%%rax),%%xmm1	\n\t	movaps	0x30(%%rax),%%xmm4	\n\t"\
					"movaps	0x40(%%rax),%%xmm2	\n\t	movaps	0x50(%%rax),%%xmm5	\n\t"\
					"mulpd	%%xmm6,%%xmm0		\n\t	mulpd	%%xmm6,%%xmm3		\n\t"\
					"mulpd	%%xmm6,%%xmm1		\n\t	mulpd	%%xmm6,%%xmm4		\n\t"\
					"mulpd	%%xmm6,%%xmm2		\n\t	mulpd	%%xmm6,%%xmm5		\n\t"\
					"movaps	%%xmm0,    (%%rax)	\n\t	movaps	%%xmm3,0x10(%%rax)	\n\t"\
					"movaps	%%xmm1,0x20(%%rax)	\n\t	movaps	%%xmm4,0x30(%%rax)	\n\t"\
					"movaps	%%xmm2,0x40(%%rax)	\n\t	movaps	%%xmm5,0x50(%%rax)	\n\t"\
					:					/* outputs: none */\
					: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
					 ,[__two26i] "m" (two26i)	\
					: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6"	/* Clobbered registers */\
				);

			/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
			j = start_index-1;

			/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
			LSHIFT96(qinv0, zshift0, lo0);	lo0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
			LSHIFT96(qinv1, zshift0, lo1);	lo1.d1 &= 0x0000000000003fff;
			LSHIFT96(qinv2, zshift0, lo2);	lo2.d1 &= 0x0000000000003fff;
			LSHIFT96(qinv3, zshift0, lo3);	lo3.d1 &= 0x0000000000003fff;

			MUL_LOHI96_PROD192(q0,lo0,prod192);	RSHIFT192(prod192,78,prod192);	lo0.d0 = prod192.d0;	lo0.d1 = prod192.d1;
			MUL_LOHI96_PROD192(q1,lo1,prod192);	RSHIFT192(prod192,78,prod192);	lo1.d0 = prod192.d0;	lo1.d1 = prod192.d1;
			MUL_LOHI96_PROD192(q2,lo2,prod192);	RSHIFT192(prod192,78,prod192);	lo2.d0 = prod192.d0;	lo2.d1 = prod192.d1;
			MUL_LOHI96_PROD192(q3,lo3,prod192);	RSHIFT192(prod192,78,prod192);	lo3.d0 = prod192.d0;	lo3.d1 = prod192.d1;

			/* hi = 0 in this instance, which simplifies things. */
			SUB96(q0, lo0, x0);	/* Put the result in lo (rather than x), to ease overflow check below */
			SUB96(q1, lo1, x1);
			SUB96(q2, lo2, x2);
			SUB96(q3, lo3, x3);

			/* Due to the common starting zshift, zshift1 can be >= 96, in which case the hi terms are nonzero and the lo terms = 0: */
			if(zshift1 > 95) {
				x4->d0 = x5->d0 = x6->d0 = x7->d0 = (uint64)1 << (zshift1 - 96);
			} else {
				LSHIFT96_PTR(qinv4, zshift1, lo4);
				LSHIFT96_PTR(qinv5, zshift1, lo5);
				LSHIFT96_PTR(qinv6, zshift1, lo6);
				LSHIFT96_PTR(qinv7, zshift1, lo7);

				MULH96_PTR_q4(  qptr4, lo4, lo4
							  , qptr5, lo5, lo5
							  , qptr6, lo6, lo6
							  , qptr7, lo7, lo7);

				SUB96_PTR(qptr4, lo4, x4);
				SUB96_PTR(qptr5, lo5, x5);
				SUB96_PTR(qptr6, lo6, x6);
				SUB96_PTR(qptr7, lo7, x7);
			}

			if((pshift >> j) & (uint32)1)
			{
				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
				if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
				if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
				if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
				if(CMPUGT96_PTR(x4, qhalf4)){ ADD96_PTR(x4, x4, x4); SUB96_PTR(x4, qptr4, x4); }else{ ADD96_PTR(x4, x4, x4); }
				if(CMPUGT96_PTR(x5, qhalf5)){ ADD96_PTR(x5, x5, x5); SUB96_PTR(x5, qptr5, x5); }else{ ADD96_PTR(x5, x5, x5); }
				if(CMPUGT96_PTR(x6, qhalf6)){ ADD96_PTR(x6, x6, x6); SUB96_PTR(x6, qptr6, x6); }else{ ADD96_PTR(x6, x6, x6); }
				if(CMPUGT96_PTR(x7, qhalf7)){ ADD96_PTR(x7, x7, x7); SUB96_PTR(x7, qptr7, x7); }else{ ADD96_PTR(x7, x7, x7); }
			}

			/* Convert x to floating form: */
			CVT_UINT78_3WORD_DOUBLE(x0 ,*fx0,*fx1,*fx2);
			CVT_UINT78_3WORD_DOUBLE(x1 ,*gx0,*gx1,*gx2);
			CVT_UINT78_3WORD_DOUBLE(x2 ,*hx0,*hx1,*hx2);
			CVT_UINT78_3WORD_DOUBLE(x3 ,*ix0,*ix1,*ix2);

			__asm__ volatile (\
				"movq	%[__fx0],%%rax	\n\t"\
				"movaps	    (%%rax),%%xmm0	\n\t	movaps	0x10(%%rax),%%xmm8	\n\t"\
				"movaps	0x20(%%rax),%%xmm2	\n\t	movaps	0x30(%%rax),%%xmm10	\n\t"\
				"movaps	0x40(%%rax),%%xmm4	\n\t	movaps	0x50(%%rax),%%xmm12	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
			);


				/* Byte offsets for various key pointers in the 2 side-by-side instruction streams: */
			#if 0
				Floating-point:				Integer:
				fq0-2		0x000/20/40		qptr4		0x280
				fqhi52		0x060			qinv4		0x2c0
				fqinv0-2	0x080/a0/c0		x4			0x300
				fx0-2		0x0e0/100/120	lo4			0x340
				flo0-2		0x140/160/180	qhalf4		0x380
				fhi0-2		0x1a0/1c0/1e0	hi4			0x3c0
				two13i		0x200			ONE96_PTR	0x400
				two26f		0x210
				two26i		0x220
				sse2_rnd	0x230
				half		0x240
			#endif
				__asm__ volatile (\
					"movq	%[__fq0],%%rsi	/* Use fq0 as the base address throughout */\n\t"\
					"movslq	%[__start_index], %%rcx		/* for(j = start_index-2; j >= 0; j--) { */\n\t"\
					"subq $2,%%rcx						\n\t"\
					"test %%rcx, %%rcx					\n\t"\
					"jl Loop8End		/* Skip if n < 0 */	\n\t"\
				"Loop8Beg:								\n\t"\
					"movaps	0x200(%%rsi),%%xmm6	/* two13i shared between both columns */	\n\t		/* SQR_LOHI96_q4(x*, lo*, hi*): */	\n\t"\
					"																						/* Load the x.d0 data: */	\n\t"\
					"/* SQR_LOHI78_3WORD_DOUBLE_q4(): */\n\t													movq	0x300(%%rsi),%%rax	/* Dereference the x0 pointer... */\n\t"\
					"/* fx0,1,2 assumed in xmm0,2,4 on loop entry */\n\t										mulq	%%rax			/* And square the low 64 bits. */\n\t"\
					"mulpd	%%xmm6,%%xmm0			\n\t	mulpd	%%xmm6 ,%%xmm8			\n\t				movq	%%rax,%%r8	/* lo0 */	\n\t"\
					"																							movq	%%rdx,%%r12	/* hi0 */	\n\t"\
					"mulpd	%%xmm6,%%xmm2			\n\t	mulpd	%%xmm6 ,%%xmm10			\n\t				movq	0x310(%%rsi),%%rax	\n\t"\
					"mulpd	%%xmm6,%%xmm4			\n\t	mulpd	%%xmm6 ,%%xmm12			\n\t				mulq	%%rax		\n\t"\
					"																							movq	%%rax,%%r9	/* lo1 */	\n\t"\
					"movaps	%%xmm0,%%xmm1			\n\t	movaps	%%xmm8 ,%%xmm9			\n\t				movq	%%rdx,%%r13	/* hi1 */	\n\t"\
					"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t				movq	0x320(%%rsi),%%rax	\n\t"\
					"addpd	%%xmm1,%%xmm1			\n\t	addpd	%%xmm9 ,%%xmm9			\n\t				mulq	%%rax		\n\t"\
					"																							movq	%%rax,%%r10	/* lo2 */	\n\t"\
					"addpd	%%xmm3,%%xmm3			\n\t	addpd	%%xmm11,%%xmm11			\n\t				movq	%%rdx,%%r14	/* hi2 */	\n\t"\
					"movaps	%%xmm1,%%xmm5			\n\t	movaps	%%xmm9 ,%%xmm13			\n\t				movq	0x330(%%rsi),%%rax	\n\t"\
					"																							mulq	%%rax		\n\t"\
					"mulpd	%%xmm0,%%xmm0			\n\t	mulpd	%%xmm8 ,%%xmm8			\n\t				movq	%%rax,%%r11	/* lo3 */	\n\t"\
					"mulpd	%%xmm2,%%xmm1			\n\t	mulpd	%%xmm10,%%xmm9			\n\t				movq	%%rdx,%%r15	/* hi3 */	\n\t"\
					"																						/* (lo0, hi0): */	\n\t"\
					"mulpd	%%xmm2,%%xmm2			\n\t	mulpd	%%xmm10,%%xmm10			\n\t				movq	0x300(%%rsi),%%rax	/* Move xlo into a-reg in prep for the mul */	\n\t"\
					"mulpd	%%xmm4,%%xmm5			\n\t	mulpd	%%xmm12,%%xmm13			\n\t				movl	0x308(%%rsi),%%edi	\n\t"\
					"																							shlq	$1,%%rdi		\n\t"\
					"mulpd	%%xmm4,%%xmm3			\n\t	mulpd	%%xmm12,%%xmm11			\n\t				mulq	%%rdi	/* 2*lo*hi in rax:rdx */	\n\t"\
					"mulpd	%%xmm4,%%xmm4			\n\t	mulpd	%%xmm12,%%xmm12			\n\t				shrq	$1,%%rdi		\n\t"\
					"																							imulq	%%rdi,%%rdi	/* hi*hi is just 64-bits, do in place */	\n\t"\
					"/* Move this part of Digit 2 computation here to free up xmm5,13: */	\n\t				addq	%%rax,%%r12	/* Add mid words, result in r12, CF bit set to indicate carryout */ \n\t"\
					"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10			\n\t				adcq	%%rdx,%%rdi	/* Add hi words + carryin, result in rdi, CF bit will be cleared */	\n\t"\
					"/* Digit 0: */					\n\t														movq	%%r12,0x300(%%rsi)	/*Dump back into x (rather than hi) */\n\t"\
					"movaps 0x210(%%rsi),%%xmm5		/* xmm5 = two26f */\n\t										movq	%%rdi,0x308(%%rsi)	\n\t"\
					"movaps 0x220(%%rsi),%%xmm15	/* xmm15 = two26i */\n\t								/* (lo1, hi1): */	\n\t"\
					"																							movq	0x310(%%rsi),%%rax	\n\t"\
					"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t				movl	0x318(%%rsi),%%edi	\n\t"\
					"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14		\n\t				shlq	$1,%%rdi		\n\t"\
					"																							mulq	%%rdi		\n\t"\
					"																							shrq	$1,%%rdi		\n\t"\
					"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t				imulq	%%rdi,%%rdi	\n\t"\
					"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5,%%xmm8			\n\t				addq	%%rax,%%r13	\n\t"\
					"																							adcq	%%rdx,%%rdi	\n\t"\
					"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t				movq	%%r13,0x310(%%rsi)	\n\t"\
					"/* Digit 1: */					\n\t														movq	%%rdi,0x318(%%rsi)	\n\t"\
					"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t			/* (lo2, hi2): */	\n\t"\
					"																							movq	0x320(%%rsi),%%rax	\n\t"\
					"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t				movl	0x328(%%rsi),%%edi	\n\t"\
					"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14		\n\t				shlq	$1,%%rdi		\n\t"\
					"																							mulq	%%rdi		\n\t"\
					"																							shrq	$1,%%rdi		\n\t"\
					"subpd	%%xmm6,%%xmm1			\n\t	subpd	%%xmm14,%%xmm9			\n\t				imulq	%%rdi,%%rdi	\n\t"\
					"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5,%%xmm9			\n\t				addq	%%rax,%%r14	\n\t"\
					"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t				adcq	%%rdx,%%rdi	\n\t"\
					"																							movq	%%r14,0x320(%%rsi)	\n\t"\
					"/* Digit 2: Require both hi and lo half of output to be >= 0, unbalanced: */\n\t			movq	%%rdi,0x328(%%rsi)	\n\t"\
					"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t			/* (lo3, hi3): */	\n\t"\
					"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t				movq	0x330(%%rsi),%%rax	\n\t"\
					"																							movl	0x338(%%rsi),%%edi	\n\t"\
					"subpd	0x240(%%rsi),%%xmm6		\n\t	subpd	0x240(%%rsi),%%xmm14	\n\t				shlq	$1,%%rdi		\n\t"\
					"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14		\n\t				mulq	%%rdi		\n\t"\
					"																							shrq	$1,%%rdi		\n\t"\
					"																							imulq	%%rdi,%%rdi	\n\t"\
					"subpd	%%xmm6,%%xmm2			\n\t	subpd	%%xmm14,%%xmm10			\n\t				addq	%%rax,%%r15	\n\t"\
					"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t				adcq	%%rdx,%%rdi	\n\t"\
					"																							movq	%%r15,0x330(%%rsi)	\n\t"\
					"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t				movq	%%rdi,0x338(%%rsi)	\n\t"\
					"/* Digit 3: */					\n\t												/* MULL96_q4((qinv*, lo*, lo*): */	\n\t"\
					"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t			/* lo0 * qinv0: */	\n\t"\
					"movaps	%%xmm3,%%xmm6			\n\t	movaps	%%xmm11,%%xmm14			\n\t				movq	0x2c0(%%rsi),%%rdi	/* qinv0 */\n\t"\
					"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14		\n\t				movq	%%rdi,%%rax	\n\t"\
					"																							mulq	%%r8	/* Use rax:rdx as product accumulator */\n\t"\
					"																							imulq	%%rdi,%%r12			/* qinv.lo * lo.hi */\n\t"\
					"subpd	%%xmm6,%%xmm3			\n\t	subpd	%%xmm14,%%xmm11			\n\t				addq	%%r12,%%rdx	\n\t"\
					"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t				imulq	0x2c8(%%rsi),%%r8	/* qinv.hi * lo.lo */\n\t"\
					"																							addq	%%r8 ,%%rdx	\n\t"\
					"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t				movq	%%rax,%%r12	/* Move low 64 bits into free reg */\n\t"\
					"/* Digit 4: */					\n\t														movl	%%edx,0x348(%%rsi)	/* Write hi 32 bits into [__lo0] */\n\t"\
					"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t			/* lo1 * qinv1: */	\n\t"\
					"movaps	%%xmm4,%%xmm6			\n\t	movaps	%%xmm12,%%xmm14			\n\t				movq	0x2d0(%%rsi),%%rdi	\n\t"\
					"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14		\n\t				movq	%%rdi,%%rax	\n\t"\
					"																							mulq	%%r9		\n\t"\
					"																							imulq	%%rdi,%%r13	\n\t"\
					"subpd	%%xmm6,%%xmm4			\n\t	subpd	%%xmm14,%%xmm12			\n\t				addq	%%r13,%%rdx	\n\t"\
					"mulpd	%%xmm5,%%xmm4			\n\t	mulpd	%%xmm5 ,%%xmm12			\n\t				imulq	0x2d8(%%rsi),%%r9	\n\t"\
					"																							addq	%%r9 ,%%rdx	\n\t"\
					"/* Digit 5 = the carry. flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6 */	\n\t				movq	%%rax,%%r13	\n\t"\
					"mulpd	%%xmm5,%%xmm6			\n\t	mulpd	%%xmm5 ,%%xmm14			\n\t				movl	%%edx,0x358(%%rsi)	\n\t"\
					"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t			/* lo2 * qinv2: */	\n\t"\
					"movaps	%%xmm3,0x1a0(%%rsi)		\n\t	movaps	%%xmm11,0x1b0(%%rsi)	\n\t				movq	0x2e0(%%rsi),%%rdi	\n\t"\
					"																							movq	%%rdi,%%rax	\n\t"\
					"movaps	%%xmm4,0x1c0(%%rsi)		\n\t	movaps	%%xmm12,0x1d0(%%rsi)	\n\t				mulq	%%r10		\n\t"\
					"/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */						\n\t				imulq	%%rdi,%%r14	\n\t"\
					"/* Digit 0: */					\n\t														addq	%%r14,%%rdx	\n\t"\
					"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t				imulq	0x2e8(%%rsi),%%r10	\n\t"\
					"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t				addq	%%r10,%%rdx	\n\t"\
					"mulpd	0x80(%%rsi),%%xmm0		\n\t	mulpd	0x90(%%rsi),%%xmm8		\n\t				movq	%%rax,%%r14	\n\t"\
					"																							movl	%%edx,0x368(%%rsi)	\n\t"\
					"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t			/* lo3 * qinv3: */	\n\t"\
					"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14		\n\t				movq	0x2f0(%%rsi),%%rdi	\n\t"\
					"																							movq	%%rdi,%%rax	\n\t"\
					"																							mulq	%%r11		\n\t"\
					"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t				imulq	%%rdi,%%r15	\n\t"\
					"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t				addq	%%r15,%%rdx	\n\t"\
					"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t				imulq	0x2f8(%%rsi),%%r11	\n\t"\
					"																							addq	%%r11,%%rdx	\n\t"\
					"/* Digit 1: */					\n\t														movq	%%rax,%%r15	\n\t"\
					"mulpd	0xa0(%%rsi),%%xmm3		\n\t	mulpd	0xb0(%%rsi),%%xmm11		\n\t				movl	%%edx,0x378(%%rsi)	\n\t"\
					"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t		/* MULH96_q4((q*, lo*, lo*): Low 64 bits of of lo0-3 in r12-15: */	\n\t"\
					"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t			/* q0 * lo0: */	\n\t"\
					"																							movq	0x280(%%rsi),%%rax	/* q0 */\n\t"\
					"mulpd	0x80(%%rsi),%%xmm1		\n\t	mulpd	0x90(%%rsi),%%xmm9		\n\t				mulq	%%r12	/* lo.lo*q.lo in rax:rdx */\n\t"\
					"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t				movl	0x288(%%rsi),%%edi	\n\t"\
					"movaps	%%xmm1,%%xmm3			\n\t	movaps	%%xmm9 ,%%xmm11			\n\t				movq	%%rdx,%%r8	/* Discard low 64 bits [rax] */\n\t"\
					"																							movq	%%r12,%%rax	\n\t"\
					"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd	$0,%%xmm11,%%xmm11		\n\t				mulq	%%rdi	/* lo.lo*q.hi in rax:rdx */\n\t"\
					"																							xorq	%%r12,%%r12	\n\t"\
					"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t				addq	%%rax,%%r8	\n\t"\
					"																							adcq	%%rdx,%%r12	\n\t"\
					"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t				movl	0x348(%%rsi),%%eax	/* Cannot do imulq __lohi,[reg] because that does a sign-extended load of __lohi */\n\t"\
					"mulpd	%%xmm15,%%xmm3			\n\t	mulpd	%%xmm15,%%xmm11			\n\t				imulq	%%rax,%%rdi			/* ...so first load __lohi into low 32 bits of a-reg, then compute lo.hi*q.hi. */\n\t"\
					"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t			addq	%%rdi,%%r12	\n\t"\
					"mulpd	0x80(%%rsi),%%xmm2		\n\t	mulpd	0x90(%%rsi),%%xmm10		\n\t				mulq	0x280(%%rsi)		/* q.lo*lo.hi in rax:rdx */\n\t"\
					"																							addq	%%rax,%%r8	\n\t"\
					"mulpd	0xa0(%%rsi),%%xmm6		\n\t	mulpd	0xb0(%%rsi),%%xmm14		\n\t				adcq	%%rdx,%%r12	\n\t"\
					"mulpd	0xc0(%%rsi),%%xmm4		\n\t	mulpd	0xd0(%%rsi),%%xmm12		\n\t	/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
					"																							movq	0x300(%%rsi),%%rax	/* h.d0 */\n\t"\
					"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t				movq	0x308(%%rsi),%%rdx	/* h.d1 */\n\t"\
					"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11			\n\t				subq	%%r8 ,%%rax	/* rax = (h-l).d0 */\n\t"\
					"addpd	%%xmm3,%%xmm2			\n\t	addpd	%%xmm11,%%xmm10			\n\t				sbbq	%%r12,%%rdx	/* rdx = (h-l).d1 */\n\t"\
					"																							sbbq	%%rdi,%%rdi	\n\t"\
					"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t				shrdq	$32,%%rdx,%%rax	/* Only keep high 96-bit output */\n\t"\
					"subpd	0x240(%%rsi),%%xmm3		\n\t	subpd	0x240(%%rsi),%%xmm11	\n\t				shrq	$32,%%rdx			\n\t"\
					"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd	$0,%%xmm11,%%xmm11		\n\t				movq	0x280(%%rsi),%%r8	\n\t"\
					"																							movq	0x288(%%rsi),%%r12	\n\t"\
					"																							andq	%%rdi,%%r8	\n\t"\
					"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t				andq	%%rdi,%%r12	\n\t"\
					"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t				addq	%%rax,%%r8	\n\t"\
					"																							adcq	%%rdx,%%r12	\n\t"\
					"/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */								\n\t			/* q1 * lo1: */	\n\t"\
					"/* Digit 0: */					\n\t														movq	0x290(%%rsi),%%rax	\n\t"\
					"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t				mulq	%%r13	\n\t"\
					"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t				movl	0x298(%%rsi),%%edi	\n\t"\
					"mulpd	    (%%rsi),%%xmm0		\n\t	mulpd	0x10(%%rsi),%%xmm8		\n\t				movq	%%rdx,%%r9	\n\t"\
					"																							movq	%%r13,%%rax	\n\t"\
					"roundpd	$0,%%xmm0,%%xmm0	\n\t	roundpd	$0,%%xmm8 ,%%xmm8		\n\t				mulq	%%rdi	\n\t"\
					"																							xorq	%%r13,%%r13	\n\t"\
					"mulpd	    %%xmm15,%%xmm0		\n\t	mulpd	    %%xmm15,%%xmm8		\n\t				addq	%%rax,%%r9	\n\t"\
					"/* Digit 1: */					\n\t														adcq	%%rdx,%%r13	\n\t"\
					"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t				movl	0x358(%%rsi),%%eax	\n\t"\
					"mulpd	0x20(%%rsi),%%xmm3		\n\t	mulpd	0x30(%%rsi),%%xmm11		\n\t				imulq	%%rax,%%rdi	\n\t"\
					"																							addq	%%rdi,%%r13	\n\t"\
					"addpd	%%xmm0,%%xmm3			\n\t	addpd	%%xmm8 ,%%xmm11			\n\t				mulq	0x290(%%rsi)	\n\t"\
					"																							addq	%%rax,%%r9	\n\t"\
					"mulpd	    (%%rsi),%%xmm1		\n\t	mulpd	0x10(%%rsi),%%xmm9		\n\t				adcq	%%rdx,%%r13	\n\t"\
					"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t	/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
					"roundpd	$0,%%xmm1,%%xmm1	\n\t	roundpd	$0,%%xmm9 ,%%xmm9		\n\t				movq	0x310(%%rsi),%%rax	\n\t"\
					"																							movq	0x318(%%rsi),%%rdx	\n\t"\
					"																							subq	%%r9 ,%%rax		\n\t"\
					"mulpd	    %%xmm15,%%xmm1		\n\t	mulpd	    %%xmm15,%%xmm9		\n\t				sbbq	%%r13,%%rdx		\n\t"\
					"																							sbbq	%%rdi,%%rdi	/* bitmask for q */\n\t"\
					"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t			shrdq	$32,%%rdx,%%rax		\n\t"\
					"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t				shrq	$32,%%rdx			\n\t"\
					"movaps	%%xmm6,%%xmm0			\n\t	movaps	%%xmm14,%%xmm8			\n\t				movq	0x290(%%rsi),%%r9	\n\t"\
					"mulpd	    (%%rsi),%%xmm2		\n\t	mulpd	0x10(%%rsi),%%xmm10		\n\t				movq	0x298(%%rsi),%%r13	\n\t"\
					"																							andq	%%rdi,%%r9	\n\t"\
					"mulpd	0x20(%%rsi),%%xmm6		\n\t	mulpd	0x30(%%rsi),%%xmm14		\n\t				andq	%%rdi,%%r13	\n\t"\
					"mulpd	0x40(%%rsi),%%xmm4		\n\t	mulpd	0x50(%%rsi),%%xmm12		\n\t				addq	%%rax,%%r9	\n\t"\
					"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t				adcq	%%rdx,%%r13	\n\t"\
					"																						/* q2 * lo2: */	\n\t"\
					"addpd	%%xmm4,%%xmm1			\n\t	addpd	%%xmm12,%%xmm9			\n\t				movq	0x2a0(%%rsi),%%rax	\n\t"\
					"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm9 ,%%xmm10			\n\t				mulq	%%r14	\n\t"\
					"subpd	0x240(%%rsi),%%xmm2		\n\t	subpd	0x240(%%rsi),%%xmm10	\n\t				movl	0x2a8(%%rsi),%%edi	\n\t"\
					"																							movq	%%rdx,%%r10	\n\t"\
					"roundpd	$0,%%xmm2,%%xmm2	\n\t	roundpd	$0,%%xmm10,%%xmm10		\n\t				movq	%%r14,%%rax	\n\t"\
					"																							mulq	%%rdi	\n\t"\
					"mulpd	    %%xmm15,%%xmm2		\n\t	mulpd	    %%xmm15,%%xmm10		\n\t				xorq	%%r14,%%r14	\n\t"\
					"																							addq	%%rax,%%r10	\n\t"\
					"/* Precompute all the needed partial products: */						\n\t				adcq	%%rdx,%%r14	\n\t"\
					"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t				movl	0x368(%%rsi),%%eax	\n\t"\
					"mulpd	0x40(%%rsi),%%xmm0		\n\t	mulpd	0x50(%%rsi),%%xmm8		\n\t				imulq	%%rax,%%rdi	\n\t"\
					"																							addq	%%rdi,%%r14	\n\t"\
					"mulpd	0x20(%%rsi),%%xmm1		\n\t	mulpd	0x30(%%rsi),%%xmm9		\n\t				mulq	0x2a0(%%rsi)	\n\t"\
					"mulpd	0x40(%%rsi),%%xmm3		\n\t	mulpd	0x50(%%rsi),%%xmm11		\n\t				addq	%%rax,%%r10	\n\t"\
					"																							adcq	%%rdx,%%r14	\n\t"\
					"/* Digit 3: */					\n\t											/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
					"addpd	%%xmm2,%%xmm1			\n\t	addpd	%%xmm10,%%xmm9			\n\t				movq	0x320(%%rsi),%%rax	\n\t"\
					"addpd	%%xmm1,%%xmm0			\n\t	addpd	%%xmm9 ,%%xmm8			\n\t				movq	0x328(%%rsi),%%rdx	\n\t"\
					"																							subq	%%r10,%%rax		\n\t"\
					"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t				sbbq	%%r14,%%rdx		\n\t"\
					"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14		\n\t				sbbq	%%rdi,%%rdi	/* bitmask for q */\n\t"\
					"																							shrdq	$32,%%rdx,%%rax		\n\t"\
					"																							shrq	$32,%%rdx			\n\t"\
					"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t				movq	0x2a0(%%rsi),%%r10	\n\t"\
					"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t				movq	0x2a8(%%rsi),%%r14	\n\t"\
					"mulpd	    %%xmm15,%%xmm6		\n\t	mulpd	    %%xmm15,%%xmm14		\n\t				andq	%%rdi,%%r10	\n\t"\
					"																							andq	%%rdi,%%r14	\n\t"\
					"/* Digit 4: */					\n\t														addq	%%rax,%%r10	\n\t"\
					"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t				adcq	%%rdx,%%r14	\n\t"\
					"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t			/* q3 * lo3: */	\n\t"\
					"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd	$0,%%xmm11,%%xmm11		\n\t				movq	0x2b0(%%rsi),%%rax	\n\t"\
					"																							mulq	%%r15	\n\t"\
					"																							movl	0x2b8(%%rsi),%%edi	\n\t"\
					"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t				movq	%%rdx,%%r11	\n\t"\
					"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t				movq	%%r15,%%rax	\n\t"\
					"/* If h < l, calculate h-l+q; otherwise h-l. */						\n\t				mulq	%%rdi	\n\t"\
					"/* Use leading 52 bits to approximate full 78-bit compare. Result in [0, q). */\n\t		xorq	%%r15,%%r15	\n\t"\
					"movaps %%xmm5,%%xmm13	/* Need a copy of two26f, stored in %%xmm5 */	\n\t				addq	%%rax,%%r11	\n\t"\
					"movaps	0x1c0(%%rsi),%%xmm2		\n\t	movaps	0x1d0(%%rsi),%%xmm10	\n\t				adcq	%%rdx,%%r15	\n\t"\
					"																							movl	0x378(%%rsi),%%eax	\n\t"\
					"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t				imulq	%%rax,%%rdi	\n\t"\
					"movaps	0x60(%%rsi),%%xmm4		\n\t	movaps	0x70(%%rsi),%%xmm12		\n\t				addq	%%rdi,%%r15	\n\t"\
					"																							mulq	0x2b0(%%rsi)	\n\t"\
					"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t				addq	%%rax,%%r11	\n\t"\
					"mulpd	    (%%rsi),%%xmm5		\n\t	mulpd	0x10(%%rsi),%%xmm13		\n\t				adcq	%%rdx,%%r15	\n\t"\
					"addpd	%%xmm1,%%xmm3			\n\t	addpd	%%xmm9 ,%%xmm11			\n\t	/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
					"																							movq	0x330(%%rsi),%%rax	\n\t"\
					"movaps	0x1a0(%%rsi),%%xmm1		\n\t	movaps	0x1b0(%%rsi),%%xmm9		\n\t				movq	0x338(%%rsi),%%rdx	\n\t"\
					"cmppd	$0x1,%%xmm3,%%xmm6		\n\t	cmppd	$0x1,%%xmm11,%%xmm14	\n\t				subq	%%r11,%%rax		\n\t"\
					"subpd	%%xmm0,%%xmm1			\n\t	subpd	%%xmm8 ,%%xmm9			\n\t				sbbq	%%r15,%%rdx		\n\t"\
					"																							sbbq	%%rdi,%%rdi	/* bitmask for q */\n\t"\
					"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t				shrdq	$32,%%rdx,%%rax		\n\t"\
					"movaps	%%xmm4,%%xmm0			\n\t	movaps	%%xmm12,%%xmm8			\n\t				shrq	$32,%%rdx			\n\t"\
					"andpd	%%xmm6,%%xmm4			\n\t	andpd	%%xmm14,%%xmm12			\n\t				movq	0x2b0(%%rsi),%%r11	\n\t"\
					"																							movq	0x2b8(%%rsi),%%r15	\n\t"\
					"addpd	%%xmm4,%%xmm2			\n\t	addpd	%%xmm12,%%xmm10			\n\t				andq	%%rdi,%%r11	\n\t"\
					"andpd	%%xmm5,%%xmm6			\n\t	andpd	%%xmm13,%%xmm14			\n\t				andq	%%rdi,%%r15	\n\t"\
					"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t				addq	%%rax,%%r11	\n\t"\
					"/* qlo26 in xmm5, qhi52 in xmm0 */															adcq	%%rdx,%%r15	\n\t"\
					"/* if((pshift >> j) & (uint64)1): */		\n\t"\
					"movq	%[__pshift],%%rax					\n\t"\
					"shrq	%%cl,%%rax	/* j already in c-reg */\n\t"\
					"andq	$0x1,%%rax							\n\t"\
				"je	twopmodq96_q4_pshiftjmp						\n\t"\
					"			/* Int64 code: if h<l carryout of low 64 bits gives hi=2^32 = 0x100000000, need to zero upper 32 bits prior to double step: */\n\t"\
					"																							movq	$-1,%%rdi	\n\t"\
					"																							shrq	$32,%%rdi	\n\t"\
					"																							andq	%%rdi,%%r12	\n\t"\
					"																							andq	%%rdi,%%r13	\n\t"\
					"																							andq	%%rdi,%%r14	\n\t"\
					"																							andq	%%rdi,%%r15	\n\t"\
					"/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */\n\t"\
					"	movaps		%%xmm0,%%xmm6											\n\t			/* x0: */	\n\t"\
					"																							addq	%%r8 ,%%r8	\n\t"\
					"																							adcq	%%r12,%%r12	\n\t"\
					"	movaps		%%xmm8 ,%%xmm14	/* cpy of qhi */						\n\t				movq	0x280(%%rsi),%%rax	\n\t"\
					"																							movq	0x288(%%rsi),%%rdx	\n\t"\
					"																							subq	%%rax,%%r8		\n\t"\
					"	addpd		%%xmm2,%%xmm2											\n\t				sbbq	%%rdx,%%r12		\n\t"\
					"																							sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
					"																							andq	%%rdi,%%rax	\n\t"\
					"	addpd		%%xmm10,%%xmm10	/* top 52 bits */						\n\t				andq	%%rdi,%%rdx	\n\t"\
					"																							addq	%%rax,%%r8	\n\t"\
					"																							adcq	%%rdx,%%r12	\n\t"\
					"	addpd		%%xmm1,%%xmm1											\n\t			/* x1: */	\n\t"\
					"																							addq	%%r9 ,%%r9	\n\t"\
					"																							adcq	%%r13,%%r13	\n\t"\
					"	addpd		%%xmm9 ,%%xmm9	/* low 26 bits */						\n\t				movq	0x290(%%rsi),%%rax	\n\t"\
					"																							movq	0x298(%%rsi),%%rdx	\n\t"\
					"/* If x > q, subtract q: */											\n\t				subq	%%rax,%%r9		\n\t"\
					"	cmppd	$0x2,%%xmm2,%%xmm6											\n\t				sbbq	%%rdx,%%r13		\n\t"\
					"																							sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
					"																							andq	%%rdi,%%rax	\n\t"\
					"	cmppd	$0x2,%%xmm10,%%xmm14/* bitmask = (qhi <= xhi) */			\n\t				andq	%%rdi,%%rdx	\n\t"\
					"																							addq	%%rax,%%r9	\n\t"\
					"																							adcq	%%rdx,%%r13	\n\t"\
					"	andpd		%%xmm6,%%xmm0											\n\t			/* x2: */	\n\t"\
					"																							addq	%%r10,%%r10	\n\t"\
					"																							adcq	%%r14,%%r14	\n\t"\
					"	andpd		%%xmm14,%%xmm8	/* qhi52 & bitmask */					\n\t				movq	0x2a0(%%rsi),%%rax	\n\t"\
					"																							movq	0x2a8(%%rsi),%%rdx	\n\t"\
					"																							subq	%%rax,%%r10		\n\t"\
					"	andpd		%%xmm6,%%xmm5											\n\t				sbbq	%%rdx,%%r14		\n\t"\
					"																							sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
					"																							andq	%%rdi,%%rax	\n\t"\
					"	andpd		%%xmm14,%%xmm13	/* qlo26 & bitmask */					\n\t				andq	%%rdi,%%rdx	\n\t"\
					"																							addq	%%rax,%%r10	\n\t"\
					"																							adcq	%%rdx,%%r14	\n\t"\
					"	subpd		%%xmm0,%%xmm2											\n\t			/* x3: */	\n\t"\
					"																							addq	%%r11,%%r11	\n\t"\
					"																							adcq	%%r15,%%r15	\n\t"\
					"	subpd		%%xmm8 ,%%xmm10	/* x mod q, top 52 bits */				\n\t				movq	0x2b0(%%rsi),%%rax	\n\t"\
					"																							movq	0x2b8(%%rsi),%%rdx	\n\t"\
					"																							subq	%%rax,%%r11		\n\t"\
					"	subpd		%%xmm5,%%xmm1											\n\t				sbbq	%%rdx,%%r15		\n\t"\
					"																							sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
					"																							andq	%%rdi,%%rax	\n\t"\
					"	subpd		%%xmm13,%%xmm9	/* x mod q, low 26 bits */				\n\t				andq	%%rdi,%%rdx	\n\t"\
					"																							addq	%%rax,%%r11	\n\t"\
					"																							adcq	%%rdx,%%r15	\n\t"\
				"twopmodq96_q4_pshiftjmp:													\n\t"\
					"/* } */																\n\t"\
					"/* Normalize the result: */											\n\t"\
					"movaps	0x210(%%rsi),%%xmm4		\n\t	movaps	0x210(%%rsi),%%xmm12	\n\t"\
					"movaps	     %%xmm15,%%xmm3		\n\t	movaps	     %%xmm15,%%xmm11	\n\t				movq	%%r8 ,0x300(%%rsi)	/* Write lo 64 bits into [__x.d0] */\n\t"\
					"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm11,%%xmm9			\n\t"\
					"mulpd	%%xmm3,%%xmm2			\n\t	mulpd	%%xmm11,%%xmm10			\n\t				movq	%%r12,0x308(%%rsi)	/* Write hi 32 bits into [__x.d1] */\n\t"\
					"movaps	%%xmm1,%%xmm0			\n\t	movaps	%%xmm9 ,%%xmm8			\n\t"\
					"roundpd	$0,%%xmm1,%%xmm1	\n\t	roundpd	$0,%%xmm9,%%xmm9		\n\t				movq	%%r9 ,0x310(%%rsi)	\n\t"\
					"subpd	%%xmm1,%%xmm0			\n\t	subpd	%%xmm9 ,%%xmm8			\n\t				movq	%%r13,0x318(%%rsi)	\n\t"\
					"mulpd	%%xmm4,%%xmm0			\n\t	mulpd	%%xmm12,%%xmm8			\n\t"\
					"mulpd	%%xmm3,%%xmm1			\n\t	mulpd	%%xmm11,%%xmm9			\n\t				movq	%%r10,0x320(%%rsi)	\n\t"\
					"addpd	%%xmm1,%%xmm2			\n\t	addpd	%%xmm9 ,%%xmm10			\n\t"\
					"movaps	%%xmm2,%%xmm1			\n\t	movaps	%%xmm10,%%xmm9			\n\t				movq	%%r14,0x328(%%rsi)	\n\t"\
					"roundpd	$0,%%xmm2,%%xmm2	\n\t	roundpd	$0,%%xmm10,%%xmm10		\n\t"\
					"subpd	%%xmm2,%%xmm1			\n\t	subpd	%%xmm10,%%xmm9			\n\t"\
					"mulpd	%%xmm4,%%xmm1			\n\t	mulpd	%%xmm12,%%xmm9			\n\t				movq	%%r11,0x330(%%rsi)	\n\t"\
					"/* Move high 2 words of result into input registers expected by start of loop body: */\n\t"\
					"movaps	%%xmm2,%%xmm4 /* fx2 */	\n\t	movaps	%%xmm10,%%xmm12	/* hx2 */	\n\t"\
					"movaps	%%xmm1,%%xmm2 /* fx1 */	\n\t	movaps	%%xmm9 ,%%xmm10	/* hx1 */	\n\t			movq	%%r15,0x338(%%rsi)	\n\t"\
					"subq	$1,%%rcx	/* j-- */		\n\t"\
					"cmpq	$0,%%rcx	/* j > 0 ?	*/	\n\t"\
					"jge	Loop8Beg	/* if (j >= 0), loop */	\n\t"\
				"Loop8End:							\n\t"\
					:					/* outputs: none */\
					: [__fq0] "m" (fq0)	/* All inputs from memory addresses here */\
					 ,[__pshift] "m" (pshift)	\
					 ,[__start_index] "m" (start_index)		\
					: "cc","memory","rax","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"\
					,"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* All user registers but rbx,xmm7 get clobbered. */\
				);

			__asm__ volatile (\
				"movq	%[__fx0],%%rax	\n\t"\
				"movaps	%%xmm0,    (%%rax)	\n\t	movaps	%%xmm8 ,0x10(%%rax)	\n\t"\
				"movaps	%%xmm2,0x20(%%rax)	\n\t	movaps	%%xmm10,0x30(%%rax)	\n\t"\
				"movaps	%%xmm4,0x40(%%rax)	\n\t	movaps	%%xmm12,0x50(%%rax)	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
			);

			/* Need to un-scale the fq-terms via remultiply by 2^26 for the postprocessing step below: */
			__asm__ volatile (\
				"movq	%[__fq0],%%rax						\n\t"\
				"movq	%[__two26f],%%rsi					\n\t"\
				"movaps	    (%%rsi),%%xmm8	/* two26f */	\n\t"\
				"movaps	    (%%rax),%%xmm0	/* fq0 */		\n\t	movaps	0x10(%%rax),%%xmm4	/* hq0 */		\n\t"\
				"movaps	0x20(%%rax),%%xmm1	/* fq1 */		\n\t	movaps	0x30(%%rax),%%xmm5	/* hq1 */		\n\t"\
				"movaps	0x40(%%rax),%%xmm2	/* fq2 */		\n\t	movaps	0x50(%%rax),%%xmm6	/* hq2 */		\n\t"\
				"mulpd	%%xmm8,%%xmm0						\n\t	mulpd	%%xmm8,%%xmm4						\n\t"\
				"mulpd	%%xmm8,%%xmm1						\n\t	mulpd	%%xmm8,%%xmm5						\n\t"\
				"mulpd	%%xmm8,%%xmm2						\n\t	mulpd	%%xmm8,%%xmm6						\n\t"\
				"movaps	%%xmm0,    (%%rax)	/* fq0/2^26 */	\n\t	movaps	%%xmm4,0x10(%%rax)	/* hq0/2^26 */	\n\t"\
				"movaps	%%xmm1,0x20(%%rax)	/* fq1/2^26 */	\n\t	movaps	%%xmm5,0x30(%%rax)	/* hq1/2^26 */	\n\t"\
				"movaps	%%xmm2,0x40(%%rax)	/* fq2/2^26 */	\n\t	movaps	%%xmm6,0x50(%%rax)	/* hq2/2^26 */	\n\t"\
				:					/* outputs: none */\
				: [__fq0] "m" (fq0)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm4","xmm5","xmm6","xmm8"		/* Clobbered registers */\
			);

			/*...Double and return.	These are specialized for the case
			where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
			*/
			CVT78_3WORD_DOUBLE_UINT96(*fx0,*fx1,*fx2, x0);
			CVT78_3WORD_DOUBLE_UINT96(*gx0,*gx1,*gx2, x1);
			CVT78_3WORD_DOUBLE_UINT96(*hx0,*hx1,*hx2, x2);
			CVT78_3WORD_DOUBLE_UINT96(*ix0,*ix1,*ix2, x3);

			/* For the 4 FP moduli, since we added k = (96-78) to 'true' value of pshift,
			need to multiply output x of the modular powering loop by 2^k (mod q).
			Since k is small (say small enougb to fit into a floating mantissa), use floating-divide to compute (x*2^k)/q:
			*/
		if(~pshift != p+78) {
			// Extra power of 2 is because in this flow we do not do the final 2*x-q step in the 'else' below:
			kmul0 = (double)((uint64)1 << (~pshift - (p+77)));
			kmul1 = kmul0;
			kmul2 = kmul0;
			kmul3 = kmul0;
			qmul0  = *fx0 + *fx1 * TWO26FLOAT;
			qmul1  = *gx0 + *gx1 * TWO26FLOAT;
			qmul2  = *hx0 + *hx1 * TWO26FLOAT;
			qmul3  = *ix0 + *ix1 * TWO26FLOAT;
			qmul0 += *fx2 * TWO26FLOAT * TWO26FLOAT;
			qmul1 += *gx2 * TWO26FLOAT * TWO26FLOAT;
			qmul2 += *hx2 * TWO26FLOAT * TWO26FLOAT;
			qmul3 += *ix2 * TWO26FLOAT * TWO26FLOAT;
			qmul0 *= kmul0;
			qmul1 *= kmul1;
			qmul2 *= kmul2;
			qmul3 *= kmul3;
			// This is the multiplicative inverse of q, not to be confused with the Montgomery-inverse:
			qdiv0  = *fq0 + *fq1 * TWO26FLOAT;
			qdiv1  = *gq0 + *gq1 * TWO26FLOAT;
			qdiv2  = *hq0 + *hq1 * TWO26FLOAT;
			qdiv3  = *iq0 + *iq1 * TWO26FLOAT;
			qdiv0 += *fq2 * TWO26FLOAT * TWO26FLOAT;
			qdiv1 += *gq2 * TWO26FLOAT * TWO26FLOAT;
			qdiv2 += *hq2 * TWO26FLOAT * TWO26FLOAT;
			qdiv3 += *iq2 * TWO26FLOAT * TWO26FLOAT;
			// Compute approx. ratio:
			qmul0 /= qdiv0;
			qmul1 /= qdiv1;
			qmul2 /= qdiv2;
			qmul3 /= qdiv3;
			qmul0 = DNINT(qmul0);
			qmul1 = DNINT(qmul1);
			qmul2 = DNINT(qmul2);
			qmul3 = DNINT(qmul3);
			// Use 96-bit variable to store q*qmul in preparation for exact-arithmetic subtract:
			MUL_LOHI64( x0.d0, (uint64)kmul0, x0.d0, tmp0);
			MUL_LOHI64( x1.d0, (uint64)kmul1, x1.d0, tmp1);
			MUL_LOHI64( x2.d0, (uint64)kmul2, x2.d0, tmp2);
			MUL_LOHI64( x3.d0, (uint64)kmul3, x3.d0, tmp3);
			x0.d1 = tmp0 + kmul0*x0.d1;
			x1.d1 = tmp1 + kmul1*x1.d1;
			x2.d1 = tmp2 + kmul2*x2.d1;
			x3.d1 = tmp3 + kmul3*x3.d1;
			MUL_LOHI64(q0.d0,(uint64)qmul0 ,lo0.d0, tmp0);
			MUL_LOHI64(q1.d0,(uint64)qmul1 ,lo1.d0, tmp1);
			MUL_LOHI64(q2.d0,(uint64)qmul2 ,lo2.d0, tmp2);
			MUL_LOHI64(q3.d0,(uint64)qmul3 ,lo3.d0, tmp3);
			lo0.d1 = tmp0 + qmul0*q0.d1;	lo0.d0 -= FERMAT;
			lo1.d1 = tmp1 + qmul1*q1.d1;	lo1.d0 -= FERMAT;
			lo2.d1 = tmp2 + qmul2*q2.d1;	lo2.d0 -= FERMAT;
			lo3.d1 = tmp3 + qmul3*q3.d1;	lo3.d0 -= FERMAT;
			SUB96(x0,lo0,x0);
			SUB96(x1,lo1,x1);
			SUB96(x2,lo2,x2);
			SUB96(x3,lo3,x3);
		} else {
			ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
			ADD96(x1,x1,x1);
			ADD96(x2,x2,x2);
			ADD96(x3,x3,x3);
			if(FERMAT) {
				q0.d0 -= 2;
				q1.d0 -= 2;
				q2.d0 -= 2;
				q3.d0 -= 2;
			}
			SUB96(x0,q0,x0);
			SUB96(x1,q1,x1);
			SUB96(x2,q2,x2);
			SUB96(x3,q3,x3);
		}

			ADD96_PTR(x4 ,x4, x4);
			ADD96_PTR(x5 ,x5, x5);
			ADD96_PTR(x6 ,x6, x6);
			ADD96_PTR(x7 ,x7, x7);

			SUB96_PTR(x4,(qptr4-FERMAT),x4);
			SUB96_PTR(x5,(qptr5-FERMAT),x5);
			SUB96_PTR(x6,(qptr6-FERMAT),x6);
			SUB96_PTR(x7,(qptr7-FERMAT),x7);

			tmp0 = CMPEQ96(x0, ONE96);
			tmp1 = CMPEQ96(x1, ONE96);
			tmp2 = CMPEQ96(x2, ONE96);
			tmp3 = CMPEQ96(x3, ONE96);
			tmp4 = CMPEQ96_PTR(x4, ONE96_PTR);
			tmp5 = CMPEQ96_PTR(x5, ONE96_PTR);
			tmp6 = CMPEQ96_PTR(x6, ONE96_PTR);
			tmp7 = CMPEQ96_PTR(x7, ONE96_PTR);

			r = tmp0;
			r += tmp1 << 1;
			r += tmp2 << 2;
			r += tmp3 << 3;
			r += tmp4 << 4;
			r += tmp5 << 5;
			r += tmp6 << 6;
			r += tmp7 << 7;
			return r;
		}

	  #endif

	  #ifdef USE_AVX	// Nov 2013: Now that have prototyped 2-TF-input/4-SSE-register modmul macro,
						// use 8-input all-float SSE2 code as basis for 16-input AVX and AVX2 modpow mode:

		uint64 twopmodq78_3WORD_DOUBLE_q16(uint64 p, uint64 k[], int init_sse2, int thr_id)
		{
		#if FAC_DEBUG
			int dbg  = (p==17715697 && k[0] = 67885719840ull);
		#endif
			const char func[] = "twopmodq78_3WORD_DOUBLE_q16";
			 int32 j = 0;	/* This needs to be signed because of the LR binary exponentiation. */
			uint32 q32[16], qinv32[16], tmp32;
			uint64 lead7, r = 0, tmp64;
			uint96 q[16], qinv[16], qhalf[16], x[16], tmp96;
			uint192 prod192;
			static uint64 psave = 0, pshift;
			static uint32 start_index, zshift, first_entry = TRUE;
			uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

			static int max_threads = 1;	// Default local-array-init is for just a single thread ... caller can re-init for > 1 threads later, if desired.
			static double *sc_arr = 0x0, *sc_ptr;
			double *dptr, dtmp;
			vec_dbl *tmp;
		  #ifdef MULTITHREAD
			static double *__r0;	// Base address for discrete per-thread local stores ...
									// *NOTE*: more convenient to use double* rather than vec_dbl* here
		  #else
			static	// Following set of pointers only static-izable in single-thread mode
		  #endif
			double *fq0[16],*fq1[16],*fq2[16],*fqhi52[16], *fqinv0[16],*fqinv1[16],*fqinv2[16], *fx0[16],*fx1[16],*fx2[16],
					*two13i, *two26f,*two26i, *two52f,*two52i;

		#if FAC_DEBUG
			if(dbg) printf("%s with p = %llu, k[] = %llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n",
						func,p,k[0x0],k[0x1],k[0x2],k[0x3],k[0x4],k[0x5],k[0x6],k[0x7],k[0x8],k[0x9],k[0xa],k[0xb],k[0xc],k[0xd],k[0xe],k[0xf]);
		#endif
			if(p != psave)
			{
			//	first_entry = FALSE;
				psave  = p;
				pshift = p + 78;
				j = leadz64(pshift);
				/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 77: */
				lead7 = ((pshift<<j) >> 57);	// leftmost bit = 1, thus lead7 in [64,127]
				if(lead7 > 77)
				{
					lead7 >>= 1;			// lead-6 bits in [39,63]
					start_index =  64-j-6;	/* Use only the leftmost 6 bits */
				}
				else
					start_index =  64-j-7;
				// lead7 in [39,77], thus zshift in [0,38]:
				zshift = 77 - lead7;
				zshift <<= 1;	/* Doubling the shift count here takes cares of the first SQR_LOHI */
								// Result in [0,76], i.e. qinv << (zshift<<1) always has at least the leading bit set.
				pshift = ~pshift;
			  #if FAC_DEBUG
				if(dbg) printf("pshift = 0x%llX\n",pshift);
			  #endif
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
				// Alloc the local-memory block:
				sc_arr = ALLOC_DOUBLE(sc_arr, 0xfc*max_threads + 4);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
				sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);	// Force vec_dbl-alignment
				ASSERT(((uintptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			#ifdef MULTITHREAD
				__r0  = sc_ptr;
				two13i = sc_ptr + 0xc0;
				two26f = sc_ptr + 0xc4;
				two26i = sc_ptr + 0xc8;
				two52f = sc_ptr + 0xcc;	// (only used for AVX2/FMA)
				two52i = sc_ptr + 0xd0;	// (only used for AVX2/FMA)
				for(j = 0; j < max_threads; ++j) {
					/* These remain fixed within each per-thread local store: */
					// Need specific-length inits here since might be building overall in AVX-512 mode:
					VEC_DBL_INIT_4((vec_dbl*)two13i, TWO13FLINV);
					VEC_DBL_INIT_4((vec_dbl*)two26f, TWO26FLOAT);
					VEC_DBL_INIT_4((vec_dbl*)two26i, TWO26FLINV);
					VEC_DBL_INIT_4((vec_dbl*)two52f, TWO52FLOAT);
					VEC_DBL_INIT_4((vec_dbl*)two52i, TWO52FLINV);
					// Move on to next thread's local store:
					two13i   += 0xfc;
					two26f   += 0xfc;
					two26i   += 0xfc;
					two52f   += 0xfc;
					two52i   += 0xfc;
				}
			#else
				/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 4 to span an AVX register: */
				fq0   [0] = sc_ptr + 0x00;
				fq1   [0] = sc_ptr + 0x10;
				fq2   [0] = sc_ptr + 0x20;
				fqhi52[0] = sc_ptr + 0x30;
				fqinv0[0] = sc_ptr + 0x40;
				fqinv1[0] = sc_ptr + 0x50;
				fqinv2[0] = sc_ptr + 0x60;
				fx0   [0] = sc_ptr + 0x70;
				fx1   [0] = sc_ptr + 0x80;
				fx2   [0] = sc_ptr + 0x90;
				for(j = 1, dptr = sc_ptr+1; j < 16; j++, dptr++) {
					fq0   [j] = dptr + 0x00;
					fq1   [j] = dptr + 0x10;
					fq2   [j] = dptr + 0x20;
					fqhi52[j] = dptr + 0x30;
					fqinv0[j] = dptr + 0x40;
					fqinv1[j] = dptr + 0x50;
					fqinv2[j] = dptr + 0x60;
					fx0   [j] = dptr + 0x70;
					fx1   [j] = dptr + 0x80;
					fx2   [j] = dptr + 0x90;
				}
			// +0xa0,b0 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
			// Consts: only get 1 AVX-reg (4 doubles) per:
				two13i = sc_ptr + 0xc0;
				two26f = sc_ptr + 0xc4;
				two26i = sc_ptr + 0xc8;
				two52f = sc_ptr + 0xcc;	// (only used for AVX2/FMA)
				two52i = sc_ptr + 0xd0;	// (only used for AVX2/FMA)
			// Total: Equivalent of 12*4 + 3 = 53 AVX slots, alloc 56
				/* Can premultiply each of the multiword-mul inputs by 1/sqrt(2^13) due to quadraticity: */
				// Need specific-length inits here since might be building overall in AVX-512 mode:
				VEC_DBL_INIT_4((vec_dbl*)two13i, TWO13FLINV);
				VEC_DBL_INIT_4((vec_dbl*)two26f, TWO26FLOAT);
				VEC_DBL_INIT_4((vec_dbl*)two26i, TWO26FLINV);
				VEC_DBL_INIT_4((vec_dbl*)two52f, TWO52FLOAT);
				VEC_DBL_INIT_4((vec_dbl*)two52i, TWO52FLINV);
			#endif
				if(init_sse2) return 0;
			}	/* end of inits */

			/* If multithreaded, set the local-store pointers needed for the current thread; */
		#ifdef MULTITHREAD
			ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
			sc_ptr = __r0 + thr_id*0xfc;
			/* Remember, these are POINTERS-TO-DOUBLES, so need an increment of 4 to span an AVX register: */
			fq0   [0] = sc_ptr + 0x00;
			fq1   [0] = sc_ptr + 0x10;
			fq2   [0] = sc_ptr + 0x20;
			fqhi52[0] = sc_ptr + 0x30;
			fqinv0[0] = sc_ptr + 0x40;
			fqinv1[0] = sc_ptr + 0x50;
			fqinv2[0] = sc_ptr + 0x60;
			fx0   [0] = sc_ptr + 0x70;
			fx1   [0] = sc_ptr + 0x80;
			fx2   [0] = sc_ptr + 0x90;
			for(j = 1, dptr = sc_ptr+1; j < 16; j++, dptr++) {
				fq0   [j] = dptr + 0x00;
				fq1   [j] = dptr + 0x10;
				fq2   [j] = dptr + 0x20;
				fqhi52[j] = dptr + 0x30;
				fqinv0[j] = dptr + 0x40;
				fqinv1[j] = dptr + 0x50;
				fqinv2[j] = dptr + 0x60;
				fx0   [j] = dptr + 0x70;
				fx1   [j] = dptr + 0x80;
				fx2   [j] = dptr + 0x90;
			}
		// +0xa0,b0 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
		// Consts: only get 1 AVX-reg (4 doubles) per:
			two13i = sc_ptr + 0xc0;
			two26f = sc_ptr + 0xc4;
			two26i = sc_ptr + 0xc8;
		/*** In Non-AVX2 mode, Byte-offsets in 8*[0xcc,0xfc) = [0x660,0x7e0), i.e. 12 YMM-sized slots,
		= ax0 + 8*[0x5c,0x8c) = ax0 + [0x2e0,0x460) = (ax0 + 0x2e0) + [0,0x180)  available for debug ***/
			two52f = sc_ptr + 0xcc;	// (only used for AVX2/FMA)
			two52i = sc_ptr + 0xd0;	// (only used for AVX2/FMA)
		/*** In AVX2 mode, Byte-offsets in 8*[0xd4,0xfc) = [0x6a0,0x7e0), i.e. 10 YMM-sized slots,
		= ax0 + 8*[0x64,0x8c) = ax0 + [0x320,0x460) = (ax0 + 0x320) + [0,0x140)  available for debug ***/
		#endif

		#ifdef MUL_LOHI64_SUBROUTINE
			#error MUL_LOHI64_SUBROUTINE defined!
		#endif
			ASSERT((p >> 63) == 0, "twopmodq78_q16: p must be < 2^63!");
			for(j = 0; j < 16; j++)
			{
				q[j].d0 = p+p;
				MUL_LOHI64(q[j].d0, k[j], q[j].d0, q[j].d1);	q[j].d0 += 1;	// Since 2*p*k even, no need to check for overflow here
				q32[j] = (uint32)q[j].d0;
				CVT_UINT78_3WORD_DOUBLE(q[j] ,*fq0[j],*fq1[j],*fq2[j]);	// Convert q to floating form
				tmp96.d0 = q[j].d0; tmp96.d1 = q[j].d1;	// Stash-in-scalar-uint96 to work around "expected expression" compiler error
				RSHIFT_FAST96(tmp96, 1, tmp96);	// qhalf = (q-1)/2, since q odd.
				qhalf[j].d0 = tmp96.d0; qhalf[j].d1 = tmp96.d1;
				// Iterative mod-inverse computation:
				/* 8 bits: */
				qinv32[j] = minv8[(q32[j] & 0xff)>>1];
				/* 16 bits: */
				tmp32 = q32[j]*qinv32[j];	qinv32[j] = qinv32[j]*((uint32)2 - tmp32);
				/* 32 bits: */
				tmp32 = q32[j]*qinv32[j];	qinv32[j] = qinv32[j]*((uint32)2 - tmp32);
				/* 64 bits: */
				qinv[j].d0 = (uint64)qinv32[j];
				tmp64 = q[j].d0*qinv[j].d0;
				qinv[j].d0 = qinv[j].d0*((uint64)2 - tmp64);
				/* 96 bits, truncated to 78: */
				MULH64(q[j].d0, qinv[j].d0, tmp64);
				qinv[j].d1 = -qinv[j].d0*(q[j].d1*qinv[j].d0 + tmp64);
				qinv[j].d1 &= 0x3fff;
				/* Convert 78-bit qinv to floating form: */
				CVT_UINT78_3WORD_DOUBLE(qinv[j], *fqinv0[j],*fqinv1[j],*fqinv2[j]);
			}

			dptr = fq0[0];
			__asm__ volatile (\
				"movq	%[__fq0],%%rax						\n\t"\
				"movq	%[__two26f],%%rsi					\n\t"\
				"vmovaps	    (%%rsi),%%ymm8	/* two26f */	\n\t"\
				"vmovaps	0x20(%%rsi),%%ymm9	/* two26i */	\n\t"\
				"vmovaps	     (%%rax),%%ymm0	/* fq0[0-3] */	\n\t	vmovaps	0x020(%%rax),%%ymm4	/* fq0[4-7] */	\n\t"\
				"vmovaps	0x080(%%rax),%%ymm1	/* fq1[0-3] */	\n\t	vmovaps	0x0a0(%%rax),%%ymm5	/* fq1[4-7] */	\n\t"\
				"vmovaps	0x100(%%rax),%%ymm2	/* fq2[0-3] */	\n\t	vmovaps	0x120(%%rax),%%ymm6	/* fq2[4-7] */	\n\t"\
				"vmovaps	%%ymm2,%%ymm3		/* cpy fq2 */	\n\t	vmovaps	%%ymm6,%%ymm7		/* cpy fq2 */	\n\t"\
				"vmulpd	%%ymm8,%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm8,%%ymm7,%%ymm7				\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3/* Hi 52 out bits */\n\t	vaddpd	%%ymm5,%%ymm7,%%ymm7	/* Hi 52*/	\n\t"\
				"vmulpd	%%ymm9,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm9,%%ymm4,%%ymm4				\n\t"\
				"vmulpd	%%ymm9,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm9,%%ymm5,%%ymm5				\n\t"\
				"vmulpd	%%ymm9,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm9,%%ymm6,%%ymm6				\n\t"\
				"vmovaps	%%ymm0,     (%%rax)	/* fq0/2^26 */	\n\t	vmovaps	%%ymm4,0x020(%%rax)	/* fq0/2^26 */	\n\t"\
				"vmovaps	%%ymm1,0x080(%%rax)	/* fq1/2^26 */	\n\t	vmovaps	%%ymm5,0x0a0(%%rax)	/* fq1/2^26 */	\n\t"\
				"vmovaps	%%ymm2,0x100(%%rax)	/* fq2/2^26 */	\n\t	vmovaps	%%ymm6,0x120(%%rax)	/* fq2/2^26 */	\n\t"\
				"vmovaps	%%ymm3,0x180(%%rax)	/* fq[hi52] */	\n\t	vmovaps	%%ymm7,0x1a0(%%rax)	/* fq[hi52] */	\n\t"\
				/* Now do 2nd set of doubles: */\
				"vmovaps	0x040(%%rax),%%ymm0	/* fq0[8-11] */	\n\t	vmovaps	0x060(%%rax),%%ymm4	/* fq0[12-15] */\n\t"\
				"vmovaps	0x0c0(%%rax),%%ymm1	/* fq1[8-11] */	\n\t	vmovaps	0x0e0(%%rax),%%ymm5	/* fq1[12-15] */\n\t"\
				"vmovaps	0x140(%%rax),%%ymm2	/* fq2[8-11] */	\n\t	vmovaps	0x160(%%rax),%%ymm6	/* fq2[12-15] */\n\t"\
				"vmovaps	%%ymm2,%%ymm3		/* cpy fq2 */	\n\t	vmovaps	%%ymm6,%%ymm7		/* cpy fq2 */	\n\t"\
				"vmulpd	%%ymm8,%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm8,%%ymm7,%%ymm7				\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3/* Hi 52 out bits */\n\t	vaddpd	%%ymm5,%%ymm7,%%ymm7	/* Hi 52*/	\n\t"\
				"vmulpd	%%ymm9,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm9,%%ymm4,%%ymm4				\n\t"\
				"vmulpd	%%ymm9,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm9,%%ymm5,%%ymm5				\n\t"\
				"vmulpd	%%ymm9,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm9,%%ymm6,%%ymm6				\n\t"\
				"vmovaps	%%ymm0,0x040(%%rax)	/* fq0/2^26 */	\n\t	vmovaps	%%ymm4,0x060(%%rax)	/* fq0/2^26 */	\n\t"\
				"vmovaps	%%ymm1,0x0c0(%%rax)	/* fq1/2^26 */	\n\t	vmovaps	%%ymm5,0x0e0(%%rax)	/* fq1/2^26 */	\n\t"\
				"vmovaps	%%ymm2,0x140(%%rax)	/* fq2/2^26 */	\n\t	vmovaps	%%ymm6,0x160(%%rax)	/* fq2/2^26 */	\n\t"\
				"vmovaps	%%ymm3,0x1c0(%%rax)	/* fq[hi52] */	\n\t	vmovaps	%%ymm7,0x1e0(%%rax)	/* fq[hi52] */	\n\t"\
				:					/* outputs: none */\
				: [__fq0] "m" (dptr)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
			);
		  #if FAC_DEBUG
			if(dbg) printf("fq = %20.15f, %20.15f, %20.15f\n", *fq0[0] * TWO26FLOAT, *fq1[0] * TWO26FLOAT, *fq2[0] * TWO26FLOAT);
		  #endif

			dptr = fqinv0[0];
		  #ifdef USE_AVX2
			// In AVX2+FMA version, overwrite middle term of usual qinv[0,1,2] sequence with low 52 bits of qinv:
			/* Moved inside Iter-1-left-shift-MULL-emulation below, since needs to come after those qinv-loads */
		  #else
			__asm__ volatile (\
				"movq	%[__fqinv0],%%rax	\n\t"\
				"movq	%[__two26i],%%rsi	\n\t"\
				"vmovaps	(%%rsi),%%ymm12		\n\t"/* TWO26FLINV */\
				"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	0x020(%%rax),%%ymm3	\n\t	vmovaps	0x040(%%rax),%%ymm6	\n\t	vmovaps	0x060(%%rax),%%ymm9	\n\t"\
				"vmovaps	0x080(%%rax),%%ymm1	\n\t	vmovaps	0x0a0(%%rax),%%ymm4	\n\t	vmovaps	0x0c0(%%rax),%%ymm7	\n\t	vmovaps	0x0e0(%%rax),%%ymm10	\n\t"\
				"vmovaps	0x100(%%rax),%%ymm2	\n\t	vmovaps	0x120(%%rax),%%ymm5	\n\t	vmovaps	0x140(%%rax),%%ymm8	\n\t	vmovaps	0x160(%%rax),%%ymm11	\n\t"\
				"vmulpd	%%ymm12,%%ymm0,%%ymm0	\n\t vmulpd	%%ymm12,%%ymm3,%%ymm3	\n\t vmulpd	%%ymm12,%%ymm6,%%ymm6	\n\t vmulpd	%%ymm12,%%ymm9 ,%%ymm9		\n\t"\
				"vmulpd	%%ymm12,%%ymm1,%%ymm1	\n\t vmulpd	%%ymm12,%%ymm4,%%ymm4	\n\t vmulpd	%%ymm12,%%ymm7,%%ymm7	\n\t vmulpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
				"vmulpd	%%ymm12,%%ymm2,%%ymm2	\n\t vmulpd	%%ymm12,%%ymm5,%%ymm5	\n\t vmulpd	%%ymm12,%%ymm8,%%ymm8	\n\t vmulpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
				"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm3,0x020(%%rax)	\n\t	vmovaps	%%ymm6,0x040(%%rax)	\n\t	vmovaps	%%ymm9 ,0x060(%%rax)	\n\t"\
				"vmovaps	%%ymm1,0x080(%%rax)	\n\t	vmovaps	%%ymm4,0x0a0(%%rax)	\n\t	vmovaps	%%ymm7,0x0c0(%%rax)	\n\t	vmovaps	%%ymm10,0x0e0(%%rax)	\n\t"\
				"vmovaps	%%ymm2,0x100(%%rax)	\n\t	vmovaps	%%ymm5,0x120(%%rax)	\n\t	vmovaps	%%ymm8,0x140(%%rax)	\n\t	vmovaps	%%ymm11,0x160(%%rax)	\n\t"\
				:					/* outputs: none */\
				: [__fqinv0] "m" (dptr)	/* All inputs from memory addresses here */\
				 ,[__two26i] "m" (two26i)	\
				: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12"	/* Clobbered registers */\
			);
		   #if FAC_DEBUG
			if(dbg) printf("fqinv[0,1,2,lo52] = %20.15f, %20.15f, %20.15f, %20.15f\n", *fqinv0[0] * TWO26FLOAT, *fqinv1[0] * TWO26FLOAT, *fqinv2[0] * TWO26FLOAT, (*fqinv0[0] + *fqinv1[0] * TWO26FLOAT) * TWO26FLOAT);
		   #endif
		  #endif

		#ifdef USE_AVX2
			// zshift is in [0,76] - The SIMD-double-based shift code below needs a shift count < 30 (roughly),
			// so fiddle our base-2^26 data here to be able to use zshift (mod 26):
			dptr = fqinv0[0];
			if(zshift < 26) {
				dtmp = 1<<zshift;
				__asm__ volatile (\
					"movq	%[__fqinv0],%%rax	\n\t"\
					"movq	%[__two26f],%%rsi	\n\t"\
					"vmovaps	(%%rsi),%%ymm3	\n\t"/* TWO26FLOAT */\
					"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	0x020(%%rax),%%ymm4		\n\t	vmovaps	0x040(%%rax),%%ymm8		\n\t	vmovaps	0x060(%%rax),%%ymm12	\n\t"\
					"vmovaps	0x080(%%rax),%%ymm1		\n\t	vmovaps	0x0a0(%%rax),%%ymm5		\n\t	vmovaps	0x0c0(%%rax),%%ymm9		\n\t	vmovaps	0x0e0(%%rax),%%ymm13	\n\t"\
					"vmovaps	%%ymm1,%%ymm2			\n\t	vmovaps	%%ymm5,%%ymm6			\n\t	vmovaps	%%ymm9,%%ymm10			\n\t	vmovaps	%%ymm13,%%ymm14			\n\t"/* cpy fhi into flo-destination reg */\
				"vfmadd132pd %%ymm3,%%ymm0,%%ymm2		\n\t vfmadd132pd %%ymm3,%%ymm4,%%ymm6	\n\t vfmadd132pd %%ymm3,%%ymm8,%%ymm10	\n\t vfmadd132pd %%ymm3,%%ymm12,%%ymm14	\n\t"/* (qinv0 + 2^26*qinv1)*/\
					"vmovaps	%%ymm2,0x080(%%rax)		\n\t	vmovaps	%%ymm6,0x0a0(%%rax)		\n\t	vmovaps	%%ymm10,0x0c0(%%rax)	\n\t	vmovaps	%%ymm14,0x0e0(%%rax)	\n\t"\
					"vmovaps	0x100(%%rax),%%ymm2		\n\t	vmovaps	0x120(%%rax),%%ymm6		\n\t	vmovaps	0x140(%%rax),%%ymm10	\n\t	vmovaps	0x160(%%rax),%%ymm14	\n\t"\
					:					/* outputs: none */\
					: [__fqinv0] "m" (dptr)	/* All inputs from memory addresses here */\
					 ,[__two26f] "m" (two26f)	\
					: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6", "xmm8","xmm9","xmm10", "xmm12","xmm13","xmm14"	/* Clobbered registers */\
				);
			} else if(zshift < 52) {
				dtmp = 1<<(zshift-26);
				// qinv <<= 26:
				__asm__ volatile (\
					"movq	%[__fqinv0],%%rax	\n\t"\
					"movq	%[__two26f],%%rsi	\n\t"\
					"vmovaps	(%%rsi),%%ymm3	\n\t"/* TWO26FLOAT */\
					"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	0x020(%%rax),%%ymm4		\n\t	vmovaps	0x040(%%rax),%%ymm8		\n\t	vmovaps	0x060(%%rax),%%ymm12	\n\t"\
					"vmovaps	0x080(%%rax),%%ymm1		\n\t	vmovaps	0x0a0(%%rax),%%ymm5		\n\t	vmovaps	0x0c0(%%rax),%%ymm9		\n\t	vmovaps	0x0e0(%%rax),%%ymm13	\n\t"\
					"vmovaps	%%ymm1,%%ymm2			\n\t	vmovaps	%%ymm5,%%ymm6			\n\t	vmovaps	%%ymm9,%%ymm10			\n\t	vmovaps	%%ymm13,%%ymm14			\n\t"/* cpy fhi into flo-destination reg */\
				"vfmadd132pd %%ymm3,%%ymm0,%%ymm2		\n\t vfmadd132pd %%ymm3,%%ymm4,%%ymm6	\n\t vfmadd132pd %%ymm3,%%ymm8,%%ymm10	\n\t vfmadd132pd %%ymm3,%%ymm12,%%ymm14	\n\t"/* (qinv0 + 2^26*qinv1)*/\
					"vmovaps	%%ymm2,0x080(%%rax)		\n\t	vmovaps	%%ymm6,0x0a0(%%rax)		\n\t	vmovaps	%%ymm10,0x0c0(%%rax)	\n\t	vmovaps	%%ymm14,0x0e0(%%rax)	\n\t"\
					/* Now shift all data leftward 1 word, losing the topmost digit, then zero the low word: */\
					"vmovaps	%%ymm1,%%ymm2			\n\t	vmovaps	%%ymm5,%%ymm6			\n\t	vmovaps	%%ymm9,%%ymm10			\n\t	vmovaps	%%ymm13,%%ymm14			\n\t"\
					"vmovaps	%%ymm0,%%ymm1			\n\t	vmovaps	%%ymm4,%%ymm5			\n\t	vmovaps	%%ymm8,%%ymm9			\n\t	vmovaps	%%ymm12,%%ymm13			\n\t"\
					"vxorpd		%%ymm0,%%ymm0,%%ymm0	\n\t	vxorpd	%%ymm4,%%ymm4,%%ymm4	\n\t	vxorpd	%%ymm8,%%ymm8,%%ymm8	\n\t	vxorpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
					:					/* outputs: none */\
					: [__fqinv0] "m" (dptr)	/* All inputs from memory addresses here */\
					 ,[__two26f] "m" (two26f)	\
					: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6", "xmm8","xmm9","xmm10", "xmm12","xmm13","xmm14"	/* Clobbered registers */\
				);
			} else if(zshift < 78) {
				dtmp = 1<<(zshift-52);
				// qinv <<= 52:
				__asm__ volatile (\
					"movq	%[__fqinv0],%%rax	\n\t"\
					"movq	%[__two26f],%%rsi	\n\t"\
					"vmovaps	(%%rsi),%%ymm3	\n\t"/* TWO26FLOAT */\
					"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	0x020(%%rax),%%ymm4		\n\t	vmovaps	0x040(%%rax),%%ymm8		\n\t	vmovaps	0x060(%%rax),%%ymm12	\n\t"\
					"vmovaps	0x080(%%rax),%%ymm1		\n\t	vmovaps	0x0a0(%%rax),%%ymm5		\n\t	vmovaps	0x0c0(%%rax),%%ymm9		\n\t	vmovaps	0x0e0(%%rax),%%ymm13	\n\t"\
					"vmovaps	%%ymm1,%%ymm2			\n\t	vmovaps	%%ymm5,%%ymm6			\n\t	vmovaps	%%ymm9,%%ymm10			\n\t	vmovaps	%%ymm13,%%ymm14			\n\t"/* cpy fhi into flo-destination reg */\
				"vfmadd132pd %%ymm3,%%ymm0,%%ymm2		\n\t vfmadd132pd %%ymm3,%%ymm4,%%ymm6	\n\t vfmadd132pd %%ymm3,%%ymm8,%%ymm10	\n\t vfmadd132pd %%ymm3,%%ymm12,%%ymm14	\n\t"/* (qinv0 + 2^26*qinv1)*/\
					"vmovaps	%%ymm2,0x080(%%rax)		\n\t	vmovaps	%%ymm6,0x0a0(%%rax)		\n\t	vmovaps	%%ymm10,0x0c0(%%rax)	\n\t	vmovaps	%%ymm14,0x0e0(%%rax)	\n\t"\
					/* Now shift all data leftward 2 word, losing the topmost 2 digits, then zero the low 2 words: */\
					"vmovaps	%%ymm0,%%ymm2			\n\t	vmovaps	%%ymm4,%%ymm6			\n\t	vmovaps	%%ymm8,%%ymm10			\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
					"vxorpd		%%ymm0,%%ymm0,%%ymm0	\n\t	vxorpd	%%ymm4,%%ymm4,%%ymm4	\n\t	vxorpd	%%ymm8,%%ymm8,%%ymm8	\n\t	vxorpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
					"vxorpd		%%ymm1,%%ymm1,%%ymm1	\n\t	vxorpd	%%ymm5,%%ymm5,%%ymm5	\n\t	vxorpd	%%ymm9,%%ymm9,%%ymm9	\n\t	vxorpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
					:					/* outputs: none */\
					: [__fqinv0] "m" (dptr)	/* All inputs from memory addresses here */\
					 ,[__two26f] "m" (two26f)	\
					: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6", "xmm8","xmm9","xmm10", "xmm12","xmm13","xmm14"	/* Clobbered registers */\
				);
			} else {
				ASSERT(0,"zshift out of range!");
			}
		//	VEC_DBL_INIT_4((vec_dbl*)fx0[0], dtmp);	This wound up using zmm0 in some of my builds, so moved this init inside the iter1 macro below

			// Since zshift is a power of two < 2^77, use a streamlined code sequence for the first iteration,
			// in which x0^2 is just 1<<zshift, i.e. we compute x^2; lo = MULL(qinv, x^2) simply as qinv<<zshift (discard bits beyond lowest 78).
			// NOTE: The ensuing first-iteration MULH(q,lo) is prone to delicate low-bits cancellation, i.e. we cannot us the low-product-omitting
			// version of that as e.g. in our FMA-based 78-bit modmul macros. Thus, base our streamlined first-iteration SIMD macro below on the
			// non-FMA version of SSE2_twopmodq78_modmul_q16():
			dptr = &dtmp;	j = start_index-1;
			SSE2_twopmodq78_modmul_q16_iter1(fq0[0],fqinv0[0],fx0[0],two26f,two26i,pshift,j,dptr);

		#else
			/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
			tmp32 = ((pshift >> (start_index-1)) & (uint32)1);	// This is same for all candidates q's
			for(j = 0; j < 16; j++) {
				/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
				LSHIFT96(qinv[j], zshift, tmp96);		tmp96.d1 &= 0x3fff;
				// Inline a (... >> 78) specialization of the RSHIFT192 macro here, since have a constant shift count:
				MUL_LOHI96_PROD192(q[j],tmp96,prod192);	tmp96.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	tmp96.d1 = (prod192.d2 >> 14);
				/* hi = 0 in this instance, which simplifies things. */
				/* Put the result in lo (rather than x), to ease overflow check below */
				SUB96(q[j], tmp96, x[j]);
				if(tmp32) {
					/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
					if(CMPUGT96(x[j], qhalf[j])){ ADD96(x[j], x[j], x[j]); SUB96(x[j], q[j], x[j]); }else{ ADD96(x[j], x[j], x[j]); }
				}
				/* Convert x to floating form: */
				CVT_UINT78_3WORD_DOUBLE(x[j] ,*fx0[j],*fx1[j],*fx2[j]);
			}
			// Load starting data into the SIMD registers:
			dptr = fx0[0];
			__asm__ volatile (\
				"movq	%[__fx0],%%rax	\n\t"\
				"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	0x020(%%rax),%%ymm4	\n\t	vmovaps	0x040(%%rax),%%ymm8		\n\t	vmovaps	0x060(%%rax),%%ymm12	\n\t"\
				"vmovaps	0x080(%%rax),%%ymm1	\n\t	vmovaps	0x0a0(%%rax),%%ymm5	\n\t	vmovaps	0x0c0(%%rax),%%ymm9		\n\t	vmovaps	0x0e0(%%rax),%%ymm13	\n\t"\
				"vmovaps	0x100(%%rax),%%ymm2	\n\t	vmovaps	0x120(%%rax),%%ymm6	\n\t	vmovaps	0x140(%%rax),%%ymm10	\n\t	vmovaps	0x160(%%rax),%%ymm14	\n\t"\
				:					/* outputs: none */\
				: [__fx0] "m" (dptr)	/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm12","xmm13","xmm14"	/* Clobbered registers */\
			);
		#endif

		#if FAC_DEBUG
		  if(dbg) {
			printf("p = %llu, k0 = %llu, start_index0 = %u, initial shift = %u\n",p,k[0],start_index,zshift);
			printf("On modpow-loop entry: start_index = %u,\n\tfx0-2 = %20.15f, %20.15f, %20.15f, %20.15f\n",start_index, *fx0[0],*fx1[0],*fx2[0]);
		  }
		#endif

			/*...x^2 mod q is returned in x. */
			/* All 3-word-double-form operands have components in the following size ranges:
				fword0,1 in [-2^25, +2^25]
				fword2   in [   -1, +2^26]
			*/

		/* Inner loop body needs 42 MOVAPS, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */
		/* !ASM_LOOP: High-level loop construct, and explicit branch to do the modular doubling: */
			for(j = start_index-2; j >= 0; j--)
			{
				SSE2_twopmodq78_modmul_q16(fq0[0],fqinv0[0],fx0[0],two26i,pshift,j);
			#if FAC_DEBUG
				// Sep 2016: Examination of .s-file (-c ==> -S -fverbose-asm) shows *ax0-2 impl below using xmm0-2:
				#warning ******** activating this inside-loop print may clobber one or more of the above xmm-data ********
				if(dbg) {
					printf("J = %d, 2x = %d:\n",j,((pshift >> j) & (uint64)1) == 1);
				#if 1
					dptr = fx0[0];
					printf("\tx0-2 = %20.5f, %20.5f, %20.5f\n",*dptr,*(dptr+0x10),*(dptr+0x20));
		// Print contents of debug-slots:
		dptr = fx0[0] + 0x64;
		printf("\tdbg-slots = %20.5f, %20.5f, %20.5f\n",*dptr,*(dptr+4),*(dptr+8),*(dptr+12));
				#elif 0
					// Print contents of fx[0]-... slots:
					dptr = fx0[0];
					printf("\th.lo26,hi52 = %20.5f, %20.5f\n",*(dptr+0x30),*(dptr+0x40));
					printf("\tl.lo26,hi52 = %20.5f, %20.5f\n",*dptr,*(dptr+0x10));
					printf("\tl.CY/flo2 ?? = %20.5f\n",*(dptr+0x20));
				#else
					// Print contents of debug-slots:
					dptr = fx0[0] + 0x80;
					printf("\tfx0-2 = %20.5f, %20.5f, %20.5f\n",*dptr,*(dptr+4),*(dptr+8));
				#endif
					exit(0);
				}
			#endif
			}	/* for(j...) */

			/* Do a final mod-doubling and return. Code here is specialized for the case
			where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2 or -(q-1)/2:
			*/
		#ifdef USE_AVX2

			dptr = fq0[0];
			__asm__ volatile (\
				"movq	%[__fq0]   ,%%rax	\n\t"\
				"movq	%[__two26f],%%rbx	\n\t"\
				"vmovaps	0x20(%%rbx),%%ymm3		\n\t	vmovaps	(%%rbx),%%ymm7				\n\t"/* two26i,f */\
			"vfmadd213pd %%ymm1,%%ymm7,%%ymm2		\n\t vfmadd213pd %%ymm5,%%ymm7,%%ymm6		\n\t vfmadd213pd %%ymm9,%%ymm7,%%ymm10			\n\t vfmadd213pd %%ymm13,%%ymm7,%%ymm14		\n\t"/* (q1 + 2^26*qinv2) = hi52 */\
				/* Mod-doubling code taken from SSE2_twopmodq78_modmul_q16, expects x.lo26,hi52 in ymm1,2, resp., so copy ymm0,4,8,12 -> ymm1,5,9,13 after computing x.hi52 (in ymm2,6,10,14): */\
				"vaddpd	%%ymm0,%%ymm0,%%ymm1		\n\t	vaddpd	%%ymm4,%%ymm4,%%ymm5		\n\t	vaddpd	%%ymm8,%%ymm8,%%ymm9			\n\t vaddpd	%%ymm12,%%ymm12,%%ymm13			\n\t"/* low 26 bits */\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10			\n\t vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"/* top 52 bits */\
				/* And subtract q.lo26,hi52 - note q.lo26 is balanced-digit normalized and divided by 2^26, so do this sub prior to final >0 normalization below: */\
			"vfnmadd231pd	 (%%rax),%%ymm7,%%ymm1	\n\t vfnmadd231pd 0x20(%%rax),%%ymm7,%%ymm5	\n\t vfnmadd231pd 0x40(%%rax),%%ymm7,%%ymm9		\n\t vfnmadd231pd 0x60(%%rax),%%ymm7,%%ymm13\n\t"/* -= qlo26 */\
				"vsubpd	0x180(%%rax),%%ymm2,%%ymm2	\n\t	vsubpd	0x1a0(%%rax),%%ymm6,%%ymm6	\n\t	vsubpd	0x1c0(%%rax),%%ymm10,%%ymm10	\n\t	vsubpd	0x1e0(%%rax),%%ymm14,%%ymm14\n\t"/* -= qhi52 */\
				/* Final-residue computation needs digits >= 0 normalization - since fq0 will get folded in with fq2 anyway, only need renormalize fq0: */\
				"vmulpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vmulpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vmulpd	%%ymm3,%%ymm9,%%ymm9			\n\t vmulpd	%%ymm3,%%ymm13,%%ymm13			\n\t"/* xlo *= two26i */\
				"vroundpd $1,%%ymm1,%%ymm0			\n\t	vroundpd $1,%%ymm5,%%ymm4			\n\t	vroundpd $1,%%ymm9,%%ymm8				\n\t vroundpd $1,%%ymm13,%%ymm12			\n\t"/* fcy */\
				"vsubpd	%%ymm0,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm4,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm8,%%ymm9,%%ymm9			\n\t vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"/* fx0 -= fcy */\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm8,%%ymm10,%%ymm10			\n\t vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"/* Add carry into xhi */\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1		\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vmulpd	%%ymm7,%%ymm9,%%ymm9			\n\t vmulpd	%%ymm7 ,%%ymm13,%%ymm13			\n\t"/* fx0 *= two26f ... result must go into m0 */\
				/* (lo26 == 1 && hi52 == 0), save result in a bitmask: */
				"vmovaps	(%%rbx),%%ymm0			\n\t"/* TWO26FLOAT */\
				"vmulpd	0x20(%%rbx),%%ymm0,%%ymm3	\n\t"/* 1.0 ... no good way to generate 1.0 on-the-fly, so use 2^26 * 2^-26 */\
				"vxorpd	     %%ymm0,%%ymm0,%%ymm0	\n\t"/* 0.0 */\
				"vcmppd $0,%%ymm1,%%ymm3,%%ymm1		\n\t	vcmppd $0,%%ymm5,%%ymm3,%%ymm5		\n\t	vcmppd $0,%%ymm9,%%ymm3,%%ymm9			\n\t vcmppd $0,%%ymm13,%%ymm3,%%ymm13		\n\t"/* low 26 bits == 1 ? */\
				"vcmppd $0,%%ymm2,%%ymm0,%%ymm2		\n\t	vcmppd $0,%%ymm6,%%ymm0,%%ymm6		\n\t	vcmppd $0,%%ymm10,%%ymm0,%%ymm10		\n\t vcmppd $0,%%ymm14,%%ymm0,%%ymm14		\n\t"/* top 52 bits == 0 ? */\
				"vandpd	%%ymm2,%%ymm1,%%ymm0		\n\t	vandpd	%%ymm6,%%ymm5,%%ymm4		\n\t	vandpd	%%ymm10,%%ymm9,%%ymm8			\n\t	vandpd	%%ymm14,%%ymm13,%%ymm12		\n\t"/* and the 2 results together */\
				"vmovmskpd	%%ymm0,%%rax			\n\t	vmovmskpd	%%ymm4,%%rbx			\n\t	vmovmskpd	%%ymm8,%%rcx				\n\t	vmovmskpd	%%ymm12,%%rdx			\n\t"/* extract result into four 4-bit bitmasks */\
				/* concatenate into 16-bit result and write: */\
				"shlq	$4 ,%%rbx		\n\t"\
				"shlq	$8 ,%%rcx		\n\t"\
				"shlq	$12,%%rdx		\n\t"\
				"addq	%%rbx,%%rax		\n\t"\
				"addq	%%rdx,%%rcx		\n\t"\
				"addq	%%rcx,%%rax		\n\t"\
				"movq	%%rax,%[__result]	\n\t"\
				:					/* outputs: none */\
				: [__fq0] "m" (dptr)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)	\
				 ,[__result] "m" (r)	\
				: "cc","memory","cl","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
			);

		#else	// AVX code uses legacy convert-back-to-int code:

			dptr = fx0[0];
			__asm__ volatile (\
				"movq	%[__fx0],%%rax	\n\t"\
				"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm4,0x020(%%rax)	\n\t	vmovaps	%%ymm8 ,0x040(%%rax)	\n\t	vmovaps	%%ymm12,0x060(%%rax)	\n\t"\
				"vmovaps	%%ymm1,0x080(%%rax)	\n\t	vmovaps	%%ymm5,0x0a0(%%rax)	\n\t	vmovaps	%%ymm9 ,0x0c0(%%rax)	\n\t	vmovaps	%%ymm13,0x0e0(%%rax)	\n\t"\
				"vmovaps	%%ymm2,0x100(%%rax)	\n\t	vmovaps	%%ymm6,0x120(%%rax)	\n\t	vmovaps	%%ymm10,0x140(%%rax)	\n\t	vmovaps	%%ymm14,0x160(%%rax)	\n\t"\
				:					/* outputs: none */\
				: [__fx0] "m" (dptr)	/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm1","xmm2","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm12","xmm13","xmm14"	/* Clobbered registers */\
			);

		  #if FAC_DEBUG
			if(dbg) printf("fx = %20.15f, %20.15f, %20.15f\n", *fx0[0], *fx1[0], *fx2[0]);
		  #endif
			for(j = 0; j < 16; j++) {
				CVT78_3WORD_DOUBLE_UINT96(*fx0[j],*fx1[j],*fx2[j], x[j]);
				// In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow:
				ADD96(x[j],x[j],x[j]);
				q[j].d0 -= FERMAT;	// += 2 of result (achieve by q -= 2, and then sub-q from x) if Fermat number
				SUB96(x[j],q[j],x[j]);
				tmp64 = CMPEQ96(x[j], ONE96);
			#if FAC_DEBUG
				if(dbg) printf("xout[0] = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x[j])]);
			#endif
				r += (tmp64 << j);
			}

		#endif	// AVX or AVX2 ?

		#if FAC_DEBUG
			if(dbg) {
				printf("xout_q16 = %llX\n",r);
				exit(0);
			}
		#endif
			return r;
		}

	  #endif	// USE_AVX ?

	  #ifdef USE_AVX512	// Oct 2016: Use 16-input float-FMA AVX2 code as basis for 32-input AVX-512 modpow mode.
						// Basic implementation - requiring just the AVX-512 Foundation instructions - is just a
						// doubled-width version of above _q16 routine, which I initially deployed for Haswell, but
						// starting with Cannonlake (expected release 2H2017) will get a major speedup via the new
						// AVX-512 IFMA (integer fused multiply-add) instruction extensions, which will likely share
						// their wide-MUL hardware with the floating-point FMA, but whose VPMADD52[L|H]UQ instructions
						// will allow us to do the same 104-bit FMA-based product arithmetic is in the _q16 code,
						// but with half as many MULs, no ADDs or ROUNDs, and of course at twice the vector width. Use
						// of the latter extensions will be triggered in my code via the USE_AVX512_I preprocessor flag.

		uint64 twopmodq78_3WORD_DOUBLE_q32(uint64 p, uint64 k[], int init_sse2, int thr_id)
		{//	p = 17732269; k[0] = k[1] = k[2] = k[3] = 133912345560ull;	// Debug
			const char func[] = "twopmodq78_3WORD_DOUBLE_q32";
			 int32 j = 0;	/* This needs to be signed because of the LR binary exponentiation. */
			uint64 lead7, r = 0, tmp64;
			static uint64 psave = 0, pshift;
			static uint32 start_index, zshift, first_entry = TRUE;
			uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case
			uint8* minv8_ptr = minv8;	// Ptr to Table of precomputed byte-inverses def'd in mi64.h
			static int max_threads = 1;	// Default local-array-init is for just a single thread ... caller can re-init for > 1 threads later, if desired.
		#ifdef USE_AVX512_I
		  #error AVX-512 IFMA instruction extensions not yet supported!
			static uint64 *sc_arr = 0x0, *sc_ptr;
			uint64 *itmp;
			vec_u64 *tmp;
		  #ifdef MULTITHREAD
			static uint64 *__r0;	// Base address for discrete per-thread local stores ...
									// *NOTE*: more convenient to use uint64* rather than vec_u64* here
		  #else
			static	// Following set of pointers only static-izable in single-thread mode
		  #endif
			uint64 *fq0[32],*fq1[32],*fq2[32],*fqhi52[32], *fqinv0[32],*fqinv1[32],*fqinv2[32], *fx0[32],*fx1[32],*fx2[32],
					*mask_lo26,*mask_lo52;
			for(j = 0; j < 32; j++) {
				ASSERT((k[j] >> 52) == 0ull, "Ks must be < 2^52!");
			}
		#else
			static double *sc_arr = 0x0, *sc_ptr;
			double *dptr, dtmp;
			vec_dbl *tmp,*tm1;
		  #ifdef MULTITHREAD
			static double *__r0;	// Base address for discrete per-thread local stores ...
		  #else
			static	// Following set of pointers only static-izable in single-thread mode
		  #endif
			double *fq0[32],*fq1[32],*fq2[32],*fqhi52[32], *fqinv0[32],*fqinv1[32],*fqinv2[32], *fx0[32],*fx1[32],*fx2[32],
					*two13i, *two26f,*two26i, *two52f,*two52i,
					kdbl[32];
			// AVX-512 Foundation lacks the needed DQ extensions, so use HLL to convert kvec entries to double:
			for(j = 0; j < 32; j++) {
				ASSERT((k[j] >> 52) == 0ull, "Ks must be < 2^52!");
				kdbl[j] = (double)k[j];
			}
		#endif

			if(p != psave)
			{
			//	first_entry = FALSE;
				psave  = p;
				pshift = p + 78;
				j = leadz64(pshift);
				/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 77: */
				lead7 = ((pshift<<j) >> 57);
				if(lead7 > 77)
				{
					lead7 >>= 1;
					start_index =  64-j-6;	/* Use only the leftmost 6 bits */
				}
				else
					start_index =  64-j-7;

				zshift = 77 - lead7;
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

				sc_arr = ALLOC_UINT64(sc_arr, 0x1c0*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
				sc_ptr = (uint64 *)ALIGN_VEC_U64(sc_arr);	// Force vec_u64-alignment
				ASSERT(((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			  #ifdef MULTITHREAD
				__r0  = sc_ptr;
				mask_lo26 = sc_ptr + 0x180;
				mask_lo52 = sc_ptr + 0x188;
				for(j = 0; j < max_threads; ++j) {
					/* These remain fixed within each per-thread local store: */
					// Need specific-length inits here since might be building overall in AVX-512 mode:
					VEC_U64_INIT((vec_u64*)mask_lo26, 0x0000000003FFFFFFull);
					VEC_U64_INIT((vec_u64*)mask_lo52, 0x000FFFFFFFFFFFFFull);
					// Move on to next thread's local store:
					mask_lo26 += 0x1c0;
					mask_lo52 += 0x1c0;
				}
			  #else
				/* Remember, these are POINTERS-TO-UINT64, so need an increment of 8 to span an AVX-512 register: */
				fq0   [0] = sc_ptr + 0x000;
				fq1   [0] = sc_ptr + 0x020;
				fq2   [0] = sc_ptr + 0x040;
				fqhi52[0] = sc_ptr + 0x060;
				fqinv0[0] = sc_ptr + 0x080;
				fqinv1[0] = sc_ptr + 0x0A0;
				fqinv2[0] = sc_ptr + 0x0C0;
				fx0   [0] = sc_ptr + 0x0E0;
				fx1   [0] = sc_ptr + 0x100;
				fx2   [0] = sc_ptr + 0x120;
				for(j = 1, itmp = sc_ptr+1; j < 32; j++, itmp++) {
					fq0   [j] = itmp + 0x000;
					fq1   [j] = itmp + 0x020;
					fq2   [j] = itmp + 0x040;
					fqhi52[j] = itmp + 0x060;
					fqinv0[j] = itmp + 0x080;
					fqinv1[j] = itmp + 0x0A0;
					fqinv2[j] = itmp + 0x0C0;
					fx0   [j] = itmp + 0x0E0;
					fx1   [j] = itmp + 0x100;
					fx2   [j] = itmp + 0x120;
				}
			// +0x140,160 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
			// Consts: only get 1 AVX-512 reg (8 uint64s) per:
				mask_lo26 = sc_ptr + 0x180;
				mask_lo52 = sc_ptr + 0x188;
			// Total: Equivalent of 12*4 + 2 = 50 AVX-512 slots, alloc 56
				VEC_U64_INIT((vec_u64*)mask_lo26, 0x0000000003FFFFFFull);
				VEC_U64_INIT((vec_u64*)mask_lo52, 0x000FFFFFFFFFFFFFull);
			  #endif

			/***************************************************/
			#else	// Default AVX-512 floating-point-FMA mode
			/***************************************************/

				sc_arr = ALLOC_DOUBLE(sc_arr, 0x1c0*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
				sc_ptr = (uint64 *)ALIGN_VEC_DBL(sc_arr);	// Force vec_u64-alignment
				ASSERT(((uintptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			  #ifdef MULTITHREAD
				__r0  = sc_ptr;
				two13i = sc_ptr + 0x180;
				two26f = sc_ptr + 0x188;
				two26i = sc_ptr + 0x190;
				two52f = sc_ptr + 0x198;
				two52i = sc_ptr + 0x1a0;
				for(j = 0; j < max_threads; ++j) {
					/* These remain fixed within each per-thread local store: */
					// Need specific-length inits here since might be building overall in AVX-512 mode:
					VEC_DBL_INIT((vec_dbl*)two13i, TWO13FLINV);
					VEC_DBL_INIT((vec_dbl*)two26f, TWO26FLOAT);
					VEC_DBL_INIT((vec_dbl*)two26i, TWO26FLINV);
					VEC_DBL_INIT((vec_dbl*)two52f, TWO52FLOAT);
					VEC_DBL_INIT((vec_dbl*)two52i, TWO52FLINV);
					// Move on to next thread's local store:
					two13i   += 0x1c0;
					two26f   += 0x1c0;
					two26i   += 0x1c0;
					two52f   += 0x1c0;
					two52i   += 0x1c0;
				}
			  #else
				/* Remember, these are POINTERS-TO-DOUBLE, so need an increment of 8 to span an AVX-512 register: */
				fq0   [0] = sc_ptr + 0x000;
				fq1   [0] = sc_ptr + 0x020;
				fq2   [0] = sc_ptr + 0x040;
				fqhi52[0] = sc_ptr + 0x060;
				fqinv0[0] = sc_ptr + 0x080;
				fqinv1[0] = sc_ptr + 0x0A0;
				fqinv2[0] = sc_ptr + 0x0C0;
				fx0   [0] = sc_ptr + 0x0E0;
				fx1   [0] = sc_ptr + 0x100;
				fx2   [0] = sc_ptr + 0x120;
				for(j = 1, dptr = sc_ptr+1; j < 32; j++, dptr++) {
					fq0   [j] = dptr + 0x000;
					fq1   [j] = dptr + 0x020;
					fq2   [j] = dptr + 0x040;
					fqhi52[j] = dptr + 0x060;
					fqinv0[j] = dptr + 0x080;
					fqinv1[j] = dptr + 0x0A0;
					fqinv2[j] = dptr + 0x0C0;
					fx0   [j] = dptr + 0x0E0;
					fx1   [j] = dptr + 0x100;
					fx2   [j] = dptr + 0x120;
				}
			// +0x140,160 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
			// Consts: only get 1 AVX-512 reg (8 doubles) per:
				two13i = sc_ptr + 0x180;
				two26f = sc_ptr + 0x188;
				two26i = sc_ptr + 0x190;
				two52f = sc_ptr + 0x198;
				two52i = sc_ptr + 0x1a0;
			// Total: Equivalent of 12*4 + 5 = 53 AVX-512 slots, alloc 56
				VEC_DBL_INIT((vec_dbl*)two13i, TWO13FLINV);
				VEC_DBL_INIT((vec_dbl*)two26f, TWO26FLOAT);
				VEC_DBL_INIT((vec_dbl*)two26i, TWO26FLINV);
				VEC_DBL_INIT((vec_dbl*)two52f, TWO52FLOAT);
				VEC_DBL_INIT((vec_dbl*)two52i, TWO52FLINV);
			  #endif

			#endif
				if(init_sse2) return 0;
			}	/* end of inits */

			/* If multithreaded, set the local-store pointers needed for the current thread; */
		#ifdef MULTITHREAD

			ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
			sc_ptr = __r0 + thr_id*0x1c0;

		  #ifdef USE_AVX512_I

			/* Remember, these are POINTERS-TO-UINT64, so need an increment of 8 to span an AVX-512 register: */
			fq0   [0] = sc_ptr + 0x000;
			fq1   [0] = sc_ptr + 0x020;
			fq2   [0] = sc_ptr + 0x040;
			fqhi52[0] = sc_ptr + 0x060;
			fqinv0[0] = sc_ptr + 0x080;
			fqinv1[0] = sc_ptr + 0x0A0;
			fqinv2[0] = sc_ptr + 0x0C0;
			fx0   [0] = sc_ptr + 0x0E0;
			fx1   [0] = sc_ptr + 0x100;
			fx2   [0] = sc_ptr + 0x120;
			for(j = 1, itmp = sc_ptr+1; j < 32; j++, itmp++) {
				fq0   [j] = itmp + 0x000;
				fq1   [j] = itmp + 0x020;
				fq2   [j] = itmp + 0x040;
				fqhi52[j] = itmp + 0x060;
				fqinv0[j] = itmp + 0x080;
				fqinv1[j] = itmp + 0x0A0;
				fqinv2[j] = itmp + 0x0C0;
				fx0   [j] = itmp + 0x0E0;
				fx1   [j] = itmp + 0x100;
				fx2   [j] = itmp + 0x120;
			}
			mask_lo26 = sc_ptr + 0x180;
			mask_lo52 = sc_ptr + 0x188;
			VEC_U64_INIT((vec_u64*)mask_lo26, 0x0000000003FFFFFFull);
			VEC_U64_INIT((vec_u64*)mask_lo52, 0x000FFFFFFFFFFFFFull);

		  #else

			/* Remember, these are POINTERS-TO-DOUBLE, so need an increment of 8 to span an AVX-512 register: */
			fq0   [0] = sc_ptr + 0x000;
			fq1   [0] = sc_ptr + 0x020;
			fq2   [0] = sc_ptr + 0x040;
			fqhi52[0] = sc_ptr + 0x060;
			fqinv0[0] = sc_ptr + 0x080;
			fqinv1[0] = sc_ptr + 0x0A0;
			fqinv2[0] = sc_ptr + 0x0C0;
			fx0   [0] = sc_ptr + 0x0E0;
			fx1   [0] = sc_ptr + 0x100;
			fx2   [0] = sc_ptr + 0x120;
			for(j = 1, dptr = sc_ptr+1; j < 32; j++, dptr++) {
				fq0   [j] = dptr + 0x000;
				fq1   [j] = dptr + 0x020;
				fq2   [j] = dptr + 0x040;
				fqhi52[j] = dptr + 0x060;
				fqinv0[j] = dptr + 0x080;
				fqinv1[j] = dptr + 0x0A0;
				fqinv2[j] = dptr + 0x0C0;
				fx0   [j] = dptr + 0x0E0;
				fx1   [j] = dptr + 0x100;
				fx2   [j] = dptr + 0x120;
			}
			two13i = sc_ptr + 0x180;
			two26f = sc_ptr + 0x188;
			two26i = sc_ptr + 0x190;
			two52f = sc_ptr + 0x198;
			two52i = sc_ptr + 0x1a0;
			VEC_DBL_INIT((vec_dbl*)two13i, TWO13FLINV);
			VEC_DBL_INIT((vec_dbl*)two26f, TWO26FLOAT);
			VEC_DBL_INIT((vec_dbl*)two26i, TWO26FLINV);
			VEC_DBL_INIT((vec_dbl*)two52f, TWO52FLOAT);
			VEC_DBL_INIT((vec_dbl*)two52i, TWO52FLINV);

		  #endif

		#endif

		#ifdef MUL_LOHI64_SUBROUTINE
			#error MUL_LOHI64_SUBROUTINE defined!
		#endif
			ASSERT((p >> 63) == 0, "twopmodq78_q32: p must be < 2^63!");

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
	#if 0
	W.r.to the vector mod-inverse computation, consider our options here:
	o Prefer pure-integer math here in order to leverage the properties of 2s-complement integer math.
	  That means first computing q = 2.p.k+1 (restricting p,k each to 52 bits for now) using floating-FMA,
	  normalizing w.r.to base-2^26 balanced-digit and storing those into the fq local-store, then converting
	  to base-2^23 unsigned-integer form for the ensuing mod-inverse computation.

	o Since (at least on KNL) mod-inverses will be stored in base-2^26 balanced-digit flaoting-point form,
	  desire the entire mod-inverse computation to be entirely in-register, to avoid having to set aside
	  special integer local-memory storage for intermediates.

	The second column [*] indicates whether the instruction in question is part of the AVX-512
	Foundation subset (i.e. available on KNL) or some future-deployment subset, as noted.
	The third column [||] indicates the instruction`s degree of parallelism with respect to
	full-register-width AVX-512.

	Instruction	*	||	Description
	VPMULLD		F	16	Multiply the packed dword signed integers in zmm2 and zmm3/m512/m32bcst
						and store the low 32 bits of each product in zmm1 under writemask k1.
	VPMULLQ		DQ	8	Multiply the packed dword signed integers in zmm2 and zmm3/m512/m32bcst
						and store the low 64 bits of each product in zmm1 under writemask k1.
	VPMULDQ		F	8	Multiply packed signed doubleword integers in zmm2 by packed signed doubleword integers
						in zmm3/m512/m64bcst, and store the quadword results in zmm1 using writemask k1.
	VPMULHUW	BW	32	Multiply the packed unsigned word integers in zmm2 and zmm3/m512, and store the
						high 16 bits of the results in zmm1 under writemask k1.
						***NOTE:*** Only available in 'W' form, i.e. for 16-bit 'word' integer inputs.
									256-bit version available under AVX2; 512-bit needs AVX-512 BW.
	VPMULUDQ	F		Multiply packed unsigned doubleword integers in zmm2 by packed unsigned doubleword ints
						in zmm3/m512/m64bcst, and store the quadword results in zmm1 under writemask k1.

	Upshot: Since this preprocessing code is executed just once per parallel-modpow call, no point trying to
	get fancy and maximize ||-bandwidth at each stage; instead just keep the double --> dword-int convrsion
	data in qword-width register subfields and operate on those. Here are iteration-specific notes, which we
	illustrate using a 78-bit example, the largest known factor of MM31,
		q = 178021379228511215367151
		= 9650<<64 + 10298917214042272751
		= 39528686<<52 + 54962424<<26 + 35992559 :

	0. Vector-load of need qinv 8-bit inits from the minv8[] byte-table def`d in mi64.h. This gather-load must needs
	   use the VPGATHERDD instruction, but no need to mask-off the high 3 bytes of each dword datum gathered
	   in the resulting vector because of the self-correcting nature of the ensuing qinv-computing iteration:

			qinv8 = minv8[q8>>1] = minv8[239>>1] = minv8[119] = 15

	1. q-data are 32-bit (of which we need the low 16, but just use q32 and VPMULLD), qinv are 8-bit, compute

			tmp16  = MULL16(q16,qinv8);				 2817
			qinv16 = MULL16(qinv8,(2-tmp16));		23311; Verify: MULL16(q16,qinv16) = 1

	2. q-data are 32-bit (3794088943), qinv-data are 16-bit (23311), compute

			tmp32  = MULL32(q32,qinv16);			2040791041
			qinv32 = MULL32(qinv16,(2-tmp32));		2472827663; Verify: MULL32(q32,qinv32) = 1

	3. Now things begin to get interesting - q-data are 64-bit (10298917214042272751), qinv-data are 32-bit.
	   Here is the 32x32 ==> 64-bit qinv-updating code from my twopmodq96 routines, with the corresponding
	   required vector instructions noted in [] the right-column comments:

			tmp64  = MULL64(q64,qinv32);			 3115510555725529089	Needs VPMULUDQ(q64.lo32,qinv32) + VPMULLD(q64.hi32,qinv32)<<32
			qinv64 = MULL64(qinv32,(2-tmp64));		13405235700914477839	Needs VPMULUDQ((2-tmp64).lo32,qinv32) + VPMULLD((2-tmp64).hi32,qinv32)<<32
											= 3121149656<<32 + 2472827663	Verify: MULL64(q64,qinv64) = 1

										Q: qinv is 32-bit in RHS here - can we truncate (2 - tmp64) to 32-bit and again use 32x32 ==> 64-bit MUL?
										A: No, because by construction tmp64 == 1 (mod 2^32), i.e. that would give (2 - tmp) = 1 and leave qinv
										unchanged at its 32-bit input value. That is, qinv.d0*((uint64)2 - tmp64) needs to be a genuine MULL64.
	4. 64 => 96-bits: q-data are 96-bit (in our integer-qinv-compute-then-truncate-to-78-bit and cvt-to-double model), qinv-data are 64-bit:

		tmp64  = MULH64(q.lo64, qinv);				 7484215760070392893
		tmp64 += MULL64(q.hi32, qinv);				  992540659696056235 (Note add is mod-2^64)
		qinv.hi32 = MULL64(-qinv, tmp64);			 5580892910975415291 (this gives 128-bits-good qinv; truncate as desired.)
	96-bit truncation of qinv = 2687464443<<64 + 3121149656<<32 + 2472827663;
	78-bit truncation of qinv =      13307<<64 + 3121149656<<32 + 2472827663 = 54508448<<52 + 37598756<<26 + 56908559 .
	Based on the above arithmetic structure, we should be able to use 32-bit MULs, check it:

		tmp32  = MULH64(q.lo64 * qinv)%2^32;		2841656381
		tmp32 += MULL32(q.hi32, qinv.lo32);			2790307755 (Note add is mod-2^32)
		qinv.hi32 = MULL32(-qinv.lo32, tmp32);		2687464443 (this = 5580892910975415291.lo32 and gives 96-bits-good qinv; truncate as desired.)

	The only one of the resulting truncated-operand MULs needing restructuring to exploit the capabilities of AVX vec-int arithmetic is the first;
	let`s break it down into 32-bit operand sub-multiplies and see how the resulting bits line up:

		(q.lo64 * qinv) = (2^32*q.md32 + q.lo32) * (2^32*qinv.md32 + qinv.lo32)

	Bitranges:																Subproducts:
										|-----------  63: 0 -----------|	q.lo32 * qinv.lo32
						|-----------  95:32 -----------|					q.md32 * qinv.lo32, q.lo32 * qinv.md32
		|----------- 127:64 -----------|									q.md32 * qinv.md32
						|--- OUTPUT ---|

	We need to compute all 4 subproducts because of the carryins to the higher-order terms, but only need 32x32 ==> 64-bit MULs for the lower-order 3.
	Here is the needed vector-operation sequence:

		A. VPMULUDQ(q.lo32, qinv.lo32), right-shift 64-bit result 32 bits to keep just the high halves;
		B. VPMULUDQ(q.md32, qinv.lo32), VPADDQ to result of [A], can ignore any carry because it would be into bit 96;
		C. VPMULUDQ(q.lo32, qinv.md32), VPADDQ to result of [B], again can ignore any carry, right-shift 64-bit result 32 bits;
		D. VPMULLD (q.md32, qinv.md32), VPADDQ to result of [C], mask off high 32 bits of sum.
	which gives the initial tmp32 ... and then tmp32 += MULL32(q.hi32, qinv.lo32) and feed result into qinv.hi32 = MULL32(-qinv.lo32, tmp32).
	=========================
	78-bit q = 178021379228511215367151 = 3794088943 + 2397903523<<32 + 9650<<64: lo26 = 35992559, hi52 = 2652725267835128 = 54962424 + 39528686<<26
	52-bit q =          295257526626031 = 4294829807 +      68744<<32           : lo26 = 66971375, hi52 = 4399679; both have bottom-byte = 239 ==> qinv8 = minv8[239>>1] = minv8[119] = 15
	=========================
	#endif
			dtmp = p+p;	dptr = &dtmp;	// dtmp = double; dptr = double*
			tmp = (vec_dbl*)fq0[0];	tm1 = (vec_dbl*)kdbl;	// Use vec_dbl* ptrs to hold address args ... asm-code cares not a whit about types
			__asm__ volatile (\
				"movq	%[__fq0] ,%%rax				\n\t"\
				"movq	%[__p2]  ,%%rbx				\n\t"\
				"movq	%[__kvec],%%rcx				\n\t"/* Each column uses 7 distinct vector registers, so increment all non-shared vec-reg indices by +7 for each additional rightward column: */\
				"vbroadcastsd	(%%rbx),%%zmm1		\n\t	vmovaps	%%zmm1,%%zmm8				\n\t	vmovaps	%%zmm1,%%zmm15				\n\t	vmovaps	%%zmm1,%%zmm22				\n\t"/* broadcast 2*p to all 8 qwords of zmm ... double-prec version is from mem-address */\
				"vmovupd	(%%rcx),%%zmm3			\n\t	vmovupd	0x40(%%rcx),%%zmm10			\n\t	vmovupd	0x80(%%rcx),%%zmm17			\n\t	vmovupd	0xc0(%%rcx),%%zmm24			\n\t"/* load the factor-candidate k's, 8 per zmm ... note kvec not nec. 16-byte aligned */\
				"movq	%[__two26f],%%rbx			\n\t"\
				"vmulpd		%%zmm1,%%zmm3,%%zmm2	\n\t	vmulpd		%%zmm8,%%zmm10,%%zmm9	\n\t	vmulpd		%%zmm15,%%zmm17,%%zmm16	\n\t	vmulpd		%%zmm22,%%zmm24,%%zmm23	\n\t"/* hi = x*y = 2p*k */\
				"vmovaps	%%zmm2,%%zmm0			\n\t	vmovaps	%%zmm9,%%zmm7				\n\t	vmovaps	%%zmm16,%%zmm14				\n\t	vmovaps	%%zmm23,%%zmm21				\n\t"/* cpy hi into lo-destination reg */\
			"vfmsub231pd	%%zmm1,%%zmm3,%%zmm0	\n\t vfmsub231pd	%%zmm8,%%zmm10,%%zmm7	\n\t vfmsub231pd	%%zmm15,%%zmm17,%%zmm14	\n\t vfmsub231pd	%%zmm22,%%zmm24,%%zmm21	\n\t"/* lo = fma(x,y,-hi; xmm1,3 FREE */\
				"movl	$1,%%esi					\n\t"\
				"vpbroadcastd	%%esi,%%zmm3		\n\t	vpbroadcastd	%%esi,%%zmm10		\n\t	vpbroadcastd	%%esi,%%zmm17		\n\t	vpbroadcastd	%%esi,%%zmm24		\n\t"/* Broadcast to all 8 int32 slots of zmm3, since KNL lacks AVX-512 VL extensions */\
				"vcvtdq2pd	%%ymm3,%%zmm3			\n\t	vcvtdq2pd	%%ymm10,%%zmm10			\n\t	vcvtdq2pd	%%ymm17,%%zmm17			\n\t	vcvtdq2pd	%%ymm24,%%zmm24			\n\t"/* zmm3 = (double)1, conversion from 4-way int32 in lower half of zmm3 */\
				"vaddpd		%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd		%%zmm10,%%zmm7,%%zmm7	\n\t	vaddpd		%%zmm17,%%zmm14,%%zmm14	\n\t	vaddpd		%%zmm24,%%zmm21,%%zmm21	\n\t"/* lo += 1 */\
			/* Normalize to properly split the resultant significant bits across lo,hi, use base-2^52 signed-int-stored-as-double: */\
				"vmovaps	%%zmm2,%%zmm1			\n\t	vmovaps	%%zmm9,%%zmm8				\n\t	vmovaps	%%zmm16,%%zmm15				\n\t	vmovaps	%%zmm23,%%zmm22				\n\t"/* cpy hi into hh-destination reg */\
				"vmulpd	  0xc0(%%rbx),%%zmm1,%%zmm1	\n\t	vmulpd	  0xc0(%%rbx),%%zmm8,%%zmm8	\n\t	vmulpd	0xc0(%%rbx),%%zmm15,%%zmm15	\n\t	vmulpd	0xc0(%%rbx),%%zmm22,%%zmm22	\n\t"/*            hi*TWO52FLINV  */\
				"vrndscalepd	$1,%%zmm1,%%zmm1	\n\t	vrndscalepd	$1,%%zmm8,%%zmm8		\n\t	vrndscalepd	$1,%%zmm15,%%zmm15		\n\t	vrndscalepd	$1,%%zmm22,%%zmm22		\n\t"/* hh = FLOOR(hi*TWO52FLINV) */\
				"vmovaps	%%zmm2,%%zmm3			\n\t	vmovaps	%%zmm9,%%zmm10				\n\t	vmovaps	%%zmm16,%%zmm17				\n\t	vmovaps	%%zmm23,%%zmm24				\n\t"/* cpy hi into cy-destination reg */\
			"vfnmadd231pd 0x80(%%rbx),%%zmm1,%%zmm3	\n\t vfnmadd231pd 0x80(%%rbx),%%zmm8,%%zmm10\n\tvfnmadd231pd 0x80(%%rbx),%%zmm15,%%zmm17\n\tvfnmadd231pd 0x80(%%rbx),%%zmm22,%%zmm24\n\t"/* cy = fma(hh,-TWO52FLOAT,hi)	Backward carry from hi into lo */\
				"vaddpd		%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd		%%zmm10,%%zmm7,%%zmm7	\n\t	vaddpd		%%zmm17,%%zmm14,%%zmm14	\n\t	vaddpd		%%zmm24,%%zmm21,%%zmm21	\n\t"/* lo += cy */\
				"vmovaps	%%zmm1,0x200(%%rax)		\n\t	vmovaps	%%zmm8,0x240(%%rax)			\n\t	vmovaps	%%zmm15,0x280(%%rax)		\n\t	vmovaps	%%zmm22,0x2c0(%%rax)		\n\t"/* store hi26 into fq2 */\
			/* Normalize to convert 78-bit q from balanced-digit [lo52,hi26] representation to nonnegative-digit [lo26,hi52]: */\
				"vmovaps	%%zmm0,%%zmm2			\n\t	vmovaps	%%zmm7,%%zmm9				\n\t	vmovaps	%%zmm14,%%zmm16				\n\t	vmovaps	%%zmm21,%%zmm23				\n\t"/* cpy lo into destination reg */\
				"vmulpd	  0x40(%%rbx),%%zmm2,%%zmm2	\n\t	vmulpd	  0x40(%%rbx),%%zmm9,%%zmm9	\n\t	vmulpd	0x40(%%rbx),%%zmm16,%%zmm16	\n\t	vmulpd	0x40(%%rbx),%%zmm23,%%zmm23	\n\t"/*            lo*TWO26FLINV  */\
				"vrndscalepd	$1,%%zmm2,%%zmm2	\n\t	vrndscalepd	$1,%%zmm9,%%zmm9		\n\t	vrndscalepd	$1,%%zmm16,%%zmm16		\n\t	vrndscalepd	$1,%%zmm23,%%zmm23		\n\t"/* cy = FLOOR(lo*TWO26FLINV)	Forward carry from lo into hi */\
				"vmovaps	%%zmm2,0x100(%%rax)		\n\t	vmovaps	%%zmm9,0x140(%%rax)			\n\t	vmovaps	%%zmm16,0x180(%%rax)		\n\t	vmovaps	%%zmm23,0x1c0(%%rax)		\n\t"/* store cy into fq1 */\
			"vfnmadd231pd     (%%rbx),%%zmm2,%%zmm0	\n\t vfnmadd231pd    (%%rbx),%%zmm9,%%zmm7	\n\tvfnmadd231pd    (%%rbx),%%zmm16,%%zmm14	\n\tvfnmadd231pd    (%%rbx),%%zmm23,%%zmm21	\n\t"/* lo26 = fma(cy,-TWO26FLOAT,lo) */\
			" vfmadd132pd     (%%rbx),%%zmm2,%%zmm1	\n\t vfmadd132pd     (%%rbx),%%zmm9,%%zmm8	\n\t vfmadd132pd    (%%rbx),%%zmm16,%%zmm15	\n\t vfmadd132pd    (%%rbx),%%zmm23,%%zmm22	\n\t"/* hi52 = fma(hi,+TWO26FLOAT,cy) */\
				"vmovaps	%%zmm0,0x000(%%rax)		\n\t	vmovaps	%%zmm7,0x040(%%rax)			\n\t	vmovaps	%%zmm14,0x080(%%rax)		\n\t	vmovaps	%%zmm21,0x0c0(%%rax)		\n\t"/* store lo26 into fq0 */\
				"vmovaps	%%zmm1,0x300(%%rax)		\n\t	vmovaps	%%zmm8,0x340(%%rax)			\n\t	vmovaps	%%zmm15,0x380(%%rax)		\n\t	vmovaps	%%zmm22,0x3c0(%%rax)		\n\t"/* store hi52 into fqhi52 */\
			/*** hi = hh ==> lo26 in zmm0, hi52 = zmm1. ***/\
			/* Assemble 32-bit pieces of q from the 26-bit floating-point ones ... AVX-512F has just double->int32, so
			convert our three 26-bit double-vec chunks to 8-way int32s (i.e. zmm -> ymm), then re-expand result to 8-way
			int64s (with high 32-bit halves all zero) in preparation for mod-inverse computation: */\
				"vcvttpd2dq		  %%zmm0,%%ymm0		\n\t	vcvttpd2dq		  %%zmm7,%%ymm7		\n\t	vcvttpd2dq		  %%zmm14,%%ymm14	\n\t	vcvttpd2dq		  %%zmm21,%%ymm21	\n\t"/* (int32)fq0 */\
				"vcvttpd2dq	0x100(%%rax),%%ymm1		\n\t	vcvttpd2dq	0x140(%%rax),%%ymm8		\n\t	vcvttpd2dq	0x180(%%rax),%%ymm15	\n\t	vcvttpd2dq	0x1c0(%%rax),%%ymm22	\n\t"/* (int32)fq1 */\
				"vcvttpd2dq	0x200(%%rax),%%ymm2		\n\t	vcvttpd2dq	0x240(%%rax),%%ymm9		\n\t	vcvttpd2dq	0x280(%%rax),%%ymm16	\n\t	vcvttpd2dq	0x2c0(%%rax),%%ymm23	\n\t"/* (int32)fq2 */\
		/*** Must use full-length zmm-regs here (even though upper 256 bits unused) due to AVX-512F full-width constraint (applying to ymm needs AVX-512VL extensions): ***/\
				"vpslld	$26,%%zmm1,%%zmm3			\n\t	vpslld	$26,%%zmm8,%%zmm10			\n\t	vpslld	$26,%%zmm15,%%zmm17			\n\t	vpslld	$26,%%zmm22,%%zmm24			\n\t"/* left-justify lo6 bits of fq1... */\
			"vpaddd	%%zmm3,%%zmm0,%%zmm0			\n\t	vpaddd	%%zmm10,%%zmm7,%%zmm7		\n\t	vpaddd	%%zmm17,%%zmm14,%%zmm14		\n\t	vpaddd	%%zmm24,%%zmm21,%%zmm21		\n\t"/* ...and add to fq0 to get q.lo32 */\
				"vpsrld	$06,%%zmm1,%%zmm1			\n\t	vpsrld	$06,%%zmm8,%%zmm8			\n\t	vpsrld	$06,%%zmm15,%%zmm15			\n\t	vpsrld	$06,%%zmm22,%%zmm22			\n\t"/* fq1>>6 leaves lo20 bits of eventual q.md32... */\
				"vpslld	$20,%%zmm2,%%zmm3			\n\t	vpslld	$20,%%zmm9,%%zmm10			\n\t	vpslld	$20,%%zmm16,%%zmm17			\n\t	vpslld	$20,%%zmm23,%%zmm24			\n\t"/* ...left-justify lo12 bits of fq2... */\
			"vpaddd	%%zmm3,%%zmm1,%%zmm1			\n\t	vpaddd	%%zmm10,%%zmm8,%%zmm8		\n\t	vpaddd	%%zmm17,%%zmm15,%%zmm15		\n\t	vpaddd	%%zmm24,%%zmm22,%%zmm22		\n\t"/* ...and add to get q.md32 */\
			"vpsrld	$12   ,%%zmm2,%%zmm2			\n\t	vpsrld	$12   ,%%zmm9,%%zmm9		\n\t	vpsrld	$12   ,%%zmm16,%%zmm16		\n\t	vpsrld	$12   ,%%zmm23,%%zmm23		\n\t"/* fq2>>12 gives q.hi32, which is at most 12-bits. */\
			/* Zero-extend each uint32 in ymm0-2 to a uint64, leaving q.lo32,md32,hi32 in qword-form in zmm0,1,2: */\
				"vpmovzxdq	%%ymm0,%%zmm0			\n\t	vpmovzxdq	%%ymm7,%%zmm7			\n\t	vpmovzxdq	%%ymm14,%%zmm14			\n\t	vpmovzxdq	%%ymm21,%%zmm21			\n\t"\
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
			/* Note VPGATHERDD may not use default-opmask k0.
			Re. inline-asm-with-opmask, this is what I first ran into when trying the 'obvious' {%%k1} syntax under both GCC and ICC,
			e.g. adding this line to an AVX-512 macro which compiles and runs fine sans any explicit opmask-register references:

				"vaddpd %%zmm3,%%zmm2,%%zmm1{%%k1} \n\t"	vaddpd %%zmm10,%%zmm9,%%zmm8{%%k2}	\n\t	vaddpd %%zmm17,%%zmm16,%%zmm15{%%k3}\n\t	vaddpd %%zmm24,%%zmm23,%%zmm22{%%k4}\n\t\

			gives "Assembler messages: Error: junk `%k1' after register". But, inside a basic (= non-extended) macro, with the
			escape-char-munge %% --> %, it compiles without issue. So it seems the culprit is GCC extended-asm. As far as a solution,
			thanks to Laurent Desnogues for turning up the following StackOverflow links.

			[1] http://stackoverflow.com/questions/21032395/masked-vector-instructions : Recommends double-curly-braces trickeration
				{{%%k1}} to make extended-asm play nice with AVX-512 opmask-register invocation. Alas, that hack is ICC-only; GCC gives
				"error: invalid 'asm': nested assembly dialect alternatives". That non-portability leads us to even-uglier hack number
			[2] http://stackoverflow.com/questions/34327831/invalid-asm-nested-assembly-dialect-alternatives : Explains that under GCC,
				double-curly-braces fail because
					"{} in GNU C inline asm already has a special meaning: providing alternatives for different
					ASM dialects. Using {{%%k1}} looks to gcc like nested alternatives, which is invalid."
				The GCC-workaround is the ugly %-escape-char-laden syntax %{%%k1%}, which, thankfully, works under both GCC and ICC.
				: */
				"movl	$-1,%%esi					\n\t"\
				"kmovw	%%esi,%%k1					\n\t	kmovw	%%esi,%%k2					\n\t	kmovw	%%esi,%%k3					\n\t	kmovw	%%esi,%%k4					\n\t"/* Init opmask k1 (Only need the low byte) */\
		"vpgatherqq (%%rbx,%%zmm4),%%zmm5%{%%k1%}	\n\tvpgatherqq (%%rbx,%%zmm11),%%zmm12%{%%k2%}\n\tvpgatherqq (%%rbx,%%zmm18),%%zmm19%{%%k3%}\n\tvpgatherqq (%%rbx,%%zmm25),%%zmm26%{%%k4%}\n\t"/* minv8 is a byte-table ... instruction sets mask-reg = 0, so each col uses separate same-valued k-reg */\
				"vpandq	%%zmm3,%%zmm5,%%zmm3		\n\t	vpandq	%%zmm10,%%zmm12,%%zmm10		\n\t	vpandq	%%zmm17,%%zmm19,%%zmm17		\n\t	vpandq	%%zmm24,%%zmm26,%%zmm24		\n\t"/* Iterative qinv-computation doesn't care about the garbage upper 7 bytes of each qword
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
				"vpmulld	%%zmm3,%%zmm4,%%zmm3	\n\t	vpmulld	%%zmm10,%%zmm11,%%zmm10		\n\t	vpmulld	%%zmm17,%%zmm18,%%zmm17		\n\t	vpmulld	%%zmm24,%%zmm25,%%zmm24		\n\t"/* qinv32 = MULL32(qinv16,(2-tmp32)); */\
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
			/* 78-bit truncation of qinv: */\
				"vpbroadcastd	%%esi,%%zmm31		\n\t"/* Mask-reg = 0x11...11 */\
				"vpsrlq	$50,%%zmm31,%%zmm31			\n\t"/* 0 all but bottom 14 bits of each qword */\
				"vpandq		%%zmm5,%%zmm31,%%zmm5	\n\t	vpandq		%%zmm12,%%zmm31,%%zmm12	\n\t	vpandq		%%zmm19,%%zmm31,%%zmm19	\n\t	vpandq		%%zmm26,%%zmm31,%%zmm26	\n\t"/* Hi 14 bits of 78-bit qinv. */\
			/* Reassemble into three 26-bit chunks, convert to double and write to local data store: */
				"vpbroadcastd	%%esi,%%zmm31		\n\t"/* Mask-reg = 0x11...11 */\
				"vpsrlq	$38,%%zmm31,%%zmm31			\n\t"/* 0 all but bottom 26 bits of each qword */\
				"vpsrlq	$26,%%zmm3,%%zmm6			\n\t	vpsrlq	$26,%%zmm10,%%zmm13			\n\t	vpsrlq	$26,%%zmm17,%%zmm20			\n\t	vpsrlq	$26,%%zmm24,%%zmm27			\n\t"/* hi6 bits of lo32 will go into md26... */\
				"vpsllq	$06,%%zmm4,%%zmm4			\n\t	vpsllq	$06,%%zmm11,%%zmm11			\n\t	vpsllq	$06,%%zmm18,%%zmm18			\n\t	vpsllq	$06,%%zmm25,%%zmm25			\n\t"/* md32<<6 can be done in place... */\
				"vpaddq		%%zmm6,%%zmm4,%%zmm4	\n\t	vpaddq		%%zmm13,%%zmm11,%%zmm11	\n\t	vpaddq		%%zmm20,%%zmm18,%%zmm18	\n\t	vpaddq		%%zmm27,%%zmm25,%%zmm25	\n\t"/* bits <26:63>, a 38-bit-wide temp. */\
				"vpandq		%%zmm3,%%zmm31,%%zmm3	\n\t	vpandq		%%zmm10,%%zmm31,%%zmm10	\n\t	vpandq		%%zmm17,%%zmm31,%%zmm17	\n\t	vpandq		%%zmm24,%%zmm31,%%zmm24	\n\t"/* qinv.lo26 */\
				"vpsrlq	$26,%%zmm4,%%zmm6			\n\t	vpsrlq	$26,%%zmm11,%%zmm13			\n\t	vpsrlq	$26,%%zmm18,%%zmm20			\n\t	vpsrlq	$26,%%zmm25,%%zmm27			\n\t"/* hi12 bits of md38 will go into hi26... */\
				"vpandq		%%zmm4,%%zmm31,%%zmm4	\n\t	vpandq		%%zmm11,%%zmm31,%%zmm11	\n\t	vpandq		%%zmm18,%%zmm31,%%zmm18	\n\t	vpandq		%%zmm25,%%zmm31,%%zmm25	\n\t"/* qinv.md26 */\
				"vpsllq	$12,%%zmm5,%%zmm5			\n\t	vpsllq	$12,%%zmm12,%%zmm12			\n\t	vpsllq	$12,%%zmm19,%%zmm19			\n\t	vpsllq	$12,%%zmm26,%%zmm26			\n\t"/* hi14<<12 can be done in place... */\
				"vpaddq		%%zmm6,%%zmm5,%%zmm5	\n\t	vpaddq		%%zmm13,%%zmm12,%%zmm12	\n\t	vpaddq		%%zmm20,%%zmm19,%%zmm19	\n\t	vpaddq		%%zmm27,%%zmm26,%%zmm26	\n\t"/* qinv.hi26 */\
			/* Down-convert each uint64 in zmm0-2 to a uint32, leaving qinv.lo32,md32,hi32 in dword-form in zmm3-5: */\
				"vpmovqd	%%zmm3,%%ymm3			\n\t	vpmovqd	%%zmm10,%%ymm10				\n\t	vpmovqd	%%zmm17,%%ymm17				\n\t	vpmovqd	%%zmm24,%%ymm24				\n\t"\
				"vpmovqd	%%zmm4,%%ymm4			\n\t	vpmovqd	%%zmm11,%%ymm11				\n\t	vpmovqd	%%zmm18,%%ymm18				\n\t	vpmovqd	%%zmm25,%%ymm25				\n\t"\
				"vpmovqd	%%zmm5,%%ymm5			\n\t	vpmovqd	%%zmm12,%%ymm12				\n\t	vpmovqd	%%zmm19,%%ymm19				\n\t	vpmovqd	%%zmm26,%%ymm26				\n\t"\
				"vcvtdq2pd	%%ymm3,%%zmm3			\n\t	vcvtdq2pd	%%ymm10,%%zmm10			\n\t	vcvtdq2pd	%%ymm17,%%zmm17			\n\t	vcvtdq2pd	%%ymm24,%%zmm24			\n\t"/* (double)fqinv0 */\
				"vcvtdq2pd	%%ymm4,%%zmm4			\n\t	vcvtdq2pd	%%ymm11,%%zmm11			\n\t	vcvtdq2pd	%%ymm18,%%zmm18			\n\t	vcvtdq2pd	%%ymm25,%%zmm25			\n\t"/* (double)fqinv1 */\
				"vcvtdq2pd	%%ymm5,%%zmm5			\n\t	vcvtdq2pd	%%ymm12,%%zmm12			\n\t	vcvtdq2pd	%%ymm19,%%zmm19			\n\t	vcvtdq2pd	%%ymm26,%%zmm26			\n\t"/* (double)fqinv2 */\
			/* Store base-2^26 double-based qinv: */\
				"vmovaps	%%zmm3,0x400(%%rax)		\n\t	vmovaps	%%zmm10,0x440(%%rax)		\n\t	vmovaps	%%zmm17,0x480(%%rax)		\n\t	vmovaps	%%zmm24,0x4c0(%%rax)		\n\t"/* store lo26 into fqinv0 */\
				"vmovaps	%%zmm4,0x500(%%rax)		\n\t	vmovaps	%%zmm11,0x540(%%rax)		\n\t	vmovaps	%%zmm18,0x580(%%rax)		\n\t	vmovaps	%%zmm25,0x5c0(%%rax)		\n\t"/* store lo52 into fqinv1 */\
				"vmovaps	%%zmm5,0x600(%%rax)		\n\t	vmovaps	%%zmm12,0x640(%%rax)		\n\t	vmovaps	%%zmm19,0x680(%%rax)		\n\t	vmovaps	%%zmm26,0x6c0(%%rax)		\n\t"/* store hi26 into fqinv2 */\
				:					/* outputs: none */\
				: [__fq0]  "m" (tmp)	/* All inputs from memory addresses here */\
				 ,[__p2]   "m" (dptr)\
				 ,[__kvec] "m" (tm1)	/* vec-of-doubles copy of input kvec */\
				 ,[__minv8]   "m" (minv8_ptr)	/* Ptr to Table of precomputed byte-inverses def'd in mi64.h */\
				 ,[__two26f] "m" (two26f)\
				: "cc","memory","rax","rbx","rcx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
			);

		#endif

			// In AVX2|512+FMA version, overwrite middle term of usual qinv[0,1,2] sequence with low 52 bits of qinv:
			dptr = fqinv0[0];
			__asm__ volatile (\
				"movq	%[__fqinv0],%%rdx	\n\t"\
				"movq	%[__two26f],%%rsi	\n\t"\
				"vmovaps	(%%rsi),%%zmm3	\n\t"/* TWO26FLOAT */\
				"vmovaps	     (%%rdx),%%zmm0		\n\t	vmovaps	0x040(%%rdx),%%zmm4		\n\t	vmovaps	0x080(%%rdx),%%zmm8		\n\t	vmovaps	0x0c0(%%rdx),%%zmm12	\n\t"\
				"vmovaps	0x100(%%rdx),%%zmm1		\n\t	vmovaps	0x140(%%rdx),%%zmm5		\n\t	vmovaps	0x180(%%rdx),%%zmm9		\n\t	vmovaps	0x1c0(%%rdx),%%zmm13	\n\t"\
				"vmovaps	%%zmm1,%%zmm2			\n\t	vmovaps	%%zmm5,%%zmm6			\n\t	vmovaps	%%zmm9,%%zmm10			\n\t	vmovaps	%%zmm13,%%zmm14			\n\t"/* cpy fhi into flo-destination reg */\
			"vfmadd132pd %%zmm3,%%zmm0,%%zmm2		\n\t vfmadd132pd %%zmm3,%%zmm4,%%zmm6	\n\t vfmadd132pd %%zmm3,%%zmm8,%%zmm10	\n\t vfmadd132pd %%zmm3,%%zmm12,%%zmm14	\n\t"/* (qinv0 + 2^26*qinv1)*/\
				"vmovaps	%%zmm2,0x100(%%rdx)		\n\t	vmovaps	%%zmm6,0x140(%%rdx)		\n\t	vmovaps	%%zmm10,0x180(%%rdx)	\n\t	vmovaps	%%zmm14,0x1c0(%%rdx)	\n\t"\
				:					/* outputs: none */\
				: [__fqinv0] "m" (dptr)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)	\
				: "cc","memory","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6", "xmm8","xmm9","xmm10", "xmm12","xmm13","xmm14"	/* Clobbered registers */\
			);

			// zshift is in [0,76]:
			if(zshift < 26) {
				dtmp = 1<<zshift;		for(j = 0; j < 32; j += 8) { VEC_DBL_INIT_8((vec_dbl*)fx0[j],dtmp); VEC_DBL_INIT_8((vec_dbl*)fx1[j], 0.0); VEC_DBL_INIT_8((vec_dbl*)fx2[j], 0.0); }
			} else if(zshift < 52) {
				dtmp = 1<<(zshift-26);	for(j = 0; j < 32; j += 8) { VEC_DBL_INIT_8((vec_dbl*)fx0[j], 0.0); VEC_DBL_INIT_8((vec_dbl*)fx1[j],dtmp); VEC_DBL_INIT_8((vec_dbl*)fx2[j], 0.0); }
			} else if(zshift < 78) {
				dtmp = 1<<(zshift-52);	for(j = 0; j < 32; j += 8) { VEC_DBL_INIT_8((vec_dbl*)fx0[j], 0.0); VEC_DBL_INIT_8((vec_dbl*)fx1[j], 0.0); VEC_DBL_INIT_8((vec_dbl*)fx2[j],dtmp); }
			} else {
				ASSERT(0,"zshift out of range!");
			}

			/*...x^2 mod q is returned in x. */
			// All 3-word-uint64-form operands have components in [0, 2^26).

		/* Inner loop body needs 42 MOVAPS, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */
		/* !ASM_LOOP: High-level loop construct, and explicit branch to do the modular doubling: */
			for(j = start_index-1; j >= 0; j--) {
				SSE2_twopmodq78_modmul_q32(fq0[0],fqinv0[0],fx0[0],two26i,pshift,j);
			}	/* for(j...) */
	#if 0
	====================================
	19 Nov 2106: Got rid of fugly special Iter1 code, but delving into resulting 1st-iter failures shows why adding a tiny const
				in UMULH computation fixes those failures:

	bc-based mont_sqr:
	a = 2^26; b = a^3
	p = 17732269; k = 133912345560; q = 2*k*p+1; q0 = q%a = 32195057; qhi52 = q/a = 70767692741
	i = 147853699419921 + 49188788*a^2
	x = 2^10

	j = 17:																					match?
	u = ?; h = ?		?
	v = 245718092949944494194688 = 17825792 + a*3661484911291964	[mull(lo,qinv)]
	l =      3861111986337213858 = 60943778 + a*57535052095			[umulh(lo,q)]			60943777 + a*57535052095, dropped a low bit!
	Here are details

	/* MULH78(q,lo,lo) --> lo =(q*lo)/2^78: flo0,1,2 enter in zmm0,1,2. NOTE: Copy of flo0,1 assumed in 0x00,0x100(rax). IN THE COMMENTS X = LO, Y = Q: */\
	/* Cost of FMA-wide-mul-based MULL78:  1 ADD, 12 MUL (8 FMA), 14 LD/ST (not counting reg-copy) */\
	/* Compare to non-FMA-version MULL78: 10 ADD, 15 MUL (0 FMA), 18 LD/ST (not counting reg-copy); thus MULs ~same, but far fewer ADD. */\
		"movq	%[__aq0],%%rdx				\n\t"\
		"vmovaps	-0x40(%%rbx),%%zmm3		\n\t"/* TWO26FLOAT */\
	"vfmadd213pd	%%zmm1,%%zmm3,%%zmm0	\n\t"/* hi52 = x1 + x2*2^26; overwrites x2 */\	3661484911291964
		"vmovaps	0x300(%%rdx),%%zmm1		\n\t"/* Load mi52 = qhi52 */\					70767692741
		"vmovaps	%%zmm0,0x200(%%rax)		\n\t"/* Save copy of hi52 in flo2 slot */\
	/* Bits <52:155> : */\
		"vmulpd		%%zmm0,%%zmm1,%%zmm2	\n\t"/* fhi1 = hi52*mi52 */\					259114839178117341979344896 [top 53 sigbits]
		"vmovaps	%%zmm2,%%zmm3			\n\t"/* cpy fhi into flo-destination reg */\
	"vfmsub231pd	%%zmm0,%%zmm1,%%zmm3	\n\t"/* ftmp = fma(hi52,mi52,-fhi1) */\			7715088428
		"vmovaps	%%zmm2,%%zmm0			\n\t"/* cpy fhi into fcy-destination reg */\
		"vmulpd	  0x80(%%rbx),%%zmm2,%%zmm2	\n\t"/*             fhi1*TWO52FLINV  */\		57535052095
		"vrndscalepd $0,%%zmm2,%%zmm2		\n\t"/* fhh = DNINT(fhi1*TWO52FLINV)	This part remains in fhi1... */\
	"vfnmadd231pd 0x40(%%rbx),%%zmm2,%%zmm0	\n\t"/* fcy = fma(fhh,-TWO52FLOAT,fhi1)	Backward carry from hiA into loA */\	2333266753355776
		"vaddpd	%%zmm3,%%zmm0,%%zmm0		\n\t"/* ftmp += fcy */\		2333274468444204
	/* fhi1 in zmm2; ftmp in zmm0; hi52 in 0x200(rax); mi52 in zmm1; zmm3 free */\
	/* Bits <26:103> - 2 separate cross-products, only need high ~53 bits of each: */\
		"vmovaps	      (%%rdx),%%zmm3	\n\t"/* Reload y0 (a.k.a. q0) and *= 2^26 */\	32195057
		"vmulpd	 0x200(%%rax),%%zmm3,%%zmm3	\n\t"/* fhh = hi52*fy0 */\						117881715423684724621948
	"vfmadd231pd      (%%rax),%%zmm1,%%zmm3	\n\t "/* fhh = fma(mi52,fx0, fhh) */\			117882976913855845597820
	/* fhi1 in zmm2; ftmp in zmm0; fhh in zmm3; zmm1 free */\
		"vmovaps		(%%rbx),%%zmm1		\n\t"/* TWO26FLINV */\
		"vmulpd		%%zmm1,%%zmm0,%%zmm0	\n\t"/* ftmp *= TWO26FLINV  */\					2333274468444204/2^26 = 34768498.964
		"vmovaps	0x80(%%rbx),%%zmm1		\n\t"/* TWO52FLINV */\
	"vfmadd231pd	%%zmm1,%%zmm3,%%zmm0	\n\t"/* ftmp*TWO26FLINV + fhh*TWO52FLINV */\	34768498.9637763 + 26175279.0362237 = 60943777.999999998, dropped a bit!
		"vrndscalepd $1,%%zmm0,%%zmm0		\n\t"/* fhi0 = floor("		") */\
	/* Digit 3 (lo0 word of MULH78 result) in zmm0: */\
	/* Digits 4,5 (lo1,2 words of MULH78 result) remain in a 52-bit double (zmm2): */\	code gives xmm0 = -6165086.0000000037 = 60943777.9999999963-2^26

	Here is bitfield governing 2x in modpow: need low (start_index = 18) bits, 010110110100000100:
	j	bit	mod_sqr sequence																		xout			code gives same?
	--	--	------------------------------------------------------------------------------------	---------------		----
	24
	23
	22
	21
	20
	19
	18
	17	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				888027481444537423	...424, 1 too large!
	16	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	3153079875297589710
	15	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				3349662228228641186
	14	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	2479112533014143301
	13	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	1818305489032730849
	12	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				2069090628149556142
	11	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	3079917231797630540
	10	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	1571057852565897090
	9	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				2134637900068426919
	8	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	2814339156985005729
	7	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				2973371148533820384
	6	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				2176732355512285298
	5	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				3519776166192070528
	4	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				2305766521945167233
	3	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				1144711116339289293
	2	1	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x = 2*x%q ; x	2014514757755243580
	1	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				740117870607009439
	0	0	x *= x ; u = x%b ; h = x/b ; v = u*i%b ; l = v*q/b ; x = h-l+q*(h<l) ; x				2374569733890875641
	...and a final mod-doubling gives x = 1, as desired.
	====================================
	#endif

			/* Do a final mod-doubling and return. Code here is specialized for the case
			where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2 or -(q-1)/2:
			*/
			dptr = fq0[0];
			__asm__ volatile (\
				"movq	%[__fq0]   ,%%rax	\n\t"\
				"movq	%[__two26f],%%rbx	\n\t"\
				"vmovaps	0x40(%%rbx),%%zmm3		\n\t	vmovaps	(%%rbx),%%zmm7				\n\t"/* two26i,f */\
			"vfmadd213pd %%zmm1,%%zmm7,%%zmm2		\n\t vfmadd213pd %%zmm5,%%zmm7,%%zmm6		\n\t vfmadd213pd %%zmm9,%%zmm7,%%zmm10			\n\t vfmadd213pd %%zmm13,%%zmm7,%%zmm14		\n\t"/* (x1 + 2^26*x2) = hi52 */\
				/* Mod-doubling code taken from SSE2_twopmodq78_modmul_q16, expects x.lo26,hi52 in zmm1,2, resp., so copy zmm0,4,8,12 -> zmm1,5,9,13 after computing x.hi52 (in zmm2,6,10,14): */\
				"vaddpd	%%zmm0,%%zmm0,%%zmm1		\n\t	vaddpd	%%zmm4,%%zmm4,%%zmm5		\n\t	vaddpd	%%zmm8 ,%%zmm8 ,%%zmm9			\n\t vaddpd	%%zmm12,%%zmm12,%%zmm13			\n\t"/* low 26 bits */\
				"vaddpd	%%zmm2,%%zmm2,%%zmm2		\n\t	vaddpd	%%zmm6,%%zmm6,%%zmm6		\n\t	vaddpd	%%zmm10,%%zmm10,%%zmm10			\n\t vaddpd	%%zmm14,%%zmm14,%%zmm14			\n\t"/* top 52 bits */\
				/* And subtract q.lo26,hi52 - note q.lo26 is balanced-digit normalized and divided by 2^26, so do this sub prior to final >0 normalization below: */\
				"vsubpd	     (%%rax),%%zmm1,%%zmm1	\n\t	vsubpd	0x040(%%rax),%%zmm5,%%zmm5	\n\t	vsubpd	0x080(%%rax),%%zmm9 ,%%zmm9		\n\t	vsubpd	0x0c0(%%rax),%%zmm13,%%zmm13\n\t"/* -= qlo26 */\
				"vsubpd	0x300(%%rax),%%zmm2,%%zmm2	\n\t	vsubpd	0x340(%%rax),%%zmm6,%%zmm6	\n\t	vsubpd	0x380(%%rax),%%zmm10,%%zmm10	\n\t	vsubpd	0x3c0(%%rax),%%zmm14,%%zmm14\n\t"/* -= qhi52 */\
				/* Final-residue computation needs digits >= 0 normalization - since fq0 will get folded in with fq2 anyway, only need renormalize fq0: */\
				"vmulpd	%%zmm3,%%zmm1,%%zmm1		\n\t	vmulpd	%%zmm3,%%zmm5,%%zmm5		\n\t	vmulpd	%%zmm3,%%zmm9,%%zmm9			\n\t vmulpd	%%zmm3,%%zmm13,%%zmm13			\n\t"/* xlo *= two26i */\
				"vrndscalepd $1,%%zmm1,%%zmm0		\n\t	vrndscalepd $1,%%zmm5,%%zmm4		\n\t	vrndscalepd $1,%%zmm9,%%zmm8			\n\t vrndscalepd $1,%%zmm13,%%zmm12			\n\t"/* fcy */\
				"vsubpd	%%zmm0,%%zmm1,%%zmm1		\n\t	vsubpd	%%zmm4,%%zmm5,%%zmm5		\n\t	vsubpd	%%zmm8,%%zmm9,%%zmm9			\n\t vsubpd	%%zmm12,%%zmm13,%%zmm13			\n\t"/* fx0 -= fcy */\
				"vaddpd	%%zmm0,%%zmm2,%%zmm2		\n\t	vaddpd	%%zmm4,%%zmm6,%%zmm6		\n\t	vaddpd	%%zmm8,%%zmm10,%%zmm10			\n\t vaddpd	%%zmm12,%%zmm14,%%zmm14			\n\t"/* Add carry into xhi */\
				"vmulpd	%%zmm7,%%zmm1,%%zmm1		\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5		\n\t	vmulpd	%%zmm7,%%zmm9,%%zmm9			\n\t vmulpd	%%zmm7 ,%%zmm13,%%zmm13			\n\t"/* fx0 *= two26f ... result must go into m0 */\
				/* (lo26 == 1 && hi52 == 0), save result in a bitmask: */
				"vmovaps	(%%rbx),%%zmm0			\n\t"/* TWO26FLOAT */\
				"vmulpd	0x40(%%rbx),%%zmm0,%%zmm3	\n\t"/* 1.0 ... no good way to generate 1.0 on-the-fly, so use 2^26 * 2^-26 */\
				"vpxorq	     %%zmm0,%%zmm0,%%zmm0	\n\t"/* 0.0 */\
				"vcmppd $0,%%zmm1,%%zmm3,%%k0	\n\t	vcmppd $0,%%zmm5,%%zmm3,%%k1	\n\t	vcmppd $0,%%zmm9,%%zmm3,%%k2	\n\t	vcmppd $0,%%zmm13,%%zmm3,%%k3	\n\t"/* low 26 bits == 1 ? In AVX512 version, mask-regs k0-3 replace dest-regs m1,5,9,13 */\
				"vcmppd $0,%%zmm2,%%zmm0,%%k4	\n\t	vcmppd $0,%%zmm6,%%zmm0,%%k5	\n\t	vcmppd $0,%%zmm10,%%zmm0,%%k6	\n\t	vcmppd $0,%%zmm14,%%zmm0,%%k7	\n\t"/* top 52 bits == 0 ? In AVX512 version, mask-regs k4-7 replace dest-regs m0,4,8,12 */\
				"kandw	%%k0,%%k4,%%k0			\n\t	kandw	%%k1,%%k5,%%k1			\n\t	kandw	%%k2,%%k6,%%k2			\n\t	kandw	%%k3,%%k7,%%k3			\n\t"/* and the 2 results together into four 8-bit bitmasks */\
				"kmovw	%%k0,%%eax				\n\t	kmovw	%%k1,%%ebx				\n\t	kmovw	%%k2,%%ecx				\n\t	kmovw	%%k3,%%edx				\n\t"\
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
				 ,[__two26f] "m" (two26f)	\
				 ,[__result] "m" (r)	\
				: "cc","memory","cl","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
			);

			return r;
		}

		uint64 twopmodq78_3WORD_DOUBLE_q64(uint64 p, uint64 k[], int init_sse2, int thr_id)
		{
			const char func[] = "twopmodq78_3WORD_DOUBLE_q64";
			 int32 j = 0;	/* This needs to be signed because of the LR binary exponentiation. */
			uint64 lead7, r = 0, tmp64;
			static uint64 psave = 0, pshift;
			static uint32 start_index, zshift, first_entry = TRUE;
			uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case
			uint8* minv8_ptr = minv8;	// Ptr to Table of precomputed byte-inverses def'd in mi64.h
			static int max_threads = 1;	// Default local-array-init is for just a single thread ... caller can re-init for > 1 threads later, if desired.
		#ifdef USE_AVX512_I
		  #error AVX-512 IFMA instruction extensions not yet supported!
			static uint64 *sc_arr = 0x0, *sc_ptr;
			uint64 *itmp;
			vec_u64 *tmp;
		  #ifdef MULTITHREAD
			static uint64 *__r0;	// Base address for discrete per-thread local stores ...
									// *NOTE*: more convenient to use uint64* rather than vec_u64* here
		  #else
			static	// Following set of pointers only static-izable in single-thread mode
		  #endif
			uint64 *fq0[64],*fq1[64],*fq2[64],*fqhi52[64], *fqinv0[64],*fqinv1[64],*fqinv2[64], *fx0[64],*fx1[64],*fx2[64],
					*mask_lo26,*mask_lo52;
			for(j = 0; j < 64; j++) {
				ASSERT((k[j] >> 52) == 0ull, "Ks must be < 2^52!");
			}
		#else
			static double *sc_arr = 0x0, *sc_ptr;
			double *dptr, dtmp;
			vec_dbl *tmp,*tm1;
		  #ifdef MULTITHREAD
			static double *__r0;	// Base address for discrete per-thread local stores ...
		  #else
			static	// Following set of pointers only static-izable in single-thread mode
		  #endif
			double *fq0[64],*fq1[64],*fq2[64],*fqhi52[64], *fqinv0[64],*fqinv1[64],*fqinv2[64], *fx0[64],*fx1[64],*fx2[64],
					*two13i, *two26f,*two26i, *two52f,*two52i,
					kdbl[64];
			// AVX-512 Foundation lacks the needed DQ extensions, so use HLL to convert kvec entries to double:
			for(j = 0; j < 64; j++) {
				ASSERT((k[j] >> 52) == 0ull, "Ks must be < 2^52!");
				kdbl[j] = (double)k[j];
			}
		#endif

			if(p != psave)
			{
			//	first_entry = FALSE;
				psave  = p;
				pshift = p + 78;
				j = leadz64(pshift);
				/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 77: */
				lead7 = ((pshift<<j) >> 57);
				if(lead7 > 77)
				{
					lead7 >>= 1;
					start_index =  64-j-6;	/* Use only the leftmost 6 bits */
				}
				else
					start_index =  64-j-7;

				zshift = 77 - lead7;
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

				sc_arr = ALLOC_UINT64(sc_arr, 0x380*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
				sc_ptr = (uint64 *)ALIGN_VEC_U64(sc_arr);	// Force vec_u64-alignment
				ASSERT(((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			  #ifdef MULTITHREAD
				__r0  = sc_ptr;
				mask_lo26 = sc_ptr + 0x300;
				mask_lo52 = sc_ptr + 0x310;
				for(j = 0; j < max_threads; ++j) {
					/* These remain fixed within each per-thread local store: */
					// Need specific-length inits here since might be building overall in AVX-512 mode:
					VEC_U64_INIT((vec_u64*)mask_lo26, 0x0000000003FFFFFFull);
					VEC_U64_INIT((vec_u64*)mask_lo52, 0x000FFFFFFFFFFFFFull);
					// Move on to next thread's local store:
					mask_lo26 += 0x380;
					mask_lo52 += 0x380;
				}
			  #else
				/* Remember, these are POINTERS-TO-UINT64, so need an increment of 8 to span an AVX-512 register: */
				fq0   [0] = sc_ptr + 0x000;
				fq1   [0] = sc_ptr + 0x040;
				fq2   [0] = sc_ptr + 0x080;
				fqhi52[0] = sc_ptr + 0x0C0;
				fqinv0[0] = sc_ptr + 0x100;
				fqinv1[0] = sc_ptr + 0x140;
				fqinv2[0] = sc_ptr + 0x180;
				fx0   [0] = sc_ptr + 0x1C0;
				fx1   [0] = sc_ptr + 0x200;
				fx2   [0] = sc_ptr + 0x240;
				for(j = 1, itmp = sc_ptr+1; j < 64; j++, itmp++) {
					fq0   [j] = itmp + 0x000;
					fq1   [j] = itmp + 0x040;
					fq2   [j] = itmp + 0x080;
					fqhi52[j] = itmp + 0x0C0;
					fqinv0[j] = itmp + 0x100;
					fqinv1[j] = itmp + 0x140;
					fqinv2[j] = itmp + 0x180;
					fx0   [j] = itmp + 0x1C0;
					fx1   [j] = itmp + 0x200;
					fx2   [j] = itmp + 0x240;
				}
			// +0x140,160 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
			// Consts: only get 1 AVX-512 reg (8 uint64s) per:
				mask_lo26 = sc_ptr + 0x300;
				mask_lo52 = sc_ptr + 0x310;
			// Total: Equivalent of 12*4 + 2 = 50 AVX-512 slots, alloc 56
				VEC_U64_INIT((vec_u64*)mask_lo26, 0x0000000003FFFFFFull);
				VEC_U64_INIT((vec_u64*)mask_lo52, 0x000FFFFFFFFFFFFFull);
			  #endif

			/***************************************************/
			#else	// Default AVX-512 floating-point-FMA mode
			/***************************************************/

				sc_arr = ALLOC_DOUBLE(sc_arr, 0x380*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
				sc_ptr = (uint64 *)ALIGN_VEC_DBL(sc_arr);	// Force vec_u64-alignment
				ASSERT(((uintptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			  #ifdef MULTITHREAD
				__r0  = sc_ptr;
				two13i = sc_ptr + 0x300;
				two26f = sc_ptr + 0x308;
				two26i = sc_ptr + 0x310;
				two52f = sc_ptr + 0x318;
				two52i = sc_ptr + 0x320;
				for(j = 0; j < max_threads; ++j) {
					/* These remain fixed within each per-thread local store: */
					// Need specific-length inits here since might be building overall in AVX-512 mode:
					VEC_DBL_INIT((vec_dbl*)two13i, TWO13FLINV);
					VEC_DBL_INIT((vec_dbl*)two26f, TWO26FLOAT);
					VEC_DBL_INIT((vec_dbl*)two26i, TWO26FLINV);
					VEC_DBL_INIT((vec_dbl*)two52f, TWO52FLOAT);
					VEC_DBL_INIT((vec_dbl*)two52i, TWO52FLINV);
					// Move on to next thread's local store:
					two13i   += 0x380;
					two26f   += 0x380;
					two26i   += 0x380;
					two52f   += 0x380;
					two52i   += 0x380;
				}
			  #else
				/* Remember, these are POINTERS-TO-DOUBLE, so need an increment of 8 to span an AVX-512 register: */
				fq0   [0] = sc_ptr + 0x000;
				fq1   [0] = sc_ptr + 0x040;
				fq2   [0] = sc_ptr + 0x080;
				fqhi52[0] = sc_ptr + 0x0C0;
				fqinv0[0] = sc_ptr + 0x100;
				fqinv1[0] = sc_ptr + 0x140;
				fqinv2[0] = sc_ptr + 0x180;
				fx0   [0] = sc_ptr + 0x1C0;
				fx1   [0] = sc_ptr + 0x200;
				fx2   [0] = sc_ptr + 0x240;
				for(j = 1, dptr = sc_ptr+1; j < 64; j++, dptr++) {
					fq0   [j] = dptr + 0x000;
					fq1   [j] = dptr + 0x040;
					fq2   [j] = dptr + 0x080;
					fqhi52[j] = dptr + 0x0C0;
					fqinv0[j] = dptr + 0x100;
					fqinv1[j] = dptr + 0x140;
					fqinv2[j] = dptr + 0x180;
					fx0   [j] = dptr + 0x1C0;
					fx1   [j] = dptr + 0x200;
					fx2   [j] = dptr + 0x240;
				}
			// +0x280,2c0 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
			// Consts: only get 1 AVX-512 reg (8 doubles) per:
				two13i = sc_ptr + 0x300;
				two26f = sc_ptr + 0x308;
				two26i = sc_ptr + 0x310;
				two52f = sc_ptr + 0x318;
				two52i = sc_ptr + 0x320;
			// Total: Equivalent of 12*8 + 5 = 101 AVX-512 slots, alloc 112
				VEC_DBL_INIT((vec_dbl*)two13i, TWO13FLINV);
				VEC_DBL_INIT((vec_dbl*)two26f, TWO26FLOAT);
				VEC_DBL_INIT((vec_dbl*)two26i, TWO26FLINV);
				VEC_DBL_INIT((vec_dbl*)two52f, TWO52FLOAT);
				VEC_DBL_INIT((vec_dbl*)two52i, TWO52FLINV);
			  #endif

			#endif
				if(init_sse2) return 0;
			}	/* end of inits */

			/* If multithreaded, set the local-store pointers needed for the current thread; */
		#ifdef MULTITHREAD

			ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
			sc_ptr = __r0 + thr_id*0x380;

		  #ifdef USE_AVX512_I

			/* Remember, these are POINTERS-TO-UINT64, so need an increment of 8 to span an AVX-512 register: */
			fq0   [0] = sc_ptr + 0x000;
			fq1   [0] = sc_ptr + 0x040;
			fq2   [0] = sc_ptr + 0x080;
			fqhi52[0] = sc_ptr + 0x0C0;
			fqinv0[0] = sc_ptr + 0x100;
			fqinv1[0] = sc_ptr + 0x140;
			fqinv2[0] = sc_ptr + 0x180;
			fx0   [0] = sc_ptr + 0x1C0;
			fx1   [0] = sc_ptr + 0x200;
			fx2   [0] = sc_ptr + 0x240;
			for(j = 1, itmp = sc_ptr+1; j < 64; j++, itmp++) {
				fq0   [j] = itmp + 0x000;
				fq1   [j] = itmp + 0x040;
				fq2   [j] = itmp + 0x080;
				fqhi52[j] = itmp + 0x0C0;
				fqinv0[j] = itmp + 0x100;
				fqinv1[j] = itmp + 0x140;
				fqinv2[j] = itmp + 0x180;
				fx0   [j] = itmp + 0x1C0;
				fx1   [j] = itmp + 0x200;
				fx2   [j] = itmp + 0x240;
			}
			mask_lo26 = sc_ptr + 0x300;
			mask_lo52 = sc_ptr + 0x308;
			VEC_U64_INIT((vec_u64*)mask_lo26, 0x0000000003FFFFFFull);
			VEC_U64_INIT((vec_u64*)mask_lo52, 0x000FFFFFFFFFFFFFull);

		  #else

			/* Remember, these are POINTERS-TO-DOUBLE, so need an increment of 8 to span an AVX-512 register: */
			fq0   [0] = sc_ptr + 0x000;
			fq1   [0] = sc_ptr + 0x040;
			fq2   [0] = sc_ptr + 0x080;
			fqhi52[0] = sc_ptr + 0x0C0;
			fqinv0[0] = sc_ptr + 0x100;
			fqinv1[0] = sc_ptr + 0x140;
			fqinv2[0] = sc_ptr + 0x180;
			fx0   [0] = sc_ptr + 0x1C0;
			fx1   [0] = sc_ptr + 0x200;
			fx2   [0] = sc_ptr + 0x240;
			for(j = 1, dptr = sc_ptr+1; j < 64; j++, dptr++) {
				fq0   [j] = dptr + 0x000;
				fq1   [j] = dptr + 0x040;
				fq2   [j] = dptr + 0x080;
				fqhi52[j] = dptr + 0x0C0;
				fqinv0[j] = dptr + 0x100;
				fqinv1[j] = dptr + 0x140;
				fqinv2[j] = dptr + 0x180;
				fx0   [j] = dptr + 0x1C0;
				fx1   [j] = dptr + 0x200;
				fx2   [j] = dptr + 0x240;
			}
			two13i = sc_ptr + 0x300;
			two26f = sc_ptr + 0x308;
			two26i = sc_ptr + 0x310;
			two52f = sc_ptr + 0x318;
			two52i = sc_ptr + 0x320;
			VEC_DBL_INIT((vec_dbl*)two13i, TWO13FLINV);
			VEC_DBL_INIT((vec_dbl*)two26f, TWO26FLOAT);
			VEC_DBL_INIT((vec_dbl*)two26i, TWO26FLINV);
			VEC_DBL_INIT((vec_dbl*)two52f, TWO52FLOAT);
			VEC_DBL_INIT((vec_dbl*)two52i, TWO52FLINV);

		  #endif

		#endif

		#ifdef MUL_LOHI64_SUBROUTINE
			#error MUL_LOHI64_SUBROUTINE defined!
		#endif
			ASSERT((p >> 63) == 0, "twopmodq78_q64: p must be < 2^63!");

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
			dtmp = p+p;	dptr = &dtmp;	// dtmp = double; dptr = double*
			tmp = (vec_dbl*)fq0[0];	tm1 = (vec_dbl*)kdbl;	// Use vec_dbl* ptrs to hold address args ... asm-code cares not a whit about types
			__asm__ volatile (\
				"movq	%[__fq0] ,%%rax				\n\t"\
				"movq	%[__p2]  ,%%rbx				\n\t"\
				"movq	%[__kvec],%%rcx				\n\t"/* Each column uses 7 distinct vector registers, so increment all non-shared vec-reg indices by +7 for each additional rightward column: */\
				"vbroadcastsd		(%%rbx),%%zmm1	\n\t"/* broadcast 2*p to all 8 qwords of zmm ... double-prec version is from mem-address */\
				"vmovupd	0x000(%%rcx),%%zmm3		\n\t	vmovupd	0x040(%%rcx),%%zmm10		\n\t	vmovupd	0x080(%%rcx),%%zmm17		\n\t	vmovupd	0x0c0(%%rcx),%%zmm24		\n\t"/* load the factor-candidate k's, 8 per zmm ... note kvec not nec. 16-byte aligned */\
				"vmulpd		%%zmm1,%%zmm3,%%zmm2	\n\t	vmulpd		%%zmm1,%%zmm10,%%zmm9	\n\t	vmulpd		%%zmm1,%%zmm17,%%zmm16	\n\t	vmulpd		%%zmm1,%%zmm24,%%zmm23	\n\t"/* hi = x*y = 2p*k */\
				"vmovaps	%%zmm2,%%zmm0			\n\t	vmovaps	%%zmm9,%%zmm7				\n\t	vmovaps	%%zmm16,%%zmm14				\n\t	vmovaps	%%zmm23,%%zmm21				\n\t"/* cpy hi into lo-destination reg */\
			"vfmsub231pd	%%zmm1,%%zmm3,%%zmm0	\n\t vfmsub231pd	%%zmm1,%%zmm10,%%zmm7	\n\t vfmsub231pd	%%zmm1,%%zmm17,%%zmm14	\n\t vfmsub231pd	%%zmm1,%%zmm24,%%zmm21	\n\t"/* lo = fma(x,y,-hi; xmm1,3 FREE */\
				"movl	$1,%%esi					\n\t"\
				"vpbroadcastd	%%esi,%%zmm3		\n\t	vpbroadcastd	%%esi,%%zmm10		\n\t	vpbroadcastd	%%esi,%%zmm17		\n\t	vpbroadcastd	%%esi,%%zmm24		\n\t"/* Broadcast to all 8 int32 slots of zmm3, since KNL lacks AVX-512 VL extensions */\
				"vcvtdq2pd	%%ymm3,%%zmm3			\n\t	vcvtdq2pd	%%ymm10,%%zmm10			\n\t	vcvtdq2pd	%%ymm17,%%zmm17			\n\t	vcvtdq2pd	%%ymm24,%%zmm24			\n\t"/* zmm3 = (double)1, conversion from 4-way int32 in lower half of zmm3 */\
				"vaddpd		%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd		%%zmm10,%%zmm7,%%zmm7	\n\t	vaddpd		%%zmm17,%%zmm14,%%zmm14	\n\t	vaddpd		%%zmm24,%%zmm21,%%zmm21	\n\t"/* lo += 1 */\
			/* Normalize to properly split the resultant significant bits across lo,hi, use base-2^52 signed-int-stored-as-double: */\
				"movq	%[__two26f],%%rbx			\n\t"\
			"vbroadcastsd	0xc0(%%rbx),%%zmm31		\n\t"/* TWO52FLINV */\
			"vbroadcastsd	0x80(%%rbx),%%zmm30		\n\t"/* TWO52FLOAT */\
			"vbroadcastsd	0x40(%%rbx),%%zmm29		\n\t"/* TWO26FLINV */\
			"vbroadcastsd		(%%rbx),%%zmm28		\n\t"/* TWO26FLOAT */\
				"vmovaps	%%zmm2,%%zmm1			\n\t	vmovaps	%%zmm9,%%zmm8				\n\t	vmovaps	%%zmm16,%%zmm15				\n\t	vmovaps	%%zmm23,%%zmm22				\n\t"/* cpy hi into hh-destination reg */\
				"vmulpd		%%zmm31,%%zmm1,%%zmm1	\n\t	vmulpd		%%zmm31,%%zmm8,%%zmm8	\n\t	vmulpd		%%zmm31,%%zmm15,%%zmm15	\n\t	vmulpd		%%zmm31,%%zmm22,%%zmm22	\n\t"/*            hi*TWO52FLINV  */\
				"vrndscalepd	$1,%%zmm1,%%zmm1	\n\t	vrndscalepd	$1,%%zmm8,%%zmm8		\n\t	vrndscalepd	$1,%%zmm15,%%zmm15		\n\t	vrndscalepd	$1,%%zmm22,%%zmm22		\n\t"/* hh = FLOOR(hi*TWO52FLINV) */\
				"vmovaps	%%zmm2,%%zmm3			\n\t	vmovaps	%%zmm9,%%zmm10				\n\t	vmovaps	%%zmm16,%%zmm17				\n\t	vmovaps	%%zmm23,%%zmm24				\n\t"/* cpy hi into cy-destination reg */\
			"vfnmadd231pd	%%zmm30,%%zmm1,%%zmm3	\n\t vfnmadd231pd	%%zmm30,%%zmm8,%%zmm10	\n\t vfnmadd231pd	%%zmm30,%%zmm15,%%zmm17	\n\t vfnmadd231pd	%%zmm30,%%zmm22,%%zmm24\n\t"/* cy = fma(hh,-TWO52FLOAT,hi)	Backward carry from hi into lo */\
				"vaddpd		%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd		%%zmm10,%%zmm7,%%zmm7	\n\t	vaddpd		%%zmm17,%%zmm14,%%zmm14	\n\t	vaddpd		%%zmm24,%%zmm21,%%zmm21	\n\t"/* lo += cy */\
				"vmovaps	%%zmm1,0x400(%%rax)		\n\t	vmovaps	%%zmm8,0x440(%%rax)			\n\t	vmovaps	%%zmm15,0x480(%%rax)		\n\t	vmovaps	%%zmm22,0x4c0(%%rax)		\n\t"/* store hi26 into fq2 */\
			/* Normalize to convert 78-bit q from balanced-digit [lo52,hi26] representation to nonnegative-digit [lo26,hi52]: */\
				"vmovaps	%%zmm0,%%zmm2			\n\t	vmovaps	%%zmm7,%%zmm9				\n\t	vmovaps	%%zmm14,%%zmm16				\n\t	vmovaps	%%zmm21,%%zmm23				\n\t"/* cpy lo into destination reg */\
				"vmulpd		%%zmm29,%%zmm2,%%zmm2	\n\t	vmulpd		%%zmm29,%%zmm9,%%zmm9	\n\t	vmulpd		%%zmm29,%%zmm16,%%zmm16	\n\t	vmulpd		%%zmm29,%%zmm23,%%zmm23	\n\t"/*            lo*TWO26FLINV  */\
				"vrndscalepd	$1,%%zmm2,%%zmm2	\n\t	vrndscalepd	$1,%%zmm9,%%zmm9		\n\t	vrndscalepd	$1,%%zmm16,%%zmm16		\n\t	vrndscalepd	$1,%%zmm23,%%zmm23		\n\t"/* cy = FLOOR(lo*TWO26FLINV)	Forward carry from lo into hi */\
				"vmovaps	%%zmm2,0x200(%%rax)		\n\t	vmovaps	%%zmm9,0x240(%%rax)			\n\t	vmovaps	%%zmm16,0x280(%%rax)		\n\t	vmovaps	%%zmm23,0x2c0(%%rax)		\n\t"/* store cy into fq1 */\
			"vfnmadd231pd	%%zmm28,%%zmm2,%%zmm0	\n\t vfnmadd231pd	%%zmm28,%%zmm9,%%zmm7	\n\t vfnmadd231pd	%%zmm28,%%zmm16,%%zmm14	\n\t vfnmadd231pd	%%zmm28,%%zmm23,%%zmm21	\n\t"/* lo26 = fma(cy,-TWO26FLOAT,lo) */\
			" vfmadd132pd	%%zmm28,%%zmm2,%%zmm1	\n\t vfmadd132pd	%%zmm28,%%zmm9,%%zmm8	\n\t vfmadd132pd	%%zmm28,%%zmm16,%%zmm15	\n\t vfmadd132pd	%%zmm28,%%zmm23,%%zmm22	\n\t"/* hi52 = fma(hi,+TWO26FLOAT,cy) */\
				"vmovaps	%%zmm0,0x000(%%rax)		\n\t	vmovaps	%%zmm7,0x040(%%rax)			\n\t	vmovaps	%%zmm14,0x080(%%rax)		\n\t	vmovaps	%%zmm21,0x0c0(%%rax)		\n\t"/* store lo26 into fq0 */\
				"vmovaps	%%zmm1,0x600(%%rax)		\n\t	vmovaps	%%zmm8,0x640(%%rax)			\n\t	vmovaps	%%zmm15,0x680(%%rax)		\n\t	vmovaps	%%zmm22,0x6c0(%%rax)		\n\t"/* store hi52 into fqhi52 */\
			/*** hi = hh ==> lo26 in zmm0, hi52 = zmm1. ***/\
			/* Assemble 32-bit pieces of q from the 26-bit floating-point ones ... AVX-512F has just double->int32, so
			convert our three 26-bit double-vec chunks to 8-way int32s (i.e. zmm -> ymm), then re-expand result to 8-way
			int64s (with high 32-bit halves all zero) in preparation for mod-inverse computation: */\
				"vcvttpd2dq		  %%zmm0,%%ymm0		\n\t	vcvttpd2dq		  %%zmm7,%%ymm7		\n\t	vcvttpd2dq		  %%zmm14,%%ymm14	\n\t	vcvttpd2dq		  %%zmm21,%%ymm21	\n\t"/* (int32)fq0 */\
				"vcvttpd2dq	0x200(%%rax),%%ymm1		\n\t	vcvttpd2dq	0x240(%%rax),%%ymm8		\n\t	vcvttpd2dq	0x280(%%rax),%%ymm15	\n\t	vcvttpd2dq	0x2c0(%%rax),%%ymm22	\n\t"/* (int32)fq1 */\
				"vcvttpd2dq	0x400(%%rax),%%ymm2		\n\t	vcvttpd2dq	0x440(%%rax),%%ymm9		\n\t	vcvttpd2dq	0x480(%%rax),%%ymm16	\n\t	vcvttpd2dq	0x4c0(%%rax),%%ymm23	\n\t"/* (int32)fq2 */\
		/*** Must use full-length zmm-regs here (even though upper 256 bits unused) due to AVX-512F full-width constraint (applying to ymm needs AVX-512VL extensions): ***/\
				"vpslld	$26,%%zmm1,%%zmm3			\n\t	vpslld	$26,%%zmm8,%%zmm10			\n\t	vpslld	$26,%%zmm15,%%zmm17			\n\t	vpslld	$26,%%zmm22,%%zmm24			\n\t"/* left-justify lo6 bits of fq1... */\
			"vpaddd	%%zmm3,%%zmm0,%%zmm0			\n\t	vpaddd	%%zmm10,%%zmm7,%%zmm7		\n\t	vpaddd	%%zmm17,%%zmm14,%%zmm14		\n\t	vpaddd	%%zmm24,%%zmm21,%%zmm21		\n\t"/* ...and add to fq0 to get q.lo32 */\
				"vpsrld	$06,%%zmm1,%%zmm1			\n\t	vpsrld	$06,%%zmm8,%%zmm8			\n\t	vpsrld	$06,%%zmm15,%%zmm15			\n\t	vpsrld	$06,%%zmm22,%%zmm22			\n\t"/* fq1>>6 leaves lo20 bits of eventual q.md32... */\
				"vpslld	$20,%%zmm2,%%zmm3			\n\t	vpslld	$20,%%zmm9,%%zmm10			\n\t	vpslld	$20,%%zmm16,%%zmm17			\n\t	vpslld	$20,%%zmm23,%%zmm24			\n\t"/* ...left-justify lo12 bits of fq2... */\
			"vpaddd	%%zmm3,%%zmm1,%%zmm1			\n\t	vpaddd	%%zmm10,%%zmm8,%%zmm8		\n\t	vpaddd	%%zmm17,%%zmm15,%%zmm15		\n\t	vpaddd	%%zmm24,%%zmm22,%%zmm22		\n\t"/* ...and add to get q.md32 */\
			"vpsrld	$12   ,%%zmm2,%%zmm2			\n\t	vpsrld	$12   ,%%zmm9,%%zmm9		\n\t	vpsrld	$12   ,%%zmm16,%%zmm16		\n\t	vpsrld	$12   ,%%zmm23,%%zmm23		\n\t"/* fq2>>12 gives q.hi32, which is at most 12-bits. */\
			/* Zero-extend each uint32 in ymm0-2 to a uint64, leaving q.lo32,md32,hi32 in qword-form in zmm0,1,2: */\
				"vpmovzxdq	%%ymm0,%%zmm0			\n\t	vpmovzxdq	%%ymm7,%%zmm7			\n\t	vpmovzxdq	%%ymm14,%%zmm14			\n\t	vpmovzxdq	%%ymm21,%%zmm21			\n\t"\
				"vpmovzxdq	%%ymm1,%%zmm1			\n\t	vpmovzxdq	%%ymm8,%%zmm8			\n\t	vpmovzxdq	%%ymm15,%%zmm15			\n\t	vpmovzxdq	%%ymm22,%%zmm22			\n\t"\
				"vpmovzxdq	%%ymm2,%%zmm2			\n\t	vpmovzxdq	%%ymm9,%%zmm9			\n\t	vpmovzxdq	%%ymm16,%%zmm16			\n\t	vpmovzxdq	%%ymm23,%%zmm23			\n\t"\
			/* Gather-load the needed byte-sized qinv initializers: */\
				"movq	%[__minv8],%%rbx			\n\t"\
				"movl	$-1,%%esi					\n\t"/* Mask-reg = 0x11...11; no longer use VPCMPEQD for this since in AVX-512
															that outputs one result bit per zmm-subfield-compare to a mask-reg */\
				"vpbroadcastd	%%esi,%%zmm3		\n\t	vpbroadcastd	%%esi,%%zmm10		\n\t	vpbroadcastd	%%esi,%%zmm17		\n\t	vpbroadcastd	%%esi,%%zmm24		\n\t"\
				"vpsrlq	$56,%%zmm3,%%zmm3			\n\t	vpsrlq	$56,%%zmm10,%%zmm10			\n\t	vpsrlq	$56,%%zmm17,%%zmm17			\n\t	vpsrlq	$56,%%zmm24,%%zmm24			\n\t"/* Each qword = 0x0000000000FF */\
				"vpandq	%%zmm3,%%zmm0,%%zmm4		\n\t	vpandq	%%zmm10,%%zmm7,%%zmm11		\n\t	vpandq	%%zmm17,%%zmm14,%%zmm18		\n\t	vpandq	%%zmm24,%%zmm21,%%zmm25		\n\t"/* AND with mask-reg leaves just low byte of each uint64 subfield */\
				"vpsrlq	$01,%%zmm4,%%zmm4			\n\t	vpsrlq	$01,%%zmm11,%%zmm11			\n\t	vpsrlq	$01,%%zmm18,%%zmm18			\n\t	vpsrlq	$01,%%zmm25,%%zmm25			\n\t"/* Lookup indices in minv byte-table */\
			/* Note VPGATHERDD may not use default-opmask k0.
			Re. inline-asm-with-opmask, this is what I first ran into when trying the 'obvious' {%%k1} syntax under both GCC and ICC,
			e.g. adding this line to an AVX-512 macro which compiles and runs fine sans any explicit opmask-register references:

				"vaddpd %%zmm3,%%zmm2,%%zmm1{%%k1} \n\t"	vaddpd %%zmm10,%%zmm9,%%zmm8{%%k2}	\n\t	vaddpd %%zmm17,%%zmm16,%%zmm15{%%k3}\n\t	vaddpd %%zmm24,%%zmm23,%%zmm22{%%k4}\n\t\

			gives "Assembler messages: Error: junk `%k1' after register". But, inside a basic (= non-extended) macro, with the
			escape-char-munge %% --> %, it compiles without issue. So it seems the culprit is GCC extended-asm. As far as a solution,
			thanks to Laurent Desnogues for turning up the following StackOverflow links.

			[1] http://stackoverflow.com/questions/21032395/masked-vector-instructions : Recommends double-curly-braces trickeration
				{{%%k1}} to make extended-asm play nice with AVX-512 opmask-register invocation. Alas, that hack is ICC-only; GCC gives
				"error: invalid 'asm': nested assembly dialect alternatives". That non-portability leads us to even-uglier hack number
			[2] http://stackoverflow.com/questions/34327831/invalid-asm-nested-assembly-dialect-alternatives : Explains that under GCC,
				double-curly-braces fail because
					"{} in GNU C inline asm already has a special meaning: providing alternatives for different
					ASM dialects. Using {{%%k1}} looks to gcc like nested alternatives, which is invalid."
				The GCC-workaround is the ugly %-escape-char-laden syntax %{%%k1%}, which, thankfully, works under both GCC and ICC.
				: */
				"movl	$-1,%%esi					\n\t"\
				"kmovw	%%esi,%%k1					\n\t	kmovw	%%esi,%%k2					\n\t	kmovw	%%esi,%%k3					\n\t	kmovw	%%esi,%%k4					\n\t"/* Init opmask k1 (Only need the low byte) */\
		"vpgatherqq (%%rbx,%%zmm4),%%zmm5%{%%k1%}	\n\tvpgatherqq (%%rbx,%%zmm11),%%zmm12%{%%k2%}\n\tvpgatherqq (%%rbx,%%zmm18),%%zmm19%{%%k3%}\n\tvpgatherqq (%%rbx,%%zmm25),%%zmm26%{%%k4%}\n\t"/* minv8 is a byte-table ... instruction sets mask-reg = 0, so each col uses separate same-valued k-reg */\
				"vpandq	%%zmm3,%%zmm5,%%zmm3		\n\t	vpandq	%%zmm10,%%zmm12,%%zmm10		\n\t	vpandq	%%zmm17,%%zmm19,%%zmm17		\n\t	vpandq	%%zmm24,%%zmm26,%%zmm24		\n\t"/* Iterative qinv-computation doesn't care about the garbage upper 7 bytes of each qword
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
				"vpmulld	%%zmm3,%%zmm4,%%zmm3	\n\t	vpmulld	%%zmm10,%%zmm11,%%zmm10		\n\t	vpmulld	%%zmm17,%%zmm18,%%zmm17		\n\t	vpmulld	%%zmm24,%%zmm25,%%zmm24		\n\t"/* qinv32 = MULL32(qinv16,(2-tmp32)); */\
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
			/* 78-bit truncation of qinv: */\
				"vpbroadcastd	%%esi,%%zmm31		\n\t"/* Mask-reg = 0x11...11 */\
				"vpsrlq	$50,%%zmm31,%%zmm31			\n\t"/* 0 all but bottom 14 bits of each qword */\
				"vpandq		%%zmm5,%%zmm31,%%zmm5	\n\t	vpandq		%%zmm12,%%zmm31,%%zmm12	\n\t	vpandq		%%zmm19,%%zmm31,%%zmm19	\n\t	vpandq		%%zmm26,%%zmm31,%%zmm26	\n\t"/* Hi 14 bits of 78-bit qinv. */\
			/* Reassemble into three 26-bit chunks, convert to double and write to local data store: */
				"vpbroadcastd	%%esi,%%zmm31		\n\t"/* Mask-reg = 0x11...11 */\
				"vpsrlq	$38,%%zmm31,%%zmm31			\n\t"/* 0 all but bottom 26 bits of each qword */\
				"vpsrlq	$26,%%zmm3,%%zmm6			\n\t	vpsrlq	$26,%%zmm10,%%zmm13			\n\t	vpsrlq	$26,%%zmm17,%%zmm20			\n\t	vpsrlq	$26,%%zmm24,%%zmm27			\n\t"/* hi6 bits of lo32 will go into md26... */\
				"vpsllq	$06,%%zmm4,%%zmm4			\n\t	vpsllq	$06,%%zmm11,%%zmm11			\n\t	vpsllq	$06,%%zmm18,%%zmm18			\n\t	vpsllq	$06,%%zmm25,%%zmm25			\n\t"/* md32<<6 can be done in place... */\
				"vpaddq		%%zmm6,%%zmm4,%%zmm4	\n\t	vpaddq		%%zmm13,%%zmm11,%%zmm11	\n\t	vpaddq		%%zmm20,%%zmm18,%%zmm18	\n\t	vpaddq		%%zmm27,%%zmm25,%%zmm25	\n\t"/* bits <26:63>, a 38-bit-wide temp. */\
				"vpandq		%%zmm3,%%zmm31,%%zmm3	\n\t	vpandq		%%zmm10,%%zmm31,%%zmm10	\n\t	vpandq		%%zmm17,%%zmm31,%%zmm17	\n\t	vpandq		%%zmm24,%%zmm31,%%zmm24	\n\t"/* qinv.lo26 */\
				"vpsrlq	$26,%%zmm4,%%zmm6			\n\t	vpsrlq	$26,%%zmm11,%%zmm13			\n\t	vpsrlq	$26,%%zmm18,%%zmm20			\n\t	vpsrlq	$26,%%zmm25,%%zmm27			\n\t"/* hi12 bits of md38 will go into hi26... */\
				"vpandq		%%zmm4,%%zmm31,%%zmm4	\n\t	vpandq		%%zmm11,%%zmm31,%%zmm11	\n\t	vpandq		%%zmm18,%%zmm31,%%zmm18	\n\t	vpandq		%%zmm25,%%zmm31,%%zmm25	\n\t"/* qinv.md26 */\
				"vpsllq	$12,%%zmm5,%%zmm5			\n\t	vpsllq	$12,%%zmm12,%%zmm12			\n\t	vpsllq	$12,%%zmm19,%%zmm19			\n\t	vpsllq	$12,%%zmm26,%%zmm26			\n\t"/* hi14<<12 can be done in place... */\
				"vpaddq		%%zmm6,%%zmm5,%%zmm5	\n\t	vpaddq		%%zmm13,%%zmm12,%%zmm12	\n\t	vpaddq		%%zmm20,%%zmm19,%%zmm19	\n\t	vpaddq		%%zmm27,%%zmm26,%%zmm26	\n\t"/* qinv.hi26 */\
			/* Down-convert each uint64 in zmm0-2 to a uint32, leaving qinv.lo32,md32,hi32 in dword-form in zmm3-5: */\
				"vpmovqd	%%zmm3,%%ymm3			\n\t	vpmovqd	%%zmm10,%%ymm10				\n\t	vpmovqd	%%zmm17,%%ymm17				\n\t	vpmovqd	%%zmm24,%%ymm24				\n\t"\
				"vpmovqd	%%zmm4,%%ymm4			\n\t	vpmovqd	%%zmm11,%%ymm11				\n\t	vpmovqd	%%zmm18,%%ymm18				\n\t	vpmovqd	%%zmm25,%%ymm25				\n\t"\
				"vpmovqd	%%zmm5,%%ymm5			\n\t	vpmovqd	%%zmm12,%%ymm12				\n\t	vpmovqd	%%zmm19,%%ymm19				\n\t	vpmovqd	%%zmm26,%%ymm26				\n\t"\
				"vcvtdq2pd	%%ymm3,%%zmm3			\n\t	vcvtdq2pd	%%ymm10,%%zmm10			\n\t	vcvtdq2pd	%%ymm17,%%zmm17			\n\t	vcvtdq2pd	%%ymm24,%%zmm24			\n\t"/* (double)fqinv0 */\
				"vcvtdq2pd	%%ymm4,%%zmm4			\n\t	vcvtdq2pd	%%ymm11,%%zmm11			\n\t	vcvtdq2pd	%%ymm18,%%zmm18			\n\t	vcvtdq2pd	%%ymm25,%%zmm25			\n\t"/* (double)fqinv1 */\
				"vcvtdq2pd	%%ymm5,%%zmm5			\n\t	vcvtdq2pd	%%ymm12,%%zmm12			\n\t	vcvtdq2pd	%%ymm19,%%zmm19			\n\t	vcvtdq2pd	%%ymm26,%%zmm26			\n\t"/* (double)fqinv2 */\
			/* Store base-2^26 double-based qinv: */\
				"vmovaps	%%zmm3,0x800(%%rax)		\n\t	vmovaps	%%zmm10,0x840(%%rax)		\n\t	vmovaps	%%zmm17,0x880(%%rax)		\n\t	vmovaps	%%zmm24,0x8c0(%%rax)		\n\t"/* store lo26 into fqinv0 */\
				"vmovaps	%%zmm4,0xa00(%%rax)		\n\t	vmovaps	%%zmm11,0xa40(%%rax)		\n\t	vmovaps	%%zmm18,0xa80(%%rax)		\n\t	vmovaps	%%zmm25,0xac0(%%rax)		\n\t"/* store lo52 into fqinv1 */\
				"vmovaps	%%zmm5,0xc00(%%rax)		\n\t	vmovaps	%%zmm12,0xc40(%%rax)		\n\t	vmovaps	%%zmm19,0xc80(%%rax)		\n\t	vmovaps	%%zmm26,0xcc0(%%rax)		\n\t"/* store hi26 into fqinv2 */\
			/****************************************************************************************/\
			/*                 And now process 2nd quartet of vector data:                          */\
			/****************************************************************************************/\
				"movq	%[__fq0] ,%%rax				\n\t"\
				"movq	%[__p2]  ,%%rbx				\n\t"\
				"movq	%[__kvec],%%rcx				\n\t"/* Each column uses 7 distinct vector registers, so increment all non-shared vec-reg indices by +7 for each additional rightward column: */\
				"vbroadcastsd		(%%rbx),%%zmm1	\n\t"/* broadcast 2*p to all 8 qwords of zmm ... double-prec version is from mem-address */\
				"vmovupd	0x100(%%rcx),%%zmm3		\n\t	vmovupd	0x140(%%rcx),%%zmm10		\n\t	vmovupd	0x180(%%rcx),%%zmm17		\n\t	vmovupd	0x1c0(%%rcx),%%zmm24		\n\t"/* load the factor-candidate k's, 8 per zmm ... note kvec not nec. 16-byte aligned */\
				"vmulpd		%%zmm1,%%zmm3,%%zmm2	\n\t	vmulpd		%%zmm1,%%zmm10,%%zmm9	\n\t	vmulpd		%%zmm1,%%zmm17,%%zmm16	\n\t	vmulpd		%%zmm1,%%zmm24,%%zmm23	\n\t"/* hi = x*y = 2p*k */\
				"vmovaps	%%zmm2,%%zmm0			\n\t	vmovaps	%%zmm9,%%zmm7				\n\t	vmovaps	%%zmm16,%%zmm14				\n\t	vmovaps	%%zmm23,%%zmm21				\n\t"/* cpy hi into lo-destination reg */\
			"vfmsub231pd	%%zmm1,%%zmm3,%%zmm0	\n\t vfmsub231pd	%%zmm1,%%zmm10,%%zmm7	\n\t vfmsub231pd	%%zmm1,%%zmm17,%%zmm14	\n\t vfmsub231pd	%%zmm1,%%zmm24,%%zmm21	\n\t"/* lo = fma(x,y,-hi; xmm1,3 FREE */\
				"movl	$1,%%esi					\n\t"\
				"vpbroadcastd	%%esi,%%zmm3		\n\t	vpbroadcastd	%%esi,%%zmm10		\n\t	vpbroadcastd	%%esi,%%zmm17		\n\t	vpbroadcastd	%%esi,%%zmm24		\n\t"/* Broadcast to all 8 int32 slots of zmm3, since KNL lacks AVX-512 VL extensions */\
				"vcvtdq2pd	%%ymm3,%%zmm3			\n\t	vcvtdq2pd	%%ymm10,%%zmm10			\n\t	vcvtdq2pd	%%ymm17,%%zmm17			\n\t	vcvtdq2pd	%%ymm24,%%zmm24			\n\t"/* zmm3 = (double)1, conversion from 4-way int32 in lower half of zmm3 */\
				"vaddpd		%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd		%%zmm10,%%zmm7,%%zmm7	\n\t	vaddpd		%%zmm17,%%zmm14,%%zmm14	\n\t	vaddpd		%%zmm24,%%zmm21,%%zmm21	\n\t"/* lo += 1 */\
			/* Normalize to properly split the resultant significant bits across lo,hi, use base-2^52 signed-int-stored-as-double: */\
				"movq	%[__two26f],%%rbx			\n\t"/* Several of these const-regs got overwritten in above qinv computation, so re-init */\
			"vbroadcastsd	0xc0(%%rbx),%%zmm31		\n\t	vbroadcastsd	0x80(%%rbx),%%zmm30	\n\t	vbroadcastsd	0x40(%%rbx),%%zmm29	\n\t	vbroadcastsd		(%%rbx),%%zmm28	\n\t"/* two52i,two52f,two26i,two26f */\
				"vmovaps	%%zmm2,%%zmm1			\n\t	vmovaps	%%zmm9,%%zmm8				\n\t	vmovaps	%%zmm16,%%zmm15				\n\t	vmovaps	%%zmm23,%%zmm22				\n\t"/* cpy hi into hh-destination reg */\
				"vmulpd		%%zmm31,%%zmm1,%%zmm1	\n\t	vmulpd		%%zmm31,%%zmm8,%%zmm8	\n\t	vmulpd		%%zmm31,%%zmm15,%%zmm15	\n\t	vmulpd		%%zmm31,%%zmm22,%%zmm22	\n\t"/*            hi*TWO52FLINV  */\
				"vrndscalepd	$1,%%zmm1,%%zmm1	\n\t	vrndscalepd	$1,%%zmm8,%%zmm8		\n\t	vrndscalepd	$1,%%zmm15,%%zmm15		\n\t	vrndscalepd	$1,%%zmm22,%%zmm22		\n\t"/* hh = FLOOR(hi*TWO52FLINV) */\
				"vmovaps	%%zmm2,%%zmm3			\n\t	vmovaps	%%zmm9,%%zmm10				\n\t	vmovaps	%%zmm16,%%zmm17				\n\t	vmovaps	%%zmm23,%%zmm24				\n\t"/* cpy hi into cy-destination reg */\
			"vfnmadd231pd	%%zmm30,%%zmm1,%%zmm3	\n\t vfnmadd231pd	%%zmm30,%%zmm8,%%zmm10	\n\t vfnmadd231pd	%%zmm30,%%zmm15,%%zmm17	\n\t vfnmadd231pd	%%zmm30,%%zmm22,%%zmm24\n\t"/* cy = fma(hh,-TWO52FLOAT,hi)	Backward carry from hi into lo */\
				"vaddpd		%%zmm3,%%zmm0,%%zmm0	\n\t	vaddpd		%%zmm10,%%zmm7,%%zmm7	\n\t	vaddpd		%%zmm17,%%zmm14,%%zmm14	\n\t	vaddpd		%%zmm24,%%zmm21,%%zmm21	\n\t"/* lo += cy */\
				"vmovaps	%%zmm1,0x500(%%rax)		\n\t	vmovaps	%%zmm8,0x540(%%rax)			\n\t	vmovaps	%%zmm15,0x580(%%rax)		\n\t	vmovaps	%%zmm22,0x5c0(%%rax)		\n\t"/* store hi26 into fq2 */\
			/* Normalize to convert 78-bit q from balanced-digit [lo52,hi26] representation to nonnegative-digit [lo26,hi52]: */\
				"vmovaps	%%zmm0,%%zmm2			\n\t	vmovaps	%%zmm7,%%zmm9				\n\t	vmovaps	%%zmm14,%%zmm16				\n\t	vmovaps	%%zmm21,%%zmm23				\n\t"/* cpy lo into destination reg */\
				"vmulpd		%%zmm29,%%zmm2,%%zmm2	\n\t	vmulpd		%%zmm29,%%zmm9,%%zmm9	\n\t	vmulpd		%%zmm29,%%zmm16,%%zmm16	\n\t	vmulpd		%%zmm29,%%zmm23,%%zmm23	\n\t"/*            lo*TWO26FLINV  */\
				"vrndscalepd	$1,%%zmm2,%%zmm2	\n\t	vrndscalepd	$1,%%zmm9,%%zmm9		\n\t	vrndscalepd	$1,%%zmm16,%%zmm16		\n\t	vrndscalepd	$1,%%zmm23,%%zmm23		\n\t"/* cy = FLOOR(lo*TWO26FLINV)	Forward carry from lo into hi */\
				"vmovaps	%%zmm2,0x300(%%rax)		\n\t	vmovaps	%%zmm9,0x340(%%rax)			\n\t	vmovaps	%%zmm16,0x380(%%rax)		\n\t	vmovaps	%%zmm23,0x3c0(%%rax)		\n\t"/* store cy into fq1 */\
			"vfnmadd231pd	%%zmm28,%%zmm2,%%zmm0	\n\t vfnmadd231pd	%%zmm28,%%zmm9,%%zmm7	\n\t vfnmadd231pd	%%zmm28,%%zmm16,%%zmm14	\n\t vfnmadd231pd	%%zmm28,%%zmm23,%%zmm21	\n\t"/* lo26 = fma(cy,-TWO26FLOAT,lo) */\
			" vfmadd132pd	%%zmm28,%%zmm2,%%zmm1	\n\t vfmadd132pd	%%zmm28,%%zmm9,%%zmm8	\n\t vfmadd132pd	%%zmm28,%%zmm16,%%zmm15	\n\t vfmadd132pd	%%zmm28,%%zmm23,%%zmm22	\n\t"/* hi52 = fma(hi,+TWO26FLOAT,cy) */\
				"vmovaps	%%zmm0,0x100(%%rax)		\n\t	vmovaps	%%zmm7,0x140(%%rax)			\n\t	vmovaps	%%zmm14,0x180(%%rax)		\n\t	vmovaps	%%zmm21,0x1c0(%%rax)		\n\t"/* store lo26 into fq0 */\
				"vmovaps	%%zmm1,0x700(%%rax)		\n\t	vmovaps	%%zmm8,0x740(%%rax)			\n\t	vmovaps	%%zmm15,0x780(%%rax)		\n\t	vmovaps	%%zmm22,0x7c0(%%rax)		\n\t"/* store hi52 into fqhi52 */\
			/*** hi = hh ==> lo26 in zmm0, hi52 = zmm1. ***/\
			/* Assemble 32-bit pieces of q from the 26-bit floating-point ones ... AVX-512F has just double->int32, so
			convert our three 26-bit double-vec chunks to 8-way int32s (i.e. zmm -> ymm), then re-expand result to 8-way
			int64s (with high 32-bit halves all zero) in preparation for mod-inverse computation: */\
				"vcvttpd2dq		  %%zmm0,%%ymm0		\n\t	vcvttpd2dq		  %%zmm7,%%ymm7		\n\t	vcvttpd2dq		  %%zmm14,%%ymm14	\n\t	vcvttpd2dq		  %%zmm21,%%ymm21	\n\t"/* (int32)fq0 */\
				"vcvttpd2dq	0x300(%%rax),%%ymm1		\n\t	vcvttpd2dq	0x340(%%rax),%%ymm8		\n\t	vcvttpd2dq	0x380(%%rax),%%ymm15	\n\t	vcvttpd2dq	0x3c0(%%rax),%%ymm22	\n\t"/* (int32)fq1 */\
				"vcvttpd2dq	0x500(%%rax),%%ymm2		\n\t	vcvttpd2dq	0x540(%%rax),%%ymm9		\n\t	vcvttpd2dq	0x580(%%rax),%%ymm16	\n\t	vcvttpd2dq	0x5c0(%%rax),%%ymm23	\n\t"/* (int32)fq2 */\
		/*** Must use full-length zmm-regs here (even though upper 256 bits unused) due to AVX-512F full-width constraint (applying to ymm needs AVX-512VL extensions): ***/\
				"vpslld	$26,%%zmm1,%%zmm3			\n\t	vpslld	$26,%%zmm8,%%zmm10			\n\t	vpslld	$26,%%zmm15,%%zmm17			\n\t	vpslld	$26,%%zmm22,%%zmm24			\n\t"/* left-justify lo6 bits of fq1... */\
			"vpaddd	%%zmm3,%%zmm0,%%zmm0			\n\t	vpaddd	%%zmm10,%%zmm7,%%zmm7		\n\t	vpaddd	%%zmm17,%%zmm14,%%zmm14		\n\t	vpaddd	%%zmm24,%%zmm21,%%zmm21		\n\t"/* ...and add to fq0 to get q.lo32 */\
				"vpsrld	$06,%%zmm1,%%zmm1			\n\t	vpsrld	$06,%%zmm8,%%zmm8			\n\t	vpsrld	$06,%%zmm15,%%zmm15			\n\t	vpsrld	$06,%%zmm22,%%zmm22			\n\t"/* fq1>>6 leaves lo20 bits of eventual q.md32... */\
				"vpslld	$20,%%zmm2,%%zmm3			\n\t	vpslld	$20,%%zmm9,%%zmm10			\n\t	vpslld	$20,%%zmm16,%%zmm17			\n\t	vpslld	$20,%%zmm23,%%zmm24			\n\t"/* ...left-justify lo12 bits of fq2... */\
			"vpaddd	%%zmm3,%%zmm1,%%zmm1			\n\t	vpaddd	%%zmm10,%%zmm8,%%zmm8		\n\t	vpaddd	%%zmm17,%%zmm15,%%zmm15		\n\t	vpaddd	%%zmm24,%%zmm22,%%zmm22		\n\t"/* ...and add to get q.md32 */\
			"vpsrld	$12   ,%%zmm2,%%zmm2			\n\t	vpsrld	$12   ,%%zmm9,%%zmm9		\n\t	vpsrld	$12   ,%%zmm16,%%zmm16		\n\t	vpsrld	$12   ,%%zmm23,%%zmm23		\n\t"/* fq2>>12 gives q.hi32, which is at most 12-bits. */\
			/* Zero-extend each uint32 in ymm0-2 to a uint64, leaving q.lo32,md32,hi32 in qword-form in zmm0,1,2: */\
				"vpmovzxdq	%%ymm0,%%zmm0			\n\t	vpmovzxdq	%%ymm7,%%zmm7			\n\t	vpmovzxdq	%%ymm14,%%zmm14			\n\t	vpmovzxdq	%%ymm21,%%zmm21			\n\t"\
				"vpmovzxdq	%%ymm1,%%zmm1			\n\t	vpmovzxdq	%%ymm8,%%zmm8			\n\t	vpmovzxdq	%%ymm15,%%zmm15			\n\t	vpmovzxdq	%%ymm22,%%zmm22			\n\t"\
				"vpmovzxdq	%%ymm2,%%zmm2			\n\t	vpmovzxdq	%%ymm9,%%zmm9			\n\t	vpmovzxdq	%%ymm16,%%zmm16			\n\t	vpmovzxdq	%%ymm23,%%zmm23			\n\t"\
			/* Gather-load the needed byte-sized qinv initializers: */\
				"movq	%[__minv8],%%rbx			\n\t"\
				"movl	$-1,%%esi					\n\t"/* Mask-reg = 0x11...11 */\
				"vpbroadcastd	%%esi,%%zmm3		\n\t	vpbroadcastd	%%esi,%%zmm10		\n\t	vpbroadcastd	%%esi,%%zmm17		\n\t	vpbroadcastd	%%esi,%%zmm24		\n\t"\
				"vpsrlq	$56,%%zmm3,%%zmm3			\n\t	vpsrlq	$56,%%zmm10,%%zmm10			\n\t	vpsrlq	$56,%%zmm17,%%zmm17			\n\t	vpsrlq	$56,%%zmm24,%%zmm24			\n\t"/* Each qword = 0x0000000000FF */\
				"vpandq	%%zmm3,%%zmm0,%%zmm4		\n\t	vpandq	%%zmm10,%%zmm7,%%zmm11		\n\t	vpandq	%%zmm17,%%zmm14,%%zmm18		\n\t	vpandq	%%zmm24,%%zmm21,%%zmm25		\n\t"/* AND with mask-reg leaves just low byte of each uint64 subfield */\
				"vpsrlq	$01,%%zmm4,%%zmm4			\n\t	vpsrlq	$01,%%zmm11,%%zmm11			\n\t	vpsrlq	$01,%%zmm18,%%zmm18			\n\t	vpsrlq	$01,%%zmm25,%%zmm25			\n\t"/* Lookup indices in minv byte-table */\
			/* Note VPGATHERDD may not use default-opmask k0: */
				"movl	$-1,%%esi					\n\t"\
				"kmovw	%%esi,%%k1					\n\t	kmovw	%%esi,%%k2					\n\t	kmovw	%%esi,%%k3					\n\t	kmovw	%%esi,%%k4					\n\t"/* Init opmask k1 (Only need the low byte) */\
		"vpgatherqq (%%rbx,%%zmm4),%%zmm5%{%%k1%}	\n\tvpgatherqq (%%rbx,%%zmm11),%%zmm12%{%%k2%}\n\tvpgatherqq (%%rbx,%%zmm18),%%zmm19%{%%k3%}\n\tvpgatherqq (%%rbx,%%zmm25),%%zmm26%{%%k4%}\n\t"/* minv8 is a byte-table ... instruction sets mask-reg = 0, so each col uses separate same-valued k-reg */\
				"vpandq	%%zmm3,%%zmm5,%%zmm3		\n\t	vpandq	%%zmm10,%%zmm12,%%zmm10		\n\t	vpandq	%%zmm17,%%zmm19,%%zmm17		\n\t	vpandq	%%zmm24,%%zmm26,%%zmm24		\n\t"/* Iterative qinv-computation doesn't care about the garbage upper 7 bytes of each qword
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
				"vpmulld	%%zmm3,%%zmm4,%%zmm3	\n\t	vpmulld	%%zmm10,%%zmm11,%%zmm10		\n\t	vpmulld	%%zmm17,%%zmm18,%%zmm17		\n\t	vpmulld	%%zmm24,%%zmm25,%%zmm24		\n\t"/* qinv32 = MULL32(qinv16,(2-tmp32)); */\
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
			/* 78-bit truncation of qinv: */\
				"vpbroadcastd	%%esi,%%zmm31		\n\t"/* Mask-reg = 0x11...11 */\
				"vpsrlq	$50,%%zmm31,%%zmm31			\n\t"/* 0 all but bottom 14 bits of each qword */\
				"vpandq		%%zmm5,%%zmm31,%%zmm5	\n\t	vpandq		%%zmm12,%%zmm31,%%zmm12	\n\t	vpandq		%%zmm19,%%zmm31,%%zmm19	\n\t	vpandq		%%zmm26,%%zmm31,%%zmm26	\n\t"/* Hi 14 bits of 78-bit qinv. */\
			/* Reassemble into three 26-bit chunks, convert to double and write to local data store: */
				"vpbroadcastd	%%esi,%%zmm31		\n\t"/* Mask-reg = 0x11...11 */\
				"vpsrlq	$38,%%zmm31,%%zmm31			\n\t"/* 0 all but bottom 26 bits of each qword */\
				"vpsrlq	$26,%%zmm3,%%zmm6			\n\t	vpsrlq	$26,%%zmm10,%%zmm13			\n\t	vpsrlq	$26,%%zmm17,%%zmm20			\n\t	vpsrlq	$26,%%zmm24,%%zmm27			\n\t"/* hi6 bits of lo32 will go into md26... */\
				"vpsllq	$06,%%zmm4,%%zmm4			\n\t	vpsllq	$06,%%zmm11,%%zmm11			\n\t	vpsllq	$06,%%zmm18,%%zmm18			\n\t	vpsllq	$06,%%zmm25,%%zmm25			\n\t"/* md32<<6 can be done in place... */\
				"vpaddq		%%zmm6,%%zmm4,%%zmm4	\n\t	vpaddq		%%zmm13,%%zmm11,%%zmm11	\n\t	vpaddq		%%zmm20,%%zmm18,%%zmm18	\n\t	vpaddq		%%zmm27,%%zmm25,%%zmm25	\n\t"/* bits <26:63>, a 38-bit-wide temp. */\
				"vpandq		%%zmm3,%%zmm31,%%zmm3	\n\t	vpandq		%%zmm10,%%zmm31,%%zmm10	\n\t	vpandq		%%zmm17,%%zmm31,%%zmm17	\n\t	vpandq		%%zmm24,%%zmm31,%%zmm24	\n\t"/* qinv.lo26 */\
				"vpsrlq	$26,%%zmm4,%%zmm6			\n\t	vpsrlq	$26,%%zmm11,%%zmm13			\n\t	vpsrlq	$26,%%zmm18,%%zmm20			\n\t	vpsrlq	$26,%%zmm25,%%zmm27			\n\t"/* hi12 bits of md38 will go into hi26... */\
				"vpandq		%%zmm4,%%zmm31,%%zmm4	\n\t	vpandq		%%zmm11,%%zmm31,%%zmm11	\n\t	vpandq		%%zmm18,%%zmm31,%%zmm18	\n\t	vpandq		%%zmm25,%%zmm31,%%zmm25	\n\t"/* qinv.md26 */\
				"vpsllq	$12,%%zmm5,%%zmm5			\n\t	vpsllq	$12,%%zmm12,%%zmm12			\n\t	vpsllq	$12,%%zmm19,%%zmm19			\n\t	vpsllq	$12,%%zmm26,%%zmm26			\n\t"/* hi14<<12 can be done in place... */\
				"vpaddq		%%zmm6,%%zmm5,%%zmm5	\n\t	vpaddq		%%zmm13,%%zmm12,%%zmm12	\n\t	vpaddq		%%zmm20,%%zmm19,%%zmm19	\n\t	vpaddq		%%zmm27,%%zmm26,%%zmm26	\n\t"/* qinv.hi26 */\
			/* Down-convert each uint64 in zmm0-2 to a uint32, leaving qinv.lo32,md32,hi32 in dword-form in zmm3-5: */\
				"vpmovqd	%%zmm3,%%ymm3			\n\t	vpmovqd	%%zmm10,%%ymm10				\n\t	vpmovqd	%%zmm17,%%ymm17				\n\t	vpmovqd	%%zmm24,%%ymm24				\n\t"\
				"vpmovqd	%%zmm4,%%ymm4			\n\t	vpmovqd	%%zmm11,%%ymm11				\n\t	vpmovqd	%%zmm18,%%ymm18				\n\t	vpmovqd	%%zmm25,%%ymm25				\n\t"\
				"vpmovqd	%%zmm5,%%ymm5			\n\t	vpmovqd	%%zmm12,%%ymm12				\n\t	vpmovqd	%%zmm19,%%ymm19				\n\t	vpmovqd	%%zmm26,%%ymm26				\n\t"\
				"vcvtdq2pd	%%ymm3,%%zmm3			\n\t	vcvtdq2pd	%%ymm10,%%zmm10			\n\t	vcvtdq2pd	%%ymm17,%%zmm17			\n\t	vcvtdq2pd	%%ymm24,%%zmm24			\n\t"/* (double)fqinv0 */\
				"vcvtdq2pd	%%ymm4,%%zmm4			\n\t	vcvtdq2pd	%%ymm11,%%zmm11			\n\t	vcvtdq2pd	%%ymm18,%%zmm18			\n\t	vcvtdq2pd	%%ymm25,%%zmm25			\n\t"/* (double)fqinv1 */\
				"vcvtdq2pd	%%ymm5,%%zmm5			\n\t	vcvtdq2pd	%%ymm12,%%zmm12			\n\t	vcvtdq2pd	%%ymm19,%%zmm19			\n\t	vcvtdq2pd	%%ymm26,%%zmm26			\n\t"/* (double)fqinv2 */\
			/* Store base-2^26 double-based qinv: */\
				"vmovaps	%%zmm3,0x900(%%rax)		\n\t	vmovaps	%%zmm10,0x940(%%rax)		\n\t	vmovaps	%%zmm17,0x980(%%rax)		\n\t	vmovaps	%%zmm24,0x9c0(%%rax)		\n\t"/* store lo26 into fqinv0 */\
				"vmovaps	%%zmm4,0xb00(%%rax)		\n\t	vmovaps	%%zmm11,0xb40(%%rax)		\n\t	vmovaps	%%zmm18,0xb80(%%rax)		\n\t	vmovaps	%%zmm25,0xbc0(%%rax)		\n\t"/* store lo52 into fqinv1 */\
				"vmovaps	%%zmm5,0xd00(%%rax)		\n\t	vmovaps	%%zmm12,0xd40(%%rax)		\n\t	vmovaps	%%zmm19,0xd80(%%rax)		\n\t	vmovaps	%%zmm26,0xdc0(%%rax)		\n\t"/* store hi26 into fqinv2 */\
				:					/* outputs: none */\
				: [__fq0]  "m" (tmp)	/* All inputs from memory addresses here */\
				 ,[__p2]   "m" (dptr)\
				 ,[__kvec] "m" (tm1)	/* vec-of-doubles copy of input kvec */\
				 ,[__minv8]   "m" (minv8_ptr)	/* Ptr to Table of precomputed byte-inverses def'd in mi64.h */\
				 ,[__two26f] "m" (two26f)\
				: "cc","memory","rax","rbx","rcx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15","xmm16","xmm17","xmm18","xmm19","xmm20","xmm21","xmm22","xmm23","xmm24","xmm25","xmm26","xmm27","xmm28","xmm29","xmm30","xmm31"	/* Clobbered registers */\
			);

		#endif

			// In AVX2|512+FMA version, overwrite middle term of usual qinv[0,1,2] sequence with low 52 bits of qinv:
			dptr = fqinv0[0];
			__asm__ volatile (\
				"movq	%[__fqinv0],%%rdx	\n\t"\
				"movq	%[__two26f],%%rsi	\n\t"\
				"vmovaps	(%%rsi),%%zmm3	\n\t"/* TWO26FLOAT */\
				"vmovaps	0x000(%%rdx),%%zmm0		\n\t	vmovaps	0x040(%%rdx),%%zmm4		\n\t	vmovaps	0x080(%%rdx),%%zmm8		\n\t	vmovaps	0x0c0(%%rdx),%%zmm12	\n\t"\
				"vmovaps	0x200(%%rdx),%%zmm1		\n\t	vmovaps	0x240(%%rdx),%%zmm5		\n\t	vmovaps	0x280(%%rdx),%%zmm9		\n\t	vmovaps	0x2c0(%%rdx),%%zmm13	\n\t"\
				"vmovaps	%%zmm1,%%zmm2			\n\t	vmovaps	%%zmm5,%%zmm6			\n\t	vmovaps	%%zmm9,%%zmm10			\n\t	vmovaps	%%zmm13,%%zmm14			\n\t"/* cpy fhi into flo-destination reg */\
			"vfmadd132pd %%zmm3,%%zmm0,%%zmm2		\n\t vfmadd132pd %%zmm3,%%zmm4,%%zmm6	\n\t vfmadd132pd %%zmm3,%%zmm8,%%zmm10	\n\t vfmadd132pd %%zmm3,%%zmm12,%%zmm14	\n\t"/* (qinv0 + 2^26*qinv1)*/\
				"vmovaps	%%zmm2,0x200(%%rdx)		\n\t	vmovaps	%%zmm6,0x240(%%rdx)		\n\t	vmovaps	%%zmm10,0x280(%%rdx)	\n\t	vmovaps	%%zmm14,0x2c0(%%rdx)	\n\t"\
			/*** 2nd quartet of vector data: ***/
				"vmovaps	0x100(%%rdx),%%zmm0		\n\t	vmovaps	0x140(%%rdx),%%zmm4		\n\t	vmovaps	0x180(%%rdx),%%zmm8		\n\t	vmovaps	0x1c0(%%rdx),%%zmm12	\n\t"\
				"vmovaps	0x300(%%rdx),%%zmm1		\n\t	vmovaps	0x340(%%rdx),%%zmm5		\n\t	vmovaps	0x380(%%rdx),%%zmm9		\n\t	vmovaps	0x3c0(%%rdx),%%zmm13	\n\t"\
				"vmovaps	%%zmm1,%%zmm2			\n\t	vmovaps	%%zmm5,%%zmm6			\n\t	vmovaps	%%zmm9,%%zmm10			\n\t	vmovaps	%%zmm13,%%zmm14			\n\t"/* cpy fhi into flo-destination reg */\
			"vfmadd132pd %%zmm3,%%zmm0,%%zmm2		\n\t vfmadd132pd %%zmm3,%%zmm4,%%zmm6	\n\t vfmadd132pd %%zmm3,%%zmm8,%%zmm10	\n\t vfmadd132pd %%zmm3,%%zmm12,%%zmm14	\n\t"/* (qinv0 + 2^26*qinv1)*/\
				"vmovaps	%%zmm2,0x300(%%rdx)		\n\t	vmovaps	%%zmm6,0x340(%%rdx)		\n\t	vmovaps	%%zmm10,0x380(%%rdx)	\n\t	vmovaps	%%zmm14,0x3c0(%%rdx)	\n\t"\
				:					/* outputs: none */\
				: [__fqinv0] "m" (dptr)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)	\
				: "cc","memory","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6", "xmm8","xmm9","xmm10", "xmm12","xmm13","xmm14"	/* Clobbered registers */\
			);

			// zshift is in [0,76]:
			if(zshift < 26) {
				dtmp = 1<<zshift;		for(j = 0; j < 64; j += 8) { VEC_DBL_INIT_8((vec_dbl*)fx0[j],dtmp); VEC_DBL_INIT_8((vec_dbl*)fx1[j], 0.0); VEC_DBL_INIT_8((vec_dbl*)fx2[j], 0.0); }
			} else if(zshift < 52) {
				dtmp = 1<<(zshift-26);	for(j = 0; j < 64; j += 8) { VEC_DBL_INIT_8((vec_dbl*)fx0[j], 0.0); VEC_DBL_INIT_8((vec_dbl*)fx1[j],dtmp); VEC_DBL_INIT_8((vec_dbl*)fx2[j], 0.0); }
			} else if(zshift < 78) {
				dtmp = 1<<(zshift-52);	for(j = 0; j < 64; j += 8) { VEC_DBL_INIT_8((vec_dbl*)fx0[j], 0.0); VEC_DBL_INIT_8((vec_dbl*)fx1[j], 0.0); VEC_DBL_INIT_8((vec_dbl*)fx2[j],dtmp); }
			} else {
				ASSERT(0,"zshift out of range!");
			}

			/*...x^2 mod q is returned in x. */
			// All 3-word-uint64-form operands have components in [0, 2^26).

		/* Inner loop body needs 42 MOVAPS, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */
		/* !ASM_LOOP: High-level loop construct, and explicit branch to do the modular doubling: */
			for(j = start_index-1; j >= 0; j--) {
				SSE2_twopmodq78_modmul_q64(fq0[0],fqinv0[0],fx0[0],two26i,pshift,j);
			}	/* for(j...) */

			/* Do a final mod-doubling and return. Code here is specialized for the case
			where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2 or -(q-1)/2:
			*/
			dptr = fq0[0];
			__asm__ volatile (\
				"movq	%[__fq0]   ,%%rax	\n\t"\
				"movq	%[__two26f],%%rbx	\n\t"\
				"vmovaps	0x40(%%rbx),%%zmm3		\n\t	vmovaps	(%%rbx),%%zmm7			\n\t"/* two26i,f */\
			"vfmadd213pd %%zmm1,%%zmm7,%%zmm2		\n\t vfmadd213pd %%zmm5,%%zmm7,%%zmm6	\n\t vfmadd213pd %%zmm9,%%zmm7,%%zmm10	\n\t vfmadd213pd %%zmm13,%%zmm7,%%zmm14	\n\t"/* (x1 + 2^26*x2) = hi52 */\
				/* Mod-doubling code taken from SSE2_twopmodq78_modmul_q16, expects x.lo26,hi52 in zmm1,2, resp., so copy zmm0,4,8,12 -> zmm1,5,9,13 after computing x.hi52 (in zmm2,6,10,14): */\
				"vaddpd	%%zmm0,%%zmm0,%%zmm1		\n\t	vaddpd	%%zmm4,%%zmm4,%%zmm5	\n\t	vaddpd	%%zmm8 ,%%zmm8 ,%%zmm9	\n\t	vaddpd	%%zmm12,%%zmm12,%%zmm13	\n\t"/* low 26 bits */\
				"vaddpd	%%zmm2,%%zmm2,%%zmm2		\n\t	vaddpd	%%zmm6,%%zmm6,%%zmm6	\n\t	vaddpd	%%zmm10,%%zmm10,%%zmm10	\n\t	vaddpd	%%zmm14,%%zmm14,%%zmm14	\n\t"/* top 52 bits */\
				/* And subtract q.lo26,hi52 - note q.lo26 is balanced-digit normalized and divided by 2^26, so do this sub prior to final >0 normalization below: */\
				"vsubpd	     (%%rax),%%zmm1,%%zmm1	\n\t vsubpd	0x040(%%rax),%%zmm5,%%zmm5	\n\t vsubpd	0x080(%%rax),%%zmm9 ,%%zmm9	\n\t vsubpd	0x0c0(%%rax),%%zmm13,%%zmm13\n\t"/* -= qlo26 */\
				"vsubpd	0x600(%%rax),%%zmm2,%%zmm2	\n\t vsubpd	0x640(%%rax),%%zmm6,%%zmm6	\n\t vsubpd	0x680(%%rax),%%zmm10,%%zmm10\n\t vsubpd	0x6c0(%%rax),%%zmm14,%%zmm14\n\t"/* -= qhi52 */\
				/* Final-residue computation needs digits >= 0 normalization - since fq0 will get folded in with fq2 anyway, only need renormalize fq0: */\
				"vmulpd	%%zmm3,%%zmm1,%%zmm1		\n\t	vmulpd	%%zmm3,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm3,%%zmm9,%%zmm9	\n\t	vmulpd	%%zmm3,%%zmm13,%%zmm13	\n\t"/* xlo *= two26i */\
				"vrndscalepd $1,%%zmm1,%%zmm0		\n\t	vrndscalepd $1,%%zmm5,%%zmm4	\n\t	vrndscalepd $1,%%zmm9,%%zmm8	\n\t	vrndscalepd $1,%%zmm13,%%zmm12	\n\t"/* fcy */\
				"vsubpd	%%zmm0,%%zmm1,%%zmm1		\n\t	vsubpd	%%zmm4,%%zmm5,%%zmm5	\n\t	vsubpd	%%zmm8,%%zmm9,%%zmm9	\n\t	vsubpd	%%zmm12,%%zmm13,%%zmm13	\n\t"/* fx0 -= fcy */\
				"vaddpd	%%zmm0,%%zmm2,%%zmm2		\n\t	vaddpd	%%zmm4,%%zmm6,%%zmm6	\n\t	vaddpd	%%zmm8,%%zmm10,%%zmm10	\n\t	vaddpd	%%zmm12,%%zmm14,%%zmm14	\n\t"/* Add carry into xhi */\
				"vmulpd	%%zmm7,%%zmm1,%%zmm1		\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm7,%%zmm9,%%zmm9	\n\t	vmulpd	%%zmm7 ,%%zmm13,%%zmm13	\n\t"/* fx0 *= two26f ... result must go into m0 */\
				/* (lo26 == 1 && hi52 == 0), save result in a bitmask: */
				"vmovaps	(%%rbx),%%zmm0			\n\t"/* TWO26FLOAT */\
				"vmulpd	0x40(%%rbx),%%zmm0,%%zmm3	\n\t"/* 1.0 ... no good way to generate 1.0 on-the-fly, so use 2^26 * 2^-26 */\
				"vpxorq	     %%zmm0,%%zmm0,%%zmm0	\n\t"/* 0.0 */\
				"vcmppd $0,%%zmm1,%%zmm3,%%k0		\n\t	vcmppd $0,%%zmm5,%%zmm3,%%k1	\n\t	vcmppd $0,%%zmm9,%%zmm3,%%k2	\n\t	vcmppd $0,%%zmm13,%%zmm3,%%k3	\n\t"/* low 26 bits == 1 ? In AVX512 version, mask-regs k0-3 replace dest-regs m1,5,9,13 */\
				"vcmppd $0,%%zmm2,%%zmm0,%%k4		\n\t	vcmppd $0,%%zmm6,%%zmm0,%%k5	\n\t	vcmppd $0,%%zmm10,%%zmm0,%%k6	\n\t	vcmppd $0,%%zmm14,%%zmm0,%%k7	\n\t"/* top 52 bits == 0 ? In AVX512 version, mask-regs k4-7 replace dest-regs m0,4,8,12 */\
				"kandw	%%k0,%%k4,%%k0				\n\t	kandw	%%k1,%%k5,%%k1			\n\t	kandw	%%k2,%%k6,%%k2			\n\t	kandw	%%k3,%%k7,%%k3			\n\t"/* and the 2 results together into four 8-bit bitmasks */\
				"kmovw	%%k0,%%eax					\n\t	kmovw	%%k1,%%ebx				\n\t	kmovw	%%k2,%%ecx				\n\t	kmovw	%%k3,%%edx				\n\t"\
				/* concatenate into 32-bit result and write: */\
				"shlq	$8 ,%%rbx		\n\t"\
				"shlq	$16,%%rcx		\n\t"\
				"shlq	$24,%%rdx		\n\t"\
				"addq	%%rbx,%%rax		\n\t"\
				"addq	%%rdx,%%rcx		\n\t"\
				"addq	%%rcx,%%rax		\n\t"\
				"movq	%%rax,%[__result]	\n\t"\
			/****************************************************************************************/\
			/* Do 2nd quartet of vector data, entering in reg-triplets 16-18, 20-22, 24-26, 28-30:  */\
			/****************************************************************************************/\
				"movq	%[__fq0]   ,%%rax	\n\t"\
				"movq	%[__two26f],%%rbx	\n\t"\
				"vmovaps	0x40(%%rbx),%%zmm3		\n\t	vmovaps	(%%rbx),%%zmm7			\n\t"/* two26i,f */\
			"vfmadd213pd %%zmm17,%%zmm7,%%zmm18		\n\t vfmadd213pd %%zmm21,%%zmm7,%%zmm22	\n\t vfmadd213pd %%zmm25,%%zmm7,%%zmm26	\n\t vfmadd213pd %%zmm29,%%zmm7,%%zmm30	\n\t"/* (x1 + 2^26*x2) = hi52 */\
				/* Mod-doubling code taken from SSE2_twopmodq78_modmul_q16, expects x.lo26,hi52 in zmm17,2, resp., so copy zmm16,4,8,12 -> zmm17,5,9,13 after computing x.hi52 (in zmm18,6,10,14): */\
				"vaddpd	%%zmm16,%%zmm16,%%zmm17		\n\t	vaddpd	%%zmm20,%%zmm20,%%zmm21	\n\t	vaddpd	%%zmm24 ,%%zmm24 ,%%zmm25	\n\t	vaddpd	%%zmm28,%%zmm28,%%zmm29	\n\t"/* low 26 bits */\
				"vaddpd	%%zmm18,%%zmm18,%%zmm18		\n\t	vaddpd	%%zmm22,%%zmm22,%%zmm22	\n\t	vaddpd	%%zmm26,%%zmm26,%%zmm26	\n\t	vaddpd	%%zmm30,%%zmm30,%%zmm30	\n\t"/* top 52 bits */\
				/* And subtract q.lo26,hi52 - note q.lo26 is balanced-digit normalized and divided by 2^26, so do this sub prior to final >0 normalization below: */\
				"vsubpd	0x100(%%rax),%%zmm17,%%zmm17	\n\t vsubpd	0x140(%%rax),%%zmm21,%%zmm21	\n\t vsubpd	0x180(%%rax),%%zmm25,%%zmm25	\n\t vsubpd	0x1c0(%%rax),%%zmm29,%%zmm29\n\t"/* -= qlo26 */\
				"vsubpd	0x700(%%rax),%%zmm18,%%zmm18	\n\t vsubpd	0x740(%%rax),%%zmm22,%%zmm22	\n\t vsubpd	0x780(%%rax),%%zmm26,%%zmm26\n\t vsubpd	0x7c0(%%rax),%%zmm30,%%zmm30\n\t"/* -= qhi52 */\
				/* Final-residue computation needs digits >= 0 normalization - since fq0 will get folded in with fq2 anyway, only need renormalize fq0: */\
				"vmulpd	%%zmm3,%%zmm17,%%zmm17		\n\t	vmulpd	%%zmm3,%%zmm21,%%zmm21	\n\t	vmulpd	%%zmm3,%%zmm25,%%zmm25	\n\t	vmulpd	%%zmm3,%%zmm29,%%zmm29	\n\t"/* xlo *= two26i */\
				"vrndscalepd $1,%%zmm17,%%zmm16		\n\t	vrndscalepd $1,%%zmm21,%%zmm20	\n\t	vrndscalepd $1,%%zmm25,%%zmm24	\n\t	vrndscalepd $1,%%zmm29,%%zmm28	\n\t"/* fcy */\
				"vsubpd	%%zmm16,%%zmm17,%%zmm17		\n\t	vsubpd	%%zmm20,%%zmm21,%%zmm21	\n\t	vsubpd	%%zmm24,%%zmm25,%%zmm25	\n\t	vsubpd	%%zmm28,%%zmm29,%%zmm29	\n\t"/* fx0 -= fcy */\
				"vaddpd	%%zmm16,%%zmm18,%%zmm18		\n\t	vaddpd	%%zmm20,%%zmm22,%%zmm22	\n\t	vaddpd	%%zmm24,%%zmm26,%%zmm26	\n\t	vaddpd	%%zmm28,%%zmm30,%%zmm30	\n\t"/* Add carry into xhi */\
				"vmulpd	%%zmm7,%%zmm17,%%zmm17		\n\t	vmulpd	%%zmm7,%%zmm21,%%zmm21	\n\t	vmulpd	%%zmm7,%%zmm25,%%zmm25	\n\t	vmulpd	%%zmm7 ,%%zmm29,%%zmm29	\n\t"/* fx0 *= two26f ... result must go into m0 */\
				/* (lo26 == 1 && hi52 == 0), save result in a bitmask: */
				"vmovaps	(%%rbx),%%zmm0			\n\t"/* TWO26FLOAT */\
				"vmulpd	0x40(%%rbx),%%zmm0,%%zmm3	\n\t"/* 1.0 ... no good way to generate 1.0 on-the-fly, so use 2^26 * 2^-26 */\
				"vpxorq	     %%zmm0,%%zmm0,%%zmm0	\n\t"/* 0.0 */\
				"vcmppd $0,%%zmm17,%%zmm3,%%k0		\n\t	vcmppd $0,%%zmm21,%%zmm3,%%k1	\n\t	vcmppd $0,%%zmm25,%%zmm3,%%k2	\n\t	vcmppd $0,%%zmm29,%%zmm3,%%k3	\n\t"/* low 26 bits == 1 ? In AVX512 version, mask-regs k0-3 replace dest-regs m17,5,9,13 */\
				"vcmppd $0,%%zmm18,%%zmm0,%%k4		\n\t	vcmppd $0,%%zmm22,%%zmm0,%%k5	\n\t	vcmppd $0,%%zmm26,%%zmm0,%%k6	\n\t	vcmppd $0,%%zmm30,%%zmm0,%%k7	\n\t"/* top 52 bits == 0 ? In AVX512 version, mask-regs k4-7 replace dest-regs m0,4,8,12 */\
				"kandw	%%k0,%%k4,%%k0				\n\t	kandw	%%k1,%%k5,%%k1			\n\t	kandw	%%k2,%%k6,%%k2			\n\t	kandw	%%k3,%%k7,%%k3			\n\t"/* and the 2 results together into four 8-bit bitmasks */\
				"kmovw	%%k0,%%eax					\n\t	kmovw	%%k1,%%ebx				\n\t	kmovw	%%k2,%%ecx				\n\t	kmovw	%%k3,%%edx				\n\t"\
				/* concatenate into 32-bit result and write into high 32-bit half of 64-bit result field: */\
				"shlq	$8 ,%%rbx		\n\t"\
				"shlq	$16,%%rcx		\n\t"\
				"shlq	$24,%%rdx		\n\t"\
				"addq	%%rbx,%%rax		\n\t"\
				"addq	%%rdx,%%rcx		\n\t"\
				"addq	%%rcx,%%rax		\n\t"\
				"shlq	$32,%%rax	\n\t	addq	%%rax,%[__result]	\n\t"\
				:					/* outputs: none */\
				: [__fq0] "m" (dptr)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)	\
				 ,[__result] "m" (r)	\
				: "cc","memory","cl","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10", "xmm12","xmm13","xmm14", "xmm16","xmm17","xmm18", "xmm20","xmm21","xmm22", "xmm24","xmm25","xmm26", "xmm28","xmm29","xmm30"	/* Clobbered registers */\
			);
			return r;
		}

	  #endif	// USE_AVX512 ?

	#endif	// X64_ASM ?

#endif	// #ifdef USE_IMCI512?
