/*******************************************************************************
*                                                                              *
*   (C) 1997-2012 by Ernst W. Mayer.                                           *
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
		const uint32 p, const uint32 pshift, const uint64 k,
		const uint32 start_index, const uint32 zshift,
		const int i	// thread id (handy for debug)
	)
	{
	#ifdef __CUDA_ARCH__
		 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
		uint64 hi64;
		uint96 q, qhalf, qinv, x, lo;
		uint192 prod192;
		// host variables - in this case the globals with these names - cannot be directly read in a device function:
		const uint96 ONE96  = {(uint64)1, (uint32)0};
		const double TWO26FLOAT = (double)0x04000000, TWO26FLINV = 1.0/TWO26FLOAT;
		double fq0, fq1, fq2;
		double fqinv0, fqinv1, fqinv2;
		double fx0, fx1, fx2;
		double flo0,flo1,flo2,fhi0,fhi1,fhi2;
		/* Experimental: Vars needed for post-powering-loop x * 2^(96-78) modmul: To test, try changing code below to e.g. pshift = p+96i */
		double qdiv, qmul, kmul;

		q.d0 = p+p;
		MUL_LOHI64(q.d0, k, q.d0, q.d1);	// Common HLL C-code IMUL macro shared by both device and host-code compile passes
		q.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */

		/* Convert q to floating form: */
		CVT_UINT78_3WORD_DOUBLE(q, fq0,fq1,fq2);

		RSHIFT_FAST96(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

		qinv.d0 = (q.d0 + q.d0 + q.d0) ^ (uint64)2;	qinv.d1 = (uint64)0;
		for(j = 0; j < 4; j++)
		{
			hi64 = q.d0*qinv.d0;
			qinv.d0 = qinv.d0*((uint64)2 - hi64);
		}

		MULH64(q.d0, qinv.d0, hi64);
		qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + hi64);
		qinv.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
//	if(i == 0)printf("In twopmodq78_gpu with p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,k,zshift,start_index);
		/* Convert qinv to floating form: */
		CVT_UINT78_3WORD_DOUBLE(qinv, fqinv0,fqinv1,fqinv2);
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
//	if(i == 0)printf("j = %u, fx2-0 = %20.15f, %20.15f, %20.15f\n", j, fx0,fx1,fx2);
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
				/* Normalize the result: */
				NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);
			}
		}
		/*...Double and return.	These are specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		*/
		CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);
		ADD96(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
		SUB96(x,q,x);
//	if(i == 0)printf("X_out = %u*2^64 + %llu\n", x.d1,x.d0);
		return (CMPEQ96(x, ONE96));
	#else	// ifndef __CUDA_ARCH__ :
		ASSERT(HERE, 0, "Device code being called in host mode!");
		return 0;
	#endif
	}

	__global__ void VecModpow(const uint32*pvec, const uint32*pshft, const uint32*zshft, const uint32*stidx, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
//	if(i == 0)printf("i = %u ... ",i);
		uint32 p           = pvec[i];
		uint32 pshift      = pshft[i];
		uint32 zshift      = zshft[i];
		uint32 start_index = stidx[i];
		uint64 k           = kvec[i];
		// Result returned in rvec[i]:
//	if(i == 0)printf("Calling twopmodq78_gpu with p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,k,zshift,start_index);
		rvec[i] = twopmodq78_gpu(p, pshift, k, start_index, zshift, i);
//	if(i == 0)printf("result = %u\n",rvec[i]);
	}

	__global__ void GPU_TF78(const uint32 p, const uint32 pshift, const uint32 zshift, const uint32 start_index, const uint64*kvec, uint8*rvec, int N)
	{
		int i = blockDim.x * blockIdx.x + threadIdx.x;	// Hardware thread ID:
		uint64 k           = kvec[i];
		// Result returned in rvec[i]:
	//if(i == 0)printf("Calling twopmodq78_gpu with p = %u, pshift = %u, k = %llu, zshift = %u, start_index = %u\n", p,pshift,k,zshift,start_index);
		rvec[i] = twopmodq78_gpu(p, pshift, k, start_index, zshift, i);
	}

#endif	// __CUDACC__


uint64 twopmodq78_3WORD_DOUBLE(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k)
{
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

	ASSERT(HERE, (p >> 63) == 0, "twopmodq78_q2 : p must be < 2^63!");
	q.d0 = p+p;
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(q.d0, k,&q.d0,&q.d1);
#else
	MUL_LOHI64(q.d0, k, q.d0, q.d1);
#endif
	q.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	ASSERT(HERE, (q.d1 >> 14) == 0, "twopmodq78 : (q.d1 >> 14) != 0");
	*checksum1 += q.d0;

	/* Convert q to floating form: */
	CVT_UINT78_3WORD_DOUBLE(q, fq0,fq1,fq2);

	RSHIFT_FAST96(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 96;

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

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);

if(~pshift != p+78) {
	ASSERT(HERE, ~pshift > (p+78), "twopmodq80 : Only support pshift >= true value!");
	ASSERT(HERE,(~pshift - (p+78)) < 32, "twopmodq80 : Only support pshift-diff < 32!");
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
	MUL_LOHI64( x.d0, (uint64)kmul, x.d0, hi64);
	x.d1 = hi64 + kmul*x.d1;
	lo = q;
	MUL_LOHI64(lo.d0, (uint64)qmul,lo.d0, hi64);
	lo.d1 = hi64 + qmul*lo.d1;
	SUB96(x,lo,x);
//printf("X_out[A] = %u*2^64 + %llu\n", x.d1,x.d0);
} else {
	ADD96(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
	SUB96(x,q,x);
//printf("X_out[B] = %u*2^64 + %llu\n", x.d1,x.d0);
}

	*checksum2 += x.d0;
	return (CMPEQ96(x, ONE96));
}

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
uint64 twopmodq78_3WORD_DOUBLE_q2(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1)
{
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead7, tmp0, tmp1, r;
	uint96 q0, q1, qinv0, qinv1, qhalf0, qhalf1, x0, x1, lo0, lo1;
	uint192 prod192;
	static uint64 psave = 0;
	static uint32 pshift;	/* Require this to be 32-bit here for 32-bit ASM support */
	static uint32 start_index, zshift, first_entry = TRUE;

#ifdef USE_SSE2

	static struct complex *sc_arr = 0x0;
	static         double *sc_ptr;
	static double *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2;
	static double *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2;
	static double *half, *two26f, *two26i, *two13i, *sse2_rnd;
	static uint64 ihalf = 0x3FDfffffffffffffull;	/* Bitfield storing 0.5*(1-epsilon) in IEEE64 format */

#else

	double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2, fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2;
	double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2;

#endif

#ifdef USE_SSE2

	if(first_entry)
	{
		sc_arr = ALLOC_COMPLEX(sc_arr, 44);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = (double *)ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		/* Remember, rhese are POINTERS-TO-DOUBLES, so need an increment of 2 to span an SSE register: */
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
		*sse2_rnd++ = 3.0*0x4000000*0x2000000;
		*sse2_rnd-- = 3.0*0x4000000*0x2000000;
		/* We init "half" = 0.5-epsilon here, because emulating FLOOR(x) via DNINT(x-half) requires
		us to always round up if x is a whole number, and our DNINT emulation can round either way if fractional part = 0.5:
		*/
		*half++ = *(double*)&ihalf;	*half-- = *(double*)&ihalf;
	}	/* first_entry */
#endif

	ASSERT(HERE, (p >> 63) == 0, "twopmodq78_q2 : p must be < 2^63!");
	q0.d0 = q1.d0 = p+p;
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(q0.d0, k0,&q0.d0,&q0.d1);
	MUL_LOHI64(q1.d0, k1,&q1.d0,&q1.d1);
#else
	MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
	MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
#endif
	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q2 : (q0.d1 >> 14) != 0");
	ASSERT(HERE, (q1.d1 >> 14) == 0, "twopmodq78_q2 : (q1.d1 >> 14) != 0");
	*checksum1 += q0.d0 + q1.d0;

	/* Convert q to floating form: */
#ifdef USE_SSE2

	CVT_UINT78_3WORD_DOUBLE(q0 ,*fq0,*fq1,*fq2);
	CVT_UINT78_3WORD_DOUBLE(q1 ,*gq0,*gq1,*gq2);

	#if defined(COMPILER_TYPE_MSVC)

		/* Scale by TWO26FLINV: */
		__asm	mov	eax, fq0
		__asm	mov	ebx, two26f
		__asm	movaps	xmm5,[ebx     ]	/* two26f */
		__asm	movaps	xmm6,[ebx+0x10]	/* two26i */
		__asm	movaps	xmm0,[eax     ]	/* fq0 */
		__asm	movaps	xmm1,[eax+0x10]	/* fq1 */
		__asm	movaps	xmm2,[eax+0x20]	/* fq2 */
		__asm	movaps	xmm3,xmm2		/* cpy fq2 */
		__asm	mulpd	xmm3,xmm5
		__asm	addpd	xmm3,xmm1		/* Hi 52 bits of q */
		__asm	mulpd	xmm0,xmm6
		__asm	mulpd	xmm1,xmm6
		__asm	mulpd	xmm2,xmm6
		__asm	movaps	[eax     ],xmm0	/* fq0/2^26 */
		__asm	movaps	[eax+0x10],xmm1	/* fq1/2^26 */
		__asm	movaps	[eax+0x20],xmm2	/* fq2/2^26 */
		__asm	movaps	[eax+0x30],xmm3	/* fqhi52 */

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			__asm__ volatile (\
				"movl	%[__fq0],%%eax		\n\t"\
				"movl	%[__two26f],%%esi	\n\t"\
				"movaps	    (%%esi),%%xmm5	\n\t"\
				"movaps	0x10(%%esi),%%xmm6	\n\t"\
				"movaps	    (%%eax),%%xmm0	\n\t"\
				"movaps	0x10(%%eax),%%xmm1	\n\t"\
				"movaps	0x20(%%eax),%%xmm2	\n\t"\
				"movaps	%%xmm2,%%xmm3		\n\t"\
				"mulpd	%%xmm5,%%xmm3		\n\t"\
				"addpd	%%xmm1,%%xmm3		\n\t"\
				"mulpd	%%xmm6,%%xmm0		\n\t"\
				"mulpd	%%xmm6,%%xmm1		\n\t"\
				"mulpd	%%xmm6,%%xmm2		\n\t"\
				"movaps	%%xmm0,    (%%eax)	\n\t"\
				"movaps	%%xmm1,0x10(%%eax)	\n\t"\
				"movaps	%%xmm2,0x20(%%eax)	\n\t"\
				"movaps	%%xmm3,0x30(%%eax)	\n\t"\
				:					/* outputs: none */\
				: [__fq0] "m" (fq0)	/* All inputs from memory addresses here */\
				,[__two26f] "m" (two26f)\
				: "cc","memory","eax","esi","xmm0","xmm1","xmm2","xmm3","xmm5","xmm6"	/* Clobbered registers */\
			);

		#else

			__asm__ volatile (\
				"movq	%[__fq0],%%rax 	\n\t"\
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

	#endif

#else
	CVT_UINT78_3WORD_DOUBLE(q0, fq0,fq1,fq2);
	CVT_UINT78_3WORD_DOUBLE(q1, gq0,gq1,gq2);
#endif

	RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST96(q1, 1, qhalf1);

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		ASSERT(HERE, ((p + 78) >> 32) == 0, "twopmodq78_q2 : (p+78) exceeds 32 bits!");
		pshift = (uint32)p + 78;

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
		j = leadz32(pshift);
		/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
		lead7 = ((pshift<<j) >> 25);
		if(lead7 > 77)
		{
			lead7 >>= 1;
			start_index =  32-j-6;	/* Use only the leftmost 6 bits */
		}
		else
			start_index =  32-j-7;

		zshift = 77 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

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

	#if defined(COMPILER_TYPE_MSVC)

		/* Scale by TWO26FLINV: */
		__asm	mov	eax, fqinv0
		__asm	mov	ebx, two26i
		__asm	movaps	xmm6,[ebx]	/* two26i */
		__asm	movaps	xmm0,[eax     ]	/* fqinv0 */
		__asm	movaps	xmm1,[eax+0x10]	/* fqinv1 */
		__asm	movaps	xmm2,[eax+0x20]	/* fqinv2 */
		__asm	mulpd	xmm0,xmm6
		__asm	mulpd	xmm1,xmm6
		__asm	mulpd	xmm2,xmm6
		__asm	movaps	[eax     ],xmm0	/* fqinv0 */
		__asm	movaps	[eax+0x10],xmm1	/* fqinv1 */
		__asm	movaps	[eax+0x20],xmm2	/* fqinv2 */

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

		__asm__ volatile (\
			"movl	%[__fqinv0],%%eax 	\n\t"\
			"movl	%[__two26i],%%esi 	\n\t"\
			"movaps	    (%%esi),%%xmm6	\n\t"\
			"movaps	    (%%eax),%%xmm0	\n\t"\
			"movaps	0x10(%%eax),%%xmm1	\n\t"\
			"movaps	0x20(%%eax),%%xmm2	\n\t"\
			"mulpd	%%xmm6,%%xmm0		\n\t"\
			"mulpd	%%xmm6,%%xmm1		\n\t"\
			"mulpd	%%xmm6,%%xmm2		\n\t"\
			"movaps	%%xmm0,    (%%eax)	\n\t"\
			"movaps	%%xmm1,0x10(%%eax)	\n\t"\
			"movaps	%%xmm2,0x20(%%eax)	\n\t"\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			: "cc","memory","eax","esi","xmm0","xmm1","xmm2","xmm6"	/* Clobbered registers */\
		);

		#else

		__asm__ volatile (\
			"movq	%[__fqinv0],%%rax 	\n\t"\
			"movq	%[__two26i],%%rsi 	\n\t"\
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

	#endif

#else
	CVT_UINT78_3WORD_DOUBLE(qinv0, fqinv0,fqinv1,fqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv1, gqinv0,gqinv1,gqinv2);
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

	#if defined(COMPILER_TYPE_MSVC)

		__asm	mov	eax, fx0
		__asm	movaps	xmm0,[eax     ]	/* fx0 */
		__asm	movaps	xmm2,[eax+0x10]	/* fx1 */
		__asm	movaps	xmm4,[eax+0x20]	/* fx2 */

	#else	/* GCC-style inline ASM: */

		/* GCC builds require us to save the xmm0/2/4 registers on each pass through the loop,
		so skip these pre-loop-entry snippets until we get the full loop-logic asm'ed:
		*/
		/* 11/10/2010: Update - Since I found it necessary to build with -O1 for the SSE2 code
		to run properly anyway, tried activating these here and commenting them out inside the
		loop body. That works - but strangely, it`s slower in my Mingw64-for-windows-7 build.
		(Slightly faster on amd64/linux).
		*/
		#if OS_BITS == 32
			__asm__ volatile (\
				"movl	%[__fx0],%%eax 	\n\t"\
				"movaps	    (%%eax),%%xmm0	\n\t"\
				"movaps	0x10(%%eax),%%xmm2	\n\t"\
				"movaps	0x20(%%eax),%%xmm4	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "cc","memory","eax","xmm0","xmm2","xmm4"	/* Clobbered registers */\
			);
		#else	// ^4-bit mode has inputs in xmm0,1,2:
			__asm__ volatile (\
				"movq	%[__fx0],%%rax 	\n\t"\
				"movaps	    (%%rax),%%xmm0	\n\t"\
				"movaps	0x10(%%rax),%%xmm1	\n\t"\
				"movaps	0x20(%%rax),%%xmm2	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm1","xmm2"	/* Clobbered registers */\
			);
		#endif

	#endif

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

#elif(defined(USE_SSE2))
/* Inner loop body needs 42 movaps, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */

  #if 0	/* Set = 1 to do timings of all-but-inner-loop */

  #else	/* Inner loop body needs 42 movaps, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */

	#if defined(COMPILER_TYPE_MSVC)

	  #define ASM_LOOP 0

	  #if ASM_LOOP	/* Pure-ASM loop control */
		__asm	mov ecx, start_index
		__asm	sub ecx, 2
		__asm	test ecx, ecx	/* int tmp = esi & esi = n*/
		__asm	jl LoopEnd		/* Skip if n < 0 */
		LoopBeg:
	  #else
		for(j = start_index-2; j >= 0; j--)
		{
	  #endif
		/*
		SQR_LOHI78_3WORD_DOUBLE_q2(
			  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2
			, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2
		);
		*/
			__asm	mov	eax, fx0
			__asm	mov	ebx, two26i
			__asm	movaps	xmm6,[ebx-0x20]	/* two13i */
			/* fx0,1,2 assumed in xmm0,2,4 on loop entry */
			__asm	mulpd	xmm0,xmm6	/* scale [fx0 . 2^-13] */
			__asm	mulpd	xmm2,xmm6	/* scale [fx1 . 2^-13] */
			__asm	mulpd	xmm4,xmm6	/* scale [fx2 . 2^-13] */
			__asm	movaps	xmm7,[ebx+0x10]	/* xmm7 = RND */
			__asm	movaps	xmm1,xmm0	/* cpy of fx0 */
			__asm	movaps	xmm3,xmm2	/* cpy of fx1 */
			__asm	addpd	xmm1,xmm1	/* 2.fx0 */
			__asm	addpd	xmm3,xmm3	/* 2.fx1 */
			__asm	movaps	xmm5,xmm1	/* cpy of 2.fx0 */
			__asm	mulpd	xmm0,xmm0	/* [  fx0 *= fx0] / 2^26 */
			__asm	mulpd	xmm1,xmm2	/* [2.fx0 *= fx1] / 2^26 */
			__asm	mulpd	xmm2,xmm2	/* [  fx1 *= fx1] / 2^26 */
			__asm	mulpd	xmm5,xmm4	/* [2.fx0 *= fx2] / 2^26 */
			__asm	mulpd	xmm3,xmm4	/* [2.fx1 *= fx2] / 2^26 */
			__asm	mulpd	xmm4,xmm4	/* [  fx2 *= fx2] / 2^26 */
			/* Digit 0:
			__fprod0  =  __fx0*__fx0;		in [0, 2^50]
			__fcy     = DNINT(__fprod0*TWO26FLINV);
			__fprod0 -= __fcy*TWO26FLOAT;	in [-2^24, +2^24]
			*/
			__asm	movaps	xmm6,xmm0		/* Init: fcy = cpy of fx0*fx0 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7		/* fcy */
			__asm	subpd	xmm0,xmm6		/* fprod0 -= fcy */
			__asm	mulpd	xmm0,[ebx-0x10]	/* fprod0 *= two26f */
			__asm	mulpd	xmm6,[ebx     ]	/* fcy    *= two26i */

			/* Digit 1:
			__fprod1  = __f2x0*__fx1 + __fcy;	in [-2^51-2^24, +2^51+2^24]
			__fcy     = DNINT(__fprod1*TWO26FLINV);
			__fprod1 -= __fcy*TWO26FLOAT;	in [-2^25, +2^25]
			*/
			__asm	addpd	xmm1,xmm6		/* fprod1 = 2.fx0*fx1 + fcy */
			__asm	movaps	xmm6,xmm1		/* fcy = cpy of fprod1 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7		/* fcy */
			__asm	subpd	xmm1,xmm6		/* fprod1 -= fcy */
			__asm	mulpd	xmm1,[ebx-0x10]	/* fprod1 *= two26f */
			__asm	mulpd	xmm6,[ebx     ]	/* fcy    *= two26i */

			/* Digit 2:
			__fprod2  = __f2x0*__fx2 + __fx1*__fx1 + __fcy;	in [-2^52-2^25, +2^52+2^25]
			__fcy     = DNINT(__fprod2*TWO26FLINV);
			__fprod2 -= __fcy*TWO26FLOAT;
			__itmp    = (__fprod2 < 0);
			__fcy    -= (double)__itmp;
			__fprod2 += (double)(__itmp << 26);
			*/
			__asm	addpd	xmm2,xmm5		/* fx1*fx1 + 2.fx0*fx2,		xmm5 FREE */
			__asm	movaps	xmm5,[ebx-0x10]	/* 2^26 */
			__asm	addpd	xmm2,xmm6		/* fprod2 = 2.fx0*fx2 + fx1*fx1 + fcy */
			__asm	movaps	xmm6,xmm2		/* fcy = cpy of fprod2 */
			__asm	subpd	xmm6,[ebx+0x20]	/* fprod2 - 0.5 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7		/* fcy = floor(fprod2) */
			__asm	subpd	xmm2,xmm6		/* fprod2 -= fcy */
			__asm	mulpd	xmm2,xmm5		/* fprod2 *= two26f */
			__asm	mulpd	xmm6,[ebx]		/* fcy    *= two26i */

			/* Digit 3:
			__fprod3  = __f2x1*__fx2 + __fcy;
			__fcy     = DNINT(__fprod3*TWO26FLINV);
			__fprod3 -= __fcy*TWO26FLOAT;
			*/
			__asm	addpd	xmm3,xmm6	/* fprod3 = 2.fx1*fx2 + fcy */
			__asm	movaps	xmm6,xmm3	/* fcy = cpy of fprod3 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm3,xmm6	/* fprod3 -= fcy */
			__asm	mulpd	xmm3,xmm5	/* fprod3 *= two26f */
			__asm	mulpd	xmm6,[ebx]	/* fcy    *= two26i */

			/* Digit 4:
			__fprod4  =  __fx2*__fx2 + __fcy;
			__fcy     = DNINT(__fprod4*TWO26FLINV);
			__fprod4 -= __fcy*TWO26FLOAT;
			*/
			__asm	addpd	xmm4,xmm6	/* fprod4 = fx2*fx2 + fcy */
			__asm	movaps	xmm6,xmm4	/* fcy = cpy of fprod4 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm4,xmm6	/* fprod4 -= fcy */
			__asm	mulpd	xmm4,xmm5	/* fprod4 *= two26f */

			/* Digit 5 = the carry.
			flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6
			*/
			__asm	mulpd	xmm6,xmm5			/* fhi2 * two26f */
			__asm	addpd	xmm4,xmm6			/* fhi, top 52 bits */
			__asm	movaps	[eax+0x30],xmm3	/* Store fhi0,1,2 to free up 3 registers */
			__asm	movaps	[eax+0x40],xmm4

		/* MULL78_3WORD_DOUBLE_q2(
			  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
			, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2);
		*/
			/* Precompute all the needed partial products: x0,1,2 = flo0,1,2; y0,1,2 = fqinv0,1,2:
			__prod0  = __x0*__y0;
			__prod1  = __x0*__y1 + __x1*__y0;
			__prod2  = __x0*__y2 + __x1*__y1 + __x2*__y0;
			*/
			__asm	mov	edx, fqinv0

			/* Digit 0:
			__fcy    = DNINT(__prod0*TWO26FLINV);
			*__lo0   = __prod0 - __fcy*TWO26FLOAT;
			*/
			__asm	movaps	xmm3,xmm0		/* xmm3 = cpy of x0 */
			__asm	movaps	xmm4,xmm0		/* xmm4 = cpy of x0 */
			__asm	mulpd	xmm0,[edx     ]	/* fprod0 = x0*y0 */
			__asm	movaps	xmm6,xmm0		/* fcy = cpy of fprod0 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7		/* fcy */
			__asm	subpd	xmm0,xmm6		/* fprod0 -= fcy */
			__asm	mulpd	xmm0,xmm5		/* fprod0 *= two26f */
			__asm	mulpd	xmm6,[ebx]		/* fcy    *= two26i */

			/* Digit 1:
			__prod1 += __fcy;
			__fcy    = DNINT(__prod1*TWO26FLINV);
			*__lo1   = __prod1 - __fcy*TWO26FLOAT;
			*/
			__asm	mulpd	xmm3,[edx+0x10]	/* x0*y1 */
			__asm	addpd	xmm3,xmm6		/* x0*y1 + fcy; xmm6 FREE */
			__asm	movaps	xmm6,xmm1		/* xmm6 = cpy of x1 */
			__asm	mulpd	xmm1,[edx     ]	/* x1*y0 */
			__asm	addpd	xmm1,xmm3		/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */
			__asm	movaps	xmm3,xmm1		/* fcy = cpy of fprod1 */
			__asm	addpd	xmm3,xmm7
			__asm	subpd	xmm3,xmm7		/* fcy */
			__asm	subpd	xmm1,xmm3		/* fprod1 -= fcy */
			__asm	mulpd	xmm1,xmm5		/* fprod1 *= two26f */
			__asm	mulpd	xmm3,[ebx]		/* fcy    *= two26i */

			/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced:
			__prod2 += __fcy;
			__fcy    = DNINT(__prod2*TWO26FLINV);
			*__lo2   = __prod2 - __fcy*TWO26FLOAT;
			*/
			__asm	mulpd	xmm2,[edx     ]	/* x2*y0 */
			__asm	mulpd	xmm6,[edx+0x10]	/* x1*y1 */
			__asm	mulpd	xmm4,[edx+0x20]	/* x0*y2 */
			__asm	addpd	xmm2,xmm6		/* x1*y1 + x2*y0; xmm6 FREE */
			__asm	addpd	xmm3,xmm4		/* x0*y2 + fcy; xmm4 FREE */
			__asm	addpd	xmm2,xmm3		/* fprod2; xmm3 FREE */
			__asm	movaps	xmm3,xmm2		/* fcy = cpy of fprod2 */
			__asm	subpd	xmm3,[ebx+0x20]	/* fprod2 - 0.5 */
			__asm	addpd	xmm3,xmm7
			__asm	subpd	xmm3,xmm7		/* fcy */
			__asm	subpd	xmm2,xmm3		/* fprod2 -= fcy */
			__asm	mulpd	xmm2,xmm5		/* fprod2 *= two26f */

		/* MULH78(q,lo,lo) --> lo = (q*lo)/2^78
		MULH78_3WORD_DOUBLE_q2(
			  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2
			, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2 );
		*/
			/* Precompute all the needed partial products: x0,1,2 = (fq0,1,2)/2^26 ; y0,1,2 = flo0,1,2: */
			__asm	mov	edx, fq0

			/* Digit 0:
			__ftmp =  __fx0*__fy0;
			__fcy  = DNINT(__ftmp*TWO26FLINV);
			*/
			__asm	movaps	xmm3,xmm0	/* xmm3 = cpy of y0 */
			__asm	movaps	xmm4,xmm0	/* xmm4 = cpy of y0 */
			__asm	mulpd	xmm0,[edx]	/* fprod0 = y0*x0 */
//	__asm	addpd	xmm0,xmm7
//	__asm	subpd	xmm0,xmm7	/* fcy */
			__asm	mulpd	xmm0,[ebx]	/* fcy    *= two26i */

			/* Digit 1:
			__ftmp = __fx0*__fy1 + __fx1*__fy0 + __fcy;
			__fcy  = DNINT(__ftmp*TWO26FLINV);
			*/
			__asm	movaps	xmm6,xmm1		/* xmm6 = cpy of y1 */
			__asm	mulpd	xmm3,[edx+0x10]	/* y0*x1 */
			__asm	addpd	xmm3,xmm0		/* y0*x1 + fcy; xmm0 FREE */
			__asm	mulpd	xmm1,[edx     ]	/* y1*x0 */
			__asm	addpd	xmm1,xmm3		/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */
//	__asm	addpd	xmm1,xmm7
//	__asm	subpd	xmm1,xmm7	/* fcy */
			__asm	mulpd	xmm1,[ebx]	/* fcy    *= two26i */

			/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced:
			__ftmp = __fx0*__fy2 + __fx1*__fy1 + __fx2*__fy0 + __fcy;
			__fcy  = DNINT(__ftmp*TWO26FLINV);
			__ftmp-= __fcy*TWO26FLOAT;
			__itmp = (__ftmp < 0);
			__fcy -= (double)__itmp;
			*/
			__asm	movaps	xmm3,xmm2		/* xmm3 = cpy of y2 */
			__asm	movaps	xmm0,xmm6		/* xmm0 = cpy of y1 */
			__asm	mulpd	xmm2,[edx     ]	/* y2*x0 */
			__asm	mulpd	xmm6,[edx+0x10]	/* y1*x1 */
			__asm	mulpd	xmm4,[edx+0x20]	/* y0*x2 */
			__asm	addpd	xmm2,xmm6		/* y1*x1 + y2*x0; xmm6 FREE */
			__asm	addpd	xmm1,xmm4		/* y0*x2 + fcy; xmm4 FREE */
			__asm	addpd	xmm2,xmm1		/* fprod2; xmm1 FREE */
			__asm	subpd	xmm2,[ebx+0x20]	/* fprod2 - 0.5 */
			__asm	addpd	xmm2,xmm7
			__asm	subpd	xmm2,xmm7	/* fcy */
			__asm	mulpd	xmm2,[ebx]	/* fcy    *= two26i */

		/* xmm1,3,4,6 FREE */
			/* Precompute all the needed partial products: */

			__asm	movaps	xmm1,xmm3		/* xmm1 = cpy of y2 */
			__asm	mulpd	xmm0,[edx+0x20]	/* y1*x2 */
			__asm	mulpd	xmm1,[edx+0x10]	/* y2*x1 */
			__asm	mulpd	xmm3,[edx+0x20]	/* y2*x2 */

		/* xmm4,6 FREE */
			/* Digit 3:
			__fprod3 = __fx1*__fy2 + __fx2*__fy1 + __fcy;
			__fcy  = DNINT(__fprod3*TWO26FLINV);
			__fhi0 = __fprod3 - __fcy*TWO26FLOAT;
			*/
			__asm	addpd	xmm1,xmm2	/* fy2*fx1 + fcy */
			__asm	addpd	xmm0,xmm1	/* fprod3 = fy1*fx2 + fy2*fx1 + fcy */
			__asm	movaps	xmm6,xmm0	/* fcy = cpy of fprod3 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm0,xmm6	/* fprod3 -= fcy */
			__asm	mulpd	xmm0,xmm5	/* fprod3 *= two26f */
			__asm	mulpd	xmm6,[ebx]	/* fcy    *= two26i */

		/* xmm1,2,4 FREE */
			/* Digit 4:
			__fprod4 = __fx2*__fy2;
			__fprod4 += __fcy;
			__fcy  = DNINT(__fprod4*TWO26FLINV);
			__fhi1 = __fprod4 - __fcy*TWO26FLOAT;
			*/
			__asm	addpd	xmm3,xmm6	/* fcy = fy2*fx2 + fcy */
			__asm	movaps	xmm1,xmm3	/* fprod4 = cpy of fcy */
			__asm	addpd	xmm3,xmm7
			__asm	subpd	xmm3,xmm7	/* fcy */
			__asm	subpd	xmm1,xmm3	/* fprod4 -= fcy */
			__asm	mulpd	xmm1,xmm5	/* fprod4 *= two26f */

		/* xmm2,4,6 FREE, mulh(q,lo) digits0,1,2 in xmm0,1,3 */

			/* If h < l, then calculate h-l+q; otherwise calculate h-l. Use leading 52 bits to approximate the full 78-bit compare.
			   Result is in [0, q). */
			__asm	movaps	xmm2,[eax+0x40]	/* fhi, top 52 bits */
			__asm	movaps	xmm6,xmm2			/* cpy of fhi */
			__asm	movaps	xmm4,[edx+0x30]	/* fq, top 52 bits */
			__asm	mulpd	xmm3,xmm5			/* flo2 * two26f */
			__asm	mulpd	xmm5,[edx]			/* fq, low 26 bits */
			__asm	addpd	xmm3,xmm1			/* flo, top 52 bits */
			__asm	movaps	xmm1,[eax+0x30]	/* fhi, low 26 bits */
			__asm	cmppd	xmm6,xmm3,0x1		/* bitmask = (fhi < flo)? */
			__asm	subpd	xmm1,xmm0			/* (fhi - flo), low 26 bits */
			__asm	subpd	xmm2,xmm3			/* (fhi - flo), top 52 bits */
			__asm	movaps	xmm0,xmm4			/* cpy of qhi52 */

			__asm	andpd	xmm4,xmm6			/* qhi52 & bitmask */
			__asm	addpd	xmm2,xmm4			/* xhi = (h-l)hi + (qhi52 & bitmask) */
			__asm	andpd	xmm6,xmm5			/* bitmask & qlo26 */
			__asm	addpd	xmm1,xmm6			/* xlo = (h-l)lo + (qlo26 & bitmask) */
												/* qlo26 in xmm5, qhi52 in xmm0 */
			/* if((pshift >> j) & (uint64)1) { */
				__asm	mov	eax, pshift
	  #if !ASM_LOOP
			__asm	mov ecx, j
	  #endif
				__asm	shr	eax, cl		/* j stored in ecx; only need lowest byte for shift count */
				__asm	and	eax, 0x1	/* (p odd)? = ((pshift >> j) & (uint64)1) */
		#if 1	// Branchless version of the conditional doubling
				__asm	xor	eax, 0x1	/* If result = 0, want a 0x0000.... double bitmask, otherwise want a 0xffff... mask */
				__asm	movd	xmm6,eax	/* xmm[0:31] <- eax */
				__asm	pshufd	xmm6,xmm6,0	/* Broadcast low 32 bits of xmm6 to all 4 slots of xmm0 (But only care about bottom 2 slots of result) */
				__asm cvtdq2pd	xmm6,xmm6	/* Convert to double */
				__asm	xorpd	xmm3,xmm3	/* All 0s */
				__asm	cmppd	xmm6,xmm3,0	/* bitmask = ((pshift >> j) odd)? */

				/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */
				__asm	movaps	xmm4,xmm2		/* cpy of xhi */
				__asm	movaps	xmm3,xmm1		/* cpy of xlo */
				__asm	andpd	xmm4,xmm6		/* xhi_ & bitmask */
				__asm	andpd	xmm3,xmm6		/* xlo_ & bitmask */
				__asm	addpd	xmm2,xmm4		/* top 52 bits */
				__asm	addpd	xmm1,xmm3		/* low 26 bits */

				/* If x > q, subtract q: */
				__asm	movaps	xmm6,xmm0		/* cpy of qhi */
				__asm	cmppd	xmm6,xmm2,0x2	/* bitmask = (qhi <= xhi)? */
				__asm	andpd	xmm0,xmm6		/* qhi52 & bitmask */
				__asm	andpd	xmm5,xmm6		/* qlo26 & bitmask */
				__asm	subpd	xmm2,xmm0		/* x % q, top 52 bits */
				__asm	subpd	xmm1,xmm5		/* x % q, low 26 bits */
		#else
			__asm	je	SHORT twopmodq78_3wdq2
				/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */
				__asm	movaps	xmm6,xmm0		/* cpy of qhi */
				__asm	addpd	xmm2,xmm2		/* top 52 bits */
				__asm	addpd	xmm1,xmm1		/* low 26 bits */

				/* If x > q, subtract q: */
				__asm	cmppd	xmm6,xmm2,0x2	/* bitmask = (qhi <= xhi)? */
				__asm	andpd	xmm0,xmm6		/* qhi52 & bitmask */
				__asm	andpd	xmm5,xmm6		/* qlo26 & bitmask */
				__asm	subpd	xmm2,xmm0		/* x % q, top 52 bits */
				__asm	subpd	xmm1,xmm5		/* x % q, low 26 bits */
			twopmodq78_3wdq2:
		#endif
			/* } */

			/* Normalize the result:
			__fcy  = DNINT(__xlo*TWO26FLINV);
			__xlo -=  __fcy*TWO26FLOAT;
			*/
			__asm	movaps	xmm4,[ebx-0x10]	/* two26f */
			__asm	movaps	xmm3,[ebx     ]	/* two26i */
			__asm	mulpd	xmm1,xmm3		/* xlo *= two26i */
			__asm	mulpd	xmm2,xmm3		/* xhi *= two26i */
			__asm	movaps	xmm0,xmm1	/* Init: fcy = cpy of ~xlo */
			__asm	addpd	xmm1,xmm7
			__asm	subpd	xmm1,xmm7	/* fcy */
			__asm	subpd	xmm0,xmm1	/* fx0 -= fcy */
			__asm	mulpd	xmm0,xmm4	/* fx0 *= two26f */
			__asm	mulpd	xmm1,xmm3	/* fcy *= two26i */
			__asm	addpd	xmm2,xmm1	/* Add carry into xhi */
			__asm	movaps	xmm1,xmm2	/* Init: fcy = cpy of ~xhi */
			__asm	addpd	xmm2,xmm7
			__asm	subpd	xmm2,xmm7	/* fx2 = fcy (no divide by 2^26) */
			__asm	subpd	xmm1,xmm2	/* fx1 -= fcy */
			__asm	mulpd	xmm1,xmm4	/* fx1 *= two26f */

			/* Move high 2 words of result into input registers expected by start of loop body: */
			__asm	movaps	xmm4,xmm2	/* fx2 */
			__asm	movaps	xmm2,xmm1	/* fx1 */

	  #if ASM_LOOP
		__asm	sub ecx, 1	/* j-- */
		__asm	cmp ecx, 0	/* j > 0 ?	*/
		__asm	jge LoopBeg	/* if (j >= 0), Loop */
		LoopEnd:
	  #else	/* Pure-ASM loop control */
		}
	  #endif

	#else	/* GCC-style inline ASM: */

		for(j = start_index-2; j >= 0; j--)
		{
		#if OS_BITS == 32

		__asm__ volatile (\
		"/* SQR_LOHI78_3WORD_DOUBLE_q2(fx, flo,fhi): */\n\t"\
			"movl	%[__fx0],%%eax		\n\t"\
			"movl	%[__two26i],%%esi	\n\t"\
			"movaps	-0x20(%%esi),%%xmm6	\n\t"/* two13i */\
			"/* fx0,1,2 assumed in xmm0,2,4 on loop entry */\n\t"\
			"mulpd	%%xmm6,%%xmm0		\n\t"/* scale [fx0 . 2^-13] */\
			"mulpd	%%xmm6,%%xmm2		\n\t"/* scale [fx1 . 2^-13] */\
			"mulpd	%%xmm6,%%xmm4		\n\t"/* scale [fx2 . 2^-13] */\
			"movaps	0x10(%%esi),%%xmm7	\n\t"/* xmm7 = RND */\
			"movaps	%%xmm0,%%xmm1		\n\t"/* cpy of fx0 */\
			"movaps	%%xmm2,%%xmm3		\n\t"/* cpy of fx1 */\
			"addpd	%%xmm1,%%xmm1		\n\t"/* 2.fx0 */\
			"addpd	%%xmm3,%%xmm3		\n\t"/* 2.fx1 */\
			"movaps	%%xmm1,%%xmm5		\n\t"/* cpy of 2.fx0 */\
			"mulpd	%%xmm0,%%xmm0		\n\t"/* [  fx0 *= fx0] / 2^26 */\
			"mulpd	%%xmm2,%%xmm1		\n\t"/* [2.fx0 *= fx1] / 2^26 */\
			"mulpd	%%xmm2,%%xmm2		\n\t"/* [  fx1 *= fx1] / 2^26 */\
			"mulpd	%%xmm4,%%xmm5		\n\t"/* [2.fx0 *= fx2] / 2^26 */\
			"mulpd	%%xmm4,%%xmm3		\n\t"/* [2.fx1 *= fx2] / 2^26 */\
			"mulpd	%%xmm4,%%xmm4		\n\t"/* [  fx2 *= fx2] / 2^26 */\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm6		\n\t"/* Init: fcy = cpy of fx0*fx0 */\
			"addpd	%%xmm7,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm6		\n\t"/* fcy */\
			"subpd	%%xmm6,%%xmm0		\n\t"/* fprod0 -= fcy */\
			"mulpd	-0x10(%%esi),%%xmm0	\n\t"/* fprod0 *= two26f */\
			"mulpd	     (%%esi),%%xmm6	\n\t"/* fcy    *= two26i */\
			"/* Digit 1: */\n\t"\
			"addpd	%%xmm6,%%xmm1		\n\t"/* fprod1 = 2.fx0*fx1 + fcy */\
			"movaps	%%xmm1,%%xmm6		\n\t"/* fcy = cpy of fprod1 */\
			"addpd	%%xmm7,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm6		\n\t"/* fcy */\
			"subpd	%%xmm6,%%xmm1		\n\t"/* fprod1 -= fcy */\
			"mulpd	-0x10(%%esi),%%xmm1	\n\t"/* fprod1 *= two26f */\
			"mulpd	     (%%esi),%%xmm6	\n\t"/* fcy    *= two26i */\
			"/* Digit 2: */\n\t"\
			"addpd	%%xmm5,%%xmm2		\n\t"/* fx1*fx1 + 2.fx0*fx2,		xmm5 FREE */\
			"movaps	-0x10(%%esi),%%xmm5	\n\t"/* 2^26 */\
			"addpd	%%xmm6,%%xmm2		\n\t"/* fprod2 = 2.fx0*fx2 + fx1*fx1 + fcy */\
			"movaps	%%xmm2,%%xmm6		\n\t"/* fcy = cpy of fprod2 */\
			"subpd	0x20(%%esi),%%xmm6	\n\t"/* fprod2 - 0.5 */\
			"addpd	%%xmm7,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm6		\n\t"/* fcy = floor(fprod2) */\
			"subpd	%%xmm6,%%xmm2		\n\t"/* fprod2 -= fcy */\
			"mulpd	%%xmm5,%%xmm2		\n\t"/* fprod2 *= two26f */\
			"mulpd	    (%%esi),%%xmm6	\n\t"/* fcy    *= two26i */\
			"/* Digit 3: */\n\t"\
			"addpd	%%xmm6,%%xmm3		\n\t"/* fprod3 = 2.fx1*fx2 + fcy */\
			"movaps	%%xmm3,%%xmm6		\n\t"/* fcy = cpy of fprod3 */\
			"addpd	%%xmm7,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm6		\n\t"/* fcy */\
			"subpd	%%xmm6,%%xmm3		\n\t"/* fprod3 -= fcy */\
			"mulpd	%%xmm5,%%xmm3		\n\t"/* fprod3 *= two26f */\
			"mulpd	    (%%esi),%%xmm6	\n\t"/* fcy    *= two26i */\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"/* fprod4 = fx2*fx2 + fcy */\
			"movaps	%%xmm4,%%xmm6		\n\t"/* fcy = cpy of fprod4 */\
			"addpd	%%xmm7,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm6		\n\t"/* fcy */\
			"subpd	%%xmm6,%%xmm4		\n\t"/* fprod4 -= fcy */\
			"mulpd	%%xmm5,%%xmm4		\n\t"/* fprod4 *= two26f */\
			"/* Digit 5 = the carry. flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6 */\n\t"\
			"mulpd	%%xmm5,%%xmm6		\n\t"/* fhi2 * two26f */\
			"addpd	%%xmm6,%%xmm4		\n\t"/* fhi, top 52 bits; xmm6 FREE */\
			"movaps	%%xmm3,0x30(%%eax)	\n\t"/* Store fhi0,1 = fprod3,4 to free up 2 more registers */\
			"movaps	%%xmm4,0x40(%%eax)	\n\t"/* Recall: fx0-2 in eax+0x[012]0, flo0-2 in eax+0x[345]0, fhi0-2 in eax+0x[678]0. */\
		"/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */\n\t"\
			"movl	%[__fqinv0],%%edx	\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm3		\n\t"/* xmm3 = cpy of x0 */\
			"movaps	%%xmm0,%%xmm4		\n\t"/* xmm4 = cpy of x0 */\
			"mulpd	    (%%edx),%%xmm0	\n\t"/* fprod0 = x0*y0 */\
			"movaps	%%xmm0,%%xmm6		\n\t"/* fcy = cpy of fprod0 */\
			"addpd	%%xmm7,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm6		\n\t"/* fcy */\
			"subpd	%%xmm6,%%xmm0		\n\t"/* fprod0 -= fcy */\
			"mulpd	%%xmm5,%%xmm0		\n\t"/* fprod0 *= two26f */\
			"mulpd	    (%%esi),%%xmm6	\n\t"/* fcy    *= two26i */\
			"/* Digit 1: */\n\t"\
			"mulpd	0x10(%%edx),%%xmm3	\n\t"/* x0*y1 */\
			"addpd	%%xmm6,%%xmm3		\n\t"/* x0*y1 + fcy; xmm6 FREE */\
			"movaps	%%xmm1,%%xmm6		\n\t"/* xmm6 = cpy of x1 */\
			"mulpd	    (%%edx),%%xmm1	\n\t"/* x1*y0 */\
			"addpd	%%xmm3,%%xmm1		\n\t"/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */\
			"movaps	%%xmm1,%%xmm3		\n\t"/* fcy = cpy of fprod1 */\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm7,%%xmm3		\n\t"/* fcy */\
			"subpd	%%xmm3,%%xmm1		\n\t"/* fprod1 -= fcy */\
			"mulpd	%%xmm5,%%xmm1		\n\t"/* fprod1 *= two26f */\
			"mulpd	    (%%esi),%%xmm3	\n\t"/* fcy    *= two26i */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"mulpd	    (%%edx),%%xmm2	\n\t"/* x2*y0 */\
			"mulpd	0x10(%%edx),%%xmm6	\n\t"/* x1*y1 */\
			"mulpd	0x20(%%edx),%%xmm4	\n\t"/* x0*y2 */\
			"addpd	%%xmm6,%%xmm2		\n\t"/* x1*y1 + x2*y0; xmm6 FREE */\
			"addpd	%%xmm4,%%xmm3		\n\t"/* x0*y2 + fcy; xmm4 FREE */\
			"addpd	%%xmm3,%%xmm2		\n\t"/* fprod2; xmm3 FREE */\
			"movaps	%%xmm2,%%xmm3		\n\t"/* fcy = cpy of fprod2 */\
			"subpd	0x20(%%esi),%%xmm3	\n\t"	/* fprod2 - 0.5 */\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm7,%%xmm3		\n\t"/* fcy */\
			"subpd	%%xmm3,%%xmm2		\n\t"/* fprod2 -= fcy */\
			"mulpd	%%xmm5,%%xmm2		\n\t"/* fprod2 *= two26f */\
		"/* MULH96(q,lo,lo) --> lo =(q*lo)/2^78 */\n\t"\
			"movl	%[__fq0],%%edx		\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm3		\n\t"/* xmm3 = cpy of y0 */\
			"movaps	%%xmm0,%%xmm4		\n\t"/* xmm4 = cpy of y0 */\
			"mulpd	    (%%edx),%%xmm0	\n\t"/* fprod0 = y0*x0 */\
			"mulpd	    (%%esi),%%xmm0	\n\t"/* fcy    *= two26i */\
			"/* Digit 1: */\n\t"\
			"movaps	%%xmm1,%%xmm6		\n\t"/* xmm6 = cpy of y1 */\
			"mulpd	0x10(%%edx),%%xmm3	\n\t"/* y0*x1 */\
			"addpd	%%xmm0,%%xmm3		\n\t"/* y0*x1 + fcy; xmm0 FREE */\
			"mulpd	    (%%edx),%%xmm1	\n\t"/* y1*x0 */\
			"addpd	%%xmm3,%%xmm1		\n\t"/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */\
			"mulpd	    (%%esi),%%xmm1	\n\t"/* fcy    *= two26i */\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"movaps	%%xmm2,%%xmm3		\n\t"/* xmm3 = cpy of y2 */\
			"movaps	%%xmm6,%%xmm0		\n\t"/* xmm0 = cpy of y1 */\
			"mulpd	    (%%edx),%%xmm2	\n\t"/* y2*x0 */\
			"mulpd	0x10(%%edx),%%xmm6	\n\t"/* y1*x1 */\
			"mulpd	0x20(%%edx),%%xmm4	\n\t"/* y0*x2 */\
			"addpd	%%xmm6,%%xmm2		\n\t"/* y1*x1 + y2*x0; xmm6 FREE */\
			"addpd	%%xmm4,%%xmm1		\n\t"/* y0*x2 + fcy; xmm4 FREE */\
			"addpd	%%xmm1,%%xmm2		\n\t"/* fprod2; xmm1 FREE */\
			"subpd	0x20(%%esi),%%xmm2	\n\t"/* fprod2 - 0.5 */\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"/* fcy */\
			"mulpd	    (%%esi),%%xmm2	\n\t"/* fcy    *= two26i */\
			"/* xmm1,3,4,6 FREE; Precompute all the needed partial products: */\n\t"\
			"movaps	%%xmm3,%%xmm1		\n\t"/* xmm1 = cpy of y2 */\
			"mulpd	0x20(%%edx),%%xmm0	\n\t"/* y1*x2 */\
			"mulpd	0x10(%%edx),%%xmm1	\n\t"/* y2*x1 */\
			"mulpd	0x20(%%edx),%%xmm3	\n\t"/* y2*x2; xmm4,6 FREE */\
			"/* Digit 3: */\n\t"\
			"addpd	%%xmm2,%%xmm1		\n\t"/* fy2*fx1 + fcy */\
			"addpd	%%xmm1,%%xmm0		\n\t"/* fprod3 = fy1*fx2 + fy2*fx1 + fcy */\
			"movaps	%%xmm0,%%xmm6		\n\t"/* fcy = cpy of fprod3 */\
			"addpd	%%xmm7,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm6		\n\t"/* fcy */\
			"subpd	%%xmm6,%%xmm0		\n\t"/* fprod3 -= fcy */\
			"mulpd	%%xmm5,%%xmm0		\n\t"/* fprod3 *= two26f */\
			"mulpd	    (%%esi),%%xmm6	\n\t"/* fcy    *= two26i; xmm1,2,4 FREE */\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm6,%%xmm3		\n\t"/* fcy = fy2*fx2 + fcy */\
			"movaps	%%xmm3,%%xmm1		\n\t"/* fprod4 = cpy of fcy */\
			"addpd	%%xmm7,%%xmm3		\n\t"\
			"subpd	%%xmm7,%%xmm3		\n\t"/* fcy */\
			"subpd	%%xmm3,%%xmm1		\n\t"/* fprod4 -= fcy */\
			"mulpd	%%xmm5,%%xmm1		\n\t"/* fprod4 *= two26f */\
			"/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\n\t"\
			"movaps	0x30(%%eax),%%xmm2	\n\t"\
			"movaps	%%xmm2,%%xmm6		\n\t"\
			"movaps	0x30(%%edx),%%xmm4	\n\t"\
			"mulpd	%%xmm5,%%xmm3		\n\t"\
			"mulpd	    (%%edx),%%xmm5	\n\t"\
			"addpd	%%xmm1,%%xmm3		\n\t"\
			"movaps	0x30(%%eax),%%xmm1	\n\t"\
			"cmppd	$0x1,%%xmm3,%%xmm6	\n\t"\
			"subpd	%%xmm0,%%xmm1		\n\t"\
			"subpd	%%xmm3,%%xmm2		\n\t"\
			"movaps	%%xmm4,%%xmm0		\n\t"\
			"andpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm4,%%xmm2		\n\t"\
			"andpd	%%xmm5,%%xmm6		\n\t"\
			"addpd	%%xmm6,%%xmm1		\n\t"\
			"/* qlo26 in xmm5, qhi52 in xmm0 */\n\t"\
			"/* if((pshift >> j) & (uint64)1) { */\n\t"\
			"movl	%[__pshift],%%eax	\n\t"\
			"movl	%[__j],%%ecx		\n\t"\
			"shrl	%%cl,%%eax		\n\t"\
			"andl	$0x1,%%eax		\n\t"\
			"/* Branchless verbxon of the conditional doubling: */\n\t"\
			"xorl	$0x1,%%eax		\n\t"\
			"movd	%%eax,%%xmm6		\n\t"\
			"pshufd	$0,%%xmm6,%%xmm6	\n\t"\
			"cvtdq2pd	%%xmm6,%%xmm6	\n\t"\
			"xorpd	%%xmm3,%%xmm3		\n\t"\
			"cmppd	$0,%%xmm3,%%xmm6	\n\t"\
			"/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */\n\t"\
			"movaps	%%xmm2,%%xmm4		\n\t"\
			"movaps	%%xmm1,%%xmm3		\n\t"\
			"andpd	%%xmm6,%%xmm4		\n\t"\
			"andpd	%%xmm6,%%xmm3		\n\t"\
			"addpd	%%xmm4,%%xmm2		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"/* If x > q, subtract q: */\n\t"\
			"movaps	%%xmm0,%%xmm6		\n\t"\
			"cmppd	$0x2,%%xmm2,%%xmm6	\n\t"\
			"andpd	%%xmm6,%%xmm0		\n\t"\
			"andpd	%%xmm6,%%xmm5		\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"/* } */\n\t"\
			"/* Normalize the result: */\n\t"\
			"movaps	-0x10(%%esi),%%xmm4	\n\t"\
			"movaps	     (%%esi),%%xmm3	\n\t"\
			"mulpd	%%xmm3,%%xmm1		\n\t"\
			"mulpd	%%xmm3,%%xmm2		\n\t"\
			"movaps	%%xmm1,%%xmm0		\n\t"\
			"addpd	%%xmm7,%%xmm1		\n\t"\
			"subpd	%%xmm7,%%xmm1		\n\t"\
			"subpd	%%xmm1,%%xmm0		\n\t"\
			"mulpd	%%xmm4,%%xmm0		\n\t"\
			"mulpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm1,%%xmm2		\n\t"\
			"movaps	%%xmm2,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm2,%%xmm1		\n\t"\
			"mulpd	%%xmm4,%%xmm1		\n\t"\
			"/* Move high 2 words of result into input registers expected by start of loop body: */\n\t"\
			"movaps	%%xmm2,%%xmm4	/* fx2 */\n\t"\
			"movaps	%%xmm1,%%xmm2	/* fx1 */\n\t"\
	"/* GCC builds require us to save the xmm0/2/4 registers on each pass through the loop: */\n\t"\
		"movl	%[__fx0],%%eax	 		\n\t"\
		"movaps	%%xmm0,     (%%eax)		\n\t"\
		"movaps	%%xmm2,0x10(%%eax)		\n\t"\
		"movaps	%%xmm4,0x20(%%eax)		\n\t"\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			 ,[__fx0]	 "m" (fx0)		\
			 ,[__fq0]	 "m" (fq0)		\
			 ,[__pshift] "m" (pshift)	\
			 ,[__j]		 "m" (j)		\
			: "cl","eax","esi","ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
		);

		#else	// The 64-bit version of the macro is timing-suboptimal because I used it as a testbed:
				// This 2-TF-input/4-xxm-register version serves as the basis for an 8-input version
				// using all 16 xmm-regs available in 64-bit mode. Doing 4 independent 4-register-using
				// instruction streams in parallel should allo for better latency hiding than we get with
				// original 2-input/8-register implementation and its translation to 4-input/16-register,
				// which can be seen in the _q4 version of this routine.

		__asm__ volatile (\
		/* SQR_LOHI78_3WORD_DOUBLE_q2(fx, flo,fhi). Inputs flo0,1,2 enter in xmm0,1,2: */\
			"movq	%[__fx0],%%rax		\n\t"\
			"movq	%[__two26i],%%rbx	\n\t"\
			"movaps	-0x20(%%rbx),%%xmm3	\n\t"/* two13i */\
			"/* fx0,1,2 assumed in xmm0,1,2 on loop entry */\n\t"\
			"mulpd	%%xmm3,%%xmm0		\n\t"/* scale [fx0 . 2^-13] */\
			"mulpd	%%xmm3,%%xmm1		\n\t"/* scale [fx1 . 2^-13] */\
			"mulpd	%%xmm3,%%xmm2		\n\t"/* scale [fx2 . 2^-13] */\
			"movaps	%%xmm0,%%xmm3		\n\t"/* cpy of fx0 */\
			"addpd	%%xmm3,%%xmm3		\n\t"/* 2.fx0 */\
			"movaps	%%xmm2,0x20(%%rax)	\n\t"/* Store fx2 to free up a register */\
			"/* Digit 0: */\n\t"\
			"mulpd	%%xmm0,%%xmm0		\n\t"/* [  fx0 *= fx0] / 2^26 */\
			"roundpd $0,%%xmm0,%%xmm2	\n\t"/* Init: FCY = dnint(fx0*fx0) */\
			"subpd	%%xmm2,%%xmm0		\n\t"/* fprod0 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fprod0 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* FCY    *= two26i */\
			"movaps	%%xmm0,    (%%rax)	\n\t"/* Store fprod0 to free up a register */\
			"/* Digit 1: */\n\t"\
			"movaps	%%xmm1,%%xmm0		\n\t"/* cpy of fx1 */\
			"mulpd	%%xmm3,%%xmm0		\n\t"/* [fx1 *= 2.fx0] / 2^26 */\
			"addpd	%%xmm2,%%xmm0		\n\t"/* fprod1 = 2.fx0*fx1 + FCY */\
			"roundpd $0,%%xmm0,%%xmm2	\n\t"/* FCY */\
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
			"roundpd $1,%%xmm1,%%xmm2	\n\t"/* FCY = floor(fprod2) */\
			"subpd	%%xmm2,%%xmm1		\n\t"/* fprod2 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm1	\n\t"/* fprod2 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* FCY    *= two26i */\
			"/* Digit 3: */\n\t"\
			"movaps	0x20(%%rax),%%xmm3	\n\t"/* Reload fx2 */\
			"mulpd	%%xmm3,%%xmm0		\n\t"/* [2.fx1 *= fx2] / 2^26 */\
			"addpd	%%xmm2,%%xmm0		\n\t"/* fprod3 = 2.fx1*fx2 + FCY */\
			"roundpd $0,%%xmm0,%%xmm2	\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t"/* fprod3 -= FCY */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fprod3 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* FCY    *= two26i */\
			"/* Digit 4: */\n\t"\
			"mulpd	%%xmm3,%%xmm3		\n\t"/* [  fx2 *= fx2] / 2^26 */\
			"addpd	%%xmm2,%%xmm3		\n\t"/* fprod4 = fx2*fx2 + fcy */\
			"roundpd $0,%%xmm3,%%xmm2	\n\t"\
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
			"roundpd $0,%%xmm0,%%xmm1	\n\t"\
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
			"roundpd $0,%%xmm1,%%xmm3	\n\t"\
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
			"roundpd $1,%%xmm2,%%xmm3	\n\t"/* fcy = floor(fprod2) */\
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
			"roundpd $1,%%xmm2,%%xmm2	\n\t"/* CY = floor(fprod2) */\
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
			"roundpd $0,%%xmm0,%%xmm2	\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t"/* fprod3 -= CY */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fprod3 *= two26f */\
			"mulpd	     (%%rbx),%%xmm2	\n\t"/* CY    *= two26i */\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm2,%%xmm3		\n\t"/* CY = y2*x2 + CY */\
			"movaps	%%xmm3,%%xmm1		\n\t"/* fprod4 = cpy of CY */\
			"roundpd $0,%%xmm3,%%xmm3	\n\t"\
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
			"movslq	%[__pshift],%%rax	\n\t"\
			"movslq	%[__j],%%rcx		\n\t"\
			"shrq	%%cl,%%rax		\n\t"\
			"andq	$0x1,%%rax		\n\t"\
		"je twopmodq78_3wdq2_gcc64	\n\t"\
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
			"roundpd $0,%%xmm1,%%xmm1	\n\t"/* fcy */\
			"subpd	%%xmm1,%%xmm0		\n\t"/* fx0 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm0	\n\t"/* fx0 *= two26f */\
			"mulpd	%%xmm3,%%xmm1		\n\t"/* fcy *= two26i */\
			"addpd	%%xmm1,%%xmm2		\n\t"/* Add carry into xhi */\
			"movaps	%%xmm2,%%xmm1		\n\t"/* Init: fcy = cpy of ~xhi */\
			"roundpd $0,%%xmm2,%%xmm2	\n\t"/* fx2 = fcy (no divide by 2^26) */\
			"subpd	%%xmm2,%%xmm1		\n\t"/* fx1 -= fcy */\
			"mulpd	-0x10(%%rbx),%%xmm1	\n\t"/* fx1 *= two26f */\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			 ,[__fx0]	 "m" (fx0)		\
			 ,[__fq0]	 "m" (fq0)		\
			 ,[__pshift] "m" (pshift)	\
			 ,[__j]		 "m" (j)		\
			: "cl","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3"	/* Clobbered registers */\
		);

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
		}	/* for(j...) */

	#endif	/* COMPILER_TYPE_ */

  #endif /* if(0) */

#endif	/* USE_SSE2 */

#ifdef DEBUG_SSE2
//	exit(0);
#endif

#ifdef USE_SSE2

	// Dump SSE2 registers to memory: */
	#if defined(COMPILER_TYPE_MSVC)

		__asm	mov	eax, fx0
		__asm	movaps	[eax     ],xmm0	/* flo0 */
		__asm	movaps	[eax+0x10],xmm2	/* flo1 */
		__asm	movaps	[eax+0x20],xmm4	/* flo2 */


	#else	/* GCC-style inline ASM: */

		/* GCC builds require us to save the xmm0/2/4 registers on each pass through the loop,
		so skip these pre-loop-entry snippets until we get the full loop-logic asm'ed: */
		#if OS_BITS == 32
			__asm__ volatile (\
				"movl	%[__fx0],%%eax 	\n\t"\
				"movaps	%%xmm0,    (%%eax)	\n\t"\
				"movaps	%%xmm2,0x10(%%eax)	\n\t"\
				"movaps	%%xmm4,0x20(%%eax)	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "cc","memory","eax","xmm0","xmm2","xmm4"	/* Clobbered registers */\
			);
		#else	// ^4-bit mode has in/outputs in xmm0,1,2:
			__asm__ volatile (\
				"movq	%[__fx0],%%rax 	\n\t"\
				"movaps	%%xmm0,    (%%rax)	\n\t"\
				"movaps	%%xmm1,0x10(%%rax)	\n\t"\
				"movaps	%%xmm2,0x20(%%rax)	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "cc","memory","rax","xmm0","xmm1","xmm2"	/* Clobbered registers */\
			);
		#endif

	#endif

#endif


	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
#ifdef USE_SSE2
	CVT78_3WORD_DOUBLE_UINT96(*fx0,*fx1,*fx2, x0);
	CVT78_3WORD_DOUBLE_UINT96(*gx0,*gx1,*gx2, x1);
#else
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
	CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
#endif
	ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
	ADD96(x1,x1,x1);

	SUB96(x0,q0,x0);
	SUB96(x1,q1,x1);

	tmp0 = CMPEQ96(x0, ONE96);	*checksum2 += x0.d0;
	tmp1 = CMPEQ96(x1, ONE96);	*checksum2 += x1.d0;
	r = tmp0;
	r += tmp1 << 1;
	return r;
}

#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
	#define YES_SSE2
#endif

/*** 4-trial-factor version ***/
uint64 twopmodq78_3WORD_DOUBLE_q4(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
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
	static uint64 psave = 0;
	static uint32 pshift;	/* Require this to be 32-bit here for 32-bit ASM support */
	static uint32 start_index, zshift, first_entry = TRUE;

#ifdef YES_SSE2

	static struct complex *sc_arr = 0x0;
	static         double *sc_ptr;
	static double *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2, *flo0,*flo1,*flo2, *fhi0,*fhi1,*fhi2;
	static double *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2, *glo0,*glo1,*glo2, *ghi0,*ghi1,*ghi2;
	static double *hq0,*hq1,*hq2,*hqhi52, *hqinv0,*hqinv1,*hqinv2, *hx0,*hx1,*hx2, *hlo0,*hlo1,*hlo2, *hhi0,*hhi1,*hhi2;
	static double *iq0,*iq1,*iq2,*iqhi52, *iqinv0,*iqinv1,*iqinv2, *ix0,*ix1,*ix2, *ilo0,*ilo1,*ilo2, *ihi0,*ihi1,*ihi2;
	static double *half, *two26f, *two26i, *two13i, *sse2_rnd;
	static uint64 ihalf = 0x3FDfffffffffffffull;	/* Bitfield storing 0.5*(1-epsilon) in IEEE64 format */

#else

	double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2,fqhalf0, fx0,fx1,fx2, flo0,flo1,flo2,flohi52, fhi0,fhi1,fhi2, ftmp0,ftmp1,ftmp2;
	double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2,fqhalf1, gx0,gx1,gx2, glo0,glo1,glo2,glohi52, ghi0,ghi1,ghi2, gtmp0,gtmp1,gtmp2;
	double hq0,hq1,hq2, hqinv0,hqinv1,hqinv2,fqhalf2, hx0,hx1,hx2, hlo0,hlo1,hlo2,hlohi52, hhi0,hhi1,hhi2, htmp0,htmp1,htmp2;
	double iq0,iq1,iq2, iqinv0,iqinv1,iqinv2,fqhalf3, ix0,ix1,ix2, ilo0,ilo1,ilo2,ilohi52, ihi0,ihi1,ihi2, itmp0,itmp1,itmp2;
	double fq_or_nil_lo26[2], fq_or_nil_hi52[2];
	double gq_or_nil_lo26[2], gq_or_nil_hi52[2];
	double hq_or_nil_lo26[2], hq_or_nil_hi52[2];
	double iq_or_nil_lo26[2], iq_or_nil_hi52[2];
	int fidx,gidx,hidx,iidx;

#endif

#ifdef YES_SSE2

	if(first_entry)
	{
		sc_arr = ALLOC_COMPLEX(sc_arr, 80);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = (double *)ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		/* Remember, rhese are POINTERS-TO-DOUBLES, so need an increment of 2 to span an SSE register: */
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
		*sse2_rnd++ = 3.0*0x4000000*0x2000000;
		*sse2_rnd-- = 3.0*0x4000000*0x2000000;
		/* We init "half" = 0.5-epsilon here, because emulating FLOOR(x) via DNINT(x-half) requires
		us to always round up if x is a whole number, and our DNINT emulation can round either way if fractional part = 0.5:
		*/
		*half++ = *(double*)&ihalf;	*half-- = *(double*)&ihalf;
	}	/* first_entry */

#endif

	ASSERT(HERE, (p >> 63) == 0, "twopmodq78_q4 : p must be < 2^63!");
	q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(q0.d0, k0,&q0.d0,&q0.d1);
	MUL_LOHI64(q1.d0, k1,&q1.d0,&q1.d1);
	MUL_LOHI64(q2.d0, k2,&q2.d0,&q2.d1);
	MUL_LOHI64(q3.d0, k3,&q3.d0,&q3.d1);
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
	ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q4 : (q0.d1 >> 14) != 0");
	ASSERT(HERE, (q1.d1 >> 14) == 0, "twopmodq78_q4 : (q1.d1 >> 14) != 0");
	ASSERT(HERE, (q2.d1 >> 14) == 0, "twopmodq78_q4 : (q2.d1 >> 14) != 0");
	ASSERT(HERE, (q3.d1 >> 14) == 0, "twopmodq78_q4 : (q3.d1 >> 14) != 0");

	*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0;

	q32_0 = (uint32)q0.d0;
	q32_1 = (uint32)q1.d0;
	q32_2 = (uint32)q2.d0;
	q32_3 = (uint32)q3.d0;

	/* Convert q to floating form: */
#ifdef YES_SSE2

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

#ifndef YES_SSE2
	fqhalf0 = qhalf0.d1*TWO64FLOAT + qhalf0.d0;
	fqhalf1 = qhalf1.d1*TWO64FLOAT + qhalf1.d0;
	fqhalf2 = qhalf2.d1*TWO64FLOAT + qhalf2.d0;
	fqhalf3 = qhalf3.d1*TWO64FLOAT + qhalf3.d0;
#endif

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		ASSERT(HERE, ((p + 78) >> 32) == 0, "twopmodq78_q2 : (p+78) exceeds 32 bits!");
		pshift = (uint32)p + 78;
		j = leadz32(pshift);
		lead7 = ((pshift<<j) >> 25);
		if(lead7 > 77)
		{
			lead7 >>= 1;
			start_index =  32-j-6;	/* Use only the leftmost 6 bits */
		}
		else
			start_index =  32-j-7;

		zshift = 77 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

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
#ifdef YES_SSE2

	CVT_UINT78_3WORD_DOUBLE(qinv0 ,*fqinv0,*fqinv1,*fqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv1 ,*gqinv0,*gqinv1,*gqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv2 ,*hqinv0,*hqinv1,*hqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv3 ,*iqinv0,*iqinv1,*iqinv2);

	#if defined(COMPILER_TYPE_MSVC) && (TRYQ == 4)

		#error twopmodq78_3WORD_DOUBLE_q4: 4-operand SSE2 modmul does not support MSVC inline assembler!

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#error twopmodq78_3WORD_DOUBLE_q4: 4-operand SSE2 modmul requires 64-bit GCC build!

		#else

		__asm__ volatile (\
			"movq	%[__fqinv0],%%rax 	\n\t"\
			"movq	%[__two26i],%%rsi 	\n\t"\
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
#ifdef YES_SSE2

	CVT_UINT78_3WORD_DOUBLE(x0 ,*fx0,*fx1,*fx2);
	CVT_UINT78_3WORD_DOUBLE(x1 ,*gx0,*gx1,*gx2);
	CVT_UINT78_3WORD_DOUBLE(x2 ,*hx0,*hx1,*hx2);
	CVT_UINT78_3WORD_DOUBLE(x3 ,*ix0,*ix1,*ix2);

	__asm__ volatile (\
		"movq	%[__fx0],%%rax 	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t	movaps	0x10(%%rax),%%xmm8	\n\t"\
		"movaps	0x20(%%rax),%%xmm2	\n\t	movaps	0x30(%%rax),%%xmm10	\n\t"\
		"movaps	0x40(%%rax),%%xmm4	\n\t	movaps	0x50(%%rax),%%xmm12	\n\t"\
		:				/* outputs: none */\
		: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
	);

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
#ifndef YES_SSE2

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
	}

#elif(defined(YES_SSE2))	/* Inner loop body needs 42 movaps, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */

	#define ASM_LOOP 0

	#if ASM_LOOP	/* Pure-ASM loop control, and branchless mod-doubling segment: */
		__asm__ volatile (\
			"movslq	%[__start_index], %%rcx		\n\t"\
			"subq $2,%%rcx						\n\t"\
			"test %%rcx, %%rcx					\n\t"\
			"jl Loop4End		/* Skip if n < 0 */	\n\t"\
		"Loop4Beg:								\n\t"\
			"/* SQR_LOHI78_3WORD_DOUBLE_q4(): */\n\t"\
			"movq	%[__fx0],%%rax				\n\t"\
			"movq	%[__two26i],%%rbx			\n\t"\
			"movaps	-0x20(%%rbx),%%xmm6	/* two13i */	\n\t"\
			"/* fx0,1,2 assumed in xmm0,2,4 on loop entry */\n\t"\
			"mulpd	%%xmm6,%%xmm0			\n\t	mulpd	%%xmm6 ,%%xmm8			\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t	mulpd	%%xmm6 ,%%xmm10			\n\t"\
			"mulpd	%%xmm6,%%xmm4			\n\t	mulpd	%%xmm6 ,%%xmm12			\n\t"\
			"movaps	0x10(%%rbx),%%xmm7		/* xmm7 = rnd_const shared between both columns*/\n\t"\
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
			"/* Move this part of Digit 2 computation here to free up xmm5,13: */	\n\t"\
			"addpd	%%xmm5,%%xmm2			\n\t	addpd	%%xmm13,%%xmm10			\n\t"\
			"/* Digit 0: */					\n\t"\
			"movaps (%%rbx),%%xmm15			\n\t	movaps	-0x10(%%rbx),%%xmm5		\n\t"\
			"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
			"mulpd	%%xmm5,%%xmm0		\n\t	mulpd	%%xmm5,%%xmm8		\n\t"\
			"mulpd	     %%xmm15,%%xmm6		\n\t	mulpd	     %%xmm15,%%xmm14	\n\t"\
			"/* Digit 1: */					\n\t"\
			"addpd	%%xmm6,%%xmm1			\n\t	addpd	%%xmm14,%%xmm9			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
			"mulpd	%%xmm5,%%xmm1		\n\t	mulpd	%%xmm5,%%xmm9		\n\t"\
			"mulpd	     %%xmm15,%%xmm6		\n\t	mulpd	     %%xmm15,%%xmm14	\n\t"\
			"/* Digit 2: */					\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t"\
			"movaps	%%xmm2,%%xmm6			\n\t	movaps	%%xmm10,%%xmm14			\n\t"\
			"subpd	0x20(%%rbx),%%xmm6		\n\t	subpd	0x20(%%rbx),%%xmm14		\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t	subpd	%%xmm14,%%xmm10			\n\t"\
			"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"/* Digit 3: */					\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
			"movaps	%%xmm3,%%xmm6			\n\t	movaps	%%xmm11,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t	subpd	%%xmm14,%%xmm11			\n\t"\
			"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"/* Digit 4: */					\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm4			\n\t	subpd	%%xmm14,%%xmm12			\n\t"\
			"mulpd	%%xmm5,%%xmm4			\n\t	mulpd	%%xmm5 ,%%xmm12			\n\t"\
			"/* Digit 5 = the carry. flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6 */	\n\t"\
			"mulpd	%%xmm5,%%xmm6			\n\t	mulpd	%%xmm5 ,%%xmm14			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t"\
			"movaps	%%xmm3,0xc0(%%rax)		\n\t	movaps	%%xmm11,0xd0(%%rax)		\n\t"\
			"movaps	%%xmm4,0xe0(%%rax)		\n\t	movaps	%%xmm12,0xf0(%%rax)		\n\t"\
			"/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */						\n\t"\
			"movq	%[__fqinv0],%%rdx		\n\t"\
			"/* Digit 0: */					\n\t"\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t"\
			"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t"\
			"mulpd	    (%%rdx),%%xmm0		\n\t	mulpd	0x10(%%rdx),%%xmm8		\n\t"\
			"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
			"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t"\
			"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t"\
			"/* Digit 1: */					\n\t"\
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
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t"\
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
			"/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */								\n\t"\
			"movq	%[__fq0],%%rdx		"
			"/* Digit 0: */					\n\t"\
			"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t"\
			"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t"\
			"mulpd	    (%%rdx),%%xmm0		\n\t	mulpd	0x10(%%rdx),%%xmm8		\n\t"\
			"addpd	%%xmm7,%%xmm0			\n\t	addpd	%%xmm7 ,%%xmm8			\n\t"\
			"subpd	%%xmm7,%%xmm0			\n\t	subpd	%%xmm7 ,%%xmm8			\n\t"\
			"mulpd	    %%xmm15,%%xmm0		\n\t	mulpd	    %%xmm15,%%xmm8		\n\t"\
			"/* Digit 1: */					\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
			"mulpd	0x20(%%rdx),%%xmm3		\n\t	mulpd	0x30(%%rdx),%%xmm11		\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t	addpd	%%xmm8 ,%%xmm11			\n\t"\
			"mulpd	    (%%rdx),%%xmm1		\n\t	mulpd	0x10(%%rdx),%%xmm9		\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
			"addpd	%%xmm7,%%xmm1			\n\t	addpd	%%xmm7 ,%%xmm9			\n\t"\
			"subpd	%%xmm7,%%xmm1			\n\t	subpd	%%xmm7 ,%%xmm9			\n\t"\
			"mulpd	    %%xmm15,%%xmm1		\n\t	mulpd	    %%xmm15,%%xmm9		\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t"\
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
			"/* Precompute all the needed partial products: */						\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t"\
			"mulpd	0x40(%%rdx),%%xmm0		\n\t	mulpd	0x50(%%rdx),%%xmm8		\n\t"\
			"mulpd	0x20(%%rdx),%%xmm1		\n\t	mulpd	0x30(%%rdx),%%xmm9		\n\t"\
			"mulpd	0x40(%%rdx),%%xmm3		\n\t	mulpd	0x50(%%rdx),%%xmm11		\n\t"\
			"/* Digit 3: */					\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t	addpd	%%xmm10,%%xmm9			\n\t"\
			"addpd	%%xmm1,%%xmm0			\n\t	addpd	%%xmm9 ,%%xmm8			\n\t"\
			"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
			"addpd	%%xmm7,%%xmm6			\n\t	addpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm7,%%xmm6			\n\t	subpd	%%xmm7 ,%%xmm14			\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t	subpd	%%xmm14,%%xmm8			\n\t"\
			"mulpd	%%xmm5,%%xmm0			\n\t	mulpd	%%xmm5 ,%%xmm8			\n\t"\
			"mulpd	    %%xmm15,%%xmm6		\n\t	mulpd	    %%xmm15,%%xmm14		\n\t"\
			"/* Digit 4: */					\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t	addpd	%%xmm14,%%xmm11			\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t	movaps	%%xmm11,%%xmm9			\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t	addpd	%%xmm7 ,%%xmm11			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t	subpd	%%xmm7 ,%%xmm11			\n\t"\
			"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t"\
			"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t"\
			"/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\n\t"\
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
			"/* qlo26 in xmm5, qhi52 in xmm0 */		\n\t"\
			"/* if((pshift >> j) & (uint64)1) { */	\n\t"\
			"movslq	%[__pshift],%%rax		\n\t"\
			"shrq	%%cl,%%rax				\n\t"\
			"andq	$0x1,%%rax				\n\t"\
			"/* Branchless version of the conditional doubling: */\n\t"\
			"xorq	$0x1,%%rax				\n\t"\
			"movd	%%rax,%%xmm6			\n\t"\
			"pshufd	$0,%%xmm6,%%xmm6		\n\t"\
			"cvtdq2pd	%%xmm6,%%xmm6		\n\t"\
			"xorpd	%%xmm3,%%xmm3			\n\t"\
			"cmppd	$0x0,%%xmm3,%%xmm6		\n\t"\
			"/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm10,xmm9; qhi52,qlo26 in xmm8,xmm13 */\n\t"\
			"movaps	%%xmm2,%%xmm4			\n\t	movaps	%%xmm10,%%xmm12			\n\t"\
			"movaps	%%xmm1,%%xmm3			\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
			"andpd	%%xmm6,%%xmm4			\n\t	andpd	%%xmm6 ,%%xmm12			\n\t"\
			"andpd	%%xmm6,%%xmm3			\n\t	andpd	%%xmm6 ,%%xmm11			\n\t"\
			"addpd	%%xmm4,%%xmm2			\n\t	addpd	%%xmm12,%%xmm10			\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
			"/* If x > q, subtract q: */	\n\t"\
			"movaps	%%xmm0,%%xmm6			\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
			"cmppd	$0x2,%%xmm2,%%xmm6		\n\t	cmppd	$0x2,%%xmm10,%%xmm14	\n\t"\
			"andpd	%%xmm6,%%xmm0			\n\t	andpd	%%xmm14,%%xmm8			\n\t"\
			"andpd	%%xmm6,%%xmm5			\n\t	andpd	%%xmm14,%%xmm13			\n\t"\
			"subpd	%%xmm0,%%xmm2			\n\t	subpd	%%xmm8 ,%%xmm10			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t	subpd	%%xmm13,%%xmm9			\n\t"\
			"/* } */						\n\t"\
			"/* Normalize the result: */	\n\t"\
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
		"movq	%[__fx0],%%rax 	\n\t"\
		"movaps	%%xmm0,    (%%rax)	\n\t	movaps	%%xmm8 ,0x10(%%rax)	\n\t"\
		"movaps	%%xmm2,0x20(%%rax)	\n\t	movaps	%%xmm10,0x30(%%rax)	\n\t"\
		"movaps	%%xmm4,0x40(%%rax)	\n\t	movaps	%%xmm12,0x50(%%rax)	\n\t"\
		:				/* outputs: none */\
		: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
	);

#endif	/* YES_SSE2 */

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
#ifdef YES_SSE2
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

	SUB96(x0,q0,x0);
	SUB96(x1,q1,x1);
	SUB96(x2,q2,x2);
	SUB96(x3,q3,x3);

	tmp0 = CMPEQ96(x0, ONE96);	*checksum2 += x0.d0;
	tmp1 = CMPEQ96(x1, ONE96);	*checksum2 += x1.d0;
	tmp2 = CMPEQ96(x2, ONE96);	*checksum2 += x2.d0;
	tmp3 = CMPEQ96(x3, ONE96);	*checksum2 += x3.d0;
	r = tmp0;
	r += tmp1 << 1;
	r += tmp2 << 2;
	r += tmp3 << 3;
	return r;
}

/*** Older "reference" version ***/
uint64 twopmodq78_3WORD_DOUBLE_q4_REF(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
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

	ASSERT(HERE, (p >> 63) == 0, "twopmodq78_q4 : p must be < 2^63!");
	q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(q0.d0, k0,&q0.d0,&q0.d1);
	MUL_LOHI64(q1.d0, k1,&q1.d0,&q1.d1);
	MUL_LOHI64(q2.d0, k2,&q2.d0,&q2.d1);
	MUL_LOHI64(q3.d0, k3,&q3.d0,&q3.d1);
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
	ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q4 : (q0.d1 >> 14) != 0");
	ASSERT(HERE, (q1.d1 >> 14) == 0, "twopmodq78_q4 : (q1.d1 >> 14) != 0");
	ASSERT(HERE, (q2.d1 >> 14) == 0, "twopmodq78_q4 : (q2.d1 >> 14) != 0");
	ASSERT(HERE, (q3.d1 >> 14) == 0, "twopmodq78_q4 : (q3.d1 >> 14) != 0");

	*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0;

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

	SUB96(x0,q0,x0);
	SUB96(x1,q1,x1);
	SUB96(x2,q2,x2);
	SUB96(x3,q3,x3);

	*checksum2 += x0.d0 + x1.d0 + x2.d0 + x3.d0;

	r = 0;
	if(CMPEQ96(x0, ONE96)) r +=  1;
	if(CMPEQ96(x1, ONE96)) r +=  2;
	if(CMPEQ96(x2, ONE96)) r +=  4;
	if(CMPEQ96(x3, ONE96)) r +=  8;
	return(r);
}

#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)

  #ifdef USE_AVX	// Nov 2013: Now that have prototyped 2-TF-input/4-SSE-register modmul macro, use 8-input all-float in both SSE2 and AVX mode

	// AVX/8-input analog of twopmodq78_3WORD_DOUBLE_q4 :
	uint64 twopmodq78_3WORD_DOUBLE_q8(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
	{
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
		static uint64 psave = 0;
		static uint32 pshift;	/* Require this to be 32-bit here for 32-bit ASM support */
		static uint32 start_index, zshift, first_entry = TRUE;
		static vec_dbl *sc_arr = 0x0;
		static double *sc_ptr;
		double *tmp;
		static double *aq0,*aq1,*aq2,*aqhi52, *aqinv0,*aqinv1,*aqinv2, *ax0,*ax1,*ax2;
		static double *bq0,*bq1,*bq2,*bqhi52, *bqinv0,*bqinv1,*bqinv2, *bx0,*bx1,*bx2;
		static double *cq0,*cq1,*cq2,*cqhi52, *cqinv0,*cqinv1,*cqinv2, *cx0,*cx1,*cx2;
		static double *dq0,*dq1,*dq2,*dqhi52, *dqinv0,*dqinv1,*dqinv2, *dx0,*dx1,*dx2;
		static double *eq0,*eq1,*eq2,*eqhi52, *eqinv0,*eqinv1,*eqinv2, *ex0,*ex1,*ex2;
		static double *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2;
		static double *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2;
		static double *hq0,*hq1,*hq2,*hqhi52, *hqinv0,*hqinv1,*hqinv2, *hx0,*hx1,*hx2;
		static vec_dbl *two26f, *two26i, *two13i;

		if(first_entry)
		{
		#ifdef USE_AVX
			sc_arr = ALLOC_VEC_DBL(sc_arr, 28);
		#else	// SSE2:
			sc_arr = ALLOC_VEC_DBL(sc_arr, 52);
		#endif
			if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
			sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);
			ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			/* Remember, rhese are POINTERS-TO-DOUBLES, so need an increment of 4 to span an AVX register: */
			tmp = (double *)sc_ptr;
			aq0    = tmp + 0x00;	bq0    = tmp + 0x01;	cq0    = tmp + 0x02;	dq0    = tmp + 0x03;	eq0    = tmp + 0x04;	fq0    = tmp + 0x05;	gq0    = tmp + 0x06;	hq0    = tmp + 0x07;
			aq1    = tmp + 0x08;	bq1    = tmp + 0x09;	cq1    = tmp + 0x0a;	dq1    = tmp + 0x0b;	eq1    = tmp + 0x0c;	fq1    = tmp + 0x0d;	gq1    = tmp + 0x0e;	hq1    = tmp + 0x0f;
			aq2    = tmp + 0x10;	bq2    = tmp + 0x11;	cq2    = tmp + 0x12;	dq2    = tmp + 0x13;	eq2    = tmp + 0x14;	fq2    = tmp + 0x15;	gq2    = tmp + 0x16;	hq2    = tmp + 0x17;
			aqhi52 = tmp + 0x18;	bqhi52 = tmp + 0x19;	cqhi52 = tmp + 0x1a;	dqhi52 = tmp + 0x1b;	eqhi52 = tmp + 0x1c;	fqhi52 = tmp + 0x1d;	gqhi52 = tmp + 0x1e;	hqhi52 = tmp + 0x1f;
			aqinv0 = tmp + 0x20;	bqinv0 = tmp + 0x21;	cqinv0 = tmp + 0x22;	dqinv0 = tmp + 0x23;	eqinv0 = tmp + 0x24;	fqinv0 = tmp + 0x25;	gqinv0 = tmp + 0x26;	hqinv0 = tmp + 0x27;
			aqinv1 = tmp + 0x28;	bqinv1 = tmp + 0x29;	cqinv1 = tmp + 0x2a;	dqinv1 = tmp + 0x2b;	eqinv1 = tmp + 0x2c;	fqinv1 = tmp + 0x2d;	gqinv1 = tmp + 0x2e;	hqinv1 = tmp + 0x2f;
			aqinv2 = tmp + 0x30;	bqinv2 = tmp + 0x31;	cqinv2 = tmp + 0x32;	dqinv2 = tmp + 0x33;	eqinv2 = tmp + 0x34;	fqinv2 = tmp + 0x35;	gqinv2 = tmp + 0x36;	hqinv2 = tmp + 0x37;
			ax0    = tmp + 0x38;	bx0    = tmp + 0x39;	cx0    = tmp + 0x3a;	dx0    = tmp + 0x3b;	ex0    = tmp + 0x3c;	fx0    = tmp + 0x3d;	gx0    = tmp + 0x3e;	hx0    = tmp + 0x3f;
			ax1    = tmp + 0x40;	bx1    = tmp + 0x41;	cx1    = tmp + 0x42;	dx1    = tmp + 0x43;	ex1    = tmp + 0x44;	fx1    = tmp + 0x45;	gx1    = tmp + 0x46;	hx1    = tmp + 0x47;
			ax2    = tmp + 0x48;	bx2    = tmp + 0x49;	cx2    = tmp + 0x4a;	dx2    = tmp + 0x4b;	ex2    = tmp + 0x4c;	fx2    = tmp + 0x4d;	gx2    = tmp + 0x4e;	hx2    = tmp + 0x4f;
			// +0x50,58 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
		#ifdef USE_AVX
			two13i = (vec_dbl *)(sc_ptr + 0x60);
			two26f = (vec_dbl *)(sc_ptr + 0x64);
			two26i = (vec_dbl *)(sc_ptr + 0x68);	// Equivalent of 27 AVX slots, alloc 28
		#else	// SSE2:
			two13i = (vec_dbl *)(sc_ptr + 0x60);
			two26f = (vec_dbl *)(sc_ptr + 0x62);
			two26i = (vec_dbl *)(sc_ptr + 0x64);	// Equivalent of 51 SSE2 slots, alloc 52
		#endif
			/* Can premultiply each of the multiword-mul inputs by 1/sqrt(2^13) due to quadraticity: */
			VEC_DBL_INIT(two13i, TWO13FLINV);
			VEC_DBL_INIT(two26f, TWO26FLOAT);
			VEC_DBL_INIT(two26i, TWO26FLINV);
		}	/* first_entry */

		ASSERT(HERE, (p >> 63) == 0, "twopmodq78_q4 : p must be < 2^63!");
		q0.d0 = q1.d0 = q2.d0 = q3.d0 = q4.d0 = q5.d0 = q6.d0 = q7.d0 = p+p;
	#ifdef MUL_LOHI64_SUBROUTINE
		#error MUL_LOHI64_SUBROUTINE defined!
	#endif
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
		MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
		MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);
		MUL_LOHI64(q4.d0, k4, q4.d0, q4.d1);
		MUL_LOHI64(q5.d0, k5, q5.d0, q5.d1);
		MUL_LOHI64(q6.d0, k6, q6.d0, q6.d1);
		MUL_LOHI64(q7.d0, k7, q7.d0, q7.d1);
		/* Since 2*p*k even, no need to check for overflow here */
		q0.d0 += 1;
		q1.d0 += 1;
		q2.d0 += 1;
		q3.d0 += 1;
		q4.d0 += 1;
		q5.d0 += 1;
		q6.d0 += 1;
		q7.d0 += 1;
		ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q8 : (q0.d1 >> 14) != 0");
		ASSERT(HERE, (q1.d1 >> 14) == 0, "twopmodq78_q8 : (q1.d1 >> 14) != 0");
		ASSERT(HERE, (q2.d1 >> 14) == 0, "twopmodq78_q8 : (q2.d1 >> 14) != 0");
		ASSERT(HERE, (q3.d1 >> 14) == 0, "twopmodq78_q8 : (q3.d1 >> 14) != 0");
		ASSERT(HERE, (q4.d1 >> 14) == 0, "twopmodq78_q8 : (q4.d1 >> 14) != 0");
		ASSERT(HERE, (q5.d1 >> 14) == 0, "twopmodq78_q8 : (q5.d1 >> 14) != 0");
		ASSERT(HERE, (q6.d1 >> 14) == 0, "twopmodq78_q8 : (q6.d1 >> 14) != 0");
		ASSERT(HERE, (q7.d1 >> 14) == 0, "twopmodq78_q8 : (q7.d1 >> 14) != 0");

		*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0 + q4.d0 + q5.d0 + q6.d0 + q7.d0;

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

		if(first_entry || p != psave)
		{
			first_entry = FALSE;
			psave  = p;
			ASSERT(HERE, ((p + 78) >> 32) == 0, "twopmodq78_q2 : (p+78) exceeds 32 bits!");
			pshift = (uint32)p + 78;
			j = leadz32(pshift);
			lead7 = ((pshift<<j) >> 25);
			if(lead7 > 77)
			{
				lead7 >>= 1;
				start_index =  32-j-6;	/* Use only the leftmost 6 bits */
			}
			else
				start_index =  32-j-7;

			zshift = 77 - lead7;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

			pshift = ~pshift;
		}

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
			"movq	%[__aqinv0],%%rax 	\n\t"\
			"movq	%[__two26i],%%rsi 	\n\t"\
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
			"movq	%[__aqinv0],%%rax 	\n\t"\
			"movq	%[__two26i],%%rsi 	\n\t"\
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
		// Inline a (... >> 78) specilaization of the RSHIFT192 macro here, since have a constant shift count:
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
			"movq	%[__ax0],%%rax 	\n\t"\
			"vmovaps	    (%%rax),%%ymm0	\n\t	vmovaps	0x20(%%rax),%%ymm8	\n\t"\
			"vmovaps	0x40(%%rax),%%ymm2	\n\t	vmovaps	0x60(%%rax),%%ymm10	\n\t"\
			"vmovaps	0x80(%%rax),%%ymm4	\n\t	vmovaps	0xa0(%%rax),%%ymm12	\n\t"\
			:				/* outputs: none */\
			: [__ax0] "m" (ax0)		/* All inputs from memory addresses here */\
			: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
		);

	#else	// SSE2:

		__asm__ volatile (\
			"movq	%[__ax0],%%rax 	\n\t"\
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
			"movq	%[__ax0],%%rax 	\n\t"\
			"vmovaps	%%ymm0,    (%%rax)	\n\t	vmovaps	%%ymm8 ,0x20(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x40(%%rax)	\n\t	vmovaps	%%ymm10,0x60(%%rax)	\n\t"\
			"vmovaps	%%ymm4,0x80(%%rax)	\n\t	vmovaps	%%ymm12,0xa0(%%rax)	\n\t"\
			:				/* outputs: none */\
			: [__ax0] "m" (ax0)		/* All inputs from memory addresses here */\
			: "cc","memory","rax","xmm0","xmm2","xmm4","xmm8","xmm10","xmm12"	/* Clobbered registers */\
		);

	#else	// SSE2:

		__asm__ volatile (\
			"movq	%[__ax0],%%rax 	\n\t"\
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
			SUB96(x0,q0,x0);			SUB96(x4,q4,x4);
			SUB96(x1,q1,x1);			SUB96(x5,q5,x5);
			SUB96(x2,q2,x2);			SUB96(x6,q6,x6);
			SUB96(x3,q3,x3);			SUB96(x7,q7,x7);
		tmp0 = CMPEQ96(x0, ONE96);		tmp4 = CMPEQ96(x4, ONE96);
		tmp1 = CMPEQ96(x1, ONE96);		tmp5 = CMPEQ96(x5, ONE96);
		tmp2 = CMPEQ96(x2, ONE96);		tmp6 = CMPEQ96(x6, ONE96);
		tmp3 = CMPEQ96(x3, ONE96);		tmp7 = CMPEQ96(x7, ONE96);
		*checksum2 += x0.d0 + x4.d0;
		*checksum2 += x1.d0 + x5.d0;
		*checksum2 += x2.d0 + x6.d0;
		*checksum2 += x3.d0 + x7.d0;
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

  #else

	/************ 64-BIT SSE2 ONLY *** 8-trial-factor, hybrid int/float version **********************/
	/* The integer asm test-code is in twopmodq96.c, also wrapped in the same-named YES_ASM #define. */
	/* The floating-point code is the above 4-factor 64-bit SSE2 code, interleaved with the integer. */
	/*************************************************************************************************/
	uint64 twopmodq78_3WORD_DOUBLE_q8(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
	{
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
		static struct complex *sc_arr = 0x0;
		static         double *sc_ptr;
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
			sc_arr = ALLOC_COMPLEX(sc_arr, 40+32);	ASSERT(HERE, sc_arr != 0x0, "FATAL: unable to allocate sc_arr!");
			sc_ptr = (double *)ALIGN_COMPLEX(sc_arr);
			ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			/* Remember, rhese are POINTERS-TO-DOUBLES, so need an increment of 2 to span an SSE register.
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
			ASSERT(HERE, (uint32)sm_ptr == ((uint32)sc_ptr +  0x280), "sm_ptr not offset as expected!");
			/* Remember, rhese are pointers-to-uint128, so need an increment of 2 to span a memory slot: */											/* Byte offsets: */
			qptr4  = (uint96*)(sm_ptr + 0x00);	qptr5  = (uint96*)(sm_ptr + 0x02);	qptr6  = (uint96*)(sm_ptr + 0x04);	qptr7  = (uint96*)(sm_ptr + 0x06);	/* 0x280 */
			qinv4  = (uint96*)(sm_ptr + 0x08);	qinv5  = (uint96*)(sm_ptr + 0x0a);	qinv6  = (uint96*)(sm_ptr + 0x0c);	qinv7  = (uint96*)(sm_ptr + 0x0e);	/* 0x2c0 */
			x4     = (uint96*)(sm_ptr + 0x10);	x5     = (uint96*)(sm_ptr + 0x12);	x6     = (uint96*)(sm_ptr + 0x14);	x7     = (uint96*)(sm_ptr + 0x16);	/* 0x300 */
			lo4    = (uint96*)(sm_ptr + 0x18);	lo5    = (uint96*)(sm_ptr + 0x1a);	lo6    = (uint96*)(sm_ptr + 0x1c);	lo7    = (uint96*)(sm_ptr + 0x1e);	/* 0x340 */
			qhalf4 = (uint96*)(sm_ptr + 0x20);	qhalf5 = (uint96*)(sm_ptr + 0x22);	qhalf6 = (uint96*)(sm_ptr + 0x24);	qhalf7 = (uint96*)(sm_ptr + 0x26);	/* 0x380 */
			hi4    = (uint96*)(sm_ptr + 0x28);	hi5    = (uint96*)(sm_ptr + 0x2a);	hi6    = (uint96*)(sm_ptr + 0x2c);	hi7    = (uint96*)(sm_ptr + 0x2e);	/* 0x3c0 */
			ONE96_PTR = (uint96*)(sm_ptr + 0x30);
			ptr64 = (uint64*)ONE96_PTR;	*ptr64++ = ONE96.d0;	*ptr64-- = ONE96.d1;
		}	/* first_entry */

		ASSERT(HERE, (p >> 63) == 0, "p must be < 2^63!");
		q0.d0 = q1.d0 = q2.d0 = q3.d0 = q4.d0 = q5.d0 = q6.d0 = q7.d0 = p+p;
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
		MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
		MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);
		MUL_LOHI64(q4.d0, k4, q4.d0, q4.d1);
		MUL_LOHI64(q5.d0, k5, q5.d0, q5.d1);
		MUL_LOHI64(q6.d0, k6, q6.d0, q6.d1);
		MUL_LOHI64(q7.d0, k7, q7.d0, q7.d1);

		q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
		q1.d0 += 1;
		q2.d0 += 1;
		q3.d0 += 1;
		q4.d0 += 1;
		q5.d0 += 1;
		q6.d0 += 1;
		q7.d0 += 1;
		ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q8 : (q0.d1 >> 14) != 0");
		ASSERT(HERE, (q1.d1 >> 14) == 0, "twopmodq78_q8 : (q1.d1 >> 14) != 0");
		ASSERT(HERE, (q2.d1 >> 14) == 0, "twopmodq78_q8 : (q2.d1 >> 14) != 0");
		ASSERT(HERE, (q3.d1 >> 14) == 0, "twopmodq78_q8 : (q3.d1 >> 14) != 0");
		ASSERT(HERE, (q4.d1 >> 14) == 0, "twopmodq78_q8 : (q4.d1 >> 14) != 0");
		ASSERT(HERE, (q5.d1 >> 14) == 0, "twopmodq78_q8 : (q5.d1 >> 14) != 0");
		ASSERT(HERE, (q6.d1 >> 14) == 0, "twopmodq78_q8 : (q6.d1 >> 14) != 0");
		ASSERT(HERE, (q7.d1 >> 14) == 0, "twopmodq78_q8 : (q7.d1 >> 14) != 0");

		*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0 + q4.d0 + q5.d0 + q6.d0 + q7.d0;

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
				"movq	%[__fqinv0],%%rax 	\n\t"\
				"movq	%[__two26i],%%rsi 	\n\t"\
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
			"movq	%[__fx0],%%rax 	\n\t"\
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
				"mulpd	%%xmm6,%%xmm0			\n\t	mulpd	%%xmm6 ,%%xmm8			\n\t				movq	%%rax,%%r8 	/* lo0 */	\n\t"\
				"																							movq	%%rdx,%%r12	/* hi0 */	\n\t"\
				"mulpd	%%xmm6,%%xmm2			\n\t	mulpd	%%xmm6 ,%%xmm10			\n\t				movq	0x310(%%rsi),%%rax	\n\t"\
				"mulpd	%%xmm6,%%xmm4			\n\t	mulpd	%%xmm6 ,%%xmm12			\n\t				mulq	%%rax		\n\t"\
				"																							movq	%%rax,%%r9 	/* lo1 */	\n\t"\
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
				"																							mulq	%%r8 	/* Use rax:rdx as product accumulator */\n\t"\
				"																							imulq	%%rdi,%%r12			/* qinv.lo * lo.hi */\n\t"\
				"subpd	%%xmm6,%%xmm3			\n\t	subpd	%%xmm14,%%xmm11			\n\t				addq	%%r12,%%rdx	\n\t"\
				"mulpd	%%xmm5,%%xmm3			\n\t	mulpd	%%xmm5 ,%%xmm11			\n\t				imulq	0x2c8(%%rsi),%%r8 	/* qinv.hi * lo.lo */\n\t"\
				"																							addq	%%r8 ,%%rdx	\n\t"\
				"mulpd	%%xmm15,%%xmm6			\n\t	mulpd	%%xmm15,%%xmm14			\n\t				movq	%%rax,%%r12	/* Move low 64 bits into free reg */\n\t"\
				"/* Digit 4: */					\n\t														movl	%%edx,0x348(%%rsi)	/* Write hi 32 bits into [__lo0] */\n\t"\
				"addpd	%%xmm6,%%xmm4			\n\t	addpd	%%xmm14,%%xmm12			\n\t			/* lo1 * qinv1: */	\n\t"\
				"movaps	%%xmm4,%%xmm6			\n\t	movaps	%%xmm12,%%xmm14			\n\t				movq	0x2d0(%%rsi),%%rdi	\n\t"\
				"roundpd	$0,%%xmm6,%%xmm6	\n\t	roundpd	$0,%%xmm14,%%xmm14		\n\t				movq	%%rdi,%%rax	\n\t"\
				"																							mulq	%%r9 		\n\t"\
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
				"movaps	%%xmm1,%%xmm3			\n\t	movaps	%%xmm9 ,%%xmm11			\n\t				movq	%%rdx,%%r8 	/* Discard low 64 bits [rax] */\n\t"\
				"																							movq	%%r12,%%rax	\n\t"\
				"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd	$0,%%xmm11,%%xmm11		\n\t				mulq	%%rdi	/* lo.lo*q.hi in rax:rdx */\n\t"\
				"																							xorq	%%r12,%%r12	\n\t"\
				"subpd	%%xmm3,%%xmm1			\n\t	subpd	%%xmm11,%%xmm9			\n\t				addq	%%rax,%%r8 	\n\t"\
				"																							adcq	%%rdx,%%r12	\n\t"\
				"mulpd	%%xmm5,%%xmm1			\n\t	mulpd	%%xmm5 ,%%xmm9			\n\t				movl	0x348(%%rsi),%%eax	/* Cannot do imulq __lohi,[reg] because that does a sign-extended load of __lohi */\n\t"\
				"mulpd	%%xmm15,%%xmm3			\n\t	mulpd	%%xmm15,%%xmm11			\n\t				imulq	%%rax,%%rdi			/* ...so first load __lohi into low 32 bits of a-reg, then compute lo.hi*q.hi. */\n\t"\
				"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t			addq	%%rdi,%%r12	\n\t"\
				"mulpd	0x80(%%rsi),%%xmm2		\n\t	mulpd	0x90(%%rsi),%%xmm10		\n\t				mulq	0x280(%%rsi)		/* q.lo*lo.hi in rax:rdx */\n\t"\
				"																							addq	%%rax,%%r8 	\n\t"\
				"mulpd	0xa0(%%rsi),%%xmm6		\n\t	mulpd	0xb0(%%rsi),%%xmm14		\n\t				adcq	%%rdx,%%r12	\n\t"\
				"mulpd	0xc0(%%rsi),%%xmm4		\n\t	mulpd	0xd0(%%rsi),%%xmm12		\n\t	/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
				"																							movq	0x300(%%rsi),%%rax	/* h.d0 */\n\t"\
				"addpd	%%xmm6,%%xmm2			\n\t	addpd	%%xmm14,%%xmm10			\n\t				movq	0x308(%%rsi),%%rdx	/* h.d1 */\n\t"\
				"addpd	%%xmm4,%%xmm3			\n\t	addpd	%%xmm12,%%xmm11			\n\t				subq	%%r8 ,%%rax	/* rax = (h-l).d0 */\n\t"\
				"addpd	%%xmm3,%%xmm2			\n\t	addpd	%%xmm11,%%xmm10			\n\t				sbbq	%%r12,%%rdx	/* rdx = (h-l).d1 */\n\t"\
				"																							sbbq	%%rdi,%%rdi	\n\t"\
				"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t				shrdq	$32,%%rdx,%%rax	/* Only keep high 96-bit output */\n\t"\
				"subpd	0x240(%%rsi),%%xmm3		\n\t	subpd	0x240(%%rsi),%%xmm11	\n\t				shrq	$32,%%rdx			\n\t"\
				"roundpd	$0,%%xmm3,%%xmm3	\n\t	roundpd	$0,%%xmm11,%%xmm11		\n\t				movq	0x280(%%rsi),%%r8 	\n\t"\
				"																							movq	0x288(%%rsi),%%r12	\n\t"\
				"																							andq	%%rdi,%%r8 	\n\t"\
				"subpd	%%xmm3,%%xmm2			\n\t	subpd	%%xmm11,%%xmm10			\n\t				andq	%%rdi,%%r12	\n\t"\
				"mulpd	%%xmm5,%%xmm2			\n\t	mulpd	%%xmm5 ,%%xmm10			\n\t				addq	%%rax,%%r8 	\n\t"\
				"																							adcq	%%rdx,%%r12	\n\t"\
				"/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */								\n\t			/* q1 * lo1: */	\n\t"\
				"/* Digit 0: */					\n\t														movq	0x290(%%rsi),%%rax	\n\t"\
				"movaps	%%xmm0,%%xmm3			\n\t	movaps	%%xmm8 ,%%xmm11			\n\t				mulq	%%r13	\n\t"\
				"movaps	%%xmm0,%%xmm4			\n\t	movaps	%%xmm8 ,%%xmm12			\n\t				movl	0x298(%%rsi),%%edi	\n\t"\
				"mulpd	    (%%rsi),%%xmm0		\n\t	mulpd	0x10(%%rsi),%%xmm8		\n\t				movq	%%rdx,%%r9 	\n\t"\
				"																							movq	%%r13,%%rax	\n\t"\
				"roundpd	$0,%%xmm0,%%xmm0	\n\t	roundpd	$0,%%xmm8 ,%%xmm8		\n\t				mulq	%%rdi	\n\t"\
				"																							xorq	%%r13,%%r13	\n\t"\
				"mulpd	    %%xmm15,%%xmm0		\n\t	mulpd	    %%xmm15,%%xmm8		\n\t				addq	%%rax,%%r9 	\n\t"\
				"/* Digit 1: */					\n\t														adcq	%%rdx,%%r13	\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t	movaps	%%xmm9 ,%%xmm14			\n\t				movl	0x358(%%rsi),%%eax	\n\t"\
				"mulpd	0x20(%%rsi),%%xmm3		\n\t	mulpd	0x30(%%rsi),%%xmm11		\n\t				imulq	%%rax,%%rdi	\n\t"\
				"																							addq	%%rdi,%%r13	\n\t"\
				"addpd	%%xmm0,%%xmm3			\n\t	addpd	%%xmm8 ,%%xmm11			\n\t				mulq	0x290(%%rsi)	\n\t"\
				"																							addq	%%rax,%%r9 	\n\t"\
				"mulpd	    (%%rsi),%%xmm1		\n\t	mulpd	0x10(%%rsi),%%xmm9		\n\t				adcq	%%rdx,%%r13	\n\t"\
				"addpd	%%xmm3,%%xmm1			\n\t	addpd	%%xmm11,%%xmm9			\n\t	/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
				"roundpd	$0,%%xmm1,%%xmm1	\n\t	roundpd	$0,%%xmm9 ,%%xmm9		\n\t				movq	0x310(%%rsi),%%rax	\n\t"\
				"																							movq	0x318(%%rsi),%%rdx	\n\t"\
				"																							subq	%%r9 ,%%rax		\n\t"\
				"mulpd	    %%xmm15,%%xmm1		\n\t	mulpd	    %%xmm15,%%xmm9		\n\t				sbbq	%%r13,%%rdx		\n\t"\
				"																							sbbq	%%rdi,%%rdi	/* bitmask for q */\n\t"\
				"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */	\n\t			shrdq	$32,%%rdx,%%rax		\n\t"\
				"movaps	%%xmm2,%%xmm3			\n\t	movaps	%%xmm10,%%xmm11			\n\t				shrq	$32,%%rdx			\n\t"\
				"movaps	%%xmm6,%%xmm0			\n\t	movaps	%%xmm14,%%xmm8			\n\t				movq	0x290(%%rsi),%%r9 	\n\t"\
				"mulpd	    (%%rsi),%%xmm2		\n\t	mulpd	0x10(%%rsi),%%xmm10		\n\t				movq	0x298(%%rsi),%%r13	\n\t"\
				"																							andq	%%rdi,%%r9 	\n\t"\
				"mulpd	0x20(%%rsi),%%xmm6		\n\t	mulpd	0x30(%%rsi),%%xmm14		\n\t				andq	%%rdi,%%r13	\n\t"\
				"mulpd	0x40(%%rsi),%%xmm4		\n\t	mulpd	0x50(%%rsi),%%xmm12		\n\t				addq	%%rax,%%r9 	\n\t"\
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
				"movslq	%[__pshift],%%rax					\n\t"\
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
				"																							addq	%%r8 ,%%r8 	\n\t"\
				"																							adcq	%%r12,%%r12	\n\t"\
				"	movaps		%%xmm8 ,%%xmm14	/* cpy of qhi */						\n\t				movq	0x280(%%rsi),%%rax	\n\t"\
				"																							movq	0x288(%%rsi),%%rdx	\n\t"\
				"																							subq	%%rax,%%r8 		\n\t"\
				"	addpd		%%xmm2,%%xmm2											\n\t				sbbq	%%rdx,%%r12		\n\t"\
				"																							sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
				"																							andq	%%rdi,%%rax	\n\t"\
				"	addpd		%%xmm10,%%xmm10	/* top 52 bits */						\n\t				andq	%%rdi,%%rdx	\n\t"\
				"																							addq	%%rax,%%r8 	\n\t"\
				"																							adcq	%%rdx,%%r12	\n\t"\
				"	addpd		%%xmm1,%%xmm1											\n\t			/* x1: */	\n\t"\
				"																							addq	%%r9 ,%%r9 	\n\t"\
				"																							adcq	%%r13,%%r13	\n\t"\
				"	addpd		%%xmm9 ,%%xmm9	/* low 26 bits */						\n\t				movq	0x290(%%rsi),%%rax	\n\t"\
				"																							movq	0x298(%%rsi),%%rdx	\n\t"\
				"/* If x > q, subtract q: */											\n\t				subq	%%rax,%%r9 		\n\t"\
				"	cmppd	$0x2,%%xmm2,%%xmm6											\n\t				sbbq	%%rdx,%%r13		\n\t"\
				"																							sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
				"																							andq	%%rdi,%%rax	\n\t"\
				"	cmppd	$0x2,%%xmm10,%%xmm14/* bitmask = (qhi <= xhi) */			\n\t				andq	%%rdi,%%rdx	\n\t"\
				"																							addq	%%rax,%%r9 	\n\t"\
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
			"movq	%[__fx0],%%rax 	\n\t"\
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
		lo0.d1 = tmp0 + qmul0*q0.d1;
		lo1.d1 = tmp1 + qmul1*q1.d1;
		lo2.d1 = tmp2 + qmul2*q2.d1;
		lo3.d1 = tmp3 + qmul3*q3.d1;
		SUB96(x0,lo0,x0);
		SUB96(x1,lo1,x1);
		SUB96(x2,lo2,x2);
		SUB96(x3,lo3,x3);
	} else {
		ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
		ADD96(x1,x1,x1);
		ADD96(x2,x2,x2);
		ADD96(x3,x3,x3);

		SUB96(x0,q0,x0);
		SUB96(x1,q1,x1);
		SUB96(x2,q2,x2);
		SUB96(x3,q3,x3);
	}

		ADD96_PTR(x4 ,x4, x4);
		ADD96_PTR(x5 ,x5, x5);
		ADD96_PTR(x6 ,x6, x6);
		ADD96_PTR(x7 ,x7, x7);

		SUB96_PTR(x4, qptr4, x4);
		SUB96_PTR(x5, qptr5, x5);
		SUB96_PTR(x6, qptr6, x6);
		SUB96_PTR(x7, qptr7, x7);

		tmp0 = CMPEQ96(x0, ONE96);			*checksum2 += x0.d0;
		tmp1 = CMPEQ96(x1, ONE96);			*checksum2 += x1.d0;
		tmp2 = CMPEQ96(x2, ONE96);			*checksum2 += x2.d0;
		tmp3 = CMPEQ96(x3, ONE96);			*checksum2 += x3.d0;
		tmp4 = CMPEQ96_PTR(x4, ONE96_PTR);	*checksum2 += x4->d0;
		tmp5 = CMPEQ96_PTR(x5, ONE96_PTR);	*checksum2 += x5->d0;
		tmp6 = CMPEQ96_PTR(x6, ONE96_PTR);	*checksum2 += x6->d0;
		tmp7 = CMPEQ96_PTR(x7, ONE96_PTR);	*checksum2 += x7->d0;

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

  #ifdef USE_AVX	// Nov 2013: Now that have prototyped 2-TF-input/4-SSE-register modmul macro, use 8-input all-float in both SSE2 and AVX mode

	// AVX/8-input analog of twopmodq78_3WORD_DOUBLE_q4 :
	uint64 twopmodq78_3WORD_DOUBLE_q16(uint64*checksum1, uint64*checksum2, uint64 p
		, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7
		, uint64 k8, uint64 k9, uint64 ka, uint64 kb, uint64 kc, uint64 kd, uint64 ke, uint64 kf
	)
	{
		 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
		uint32 q32_0, qinv32_0, tmp32_0, q32_8, qinv32_8, tmp32_8
			 , q32_1, qinv32_1, tmp32_1, q32_9, qinv32_9, tmp32_9
			 , q32_2, qinv32_2, tmp32_2, q32_a, qinv32_a, tmp32_a
			 , q32_3, qinv32_3, tmp32_3, q32_b, qinv32_b, tmp32_b
			 , q32_4, qinv32_4, tmp32_4, q32_c, qinv32_c, tmp32_c
			 , q32_5, qinv32_5, tmp32_5, q32_d, qinv32_d, tmp32_d
			 , q32_6, qinv32_6, tmp32_6, q32_e, qinv32_e, tmp32_e
			 , q32_7, qinv32_7, tmp32_7, q32_f, qinv32_f, tmp32_f;
		uint64 lead7, r = 0
			, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
			, tmp8, tmp9, tmpa, tmpb, tmpc, tmpd, tmpe, tmpf;
		uint96 q0, qinv0, qhalf0, x0, lo0, q8, qinv8, qhalf8, x8, lo8
			 , q1, qinv1, qhalf1, x1, lo1, q9, qinv9, qhalf9, x9, lo9
			 , q2, qinv2, qhalf2, x2, lo2, qa, qinva, qhalfa, xa, loa
			 , q3, qinv3, qhalf3, x3, lo3, qb, qinvb, qhalfb, xb, lob
			 , q4, qinv4, qhalf4, x4, lo4, qc, qinvc, qhalfc, xc, loc
			 , q5, qinv5, qhalf5, x5, lo5, qd, qinvd, qhalfd, xd, lod
			 , q6, qinv6, qhalf6, x6, lo6, qe, qinve, qhalfe, xe, loe
			 , q7, qinv7, qhalf7, x7, lo7, qf, qinvf, qhalff, xf, lof;
		uint192 prod192;
		static uint64 psave = 0;
		static uint32 pshift;	/* Require this to be 32-bit here for 32-bit ASM support */
		static uint32 start_index, zshift, first_entry = TRUE;
		static vec_dbl *sc_arr = 0x0;
		static double *sc_ptr;
		double *tmp;
		static double *aq0,*aq1,*aq2,*aqhi52, *aqinv0,*aqinv1,*aqinv2, *ax0,*ax1,*ax2;
		static double *bq0,*bq1,*bq2,*bqhi52, *bqinv0,*bqinv1,*bqinv2, *bx0,*bx1,*bx2;
		static double *cq0,*cq1,*cq2,*cqhi52, *cqinv0,*cqinv1,*cqinv2, *cx0,*cx1,*cx2;
		static double *dq0,*dq1,*dq2,*dqhi52, *dqinv0,*dqinv1,*dqinv2, *dx0,*dx1,*dx2;
		static double *eq0,*eq1,*eq2,*eqhi52, *eqinv0,*eqinv1,*eqinv2, *ex0,*ex1,*ex2;
		static double *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2;
		static double *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2;
		static double *hq0,*hq1,*hq2,*hqhi52, *hqinv0,*hqinv1,*hqinv2, *hx0,*hx1,*hx2;
		static double *iq0,*iq1,*iq2,*iqhi52, *iqinv0,*iqinv1,*iqinv2, *ix0,*ix1,*ix2;
		static double *jq0,*jq1,*jq2,*jqhi52, *jqinv0,*jqinv1,*jqinv2, *jx0,*jx1,*jx2;
		static double *kq0,*kq1,*kq2,*kqhi52, *kqinv0,*kqinv1,*kqinv2, *kx0,*kx1,*kx2;
		static double *lq0,*lq1,*lq2,*lqhi52, *lqinv0,*lqinv1,*lqinv2, *lx0,*lx1,*lx2;
		static double *mq0,*mq1,*mq2,*mqhi52, *mqinv0,*mqinv1,*mqinv2, *mx0,*mx1,*mx2;
		static double *nq0,*nq1,*nq2,*nqhi52, *nqinv0,*nqinv1,*nqinv2, *nx0,*nx1,*nx2;
		static double *oq0,*oq1,*oq2,*oqhi52, *oqinv0,*oqinv1,*oqinv2, *ox0,*ox1,*ox2;
		static double *pq0,*pq1,*pq2,*pqhi52, *pqinv0,*pqinv1,*pqinv2, *px0,*px1,*px2;
		static vec_dbl *two26f, *two26i, *two13i;

		if(first_entry)
		{
			sc_arr = ALLOC_VEC_DBL(sc_arr, 52);
			if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
			sc_ptr = (double *)ALIGN_VEC_DBL(sc_arr);
			ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
			/* Remember, rhese are POINTERS-TO-DOUBLES, so need an increment of 4 to span an AVX register: */
			tmp = (double *)sc_ptr;
		// 1st octet:
			aq0    = tmp + 0x00;	bq0    = tmp + 0x01;	cq0    = tmp + 0x02;	dq0    = tmp + 0x03;	eq0    = tmp + 0x04;	fq0    = tmp + 0x05;	gq0    = tmp + 0x06;	hq0    = tmp + 0x07;
			aq1    = tmp + 0x10;	bq1    = tmp + 0x11;	cq1    = tmp + 0x12;	dq1    = tmp + 0x13;	eq1    = tmp + 0x14;	fq1    = tmp + 0x15;	gq1    = tmp + 0x16;	hq1    = tmp + 0x17;
			aq2    = tmp + 0x20;	bq2    = tmp + 0x21;	cq2    = tmp + 0x22;	dq2    = tmp + 0x23;	eq2    = tmp + 0x24;	fq2    = tmp + 0x25;	gq2    = tmp + 0x26;	hq2    = tmp + 0x27;
			aqhi52 = tmp + 0x30;	bqhi52 = tmp + 0x31;	cqhi52 = tmp + 0x32;	dqhi52 = tmp + 0x33;	eqhi52 = tmp + 0x34;	fqhi52 = tmp + 0x35;	gqhi52 = tmp + 0x36;	hqhi52 = tmp + 0x37;
			aqinv0 = tmp + 0x40;	bqinv0 = tmp + 0x41;	cqinv0 = tmp + 0x42;	dqinv0 = tmp + 0x43;	eqinv0 = tmp + 0x44;	fqinv0 = tmp + 0x45;	gqinv0 = tmp + 0x46;	hqinv0 = tmp + 0x47;
			aqinv1 = tmp + 0x50;	bqinv1 = tmp + 0x51;	cqinv1 = tmp + 0x52;	dqinv1 = tmp + 0x53;	eqinv1 = tmp + 0x54;	fqinv1 = tmp + 0x55;	gqinv1 = tmp + 0x56;	hqinv1 = tmp + 0x57;
			aqinv2 = tmp + 0x60;	bqinv2 = tmp + 0x61;	cqinv2 = tmp + 0x62;	dqinv2 = tmp + 0x63;	eqinv2 = tmp + 0x64;	fqinv2 = tmp + 0x65;	gqinv2 = tmp + 0x66;	hqinv2 = tmp + 0x67;
			ax0    = tmp + 0x70;	bx0    = tmp + 0x71;	cx0    = tmp + 0x72;	dx0    = tmp + 0x73;	ex0    = tmp + 0x74;	fx0    = tmp + 0x75;	gx0    = tmp + 0x76;	hx0    = tmp + 0x77;
			ax1    = tmp + 0x80;	bx1    = tmp + 0x81;	cx1    = tmp + 0x82;	dx1    = tmp + 0x83;	ex1    = tmp + 0x84;	fx1    = tmp + 0x85;	gx1    = tmp + 0x86;	hx1    = tmp + 0x87;
			ax2    = tmp + 0x90;	bx2    = tmp + 0x91;	cx2    = tmp + 0x92;	dx2    = tmp + 0x93;	ex2    = tmp + 0x94;	fx2    = tmp + 0x95;	gx2    = tmp + 0x96;	hx2    = tmp + 0x97;
		// 2nd octet:
			iq0    = tmp + 0x08;	jq0    = tmp + 0x09;	kq0    = tmp + 0x0a;	lq0    = tmp + 0x0b;	mq0    = tmp + 0x0c;	nq0    = tmp + 0x0d;	oq0    = tmp + 0x0e;	pq0    = tmp + 0x0f;
			iq1    = tmp + 0x18;	jq1    = tmp + 0x19;	kq1    = tmp + 0x1a;	lq1    = tmp + 0x1b;	mq1    = tmp + 0x1c;	nq1    = tmp + 0x1d;	oq1    = tmp + 0x1e;	pq1    = tmp + 0x1f;
			iq2    = tmp + 0x28;	jq2    = tmp + 0x29;	kq2    = tmp + 0x2a;	lq2    = tmp + 0x2b;	mq2    = tmp + 0x2c;	nq2    = tmp + 0x2d;	oq2    = tmp + 0x2e;	pq2    = tmp + 0x2f;
			iqhi52 = tmp + 0x38;	jqhi52 = tmp + 0x39;	kqhi52 = tmp + 0x3a;	lqhi52 = tmp + 0x3b;	mqhi52 = tmp + 0x3c;	nqhi52 = tmp + 0x3d;	oqhi52 = tmp + 0x3e;	pqhi52 = tmp + 0x3f;
			iqinv0 = tmp + 0x48;	jqinv0 = tmp + 0x49;	kqinv0 = tmp + 0x4a;	lqinv0 = tmp + 0x4b;	mqinv0 = tmp + 0x4c;	nqinv0 = tmp + 0x4d;	oqinv0 = tmp + 0x4e;	pqinv0 = tmp + 0x4f;
			iqinv1 = tmp + 0x58;	jqinv1 = tmp + 0x59;	kqinv1 = tmp + 0x5a;	lqinv1 = tmp + 0x5b;	mqinv1 = tmp + 0x5c;	nqinv1 = tmp + 0x5d;	oqinv1 = tmp + 0x5e;	pqinv1 = tmp + 0x5f;
			iqinv2 = tmp + 0x68;	jqinv2 = tmp + 0x69;	kqinv2 = tmp + 0x6a;	lqinv2 = tmp + 0x6b;	mqinv2 = tmp + 0x6c;	nqinv2 = tmp + 0x6d;	oqinv2 = tmp + 0x6e;	pqinv2 = tmp + 0x6f;
			ix0    = tmp + 0x78;	jx0    = tmp + 0x79;	kx0    = tmp + 0x7a;	lx0    = tmp + 0x7b;	mx0    = tmp + 0x7c;	nx0    = tmp + 0x7d;	ox0    = tmp + 0x7e;	px0    = tmp + 0x7f;
			ix1    = tmp + 0x88;	jx1    = tmp + 0x89;	kx1    = tmp + 0x8a;	lx1    = tmp + 0x8b;	mx1    = tmp + 0x8c;	nx1    = tmp + 0x8d;	ox1    = tmp + 0x8e;	px1    = tmp + 0x8f;
			ix2    = tmp + 0x98;	jx2    = tmp + 0x99;	kx2    = tmp + 0x9a;	lx2    = tmp + 0x9b;	mx2    = tmp + 0x9c;	nx2    = tmp + 0x9d;	ox2    = tmp + 0x9e;	px2    = tmp + 0x9f;
			// +0xa0,b0 - Insert another 2 pairs of padding slots here for high-product-words register spills (we spill 2 of 3 words)
			two13i = (vec_dbl *)(sc_ptr + 0xc0);
			two26f = (vec_dbl *)(sc_ptr + 0xc4);
			two26i = (vec_dbl *)(sc_ptr + 0xc8);	// Equivalent of 12*4 + 3 = 51 AVX slots, alloc 52
			/* Can premultiply each of the multiword-mul inputs by 1/sqrt(2^13) due to quadraticity: */
			VEC_DBL_INIT(two13i, TWO13FLINV);
			VEC_DBL_INIT(two26f, TWO26FLOAT);
			VEC_DBL_INIT(two26i, TWO26FLINV);
		}	/* first_entry */

		ASSERT(HERE, (p >> 63) == 0, "twopmodq78_q4 : p must be < 2^63!");
		q0.d0 = q1.d0 = q2.d0 = q3.d0 = q4.d0 = q5.d0 = q6.d0 = q7.d0 = q8.d0 = q9.d0 = qa.d0 = qb.d0 = qc.d0 = qd.d0 = qe.d0 = qf.d0 = p+p;
	#ifdef MUL_LOHI64_SUBROUTINE
		#error MUL_LOHI64_SUBROUTINE defined!
	#endif
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);		MUL_LOHI64(q8.d0, k8, q8.d0, q8.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);		MUL_LOHI64(q9.d0, k9, q9.d0, q9.d1);
		MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);		MUL_LOHI64(qa.d0, ka, qa.d0, qa.d1);
		MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);		MUL_LOHI64(qb.d0, kb, qb.d0, qb.d1);
		MUL_LOHI64(q4.d0, k4, q4.d0, q4.d1);		MUL_LOHI64(qc.d0, kc, qc.d0, qc.d1);
		MUL_LOHI64(q5.d0, k5, q5.d0, q5.d1);		MUL_LOHI64(qd.d0, kd, qd.d0, qd.d1);
		MUL_LOHI64(q6.d0, k6, q6.d0, q6.d1);		MUL_LOHI64(qe.d0, ke, qe.d0, qe.d1);
		MUL_LOHI64(q7.d0, k7, q7.d0, q7.d1);		MUL_LOHI64(qf.d0, kf, qf.d0, qf.d1);
		/* Since 2*p*k even, no need to check for overflow here */
		q0.d0 += 1;		q8.d0 += 1;
		q1.d0 += 1;		q9.d0 += 1;
		q2.d0 += 1;		qa.d0 += 1;
		q3.d0 += 1;		qb.d0 += 1;
		q4.d0 += 1;		qc.d0 += 1;
		q5.d0 += 1;		qd.d0 += 1;
		q6.d0 += 1;		qe.d0 += 1;
		q7.d0 += 1;		qf.d0 += 1;
		ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q8 : (q0.d1 >> 14) != 0");

		*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0 + q4.d0 + q5.d0 + q6.d0 + q7.d0 + q8.d0 + q9.d0 + qa.d0 + qb.d0 + qc.d0 + qd.d0 + qe.d0 + qf.d0;

		q32_0 = (uint32)q0.d0;		q32_8 = (uint32)q8.d0;
		q32_1 = (uint32)q1.d0;		q32_9 = (uint32)q9.d0;
		q32_2 = (uint32)q2.d0;		q32_a = (uint32)qa.d0;
		q32_3 = (uint32)q3.d0;		q32_b = (uint32)qb.d0;
		q32_4 = (uint32)q4.d0;		q32_c = (uint32)qc.d0;
		q32_5 = (uint32)q5.d0;		q32_d = (uint32)qd.d0;
		q32_6 = (uint32)q6.d0;		q32_e = (uint32)qe.d0;
		q32_7 = (uint32)q7.d0;		q32_f = (uint32)qf.d0;

		/* Convert q to floating form: */
		CVT_UINT78_3WORD_DOUBLE(q0 ,*aq0,*aq1,*aq2);
		CVT_UINT78_3WORD_DOUBLE(q1 ,*bq0,*bq1,*bq2);
		CVT_UINT78_3WORD_DOUBLE(q2 ,*cq0,*cq1,*cq2);
		CVT_UINT78_3WORD_DOUBLE(q3 ,*dq0,*dq1,*dq2);
		CVT_UINT78_3WORD_DOUBLE(q4 ,*eq0,*eq1,*eq2);
		CVT_UINT78_3WORD_DOUBLE(q5 ,*fq0,*fq1,*fq2);
		CVT_UINT78_3WORD_DOUBLE(q6 ,*gq0,*gq1,*gq2);
		CVT_UINT78_3WORD_DOUBLE(q7 ,*hq0,*hq1,*hq2);
		CVT_UINT78_3WORD_DOUBLE(q8 ,*iq0,*iq1,*iq2);
		CVT_UINT78_3WORD_DOUBLE(q9 ,*jq0,*jq1,*jq2);
		CVT_UINT78_3WORD_DOUBLE(qa ,*kq0,*kq1,*kq2);
		CVT_UINT78_3WORD_DOUBLE(qb ,*lq0,*lq1,*lq2);
		CVT_UINT78_3WORD_DOUBLE(qc ,*mq0,*mq1,*mq2);
		CVT_UINT78_3WORD_DOUBLE(qd ,*nq0,*nq1,*nq2);
		CVT_UINT78_3WORD_DOUBLE(qe ,*oq0,*oq1,*oq2);
		CVT_UINT78_3WORD_DOUBLE(qf ,*pq0,*pq1,*pq2);

		__asm__ volatile (\
			"movq	%[__aq0],%%rax						\n\t"\
			"movq	%[__two26f],%%rsi					\n\t"\
			"vmovaps	    (%%rsi),%%ymm8	/* two26f */	\n\t"\
			"vmovaps	0x20(%%rsi),%%ymm9	/* two26i */	\n\t"\
			"vmovaps	     (%%rax),%%ymm0	/* aq0 */		\n\t	vmovaps	0x020(%%rax),%%ymm4	/* eq0 */		\n\t"\
			"vmovaps	0x080(%%rax),%%ymm1	/* aq1 */		\n\t	vmovaps	0x0a0(%%rax),%%ymm5	/* eq1 */		\n\t"\
			"vmovaps	0x100(%%rax),%%ymm2	/* aq2 */		\n\t	vmovaps	0x120(%%rax),%%ymm6	/* eq2 */		\n\t"\
			"vmovaps	%%ymm2,%%ymm3		/* cpy aq2 */	\n\t	vmovaps	%%ymm6,%%ymm7		/* cpy eq2 */	\n\t"\
			"vmulpd	%%ymm8,%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm8,%%ymm7,%%ymm7				\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3/* Hi 52 out bits */\n\t	vaddpd	%%ymm5,%%ymm7,%%ymm7	/* Hi 52*/	\n\t"\
			"vmulpd	%%ymm9,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm9,%%ymm4,%%ymm4				\n\t"\
			"vmulpd	%%ymm9,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm9,%%ymm5,%%ymm5				\n\t"\
			"vmulpd	%%ymm9,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm9,%%ymm6,%%ymm6				\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	/* aq0/2^26 */	\n\t	vmovaps	%%ymm4,0x020(%%rax)	/* eq0/2^26 */	\n\t"\
			"vmovaps	%%ymm1,0x080(%%rax)	/* aq1/2^26 */	\n\t	vmovaps	%%ymm5,0x0a0(%%rax)	/* eq1/2^26 */	\n\t"\
			"vmovaps	%%ymm2,0x100(%%rax)	/* aq2/2^26 */	\n\t	vmovaps	%%ymm6,0x120(%%rax)	/* eq2/2^26 */	\n\t"\
			"vmovaps	%%ymm3,0x180(%%rax)	/* aq[hi52] */	\n\t	vmovaps	%%ymm7,0x1a0(%%rax)	/* eq[hi52] */	\n\t"\
			/* Now do 2nd set of doubles: */\
			"vmovaps	0x040(%%rax),%%ymm0	/* iq0 */		\n\t	vmovaps	0x060(%%rax),%%ymm4	/* mq0 */		\n\t"\
			"vmovaps	0x0c0(%%rax),%%ymm1	/* iq1 */		\n\t	vmovaps	0x0e0(%%rax),%%ymm5	/* mq1 */		\n\t"\
			"vmovaps	0x140(%%rax),%%ymm2	/* iq2 */		\n\t	vmovaps	0x160(%%rax),%%ymm6	/* mq2 */		\n\t"\
			"vmovaps	%%ymm2,%%ymm3		/* cpy iq2 */	\n\t	vmovaps	%%ymm6,%%ymm7		/* cpy mq2 */	\n\t"\
			"vmulpd	%%ymm8,%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm8,%%ymm7,%%ymm7				\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3/* Hi 52 out bits */\n\t	vaddpd	%%ymm5,%%ymm7,%%ymm7	/* Hi 52*/	\n\t"\
			"vmulpd	%%ymm9,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm9,%%ymm4,%%ymm4				\n\t"\
			"vmulpd	%%ymm9,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm9,%%ymm5,%%ymm5				\n\t"\
			"vmulpd	%%ymm9,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm9,%%ymm6,%%ymm6				\n\t"\
			"vmovaps	%%ymm0,0x040(%%rax)	/* iq0/2^26 */	\n\t	vmovaps	%%ymm4,0x060(%%rax)	/* mq0/2^26 */	\n\t"\
			"vmovaps	%%ymm1,0x0c0(%%rax)	/* iq1/2^26 */	\n\t	vmovaps	%%ymm5,0x0e0(%%rax)	/* mq1/2^26 */	\n\t"\
			"vmovaps	%%ymm2,0x140(%%rax)	/* iq2/2^26 */	\n\t	vmovaps	%%ymm6,0x160(%%rax)	/* mq2/2^26 */	\n\t"\
			"vmovaps	%%ymm3,0x1c0(%%rax)	/* iq[hi52] */	\n\t	vmovaps	%%ymm7,0x1e0(%%rax)	/* mq[hi52] */	\n\t"\
			:					/* outputs: none */\
			: [__aq0] "m" (aq0)	/* All inputs from memory addresses here */\
			 ,[__two26f] "m" (two26f)\
			: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
		);

		/* = (q-1)/2, since q odd. */
		RSHIFT_FAST96(q0, 1, qhalf0);		RSHIFT_FAST96(q8, 1, qhalf8);
		RSHIFT_FAST96(q1, 1, qhalf1);		RSHIFT_FAST96(q9, 1, qhalf9);
		RSHIFT_FAST96(q2, 1, qhalf2);		RSHIFT_FAST96(qa, 1, qhalfa);
		RSHIFT_FAST96(q3, 1, qhalf3);		RSHIFT_FAST96(qb, 1, qhalfb);
		RSHIFT_FAST96(q4, 1, qhalf4);		RSHIFT_FAST96(qc, 1, qhalfc);
		RSHIFT_FAST96(q5, 1, qhalf5);		RSHIFT_FAST96(qd, 1, qhalfd);
		RSHIFT_FAST96(q6, 1, qhalf6);		RSHIFT_FAST96(qe, 1, qhalfe);
		RSHIFT_FAST96(q7, 1, qhalf7);		RSHIFT_FAST96(qf, 1, qhalff);

		if(first_entry || p != psave)
		{
			first_entry = FALSE;
			psave  = p;
			ASSERT(HERE, ((p + 78) >> 32) == 0, "twopmodq78_q2 : (p+78) exceeds 32 bits!");
			pshift = (uint32)p + 78;
			j = leadz32(pshift);
			lead7 = ((pshift<<j) >> 25);
			if(lead7 > 77)
			{
				lead7 >>= 1;
				start_index =  32-j-6;	/* Use only the leftmost 6 bits */
			}
			else
				start_index =  32-j-7;

			zshift = 77 - lead7;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

			pshift = ~pshift;
		}

		qinv32_0 = (q32_0 + q32_0 + q32_0) ^ (uint32)2;	qinv32_8 = (q32_8 + q32_8 + q32_8) ^ (uint32)2;
		qinv32_1 = (q32_1 + q32_1 + q32_1) ^ (uint32)2;	qinv32_9 = (q32_9 + q32_9 + q32_9) ^ (uint32)2;
		qinv32_2 = (q32_2 + q32_2 + q32_2) ^ (uint32)2;	qinv32_a = (q32_a + q32_a + q32_a) ^ (uint32)2;
		qinv32_3 = (q32_3 + q32_3 + q32_3) ^ (uint32)2;	qinv32_b = (q32_b + q32_b + q32_b) ^ (uint32)2;
		qinv32_4 = (q32_4 + q32_4 + q32_4) ^ (uint32)2;	qinv32_c = (q32_c + q32_c + q32_c) ^ (uint32)2;
		qinv32_5 = (q32_5 + q32_5 + q32_5) ^ (uint32)2;	qinv32_d = (q32_d + q32_d + q32_d) ^ (uint32)2;
		qinv32_6 = (q32_6 + q32_6 + q32_6) ^ (uint32)2;	qinv32_e = (q32_e + q32_e + q32_e) ^ (uint32)2;
		qinv32_7 = (q32_7 + q32_7 + q32_7) ^ (uint32)2;	qinv32_f = (q32_f + q32_f + q32_f) ^ (uint32)2;
		/* 4 bits: */
		tmp32_0 = q32_0*qinv32_0;						tmp32_8 = q32_8*qinv32_8;
		tmp32_1 = q32_1*qinv32_1;						tmp32_9 = q32_9*qinv32_9;
		tmp32_2 = q32_2*qinv32_2;						tmp32_a = q32_a*qinv32_a;
		tmp32_3 = q32_3*qinv32_3;						tmp32_b = q32_b*qinv32_b;
		tmp32_4 = q32_4*qinv32_4;						tmp32_c = q32_c*qinv32_c;
		tmp32_5 = q32_5*qinv32_5;						tmp32_d = q32_d*qinv32_d;
		tmp32_6 = q32_6*qinv32_6;						tmp32_e = q32_e*qinv32_e;
		tmp32_7 = q32_7*qinv32_7;						tmp32_f = q32_f*qinv32_f;
		qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);		qinv32_8 = qinv32_8*((uint32)2 - tmp32_8);
		qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);		qinv32_9 = qinv32_9*((uint32)2 - tmp32_9);
		qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);		qinv32_a = qinv32_a*((uint32)2 - tmp32_a);
		qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);		qinv32_b = qinv32_b*((uint32)2 - tmp32_b);
		qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);		qinv32_c = qinv32_c*((uint32)2 - tmp32_c);
		qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);		qinv32_d = qinv32_d*((uint32)2 - tmp32_d);
		qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);		qinv32_e = qinv32_e*((uint32)2 - tmp32_e);
		qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);		qinv32_f = qinv32_f*((uint32)2 - tmp32_f);
		/* 8 bits: */
		tmp32_0 = q32_0*qinv32_0;						tmp32_8 = q32_8*qinv32_8;
		tmp32_1 = q32_1*qinv32_1;						tmp32_9 = q32_9*qinv32_9;
		tmp32_2 = q32_2*qinv32_2;						tmp32_a = q32_a*qinv32_a;
		tmp32_3 = q32_3*qinv32_3;						tmp32_b = q32_b*qinv32_b;
		tmp32_4 = q32_4*qinv32_4;						tmp32_c = q32_c*qinv32_c;
		tmp32_5 = q32_5*qinv32_5;						tmp32_d = q32_d*qinv32_d;
		tmp32_6 = q32_6*qinv32_6;						tmp32_e = q32_e*qinv32_e;
		tmp32_7 = q32_7*qinv32_7;						tmp32_f = q32_f*qinv32_f;
		qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);		qinv32_8 = qinv32_8*((uint32)2 - tmp32_8);
		qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);		qinv32_9 = qinv32_9*((uint32)2 - tmp32_9);
		qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);		qinv32_a = qinv32_a*((uint32)2 - tmp32_a);
		qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);		qinv32_b = qinv32_b*((uint32)2 - tmp32_b);
		qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);		qinv32_c = qinv32_c*((uint32)2 - tmp32_c);
		qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);		qinv32_d = qinv32_d*((uint32)2 - tmp32_d);
		qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);		qinv32_e = qinv32_e*((uint32)2 - tmp32_e);
		qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);		qinv32_f = qinv32_f*((uint32)2 - tmp32_f);
		/* 16 bits: */
		tmp32_0 = q32_0*qinv32_0;						tmp32_8 = q32_8*qinv32_8;
		tmp32_1 = q32_1*qinv32_1;						tmp32_9 = q32_9*qinv32_9;
		tmp32_2 = q32_2*qinv32_2;						tmp32_a = q32_a*qinv32_a;
		tmp32_3 = q32_3*qinv32_3;						tmp32_b = q32_b*qinv32_b;
		tmp32_4 = q32_4*qinv32_4;						tmp32_c = q32_c*qinv32_c;
		tmp32_5 = q32_5*qinv32_5;						tmp32_d = q32_d*qinv32_d;
		tmp32_6 = q32_6*qinv32_6;						tmp32_e = q32_e*qinv32_e;
		tmp32_7 = q32_7*qinv32_7;						tmp32_f = q32_f*qinv32_f;
		qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);		qinv32_8 = qinv32_8*((uint32)2 - tmp32_8);
		qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);		qinv32_9 = qinv32_9*((uint32)2 - tmp32_9);
		qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);		qinv32_a = qinv32_a*((uint32)2 - tmp32_a);
		qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);		qinv32_b = qinv32_b*((uint32)2 - tmp32_b);
		qinv32_4 = qinv32_4*((uint32)2 - tmp32_4);		qinv32_c = qinv32_c*((uint32)2 - tmp32_c);
		qinv32_5 = qinv32_5*((uint32)2 - tmp32_5);		qinv32_d = qinv32_d*((uint32)2 - tmp32_d);
		qinv32_6 = qinv32_6*((uint32)2 - tmp32_6);		qinv32_e = qinv32_e*((uint32)2 - tmp32_e);
		qinv32_7 = qinv32_7*((uint32)2 - tmp32_7);		qinv32_f = qinv32_f*((uint32)2 - tmp32_f);
		/* 32 bits: */
		qinv0.d0 = (uint64)qinv32_0;					qinv8.d0 = (uint64)qinv32_8;
		qinv1.d0 = (uint64)qinv32_1;					qinv9.d0 = (uint64)qinv32_9;
		qinv2.d0 = (uint64)qinv32_2;					qinva.d0 = (uint64)qinv32_a;
		qinv3.d0 = (uint64)qinv32_3;					qinvb.d0 = (uint64)qinv32_b;
		qinv4.d0 = (uint64)qinv32_4;					qinvc.d0 = (uint64)qinv32_c;
		qinv5.d0 = (uint64)qinv32_5;					qinvd.d0 = (uint64)qinv32_d;
		qinv6.d0 = (uint64)qinv32_6;					qinve.d0 = (uint64)qinv32_e;
		qinv7.d0 = (uint64)qinv32_7;					qinvf.d0 = (uint64)qinv32_f;
		tmp0 = q0.d0*qinv0.d0;							tmp8 = q8.d0*qinv8.d0;
		tmp1 = q1.d0*qinv1.d0;							tmp9 = q9.d0*qinv9.d0;
		tmp2 = q2.d0*qinv2.d0;							tmpa = qa.d0*qinva.d0;
		tmp3 = q3.d0*qinv3.d0;							tmpb = qb.d0*qinvb.d0;
		tmp4 = q4.d0*qinv4.d0;							tmpc = qc.d0*qinvc.d0;
		tmp5 = q5.d0*qinv5.d0;							tmpd = qd.d0*qinvd.d0;
		tmp6 = q6.d0*qinv6.d0;							tmpe = qe.d0*qinve.d0;
		tmp7 = q7.d0*qinv7.d0;							tmpf = qf.d0*qinvf.d0;
		qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);			qinv8.d0 = qinv8.d0*((uint64)2 - tmp8);
		qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);			qinv9.d0 = qinv9.d0*((uint64)2 - tmp9);
		qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);			qinva.d0 = qinva.d0*((uint64)2 - tmpa);
		qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);			qinvb.d0 = qinvb.d0*((uint64)2 - tmpb);
		qinv4.d0 = qinv4.d0*((uint64)2 - tmp4);			qinvc.d0 = qinvc.d0*((uint64)2 - tmpc);
		qinv5.d0 = qinv5.d0*((uint64)2 - tmp5);			qinvd.d0 = qinvd.d0*((uint64)2 - tmpd);
		qinv6.d0 = qinv6.d0*((uint64)2 - tmp6);			qinve.d0 = qinve.d0*((uint64)2 - tmpe);
		qinv7.d0 = qinv7.d0*((uint64)2 - tmp7);			qinvf.d0 = qinvf.d0*((uint64)2 - tmpf);
		/* 64 bits: */
		MULH64(q0.d0, qinv0.d0, tmp0);					MULH64(q8.d0, qinv8.d0, tmp8);
		MULH64(q1.d0, qinv1.d0, tmp1);					MULH64(q9.d0, qinv9.d0, tmp9);
		MULH64(q2.d0, qinv2.d0, tmp2);					MULH64(qa.d0, qinva.d0, tmpa);
		MULH64(q3.d0, qinv3.d0, tmp3);					MULH64(qb.d0, qinvb.d0, tmpb);
		MULH64(q4.d0, qinv4.d0, tmp4);					MULH64(qc.d0, qinvc.d0, tmpc);
		MULH64(q5.d0, qinv5.d0, tmp5);					MULH64(qd.d0, qinvd.d0, tmpd);
		MULH64(q6.d0, qinv6.d0, tmp6);					MULH64(qe.d0, qinve.d0, tmpe);
		MULH64(q7.d0, qinv7.d0, tmp7);					MULH64(qf.d0, qinvf.d0, tmpf);
		qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);	qinv8.d1 = -qinv8.d0*(q8.d1*qinv8.d0 + tmp8);
		qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);	qinv9.d1 = -qinv9.d0*(q9.d1*qinv9.d0 + tmp9);
		qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + tmp2);	qinva.d1 = -qinva.d0*(qa.d1*qinva.d0 + tmpa);
		qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + tmp3);	qinvb.d1 = -qinvb.d0*(qb.d1*qinvb.d0 + tmpb);
		qinv4.d1 = -qinv4.d0*(q4.d1*qinv4.d0 + tmp4);	qinvc.d1 = -qinvc.d0*(qc.d1*qinvc.d0 + tmpc);
		qinv5.d1 = -qinv5.d0*(q5.d1*qinv5.d0 + tmp5);	qinvd.d1 = -qinvd.d0*(qd.d1*qinvd.d0 + tmpd);
		qinv6.d1 = -qinv6.d0*(q6.d1*qinv6.d0 + tmp6);	qinve.d1 = -qinve.d0*(qe.d1*qinve.d0 + tmpe);
		qinv7.d1 = -qinv7.d0*(q7.d1*qinv7.d0 + tmp7);	qinvf.d1 = -qinvf.d0*(qf.d1*qinvf.d0 + tmpf);
		qinv0.d1 &= 0x3fff;								qinv8.d1 &= 0x3fff;
		qinv1.d1 &= 0x3fff;								qinv9.d1 &= 0x3fff;
		qinv2.d1 &= 0x3fff;								qinva.d1 &= 0x3fff;
		qinv3.d1 &= 0x3fff;								qinvb.d1 &= 0x3fff;
		qinv4.d1 &= 0x3fff;								qinvc.d1 &= 0x3fff;
		qinv5.d1 &= 0x3fff;								qinvd.d1 &= 0x3fff;
		qinv6.d1 &= 0x3fff;								qinve.d1 &= 0x3fff;
		qinv7.d1 &= 0x3fff;								qinvf.d1 &= 0x3fff;
		/* Convert qinv to floating form: */
		CVT_UINT78_3WORD_DOUBLE(qinv0 ,*aqinv0,*aqinv1,*aqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv1 ,*bqinv0,*bqinv1,*bqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv2 ,*cqinv0,*cqinv1,*cqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv3 ,*dqinv0,*dqinv1,*dqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv4 ,*eqinv0,*eqinv1,*eqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv5 ,*fqinv0,*fqinv1,*fqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv6 ,*gqinv0,*gqinv1,*gqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv7 ,*hqinv0,*hqinv1,*hqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv8 ,*iqinv0,*iqinv1,*iqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinv9 ,*jqinv0,*jqinv1,*jqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinva ,*kqinv0,*kqinv1,*kqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinvb ,*lqinv0,*lqinv1,*lqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinvc ,*mqinv0,*mqinv1,*mqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinvd ,*nqinv0,*nqinv1,*nqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinve ,*oqinv0,*oqinv1,*oqinv2);
		CVT_UINT78_3WORD_DOUBLE(qinvf ,*pqinv0,*pqinv1,*pqinv2);

		__asm__ volatile (\
			"movq	%[__aqinv0],%%rax 	\n\t"\
			"movq	%[__two26i],%%rsi 	\n\t"\
			"vmovaps	(%%rsi),%%ymm12		\n\t"\
			"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	0x020(%%rax),%%ymm3	\n\t	vmovaps	0x040(%%rax),%%ymm6	\n\t	vmovaps	0x060(%%rax),%%ymm9	\n\t"\
			"vmovaps	0x080(%%rax),%%ymm1	\n\t	vmovaps	0x0a0(%%rax),%%ymm4	\n\t	vmovaps	0x0c0(%%rax),%%ymm7	\n\t	vmovaps	0x0e0(%%rax),%%ymm10	\n\t"\
			"vmovaps	0x100(%%rax),%%ymm2	\n\t	vmovaps	0x120(%%rax),%%ymm5	\n\t	vmovaps	0x140(%%rax),%%ymm8	\n\t	vmovaps	0x160(%%rax),%%ymm11	\n\t"\
			"vmulpd	%%ymm12,%%ymm0,%%ymm0	\n\t vmulpd	%%ymm12,%%ymm3,%%ymm3	\n\t vmulpd	%%ymm12,%%ymm6,%%ymm6	\n\t vmulpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
			"vmulpd	%%ymm12,%%ymm1,%%ymm1	\n\t vmulpd	%%ymm12,%%ymm4,%%ymm4	\n\t vmulpd	%%ymm12,%%ymm7,%%ymm7	\n\t vmulpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
			"vmulpd	%%ymm12,%%ymm2,%%ymm2	\n\t vmulpd	%%ymm12,%%ymm5,%%ymm5	\n\t vmulpd	%%ymm12,%%ymm8,%%ymm8	\n\t vmulpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm3,0x020(%%rax)	\n\t	vmovaps	%%ymm6,0x040(%%rax)	\n\t	vmovaps	%%ymm9 ,0x060(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x080(%%rax)	\n\t	vmovaps	%%ymm4,0x0a0(%%rax)	\n\t	vmovaps	%%ymm7,0x0c0(%%rax)	\n\t	vmovaps	%%ymm10,0x0e0(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x100(%%rax)	\n\t	vmovaps	%%ymm5,0x120(%%rax)	\n\t	vmovaps	%%ymm8,0x140(%%rax)	\n\t	vmovaps	%%ymm11,0x160(%%rax)	\n\t"\
			:					/* outputs: none */\
			: [__aqinv0] "m" (aqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			: "cc","memory","rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12"	/* Clobbered registers */\
		);

		/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
		j = start_index-1;

		/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
		LSHIFT96(qinv0, zshift, lo0);		LSHIFT96(qinv8, zshift, lo8);
		LSHIFT96(qinv1, zshift, lo1);		LSHIFT96(qinv9, zshift, lo9);
		LSHIFT96(qinv2, zshift, lo2);		LSHIFT96(qinva, zshift, loa);
		LSHIFT96(qinv3, zshift, lo3);		LSHIFT96(qinvb, zshift, lob);
		LSHIFT96(qinv4, zshift, lo4);		LSHIFT96(qinvc, zshift, loc);
		LSHIFT96(qinv5, zshift, lo5);		LSHIFT96(qinvd, zshift, lod);
		LSHIFT96(qinv6, zshift, lo6);		LSHIFT96(qinve, zshift, loe);
		LSHIFT96(qinv7, zshift, lo7);		LSHIFT96(qinvf, zshift, lof);
		lo0.d1 &= 0x3fff;					lo8.d1 &= 0x3fff;
		lo1.d1 &= 0x3fff;					lo9.d1 &= 0x3fff;
		lo2.d1 &= 0x3fff;					loa.d1 &= 0x3fff;
		lo3.d1 &= 0x3fff;					lob.d1 &= 0x3fff;
		lo4.d1 &= 0x3fff;					loc.d1 &= 0x3fff;
		lo5.d1 &= 0x3fff;					lod.d1 &= 0x3fff;
		lo6.d1 &= 0x3fff;					loe.d1 &= 0x3fff;
		lo7.d1 &= 0x3fff;					lof.d1 &= 0x3fff;
		// These share the prod192 variable, so cannot be done in side-by-side fashion:
		// Inline a (... >> 78) specilaization of the RSHIFT192 macro here, since have a constant shift count:
		MUL_LOHI96_PROD192(q0,lo0,prod192);	lo0.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo0.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q1,lo1,prod192);	lo1.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo1.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q2,lo2,prod192);	lo2.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo2.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q3,lo3,prod192);	lo3.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo3.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q4,lo4,prod192);	lo4.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo4.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q5,lo5,prod192);	lo5.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo5.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q6,lo6,prod192);	lo6.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo6.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q7,lo7,prod192);	lo7.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo7.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q8,lo8,prod192);	lo8.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo8.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(q9,lo9,prod192);	lo9.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lo9.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(qa,loa,prod192);	loa.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	loa.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(qb,lob,prod192);	lob.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lob.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(qc,loc,prod192);	loc.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	loc.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(qd,lod,prod192);	lod.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lod.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(qe,loe,prod192);	loe.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	loe.d1 = (prod192.d2 >> 14);
		MUL_LOHI96_PROD192(qf,lof,prod192);	lof.d0 = (prod192.d1 >> 14) + (prod192.d2 << 50);	lof.d1 = (prod192.d2 >> 14);

		/* hi = 0 in this instance, which simplifies things. */
		/* Put the result in lo (rather than x), to ease overflow check below */
		SUB96(q0, lo0, x0);		SUB96(q8, lo8, x8);
		SUB96(q1, lo1, x1);		SUB96(q9, lo9, x9);
		SUB96(q2, lo2, x2);		SUB96(qa, loa, xa);
		SUB96(q3, lo3, x3);		SUB96(qb, lob, xb);
		SUB96(q4, lo4, x4);		SUB96(qc, loc, xc);
		SUB96(q5, lo5, x5);		SUB96(qd, lod, xd);
		SUB96(q6, lo6, x6);		SUB96(qe, loe, xe);
		SUB96(q7, lo7, x7);		SUB96(qf, lof, xf);

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
			if(CMPUGT96(x8, qhalf8)){ ADD96(x8, x8, x8); SUB96(x8, q8, x8); }else{ ADD96(x8, x8, x8); }
			if(CMPUGT96(x9, qhalf9)){ ADD96(x9, x9, x9); SUB96(x9, q9, x9); }else{ ADD96(x9, x9, x9); }
			if(CMPUGT96(xa, qhalfa)){ ADD96(xa, xa, xa); SUB96(xa, qa, xa); }else{ ADD96(xa, xa, xa); }
			if(CMPUGT96(xb, qhalfb)){ ADD96(xb, xb, xb); SUB96(xb, qb, xb); }else{ ADD96(xb, xb, xb); }
			if(CMPUGT96(xc, qhalfc)){ ADD96(xc, xc, xc); SUB96(xc, qc, xc); }else{ ADD96(xc, xc, xc); }
			if(CMPUGT96(xd, qhalfd)){ ADD96(xd, xd, xd); SUB96(xd, qd, xd); }else{ ADD96(xd, xd, xd); }
			if(CMPUGT96(xe, qhalfe)){ ADD96(xe, xe, xe); SUB96(xe, qe, xe); }else{ ADD96(xe, xe, xe); }
			if(CMPUGT96(xf, qhalff)){ ADD96(xf, xf, xf); SUB96(xf, qf, xf); }else{ ADD96(xf, xf, xf); }
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
		CVT_UINT78_3WORD_DOUBLE(x8 ,*ix0,*ix1,*ix2);
		CVT_UINT78_3WORD_DOUBLE(x9 ,*jx0,*jx1,*jx2);
		CVT_UINT78_3WORD_DOUBLE(xa ,*kx0,*kx1,*kx2);
		CVT_UINT78_3WORD_DOUBLE(xb ,*lx0,*lx1,*lx2);
		CVT_UINT78_3WORD_DOUBLE(xc ,*mx0,*mx1,*mx2);
		CVT_UINT78_3WORD_DOUBLE(xd ,*nx0,*nx1,*nx2);
		CVT_UINT78_3WORD_DOUBLE(xe ,*ox0,*ox1,*ox2);
		CVT_UINT78_3WORD_DOUBLE(xf ,*px0,*px1,*px2);

		__asm__ volatile (\
			"movq	%[__ax0],%%rax 	\n\t"\
			"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	0x020(%%rax),%%ymm4	\n\t	vmovaps	0x040(%%rax),%%ymm8		\n\t	vmovaps	0x060(%%rax),%%ymm12	\n\t"\
			"vmovaps	0x080(%%rax),%%ymm1	\n\t	vmovaps	0x0a0(%%rax),%%ymm5	\n\t	vmovaps	0x0c0(%%rax),%%ymm9		\n\t	vmovaps	0x0e0(%%rax),%%ymm13	\n\t"\
			"vmovaps	0x100(%%rax),%%ymm2	\n\t	vmovaps	0x120(%%rax),%%ymm6	\n\t	vmovaps	0x140(%%rax),%%ymm10	\n\t	vmovaps	0x160(%%rax),%%ymm14	\n\t"\
			:					/* outputs: none */\
			: [__ax0] "m" (ax0)	/* All inputs from memory addresses here */\
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm12","xmm13","xmm14"	/* Clobbered registers */\
		);

		/*...x^2 mod q is returned in x. */
		/* All 3-word-double-form operands have components in the following size ranges:
			fword0,1 in [-2^25, +2^25]
			fword2   in [   -1, +2^26]
		*/

	/* Inner loop body needs 42 MOVAPS, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */
	/* !ASM_LOOP: High-level loop construct, and explicit branch to do the modular doubling: */

		for(j = start_index-2; j >= 0; j--)
		{
			SSE2_twopmodq78_modmul_q16(aq0,aqinv0,ax0,two26i,pshift,j);
		}	/* for(j...) */

		__asm__ volatile (\
			"movq	%[__ax0],%%rax 	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm4,0x020(%%rax)	\n\t	vmovaps	%%ymm8 ,0x040(%%rax)	\n\t	vmovaps	%%ymm12,0x060(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x080(%%rax)	\n\t	vmovaps	%%ymm5,0x0a0(%%rax)	\n\t	vmovaps	%%ymm9 ,0x0c0(%%rax)	\n\t	vmovaps	%%ymm13,0x0e0(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x100(%%rax)	\n\t	vmovaps	%%ymm6,0x120(%%rax)	\n\t	vmovaps	%%ymm10,0x140(%%rax)	\n\t	vmovaps	%%ymm14,0x160(%%rax)	\n\t"\
			:					/* outputs: none */\
			: [__ax0] "m" (ax0)	/* All inputs from memory addresses here */\
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm4","xmm5","xmm6","xmm8","xmm9","xmm10","xmm12","xmm13","xmm14"	/* Clobbered registers */\
		);

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
		CVT78_3WORD_DOUBLE_UINT96(*ix0,*ix1,*ix2, x8);
		CVT78_3WORD_DOUBLE_UINT96(*jx0,*jx1,*jx2, x9);
		CVT78_3WORD_DOUBLE_UINT96(*kx0,*kx1,*kx2, xa);
		CVT78_3WORD_DOUBLE_UINT96(*lx0,*lx1,*lx2, xb);
		CVT78_3WORD_DOUBLE_UINT96(*mx0,*mx1,*mx2, xc);
		CVT78_3WORD_DOUBLE_UINT96(*nx0,*nx1,*nx2, xd);
		CVT78_3WORD_DOUBLE_UINT96(*ox0,*ox1,*ox2, xe);
		CVT78_3WORD_DOUBLE_UINT96(*px0,*px1,*px2, xf);
		/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
		ADD96(x0,x0,x0);				ADD96(x8,x8,x8);
		ADD96(x1,x1,x1);				ADD96(x9,x9,x9);
		ADD96(x2,x2,x2);				ADD96(xa,xa,xa);
		ADD96(x3,x3,x3);				ADD96(xb,xb,xb);
		ADD96(x4,x4,x4);				ADD96(xc,xc,xc);
		ADD96(x5,x5,x5);				ADD96(xd,xd,xd);
		ADD96(x6,x6,x6);				ADD96(xe,xe,xe);
		ADD96(x7,x7,x7);				ADD96(xf,xf,xf);
			SUB96(x0,q0,x0);				SUB96(x8,q8,x8);
			SUB96(x1,q1,x1);				SUB96(x9,q9,x9);
			SUB96(x2,q2,x2);				SUB96(xa,qa,xa);
			SUB96(x3,q3,x3);				SUB96(xb,qb,xb);
			SUB96(x4,q4,x4);				SUB96(xc,qc,xc);
			SUB96(x5,q5,x5);				SUB96(xd,qd,xd);
			SUB96(x6,q6,x6);				SUB96(xe,qe,xe);
			SUB96(x7,q7,x7);				SUB96(xf,qf,xf);
		tmp0 = CMPEQ96(x0, ONE96);		tmp8 = CMPEQ96(x8, ONE96);
		tmp1 = CMPEQ96(x1, ONE96);		tmp9 = CMPEQ96(x9, ONE96);
		tmp2 = CMPEQ96(x2, ONE96);		tmpa = CMPEQ96(xa, ONE96);
		tmp3 = CMPEQ96(x3, ONE96);		tmpb = CMPEQ96(xb, ONE96);
		tmp4 = CMPEQ96(x4, ONE96);		tmpc = CMPEQ96(xc, ONE96);
		tmp5 = CMPEQ96(x5, ONE96);		tmpd = CMPEQ96(xd, ONE96);
		tmp6 = CMPEQ96(x6, ONE96);		tmpe = CMPEQ96(xe, ONE96);
		tmp7 = CMPEQ96(x7, ONE96);		tmpf = CMPEQ96(xf, ONE96);

		*checksum2 += x0.d0 + x1.d0 + x2.d0 + x3.d0 + x4.d0 + x5.d0 + x6.d0 + x7.d0 + x8.d0 + x9.d0 + xa.d0 + xb.d0 + xc.d0 + xd.d0 + xe.d0 + xf.d0;

		r += (tmp0 << 0) + (tmp8 << 0x8);
		r += (tmp1 << 1) + (tmp9 << 0x9);
		r += (tmp2 << 2) + (tmpa << 0xa);
		r += (tmp3 << 3) + (tmpb << 0xb);
		r += (tmp4 << 4) + (tmpc << 0xc);
		r += (tmp5 << 5) + (tmpd << 0xd);
		r += (tmp6 << 6) + (tmpe << 0xe);
		r += (tmp7 << 7) + (tmpf << 0xf);
		return r;
	}

  #endif

#endif

