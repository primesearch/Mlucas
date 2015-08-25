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
#include "align.h"

#define FAC_DEBUG 0
#if FAC_DEBUG
	char char_buf[1024], str0[64], str1[64];
#endif

/***********************************************************************************/
/*** 96-BIT INPUTS *****************************************************************/
/***********************************************************************************/
/*
Function to find 2^(-p) mod q, where p and q are 64 and 96-bit unsigned integers, respectively.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^96 (i.e. our MODQ
operation effects multiply modulo 2^96).

The key 3-operation sequence here is as follows:

	SQR_LOHI96(x,lo,hi);	// Input   x has 96 bits; outputs lo & hi have 96 bits
	MULL96(lo,qinv,lo);		// Inputs lo & qinv, and output (overwrites lo) have 96 bits
	MULH96(q,lo,lo);	// Inputs  q &   lo, and output (overwrites lo) have 96 bits
*/
uint96 twopmodq96(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k)
{
#if FAC_DEBUG
	int dbg = (p == 0);
	uint128 y;
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead7,  hi64;
	uint96 q, qhalf, qinv, x, lo, hi;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq96:\n");
#endif
	ASSERT(HERE, (p >> 63) == 0, "p must be < 2^63!");
	q.d0 = p+p;
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(q.d0, k,&q.d0,&q.d1);
#else
	MUL_LOHI64(q.d0, k, q.d0, q.d1);
#endif
	q.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	*checksum1 += q.d0;

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
		/* Extract leftmost 7 bits of pshift (if > 95, use the leftmost 6) and subtract from 96: */
		lead7 = ((pshift<<j) >> 57);	/* In [64,127] */
		if(lead7 > 95)
		{
			lead7 >>= 1;	/* Guarantees that lead7 in [48,95] */
			start_index =  64-j-6;	/* Use only the leftmost 6 bits */
		}
		else
			start_index =  64-j-7;

#if FAC_DEBUG
	if(dbg)	printf("lead7 = %u\n", (uint32)lead7);
#endif
		zshift = 95 - lead7;	/* In [0,47] */
		zshift <<= 1;			/* In [0,94]; Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

	/*
	!    Find modular inverse (mod 2^96) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q.d0 & (uint64)1) == 1, "twopmodq96 : (q.d0 & (uint64)1) == 1");
#endif
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
#if FAC_DEBUG
	MULL96(q, qinv, x);			MULL128(q, qinv, y);
	SUB96 (TWO96, x, x);		SUB128 (TWO128, y, y);
	MULL96(qinv, x, x);			MULL128(qinv, y, y);
	ASSERT(HERE, x.d1 == (y.d1 & 0x00000000ffffffff) && x.d0 == y.d0, "x.d1 == (y.d1 & 0x00000000ffffffff) && x.d0 == y.d0");
#endif
	/* qinv has 96 bits, but only the upper 32 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + __MULH64(q.d0, qinv.d0));
#else
	MULH64(q.d0, qinv.d0, hi64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + hi64);
#endif
	qinv.d1 &= 0x00000000ffffffff;	/* Only want the lower 32 bits here */

#if FAC_DEBUG
	ASSERT(HERE, qinv.d1 == x.d1 && qinv.d0 == x.d0, "twopmodq96 : qinv.d1 == x.d1 && qinv.d0 == x.d0");
	if(dbg) printf("q    = %s\n", &char_buf[convert_uint96_base10_char(char_buf, q   )]);
	if(dbg) printf("qinv = %s\n", &char_buf[convert_uint96_base10_char(char_buf, qinv)]);
#endif

	/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("zshift  = %u\n", zshift);
#endif
	LSHIFT96(qinv, zshift, lo);

#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo)]);
#endif

	MULH96(q,lo,lo);

#if FAC_DEBUG
	if(dbg) printf("q*lo/2^96 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB96(q, lo, x);	/* Put the result in lo (rather than x), to ease overflow check below */

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x)]);
#endif

	if((pshift >> j) & (uint64)1)
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT96(x, q), "twopmodq96 : CMPULT96(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT96(x, qhalf)){ ADD96(x, x, x); SUB96(x, q, x); }else{ ADD96(x, x, x); }
	}

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint96_base10_char(char_buf, x)]);
#endif
#if FAC_DEBUG
	if(CMPULT96(q, x)){ sprintf(char_buf, "twopmodq96 : (x0 = %s) >= (q = %s)", &str0[convert_uint96_base10_char(str0, x)], &str1[convert_uint96_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI96(x,lo,hi);	/* x=x*x;l=x%b;h=(x-l)/b; */
#if FAC_DEBUG
	if(dbg) printf("lo= %s\n", &char_buf[convert_uint96_base10_char(char_buf,lo)]);
	if(dbg) printf("hi= %s\n", &char_buf[convert_uint96_base10_char(char_buf,hi)]);
#endif
		MULL96(lo,qinv,lo);		/* l=i*l%b; */
#if FAC_DEBUG
	if(dbg) printf("lo = MULL96(lo,qinv) = %s\n", &char_buf[convert_uint96_base10_char(char_buf,lo)]);
#endif
		MULH96(q,lo,lo);		/* x=q*l%b;l=(l-x)/b; */
								/* x=(h-l)%q */
#if FAC_DEBUG
	if(dbg) printf("lo = MULH96(q,lo) = %s\n", &char_buf[convert_uint96_base10_char(char_buf,lo)]);
#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT96(hi, lo))
		{
			SUB96(q, lo, lo);
			ADD96(lo, hi, x);
#if FAC_DEBUG
	if(dbg) printf("q-l   = %10u, %20llu\n", lo.d1, lo.d0);
	if(dbg) printf("q-l+h = %10u, %20llu\n",  x.d1,  x.d0);
#endif
		}
		else
		{
			SUB96(hi, lo, x);
#if FAC_DEBUG
	if(dbg) printf("q=h-l = %10u, %20llu\n",  x.d1,  x.d0);
#endif
		}

#if FAC_DEBUG
if(dbg)printf("j = %2d, x = %s",j, &char_buf[convert_uint96_base10_char(char_buf, x)]);
#endif

		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT96(x, q), "twopmodq96 : CMPULT96(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x, qhalf)){ ADD96(x, x, x); SUB96(x, q, x); }else{ ADD96(x, x, x); }

#if FAC_DEBUG
	if(dbg) printf("*2= %s", &char_buf[convert_uint96_base10_char(char_buf, x)]);
#endif
		}
#if FAC_DEBUG
	if(dbg) printf("\n");
#endif
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD96(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
	SUB96(x,q,x);

#if FAC_DEBUG
if(dbg)printf("xout = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x)]);
#endif

	*checksum2 += x.d0;
	return x;
}

/******************************/
/*** 4-trial-factor version ***/
/******************************/
#undef YES_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
#endif

#ifdef YES_ASM

	uint64 twopmodq96_q4(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
	{
	#if FAC_DEBUG
		int dbg = (p != 0);
	/*	int dbg = STREQ(&char_buf[convert_uint64_base10_char(char_buf, p)], "0");*/
	#endif
		 int32 j;
		uint64 lead7, tmp0, tmp1, tmp2, tmp3, r;
		static uint64 psave = 0, pshift;
		static uint32 start_index, zshift, first_entry = TRUE;
		static uint64 *sm_arr = 0x0, *sm_ptr, *ptr64;
		static uint96 q0,q1,q2,q3
			,*qptr0,*qptr1,*qptr2,*qptr3
			,*qinv0,*qinv1,*qinv2,*qinv3
			,*qhalf0,*qhalf1,*qhalf2,*qhalf3
			,*x0,*x1,*x2,*x3
			,*lo0,*lo1,*lo2,*lo3
			,*hi0,*hi1,*hi2,*hi3
			,*ONE96_PTR;

		if(first_entry || p != psave)
		{
			first_entry = FALSE;
			psave  = p;
			pshift = p + 96;

			j = leadz64(pshift);
			/* Extract leftmost 7 bits of pshift (if > 95, use the leftmost 6) and subtract from 96: */
			lead7 = ((pshift<<j) >> 57);
			if(lead7 > 95) {
				lead7 >>= 1;
				start_index =  64-j-6;	/* Use only the leftmost 6 bits */
			} else {
				start_index =  64-j-7;
			}
		#if FAC_DEBUG
			if(dbg)	printf("twopmodq96_q4: lead7 = %u\n", (uint32)lead7);
		#endif
			zshift = 95 - lead7;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

			pshift = ~pshift;

			sm_arr = ALLOC_UINT64(sm_arr, 32);	if(!sm_arr){ ASSERT(HERE, 0, "FATAL: unable to allocate sm_arr!"); }
			sm_ptr = (uint64*)ALIGN_UINT64(sm_arr);
			ASSERT(HERE, ((uint64)sm_ptr & 0xf) == 0, "sm_ptr not 16-byte aligned!");
			/* Remember, rhese are pointers-to-uint128, so need an increment of 2 to span a memory slot: */
			qptr0  = (uint96*)(sm_ptr + 0x00);	qptr1  = (uint96*)(sm_ptr + 0x02);	qptr2  = (uint96*)(sm_ptr + 0x04);	qptr3  = (uint96*)(sm_ptr + 0x06);
			qinv0  = (uint96*)(sm_ptr + 0x08);	qinv1  = (uint96*)(sm_ptr + 0x0a);	qinv2  = (uint96*)(sm_ptr + 0x0c);	qinv3  = (uint96*)(sm_ptr + 0x0e);
			x0     = (uint96*)(sm_ptr + 0x10);	x1     = (uint96*)(sm_ptr + 0x12);	x2     = (uint96*)(sm_ptr + 0x14);	x3     = (uint96*)(sm_ptr + 0x16);
			lo0    = (uint96*)(sm_ptr + 0x18);	lo1    = (uint96*)(sm_ptr + 0x1a);	lo2    = (uint96*)(sm_ptr + 0x1c);	lo3    = (uint96*)(sm_ptr + 0x1e);
			qhalf0 = (uint96*)(sm_ptr + 0x20);	qhalf1 = (uint96*)(sm_ptr + 0x22);	qhalf2 = (uint96*)(sm_ptr + 0x24);	qhalf3 = (uint96*)(sm_ptr + 0x26);
			hi0    = (uint96*)(sm_ptr + 0x28);	hi1    = (uint96*)(sm_ptr + 0x2a);	hi2    = (uint96*)(sm_ptr + 0x2c);	hi3    = (uint96*)(sm_ptr + 0x2e);
			ONE96_PTR = (uint96*)(sm_ptr + 0x30);
			ptr64 = (uint64*)ONE96_PTR;	*ptr64++ = ONE96.d0;	*ptr64-- = ONE96.d1;
		}

	#if FAC_DEBUG
		if(dbg)printf("twopmodq96_q4:\n");
	#endif
		ASSERT(HERE, (p >> 63) == 0, "p must be < 2^63!");
		q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
		MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
		MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);
	
		q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
		q1.d0 += 1;
		q2.d0 += 1;
		q3.d0 += 1;

		*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0;

		ptr64 = (uint64*)qptr0;
		*ptr64++ = q0.d0;	*ptr64++ = q0.d1;
		*ptr64++ = q1.d0;	*ptr64++ = q1.d1;
		*ptr64++ = q2.d0;	*ptr64++ = q2.d1;
		*ptr64++ = q3.d0;	*ptr64++ = q3.d1;

		RSHIFT_FAST96_PTR(qptr0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
		RSHIFT_FAST96_PTR(qptr1, 1, qhalf1);
		RSHIFT_FAST96_PTR(qptr2, 1, qhalf2);
		RSHIFT_FAST96_PTR(qptr3, 1, qhalf3);

		/*
		!    Find modular inverse (mod 2^96) of q in preparation for modular multiply.
		*/
		/* q must be odd for Montgomery-style modmul to work: */
	#if FAC_DEBUG
		ASSERT(HERE, (q0.d0 & (uint64)1) == 1, "twopmodq96_q4 : (q0.d0 & (uint64)1) == 1");
	#endif
		qinv0->d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0->d1 = (uint64)0;
		qinv1->d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1->d1 = (uint64)0;
		qinv2->d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2->d1 = (uint64)0;
		qinv3->d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3->d1 = (uint64)0;

		/* Newton iteration involves repeated steps of form

			qinv = qinv*(2 - q*qinv);

		Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
		defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
		*/
		for(j = 0; j < 4; j++)
		{
			tmp0 = qptr0->d0 * qinv0->d0;
			tmp1 = qptr1->d0 * qinv1->d0;
			tmp2 = qptr2->d0 * qinv2->d0;
			tmp3 = qptr3->d0 * qinv3->d0;

			qinv0->d0 = qinv0->d0 * ((uint64)2 - tmp0);
			qinv1->d0 = qinv1->d0 * ((uint64)2 - tmp1);
			qinv2->d0 = qinv2->d0 * ((uint64)2 - tmp2);
			qinv3->d0 = qinv3->d0 * ((uint64)2 - tmp3);
		}

		/* Now that have bottom 64 bits of qinv, do one more Newton iteration
		using full 96-bit operands. See twopmodq96 for details on streamlining here.
		*/
		/* qinv has 96 bits, but only the upper 64 get modified here. */
		MULH64(qptr0->d0, qinv0->d0, tmp0);
		MULH64(qptr1->d0, qinv1->d0, tmp1);
		MULH64(qptr2->d0, qinv2->d0, tmp2);
		MULH64(qptr3->d0, qinv3->d0, tmp3);

		qinv0->d1 = -qinv0->d0 * (qptr0->d1 * qinv0->d0 + tmp0);
		qinv1->d1 = -qinv1->d0 * (qptr1->d1 * qinv1->d0 + tmp1);
		qinv2->d1 = -qinv2->d0 * (qptr2->d1 * qinv2->d0 + tmp2);
		qinv3->d1 = -qinv3->d0 * (qptr3->d1 * qinv3->d0 + tmp3);

		qinv0->d1 &= 0x00000000ffffffff;	/*  Only want the lower 32 bits here  */
		qinv1->d1 &= 0x00000000ffffffff;
		qinv2->d1 &= 0x00000000ffffffff;
		qinv3->d1 &= 0x00000000ffffffff;

	#if FAC_DEBUG
		if(dbg)
		{
			printf("q0 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, qptr0)]);
			printf("q1 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, qptr1)]);
			printf("q2 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, qptr2)]);
			printf("q3 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, qptr3)]);
			printf("\n");
			printf("qinv0 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, qinv0)]);
			printf("qinv1 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, qinv1)]);
			printf("qinv2 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, qinv2)]);
			printf("qinv3 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, qinv3)]);
		}
	#endif

		/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
		j = start_index-1;

	#if FAC_DEBUG
		if(dbg) printf("zshift  = %u\n", zshift);
	#endif
		/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
		LSHIFT96_PTR(qinv0, zshift, lo0);
		LSHIFT96_PTR(qinv1, zshift, lo1);
		LSHIFT96_PTR(qinv2, zshift, lo2);
		LSHIFT96_PTR(qinv3, zshift, lo3);
	#if FAC_DEBUG
		if(dbg) printf("lo0 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, lo0)]);
	#endif

		MULH96_PTR_q4(
		  qptr0, lo0, lo0
		, qptr1, lo1, lo1
		, qptr2, lo2, lo2
		, qptr3, lo3, lo3);
	#if FAC_DEBUG
		if(dbg) printf("q0*lo0/2^96 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, lo0)]);
	#endif

		/* hi = 0 in this instance, which simplifies things. */
		SUB96_PTR(qptr0, lo0, x0);
		SUB96_PTR(qptr1, lo1, x1);
		SUB96_PTR(qptr2, lo2, x2);
		SUB96_PTR(qptr3, lo3, x3);

		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT96_PTR(x0, qptr0), "twopmodq96_q4 : CMPULT96(x0,q0)");
		#endif
			if(CMPUGT96_PTR(x0, qhalf0)){ ADD96_PTR(x0, x0, x0); SUB96_PTR(x0, qptr0, x0); }else{ ADD96_PTR(x0, x0, x0); }
			if(CMPUGT96_PTR(x1, qhalf1)){ ADD96_PTR(x1, x1, x1); SUB96_PTR(x1, qptr1, x1); }else{ ADD96_PTR(x1, x1, x1); }
			if(CMPUGT96_PTR(x2, qhalf2)){ ADD96_PTR(x2, x2, x2); SUB96_PTR(x2, qptr2, x2); }else{ ADD96_PTR(x2, x2, x2); }
			if(CMPUGT96_PTR(x3, qhalf3)){ ADD96_PTR(x3, x3, x3); SUB96_PTR(x3, qptr3, x3); }else{ ADD96_PTR(x3, x3, x3); }
		#if FAC_DEBUG
			if(dbg) printf("x0 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, x0)]);
			ASSERT(HERE, CMPULT96_PTR(x0, qptr0), "twopmodq96_q4 : CMPULT96(x0,q0)");
		#endif
		}

		__asm__ volatile (\
			"movq	%[__q0],%%rsi	/* Use qptr0 as the base address throughout */\n\t"\
		"/* Load the x.d0's: */	\n\t"\
			"movq	0x80(%%rsi),%%r8 	/* Dereference the x0 pointer... */\n\t"\
			"movq	0x90(%%rsi),%%r9 	\n\t"\
			"movq	0xa0(%%rsi),%%r10	\n\t"\
			"movq	0xb0(%%rsi),%%r11	\n\t"\
		"/* Pure-ASM loop control: for(j = start_index-2; j >= 0; j--) */\n\t"\
			"movslq	%[__start_index], %%rcx		\n\t"\
			"subq $2,%%rcx						\n\t"\
			"test %%rcx, %%rcx					\n\t"\
			"jl LoopEnd		/* Skip if n < 0 */	\n\t"\
		"LoopBeg:								\n\t"\
		"/* SQR_LOHI96_q4(x*, lo*, hi*): */	\n\t"\
			"movq	%%r8 ,%%rax	/* Low 64 bits of x0-3 in r8-11. */\n\t"\
			"mulq	%%rax		\n\t"\
			"movq	%%rax,%%r8 	/* lo0 */	\n\t"\
			"movq	%%rdx,%%r12	/* hi0 */	\n\t"\
			"movq	%%r9 ,%%rax	\n\t"\
			"mulq	%%rax		\n\t"\
			"movq	%%rax,%%r9 	/* lo1 */	\n\t"\
			"movq	%%rdx,%%r13	/* hi1 */	\n\t"\
			"movq	%%r10,%%rax	\n\t"\
			"mulq	%%rax		\n\t"\
			"movq	%%rax,%%r10	/* lo2 */	\n\t"\
			"movq	%%rdx,%%r14	/* hi2 */	\n\t"\
			"movq	%%r11,%%rax	\n\t"\
			"mulq	%%rax		\n\t"\
			"movq	%%rax,%%r11	/* lo3 */	\n\t"\
			"movq	%%rdx,%%r15	/* hi3 */	\n\t"\
			"\n\t"\
		"/* Load each x.d0 in turn into the a-register and multiply it by the doubled high word (which is stored in a uint32): */	\n\t"\
			"\n\t"\
		"/* (lo0, hi0): */	\n\t"\
			"movq	0x80(%%rsi),%%rax	/* Move xlo into a-reg in prep for the mul */	\n\t"\
			"movl	0x88(%%rsi),%%edi	\n\t"\
			"shlq	$1,%%rdi		\n\t"\
			"mulq	%%rdi	/* 2*lo*hi in rax:rdx */	\n\t"\
			"shrq	$1,%%rdi		\n\t"\
			"imulq	%%rdi,%%rdi	/* hi*hi is just 64-bits, do in place */	\n\t"\
			"/* Assemble the 3-word result from (r8 ,r12,  0) + (  0,rax,rdx) + (  0,  0,rdi). Only carry is from (r12+rax), into (rdx+rdi): */	\n\t"\
			"addq	%%rax,%%r12	/* Add mid words, result in r12, CF bit set to indicate carryout */ \n\t"\
			"adcq	%%rdx,%%rdi	/* Add hi words + carryin, result in rdi, CF bit will be cleared */	\n\t"\
			"/* x0^2 in (r8 ,r12,rdi). Dump high half of result to 96-bit local hi0: */	\n\t"\
			"movq	%%r12,0x80(%%rsi)	/************* Dump back into x (rather than hi) for now *****************/\n\t"\
			"movq	%%rdi,0x88(%%rsi)	\n\t"\
			"\n\t"\
		"/* (lo1, hi1): */	\n\t"\
			"movq	0x90(%%rsi),%%rax	\n\t"\
			"movl	0x98(%%rsi),%%edi	\n\t"\
			"shlq	$1,%%rdi		\n\t"\
			"mulq	%%rdi		\n\t"\
			"shrq	$1,%%rdi		\n\t"\
			"imulq	%%rdi,%%rdi	\n\t"\
			"addq	%%rax,%%r13	\n\t"\
			"adcq	%%rdx,%%rdi	\n\t"\
			"/* x1^2 in (r9 ,r13,rdi). Dump high half of result to 96-bit local hi1: */	\n\t"\
			"movq	%%r13,0x90(%%rsi)	\n\t"\
			"movq	%%rdi,0x98(%%rsi)	\n\t"\
			"\n\t"\
		"/* (lo2, hi2): */	\n\t"\
			"movq	0xa0(%%rsi),%%rax	\n\t"\
			"movl	0xa8(%%rsi),%%edi	\n\t"\
			"shlq	$1,%%rdi		\n\t"\
			"mulq	%%rdi		\n\t"\
			"shrq	$1,%%rdi		\n\t"\
			"imulq	%%rdi,%%rdi	\n\t"\
			"addq	%%rax,%%r14	\n\t"\
			"adcq	%%rdx,%%rdi	\n\t"\
			"/* x2^2 in (r10,r14,rdi). Dump high half of result to 96-bit local hi2: */	\n\t"\
			"movq	%%r14,0xa0(%%rsi)	\n\t"\
			"movq	%%rdi,0xa8(%%rsi)	\n\t"\
			"\n\t"\
		"/* (lo3, hi3): */	\n\t"\
			"movq	0xb0(%%rsi),%%rax	\n\t"\
			"movl	0xb8(%%rsi),%%edi	\n\t"\
			"shlq	$1,%%rdi		\n\t"\
			"mulq	%%rdi		\n\t"\
			"shrq	$1,%%rdi		\n\t"\
			"imulq	%%rdi,%%rdi	\n\t"\
			"addq	%%rax,%%r15	\n\t"\
			"adcq	%%rdx,%%rdi	\n\t"\
			"/* x3^2 in (r11,r15,rdi). Dump high half of result to 96-bit local hi3: */	\n\t"\
			"movq	%%r15,0xb0(%%rsi)	\n\t"\
			"movq	%%rdi,0xb8(%%rsi)	\n\t"\
			"\n\t"\
		"/* MULL96_q4((qinv*, lo*, lo*): */	\n\t"\
			"\n\t"\
		"/* lo0 * qinv0: */	\n\t"\
			"movq	0x40(%%rsi),%%rdi	/* qinv0 */\n\t"\
			"movq	%%rdi,%%rax	\n\t"\
			"mulq	%%r8 	/* Use rax:rdx as product accumulator */\n\t"\
			"imulq	%%rdi,%%r12			/* qinv.lo * lo.hi. Hi 32 bits of r12 not cleared, but we keep only the lo 32 bits of the result here. */\n\t"\
			"addq	%%r12,%%rdx	\n\t"\
			"imulq	0x48(%%rsi),%%r8 	/* qinv.hi * lo.lo */\n\t"\
			"addq	%%r8 ,%%rdx	\n\t"\
			"movq	%%rax,%%r12	/* Move low 64 bits into free reg */\n\t"\
			"movl	%%edx,0xc8(%%rsi)	/* Write hi 32 bits into [__lo0] */\n\t"\
			"\n\t"\
		"/* lo1 * qinv1: */	\n\t"\
			"movq	0x50(%%rsi),%%rdi	\n\t"\
			"movq	%%rdi,%%rax	\n\t"\
			"mulq	%%r9 		\n\t"\
			"imulq	%%rdi,%%r13	\n\t"\
			"addq	%%r13,%%rdx	\n\t"\
			"imulq	0x58(%%rsi),%%r9	\n\t"\
			"addq	%%r9 ,%%rdx	\n\t"\
			"movq	%%rax,%%r13	\n\t"\
			"movl	%%edx,0xd8(%%rsi)	\n\t"\
		"/* lo2 * qinv2: */	\n\t"\
			"movq	0x60(%%rsi),%%rdi	\n\t"\
			"movq	%%rdi,%%rax	\n\t"\
			"mulq	%%r10		\n\t"\
			"imulq	%%rdi,%%r14	\n\t"\
			"addq	%%r14,%%rdx	\n\t"\
			"imulq	0x68(%%rsi),%%r10	\n\t"\
			"addq	%%r10,%%rdx	\n\t"\
			"movq	%%rax,%%r14	\n\t"\
			"movl	%%edx,0xe8(%%rsi)	\n\t"\
		"/* lo3 * qinv3: */	\n\t"\
			"movq	0x70(%%rsi),%%rdi	\n\t"\
			"movq	%%rdi,%%rax	\n\t"\
			"mulq	%%r11		\n\t"\
			"imulq	%%rdi,%%r15	\n\t"\
			"addq	%%r15,%%rdx	\n\t"\
			"imulq	0x78(%%rsi),%%r11	\n\t"\
			"addq	%%r11,%%rdx	\n\t"\
			"movq	%%rax,%%r15	\n\t"\
			"movl	%%edx,0xf8(%%rsi)	\n\t"\
			"\n\t"\
		"/* MULH96_q4((q*, lo*, lo*): Low 64 bits of of lo0-3 in r12-15: */	\n\t"\
			"\n\t"\
		"/* q0 * lo0: */	\n\t"\
			"movq	    (%%rsi),%%rax	/* q0 */\n\t"\
			"mulq	%%r12	/* lo.lo*q.lo in rax:rdx */\n\t"\
			"movl	0x08(%%rsi),%%edi	\n\t"\
			"movq	%%rdx,%%r8 	/* Discard low 64 bits [rax] */\n\t"\
			"movq	%%r12,%%rax	\n\t"\
			"mulq	%%rdi	/* lo.lo*q.hi in rax:rdx */\n\t"\
			"xorq	%%r12,%%r12	\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"movl	0xc8(%%rsi),%%eax	/* Can't do imulq __lohi,[reg] because that does a sign-extended load of __lohi */\n\t"\
			"imulq	%%rax,%%rdi			/* ...so first load __lohi into low 32 bits of a-reg, then compute lo.hi*q.hi. */\n\t"\
			"addq	%%rdi,%%r12	\n\t"\
			"mulq	    (%%rsi)		/* q.lo*lo.hi in rax:rdx */\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
			"/* Compare upper 64 bits against hi output of squaring, which is also stored in not-yet-right-32-shifted 128-bit form: */\n\t"\
			"movq	0x80(%%rsi),%%rax	/* h.d0 */\n\t"\
			"movq	0x88(%%rsi),%%rdx	/* h.d1 */\n\t"\
			"subq	%%r8 ,%%rax	/* rax = (h-l).d0 */\n\t"\
			"sbbq	%%r12,%%rdx	/* rdx = (h-l).d1 */\n\t"\
			"/* If there was a borrow out of the hi-word subtract (h-l), use it to create a bitmask of all ones in the free rdi-register: */\n\t"\
			"sbbq	%%rdi,%%rdi	\n\t"\
			"/* Right-shift the 2-word result 32 bits: */\n\t"\
			"shrdq	$32,%%rdx,%%rax	/* Only keep high 96-bit output */\n\t"\
			"/* Initial impulse here is to use arithmetic rshift (sarq) of hi word, since if h<l, carryout of low 64 bits -> hi=2^32 = 0x100000000, */\n\t"\
			"/* which would need to be explicitly zeroed prior to the conditional doubling step. But arith-shift needs the hi bit of rdx to reflect */\n\t"\
			"/* any carryout, which will not be the case if (h-l) has 96 sig. bits (i.e. 128 bits in the present << 32 form. To properly handle that, */\n\t"\
			"/* need either another 2-word shift 'shrdq $32,%%rdi,%%rdx' here or a later explicit zeroing of the hi 32 bits ... prefer the latter. */\n\t"\
			"shrq	$32,%%rdx			\n\t"\
			"/* Apply mask to lo/hi words of q and add result to (h-l): */\n\t"\
			"movq	    (%%rsi),%%r8 	\n\t"\
			"movq	0x08(%%rsi),%%r12	\n\t"\
			"andq	%%rdi,%%r8 	\n\t"\
			"andq	%%rdi,%%r12	\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"\n\t"\
		"/* q1 * lo1: */	\n\t"\
			"movq	0x10(%%rsi),%%rax	\n\t"\
			"mulq	%%r13	\n\t"\
			"movl	0x18(%%rsi),%%edi	\n\t"\
			"movq	%%rdx,%%r9 	\n\t"\
			"movq	%%r13,%%rax	\n\t"\
			"mulq	%%rdi	\n\t"\
			"xorq	%%r13,%%r13	\n\t"\
			"addq	%%rax,%%r9 	\n\t"\
			"adcq	%%rdx,%%r13	\n\t"\
			"movl	0xd8(%%rsi),%%eax	\n\t"\
			"imulq	%%rax,%%rdi	\n\t"\
			"addq	%%rdi,%%r13	\n\t"\
			"mulq	0x10(%%rsi)	\n\t"\
			"addq	%%rax,%%r9 	\n\t"\
			"adcq	%%rdx,%%r13	\n\t"\
			"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
			"movq	0x90(%%rsi),%%rax	\n\t"\
			"movq	0x98(%%rsi),%%rdx	\n\t"\
			"subq	%%r9 ,%%rax		\n\t"\
			"sbbq	%%r13,%%rdx		\n\t"\
			"sbbq	%%rdi,%%rdi	/* bitmask for q */\n\t"\
			"shrdq	$32,%%rdx,%%rax		\n\t"\
			"shrq	$32,%%rdx			\n\t"\
			"movq	0x10(%%rsi),%%r9 	\n\t"\
			"movq	0x18(%%rsi),%%r13	\n\t"\
			"andq	%%rdi,%%r9 	\n\t"\
			"andq	%%rdi,%%r13	\n\t"\
			"addq	%%rax,%%r9 	\n\t"\
			"adcq	%%rdx,%%r13	\n\t"\
			"\n\t"\
		"/* q2 * lo2: */	\n\t"\
			"movq	0x20(%%rsi),%%rax	\n\t"\
			"mulq	%%r14	\n\t"\
			"movl	0x28(%%rsi),%%edi	\n\t"\
			"movq	%%rdx,%%r10	\n\t"\
			"movq	%%r14,%%rax	\n\t"\
			"mulq	%%rdi	\n\t"\
			"xorq	%%r14,%%r14	\n\t"\
			"addq	%%rax,%%r10	\n\t"\
			"adcq	%%rdx,%%r14	\n\t"\
			"movl	0xe8(%%rsi),%%eax	\n\t"\
			"imulq	%%rax,%%rdi	\n\t"\
			"addq	%%rdi,%%r14	\n\t"\
			"mulq	0x20(%%rsi)	\n\t"\
			"addq	%%rax,%%r10	\n\t"\
			"adcq	%%rdx,%%r14	\n\t"\
			"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
			"movq	0xa0(%%rsi),%%rax	\n\t"\
			"movq	0xa8(%%rsi),%%rdx	\n\t"\
			"subq	%%r10,%%rax		\n\t"\
			"sbbq	%%r14,%%rdx		\n\t"\
			"sbbq	%%rdi,%%rdi	/* bitmask for q */\n\t"\
			"shrdq	$32,%%rdx,%%rax		\n\t"\
			"shrq	$32,%%rdx			\n\t"\
			"movq	0x20(%%rsi),%%r10	\n\t"\
			"movq	0x28(%%rsi),%%r14	\n\t"\
			"andq	%%rdi,%%r10	\n\t"\
			"andq	%%rdi,%%r14	\n\t"\
			"addq	%%rax,%%r10	\n\t"\
			"adcq	%%rdx,%%r14	\n\t"\
			"\n\t"\
		"/* q3 * lo3: */	\n\t"\
			"movq	0x30(%%rsi),%%rax	\n\t"\
			"mulq	%%r15	\n\t"\
			"movl	0x38(%%rsi),%%edi	\n\t"\
			"movq	%%rdx,%%r11	\n\t"\
			"movq	%%r15,%%rax	\n\t"\
			"mulq	%%rdi	\n\t"\
			"xorq	%%r15,%%r15	\n\t"\
			"addq	%%rax,%%r11	\n\t"\
			"adcq	%%rdx,%%r15	\n\t"\
			"movl	0xf8(%%rsi),%%eax	\n\t"\
			"imulq	%%rax,%%rdi	\n\t"\
			"addq	%%rdi,%%r15	\n\t"\
			"mulq	0x30(%%rsi)	\n\t"\
			"addq	%%rax,%%r11	\n\t"\
			"adcq	%%rdx,%%r15	\n\t"\
			"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
			"movq	0xb0(%%rsi),%%rax	\n\t"\
			"movq	0xb8(%%rsi),%%rdx	\n\t"\
			"subq	%%r11,%%rax		\n\t"\
			"sbbq	%%r15,%%rdx		\n\t"\
			"sbbq	%%rdi,%%rdi	/* bitmask for q */\n\t"\
			"shrdq	$32,%%rdx,%%rax		\n\t"\
			"shrq	$32,%%rdx			\n\t"\
			"movq	0x30(%%rsi),%%r11	\n\t"\
			"movq	0x38(%%rsi),%%r15	\n\t"\
			"andq	%%rdi,%%r11	\n\t"\
			"andq	%%rdi,%%r15	\n\t"\
			"addq	%%rax,%%r11	\n\t"\
			"adcq	%%rdx,%%r15	\n\t"\
			"\n\t"\
	"/* If current bit of pshift == 1, double each output modulo q: */	\n\t"\
			"/* if((pshift >> j) & (uint64)1) { */	\n\t"\
			"movl	%[__pshift],%%eax		\n\t"\
			"shrq	%%cl,%%rax				\n\t"\
			"andq	$0x1,%%rax				\n\t"\
		"je	twopmodq96_q4_pshiftjmp			\n\t"\
			"\n\t"\
		"/* if h<l carryout of low 64 bits gives hi=2^32 = 0x100000000, need to zero upper 32 bits prior to double step: */\n\t"\
			"movq	$-1,%%rdi	\n\t"\
			"shrq	$32,%%rdi	\n\t"\
			"andq	%%rdi,%%r12	\n\t"\
			"andq	%%rdi,%%r13	\n\t"\
			"andq	%%rdi,%%r14	\n\t"\
			"andq	%%rdi,%%r15	\n\t"\
		"/* x0: */	\n\t"\
			"addq	%%r8 ,%%r8 	\n\t"\
			"adcq	%%r12,%%r12	\n\t"\
			"movq	    (%%rsi),%%rax	\n\t"\
			"movq	0x08(%%rsi),%%rdx	\n\t"\
			"subq	%%rax,%%r8 		\n\t"\
			"sbbq	%%rdx,%%r12		\n\t"\
			"sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
			"andq	%%rdi,%%rax	\n\t"\
			"andq	%%rdi,%%rdx	\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"\n\t"\
		"/* x1: */	\n\t"\
			"addq	%%r9 ,%%r9 	\n\t"\
			"adcq	%%r13,%%r13	\n\t"\
			"movq	0x10(%%rsi),%%rax	\n\t"\
			"movq	0x18(%%rsi),%%rdx	\n\t"\
			"subq	%%rax,%%r9 		\n\t"\
			"sbbq	%%rdx,%%r13		\n\t"\
			"sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
			"andq	%%rdi,%%rax	\n\t"\
			"andq	%%rdi,%%rdx	\n\t"\
			"addq	%%rax,%%r9 	\n\t"\
			"adcq	%%rdx,%%r13	\n\t"\
			"\n\t"\
		"/* x2: */	\n\t"\
			"addq	%%r10,%%r10	\n\t"\
			"adcq	%%r14,%%r14	\n\t"\
			"movq	0x20(%%rsi),%%rax	\n\t"\
			"movq	0x28(%%rsi),%%rdx	\n\t"\
			"subq	%%rax,%%r10		\n\t"\
			"sbbq	%%rdx,%%r14		\n\t"\
			"sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
			"andq	%%rdi,%%rax	\n\t"\
			"andq	%%rdi,%%rdx	\n\t"\
			"addq	%%rax,%%r10	\n\t"\
			"adcq	%%rdx,%%r14	\n\t"\
			"\n\t"\
		"/* x3: */	\n\t"\
			"addq	%%r11,%%r11	\n\t"\
			"adcq	%%r15,%%r15	\n\t"\
			"movq	0x30(%%rsi),%%rax	\n\t"\
			"movq	0x38(%%rsi),%%rdx	\n\t"\
			"subq	%%rax,%%r11		\n\t"\
			"sbbq	%%rdx,%%r15		\n\t"\
			"sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t"\
			"andq	%%rdi,%%rax	\n\t"\
			"andq	%%rdi,%%rdx	\n\t"\
			"addq	%%rax,%%r11	\n\t"\
			"adcq	%%rdx,%%r15	\n\t"\
		"twopmodq96_q4_pshiftjmp:					\n\t"\
			"/* } endif((pshift >> j) & (uint64)1) */						\n\t"\
			"movq	%%r8 ,0x80(%%rsi)	/* Write lo 64 bits into [__x.d0] */\n\t"\
			"movq	%%r12,0x88(%%rsi)	/* Write hi 32 bits into [__x.d1] */\n\t"\
			"movq	%%r9 ,0x90(%%rsi)	\n\t"\
			"movq	%%r13,0x98(%%rsi)	\n\t"\
			"movq	%%r10,0xa0(%%rsi)	\n\t"\
			"movq	%%r14,0xa8(%%rsi)	\n\t"\
			"movq	%%r11,0xb0(%%rsi)	\n\t"\
			"movq	%%r15,0xb8(%%rsi)	\n\t"\
			"subq	$1,%%rcx	/* j-- */		\n\t"\
			"cmpq	$0,%%rcx	/* compare j vs 0 */\n\t"\
			"jge	LoopBeg	/* if (j >= 0), Loop */	\n\t"\
		"LoopEnd:							\n\t"\
			:	/* outputs: none */\
			: [__q0] "m" (qptr0)	/* All inputs from memory addresses here */\
			 ,[__pshift] "m" (pshift)	\
			 ,[__start_index] "m" (start_index)	\
			: "cc","memory","rax","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"		/* Clobbered registers */\
		);
		/* endfor() */

		/*...Double and return.	These are specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		*/
		ADD96_PTR(x0 ,x0, x0);
		ADD96_PTR(x1 ,x1, x1);
		ADD96_PTR(x2 ,x2, x2);
		ADD96_PTR(x3 ,x3, x3);

		SUB96_PTR(x0, qptr0, x0);
		SUB96_PTR(x1, qptr1, x1);
		SUB96_PTR(x2, qptr2, x2);
		SUB96_PTR(x3, qptr3, x3);

		tmp0 = CMPEQ96_PTR(x0, ONE96_PTR);	*checksum2 += x0->d0;
		tmp1 = CMPEQ96_PTR(x1, ONE96_PTR);	*checksum2 += x1->d0;
		tmp2 = CMPEQ96_PTR(x2, ONE96_PTR);	*checksum2 += x2->d0;
		tmp3 = CMPEQ96_PTR(x3, ONE96_PTR);	*checksum2 += x3->d0;
	#if FAC_DEBUG
		if(dbg)
		{
			printf("xout0 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, x0)]);
			printf("xout1 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, x1)]);
			printf("xout2 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, x2)]);
			printf("xout3 = %s\n", &char_buf[convert_uint96ptr_base10_char(char_buf, x3)]);
			exit(0);
		}
	#endif

		r = tmp0;
		r += tmp1 << 1;
		r += tmp2 << 2;
		r += tmp3 << 3;
		return r;
	}

#else	/* Reference version using multiword macros: */

	uint64 twopmodq96_q4(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
	{
	#if FAC_DEBUG
		int dbg = (p != 0);
	/*	int dbg = STREQ(&char_buf[convert_uint64_base10_char(char_buf, p)], "0");*/
	#endif
		 int32 j;
	#undef	NO_BRANCH
	#define	NO_BRANCH	0
	#if	NO_BRANCH
		uint64 fq_or_nil_lo64[2],gq_or_nil_lo64[2],hq_or_nil_lo64[2],iq_or_nil_lo64[2];
		uint32 fq_or_nil_hi32[2],gq_or_nil_hi32[2],hq_or_nil_hi32[2],iq_or_nil_hi32[2];
		int fidx,gidx,hidx,iidx;
	#endif
		uint64 lead7, tmp0, tmp1, tmp2, tmp3, r;
		static uint64 psave = 0, pshift;
		static uint32 start_index, zshift, first_entry = TRUE;
		uint96 q0, q1, q2, q3
			, qinv0, qinv1, qinv2, qinv3
			, qhalf0, qhalf1, qhalf2, qhalf3
			, x0, x1, x2, x3
			, lo0, lo1, lo2, lo3
			, hi0, hi1, hi2, hi3;

		if(first_entry || p != psave)
		{
			first_entry = FALSE;
			psave  = p;
			pshift = p + 96;

			j = leadz64(pshift);
			/* Extract leftmost 7 bits of pshift (if > 95, use the leftmost 6) and subtract from 96: */
			lead7 = ((pshift<<j) >> 57);
			if(lead7 > 95) {
				lead7 >>= 1;
				start_index =  64-j-6;	/* Use only the leftmost 6 bits */
			} else {
				start_index =  64-j-7;
			}
		#if FAC_DEBUG
			if(dbg)	printf("twopmodq96_q4: lead7 = %u\n", (uint32)lead7);
		#endif
			zshift = 95 - lead7;
			zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

			pshift = ~pshift;
		}

	#if FAC_DEBUG
	if(dbg)printf("twopmodq96_q4:\n");
	#endif
		ASSERT(HERE, (p >> 63) == 0, "p must be < 2^63!");
		q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;
		MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
		MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
		MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
		MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);
	
		q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
		q1.d0 += 1;
		q2.d0 += 1;
		q3.d0 += 1;

		*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0;

	#if	NO_BRANCH
		fq_or_nil_lo64[0] = 0ull; fq_or_nil_lo64[1] = q0.d0; fq_or_nil_hi32[0] = 0; fq_or_nil_hi32[1] = (uint32)q0.d1;
		gq_or_nil_lo64[0] = 0ull; gq_or_nil_lo64[1] = q1.d0; gq_or_nil_hi32[0] = 0; gq_or_nil_hi32[1] = (uint32)q1.d1;
		hq_or_nil_lo64[0] = 0ull; hq_or_nil_lo64[1] = q2.d0; hq_or_nil_hi32[0] = 0; hq_or_nil_hi32[1] = (uint32)q2.d1;
		iq_or_nil_lo64[0] = 0ull; iq_or_nil_lo64[1] = q3.d0; iq_or_nil_hi32[0] = 0; iq_or_nil_hi32[1] = (uint32)q3.d1;
	#endif

		RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
		RSHIFT_FAST96(q1, 1, qhalf1);
		RSHIFT_FAST96(q2, 1, qhalf2);
		RSHIFT_FAST96(q3, 1, qhalf3);

		/*
		!    Find modular inverse (mod 2^96) of q in preparation for modular multiply.
		*/
		/* q must be odd for Montgomery-style modmul to work: */
	#if FAC_DEBUG
		ASSERT(HERE, (q0.d0 & (uint64)1) == 1, "twopmodq96_q4 : (q0.d0 & (uint64)1) == 1");
	#endif
		qinv0.d0  = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1  = (uint64)0;
		qinv1.d0  = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1  = (uint64)0;
		qinv2.d0  = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d1  = (uint64)0;
		qinv3.d0  = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d1  = (uint64)0;

		/* Newton iteration involves repeated steps of form

			qinv = qinv*(2 - q*qinv);

		Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
		defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
		*/
		for(j = 0; j < 4; j++)
		{
			tmp0 = q0.d0 * qinv0.d0;
			tmp1 = q1.d0 * qinv1.d0;
			tmp2 = q2.d0 * qinv2.d0;
			tmp3 = q3.d0 * qinv3.d0;

			qinv0.d0 = qinv0.d0 * ((uint64)2 - tmp0);
			qinv1.d0 = qinv1.d0 * ((uint64)2 - tmp1);
			qinv2.d0 = qinv2.d0 * ((uint64)2 - tmp2);
			qinv3.d0 = qinv3.d0 * ((uint64)2 - tmp3);
		}

		/* Now that have bottom 64 bits of qinv, do one more Newton iteration
		using full 96-bit operands. See twopmodq96 for details on streamlining here.
		*/
		/* qinv has 96 bits, but only the upper 64 get modified here. */
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

		qinv0.d1 &= 0x00000000ffffffff;	/* Only want the lower 32 bits here */
		qinv1.d1 &= 0x00000000ffffffff;
		qinv2.d1 &= 0x00000000ffffffff;
		qinv3.d1 &= 0x00000000ffffffff;

	  #if FAC_DEBUG
		if(dbg)
		{
			printf("q0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, q0)]);
			printf("q1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, q1)]);
			printf("q2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, q2)]);
			printf("q3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, q3)]);
			printf("\n");
			printf("qinv0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, qinv0)]);
			printf("qinv1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, qinv1)]);
			printf("qinv2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, qinv2)]);
			printf("qinv3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, qinv3)]);
		}
	  #endif
	#endif

		/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
		j = start_index-1;

	#if FAC_DEBUG
		if(dbg) printf("zshift  = %u\n", zshift);
	#endif
		/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
		LSHIFT96(qinv0, zshift, lo0);
		LSHIFT96(qinv1, zshift, lo1);
		LSHIFT96(qinv2, zshift, lo2);
		LSHIFT96(qinv3, zshift, lo3);
	#if FAC_DEBUG
		if(dbg) printf("lo0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo0)]);
	#endif

		MULH96_q4(
		  q0, lo0, lo0
		, q1, lo1, lo1
		, q2, lo2, lo2
		, q3, lo3, lo3);
	#if FAC_DEBUG
		if(dbg) printf("q0*lo0/2^96 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo0)]);
	#endif

		/* hi = 0 in this instance, which simplifies things. */
		SUB96(q0, lo0, x0);
		SUB96(q1, lo1, x1);
		SUB96(q2, lo2, x2);
		SUB96(q3, lo3, x3);

		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT96(x0, q0), "twopmodq96_q4 : CMPULT96(x0,q0)");
		#endif
			if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
			if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
			if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
			if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
		#if FAC_DEBUG
			if(dbg) printf("x0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x0)]);
			ASSERT(HERE, CMPULT96(x0, q0), "twopmodq96_q4 : CMPULT96(x0,q0)");
			if(CMPULT96(q0, x0)){
				sprintf(char_buf, "twopmodq96_q4 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint96_base10_char(str0, x0)], &str1[convert_uint96_base10_char(str1, q0)] );
			//	DBG_WARN(HERE, char_buf, STATFILE, !restart);
			}
		#endif
		}

		for(j = start_index-2; j >= 0; j--)
		{
		#if FAC_DEBUG
			if(dbg)printf("j = %2d:\n", j);
		#endif
			/*...x^2 mod q is returned in x. */
			SQR_LOHI96_q4(
			  x0, lo0, hi0
			, x1, lo1, hi1
			, x2, lo2, hi2
			, x3, lo3, hi3);
		#if FAC_DEBUG
			if(dbg) {
				printf("Multiply output #1.lo0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo0)]);
				printf("Multiply output #1.lo1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo1)]);
				printf("Multiply output #1.lo2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo2)]);
				printf("Multiply output #1.lo3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo3)]);

				printf("Multiply output #1.hi0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, hi0)]);
				printf("Multiply output #1.hi1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, hi1)]);
				printf("Multiply output #1.hi2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, hi2)]);
				printf("Multiply output #1.hi3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, hi3)]);
			}
		#endif

			MULL96_q4(
			  lo0, qinv0, lo0
			, lo1, qinv1, lo1
			, lo2, qinv2, lo2
			, lo3, qinv3, lo3);
		#if FAC_DEBUG
			if(dbg) {
				printf("Multiply output #2.0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo0)]);
				printf("Multiply output #2.1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo1)]);
				printf("Multiply output #2.2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo2)]);
				printf("Multiply output #2.3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo3)]);
			}
		#endif

			MULH96_q4(
			  q0, lo0, lo0
			, q1, lo1, lo1
			, q2, lo2, lo2
			, q3, lo3, lo3);
		#if FAC_DEBUG
			if(dbg) {
				printf("Multiply output #3.0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo0)]);
				printf("Multiply output #3.1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo1)]);
				printf("Multiply output #3.2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo2)]);
				printf("Multiply output #3.3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, lo3)]);
			}
		#endif

			/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */

	#if NO_BRANCH	/* Branchless versus branched post-multiply (h-l) and double-if-current-bit-set switch */

			/* First compare upper 32 bits of the respective 96-bit outputs: */
			fidx = (hi0.d1 < lo0.d1);
			gidx = (hi1.d1 < lo1.d1);
			hidx = (hi2.d1 < lo2.d1);
			iidx = (hi3.d1 < lo3.d1);
			/* If high-32 compare = 0, the upper pieces might still be equal, in which case we need the result of the low-64-bit compare.
			We want to avoid a branch here, so use the high-part compare as a mask, which discards (mask = 0) the == and low-64 compare if the
			high-part cpmpare returned 'true' and keeps (mask = -1) them if 'false': */
			fidx += (fidx-1) & ((hi0.d1 == lo0.d1) && (hi0.d0 < lo0.d0));
			gidx += (gidx-1) & ((hi1.d1 == lo1.d1) && (hi1.d0 < lo1.d0));
			hidx += (hidx-1) & ((hi2.d1 == lo2.d1) && (hi2.d0 < lo2.d0));
			iidx += (iidx-1) & ((hi3.d1 == lo3.d1) && (hi3.d0 < lo3.d0));

			/* Double-and-mod conditional below needs original value of hi* terms */
			x0.d1 = hi0.d1 - lo0.d1;
			x1.d1 = hi1.d1 - lo1.d1;
			x2.d1 = hi2.d1 - lo2.d1;
			x3.d1 = hi3.d1 - lo3.d1;

			x0.d0 = hi0.d0 - lo0.d0;
			x1.d0 = hi1.d0 - lo1.d0;
			x2.d0 = hi2.d0 - lo2.d0;
			x3.d0 = hi3.d0 - lo3.d0;

			/* Adjust high word for any borrow: */
			x0.d1 -= (x0.d0 > hi0.d0);
			x1.d1 -= (x1.d0 > hi1.d0);
			x2.d1 -= (x2.d0 > hi2.d0);
			x3.d1 -= (x3.d0 > hi3.d0);

			/* The *idx value tells us whether to add q or 0 to the corresponding h-l result: */
			x0.d0 += fq_or_nil_lo64[fidx];
			x1.d0 += gq_or_nil_lo64[gidx];
			x2.d0 += hq_or_nil_lo64[hidx];
			x3.d0 += iq_or_nil_lo64[iidx];

			x0.d1 += x0.d0 < fq_or_nil_lo64[fidx];	// Carries into high word
			x1.d1 += x1.d0 < gq_or_nil_lo64[gidx];
			x2.d1 += x2.d0 < hq_or_nil_lo64[hidx];
			x3.d1 += x3.d0 < iq_or_nil_lo64[iidx];

			x0.d1 += fq_or_nil_hi32[fidx];
			x1.d1 += gq_or_nil_hi32[gidx];
			x2.d1 += hq_or_nil_hi32[hidx];
			x3.d1 += iq_or_nil_hi32[iidx];

			/* Branchless version of the double-and-subtract-q-if-(2x >= q) sequence:
				x = x + x - ((-(x > qhalf)) & q);
			*/
			if((pshift >> j) & (uint64)1)
			{
			#if FAC_DEBUG
				ASSERT(HERE, CMPULT96(x0, q0), "twopmodq96_q4 : CMPULT96(x0,q0)");
				ASSERT(HERE, CMPULT96(x1, q1), "twopmodq96_q4 : CMPULT96(x1,q1)");
				ASSERT(HERE, CMPULT96(x2, q2), "twopmodq96_q4 : CMPULT96(x2,q2)");
				ASSERT(HERE, CMPULT96(x3, q3), "twopmodq96_q4 : CMPULT96(x3,q3)");
			#endif
				/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
				if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
				if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
				if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
				if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
			#if 0
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

				fx0 -= fq_or_nil_lo26[fidx];
				gx0 -= gq_or_nil_lo26[gidx];
				hx0 -= hq_or_nil_lo26[hidx];
				ix0 -= iq_or_nil_lo26[iidx];
			#endif
			#if FAC_DEBUG
				if(dbg) printf("*2= %s", &char_buf[convert_uint96_base10_char(char_buf, x0)]);
			#endif
			}
		#if FAC_DEBUG
			if(dbg) printf("\n");
		#endif

	#else

			/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */
			if(CMPULT96(hi0, lo0)) { SUB96(q0, lo0, lo0);	ADD96(lo0, hi0, x0); } else { SUB96(hi0, lo0, x0); }
			if(CMPULT96(hi1, lo1)) { SUB96(q1, lo1, lo1);	ADD96(lo1, hi1, x1); } else { SUB96(hi1, lo1, x1); }
			if(CMPULT96(hi2, lo2)) { SUB96(q2, lo2, lo2);	ADD96(lo2, hi2, x2); } else { SUB96(hi2, lo2, x2); }
			if(CMPULT96(hi3, lo3)) { SUB96(q3, lo3, lo3);	ADD96(lo3, hi3, x3); } else { SUB96(hi3, lo3, x3); }
		#if FAC_DEBUG
			printf("(h-l) mod q #0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x0)]);
			printf("(h-l) mod q #1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x1)]);
			printf("(h-l) mod q #2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x2)]);
			printf("(h-l) mod q #3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x3)]);
		#endif

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if((pshift >> j) & (uint64)1)
			{
			#if FAC_DEBUG
				ASSERT(HERE, CMPULT96(x0, q0), "twopmodq96_q4 : CMPULT96(x0,q0)");
				ASSERT(HERE, CMPULT96(x1, q1), "twopmodq96_q4 : CMPULT96(x1,q1)");
				ASSERT(HERE, CMPULT96(x2, q2), "twopmodq96_q4 : CMPULT96(x2,q2)");
				ASSERT(HERE, CMPULT96(x3, q3), "twopmodq96_q4 : CMPULT96(x3,q3)");
			#endif
				if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
				if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
				if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
				if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
			#if FAC_DEBUG
				if(dbg) {
					printf("x0 = %s, *2 .\n", &char_buf[convert_uint96_base10_char(char_buf, x0)]);
					printf("x1 = %s, *2 .\n", &char_buf[convert_uint96_base10_char(char_buf, x1)]);
					printf("x2 = %s, *2 .\n", &char_buf[convert_uint96_base10_char(char_buf, x2)]);
					printf("x3 = %s, *2 .\n", &char_buf[convert_uint96_base10_char(char_buf, x3)]);
				}
			#endif
			} else {
			#if FAC_DEBUG
				if(dbg) {
					printf("x0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x0)]);
					printf("x1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x1)]);
					printf("x2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x2)]);
					printf("x3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x3)]);
				}
			#endif
			}
		#if FAC_DEBUG
			if(dbg) printf("\n");
		#endif

	#endif	/* Branchless versus branched post-multiply (h-l) and doubleif-vurrent-bit-set switch */
		}	/* endfor() */
		/*...Double and return.	These are specialized for the case
		where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
		*/
		ADD96(x0 ,x0, x0);
		ADD96(x1 ,x1, x1);
		ADD96(x2 ,x2, x2);
		ADD96(x3 ,x3, x3);

		SUB96(x0, q0, x0);
		SUB96(x1, q1, x1);
		SUB96(x2, q2, x2);
		SUB96(x3, q3, x3);

		tmp0 = CMPEQ96(x0, ONE96);	*checksum2 += x0.d0;
		tmp1 = CMPEQ96(x1, ONE96);	*checksum2 += x1.d0;
		tmp2 = CMPEQ96(x2, ONE96);	*checksum2 += x2.d0;
		tmp3 = CMPEQ96(x3, ONE96);	*checksum2 += x3.d0;
	#if FAC_DEBUG
		if(dbg)
		{
			printf("xout0 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x0)]);
			printf("xout1 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x1)]);
			printf("xout2 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x2)]);
			printf("xout3 = %s\n", &char_buf[convert_uint96_base10_char(char_buf, x3)]);
			exit(0);
		}
	#endif

		r = tmp0;
		r += tmp1 << 1;
		r += tmp2 << 2;
		r += tmp3 << 3;
		return r;
	}

#endif


/*** 8-trial-factor version ***/
uint64 twopmodq96_q8(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint64_base10_char(char_buf, p)], "0");
#endif
	 int32 j;
	uint64 lead7, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r;
	uint96 q0, q1, q2, q3, q4, q5, q6, q7
		, qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, qhalf0, qhalf1, qhalf2, qhalf3, qhalf4, qhalf5, qhalf6, qhalf7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
		, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
	if(dbg)printf("twopmodq96_q8:\n");
#endif
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

	*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0 + q4.d0 + q5.d0 + q6.d0 + q7.d0;

	RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
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
		pshift = p + 96;

		j = leadz64(pshift);
		/* Extract leftmost 7 bits of pshift (if > 95, use the leftmost 6) and subtract from 96: */
		lead7 = ((pshift<<j) >> 57);
		if(lead7 > 95)
		{
			lead7 >>= 1;
			start_index =  64-j-6;	/* Use only the leftmost 6 bits */
		}
		else
			start_index =  64-j-7;

#if FAC_DEBUG
	if(dbg)	printf("lead7 = %u\n", (uint32)lead7);
#endif
		zshift = 95 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

	ASSERT(HERE, (p >> 63) == 0, "p must be < 2^63!");
	q0.d0 = q1.d0 = q2.d0 = q3.d0 = q4.d0 = q5.d0 = q6.d0 = q7.d0 = p+p;
	q0.d1 = q1.d1 = q2.d1 = q3.d1 = q4.d1 = q5.d1 = q6.d1 = q7.d1 = 0;
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

	/*
	!    Find modular inverse (mod 2^96) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q0.d0 & (uint64)1) == 1, "twopmodq96_q8 : (q0.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q1.d0 & (uint64)1) == 1, "twopmodq96_q8 : (q1.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q2.d0 & (uint64)1) == 1, "twopmodq96_q8 : (q2.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q3.d0 & (uint64)1) == 1, "twopmodq96_q8 : (q3.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q4.d0 & (uint64)1) == 1, "twopmodq96_q8 : (q4.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q5.d0 & (uint64)1) == 1, "twopmodq96_q8 : (q5.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q6.d0 & (uint64)1) == 1, "twopmodq96_q8 : (q6.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q7.d0 & (uint64)1) == 1, "twopmodq96_q8 : (q7.d0 & (uint64)1) == 1");
#endif
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d1 = (uint64)0;
	qinv4.d0 = (q4.d0 + q4.d0 + q4.d0) ^ (uint64)2;	qinv4.d1 = (uint64)0;
	qinv5.d0 = (q5.d0 + q5.d0 + q5.d0) ^ (uint64)2;	qinv5.d1 = (uint64)0;
	qinv6.d0 = (q6.d0 + q6.d0 + q6.d0) ^ (uint64)2;	qinv6.d1 = (uint64)0;
	qinv7.d0 = (q7.d0 + q7.d0 + q7.d0) ^ (uint64)2;	qinv7.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*(2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
	defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
	*/
	for(j = 0; j < 4; j++)
	{
		tmp0 = q0.d0*qinv0.d0;
		tmp1 = q1.d0*qinv1.d0;
		tmp2 = q2.d0*qinv2.d0;
		tmp3 = q3.d0*qinv3.d0;
		tmp4 = q4.d0*qinv4.d0;
		tmp5 = q5.d0*qinv5.d0;
		tmp6 = q6.d0*qinv6.d0;
		tmp7 = q7.d0*qinv7.d0;

		qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
		qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
		qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);
		qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);
		qinv4.d0 = qinv4.d0*((uint64)2 - tmp4);
		qinv5.d0 = qinv5.d0*((uint64)2 - tmp5);
		qinv6.d0 = qinv6.d0*((uint64)2 - tmp6);
		qinv7.d0 = qinv7.d0*((uint64)2 - tmp7);
	}

	/* Now that have bottom 64 bits of qinv, do one more Newton iteration
	using full 96-bit operands. See twopmodq96 for details on streamlining here.
	*/
	/* qinv has 96 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + __MULH64(q2.d0, qinv2.d0));
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + __MULH64(q3.d0, qinv3.d0));
	qinv4.d1 = -qinv4.d0*(q4.d1*qinv4.d0 + __MULH64(q4.d0, qinv4.d0));
	qinv5.d1 = -qinv5.d0*(q5.d1*qinv5.d0 + __MULH64(q5.d0, qinv5.d0));
	qinv6.d1 = -qinv6.d0*(q6.d1*qinv6.d0 + __MULH64(q6.d0, qinv6.d0));
	qinv7.d1 = -qinv7.d0*(q7.d1*qinv7.d0 + __MULH64(q7.d0, qinv7.d0));
#else
	MULH64(q0.d0, qinv0.d0, tmp0);
	MULH64(q1.d0, qinv1.d0, tmp1);
	MULH64(q2.d0, qinv2.d0, tmp2);
	MULH64(q3.d0, qinv3.d0, tmp3);
	MULH64(q4.d0, qinv4.d0, tmp4);
	MULH64(q5.d0, qinv5.d0, tmp5);
	MULH64(q6.d0, qinv6.d0, tmp6);
	MULH64(q7.d0, qinv7.d0, tmp7);

	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + tmp2);
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + tmp3);
	qinv4.d1 = -qinv4.d0*(q4.d1*qinv4.d0 + tmp4);
	qinv5.d1 = -qinv5.d0*(q5.d1*qinv5.d0 + tmp5);
	qinv6.d1 = -qinv6.d0*(q6.d1*qinv6.d0 + tmp6);
	qinv7.d1 = -qinv7.d0*(q7.d1*qinv7.d0 + tmp7);
#endif
	qinv0.d1 &= 0x00000000ffffffff;	/* Only want the lower 32 bits here */
	qinv1.d1 &= 0x00000000ffffffff;
	qinv2.d1 &= 0x00000000ffffffff;
	qinv3.d1 &= 0x00000000ffffffff;
	qinv4.d1 &= 0x00000000ffffffff;
	qinv5.d1 &= 0x00000000ffffffff;
	qinv6.d1 &= 0x00000000ffffffff;
	qinv7.d1 &= 0x00000000ffffffff;

	/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT96(qinv0, zshift, lo0);
	LSHIFT96(qinv1, zshift, lo1);
	LSHIFT96(qinv2, zshift, lo2);
	LSHIFT96(qinv3, zshift, lo3);
	LSHIFT96(qinv4, zshift, lo4);
	LSHIFT96(qinv5, zshift, lo5);
	LSHIFT96(qinv6, zshift, lo6);
	LSHIFT96(qinv7, zshift, lo7);

	MULH96_q8(
	  q0, lo0, lo0
	, q1, lo1, lo1
	, q2, lo2, lo2
	, q3, lo3, lo3
	, q4, lo4, lo4
	, q5, lo5, lo5
	, q6, lo6, lo6
	, q7, lo7, lo7);

	/* hi = 0 in this instance, which simplifies things. */
	SUB96(q0, lo0, x0);
	SUB96(q1, lo1, x1);
	SUB96(q2, lo2, x2);
	SUB96(q3, lo3, x3);
	SUB96(q4, lo4, x4);
	SUB96(q5, lo5, x5);
	SUB96(q6, lo6, x6);
	SUB96(q7, lo7, x7);

	if((pshift >> j) & (uint64)1)
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT96(x0, q0), "twopmodq96_q8 : CMPULT96(x0,q0)");
		ASSERT(HERE, CMPULT96(x1, q1), "twopmodq96_q8 : CMPULT96(x1,q1)");
		ASSERT(HERE, CMPULT96(x2, q2), "twopmodq96_q8 : CMPULT96(x2,q2)");
		ASSERT(HERE, CMPULT96(x3, q3), "twopmodq96_q8 : CMPULT96(x3,q3)");
		ASSERT(HERE, CMPULT96(x4, q4), "twopmodq96_q8 : CMPULT96(x4,q4)");
		ASSERT(HERE, CMPULT96(x5, q5), "twopmodq96_q8 : CMPULT96(x5,q5)");
		ASSERT(HERE, CMPULT96(x6, q6), "twopmodq96_q8 : CMPULT96(x6,q6)");
		ASSERT(HERE, CMPULT96(x7, q7), "twopmodq96_q8 : CMPULT96(x7,q7)");
	#endif
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

#if FAC_DEBUG
	if(CMPULT96(q0, x0)){ sprintf(char_buf, "twopmodq96_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint96_base10_char(str0, x0)], &str1[convert_uint96_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q1, x1)){ sprintf(char_buf, "twopmodq96_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint96_base10_char(str0, x1)], &str1[convert_uint96_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q2, x2)){ sprintf(char_buf, "twopmodq96_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint96_base10_char(str0, x2)], &str1[convert_uint96_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q3, x3)){ sprintf(char_buf, "twopmodq96_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint96_base10_char(str0, x3)], &str1[convert_uint96_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q4, x4)){ sprintf(char_buf, "twopmodq96_q8 : (x4 = %s) >= (q4 = %s)", &str0[convert_uint96_base10_char(str0, x4)], &str1[convert_uint96_base10_char(str1, q4)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q5, x5)){ sprintf(char_buf, "twopmodq96_q8 : (x5 = %s) >= (q5 = %s)", &str0[convert_uint96_base10_char(str0, x5)], &str1[convert_uint96_base10_char(str1, q5)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q6, x6)){ sprintf(char_buf, "twopmodq96_q8 : (x6 = %s) >= (q6 = %s)", &str0[convert_uint96_base10_char(str0, x6)], &str1[convert_uint96_base10_char(str1, q6)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q7, x7)){ sprintf(char_buf, "twopmodq96_q8 : (x7 = %s) >= (q7 = %s)", &str0[convert_uint96_base10_char(str0, x7)], &str1[convert_uint96_base10_char(str1, q7)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI96_q8(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3
		, x4, lo4, hi4
		, x5, lo5, hi5
		, x6, lo6, hi6
		, x7, lo7, hi7);

		MULL96_q8(
		  lo0, qinv0, lo0
		, lo1, qinv1, lo1
		, lo2, qinv2, lo2
		, lo3, qinv3, lo3
		, lo4, qinv4, lo4
		, lo5, qinv5, lo5
		, lo6, qinv6, lo6
		, lo7, qinv7, lo7);

		MULH96_q8(
		  q0, lo0, lo0
		, q1, lo1, lo1
		, q2, lo2, lo2
		, q3, lo3, lo3
		, q4, lo4, lo4
		, q5, lo5, lo5
		, q6, lo6, lo6
		, q7, lo7, lo7);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT96(hi0, lo0)) { SUB96(q0, lo0, lo0);	ADD96(lo0, hi0, x0); } else { SUB96(hi0, lo0, x0); }
		if(CMPULT96(hi1, lo1)) { SUB96(q1, lo1, lo1);	ADD96(lo1, hi1, x1); } else { SUB96(hi1, lo1, x1); }
		if(CMPULT96(hi2, lo2)) { SUB96(q2, lo2, lo2);	ADD96(lo2, hi2, x2); } else { SUB96(hi2, lo2, x2); }
		if(CMPULT96(hi3, lo3)) { SUB96(q3, lo3, lo3);	ADD96(lo3, hi3, x3); } else { SUB96(hi3, lo3, x3); }
		if(CMPULT96(hi4, lo4)) { SUB96(q4, lo4, lo4);	ADD96(lo4, hi4, x4); } else { SUB96(hi4, lo4, x4); }
		if(CMPULT96(hi5, lo5)) { SUB96(q5, lo5, lo5);	ADD96(lo5, hi5, x5); } else { SUB96(hi5, lo5, x5); }
		if(CMPULT96(hi6, lo6)) { SUB96(q6, lo6, lo6);	ADD96(lo6, hi6, x6); } else { SUB96(hi6, lo6, x6); }
		if(CMPULT96(hi7, lo7)) { SUB96(q7, lo7, lo7);	ADD96(lo7, hi7, x7); } else { SUB96(hi7, lo7, x7); }

		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT96(x0, q0), "twopmodq96_q8 : CMPULT96(x0,q0)");
			ASSERT(HERE, CMPULT96(x1, q1), "twopmodq96_q8 : CMPULT96(x1,q1)");
			ASSERT(HERE, CMPULT96(x2, q2), "twopmodq96_q8 : CMPULT96(x2,q2)");
			ASSERT(HERE, CMPULT96(x3, q3), "twopmodq96_q8 : CMPULT96(x3,q3)");
			ASSERT(HERE, CMPULT96(x4, q4), "twopmodq96_q8 : CMPULT96(x4,q4)");
			ASSERT(HERE, CMPULT96(x5, q5), "twopmodq96_q8 : CMPULT96(x5,q5)");
			ASSERT(HERE, CMPULT96(x6, q6), "twopmodq96_q8 : CMPULT96(x6,q6)");
			ASSERT(HERE, CMPULT96(x7, q7), "twopmodq96_q8 : CMPULT96(x7,q7)");
		#endif
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
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD96(x0 ,x0, x0);
	ADD96(x1 ,x1, x1);
	ADD96(x2 ,x2, x2);
	ADD96(x3 ,x3, x3);
	ADD96(x4 ,x4, x4);
	ADD96(x5 ,x5, x5);
	ADD96(x6 ,x6, x6);
	ADD96(x7 ,x7, x7);

	SUB96(x0, q0, x0);
	SUB96(x1, q1, x1);
	SUB96(x2, q2, x2);
	SUB96(x3, q3, x3);
	SUB96(x4, q4, x4);
	SUB96(x5, q5, x5);
	SUB96(x6, q6, x6);
	SUB96(x7, q7, x7);

	*checksum2 += x0.d0 + x1.d0 + x2.d0 + x3.d0 + x4.d0 + x5.d0 + x6.d0 + x7.d0;

	r = 0;
	if(CMPEQ96(x0, ONE96)) r +=  1;
	if(CMPEQ96(x1, ONE96)) r +=  2;
	if(CMPEQ96(x2, ONE96)) r +=  4;
	if(CMPEQ96(x3, ONE96)) r +=  8;
	if(CMPEQ96(x4, ONE96)) r += 16;
	if(CMPEQ96(x5, ONE96)) r += 32;
	if(CMPEQ96(x6, ONE96)) r += 64;
	if(CMPEQ96(x7, ONE96)) r +=128;
	return(r);
}
#if 0

You might be able to save a cycle by using CMOVcc instead of bitmasking. Before:

Code:

; eax = hi, ebx = lo, ecx = q
1-2 sub eax, ebx ; hi - lo
2-4 sbb edx, edx ; mask
4-5 and edx, ecx ; q or 0
5-6 add eax, edx ; x

With CMOVcc:

Code:

; eax = hi, ebx = lo, ecx = q
1-2 sub eax, ebx ; hi - lo
2-3 lea edx, [eax+ecx] ; hi - lo + q
3-5 cmovc eax, edx ; x

or

Code:

; eax = hi, ebx = lo, ecx = q
1-2 lea edx, [eax+ecx] ; hi + q
2-3 sub edx, ebx       ; hi + q - lo
1-2 sub eax, ebx       ; hi - lo
3-5 cmovc eax, edx     ; x

The numbers on the left are the cycles when the instruction is expected to start, and when its result is ready to be used, in recent Intel CPUs.

#endif

#if 0

/* Byte offsets for various key pointers in the 2 side-by-side instruction streams: */
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

			movq	%[__q0],%%rsi	/* Use qptr0 as the base address throughout */\n\t@\


	/* If current bit of pshift == 1, double each output modulo q: */	\n\t@\
			/* if((pshift >> j) & (uint64)1) { */	\n\t@\
			movl	%[__pshift],%%eax		\n\t@\
			movl	%[__j],%%ecx			\n\t@\
			shrq	%%cl,%%rax				\n\t@\
			andq	$0x1,%%rax				\n\t@\
		je	twopmodq96_q4_pshiftjmp			\n\t@\
			\n\t@\
		/* if h<l carryout of low 64 bits gives hi=2^32 = 0x100000000, need to zero upper 32 bits prior to double step: */\n\t@\
			movq	$-1,%%rdi	\n\t@\
			shrq	$32,%%rdi	\n\t@\
			andq	%%rdi,%%r12	\n\t@\
			andq	%%rdi,%%r13	\n\t@\
			andq	%%rdi,%%r14	\n\t@\
			andq	%%rdi,%%r15	\n\t@\
		/* x0: */	\n\t@\
			addq	%%r8 ,%%r8 	\n\t@\
			adcq	%%r12,%%r12	\n\t@\
			movq	0x280(%%rsi),%%rax	\n\t@\
			movq	0x288(%%rsi),%%rdx	\n\t@\
			subq	%%rax,%%r8 		\n\t@\
			sbbq	%%rdx,%%r12		\n\t@\
			sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t@\
			andq	%%rdi,%%rax	\n\t@\
			andq	%%rdi,%%rdx	\n\t@\
			addq	%%rax,%%r8 	\n\t@\
			adcq	%%rdx,%%r12	\n\t@\
			\n\t@\
		/* x1: */	\n\t@\
			addq	%%r9 ,%%r9 	\n\t@\
			adcq	%%r13,%%r13	\n\t@\
			movq	0x290(%%rsi),%%rax	\n\t@\
			movq	0x298(%%rsi),%%rdx	\n\t@\
			subq	%%rax,%%r9 		\n\t@\
			sbbq	%%rdx,%%r13		\n\t@\
			sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t@\
			andq	%%rdi,%%rax	\n\t@\
			andq	%%rdi,%%rdx	\n\t@\
			addq	%%rax,%%r9 	\n\t@\
			adcq	%%rdx,%%r13	\n\t@\
			\n\t@\
		/* x2: */	\n\t@\
			addq	%%r10,%%r10	\n\t@\
			adcq	%%r14,%%r14	\n\t@\
			movq	0x2a0(%%rsi),%%rax	\n\t@\
			movq	0x2a8(%%rsi),%%rdx	\n\t@\
			subq	%%rax,%%r10		\n\t@\
			sbbq	%%rdx,%%r14		\n\t@\
			sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t@\
			andq	%%rdi,%%rax	\n\t@\
			andq	%%rdi,%%rdx	\n\t@\
			addq	%%rax,%%r10	\n\t@\
			adcq	%%rdx,%%r14	\n\t@\
			\n\t@\
		/* x3: */	\n\t@\
			addq	%%r11,%%r11	\n\t@\
			adcq	%%r15,%%r15	\n\t@\
			movq	0x2b0(%%rsi),%%rax	\n\t@\
			movq	0x2b8(%%rsi),%%rdx	\n\t@\
			subq	%%rax,%%r11		\n\t@\
			sbbq	%%rdx,%%r15		\n\t@\
			sbbq	%%rdi,%%rdi	/* If 2x-q < 0, re-add q: */\n\t@\
			andq	%%rdi,%%rax	\n\t@\
			andq	%%rdi,%%rdx	\n\t@\
			addq	%%rax,%%r11	\n\t@\
			adcq	%%rdx,%%r15	\n\t@\
		twopmodq96_q4_pshiftjmp:					\n\t@\
			/* } endif((pshift >> j) & (uint64)1) */						\n\t@\
			movq	%%r8 ,0x300(%%rsi)	/* Write lo 64 bits into [__x.d0] */\n\t@\
			movq	%%r12,0x308(%%rsi)	/* Write hi 32 bits into [__x.d1] */\n\t@\
			movq	%%r9 ,0x310(%%rsi)	\n\t@\
			movq	%%r13,0x318(%%rsi)	\n\t@\
			movq	%%r10,0x320(%%rsi)	\n\t@\
			movq	%%r14,0x328(%%rsi)	\n\t@\
			movq	%%r11,0x330(%%rsi)	\n\t@\
			movq	%%r15,0x338(%%rsi)	\n\t@\

Alex`s CADO code snip:

unsigned long tr, t = a[0];
    __asm__ (
      "sub %2, %1\n\t"  /* t -= b ( = a - b) */
      "lea (%1,%3,1), %0\n\t" /* tr = t + m ( = a - b + m) */
      "cmovnc %1, %0\n\t" /* if (a >= b) tr = t */
      : "=&r" (tr), "+&r" (t)
      : "g" (b[0]), "r" (m[0])
      : "cc"


Current code for 96-bit UMULH:

	"/* MULH96_q4((q*, lo*, lo*): Low 64 bits of of lo0-3 in r12-15: */	\n\t"\
			"\n\t"\
		"/* q0 * lo0: */	\n\t"\
			"movq	    (%%rsi),%%rax	/* q0 */\n\t"\
			"mulq	%%r12	/* lo.lo*q.lo in rax:rdx */\n\t"\
			"movl	0x08(%%rsi),%%edi	\n\t"\
			"movq	%%rdx,%%r8 	/* Discard low 64 bits [rax] */\n\t"\
			"movq	%%r12,%%rax	\n\t"\
			"mulq	%%rdi	/* lo.lo*q.hi in rax:rdx */\n\t"\
			"xorq	%%r12,%%r12	\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"movl	0xc8(%%rsi),%%eax	/* Can`t do imulq __lohi,[reg] because that does a sign-extended load of __lohi */\n\t"\
			"imulq	%%rax,%%rdi			/* ...so first load __lohi into low 32 bits of a-reg, then compute lo.hi*q.hi. */\n\t"\
			"addq	%%rdi,%%r12	\n\t"\
			"mulq	    (%%rsi)		/* q.lo*lo.hi in rax:rdx */\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
			"/* Compare upper 64 bits against hi output of squaring, which is also stored in not-yet-right-32-shifted 128-bit form: */\n\t"\
			"movq	0x80(%%rsi),%%rax	/* h.d0 */\n\t"\
			"movq	0x88(%%rsi),%%rdx	/* h.d1 */\n\t"\
			"subq	%%r8 ,%%rax	/* rax = (h-l).d0 */\n\t"\
			"sbbq	%%r12,%%rdx	/* rdx = (h-l).d1 */\n\t"\
			"/* If there was a borrow out of the hi-word subtract (h-l), use it to create a bitmask of all ones in the free rdi-register: */\n\t"\
			"sbbq	%%rdi,%%rdi	\n\t"\
			"/* Right-shift the 2-word result 32 bits: */\n\t"\
			"shrdq	$32,%%rdx,%%rax	/* Only keep high 96-bit output */\n\t"\
			"/* Initial impulse here is to use arithmetic rshift (sarq) of hi word, since if h<l, carryout of low 64 bits -> hi=2^32 = 0x100000000, */\n\t"\
			"/* which would need to be explicitly zeroed prior to the conditional doubling step. But arith-shift needs the hi bit of rdx to reflect */\n\t"\
			"/* any carryout, which will not be the case if (h-l) has 96 sig. bits (i.e. 128 bits in the present << 32 form. To properly handle that, */\n\t"\
			"/* need either another 2-word shift `shrdq $32,%%rdi,%%rdx` here or a later explicit zeroing of the hi 32 bits ... prefer the latter. */\n\t"\
			"shrq	$32,%%rdx			\n\t"\
			"/* Apply mask to lo/hi words of q and add result to (h-l): */\n\t"\
			"movq	    (%%rsi),%%r8 	\n\t"\
			"movq	0x08(%%rsi),%%r12	\n\t"\
			"andq	%%rdi,%%r8 	\n\t"\
			"andq	%%rdi,%%r12	\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"\n\t"\

Try speeding the UMULHs in the 96-bit modmul by modifying the 2-word representation to

[2] x = a+b*2^32, y = c+d*2^32, a/c 32-bit, b/d 64-bit:

x*y = a*b[64-bit] + (a*d+b*c)[96-97-bit]*2^32 + b*d[128-bit]*2^64, extracting hi 96: a*b[none]+(a*d+b*c)[hi32-33]+b*d[hi96]

Notes:
- Don`t need to multiply the low 32-bit words at all;
- Can do 2 32-bit umulh to get the hi32 bits of a*d and b*c, then add those in a 64-bit register;
- Need 1 full-length 128-bit MUL to get b*d, of which we discard the lo32 bits via right-shft, then add-with-carry the 33-bit partial sum (a*d + b*c).hi .

	"/* MULH96_q4((q*, lo*, lo*): Low 64 bits of of lo0-3 in r12-15: */	\n\t"\
			"\n\t"\
		"/* q0 * lo0: */	\n\t"\
			"movq	    (%%rsi),%%rax	/* q0 */\n\t"\
			"shrq	$32,%%rax	/* Discard lo32 bits of q, keep mid32 */\n\t"\
			"movq	%%r12,%%rdx	\n\t"\
			"mull	%%edx	/* lo.lo32*q.mid32 in eax:edx */\n\t"\
			"movq	%%rdx,%%r8 	/* Store hi32 bits of 64-bit product in lower half of a 64-bit register */\n\t"\
			"movq	    (%%rsi),%%rax	/* q0 */\n\t"\
			"movq	%%r12,%%rdx	\n\t"\
			"shrq	$32,%%rdx	/* Discard lo32 bits of x, keep mid32 */\n\t"\
			"mull	%%edx	/* q.lo32*lo.mid32 in eax:edx */\n\t"\
			"addq	%%rdx,%%r8 	/* Add 32-bit (a*d + b*c).hi32 terms in a 64-bit register */\n\t"\
			"movq	        (%%rsi),%%rax	/* q0.lo64 */\n\t"\
			"shrdq	$32,0x08(%%rsi),%%rax	/* Rightward shift-in hi32 bits to get q0.hi64 */\n\t"\
			"shrdq	$32,0xc8(%%rsi),%%r12	/* Rightward shift-in hi32 bits to get lo.hi64 */ */\n\t"\
			"mulq	%%r12	/* q.hi64*lo.hi64 in rax:rdx */\n\t"\
			"shrdq	$32,%%rdx,%%rax	/* Rightward shift 128-bit hi-part product to keep just the hi96 bits */\n\t"\
			"xorq	%%r12,%%r12	/* zero the r12 register which will hold the hi32 bits of the UMULH result */\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"shrlq	$32,%%r8 ,%%r12	/* Lefttward shift 96-bit UMULH result to maintain compatibility with existing code sequence below */\n\t"\

Compare cost of mul-sequence above:

Old: 4 mul (1 imulq, 3 mulq), 5 add (3 addq, 2 adcq), 0 shift, 1 xor, 5 mov
New: 3 mul (2 mull, 1 mulq), 3 add (2 addq, 1 adcq), 6 shift (2 single-reg, 4 double-reg), 1 xor, 6 mov

New scheme has 1 fewer mul, 2 fewer add, but 1 more mov and a whopping 6 more shift!

			"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
			"/* Compare upper 64 bits against hi output of squaring, which is also stored in not-yet-right-32-shifted 128-bit form: */\n\t"\
			"movq	0x80(%%rsi),%%rax	/* h.d0 */\n\t"\
			"movq	0x88(%%rsi),%%rdx	/* h.d1 */\n\t"\
			"subq	%%r8 ,%%rax	/* rax = (h-l).d0 */\n\t"\
			"sbbq	%%r12,%%rdx	/* rdx = (h-l).d1 */\n\t"\
			"/* If there was a borrow out of the hi-word subtract (h-l), use it to create a bitmask of all ones in the free rdi-register: */\n\t"\
			"sbbq	%%rdi,%%rdi	\n\t"\
			"/* Right-shift the 2-word result 32 bits: */\n\t"\
			"shrdq	$32,%%rdx,%%rax	/* Only keep high 96-bit output */\n\t"\
			"/* Initial impulse here is to use arithmetic rshift (sarq) of hi word, since if h<l, carryout of low 64 bits -> hi=2^32 = 0x100000000, */\n\t"\
			"/* which would need to be explicitly zeroed prior to the conditional doubling step. But arith-shift needs the hi bit of rdx to reflect */\n\t"\
			"/* any carryout, which will not be the case if (h-l) has 96 sig. bits (i.e. 128 bits in the present << 32 form. To properly handle that, */\n\t"\
			"/* need either another 2-word shift `shrdq $32,%%rdi,%%rdx` here or a later explicit zeroing of the hi 32 bits ... prefer the latter. */\n\t"\
			"shrq	$32,%%rdx			\n\t"\
			"/* Apply mask to lo/hi words of q and add result to (h-l): */\n\t"\
			"movq	    (%%rsi),%%r8 	\n\t"\
			"movq	0x08(%%rsi),%%r12	\n\t"\
			"andq	%%rdi,%%r8 	\n\t"\
			"andq	%%rdi,%%r12	\n\t"\
			"addq	%%rax,%%r8 	\n\t"\
			"adcq	%%rdx,%%r12	\n\t"\
			"\n\t"\

#endif

