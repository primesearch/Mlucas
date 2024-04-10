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

#if FAC_DEBUG
	char char_buf[1024], str0[64], str1[64];
#endif

/***********************************************************************************/
/*** 96-BIT INPUTS, 128-BIT MULL OPERATION *****************************************/
/***********************************************************************************/
/*
Function to find 2^(-p) mod q, where p and q are 64 and 128-bit unsigned integers,
respectively, and at most the least-significant 96 bits of q are nonzero.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^128 (i.e. our MODQ
operation effects multiply modulo 2^128).

The key 3-operation sequence here is as follows:

	SQR_LOHI128_96(x,lo,hi);// Input   x has 96 bits; output lo has 128 bits; hi has 64 bits
	MULL128(lo,qinv,lo);	// Inputs lo & qinv, and output (overwrites lo) have 128 bits
	MULH128(q,lo,lo);	// Inputs  q &   lo, and output (overwrites lo) have 128 bits
							//    (but only lower 96 bits of q and output are nonzero).
*/
uint64 twopmodq128_96(uint64 p, uint64 k)
{
#if FAC_DEBUG
	int dbg = (p == 0);
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 hi64;
	uint128 q, qinv, x, lo;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

#if FAC_DEBUG
if(dbg)printf("twopmodq128_96:\n");
#endif

	ASSERT((p >> 63) == 0, "p must be < 2^63!");
	q.d0 = p+p;	q.d1 = 0;
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(q.d0, k,&q.d0,&q.d1);
#else
	MUL_LOHI64(q.d0, k, q.d0, q.d1);
#endif
	q.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	ASSERT((q.d1 >> 32) == 0, "(q.d1 >> 32) != 0");

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 128;
		/*
		!    find number of leading zeros in p, use it to find the position of the leftmost
		!    ones bit, and subtract 7 to account for the fact that we can do the powering for the leftmost
		!    7 bits via a simple shift.
		*/
	#if 0
		start_index = 64-leadz64(pshift)-7;	/* Leftward bit at which to start the l-r binary powering, assuming
							 the leftmost 7 bits have already been processed via a shift (see next). */
		zshift = 127-ibits64(pshift,start_index,7);	/* Since 7 bits with leftmost bit = 1 is guaranteed
								to be in [64,127], the shift count here is in [0, 63].
								That means that zstart < 2^64. Together with the fact that
								squaring a power of two gives another power of two, we can
								simplify the modmul code sequence for the first iteration.
								Every little bit counts (literally in this case :), right? */
	#else
		j = leadz64(pshift);
		start_index =  64-j-7;
		zshift = 127 - ((pshift<<j) >> 57);
	#endif
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

#if FAC_DEBUG
	if(dbg)	printf("start_index = %u\n", (uint32)start_index);
#endif
	/*
	!    Find modular inverse (mod 2^128) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT((q.d0 & (uint64)1) == 1, "twopmodq128_96 : (q.d0 & (uint64)1) == 1");
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

	/* qinv.d1 = 0 and q*qinv will == 1 mod 2^64 here; can take advantage of that fact to speed the iteration
	(i.e. q*qinv = (q.d1*2^64 + q.d0)*qinv.d0 = q.d1*qinv.d0*2^64 + q.d0*qinv.d0 == 1 mod 2^64,
	i.e. do a MULL64(q.d1,qinv.d0) + MULH64(q.d0,qinv.d0) to get bits 64-127 of q*qinv = r, set lower half = 1,
	then simply negate upper half to get 2-q*qinv (mod 2^128), then use that lower half still = 1 to
	simplify the qinv*r MULL128 operation, i.e. qinv*r = qinv.d0*(r.d1*2^64 + 1) == (MULL64(qinv.d0*r.d1)*2^64 + qinv.d0) mod 2^128.
	This needs a total of just 3 MUL instructions and 2 ALUs, compared to 8 MULs and 8 ALUs for the original sequence.
	*/
#if FAC_DEBUG
	MULL128(q, qinv, x);
	SUB128 (TWO128, x, x);
	MULL128(qinv, x, x);
#endif
	/* qinv has 128 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + __MULH64(q.d0, qinv.d0));
#else
	MULH64(q.d0, qinv.d0, hi64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + hi64);
#endif

#if FAC_DEBUG
	ASSERT(qinv.d1 == x.d1 && qinv.d0 == x.d0, "twopmodq128_96 : qinv.d1 == x.d1 && qinv.d0 == x.d0");
#endif

#if FAC_DEBUG
	if(dbg) printf("q    = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q   )]);
	if(dbg) printf("qinv = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv)]);
#endif

	/* Since zstart is a power of two < 2^128, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

#if FAC_DEBUG
	if(dbg) printf("zshift  = %u\n", zshift);
#endif

	/* MULL128(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT128(qinv, zshift, lo);

#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	MULH128(q,lo,lo);

#if FAC_DEBUG
	if(dbg) printf("q*lo/2^128 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	/* hi64 = 0 in this instance, which simplifies things. */
	SUB128(q, lo, x);

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

	if((pshift >> j) & (uint64)1)
	{
	#if FAC_DEBUG
		ASSERT(CMPULT128(x,q), "twopmodq128_96 : CMPULT128(x,q)");
	#endif
		ADD128(x,x,x);
		if(CMPULE128(q,x)) SUB128(x,q,x);
	}

#if FAC_DEBUG
	if(CMPULT128(q, x)){ sprintf(char_buf, "twopmodq128_96 : (x0 = %s) >= (q = %s)", &str0[convert_uint128_base10_char(str0, x)], &str1[convert_uint128_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI128_96(x,lo,hi64);	/* lo has 128 bits, hi64 has 64. */
		MULL128(lo,qinv,lo);		/* lo has 128 bits */
		MULH128(q,lo,lo);		/* On output, lo has  96 bits */

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(lo.d1 != 0 || hi64 < lo.d0)
		{
			SUB128(q, lo, x);
			x.d0 +=  hi64;
			x.d1 += (x.d0 < hi64);	/* Had a carry. */
		}
		else
		{
			x.d0 =  hi64 - lo.d0;
			x.d1 = (uint64)0;
		}

#if FAC_DEBUG
if(dbg)printf("j = %2d, x = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(CMPULT128(x,q), "twopmodq128_96 : CMPULT128(x,q)");
		#endif
			ADD128(x,x,x);	/* Since we're using 128-bit arithmetic for the add, x+x cannot overflow. */
			if(CMPULE128(q,x)) SUB128(x,q,x);
#if FAC_DEBUG
	if(dbg) printf("*2= %s", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
		}
#if FAC_DEBUG
if(dbg)printf("\n");
#endif
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD128(x,x,x);
	//*** Aug 2022: see note at end of twopmodq96() for why we can add-sans-carry the constant (0 or 2) here:
	x.d0 += FERMAT;
	SUB128(x,q,x);

#if FAC_DEBUG
if(dbg)printf("x0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
	return (uint64)CMPEQ128(x, ONE128) ;
}

/*** 4-trial-factor version ***/
uint64 twopmodq128_96_q4(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
#if FAC_DEBUG
	int dbg = (p == 0);
	uint128 x,y,lo,hi;
	uint64 hi64;
#endif
	 int32 j;
	uint64 hi0, hi1, hi2, hi3, r;
	uint128 q0, q1, q2, q3
		, qinv0, qinv1, qinv2, qinv3
		, x0, x1, x2, x3
		, lo0, lo1, lo2, lo3;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

#if FAC_DEBUG
if(dbg)printf("twopmodq128_96_q4:\n");
#endif
	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 128;
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 7 to account for the fact that we can do the powering for the leftmost
	!    7 bits via a simple shift.
	*/
		start_index = 64-leadz64(pshift)-7;	/* Leftward bit at which to start the l-r binary powering, assuming
							 the leftmost 7 bits have already been processed via a shift (see next). */
		zshift = 127-ibits64(pshift,start_index,7);	/* Since 7 bits with leftmost bit = 1 is guaranteed
								to be in [64,127], the shift count here is in [0, 63].
								That means that zstart < 2^64. Together with the fact that
								squaring a power of two gives another power of two, we can
								simplify the modmul code sequence for the first iteration.
								Every little bit counts (literally in this case :), right? */
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

#if FAC_DEBUG
	if(dbg)	printf("start_index = %u\n", (uint32)start_index);
#endif

	ASSERT((p >> 63) == 0, "p must be < 2^63!");
	q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;
	q0.d1 = q1.d1 = q2.d1 = q3.d1 = 0;
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
	ASSERT((q0.d1 >> 32) == 0, "(q0.d1 >> 32) != 0");
	ASSERT((q1.d1 >> 32) == 0, "(q1.d1 >> 32) != 0");
	ASSERT((q2.d1 >> 32) == 0, "(q2.d1 >> 32) != 0");
	ASSERT((q3.d1 >> 32) == 0, "(q3.d1 >> 32) != 0");

	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	q2.d0 += 1;
	q3.d0 += 1;

	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT((q0.d0 & (uint64)1) == 1, "twopmodq128_96_q4 : (q0.d0 & (uint64)1) == 1");
	ASSERT((q1.d0 & (uint64)1) == 1, "twopmodq128_96_q4 : (q1.d0 & (uint64)1) == 1");
	ASSERT((q2.d0 & (uint64)1) == 1, "twopmodq128_96_q4 : (q2.d0 & (uint64)1) == 1");
	ASSERT((q3.d0 & (uint64)1) == 1, "twopmodq128_96_q4 : (q3.d0 & (uint64)1) == 1");
#endif
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d1 = (uint64)0;

	for(j = 0; j < 4; j++)
	{
		hi0 = q0.d0*qinv0.d0;
		hi1 = q1.d0*qinv1.d0;
		hi2 = q2.d0*qinv2.d0;
		hi3 = q3.d0*qinv3.d0;

		qinv0.d0 = qinv0.d0*((uint64)2 - hi0);
		qinv1.d0 = qinv1.d0*((uint64)2 - hi1);
		qinv2.d0 = qinv2.d0*((uint64)2 - hi2);
		qinv3.d0 = qinv3.d0*((uint64)2 - hi3);
	}

	/* qinv has 128 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + __MULH64(q2.d0, qinv2.d0));
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + __MULH64(q3.d0, qinv3.d0));
#else
	MULH64(q0.d0, qinv0.d0, hi0);
	MULH64(q1.d0, qinv1.d0, hi1);
	MULH64(q2.d0, qinv2.d0, hi2);
	MULH64(q3.d0, qinv3.d0, hi3);

	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + hi0);
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + hi1);
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + hi2);
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + hi3);
#endif

#if FAC_DEBUG
	if(dbg) printf("q    = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q0   )]);
	if(dbg) printf("qinv = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv0)]);
#endif

	/* Since zstart is a power of two < 2^128, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

#if FAC_DEBUG
	if(dbg) printf("zshift  = %u\n", zshift);
#endif

	/* MULL128(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT128(qinv0, zshift, lo0);
	LSHIFT128(qinv1, zshift, lo1);
	LSHIFT128(qinv2, zshift, lo2);
	LSHIFT128(qinv3, zshift, lo3);

#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo0)]);
#endif

#if(USE_128x96 > 0)
	/* Need to be careful about the order of the 2 inputs here,
	as MULH128x96 assumes the 2nd input is the one which is < 2^96: */
	MULH128x96_q4(
	  lo0, q0, lo0
	, lo1, q1, lo1
	, lo2, q2, lo2
	, lo3, q3, lo3);
#else
	MULH128_q4(
	  q0, lo0, lo0
	, q1, lo1, lo1
	, q2, lo2, lo2
	, q3, lo3, lo3);
#endif

#if FAC_DEBUG
	if(dbg) printf("q*lo/2^128 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo0)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB128(q0, lo0, x0);
	SUB128(q1, lo1, x1);
	SUB128(q2, lo2, x2);
	SUB128(q3, lo3, x3);

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x0)]);
#endif

	if((pshift >> j) & (uint64)1)
	{
	#if FAC_DEBUG
		ASSERT(CMPULT128(x0, q0), "twopmodq128_96_q4 : CMPULT128(x0, q0)");
		ASSERT(CMPULT128(x1, q1), "twopmodq128_96_q4 : CMPULT128(x1, q1)");
		ASSERT(CMPULT128(x2, q2), "twopmodq128_96_q4 : CMPULT128(x2, q2)");
		ASSERT(CMPULT128(x3, q3), "twopmodq128_96_q4 : CMPULT128(x3, q3)");
	#endif
		ADD128(x0, x0, x0);
		ADD128(x1, x1, x1);
		ADD128(x2, x2, x2);
		ADD128(x3, x3, x3);

		if(CMPULE128(q0, x0)) SUB128(x0, q0, x0);
		if(CMPULE128(q1, x1)) SUB128(x1, q1, x1);
		if(CMPULE128(q2, x2)) SUB128(x2, q2, x2);
		if(CMPULE128(q3, x3)) SUB128(x3, q3, x3);
	}

#if FAC_DEBUG
	if(CMPULT128(q0, x0)){ sprintf(char_buf, "twopmodq128_96_q4 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint128_base10_char(str0, x0)], &str1[convert_uint128_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q1, x1)){ sprintf(char_buf, "twopmodq128_96_q4 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint128_base10_char(str0, x1)], &str1[convert_uint128_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q2, x2)){ sprintf(char_buf, "twopmodq128_96_q4 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint128_base10_char(str0, x2)], &str1[convert_uint128_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q3, x3)){ sprintf(char_buf, "twopmodq128_96_q4 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint128_base10_char(str0, x3)], &str1[convert_uint128_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x0)]);
#endif

	for(j = start_index-2; j >= 0; j--)
	{
#if FAC_DEBUG
	x=x0;
	SQR_LOHI128(x,lo,hi);
#endif
		/* Haven't gotten IA64 version of this working properly yet:
		SQR_LOHI_INPLACE128_96_q4(
		  x0, hi0
		, x1, hi1
		, x2, hi2
		, x3, hi3);
		*/
		SQR_LOHI128_96_q4(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3);

		/* For unknown reasons, the 8-operand version of MULL128 was slower than one-at-a-time. */
#if FAC_DEBUG
	ASSERT(CMPEQ128(lo,lo0), "twopmodq128_96_q4 : CMPEQ128(SQR_LO)");
	ASSERT(hi.d1 == 0      , "twopmodq128_96_q4 : hi.d1 != 0");
	hi64 = hi.d0;
	ASSERT(hi64 == hi0     , "twopmodq128_96_q4 : CMPEQ128(SQR_HI)");
	x=lo0;y=qinv0;
	MULL128(x,y,lo);
#endif
		MULL128_INPLACE_q4(
		  lo0, qinv0
		, lo1, qinv1
		, lo2, qinv2
		, lo3, qinv3);
	#if(USE_128x96 > 0)
	/* Need to be careful about the order of the 2 inputs here,
	as MULH128x96 assumes the 2nd input is the one which is < 2^96: */
#if FAC_DEBUG
	ASSERT(CMPEQ128(lo,lo0), "twopmodq128_96_q4 : CMPEQ128(MULL128)");
	x=lo0;y=q0;
	MULH128(x,y,lo);
	MULH128(y,x,hi);
	ASSERT(CMPEQ128(lo, hi), "twopmodq128_96_q4 : MULH(X,Y) != MULH(Y,X)");
#endif
		MULH128x96_q4(
		  lo0, q0, lo0
		, lo1, q1, lo1
		, lo2, q2, lo2
		, lo3, q3, lo3);
	#else
		MULH128_q4(
		  q0, lo0, lo0
		, q1, lo1, lo1
		, q2, lo2, lo2
		, q3, lo3, lo3);
	#endif
#if FAC_DEBUG
	ASSERT(CMPEQ128(lo,lo0), "twopmodq128_96_q4 : CMPEQ128(MULH())");
	/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
	if(lo.d1 != 0 || hi64 < lo.d0)
	{
		SUB128(q0, lo, x);
		x.d0 +=  hi64;
		x.d1 += (x.d0 < hi64);	/* Had a carry. */
	}
	else
	{
		x.d0 =  hi64 - lo.d0;
		x.d1 = (uint64)0;
	}
#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(lo0.d1 != 0 || hi0 < lo0.d0){ SUB128(q0, lo0, x0);	x0.d0 +=  hi0; x0.d1 += (x0.d0 < hi0); } else { x0.d0 =  hi0 - lo0.d0; x0.d1 = (uint64)0; }
		if(lo1.d1 != 0 || hi1 < lo1.d0){ SUB128(q1, lo1, x1);	x1.d0 +=  hi1; x1.d1 += (x1.d0 < hi1); } else { x1.d0 =  hi1 - lo1.d0; x1.d1 = (uint64)0; }
		if(lo2.d1 != 0 || hi2 < lo2.d0){ SUB128(q2, lo2, x2);	x2.d0 +=  hi2; x2.d1 += (x2.d0 < hi2); } else { x2.d0 =  hi2 - lo2.d0; x2.d1 = (uint64)0; }
		if(lo3.d1 != 0 || hi3 < lo3.d0){ SUB128(q3, lo3, x3);	x3.d0 +=  hi3; x3.d1 += (x3.d0 < hi3); } else { x3.d0 =  hi3 - lo3.d0; x3.d1 = (uint64)0; }
#if FAC_DEBUG
	ASSERT(CMPEQ128( x, x0), "twopmodq128_96_q4 : CMPEQ128(MULH())");
#endif

#if FAC_DEBUG
if(dbg)printf("j = %2d, x = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x0)]);
#endif

		if((pshift >> j) & (uint64)1)
		{
#if FAC_DEBUG
	ADD128(x,x,x);	/* Since we're using 128-bit arithmetic for the add, x+x cannot overflow. */
	if(CMPULE128(q0,x)) SUB128(x,q0,x);
#endif
		#if FAC_DEBUG
			ASSERT(CMPULT128(x0, q0), "twopmodq128_96_q4 : CMPULT128(x0, q0)");
			ASSERT(CMPULT128(x1, q1), "twopmodq128_96_q4 : CMPULT128(x1, q1)");
			ASSERT(CMPULT128(x2, q2), "twopmodq128_96_q4 : CMPULT128(x2, q2)");
			ASSERT(CMPULT128(x3, q3), "twopmodq128_96_q4 : CMPULT128(x3, q3)");
		#endif
			ADD128(x0, x0, x0);
			ADD128(x1, x1, x1);
			ADD128(x2, x2, x2);
			ADD128(x3, x3, x3);

			if(CMPULE128(q0, x0)) SUB128(x0, q0, x0);
			if(CMPULE128(q1, x1)) SUB128(x1, q1, x1);
			if(CMPULE128(q2, x2)) SUB128(x2, q2, x2);
			if(CMPULE128(q3, x3)) SUB128(x3, q3, x3);
#if FAC_DEBUG
	if(dbg) printf("*2= %s", &char_buf[convert_uint128_base10_char(char_buf, x0)]);
#endif
		}
#if FAC_DEBUG
if(dbg)printf("\n");
	ASSERT(CMPEQ128( x, x0), "twopmodq128_96_q4 : CMPEQ128(MULH())");
#endif
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	ADD128(x0 ,x0, x0);
	ADD128(x1 ,x1, x1);
	ADD128(x2 ,x2, x2);
	ADD128(x3 ,x3, x3);
	//*** Aug 2022: see note at end of twopmodq96() for why we can add-sans-carry the constant (0 or 2) here:
	x0.d0 += FERMAT;
	x1.d0 += FERMAT;
	x2.d0 += FERMAT;
	x3.d0 += FERMAT;
	SUB128(x0, q0, x0);
	SUB128(x1, q1, x1);
	SUB128(x2, q2, x2);
	SUB128(x3, q3, x3);

#if FAC_DEBUG
if(dbg)printf("x0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x0)]);
#endif

	/* Only do the full 128-bit (Xj== 1) check if the bottom 64 bits of Xj == 1: */
	r = 0;
	if(x0.d0 == 1) r += ((uint64)CMPEQ128(x0, ONE128) << 0);
	if(x1.d0 == 1) r += ((uint64)CMPEQ128(x1, ONE128) << 1);
	if(x2.d0 == 1) r += ((uint64)CMPEQ128(x2, ONE128) << 2);
	if(x3.d0 == 1) r += ((uint64)CMPEQ128(x3, ONE128) << 3);
	return(r);
}

/*** 8-trial-factor version ***/
uint64 twopmodq128_96_q8(uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint64_base10_char(char_buf, p)], "0");
#endif
	 int32 j;
	uint64 hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7, r;
	uint128 q0, q1, q2, q3, q4, q5, q6, q7
		, qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	uint32 FERMAT = isPow2_64(p)<<1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

#if FAC_DEBUG
if(dbg)printf("twopmodq128_96_q8:\n");
#endif
	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 128;
		/*
		!    find number of leading zeros in p, use it to find the position of the leftmost
		!    ones bit, and subtract 7 to account for the fact that we can do the powering for the leftmost
		!    7 bits via a simple shift.
		*/
	#if 0
		start_index = 64-leadz64(pshift)-7;	/* Leftward bit at which to start the l-r binary powering, assuming
							 the leftmost 7 bits have already been processed via a shift (see next). */
		zshift = 127-ibits64(pshift,start_index,7);	/* Since 7 bits with leftmost bit = 1 is guaranteed
								to be in [64,127], the shift count here is in [0, 63].
								That means that zstart < 2^64. Together with the fact that
								squaring a power of two gives another power of two, we can
								simplify the modmul code sequence for the first iteration.
								Every little bit counts (literally in this case :), right? */
	#else
		j = leadz64(pshift);
		start_index =  64-j-7;
		zshift = 127 - ((pshift<<j) >> 57);
	#endif
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

	ASSERT((p >> 63) == 0, "p must be < 2^63!");
	q0.d0 = q1.d0 = q2.d0 = q3.d0 = q4.d0 = q5.d0 = q6.d0 = q7.d0 = p+p;
	q0.d1 = q1.d1 = q2.d1 = q3.d1 = q4.d1 = q5.d1 = q6.d1 = q7.d1 = 0;
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(q0.d0, k0,&q0.d0,&q0.d1);
	MUL_LOHI64(q1.d0, k1,&q1.d0,&q1.d1);
	MUL_LOHI64(q2.d0, k2,&q2.d0,&q2.d1);
	MUL_LOHI64(q3.d0, k3,&q3.d0,&q3.d1);
	MUL_LOHI64(q4.d0, k4,&q4.d0,&q4.d1);
	MUL_LOHI64(q5.d0, k5,&q5.d0,&q5.d1);
	MUL_LOHI64(q6.d0, k6,&q6.d0,&q6.d1);
	MUL_LOHI64(q7.d0, k7,&q7.d0,&q7.d1);
#else
	MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
	MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
	MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
	MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);
	MUL_LOHI64(q4.d0, k4, q4.d0, q4.d1);
	MUL_LOHI64(q5.d0, k5, q5.d0, q5.d1);
	MUL_LOHI64(q6.d0, k6, q6.d0, q6.d1);
	MUL_LOHI64(q7.d0, k7, q7.d0, q7.d1);
#endif

	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	q2.d0 += 1;
	q3.d0 += 1;
	q4.d0 += 1;
	q5.d0 += 1;
	q6.d0 += 1;
	q7.d0 += 1;
	ASSERT((q0.d1 >> 32) == 0, "(q0.d1 >> 32) != 0");
	ASSERT((q1.d1 >> 32) == 0, "(q1.d1 >> 32) != 0");
	ASSERT((q2.d1 >> 32) == 0, "(q2.d1 >> 32) != 0");
	ASSERT((q3.d1 >> 32) == 0, "(q3.d1 >> 32) != 0");
	ASSERT((q4.d1 >> 32) == 0, "(q4.d1 >> 32) != 0");
	ASSERT((q5.d1 >> 32) == 0, "(q5.d1 >> 32) != 0");
	ASSERT((q6.d1 >> 32) == 0, "(q6.d1 >> 32) != 0");
	ASSERT((q7.d1 >> 32) == 0, "(q7.d1 >> 32) != 0");

	/*
	!    Find modular inverse (mod 2^128) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT((q0.d0 & (uint64)1) == 1, "twopmodq128_96_q8 : (q0.d0 & (uint64)1) == 1");
	ASSERT((q1.d0 & (uint64)1) == 1, "twopmodq128_96_q8 : (q1.d0 & (uint64)1) == 1");
	ASSERT((q2.d0 & (uint64)1) == 1, "twopmodq128_96_q8 : (q2.d0 & (uint64)1) == 1");
	ASSERT((q3.d0 & (uint64)1) == 1, "twopmodq128_96_q8 : (q3.d0 & (uint64)1) == 1");
	ASSERT((q4.d0 & (uint64)1) == 1, "twopmodq128_96_q8 : (q4.d0 & (uint64)1) == 1");
	ASSERT((q5.d0 & (uint64)1) == 1, "twopmodq128_96_q8 : (q5.d0 & (uint64)1) == 1");
	ASSERT((q6.d0 & (uint64)1) == 1, "twopmodq128_96_q8 : (q6.d0 & (uint64)1) == 1");
	ASSERT((q7.d0 & (uint64)1) == 1, "twopmodq128_96_q8 : (q7.d0 & (uint64)1) == 1");
#endif
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d1 = (uint64)0;
	qinv4.d0 = (q4.d0 + q4.d0 + q4.d0) ^ (uint64)2;	qinv4.d1 = (uint64)0;
	qinv5.d0 = (q5.d0 + q5.d0 + q5.d0) ^ (uint64)2;	qinv5.d1 = (uint64)0;
	qinv6.d0 = (q6.d0 + q6.d0 + q6.d0) ^ (uint64)2;	qinv6.d1 = (uint64)0;
	qinv7.d0 = (q7.d0 + q7.d0 + q7.d0) ^ (uint64)2;	qinv7.d1 = (uint64)0;

	for(j = 0; j < 4; j++)
	{
		hi0 = q0.d0*qinv0.d0;
		hi1 = q1.d0*qinv1.d0;
		hi2 = q2.d0*qinv2.d0;
		hi3 = q3.d0*qinv3.d0;
		hi4 = q4.d0*qinv4.d0;
		hi5 = q5.d0*qinv5.d0;
		hi6 = q6.d0*qinv6.d0;
		hi7 = q7.d0*qinv7.d0;

		qinv0.d0 = qinv0.d0*((uint64)2 - hi0);
		qinv1.d0 = qinv1.d0*((uint64)2 - hi1);
		qinv2.d0 = qinv2.d0*((uint64)2 - hi2);
		qinv3.d0 = qinv3.d0*((uint64)2 - hi3);
		qinv4.d0 = qinv4.d0*((uint64)2 - hi4);
		qinv5.d0 = qinv5.d0*((uint64)2 - hi5);
		qinv6.d0 = qinv6.d0*((uint64)2 - hi6);
		qinv7.d0 = qinv7.d0*((uint64)2 - hi7);
	}

	/* qinv has 128 bits, but only the upper 64 get modified here. */
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
	MULH64(q0.d0, qinv0.d0, hi0);
	MULH64(q1.d0, qinv1.d0, hi1);
	MULH64(q2.d0, qinv2.d0, hi2);
	MULH64(q3.d0, qinv3.d0, hi3);
	MULH64(q4.d0, qinv4.d0, hi4);
	MULH64(q5.d0, qinv5.d0, hi5);
	MULH64(q6.d0, qinv6.d0, hi6);
	MULH64(q7.d0, qinv7.d0, hi7);

	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + hi0);
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + hi1);
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + hi2);
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + hi3);
	qinv4.d1 = -qinv4.d0*(q4.d1*qinv4.d0 + hi4);
	qinv5.d1 = -qinv5.d0*(q5.d1*qinv5.d0 + hi5);
	qinv6.d1 = -qinv6.d0*(q6.d1*qinv6.d0 + hi6);
	qinv7.d1 = -qinv7.d0*(q7.d1*qinv7.d0 + hi7);
#endif

	/* Since zstart is a power of two < 2^128, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL128(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT128(qinv0, zshift, lo0);
	LSHIFT128(qinv1, zshift, lo1);
	LSHIFT128(qinv2, zshift, lo2);
	LSHIFT128(qinv3, zshift, lo3);
	LSHIFT128(qinv4, zshift, lo4);
	LSHIFT128(qinv5, zshift, lo5);
	LSHIFT128(qinv6, zshift, lo6);
	LSHIFT128(qinv7, zshift, lo7);

#if(USE_128x96 > 0)
	/* Need to be careful about the order of the 2 inputs here,
	as MULH128x96 assumes the 2nd input is the one which is < 2^96: */
	MULH128x96_q8(
	  lo0, q0, lo0
	, lo1, q1, lo1
	, lo2, q2, lo2
	, lo3, q3, lo3
	, lo4, q4, lo4
	, lo5, q5, lo5
	, lo6, q6, lo6
	, lo7, q7, lo7);
#else
	MULH128_q8(
	  q0, lo0, lo0
	, q1, lo1, lo1
	, q2, lo2, lo2
	, q3, lo3, lo3
	, q4, lo4, lo4
	, q5, lo5, lo5
	, q6, lo6, lo6
	, q7, lo7, lo7);
#endif
	/* hi = 0 in this instance, which simplifies things. */
	SUB128(q0, lo0, x0);
	SUB128(q1, lo1, x1);
	SUB128(q2, lo2, x2);
	SUB128(q3, lo3, x3);
	SUB128(q4, lo4, x4);
	SUB128(q5, lo5, x5);
	SUB128(q6, lo6, x6);
	SUB128(q7, lo7, x7);

	if((pshift >> j) & (uint64)1)
	{
	#if FAC_DEBUG
		ASSERT(CMPULT128(x0, q0), "twopmodq128_96_q8 : CMPULT128(x0, q0)");
		ASSERT(CMPULT128(x1, q1), "twopmodq128_96_q8 : CMPULT128(x1, q1)");
		ASSERT(CMPULT128(x2, q2), "twopmodq128_96_q8 : CMPULT128(x2, q2)");
		ASSERT(CMPULT128(x3, q3), "twopmodq128_96_q8 : CMPULT128(x3, q3)");
		ASSERT(CMPULT128(x4, q4), "twopmodq128_96_q8 : CMPULT128(x4, q4)");
		ASSERT(CMPULT128(x5, q5), "twopmodq128_96_q8 : CMPULT128(x5, q5)");
		ASSERT(CMPULT128(x6, q6), "twopmodq128_96_q8 : CMPULT128(x6, q6)");
		ASSERT(CMPULT128(x7, q7), "twopmodq128_96_q8 : CMPULT128(x7, q7)");
	#endif
		ADD128(x0, x0, x0);
		ADD128(x1, x1, x1);
		ADD128(x2, x2, x2);
		ADD128(x3, x3, x3);
		ADD128(x4, x4, x4);
		ADD128(x5, x5, x5);
		ADD128(x6, x6, x6);
		ADD128(x7, x7, x7);

		if(CMPULE128(q0, x0)) SUB128(x0, q0, x0);
		if(CMPULE128(q1, x1)) SUB128(x1, q1, x1);
		if(CMPULE128(q2, x2)) SUB128(x2, q2, x2);
		if(CMPULE128(q3, x3)) SUB128(x3, q3, x3);
		if(CMPULE128(q4, x4)) SUB128(x4, q4, x4);
		if(CMPULE128(q5, x5)) SUB128(x5, q5, x5);
		if(CMPULE128(q6, x6)) SUB128(x6, q6, x6);
		if(CMPULE128(q7, x7)) SUB128(x7, q7, x7);
	}

#if FAC_DEBUG
	if(CMPULT128(q0, x0)){ sprintf(char_buf, "twopmodq128_96_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint128_base10_char(str0, x0)], &str1[convert_uint128_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q1, x1)){ sprintf(char_buf, "twopmodq128_96_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint128_base10_char(str0, x1)], &str1[convert_uint128_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q2, x2)){ sprintf(char_buf, "twopmodq128_96_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint128_base10_char(str0, x2)], &str1[convert_uint128_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q3, x3)){ sprintf(char_buf, "twopmodq128_96_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint128_base10_char(str0, x3)], &str1[convert_uint128_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q4, x4)){ sprintf(char_buf, "twopmodq128_96_q8 : (x4 = %s) >= (q4 = %s)", &str0[convert_uint128_base10_char(str0, x4)], &str1[convert_uint128_base10_char(str1, q4)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q5, x5)){ sprintf(char_buf, "twopmodq128_96_q8 : (x5 = %s) >= (q5 = %s)", &str0[convert_uint128_base10_char(str0, x5)], &str1[convert_uint128_base10_char(str1, q5)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q6, x6)){ sprintf(char_buf, "twopmodq128_96_q8 : (x6 = %s) >= (q6 = %s)", &str0[convert_uint128_base10_char(str0, x6)], &str1[convert_uint128_base10_char(str1, q6)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q7, x7)){ sprintf(char_buf, "twopmodq128_96_q8 : (x7 = %s) >= (q7 = %s)", &str0[convert_uint128_base10_char(str0, x7)], &str1[convert_uint128_base10_char(str1, q7)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
#if FAC_DEBUG
if(dbg)printf("A: l = %20llu + 2^64* %20llu\n",lo0.d0,lo0.d1);
#endif
		/* Haven't gotten IA64 version of this working properly yet:
		SQR_LOHI_INPLACE128_96_q8(
		  x0, hi0
		, x1, hi1
		, x2, hi2
		, x3, hi3
		, x4, hi4
		, x5, hi5
		, x6, hi6
		, x7, hi7);
		*/
		SQR_LOHI128_96_q8(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3
		, x4, lo4, hi4
		, x5, lo5, hi5
		, x6, lo6, hi6
		, x7, lo7, hi7);
#if FAC_DEBUG
if(dbg)printf("B: l = %20llu + 2^64* %20llu\n",lo0.d0,lo0.d1);
#endif

		/* For unknown reasons, the 8-operand version of MULL128 was slower than one-at-a-time. */
		MULL128_INPLACE_q8(
		  lo0, qinv0
		, lo1, qinv1
		, lo2, qinv2
		, lo3, qinv3
		, lo4, qinv4
		, lo5, qinv5
		, lo6, qinv6
		, lo7, qinv7);
#if FAC_DEBUG
if(dbg)printf("C: l = %20llu + 2^64* %20llu\n",lo0.d0,lo0.d1);
#endif
	#if(USE_128x96 > 0)
		/* Need to be careful about the order of the 2 inputs here,
		as MULH128x96 assumes the 2nd input is the one which is < 2^96: */
		MULH128x96_q8(
		  lo0, q0, lo0
		, lo1, q1, lo1
		, lo2, q2, lo2
		, lo3, q3, lo3
		, lo4, q4, lo4
		, lo5, q5, lo5
		, lo6, q6, lo6
		, lo7, q7, lo7);
	#else
		MULH128_q8(
		  q0, lo0, lo0
		, q1, lo1, lo1
		, q2, lo2, lo2
		, q3, lo3, lo3
		, q4, lo4, lo4
		, q5, lo5, lo5
		, q6, lo6, lo6
		, q7, lo7, lo7);
	#endif
#if FAC_DEBUG
if(dbg)printf("D: l = %20llu + 2^64* %20llu\n",lo0.d0,lo0.d1);
#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
#if FAC_DEBUG
if(dbg)printf("On entry to (h<l): hi = %20llu\n",hi0);
if(dbg)printf("l = %20llu + 2^64* %20llu\n",lo0.d0,lo0.d1);
if(dbg)printf("x = %20llu + 2^64* %20llu\n",x0.d0,x0.d1);
#endif
		if(lo0.d1 != 0 || hi0 < lo0.d0){ SUB128(q0, lo0, x0);	x0.d0 +=  hi0; x0.d1 += (x0.d0 < hi0); } else { x0.d0 =  hi0 - lo0.d0; x0.d1 = (uint64)0; }
		if(lo1.d1 != 0 || hi1 < lo1.d0){ SUB128(q1, lo1, x1);	x1.d0 +=  hi1; x1.d1 += (x1.d0 < hi1); } else { x1.d0 =  hi1 - lo1.d0; x1.d1 = (uint64)0; }
		if(lo2.d1 != 0 || hi2 < lo2.d0){ SUB128(q2, lo2, x2);	x2.d0 +=  hi2; x2.d1 += (x2.d0 < hi2); } else { x2.d0 =  hi2 - lo2.d0; x2.d1 = (uint64)0; }
		if(lo3.d1 != 0 || hi3 < lo3.d0){ SUB128(q3, lo3, x3);	x3.d0 +=  hi3; x3.d1 += (x3.d0 < hi3); } else { x3.d0 =  hi3 - lo3.d0; x3.d1 = (uint64)0; }
		if(lo4.d1 != 0 || hi4 < lo4.d0){ SUB128(q4, lo4, x4);	x4.d0 +=  hi4; x4.d1 += (x4.d0 < hi4); } else { x4.d0 =  hi4 - lo4.d0; x4.d1 = (uint64)0; }
		if(lo5.d1 != 0 || hi5 < lo5.d0){ SUB128(q5, lo5, x5);	x5.d0 +=  hi5; x5.d1 += (x5.d0 < hi5); } else { x5.d0 =  hi5 - lo5.d0; x5.d1 = (uint64)0; }
		if(lo6.d1 != 0 || hi6 < lo6.d0){ SUB128(q6, lo6, x6);	x6.d0 +=  hi6; x6.d1 += (x6.d0 < hi6); } else { x6.d0 =  hi6 - lo6.d0; x6.d1 = (uint64)0; }
		if(lo7.d1 != 0 || hi7 < lo7.d0){ SUB128(q7, lo7, x7);	x7.d0 +=  hi7; x7.d1 += (x7.d0 < hi7); } else { x7.d0 =  hi7 - lo7.d0; x7.d1 = (uint64)0; }
#if FAC_DEBUG
if(dbg)printf("j = %2d, Res = %20llu + 2^64* %20llu",j,x0.d0,x0.d1);
#endif

		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(CMPULT128(x0, q0), "twopmodq128_96_q8 : CMPULT128(x0, q0)");
			ASSERT(CMPULT128(x1, q1), "twopmodq128_96_q8 : CMPULT128(x1, q1)");
			ASSERT(CMPULT128(x2, q2), "twopmodq128_96_q8 : CMPULT128(x2, q2)");
			ASSERT(CMPULT128(x3, q3), "twopmodq128_96_q8 : CMPULT128(x3, q3)");
			ASSERT(CMPULT128(x4, q4), "twopmodq128_96_q8 : CMPULT128(x4, q4)");
			ASSERT(CMPULT128(x5, q5), "twopmodq128_96_q8 : CMPULT128(x5, q5)");
			ASSERT(CMPULT128(x6, q6), "twopmodq128_96_q8 : CMPULT128(x6, q6)");
			ASSERT(CMPULT128(x7, q7), "twopmodq128_96_q8 : CMPULT128(x7, q7)");
		#endif
			ADD128(x0, x0, x0);
			ADD128(x1, x1, x1);
			ADD128(x2, x2, x2);
			ADD128(x3, x3, x3);
			ADD128(x4, x4, x4);
			ADD128(x5, x5, x5);
			ADD128(x6, x6, x6);
			ADD128(x7, x7, x7);

			if(CMPULE128(q0, x0)) SUB128(x0, q0, x0);
			if(CMPULE128(q1, x1)) SUB128(x1, q1, x1);
			if(CMPULE128(q2, x2)) SUB128(x2, q2, x2);
			if(CMPULE128(q3, x3)) SUB128(x3, q3, x3);
			if(CMPULE128(q4, x4)) SUB128(x4, q4, x4);
			if(CMPULE128(q5, x5)) SUB128(x5, q5, x5);
			if(CMPULE128(q6, x6)) SUB128(x6, q6, x6);
			if(CMPULE128(q7, x7)) SUB128(x7, q7, x7);
#if FAC_DEBUG
if(dbg)printf(" *2 = %20llu + 2^64* %20llu",x0.d0,x0.d1);
#endif
		}
#if FAC_DEBUG
if(dbg)printf("\n");
#endif
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	ADD128(x0 ,x0, x0);
	ADD128(x1 ,x1, x1);
	ADD128(x2 ,x2, x2);
	ADD128(x3 ,x3, x3);
	ADD128(x4 ,x4, x4);
	ADD128(x5 ,x5, x5);
	ADD128(x6 ,x6, x6);
	ADD128(x7 ,x7, x7);
	//*** Aug 2022: see note at end of twopmodq96() for why we can add-sans-carry the constant (0 or 2) here:
	x0.d0 += FERMAT;
	x1.d0 += FERMAT;
	x2.d0 += FERMAT;
	x3.d0 += FERMAT;
	x4.d0 += FERMAT;
	x5.d0 += FERMAT;
	x6.d0 += FERMAT;
	x7.d0 += FERMAT;
	SUB128(x0, q0, x0);
	SUB128(x1, q1, x1);
	SUB128(x2, q2, x2);
	SUB128(x3, q3, x3);
	SUB128(x4, q4, x4);
	SUB128(x5, q5, x5);
	SUB128(x6, q6, x6);
	SUB128(x7, q7, x7);

#if FAC_DEBUG
if(dbg)printf("x0 = %20llu + 2^64* %20llu\n",x0.d0, x0.d1);
#endif

	/* Only do the full 128-bit (Xj== 1) check if the bottom 64 bits of Xj == 1: */
	r = 0;
	if(x0.d0 == 1) r += ((uint64)CMPEQ128(x0, ONE128) << 0);
	if(x1.d0 == 1) r += ((uint64)CMPEQ128(x1, ONE128) << 1);
	if(x2.d0 == 1) r += ((uint64)CMPEQ128(x2, ONE128) << 2);
	if(x3.d0 == 1) r += ((uint64)CMPEQ128(x3, ONE128) << 3);
	if(x4.d0 == 1) r += ((uint64)CMPEQ128(x4, ONE128) << 4);
	if(x5.d0 == 1) r += ((uint64)CMPEQ128(x5, ONE128) << 5);
	if(x6.d0 == 1) r += ((uint64)CMPEQ128(x6, ONE128) << 6);
	if(x7.d0 == 1) r += ((uint64)CMPEQ128(x7, ONE128) << 7);
	return(r);
}

