#include "factor.h"
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
uint96 twopmodq96(uint64 p, uint96 q)
{
#if FAC_DEBUG
	int dbg = (p == 16219289);
	uint128 y;
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead7,  hi64;
	uint96 qhalf, qinv, x, lo, hi;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq96:\n");
#endif

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
	/* qinv has 128 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + __MULH64(q.d0, qinv.d0));
#else
	MULH64(q.d0, qinv.d0, hi64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + hi64);
#endif
	qinv.d1 &= 0x00000000ffffffff;	/* Only want the lower 32 bits here */

#if FAC_DEBUG
	ASSERT(HERE, qinv.d1 == x.d1 && qinv.d0 == x.d0, "twopmodq96 : qinv.d1 == x.d1 && qinv.d0 == x.d0");
	if(dbg) printf("q    = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q   )]);
	if(dbg) printf("qinv = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv)]);
#endif

	/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("zshift  = %u\n", zshift);
#endif
	LSHIFT96(qinv, zshift, lo);

#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	MULH96(q,lo,lo);

#if FAC_DEBUG
	if(dbg) printf("q*lo/2^96 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB96(q, lo, x);	/* Put the result in lo (rather than x), to ease overflow check below */

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
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
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
#if FAC_DEBUG
	if(CMPULT96(q, x)){ sprintf(char_buf, "twopmodq96 : (x0 = %s) >= (q = %s)", &str0[convert_uint128_base10_char(str0, x)], &str1[convert_uint128_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI96(x,lo,hi);
		MULL96(lo,qinv,lo);
		MULH96(q,lo,lo);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT96(hi, lo))
		{
			SUB96(q, lo, lo);
			ADD96(lo, hi, x);
		}
		else
		{
			SUB96(hi, lo, x);
		}

#if FAC_DEBUG
if(dbg)printf("j = %2d, x = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT96(x, q), "twopmodq96 : CMPULT96(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x, qhalf)){ ADD96(x, x, x); SUB96(x, q, x); }else{ ADD96(x, x, x); }

#if FAC_DEBUG
	if(dbg) printf("*2= %s", &char_buf[convert_uint128_base10_char(char_buf, x)]);
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
if(dbg)printf("x0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

	return x;
}

/*** 4-trial-factor version ***/
uint64 twopmodq96_q4(uint64 p, uint96 q0, uint96 q1, uint96 q2, uint96 q3)
{
#if FAC_DEBUG
	int dbg = (p == 16219289);
/*	int dbg = STREQ(&char_buf[convert_uint64_base10_char(char_buf, p)], "0");*/
#endif
	 int32 j;
	uint64 lead7, tmp0, tmp1, tmp2, tmp3, r;
	uint96 qinv0, qinv1, qinv2, qinv3
		, qhalf0, qhalf1, qhalf2, qhalf3
		, x0, x1, x2, x3
		, lo0, lo1, lo2, lo3
		, hi0, hi1, hi2, hi3;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq96_q4:\n");
#endif
	RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST96(q1, 1, qhalf1);
	RSHIFT_FAST96(q2, 1, qhalf2);
	RSHIFT_FAST96(q3, 1, qhalf3);

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

	/*
	!    Find modular inverse (mod 2^96) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q0.d0 & (uint64)1) == 1, "twopmodq96_q4 : (q0.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q1.d0 & (uint64)1) == 1, "twopmodq96_q4 : (q1.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q2.d0 & (uint64)1) == 1, "twopmodq96_q4 : (q2.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q3.d0 & (uint64)1) == 1, "twopmodq96_q4 : (q3.d0 & (uint64)1) == 1");
#endif
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d1 = (uint64)0;

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

		qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
		qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
		qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);
		qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);
	}

	/* Now that have bottom 64 bits of qinv, do one more Newton iteration
	using full 96-bit operands. See twopmodq96 for details on streamlining here.
	*/
	/* qinv has 128 bits, but only the upper 64 get modified here. */
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
	qinv0.d1 &= 0x00000000ffffffff;	/* Only want the lower 32 bits here */
	qinv1.d1 &= 0x00000000ffffffff;
	qinv2.d1 &= 0x00000000ffffffff;
	qinv3.d1 &= 0x00000000ffffffff;

#if FAC_DEBUG
	if(dbg)
	{
		printf("q   0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q0   )]);
		printf("qinv0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv0)]);
		printf("q   1 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q1   )]);
		printf("qinv1 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv1)]);
		printf("q   2 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q2   )]);
		printf("qinv2 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv2)]);
		printf("q   3 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q3   )]);
		printf("qinv3 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv3)]);
	}
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
	if(dbg) printf("lo0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo0)]);
	if(dbg) printf("lo1 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo1)]);
	if(dbg) printf("lo2 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo2)]);
	if(dbg) printf("lo3 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo3)]);
#endif

	MULH96_q4(
	  q0, lo0, lo0
	, q1, lo1, lo1
	, q2, lo2, lo2
	, q3, lo3, lo3);
#if FAC_DEBUG
	if(dbg) printf("q0*lo0/2^96 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo0)]);
	if(dbg) printf("q1*lo1/2^96 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo1)]);
	if(dbg) printf("q2*lo2/2^96 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo2)]);
	if(dbg) printf("q3*lo3/2^96 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo3)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB96(q0, lo0, x0);
	SUB96(q1, lo1, x1);
	SUB96(q2, lo2, x2);
	SUB96(q3, lo3, x3);
#if FAC_DEBUG
	if(dbg) printf("x0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x0)]);
	if(dbg) printf("x1 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x1)]);
	if(dbg) printf("x2 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x2)]);
	if(dbg) printf("x3 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x3)]);
#endif

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
	}

#if FAC_DEBUG
	if(CMPULT96(q0, x0)){ sprintf(char_buf, "twopmodq96_q4 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint128_base10_char(str0, x0)], &str1[convert_uint128_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q1, x1)){ sprintf(char_buf, "twopmodq96_q4 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint128_base10_char(str0, x1)], &str1[convert_uint128_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q2, x2)){ sprintf(char_buf, "twopmodq96_q4 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint128_base10_char(str0, x2)], &str1[convert_uint128_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q3, x3)){ sprintf(char_buf, "twopmodq96_q4 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint128_base10_char(str0, x3)], &str1[convert_uint128_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI96_q4(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3);

		MULL96_q4(
		  lo0, qinv0, lo0
		, lo1, qinv1, lo1
		, lo2, qinv2, lo2
		, lo3, qinv3, lo3);

		MULH96_q4(
		  q0, lo0, lo0
		, q1, lo1, lo1
		, q2, lo2, lo2
		, q3, lo3, lo3);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT96(hi0, lo0)) { SUB96(q0, lo0, lo0);	ADD96(lo0, hi0, x0); } else { SUB96(hi0, lo0, x0); }
		if(CMPULT96(hi1, lo1)) { SUB96(q1, lo1, lo1);	ADD96(lo1, hi1, x1); } else { SUB96(hi1, lo1, x1); }
		if(CMPULT96(hi2, lo2)) { SUB96(q2, lo2, lo2);	ADD96(lo2, hi2, x2); } else { SUB96(hi2, lo2, x2); }
		if(CMPULT96(hi3, lo3)) { SUB96(q3, lo3, lo3);	ADD96(lo3, hi3, x3); } else { SUB96(hi3, lo3, x3); }

#if FAC_DEBUG
if(dbg)
{
	printf("j = %2d, x0 = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x0)]);
	printf("j = %2d, x1 = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x1)]);
	printf("j = %2d, x2 = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x2)]);
	printf("j = %2d, x3 = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x3)]);
}
#endif

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
		}
	}

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

#if FAC_DEBUG
if(dbg)
{
	printf("x0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x0)]);
	printf("x1 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x1)]);
	printf("x2 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x2)]);
	printf("x3 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x3)]);
}
#endif

	r = 0;
	if(CMPEQ96(x0, ONE96)) r +=  1;
	if(CMPEQ96(x1, ONE96)) r +=  2;
	if(CMPEQ96(x2, ONE96)) r +=  4;
	if(CMPEQ96(x3, ONE96)) r +=  8;
	return(r);
}

/*** 8-trial-factor version ***/
uint64 twopmodq96_q8(uint64 p, uint96 q0, uint96 q1, uint96 q2, uint96 q3, uint96 q4, uint96 q5, uint96 q6, uint96 q7)
{
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint64_base10_char(char_buf, p)], "0");
#endif
	 int32 j;
	uint64 lead7, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r;
	uint96 qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, qhalf0, qhalf1, qhalf2, qhalf3, qhalf4, qhalf5, qhalf6, qhalf7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
		, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq96_q8:\n");
#endif
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
	if(CMPULT96(q0, x0)){ sprintf(char_buf, "twopmodq96_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint128_base10_char(str0, x0)], &str1[convert_uint128_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q1, x1)){ sprintf(char_buf, "twopmodq96_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint128_base10_char(str0, x1)], &str1[convert_uint128_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q2, x2)){ sprintf(char_buf, "twopmodq96_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint128_base10_char(str0, x2)], &str1[convert_uint128_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q3, x3)){ sprintf(char_buf, "twopmodq96_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint128_base10_char(str0, x3)], &str1[convert_uint128_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q4, x4)){ sprintf(char_buf, "twopmodq96_q8 : (x4 = %s) >= (q4 = %s)", &str0[convert_uint128_base10_char(str0, x4)], &str1[convert_uint128_base10_char(str1, q4)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q5, x5)){ sprintf(char_buf, "twopmodq96_q8 : (x5 = %s) >= (q5 = %s)", &str0[convert_uint128_base10_char(str0, x5)], &str1[convert_uint128_base10_char(str1, q5)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q6, x6)){ sprintf(char_buf, "twopmodq96_q8 : (x6 = %s) >= (q6 = %s)", &str0[convert_uint128_base10_char(str0, x6)], &str1[convert_uint128_base10_char(str1, q6)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT96(q7, x7)){ sprintf(char_buf, "twopmodq96_q8 : (x7 = %s) >= (q7 = %s)", &str0[convert_uint128_base10_char(str0, x7)], &str1[convert_uint128_base10_char(str1, q7)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
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

