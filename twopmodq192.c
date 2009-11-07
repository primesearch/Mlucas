#include "factor.h"

#if FAC_DEBUG
	char char_buf[1024], str0[64], str1[64];
#endif

/***********************************************************************************/
/***192-BIT INPUTS *****************************************************************/
/***********************************************************************************/
/*
Function to find 2^(-p) mod q, where p and q are both 192-bit unsigned integers.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^192 (i.e. our MODQ
operation effects multiply modulo 2^192).

The key 3-operation inner-loop sequence here is as follows:

loop{
	SQR_LOHI192(x,lo,hi);	// lo + hi*2^192 = x^2. Input x has 192 bits; lo/hi have 192 bits
	MULL192(lo,qinv,lo);	// lo = lo*qinv, where qinv is the mod-2^192 inverse of q. High half of result discarded.
	MULH192(q,lo,lo);		// lo = UMULH(lo, q). Low half of result discarded.
	// If hi < lo, then calculate q - lo + hi < q; otherwise calculate h-l.
}
*/
uint192 twopmodq192(uint192 p, uint192 q)
{
#if FAC_DEBUG
	int dbg = TRUE;//STREQ(&char_buf[convert_uint192_base10_char(char_buf, p)], "0");
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead8, lo64;
	uint192 qhalf, qinv, x, lo, hi;
	static uint192 psave = {0ull,0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)
	printf("twopmodq192:\n");
#endif

	RSHIFT_FAST192(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || !CMPEQ192(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		x.d0 = 192; x.d1 = x.d2 = 0;
		ADD192(p, x, pshift);
#if FAC_DEBUG
	if(dbg)
	{
		printf("p = %s\n", &char_buf[convert_uint192_base10_char(char_buf, p     )]);
		printf("p+= %s\n", &char_buf[convert_uint192_base10_char(char_buf, pshift)]);
	}
#endif
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 7 or 8 to account for the fact that we can do the powering for the leftmost
	!    7 or 8 bits (depending on whether the leftmost 8 > 191 or not) via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming
		the leftmost 7/8 bits have already been processed via a shift.

		Since 7 bits with leftmost bit = 1 is guaranteed
		to be in [64,127], the shift count here is in [0, 63].
		That means that zstart < 2^64. Together with the fact that
		squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration.
		Every little bit counts (literally in this case :), right?
		*/
		if(pshift.d2)
		{
			j = leadz64(pshift.d2);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = (((pshift.d2<<j) + (pshift.d1>>(64-j))) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index = 192-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index = 192-j-8;
		}
		else if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index = 128-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index = 128-j-8;
		}
		else
		{
			j = leadz64(pshift.d0);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = ((pshift.d0<<j) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index =  64-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index =  64-j-8;
		}

#if FAC_DEBUG
	if(dbg)	printf("lead8 = %u\n", (uint32)lead8);
#endif
		zshift = 191 - lead8;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift.d2 = ~pshift.d2;	pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/*
	!    Find modular inverse (mod 2^192) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q.d0 & (uint64)1) == 1, "twopmodq192 : (q.d0 & (uint64)1) == 1");
#endif
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv.d0 = (q.d0 + q.d0 + q.d0) ^ (uint64)2;	qinv.d2 = qinv.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*(2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
	defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
	*/
	for(j = 0; j < 4; j++)
	{
		lo64 = q.d0*qinv.d0;
		qinv.d0 = qinv.d0*((uint64)2 - lo64);
	}

	/* Now that have bottom 64 bits of qinv, do two more Newton iterations using full 192-bit operands: */
	MULL192(q, qinv, x);
	SUB192 (TWO192, x, x);
	MULL192(qinv, x, qinv);

	MULL192(q, qinv, x);
	SUB192 (TWO192, x, x);
	MULL192(qinv, x, qinv);

#if FAC_DEBUG
	if(dbg)
	{
		printf("q    = %s\n", &char_buf[convert_uint192_base10_char(char_buf, q   )]);
		printf("qinv = %s\n", &char_buf[convert_uint192_base10_char(char_buf, qinv)]);
	}
#endif
	/* Since zstart is a power of two < 2^192, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL192(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("j = start_index - 1 = %u\n", j);
	if(dbg) printf("zshift  = %u\n", zshift);
#endif
	LSHIFT192(qinv, zshift, lo);
#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint192_base10_char(char_buf, lo)]);
#endif
	MULH192(q,lo,lo);
#if FAC_DEBUG
	if(dbg) printf("q*lo/2^192 = %s\n", &char_buf[convert_uint192_base10_char(char_buf, lo)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB192(q, lo, x);

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
	if(TEST_BIT192(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT192(x,q), "twopmodq192 : CMPULT192(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT192(x, qhalf)){ ADD192(x, x, x); SUB192(x, q, x); }else{ ADD192(x, x, x); }
#if FAC_DEBUG
if(dbg) printf("2x= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
	}

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
#if FAC_DEBUG
	if(CMPULT192(q, x)){ sprintf(char_buf, "twopmodq192 : (x0 = %s) >= (q = %s)", &str0[convert_uint192_base10_char(str0, x)], &str1[convert_uint192_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI192(x,lo,hi);
		MULL192(lo,qinv,lo);
		MULH192(q,lo,lo);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT192(hi, lo))
		{
			SUB192(q, lo, lo);
			ADD192(lo, hi, x);
		}
		else
		{
			SUB192(hi, lo, x);
		}
#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif

		if(TEST_BIT192(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT192(x,q), "twopmodq192 : CMPULT192(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT192(x, qhalf)){ ADD192(x, x, x); SUB192(x, q, x); }else{ ADD192(x, x, x); }
#if FAC_DEBUG
	if(dbg) printf("2x= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD192(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^191, so x + x cannot overflow. */
#if FAC_DEBUG
	if(dbg) printf("Final x = %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
	SUB192(x,q,x);
#if FAC_DEBUG
	if(dbg) printf("Final x-q=%s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif

	return x;
}

/*** 4-trial-factor version ***/
uint64 twopmodq192_q4(uint192 p, uint192 q0, uint192 q1, uint192 q2, uint192 q3)
{
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, lead8, r;
	uint192 qinv0, qinv1, qinv2, qinv3
			, qhalf0, qhalf1, qhalf2, qhalf3
			, x0, x1, x2, x3
			, lo0, lo1, lo2, lo3
			, hi0, hi1, hi2, hi3;
	uint192 x;
	static uint192 psave = {0ull,0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	RSHIFT_FAST192(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST192(q1, 1, qhalf1);
	RSHIFT_FAST192(q2, 1, qhalf2);
	RSHIFT_FAST192(q3, 1, qhalf3);

	if(first_entry || !CMPEQ192(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		x.d0 = 192; x.d1 = x.d2 = 0;
		ADD192(p, x, pshift);
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 7 or 8 to account for the fact that we can do the powering for the leftmost
	!    7 or 8 bits (depending on whether the leftmost 8 > 191 or not) via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming
		the leftmost 7/8 bits have already been processed via a shift.

		Since 7 bits with leftmost bit = 1 is guaranteed
		to be in [64,127], the shift count here is in [0, 63].
		That means that zstart < 2^64. Together with the fact that
		squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration.
		Every little bit counts (literally in this case :), right?
		*/
		if(pshift.d2)
		{
			j = leadz64(pshift.d2);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = (((pshift.d2<<j) + (pshift.d1>>(64-j))) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index = 192-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index = 192-j-8;
		}
		else if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index = 128-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index = 128-j-8;
		}
		else
		{
			j = leadz64(pshift.d0);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = ((pshift.d0<<j) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index =  64-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index =  64-j-8;
		}

  		zshift = 191 - lead8;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift.d2 = ~pshift.d2;	pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/*
	!    Find modular inverse (mod 2^192) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q0.d0 & (uint64)1) == 1, "twopmodq192_q4 : (q0.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q1.d0 & (uint64)1) == 1, "twopmodq192_q4 : (q1.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q2.d0 & (uint64)1) == 1, "twopmodq192_q4 : (q2.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q3.d0 & (uint64)1) == 1, "twopmodq192_q4 : (q3.d0 & (uint64)1) == 1");
#endif
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d2 = qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d2 = qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d2 = qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d2 = qinv3.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*(2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
	defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
	*/
	for(j = 0; j < 4; j++)
	{
		lo64_0 = q0.d0*qinv0.d0;
		lo64_1 = q1.d0*qinv1.d0;
		lo64_2 = q2.d0*qinv2.d0;
		lo64_3 = q3.d0*qinv3.d0;

		qinv0.d0 = qinv0.d0*((uint64)2 - lo64_0);
		qinv1.d0 = qinv1.d0*((uint64)2 - lo64_1);
		qinv2.d0 = qinv2.d0*((uint64)2 - lo64_2);
		qinv3.d0 = qinv3.d0*((uint64)2 - lo64_3);
	}

	/* Now that have bottom 64 bits of qinv, do two more Newton iterations using full 192-bit operands: */
	MULL192(q0, qinv0, x0);	SUB192 (TWO192, x0, x0);	MULL192(qinv0, x0, qinv0);
	MULL192(q1, qinv1, x1);	SUB192 (TWO192, x1, x1);	MULL192(qinv1, x1, qinv1);
	MULL192(q2, qinv2, x2);	SUB192 (TWO192, x2, x2);	MULL192(qinv2, x2, qinv2);
	MULL192(q3, qinv3, x3);	SUB192 (TWO192, x3, x3);	MULL192(qinv3, x3, qinv3);

	MULL192(q0, qinv0, x0);	SUB192 (TWO192, x0, x0);	MULL192(qinv0, x0, qinv0);
	MULL192(q1, qinv1, x1);	SUB192 (TWO192, x1, x1);	MULL192(qinv1, x1, qinv1);
	MULL192(q2, qinv2, x2);	SUB192 (TWO192, x2, x2);	MULL192(qinv2, x2, qinv2);
	MULL192(q3, qinv3, x3);	SUB192 (TWO192, x3, x3);	MULL192(qinv3, x3, qinv3);

	/* Since zstart is a power of two < 2^192, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL192(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv.
	hi = 0 in this instance, which simplifies things in the final subtract step.
	*/
	LSHIFT192(qinv0, zshift, lo0);	MULH192(q0,lo0,lo0);	SUB192(q0, lo0, x0);
	LSHIFT192(qinv1, zshift, lo1);	MULH192(q1,lo1,lo1);	SUB192(q1, lo1, x1);
	LSHIFT192(qinv2, zshift, lo2);	MULH192(q2,lo2,lo2);	SUB192(q2, lo2, x2);
	LSHIFT192(qinv3, zshift, lo3);	MULH192(q3,lo3,lo3);	SUB192(q3, lo3, x3);

	if(TEST_BIT192(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT192(x0, q0), "twopmodq192_q4 : CMPULT192(x0, q0)");
		ASSERT(HERE, CMPULT192(x1, q1), "twopmodq192_q4 : CMPULT192(x1, q1)");
		ASSERT(HERE, CMPULT192(x2, q2), "twopmodq192_q4 : CMPULT192(x2, q2)");
		ASSERT(HERE, CMPULT192(x3, q3), "twopmodq192_q4 : CMPULT192(x3, q3)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
		if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
		if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
		if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
	}

#if FAC_DEBUG
	if(CMPULT192(q0, x0)){ sprintf(char_buf, "twopmodq192_q4 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint192_base10_char(str0, x0)], &str1[convert_uint192_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q1, x1)){ sprintf(char_buf, "twopmodq192_q4 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint192_base10_char(str0, x1)], &str1[convert_uint192_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q2, x2)){ sprintf(char_buf, "twopmodq192_q4 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint192_base10_char(str0, x2)], &str1[convert_uint192_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q3, x3)){ sprintf(char_buf, "twopmodq192_q4 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint192_base10_char(str0, x3)], &str1[convert_uint192_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
	/*...x^2 mod q is returned in x. */
	#ifdef PIPELINE_MUL192
		SQR_LOHI192_q4(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3);

		MULL192_q4(
		  lo0, qinv0, lo0
		, lo1, qinv1, lo1
		, lo2, qinv2, lo2
		, lo3, qinv3, lo3);

		MULH192_q4(
		  q0, lo0, lo0
		, q1, lo1, lo1
		, q2, lo2, lo2
		, q3, lo3, lo3);
	#else
	/* Strangely, on alpha/ia64 nopipe is better than the pipelined ..._q4 version: */
		SQR_LOHI192(x0,lo0,hi0);	MULL192(lo0,qinv0,lo0);		MULH192(q0,lo0,lo0);
		SQR_LOHI192(x1,lo1,hi1);	MULL192(lo1,qinv1,lo1);		MULH192(q1,lo1,lo1);
		SQR_LOHI192(x2,lo2,hi2);	MULL192(lo2,qinv2,lo2);		MULH192(q2,lo2,lo2);
		SQR_LOHI192(x3,lo3,hi3);	MULL192(lo3,qinv3,lo3);		MULH192(q3,lo3,lo3);
	#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT192(hi0, lo0)) { SUB192(q0, lo0, lo0);	ADD192(lo0, hi0, x0); } else { SUB192(hi0, lo0, x0); }
		if(CMPULT192(hi1, lo1)) { SUB192(q1, lo1, lo1);	ADD192(lo1, hi1, x1); } else { SUB192(hi1, lo1, x1); }
		if(CMPULT192(hi2, lo2)) { SUB192(q2, lo2, lo2);	ADD192(lo2, hi2, x2); } else { SUB192(hi2, lo2, x2); }
		if(CMPULT192(hi3, lo3)) { SUB192(q3, lo3, lo3);	ADD192(lo3, hi3, x3); } else { SUB192(hi3, lo3, x3); }

		if(TEST_BIT192(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT192(x0,q0), "twopmodq192_q4 : CMPULT192(x0,q0)");
			ASSERT(HERE, CMPULT192(x1,q1), "twopmodq192_q4 : CMPULT192(x1,q1)");
			ASSERT(HERE, CMPULT192(x2,q2), "twopmodq192_q4 : CMPULT192(x2,q2)");
			ASSERT(HERE, CMPULT192(x3,q3), "twopmodq192_q4 : CMPULT192(x3,q3)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
			if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
			if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
			if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD192(x0,x0,x0);	SUB192(x0,q0,x0);
	ADD192(x1,x1,x1);	SUB192(x1,q1,x1);
	ADD192(x2,x2,x2);	SUB192(x2,q2,x2);
	ADD192(x3,x3,x3);	SUB192(x3,q3,x3);

	/* Only do the full 192-bit (Xj== 1) check if the bottom 64 bits of Xj == 1: */
	r = 0;
	if(x0.d0 == 1) r += ((uint64)CMPEQ192(x0, ONE192) << 0);
	if(x1.d0 == 1) r += ((uint64)CMPEQ192(x1, ONE192) << 1);
	if(x2.d0 == 1) r += ((uint64)CMPEQ192(x2, ONE192) << 2);
	if(x3.d0 == 1) r += ((uint64)CMPEQ192(x3, ONE192) << 3);
	return(r);
}

/*** 8-trial-factor version ***/
uint64 twopmodq192_q8(uint192 p, uint192 q0, uint192 q1, uint192 q2, uint192 q3, uint192 q4, uint192 q5, uint192 q6, uint192 q7)
{
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, lo64_4, lo64_5, lo64_6, lo64_7, lead8, r;
	uint192 qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
			, qhalf0, qhalf1, qhalf2, qhalf3, qhalf4, qhalf5, qhalf6, qhalf7
			, x0, x1, x2, x3, x4, x5, x6, x7
			, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
			, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7;
	uint192 x;
	static uint192 psave = {0ull,0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	RSHIFT_FAST192(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST192(q1, 1, qhalf1);
	RSHIFT_FAST192(q2, 1, qhalf2);
	RSHIFT_FAST192(q3, 1, qhalf3);
	RSHIFT_FAST192(q4, 1, qhalf4);
	RSHIFT_FAST192(q5, 1, qhalf5);
	RSHIFT_FAST192(q6, 1, qhalf6);
	RSHIFT_FAST192(q7, 1, qhalf7);

	if(first_entry || !CMPEQ192(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		x.d0 = 192; x.d1 = x.d2 = 0;
		ADD192(p, x, pshift);
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 7 or 8 to account for the fact that we can do the powering for the leftmost
	!    7 or 8 bits (depending on whether the leftmost 8 > 191 or not) via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming
		the leftmost 7/8 bits have already been processed via a shift.

		Since 7 bits with leftmost bit = 1 is guaranteed
		to be in [64,127], the shift count here is in [0, 63].
		That means that zstart < 2^64. Together with the fact that
		squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration.
		Every little bit counts (literally in this case :), right?
		*/
		if(pshift.d2)
		{
			j = leadz64(pshift.d2);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = (((pshift.d2<<j) + (pshift.d1>>(64-j))) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index = 192-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index = 192-j-8;
		}
		else if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index = 128-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index = 128-j-8;
		}
		else
		{
			j = leadz64(pshift.d0);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = ((pshift.d0<<j) >> 56);
			if(lead8 > 191)
			{
				lead8 >>= 1;
				start_index =  64-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index =  64-j-8;
		}

  		zshift = 191 - lead8;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift.d2 = ~pshift.d2;	pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/*
	!    Find modular inverse (mod 2^192) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q0.d0 & (uint64)1) == 1, "twopmodq192_q8 : (q0.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q1.d0 & (uint64)1) == 1, "twopmodq192_q8 : (q1.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q2.d0 & (uint64)1) == 1, "twopmodq192_q8 : (q2.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q3.d0 & (uint64)1) == 1, "twopmodq192_q8 : (q3.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q4.d0 & (uint64)1) == 1, "twopmodq192_q8 : (q4.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q5.d0 & (uint64)1) == 1, "twopmodq192_q8 : (q5.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q6.d0 & (uint64)1) == 1, "twopmodq192_q8 : (q6.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q7.d0 & (uint64)1) == 1, "twopmodq192_q8 : (q7.d0 & (uint64)1) == 1");
#endif
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d2 = qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d2 = qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d2 = qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d2 = qinv3.d1 = (uint64)0;
	qinv4.d0 = (q4.d0 + q4.d0 + q4.d0) ^ (uint64)2;	qinv4.d2 = qinv4.d1 = (uint64)0;
	qinv5.d0 = (q5.d0 + q5.d0 + q5.d0) ^ (uint64)2;	qinv5.d2 = qinv5.d1 = (uint64)0;
	qinv6.d0 = (q6.d0 + q6.d0 + q6.d0) ^ (uint64)2;	qinv6.d2 = qinv6.d1 = (uint64)0;
	qinv7.d0 = (q7.d0 + q7.d0 + q7.d0) ^ (uint64)2;	qinv7.d2 = qinv7.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*(2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
	defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
	*/
	for(j = 0; j < 4; j++)
	{
		lo64_0 = q0.d0*qinv0.d0;
		lo64_1 = q1.d0*qinv1.d0;
		lo64_2 = q2.d0*qinv2.d0;
		lo64_3 = q3.d0*qinv3.d0;
		lo64_4 = q4.d0*qinv4.d0;
		lo64_5 = q5.d0*qinv5.d0;
		lo64_6 = q6.d0*qinv6.d0;
		lo64_7 = q7.d0*qinv7.d0;

		qinv0.d0 = qinv0.d0*((uint64)2 - lo64_0);
		qinv1.d0 = qinv1.d0*((uint64)2 - lo64_1);
		qinv2.d0 = qinv2.d0*((uint64)2 - lo64_2);
		qinv3.d0 = qinv3.d0*((uint64)2 - lo64_3);
		qinv4.d0 = qinv4.d0*((uint64)2 - lo64_4);
		qinv5.d0 = qinv5.d0*((uint64)2 - lo64_5);
		qinv6.d0 = qinv6.d0*((uint64)2 - lo64_6);
		qinv7.d0 = qinv7.d0*((uint64)2 - lo64_7);
	}

	/* Now that have bottom 64 bits of qinv, do two more Newton iterations using full 192-bit operands: */
	MULL192(q0, qinv0, x0);	SUB192 (TWO192, x0, x0);	MULL192(qinv0, x0, qinv0);
	MULL192(q1, qinv1, x1);	SUB192 (TWO192, x1, x1);	MULL192(qinv1, x1, qinv1);
	MULL192(q2, qinv2, x2);	SUB192 (TWO192, x2, x2);	MULL192(qinv2, x2, qinv2);
	MULL192(q3, qinv3, x3);	SUB192 (TWO192, x3, x3);	MULL192(qinv3, x3, qinv3);
	MULL192(q4, qinv4, x4);	SUB192 (TWO192, x4, x4);	MULL192(qinv4, x4, qinv4);
	MULL192(q5, qinv5, x5);	SUB192 (TWO192, x5, x5);	MULL192(qinv5, x5, qinv5);
	MULL192(q6, qinv6, x6);	SUB192 (TWO192, x6, x6);	MULL192(qinv6, x6, qinv6);
	MULL192(q7, qinv7, x7);	SUB192 (TWO192, x7, x7);	MULL192(qinv7, x7, qinv7);

	MULL192(q0, qinv0, x0);	SUB192 (TWO192, x0, x0);	MULL192(qinv0, x0, qinv0);
	MULL192(q1, qinv1, x1);	SUB192 (TWO192, x1, x1);	MULL192(qinv1, x1, qinv1);
	MULL192(q2, qinv2, x2);	SUB192 (TWO192, x2, x2);	MULL192(qinv2, x2, qinv2);
	MULL192(q3, qinv3, x3);	SUB192 (TWO192, x3, x3);	MULL192(qinv3, x3, qinv3);
	MULL192(q4, qinv4, x4);	SUB192 (TWO192, x4, x4);	MULL192(qinv4, x4, qinv4);
	MULL192(q5, qinv5, x5);	SUB192 (TWO192, x5, x5);	MULL192(qinv5, x5, qinv5);
	MULL192(q6, qinv6, x6);	SUB192 (TWO192, x6, x6);	MULL192(qinv6, x6, qinv6);
	MULL192(q7, qinv7, x7);	SUB192 (TWO192, x7, x7);	MULL192(qinv7, x7, qinv7);

	/* Since zstart is a power of two < 2^192, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL192(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv.
	hi = 0 in this instance, which simplifies things in the final subtract step.
	*/
	LSHIFT192(qinv0, zshift, lo0);	MULH192(q0,lo0,lo0);	SUB192(q0, lo0, x0);
	LSHIFT192(qinv1, zshift, lo1);	MULH192(q1,lo1,lo1);	SUB192(q1, lo1, x1);
	LSHIFT192(qinv2, zshift, lo2);	MULH192(q2,lo2,lo2);	SUB192(q2, lo2, x2);
	LSHIFT192(qinv3, zshift, lo3);	MULH192(q3,lo3,lo3);	SUB192(q3, lo3, x3);
	LSHIFT192(qinv4, zshift, lo4);	MULH192(q4,lo4,lo4);	SUB192(q4, lo4, x4);
	LSHIFT192(qinv5, zshift, lo5);	MULH192(q5,lo5,lo5);	SUB192(q5, lo5, x5);
	LSHIFT192(qinv6, zshift, lo6);	MULH192(q6,lo6,lo6);	SUB192(q6, lo6, x6);
	LSHIFT192(qinv7, zshift, lo7);	MULH192(q7,lo7,lo7);	SUB192(q7, lo7, x7);

	if(TEST_BIT192(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT192(x0, q0), "twopmodq192_q8 : CMPULT192(x0, q0)");
		ASSERT(HERE, CMPULT192(x1, q1), "twopmodq192_q8 : CMPULT192(x1, q1)");
		ASSERT(HERE, CMPULT192(x2, q2), "twopmodq192_q8 : CMPULT192(x2, q2)");
		ASSERT(HERE, CMPULT192(x3, q3), "twopmodq192_q8 : CMPULT192(x3, q3)");
		ASSERT(HERE, CMPULT192(x4, q4), "twopmodq192_q8 : CMPULT192(x4, q4)");
		ASSERT(HERE, CMPULT192(x5, q5), "twopmodq192_q8 : CMPULT192(x5, q5)");
		ASSERT(HERE, CMPULT192(x6, q6), "twopmodq192_q8 : CMPULT192(x6, q6)");
		ASSERT(HERE, CMPULT192(x7, q7), "twopmodq192_q8 : CMPULT192(x7, q7)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
		if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
		if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
		if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
		if(CMPUGT192(x4, qhalf4)){ ADD192(x4, x4, x4); SUB192(x4, q4, x4); }else{ ADD192(x4, x4, x4); }
		if(CMPUGT192(x5, qhalf5)){ ADD192(x5, x5, x5); SUB192(x5, q5, x5); }else{ ADD192(x5, x5, x5); }
		if(CMPUGT192(x6, qhalf6)){ ADD192(x6, x6, x6); SUB192(x6, q6, x6); }else{ ADD192(x6, x6, x6); }
		if(CMPUGT192(x7, qhalf7)){ ADD192(x7, x7, x7); SUB192(x7, q7, x7); }else{ ADD192(x7, x7, x7); }
	}

#if FAC_DEBUG
	if(CMPULT192(q0, x0)){ sprintf(char_buf, "twopmodq192_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint192_base10_char(str0, x0)], &str1[convert_uint192_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q1, x1)){ sprintf(char_buf, "twopmodq192_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint192_base10_char(str0, x1)], &str1[convert_uint192_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q2, x2)){ sprintf(char_buf, "twopmodq192_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint192_base10_char(str0, x2)], &str1[convert_uint192_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q3, x3)){ sprintf(char_buf, "twopmodq192_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint192_base10_char(str0, x3)], &str1[convert_uint192_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q4, x4)){ sprintf(char_buf, "twopmodq192_q8 : (x4 = %s) >= (q4 = %s)", &str0[convert_uint192_base10_char(str0, x4)], &str1[convert_uint192_base10_char(str1, q4)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q5, x5)){ sprintf(char_buf, "twopmodq192_q8 : (x5 = %s) >= (q5 = %s)", &str0[convert_uint192_base10_char(str0, x5)], &str1[convert_uint192_base10_char(str1, q5)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q6, x6)){ sprintf(char_buf, "twopmodq192_q8 : (x6 = %s) >= (q6 = %s)", &str0[convert_uint192_base10_char(str0, x6)], &str1[convert_uint192_base10_char(str1, q6)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT192(q7, x7)){ sprintf(char_buf, "twopmodq192_q8 : (x7 = %s) >= (q7 = %s)", &str0[convert_uint192_base10_char(str0, x7)], &str1[convert_uint192_base10_char(str1, q7)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
	/*...x^2 mod q is returned in x. */
	#ifdef PIPELINE_MUL192
		SQR_LOHI192_q8(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3
		, x4, lo4, hi4
		, x5, lo5, hi5
		, x6, lo6, hi6
		, x7, lo7, hi7);

		MULL192_q8(
		  lo0, qinv0, lo0
		, lo1, qinv1, lo1
		, lo2, qinv2, lo2
		, lo3, qinv3, lo3
		, lo4, qinv4, lo4
		, lo5, qinv5, lo5
		, lo6, qinv6, lo6
		, lo7, qinv7, lo7);

		MULH192_q8(
		  q0, lo0, lo0
		, q1, lo1, lo1
		, q2, lo2, lo2
		, q3, lo3, lo3
		, q4, lo4, lo4
		, q5, lo5, lo5
		, q6, lo6, lo6
		, q7, lo7, lo7);
	#else
	/* Strangely, on alpha/ia64 nopipe is better than the pipelined ..._q8 version: */
		SQR_LOHI192(x0,lo0,hi0);	MULL192(lo0,qinv0,lo0);		MULH192(q0,lo0,lo0);
		SQR_LOHI192(x1,lo1,hi1);	MULL192(lo1,qinv1,lo1);		MULH192(q1,lo1,lo1);
		SQR_LOHI192(x2,lo2,hi2);	MULL192(lo2,qinv2,lo2);		MULH192(q2,lo2,lo2);
		SQR_LOHI192(x3,lo3,hi3);	MULL192(lo3,qinv3,lo3);		MULH192(q3,lo3,lo3);
		SQR_LOHI192(x4,lo4,hi4);	MULL192(lo4,qinv4,lo4);		MULH192(q4,lo4,lo4);
		SQR_LOHI192(x5,lo5,hi5);	MULL192(lo5,qinv5,lo5);		MULH192(q5,lo5,lo5);
		SQR_LOHI192(x6,lo6,hi6);	MULL192(lo6,qinv6,lo6);		MULH192(q6,lo6,lo6);
		SQR_LOHI192(x7,lo7,hi7);	MULL192(lo7,qinv7,lo7);		MULH192(q7,lo7,lo7);
	#endif

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT192(hi0, lo0)) { SUB192(q0, lo0, lo0);	ADD192(lo0, hi0, x0); } else { SUB192(hi0, lo0, x0); }
		if(CMPULT192(hi1, lo1)) { SUB192(q1, lo1, lo1);	ADD192(lo1, hi1, x1); } else { SUB192(hi1, lo1, x1); }
		if(CMPULT192(hi2, lo2)) { SUB192(q2, lo2, lo2);	ADD192(lo2, hi2, x2); } else { SUB192(hi2, lo2, x2); }
		if(CMPULT192(hi3, lo3)) { SUB192(q3, lo3, lo3);	ADD192(lo3, hi3, x3); } else { SUB192(hi3, lo3, x3); }
		if(CMPULT192(hi4, lo4)) { SUB192(q4, lo4, lo4);	ADD192(lo4, hi4, x4); } else { SUB192(hi4, lo4, x4); }
		if(CMPULT192(hi5, lo5)) { SUB192(q5, lo5, lo5);	ADD192(lo5, hi5, x5); } else { SUB192(hi5, lo5, x5); }
		if(CMPULT192(hi6, lo6)) { SUB192(q6, lo6, lo6);	ADD192(lo6, hi6, x6); } else { SUB192(hi6, lo6, x6); }
		if(CMPULT192(hi7, lo7)) { SUB192(q7, lo7, lo7);	ADD192(lo7, hi7, x7); } else { SUB192(hi7, lo7, x7); }

		if(TEST_BIT192(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT192(x0,q0), "twopmodq192_q8 : CMPULT192(x0,q0)");
			ASSERT(HERE, CMPULT192(x1,q1), "twopmodq192_q8 : CMPULT192(x1,q1)");
			ASSERT(HERE, CMPULT192(x2,q2), "twopmodq192_q8 : CMPULT192(x2,q2)");
			ASSERT(HERE, CMPULT192(x3,q3), "twopmodq192_q8 : CMPULT192(x3,q3)");
			ASSERT(HERE, CMPULT192(x4,q4), "twopmodq192_q8 : CMPULT192(x4,q4)");
			ASSERT(HERE, CMPULT192(x5,q5), "twopmodq192_q8 : CMPULT192(x5,q5)");
			ASSERT(HERE, CMPULT192(x6,q6), "twopmodq192_q8 : CMPULT192(x6,q6)");
			ASSERT(HERE, CMPULT192(x7,q7), "twopmodq192_q8 : CMPULT192(x7,q7)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
			if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
			if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
			if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
			if(CMPUGT192(x4, qhalf4)){ ADD192(x4, x4, x4); SUB192(x4, q4, x4); }else{ ADD192(x4, x4, x4); }
			if(CMPUGT192(x5, qhalf5)){ ADD192(x5, x5, x5); SUB192(x5, q5, x5); }else{ ADD192(x5, x5, x5); }
			if(CMPUGT192(x6, qhalf6)){ ADD192(x6, x6, x6); SUB192(x6, q6, x6); }else{ ADD192(x6, x6, x6); }
			if(CMPUGT192(x7, qhalf7)){ ADD192(x7, x7, x7); SUB192(x7, q7, x7); }else{ ADD192(x7, x7, x7); }
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD192(x0,x0,x0);	SUB192(x0,q0,x0);
	ADD192(x1,x1,x1);	SUB192(x1,q1,x1);
	ADD192(x2,x2,x2);	SUB192(x2,q2,x2);
	ADD192(x3,x3,x3);	SUB192(x3,q3,x3);
	ADD192(x4,x4,x4);	SUB192(x4,q4,x4);
	ADD192(x5,x5,x5);	SUB192(x5,q5,x5);
	ADD192(x6,x6,x6);	SUB192(x6,q6,x6);
	ADD192(x7,x7,x7);	SUB192(x7,q7,x7);

	/* Only do the full 192-bit (Xj== 1) check if the bottom 64 bits of Xj == 1: */
	r = 0;
	if(x0.d0 == 1) r += ((uint64)CMPEQ192(x0, ONE192) << 0);
	if(x1.d0 == 1) r += ((uint64)CMPEQ192(x1, ONE192) << 1);
	if(x2.d0 == 1) r += ((uint64)CMPEQ192(x2, ONE192) << 2);
	if(x3.d0 == 1) r += ((uint64)CMPEQ192(x3, ONE192) << 3);
	if(x4.d0 == 1) r += ((uint64)CMPEQ192(x4, ONE192) << 4);
	if(x5.d0 == 1) r += ((uint64)CMPEQ192(x5, ONE192) << 5);
	if(x6.d0 == 1) r += ((uint64)CMPEQ192(x6, ONE192) << 6);
	if(x7.d0 == 1) r += ((uint64)CMPEQ192(x7, ONE192) << 7);
	return(r);
}

