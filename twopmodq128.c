#include "factor.h"

#if FAC_DEBUG
	char char_buf[1024], str0[64], str1[64];
#endif

/***********************************************************************************/
/***128-BIT INPUTS *****************************************************************/
/***********************************************************************************/
/*
Function to find 2^(-p) mod q, where p is a 64-bit and q a 128-bit unsigned integer.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^128 (i.e. our MODQ
operation effects multiply modulo 2^128).

The key 3-operation sequence here is as follows:

	SQR_LOHI128(x,lo,hi);	// Input   x has 128 bits; lo/hi have 128 bits
	MULL128(lo,qinv,lo);	// Inputs lo & qinv, and output (overwrites lo) have 128 bits
	MULH128(q,lo,lo);	// Inputs  q &   lo, and output (overwrites lo) have 128 bits
*/
uint128 twopmodq128(uint64 p, uint128 q)
{
#if FAC_DEBUG
	int dbg = (p == 0);
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64;
	uint128 qhalf, qinv, x, lo, hi;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq128:\n");
#endif

	RSHIFT_FAST128(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

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

	/*
	!    Find modular inverse (mod 2^128) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q.d0 & (uint64)1) == 1, "twopmodq128 : (q.d0 & (uint64)1) == 1");
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
		lo64 = q.d0*qinv.d0;
		qinv.d0 = qinv.d0*((uint64)2 - lo64);
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
	MULH64(q.d0, qinv.d0, lo64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + lo64);
#endif

#if FAC_DEBUG
	ASSERT(HERE, qinv.d1 == x.d1 && qinv.d0 == x.d0, "twopmodq128 : qinv.d1 == x.d1 && qinv.d0 == x.d0");
	if(dbg) printf("q    = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q   )]);
	if(dbg) printf("qinv = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv)]);
#endif

	/* Since zstart is a power of two < 2^128, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL128(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("zshift  = %u\n", zshift);
#endif
	LSHIFT128(qinv, zshift, lo);

#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	MULH128(q,lo,lo);

#if FAC_DEBUG
	if(dbg) printf("q*lo/2^128 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB128(q, lo, x);	/* Put the result in lo (rather than x), to ease overflow check below */

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

	if((pshift >> j) & (uint64)1)
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT128(x, q), "twopmodq128 : CMPULT128(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT128(x, qhalf)){ ADD128(x, x, x); SUB128(x, q, x); }else{ ADD128(x, x, x); }
	}

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
#if FAC_DEBUG
	if(CMPULT128(q, x)){ sprintf(char_buf, "twopmodq128 : (x0 = %s) >= (q = %s)", &str0[convert_uint128_base10_char(str0, x)], &str1[convert_uint128_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI128(x,lo,hi);
		MULL128(lo,qinv,lo);
		MULH128(q,lo,lo);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT128(hi, lo))
		{
			SUB128(q, lo, lo);
			ADD128(lo, hi, x);
		}
		else
		{
			SUB128(hi, lo, x);
		}

#if FAC_DEBUG
if(dbg)printf("j = %2d, x = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT128(x, q), "twopmodq128 : CMPULT128(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT128(x, qhalf)){ ADD128(x, x, x); SUB128(x, q, x); }else{ ADD128(x, x, x); }

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
	ADD128(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^127, so x + x cannot overflow. */
	SUB128(x,q,x);

#if FAC_DEBUG
if(dbg)printf("x0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

	return x;
}

/*
!...Function to find 2^(-p) mod q, where p and q are both 128-bit unsigned integers.
*/
uint128 twopmodq128x2(uint64 *p, uint128 q)
{
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_mi64_base10_char(char_buf, p, 2)], "0");
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64;
	uint128 qhalf, qinv, x, lo, hi;
	static uint128 psave = {0ull, 0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq128x2:\n");
#endif

	RSHIFT_FAST128(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || !mi64_cmpeq(p, (uint64*)&psave, 2))
	{
		first_entry = FALSE;
		psave.d0 = p[0];	psave.d1 = p[1];
		pshift.d0 = p[0] + 128;	pshift.d1 = p[1] + (pshift.d0 < 128);
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 7 to account for the fact that we can do the powering for the leftmost
	!    7 bits via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming
		the leftmost 7 bits have already been processed via a shift.

		Since 7 bits with leftmost bit = 1 is guaranteed
		to be in [64,127], the shift count here is in [0, 63].
		That means that zstart < 2^64. Together with the fact that
		squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration.
		Every little bit counts (literally in this case :), right?
		*/
		if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			start_index = 128-j-7;
			/* Extract leftmost 7 bits of pshift and subtract from 127: */
			zshift = 127 - (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 57);
		}
		else
		{
			j = leadz64(pshift.d0);
			start_index =  64-j-7;
			zshift = 127 - ((pshift.d0<<j) >> 57);
		}


		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/*
	!    Find modular inverse (mod 2^128) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q.d0 & (uint64)1) == 1, "twopmodq128x2 : (q.d0 & (uint64)1) == 1");
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
		lo64 = q.d0*qinv.d0;
		qinv.d0 = qinv.d0*((uint64)2 - lo64);
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
	MULH64(q.d0, qinv.d0, lo64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + lo64);
#endif

#if FAC_DEBUG
	ASSERT(HERE, qinv.d1 == x.d1 && qinv.d0 == x.d0, "twopmodq128x2 : qinv.d1 == x.d1 && qinv.d0 == x.d0");
	if(dbg) printf("q    = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q   )]);
	if(dbg) printf("qinv = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv)]);
#endif
	/* Since zstart is a power of two < 2^128, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL128(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("zshift  = %u\n", zshift);
#endif
	LSHIFT128(qinv, zshift, lo);

#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	MULH128(q,lo,lo);

#if FAC_DEBUG
	if(dbg) printf("q*lo/2^128 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB128(q, lo, x);

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

	if(TEST_BIT128(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT128(x,q), "twopmodq128x2 : CMPULT128(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT128(x, qhalf)){ ADD128(x, x, x); SUB128(x, q, x); }else{ ADD128(x, x, x); }
	}

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
#if FAC_DEBUG
	if(CMPULT128(q, x)){ sprintf(char_buf, "twopmodq128x2 : (x0 = %s) >= (q = %s)", &str0[convert_uint128_base10_char(str0, x)], &str1[convert_uint128_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI128(x,lo,hi);
		MULL128(lo,qinv,lo);
		MULH128(q,lo,lo);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT128(hi, lo))
		{
			SUB128(q, lo, lo);
			ADD128(lo, hi, x);
		}
		else
		{
			SUB128(hi, lo, x);
		}
#if FAC_DEBUG
	if(dbg) printf("x = %s", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

		if(TEST_BIT128(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT128(x,q), "twopmodq128x2 : CMPULT128(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT128(x, qhalf)){ ADD128(x, x, x); SUB128(x, q, x); }else{ ADD128(x, x, x); }
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
	ADD128(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^127, so x + x cannot overflow. */
#if FAC_DEBUG
	if(dbg) printf("Final x = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
	SUB128(x,q,x);
#if FAC_DEBUG
	if(dbg) printf("Final x*= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

	return x;
}

/*** 4-trial-factor version ***/
uint64 twopmodq128_q4(uint128 p, uint128 q0, uint128 q1, uint128 q2, uint128 q3)
{
	 int32 j;
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, r;
	uint128 qinv0, qinv1, qinv2, qinv3
		, qhalf0, qhalf1, qhalf2, qhalf3
		, x0, x1, x2, x3
		, lo0, lo1, lo2, lo3
		, hi0, hi1, hi2, hi3;
	static uint128 psave = {0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	RSHIFT_FAST128(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST128(q1, 1, qhalf1);
	RSHIFT_FAST128(q2, 1, qhalf2);
	RSHIFT_FAST128(q3, 1, qhalf3);

	if(first_entry || !CMPEQ128(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		pshift.d0 = p.d0 + 128;	pshift.d1 = p.d1 + (pshift.d0 < 128);

		if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			start_index = 128-j-7;
			/* Extract leftmost 7 bits of pshift and subtract from 127: */
			zshift = 127 - (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 57);
		}
		else
		{
			j = leadz64(pshift.d0);
			start_index =  64-j-7;
			zshift = 127 - ((pshift.d0<<j) >> 57);
		}

		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q0.d0 & (uint64)1) == 1, "twopmodq128_q4 : (q0.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q1.d0 & (uint64)1) == 1, "twopmodq128_q4 : (q1.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q2.d0 & (uint64)1) == 1, "twopmodq128_q4 : (q2.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q3.d0 & (uint64)1) == 1, "twopmodq128_q4 : (q3.d0 & (uint64)1) == 1");
#endif
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d1 = (uint64)0;

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

	/* qinv has 128 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + __MULH64(q2.d0, qinv2.d0));
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + __MULH64(q3.d0, qinv3.d0));
#else
	MULH64(q0.d0, qinv0.d0, lo64_0);
	MULH64(q1.d0, qinv1.d0, lo64_1);
	MULH64(q2.d0, qinv2.d0, lo64_2);
	MULH64(q3.d0, qinv3.d0, lo64_3);

	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + lo64_0);
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + lo64_1);
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + lo64_2);
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + lo64_3);
#endif

	/* Since zstart is a power of two < 2^128, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL128(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT128(qinv0, zshift, lo0);
	LSHIFT128(qinv1, zshift, lo1);
	LSHIFT128(qinv2, zshift, lo2);
	LSHIFT128(qinv3, zshift, lo3);

	MULH128_q4(
	  q0, lo0, lo0
	, q1, lo1, lo1
	, q2, lo2, lo2
	, q3, lo3, lo3);

	/* hi = 0 in this instance, which simplifies things. */
	SUB128(q0, lo0, x0);
	SUB128(q1, lo1, x1);
	SUB128(q2, lo2, x2);
	SUB128(q3, lo3, x3);

	if(TEST_BIT128(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT128(x0, q0), "twopmodq128_q4 : CMPULT128(x0, q0)");
		ASSERT(HERE, CMPULT128(x1, q1), "twopmodq128_q4 : CMPULT128(x1, q1)");
		ASSERT(HERE, CMPULT128(x2, q2), "twopmodq128_q4 : CMPULT128(x2, q2)");
		ASSERT(HERE, CMPULT128(x3, q3), "twopmodq128_q4 : CMPULT128(x3, q3)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT128(x0, qhalf0)){ ADD128(x0, x0, x0); SUB128(x0, q0, x0); }else{ ADD128(x0, x0, x0); }
		if(CMPUGT128(x1, qhalf1)){ ADD128(x1, x1, x1); SUB128(x1, q1, x1); }else{ ADD128(x1, x1, x1); }
		if(CMPUGT128(x2, qhalf2)){ ADD128(x2, x2, x2); SUB128(x2, q2, x2); }else{ ADD128(x2, x2, x2); }
		if(CMPUGT128(x3, qhalf3)){ ADD128(x3, x3, x3); SUB128(x3, q3, x3); }else{ ADD128(x3, x3, x3); }
	}

#if FAC_DEBUG
	if(CMPULT128(q0, x0)){ sprintf(char_buf, "twopmodq128_q4 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint128_base10_char(str0, x0)], &str1[convert_uint128_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q1, x1)){ sprintf(char_buf, "twopmodq128_q4 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint128_base10_char(str0, x1)], &str1[convert_uint128_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q2, x2)){ sprintf(char_buf, "twopmodq128_q4 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint128_base10_char(str0, x2)], &str1[convert_uint128_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q3, x3)){ sprintf(char_buf, "twopmodq128_q4 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint128_base10_char(str0, x3)], &str1[convert_uint128_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
	#if THREE_OP128
		/* Fused version of all 3 of the above function calls. Surprisingly, on Alpha this was significantly slower
		than the 3-function version. */
		THREE_OP128_q4(
		  x0, qinv0, q0, lo0, hi0
		, x1, qinv1, q1, lo1, hi1
		, x2, qinv2, q2, lo2, hi2
		, x3, qinv3, q3, lo3, hi3);
	#else
		/* Haven't gotten IA64 version of this working properly yet:
		SQR_LOHI_INPLACE128_q4(
		  x0, hi0
		, x1, hi1
		, x2, hi2
		, x3, hi3);
		*/
		SQR_LOHI128_q4(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3);
	/*********** TRY 127-BIT VERSION: *************
		SQR_LOHI127_q4(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3);
	**********************************************/
		/* For unknown reasons, the 8-operand version of MULL128 was slower than one-at-a-time. */

		MULL128_INPLACE_q4(
		  lo0, qinv0
		, lo1, qinv1
		, lo2, qinv2
		, lo3, qinv3);
		/*
		MULL128_q4(
		  lo0, qinv0, lo0
		, lo1, qinv1, lo1
		, lo2, qinv2, lo2
		, lo3, qinv3, lo3);
		//
		MULL128(lo0, qinv0, lo0);
		MULL128(lo1, qinv1, lo1);
		MULL128(lo2, qinv2, lo2);
		MULL128(lo3, qinv3, lo3);
		*/

		MULH128_q4(
		  lo0, q0, lo0
		, lo1, q1, lo1
		, lo2, q2, lo2
		, lo3, q3, lo3);
		/*
		MULH128_q4(
		  lo0, q0, lo0
		, lo1, q1, lo1
		, lo2, q2, lo2
		, lo3, q3, lo3);
		//
		MULH128(lo0, q0, lo0);
		MULH128(lo1, q1, lo1);
		MULH128(lo2, q2, lo2);
		MULH128(lo3, q3, lo3);
		*/
	#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT128(hi0, lo0)) { SUB128(q0, lo0, lo0);	ADD128(lo0, hi0, x0); } else { SUB128(hi0, lo0, x0); }
		if(CMPULT128(hi1, lo1)) { SUB128(q1, lo1, lo1);	ADD128(lo1, hi1, x1); } else { SUB128(hi1, lo1, x1); }
		if(CMPULT128(hi2, lo2)) { SUB128(q2, lo2, lo2);	ADD128(lo2, hi2, x2); } else { SUB128(hi2, lo2, x2); }
		if(CMPULT128(hi3, lo3)) { SUB128(q3, lo3, lo3);	ADD128(lo3, hi3, x3); } else { SUB128(hi3, lo3, x3); }

		if(TEST_BIT128(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT128(x0,q0), "twopmodq128_q4 : CMPULT128(x0,q0)");
			ASSERT(HERE, CMPULT128(x1,q1), "twopmodq128_q4 : CMPULT128(x1,q1)");
			ASSERT(HERE, CMPULT128(x2,q2), "twopmodq128_q4 : CMPULT128(x2,q2)");
			ASSERT(HERE, CMPULT128(x3,q3), "twopmodq128_q4 : CMPULT128(x3,q3)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT128(x0, qhalf0)){ ADD128(x0, x0, x0); SUB128(x0, q0, x0); }else{ ADD128(x0, x0, x0); }
			if(CMPUGT128(x1, qhalf1)){ ADD128(x1, x1, x1); SUB128(x1, q1, x1); }else{ ADD128(x1, x1, x1); }
			if(CMPUGT128(x2, qhalf2)){ ADD128(x2, x2, x2); SUB128(x2, q2, x2); }else{ ADD128(x2, x2, x2); }
			if(CMPUGT128(x3, qhalf3)){ ADD128(x3, x3, x3); SUB128(x3, q3, x3); }else{ ADD128(x3, x3, x3); }
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD128(x0 ,x0, x0);
	ADD128(x1 ,x1, x1);
	ADD128(x2 ,x2, x2);
	ADD128(x3 ,x3, x3);

	SUB128(x0, q0, x0);
	SUB128(x1, q1, x1);
	SUB128(x2, q2, x2);
	SUB128(x3, q3, x3);

	r = 0;
	if(CMPEQ128(x0, ONE128)) r +=  1;
	if(CMPEQ128(x1, ONE128)) r +=  2;
	if(CMPEQ128(x2, ONE128)) r +=  4;
	if(CMPEQ128(x3, ONE128)) r +=  8;
	return(r);
}

/*** 8-trial-factor version ***/
uint64 twopmodq128_q8(uint128 p, uint128 q0, uint128 q1, uint128 q2, uint128 q3, uint128 q4, uint128 q5, uint128 q6, uint128 q7)
{
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint128_base10_char(char_buf, p)], "0");
#endif
	 int32 j;
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, lo64_4, lo64_5, lo64_6, lo64_7, r;
	uint128 qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, qhalf0, qhalf1, qhalf2, qhalf3, qhalf4, qhalf5, qhalf6, qhalf7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
		, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7;
	static uint128 psave = {0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq128_q8:\n");
#endif
	RSHIFT_FAST128(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST128(q1, 1, qhalf1);
	RSHIFT_FAST128(q2, 1, qhalf2);
	RSHIFT_FAST128(q3, 1, qhalf3);
	RSHIFT_FAST128(q4, 1, qhalf4);
	RSHIFT_FAST128(q5, 1, qhalf5);
	RSHIFT_FAST128(q6, 1, qhalf6);
	RSHIFT_FAST128(q7, 1, qhalf7);

	if(first_entry || !CMPEQ128(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		pshift.d0 = p.d0 + 128;	pshift.d1 = p.d1 + (pshift.d0 < 128);

		if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			start_index = 128-j-7;
			/* Extract leftmost 7 bits of pshift and subtract from 127: */
			zshift = 127 - (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 57);
		}
		else
		{
			j = leadz64(pshift.d0);
			start_index =  64-j-7;
			zshift = 127 - ((pshift.d0<<j) >> 57);
		}

		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/*
	!    Find modular inverse (mod 2^128) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q0.d0 & (uint64)1) == 1, "twopmodq128_q8 : (q0.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q1.d0 & (uint64)1) == 1, "twopmodq128_q8 : (q1.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q2.d0 & (uint64)1) == 1, "twopmodq128_q8 : (q2.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q3.d0 & (uint64)1) == 1, "twopmodq128_q8 : (q3.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q4.d0 & (uint64)1) == 1, "twopmodq128_q8 : (q4.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q5.d0 & (uint64)1) == 1, "twopmodq128_q8 : (q5.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q6.d0 & (uint64)1) == 1, "twopmodq128_q8 : (q6.d0 & (uint64)1) == 1");
	ASSERT(HERE, (q7.d0 & (uint64)1) == 1, "twopmodq128_q8 : (q7.d0 & (uint64)1) == 1");
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
	MULH64(q0.d0, qinv0.d0, lo64_0);
	MULH64(q1.d0, qinv1.d0, lo64_1);
	MULH64(q2.d0, qinv2.d0, lo64_2);
	MULH64(q3.d0, qinv3.d0, lo64_3);
	MULH64(q4.d0, qinv4.d0, lo64_4);
	MULH64(q5.d0, qinv5.d0, lo64_5);
	MULH64(q6.d0, qinv6.d0, lo64_6);
	MULH64(q7.d0, qinv7.d0, lo64_7);

	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + lo64_0);
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + lo64_1);
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + lo64_2);
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + lo64_3);
	qinv4.d1 = -qinv4.d0*(q4.d1*qinv4.d0 + lo64_4);
	qinv5.d1 = -qinv5.d0*(q5.d1*qinv5.d0 + lo64_5);
	qinv6.d1 = -qinv6.d0*(q6.d1*qinv6.d0 + lo64_6);
	qinv7.d1 = -qinv7.d0*(q7.d1*qinv7.d0 + lo64_7);
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

	MULH128_q8(
	  q0, lo0, lo0
	, q1, lo1, lo1
	, q2, lo2, lo2
	, q3, lo3, lo3
	, q4, lo4, lo4
	, q5, lo5, lo5
	, q6, lo6, lo6
	, q7, lo7, lo7);

	/* hi = 0 in this instance, which simplifies things. */
	/* Put the result in lo (rather than x), to ease overflow check below */
	SUB128(q0, lo0, x0);
	SUB128(q1, lo1, x1);
	SUB128(q2, lo2, x2);
	SUB128(q3, lo3, x3);
	SUB128(q4, lo4, x4);
	SUB128(q5, lo5, x5);
	SUB128(q6, lo6, x6);
	SUB128(q7, lo7, x7);

	if(TEST_BIT128(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT128(x0, q0), "twopmodq128_q8 : CMPULT128(x0, q0)");
		ASSERT(HERE, CMPULT128(x1, q1), "twopmodq128_q8 : CMPULT128(x1, q1)");
		ASSERT(HERE, CMPULT128(x2, q2), "twopmodq128_q8 : CMPULT128(x2, q2)");
		ASSERT(HERE, CMPULT128(x3, q3), "twopmodq128_q8 : CMPULT128(x3, q3)");
		ASSERT(HERE, CMPULT128(x4, q4), "twopmodq128_q8 : CMPULT128(x4, q4)");
		ASSERT(HERE, CMPULT128(x5, q5), "twopmodq128_q8 : CMPULT128(x5, q5)");
		ASSERT(HERE, CMPULT128(x6, q6), "twopmodq128_q8 : CMPULT128(x6, q6)");
		ASSERT(HERE, CMPULT128(x7, q7), "twopmodq128_q8 : CMPULT128(x7, q7)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT128(x0, qhalf0)){ ADD128(x0, x0, x0); SUB128(x0, q0, x0); }else{ ADD128(x0, x0, x0); }
		if(CMPUGT128(x1, qhalf1)){ ADD128(x1, x1, x1); SUB128(x1, q1, x1); }else{ ADD128(x1, x1, x1); }
		if(CMPUGT128(x2, qhalf2)){ ADD128(x2, x2, x2); SUB128(x2, q2, x2); }else{ ADD128(x2, x2, x2); }
		if(CMPUGT128(x3, qhalf3)){ ADD128(x3, x3, x3); SUB128(x3, q3, x3); }else{ ADD128(x3, x3, x3); }
		if(CMPUGT128(x4, qhalf4)){ ADD128(x4, x4, x4); SUB128(x4, q4, x4); }else{ ADD128(x4, x4, x4); }
		if(CMPUGT128(x5, qhalf5)){ ADD128(x5, x5, x5); SUB128(x5, q5, x5); }else{ ADD128(x5, x5, x5); }
		if(CMPUGT128(x6, qhalf6)){ ADD128(x6, x6, x6); SUB128(x6, q6, x6); }else{ ADD128(x6, x6, x6); }
		if(CMPUGT128(x7, qhalf7)){ ADD128(x7, x7, x7); SUB128(x7, q7, x7); }else{ ADD128(x7, x7, x7); }
	}

#if FAC_DEBUG
	if(CMPULT128(q0, x0)){ sprintf(char_buf, "twopmodq128_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint128_base10_char(str0, x0)], &str1[convert_uint128_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q1, x1)){ sprintf(char_buf, "twopmodq128_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint128_base10_char(str0, x1)], &str1[convert_uint128_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q2, x2)){ sprintf(char_buf, "twopmodq128_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint128_base10_char(str0, x2)], &str1[convert_uint128_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q3, x3)){ sprintf(char_buf, "twopmodq128_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint128_base10_char(str0, x3)], &str1[convert_uint128_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q4, x4)){ sprintf(char_buf, "twopmodq128_q8 : (x4 = %s) >= (q4 = %s)", &str0[convert_uint128_base10_char(str0, x4)], &str1[convert_uint128_base10_char(str1, q4)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q5, x5)){ sprintf(char_buf, "twopmodq128_q8 : (x5 = %s) >= (q5 = %s)", &str0[convert_uint128_base10_char(str0, x5)], &str1[convert_uint128_base10_char(str1, q5)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q6, x6)){ sprintf(char_buf, "twopmodq128_q8 : (x6 = %s) >= (q6 = %s)", &str0[convert_uint128_base10_char(str0, x6)], &str1[convert_uint128_base10_char(str1, q6)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(CMPULT128(q7, x7)){ sprintf(char_buf, "twopmodq128_q8 : (x7 = %s) >= (q7 = %s)", &str0[convert_uint128_base10_char(str0, x7)], &str1[convert_uint128_base10_char(str1, q7)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
#if FAC_DEBUG
if(dbg)printf("A: x = %20llu + 2^64* %20llu\n",x0.d0,x0.d1);
#endif
	#if THREE_OP128
		/* Fused version of all 3 of the above function calls. Surprisingly, on Alpha this was significantly slower
		than the 3-function version. */
		THREE_OP128_q8(
		  x0, qinv0, q0, lo0, hi0
		, x1, qinv1, q1, lo1, hi1
		, x2, qinv2, q2, lo2, hi2
		, x3, qinv3, q3, lo3, hi3
		, x4, qinv4, q4, lo4, hi4
		, x5, qinv5, q5, lo5, hi5
		, x6, qinv6, q6, lo6, hi6
		, x7, qinv7, q7, lo7, hi7);
	#else
		/* Haven't gotten IA64 version of this working properly yet:
		SQR_LOHI_INPLACE128_q8(
		  x0, hi0
		, x1, hi1
		, x2, hi2
		, x3, hi3
		, x4, hi4
		, x5, hi5
		, x6, hi6
		, x7, hi7);
		*/
	/*********** TRY 127-BIT VERSION: *************
		SQR_LOHI127_q8(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3
		, x4, lo4, hi4
		, x5, lo5, hi5
		, x6, lo6, hi6
		, x7, lo7, hi7);
	**********************************************/
		SQR_LOHI128_q8(
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
if(dbg)printf("B: h = %20llu + 2^64* %20llu\n",hi0.d0,hi0.d1);
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
		/*
		MULL128_q8(
		  lo0, qinv0, lo0
		, lo1, qinv1, lo1
		, lo2, qinv2, lo2
		, lo3, qinv3, lo3
		, lo4, qinv4, lo4
		, lo5, qinv5, lo5
		, lo6, qinv6, lo6
		, lo7, qinv7, lo7);
		//
		MULL128(lo0, qinv0, lo0);
		MULL128(lo1, qinv1, lo1);
		MULL128(lo2, qinv2, lo2);
		MULL128(lo3, qinv3, lo3);
		MULL128(lo4, qinv4, lo4);
		MULL128(lo5, qinv5, lo5);
		MULL128(lo6, qinv6, lo6);
		MULL128(lo7, qinv7, lo7);
		*/
#if FAC_DEBUG
if(dbg)printf("C: l = %20llu + 2^64* %20llu\n",lo0.d0,lo0.d1);
#endif

		MULH128_q8(
		  q0, lo0, lo0
		, q1, lo1, lo1
		, q2, lo2, lo2
		, q3, lo3, lo3
		, q4, lo4, lo4
		, q5, lo5, lo5
		, q6, lo6, lo6
		, q7, lo7, lo7);
		/*
		MULH128_q8(
		  lo0, q0, lo0
		, lo1, q1, lo1
		, lo2, q2, lo2
		, lo3, q3, lo3
		, lo4, q4, lo4
		, lo5, q5, lo5
		, lo6, q6, lo6
		, lo7, q7, lo7);
		//
		MULH128(lo0, q0, lo0);
		MULH128(lo1, q1, lo1);
		MULH128(lo2, q2, lo2);
		MULH128(lo3, q3, lo3);
		MULH128(lo4, q4, lo4);
		MULH128(lo5, q5, lo5);
		MULH128(lo6, q6, lo6);
		MULH128(lo7, q7, lo7);
		*/
#if FAC_DEBUG
if(dbg)printf("D: l = %20llu + 2^64* %20llu\n",lo0.d0,lo0.d1);
#endif
	#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT128(hi0, lo0)) { SUB128(q0, lo0, lo0);	ADD128(lo0, hi0, x0); } else { SUB128(hi0, lo0, x0); }
		if(CMPULT128(hi1, lo1)) { SUB128(q1, lo1, lo1);	ADD128(lo1, hi1, x1); } else { SUB128(hi1, lo1, x1); }
		if(CMPULT128(hi2, lo2)) { SUB128(q2, lo2, lo2);	ADD128(lo2, hi2, x2); } else { SUB128(hi2, lo2, x2); }
		if(CMPULT128(hi3, lo3)) { SUB128(q3, lo3, lo3);	ADD128(lo3, hi3, x3); } else { SUB128(hi3, lo3, x3); }
		if(CMPULT128(hi4, lo4)) { SUB128(q4, lo4, lo4);	ADD128(lo4, hi4, x4); } else { SUB128(hi4, lo4, x4); }
		if(CMPULT128(hi5, lo5)) { SUB128(q5, lo5, lo5);	ADD128(lo5, hi5, x5); } else { SUB128(hi5, lo5, x5); }
		if(CMPULT128(hi6, lo6)) { SUB128(q6, lo6, lo6);	ADD128(lo6, hi6, x6); } else { SUB128(hi6, lo6, x6); }
		if(CMPULT128(hi7, lo7)) { SUB128(q7, lo7, lo7);	ADD128(lo7, hi7, x7); } else { SUB128(hi7, lo7, x7); }
#if FAC_DEBUG
if(dbg)printf("j = %2d, Res = %20llu + 2^64* %20llu",j,x0.d0,x0.d1);
#endif

		if(TEST_BIT128(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT128(x0,q0), "twopmodq128_q8 : CMPULT128(x0,q0)");
			ASSERT(HERE, CMPULT128(x1,q1), "twopmodq128_q8 : CMPULT128(x1,q1)");
			ASSERT(HERE, CMPULT128(x2,q2), "twopmodq128_q8 : CMPULT128(x2,q2)");
			ASSERT(HERE, CMPULT128(x3,q3), "twopmodq128_q8 : CMPULT128(x3,q3)");
			ASSERT(HERE, CMPULT128(x4,q4), "twopmodq128_q8 : CMPULT128(x4,q4)");
			ASSERT(HERE, CMPULT128(x5,q5), "twopmodq128_q8 : CMPULT128(x5,q5)");
			ASSERT(HERE, CMPULT128(x6,q6), "twopmodq128_q8 : CMPULT128(x6,q6)");
			ASSERT(HERE, CMPULT128(x7,q7), "twopmodq128_q8 : CMPULT128(x7,q7)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT128(x0, qhalf0)){ ADD128(x0, x0, x0); SUB128(x0, q0, x0); }else{ ADD128(x0, x0, x0); }
			if(CMPUGT128(x1, qhalf1)){ ADD128(x1, x1, x1); SUB128(x1, q1, x1); }else{ ADD128(x1, x1, x1); }
			if(CMPUGT128(x2, qhalf2)){ ADD128(x2, x2, x2); SUB128(x2, q2, x2); }else{ ADD128(x2, x2, x2); }
			if(CMPUGT128(x3, qhalf3)){ ADD128(x3, x3, x3); SUB128(x3, q3, x3); }else{ ADD128(x3, x3, x3); }
			if(CMPUGT128(x4, qhalf4)){ ADD128(x4, x4, x4); SUB128(x4, q4, x4); }else{ ADD128(x4, x4, x4); }
			if(CMPUGT128(x5, qhalf5)){ ADD128(x5, x5, x5); SUB128(x5, q5, x5); }else{ ADD128(x5, x5, x5); }
			if(CMPUGT128(x6, qhalf6)){ ADD128(x6, x6, x6); SUB128(x6, q6, x6); }else{ ADD128(x6, x6, x6); }
			if(CMPUGT128(x7, qhalf7)){ ADD128(x7, x7, x7); SUB128(x7, q7, x7); }else{ ADD128(x7, x7, x7); }
#if FAC_DEBUG
if(dbg)printf(" *2 = %20llu + 2^64* %20llu",x0.d0,x0.d1);
#endif
		}
#if FAC_DEBUG
if(dbg)printf("\n");
#endif
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD128(x0 ,x0, x0);
	ADD128(x1 ,x1, x1);
	ADD128(x2 ,x2, x2);
	ADD128(x3 ,x3, x3);
	ADD128(x4 ,x4, x4);
	ADD128(x5 ,x5, x5);
	ADD128(x6 ,x6, x6);
	ADD128(x7 ,x7, x7);

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

	r = 0;
	if(CMPEQ128(x0, ONE128)) r +=  1;
	if(CMPEQ128(x1, ONE128)) r +=  2;
	if(CMPEQ128(x2, ONE128)) r +=  4;
	if(CMPEQ128(x3, ONE128)) r +=  8;
	if(CMPEQ128(x4, ONE128)) r += 16;
	if(CMPEQ128(x5, ONE128)) r += 32;
	if(CMPEQ128(x6, ONE128)) r += 64;
	if(CMPEQ128(x7, ONE128)) r +=128;
	return(r);
}

