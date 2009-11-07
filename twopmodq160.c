#include "factor.h"

#if FAC_DEBUG
	char char_buf[1024], str0[64], str1[64];
#endif

/***********************************************************************************/
/***160-BIT INPUTS *****************************************************************/
/***********************************************************************************/
/*
Function to find 2^(-p) mod q, where p and q are both 192-bit unsigned integers,
and both p and q are < 2^160, i.e. their upper 32 bits are zero.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^160 (i.e. our MODQ
operation effects multiply modulo 2^160).

The key 3-operation sequence here is as follows:

	SQR_LOHI160(x,lo,hi);
	MULL160(lo,qinv,lo);
	MULH160(q,lo,t);
*/
uint160 twopmodq160(uint160 p, uint160 q)
{
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint192_base10_char(char_buf, p)], "0");
	uint160 y,z;
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead8, lo64;
	uint160 qhalf, qinv, x, lo, hi;
	static uint160 psave = {0ull,0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq160:\n");
#endif

	RSHIFT_FAST160(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || !CMPEQ160(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		x.d0 = 160; x.d1 = x.d2 = 0;
		ADD160(p, x, pshift);
#if FAC_DEBUG
	if(dbg)
	{
		printf("p = %s\n", &char_buf[convert_uint192_base10_char(char_buf, p     )]);
		printf("p+= %s\n", &char_buf[convert_uint192_base10_char(char_buf, pshift)]);
	}
#endif
	/*
	find number of leading zeros in p, use it to find the position of the leftmost ones bit,
	and subtract 7 or 8 to account for the fact that we can do the powering for the leftmost
	7 or 8 bits (depending on whether the leftmost 8 > 159 or not) via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering: */
		if(pshift.d2)
		{
			j = leadz64(pshift.d2);	/* Remember, pshift is stored in a 192-bit... */
		#if FAC_DEBUG
			ASSERT(HERE, j >= 32,"twopmodq160: j >= 32");
		#endif
			/* Extract leftmost 8 bits of pshift (if > 159, use the leftmost 7) and subtract from 160: */
			lead8 = (((pshift.d2<<j) + (pshift.d1>>(64-j))) >> 56);	/* lead8 in [128, 255] */
			if(lead8 > 159)
			{
				lead8 >>= 1;	/* lead8 in [80, 128] */
				start_index = 192-j-7;	/* Use only the leftmost 7 bits */
			}
			else
				start_index = 192-j-8;
		}
		else if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			/* Extract leftmost 8 bits of pshift (if > 159, use the leftmost 7) and subtract from 160: */
			lead8 = (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 56);
			if(lead8 > 159)
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
			/* Extract leftmost 8 bits of pshift (if > 159, use the leftmost 7) and subtract from 160: */
			lead8 = ((pshift.d0<<j) >> 56);
			if(lead8 > 159)
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
		zshift = 159 - lead8;	/* lead8 in [80, 128], so zshift in [31, 79] */
		zshift <<= 1;			/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift.d2 = ~pshift.d2;	pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
		pshift.d2 &= 0x00000000ffffffff;	/* Don't really care about upper 32 bits of pshift, but
											if we don't zero them it'll trigger assertions in debug mode. */
	}

	/*
	!    Find modular inverse (mod 2^160) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q.d0 & (uint64)1) == 1, "twopmodq160 : (q.d0 & (uint64)1) == 1");
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

	/* Now that have bottom 64 bits of qinv, do two more Newton iterations using full 160-bit operands: */
	MULL160(q, qinv, x);
#if FAC_DEBUG
	MULL192(q, qinv, y);	y.d2 &= 0x00000000ffffffff;
	ASSERT(HERE, CMPEQ192(x, y), "twopmodq160: CMPEQ192(x, y)");
	SUB160 (TWO160, y, y);
	MULL192(y, qinv, y);	y.d2 &= 0x00000000ffffffff;
#endif
	SUB160 (TWO160, x, x);
	MULL160(qinv, x, qinv);
#if FAC_DEBUG
	ASSERT(HERE, CMPEQ192(qinv, y), "twopmodq160: CMPEQ192(qinv, y)");
#endif

	MULL160(q, qinv, x);
	SUB160 (TWO160, x, x);
	MULL160(qinv, x, qinv);

#if FAC_DEBUG
	if(dbg)
	{
		printf("q    = %s\n", &char_buf[convert_uint192_base10_char(char_buf, q   )]);
		printf("qinv = %s\n", &char_buf[convert_uint192_base10_char(char_buf, qinv)]);
	}
#endif
	/* Since zstart is a power of two < 2^160, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL160(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("j = start_index - 1 = %u\n", j);
	if(dbg) printf("zshift  = %u\n", zshift);
#endif
	LSHIFT160(qinv, zshift, lo);
#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint192_base10_char(char_buf, lo)]);
if(p.d0==629)
printf("");
	LSHIFT_FAST192(lo,32,y);	/* y = lo*2^32 */
	MULH192(q,y,y);
#endif
	MULH160(q,lo,lo);
#if FAC_DEBUG
	ASSERT(HERE, CMPEQ192(lo,y), "twopmodq160: CMPEQ192(lo,y)");
	if(dbg) printf("q*lo/2^160 = %s\n", &char_buf[convert_uint192_base10_char(char_buf, lo)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB160(q, lo, x);

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint192_base10_char(char_buf, lo)]);
#endif
	if(TEST_BIT160(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT160(x,q), "twopmodq160 : CMPULT160(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT160(x, qhalf)){ ADD160(x, x, x); SUB160(x, q, x); }else{ ADD160(x, x, x); }
#if FAC_DEBUG
if(dbg) printf("2x= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
	}

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
#if FAC_DEBUG
	if(CMPULT192(q, x)){ sprintf(char_buf, "twopmodq160 : (x0 = %s) >= (q = %s)", &str0[convert_uint192_base10_char(str0, x)], &str1[convert_uint192_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI160(x,lo,hi);
	#if FAC_DEBUG
		SQR_LOHI192(x, y, z);
		LSHIFT_FAST192(z,32,z);	z.d0 += (y.d2 >> 32);	/* x^2/2^160 */
		y.d2 &= 0x00000000ffffffff;							/* x^2%2^160 */
		ASSERT(HERE, CMPEQ192(lo,y), "twopmodq160: SQR_LOHI160: CMPEQ192(lo,y)");
		ASSERT(HERE, CMPEQ192(hi,z), "twopmodq160: SQR_LOHI160: CMPEQ192(hi,z)");
		y = lo;
		MULL192(y, qinv, y);	y.d2 &= 0x00000000ffffffff;
	#endif
		MULL160(lo,qinv,lo);
	#if FAC_DEBUG
		ASSERT(HERE, CMPEQ192(lo,y), "twopmodq160: MULL160: CMPEQ192(lo,y)");
		LSHIFT_FAST192(lo,32,y);	/* y = lo*2^32 */
		MULH192(q,y,y);
	#endif
		MULH160(q,lo,lo);
	#if FAC_DEBUG
		ASSERT(HERE, CMPEQ192(lo,y), "twopmodq160: MULH160: CMPEQ192(lo,y)");
	#endif

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT160(hi, lo))
		{
			SUB160(q, lo, lo);
			ADD160(lo, hi, x);
		}
		else
		{
			SUB160(hi, lo, x);
		}
#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif

		if(TEST_BIT160(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT160(x,q), "twopmodq160 : CMPULT160(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT160(x, qhalf)){ ADD160(x, x, x); SUB160(x, q, x); }else{ ADD160(x, x, x); }
#if FAC_DEBUG
	if(dbg) printf("2x= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD160(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^159, so x + x cannot overflow. */
#if FAC_DEBUG
	if(dbg) printf("Final x = %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
	SUB160(x,q,x);
#if FAC_DEBUG
	if(dbg) printf("Final x-q=%s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif

	return x;
}

