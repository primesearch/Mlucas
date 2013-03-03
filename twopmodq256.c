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
#include "imul256_macro.h"

#if FAC_DEBUG
	char str0[64], str1[64];
#endif

/***********************************************************************************/
/***256-BIT INPUTS *****************************************************************/
/***********************************************************************************/
/*
!...Function to find 2^(-p) mod q, where p and q are both 256-bit unsigned integers.
*/
uint256 twopmodq256(uint256 p, uint256 q)
{
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint256_base10_char(char_buf, p)], "0");
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead8, lo64;
	uint256 qhalf, qinv, x, lo, hi;
	static uint256 psave = {0ull, 0ull, 0ull, 0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq256:\n");
#endif

	RSHIFT_FAST256(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || !CMPEQ256(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		x.d0 = 256; x.d1 = x.d2 = x.d3 = 0;
		ADD256(p, x, pshift);
#if FAC_DEBUG
	if(dbg)
	{
		printf("p = %s\n", &char_buf[convert_uint256_base10_char(char_buf, p     )]);
		printf("p+= %s\n", &char_buf[convert_uint256_base10_char(char_buf, pshift)]);
	}
#endif
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 8 to account for the fact that we can do the powering for the leftmost
	!    8 bits via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming
		the leftmost 8 bits have already been processed via a shift.

		Since 8 bits with leftmost bit = 1 is guaranteed
		to be in [128,255], the shift count here is in [0,127].
		That means that zstart < 2^127. Together with the fact that
		squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration.
		Every little bit counts (literally in this case :), right?
		*/
		if(pshift.d3)
		{
			j = leadz64(pshift.d3);
			/* Extract leftmost 8 bits of pshift and subtract from 256: */
			lead8 = (((pshift.d3<<j) + (pshift.d2>>(64-j))) >> 56);
			start_index = 256-j-8;
		}
		else if(pshift.d2)
		{
			j = leadz64(pshift.d2);
			/* Extract leftmost 8 bits of pshift and subtract from 192: */
			lead8 = (((pshift.d2<<j) + (pshift.d1>>(64-j))) >> 56);
			start_index = 192-j-8;
		}
		else if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			/* Extract leftmost 8 bits of pshift and subtract from 128: */
			lead8 = (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 56);
			start_index = 128-j-8;
		}
		else
		{
			j = leadz64(pshift.d0);
			/* Extract leftmost 8 bits of pshift (if > 191, use the leftmost 7) and subtract from 192: */
			lead8 = ((pshift.d0<<j) >> 56);
			start_index =  64-j-8;
		}

#if FAC_DEBUG
	if(dbg)	printf("lead8 = %u\n", (uint32)lead8);
#endif
		zshift = 255 - lead8;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift.d3 = ~pshift.d3;	pshift.d2 = ~pshift.d2;	pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/*
	!    Find modular inverse (mod 2^256) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q.d0 & (uint64)1) == 1, "twopmodq256 : (q.d0 & (uint64)1) == 1");
#endif
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv.d0 = (q.d0 + q.d0 + q.d0) ^ (uint64)2;	qinv.d3 = qinv.d2 = qinv.d1 = (uint64)0;

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

	/* Now that have bottom 64 bits of qinv, do two more Newton iterations using full 256-bit operands: */
	MULL256(q, qinv, x);
	SUB256 (TWO256, x, x);
	MULL256(qinv, x, qinv);

	MULL256(q, qinv, x);
	SUB256 (TWO256, x, x);
	MULL256(qinv, x, qinv);

#if FAC_DEBUG
	if(dbg)
	{
		printf("q    = %s\n", &char_buf[convert_uint256_base10_char(char_buf, q   )]);
		printf("qinv = %s\n", &char_buf[convert_uint256_base10_char(char_buf, qinv)]);
	}
#endif
	/* Since zstart is a power of two < 2^256, use a streamlined code sequence for the first iteration: */
	ASSERT(HERE, start_index>=2, "twopmodq256 : start_index < 2!");
	j = start_index-1;

	/* MULL256(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("j = start_index - 1 = %u\n", j);
	if(dbg) printf("zshift  = %u\n", zshift);
#endif
	LSHIFT256(qinv, zshift, lo);
#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint256_base10_char(char_buf, lo)]);
#endif
	MULH256(q,lo,lo);
#if FAC_DEBUG
	if(dbg) printf("q*lo/2^256 = %s\n", &char_buf[convert_uint256_base10_char(char_buf, lo)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB256(q, lo, x);

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint256_base10_char(char_buf, lo)]);
#endif
	if(TEST_BIT256(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT256(x,q), "twopmodq256 : CMPULT256(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT256(x, qhalf)){ ADD256(x, x, x); SUB256(x, q, x); }else{ ADD256(x, x, x); }
#if FAC_DEBUG
if(dbg) printf("2x= %s\n", &char_buf[convert_uint256_base10_char(char_buf, x)]);
#endif
	}

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint256_base10_char(char_buf, x)]);
#endif
#if FAC_DEBUG
	if(CMPULT256(q, x)){ sprintf(char_buf, "twopmodq256 : (x0 = %s) >= (q = %s)", &str0[convert_uint256_base10_char(str0, x)], &str1[convert_uint256_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI256(x,lo,hi);
		MULL256(lo,qinv,lo);
		MULH256(q,lo,lo);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT256(hi, lo))
		{
			SUB256(q, lo, lo);
			ADD256(lo, hi, x);
		}
		else
		{
			SUB256(hi, lo, x);
		}
#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint256_base10_char(char_buf, x)]);
#endif

		if(TEST_BIT256(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT256(x,q), "twopmodq256 : CMPULT256(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT256(x, qhalf)){ ADD256(x, x, x); SUB256(x, q, x); }else{ ADD256(x, x, x); }
#if FAC_DEBUG
	if(dbg) printf("2x= %s\n", &char_buf[convert_uint256_base10_char(char_buf, x)]);
#endif
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD256(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^255, so x + x cannot overflow. */
#if FAC_DEBUG
	if(dbg) printf("Final x = %s\n", &char_buf[convert_uint256_base10_char(char_buf, x)]);
#endif
	SUB256(x,q,x);
#if FAC_DEBUG
	if(dbg) printf("Final x-q=%s\n", &char_buf[convert_uint256_base10_char(char_buf, x)]);
#endif

	return x;
}

