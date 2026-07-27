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

#ifndef FAC_DEBUG
	#define FAC_DEBUG	0
#elif FAC_DEBUG < 0 || FAC_DEBUG > 2
	#error If FAC_DEBUG defined it must be assigned a value of 0 (off) or 1 (on).
#endif
#if FAC_DEBUG
	#warning FAC_DEBUG = 1: Enabling debug-printing.
	char char_buf[1024], str0[64], str1[64];
#endif

#define MONT_MUL128(__x,__y,__q,__qinv,__z)\
{\
	uint128 lo,hi;					\
	MUL_LOHI128(__x,__y,lo,hi);		\
	MULL128(__qinv,lo,lo);			\
	MULH128(__q,lo,lo);			\
	/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */\
	if(CMPULT128(hi, lo)) {	\
		SUB128(__q, lo, lo);	ADD128(lo, hi, __z);	\
	} else {	\
		SUB128(hi, lo, __z);	\
	}	\
};

/***********************************************************************************/
/***128-BIT INPUTS *****************************************************************/
/***********************************************************************************/
// Conventional positive-power version of twopmodq128, returns true mod:
uint128 twopmmodq128(uint128 p, uint128 q)
{
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint128_base10_char(char_buf, q)], "0");	// Replace "0" with "[desired decimal-form debug modulus]"
	int32 pow;
#endif
	 int32 j;	// j needs to be signed because of the LR binary exponentiation
	uint32 curr_bit, leadb, start_index, nshift;
	uint64 lo64;
	uint128 pshift, qhalf, qinv, x, lo,hi, rsqr;
	// char_buf is local, g_cstr is globally available:
#if FAC_DEBUG
	if(dbg) printf("twopmmodq128: computing 2^%s (mod %s)\n",&char_buf[convert_uint128_base10_char(char_buf,p)],&g_cstr[convert_uint128_base10_char(g_cstr,q)]);
#endif
	RSHIFT_FAST128(q, 1, qhalf);	/* = (q-1)/2, since q odd. */
	// If p <= 128, directly compute 2^p (mod q):
	if(p.d1 == 0 && p.d0 <= 128) {
		// Lshift (1 << j) to align with leading bit of q, then do (p - j) repeated mod-doublings:
		x.d0 = 1ull; x.d1 = 0ull;
		j = leadz128(q);	// q >= 2^(128 - j - 1)
		j = (128 - j - 1);
		if(j > p.d0) {
			LSHIFT128(x,(uint32)p.d0,x);
		} else {
			LSHIFT128(x,j,x);
		}
		for( ; j < p.d0; j++) {
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT128(x, qhalf)){ ADD128(x, x, x); SUB128(x, q, x); }else{ ADD128(x, x, x); }
		}
		return x;
	}
	// If get here, p > 128: set up for Montgomery-mul-based powering loop:
	nshift = trailz128(q);
	if(nshift) {
		x.d0 = (uint64)nshift; x.d1 = 0ull; SUB128(p,x,p);	// p >= nshift guaranteed here:
		RSHIFT128(q,nshift,q);	// Right-shift dividend by (nshift) bits; for 2^p this means subtracting nshift from p
	#if FAC_DEBUG
		if(dbg) printf("Removed power-of-2 from q: q' = (q >> %u) = %s\n",nshift,&char_buf[convert_uint128_base10_char(char_buf,q)]);
	#endif
	}
	// Extract leftmost 8 bits of (p - 128); if > 128, use leftmost 7 instead:
	x.d0 = 128ull; x.d1 = 0ull; SUB128(p,x,pshift); j = leadz128(pshift);
	LSHIFT128(pshift,j,x);	leadb = x.d1 >> 56;	// leadb = (pshift<<j) >> 57; no (pshift = ~pshift) step in positive-power algorithm!
	if(leadb > 128) {
		start_index = 128-7-j;
		leadb >>= 1;
	} else {
		start_index = 128-8-j;
	}
#if FAC_DEBUG
	if(dbg) {
		printf("leadb = %u\n",j);
		printf("pshift = p - %u = %s\n",128,&char_buf[convert_uint128_base10_char(char_buf,pshift)]);
	}
#endif
	// Find inverse (mod 2^128) of q; q must be odd for Montgomery-style modmul to work:
	ASSERT((q.d0 & (uint64)1) == 1, "twopmmodq128 : q must be odd for Montgomery-style modmul!");
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv.d0 = (q.d0 + q.d0 + q.d0) ^ (uint64)2;	qinv.d1 = (uint64)0;
	/* Compute qinv  = q^-1 (mod R = 2^128) via Newton iteration qinv = qinv*(2 - q*qinv), starting with
	5-bits-good approximation qinv_0 = 3*q ^ 2. Number-good-bits doubles each iteration until >= lg2(R) = 128:
	*/
	for(j = 0; j < 4; j++) 	{
		lo64 = q.d0*qinv.d0;
		qinv.d0 = qinv.d0*((uint64)2 - lo64);
	}
	/* qinv.d1 = 0 and q*qinv will == 1 mod 2^64 here; can take advantage of that fact to speed the iteration
	(i.e. q*qinv = (q.d1*2^64 + q.d0)*qinv.d0 = q.d1*qinv.d0*2^64 + q.d0*qinv.d0 == 1 mod 2^64,
	so do a MULL64(q.d1,qinv.d0) + MULH64(q.d0,qinv.d0) to get bits 64-127 of q*qinv = r, set lower half = 1,
	then simply negate upper half to get 2-q*qinv (mod 2^128), then use that lower half still = 1 to
	simplify the qinv*r MULL128 operation, i.e. qinv*r = qinv.d0*(r.d1*2^64 + 1) == (MULL64(qinv.d0*r.d1)*2^64 + qinv.d0) mod 2^128.
	This needs a total of just 3 MUL instructions and 2 ALUs, compared to 8 MULs and 8 ALUs for the original sequence.
	*/
	// qinv is 128 bits wide, but only the upper 64 get modified here:
#ifdef MUL_LOHI64_SUBROUTINE
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + __MULH64(q.d0, qinv.d0));
#else
	MULH64(q.d0, qinv.d0, lo64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + lo64);
#endif
	/* Initialize binary powering = R*x (mod q), where R = binary radix (2^128 here);
	at present don't care about optimizing this rarely-used function. */
	// First compute R^2 (mod q) in prep. for Mont-mul with initial seed:
	uint64 vtmp[5] = {0ull,0ull,0ull,0ull,1ull};	// R^2 = 2^256
	mi64_div_binary((const uint64*)vtmp, (const uint64*)&q, 5,2, 0x0, (uint32*)&j, (uint64*)&rsqr);

	// If leadb = 128, x = 2^128 = R, thus rsqr holds our desired starting value for x:
	if(leadb == 128)
		x = rsqr;
	else {
		x.d0 = 1ull; LSHIFT128(x,leadb,x);	// x <<= leadb;
		MONT_MUL128(x,rsqr, q,qinv, x);	// x*R (mod q) = MONT_MUL(x,R^2 (mod q),q,qinv)
 	}

#if FAC_DEBUG
	if(dbg) {
		printf("qinv = %s\n",&char_buf[convert_uint128_base10_char(char_buf,qinv)]);
		printf("leadb = %u, x0 = %s\n",leadb,&char_buf[convert_uint128_base10_char(char_buf,x)]);
		pow = leadb + 128;
		printf("Initial power = 2^(%u+128) = 2^%u mod q' = %s\n",leadb,pow,&char_buf[convert_uint128_base10_char(char_buf,x)]);
		printf("Looping over %u remaining bits in power:\n",start_index);
	}
#endif
	// LR binary powering loop:
	for(j = start_index-1; j >= 0; j--) {
		curr_bit = TEST_BIT128(pshift, j);
		SQR_LOHI128(x,lo,hi);	// x^2 mod q is returned in x
		MULL128(lo,qinv,lo);
		MULH128(q,lo,lo);
	#if FAC_DEBUG
		if(dbg) { pow = 2*pow + curr_bit - 128; printf("\tJ = %2u: [bit = %u]pow = %u, x = %s\n",j,curr_bit,pow,&char_buf[convert_uint128_base10_char(char_buf,x)]); }
	#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT128(hi, lo)) {
			SUB128(q, lo, lo);	ADD128(lo, hi, x);
		} else {
			SUB128(hi, lo, x);
		}
		if(curr_bit) {	// Combines overflow-on-add and need-to-subtract-q-from-sum checks:
			if(CMPUGT128(x, qhalf)){ ADD128(x, x, x); SUB128(x, q, x); }else{ ADD128(x, x, x); }
		}
	}
	// Since pre-subtracted lg2(R) = 128 from computed powermod exponent, no need to un-scale the loop output.
#if FAC_DEBUG
	if(dbg) printf("pow = %u, x = %s\n",pow,&char_buf[convert_uint128_base10_char(char_buf,x)]);
#endif
	// If we applied an initial right-justify shift to the modulus, restore the shift to the
	// current (partial) remainder and re-add the off-shifted part of the true remainder.
	if(nshift) {
		LSHIFT128(x,nshift,x);
	#if FAC_DEBUG
		if(dbg) printf("Restoring power-of-2: pow = %u, x *= 2^%u = %s\n",pow+nshift,nshift, &char_buf[convert_uint128_base10_char(char_buf, x)]);
	#endif
	}
#if FAC_DEBUG
	if(dbg) printf("xout = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
	return x;
}

/*
Function to find 2^(-p) mod q, where p is a 64-bit and q a 128-bit unsigned integer.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^128 (i.e. our MODQ
operation effects multiply modulo 2^128).

The key 3-operation sequence here is as follows:

	SQR_LOHI128(x,lo,hi);	// Input   x has 128 bits; lo/hi have 128 bits
	MULL128(lo,qinv,lo);	// Inputs lo & qinv, and output (overwrites lo) have 128 bits
	MULH128(q,lo,lo);	// Inputs  q &   lo, and output (overwrites lo) have 128 bits
*/
uint128 twopmodq128(uint128 p, uint128 q)
{
#if FAC_DEBUG
	int dbg = 0;//STREQ(&char_buf[convert_uint128_base10_char(char_buf, p)], "0");
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64;
	uint128 qhalf, qinv, x, lo, hi;
	static uint128 psave = {0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	uint32 FERMAT;
	if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

#if FAC_DEBUG
if(dbg) printf("twopmodq128:\n");
#endif

	RSHIFT_FAST128(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || !CMPEQ128(p, psave))
	{
		first_entry = FALSE;
		psave.d0 = p.d0;	psave.d1 = p.d1;
		pshift.d0 = p.d0 + 128;	pshift.d1 = p.d1 + (pshift.d0 < 128);
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
		if(pshift.d1) {
			j = leadz64(pshift.d1);
			start_index = 128-j-7;
			/* Extract leftmost 7 bits of pshift and subtract from 127: */
			zshift = 127 - (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 57);
		} else {
			j = leadz64(pshift.d0);
			start_index =  64-j-7;
			zshift = 127 - ((pshift.d0<<j) >> 57);
		}
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	#if FAC_DEBUG
		if(dbg) printf("leadb = %u\n",j);
		if(dbg) printf("zshift  = %u\n", zshift);
		if(dbg) printf("pshift = p + %u = %s\n",128,&char_buf[convert_uint128_base10_char(char_buf,pshift)]);
	#endif
		pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/*
	!    Find modular inverse (mod 2^128) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT((q.d0 & (uint64)1) == 1, "twopmodq128 : q must be odd for Montgomery-style modmul!");
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
	ASSERT(qinv.d1 == x.d1 && qinv.d0 == x.d0, "twopmodq128 : qinv.d1 == x.d1 && qinv.d0 == x.d0");
	if(dbg) printf("q    = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q   )]);
	if(dbg) printf("qinv = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv)]);
#endif

	/* Since zstart is a power of two < 2^128, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL128(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
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

	if(TEST_BIT128(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(CMPULT128(x, q), "twopmodq128 : CMPULT128(x,q)");
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
if(dbg) printf("j = %2d, x = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

		if(TEST_BIT128(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(CMPULT128(x, q), "twopmodq128 : CMPULT128(x,q)");
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
	//*** Aug 2022: see note at end of twopmodq96() for why we can add-sans-carry the constant (0 or 2) here:
	x.d0 += FERMAT;
	SUB128(x,q,x);

#if FAC_DEBUG
if(dbg) printf("x0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
	return x;
}

/*
!...Function to find 2^(-p) mod q, where p and q are both 128-bit unsigned integers.
*/
uint64 twopmodq128x2(uint64 *p_in, uint64 k)
{
	ASSERT(p_in != 0x0, "Null p_in pointer!");
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_mi64_base10_char(char_buf, p_in, 2, 0)], "0");
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64;
	uint128 p, q, qhalf, qinv, x, lo, hi;
	static uint128 psave = {0ull, 0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
if(dbg) printf("twopmodq128x2:\n");
#endif

	p.d0 = p_in[0]; p.d1 = p_in[1];
	uint32 FERMAT;
	if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	// Use x as tmp to hold 2*p:
	ADD128(p,p, x);
	ASSERT(!mi64_mul_scalar((uint64 *)&x, k, (uint64 *)&q, 2), "q must be < 2^128!");
	q.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */

	RSHIFT_FAST128(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || !CMPEQ128(p, psave))
	{
		first_entry = FALSE;
		psave = p;
		pshift.d0 = p_in[0] + 128;	pshift.d1 = p_in[1] + (pshift.d0 < 128);
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
			pshift.d1 = ~pshift.d1;
		}
		else
		{
			j = leadz64(pshift.d0);
			start_index =  64-j-7;
			zshift = 127 - ((pshift.d0<<j) >> 57);
		}
		pshift.d0 = ~pshift.d0;

		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	}

	/*
	!    Find modular inverse (mod 2^128) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT((q.d0 & (uint64)1) == 1, "twopmodq128x2 : q must be odd for Montgomery-style modmul!");
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
	ASSERT(qinv.d1 == x.d1 && qinv.d0 == x.d0, "twopmodq128x2 : qinv.d1 == x.d1 && qinv.d0 == x.d0");
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
		ASSERT(CMPULT128(x,q), "twopmodq128x2 : CMPULT128(x,q)");
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
			ASSERT(CMPULT128(x,q), "twopmodq128x2 : CMPULT128(x,q)");
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
	//*** Aug 2022: see note at end of twopmodq96() for why we can add-sans-carry the constant (0 or 2) here:
	x.d0 += FERMAT;
	SUB128(x,q,x);
#if FAC_DEBUG
	if(dbg) printf("Final x*= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
	return (uint64)CMPEQ128(x, ONE128) ;
}

// Second version of above, which takes factor candidate q in uint128 form:
uint64 twopmodq128x2B(uint64 *p_in, uint128 q)
{
	ASSERT(p_in != 0x0, "Null p_in pointer!");
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64;
	uint128 p, qhalf, qinv, x, lo, hi;
	static uint128 psave = {0ull, 0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	p.d0 = p_in[0]; p.d1 = p_in[1];
	uint32 FERMAT;
	if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	RSHIFT_FAST128(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || !CMPEQ128(p, psave)) {
		first_entry = FALSE;
		psave = p;
		pshift.d0 = p_in[0] + 128;	pshift.d1 = p_in[1] + (pshift.d0 < 128);
		if(pshift.d1) {
			j = leadz64(pshift.d1);
			start_index = 128-j-7;
			/* Extract leftmost 7 bits of pshift and subtract from 127: */
			zshift = 127 - (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 57);
			pshift.d1 = ~pshift.d1;
		} else {
			j = leadz64(pshift.d0);
			start_index =  64-j-7;
			zshift = 127 - ((pshift.d0<<j) >> 57);
		}
		pshift.d0 = ~pshift.d0;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	}

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT((q.d0 & (uint64)1) == 1, "twopmodq128x2B: q must be odd for Montgomery-style modmul!");
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv.d0 = (q.d0 + q.d0 + q.d0) ^ (uint64)2;	qinv.d1 = (uint64)0;
	for(j = 0; j < 4; j++) {
		lo64 = q.d0*qinv.d0;
		qinv.d0 = qinv.d0*((uint64)2 - lo64);
	}

	/* qinv has 128 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + __MULH64(q.d0, qinv.d0));
#else
	MULH64(q.d0, qinv.d0, lo64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + lo64);
#endif

	/* Since zstart is a power of two < 2^128, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL128(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT128(qinv, zshift, lo);
	MULH128(q,lo,lo);
	/* hi = 0 in this instance, which simplifies things. */
	SUB128(q, lo, x);

	if(TEST_BIT128(pshift, j)) {
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT128(x, qhalf)){ ADD128(x, x, x); SUB128(x, q, x); }else{ ADD128(x, x, x); }
	}

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI128(x,lo,hi);
		MULL128(lo,qinv,lo);
		MULH128(q,lo,lo);
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT128(hi, lo)) {
			SUB128(q, lo, lo);
			ADD128(lo, hi, x);
		} else {
			SUB128(hi, lo, x);
		}

		if(TEST_BIT128(pshift, j)) {
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT128(x, qhalf)){ ADD128(x, x, x); SUB128(x, q, x); }else{ ADD128(x, x, x); }
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD128(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^127, so x + x cannot overflow. */
	//*** Aug 2022: see note at end of twopmodq96() for why we can add-sans-carry the constant (0 or 2) here:
	x.d0 += FERMAT;
	SUB128(x,q,x);
#if FAC_DEBUG
	printf("Final x*= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
	return (uint64)CMPEQ128(x, ONE128) ;
}

/*** 4-trial-factor version ***/
uint64 twopmodq128_q4(uint64* p_in, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
	ASSERT(p_in != 0x0, "Null p_in pointer!");
	 int32 j;
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, lead7, r;
	uint128 p, q0, q1, q2, q3
		, qinv0, qinv1, qinv2, qinv3
		, qhalf0, qhalf1, qhalf2, qhalf3
		, x0, x1, x2, x3
		, lo0, lo1, lo2, lo3
		, hi0, hi1, hi2, hi3;
	static uint128 psave = {0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	p.d0 = p_in[0]; p.d1 = p_in[1];
	uint32 FERMAT;
	if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	// Use x0 as tmp to hold 2*p:
	ADD128(p,p, x0);
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k0, (uint64 *)&q0, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k1, (uint64 *)&q1, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k2, (uint64 *)&q2, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k3, (uint64 *)&q3, 2), "q must be < 2^128!");

	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	q2.d0 += 1;
	q3.d0 += 1;

	if(first_entry || !CMPEQ128(p, psave))
	{
		first_entry = FALSE;
		psave = p;
		pshift.d0 = p_in[0] + 128;	pshift.d1 = p_in[1] + (pshift.d0 < 128);
		/* Extract leftmost 7 bits of pshift and subtract from 128: */
		if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			lead7 = (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 57);
			start_index = 128-j-7;
			pshift.d1 = ~pshift.d1;
		}
		else
		{
			j = leadz64(pshift.d0);
			/* Extract leftmost 7 bits of pshift and subtract from 128: */
			lead7 = ((pshift.d0<<j) >> 57);
			start_index =  64-j-7;
		}
		pshift.d0 = ~pshift.d0;
  		zshift = 127 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	}
//=================
	RSHIFT_FAST128(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST128(q1, 1, qhalf1);
	RSHIFT_FAST128(q2, 1, qhalf2);
	RSHIFT_FAST128(q3, 1, qhalf3);

	// Find modular inverse (mod 2^128) of q in preparation for modular multiply:
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
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT128(x0, qhalf0)){ ADD128(x0, x0, x0); SUB128(x0, q0, x0); }else{ ADD128(x0, x0, x0); }
		if(CMPUGT128(x1, qhalf1)){ ADD128(x1, x1, x1); SUB128(x1, q1, x1); }else{ ADD128(x1, x1, x1); }
		if(CMPUGT128(x2, qhalf2)){ ADD128(x2, x2, x2); SUB128(x2, q2, x2); }else{ ADD128(x2, x2, x2); }
		if(CMPUGT128(x3, qhalf3)){ ADD128(x3, x3, x3); SUB128(x3, q3, x3); }else{ ADD128(x3, x3, x3); }
	}

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
	//*** Aug 2022: see note at end of twopmodq96() for why we can add-sans-carry the constant (0 or 2) here:
	x0.d0 += FERMAT;
	x1.d0 += FERMAT;
	x2.d0 += FERMAT;
	x3.d0 += FERMAT;
	SUB128(x0, q0, x0);
	SUB128(x1, q1, x1);
	SUB128(x2, q2, x2);
	SUB128(x3, q3, x3);

	/* Only do the full 128-bit (Xj== 1) check if the bottom 64 bits of Xj == 1: */
	r = 0;
	if(x0.d0 == 1) r += ((uint64)CMPEQ128(x0, ONE128) << 0);
	if(x1.d0 == 1) r += ((uint64)CMPEQ128(x1, ONE128) << 1);
	if(x2.d0 == 1) r += ((uint64)CMPEQ128(x2, ONE128) << 2);
	if(x3.d0 == 1) r += ((uint64)CMPEQ128(x3, ONE128) << 3);
	return(r);
}

/*** 8-trial-factor version ***/
uint64 twopmodq128_q8(uint64 *p_in, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
	ASSERT(p_in != 0x0, "Null p_in pointer!");
#if FAC_DEBUG
	int dbg = 0;
#endif
	 int32 j;
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, lo64_4, lo64_5, lo64_6, lo64_7, lead7, r;
	uint128 p, q0, q1, q2, q3, q4, q5, q6, q7
		,  qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, qhalf0, qhalf1, qhalf2, qhalf3, qhalf4, qhalf5, qhalf6, qhalf7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
		, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7;
	static uint128 psave = {0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
#if FAC_DEBUG
if(dbg) printf("twopmodq128_q8:\n");
#endif

	p.d0 = p_in[0]; p.d1 = p_in[1];
	uint32 FERMAT;
	if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	// Use x0 as tmp to hold 2*p:
	ADD128(p,p, x0);
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k0, (uint64 *)&q0, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k1, (uint64 *)&q1, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k2, (uint64 *)&q2, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k3, (uint64 *)&q3, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k4, (uint64 *)&q4, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k5, (uint64 *)&q5, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k6, (uint64 *)&q6, 2), "q must be < 2^128!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k7, (uint64 *)&q7, 2), "q must be < 2^128!");

	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	q2.d0 += 1;
	q3.d0 += 1;
	q4.d0 += 1;
	q5.d0 += 1;
	q6.d0 += 1;
	q7.d0 += 1;

	if(first_entry || !CMPEQ128(p, psave))
	{
		first_entry = FALSE;
		psave = p;
		pshift.d0 = p_in[0] + 128;	pshift.d1 = p_in[1] + (pshift.d0 < 128);
		/* Extract leftmost 7 bits of pshift and subtract from 128: */
		if(pshift.d1)
		{
			j = leadz64(pshift.d1);
			lead7 = (((pshift.d1<<j) + (pshift.d0>>(64-j))) >> 57);
			start_index = 128-j-7;
			pshift.d1 = ~pshift.d1;
		}
		else
		{
			j = leadz64(pshift.d0);
			/* Extract leftmost 7 bits of pshift and subtract from 128: */
			lead7 = ((pshift.d0<<j) >> 57);
			start_index =  64-j-7;
		}
		pshift.d0 = ~pshift.d0;
  		zshift = 127 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	}
//=================

	RSHIFT_FAST128(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST128(q1, 1, qhalf1);
	RSHIFT_FAST128(q2, 1, qhalf2);
	RSHIFT_FAST128(q3, 1, qhalf3);
	RSHIFT_FAST128(q4, 1, qhalf4);
	RSHIFT_FAST128(q5, 1, qhalf5);
	RSHIFT_FAST128(q6, 1, qhalf6);
	RSHIFT_FAST128(q7, 1, qhalf7);

	// Find modular inverse (mod 2^128) of q in preparation for modular multiply:
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

	for(j = start_index-2; j >= 0; j--)
	{
#if FAC_DEBUG
if(dbg) printf("A: x = %20" PRIu64 " + 2^64* %20" PRIu64 "\n",x0.d0,x0.d1);
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
if(dbg) printf("B: l = %20" PRIu64 " + 2^64* %20" PRIu64 "\n",lo0.d0,lo0.d1);
if(dbg) printf("B: h = %20" PRIu64 " + 2^64* %20" PRIu64 "\n",hi0.d0,hi0.d1);
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
if(dbg) printf("C: l = %20" PRIu64 " + 2^64* %20" PRIu64 "\n",lo0.d0,lo0.d1);
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
if(dbg) printf("D: l = %20" PRIu64 " + 2^64* %20" PRIu64 "\n",lo0.d0,lo0.d1);
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
if(dbg) printf("j = %2d, Res = %20" PRIu64 " + 2^64* %20" PRIu64,j,x0.d0,x0.d1);
#endif

		if(TEST_BIT128(pshift, j))
		{
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
if(dbg) printf(" *2 = %20" PRIu64 " + 2^64* %20" PRIu64,x0.d0,x0.d1);
#endif
		}
#if FAC_DEBUG
if(dbg) printf("\n");
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
if(dbg) printf("x0 = %20" PRIu64 " + 2^64* %20" PRIu64 "\n",x0.d0, x0.d1);
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

