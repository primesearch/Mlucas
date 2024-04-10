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
	#warning FAC_DEBUG = 1: Enabling dbg-printing.
	char char_buf[1024], str0[64], str1[64];
#endif

#define MONT_MUL192(__x,__y,__q,__qinv,__z)\
{\
	uint192 lo,hi;					\
	MUL_LOHI192(__x,__y,lo,hi);		\
	MULL192(__qinv,lo,lo);			\
	MULH192(__q,lo,lo);			\
	/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */\
	if(CMPULT192(hi, lo)) {	\
		SUB192(__q, lo, lo);	ADD192(lo, hi, __z);	\
	} else {	\
		SUB192(hi, lo, __z);	\
	}	\
};

/***********************************************************************************/
/***192-BIT INPUTS *****************************************************************/
/***********************************************************************************/
// Conventional positive-power version of twopmodq192, returns true mod:
uint192 twopmmodq192(uint192 p, uint192 q)
{
	 int32 j, pow;	// j needs to be signed because of the LR binary exponentiation
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint192_base10_char(char_buf, q)], "0");	// Replace "0" with "[desired decimal-form debug modulus]"
#endif
	uint32 curr_bit, leadb, start_index, nshift;
	uint64 lo64;
	uint192 pshift, qhalf, qinv, x, lo,hi, rsqr;
	// char_buf is local, cstr is globall available:
#if FAC_DEBUG
	if(dbg) printf("twopmmodq192: computing 2^%s (mod %s)\n",&char_buf[convert_uint192_base10_char(char_buf,p)],&cstr[convert_uint192_base10_char(cstr,q)]);
#endif
	RSHIFT_FAST192(q, 1, qhalf);	/* = (q-1)/2, since q odd. */
	// If p <= 192, directly compute 2^p (mod q):
	if(p.d2 == 0 && p.d1 == 0 && p.d0 <= 192) {
		// Lshift (1 << j) to align with leading bit of q, then do (p - j) repeated mod-doublings:
		x.d0 = 1ull; x.d1 = x.d2 = 0ull;
		j = leadz192(q);	// q >= 2^(192 - j - 1)
		j = (192 - j - 1);
		if(j > p.d0) {
			LSHIFT192(x,(uint32)p.d0,x);
		} else {
			LSHIFT192(x,j,x);
		}
		for( ; j < p.d0; j++) {
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT192(x, qhalf)){ ADD192(x, x, x); SUB192(x, q, x); }else{ ADD192(x, x, x); }
		}
		return x;
	}
	// If get here, p > 192: set up for Montgomery-mul-based powering loop:
	nshift = trailz192(q);
	if(nshift) {
		x.d0 = (uint64)nshift; x.d1 = x.d2 = 0ull; SUB192(p,x,p);	// p >= nshift guaranteed here:
		RSHIFT192(q,nshift,q);	// Right-shift dividend by (nshift) bits; for 2^p this means subtracting nshift from p
	#if FAC_DEBUG
		if(dbg) printf("Removed power-of-2 from q: q' = (q >> %u) = %s\n",nshift,&char_buf[convert_uint192_base10_char(char_buf,q)]);
	#endif
	}
	// Extract leftmost 8 bits of (p - 192); if > 192, use leftmost 7 instead:
	x.d0 = 192ull; x.d1 = x.d2 = 0ull; SUB192(p,x,pshift); j = leadz192(pshift);
	LSHIFT192(pshift,j,x);	leadb = x.d2 >> 56;	// leadb = (pshift<<j) >> 57; no (pshift = ~pshift) step in positive-power algorithm!
	if(leadb > 192) {
		start_index = 192-7-j;
		leadb >>= 1;
	} else {
		start_index = 192-8-j;
	}
#if FAC_DEBUG
	if(dbg) {
		printf("leadb = %u\n",j);
		printf("pshift = p - %u = %s\n",192,&char_buf[convert_uint192_base10_char(char_buf,pshift)]);
	}
#endif
	// Find inverse (mod 2^192) of q; q must be odd for Montgomery-style modmul to work:
	ASSERT((q.d0 & (uint64)1) == 1, "twopmmodq192 : q must be odd for Montgomery-style modmul to work");
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv.d0 = (q.d0 + q.d0 + q.d0) ^ (uint64)2;	qinv.d1 = qinv.d2 = 0ull;
	/* Compute qinv  = q^-1 (mod R = 2^192) via Newton iteration qinv = qinv*(2 - q*qinv), starting with
	5-bits-good approximation qinv_0 = 3*q ^ 2. Number-good-bits doubles each iteration until >= lg2(R) = 192:
	*/
	for(j = 0; j < 4; j++) 	{
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

	/* Initialize binary powering = R*x (mod q), where R = binary radix (2^192 here);
	at present don't care about optimizing this rarely-used function. */
	// First compute R^2 (mod q) in prep. for Mont-mul with initial seed:
	uint64 vtmp[7] = {0ull,0ull,0ull,0ull,0ull,0ull,1ull};	// R^2 = 2^384
	mi64_div_binary((const uint64*)vtmp, (const uint64*)&q, 7,3, 0x0, (uint32*)&j, (uint64*)&rsqr);

	// If leadb = 192, x = 2^192 = R, thus rsqr holds our desired starting value for x:
	if(leadb == 192)
		x = rsqr;
	else {
		x.d0 = 1ull; LSHIFT192(x,leadb,x);	// x <<= leadb;
		MONT_MUL192(x,rsqr, q,qinv, x);	// x*R (mod q) = MONT_MUL(x,R^2 (mod q),q,qinv)
 	}

#if FAC_DEBUG
	if(dbg) {
		printf("qinv = %s\n",&char_buf[convert_uint192_base10_char(char_buf,qinv)]);
		printf("leadb = %u, x0 = %s\n",leadb,&char_buf[convert_uint192_base10_char(char_buf,x)]);
		pow = leadb + 192;
		printf("Initial power = 2^(%u+192) = 2^%u mod q' = %s\n",leadb,pow,&char_buf[convert_uint192_base10_char(char_buf,x)]);
		printf("Looping over %u remaining bits in power:\n",start_index);
	}
#endif
	// LR binary powering loop:
	for(j = start_index-1; j >= 0; j--) {
		curr_bit = TEST_BIT192(pshift, j);
		SQR_LOHI192(x,lo,hi);	// x^2 mod q is returned in x
		MULL192(lo,qinv,lo);
		MULH192(q,lo,lo);
	#if FAC_DEBUG
		if(dbg) { pow = 2*pow + curr_bit - 192; printf("\tJ = %2u: [bit = %u]pow = %u, x = %s\n",j,curr_bit,pow,&char_buf[convert_uint192_base10_char(char_buf,x)]); }
	#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT192(hi, lo)) {
			SUB192(q, lo, lo);	ADD192(lo, hi, x);
		} else {
			SUB192(hi, lo, x);
		}
		if(curr_bit) {	// Combines overflow-on-add and need-to-subtract-q-from-sum checks:
			if(CMPUGT192(x, qhalf)){ ADD192(x, x, x); SUB192(x, q, x); }else{ ADD192(x, x, x); }
		}
	}
	// Since pre-subtracted lg2(R) = 192 from computed powermod exponent, no need to un-scale the loop output.
#if FAC_DEBUG
	if(dbg) printf("pow = %u, x = %s\n",pow,&char_buf[convert_uint192_base10_char(char_buf,x)]);
#endif
	// If we applied an initial right-justify shift to the modulus, restore the shift to the
	// current (partial) remainder and re-add the off-shifted part of the true remainder.
	if(nshift) {
		LSHIFT192(x,nshift,x);
	#if FAC_DEBUG
		if(dbg) printf("Restoring power-of-2: pow = %u, x *= 2^%u = %s\n",pow+nshift,nshift, &char_buf[convert_uint192_base10_char(char_buf, x)]);
	#endif
	}
#if FAC_DEBUG
	if(dbg) printf("xout = %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
#endif
	return x;
}

/*
Function to find 2^(-p) mod q, where p and q are both 192-bit unsigned integers.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^192 (i.e. our MODQ
operation effects multiply modulo 2^192).

The key 3-operation sequence here is as follows:

	SQR_LOHI192(x,lo,hi);	// Input   x has 192 bits; lo/hi have 192 bits
	MULL192(lo,qinv,lo);	// Inputs lo & qinv, and output (overwrites lo) have 192 bits
	MULH192(q,lo,lo);	// Inputs  q &   lo, and output (overwrites lo) have 192 bits
*/
uint192 twopmodq192(uint192 p, uint192 q)
{
#if FAC_DEBUG
	int dbg = 0;//STREQ(&char_buf[convert_uint192_base10_char(char_buf, q)], "569998349628599779154827250713823");
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead8, lo64;
	uint192 qhalf, qinv, x, lo, hi;
	static uint192 psave = {0ull,0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	uint32 FERMAT;
	if(p.d2 != 0ull)
		FERMAT = isPow2_64(p.d2) && (p.d1 == 0ull) && (p.d0 == 0ull);
	else if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	RSHIFT_FAST192(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	#if FAC_DEBUG
		// Now compute k = (q-1)/2p:
		mi64_div((uint64*)&qhalf,(uint64*)&p, 3,3, (uint64*)&x, (uint64*)&lo);	// x contains k; lo = (q-1)/2 % p
	//	dbg = (x.d0 == 488) && (x.d1 == 0 && x.d2 == 0);
	if(dbg) {
		ASSERT(mi64_iszero((uint64*)&lo, 3), "k must divide (q-1)/2!");
		printf("twopmodq192:\n");
	}
	#endif

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

		zshift = 191 - lead8;	/* zshift < 2^64  */
		/* Doubling the shift count here takes cares of the first SQR_LOHI */
		zshift <<= 1;			/* zshift < 2^128 */
	#if FAC_DEBUG
		if(dbg)	printf("lead8 = %u, zshift = %u\n", (uint32)lead8,zshift);
	#endif
		pshift.d2 = ~pshift.d2;	pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	/*
	!    Find modular inverse (mod 2^192) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT((q.d0 & (uint64)1) == 1, "twopmodq192 : q must be odd for Montgomery-style modmul to work!");
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
	if(dbg) {
		printf("q    = %s\n", &char_buf[convert_uint192_base10_char(char_buf, q   )]);
		printf("qinv = %s\n", &char_buf[convert_uint192_base10_char(char_buf, qinv)]);
	}
#endif
	/* Since zstart is a power of two < 2^192 (in fact < 2^128), use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL192(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("j = start_index - 1 = %u\n", j);
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
/*
Q: Can we also simplify the high-half mul here? By construction, q*qinv = u*2^192 + 1, where u < 2^192 and
is "discarded" in the computation of qinv, though it remains implicitly present every time we do this sequence
in the loop below:
		lo = MULL192(lo,qinv);
		lo = MULH192(q,lo).
If the input lo is a power of 2 which is < 2^192, can we avoid doing the MULs here and instead do some shifts
and masks, using the precomputed and stored value u := MULH192(q,qinv) ?
Let's have a look: Let lo = 2^x.

Then lo*qinv = qinv*lo = B*2^192 + A = (qinv << x) and taking the low 192 bits of this gives
MULL192(lo,qinv) = A = (qinv << x) & (2^192-1), i.e. do a 192-bit left off-shift of qinv by x bits. We discard the upper half B, although it's basically free.
In bitfield-form this looks like:

  qinv*lo = |000 (192-x bits) 000||-------------- qinv (192 bits) ---------------||00000000 (x bits) 000000|
          = |--------- MULH192(lo,qinv) (discarded) -------||-------------- MULL192(lo,qinv) --------------|

Now if we delay the modulo implied by the off-shift, the full 3-way product from which we will eventually extract the MULH192 output is
q*qinv*lo = ((u*2^192 + 1) << x). In bitfield-form this looks like:

q*qinv*lo = |000 (192-x bits) 000||-------------------------------------- q*qinv (384 bits) -------------------------------------||00000000 (x bits) 000000|
          = |000 (192-x bits) 000||-------------- u (192 bits) ------------------||0000000000000000 (192 bits) 000000000000000001||00000000 (x bits) 000000|

???????
          = |--------- MULH192(lo,qinv) (discarded) -------||-------------- MULL192(lo,qinv) --------------|
*/

	if(TEST_BIT192(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(CMPULT192(x,q), "twopmodq192 : CMPULT192(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT192(x, qhalf)){ ADD192(x, x, x); SUB192(x, q, x); }else{ ADD192(x, x, x); }
	#if FAC_DEBUG
		if(dbg) printf(", *2= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
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
		if(dbg) printf("j = %d, x = %s\n", j,&char_buf[convert_uint192_base10_char(char_buf, x)]);
	#endif

		if(TEST_BIT192(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(CMPULT192(x,q), "twopmodq192 : CMPULT192(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT192(x, qhalf)){ ADD192(x, x, x); SUB192(x, q, x); }else{ ADD192(x, x, x); }
		#if FAC_DEBUG
			if(dbg) printf(", *2= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x)]);
		#endif
		}
	}
	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD192(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^159, so x + x cannot overflow. */
	// For Fn with n > 50-or-so it is not uncommon to have q = b*2^64 + 1, thus need to check for borrow
	lo.d2 = lo.d1 = 0ull; lo.d0 = FERMAT;	SUB192(q,lo,q);	// Since carry may propagate > 1 word (e.g. F133 testfac), use general SUB
	SUB192(x,q,x);
	return x;
}

/*** 4-trial-factor version ***/
uint64 twopmodq192_q4(uint64 *p_in, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
#if FAC_DEBUG
	int dbg = 0;//(k0 == 460441ull);
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, lead8, r;
	uint192 p, q0, q1, q2, q3
			, qinv0, qinv1, qinv2, qinv3
			, qhalf0, qhalf1, qhalf2, qhalf3
			, x0, x1, x2, x3
			, lo0, lo1, lo2, lo3
			, hi0, hi1, hi2, hi3;
	uint192 x;
	static uint192 psave = {0ull,0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	p.d0 = p_in[0]; p.d1 = p_in[1]; p.d2 = p_in[2];
	uint32 FERMAT;
	if(p.d2 != 0ull)
		FERMAT = isPow2_64(p.d2) && (p.d1 == 0ull) && (p.d0 == 0ull);
	else if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	// Use x0 as tmp to hold 2*p:
	ADD192(p,p, x0);
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k0, (uint64 *)&q0, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k1, (uint64 *)&q1, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k2, (uint64 *)&q2, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k3, (uint64 *)&q3, 3), "q must be < 2^192!");
	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	q2.d0 += 1;
	q3.d0 += 1;
#if FAC_DEBUG
	if(dbg) {	// Compare vs 1-operand version:
		printf("twopmodq192_q4(p,q) with p = %s, q0 = %s\n",&str0[convert_uint192_base10_char(str0,p)],&char_buf[convert_uint192_base10_char(char_buf,q0)]);
		x0 = twopmodq192(p,q0);
		printf("Reference: twopmodq192(p,q) = %s\n", &char_buf[convert_uint192_base10_char(char_buf,x0)]);
	}
#endif
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
	#if FAC_DEBUG
		if(dbg) {
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

  		zshift = 191 - lead8;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	#if FAC_DEBUG
		if(dbg)	printf("lead8 = %u, zshift = %u\n", (uint32)lead8,zshift);
	#endif
		pshift.d2 = ~pshift.d2;	pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	// Find modular inverse (mod 2^128) of q in preparation for modular multiply:
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d2 = qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d2 = qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d2 = qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d2 = qinv3.d1 = (uint64)0;

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

	// Get next 64 bits of qinv via   qinv.d1 = MULL_64(-qinv.d0, MULL_64(q.d1, qinv.d0) + UMULH_64(q.d0, qinv.d0)):
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

	// Now that have bottom 128 bits of qinv, do one more Newton iteration using full 192-bit operands:
	MULL192(q0, qinv0, x0);	SUB192 (TWO192, x0, x0);	MULL192(qinv0, x0, qinv0);
	MULL192(q1, qinv1, x1);	SUB192 (TWO192, x1, x1);	MULL192(qinv1, x1, qinv1);
	MULL192(q2, qinv2, x2);	SUB192 (TWO192, x2, x2);	MULL192(qinv2, x2, qinv2);
	MULL192(q3, qinv3, x3);	SUB192 (TWO192, x3, x3);	MULL192(qinv3, x3, qinv3);

#if FAC_DEBUG
	if(dbg) {
		printf("q    = %s\n", &char_buf[convert_uint192_base10_char(char_buf, q0   )]);
		printf("qinv = %s\n", &char_buf[convert_uint192_base10_char(char_buf, qinv0)]);
	}
#endif
	/* Since zstart is a power of two < 2^192 (in fact < 2^128), use a streamlined code sequence for the first iteration: */
	j = start_index-1;

#if FAC_DEBUG
	if(dbg) printf("j = start_index - 1 = %u\n", j);
	LSHIFT192(qinv0, zshift, lo0);
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint192_base10_char(char_buf, lo0)]);
#endif
	/* MULL192(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv.
	hi = 0 in this instance, which simplifies things in the final subtract step.
	*/
	LSHIFT192(qinv0, zshift, lo0);	MULH192(q0,lo0,lo0);	SUB192(q0, lo0, x0);
	LSHIFT192(qinv1, zshift, lo1);	MULH192(q1,lo1,lo1);	SUB192(q1, lo1, x1);
	LSHIFT192(qinv2, zshift, lo2);	MULH192(q2,lo2,lo2);	SUB192(q2, lo2, x2);
	LSHIFT192(qinv3, zshift, lo3);	MULH192(q3,lo3,lo3);	SUB192(q3, lo3, x3);
#if FAC_DEBUG
	if(dbg) printf("q*lo/2^192 = %s\n", &char_buf[convert_uint192_base10_char(char_buf, lo0)]);
	if(dbg) printf("x = %s\n", &char_buf[convert_uint192_base10_char(char_buf, x0)]);
#endif
	if(TEST_BIT192(pshift, j))
	{
	#if FAC_DEBUG
		ASSERT(CMPULT192(x0,q0), "twopmodq192_q4: CMPULT192(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
		if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
		if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
		if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
	#if FAC_DEBUG
		if(dbg) printf(", *2= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x0)]);
	#endif
	}

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x0)]);
#endif
#if FAC_DEBUG
	if(CMPULT192(q0,x0)) {
		printf("(x0 = %s) >= (q = %s)\n", &str0[convert_uint192_base10_char(str0, x0)], &str1[convert_uint192_base10_char(str1, q0)] );
	}
#endif

	for(j = start_index-2; j >= 0; j--)
	{
	/*...x^2 mod q is returned in x. */
	#if PIPELINE_MUL192
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
	#if FAC_DEBUG
		if(dbg) printf("j = %d, x = %s\n", j,&char_buf[convert_uint192_base10_char(char_buf,x0)]);
	#endif

		if(TEST_BIT192(pshift, j))
		{
		#if FAC_DEBUG
			ASSERT(CMPULT192(x0,q0), "twopmodq192_q4 : CMPULT192(x,q)");
		#endif
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
			if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
			if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
			if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
		#if FAC_DEBUG
			if(dbg) printf(", *2= %s\n", &char_buf[convert_uint192_base10_char(char_buf, x0)]);
		#endif
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD192(x0 ,x0, x0);
	ADD192(x1 ,x1, x1);
	ADD192(x2 ,x2, x2);
	ADD192(x3 ,x3, x3);
	// For Fn with n > 50-or-so it is not uncommon to have q = b*2^64 + 1, thus need to check for borrow
	if(FERMAT) {
		lo0.d2 = lo0.d1 = 0ull; lo0.d0 = FERMAT;	// Since carry may propagate > 1 word (e.g. F133 testfac), use general SUB
		SUB192(q0,lo0,q0);
		SUB192(q1,lo0,q1);
		SUB192(q2,lo0,q2);
		SUB192(q3,lo0,q3);
	}
	SUB192(x0, q0, x0);
	SUB192(x1, q1, x1);
	SUB192(x2, q2, x2);
	SUB192(x3, q3, x3);

	/* Only do the full 192-bit (Xj== 1) check if the bottom 64 bits of Xj == 1: */
	r = 0;
	if(x0.d0 == 1) r += ((uint64)CMPEQ192(x0, ONE192) << 0);
	if(x1.d0 == 1) r += ((uint64)CMPEQ192(x1, ONE192) << 1);
	if(x2.d0 == 1) r += ((uint64)CMPEQ192(x2, ONE192) << 2);
	if(x3.d0 == 1) r += ((uint64)CMPEQ192(x3, ONE192) << 3);
	return(r);
}

/*** 4-trial-factor, using specialized O(n) version of UMULH,
mi64_mul_vector_hi_half for moduli q = 2.k.M(p) + 1, where M(p) is a Mersenne prime
***/
// MSVC does a horrible job with the uint192-specific macro, whereas GCC hates (in terms of producing horribly slow code) the cast-based version:
#ifdef COMPILER_TYPE_MSVC
	// Use mi64 version of hi-mul, with reinterpret-casts of uint192 <--> uint64[]:
	#define MULH192_QMMP(vin,p,k,vout,len) mi64_mul_vector_hi_qmmp((uint64 *)&vin, p, k, (uint64 *)&vout, 192);
//	#define MULH192_QMMP(vin,p,k,vout,len) mi64_mul_vector_hi_fast((uint64 *)&vin, p, k, (uint64 *)&vout, len);
#else
	// uint192-specific macro version ... assunmes vin != vout:
	#define MULH192_QMMP(__vin,__p,__k,__vout,__len)\
	{\
	uint32 __i, __nwshift, __rembits, __m64bits;\
	uint64 __lo, __hi, __bw, __cy, __cw, __cz, __nshift, *__zptr, __k2m1 = (__k << 1) - 1;\
	/* 1. compute z' = (2k-1).y via vector-scalar mul, the carryout word cw = ((2k-1).Y >> B); */\
		/* cw = mi64_mul_scalar(y,k2m1,z,len):	// z' = (2k-1).y */\
		MUL_LOHI64(__k2m1, __vin.d0, __vout.d0, __cy);\
		MUL_LOHI64(__k2m1, __vin.d1, __lo, __hi);\
		__vout.d1 = __lo + __cy;\
		__cy = __hi + (__vout.d1 < __lo);\
		MUL_LOHI64(__k2m1, __vin.d2, __lo, __hi);\
		__vout.d2 = __lo + __cy;\
		__cw = __hi + (__vout.d2 < __lo);	/* carryout into cw */\
		__lo = __vout.d2;	/* bw0 = z[len-1]; */\
/*if(__k==900) {printf("Macro: bw0 = %20llu, cw = %20llu, z` = %s\n", __lo,__cw,&char_buf[convert_uint192_base10_char(char_buf,__vout)]);}*/\
	/* 2. compute low n words of z = z' + y via vector-vector add, any carryout of that gets added to a 2nd copy of cw, cz: */\
		/* mi64_add(y,z,z, len):	// z = z' + y */\
		__vout.d0 = __vin.d0 + __vout.d0;\
		__cy = (__vout.d0 < __vin.d0);\
		/* Although &vin != &vout, still need 2-step add with 2 separate carry checks for d1 and d2-terms since any 2-term sum may = 2^64 exactly */\
		__hi = __vin.d1 + __cy;\
		__cy = (__hi < __vin.d1);\
		__vout.d1 += __hi;\
		__cy += (__vout.d1 < __hi);\
		/* d2 terms: */\
		__hi = __vin.d2 + __cy;\
		__cy = (__hi < __vin.d2);\
		__vout.d2 += __hi;\
		__cy += (__vout.d2 < __hi);\
/*if(__k==900) {printf("Macro: __vout.d2 [out] = %20llu\n", __vout.d2);}*/\
		__cz = __cw + __cy;	/* cz = cw + mi64_add(y,z,z, len);	// z = z' + y */\
/*if(__k==900) {printf("Macro: cz = %20llu, z = %s\n", __cz,&char_buf[convert_uint192_base10_char(char_buf,__vout)]);}*/\
\
	/* 3. compute low n words of z >> (b-p), then separately shift in cz from the left, via (2^b*cz) >> (b-p) = (cz << p). */\
		/* bw1 = mi64_shrl(z,z,nshift,len);	// low n words of z >> (b-p); high 64 bits of off-shifted portion saved in bw1 */\
		__nshift = (__len<<6) - __p;\
		__nwshift = __len - (__p >> 6); __rembits = (__nshift & 63);\
		__hi = 0ull;	/* hi stands in for bw1 */\
		/* hi plays the part of "y[-1]": Must do this *before* the loop because of in-place possibility: */\
		if(__nwshift)\
		{\
			__hi = __vout.d0;\
			/* Take care of the whole-word part of the shift: */\
			__vout.d0 = __vout.d1;\
			__vout.d1 = __vout.d2;\
			__vout.d2 = 0ull;\
		}\
		/* If __nshift not an exact multiple of the wordlength, take care of remaining shift bits: */\
		if(__rembits)\
		{\
			__m64bits = (64-__rembits);\
			__hi = (__hi >> __rembits) + (__vout.d0 << __m64bits);\
			/* Process all but the most-significant element, in reverse order: */\
			__vout.d0 = (__vout.d0 >> __rembits) + (__vout.d1 << __m64bits);\
			__vout.d1 = (__vout.d1 >> __rembits) + (__vout.d2 << __m64bits);\
			/* Most-significant element gets zeros shifted in from the left: */\
			__vout.d2 >>= __rembits;\
		}\
/*if(__k==900) {printf("Macro: bw1 = %20llu, z>> = %s\n", __hi,&char_buf[convert_uint192_base10_char(char_buf,__vout)]);}*/\
\
		/* Check for borrow-on-subtract of to-be-off-shifted sections: */\
		__bw = (__lo > __hi);\
		/* Add (cw << (b-p)) to result: */\
		__rembits = (__p&63);\
		__nshift = (__p >> 6);	/* i = index of low word into which cw will get added */\
		__zptr = (uint64 *)&__vout + __nshift;	/* z+i */\
		/* If (b-p) == 0 (mod 64) all of cz goes into z[nshift], with nshift = (b-p)/64: */\
		if(__rembits == 0) {\
			/* mi64_add_scalar(__zptr, __cz, __zptr, __rembits);	// mi64_add_scalar(&z[nshift],cz,&z[nshift],len-i): */\
			/* In-place: Only need to proceed until carry peters out: */\
			__cy = 0ull;\
			__rembits = 3-__nshift;	/* use this to store loop upper-bound */\
			for(__i = 0; __i < __rembits; __i++) {\
				*__zptr += __cy;\
				__cy = (*__zptr < __cy);\
				if(!__cy) { break; }\
			}\
		/* Otherwise cz gets split between z[nshift] and z[nshift+1]: */\
		} else {\
			/* low 64-(p%64) rembits of cz = (cz << rembits) go into z[i]: */\
			__cy = (__cz << __rembits);\
			*__zptr += __cy; __cy = (*__zptr++ < __cy);\
			/* high (p%64) bits of cw = (cw >> bits) go into z[i+1]: */\
			/* mi64_add_scalar(__zptr, (__cz >> (64-__rembits)) + __cy, __zptr, 2-__nshift); */\
			__m64bits = (64-__rembits);\
			__rembits = 2-__nshift;	/* use this to store loop upper-bound */\
			for(__i = 0; __i < __rembits; __i++) {\
				*__zptr += __cy;\
				__cy = (*__zptr < __cy);\
				if(!__cy) { break; }\
			}\
		}\
\
	/* 4. subtract scalar (bw + cw) from resulting vector to effect ... - (2k-1).Y step in [*]: */\
		/* __zptr = (uint64 *)&__vout;	mi64_sub_scalar(__zptr, (__bw + __cw), __zptr, 3); */\
		__bw += __cw;\
		while(__bw)\
		{\
			__cy = __vout.d0 - __bw;\
			__bw = (__cy > __vout.d0);\
			__vout.d0 = __cy;\
			if(!__bw) break;\
			__cy = __vout.d1 - __bw;\
			__bw = (__cy > __vout.d1);\
			__vout.d1 = __cy;\
			if(!__bw) break;\
			__cy = __vout.d2 - __bw;\
			__bw = (__cy > __vout.d2);\
			__vout.d2 = __cy;\
			ASSERT(!__bw, "bw != 0");\
		}\
	}

#endif

// Specialized version of above _q4 modpow routine, for double-Mersenne numbers M(M(p)):
uint64 twopmodq192_q4_qmmp(uint64 *p_in, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
#if FAC_DEBUG
	int dbg = 0;//(k0==900 || k1==900 || k2==900 || k3==900);
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, lead8, r;
	uint192 p, q0, q1, q2, q3
			, qinv0, qinv1, qinv2, qinv3
			, qhalf0, qhalf1, qhalf2, qhalf3
			, x0, x1, x2, x3
			, lo0, lo1, lo2, lo3
			, hi0, hi1, hi2, hi3;
	uint192 x;
	static uint64 mmpsave;
	static uint192 psave = {0ull,0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#if FAC_DEBUG
//	if(dbg) { printf("k0-3 = %u %u %u %u\n", k0,k1,k2,k3); exit(0); }
#endif
	p.d0 = p_in[0]; p.d1 = p_in[1]; p.d2 = p_in[2];
	// Use x0 as tmp to hold 2*p:
	ADD192(p,p, x0);
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k0, (uint64 *)&q0, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k1, (uint64 *)&q1, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k2, (uint64 *)&q2, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k3, (uint64 *)&q3, 3), "q must be < 2^192!");

	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	q2.d0 += 1;
	q3.d0 += 1;

	RSHIFT_FAST192(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST192(q1, 1, qhalf1);
	RSHIFT_FAST192(q2, 1, qhalf2);
	RSHIFT_FAST192(q3, 1, qhalf3);

	if(first_entry || !CMPEQ192(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		mmpsave = 192 - mi64_leadz(p_in,3);
		// Check that it's really a double-Mersenne: Adding one, right-shift by mmpsave = #bits give 1:
		mi64_add_scalar(p_in, 1ull, (uint64*)&x, 3);
		mi64_shrl((uint64*)&x, (uint64*)&x, mmpsave, 3,3);
		--x.d0;	ASSERT(mi64_iszero((uint64*)&x, 3), "MMp check failed!");
		x.d0 = 192; x.d1 = x.d2 = 0;
		ADD192(p, x, pshift);
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

	// Find modular inverse (mod 2^128) of q in preparation for modular multiply:
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d2 = qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d2 = qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d2 = qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d2 = qinv3.d1 = (uint64)0;

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

	// Get next 64 bits of qinv via   qinv.d1 = MULL_64(-qinv.d0, MULL_64(q.d1, qinv.d0) + UMULH_64(q.d0, qinv.d0)):
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

	// Now that have bottom 128 bits of qinv, do one more Newton iteration using full 192-bit operands:
	MULL192(q0, qinv0, x0);	SUB192 (TWO192, x0, x0);	MULL192(qinv0, x0, qinv0);
	MULL192(q1, qinv1, x1);	SUB192 (TWO192, x1, x1);	MULL192(qinv1, x1, qinv1);
	MULL192(q2, qinv2, x2);	SUB192 (TWO192, x2, x2);	MULL192(qinv2, x2, qinv2);
	MULL192(q3, qinv3, x3);	SUB192 (TWO192, x3, x3);	MULL192(qinv3, x3, qinv3);

	/* Since zstart is a power of two < 2^192, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL192(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv.
	Put MULL output into x, because MULH192_QMMP requires a separate in/out arrays.
	hi = 0 in this instance, which simplifies things in the final subtract step.
	*/
	LSHIFT192(qinv0, zshift, x0);	MULH192_QMMP(x0,mmpsave,k0,lo0,3);	SUB192(q0, lo0, x0);
	LSHIFT192(qinv1, zshift, x1);	MULH192_QMMP(x1,mmpsave,k1,lo1,3);	SUB192(q1, lo1, x1);
	LSHIFT192(qinv2, zshift, x2);	MULH192_QMMP(x2,mmpsave,k2,lo2,3);	SUB192(q2, lo2, x2);
	LSHIFT192(qinv3, zshift, x3);	MULH192_QMMP(x3,mmpsave,k3,lo3,3);	SUB192(q3, lo3, x3);
#if FAC_DEBUG
if(dbg) {
	printf("pshift = %s\n", &char_buf[convert_uint192_base10_char(char_buf, pshift)]);
	printf("qinv   = %s\n", &char_buf[convert_uint192_base10_char(char_buf, qinv3)]);
	printf("umulh  = %s\n", &char_buf[convert_uint192_base10_char(char_buf, lo3)]);
	printf("x`     = %s\n", &char_buf[convert_uint192_base10_char(char_buf, x3)]);
//	exit(0);
}
#endif

	if(TEST_BIT192(pshift, j))
	{
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
		if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
		if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
		if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
	}

	for(j = start_index-2; j >= 0; j--)
	{
	/*...x^2 mod q is returned in x. */
	#if PIPELINE_MUL192
		SQR_LOHI192_q4(
		  x0, lo0, hi0
		, x1, lo1, hi1
		, x2, lo2, hi2
		, x3, lo3, hi3);
		// Put MULL output into x, because MULH192_QMMP requires a separate in/out arrays:
		MULL192_q4(
		  lo0, qinv0, x0
		, lo1, qinv1, x1
		, lo2, qinv2, x2
		, lo3, qinv3, x3);
	#else
	/* Strangely, on alpha/ia64 nopipe is better than the pipelined ..._q4 version: */
		SQR_LOHI192(x0,lo0,hi0);
		SQR_LOHI192(x1,lo1,hi1);
		SQR_LOHI192(x2,lo2,hi2);
		SQR_LOHI192(x3,lo3,hi3);
	#if FAC_DEBUG
		if(dbg) printf("j = %d, lo192 = %s\n", j,&char_buf[convert_uint192_base10_char(char_buf,lo3)]);
		if(dbg) printf("j = %d, hi192 = %s\n", j,&char_buf[convert_uint192_base10_char(char_buf,hi3)]);
	#endif

		// Put MULL output into x, because MULH192_QMMP requires a separate in/out arrays:
		MULL192(lo0,qinv0,x0);
		MULL192(lo1,qinv1,x1);
		MULL192(lo2,qinv2,x2);
		MULL192(lo3,qinv3,x3);
	#if FAC_DEBUG
		if(dbg){printf("j = %d, mull  = %s\n", j,&char_buf[convert_uint192_base10_char(char_buf, x3)]);}
	#endif
	#endif

	#if 1
		MULH192_QMMP(x0,mmpsave,k0,lo0,3);	//MULH192(x0,q0,x);	ASSERT(CMPEQ192(lo0, x), "MULH192_QMMP fail!");
		MULH192_QMMP(x1,mmpsave,k1,lo1,3);	//MULH192(x1,q1,x);	ASSERT(CMPEQ192(lo1, x), "MULH192_QMMP fail!");
		MULH192_QMMP(x2,mmpsave,k2,lo2,3);	//MULH192(x2,q2,x);	ASSERT(CMPEQ192(lo2, x), "MULH192_QMMP fail!");
		MULH192_QMMP(x3,mmpsave,k3,lo3,3);	//MULH192(x3,q3,x);	ASSERT(CMPEQ192(lo3, x), "MULH192_QMMP fail!");
	#else
		MULH192(x0,q0,lo0);
		MULH192(x1,q1,lo1);
		MULH192(x2,q2,lo2);
		MULH192(x3,q3,lo3);
	#endif
	#if FAC_DEBUG
		if(dbg) printf("j = %d, x = %s\n", j,&char_buf[convert_uint192_base10_char(char_buf,lo3)]);
	#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT192(hi0, lo0)) { SUB192(q0, lo0, lo0);	ADD192(lo0, hi0, x0); } else { SUB192(hi0, lo0, x0); }
		if(CMPULT192(hi1, lo1)) { SUB192(q1, lo1, lo1);	ADD192(lo1, hi1, x1); } else { SUB192(hi1, lo1, x1); }
		if(CMPULT192(hi2, lo2)) { SUB192(q2, lo2, lo2);	ADD192(lo2, hi2, x2); } else { SUB192(hi2, lo2, x2); }
		if(CMPULT192(hi3, lo3)) { SUB192(q3, lo3, lo3);	ADD192(lo3, hi3, x3); } else { SUB192(hi3, lo3, x3); }

		if(TEST_BIT192(pshift, j))
		{
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
			if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
			if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
			if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
		}
	}

	//...Double and return:
	if(CMPUGT192(x0, qhalf0)){ ADD192(x0, x0, x0); SUB192(x0, q0, x0); }else{ ADD192(x0, x0, x0); }
	if(CMPUGT192(x1, qhalf1)){ ADD192(x1, x1, x1); SUB192(x1, q1, x1); }else{ ADD192(x1, x1, x1); }
	if(CMPUGT192(x2, qhalf2)){ ADD192(x2, x2, x2); SUB192(x2, q2, x2); }else{ ADD192(x2, x2, x2); }
	if(CMPUGT192(x3, qhalf3)){ ADD192(x3, x3, x3); SUB192(x3, q3, x3); }else{ ADD192(x3, x3, x3); }
#if FAC_DEBUG
	if(dbg) {
		printf("Final x = %s\n",&char_buf[convert_uint192_base10_char(char_buf, x3)]);
		exit(0);
	}
#endif

	/* Only do the full 192-bit (Xj== 1) check if the bottom 64 bits of Xj == 1: */
	r = 0;
	if(x0.d0 == 1) r += ((uint64)CMPEQ192(x0, ONE192) << 0);
	if(x1.d0 == 1) r += ((uint64)CMPEQ192(x1, ONE192) << 1);
	if(x2.d0 == 1) r += ((uint64)CMPEQ192(x2, ONE192) << 2);
	if(x3.d0 == 1) r += ((uint64)CMPEQ192(x3, ONE192) << 3);
	return(r);
}

/*** 8-trial-factor version ***/
uint64 twopmodq192_q8(uint64 *p_in, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lo64_0, lo64_1, lo64_2, lo64_3, lo64_4, lo64_5, lo64_6, lo64_7, lead8, r;
	uint192 p, q0, q1, q2, q3, q4, q5, q6, q7
			, qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
			, qhalf0, qhalf1, qhalf2, qhalf3, qhalf4, qhalf5, qhalf6, qhalf7
			, x0, x1, x2, x3, x4, x5, x6, x7
			, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
			, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7;
	uint192 x;
	static uint192 psave = {0ull,0ull,0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	p.d0 = p_in[0]; p.d1 = p_in[1]; p.d2 = p_in[2];
	uint32 FERMAT;
	if(p.d2 != 0ull)
		FERMAT = isPow2_64(p.d2) && (p.d1 == 0ull) && (p.d0 == 0ull);
	else if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	// Use x0 as tmp to hold 2*p:
	ADD192(p,p, x0);
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k0, (uint64 *)&q0, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k1, (uint64 *)&q1, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k2, (uint64 *)&q2, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k3, (uint64 *)&q3, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k4, (uint64 *)&q4, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k5, (uint64 *)&q5, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k6, (uint64 *)&q6, 3), "q must be < 2^192!");
	ASSERT(!mi64_mul_scalar((uint64 *)&x0, k7, (uint64 *)&q7, 3), "q must be < 2^192!");

	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	q2.d0 += 1;
	q3.d0 += 1;
	q4.d0 += 1;
	q5.d0 += 1;
	q6.d0 += 1;
	q7.d0 += 1;

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

	// Find modular inverse (mod 2^128) of q in preparation for modular multiply:
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d2 = qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d2 = qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d2 = qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d2 = qinv3.d1 = (uint64)0;
	qinv4.d0 = (q4.d0 + q4.d0 + q4.d0) ^ (uint64)2;	qinv4.d2 = qinv4.d1 = (uint64)0;
	qinv5.d0 = (q5.d0 + q5.d0 + q5.d0) ^ (uint64)2;	qinv5.d2 = qinv5.d1 = (uint64)0;
	qinv6.d0 = (q6.d0 + q6.d0 + q6.d0) ^ (uint64)2;	qinv6.d2 = qinv6.d1 = (uint64)0;
	qinv7.d0 = (q7.d0 + q7.d0 + q7.d0) ^ (uint64)2;	qinv7.d2 = qinv7.d1 = (uint64)0;

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

	for(j = start_index-2; j >= 0; j--)
	{
	/*...x^2 mod q is returned in x. */
	#if PIPELINE_MUL192
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
	ADD192(x0 ,x0, x0);
	ADD192(x1 ,x1, x1);
	ADD192(x2 ,x2, x2);
	ADD192(x3 ,x3, x3);
	ADD192(x4 ,x4, x4);
	ADD192(x5 ,x5, x5);
	ADD192(x6 ,x6, x6);
	ADD192(x7 ,x7, x7);
	if(FERMAT) {
		lo0.d2 = lo0.d1 = 0ull; lo0.d0 = FERMAT;	// Since carry may propagate > 1 word (e.g. F133 testfac), use general SUB
		SUB192(q0,lo0,q0);
		SUB192(q1,lo0,q1);
		SUB192(q2,lo0,q2);
		SUB192(q3,lo0,q3);
		SUB192(q4,lo0,q4);
		SUB192(q5,lo0,q5);
		SUB192(q6,lo0,q6);
		SUB192(q7,lo0,q7);
	}
	SUB192(x0, q0, x0);
	SUB192(x1, q1, x1);
	SUB192(x2, q2, x2);
	SUB192(x3, q3, x3);
	SUB192(x4, q4, x4);
	SUB192(x5, q5, x5);
	SUB192(x6, q6, x6);
	SUB192(x7, q7, x7);

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

/*
Good:
pshift = 6277101735386680763835789423049210091073826769276946612032
qinv   = 4249870423643896013758727800263806032353769591780234971465
x`     = 101785508172855884961760149503358
j = 88, lo192 = 2307869884805901223499303143938014440754487531554898070020
j = 88, hi192 = 1650489
j = 88, mull  = 614654487165510378251812999766475911146375692465063739172
j = 88, x = 27928909160267978861999257736012
j = 87, x = 83737795926253879717701166788654
j = 86, x = 38657137461790946051076066555840
j = 85, x = 39309390462328548112272409919459
j = 84, x = 225543367269777476900523971906596
j = 83, x = 16937395093467294434228161485550
j = 82, x = 52054413785730228403548779747850
j = 81, x = 100196111389013757184049357920276
j = 80, x = 37765253136558824295972530487457
j = 79, x = 59565409497751052946259494818941
j = 78, x = 115270234154731440654023261356344
j = 77, x = 54180909019591871072874929271682
j = 76, x = 39623315280572538344995178094014
j = 75, x = 127145906491756825176102542178185
j = 74, x = 84424154639354248486199519053285
j = 73, x = 136791523179468741731124936440281
j = 72, x = 42326871110060381298834509744979
j = 71, x = 238651999924030502645753166275826
j = 70, x = 181550406295368648469118554709825
j = 69, x = 156514469790176095389204206948931
j = 68, x = 22541504766401649939042285417434
j = 67, x = 76903282459274333272018893322306
j = 66, x = 13664852293193957489684078746302
j = 65, x = 54532416259590898580662714671125
j = 64, x = 191798668298230317338054211697250
j = 63, x = 216938377281815535076503447528405
j = 62, x = 107985826597111309963874785182867
j = 61, x = 145555090635103908193440750689944
j = 60, x = 159787122695578918135522978597759
j = 59, x = 243117578510900771695075460784442
j = 58, x = 138847375315209812118571176471504
j = 57, x = 201593799535256183463986370966423
j = 56, x = 140772200362712948753593297959442
j = 55, x = 215999525398574923206235433818409
j = 54, x = 214180714999124008695199390982083
j = 53, x = 89210777703487099348094768947053
j = 52, x = 190560440449678810075430215482433
j = 51, x = 150171250381073615316190846148993
j = 50, x = 177110584909956222590648279575193
j = 49, x = 164830198370292248371394492405827
j = 48, x = 245674334884056672439157626974142
j = 47, x = 118821749248607011907949354016590
j = 46, x = 141755767500945933367594210626243
j = 45, x = 3635565132091631480272315531524
j = 44, x = 135978882207900502637476231363668
j = 43, x = 93456468784858326616477107609949
j = 42, x = 85154905018628079927611723975151
j = 41, x = 184694205220174018346251227841881
j = 40, x = 23449771345365861335887873276200
j = 39, x = 142048668539424360842648993229672
j = 38, x = 242549386252727665781470673296842
j = 37, x = 122754736292091812498895209817825
j = 36, x = 139263608183566868784018350525057
j = 35, x = 214041366752509840863768815987335
j = 34, x = 144029448046935736767088941802387
j = 33, x = 55311893908534004532676554475296
j = 32, x = 103919833466840951347514126478772
j = 31, x = 254262641105532406982548512839416
j = 30, x = 243806173561487226374780153801324
j = 29, x = 43530383534983198248129726665339
j = 28, x = 223499957507363569113558469300943
j = 27, x = 76144763137639228360244902381965
j = 26, x = 177578767097520325602590693481500
j = 25, x = 32961836846102346154119094783215
j = 24, x = 198612709730785426222704521788610
j = 23, x = 46929462518342229334291051145023
j = 22, x = 42576567486656682210463266984071
j = 21, x = 46649247508444037663665912072189
j = 20, x = 218545104774736963855780738022385
j = 19, x = 125449400505904927102946798858318
j = 18, x = 60987305371745659000513543423487
j = 17, x = 249751889511166889686357125400032
j = 16, x = 78544050505203299644996636757348
j = 15, x = 113776067172341897183975731737641
j = 14, x = 264570876537333033680917613042333
j = 13, x = 32326339293057949532140044876796
j = 12, x = 257197893627768616747977113332817
j = 11, x = 189816187977296012484719239912551
j = 10, x = 34957015492930785080143467015123
j = 9, x = 138426901347332711626049083801894
j = 8, x = 82879317131706515433715666840394
j = 7, x = 99710159123210341697421086358890
j = 6, x = 124502691823884579280566892058966
j = 5, x = 202710371155775592638518110720051
j = 4, x = 238168071194091464328775678969401
j = 3, x = 80511908756680588114691729279299
j = 2, x = 1007885424503737279669410203767
j = 1, x = 112495345681500264232012357158068
j = 0, x = 142610692525675807668379115356775
Final x = 1

Macro:
*/

