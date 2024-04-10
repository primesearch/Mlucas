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

#warning To-Do: Add support for PIPELINE_MUL256

#include "factor.h"
#include "imul256_macro.h"

#ifndef FAC_DEBUG
	#define FAC_DEBUG	0
#elif FAC_DEBUG < 0 || FAC_DEBUG > 2
	#error If FAC_DEBUG defined it must be assigned a value of 0 (off) or 1 (on).
#endif
#if FAC_DEBUG
	#warning FAC_DEBUG = 1: Enabling dbg-printing.
	char char_buf[1024], str0[64], str1[64];
#endif

// Specialized versions of MULL128, in which we multiply the
// low half of the 256-bit input x with the high half of the 256-bit input y:
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULL128_256lohihalf(__x, __y, __lo)\
    {\
		uint64 __w0,__w1,__a,__c;\
		\
		__a = __MULL64(__x.d0,__y.d3);\
		__c = __MULL64(__y.d2,__x.d1);\
		MUL_LOHI64(__x.d0,__y.d2,&__w0,&__w1);\
		__w1 += __a;\
		__w1 += __c;\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
    }

#else

    #define MULL128_256lohihalf(__x, __y, __lo)\
    {\
		uint64 __w0,__w1,__a,__c;\
		\
		MULL64(__x.d0,__y.d3,__a);\
		MULL64(__y.d2,__x.d1,__c);\
		MUL_LOHI64(__x.d0,__y.d2, __w0, __w1);\
		__w1 += __a;\
		__w1 += __c;\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
    }

#endif

#define MONT_MUL256(__x,__y,__q,__qinv,__z)\
{\
	uint256 lo,hi;					\
	MUL_LOHI256(__x,__y,lo,hi);		\
	MULL256(__qinv,lo,lo);			\
	MULH256(__q,lo,lo);			\
	/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */\
	if(CMPULT256(hi, lo)) {	\
		SUB256(__q, lo, lo);	ADD256(lo, hi, __z);	\
	} else {	\
		SUB256(hi, lo, __z);	\
	}	\
};

/***********************************************************************************/
/***256-BIT INPUTS *****************************************************************/
/***********************************************************************************/
// Conventional positive-power version of twopmodq256, returns true mod:
uint256 twopmmodq256(uint256 p, uint256 q)
{
	 int32 j, pow;	// j needs to be signed because of the LR binary exponentiation
#if FAC_DEBUG
	int dbg = STREQ(&char_buf[convert_uint256_base10_char(char_buf, q)], "0");	// Replace "0" with "[desired decimal-form debug modulus]"
#endif
	uint32 curr_bit, leadb, start_index, nshift;
	uint64 lo64;
	uint256 pshift, qhalf, qinv, x, lo,hi, rsqr;
	// char_buf is local, cstr is globall available:
#if FAC_DEBUG
	if(dbg) printf("twopmmodq256: computing 2^%s (mod %s)\n",&char_buf[convert_uint256_base10_char(char_buf,p)],&cstr[convert_uint256_base10_char(cstr,q)]);
#endif
	RSHIFT_FAST256(q, 1, qhalf);	/* = (q-1)/2, since q odd. */
	// If p <= 256, directly compute 2^p (mod q):
	if(p.d3 == 0 && p.d2 == 0 && p.d1 == 0 && p.d0 <= 256) {
		// Lshift (1 << j) to align with leading bit of q, then do (p - j) repeated mod-doublings:
		x.d0 = 1ull; x.d1 = x.d2 = x.d3 = 0ull;
		j = leadz256(q);	// q >= 2^(256 - j - 1)
		j = (256 - j - 1);
		if(j > p.d0) {
			LSHIFT256(x,(uint32)p.d0,x);
		} else {
			LSHIFT256(x,j,x);
		}
		for( ; j < p.d0; j++) {
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT256(x, qhalf)){ ADD256(x, x, x); SUB256(x, q, x); }else{ ADD256(x, x, x); }
		}
		return x;
	}
	// If get here, p > 256: set up for Montgomery-mul-based powering loop:
	nshift = trailz256(q);
	if(nshift) {
		x.d0 = (uint64)nshift; x.d1 = x.d2 = x.d3 = 0ull; SUB256(p,x,p);	// p >= nshift guaranteed here:
		RSHIFT256(q,nshift,q);	// Right-shift dividend by (nshift) bits; for 2^p this means subtracting nshift from p
	#if FAC_DEBUG
		if(dbg) printf("Removed power-of-2 from q: q' = (q >> %u) = %s\n",nshift,&char_buf[convert_uint256_base10_char(char_buf,q)]);
	#endif
	}
	// Extract leftmost 8 bits of (p - 256):
	x.d0 = 256ull; x.d1 = x.d2 = x.d3 = 0ull; SUB256(p,x,pshift); j = leadz256(pshift);
	LSHIFT256(pshift,j,x);	leadb = x.d3 >> 56;	// leadb = (pshift<<j) >> 57; no (pshift = ~pshift) step in positive-power algorithm!
	start_index = 256-8-j;
#if FAC_DEBUG
	if(dbg) {
		printf("leadb = %u\n",j);
		printf("pshift = p - %u = %s\n",256,&char_buf[convert_uint256_base10_char(char_buf,pshift)]);
	}
#endif
	// Find inverse (mod 2^256) of q; q must be odd for Montgomery-style modmul to work:
	ASSERT((q.d0 & (uint64)1) == 1, "twopmmodq256 : q must be odd for Montgomery-style modmul to work");
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv.d0 = (q.d0 + q.d0 + q.d0) ^ (uint64)2;	qinv.d1 = qinv.d2 = qinv.d3 = 0ull;
	/* Compute qinv  = q^-1 (mod R = 2^256) via Newton iteration qinv = qinv*(2 - q*qinv), starting with
	5-bits-good approximation qinv_0 = 3*q ^ 2. Number-good-bits doubles each iteration until >= lg2(R) = 256:
	*/
	for(j = 0; j < 4; j++) 	{
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

	/* Initialize binary powering = R*x (mod q), where R = binary radix (2^256 here);
	at present don't care about optimizing this rarely-used function. */
	// First compute R^2 (mod q) in prep. for Mont-mul with initial seed:
	uint64 vtmp[9] = {0ull,0ull,0ull,0ull,0ull,0ull,0ull,0ull,1ull};	// R^2 = 2^384
	mi64_div_binary((const uint64*)vtmp, (const uint64*)&q, 9,4, 0x0, (uint32*)&j, (uint64*)&rsqr);

	// If leadb = 256, x = 2^256 = R, thus rsqr holds our desired starting value for x:
	if(leadb == 256)
		x = rsqr;
	else {
		x.d0 = 1ull; LSHIFT256(x,leadb,x);	// x <<= leadb;
		MONT_MUL256(x,rsqr, q,qinv, x);	// x*R (mod q) = MONT_MUL(x,R^2 (mod q),q,qinv)
 	}

#if FAC_DEBUG
	if(dbg) {
		printf("qinv = %s\n",&char_buf[convert_uint256_base10_char(char_buf,qinv)]);
		printf("leadb = %u, x0 = %s\n",leadb,&char_buf[convert_uint256_base10_char(char_buf,x)]);
		pow = leadb + 256;
		printf("Initial power = 2^(%u+256) = 2^%u mod q' = %s\n",leadb,pow,&char_buf[convert_uint256_base10_char(char_buf,x)]);
		printf("Looping over %u remaining bits in power:\n",start_index);
	}
#endif
	// LR binary powering loop:
	for(j = start_index-1; j >= 0; j--) {
		curr_bit = TEST_BIT256(pshift, j);
		SQR_LOHI256(x,lo,hi);	// x^2 mod q is returned in x
		MULL256(lo,qinv,lo);
		MULH256(q,lo,lo);
	#if FAC_DEBUG
		if(dbg) { pow = 2*pow + curr_bit - 256; printf("\tJ = %2u: [bit = %u]pow = %u, x = %s\n",j,curr_bit,pow,&char_buf[convert_uint256_base10_char(char_buf,x)]); }
	#endif
		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT256(hi, lo)) {
			SUB256(q, lo, lo);	ADD256(lo, hi, x);
		} else {
			SUB256(hi, lo, x);
		}
		if(curr_bit) {	// Combines overflow-on-add and need-to-subtract-q-from-sum checks:
			if(CMPUGT256(x, qhalf)){ ADD256(x, x, x); SUB256(x, q, x); }else{ ADD256(x, x, x); }
		}
	}
	// Since pre-subtracted lg2(R) = 256 from computed powermod exponent, no need to un-scale the loop output.
#if FAC_DEBUG
	if(dbg) printf("pow = %u, x = %s\n",pow,&char_buf[convert_uint256_base10_char(char_buf,x)]);
#endif
	// If we applied an initial right-justify shift to the modulus, restore the shift to the
	// current (partial) remainder and re-add the off-shifted part of the true remainder.
	if(nshift) {
		LSHIFT256(x,nshift,x);
	#if FAC_DEBUG
		if(dbg) printf("Restoring power-of-2: pow = %u, x *= 2^%u = %s\n",pow+nshift,nshift, &char_buf[convert_uint256_base10_char(char_buf, x)]);
	#endif
	}
#if FAC_DEBUG
	if(dbg) printf("xout = %s\n", &char_buf[convert_uint256_base10_char(char_buf, x)]);
#endif
	return x;
}

/*
Function to find 2^(-p) mod q, where p and q are both 256-bit unsigned integers.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^256 (i.e. our MODQ
operation effects multiply modulo 2^256).

The key 3-operation sequence here is as follows:

	SQR_LOHI256(x,lo,hi);	// Input   x has 256 bits; lo/hi have 256 bits
	MULL256(lo,qinv,lo);	// Inputs lo & qinv, and output (overwrites lo) have 256 bits
	MULH256(q,lo,lo);	// Inputs  q &   lo, and output (overwrites lo) have 256 bits
*/
uint256 twopmodq256(uint256 p, uint256 q)
{
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead8, lo64;
	uint256 qhalf, qinv, x, lo, hi;
	static uint256 psave = {0ull, 0ull, 0ull, 0ull}, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	uint32 FERMAT;
	if(p.d3 != 0ull)
		FERMAT = isPow2_64(p.d3) && (p.d2 == 0ull) && (p.d1 == 0ull) && (p.d0 == 0ull);
	else if(p.d2 != 0ull)
		FERMAT = isPow2_64(p.d2) && (p.d1 == 0ull) && (p.d0 == 0ull);
	else if(p.d1 != 0ull)
		FERMAT = isPow2_64(p.d1) && (p.d0 == 0ull);
	else
		FERMAT = isPow2_64(p.d0);
	FERMAT <<= 1;	// *2 is b/c need to add 2 to the usual Mers-mod residue in the Fermat case

	RSHIFT_FAST256(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || !CMPEQ256(p, psave))
	{
		first_entry = FALSE;
		psave  = p;
		x.d0 = 256; x.d1 = x.d2 = x.d3 = 0;
		ADD256(p, x, pshift);
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

		zshift = 255 - lead8;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		pshift.d3 = ~pshift.d3;	pshift.d2 = ~pshift.d2;	pshift.d1 = ~pshift.d1;	pshift.d0 = ~pshift.d0;
	}

	// Find inverse (mod 2^256) of q in preparation for modular multiply: q must be odd for Montgomery-style modmul to work:
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

	/* Now that have bottom 64 bits of qinv, do two more Newton iterations using 128-bit and full 256-bit operands, resp: */
	// 64 -> 128 bits: CF. twopmodq128.c: qinv has 128 bits, but only the upper 64 get modified here:
#ifdef MUL_LOHI64_SUBROUTINE
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + __MULH64(q.d0, qinv.d0));
#else
	MULH64(q.d0, qinv.d0, lo64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + lo64);
#endif
	// 128 -> 256 bits:
#if 1
	MULH128(q, qinv, x);
	// NB: input order here matters! Low-half input must be first in arglist:
	MULL128_256lohihalf(qinv, q, lo);	//      MULL128([hi half of q], [lo half of qinv])
	ADD128(x, lo, x);					// x += MULL128_256hilohalf result [only need low 128 bits]
	MULL128(qinv, x, hi);				// MULL128(qinv, x)
	qinv.d2 = ~hi.d0;	qinv.d3 = ~hi.d1;			// [hi half of qinv] = -MULL128(qinv, x) ...
	qinv.d2 += 1ull; qinv.d3 += (qinv.d2 == 0ull);	// ... Negation here uses multiword analog of the identity -x = ~x + 1 .
#else
	MULL256(q, qinv, x);
	SUB256 (TWO256, x, x);
	MULL256(qinv, x, qinv);
#endif

	/* Since zstart is a power of two < 2^256, use a streamlined code sequence for the first iteration: */
	ASSERT(start_index>=2, "twopmodq256 : start_index < 2!");
	j = start_index-1;

	/* MULL256(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT256(qinv, zshift, lo);
	MULH256(q,lo,lo);
	SUB256(q, lo, x);	// hi = 0 in this instance, which simplifies things.

	if(TEST_BIT256(pshift, j)) {
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT256(x, qhalf)){ ADD256(x, x, x); SUB256(x, q, x); }else{ ADD256(x, x, x); }
	}

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		SQR_LOHI256(x,lo,hi);
		MULL256(lo,qinv,lo);
		MULH256(q,lo,lo);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPULT256(hi, lo)) {
			SUB256(q, lo, lo);
			ADD256(lo, hi, x);
		} else {
			SUB256(hi, lo, x);
		}

		if(TEST_BIT256(pshift, j)) {
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT256(x, qhalf)){ ADD256(x, x, x); SUB256(x, q, x); }else{ ADD256(x, x, x); }
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	ADD256(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^255, so x + x cannot overflow. */
	// For Fn with n > 50-or-so it is not uncommon to have q = b*2^64 + 1, thus need to check for borrow
	lo.d3 = lo.d2 = lo.d1 = 0ull; lo.d0 = FERMAT;	SUB256(q,lo,q);	// Since carry may propagate > 1 word, use general SUB
	SUB256(x,q,x);
	return x;
}
/* Sep 2015: timing comparisons on my clunky old 32-bit Core Duo:

192-bit:
	F(147) has 1 factors in range k = [245044800, 261381120], passes 0-15
	Performed 580390 trial divides
	real	0m28.437s
	user	0m53.595s

	~170 modmul per q --> ~1075 cycle/modmul [lumps sieve time in w/modmul]

256-bit:
	F(195) has 1 factors in range k = [97190681427840, 97190697764160], passes 0-15
	Performed 579792 trial divides
	real	1m9.458s
	user	2m10.487s

	~235 modmul per q --> ~1900 cycle/modmul [lumps sieve time in w/modmul]
*/
