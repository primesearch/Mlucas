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

#undef YES_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
#endif

#if FAC_DEBUG
#error !!!
	char char_buf[STR_MAX_LEN], str0[64], str1[64];
#endif

/***********************************************************************************/
/* Calculate x^2%q using Montgomery-style divideless modular multiply algorithm,   */
/* for various-sized operands (ranging from <= 63 bits to 192 bits), in 1,4 and 8- */
/* operands-at-a-time versions to exploit pipelining:                              */
/***********************************************************************************/

/* First, a couple of small functions for testing the modmul: */

uint64 test_modsqr64(uint64 x, uint64 q)
{
	uint32 j;
	uint64 qinv,t,hi,lo;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q & 0x1, "q must be odd!");
	qinv = (q+q+q) ^ (uint64)2;
	for(j = 0; j < 4; j++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

#ifdef MUL_LOHI64_SUBROUTINE
	SQR_LOHI64(x,&lo,&hi);
#else
	SQR_LOHI64(x,lo,hi);
#endif

/*...x*y mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/
	lo *= qinv;
#ifdef MUL_LOHI64_SUBROUTINE
	lo = __MULH64(q,lo);
#else
	MULH64(q,lo,lo);
#endif

	t = hi - lo;

	if(t > hi)
	{
		t += q;	/* had a borrow */
	}

	return t;
}

/***********************************************************************************/

uint96 test_modsqr96(uint96 x, uint96 q)
{
	uint32 j;
	uint96 qinv,t,hi,lo;
#if 1
	uint64 __l,__m,__a,__b;
	uint32 __tt = x.d1, __hl32,__hh32;
	MUL64x32(x.d0,__tt, __a, __b);
	SQR_LOHI64(x.d0,     __l, __m);
	MUL_LOHI32(__tt,__tt,__hl32,__hh32);
#endif

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q.d0 & 0x1, "q must be odd!");
	/* Init qinv = q. Since we're only interested in the bottom 3 bits of q, can use 64-bit math for that:*/
	qinv.d0 = q.d0;	qinv.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*((uint64)2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 3,6,12,24,48,96. */
	for(j = 0; j < 4; j++)
	{
		qinv.d0 = qinv.d0*((uint64)2 - q.d0*qinv.d0);
	}

	/* Since #SBs < 64 until last iteration, only need 96-bit mod-2^96 math on that one. */
	MULL96(q, qinv, t);
	SUB96 (TWO96, t, t);
	MULL96(qinv, t, qinv);

/*...x*y mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/
	SQR_LOHI96(x,lo,hi);
	MULL96(lo,qinv,lo);
	MULH96(q,lo,lo);

	SUB96(hi,lo,t);

	if(CMPULT96(hi,t))
	{
		ADD96(t,q,t);	/* had a borrow */
	}

	return t;
}

/***********************************************************************************/

uint128 test_modsqr128(uint128 x, uint128 q)
{
	uint32 j;
	uint128 qinv,t,hi,lo;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q.d0 & 0x1, "q must be odd!");
	/* Init qinv = q. Since we're only interested in the bottom 3 bits of q, can use 64-bit math for that:*/
	qinv.d0 = (q.d0+q.d0+q.d0) ^ (uint64)2;
	qinv.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*((uint64)2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed. */
	for(j = 0; j < 4; j++)
	{
		qinv.d0 = qinv.d0*((uint64)2 - q.d0*qinv.d0);
	}

	/* Since #SBs <= 64 until last iteration, only need 128-bit mod-2^128 math on that one. */
	MULL128(q, qinv, t);
	SUB128 (TWO128, t, t);
	MULL128(qinv, t, qinv);

/*...x*y mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/
	SQR_LOHI128(x,lo,hi);
	MULL128(lo,qinv,lo);
	MULH128(q,lo,lo);

	SUB128(hi,lo,t);

	if(CMPULT128(hi,t))
	{
		ADD128(t,q,t);	/* had a borrow */
	}

	return t;
}

/***********************************************************************************/

uint128 test_modsqr128_96(uint128 x, uint128 q)
{
	uint32 j;
	uint64 hi64;
	uint128 qinv,t,lo;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q.d0 & 0x1, "q must be odd!");
	/* Init qinv = q. Since we're only interested in the bottom 3 bits of q, can use 64-bit math for that:*/
	qinv.d0 = (q.d0+q.d0+q.d0) ^ (uint64)2;
	qinv.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*((uint64)2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed. */
	for(j = 0; j < 4; j++)
	{
		qinv.d0 = qinv.d0*((uint64)2 - q.d0*qinv.d0);
	}

	/* Since #SBs <= 64 until last iteration, only need 128-bit mod-2^128 math on that one. */
	MULL128(q, qinv, t);
	SUB128 (TWO128, t, t);
	MULL128(qinv, t, qinv);	/* qinv has 128 bits */

	/*...x*y mod q is returned in x. */
	SQR_LOHI128_96(x,lo,hi64);	/* lo has 128 bits, hi64 has 64. */
	MULL128(lo,qinv,lo);		/* lo has 128 bits */
	MULH128(q,lo,lo);		/* lo has  96 bits */

	/* Calculate h-l+q. Since the odds of h > l are at most 1 in 2^32, we assume q-l+h < q, and calculate that sans branches. */
	SUB128(q,lo,t);
	t.d0 +=  hi64;
	t.d1 += (t.d0 < hi64);	/* Had a carry. */

	return t;
}


/*
!   We calculate 2^(-p) rather than 2^p since Montgomery's MULH64-based modular
!   multiply algorithm, used here, lends itself to negative modular exponents.
!   Uses a modified Left-to-Right binary method of exponentiation similar to
!   that used in evaluating Lucas sequences.
!
!   Example: perform a Fermat base-2 compositeness test on q = 1000000007
!   (a known base-2 pseudoprime, i.e. the result should = 1). Using the inverse-power
!   scheme, the powering 2^(-(q-1)) proceeds as follows: we first calculate the
!   biased exponent q-1+64, since Montgomery's modular multiply using base 2^64
!   yields a result that needs to be multiplied by 2^64 to get the desired result.
!   This scheme uses the property that (using base 2^64 for the M-multiply) if
!   the multiplicands x and y are scaled by 2^64, i.e. we do the Montgomery multiply
!   of x*2^64 and y*2^64, the output is also scaled by 2^64, i.e. = (x*y mod q)*2^64.
!   This is why we add 64 to the exponent of the binary powering - if we scale the
!   initial inputs by 2^64, this scaling factor is carried through to the end of the
!   powering calculation, which thus yields 2^-(input power)*2^64. Thus, setting
!   (input power) = E + 64, where E = (desired exponent of 2), we wind up with
!   2^-(E + 64) * 2^64 = 2^-E mod q. For our example, E = q-1 and so we have
!
!   (input power) = q-1+64 1000000070, which in binary is = 111011100110101100101001000110.
!
!   We use a left-to-right binary exponentiation, and get the leading (leftmost)
!   6 bits basically for free, using just a shift. These are 111011 = 59, i.e. we
!   start with 2^-59 == 2^(64-59) = 2^5 modulo 2^64. The remaining bits then get
!   processed as follows (note that the input to each modular multiply is ,
!   so instead of doing
!
!	- square inputs, divide by 2 if the current bit = 1,
!
!   we do the following, which uses only modular muls by 2, not divides:
!
!	- if current bit = 0, multiply one input by 2
!
!   Bit				Power computed using squarings and mults (the latter if the next bit to be processed = 1):
!   -	2^64*2^(-        59)	Start with 2^(-60) (i.e. 2^4). At each stage we'll get the desired (negative) exponent minus one:
!   1	2^64*2^(-       119)	-120		square       2^( -60)
!   0	2^64*2^(-       238)	-239		mul 2^(-119)*2^(-120)
!   0	2^64*2^(-       476)	-477		mul 2^(-238)*2^(-239)
!   1	2^64*2^(-       953)	-954		square       2^(-477)
!   1	2^64*2^(-      1907) 	...
!   0	2^64*2^(-      3814) 	i.e. if current bit = 1, square inputs, otherwise square and double.
!   1	2^64*2^(-      7629)
!   0	2^64*2^(-     15258)
!   1	2^64*2^(-     30517)
!   1	2^64*2^(-     61035)
!   0	2^64*2^(-    122070)
!   0	2^64*2^(-    244140)
!   1	2^64*2^(-    488281)
!   0	2^64*2^(-    976562)
!   1	2^64*2^(-   1953125)
!   0	2^64*2^(-   3906250)
!   0	2^64*2^(-   7812500)
!   1	2^64*2^(-  15625001)
!   0	2^64*2^(-  31250002)
!   0	2^64*2^(-  62500004)
!   0	2^64*2^(- 125000008)
!   1	2^64*2^(- 250000017)
!   1	2^64*2^(- 500000035)
!   0	2^64*2^(-1000000070)	-1000000071, one less than the desired result, so double and return.
*/

/***********************************************************************************/
/*** 63-BIT INPUTS *****************************************************************/
/***********************************************************************************/

uint64 twopmodq63(uint64 p, uint64 q)
{
	 int32 j;
	uint64 qinv, x, lo, hi, r;
	static uint64 psave = 0, pshift, zstart;
	static uint32 start_index, first_entry = TRUE;

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 6 to account for the fact that we can do the powering for the leftmost
	!    6 bits via a simple shift.
	*/
		start_index = 64-leadz64(pshift)-6;	/* Leftward bit at which to start the l-r binary powering, assuming
							 the leftmost 6 bits have already been processed via a shift (see next). */
		zstart = ((uint64)1) << (63-ibits64(pshift,start_index,6));	/* no need to mod this, since ibits64(...) > 0. */
		pshift = ~pshift;
	}

	/*
	Find modular inverse (mod 2^64) of q in preparation for modular multiply.
	We use the simple and elegant iterative inversion method of Montgomery,
	which amounts to a modular analogue of Newton's method for iterative inversion:

	     Zinv = Z                   ! Z*Zinv == 1 (mod 2^3)
	     Zinv = Zinv*(2 - Z*Zinv)   ! Z*Zinv == 1 (mod 2^6)
	     Zinv = Zinv*(2 - Z*Zinv)   ! Z*Zinv == 1 (mod 2^12)
	     Zinv = Zinv*(2 - Z*Zinv)   ! Z*Zinv == 1 (mod 2^24)
	     Zinv = Zinv*(2 - Z*Zinv)   ! Z*Zinv == 1 (mod 2^48)
	     Zinv = Zinv*(2 - Z*Zinv)   ! Z*Zinv == 1 (mod 2^64)

	(all integer arithmetic, mod 2^64).  The number of correct bits (at the low end)
	doubles at each step. We choose a different starting value of Zinv, XOR(3*Z, 2),
	so the first congruence holds modulo 2^4, thus requiring just four iterations.
	(One can also initialize a table of inverses mod 2^8 or 2^16 and cut the number
	of iterations to 3 or 2, respectively, but that seems more trouble than it's worth.)

	All this works just fine if Z has > 64 bits (working with the bottom 64 bits of Z),
	but the second phase of the Montgomery algorithm, i.e. the MULH and SQR_HIGH operations,
	need modification. Specifically, if multiplicands x and y and modulus q have (as many as)
	64 + E (read E as standing for 'excess') bits, then the upper part of the square x^2
	has as many as 64 +2*E bits, and the result of MULH_EXACT(q, lo*qinv) (where lo*qinv is
	modulo 2^64) has as many as 64 + E bits. If E is small, we can calculate these quantities
	using a call to SQR_LOHI (square modulo 2^128) and SQR_HIGH (upper 64 bits of product modulo 2^128)
	which use hardware multiply functionality, plus some adds and so forth, i.e. using no
	more hardware integer multiplies than for 64-bit inputs. Some examples:

	E = 1 (65-bit inputs, q has precisely 65 bits, no more, no less): let x = X + A*2^64, then

	x^2 = X^2 + 2*A*X*2^64 + A^2*2^128, where A is at most one bit.
	This can be computed via the following sequence:

	SQR_LOHI64(X, lo, hi);	# lo has lower 64 bits of 130-bit square
	if(A > 0)
	{
		Y  = (-A & X);
		hi+= Y; if(hi<Y) ++A;
		hi+= Y; if(hi<Y) ++A;
	}

	The subsequent call to MULH64(q, lo*qinv) generates a product of at most 129 bits,
	so can be gotten via the simpler sequence (where q = Q + 2^64)

	uint64 lq = lo*qinv, B = 0;
	ho = MULH64(Q, lq);
	ho += lq; if(ho < lq) ++B;	# Upper (65th) bit of q guaranteed to be 1.

	And now we need to get SQR_HIGH_EXACT(x) - MULH_EXACT(q, lo*qinv), normalizing as we go.
	The full implementation may be found in the 65-bit version of the function.

	*/
	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q & 0x1, "q must be odd!");

	qinv = (q+q+q) ^ (uint64)2;
	for(j = 0; j < 4; j++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

	/*...Initialize the binary powering...*/
	x = zstart;
	DBG_ASSERT(HERE, x < q, "twopmodq63 : x0 < q");

#if FAC_DEBUG
	fprintf(stderr, "twopmodq63 : x0 = %s, q = %s\n", &str0[convert_uint64_base10_char(str0, x)], &str1[convert_uint64_base10_char(str1, q)] );
#endif

	for(j = start_index-1; j >= 0; j--)
	{
	#ifdef MUL_LOHI64_SUBROUTINE
		SQR_LOHI64(x,&lo,&hi);
	#else
		SQR_LOHI64(x,lo,hi);
	#endif

	/*...x^2 mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/

		MULL64(qinv,lo,lo);
	#ifdef MUL_LOHI64_SUBROUTINE
		lo = __MULH64(q,lo);
	#else
		MULH64(q,lo,lo);
	#endif
		x = hi - lo + q;
		if(x >= q) x -= q;

#if FAC_DEBUG
	fprintf(stderr, "twopmodq63 : j = %d, double = %1u, x = %s\n", j, (uint32)((pshift >> j) & (uint64)1), &str0[convert_uint64_base10_char(str0, x)] );
#endif
		if((pshift >> j) & (uint64)1)
		{
			x = x + x;
			if(x >= q) x -= q;
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q
	implies divisibility, in which case x = (q+1)/2.
	*/
	r = x+x-q;	/* In the case of interest, x = (q+1)/2 < 2^63, so x + x cannot overflow. */
	return r;
}

/*** 4-trial-factor version ***/
uint64 twopmodq63_q4(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
	 int32 j;
	uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3
		, qinv0, qinv1, qinv2, qinv3
		, x0, x1, x2, x3
		, y0, y1, y2, y3
		, lo0, lo1, lo2, lo3
		, hi0, hi1, hi2, hi3, r;
	static uint64 psave = 0, pshift, zstart;
	static uint32 start_index, first_entry = TRUE;
#ifdef MOD_INI_Q4
	uint32 ql0,ql1,ql2,ql3,qh0,qh1,qh2,qh3,qil0,qil1,qil2,qil3,qih0,qih1,qih2,qih3;
#endif

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 6 to account for the scaling factor of 2^64...
	*/
		start_index = 64-leadz64(pshift)-6;	/* Leftward bit at which to start the l-r binary powering, assuming
								 the leftmost 6 bits have already been processed via a shift (see next). */
		zstart = ((uint64)1) << (63-ibits64(pshift,start_index,6));	/* no need to mod this, since ibits64(...) > 0. */
		pshift = ~pshift;
	}
	*checksum1 += q0 + q1 + q2 + q3;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 , "even modulus!");
	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;

	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
	}

#ifdef MOD_INI_Q4
	MOD_INI_Q4(
	 q0,qinv0
	,q1,qinv1
	,q2,qinv2
	,q3,qinv3
	);
#endif

/*...Initialize the binary powering...*/

	x0 = x1 = x2 = x3 = zstart;
#if FAC_DEBUG
  #if 0	/* These appear to be benign: */
	if(x0 >= q0){ sprintf(char_buf, "twopmodq63_q4 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint64_base10_char(str0, x0)], &str1[convert_uint64_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x1 >= q1){ sprintf(char_buf, "twopmodq63_q4 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint64_base10_char(str0, x1)], &str1[convert_uint64_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x2 >= q2){ sprintf(char_buf, "twopmodq63_q4 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint64_base10_char(str0, x2)], &str1[convert_uint64_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x3 >= q3){ sprintf(char_buf, "twopmodq63_q4 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint64_base10_char(str0, x3)], &str1[convert_uint64_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
  #endif
#endif

	for(j = start_index-1; j >= 0; j--)
	{
	#ifdef MOD_INI_Q4
		MOD_SQR_Q4(
		 x0,hi0,y0
		,x1,hi1,y1
		,x2,hi2,y2
		,x3,hi3,y3
		);
	#else
	  #ifdef MUL_LOHI64_SUBROUTINE
		SQR_LOHI64(x0,&lo0,&hi0);
		SQR_LOHI64(x1,&lo1,&hi1);
		SQR_LOHI64(x2,&lo2,&hi2);
		SQR_LOHI64(x3,&lo3,&hi3);
	  #else
		SQR_LOHI64(x0,lo0,hi0);
		SQR_LOHI64(x1,lo1,hi1);
		SQR_LOHI64(x2,lo2,hi2);
		SQR_LOHI64(x3,lo3,hi3);
	 #endif

	/*...x^2 mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/
	/*
		lo0 *= qinv0;
		lo1 *= qinv1;
		lo2 *= qinv2;
		lo3 *= qinv3;
	*/
		MULL64(lo0,qinv0,lo0);
		MULL64(lo1,qinv1,lo1);
		MULL64(lo2,qinv2,lo2);
		MULL64(lo3,qinv3,lo3);

	  #ifdef MUL_LOHI64_SUBROUTINE
		y0 = __MULH64(q0,lo0);
		y1 = __MULH64(q1,lo1);
		y2 = __MULH64(q2,lo2);
		y3 = __MULH64(q3,lo3);
	  #else
		MULH64(q0,lo0,y0);
		MULH64(q1,lo1,y1);
		MULH64(q2,lo2,y2);
		MULH64(q3,lo3,y3);
	  #endif
	#endif
		x0 = hi0 - y0 + q0;
		x1 = hi1 - y1 + q1;
		x2 = hi2 - y2 + q2;
		x3 = hi3 - y3 + q3;

	#ifdef NOBRANCH
		x0 -= q0 & -(x0 >= q0);
		x1 -= q1 & -(x1 >= q1);
		x2 -= q2 & -(x2 >= q2);
		x3 -= q3 & -(x3 >= q3);
	#else
		if(x0 >= q0) x0 -= q0;
		if(x1 >= q1) x1 -= q1;
		if(x2 >= q2) x2 -= q2;
		if(x3 >= q3) x3 -= q3;
	#endif

		if((pshift >> j) & (uint64)1)
		{
			x0 = x0 + x0;
			x1 = x1 + x1;
			x2 = x2 + x2;
			x3 = x3 + x3;
		#ifdef NOBRANCH
			x0 -= q0 & -(x0 >= q0);
			x1 -= q1 & -(x1 >= q1);
			x2 -= q2 & -(x2 >= q2);
			x3 -= q3 & -(x3 >= q3);
		#else
			if(x0 >= q0) x0 -= q0;
			if(x1 >= q1) x1 -= q1;
			if(x2 >= q2) x2 -= q2;
			if(x3 >= q3) x3 -= q3;
		#endif
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	*checksum2 += x0 + x1 + x2 + x3;

	r = 0;
	if(x0+x0-q0 == 1)r += 1;
	if(x1+x1-q1 == 1)r += 2;
	if(x2+x2-q2 == 1)r += 4;
	if(x3+x3-q3 == 1)r += 8;
	return r;
}

/*** 8-trial-factor version ***/
uint64 twopmodq63_q8(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
	 int32 j;
	uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3, q4 = 1+(p<<1)*k4, q5 = 1+(p<<1)*k5, q6 = 1+(p<<1)*k6, q7 = 1+(p<<1)*k7
		, qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
		, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7, r;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 6 to account for the scaling factor of 2^64...
	*/
		start_index = 64-leadz64(pshift)-6;	/* Leftward bit at which to start the l-r binary powering, assuming
							 the leftmost 6 bits have already been processed via a shift (see next). */
		zshift = 63-ibits64(pshift,start_index,6);	/* Since 6 bits with leftmost bit = 1 is guaranteed
								to be in [32,63], the shift count here is in [0, 31].
								That means that zstart < 2^32. Together with the fact that
								squaring a power of two gives another power of two, we can
								simplify the modmul code sequence for the first iteration.
								Every little bit counts (literally in this case :), right? */
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}
	*checksum1 += q0 + q1 + q2 + q3 + q4 + q5 + q6 + q7;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 && q4 & 1 && q5 & 1 && q6 & 1 && q7 & 1 , "even modulus!");
	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;
	qinv4 = (q4+q4+q4) ^ (uint64)2;
	qinv5 = (q5+q5+q5) ^ (uint64)2;
	qinv6 = (q6+q6+q6) ^ (uint64)2;
	qinv7 = (q7+q7+q7) ^ (uint64)2;

	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
		qinv4 = qinv4*((uint64)2 - q4*qinv4);
		qinv5 = qinv5*((uint64)2 - q5*qinv5);
		qinv6 = qinv6*((uint64)2 - q6*qinv6);
		qinv7 = qinv7*((uint64)2 - q7*qinv7);
	}

	/* Since zstart is a power of two < 2^64, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL64(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	lo0 = qinv0 << zshift;
	lo1 = qinv1 << zshift;
	lo2 = qinv2 << zshift;
	lo3 = qinv3 << zshift;
	lo4 = qinv4 << zshift;
	lo5 = qinv5 << zshift;
	lo6 = qinv6 << zshift;
	lo7 = qinv7 << zshift;

#ifdef MUL_LOHI64_SUBROUTINE
	lo0 = __MULH64(q0,lo0);
	lo1 = __MULH64(q1,lo1);
	lo2 = __MULH64(q2,lo2);
	lo3 = __MULH64(q3,lo3);
	lo4 = __MULH64(q4,lo4);
	lo5 = __MULH64(q5,lo5);
	lo6 = __MULH64(q6,lo6);
	lo7 = __MULH64(q7,lo7);
#else
	MULH64(q0,lo0,lo0);
	MULH64(q1,lo1,lo1);
	MULH64(q2,lo2,lo2);
	MULH64(q3,lo3,lo3);
	MULH64(q4,lo4,lo4);
	MULH64(q5,lo5,lo5);
	MULH64(q6,lo6,lo6);
	MULH64(q7,lo7,lo7);
#endif
	/* hi = 0 in this instance, which simplifies things. */
	x0 = q0 - lo0;
	x1 = q1 - lo1;
	x2 = q2 - lo2;
	x3 = q3 - lo3;
	x4 = q4 - lo4;
	x5 = q5 - lo5;
	x6 = q6 - lo6;
	x7 = q7 - lo7;

	if((pshift >> j) & (uint64)1)
	{
		x0 = x0 + x0;
		x1 = x1 + x1;
		x2 = x2 + x2;
		x3 = x3 + x3;
		x4 = x4 + x4;
		x5 = x5 + x5;
		x6 = x6 + x6;
		x7 = x7 + x7;
#ifdef NOBRANCH
		x0 -= q0 & -(x0 >= q0);
		x1 -= q1 & -(x1 >= q1);
		x2 -= q2 & -(x2 >= q2);
		x3 -= q3 & -(x3 >= q3);
		x4 -= q4 & -(x4 >= q4);
		x5 -= q5 & -(x5 >= q5);
		x6 -= q6 & -(x6 >= q6);
		x7 -= q7 & -(x7 >= q7);
#else
		if(x0 >= q0) x0 -= q0;
		if(x1 >= q1) x1 -= q1;
		if(x2 >= q2) x2 -= q2;
		if(x3 >= q3) x3 -= q3;
		if(x4 >= q4) x4 -= q4;
		if(x5 >= q5) x5 -= q5;
		if(x6 >= q6) x6 -= q6;
		if(x7 >= q7) x7 -= q7;
#endif
	}

#if FAC_DEBUG
  #if 0	/* These appear to be benign: */
	if(x0 >= q0){ sprintf(char_buf, "twopmodq63_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint64_base10_char(str0, x0)], &str1[convert_uint64_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x1 >= q1){ sprintf(char_buf, "twopmodq63_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint64_base10_char(str0, x1)], &str1[convert_uint64_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x2 >= q2){ sprintf(char_buf, "twopmodq63_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint64_base10_char(str0, x2)], &str1[convert_uint64_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x3 >= q3){ sprintf(char_buf, "twopmodq63_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint64_base10_char(str0, x3)], &str1[convert_uint64_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x4 >= q4){ sprintf(char_buf, "twopmodq63_q8 : (x4 = %s) >= (q4 = %s)", &str0[convert_uint64_base10_char(str0, x4)], &str1[convert_uint64_base10_char(str1, q4)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x5 >= q5){ sprintf(char_buf, "twopmodq63_q8 : (x5 = %s) >= (q5 = %s)", &str0[convert_uint64_base10_char(str0, x5)], &str1[convert_uint64_base10_char(str1, q5)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x6 >= q6){ sprintf(char_buf, "twopmodq63_q8 : (x6 = %s) >= (q6 = %s)", &str0[convert_uint64_base10_char(str0, x6)], &str1[convert_uint64_base10_char(str1, q6)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x7 >= q7){ sprintf(char_buf, "twopmodq63_q8 : (x7 = %s) >= (q7 = %s)", &str0[convert_uint64_base10_char(str0, x7)], &str1[convert_uint64_base10_char(str1, q7)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
  #endif
#endif

	for(j = start_index-2; j >= 0; j--)
	{
#ifdef MUL_LOHI64_SUBROUTINE

		SQR_LOHI64(x0,&lo0,&hi0);	/*FSQR_LOHI64(x0,lo7,hi7);	if(lo0 != lo7 || hi0 != hi7){ printf("FSQR_LOHI mismatch! x = %s, sqr(hi,lo) = (%s, %s), fsqr = (%s, %s)\n", x0,lo0,hi0,lo7,hi7); } */
		SQR_LOHI64(x1,&lo1,&hi1);
		SQR_LOHI64(x2,&lo2,&hi2);
		SQR_LOHI64(x3,&lo3,&hi3);
		SQR_LOHI64(x4,&lo4,&hi4);
		SQR_LOHI64(x5,&lo5,&hi5);
		SQR_LOHI64(x6,&lo6,&hi6);
		SQR_LOHI64(x7,&lo7,&hi7);
		/*
		FSQR_LOHI64(x0,lo0,hi0);
		FSQR_LOHI64(x1,lo1,hi1);
		FSQR_LOHI64(x2,lo2,hi2);
		FSQR_LOHI64(x3,lo3,hi3);
		FSQR_LOHI64(x4,lo4,hi4);
		FSQR_LOHI64(x5,lo5,hi5);
		FSQR_LOHI64(x6,lo6,hi6);
		FSQR_LOHI64(x7,lo7,hi7);
		*/
#else
								/* Uncomment these to test out various versions of the FSQR... routines: */
		SQR_LOHI64(x0,lo0,hi0);	/*FSQR_LOHI64(x0,lo7,hi7);	if(lo0 != lo7 || hi0 != hi7){ printf("FSQR_LOHI mismatch! x = %s, sqr(lo,hi) = (%s, %s), fsqr = (%s, %s)\n", x0,lo0,hi0,lo7,hi7); abort(); } */
		SQR_LOHI64(x1,lo1,hi1);	/*FSQR_LOHI64(x1,lo7,hi7);	if(lo1 != lo7 || hi1 != hi7){ printf("FSQR_LOHI mismatch! x = %s, sqr(lo,hi) = (%s, %s), fsqr = (%s, %s)\n", x1,lo1,hi1,lo7,hi7); abort(); } */
		SQR_LOHI64(x2,lo2,hi2);	/*FSQR_LOHI64(x2,lo7,hi7);	if(lo2 != lo7 || hi2 != hi7){ printf("FSQR_LOHI mismatch! x = %s, sqr(lo,hi) = (%s, %s), fsqr = (%s, %s)\n", x2,lo2,hi2,lo7,hi7); abort(); } */
		SQR_LOHI64(x3,lo3,hi3);	/*FSQR_LOHI64(x3,lo7,hi7);	if(lo3 != lo7 || hi3 != hi7){ printf("FSQR_LOHI mismatch! x = %s, sqr(lo,hi) = (%s, %s), fsqr = (%s, %s)\n", x3,lo3,hi3,lo7,hi7); abort(); } */
		SQR_LOHI64(x4,lo4,hi4);	/*FSQR_LOHI64(x4,lo7,hi7);	if(lo4 != lo7 || hi4 != hi7){ printf("FSQR_LOHI mismatch! x = %s, sqr(lo,hi) = (%s, %s), fsqr = (%s, %s)\n", x4,lo4,hi4,lo7,hi7); abort(); } */
		SQR_LOHI64(x5,lo5,hi5);	/*FSQR_LOHI64(x5,lo7,hi7);	if(lo5 != lo7 || hi5 != hi7){ printf("FSQR_LOHI mismatch! x = %s, sqr(lo,hi) = (%s, %s), fsqr = (%s, %s)\n", x5,lo5,hi5,lo7,hi7); abort(); } */
		SQR_LOHI64(x6,lo6,hi6);	/*FSQR_LOHI64(x6,lo7,hi7);	if(lo6 != lo7 || hi6 != hi7){ printf("FSQR_LOHI mismatch! x = %s, sqr(lo,hi) = (%s, %s), fsqr = (%s, %s)\n", x6,lo6,hi6,lo7,hi7); abort(); } */
		SQR_LOHI64(x7,lo7,hi7);
		/*
		FSQR_LOHI_q8(
		  x0,lo0,hi0
		, x1,lo1,hi1
		, x2,lo2,hi2
		, x3,lo3,hi3
		, x4,lo4,hi4
		, x5,lo5,hi5
		, x6,lo6,hi6
		, x7,lo7,hi7);
		//
		FSQR_LOHI25(x0,lo0,hi0);
		FSQR_LOHI25(x1,lo1,hi1);
		FSQR_LOHI25(x2,lo2,hi2);
		FSQR_LOHI25(x3,lo3,hi3);
		FSQR_LOHI25(x4,lo4,hi4);
		FSQR_LOHI25(x5,lo5,hi5);
		FSQR_LOHI25(x6,lo6,hi6);
		FSQR_LOHI25(x7,lo7,hi7);
		*/
#endif
	/*...x^2 mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/

		MULL64(lo0,qinv0,lo0);
		MULL64(lo1,qinv1,lo1);
		MULL64(lo2,qinv2,lo2);
		MULL64(lo3,qinv3,lo3);
		MULL64(lo4,qinv4,lo4);
		MULL64(lo5,qinv5,lo5);
		MULL64(lo6,qinv6,lo6);
		MULL64(lo7,qinv7,lo7);

		x0 = hi0 + q0;
		x1 = hi1 + q1;
		x2 = hi2 + q2;
		x3 = hi3 + q3;
		x4 = hi4 + q4;
		x5 = hi5 + q5;
		x6 = hi6 + q6;
		x7 = hi7 + q7;

#ifdef MUL_LOHI64_SUBROUTINE
		lo0 = __MULH64(q0,lo0);
		lo1 = __MULH64(q1,lo1);
		lo2 = __MULH64(q2,lo2);
		lo3 = __MULH64(q3,lo3);
		lo4 = __MULH64(q4,lo4);
		lo5 = __MULH64(q5,lo5);
		lo6 = __MULH64(q6,lo6);
		lo7 = __MULH64(q7,lo7);
#else
		MULH64(q0,lo0,lo0);
		MULH64(q1,lo1,lo1);
		MULH64(q2,lo2,lo2);
		MULH64(q3,lo3,lo3);
		MULH64(q4,lo4,lo4);
		MULH64(q5,lo5,lo5);
		MULH64(q6,lo6,lo6);
		MULH64(q7,lo7,lo7);
#endif

		x0 -= lo0;
		x1 -= lo1;
		x2 -= lo2;
		x3 -= lo3;
		x4 -= lo4;
		x5 -= lo5;
		x6 -= lo6;
		x7 -= lo7;

#ifdef NOBRANCH
		x0 -= q0 & -(x0 >= q0);
		x1 -= q1 & -(x1 >= q1);
		x2 -= q2 & -(x2 >= q2);
		x3 -= q3 & -(x3 >= q3);
		x4 -= q4 & -(x4 >= q4);
		x5 -= q5 & -(x5 >= q5);
		x6 -= q6 & -(x6 >= q6);
		x7 -= q7 & -(x7 >= q7);
#else
		if(x0 >= q0) x0 -= q0;
		if(x1 >= q1) x1 -= q1;
		if(x2 >= q2) x2 -= q2;
		if(x3 >= q3) x3 -= q3;
		if(x4 >= q4) x4 -= q4;
		if(x5 >= q5) x5 -= q5;
		if(x6 >= q6) x6 -= q6;
		if(x7 >= q7) x7 -= q7;
#endif

		if((pshift >> j) & (uint64)1)
		{
			x0 = x0 + x0;
			x1 = x1 + x1;
			x2 = x2 + x2;
			x3 = x3 + x3;
			x4 = x4 + x4;
			x5 = x5 + x5;
			x6 = x6 + x6;
			x7 = x7 + x7;
#ifdef NOBRANCH
			x0 -= q0 & -(x0 >= q0);
			x1 -= q1 & -(x1 >= q1);
			x2 -= q2 & -(x2 >= q2);
			x3 -= q3 & -(x3 >= q3);
			x4 -= q4 & -(x4 >= q4);
			x5 -= q5 & -(x5 >= q5);
			x6 -= q6 & -(x6 >= q6);
			x7 -= q7 & -(x7 >= q7);
#else
			if(x0 >= q0) x0 -= q0;
			if(x1 >= q1) x1 -= q1;
			if(x2 >= q2) x2 -= q2;
			if(x3 >= q3) x3 -= q3;
			if(x4 >= q4) x4 -= q4;
			if(x5 >= q5) x5 -= q5;
			if(x6 >= q6) x6 -= q6;
			if(x7 >= q7) x7 -= q7;
#endif
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	*checksum2 += x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7;

	r = 0;
	if(x0+x0-q0 == 1)r +=  1;
	if(x1+x1-q1 == 1)r +=  2;
	if(x2+x2-q2 == 1)r +=  4;
	if(x3+x3-q3 == 1)r +=  8;
	if(x4+x4-q4 == 1)r += 16;
	if(x5+x5-q5 == 1)r += 32;
	if(x6+x6-q6 == 1)r += 64;
	if(x7+x7-q7 == 1)r +=128;
	return r;
}

#if 1	/************ experimental code below, needs to be modified from q to k-form ******************/
/* Does an 8-fold base-2 PRP test on the prime candidates q0-7. */
uint64 twopmodq63_x8(uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint64 q4, uint64 q5, uint64 q6, uint64 q7)
{
	 int32 j, retval = 0;
	uint32 start_index;
	uint64 lead0, pshift0, qinv0, zshift0, x0, lo0, hi0;
	uint64 lead1, pshift1, qinv1, zshift1, x1, lo1, hi1;
	uint64 lead2, pshift2, qinv2, zshift2, x2, lo2, hi2;
	uint64 lead3, pshift3, qinv3, zshift3, x3, lo3, hi3;
	uint64 lead4, pshift4, qinv4, zshift4, x4, lo4, hi4;
	uint64 lead5, pshift5, qinv5, zshift5, x5, lo5, hi5;
	uint64 lead6, pshift6, qinv6, zshift6, x6, lo6, hi6;
	uint64 lead7, pshift7, qinv7, zshift7, x7, lo7, hi7;

	DBG_ASSERT(HERE, (q0 < q1) && (q1 < q2) && (q2 < q3) && (q3 < q4) && (q4 < q5) && (q5 < q6) && (q6 < q7), "twopmodq63_x8: Inputs nonmonotone!");

	pshift0 = q0 + 63;
	pshift1 = q1 + 63;
	pshift2 = q2 + 63;
	pshift3 = q3 + 63;
	pshift4 = q4 + 63;
	pshift5 = q5 + 63;
	pshift6 = q6 + 63;
	pshift7 = q7 + 63;

	/* Find number of leading zeros in p, use it to find the position of the leftmost ones bit: */
	j = leadz64(pshift0);
	if( leadz64(pshift7) != j )	/* Fused 8-fold algo needs all p's to have same bitlength */
	{
		retval  = ((uint64)twopmodq63(q0-1, q0) == 1ull);
		retval += ((uint64)twopmodq63(q1-1, q1) == 1ull) << 1;
		retval += ((uint64)twopmodq63(q2-1, q2) == 1ull) << 2;
		retval += ((uint64)twopmodq63(q3-1, q3) == 1ull) << 3;
		retval += ((uint64)twopmodq63(q4-1, q4) == 1ull) << 4;
		retval += ((uint64)twopmodq63(q5-1, q5) == 1ull) << 5;
		retval += ((uint64)twopmodq63(q6-1, q6) == 1ull) << 6;
		retval += ((uint64)twopmodq63(q7-1, q7) == 1ull) << 7;
		return retval;
	}

	/* Extract leftmost 5 bits of pshift and subtract from 32: */
	lead0 = ((pshift0<<j) >> 58);
	lead1 = ((pshift1<<j) >> 58);
	lead2 = ((pshift2<<j) >> 58);
	lead3 = ((pshift3<<j) >> 58);
	lead4 = ((pshift4<<j) >> 58);
	lead5 = ((pshift5<<j) >> 58);
	lead6 = ((pshift6<<j) >> 58);
	lead7 = ((pshift7<<j) >> 58);

	start_index = 64-j-6;	/* Leftward bit at which to start the l-r binary powering, assuming
				 the leftmost 5 bits have already been processed via a shift (see next). */

	/* Doubling the shift count here takes cares of the first SQR_LOHI */
	zshift0 = 63 - lead0;	zshift0 <<= 1;	pshift0 = ~pshift0;
	zshift1 = 63 - lead1;	zshift1 <<= 1;	pshift1 = ~pshift1;
	zshift2 = 63 - lead2;	zshift2 <<= 1;	pshift2 = ~pshift2;
	zshift3 = 63 - lead3;	zshift3 <<= 1;	pshift3 = ~pshift3;
	zshift4 = 63 - lead4;	zshift4 <<= 1;	pshift4 = ~pshift4;
	zshift5 = 63 - lead5;	zshift5 <<= 1;	pshift5 = ~pshift5;
	zshift6 = 63 - lead6;	zshift6 <<= 1;	pshift6 = ~pshift6;
	zshift7 = 63 - lead7;	zshift7 <<= 1;	pshift7 = ~pshift7;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 && q4 & 1 && q5 & 1 && q6 & 1 && q7 & 1 , "even modulus!");

	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;
	qinv4 = (q4+q4+q4) ^ (uint64)2;
	qinv5 = (q5+q5+q5) ^ (uint64)2;
	qinv6 = (q6+q6+q6) ^ (uint64)2;
	qinv7 = (q7+q7+q7) ^ (uint64)2;

	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
		qinv4 = qinv4*((uint64)2 - q4*qinv4);
		qinv5 = qinv5*((uint64)2 - q5*qinv5);
		qinv6 = qinv6*((uint64)2 - q6*qinv6);
		qinv7 = qinv7*((uint64)2 - q7*qinv7);
	}

	/* Since zstart is a power of two < 2^64, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL64(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	lo0 = qinv0 << zshift0;
	lo1 = qinv1 << zshift1;
	lo2 = qinv2 << zshift2;
	lo3 = qinv3 << zshift3;
	lo4 = qinv4 << zshift4;
	lo5 = qinv5 << zshift5;
	lo6 = qinv6 << zshift6;
	lo7 = qinv7 << zshift7;

	/* lo = MULH64(q, lo): */
#ifdef MUL_LOHI64_SUBROUTINE
	lo0 = __MULH64(q0,lo0);
	lo1 = __MULH64(q1,lo1);
	lo2 = __MULH64(q2,lo2);
	lo3 = __MULH64(q3,lo3);
	lo4 = __MULH64(q4,lo4);
	lo5 = __MULH64(q5,lo5);
	lo6 = __MULH64(q6,lo6);
	lo7 = __MULH64(q7,lo7);
#else
	MULH64(q0,lo0,lo0);
	MULH64(q1,lo1,lo1);
	MULH64(q2,lo2,lo2);
	MULH64(q3,lo3,lo3);
	MULH64(q4,lo4,lo4);
	MULH64(q5,lo5,lo5);
	MULH64(q6,lo6,lo6);
	MULH64(q7,lo7,lo7);
#endif
	/* hi = 0 in this instance, which simplifies things. */
	x0 = q0 - lo0;
	x1 = q1 - lo1;
	x2 = q2 - lo2;
	x3 = q3 - lo3;
	x4 = q4 - lo4;
	x5 = q5 - lo5;
	x6 = q6 - lo6;
	x7 = q7 - lo7;

	if((pshift0 >> j) & (uint64)1){ x0 = x0 + x0;	x0 -= q0 & -(x0 >= q0); }
	if((pshift1 >> j) & (uint64)1){ x1 = x1 + x1;	x1 -= q1 & -(x1 >= q1); }
	if((pshift2 >> j) & (uint64)1){ x2 = x2 + x2;	x2 -= q2 & -(x2 >= q2); }
	if((pshift3 >> j) & (uint64)1){ x3 = x3 + x3;	x3 -= q3 & -(x3 >= q3); }
	if((pshift4 >> j) & (uint64)1){ x4 = x4 + x4;	x4 -= q4 & -(x4 >= q4); }
	if((pshift5 >> j) & (uint64)1){ x5 = x5 + x5;	x5 -= q5 & -(x5 >= q5); }
	if((pshift6 >> j) & (uint64)1){ x6 = x6 + x6;	x6 -= q6 & -(x6 >= q6); }
	if((pshift7 >> j) & (uint64)1){ x7 = x7 + x7;	x7 -= q7 & -(x7 >= q7); }

#if FAC_DEBUG
	fprintf(stderr, "twopmodq63_q8 : x0 = %s, q = %s\n", &str0[convert_uint64_base10_char(str0, x0)], &str1[convert_uint64_base10_char(str1, q0)] );
  #if 0	/* These appear to be benign: */
	if(x0 >= q0){ sprintf(char_buf, "twopmodq63_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint64_base10_char(str0, x0)], &str1[convert_uint64_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x1 >= q1){ sprintf(char_buf, "twopmodq63_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint64_base10_char(str0, x1)], &str1[convert_uint64_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x2 >= q2){ sprintf(char_buf, "twopmodq63_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint64_base10_char(str0, x2)], &str1[convert_uint64_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x3 >= q3){ sprintf(char_buf, "twopmodq63_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint64_base10_char(str0, x3)], &str1[convert_uint64_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x4 >= q4){ sprintf(char_buf, "twopmodq63_q8 : (x4 = %s) >= (q4 = %s)", &str0[convert_uint64_base10_char(str0, x4)], &str1[convert_uint64_base10_char(str1, q4)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x5 >= q5){ sprintf(char_buf, "twopmodq63_q8 : (x5 = %s) >= (q5 = %s)", &str0[convert_uint64_base10_char(str0, x5)], &str1[convert_uint64_base10_char(str1, q5)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x6 >= q6){ sprintf(char_buf, "twopmodq63_q8 : (x6 = %s) >= (q6 = %s)", &str0[convert_uint64_base10_char(str0, x6)], &str1[convert_uint64_base10_char(str1, q6)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x7 >= q7){ sprintf(char_buf, "twopmodq63_q8 : (x7 = %s) >= (q7 = %s)", &str0[convert_uint64_base10_char(str0, x7)], &str1[convert_uint64_base10_char(str1, q7)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
  #endif
#endif

	for(j = start_index-2; j >= 0; j--)
	{
	#ifdef MUL_LOHI64_SUBROUTINE
		SQR_LOHI64(x0,&lo0,&hi0);
		SQR_LOHI64(x1,&lo1,&hi1);
		SQR_LOHI64(x2,&lo2,&hi2);
		SQR_LOHI64(x3,&lo3,&hi3);
		SQR_LOHI64(x4,&lo4,&hi4);
		SQR_LOHI64(x5,&lo5,&hi5);
		SQR_LOHI64(x6,&lo6,&hi6);
		SQR_LOHI64(x7,&lo7,&hi7);
	#else
		SQR_LOHI64(x0,lo0,hi0);
		SQR_LOHI64(x1,lo1,hi1);
		SQR_LOHI64(x2,lo2,hi2);
		SQR_LOHI64(x3,lo3,hi3);
		SQR_LOHI64(x4,lo4,hi4);
		SQR_LOHI64(x5,lo5,hi5);
		SQR_LOHI64(x6,lo6,hi6);
		SQR_LOHI64(x7,lo7,hi7);
	#endif

	/*...x^2 mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/

		MULL64(lo0,qinv0,lo0);
		MULL64(lo1,qinv1,lo1);
		MULL64(lo2,qinv2,lo2);
		MULL64(lo3,qinv3,lo3);
		MULL64(lo4,qinv4,lo4);
		MULL64(lo5,qinv5,lo5);
		MULL64(lo6,qinv6,lo6);
		MULL64(lo7,qinv7,lo7);

		x0 = hi0 + q0;
		x1 = hi1 + q1;
		x2 = hi2 + q2;
		x3 = hi3 + q3;
		x4 = hi4 + q4;
		x5 = hi5 + q5;
		x6 = hi6 + q6;
		x7 = hi7 + q7;

#ifdef MUL_LOHI64_SUBROUTINE
		lo0 = __MULH64(q0,lo0);
		lo1 = __MULH64(q1,lo1);
		lo2 = __MULH64(q2,lo2);
		lo3 = __MULH64(q3,lo3);
		lo4 = __MULH64(q4,lo4);
		lo5 = __MULH64(q5,lo5);
		lo6 = __MULH64(q6,lo6);
		lo7 = __MULH64(q7,lo7);
#else
		MULH64(q0,lo0,lo0);
		MULH64(q1,lo1,lo1);
		MULH64(q2,lo2,lo2);
		MULH64(q3,lo3,lo3);
		MULH64(q4,lo4,lo4);
		MULH64(q5,lo5,lo5);
		MULH64(q6,lo6,lo6);
		MULH64(q7,lo7,lo7);
#endif

		x0 -= lo0;
		x1 -= lo1;
		x2 -= lo2;
		x3 -= lo3;
		x4 -= lo4;
		x5 -= lo5;
		x6 -= lo6;
		x7 -= lo7;

#if NOBRANCH
		x0 -= q0 & -(x0 >= q0);
		x1 -= q1 & -(x1 >= q1);
		x2 -= q2 & -(x2 >= q2);
		x3 -= q3 & -(x3 >= q3);
		x4 -= q4 & -(x4 >= q4);
		x5 -= q5 & -(x5 >= q5);
		x6 -= q6 & -(x6 >= q6);
		x7 -= q7 & -(x7 >= q7);
#else
		if(x0 >= q0) x0 -= q0;
		if(x1 >= q1) x1 -= q1;
		if(x2 >= q2) x2 -= q2;
		if(x3 >= q3) x3 -= q3;
		if(x4 >= q4) x4 -= q4;
		if(x5 >= q5) x5 -= q5;
		if(x6 >= q6) x6 -= q6;
		if(x7 >= q7) x7 -= q7;
#endif

#if FAC_DEBUG
	fprintf(stderr, "twopmodq63_q8 : j = %d, double = %1u, x = %s\n", j, (uint32)((pshift0 >> j) & (uint64)1), &str0[convert_uint64_base10_char(str0, x0)] );
#endif
		if((pshift0 >> j) & (uint64)1){ x0 = x0 + x0;	x0 -= q0 & -(x0 >= q0); }
		if((pshift1 >> j) & (uint64)1){ x1 = x1 + x1;	x1 -= q1 & -(x1 >= q1); }
		if((pshift2 >> j) & (uint64)1){ x2 = x2 + x2;	x2 -= q2 & -(x2 >= q2); }
		if((pshift3 >> j) & (uint64)1){ x3 = x3 + x3;	x3 -= q3 & -(x3 >= q3); }
		if((pshift4 >> j) & (uint64)1){ x4 = x4 + x4;	x4 -= q4 & -(x4 >= q4); }
		if((pshift5 >> j) & (uint64)1){ x5 = x5 + x5;	x5 -= q5 & -(x5 >= q5); }
		if((pshift6 >> j) & (uint64)1){ x6 = x6 + x6;	x6 -= q6 & -(x6 >= q6); }
		if((pshift7 >> j) & (uint64)1){ x7 = x7 + x7;	x7 -= q7 & -(x7 >= q7); }
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	retval += ((x0 + x0 - q0) == 1)     ;
	retval += ((x1 + x1 - q1) == 1) << 1;
	retval += ((x2 + x2 - q2) == 1) << 2;
	retval += ((x3 + x3 - q3) == 1) << 3;
	retval += ((x4 + x4 - q4) == 1) << 4;
	retval += ((x5 + x5 - q5) == 1) << 5;
	retval += ((x6 + x6 - q6) == 1) << 6;
	retval += ((x7 + x7 - q7) == 1) << 7;
	return(retval);
}
#endif	/*************/

/***********************************************************************************/
/*** 64-BIT INPUTS *****************************************************************/
/***********************************************************************************/

// Conventional positive-power version of twopmodq, returns true mod:
uint64 twopmmodq64(uint64 p, uint64 q)
{
	 int32 j,nshift;
	uint64 qhalf, qinv, x, rsqr;
	static uint64 psave = 0, zstart;
	static uint32 start_index, first_entry = TRUE;

	qhalf  = q>>1;	/* = (q-1)/2, since q odd. */

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		start_index = 64-leadz64(p)-6;
		zstart = (uint64)1 << ibits64(p,start_index,6);
	}

	nshift = trailz64(q);
	if(nshift)
	{
		ASSERT(HERE, p >= nshift, "twopmodq64 : Must add code to explicitly save offshifted low bits of modulus!");
		q >>= nshift;
		p -= nshift;	// Must also right-shift dividend by (nshift) bits; for 2^p this means subtracting nshift from p
	}
	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q & 0x1, "q must be odd!");
	qinv = (q+q+q) ^ (uint64)2;
	for(j = 0; j < 4; j++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

	/*...Initialize the binary powering...*/
	x = zstart;

	rsqr = radix_power64(q,qinv,2);	// Compute R^2 (mod q) in prep. for Mont-mul with initial seed...

	MONT_MUL64(x,rsqr, q,qinv, x);	// x*R (mod q) = MONT_MUL(x,R^2 (mod q),q,qinv)
	for(j = start_index-1; j >= 0; j--)
	{
		MONT_SQR64(x,q,qinv,x);
		if((p >> j) & (uint64)1)
		{
			if(x > qhalf)	/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			{
				x = x + x;
				x -= q;
			}
			else
				x = x + x;
		}
	}
	MONT_UNITY_MUL64(x,q,qinv,x);	// Un-scale the loop output

	// If we applied an initial right-justify shift to the modulus, restore the shift to the
	// current (partial) remainder and re-add the off-shifted part of the true remainder.
	if(nshift)
	{
		x = (x << nshift);// + rem_save;
	}

	return x;
}

#ifndef YES_ASM	/* Use x86_64-optimized asm version if available */

  #ifdef MUL_LOHI64_SUBROUTINE

	#define MONT_SQR64_q4(__x0,__x1,__x2,__x3,__q0,__q1,__q2,__q3,__qinv0,__qinv1,__qinv2,__qinv3,__z0,__z1,__z2,__z3)\
	{\
		uint64 lo0,lo1,lo2,lo3,hi0,hi1,hi2,hi3;					\
		SQR_LOHI64(__x0,&lo0,&hi0);								\
		SQR_LOHI64(__x1,&lo1,&hi1);								\
		SQR_LOHI64(__x2,&lo2,&hi2);								\
		SQR_LOHI64(__x3,&lo3,&hi3);								\
		MULL64(__qinv0,lo0,lo0);								\
		MULL64(__qinv1,lo1,lo1);								\
		MULL64(__qinv2,lo2,lo2);								\
		MULL64(__qinv3,lo3,lo3);								\
		lo0 = __MULH64(__q0,lo0);								\
		lo1 = __MULH64(__q1,lo1);								\
		lo2 = __MULH64(__q2,lo2);								\
		lo3 = __MULH64(__q3,lo3);								\
		/* did we have a borrow from (hi-lo)? */				\
		__z0 = hi0 - lo0 + ((-(int64)(hi0 < lo0)) & __q0);		\
		__z1 = hi1 - lo1 + ((-(int64)(hi1 < lo1)) & __q1);		\
		__z2 = hi2 - lo2 + ((-(int64)(hi2 < lo2)) & __q2);		\
		__z3 = hi3 - lo3 + ((-(int64)(hi3 < lo3)) & __q3);		\
	}

  #else

	#define MONT_SQR64_q4(__x0,__x1,__x2,__x3,__q0,__q1,__q2,__q3,__qinv0,__qinv1,__qinv2,__qinv3,__z0,__z1,__z2,__z3)\
	{\
		uint64 lo0,lo1,lo2,lo3,hi0,hi1,hi2,hi3;					\
		SQR_LOHI64(__x0, lo0, hi0);								\
		SQR_LOHI64(__x1, lo1, hi1);								\
		SQR_LOHI64(__x2, lo2, hi2);								\
		SQR_LOHI64(__x3, lo3, hi3);								\
		MULL64(__qinv0,lo0,lo0);								\
		MULL64(__qinv1,lo1,lo1);								\
		MULL64(__qinv2,lo2,lo2);								\
		MULL64(__qinv3,lo3,lo3);								\
		MULH64(__q0,lo0,lo0);									\
		MULH64(__q1,lo1,lo1);									\
		MULH64(__q2,lo2,lo2);									\
		MULH64(__q3,lo3,lo3);									\
		/* did we have a borrow from (hi-lo)? */				\
		__z0 = hi0 - lo0 + ((-(int64)(hi0 < lo0)) & __q0);		\
		__z1 = hi1 - lo1 + ((-(int64)(hi1 < lo1)) & __q1);		\
		__z2 = hi2 - lo2 + ((-(int64)(hi2 < lo2)) & __q2);		\
		__z3 = hi3 - lo3 + ((-(int64)(hi3 < lo3)) & __q3);		\
	}

  #endif

#endif


// This variant returns the 4 true mods, overwriting the inputs
void twopmmodq64_q4(uint64 p, uint64 *i0, uint64 *i1, uint64 *i2, uint64 *i3, uint64 *qi0, uint64 *qi1, uint64 *qi2, uint64 *qi3)
{
	int32 j;
	uint64 q0 = *i0, q1 = *i1, q2 = *i2, q3 = *i3
	, qinv0, qinv1, qinv2, qinv3
	, x0, x1, x2, x3
	, y0, y1, y2, y3;
	static uint64 psave = 0, zstart;
	static uint32 start_index, first_entry = TRUE;
	
	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		start_index = 64-leadz64(p)-6;
		zstart = (uint64)1 << ibits64(p,start_index,6);
	}

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q1 > 1 && q1 > 1 && q2 > 1 && q3 > 1 , "modulus must be > 1!");
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 , "even modulus!");

	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;
	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
	}
	x0 = x1 = x2 = x3 = zstart;
	y0 = radix_power64(q0,qinv0,2);	// Compute R^2 (mod q) in prep. for Mont-mul with initial seed...
	y1 = radix_power64(q1,qinv1,2);
	y2 = radix_power64(q2,qinv2,2);
	y3 = radix_power64(q3,qinv3,2);
	MONT_MUL64(x0,y0, q0,qinv0, x0);	// x*R (mod q) = MONT_MUL(x,R^2 (mod q),q,qinv)
	MONT_MUL64(x1,y1, q1,qinv1, x1);
	MONT_MUL64(x2,y2, q2,qinv2, x2);
	MONT_MUL64(x3,y3, q3,qinv3, x3);

  #ifndef YES_ASM

	for(j = start_index-1; j >= 0; j--)
	{
		MONT_SQR64_q4(x0,x1,x2,x3,q0,q1,q2,q3,qinv0,qinv1,qinv2,qinv3,y0,y1,y2,y3);
		if((p >> j) & (uint64)1)
		{
		#ifdef NOBRANCH
			x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
			x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
			x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
			x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
		#else
			x0 = y0 + y0;	if(x0 >= q0 || x0 < y0) x0 -= q0;
			x1 = y1 + y1;	if(x1 >= q1 || x1 < y1) x1 -= q1;
			x2 = y2 + y2;	if(x2 >= q2 || x2 < y2) x2 -= q2;
			x3 = y3 + y3;	if(x3 >= q3 || x3 < y3) x3 -= q3;
		#endif
		} else {
			x0 = y0;
			x1 = y1;
			x2 = y2;
			x3 = y3;
		}
	}

  #else

	j = start_index-1;
//printf("twopmodq64_q4 : x1 = %s\n", &str0[convert_uint64_base10_char(str0, x1)] );
//for(j = start_index-1; j >= 0; j--) {
	__asm__ volatile (\
	"/* Load the q0|1|2 into rbx|rsi|rdi, keeping rcx free for loop index and rax:rdx for double-width MULs: */	\n\t"\
		"movq	%[__q0],%%rbx	\n\t"\
		"movq	%[__q1],%%rsi	\n\t"\
		"movq	%[__q2],%%rdi	\n\t"\
		"/* Must load q3 as-needed due to lack of user register to store it */\n\t"\
	"/* Load the x's into r12-15: */	\n\t"\
		"movq	%[__x0],%%r12	\n\t"\
		"movq	%[__x1],%%r13	\n\t"\
		"movq	%[__x2],%%r14	\n\t"\
		"movq	%[__x3],%%r15	\n\t"\
	"/* Pure-ASM loop control: for(j = start_index-2; j >= 0; j--) */\n\t"\
		"movslq	%[__start_index], %%rcx		\n\t"\
		"subq $1,%%rcx						\n\t"\
		"test %%rcx, %%rcx					\n\t"\
		"jl LoopEnd4a		/* Skip if n < 0 */	\n\t"\
	"LoopBeg4a:								\n\t"\
	"/* SQR_LOHI_q4(x*, lo*, hi*): */	\n\t"\
		"movq	%%r12,%%rax	/* x0-3 in r8-11. */\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r8 	/* lo0 */	\n\t"\
		"movq	%%rdx,%%r12	/* hi0 */	\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r9 	/* lo1 */	\n\t"\
		"movq	%%rdx,%%r13	/* hi1 */	\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r10	/* lo2 */	\n\t"\
		"movq	%%rdx,%%r14	/* hi2 */	\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r11	/* lo3 */	\n\t"\
		"movq	%%rdx,%%r15	/* hi3 */	\n\t"\
		"\n\t"\
		"\n\t"\
	"/* MULL_q4((qinv*, lo*, lo*): */	\n\t"\
		"\n\t"\
		"imulq	%[__qinv0],%%r8 	\n\t"\
		"imulq	%[__qinv1],%%r9 	\n\t"\
		"imulq	%[__qinv2],%%r10	\n\t"\
		"imulq	%[__qinv3],%%r11	\n\t"\
		"\n\t"\
	"/* UMULH_q4((q*, lo*, lo*): lo0-3 in r8-11: */	\n\t"\
		"\n\t"\
		"movq	%%r8 ,%%rax	/* lo0 */\n\t"\
		"mulq	%%rbx	/* q0 in rbx; q0*lo0 in rax:rdx */\n\t"\
		"movq	%%rdx,%%r8 	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
		"movq	%%r9 ,%%rax	/* lo1 */\n\t"\
		"mulq	%%rsi	/* q1 in rsi; q1*lo1 in rax:rdx */\n\t"\
		"movq	%%rdx,%%r9 	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
		"movq	%%r10,%%rax	/* lo2 */\n\t"\
		"mulq	%%rdi	/* q2 in rdi; q2*lo2 in rax:rdx */\n\t"\
		"movq	%%rdx,%%r10	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
		"movq	%%r11,%%rax	/* lo3 */\n\t"\
		"mulq	%[__q3]		\n\t"\
		"movq	%%rdx,%%r11	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
	"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
		"subq	%%r8 ,%%r12			\n\t" /* r12 = (h-l) */\
		"leaq (%%r12,%%rbx),%%r8 	\n\t" /* r8  = (h-l)+q */\
		"cmovcq %%r8 ,%%r12	\n\t" /* if carry = true (i.e. h > l), copy source (r8 = h-l+q) into dest (r12), else leave dest = h-l. */\
		"\n\t"\
		"subq	%%r9 ,%%r13			\n\t"\
		"leaq (%%r13,%%rsi),%%r9 	\n\t"\
		"cmovcq %%r9 ,%%r13	\n\t"\
		"\n\t"\
		"subq	%%r10,%%r14			\n\t"\
		"leaq (%%r14,%%rdi),%%r10	\n\t"\
		"cmovcq %%r10,%%r14	\n\t"\
		"\n\t"\
		"movq	%[__q3],%%rdx	/* Re-use q3 several times below, so store in free register */\n\t"\
		"subq	%%r11,%%r15			\n\t"\
		"leaq (%%r15,%%rdx),%%r11	\n\t"\
		"cmovcq %%r11,%%r15	\n\t"\
		"\n\t"\
"/* If current bit of p == 1, double each output modulo q: */	\n\t"\
		"/* if((p >> j) & (uint64)1) { */	\n\t"\
		"movl	%[__p],%%eax	/* Need to follow this with load-j-into-ecx if use HLL loop control in debug mode */\n\t"\
		"shrq	%%cl,%%rax				\n\t"\
		"andq	$0x1,%%rax				\n\t"\
	"je	twopmodq64_q4_pjmp			\n\t"\
		"\n\t"\
		"movq	%%r12,%%r8 	/* r8  <- Copy of x */\n\t"\
		"movq	%%r12,%%rax	/* rax <- Copy of x */\n\t"\
		"leaq (%%r12,%%r12),%%r12	/* x+x */\n\t"\
		"subq	%%rbx,%%r8	/* r8  <- x-q */\n\t"\
		"addq	%%rax,%%r8 	/* r8 <- 2x-q */\n\t"\
		"cmovcq %%r8 ,%%r12	\n\t"/* if carry (i.e. x+x needed modding), copy source (r8 = x+x-q) into dest (r12), else leave dest = x+x. */\
		"\n\t"\
		"movq	%%r13,%%r9 	\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"leaq (%%r13,%%r13),%%r13	\n\t"\
		"subq	%%rsi,%%r9	\n\t"\
		"addq	%%rax,%%r9 	\n\t"\
		"cmovcq %%r9 ,%%r13	\n\t"\
		"\n\t"\
		"movq	%%r14,%%r10	\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"leaq (%%r14,%%r14),%%r14	\n\t"\
		"subq	%%rdi,%%r10	\n\t"\
		"addq	%%rax,%%r10	\n\t"\
		"cmovcq %%r10,%%r14	\n\t"\
		"\n\t"\
		"movq	%%r15,%%r11	\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"leaq (%%r15,%%r15),%%r15	\n\t"\
		"subq	%[__q3],%%r11	/* Weird...using rdx in place of __q3 here made timings 1-2 percent  worse. */\n\t"\
		"addq	%%rax,%%r11	\n\t"\
		"cmovcq %%r11,%%r15	\n\t"\
		"\n\t"\
	"twopmodq64_q4_pjmp:					\n\t"\
		"/* } endif((p >> j) & (uint64)1) */						\n\t"\
		"subq	$1,%%rcx	/* j-- */		\n\t"\
		"cmpq	$0,%%rcx	/* compare j vs 0 */\n\t"\
		"jge	LoopBeg4a	/* if (j >= 0), Loop */	\n\t"\
	"LoopEnd4a:							\n\t"\
		"movq	%%r12,%[__x0]	\n\t"\
		"movq	%%r13,%[__x1]	\n\t"\
		"movq	%%r14,%[__x2]	\n\t"\
		"movq	%%r15,%[__x3]	\n\t"\
		:	/* outputs: none */\
		: [__q0] "m" (q0)	/* All inputs from memory addresses here */\
		 ,[__q1] "m" (q1)	\
		 ,[__q2] "m" (q2)	\
		 ,[__q3] "m" (q3)	\
		 ,[__qinv0] "m" (qinv0)	\
		 ,[__qinv1] "m" (qinv1)	\
		 ,[__qinv2] "m" (qinv2)	\
		 ,[__qinv3] "m" (qinv3)	\
		 ,[__x0] "m" (x0)	\
		 ,[__x1] "m" (x1)	\
		 ,[__x2] "m" (x2)	\
		 ,[__x3] "m" (x3)	\
		 ,[__p] "m" (p)	\
		 ,[__j] "m" (j)	/* Only need this if debug and explicit loop enabled, but try with/sans for each version of the asm, pivk faster one. */\
		 ,[__start_index] "m" (start_index)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"		/* Clobbered registers */\
		);

  #endif	// YES_ASM

	MONT_UNITY_MUL64(x0,q0,qinv0,*i0);	// Un-scale the loop output
	MONT_UNITY_MUL64(x1,q1,qinv1,*i1);
	MONT_UNITY_MUL64(x2,q2,qinv2,*i2);
	MONT_UNITY_MUL64(x3,q3,qinv3,*i3);

	// User needs mod-inverses returned:
	if(qi0) {
		*qi0 = qinv0;
		*qi1 = qinv1;
		*qi2 = qinv2;
		*qi3 = qinv3;
	}
	return;
}

// Negative-power version, needing no explicit radix-mod scalings:
uint64 twopmodq64(uint64 p, uint64 q)
{
	 int32 j;
	uint64 qhalf, qinv, x;
	static uint64 psave = 0, pshift, zstart;
	static uint32 start_index, first_entry = TRUE;

	qhalf  = q>>1;	/* = (q-1)/2, since q odd. */

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
		start_index = 64-leadz64(pshift)-6;

		zstart = ((uint64)1) << (63-ibits64(pshift,start_index,6));
		pshift = ~pshift;
	}

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q & 0x1, "q must be odd!");
	qinv = (q+q+q) ^ (uint64)2;
	for(j = 0; j < 4; j++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

	/*...Initialize the binary powering...*/

	x = zstart;
#if FAC_DEBUG
/*	ASSERT(HERE, x < q, "twopmodq64 : x0 < q");	*/
  #if 0	/* These appear to be benign: */
	if(x >= q){ sprintf(char_buf, "twopmodq64 : (x0 = %s) >= (q = %s)", &str0[convert_uint64_base10_char(str0, x)], &str1[convert_uint64_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
  #endif
#endif

	for(j = start_index-1; j >= 0; j--)
	{
		MONT_SQR64(x,q,qinv,x);
		if((pshift >> j) & (uint64)1)
		{
			if(x > qhalf)	/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			{
				x = x + x;
				x -= q;
			}
			else
				x = x + x;
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
	x = x+x-q;	/* In the case of interest, x = (q+1)/2 < 2^63, so x + x cannot overflow. */
	return x;
}


#ifndef YES_ASM	// YES_ASM versions of the following routines are in twopmodq64_test.c:

/******************************/
/*** 4-trial-factor version ***/
/******************************/

// Negative-power version, needing no explicit radix-mod scalings:
uint64 twopmodq64_q4(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
	 int32 j;
	uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3
		, qinv0, qinv1, qinv2, qinv3
		, x0, x1, x2, x3
		, y0, y1, y2, y3
		, lo0, lo1, lo2, lo3, r;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
		start_index = 64-leadz64(pshift)-6;

		zshift = 63-ibits64(pshift,start_index,6);	/* Since 6 bits with leftmost bit = 1 is guaranteed
								to be in [32,63], the shift count here is in [0, 31].
								That means that zstart < 2^32. Together with the fact that
								squaring a power of two gives another power of two, we can
								simplify the modmul code sequence for the first iteration.
								Every little bit counts (literally in this case :), right? */
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}
	*checksum1 += q0 + q1 + q2 + q3;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 , "even modulus!");

	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;

	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
	}

	/* Since zstart is a power of two < 2^128, use a
	streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL64(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	lo0 = qinv0 << zshift;
	lo1 = qinv1 << zshift;
	lo2 = qinv2 << zshift;
	lo3 = qinv3 << zshift;

#ifdef MUL_LOHI64_SUBROUTINE
	lo0 = __MULH64(q0,lo0);
	lo1 = __MULH64(q1,lo1);
	lo2 = __MULH64(q2,lo2);
	lo3 = __MULH64(q3,lo3);
#else
	MULH64(q0,lo0,lo0);
	MULH64(q1,lo1,lo1);
	MULH64(q2,lo2,lo2);
	MULH64(q3,lo3,lo3);
#endif
	/* hi = 0 in this instance, which simplifies things. */
	y0 = q0 - lo0;
	y1 = q1 - lo1;
	y2 = q2 - lo2;
	y3 = q3 - lo3;

	if((pshift >> j) & (uint64)1)
	{
	#ifdef NOBRANCH
		x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
		x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
		x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
		x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
	#else
		x0 = y0 + y0;	if(x0 >= q0 || x0 < y0) x0 -= q0;
		x1 = y1 + y1;	if(x1 >= q1 || x1 < y1) x1 -= q1;
		x2 = y2 + y2;	if(x2 >= q2 || x2 < y2) x2 -= q2;
		x3 = y3 + y3;	if(x3 >= q3 || x3 < y3) x3 -= q3;
	#endif
	}
	else
	{
		x0 = y0;
		x1 = y1;
		x2 = y2;
		x3 = y3;
	}

#if FAC_DEBUG
  #if 0	/* These appear to be benign: */
	if(x0 >= q0){ sprintf(char_buf, "twopmodq64_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint64_base10_char(str0, x0)], &str1[convert_uint64_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x1 >= q1){ sprintf(char_buf, "twopmodq64_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint64_base10_char(str0, x1)], &str1[convert_uint64_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x2 >= q2){ sprintf(char_buf, "twopmodq64_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint64_base10_char(str0, x2)], &str1[convert_uint64_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x3 >= q3){ sprintf(char_buf, "twopmodq64_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint64_base10_char(str0, x3)], &str1[convert_uint64_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
  #endif
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		MONT_SQR64_q4(x0,x1,x2,x3,q0,q1,q2,q3,qinv0,qinv1,qinv2,qinv3,y0,y1,y2,y3);
		if((pshift >> j) & (uint64)1)
		{
		#ifdef NOBRANCH
			x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
			x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
			x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
			x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
		#else
			x0 = y0 + y0;	if(x0 >= q0 || x0 < y0) x0 -= q0;
			x1 = y1 + y1;	if(x1 >= q1 || x1 < y1) x1 -= q1;
			x2 = y2 + y2;	if(x2 >= q2 || x2 < y2) x2 -= q2;
			x3 = y3 + y3;	if(x3 >= q3 || x3 < y3) x3 -= q3;
		#endif
		}
		else
		{
			x0 = y0;
			x1 = y1;
			x2 = y2;
			x3 = y3;
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
	*checksum2 += x0 + x1 + x2 + x3;

	r = 0;
	if(x0+x0-q0 == 1) r +=  1;
	if(x1+x1-q1 == 1) r +=  2;
	if(x2+x2-q2 == 1) r +=  4;
	if(x3+x3-q3 == 1) r +=  8;
	return r;
}


/******************************/
/*** 8-trial-factor version ***/
/******************************/

uint64 twopmodq64_q8(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
	 int32 j;
	uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3, q4 = 1+(p<<1)*k4, q5 = 1+(p<<1)*k5, q6 = 1+(p<<1)*k6, q7 = 1+(p<<1)*k7
		, qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, y0, y1, y2, y3, y4, y5, y6, y7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
		, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7, r;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
		start_index = 64-leadz64(pshift)-6;

		zshift = 63-ibits64(pshift,start_index,6);	/* Since 6 bits with leftmost bit = 1 is guaranteed
								to be in [32,63], the shift count here is in [0, 31].
								That means that zstart < 2^32. Together with the fact that
								squaring a power of two gives another power of two, we can
								simplify the modmul code sequence for the first iteration.
								Every little bit counts (literally in this case :), right? */
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}
	*checksum1 += q0 + q1 + q2 + q3 + q4 + q5 + q6 + q7;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 && q4 & 1 && q5 & 1 && q6 & 1 && q7 & 1 , "even modulus!");

	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;
	qinv4 = (q4+q4+q4) ^ (uint64)2;
	qinv5 = (q5+q5+q5) ^ (uint64)2;
	qinv6 = (q6+q6+q6) ^ (uint64)2;
	qinv7 = (q7+q7+q7) ^ (uint64)2;

	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
		qinv4 = qinv4*((uint64)2 - q4*qinv4);
		qinv5 = qinv5*((uint64)2 - q5*qinv5);
		qinv6 = qinv6*((uint64)2 - q6*qinv6);
		qinv7 = qinv7*((uint64)2 - q7*qinv7);
	}

	/* Since zstart is a power of two < 2^128, use a
	streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL64(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	lo0 = qinv0 << zshift;
	lo1 = qinv1 << zshift;
	lo2 = qinv2 << zshift;
	lo3 = qinv3 << zshift;
	lo4 = qinv4 << zshift;
	lo5 = qinv5 << zshift;
	lo6 = qinv6 << zshift;
	lo7 = qinv7 << zshift;

#ifdef MUL_LOHI64_SUBROUTINE
	lo0 = __MULH64(q0,lo0);
	lo1 = __MULH64(q1,lo1);
	lo2 = __MULH64(q2,lo2);
	lo3 = __MULH64(q3,lo3);
	lo4 = __MULH64(q4,lo4);
	lo5 = __MULH64(q5,lo5);
	lo6 = __MULH64(q6,lo6);
	lo7 = __MULH64(q7,lo7);
#else
	MULH64(q0,lo0,lo0);
	MULH64(q1,lo1,lo1);
	MULH64(q2,lo2,lo2);
	MULH64(q3,lo3,lo3);
	MULH64(q4,lo4,lo4);
	MULH64(q5,lo5,lo5);
	MULH64(q6,lo6,lo6);
	MULH64(q7,lo7,lo7);
#endif
	/* hi = 0 in this instance, which simplifies things. */
	y0 = q0 - lo0;
	y1 = q1 - lo1;
	y2 = q2 - lo2;
	y3 = q3 - lo3;
	y4 = q4 - lo4;
	y5 = q5 - lo5;
	y6 = q6 - lo6;
	y7 = q7 - lo7;

	if((pshift >> j) & (uint64)1)
	{
	#ifdef NOBRANCH
		x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
		x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
		x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
		x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
		x4 = y4 + y4;	x4 -= q4 & -(x4 >= q4 || x4 < y4);
		x5 = y5 + y5;	x5 -= q5 & -(x5 >= q5 || x5 < y5);
		x6 = y6 + y6;	x6 -= q6 & -(x6 >= q6 || x6 < y6);
		x7 = y7 + y7;	x7 -= q7 & -(x7 >= q7 || x7 < y7);
	#else
		x0 = y0 + y0;	if(x0 >= q0 || x0 < y0) x0 -= q0;
		x1 = y1 + y1;	if(x1 >= q1 || x1 < y1) x1 -= q1;
		x2 = y2 + y2;	if(x2 >= q2 || x2 < y2) x2 -= q2;
		x3 = y3 + y3;	if(x3 >= q3 || x3 < y3) x3 -= q3;
		x4 = y4 + y4;	if(x4 >= q4 || x4 < y4) x4 -= q4;
		x5 = y5 + y5;	if(x5 >= q5 || x5 < y5) x5 -= q5;
		x6 = y6 + y6;	if(x6 >= q6 || x6 < y6) x6 -= q6;
		x7 = y7 + y7;	if(x7 >= q7 || x7 < y7) x7 -= q7;
	#endif
	}
	else
	{
		x0 = y0;
		x1 = y1;
		x2 = y2;
		x3 = y3;
		x4 = y4;
		x5 = y5;
		x6 = y6;
		x7 = y7;
	}

#if FAC_DEBUG
  #if 0	/* These appear to be benign: */
	if(x0 >= q0){ sprintf(char_buf, "twopmodq64_q8 : (x0 = %s) >= (q0 = %s)", &str0[convert_uint64_base10_char(str0, x0)], &str1[convert_uint64_base10_char(str1, q0)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x1 >= q1){ sprintf(char_buf, "twopmodq64_q8 : (x1 = %s) >= (q1 = %s)", &str0[convert_uint64_base10_char(str0, x1)], &str1[convert_uint64_base10_char(str1, q1)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x2 >= q2){ sprintf(char_buf, "twopmodq64_q8 : (x2 = %s) >= (q2 = %s)", &str0[convert_uint64_base10_char(str0, x2)], &str1[convert_uint64_base10_char(str1, q2)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x3 >= q3){ sprintf(char_buf, "twopmodq64_q8 : (x3 = %s) >= (q3 = %s)", &str0[convert_uint64_base10_char(str0, x3)], &str1[convert_uint64_base10_char(str1, q3)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x4 >= q4){ sprintf(char_buf, "twopmodq64_q8 : (x4 = %s) >= (q4 = %s)", &str0[convert_uint64_base10_char(str0, x4)], &str1[convert_uint64_base10_char(str1, q4)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x5 >= q5){ sprintf(char_buf, "twopmodq64_q8 : (x5 = %s) >= (q5 = %s)", &str0[convert_uint64_base10_char(str0, x5)], &str1[convert_uint64_base10_char(str1, q5)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x6 >= q6){ sprintf(char_buf, "twopmodq64_q8 : (x6 = %s) >= (q6 = %s)", &str0[convert_uint64_base10_char(str0, x6)], &str1[convert_uint64_base10_char(str1, q6)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
	if(x7 >= q7){ sprintf(char_buf, "twopmodq64_q8 : (x7 = %s) >= (q7 = %s)", &str0[convert_uint64_base10_char(str0, x7)], &str1[convert_uint64_base10_char(str1, q7)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
  #endif
#endif

	for(j = start_index-2; j >= 0; j--)
	{
	#ifdef MUL_LOHI64_SUBROUTINE
		SQR_LOHI64(x0,&lo0,&hi0);
		SQR_LOHI64(x1,&lo1,&hi1);
		SQR_LOHI64(x2,&lo2,&hi2);
		SQR_LOHI64(x3,&lo3,&hi3);
		SQR_LOHI64(x4,&lo4,&hi4);
		SQR_LOHI64(x5,&lo5,&hi5);
		SQR_LOHI64(x6,&lo6,&hi6);
		SQR_LOHI64(x7,&lo7,&hi7);
	#else
		SQR_LOHI64(x0,lo0,hi0);
		SQR_LOHI64(x1,lo1,hi1);
		SQR_LOHI64(x2,lo2,hi2);
		SQR_LOHI64(x3,lo3,hi3);
		SQR_LOHI64(x4,lo4,hi4);
		SQR_LOHI64(x5,lo5,hi5);
		SQR_LOHI64(x6,lo6,hi6);
		SQR_LOHI64(x7,lo7,hi7);
	#endif
	/*...x^2 mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/

		MULL64(lo0,qinv0,lo0);
		MULL64(lo1,qinv1,lo1);
		MULL64(lo2,qinv2,lo2);
		MULL64(lo3,qinv3,lo3);
		MULL64(lo4,qinv4,lo4);
		MULL64(lo5,qinv5,lo5);
		MULL64(lo6,qinv6,lo6);
		MULL64(lo7,qinv7,lo7);

	#ifdef MUL_LOHI64_SUBROUTINE
		lo0 = __MULH64(q0,lo0);
		lo1 = __MULH64(q1,lo1);
		lo2 = __MULH64(q2,lo2);
		lo3 = __MULH64(q3,lo3);
		lo4 = __MULH64(q4,lo4);
		lo5 = __MULH64(q5,lo5);
		lo6 = __MULH64(q6,lo6);
		lo7 = __MULH64(q7,lo7);
	#else
		MULH64(q0,lo0,lo0);
		MULH64(q1,lo1,lo1);
		MULH64(q2,lo2,lo2);
		MULH64(q3,lo3,lo3);
		MULH64(q4,lo4,lo4);
		MULH64(q5,lo5,lo5);
		MULH64(q6,lo6,lo6);
		MULH64(q7,lo7,lo7);
	#endif
		y0 = hi0 - lo0;
		y1 = hi1 - lo1;
		y2 = hi2 - lo2;
		y3 = hi3 - lo3;
		y4 = hi4 - lo4;
		y5 = hi5 - lo5;
		y6 = hi6 - lo6;
		y7 = hi7 - lo7;

	#ifdef NOBRANCH
		y0 += q0 & -(y0 > hi0);
		y1 += q1 & -(y1 > hi1);
		y2 += q2 & -(y2 > hi2);
		y3 += q3 & -(y3 > hi3);
		y4 += q4 & -(y4 > hi4);
		y5 += q5 & -(y5 > hi5);
		y6 += q6 & -(y6 > hi6);
		y7 += q7 & -(y7 > hi7);
	#else
		if(y0 > hi0) y0 += q0;
		if(y1 > hi1) y1 += q1;
		if(y2 > hi2) y2 += q2;
		if(y3 > hi3) y3 += q3;
		if(y4 > hi4) y4 += q4;
		if(y5 > hi5) y5 += q5;
		if(y6 > hi6) y6 += q6;
		if(y7 > hi7) y7 += q7;
	#endif

		if((pshift >> j) & (uint64)1)
		{
		#ifdef NOBRANCH
			x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
			x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
			x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
			x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
			x4 = y4 + y4;	x4 -= q4 & -(x4 >= q4 || x4 < y4);
			x5 = y5 + y5;	x5 -= q5 & -(x5 >= q5 || x5 < y5);
			x6 = y6 + y6;	x6 -= q6 & -(x6 >= q6 || x6 < y6);
			x7 = y7 + y7;	x7 -= q7 & -(x7 >= q7 || x7 < y7);
		#else
			x0 = y0 + y0;	if(x0 >= q0 || x0 < y0) x0 -= q0;
			x1 = y1 + y1;	if(x1 >= q1 || x1 < y1) x1 -= q1;
			x2 = y2 + y2;	if(x2 >= q2 || x2 < y2) x2 -= q2;
			x3 = y3 + y3;	if(x3 >= q3 || x3 < y3) x3 -= q3;
			x4 = y4 + y4;	if(x4 >= q4 || x4 < y4) x4 -= q4;
			x5 = y5 + y5;	if(x5 >= q5 || x5 < y5) x5 -= q5;
			x6 = y6 + y6;	if(x6 >= q6 || x6 < y6) x6 -= q6;
			x7 = y7 + y7;	if(x7 >= q7 || x7 < y7) x7 -= q7;
		#endif
		}
		else
		{
			x0 = y0;
			x1 = y1;
			x2 = y2;
			x3 = y3;
			x4 = y4;
			x5 = y5;
			x6 = y6;
			x7 = y7;
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	*checksum2 += x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7;

	r = 0;
	if(x0+x0-q0 == 1) r +=  1;
	if(x1+x1-q1 == 1) r +=  2;
	if(x2+x2-q2 == 1) r +=  4;
	if(x3+x3-q3 == 1) r +=  8;
	if(x4+x4-q4 == 1) r += 16;
	if(x5+x5-q5 == 1) r += 32;
	if(x6+x6-q6 == 1) r += 64;
	if(x7+x7-q7 == 1) r +=128;
	return r;
}

#endif /* endifndef(YES_ASM) */


/***********************************************************************************/
/*** 65-BIT INPUTS: it is assumed here that q has a hidden 65th bit! ***************/
/***********************************************************************************/

uint64 twopmodq65(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k)
{
#if FAC_DEBUG
	int dbg = (p == 0);
#endif
	 int32 j;
	uint64 q, qinv, x, y, lo, hi, r;
	uint64 A, B, t;
	static uint64 psave = 0, pshift, zstart;
	static uint32 start_index, first_entry = TRUE;

#if FAC_DEBUG
if(dbg)printf("twopmodq65:\n");
#endif
	// Assume q is 65-bits here, so check that during construction of q = 2.k.p+1:
	q = k*p;	ASSERT(HERE, q+q < q, "q not 65 bits!");
	q = (q << 1) + 1;
	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
		start_index = 64-leadz64(pshift)-6;

		zstart = ((uint64)1) << (63-ibits64(pshift,start_index,6));
		pshift = ~pshift;
	}
	*checksum1 += q;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q & 0x1, "q must be odd!");
	qinv = (q+q+q) ^ (uint64)2;
	for(j = 0; j < 4; j++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

	/*...Initialize the binary powering...*/
	x = zstart;
	A = 0;

#if FAC_DEBUG
	if(dbg) printf("twopmodq65: jhi = %d\n", start_index-1);
#endif
	for(j = start_index-1; j >= 0; j--)
	{
#ifdef MUL_LOHI64_SUBROUTINE
		SQR_LOHI64(x,&lo,&hi);
#else
		SQR_LOHI64(x,lo,hi);
#endif

		if(A > 0)
		{
			y  = (-A & x);
			hi+= y; if(hi<y) ++A;
			hi+= y; if(hi<y) ++A;
		}

	/*...x^2 mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/

		lo *= qinv;

#ifdef MUL_LOHI64_SUBROUTINE
		y = __MULH64(q,lo);
#else
		MULH64(q,lo,y);	/* Need original lo for a few more steps here, so store MULH result in y. */
#endif
		B = 0;
		y += lo; if(y < lo) ++B;	/* Upper (65th) bit of q guaranteed to be 1. */
	/*
	!	SQR_HIGH_EXACT(x) now in hi + A*2^64,  MULH_EXACT(q, lo*qinv) in y + B*2^64.
	!	Get SQR_HIGH_EXACT(x) - MULH_EXACT(q, lo*qinv), normalizing as we go.
	*/
#if FAC_DEBUG
	if(dbg) printf("twopmodq65: while(%llu++ < %llu) || (%llu+=%llu < %llu)\n",A,B,hi,q,y);
#endif
		while(A < B || (A == B && hi < y))	/* SQR_HIGH_EXACT(x) < MULH_EXACT(q, lo*qinv); add q until >= . */
		{
			++A; hi += q;
			if(hi < q) ++A;	/* had a carry */
		}

	/* Do the subtraction. Result is in (hi,A). */

	#if FAC_DEBUG
		/*ASSERT(HERE, A > B || (A == B && hi >= y), "twopmodq65 : A > B || (A == B && hi >= y)"); */
	#endif

		A -= B; x = hi; hi -= y;
		if(hi > x) --A;	/* had a borrow */

	/* ...and normalize. Result is in (x, A).. */

		x = hi;
#if FAC_DEBUG
	if(dbg) printf("twopmodq65: while(A=%llu-- > 1) || (%llu-=%llu >=%llu)\n",A,x,q,q);
#endif
		while(A > 1 || (A == 1 && x >= q))
		{
			--A; x -= q;

#if FAC_DEBUG
	if(dbg) printf("twopmodq65: A = %llu, x-q = %llu, q = %llu, hi = %llu\n",A,x,q,hi);
	if(dbg) printf("twopmodq65: (x > hi) = %llu\n",(x > hi));
	if(dbg) printf("twopmodq65: (x <=hi) = %llu\n",(x <=hi));
	if(dbg) printf("twopmodq65: (hi < x) = %llu\n",(hi < x));
	if(dbg) printf("twopmodq65: (hi <=x) = %llu\n",(hi <=x));
	if(dbg) printf("twopmodq65: (hi >=x) = %llu\n",(hi >=x));
	if(dbg) printf("twopmodq65: (hi -x ) = %llu - %llu = %llu\n",hi,x,(hi - x));
	if(dbg) printf("twopmodq65: (hi.-x.) = %lf  - %lf  = %lf \n",(double)hi,(double)x,((double)hi - (double)x));
#endif
			/* had a borrow: */
			A -= (x > hi);	/* 12/16/2003: The HP C compiler for IA64 under HPUX
								f***s up an (x > hi) compare for certain operands,
								so replace that with (hi <= x), which appears to get
								handled OK. */
#if FAC_DEBUG
	ASSERT(HERE, ((double)x > (double)hi) == (x > hi),"((double)x > (double)hi) == (x > hi)");
	ASSERT(HERE, (int64)A >=0,"(int64)A >=0");
#endif
		}

		if((pshift >> j) & (uint64)1)
		{
			/* Add x to itself... */
			y = x + x; B = A + A;
			if(y < x) ++B;		/* had a carry */

			/* ...and normalize the result. */
#if FAC_DEBUG
	if(dbg) printf("twopmodq65: while(B=%llu-- > 1) || (%llu-=%llu >=%llu)\n",B,y,q,q);
#endif
			while(B > 1 || (B == 1 && y >= q))
			{
				--B; t = y; y -= q;
				if(y > t) --B;		/* had a borrow */
			}
			x = y; A = B;
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	/* Add x to itself... */
	y = x + x; B = A + A;
	if(y < x) ++B;		/* had a carry */
	/* ...and normalize the result. */

#if FAC_DEBUG
	if(dbg) printf("twopmodq65:#while(%llu-- > 1) || (%llu-=%llu >=%llu)\n",B,y,q,q);
#endif
	while(B > 1 || (B == 1 && y >= q))
	{
		--B; t = y; y -= q;
		if(y > t) --B;		/* had a borrow */
	}

	/*if(y == (uint64)1) ASSERT(HERE, B == 0, "twopmodq65 : B == 0");*/
	*checksum2 += x;

	r = 0;
	if(y == 1)
	{
		r += 1;
	}
	return r;
}

/*** 4-trial-factor version ***/
uint64 twopmodq65_q4(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
	 int32 j;
	uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3
		, qinv0, qinv1, qinv2, qinv3
		, x0, x1, x2, x3
		, y0, y1, y2, y3
		, A0, A1, A2, A3
		, B0, B1, B2, B3
		, lo0, lo1, lo2, lo3
		, hi0, hi1, hi2, hi3, r, t;
	static uint64 psave = 0, pshift, zstart;
	static uint32 start_index, first_entry = TRUE;
#ifdef MOD_INI_Q4
	uint32 ql0,ql1,ql2,ql3,qh0,qh1,qh2,qh3,qil0,qil1,qil2,qil3,qih0,qih1,qih2,qih3;
#endif

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
		start_index = 64-leadz64(pshift)-6;

		zstart = ((uint64)1) << (63-ibits64(pshift,start_index,6));
		pshift = ~pshift;
	}
	*checksum1 += q0 + q1 + q2 + q3;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 , "even modulus!");

	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;

	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
	}

#ifdef MOD_INI_Q4
	MOD_INI_Q4(
	 q0,qinv0
	,q1,qinv1
	,q2,qinv2
	,q3,qinv3
	);
#endif

/*...Initialize the binary powering...*/

	x0 = x1 = x2 = x3 = zstart;
	A0 = A1 = A2 = A3 = 0;

	for(j = start_index-1; j >= 0; j--)
	{

#ifdef MOD_INI_Q4
		MOD_SQR_Q4(
		 x0,hi0,y0
		,x1,hi1,y1
		,x2,hi2,y2
		,x3,hi3,y3
		);
#else
	#ifdef MUL_LOHI64_SUBROUTINE
		SQR_LOHI64(x0,&lo0,&hi0);
		SQR_LOHI64(x1,&lo1,&hi1);
		SQR_LOHI64(x2,&lo2,&hi2);
		SQR_LOHI64(x3,&lo3,&hi3);
	#else
		SQR_LOHI64(x0,lo0,hi0);
		SQR_LOHI64(x1,lo1,hi1);
		SQR_LOHI64(x2,lo2,hi2);
		SQR_LOHI64(x3,lo3,hi3);
	#endif

	/*...x^2 mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/

		MULL64(lo0,qinv0,lo0);
		MULL64(lo1,qinv1,lo1);
		MULL64(lo2,qinv2,lo2);
		MULL64(lo3,qinv3,lo3);

	#ifdef MUL_LOHI64_SUBROUTINE
		y0 = __MULH64(q0,lo0);
		y1 = __MULH64(q1,lo1);
		y2 = __MULH64(q2,lo2);
		y3 = __MULH64(q3,lo3);
	#else
		MULH64(q0,lo0,y0);
		MULH64(q1,lo1,y1);
		MULH64(q2,lo2,y2);
		MULH64(q3,lo3,y3);
	#endif
#endif

		/* Use B0-3 as temporaries here... */
#ifdef NOBRANCH
		if(A0 > 0){ B0 = (-A0 & x0); hi0 += B0; A0 += (hi0 < B0); hi0 += B0; A0 += (hi0 < B0); }
		if(A1 > 0){ B1 = (-A1 & x1); hi1 += B1; A1 += (hi1 < B1); hi1 += B1; A1 += (hi1 < B1); }
		if(A2 > 0){ B2 = (-A2 & x2); hi2 += B2; A2 += (hi2 < B2); hi2 += B2; A2 += (hi2 < B2); }
		if(A3 > 0){ B3 = (-A3 & x3); hi3 += B3; A3 += (hi3 < B3); hi3 += B3; A3 += (hi3 < B3); }
#else
		if(A0 > 0){ B0 = (-A0 & x0); hi0 += B0; if(hi0 < B0) ++A0; hi0 += B0; if(hi0 < B0) ++A0; }
		if(A1 > 0){ B1 = (-A1 & x1); hi1 += B1; if(hi1 < B1) ++A1; hi1 += B1; if(hi1 < B1) ++A1; }
		if(A2 > 0){ B2 = (-A2 & x2); hi2 += B2; if(hi2 < B2) ++A2; hi2 += B2; if(hi2 < B2) ++A2; }
		if(A3 > 0){ B3 = (-A3 & x3); hi3 += B3; if(hi3 < B3) ++A3; hi3 += B3; if(hi3 < B3) ++A3; }
#endif
		/* ...And then give them their proper value here: */
		y0 += lo0; B0 = (y0 < lo0);
		y1 += lo1; B1 = (y1 < lo1);
		y2 += lo2; B2 = (y2 < lo2);
		y3 += lo3; B3 = (y3 < lo3);
	/*
	!	SQR_HIGH_EXACT(x) now in hi + A*2^64,  MULH_EXACT(q, lo*qinv) in y + B*2^64. Get SQR_HIGH_EXACT(x) - MULH_EXACT(q, lo*qinv), normalizing
	!	as we go. Since MULH_EXACT(q, lo*qinv) < q, we need to add at most one q to SQR_HIGH_EXACT(x) to make the difference nonnegative.
	!	On the other side of the ledger, since SQR_HIGH_EXACT(x) is < 2q, the difference will be < 2q or < 3q, respectively,
	!	depending on whether we add the q to SQR_HIGH_EXACT(x) conditionally or unconditionally.
	*/
		if(A0 < B0 || (A0 == B0 && hi0 < y0)){ ++A0; hi0 += q0; if(hi0 < q0) ++A0; }
		if(A1 < B1 || (A1 == B1 && hi1 < y1)){ ++A1; hi1 += q1; if(hi1 < q1) ++A1; }
		if(A2 < B2 || (A2 == B2 && hi2 < y2)){ ++A2; hi2 += q2; if(hi2 < q2) ++A2; }
		if(A3 < B3 || (A3 == B3 && hi3 < y3)){ ++A3; hi3 += q3; if(hi3 < q3) ++A3; }

	/* Do the subtraction. Result is in (hi,A). */

		A0 -= B0; x0 = hi0; hi0 -= y0; if(hi0 > x0) --A0;
		A1 -= B1; x1 = hi1; hi1 -= y1; if(hi1 > x1) --A1;
		A2 -= B2; x2 = hi2; hi2 -= y2; if(hi2 > x2) --A2;
		A3 -= B3; x3 = hi3; hi3 -= y3; if(hi3 > x3) --A3;

	/* ...and normalize. Result is in (x, A), and is < 2q if the previous normalization was done conditionally, < 3q otherwise. */

		x0 = hi0; if(A0 > 1 || (A0 == 1 && x0 >= q0)){ --A0; x0 -= q0; if(x0 > hi0) --A0; }
		x1 = hi1; if(A1 > 1 || (A1 == 1 && x1 >= q1)){ --A1; x1 -= q1; if(x1 > hi1) --A1; }
		x2 = hi2; if(A2 > 1 || (A2 == 1 && x2 >= q2)){ --A2; x2 -= q2; if(x2 > hi2) --A2; }
		x3 = hi3; if(A3 > 1 || (A3 == 1 && x3 >= q3)){ --A3; x3 -= q3; if(x3 > hi3) --A3; }

	/* Add x to itself and normalize the result. */

		if((pshift >> j) & (uint64)1)
		{
			y0 = x0 + x0; B0 = A0 + A0; if(y0 < x0) ++B0; if(B0 > 1 || (B0 == 1 && y0 >= q0)){ --B0; t = y0; y0 -= q0; if(y0 > t) --B0; }
			y1 = x1 + x1; B1 = A1 + A1; if(y1 < x1) ++B1; if(B1 > 1 || (B1 == 1 && y1 >= q1)){ --B1; t = y1; y1 -= q1; if(y1 > t) --B1; }
			y2 = x2 + x2; B2 = A2 + A2; if(y2 < x2) ++B2; if(B2 > 1 || (B2 == 1 && y2 >= q2)){ --B2; t = y2; y2 -= q2; if(y2 > t) --B2; }
			y3 = x3 + x3; B3 = A3 + A3; if(y3 < x3) ++B3; if(B3 > 1 || (B3 == 1 && y3 >= q3)){ --B3; t = y3; y3 -= q3; if(y3 > t) --B3; }

			x0 = y0; A0 = B0;
			x1 = y1; A1 = B1;
			x2 = y2; A2 = B2;
			x3 = y3; A3 = B3;
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	y0 = x0 + x0; B0 = A0 + A0; if(y0 < x0) ++B0; if(B0 > 1 || (B0 == 1 && y0 >= q0)){ --B0; t = y0; y0 -= q0; if(y0 > t) --B0; }
	y1 = x1 + x1; B1 = A1 + A1; if(y1 < x1) ++B1; if(B1 > 1 || (B1 == 1 && y1 >= q1)){ --B1; t = y1; y1 -= q1; if(y1 > t) --B1; }
	y2 = x2 + x2; B2 = A2 + A2; if(y2 < x2) ++B2; if(B2 > 1 || (B2 == 1 && y2 >= q2)){ --B2; t = y2; y2 -= q2; if(y2 > t) --B2; }
	y3 = x3 + x3; B3 = A3 + A3; if(y3 < x3) ++B3; if(B3 > 1 || (B3 == 1 && y3 >= q3)){ --B3; t = y3; y3 -= q3; if(y3 > t) --B3; }

	*checksum2 += x0 + x1 + x2 + x3;

	r = 0;
	if(y0 == 1) r +=  1;
	if(y1 == 1) r +=  2;
	if(y2 == 1) r +=  4;
	if(y3 == 1) r +=  8;
	return r;
}

/*** 8-trial-factor version ***/
uint64 twopmodq65_q8(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
	 int32 j;
	uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3, q4 = 1+(p<<1)*k4, q5 = 1+(p<<1)*k5, q6 = 1+(p<<1)*k6, q7 = 1+(p<<1)*k7
		, qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, y0, y1, y2, y3, y4, y5, y6, y7
		, A0, A1, A2, A3, A4, A5, A6, A7
		, B0, B1, B2, B3, B4, B5, B6, B7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7
		, hi0, hi1, hi2, hi3, hi4, hi5, hi6, hi7, r, t;
	static uint64 psave = 0, pshift, zstart;
	static uint32 start_index, first_entry = TRUE;

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
		start_index = 64-leadz64(pshift)-6;

		zstart = ((uint64)1) << (63-ibits64(pshift,start_index,6));
		pshift = ~pshift;
	}
	*checksum1 += q0 + q1 + q2 + q3 + q4 + q5 + q6 + q7;

	/* q must be odd for Montgomery-style modmul to work: */
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 && q4 & 1 && q5 & 1 && q6 & 1 && q7 & 1 , "even modulus!");

	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;
	qinv4 = (q4+q4+q4) ^ (uint64)2;
	qinv5 = (q5+q5+q5) ^ (uint64)2;
	qinv6 = (q6+q6+q6) ^ (uint64)2;
	qinv7 = (q7+q7+q7) ^ (uint64)2;

	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
		qinv4 = qinv4*((uint64)2 - q4*qinv4);
		qinv5 = qinv5*((uint64)2 - q5*qinv5);
		qinv6 = qinv6*((uint64)2 - q6*qinv6);
		qinv7 = qinv7*((uint64)2 - q7*qinv7);
	}

	/*...Initialize the binary powering...*/

	x0 = x1 = x2 = x3 = x4 = x5 = x6 = x7 = zstart;
	A0 = A1 = A2 = A3 = A4 = A5 = A6 = A7 = 0;

	for(j = start_index-1; j >= 0; j--)
	{
#ifdef MUL_LOHI64_SUBROUTINE
		SQR_LOHI64(x0,&lo0,&hi0);
		SQR_LOHI64(x1,&lo1,&hi1);
		SQR_LOHI64(x2,&lo2,&hi2);
		SQR_LOHI64(x3,&lo3,&hi3);
		SQR_LOHI64(x4,&lo4,&hi4);
		SQR_LOHI64(x5,&lo5,&hi5);
		SQR_LOHI64(x6,&lo6,&hi6);
		SQR_LOHI64(x7,&lo7,&hi7);
#else
		SQR_LOHI64(x0,lo0,hi0);
		SQR_LOHI64(x1,lo1,hi1);
		SQR_LOHI64(x2,lo2,hi2);
		SQR_LOHI64(x3,lo3,hi3);
		SQR_LOHI64(x4,lo4,hi4);
		SQR_LOHI64(x5,lo5,hi5);
		SQR_LOHI64(x6,lo6,hi6);
		SQR_LOHI64(x7,lo7,hi7);
#endif

	/*...x^2 mod q is returned in x. On MIPS, we discard the lower half of DMULTU(q,x*y*qinv).	*/

		MULL64(lo0,qinv0,lo0);
		MULL64(lo1,qinv1,lo1);
		MULL64(lo2,qinv2,lo2);
		MULL64(lo3,qinv3,lo3);
		MULL64(lo4,qinv4,lo4);
		MULL64(lo5,qinv5,lo5);
		MULL64(lo6,qinv6,lo6);
		MULL64(lo7,qinv7,lo7);

	#ifdef MUL_LOHI64_SUBROUTINE
		y0 = __MULH64(q0,lo0);
		y1 = __MULH64(q1,lo1);
		y2 = __MULH64(q2,lo2);
		y3 = __MULH64(q3,lo3);
		y4 = __MULH64(q4,lo4);
		y5 = __MULH64(q5,lo5);
		y6 = __MULH64(q6,lo6);
		y7 = __MULH64(q7,lo7);
	#else
		MULH64(q0,lo0,y0);
		MULH64(q1,lo1,y1);
		MULH64(q2,lo2,y2);
		MULH64(q3,lo3,y3);
		MULH64(q4,lo4,y4);
		MULH64(q5,lo5,y5);
		MULH64(q6,lo6,y6);
		MULH64(q7,lo7,y7);
	#endif

		/* Use B0-3 as temporaries here... */
#ifdef NOBRANCH
		if(A0 > 0){ B0 = (-A0 & x0); hi0 += B0; A0 += (hi0 < B0); hi0 += B0; A0 += (hi0 < B0); }
		if(A1 > 0){ B1 = (-A1 & x1); hi1 += B1; A1 += (hi1 < B1); hi1 += B1; A1 += (hi1 < B1); }
		if(A2 > 0){ B2 = (-A2 & x2); hi2 += B2; A2 += (hi2 < B2); hi2 += B2; A2 += (hi2 < B2); }
		if(A3 > 0){ B3 = (-A3 & x3); hi3 += B3; A3 += (hi3 < B3); hi3 += B3; A3 += (hi3 < B3); }
		if(A4 > 0){ B4 = (-A4 & x4); hi4 += B4; A4 += (hi4 < B4); hi4 += B4; A4 += (hi4 < B4); }
		if(A5 > 0){ B5 = (-A5 & x5); hi5 += B5; A5 += (hi5 < B5); hi5 += B5; A5 += (hi5 < B5); }
		if(A6 > 0){ B6 = (-A6 & x6); hi6 += B6; A6 += (hi6 < B6); hi6 += B6; A6 += (hi6 < B6); }
		if(A7 > 0){ B7 = (-A7 & x7); hi7 += B7; A7 += (hi7 < B7); hi7 += B7; A7 += (hi7 < B7); }
#else
		if(A0 > 0){ B0 = (-A0 & x0); hi0 += B0; if(hi0 < B0) ++A0; hi0 += B0; if(hi0 < B0) ++A0; }
		if(A1 > 0){ B1 = (-A1 & x1); hi1 += B1; if(hi1 < B1) ++A1; hi1 += B1; if(hi1 < B1) ++A1; }
		if(A2 > 0){ B2 = (-A2 & x2); hi2 += B2; if(hi2 < B2) ++A2; hi2 += B2; if(hi2 < B2) ++A2; }
		if(A3 > 0){ B3 = (-A3 & x3); hi3 += B3; if(hi3 < B3) ++A3; hi3 += B3; if(hi3 < B3) ++A3; }
		if(A4 > 0){ B4 = (-A4 & x4); hi4 += B4; if(hi4 < B4) ++A4; hi4 += B4; if(hi4 < B4) ++A4; }
		if(A5 > 0){ B5 = (-A5 & x5); hi5 += B5; if(hi5 < B5) ++A5; hi5 += B5; if(hi5 < B5) ++A5; }
		if(A6 > 0){ B6 = (-A6 & x6); hi6 += B6; if(hi6 < B6) ++A6; hi6 += B6; if(hi6 < B6) ++A6; }
		if(A7 > 0){ B7 = (-A7 & x7); hi7 += B7; if(hi7 < B7) ++A7; hi7 += B7; if(hi7 < B7) ++A7; }
#endif
		/* ...And then give them their proper value here: */
		y0 += lo0; B0 = (y0 < lo0);
		y1 += lo1; B1 = (y1 < lo1);
		y2 += lo2; B2 = (y2 < lo2);
		y3 += lo3; B3 = (y3 < lo3);
		y4 += lo4; B4 = (y4 < lo4);
		y5 += lo5; B5 = (y5 < lo5);
		y6 += lo6; B6 = (y6 < lo6);
		y7 += lo7; B7 = (y7 < lo7);
	/*
	!	SQR_HIGH_EXACT(x) now in hi + A*2^64,  MULH_EXACT(q, lo*qinv) in y + B*2^64. Get SQR_HIGH_EXACT(x) - MULH_EXACT(q, lo*qinv), normalizing
	!	as we go. Since MULH_EXACT(q, lo*qinv) < q, we need to add at most one q to SQR_HIGH_EXACT(x) to make the difference nonnegative.
	!	On the other side of the ledger, since SQR_HIGH_EXACT(x) is < 2q, the difference will be < 2q or < 3q, respectively,
	!	depending on whether we add the q to SQR_HIGH_EXACT(x) conditionally or unconditionally.
	*/
		if(A0 < B0 || (A0 == B0 && hi0 < y0)){ ++A0; hi0 += q0; if(hi0 < q0) ++A0; }
		if(A1 < B1 || (A1 == B1 && hi1 < y1)){ ++A1; hi1 += q1; if(hi1 < q1) ++A1; }
		if(A2 < B2 || (A2 == B2 && hi2 < y2)){ ++A2; hi2 += q2; if(hi2 < q2) ++A2; }
		if(A3 < B3 || (A3 == B3 && hi3 < y3)){ ++A3; hi3 += q3; if(hi3 < q3) ++A3; }
		if(A4 < B4 || (A4 == B4 && hi4 < y4)){ ++A4; hi4 += q4; if(hi4 < q4) ++A4; }
		if(A5 < B5 || (A5 == B5 && hi5 < y5)){ ++A5; hi5 += q5; if(hi5 < q5) ++A5; }
		if(A6 < B6 || (A6 == B6 && hi6 < y6)){ ++A6; hi6 += q6; if(hi6 < q6) ++A6; }
		if(A7 < B7 || (A7 == B7 && hi7 < y7)){ ++A7; hi7 += q7; if(hi7 < q7) ++A7; }

	/* Do the subtraction. Result is in (hi,A). */

		A0 -= B0; x0 = hi0; hi0 -= y0; if(hi0 > x0) --A0;
		A1 -= B1; x1 = hi1; hi1 -= y1; if(hi1 > x1) --A1;
		A2 -= B2; x2 = hi2; hi2 -= y2; if(hi2 > x2) --A2;
		A3 -= B3; x3 = hi3; hi3 -= y3; if(hi3 > x3) --A3;
		A4 -= B4; x4 = hi4; hi4 -= y4; if(hi4 > x4) --A4;
		A5 -= B5; x5 = hi5; hi5 -= y5; if(hi5 > x5) --A5;
		A6 -= B6; x6 = hi6; hi6 -= y6; if(hi6 > x6) --A6;
		A7 -= B7; x7 = hi7; hi7 -= y7; if(hi7 > x7) --A7;

	/* ...and normalize. Result is in (x, A), and is < 2q if the previous normalization was done conditionally, < 3q otherwise. */

		x0 = hi0; if(A0 > 1 || (A0 == 1 && x0 >= q0)){ --A0; x0 -= q0; if(x0 > hi0) --A0; }
		x1 = hi1; if(A1 > 1 || (A1 == 1 && x1 >= q1)){ --A1; x1 -= q1; if(x1 > hi1) --A1; }
		x2 = hi2; if(A2 > 1 || (A2 == 1 && x2 >= q2)){ --A2; x2 -= q2; if(x2 > hi2) --A2; }
		x3 = hi3; if(A3 > 1 || (A3 == 1 && x3 >= q3)){ --A3; x3 -= q3; if(x3 > hi3) --A3; }
		x4 = hi4; if(A4 > 1 || (A4 == 1 && x4 >= q4)){ --A4; x4 -= q4; if(x4 > hi4) --A4; }
		x5 = hi5; if(A5 > 1 || (A5 == 1 && x5 >= q5)){ --A5; x5 -= q5; if(x5 > hi5) --A5; }
		x6 = hi6; if(A6 > 1 || (A6 == 1 && x6 >= q6)){ --A6; x6 -= q6; if(x6 > hi6) --A6; }
		x7 = hi7; if(A7 > 1 || (A7 == 1 && x7 >= q7)){ --A7; x7 -= q7; if(x7 > hi7) --A7; }

	/* Add x to itself and normalize the result. */

		if((pshift >> j) & (uint64)1)
		{
			y0 = x0 + x0; B0 = A0 + A0; if(y0 < x0) ++B0; if(B0 > 1 || (B0 == 1 && y0 >= q0)){ --B0; t = y0; y0 -= q0; if(y0 > t) --B0; }
			y1 = x1 + x1; B1 = A1 + A1; if(y1 < x1) ++B1; if(B1 > 1 || (B1 == 1 && y1 >= q1)){ --B1; t = y1; y1 -= q1; if(y1 > t) --B1; }
			y2 = x2 + x2; B2 = A2 + A2; if(y2 < x2) ++B2; if(B2 > 1 || (B2 == 1 && y2 >= q2)){ --B2; t = y2; y2 -= q2; if(y2 > t) --B2; }
			y3 = x3 + x3; B3 = A3 + A3; if(y3 < x3) ++B3; if(B3 > 1 || (B3 == 1 && y3 >= q3)){ --B3; t = y3; y3 -= q3; if(y3 > t) --B3; }
			y4 = x4 + x4; B4 = A4 + A4; if(y4 < x4) ++B4; if(B4 > 1 || (B4 == 1 && y4 >= q4)){ --B4; t = y4; y4 -= q4; if(y4 > t) --B4; }
			y5 = x5 + x5; B5 = A5 + A5; if(y5 < x5) ++B5; if(B5 > 1 || (B5 == 1 && y5 >= q5)){ --B5; t = y5; y5 -= q5; if(y5 > t) --B5; }
			y6 = x6 + x6; B6 = A6 + A6; if(y6 < x6) ++B6; if(B6 > 1 || (B6 == 1 && y6 >= q6)){ --B6; t = y6; y6 -= q6; if(y6 > t) --B6; }
			y7 = x7 + x7; B7 = A7 + A7; if(y7 < x7) ++B7; if(B7 > 1 || (B7 == 1 && y7 >= q7)){ --B7; t = y7; y7 -= q7; if(y7 > t) --B7; }

			x0 = y0; A0 = B0;
			x1 = y1; A1 = B1;
			x2 = y2; A2 = B2;
			x3 = y3; A3 = B3;
			x4 = y4; A4 = B4;
			x5 = y5; A5 = B5;
			x6 = y6; A6 = B6;
			x7 = y7; A7 = B7;
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	y0 = x0 + x0; B0 = A0 + A0; if(y0 < x0) ++B0; if(B0 > 1 || (B0 == 1 && y0 >= q0)){ --B0; t = y0; y0 -= q0; if(y0 > t) --B0; }
	y1 = x1 + x1; B1 = A1 + A1; if(y1 < x1) ++B1; if(B1 > 1 || (B1 == 1 && y1 >= q1)){ --B1; t = y1; y1 -= q1; if(y1 > t) --B1; }
	y2 = x2 + x2; B2 = A2 + A2; if(y2 < x2) ++B2; if(B2 > 1 || (B2 == 1 && y2 >= q2)){ --B2; t = y2; y2 -= q2; if(y2 > t) --B2; }
	y3 = x3 + x3; B3 = A3 + A3; if(y3 < x3) ++B3; if(B3 > 1 || (B3 == 1 && y3 >= q3)){ --B3; t = y3; y3 -= q3; if(y3 > t) --B3; }
	y4 = x4 + x4; B4 = A4 + A4; if(y4 < x4) ++B4; if(B4 > 1 || (B4 == 1 && y4 >= q4)){ --B4; t = y4; y4 -= q4; if(y4 > t) --B4; }
	y5 = x5 + x5; B5 = A5 + A5; if(y5 < x5) ++B5; if(B5 > 1 || (B5 == 1 && y5 >= q5)){ --B5; t = y5; y5 -= q5; if(y5 > t) --B5; }
	y6 = x6 + x6; B6 = A6 + A6; if(y6 < x6) ++B6; if(B6 > 1 || (B6 == 1 && y6 >= q6)){ --B6; t = y6; y6 -= q6; if(y6 > t) --B6; }
	y7 = x7 + x7; B7 = A7 + A7; if(y7 < x7) ++B7; if(B7 > 1 || (B7 == 1 && y7 >= q7)){ --B7; t = y7; y7 -= q7; if(y7 > t) --B7; }

	*checksum2 += x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7;

	r = 0;
	if(y0 == 1) r +=  1;
	if(y1 == 1) r +=  2;
	if(y2 == 1) r +=  4;
	if(y3 == 1) r +=  8;
	if(y4 == 1) r += 16;
	if(y5 == 1) r += 32;
	if(y6 == 1) r += 64;
	if(y7 == 1) r +=128;
	return r;
}

