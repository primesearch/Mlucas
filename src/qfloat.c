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

#include "qfloat.h"

/* Portable module providing basic 128-bit float emulation, using IEEE64 format for the most-significant
64-bits and simply adding 64-bits to the significand length, i.e. 1 bit for sign, 11 for scaled binary
exponent, and 117 for the significand, which includes 1 hidden bit.

*** ASSUMES: ***
1. The operands are normal - some opeartions may work with denormals, but not all;
2. The endianness is the same for floating-point numbers as for integers.

No support for NaNs of the various kinds is provided. Aagin, some ops may work fine with NaNs, but most won't.
*/
//#undef X87_ASM

#if 0
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC))
	#define X87_ASM
	#if FP_MANTISSA_BITS_DOUBLE != 64
		#error x86_64 asm requires FP_MANTISSA_BITS_DOUBLE == 64!
	#endif
#endif
#endif

// Set nonzero to enable timing report of various qfloat operations:
#define TIMING_TEST	0

/*
16 Sep 2012: First precise timing tests of basic qfloat operations:
Operation	#cycles		Comments													09/25/2012	Notes
=========	=======		============================								========	========================
qfdbl          3.12																	    3.04
qfmul_pow2     0.01																	   24.78	Timing ballooned after added denormal handling
qfmul         88.11		Not terrible, but each of these 3 core arithmetic ops		   74.59
qfadd         65.30		needs to be speeded by at least 1.5-2x; try some x86_64		   57.59
qfsub        101.44		assembly code, using mix of CF-usage and sse 64-bit math.	   74.89	Cut cycles off this via inline-ASM of x86 BSR instruction
qfsqrt      1112.69																	  558.46
qfinv        712.05																	  347.39	Init with 64-bit x87 result, 1 N-R iter
qfdiv        807.13																	  426.03	Init with 64-bit x87 result, 1 N-R iter
qfsin       4117.68																	 3583.69	Replace /(Pi/2) with *(2/Pi)
qfcos       3967.09																	 3662.32	Replace /(Pi/2) with *(2/Pi)
qftan      11386.77		impl simply as sin/cos										 4798.56	Cut cycles via sin*isqrt(1-sin^2)
qfcot			n/a																	 4923.35
qfatan     34836.71		Cut to ~25000 using faster tan(x) ... Try AGM-based?		 5479.69	Init with 64-bit x87 result, 1 N-R iter cuts time; further adding qftan_and_cos() saves another 30%
qflog      14013.66																	 4805.48	Init with 64-bit x87 result, 1 N-R iter
qfexp       7347.76		Cut to ~6400 using nested-product and INC in place of ADD	 4113.49	Use cosh/sinh trick saves ~2000 cycles
qfsinh			n/a																	 3178.65
qfcosh			n/a																	 3304.67
qftanh			n/a																	 4159.45
qfnint        27.89		This returns result as 2s-complement uint128				   36.70	These both slower due to cleaner range checking
qfint         26.52																	   69.77

Use of all the various speedups cuts time of 1/x, sqrt(x), x/y, exp(x), tan(x) by ~2x; log(x) by ~3x; log(x) by ~2x; atan(x) by a massive 6x!
*/

/* Can turn on various qfloat-specific debugging features here: */
#ifndef QFDEBUG
	#define QFDEBUG		0
#endif

#if QFDEBUG
	#define QFONLY		0	/* Set to 1 to compile qfloat package in standalone mode - required if QFDEBUG is set! */
#else
	#define QFONLY		0	/* Set to 1 to compile qfloat package in standalone mode. */
#endif

#define DEBUG_QFMUL	0

/* Needed qfloat constants. (Remember, the high half of any qfloat has the same
	form as the IEEE double approximation to the qfloat, except for any differences
	due to rounding of the lowest bit in the latter.)
	Here's how to do this for any desired constant using the Unix 'bc' utility:

	0) Invoke 'bc -l' (i.e. bc in floating-point mode). Set (say) 75 decimal digits
	of precision and specify binary output, as follows:

	scale=75
	obase=2

	1) Generate the desired constant. We'll use 2*pi as out example, and we'll get
	this value as 8*arctan(1), which in bc syntax is:

	8*a(1)

	The output is

	110.010010000111111011010101000100010000101101000110000100011010011000\
	1001100011001100010100010111000000011011100000111001101000100101001000\
	0001001001110000010001000101001100111110011000111010000000010000010111\
	01111101010011000111011000100111001101011

	2) If the number in question is nonnegative, add the number of digits to the left of the binary point to the
	hex constant 0x3FE (= 1022 in decimal). In our example, we add 3, getting 0x401. If the number is negative,
	add to 0xBFE. The result is the 12-bit sign/exponent field of the qfloat. If the binary representation of the
	number has no nonzero digits to the left of the binary point, count the number of leading zeros to the right
	of the binary point and subtract this number from the appropriate one of the above two hex constants.
	(If the result is <= 0, your constant has underflowed, and special steps must be taken to represent it properly -
	see any reference on the IEEE floating-point standard for how to do this; most computer architecture texts have
	a section on this.)

	3) Delete the binary point and any leading zeros from the bc printout, and also delete the leading (leftmost)
	ones bit - this is the so-called "hidden bit" of the IEEE representation. Invoke another bc run, this time
	without the -l flag (i.e. now work in whole-number mode), and set the input base to be 2 by entering

	ibase=2

	Now cut the leading 52 bits left over after the above hidden-bit deletion and paste them into the bc shell.
	Then cut and paste the next 64 bits. If the leading bit of the remaining bits = 1, add one to the 64 bits
	that form the lower half of our qfloat. In our example, the leading 52 bits (sans hidden bit) are

	1001001000011111101101010100010001000010110100011000

	and the bc decimal output is

	2570638124657944 .

	The next 64 bits are

	0100011010011000100110001100110001010001011100000001101110000011

	or, in decimal:

	5086983782422027139 .

	The leftover bits are

	1001101000100101001000...,

	which have a leading ones bit, so we round by adding one to 5086983782422027139, getting 5086983782422027140 .


	4) This next step is separate from (3) since I've yet to find a bc implementation that allows one to switch
	both the input base and the output base from their defaults of 10 and correctly print things. That means that
	to convert binary to hex (which is what we are doing) we need to do 2 steps: (i) binary to decimal;
	(ii) decimal to hex. So, kill your bc run of step (3) and start another, again working in whole-number mode.
	Set the output format to hex via

	obase=16

	and paste in the two decimal outputs of step (3). In our example, the upper 52 bits of our significand are,
	in decimal:

	2570638124657944 .

	and the bc hex output is

	921FB54442D18 .

	Left-zero-pad (if necessary) the upper-52-bit result to be 13 hex characters long and concatenate
	the resulting string with the sign/exponent string you got in step (2). Thus, 0x401921FB54442D18
	is the high half of the qfloat approximation to 2*pi. The next 64 bits have decimal form

	5086983782422027140

	which in hex is

	469898CC51701B84 .

	This (left-zero-padded and appended to "Ox") is the hex representation of the lower half of our qfloat.
*/

const struct qfloat QZRO    = {0x0000000000000000ull, 0x0000000000000000ull};
const struct qfloat QONE    = {0x3FF0000000000000ull, 0x0000000000000000ull};
const struct qfloat QEPS    = {0x3890000000000000ull, 0x0000000000000000ull};	/* Set QEPS to 2^-118 */
const struct qfloat QHALF   = {0x3FE0000000000000ull, 0x0000000000000000ull};
const struct qfloat QTWO    = {0x4000000000000000ull, 0x0000000000000000ull};
const struct qfloat QTHREE  = {0x4008000000000000ull, 0x0000000000000000ull};
const struct qfloat Q2PI    = {0x401921FB54442D18ull, 0x469898CC51701B84ull};	// 2*Pi
const struct qfloat QPI     = {0x400921FB54442D18ull, 0x469898CC51701B84ull};	// Pi
const struct qfloat QPIHALF = {0x3FF921FB54442D18ull, 0x469898CC51701B84ull};	// Pi/2
const struct qfloat QPI4TH  = {0x3FE921FB54442D18ull, 0x469898CC51701B84ull};	// Pi/4
const struct qfloat QIPIHLF = {0x3FE45F306DC9C882ull, 0xA53F84EAFA3EA69Cull};	// 2/Pi
const struct qfloat QLN2    = {0x3FE62E42FEFA39EFull, 0x35793C7673007E5Full};	// log(2)
const struct qfloat QEXP    = {0x4005BF0A8B145769ull, 0x5355FB8AC404E7A8ull};	// E
const struct qfloat QSQRT2  = {0x3FF6A09E667F3BCCull, 0x908B2FB1366EA958ull};	// NB: Use Pari to output sqrt(2), then display in hex using 'bc -l' and 'obase=16' gives
const struct qfloat QISRT2  = {0x3FE6A09E667F3BCCull, 0x908B2FB1366EA958ull};	// 1.6A09E667F3BCC908B2FB1366EA957D4, which is already normalized with hidden bit left of the '.'
const struct qfloat QI64AGM = {0x403D1FB7DBFA8B50ull, 0x9906EE26BBF22BBBull};	/* 1/AGM(1,1/2^64) for log(x), @full precision */
const struct qfloat QSNAN   = {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull};	/* signaling NaN */

// Array of precomputed small-integer inverses
// !$%^%$@ GCC/Clang would not allow me to set the first 2 elts = QSNAN, QONE: "error: initializer element is not a compile-time constant"
static const struct qfloat QNINV[] = {
	{0xFFFFFFFFFFFFFFFFull,0xFFFFFFFFFFFFFFFFull},	// n =  0; set inv = QSNAN
	{0x3FF0000000000000ull,0x0000000000000000ull},	// n =  1
	{0x3FE0000000000000ull,                 0ull},	// n =  2
	{0x3FD5555555555555ull,0x5555555555555555ull},	// n =  3
	{0x3FD0000000000000ull,                 0ull},	// n =  4
	{0x3FC9999999999999ull,0x9999999999999999ull},	// n =  5
	{0x3FC5555555555555ull,0x5555555555555555ull},	// n =  6
	{0x3FC2492492492492ull,0x4924924924924925ull},	// n =  7
	{0x3FC0000000000000ull,                 0ull},	// n =  8
	{0x3FBC71C71C71C71Cull,0x71C71C71C71C71C7ull},	// n =  9
	{0x3FB9999999999999ull,0x9999999999999999ull},	// n = 10
	{0x3FB745D1745D1745ull,0xD1745D1745D1745Dull},	// n = 11
	{0x3FB5555555555555ull,0x5555555555555555ull},	// n = 12
	{0x3FB3B13B13B13B13ull,0xB13B13B13B13B13Bull},	// n = 13
	{0x3FB2492492492492ull,0x4924924924924924ull},	// n = 14
	{0x3FB1111111111111ull,0x1111111111111111ull},	// n = 15
	{0x3FB0000000000000ull,                 0ull},	// n = 16
	{0x3FAE1E1E1E1E1E1Eull,0x1E1E1E1E1E1E1E1Eull},	// n = 17
	{0x3FAC71C71C71C71Cull,0x71C71C71C71C71C7ull},	// n = 18
	{0x3FAAF286BCA1AF28ull,0x6BCA1AF286BCA1AFull},	// n = 19
	{0x3FA9999999999999ull,0x9999999999999999ull},	// n = 20
	{0x3FA8618618618618ull,0x6186186186186186ull},	// n = 21
	{0x3FA745D1745D1745ull,0xD1745D1745D1745Dull},	// n = 22
	{0x3FA642C8590B2164ull,0x2C8590B21642C860ull},	// n = 23
	{0x3FA5555555555555ull,0x5555555555555555ull},	// n = 24
	{0x3FA47AE147AE147Aull,0xE147AE147AE147AEull},	// n = 25
	{0x3FA3B13B13B13B13ull,0xB13B13B13B13B13Bull},	// n = 26
	{0x3FA2F684BDA12F68ull,0x4BDA12F684BDA12Full},	// n = 27
	{0x3FA2492492492492ull,0x4924924924924924ull},	// n = 28
	{0x3FA1A7B9611A7B96ull,0x11A7B9611A7B9612ull},	// n = 29
	{0x3FA1111111111111ull,0x1111111111111111ull},	// n = 30
	{0x3FA0842108421084ull,0x2108421084210842ull},	// n = 31
	{0x3FA0000000000000ull,                 0ull},	// n = 32
	{0x3F9F07C1F07C1F07ull,0xC1F07C1F07C1F07Cull},	// n = 33
	{0x3F9E1E1E1E1E1E1Eull,0x1E1E1E1E1E1E1E1Eull},	// n = 34
	{0x3F9D41D41D41D41Dull,0x41D41D41D41D41D4ull},	// n = 35
	{0x3F9C71C71C71C71Cull,0x71C71C71C71C71C7ull},	// n = 36
	{0x3F9BACF914C1BACFull,0x914C1BACF914C1BBull},	// n = 37
	{0x3F9AF286BCA1AF28ull,0x6BCA1AF286BCA1AFull},	// n = 38
	{0x3F9A41A41A41A41Aull,0x41A41A41A41A41A4ull},	// n = 39
	{0x3F99999999999999ull,0x9999999999999999ull},	// n = 40
};

#if QFONLY	/* If set QFDEBUG to true, must compile in standalone mode. */
int main()
{
	return qtest();
}
#endif


/* qfloat negation. */
struct qfloat qfneg(struct qfloat q)
{
	q.hi ^= MASK_SIGN;
	return q;
}


/* qfloat absolute value. */
struct qfloat qfabs(struct qfloat q)
{
	q.hi &= ~MASK_SIGN;
	return q;
}

uint32 qfiszero(const struct qfloat q)
{
	// Cannot use qfcmpeq(q, QZRO) here due to +-0 cases:
	return ((q.hi & ~MASK_SIGN) == 0) && (q.lo == 0);
}

/* qfloat comparison: q1 == q2. 0.0 and -0.0 are treated as equal; cf. http://en.wikipedia.org/wiki/Signed_zero */
uint32 qfcmpeq(struct qfloat q1, struct qfloat q2)
{
	// +-0.0 compare as equal:
	if( qfiszero(q1) && qfiszero(q2) ) {
		return TRUE;
	}
	return (q1.hi == q2.hi && q1.lo == q2.lo);
}

/* qfloat comparison: q1 != q2 */
uint32 qfcmpne(struct qfloat q1, struct qfloat q2)
{
	return !qfcmpeq(q1,q2);
}

/* qfloat comparison: q1 < q2 */
uint32 qfcmplt(struct qfloat q1, struct qfloat q2)
{
	uint32 sgn1 = (int64)q1.hi < 0, sgn2 = (int64)q2.hi < 0, both_signs = (sgn1<<1) + sgn2;
	// +-0.0 compare as equal:
	if( qfiszero(q1) && qfiszero(q2) ) {
		return FALSE;
	}
	switch(both_signs)
	{
	case(0) :	/* Both q1 and q2 nonnegative */
		return (q1.hi < q2.hi || (q1.hi == q2.hi && q1.lo < q2.lo));
	case(1) :	/* q1 >= 0, q2 <= 0, but not both zero: */
		return 0;
	case(2) :	/* q1 <= 0, q2 >= 0, but not both zero: */
		return 1;
	case(3) :	/* Both q1 and q2 negative, in which case a more-negative q1 looks larger w.r.to the unsigned compare */
		return (q1.hi > q2.hi || (q1.hi == q2.hi && q1.lo > q2.lo));
	default:
		ASSERT(0,"ERROR 98 in qfloat.c");
	}
	return 0;	/* Just to get the compiler to shut up ... this should never be reached. */
}

/* qfloat comparison: q1 > q2 */
uint32 qfcmpgt(struct qfloat q1, struct qfloat q2)
{
	return qfcmplt(q2,q1);
}

/* qfloat comparison: q1 <= q2 */
uint32 qfcmple(struct qfloat q1, struct qfloat q2)
{
	return !qfcmpgt(q1,q2);
}

/* qfloat comparison: q1 >= q2 */
uint32 qfcmpge(struct qfloat q1, struct qfloat q2)
{
	return !qfcmplt(q1,q2);
}

/*******************************************************************************/
/*                  Type Conversion Utilities                                  */
/*******************************************************************************/

// Returns qfdbl result as a 64-bit integer bitfield
uint64 qfdbl_as_uint64(const struct qfloat q)
{
	return q.hi + (q.lo >> 63);
}

/*
Return IEEE64-compliant floating double approximation to a qfloat.
Since the high part of a qfloat is already in IEEE double form, we only
need to round in the high bit of the low part and cast the result to a double.
*/
double qfdbl(struct qfloat q)
{
	uint64 hi;

	hi   = q.hi + (q.lo >> 63);
	/* If adding the carry rippled into the bottom bit of the exponent field, need to right-shift
	the significand one place. Note that this can only happen if bits <0:52> of hi were all flipped
	to 0, so don't need to worry about rounding the bit that gets shifted off. In fact, we don't even
	need to explicitly restore any hidden bit or do an actual right-shift while doing this, since:

	a) If number was normalized, the 1 that got carried out of the bottom 52 bits
		gets added to the exponent, which is what we wanted to do anyway;

	b) If number was denormalized, MSB of significand now gets treated as exponent
		field of 1, which again is what we wanted to do anyway.
	*/
	return u64_to_f64(hi);
}

// Same as above, but deliberately round the LSB the wrong way, e.g. for sensitivity analysis of a const doubles:
double qfdbl_wrong_way_rnd(struct qfloat q)
{
	uint64 hi,lo;
	lo = (q.lo >> 63) ^ 0x1;
	hi = q.hi + lo;
	return u64_to_f64(hi);
}

// For generic doubles (i.e. we do not know what the less-significant bits were which were rounded in),
// simply provide a function which toggles the LSB:
double dbl_flip_lsb(double d)
{
	uint64 c = f64_to_u64(d);
	BIT_FLIP(c,0);
	return u64_to_f64(c);
}

#ifdef X87_ASM
/* qfloat --> long double conversion utility.
Example: x = QLN2:
	hi = 0x3fe62e42fefa39ef,
	lo = 0x35793c7673007e5f, low 64 bits are all-extended-mantissa bits, separate into hi11 and lo53:
	= 00110101011 11001001111000111011001110011000000000111111001011111, lo53 have MSB = 1, so add 1 to hi11:
	= 001 1010 1100 [rest discarded], 0-pad 1 bit at top, vie rest as 3-digit hex,
	upper 11 bits (after rounding) = 0001 1010 1100 = 1AC
==> Top 64 mant-bits = (162E42FEFA39EF)*2^B + 1AC = B17217F7D1CF79AC = 12786308645202655660_10 .
	Input exp-bits = 0x3FE, get widened to 0x3FFE .
*/
long double qfldbl(struct qfloat x)
{
	long double ld;
	uint64 *ld_ptr = (uint64 *)&ld, nonhidden;
	int32 exp = (int32)((x.hi & ~MASK_SIGN)>>52);
	ASSERT(sizeof(long double) == 16, "QFLDBL assumes 16-byte long double type!");
	// Denormal check:
	ASSERT((exp != 0) && (exp != 0x7ff), "QFLDBL requires normal input!");
	exp -= (int32)0x400;	// x87 80-bit reg-format has 4 more bits in exp, centered around 0x4000 rather than 0x400
	nonhidden = ((x.hi & MASK_MANT)<<11) + (x.lo>>53) + ((x.lo>>52)&0x1);
	// Rounding of the off-shifted portion may cause nonhidden-bit summation to overflow into sign bit:
	if(MASK_SIGN == nonhidden) {
		*ld_ptr = MASK_SIGN;	++exp;		// which requires special handling
	} else {
		*ld_ptr = MASK_SIGN + nonhidden;	// 64-bit x87 mantissa, MASK_SIGN used here to instantiate (non)hidden top bit
	}
	*(ld_ptr+1) = (uint64)((int32)0x4000 + exp) + ((x.hi & MASK_SIGN)>>48);	// Widen exp-field and restore sign
	return ld;
}

/* long double --> qfloat conversion utility */
struct qfloat ldbl_to_q(long double x)
{
	struct qfloat q;
	long double ld = x;
	uint64 *ld_ptr = (uint64 *)&ld, x87_mant, x87_sexp;	// Note high 48 bits of x87_sexp are uninited
	int32 exp;
	DBG_ASSERT(sizeof(long double) == 16, "LDBL_TO_Q assumes 16-byte long double type!");
	x87_mant = *ld_ptr; x87_sexp = *(ld_ptr+1);
	if(!x87_mant) return QZRO;
	// Denormal check:
	exp = (int32)(((x87_sexp<<48) & ~MASK_SIGN)>>48) - (int32)0x4000;
	ASSERT(ABS(exp) <= 0x3ff, "LDBL_TO_Q requires double-compatible normal input!");
	q.hi = ((x87_sexp<<48) & MASK_SIGN) + ((uint64)((int32)0x400 + exp)<<52) + ((x87_mant>>11) & MASK_MANT);
	q.lo = (x87_mant<<53);
	return q;
}
#endif


/* Convert IEEE64 double to qfloat. */
struct qfloat dbl_to_q(double d)
{
	struct qfloat q;

	q.hi = f64_to_u64(d);	/* Copy bit pattern of double into a uint64. */
	q.lo = (uint64)0;
	return q;
}


/* Convert 64-bit signed int to qfloat. */
struct qfloat i64_to_q(int64 i64)
{
	struct qfloat q;
	int32 lz, shift;
	uint64 sexp;
	if(!i64) return QZRO;

	sexp = i64 & MASK_SIGN;
	if(sexp) i64 = -i64;
	lz = leadz64(i64);	/* Find # of leading zeros in |i64|. */

	/* Put leading ones bit of mantissa into bit position 52 of <0:63> and add to sign/exponent. */
	shift = lz - (int)11;	/* if lz = 11, int is already correctly positioned. If <> 11,
							i.e. shift <> 0, need to right/left-shift mantissa */

	sexp += ((uint64)1074 - shift) << 52;	/* Ex: i64 = 3, with lz = 62 and shift = 51, should yield exponent = 0x400 = 1024
											(prior to left-shift by 52 places); each bit less in lz means one more in exponent.
											Since the leftmost mantissa bit = 1 and gets added to the exponent, subtract an
											extra one, i.e. 1075 - 1 - shift = 1074 - shift. */
	if(shift < 0)
	{
		q.hi = sexp + (i64 >> (-shift));
		q.lo = i64 << (64+shift);
	}
	else	/* need to left-shift mantissa */
	{
		q.hi = sexp + (i64 << shift);
		q.lo = (uint64)0;
	}

	return q;
}

/* Convert 128-bit signed int to qfloat. We formally pass a uint128, but signedness is inferred from 2s-comp representation: */
struct qfloat i128_to_q(uint128 i)
{
	struct qfloat q;
	int32 lz, shift, rshift,lshift;
	uint64 sexp, offword;
	if(!i.d1 && !i.d0) return QZRO;

	sexp = i.d1 & MASK_SIGN;
	if(sexp) {	// General formula for 2s-comp arithmetic negation is -i = ~i + 1
		i.d0 = ~i.d0; i.d1 = ~i.d1;
		ADD128(i, ONE128, i);
	}
	lz = LEADZ128(i);	/* Find # of leading zeros in |i|. */

	/* Put leading ones bit of mantissa into bit position 52 of q.d1 and add to sign/exponent. */
	shift = lz - (int)11;	/* if lz = 11, int is already correctly positioned. If <> 11,
							i.e. shift <> 0, need to right/left-shift mantissa */

	sexp += ((uint64)1074 - shift) << 52;	/* Ex: i = 3, with lz = 62 and shift = 51, should yield exponent = 0x400 = 1024
											(prior to left-shift by 52 places); each bit less in lz means one more in exponent.
											Since the leftmost mantissa bit = 1 and gets added to the exponent, subtract an
											extra one, i.e. 1075 - 1 - shift = 1074 - shift. */
	if(shift < 0)
	{
		// Unlike the 64-bit case, need to properly round any off-shifted bits, but defer for later due to complications
		// this adds (mainly, that the round-carry addition can ripple-carry all the way up, requiring an added mantissa-shift).
		rshift = -shift; lshift = (64+shift);
		q.hi = sexp + (i.d1 >> rshift);
		offword = (i.d0 << lshift) >> 63;	// MSB of off-shifted portion
		q.lo = (i.d1 << lshift) + (i.d0 >> rshift) + offword;
		ASSERT(q.lo >= offword, "Ripple-carry!");	// For now, just check for ripple-carry. Proper handling will come later.
	}
	else	/* need to left-shift mantissa */
	{
		LSHIFT128(i,shift,i);
		q.hi = sexp + i.d1;
		q.lo = i.d0;
	}

	return q;
}


/* Multiply by an integer power of 2 is especially simple: */
struct qfloat qfmul_pow2(struct qfloat q, int32 pow)
{
	int64 exp0,exp1, sgn = q.hi & MASK_SIGN, lz, shft;
	struct qfloat qt = q;
	if(!pow) return qt;	// *2^0 is a no-op
	// Extract 11-bit exponent field:
	exp0 = ((q.hi - sgn) >> 52);
	exp1 = exp0 + pow;	// Exponent field after shift, in absence of any denormal considerations
	// Underflow: Set exp = 0 and rshift mantissa by remaining bits
	if(exp1 <= 0) {
		WARN(HERE,"DENORM in qfmul_pow2", "",0);
		// If denorm to begin with, we right-shift the mantissa a further (pow) bits:
		if(exp0 == 0) {
			qt.hi &= MASK_MANT;	// Mantissa-only
			QRSHIFT(qt, -pow, qt);
			qt.hi += sgn;		// Restore sign
		} else {
		// If normal input, restore hidden bit, right-shift result -pow - (exp0-1) places
		// (corr. to setting the exp = 1) and set exp = 0 to signal denomral result:
			qt.hi = (qt.hi & MASK_MANT) + TWOE52;	// Mantissa-only
			shft = -pow - (exp0-1);
			QRSHIFT(qt, shft, qt);
			qt.hi += sgn;		// Restore sign
		}
	} else if(exp1 >> 11) {	// Overflow: exp+pow carried into sign-bit slpt:
		ASSERT(0,"OVERFLOW in qfmul_pow2");
	} else {	// Result is normal
		if(exp0) {
			// If normal input, update exponent field and return:
			qt.hi += ((uint64)pow << 52);
		} else {
			// If denorm to begin with, we left-shift the mantissa by the smaller of
			// (pow, [shift count needed to move leftmost mantissa bit into exp-field]) bits:
			qt.hi &= MASK_MANT;	// Mantissa-only
			if(!qfiszero(qt)) {
				lz = QLEADZ(qt) - 11;	// Number of leading zero bits in denormalized mantissa (i.e. shift count needed to move leading bit into hidden-bit location)
				if(pow > lz) {	// Result will be normal
					QLSHIFT(qt, lz, qt);
					ASSERT((qt.hi>>52) == 1, "Bad mantissa left-shift count in qfmul_pow2!");
					qt.hi += (((uint64)pow-lz)<<52) - TWOE52;	// Don't fold -TWOE52 in via (pow-lz-1)<<52, since may have pow = lz here.
				} else {	// Result still denormal
					QLSHIFT(qt, pow, qt);
					ASSERT((qt.hi>>52) == 0, "Bad mantissa left-shift count in qfmul_pow2!");
				}
			}
			qt.hi += sgn;		// Restore sign
		}
	}
	return(qt);
}


/* Nearest integer (round-to-nearest) of a qfloat. The result is stored in a uint128,
interpreted in 2s-comp form. The allowable range of the result is [-2^127, +2^127-1].
*/
uint128 qfnint(struct qfloat q)
{
	const struct qfloat two127 = {0x47E0000000000000ull, 0x0000000000000000ull};
	int32 exp, sign, rshift, lshift;
	uint64 offword, carry;
	uint128 i;
	i.d1 = q.hi; i.d0 = q.lo;
	ASSERT(qfcmpge(q, qfneg(two127)) && qfcmplt(q, two127), "QFNINT input out of range!");
	/* Separate upper part of the significand from the sign/exponent fields: */
	sign = (int32)(i.d1 >> 63);
	exp  = (int32)(i.d1 >> 52) & MASK_EXP;
	i.d1 =         i.d1        & MASK_MANT;

	/* If number is normalized, restore hidden bit. Otherwise, number is so small that we can immediately return zero. */
	if(exp) {
		i.d1 += TWOE52;
	} else {
		return NIL128;
	}

	/* Get right-shift count. E.g. 1 has exp = 0x3FF and rshift = 116, 2 has exp = 0x400 and rshift = 115,
	so rshift = 0x400 + 115 - exp = 0x473 - exp.
	*/
	rshift = 0x473 - exp;	// Note (qfloat)2^127 has exp = 0x3ff + 0x07f = 0x47e --> rshift = -11.
	carry = (uint64)0;

	if(rshift <= 0)		/* If rshift <= 0, require rshift strictly > -11, unless rshift = -11 and result = -2^127. */
	{
		lshift =     - rshift;
		rshift = (64 + rshift) & 63;	/* shift count of 64 must get aliased to 0. */
		i.d1 = (i.d1 << lshift) + (i.d0 >> rshift);
		i.d0 = (i.d0 << lshift);
	}
	else if(rshift <= 53)	/* Hi part partially shifted into lo, lo part partially shifted off. */
	{
		lshift = 64 - rshift;
		offword = (i.d0 << lshift) >> 63;	/* Round by adding the MSB of off-shifted low word into next-higher word. */
		i.d0 = (i.d0 >> rshift) + (i.d1 << lshift);
		i.d0 += offword;
		if(i.d0 < offword)	/* Had a carry from the rounded bit being added into the right-shifted low part. */
		{
			carry = (uint64)1;
		}
		i.d1 = (i.d1 >> rshift) + carry;
	}
	else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into i.d1, lo part partially shifted off. */
	{
		lshift = 64 - rshift;
		offword = (i.d0 << lshift) >> 63;
		i.d0 = (i.d0 >> rshift) + (i.d1 << lshift) + offword;
		i.d1 = (uint64)0;
	}
	else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into i.d1, lo part partially shifted off. */
	{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
		offword = i.d0 >> 63;
		i.d0 = i.d1 + offword;
		i.d1 = (uint64)0;
	}
	else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
	{
		rshift -= 64;
		lshift = 64 - rshift;
		offword = (i.d1 << lshift) >> 63;
		i.d0 = (i.d1 >> rshift) + offword;
		i.d1 = (uint64)0;
	}
	else			/* Entire significand shifted off, no need to round result. */
	{
		return NIL128;
	}

	/* If sign bit was set, negate result by arithmetic negation of lo word, logical negation of hi word. */
	if(sign)
	{
		i.d1 = ~i.d1;
		i.d0 = -i.d0;
	}

	return i;
}

/* Cast to integer (round-toward-zero) of a qfloat. The result is stored in a uint128,
interpreted in 2s-comp form. The allowable range of the result is [-2^127, +2^127-1].
*/
uint128 qfint(struct qfloat q)
{
	const struct qfloat two127 = {0x47E0000000000000ull, 0x0000000000000000ull};
	int32 exp, sign, rshift, lshift;
	uint128 i;
	i.d1 = q.hi; i.d0 = q.lo;
	ASSERT(qfcmpge(q, qfneg(two127)) && qfcmplt(q, two127), "QFNINT input out of range!");

	/* Separate upper part of the significand from the sign/exponent fields: */
	sign = (int32)(i.d1 >> 63);
	exp  = (int32)(i.d1 >> 52) & MASK_EXP;
	i.d1 =         i.d1        & MASK_MANT;

	/* If number is normalized, restore hidden bit. Otherwise, number is so small that we can immediately return zero. */
	if(exp) {
		i.d1 += TWOE52;
	} else {
		return NIL128;
	}

	/* Get right-shift count. E.g. 1 has exp = 0x3FF and rshift = 116, 2 has exp = 0x400 and rshift = 115,
	so rshift = 0x400 + 115 - exp = 0x473 - exp.
	*/
	rshift = 0x473 - exp;

	if(rshift <= 0)		/* If rshift <= 0, require rshift strictly > -11, unless rshift = -11 and result = -2^127. */
	{
		if(rshift <= -11)
		{
			if(rshift == -11 && (!sign || (i.d1 << -rshift) != MASK_SIGN || i.d0 != (uint64)0))
			{
				ASSERT(0,"ERROR: qfloat is too large to convert to 128-bit integer.");
			}
		}
		lshift =     - rshift;
		rshift = (64 + rshift) & 63;	/* shift count of 64 must get aliased to 0. */
		i.d1 = (i.d1 << lshift) + (i.d0 >> rshift);
		i.d0 = (i.d0 << lshift);
	}
	else if(rshift <= 53)	/* Hi part partially shifted into lo, lo part partially shifted off. */
	{
		lshift = 64 - rshift;
		i.d0 = (i.d0 >> rshift) + (i.d1 << lshift);
		i.d1 = (i.d1 >> rshift);
	}
	else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into i.d1, lo part partially shifted off. */
	{
		lshift = 64 - rshift;
		i.d0 = (i.d0 >> rshift) + (i.d1 << lshift);
		i.d1 = (uint64)0;
	}
	else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into i.d1, lo part partially shifted off. */
	{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
		i.d0 = i.d1;
		i.d1 = (uint64)0;
	}
	else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
	{
		rshift -= 64;
		lshift = 64 - rshift;
		i.d0 = (i.d1 >> rshift);
		i.d1 = (uint64)0;
	}
	else			/* Entire significand shifted off, no need to round result. */
	{
		return NIL128;
	}

	/* If sign bit was set, negate result by arithmetic negation of lo word, logical negation of hi word. */
	if(sign)
	{
		i.d1 = ~i.d1;
		i.d0 = -i.d0;
	}

	return i;
}

/******************** I/O functions ******************************/

/* Convert a qfloat to normalized scientific/exponential notation, with respective subfields char counts listed beneath:

[sign][mantissa in range 1.0 <= m < 10.0, with 36 decimal digits of displayed significance] E[sign][power of 10] .
  1                               37 (incl. decimal point)                                  2  1        3

Total = 1+36+2+1+3 = 43 characters, so return a length-45 character string, incl. terminating nullchar.

Algorithm proceeds by copying mantissa (with hidden bit made explicit for normal inputs) to an mi64 array of length 2^11 bits
or 2^5 64-bit unsigned integer words, with the binary radix point taken to divide the upper and lower halves of this array.Thus
the number 1.0 ends up as the array (where we write the MSW at left and significance decreases rightward)

1.0 = [31 words = 0] , 1 <.> , [32 trailing words = 0] .

Thus, our initial bitstring is of form [117-bit mantissa] << (exp-116), for both normal and denormal inputs.

Example: qfloat Pi = 0x400921FB54442D18 469898CC51701B84, where the space separates the left (hi) and right(lo) 64-bit words.

Separating sign/exp, restoring the hidden bit and left-shifting 0x400-115 = 909 places gives mantissa field
 = 11.0010010000111111011010101000100010000101101000110000100011010011000100110001100110001010001011100000001101110000100 ,
which has 2 nozero bits left of the radix point, the remaining 115 mantissa bits to the right of the radix point, and
1024-2 = 1022 leading binary zeros and 1024-115 = 909 trailing binary zeros bracketing the 117 mantissa bits.

After this initial mantissa shift, we next adjust this initial bitstring by successive muls-by-10 (if the initial value < 1.0)
or divs-by 10 (if >= 10.0) in order to produce a result lying in the range [1.0, 10). The count of muls-or-divs-by-10 needed
to accomplish this yields the power-of-10 exponent appearing in our final scientific-notation result, as E+[count] if we used
muls-by-10, and E-[count] is we used divs. (By convention if we needed no muls or divs we write E+000).

In our Pi example, we need no such power-of-10 adjustment.

We then process our mantissa bits by dumping the current ones to the left of the radix point into a decimal digit, zeroing them,
multiplying the remaining bitstring by 10, dumping the resulting left-of0radix-point bits into the next-lower decimal digit, etc.
We do this until we have filled our 36-decimal-digit display mantissa. The 36th is in fact at most semi-significant, but include
it so user can infer properly rounded 35th digit.

Again for our Pi example, here are the resulting decimal digits and leftover bitstrings:

3	0011.0010010000111111011010101000100010000101101000110000100011010011000100110001100110001010001011100000001101110000100
1	0001.0110101001111010001010010101010100111000010111100101100000111110101111101111111101100101110011000010001001100101000
4	0100.0010100011000101100111010101010000110011101011110111001001110011011101011111100111111001111110010101011111110010000
1	0001.1001011110111000001001010100101000000100110110100111100010000010100110111100001111000011101111010110111101110100000
5	0101.1110110100110001011101001110010000110000100010001011010100011010000101011010010110100101011001100101101010001000000
9	1001.0100001111101110100100001110100111100101010101110001001100000100110110000111100001110101111111111000100101010000000
2	0010.1010011101010001101010010010001011110101011001101011111000110000011101001011010010011011111110110101110100100000000
6	0110.1000100100110000100110110101110110010110000000110110110111100100100011110000111000010111110100011010001101000000000
5	0101.0101101111100110000100011010011111011100001000100100101011101101100101101000110011101110001100000110000010000000000
3	0011.1001011011111100101100001000111010011001010101101110110101000111111000011000000101001101111000111100010100000000000

5	0101.1110010111011110111001011001000111111101011001010100010011001110110011110000110100001010111001011011001000000000000
8	1000.1111101010110100111101111011001111100101111101001011000000010100000101101000001001101100111110001111010000000000000
9	1001.1100101100010001101011010000011011111011100011101110000011001000111000010001100001000001101110011000100000000000000
7	0111.1110111010110000110000100100010111010011100101001100011111011000110010101111001010010001001111110101000000000000000
9	1001.0101001011100111100101101011101001000011110011111100111001110111111011010111100110101100011110010010000000000000000
3	0011.0011110100001011111000110100011010100110000111100001000010101111010001101100000010111100101110110100000000000000000
2	0010.0110001001110110111000001100001001111101001011001010011011011000110000111000011101011111010100001000000000000000000
3	0011.1101100010100100110001111001100011100011101111101000010001110111101000110100100110111001001001010000000000000000000
8	1000.0111011001101111110010111111100011100101011100010010110010101100011000001110000100111011011100100000000000000000000
4	0100.1010000001011101111101111011100011110110011010111011111010111011110010001100110001010010011101000000000000000000000

6	0110.0100001110101011101011010011100110100000001101010111001101010101110101111111101100111000100010000000000000000000000
2	0010.1010010010110100110001000100000001000010000101101000000101011010011011111101000000110101010100000000000000000000000
6	0110.0110111100001111101010101000001010010100111000010000110110001000010111100010001000010101001000000000000000000000000
4	0100.0101011010011100101010010001100111010000110010101000011101010011101011010101010011010011010000000000000000000000000
3	0011.0110001000011110100110110000001000100111111010010100100101000100110001010101000001000000100000000000000000000000000
3	0011.1101010100110010000011100001010110001111000111001101110010101111101101010010001010000101000000000000000000000000000
8	1000.0101001111110100100011001101011110010111001000001001111011011101000100110101100100110010000000000000000000000000000
3	0011.0100011110001101100000000110101111100111010001100011010010100010110000010111101111110100000000000000000000000000000
2	0010.1100101110000111000001000011011100001000101111100000111001011011100011101101011110001000000000000000000000000000000
7	0111.1111001101000110001010100010011001010111011011001000111110010011100101000110101101010000000000000000000000000000000

9	1001.1000000010111101101001010111111101101010001111011001101111000011110011000011000100100000000000000000000000000000000
5	0101.0000011101101000011101101111101000100110011010000001010110100101111110011110101101000000000000000000000000000000000
0	0000.0100101000010100101001011100010110000000000100001101100001111011110000110011000010000000000000000000000000000000000
2	0010.1110010011001110011110011011011100000000101010000111010011010101100111111110010100000000000000000000000000000000000
8	1000.1111000000010000110000010010011000000110100101001001000001011000001111101111001000000000000000000000000000000000000

and we should round up the trailing 8 to a 9 because the leading bit of the remaining stuff to the right of the radix point = 1.
If we compute the next few digits we get

9	1001.0110000010100111100010110111110001000001110011011010001101110010011101010111010000000000000000000000000000000000000
3	0011.1100011010001011011100101101101010010010000010000110001001111000100101101000100000000000000000000000000000000000000

of which at most the '9' has any significance, and indeed the corresponding 7 digits of exact Pi = ...9502884 . We return
...950289, which is correct in the sense of a 35-digit-rounded result ...95029 .
*/
char* qf2str(struct qfloat q)
{
	static char os[45];
	int i, len = 32, exp, sgn = (int64)q.hi < 0, shft, pow10 = 0, mod10;
	static uint64 *u = 0x0;	// Scratch array
	const char pm[2] = {'+','-'};
	os[0] = pm[sgn];	os[1] = '\0';
	if(!u) {	// Allocate scratch space
		u = (uint64 *)calloc(32, sizeof(uint64));
	}
	mi64_clear(u,32);
	if(qfiszero(q)) {	// +-0.0
		strcat(os,"0.00000000000000000000000000000000000E+000");
		return os;	// *2^0 is a no-op
	}
	// Extract 11-bit exponent field and use to compute binary mantissa-shift count:
	exp = (int)((q.hi >> 52) & MASK_EXP);
	u[0] = q.lo; u[1] = q.hi & MASK_MANT;
	if(exp) u[1] += TWOE52;	// If normal, restore hidden bit
	shft = exp - 115;	// = exp - [number of explicit mantissa bits]
	if(shft > 0) {
		mi64_shl(u,u,shft,len);
	}
	if(exp > 0x401) {			// Need 1 or more div-by-10
		while(mi64_getlen(u,len) > 17 || u[16] >= 10) {
			mi64_div_by_scalar64(u, (uint64)10, len, u);	// For now, just discard remainder
			++pow10;
		}
	} else if(exp < 0x3FF) {	// Need 1 or more mul-by-10
		while(mi64_getlen(u,len) < 17 || !u[16]) {
			mi64_mul_scalar     (u, (uint64)10, u, len);	// No need to check for carryout
			--pow10;
		}
	}
	ASSERT(mi64_getlen(u,len) == 17 && u[16] && u[16] < 10, "QF2STRING: Normalization error!");
	os[1] = u[16] + '0';	// Put MSD to left of decimal point
	os[2] = '.';
	for(i = 3; i < 38; ++i) {
		u[16] = mi64_mul_scalar(u+14, (uint64)10, u+14, 2);	// Only use 128-bits to right of radix point
		os[i] = u[16] + '0';
	}
	// Now do decimal-form exponent:
	os[i++] = ' ';
	os[i++] = 'E';
	if(pow10 < 0) {
		os[i++] = '-';
	} else {
		os[i++] = '+';
	}
	pow10 = ABS(pow10);
	mod10 = pow10%10; os[i+2] = '0' + mod10;	pow10 /= 10;
	mod10 = pow10%10; os[i+1] = '0' + mod10;	pow10 /= 10;
	mod10 = pow10%10; os[i+0] = '0' + mod10;	pow10 /= 10;
	os[44] = '\0';
	return os;
}

/************************************* Arithmetic functions ****************************/

/*
Qfloat multiply algorithm. Needs only IEEE64-floating-compliant multiply and a macro to generate the
128-bit product of two unsigned 64-bit ints.

Rearrange the inputs as x = b + a*2^53, y = d + c*2^53, where a and c are 53 bits long,
and b and d are 64 bits long, in contrast with the obvious approach, which keeps 64 bits into the low part
of the 2-word representation and 53 bits into the high part. The multiply then looks like

        x*y = a*c       *2^106   (upper 128 bits)
            +(a*d + c*b)*2^53    (next-lower 117 bits, right-shifted 53 bits w.r.to a*c)
            + b*d                (lowest-order 106 bits, right-shifted 106 bits w.r.to a*c).

Notice what happens: superficially this looks like the same amount of work to calculate,
but hold on - since the middle (a*d + c*b) term of the above product is right-shifted 53 bits w.r.to
the leading order (a*c) term, only the upper 64 bits of a*d and c*b (considered separately) overlap
the a*c term, i.e. they align precisely with the lower half of a*c. Even better, since we only need
the upper 53 bits of these 3 64-bit aligned terms (i.e. MULL64(a*c), MULH64(a*d) and MULH64(c*b)) and since
MULH64(a*d) and MULH64(c*b) have only 53 bits to begin with, we can simply use an FMUL to emulate these
two product operations. Unless we absolutely need a result in which all 117 bits are correct, we don't
even care if the bottom bit of the FMUL result is wrong, i.e. need no expensive and messy error correction
operations to correct this bit. Our total cost is now just one 128-bit integer product (b*d), 2 FMUL,
and 2 adds with a potential carry into the upper 64 bits of a*c. The latter can be simplified by first
rounding MULL64(a*c) to its upper 53 bits, adding in FMUL(a*d) and FMUL(c*b) (properly scaled), and
combining the result with the lower 11 bits of MULH64(a*c). Sweet!

TOTAL COST (uisng Alpha multiply model): 6 FLTCVT, 4 FMUL, 2 IMUL + 50-70 IOPS (shift, add, etc.)
This is a lot of total operations, but the most-expensive operation, integer multiply, has been drastically
reduced in number from the all-integer model. Assuming a typical RISC multiple-instruction-issue model and
two 64-bit integer logic units which can operate in parallel on simple IOPS, we estimate 30-40 cycles per
qfloat multiply. That's a lot slower than the one cycle for multiplying two doubles, but still fast enough
to be useful. Also, were we to use an all-floating implementation, we'd need roughly a length-16 floating
convolution, followed by a carry step, and a bunch of floating-int conversions, which I estimate would be
~4x slower. (The latter approach is of course asymptotically better since we can use FFT, but that doesn't
begin to pay off until runlengths become much larger than the # of bits in a qfloat.)

IMPORTANT NOTE: since the input mantissas are in [2^116, 2^117)
(and this neglects the possibility of denormalized input operands, which requires additional special casing),
the product is in [2^232, 2^234), i.e. may have either 233 or 234 bits. We detect which by checking whether the
result of the MUL_LOHI64(b, d) has 127 or 128 bits, and then left-justify subsequent operations accordingly,
in order not to leave the MSB of the result mantissa vacant.


Example: in our qfloat format, the 117-bit approximation to e looks like

2^2 * 0.10101101111110000101010001011000101000101011101101001 0101001101010101111110111000101011000100000001001110011110101000

Rearranging the significand as b + a*2^53 gives

a = 1010110111111000010101000101100010100010101110110100101010011010 = 12535862302449814170
b = 10101111110111000101011000100000001001110011110101000            = 6187547923638184

Squaring, we get

a*a.hi = 8519001675203524399
a*a.lo = 16107666474326910116
a*b.hi = 4204874770886292	(rightmost digit should be 3 after rounding in low part).

Writing these in binary and aligning according to their relative power-of-2 multipliers:
																x
a*a*2^53 = 0111011000111001100100101110001101010011011101101011011100101111 1101111110001001111011010101000011111010101110010110010010100100
a*b      =                                                                  01110111100000100111110110011000010110010111010010100


*/
struct qfloat qfmul(struct qfloat q1, struct qfloat q2)
{
	struct qfloat return_val;
	uint64 sexp1, sexp2, hi1, hi2, lo1, lo2;
	uint64 a,b,c,d,lo,hi,bvac;
#ifdef USE_FMUL_FOR_LOW_WORD
	const double fmul_scale = 1.0/((double)0x00400000*0x80000000);	/* 2^(-53) */
	double da,db,dc,dd, adhi,bchi;
#else
	uint64 adhi,bchi;
#endif

	/* If either operand is zero, return zero. We mask off the sign bit so -0 is treated like 0. */
	if((!(q1.hi & ~MASK_SIGN) && !q1.lo) || (!(q2.hi & ~MASK_SIGN) && !q2.lo))
	{
		return QZRO;
	}

	/* Separate upper parts of the 2 significands from the sign/exponent fields: */

	hi1 = q1.hi & MASK_MANT;
	hi2 = q2.hi & MASK_MANT;

	sexp1 = q1.hi - hi1;
	sexp2 = q2.hi - hi2;
	/* If number is normalized, restore hidden bit. */
	if(sexp1 & ~MASK_SIGN)
	{
		hi1 += TWOE52;
	}
	else
	{
		printf("Multiply by denormalized operand not supported!");
	//	ASSERT(0,"ERROR in qfloat.c : qfmul");
		return QZRO;
	}

	if(sexp2 & ~MASK_SIGN)
	{
		hi2 += TWOE52;
	}
	else
	{
		printf("Multiply by denormalized operand not supported!");
	//	ASSERT(0,"ERROR in qfloat.c : qfmul");
		return QZRO;
	}

	lo1 = q1.lo;
	lo2 = q2.lo;

	/* Need to subtract scaling factor prior to adding sexp1 and sexp2 together, so don't overflow into the sign bit.
	Also check if resulting exponent sum is < 0x400, indicating a denormalized MUL result, which gets flushed to 0:
	*/
	if( ((sexp1 & ~MASK_SIGN) + (sexp2 & ~MASK_SIGN)) < 0x4000000000000000ull) {
	#if QFDEBUG
		WARN(HERE, "DENORM result in QFMUL ... flushing to 0.\n", "", 0);
		ASSERT(fabs(qfdbl(q1)*qfdbl(q2)) < 1e-300, "Incorrect DENORM result in QFMUL!");
	#endif
		return QZRO;
	}
	sexp1 -= 0x3FF0000000000000ull;	/* Hex constant = ((uint64)1023 << 52) */

	a = (hi1<<11) + (lo1>>53);
	c = (hi2<<11) + (lo2>>53);

#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(a,c,&lo,&hi);	/* start the integer product ASAP, since it will generally have a large latency. */
#else
	MUL_LOHI64(a,c,lo,hi);
#endif

	b = (lo1 << 11) >> 11;
	d = (lo2 << 11) >> 11;

#ifndef USE_FMUL_FOR_LOW_WORD

  #ifdef MUL_LOHI64_SUBROUTINE
	adhi = __MULH64(a,d);
	bchi = __MULH64(b,c);
  #else
	MULH64(a,d,adhi);
	MULH64(b,c,bchi);
  #endif

#else	// USE_FMUL_FOR_LOW_WORD

	/* Make these as accurate as possible by rounding rather than truncating.
	If operand overflows into bit position 53 (of <0:63>) on rounding, the high part was = 2^53 - 1
	and is now = 2^53, which can be exactly represented as a double, i.e. we lose no precision if this occurs.
	*/
	da = (double)(hi1 + (lo1 >> 63)) * fmul_scale;
	dc = (double)(hi2 + (lo2 >> 63)) * fmul_scale;
	db = (double)(b);
	dd = (double)(d);

  /* DEBUG: make sure we didn't lose any bits of b or d during the conversion to float. */
  #if QFDEBUG
	if(!((uint64)db == b)) ASSERT(0,"ERROR 120 in qfloat.c");
	if(!((uint64)dd == d)) ASSERT(0,"ERROR 121 in qfloat.c");
  #endif

	adhi = da*dd;
	bchi = db*dc;

#endif	// USE_FMUL_FOR_LOW_WORD

	bvac = (uint64)leadz64(hi);
	ASSERT((bvac < 2), "ERROR 130 in qfloat.c");	/* Make sure at most the leftmost bit of high part is vacant. This check
						needs to remain in place until support for denormalized oprands is added. */
	/*
	Now need to right-shift MUL_LOHI result (12-bvac) places, FMUL results (1-bvac) place, and add together.
	(We could contemplate FADDing the 2 FMUL results together, which saves a float->int conversion
	but risks losing a bit of precision at the lower end. If FLTCVT is slow, it's probaby a good trade.)
	First add the 2 FMUL results to (LO >> (11-bvac)), which preserves an extra low bit, then right-shift
	one place, round and add to the lower (12-bvac) bits of HI.
	*/
	/*
	lo = (hi << 52) + (lo >> 12) + (lo & 1) + (uint64)(adhi + bchi);
 	hi = (hi >> 12);
	printf("mul_hi = %16" PRIX64 " = %20" PRIu64 "\n", hi, hi);
	printf("mul_lo = %16" PRIX64 " = %20" PRIu64 "\n", lo, lo);
	*/

	return_val.hi = (hi >> (11-bvac));
	ASSERT((return_val.hi >> 52) == 1, "ERROR 140 in qfloat.c");
	return_val.lo = (hi << (53+bvac)) + (lo >> (11-bvac)) + ((lo >> (10-bvac)) & (uint64)1) + (((uint64)adhi + (uint64)bchi) << bvac);
                                                            /* ^^^^rounding is here^^^^^ */   /* Maximize accuracy by converting to int prior to add. */
	if(return_val.lo < (hi << (53+bvac)))	/* had a carry out of lo part. */
	{
		return_val.hi++;

		/* If adding the carry rippled into the bottom bit of the exponent field, need to right-shift the significand one place.
		(This actually happens quite frequently, e.g. in qfinv, when we multiply x by its estimated inverse, and the product may
		be extremely close to unity, which may show up as return_val.hi = 0x001FFFFFFFFFFFFF prior to adding carrylo, and then
		we get 0x0020000000000000 after. Note that we don't actually need to right-shift the hi part in this case, since the
		bottom bits are all zero and the carry needs to be added to the exponent field to account for the right-shift, anyway.
		*/
		if(return_val.hi == 0x0020000000000000ull)	/* 	Note that this can only happen if bits <0:52> of hi were all flipped to 0... */
		{
			return_val.lo = (return_val.lo >> 1) + (return_val.lo & 1);	/* ...so don't need to worry about carrying in the lo bit of hi here. */
		}
	}
	return_val.hi = (sexp1 + sexp2 - (bvac << 52)) + return_val.hi;

#if QFDEBUG
	double qres = qfdbl(return_val), dres = qfdbl(q1)*qfdbl(q2);
	double esum = fabs(qres + dres);
	if(esum > 1e-15) {
		double rerr = fabs( (qres - dres)/esum );
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFMUL!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFMUL!\n");
	}
#endif

	return(return_val);
}

/*
Increment and decrement:
*/
struct qfloat qfinc(struct qfloat x)
{
	struct qfloat q = x;
	uint32 rshift, lshift;
	uint64 exp0,hi0,hi1,offword;
	if((int64)q.hi < 0) {
		return qfneg(qfdec(q));	// For x < 0, Use x + 1 = 1 - |x| = -(|x| - 1)
	} else {
		// Check for denormal (over and underflow):
		ASSERT(((q.hi>>52) + 1) >= 2, "Denormal not supported!");

		if(q.hi < QONE.hi) {
		// This is just the significand-add section of qfsum with the following argument value specializations:
		//		lo0 = QONE.lo = 0, hi0 = (QONE.hi & MASK_MANT) + TWOE52 = TWOE52
		//		lo1 =        q.lo, hi1 =                  q.hi + TWOE52  (i.e. q.hi with hidden bit restored)
			exp0 = QONE.hi;	// No nonzero significand bits needing masking-off here
			hi0 = TWOE52;
			rshift = (QONE.hi>>52) - (q.hi>>52);	// Must shift before - to avoid borrow-into-exponent
			lshift = 64 - rshift;
			// Mask off sign/exp-bits (sign is moot here, since already checked for < 0) and restore hidden bit:
			q.hi = TWOE52 + (q.hi & MASK_MANT);
			// When adding x < 1 to 1, only chance of carry into exponent field is if x = (1 - [tiny]) and gets rounded to 1:
			/***** <= 53 changed to < 53 here: ***********************/
			if(rshift < 53)	/* Hi part partially shifted into lo, lo part partially shifted off. */
			{
				offword = (q.lo << lshift) >> 63;	/* Round by adding the MSB of off-shifted low word into next-higher word. */
				q.lo = (q.lo >> rshift) + (q.hi << lshift);
				q.lo += offword;
				q.hi = (q.hi >> rshift) + (q.lo < offword);
				hi0 += q.hi;
			}
			else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into q.hi, lo part partially shifted off. */
			{
				offword = (q.lo << lshift) >> 63;
				q.lo = (q.lo >> rshift) + (q.hi << lshift) + offword;
			}
			else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into q.hi, lo part partially shifted off. */
			{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
				offword = q.lo >> 63;
				q.lo = q.hi + offword;
			}
			else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
			{
				rshift -= 64;
				lshift += 64;	/* lshift treated as unsigned, so put back into range. This makes q.hi behave like q.lo does in the above segments. */
				offword = (q.hi << lshift) >> 63;
				q.lo = (q.hi >> rshift) + offword;
			}
			else			/* Entire summand completely shifted off,  no need to round result. */
			{
				return QONE;
			}
		} else {	// If x >= 1, right-shift 1 to align with x and add
		// This is just the significand-add section of qfsum with the following argument value specializations:
		//		lo0 =        q.lo, hi0 =                                    q.hi .
		//		lo1 = QONE.lo = 0, hi1 = (QONE.hi & MASK_MANT) + TWOE52 = TWOE52
			exp0 = q.hi & ~MASK_MANT;
			hi0 = TWOE52 + (q.hi & MASK_MANT);
			hi1 = TWOE52;
			rshift = (q.hi>>52) - (QONE.hi>>52);	// Must shift before - to avoid borrow-into-exponent
			lshift = 64 - rshift;

			if(rshift == 0)		/* Operands perfectly aligned */
			{
				hi0 += hi1;
			}
			/***** <= 53 changed to < 53 here: ***********************/
			else if(rshift < 53)	/* Hi part partially shifted into lo, lo part partially shifted off. */
			{	// lshift in [12,63]
				// All off-shifted bits = 0 here:
				hi0 += (hi1 >> rshift);
			}
			else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
			{	// lshift in [1,11]
				q.lo += (hi1 << lshift);
				hi0 += (q.lo < (hi1 << lshift));
			}
			else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
			{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
				q.lo += hi1;
				hi0 += (q.lo < hi1);
			}
			else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
			{
				rshift -= 64;
				lshift += 64;	/* lshift treated as unsigned, so put back into range. This makes hi1 behave like 0 does in the above segments. */
				offword = (hi1 << lshift) >> 63;
				offword += (hi1 >> rshift);
				q.lo += offword;
				hi0 += (q.lo < offword);
			}
			else			/* Entire summand completely shifted off,  no need to round result. */
			{
				return q;
			}
		}
		/* If carry a one into 53-bit <53:53> of the high part, right-shift sum one place and add one to result exponent.
		Potentially need to do this step twice, if round bit from low word ripples all the way back into 53-bit, but we'll
		conveniently ignore that complication.
		*/
		rshift = hi0 >> 52;	/* rshift = 3, 2, 1 or 0 */
		if(rshift > 1)
		{
			q.hi = exp0 + (hi0 >> 1);			/* Right-shift high word one place, then simply add result to exponent field,
														exp. gets incremented automatically if 53-bit was set. */
			q.lo = (q.lo >> 1) + (hi0 << 63);
		}
		else if(rshift == 1)
		{
			q.hi = exp0 + (hi0 & MASK_MANT);	/* Mask off MSB of result significand prior to adding to exponent field. */
		}
		else	/* Both inputs were denormalized, result is also denormalized. Subtract one from exponent field prior to returning. */
		{
			q.hi = exp0 - TWOE52 + hi0;
		}
	}
#if QFDEBUG
	ASSERT(qfcmpeq(q, qfadd(x,QONE)), "qfinc fails!");
#endif
	return q;
}

struct qfloat qfdec(struct qfloat q)
{
	ASSERT(0, "qfdec not supported yet!");
	return qfsub(q, QONE);
}


/*
Top-level add and subtract routines seen by the caller - these examine the
signs of the inputs, send the proper combination of +-q1 and +-q2 to the
low-level qfsum and qfdif functions, and properly sign the result.
*/
struct qfloat qfadd	(struct qfloat q1, struct qfloat q2)
{
	struct qfloat q;
	uint32 sgn1 = (int64)q1.hi < 0, sgn2 = (int64)q2.hi < 0, both_signs = (sgn1<<1) + sgn2;

	if(both_signs == 3)		/* Both inputs negative */
	{
		q1.hi ^= MASK_SIGN;	/* Send -q1 and -q2 to low-level add routine... */
		q2.hi ^= MASK_SIGN;
		q = qfsum(q1, q2); q.hi ^= MASK_SIGN;	/* ...and return negated abs-value sum. */
	}
	else if(both_signs == 2)	/* q1 negative */
	{
		q1.hi ^= MASK_SIGN; /* Send q2 and -q1 to low-level sub routine. */
		q = qfdif(q2, q1);
	}
	else if(both_signs == 1)	/* q2 negative */
	{
		q2.hi ^= MASK_SIGN; /* Send q1 and -q2 to low-level sub routine. */
		q = qfdif(q1, q2);
	}
	else if(both_signs == 0)	/* Both inputs positive */
	{
		q = qfsum(q1, q2);
	}
	else
	{
		ASSERT(0,"ERROR: unrecognized sign combination in QFADD");
		q = QZRO; // silence warning
	}
#if QFDEBUG
	double qres = qfdbl(q), dres = (1-2.0*sgn1)*qfdbl(q1) + (1-2.0*sgn2)*qfdbl(q2);	// Must cast sgn1,2 to double prior to 1-...
	double esum = fabs(qres + dres);
	if(esum > 1e-15) {
		double rerr = fabs( (qres - dres)/esum );
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFADD!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFADD!\n");
	}
#endif

	return q;
}

struct qfloat qfsub	(struct qfloat q1, struct qfloat q2)
{
	struct qfloat q;
	uint32 sgn1 = (int64)q1.hi < 0, sgn2 = (int64)q2.hi < 0, both_signs = (sgn1<<1) + sgn2;

	if(both_signs == 3)		/* Both inputs negative */
	{
		q1.hi ^= MASK_SIGN;
		q2.hi ^= MASK_SIGN;
		q = qfdif(q2, q1);	/* Return -x - (-y) computed as (y - x), where x,y are absvals of inputs. */
	}
	else if(both_signs == 2)	/* q1 negative */
	{
		q1.hi ^= MASK_SIGN; /* Send q2 and -q1 to low-level add routine. */
		q = qfsum(q1, q2); q.hi ^= MASK_SIGN;	/* ...and return negated abs-value-sum. */
	}
	else if(both_signs == 1)	/* q2 negative */
	{
		q2.hi ^= MASK_SIGN; /* Send q1 and -q2 to low-level add routine. */
		q = qfsum(q1, q2);
	}
	else if(both_signs == 0)	/* Both inputs positive */
	{
		q = qfdif(q1, q2);
	}
	else
	{
		ASSERT(0,"ERROR: unrecognized sign combination in QFSUB");
		q = QZRO; // silence warning
	}
#if QFDEBUG
	double qres = qfdbl(q), dres = (1-2.0*sgn1)*qfdbl(q1) - (1-2.0*sgn2)*qfdbl(q2);	// Must cast sgn1,2 to double prior to 1-...
	double esum = fabs(qres + dres);
	if(esum > 1e-15) {
		double rerr = fabs( (qres - dres)/esum );
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFSUB!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFSUB!\n");
	}
#endif

	return q;
}

// Absolute-difference:
struct qfloat qfabsdiff(struct qfloat q1, struct qfloat q2)
{
	return qfabs(qfsub(q1,q2));
}

/*
Qfloat addition algorithm. Assumes both summands are nonnegative. A fair amount of code,
but should need only about 10 cycles in a nicely optimized machine implementation.

Example: in our qfloat format,

pi = 2^2 * 0.1100100100001111110110101010001000100001011010001100	0010001101001100010011000110011000101000101110000000110111000010,
 e = 2^2 * 0.1010110111111000010101000101100010100010101110110100	1010100110101010111111011100010101100010000000100111001111010100,

and the exact sum pi + e looks like

     2^2 * 1.0111011100001000001011101111101011000100001001000000	1100110011110111010010100010101110001010101110101000000110010110
   = 2^3 * 0.1011101110000100000101110111110101100010000100100000	0110011001111011101001010001010111000101010111010100000011001011.

Example 2: the exact sum pi + 1024*e looks like

     2^12 * 0.0000000000110010010000111111011010101000100010000101	1010001100001000110100110001001100011001100010100010111000000011
   + 2^12 * 0.1010110111111000010101000101100010100010101110110100	1010100110101010111111011100010101100010000000100111001111010100

   = 2^12 * 0.1010111000101010100110000100111101001011010000111010	0100110010110011110100001101100001111011100011001010000111010111.
*/
struct qfloat qfsum(struct qfloat q1, struct qfloat q2)
{
	struct qfloat return_val;
	struct qfloat *ptr0, *ptr1;
	uint32 rshift, lshift;
	uint64 exp0, exp1, hi0, hi1, lo0, lo1, offword;

	/* Make sure both inputs are nonnegative. */
	DBG_ASSERT(((int64)q1.hi >= 0 && (int64)q2.hi >= 0),"ERROR 160 in qfloat.c");

	/* Larger of the two operands gets index 0 in our local length-2 arrays: */
	if(qfcmple(q2, q1))
	{
		ptr0 = &q1;	ptr1 = &q2;
	}
	else
	{
		ptr0 = &q2;	ptr1 = &q1;
	}

	/* Separate upper parts of the 2 significands from the sign/exponent fields: */

	hi0  = ptr0->hi & MASK_MANT;
	hi1  = ptr1->hi & MASK_MANT;

	/* Sign = zero, so get exponent by simply subtracting off low bits. */
	exp0 = ptr0->hi - hi0;
	exp1 = ptr1->hi - hi1;

	lo0  = ptr0->lo;
	lo1  = ptr1->lo;

	/* If number is normalized, restore hidden bit. Otherwise, add one to the exponent field to get the true exponent. */
	if(exp0) {
		hi0  += TWOE52;
	} else {
		exp0 += TWOE52;
	}

	if(exp1) {
		hi1  += TWOE52;
	} else {
		exp1 += TWOE52;
	}

	/* ...then right-shift the smaller summand and add. */
	rshift = (exp0 - exp1) >> 52;
	lshift = 64 - rshift;

	if(rshift == 0)		/* Operands perfectly aligned */
	{
		lo0 += lo1;
		hi0 += hi1 + (lo0 < lo1);	/* carry? */
	}
	else if(rshift <= 53)	/* Hi part partially shifted into lo, lo part partially shifted off. */
	{
		offword = (lo1 << lshift) >> 63;	/* Round by adding the MSB of off-shifted low word into next-higher word. */
		lo1 = (lo1 >> rshift) + (hi1 << lshift);
		lo1 += offword;
		hi1 = (hi1 >> rshift) + (lo1 < offword);
		lo0 += lo1;
		hi0 += hi1 + (lo0 < lo1);
	}
	else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
	{
		offword = (lo1 << lshift) >> 63;
		lo1 = (lo1 >> rshift) + (hi1 << lshift) + offword;
		lo0 += lo1;
		hi0 += (lo0 < lo1);
	}
	else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
	{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
		offword = lo1 >> 63;
		lo1 = hi1 + offword;
		lo0 += lo1;
		hi0 += (lo0 < lo1);
	}
	else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
	{
		rshift -= 64;
		lshift += 64;	/* lshift treated as unsigned, so put back into range. This makes hi1 behave like lo1 does in the above segments. */
		offword = (hi1 << lshift) >> 63;
		lo1 = (hi1 >> rshift) + offword;
		lo0 += lo1;
		hi0 += (lo0 < lo1);
	}
	else			/* Entire summand completely shifted off,  no need to round result. */
	{
		return *ptr0;
	}

	/* If carry a one into 53-bit <53:53> of the high part, right-shift sum one place and add one to result exponent.
	Potentially need to do this step twice, if round bit from low word ripples all the way back into 53-bit, but we'll
	conveniently ignore that complication.
	*/
	rshift = hi0 >> 52;	/* rshift = 3, 2, 1 or 0 */
	if(rshift > 1)
	{
		return_val.hi = exp0 + (hi0 >> 1);			/* Right-shift high word one place, then simply add result to exponent field,
													exp. gets incremented automatically if 53-bit was set. */
		return_val.lo = (lo0 >> 1) + (hi0 << 63);
	}
	else if(rshift == 1)
	{
		return_val.hi = exp0 + (hi0 & MASK_MANT);	/* Mask off MSB of result significand prior to adding to exponent field. */
		return_val.lo = lo0;
	}
	else	/* Both inputs were denormalized, result is also denormalized. Subtract one from exponent field prior to returning. */
	{
		return_val.hi = exp0 - TWOE52 + hi0;
		return_val.lo = lo0;
	}

#if QFDEBUG
	double qres = qfdbl(return_val), dres = qfdbl(q1)+qfdbl(q2);
	double esum = fabs(qres + dres);
	if(esum > 1e-15) {
		double rerr = fabs( (qres - dres)/esum );
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFSUM!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFSUM!\n");
	}
#endif
	return(return_val);
}


/*
Qfloat subtraction algorithm. Assumes both inputs are nonnegative.

Example: in our qfloat format,

pi = 2^2 * 0.1100100100001111110110101010001000100001011010001100	0010001101001100010011000110011000101000101110000000110111000010,
 e = 2^2 * 0.1010110111111000010101000101100010100010101110110100	1010100110101010111111011100010101100010000000100111001111010100,

and the difference pi - e looks like (with a borrow from the low-part subtract)

	 2^2 * 0.0001101100010111100001100100100101111110101011010111   0111100110100001010011101010000011000110101101011001100111110010
   = 2^-1* 0.1101100010111100001100100100101111110101011010111011   1100110100001010011101010000011000110101101011001100111101110000
In hex qfloat format, this is 0x399D8BC324BF56BB  0xCD0A750635ACCF70 .
*/
struct qfloat qfdif(struct qfloat q1, struct qfloat q2)
{
	struct qfloat return_val;
	struct qfloat *ptr0, *ptr1;
	uint32 rshift, lshift, lshift_max;
	uint64 exp0, exp1, hi0, hi1, lo0, lo1, offword;

	/* Make sure both inputs are nonnegative. */
	DBG_ASSERT(((int64)q1.hi >= 0) && ((int64)q2.hi >= 0),"ERROR 170 in qfloat.c");

	/* Larger of the two operands gets index 0 in our local length-2 arrays: */
	if(qfcmple(q2, q1))
	{
		/* One more reason to make function arguments qfloat rather than qfloat * : if q1, q2 were ptrs
		and pointed to the same thing, ptr0 and ptr1 would also point to that one address, which would
		trigger the sign-flip conditional prior to returning, causing a return value of -0 here.
		*/
		ptr0 = &q1;	ptr1 = &q2;
	}
	else
	{
		ptr0 = &q2;	ptr1 = &q1;
	}

	/* Separate upper parts of the 2 significands from the sign/exponent fields: */
	// Right-col annotations are for corner case      q1 = 3FF0000000000000                0 -
	//                                                q2 = 3FEFFFFFFFFFFFFF FFFFFFFFFFFFFFFF
	// Which yielded the "inverted" result +NaN up though my debug of Sep 2012:
	hi0  = ptr0->hi & MASK_MANT;	//             0
	hi1  = ptr1->hi & MASK_MANT;	// FFFFFFFFFFFFF

	/* Sign assumed zero, so get exponent by simply subtracting off low bits. */
	exp0 = ptr0->hi - hi0;	// 3FF0000000000000
	exp1 = ptr1->hi - hi1;	// 3FE0000000000000

	lo0  = ptr0->lo;	//                0
	lo1  = ptr1->lo;	// FFFFFFFFFFFFFFFF

	/* If number is normalized, restore hidden bit.
	Otherwise, add one to the exponent field to get the true exponent. */
	if(exp0) {
		hi0  += TWOE52;	// 10000000000000
	} else {
		exp0 += TWOE52;
	}

	if(exp1) {
		hi1  += TWOE52;	// 1FFFFFFFFFFFFF
	} else {
		exp1 += TWOE52;
	}

	/* ...then right-shift the smaller summand and add. */

	rshift = (exp0 - exp1) >> 52;	//  1
	lshift = 64 - rshift;			// 63

	if(rshift == 0)		/* Operands perfectly aligned */
	{
		lo0 -= lo1;
		hi0 -= hi1 + (lo0 > ptr0->lo);		/* Had a borrow. */
		/* Since ptr0 >= ptr1, hi0 >= 0, no check for borrow-out needed here: */
	}
	else if(rshift <= 53)		/* Hi part partially shifted into lo, lo part partially shifted off. */
	{
		/* Round by adding the MSB of off-shifted bits into resulting right-justified low word, and shift in the low bits of hi1: */
		offword = (lo1 << lshift) >> 63;			// 1
		lo1 = (lo1 >> rshift) + (hi1 << lshift);	// 7FFFFFFFFFFFFFFF + 8000000000000000 = FFFFFFFFFFFFFFFF
		lo1 += offword;								// 0
		hi1 = (hi1 >> rshift) + (lo1 < offword);	/* Had a carry from the rounded bit being added into the right-shifted low part. */
													// (1FFFFFFFFFFFFF >> 1) + 1 = FFFFFFFFFFFFF + 1 = 10000000000000
		lo0 -= lo1;	// 0 - 0 = 0
		hi0 -= hi1 + (lo0 > ptr0->lo);				/* Had a borrow. */
													// 10000000000000 - 10000000000000 + 0 = 0 [before 9/10/12 bugfix gave 1]
	}
	else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
	{
		offword = (lo1 << lshift) >> 63;
		lo1 = (lo1 >> rshift) + (hi1 << lshift) + offword;

		lo0 -= lo1;
		hi0 -= (uint64)(lo0 > ptr0->lo);		/* Had a borrow. */
	}
	else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
	{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
		offword = lo1 >> 63;
		lo1 = hi1 + offword;

		lo0 -= lo1;
		hi0 -= (uint64)(lo0 > ptr0->lo);		/* Had a borrow. */
	}
	else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
	{
		rshift -= 64;
		lshift += 64;	/* lshift treated as unsigned, so put back into range. This makes hi1 behave like lo1 does in the above segments. */

		offword = (hi1 << lshift) >> 63;
		lo1 = (hi1 >> rshift) + offword;

		lo0 -= lo1;
		hi0 -= (uint64)(lo0 > ptr0->lo);	/* Had a borrow. */
	}
	else			/* Entire summand completely shifted off, no need to round result. */
	{
		if(ptr0 == &q2) {			// q1 - [rounds-to-0] = +q1
			ptr0->hi ^= MASK_SIGN;	// [rounds-to-0] - q2 = -q2
		}
		return *ptr0;
	}

	/* Now left-shift result so MSB is in bit <52:52> of high word. */
	if(hi0) {
		lshift = leadz64(hi0) - 11;	// [Before 9/10/12 fix took this branch, gave -11]
	} else {
		lshift = leadz64(lo0) + 53;	// 117
	}

	if(lshift && (exp0 > TWOE52))	/* Left-justify mantissa if subtract caused us to lose bits at the left end AND larger input was normalized. */
	{
		lshift_max = (exp0 - TWOE52) >> 52;	// (3FF0000000000000 - 0010000000000000)>>52 = (3FE0000000000000>>52 = 3FE = 1022_10
		/* If result underflows, left-shift mantissa just enough to make exp field = 1. */
		lshift = MIN(lshift, lshift_max);	// = lshift = 117
		rshift = 64 - lshift;				/* Don't have to worry about this being  = 64, since lshift > 0 here. */

		if(lshift < 53)	/* Hi part is nonzero, gets left-shifted <= 52 bits, upper (lshift) bits of lo shifted in to fill vacated bits. */
		{
			hi0 = (hi0 << lshift) + (lo0 >> rshift);
			lo0 = (lo0 << lshift);
			exp0 -= (uint64)lshift << 52;
		}
		else	/* Hi part is zero. Assuming lo part has lzlo lead zeros, right-shift lo (11-lzlo) bits and put that into hi. */
		{
		#if QFDEBUG
			printf("WARNING: catastrophic loss of precision in subtract:\n %16" PRIX64 " %16" PRIX64 " -\n %16" PRIX64 " %16" PRIX64 "\n", ptr0->hi, ptr0->lo, ptr1->hi, ptr1->lo);
		#endif
		//	return QZRO; *** Taking the easy way out breaks older already-tested stuff in qtest() ***
			if((int32)rshift > 0)	/* Lo part has > 53 SB, upper 53 of which get put into hi part. */
			{						/* That is, lo gets right-shifted <= 11 bits, result put into hi. Lshift now in [52, 63), i.e. in range. */
				hi0 = (lo0 >> rshift);
				lo0 = (lo0 << lshift);
				exp0 -= (uint64)lshift << 52;
			}
			else	// (int32)rshift = -53		/* Lo part has <= 53 SB, all of which get put into hi part. */
			{
				if(lshift < 117)	/* Lo part has at least one SB */
				{
					hi0 = (lo0 << -rshift);
					lo0 = (uint64)0;
					exp0 -= (uint64)lshift << 52;
				}
				else				/* Result = 0 */
				{
					exp0 = 	TWOE52;	/* Set exponent = 1, since we're going to subtract one prior to returning. */
				}
			}
		}
	}

	return_val.hi = exp0 - TWOE52 + hi0;	/* exponent field unchanged if 52-bit = 1, decremented by 1 if both = 0. */
	if(ptr0 == &q2) {
		return_val.hi ^= MASK_SIGN;	/* If flipped order of inputs, negate result. */
	}
	return_val.lo = lo0;
#if QFDEBUG
	double qres = qfdbl(return_val), dres = qfdbl(q1)-qfdbl(q2);
	double esum = fabs(qres + dres);
	if(esum > 1e-15) {
		double rerr = fabs( (qres - dres)/esum );
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFDIF!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFDIF!\n");
	}
#endif

	return(return_val);
}


/*
Newtonian iterative inverse scheme: given nonzero x, whose multiplicative
inverse we wish to find, we begin with an initial guess y and iterate as follows:

		y = y*(2 - x*y) .

Assuming the initial guess is reasonably close to 1/x, the number of significant
bits in y should roughly double on each iteration. If x is a qfloat and our initial
guess is the double-precision-accurate inverse of x, we need just 2 such iterations.
*/
struct qfloat qfinv(struct qfloat x)
{
	struct qfloat xinv, xyprod;
#ifdef X87_ASM
	long double ld;
  #if QFDEBUG
 	uint64 *ld_ptr;
  #endif
#else
	double xinv_dble;
#endif
#if QFDEBUG
	double qres, dres, esum, rerr;
#endif

	/* Make sure x is properly normalized. This also catches potential divides-by-zero. */
	if((x.hi & ~(MASK_SIGN + MASK_MANT)) == (uint64)0)
	{
		ASSERT(0,"ERROR: divide by denormalized input not supported.");
	}
#ifdef X87_ASM
	ld = qfldbl(x);
  #if QFDEBUG
 	ld_ptr = &ld;
  #endif
	asm ("fldt %1;"	/* push (log double)x onto FP stack*/
		 "fld1;"	/* push +1.0 onto FP stack, so st0 = 1.0, st1 = x */
		 "fdivp %%st(1);"	/* st0 = 1.0/x */
		 "fstpt %0"	/* (long double)x = 1.0/x */
		 : "+m"(ld) : "m"(ld));
	xinv = ldbl_to_q(ld);
#else
	/* Get double-precision approximation to 1/x: */
	xinv_dble = fisqrtest(qfdbl(x), 53);
	xinv_dble = finvest(qfdbl(x), 53);
	xinv = dbl_to_q(xinv_dble);
#endif

	/* Do 1 or 2 Newton iterations, depending on initial precision: */
	// Iter #1:
	xyprod = qfmul(x, xinv);
	xyprod = qfsub(QTWO, xyprod);
	xinv   = qfmul(xinv, xyprod);

#ifndef X87_ASM	// With higher-precision x87 init, only need one N-R iteration
	// Iter #2:
	xyprod = qfmul(x, xinv);
	xyprod = qfsub(QTWO, xyprod);
	xinv   = qfmul(xinv, xyprod);
#endif

#if QFDEBUG
	qres = qfdbl(xinv), dres = 1.0/qfdbl(x);
	esum = fabs(qres + dres);
	if(esum > 1e-15) {
		rerr = fabs( (qres - dres)/esum );
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFINV!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFINV!\n");
	}
#endif

	return xinv;
}


/*
Divide q1/q2, using iterative inverse and a multiply.
*/
struct qfloat qfdiv(struct qfloat q1, struct qfloat q2)
{
	struct qfloat qinv = qfmul(q1, qfinv(q2));

#if QFDEBUG
	double qres = qfdbl(qinv), dres = qfdbl(q1)/qfdbl(q2);
	double esum = fabs(qres + dres);
	if(esum > 1e-15) {
		double rerr = fabs( (qres - dres)/esum );
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFDIV!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFDIV!\n");
	}
#endif
	return qinv;
}

/*
Return a quotient p/q of 64-bit ints in qfloat form.
*/
struct qfloat qf_rational_quotient(int64 p, int64 q)
{
	return qfdiv(i64_to_q(p),i64_to_q(q));
}

/*
Newtonian iterative square root scheme: given nonnegative x, whose square root
we wish to find, we begin with an initial guess to 1/sqrt(x) and iterate as follows:

		y = y*(3 - x*y^2)/2 .

Assuming the initial guess is reasonably close, the iteration should converge to 1/sqrt(x),
with the number of significant bits in y roughly doubling on each iteration. If x is a qfloat and
our initial guess is the double-precision-accurate approximation to sqrt(x), we need just 2 such
iterations.
*/
struct qfloat qfsqrt(struct qfloat x)
{
	struct qfloat xisrt, xysqr;
#ifdef X87_ASM
	long double ld;
  #if QFDEBUG
 	uint64 *ld_ptr;
  #endif
#else
	double xisrt_dble;
#endif
#if QFDEBUG
	double qres, dres, esum, rerr;
#endif

	/* Make sure x is nonnegative. This also catches potential divides-by-zero. */
	ASSERT(!(x.hi >> 63),"ERROR: sqrt of a negative number not supported.");
	if(qfcmpeq(x, QZRO)) return QZRO;
#ifdef X87_ASM
	ld = qfldbl(x);
  #if QFDEBUG
 	ld_ptr = &ld;
  #endif
	asm ("fldt %1;"
		 "fst %%st(1);"
		 "fsqrt;"
		 "fdivp %%st(1);"
		 "fstpt %0" : "+m"(ld) : "m"(ld));
	xisrt = ldbl_to_q(ld);
#else
	/* Get double-precision approximation to 1/sqrt(x): */
	xisrt_dble = fisqrtest(qfdbl(x), 53);
//	xisrt_dble = (double)1.0/sqrt(qfdbl(x));
	xisrt = dbl_to_q(xisrt_dble);
#endif

	/* Do 1 or 2 Newton iterations, depending on initial precision: */
	// Iter #1:
	xysqr = qfmul(x, qfmul(xisrt, xisrt));
	xysqr = qfsub(QTHREE, xysqr);
	xisrt = qfmul_pow2(qfmul(xisrt, xysqr), -1);

#ifndef X87_ASM	// With higher-precision x87 init, only need one N-R iteration
	// Iter #2:
	xysqr = qfmul(x, qfmul(xisrt, xisrt));
	xysqr = qfsub(QTHREE, xysqr);
	xisrt = qfmul_pow2(qfmul(xisrt, xysqr), -1);
#endif

#if QFDEBUG
	qres = qfdbl(qisrt), dres = 1.0/sqrt(qfdbl(q));
	esum = fabs(qres + dres);
	if(esum > 1e-15) {
		rerr = fabs( (qres - dres)/esum );
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFSQRT!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFSQRT!\n");
	}
#endif
	/* Multiply 1/sqrt(x) by x to get sqrt(x). */
	return qfmul(x, xisrt);
}

// Some applications need the inverse sqrt, so can skip the final multiply-by-operand step:
struct qfloat qisqrt(struct qfloat q)
{
	double qisrt_dble;
	struct qfloat qisrt, xysqr;

	/* Make sure q is nonnegative. This also catches potential divides-by-zero. */
	if(q.hi >> 63)
	{
		ASSERT(0,"ERROR: sqrt of a negative number not supported.");
	}
	else if(qfcmpeq(q, QZRO))
	{
		return QZRO;
	}

	/* Get double-precision approximation to 1/sqrt(x): */
	qisrt_dble = fisqrtest(qfdbl(q), 53);
	qisrt = dbl_to_q(qisrt_dble);

	/* Do 2 Newton iterations: */

	xysqr = qfmul(q, qfmul(qisrt, qisrt));
	xysqr = qfsub(QTHREE, xysqr);
	qisrt = qfmul_pow2(qfmul(qisrt, xysqr), -1);

	xysqr = qfmul(q, qfmul(qisrt, qisrt));
	xysqr = qfsub(QTHREE, xysqr);
	qisrt = qfmul_pow2(qfmul(qisrt, xysqr), -1);

	return qisrt;
}

/***********************************************************************************************/
/************** Elementary transcendentals use AGM for speed wherever convenient ***************/
/***********************************************************************************************/

// AGM-iteration utility:
struct qfloat qfagm(struct qfloat x, struct qfloat y)
{
	int i;	uint64 idiff;
	struct qfloat a = x, b = y, tmp;
	// Limit #iterations to 20 (note we *cannot* do 1st few in lower-prec); allow low halves to differ slightly to avoid
	// excess work related from trying to purify away that last tiny bit of roundoff error:
#if YOU_DONT_BELEIV_ME_ABOUT_THE_PRECISION_THING
	double da = qfdbl(a),db = qfdbl(b),dtmp,rerr;
	for(i = 0; i < 10; ++i) {
		dtmp = da*db;
		da = 0.5*(da + db);
		db = sqrt(dtmp);
		rerr = (da - db)/(da + db);
		rerr = fabs( rerr );
		if( rerr < 1e-12 ) {
			break;
		}
	}
//	a = dbl_to_q(da);	b = dbl_to_q(db);
#endif
	for(i = 0; i < 20; ++i) {
		tmp = qfmul(a,b);
		a = qfmul_pow2(qfadd(a,b), -1);
		b = qfsqrt(tmp);
		if(a.hi == b.hi) {
		//	if(a.lo>>3) == (b.lo>>3)) break;	// More robust that !qfcmpeq(a,b), but still suffers from "long carry, small diff" errors
			idiff = (a.lo-b.lo);	// Unsigned difference...
			if(idiff>>63) {			// Use absolute diff for closeness check
				idiff = -idiff;
			}
			if(idiff < 8) break;
		}
	}
	ASSERT((i < 20), "Failure to converge in QFAGM!");
	return a;
}

/* Natural logarithm: */
/*
Algo 0:
First scale: find smallest k s.t. y = (2^k)*x >= 2^64, then use that for y large,

            Pi/2
log(y) = ---------- + O(1/y^2) .
         AGM(1,4/y)

Thus for y >= 2^64 our rel.err. is negligible relative to qfloat precision.

Example: Compute log(x) with x = 2^64 this way. using Pari:

a=1.0;x=2^64;b=4.0/x;b*x

t=a*b;a=(a+b)/2;a
b=sqrt(t);b
[repeat until converged]

0.50000000000000000010842021724855044340	4.6566128730773925781250000000000000000 E-10
0.25000000023283064370807973753052522170	1.5258789062500000001654361225106055350 E-5
0.12500762951094657185404069594587516388	0.0019531250009094947018788073563711304115
0.063480377255928033277959751651123147144	0.015625476840796292956088241362414450989
0.039552927048362163117023996506768799066	0.031494621202000750749563326632941336026
0.035523774125181456933293661569855067546	0.035294538597615013392143169346714794852
0.035409156361398235162718415458284931199	0.035408970854773150808257662478105142829
0.035409063608085692985488038968195037014	0.035409063607964210289944559157030849704
0.035409063608024951637716299062612943359	0.035409063608024951637716246964209975548
0.035409063608024951637716273013411459454	0.035409063608024951637716273013411459454

Pi/2/a
%27 = 44.361419555836499802702855773323300358
? Pi/2/a-log(x)
%28 = 1.1284745767893960076 E-36

qfloat: log(2^64) = {0x40462E42FEFA39EF,35793C7673007E5F}
*/
/*
Algo A: Define a small parameter eps (= 2^64 in our case), then use that

                     1              1
log(x) = Pi/2 * [ ---------- + ------------ ] + O(1/eps^2) .
                  AGM(1,eps)   AGM(1,eps*x)

Thus for eps < 1/2^64 our rel.err. is negligible relative to qfloat precision.
This runs into difficulty for large args (e.g. 2^64), so prefer Algo 0.
*/
/*
Algo B:
If x is regarded as fixed, given a function F(y) for which we have an effective algorithm, use Newton`s method
to compute inverse function G(y) at the given value x by applying N-R to the function F(y) - x:

	                F(y_n) - x
	y_{n+1} = y_n - ----------
	                  F`(y_n)

converges to the desired inverse function, y = G(x), under certain conditions.

Example: to compute G(x) = log(x) define F(y) = exp(y), use DP to compute initial seed y0 = log(x) to ~53 bits,
then do 1 or 2 of the above N-R inverse-function iterations in qfloat mode using F(y) = exp(y).
*/
struct qfloat qflog(struct qfloat x)
{
	struct qfloat y;
#ifdef X87_ASM
	long double ld;
  #if QFDEBUG
 	uint64 *ld_ptr;
  #endif
#else
	double lnx;
#endif
#if QFDEBUG
	double qres, dres, esum, rerr;
#endif
	struct qfloat expy;

	ASSERT(qfcmplt(QZRO,x), "Arg must be > 0 in QFLOG!");

#if 0	// Algo 0

	uint32 efield,k;
	// Find smallest k s.t. y = (2^k)*x >= 2^64:
	/* Extract 11-bit exponent field and add sign-extended power-of-2 exponent: */
	efield = ((x.hi >> 52) & MASK_EXP);	ASSERT(efield,"Denormalized numbers not currently supported in QFLOG");
	// 1.0 has efield = 0x3FF, so use that to compute k:
	k = 64 - (efield - 0x3FF);
	y = qfmul_pow2(x, (uint64)k);
	y = qfmul_pow2(qfinv(y), +2);	// 4/y
	y = qfagm(QONE,y);		// AGM(1.0, 4/y)
	y = qfdiv(QPIHALF,y);	// log(y) = (Pi/2) / AGM(1.0, 4/y)
	// log(x) = log(y/2^k) = log(y) - k*log2
	y = qfsub( y, qfmul( QLN2, i64_to_q((uint64)k ) ) );

#elif 0	// Algo A

	struct qfloat a,b;	//, qeps = qfmul_pow2(QONE,-64);	// qeps = 1/2^64
	a = QI64AGM;	// qfagm(QONE,qeps); a = qfinv(a);	// 1/AGM(1,eps); prefer precomputed version
	b = qfagm(QONE,qfmul_pow2(x,-64));	// AGM(1,eps*x);
	b = qfinv(b);
	y = qfmul(QPIHALF,qfsub(a,b));		// log(x)

#else	// Algo B

  #ifdef X87_ASM

	ld = qfldbl(x);

   // x87 supports scaled base-2 logarithm, y*l2(x); Use that ln(x) = l2(x)*ln(2), call fy2lx with scaling multiplier y = ln(2):
	asm ("fldln2;"
		 "fldt %1;"
		 "fyl2x;"
		 "fstpt %0" : "+m"(ld) : "m"(ld));
	y = ldbl_to_q(ld);
  #else
	lnx = log(qfdbl(x));
	y.hi = f64_to_u64(lnx); y.lo = 0ull;
  #endif

	// Iter #1:
	expy = qfexp(y);
	y = qfsub(y, qfdiv( qfsub(expy,x) , expy ) );	/*** If could get ~4 more bits in initial guess, could do just 1 iter ***/

  #ifndef X87_ASM	// With higher-precision x87 init, only need one N-R iteration
	// Iter #2:
	expy = qfexp(y);
	y = qfsub(y, qfdiv( qfsub(expy,x) , expy ) );
  #endif

#endif

#if QFDEBUG
	qres = qfdbl(y), dres = log(qfdbl(x));
	esum = fabs(qres + dres);
	if(esum > 1e-15) {
		rerr = fabs( (qres - dres)/esum );
	fprintf(stderr,"QFLOG: x = %20.10e, relerr = %20.10e\n", qfdbl(x), rerr);
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFLOG!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFLOG!\n");
	}
#endif
	return y;
}

/* Algo A uses simple argument-scaled Taylor series summation, followed by repeated squaring to un-scale the output.
This has 2 variants: The usual sum, and the "inside out" nested-product evaluation which requires #terms to be predomputed:

	exp(x) ~= 1 + x*(1 + (x/2)*(1 + (x/3)*(...)))

Algo B uses the basic complex exponential:

	exp(iy) = cos(y) + i*sin(y)

Setting y = -ix (x our real input argument), this is

	exp(x) = cos(ix) - i*sin(ix) .

How is this useful? Well, our sin() and cos() summations are faster than exp() because they involve half as many terms.
We can easily modify our real-argument cosine summation to instead return cos(ix):

	cos(ix) = 1 - (ix)^2/2 + (ix)^4/4! - (ix)^6/6! + ...

			= 1 + x^2/2 + x^4/4! + x^6/6! + ... = cosh(x), which is evaluatable just like cos(x) but with all + signs.

We could similarly modify the sine summation to get

	sin(ix) = (ix) - (ix)^3/3! + (ix)^5/5! - ...

			= i*[ x + x^3/3! + x^5/5! + ... ] , whence i*sin(ix) = -[ x + x^3/3! + x^5/5! + ... ] = -sinh(x),

but that doubles the work, negating the savings from the half-terms-summation property; instead use that the identity
(sin^2 + cos^2 = 1) holds equally well for real and complex argument (or in our case, pure-imaginary arguments) to
obtain sin(ix) from cos(ix) via

	sin(ix) = sqrt(1 - cos^2(ix)) = i*sqrt(cos^2(ix) - 1)   ==>   i*sin(ix) = -sqrt(cos^2(ix) - 1), where the arg of sqrt() is real.

Algo X [deprecated] uses Newtonian inverse-function iteration: to compute G(x) = exp(x) define F(y) = log(y) (computed via AGM),
then compute inverse function G(y) at the given value x by applying N-R to the function F(y) - x:

	                F(y_n) - x
	y_{n+1} = y_n - ---------- .
	                  F`(y_n)

To init, use DP to compute initial seed y0 = exp(x) to ~53 bits,
then do 1 or 2 of the above N-R inverse-function iterations in qfloat mode using F`(y) = 1/y.
*/
struct qfloat qfexp(struct qfloat x)
{
	 int32 pow2;
	double darg = qfdbl(x);
	struct qfloat xabs, y, sinh;

	// If arg < 0 compute 1/exp(-x) to avoid cancellation in the summation
	xabs = qfabs(x);
	if(xabs.hi < QEPS.hi) {
		return QONE;
	}
	// Normalize abs(argument) to be < 2 by manipulating exp-field. The 3FD threshold here is the largest
	// subtrahend which allows inputs of arbitrary size to be handled using at most 30 terms in the summation:
	pow2 = (xabs.hi >> 52) - 0x3fd;
	if(abs(pow2) > 9) {	// If arg > +- 512, check for over/underflow which occurs for |arg| ~> 700
		if(darg > 700) {
			fprintf(stderr,"QFEXP: xabs.hi = %16" PRIX64 ", pow2 = %u, darg = %10.10e\n",xabs.hi,pow2,darg);
			ASSERT(0,"expo overflow!");
		} else if(darg < -700) {
			return QZRO;	// expo underflow flushes to zero
		}
	}
	if(pow2 > 1) {
		xabs.hi -= ((uint64)(pow2-1) << 52);
	}

#if 0	// Algo X:
	#error *** need to add x87 high-prec exp() init ***
	double expx = exp(qfdbl(x));
	struct qfloat y, logy;
	y.hi = f64_to_u64(expx); y.lo = 0ull;
	// Iter #1:
	logy = qflog(y);
	y = qfsub(y, qfmul( qfsub(logy,x) , y ) );	/*** If could get ~4 more bits in initial guess, could do just 1 iter ***/
	// Iter #2:
	logy = qflog(y);
	y = qfsub(y, qfmul( qfsub(logy,x) , y ) );
	return y;

#elif 1	// Algo B:

	y = qfcosh(xabs);
	sinh = qfsqrt( qfsub(qfmul(y,y),QONE) );
	y = qfadd(y,sinh);

#elif 1	// Algo A:

	uint32 i;
	uint32 nterm,nterm_idx;
	const uint8 nterm_arr[64] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,6,6,6,7,7,8,8,9,9,10,10,11,12,13,15,17,19,22,25,30,30,30,30,30,30,30};
	struct qfloat curr_term, mult;
	/*
	Curve fit of #terms vs (0x3fd-capped) exp-field:
	Exp-range	#terms (excluding the 1 + ...)
	---------	------	hi_exp of range - 3C5 / #entries in lookup table
	3FE			30		57	/	1	[Use 6 of these to pad the byte-based lookup table to 64 entries]
	3FD			25		56	/	1
	3FC			22		55	/	1
	3FB			19		54	/	1
	3FA			17		53	/	1
	3F9			15		52	/	1
	3F8			13		51	/	1
	3F7			12		50	/	1
	3F6			11		49	/	1
	3F5-3F4		10		48	/	2
	3F3-3F2		9		46	/	2
	3F1-3F0		8		44	/	2
	3EF-3ED		7		42	/	2
	3EC-3E9		6		39	/	3
	3E8-3E3		5		35	/	4
	3E2-3D9		4		29	/	6
	3D8-3C6		3		19	/	10
	<= 3C5		2		0	/	19
	<= 38B		1
	*/
	i = (xabs.hi>>52);
	if(i < 0x3C5) {
		nterm = 1;
	} else {
		nterm_idx = i - 0x3C5;
	//	printf("Input exp-field = %3X, nterm_idx = %d; \n",i,nterm_idx);
		ASSERT(nterm_idx < 64, "nterm_idx ou of range!");
		nterm = nterm_arr[nterm_idx];
	}

	// 2 variants of summation:
  #define NESTED_PRODUCT	1
  #if NESTED_PRODUCT

	/* "inside-out" evaluation of N-term partial sum: */
	curr_term = qfadd(QONE, qfmul(xabs, QNINV[nterm]));	// Innermost term: 1 + x/N
	for(i = nterm-1; i > 1; i--)
	{
		mult      = qfmul(xabs, QNINV[i]);
		curr_term = qfinc(qfmul(curr_term, mult));	// The nested form allows us to replace generic ADD with faster INCrement
	}
	y = qfadd(QONE, qfmul(xabs, curr_term));

  #else
	/* Initialize partial sum to 1 + x... */
	curr_term = xabs;
	y = qfadd(QONE, curr_term);
	/* get next term of series by multiplying current one by x/i */
	for(i = 2; i <= nterm; i++)
	{
		mult      = qfmul(xabs, QNINV[i]);
		curr_term = qfmul(curr_term, mult);
		y         = qfadd(curr_term, y);
	}
   #if QFDEBUG
	i = (y.hi - curr_term.hi)>>52;
	printf("Input exp-field = %3X, Used %u-term summation, exp_diff = %d; \n",xabs.hi>>52,nterm,i);
   #endif

  #endif

#endif

	/*** Algos A & B share this: ***/
	// If arg < 0 compute 1/exp(-x) to avoid cancellation in the summation:
	if((int64)x.hi < 0) {
		y = qfinv(y);
	}
	// If normalized exponent, recover original via successive squarings:
	while(pow2 > 1) {
		--pow2;		y = qfmul(y,y);
	}

#if QFDEBUG
	double qres = qfdbl(y), dres = exp(qfdbl(x));
	double esum = fabs(qres + dres);
	if(esum > 1e-15) {
		double rerr = fabs( (qres - dres)/esum );
	fprintf(stderr,"QFEXP: x = %20.10e, relerr = %20.10e\n", qfdbl(x), rerr);
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFEXP!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFEXP!\n");
	}
#endif
	return y;
}

/* Arctangent: Use Newtonian inverse-function iteration: to compute G(x) = atan(x) define F(y) = tan(y),
then compute inverse function G(y) at the given value x by applying N-R to the function F(y) - x:

	                F(y_n) - x
	y_{n+1} = y_n - ---------- .
	                  F`(y_n)

To init, use DP to compute initial seed y0 = tan(x) to ~53 bits [64 bits ia inline ASM and long double on x87]
then do 1 or 2 of the above N-R inverse-function iterations in qfloat mode using F`(y) = 1/cos^2(y).
*/
struct qfloat qfatan(struct qfloat x)
{
	struct qfloat y, cos2y, tany;
#ifdef X87_ASM
	long double ld = qfldbl(x);
  #if QFDEBUG
 	uint64 *ld_ptr = &ld;
  #endif
	asm ("fldt %1;"
		 "fld1;"
		 "fpatan;"
		 "fstpt %0" : "+m"(ld) : "m"(ld));
	y = ldbl_to_q(ld);
#else
	double darct = atan(qfdbl(x));
	y.hi = f64_to_u64(darct); y.lo = 0ull;
#endif

	// Iter #1:
	cos2y = y;	// qftan_and_cos() Overwrites arg with cos(y), so first copy y into cos2y, since we still need original y
	tany = qftan_and_cos(&cos2y);
	cos2y = qfmul(cos2y,cos2y);	// cos^2(y)
	y = qfsub(y, qfmul( qfsub(tany,x) , cos2y ) );

#ifndef X87_ASM	// With higher-precision x87 init, only need one N-R iteration
	// Iter #2:
	tany = qftan(y);
	cos2y = qfcos(y);	cos2y = qfmul(cos2y,cos2y);	// cos^2(y)
	y = qfsub(y, qfmul( qfsub(tany,x) , cos2y ) );
#endif

#if QFDEBUG
	double qres = qfdbl(y), dres = atan(qfdbl(x));
	double esum = fabs(qres + dres);
	if(esum > 1e-15) {
		double rerr = fabs( (qres - dres)/esum );
	printf("QFATAN: x = %20.10e, relerr = %20.10e\n", qfdbl(x), rerr);
		if( rerr > 1e-12 ) {
			WARN(HERE, "High Error Level in QFATAN!\n", "", 0);
		}
		ASSERT(rerr < 1e-2, "Fatal ROE in QFATAN!\n");
	}
#endif
	return y;
}

/* Factorial: */
struct qfloat qffact(uint32 n)
{
	int i;
	struct qfloat nf = QONE;
	if(n > 170) {
		return QSNAN;	// 170! ~= 1e307, so 171! exceeds IEEE64 exp-field, also used by qfloat type
	}
	for(i = 2; i <= n; i++) {
		nf = qfmul(nf, i64_to_q((uint64)i));
	}
	return nf;
}

/*
Tangent function. We use sin/cos since the Taylor series for tan(x) converges vastly more slowly than either of those,
due to the poles where cos(x) = 0 [pari/GP result below]:

default(seriesprecision, 30)
tan(x) = x + 1/3*x^3 + 2/3.5*x^5 + 17/3^2.5.7*x^7 + 2.31/3^4.5.7*x^9 + 2.691/3^4.5^2.7.11*x^11 + 2^2.43.127/3^5.5^2.7.11.13*x^13
+ 257.3617/3^6.5^3.7^2.11.13*x^15 + 2.73.43867/3^6.5^3.7^2.11.13.17*x^17 + 2.31.41.283.617/3^8.5^3.7^2.11.13.17.19*x^19
+ 2^2.89.131.593.683/3^9.5^4.7^3.11.13.17.19*x^21 + 2.103.241.2294797/3^9.5^4.7^3.11^2.13.17.19.23*x^23
+ 2^2.2731.8191.657931/3^10.5^6.7^3.11^2.13.17.19.23*x^25 + 2^2.43.113.127.9349.362903/3^13.5^6.7^3.11^2.13^2.17.19.23*x^27
+ 2^3.151.331.1721.1001259881/3^13.5^6.7^4.11^2.13^2.17.19.23.29*x^29 + O(x^31)

Compare the coefficients to the 1/n! (n odd) of the sine series;

n	n!/tan_n
--	---------
1	1
3	2
5	16
7	272
9	7936
11	353792
13	22368256
15	1903757312
17	209865342976
19	29088885112832
21	4951498053124096
23	1015423886506852352
25	246921480190207983616
27	70251601603943959887872
29	23119184187809597841473536

The tan(x) coefficients thus decay woefully more slowly than the inverse-factorial of sin(x).
*/
struct qfloat qftan(struct qfloat x)
{
#if 0
	return qfdiv( qfsin(x), qfcos(x) );
#else
	struct qfloat sin = qfsin(x), sec = qisqrt(qfsub(QONE, qfmul(sin,sin)));
	return qfmul( sin, sec );
#endif
}

// Return tan(x) as function result, overwrite x with sin(x)
struct qfloat qftan_and_sin(struct qfloat *x)
{
	struct qfloat sin = qfsin(*x), sec = qisqrt(qfsub(QONE, qfmul(sin,sin)));
	*x = sin;
	return qfmul( sin, sec );
}

// Return tan(x) as function result, overwrite x with cos(x)
struct qfloat qftan_and_cos(struct qfloat *x)
{
	struct qfloat cos = qfcos(*x), sin = qfsqrt(qfsub(QONE, qfmul(cos,cos)));
	*x = cos;
	return qfdiv( sin, cos );
}

struct qfloat qfcot(struct qfloat q)
{
#if 0
	return qfdiv( qfcos(q), qfsin(q) );
#else
	struct qfloat cos = qfcos(q), csc = qisqrt(qfsub(QONE, qfmul(cos,cos)));
	return qfmul( cos, csc );
#endif
}

/* Top-level sine and cosine routines seen by the caller - these examine the
   sign and magnitude of the input and map that to the proper call to either
   +-sin(arg) or +-cos(arg), where arg is in [0, pi/2).

   Use the following symmetries, all based on the angle-sum identities

	cos(x+a) = cos(x)*cos(a) - sin(x)*sin(a);	sin(x+a) = sin(x)*cos(a) + cos(x)*sin(a);

   1) exp(I*(x + 2*k*pi)) = exp(I*x);						(full-rotation symmetry, a = 2*k*pi)
   2) cos(x + pi  ) = -cos(x);      sin(x + pi  ) = -sin(x);	(half-rotation symmetry, a = pi  )
   3) cos(x + pi/2) = -sin(x);      sin(x + pi/2) = +cos(x);	(symmetry about y-axis,  a = pi/2)

	Let i = int(theta/(pi/2) (mod 4) = 0-3 be the quadrant index, and gamma := theta - i*pi/2 in [0, pi/2). Then:

		i	cos(theta)		sin(theta)
		-	-----------		-----------
		0	+cos(gamma)		+sin(gamma)
		1	-sin(gamma)		+cos(gamma)
		2	-cos(gamma)		-sin(gamma)
		3	+sin(gamma)		-cos(gamma)

	We can even go one better by using the fact that all roots can be mapped to data in the first octant.

	Let i = int(theta/(pi/4) (mod 8) = 0-7 be the octant index, and gamma := theta - i*pi/4 in [0, pi/4). Then:

		i	cos(theta)			sin(theta)
		-	------------------	------------------
		0	+cos(       gamma)	+sin(       gamma)
		1	+sin(  pi/2-gamma)	+cos(  pi/2-gamma)
		2	-sin(       gamma)	+cos(       gamma)
		3	-cos(  pi  -gamma)	+sin(  pi  -gamma)
		4	-cos(       gamma)	-sin(       gamma)
		5	-sin(3*pi/2-gamma)	-cos(3*pi/2-gamma)
		6	+sin(       gamma)	-cos(       gamma)
		7	+cos(2*pi  -gamma)	-sin(2*pi  -gamma)

	We implement the latter by starting with the quadrant-based mapping, then mapping the resulting
	args in [0, pi/2) to [0, pi/4) in the low-level qfcs1 and qfsn1 routines using the symmetries

		cos(pi/2 - x) = sin(x),	    sin(pi/2 - x) = cos(x).
*/
struct qfloat qfcos	(struct qfloat q)
{
	struct qfloat qt;
	uint128 i;
	uint32 quad;
	qt = qfabs(q);				// Need deal only with nonnegative numbers.
	qt = qfmul(qt, QIPIHLF);	// Get q/(pi/2)...
	i  = qfint(qt);				// ...And truncate that to the next-smaller integer (128-bit int in this case).
	// For the quadrant, we only need the result modulo 4:
	quad = i.d0 & (uint64)3;
	ASSERT(!i.d1 && (int64)i.d0 >= 0,"QFCOS: quadrant error");
	qt = i64_to_q((int64)i.d0);
	// The calling argument is q mod pi/2:
	q = qfsub(q, qfmul(qt, QPIHALF));
	if(quad == 0)
		q = qfcs1(q);
	else if(quad == 1)
		q = qfneg(qfsn1(q));
	else if(quad == 2)
		q = qfneg(qfcs1(q));
	else if(quad == 3)
		q = qfsn1(q);
	return q;
}

struct qfloat qfsin	(struct qfloat q)
{
	struct qfloat qt;
	uint128 i;
	uint32 quad, sign = qfcmplt(q, QZRO);
	qt = qfabs(q);				// Need deal only with nonnegative numbers.
	qt = qfmul(qt, QIPIHLF);	// Get q/(pi/2)...
	i  = qfint(qt);				// ...And truncate that to the next-smaller integer (128-bit int in this case).
	// For the quadrant, we only need the result modulo 4:
	quad = i.d0 & (uint64)3;
	ASSERT(!i.d1 && (int64)i.d0 >= 0,"QFSIN: quadrant error");
	qt = i64_to_q((int64)i.d0);
	// The calling argument is q mod pi/2:
	q = qfsub(q, qfmul(qt, QPIHALF));
	if(quad == 0)
		q = qfsn1(q);
	else if(quad == 1)
		q = qfcs1(q);
	else if(quad == 2)
		q = qfneg(qfsn1(q));
	else if(quad == 3)
		q = qfneg(qfcs1(q));

	if(sign)
		return qfneg(q);
	else
		return q;
}

/***********************************************************************************/
/* low-level trig routines, whose args are presumed normalized to lie in [0,Pi/2). */
/***********************************************************************************/

/*
Cosine function. Algorithm is not sophisticated - given qfloat x, use the Taylor series about 0.
Requires an argument in [0, pi/2), and then maps that to [0, pi/4) using the symmetry

		cos(pi/2 - x) = sin(x);

Actually, the Taylor series itself doesn't require any limitations whatsoever on the input
argument, but for the sake of fast convergence and high accuracy, we desire smallish arguments.
Since this routine is designed strictly to be called from the high-level qfcos routine, which
maps an arbitrary qfloat argument into [0, pi/2), the restriction on inputs is no problem.

Jan 2002: Add fixed-order polynomial exoansions based on Chebyshev-basis approximation to cos,sin(x)
for x in [-pi/4,pi/4], using first 24 T_j(x). Here are the resulting max abs-error values relative to
the existing Taylor-series-based expansion:

maxerr =     1.811577853939244e-33 at x =    0.001315541923691
maxerr =     3.129636159629258e-34 at x =    0.753856573155407

*/
// Jan 2020: This is ~10% faster than Taylor-series approach on my CoreDuo MacBook
// but slightly less accurate, so leave the Taylor series way as the default:
#undef USE_CHEB_EXPANSION
#define USE_CHEB_EXPANSION	0	// 0 = Taylor series, > 0 = Chebyshev-derived near-minimax poly, 1 = branchless, 2 = branched
#if (USE_CHEB_EXPANSION < 0) || (USE_CHEB_EXPANSION > 2)
	#error Allowable values of USE_CHEB_EXPANSION = 0,1,2!
#endif
struct qfloat qfcs1(struct qfloat q)
{
#if (USE_CHEB_EXPANSION == 1)	// Branchless algorithm:
	const struct qfloat g[2] = {QZRO,QPIHALF}, qeps = {0x3930000000000000ull, 0x0000000000000000ull};;
	struct qfloat q2, sum;
#elif (USE_CHEB_EXPANSION == 2)
	int i;
	struct qfloat qsqr, sum;
	const struct qfloat coeff[12] = {
		{0x3FEFFFFFFFFFFFFFull,0xFFFFFFFFFFFFFED7ull},
		{0xBFDFFFFFFFFFFFFFull,0xFFFFFFFFFFFD69B3ull},
		{0x3FA5555555555555ull,0x55555555502458ABull},
		{0xBF56C16C16C16C16ull,0xC16C16AFDFC52F06ull},
		{0x3EFA01A01A01A01Aull,0x019FDA92EFCFD33Cull},
		{0xBE927E4FB7789F5Cull,0x7264FD63BD813E8Bull},
		{0x3E21EED8EFF8D896ull,0x33A732C1DF58EB2Eull},
		{0xBDA93974A8C07706ull,0x64DDC3EB89ACE555ull},
		{0x3D2AE7F3E725CD67ull,0xCF3B27A90337C664ull},
		{0xBCA682784CBB96C1ull,0xAFC0F5C914D95396ull},
		{0x3C1E53FB4C09439Aull,0x5F593CD59E826B71ull},
		{0xBB90B0F048BCA62Full,0x4FF38D712E5D4A81ull}
	};
#else
	int i,j;
	uint32 e_sum, e_new;
	struct qfloat qsqr, sum, curr_term, q1, mult;
	const int qsz = sizeof(struct qfloat);
	static struct qfloat*denoms = 0x0;
	static int first_entry = TRUE;
	if(first_entry) {
		first_entry = FALSE;
		denoms = (struct qfloat *)malloc(20*qsz);	ASSERT(denoms != NULL, "alloc failed!");
		for(i = 4; i < 38; i += 4)	// Must limit largest index into QNINV[] to 40, hence (i+2) < 40
		{
			j = (i>>1)-1;
			denoms[j-1] = qfmul(QNINV[i-1], QNINV[i  ]);
			denoms[j  ] = qfmul(QNINV[i+1], QNINV[i+2]);
		}
	}
#endif
	/* Make sure argument is in range... */
	DBG_ASSERT((qfcmple(qfneg(QEPS), q) && qfcmplt(q, qfadd(QPIHALF, QEPS))), "ERROR 200 in qfloat.c");

#if (USE_CHEB_EXPANSION == 1)	// Branchless algorithm:

	// q < pi/4: compute cos(q); q >= pi/4: compute sin(pi/2 - q). Break this into the following steps:
	int cos_or_sin = qfcmpge(q,QPI4TH);	// [0] q < pi/4: cos_or_sin = 0; q >= pi/4: cos_or_sin = 1
	q2 = qfsub(q,zero_or_pi2[cos_or_sin]);	// [1] q < pi/4: q2 = q - 0; q >= pi/4: q2 = q - pi/2
	q2.hi ^= (uint64)cos_or_sin << 63;		// [2] q < pi/4: q2 = q - 0; q >= pi/4: q2 = pi/2 - q
	sum = qfcos_or_sin1(q2,cos_or_sin);	// [3] q < pi/4: s2 = cos(q); q >= pi/4: s2 = sin(pi/2 - q).

#else	// Legacy algorithm w/branch

	if(qfcmpeq(QZRO, q)) {
		return QONE;
	} else if(qfcmpgt(q,QPI4TH)) {		// Argument > pi/4; return sin(pi/2 - q). ***NOTE:*** Use > here and >= in the corresponding
		return qfsn1(qfsub(QPIHALF, q));// qfsn1 code; if use >= for both, execution will bounce back and forth in an infinite recursive
	}									// loop (==> "EXC_BAD_ACCESS, Could not access memory") when q == QPI4TH.

  #if !USE_CHEB_EXPANSION
	/* Initialize partial sum to 1 - x^2/2... */
	qsqr= qfmul(q, q);
	curr_term = qfmul_pow2(qsqr, -1);
	sum = qfsub(QONE, curr_term);

	/* get next term of series by multiplying current one by x^2/(i*(i-1)) */
	for(i = 4; i < 38; i += 4)	// Must limit largest index into QNINV[] to 40, hence (i+2) < 40
	{
		j = (i>>1)-1;
		mult      = qfmul(qsqr, denoms[j-1]);
		curr_term = qfmul(curr_term, mult);
		q1        = curr_term;

		mult      = qfmul(qsqr, denoms[j  ]);
		curr_term = qfmul(curr_term, mult);
		q1        = qfsub(q1, curr_term);	/* First calculate the difference of the current 2 terms... */
		sum       = qfadd(sum, q1);			/* ...and only then update the partial sum. */

		/* If the current term differs from the sum by more than a factor of 2^110, stop. */
		e_sum = (uint32)((sum      .hi & ~MASK_SIGN) >> 52);
		e_new = (uint32)((curr_term.hi & ~MASK_SIGN) >> 52);
		if((int32)(e_sum - e_new) > 115)
		{
		#if QFDEBUG
			printf("End cos(x) summation after %d terms.\n", i/2);
		#endif
			break;
		}
	}
	ASSERT(((int32)(e_sum - e_new) > 115), "Unconverged cos(x) summation!");

  #elif (USE_CHEB_EXPANSION == 2)

	/* Compute first-12-even-bases Chebyshev approximation,
		c0 + c2.x^2 + c4.x^4 + c6.x^6 + c8.x^8 + c10.x^10 + c12.x^12 + c14.x^14 + c16.x^16 + c18.x^18 + c20.x^20 + c22.x^22,
	with the multiply-accumulation structured like so, preserving maximal accuracy by starting at the small (right) end and working leftward:
		c0 + x^2.(c2 + x^2.(c4 + x^2.(c6 + x^2.(c8 + x^2.(c10 + x^2.(c12 + x^2.(c14 + x^2.(c16 + x^2.(c18 + x^2.(c20 + x^2.c22))))))))))
	*/
	// Initialize accumulator to c20 + c22.x^2:
	qsqr= qfmul(q, q);
	sum = qfmul(qsqr,coeff[11]);
	sum = qfadd(coeff[10], sum);
	for(i = 9; i >= 0; i--)	{
		sum = qfmul(qsqr,sum);
		sum = qfadd(coeff[i], sum);
	}

  #endif
#endif
	return sum;
}

/*
Sine function. Algorithm is not sophisticated - given qfloat x, use the Taylor series about 0.
Requires an argument in [0, pi/2), and then maps that to [0, pi/4) using the symmetry

	    sin(pi/2 - x) = cos(x).

Actually, the Taylor series itself doesn't require any limitations whatsoever on the input
argument, but for the sake of fast convergence and high accuracy, we desire smallish arguments.
Since this routine is designed strictly to be called from the high-level qfsin routine, which
maps an arbitrary qfloat argument into [0, pi/2), the restriction on inputs is no problem.
*/
struct qfloat qfsn1(struct qfloat q)
{
#if (USE_CHEB_EXPANSION == 1)	// Branchless algorithm:
	const struct qfloat zero_or_pi2[2] = {QZRO,QPIHALF}, qeps = {0x3930000000000000ull, 0x0000000000000000ull};;
	struct qfloat q2, sum;
#elif (USE_CHEB_EXPANSION == 2)
	int i;
	struct qfloat qsqr, sum;
	const struct qfloat coeff[12] = {
		{0x3FF0000000000000ull,0x000000000000019Full},
		{0xBFC5555555555555ull,0x55555555555D3720ull},
		{0x3F81111111111111ull,0x11111111274A67C2ull},
		{0xBF2A01A01A01A01Aull,0x01A01A7441697770ull},
		{0x3EC71DE3A556C733ull,0x8FAC077F9BABDB73ull},
		{0xBE5AE64567F544E3ull,0x944939AFB8657D0Cull},
		{0x3DE6124613A86D13ull,0x19936FC2A332D1C8ull},
		{0xBD6AE7F3E733D377ull,0x7C466FD62F4FAF8Cull},
		{0x3CE952C770628C4Dull,0x785AD699F7C50A9Eull},
		{0xBC62F49B7DC66453ull,0x343E4FE1AC28ECB9ull},
		{0x3BD71BCED110E8CEull,0x0B2F6B7014428EB2ull},
		{0xBB4773EEB2D02E27ull,0x8A69E1E69474DF98ull}
	};
#else
	int i,j;
	uint32 e_sum, e_new;
	struct qfloat qsqr, sum, curr_term, q1, mult;
	const int qsz = sizeof(struct qfloat);
	static struct qfloat*denoms = 0x0;
	static int first_entry = TRUE;
	if(first_entry) {
		first_entry = FALSE;
		denoms = (struct qfloat *)malloc(20*qsz);	ASSERT(denoms != NULL, "alloc failed!");
		for(i = 3; i < 38; i += 4)	// Must limit largest index into QNINV[] to 40, hence (i+2) < 40
		{
			j = (i>>1);
			denoms[j-1] = qfmul(QNINV[i-1], QNINV[i  ]);
			denoms[j  ] = qfmul(QNINV[i+1], QNINV[i+2]);
		}
	}
#endif
	/* Make sure argument is in range... */
	DBG_ASSERT((qfcmple(qfneg(QEPS), q) && qfcmplt(q, qfadd(QPIHALF, QEPS))), "ERROR 210 in qfloat.c");

#if (USE_CHEB_EXPANSION == 1)	// Branchless algorithm:

	// q < pi/4: compute sin(q); q >= pi/4: compute cos(pi/2 - q). Break this into the following steps:
	int sin_or_cos = qfcmpge(q,QPI4TH);	// [0] q < pi/4: sin_or_cos = 0; q >= pi/4: sin_or_cos = 1
	q2 = qfsub(q,zero_or_pi2[sin_or_cos]);	// [1] q < pi/4: q2 = q - 0; q >= pi/4: q2 = q - pi/2
	q2.hi ^= (uint64)sin_or_cos << 63;		// [2] q < pi/4: q2 = q - 0; q >= pi/4: q2 = pi/2 - q
	sum = qfcos_or_sin1(q2,1-sin_or_cos);	// [3] q < pi/4: s2 = sin(q); q >= pi/4: s2 = cos(pi/2 - q).

#else	// Legacy algorithm w/branch

	if(qfcmpeq(QZRO, q)) {
		return QZRO;
	} else if(qfcmpge(q,QPI4TH)) {		// Argument >= pi/4; return cos(pi/2 - q). ***NOTE:*** Use >= here and > in the corresponding
		return qfcs1(qfsub(QPIHALF, q));// qfcs1 code; if use >= for both, execution will bounce back and forth in an infinite recursive
	}									// loop (==> "EXC_BAD_ACCESS, Could not access memory") when q == QPI4TH.

  #if !USE_CHEB_EXPANSION
	/* Initialize partial sum to x... */
	qsqr= qfmul(q, q);
	curr_term = sum = q;

	/* get next term of series by multiplying current one by x^2/(i*(i-1)) */
	for(i = 3; i < 38; i += 4)	// Must limit largest index into QNINV[] to 40, hence (i+2) < 40
	{
		j = (i>>1);
		mult      = qfmul(qsqr, denoms[j-1]);
		curr_term = qfmul(curr_term, mult);
		q1        = curr_term;

		mult      = qfmul(qsqr, denoms[j  ]);
		curr_term = qfmul(curr_term, mult);
		q1        = qfsub(q1, curr_term);	/* First calculate the difference of the current 2 terms... */
		sum       = qfsub(sum, q1);			/* ...and only then update the partial sum. */

		/* If the current term differs from the sum by more than a factor of 2^110, stop. */
		e_sum = (uint32)((sum      .hi & ~MASK_SIGN) >> 52);
		e_new = (uint32)((curr_term.hi & ~MASK_SIGN) >> 52);
		if((int32)(e_sum - e_new) > 115)
		{
		#if QFDEBUG
			printf("End sin(x) summation after %d terms.\n", i/2);
		#endif
			break;
		}
	}
	ASSERT(((int32)(e_sum - e_new) > 115), "Unconverged sin(x) summation!");

  #elif (USE_CHEB_EXPANSION == 2)

	// Compute first-12-odd-bases Chebyshev approximation;
	// Initialize accumulator to c21 + c23.x^2:
	qsqr= qfmul(q, q);
	sum = qfmul(qsqr,coeff[11]);
	sum = qfadd(coeff[10], sum);
	for(i = 9; i >= 0; i--)	{
		sum = qfmul(qsqr,sum);
		sum = qfadd(coeff[i], sum);
	}
	sum = qfmul(q,sum);	// Need a final mul-by-x in the sin(x) case

  #endif
#endif
	return sum;
}

struct qfloat qfcos_or_sin1(struct qfloat q, int cos_or_sin)
{
	int i,i0;
	struct qfloat qsqr, sum;
	const struct qfloat coeff[24] = {	// Now interleave the even (cosine) and odd (sine) Chebyshev coefficients:
		{0x3FEFFFFFFFFFFFFFull,0xFFFFFFFFFFFFFED7ull},
		{0x3FF0000000000000ull,0x000000000000019Full},
		{0xBFDFFFFFFFFFFFFFull,0xFFFFFFFFFFFD69B3ull},
		{0xBFC5555555555555ull,0x55555555555D3720ull},
		{0x3FA5555555555555ull,0x55555555502458ABull},
		{0x3F81111111111111ull,0x11111111274A67C2ull},
		{0xBF56C16C16C16C16ull,0xC16C16AFDFC52F06ull},
		{0xBF2A01A01A01A01Aull,0x01A01A7441697770ull},
		{0x3EFA01A01A01A01Aull,0x019FDA92EFCFD33Cull},
		{0x3EC71DE3A556C733ull,0x8FAC077F9BABDB73ull},
		{0xBE927E4FB7789F5Cull,0x7264FD63BD813E8Bull},
		{0xBE5AE64567F544E3ull,0x944939AFB8657D0Cull},
		{0x3E21EED8EFF8D896ull,0x33A732C1DF58EB2Eull},
		{0x3DE6124613A86D13ull,0x19936FC2A332D1C8ull},
		{0xBDA93974A8C07706ull,0x64DDC3EB89ACE555ull},
		{0xBD6AE7F3E733D377ull,0x7C466FD62F4FAF8Cull},
		{0x3D2AE7F3E725CD67ull,0xCF3B27A90337C664ull},
		{0x3CE952C770628C4Dull,0x785AD699F7C50A9Eull},
		{0xBCA682784CBB96C1ull,0xAFC0F5C914D95396ull},
		{0xBC62F49B7DC66453ull,0x343E4FE1AC28ECB9ull},
		{0x3C1E53FB4C09439Aull,0x5F593CD59E826B71ull},
		{0x3BD71BCED110E8CEull,0x0B2F6B7014428EB2ull},
		{0xBB90B0F048BCA62Full,0x4FF38D712E5D4A81ull},
		{0xBB4773EEB2D02E27ull,0x8A69E1E69474DF98ull}
	};
	// GCC/Clang would not allow me to set the first 2 elts = QONE,QZRO: "error: initializer element is not a compile-time constant":
	static struct qfloat one_or_q[2] = {
		{0x3FF0000000000000ull, 0x0000000000000000ull},	// = QONE
		{0x0000000000000000ull, 0x0000000000000000ull}	// = QZRO, used here as a compile-time placeholder
	};
	// Now insert the runtie value of the angular arg in the 2nd slot:
	one_or_q[1] = q;
	// Initialize accumulator to c20 + c22.x^2 [cosine case] or c21 + c23.x^2 [sine case]:
	i0 = 22 + cos_or_sin;
	qsqr= qfmul(q, q);
	sum = qfmul(qsqr,coeff[i0]);
	sum = qfadd(coeff[i0-2], sum);
	for(i = i0-4; i >= 0; i -= 2) {
		sum = qfmul(qsqr,sum);
		sum = qfadd(coeff[i], sum);
	}
	sum = qfmul(one_or_q[cos_or_sin],sum);	// Need a final mul-by-x in the sin(x) case
	return sum;
}

/***********************************************************************************/
/*                           Hyperbolic-trig functions                             */
/***********************************************************************************/

struct qfloat qfcosh(struct qfloat q)
{
	// If |x| > 2 the series needs more terms than we want to deal with, so revert to using exp(x)-construction,
	// since the exponential function allows scaling of the argument to always be O(1). That will then call this
	// function in turn to compute exp(+y) with a scaled argument y. We can be sure the recursion terminates here
	// because the scaled argument is guaranteed < 2, causing the direct summation code here to be executed:
	if((q.hi & ~MASK_SIGN) > 0x4000000000000000ull) {	// Fast |q| >= 2 check
		struct qfloat e = qfexp(q);
		return qfmul_pow2( qfadd( e,qfinv(e) ), -1 );
	} else {
		// Use positive-terms analog of cosine Taylor series summation:
		uint32 i,j, e_sum, e_new;
		struct qfloat qsqr, sum, curr_term, q1, mult;
		const int qsz = sizeof(struct qfloat);
		static struct qfloat*denoms = 0x0;
		static int first_entry = TRUE;
		if(first_entry) {
			first_entry = FALSE;
			denoms = (struct qfloat *)malloc(20*qsz);	ASSERT(denoms != NULL, "alloc failed!");
			for(i = 4; i < 38; i += 4)	// Limit largest index into QNINV[] to 40, hence (i+2) < 40 (--> i_max = 36 here)
			{
				j = (i>>1)-1;
				denoms[j-1] = qfmul(QNINV[i-1], QNINV[i  ]);
				denoms[j  ] = qfmul(QNINV[i+1], QNINV[i+2]);
			}
			// Compute one added denoms[] entry for neglected-term smallness check:
			j = (i>>1)-1;	// i = 40 on exit from above loop, so j = 19 here
			denoms[j-1] = qfmul(QNINV[i-1], QNINV[i  ]);
		}

		/* Initialize partial sum to 1 + x^2/2... */
		qsqr= qfmul(q, q);
		curr_term = qfmul_pow2(qsqr, -1);
		sum = qfadd(QONE, curr_term);

		/* get next term of series by multiplying current one by x^2/(i*(i-1)) */
		for(i = 4; i < 38; i += 4)	// Must limit largest index into QNINV[] to 40, hence (i+2) < 40
		{
			j = (i>>1)-1;
			mult      = qfmul(qsqr, denoms[j-1]);
			curr_term = qfmul(curr_term, mult);
			q1        = curr_term;

			mult      = qfmul(qsqr, denoms[j  ]);
			curr_term = qfmul(curr_term, mult);
			q1        = qfadd(q1, curr_term);	/* First calculate the sum of the current 2 terms... */
			sum       = qfadd(sum, q1);			/* ...and only then update the partial sum. */

			/* If the current term differs from the sum by more than a factor of 2^110, stop. */
			e_sum = (uint32)((sum      .hi & ~MASK_SIGN) >> 52);
			e_new = (uint32)((curr_term.hi & ~MASK_SIGN) >> 52);
			if((int32)(e_sum - e_new) > 115)
			{
			#if QFDEBUG
				printf("End cosh(x) summation after %d terms.\n", i/2);
			#endif
				break;
			}
		}
		// If used maximum #terms, do a more refined test here which checks if the *next* term is negligible:
		if(i >= 38) {	// i = 40 for cosh() series
			j = (i>>1)-1;
			mult      = qfmul(qsqr, denoms[j-1]);
			curr_term = qfmul(curr_term, mult);
			e_sum = (uint32)((sum      .hi & ~MASK_SIGN) >> 52);
			e_new = (uint32)((curr_term.hi & ~MASK_SIGN) >> 52);
			ASSERT(((int32)(e_sum - e_new) > 115), "Unconverged cosh(x) summation!");
		}
		return sum;
	}	// |x| >= 2 ?
}

struct qfloat qfsinh(struct qfloat q)
{
	// If |x| > 2 the series needs more terms than we want to deal with, revert to exp(x)-construction:
	if((q.hi & ~MASK_SIGN) > 0x4000000000000000ull) {	// Fast |q| >= 2 check
		struct qfloat e = qfexp(q);
		return qfmul_pow2( qfsub( e,qfinv(e) ), -1 );
	} else {
		// Use positive-terms analog of sine Taylor series summation:
		uint32 i,j, e_sum, e_new;
		struct qfloat qsqr, sum, curr_term, q1, mult;
		const int qsz = sizeof(struct qfloat);
		static struct qfloat*denoms = 0x0;
		static int first_entry = TRUE;
		if(first_entry) {
			first_entry = FALSE;
			denoms = (struct qfloat *)malloc(20*qsz);	ASSERT(denoms != NULL, "alloc failed!");
			for(i = 3; i < 38; i += 4)	// Must limit largest index into QNINV[] to 40, hence (i+2) < 40
			{
				j = (i>>1);
				denoms[j-1] = qfmul(QNINV[i-1], QNINV[i  ]);
				denoms[j  ] = qfmul(QNINV[i+1], QNINV[i+2]);
			}
			// Compute one added denoms[] entry for neglected-term smallness check:
			j = (i>>1);	// i = 39 on exit from above loop, so j = 19 here
			denoms[j-1] = qfmul(QNINV[i-1], QNINV[i  ]);
		}

		/* Initialize partial sum to x... */
		qsqr= qfmul(q, q);
		curr_term = sum = q;

		/* get next term of series by multiplying current one by x^2/(i*(i-1)) */
		for(i = 3; i < 38; i += 4)	// Must limit largest index into QNINV[] to 40, hence (i+2) < 40
		{
			j = (i>>1);
			mult      = qfmul(qsqr, denoms[j-1]);
			curr_term = qfmul(curr_term, mult);
			q1        = curr_term;

			mult      = qfmul(qsqr, denoms[j  ]);
			curr_term = qfmul(curr_term, mult);
			q1        = qfadd(q1, curr_term);	/* First calculate the sum of the current 2 terms... */
			sum       = qfadd(sum, q1);			/* ...and only then update the partial sum. */

			/* If the current term differs from the sum by more than a factor of 2^110, stop. */
			e_sum = (uint32)((sum      .hi & ~MASK_SIGN) >> 52);
			e_new = (uint32)((curr_term.hi & ~MASK_SIGN) >> 52);
			if((int32)(e_sum - e_new) > 115)
			{
			#if QFDEBUG
				printf("End sinh(x) summation after %d terms.\n", i/2);
			#endif
				break;
			}
		}
		// If used maximum #terms, do a more refined test here which checks if the *next* term is negligible:
		if(i >= 38) {	// i = 39 for cosh() series
			j = (i>>1);
			mult      = qfmul(qsqr, denoms[j-1]);
			curr_term = qfmul(curr_term, mult);
			e_sum = (uint32)((sum      .hi & ~MASK_SIGN) >> 52);
			e_new = (uint32)((curr_term.hi & ~MASK_SIGN) >> 52);
			ASSERT(((int32)(e_sum - e_new) > 115), "Unconverged sinh(x) summation!");
		}
		return sum;
	}	// |x| >= 2 ?
}


struct qfloat qftanh(struct qfloat x)
{
#if 0
	return qfdiv( qfsinh(x), qfcosh(x) );
#else
	struct qfloat sinh = qfsinh(x), sech = qisqrt(qfadd(QONE, qfmul(sinh,sinh)));
	return qfmul( sinh, sech );
#endif
}

// Return tanh(x) as function result, overwrite x with sinh(x)
struct qfloat qftanh_and_sinh(struct qfloat *x)
{
	struct qfloat sinh = qfsinh(*x), sech = qisqrt(qfadd(QONE, qfmul(sinh,sinh)));
	*x = sinh;
	return qfmul( sinh, sech );
}

// Return tanh(x) as function result, overwrite x with cosh(x)
struct qfloat qftanh_and_cosh(struct qfloat *x)
{
	struct qfloat cosh = qfcosh(*x), sinh = qfsqrt(qfsub(qfmul(cosh,cosh), QONE));
	*x = cosh;
	return qfdiv( sinh, cosh );
}

/*=========================== test routine: ============================*/

int qtest(void)
{
	double c, d;
	int64 hidiff;
	uint128 i128;
	struct qfloat p, q,qref,qerr, r;
	double derr, pi = 3.1415926535897932384;
#if TIMING_TEST
	uint32 i;
	int64 lodiff;
	/*...time-related stuff	*/
	clock_t clock1, clock2;
	const uint32 titers = 100000;
	const double CPU_FREQUENCY = 2000000000.0;	// One would think CLOCKS_PER_SEC would encode exactly what its name implies, but noooo......
	double td = 0, tdiff, cycles, cycles_for_qfdbl;
#endif

#ifdef MUL_LOHI64_SUBROUTINE
	printf("INFO: qfloat routines using subroutine form of MUL_LOHI\n");
#endif

#ifdef X87_ASM
	// Test long-double interface used for high-precision inits of NR seeds for certain transcendental functions:
	long double ld;
	uint64 *ld_ptr = &ld, x87_mant, x87_sexp;

	// Test I/O functions:
	ASSERT(STREQ( qf2str(QPI), "+3.14159265358979323846264338327950289 E+000" ), "I/O test failed!");

	asm ("fldln2;"
		 "fstpt %0" : "=m"(ld) : );
	x87_mant = *ld_ptr; x87_sexp = *(ld_ptr+1) & 0x000000000000FFFFull;	// Mask off high 48 bits of x87_sexp field, as these are uninited
	if(x87_mant != 0xB17217F7D1CF79ACull) {
		printf("ln2 = %30.20Le\n", ld);
		printf("x87_mant = %16" PRIx64 ", expected 0xB17217F7D1CF79ACull\n", x87_mant);	// x87_mant = b17217f7d1cf79ac, left-shift one place to off-shift hidden bit
		WARN(HERE, "Ln2 long-double mantissa conversion error", "", 0);
	}
//	ASSERT(x87_mant == 0xB17217F7D1CF79ACull, "Ln2 long-double mantissa conversion error");

//	printf("x87_sexp = %16" PRIx64 "\n", x87_sexp);	// x87_sexp = 3ffe, clear high 4 bits to get qfloat/double-compatible exp-field
	ASSERT(x87_sexp == 0x0000000000003FFEull, "Ln2 long-double exponent conversion error");

	asm ("fld1;"
		 "fadd %%st(0), %%st(0);"
		 "fsqrt;"
		 "fchs;"
		 "fstpt %0" : "=m"(ld) : );
	x87_mant = *ld_ptr; x87_sexp = *(ld_ptr+1) & 0x000000000000FFFFull;

	if(x87_mant != 0xB504F333F9DE6484ull) {
		printf("-Sqrt2 = %30.20Le\n", ld);
		printf("x87_mant = %16" PRIx64 ", expected 0xB504F333F9DE6484ull\n", x87_mant);
		WARN(HERE, "-Sqrt2 long-double mantissa conversion error", "", 0);
	}
//	ASSERT(x87_mant == 0xB504F333F9DE6484ull, "-Sqrt2 long-double mantissa conversion error");

//	printf("x87_sexp = %16" PRIx64 "\n", x87_sexp);
	ASSERT(x87_sexp == 0x000000000000BFFFull, "-Sqrt2 long-double exponent conversion error");

#endif

	ASSERT((ABS((int64)0x1234567890ABCDEFull) == 0x1234567890ABCDEFull), "ERROR 10 in qfloat.c");

	/*********** TEST THE TYPE CONVERSIONS **************/
#if TIMING_TEST
	printf	("Performing timing tests of basic qfloat operations:\n");
	clock1 = clock();
	for(i = 0; i < titers; ++i) {
		td += qfdbl(QISRT2);
		td += qfdbl(QPI);
		td += qfdbl(QLN2);
		td += qfdbl(QEXP);
	}
	clock2 = clock();
	ASSERT(td != 0.0, "!");
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= 4.0*(double)titers;
	printf	("qfdbl   : cycles/operation = %10.2f\n",cycles);
	cycles_for_qfdbl = cycles;
#endif
	c = 0.0;	d = qfdbl(QZRO);
#if QFDEBUG
		printf("dble(0.0) = %16" PRIX64 "  %16" PRIX64 "\n",f64_to_i64(c), f64_to_i64(d));
#endif
	hidiff = f64_to_i64(c) - f64_to_i64(d);	if(!(hidiff == (int64)0)) ASSERT(0,"ERROR 12 in qfloat.c");

	c = 1.0;	d = qfdbl(QONE);
#if QFDEBUG
		printf("dble(1.0) = %16" PRIX64 "  %16" PRIX64 "\n",f64_to_i64(c), f64_to_i64(d));
#endif
	hidiff = f64_to_i64(c) - f64_to_i64(d);	if(!(hidiff == (int64)0)) ASSERT(0,"ERROR 14 in qfloat.c");

	c = 2.0;	d = qfdbl(QTWO);
#if QFDEBUG
		printf("dble(2.0) = %16" PRIX64 "  %16" PRIX64 "\n",f64_to_i64(c), f64_to_i64(d));
#endif
	hidiff = f64_to_i64(c) - f64_to_i64(d);	if(!(hidiff == (int64)0)) ASSERT(0,"ERROR 16 in qfloat.c");

	c =-2.0;	d = qfdbl(qfneg(QTWO));
#if QFDEBUG
		printf("dble(-2.0)= %16" PRIX64 "  %16" PRIX64 "\n",f64_to_i64(c), f64_to_i64(d));
#endif
	hidiff = f64_to_i64(c) - f64_to_i64(d);	if(!(hidiff == (int64)0)) ASSERT(0,"ERROR 18 in qfloat.c");

	c = 2*pi;	d = qfdbl(Q2PI);
#if QFDEBUG
		printf("dble(2pi) = %16" PRIX64 "  %16" PRIX64 "\n",f64_to_i64(c), f64_to_i64(d));
#endif
	hidiff = f64_to_i64(c) - f64_to_i64(d);	if(!(ABS(hidiff) < (int64)2)) ASSERT(0,"ERROR 20 in qfloat.c");

	c =log(2.0);d = qfdbl(QLN2);
#if QFDEBUG
		printf("dble(ln2) = %16" PRIX64 "  %16" PRIX64 "\n",f64_to_i64(c), f64_to_i64(d));
#endif
	hidiff = f64_to_i64(c) - f64_to_i64(d);	if(!(ABS(hidiff) < (int64)2)) ASSERT(0,"ERROR 22 in qfloat.c");

	c = exp(1.0);
	d = qfdbl(QEXP);
#if QFDEBUG
		printf("dble(exp) = %16" PRIX64 "  %16" PRIX64 "\n",f64_to_i64(c), f64_to_i64(d));
#endif
	hidiff = f64_to_i64(c) - f64_to_i64(d);	if(!(ABS(hidiff) < (int64)2)) ASSERT(0,"ERROR 24 in qfloat.c");

	c = -c;		d = qfdbl(qfneg(QEXP));
#if QFDEBUG
		printf("dble(-exp)= %16" PRIX64 "  %16" PRIX64 "\n",f64_to_i64(c), f64_to_i64(d));
#endif
	hidiff = f64_to_i64(c) - f64_to_i64(d);	if(!(ABS(hidiff) < (int64)2)) ASSERT(0,"ERROR 26 in qfloat.c");

	/*********** TEST THE MULTIPLY ALGORITHM ************/
#if TIMING_TEST
	clock1 = clock();
	// For some reason the usual qfdbl(...) here ends up running a tad faster than the above pure qfdbl-loop,
	hidiff = 0;						// so instead just accum the hi result part in an int64.
	for(i = 0; i < titers; ++i) {
		hidiff += qfmul_pow2(QLN2,+1).hi;
		hidiff += qfmul_pow2(QLN2,+1).hi;
		hidiff += qfmul_pow2(QLN2,+1).hi;
		hidiff += qfmul_pow2(QLN2,+1).hi;
	}
	ASSERT(hidiff, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= 4.0*(double)titers;
	printf	("qfmul_pow2: cycles/operation = %8.2f\n",cycles);
#endif
#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfmul(QEXP,QEXP));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfmul   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* e*e: 	0x401D8E64B8D4DDAD, 0xCC33A3BA206B68AC	*/
	q = qfmul(QEXP,QEXP);
#if QFDEBUG
		printf("      e*e  = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
		printf("dble( e*e) = %25.16e\n",qfdbl(q));
#endif
	qref.hi = 0x401D8E64B8D4DDADull;	qref.lo = 0xCC33A3BA206B68ACull;
	// This is better than the separate hi/lo-word test, since it allows for the ROE to be of either sign:
	qerr = qfabs(qfsub(q,qref));			// Div-by-eps same as mul-by-by-2^118
	derr = qfdbl( qfmul_pow2(qerr,+118) );	// The threshold here typically needs to be ~16*[magnitude of output]
	ASSERT(derr < 64.0 ,"ERROR in QFMUL error-level check!");

	/* ln2*e:	0x3FFE258ECC242F82, 0x5DEC567E6A0E1111	*/
	q = qfmul(QLN2,QEXP);
#if QFDEBUG
		printf("     L2*e  = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
		printf("dble(L2*e) = %25.16e\n",qfdbl(q));
#endif
	qref.hi = 0x3FFE258ECC242F82ull;	qref.lo = 0x5DEC567E6A0E1111ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr < 64.0 ,"ERROR in QFMUL error-level check!");

	/* ln2^2:	0x3FDEBFBDFF82C58E, 0xA86F16B06EC97360	*/
	q = qfmul(QLN2,QLN2);
#if QFDEBUG
		printf("     L2^2  = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
		printf("dble(L2^2) = %25.16e\n",qfdbl(q));
#endif
	qref.hi = 0x3FDEBFBDFF82C58Eull;	qref.lo = 0xA86F16B06EC97360ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr < 64.0 ,"ERROR in QFMUL error-level check!");

	/* ln2*2pi:	0x40116BB24190A0B6, 0xE765BE0D06135E60	*/
	q = qfmul(QLN2,Q2PI);
#if QFDEBUG
		printf("     Ln2*pi = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
		printf("dble(Ln2*pi)= %25.16e\n",qfdbl(q));
#endif
	qref.hi = 0x40116BB24190A0B6ull;	qref.lo = 0xE765BE0D06135E60ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr < 64.0 ,"ERROR in QFMUL error-level check!");

	/* 2pi*e:	0x403114580B45D474, 0x9E6108579A2D0CA7	*/
	q = qfmul(Q2PI,QEXP);
#if QFDEBUG
		printf("     pi*e  = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
		printf("dble(pi*e) = %25.16e\n",qfdbl(q));
#endif
	qref.hi = 0x403114580B45D474ull;	qref.lo = 0x9E6108579A2D0CA7ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr < 128.0 ,"ERROR in QFMUL error-level check!");

	/* 2pi*2pi:	0x4043BD3CC9BE45DE, 0x5A4ADC4D9B301183	*/
	q = qfmul(Q2PI,Q2PI);
#if QFDEBUG
		printf("  (2*pi)^2 = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
		printf("dble(2pi^2)= %25.16e\n",qfdbl(q));
		printf("dble(2pi^2)= %25.16e\n",pi*pi);
#endif
	qref.hi = 0x4043BD3CC9BE45DEull;	qref.lo = 0x5A4ADC4D9B301183ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr < 128.0 ,"ERROR in QFMUL error-level check!");

	/*********** TEST THE ADDITION ALGORITHM ************/
#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfadd(QEXP,QEXP));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfadd   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* 2*pi+e:	0x402200C04CE72C66, 0x7821CB48D9B947AC	*/
	q = qfadd(QEXP,Q2PI);
#if QFDEBUG
		printf("  2*pi + e = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x402200C04CE72C66ull;	qref.lo = 0x7821CB48D9B947ACull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr < 64.0 ,"ERROR in QFMUL error-level check!");

	/********** TEST THE SUBTRACTION ALGORITHM **********/
#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfsub(QEXP,QEXP));
	}
	ASSERT(td == 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfsub   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* Both inputs   normalized, output   normalized, with just one significant bit. */
	q.hi = 0x3FEFFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFFull;
	q = qfsub(q, q);
#if QFDEBUG
		printf("result1 = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = qref.lo = 0x0000000000000000ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr == 0.0 ,"ERROR in QFSUB error-level check!");

	p.hi = 0x3FEFFFFFFFFFFFFFull;	p.lo = 0xFFFFFFFFFFFFFFFFull;
	q.hi = 0x3FEFFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFEull;
	q = qfsub(p, q);
#if QFDEBUG
		printf("result2 = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x38A0000000000000ull;	qref.lo = 0x0000000000000000ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr == 0.0 ,"ERROR in QFSUB error-level check!");

	/* Both inputs   normalized, output denormalized, with just one significant bit. */
	p.hi = 0x00FFFFFFFFFFFFFFull;	p.lo = 0xFFFFFFFFFFFFFFFFull;
	q.hi = 0x00FFFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFEull;
	q = qfsub(p, q);
#if QFDEBUG
		printf("result3 = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x0000000000000000ull;	qref.lo = 0x0000000000004000ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr == 0.0 ,"ERROR in QFSUB error-level check!");

	/* Both inputs denormalized, output zero */
	q.hi = 0x000FFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFFull;
	q = qfsub(q, q);
#if QFDEBUG
		printf("result4 = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = qref.lo = 0ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr == 0.0 ,"ERROR in QFSUB error-level check!");

	/* Both inputs denormalized, output denormalized, with just one significant bit. */
	p.hi = 0x000FFFFFFFFFFFFFull;	p.lo = 0xFFFFFFFFFFFFFFFFull;
	q.hi = 0x000FFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFEull;
	q = qfsub(p, q);
#if QFDEBUG
		printf("result5 = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0ull;	qref.lo = 1ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr == 0.0 ,"ERROR in QFSUB error-level check!");

	/* 2*pi-e:	0x400C84EC1D7402C7, 0x39DB360DDEDB4F60	*/
	q = qfsub(Q2PI,QEXP);
#if QFDEBUG
		printf("    2pi- e = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x400C84EC1D7402C7ull;	qref.lo = 0x39DB360DDEDB4F60ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr == 0.0 ,"ERROR in QFSUB error-level check!");

	/* e-2*pi:	0xC00C84EC1D7402C7, 0x39DB360DDEDB4F60	*/
	r = qfsub(QEXP,Q2PI);
#if QFDEBUG
		printf("     e-2pi = %16" PRIX64 "  %16" PRIX64 "\n",r.hi,r.lo);
#endif
	if(!(qfcmpeq(r, qfneg(q)))) ASSERT(0,"ERROR 54 in qfloat.c");

	/*********** TEST THE SQUARE ROOT ALGORITHM ************/
#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfsqrt(QEXP));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfsqrt  : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* sqrt(2):	0x3FF6A09E667F3BCC, 0x908B2FB1366EA958, qfsqrt gives ...956. */
	q = qfsqrt(QTWO);
#if QFDEBUG
		printf("sqrt(2) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x3FF6A09E667F3BCCull;	qref.lo = 0x908B2FB1366EA958ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFSQRT error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFSQRT error-level check!");

	/*********** TEST THE INVERSION AND DIVISION ALGORITHMS ************/
#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfinv(QEXP));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfinv   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* 1/(2*pi):0x3FC45F306DC9C882, 0xA53F84EAFA3EA69B(B81B...), qfinv gives ...698. */
	q = qfinv(Q2PI);
#if QFDEBUG
		printf("1/(2*pi) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x3FC45F306DC9C882ull;	qref.lo = 0xA53F84EAFA3EA69Bull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFINV error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFINV error-level check!");

	/* 1/e:		0x3FD78B56362CEF37, 0xC6AEB7B1E0A4153E(4376...), qfinv gives ...53C. */
	q = qfinv(QEXP);
#if QFDEBUG
		printf("1/e      = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x3FD78B56362CEF37ull;	qref.lo = 0xC6AEB7B1E0A4153Eull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFINV error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFINV error-level check!");

	/* 1/ln2:	0x3FF71547652B82FE, 0x1777D0FFDA0D23A7(D11D...), qfinv gives ...3A6. */
	q = qfinv(QLN2);
#if QFDEBUG
		printf("1/ln(2)  = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x3FF71547652B82FEull;	qref.lo = 0x1777D0FFDA0D23A7ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFINV error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFINV error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfdiv(QEXP,QPI));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfdiv   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* 2*pi/ln2:0x40222123045B5DEB, 0x9C5398CE82C06E4B(80DB...), qfdiv gives ...E4A. */
	q = qfdiv(Q2PI, QLN2);
#if QFDEBUG
		printf("2*pi/ln(2)  = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x40222123045B5DEBull;	qref.lo = 0x9C5398CE82C06E4Bull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 128.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFDIV error-level check!", "", 0);
	}
//	ASSERT(derr < 128.0 ,"ERROR in QFDIV error-level check!");

	/*********** TEST THE TRANSCENDENTAL FUNCTIONS ************/
#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfsn1(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfsin   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* sin(Pi/4):	Compare to precomputed 1/sqrt2: */
	q = qfsn1(QPI4TH);
	qref.hi = QISRT2.hi;	qref.lo = QISRT2.lo;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFSIN error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFSIN error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfcs1(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfcos   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* cos(Pi/4):	Compare to precomputed 1/sqrt2: */
	q = qfcs1(QPI4TH);
	qref.hi = QISRT2.hi;	qref.lo = QISRT2.lo;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFCOS error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFCOS error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qftan(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qftan   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* tan(Pi/4):	Compare to 1: */
	q = qftan(QPI4TH);
#if QFDEBUG
		printf("qtfan(PI/4) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = QONE.hi;	qref.lo = QONE.lo;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFTAN error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFTAN error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfcot(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfcot   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* cot(Pi/4):	Compare to 1: */
	q = qfcot(QPI4TH);
#if QFDEBUG
		printf("qfcot(PI/4) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = QONE.hi;	qref.lo = QONE.lo;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFCOT error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFCOT error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfatan(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfatan  : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* atan(1):	Compare to precomputed Pi/4: */
	q = qfatan(QONE);
#if QFDEBUG
		printf("qatan(1) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = QPI4TH.hi;	qref.lo = QPI4TH.lo;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFATAN error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFATAN error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qflog(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qflog   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* log(2):	Compare to precomputed QLN2 = {0x3FE62E42FEFA39EFull, 0x35793C7673007E5Full}: */
	q = qflog(QTWO);
#if QFDEBUG
		printf("qlog(2) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x3FE62E42FEFA39EFull;	qref.lo = 0x35793C7673007E5Full;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	ASSERT(derr < 1100.0,"ERROR in QFLOG error-level check!");	// AGM-based log is fast but error-prone

	/* log(2^64):	Compare to precomputed log(2^64) = (same as log(2) but exp-field += 6): */
	q = qfmul_pow2(QONE,+64);
	q = qflog(q);

#if QFDEBUG
		printf("qlog(2^64) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x40462E42FEFA39EFull;	qref.lo = 0x35793C7673007E5Full;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 1100.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFLOG error-level check!", "", 0);
	}
//	ASSERT(derr < 1100.0 ,"ERROR in QFLOG error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfexp(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfexp   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* exp(1):	0x4005BF0A8B145769, 0x5355FB8AC404E7A7(9E3B...), qfexp gives ...4E7A7, ~116 bits of accuracy. */
	q = qfexp(QONE);
#if QFDEBUG
		printf("qexp(1) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x4005BF0A8B145769ull;	qref.lo = 0x5355FB8AC404E7A7ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr > 64.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFEXP error-level check!", "", 0);
	}
//	ASSERT(derr <= 64.0 ,"ERROR in QFEXP error-level check!");

	/* Sine and cosine are somewhat roundoff-error prone, so raise the error limit slightly. */
	/* cos(1):	0x3FE14A280FB5068B, 0x923848CDB2ED0E37(A534...), qfcs1 gives ...D0E38, ~116 bits of accuracy */
	q = qfcs1(QONE);
#if QFDEBUG
		printf("qcs1(1) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x3FE14A280FB5068Bull;	qref.lo = 0x923848CDB2ED0E37ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFCS1 error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFCS1 error-level check!");

	r = qfcos(QONE);
#if QFDEBUG
		printf("qcos(1) = %16" PRIX64 "  %16" PRIX64 "\n",r.hi,r.lo);
#endif
	if(!(qfcmpeq(r, q))) ASSERT(0,"ERROR 70 in qfloat.c");

	/* sin(1):	0x3FEAED548F090CEE, 0x0418DD3D2138A1E7(8651...), qfsn1 gives ...8A1E9, ~115 bits of accuracy */
	q = qfsn1(QONE);
#if QFDEBUG
		printf("qsn1(1) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x3FEAED548F090CEEull;	qref.lo = 0x0418DD3D2138A1E7ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFSN1 error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFSN1 error-level check!");

	r = qfsin(QONE);
#if QFDEBUG
		printf("qsin(1) = %16" PRIX64 "  %16" PRIX64 "\n",r.hi,r.lo);
#endif
	if(!(qfcmpeq(r, q))) ASSERT(0,"ERROR 74 in qfloat.c");

	/* cos(100):0x3FEB981DBF665FDF, 0x63F433736617A041(5D8A...), qfcos gives ...7A023, ~114 bits of accuracy */
	q = qfcos(i64_to_q((int64)100));
#if QFDEBUG
		printf("qcos(100) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0x3FEB981DBF665FDFull;	qref.lo = 0x63F433736617A041ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 128.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFCOS error-level check!", "", 0);
	}
//	ASSERT(derr < 128.0 ,"ERROR in QFCOS error-level check!");

	/* sin(100):0xBFE03425B78C4DB8, 0x0708F6155D083EB2(1C6B...), qfsin gives ...83EE5, ~109 bits of accuracy */
	q = qfsin(i64_to_q((int64)100));
#if QFDEBUG
		printf("qsin(100) = %16" PRIX64 "  %16" PRIX64 "\n",q.hi,q.lo);
#endif
	qref.hi = 0xBFE03425B78C4DB8ull;	qref.lo = 0x0708F6155D083EB2ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 128.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFSIN error-level check!", "", 0);
	}
//	ASSERT(derr < 128.0 ,"ERROR in QFSIN error-level check!");

	/*********** Test the hyperbolic-trigs: **********************/
#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfsinh(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfsinh  : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* sinh(ln(2)):	Compare to 0.75: */
	q = qfsinh(QLN2);
	qref.hi = 0x3FE8000000000000ull;	qref.lo = 0x0000000000000000ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFSINH error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFSINH error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qfcosh(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfcosh  : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* cosh(ln(2)):	Compare to 1.25: */
	q = qfcosh(QLN2);
	qref.hi = 0x3FF4000000000000ull;	qref.lo = 0x0000000000000000ull;
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFCOSH error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFCOSH error-level check!");

#if TIMING_TEST
	clock1 = clock();
	td = 0.0;
	for(i = 0; i < titers; ++i) {
		td += qfdbl(qftanh(QLN2));
	}
	ASSERT(td != 0.0, "!");
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qftanh  : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	/* tanh(ln(2)):	Compare to 3/5: */
	q = qftanh(QLN2);
	qref = qfmul( i64_to_q(3ull), QNINV[5] );
	qerr = qfabs(qfsub(q,qref));	derr = qfdbl( qfmul_pow2(qerr,+118) );
	if(derr >= 16.0) {
		printf("derr = %10.5f\n", derr);
		WARN(HERE, "ERROR in QFTANH error-level check!", "", 0);
	}
//	ASSERT(derr < 16.0 ,"ERROR in QFTANH error-level check!");

	/*********** TEST THE INT --> QFLOAT and ROUND-TOWARD-ZERO AND ROUND-TO-NEAREST FUNCTIONS ************/
	ASSERT(CMPEQ128( qfint(qfneg( i64_to_q(  0ull))), NIL128 ), "error!");
	ASSERT(CMPEQ128( qfint(qfneg(i128_to_q(NIL128))), NIL128 ), "error!");
#if TIMING_TEST
	clock1 = clock();
	hidiff = lodiff = 0ull;
	for(i = 0; i < titers; ++i) {
		i128 = qfnint(QLN2);
		hidiff += i128.d1;
		lodiff += i128.d0;
	}
	ASSERT(!hidiff && (lodiff == titers), "!");	// NINT(ln2) = 1, titers times
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfnint  : cycles/operation = %10.2f\n",cycles);
#endif
	q = qfmul_pow2(QONE, -1);
	i128 = qfnint(q);
#if QFDEBUG
		printf("qfnint(0.5) = %16" PRIX64 "  %16" PRIX64 "\n",i128.d1,i128.d0);
#endif
	ASSERT((!i128.d1 && i128.d0 == (uint64)1),"ERROR 80 in qfloat.c");

#if TIMING_TEST
	clock1 = clock();
	hidiff = lodiff = 0ull;
	for(i = 0; i < titers; ++i) {
		i128 = qfint(QLN2);
		hidiff += i128.d1;
		lodiff += i128.d0 + qfint(QPI).d0;
	}
	ASSERT(!hidiff && (lodiff == 3*titers), "!");	// INT(ln2) = 0 and INT(pi) = 3, summed (titers) times
	clock2 = clock();
	tdiff = (double)(clock2 - clock1);
	cycles = tdiff*CPU_FREQUENCY/CLOCKS_PER_SEC;
	cycles /= (double)titers;
	printf	("qfint   : cycles/operation = %10.2f\n",cycles - cycles_for_qfdbl);
#endif
	i128 = qfnint(QHALF);	ASSERT((!i128.d1 && i128.d0 == (uint64)1),"ERROR 82 in qfloat.c");
	i128 = qfint(QHALF);	ASSERT((!i128.d1 && i128.d0 == (uint64)0),"ERROR 83 in qfloat.c");
	i128 = qfnint(QEXP);	ASSERT((!i128.d1 && i128.d0 == (uint64)3),"ERROR 84 in qfloat.c");
	i128 = qfint(QEXP);		ASSERT((!i128.d1 && i128.d0 == (uint64)2),"ERROR 85 in qfloat.c");
	i128 = qfnint(Q2PI);	ASSERT((!i128.d1 && i128.d0 == (uint64)6),"ERROR 86 in qfloat.c");
	i128 = qfint(Q2PI);		ASSERT((!i128.d1 && i128.d0 == (uint64)6),"ERROR 87 in qfloat.c");
	q = qfmul_pow2(Q2PI, 20);
	i128 = qfnint(q);		ASSERT((!i128.d1 && i128.d0 == (uint64)6588397),"ERROR 90 in qfloat.c");

	q = qfmul_pow2(QPI, 125);	/* This gives pi*2^125, which should still fit into a signed 128-bit int. */
	i128 = qfnint(q);
	ASSERT((i128.d1 = (uint64)0x6487ED5110B4611Aull && i128.d0 == (uint64)0x62633145C06E1000ull),"ERROR 92 in qfloat.c");

#if TIMING_TEST
	exit(0);
#endif
	return 0;
}

