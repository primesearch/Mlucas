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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef qfloat_h_included
#define qfloat_h_included

#include "Mdata.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
Implements a 128-bit floating-point emulation.

***NOTE: Any program that includes this header must call qtest() prior to calling
any other qfloat functions!***
*/

/* Integer constants needed to manipulate upper half of a qfloat: */

/* XOR high part of a qfloat with this to flip sign. */
#define MASK_SIGN	(uint64)0x8000000000000000ull
/* Use to mask off sign bit of a right-justied sign/exponent combo. */
#define MASK_EXP	(uint64)0x00000000000007FFull
/* Use to mask off significand. */
#define MASK_MANT	(uint64)0x000FFFFFFFFFFFFFull
/* Hidden bit, for normalized numbers. */
#define TWOE52		(uint64)0x0010000000000000ull

struct qfloat{
	/*
	2004: In the current version, for reasons of debug and simpler code development,
	the Upper 64 bits contain a bit pattern identical to an IEEE754 64-bit float.
	That is:
	hi<63:63> - Sign bit S. 1 for negative number, 0 for nonnegative.
	hi<52:62> - 11-bit biased binary exponent E, treated as unsigned int in [0,2047].
	hi< 0:51> - High 52 bits (excluding hidden bit) of (53+64)-bit unsigned significand
			Mhi, binary point assumed at left of uppermost (hidden) bit, which is assumed
			to be a one, unless the exponent field E is zero, which indicates underflow
			has occurred (see below.)
	lo< 0:63> - Low 64 bits of significand Mlo, also treated as unsigned.

	That is, if the biased exponent E > 0, the (normalized) number represented is

		(-1)^S * 2^(E-1023) * 2^(-53)*[2^52 + Mhi + 2^(-64)*Mlo] .

	If the biased exponent E = 0, the (denormalized) number represented is

		(-1)^S * 2^(E-1023) * 2^(-53)*[       Mhi + 2^(-64)*Mlo] .

	TODO: Allow exponent fields with > 11 bits. (IEEE754 mandates >= 15 bits.)
	In fact, there is a practical advantage to having < 53 high-word mantissa bits,
	namely that, using a floating-mul (FMUL) emulation of the MULH64 instruction
	(with a regular MULL64 alongside it to get the lowest 64 bits of the product), we can
	exactly calculate product of up to a 52-bit x 64-bit unsigned integer - the FMUL
	generates a 53-bit result and the MULL64 a 64-bit result, but we need one overlap bit
	between the two to effect error correction of the least-significant bit of the
	floating-point product. This floating-point MULH64 emulation allows the quad multiply
	to be fast even on hardware lacking an intrinsic MULH64-like instruction and/or where
	integer multiply is slow.

	TODO: Implement NaN and SNaN support, together with exception-return value codes.
	*/
	uint64 hi;
	uint64 lo;
};

/* Useful constants in qfloat form: */
extern const struct qfloat QZRO;
extern const struct qfloat QONE;
extern const struct qfloat QEPS;
extern const struct qfloat QHALF;
extern const struct qfloat QTWO;
extern const struct qfloat QTHREE;
extern const struct qfloat Q2PI;
extern const struct qfloat QPI;
extern const struct qfloat QPIHALF;
extern const struct qfloat QIPIHLF;
extern const struct qfloat QPI4TH;
extern const struct qfloat QLN2;
extern const struct qfloat QEXP;
extern const struct qfloat QSQRT2;
extern const struct qfloat QISRT2;
extern const struct qfloat QI64AGM;
extern const struct qfloat QSNAN;	// signaling NaN

/*********************/
/* Qfloat functions: */
/*********************/
// Self-tests:
int    qtest		(void);
// I/O:
char* qf2str(struct qfloat q);
// Negation, absolute value:
struct qfloat qfneg	(const struct qfloat q);
struct qfloat qfabs	(const struct qfloat q);
// Comparisons:
uint32 qfiszero		(const struct qfloat q);
uint32 qfcmpeq		(const struct qfloat q1, struct qfloat q2);
uint32 qfcmpne		(const struct qfloat q1, struct qfloat q2);
uint32 qfcmplt		(const struct qfloat q1, struct qfloat q2);
uint32 qfcmpgt		(const struct qfloat q1, struct qfloat q2);
uint32 qfcmple		(const struct qfloat q1, struct qfloat q2);
uint32 qfcmpge		(const struct qfloat q1, struct qfloat q2);
// qfloat --> double / long-double conversions:
uint64 qfdbl_as_uint64(const struct qfloat q);
double qfdbl		(const struct qfloat q);
double qfdbl_wrong_way_rnd(struct qfloat q);
double dbl_flip_lsb	(double d);
long double qfldbl	(const struct qfloat q);	// This should only be used for x87 code.
// double / long-double --> qfloat conversions:
struct qfloat dbl_to_q	(double d);
struct qfloat ldbl_to_q	(long double ld);		// This should only be used for x87 code.
// int --> qfloat conversion:
struct qfloat i64_to_q	(int64 i64);
struct qfloat i128_to_q	(uint128 i);
// Multiply by power of 2 ("binary shift", but with added exponent/significand handling needed)
struct qfloat qfmul_pow2(struct qfloat q, int32 pow);
// Round-to-nearest / Round toward zero - Note these return result as a twos-comp 128-bit *integer*:
uint128 qfnint	(struct qfloat q);
uint128 qfint	(struct qfloat q);
// q1*q2:
struct qfloat qfmul	(struct qfloat q1, struct qfloat q2);
// q1/q2:
struct qfloat qfdiv	(struct qfloat q1, struct qfloat q2);
// qfloat approximation to rational p/q:
struct qfloat qf_rational_quotient(int64 p, int64 q);

// Increment & decrement
struct qfloat qfinc(struct qfloat q);
struct qfloat qfdec(struct qfloat q);

/* Top-level add, subtract and absolute-difference routines seen by the caller - these examine the
   signs of the inputs, send the proper combination of +-q1 and +-q2 to the
   low-level qfsum and qfdif functions, and properly sign the result. */
struct qfloat qfadd	(struct qfloat q1, struct qfloat q2);
struct qfloat qfsub	(struct qfloat q1, struct qfloat q2);
struct qfloat qfabsdiff(struct qfloat q1, struct qfloat q2);
/* Low-level add and subtract routines, which assume both inputs are nonnegative.
NOTE: these are intended to be called ONLY from qfadd and qfsub. */
struct qfloat qfsum	(struct qfloat q1, struct qfloat q2);
struct qfloat qfdif	(struct qfloat q1, struct qfloat q2);
/* Multiplicative inverse: */
struct qfloat qfinv	(struct qfloat q);
/* Square root and inverse thereof: */
struct qfloat qfsqrt(struct qfloat q);
struct qfloat qisqrt(struct qfloat q);
/* AGM-iteration utility: */
struct qfloat qfagm	(struct qfloat x, struct qfloat y);
/* Natural exponential and log, base-10 log: */
struct qfloat qfexp	(struct qfloat q);
struct qfloat qflog	(struct qfloat q);
struct qfloat qflog10(struct qfloat q);
struct qfloat qfatan(struct qfloat x);
/* Factorial: */
struct qfloat qffact(uint32 n);

/* Top-level trigonometric routines seen by the caller - these examine the
   sign and magnitude of the input and map that to the proper call to either
   +-qfsn1(arg) or +-qfcs1(arg), where arg is in [0, pi/2). */
struct qfloat qfcos	(struct qfloat q);
struct qfloat qfsin	(struct qfloat q);
struct qfloat qftan (struct qfloat q);
struct qfloat qfcot (struct qfloat q);
struct qfloat qftan_and_sin(struct qfloat *x);
struct qfloat qftan_and_cos(struct qfloat *x);
struct qfloat qfcosh(struct qfloat q);
struct qfloat qfsinh(struct qfloat q);
struct qfloat qftanh(struct qfloat q);
struct qfloat qftanh_and_sinh(struct qfloat *x);
struct qfloat qftanh_and_cosh(struct qfloat *x);

/* low-level sine and cosine routines - these require an argument in [0, pi/2). */
struct qfloat qfcs1	(struct qfloat q);
struct qfloat qfsn1	(struct qfloat q);
// The above 2 both feed either q or pi/2-q to this one, which requires arg in [0, pi/4):
struct qfloat qfcos_or_sin1(struct qfloat q, int cos_or_sin);

// Utility macros:
#define QLEADZ(__x)	( leadz64(__x.hi) + ((-(sint32)(__x.hi == 0)) && leadz64(__x.lo)) )

/* Left-shift: This should move at most a bit into the lowest exp-field, so check that on output: */
#define QLSHIFT(__x, __n, __y)\
{\
	/* Make sure sign/exp fields have been cleared and shift count >= 0: */\
	ASSERT((__x.hi>>52) == 0,"QLSHIFT: sign/exp fields not zero!");\
	ASSERT((int64)__n >= 0,"QLSHIFT: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.hi = ((uint64)__x.hi);\
		__y.lo = ((uint64)__x.lo);\
	}\
	else if(__n < 64)\
	{\
		__y.hi = ((uint64)__x.hi << __n) + ((uint64)__x.lo >> (64-__n));\
		__y.lo = ((uint64)__x.lo << __n);\
	}\
	else if(__n < 128)\
	{\
		__y.hi = ((uint64)__x.lo << (__n-64));\
		__y.lo = (uint64)0;\
	}\
	else\
	{\
		__y.hi = (uint64)0;\
		__y.lo = (uint64)0;\
	}\
	/* Make sure exp field at most 1 after shift: */\
	ASSERT((__x.hi>>52) <= 1,"QLSHIFT: exp field out of range on output!");\
}

/* (Logical) Right-shift: */
#define QRSHIFT(__x, __n, __y)\
{\
	/* Make sure sign/exp fields have been cleared and shift count >= 0: */\
	ASSERT((__x.hi>>52) == 0,"QRSHIFT:  sign/exp fields not zero!");\
	ASSERT((int64)(__n) >= 0,"QRSHIFT: (int64)(__n) >= 0 !");\
	/* Need to handle zero shift count separately: */\
	if((__n) == 0)\
	{\
		__y.lo = ((uint64)__x.lo);\
		__y.hi = ((uint64)__x.hi);\
	}\
	else if((__n) < 64)\
	{\
		__y.lo = ((uint64)__x.lo >> (__n)) + ((uint64)__x.hi << (64-(__n)));\
		__y.hi = ((uint64)__x.hi >> (__n));\
	}\
	else if((__n) < 128)\
	{\
		__y.lo = ((uint64)__x.hi >> ((__n)-64));\
		__y.hi = (uint64)0;\
	}\
	else\
	{\
		__y.lo = (uint64)0;\
		__y.hi = (uint64)0;\
	}\
}

#ifdef __cplusplus
}
#endif

#endif	/* qfloat_h_included */

