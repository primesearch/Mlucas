/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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

#include "Mdata.h"

/*
Implements a 128-bit floating-point emulation.

***NOTE: Any program that includes this header must call qtest() prior to calling
any other qfloat functions!***
*/

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef qfloat_h_included
#define qfloat_h_included

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
extern const struct qfloat QTWO;
extern const struct qfloat QTHREE;
extern const struct qfloat Q2PI;
extern const struct qfloat QPI;
extern const struct qfloat QPIHALF;
extern const struct qfloat QLN2;
extern const struct qfloat QEXP;
extern const struct qfloat QSQRT2;

/* Qfloat functions: */
int    qtest		(void);
struct qfloat qfneg	(struct qfloat q);
struct qfloat qfabs	(struct qfloat q);
uint32 qfcmpeq		(struct qfloat q1, struct qfloat q2);
uint32 qfcmplt		(struct qfloat q1, struct qfloat q2);
uint32 qfcmple		(struct qfloat q1, struct qfloat q2);
double qfdbl		(struct qfloat q);
struct qfloat dbl_to_q	(double d);
struct qfloat i64_to_q	(int64 i64);
struct qfloat qfmul_pow2(struct qfloat q, int32 pow);
struct qfloat qfnint	(struct qfloat q);
struct qfloat qfint	(struct qfloat q);
struct qfloat qfmul	(struct qfloat q1, struct qfloat q2);
struct qfloat qfdiv	(struct qfloat q1, struct qfloat q2);
struct qfloat qf_rational_quotient(int64 p, int64 q);

/* Top-level add and subtract routines seen by the caller - these examine the
   signs of the inputs, send the proper combination of +-q1 and +-q2 to the
   low-level qfsum and qfdif functions, and properly sign the result. */
struct qfloat qfadd	(struct qfloat q1, struct qfloat q2);
struct qfloat qfsub	(struct qfloat q1, struct qfloat q2);
/* Low-level add and subtract routines, which assume both inputs are nonnegative.
NOTE: these are intended to be called ONLY from qfadd and qfsub. */
struct qfloat qfsum	(struct qfloat q1, struct qfloat q2);
struct qfloat qfdif	(struct qfloat q1, struct qfloat q2);
/* Multiplicative inverse: */
struct qfloat qfinv	(struct qfloat q);
/* Square root: */
struct qfloat qfsqrt(struct qfloat q);
/* Exponential: */
struct qfloat qfexp	(struct qfloat q);
/* Top-level sine and cosine routines seen by the caller - these examine the
   sign and magnitude of the input and map that to the proper call to either
   +-qfsn1(arg) or +-qfcs1(arg), where arg is in [0, pi/2). */
struct qfloat qfcos	(struct qfloat q);
struct qfloat qfsin	(struct qfloat q);
/* low-level sine and cosine routines - these require an argument in [0, pi/2). */
struct qfloat qfcs1	(struct qfloat q);
struct qfloat qfsn1	(struct qfloat q);

#endif	/* qfloat_h_included */
