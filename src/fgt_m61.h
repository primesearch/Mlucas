/*******************************************************************************
*                                                                              *
*   (C) 1997-2019 by Ernst W. Mayer.                                           *
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
#ifndef fgt_m61_h_included
#define fgt_m61_h_included

#include "util.h"

// Our modulus q = 2^61 - 1:

/***************/
// NB: Since args to these reduce macros will more often than not be expressions (e.g. qreduce(x - y + q4)),
// start each by copying arg into a local uint64, to ensure that any input expression only gets evaluated once:
/*
Returns x (mod q), but in the sense of a possible partial modular reduction: Outputs are in [0, B], where B = q+7.
Note: if x = q, QREDUCE returns q, not zero.
*/
#define qreduce(x)	\
	({ uint64 tmp = x;	\
		tmp = (tmp >> 61) + (tmp & 0x1FFFFFFFFFFFFFFFull);	\
		tmp;	})

// ...or this if you want to finish reducing a qreduce() output:
#define qreduce_finish(x)	\
	({ uint64 tmp = x;	\
		tmp -= (-(uint64)(tmp >= 0x1FFFFFFFFFFFFFFFull)) & 0x1FFFFFFFFFFFFFFFull;	\
		tmp;	})

// Use this if you require a guaranteed-full reduction of x (mod q)...
#define qreduce_full(x)	\
	({ uint64 tmp = x;	\
		tmp = (tmp >> 61) + (tmp & 0x1FFFFFFFFFFFFFFFull);	\
		tmp -= (-(uint64)(tmp >= 0x1FFFFFFFFFFFFFFFull)) & 0x1FFFFFFFFFFFFFFFull;	\
		tmp;	})

/***************/

/*
Returns sqrt(1/2)*x (mod q).
sqrt(1/2) == 2^30 mod q, so the multiply can be effected via 2 shifts, an AND, and an add.
For normalized inputs (< q), Output is in [0, B30], where B30 = q + 7*2^30 = 2^61 + 2^33 - 2^30 - 1.
*/
#define mul_i2(x)	(((x) << 30) & 0x1FFFFFFFFFFFFFFFull) + ((x) >> 31)

/***************/

/*
Returns sqrt(2)*x (mod q).
sqrt(2) == 2^31 mod q, so the multiply can be effected via 2 shifts, an AND, and an add.
Outputs are in [0, B31], where B31 = q + 7*2^31 = 2^61 + 2^34 - 2^31 - 1.
*/
#define mul_s2(x)	(((x) << 31) & 0x1FFFFFFFFFFFFFFFull) + ((x) >> 30)

/***************/

/*
Returns 2^n * x (mod q). x is a uint64; The shift count n is assumed to be any kind of int, with value in [0,61].

If x only partially normalized (i.e. in [0, b]) on entry and n = 0, result is fully normalized, i.e. xout in [0,q].
If x unnormalized on entry and n = 0, the result is partially normalized, i.e. xout in [0,b].
The special case n = 61 leaves x unchanged.

For general operands x in [0,2^64-1] and n in [0,60], ((x << n) & q) is in [0, q - (2^n - 1)] = [0, 2^61 - 2^n]
and (x >> (61-n)) is in [0, 2^(3+n) - 1]. The sum is bounded above by 2^61 - 2^n + 2^(3+n) - 1 = q + 2^(3+n) - 2^n.

OK, let`s do some crude estimation for non-normalized inputs:

The sum is maximized for x = 2^64-1 and n = 60, giving 2^63 - 1 + 2^60 = 9*2^60 - 1 ~= 4.5*q,
i.e. inputs approximately in [0,8q] yield outputs approximately in [0,5q].
For x = 2^64-1 and n = 59, the sum is bounded by ~2.75*q, etc., approaching q+7 from above.

x = 2^63-1 and n = 60 gives q + 2^62 - 2^60 ~= 2.5*q .
x = 2^63-1 and n = 59 gives q + 2^61 - 2^59 ~= 1.75*q . This case is important in the between-forward-and-inverse-FFT
							pair_square step, where we multiply inputs in [0,4q] by the modular inverse of 4 == 2^59.
x = 2^62-1 and n = 60 gives 2^61 - 1 + 2^60 = 3*2^60 - 1 ~= 1.5*q, i.e. inputs approximately in [0,2q]
																	yield outputs approximately in [0,2q].

NEGATIVE POWERS OF 2:

The modular analog of 1/2 (call it w) satisfies 2*w == 1 (mod q), thus w = (q+1)/2 = 2^60. More generally,
any negative-integer power of 2 (mod q) satisfies 2^(-p) == 2^(61-p), with p < 61. We obtain the same
result by simply analogizing the mul_pow2_modq macro to negative powers, and thus can effect multiply
by 2^(-p) by simply calling the mul_pow2_modq macro with power-of-2 argument (61-p).

Thus e.g. to effect a modular x*(1/2) we call mul_pow2_modq(x,60).
*/
#define mul_pow2_modq(x,n)	(((x) << n) & 0x1FFFFFFFFFFFFFFFull) + ((x) >> (61-n))

/****** Prototypes for functions defined in fgt_m61.c are collected in util.h *******/

#endif	/* fgt_m61_h_included */
