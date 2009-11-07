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

#include "factor.h"

//#define DBG_ASSERT ASSERT

#ifdef USE_SSE2
	#include "align.h"

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

#endif

//#define FAC_DEBUG	1
#if FAC_DEBUG
	char char_buf[1024], str0[64], str1[64];
#endif

/* Useful Macros whose scope is this file only: */

/* Converts a 78-bit unsigned input __x (stored in a uint96)
to balanced-digit floating-point form. Outputs have the following size ranges:

	fword0,1 in [-2^25, +2^25]
	fword2   in [   -1, +2^26]
*/
/* Subroutine version is handy for debugging: */
void cvt_uint78_3word_double(uint96 __x, double*__fword0, double*__fword1, double*__fword2)
{
	uint64 __tmp64;
	 int64 __itmp, __cy;

	DBG_ASSERT(HERE, (__x.d1 >> 14) == 0, "cvt_uint78_3word_double : Input > 78-bit limit!");

	/* Digit 0: */
	__tmp64 = __x.d0;
	__itmp = __tmp64 & 0x0000000003ffffff;
	/* Is the current digit >= (base/2)? */
	__cy = (int64)((uint64)__itmp>>25);	/* Cast to unsigned to ensure logical right-shift */
	/* If yes, balance it by subtracting the base: */
	/* RHS terms must be signed to prevent integer underflow-on-subtract: */
	*__fword0 = (double)(__itmp - (__cy<<26));

	/* Digit 1: */
	__tmp64 = (__tmp64 >> 26) + __cy;\
	__itmp = __tmp64 & 0x0000000003ffffff;
	/* Is the current digit >= (base/2)? */
	__cy = (int64)((uint64)__itmp>>25);	/* Cast to unsigned to ensure logical right-shift */
	/* If yes, balance it by subtracting the base: */
	/* RHS terms must be signed to prevent integer underflow-on-subtract: */
	*__fword1 = (double)(__itmp - (__cy<<26));

	/* Digit 2: */
	__tmp64 = (__tmp64 >> 26) + __cy;\
	__tmp64 += (__x.d1 << 12);	/* 12 = (64 - 2*26) */
	/* No balanced-digit normalization of MSW: */
	*__fword2 = (double)__tmp64;

	DBG_ASSERT(HERE, *__fword2 <= TWO26FLOAT, "cvt_uint78_3word_double : MSW > TWO26FLOAT");
}

#define CVT_UINT78_3WORD_DOUBLE(__x, __fword0, __fword1, __fword2)\
{\
	uint64 __tmp64;\
	 int64 __itmp, __cy;\
	\
	DBG_ASSERT(HERE, (__x.d1 >> 14) == 0, "CVT_UINT78_3WORD_DOUBLE : Input > 78-bit limit!");\
	\
	/* Digit 0: */\
	__tmp64 = __x.d0;\
	__itmp = __tmp64 & 0x0000000003ffffff;\
	/* Is the current digit >= (base/2)? */\
	__cy = (int64)((uint64)__itmp>>25);	/* Cast to unsigned to ensure logical right-shift */\
	/* If yes, balance it by subtracting the base: */\
	/* RHS terms must be signed to prevent integer underflow-on-subtract: */\
	__fword0 = (double)(__itmp - (__cy<<26));\
	\
	/* Digit 1: */\
	__tmp64 = (__tmp64 >> 26) + __cy;\
	__itmp = __tmp64 & 0x0000000003ffffff;\
	/* Is the current digit >= (base/2)? */\
	__cy = (int64)((uint64)__itmp>>25);	/* Cast to unsigned to ensure logical right-shift */\
	/* If yes, balance it by subtracting the base: */\
	/* RHS terms must be signed to prevent integer underflow-on-subtract: */\
	__fword1 = (double)(__itmp - (__cy<<26));\
	\
	/* Digit 2: */\
	__tmp64 = (__tmp64 >> 26) + __cy;\
	__tmp64 += (__x.d1 << 12);	/* 12 = (64 - 2*26) */\
	/* No balanced-digit normalization of MSW: */\
	__fword2 = (double)__tmp64;\
	\
	DBG_ASSERT(HERE, __fword2 <= TWO26FLOAT, "CVT_UINT78_3WORD_DOUBLE : MSW > TWO26FLOAT");\
}

/* Converts a 78-bit unsigned input __x (stored in balanced-digit
floating-point form) to a uint96. Assumes the FP input is properly normalized.
*/
/* Subroutine version is handy for debugging: */
void cvt78_3word_double_uint96(double __fword0, double __fword1, double __fword2, uint96*__x)
{
	int64 __itmp, __cy;

	/* Cast current digit to int64 form, subtracting any borrow from previous digit: */
	__itmp = (int64)__fword0;
	if(__itmp < 0)	/* If current digit < 0, add the base and set carry = -1	*/
	{
		__itmp += TWO26FLOAT;
		DBG_ASSERT(HERE, __itmp >= 0, "cvt78_3word_double_uint96 : Normalized digit still < 0!");
		__cy = -1;
	}
	else
	{
		__cy = 0;
	}
	__x->d0 = (uint64)__itmp;

	/* Digit 1: */
	__itmp = (int64)__fword1 +  __cy;
	if(__itmp < 0)
	{
		__itmp += TWO26FLOAT;
		DBG_ASSERT(HERE, __itmp >= 0, "cvt78_3word_double_uint96 : Normalized digit still < 0!");
		__cy = -1;
	}
	else
	{
		__cy = 0;
	}
	__x->d0 += ((uint64)__itmp << 26);

	/* Digit 2: */
	__itmp = (int64)__fword2 +  __cy;
	if(__itmp < 0)
	{
		__itmp += TWO26FLOAT;
		DBG_ASSERT(HERE, __itmp >= 0, "cvt78_3word_double_uint96 : Normalized digit still < 0!");
		__cy = -1;
	}
	else
	{
		__cy = 0;
	}
	__x->d0 += ((uint64)__itmp << 52);
	__x->d1  = ((uint64)__itmp >> 12);	/* Only case where we really need the (uint64) cast */

	DBG_ASSERT(HERE, (__x->d1 >> 14) == 0, "cvt78_3word_double_uint96 : Output > 78-bit limit!");
	DBG_ASSERT(HERE,  __cy           == 0, "cvt78_3word_double_uint96 : Nonzero exit carry!");
}

#define CVT78_3WORD_DOUBLE_UINT96(__fword0, __fword1, __fword2, __x)\
{\
	int64 __itmp, __cy;\
	\
	/* Cast current digit to int64 form, subtracting any borrow from previous digit: */\
	__itmp = (int64)__fword0;\
	if(__itmp < 0)	/* If current digit < 0, add the base and set carry = -1	*/\
	{\
		__itmp += TWO26FLOAT;\
		DBG_ASSERT(HERE, __itmp >= 0, "CVT78_3WORD_DOUBLE_UINT96 : Normalized digit still < 0!");\
		__cy = -1;\
	}\
	else\
	{\
		__cy = 0;\
	}\
	__x.d0 = (uint64)__itmp;\
\
	/* Digit 1: */\
	__itmp = (int64)__fword1 +  __cy;\
	if(__itmp < 0)\
	{\
		__itmp += TWO26FLOAT;\
		DBG_ASSERT(HERE, __itmp >= 0, "CVT78_3WORD_DOUBLE_UINT96 : Normalized digit still < 0!");\
		__cy = -1;\
	}\
	else\
	{\
		__cy = 0;\
	}\
	__x.d0 += ((uint64)__itmp << 26);\
\
	/* Digit 2: */\
	__itmp = (int64)__fword2 +  __cy;\
	if(__itmp < 0)\
	{\
		__itmp += TWO26FLOAT;\
		DBG_ASSERT(HERE, __itmp >= 0, "CVT78_3WORD_DOUBLE_UINT96 : Normalized digit still < 0!");\
		__cy = -1;\
	}\
	else\
	{\
		__cy = 0;\
	}\
	__x.d0 += ((uint64)__itmp << 52);\
	__x.d1  = ((uint64)__itmp >> 12) & 0x0000000000003fff;	/* Only case where we really need the (uint64) cast */\
	\
	DBG_ASSERT(HERE, (__x.d1 >> 14) == 0, "CVT78_3WORD_DOUBLE_UINT96 : Output > 78-bit limit!");\
	DBG_ASSERT(HERE,  __cy          == 0, "CVT78_3WORD_DOUBLE_UINT96 : Nonzero exit carry!");\
}

/* Takes a 78-bit unsigned input __x stored in balanced-digit floating-point form
and renormalizes with respect to the balanced-digit base.
*/
#define NORMALIZE78_3WORD_DOUBLE(__x0, __x1, __x2)\
{\
	double __fcy;\
	\
	/* Digit 0: */\
	__fcy = NINT(__x0*TWO26FLINV);\
	__x0 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__x1 += __fcy;\
	__fcy = NINT(__x1*TWO26FLINV);\
	__x1 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__x2 += __fcy;\
	\
	DBG_ASSERT(HERE, __x2 <= TWO26FLOAT, "NORMALIZE_UINT78_3WORD_DOUBLE : MSW > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x2 >= 0        , "NORMALIZE_UINT78_3WORD_DOUBLE : MSW < 0!");\
}

/* Takes a 156-bit unsigned input __x stored in balanced-digit floating-point form
and renormalizes with respect to an 78-bit = (26,26,26)-bit balanced-digit base.
Because we expect that we may wind up using the upper and lower halves of the result
separately, we require the MSW of each to be nonnegative, i.e. we don't balance __x2.
*/
#define NORMALIZE156_6WORD_DOUBLE(__x0, __x1, __x2, __x3, __x4, __x5)\
{\
	double __fcy;\
	\
	/* Digit 0: */\
	__fcy = NINT(__x0*TWO26FLINV);\
	__x0 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__x1 += __fcy;\
	__fcy = NINT(__x1*TWO26FLINV);\
	__x1 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__x2 += __fcy;\
	__fcy = NINT(__x2*TWO26FLINV);\
	__x2 -= __fcy*TWO26FLOAT;\
	if(__x2 < 0)\
	{\
		__x2 += TWO26FLOAT;\
		__fcy--;\
	}\
	\
	/* Digit 3: */\
	__x3 += __fcy;\
	__fcy = NINT(__x3*TWO26FLINV);\
	__x3 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__x4 += __fcy;\
	__fcy = NINT(__x4*TWO26FLINV);\
	__x4 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__x5 += __fcy;\
	\
	DBG_ASSERT(HERE, __x2 >= 0         , "NORMALIZE_UINT156_6WORD_DOUBLE : _x2 < 0!");\
	DBG_ASSERT(HERE, __x2 <= TWO26FLOAT, "NORMALIZE_UINT156_6WORD_DOUBLE : _x2 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x5 >= 0         , "NORMALIZE_UINT156_6WORD_DOUBLE : MSW < 0!");\
	DBG_ASSERT(HERE, __x5 <= TWO26FLOAT, "NORMALIZE_UINT156_6WORD_DOUBLE : MSW > TWO26FLOAT");\
}

/**********************************************************************************/
/* Balanced-digit FP arithmetic with base 2^26 allows MULTIPLY inputs up to 2^78. */
/**********************************************************************************/

/*...Square of a 78-bit input __x .
Lower and upper halves of __x^2 are returned in __lo and __hi, respectively.
Current version needs 14 FMUL, 22 FADD.

Because we expect that we may wind up using the upper and lower halves of the result
separately, we require the MSW of each to be nonnegative, i.e. we don't balance __x2.

ASSUMES:
	- None of the input and output addresses coincide;
*/
/* Subroutine version is handy for debugging: */
void sqr_lohi78_3word_double(double __x0,double __x1,double __x2, double*__prod0,double*__prod1,double*__prod2,double*__prod3,double*__prod4,double*__prod5)
{
	double __2x0 = __x0 + __x0, __2x1 = __x1 + __x1;\
	double __fcy;
	uint32 __itmp;

	DBG_ASSERT(HERE, __x0 < TWO26FLOAT, "sqr_lohi78_3word_double : x0 > TWO26FLOAT");
	DBG_ASSERT(HERE, __x1 < TWO26FLOAT, "sqr_lohi78_3word_double : x1 > TWO26FLOAT");
	DBG_ASSERT(HERE, __x2 < TWO26FLOAT, "sqr_lohi78_3word_double : x2 > TWO26FLOAT");

	/* Digit 0: */
	*__prod0  =  __x0*__x0;
	__fcy    = NINT(*__prod0 * TWO26FLINV);
	*__prod0 -= __fcy*TWO26FLOAT;

	/* Digit 1: */
	*__prod1  = __2x0*__x1 + __fcy;
	__fcy    = NINT(*__prod1 * TWO26FLINV);
	*__prod1 -= __fcy*TWO26FLOAT;

	/* Digit 2: */
	*__prod2  = __2x0*__x2 + __x1*__x1 + __fcy;
	__fcy    = NINT(*__prod2 * TWO26FLINV);
	*__prod2 -= __fcy*TWO26FLOAT;
	/* Branchless sequence to unbalance the *__prod2 term: */
	__itmp   = (*__prod2 < 0);
	__fcy   -= (double)__itmp;
	*__prod2 += (double)(__itmp << 26);

	/* Digit 3: */
	*__prod3  = __2x1*__x2 + __fcy;
	__fcy    = NINT(*__prod3 * TWO26FLINV);
	*__prod3 -= __fcy*TWO26FLOAT;

	/* Digit 4: */
	*__prod4  =  __x2*__x2 + __fcy;
	__fcy    = NINT(*__prod4 * TWO26FLINV);
	*__prod4 -= __fcy*TWO26FLOAT;

	/* Digit 5: */
	*__prod5  = __fcy;

	DBG_ASSERT(HERE, *__prod2 >= 0         , "sqr_lohi78_3word_double : _x2 < 0!");
	DBG_ASSERT(HERE, *__prod2 <= TWO26FLOAT, "sqr_lohi78_3word_double : _x2 > TWO26FLOAT");
	DBG_ASSERT(HERE, *__prod5 >= 0         , "sqr_lohi78_3word_double : MSW < 0!");
	DBG_ASSERT(HERE, *__prod5 <= TWO26FLOAT, "sqr_lohi78_3word_double : MSW > TWO26FLOAT");
}

#define SQR_LOHI78_3WORD_DOUBLE(__x0,__x1,__x2, __prod0,__prod1,__prod2,__prod3,__prod4,__prod5)\
{\
	double __2x0 = __x0 + __x0, __2x1 = __x1 + __x1;\
	double __fcy;\
	uint32 __itmp;\
	\
	DBG_ASSERT(HERE, __x0 < TWO26FLOAT, "SQR_LOHI78_3WORD_DOUBLE : x0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x1 < TWO26FLOAT, "SQR_LOHI78_3WORD_DOUBLE : x1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x2 < TWO26FLOAT, "SQR_LOHI78_3WORD_DOUBLE : x2 > TWO26FLOAT");\
	\
	/* Digit 0: */\
	__prod0  =  __x0*__x0;\
	__fcy    = NINT(__prod0*TWO26FLINV);\
	__prod0 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__prod1  = __2x0*__x1 + __fcy;\
	__fcy    = NINT(__prod1*TWO26FLINV);\
	__prod1 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__prod2  = __2x0*__x2 + __x1*__x1 + __fcy;\
	__fcy    = NINT(__prod2*TWO26FLINV);\
	__prod2 -= __fcy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __prod2 term: */\
	__itmp   = (__prod2 < 0);\
	__fcy   -= (double)__itmp;\
	__prod2 += (double)(__itmp << 26);\
	\
	/* Digit 3: */\
	__prod3  = __2x1*__x2 + __fcy;\
	__fcy    = NINT(__prod3*TWO26FLINV);\
	__prod3 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__prod4  =  __x2*__x2 + __fcy;\
	__fcy    = NINT(__prod4*TWO26FLINV);\
	__prod4 -= __fcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__prod5  = __fcy;\
	\
	DBG_ASSERT(HERE, __prod2 >= 0         , "SQR_LOHI78_3WORD_DOUBLE : _x2 < 0!");\
	DBG_ASSERT(HERE, __prod2 <= TWO26FLOAT, "SQR_LOHI78_3WORD_DOUBLE : _x2 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __prod5 >= 0         , "SQR_LOHI78_3WORD_DOUBLE : MSW < 0!");\
	DBG_ASSERT(HERE, __prod5 <= TWO26FLOAT, "SQR_LOHI78_3WORD_DOUBLE : MSW > TWO26FLOAT");\
}

#define SQR_LOHI78_3WORD_DOUBLE_q2(\
  __fx0,__fx1,__fx2, __fprod0,__fprod1,__fprod2,__fprod3,__fprod4,__fprod5\
, __gx0,__gx1,__gx2, __gprod0,__gprod1,__gprod2,__gprod3,__gprod4,__gprod5\
)\
{\
	double __f2x0 = __fx0 + __fx0, __f2x1 = __fx1 + __fx1;\
	double __g2x0 = __gx0 + __gx0, __g2x1 = __gx1 + __gx1;\
	double __fcy, __gcy;\
	uint32 __itmp, __jtmp;\
	\
	/* Digit 0: */\
	__fprod0  =  __fx0*__fx0;\
	__fcy     = NINT(__fprod0*TWO26FLINV);\
	__fprod0 -= __fcy*TWO26FLOAT;\
	\
	__gprod0  =  __gx0*__gx0;\
	__gcy     = NINT(__gprod0*TWO26FLINV);\
	__gprod0 -= __gcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__fprod1  = __f2x0*__fx1 + __fcy;\
	__fcy     = NINT(__fprod1*TWO26FLINV);\
	__fprod1 -= __fcy*TWO26FLOAT;\
	\
	__gprod1  = __g2x0*__gx1 + __gcy;\
	__gcy     = NINT(__gprod1*TWO26FLINV);\
	__gprod1 -= __gcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__fprod2  = __f2x0*__fx2 + __fx1*__fx1 + __fcy;\
	__fcy     = NINT(__fprod2*TWO26FLINV);\
	__fprod2 -= __fcy*TWO26FLOAT;\
	__itmp    = (__fprod2 < 0);\
	__fcy    -= (double)__itmp;\
	__fprod2 += (double)(__itmp << 26);\
	\
	__gprod2  = __g2x0*__gx2 + __gx1*__gx1 + __gcy;\
	__gcy     = NINT(__gprod2*TWO26FLINV);\
	__gprod2 -= __gcy*TWO26FLOAT;\
	__jtmp    = (__gprod2 < 0);\
	__gcy    -= (double)__jtmp;\
	__gprod2 += (double)(__jtmp << 26);\
	\
	/* Digit 3: */\
	__fprod3  = __f2x1*__fx2 + __fcy;\
	__fcy     = NINT(__fprod3*TWO26FLINV);\
	__fprod3 -= __fcy*TWO26FLOAT;\
	\
	__gprod3  = __g2x1*__gx2 + __gcy;\
	__gcy     = NINT(__gprod3*TWO26FLINV);\
	__gprod3 -= __gcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__fprod4  =  __fx2*__fx2 + __fcy;\
	__fcy     = NINT(__fprod4*TWO26FLINV);\
	__fprod4 -= __fcy*TWO26FLOAT;\
	\
	__gprod4  =  __gx2*__gx2 + __gcy;\
	__gcy     = NINT(__gprod4*TWO26FLINV);\
	__gprod4 -= __gcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__fprod5  = __fcy;\
	__gprod5  = __gcy;\
}

#define SQR_LOHI78_3WORD_DOUBLE_q4(\
  __fx0,__fx1,__fx2, __fprod0,__fprod1,__fprod2,__fprod3,__fprod4,__fprod5\
, __gx0,__gx1,__gx2, __gprod0,__gprod1,__gprod2,__gprod3,__gprod4,__gprod5\
, __hx0,__hx1,__hx2, __hprod0,__hprod1,__hprod2,__hprod3,__hprod4,__hprod5\
, __ix0,__ix1,__ix2, __iprod0,__iprod1,__iprod2,__iprod3,__iprod4,__iprod5\
)\
{\
	double __f2x0 = __fx0 + __fx0, __f2x1 = __fx1 + __fx1;\
	double __g2x0 = __gx0 + __gx0, __g2x1 = __gx1 + __gx1;\
	double __h2x0 = __hx0 + __hx0, __h2x1 = __hx1 + __hx1;\
	double __i2x0 = __ix0 + __ix0, __i2x1 = __ix1 + __ix1;\
	double __fcy, __gcy, __hcy, __icy;\
	uint32 __itmp, __jtmp, __ktmp, __ltmp;\
	\
	/* Digit 0: */\
	__fprod0  =  __fx0*__fx0;							__gprod0  =  __gx0*__gx0;							__hprod0  =  __hx0*__hx0;							__iprod0  =  __ix0*__ix0;\
	__fcy     = (__fprod0*TWO26FLINV);	__gcy     = (__gprod0*TWO26FLINV);	__hcy     = (__hprod0*TWO26FLINV);	__icy     = (__iprod0*TWO26FLINV);\
	__fprod0 -= __fcy*TWO26FLOAT;						__gprod0 -= __gcy*TWO26FLOAT;						__hprod0 -= __hcy*TWO26FLOAT;						__iprod0 -= __icy*TWO26FLOAT;\
	/* Digit 1: */\
	__fprod1  = __f2x0*__fx1 + __fcy;					__gprod1  = __g2x0*__gx1 + __gcy;					__hprod1  = __h2x0*__hx1 + __hcy;					__iprod1  = __i2x0*__ix1 + __icy;\
	__fcy     = (__fprod1*TWO26FLINV);	__gcy     = (__gprod1*TWO26FLINV);	__hcy     = (__hprod1*TWO26FLINV);	__icy     = (__iprod1*TWO26FLINV);\
	__fprod1 -= __fcy*TWO26FLOAT;						__gprod1 -= __gcy*TWO26FLOAT;						__hprod1 -= __hcy*TWO26FLOAT;						__iprod1 -= __icy*TWO26FLOAT;\
	/* Digit 2: */\
	__fprod2  = __f2x0*__fx2 + __fx1*__fx1 + __fcy;		__gprod2  = __g2x0*__gx2 + __gx1*__gx1 + __gcy;		__hprod2  = __h2x0*__hx2 + __hx1*__hx1 + __hcy;		__iprod2  = __i2x0*__ix2 + __ix1*__ix1 + __icy;\
	__fcy     = (__fprod2*TWO26FLINV);	__gcy     = (__gprod2*TWO26FLINV);	__hcy     = (__hprod2*TWO26FLINV);	__icy     = (__iprod2*TWO26FLINV);\
	__fprod2 -= __fcy*TWO26FLOAT;						__gprod2 -= __gcy*TWO26FLOAT;						__hprod2 -= __hcy*TWO26FLOAT;						__iprod2 -= __icy*TWO26FLOAT;\
	__itmp    = (__fprod2 < 0);							__jtmp    = (__gprod2 < 0);							__ktmp    = (__hprod2 < 0);							__ltmp    = (__iprod2 < 0);\
	__fcy    -= (double)__itmp;							__gcy    -= (double)__jtmp;							__hcy    -= (double)__ktmp;							__icy    -= (double)__ltmp;\
	__fprod2 += (double)(__itmp << 26);					__gprod2 += (double)(__jtmp << 26);					__hprod2 += (double)(__ktmp << 26);					__iprod2 += (double)(__ltmp << 26);\
	/* Digit 3: */\
	__fprod3  = __f2x1*__fx2 + __fcy;					__gprod3  = __g2x1*__gx2 + __gcy;					__hprod3  = __h2x1*__hx2 + __hcy;					__iprod3  = __i2x1*__ix2 + __icy;\
	__fcy     = (__fprod3*TWO26FLINV);	__gcy     = (__gprod3*TWO26FLINV);	__hcy     = (__hprod3*TWO26FLINV);	__icy     = (__iprod3*TWO26FLINV);\
	__fprod3 -= __fcy*TWO26FLOAT;						__gprod3 -= __gcy*TWO26FLOAT;						__hprod3 -= __hcy*TWO26FLOAT;						__iprod3 -= __icy*TWO26FLOAT;\
	/* Digit 4: */\
	__fprod4  =  __fx2*__fx2 + __fcy;					__gprod4  =  __gx2*__gx2 + __gcy;					__hprod4  =  __hx2*__hx2 + __hcy;					__iprod4  =  __ix2*__ix2 + __icy;\
	__fcy     = (__fprod4*TWO26FLINV);	__gcy     = (__gprod4*TWO26FLINV);	__hcy     = (__hprod4*TWO26FLINV);	__icy     = (__iprod4*TWO26FLINV);\
	__fprod4 -= __fcy*TWO26FLOAT;						__gprod4 -= __gcy*TWO26FLOAT;						__hprod4 -= __hcy*TWO26FLOAT;						__iprod4 -= __icy*TWO26FLOAT;\
	/* Digit 5: */\
	__fprod5  = __fcy;									__gprod5  = __gcy;									__hprod5  = __hcy;									__iprod5  = __icy;\
}

/* Lower half of __x * __y .
Current version needs 12 FMUL, 14 FADD.

Because we may desire to overwrite one of the two sets of inputs with the outputs,
we code so that any or all of __X, __Y and __LO may have the same addresses.
*/
void mull78_3word_double(double __x0,double __x1,double __x2, double __y0,double __y1,double __y2, double*__lo0,double*__lo1,double*__lo2)
{
	double __fcy, __prod0, __prod1, __prod2;

	DBG_ASSERT(HERE, __x0 < TWO26FLOAT, "mull78_3word_double : x0 > TWO26FLOAT");
	DBG_ASSERT(HERE, __x1 < TWO26FLOAT, "mull78_3word_double : x1 > TWO26FLOAT");
	DBG_ASSERT(HERE, __x2 < TWO26FLOAT, "mull78_3word_double : x2 > TWO26FLOAT");

	DBG_ASSERT(HERE, __y0 < TWO26FLOAT, "mull78_3word_double : y0 > TWO26FLOAT");
	DBG_ASSERT(HERE, __y1 < TWO26FLOAT, "mull78_3word_double : y1 > TWO26FLOAT");
	DBG_ASSERT(HERE, __y2 < TWO26FLOAT, "mull78_3word_double : y2 > TWO26FLOAT");

	/* Precompute all the needed partial products: */
	__prod0  = __x0*__y0;
	__prod1  = __x0*__y1 + __x1*__y0;
	__prod2  = __x0*__y2 + __x1*__y1 + __x2*__y0;

	/* Digit 0: */
	__fcy    = NINT(__prod0*TWO26FLINV);
	*__lo0   = __prod0 - __fcy*TWO26FLOAT;

	/* Digit 1: */
	__prod1 += __fcy;
	__fcy    = NINT(__prod1*TWO26FLINV);
	*__lo1   = __prod1 - __fcy*TWO26FLOAT;

	/* Digit 2: */
	__prod2 += __fcy;
	__fcy    = NINT(__prod2*TWO26FLINV);
	*__lo2   = __prod2 - __fcy*TWO26FLOAT;
	/* Require output to be nonnegative, so leave MSW unbalanced: */
	if(*__lo2 < 0)
	{
		*__lo2 += TWO26FLOAT;
	}

	DBG_ASSERT(HERE, *__lo2 >= 0, "mull78_3word_double : MSW < 0!");
}

#define MULL78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2, __lo0,__lo1,__lo2)\
{\
	double __fcy, __prod0, __prod1, __prod2;\
	uint32 __itmp;\
	\
	DBG_ASSERT(HERE, __x0 < TWO26FLOAT, "MULL78_3WORD_DOUBLE : x0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x1 < TWO26FLOAT, "MULL78_3WORD_DOUBLE : x1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x2 < TWO26FLOAT, "MULL78_3WORD_DOUBLE : x2 > TWO26FLOAT");\
	\
	DBG_ASSERT(HERE, __y0 < TWO26FLOAT, "MULL78_3WORD_DOUBLE : y0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __y1 < TWO26FLOAT, "MULL78_3WORD_DOUBLE : y1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __y2 < TWO26FLOAT, "MULL78_3WORD_DOUBLE : y2 > TWO26FLOAT");\
	\
	/* Precompute all the needed partial products: */\
	__prod0  = __x0*__y0;\
	__prod1  = __x0*__y1 + __x1*__y0;\
	__prod2  = __x0*__y2 + __x1*__y1 + __x2*__y0;\
	\
	/* Digit 0: */\
	__fcy    = NINT(__prod0*TWO26FLINV);\
	__lo0    = __prod0 - __fcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__prod1 += __fcy;\
	__fcy    = NINT(__prod1*TWO26FLINV);\
	__lo1    = __prod1 - __fcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__prod2 += __fcy;\
	__fcy    = NINT(__prod2*TWO26FLINV);\
	__lo2    = __prod2 - __fcy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __lo2 term: */\
	__itmp   = (__lo2 < 0);\
	__lo2   += (double)(__itmp << 26);\
	/* Require output to be nonnegative, so leave MSW unbalanced: */\
	DBG_ASSERT(HERE, __lo2 >= 0, "MULL78_3WORD_DOUBLE : MSW < 0!");\
}

#define MULL78_3WORD_DOUBLE_q2(\
  __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __flo0,__flo1,__flo2\
, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __glo0,__glo1,__glo2\
)\
{\
	double __fcy, __fprod0, __fprod1, __fprod2;\
	double __gcy, __gprod0, __gprod1, __gprod2;\
	uint32 __itmp, __jtmp;\
	\
	/* Precompute all the needed partial products: */\
	__fprod0  = __fx0*__fy0;\
	__fprod1  = __fx0*__fy1 + __fx1*__fy0;\
	__fprod2  = __fx0*__fy2 + __fx1*__fy1 + __fx2*__fy0;\
	\
	__gprod0  = __gx0*__gy0;\
	__gprod1  = __gx0*__gy1 + __gx1*__gy0;\
	__gprod2  = __gx0*__gy2 + __gx1*__gy1 + __gx2*__gy0;\
	\
	/* Digit 0: */\
	__fcy     = NINT(__fprod0*TWO26FLINV);\
	__flo0    = __fprod0 - __fcy*TWO26FLOAT;\
	\
	__gcy     = NINT(__gprod0*TWO26FLINV);\
	__glo0    = __gprod0 - __gcy*TWO26FLOAT;\
	\
	/* Digit 1: */\
	__fprod1 += __fcy;\
	__fcy     = NINT(__fprod1*TWO26FLINV);\
	__flo1    = __fprod1 - __fcy*TWO26FLOAT;\
	\
	__gprod1 += __gcy;\
	__gcy     = NINT(__gprod1*TWO26FLINV);\
	__glo1    = __gprod1 - __gcy*TWO26FLOAT;\
	\
	/* Digit 2: */\
	__fprod2 += __fcy;\
	__fcy     = NINT(__fprod2*TWO26FLINV);\
	__flo2    = __fprod2 - __fcy*TWO26FLOAT;\
	\
	__gprod2 += __gcy;\
	__gcy     = NINT(__gprod2*TWO26FLINV);\
	__glo2    = __gprod2 - __gcy*TWO26FLOAT;\
	\
	/* Branchless sequence to unbalance the __fprod2 term: */\
	__itmp    = (__flo2 < 0);\
	__flo2 += (double)(__itmp << 26);\
	\
	__jtmp    = (__glo2 < 0);\
	__glo2 += (double)(__jtmp << 26);\
}

#define MULL78_3WORD_DOUBLE_q4(\
  __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __flo0,__flo1,__flo2\
, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __glo0,__glo1,__glo2\
, __hx0,__hx1,__hx2, __hy0,__hy1,__hy2, __hlo0,__hlo1,__hlo2\
, __ix0,__ix1,__ix2, __iy0,__iy1,__iy2, __ilo0,__ilo1,__ilo2\
)\
{\
	double __fprod0, __fprod1, __fprod2;\
	double __gprod0, __gprod1, __gprod2;\
	double __hprod0, __hprod1, __hprod2;\
	double __iprod0, __iprod1, __iprod2;\
	double __fcy, __gcy, __hcy, __icy;\
	uint32 __itmp, __jtmp, __ktmp, __ltmp;\
	\
	/* Precompute all the needed partial products: */\
	__fprod0  = __fx0*__fy0;							__gprod0  = __gx0*__gy0;							__hprod0  = __hx0*__hy0;							__iprod0  = __ix0*__iy0;							\
	__fprod1  = __fx0*__fy1 + __fx1*__fy0;				__gprod1  = __gx0*__gy1 + __gx1*__gy0;				__hprod1  = __hx0*__hy1 + __hx1*__hy0;				__iprod1  = __ix0*__iy1 + __ix1*__iy0;				\
	__fprod2  = __fx0*__fy2 + __fx1*__fy1 + __fx2*__fy0;__gprod2  = __gx0*__gy2 + __gx1*__gy1 + __gx2*__gy0;__hprod2  = __hx0*__hy2 + __hx1*__hy1 + __hx2*__hy0;__iprod2  = __ix0*__iy2 + __ix1*__iy1 + __ix2*__iy0;\
	/* Digit 0: */\
	__fcy     = (__fprod0*TWO26FLINV);	__gcy     = (__gprod0*TWO26FLINV);	__hcy     = (__hprod0*TWO26FLINV);	__icy     = (__iprod0*TWO26FLINV);\
	__flo0    = __fprod0 - __fcy*TWO26FLOAT;			__glo0    = __gprod0 - __gcy*TWO26FLOAT;			__hlo0    = __hprod0 - __hcy*TWO26FLOAT;			__ilo0    = __iprod0 - __icy*TWO26FLOAT;\
	/* Digit 1: */\
	__fprod1 += __fcy;									__gprod1 += __gcy;									__hprod1 += __hcy;									__iprod1 += __icy;\
	__fcy     = (__fprod1*TWO26FLINV);	__gcy     = (__gprod1*TWO26FLINV);	__hcy     = (__hprod1*TWO26FLINV);	__icy     = (__iprod1*TWO26FLINV);\
	__flo1    = __fprod1 - __fcy*TWO26FLOAT;			__glo1    = __gprod1 - __gcy*TWO26FLOAT;			__hlo1    = __hprod1 - __hcy*TWO26FLOAT;			__ilo1    = __iprod1 - __icy*TWO26FLOAT;\
	/* Digit 2: */\
	__fprod2 += __fcy;									__gprod2 += __gcy;									__hprod2 += __hcy;									__iprod2 += __icy;\
	__fcy     = (__fprod2*TWO26FLINV);	__gcy     = (__gprod2*TWO26FLINV);	__hcy     = (__hprod2*TWO26FLINV);	__icy     = (__iprod2*TWO26FLINV);\
	__flo2    = __fprod2 - __fcy*TWO26FLOAT;			__glo2    = __gprod2 - __gcy*TWO26FLOAT;			__hlo2    = __hprod2 - __hcy*TWO26FLOAT;			__ilo2    = __iprod2 - __icy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __fprod2 term: */\
	__itmp    = (__flo2 < 0);							__jtmp    = (__glo2 < 0);							__ktmp    = (__hlo2 < 0);							__ltmp    = (__ilo2 < 0);\
	__flo2 += (double)(__itmp << 26);					__glo2 += (double)(__jtmp << 26);					__hlo2 += (double)(__ktmp << 26);					__ilo2 += (double)(__ltmp << 26);\
}

/* Upper half of __x * __y .
Current version needs 16 FMUL, 20 FADD.

Because we may desire to overwrite one of the two sets of inputs with the outputs,
we code so that any or all of __X, __Y and __LO may have the same addresses.
*/
#define MULH78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2, __hi0,__hi1,__hi2)\
{\
	double __fcy, __tmp, __prod3, __prod4;\
	uint32 __itmp;\
	\
	DBG_ASSERT(HERE, __x0 < TWO26FLOAT, "MULH78_3WORD_DOUBLE : x0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x1 < TWO26FLOAT, "MULH78_3WORD_DOUBLE : x1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __x2 < TWO26FLOAT, "MULH78_3WORD_DOUBLE : x2 > TWO26FLOAT");\
	\
	DBG_ASSERT(HERE, __y0 < TWO26FLOAT, "MULH78_3WORD_DOUBLE : y0 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __y1 < TWO26FLOAT, "MULH78_3WORD_DOUBLE : y1 > TWO26FLOAT");\
	DBG_ASSERT(HERE, __y2 < TWO26FLOAT, "MULH78_3WORD_DOUBLE : y2 > TWO26FLOAT");\
	\
	/* Digit 0: */\
	__tmp  =  __x0*__y0;\
	__fcy  = NINT(__tmp*TWO26FLINV);\
	\
	/* Digit 1: */\
	__tmp  = __x0*__y1 + __x1*__y0 + __fcy;\
	__fcy  = NINT(__tmp*TWO26FLINV);\
	\
	/* Digit 2: */\
	__tmp  = __x0*__y2 + __x1*__y1 + __x2*__y0 + __fcy;\
	__fcy  = NINT(__tmp*TWO26FLINV);\
	__tmp -= __fcy*TWO26FLOAT;\
	/* Branchless sequence to unbalance the __prod2 term: */\
	__itmp = (__tmp < 0);\
	__fcy -= (double)__itmp;\
	/* Require low half to be nonnegative, so leave this term unbalanced: */\
	/*if(__tmp < 0)	*/\
	/*{				*/\
	/*	__fcy--;	*/\
	/*}				*/\
	\
	/* At this point the possibility of same-address in-and-outputs comes into play: */\
	/* Precompute all the needed partial products: */\
	__prod3 = __x1*__y2 + __x2*__y1 + __fcy;\
	__prod4 = __x2*__y2;\
	\
	/* Digit 3: */\
	__fcy    = NINT(__prod3*TWO26FLINV);\
	__hi0    = __prod3 - __fcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__prod4 += __fcy;\
	__fcy    = NINT(__prod4*TWO26FLINV);\
	__hi1    = __prod4 - __fcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__hi2    = __fcy;\
	\
	DBG_ASSERT(HERE, __hi2 >= 0, "MULH78_3WORD_DOUBLE : MSW < 0!");\
}

#define MULH78_3WORD_DOUBLE_q2(\
  __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __fhi0,__fhi1,__fhi2\
, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __ghi0,__ghi1,__ghi2\
  )\
{\
	double __fcy, __ftmp, __fprod3, __fprod4;\
	double __gcy, __gtmp, __gprod3, __gprod4;\
	uint32 __itmp,__jtmp;\
	\
	/* Digit 0: */\
	__ftmp =  __fx0*__fy0;\
	__fcy  = NINT(__ftmp*TWO26FLINV);\
	\
	__gtmp =  __gx0*__gy0;\
	__gcy  = NINT(__gtmp*TWO26FLINV);\
	\
	/* Digit 1: */\
	__ftmp = __fx0*__fy1 + __fx1*__fy0 + __fcy;\
	__fcy  = NINT(__ftmp*TWO26FLINV);\
	\
	__gtmp = __gx0*__gy1 + __gx1*__gy0 + __gcy;\
	__gcy  = NINT(__gtmp*TWO26FLINV);\
	\
	/* Digit 2: */\
	__ftmp = __fx0*__fy2 + __fx1*__fy1 + __fx2*__fy0 + __fcy;\
	__fcy  = NINT(__ftmp*TWO26FLINV);\
	__ftmp-= __fcy*TWO26FLOAT;\
	__itmp = (__ftmp < 0);\
	__fcy -= (double)__itmp;\
	\
	__gtmp = __gx0*__gy2 + __gx1*__gy1 + __gx2*__gy0 + __gcy;\
	__gcy  = NINT(__gtmp*TWO26FLINV);\
	__gtmp-= __gcy*TWO26FLOAT;\
	__jtmp = (__gtmp < 0);\
	__gcy -= (double)__jtmp;\
	\
	/* At this point the possibility of same-address in-and-outputs comes into play: */\
	/* Precompute all the needed partial products: */\
	__fprod3 = __fx1*__fy2 + __fx2*__fy1 + __fcy;\
	__fprod4 = __fx2*__fy2;\
	\
	__gprod3 = __gx1*__gy2 + __gx2*__gy1 + __gcy;\
	__gprod4 = __gx2*__gy2;\
	\
	/* Digit 3: */\
	__fcy  = NINT(__fprod3*TWO26FLINV);\
	__fhi0 = __fprod3 - __fcy*TWO26FLOAT;\
	\
	__gcy  = NINT(__gprod3*TWO26FLINV);\
	__ghi0 = __gprod3 - __gcy*TWO26FLOAT;\
	\
	/* Digit 4: */\
	__fprod4 += __fcy;\
	__fcy  = NINT(__fprod4*TWO26FLINV);\
	__fhi1 = __fprod4 - __fcy*TWO26FLOAT;\
	\
	__gprod4 += __gcy;\
	__gcy  = NINT(__gprod4*TWO26FLINV);\
	__ghi1 = __gprod4 - __gcy*TWO26FLOAT;\
	\
	/* Digit 5: */\
	__fhi2 = __fcy;\
	__ghi2 = __gcy;\
}

#define MULH78_3WORD_DOUBLE_q4(\
  __fx0,__fx1,__fx2, __fy0,__fy1,__fy2, __fhi0,__fhi1,__fhi2\
, __gx0,__gx1,__gx2, __gy0,__gy1,__gy2, __ghi0,__ghi1,__ghi2\
, __hx0,__hx1,__hx2, __hy0,__hy1,__hy2, __hhi0,__hhi1,__hhi2\
, __ix0,__ix1,__ix2, __iy0,__iy1,__iy2, __ihi0,__ihi1,__ihi2\
  )\
{\
	double __ftmp, __fprod3, __fprod4;\
	double __gtmp, __gprod3, __gprod4;\
	double __htmp, __hprod3, __hprod4;\
	double __itmp, __iprod3, __iprod4;\
	double __fcy, __gcy, __hcy, __icy;\
	\
	/* Digit 0: */\
	__ftmp =  __fx0*__fy0;								__gtmp =  __gx0*__gy0;								__htmp =  __hx0*__hy0;								__itmp =  __ix0*__iy0;\
	__fcy  = (__ftmp*TWO26FLINV);		__gcy  = (__gtmp*TWO26FLINV);		__hcy  = (__htmp*TWO26FLINV);		__icy  = (__itmp*TWO26FLINV);\
	/* Digit 1: */\
	__ftmp = __fx0*__fy1 + __fx1*__fy0 + __fcy;			__gtmp = __gx0*__gy1 + __gx1*__gy0 + __gcy;			__htmp = __hx0*__hy1 + __hx1*__hy0 + __hcy;			__itmp = __ix0*__iy1 + __ix1*__iy0 + __icy;\
	__fcy  = (__ftmp*TWO26FLINV);		__gcy  = (__gtmp*TWO26FLINV);		__hcy  = (__htmp*TWO26FLINV);		__icy  = (__itmp*TWO26FLINV);\
	/* Digit 2: */\
	__ftmp = __fx0*__fy2+__fx1*__fy1+__fx2*__fy0+__fcy;	__gtmp = __gx0*__gy2+__gx1*__gy1+__gx2*__gy0+__gcy;	__htmp = __hx0*__hy2+__hx1*__hy1+__hx2*__hy0+__hcy;	__itmp = __ix0*__iy2+__ix1*__iy1+__ix2*__iy0+__icy;\
	__fcy  = (__ftmp*TWO26FLINV);		__gcy  = (__gtmp*TWO26FLINV+RND_A) - RND_B;			__hcy  = (__htmp*TWO26FLINV+RND_A) - RND_B;			__icy  = (__itmp*TWO26FLINV+RND_A) - RND_B;\
	__ftmp-= __fcy*TWO26FLOAT;							__gtmp-= __gcy*TWO26FLOAT;							__htmp-= __hcy*TWO26FLOAT;							__itmp-= __icy*TWO26FLOAT;\
	__fcy -= (double)(__ftmp < 0);						__gcy -= (double)(__gtmp < 0);						__hcy -= (double)(__htmp < 0);						__icy -= (double)(__itmp < 0);\
	/* At this point the possibility of same-address in-and-outputs comes into play: */\
	/* Precompute all the needed partial products: */\
	__fprod3 = __fx1*__fy2 + __fx2*__fy1 + __fcy;		__gprod3 = __gx1*__gy2 + __gx2*__gy1 + __gcy;		__hprod3 = __hx1*__hy2 + __hx2*__hy1 + __hcy;		__iprod3 = __ix1*__iy2 + __ix2*__iy1 + __icy;\
	__fprod4 = __fx2*__fy2;								__gprod4 = __gx2*__gy2;								__hprod4 = __hx2*__hy2;								__iprod4 = __ix2*__iy2;\
	/* Digit 3: */\
	__fcy  = (__fprod3*TWO26FLINV);		__gcy  = (__gprod3*TWO26FLINV);		__hcy  = (__hprod3*TWO26FLINV);		__icy  = (__iprod3*TWO26FLINV);\
	__fhi0 = __fprod3 - __fcy*TWO26FLOAT;				__ghi0 = __gprod3 - __gcy*TWO26FLOAT;				__hhi0 = __hprod3 - __hcy*TWO26FLOAT;				__ihi0 = __iprod3 - __icy*TWO26FLOAT;\
	/* Digit 4: */\
	__fprod4 += __fcy;									__gprod4 += __gcy;									__hprod4 += __hcy;									__iprod4 += __icy;\
	__fcy  = (__fprod4*TWO26FLINV);		__gcy  = (__gprod4*TWO26FLINV);		__hcy  = (__hprod4*TWO26FLINV);		__icy  = (__iprod4*TWO26FLINV);\
	__fhi1 = __fprod4 - __fcy*TWO26FLOAT;				__ghi1 = __gprod4 - __gcy*TWO26FLOAT;				__hhi1 = __hprod4 - __hcy*TWO26FLOAT;				__ihi1 = __iprod4 - __icy*TWO26FLOAT;\
	/* Digit 5: */\
	__fhi2 = __fcy;										__ghi2 = __gcy;										__hhi2 = __hcy;										__ihi2 = __icy;\
}

#define	CMPLT78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2)	(__x2 < __y2 || (__x2 == __y2 && __x1 < __y1) || (__x2 == __y2 && __x1 == __y1 && __x0 < __y0))

#define	CMPGT78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2)	(__x2 > __y2 || (__x2 == __y2 && __x1 > __y1) || (__x2 == __y2 && __x1 == __y1 && __x0 > __y0))

#define ADD78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2, __z0,__z1,__z2)\
{\
	__z0 = __x0 + __y0;\
	__z1 = __x1 + __y1;\
	__z2 = __x2 + __y2;\
}

#define SUB78_3WORD_DOUBLE(__x0,__x1,__x2, __y0,__y1,__y2, __z0,__z1,__z2)\
{\
	__z0 = __x0 - __y0;\
	__z1 = __x1 - __y1;\
	__z2 = __x2 - __y2;\
}


/***********************************************************************************/
/*** 78-BIT INPUTS, USING FLOATING-POINT-DOUBLE ARITHMETIC FOR MODMUL **************/
/***********************************************************************************/
/*
Function to find 2^(-p) mod q, where p and q are 64 and 78-bit unsigned integers,
respectively, with q stored in a uint96. Uses a floating-double-based Montgomery-style
modmul with a power-of-2 modulus = TWO26FLOAT^3 (i.e. our MODQ operation effects multiply
modulo TWO26FLOAT^3). Using balanced-digit representation in the floating-point form the
largest base we can use with our 3-word-input convolution algorithm is TWO26FLOAT = 2^26,
with the MSW of each input being non-balanced, i.e. strictly nonnegative.

The key 3-operation sequence here is as follows:

	SQR_LOHI78(x,lo,hi);	// Input   x has 78 bits; outputs lo & hi have 78 bits
	MULL78(lo,qinv,lo);		// Inputs lo & qinv, and output (overwrites lo) have 78 bits
	MULH78(q,lo,lo);		// Inputs  q &   lo, and output (overwrites lo) have 78 bits
*/
uint64 twopmodq78_3WORD_DOUBLE(uint64 p, uint96 q)
{
#if FAC_DEBUG
	int dbg = (p > 0);
	uint96 x96, y96, hi;
	uint128 y;
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead7, hi64;
	uint96 qhalf, qinv, x, lo;
	uint192 prod192;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	double fq0, fq1, fq2;
	double fqinv0, fqinv1, fqinv2;
	double fx0, fx1, fx2;
	double flo0,flo1,flo2,fhi0,fhi1,fhi2;

	ASSERT(HERE, (q.d1 >> 14) == 0, "twopmodq78 : (q.d1 >> 14) != 0");

#if FAC_DEBUG
if(dbg)printf("twopmodq78:\n");
	ASSERT(HERE, (TWO26FLOAT * TWO26FLINV) == 1.0, "twopmodq78 : (TWO26FLOAT * TWO26FLINV) != 1.0");
#endif

	/* Convert q to floating form: */
	CVT_UINT78_3WORD_DOUBLE(q, fq0,fq1,fq2);

	RSHIFT_FAST96(q, 1, qhalf);	/* = (q-1)/2, since q odd. */

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 78;

	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 7 or 8 to account for the fact that we can do the powering for the
	!    leftmost 7 or 8 bits (depending on whether the leftmost 7 > 95 or not) via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming
		the leftmost 6/7 bits have already been processed via a shift.

		Since 7 bits with leftmost bit = 1 is guaranteed
		to be in [64,127], the shift count here is in [0, 63].
		That means that zstart < 2^64. Together with the fact that
		squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration.
		Every little bit counts (literally in this case :), right?
		*/
		j = leadz64(pshift);
		/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
		lead7 = ((pshift<<j) >> 57);
		if(lead7 > 77)
		{
			lead7 >>= 1;
			start_index =  64-j-6;	/* Use only the leftmost 6 bits */
		}
		else
			start_index =  64-j-7;

#if FAC_DEBUG
	if(dbg)	printf("lead7 = %u\n", (uint32)lead7);
#endif
		zshift = 77 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

	/*
	!    Find modular inverse (mod 2^78) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q.d0 & (uint64)1) == 1, "twopmodq78 : (q.d0 & (uint64)1) == 1");
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
		hi64 = q.d0*qinv.d0;
		qinv.d0 = qinv.d0*((uint64)2 - hi64);
	}

	/* Now that have bottom 64 bits of qinv, do one more Newton iteration using full 96-bit operands.

	Since qinv.d1 = 0 and q.d0*qinv.d0 == 1 mod 2^64 here; can take advantage of this to speed the iteration:

		q*qinv	= (q.d1*2^64 + q.d0)*qinv.d0
				= q.d1*qinv.d0*2^64 + q.d0*qinv.d0, but this == 1 mod 2^64, i.e.
				= (q.d1*qinv.d0 + MULH64(q.d0,qinv.d0))*2^64 + 1 ,

	i.e. do a MULL64(q.d1,qinv.d0) + MULH64(q.d0,qinv.d0) to get bits 64-127 of q*qinv, set bits 0-63 = 1,
	simply negate bits 64-127 to get r = 2-q*qinv (mod 2^128), then use that (bits 0-63 of r) = 1 to
	simplify the MULL128(qinv,r) operation, i.e.

		qinv*r = qinv.d0*(r.d1*2^64 + 1) == (MULL64(qinv.d0,r.d1)*2^64 + qinv.d0) mod 2^128 .

	In the 96-bit case we then zero the upper 32 bits of the 128-bit result.
	This needs a total of just 3 MULs and 2 ALUs, compared to 8 MULs and 8 ALUs for the naive sequence.
	(May be able to squeeze a bit more out by using 32-bit MULLs rather than 64-bit MULLs in some places,
	but since the main action is in the j-loop below, don't bother with further speedups here.
	*/
#if FAC_DEBUG
	MULL96(q, qinv, x);			MULL128(q, qinv, y);
	SUB96 (TWO96, x, x);		SUB128 (TWO128, y, y);
	MULL96(qinv, x, x);			MULL128(qinv, y, y);
	ASSERT(HERE, x.d1 == (y.d1 & 0x00000000ffffffff) && x.d0 == y.d0, "x.d1 == (y.d1 & 0x00000000ffffffff) && x.d0 == y.d0");
#endif
	/* qinv has 128 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + __MULH64(q.d0, qinv.d0));
#else
	MULH64(q.d0, qinv.d0, hi64);
	qinv.d1 = -qinv.d0*(q.d1*qinv.d0 + hi64);
#endif
	qinv.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */

#if FAC_DEBUG
	ASSERT(HERE, qinv.d1 == (x.d1 & 0x0000000000003fff) && qinv.d0 == x.d0, "twopmodq78 : qinv.d1 == x.d1 && qinv.d0 == x.d0");
	if(dbg) printf("q    = %s\n", &char_buf[convert_uint128_base10_char(char_buf, q   )]);
	if(dbg) printf("qinv = %s\n", &char_buf[convert_uint128_base10_char(char_buf, qinv)]);
#endif

	/* Convert qinv to floating form: */
/*	cvt_uint78_3word_double(qinv, &fqinv0,&fqinv1,&fqinv2);	*/
	CVT_UINT78_3WORD_DOUBLE(qinv, fqinv0,fqinv1,fqinv2);

	/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if FAC_DEBUG
	if(dbg) printf("zshift  = %u\n", zshift);
#endif
	LSHIFT96(qinv, zshift, lo);
	lo.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */

#if FAC_DEBUG
	if(dbg) printf("lo = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	MUL_LOHI96_PROD192(q,lo,prod192);
	RSHIFT192(prod192,78,prod192);
	lo.d0 = prod192.d0;	lo.d1 = prod192.d1;

#if FAC_DEBUG
	if(dbg) printf("q*lo/2^78 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, lo)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	SUB96(q, lo, x);	/* Put the result in lo (rather than x), to ease overflow check below */

#if FAC_DEBUG
	if(dbg) printf("x = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif

	if((pshift >> j) & (uint64)1)
	{
	#if FAC_DEBUG
		ASSERT(HERE, CMPULT96(x, q), "twopmodq80 : CMPULT96(x,q)");
	#endif
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT96(x, qhalf)){ ADD96(x, x, x); SUB96(x, q, x); }else{ ADD96(x, x, x); }
	}

#if FAC_DEBUG
	if(dbg) printf("x0= %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
#endif
#if FAC_DEBUG
	if(CMPULT96(q, x)){ sprintf(char_buf, "twopmodq80 : (x0 = %s) >= (q = %s)", &str0[convert_uint128_base10_char(str0, x)], &str1[convert_uint128_base10_char(str1, q)] );	DBG_WARN(HERE, char_buf, STATFILE, !restart); }
#endif

	/* Convert x to floating form: */
	CVT_UINT78_3WORD_DOUBLE(x, fx0,fx1,fx2);
#if FAC_DEBUG
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x96);
	ASSERT(HERE, CMPEQ96(x, x96), "twopmodq78 : x != fx");
#endif

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		/* SQR_LOHI96(x,lo,hi); */
//SQR_LOHI78_3WORD_DOUBLE(fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2);
		sqr_lohi78_3word_double(fx0,fx1,fx2, &flo0,&flo1,&flo2, &fhi0,&fhi1,&fhi2);
#if FAC_DEBUG
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);		SQR_LOHI96(x,lo,x96);
	/* Split SQR_LOHI96 output into two 78-bit pieces: */
	/* lo has 96 bits but we keep only the low 78 bits of that,
	so upper 18 bits of lo get moved into LSB of hi: */
	hi.d0 = (lo.d1 >> 14) + (x96.d0 << 18);
	lo.d1 &= 0x0000000000003fff;
	/* x96 has at most 2*78 - 96 = 60 bits, so >> (64-18) = 46 gives
	hi.d1 with at most 14 bits, as desired: */
	hi.d1 = (x96.d0 >> 46);
	CVT78_3WORD_DOUBLE_UINT96(flo0,flo1,flo2, x96);
	ASSERT(HERE, CMPEQ96(lo, x96), "twopmodq78A: lo != flo");
	cvt78_3word_double_uint96(fhi0,fhi1,fhi2,&y96);
	ASSERT(HERE, CMPEQ96(hi, y96), "twopmodq78A: hi != fhi");
#endif
		/* MULL96(lo,qinv,lo); */
		MULL78_3WORD_DOUBLE(flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2);
#if FAC_DEBUG
	MULL96(lo,qinv,lo);	lo.d1 &= 0x0000000000003fff;
	CVT78_3WORD_DOUBLE_UINT96(flo0,flo1,flo2, x96);
	ASSERT(HERE, CMPEQ96(lo, x96), "twopmodq78B: lo != flo");
#endif
		/* MULH96(q,lo,lo); */
		MULH78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
#if FAC_DEBUG
	CVT78_3WORD_DOUBLE_UINT96(flo0,flo1,flo2, x96);
	MUL_LOHI96_PROD192(q,lo,prod192);
	RSHIFT192(prod192,78,prod192);
	lo.d0 = prod192.d0;	lo.d1 = prod192.d1;
	ASSERT(HERE, CMPEQ96(lo, x96), "twopmodq78C: lo != flo");
#endif

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
		{
			SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
			ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
		}
		/* Normalize the result: */
		NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

#if FAC_DEBUG
if(dbg)
{
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);
	printf("j = %2d, x = %s",j, &char_buf[convert_uint128_base10_char(char_buf, x)]);
}
#endif

		if((pshift >> j) & (uint64)1)
		{
		#if FAC_DEBUG
			ASSERT(HERE, CMPULT96(x, q), "twopmodq80 : CMPULT96(x,q)");
		#endif
			/* Convert fx to uint96 for purpose of this particular comparison: */
			CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x, qhalf))
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
			}
			/* Normalize the result: */
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

#if FAC_DEBUG
if(dbg)
{
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);
	printf("*2= %s", &char_buf[convert_uint128_base10_char(char_buf, x)]);
}
#endif
		}
#if FAC_DEBUG
if(dbg)
{
	printf("\n");
}
#endif
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x);
	ADD96(x,x,x);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
	SUB96(x,q,x);

#if FAC_DEBUG
if(dbg)
{	printf("x0 = %s\n", &char_buf[convert_uint128_base10_char(char_buf, x)]);
	x96 = twopmodq96(p, q);
	ASSERT(HERE, CMPEQ96(x, x96), "twopmodq80 : Result differs from that returned by twopmodq96!");
}
#endif

	return (CMPEQ96(x, ONE96));
}

/*
This routine needs the SSE rounding mode to be set to round-toward-minus-infinity, so we need to twiddle the MXCSR register -
from http://softpixel.com/~cwright/programming/simd/sse.php :

The MXCSR register is a 32-bit register containing flags for control and status information regarding SSE instructions.
As of SSE3, only bits 0-15 have been defined.

Mnemonic	Bit Location	Description				[EWM: Default value on MSVC/ia32 = 0x1FA0, so the bits marked with [x] are set:]
--------	------------	---------------------	-----
	FZ		bit 15			Flush To Zero
	R+		bit<14:13> = 10	Round Positive
	R-		bit<14:13> = 01	Round Negative
	RZ		bit<14:13> = 11	Round To Zero
	RN		bit<14:13> = 00	Round To Nearest		[x]
	PM		bit 12			Precision Mask			[x]
	UM		bit 11			Underflow Mask			[x]
	OM		bit 10			Overflow Mask			[x]
	ZM		bit 9			Divide By Zero Mask		[x]
	DM		bit 8			Denormal Mask			[x]
	IM		bit 7			Invalid Operation Mask	[x]
	DAZ		bit 6			Denormals Are Zero
	PE		bit 5			Precision Flag			[x]
	UE		bit 4			Underflow Flag
	OE		bit 3			Overflow Flag
	ZE		bit 2			Divide By Zero Flag
	DE		bit 1			Denormal Flag
	IE		bit 0			Invalid Operation Flag

So to enable round-toward-minus-infinity, we set bit 13.

The problem with this is that we need to do it and undo it every pass through the squaring sequence:

	__asm	stmxcsr	mxcsr_ptr
	__asm	xor		mxcsr_ptr, 0x2000
	__asm	ldmxcsr	mxcsr_ptr
		__asm	movaps	xmm6,xmm2	// fcy = cpy of fprod2
		__asm cvtpd2dq	xmm6,xmm6	// Convert to 32-bit int, rounding toward -oo
		__asm cvtdq2pd	xmm6,xmm6	// ...and back to double.
	__asm	xor		mxcsr_ptr, 0x2000
	__asm	ldmxcsr	mxcsr_ptr

So we instead emulate [round toward -oo](x) via NINT(x-(0.5-epsilon)).
*/

/*** 2-trial-factor version ***/
uint64 twopmodq78_3WORD_DOUBLE_q2(uint64 p, uint96 q0, uint96 q1)
{
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead7, tmp0, tmp1, r;
	uint96 qinv0, qinv1
		, qhalf0, qhalf1
		, x0, x1
		, lo0, lo1;
	uint192 prod192;
	static uint64 psave = 0;
	static uint32 pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

#ifdef USE_SSE2

	static struct complex *sc_arr = 0x0;
	static         double *sc_ptr;
	static double *fq0,*fq1,*fq2,*fqhi52, *fqinv0,*fqinv1,*fqinv2, *fx0,*fx1,*fx2, *flo0,*flo1,*flo2, *fhi0,*fhi1,*fhi2;
	static double *gq0,*gq1,*gq2,*gqhi52, *gqinv0,*gqinv1,*gqinv2, *gx0,*gx1,*gx2, *glo0,*glo1,*glo2, *ghi0,*ghi1,*ghi2;
	static double *half, *two26f, *two26i, *two13i, *sse2_rnd;
	static uint64 ihalf = 0x3FDfffffffffffffull;	/* Bitfield storing 0.5*(1-epsilon) in IEEE64 format */

#else

	double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2, fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2;
	double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2;

#endif

#ifdef USE_SSE2

	if(first_entry)
	{
		sc_arr = ALLOC_COMPLEX(sc_arr, 24);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = (double *)ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

		fq0    = sc_ptr + 0x00;		gq0    = sc_ptr + 0x01;
		fq1    = sc_ptr + 0x02;		gq1    = sc_ptr + 0x03;
		fq2    = sc_ptr + 0x04;		gq2    = sc_ptr + 0x05;
		fqhi52 = sc_ptr + 0x06;		gqhi52 = sc_ptr + 0x07;
		fqinv0 = sc_ptr + 0x08;		gqinv0 = sc_ptr + 0x09;
		fqinv1 = sc_ptr + 0x0a;		gqinv1 = sc_ptr + 0x0b;
		fqinv2 = sc_ptr + 0x0c;		gqinv2 = sc_ptr + 0x0d;
		fx0    = sc_ptr + 0x0e;		gx0    = sc_ptr + 0x0f;
		fx1    = sc_ptr + 0x10;		gx1    = sc_ptr + 0x11;
		fx2    = sc_ptr + 0x12;		gx2    = sc_ptr + 0x13;
		flo0   = sc_ptr + 0x14;		glo0   = sc_ptr + 0x15;
		flo1   = sc_ptr + 0x16;		glo1   = sc_ptr + 0x17;
		flo2   = sc_ptr + 0x18;		glo2   = sc_ptr + 0x19;
		fhi0   = sc_ptr + 0x1a;		ghi0   = sc_ptr + 0x1b;
		fhi1   = sc_ptr + 0x1c;		ghi1   = sc_ptr + 0x1d;
		fhi2   = sc_ptr + 0x1e;		ghi2   = sc_ptr + 0x1f;
		two13i = sc_ptr + 0x20;
		two26f = sc_ptr + 0x22;
		two26i = sc_ptr + 0x24;
		sse2_rnd=sc_ptr + 0x26;
		half   = sc_ptr + 0x28;
		/* Can premultiply each of the multiword-mul inputs by 1/sqrt(2^13) due to quadraticity: */
		*two13i++ = TWO13FLINV;		*two13i-- = TWO13FLINV;
		*two26f++ = TWO26FLOAT;		*two26f-- = TWO26FLOAT;
		*two26i++ = TWO26FLINV;		*two26i-- = TWO26FLINV;
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		*sse2_rnd++ = 3.0*0x4000000*0x2000000;
		*sse2_rnd-- = 3.0*0x4000000*0x2000000;
		/* We init "half" = 0.5-epsilon here, because emulating FLOOR(x) via NINT(x-half) requires
		us to always round up if x is a whole number, and our NINT emulation can round either way if fractional part = 0.5:
		*/
		*half++ = *(double*)&ihalf;	*half-- = *(double*)&ihalf;
	}	/* first_entry */
#endif

	ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q2 : (q0.d1 >> 14) != 0");
	ASSERT(HERE, (q1.d1 >> 14) == 0, "twopmodq78_q2 : (q1.d1 >> 14) != 0");

	/* Convert q to floating form: */
#ifdef USE_SSE2

	CVT_UINT78_3WORD_DOUBLE(q0 ,*fq0,*fq1,*fq2);
	CVT_UINT78_3WORD_DOUBLE(q1 ,*gq0,*gq1,*gq2);

	#if defined(COMPILER_TYPE_MSVC)

		/* Scale by TWO26FLINV: */
		__asm	mov	eax, fq0
		__asm	mov	ebx, two26f
		__asm	movaps	xmm5,[ebx      ]	/* two26f */
		__asm	movaps	xmm6,[ebx+0x010]	/* two26i */
		__asm	movaps	xmm0,[eax      ]	/* fq0 */
		__asm	movaps	xmm1,[eax+0x010]	/* fq1 */
		__asm	movaps	xmm2,[eax+0x020]	/* fq2 */
		__asm	movaps	xmm3,xmm2		/* cpy fq2 */
		__asm	mulpd	xmm3,xmm5
		__asm	addpd	xmm3,xmm1		/* Hi 52 bits of q */
		__asm	mulpd	xmm0,xmm6
		__asm	mulpd	xmm1,xmm6
		__asm	mulpd	xmm2,xmm6
		__asm	movaps	[eax      ],xmm0	/* fq0/2^26 */
		__asm	movaps	[eax+0x010],xmm1	/* fq1/2^26 */
		__asm	movaps	[eax+0x020],xmm2	/* fq2/2^26 */
		__asm	movaps	[eax+0x030],xmm3	/* fqhi52 */

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			__asm__ volatile (\
				"movl	%[__fq0],%%eax		\n\t"\
				"movl	%[__two26f],%%esi	\n\t"\
				"movaps	     (%%esi),%%xmm5	\n\t"\
				"movaps	0x010(%%esi),%%xmm6	\n\t"\
				"movaps	     (%%eax),%%xmm0	\n\t"\
				"movaps	0x010(%%eax),%%xmm1	\n\t"\
				"movaps	0x020(%%eax),%%xmm2	\n\t"\
				"movaps	%%xmm2,%%xmm3		\n\t"\
				"mulpd	%%xmm5,%%xmm3		\n\t"\
				"addpd	%%xmm1,%%xmm3		\n\t"\
				"mulpd	%%xmm6,%%xmm0		\n\t"\
				"mulpd	%%xmm6,%%xmm1		\n\t"\
				"mulpd	%%xmm6,%%xmm2		\n\t"\
				"movaps	%%xmm0,     (%%eax)	\n\t"\
				"movaps	%%xmm1,0x010(%%eax)	\n\t"\
				"movaps	%%xmm2,0x020(%%eax)	\n\t"\
				"movaps	%%xmm3,0x030(%%eax)	\n\t"\
				:					/* outputs: none */\
				: [__fq0] "m" (fq0)	/* All inputs from memory addresses here */\
				,[__two26f] "m" (two26f)\
				: "eax","esi"		/* Clobbered registers */\
			);

		#else

			__asm__ volatile (\
				"movq	%[__fq0],%%rax 	\n\t"\
				"movq	%[__two26f],%%rsi \n\t"\
				"movaps	     (%%rsi),%%xmm5	\n\t"\
				"movaps	0x010(%%rsi),%%xmm6	\n\t"\
				"movaps	     (%%rax),%%xmm0	\n\t"\
				"movaps	0x010(%%rax),%%xmm1	\n\t"\
				"movaps	0x020(%%rax),%%xmm2	\n\t"\
				"movaps	%%xmm2,%%xmm3		\n\t"\
				"mulpd	%%xmm5,%%xmm3		\n\t"\
				"addpd	%%xmm1,%%xmm3		\n\t"\
				"mulpd	%%xmm6,%%xmm0		\n\t"\
				"mulpd	%%xmm6,%%xmm1		\n\t"\
				"mulpd	%%xmm6,%%xmm2		\n\t"\
				"movaps	%%xmm0,     (%%rax)	\n\t"\
				"movaps	%%xmm1,0x010(%%rax)	\n\t"\
				"movaps	%%xmm2,0x020(%%rax)	\n\t"\
				"movaps	%%xmm3,0x030(%%rax)	\n\t"\
				:					/* outputs: none */\
				: [__fq0] "m" (fq0)	/* All inputs from memory addresses here */\
				 ,[__two26f] "m" (two26f)\
				: "rax","rsi","xmm0","xmm1","xmm2","xmm3","xmm5","xmm6"		/* Clobbered registers */\
			);

		#endif

	#endif

#else
	CVT_UINT78_3WORD_DOUBLE(q0, fq0,fq1,fq2);
	CVT_UINT78_3WORD_DOUBLE(q1, gq0,gq1,gq2);
#endif

	RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST96(q1, 1, qhalf1);

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = (uint32)p + 78;

	/*
	!    find number of leading zeros in p, use it to find the position of the leftmost
	!    ones bit, and subtract 7 or 8 to account for the fact that we can do the powering for the
	!    leftmost 7 or 8 bits (depending on whether the leftmost 7 > 95 or not) via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming the leftmost 6/7 bits have already been processed via a shift.

		Since 7 bits with leftmost bit = 1 is guaranteed to be in [64,127], the shift count here is in [0, 63].
		That means that zstart < 2^64. Together with the fact that squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration. Every little bit counts (literally in this case :), right?
		*/
		j = leadz32(pshift);
		/* Extract leftmost 7 bits of pshift (if > 77, use the leftmost 6) and subtract from 96: */
		lead7 = ((pshift<<j) >> 25);
		if(lead7 > 77)
		{
			lead7 >>= 1;
			start_index =  32-j-6;	/* Use only the leftmost 6 bits */
		}
		else
			start_index =  32-j-7;

		zshift = 77 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

	/*
	!    Find modular inverse (mod 2^78) of q in preparation for modular multiply.
	*/
	/* q must be odd for Montgomery-style modmul to work: */
	/* Init qinv = q. We're really only interested in the bottom 2 bits of q. */
	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*(2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
	defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
	*/
	for(j = 0; j < 4; j++)
	{
		tmp0 = q0.d0*qinv0.d0;
		tmp1 = q1.d0*qinv1.d0;

		qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
		qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
	}

	/* Now that have bottom 64 bits of qinv, do one more Newton iteration using full 96-bit operands.
	qinv has 128 bits, but only the upper 64 get modified here. */
#ifdef MUL_LOHI64_SUBROUTINE
	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
#else
	MULH64(q0.d0, qinv0.d0, tmp0);
	MULH64(q1.d0, qinv1.d0, tmp1);

	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);
#endif
	qinv0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
	qinv1.d1 &= 0x0000000000003fff;

	/* Convert qinv to floating form: */
#ifdef USE_SSE2
	CVT_UINT78_3WORD_DOUBLE(qinv0 ,*fqinv0,*fqinv1,*fqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv1 ,*gqinv0,*gqinv1,*gqinv2);

	#if defined(COMPILER_TYPE_MSVC)

		/* Scale by TWO26FLINV: */
		__asm	mov	eax, fqinv0
		__asm	mov	ebx, two26i
		__asm	movaps	xmm6,[ebx]	/* two26i */
		__asm	movaps	xmm0,[eax      ]	/* fqinv0 */
		__asm	movaps	xmm1,[eax+0x010]	/* fqinv1 */
		__asm	movaps	xmm2,[eax+0x020]	/* fqinv2 */
		__asm	mulpd	xmm0,xmm6
		__asm	mulpd	xmm1,xmm6
		__asm	mulpd	xmm2,xmm6
		__asm	movaps	[eax      ],xmm0	/* fqinv0 */
		__asm	movaps	[eax+0x010],xmm1	/* fqinv1 */
		__asm	movaps	[eax+0x020],xmm2	/* fqinv2 */

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

		__asm__ volatile (\
			"movl	%[__fqinv0],%%eax 	\n\t"\
			"movl	%[__two26i],%%esi 	\n\t"\
			"movaps	     (%%esi),%%xmm6	\n\t"\
			"movaps	     (%%eax),%%xmm0	\n\t"\
			"movaps	0x010(%%eax),%%xmm1	\n\t"\
			"movaps	0x020(%%eax),%%xmm2	\n\t"\
			"mulpd	%%xmm6,%%xmm0		\n\t"\
			"mulpd	%%xmm6,%%xmm1		\n\t"\
			"mulpd	%%xmm6,%%xmm2		\n\t"\
			"movaps	%%xmm0,     (%%eax)	\n\t"\
			"movaps	%%xmm1,0x010(%%eax)	\n\t"\
			"movaps	%%xmm2,0x020(%%eax)	\n\t"\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			: "eax","esi"	/* Clobbered registers */\
		);

		#else

		__asm__ volatile (\
			"movq	%[__fqinv0],%%rax 	\n\t"\
			"movq	%[__two26i],%%rsi 	\n\t"\
			"movaps	     (%%rsi),%%xmm6	\n\t"\
			"movaps	     (%%rax),%%xmm0	\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t"\
			"movaps	0x020(%%rax),%%xmm2	\n\t"\
			"mulpd	%%xmm6,%%xmm0		\n\t"\
			"mulpd	%%xmm6,%%xmm1		\n\t"\
			"mulpd	%%xmm6,%%xmm2		\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t"\
			"movaps	%%xmm2,0x020(%%rax)	\n\t"\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			,[__two26i] "m" (two26i)	\
			: "rax","rsi","xmm0","xmm1","xmm2","xmm6"	/* Clobbered registers */\
		);

		#endif

	#endif

  #ifdef DEBUG_SSE2
	fprintf(stderr, "[f|g]qinv0 = %20.0f, %20.0f\n",*fqinv0*TWO26FLOAT,*gqinv0*TWO26FLOAT);
	fprintf(stderr, "[f|g]qinv1 = %20.0f, %20.0f\n",*fqinv1*TWO26FLOAT,*gqinv1*TWO26FLOAT);
	fprintf(stderr, "[f|g]qinv2 = %20.0f, %20.0f\n",*fqinv2*TWO26FLOAT,*gqinv2*TWO26FLOAT);
	fprintf(stderr, "\n");
	/******************************************************/
  #endif
#else
	CVT_UINT78_3WORD_DOUBLE(qinv0, fqinv0,fqinv1,fqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv1, gqinv0,gqinv1,gqinv2);
#endif

	/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT96(qinv0, zshift, lo0);	lo0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
	LSHIFT96(qinv1, zshift, lo1);	lo1.d1 &= 0x0000000000003fff;

	MUL_LOHI96_PROD192(q0,lo0,prod192);	RSHIFT192(prod192,78,prod192);	lo0.d0 = prod192.d0;	lo0.d1 = prod192.d1;
	MUL_LOHI96_PROD192(q1,lo1,prod192);	RSHIFT192(prod192,78,prod192);	lo1.d0 = prod192.d0;	lo1.d1 = prod192.d1;

	/* hi = 0 in this instance, which simplifies things. */
	SUB96(q0, lo0, x0);	/* Put the result in lo (rather than x), to ease overflow check below */
	SUB96(q1, lo1, x1);

	if((pshift >> j) & (uint32)1)
	{
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
		if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
	}

	/* Convert x to floating form: */
#ifdef USE_SSE2

	CVT_UINT78_3WORD_DOUBLE(x0 ,*fx0,*fx1,*fx2);
	CVT_UINT78_3WORD_DOUBLE(x1 ,*gx0,*gx1,*gx2);

	#if defined(COMPILER_TYPE_MSVC)

		__asm	mov	eax, fx0
		__asm	movaps	xmm0,[eax      ]	/* fx0 */
		__asm	movaps	xmm2,[eax+0x010]	/* fx1 */
		__asm	movaps	xmm4,[eax+0x020]	/* fx2 */

	#else	/* GCC-style inline ASM: */

		/* GCC builds require us to save the xmm0/2/4 registers on each pass thru the loop,
		so skip these pre-loop-entry snippets until we get the full loop-logic asm'ed: */
		#if 0
		  #if OS_BITS == 32
			__asm__ volatile (\
				"movl	%[__fx0],%%eax 	\n\t"\
				"movaps	     (%%eax),%%xmm0	\n\t"\
				"movaps	0x010(%%eax),%%xmm2	\n\t"\
				"movaps	0x020(%%eax),%%xmm4	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "eax"	/* Clobbered registers */\
			);
		  #else
			__asm__ volatile (\
				"movq	%[__fx0],%%rax 	\n\t"\
				"movaps	     (%%rax),%%xmm0	\n\t"\
				"movaps	0x010(%%rax),%%xmm2	\n\t"\
				"movaps	0x020(%%rax),%%xmm4	\n\t"\
				:				/* outputs: none */\
				: [__fx0] "m" (fx0)		/* All inputs from memory addresses here */\
				: "rax","xmm0","xmm2","xmm4"	/* Clobbered registers */\
			);
		  #endif
		#endif
	#endif

#else

	CVT_UINT78_3WORD_DOUBLE(x0, fx0,fx1,fx2);
	CVT_UINT78_3WORD_DOUBLE(x1, gx0,gx1,gx2);

#endif

	/*...x^2 mod q is returned in x. */
	/* All 3-word-double-form operands have components in the following size ranges:
		fword0,1 in [-2^25, +2^25]
		fword2   in [   -1, +2^26]
	*/
#ifndef USE_SSE2

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		/* SQR_LOHI96(x,lo,hi); */
		SQR_LOHI78_3WORD_DOUBLE_q2(
			  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2
			, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2
		);
		/* MULL96(lo,qinv,lo); */
		MULL78_3WORD_DOUBLE_q2(
			  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
			, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
		);
		/* MULH96(q,lo,lo); */
		MULH78_3WORD_DOUBLE_q2(
			  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2
			, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2
		);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
		{
			SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
			ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
		}
		/* Normalize the result: */
		NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

		if(CMPLT78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2))
		{
			SUB78_3WORD_DOUBLE(gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2);
			ADD78_3WORD_DOUBLE(glo0,glo1,glo2, ghi0,ghi1,ghi2, gx0,gx1,gx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2, gx0,gx1,gx2);
		}
		/* Normalize the result: */
		NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

		if((pshift >> j) & (uint64)1)
		{
			/* Convert fx to uint96 for purpose of this particular comparison: */
			CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
			CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x0, qhalf0))
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
			}
			/* Normalize the result: */
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x1, qhalf1))
			{
				ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
				SUB78_3WORD_DOUBLE(gx0,gx1,gx2, gq0,gq1,gq2, gx0,gx1,gx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
			}
			/* Normalize the result: */
			NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);
		}

	#ifdef DEBUG_THIS
		fprintf(stderr, "MULH78_3WORD_DOUBLE_q2: [f|g]x0 = %20.5f, %20.5f\n",*fx0,*gx0);
		fprintf(stderr, "MULH78_3WORD_DOUBLE_q2: [f|g]x1 = %20.5f, %20.5f\n",*fx1,*gx1);
		fprintf(stderr, "MULH78_3WORD_DOUBLE_q2: [f|g]x2 = %20.5f, %20.5f\n",*fx2,*gx2);
		fprintf(stderr, "\n");
	#endif

	}	/* for(j...) */

#elif(defined(USE_SSE2))	/* Inner loop body needs 42 movaps, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */

	#if 0	/* Set = 1 to do timings of all-but-inner-loop */
	#else	/* Inner loop body needs 42 movaps, 76 ADD/SUBPD, 52 MULPD, 13 MISC/ALU (ANDPD, XORPD, CMPPD, etc) */

		#if defined(COMPILER_TYPE_MSVC)

		__asm	mov ecx, start_index
		__asm	sub ecx, 2
		__asm	test ecx, ecx	/* int tmp = esi & esi = n*/
		__asm	jl LoopEnd		/* Skip if n < 0 */
		LoopBeg:
		/*
		SQR_LOHI78_3WORD_DOUBLE_q2(
			  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2
			, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2
		);
		*/
			__asm	mov	eax, fx0
			__asm	mov	ebx, two26i
			__asm	movaps	xmm6,[ebx-0x020]	/* two13i */
			/* fx0,1,2 assumed in xmm0,2,4 on loop entry */
			__asm	mulpd	xmm0,xmm6	/* scale [fx0 . 2^-13] */
			__asm	mulpd	xmm2,xmm6	/* scale [fx1 . 2^-13] */
			__asm	mulpd	xmm4,xmm6	/* scale [fx2 . 2^-13] */
			__asm	movaps	xmm7,[ebx+0x10]	/* xmm7 = RND */
			__asm	movaps	xmm1,xmm0	/* cpy of fx0 */
			__asm	movaps	xmm3,xmm2	/* cpy of fx1 */
			__asm	addpd	xmm1,xmm1	/* 2.fx0 */
			__asm	addpd	xmm3,xmm3	/* 2.fx1 */
			__asm	movaps	xmm5,xmm1	/* cpy of 2.fx0 */
			__asm	mulpd	xmm0,xmm0	/* [  fx0 *= fx0] / 2^26 */
			__asm	mulpd	xmm1,xmm2	/* [2.fx0 *= fx1] / 2^26 */
			__asm	mulpd	xmm2,xmm2	/* [  fx1 *= fx1] / 2^26 */
			__asm	mulpd	xmm5,xmm4	/* [2.fx0 *= fx2] / 2^26 */
			__asm	mulpd	xmm3,xmm4	/* [2.fx1 *= fx2] / 2^26 */
			__asm	mulpd	xmm4,xmm4	/* [  fx2 *= fx2] / 2^26 */
			/* Digit 0:
			__fprod0  =  __fx0*__fx0;		in [0, 2^50]
			__fcy     = NINT(__fprod0*TWO26FLINV);
			__fprod0 -= __fcy*TWO26FLOAT;	in [-2^24, +2^24]
			*/
			__asm	movaps	xmm6,xmm0	/* Init: fcy = cpy of fx0*fx0 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm0,xmm6	/* fprod0 -= fcy */
			__asm	mulpd	xmm0,[ebx-0x10]	/* fprod0 *= two26f */
			__asm	mulpd	xmm6,[ebx     ]	/* fcy    *= two26i */

			/* Digit 1:
			__fprod1  = __f2x0*__fx1 + __fcy;	in [-2^51-2^24, +2^51+2^24]
			__fcy     = NINT(__fprod1*TWO26FLINV);
			__fprod1 -= __fcy*TWO26FLOAT;	in [-2^25, +2^25]
			*/
			__asm	addpd	xmm1,xmm6	/* fprod1 = 2.fx0*fx1 + fcy */
			__asm	movaps	xmm6,xmm1	/* fcy = cpy of fprod1 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm1,xmm6	/* fprod1 -= fcy */
			__asm	mulpd	xmm1,[ebx-0x10]	/* fprod1 *= two26f */
			__asm	mulpd	xmm6,[ebx     ]	/* fcy    *= two26i */

			/* Digit 2:
			__fprod2  = __f2x0*__fx2 + __fx1*__fx1 + __fcy;	in [-2^52-2^25, +2^52+2^25]
			__fcy     = NINT(__fprod2*TWO26FLINV);
			__fprod2 -= __fcy*TWO26FLOAT;
			__itmp    = (__fprod2 < 0);
			__fcy    -= (double)__itmp;
			__fprod2 += (double)(__itmp << 26);
			*/
			__asm	addpd	xmm2,xmm5	/* fx1*fx1 + 2.fx0*fx2,		xmm5 FREE */
			__asm	movaps	xmm5,[ebx-0x10]	/* 2^26 */
			__asm	addpd	xmm2,xmm6	/* fprod2 = 2.fx0*fx2 + fx1*fx1 + fcy */
			__asm	movaps	xmm6,xmm2	/* fcy = cpy of fprod2 */
			__asm	subpd	xmm6,[ebx+0x020]	/* fprod2 - 0.5 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy = floor(fprod2) */
			__asm	subpd	xmm2,xmm6	/* fprod2 -= fcy */
			__asm	mulpd	xmm2,xmm5	/* fprod2 *= two26f */
			__asm	mulpd	xmm6,[ebx]	/* fcy    *= two26i */

			/* Digit 3:
			__fprod3  = __f2x1*__fx2 + __fcy;
			__fcy     = NINT(__fprod3*TWO26FLINV);
			__fprod3 -= __fcy*TWO26FLOAT;
			*/
			__asm	addpd	xmm3,xmm6	/* fprod3 = 2.fx1*fx2 + fcy */
			__asm	movaps	xmm6,xmm3	/* fcy = cpy of fprod3 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm3,xmm6	/* fprod3 -= fcy */
			__asm	mulpd	xmm3,xmm5	/* fprod3 *= two26f */
			__asm	mulpd	xmm6,[ebx]	/* fcy    *= two26i */

			/* Digit 4:
			__fprod4  =  __fx2*__fx2 + __fcy;
			__fcy     = NINT(__fprod4*TWO26FLINV);
			__fprod4 -= __fcy*TWO26FLOAT;
			*/
			__asm	addpd	xmm4,xmm6	/* fprod4 = fx2*fx2 + fcy */
			__asm	movaps	xmm6,xmm4	/* fcy = cpy of fprod4 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm4,xmm6	/* fprod4 -= fcy */
			__asm	mulpd	xmm4,xmm5	/* fprod4 *= two26f */

			/* Digit 5 = the carry.
			flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6
			*/
			__asm	mulpd	xmm6,xmm5			/* fhi2 * two26f */
			__asm	addpd	xmm4,xmm6			/* fhi, top 52 bits */
			__asm	movaps	[eax+0x060],xmm3	/* Store fhi0,1,2 to free up 3 registers */
			__asm	movaps	[eax+0x070],xmm4

		/* MULL78_3WORD_DOUBLE_q2(
			  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
			, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2);
		*/
			/* Precompute all the needed partial products: x0,1,2 = flo0,1,2; y0,1,2 = fqinv0,1,2:
			__prod0  = __x0*__y0;
			__prod1  = __x0*__y1 + __x1*__y0;
			__prod2  = __x0*__y2 + __x1*__y1 + __x2*__y0;
			*/
			__asm	mov	edx, fqinv0

			/* Digit 0:
			__fcy    = NINT(__prod0*TWO26FLINV);
			*__lo0   = __prod0 - __fcy*TWO26FLOAT;
			*/
			__asm	movaps	xmm3,xmm0		/* xmm3 = cpy of x0 */
			__asm	movaps	xmm4,xmm0		/* xmm4 = cpy of x0 */
			__asm	mulpd	xmm0,[edx     ]	/* fprod0 = x0*y0 */
			__asm	movaps	xmm6,xmm0	/* fcy = cpy of fprod0 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm0,xmm6	/* fprod0 -= fcy */
			__asm	mulpd	xmm0,xmm5	/* fprod0 *= two26f */
			__asm	mulpd	xmm6,[ebx]	/* fcy    *= two26i */

			/* Digit 1:
			__prod1 += __fcy;
			__fcy    = NINT(__prod1*TWO26FLINV);
			*__lo1   = __prod1 - __fcy*TWO26FLOAT;
			*/
			__asm	mulpd	xmm3,[edx+0x10]	/* x0*y1 */
			__asm	addpd	xmm3,xmm6		/* x0*y1 + fcy; xmm6 FREE */
			__asm	movaps	xmm6,xmm1		/* xmm6 = cpy of x1 */
			__asm	mulpd	xmm1,[edx     ]	/* x1*y0 */
			__asm	addpd	xmm1,xmm3		/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */
			__asm	movaps	xmm3,xmm1		/* fcy = cpy of fprod1 */
			__asm	addpd	xmm3,xmm7
			__asm	subpd	xmm3,xmm7	/* fcy */
			__asm	subpd	xmm1,xmm3	/* fprod1 -= fcy */
			__asm	mulpd	xmm1,xmm5	/* fprod1 *= two26f */
			__asm	mulpd	xmm3,[ebx]	/* fcy    *= two26i */

			/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced:
			__prod2 += __fcy;
			__fcy    = NINT(__prod2*TWO26FLINV);
			*__lo2   = __prod2 - __fcy*TWO26FLOAT;
			*/
			__asm	mulpd	xmm2,[edx     ]	/* x2*y0 */
			__asm	mulpd	xmm6,[edx+0x10]	/* x1*y1 */
			__asm	mulpd	xmm4,[edx+0x20]	/* x0*y2 */
			__asm	addpd	xmm2,xmm6		/* x1*y1 + x2*y0; xmm6 FREE */
			__asm	addpd	xmm3,xmm4		/* x0*y2 + fcy; xmm4 FREE */
			__asm	addpd	xmm2,xmm3		/* fprod2; xmm3 FREE */
			__asm	movaps	xmm3,xmm2		/* fcy = cpy of fprod2 */
			__asm	subpd	xmm3,[ebx+0x020]	/* fprod2 - 0.5 */
			__asm	addpd	xmm3,xmm7
			__asm	subpd	xmm3,xmm7	/* fcy */
			__asm	subpd	xmm2,xmm3	/* fprod2 -= fcy */
			__asm	mulpd	xmm2,xmm5	/* fprod2 *= two26f */

			/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78
			MULH78_3WORD_DOUBLE_q2(
				  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2
				, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2 );
			*/
			/* Precompute all the needed partial products: x0,1,2 = (fq0,1,2)/2^26 ; y0,1,2 = flo0,1,2: */
			__asm	mov	edx, fq0

			/* Digit 0:
			__ftmp =  __fx0*__fy0;
			__fcy  = NINT(__ftmp*TWO26FLINV);
			*/
			__asm	movaps	xmm3,xmm0	/* xmm3 = cpy of y0 */
			__asm	movaps	xmm4,xmm0	/* xmm4 = cpy of y0 */
			__asm	mulpd	xmm0,[edx]	/* fprod0 = y0*x0 */
			__asm	addpd	xmm0,xmm7
			__asm	subpd	xmm0,xmm7	/* fcy */
			__asm	mulpd	xmm0,[ebx]	/* fcy    *= two26i */

			/* Digit 1:
			__ftmp = __fx0*__fy1 + __fx1*__fy0 + __fcy;
			__fcy  = NINT(__ftmp*TWO26FLINV);
			*/
			__asm	movaps	xmm6,xmm1		/* xmm6 = cpy of y1 */
			__asm	mulpd	xmm3,[edx+0x10]	/* y0*x1 */
			__asm	addpd	xmm3,xmm0		/* y0*x1 + fcy; xmm0 FREE */
			__asm	mulpd	xmm1,[edx     ]	/* y1*x0 */
			__asm	addpd	xmm1,xmm3		/* fprod1 = x0*y1 + x1*y0; xmm3 FREE */
			__asm	addpd	xmm1,xmm7
			__asm	subpd	xmm1,xmm7	/* fcy */
			__asm	mulpd	xmm1,[ebx]	/* fcy    *= two26i */

			/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced:
			__ftmp = __fx0*__fy2 + __fx1*__fy1 + __fx2*__fy0 + __fcy;
			__fcy  = NINT(__ftmp*TWO26FLINV);
			__ftmp-= __fcy*TWO26FLOAT;
			__itmp = (__ftmp < 0);
			__fcy -= (double)__itmp;
			*/
			__asm	movaps	xmm3,xmm2		/* xmm3 = cpy of y2 */
			__asm	movaps	xmm0,xmm6		/* xmm0 = cpy of y1 */
			__asm	mulpd	xmm2,[edx     ]	/* y2*x0 */
			__asm	mulpd	xmm6,[edx+0x10]	/* y1*x1 */
			__asm	mulpd	xmm4,[edx+0x20]	/* y0*x2 */
			__asm	addpd	xmm2,xmm6		/* y1*x1 + y2*x0; xmm6 FREE */
			__asm	addpd	xmm1,xmm4		/* y0*x2 + fcy; xmm4 FREE */
			__asm	addpd	xmm2,xmm1		/* fprod2; xmm1 FREE */
			__asm	subpd	xmm2,[ebx+0x020]	/* fprod2 - 0.5 */
			__asm	addpd	xmm2,xmm7
			__asm	subpd	xmm2,xmm7	/* fcy */
			__asm	mulpd	xmm2,[ebx]	/* fcy    *= two26i */

		/* xmm1,3,4,6 FREE */
			/* Precompute all the needed partial products: */

			__asm	movaps	xmm1,xmm3		/* xmm1 = cpy of y2 */
			__asm	mulpd	xmm0,[edx+0x20]	/* y1*x2 */
			__asm	mulpd	xmm1,[edx+0x10]	/* y2*x1 */
			__asm	mulpd	xmm3,[edx+0x20]	/* y2*x2 */

		/* xmm4,6 FREE */
			/* Digit 3:
			__fprod3 = __fx1*__fy2 + __fx2*__fy1 + __fcy;
			__fcy  = NINT(__fprod3*TWO26FLINV);
			__fhi0 = __fprod3 - __fcy*TWO26FLOAT;
			*/
			__asm	addpd	xmm1,xmm2	/* fy2*fx1 + fcy */
			__asm	addpd	xmm0,xmm1	/* fprod3 = fy1*fx2 + fy2*fx1 + fcy */
			__asm	movaps	xmm6,xmm0	/* fcy = cpy of fprod3 */
			__asm	addpd	xmm6,xmm7
			__asm	subpd	xmm6,xmm7	/* fcy */
			__asm	subpd	xmm0,xmm6	/* fprod3 -= fcy */
			__asm	mulpd	xmm0,xmm5	/* fprod3 *= two26f */
			__asm	mulpd	xmm6,[ebx]	/* fcy    *= two26i */

		/* xmm1,2,4 FREE */
			/* Digit 4:
			__fprod4 = __fx2*__fy2;
			__fprod4 += __fcy;
			__fcy  = NINT(__fprod4*TWO26FLINV);
			__fhi1 = __fprod4 - __fcy*TWO26FLOAT;
			*/
			__asm	addpd	xmm3,xmm6	/* fcy = fy2*fx2 + fcy */
			__asm	movaps	xmm1,xmm3	/* fprod4 = cpy of fcy */
			__asm	addpd	xmm3,xmm7
			__asm	subpd	xmm3,xmm7	/* fcy */
			__asm	subpd	xmm1,xmm3	/* fprod4 -= fcy */
			__asm	mulpd	xmm1,xmm5	/* fprod4 *= two26f */

		/* xmm2,4,6 FREE, mulh(q,lo) digits0,1,2 in xmm0,1,3 */

			/* If h < l, then calculate h-l+q; otherwise calculate h-l. Use leading 52 bits to approximate the full 78-bit compare.
			   Result is in [0, q). */
			__asm	movaps	xmm2,[eax+0x070]	/* fhi, top 52 bits */
			__asm	movaps	xmm6,xmm2			/* cpy of fhi */
			__asm	movaps	xmm4,[edx+0x030]	/* fq, top 52 bits */
			__asm	mulpd	xmm3,xmm5			/* flo2 * two26f */
			__asm	mulpd	xmm5,[edx]			/* fq, low 26 bits */
			__asm	addpd	xmm3,xmm1			/* flo, top 52 bits */
			__asm	movaps	xmm1,[eax+0x060]	/* fhi, low 26 bits */
			__asm	cmppd	xmm6,xmm3,0x1		/* bitmask = (fhi < flo)? */
			__asm	subpd	xmm1,xmm0			/* (fhi - flo), low 26 bits */
			__asm	subpd	xmm2,xmm3			/* (fhi - flo), top 52 bits */
			__asm	movaps	xmm0,xmm4			/* cpy of qhi52 */

			__asm	andpd	xmm4,xmm6			/* qhi52 & bitmask */
			__asm	addpd	xmm2,xmm4			/* xhi = (h-l)hi + (qhi52 & bitmask) */
			__asm	andpd	xmm6,xmm5			/* bitmask & qlo26 */
			__asm	addpd	xmm1,xmm6			/* xlo = (h-l)lo + (qlo26 & bitmask) */
												/* qlo26 in xmm5, qhi52 in xmm0 */
			/* if((pshift >> j) & (uint64)1) { */
				__asm	mov	eax, pshift
				__asm	shr	eax, cl		/* j stored in ecx; only need lowest byte for shift count */
				__asm	and	eax, 0x1	/* (p odd)? = ((pshift >> j) & (uint64)1) */
		#if 1	// Branchless version of the conditional doubling
				__asm	xor	eax, 0x1	/* If result = 0, want a 0x0000.... double bitmask, otherwise want a 0xffff... mask */
				__asm	movd	xmm6,eax	/* xmm[0:31] <- eax */
				__asm	pshufd	xmm6,xmm6,0	/* Broadcast low 32 bits of xmm6 to all 4 slots of xmm0 (But only care about bottom 2 slots of result) */
				__asm cvtdq2pd	xmm6,xmm6	/* Convert to double */
				__asm	xorpd	xmm3,xmm3	/* All 0s */
				__asm	cmppd	xmm6,xmm3,0x0	/* bitmask = ((pshift >> j) odd)? */

				/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */
				__asm	movaps	xmm4,xmm2		/* cpy of xhi */
				__asm	movaps	xmm3,xmm1		/* cpy of xlo */
				__asm	andpd	xmm4,xmm6		/* xhi_ & bitmask */
				__asm	andpd	xmm3,xmm6		/* xlo_ & bitmask */
				__asm	addpd	xmm2,xmm4		/* top 52 bits */
				__asm	addpd	xmm1,xmm3		/* low 26 bits */

				/* If x > q, subtract q: */
				__asm	movaps	xmm6,xmm0		/* cpy of qhi */
				__asm	cmppd	xmm6,xmm2,0x2	/* bitmask = (qhi <= xhi)? */
				__asm	andpd	xmm0,xmm6		/* qhi52 & bitmask */
				__asm	andpd	xmm5,xmm6		/* qlo26 & bitmask */
				__asm	subpd	xmm2,xmm0		/* x % q, top 52 bits */
				__asm	subpd	xmm1,xmm5		/* x % q, low 26 bits */
		#else
			__asm	je	SHORT twopmodq78_3wdq2
				/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */
				__asm	movaps	xmm6,xmm0		/* cpy of qhi */
				__asm	addpd	xmm2,xmm2		/* top 52 bits */
				__asm	addpd	xmm1,xmm1		/* low 26 bits */

				/* If x > q, subtract q: */
				__asm	cmppd	xmm6,xmm2,0x2	/* bitmask = (qhi <= xhi)? */
				__asm	andpd	xmm0,xmm6		/* qhi52 & bitmask */
				__asm	andpd	xmm5,xmm6		/* qlo26 & bitmask */
				__asm	subpd	xmm2,xmm0		/* x % q, top 52 bits */
				__asm	subpd	xmm1,xmm5		/* x % q, low 26 bits */
			twopmodq78_3wdq2:
		#endif
			/* } */

			/* Normalize the result:
			__fcy  = NINT(__xlo*TWO26FLINV);
			__xlo -=  __fcy*TWO26FLOAT;
			*/
			__asm	movaps	xmm4,[ebx-0x010]	/* two26f */
			__asm	movaps	xmm3,[ebx      ]	/* two26i */

			__asm	mulpd	xmm1,xmm3		/* xlo *= two26i */
			__asm	mulpd	xmm2,xmm3		/* xhi *= two26i */
			__asm	movaps	xmm0,xmm1	/* Init: fcy = cpy of ~xlo */
			__asm	addpd	xmm1,xmm7
			__asm	subpd	xmm1,xmm7	/* fcy */
			__asm	subpd	xmm0,xmm1	/* fx0 -= fcy */
			__asm	mulpd	xmm0,xmm4	/* fx0 *= two26f */
			__asm	mulpd	xmm1,xmm3	/* fcy *= two26i */

			__asm	addpd	xmm2,xmm1	/* Add carry into xhi */
			__asm	movaps	xmm1,xmm2	/* Init: fcy = cpy of ~xhi */
			__asm	addpd	xmm2,xmm7
			__asm	subpd	xmm2,xmm7	/* fx2 = fcy (no divide by 2^26) */
			__asm	subpd	xmm1,xmm2	/* fx1 -= fcy */
			__asm	mulpd	xmm1,xmm4	/* fx1 *= two26f */

			/* Move high 2 words of result into input registers expected by start of loop body: */
			__asm	movaps	xmm4,xmm2	/* fx2 */
			__asm	movaps	xmm2,xmm1	/* fx1 */

		#ifdef DEBUG_SSE2
		/******************************************************/
		// Dump SSE2 registers to memory and compare to pure-int result: */
		__asm	mov	eax, fx0
		__asm	movaps	[eax      ],xmm0	/* flo0 */
		__asm	movaps	[eax+0x010],xmm2	/* flo1 */
		__asm	movaps	[eax+0x020],xmm4	/* flo2 */
		/******************************************************/
		#endif

		__asm	sub ecx, 1	/* j-- */
		__asm	cmp ecx, 0	/* j > 0 ?	*/
		__asm	jge LoopBeg	/* if (j >= 0), Loop */
		LoopEnd:

	  #else	/* GCC-style inline ASM: */

		for(j = start_index-2; j >= 0; j--)
		{
		#if OS_BITS == 32

		__asm__ volatile (\
			"/* SQR_LOHI78_3WORD_DOUBLE_q2(fx, flo,fhi): */\n\t"\
			"movl	%[__fx0],%%eax\n\t"\
			"movl	%[__two26i],%%esi\n\t"\
			"movaps	-0x20(%%esi),%%xmm6\n\t"\
			"/* fx0,1,2 assumed in xmm0,2,4 on loop entry */\n\t"\
	"/* GCC builds require us to reload the xmm0/2/4 registers on each pass thru the loop: */\n\t"\
		"movaps	     (%%eax),%%xmm0	\n\t"\
		"movaps	0x010(%%eax),%%xmm2	\n\t"\
		"movaps	0x020(%%eax),%%xmm4	\n\t"\
			"mulpd	%%xmm6,%%xmm0\n\t"\
			"mulpd	%%xmm6,%%xmm2\n\t"\
			"mulpd	%%xmm6,%%xmm4\n\t"\
			"movaps	0x10(%%esi),%%xmm7\n\t"\
			"movaps	%%xmm0,%%xmm1\n\t"\
			"movaps	%%xmm2,%%xmm3\n\t"\
			"addpd	%%xmm1,%%xmm1\n\t"\
			"addpd	%%xmm3,%%xmm3\n\t"\
			"movaps	%%xmm1,%%xmm5\n\t"\
			"mulpd	%%xmm0,%%xmm0\n\t"\
			"mulpd	%%xmm2,%%xmm1\n\t"\
			"mulpd	%%xmm2,%%xmm2\n\t"\
			"mulpd	%%xmm4,%%xmm5\n\t"\
			"mulpd	%%xmm4,%%xmm3\n\t"\
			"mulpd	%%xmm4,%%xmm4\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm0\n\t"\
			"mulpd	-0x10(%%esi),%%xmm0\n\t"\
			"mulpd	     (%%esi),%%xmm6\n\t"\
			"/* Digit 1: */\n\t"\
			"addpd	%%xmm6,%%xmm1\n\t"\
			"movaps	%%xmm1,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm1\n\t"\
			"mulpd	-0x10(%%esi),%%xmm1\n\t"\
			"mulpd	     (%%esi),%%xmm6\n\t"\
			"/* Digit 2: */\n\t"\
			"addpd	%%xmm5,%%xmm2\n\t"\
			"movaps	-0x10(%%esi),%%xmm5\n\t"\
			"addpd	%%xmm6,%%xmm2\n\t"\
			"movaps	%%xmm2,%%xmm6\n\t"\
			"subpd	0x20(%%esi),%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm2\n\t"\
			"mulpd	%%xmm5,%%xmm2\n\t"\
			"mulpd	    (%%esi),%%xmm6\n\t"\
			"/* Digit 3: */\n\t"\
			"addpd	%%xmm6,%%xmm3\n\t"\
			"movaps	%%xmm3,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm3\n\t"\
			"mulpd	%%xmm5,%%xmm3\n\t"\
			"mulpd	    (%%esi),%%xmm6\n\t"\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm6,%%xmm4\n\t"\
			"movaps	%%xmm4,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm4\n\t"\
			"mulpd	%%xmm5,%%xmm4\n\t"\
			"/* Digit 5 = the carry. flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6 */\n\t"\
			"mulpd	%%xmm5,%%xmm6\n\t"\
			"addpd	%%xmm6,%%xmm4\n\t"\
			"movaps	%%xmm3,0x60(%%eax)\n\t"\
			"movaps	%%xmm4,0x70(%%eax)\n\t"\
			"/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */\n\t"\
			"movl	%[__fqinv0],%%edx\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm3\n\t"\
			"movaps	%%xmm0,%%xmm4\n\t"\
			"mulpd	    (%%edx),%%xmm0\n\t"\
			"movaps	%%xmm0,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm0\n\t"\
			"mulpd	%%xmm5,%%xmm0\n\t"\
			"mulpd	    (%%esi),%%xmm6\n\t"\
			"/* Digit 1: */\n\t"\
			"mulpd	0x10(%%edx),%%xmm3\n\t"\
			"addpd	%%xmm6,%%xmm3\n\t"\
			"movaps	%%xmm1,%%xmm6\n\t"\
			"mulpd	    (%%edx),%%xmm1\n\t"\
			"addpd	%%xmm3,%%xmm1\n\t"\
			"movaps	%%xmm1,%%xmm3\n\t"\
			"addpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm3,%%xmm1\n\t"\
			"mulpd	%%xmm5,%%xmm1\n\t"\
			"mulpd	    (%%esi),%%xmm3\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"mulpd	    (%%edx),%%xmm2\n\t"\
			"mulpd	0x10(%%edx),%%xmm6\n\t"\
			"mulpd	0x20(%%edx),%%xmm4\n\t"\
			"addpd	%%xmm6,%%xmm2\n\t"\
			"addpd	%%xmm4,%%xmm3\n\t"\
			"addpd	%%xmm3,%%xmm2\n\t"\
			"movaps	%%xmm2,%%xmm3\n\t"\
			"subpd	0x20(%%esi),%%xmm3\n\t"\
			"addpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm3,%%xmm2\n\t"\
			"mulpd	%%xmm5,%%xmm2\n\t"\
			"/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */\n\t"\
			"movl	%[__fq0],%%edx\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm3\n\t"\
			"movaps	%%xmm0,%%xmm4\n\t"\
			"mulpd	    (%%edx),%%xmm0\n\t"\
			"addpd	%%xmm7,%%xmm0\n\t"\
			"subpd	%%xmm7,%%xmm0\n\t"\
			"mulpd	    (%%esi),%%xmm0\n\t"\
			"/* Digit 1: */\n\t"\
			"movaps	%%xmm1,%%xmm6\n\t"\
			"mulpd	0x10(%%edx),%%xmm3\n\t"\
			"addpd	%%xmm0,%%xmm3\n\t"\
			"mulpd	    (%%edx),%%xmm1\n\t"\
			"addpd	%%xmm3,%%xmm1\n\t"\
			"addpd	%%xmm7,%%xmm1\n\t"\
			"subpd	%%xmm7,%%xmm1\n\t"\
			"mulpd	    (%%esi),%%xmm1\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"movaps	%%xmm2,%%xmm3\n\t"\
			"movaps	%%xmm6,%%xmm0\n\t"\
			"mulpd	    (%%edx),%%xmm2\n\t"\
			"mulpd	0x10(%%edx),%%xmm6\n\t"\
			"mulpd	0x20(%%edx),%%xmm4\n\t"\
			"addpd	%%xmm6,%%xmm2\n\t"\
			"addpd	%%xmm4,%%xmm1\n\t"\
			"addpd	%%xmm1,%%xmm2\n\t"\
			"subpd	0x20(%%esi),%%xmm2\n\t"\
			"addpd	%%xmm7,%%xmm2\n\t"\
			"subpd	%%xmm7,%%xmm2\n\t"\
			"mulpd	    (%%esi),%%xmm2\n\t"\
			"/* Precompute all the needed partial products: */\n\t"\
			"movaps	%%xmm3,%%xmm1\n\t"\
			"mulpd	0x20(%%edx),%%xmm0\n\t"\
			"mulpd	0x10(%%edx),%%xmm1\n\t"\
			"mulpd	0x20(%%edx),%%xmm3\n\t"\
			"/* Digit 3: */\n\t"\
			"addpd	%%xmm2,%%xmm1\n\t"\
			"addpd	%%xmm1,%%xmm0\n\t"\
			"movaps	%%xmm0,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm0\n\t"\
			"mulpd	%%xmm5,%%xmm0\n\t"\
			"mulpd	    (%%esi),%%xmm6\n\t"\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm6,%%xmm3\n\t"\
			"movaps	%%xmm3,%%xmm1\n\t"\
			"addpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm3,%%xmm1\n\t"\
			"mulpd	%%xmm5,%%xmm1\n\t"\
			"/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\n\t"\
			"movaps	0x70(%%eax),%%xmm2\n\t"\
			"movaps	%%xmm2,%%xmm6\n\t"\
			"movaps	0x30(%%edx),%%xmm4\n\t"\
			"mulpd	%%xmm5,%%xmm3\n\t"\
			"mulpd	    (%%edx),%%xmm5\n\t"\
			"addpd	%%xmm1,%%xmm3\n\t"\
			"movaps	0x60(%%eax),%%xmm1\n\t"\
			"cmppd	$0x1,%%xmm3,%%xmm6\n\t"\
			"subpd	%%xmm0,%%xmm1\n\t"\
			"subpd	%%xmm3,%%xmm2\n\t"\
			"movaps	%%xmm4,%%xmm0\n\t"\
			"andpd	%%xmm6,%%xmm4\n\t"\
			"addpd	%%xmm4,%%xmm2\n\t"\
			"andpd	%%xmm5,%%xmm6\n\t"\
			"addpd	%%xmm6,%%xmm1\n\t"\
			"/* qlo26 in xmm5, qhi52 in xmm0 */\n\t"\
			"/* if((pshift >> j) & (uint64)1) { */\n\t"\
			"movl	%[__pshift],%%eax\n\t"\
			"movl	%[__j],%%ecx\n\t"\
			"shrl	%%cl,%%eax\n\t"\
			"andl	$0x1,%%eax\n\t"\
			"/* Branchless version of the conditional doubling: */\n\t"\
			"xorl	$0x1,%%eax\n\t"\
			"movd	%%eax,%%xmm6\n\t"\
			"pshufd	$0,%%xmm6,%%xmm6\n\t"\
			"cvtdq2pd	%%xmm6,%%xmm6\n\t"\
			"xorpd	%%xmm3,%%xmm3\n\t"\
			"cmppd	$0x0,%%xmm3,%%xmm6\n\t"\
			"/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */\n\t"\
			"movaps	%%xmm2,%%xmm4\n\t"\
			"movaps	%%xmm1,%%xmm3\n\t"\
			"andpd	%%xmm6,%%xmm4\n\t"\
			"andpd	%%xmm6,%%xmm3\n\t"\
			"addpd	%%xmm4,%%xmm2\n\t"\
			"addpd	%%xmm3,%%xmm1\n\t"\
			"/* If x > q, subtract q: */\n\t"\
			"movaps	%%xmm0,%%xmm6\n\t"\
			"cmppd	$0x2,%%xmm2,%%xmm6\n\t"\
			"andpd	%%xmm6,%%xmm0\n\t"\
			"andpd	%%xmm6,%%xmm5\n\t"\
			"subpd	%%xmm0,%%xmm2\n\t"\
			"subpd	%%xmm5,%%xmm1\n\t"\
			"/* } */\n\t"\
			"/* Normalize the result: */\n\t"\
			"movaps	-0x10(%%esi),%%xmm4\n\t"\
			"movaps	     (%%esi),%%xmm3\n\t"\
			"mulpd	%%xmm3,%%xmm1\n\t"\
			"mulpd	%%xmm3,%%xmm2\n\t"\
			"movaps	%%xmm1,%%xmm0\n\t"\
			"addpd	%%xmm7,%%xmm1\n\t"\
			"subpd	%%xmm7,%%xmm1\n\t"\
			"subpd	%%xmm1,%%xmm0\n\t"\
			"mulpd	%%xmm4,%%xmm0\n\t"\
			"mulpd	%%xmm3,%%xmm1\n\t"\
			"addpd	%%xmm1,%%xmm2\n\t"\
			"movaps	%%xmm2,%%xmm1\n\t"\
			"addpd	%%xmm7,%%xmm2\n\t"\
			"subpd	%%xmm7,%%xmm2\n\t"\
			"subpd	%%xmm2,%%xmm1\n\t"\
			"mulpd	%%xmm4,%%xmm1\n\t"\
			"/* Move high 2 words of result into input registers expected by start of loop body: */\n\t"\
			"movaps	%%xmm2,%%xmm4	/* fx2 */\n\t"\
			"movaps	%%xmm1,%%xmm2	/* fx1 */\n\t"\
	"/* GCC builds require us to save the xmm0/2/4 registers on each pass thru the loop: */\n\t"\
		"movl	%[__fx0],%%eax 	\n\t"\
		"movaps	%%xmm0,     (%%eax)	\n\t"\
		"movaps	%%xmm2,0x010(%%eax)	\n\t"\
		"movaps	%%xmm4,0x020(%%eax)	\n\t"\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			 ,[__fx0]	 "m" (fx0)		\
			 ,[__fq0]	 "m" (fq0)		\
			 ,[__pshift] "m" (pshift)	\
			 ,[__j]		 "m" (j)		\
			: "cl","eax","esi","ecx","edx"	/* Clobbered registers */\
		);

		#else

		__asm__ volatile (\
			"/* SQR_LOHI78_3WORD_DOUBLE_q2(fx, flo,fhi): */\n\t"\
			"movq	%[__fx0],%%rax\n\t"\
			"movq	%[__two26i],%%rsi\n\t"\
			"movaps	-0x20(%%rsi),%%xmm6\n\t"\
			"/* fx0,1,2 assumed in xmm0,2,4 on loop entry */\n\t"\
	"/* GCC builds require us to reload the xmm0/2/4 registers on each pass thru the loop: */\n\t"\
		"movaps	     (%%rax),%%xmm0	\n\t"\
		"movaps	0x010(%%rax),%%xmm2	\n\t"\
		"movaps	0x020(%%rax),%%xmm4	\n\t"\
			"mulpd	%%xmm6,%%xmm0\n\t"\
			"mulpd	%%xmm6,%%xmm2\n\t"\
			"mulpd	%%xmm6,%%xmm4\n\t"\
			"movaps	0x10(%%rsi),%%xmm7\n\t"\
			"movaps	%%xmm0,%%xmm1\n\t"\
			"movaps	%%xmm2,%%xmm3\n\t"\
			"addpd	%%xmm1,%%xmm1\n\t"\
			"addpd	%%xmm3,%%xmm3\n\t"\
			"movaps	%%xmm1,%%xmm5\n\t"\
			"mulpd	%%xmm0,%%xmm0\n\t"\
			"mulpd	%%xmm2,%%xmm1\n\t"\
			"mulpd	%%xmm2,%%xmm2\n\t"\
			"mulpd	%%xmm4,%%xmm5\n\t"\
			"mulpd	%%xmm4,%%xmm3\n\t"\
			"mulpd	%%xmm4,%%xmm4\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm0\n\t"\
			"mulpd	-0x10(%%rsi),%%xmm0\n\t"\
			"mulpd	     (%%rsi),%%xmm6\n\t"\
			"/* Digit 1: */\n\t"\
			"addpd	%%xmm6,%%xmm1\n\t"\
			"movaps	%%xmm1,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm1\n\t"\
			"mulpd	-0x10(%%rsi),%%xmm1\n\t"\
			"mulpd	     (%%rsi),%%xmm6\n\t"\
			"/* Digit 2: */\n\t"\
			"addpd	%%xmm5,%%xmm2\n\t"\
			"movaps	-0x10(%%rsi),%%xmm5\n\t"\
			"addpd	%%xmm6,%%xmm2\n\t"\
			"movaps	%%xmm2,%%xmm6\n\t"\
			"subpd	0x20(%%rsi),%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm2\n\t"\
			"mulpd	%%xmm5,%%xmm2\n\t"\
			"mulpd	    (%%rsi),%%xmm6\n\t"\
			"/* Digit 3: */\n\t"\
			"addpd	%%xmm6,%%xmm3\n\t"\
			"movaps	%%xmm3,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm3\n\t"\
			"mulpd	%%xmm5,%%xmm3\n\t"\
			"mulpd	    (%%rsi),%%xmm6\n\t"\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm6,%%xmm4\n\t"\
			"movaps	%%xmm4,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm4\n\t"\
			"mulpd	%%xmm5,%%xmm4\n\t"\
			"/* Digit 5 = the carry. flo0,1,2 in xmm0,1,2; fhi0,1,2 in xmm3,4,6 */\n\t"\
			"mulpd	%%xmm5,%%xmm6\n\t"\
			"addpd	%%xmm6,%%xmm4\n\t"\
			"movaps	%%xmm3,0x60(%%rax)\n\t"\
			"movaps	%%xmm4,0x70(%%rax)\n\t"\
			"/* flo = MULL78_3WORD_DOUBLE_q2(flo, fqinv): */\n\t"\
			"movq	%[__fqinv0],%%rdx\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm3\n\t"\
			"movaps	%%xmm0,%%xmm4\n\t"\
			"mulpd	    (%%rdx),%%xmm0\n\t"\
			"movaps	%%xmm0,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm0\n\t"\
			"mulpd	%%xmm5,%%xmm0\n\t"\
			"mulpd	    (%%rsi),%%xmm6\n\t"\
			"/* Digit 1: */\n\t"\
			"mulpd	0x10(%%rdx),%%xmm3\n\t"\
			"addpd	%%xmm6,%%xmm3\n\t"\
			"movaps	%%xmm1,%%xmm6\n\t"\
			"mulpd	    (%%rdx),%%xmm1\n\t"\
			"addpd	%%xmm3,%%xmm1\n\t"\
			"movaps	%%xmm1,%%xmm3\n\t"\
			"addpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm3,%%xmm1\n\t"\
			"mulpd	%%xmm5,%%xmm1\n\t"\
			"mulpd	    (%%rsi),%%xmm3\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"mulpd	    (%%rdx),%%xmm2\n\t"\
			"mulpd	0x10(%%rdx),%%xmm6\n\t"\
			"mulpd	0x20(%%rdx),%%xmm4\n\t"\
			"addpd	%%xmm6,%%xmm2\n\t"\
			"addpd	%%xmm4,%%xmm3\n\t"\
			"addpd	%%xmm3,%%xmm2\n\t"\
			"movaps	%%xmm2,%%xmm3\n\t"\
			"subpd	0x20(%%rsi),%%xmm3\n\t"\
			"addpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm3,%%xmm2\n\t"\
			"mulpd	%%xmm5,%%xmm2\n\t"\
			"/* MULH96(q,lo,lo) --> lo = (q*lo)/2^78 */\n\t"\
			"movq	%[__fq0],%%rdx\n\t"\
			"/* Digit 0: */\n\t"\
			"movaps	%%xmm0,%%xmm3\n\t"\
			"movaps	%%xmm0,%%xmm4\n\t"\
			"mulpd	    (%%rdx),%%xmm0\n\t"\
			"addpd	%%xmm7,%%xmm0\n\t"\
			"subpd	%%xmm7,%%xmm0\n\t"\
			"mulpd	    (%%rsi),%%xmm0\n\t"\
			"/* Digit 1: */\n\t"\
			"movaps	%%xmm1,%%xmm6\n\t"\
			"mulpd	0x10(%%rdx),%%xmm3\n\t"\
			"addpd	%%xmm0,%%xmm3\n\t"\
			"mulpd	    (%%rdx),%%xmm1\n\t"\
			"addpd	%%xmm3,%%xmm1\n\t"\
			"addpd	%%xmm7,%%xmm1\n\t"\
			"subpd	%%xmm7,%%xmm1\n\t"\
			"mulpd	    (%%rsi),%%xmm1\n\t"\
			"/* Digit 2: Require output to be nonnegative, so leave MSW unbalanced: */\n\t"\
			"movaps	%%xmm2,%%xmm3\n\t"\
			"movaps	%%xmm6,%%xmm0\n\t"\
			"mulpd	    (%%rdx),%%xmm2\n\t"\
			"mulpd	0x10(%%rdx),%%xmm6\n\t"\
			"mulpd	0x20(%%rdx),%%xmm4\n\t"\
			"addpd	%%xmm6,%%xmm2\n\t"\
			"addpd	%%xmm4,%%xmm1\n\t"\
			"addpd	%%xmm1,%%xmm2\n\t"\
			"subpd	0x20(%%rsi),%%xmm2\n\t"\
			"addpd	%%xmm7,%%xmm2\n\t"\
			"subpd	%%xmm7,%%xmm2\n\t"\
			"mulpd	    (%%rsi),%%xmm2\n\t"\
			"/* Precompute all the needed partial products: */\n\t"\
			"movaps	%%xmm3,%%xmm1\n\t"\
			"mulpd	0x20(%%rdx),%%xmm0\n\t"\
			"mulpd	0x10(%%rdx),%%xmm1\n\t"\
			"mulpd	0x20(%%rdx),%%xmm3\n\t"\
			"/* Digit 3: */\n\t"\
			"addpd	%%xmm2,%%xmm1\n\t"\
			"addpd	%%xmm1,%%xmm0\n\t"\
			"movaps	%%xmm0,%%xmm6\n\t"\
			"addpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm7,%%xmm6\n\t"\
			"subpd	%%xmm6,%%xmm0\n\t"\
			"mulpd	%%xmm5,%%xmm0\n\t"\
			"mulpd	    (%%rsi),%%xmm6\n\t"\
			"/* Digit 4: */\n\t"\
			"addpd	%%xmm6,%%xmm3\n\t"\
			"movaps	%%xmm3,%%xmm1\n\t"\
			"addpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm7,%%xmm3\n\t"\
			"subpd	%%xmm3,%%xmm1\n\t"\
			"mulpd	%%xmm5,%%xmm1\n\t"\
			"/* If h < l, calculate h-l+q; otherwise h-l. Use leading 52 bits to approximate the full 78-bit compare. Result is in [0, q). */\n\t"\
			"movaps	0x70(%%rax),%%xmm2\n\t"\
			"movaps	%%xmm2,%%xmm6\n\t"\
			"movaps	0x30(%%rdx),%%xmm4\n\t"\
			"mulpd	%%xmm5,%%xmm3\n\t"\
			"mulpd	    (%%rdx),%%xmm5\n\t"\
			"addpd	%%xmm1,%%xmm3\n\t"\
			"movaps	0x60(%%rax),%%xmm1\n\t"\
			"cmppd	$0x1,%%xmm3,%%xmm6\n\t"\
			"subpd	%%xmm0,%%xmm1\n\t"\
			"subpd	%%xmm3,%%xmm2\n\t"\
			"movaps	%%xmm4,%%xmm0\n\t"\
			"andpd	%%xmm6,%%xmm4\n\t"\
			"addpd	%%xmm4,%%xmm2\n\t"\
			"andpd	%%xmm5,%%xmm6\n\t"\
			"addpd	%%xmm6,%%xmm1\n\t"\
			"/* qlo26 in xmm5, qhi52 in xmm0 */\n\t"\
			"/* if((pshift >> j) & (uint64)1) { */\n\t"\
			"movq	%[__pshift],%%rax\n\t"\
			"movq	%[__j],%%rcx\n\t"\
			"shrq	%%cl,%%rax\n\t"\
			"andq	$0x1,%%rax\n\t"\
			"/* Branchless version of the conditional doubling: */\n\t"\
			"xorq	$0x1,%%rax\n\t"\
			"movd	%%eax,%%xmm6\n\t"\
			"pshufd	$0,%%xmm6,%%xmm6\n\t"\
			"cvtdq2pd	%%xmm6,%%xmm6\n\t"\
			"xorpd	%%xmm3,%%xmm3\n\t"\
			"cmppd	$0x0,%%xmm3,%%xmm6\n\t"\
			"/* Double: Use that high part < 2^52 (strictly >= 0): xhi52,xlo26 in xmm2,xmm1; qhi52,qlo26 in xmm0,xmm5 */\n\t"\
			"movaps	%%xmm2,%%xmm4\n\t"\
			"movaps	%%xmm1,%%xmm3\n\t"\
			"andpd	%%xmm6,%%xmm4\n\t"\
			"andpd	%%xmm6,%%xmm3\n\t"\
			"addpd	%%xmm4,%%xmm2\n\t"\
			"addpd	%%xmm3,%%xmm1\n\t"\
			"/* If x > q, subtract q: */\n\t"\
			"movaps	%%xmm0,%%xmm6\n\t"\
			"cmppd	$0x2,%%xmm2,%%xmm6\n\t"\
			"andpd	%%xmm6,%%xmm0\n\t"\
			"andpd	%%xmm6,%%xmm5\n\t"\
			"subpd	%%xmm0,%%xmm2\n\t"\
			"subpd	%%xmm5,%%xmm1\n\t"\
			"/* } */\n\t"\
			"/* Normalize the result: */\n\t"\
			"movaps	-0x10(%%rsi),%%xmm4\n\t"\
			"movaps	     (%%rsi),%%xmm3\n\t"\
			"mulpd	%%xmm3,%%xmm1\n\t"\
			"mulpd	%%xmm3,%%xmm2\n\t"\
			"movaps	%%xmm1,%%xmm0\n\t"\
			"addpd	%%xmm7,%%xmm1\n\t"\
			"subpd	%%xmm7,%%xmm1\n\t"\
			"subpd	%%xmm1,%%xmm0\n\t"\
			"mulpd	%%xmm4,%%xmm0\n\t"\
			"mulpd	%%xmm3,%%xmm1\n\t"\
			"addpd	%%xmm1,%%xmm2\n\t"\
			"movaps	%%xmm2,%%xmm1\n\t"\
			"addpd	%%xmm7,%%xmm2\n\t"\
			"subpd	%%xmm7,%%xmm2\n\t"\
			"subpd	%%xmm2,%%xmm1\n\t"\
			"mulpd	%%xmm4,%%xmm1\n\t"\
			"/* Move high 2 words of result into input registers expected by start of loop body: */\n\t"\
			"movaps	%%xmm2,%%xmm4	/* fx2 */\n\t"\
			"movaps	%%xmm1,%%xmm2	/* fx1 */\n\t"\
	"/* GCC builds require us to save the xmm0/2/4 registers on each pass thru the loop: */\n\t"\
		"movq	%[__fx0],%%rax 	\n\t"\
		"movaps	%%xmm0,     (%%rax)	\n\t"\
		"movaps	%%xmm2,0x010(%%rax)	\n\t"\
		"movaps	%%xmm4,0x020(%%rax)	\n\t"\
			:					/* outputs: none */\
			: [__fqinv0] "m" (fqinv0)	/* All inputs from memory addresses here */\
			 ,[__two26i] "m" (two26i)	\
			 ,[__fx0]	 "m" (fx0)		\
			 ,[__fq0]	 "m" (fq0)		\
			 ,[__pshift] "m" (pshift)	\
			 ,[__j]		 "m" (j)		\
			: "cl","eax","rax","rsi","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
		);

		#endif	/* OS_BITS */

		}	/* for(j...) */

	  #endif	/* COMPILER_TYPE_ */

	  #ifdef DEBUG_SSE2
		fprintf(stderr, "[f|g]x0  = %20.0f, %20.0f\n",*fx0,*gx0);
		fprintf(stderr, "[f|g]x1  = %20.0f, %20.0f\n",*fx1,*gx1);
		fprintf(stderr, "[f|g]x2  = %20.0f, %20.0f\n",*fx2,*gx2);
		fprintf(stderr, "\n");
		/******************************************************/
	  #endif

	#endif /* if(0) */

#endif	/* USE_SSE2 */

#ifdef DEBUG_SSE2
//	exit(0);
#endif

#ifdef USE_SSE2

	// Dump SSE2 registers to memory: */
	#if defined(COMPILER_TYPE_MSVC)

		__asm	mov	eax, fx0
		__asm	movaps	[eax      ],xmm0	/* flo0 */
		__asm	movaps	[eax+0x010],xmm2	/* flo1 */
		__asm	movaps	[eax+0x020],xmm4	/* flo2 */

	#endif

#endif


	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
#ifdef USE_SSE2
	CVT78_3WORD_DOUBLE_UINT96(*fx0,*fx1,*fx2, x0);
	CVT78_3WORD_DOUBLE_UINT96(*gx0,*gx1,*gx2, x1);
#else
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
	CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
#endif


	ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
	ADD96(x1,x1,x1);

	SUB96(x0,q0,x0);
	SUB96(x1,q1,x1);

	r = 0;
	if(CMPEQ96(x0, ONE96)) r +=  1;
	if(CMPEQ96(x1, ONE96)) r +=  2;
	return(r);
}

/*** 4-trial-factor version ***/
uint64 twopmodq78_3WORD_DOUBLE_q4(uint64 p, uint96 q0, uint96 q1, uint96 q2, uint96 q3)
{
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint64 lead7, tmp0, tmp1, tmp2, tmp3, r;
	uint96 qinv0, qhalf0, x0, lo0
		 , qinv1, qhalf1, x1, lo1
		 , qinv2, qhalf2, x2, lo2
		 , qinv3, qhalf3, x3, lo3;
	uint192 prod192;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2, fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2;
	double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2;
	double hq0,hq1,hq2, hqinv0,hqinv1,hqinv2, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hhi1,hhi2;
	double iq0,iq1,iq2, iqinv0,iqinv1,iqinv2, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ihi1,ihi2;

	ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q4 : (q0.d1 >> 14) != 0");
	ASSERT(HERE, (q1.d1 >> 14) == 0, "twopmodq78_q4 : (q1.d1 >> 14) != 0");
	ASSERT(HERE, (q2.d1 >> 14) == 0, "twopmodq78_q4 : (q2.d1 >> 14) != 0");
	ASSERT(HERE, (q3.d1 >> 14) == 0, "twopmodq78_q4 : (q3.d1 >> 14) != 0");

	/* Convert q to floating form: */
	CVT_UINT78_3WORD_DOUBLE(q0, fq0,fq1,fq2);
	CVT_UINT78_3WORD_DOUBLE(q1, gq0,gq1,gq2);
	CVT_UINT78_3WORD_DOUBLE(q2, hq0,hq1,hq2);
	CVT_UINT78_3WORD_DOUBLE(q3, iq0,iq1,iq2);

	RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST96(q1, 1, qhalf1);
	RSHIFT_FAST96(q2, 1, qhalf2);
	RSHIFT_FAST96(q3, 1, qhalf3);

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 78;
		j = leadz64(pshift);
		lead7 = ((pshift<<j) >> 57);
		if(lead7 > 77)
		{
			lead7 >>= 1;
			start_index =  64-j-6;	/* Use only the leftmost 6 bits */
		}
		else
			start_index =  64-j-7;

		zshift = 77 - lead7;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}

	qinv0.d0 = (q0.d0 + q0.d0 + q0.d0) ^ (uint64)2;	qinv0.d1 = (uint64)0;
	qinv1.d0 = (q1.d0 + q1.d0 + q1.d0) ^ (uint64)2;	qinv1.d1 = (uint64)0;
	qinv2.d0 = (q2.d0 + q2.d0 + q2.d0) ^ (uint64)2;	qinv2.d1 = (uint64)0;
	qinv3.d0 = (q3.d0 + q3.d0 + q3.d0) ^ (uint64)2;	qinv3.d1 = (uint64)0;

	for(j = 0; j < 4; j++)
	{
		tmp0 = q0.d0*qinv0.d0;
		tmp1 = q1.d0*qinv1.d0;
		tmp2 = q2.d0*qinv2.d0;
		tmp3 = q3.d0*qinv3.d0;

		qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
		qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
		qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);
		qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);
	}

#ifdef MUL_LOHI64_SUBROUTINE
	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + __MULH64(q0.d0, qinv0.d0));
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + __MULH64(q1.d0, qinv1.d0));
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + __MULH64(q2.d0, qinv2.d0));
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + __MULH64(q3.d0, qinv3.d0));
#else
	MULH64(q0.d0, qinv0.d0, tmp0);
	MULH64(q1.d0, qinv1.d0, tmp1);
	MULH64(q2.d0, qinv2.d0, tmp2);
	MULH64(q3.d0, qinv3.d0, tmp3);

	qinv0.d1 = -qinv0.d0*(q0.d1*qinv0.d0 + tmp0);
	qinv1.d1 = -qinv1.d0*(q1.d1*qinv1.d0 + tmp1);
	qinv2.d1 = -qinv2.d0*(q2.d1*qinv2.d0 + tmp2);
	qinv3.d1 = -qinv3.d0*(q3.d1*qinv3.d0 + tmp3);
#endif
	qinv0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
	qinv1.d1 &= 0x0000000000003fff;
	qinv2.d1 &= 0x0000000000003fff;
	qinv3.d1 &= 0x0000000000003fff;

	/* Convert qinv to floating form: */
	CVT_UINT78_3WORD_DOUBLE(qinv0, fqinv0,fqinv1,fqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv1, gqinv0,gqinv1,gqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv2, hqinv0,hqinv1,hqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv3, iqinv0,iqinv1,iqinv2);

	/* Since zstart is a power of two < 2^96, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL96(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT96(qinv0, zshift, lo0);	lo0.d1 &= 0x0000000000003fff;	/* Only want the lower 14 bits here */
	LSHIFT96(qinv1, zshift, lo1);	lo1.d1 &= 0x0000000000003fff;
	LSHIFT96(qinv2, zshift, lo2);	lo2.d1 &= 0x0000000000003fff;
	LSHIFT96(qinv3, zshift, lo3);	lo3.d1 &= 0x0000000000003fff;

	MUL_LOHI96_PROD192(q0,lo0,prod192);	RSHIFT192(prod192,78,prod192);	lo0.d0 = prod192.d0;	lo0.d1 = prod192.d1;
	MUL_LOHI96_PROD192(q1,lo1,prod192);	RSHIFT192(prod192,78,prod192);	lo1.d0 = prod192.d0;	lo1.d1 = prod192.d1;
	MUL_LOHI96_PROD192(q2,lo2,prod192);	RSHIFT192(prod192,78,prod192);	lo2.d0 = prod192.d0;	lo2.d1 = prod192.d1;
	MUL_LOHI96_PROD192(q3,lo3,prod192);	RSHIFT192(prod192,78,prod192);	lo3.d0 = prod192.d0;	lo3.d1 = prod192.d1;

	/* hi = 0 in this instance, which simplifies things. */
	SUB96(q0, lo0, x0);	/* Put the result in lo (rather than x), to ease overflow check below */
	SUB96(q1, lo1, x1);
	SUB96(q2, lo2, x2);
	SUB96(q3, lo3, x3);

	if((pshift >> j) & (uint64)1)
	{
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(CMPUGT96(x0, qhalf0)){ ADD96(x0, x0, x0); SUB96(x0, q0, x0); }else{ ADD96(x0, x0, x0); }
		if(CMPUGT96(x1, qhalf1)){ ADD96(x1, x1, x1); SUB96(x1, q1, x1); }else{ ADD96(x1, x1, x1); }
		if(CMPUGT96(x2, qhalf2)){ ADD96(x2, x2, x2); SUB96(x2, q2, x2); }else{ ADD96(x2, x2, x2); }
		if(CMPUGT96(x3, qhalf3)){ ADD96(x3, x3, x3); SUB96(x3, q3, x3); }else{ ADD96(x3, x3, x3); }
	}

	/* Convert x to floating form: */
	CVT_UINT78_3WORD_DOUBLE(x0, fx0,fx1,fx2);
	CVT_UINT78_3WORD_DOUBLE(x1, gx0,gx1,gx2);
	CVT_UINT78_3WORD_DOUBLE(x2, hx0,hx1,hx2);
	CVT_UINT78_3WORD_DOUBLE(x3, ix0,ix1,ix2);

	for(j = start_index-2; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		/* SQR_LOHI96(x,lo,hi); */
		SQR_LOHI78_3WORD_DOUBLE_q4(
			  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2
			, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2
			, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hhi1,hhi2
			, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ihi1,ihi2
		);
		/* MULL96(lo,qinv,lo); */
		MULL78_3WORD_DOUBLE_q4(
			  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
			, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
			, hlo0,hlo1,hlo2, hqinv0,hqinv1,hqinv2, hlo0,hlo1,hlo2
			, ilo0,ilo1,ilo2, iqinv0,iqinv1,iqinv2, ilo0,ilo1,ilo2
		);
		/* MULH96(q,lo,lo); */
		MULH78_3WORD_DOUBLE_q4(
			  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2
			, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2
			, hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlo1,hlo2
			, iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilo1,ilo2
		);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
		{
			SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
			ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
		}
		NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

		if(CMPLT78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2))
		{
			SUB78_3WORD_DOUBLE(gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2);
			ADD78_3WORD_DOUBLE(glo0,glo1,glo2, ghi0,ghi1,ghi2, gx0,gx1,gx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2, gx0,gx1,gx2);
		}
		NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

		if(CMPLT78_3WORD_DOUBLE(hhi0,hhi1,hhi2, hlo0,hlo1,hlo2))
		{
			SUB78_3WORD_DOUBLE(hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlo1,hlo2);
			ADD78_3WORD_DOUBLE(hlo0,hlo1,hlo2, hhi0,hhi1,hhi2, hx0,hx1,hx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(hhi0,hhi1,hhi2, hlo0,hlo1,hlo2, hx0,hx1,hx2);
		}
		NORMALIZE78_3WORD_DOUBLE(hx0,hx1,hx2);

		if(CMPLT78_3WORD_DOUBLE(ihi0,ihi1,ihi2, ilo0,ilo1,ilo2))
		{
			SUB78_3WORD_DOUBLE(iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilo1,ilo2);
			ADD78_3WORD_DOUBLE(ilo0,ilo1,ilo2, ihi0,ihi1,ihi2, ix0,ix1,ix2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(ihi0,ihi1,ihi2, ilo0,ilo1,ilo2, ix0,ix1,ix2);
		}
		NORMALIZE78_3WORD_DOUBLE(ix0,ix1,ix2);

		if((pshift >> j) & (uint64)1)
		{
			/* Convert fx to uint96 for purpose of this particular comparison: */
			CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
			CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
			CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
			CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x0, qhalf0))
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
			}
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x1, qhalf1))
			{
				ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
				SUB78_3WORD_DOUBLE(gx0,gx1,gx2, gq0,gq1,gq2, gx0,gx1,gx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
			}
			NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x2, qhalf2))
			{
				ADD78_3WORD_DOUBLE(hx0,hx1,hx2, hx0,hx1,hx2, hx0,hx1,hx2);
				SUB78_3WORD_DOUBLE(hx0,hx1,hx2, hq0,hq1,hq2, hx0,hx1,hx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(hx0,hx1,hx2, hx0,hx1,hx2, hx0,hx1,hx2);
			}
			NORMALIZE78_3WORD_DOUBLE(hx0,hx1,hx2);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x3, qhalf3))
			{
				ADD78_3WORD_DOUBLE(ix0,ix1,ix2, ix0,ix1,ix2, ix0,ix1,ix2);
				SUB78_3WORD_DOUBLE(ix0,ix1,ix2, iq0,iq1,iq2, ix0,ix1,ix2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(ix0,ix1,ix2, ix0,ix1,ix2, ix0,ix1,ix2);
			}
			NORMALIZE78_3WORD_DOUBLE(ix0,ix1,ix2);
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
	CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
	CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
	CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);

	ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
	ADD96(x1,x1,x1);
	ADD96(x2,x2,x2);
	ADD96(x3,x3,x3);

	SUB96(x0,q0,x0);
	SUB96(x1,q1,x1);
	SUB96(x2,q2,x2);
	SUB96(x3,q3,x3);

	r = 0;
	if(CMPEQ96(x0, ONE96)) r +=  1;
	if(CMPEQ96(x1, ONE96)) r +=  2;
	if(CMPEQ96(x2, ONE96)) r +=  4;
	if(CMPEQ96(x3, ONE96)) r +=  8;
	return(r);
}

