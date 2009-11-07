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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef types_h_included
#define types_h_included

/* Include any needed level-0 header files: */
#include "platform.h"

/* Typedefs */

/*...useful utility parameters */

#undef TRUE
#define TRUE	1

#undef FALSE
#define FALSE	0

/* Basic integer types - we assume char/short/int mean 8/16/32 bits, respectively,
but this assumption gets checked at the start of program execution,
so we're not flying blind:
*/
#undef	 int8
#undef	sint8
#undef	uint8

#undef	 int16
#undef	sint16
#undef	uint16

#undef	 int32
#undef	sint32
#undef	uint32

#undef	 int64
#undef	sint64
#undef	uint64

#undef	 int64c
#undef	sint64c
#undef	uint64c

typedef          char		 int8;
typedef          char		sint8;
typedef unsigned char		uint8;

typedef          short		 int16;
typedef          short		sint16;
typedef unsigned short		uint16;

typedef          int		 int32;
typedef          int		sint32;
typedef unsigned int		uint32;

/* 64-bit int: */
/* MSVC doesn't like 'long long', and of course MS has their own
completely non-portable substitute:
*/
#ifdef OS_TYPE_WINDOWS
	typedef   signed __int64	 int64;
	typedef   signed __int64	sint64;
	typedef unsigned __int64	uint64;
	typedef const  signed __int64	 int64c;
	typedef const  signed __int64	sint64c;
	typedef const unsigned __int64	uint64c;
#else
	typedef          long long	 int64;
	typedef          long long	sint64;
	typedef unsigned long long	uint64;
	typedef const          long long	 int64c;
	typedef const          long long	sint64c;
	typedef const unsigned long long	uint64c;
#endif

/*******************************************************************************
   Some useful utility macros:
*******************************************************************************/

#undef	HERE
#define	HERE	__LINE__, __FILE__

/* Slow version of nearest-int; for doubles use faster trick,
but this is useful for reference purposes and error checking:
*/
#undef  NINT
#define NINT(x) floor(x + 0.5)

/* Fast double-float NINT. RND_A & RND_B declared in Mdata.h, defined in util.c:check_nbits_in_types .

***NOTE:*** The util.c functions set_fpu_params() and check_nbits_in_types()
MUST BE CALLED (in that order) AT PROGRAM INVOCATION THIS MACRO TO WORK PROPERLY!!!
*/
#undef  DNINT
#ifdef COMPILER_TYPE_ICC
	#define DNINT(x) rint((x))
#else
	#define DNINT(x) ((x) + RND_A) - RND_B
#endif

/* NOTE: when calling the MAX/MIN/ABS macros, do NOT allow either argument to be a
function call, since that may result in the function being called twice. */
#undef  MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#undef  MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#undef  ABS					/* 64-bit absval is not portable, so alias it. */
#define ABS(a) ((a) < 0 ? -(a) : (a))

#define IS_ODD(a)	( (int)(a) & 1)

#define IS_EVEN(a)	(~(int)(a) & 1)

#define	STREQ(s1,s2)	(!strcmp(s1,s2))

#define	STREQN(s1,s2,n)	(!strncmp(s1,s2,n))

#define	STRNEQ(s1,s2)	( strcmp(s1,s2))

#define	STRNEQN(s1,s2,n)( strncmp(s1,s2,n))

/* Complex multiplication C = A*B, input as pairs of doubles.
Designed so any or all of A, B, C may point to the same memory location. */
#define CMUL(ar,ai,br,bi,cr,ci)\
{\
	double __tmp = ar;\
	ci = __tmp*bi + ai*br;\
	cr = __tmp*br - ai*bi;\
}

/* Pairs-of-doubles structs: we'd like to give these all a GCC-style
"_attribute__ ((aligned (16)))" alignment flag, but as that's not portable,
we hope other compilers will be smart enough to properly align these
without external prompting.

Macros that manipulate complex (real, imag) pairs will have prefix CMPLX_
and will have declarations residing in the header file cmplx.h ;

Macros that manipulate non-complex pairs of doubles will have prefix VD_
and will have declarations residing in the header file vec_double.h ;

There will likely be significant overlap in the low-level functionality
(e.g. SSE2 intrinsics and data types) used by these 2 classes of macros;
which is why we union-ize the 128-bit structs storing the relevant data.
*/
#if 0	/*	#ifdef COMPILER_TYPE_GCC	*/

	#undef complex
	struct complex{
		double re;
		double im;
	} _attribute__ ((aligned (16)));	/* Gives 'error: parse error before '(' token' */

	#undef lohi_double_pair
	struct lohi_double_pair
	{
		double lo;
		double hi;
	} _attribute__ ((aligned (16)));

	#undef ddouble
	struct ddouble
	{
		double d[2];
	} _attribute__ ((aligned (16)));

#else

	#undef complex
	struct complex{
		double re;
		double im;
	};

	#undef lohi_double_pair
	struct lohi_double_pair
	{
		double lo;
		double hi;
	};

	#undef ddouble
	struct ddouble
	{
		double d[2];
	};

#endif

/* Unionize the paired doubles! Complex couples of the world unite, and stuff... */
union vec_double
{
	struct complex			 cmplx;
	struct lohi_double_pair lo_hi;
	struct ddouble			 dd;
};

/* 128-bit vector data types for AltiVec, Cell and SSE2/3:
Alas, we can't use the above typedefs here, i.e. "vector uint32" won't work:
*/
#if(CPU_HAS_ALTIVEC || CPU_IS_CELL)
	typedef	vector unsigned char	vec_uint8X16 ;
	typedef	vector unsigned short	vec_uint16X8 ;
	typedef	vector unsigned int		vec_uint32X4 ;
#endif
/* Key difference is that the Cell SPU also supports vector floating doubles, which AltiVec does not.
Note that Cell also supports vector long long ints, but as there is currently no hardware arithmetic
support for same, there's little point in using them.
*/
#if(CPU_IS_CELL)
	typedef	vector double			vec_double ;
#endif

/* For 96-bit ints, define a uint64+32 basic type, then declare a typedef
of that to uint192 so we can declare it like a uint64 (i.e. without explicitly
using the 'struct ...' declarator everywhere) in our code:
*/
/* 5/04/2005:
	OBSOLETE - prefer all these multiword ints to be a multiple of 64 bits long,
	e.g. uint96s are really uint128s with upper 32 bits zero, and so forth.

#undef uint64p32
struct uint64p32{
	uint64 d0;
	uint32 d1;
};

#undef uint96
typedef	struct uint64p32	uint96;
*/

/* Proceed analogously for 128-bit ints: */
#undef uint64x2
struct uint64x2{
	uint64 d0;
	uint64 d1;
};

#undef uint96
#undef uint128
typedef	struct uint64x2		uint96;
typedef	struct uint64x2		uint128;


/* For 192-bit ints, want to be able to access the low 2 words either as uint192.d0,d1
or as a uint128, so define uint64x3 and uint128+64 basic types, declare union192 as a
union of these, and also declare a typedef of uint64x3 to uint192 so we can declare it like
a uint64 (and with the same kind of subfield accessor names as a uint128) in our code:
*/
#undef uint64x3
struct uint64x3{
	uint64 d0;
	uint64 d1;
	uint64 d2;
};

#undef uint160
#undef uint192
typedef	struct uint64x3		uint160;
typedef	struct uint64x3		uint192;

#undef uint128p64
struct uint128p64{
	uint128 lo128;
	uint64  hi64;
};

#undef union192
union  union192{
	       uint192		u192;	/* By declaring this type first in the union, force inits in this form */
	struct uint128p64	u128p64;
};

#undef unio192
typedef	union union192		unio192;

/* 256-bit ints: */
#undef uint64x4
struct uint64x4{
	uint64 d0;
	uint64 d1;
	uint64 d2;
	uint64 d3;
};

#undef uint256
typedef	struct uint64x4		uint256;

/* Useful extern constants to export (defined in types.c): */

extern const uint96 NIL96;
extern const uint96 ONE96;
extern const uint96 TWO96;

extern const uint128 NIL128;
extern const uint128 ONE128;
extern const uint128 TWO128;

extern const uint160 NIL160;
extern const uint160 ONE160;
extern const uint160 TWO160;

extern const uint192 NIL192;
extern const uint192 ONE192;
extern const uint192 TWO192;

extern const uint256 NIL256;
extern const uint256 ONE256;
extern const uint256 TWO256;

#endif	/* types_h_included */
