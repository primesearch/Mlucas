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
#ifndef imul_macro1_h_included
#define imul_macro1_h_included

#include "imul_macro0.h"

#undef YES_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Macros for 96/128/160/192/256-bit unsigned integer add & subtract. No data type checking
or checking for carry/borrow out of the high word is performed, i.e. arithmetic
is modulo 2^96/128/160/192/256, respectively. In some cases (e.g. add/sub/shift/compare)
we alias 96/160-bit ops to the analogous 128/192-bit op for mnemonic convenience,
in others (esp. multiply) the 96/160-bit variant
genuinely takes advantage of the fact that operands are 96/160 bits long in a fashion
that is advantageous on at least an appreciable set of CPUs.
*/

/* Multiword Add: */
#define ADD128(__x, __y, __sum)\
{\
	uint64 __t;	/* Need this in case any or all of x, y, sum point to the same address. */\
	__t      = __x.d0;\
	__sum.d0 = __x.d0 + __y.d0;\
	__sum.d1 = __x.d1 + __y.d1 + (__sum.d0 < __t); /* Overflow into high word is checked here. */\
}
#define ADD96(__x, __y, __sum)\
{\
	ADD128(__x, __y, __sum);\
}

#define ADD192(__x, __y, __sum)\
{\
	uint64 __s,__t, __cy0;	/* Need these in case any or all of x, y, sum point to the same address. */\
	__s      = __x.d0;\
	__t      = __x.d1;\
	__sum.d0 = __x.d0 + __y.d0;	__cy0 = (__sum.d0 < __s);\
	__sum.d1 = __x.d1 + __y.d1;\
	__sum.d2 = __x.d2 + __y.d2 + (__sum.d1 < __t);\
	/* Need a separate add-and-overflow check for the carry out of the low word: */\
	__sum.d1 += __cy0;\
	__sum.d2 += (__sum.d1 < __cy0);\
}

#define ADD192_PTR(__x, __y, __sum)\
{\
	uint64 __s,__t, __cy0;	/* Need these in case any or all of x, y, sum point to the same address. */\
	__s      = __x->d0;\
	__t      = __x->d1;\
	__sum->d0 = __x->d0 + __y->d0;	__cy0 = (__sum->d0 < __s);\
	__sum->d1 = __x->d1 + __y->d1;\
	__sum->d2 = __x->d2 + __y->d2 + (__sum->d1 < __t);\
	/* Need a separate add-and-overflow check for the carry out of the low word: */\
	__sum->d1 += __cy0;\
	__sum->d2 += (__sum->d1 < __cy0);\
}

#define ADD160(__x, __y, __sum)\
{\
	DBG_ASSERT((__x.d2 >> 32) == 0,"ADD160: (__x.d2 >> 32) == 0");\
	DBG_ASSERT((__y.d2 >> 32) == 0,"ADD160: (__y.d2 >> 32) == 0");\
	ADD192(__x, __y, __sum);\
	__sum.d2 &= 0x00000000ffffffff;	/* In case of add need to take care to get proper mod-2^160 result */\
}

#define ADD256(__x, __y, __sum)\
{\
	uint64 __tmp, __cy;						\
											\
	__tmp = __x.d0+ __y.d0;					\
	__cy  = (__tmp < __y.d0);				\
	__sum.d0 = __tmp;						\
											\
	__tmp = __x.d1 + __cy;					\
	__cy  = (__tmp < __x.d1);				\
	__tmp = __tmp + __y.d1;					\
	__cy += (__tmp < __y.d1);				\
	__sum.d1 = __tmp;						\
											\
	__tmp = __x.d2 + __cy;					\
	__cy  = (__tmp < __x.d2);				\
	__tmp = __tmp + __y.d2;					\
	__cy += (__tmp < __y.d2);				\
	__sum.d2 = __tmp;						\
											\
	__sum.d3 = __x.d3 + __y.d3 + __cy;		\
}

#define ADD256_PTR(__x, __y, __sum)\
{\
	uint64 __tmp, __cy;						\
											\
	__tmp = __x->d0+ __y->d0;					\
	__cy  = (__tmp < __y->d0);				\
	__sum->d0 = __tmp;						\
											\
	__tmp = __x->d1 + __cy;					\
	__cy  = (__tmp < __x->d1);				\
	__tmp = __tmp + __y->d1;					\
	__cy += (__tmp < __y->d1);				\
	__sum->d1 = __tmp;						\
											\
	__tmp = __x->d2 + __cy;					\
	__cy  = (__tmp < __x->d2);				\
	__tmp = __tmp + __y->d2;					\
	__cy += (__tmp < __y->d2);				\
	__sum->d2 = __tmp;						\
											\
	__sum->d3 = __x->d3 + __y->d3 + __cy;		\
}

/* Subtract: */
#define SUB128(__x, __y, __dif)\
{\
	uint64 __t;	/* Need this in case any or all of x, y, dif point to the same address. */\
	__t      = __x.d0;\
	__dif.d0 = __x.d0 - __y.d0;\
	__dif.d1 = __x.d1 - __y.d1 - (__dif.d0 > __t); /* Borrow from high word is checked here. */\
}
#define SUB96(__x, __y, __dif)\
{\
	SUB128(__x, __y, __dif);\
}

#define SUB128B(__xhi, __xlo, __y, __dif)\
{\
	uint64 __t;	/* Need this in case any or all of x, y, dif point to the same address. */\
	__t      = __xlo;\
	__dif.d0 = __xlo - __y.d0;\
	__dif.d1 = __xhi - __y.d1 - (__dif.d0 > __t); /* Borrow from high word is checked here. */\
}
#define SUB96B(__xhi, __xlo, __y, __dif)\
{\
	SUB128B(__xhi, __xlo, __y, __dif);\
}

#define SUB128C(__x, __yhi, __ylo, __dif)\
{\
	uint64 __t;	/* Need this in case any or all of x, y, dif point to the same address. */\
	__t      = __x.d0;\
	__dif.d0 = __x.d0 - __ylo;\
	__dif.d1 = __x.d1 - __yhi - (__dif.d0 > __t); /* Borrow from high word is checked here. */\
}
#define SUB96C(__x, __yhi, __ylo, __dif)\
{\
	SUB128C(__x, __yhi, __ylo, __dif);\
}

#define SUB192(__x, __y, __dif)\
{\
	uint64 __s,__t, __cy0;	/* Need these in case any or all of x, y, sum point to the same address. */\
	__s      = __x.d0;\
	__t      = __x.d1;\
	__dif.d0 = __x.d0 - __y.d0;	__cy0 = (__dif.d0 > __s);\
	__dif.d1 = __x.d1 - __y.d1;\
	__dif.d2 = __x.d2 - __y.d2 - (__dif.d1 > __t);\
	/* Need a separate sub-and-underflow check for the borrow out of the low word: */\
	__t      = __dif.d1;	/* Need an extra temp here due to asymmetry of subtract */\
	__dif.d1 -= __cy0;\
	__dif.d2 -= (__dif.d1 > __t);\
}

#define SUB160(__x, __y, __dif)\
{\
	DBG_ASSERT((__x.d2 >> 32) == 0,"SUB160: (__x.d2 >> 32) == 0");\
	DBG_ASSERT((__y.d2 >> 32) == 0,"SUB160: (__y.d2 >> 32) == 0");\
	SUB192(__x, __y, __dif);\
	__dif.d2 &= 0x00000000ffffffff;	/* In case of add need to take care to get proper mod-2^160 result */\
}


#define SUB256(__x, __y, __dif)\
{\
	/* Need an extra temp here due to asymmetry of subtract: */\
	uint64 __tmp, __tmp2, __cy;				\
											\
	__tmp2= __x.d0 - __y.d0;				\
	__cy  = (__tmp2 > __x.d0);				\
	__dif.d0 = __tmp2;						\
											\
	__tmp = __x.d1 - __cy;					\
	__cy  = (__tmp > __x.d1);				\
	__tmp2= __tmp - __y.d1;					\
	__cy += (__tmp2 > __tmp);				\
	__dif.d1 = __tmp2;						\
											\
	__tmp = __x.d2 - __cy;					\
	__cy  = (__tmp > __x.d2);				\
	__tmp2= __tmp - __y.d2;					\
	__cy += (__tmp2 > __tmp);				\
	__dif.d2 = __tmp2;						\
											\
	__dif.d3 = __x.d3 - __y.d3 - __cy;		\
}


/* Macros for 96/128/160/192/256-bit unsigned integer compare: */
/* Unsigned < returns 1 if x < y, 0 otherwise: */
#define	CMPULT128(__x, __y)				((uint64)__x.d1 < (uint64)__y.d1 || ((uint64)__x.d1 == (uint64)__y.d1 && (uint64)__x.d0 < (uint64)__y.d0))
#define	CMPULT96(__x, __y)				((uint32)__x.d1 < (uint32)__y.d1 || ((uint32)__x.d1 == (uint32)__y.d1 && (uint64)__x.d0 < (uint64)__y.d0))

#define	CMPULT128B(__xhi, __xlo, __y)	((uint32)__xhi < (uint32)__y.d1 || ((uint32)__xhi == (uint32)__y.d1 && (uint64)__xlo < (uint64)__y.d0))
#define	CMPULT96B(__xhi, __xlo, __y)	((uint32)__xhi < (uint32)__y.d1 || ((uint32)__xhi == (uint32)__y.d1 && (uint64)__xlo < (uint64)__y.d0))

#define	CMPULT128C(__x, __yhi, __ylo)	((uint32)__x.d1 < (uint32)__yhi || ((uint32)__x.d1 == (uint32)__yhi && (uint64)__x.d0 < (uint64)__ylo))
#define	CMPULT96C(__x, __yhi, __ylo)	((uint32)__x.d1 < (uint32)__yhi || ((uint32)__x.d1 == (uint32)__yhi && (uint64)__x.d0 < (uint64)__ylo))

#define	CMPULT192(__x, __y)				((uint64)__x.d2 < (uint64)__y.d2 || ((uint64)__x.d2 == (uint64)__y.d2 && (uint64)__x.d1 < (uint64)__y.d1) || ((uint64)__x.d2 == (uint64)__y.d2 && (uint64)__x.d1 == (uint64)__y.d1 && (uint64)__x.d0 < (uint64)__y.d0))
#define	CMPULT160(__x, __y)				CMPULT192(__x, __y)

#define	CMPULT256(__x, __y)				((uint64)__x.d3 < (uint64)__y.d3 || ((uint64)__x.d3 == (uint64)__y.d3 && (uint64)__x.d2 < (uint64)__y.d2) || ((uint64)__x.d3 == (uint64)__y.d3 && (uint64)__x.d2 == (uint64)__y.d2 && (uint64)__x.d1 < (uint64)__y.d1) || ((uint64)__x.d3 == (uint64)__y.d3 && (uint64)__x.d2 == (uint64)__y.d2 && (uint64)__x.d1 == (uint64)__y.d1 && (uint64)__x.d0 < (uint64)__y.d0))

/* Unsigned > returns 1 if x > y, 0 otherwise: */
#define CMPUGT128(__x, __y)				CMPULT128(__y, __x)
#define CMPUGT96(__x, __y)				CMPULT96(__y, __x)

#define	CMPUGT128B(__xhi, __xlo, __y)	CMPULT128C(__y, __xhi, __xlo)
#define	CMPUGT96B(__xhi, __xlo, __y)	CMPULT96C(__y, __xhi, __xlo)

#define	CMPUGT128C(__x, __yhi, __ylo)	CMPULT128B(__yhi, __ylo, __x)
#define	CMPUGT96C(__x, __yhi, __ylo)	CMPULT96B(__yhi, __ylo, __x)

#define	CMPUGT192(__x, __y)				CMPULT192(__y, __x)
#define	CMPUGT160(__x, __y)				CMPULT192(__y, __x)

#define	CMPUGT256(__x, __y)				CMPULT256(__y, __x)

/* Unsigned <= returns 1 if x <= y, 0 otherwise: */
#define	CMPULE128(__x, __y)				!CMPUGT128(__x, __y)
#define	CMPULE96(__x, __y)				!CMPUGT96(__x, __y)

#define	CMPULE128B(__xhi, __xlo, __y)	!CMPUGT128B(__xhi, __xlo, __y)

#define	CMPULE192(__x, __y)				!CMPUGT192(__x, __y)
#define	CMPULE160(__x, __y)				CMPULE192(__x, __y)

#define	CMPULE256(__x, __y)				!CMPUGT256(__x, __y)

/* Unsigned >= returns 1 if x >= y, 0 otherwise: */
#define	CMPUGE128(__x, __y)				!CMPULT128(__x, __y)
#define	CMPUGE96(__x, __y)				!CMPULT96(__x, __y)

#define	CMPUGE128B(__xhi, __xlo, __y)	!CMPULT128B(__xhi, __xlo, __y)

#define	CMPUGE192(__x, __y)				!CMPULT192(__x, __y)
#define	CMPUGE160(__x, __y)				CMPUGE192(__x, __y)

#define	CMPUGE256(__x, __y)				!CMPULT256(__x, __y)

/* == returns 1 if x == y, 0 otherwise: */
#define	CMPEQ128(__x, __y)				(__x.d1 == __y.d1 && __x.d0 == __y.d0)
#define	CMPEQ96(__x, __y)				CMPEQ128(__x, __y)

#define	CMPEQ128B(__xhi, __xlo, __y)	(__xhi == __y.d1 && __xlo == __y.d0)

#define	CMPEQ192(__x, __y)				(__x.d2 == __y.d2 && __x.d1 == __y.d1 && __x.d0 == __y.d0)
#define	CMPEQ160(__x, __y)				CMPEQ192(__x, __y)

#define	CMPEQ256(__x, __y)				(__x.d3 == __y.d3 && __x.d2 == __y.d2 && __x.d1 == __y.d1 && __x.d0 == __y.d0)

/* Binary predicates for use of stdlib qsort(): Only support unsigned for multiword ints (for now): */
int ncmp_uint128(const void * a, const void * b);
int ncmp_uint192(const void * a, const void * b);
int ncmp_uint256(const void * a, const void * b);

/*
Macros for 96/128/160/192/256-bit unsigned integer shift. These are constructed so the
input (__x) and output (__y) arguments can point to the same object, if desired.
If the shift count (__n) is >= the width of the integer type, 0 is returned.
*/
/* Left-shifts: */
#define LSHIFT128(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >= 0,"LSHIFT128: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.d1 = ((uint64)__x.d1);\
		__y.d0 = ((uint64)__x.d0);\
	}\
	else if(__n < 64)\
	{\
		__y.d1 = ((uint64)__x.d1 << __n) + ((uint64)__x.d0 >> (64-__n));\
		__y.d0 = ((uint64)__x.d0 << __n);\
	}\
	else if(__n < 128)\
	{\
		__y.d1 = ((uint64)__x.d0 << (__n-64));\
		__y.d0 = (uint64)0;\
	}\
	else\
	{\
		__y.d1 = (uint64)0;\
		__y.d0 = (uint64)0;\
	}\
}
#define LSHIFT96(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >= 0,"LSHIFT96: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.d1 = ((uint64)__x.d1);\
		__y.d0 = ((uint64)__x.d0);\
	}\
	else if(__n < 32)\
	{\
		__y.d1 = ((uint32)__x.d1 << __n) + (uint32)((uint64)__x.d0 >> (64-__n));\
		__y.d0 = ((uint64)__x.d0 << __n);\
	}\
	else if(__n < 64)\
	{\
		__y.d1 =                           (uint32)((uint64)__x.d0 >> (64-__n));\
		__y.d0 = ((uint64)__x.d0 << __n);\
	}\
	else if(__n < 96)\
	{\
		__y.d1 = (uint32)((uint64)__x.d0 << (__n-64));\
		__y.d0 = (uint64)0;\
	}\
	else\
	{\
		__y.d1 = (uint32)0;\
		__y.d0 = (uint64)0;\
	}\
}

#define LSHIFT192(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >= 0,"LSHIFT192: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.d2 = ((uint64)__x.d2);\
		__y.d1 = ((uint64)__x.d1);\
		__y.d0 = ((uint64)__x.d0);\
	}\
	else if(__n < 64)\
	{\
		__y.d2 = ((uint64)__x.d2 << __n) + ((uint64)__x.d1 >> (64-__n));\
		__y.d1 = ((uint64)__x.d1 << __n) + ((uint64)__x.d0 >> (64-__n));\
		__y.d0 = ((uint64)__x.d0 << __n);\
	}\
	else if(__n == 64)/* Need an extra n==64 case here due to >> (128-n) below: */\
	{\
		__y.d2 = (uint64)__x.d1;\
		__y.d1 = (uint64)__x.d0;\
		__y.d0 = (uint64)0;\
	}\
	else if(__n < 128)\
	{\
		__y.d2 = ((uint64)__x.d1 << (__n-64)) + ((uint64)__x.d0 >> (128-__n));\
		__y.d1 = ((uint64)__x.d0 << (__n-64));\
		__y.d0 = (uint64)0;\
	}\
	else if(__n < 192)\
	{\
		__y.d2 = ((uint64)__x.d0 << (__n-128));\
		__y.d1 = (uint64)0;\
		__y.d0 = (uint64)0;\
	}\
	else\
	{\
		__y.d2 = (uint64)0;\
		__y.d1 = (uint64)0;\
		__y.d0 = (uint64)0;\
	}\
}
#define LSHIFT160(__x, __n, __y)\
{\
	DBG_ASSERT(((uint64)__x.d2 >> 32) == 0,"LSHIFT160: ((uint64)__x.d2 >> 32) == 0");\
	LSHIFT192(__x,__n, __y);\
	__y.d2 &= 0x00000000ffffffff;\
}

#define LSHIFT256(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >= 0,"LSHIFT256: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.d3 = ((uint64)__x.d3);\
		__y.d2 = ((uint64)__x.d2);\
		__y.d1 = ((uint64)__x.d1);\
		__y.d0 = ((uint64)__x.d0);\
	}\
	else if(__n < 64)\
	{\
		__y.d3 = ((uint64)__x.d3 << __n) + ((uint64)__x.d2 >> (64-__n));\
		__y.d2 = ((uint64)__x.d2 << __n) + ((uint64)__x.d1 >> (64-__n));\
		__y.d1 = ((uint64)__x.d1 << __n) + ((uint64)__x.d0 >> (64-__n));\
		__y.d0 = ((uint64)__x.d0 << __n);\
	}\
	else if(__n == 64)/* Need an extra n==64 case here due to >> (128-n) below: */\
	{\
		__y.d3 = (uint64)__x.d2;\
		__y.d2 = (uint64)__x.d1;\
		__y.d1 = (uint64)__x.d0;\
		__y.d0 = (uint64)0;\
	}\
	else if(__n < 128)\
	{\
		__y.d3 = ((uint64)__x.d2 << (__n-64)) + ((uint64)__x.d1 >> (128-__n));\
		__y.d2 = ((uint64)__x.d1 << (__n-64)) + ((uint64)__x.d0 >> (128-__n));\
		__y.d1 = ((uint64)__x.d0 << (__n-64));\
		__y.d0 = (uint64)0;\
	}\
	else if(__n ==128)/* Need an extra n==128 case here due to >> (192-n) below: */\
	{\
		__y.d3 = (uint64)__x.d1;\
		__y.d2 = (uint64)__x.d0;\
		__y.d1 = (uint64)0;\
		__y.d0 = (uint64)0;\
	}\
	else if(__n < 192)\
	{\
		__y.d3 = ((uint64)__x.d1 << (__n-128)) + ((uint64)__x.d0 >> (192-__n));\
		__y.d2 = ((uint64)__x.d0 << (__n-128));\
		__y.d1 = (uint64)0;\
		__y.d0 = (uint64)0;\
	}\
	else if(__n < 256)\
	{\
		__y.d3 = ((uint64)__x.d0 << (__n-192));\
		__y.d2 = (uint64)0;\
		__y.d1 = (uint64)0;\
		__y.d0 = (uint64)0;\
	}\
	else\
	{\
		__y.d3 = (uint64)0;\
		__y.d2 = (uint64)0;\
		__y.d1 = (uint64)0;\
		__y.d0 = (uint64)0;\
	}\
}


/* (Logical) Right-shifts: */
#define RSHIFT128(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >= 0,"RSHIFT128: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.d0 = ((uint64)__x.d0);\
		__y.d1 = ((uint64)__x.d1);\
	}\
	else if(__n < 64)\
	{\
		__y.d0 = ((uint64)__x.d0 >> __n) + ((uint64)__x.d1 << (64-__n));\
		__y.d1 = ((uint64)__x.d1 >> __n);\
	}\
	else if(__n < 128)\
	{\
		__y.d0 = ((uint64)__x.d1 >> (__n-64));\
		__y.d1 = (uint64)0;\
	}\
	else\
	{\
		__y.d0 = (uint64)0;\
		__y.d1 = (uint64)0;\
	}\
}
#define RSHIFT96(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >= 0,"RSHIFT96: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.d0 = ((uint64)__x.d0);\
		__y.d1 = ((uint32)__x.d1);\
	}\
	else if(__n < 32)\
	{\
		__y.d0 = ((uint64)__x.d0 >> __n) + ((uint64)__x.d1 << (64-__n));\
		__y.d1 = ((uint32)__x.d1 >> __n);\
	}\
	else if(__n < 64)\
	{\
		__y.d0 = ((uint64)__x.d0 >> __n) + ((uint64)__x.d1 << (64-__n));\
		__y.d1 = (uint32)0;\
	}\
	else if(__n < 96)\
	{\
		__y.d0 = ((uint64)__x.d1 >> (__n-64));\
		__y.d1 = (uint32)0;\
	}\
	else\
	{\
		__y.d0 = (uint64)0;\
		__y.d1 = (uint32)0;\
	}\
}

/* The store-mod-64-shift-counts-in-locals stuff here is to get GCC to finally STFU by silencing the following spurious warnings:
	warning: right shift count >= width of type
	warning: right shift count is negative
*/
#define RSHIFT192(__x, __n, __y)\
{\
	int __lsh,__rsh;\
	DBG_ASSERT((int64)__n >= 0,"RSHIFT192: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.d0 = ((uint64)__x.d0);\
		__y.d1 = ((uint64)__x.d1);\
		__y.d2 = ((uint64)__x.d2);\
	}\
	else if(__n < 64)\
	{\
		__lsh = 64-__n;\
		__rsh = __n&63;\
		__y.d0 = ((uint64)__x.d0 >> __rsh) + ((uint64)__x.d1 << __lsh);\
		__y.d1 = ((uint64)__x.d1 >> __rsh) + ((uint64)__x.d2 << __lsh);\
		__y.d2 = ((uint64)__x.d2 >> __rsh);\
	}\
	else if(__n == 64)/* Need an extra n==64 case here due to << (128-n) below: */\
	{\
		__y.d0 = (uint64)__x.d1;\
		__y.d1 = (uint64)__x.d2;\
		__y.d2 = (uint64)0;\
	}\
	else if(__n < 128)\
	{\
		__lsh = 128-__n;\
		__rsh = (__n- 64)&63;\
		__y.d0 = ((uint64)__x.d1 >> __rsh) + ((uint64)__x.d2 << __lsh);\
		__y.d1 = ((uint64)__x.d2 >> __rsh);\
		__y.d2 = (uint64)0;\
	}\
	else if(__n < 192)\
	{\
		__rsh = (__n-128)&63;\
		__y.d0 = ((uint64)__x.d2 >> __rsh);\
		__y.d1 = (uint64)0;\
		__y.d2 = (uint64)0;\
	}\
	else\
	{\
		__y.d0 = (uint64)0;\
		__y.d1 = (uint64)0;\
		__y.d2 = (uint64)0;\
	}\
}
#define RSHIFT160(__x, __n, __y)\
{\
	DBG_ASSERT(((uint64)__x.d2 >> 32) == 0,"RSHIFT160: ((uint64)__x.d2 >> 32) == 0");\
	RSHIFT192(__x,__n, __y);\
	DBG_ASSERT((uint64)__y.d2 <= (uint64)__x.d2,"RSHIFT160: (uint64)__y.d2 <= (uint64)__x.d2");\
}

#define RSHIFT256(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >= 0,"RSHIFT256: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y.d0 = ((uint64)__x.d0);\
		__y.d1 = ((uint64)__x.d1);\
		__y.d2 = ((uint64)__x.d2);\
		__y.d3 = ((uint64)__x.d3);\
	}\
	else if(__n < 64)\
	{\
		__y.d0 = ((uint64)__x.d0 >> __n) + ((uint64)__x.d1 << (64-__n));\
		__y.d1 = ((uint64)__x.d1 >> __n) + ((uint64)__x.d2 << (64-__n));\
		__y.d2 = ((uint64)__x.d2 >> __n) + ((uint64)__x.d3 << (64-__n));\
		__y.d3 = ((uint64)__x.d3 >> __n);\
	}\
	else if(__n == 64)/* Need an extra n==64 case here due to << (128-n) below: */\
	{\
		__y.d0 = (uint64)__x.d1;\
		__y.d1 = (uint64)__x.d2;\
		__y.d2 = (uint64)__x.d3;\
		__y.d3 = (uint64)0;\
	}\
	else if(__n < 128)\
	{\
		__y.d0 = ((uint64)__x.d1 >> (__n-64)) + ((uint64)__x.d2 << (128-__n));\
		__y.d1 = ((uint64)__x.d2 >> (__n-64)) + ((uint64)__x.d3 << (128-__n));\
		__y.d2 = ((uint64)__x.d3 >> (__n-64));\
		__y.d3 = (uint64)0;\
	}\
	else if(__n ==128)/* Need an extra n==128 case here due to << (192-n) below: */\
	{\
		__y.d0 = (uint64)__x.d2;\
		__y.d1 = (uint64)__x.d3;\
		__y.d2 = (uint64)0;\
		__y.d3 = (uint64)0;\
	}\
	else if(__n < 192)\
	{\
		__y.d0 = ((uint64)__x.d2 >> (__n-128)) + ((uint64)__x.d3 << (192-__n));\
		__y.d1 = ((uint64)__x.d3 >> (__n-128));\
		__y.d2 = (uint64)0;\
		__y.d3 = (uint64)0;\
	}\
	else if(__n < 256)\
	{\
		__y.d0 = ((uint64)__x.d3 >> (__n-192));\
		__y.d1 = (uint64)0;\
		__y.d2 = (uint64)0;\
		__y.d3 = (uint64)0;\
	}\
	else\
	{\
		__y.d0 = (uint64)0;\
		__y.d1 = (uint64)0;\
		__y.d2 = (uint64)0;\
		__y.d3 = (uint64)0;\
	}\
}

/* Special versions that assume the shift count is strictly in [1, 63] */
/* Left-shifts: */
#define LSHIFT_FAST128(__x, __n, __y)\
{\
	DBG_ASSERT((uint64)__n != 0,"LSHIFT_FAST128: (uint64)__n != 0");\
	DBG_ASSERT((uint64)__n < 64,"LSHIFT_FAST128: (uint64)__n < 64");\
	__y.d1 = ((uint64)__x.d1 << __n) + ((uint64)__x.d0 >> (64-__n));\
	__y.d0 = ((uint64)__x.d0 << __n);\
}
#define LSHIFT_FAST96(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >  0,"LSHIFT96: (int64)__n >  0");\
	DBG_ASSERT((int64)__n < 32,"LSHIFT96: (int64)__n < 32");\
	__y.d1 = ((uint32)__x.d1 << __n) + (uint32)((uint64)__x.d0 >> (64-__n));\
	__y.d0 = ((uint64)__x.d0 << __n);\
}

#define LSHIFT_FAST192(__x, __n, __y)\
{\
	DBG_ASSERT((uint64)__n != 0,"LSHIFT_FAST192: (uint64)__n != 0");\
	DBG_ASSERT((uint64)__n < 64,"LSHIFT_FAST192: (uint64)__n < 64");\
	__y.d2 = ((uint64)__x.d2 << __n) + ((uint64)__x.d1 >> (64-__n));\
	__y.d1 = ((uint64)__x.d1 << __n) + ((uint64)__x.d0 >> (64-__n));\
	__y.d0 = ((uint64)__x.d0 << __n);\
}
#define LSHIFT_FAST160(__x, __n, __y)\
{\
	DBG_ASSERT(((uint64)__x.d2 >> 32) == 0,"LSHIFT_FAST160: ((uint64)__x.d2 >> 32) == 0");\
	LSHIFT_FAST192(__x,__n, __y);\
	__y.d2 &= 0x00000000ffffffff;\
}

#define LSHIFT_FAST256(__x, __n, __y)\
{\
	DBG_ASSERT((uint64)__n != 0,"LSHIFT_FAST256: (uint64)__n != 0");\
	DBG_ASSERT((uint64)__n < 64,"LSHIFT_FAST256: (uint64)__n < 64");\
	__y.d3 = ((uint64)__x.d3 << __n) + ((uint64)__x.d2 >> (64-__n));\
	__y.d2 = ((uint64)__x.d2 << __n) + ((uint64)__x.d1 >> (64-__n));\
	__y.d1 = ((uint64)__x.d1 << __n) + ((uint64)__x.d0 >> (64-__n));\
	__y.d0 = ((uint64)__x.d0 << __n);\
}

/* (Logical) Right-shifts: */
#define RSHIFT_FAST128(__x, __n, __y)\
{\
	DBG_ASSERT((uint64)__n != 0,"RSHIFT_FAST128: (uint64)__n != 0");\
	DBG_ASSERT((uint64)__n < 64,"RSHIFT_FAST128: (uint64)__n < 64");\
	__y.d0 = ((uint64)__x.d0 >> __n) + ((uint64)__x.d1 << (64-__n));\
	__y.d1 = ((uint64)__x.d1 >> __n);\
}
#define RSHIFT_FAST96(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >  0,"RSHIFT96: (int64)__n >  0");\
	DBG_ASSERT((int64)__n < 32,"RSHIFT96: (int64)__n < 32");\
	__y.d0 = ((uint64)__x.d0 >> __n) + ((uint64)__x.d1 << (64-__n));\
	__y.d1 = ((uint32)__x.d1 >> __n);\
}

#define RSHIFT_FAST192(__x, __n, __y)\
{\
	DBG_ASSERT((uint64)__n != 0,"RSHIFT_FAST192: (uint64)__n != 0");\
	DBG_ASSERT((uint64)__n < 64,"RSHIFT_FAST192: (uint64)__n < 64");\
	__y.d0 = ((uint64)__x.d0 >> __n) + ((uint64)__x.d1 << (64-__n));\
	__y.d1 = ((uint64)__x.d1 >> __n) + ((uint64)__x.d2 << (64-__n));\
	__y.d2 = ((uint64)__x.d2 >> __n);\
}
#define RSHIFT_FAST160(__x, __n, __y)\
{\
	DBG_ASSERT(((uint64)__x.d2 >> 32) == 0,"RSHIFT_FAST160: ((uint64)__x.d2 >> 32) == 0");\
	RSHIFT_FAST192(__x,__n, __y);\
	__y.d2 &= 0x00000000ffffffff;\
}

#define RSHIFT_FAST256(__x, __n, __y)\
{\
	DBG_ASSERT((uint64)__n != 0,"RSHIFT_FAST256: (uint64)__n != 0");\
	DBG_ASSERT((uint64)__n < 64,"RSHIFT_FAST256: (uint64)__n < 64");\
	__y.d0 = ((uint64)__x.d0 >> __n) + ((uint64)__x.d1 << (64-__n));\
	__y.d1 = ((uint64)__x.d1 >> __n) + ((uint64)__x.d2 << (64-__n));\
	__y.d2 = ((uint64)__x.d2 >> __n) + ((uint64)__x.d3 << (64-__n));\
	__y.d3 = ((uint64)__x.d3 >> __n);\
}


/* Macros for 96/128/160/192-bit unsigned integer leading-zeros counting. For > 192-bit, call mi64_leadz instead.
Cast the result of the high-part-equals-zero test to a signed 32-bit (-1) because leadz* returns a 32-bit value for inputs of all sizes: */
#define LEADZ128(__x)	( leadz64(__x.d1) + ((-(sint32)(__x.d1 == 0)) && leadz64(__x.d0)) )
#define LEADZ96(__x)	( leadz32(__x.d1) + ((-(sint32)(__x.d1 == 0)) && leadz64(__x.d0)) )

#define LEADZ192(__x)	( leadz64(__x.d2) + ((-(sint32)(__x.d2 == 0)) && leadz64(__x.d1))+ ((-(sint32)(__x.d2 == 0 && __x.d1 == 0)) && leadz64(__x.d0)) )
#define LEADZ160(__x)	LEADZ192(__x)

/* For larger operands it's more convenient to use a conditional + leadz64: */
#define LEADZ256(__x)	leadz256(__x)

/*** Special pointer-based versions of key 96-bit macros, for the pure-ASM int64 code *****/
#define	CMPULT96_PTR(__x, __y)	((uint32)__x->d1 < (uint32)__y->d1 || ((uint32)__x->d1 == (uint32)__y->d1 && (uint64)__x->d0 < (uint64)__y->d0))
#define CMPUGT96_PTR(__x, __y)	CMPULT96_PTR(__y, __x)
#define	CMPEQ96_PTR(__x, __y)	((uint32)__x->d1 == (uint32)__y->d1 && (uint64)__x->d0 == (uint64)__y->d0)

#define RSHIFT_FAST96_PTR(__x, __n, __y)\
{\
	__y->d0 = ((uint64)__x->d0 >> __n) + ((uint64)__x->d1 << (64-__n));\
	__y->d1 = ((uint32)__x->d1 >> __n);\
}

#define LSHIFT96_PTR(__x, __n, __y)\
{\
	DBG_ASSERT((int64)__n >= 0,"LSHIFT96_PTR: (int64)__n >= 0");\
	/* Need to handle zero shift count separately: */\
	if(__n == 0)\
	{\
		__y->d1 = ((uint64)__x->d1);\
		__y->d0 = ((uint64)__x->d0);\
	}\
	else if(__n < 32)\
	{\
		__y->d1 = ((uint32)__x->d1 << __n) + (uint32)((uint64)__x->d0 >> (64-__n));\
		__y->d0 = ((uint64)__x->d0 << __n);\
	}\
	else if(__n < 64)\
	{\
		__y->d1 =                           (uint32)((uint64)__x->d0 >> (64-__n));\
		__y->d0 = ((uint64)__x->d0 << __n);\
	}\
	else if(__n < 96)\
	{\
		__y->d1 = (uint32)((uint64)__x->d0 << (__n-64));\
		__y->d0 = (uint64)0;\
	}\
	else\
	{\
		__y->d1 = (uint32)0;\
		__y->d0 = (uint64)0;\
	}\
}

#define ADD96_PTR(__x, __y, __sum)\
{\
	uint64 __t;	/* Need this in case any or all of x, y, sum point to the same address. */\
	__t       = __x->d0;\
	__sum->d0 = __x->d0 + __y->d0;\
	__sum->d1 = __x->d1 + __y->d1 + (__sum->d0 < __t); /* Overflow into high word is checked here. */\
}

#define SUB96_PTR(__x, __y, __dif)\
{\
	uint64 __t;	/* Need this in case any or all of x, y, dif point to the same address. */\
	__t       = __x->d0;\
	__dif->d0 = __x->d0 - __y->d0;\
	__dif->d1 = __x->d1 - __y->d1 - (__dif->d0 > __t); /* Borrow from high word is checked here. */\
}


#define MULH96_PTR_q4(\
  __x0, __y0, __hi0\
, __x1, __y1, __hi1\
, __x2, __y2, __hi2\
, __x3, __y3, __hi3)\
{\
	uint32 __l32_0, __h32_0;\
	uint32 __l32_1, __h32_1;\
	uint32 __l32_2, __h32_2;\
	uint32 __l32_3, __h32_3;\
	uint64 __a0,__b0,__c0,__d0,__m0,__h0;\
	uint64 __a1,__b1,__c1,__d1,__m1,__h1;\
	uint64 __a2,__b2,__c2,__d2,__m2,__h2;\
	uint64 __a3,__b3,__c3,__d3,__m3,__h3;\
	\
	MUL_LOHI32(__x0->d1, __y0->d1, __l32_0, __h32_0);\
	MUL_LOHI32(__x1->d1, __y1->d1, __l32_1, __h32_1);\
	MUL_LOHI32(__x2->d1, __y2->d1, __l32_2, __h32_2);\
	MUL_LOHI32(__x3->d1, __y3->d1, __l32_3, __h32_3);\
	\
	__h0 = __l32_0 + ((uint64)__h32_0 << 32);\
	__h1 = __l32_1 + ((uint64)__h32_1 << 32);\
	__h2 = __l32_2 + ((uint64)__h32_2 << 32);\
	__h3 = __l32_3 + ((uint64)__h32_3 << 32);\
	\
	MULH64(__x0->d0,__y0->d0,       __m0);\
	MULH64(__x1->d0,__y1->d0,       __m1);\
	MULH64(__x2->d0,__y2->d0,       __m2);\
	MULH64(__x3->d0,__y3->d0,       __m3);\
	\
	MUL64x32(__x0->d0,__y0->d1, __a0, __b0);\
	MUL64x32(__x1->d0,__y1->d1, __a1, __b1);\
	MUL64x32(__x2->d0,__y2->d1, __a2, __b2);\
	MUL64x32(__x3->d0,__y3->d1, __a3, __b3);\
	\
	MUL64x32(__y0->d0,__x0->d1, __c0, __d0);\
	MUL64x32(__y1->d0,__x1->d1, __c1, __d1);\
	MUL64x32(__y2->d0,__x2->d1, __c2, __d2);\
	MUL64x32(__y3->d0,__x3->d1, __c3, __d3);\
	\
	__m0 += __a0;\
	__m1 += __a1;\
	__m2 += __a2;\
	__m3 += __a3;\
	\
	__h0 += __b0 + (__m0 < __a0);\
	__h1 += __b1 + (__m1 < __a1);\
	__h2 += __b2 + (__m2 < __a2);\
	__h3 += __b3 + (__m3 < __a3);\
	\
	__m0 += __c0;\
	__m1 += __c1;\
	__m2 += __c2;\
	__m3 += __c3;\
	\
	__h0 += __d0 + (__m0 < __c0);\
	__h1 += __d1 + (__m1 < __c1);\
	__h2 += __d2 + (__m2 < __c2);\
	__h3 += __d3 + (__m3 < __c3);\
	\
	/* Now put upper half of the result into __hi: */\
	__hi0->d0 = (__m0>>32) + (__h0<<32);	__hi0->d1 = (uint32)(__h0>>32);\
	__hi1->d0 = (__m1>>32) + (__h1<<32);	__hi1->d1 = (uint32)(__h1>>32);\
	__hi2->d0 = (__m2>>32) + (__h2<<32);	__hi2->d1 = (uint32)(__h2>>32);\
	__hi3->d0 = (__m3>>32) + (__h3<<32);	__hi3->d1 = (uint32)(__h3>>32);\
}


/****************************************************************************************************************************
*																															*
* Macros for 96/128/160/192/256-bit unsigned integer multiply.																		*
*																															*
* Doing 2^p mod q using a Montgomery-style mod needs a SQR_LOHI, MULL64 and MULH for each bit of p processed,				*
* so a rough opcount (and exact MUL opcount) can be gotten by summing the individual opcounts for these three operations.	*
*																															*
* For example, using the Alpha-style 64-bit integer MUL model:																*
*	Using uint96_80   math, total opcount for a modsquare is  ?+ ?+ ? = ?? MUL,  ?+ ?+?? = ?? ALU ops. 						*
*	Using uint96      math, total opcount for a modsquare is  5+ 4+ 6 = 15 MUL,  9+ 3+12 = 23 ALU ops. 						*
*	Using uint128_96  math, total opcount is                  5+ 4+ 7 = 16 MUL,  5+ 2+12 = 19 ALU ops.                    	*
*	Using uint128     math, total opcount is                  6+ 4+ 7 = 17 MUL, 12+ 2+12 = 26 ALU ops.                    	*
*	Using uint160     math, total opcount is                  9+ 9+11 = 29 MUL, 23+10+22 = 55 ALU ops.                    	*
*	Using uint192     math, total opcount is                 12+ 9+17 = 38 MUL, 36+10+36 = 82 ALU ops.                    	*
*																															*
****************************************************************************************************************************/

/*****************************************************/
/*                     96 x 96                       */
/*****************************************************/

/* 192-bit square of uint96 *x. Result is returned in a pair of uint96s.

On Alpha, this needs a total of 5 MUL instructions and 9 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define SQR_LOHI96(__x, __lo, __hi)\
    {\
		uint64 __l,__m,__h,__a,__b,__t;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"SQR_LOHI96: (__x.d1 >> 32) == 0");\
		__t   = (uint64)(__x.d1);\
		__h   = __t*__t;\
		SQR_LOHI64(__x.d0,    &__l,&__m);\
		MUL64x32(__x.d0,__t,&__a,&__b);\
		__b = (__b << 1) + (__a >> 63);	__a <<= 1;	/* Double the cross term */\
		__m += __a;\
		__h += __b + (__m < __a); /* Overflow into high word is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __l;						__lo.d1 =  __m & 0x00000000ffffffff;\
		__hi.d0 = (__m>>32) + (__h<<32);	__hi.d1 = (__h>>32);\
    }
	/* Want a 95-bit version in this case, since the "pipelined" version uses the
	single-operand versions: */
    #define SQR_LOHI95(__x, __lo, __hi)\
    {\
		uint64 __l,__m,__h,__a,__b,__t;\
		\
		DBG_ASSERT((__x.d1 >> 31) == 0,"SQR_LOHI95: (__x.d1 >> 31) == 0");\
		__t   = (uint64)(__x.d1);\
		__h   = __t*__t;\
		SQR_LOHI64(__x.d0,         &__l,&__m);\
		MUL64x32(__x.d0,__t << 1,&__a,&__b);\
		__m += __a;\
		__h += __b + (__m < __a); /* Overflow into high word is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __l;						__lo.d1 =  __m & 0x00000000ffffffff;\
		__hi.d0 = (__m>>32) + (__h<<32);	__hi.d1 = (__h>>32);\
    }

    /* 4-operand-pipelined version: */
    #define SQR_LOHI96_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
		SQR_LOHI96(__x0, __y0, __lo0);\
		SQR_LOHI96(__x1, __y1, __lo1);\
		SQR_LOHI96(__x2, __y2, __lo2);\
		SQR_LOHI96(__x3, __y3, __lo3);\
	}
    /* 4-operand-pipelined version, specialized to 95-bit inputs: */
    #define SQR_LOHI95_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
		SQR_LOHI95(__x0, __y0, __lo0);\
		SQR_LOHI95(__x1, __y1, __lo1);\
		SQR_LOHI95(__x2, __y2, __lo2);\
		SQR_LOHI95(__x3, __y3, __lo3);\
	}

    /* 8-operand-pipelined version: */
    #define SQR_LOHI96_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
		SQR_LOHI96(__x0, __y0, __lo0);\
		SQR_LOHI96(__x1, __y1, __lo1);\
		SQR_LOHI96(__x2, __y2, __lo2);\
		SQR_LOHI96(__x3, __y3, __lo3);\
		SQR_LOHI96(__x4, __y4, __lo4);\
		SQR_LOHI96(__x5, __y5, __lo5);\
		SQR_LOHI96(__x6, __y6, __lo6);\
		SQR_LOHI96(__x7, __y7, __lo7);\
	}
    /* 8-operand-pipelined version, specialized to 95-bit inputs: */
    #define SQR_LOHI95_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
		SQR_LOHI95(__x0, __y0, __lo0);\
		SQR_LOHI95(__x1, __y1, __lo1);\
		SQR_LOHI95(__x2, __y2, __lo2);\
		SQR_LOHI95(__x3, __y3, __lo3);\
		SQR_LOHI95(__x4, __y4, __lo4);\
		SQR_LOHI95(__x5, __y5, __lo5);\
		SQR_LOHI95(__x6, __y6, __lo6);\
		SQR_LOHI95(__x7, __y7, __lo7);\
	}

#else

    #define SQR_LOHI96(__x, __lo, __hi)\
    {\
		uint64 __l,__m,__h,__a,__b;\
		uint32 __tt = __x.d1, __hl32,__hh32;\
		DBG_ASSERT((__x.d1 >> 32) == 0,"SQR_LOHI96: (__x.d1 >> 32) == 0");\
		MUL64x32(__x.d0,__tt, __a, __b);\
		SQR_LOHI64(__x.d0,     __l, __m);\
		MUL_LOHI32(__tt,__tt,__hl32,__hh32);\
		__b = (__b << 1) + (__a >> 63);	__a <<= 1;	/* Double the cross term */\
		__h = (uint64)(__hl32) + ((uint64)__hh32 << 32);\
		__m += __a;\
		__h += __b + (__m < __a); /* Overflow into high word is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __l;						__lo.d1 = (uint32)__m;\
		__hi.d0 = (__m>>32) + (__h<<32);	__hi.d1 = (__h>>32);\
    }

    /* 4-operand-pipelined version: designed so low half of output may be returned
    in input argument x, if x and lo happen to point to the same addresses:
    */
  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI96_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI96_q4: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI96_q4: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI96_q4: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI96_q4: (__x3.d1 >> 32) == 0");\
		\
		SQR_LOHI    (__x0.d0,                             __a0 ,     __b0 );\
		SQR_LOHI    (__x1.d0,                             __a1 ,     __b1 );\
		SQR_LOHI    (__x2.d0,                             __a2 ,     __b2 );\
		SQR_LOHI    (__x3.d0,                             __a3 ,     __b3 );\
		\
		MUL_LOHI64_ADD(__x0.d0, (uint64)__x0.d1 << 1,     __b0,     __b0 , __hi0.d0 );\
		MUL_LOHI64_ADD(__x1.d0, (uint64)__x1.d1 << 1,     __b1,     __b1 , __hi1.d0 );\
		MUL_LOHI64_ADD(__x2.d0, (uint64)__x2.d1 << 1,     __b2,     __b2 , __hi2.d0 );\
		MUL_LOHI64_ADD(__x3.d0, (uint64)__x3.d1 << 1,     __b3,     __b3 , __hi3.d0 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		__hi0.d0 += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1.d0 += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2.d0 += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3.d0 += (uint64)__x3.d1*(uint64)__x3.d1;\
		\
		__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d1 = __b3 & 0x00000000ffffffff;\
		\
		__hi0.d1 = (__hi0.d0 >> 32);\
		__hi1.d1 = (__hi1.d0 >> 32);\
		__hi2.d1 = (__hi2.d0 >> 32);\
		__hi3.d1 = (__hi3.d0 >> 32);\
		\
		__hi0.d0 = (__hi0.d0 << 32) + (__b0 >> 32);\
		__hi1.d0 = (__hi1.d0 << 32) + (__b1 >> 32);\
		__hi2.d0 = (__hi2.d0 << 32) + (__b2 >> 32);\
		__hi3.d0 = (__hi3.d0 << 32) + (__b3 >> 32);\
    }
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI95_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
	{\
		SQR_LOHI96_q4(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3)\
	}

  #elif(defined(MULH64_FAST))	/* Fast 64-bit MUL but no integer fused MADD */

    #define SQR_LOHI96_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__s0,__s1,__s2,__s3;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI96_q4: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI96_q4: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI96_q4: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI96_q4: (__x3.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,                    __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                    __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                    __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                    __a3 ,     __b3 );\
		\
		MUL_LOHI64(__x0.d0, (uint64)__x0.d1 << 1 ,     __s0 ,    __hi0.d0 );\
		MUL_LOHI64(__x1.d0, (uint64)__x1.d1 << 1 ,     __s1 ,    __hi1.d0 );\
		MUL_LOHI64(__x2.d0, (uint64)__x2.d1 << 1 ,     __s2 ,    __hi2.d0 );\
		MUL_LOHI64(__x3.d0, (uint64)__x3.d1 << 1 ,     __s3 ,    __hi3.d0 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		__hi0.d0+= (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1.d0+= (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2.d0+= (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3.d0+= (uint64)__x3.d1*(uint64)__x3.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__b0 = __b0 + __s0;	__hi0.d0 += (__b0 < __s0);\
		__b1 = __b1 + __s1;	__hi1.d0 += (__b1 < __s1);\
		__b2 = __b2 + __s2;	__hi2.d0 += (__b2 < __s2);\
		__b3 = __b3 + __s3;	__hi3.d0 += (__b3 < __s3);\
		\
		__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d1 = __b3 & 0x00000000ffffffff;\
		\
		__hi0.d1 = (__hi0.d0 >> 32);\
		__hi1.d1 = (__hi1.d0 >> 32);\
		__hi2.d1 = (__hi2.d0 >> 32);\
		__hi3.d1 = (__hi3.d0 >> 32);\
		\
		__hi0.d0 = (__hi0.d0 << 32) + (__b0 >> 32);\
		__hi1.d0 = (__hi1.d0 << 32) + (__b1 >> 32);\
		__hi2.d0 = (__hi2.d0 << 32) + (__b2 >> 32);\
		__hi3.d0 = (__hi3.d0 << 32) + (__b3 >> 32);\
    }
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI95_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
	{\
		SQR_LOHI96_q4(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3)\
	}


  #else

    /* If slow MULH64, use specialized 64x32-bit cross-product MULs and extra ALU ops
    if inputs are genuinely 96-bit. In this case it pays to use a specialized fast version
    for inputs that are 95 bits or less, which saves 3 ALU ops per input: */
    #define SQR_LOHI96_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64 __l0,__m0,__h0,__a0,__b0;\
		uint64 __l1,__m1,__h1,__a1,__b1;\
		uint64 __l2,__m2,__h2,__a2,__b2;\
		uint64 __l3,__m3,__h3,__a3,__b3;\
		uint32 __tt0 = (uint32)__x0.d1, __hl32_0,__hh32_0;\
		uint32 __tt1 = (uint32)__x1.d1, __hl32_1,__hh32_1;\
		uint32 __tt2 = (uint32)__x2.d1, __hl32_2,__hh32_2;\
		uint32 __tt3 = (uint32)__x3.d1, __hl32_3,__hh32_3;\
		\
		MUL64x32(__x0.d0,__tt0, __a0, __b0);\
		MUL64x32(__x1.d0,__tt1, __a1, __b1);\
		MUL64x32(__x2.d0,__tt2, __a2, __b2);\
		MUL64x32(__x3.d0,__tt3, __a3, __b3);\
		\
		SQR_LOHI64(__x0.d0,     __l0, __m0);\
		SQR_LOHI64(__x1.d0,     __l1, __m1);\
		SQR_LOHI64(__x2.d0,     __l2, __m2);\
		SQR_LOHI64(__x3.d0,     __l3, __m3);\
		\
		MUL_LOHI32(__tt0,__tt0,__hl32_0,__hh32_0);\
		MUL_LOHI32(__tt1,__tt1,__hl32_1,__hh32_1);\
		MUL_LOHI32(__tt2,__tt2,__hl32_2,__hh32_2);\
		MUL_LOHI32(__tt3,__tt3,__hl32_3,__hh32_3);\
		/* Double the cross term */\
		__b0 = (__b0 << 1) + (__a0 >> 63);	__a0 <<= 1;\
		__b1 = (__b1 << 1) + (__a1 >> 63);	__a1 <<= 1;\
		__b2 = (__b2 << 1) + (__a2 >> 63);	__a2 <<= 1;\
		__b3 = (__b3 << 1) + (__a3 >> 63);	__a3 <<= 1;\
		\
		__h0 = (uint64)(__hl32_0) + ((uint64)__hh32_0 << 32);\
		__h1 = (uint64)(__hl32_1) + ((uint64)__hh32_1 << 32);\
		__h2 = (uint64)(__hl32_2) + ((uint64)__hh32_2 << 32);\
		__h3 = (uint64)(__hl32_3) + ((uint64)__hh32_3 << 32);\
		\
		__m0 += __a0;\
		__m1 += __a1;\
		__m2 += __a2;\
		__m3 += __a3;\
		/* Overflow into high word is checked here. */\
		__h0 += __b0 + (__m0 < __a0);\
		__h1 += __b1 + (__m1 < __a1);\
		__h2 += __b2 + (__m2 < __a2);\
		__h3 += __b3 + (__m3 < __a3);\
		/* Now split the result between __lo and __hi: */\
		__lo0.d0 =  __l0;					__lo0.d1 = (uint32)__m0;\
		__lo1.d0 =  __l1;					__lo1.d1 = (uint32)__m1;\
		__lo2.d0 =  __l2;					__lo2.d1 = (uint32)__m2;\
		__lo3.d0 =  __l3;					__lo3.d1 = (uint32)__m3;\
		\
		__hi0.d0 = (__m0>>32) + (__h0<<32);	__hi0.d1 = (__h0>>32);\
		__hi1.d0 = (__m1>>32) + (__h1<<32);	__hi1.d1 = (__h1>>32);\
		__hi2.d0 = (__m2>>32) + (__h2<<32);	__hi2.d1 = (__h2>>32);\
		__hi3.d0 = (__m3>>32) + (__h3<<32);	__hi3.d1 = (__h3>>32);\
    }

	/* Special version for x < 2^95: */
    #define SQR_LOHI95_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__s0,__s1,__s2,__s3;\
		\
		DBG_ASSERT((__x0.d1 >> 31) == 0,"SQR_LOHI95_q4: (__x0.d1 >> 31) == 0");\
		DBG_ASSERT((__x1.d1 >> 31) == 0,"SQR_LOHI95_q4: (__x1.d1 >> 31) == 0");\
		DBG_ASSERT((__x2.d1 >> 31) == 0,"SQR_LOHI95_q4: (__x2.d1 >> 31) == 0");\
		DBG_ASSERT((__x3.d1 >> 31) == 0,"SQR_LOHI95_q4: (__x3.d1 >> 31) == 0");\
		\
		SQR_LOHI64(__x0.d0,                    __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                    __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                    __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                    __a3 ,     __b3 );\
		\
		MUL64x32(__x0.d0, __x0.d1 << 1 ,     __s0 , __hi0.d0 );\
		MUL64x32(__x1.d0, __x1.d1 << 1 ,     __s1 , __hi1.d0 );\
		MUL64x32(__x2.d0, __x2.d1 << 1 ,     __s2 , __hi2.d0 );\
		MUL64x32(__x3.d0, __x3.d1 << 1 ,     __s3 , __hi3.d0 );\
		\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		__hi0.d0+= (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1.d0+= (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2.d0+= (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3.d0+= (uint64)__x3.d1*(uint64)__x3.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__b0 = __b0 + __s0;	__hi0.d0 += (__b0 < __s0);\
		__b1 = __b1 + __s1;	__hi1.d0 += (__b1 < __s1);\
		__b2 = __b2 + __s2;	__hi2.d0 += (__b2 < __s2);\
		__b3 = __b3 + __s3;	__hi3.d0 += (__b3 < __s3);\
		\
		__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d1 = __b3 & 0x00000000ffffffff;\
		\
		__hi0.d1 = (__hi0.d0 >> 32);\
		__hi1.d1 = (__hi1.d0 >> 32);\
		__hi2.d1 = (__hi2.d0 >> 32);\
		__hi3.d1 = (__hi3.d0 >> 32);\
		\
		__hi0.d0 = (__hi0.d0 << 32) + (__b0 >> 32);\
		__hi1.d0 = (__hi1.d0 << 32) + (__b1 >> 32);\
		__hi2.d0 = (__hi2.d0 << 32) + (__b2 >> 32);\
		__hi3.d0 = (__hi3.d0 << 32) + (__b3 >> 32);\
    }
  #endif

    /* 8-operand-pipelined version: designed so low half of output may be returned
    in input argument x, if x and lo happen to point to the same addresses:
    */
  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI96_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x7.d1 >> 32) == 0");\
		\
		SQR_LOHI    (__x0.d0,                             __a0 ,     __b0 );\
		SQR_LOHI    (__x1.d0,                             __a1 ,     __b1 );\
		SQR_LOHI    (__x2.d0,                             __a2 ,     __b2 );\
		SQR_LOHI    (__x3.d0,                             __a3 ,     __b3 );\
		SQR_LOHI    (__x4.d0,                             __a4 ,     __b4 );\
		SQR_LOHI    (__x5.d0,                             __a5 ,     __b5 );\
		SQR_LOHI    (__x6.d0,                             __a6 ,     __b6 );\
		SQR_LOHI    (__x7.d0,                             __a7 ,     __b7 );\
		\
		MUL_LOHI64_ADD(__x0.d0, (uint64)__x0.d1 << 1,     __b0,     __b0 , __hi0.d0 );\
		MUL_LOHI64_ADD(__x1.d0, (uint64)__x1.d1 << 1,     __b1,     __b1 , __hi1.d0 );\
		MUL_LOHI64_ADD(__x2.d0, (uint64)__x2.d1 << 1,     __b2,     __b2 , __hi2.d0 );\
		MUL_LOHI64_ADD(__x3.d0, (uint64)__x3.d1 << 1,     __b3,     __b3 , __hi3.d0 );\
		MUL_LOHI64_ADD(__x4.d0, (uint64)__x4.d1 << 1,     __b4,     __b4 , __hi4.d0 );\
		MUL_LOHI64_ADD(__x5.d0, (uint64)__x5.d1 << 1,     __b5,     __b5 , __hi5.d0 );\
		MUL_LOHI64_ADD(__x6.d0, (uint64)__x6.d1 << 1,     __b6,     __b6 , __hi6.d0 );\
		MUL_LOHI64_ADD(__x7.d0, (uint64)__x7.d1 << 1,     __b7,     __b7 , __hi7.d0 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		__hi0.d0 += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1.d0 += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2.d0 += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3.d0 += (uint64)__x3.d1*(uint64)__x3.d1;\
		__hi4.d0 += (uint64)__x4.d1*(uint64)__x4.d1;\
		__hi5.d0 += (uint64)__x5.d1*(uint64)__x5.d1;\
		__hi6.d0 += (uint64)__x6.d1*(uint64)__x6.d1;\
		__hi7.d0 += (uint64)__x7.d1*(uint64)__x7.d1;\
		\
		__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d1 = __b3 & 0x00000000ffffffff;\
		__lo4.d1 = __b4 & 0x00000000ffffffff;\
		__lo5.d1 = __b5 & 0x00000000ffffffff;\
		__lo6.d1 = __b6 & 0x00000000ffffffff;\
		__lo7.d1 = __b7 & 0x00000000ffffffff;\
		\
		__hi0.d1 = (__hi0.d0 >> 32);\
		__hi1.d1 = (__hi1.d0 >> 32);\
		__hi2.d1 = (__hi2.d0 >> 32);\
		__hi3.d1 = (__hi3.d0 >> 32);\
		__hi4.d1 = (__hi4.d0 >> 32);\
		__hi5.d1 = (__hi5.d0 >> 32);\
		__hi6.d1 = (__hi6.d0 >> 32);\
		__hi7.d1 = (__hi7.d0 >> 32);\
		\
		__hi0.d0 = (__hi0.d0 << 32) + (__b0 >> 32);\
		__hi1.d0 = (__hi1.d0 << 32) + (__b1 >> 32);\
		__hi2.d0 = (__hi2.d0 << 32) + (__b2 >> 32);\
		__hi3.d0 = (__hi3.d0 << 32) + (__b3 >> 32);\
		__hi4.d0 = (__hi4.d0 << 32) + (__b4 >> 32);\
		__hi5.d0 = (__hi5.d0 << 32) + (__b5 >> 32);\
		__hi6.d0 = (__hi6.d0 << 32) + (__b6 >> 32);\
		__hi7.d0 = (__hi7.d0 << 32) + (__b7 >> 32);\
    }
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI95_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
	{\
		SQR_LOHI96_q8(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3\
		, __x4, __lo4, __hi4\
		, __x5, __lo5, __hi5\
		, __x6, __lo6, __hi6\
		, __x7, __lo7, __hi7)\
	}

  #elif(defined(MULH64_FAST))	/* Fast 64-bit MUL but no integer fused MADD */

    #define SQR_LOHI96_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x7.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,                    __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                    __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                    __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                    __a3 ,     __b3 );\
		SQR_LOHI64(__x4.d0,                    __a4 ,     __b4 );\
		SQR_LOHI64(__x5.d0,                    __a5 ,     __b5 );\
		SQR_LOHI64(__x6.d0,                    __a6 ,     __b6 );\
		SQR_LOHI64(__x7.d0,                    __a7 ,     __b7 );\
		\
		MUL_LOHI64(__x0.d0, (uint64)__x0.d1 << 1 ,     __s0 , __hi0.d0 );\
		MUL_LOHI64(__x1.d0, (uint64)__x1.d1 << 1 ,     __s1 , __hi1.d0 );\
		MUL_LOHI64(__x2.d0, (uint64)__x2.d1 << 1 ,     __s2 , __hi2.d0 );\
		MUL_LOHI64(__x3.d0, (uint64)__x3.d1 << 1 ,     __s3 , __hi3.d0 );\
		MUL_LOHI64(__x4.d0, (uint64)__x4.d1 << 1 ,     __s4 , __hi4.d0 );\
		MUL_LOHI64(__x5.d0, (uint64)__x5.d1 << 1 ,     __s5 , __hi5.d0 );\
		MUL_LOHI64(__x6.d0, (uint64)__x6.d1 << 1 ,     __s6 , __hi6.d0 );\
		MUL_LOHI64(__x7.d0, (uint64)__x7.d1 << 1 ,     __s7 , __hi7.d0 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		__hi0.d0+= (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1.d0+= (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2.d0+= (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3.d0+= (uint64)__x3.d1*(uint64)__x3.d1;\
		__hi4.d0+= (uint64)__x4.d1*(uint64)__x4.d1;\
		__hi5.d0+= (uint64)__x5.d1*(uint64)__x5.d1;\
		__hi6.d0+= (uint64)__x6.d1*(uint64)__x6.d1;\
		__hi7.d0+= (uint64)__x7.d1*(uint64)__x7.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__b0 = __b0 + __s0;	__hi0.d0 += (__b0 < __s0);\
		__b1 = __b1 + __s1;	__hi1.d0 += (__b1 < __s1);\
		__b2 = __b2 + __s2;	__hi2.d0 += (__b2 < __s2);\
		__b3 = __b3 + __s3;	__hi3.d0 += (__b3 < __s3);\
		__b4 = __b4 + __s4;	__hi4.d0 += (__b4 < __s4);\
		__b5 = __b5 + __s5;	__hi5.d0 += (__b5 < __s5);\
		__b6 = __b6 + __s6;	__hi6.d0 += (__b6 < __s6);\
		__b7 = __b7 + __s7;	__hi7.d0 += (__b7 < __s7);\
		\
		__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d1 = __b3 & 0x00000000ffffffff;\
		__lo4.d1 = __b4 & 0x00000000ffffffff;\
		__lo5.d1 = __b5 & 0x00000000ffffffff;\
		__lo6.d1 = __b6 & 0x00000000ffffffff;\
		__lo7.d1 = __b7 & 0x00000000ffffffff;\
		\
		__hi0.d1 = (__hi0.d0 >> 32);\
		__hi1.d1 = (__hi1.d0 >> 32);\
		__hi2.d1 = (__hi2.d0 >> 32);\
		__hi3.d1 = (__hi3.d0 >> 32);\
		__hi4.d1 = (__hi4.d0 >> 32);\
		__hi5.d1 = (__hi5.d0 >> 32);\
		__hi6.d1 = (__hi6.d0 >> 32);\
		__hi7.d1 = (__hi7.d0 >> 32);\
		\
		__hi0.d0 = (__hi0.d0 << 32) + (__b0 >> 32);\
		__hi1.d0 = (__hi1.d0 << 32) + (__b1 >> 32);\
		__hi2.d0 = (__hi2.d0 << 32) + (__b2 >> 32);\
		__hi3.d0 = (__hi3.d0 << 32) + (__b3 >> 32);\
		__hi4.d0 = (__hi4.d0 << 32) + (__b4 >> 32);\
		__hi5.d0 = (__hi5.d0 << 32) + (__b5 >> 32);\
		__hi6.d0 = (__hi6.d0 << 32) + (__b6 >> 32);\
		__hi7.d0 = (__hi7.d0 << 32) + (__b7 >> 32);\
    }
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI95_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
	{\
		SQR_LOHI96_q8(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3\
		, __x4, __lo4, __hi4\
		, __x5, __lo5, __hi5\
		, __x6, __lo6, __hi6\
		, __x7, __lo7, __hi7)\
	}

  #else

    /* If slow MULH64, use specialized 64x32-bit cross-product MULs and extra ALU ops
    if inputs are genuinely 96-bit. In this case it pays to use a specialized fast version
    for inputs that are 95 bits or less, which saves 3 ALU ops per input: */
    #define SQR_LOHI96_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 32) == 0,"SQR_LOHI96_q8: (__x7.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,              __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,              __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,              __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,              __a3 ,     __b3 );\
		SQR_LOHI64(__x4.d0,              __a4 ,     __b4 );\
		SQR_LOHI64(__x5.d0,              __a5 ,     __b5 );\
		SQR_LOHI64(__x6.d0,              __a6 ,     __b6 );\
		SQR_LOHI64(__x7.d0,              __a7 ,     __b7 );\
		\
		MUL64x32(__x0.d0, __x0.d1,     __s0 , __hi0.d0 );\
		MUL64x32(__x1.d0, __x1.d1,     __s1 , __hi1.d0 );\
		MUL64x32(__x2.d0, __x2.d1,     __s2 , __hi2.d0 );\
		MUL64x32(__x3.d0, __x3.d1,     __s3 , __hi3.d0 );\
		MUL64x32(__x4.d0, __x4.d1,     __s4 , __hi4.d0 );\
		MUL64x32(__x5.d0, __x5.d1,     __s5 , __hi5.d0 );\
		MUL64x32(__x6.d0, __x6.d1,     __s6 , __hi6.d0 );\
		MUL64x32(__x7.d0, __x7.d1,     __s7 , __hi7.d0 );\
		\
		/* Double the cross product term: */\
		__hi0.d0 = (__hi0.d0 << 1) + (__s0 >> 63);\
		__hi1.d0 = (__hi1.d0 << 1) + (__s1 >> 63);\
		__hi2.d0 = (__hi2.d0 << 1) + (__s2 >> 63);\
		__hi3.d0 = (__hi3.d0 << 1) + (__s3 >> 63);\
		__hi4.d0 = (__hi4.d0 << 1) + (__s4 >> 63);\
		__hi5.d0 = (__hi5.d0 << 1) + (__s5 >> 63);\
		__hi6.d0 = (__hi6.d0 << 1) + (__s6 >> 63);\
		__hi7.d0 = (__hi7.d0 << 1) + (__s7 >> 63);\
		\
		__s0 <<= 1;\
		__s1 <<= 1;\
		__s2 <<= 1;\
		__s3 <<= 1;\
		__s4 <<= 1;\
		__s5 <<= 1;\
		__s6 <<= 1;\
		__s7 <<= 1;\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		__hi0.d0+= (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1.d0+= (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2.d0+= (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3.d0+= (uint64)__x3.d1*(uint64)__x3.d1;\
		__hi4.d0+= (uint64)__x4.d1*(uint64)__x4.d1;\
		__hi5.d0+= (uint64)__x5.d1*(uint64)__x5.d1;\
		__hi6.d0+= (uint64)__x6.d1*(uint64)__x6.d1;\
		__hi7.d0+= (uint64)__x7.d1*(uint64)__x7.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__b0 = __b0 + __s0;	__hi0.d0 += (__b0 < __s0);\
		__b1 = __b1 + __s1;	__hi1.d0 += (__b1 < __s1);\
		__b2 = __b2 + __s2;	__hi2.d0 += (__b2 < __s2);\
		__b3 = __b3 + __s3;	__hi3.d0 += (__b3 < __s3);\
		__b4 = __b4 + __s4;	__hi4.d0 += (__b4 < __s4);\
		__b5 = __b5 + __s5;	__hi5.d0 += (__b5 < __s5);\
		__b6 = __b6 + __s6;	__hi6.d0 += (__b6 < __s6);\
		__b7 = __b7 + __s7;	__hi7.d0 += (__b7 < __s7);\
		\
		__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d1 = __b3 & 0x00000000ffffffff;\
		__lo4.d1 = __b4 & 0x00000000ffffffff;\
		__lo5.d1 = __b5 & 0x00000000ffffffff;\
		__lo6.d1 = __b6 & 0x00000000ffffffff;\
		__lo7.d1 = __b7 & 0x00000000ffffffff;\
		\
		__hi0.d1 = (__hi0.d0 >> 32);\
		__hi1.d1 = (__hi1.d0 >> 32);\
		__hi2.d1 = (__hi2.d0 >> 32);\
		__hi3.d1 = (__hi3.d0 >> 32);\
		__hi4.d1 = (__hi4.d0 >> 32);\
		__hi5.d1 = (__hi5.d0 >> 32);\
		__hi6.d1 = (__hi6.d0 >> 32);\
		__hi7.d1 = (__hi7.d0 >> 32);\
		\
		__hi0.d0 = (__hi0.d0 << 32) + (__b0 >> 32);\
		__hi1.d0 = (__hi1.d0 << 32) + (__b1 >> 32);\
		__hi2.d0 = (__hi2.d0 << 32) + (__b2 >> 32);\
		__hi3.d0 = (__hi3.d0 << 32) + (__b3 >> 32);\
		__hi4.d0 = (__hi4.d0 << 32) + (__b4 >> 32);\
		__hi5.d0 = (__hi5.d0 << 32) + (__b5 >> 32);\
		__hi6.d0 = (__hi6.d0 << 32) + (__b6 >> 32);\
		__hi7.d0 = (__hi7.d0 << 32) + (__b7 >> 32);\
    }

	/* Special version for x < 2^95: */
    #define SQR_LOHI95_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7;\
		\
		DBG_ASSERT((__x0.d1 >> 31) == 0,"SQR_LOHI95_q8: (__x0.d1 >> 31) == 0");\
		DBG_ASSERT((__x1.d1 >> 31) == 0,"SQR_LOHI95_q8: (__x1.d1 >> 31) == 0");\
		DBG_ASSERT((__x2.d1 >> 31) == 0,"SQR_LOHI95_q8: (__x2.d1 >> 31) == 0");\
		DBG_ASSERT((__x3.d1 >> 31) == 0,"SQR_LOHI95_q8: (__x3.d1 >> 31) == 0");\
		DBG_ASSERT((__x4.d1 >> 31) == 0,"SQR_LOHI95_q8: (__x4.d1 >> 31) == 0");\
		DBG_ASSERT((__x5.d1 >> 31) == 0,"SQR_LOHI95_q8: (__x5.d1 >> 31) == 0");\
		DBG_ASSERT((__x6.d1 >> 31) == 0,"SQR_LOHI95_q8: (__x6.d1 >> 31) == 0");\
		DBG_ASSERT((__x7.d1 >> 31) == 0,"SQR_LOHI95_q8: (__x7.d1 >> 31) == 0");\
		\
		SQR_LOHI64(__x0.d0,                    __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                    __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                    __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                    __a3 ,     __b3 );\
		SQR_LOHI64(__x4.d0,                    __a4 ,     __b4 );\
		SQR_LOHI64(__x5.d0,                    __a5 ,     __b5 );\
		SQR_LOHI64(__x6.d0,                    __a6 ,     __b6 );\
		SQR_LOHI64(__x7.d0,                    __a7 ,     __b7 );\
		\
		MUL64x32(__x0.d0, __x0.d1 << 1 ,     __s0 , __hi0.d0 );\
		MUL64x32(__x1.d0, __x1.d1 << 1 ,     __s1 , __hi1.d0 );\
		MUL64x32(__x2.d0, __x2.d1 << 1 ,     __s2 , __hi2.d0 );\
		MUL64x32(__x3.d0, __x3.d1 << 1 ,     __s3 , __hi3.d0 );\
		MUL64x32(__x4.d0, __x4.d1 << 1 ,     __s4 , __hi4.d0 );\
		MUL64x32(__x5.d0, __x5.d1 << 1 ,     __s5 , __hi5.d0 );\
		MUL64x32(__x6.d0, __x6.d1 << 1 ,     __s6 , __hi6.d0 );\
		MUL64x32(__x7.d0, __x7.d1 << 1 ,     __s7 , __hi7.d0 );\
		\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		__hi0.d0+= (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1.d0+= (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2.d0+= (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3.d0+= (uint64)__x3.d1*(uint64)__x3.d1;\
		__hi4.d0+= (uint64)__x4.d1*(uint64)__x4.d1;\
		__hi5.d0+= (uint64)__x5.d1*(uint64)__x5.d1;\
		__hi6.d0+= (uint64)__x6.d1*(uint64)__x6.d1;\
		__hi7.d0+= (uint64)__x7.d1*(uint64)__x7.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__b0 = __b0 + __s0;	__hi0.d0 += (__b0 < __s0);\
		__b1 = __b1 + __s1;	__hi1.d0 += (__b1 < __s1);\
		__b2 = __b2 + __s2;	__hi2.d0 += (__b2 < __s2);\
		__b3 = __b3 + __s3;	__hi3.d0 += (__b3 < __s3);\
		__b4 = __b4 + __s4;	__hi4.d0 += (__b4 < __s4);\
		__b5 = __b5 + __s5;	__hi5.d0 += (__b5 < __s5);\
		__b6 = __b6 + __s6;	__hi6.d0 += (__b6 < __s6);\
		__b7 = __b7 + __s7;	__hi7.d0 += (__b7 < __s7);\
		\
		__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d1 = __b3 & 0x00000000ffffffff;\
		__lo4.d1 = __b4 & 0x00000000ffffffff;\
		__lo5.d1 = __b5 & 0x00000000ffffffff;\
		__lo6.d1 = __b6 & 0x00000000ffffffff;\
		__lo7.d1 = __b7 & 0x00000000ffffffff;\
		\
		__hi0.d1 = (__hi0.d0 >> 32);\
		__hi1.d1 = (__hi1.d0 >> 32);\
		__hi2.d1 = (__hi2.d0 >> 32);\
		__hi3.d1 = (__hi3.d0 >> 32);\
		__hi4.d1 = (__hi4.d0 >> 32);\
		__hi5.d1 = (__hi5.d0 >> 32);\
		__hi6.d1 = (__hi6.d0 >> 32);\
		__hi7.d1 = (__hi7.d0 >> 32);\
		\
		__hi0.d0 = (__hi0.d0 << 32) + (__b0 >> 32);\
		__hi1.d0 = (__hi1.d0 << 32) + (__b1 >> 32);\
		__hi2.d0 = (__hi2.d0 << 32) + (__b2 >> 32);\
		__hi3.d0 = (__hi3.d0 << 32) + (__b3 >> 32);\
		__hi4.d0 = (__hi4.d0 << 32) + (__b4 >> 32);\
		__hi5.d0 = (__hi5.d0 << 32) + (__b5 >> 32);\
		__hi6.d0 = (__hi6.d0 << 32) + (__b6 >> 32);\
		__hi7.d0 = (__hi7.d0 << 32) + (__b7 >> 32);\
    }
  #endif

#endif

/* 192-bit product of uint96 *x and *y. Result is returned in a pair of uint96s.

On Alpha, this needs a total of 7 MUL, 12 ALU op.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MUL_LOHI96(__x, __y, __lo, __hi)\
    {\
		uint64 __l,__m,__h,__a,__b,__c,__d,__s,__t;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MUL_LOHI96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MUL_LOHI96: (__y.d1 >> 32) == 0");\
		__s   = (uint64)(__x.d1);\
		__t   = (uint64)(__y.d1);\
		__h   = __s*__t;\
		MUL_LOHI64(__x.d0, __y.d0,&__l,&__m);\
		MUL64x32(__x.d0, __t   ,&__a,&__b);\
		MUL64x32(__y.d0, __s   ,&__c,&__d);\
		__m += __a;\
		__h += __b + (__m < __a); /* Overflow into high word is checked here. */\
		__m += __c;\
		__h += __d + (__m < __c); /* Overflow into high word is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __l;						__lo.d1 =  __m;\
		__hi.d0 = (__m>>32) + (__h<<32);	__hi.d1 = (__h>>32);\
    }

#else

    #define MUL_LOHI96(__x, __y, __lo, __hi)\
    {\
		uint64 __l,__m,__h,__a,__b,__c,__d,__s,__t;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MUL_LOHI96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MUL_LOHI96: (__y.d1 >> 32) == 0");\
		__s   = (uint64)(__x.d1);\
		__t   = (uint64)(__y.d1);\
		__h   = __s*__t;\
		MUL_LOHI64(__x.d0,__y.d0, __l, __m);\
		MUL64x32(__x.d0,__t   , __a, __b);\
		MUL64x32(__y.d0,__s   , __c, __d);\
		__m += __a;\
		__h += __b + (__m < __a); /* Overflow into high word is checked here. */\
		__m += __c;\
		__h += __d + (__m < __c); /* Overflow into high word is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __l;						__lo.d1 =  __m;\
		__hi.d0 = (__m>>32) + (__h<<32);	__hi.d1 = (__h>>32);\
    }

#endif

/* 192-bit product of uint96 *x and *y. Result is returned in a uint192.

On Alpha, this needs a total of 7 MUL, 12 ALU op.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MUL_LOHI96_PROD192(__x, __y, __prod192)\
    {\
		uint64 __l,__m,__h,__a,__b,__c,__d,__s,__t;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MUL_LOHI96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MUL_LOHI96: (__y.d1 >> 32) == 0");\
		__s   = (uint64)(__x.d1);\
		__t   = (uint64)(__y.d1);\
		__h   = __s*__t;\
		MUL_LOHI64(__x.d0, __y.d0,&__l,&__m);\
		MUL64x32(__x.d0, __t   ,&__a,&__b);\
		MUL64x32(__y.d0, __s   ,&__c,&__d);\
		__m += __a;\
		__h += __b + (__m < __a); /* Overflow into high word is checked here. */\
		__m += __c;\
		__h += __d + (__m < __c); /* Overflow into high word is checked here. */\
		/* Now store the result in __prod192: */\
		__prod192.d0 =  __l;	__prod192.d1 =  __m;	__prod192.d2 = __h;\
    }

#else

    #define MUL_LOHI96_PROD192(__x, __y, __prod192)\
    {\
		uint64 __l,__m,__h,__a,__b,__c,__d,__s,__t;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MUL_LOHI96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MUL_LOHI96: (__y.d1 >> 32) == 0");\
		__s   = (uint64)(__x.d1);\
		__t   = (uint64)(__y.d1);\
		__h   = __s*__t;\
		MUL_LOHI64(__x.d0,__y.d0, __l, __m);\
		MUL64x32(__x.d0,__t   , __a, __b);\
		MUL64x32(__y.d0,__s   , __c, __d);\
		__m += __a;\
		__h += __b + (__m < __a); /* Overflow into high word is checked here. */\
		__m += __c;\
		__h += __d + (__m < __c); /* Overflow into high word is checked here. */\
		/* Now store the result in __prod192: */\
		__prod192.d0 =  __l;	__prod192.d1 =  __m;	__prod192.d2 = __h;\
    }

#endif

/* 96-bit product modulo 2^96 of uint96 *x and *y. Result is returned in a uint96.
   Designed so output may be returned in input argument x,
   if x and lo happen to point to the same addresses.

   On Alpha, this needs a total of 4 MUL (2 MULL32, 1 MULL64, 1 MULH64), 3 ALU.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULL96(__x, __y, __lo)\
    {\
		uint64 __l,__m;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MULL96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MULL96: (__y.d1 >> 32) == 0");\
		MUL_LOHI64(__x.d0,__y.d0,&__l,&__m);\
		__m += __MULL32(__x.d1,__y.d0) + __MULL32(__y.d1,__x.d0);	/* Only need the bottom 32 bits of each product here */\
		__lo.d0 =  __l;	__lo.d1 = __m & 0x00000000ffffffff;\
    }

    /* 4-operand-pipelined version: */
    #define MULL96_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
		MULL96(__x0, __y0, __lo0);\
		MULL96(__x1, __y1, __lo1);\
		MULL96(__x2, __y2, __lo2);\
		MULL96(__x3, __y3, __lo3);\
	}

    /* 8-operand-pipelined version: */
    #define MULL96_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
		MULL96(__x0, __y0, __lo0);\
		MULL96(__x1, __y1, __lo1);\
		MULL96(__x2, __y2, __lo2);\
		MULL96(__x3, __y3, __lo3);\
		MULL96(__x4, __y4, __lo4);\
		MULL96(__x5, __y5, __lo5);\
		MULL96(__x6, __y6, __lo6);\
		MULL96(__x7, __y7, __lo7);\
	}

#else

    #define MULL96(__x, __y, __lo)\
    {\
		uint64 __l,__m;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MULL96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MULL96: (__y.d1 >> 32) == 0");\
		MUL_LOHI64(__x.d0,__y.d0, __l, __m);\
		__m += __MULL32(__x.d1,__y.d0) + __MULL32(__y.d1,__x.d0);	/* Only need the bottom 32 bits of each product here */\
		__lo.d0 =  __l;	__lo.d1 = __m & 0x00000000ffffffff;\
    }

    /* 4-operand-pipelined version: designed so low half of output may be returned
    in either of the input arguments, if it and lo happen to point to the same address:
    */
    #define MULL96_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"MULL96_q4: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"MULL96_q4: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"MULL96_q4: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"MULL96_q4: (__x3.d1 >> 32) == 0");\
		\
		DBG_ASSERT((__y0.d1 >> 32) == 0,"MULL96_q4: (__y0.d1 >> 32) == 0");\
		DBG_ASSERT((__y1.d1 >> 32) == 0,"MULL96_q4: (__y1.d1 >> 32) == 0");\
		DBG_ASSERT((__y2.d1 >> 32) == 0,"MULL96_q4: (__y2.d1 >> 32) == 0");\
		DBG_ASSERT((__y3.d1 >> 32) == 0,"MULL96_q4: (__y3.d1 >> 32) == 0");\
		\
		MUL_LOHI64(__x0.d0,__y0.d0, __a0, __b0);\
		MUL_LOHI64(__x1.d0,__y1.d0, __a1, __b1);\
		MUL_LOHI64(__x2.d0,__y2.d0, __a2, __b2);\
		MUL_LOHI64(__x3.d0,__y3.d0, __a3, __b3);\
		\
		__b0 += __MULL32(__x0.d1,__y0.d0) + __MULL32(__y0.d1,__x0.d0);\
		__b1 += __MULL32(__x1.d1,__y1.d0) + __MULL32(__y1.d1,__x1.d0);\
		__b2 += __MULL32(__x2.d1,__y2.d0) + __MULL32(__y2.d1,__x2.d0);\
		__b3 += __MULL32(__x3.d1,__y3.d0) + __MULL32(__y3.d1,__x3.d0);\
		\
		__lo0.d0 =  __a0;	__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d0 =  __a1;	__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d0 =  __a2;	__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d0 =  __a3;	__lo3.d1 = __b3 & 0x00000000ffffffff;\
	}

    /* 8-operand-pipelined version: designed so low half of output may be returned
    in either of the input arguments, if it and lo happen to point to the same address:
    */
    #define MULL96_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"MULL96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"MULL96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"MULL96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"MULL96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 32) == 0,"MULL96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 32) == 0,"MULL96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 32) == 0,"MULL96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 32) == 0,"MULL96_q8: (__x7.d1 >> 32) == 0");\
		\
		DBG_ASSERT((__y0.d1 >> 32) == 0,"MULL96_q8: (__y0.d1 >> 32) == 0");\
		DBG_ASSERT((__y1.d1 >> 32) == 0,"MULL96_q8: (__y1.d1 >> 32) == 0");\
		DBG_ASSERT((__y2.d1 >> 32) == 0,"MULL96_q8: (__y2.d1 >> 32) == 0");\
		DBG_ASSERT((__y3.d1 >> 32) == 0,"MULL96_q8: (__y3.d1 >> 32) == 0");\
		DBG_ASSERT((__y4.d1 >> 32) == 0,"MULL96_q8: (__y4.d1 >> 32) == 0");\
		DBG_ASSERT((__y5.d1 >> 32) == 0,"MULL96_q8: (__y5.d1 >> 32) == 0");\
		DBG_ASSERT((__y6.d1 >> 32) == 0,"MULL96_q8: (__y6.d1 >> 32) == 0");\
		DBG_ASSERT((__y7.d1 >> 32) == 0,"MULL96_q8: (__y7.d1 >> 32) == 0");\
		\
		MUL_LOHI64(__x0.d0,__y0.d0, __a0, __b0);\
		MUL_LOHI64(__x1.d0,__y1.d0, __a1, __b1);\
		MUL_LOHI64(__x2.d0,__y2.d0, __a2, __b2);\
		MUL_LOHI64(__x3.d0,__y3.d0, __a3, __b3);\
		MUL_LOHI64(__x4.d0,__y4.d0, __a4, __b4);\
		MUL_LOHI64(__x5.d0,__y5.d0, __a5, __b5);\
		MUL_LOHI64(__x6.d0,__y6.d0, __a6, __b6);\
		MUL_LOHI64(__x7.d0,__y7.d0, __a7, __b7);\
		\
		__b0 += __MULL32(__x0.d1,__y0.d0) + __MULL32(__y0.d1,__x0.d0);\
		__b1 += __MULL32(__x1.d1,__y1.d0) + __MULL32(__y1.d1,__x1.d0);\
		__b2 += __MULL32(__x2.d1,__y2.d0) + __MULL32(__y2.d1,__x2.d0);\
		__b3 += __MULL32(__x3.d1,__y3.d0) + __MULL32(__y3.d1,__x3.d0);\
		__b4 += __MULL32(__x4.d1,__y4.d0) + __MULL32(__y4.d1,__x4.d0);\
		__b5 += __MULL32(__x5.d1,__y5.d0) + __MULL32(__y5.d1,__x5.d0);\
		__b6 += __MULL32(__x6.d1,__y6.d0) + __MULL32(__y6.d1,__x6.d0);\
		__b7 += __MULL32(__x7.d1,__y7.d0) + __MULL32(__y7.d1,__x7.d0);\
		\
		__lo0.d0 =  __a0;	__lo0.d1 = __b0 & 0x00000000ffffffff;\
		__lo1.d0 =  __a1;	__lo1.d1 = __b1 & 0x00000000ffffffff;\
		__lo2.d0 =  __a2;	__lo2.d1 = __b2 & 0x00000000ffffffff;\
		__lo3.d0 =  __a3;	__lo3.d1 = __b3 & 0x00000000ffffffff;\
		__lo4.d0 =  __a4;	__lo4.d1 = __b4 & 0x00000000ffffffff;\
		__lo5.d0 =  __a5;	__lo5.d1 = __b5 & 0x00000000ffffffff;\
		__lo6.d0 =  __a6;	__lo6.d1 = __b6 & 0x00000000ffffffff;\
		__lo7.d0 =  __a7;	__lo7.d1 = __b7 & 0x00000000ffffffff;\
	}

#endif

/* Upper half of 192-bit product of uint96 x and y. Result is returned in a uint96.

On Alpha, this needs a total of 6 MUL instructions and 12 ALU ops.

NOTE: if for operands < 2^96 we instead use a MULL128 (mul modulo 2^128), can get the
corresponding MULH via a single Alpha UMULH instruction, applied to the high 64 bits
of each input. For random inputs the result will be incorrect in ~1/2^64 of cases due to
neglect of the lower bits, but that seems well below the likely level of hardware error.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULH96(__x, __y, __hi)\
    {\
		uint64 __m,__h,__a,__b,__c,__d,__s,__t;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MULH96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MULH96: (__y.d1 >> 32) == 0");\
		__s   = (uint64)(__x.d1);\
		__t   = (uint64)(__y.d1);\
		__h   = __s*__t;\
		MULH64(__x.d0,__y.d0, __m     );\
		MUL64x32(__x.d0, __t     ,&__a,&__b);\
		MUL64x32(__y.d0, __s     ,&__c,&__d);\
		__m += __a;\
		__h += __b + (__m < __a); /* Overflow into high word is checked here. */\
		__m += __c;\
		__h += __d + (__m < __c); /* Overflow into high word is checked here. */\
		/* Now put upper half of the result into __hi: */\
		__hi.d0 = (__m>>32) + (__h<<32);	__hi.d1 = (uint32)(__h>>32);\
    }

    /* 4-operand-pipelined version: */
    #define MULH96_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		MULH96(__x0, __lo0, __hi0);\
		MULH96(__x1, __lo1, __hi1);\
		MULH96(__x2, __lo2, __hi2);\
		MULH96(__x3, __lo3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define MULH96_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		MULH96(__x0, __lo0, __hi0);\
		MULH96(__x1, __lo1, __hi1);\
		MULH96(__x2, __lo2, __hi2);\
		MULH96(__x3, __lo3, __hi3);\
		MULH96(__x4, __lo4, __hi4);\
		MULH96(__x5, __lo5, __hi5);\
		MULH96(__x6, __lo6, __hi6);\
		MULH96(__x7, __lo7, __hi7);\
	}

#else

/*
For operands < 2^80 using a MULL96, can get the corresponding MULH (< 2^64)
via fewer MUL instructions by only considering the high 64 bits of each input:

	x*y = (xlo + 2^16*xhi) * (ylo + 2^16*yhi)		(xlo, ylo < 2^16;  xhi, yhi < 2^64)
		= (xlo*ylo) + 2^16*(xhi*ylo + yhi*xlo) + 2^32*(xhi*yhi) .
	The intermediate product terms xhi*ylo and yhi*xlo are each < 2^80,
	so the upper 64 bits of each overlap the lower half of the 128-bit
	product (xhi*yhi). We can use a pair of MULH64s applied to (xhi, ylo<<48)
	and (yhi, xlo<<48) to exactly calculate these 64-bit carryins, a MULL64
	to get the lower half of (xhi*yhi), and 14 ALU ops to calculate the
	carry into the upper half of (xhi*yhi).

For random inputs the result will be incorrect in some tiny fraction of cases due to
neglect of the (xlo*ylo) term, but that again will be < the likely level of hardware error.

On MUL-friendly architectures this will likely be no better (just 2 fewer MULs,
offset by 2 more ALU ops), but e.g. on AltiVec we can consider using SIMD routines
to get the 16x64==>80-bit intermediate products.
*/
	/* 4 MUL, 14 ALU: */
    #define MULH96_80(__x, __y, __hi)\
    {\
		uint64 __a,__b,__xlo,__ylo,__xhi,__yhi,__lo;\
		\
		DBG_ASSERT((__x.d1 >> 16) == 0,"MULH96_80 : (__x.d1 >> 16) == 0");\
		DBG_ASSERT((__y.d1 >> 16) == 0,"MULH96_80 : (__y.d1 >> 16) == 0");\
		__xhi =(__x.d1 << 48) + (__x.d0 >> 16);\
		__yhi =(__y.d1 << 48) + (__y.d0 >> 16);\
		__xlo = __x.d0 << 48;	/* xlo << 48 */\
		__ylo = __y.d0 << 48;	/* ylo << 48 */\
		__lo  = __xhi*__yhi;	/* Start this MUL first because we need it several cycles prior to the upper half of same */\
		MULH64(__xhi, __ylo, __a);	/* MULH64(xhi, ylo<<48) */\
		MULH64(__yhi, __xlo, __b);	/* MULH64(yhi, xlo<<48) */\
		MULH64(__xhi, __yhi, __hi);\
		__lo += __a;\
		__hi += (__lo < __a); /* Overflow into high word is checked here. */\
		__lo += __b;\
		__hi += (__lo < __b); /* Overflow into high word is checked here. */\
    }

	/* 6 MUL, 12 ALU: */
    #define MULH96(__x, __y, __hi)\
    {\
		uint64 __m,__h,__aa,__bb,__cc,__dd,__s,__t;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MULH96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MULH96: (__y.d1 >> 32) == 0");\
		__s   = (uint64)(__x.d1);\
		__t   = (uint64)(__y.d1);\
		__h   = __s*__t;\
		MULH64(__x.d0,__y.d0, __m     );\
		MUL64x32(__x.d0,__t   , __aa, __bb);\
		MUL64x32(__y.d0,__s   , __cc, __dd);\
		__m += __aa;\
		__h += __bb + (__m < __aa); /* Overflow into high word is checked here. */\
		__m += __cc;\
		__h += __dd + (__m < __cc); /* Overflow into high word is checked here. */\
		/* Now put upper half of the result into __hi: */\
		__hi.d0 = (__m>>32) + (__h<<32);	__hi.d1 = (uint32)(__h>>32);\
    }

    /* 4-operand-pipelined version: designed so low half of output may be returned
    in either of the input arguments, if it and lo happen to point to the same address:
    */
    #define MULH96_q4(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3)\
    {\
    	uint32 __l32_0, __h32_0;\
    	uint32 __l32_1, __h32_1;\
    	uint32 __l32_2, __h32_2;\
    	uint32 __l32_3, __h32_3;\
		uint64 __a0,__b0,__c0,__d0,__m0,__h0;\
		uint64 __a1,__b1,__c1,__d1,__m1,__h1;\
		uint64 __a2,__b2,__c2,__d2,__m2,__h2;\
		uint64 __a3,__b3,__c3,__d3,__m3,__h3;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"MULH96_q4: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"MULH96_q4: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"MULH96_q4: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"MULH96_q4: (__x3.d1 >> 32) == 0");\
		\
		DBG_ASSERT((__y0.d1 >> 32) == 0,"MULH96_q4: (__y0.d1 >> 32) == 0");\
		DBG_ASSERT((__y1.d1 >> 32) == 0,"MULH96_q4: (__y1.d1 >> 32) == 0");\
		DBG_ASSERT((__y2.d1 >> 32) == 0,"MULH96_q4: (__y2.d1 >> 32) == 0");\
		DBG_ASSERT((__y3.d1 >> 32) == 0,"MULH96_q4: (__y3.d1 >> 32) == 0");\
		\
		MUL_LOHI32(__x0.d1, __y0.d1, __l32_0, __h32_0);\
		MUL_LOHI32(__x1.d1, __y1.d1, __l32_1, __h32_1);\
		MUL_LOHI32(__x2.d1, __y2.d1, __l32_2, __h32_2);\
		MUL_LOHI32(__x3.d1, __y3.d1, __l32_3, __h32_3);\
		\
		__h0 = __l32_0 + ((uint64)__h32_0 << 32);\
		__h1 = __l32_1 + ((uint64)__h32_1 << 32);\
		__h2 = __l32_2 + ((uint64)__h32_2 << 32);\
		__h3 = __l32_3 + ((uint64)__h32_3 << 32);\
		\
		MULH64(__x0.d0,__y0.d0,       __m0);\
		MULH64(__x1.d0,__y1.d0,       __m1);\
		MULH64(__x2.d0,__y2.d0,       __m2);\
		MULH64(__x3.d0,__y3.d0,       __m3);\
		\
		MUL64x32(__x0.d0,__y0.d1, __a0, __b0);\
		MUL64x32(__x1.d0,__y1.d1, __a1, __b1);\
		MUL64x32(__x2.d0,__y2.d1, __a2, __b2);\
		MUL64x32(__x3.d0,__y3.d1, __a3, __b3);\
		\
		MUL64x32(__y0.d0,__x0.d1, __c0, __d0);\
		MUL64x32(__y1.d0,__x1.d1, __c1, __d1);\
		MUL64x32(__y2.d0,__x2.d1, __c2, __d2);\
		MUL64x32(__y3.d0,__x3.d1, __c3, __d3);\
		\
		__m0 += __a0;\
		__m1 += __a1;\
		__m2 += __a2;\
		__m3 += __a3;\
		\
		__h0 += __b0 + (__m0 < __a0);\
		__h1 += __b1 + (__m1 < __a1);\
		__h2 += __b2 + (__m2 < __a2);\
		__h3 += __b3 + (__m3 < __a3);\
		\
		__m0 += __c0;\
		__m1 += __c1;\
		__m2 += __c2;\
		__m3 += __c3;\
		\
		__h0 += __d0 + (__m0 < __c0);\
		__h1 += __d1 + (__m1 < __c1);\
		__h2 += __d2 + (__m2 < __c2);\
		__h3 += __d3 + (__m3 < __c3);\
		\
		/* Now put upper half of the result into __hi: */\
		__hi0.d0 = (__m0>>32) + (__h0<<32);	__hi0.d1 = (uint32)(__h0>>32);\
		__hi1.d0 = (__m1>>32) + (__h1<<32);	__hi1.d1 = (uint32)(__h1>>32);\
		__hi2.d0 = (__m2>>32) + (__h2<<32);	__hi2.d1 = (uint32)(__h2>>32);\
		__hi3.d0 = (__m3>>32) + (__h3<<32);	__hi3.d1 = (uint32)(__h3>>32);\
	}

    /* 8-operand-pipelined version: designed so low half of output may be returned
    in either of the input arguments, if it and lo happen to point to the same address:
    */
    #define MULH96_q8(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3\
    , __x4, __y4, __hi4\
    , __x5, __y5, __hi5\
    , __x6, __y6, __hi6\
    , __x7, __y7, __hi7)\
    {\
		uint64 __a0,__b0,__c0,__d0,__m0,__h0;\
		uint64 __a1,__b1,__c1,__d1,__m1,__h1;\
		uint64 __a2,__b2,__c2,__d2,__m2,__h2;\
		uint64 __a3,__b3,__c3,__d3,__m3,__h3;\
		uint64 __a4,__b4,__c4,__d4,__m4,__h4;\
		uint64 __a5,__b5,__c5,__d5,__m5,__h5;\
		uint64 __a6,__b6,__c6,__d6,__m6,__h6;\
		uint64 __a7,__b7,__c7,__d7,__m7,__h7;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"MULH96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"MULH96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"MULH96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"MULH96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 32) == 0,"MULH96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 32) == 0,"MULH96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 32) == 0,"MULH96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 32) == 0,"MULH96_q8: (__x7.d1 >> 32) == 0");\
		\
		DBG_ASSERT((__y0.d1 >> 32) == 0,"MULH96_q8: (__y0.d1 >> 32) == 0");\
		DBG_ASSERT((__y1.d1 >> 32) == 0,"MULH96_q8: (__y1.d1 >> 32) == 0");\
		DBG_ASSERT((__y2.d1 >> 32) == 0,"MULH96_q8: (__y2.d1 >> 32) == 0");\
		DBG_ASSERT((__y3.d1 >> 32) == 0,"MULH96_q8: (__y3.d1 >> 32) == 0");\
		DBG_ASSERT((__y4.d1 >> 32) == 0,"MULH96_q8: (__y4.d1 >> 32) == 0");\
		DBG_ASSERT((__y5.d1 >> 32) == 0,"MULH96_q8: (__y5.d1 >> 32) == 0");\
		DBG_ASSERT((__y6.d1 >> 32) == 0,"MULH96_q8: (__y6.d1 >> 32) == 0");\
		DBG_ASSERT((__y7.d1 >> 32) == 0,"MULH96_q8: (__y7.d1 >> 32) == 0");\
		\
		__h0 = (uint64)__x0.d1*(uint64)__y0.d1;\
		__h1 = (uint64)__x1.d1*(uint64)__y1.d1;\
		__h2 = (uint64)__x2.d1*(uint64)__y2.d1;\
		__h3 = (uint64)__x3.d1*(uint64)__y3.d1;\
		__h4 = (uint64)__x4.d1*(uint64)__y4.d1;\
		__h5 = (uint64)__x5.d1*(uint64)__y5.d1;\
		__h6 = (uint64)__x6.d1*(uint64)__y6.d1;\
		__h7 = (uint64)__x7.d1*(uint64)__y7.d1;\
		\
		MULH64(__x0.d0,__y0.d0,       __m0);\
		MULH64(__x1.d0,__y1.d0,       __m1);\
		MULH64(__x2.d0,__y2.d0,       __m2);\
		MULH64(__x3.d0,__y3.d0,       __m3);\
		MULH64(__x4.d0,__y4.d0,       __m4);\
		MULH64(__x5.d0,__y5.d0,       __m5);\
		MULH64(__x6.d0,__y6.d0,       __m6);\
		MULH64(__x7.d0,__y7.d0,       __m7);\
		\
		MUL64x32(__x0.d0,__y0.d1, __a0, __b0);\
		MUL64x32(__x1.d0,__y1.d1, __a1, __b1);\
		MUL64x32(__x2.d0,__y2.d1, __a2, __b2);\
		MUL64x32(__x3.d0,__y3.d1, __a3, __b3);\
		MUL64x32(__x4.d0,__y4.d1, __a4, __b4);\
		MUL64x32(__x5.d0,__y5.d1, __a5, __b5);\
		MUL64x32(__x6.d0,__y6.d1, __a6, __b6);\
		MUL64x32(__x7.d0,__y7.d1, __a7, __b7);\
		\
		MUL64x32(__y0.d0,__x0.d1, __c0, __d0);\
		MUL64x32(__y1.d0,__x1.d1, __c1, __d1);\
		MUL64x32(__y2.d0,__x2.d1, __c2, __d2);\
		MUL64x32(__y3.d0,__x3.d1, __c3, __d3);\
		MUL64x32(__y4.d0,__x4.d1, __c4, __d4);\
		MUL64x32(__y5.d0,__x5.d1, __c5, __d5);\
		MUL64x32(__y6.d0,__x6.d1, __c6, __d6);\
		MUL64x32(__y7.d0,__x7.d1, __c7, __d7);\
		\
		__m0 += __a0;\
		__m1 += __a1;\
		__m2 += __a2;\
		__m3 += __a3;\
		__m4 += __a4;\
		__m5 += __a5;\
		__m6 += __a6;\
		__m7 += __a7;\
		\
		__h0 += __b0 + (__m0 < __a0);\
		__h1 += __b1 + (__m1 < __a1);\
		__h2 += __b2 + (__m2 < __a2);\
		__h3 += __b3 + (__m3 < __a3);\
		__h4 += __b4 + (__m4 < __a4);\
		__h5 += __b5 + (__m5 < __a5);\
		__h6 += __b6 + (__m6 < __a6);\
		__h7 += __b7 + (__m7 < __a7);\
		\
		__m0 += __c0;\
		__m1 += __c1;\
		__m2 += __c2;\
		__m3 += __c3;\
		__m4 += __c4;\
		__m5 += __c5;\
		__m6 += __c6;\
		__m7 += __c7;\
		\
		__h0 += __d0 + (__m0 < __c0);\
		__h1 += __d1 + (__m1 < __c1);\
		__h2 += __d2 + (__m2 < __c2);\
		__h3 += __d3 + (__m3 < __c3);\
		__h4 += __d4 + (__m4 < __c4);\
		__h5 += __d5 + (__m5 < __c5);\
		__h6 += __d6 + (__m6 < __c6);\
		__h7 += __d7 + (__m7 < __c7);\
		\
		/* Now put upper half of the result into __hi: */\
		__hi0.d0 = (__m0>>32) + (__h0<<32);	__hi0.d1 = (uint32)(__h0>>32);\
		__hi1.d0 = (__m1>>32) + (__h1<<32);	__hi1.d1 = (uint32)(__h1>>32);\
		__hi2.d0 = (__m2>>32) + (__h2<<32);	__hi2.d1 = (uint32)(__h2>>32);\
		__hi3.d0 = (__m3>>32) + (__h3<<32);	__hi3.d1 = (uint32)(__h3>>32);\
		__hi4.d0 = (__m4>>32) + (__h4<<32);	__hi4.d1 = (uint32)(__h4>>32);\
		__hi5.d0 = (__m5>>32) + (__h5<<32);	__hi5.d1 = (uint32)(__h5>>32);\
		__hi6.d0 = (__m6>>32) + (__h6<<32);	__hi6.d1 = (uint32)(__h6>>32);\
		__hi7.d0 = (__m7>>32) + (__h7<<32);	__hi7.d1 = (uint32)(__h7>>32);\
	}

#endif

/*****************************************************/
/*                    128 x 128                      */
/*****************************************************/

/* 256-bit square of uint128 *x. Result is returned in a pair of uint128s.

Let x = x0 + x1*2^64,  with x0, x1 < 2^64. Then,

 x^2 = x0^2 + 2*x0*x1*2^64 + x1^2*2^128 .

In terms of the 4 output coefficients, here is what goes into each of those:

	w0 = (x0^2).lo
	w1 = (x0^2).hi + (2*x0*x1).lo
	w2 =             (2*x0*x1).hi + (x1^2).lo
	w3 =                            (x1^2).hi

On Alpha, this needs a total of 6 MUL instructions and 12 ALU ops.

Note that if x < 2^127 (i.e. x.d1 has at most 63 bits), we can feed 2*x.d1 to the
MUL_LOHI and save ourselves one of the two add-and-check-for-carry-hi steps,
for an improved total of 6 MUL and 7 ALU ops.

On Itanium, we can use FMA to reduce the opcount to 6 MUL + 8 ALU op for the 128-bit
version (we only do this for the pipelined 8-operand version) and just 6 MUL + 1 ALU op
for the 127-bit version.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define SQR_LOHI128(__x, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__a,__b;\
		\
		SQR_LOHI64(__x.d0,       &__w0,&__w1);\
		MUL_LOHI64(__x.d0,__x.d1,&__a ,&__b );\
		SQR_LOHI64(__x.d1,       &__w2,&__w3);\
		/* Need to add 2*a*b, so add a*b twice: */\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }

    #define SQR_LOHI127(__x, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__a,__b;\
		\
		SQR_LOHI64(__x.d0,            &__w0,&__w1);\
		MUL_LOHI64(__x.d0,__x.d1 << 1,&__a ,&__b );\
		SQR_LOHI64(__x.d1,            &__w2,&__w3);\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }

    /* 4-operand-pipelined version: */
    #define SQR_LOHI128_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		SQR_LOHI128(__x0, __lo0, __hi0);\
		SQR_LOHI128(__x1, __lo1, __hi1);\
		SQR_LOHI128(__x2, __lo2, __hi2);\
		SQR_LOHI128(__x3, __lo3, __hi3);\
	}
    /* 4-operand-pipelined version, specialized to 127-bit inputs: */
    #define SQR_LOHI127_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		SQR_LOHI127(__x0, __lo0, __hi0);\
		SQR_LOHI127(__x1, __lo1, __hi1);\
		SQR_LOHI127(__x2, __lo2, __hi2);\
		SQR_LOHI127(__x3, __lo3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define SQR_LOHI128_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		SQR_LOHI128(__x0, __lo0, __hi0);\
		SQR_LOHI128(__x1, __lo1, __hi1);\
		SQR_LOHI128(__x2, __lo2, __hi2);\
		SQR_LOHI128(__x3, __lo3, __hi3);\
		SQR_LOHI128(__x4, __lo4, __hi4);\
		SQR_LOHI128(__x5, __lo5, __hi5);\
		SQR_LOHI128(__x6, __lo6, __hi6);\
		SQR_LOHI128(__x7, __lo7, __hi7);\
	}
    /* 8-operand-pipelined version, specialized to 127-bit inputs: */
    #define SQR_LOHI127_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		SQR_LOHI127(__x0, __lo0, __hi0);\
		SQR_LOHI127(__x1, __lo1, __hi1);\
		SQR_LOHI127(__x2, __lo2, __hi2);\
		SQR_LOHI127(__x3, __lo3, __hi3);\
		SQR_LOHI127(__x4, __lo4, __hi4);\
		SQR_LOHI127(__x5, __lo5, __hi5);\
		SQR_LOHI127(__x6, __lo6, __hi6);\
		SQR_LOHI127(__x7, __lo7, __hi7);\
	}

#else

    #define SQR_LOHI128(__x, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__a,__b;\
		\
		SQR_LOHI64(__x.d0,        __w0, __w1);\
		MUL_LOHI64(__x.d0,__x.d1, __a , __b );\
		SQR_LOHI64(__x.d1,        __w2, __w3);\
		/* Need to add 2*a*b, so add a*b twice: */\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }

  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI127(__x, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__a,__b;\
		\
		SQR_LOHI    (__x.d0,                   __w0, __a );\
		MUL_LOHI64_ADD(__x.d0,__x.d1 << 1, __a , __w1, __b );\
		SQR_LOHI64_ADD(__x.d1,             __b , __w2, __w3);\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }
  #else
    #define SQR_LOHI127(__x, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__a,__b;\
		\
		SQR_LOHI64(__x.d0,             __w0, __w1);\
		MUL_LOHI64(__x.d0,__x.d1 << 1, __a , __b );\
		SQR_LOHI64(__x.d1,             __w2, __w3);\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }
  #endif

    /* 4-operand-pipelined version: designed so low 128 bits of output may be returned
    in input arguments, if x and lo happen to point to the same addresses:
    */
  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI128_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__wa0,__wa1,__wa2,__wa3,\
				__wb0,__wb1,__wb2,__wb3,\
				__wc0,__wc1,__wc2,__wc3,\
				__wd0,__wd1,__wd2,__wd3;\
		\
		SQR_LOHI64(__x0.d0,        __wa0, __wb0);\
		SQR_LOHI64(__x1.d0,        __wa1, __wb1);\
		SQR_LOHI64(__x2.d0,        __wa2, __wb2);\
		SQR_LOHI64(__x3.d0,        __wa3, __wb3);\
		\
		MUL_LOHI64(__x0.d0,(uint64)__x0.d1, __a0 , __b0);\
		MUL_LOHI64(__x1.d0,(uint64)__x1.d1, __a1 , __b1);\
		MUL_LOHI64(__x2.d0,(uint64)__x2.d1, __a2 , __b2);\
		MUL_LOHI64(__x3.d0,(uint64)__x3.d1, __a3 , __b3);\
		/* Left-shift both halves of this last 128-bit product one place to effect the mul-by-2: */\
		__hi0.d1 = __b0 >> 63;	__b0 = (__b0 << 1) + (__a0 >> 63);	__a0 <<= 1;\
		__hi1.d1 = __b1 >> 63;	__b1 = (__b1 << 1) + (__a1 >> 63);	__a1 <<= 1;\
		__hi2.d1 = __b2 >> 63;	__b2 = (__b2 << 1) + (__a2 >> 63);	__a2 <<= 1;\
		__hi3.d1 = __b3 >> 63;	__b3 = (__b3 << 1) + (__a3 >> 63);	__a3 <<= 1;\
		\
		/* ADD 2*A TO Wb, PUT CARRYOUT INTO HI.d0: */\
		__wb0 += __a0;	__hi0.d0 = (__wb0 < __a0);\
		__wb1 += __a1;	__hi1.d0 = (__wb1 < __a1);\
		__wb2 += __a2;	__hi2.d0 = (__wb2 < __a2);\
		__wb3 += __a3;	__hi3.d0 = (__wb3 < __a3);\
		\
		/* FUSE ADD OF 2*B IN WITH FINAL SQR_LOHI: */\
		SQR_LOHI64_ADD(__x0.d1, __b0 , __wc0, __wd0);\
		SQR_LOHI64_ADD(__x1.d1, __b1 , __wc1, __wd1);\
		SQR_LOHI64_ADD(__x2.d1, __b2 , __wc2, __wd2);\
		SQR_LOHI64_ADD(__x3.d1, __b3 , __wc3, __wd3);\
		\
		/* Now split the result between __lo and __hi. since hi.d0's are at most 1	\
		here, there should not be any carries resulting from the add of the wc's: */\
		__lo0.d0 = __wa0;	__lo0.d1 = __wb0;	__hi0.d0 +=  __wc0;	__hi0.d1 = __wd0;\
		__lo1.d0 = __wa1;	__lo1.d1 = __wb1;	__hi1.d0 +=  __wc1;	__hi1.d1 = __wd1;\
		__lo2.d0 = __wa2;	__lo2.d1 = __wb2;	__hi2.d0 +=  __wc2;	__hi2.d1 = __wd2;\
		__lo3.d0 = __wa3;	__lo3.d1 = __wb3;	__hi3.d0 +=  __wc3;	__hi3.d1 = __wd3;\
    }
  #else
    #define SQR_LOHI128_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__wa0,__wa1,__wa2,__wa3,\
				__wb0,__wb1,__wb2,__wb3,\
				__wc0,__wc1,__wc2,__wc3,\
				__wd0,__wd1,__wd2,__wd3;\
		\
		SQR_LOHI64(__x0.d0,        __wa0, __wb0);\
		SQR_LOHI64(__x1.d0,        __wa1, __wb1);\
		SQR_LOHI64(__x2.d0,        __wa2, __wb2);\
		SQR_LOHI64(__x3.d0,        __wa3, __wb3);\
		/* On IA64, can add the __wb's in with the __a's here: */\
		MUL_LOHI64(__x0.d0,__x0.d1, __a0 , __b0);\
		MUL_LOHI64(__x1.d0,__x1.d1, __a1 , __b1);\
		MUL_LOHI64(__x2.d0,__x2.d1, __a2 , __b2);\
		MUL_LOHI64(__x3.d0,__x3.d1, __a3 , __b3);\
		/* On IA64, can add the __b's in with the __wc's here: */\
		SQR_LOHI64(__x0.d1,        __wc0, __wd0);\
		SQR_LOHI64(__x1.d1,        __wc1, __wd1);\
		SQR_LOHI64(__x2.d1,        __wc2, __wd2);\
		SQR_LOHI64(__x3.d1,        __wc3, __wd3);\
		\
		/* Need to add 2*a*b, so add a*b twice: */\
		__wb0 += __a0; __wc0 += __b0 + (__wb0 < __a0); __wd0 += (__wc0 < __b0); __wb0 += __a0; __wc0 += __b0 + (__wb0 < __a0); __wd0 += (__wc0 < __b0);\
		__wb1 += __a1; __wc1 += __b1 + (__wb1 < __a1); __wd1 += (__wc1 < __b1); __wb1 += __a1; __wc1 += __b1 + (__wb1 < __a1); __wd1 += (__wc1 < __b1);\
		__wb2 += __a2; __wc2 += __b2 + (__wb2 < __a2); __wd2 += (__wc2 < __b2); __wb2 += __a2; __wc2 += __b2 + (__wb2 < __a2); __wd2 += (__wc2 < __b2);\
		__wb3 += __a3; __wc3 += __b3 + (__wb3 < __a3); __wd3 += (__wc3 < __b3); __wb3 += __a3; __wc3 += __b3 + (__wb3 < __a3); __wd3 += (__wc3 < __b3);\
		\
		/* Now split the result between __lo and __hi: */\
		__lo0.d0 = __wa0;	__lo0.d1 = __wb0;	__hi0.d0 =  __wc0;	__hi0.d1 = __wd0;\
		__lo1.d0 = __wa1;	__lo1.d1 = __wb1;	__hi1.d0 =  __wc1;	__hi1.d1 = __wd1;\
		__lo2.d0 = __wa2;	__lo2.d1 = __wb2;	__hi2.d0 =  __wc2;	__hi2.d1 = __wd2;\
		__lo3.d0 = __wa3;	__lo3.d1 = __wb3;	__hi3.d0 =  __wc3;	__hi3.d1 = __wd3;\
    }
  #endif

  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI127_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3;\
		\
		SQR_LOHI    (__x0.d0,                          __a0,     __b0 );\
		SQR_LOHI    (__x1.d0,                          __a1,     __b1 );\
		SQR_LOHI    (__x2.d0,                          __a2,     __b2 );\
		SQR_LOHI    (__x3.d0,                          __a3,     __b3 );\
		\
		MUL_LOHI64_ADD(__x0.d0, __x0.d1 << 1, __b0 ,     __b0, __hi0.d0 );\
		MUL_LOHI64_ADD(__x1.d0, __x1.d1 << 1, __b1 ,     __b1, __hi1.d0 );\
		MUL_LOHI64_ADD(__x2.d0, __x2.d1 << 1, __b2 ,     __b2, __hi2.d0 );\
		MUL_LOHI64_ADD(__x3.d0, __x3.d1 << 1, __b3 ,     __b3, __hi3.d0 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		SQR_LOHI64_ADD(__x0.d1,           __hi0.d0 , __hi0.d0, __hi0.d1 );\
		SQR_LOHI64_ADD(__x1.d1,           __hi1.d0 , __hi1.d0, __hi1.d1 );\
		SQR_LOHI64_ADD(__x2.d1,           __hi2.d0 , __hi2.d0, __hi2.d1 );\
		SQR_LOHI64_ADD(__x3.d1,           __hi3.d0 , __hi3.d0, __hi3.d1 );\
		/* Delay assignment of lo.d1 terms until here, in case x & lo are the same: */\
		__lo0.d1 = __b0;\
		__lo1.d1 = __b1;\
		__lo2.d1 = __b2;\
		__lo3.d1 = __b3;\
    }
  #else
    #define SQR_LOHI127_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__s0,__s1,__s2,__s3,\
				__t0,__t1,__t2,__t3;\
		\
		SQR_LOHI64(__x0.d0,                   __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                   __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                   __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                   __a3 ,     __b3 );\
		\
		MUL_LOHI64(__x0.d0, __x0.d1 << 1,     __s0 ,     __t0 );\
		MUL_LOHI64(__x1.d0, __x1.d1 << 1,     __s1 ,     __t1 );\
		MUL_LOHI64(__x2.d0, __x2.d1 << 1,     __s2 ,     __t2 );\
		MUL_LOHI64(__x3.d0, __x3.d1 << 1,     __s3 ,     __t3 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		SQR_LOHI64(__x0.d1,               __hi0.d0 , __hi0.d1 );\
		SQR_LOHI64(__x1.d1,               __hi1.d0 , __hi1.d1 );\
		SQR_LOHI64(__x2.d1,               __hi2.d0 , __hi2.d1 );\
		SQR_LOHI64(__x3.d1,               __hi3.d0 , __hi3.d1 );\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__lo0.d1 = __b0 + __s0;	__hi0.d0 += __t0 + (__lo0.d1 < __b0);	__hi0.d1 += (__hi0.d0 < __t0);\
		__lo1.d1 = __b1 + __s1;	__hi1.d0 += __t1 + (__lo1.d1 < __b1);	__hi1.d1 += (__hi1.d0 < __t1);\
		__lo2.d1 = __b2 + __s2;	__hi2.d0 += __t2 + (__lo2.d1 < __b2);	__hi2.d1 += (__hi2.d0 < __t2);\
		__lo3.d1 = __b3 + __s3;	__hi3.d0 += __t3 + (__lo3.d1 < __b3);	__hi3.d1 += (__hi3.d0 < __t3);\
    }
  #endif	/* #endif(ia64) */

    /* 8-operand-pipelined version: designed so low 128 bits of output may be returned
    in input arguments, if x and lo happen to point to the same addresses:
    */
  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI128_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__wa0,__wa1,__wa2,__wa3,__wa4,__wa5,__wa6,__wa7,\
				__wb0,__wb1,__wb2,__wb3,__wb4,__wb5,__wb6,__wb7,\
				__wc0,__wc1,__wc2,__wc3,__wc4,__wc5,__wc6,__wc7,\
				__wd0,__wd1,__wd2,__wd3,__wd4,__wd5,__wd6,__wd7;\
		\
		SQR_LOHI64(__x0.d0,        __wa0, __wb0);\
		SQR_LOHI64(__x1.d0,        __wa1, __wb1);\
		SQR_LOHI64(__x2.d0,        __wa2, __wb2);\
		SQR_LOHI64(__x3.d0,        __wa3, __wb3);\
		SQR_LOHI64(__x4.d0,        __wa4, __wb4);\
		SQR_LOHI64(__x5.d0,        __wa5, __wb5);\
		SQR_LOHI64(__x6.d0,        __wa6, __wb6);\
		SQR_LOHI64(__x7.d0,        __wa7, __wb7);\
		\
		MUL_LOHI64(__x0.d0,__x0.d1, __a0 , __b0);\
		MUL_LOHI64(__x1.d0,__x1.d1, __a1 , __b1);\
		MUL_LOHI64(__x2.d0,__x2.d1, __a2 , __b2);\
		MUL_LOHI64(__x3.d0,__x3.d1, __a3 , __b3);\
		MUL_LOHI64(__x4.d0,__x4.d1, __a4 , __b4);\
		MUL_LOHI64(__x5.d0,__x5.d1, __a5 , __b5);\
		MUL_LOHI64(__x6.d0,__x6.d1, __a6 , __b6);\
		MUL_LOHI64(__x7.d0,__x7.d1, __a7 , __b7);\
		/* Left-shift both halves of this last 128-bit product one place to effect the mul-by-2: */\
		__hi0.d1 = __b0 >> 63;	__b0 = (__b0 << 1) + (__a0 >> 63);	__a0 <<= 1;\
		__hi1.d1 = __b1 >> 63;	__b1 = (__b1 << 1) + (__a1 >> 63);	__a1 <<= 1;\
		__hi2.d1 = __b2 >> 63;	__b2 = (__b2 << 1) + (__a2 >> 63);	__a2 <<= 1;\
		__hi3.d1 = __b3 >> 63;	__b3 = (__b3 << 1) + (__a3 >> 63);	__a3 <<= 1;\
		__hi4.d1 = __b4 >> 63;	__b4 = (__b4 << 1) + (__a4 >> 63);	__a4 <<= 1;\
		__hi5.d1 = __b5 >> 63;	__b5 = (__b5 << 1) + (__a5 >> 63);	__a5 <<= 1;\
		__hi6.d1 = __b6 >> 63;	__b6 = (__b6 << 1) + (__a6 >> 63);	__a6 <<= 1;\
		__hi7.d1 = __b7 >> 63;	__b7 = (__b7 << 1) + (__a7 >> 63);	__a7 <<= 1;\
		\
		/* ADD 2*A TO Wb, PUT CARRYOUT INTO HI.d0: */\
		__wb0 += __a0;	__hi0.d0 = (__wb0 < __a0);\
		__wb1 += __a1;	__hi1.d0 = (__wb1 < __a1);\
		__wb2 += __a2;	__hi2.d0 = (__wb2 < __a2);\
		__wb3 += __a3;	__hi3.d0 = (__wb3 < __a3);\
		__wb4 += __a4;	__hi4.d0 = (__wb4 < __a4);\
		__wb5 += __a5;	__hi5.d0 = (__wb5 < __a5);\
		__wb6 += __a6;	__hi6.d0 = (__wb6 < __a6);\
		__wb7 += __a7;	__hi7.d0 = (__wb7 < __a7);\
		\
		/* FUSE ADD OF 2*B IN WITH FINAL SQR_LOHI: */\
		SQR_LOHI64_ADD(__x0.d1, __b0 , __wc0, __wd0);\
		SQR_LOHI64_ADD(__x1.d1, __b1 , __wc1, __wd1);\
		SQR_LOHI64_ADD(__x2.d1, __b2 , __wc2, __wd2);\
		SQR_LOHI64_ADD(__x3.d1, __b3 , __wc3, __wd3);\
		SQR_LOHI64_ADD(__x4.d1, __b4 , __wc4, __wd4);\
		SQR_LOHI64_ADD(__x5.d1, __b5 , __wc5, __wd5);\
		SQR_LOHI64_ADD(__x6.d1, __b6 , __wc6, __wd6);\
		SQR_LOHI64_ADD(__x7.d1, __b7 , __wc7, __wd7);\
		\
		/* Now split the result between __lo and __hi. since hi.d0's are at most 1	\
		here, there should not be any carries resulting from the add of the wc's: */\
		__lo0.d0 = __wa0;	__lo0.d1 = __wb0;	__hi0.d0 +=  __wc0;	__hi0.d1 = __wd0;\
		__lo1.d0 = __wa1;	__lo1.d1 = __wb1;	__hi1.d0 +=  __wc1;	__hi1.d1 = __wd1;\
		__lo2.d0 = __wa2;	__lo2.d1 = __wb2;	__hi2.d0 +=  __wc2;	__hi2.d1 = __wd2;\
		__lo3.d0 = __wa3;	__lo3.d1 = __wb3;	__hi3.d0 +=  __wc3;	__hi3.d1 = __wd3;\
		__lo4.d0 = __wa4;	__lo4.d1 = __wb4;	__hi4.d0 +=  __wc4;	__hi4.d1 = __wd4;\
		__lo5.d0 = __wa5;	__lo5.d1 = __wb5;	__hi5.d0 +=  __wc5;	__hi5.d1 = __wd5;\
		__lo6.d0 = __wa6;	__lo6.d1 = __wb6;	__hi6.d0 +=  __wc6;	__hi6.d1 = __wd6;\
		__lo7.d0 = __wa7;	__lo7.d1 = __wb7;	__hi7.d0 +=  __wc7;	__hi7.d1 = __wd7;\
    }
  #else
    #define SQR_LOHI128_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__wa0,__wa1,__wa2,__wa3,__wa4,__wa5,__wa6,__wa7,\
				__wb0,__wb1,__wb2,__wb3,__wb4,__wb5,__wb6,__wb7,\
				__wc0,__wc1,__wc2,__wc3,__wc4,__wc5,__wc6,__wc7,\
				__wd0,__wd1,__wd2,__wd3,__wd4,__wd5,__wd6,__wd7;\
		\
		SQR_LOHI64(__x0.d0,        __wa0, __wb0);\
		SQR_LOHI64(__x1.d0,        __wa1, __wb1);\
		SQR_LOHI64(__x2.d0,        __wa2, __wb2);\
		SQR_LOHI64(__x3.d0,        __wa3, __wb3);\
		SQR_LOHI64(__x4.d0,        __wa4, __wb4);\
		SQR_LOHI64(__x5.d0,        __wa5, __wb5);\
		SQR_LOHI64(__x6.d0,        __wa6, __wb6);\
		SQR_LOHI64(__x7.d0,        __wa7, __wb7);\
		/* On IA64, can add the __wb's in with the __a's here: */\
		MUL_LOHI64(__x0.d0,__x0.d1, __a0 , __b0);\
		MUL_LOHI64(__x1.d0,__x1.d1, __a1 , __b1);\
		MUL_LOHI64(__x2.d0,__x2.d1, __a2 , __b2);\
		MUL_LOHI64(__x3.d0,__x3.d1, __a3 , __b3);\
		MUL_LOHI64(__x4.d0,__x4.d1, __a4 , __b4);\
		MUL_LOHI64(__x5.d0,__x5.d1, __a5 , __b5);\
		MUL_LOHI64(__x6.d0,__x6.d1, __a6 , __b6);\
		MUL_LOHI64(__x7.d0,__x7.d1, __a7 , __b7);\
		/* On IA64, can add the __b's in with the __wc's here: */\
		SQR_LOHI64(__x0.d1,        __wc0, __wd0);\
		SQR_LOHI64(__x1.d1,        __wc1, __wd1);\
		SQR_LOHI64(__x2.d1,        __wc2, __wd2);\
		SQR_LOHI64(__x3.d1,        __wc3, __wd3);\
		SQR_LOHI64(__x4.d1,        __wc4, __wd4);\
		SQR_LOHI64(__x5.d1,        __wc5, __wd5);\
		SQR_LOHI64(__x6.d1,        __wc6, __wd6);\
		SQR_LOHI64(__x7.d1,        __wc7, __wd7);\
		\
		/* Need to add 2*a*b, so add a*b twice: */\
		__wb0 += __a0; __wc0 += __b0 + (__wb0 < __a0); __wd0 += (__wc0 < __b0); __wb0 += __a0; __wc0 += __b0 + (__wb0 < __a0); __wd0 += (__wc0 < __b0);\
		__wb1 += __a1; __wc1 += __b1 + (__wb1 < __a1); __wd1 += (__wc1 < __b1); __wb1 += __a1; __wc1 += __b1 + (__wb1 < __a1); __wd1 += (__wc1 < __b1);\
		__wb2 += __a2; __wc2 += __b2 + (__wb2 < __a2); __wd2 += (__wc2 < __b2); __wb2 += __a2; __wc2 += __b2 + (__wb2 < __a2); __wd2 += (__wc2 < __b2);\
		__wb3 += __a3; __wc3 += __b3 + (__wb3 < __a3); __wd3 += (__wc3 < __b3); __wb3 += __a3; __wc3 += __b3 + (__wb3 < __a3); __wd3 += (__wc3 < __b3);\
		__wb4 += __a4; __wc4 += __b4 + (__wb4 < __a4); __wd4 += (__wc4 < __b4); __wb4 += __a4; __wc4 += __b4 + (__wb4 < __a4); __wd4 += (__wc4 < __b4);\
		__wb5 += __a5; __wc5 += __b5 + (__wb5 < __a5); __wd5 += (__wc5 < __b5); __wb5 += __a5; __wc5 += __b5 + (__wb5 < __a5); __wd5 += (__wc5 < __b5);\
		__wb6 += __a6; __wc6 += __b6 + (__wb6 < __a6); __wd6 += (__wc6 < __b6); __wb6 += __a6; __wc6 += __b6 + (__wb6 < __a6); __wd6 += (__wc6 < __b6);\
		__wb7 += __a7; __wc7 += __b7 + (__wb7 < __a7); __wd7 += (__wc7 < __b7); __wb7 += __a7; __wc7 += __b7 + (__wb7 < __a7); __wd7 += (__wc7 < __b7);\
		\
		/* Now split the result between __lo and __hi: */\
		__lo0.d0 = __wa0;	__lo0.d1 = __wb0;	__hi0.d0 =  __wc0;	__hi0.d1 = __wd0;\
		__lo1.d0 = __wa1;	__lo1.d1 = __wb1;	__hi1.d0 =  __wc1;	__hi1.d1 = __wd1;\
		__lo2.d0 = __wa2;	__lo2.d1 = __wb2;	__hi2.d0 =  __wc2;	__hi2.d1 = __wd2;\
		__lo3.d0 = __wa3;	__lo3.d1 = __wb3;	__hi3.d0 =  __wc3;	__hi3.d1 = __wd3;\
		__lo4.d0 = __wa4;	__lo4.d1 = __wb4;	__hi4.d0 =  __wc4;	__hi4.d1 = __wd4;\
		__lo5.d0 = __wa5;	__lo5.d1 = __wb5;	__hi5.d0 =  __wc5;	__hi5.d1 = __wd5;\
		__lo6.d0 = __wa6;	__lo6.d1 = __wb6;	__hi6.d0 =  __wc6;	__hi6.d1 = __wd6;\
		__lo7.d0 = __wa7;	__lo7.d1 = __wb7;	__hi7.d0 =  __wc7;	__hi7.d1 = __wd7;\
    }
  #endif

  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI127_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7;\
		\
		SQR_LOHI    (__x0.d0,                          __a0,     __b0 );\
		SQR_LOHI    (__x1.d0,                          __a1,     __b1 );\
		SQR_LOHI    (__x2.d0,                          __a2,     __b2 );\
		SQR_LOHI    (__x3.d0,                          __a3,     __b3 );\
		SQR_LOHI    (__x4.d0,                          __a4,     __b4 );\
		SQR_LOHI    (__x5.d0,                          __a5,     __b5 );\
		SQR_LOHI    (__x6.d0,                          __a6,     __b6 );\
		SQR_LOHI    (__x7.d0,                          __a7,     __b7 );\
		\
		MUL_LOHI64_ADD(__x0.d0, __x0.d1 << 1, __b0 ,     __b0, __hi0.d0 );\
		MUL_LOHI64_ADD(__x1.d0, __x1.d1 << 1, __b1 ,     __b1, __hi1.d0 );\
		MUL_LOHI64_ADD(__x2.d0, __x2.d1 << 1, __b2 ,     __b2, __hi2.d0 );\
		MUL_LOHI64_ADD(__x3.d0, __x3.d1 << 1, __b3 ,     __b3, __hi3.d0 );\
		MUL_LOHI64_ADD(__x4.d0, __x4.d1 << 1, __b4 ,     __b4, __hi4.d0 );\
		MUL_LOHI64_ADD(__x5.d0, __x5.d1 << 1, __b5 ,     __b5, __hi5.d0 );\
		MUL_LOHI64_ADD(__x6.d0, __x6.d1 << 1, __b6 ,     __b6, __hi6.d0 );\
		MUL_LOHI64_ADD(__x7.d0, __x7.d1 << 1, __b7 ,     __b7, __hi7.d0 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		SQR_LOHI64_ADD(__x0.d1,           __hi0.d0 , __hi0.d0, __hi0.d1 );\
		SQR_LOHI64_ADD(__x1.d1,           __hi1.d0 , __hi1.d0, __hi1.d1 );\
		SQR_LOHI64_ADD(__x2.d1,           __hi2.d0 , __hi2.d0, __hi2.d1 );\
		SQR_LOHI64_ADD(__x3.d1,           __hi3.d0 , __hi3.d0, __hi3.d1 );\
		SQR_LOHI64_ADD(__x4.d1,           __hi4.d0 , __hi4.d0, __hi4.d1 );\
		SQR_LOHI64_ADD(__x5.d1,           __hi5.d0 , __hi5.d0, __hi5.d1 );\
		SQR_LOHI64_ADD(__x6.d1,           __hi6.d0 , __hi6.d0, __hi6.d1 );\
		SQR_LOHI64_ADD(__x7.d1,           __hi7.d0 , __hi7.d0, __hi7.d1 );\
		/* Delay assignment of lo.d1 terms until here, in case x & lo are the same: */\
		__lo0.d1 = __b0;\
		__lo1.d1 = __b1;\
		__lo2.d1 = __b2;\
		__lo3.d1 = __b3;\
		__lo4.d1 = __b4;\
		__lo5.d1 = __b5;\
		__lo6.d1 = __b6;\
		__lo7.d1 = __b7;\
    }
  #else
    #define SQR_LOHI127_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7,\
				__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7;\
		\
		SQR_LOHI64(__x0.d0,                   __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                   __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                   __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                   __a3 ,     __b3 );\
		SQR_LOHI64(__x4.d0,                   __a4 ,     __b4 );\
		SQR_LOHI64(__x5.d0,                   __a5 ,     __b5 );\
		SQR_LOHI64(__x6.d0,                   __a6 ,     __b6 );\
		SQR_LOHI64(__x7.d0,                   __a7 ,     __b7 );\
		\
		MUL_LOHI64(__x0.d0, __x0.d1 << 1,     __s0 ,     __t0 );\
		MUL_LOHI64(__x1.d0, __x1.d1 << 1,     __s1 ,     __t1 );\
		MUL_LOHI64(__x2.d0, __x2.d1 << 1,     __s2 ,     __t2 );\
		MUL_LOHI64(__x3.d0, __x3.d1 << 1,     __s3 ,     __t3 );\
		MUL_LOHI64(__x4.d0, __x4.d1 << 1,     __s4 ,     __t4 );\
		MUL_LOHI64(__x5.d0, __x5.d1 << 1,     __s5 ,     __t5 );\
		MUL_LOHI64(__x6.d0, __x6.d1 << 1,     __s6 ,     __t6 );\
		MUL_LOHI64(__x7.d0, __x7.d1 << 1,     __s7 ,     __t7 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		SQR_LOHI64(__x0.d1,               __hi0.d0 , __hi0.d1 );\
		SQR_LOHI64(__x1.d1,               __hi1.d0 , __hi1.d1 );\
		SQR_LOHI64(__x2.d1,               __hi2.d0 , __hi2.d1 );\
		SQR_LOHI64(__x3.d1,               __hi3.d0 , __hi3.d1 );\
		SQR_LOHI64(__x4.d1,               __hi4.d0 , __hi4.d1 );\
		SQR_LOHI64(__x5.d1,               __hi5.d0 , __hi5.d1 );\
		SQR_LOHI64(__x6.d1,               __hi6.d0 , __hi6.d1 );\
		SQR_LOHI64(__x7.d1,               __hi7.d0 , __hi7.d1 );\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__lo0.d1 = __b0 + __s0;	__hi0.d0 += __t0 + (__lo0.d1 < __b0);	__hi0.d1 += (__hi0.d0 < __t0);\
		__lo1.d1 = __b1 + __s1;	__hi1.d0 += __t1 + (__lo1.d1 < __b1);	__hi1.d1 += (__hi1.d0 < __t1);\
		__lo2.d1 = __b2 + __s2;	__hi2.d0 += __t2 + (__lo2.d1 < __b2);	__hi2.d1 += (__hi2.d0 < __t2);\
		__lo3.d1 = __b3 + __s3;	__hi3.d0 += __t3 + (__lo3.d1 < __b3);	__hi3.d1 += (__hi3.d0 < __t3);\
		__lo4.d1 = __b4 + __s4;	__hi4.d0 += __t4 + (__lo4.d1 < __b4);	__hi4.d1 += (__hi4.d0 < __t4);\
		__lo5.d1 = __b5 + __s5;	__hi5.d0 += __t5 + (__lo5.d1 < __b5);	__hi5.d1 += (__hi5.d0 < __t5);\
		__lo6.d1 = __b6 + __s6;	__hi6.d0 += __t6 + (__lo6.d1 < __b6);	__hi6.d1 += (__hi6.d0 < __t6);\
		__lo7.d1 = __b7 + __s7;	__hi7.d0 += __t7 + (__lo7.d1 < __b7);	__hi7.d1 += (__hi7.d0 < __t7);\
    }
  #endif	/* #endif(ia64) */

#endif

/* 256-bit square of uint128 *x, specialized to the case where x has no more than 96 bits.
Result is returned in a uint128 (low 128 bits) and a uint64 (upper 64 bits.)
On Alpha, this needs a total of 5 MUL instructions and 5 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define SQR_LOHI128_96(__x, __lo128, __hi64)\
    {\
		uint64 __w0,__w1,__w2,__a,__b;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"SQR_LOHI128_96: (__x.d1 >> 32) == 0");\
		SQR_LOHI64(__x.d0,            &__w0,&__w1);\
		/* Need to add 2*a*b, so simply double b (which has at most 32 bits) prior to the MUL_LOHI: */\
		MUL_LOHI64(__x.d0,__x.d1 << 1,&__a ,&__b );\
		__w2  = __x.d1*__x.d1;\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo128.d0 =  __w0;	__lo128.d1 = __w1;\
		__hi64     =  __w2;\
    }

    /* 4-operand-pipelined version: */
    #define SQR_LOHI128_96_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		SQR_LOHI128_96(__x0, __lo0, __hi0);\
		SQR_LOHI128_96(__x1, __lo1, __hi1);\
		SQR_LOHI128_96(__x2, __lo2, __hi2);\
		SQR_LOHI128_96(__x3, __lo3, __hi3);\
	}
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI128_95_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
	{\
		SQR_LOHI128_96_q4(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3)\
	}

    /* 8-operand-pipelined version: */
    #define SQR_LOHI128_96_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		SQR_LOHI128_96(__x0, __lo0, __hi0);\
		SQR_LOHI128_96(__x1, __lo1, __hi1);\
		SQR_LOHI128_96(__x2, __lo2, __hi2);\
		SQR_LOHI128_96(__x3, __lo3, __hi3);\
		SQR_LOHI128_96(__x4, __lo4, __hi4);\
		SQR_LOHI128_96(__x5, __lo5, __hi5);\
		SQR_LOHI128_96(__x6, __lo6, __hi6);\
		SQR_LOHI128_96(__x7, __lo7, __hi7);\
	}
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI128_95_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
	{\
		SQR_LOHI128_96_q8(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3\
		, __x4, __lo4, __hi4\
		, __x5, __lo5, __hi5\
		, __x6, __lo6, __hi6\
		, __x7, __lo7, __hi7)\
	}

#else

  #if 0	/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI128_96(__x, __lo128, __hi64)\
    {\
		uint64 __w0,__w1,__w2,__a,__b;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"SQR_LOHI128_96: (__x.d1 >> 32) == 0");\
		SQR_LOHI64(__x.d0, __w0, __w1);\
		/* Need to add 2*a*b, so simply double b (which has at most 32 bits) prior to the MUL_LOHI: */\
		MUL_LOHI64_ADD(__x.d0, __x.d1 << 1, __w1, __a , __b );\
		__w2  = __x.d1*__x.d1;\
		__lo128.d0 =  __w0;	__lo128.d1 = __a;	__hi64     =  __b + __w2;\
    }
  #else
    #define SQR_LOHI128_96(__x, __lo128, __hi64)\
    {\
		uint64 __w0,__w1,__w2,__a,__b;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"SQR_LOHI128_96: (__x.d1 >> 32) == 0");\
		SQR_LOHI64(__x.d0,           __w0, __w1);\
		/* Need to add 2*a*b, so simply double b (which has at most 32 bits) prior to the MUL_LOHI: */\
		MUL_LOHI64(__x.d0, __x.d1 << 1, __a , __b );\
		__w2  = __x.d1*__x.d1;\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		/* Now split the result between __lo and __hi: */\
		__lo128.d0 =  __w0;	__lo128.d1 = __w1;	__hi64     =  __w2;\
    }
  #endif

    /* 8-operand-pipelined version: designed so low 128 bits of output may be returned
    in input arguments, if x and lo happen to point to the same addresses:
    */
  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI128_96_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x3.d1 >> 32) == 0");\
		\
		SQR_LOHI    (__x0.d0,                             __a0 ,     __b0 );\
		SQR_LOHI    (__x1.d0,                             __a1 ,     __b1 );\
		SQR_LOHI    (__x2.d0,                             __a2 ,     __b2 );\
		SQR_LOHI    (__x3.d0,                             __a3 ,     __b3 );\
		\
		MUL_LOHI64_ADD(__x0.d0, __x0.d1 << 1,     __b0,     __b0 ,    __hi0 );\
		MUL_LOHI64_ADD(__x1.d0, __x1.d1 << 1,     __b1,     __b1 ,    __hi1 );\
		MUL_LOHI64_ADD(__x2.d0, __x2.d1 << 1,     __b2,     __b2 ,    __hi2 );\
		MUL_LOHI64_ADD(__x3.d0, __x3.d1 << 1,     __b3,     __b3 ,    __hi3 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		__hi0 += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1 += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2 += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3 += (uint64)__x3.d1*(uint64)__x3.d1;\
		/* Delay assignment of lo.d1 terms until here, in case x & lo are the same: */\
		__lo0.d1 = __b0;\
		__lo1.d1 = __b1;\
		__lo2.d1 = __b2;\
		__lo3.d1 = __b3;\
    }
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI128_95_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
	{\
		SQR_LOHI128_96_q4(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3)\
	}

  #elif(defined(MULH64_FAST))

    #define SQR_LOHI128_96_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__s0,__s1,__s2,__s3;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x3.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,                    __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                    __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                    __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                    __a3 ,     __b3 );\
		\
		MUL_LOHI64(__x0.d0, __x0.d1 << 1 ,     __s0 ,    __hi0 );\
		MUL_LOHI64(__x1.d0, __x1.d1 << 1 ,     __s1 ,    __hi1 );\
		MUL_LOHI64(__x2.d0, __x2.d1 << 1 ,     __s2 ,    __hi2 );\
		MUL_LOHI64(__x3.d0, __x3.d1 << 1 ,     __s3 ,    __hi3 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		__hi0   += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1   += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2   += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3   += (uint64)__x3.d1*(uint64)__x3.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__lo0.d1 = __b0 + __s0;	__hi0 += (__lo0.d1 < __b0);\
		__lo1.d1 = __b1 + __s1;	__hi1 += (__lo1.d1 < __b1);\
		__lo2.d1 = __b2 + __s2;	__hi2 += (__lo2.d1 < __b2);\
		__lo3.d1 = __b3 + __s3;	__hi3 += (__lo3.d1 < __b3);\
    }
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI128_95_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
	{\
		SQR_LOHI128_96_q4(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3)\
	}

  #else

    /* If slow MULH64, use specialized 64x32-bit cross-product MULs and extra ALU ops
    if inputs are genuinely 96-bit. In this case it pays to use a specialized fast version
    for inputs that are 95 bits or less, which saves 3 ALU ops per input: */
    #define SQR_LOHI128_96_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__s0,__s1,__s2,__s3;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI128_96_q4: (__x3.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,              __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,              __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,              __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,              __a3 ,     __b3 );\
		\
		MUL64x32(__x0.d0, __x0.d1,     __s0 ,    __hi0 );\
		MUL64x32(__x1.d0, __x1.d1,     __s1 ,    __hi1 );\
		MUL64x32(__x2.d0, __x2.d1,     __s2 ,    __hi2 );\
		MUL64x32(__x3.d0, __x3.d1,     __s3 ,    __hi3 );\
		\
		/* Double the cross product term: */\
		__hi0 = (__hi0 << 1) + (__s0 >> 63);\
		__hi1 = (__hi1 << 1) + (__s1 >> 63);\
		__hi2 = (__hi2 << 1) + (__s2 >> 63);\
		__hi3 = (__hi3 << 1) + (__s3 >> 63);\
		\
		__s0 <<= 1;\
		__s1 <<= 1;\
		__s2 <<= 1;\
		__s3 <<= 1;\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		__hi0   += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1   += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2   += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3   += (uint64)__x3.d1*(uint64)__x3.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__lo0.d1 = __b0 + __s0;	__hi0 += (__lo0.d1 < __b0);\
		__lo1.d1 = __b1 + __s1;	__hi1 += (__lo1.d1 < __b1);\
		__lo2.d1 = __b2 + __s2;	__hi2 += (__lo2.d1 < __b2);\
		__lo3.d1 = __b3 + __s3;	__hi3 += (__lo3.d1 < __b3);\
    }

	/* Special version for x < 2^95: */
    #define SQR_LOHI128_95_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__s0,__s1,__s2,__s3;\
		\
		DBG_ASSERT((__x0.d1 >> 31) == 0,"SQR_LOHI128_95_q4: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 31) == 0,"SQR_LOHI128_95_q4: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 31) == 0,"SQR_LOHI128_95_q4: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 31) == 0,"SQR_LOHI128_95_q4: (__x3.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,                    __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                    __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                    __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                    __a3 ,     __b3 );\
		\
		MUL64x32(__x0.d0, __x0.d1 << 1 ,     __s0 ,    __hi0 );\
		MUL64x32(__x1.d0, __x1.d1 << 1 ,     __s1 ,    __hi1 );\
		MUL64x32(__x2.d0, __x2.d1 << 1 ,     __s2 ,    __hi2 );\
		MUL64x32(__x3.d0, __x3.d1 << 1 ,     __s3 ,    __hi3 );\
		\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		\
		__hi0   += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1   += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2   += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3   += (uint64)__x3.d1*(uint64)__x3.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__lo0.d1 = __b0 + __s0;	__hi0 += (__lo0.d1 < __b0);\
		__lo1.d1 = __b1 + __s1;	__hi1 += (__lo1.d1 < __b1);\
		__lo2.d1 = __b2 + __s2;	__hi2 += (__lo2.d1 < __b2);\
		__lo3.d1 = __b3 + __s3;	__hi3 += (__lo3.d1 < __b3);\
    }
  #endif

    /* 8-operand-pipelined version: designed so low 128 bits of output may be returned
    in input arguments, if x and lo happen to point to the same addresses:
    */
  #if 0/*(defined(CPU_IS_IA64))*/
    /* On Itanium, take advantage of fused integer mul/add: */
    #define SQR_LOHI128_96_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x7.d1 >> 32) == 0");\
		\
		SQR_LOHI    (__x0.d0,                             __a0 ,     __b0 );\
		SQR_LOHI    (__x1.d0,                             __a1 ,     __b1 );\
		SQR_LOHI    (__x2.d0,                             __a2 ,     __b2 );\
		SQR_LOHI    (__x3.d0,                             __a3 ,     __b3 );\
		SQR_LOHI    (__x4.d0,                             __a4 ,     __b4 );\
		SQR_LOHI    (__x5.d0,                             __a5 ,     __b5 );\
		SQR_LOHI    (__x6.d0,                             __a6 ,     __b6 );\
		SQR_LOHI    (__x7.d0,                             __a7 ,     __b7 );\
		\
		MUL_LOHI64_ADD(__x0.d0, __x0.d1 << 1,     __b0,     __b0 ,    __hi0 );\
		MUL_LOHI64_ADD(__x1.d0, __x1.d1 << 1,     __b1,     __b1 ,    __hi1 );\
		MUL_LOHI64_ADD(__x2.d0, __x2.d1 << 1,     __b2,     __b2 ,    __hi2 );\
		MUL_LOHI64_ADD(__x3.d0, __x3.d1 << 1,     __b3,     __b3 ,    __hi3 );\
		MUL_LOHI64_ADD(__x4.d0, __x4.d1 << 1,     __b4,     __b4 ,    __hi4 );\
		MUL_LOHI64_ADD(__x5.d0, __x5.d1 << 1,     __b5,     __b5 ,    __hi5 );\
		MUL_LOHI64_ADD(__x6.d0, __x6.d1 << 1,     __b6,     __b6 ,    __hi6 );\
		MUL_LOHI64_ADD(__x7.d0, __x7.d1 << 1,     __b7,     __b7 ,    __hi7 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		__hi0 += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1 += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2 += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3 += (uint64)__x3.d1*(uint64)__x3.d1;\
		__hi4 += (uint64)__x4.d1*(uint64)__x4.d1;\
		__hi5 += (uint64)__x5.d1*(uint64)__x5.d1;\
		__hi6 += (uint64)__x6.d1*(uint64)__x6.d1;\
		__hi7 += (uint64)__x7.d1*(uint64)__x7.d1;\
		/* Delay assignment of lo.d1 terms until here, in case x & lo are the same: */\
		__lo0.d1 = __b0;\
		__lo1.d1 = __b1;\
		__lo2.d1 = __b2;\
		__lo3.d1 = __b3;\
		__lo4.d1 = __b4;\
		__lo5.d1 = __b5;\
		__lo6.d1 = __b6;\
		__lo7.d1 = __b7;\
    }
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI128_95_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
	{\
		SQR_LOHI128_96_q8(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3\
		, __x4, __lo4, __hi4\
		, __x5, __lo5, __hi5\
		, __x6, __lo6, __hi6\
		, __x7, __lo7, __hi7)\
	}

  #elif(defined(MULH64_FAST))

    #define SQR_LOHI128_96_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x7.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,                    __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                    __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                    __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                    __a3 ,     __b3 );\
		SQR_LOHI64(__x4.d0,                    __a4 ,     __b4 );\
		SQR_LOHI64(__x5.d0,                    __a5 ,     __b5 );\
		SQR_LOHI64(__x6.d0,                    __a6 ,     __b6 );\
		SQR_LOHI64(__x7.d0,                    __a7 ,     __b7 );\
		\
		MUL_LOHI64(__x0.d0, __x0.d1 << 1 ,     __s0 ,    __hi0 );\
		MUL_LOHI64(__x1.d0, __x1.d1 << 1 ,     __s1 ,    __hi1 );\
		MUL_LOHI64(__x2.d0, __x2.d1 << 1 ,     __s2 ,    __hi2 );\
		MUL_LOHI64(__x3.d0, __x3.d1 << 1 ,     __s3 ,    __hi3 );\
		MUL_LOHI64(__x4.d0, __x4.d1 << 1 ,     __s4 ,    __hi4 );\
		MUL_LOHI64(__x5.d0, __x5.d1 << 1 ,     __s5 ,    __hi5 );\
		MUL_LOHI64(__x6.d0, __x6.d1 << 1 ,     __s6 ,    __hi6 );\
		MUL_LOHI64(__x7.d0, __x7.d1 << 1 ,     __s7 ,    __hi7 );\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		__hi0   += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1   += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2   += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3   += (uint64)__x3.d1*(uint64)__x3.d1;\
		__hi4   += (uint64)__x4.d1*(uint64)__x4.d1;\
		__hi5   += (uint64)__x5.d1*(uint64)__x5.d1;\
		__hi6   += (uint64)__x6.d1*(uint64)__x6.d1;\
		__hi7   += (uint64)__x7.d1*(uint64)__x7.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__lo0.d1 = __b0 + __s0;	__hi0 += (__lo0.d1 < __b0);\
		__lo1.d1 = __b1 + __s1;	__hi1 += (__lo1.d1 < __b1);\
		__lo2.d1 = __b2 + __s2;	__hi2 += (__lo2.d1 < __b2);\
		__lo3.d1 = __b3 + __s3;	__hi3 += (__lo3.d1 < __b3);\
		__lo4.d1 = __b4 + __s4;	__hi4 += (__lo4.d1 < __b4);\
		__lo5.d1 = __b5 + __s5;	__hi5 += (__lo5.d1 < __b5);\
		__lo6.d1 = __b6 + __s6;	__hi6 += (__lo6.d1 < __b6);\
		__lo7.d1 = __b7 + __s7;	__hi7 += (__lo7.d1 < __b7);\
    }
    /* In this version it doesn't matter whether inputs are 95 or 96-bit: */
    #define SQR_LOHI128_95_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
	{\
		SQR_LOHI128_96_q8(\
		  __x0, __lo0, __hi0\
		, __x1, __lo1, __hi1\
		, __x2, __lo2, __hi2\
		, __x3, __lo3, __hi3\
		, __x4, __lo4, __hi4\
		, __x5, __lo5, __hi5\
		, __x6, __lo6, __hi6\
		, __x7, __lo7, __hi7)\
	}

  #else

    /* If slow MULH64, use specialized 64x32-bit cross-product MULs and extra ALU ops
    if inputs are genuinely 96-bit. In this case it pays to use a specialized fast version
    for inputs that are 95 bits or less, which saves 3 ALU ops per input: */
    #define SQR_LOHI128_96_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7;\
		\
		DBG_ASSERT((__x0.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 32) == 0,"SQR_LOHI128_96_q8: (__x7.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,              __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,              __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,              __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,              __a3 ,     __b3 );\
		SQR_LOHI64(__x4.d0,              __a4 ,     __b4 );\
		SQR_LOHI64(__x5.d0,              __a5 ,     __b5 );\
		SQR_LOHI64(__x6.d0,              __a6 ,     __b6 );\
		SQR_LOHI64(__x7.d0,              __a7 ,     __b7 );\
		\
		MUL64x32(__x0.d0, __x0.d1,     __s0 ,    __hi0 );\
		MUL64x32(__x1.d0, __x1.d1,     __s1 ,    __hi1 );\
		MUL64x32(__x2.d0, __x2.d1,     __s2 ,    __hi2 );\
		MUL64x32(__x3.d0, __x3.d1,     __s3 ,    __hi3 );\
		MUL64x32(__x4.d0, __x4.d1,     __s4 ,    __hi4 );\
		MUL64x32(__x5.d0, __x5.d1,     __s5 ,    __hi5 );\
		MUL64x32(__x6.d0, __x6.d1,     __s6 ,    __hi6 );\
		MUL64x32(__x7.d0, __x7.d1,     __s7 ,    __hi7 );\
		\
		/* Double the cross product term: */\
		__hi0 = (__hi0 << 1) + (__s0 >> 63);\
		__hi1 = (__hi1 << 1) + (__s1 >> 63);\
		__hi2 = (__hi2 << 1) + (__s2 >> 63);\
		__hi3 = (__hi3 << 1) + (__s3 >> 63);\
		__hi4 = (__hi4 << 1) + (__s4 >> 63);\
		__hi5 = (__hi5 << 1) + (__s5 >> 63);\
		__hi6 = (__hi6 << 1) + (__s6 >> 63);\
		__hi7 = (__hi7 << 1) + (__s7 >> 63);\
		\
		__s0 <<= 1;\
		__s1 <<= 1;\
		__s2 <<= 1;\
		__s3 <<= 1;\
		__s4 <<= 1;\
		__s5 <<= 1;\
		__s6 <<= 1;\
		__s7 <<= 1;\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		__hi0   += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1   += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2   += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3   += (uint64)__x3.d1*(uint64)__x3.d1;\
		__hi4   += (uint64)__x4.d1*(uint64)__x4.d1;\
		__hi5   += (uint64)__x5.d1*(uint64)__x5.d1;\
		__hi6   += (uint64)__x6.d1*(uint64)__x6.d1;\
		__hi7   += (uint64)__x7.d1*(uint64)__x7.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__lo0.d1 = __b0 + __s0;	__hi0 += (__lo0.d1 < __b0);\
		__lo1.d1 = __b1 + __s1;	__hi1 += (__lo1.d1 < __b1);\
		__lo2.d1 = __b2 + __s2;	__hi2 += (__lo2.d1 < __b2);\
		__lo3.d1 = __b3 + __s3;	__hi3 += (__lo3.d1 < __b3);\
		__lo4.d1 = __b4 + __s4;	__hi4 += (__lo4.d1 < __b4);\
		__lo5.d1 = __b5 + __s5;	__hi5 += (__lo5.d1 < __b5);\
		__lo6.d1 = __b6 + __s6;	__hi6 += (__lo6.d1 < __b6);\
		__lo7.d1 = __b7 + __s7;	__hi7 += (__lo7.d1 < __b7);\
    }

	/* Special version for x < 2^95: */
    #define SQR_LOHI128_95_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7;\
		\
		DBG_ASSERT((__x0.d1 >> 31) == 0,"SQR_LOHI128_96_q8: (__x0.d1 >> 32) == 0");\
		DBG_ASSERT((__x1.d1 >> 31) == 0,"SQR_LOHI128_96_q8: (__x1.d1 >> 32) == 0");\
		DBG_ASSERT((__x2.d1 >> 31) == 0,"SQR_LOHI128_96_q8: (__x2.d1 >> 32) == 0");\
		DBG_ASSERT((__x3.d1 >> 31) == 0,"SQR_LOHI128_96_q8: (__x3.d1 >> 32) == 0");\
		DBG_ASSERT((__x4.d1 >> 31) == 0,"SQR_LOHI128_96_q8: (__x4.d1 >> 32) == 0");\
		DBG_ASSERT((__x5.d1 >> 31) == 0,"SQR_LOHI128_96_q8: (__x5.d1 >> 32) == 0");\
		DBG_ASSERT((__x6.d1 >> 31) == 0,"SQR_LOHI128_96_q8: (__x6.d1 >> 32) == 0");\
		DBG_ASSERT((__x7.d1 >> 31) == 0,"SQR_LOHI128_96_q8: (__x7.d1 >> 32) == 0");\
		\
		SQR_LOHI64(__x0.d0,                    __a0 ,     __b0 );\
		SQR_LOHI64(__x1.d0,                    __a1 ,     __b1 );\
		SQR_LOHI64(__x2.d0,                    __a2 ,     __b2 );\
		SQR_LOHI64(__x3.d0,                    __a3 ,     __b3 );\
		SQR_LOHI64(__x4.d0,                    __a4 ,     __b4 );\
		SQR_LOHI64(__x5.d0,                    __a5 ,     __b5 );\
		SQR_LOHI64(__x6.d0,                    __a6 ,     __b6 );\
		SQR_LOHI64(__x7.d0,                    __a7 ,     __b7 );\
		\
		MUL64x32(__x0.d0, __x0.d1 << 1 ,     __s0 ,    __hi0 );\
		MUL64x32(__x1.d0, __x1.d1 << 1 ,     __s1 ,    __hi1 );\
		MUL64x32(__x2.d0, __x2.d1 << 1 ,     __s2 ,    __hi2 );\
		MUL64x32(__x3.d0, __x3.d1 << 1 ,     __s3 ,    __hi3 );\
		MUL64x32(__x4.d0, __x4.d1 << 1 ,     __s4 ,    __hi4 );\
		MUL64x32(__x5.d0, __x5.d1 << 1 ,     __s5 ,    __hi5 );\
		MUL64x32(__x6.d0, __x6.d1 << 1 ,     __s6 ,    __hi6 );\
		MUL64x32(__x7.d0, __x7.d1 << 1 ,     __s7 ,    __hi7 );\
		\
		/* Delay assignment of lo.d0 terms until here, in case x & lo are the same: */\
		__lo0.d0 = __a0;\
		__lo1.d0 = __a1;\
		__lo2.d0 = __a2;\
		__lo3.d0 = __a3;\
		__lo4.d0 = __a4;\
		__lo5.d0 = __a5;\
		__lo6.d0 = __a6;\
		__lo7.d0 = __a7;\
		\
		__hi0   += (uint64)__x0.d1*(uint64)__x0.d1;\
		__hi1   += (uint64)__x1.d1*(uint64)__x1.d1;\
		__hi2   += (uint64)__x2.d1*(uint64)__x2.d1;\
		__hi3   += (uint64)__x3.d1*(uint64)__x3.d1;\
		__hi4   += (uint64)__x4.d1*(uint64)__x4.d1;\
		__hi5   += (uint64)__x5.d1*(uint64)__x5.d1;\
		__hi6   += (uint64)__x6.d1*(uint64)__x6.d1;\
		__hi7   += (uint64)__x7.d1*(uint64)__x7.d1;\
		\
		/* __a + 2^64*__b now stores 2*a*b, so only need to add once: */\
		__lo0.d1 = __b0 + __s0;	__hi0 += (__lo0.d1 < __b0);\
		__lo1.d1 = __b1 + __s1;	__hi1 += (__lo1.d1 < __b1);\
		__lo2.d1 = __b2 + __s2;	__hi2 += (__lo2.d1 < __b2);\
		__lo3.d1 = __b3 + __s3;	__hi3 += (__lo3.d1 < __b3);\
		__lo4.d1 = __b4 + __s4;	__hi4 += (__lo4.d1 < __b4);\
		__lo5.d1 = __b5 + __s5;	__hi5 += (__lo5.d1 < __b5);\
		__lo6.d1 = __b6 + __s6;	__hi6 += (__lo6.d1 < __b6);\
		__lo7.d1 = __b7 + __s7;	__hi7 += (__lo7.d1 < __b7);\
    }
  #endif

#endif	/* endif( MUL_LOHI64_SUBROUTINE ) */

/* Upper half of 256-bit product of uint128 x and y, where y < 2^96.
Result is returned in a uint128.
On Alpha, this needs a total of 7 64-bit MULs (but 4 of these are via
MUL64x32s, i.e. are significantly cheaper than full-blown MUL_LOHIs on
32-bit architectures) and 12 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULH128x96(__x, __y, __hi)\
    {\
		uint64 __w1,__w2,__w3,__a,__b,__c,__d,__cy;\
		\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MULH128x96: (__y.d1 >> 32) == 0");\
		\
		MULH64(    __x.d0,__y.d0,       __w1);\
		MUL64x32(  __x.d0,__y.d1,&__a ,&__b );\
		MUL_LOHI64(__y.d0,__x.d1,&__c ,&__d );\
		MUL64x32(  __x.d1,__y.d1,&__w2,&__w3);\
		/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
		__a  += __c;\
		__b  += __d + (__a < __c);\
		DBG_ASSERT((__b >= __d),"MULH128x96: unexpected carryout of __b");\
		/* Now add [w1,w2,w3] + [a,b,0]: */\
		__w1 += __a;\
		__cy  = (__w1 < __a);\
		__w2 += __cy;\
		__w3 += (__w2 < __cy);\
		__w2 += __b;\
		__w3 += (__w2 < __b);\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }

    /* 4-operand-pipelined version: */
    #define MULH128x96_q4(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3)\
    {\
		MULH128x96(__x0, __y0, __hi0);\
		MULH128x96(__x1, __y1, __hi1);\
		MULH128x96(__x2, __y2, __hi2);\
		MULH128x96(__x3, __y3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define MULH128x96_q8(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3\
    , __x4, __y4, __hi4\
    , __x5, __y5, __hi5\
    , __x6, __y6, __hi6\
    , __x7, __y7, __hi7)\
    {\
		MULH128x96(__x0, __y0, __hi0);\
		MULH128x96(__x1, __y1, __hi1);\
		MULH128x96(__x2, __y2, __hi2);\
		MULH128x96(__x3, __y3, __hi3);\
		MULH128x96(__x4, __y4, __hi4);\
		MULH128x96(__x5, __y5, __hi5);\
		MULH128x96(__x6, __y6, __hi6);\
		MULH128x96(__x7, __y7, __hi7);\
	}

#else

    #define MULH128x96(__x, __y, __hi)\
    {\
		uint64 __w1,__w2,__w3,__a,__b,__c,__d,__cy;\
		\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MULH128x96: (__y.d1 >> 32) == 0");\
		\
		MULH64(    __x.d0,__y.d0,       __w1);\
		MUL64x32(  __x.d0,__y.d1, __a , __b );\
		MUL_LOHI64(__y.d0,__x.d1, __c , __d );\
		MUL64x32(  __x.d1,__y.d1, __w2, __w3);\
		/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
		__a  += __c;\
		__b  += __d + (__a < __c);\
		DBG_ASSERT((__b >= __d),"MULH128x96: unexpected carryout of __b");\
		/* Now add [w1,w2,w3] + [a,b,0]: */\
		__w1 += __a;\
		__cy  = (__w1 < __a);\
		__w2 += __cy;\
		__w3 += (__w2 < __cy);\
		__w2 += __b;\
		__w3 += (__w2 < __b);\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }

    /* 4-operand-pipelined version: */
    #define MULH128x96_q4(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3)\
    {\
		uint64 __t0,__a0,__b0,__c0,__d0,__cy0;\
		uint64 __t1,__a1,__b1,__c1,__d1,__cy1;\
		uint64 __t2,__a2,__b2,__c2,__d2,__cy2;\
		uint64 __t3,__a3,__b3,__c3,__d3,__cy3;\
		\
		DBG_ASSERT((__y0.d1 >> 32) == 0,"MULH128x96_q4: (__y0.d1 >> 32) == 0");\
		DBG_ASSERT((__y1.d1 >> 32) == 0,"MULH128x96_q4: (__y1.d1 >> 32) == 0");\
		DBG_ASSERT((__y2.d1 >> 32) == 0,"MULH128x96_q4: (__y2.d1 >> 32) == 0");\
		DBG_ASSERT((__y3.d1 >> 32) == 0,"MULH128x96_q4: (__y3.d1 >> 32) == 0");\
		\
		MULH64(    __x0.d0,__y0.d0, __t0);\
		MULH64(    __x1.d0,__y1.d0, __t1);\
		MULH64(    __x2.d0,__y2.d0, __t2);\
		MULH64(    __x3.d0,__y3.d0, __t3);\
		\
		MUL64x32(  __x0.d0,__y0.d1, __a0, __b0);\
		MUL64x32(  __x1.d0,__y1.d1, __a1, __b1);\
		MUL64x32(  __x2.d0,__y2.d1, __a2, __b2);\
		MUL64x32(  __x3.d0,__y3.d1, __a3, __b3);\
		\
		MUL_LOHI64(__y0.d0,__x0.d1, __c0, __d0);\
		MUL_LOHI64(__y1.d0,__x1.d1, __c1, __d1);\
		MUL_LOHI64(__y2.d0,__x2.d1, __c2, __d2);\
		MUL_LOHI64(__y3.d0,__x3.d1, __c3, __d3);\
		\
		MUL64x32(  __x0.d1,__y0.d1, __hi0.d0, __hi0.d1);\
		MUL64x32(  __x1.d1,__y1.d1, __hi1.d0, __hi1.d1);\
		MUL64x32(  __x2.d1,__y2.d1, __hi2.d0, __hi2.d1);\
		MUL64x32(  __x3.d1,__y3.d1, __hi3.d0, __hi3.d1);\
		\
		/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
		__a0 += __c0;\
		__a1 += __c1;\
		__a2 += __c2;\
		__a3 += __c3;\
		\
		__b0 += __d0 + (__a0 < __c0);\
		__b1 += __d1 + (__a1 < __c1);\
		__b2 += __d2 + (__a2 < __c2);\
		__b3 += __d3 + (__a3 < __c3);\
		\
		/* Now add [w1,w2,w3] + [a,b,0]: */\
		__t0 += __a0;\
		__t1 += __a1;\
		__t2 += __a2;\
		__t3 += __a3;\
		\
		__cy0 = (__t0 < __a0);\
		__cy1 = (__t1 < __a1);\
		__cy2 = (__t2 < __a2);\
		__cy3 = (__t3 < __a3);\
		\
		__hi0.d0 += __cy0;\
		__hi1.d0 += __cy1;\
		__hi2.d0 += __cy2;\
		__hi3.d0 += __cy3;\
		\
		__hi0.d1 += (__hi0.d0 < __cy0);\
		__hi1.d1 += (__hi1.d0 < __cy1);\
		__hi2.d1 += (__hi2.d0 < __cy2);\
		__hi3.d1 += (__hi3.d0 < __cy3);\
		\
		__hi0.d0 += __b0;\
		__hi1.d0 += __b1;\
		__hi2.d0 += __b2;\
		__hi3.d0 += __b3;\
		\
		__hi0.d1 += (__hi0.d0 < __b0);\
		__hi1.d1 += (__hi1.d0 < __b1);\
		__hi2.d1 += (__hi2.d0 < __b2);\
		__hi3.d1 += (__hi3.d0 < __b3);\
    }

    /* 8-operand-pipelined version: */
    #define MULH128x96_q8(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3\
    , __x4, __y4, __hi4\
    , __x5, __y5, __hi5\
    , __x6, __y6, __hi6\
    , __x7, __y7, __hi7)\
    {\
		uint64 __t0,__a0,__b0,__c0,__d0,__cy0;\
		uint64 __t1,__a1,__b1,__c1,__d1,__cy1;\
		uint64 __t2,__a2,__b2,__c2,__d2,__cy2;\
		uint64 __t3,__a3,__b3,__c3,__d3,__cy3;\
		uint64 __t4,__a4,__b4,__c4,__d4,__cy4;\
		uint64 __t5,__a5,__b5,__c5,__d5,__cy5;\
		uint64 __t6,__a6,__b6,__c6,__d6,__cy6;\
		uint64 __t7,__a7,__b7,__c7,__d7,__cy7;\
		\
		DBG_ASSERT((__y0.d1 >> 32) == 0,"MULH128x96_q8: (__y0.d1 >> 32) == 0");\
		DBG_ASSERT((__y1.d1 >> 32) == 0,"MULH128x96_q8: (__y1.d1 >> 32) == 0");\
		DBG_ASSERT((__y2.d1 >> 32) == 0,"MULH128x96_q8: (__y2.d1 >> 32) == 0");\
		DBG_ASSERT((__y3.d1 >> 32) == 0,"MULH128x96_q8: (__y3.d1 >> 32) == 0");\
		DBG_ASSERT((__y4.d1 >> 32) == 0,"MULH128x96_q8: (__y4.d1 >> 32) == 0");\
		DBG_ASSERT((__y5.d1 >> 32) == 0,"MULH128x96_q8: (__y5.d1 >> 32) == 0");\
		DBG_ASSERT((__y6.d1 >> 32) == 0,"MULH128x96_q8: (__y6.d1 >> 32) == 0");\
		DBG_ASSERT((__y7.d1 >> 32) == 0,"MULH128x96_q8: (__y7.d1 >> 32) == 0");\
		\
		MULH64(    __x0.d0,__y0.d0, __t0);\
		MULH64(    __x1.d0,__y1.d0, __t1);\
		MULH64(    __x2.d0,__y2.d0, __t2);\
		MULH64(    __x3.d0,__y3.d0, __t3);\
		MULH64(    __x4.d0,__y4.d0, __t4);\
		MULH64(    __x5.d0,__y5.d0, __t5);\
		MULH64(    __x6.d0,__y6.d0, __t6);\
		MULH64(    __x7.d0,__y7.d0, __t7);\
		\
		MUL64x32(  __x0.d0,__y0.d1, __a0, __b0);\
		MUL64x32(  __x1.d0,__y1.d1, __a1, __b1);\
		MUL64x32(  __x2.d0,__y2.d1, __a2, __b2);\
		MUL64x32(  __x3.d0,__y3.d1, __a3, __b3);\
		MUL64x32(  __x4.d0,__y4.d1, __a4, __b4);\
		MUL64x32(  __x5.d0,__y5.d1, __a5, __b5);\
		MUL64x32(  __x6.d0,__y6.d1, __a6, __b6);\
		MUL64x32(  __x7.d0,__y7.d1, __a7, __b7);\
		\
		MUL_LOHI64(__y0.d0,__x0.d1, __c0, __d0);\
		MUL_LOHI64(__y1.d0,__x1.d1, __c1, __d1);\
		MUL_LOHI64(__y2.d0,__x2.d1, __c2, __d2);\
		MUL_LOHI64(__y3.d0,__x3.d1, __c3, __d3);\
		MUL_LOHI64(__y4.d0,__x4.d1, __c4, __d4);\
		MUL_LOHI64(__y5.d0,__x5.d1, __c5, __d5);\
		MUL_LOHI64(__y6.d0,__x6.d1, __c6, __d6);\
		MUL_LOHI64(__y7.d0,__x7.d1, __c7, __d7);\
		\
		MUL64x32(  __x0.d1,__y0.d1, __hi0.d0, __hi0.d1);\
		MUL64x32(  __x1.d1,__y1.d1, __hi1.d0, __hi1.d1);\
		MUL64x32(  __x2.d1,__y2.d1, __hi2.d0, __hi2.d1);\
		MUL64x32(  __x3.d1,__y3.d1, __hi3.d0, __hi3.d1);\
		MUL64x32(  __x4.d1,__y4.d1, __hi4.d0, __hi4.d1);\
		MUL64x32(  __x5.d1,__y5.d1, __hi5.d0, __hi5.d1);\
		MUL64x32(  __x6.d1,__y6.d1, __hi6.d0, __hi6.d1);\
		MUL64x32(  __x7.d1,__y7.d1, __hi7.d0, __hi7.d1);\
		\
		/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
		__a0 += __c0;\
		__a1 += __c1;\
		__a2 += __c2;\
		__a3 += __c3;\
		__a4 += __c4;\
		__a5 += __c5;\
		__a6 += __c6;\
		__a7 += __c7;\
		\
		__b0 += __d0 + (__a0 < __c0);\
		__b1 += __d1 + (__a1 < __c1);\
		__b2 += __d2 + (__a2 < __c2);\
		__b3 += __d3 + (__a3 < __c3);\
		__b4 += __d4 + (__a4 < __c4);\
		__b5 += __d5 + (__a5 < __c5);\
		__b6 += __d6 + (__a6 < __c6);\
		__b7 += __d7 + (__a7 < __c7);\
		\
		/* Now add [w1,w2,w3] + [a,b,0]: */\
		__t0 += __a0;\
		__t1 += __a1;\
		__t2 += __a2;\
		__t3 += __a3;\
		__t4 += __a4;\
		__t5 += __a5;\
		__t6 += __a6;\
		__t7 += __a7;\
		\
		__cy0 = (__t0 < __a0);\
		__cy1 = (__t1 < __a1);\
		__cy2 = (__t2 < __a2);\
		__cy3 = (__t3 < __a3);\
		__cy4 = (__t4 < __a4);\
		__cy5 = (__t5 < __a5);\
		__cy6 = (__t6 < __a6);\
		__cy7 = (__t7 < __a7);\
		\
		__hi0.d0 += __cy0;\
		__hi1.d0 += __cy1;\
		__hi2.d0 += __cy2;\
		__hi3.d0 += __cy3;\
		__hi4.d0 += __cy4;\
		__hi5.d0 += __cy5;\
		__hi6.d0 += __cy6;\
		__hi7.d0 += __cy7;\
		\
		__hi0.d1 += (__hi0.d0 < __cy0);\
		__hi1.d1 += (__hi1.d0 < __cy1);\
		__hi2.d1 += (__hi2.d0 < __cy2);\
		__hi3.d1 += (__hi3.d0 < __cy3);\
		__hi4.d1 += (__hi4.d0 < __cy4);\
		__hi5.d1 += (__hi5.d0 < __cy5);\
		__hi6.d1 += (__hi6.d0 < __cy6);\
		__hi7.d1 += (__hi7.d0 < __cy7);\
		\
		__hi0.d0 += __b0;\
		__hi1.d0 += __b1;\
		__hi2.d0 += __b2;\
		__hi3.d0 += __b3;\
		__hi4.d0 += __b4;\
		__hi5.d0 += __b5;\
		__hi6.d0 += __b6;\
		__hi7.d0 += __b7;\
		\
		__hi0.d1 += (__hi0.d0 < __b0);\
		__hi1.d1 += (__hi1.d0 < __b1);\
		__hi2.d1 += (__hi2.d0 < __b2);\
		__hi3.d1 += (__hi3.d0 < __b3);\
		__hi4.d1 += (__hi4.d0 < __b4);\
		__hi5.d1 += (__hi5.d0 < __b5);\
		__hi6.d1 += (__hi6.d0 < __b6);\
		__hi7.d1 += (__hi7.d0 < __b7);\
    }

#endif

/* 256-bit product of uint128 *x and *y. Result is returned in a uint256.
On Alpha, this needs a total of 8 MUL instructions.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MUL_LOHI128(__x, __y, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__a,__b,__c,__d;\
		\
		MUL_LOHI64(__x.d0,__y.d0,&__w0,&__w1);\
		MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );\
		MUL_LOHI64(__y.d0,__x.d1,&__c ,&__d );\
		MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		__w1 += __c;\
		__w2 += __d + (__w1 < __c); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __d); /* Overflow into word3 is checked here. */\
		/* Now store the result: */\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
	}

#else

    #define MUL_LOHI128(__x, __y, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__a,__b,__c,__d;\
		\
		MUL_LOHI64(__x.d0,__y.d0, __w0, __w1);\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );\
		MUL_LOHI64(__y.d0,__x.d1, __c , __d );\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		__w1 += __c;\
		__w2 += __d + (__w1 < __c); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __d); /* Overflow into word3 is checked here. */\
		/* Now store the result: */\
		/* Now store the result: */\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
	}


#endif

/* 128-bit product modulo 2^128 of uint128 *x and *y. Result is returned in a uint128.
On Alpha, this needs a total of 4 MUL and 2 ALU op,
and is the same whether x and y have <= 96 or <= 128 bits.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULL128(__x, __y, __lo)\
    {\
		uint64 __w0,__w1,__a,__c;\
		\
		__a = __MULL64(__x.d0,__y.d1);\
		__c = __MULL64(__y.d0,__x.d1);\
		MUL_LOHI64(__x.d0,__y.d0,&__w0,&__w1);\
		__w1 += __a;\
		__w1 += __c;\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
    }

    /* 4-operand-pipelined version: */
    #define MULL128_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
    	MULL128(__x0, __y0, __lo0);\
    	MULL128(__x1, __y1, __lo1);\
    	MULL128(__x2, __y2, __lo2);\
    	MULL128(__x3, __y3, __lo3);\
	}
    /* 4-operand-pipelined version, returning result in the first input: */
    #define MULL128_INPLACE_q4(\
      __x0, __y0\
    , __x1, __y1\
    , __x2, __y2\
    , __x3, __y3)\
    {\
    	MULL128(__x0, __y0, __x0);\
    	MULL128(__x1, __y1, __x1);\
    	MULL128(__x2, __y2, __x2);\
    	MULL128(__x3, __y3, __x3);\
	}

    /* 8-operand-pipelined version: */
    #define MULL128_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
    	MULL128(__x0, __y0, __lo0);\
    	MULL128(__x1, __y1, __lo1);\
    	MULL128(__x2, __y2, __lo2);\
    	MULL128(__x3, __y3, __lo3);\
    	MULL128(__x4, __y4, __lo4);\
    	MULL128(__x5, __y5, __lo5);\
    	MULL128(__x6, __y6, __lo6);\
    	MULL128(__x7, __y7, __lo7);\
	}
    /* 8-operand-pipelined version, returning result in the first input: */
    #define MULL128_INPLACE_q8(\
      __x0, __y0\
    , __x1, __y1\
    , __x2, __y2\
    , __x3, __y3\
    , __x4, __y4\
    , __x5, __y5\
    , __x6, __y6\
    , __x7, __y7)\
    {\
    	MULL128(__x0, __y0, __x0);\
    	MULL128(__x1, __y1, __x1);\
    	MULL128(__x2, __y2, __x2);\
    	MULL128(__x3, __y3, __x3);\
    	MULL128(__x4, __y4, __x4);\
    	MULL128(__x5, __y5, __x5);\
    	MULL128(__x6, __y6, __x6);\
    	MULL128(__x7, __y7, __x7);\
	}

#else

    #define MULL128(__x, __y, __lo)\
    {\
		uint64 __w0,__w1,__a,__c;\
		\
		MULL64(__x.d0,__y.d1,__a);\
		MULL64(__y.d0,__x.d1,__c);\
		MUL_LOHI64(__x.d0,__y.d0, __w0, __w1);\
		__w1 += __a;\
		__w1 += __c;\
		__lo.d0 =  __w0;	__lo.d1 = __w1;\
    }

    /* 4-operand-pipelined version: */
    #define MULL128_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
		uint64	__v0,__v1,__v2,__v3,\
				__w0,__w1,__w2,__w3;\
		\
		MULL64(__x0.d0,__y0.d1,__v0);\
		MULL64(__x1.d0,__y1.d1,__v1);\
		MULL64(__x2.d0,__y2.d1,__v2);\
		MULL64(__x3.d0,__y3.d1,__v3);\
		\
		MULL64(__y0.d0,__x0.d1,__w0);\
		MULL64(__y1.d0,__x1.d1,__w1);\
		MULL64(__y2.d0,__x2.d1,__w2);\
		MULL64(__y3.d0,__x3.d1,__w3);\
		\
		__lo0.d1 = __v0 + __w0;\
		__lo1.d1 = __v1 + __w1;\
		__lo2.d1 = __v2 + __w2;\
		__lo3.d1 = __v3 + __w3;\
		\
		MUL_LOHI64(__x0.d0,__y0.d0, __v0, __w0);\
		MUL_LOHI64(__x1.d0,__y1.d0, __v1, __w1);\
		MUL_LOHI64(__x2.d0,__y2.d0, __v2, __w2);\
		MUL_LOHI64(__x3.d0,__y3.d0, __v3, __w3);\
		\
		__lo0.d0 = __v0;\
		__lo1.d0 = __v1;\
		__lo2.d0 = __v2;\
		__lo3.d0 = __v3;\
		\
		__lo0.d1 += __w0;\
		__lo1.d1 += __w1;\
		__lo2.d1 += __w2;\
		__lo3.d1 += __w3;\
    }
    /* 4-operand-pipelined version, returning result in the first input: */
    #define MULL128_INPLACE_q4(\
      __x0, __y0\
    , __x1, __y1\
    , __x2, __y2\
    , __x3, __y3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__v0,__v1,__v2,__v3,\
				__w0,__w1,__w2,__w3;\
		\
		MULL64(__x0.d0,__y0.d1,__v0);\
		MULL64(__x1.d0,__y1.d1,__v1);\
		MULL64(__x2.d0,__y2.d1,__v2);\
		MULL64(__x3.d0,__y3.d1,__v3);\
		\
		MULL64(__y0.d0,__x0.d1,__w0);\
		MULL64(__y1.d0,__x1.d1,__w1);\
		MULL64(__y2.d0,__x2.d1,__w2);\
		MULL64(__y3.d0,__x3.d1,__w3);\
		\
		__a0 = __v0 + __w0;\
		__a1 = __v1 + __w1;\
		__a2 = __v2 + __w2;\
		__a3 = __v3 + __w3;\
		\
		MUL_LOHI64(__x0.d0,__y0.d0, __v0, __w0);\
		MUL_LOHI64(__x1.d0,__y1.d0, __v1, __w1);\
		MUL_LOHI64(__x2.d0,__y2.d0, __v2, __w2);\
		MUL_LOHI64(__x3.d0,__y3.d0, __v3, __w3);\
		\
		__x0.d0 = __v0;\
		__x1.d0 = __v1;\
		__x2.d0 = __v2;\
		__x3.d0 = __v3;\
		\
		__x0.d1 = __w0 + __a0;\
		__x1.d1 = __w1 + __a1;\
		__x2.d1 = __w2 + __a2;\
		__x3.d1 = __w3 + __a3;\
   }


    /* 8-operand-pipelined version: */
    #define MULL128_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
		uint64	__v0,__v1,__v2,__v3,__v4,__v5,__v6,__v7,\
				__w0,__w1,__w2,__w3,__w4,__w5,__w6,__w7;\
		\
		MULL64(__x0.d0,__y0.d1,__v0);\
		MULL64(__x1.d0,__y1.d1,__v1);\
		MULL64(__x2.d0,__y2.d1,__v2);\
		MULL64(__x3.d0,__y3.d1,__v3);\
		MULL64(__x4.d0,__y4.d1,__v4);\
		MULL64(__x5.d0,__y5.d1,__v5);\
		MULL64(__x6.d0,__y6.d1,__v6);\
		MULL64(__x7.d0,__y7.d1,__v7);\
		\
		MULL64(__y0.d0,__x0.d1,__w0);\
		MULL64(__y1.d0,__x1.d1,__w1);\
		MULL64(__y2.d0,__x2.d1,__w2);\
		MULL64(__y3.d0,__x3.d1,__w3);\
		MULL64(__y4.d0,__x4.d1,__w4);\
		MULL64(__y5.d0,__x5.d1,__w5);\
		MULL64(__y6.d0,__x6.d1,__w6);\
		MULL64(__y7.d0,__x7.d1,__w7);\
		\
		__lo0.d1 = __v0 + __w0;\
		__lo1.d1 = __v1 + __w1;\
		__lo2.d1 = __v2 + __w2;\
		__lo3.d1 = __v3 + __w3;\
		__lo4.d1 = __v4 + __w4;\
		__lo5.d1 = __v5 + __w5;\
		__lo6.d1 = __v6 + __w6;\
		__lo7.d1 = __v7 + __w7;\
		\
		MUL_LOHI64(__x0.d0,__y0.d0, __v0, __w0);\
		MUL_LOHI64(__x1.d0,__y1.d0, __v1, __w1);\
		MUL_LOHI64(__x2.d0,__y2.d0, __v2, __w2);\
		MUL_LOHI64(__x3.d0,__y3.d0, __v3, __w3);\
		MUL_LOHI64(__x4.d0,__y4.d0, __v4, __w4);\
		MUL_LOHI64(__x5.d0,__y5.d0, __v5, __w5);\
		MUL_LOHI64(__x6.d0,__y6.d0, __v6, __w6);\
		MUL_LOHI64(__x7.d0,__y7.d0, __v7, __w7);\
		\
		__lo0.d0 = __v0;\
		__lo1.d0 = __v1;\
		__lo2.d0 = __v2;\
		__lo3.d0 = __v3;\
		__lo4.d0 = __v4;\
		__lo5.d0 = __v5;\
		__lo6.d0 = __v6;\
		__lo7.d0 = __v7;\
		\
		__lo0.d1 += __w0;\
		__lo1.d1 += __w1;\
		__lo2.d1 += __w2;\
		__lo3.d1 += __w3;\
		__lo4.d1 += __w4;\
		__lo5.d1 += __w5;\
		__lo6.d1 += __w6;\
		__lo7.d1 += __w7;\
    }
    /* 8-operand-pipelined version, returning result in the first input: */
    #define MULL128_INPLACE_q8(\
      __x0, __y0\
    , __x1, __y1\
    , __x2, __y2\
    , __x3, __y3\
    , __x4, __y4\
    , __x5, __y5\
    , __x6, __y6\
    , __x7, __y7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__v0,__v1,__v2,__v3,__v4,__v5,__v6,__v7,\
				__w0,__w1,__w2,__w3,__w4,__w5,__w6,__w7;\
		\
		MULL64(__x0.d0,__y0.d1,__v0);\
		MULL64(__x1.d0,__y1.d1,__v1);\
		MULL64(__x2.d0,__y2.d1,__v2);\
		MULL64(__x3.d0,__y3.d1,__v3);\
		MULL64(__x4.d0,__y4.d1,__v4);\
		MULL64(__x5.d0,__y5.d1,__v5);\
		MULL64(__x6.d0,__y6.d1,__v6);\
		MULL64(__x7.d0,__y7.d1,__v7);\
		\
		MULL64(__y0.d0,__x0.d1,__w0);\
		MULL64(__y1.d0,__x1.d1,__w1);\
		MULL64(__y2.d0,__x2.d1,__w2);\
		MULL64(__y3.d0,__x3.d1,__w3);\
		MULL64(__y4.d0,__x4.d1,__w4);\
		MULL64(__y5.d0,__x5.d1,__w5);\
		MULL64(__y6.d0,__x6.d1,__w6);\
		MULL64(__y7.d0,__x7.d1,__w7);\
		\
		__a0 = __v0 + __w0;\
		__a1 = __v1 + __w1;\
		__a2 = __v2 + __w2;\
		__a3 = __v3 + __w3;\
		__a4 = __v4 + __w4;\
		__a5 = __v5 + __w5;\
		__a6 = __v6 + __w6;\
		__a7 = __v7 + __w7;\
		\
		MUL_LOHI64(__x0.d0,__y0.d0, __v0, __w0);\
		MUL_LOHI64(__x1.d0,__y1.d0, __v1, __w1);\
		MUL_LOHI64(__x2.d0,__y2.d0, __v2, __w2);\
		MUL_LOHI64(__x3.d0,__y3.d0, __v3, __w3);\
		MUL_LOHI64(__x4.d0,__y4.d0, __v4, __w4);\
		MUL_LOHI64(__x5.d0,__y5.d0, __v5, __w5);\
		MUL_LOHI64(__x6.d0,__y6.d0, __v6, __w6);\
		MUL_LOHI64(__x7.d0,__y7.d0, __v7, __w7);\
		\
		__x0.d0 = __v0;\
		__x1.d0 = __v1;\
		__x2.d0 = __v2;\
		__x3.d0 = __v3;\
		__x4.d0 = __v4;\
		__x5.d0 = __v5;\
		__x6.d0 = __v6;\
		__x7.d0 = __v7;\
		\
		__x0.d1 = __w0 + __a0;\
		__x1.d1 = __w1 + __a1;\
		__x2.d1 = __w2 + __a2;\
		__x3.d1 = __w3 + __a3;\
		__x4.d1 = __w4 + __a4;\
		__x5.d1 = __w5 + __a5;\
		__x6.d1 = __w6 + __a6;\
		__x7.d1 = __w7 + __a7;\
    }

#endif

/* Upper half of 256-bit product of uint128 *x and *y. Result is returned in a uint128.
On Alpha, this needs a total of 7 MUL instructions and 12 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULH128(__x, __y, __hi)\
    {\
	uint64 __w1,__w2,__w3,__a,__b,__c,__d,__cy;\
	\
	MULH64(    __x.d0,__y.d0,       __w1);\
	MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );\
	MUL_LOHI64(__y.d0,__x.d1,&__c ,&__d );\
	MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);\
	/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
	__a  += __c;\
	__b  += __d + (__a < __c);\
	DBG_ASSERT((__b >= __d),"MULH128: unexpected carryout of __b");\
	/* Now add [w1,w2,w3] + [a,b,0]: */\
	__w1 += __a;\
	__cy  = (__w1 < __a);\
	__w2 += __cy;\
	__w3 += (__w2 < __cy);\
	__w2 += __b;\
	__w3 += (__w2 < __b);\
	__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }

    /* 4-operand-pipelined version: */
    #define MULH128_q4(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3)\
    {\
    	MULH128(__x0, __y0, __hi0);\
    	MULH128(__x1, __y1, __hi1);\
    	MULH128(__x2, __y2, __hi2);\
    	MULH128(__x3, __y3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define MULH128_q8(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3\
    , __x4, __y4, __hi4\
    , __x5, __y5, __hi5\
    , __x6, __y6, __hi6\
    , __x7, __y7, __hi7)\
    {\
    	MULH128(__x0, __y0, __hi0);\
    	MULH128(__x1, __y1, __hi1);\
    	MULH128(__x2, __y2, __hi2);\
    	MULH128(__x3, __y3, __hi3);\
    	MULH128(__x4, __y4, __hi4);\
    	MULH128(__x5, __y5, __hi5);\
    	MULH128(__x6, __y6, __hi6);\
    	MULH128(__x7, __y7, __hi7);\
	}

/* Routine to perform fused version of key 3-operation sequence in 128-bit modmul:

	SQR_LOHI128(x,lo,hi);
	MULL128(lo,qinv,lo);
	MULH128(q,lo,lo);
*/
    /* 4-operand-pipelined version: */
    #define THREE_OP128_q4(\
      __x0, __qinv0, __q0, __lo0, __hi0\
    , __x1, __qinv1, __q1, __lo1, __hi1\
    , __x2, __qinv2, __q2, __lo2, __hi2\
    , __x3, __qinv3, __q3, __lo3, __hi3)\
    {\
		uint64	__a0,__b0,__t0,__wa0,__wb0,__wc0,__wd0,__cy0;\
		uint64	__a1,__b1,__t1,__wa1,__wb1,__wc1,__wd1,__cy1;\
		uint64	__a2,__b2,__t2,__wa2,__wb2,__wc2,__wd2,__cy2;\
		uint64	__a3,__b3,__t3,__wa3,__wb3,__wc3,__wd3,__cy3;\
		\
		/* SQR_LOHI128(x,lo,hi) IS HERE: */\
		SQR_LOHI64(__x0.d0,        &__wa0,&__wb0);\
		SQR_LOHI64(__x1.d0,        &__wa1,&__wb1);\
		SQR_LOHI64(__x2.d0,        &__wa2,&__wb2);\
		SQR_LOHI64(__x3.d0,        &__wa3,&__wb3);\
		/* On IA64, can add the __wb's in with the __a's here: */\
		MUL_LOHI64(__x0.d0,__x0.d1,&__a0 ,&__b0);\
		MUL_LOHI64(__x1.d0,__x1.d1,&__a1 ,&__b1);\
		MUL_LOHI64(__x2.d0,__x2.d1,&__a2 ,&__b2);\
		MUL_LOHI64(__x3.d0,__x3.d1,&__a3 ,&__b3);\
		/* On IA64, can add the __b's in with the __wc's here: */\
		SQR_LOHI64(__x0.d1,        &__wc0,&__wd0);\
		SQR_LOHI64(__x1.d1,        &__wc1,&__wd1);\
		SQR_LOHI64(__x2.d1,        &__wc2,&__wd2);\
		SQR_LOHI64(__x3.d1,        &__wc3,&__wd3);\
		\
		/* Need to add 2*a*b, so add a*b twice: */\
		__wb0 += __a0;	__wc0 += __b0 + (__wb0 < __a0);	__wd0 += (__wc0 < __b0);	__wb0 += __a0;	__wc0 += __b0 + (__wb0 < __a0);	__wd0 += (__wc0 < __b0);\
		__wb1 += __a1;	__wc1 += __b1 + (__wb1 < __a1);	__wd1 += (__wc1 < __b1);	__wb1 += __a1;	__wc1 += __b1 + (__wb1 < __a1);	__wd1 += (__wc1 < __b1);\
		__wb2 += __a2;	__wc2 += __b2 + (__wb2 < __a2);	__wd2 += (__wc2 < __b2);	__wb2 += __a2;	__wc2 += __b2 + (__wb2 < __a2);	__wd2 += (__wc2 < __b2);\
		__wb3 += __a3;	__wc3 += __b3 + (__wb3 < __a3);	__wd3 += (__wc3 < __b3);	__wb3 += __a3;	__wc3 += __b3 + (__wb3 < __a3);	__wd3 += (__wc3 < __b3);\
		\
		/* Now store the high 128 bits of the result in __hi: */\
		__hi0.d0 =  __wc0;	__hi0.d1 = __wd0;\
		__hi1.d0 =  __wc1;	__hi1.d1 = __wd1;\
		__hi2.d0 =  __wc2;	__hi2.d1 = __wd2;\
		__hi3.d0 =  __wc3;	__hi3.d1 = __wd3;\
		\
		/* MULL128(lo,qinv,lo) IS HERE: wc, wd now free, so use as 64-bit temporaries: */\
		\
		__a0 = __wa0*__qinv0.d1 + __qinv0.d0*__wb0;\
		__a1 = __wa1*__qinv1.d1 + __qinv1.d0*__wb1;\
		__a2 = __wa2*__qinv2.d1 + __qinv2.d0*__wb2;\
		__a3 = __wa3*__qinv3.d1 + __qinv3.d0*__wb3;\
		\
		MUL_LOHI64(__wa0,__qinv0.d0,&__wc0,&__wd0);\
		MUL_LOHI64(__wa1,__qinv1.d0,&__wc1,&__wd1);\
		MUL_LOHI64(__wa2,__qinv2.d0,&__wc2,&__wd2);\
		MUL_LOHI64(__wa3,__qinv3.d0,&__wc3,&__wd3);\
		\
		__wa0 = __wc0;\
		__wa1 = __wc1;\
		__wa2 = __wc2;\
		__wa3 = __wc3;\
		\
		__wb0 = __wd0 + __a0;\
		__wb1 = __wd1 + __a1;\
		__wb2 = __wd2 + __a2;\
		__wb3 = __wd3 + __a3;\
		\
		/* MULH128(q,lo,lo) IS HERE: store result in __lo: */\
		\
		MULH64(__wa0,__q0.d0, __t0);\
		MULH64(__wa1,__q1.d0, __t1);\
		MULH64(__wa2,__q2.d0, __t2);\
		MULH64(__wa3,__q3.d0, __t3);\
		\
		MUL_LOHI64(__wa0,__q0.d1, &__a0, &__b0);\
		MUL_LOHI64(__wa1,__q1.d1, &__a1, &__b1);\
		MUL_LOHI64(__wa2,__q2.d1, &__a2, &__b2);\
		MUL_LOHI64(__wa3,__q3.d1, &__a3, &__b3);\
		\
		MUL_LOHI64(__q0.d0,__wb0,&__wc0,&__wd0);\
		MUL_LOHI64(__q1.d0,__wb1,&__wc1,&__wd1);\
		MUL_LOHI64(__q2.d0,__wb2,&__wc2,&__wd2);\
		MUL_LOHI64(__q3.d0,__wb3,&__wc3,&__wd3);\
		\
		MUL_LOHI64(__wb0,__q0.d1,&__lo0.d0,&__lo0.d1);\
		MUL_LOHI64(__wb1,__q1.d1,&__lo1.d0,&__lo1.d1);\
		MUL_LOHI64(__wb2,__q2.d1,&__lo2.d0,&__lo2.d1);\
		MUL_LOHI64(__wb3,__q3.d1,&__lo3.d0,&__lo3.d1);\
		\
		/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
		__a0 += __wc0;\
		__a1 += __wc1;\
		__a2 += __wc2;\
		__a3 += __wc3;\
		\
		__b0 += __wd0 + (__a0 < __wc0);\
		__b1 += __wd1 + (__a1 < __wc1);\
		__b2 += __wd2 + (__a2 < __wc2);\
		__b3 += __wd3 + (__a3 < __wc3);\
		\
		/* Now add [w1,w2,w3] + [a,b,0]: */\
		__t0 += __a0;\
		__t1 += __a1;\
		__t2 += __a2;\
		__t3 += __a3;\
		\
		__cy0 = (__t0 < __a0);\
		__cy1 = (__t1 < __a1);\
		__cy2 = (__t2 < __a2);\
		__cy3 = (__t3 < __a3);\
		\
		__lo0.d0 += __cy0;\
		__lo1.d0 += __cy1;\
		__lo2.d0 += __cy2;\
		__lo3.d0 += __cy3;\
		\
		__lo0.d1 += (__lo0.d0 < __cy0);\
		__lo1.d1 += (__lo1.d0 < __cy1);\
		__lo2.d1 += (__lo2.d0 < __cy2);\
		__lo3.d1 += (__lo3.d0 < __cy3);\
		\
		__lo0.d0 += __b0;\
		__lo1.d0 += __b1;\
		__lo2.d0 += __b2;\
		__lo3.d0 += __b3;\
		\
		__lo0.d1 += (__lo0.d0 < __b0);\
		__lo1.d1 += (__lo1.d0 < __b1);\
		__lo2.d1 += (__lo2.d0 < __b2);\
		__lo3.d1 += (__lo3.d0 < __b3);\
    }

    /* 8-operand-pipelined version: */
    #define THREE_OP128_q8(\
      __x0, __qinv0, __q0, __lo0, __hi0\
    , __x1, __qinv1, __q1, __lo1, __hi1\
    , __x2, __qinv2, __q2, __lo2, __hi2\
    , __x3, __qinv3, __q3, __lo3, __hi3\
    , __x4, __qinv4, __q4, __lo4, __hi4\
    , __x5, __qinv5, __q5, __lo5, __hi5\
    , __x6, __qinv6, __q6, __lo6, __hi6\
    , __x7, __qinv7, __q7, __lo7, __hi7)\
    {\
		uint64	__a0,__b0,__t0,__wa0,__wb0,__wc0,__wd0,__cy0;\
		uint64	__a1,__b1,__t1,__wa1,__wb1,__wc1,__wd1,__cy1;\
		uint64	__a2,__b2,__t2,__wa2,__wb2,__wc2,__wd2,__cy2;\
		uint64	__a3,__b3,__t3,__wa3,__wb3,__wc3,__wd3,__cy3;\
		uint64	__a4,__b4,__t4,__wa4,__wb4,__wc4,__wd4,__cy4;\
		uint64	__a5,__b5,__t5,__wa5,__wb5,__wc5,__wd5,__cy5;\
		uint64	__a6,__b6,__t6,__wa6,__wb6,__wc6,__wd6,__cy6;\
		uint64	__a7,__b7,__t7,__wa7,__wb7,__wc7,__wd7,__cy7;\
		\
		/* SQR_LOHI128(x,lo,hi) IS HERE: */\
		SQR_LOHI64(__x0.d0,        &__wa0,&__wb0);\
		SQR_LOHI64(__x1.d0,        &__wa1,&__wb1);\
		SQR_LOHI64(__x2.d0,        &__wa2,&__wb2);\
		SQR_LOHI64(__x3.d0,        &__wa3,&__wb3);\
		SQR_LOHI64(__x4.d0,        &__wa4,&__wb4);\
		SQR_LOHI64(__x5.d0,        &__wa5,&__wb5);\
		SQR_LOHI64(__x6.d0,        &__wa6,&__wb6);\
		SQR_LOHI64(__x7.d0,        &__wa7,&__wb7);\
		/* On IA64, can add the __wb's in with the __a's here: */\
		MUL_LOHI64(__x0.d0,__x0.d1,&__a0 ,&__b0);\
		MUL_LOHI64(__x1.d0,__x1.d1,&__a1 ,&__b1);\
		MUL_LOHI64(__x2.d0,__x2.d1,&__a2 ,&__b2);\
		MUL_LOHI64(__x3.d0,__x3.d1,&__a3 ,&__b3);\
		MUL_LOHI64(__x4.d0,__x4.d1,&__a4 ,&__b4);\
		MUL_LOHI64(__x5.d0,__x5.d1,&__a5 ,&__b5);\
		MUL_LOHI64(__x6.d0,__x6.d1,&__a6 ,&__b6);\
		MUL_LOHI64(__x7.d0,__x7.d1,&__a7 ,&__b7);\
		/* On IA64, can add the __b's in with the __wc's here: */\
		SQR_LOHI64(__x0.d1,        &__wc0,&__wd0);\
		SQR_LOHI64(__x1.d1,        &__wc1,&__wd1);\
		SQR_LOHI64(__x2.d1,        &__wc2,&__wd2);\
		SQR_LOHI64(__x3.d1,        &__wc3,&__wd3);\
		SQR_LOHI64(__x4.d1,        &__wc4,&__wd4);\
		SQR_LOHI64(__x5.d1,        &__wc5,&__wd5);\
		SQR_LOHI64(__x6.d1,        &__wc6,&__wd6);\
		SQR_LOHI64(__x7.d1,        &__wc7,&__wd7);\
		\
		/* Need to add 2*a*b, so add a*b twice: */\
		__wb0 += __a0;	__wc0 += __b0 + (__wb0 < __a0);	__wd0 += (__wc0 < __b0);	__wb0 += __a0;	__wc0 += __b0 + (__wb0 < __a0);	__wd0 += (__wc0 < __b0);\
		__wb1 += __a1;	__wc1 += __b1 + (__wb1 < __a1);	__wd1 += (__wc1 < __b1);	__wb1 += __a1;	__wc1 += __b1 + (__wb1 < __a1);	__wd1 += (__wc1 < __b1);\
		__wb2 += __a2;	__wc2 += __b2 + (__wb2 < __a2);	__wd2 += (__wc2 < __b2);	__wb2 += __a2;	__wc2 += __b2 + (__wb2 < __a2);	__wd2 += (__wc2 < __b2);\
		__wb3 += __a3;	__wc3 += __b3 + (__wb3 < __a3);	__wd3 += (__wc3 < __b3);	__wb3 += __a3;	__wc3 += __b3 + (__wb3 < __a3);	__wd3 += (__wc3 < __b3);\
		__wb4 += __a4;	__wc4 += __b4 + (__wb4 < __a4);	__wd4 += (__wc4 < __b4);	__wb4 += __a4;	__wc4 += __b4 + (__wb4 < __a4);	__wd4 += (__wc4 < __b4);\
		__wb5 += __a5;	__wc5 += __b5 + (__wb5 < __a5);	__wd5 += (__wc5 < __b5);	__wb5 += __a5;	__wc5 += __b5 + (__wb5 < __a5);	__wd5 += (__wc5 < __b5);\
		__wb6 += __a6;	__wc6 += __b6 + (__wb6 < __a6);	__wd6 += (__wc6 < __b6);	__wb6 += __a6;	__wc6 += __b6 + (__wb6 < __a6);	__wd6 += (__wc6 < __b6);\
		__wb7 += __a7;	__wc7 += __b7 + (__wb7 < __a7);	__wd7 += (__wc7 < __b7);	__wb7 += __a7;	__wc7 += __b7 + (__wb7 < __a7);	__wd7 += (__wc7 < __b7);\
		\
		/* Now store the high 128 bits of the result in __hi: */\
		__hi0.d0 =  __wc0;	__hi0.d1 = __wd0;\
		__hi1.d0 =  __wc1;	__hi1.d1 = __wd1;\
		__hi2.d0 =  __wc2;	__hi2.d1 = __wd2;\
		__hi3.d0 =  __wc3;	__hi3.d1 = __wd3;\
		__hi4.d0 =  __wc4;	__hi4.d1 = __wd4;\
		__hi5.d0 =  __wc5;	__hi5.d1 = __wd5;\
		__hi6.d0 =  __wc6;	__hi6.d1 = __wd6;\
		__hi7.d0 =  __wc7;	__hi7.d1 = __wd7;\
		\
		/* MULL128(lo,qinv,lo) IS HERE: wc, wd now free, so use as 64-bit temporaries: */\
		\
		__a0 = __wa0*__qinv0.d1 + __qinv0.d0*__wb0;\
		__a1 = __wa1*__qinv1.d1 + __qinv1.d0*__wb1;\
		__a2 = __wa2*__qinv2.d1 + __qinv2.d0*__wb2;\
		__a3 = __wa3*__qinv3.d1 + __qinv3.d0*__wb3;\
		__a4 = __wa4*__qinv4.d1 + __qinv4.d0*__wb4;\
		__a5 = __wa5*__qinv5.d1 + __qinv5.d0*__wb5;\
		__a6 = __wa6*__qinv6.d1 + __qinv6.d0*__wb6;\
		__a7 = __wa7*__qinv7.d1 + __qinv7.d0*__wb7;\
		\
		MUL_LOHI64(__wa0,__qinv0.d0,&__wc0,&__wd0);\
		MUL_LOHI64(__wa1,__qinv1.d0,&__wc1,&__wd1);\
		MUL_LOHI64(__wa2,__qinv2.d0,&__wc2,&__wd2);\
		MUL_LOHI64(__wa3,__qinv3.d0,&__wc3,&__wd3);\
		MUL_LOHI64(__wa4,__qinv4.d0,&__wc4,&__wd4);\
		MUL_LOHI64(__wa5,__qinv5.d0,&__wc5,&__wd5);\
		MUL_LOHI64(__wa6,__qinv6.d0,&__wc6,&__wd6);\
		MUL_LOHI64(__wa7,__qinv7.d0,&__wc7,&__wd7);\
		\
		__wa0 = __wc0;\
		__wa1 = __wc1;\
		__wa2 = __wc2;\
		__wa3 = __wc3;\
		__wa4 = __wc4;\
		__wa5 = __wc5;\
		__wa6 = __wc6;\
		__wa7 = __wc7;\
		\
		__wb0 = __wd0 + __a0;\
		__wb1 = __wd1 + __a1;\
		__wb2 = __wd2 + __a2;\
		__wb3 = __wd3 + __a3;\
		__wb4 = __wd4 + __a4;\
		__wb5 = __wd5 + __a5;\
		__wb6 = __wd6 + __a6;\
		__wb7 = __wd7 + __a7;\
		\
		/* MULH128(q,lo,lo) IS HERE: store result in __lo: */\
		\
		MULH64(__wa0,__q0.d0, __t0);\
		MULH64(__wa1,__q1.d0, __t1);\
		MULH64(__wa2,__q2.d0, __t2);\
		MULH64(__wa3,__q3.d0, __t3);\
		MULH64(__wa4,__q4.d0, __t4);\
		MULH64(__wa5,__q5.d0, __t5);\
		MULH64(__wa6,__q6.d0, __t6);\
		MULH64(__wa7,__q7.d0, __t7);\
		\
		MUL_LOHI64(__wa0,__q0.d1, &__a0, &__b0);\
		MUL_LOHI64(__wa1,__q1.d1, &__a1, &__b1);\
		MUL_LOHI64(__wa2,__q2.d1, &__a2, &__b2);\
		MUL_LOHI64(__wa3,__q3.d1, &__a3, &__b3);\
		MUL_LOHI64(__wa4,__q4.d1, &__a4, &__b4);\
		MUL_LOHI64(__wa5,__q5.d1, &__a5, &__b5);\
		MUL_LOHI64(__wa6,__q6.d1, &__a6, &__b6);\
		MUL_LOHI64(__wa7,__q7.d1, &__a7, &__b7);\
		\
		MUL_LOHI64(__q0.d0,__wb0,&__wc0,&__wd0);\
		MUL_LOHI64(__q1.d0,__wb1,&__wc1,&__wd1);\
		MUL_LOHI64(__q2.d0,__wb2,&__wc2,&__wd2);\
		MUL_LOHI64(__q3.d0,__wb3,&__wc3,&__wd3);\
		MUL_LOHI64(__q4.d0,__wb4,&__wc4,&__wd4);\
		MUL_LOHI64(__q5.d0,__wb5,&__wc5,&__wd5);\
		MUL_LOHI64(__q6.d0,__wb6,&__wc6,&__wd6);\
		MUL_LOHI64(__q7.d0,__wb7,&__wc7,&__wd7);\
		\
		MUL_LOHI64(__wb0,__q0.d1,&__lo0.d0,&__lo0.d1);\
		MUL_LOHI64(__wb1,__q1.d1,&__lo1.d0,&__lo1.d1);\
		MUL_LOHI64(__wb2,__q2.d1,&__lo2.d0,&__lo2.d1);\
		MUL_LOHI64(__wb3,__q3.d1,&__lo3.d0,&__lo3.d1);\
		MUL_LOHI64(__wb4,__q4.d1,&__lo4.d0,&__lo4.d1);\
		MUL_LOHI64(__wb5,__q5.d1,&__lo5.d0,&__lo5.d1);\
		MUL_LOHI64(__wb6,__q6.d1,&__lo6.d0,&__lo6.d1);\
		MUL_LOHI64(__wb7,__q7.d1,&__lo7.d0,&__lo7.d1);\
		\
		/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
		__a0 += __wc0;\
		__a1 += __wc1;\
		__a2 += __wc2;\
		__a3 += __wc3;\
		__a4 += __wc4;\
		__a5 += __wc5;\
		__a6 += __wc6;\
		__a7 += __wc7;\
		\
		__b0 += __wd0 + (__a0 < __wc0);\
		__b1 += __wd1 + (__a1 < __wc1);\
		__b2 += __wd2 + (__a2 < __wc2);\
		__b3 += __wd3 + (__a3 < __wc3);\
		__b4 += __wd4 + (__a4 < __wc4);\
		__b5 += __wd5 + (__a5 < __wc5);\
		__b6 += __wd6 + (__a6 < __wc6);\
		__b7 += __wd7 + (__a7 < __wc7);\
		\
		/* Now add [w1,w2,w3] + [a,b,0]: */\
		__t0 += __a0;\
		__t1 += __a1;\
		__t2 += __a2;\
		__t3 += __a3;\
		__t4 += __a4;\
		__t5 += __a5;\
		__t6 += __a6;\
		__t7 += __a7;\
		\
		__cy0 = (__t0 < __a0);\
		__cy1 = (__t1 < __a1);\
		__cy2 = (__t2 < __a2);\
		__cy3 = (__t3 < __a3);\
		__cy4 = (__t4 < __a4);\
		__cy5 = (__t5 < __a5);\
		__cy6 = (__t6 < __a6);\
		__cy7 = (__t7 < __a7);\
		\
		__lo0.d0 += __cy0;\
		__lo1.d0 += __cy1;\
		__lo2.d0 += __cy2;\
		__lo3.d0 += __cy3;\
		__lo4.d0 += __cy4;\
		__lo5.d0 += __cy5;\
		__lo6.d0 += __cy6;\
		__lo7.d0 += __cy7;\
		\
		__lo0.d1 += (__lo0.d0 < __cy0);\
		__lo1.d1 += (__lo1.d0 < __cy1);\
		__lo2.d1 += (__lo2.d0 < __cy2);\
		__lo3.d1 += (__lo3.d0 < __cy3);\
		__lo4.d1 += (__lo4.d0 < __cy4);\
		__lo5.d1 += (__lo5.d0 < __cy5);\
		__lo6.d1 += (__lo6.d0 < __cy6);\
		__lo7.d1 += (__lo7.d0 < __cy7);\
		\
		__lo0.d0 += __b0;\
		__lo1.d0 += __b1;\
		__lo2.d0 += __b2;\
		__lo3.d0 += __b3;\
		__lo4.d0 += __b4;\
		__lo5.d0 += __b5;\
		__lo6.d0 += __b6;\
		__lo7.d0 += __b7;\
		\
		__lo0.d1 += (__lo0.d0 < __b0);\
		__lo1.d1 += (__lo1.d0 < __b1);\
		__lo2.d1 += (__lo2.d0 < __b2);\
		__lo3.d1 += (__lo3.d0 < __b3);\
		__lo4.d1 += (__lo4.d0 < __b4);\
		__lo5.d1 += (__lo5.d0 < __b5);\
		__lo6.d1 += (__lo6.d0 < __b6);\
		__lo7.d1 += (__lo7.d0 < __b7);\
    }

#else

    #define MULH128(__x, __y, __hi)\
    {\
		uint64 __w1,__w2,__w3,__a,__b,__c,__d;\
		\
		MULH64(__x.d0,__y.d0,       __w1);\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );\
		MUL_LOHI64(__y.d0,__x.d1, __c , __d );\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __b); /* Overflow into word3 is checked here. */\
		__w1 += __c;\
		__w2 += __d + (__w1 < __c); /* Overflow into word2 is checked here. */\
		__w3 +=       (__w2 < __d); /* Overflow into word3 is checked here. */\
		__hi.d0 =  __w2;	__hi.d1 = __w3;\
    }

    /* 4-operand-pipelined version: */
    #define MULH128_q4(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3)\
    {\
		uint64	__a0,__a1,__a2,__a3,\
				__b0,__b1,__b2,__b3,\
				__c0,__c1,__c2,__c3,\
				__d0,__d1,__d2,__d3,\
				__t0,__t1,__t2,__t3;\
		\
		MULH64(__x0.d0,__y0.d0, __t0);\
		MULH64(__x1.d0,__y1.d0, __t1);\
		MULH64(__x2.d0,__y2.d0, __t2);\
		MULH64(__x3.d0,__y3.d0, __t3);\
		\
		MUL_LOHI64(__x0.d0,__y0.d1, __a0, __b0);\
		MUL_LOHI64(__x1.d0,__y1.d1, __a1, __b1);\
		MUL_LOHI64(__x2.d0,__y2.d1, __a2, __b2);\
		MUL_LOHI64(__x3.d0,__y3.d1, __a3, __b3);\
		\
		MUL_LOHI64(__y0.d0,__x0.d1, __c0, __d0);\
		MUL_LOHI64(__y1.d0,__x1.d1, __c1, __d1);\
		MUL_LOHI64(__y2.d0,__x2.d1, __c2, __d2);\
		MUL_LOHI64(__y3.d0,__x3.d1, __c3, __d3);\
		\
		MUL_LOHI64(__x0.d1,__y0.d1, __hi0.d0, __hi0.d1);\
		MUL_LOHI64(__x1.d1,__y1.d1, __hi1.d0, __hi1.d1);\
		MUL_LOHI64(__x2.d1,__y2.d1, __hi2.d0, __hi2.d1);\
		MUL_LOHI64(__x3.d1,__y3.d1, __hi3.d0, __hi3.d1);\
		\
		__t0 += __a0;\
		__t1 += __a1;\
		__t2 += __a2;\
		__t3 += __a3;\
		\
		__hi0.d0 += __b0 + (__t0 < __a0);\
		__hi1.d0 += __b1 + (__t1 < __a1);\
		__hi2.d0 += __b2 + (__t2 < __a2);\
		__hi3.d0 += __b3 + (__t3 < __a3);\
		\
		__hi0.d1 +=        (__hi0.d0 < __b0);\
		__hi1.d1 +=        (__hi1.d0 < __b1);\
		__hi2.d1 +=        (__hi2.d0 < __b2);\
		__hi3.d1 +=        (__hi3.d0 < __b3);\
		\
		__t0 += __c0;\
		__t1 += __c1;\
		__t2 += __c2;\
		__t3 += __c3;\
		\
		__hi0.d0 += __d0 + (__t0 < __c0);\
		__hi1.d0 += __d1 + (__t1 < __c1);\
		__hi2.d0 += __d2 + (__t2 < __c2);\
		__hi3.d0 += __d3 + (__t3 < __c3);\
		\
		__hi0.d1 +=        (__hi0.d0 < __d0);\
		__hi1.d1 +=        (__hi1.d0 < __d1);\
		__hi2.d1 +=        (__hi2.d0 < __d2);\
		__hi3.d1 +=        (__hi3.d0 < __d3);\
    }

/************* TODO: OK with y == hi, but not with x == hi! *************/
    /* 8-operand-pipelined version: */
    #define MULH128_q8(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3\
    , __x4, __y4, __hi4\
    , __x5, __y5, __hi5\
    , __x6, __y6, __hi6\
    , __x7, __y7, __hi7)\
    {\
		uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
				__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
				__c0,__c1,__c2,__c3,__c4,__c5,__c6,__c7,\
				__d0,__d1,__d2,__d3,__d4,__d5,__d6,__d7,\
				__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7;\
		\
		MULH64(__x0.d0,__y0.d0, __t0);\
		MULH64(__x1.d0,__y1.d0, __t1);\
		MULH64(__x2.d0,__y2.d0, __t2);\
		MULH64(__x3.d0,__y3.d0, __t3);\
		MULH64(__x4.d0,__y4.d0, __t4);\
		MULH64(__x5.d0,__y5.d0, __t5);\
		MULH64(__x6.d0,__y6.d0, __t6);\
		MULH64(__x7.d0,__y7.d0, __t7);\
		\
		MUL_LOHI64(__x0.d0,__y0.d1, __a0, __b0);\
		MUL_LOHI64(__x1.d0,__y1.d1, __a1, __b1);\
		MUL_LOHI64(__x2.d0,__y2.d1, __a2, __b2);\
		MUL_LOHI64(__x3.d0,__y3.d1, __a3, __b3);\
		MUL_LOHI64(__x4.d0,__y4.d1, __a4, __b4);\
		MUL_LOHI64(__x5.d0,__y5.d1, __a5, __b5);\
		MUL_LOHI64(__x6.d0,__y6.d1, __a6, __b6);\
		MUL_LOHI64(__x7.d0,__y7.d1, __a7, __b7);\
		\
		MUL_LOHI64(__y0.d0,__x0.d1, __c0, __d0);\
		MUL_LOHI64(__y1.d0,__x1.d1, __c1, __d1);\
		MUL_LOHI64(__y2.d0,__x2.d1, __c2, __d2);\
		MUL_LOHI64(__y3.d0,__x3.d1, __c3, __d3);\
		MUL_LOHI64(__y4.d0,__x4.d1, __c4, __d4);\
		MUL_LOHI64(__y5.d0,__x5.d1, __c5, __d5);\
		MUL_LOHI64(__y6.d0,__x6.d1, __c6, __d6);\
		MUL_LOHI64(__y7.d0,__x7.d1, __c7, __d7);\
		\
		MUL_LOHI64(__x0.d1,__y0.d1, __hi0.d0, __hi0.d1);\
		MUL_LOHI64(__x1.d1,__y1.d1, __hi1.d0, __hi1.d1);\
		MUL_LOHI64(__x2.d1,__y2.d1, __hi2.d0, __hi2.d1);\
		MUL_LOHI64(__x3.d1,__y3.d1, __hi3.d0, __hi3.d1);\
		MUL_LOHI64(__x4.d1,__y4.d1, __hi4.d0, __hi4.d1);\
		MUL_LOHI64(__x5.d1,__y5.d1, __hi5.d0, __hi5.d1);\
		MUL_LOHI64(__x6.d1,__y6.d1, __hi6.d0, __hi6.d1);\
		MUL_LOHI64(__x7.d1,__y7.d1, __hi7.d0, __hi7.d1);\
		\
		__t0 += __a0;\
		__t1 += __a1;\
		__t2 += __a2;\
		__t3 += __a3;\
		__t4 += __a4;\
		__t5 += __a5;\
		__t6 += __a6;\
		__t7 += __a7;\
		\
		__hi0.d0 += __b0 + (__t0 < __a0);\
		__hi1.d0 += __b1 + (__t1 < __a1);\
		__hi2.d0 += __b2 + (__t2 < __a2);\
		__hi3.d0 += __b3 + (__t3 < __a3);\
		__hi4.d0 += __b4 + (__t4 < __a4);\
		__hi5.d0 += __b5 + (__t5 < __a5);\
		__hi6.d0 += __b6 + (__t6 < __a6);\
		__hi7.d0 += __b7 + (__t7 < __a7);\
		\
		__hi0.d1 +=        (__hi0.d0 < __b0);\
		__hi1.d1 +=        (__hi1.d0 < __b1);\
		__hi2.d1 +=        (__hi2.d0 < __b2);\
		__hi3.d1 +=        (__hi3.d0 < __b3);\
		__hi4.d1 +=        (__hi4.d0 < __b4);\
		__hi5.d1 +=        (__hi5.d0 < __b5);\
		__hi6.d1 +=        (__hi6.d0 < __b6);\
		__hi7.d1 +=        (__hi7.d0 < __b7);\
		\
		__t0 += __c0;\
		__t1 += __c1;\
		__t2 += __c2;\
		__t3 += __c3;\
		__t4 += __c4;\
		__t5 += __c5;\
		__t6 += __c6;\
		__t7 += __c7;\
		\
		__hi0.d0 += __d0 + (__t0 < __c0);\
		__hi1.d0 += __d1 + (__t1 < __c1);\
		__hi2.d0 += __d2 + (__t2 < __c2);\
		__hi3.d0 += __d3 + (__t3 < __c3);\
		__hi4.d0 += __d4 + (__t4 < __c4);\
		__hi5.d0 += __d5 + (__t5 < __c5);\
		__hi6.d0 += __d6 + (__t6 < __c6);\
		__hi7.d0 += __d7 + (__t7 < __c7);\
		\
		__hi0.d1 +=        (__hi0.d0 < __d0);\
		__hi1.d1 +=        (__hi1.d0 < __d1);\
		__hi2.d1 +=        (__hi2.d0 < __d2);\
		__hi3.d1 +=        (__hi3.d0 < __d3);\
		__hi4.d1 +=        (__hi4.d0 < __d4);\
		__hi5.d1 +=        (__hi5.d0 < __d5);\
		__hi6.d1 +=        (__hi6.d0 < __d6);\
		__hi7.d1 +=        (__hi7.d0 < __d7);\
    }

/* Routine to perform fused version of key 3-operation sequence in 128-bit modmul:

	SQR_LOHI128(x,lo,hi);
	MULL128(lo,qinv,lo);
	MULH128(q,lo,lo);
*/
    /* 4-operand-pipelined version: */
    #define THREE_OP128_q4(\
      __x0, __qinv0, __q0, __lo0, __hi0\
    , __x1, __qinv1, __q1, __lo1, __hi1\
    , __x2, __qinv2, __q2, __lo2, __hi2\
    , __x3, __qinv3, __q3, __lo3, __hi3)\
    {\
		uint64	__a0,__b0,__t0,__wa0,__wb0,__wc0,__wd0,__cy0;\
		uint64	__a1,__b1,__t1,__wa1,__wb1,__wc1,__wd1,__cy1;\
		uint64	__a2,__b2,__t2,__wa2,__wb2,__wc2,__wd2,__cy2;\
		uint64	__a3,__b3,__t3,__wa3,__wb3,__wc3,__wd3,__cy3;\
		\
		/* SQR_LOHI128(x,lo,hi) IS HERE: */\
		SQR_LOHI64(__x0.d0,         __wa0, __wb0);\
		SQR_LOHI64(__x1.d0,         __wa1, __wb1);\
		SQR_LOHI64(__x2.d0,         __wa2, __wb2);\
		SQR_LOHI64(__x3.d0,         __wa3, __wb3);\
		/* On IA64, can add the __wb's in with the __a's here: */\
		MUL_LOHI64(__x0.d0,__x0.d1, __a0 , __b0);\
		MUL_LOHI64(__x1.d0,__x1.d1, __a1 , __b1);\
		MUL_LOHI64(__x2.d0,__x2.d1, __a2 , __b2);\
		MUL_LOHI64(__x3.d0,__x3.d1, __a3 , __b3);\
		/* On IA64, can add the __b's in with the __wc's here: */\
		SQR_LOHI64(__x0.d1,         __wc0, __wd0);\
		SQR_LOHI64(__x1.d1,         __wc1, __wd1);\
		SQR_LOHI64(__x2.d1,         __wc2, __wd2);\
		SQR_LOHI64(__x3.d1,         __wc3, __wd3);\
		\
		/* Need to add 2*a*b, so add a*b twice: */\
		__wb0 += __a0;	__wc0 += __b0 + (__wb0 < __a0);	__wd0 += (__wc0 < __b0);	__wb0 += __a0;	__wc0 += __b0 + (__wb0 < __a0);	__wd0 += (__wc0 < __b0);\
		__wb1 += __a1;	__wc1 += __b1 + (__wb1 < __a1);	__wd1 += (__wc1 < __b1);	__wb1 += __a1;	__wc1 += __b1 + (__wb1 < __a1);	__wd1 += (__wc1 < __b1);\
		__wb2 += __a2;	__wc2 += __b2 + (__wb2 < __a2);	__wd2 += (__wc2 < __b2);	__wb2 += __a2;	__wc2 += __b2 + (__wb2 < __a2);	__wd2 += (__wc2 < __b2);\
		__wb3 += __a3;	__wc3 += __b3 + (__wb3 < __a3);	__wd3 += (__wc3 < __b3);	__wb3 += __a3;	__wc3 += __b3 + (__wb3 < __a3);	__wd3 += (__wc3 < __b3);\
		\
		/* Now store the high 128 bits of the result in __hi: */\
		__hi0.d0 =  __wc0;	__hi0.d1 = __wd0;\
		__hi1.d0 =  __wc1;	__hi1.d1 = __wd1;\
		__hi2.d0 =  __wc2;	__hi2.d1 = __wd2;\
		__hi3.d0 =  __wc3;	__hi3.d1 = __wd3;\
		\
		/* MULL128(lo,qinv,lo) IS HERE: wc, wd now free, so use as 64-bit temporaries: */\
		\
		__a0 = __wa0*__qinv0.d1 + __qinv0.d0*__wb0;\
		__a1 = __wa1*__qinv1.d1 + __qinv1.d0*__wb1;\
		__a2 = __wa2*__qinv2.d1 + __qinv2.d0*__wb2;\
		__a3 = __wa3*__qinv3.d1 + __qinv3.d0*__wb3;\
		\
		MUL_LOHI64(__wa0,__qinv0.d0, __wc0, __wd0);\
		MUL_LOHI64(__wa1,__qinv1.d0, __wc1, __wd1);\
		MUL_LOHI64(__wa2,__qinv2.d0, __wc2, __wd2);\
		MUL_LOHI64(__wa3,__qinv3.d0, __wc3, __wd3);\
		\
		__wa0 = __wc0;\
		__wa1 = __wc1;\
		__wa2 = __wc2;\
		__wa3 = __wc3;\
		\
		__wb0 = __wd0 + __a0;\
		__wb1 = __wd1 + __a1;\
		__wb2 = __wd2 + __a2;\
		__wb3 = __wd3 + __a3;\
		\
		/* MULH128(q,lo,lo) IS HERE: store result in __lo: */\
		\
		MULH64(__wa0,__q0.d0, __t0);\
		MULH64(__wa1,__q1.d0, __t1);\
		MULH64(__wa2,__q2.d0, __t2);\
		MULH64(__wa3,__q3.d0, __t3);\
		\
		MUL_LOHI64(__wa0,__q0.d1,  __a0,  __b0);\
		MUL_LOHI64(__wa1,__q1.d1,  __a1,  __b1);\
		MUL_LOHI64(__wa2,__q2.d1,  __a2,  __b2);\
		MUL_LOHI64(__wa3,__q3.d1,  __a3,  __b3);\
		\
		MUL_LOHI64(__q0.d0,__wb0, __wc0, __wd0);\
		MUL_LOHI64(__q1.d0,__wb1, __wc1, __wd1);\
		MUL_LOHI64(__q2.d0,__wb2, __wc2, __wd2);\
		MUL_LOHI64(__q3.d0,__wb3, __wc3, __wd3);\
		\
		MUL_LOHI64(__wb0,__q0.d1, __lo0.d0, __lo0.d1);\
		MUL_LOHI64(__wb1,__q1.d1, __lo1.d0, __lo1.d1);\
		MUL_LOHI64(__wb2,__q2.d1, __lo2.d0, __lo2.d1);\
		MUL_LOHI64(__wb3,__q3.d1, __lo3.d0, __lo3.d1);\
		\
		/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
		__a0 += __wc0;\
		__a1 += __wc1;\
		__a2 += __wc2;\
		__a3 += __wc3;\
		\
		__b0 += __wd0 + (__a0 < __wc0);\
		__b1 += __wd1 + (__a1 < __wc1);\
		__b2 += __wd2 + (__a2 < __wc2);\
		__b3 += __wd3 + (__a3 < __wc3);\
		\
		/* Now add [w1,w2,w3] + [a,b,0]: */\
		__t0 += __a0;\
		__t1 += __a1;\
		__t2 += __a2;\
		__t3 += __a3;\
		\
		__cy0 = (__t0 < __a0);\
		__cy1 = (__t1 < __a1);\
		__cy2 = (__t2 < __a2);\
		__cy3 = (__t3 < __a3);\
		\
		__lo0.d0 += __cy0;\
		__lo1.d0 += __cy1;\
		__lo2.d0 += __cy2;\
		__lo3.d0 += __cy3;\
		\
		__lo0.d1 += (__lo0.d0 < __cy0);\
		__lo1.d1 += (__lo1.d0 < __cy1);\
		__lo2.d1 += (__lo2.d0 < __cy2);\
		__lo3.d1 += (__lo3.d0 < __cy3);\
		\
		__lo0.d0 += __b0;\
		__lo1.d0 += __b1;\
		__lo2.d0 += __b2;\
		__lo3.d0 += __b3;\
		\
		__lo0.d1 += (__lo0.d0 < __b0);\
		__lo1.d1 += (__lo1.d0 < __b1);\
		__lo2.d1 += (__lo2.d0 < __b2);\
		__lo3.d1 += (__lo3.d0 < __b3);\
    }

    /* 8-operand-pipelined version: */
    #define THREE_OP128_q8(\
      __x0, __qinv0, __q0, __lo0, __hi0\
    , __x1, __qinv1, __q1, __lo1, __hi1\
    , __x2, __qinv2, __q2, __lo2, __hi2\
    , __x3, __qinv3, __q3, __lo3, __hi3\
    , __x4, __qinv4, __q4, __lo4, __hi4\
    , __x5, __qinv5, __q5, __lo5, __hi5\
    , __x6, __qinv6, __q6, __lo6, __hi6\
    , __x7, __qinv7, __q7, __lo7, __hi7)\
    {\
		uint64	__a0,__b0,__t0,__wa0,__wb0,__wc0,__wd0,__cy0;\
		uint64	__a1,__b1,__t1,__wa1,__wb1,__wc1,__wd1,__cy1;\
		uint64	__a2,__b2,__t2,__wa2,__wb2,__wc2,__wd2,__cy2;\
		uint64	__a3,__b3,__t3,__wa3,__wb3,__wc3,__wd3,__cy3;\
		uint64	__a4,__b4,__t4,__wa4,__wb4,__wc4,__wd4,__cy4;\
		uint64	__a5,__b5,__t5,__wa5,__wb5,__wc5,__wd5,__cy5;\
		uint64	__a6,__b6,__t6,__wa6,__wb6,__wc6,__wd6,__cy6;\
		uint64	__a7,__b7,__t7,__wa7,__wb7,__wc7,__wd7,__cy7;\
		\
		/* SQR_LOHI128(x,lo,hi) IS HERE: */\
		SQR_LOHI64(__x0.d0,        __wa0, __wb0);\
		SQR_LOHI64(__x1.d0,        __wa1, __wb1);\
		SQR_LOHI64(__x2.d0,        __wa2, __wb2);\
		SQR_LOHI64(__x3.d0,        __wa3, __wb3);\
		SQR_LOHI64(__x4.d0,        __wa4, __wb4);\
		SQR_LOHI64(__x5.d0,        __wa5, __wb5);\
		SQR_LOHI64(__x6.d0,        __wa6, __wb6);\
		SQR_LOHI64(__x7.d0,        __wa7, __wb7);\
		/* On IA64, can add the __wb's in with the __a's here: */\
		MUL_LOHI64(__x0.d0,__x0.d1, __a0 , __b0);\
		MUL_LOHI64(__x1.d0,__x1.d1, __a1 , __b1);\
		MUL_LOHI64(__x2.d0,__x2.d1, __a2 , __b2);\
		MUL_LOHI64(__x3.d0,__x3.d1, __a3 , __b3);\
		MUL_LOHI64(__x4.d0,__x4.d1, __a4 , __b4);\
		MUL_LOHI64(__x5.d0,__x5.d1, __a5 , __b5);\
		MUL_LOHI64(__x6.d0,__x6.d1, __a6 , __b6);\
		MUL_LOHI64(__x7.d0,__x7.d1, __a7 , __b7);\
		/* On IA64, can add the __b's in with the __wc's here: */\
		SQR_LOHI64(__x0.d1,        __wc0, __wd0);\
		SQR_LOHI64(__x1.d1,        __wc1, __wd1);\
		SQR_LOHI64(__x2.d1,        __wc2, __wd2);\
		SQR_LOHI64(__x3.d1,        __wc3, __wd3);\
		SQR_LOHI64(__x4.d1,        __wc4, __wd4);\
		SQR_LOHI64(__x5.d1,        __wc5, __wd5);\
		SQR_LOHI64(__x6.d1,        __wc6, __wd6);\
		SQR_LOHI64(__x7.d1,        __wc7, __wd7);\
		\
		/* Need to add 2*a*b, so add a*b twice: */\
		__wb0 += __a0;	__wc0 += __b0 + (__wb0 < __a0);	__wd0 += (__wc0 < __b0);	__wb0 += __a0;	__wc0 += __b0 + (__wb0 < __a0);	__wd0 += (__wc0 < __b0);\
		__wb1 += __a1;	__wc1 += __b1 + (__wb1 < __a1);	__wd1 += (__wc1 < __b1);	__wb1 += __a1;	__wc1 += __b1 + (__wb1 < __a1);	__wd1 += (__wc1 < __b1);\
		__wb2 += __a2;	__wc2 += __b2 + (__wb2 < __a2);	__wd2 += (__wc2 < __b2);	__wb2 += __a2;	__wc2 += __b2 + (__wb2 < __a2);	__wd2 += (__wc2 < __b2);\
		__wb3 += __a3;	__wc3 += __b3 + (__wb3 < __a3);	__wd3 += (__wc3 < __b3);	__wb3 += __a3;	__wc3 += __b3 + (__wb3 < __a3);	__wd3 += (__wc3 < __b3);\
		__wb4 += __a4;	__wc4 += __b4 + (__wb4 < __a4);	__wd4 += (__wc4 < __b4);	__wb4 += __a4;	__wc4 += __b4 + (__wb4 < __a4);	__wd4 += (__wc4 < __b4);\
		__wb5 += __a5;	__wc5 += __b5 + (__wb5 < __a5);	__wd5 += (__wc5 < __b5);	__wb5 += __a5;	__wc5 += __b5 + (__wb5 < __a5);	__wd5 += (__wc5 < __b5);\
		__wb6 += __a6;	__wc6 += __b6 + (__wb6 < __a6);	__wd6 += (__wc6 < __b6);	__wb6 += __a6;	__wc6 += __b6 + (__wb6 < __a6);	__wd6 += (__wc6 < __b6);\
		__wb7 += __a7;	__wc7 += __b7 + (__wb7 < __a7);	__wd7 += (__wc7 < __b7);	__wb7 += __a7;	__wc7 += __b7 + (__wb7 < __a7);	__wd7 += (__wc7 < __b7);\
		\
		/* Now store the high 128 bits of the result in __hi: */\
		__hi0.d0 =  __wc0;	__hi0.d1 = __wd0;\
		__hi1.d0 =  __wc1;	__hi1.d1 = __wd1;\
		__hi2.d0 =  __wc2;	__hi2.d1 = __wd2;\
		__hi3.d0 =  __wc3;	__hi3.d1 = __wd3;\
		__hi4.d0 =  __wc4;	__hi4.d1 = __wd4;\
		__hi5.d0 =  __wc5;	__hi5.d1 = __wd5;\
		__hi6.d0 =  __wc6;	__hi6.d1 = __wd6;\
		__hi7.d0 =  __wc7;	__hi7.d1 = __wd7;\
		\
		/* MULL128(lo,qinv,lo) IS HERE: wc, wd now free, so use as 64-bit temporaries: */\
		\
		__a0 = __wa0*__qinv0.d1 + __qinv0.d0*__wb0;\
		__a1 = __wa1*__qinv1.d1 + __qinv1.d0*__wb1;\
		__a2 = __wa2*__qinv2.d1 + __qinv2.d0*__wb2;\
		__a3 = __wa3*__qinv3.d1 + __qinv3.d0*__wb3;\
		__a4 = __wa4*__qinv4.d1 + __qinv4.d0*__wb4;\
		__a5 = __wa5*__qinv5.d1 + __qinv5.d0*__wb5;\
		__a6 = __wa6*__qinv6.d1 + __qinv6.d0*__wb6;\
		__a7 = __wa7*__qinv7.d1 + __qinv7.d0*__wb7;\
		\
		MUL_LOHI64(__wa0,__qinv0.d0, __wc0, __wd0);\
		MUL_LOHI64(__wa1,__qinv1.d0, __wc1, __wd1);\
		MUL_LOHI64(__wa2,__qinv2.d0, __wc2, __wd2);\
		MUL_LOHI64(__wa3,__qinv3.d0, __wc3, __wd3);\
		MUL_LOHI64(__wa4,__qinv4.d0, __wc4, __wd4);\
		MUL_LOHI64(__wa5,__qinv5.d0, __wc5, __wd5);\
		MUL_LOHI64(__wa6,__qinv6.d0, __wc6, __wd6);\
		MUL_LOHI64(__wa7,__qinv7.d0, __wc7, __wd7);\
		\
		__wa0 = __wc0;\
		__wa1 = __wc1;\
		__wa2 = __wc2;\
		__wa3 = __wc3;\
		__wa4 = __wc4;\
		__wa5 = __wc5;\
		__wa6 = __wc6;\
		__wa7 = __wc7;\
		\
		__wb0 = __wd0 + __a0;\
		__wb1 = __wd1 + __a1;\
		__wb2 = __wd2 + __a2;\
		__wb3 = __wd3 + __a3;\
		__wb4 = __wd4 + __a4;\
		__wb5 = __wd5 + __a5;\
		__wb6 = __wd6 + __a6;\
		__wb7 = __wd7 + __a7;\
		\
		/* MULH128(q,lo,lo) IS HERE: store result in __lo: */\
		\
		MULH64(__wa0,__q0.d0, __t0);\
		MULH64(__wa1,__q1.d0, __t1);\
		MULH64(__wa2,__q2.d0, __t2);\
		MULH64(__wa3,__q3.d0, __t3);\
		MULH64(__wa4,__q4.d0, __t4);\
		MULH64(__wa5,__q5.d0, __t5);\
		MULH64(__wa6,__q6.d0, __t6);\
		MULH64(__wa7,__q7.d0, __t7);\
		\
		MUL_LOHI64(__wa0,__q0.d1, __a0, __b0);\
		MUL_LOHI64(__wa1,__q1.d1, __a1, __b1);\
		MUL_LOHI64(__wa2,__q2.d1, __a2, __b2);\
		MUL_LOHI64(__wa3,__q3.d1, __a3, __b3);\
		MUL_LOHI64(__wa4,__q4.d1, __a4, __b4);\
		MUL_LOHI64(__wa5,__q5.d1, __a5, __b5);\
		MUL_LOHI64(__wa6,__q6.d1, __a6, __b6);\
		MUL_LOHI64(__wa7,__q7.d1, __a7, __b7);\
		\
		MUL_LOHI64(__q0.d0,__wb0, __wc0, __wd0);\
		MUL_LOHI64(__q1.d0,__wb1, __wc1, __wd1);\
		MUL_LOHI64(__q2.d0,__wb2, __wc2, __wd2);\
		MUL_LOHI64(__q3.d0,__wb3, __wc3, __wd3);\
		MUL_LOHI64(__q4.d0,__wb4, __wc4, __wd4);\
		MUL_LOHI64(__q5.d0,__wb5, __wc5, __wd5);\
		MUL_LOHI64(__q6.d0,__wb6, __wc6, __wd6);\
		MUL_LOHI64(__q7.d0,__wb7, __wc7, __wd7);\
		\
		MUL_LOHI64(__wb0,__q0.d1, __lo0.d0, __lo0.d1);\
		MUL_LOHI64(__wb1,__q1.d1, __lo1.d0, __lo1.d1);\
		MUL_LOHI64(__wb2,__q2.d1, __lo2.d0, __lo2.d1);\
		MUL_LOHI64(__wb3,__q3.d1, __lo3.d0, __lo3.d1);\
		MUL_LOHI64(__wb4,__q4.d1, __lo4.d0, __lo4.d1);\
		MUL_LOHI64(__wb5,__q5.d1, __lo5.d0, __lo5.d1);\
		MUL_LOHI64(__wb6,__q6.d1, __lo6.d0, __lo6.d1);\
		MUL_LOHI64(__wb7,__q7.d1, __lo7.d0, __lo7.d1);\
		\
		/* First add [a,b] + [c,d] : since b and d <= 2^64 - 2, can add carryout of a+c sans ripple-carry check: */\
		__a0 += __wc0;\
		__a1 += __wc1;\
		__a2 += __wc2;\
		__a3 += __wc3;\
		__a4 += __wc4;\
		__a5 += __wc5;\
		__a6 += __wc6;\
		__a7 += __wc7;\
		\
		__b0 += __wd0 + (__a0 < __wc0);\
		__b1 += __wd1 + (__a1 < __wc1);\
		__b2 += __wd2 + (__a2 < __wc2);\
		__b3 += __wd3 + (__a3 < __wc3);\
		__b4 += __wd4 + (__a4 < __wc4);\
		__b5 += __wd5 + (__a5 < __wc5);\
		__b6 += __wd6 + (__a6 < __wc6);\
		__b7 += __wd7 + (__a7 < __wc7);\
		\
		/* Now add [w1,w2,w3] + [a,b,0]: */\
		__t0 += __a0;\
		__t1 += __a1;\
		__t2 += __a2;\
		__t3 += __a3;\
		__t4 += __a4;\
		__t5 += __a5;\
		__t6 += __a6;\
		__t7 += __a7;\
		\
		__cy0 = (__t0 < __a0);\
		__cy1 = (__t1 < __a1);\
		__cy2 = (__t2 < __a2);\
		__cy3 = (__t3 < __a3);\
		__cy4 = (__t4 < __a4);\
		__cy5 = (__t5 < __a5);\
		__cy6 = (__t6 < __a6);\
		__cy7 = (__t7 < __a7);\
		\
		__lo0.d0 += __cy0;\
		__lo1.d0 += __cy1;\
		__lo2.d0 += __cy2;\
		__lo3.d0 += __cy3;\
		__lo4.d0 += __cy4;\
		__lo5.d0 += __cy5;\
		__lo6.d0 += __cy6;\
		__lo7.d0 += __cy7;\
		\
		__lo0.d1 += (__lo0.d0 < __cy0);\
		__lo1.d1 += (__lo1.d0 < __cy1);\
		__lo2.d1 += (__lo2.d0 < __cy2);\
		__lo3.d1 += (__lo3.d0 < __cy3);\
		__lo4.d1 += (__lo4.d0 < __cy4);\
		__lo5.d1 += (__lo5.d0 < __cy5);\
		__lo6.d1 += (__lo6.d0 < __cy6);\
		__lo7.d1 += (__lo7.d0 < __cy7);\
		\
		__lo0.d0 += __b0;\
		__lo1.d0 += __b1;\
		__lo2.d0 += __b2;\
		__lo3.d0 += __b3;\
		__lo4.d0 += __b4;\
		__lo5.d0 += __b5;\
		__lo6.d0 += __b6;\
		__lo7.d0 += __b7;\
		\
		__lo0.d1 += (__lo0.d0 < __b0);\
		__lo1.d1 += (__lo1.d0 < __b1);\
		__lo2.d1 += (__lo2.d0 < __b2);\
		__lo3.d1 += (__lo3.d0 < __b3);\
		__lo4.d1 += (__lo4.d0 < __b4);\
		__lo5.d1 += (__lo5.d0 < __b5);\
		__lo6.d1 += (__lo6.d0 < __b6);\
		__lo7.d1 += (__lo7.d0 < __b7);\
    }

#endif

/* Upper half of 256-bit product of uint128 *x and *y, specialized to x,y < 2^96.
Result is returned in a uint64.

If we write x = a + b*2^32 and y = c + d*2^32 (i.e. a and c are 32 bits, b and d are 64 bits), the multiply
gives x*y = a*b + (a*d + b*c)*2^32 + b*d*2^64. a*d and b*c are both 96 bits, so there is a roughly 50%
chance that their sum will cause a carry bit to cascade into the 128th bit position. That means that we can't
simply neglect all the low-order stuff and simply do a MULH of x and y, but rather that we must do some
error correction in order to reconstruct the carryin from the lower 128 bits of the exact full-length product.
If we are willing to accept a roughly 1 in 2^32 chance of failure, we can simply calculate a*(d>>32) + (b>>32)*c
and see if, when added to the lower 64 half of the 128-bit product b*d, that sum overflows the 64th bit position.
That needs 2 32x32->64-bit MULL64s and a 64x64->64-bit MULL64. Thus, on Alpha, this needs a total of 4 MUL instructions.

	__x = a+b.B+x.B^2,   __y = c+d.B+y.B^2, where B = 2^32.
Then,
	__x*__y = a.c + (a.d+b.c).B + (a.y+b.d+x.c).B^2 + (b.y+x.d).B^3 + x.y.B^4

If (a.y+b.d+x.c) overflows 64 bits, the overflow bit must be added to the MULH128 result.
Similarly, (b.y+x.d)>>32 must be added to the MULH128 result.
*/
#if 0/*#ifdef MUL_LOHI64_SUBROUTINE */

    #define MULH128_96(__x, __y, __hi64)\
    {\
		uint64 __a,__b,__c,__d,__lo64;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MULH128_96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MULH128_96: (__y.d1 >> 32) == 0");\
		\
		__a  = (__x.d0) & (uint64)0x00000000ffffffff;\
		__b  = (__y.d0) & (uint64)0x00000000ffffffff;\
		__c  = ((__x.d0)>>32) + ((__x.d1)<<32);\
		__d  = ((__y.d0)>>32) + ((__y.d1)<<32);\
		__a *= __y.d1;\
		__b *= __x.d1;\
		__a += __b;\
		MUL_LOHI64(__c,__d,&__lo64,&__hi64);\
		__hi64 += (__a    < __b);/* Error correction step 1 is here. */\
		__lo64 += __a;\
		__hi64 += (__lo64 < __a);/* Error correction step 2 is here. */\
    }

#elif 0/*#else*/

    #define MULH128_96(__x, __y, __hi64)\
    {\
		uint64 __a,__b,__c,__d,__lo64;\
		\
		DBG_ASSERT((__x.d1 >> 32) == 0,"MULH128_96: (__x.d1 >> 32) == 0");\
		DBG_ASSERT((__y.d1 >> 32) == 0,"MULH128_96: (__y.d1 >> 32) == 0");\
		\
		__a  = (__x.d0) & (uint64)0x00000000ffffffff;\
		__b  = (__y.d0) & (uint64)0x00000000ffffffff;\
		__c  = ((__x.d0)>>32) + ((__x.d1)<<32);\
		__d  = ((__y.d0)>>32) + ((__y.d1)<<32);\
		__a *= __y.d1;\
		__b *= __x.d1;\
		__a += __b;\
		MUL_LOHI64(__c,__d, __lo64, __hi64);\
		__hi64 += (__a    < __b);/* Error correction step 1 is here. */\
		__lo64 += __a;\
		__hi64 += (__lo64 < __a);/* Error correction step 2 is here. */\
    }

#endif


/*****************************************************/
/*                    160 x 160                      */
/*****************************************************/

/* 320-bit square of uint160 *x. Result is returned in a pair of uint160s.

Let x = x0 + x1*2^64 + x2*2^128,  with x0, x1 < 2^64, x2 < 2^32. Then,

 x^2 = x0^2 + 2*x0*x1*2^64 + (x1^2 + 2*x0*x2)*2^128 + 2*x1*x2*2^192 + x2^2*2^256 .

In terms of the 5 output coefficients, here is what goes into each of those:

	w0 = (x0^2).lo
	w1 = (x0^2).hi + (2*x0*x1).lo
	w2 =             (2*x0*x1).hi + (x1^2 + 2*x0*x2).lo
	w3 =                            (x1^2 + 2*x0*x2).hi + (2*x1*x2).lo
	w4 =                                                  (2*x1*x2).hi + x2^2 .

On Alpha, this needs a total of 9 MUL instructions and 23 ALU ops, plus a few more
ALU ops to split the 5 64-bit outputs into a pair of uint160s.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define SQR_LOHI160(__x, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__w4,__a,__b,__t;\
		\
		DBG_ASSERT((__x.d2 >> 32) == 0,"SQR_LOHI160: (__x.d2 >> 32) == 0");\
		/* First calculate high partial products and put into w3 and w4: */\
		__t  = __x.d2;\
		__w4 = __t * __t;						/*   x2^2 */\
		__t  = __t + __t;						/* Store 2*x2 in a temporary */\
		MUL_LOHI64(__x.d1,__t   ,&__w3,&__a );	/* 2*x1*x2 */\
		__w4 += __a;\
		/* Now do low PPs: */\
		SQR_LOHI64(__x.d0,       &__w0,&__w1);	/*   x0^2 ; w0 done. */\
		MUL_LOHI64(__x.d0,__x.d1,&__a ,&__b );	/*   x0*x1 */\
		/* Double x0*x1 and add to x0^2: */\
		__w3 += __b >> 63;\
		__w2 = (__b << 1) + (__a >> 63);\
		__a <<= 1;\
		__w1 += __a;	__w2 += (__w1 < __a);	/*          w1 done. */\
		/* Now calculate (x1^2 + 2*x0*x2) and add into middle part (w2 and w3): */\
		SQR_LOHI64(__x.d1,       &__a ,&__b );	/*   x1^2  */\
		__w2 += __a;	__w3 += (__w2 < __a);\
		__w3 += __b;	__w4 += (__w3 < __b);\
		MUL_LOHI64(__x.d0,__t   ,&__a ,&__b );	/* 2*x0*x2 */\
		__w2 += __a;	__w3 += (__w2 < __a);\
		__w3 += __b;	__w4 += (__w3 < __b);	/* Chance of overflow here is just 1 in 2^31. */\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 =  __w0;					__lo.d1 = __w1;						__lo.d2 =  __w2 & 0x00000000ffffffff;\
		__hi.d0 = (__w2>>32) + (__w3<<32);	__hi.d1 = (__w3>>32) + (__w4<<32);	__hi.d2 = (__w4>>32);\
    }

    /* 4-operand-pipelined version: */
    #define SQR_LOHI160_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
    	SQR_LOHI160(__x0, __lo0, __hi0);\
    	SQR_LOHI160(__x1, __lo1, __hi1);\
    	SQR_LOHI160(__x2, __lo2, __hi2);\
    	SQR_LOHI160(__x3, __lo3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define SQR_LOHI160_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
    	SQR_LOHI160(__x0, __lo0, __hi0);\
    	SQR_LOHI160(__x1, __lo1, __hi1);\
    	SQR_LOHI160(__x2, __lo2, __hi2);\
    	SQR_LOHI160(__x3, __lo3, __hi3);\
    	SQR_LOHI160(__x4, __lo4, __hi4);\
    	SQR_LOHI160(__x5, __lo5, __hi5);\
    	SQR_LOHI160(__x6, __lo6, __hi6);\
    	SQR_LOHI160(__x7, __lo7, __hi7);\
	}

#else

    #define SQR_LOHI160(__x, __lo, __hi)\
    {\
		uint64 __w0,__w1,__w2,__w3,__w4,__a,__b,__t;\
		\
		DBG_ASSERT((__x.d2 >> 32) == 0,"SQR_LOHI160: (__x.d2 >> 32) == 0");\
		/* First calculate high partial products and put into w3 and w4: */\
		__t  = __x.d2;\
		__w4 = __t * __t;						/*   x2^2 */\
		__t  = __t + __t;						/* Store 2*x2 in a temporary */\
		MUL_LOHI64(__x.d1,__t   , __w3, __a );	/* 2*x1*x2 */\
		__w4 += __a;\
		/* Now do low PPs: */\
		SQR_LOHI64(__x.d0,        __w0, __w1);	/*   x0^2 ; w0 done. */\
		MUL_LOHI64(__x.d0,__x.d1, __a , __b );	/*   x0*x1 */\
		/* Double x0*x1 and add to x0^2: */\
		__w3 += __b >> 63;\
		__w2 = (__b << 1) + (__a >> 63);\
		__a <<= 1;\
		__w1 += __a;	__w2 += (__w1 < __a);	/*          w1 done. */\
		/* Now calculate (x1^2 + 2*x0*x2) and add into middle part (w2 and w3): */\
		SQR_LOHI64(__x.d1,        __a , __b );	/*   x1^2  */\
		__w2 += __a;	__w3 += (__w2 < __a);\
		__w3 += __b;	__w4 += (__w3 < __b);\
		MUL_LOHI64(__x.d0,__t   , __a , __b );	/* 2*x0*x2 */\
		__w2 += __a;	__w3 += (__w2 < __a);\
		__w3 += __b;	__w4 += (__w3 < __b);	/* Chance of overflow here is just 1 in 2^31. */\
		__lo.d0 =  __w0;					__lo.d1 = __w1;						__lo.d2 =  __w2 & 0x00000000ffffffff;\
		__hi.d0 = (__w2>>32) + (__w3<<32);	__hi.d1 = (__w3>>32) + (__w4<<32);	__hi.d2 = (__w4>>32);\
    }

    /* 4-operand-pipelined version: */
    #define SQR_LOHI160_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
	uint64	__a0,__a1,__a2,__a3,\
			__b0,__b1,__b2,__b3,\
			__t0,__t1,__t2,__t3,\
			__wa0,__wa1,__wa2,__wa3,\
			__wb0,__wb1,__wb2,__wb3,\
			__wc0,__wc1,__wc2,__wc3,\
			__wd0,__wd1,__wd2,__wd3,\
			__we0,__we1,__we2,__we3;\
		\
		DBG_ASSERT((__x0.d2 >> 32) == 0,"SQR_LOHI160_q4: (__x0.d2 >> 32) == 0");\
		DBG_ASSERT((__x1.d2 >> 32) == 0,"SQR_LOHI160_q4: (__x1.d2 >> 32) == 0");\
		DBG_ASSERT((__x2.d2 >> 32) == 0,"SQR_LOHI160_q4: (__x2.d2 >> 32) == 0");\
		DBG_ASSERT((__x3.d2 >> 32) == 0,"SQR_LOHI160_q4: (__x3.d2 >> 32) == 0");\
		\
		__t0  = __x0.d2;\
		__t1  = __x1.d2;\
		__t2  = __x2.d2;\
		__t3  = __x3.d2;\
		\
		__we0 = __t0 * __t0;\
		__we1 = __t1 * __t1;\
		__we2 = __t2 * __t2;\
		__we3 = __t3 * __t3;\
		\
		__t0  = __t0 + __t0;\
		__t1  = __t1 + __t1;\
		__t2  = __t2 + __t2;\
		__t3  = __t3 + __t3;\
		\
		MUL_LOHI64(__x0.d1,__t0   ,__wd0, __a0 );\
		MUL_LOHI64(__x1.d1,__t1   ,__wd1, __a1 );\
		MUL_LOHI64(__x2.d1,__t2   ,__wd2, __a2 );\
		MUL_LOHI64(__x3.d1,__t3   ,__wd3, __a3 );\
		\
		__we0 += __a0;\
		__we1 += __a1;\
		__we2 += __a2;\
		__we3 += __a3;\
		\
		/* Now do low PPs: */\
		SQR_LOHI64(__x0.d0,        __wa0, __wb0);\
		SQR_LOHI64(__x1.d0,        __wa1, __wb1);\
		SQR_LOHI64(__x2.d0,        __wa2, __wb2);\
		SQR_LOHI64(__x3.d0,        __wa3, __wb3);\
		\
		MUL_LOHI64(__x0.d0,__x0.d1,__a0 , __b0 );\
		MUL_LOHI64(__x1.d0,__x1.d1,__a1 , __b1 );\
		MUL_LOHI64(__x2.d0,__x2.d1,__a2 , __b2 );\
		MUL_LOHI64(__x3.d0,__x3.d1,__a3 , __b3 );\
		\
		/* Double x0*x1 and add to x0^2: */\
		__wd0 += __b0 >> 63;\
		__wd1 += __b1 >> 63;\
		__wd2 += __b2 >> 63;\
		__wd3 += __b3 >> 63;\
		\
		__wc0 = (__b0 << 1) + (__a0 >> 63);\
		__wc1 = (__b1 << 1) + (__a1 >> 63);\
		__wc2 = (__b2 << 1) + (__a2 >> 63);\
		__wc3 = (__b3 << 1) + (__a3 >> 63);\
		\
		__a0 <<= 1;\
		__a1 <<= 1;\
		__a2 <<= 1;\
		__a3 <<= 1;\
		\
		__wb0 += __a0;	__wc0 += (__wb0 < __a0);\
		__wb1 += __a1;	__wc1 += (__wb1 < __a1);\
		__wb2 += __a2;	__wc2 += (__wb2 < __a2);\
		__wb3 += __a3;	__wc3 += (__wb3 < __a3);\
		\
		/* Now calculate (x1^2 + 2*x0*x2) and add into middle part (wc and wd): */\
		SQR_LOHI64(__x0.d1,       __a0 , __b0 );\
		SQR_LOHI64(__x1.d1,       __a1 , __b1 );\
		SQR_LOHI64(__x2.d1,       __a2 , __b2 );\
		SQR_LOHI64(__x3.d1,       __a3 , __b3 );\
		\
		__wc0 += __a0;	__wd0 += (__wc0 < __a0);\
		__wc1 += __a1;	__wd1 += (__wc1 < __a1);\
		__wc2 += __a2;	__wd2 += (__wc2 < __a2);\
		__wc3 += __a3;	__wd3 += (__wc3 < __a3);\
		\
		__wd0 += __b0;	__we0 += (__wd0 < __b0);\
		__wd1 += __b1;	__we1 += (__wd1 < __b1);\
		__wd2 += __b2;	__we2 += (__wd2 < __b2);\
		__wd3 += __b3;	__we3 += (__wd3 < __b3);\
		\
		MUL_LOHI64(__x0.d0,__t0   ,__a0 , __b0 );\
		MUL_LOHI64(__x1.d0,__t1   ,__a1 , __b1 );\
		MUL_LOHI64(__x2.d0,__t2   ,__a2 , __b2 );\
		MUL_LOHI64(__x3.d0,__t3   ,__a3 , __b3 );\
		\
		__wc0 += __a0;	__wd0 += (__wc0 < __a0);\
		__wc1 += __a1;	__wd1 += (__wc1 < __a1);\
		__wc2 += __a2;	__wd2 += (__wc2 < __a2);\
		__wc3 += __a3;	__wd3 += (__wc3 < __a3);\
		\
		__wd0 += __b0;	__we0 += (__wd0 < __b0);\
		__wd1 += __b1;	__we1 += (__wd1 < __b1);\
		__wd2 += __b2;	__we2 += (__wd2 < __b2);\
		__wd3 += __b3;	__we3 += (__wd3 < __b3);\
		\
		/* Now split the result between __lo and __hi: */\
		__lo0.d0 =  __wa0;						__lo0.d1 = __wb0;						__lo0.d2 =  __wc0 & 0x00000000ffffffff;\
		__lo1.d0 =  __wa1;						__lo1.d1 = __wb1;						__lo1.d2 =  __wc1 & 0x00000000ffffffff;\
		__lo2.d0 =  __wa2;						__lo2.d1 = __wb2;						__lo2.d2 =  __wc2 & 0x00000000ffffffff;\
		__lo3.d0 =  __wa3;						__lo3.d1 = __wb3;						__lo3.d2 =  __wc3 & 0x00000000ffffffff;\
		\
		__hi0.d0 = (__wc0>>32) + (__wd0<<32);	__hi0.d1 = (__wd0>>32) + (__we0<<32);	__hi0.d2 = (__we0>>32);\
		__hi1.d0 = (__wc1>>32) + (__wd1<<32);	__hi1.d1 = (__wd1>>32) + (__we1<<32);	__hi1.d2 = (__we1>>32);\
		__hi2.d0 = (__wc2>>32) + (__wd2<<32);	__hi2.d1 = (__wd2>>32) + (__we2<<32);	__hi2.d2 = (__we2>>32);\
		__hi3.d0 = (__wc3>>32) + (__wd3<<32);	__hi3.d1 = (__wd3>>32) + (__we3<<32);	__hi3.d2 = (__we3>>32);\
    }

    /* 8-operand-pipelined version: */
    #define SQR_LOHI160_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
	uint64	__a0,__a1,__a2,__a3,__a4,__a5,__a6,__a7,\
			__b0,__b1,__b2,__b3,__b4,__b5,__b6,__b7,\
			__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,\
			__wa0,__wa1,__wa2,__wa3,__wa4,__wa5,__wa6,__wa7,\
			__wb0,__wb1,__wb2,__wb3,__wb4,__wb5,__wb6,__wb7,\
			__wc0,__wc1,__wc2,__wc3,__wc4,__wc5,__wc6,__wc7,\
			__wd0,__wd1,__wd2,__wd3,__wd4,__wd5,__wd6,__wd7,\
			__we0,__we1,__we2,__we3,__we4,__we5,__we6,__we7;\
		\
		DBG_ASSERT((__x0.d2 >> 32) == 0,"SQR_LOHI160_q8: (__x0.d2 >> 32) == 0");\
		DBG_ASSERT((__x1.d2 >> 32) == 0,"SQR_LOHI160_q8: (__x1.d2 >> 32) == 0");\
		DBG_ASSERT((__x2.d2 >> 32) == 0,"SQR_LOHI160_q8: (__x2.d2 >> 32) == 0");\
		DBG_ASSERT((__x3.d2 >> 32) == 0,"SQR_LOHI160_q8: (__x3.d2 >> 32) == 0");\
		DBG_ASSERT((__x4.d2 >> 32) == 0,"SQR_LOHI160_q8: (__x4.d2 >> 32) == 0");\
		DBG_ASSERT((__x5.d2 >> 32) == 0,"SQR_LOHI160_q8: (__x5.d2 >> 32) == 0");\
		DBG_ASSERT((__x6.d2 >> 32) == 0,"SQR_LOHI160_q8: (__x6.d2 >> 32) == 0");\
		DBG_ASSERT((__x7.d2 >> 32) == 0,"SQR_LOHI160_q8: (__x7.d2 >> 32) == 0");\
		\
		__t0  = __x0.d2;\
		__t1  = __x1.d2;\
		__t2  = __x2.d2;\
		__t3  = __x3.d2;\
		__t4  = __x4.d2;\
		__t5  = __x5.d2;\
		__t6  = __x6.d2;\
		__t7  = __x7.d2;\
		\
		__we0 = __t0 * __t0;\
		__we1 = __t1 * __t1;\
		__we2 = __t2 * __t2;\
		__we3 = __t3 * __t3;\
		__we4 = __t4 * __t4;\
		__we5 = __t5 * __t5;\
		__we6 = __t6 * __t6;\
		__we7 = __t7 * __t7;\
		\
		__t0  = __t0 + __t0;\
		__t1  = __t1 + __t1;\
		__t2  = __t2 + __t2;\
		__t3  = __t3 + __t3;\
		__t4  = __t4 + __t4;\
		__t5  = __t5 + __t5;\
		__t6  = __t6 + __t6;\
		__t7  = __t7 + __t7;\
		\
		MUL_LOHI64(__x0.d1,__t0   ,__wd0, __a0 );\
		MUL_LOHI64(__x1.d1,__t1   ,__wd1, __a1 );\
		MUL_LOHI64(__x2.d1,__t2   ,__wd2, __a2 );\
		MUL_LOHI64(__x3.d1,__t3   ,__wd3, __a3 );\
		MUL_LOHI64(__x4.d1,__t4   ,__wd4, __a4 );\
		MUL_LOHI64(__x5.d1,__t5   ,__wd5, __a5 );\
		MUL_LOHI64(__x6.d1,__t6   ,__wd6, __a6 );\
		MUL_LOHI64(__x7.d1,__t7   ,__wd7, __a7 );\
		\
		__we0 += __a0;\
		__we1 += __a1;\
		__we2 += __a2;\
		__we3 += __a3;\
		__we4 += __a4;\
		__we5 += __a5;\
		__we6 += __a6;\
		__we7 += __a7;\
		\
		/* Now do low PPs: */\
		SQR_LOHI64(__x0.d0,        __wa0, __wb0);\
		SQR_LOHI64(__x1.d0,        __wa1, __wb1);\
		SQR_LOHI64(__x2.d0,        __wa2, __wb2);\
		SQR_LOHI64(__x3.d0,        __wa3, __wb3);\
		SQR_LOHI64(__x4.d0,        __wa4, __wb4);\
		SQR_LOHI64(__x5.d0,        __wa5, __wb5);\
		SQR_LOHI64(__x6.d0,        __wa6, __wb6);\
		SQR_LOHI64(__x7.d0,        __wa7, __wb7);\
		\
		MUL_LOHI64(__x0.d0,__x0.d1,__a0 , __b0 );\
		MUL_LOHI64(__x1.d0,__x1.d1,__a1 , __b1 );\
		MUL_LOHI64(__x2.d0,__x2.d1,__a2 , __b2 );\
		MUL_LOHI64(__x3.d0,__x3.d1,__a3 , __b3 );\
		MUL_LOHI64(__x4.d0,__x4.d1,__a4 , __b4 );\
		MUL_LOHI64(__x5.d0,__x5.d1,__a5 , __b5 );\
		MUL_LOHI64(__x6.d0,__x6.d1,__a6 , __b6 );\
		MUL_LOHI64(__x7.d0,__x7.d1,__a7 , __b7 );\
		\
		/* Double x0*x1 and add to x0^2: */\
		__wd0 += __b0 >> 63;\
		__wd1 += __b1 >> 63;\
		__wd2 += __b2 >> 63;\
		__wd3 += __b3 >> 63;\
		__wd4 += __b4 >> 63;\
		__wd5 += __b5 >> 63;\
		__wd6 += __b6 >> 63;\
		__wd7 += __b7 >> 63;\
		\
		__wc0 = (__b0 << 1) + (__a0 >> 63);\
		__wc1 = (__b1 << 1) + (__a1 >> 63);\
		__wc2 = (__b2 << 1) + (__a2 >> 63);\
		__wc3 = (__b3 << 1) + (__a3 >> 63);\
		__wc4 = (__b4 << 1) + (__a4 >> 63);\
		__wc5 = (__b5 << 1) + (__a5 >> 63);\
		__wc6 = (__b6 << 1) + (__a6 >> 63);\
		__wc7 = (__b7 << 1) + (__a7 >> 63);\
		\
		__a0 <<= 1;\
		__a1 <<= 1;\
		__a2 <<= 1;\
		__a3 <<= 1;\
		__a4 <<= 1;\
		__a5 <<= 1;\
		__a6 <<= 1;\
		__a7 <<= 1;\
		\
		__wb0 += __a0;	__wc0 += (__wb0 < __a0);\
		__wb1 += __a1;	__wc1 += (__wb1 < __a1);\
		__wb2 += __a2;	__wc2 += (__wb2 < __a2);\
		__wb3 += __a3;	__wc3 += (__wb3 < __a3);\
		__wb4 += __a4;	__wc4 += (__wb4 < __a4);\
		__wb5 += __a5;	__wc5 += (__wb5 < __a5);\
		__wb6 += __a6;	__wc6 += (__wb6 < __a6);\
		__wb7 += __a7;	__wc7 += (__wb7 < __a7);\
		\
		/* Now calculate (x1^2 + 2*x0*x2) and add into middle part (wc and wd): */\
		SQR_LOHI64(__x0.d1,       __a0 , __b0 );\
		SQR_LOHI64(__x1.d1,       __a1 , __b1 );\
		SQR_LOHI64(__x2.d1,       __a2 , __b2 );\
		SQR_LOHI64(__x3.d1,       __a3 , __b3 );\
		SQR_LOHI64(__x4.d1,       __a4 , __b4 );\
		SQR_LOHI64(__x5.d1,       __a5 , __b5 );\
		SQR_LOHI64(__x6.d1,       __a6 , __b6 );\
		SQR_LOHI64(__x7.d1,       __a7 , __b7 );\
		\
		__wc0 += __a0;	__wd0 += (__wc0 < __a0);\
		__wc1 += __a1;	__wd1 += (__wc1 < __a1);\
		__wc2 += __a2;	__wd2 += (__wc2 < __a2);\
		__wc3 += __a3;	__wd3 += (__wc3 < __a3);\
		__wc4 += __a4;	__wd4 += (__wc4 < __a4);\
		__wc5 += __a5;	__wd5 += (__wc5 < __a5);\
		__wc6 += __a6;	__wd6 += (__wc6 < __a6);\
		__wc7 += __a7;	__wd7 += (__wc7 < __a7);\
		\
		__wd0 += __b0;	__we0 += (__wd0 < __b0);\
		__wd1 += __b1;	__we1 += (__wd1 < __b1);\
		__wd2 += __b2;	__we2 += (__wd2 < __b2);\
		__wd3 += __b3;	__we3 += (__wd3 < __b3);\
		__wd4 += __b4;	__we4 += (__wd4 < __b4);\
		__wd5 += __b5;	__we5 += (__wd5 < __b5);\
		__wd6 += __b6;	__we6 += (__wd6 < __b6);\
		__wd7 += __b7;	__we7 += (__wd7 < __b7);\
		\
		MUL_LOHI64(__x0.d0,__t0   ,__a0 , __b0 );\
		MUL_LOHI64(__x1.d0,__t1   ,__a1 , __b1 );\
		MUL_LOHI64(__x2.d0,__t2   ,__a2 , __b2 );\
		MUL_LOHI64(__x3.d0,__t3   ,__a3 , __b3 );\
		MUL_LOHI64(__x4.d0,__t4   ,__a4 , __b4 );\
		MUL_LOHI64(__x5.d0,__t5   ,__a5 , __b5 );\
		MUL_LOHI64(__x6.d0,__t6   ,__a6 , __b6 );\
		MUL_LOHI64(__x7.d0,__t7   ,__a7 , __b7 );\
		\
		__wc0 += __a0;	__wd0 += (__wc0 < __a0);\
		__wc1 += __a1;	__wd1 += (__wc1 < __a1);\
		__wc2 += __a2;	__wd2 += (__wc2 < __a2);\
		__wc3 += __a3;	__wd3 += (__wc3 < __a3);\
		__wc4 += __a4;	__wd4 += (__wc4 < __a4);\
		__wc5 += __a5;	__wd5 += (__wc5 < __a5);\
		__wc6 += __a6;	__wd6 += (__wc6 < __a6);\
		__wc7 += __a7;	__wd7 += (__wc7 < __a7);\
		\
		__wd0 += __b0;	__we0 += (__wd0 < __b0);\
		__wd1 += __b1;	__we1 += (__wd1 < __b1);\
		__wd2 += __b2;	__we2 += (__wd2 < __b2);\
		__wd3 += __b3;	__we3 += (__wd3 < __b3);\
		__wd4 += __b4;	__we4 += (__wd4 < __b4);\
		__wd5 += __b5;	__we5 += (__wd5 < __b5);\
		__wd6 += __b6;	__we6 += (__wd6 < __b6);\
		__wd7 += __b7;	__we7 += (__wd7 < __b7);\
		\
		/* Now split the result between __lo and __hi: */\
		__lo0.d0 =  __wa0;						__lo0.d1 = __wb0;						__lo0.d2 =  __wc0 & 0x00000000ffffffff;\
		__lo1.d0 =  __wa1;						__lo1.d1 = __wb1;						__lo1.d2 =  __wc1 & 0x00000000ffffffff;\
		__lo2.d0 =  __wa2;						__lo2.d1 = __wb2;						__lo2.d2 =  __wc2 & 0x00000000ffffffff;\
		__lo3.d0 =  __wa3;						__lo3.d1 = __wb3;						__lo3.d2 =  __wc3 & 0x00000000ffffffff;\
		__lo4.d0 =  __wa4;						__lo4.d1 = __wb4;						__lo0.d2 =  __wc0 & 0x00000000ffffffff;\
		__lo5.d0 =  __wa5;						__lo5.d1 = __wb5;						__lo1.d2 =  __wc1 & 0x00000000ffffffff;\
		__lo6.d0 =  __wa6;						__lo6.d1 = __wb6;						__lo2.d2 =  __wc2 & 0x00000000ffffffff;\
		__lo7.d0 =  __wa7;						__lo7.d1 = __wb7;						__lo3.d2 =  __wc3 & 0x00000000ffffffff;\
		\
		__hi0.d0 = (__wc0>>32) + (__wd0<<32);	__hi0.d1 = (__wd0>>32) + (__we0<<32);	__hi0.d2 = (__we0>>32);\
		__hi1.d0 = (__wc1>>32) + (__wd1<<32);	__hi1.d1 = (__wd1>>32) + (__we1<<32);	__hi1.d2 = (__we1>>32);\
		__hi2.d0 = (__wc2>>32) + (__wd2<<32);	__hi2.d1 = (__wd2>>32) + (__we2<<32);	__hi2.d2 = (__we2>>32);\
		__hi3.d0 = (__wc3>>32) + (__wd3<<32);	__hi3.d1 = (__wd3>>32) + (__we3<<32);	__hi3.d2 = (__we3>>32);\
		__hi4.d0 = (__wc4>>32) + (__wd4<<32);	__hi4.d1 = (__wd4>>32) + (__we4<<32);	__hi4.d2 = (__we4>>32);\
		__hi5.d0 = (__wc5>>32) + (__wd5<<32);	__hi5.d1 = (__wd5>>32) + (__we5<<32);	__hi5.d2 = (__we5>>32);\
		__hi6.d0 = (__wc6>>32) + (__wd6<<32);	__hi6.d1 = (__wd6>>32) + (__we6<<32);	__hi6.d2 = (__we6>>32);\
		__hi7.d0 = (__wc7>>32) + (__wd7<<32);	__hi7.d1 = (__wd7>>32) + (__we7<<32);	__hi7.d2 = (__we7>>32);\
    }

#endif

/* 160-bit product modulo 2^160 of uint160 *x and *y. Result is returned in a uint160.
   Designed so output may be returned in input argument x,
   if x and lo happen to point to the same addresses.

Let x = x0 + x1*2^64 + x2*2^128,  with x0, x1 < 2^64, x2 < 2^32,
and y = y0 + y1*2^64 + y2*2^128,  with y0, y1 < 2^64, y2 < 2^32. Then,

 x*y = x0*y0 + (x0*y1 + x1*y0)*2^64 + (x1*y1 + x2*y0 + y2*x0)*2^128 + ...

In terms of the 3 output coefficients, here is what goes into each of those:

	w0 = (x0*y0).lo
	w1 = (x0*y0).hi + (x0*y1 + x1*y0).lo
	w2 =              (x0*y1 + x1*y0).hi + (x1*y1 + x2*y0 + y2*x0).lo ,

where we need only the bottom 32 bits of the (x0*y1 + x1*y0).hi and
(x1*y1 + x2*y0 + y2*x0).lo terms, i.e. can use a trio of MULLs for the latter.

On Alpha, this needs a total of 9 MUL (3 MULL32, 3 MULL64, 3 MULH64), 10 ALU.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULL160(__x, __y, __lo)\
    {\
		uint64 __w0,__w1,__w2,__a,__b,__c,__d;\
		\
		MUL_LOHI64(__x.d0,__y.d0,&__w0,&__w1);	/*   x0*y0 */\
		/* Only need the bottom 32 bits of (x1*y1 + x2*y0 + y2*x0).lo: */\
		__w2 = __MULL32(__x.d1,__y.d1)+__MULL32(__x.d0,__y.d2)+__MULL32(__y.d0,__x.d2);\
		\
		MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0,&__c ,&__d );	/*   x1*y0 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__w2 += (__w1 < __a);\
		__w2 += __b;\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__w2 += (__w1 < __c);\
		__w2 += __d;\
		\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2 & 0x00000000ffffffff;\
    }

    /* 4-operand-pipelined version: */

#else

    #define MULL160(__x, __y, __lo)\
    {\
		uint64 __w0,__w1,__w2,__a,__b,__c,__d;\
		\
		MUL_LOHI64(__x.d0,__y.d0, __w0, __w1);	/*   x0*y0 */\
		/* Only need the bottom 32 bits of (x1*y1 + x2*y0 + y2*x0).lo: */\
		__w2 = __MULL32(__x.d1,__y.d1)+__MULL32(__x.d0,__y.d2)+__MULL32(__y.d0,__x.d2);\
		\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0, __c , __d );	/*   x1*y0 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__w2 += (__w1 < __a);\
		__w2 += __b;\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__w2 += (__w1 < __c);\
		__w2 += __d;\
		\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2 & 0x00000000ffffffff;\
    }

#endif

    /* 4-operand-pipelined version: */
    #define MULL160_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
		MULL160(__x0, __y0, __lo0);\
		MULL160(__x1, __y1, __lo1);\
		MULL160(__x2, __y2, __lo2);\
		MULL160(__x3, __y3, __lo3);\
	}

    /* 8-operand-pipelined version: */
    #define MULL160_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
		MULL160(__x0, __y0, __lo0);\
		MULL160(__x1, __y1, __lo1);\
		MULL160(__x2, __y2, __lo2);\
		MULL160(__x3, __y3, __lo3);\
		MULL160(__x4, __y4, __lo4);\
		MULL160(__x5, __y5, __lo5);\
		MULL160(__x6, __y6, __lo6);\
		MULL160(__x7, __y7, __lo7);\
	}

/* Upper 160 bits of product of uint160 x and y. Result is returned in a uint160.

Let x = x0 + x1*2^64 + x2*2^128,  with x0, x1 < 2^64, x2 < 2^32,
and y = y0 + y1*2^64 + y2*2^128,  with y0, y1 < 2^64, y2 < 2^32. Then,

 x*y = x0*y0 + (x0*y1 + x1*y0)*2^64 + (x1*y1 + x2*y0 + y2*x0)*2^128 + (x2*y1 + y2*x1)*2^192 + x2*y2*2^256.

In terms of the 5 full-length-product output coefficients, here is what goes into each of those:

	w0 = (x0*y0).lo
	w1 = (x0*y0).hi + (x0*y1 + x1*y0).lo
	w2 =              (x0*y1 + x1*y0).hi + (x1*y1 + x2*y0 + y2*x0).lo
	w3 =                                   (x1*y1 + x2*y0 + y2*x0).hi + (x2*y1 + y2*x1).lo
	w4 =                                                                (x2*y1 + y2*x1).hi + (x2*y2).lo .

On Alpha, this needs a total of 11 MUL instructions and 22 ALU ops.

On 32-bit hardware, take advantage of the fact that x2 and y2 are only 32 bits wide.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULH160(__x, __y, __hi)\
    {\
		uint64 __w1,__w2,__w3,__w4,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		DBG_ASSERT((__x.d2 >> 32) == 0,"MULH160: (__x.d2 >> 32) == 0");\
		DBG_ASSERT((__y.d2 >> 32) == 0,"MULH160: (__y.d2 >> 32) == 0");\
		\
		__w4 = __x.d2*__y.d2;				/*   x2*y2 */\
		MULH64(__x.d0,__y.d0,       __w1);	/*   x0*y0.hi */\
		MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);	/*   x1*y1 */\
		\
		MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0,&__c ,&__d );	/*   x1*y0 */\
		\
		MUL64x32(__x.d0,__y.d2,&__e ,&__f );	/*   x0*y2 */\
		MUL64x32(__x.d2,__y.d0,&__g ,&__h );	/*   x2*y0 */\
		\
		MUL64x32(__x.d1,__y.d2,&__i ,&__j );	/*   x1*y2 */\
		MUL64x32(__x.d2,__y.d1,&__k ,&__l );	/*   x2*y1 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__w2 += (__w1 < __a);\
		__w2 += __b;	__w3 += (__w2 < __b);\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__w2 += (__w1 < __c);\
		__w2 += __d;	__w3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-4: */\
		__w2 += __e;	__w3 += (__w2 < __e);\
		__w3 += __f;	__w4 += (__w3 < __f);\
		\
		/* Add x2*y0 to w2-4: */\
		__w2 += __g;	__w3 += (__w2 < __g);\
		__w3 += __h;	__w4 += (__w3 < __h);\
		\
		/* Add x1*y2 to w3-4: */\
		__w3 += __i;	__w4 += (__w3 < __i);\
		__w4 += __j;\
		\
		/* Add x2*y1 to w3-4: */\
		__w3 += __k;	__w4 += (__w3 < __k);\
		__w4 += __l;\
		\
		__hi.d0 = (__w2>>32) + (__w3<<32);	__hi.d1 = (__w3>>32) + (__w4<<32);	__hi.d2 = (__w4>>32);\
    }

    /* 4-operand-pipelined version: */
    #define MULH160_q4(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3)\
    {\
    	MULH160(__x0, __y0, __hi0);\
    	MULH160(__x1, __y1, __hi1);\
    	MULH160(__x2, __y2, __hi2);\
    	MULH160(__x3, __y3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define MULH160_q8(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3\
    , __x4, __y4, __hi4\
    , __x5, __y5, __hi5\
    , __x6, __y6, __hi6\
    , __x7, __y7, __hi7)\
    {\
    	MULH160(__x0, __y0, __hi0);\
    	MULH160(__x1, __y1, __hi1);\
    	MULH160(__x2, __y2, __hi2);\
    	MULH160(__x3, __y3, __hi3);\
    	MULH160(__x4, __y4, __hi4);\
    	MULH160(__x5, __y5, __hi5);\
    	MULH160(__x6, __y6, __hi6);\
    	MULH160(__x7, __y7, __hi7);\
	}

#else

    #define MULH160(__x, __y, __hi)\
    {\
		uint64 __w1,__w2,__w3,__w4,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		DBG_ASSERT((__x.d2 >> 32) == 0,"MULH160: (__x.d2 >> 32) == 0");\
		DBG_ASSERT((__y.d2 >> 32) == 0,"MULH160: (__y.d2 >> 32) == 0");\
		\
		__w4 = __x.d2*__y.d2;				/*   x2*y2 */\
		MULH64(__x.d0,__y.d0,       __w1);	/*   x0*y0.hi */\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);	/*   x1*y1 */\
		\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0, __c , __d );	/*   x1*y0 */\
		\
		MUL64x32(__x.d0,__y.d2, __e , __f );	/*   x0*y2 */\
		MUL64x32(__x.d2,__y.d0, __g , __h );	/*   x2*y0 */\
		\
		MUL64x32(__x.d1,__y.d2, __i , __j );	/*   x1*y2 */\
		MUL64x32(__x.d2,__y.d1, __k , __l );	/*   x2*y1 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__w2 += (__w1 < __a);\
		__w2 += __b;	__w3 += (__w2 < __b);\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__w2 += (__w1 < __c);\
		__w2 += __d;	__w3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-4: */\
		__w2 += __e;	__w3 += (__w2 < __e);\
		__w3 += __f;	__w4 += (__w3 < __f);\
		\
		/* Add x2*y0 to w2-4: */\
		__w2 += __g;	__w3 += (__w2 < __g);\
		__w3 += __h;	__w4 += (__w3 < __h);\
		\
		/* Add x1*y2 to w3-4: */\
		__w3 += __i;	__w4 += (__w3 < __i);\
		__w4 += __j;\
		\
		/* Add x2*y1 to w3-4: */\
		__w3 += __k;	__w4 += (__w3 < __k);\
		__w4 += __l;\
		\
		__hi.d0 = (__w2>>32) + (__w3<<32);	__hi.d1 = (__w3>>32) + (__w4<<32);	__hi.d2 = (__w4>>32);\
    }

    /* 4-operand-pipelined version: */
    #define MULH160_q4(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3)\
    {\
    	MULH160(__x0, __y0, __hi0);\
    	MULH160(__x1, __y1, __hi1);\
    	MULH160(__x2, __y2, __hi2);\
    	MULH160(__x3, __y3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define MULH160_q8(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3\
    , __x4, __y4, __hi4\
    , __x5, __y5, __hi5\
    , __x6, __y6, __hi6\
    , __x7, __y7, __hi7)\
    {\
    	MULH160(__x0, __y0, __hi0);\
    	MULH160(__x1, __y1, __hi1);\
    	MULH160(__x2, __y2, __hi2);\
    	MULH160(__x3, __y3, __hi3);\
    	MULH160(__x4, __y4, __hi4);\
    	MULH160(__x5, __y5, __hi5);\
    	MULH160(__x6, __y6, __hi6);\
    	MULH160(__x7, __y7, __hi7);\
	}

#endif


/*****************************************************/
/*                    192 x 192                      */
/*****************************************************/

/* Full 384 bits of the product of uint192 x and y.

Let x = x0 + x1*2^64 + x2*2^128,  with x0, x1, x2 < 2^64,
and y = y0 + y1*2^64 + y2*2^128,  with y0, y1, y2 < 2^64.            Then,

 x*y = x0*y0 + (x0*y1 + x1*y0)*2^64 + (x0*y2 + x1*y1 + x2*y0)*2^128 + (x2*y1 + x1*y2)*2^192 + (x2*y2)*2^256 .

In terms of the 6 output coefficients, here is what goes into each of those:

	w0 = (x0*y0).lo
	w1 = (x0*y0).hi + (x0*y1 + x1*y0).lo
	w2 =              (x0*y1 + x1*y0).hi + (x0*y2 + x1*y1 + x2*y0).lo
	w3 =                                   (x0*y2 + x1*y1 + x2*y0).hi + (x2*y1 + x1*y2).lo
	w4 =                                                                (x2*y1 + x1*y2).hi + (x2*y2).lo
	w5 =                                                                                     (x2*y2).hi .

On Alpha, this needs a total of 18 MUL instructions and 36 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MUL_LOHI192(__x, __y, __lo, __hi)\
    {\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5\
						,__cy2,__cy3,__cy4,__cy5\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		MUL_LOHI64(__x.d0,__y.d0,&__w0,&__w1);	/*   x0*y0 */\
		MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);	/*   x1*y1 */\
		MUL_LOHI64(__x.d2,__y.d2,&__w4,&__w5);	/*   x2*y2 */\
		\
		MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0,&__c ,&__d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2,&__e ,&__f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0,&__g ,&__h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x.d1,__y.d2,&__i ,&__j );	/*   x1*y2 */\
		MUL_LOHI64(__x.d2,__y.d1,&__k ,&__l );	/*   x2*y1 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-4: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-4: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x1*y2 to w3-5: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x2*y1 to w3-5: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

#else

    #define MUL_LOHI192(__x, __y, __lo, __hi)\
    {\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5\
						,__cy2,__cy3,__cy4,__cy5\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		MUL_LOHI64(__x.d0,__y.d0, __w0, __w1);	/*   x0*y0 */\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);	/*   x1*y1 */\
		MUL_LOHI64(__x.d2,__y.d2, __w4, __w5);	/*   x2*y2 */\
		\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0, __c , __d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2, __e , __f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0, __g , __h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x.d1,__y.d2, __i , __j );	/*   x1*y2 */\
		MUL_LOHI64(__x.d2,__y.d1, __k , __l );	/*   x2*y1 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-4: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-4: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x1*y2 to w3-5: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x2*y1 to w3-5: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

#endif


/* 384-bit square of uint192 x. Result is returned in a pair of uint192.

On Alpha, this needs a total of 12 MUL instructions and 42 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define SQR_LOHI192(__x, __lo, __hi)\
    {\
	/* In the output assembly below, we use MUL outputs in the following order: */	\
	/*     a,b,w3,w2,w1    c,d,w4    e,f,w5    */									\
	/* IDEA: rearrange and interleave the MULs to try to take advantage of that? */	\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5\
						,__cy2,__cy3,__cy4,__cy5\
		,__a,__b,__e,__f,__i,__j;\
		\
		SQR_LOHI64(__x.d0,       &__w0,&__w1);	/*   x0^2  */\
		MUL_LOHI64(__x.d0,__x.d1,&__a ,&__b );	/*   x0*x1 */\
		SQR_LOHI64(__x.d1,       &__w2,&__w3);	/*   x1^2  */\
		MUL_LOHI64(__x.d0,__x.d2,&__e ,&__f );	/*   x0*x2 */\
		SQR_LOHI64(__x.d2,       &__w4,&__w5);	/*   x2^2  */\
		MUL_LOHI64(__x.d1,__x.d2,&__i ,&__j );	/*   x1*x2 */\
		\
		/* Now add cross terms: */\
		/* Add x0*x1 twice to w1-3: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		__w1 += __a;	__cy2 += (__w1 < __a);\
		__w2 += __b;	__cy3 += (__w2 < __b);\
		\
		/* Add x0*x2 twice to w2-4: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4 += (__w3 < __f);\
		\
		/* Add x1*y2 twice to w3-5: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5 += (__w4 < __j);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

    #define SQR_LOHI192_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
    	SQR_LOHI192(__x0, __lo0, __hi0);\
    	SQR_LOHI192(__x1, __lo1, __hi1);\
    	SQR_LOHI192(__x2, __lo2, __hi2);\
    	SQR_LOHI192(__x3, __lo3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define SQR_LOHI192_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
    	SQR_LOHI192(__x0, __lo0, __hi0);\
    	SQR_LOHI192(__x1, __lo1, __hi1);\
    	SQR_LOHI192(__x2, __lo2, __hi2);\
    	SQR_LOHI192(__x3, __lo3, __hi3);\
    	SQR_LOHI192(__x4, __lo4, __hi4);\
    	SQR_LOHI192(__x5, __lo5, __hi5);\
    	SQR_LOHI192(__x6, __lo6, __hi6);\
    	SQR_LOHI192(__x7, __lo7, __hi7);\
	}

#else

  #if 0	// In generic-C mode, prefer faster arrangement of this below

    #define SQR_LOHI192(__x, __lo, __hi)\
    {\
	/* In the output assembly below, we use MUL outputs in the following order: */	\
	/*     w0,w1,w2,a,b    c,d,w4    e,f,w5    */									\
	/* IDEA: rearrange and interleave the MULs to try to take advantage of that? */	\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5\
						,__cy2,__cy3,__cy4,__cy5\
		,__a,__b,__e,__f,__i,__j;\
		\
		SQR_LOHI64(__x.d0,        __w0, __w1);	/*   x0^2  */\
		MUL_LOHI64(__x.d0,__x.d1, __a , __b );	/*   x0*x1 */\
		SQR_LOHI64(__x.d1,        __w2, __w3);	/*   x1^2  */\
		MUL_LOHI64(__x.d0,__x.d2, __e , __f );	/*   x0*x2 */\
		SQR_LOHI64(__x.d2,        __w4, __w5);	/*   x2^2  */\
		MUL_LOHI64(__x.d1,__x.d2, __i , __j );	/*   x1*x2 */\
		\
		/* Now add cross terms: */\
		/* Add x0*x1 twice to w1-3: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		__w1 += __a;	__cy2 += (__w1 < __a);\
		__w2 += __b;	__cy3 += (__w2 < __b);\
		\
		/* Add x0*x2 twice to w2-4: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4 += (__w3 < __f);\
		\
		/* Add x1*y2 twice to w3-5: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5 += (__w4 < __j);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2; __hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

  #elif(!defined(YES_ASM))	// Standard-C version #2:

    #define SQR_LOHI192(__x, __lo, __hi)\
    {\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5, __a,__b,__e,__f,__i,__j;\
		/* Arrange MUL results in layers, like so: */\
		/* w5:w4:w3:w2:w1:w0 = x2*x2 | x1*x1 | x0*x0 */\
		/*    _j:_i:_b:_a    =     x2*x1    x1*x0    */\
		/*       _f:_e       =         x0*x2         */\
		/* Evaluate the cross-terms first, then right-shift those 1 place while we compute the w0-w5 terms: */\
		MUL_LOHI64(__x.d0,__x.d1, __a , __b );	/*   x0*x1 */\
		MUL_LOHI64(__x.d0,__x.d2, __e , __f );	/*   x0*x2 */\
		MUL_LOHI64(__x.d1,__x.d2, __i , __j );	/*   x1*x2 */\
		/* Add last pair (bottom layer) cross terms: */\
		__b += __e;	__i += __f + (__b < __e);	/* f is a UMULH output, so f + CY cannot overflow */\
		__hi.d2 = ((int64)__j < 0);\
		/* Interleave left-shift of the 4-word result with squaring MULs: */\
		/*** I wish I could say that interleaving the MULs and ALU stuff gives a speedup, but on x86_64 I get same speed if I do */\
		/* all the MULs right at the start ... likely the microcode engine ends up giving similarly runtime-restructured code. ***/\
		SQR_LOHI64(__x.d0,        __w0, __w1);	__lo.d0 = __w0;	/*   x0^2  */\
		__j = (__j << 1) + ((int64)__i < 0);\
		__i = (__i << 1) + ((int64)__b < 0);\
		SQR_LOHI64(__x.d1,        __w2, __w3);	/*   x1^2  */\
		__b = (__b << 1) + ((int64)__a < 0);\
		__a = (__a << 1);\
		SQR_LOHI64(__x.d2,        __w4, __w5);	/*   x2^2  */\
		/* Final add/carry step: */\
		__w1 += __a;\
		__w2 += __b + (__w1 < __a);	__lo.d1 = __w1;\
		__w3 += __i + (__w2 < __b);	__lo.d2 = __w2;\
		__w4 += __j + (__w3 < __i);	__hi.d0 = __w3;\
	__hi.d2 += __w5 + (__w4 < __j);	__hi.d1 = __w4;\
		/* Result now split between __lo and __hi. Cost: 15 ADD, 4 SHL, 9 CMP, versus original-arrangement */\
		/* cost of 25 ADD, 13 CMP. On x86_64 get ~6% overall savings (i.e. this routine probably ~20% faster) */\
    }

  #elif 1	// x86_64 ASM version (non-SSE):

	/* x86_64 ASM implementation of the 192-bit SQR_LOHI macros: */
    #define SQR_LOHI192(__x, __lo, __hi)\
    {\
	__asm__ volatile (\
		/* Evaluate the cross-terms first, then right-shift those 1 place while we compute the w0-w5 terms: */\
	"movq	%[__x0],%%rax	\n\t	mulq	%[__x1]	\n\t	movq	%%rax,%%r11	\n\t	movq	%%rdx,%%r12	\n\t"/* _a:_b = x0*x1 */\
	"movq	%[__x0],%%rax	\n\t	mulq	%[__x2]	\n\t	movq	%%rax,%%rdi	\n\t	movq	%%rdx,%%r13	\n\t"/* _e:_f = x0*x2 */\
	"movq	%[__x1],%%rax	\n\t	mulq	%[__x2]	\n\t		movq %%rdx,%%r14  \n\t	movq	%%rdx,%%r15	\n\t"/* _i:_j = x1*x2; leave _i in rax, make 2 copies of rdx */\
		/* Add last pair (bottom layer) cross terms: */\
		"addq	%%rdi,%%r12		\n\t"/* _b += _e */\
		"adcq	%%rax,%%r13		\n\t"/*	_i += _f + (_b < _e); _f is a UMULH output, so f + CY cannot overflow */\
		"shrq	$63,%%r15		\n\t"/* _w5 = (_j >> 63) */\
		/* Interleave left-shift of the 4-word result with squaring MULs: */\
	"movq	%[__x0],%%rax	\n\t	mulq	%%rax	\n\t	movq	%%rax,%[__lo0] \n\t	movq	%%rdx,%%r10	\n\t"/* w0:w1 = x0^2 */\
		"shldq	$1,%%r13,%%r14	\n\t"/* _j = (_j << 1) + ((int64)_i < 0) */\
		"shldq	$1,%%r12,%%r13	\n\t"/* _i = (_i << 1) + ((int64)_b < 0) */\
	"movq	%[__x1],%%rax	\n\t	mulq	%%rax	\n\t	movq	%%rax,%%rbx	\n\t	movq	%%rdx,%%rdi	\n\t"/* w2:w3 = x1^2 */\
		"shldq	$1,%%r11,%%r12	\n\t"/* _b = (_b << 1) + ((int64)_a < 0) */\
		"shlq	$1      ,%%r11	\n\t"/* _a = (_a << 1) */\
	"movq	%[__x2],%%rax	\n\t	mulq	%%rax	\n\t"/* w4:rdx = x2^2; rdx will get added to w5 below */\
		"addq	%%r10,%%r11		\n\t"/* _w1 += _a */\
		"adcq	%%rbx,%%r12		\n\t	movq	%%r11,%[__lo1] \n\t"/*	_w2 += _b +  (_w1 < _a);	_lo.d1 = _w1 */\
		"adcq	%%rdi,%%r13		\n\t	movq	%%r12,%[__lo2] \n\t"/*	_w3 += _i +  (_w2 < _b);	_lo.d2 = _w2 */\
		"adcq	%%rax,%%r14		\n\t	movq	%%r13,%[__hi0] \n\t"/*	_w4 += _j +  (_w3 < _i);	_hi.d0 = _w3 */\
		"adcq	%%rdx,%%r15		\n\t	movq	%%r14,%[__hi1] \n\t"/*	_w5 += rdx + (_w4 < _j);	_hi.d1 = _w4 */\
	"movq	%%r15,%[__hi2]	\n\t"/* _hi.d2 = _w5 */\
		: [__lo0] "=m" (__lo.d0)	/* This macro would only compile clean with outputs like so ...  "g" form below gives "hi.d0-2 uninitialized" warnings */\
		 ,[__lo1] "=m" (__lo.d1)	\
		 ,[__lo2] "=m" (__lo.d2)	\
		 ,[__hi0] "=m" (__hi.d0)	\
		 ,[__hi1] "=m" (__hi.d1)	\
		 ,[__hi2] "=m" (__hi.d2)	\
		: [__x0] "g" (__x.d0)	/* Inputs from memory/register here */\
		 ,[__x1] "g" (__x.d1)	\
		 ,[__x2] "g" (__x.d2)	\
		: "cc","memory","rax","rbx","rdx","rdi","r10","r11","r12","r13","r14","r15"	/* Clobbered registers [all but rcx] */\
	);\
	}

	/*
	To-Do: Try combining scalar MUL and SSE packed-quadword ADD/CMP:
	P[ADD|SUB]Q
	PCMP[EQ|GT]Q - Note only == and > supported, others must be synthesized from these. True = 0xffff...ffff --> subtract-to-carry
	PUNPCK[H|L]QDQ
	UNPCK[H|L]PD — Redundant Unpack, formally defined for doubles
	*/

  #endif	// YES_ASM

  #if 0
    /* 4-operand-pipelined version: */
    #define SQR_LOHI192_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
    	SQR_LOHI192(__x0, __lo0, __hi0);\
    	SQR_LOHI192(__x1, __lo1, __hi1);\
    	SQR_LOHI192(__x2, __lo2, __hi2);\
    	SQR_LOHI192(__x3, __lo3, __hi3);\
	}
  #else
    #define SQR_LOHI192_q4(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3)\
    {\
	/* In the output assembly below, we use MUL outputs in the following order: */	\
	/*     w0,w1,w2,a,b    c,d,w4    e,f,w5    */									\
	/* IDEA: rearrange and interleave the MULs to try to take advantage of that? */	\
		uint64 __wA0, __wB0, __wC0, __wD0, __wE0, __wF0;\
		uint64 __wA1, __wB1, __wC1, __wD1, __wE1, __wF1;\
		uint64 __wA2, __wB2, __wC2, __wD2, __wE2, __wF2;\
		uint64 __wA3, __wB3, __wC3, __wD3, __wE3, __wF3;\
		\
		uint64 __cyC0,__cyD0,__cyE0,__cyF0;\
		uint64 __cyC1,__cyD1,__cyE1,__cyF1;\
		uint64 __cyC2,__cyD2,__cyE2,__cyF2;\
		uint64 __cyC3,__cyD3,__cyE3,__cyF3;\
		\
		uint64 __a0,__b0,__e0,__f0,__i0,__j0;\
		uint64 __a1,__b1,__e1,__f1,__i1,__j1;\
		uint64 __a2,__b2,__e2,__f2,__i2,__j2;\
		uint64 __a3,__b3,__e3,__f3,__i3,__j3;\
		\
		SQR_LOHI64(__x0.d0,        __wA0, __wB0);\
		SQR_LOHI64(__x1.d0,        __wA1, __wB1);\
		SQR_LOHI64(__x2.d0,        __wA2, __wB2);\
		SQR_LOHI64(__x3.d0,        __wA3, __wB3);\
		\
		MUL_LOHI64(__x0.d0,__x0.d1, __a0, __b0 );\
		MUL_LOHI64(__x1.d0,__x1.d1, __a1, __b1 );\
		MUL_LOHI64(__x2.d0,__x2.d1, __a2, __b2 );\
		MUL_LOHI64(__x3.d0,__x3.d1, __a3, __b3 );\
		\
		SQR_LOHI64(__x0.d1,        __wC0, __wD0);\
		SQR_LOHI64(__x1.d1,        __wC1, __wD1);\
		SQR_LOHI64(__x2.d1,        __wC2, __wD2);\
		SQR_LOHI64(__x3.d1,        __wC3, __wD3);\
		\
		MUL_LOHI64(__x0.d0,__x0.d2, __e0, __f0 );\
		MUL_LOHI64(__x1.d0,__x1.d2, __e1, __f1 );\
		MUL_LOHI64(__x2.d0,__x2.d2, __e2, __f2 );\
		MUL_LOHI64(__x3.d0,__x3.d2, __e3, __f3 );\
		\
		SQR_LOHI64(__x0.d2,        __wE0, __wF0);\
		SQR_LOHI64(__x1.d2,        __wE1, __wF1);\
		SQR_LOHI64(__x2.d2,        __wE2, __wF2);\
		SQR_LOHI64(__x3.d2,        __wE3, __wF3);\
		\
		MUL_LOHI64(__x0.d1,__x0.d2, __i0, __j0 );\
		MUL_LOHI64(__x1.d1,__x1.d2, __i1, __j1 );\
		MUL_LOHI64(__x2.d1,__x2.d2, __i2, __j2 );\
		MUL_LOHI64(__x3.d1,__x3.d2, __i3, __j3 );\
		\
		/* Now add cross terms: */\
		/* Add x0*x1 twice to w1-3: */\
		__wB0 += __a0;	__cyC0  = (__wB0 < __a0);	__wC0 += __b0;	__cyD0  = (__wC0 < __b0);	__wB0 += __a0;	__cyC0 += (__wB0 < __a0);	__wC0 += __b0;	__cyD0 += (__wC0 < __b0);\
		__wB1 += __a1;	__cyC1  = (__wB1 < __a1);	__wC1 += __b1;	__cyD1  = (__wC1 < __b1);	__wB1 += __a1;	__cyC1 += (__wB1 < __a1);	__wC1 += __b1;	__cyD1 += (__wC1 < __b1);\
		__wB2 += __a2;	__cyC2  = (__wB2 < __a2);	__wC2 += __b2;	__cyD2  = (__wC2 < __b2);	__wB2 += __a2;	__cyC2 += (__wB2 < __a2);	__wC2 += __b2;	__cyD2 += (__wC2 < __b2);\
		__wB3 += __a3;	__cyC3  = (__wB3 < __a3);	__wC3 += __b3;	__cyD3  = (__wC3 < __b3);	__wB3 += __a3;	__cyC3 += (__wB3 < __a3);	__wC3 += __b3;	__cyD3 += (__wC3 < __b3);\
		\
		/* Add x0*x2 twice to w2-4: */\
		__wC0 += __e0;	__cyD0 += (__wC0 < __e0);	__wD0 += __f0;	__cyE0  = (__wD0 < __f0);	__wC0 += __e0;	__cyD0 += (__wC0 < __e0);	__wD0 += __f0;	__cyE0 += (__wD0 < __f0);\
		__wC1 += __e1;	__cyD1 += (__wC1 < __e1);	__wD1 += __f1;	__cyE1  = (__wD1 < __f1);	__wC1 += __e1;	__cyD1 += (__wC1 < __e1);	__wD1 += __f1;	__cyE1 += (__wD1 < __f1);\
		__wC2 += __e2;	__cyD2 += (__wC2 < __e2);	__wD2 += __f2;	__cyE2  = (__wD2 < __f2);	__wC2 += __e2;	__cyD2 += (__wC2 < __e2);	__wD2 += __f2;	__cyE2 += (__wD2 < __f2);\
		__wC3 += __e3;	__cyD3 += (__wC3 < __e3);	__wD3 += __f3;	__cyE3  = (__wD3 < __f3);	__wC3 += __e3;	__cyD3 += (__wC3 < __e3);	__wD3 += __f3;	__cyE3 += (__wD3 < __f3);\
		\
		/* Add x1*y2 twice to w2-5: */\
		__wD0 += __i0;	__cyE0 += (__wD0 < __i0);	__wE0 += __j0;	__cyF0  = (__wE0 < __j0);	__wD0 += __i0;	__cyE0 += (__wD0 < __i0);	__wE0 += __j0;	__cyF0 += (__wE0 < __j0);\
		__wD1 += __i1;	__cyE1 += (__wD1 < __i1);	__wE1 += __j1;	__cyF1  = (__wE1 < __j1);	__wD1 += __i1;	__cyE1 += (__wD1 < __i1);	__wE1 += __j1;	__cyF1 += (__wE1 < __j1);\
		__wD2 += __i2;	__cyE2 += (__wD2 < __i2);	__wE2 += __j2;	__cyF2  = (__wE2 < __j2);	__wD2 += __i2;	__cyE2 += (__wD2 < __i2);	__wE2 += __j2;	__cyF2 += (__wE2 < __j2);\
		__wD3 += __i3;	__cyE3 += (__wD3 < __i3);	__wE3 += __j3;	__cyF3  = (__wE3 < __j3);	__wD3 += __i3;	__cyE3 += (__wD3 < __i3);	__wE3 += __j3;	__cyF3 += (__wE3 < __j3);\
		\
		/* Now process carries: */\
		__wC0 += __cyC0;	__cyD0 += (__wC0 < __cyC0);	__wD0 += __cyD0;	__cyE0 += (__wD0 < __cyD0);	__wE0 += __cyE0;	__cyF0 += (__wE0 < __cyE0);	__wF0 += __cyF0;\
		__wC1 += __cyC1;	__cyD1 += (__wC1 < __cyC1);	__wD1 += __cyD1;	__cyE1 += (__wD1 < __cyD1);	__wE1 += __cyE1;	__cyF1 += (__wE1 < __cyE1);	__wF1 += __cyF1;\
		__wC2 += __cyC2;	__cyD2 += (__wC2 < __cyC2);	__wD2 += __cyD2;	__cyE2 += (__wD2 < __cyD2);	__wE2 += __cyE2;	__cyF2 += (__wE2 < __cyE2);	__wF2 += __cyF2;\
		__wC3 += __cyC3;	__cyD3 += (__wC3 < __cyC3);	__wD3 += __cyD3;	__cyE3 += (__wD3 < __cyD3);	__wE3 += __cyE3;	__cyF3 += (__wE3 < __cyE3);	__wF3 += __cyF3;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo0.d0 = __wA0; __lo0.d1 = __wB0; __lo0.d2 = __wC0; __hi0.d0 = __wD0;	__hi0.d1 = __wE0;	__hi0.d2 = __wF0;\
		__lo1.d0 = __wA1; __lo1.d1 = __wB1; __lo1.d2 = __wC1; __hi1.d0 = __wD1;	__hi1.d1 = __wE1;	__hi1.d2 = __wF1;\
		__lo2.d0 = __wA2; __lo2.d1 = __wB2; __lo2.d2 = __wC2; __hi2.d0 = __wD2;	__hi2.d1 = __wE2;	__hi2.d2 = __wF2;\
		__lo3.d0 = __wA3; __lo3.d1 = __wB3; __lo3.d2 = __wC3; __hi3.d0 = __wD3;	__hi3.d1 = __wE3;	__hi3.d2 = __wF3;\
	}
  #endif

  #if 0
    /* 8-operand-pipelined version: */
    #define SQR_LOHI192_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
    	SQR_LOHI192(__x0, __lo0, __hi0);\
    	SQR_LOHI192(__x1, __lo1, __hi1);\
    	SQR_LOHI192(__x2, __lo2, __hi2);\
    	SQR_LOHI192(__x3, __lo3, __hi3);\
    	SQR_LOHI192(__x4, __lo4, __hi4);\
    	SQR_LOHI192(__x5, __lo5, __hi5);\
    	SQR_LOHI192(__x6, __lo6, __hi6);\
    	SQR_LOHI192(__x7, __lo7, __hi7);\
	}
  #else
    #define SQR_LOHI192_q8(\
      __x0, __lo0, __hi0\
    , __x1, __lo1, __hi1\
    , __x2, __lo2, __hi2\
    , __x3, __lo3, __hi3\
    , __x4, __lo4, __hi4\
    , __x5, __lo5, __hi5\
    , __x6, __lo6, __hi6\
    , __x7, __lo7, __hi7)\
    {\
	/* In the output assembly below, we use MUL outputs in the following order: */	\
	/*     w0,w1,w2,a,b    c,d,w4    e,f,w5    */									\
	/* IDEA: rearrange and interleave the MULs to try to take advantage of that? */	\
		uint64 __wA0, __wB0, __wC0, __wD0, __wE0, __wF0;\
		uint64 __wA1, __wB1, __wC1, __wD1, __wE1, __wF1;\
		uint64 __wA2, __wB2, __wC2, __wD2, __wE2, __wF2;\
		uint64 __wA3, __wB3, __wC3, __wD3, __wE3, __wF3;\
		uint64 __wA4, __wB4, __wC4, __wD4, __wE4, __wF4;\
		uint64 __wA5, __wB5, __wC5, __wD5, __wE5, __wF5;\
		uint64 __wA6, __wB6, __wC6, __wD6, __wE6, __wF6;\
		uint64 __wA7, __wB7, __wC7, __wD7, __wE7, __wF7;\
		\
		uint64 __cyC0,__cyD0,__cyE0,__cyF0;\
		uint64 __cyC1,__cyD1,__cyE1,__cyF1;\
		uint64 __cyC2,__cyD2,__cyE2,__cyF2;\
		uint64 __cyC3,__cyD3,__cyE3,__cyF3;\
		uint64 __cyC4,__cyD4,__cyE4,__cyF4;\
		uint64 __cyC5,__cyD5,__cyE5,__cyF5;\
		uint64 __cyC6,__cyD6,__cyE6,__cyF6;\
		uint64 __cyC7,__cyD7,__cyE7,__cyF7;\
		\
		uint64 __a0,__b0,__e0,__f0,__i0,__j0;\
		uint64 __a1,__b1,__e1,__f1,__i1,__j1;\
		uint64 __a2,__b2,__e2,__f2,__i2,__j2;\
		uint64 __a3,__b3,__e3,__f3,__i3,__j3;\
		uint64 __a4,__b4,__e4,__f4,__i4,__j4;\
		uint64 __a5,__b5,__e5,__f5,__i5,__j5;\
		uint64 __a6,__b6,__e6,__f6,__i6,__j6;\
		uint64 __a7,__b7,__e7,__f7,__i7,__j7;\
		\
		SQR_LOHI64(__x0.d0,        __wA0, __wB0);\
		SQR_LOHI64(__x1.d0,        __wA1, __wB1);\
		SQR_LOHI64(__x2.d0,        __wA2, __wB2);\
		SQR_LOHI64(__x3.d0,        __wA3, __wB3);\
		SQR_LOHI64(__x4.d0,        __wA4, __wB4);\
		SQR_LOHI64(__x5.d0,        __wA5, __wB5);\
		SQR_LOHI64(__x6.d0,        __wA6, __wB6);\
		SQR_LOHI64(__x7.d0,        __wA7, __wB7);\
		\
		MUL_LOHI64(__x0.d0,__x0.d1, __a0, __b0 );\
		MUL_LOHI64(__x1.d0,__x1.d1, __a1, __b1 );\
		MUL_LOHI64(__x2.d0,__x2.d1, __a2, __b2 );\
		MUL_LOHI64(__x3.d0,__x3.d1, __a3, __b3 );\
		MUL_LOHI64(__x4.d0,__x4.d1, __a4, __b4 );\
		MUL_LOHI64(__x5.d0,__x5.d1, __a5, __b5 );\
		MUL_LOHI64(__x6.d0,__x6.d1, __a6, __b6 );\
		MUL_LOHI64(__x7.d0,__x7.d1, __a7, __b7 );\
		\
		SQR_LOHI64(__x0.d1,        __wC0, __wD0);\
		SQR_LOHI64(__x1.d1,        __wC1, __wD1);\
		SQR_LOHI64(__x2.d1,        __wC2, __wD2);\
		SQR_LOHI64(__x3.d1,        __wC3, __wD3);\
		SQR_LOHI64(__x4.d1,        __wC4, __wD4);\
		SQR_LOHI64(__x5.d1,        __wC5, __wD5);\
		SQR_LOHI64(__x6.d1,        __wC6, __wD6);\
		SQR_LOHI64(__x7.d1,        __wC7, __wD7);\
		\
		MUL_LOHI64(__x0.d0,__x0.d2, __e0, __f0 );\
		MUL_LOHI64(__x1.d0,__x1.d2, __e1, __f1 );\
		MUL_LOHI64(__x2.d0,__x2.d2, __e2, __f2 );\
		MUL_LOHI64(__x3.d0,__x3.d2, __e3, __f3 );\
		MUL_LOHI64(__x4.d0,__x4.d2, __e4, __f4 );\
		MUL_LOHI64(__x5.d0,__x5.d2, __e5, __f5 );\
		MUL_LOHI64(__x6.d0,__x6.d2, __e6, __f6 );\
		MUL_LOHI64(__x7.d0,__x7.d2, __e7, __f7 );\
		\
		SQR_LOHI64(__x0.d2,        __wE0, __wF0);\
		SQR_LOHI64(__x1.d2,        __wE1, __wF1);\
		SQR_LOHI64(__x2.d2,        __wE2, __wF2);\
		SQR_LOHI64(__x3.d2,        __wE3, __wF3);\
		SQR_LOHI64(__x4.d2,        __wE4, __wF4);\
		SQR_LOHI64(__x5.d2,        __wE5, __wF5);\
		SQR_LOHI64(__x6.d2,        __wE6, __wF6);\
		SQR_LOHI64(__x7.d2,        __wE7, __wF7);\
		\
		MUL_LOHI64(__x0.d1,__x0.d2, __i0, __j0 );\
		MUL_LOHI64(__x1.d1,__x1.d2, __i1, __j1 );\
		MUL_LOHI64(__x2.d1,__x2.d2, __i2, __j2 );\
		MUL_LOHI64(__x3.d1,__x3.d2, __i3, __j3 );\
		MUL_LOHI64(__x4.d1,__x4.d2, __i4, __j4 );\
		MUL_LOHI64(__x5.d1,__x5.d2, __i5, __j5 );\
		MUL_LOHI64(__x6.d1,__x6.d2, __i6, __j6 );\
		MUL_LOHI64(__x7.d1,__x7.d2, __i7, __j7 );\
		\
		/* Now add cross terms: */\
		/* Add x0*x1 twice to w1-3: */\
		__wB0 += __a0;	__cyC0  = (__wB0 < __a0);	__wC0 += __b0;	__cyD0  = (__wC0 < __b0);	__wB0 += __a0;	__cyC0 += (__wB0 < __a0);	__wC0 += __b0;	__cyD0 += (__wC0 < __b0);\
		__wB1 += __a1;	__cyC1  = (__wB1 < __a1);	__wC1 += __b1;	__cyD1  = (__wC1 < __b1);	__wB1 += __a1;	__cyC1 += (__wB1 < __a1);	__wC1 += __b1;	__cyD1 += (__wC1 < __b1);\
		__wB2 += __a2;	__cyC2  = (__wB2 < __a2);	__wC2 += __b2;	__cyD2  = (__wC2 < __b2);	__wB2 += __a2;	__cyC2 += (__wB2 < __a2);	__wC2 += __b2;	__cyD2 += (__wC2 < __b2);\
		__wB3 += __a3;	__cyC3  = (__wB3 < __a3);	__wC3 += __b3;	__cyD3  = (__wC3 < __b3);	__wB3 += __a3;	__cyC3 += (__wB3 < __a3);	__wC3 += __b3;	__cyD3 += (__wC3 < __b3);\
		__wB4 += __a4;	__cyC4  = (__wB4 < __a4);	__wC4 += __b4;	__cyD4  = (__wC4 < __b4);	__wB4 += __a4;	__cyC4 += (__wB4 < __a4);	__wC4 += __b4;	__cyD4 += (__wC4 < __b4);\
		__wB5 += __a5;	__cyC5  = (__wB5 < __a5);	__wC5 += __b5;	__cyD5  = (__wC5 < __b5);	__wB5 += __a5;	__cyC5 += (__wB5 < __a5);	__wC5 += __b5;	__cyD5 += (__wC5 < __b5);\
		__wB6 += __a6;	__cyC6  = (__wB6 < __a6);	__wC6 += __b6;	__cyD6  = (__wC6 < __b6);	__wB6 += __a6;	__cyC6 += (__wB6 < __a6);	__wC6 += __b6;	__cyD6 += (__wC6 < __b6);\
		__wB7 += __a7;	__cyC7  = (__wB7 < __a7);	__wC7 += __b7;	__cyD7  = (__wC7 < __b7);	__wB7 += __a7;	__cyC7 += (__wB7 < __a7);	__wC7 += __b7;	__cyD7 += (__wC7 < __b7);\
		\
		/* Add x0*x2 twice to w2-4: */\
		__wC0 += __e0;	__cyD0 += (__wC0 < __e0);	__wD0 += __f0;	__cyE0  = (__wD0 < __f0);	__wC0 += __e0;	__cyD0 += (__wC0 < __e0);	__wD0 += __f0;	__cyE0 += (__wD0 < __f0);\
		__wC1 += __e1;	__cyD1 += (__wC1 < __e1);	__wD1 += __f1;	__cyE1  = (__wD1 < __f1);	__wC1 += __e1;	__cyD1 += (__wC1 < __e1);	__wD1 += __f1;	__cyE1 += (__wD1 < __f1);\
		__wC2 += __e2;	__cyD2 += (__wC2 < __e2);	__wD2 += __f2;	__cyE2  = (__wD2 < __f2);	__wC2 += __e2;	__cyD2 += (__wC2 < __e2);	__wD2 += __f2;	__cyE2 += (__wD2 < __f2);\
		__wC3 += __e3;	__cyD3 += (__wC3 < __e3);	__wD3 += __f3;	__cyE3  = (__wD3 < __f3);	__wC3 += __e3;	__cyD3 += (__wC3 < __e3);	__wD3 += __f3;	__cyE3 += (__wD3 < __f3);\
		__wC4 += __e4;	__cyD4 += (__wC4 < __e4);	__wD4 += __f4;	__cyE4  = (__wD4 < __f4);	__wC4 += __e4;	__cyD4 += (__wC4 < __e4);	__wD4 += __f4;	__cyE4 += (__wD4 < __f4);\
		__wC5 += __e5;	__cyD5 += (__wC5 < __e5);	__wD5 += __f5;	__cyE5  = (__wD5 < __f5);	__wC5 += __e5;	__cyD5 += (__wC5 < __e5);	__wD5 += __f5;	__cyE5 += (__wD5 < __f5);\
		__wC6 += __e6;	__cyD6 += (__wC6 < __e6);	__wD6 += __f6;	__cyE6  = (__wD6 < __f6);	__wC6 += __e6;	__cyD6 += (__wC6 < __e6);	__wD6 += __f6;	__cyE6 += (__wD6 < __f6);\
		__wC7 += __e7;	__cyD7 += (__wC7 < __e7);	__wD7 += __f7;	__cyE7  = (__wD7 < __f7);	__wC7 += __e7;	__cyD7 += (__wC7 < __e7);	__wD7 += __f7;	__cyE7 += (__wD7 < __f7);\
		\
		/* Add x1*y2 twice to w3-5: */\
		__wD0 += __i0;	__cyE0 += (__wD0 < __i0);	__wE0 += __j0;	__cyF0  = (__wE0 < __j0);	__wD0 += __i0;	__cyE0 += (__wD0 < __i0);	__wE0 += __j0;	__cyF0 += (__wE0 < __j0);\
		__wD1 += __i1;	__cyE1 += (__wD1 < __i1);	__wE1 += __j1;	__cyF1  = (__wE1 < __j1);	__wD1 += __i1;	__cyE1 += (__wD1 < __i1);	__wE1 += __j1;	__cyF1 += (__wE1 < __j1);\
		__wD2 += __i2;	__cyE2 += (__wD2 < __i2);	__wE2 += __j2;	__cyF2  = (__wE2 < __j2);	__wD2 += __i2;	__cyE2 += (__wD2 < __i2);	__wE2 += __j2;	__cyF2 += (__wE2 < __j2);\
		__wD3 += __i3;	__cyE3 += (__wD3 < __i3);	__wE3 += __j3;	__cyF3  = (__wE3 < __j3);	__wD3 += __i3;	__cyE3 += (__wD3 < __i3);	__wE3 += __j3;	__cyF3 += (__wE3 < __j3);\
		__wD4 += __i4;	__cyE4 += (__wD4 < __i4);	__wE4 += __j4;	__cyF4  = (__wE4 < __j4);	__wD4 += __i4;	__cyE4 += (__wD4 < __i4);	__wE4 += __j4;	__cyF4 += (__wE4 < __j4);\
		__wD5 += __i5;	__cyE5 += (__wD5 < __i5);	__wE5 += __j5;	__cyF5  = (__wE5 < __j5);	__wD5 += __i5;	__cyE5 += (__wD5 < __i5);	__wE5 += __j5;	__cyF5 += (__wE5 < __j5);\
		__wD6 += __i6;	__cyE6 += (__wD6 < __i6);	__wE6 += __j6;	__cyF6  = (__wE6 < __j6);	__wD6 += __i6;	__cyE6 += (__wD6 < __i6);	__wE6 += __j6;	__cyF6 += (__wE6 < __j6);\
		__wD7 += __i7;	__cyE7 += (__wD7 < __i7);	__wE7 += __j7;	__cyF7  = (__wE7 < __j7);	__wD7 += __i7;	__cyE7 += (__wD7 < __i7);	__wE7 += __j7;	__cyF7 += (__wE7 < __j7);\
		\
		/* Now process carries: */\
		__wC0 += __cyC0;	__cyD0 += (__wC0 < __cyC0);	__wD0 += __cyD0;	__cyE0 += (__wD0 < __cyD0);	__wE0 += __cyE0;	__cyF0 += (__wE0 < __cyE0);	__wF0 += __cyF0;\
		__wC1 += __cyC1;	__cyD1 += (__wC1 < __cyC1);	__wD1 += __cyD1;	__cyE1 += (__wD1 < __cyD1);	__wE1 += __cyE1;	__cyF1 += (__wE1 < __cyE1);	__wF1 += __cyF1;\
		__wC2 += __cyC2;	__cyD2 += (__wC2 < __cyC2);	__wD2 += __cyD2;	__cyE2 += (__wD2 < __cyD2);	__wE2 += __cyE2;	__cyF2 += (__wE2 < __cyE2);	__wF2 += __cyF2;\
		__wC3 += __cyC3;	__cyD3 += (__wC3 < __cyC3);	__wD3 += __cyD3;	__cyE3 += (__wD3 < __cyD3);	__wE3 += __cyE3;	__cyF3 += (__wE3 < __cyE3);	__wF3 += __cyF3;\
		__wC4 += __cyC4;	__cyD4 += (__wC4 < __cyC4);	__wD4 += __cyD4;	__cyE4 += (__wD4 < __cyD4);	__wE4 += __cyE4;	__cyF4 += (__wE4 < __cyE4);	__wF4 += __cyF4;\
		__wC5 += __cyC5;	__cyD5 += (__wC5 < __cyC5);	__wD5 += __cyD5;	__cyE5 += (__wD5 < __cyD5);	__wE5 += __cyE5;	__cyF5 += (__wE5 < __cyE5);	__wF5 += __cyF5;\
		__wC6 += __cyC6;	__cyD6 += (__wC6 < __cyC6);	__wD6 += __cyD6;	__cyE6 += (__wD6 < __cyD6);	__wE6 += __cyE6;	__cyF6 += (__wE6 < __cyE6);	__wF6 += __cyF6;\
		__wC7 += __cyC7;	__cyD7 += (__wC7 < __cyC7);	__wD7 += __cyD7;	__cyE7 += (__wD7 < __cyD7);	__wE7 += __cyE7;	__cyF7 += (__wE7 < __cyE7);	__wF7 += __cyF7;\
		\
		/* Now split the result between __lo and __hi: */\
		__lo0.d0 = __wA0; __lo0.d1 = __wB0; __lo0.d2 = __wC0; __hi0.d0 = __wD0;	__hi0.d1 = __wE0;	__hi0.d2 = __wF0;\
		__lo1.d0 = __wA1; __lo1.d1 = __wB1; __lo1.d2 = __wC1; __hi1.d0 = __wD1;	__hi1.d1 = __wE1;	__hi1.d2 = __wF1;\
		__lo2.d0 = __wA2; __lo2.d1 = __wB2; __lo2.d2 = __wC2; __hi2.d0 = __wD2;	__hi2.d1 = __wE2;	__hi2.d2 = __wF2;\
		__lo3.d0 = __wA3; __lo3.d1 = __wB3; __lo3.d2 = __wC3; __hi3.d0 = __wD3;	__hi3.d1 = __wE3;	__hi3.d2 = __wF3;\
		__lo4.d0 = __wA4; __lo4.d1 = __wB4; __lo4.d2 = __wC4; __hi4.d0 = __wD4;	__hi4.d1 = __wE4;	__hi4.d2 = __wF4;\
		__lo5.d0 = __wA5; __lo5.d1 = __wB5; __lo5.d2 = __wC5; __hi5.d0 = __wD5;	__hi5.d1 = __wE5;	__hi5.d2 = __wF5;\
		__lo6.d0 = __wA6; __lo6.d1 = __wB6; __lo6.d2 = __wC6; __hi6.d0 = __wD6;	__hi6.d1 = __wE6;	__hi6.d2 = __wF6;\
		__lo7.d0 = __wA7; __lo7.d1 = __wB7; __lo7.d2 = __wC7; __hi7.d0 = __wD7;	__hi7.d1 = __wE7;	__hi7.d2 = __wF7;\
	}
  #endif

#endif

/* Lower 192 bits of the product of uint192 x and y.

Let x = x0 + x1*2^64 + x2*2^128,  with x0, x1, x2 < 2^64,
and y = y0 + y1*2^64 + y2*2^128,  with y0, y1, y2 < 2^64.            Then,

 x*y = x0*y0 + (x0*y1 + x1*y0)*2^64 + (x1*y1 + x2*y0)*2^128 + (x2*y1 + x1*y2)*2^192 + (x2*y2)*2^256 .

In terms of the 5 output coefficients, here is what goes into each of those:

	w0 = (x0*y0).lo
	w1 = (x0*y0).hi + (x0*y1 + x1*y0).lo
	w2 =              (x0*y1 + x1*y0).hi + (x0*y2 + x1*y1 + x2*y0).lo
	w3 =                                   (x0*y2 + x1*y1 + x2*y0).hi + (x2*y1 + x1*y2).lo
	w4 =                                                                (x2*y1 + x1*y2).hi + (x2*y2).lo
	w5 =                                                                                   + (x2*y2).hi .

In this case, we only need the items that go into w0-2.
On Alpha, this needs a total of  9 MUL instructions and 10 ALU ops.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULL192(__x, __y, __lo)\
    {\
		uint64 __w0,__w1,__w2,__a,__b,__c,__d,__e,__g;\
		\
		MUL_LOHI64(__x.d0,__y.d0,&__w0,&__w1);	/*   x0*y0 */\
		__w2 =   __x.d1*__y.d1            ;	/*   x1*y1 */\
		\
		MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0,&__c ,&__d );	/*   x1*y0 */\
		\
		__e  =   __x.d0*__y.d2            ;	/*   x0*y2 */\
		__g  =   __x.d2*__y.d0            ;	/*   x2*y0 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__w2 += (__w1 < __a);\
		__w2 += __b;\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__w2 += (__w1 < __c);\
		__w2 += __d;\
		\
		/* Add x0*y2 to w2-4: */\
		__w2 += __e;\
		\
		/* Add x2*y0 to w2-4: */\
		__w2 += __g;\
		\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2;\
    }

    /* 4-operand-pipelined version: */
    #define MULL192_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
    	MULL192(__x0, __y0, __lo0);\
    	MULL192(__x1, __y1, __lo1);\
    	MULL192(__x2, __y2, __lo2);\
    	MULL192(__x3, __y3, __lo3);\
	}

    /* 8-operand-pipelined version: */
    #define MULL192_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
    	MULL192(__x0, __y0, __lo0);\
    	MULL192(__x1, __y1, __lo1);\
    	MULL192(__x2, __y2, __lo2);\
    	MULL192(__x3, __y3, __lo3);\
    	MULL192(__x4, __y4, __lo4);\
    	MULL192(__x5, __y5, __lo5);\
    	MULL192(__x6, __y6, __lo6);\
    	MULL192(__x7, __y7, __lo7);\
	}

#else

  #if defined(YES_ASM)

	/* x86_64 ASM implementation of the 192-bit MUL macros: */
    #define MULL192(__x, __y, __lo)\
    {\
	__asm__ volatile (\
	"movq	%[__x0],%%rax	\n\t	mulq	%[__y0]	\n\t	movq	%%rax,%%r10	\n\t	movq	%%rdx,%%r11	\n\t"/* w0:w1 = x0*y0\ */\
	"movq	%[__x1],%%r12	\n\t	imulq	%[__y1],%%r12	\n\t"/* w2 = (x1*y1).lo */\
	"movq	%[__x0],%%rax	\n\t	mulq	%[__y1]	\n\t"/* _a:_b = x0*y1 */\
		"addq	%%rax,%%r11		\n\t"/* __w1 += __a */\
		"adcq	%%rdx,%%r12		\n\t"/* __w2 += __b + (__w1 < __a) */\
	"movq	%[__x1],%%rax	\n\t	mulq	%[__y0]	\n\t"/* _c:_d = x1*y0 */\
		"addq	%%rax,%%r11		\n\t"/* __w1 += __c */\
		"adcq	%%rdx,%%r12		\n\t"/* __w2 += __d + (__w1 < __c) */\
	"movq	%%r11,%[__lo1]	\n\t"/* __lo.d1 = __w1 */\
	"movq	%[__x0],%%rdi	\n\t	imulq	%[__y2],%%rdi	\n\t"/* _e = (x0*y2).lo */\
		"addq	%%rdi,%%r12		\n\t"/* __w2 += __e */\
	"movq	%[__x2],%%rdi	\n\t	imulq	%[__y0],%%rdi	\n\t"/* _g = (x2*y0).lo */\
		"addq	%%rdi,%%r12		\n\t"/* __w2 += __g */\
	"movq	%%r10,%[__lo0]	\n\t"/* __lo.d0 = __w0; must wait until end because MULL designed for in-place capability, i.e. lo0 may overlap x0 or y0 */\
	"movq	%%r12,%[__lo2]	\n\t"/* __lo.d2 = __w2 */\
		: [__lo0] "=m" (__lo.d0)	/* Outputs are 3 low words of 6-word product */\
		 ,[__lo1] "=m" (__lo.d1)	\
		 ,[__lo2] "=m" (__lo.d2)	\
		: [__x0] "m" (__x.d0)	/* Inputs from memory/register here */\
		 ,[__x1] "m" (__x.d1)	\
		 ,[__x2] "m" (__x.d2)	\
		 ,[__y0] "m" (__y.d0)	\
		 ,[__y1] "m" (__y.d1)	\
		 ,[__y2] "m" (__y.d2)	\
		: "cc","memory","rax","rdx","rdi","r10","r11","r12"	/* Clobbered registers [all but rcx,rsi,r8-9,r13-15] */\
	);\
	}

  #else

    #define MULL192(__x, __y, __lo)\
    {\
		uint64 __w0,__w1,__w2,__a,__b,__c,__d,__e,__g;\
		\
		MUL_LOHI64(__x.d0,__y.d0, __w0, __w1);	/*   x0*y0 - Use a tmp for low output word to allow in-placeness */\
		__w2 =   __x.d1*__y.d1            ;	/*   x1*y1 */\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0, __c , __d );	/*   x1*y0 */\
		__e  =   __x.d0*__y.d2            ;	/*   x0*y2 */\
		__g  =   __x.d2*__y.d0            ;	/*   x2*y0 */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__w2 += (__w1 < __a) + __b;	/* Add CY + __b directly into w2 since __b is a UMULH output, i.e. (CY + __b) < 2^64 */\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__w2 += (__w1 < __c) + __d;\
		/* Add mull(x0,y2) + mull(x2,y0)to w2: */\
		__lo.d0 = __w0; __lo.d1 = __w1; __lo.d2 = __w2 + __e + __g;\
    }

  #endif // YES_ASM

    /* 4-operand-pipelined version: */
    #define MULL192_q4(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3)\
    {\
    	MULL192(__x0, __y0, __lo0);\
    	MULL192(__x1, __y1, __lo1);\
    	MULL192(__x2, __y2, __lo2);\
    	MULL192(__x3, __y3, __lo3);\
	}

    /* 8-operand-pipelined version: */
    #define MULL192_q8(\
      __x0, __y0, __lo0\
    , __x1, __y1, __lo1\
    , __x2, __y2, __lo2\
    , __x3, __y3, __lo3\
    , __x4, __y4, __lo4\
    , __x5, __y5, __lo5\
    , __x6, __y6, __lo6\
    , __x7, __y7, __lo7)\
    {\
    	MULL192(__x0, __y0, __lo0);\
    	MULL192(__x1, __y1, __lo1);\
    	MULL192(__x2, __y2, __lo2);\
    	MULL192(__x3, __y3, __lo3);\
    	MULL192(__x4, __y4, __lo4);\
    	MULL192(__x5, __y5, __lo5);\
    	MULL192(__x6, __y6, __lo6);\
    	MULL192(__x7, __y7, __lo7);\
	}

#endif

/* Upper 192 bits of the product of uint192 x and y.

Let x = x0 + x1*2^64 + x2*2^128,  with x0, x1, x2 < 2^64,
and y = y0 + y1*2^64 + y2*2^128,  with y0, y1, y2 < 2^64.

Then,

 x*y = x0*y0
     + (x0*y1 + x1*y0)*2^64
     + (x1*y1 + x2*y0)*2^128
     + (x2*y1 + x1*y2)*2^192
     + (x2*y2)*2^256 .

In terms of the 5 output coefficients, here is what goes into each of those:

	w0 = (x0*y0).lo
	w1 = (x0*y0).hi + (x0*y1 + x1*y0).lo
	w2 =              (x0*y1 + x1*y0).hi + (x0*y2 + x1*y1 + x2*y0).lo
	w3 =                                   (x0*y2 + x1*y1 + x2*y0).hi + (x2*y1 + x1*y2).lo
	w4 =                                                                (x2*y1 + x1*y2).hi + (x2*y2).lo
	w5 =                                                                                   + (x2*y2).hi .

On Alpha, this needs a total of:

	17 MUL, 36 ALU ops - if we require a 100% correct computation of the carry into w3
	                     from the lower-order terms;

	14 MUL, 30 ALU ops - if we can afford a probabilistic computation of carry into w3
	                     from the lower-order terms. (Savings probably not worth it).

! 01/23/2004: I thought I would generally be able to neglect all the w0 and w1 terms,
! but I almost immediately ran into the following case: which occurs during the mod-inverse
! computation preceding the modpow loop used in checking the known factor of M1283:
!
! X =                  14*2^128 +  3291757557782450881*2^64 +  3893270457587058239
! Y = 2610710867498977111*2^128 + 18111005741868530922*2^64 + 10462639744819210687 .
!
! This occurs during the following sequence just prior to entering the for-loop:
!
!	LSHIFT192(qinv, zshift, lo);
!	MULH192(q,lo,lo);
!
! In the above, Y contains the left-shifted version of qinv, and X contains the factor q.
! By the definition of qinv, q*(qinv*2^zshift) yields a result whose lower (192 + zshift) bits
! are identically zero, except for a lone ones bit in the (zshift)th position. That means that
! in calculating this we expect to see instances where adding two 64-bit operands will yield 2^64,
! i.e. all zeros in the current word, with a carry into the next-higher word. We thus need to be
! really careful in calculating carries - can't take shortcuts and ignore the w0 and w1 outputs.
! For the above example, if neglect w0 and w1, wind up with w2 = 18446744073709551615 = 2^64 - 1
! and carry into the upper half = 1, rather than the true w2 = 0 and carry = 2.

n = 147
n++; p = 2^n-1; for (i = 1, 1000, if (isprime(p+2*i), print(p+2*i)))

./Mfactor -kmax 1000000 -m

lg p	p											factor k	pass
----	----------------------------------------	--------	----
128			 340282366920938463463374607431768211841	       4	1
129			 680564733841876926926749214863536422987	 6555772	13
130			1361129467683753853853498429727072846803	       1	0
131			2722258935367507707706996859454145692553	    1767	7
132			5444517870735015415413993718908291383419	     297	14
133		   10889035741470030830827987437816582767771	420,9904	15,1
134		   21778071482940061661655974875633165533913	      47	11
135		   43556142965880123323311949751266331067181	       4	1
136		   87112285931760246646623899502532662134269	 2657592	1
137		  174224571863520493293247799005065324267283	 1172021	9
138		  348449143727040986586495598010130648531089	  255904	0
139		  696898287454081973172991196020261297062043	 2901633	7
140		 1393796574908163946345982392040522594124799	 2532965	0
141		 2787593149816327892691964784081045188248599	  124200	15
142		 5575186299632655785383929568162090376498137	       3	0
143		11150372599265311570767859136324180752990493	52567,31275	1,3
144		22300745198530623141535718272648361505982037	   24063	0
145		44601490397061246283071436545296723011960969	     440	3
146		89202980794122492566142873090593446023922217	      39	9
147	   178405961588244985132285746181186892047844053	       3	0

!
! However, this kind of thing is only likely to happen in the mod-inverse computation and the
! q*(qinv << shift) multiply that we do just prior to entering the modmul-based powering loop, to save
! one full-cost modmul. So we define two versions of the MULH routine: an exact one which will be used
! only for the initial mul, and one which approximates the (presumably quasirandom) lower bits of the mul
! so as to get the proper carry into the upper half with very high probability, which we use everywhere else.
! **** The fast variant should only ever be used with spot-checking enabled *****

08/13/2012: Addendum: In the context of the Montgomery-mul mod-inverse computation we know the
low-half product terms should sum to unity, that is, the above w2 "carry layer" terms, if computed
exactly with carryins from the lower-order products, should sum to 0 (mod 2^64). Thus since the
truncated-carry-layer approximation gives w2 = 2^64 -1, we know that we missed a carry here, and
simply add 1 to the approximate-computed carry to obtain the correct result.
*/
#ifdef MUL_LOHI64_SUBROUTINE

    #define MULH192(__x, __y, __hi)\
    {\
		uint64 __w1, __w2, __w3, __w4, __w5\
		/* Carries are numbered in the sense of "__cyX = carry OUT OF word __wX": */\
					,__cy2,__cy3,__cy4,__cy5\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		MULH64(__x.d0,__y.d0,       __w1);	/*   x0*y0.hi */\
		MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);	/*   x1*y1 */\
		MUL_LOHI64(__x.d2,__y.d2,&__w4,&__w5);	/*   x2*y2 */\
		\
		MUL_LOHI64(__x.d0,__y.d1,&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0,&__c ,&__d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2,&__e ,&__f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0,&__g ,&__h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x.d1,__y.d2,&__i ,&__j );	/*   x1*y2 */\
		MUL_LOHI64(__x.d2,__y.d1,&__k ,&__l );	/*   x2*y1 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-4: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-4: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x1*y2 to w3-5: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x2*y1 to w3-5: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		__hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

  #if 0

    #define MULH192_FAST(__x, __y, __hi)	MULH192(__x, __y, __hi)\

  #elif defined(YES_ASM)

    #define MULH192_FAST(__x, __y, __hi)\
    {\
		uint64 __w2,__w3, __w4, __w5,__cy3,__cy4,__cy5,__bd_lo,__e,__f,__g,__h,__i,__j,__k,__l;\
		double __fprod, __scale = TWO64FLINV*TWO64FLINV;\
		\
		/* Compute (x0*y1 + x1*y0)>>64 via floating-double approximation: */\
		__fprod = ((double)__x.d0*(double)__y.d1 + (double)__y.d0*(double)__x.d1)*__scale;\
		MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);	/*   x1*y1 */\
		__cy3 = (uint64)__fprod;	/* Contribution to __cy3 of the above floating-point-approximated terms */\
		__fprod = (__fprod - (double)__cy3)*TWO64FLOAT;	/*** To-Do: Speed these float <--> int conversions ***/\
		__bd_lo = (uint64)__fprod;\
		/* Now start adding cross terms: */\
		__w2 += __bd_lo;__cy3 += (__w2 < __bd_lo);\
		MUL_LOHI64(__x.d0,__y.d2,&__e ,&__f );	/*   x0*y2 */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		MUL_LOHI64(__x.d2,__y.d0,&__g ,&__h );	/*   x2*y0 */\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		MUL_LOHI64(__x.d1,__y.d2,&__i ,&__j );	/*   x1*y2 */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		MUL_LOHI64(__x.d2,__y.d2,&__w4,&__w5);	/*   x2*y2 */\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		MUL_LOHI64(__x.d2,__y.d1,&__k ,&__l );	/*   x2*y1 */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		/* Now process last set of carries: */\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		__hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

  #endif

	/* Lower-half-approximated version of carry-approximated MULH192, for use in my bit-doubling fast modular inversion
	scheme for computing the Montgomery-mul mod-inverse. Here we know the carry out of the low-half terms = 0,
	so don't even try to compute it, either exactly or approximately. That is, compute the terms of the multiply
	rhombus, with zero carryin to w2, and caller supplying the expected value of w2 for purpose of error correction:

		w2 = (x0*y1 + x1*y0).hi + (x0*y2 + x1*y1 + x2*y0).lo
		w3 =                      (x0*y2 + x1*y1 + x2*y0).hi + (x2*y1 + x1*y2).lo
		w4 =                                                   (x2*y1 + x1*y2).hi + (x2*y2).lo
		w5 =                                                                      + (x2*y2).hi .

	Use floating-double approximation of the (x0*y1 + x1*y0).hi term; on x86 the (x0*y2 + x1*y1 + x2*y0).lo is essentially
	free because the MUL instruction needed for exact computation of the .hi halves also produces the low halves, we just
	need to add those with carry.

	At this length the truncated-rhombus version doesn't save us much work - 2 slow MUL_LOHI64s get replaced by 5 fast FMULs -
	but is important for proof-of-concept purposes.
	*/
    #define MULH192_TRUNC(__x, __y, __w2_exact, __hi)\
    {\
		uint64 __w2,__w3, __w4, __w5,__cy3,__cy4,__cy5,__bd_lo,__e,__f,__g,__h,__i,__j,__k,__l;\
		double __fprod, __scale = TWO64FLINV*TWO64FLINV;\
		\
		/* Compute (x0*y1 + x1*y0)>>64 via floating-double approximation: */\
		__fprod = ((double)__x.d0*(double)__y.d1 + (double)__y.d0*(double)__x.d1)*__scale;\
		MUL_LOHI64(__x.d1,__y.d1,&__w2,&__w3);	/*   x1*y1 */\
		__cy3 = (uint64)__fprod;	/* Contribution to __cy3 of the above floating-point-approximated terms */\
		__fprod = (__fprod - (double)__cy3)*TWO64FLOAT;	/*** To-Do: Speed these float <--> int conversions ***/\
		__bd_lo = (uint64)__fprod;\
		/* Now start adding cross terms: */\
		__w2 += __bd_lo;__cy3 += (__w2 < __bd_lo);\
		MUL_LOHI64(__x.d0,__y.d2,&__e ,&__f );	/*   x0*y2 */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		MUL_LOHI64(__x.d2,__y.d0,&__g ,&__h );	/*   x2*y0 */\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		MUL_LOHI64(__x.d1,__y.d2,&__i ,&__j );	/*   x1*y2 */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
	/* If result is close to but < [expected value], know that carry-layer approximation dropped a carry: */\
	__cy3 += ((int64)(__w2 - __w2_exact) < 0);\
		\
		MUL_LOHI64(__x.d2,__y.d2,&__w4,&__w5);	/*   x2*y2 */\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		MUL_LOHI64(__x.d2,__y.d1,&__k ,&__l );	/*   x2*y1 */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		/* Now process last set of carries: */\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		__hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
	}

#else

    #define MULH192(__x, __y, __hi)\
    {\
		uint64 __w1, __w2, __w3, __w4, __w5\
					,__cy2,__cy3,__cy4,__cy5\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		MULH64(__x.d0,__y.d0,       __w1);	/*   x0*y0.hi */\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);	/*   x1*y1 */\
		MUL_LOHI64(__x.d2,__y.d2, __w4, __w5);	/*   x2*y2 */\
		\
		MUL_LOHI64(__x.d0,__y.d1, __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x.d1,__y.d0, __c , __d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x.d0,__y.d2, __e , __f );	/*   x0*y2 */\
		MUL_LOHI64(__x.d2,__y.d0, __g , __h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x.d1,__y.d2, __i , __j );	/*   x1*y2 */\
		MUL_LOHI64(__x.d2,__y.d1, __k , __l );	/*   x2*y1 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-3: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-3: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-4: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-4: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x1*y2 to w3-5: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x2*y1 to w3-5: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		__hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

  #if 0

    #define MULH192_FAST(__x, __y, __hi)	MULH192(__x, __y, __hi)

  #elif defined(YES_ASM)

   #if 1
    #define MULH192_FAST(__x, __y, __hi)	MULH192(__x, __y, __hi)
   #else
	/* x86_64 ASM implementation of the 192-bit MULH macro: */
	#error Using x86_64 asm version of MULH192_FAST - needs debug!
    #define MULH192_FAST(__x, __y, __hi)\
    {\
	__asm__ volatile (\
		"xorq   %%rbx, %%rbx    \n\t"/* Tie rbx = 0 ... the 'b' is for 'bupkis' :) */\
		"xorq	%%r8 ,%%r8 		\n\t"/* cy3 = 0 */\
		"xorq	%%r9 ,%%r9 		\n\t"/* cy4 = 0 */\
		/* Compute (x0*y1 + x1*y0): */\
	"movq	%[__x0],%%rax	\n\t	mulq	%[__y1]	\n\t	movq	%%rdx,%%rdi	\n\t"/* rdi = MULH64(x0,y1) */\
	"movq	%[__x1],%%rax	\n\t	mulq	%[__y0]	\n\t	movq	%%rdx,%%rsi	\n\t"/* rsi = MULH64(x1,y0) */\
		"addq	%%rsi,%%rdi		\n\t"/* bd_lo in rdi */\
		"adcq	%%rbx,%%r8 		\n\t"/* cy3 = (x0*y1 + x1*y0).hi */\
	"movq	%[__x1],%%rax	\n\t	mulq	%[__y1]	\n\t	movq %%rax,%%r12 \n\t movq %%rdx,%%r13	\n\t"/* w2:w3 = x1*y1 */\
		/* Now start adding cross terms: */\
		"addq	%%rdi,%%r12		\n\t"/* __w2 += __bd_lo */\
		"adcq	%%rbx,%%r8 		\n\t"/* __cy3 += (__w2 < __bd_lo) */\
	"movq	%[__x0],%%rax	\n\t	mulq	%[__y2]	\n\t"/* _e:_f = x0*y2 */\
		"addq	%%rax,%%r12		\n\t"/* __w2 += __e */\
		"adcq	%%rbx,%%r8 		\n\t"/* __cy3 += (__w2 < __e) */\
		"addq	%%rdx,%%r13		\n\t"/* __w3 += __f */\
		"adcq	%%rbx,%%r9 		\n\t"/* __cy4 += (__w3 < __f) */\
	"movq	%[__x2],%%rax	\n\t	mulq	%[__y0]	\n\t"/* _g:_h = x2*y0 */\
		"addq	%%rax,%%r12		\n\t"/* __w2 += __g */\
		"adcq	%%rbx,%%r8 		\n\t"/* __cy3 += (__w2 < __g) */\
		"addq	%%rdx,%%r13		\n\t"/* __w3 += __h */\
		"adcq	%%rbx,%%r9 		\n\t"/* __cy4 += (__w3 < __h) */\
	"movq	%[__x2],%%rax	\n\t	mulq	%[__y2]	\n\t	movq	%%rax,%%r14	\n\t	movq	%%rdx,%%r15	\n\t"/* w4:w5 = x2*y2 */\
	"movq	%[__x1],%%rax	\n\t	mulq	%[__y2]	\n\t"/* _i:_j = x1*y2 */\
		"addq	%%rax,%%r13		\n\t"/* __w3 += __i */\
		"adcq	%%rbx,%%r9 		\n\t"/* __cy4 += (__w3 < __i) */\
		"addq	%%rdx,%%r14		\n\t"/* __w4 += __j */\
		"adcq	%%rbx,%%r15		\n\t"/* __cy5 += (__w4 < __j); simply add directly into w5 */\
	"movq	%[__x2],%%rax	\n\t	mulq	%[__y1]	\n\t"/* _k:_l = x2*y1 */\
		"addq	%%rax,%%r13		\n\t"/* __w3 += __k */\
		"adcq	%%rbx,%%r9 		\n\t"/* __cy4 += (__w3 < __k) */\
		"addq	%%rdx,%%r14		\n\t"/* __w4 += __l */\
		"adcq	%%rbx,%%r15		\n\t"/* __cy5 += (__w4 < __l); simply add directly into w5 */\
		/* Now process last set of carries: */\
		"addq	%%r8 ,%%r13		\n\t"/* __w3 += __cy3 */\
		"adcq	%%rbx,%%r9 		\n\t"/* __cy4 += (__w3 < __cy3) */\
		"addq	%%r9 ,%%r14		\n\t"/* __w4 += __cy4 */\
		"adcq	%%rbx,%%r15		\n\t"/* __cy5 += (__w4 < __cy4); simply add directly into w5 */\
		\
		"movq	%%r13,%[__hi0]	\n\t"/* __hi.d0 = __w3 */\
		"movq	%%r14,%[__hi1]	\n\t"/* __hi.d1 = __w4 */\
		"movq	%%r15,%[__hi2]	\n\t"/* __hi.d2 = __w5 */\
		: /* outputs: none */\
		: [__x0] "g" (__x.d0)	/* All inputs from memory/register here */\
		 ,[__x1] "g" (__x.d1)	\
		 ,[__x2] "g" (__x.d2)	\
		 ,[__y0] "g" (__y.d0)	\
		 ,[__y1] "g" (__y.d1)	\
		 ,[__y2] "g" (__y.d2)	\
		 ,[__hi0] "g" (__hi.d0)	\
		 ,[__hi1] "g" (__hi.d1)	\
		 ,[__hi2] "g" (__hi.d2)	\
		: "cc","memory","rax","rbx","rdx","rdi","rsi","r8","r9","r12","r13","r14","r15"	/* Clobbered registers [all but rcx,rsi,r10,r11] */\
	);\
	}
   #endif	// if(0|1)

  #else

    #define MULH192_FAST(__x, __y, __hi)\
    {\
		uint64 __w2,__w3, __w4, __w5,__cy3,__cy4,__cy5,__bd_lo,__e,__f,__g,__h,__i,__j,__k,__l;\
		double __fprod, __scale = TWO64FLINV*TWO64FLINV;\
		\
		/* Compute (x0*y1 + x1*y0)>>64 via floating-double approximation: */\
		__fprod = ((double)__x.d0*(double)__y.d1 + (double)__y.d0*(double)__x.d1)*__scale;\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);	/*   x1*y1 */\
		__cy3 = (uint64)__fprod;	/* Contribution to __cy3 of the above floating-point-approximated terms */\
		__fprod = (__fprod - (double)__cy3)*TWO64FLOAT;	/*** To-Do: Speed these float <--> int conversions ***/\
		__bd_lo = (uint64)__fprod;\
		/* Now start adding cross terms: */\
		__w2 += __bd_lo;__cy3 += (__w2 < __bd_lo);\
		MUL_LOHI64(__x.d0,__y.d2, __e , __f );	/*   x0*y2 */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		MUL_LOHI64(__x.d2,__y.d0, __g , __h );	/*   x2*y0 */\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		MUL_LOHI64(__x.d1,__y.d2, __i , __j );	/*   x1*y2 */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		MUL_LOHI64(__x.d2,__y.d2, __w4, __w5);	/*   x2*y2 */\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		MUL_LOHI64(__x.d2,__y.d1, __k , __l );	/*   x2*y1 */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		/* Now process last set of carries: */\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		__hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

#endif

    #define MULH192_TRUNC(__x, __y, __w2_exact, __hi)\
    {\
		uint64 __w2,__w3, __w4, __w5,__cy3,__cy4,__cy5,__bd_lo,__e,__f,__g,__h,__i,__j,__k,__l;\
		double __fprod, __scale = TWO64FLINV*TWO64FLINV;\
		\
		/* Compute (x0*y1 + x1*y0)>>64 via floating-double approximation: */\
		__fprod = ((double)__x.d0*(double)__y.d1 + (double)__y.d0*(double)__x.d1)*__scale;\
		MUL_LOHI64(__x.d1,__y.d1, __w2, __w3);	/*   x1*y1 */\
		__cy3 = (uint64)__fprod;	/* Contribution to __cy3 of the above floating-point-approximated terms */\
		__fprod = (__fprod - (double)__cy3)*TWO64FLOAT;	/*** To-Do: Speed these float <--> int conversions ***/\
		__bd_lo = (uint64)__fprod;\
		/* Now start adding cross terms: */\
		__w2 += __bd_lo;__cy3 += (__w2 < __bd_lo);\
		MUL_LOHI64(__x.d0,__y.d2, __e , __f );	/*   x0*y2 */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		MUL_LOHI64(__x.d2,__y.d0, __g , __h );	/*   x2*y0 */\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		MUL_LOHI64(__x.d1,__y.d2, __i , __j );	/*   x1*y2 */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
	/* If result is close to but < [expected value], know that carry-layer approximation dropped a carry: */\
	__cy3 += ((int64)(__w2 - __w2_exact) < 0);\
		\
		MUL_LOHI64(__x.d2,__y.d2, __w4, __w5);	/*   x2*y2 */\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		MUL_LOHI64(__x.d2,__y.d1, __k , __l );	/*   x2*y1 */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		/* Now process last set of carries: */\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;\
		\
		__hi.d0 = __w3;	__hi.d1 = __w4;	__hi.d2 = __w5;\
    }

#endif

    /* 4-operand-pipelined version: */
    #define MULH192_q4(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3)\
    {\
    	MULH192_FAST(__x0, __y0, __hi0);\
    	MULH192_FAST(__x1, __y1, __hi1);\
    	MULH192_FAST(__x2, __y2, __hi2);\
    	MULH192_FAST(__x3, __y3, __hi3);\
	}

    /* 8-operand-pipelined version: */
    #define MULH192_q8(\
      __x0, __y0, __hi0\
    , __x1, __y1, __hi1\
    , __x2, __y2, __hi2\
    , __x3, __y3, __hi3\
    , __x4, __y4, __hi4\
    , __x5, __y5, __hi5\
    , __x6, __y6, __hi6\
    , __x7, __y7, __hi7)\
    {\
    	MULH192_FAST(__x0, __y0, __hi0);\
    	MULH192_FAST(__x1, __y1, __hi1);\
    	MULH192_FAST(__x2, __y2, __hi2);\
    	MULH192_FAST(__x3, __y3, __hi3);\
    	MULH192_FAST(__x4, __y4, __hi4);\
    	MULH192_FAST(__x5, __y5, __hi5);\
    	MULH192_FAST(__x6, __y6, __hi6);\
    	MULH192_FAST(__x7, __y7, __hi7);\
	}

#ifdef __cplusplus
}
#endif

#endif	/* imul_macro1_h_included */

