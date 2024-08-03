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
#ifndef types_h_included
#define types_h_included

/* Include any needed level-0 header files: */
#include "platform.h"
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Typedefs */

/*...useful utility parameters */

#undef TRUE
#define TRUE	true

#undef FALSE
#define FALSE	false

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

typedef          int8_t		 int8;
typedef          int8_t		sint8;
typedef uint8_t			uint8;

typedef          int16_t	 int16;
typedef          int16_t	sint16;
typedef uint16_t		uint16;

typedef          int32_t	 int32;
typedef          int32_t	sint32;
typedef uint32_t		uint32;

typedef          int64_t	 int64;
typedef          int64_t	sint64;
typedef uint64_t		uint64;

/*
#ifdef int32_t
	#warning int32_t already defined!
#else
	typedef  int32		 int32_t;
	typedef uint32		uint32_t;
#endif

EWM: nvcc gives 'error: invalid redeclaration of type name "int64_t"' for the typedef following the #else here, i.e. ignores the #ifdef. WTF?

#ifdef int64_t
	#warning int64_t already defined!
#else
	typedef  int64		 int64_t;
	typedef uint64		uint64_t;
#endif
*/
/*******************************************************************************
   Some useful utility macros:
*******************************************************************************/

#undef	HERE
#define	HERE	__LINE__, __FILE__

/* For my snprint-into-fixed-length-buffer calls, GCC 7.1 and later emit this kind of warning:
	"warning: ‘%s’ directive output may be truncated writing up to 1023 bytes into a region of size between [foo] and [bar] [-Wformat-truncation=]
Workarounds, including the slick typedef-based one below, are discussed here:
	https://stackoverflow.com/questions/51534284/how-to-circumvent-format-truncation-warning-in-gcc
*/
#define snprintf_nowarn(...) (snprintf(__VA_ARGS__) < 0 ? abort() : (void)0)

/* Array-bounds check, return TRUE if y-array has either endpoint in x-array's range or v.v.
(need both ways to handle e.g. one array entirely contained within the other. Consider the
following cartoon illustrating the possibilities:

                         *---------------------*
                         x1                   x2

                 *----------------------------------*
                y1                                  y2

The arrays are nonoverlapping iff  y2 < x1 or y1 > x2, so non-overlap is assured
logical converse holds, that is if (y2 >= x1) && (y1 <= x2).
*/
#define ARRAYS_DISJOINT(xarr,lenx,yarr,leny)	((yarr+leny <= xarr) || (yarr >= xarr+lenx))
#define  ARRAYS_OVERLAP(xarr,lenx,yarr,leny)	!ARRAYS_DISJOINT(xarr,lenx,yarr,leny)

/* Original version of the above:
#define ARRAYS_OVERLAP(x, lenx, y, leny)	( (x <= y) && (x+lenx) > y ) || ( (x > y) && (y+leny) > x )
*/

// 32 and 64-bit 2s-comp integer mod-add/sub macros, compute z = (x +- y)%q. Cast-enforces required uint32 typing.
// Assumes inputs properly normalized x,y < q. Allow in-place: any or all of x,y,z may refer to same operand.
#define MOD_ADD32(__x, __y, __q, __z) {\
	uint32 cy,tmp;\
	/* Need tmp for initial sum, since e.g. if x,z refer to same var the cy-check z < x always comes up 'false';
	Inputs can actually be typed signed-int32 in calling code - only the compares care about that, force uint32 there. */\
	tmp = (uint32)__x + (uint32)__y;	cy = tmp < (uint32)__x;	/* Since inputs assumed normalized (< q), cy = 1 implies q > 2^31, thus */\
	__z = tmp - (uint32)__q;	cy -= (uint32)__z > tmp;	/* tmp - q guaranteed to underflow -> no need to restore-add q in this case. */\
	__z = (uint32)__z + (cy & (uint32)__q);			/* Only possible values of cy = 0,-1 here; and need to restore-add q if cy = -1. */\
}
#define MOD_SUB32(__x, __y, __q, __z) MOD_ADD32(__x, __q - __y, __q, __z)

#define MOD_ADD64(__x, __y, __q, __z) {\
	uint64 cy,tmp;\
	/* Need tmp for initial sum, since e.g. if x,z refer to same var the cy-check z < x always comes up 'false';
	Inputs can actually be typed signed-int64 in calling code - only the compares care about that, force uint64 there. */\
	tmp = (uint64)__x + (uint64)__y;	cy = tmp < (uint64)__x;	/* Since inputs assumed normalized (< q), cy = 1 implies q > 2^63, thus */\
	__z = tmp - (uint64)__q;	cy -= (uint64)__z > tmp;	/* tmp - q guaranteed to underflow -> no need to restore-add q in this case. */\
	__z = (uint64)__z + (cy & (uint64)__q);			/* Only possible values of cy = 0,-1 here; and need to restore-add q if cy = -1. */\
}
#define MOD_SUB64(__x, __y, __q, __z) MOD_ADD64(__x, __q - __y, __q, __z)

/* Slow version of nearest-int; for doubles use faster trick,
but this is useful for reference purposes and error checking:
*/
#undef  NINT
#define NINT(x) floor(x + 0.5)

/* Fast double-float NINT. For the hand-rolled versiom, RND_A & RND_B declared in Mdata.h, defined in util.c:check_nbits_in_types.
Since the add/sub-magic-constant version depends on the compiler not optimizing things away, prefer an efficient intrinsic
whenever available:
*/
#undef  DNINT
/* Consider broadening these platform checks to "Is C99 Standard supported?" if validate no adverse performance impact */
#ifdef USE_RINT

  #if defined(VERBOSE_HEADERS) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_NVCC))
	#warning Using rint() for DNINT
  #endif
	/* E.g. CUDA kernel code needs rint(), not lrint() */
	#define DNINT(x)  rint((x))

#elif(defined(COMPILER_TYPE_ICC) || defined(COMPILER_TYPE_SUNC))

  #if defined(VERBOSE_HEADERS) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_NVCC))
	#warning Using rint() for DNINT
  #endif
	#define DNINT(x)  rint((x))

// Mustn't use lrint in 32-bit mode, since that dumps result into a long, which is only 32 bits:
#elif(defined(COMPILER_TYPE_GCC)) && (OS_BITS == 64)

  #if defined(VERBOSE_HEADERS) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_NVCC))
	#warning Using lrint() for DNINT
  #endif
	#define DNINT(x) lrint((x))

#elif(defined(COMPILER_TYPE_GCC)) && (OS_BITS == 32)

  #if defined(VERBOSE_HEADERS) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_NVCC))
	#warning Using llrint() for DNINT
  #endif
	#define DNINT(x) llrint((x))

#else
	/***NOTE:*** The util.c functions set_fpu_params() and check_nbits_in_types()
	MUST BE CALLED (in that order) AT PROGRAM INVOCATION THIS MACRO TO WORK PROPERLY!!!
	*/
  #if defined(VERBOSE_HEADERS) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_NVCC))
	#warning Using +-rnd_const DNINT emulation
  #endif
	#define DNINT(x) ((x) + RND_A) - RND_B

#endif

/* NOTE: when calling the MAX/MIN/ABS macros, do NOT allow either argument to be a
function call, since that may result in the function being called twice. */
#undef  MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#undef  MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#undef  ABS					/* 64-bit absval is not portable, so alias it. */
#define ABS(a)	((a) < 0 ? -(a) : (a))

#define IS_ODD(a)	( (int)(a) & 1)

#define IS_EVEN(a)	(~(int)(a) & 1)

// Apply sign bit b to argument x (0 = no change in sign; 1 = sign is flipped).
// No check is performed whether b is strictly 0 or 1, as required for proper behavior:
#define SGN(x,b)	((b) == 1 ? -(x) : (x))

#define BIT_SET(x,b)	( (x) |= (1 << (b)) )

#define BIT_SETC(x,b,condition)	( (x) |= ((condition) << (b)) )	// SETC = "Set conditional", bit set based on truth value of condition

#define BIT_FLIP(x,b)	( (x) ^= (1 << (b)) )

#define BIT_CLR(x,b)	( (x) &= ~(1 << (b)) )

#define BIT_TEST(x,b)	( ((x) >> (b)) & 1 )

#define	STREQ(s1,s2)	(!strcmp(s1,s2))
#define	STREQN(s1,s2,n)	(!strncmp(s1,s2,n))
#define	STRNEQ(s1,s2)	( strcmp(s1,s2))
#define	STRNEQN(s1,s2,n)( strncmp(s1,s2,n))
/* Nov 2021: Add case-insensitive versions of the above - here a hand-rolled ASCII-char-only version of STRNEQ_NOCASE
adapted (to handle null a and/or b) from https://stackoverflow.com/questions/5820810/case-insensitive-string-comparison-in-c,
in case it ever proves useful to have an own-rolled for portability:
	int STRNEQ_NOCASE(char const *a, char const *b)
	{						// If one input NULL, compare as not-equal; if both NULL compare as equal. If prefer
		int d = (a != b);	// the convention that either-NULL should compare not-equal, instead init d = (!a || !b).
		for (; a && b; a++, b++) {
			int d = tolower((unsigned char)*a) - tolower((unsigned char)*b);
			if (d != 0 || !*a)	// If all leading shared-chars are case-insensitive-equal but strlen(a) < strlen(b),
				break;			// will hit d != 0 at '\0' of a; if strlen(a) == strlen(b), will exit on hitting !*a.
		}
		return d;
	}
*/
#define	STREQ_NOCASE(s1,s2)		(!strcasecmp(s1,s2))
#define	STREQN_NOCASE(s1,s2,n)	(!strncasecmp(s1,s2,n))
#define	STRNEQ_NOCASE(s1,s2)	( strcasecmp(s1,s2))
#define	STRNEQN_NOCASE(s1,s2,n)	( strncasecmp(s1,s2,n))

/* Complex multiplication C = A*B, input as pairs of doubles.
Designed so any or all of A, B, C may point to the same memory location. */
#define CMUL(ar,ai,br,bi,cr,ci)\
{\
	double __tmp = ar;\
	ci = __tmp*bi + ai*br;\
	cr = __tmp*br - ai*bi;\
}

/* Pairs-of-doubles structs: we'd like to give these all a GCC-style
"__attribute__ ((aligned (16)))" alignment flag, but as that's not portable,
we hope other compilers will be smart enough to properly align these
without external prompting.

NOTE: GCC alas does not respect alignment specifiers for > 16-byte boundaries,
which means e.g. for the 32-byte alignment required by AVX aligned-MOV instructions
we must manually enforce such alignment using the ALLOC/ALIGN macro pairs in align.h.

Macros that manipulate complex (real, imag) pairs will have prefix CMPLX_
and will have declarations residing in the header file cmplx.h ;

Macros that manipulate non-complex pairs of doubles will have prefix VD_
and will have declarations residing in the header file vec_double.h ;

There will likely be significant overlap in the low-level functionality
(e.g. SSE2 intrinsics and data types) used by these 2 classes of macros;
which is why we union-ize the 128-bit structs storing the relevant data.
*/
#ifdef COMPILER_TYPE_GCC

	#undef float_x4
	struct float_x4
	{
		float d0;
		float d1;
		float d2;
		float d3;
	} __attribute__ ((aligned (16)));

	#undef float_x8
	struct float_x8
	{
		float d0;
		float d1;
		float d2;
		float d3;
		float d4;
		float d5;
		float d6;
		float d7;
	} __attribute__ ((aligned (16)));

	#undef float_x16
	struct float_x16
	{
		float d0;
		float d1;
		float d2;
		float d3;
		float d4;
		float d5;
		float d6;
		float d7;
		float d8;
		float d9;
		float da;
		float db;
		float dc;
		float dd;
		float de;
		float df;
	} __attribute__ ((aligned (16)));

	#undef complex
	struct complex{
		double re;
		double im;
	} __attribute__ ((aligned (16)));

	#undef double_x2
	struct double_x2
	{
		double d0;
		double d1;
	} __attribute__ ((aligned (16)));

	#undef double_x4
	struct double_x4
	{
		double d0;
		double d1;
		double d2;
		double d3;
	} __attribute__ ((aligned (16)));

	#undef double_x8
	struct double_x8
	{
		double d0;
		double d1;
		double d2;
		double d3;
		double d4;
		double d5;
		double d6;
		double d7;
	} __attribute__ ((aligned (16)));

#else

	#undef float_x4
	struct float_x4
	{
		float d0;
		float d1;
		float d2;
		float d3;
	}

	#undef float_x8
	struct float_x8
	{
		float d0;
		float d1;
		float d2;
		float d3;
		float d4;
		float d5;
		float d6;
		float d7;
	}

	#undef float_x16
	struct float_x16
	{
		float d0;
		float d1;
		float d2;
		float d3;
		float d4;
		float d5;
		float d6;
		float d7;
		float d8;
		float d9;
		float da;
		float db;
		float dc;
		float dd;
		float de;
		float df;
	}

	#undef complex
	struct complex{
		double re;
		double im;
	};

	#undef double_x2
	struct double_x2
	{
		double d0;
		double d1;
	};

	#undef double_x4
	struct double_x4
	{
		double d0;
		double d1;
		double d2;
		double d3;
	};

	#undef double_x8
	struct double_x8
	{
		double d0;
		double d1;
		double d2;
		double d3;
		double d4;
		double d5;
		double d6;
		double d7;
	};

#endif

// Basic macro used to assign same double initializer (val) to all subfields of a vec_dbl:
#ifdef USE_AVX512

	typedef struct float_x16	vec_flt;
	#define VEC_FLT_INIT(vflt_ptr, val)	( (vflt_ptr)->d0 = (vflt_ptr)->d1 = (vflt_ptr)->d2 = (vflt_ptr)->d3 = (vflt_ptr)->d4 = (vflt_ptr)->d5 = (vflt_ptr)->d6 = (vflt_ptr)->d7 = (vflt_ptr)->d8 = (vflt_ptr)->d9 = (vflt_ptr)->da = (vflt_ptr)->db = (vflt_ptr)->dc = (vflt_ptr)->dd = (vflt_ptr)->de = (vflt_ptr)->df = val )

	typedef struct double_x8	vec_dbl;
	#define VEC_DBL_INIT(vdbl_ptr, val)	( (vdbl_ptr)->d0 = (vdbl_ptr)->d1 = (vdbl_ptr)->d2 = (vdbl_ptr)->d3 = (vdbl_ptr)->d4 = (vdbl_ptr)->d5 = (vdbl_ptr)->d6 = (vdbl_ptr)->d7 = val )

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	typedef struct float_x8	vec_flt;
	#define VEC_FLT_INIT(vflt_ptr, val)	( (vflt_ptr)->d0 = (vflt_ptr)->d1 = (vflt_ptr)->d2 = (vflt_ptr)->d3 = (vflt_ptr)->d4 = (vflt_ptr)->d5 = (vflt_ptr)->d6 = (vflt_ptr)->d7 = val )

	typedef struct double_x4	vec_dbl;
	#define VEC_DBL_INIT(vdbl_ptr, val)	( (vdbl_ptr)->d0 = (vdbl_ptr)->d1 = (vdbl_ptr)->d2 = (vdbl_ptr)->d3 = val )

#elif defined(USE_SSE2)

	typedef struct float_x4	vec_flt;
	#define VEC_FLT_INIT(vflt_ptr, val)	( (vflt_ptr)->d0 = (vflt_ptr)->d1 = (vflt_ptr)->d2 = (vflt_ptr)->d3 = val )

	typedef struct double_x2	vec_dbl;
	#define VEC_DBL_INIT(vdbl_ptr, val)	( (vdbl_ptr)->d0 = (vdbl_ptr)->d1 = val )

#elif defined(__CUDACC__)

	#ifdef __CUDA_ARCH__	// This only def'd for the device-code compilation pass
		#if __CUDA_ARCH__ <= 120
			#error CUDA: This code requires double precision, please add `-arch compute_13` (or greater) to your compile flags and build on an arch of compute capability >= 1.3
		#endif
	#endif

#endif

// Also need specific-vector-double-length init macros for files containing mix of SSE2/AVX/etc SIMD code:
#ifdef USE_SSE2
	#define VEC_DBL_INIT_8(vdbl_ptr, val)	( (vdbl_ptr)->d0 = (vdbl_ptr)->d1 = (vdbl_ptr)->d2 = (vdbl_ptr)->d3 = (vdbl_ptr)->d4 = (vdbl_ptr)->d5 = (vdbl_ptr)->d6 = (vdbl_ptr)->d7 = val )
	#define VEC_DBL_INIT_4(vdbl_ptr, val)	( (vdbl_ptr)->d0 = (vdbl_ptr)->d1 = (vdbl_ptr)->d2 = (vdbl_ptr)->d3 = val )
	#define VEC_DBL_INIT_2(vdbl_ptr, val)	( (vdbl_ptr)->d0 = (vdbl_ptr)->d1 = val )
#endif

/* 128-bit vector data types for AltiVec, Cell and SSE2/3:
Alas, we can't use the above typedefs here, i.e. "vector uint32" won't work:
*/
#if(CPU_HAS_ALTIVEC || CPU_IS_CELL)
#error No support for this arch!
	typedef	vector unsigned char	vec_uint8X16 ;
	typedef	vector unsigned short	vec_uint16X8 ;
	typedef	vector unsigned int		vec_uint32X4 ;
#endif
/* Key difference is that the Cell SPU also supports vector floating doubles, which AltiVec does not.
Note that Cell also supports vector long long ints, but as there is currently no hardware arithmetic
support for same, there's little point in using them.
*/
#if(CPU_IS_CELL)
#error No support for this arch!
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
#undef uint64_32
struct uint64_32{
	uint64 d0;
	uint32 d1;
};

#undef uint96
typedef	struct uint64_32		uint96;

#undef uint64x2
struct uint64x2{
	uint64 d0;
	uint64 d1;
};

#undef uint128
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

/* 128-bit meta-int consisting of 4 uint32s: */
#undef uint32x4

	struct uint32x4{
		uint32 d0;
		uint32 d1;
		uint32 d2;
		uint32 d3;
	}
  #ifdef COMPILER_TYPE_GCC
	__attribute__ ((aligned (16)));
  #else
	;
  #endif

// *** Don't use the above to define a uint128 since we already use 64x2 struct for that ***

/* 256-bit meta-int consisting of 8 uint32s: */
#undef uint32x8

	struct uint32x8{
		uint32 d0;
		uint32 d1;
		uint32 d2;
		uint32 d3;
		uint32 d4;
		uint32 d5;
		uint32 d6;
		uint32 d7;
	}
  #ifdef COMPILER_TYPE_GCC
	__attribute__ ((aligned (16)));
  #else
	;
  #endif

// *** Don't use the above to define a uint256 since we already use 64x4 struct for that ***

/* 512-bit meta-int consisting of 16 uint32s: */
#undef uint32x16

	struct uint32x16{
		uint32 d0;
		uint32 d1;
		uint32 d2;
		uint32 d3;
		uint32 d4;
		uint32 d5;
		uint32 d6;
		uint32 d7;
		uint32 d8;
		uint32 d9;
		uint32 dA;
		uint32 dB;
		uint32 dC;
		uint32 dD;
		uint32 dE;
		uint32 dF;
	}
  #ifdef COMPILER_TYPE_GCC
	__attribute__ ((aligned (16)));
  #else
	;
  #endif

/* 256-bit int consisting of 4 uint64s: */
#undef uint64x4

	struct uint64x4{
		uint64 d0;
		uint64 d1;
		uint64 d2;
		uint64 d3;
	}
  #ifdef COMPILER_TYPE_GCC
	__attribute__ ((aligned (16)));
  #else
	;
  #endif

#undef uint256
typedef	struct uint64x4		uint256;

/* 512-bit int consisting of 8 uint64s: */
#undef uint64x8

	struct uint64x8{
		uint64 d0;
		uint64 d1;
		uint64 d2;
		uint64 d3;
		uint64 d4;
		uint64 d5;
		uint64 d6;
		uint64 d7;
	}
  #ifdef COMPILER_TYPE_GCC
	__attribute__ ((aligned (16)));
  #else
	;
  #endif

#undef uint512
typedef	struct uint64x8		uint512;

// Basic macro used to assign same uint64 initializer (val) to all subfields of a vec_u64:
#ifdef USE_AVX512

	typedef struct uint64x8	vec_u64;
	#define VEC_U64_INIT(vu64_ptr, val)	( (vu64_ptr)->d0 = (vu64_ptr)->d1 = (vu64_ptr)->d2 = (vu64_ptr)->d3 = (vu64_ptr)->d4 = (vu64_ptr)->d5 = (vu64_ptr)->d6 = (vu64_ptr)->d7 = val )

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	typedef struct uint64x4	vec_u64;
	#define VEC_U64_INIT(vu64_ptr, val)	( (vu64_ptr)->d0 = (vu64_ptr)->d1 = (vu64_ptr)->d2 = (vu64_ptr)->d3 = val )

#elif defined(USE_SSE2)

	typedef struct uint64x2	vec_u64;
	#define VEC_U64_INIT(vu64_ptr, val)	( (vu64_ptr)->d0 = (vu64_ptr)->d1 = val )

#elif defined(__CUDACC__)

	#ifdef __CUDA_ARCH__	// This only def'd for the device-code compilation pass
		#if __CUDA_ARCH__ <= 120
			#error CUDA: This code requires uint64 precision, please add `-arch compute_13` (or greater) to your compile flags and build on an arch of compute capability >= 1.3
		#endif
	#endif

#endif

// Also need specific-vector-uint64-length init macros for files containing mix of SSE2/AVX/etc SIMD code:
#ifdef USE_SSE2
	#define VEC_U64_INIT_8(vu64_ptr, val)	( (vu64_ptr)->d0 = (vu64_ptr)->d1 = (vu64_ptr)->d2 = (vu64_ptr)->d3 = (vu64_ptr)->d4 = (vu64_ptr)->d5 = (vu64_ptr)->d6 = (vu64_ptr)->d7 = val )
	#define VEC_U64_INIT_4(vu64_ptr, val)	( (vu64_ptr)->d0 = (vu64_ptr)->d1 = (vu64_ptr)->d2 = (vu64_ptr)->d3 = val )
	#define VEC_U64_INIT_2(vu64_ptr, val)	( (vu64_ptr)->d0 = (vu64_ptr)->d1 = val )
#endif

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

extern const uint512 NIL512;
extern const uint512 ONE512;
extern const uint512 TWO512;

// Case-insensitive analog of strstr - this needs the ctype.h header:
char* stristr(const char* haystack, const char* needle);

/* Binary predicates for use of stdlib qsort(): */
int ncmp_int   (const void * a, const void * b);	// Default-int compare predicate
int ncmp_uint32(const void * a, const void * b);	// Mnemonic: "Numeric CoMPare of UINT32 data"
int ncmp_sint32(const void * a, const void * b);
int ncmp_uint64(const void * a, const void * b);
int ncmp_sint64(const void * a, const void * b);

#ifdef __cplusplus
}
#endif

#endif	/* types_h_included */

