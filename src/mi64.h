/*******************************************************************************
*                                                                              *
*   (C) 1997-2016 by Ernst W. Mayer.                                           *
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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef mi64_h_included
#define mi64_h_included

#include "types.h"
#include "imul_macro.h"
#include "util.h"

#ifdef __cplusplus
extern "C" {
#endif

#undef YES_ASM
#ifdef USE_GPU	// GPU-compiles disable use of x86_64 inline-asm routines:
	#define NO_ASM
#endif
// On x86_64, compile-time NO_ASM flag allows us to override the normally auto-set YES_ASM flag:
#ifndef NO_ASM
  #if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
  #endif
#endif

// Table of precomputed low-byte (mod 2^8) inverses used to init Montgomery inversions:
static const uint8 minv8[128] = {														//  Index:
	0x01,0xAB,0xCD,0xB7,0x39,0xA3,0xC5,0xEF,0xF1,0x1B,0x3D,0xA7,0x29,0x13,0x35,0xDF,	//   0- 15
	0xE1,0x8B,0xAD,0x97,0x19,0x83,0xA5,0xCF,0xD1,0xFB,0x1D,0x87,0x09,0xF3,0x15,0xBF,	//  16- 31
	0xC1,0x6B,0x8D,0x77,0xF9,0x63,0x85,0xAF,0xB1,0xDB,0xFD,0x67,0xE9,0xD3,0xF5,0x9F,	//  32- 47
	0xA1,0x4B,0x6D,0x57,0xD9,0x43,0x65,0x8F,0x91,0xBB,0xDD,0x47,0xC9,0xB3,0xD5,0x7F,	//  48- 63
	0x81,0x2B,0x4D,0x37,0xB9,0x23,0x45,0x6F,0x71,0x9B,0xBD,0x27,0xA9,0x93,0xB5,0x5F,	//  64- 79
	0x61,0x0B,0x2D,0x17,0x99,0x03,0x25,0x4F,0x51,0x7B,0x9D,0x07,0x89,0x73,0x95,0x3F,	//  80- 95
	0x41,0xEB,0x0D,0xF7,0x79,0xE3,0x05,0x2F,0x31,0x5B,0x7D,0xE7,0x69,0x53,0x75,0x1F,	//  96-111
	0x21,0xCB,0xED,0xD7,0x59,0xC3,0xE5,0x0F,0x11,0x3B,0x5D,0xC7,0x49,0x33,0x55,0xFF		// 112-127
};

/*******************************************************************************
   Function prototypes. The corresponding function definitions will either
   be in a {function name}.c file or (for cases where a .c file contains
   multiple function definitions) in the given .c file:
*******************************************************************************/

#ifdef __CUDA_ARCH__
	#define DEV __device__
#else
	#define DEV /* */
#endif

/* RNG, Set-X-equal-to-Y, Set-X-equal-to-scalar-A, Set-X-equal-to-0, set selected bit: */
DEV void	mi64_rand				(uint64 x[], uint32 len);
DEV void	mi64_set_eq				(uint64 x[], const uint64 y[], uint32 len);
DEV void	mi64_set_eq_scalar		(uint64 x[], const uint64   a, uint32 len);
DEV void	mi64_clear				(uint64 x[], uint32 len);
DEV void	mi64_set_bit			(uint64 x[], uint32 bit);
DEV int		mi64_test_bit			(const uint64 x[], uint32 bit);

/* bitwise shift y = (x <<|>> nbits); last 2 are specialized less-than-64-bit left and rightward shifts: */
DEV uint64	mi64_shl 				(const uint64 x[], uint64 y[], uint32 nbits, uint32 len);
DEV uint64	mi64_shrl				(const uint64 x[], uint64 y[], uint32 nbits, uint32 len);
DEV uint64	mi64_shl_short_ref		(const uint64 x[], uint64 y[], uint32 nbits, uint32 len);
DEV uint64	mi64_shl_short			(const uint64 x[], uint64 y[], uint32 nbits, uint32 len);
DEV uint64	mi64_shrl_short_ref		(const uint64 x[], uint64 y[], uint32 nbits, uint32 len);
DEV uint64	mi64_shrl_short			(const uint64 x[], uint64 y[], uint32 nbits, uint32 len);

/* unsigned compare: these are all built from just two elemental compares: < and == */
DEV uint32	mi64_cmpult				(const uint64 x[], const uint64 y[], uint32 len);
DEV uint32	mi64_cmp_eq				(const uint64 x[], const uint64 y[], uint32 len);
DEV uint32	mi64_cmpult_scalar		(const uint64 x[], uint64 a, uint32 len);
DEV uint32	mi64_cmpugt_scalar		(const uint64 x[], uint64 a, uint32 len);
DEV uint32	mi64_cmp_eq_scalar		(const uint64 x[], uint64 a, uint32 len);
/* 11/14/2008: replace with macros:
uint32	mi64_cmpugt				(const uint64 x[], const uint64 y[], uint32 len);
uint32	mi64_cmpule				(const uint64 x[], const uint64 y[], uint32 len);
uint32	mi64_cmpuge				(const uint64 x[], const uint64 y[], uint32 len);
*/
#define	mi64_cmpugt(x, y, len)	 mi64_cmpult(y, x, len)
#define	mi64_cmpule(x, y, len)	!mi64_cmpugt(x, y, len)
#define mi64_cmpuge(x, y, len)	!mi64_cmpult(x, y, len)

/* POPC, leading & trailing binary zeros: */
DEV uint32	mi64_popcount			(const uint64 x[], uint32 len);
DEV int		mi64_ith_set_bit		(const uint64 x[], uint32 bit, uint32 len);
DEV uint32	mi64_trailz				(const uint64 x[], uint32 len);
DEV uint32	mi64_leadz				(const uint64 x[], uint32 len);

/* Extract leading 64 significant bits, or 128-most-significant-bits-starting-at-a-specified-bit-position. */
DEV uint32	mi64_extract_lead64		(const uint64 x[], uint32 len, uint64*result);
	double	mi64_cvt_double			(const uint64 x[], uint32 len);
DEV void 	mi64_extract_lead128	(const uint64 x[], uint32 len, uint32 nshift, uint64 lead_x[]);

/* Compare-to-zero and find-leading-nonzero-element: */
DEV uint32	mi64_iszero				(const uint64 x[], uint32 len);
DEV uint32	mi64_getlen				(const uint64 x[], uint32 len);
DEV void	mi64_setlen				(uint64 x[], uint32 oldlen, uint32 newlen);
DEV uint32	mi64_isPow2				(const uint64 x[], uint32 len, uint32*pow2);

/* add/subtract (vector +- vector) : */
DEV uint64	mi64_add_ref			(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
DEV uint64	mi64_add_cyin			(const uint64 x[], const uint64 y[], uint64 z[], uint32 len, uint64 cyin);
DEV uint64	mi64_add				(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
DEV uint64	mi64_sub_ref			(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
DEV uint64	mi64_sub				(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
DEV uint64	mi64_sub_bwin			(const uint64 x[], const uint64 y[], uint64 z[], uint32 len, uint64 bwin);

/* arithmetic and logical negation: */
DEV void	mi64_nega				(const uint64 x[], uint64 y[], uint32 len);
DEV void	mi64_negl				(const uint64 x[], uint64 y[], uint32 len);

/* add/subtract (vector +- scalar) : */
DEV uint64	mi64_add_scalar			(const uint64 x[], uint64 a, uint64 y[], uint32 len);
DEV uint64	mi64_sub_scalar			(const uint64 x[], uint64 a, uint64 y[], uint32 len);

/* unsigned multiply (vector * scalar), a*X = Y : */
DEV uint64	mi64_mul_scalar			(const uint64 x[], uint64 a, uint64 y[], uint32 len);
/* unsigned multiply (vector * scalar) and add result to vector, a*X + Y = Z : */
DEV uint64	mi64_mul_scalar_add_vec2(const uint64 x[], uint64 a, const uint64 y[], uint64 z[], uint32 len);

/* unsigned multiply (vector * vector) : */
#ifdef __CUDA_ARCH__
	/* GPU versions of these need to supply their own scratch array, in u[]: */
DEV void	mi64_mul_vector			(const uint64 x[], uint32 lenX, const uint64 y[], uint32 lenY, uint64 z[], uint32 *lenZ, uint64 u[]);
DEV void	mi64_sqr_vector			(const uint64 x[], uint64 z[], uint32 len, uint64 u[]);
DEV void	mi64_mul_vector_lo_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len, uint64 u[]);
DEV void	mi64_mul_vector_hi_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len, uint64 u[], uint64 v[]);
DEV void	mi64_mul_vector_hi_qferm(const uint64 y[], const uint64 p, const uint64 k, uint64 z[], uint32 bits, uint64 u[]);
// Specialized O(n) version of mi64_mul_vector_hi_half for moduli q = 2.k.M(p) + 1, where M(p) is a Mersenne prime
DEV void	mi64_mul_vector_hi_qmmp	(const uint64 y[], const uint64 p, const uint64 k, uint64 z[], uint32 bits, uint64 u[], uint64 v[]);
#else
DEV void	mi64_mul_vector			(const uint64 x[], uint32 lenX, const uint64 y[], uint32 lenY, uint64 z[], uint32 *lenZ);
DEV void	mi64_sqr_vector			(const uint64 x[], uint64 z[], uint32 len);
DEV void	mi64_mul_vector_lo_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
DEV void	mi64_mul_vector_hi_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
DEV void	mi64_mul_vector_hi_qferm(const uint64 y[], const uint64 p, const uint64 k, uint64 z[], uint32 bits);
// Specialized O(n) version of mi64_mul_vector_hi_half for moduli q = 2.k.M(p) + 1, where M(p) is a Mersenne prime
DEV void	mi64_mul_vector_hi_qmmp	(const uint64 y[], const uint64 p, const uint64 k, uint64 z[], uint32 bits);
#endif
//void	mi64_mul_vector_hi_fast	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
DEV void	mi64_mul_vector_hi_fast	(const uint64 y[], const uint64 p, const uint64 k, uint64 z[], uint32 len);

/* Routines needed for multiword uint64<==>double conversion and FFT-based mul: */
	uint32	mi64_cvt_uint64_double	(const uint64 x[], const uint64 y[], uint32 cy, uint32 len, double a[]);
	uint32	mi64_cvt_double_uint64	(const double a[], uint32   n, uint64 x[], uint64 y[]);

/* returns 1 if the multiword unsigned integer p is a base-z Fermat pseudoprime, 0 otherwise: */
	uint32	mi64_pprimeF			(const uint64 p[], uint64 z, uint32 len);

/* Fast division based on Montgomery multiply: */
	uint64	radix_power64(const uint64 q, const uint64 qinv, uint32 n);
DEV int		mi64_div				(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint64 r[]);	// Wrapper for the next 2 routines
	int		mi64_div_mont			(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint64 r[]);
DEV uint64	mi64_div_by_scalar64	(const uint64 x[], uint64 a, uint32 lenX, uint64 q[]);
	// x declared non-const in folded versions to permit 0-padding:
DEV uint64	mi64_div_by_scalar64_u2	(uint64 x[], uint64 a, uint32 lenX, uint64 q[]);	// 2-way interleaved-|| loop
DEV uint64	mi64_div_by_scalar64_u4	(uint64 x[], uint64 a, uint32 lenX, uint64 q[]);	// 4-way interleaved-|| loop
/* Slow bit-at-a-time division to obtain quotient q = x/y and/or remainder r = x%y: */
	int		mi64_div_binary			(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint64 r[]);

/* Fast is-divisible-by-scalar, div-by-scalar y]: */
DEV int		mi64_is_div_by_scalar32 (const uint32 x[], uint32 a, uint32 len);
DEV uint32	mi64_is_div_by_scalar32_x4(const uint32 x[], uint32 tryq0, uint32 tryq1, uint32 tryq2, uint32 tryq3, uint32 len);
DEV uint32	mi64_is_div_by_scalar32_x8(const uint32 x[], uint32 tryq0, uint32 tryq1, uint32 tryq2, uint32 tryq3, uint32 tryq4, uint32 tryq5, uint32 tryq6, uint32 tryq7, uint32 len);
DEV int		mi64_is_div_by_scalar64	(const uint64 x[], uint64 a, uint32 len);
DEV int		mi64_is_div_by_scalar64_x4(const uint64 x[], uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint32 len);
DEV int		mi64_is_div_by_scalar64_u2	(const uint64 x[], uint64 a, uint32 len);	// 2-way interleaved-|| loop
DEV int		mi64_is_div_by_scalar64_u4	(const uint64 x[], uint64 a, uint32 len);	// 4-way interleaved-|| loop
DEV uint32	mi64_div_y32			(uint64 x[], uint32 y, uint64 q[], uint32 len);

/* Basic I/O routines: */
/* This is an arbitrary-string-length core routine which can be called directly, requires caller to supply allocated-length of input string: */
	int	__convert_mi64_base10_char(char char_buf[], uint32 n_alloc_chars, const uint64 x[], uint32 len, uint32 wrap_every);
/* This is a wrapper for the above, which uses default string allocated-length = STR_MAX_LEN: */
	int	  convert_mi64_base10_char(char char_buf[], const uint64 x[], uint32 len, uint32 wrap_every);

// Leading-zero-printing analog of the above 2 functions:
	int	__convert_mi64_base10_char_print_lead0(char char_buf[], uint32 n_alloc_chars, const uint64 x[], uint32 len, uint32 ndigit, uint32 wrap_every);
	int	  convert_mi64_base10_char_print_lead0(char char_buf[], const uint64 x[], uint32 len, uint32 ndigit, uint32 wrap_every);

	uint64 *convert_base10_char_mi64(const char*char_buf, uint32 *len);

/* Modular powering: returns 1 if 2^p == 1 (mod q), and (optionally) returns the full residue in res[]. k is Mersenne-factor-specific index, used only for debugging: */
#ifdef __CUDACC__
__global__ void GPU_TF_mi64(
	const uint64*p, const uint64*pshift, uint64*gpu_thread_local, const uint32 lenP, const uint32 lenQ, const uint32 FERMAT,
	const uint32 pow2, const uint32 zshift, const uint32 start_index, const uint64 k0, uint32*kvec, const uint32 N
);
// Simple GPU-ized version of mi64_twopmodq:
__device__ uint32 mi64_twopmodq_gpu(
	const uint64*p, const uint64*pshift, uint64*gpu_thread_local, const uint32 lenP, const uint32 lenQ,
	const uint32 FERMAT, const uint32 pow2, const uint64 k, const uint32 start_index, const uint32 zshift,
	const int i	// thread id (handy for debug)
);
#endif
// These are used by factor.c, so need to be declared in both host and device-code compile pass:
	uint32	mi64_twopmodq			(const uint64 p[], uint32 len_p, const uint64 k, uint64 q[], uint32 len_q, uint64*res);
	// Specialized version of mi64_twopmodq for moduli q = 2.k.M(p) + 1, where M(p) is a Mersenne prime
	uint32	mi64_twopmodq_qmmp		(const uint64 p, const uint64 k, uint64*res);

/* Simple 32x32-bit Montgomery modular multiply - assumes inputs are either 32-bit ints or 64-bit ints < 2^32: */
#ifdef MUL_LOHI64_SUBROUTINE

	#define MONT_MUL32(__x,__y,__q,__qinv,__z)\
	{\
		uint32 lo,hi;					\
		MUL_LOHI32(__x,__y,&lo,&hi);	\
		MULL32(__qinv,lo,lo);			\
		lo = __MULH32(__q,lo);			\
		__z = hi - lo + ((-(int32)(hi < lo)) & __q);	/* did we have a borrow from (hi-lo)? */\
	}

#else

	#define MONT_MUL32(__x,__y,__q,__qinv,__z)\
	{\
		uint32 lo,hi;					\
		MUL_LOHI32(__x,__y, lo, hi);	\
		MULL32(__qinv,lo,lo);			\
		lo = __MULH32(__q,lo);			\
		__z = hi - lo + ((-(int32)(hi < lo)) & __q);	/* did we have a borrow from (hi-lo)? */\
	}

#endif

/* 48x48-bit: */
#ifdef MUL_LOHI64_SUBROUTINE

	#define MONT_MUL48(__x,__y,__q,__qinv,__z)\
	{\
		uint64 lo,hi;					\
		MUL_LOHI64(__x<<16,__y,&lo,&hi);	lo >>= 16;\
		MULL64(__qinv,lo,lo);	lo &= 0x0000FFFFFFFFFFFFull;\
		lo = __MULH64(__q<<16,lo);			\
		__z = hi - lo + ((-(int64)(hi < lo)) & __q);	/* did we have a borrow from (hi-lo)? */\
	}

#else

  #ifndef YES_ASM
	#define MONT_MUL48(__x,__y,__q,__qinv,__z)\
	{\
		uint64 lo,hi;					\
		MUL_LOHI64(__x<<16,__y, lo, hi);	lo >>= 16;\
		MULL64(__qinv,lo,lo);	lo &= 0x0000FFFFFFFFFFFFull;\
		/* Bizarre: lshifting lo instead of q in the __MULH64 arglist obviates need for above masking - and we take
		advantage of that in the ASM version below - but runs noticeably slower in this high-level C implementation: */\
		lo = __MULH64(__q<<16,lo);			\
		__z = hi - lo + ((-(int64)(hi < lo)) & __q);	/* did we have a borrow from (hi-lo)? */\
	}
  #else
	#define MONT_MUL48(__Xx,__Xy,__Xq,__Xqinv,__Xz)\
	{\
		__asm__ volatile (\
		"movq	%[__x],%%rax	\n\t"\
		"shlq	$16,%%rax		\n\t"/* __x << 16 */\
		"mulq	%[__y]			\n\t"/* MUL_LOHI64(__x,__y) -- outputs lo in rax, hi in rdx */\
		"shrq	$16,%%rax		\n\t"/* __lo >>= 16 (logical rshift) */\
		"movq	%%rdx,%%rdi		\n\t"/* move hi out of rdx in prep for next pair of MULs */\
		"imulq	%[__qinv],%%rax	\n\t"/* MULL64(__qinv,lo,lo) */\
		"shlq	$16,%%rax		\n\t"/* lo << 16 */\
		"mulq	%[__q]			\n\t"/* MULH64(__q,lo,lo) -- lo output in rdx; discard low half in rax */\
		"subq	%%rdx,%%rdi		\n\t"/* rdi = hi - lo, sets CY if there was a borrow */\
		"sbbq	%%rdx,%%rdx		\n\t"/* rdx = -CY */\
		"andq	%[__q],%%rdx	\n\t"/* ((-(int64)(hi < lo)) & __q) */\
		"addq	%%rdx,%%rdi		\n\t"/* __z in rdi */\
		"movq	%%rdi,%[__z]	\n\t"\
		:	/* outputs: none */\
		: [__x] "m" (__Xx)	/* All inputs from memory addresses here */\
		 ,[__y] "m" (__Xy)	\
		 ,[__q] "m" (__Xq)	\
		 ,[__qinv] "m" (__Xqinv)	\
		 ,[__z] "m" (__Xz)	\
		: "cc","memory","rax","rdx","rdi"	/* Clobbered registers */\
		);\
	}
  #endif


#endif

/* Various flavors of 64x64-bit Montgomery modular multiply: */
#ifdef MUL_LOHI64_SUBROUTINE

	// Our performance target here is a 64-bit CPU where 64 x 64 ==> 128-bit
	// hardware SQR no faster than general MUL, so SQR-macro here simply wraps the MUL one:
	#define MONT_SQR64(__x,__q,__qinv,__z)\
	{\
		MONT_MUL64(__x,__x,__q,__qinv,__z);\
	}

	#define MONT_MUL64(__x,__y,__q,__qinv,__z)\
	{\
		uint64 lo,hi;					\
		MUL_LOHI64(__x,__y,&lo,&hi);	\
		MULL64(__qinv,lo,lo);			\
		lo = __MULH64(__q,lo);			\
		__z = hi - lo + ((-(int64)(hi < lo)) & __q);	/* did we have a borrow from (hi-lo)? */\
	}

	// "Downshift MUL" by unity to obtain x * B^(-1) mod q:
	#define MONT_UNITY_MUL64(__x,__q,__qinv,__z)\
	{\
		uint64 lo = __x;				\
		MULL64(__qinv,lo,lo);			\
		lo = __MULH64(__q,lo);			\
		__z = -lo + ((-(int64)(lo != 0)) & __q);	/* did we have a borrow from (0-lo)? */\
	}

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

	// Our performance target here is a 64-bit CPU where 64 x 64 ==> 128-bit
	// hardware SQR no faster than general MUL, so SQR-macro here simply wraps the MUL one:
	#define MONT_SQR64(__x,__q,__qinv,__z)\
	{\
		MONT_MUL64(__x,__x,__q,__qinv,__z);\
	}

  #ifndef YES_ASM
	#define MONT_MUL64(__x,__y,__q,__qinv,__z)\
	{\
		uint64 lo,hi;					\
		MUL_LOHI64(__x,__y,lo,hi);		\
		MULL64(__qinv,lo,lo);			\
		MULH64(__q,lo,lo);			\
		__z = hi - lo + ((-(int64)(hi < lo)) & __q);	/* did we have a borrow from (hi-lo)? */\
	}
  #else
	#define MONT_MUL64(__Xx,__Xy,__Xq,__Xqinv,__Xz)\
	{\
		__asm__ volatile (\
		"movq	%[__x],%%rax	\n\t"\
		"mulq	%[__y]			\n\t"/* MUL_LOHI64(__x,__y) -- outputs lo in rax, hi in rdx */\
		"movq	%%rdx,%%rdi		\n\t"/* move hi out of rdx in prep for next pair of MULs */\
		"imulq	%[__qinv],%%rax	\n\t"/* MULL64(__qinv,lo,lo) */\
		"mulq	%[__q]			\n\t"/* MULH64(__q,lo,lo) -- lo output in rdx; discard low half in rax */\
		"subq	%%rdx,%%rdi		\n\t"/* rdi = hi - lo, sets CY if there was a borrow */\
		"sbbq	%%rdx,%%rdx		\n\t"/* rdx = -CY */\
		"andq	%[__q],%%rdx	\n\t"/* ((-(int64)(hi < lo)) & __q) */\
		"addq	%%rdx,%%rdi		\n\t"/* __z in rdi */\
		"movq	%%rdi,%[__z]	\n\t"\
		:	/* outputs: none */\
		: [__x] "m" (__Xx)	/* All inputs from memory addresses here */\
		 ,[__y] "m" (__Xy)	\
		 ,[__q] "m" (__Xq)	\
		 ,[__qinv] "m" (__Xqinv)	\
		 ,[__z] "m" (__Xz)	\
		: "cc","memory","rax","rdx","rdi"	/* Clobbered registers */\
		);\
	}
  #endif

	// "Downshift MUL" by unity to obtain x * B^(-1) mod q - on x86_64 I found no speedup from ASM-izing  this macro:
  #ifndef YES_ASM
	#define MONT_UNITY_MUL64(__x,__q,__qinv,__z)\
	{\
		uint64 lo = __x;				\
		MULL64(__qinv,lo,lo);			\
		MULH64(__q,lo,lo);				\
		__z = -lo + ((-(int64)(lo != 0)) & __q);	/* did we have a borrow from (0-lo)? */\
	}
  #else
	#define MONT_UNITY_MUL64(__Xx,__Xq,__Xqinv,__Xz)\
	{\
		__asm__ volatile (\
		"movq	%[__x],%%rax	\n\t"\
		"imulq	%[__qinv],%%rax	\n\t"/* MULL64(__qinv,lo,lo) */\
		"mulq	%[__q]			\n\t"/* MULH64(__q,lo,lo) -- lo output in rdx; discard low half in rax */\
		"negq	%%rdx			\n\t"/* rdx = hi(=0) - lo, NEG sets CY if lo 1= 0 (i.e. there would a borrow from 0-lo) */\
		"sbbq	%%rdi,%%rdi		\n\t"/* rdi = -CY */\
		"andq	%[__q],%%rdi	\n\t"/* ((-(int64)(hi < lo)) & __q) */\
		"addq	%%rdx,%%rdi		\n\t"/* __z in rdi */\
		"movq	%%rdi,%[__z]	\n\t"\
		:	/* outputs: none */\
		: [__x] "m" (__Xx)	/* All inputs from memory addresses here */\
		 ,[__q] "m" (__Xq)	\
		 ,[__qinv] "m" (__Xqinv)	\
		 ,[__z] "m" (__Xz)	\
		: "cc","memory","rax","rdx","rdi"	/* Clobbered registers */\
		);\
	}
  #endif

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

// And the multiword analogs of the above:

	/* x = input vector, lo:hi = scratch vecs, q = modulus, qinv = radix-inverse of q, z = output vec, n = #words in x/lo/hi/q/qinv */
	#define MONT_SQR_N(__x,__lo,__q,__qinv, __z,__n)\
	{\
		uint64 *__hi = __lo + __n;	/* Assume user has set aside enough storage in __lo for full double-wide product */\
		mi64_sqr_vector(__x,__lo,__n);					/* lo:hi = SQR_LOHI(x) */\
		mi64_mul_vector_lo_half(__lo,__qinv,__lo,__n);	/* lo =  MULL(lo,qinv) */\
		mi64_mul_vector_hi_half(__lo,__q   ,__lo,__n);	/* lo = UMULH(lo,q   ) */\
		if(mi64_sub(__hi,__lo,__z,__n)) {	/* did we have a borrow from (hi-lo)? */\
			mi64_add(__z,__q,__z,__n);		\
		}\
	}

//char s0[STR_MAX_LEN];
	/* Same as SQR but here have 2 distinct multiplicand-input vectors x and y: */
	#define MONT_MUL_N(__x,__y,__lo,__q,__qinv, __z,__n)\
	{\
/*printf("MONT_MUL_N: Xin = %s, ndim = %u\n", &s0[convert_mi64_base10_char(s0, __x, __n, 0)], __n);*/\
/*printf("MONT_MUL_N: Yin = %s\n", &s0[convert_mi64_base10_char(s0, __y, __n, 0)]);*/\
/*printf("MONT_MUL_N: Q    = %s\n", &s0[convert_mi64_base10_char(s0, __q, __n, 0)]);*/\
/*printf("MONT_MUL_N: QINV = %s\n", &s0[convert_mi64_base10_char(s0, __qinv, __n, 0)]);*/\
		uint64 *__hi = __lo + __n;	/* Assume user has set aside enough storage in __lo for full double-wide product */\
		uint32 pdim;				\
		mi64_mul_vector(__x,__n,__y,__n,__lo,&pdim);	/* lo:hi = MUL_LOHI(x,y) ***** result 2^64 too low! ****/\
/*printf("MONT_MUL_N: X.Y.lo = %s\n", &s0[convert_mi64_base10_char(s0, __lo, __n, 0)]);*/\
/*printf("MONT_MUL_N: X.Y.hi = %s, total result len = %u\n", &s0[convert_mi64_base10_char(s0, __hi, __n, 0)], pdim);*/\
		mi64_mul_vector_lo_half(__lo,__qinv,__lo,__n);	/* lo =  MULL(lo,qinv) */\
/*printf("MONT_MUL_N: LO *= QINV = %s\n", &s0[convert_mi64_base10_char(s0, __lo, __n, 0)]);*/\
		mi64_mul_vector_hi_half(__lo,__q   ,__lo,__n);	/* lo = UMULH(lo,q   ) */\
/*printf("MONT_MUL_N: MULH(LO,Q) = %s\n", &s0[convert_mi64_base10_char(s0, __lo, __n, 0)]);*/\
		if(mi64_sub(__hi,__lo,__z,__n)) {	/* did we have a borrow from (hi-lo)? */\
			mi64_add(__z,__q,__z,__n);		\
/*printf("MONT_MUL_N: H-L+Q = %s\n", &s0[convert_mi64_base10_char(s0, __z, __n, 0)]);*/\
		} else {\
/*printf("MONT_MUL_N: H-L = %s\n", &s0[convert_mi64_base10_char(s0, __z, __n, 0)]);*/\
		}\
	}

	/* Montgomery MUL-by-unity; x/2^(n*64) (mod q) returned in z: */
	#define MONT_UNITY_MUL_N(__x,__q,__qinv, __z,__n)\
	{\
		mi64_mul_vector_lo_half(__x,__qinv,__z,__n);	/* z =  MULL(x,qinv) */\
		mi64_mul_vector_hi_half(__z,__q   ,__z,__n);	/* z = UMULH(z,q   ) */\
		if(!mi64_iszero(__z,__n)) {		/* did we have a borrow from (hi-lo)? Since hi == 0 here, that means (lo != 0)? */\
			mi64_sub(__q,__z,__z,__n);	\
		}\
	}

#ifdef MUL_LOHI64_SUBROUTINE

	#define mi64_mul_4word(__x, __y, __out)\
	{\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5, __w6, __w7\
						,__cy2,__cy3,__cy4,__cy5,__cy6,__cy7\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l\
		,__m,__n,__o,__p,__q,__r,__s,__t,__u,__v,__w,__z;\
		\
		MUL_LOHI64(__x[0],__y[0],&__w0,&__w1);	/*   x0*y0 */\
		MUL_LOHI64(__x[1],__y[1],&__w2,&__w3);	/*   x1*y1 */\
		MUL_LOHI64(__x[2],__y[2],&__w4,&__w5);	/*   x2*y2 */\
		MUL_LOHI64(__x[3],__y[3],&__w6,&__w7);	/*   x3*y3 */\
		\
		MUL_LOHI64(__x[0],__y[1],&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x[1],__y[0],&__c ,&__d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x[0],__y[2],&__e ,&__f );	/*   x0*y2 */\
		MUL_LOHI64(__x[2],__y[0],&__g ,&__h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x[0],__y[3],&__i ,&__j );	/*   x0*y3 */\
		MUL_LOHI64(__x[1],__y[2],&__k ,&__l );	/*   x1*y2 */\
		MUL_LOHI64(__x[2],__y[1],&__m ,&__n );	/*   x2*y1 */\
		MUL_LOHI64(__x[3],__y[0],&__o ,&__p );	/*   x3*y0 */\
		\
		MUL_LOHI64(__x[1],__y[3],&__q ,&__r );	/*   x1*y3 */\
		MUL_LOHI64(__x[3],__y[1],&__s ,&__t );	/*   x3*y1 */\
		\
		MUL_LOHI64(__x[2],__y[3],&__u ,&__v );	/*   x2*y3 */\
		MUL_LOHI64(__x[3],__y[2],&__w ,&__z );	/*   x3*y2 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x0*y3 to w3-4: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x1*y2 to w3-4: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		/* Add x2*y1 to w3-4: */\
		__w3 += __m;	__cy4 += (__w3 < __m);\
		__w4 += __n;	__cy5 += (__w4 < __n);\
		\
		/* Add x3*y0 to w3-4: */\
		__w3 += __o;	__cy4 += (__w3 < __o);\
		__w4 += __p;	__cy5 += (__w4 < __p);\
		\
		/* Add x1*y3 to w4-5: */\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6  = (__w5 < __r);\
		\
		/* Add x3*y1 to w4-5: */\
		__w4 += __s;	__cy5 += (__w4 < __s);\
		__w5 += __t;	__cy6 += (__w5 < __t);\
		\
		/* Add x2*y3 to w5-6: */\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7  = (__w6 < __v);\
		\
		/* Add x3*y2 to w5-6: */\
		__w5 += __w;	__cy6 += (__w5 < __w);\
		__w6 += __z;	__cy7 += (__w6 < __z);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;	__cy6 += (__w5 < __cy5);\
		__w6 += __cy6;	__cy7 += (__w6 < __cy6);\
		__w7 += __cy7;\
		\
		/* Now split the result between __lo and __hi: */\
		__out[0] = __w0; __out[1] = __w1; __out[2] = __w2; __out[3] = __w3;\
		__out[4] = __w4; __out[5] = __w5; __out[6] = __w6; __out[7] = __w7;\
	}
	
	#define mi64_mul_lo_half_4word(__x, __y, __lo)\
	{\
		uint64 __w0, __w1, __w2, __w3\
						,__cy2,__cy3\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		MUL_LOHI64(__x[0],__y[0],&__w0,&__w1);	/*   x0*y0 */\
		MUL_LOHI64(__x[1],__y[1],&__w2,&__w3);	/*   x1*y1 */\
		\
		MUL_LOHI64(__x[0],__y[1],&__a ,&__b );	/*   x0*y1 */\
		MUL_LOHI64(__x[1],__y[0],&__c ,&__d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x[0],__y[2],&__e ,&__f );	/*   x0*y2 */\
		MUL_LOHI64(__x[2],__y[0],&__g ,&__h );	/*   x2*y0 */\
		\
		__i = __x[0]*__y[3];				/* (x0*y3).lo */\
		__j = __x[1]*__y[2];				/* (x1*y2).lo */\
		__k = __x[2]*__y[1];				/* (x2*y1).lo */\
		__l = __x[3]*__y[0];				/* (x3*y0).lo */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;\
		\
		/* Add (x0*y3 + x1*y2 + x2*y1 + x3*y0).lo to w3: */\
		__w3 += __i+__j+__k+__l;\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;\
		\
		/* Now return the result in __lo: */\
		__lo[0] = __w0; __lo[1] = __w1; __lo[2] = __w2; __lo[3] = __w3;\
	}

#else

	#define mi64_mul_4word(__x, __y, __out)\
	{\
		uint64 __w0, __w1, __w2, __w3, __w4, __w5, __w6, __w7\
						,__cy2,__cy3,__cy4,__cy5,__cy6,__cy7\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l\
		,__m,__n,__o,__p,__q,__r,__s,__t,__u,__v,__w,__z;\
		\
		MUL_LOHI64(__x[0],__y[0], __w0, __w1);	/*   x0*y0 */\
		MUL_LOHI64(__x[1],__y[1], __w2, __w3);	/*   x1*y1 */\
		MUL_LOHI64(__x[2],__y[2], __w4, __w5);	/*   x2*y2 */\
		MUL_LOHI64(__x[3],__y[3], __w6, __w7);	/*   x3*y3 */\
		\
		MUL_LOHI64(__x[0],__y[1], __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x[1],__y[0], __c , __d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x[0],__y[2], __e , __f );	/*   x0*y2 */\
		MUL_LOHI64(__x[2],__y[0], __g , __h );	/*   x2*y0 */\
		\
		MUL_LOHI64(__x[0],__y[3], __i , __j );	/*   x0*y3 */\
		MUL_LOHI64(__x[1],__y[2], __k , __l );	/*   x1*y2 */\
		MUL_LOHI64(__x[2],__y[1], __m , __n );	/*   x2*y1 */\
		MUL_LOHI64(__x[3],__y[0], __o , __p );	/*   x3*y0 */\
		\
		MUL_LOHI64(__x[1],__y[3], __q , __r );	/*   x1*y3 */\
		MUL_LOHI64(__x[3],__y[1], __s , __t );	/*   x3*y1 */\
		\
		MUL_LOHI64(__x[2],__y[3], __u , __v );	/*   x2*y3 */\
		MUL_LOHI64(__x[3],__y[2], __w , __z );	/*   x3*y2 */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;	__cy4  = (__w3 < __f);\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;	__cy4 += (__w3 < __h);\
		\
		/* Add x0*y3 to w3-4: */\
		__w3 += __i;	__cy4 += (__w3 < __i);\
		__w4 += __j;	__cy5  = (__w4 < __j);\
		\
		/* Add x1*y2 to w3-4: */\
		__w3 += __k;	__cy4 += (__w3 < __k);\
		__w4 += __l;	__cy5 += (__w4 < __l);\
		\
		/* Add x2*y1 to w3-4: */\
		__w3 += __m;	__cy4 += (__w3 < __m);\
		__w4 += __n;	__cy5 += (__w4 < __n);\
		\
		/* Add x3*y0 to w3-4: */\
		__w3 += __o;	__cy4 += (__w3 < __o);\
		__w4 += __p;	__cy5 += (__w4 < __p);\
		\
		/* Add x1*y3 to w4-5: */\
		__w4 += __q;	__cy5 += (__w4 < __q);\
		__w5 += __r;	__cy6  = (__w5 < __r);\
		\
		/* Add x3*y1 to w4-5: */\
		__w4 += __s;	__cy5 += (__w4 < __s);\
		__w5 += __t;	__cy6 += (__w5 < __t);\
		\
		/* Add x2*y3 to w5-6: */\
		__w5 += __u;	__cy6 += (__w5 < __u);\
		__w6 += __v;	__cy7  = (__w6 < __v);\
		\
		/* Add x3*y2 to w5-6: */\
		__w5 += __w;	__cy6 += (__w5 < __w);\
		__w6 += __z;	__cy7 += (__w6 < __z);\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;	__cy4 += (__w3 < __cy3);\
		__w4 += __cy4;	__cy5 += (__w4 < __cy4);\
		__w5 += __cy5;	__cy6 += (__w5 < __cy5);\
		__w6 += __cy6;	__cy7 += (__w6 < __cy6);\
		__w7 += __cy7;\
		\
		/* Now split the result between __lo and __hi: */\
		__out[0] = __w0; __out[1] = __w1; __out[2] = __w2; __out[3] = __w3;\
		__out[4] = __w4; __out[5] = __w5; __out[6] = __w6; __out[7] = __w7;\
	}
	
	#define mi64_mul_lo_half_4word(__x, __y, __lo)\
	{\
		uint64 __w0, __w1, __w2, __w3\
						,__cy2,__cy3\
		,__a,__b,__c,__d,__e,__f,__g,__h,__i,__j,__k,__l;\
		\
		MUL_LOHI64(__x[0],__y[0], __w0, __w1);	/*   x0*y0 */\
		MUL_LOHI64(__x[1],__y[1], __w2, __w3);	/*   x1*y1 */\
		\
		MUL_LOHI64(__x[0],__y[1], __a , __b );	/*   x0*y1 */\
		MUL_LOHI64(__x[1],__y[0], __c , __d );	/*   x1*y0 */\
		\
		MUL_LOHI64(__x[0],__y[2], __e , __f );	/*   x0*y2 */\
		MUL_LOHI64(__x[2],__y[0], __g , __h );	/*   x2*y0 */\
		\
		__i = __x[0]*__y[3];				/* (x0*y3).lo */\
		__j = __x[1]*__y[2];				/* (x1*y2).lo */\
		__k = __x[2]*__y[1];				/* (x2*y1).lo */\
		__l = __x[3]*__y[0];				/* (x3*y0).lo */\
		\
		/* Now add cross terms: */\
		/* Add x0*y1 to w1-2: */\
		__w1 += __a;	__cy2  = (__w1 < __a);\
		__w2 += __b;	__cy3  = (__w2 < __b);\
		\
		/* Add x1*y0 to w1-2: */\
		__w1 += __c;	__cy2 += (__w1 < __c);\
		__w2 += __d;	__cy3 += (__w2 < __d);\
		\
		/* Add x0*y2 to w2-3: */\
		__w2 += __e;	__cy3 += (__w2 < __e);\
		__w3 += __f;\
		\
		/* Add x2*y0 to w2-3: */\
		__w2 += __g;	__cy3 += (__w2 < __g);\
		__w3 += __h;\
		\
		/* Add (x0*y3 + x1*y2 + x2*y1 + x3*y0).lo to w3: */\
		__w3 += __i+__j+__k+__l;\
		\
		/* Now process carries: */\
		__w2 += __cy2;	__cy3 += (__w2 < __cy2);\
		__w3 += __cy3;\
		\
		/* Now return the result in __lo: */\
		__lo[0] = __w0; __lo[1] = __w1; __lo[2] = __w2; __lo[3] = __w3;\
	}

#endif

#undef DEV

#ifdef __cplusplus
}
#endif

#endif	/* mi64_h_included */

