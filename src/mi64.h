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
	#warning mi64.h: Defining YES_ASM.
	#define YES_ASM
  #endif
#endif

// Table of leading 2048 bits of Pi treated as a no-binary-point bitstring 1100100100001111110110... in
// to base-2^64 form, i.e. the final word 0x648... contains the leftmost 64 bits of the above expansion:
static const uint64 pi_bits[32] = {														// Input index range by row:/
	0xD42A90D5EF8E5D32ull,0xD6998B8682283D19ull,0x0AB9472D45556216ull,0x8AE9130C4C7D0288ull,
	0x1CCAA4BE754AB572ull,0xEF15E5FB4AAC0B8Cull,0xDAE2AEF837A62964ull,0xCD93C1D17603D147ull,
	0xF1CF3B960C074301ull,0x19482F23171B671Dull,0x78BA3604650C10BEull,0xB3861AA7255E4C02ull,
	0xCF6A9483B84B4B36ull,0x0E3179AB1042A95Dull,0xC1B2AE91EE51D6CBull,0x348B1FD47E9267AFull,
	0xCC6D241B0E2AE9CDull,0xE1003E5C50B1DF82ull,0x24943328F6722D9Eull,0xD74F9208BE258FF3ull,
	0xF71C35FDAD44CFD2ull,0x85FFAE5B7A035BF6ull,0x7A262174D31BF6B5ull,0xF242DABB312F3F63ull,
	0xA7F09AB6B6A8E122ull,0x98158536F92F8A1Bull,0xF7CA8CD9E69D218Dull,0x28A5043CC71A026Eull,
	0x0105DF531D89CD91ull,0x948127044533E63Aull,0x62633145C06E0E68ull,0x6487ED5110B4611Aull
};

// Table of precomputed bit-reversed bytes used for bytewise bit-reversal utilities:
static const uint8 brev8[256] = {														// Input index range by row:/
	0x00,0x80,0x40,0xC0,0x20,0xA0,0x60,0xE0,0x10,0x90,0x50,0xD0,0x30,0xB0,0x70,0xF0,//  0- f
	0x08,0x88,0x48,0xC8,0x28,0xA8,0x68,0xE8,0x18,0x98,0x58,0xD8,0x38,0xB8,0x78,0xF8,// 10-1f
	0x04,0x84,0x44,0xC4,0x24,0xA4,0x64,0xE4,0x14,0x94,0x54,0xD4,0x34,0xB4,0x74,0xF4,// 20-2f
	0x0C,0x8C,0x4C,0xCC,0x2C,0xAC,0x6C,0xEC,0x1C,0x9C,0x5C,0xDC,0x3C,0xBC,0x7C,0xFC,// 30-3f
	0x02,0x82,0x42,0xC2,0x22,0xA2,0x62,0xE2,0x12,0x92,0x52,0xD2,0x32,0xB2,0x72,0xF2,// 40-4f
	0x0A,0x8A,0x4A,0xCA,0x2A,0xAA,0x6A,0xEA,0x1A,0x9A,0x5A,0xDA,0x3A,0xBA,0x7A,0xFA,// 50-5f
	0x06,0x86,0x46,0xC6,0x26,0xA6,0x66,0xE6,0x16,0x96,0x56,0xD6,0x36,0xB6,0x76,0xF6,// 60-6f
	0x0E,0x8E,0x4E,0xCE,0x2E,0xAE,0x6E,0xEE,0x1E,0x9E,0x5E,0xDE,0x3E,0xBE,0x7E,0xFE,// 70-7f
	0x01,0x81,0x41,0xC1,0x21,0xA1,0x61,0xE1,0x11,0x91,0x51,0xD1,0x31,0xB1,0x71,0xF1,// 80-8f
	0x09,0x89,0x49,0xC9,0x29,0xA9,0x69,0xE9,0x19,0x99,0x59,0xD9,0x39,0xB9,0x79,0xF9,// 90-9f
	0x05,0x85,0x45,0xC5,0x25,0xA5,0x65,0xE5,0x15,0x95,0x55,0xD5,0x35,0xB5,0x75,0xF5,// a0-af
	0x0D,0x8D,0x4D,0xCD,0x2D,0xAD,0x6D,0xED,0x1D,0x9D,0x5D,0xDD,0x3D,0xBD,0x7D,0xFD,// b0-bf
	0x03,0x83,0x43,0xC3,0x23,0xA3,0x63,0xE3,0x13,0x93,0x53,0xD3,0x33,0xB3,0x73,0xF3,// c0-cf
	0x0B,0x8B,0x4B,0xCB,0x2B,0xAB,0x6B,0xEB,0x1B,0x9B,0x5B,0xDB,0x3B,0xBB,0x7B,0xFB,// d0-df
	0x07,0x87,0x47,0xC7,0x27,0xA7,0x67,0xE7,0x17,0x97,0x57,0xD7,0x37,0xB7,0x77,0xF7,// e0-ef
	0x0F,0x8F,0x4F,0xCF,0x2F,0xAF,0x6F,0xEF,0x1F,0x9F,0x5F,0xDF,0x3F,0xBF,0x7F,0xFF	// f0-ff
};

// Table of precomputed low-byte inverses used to init Montgomery inversion of odd input q,
// such that init qinv = minv8[(q&0xff)>>1] yields ginv satisfying q*qinv == 1 (mod 2^8):
static const uint8 minv8[128] = {														//  q-values:
	0x01,0xAB,0xCD,0xB7,0x39,0xA3,0xC5,0xEF,0xF1,0x1B,0x3D,0xA7,0x29,0x13,0x35,0xDF,	//   1: 31:2
	0xE1,0x8B,0xAD,0x97,0x19,0x83,0xA5,0xCF,0xD1,0xFB,0x1D,0x87,0x09,0xF3,0x15,0xBF,	//  33: 63:2
	0xC1,0x6B,0x8D,0x77,0xF9,0x63,0x85,0xAF,0xB1,0xDB,0xFD,0x67,0xE9,0xD3,0xF5,0x9F,	//  65: 95:2
	0xA1,0x4B,0x6D,0x57,0xD9,0x43,0x65,0x8F,0x91,0xBB,0xDD,0x47,0xC9,0xB3,0xD5,0x7F,	//  97:127:2
	0x81,0x2B,0x4D,0x37,0xB9,0x23,0x45,0x6F,0x71,0x9B,0xBD,0x27,0xA9,0x93,0xB5,0x5F,	// 129:159:2
	0x61,0x0B,0x2D,0x17,0x99,0x03,0x25,0x4F,0x51,0x7B,0x9D,0x07,0x89,0x73,0x95,0x3F,	// 161:191:2
	0x41,0xEB,0x0D,0xF7,0x79,0xE3,0x05,0x2F,0x31,0x5B,0x7D,0xE7,0x69,0x53,0x75,0x1F,	// 193:223:2
	0x21,0xCB,0xED,0xD7,0x59,0xC3,0xE5,0x0F,0x11,0x3B,0x5D,0xC7,0x49,0x33,0x55,0xFF		// 225:255:2
};
/* October 2020: Traded some e-mails/code/ideas with jeff Hurchalla re. n-bits-good initial-iterate formulae.
Via Montgomery we have, for odd q, the 5-bits-good qinv_5 = (3*q)^2, which needs just 2 ADD and 1 XOR.

On 10/7/20 3:55 AM, Jeff Hurchalla wrote:

I'll paste a formula for 6 bits I think you'll find interesting.  This comes from early this year when
I wrote a brute force program that gave me

	inv_q_6bits = ((q^12)^(q<<4))+((q|8)<<2);

The program found ~50 similar formulas (with the restriction they could use only 1 cycle operations like
add/sub/xor/and/etc), and this was the simplest.  It took half a day to complete, while the 5 bit brute
force program took under a second (with the already-known formula of (3*q)^2 as its best result).
  I kind of suspect beyond 6 bits that a formula based on simple operations would need to to store and
reuse intermediate calculations multiple times to be optimal, and that goes beyond what I tried to do,
and probably beyond what's computationally feasible anyway for brute force.

Getting back to the formula, on x86 I believe it would take 3 cycles to complete and need 5 instructions
total (we can use a single LEA instruction to combine the "<<2" with the preceding add).  I think ARM can
also combine that add/shift in 1 instruction, and also that it can combine the "<<4" with the preceding
xor (though I haven't ever tried it) - which would mean 4 instructions total, and again 3 cycles latency.

As far as I can tell it's not actually useful for anything!  Maybe a fun party trick, or good for some
48 bit hardware that I've never heard about.

I did have hopes for using one stage of a newton-like refinement with cubic convergance rather than
quadratic, so that the initial 6 good bits could become 18 good bits after only a few instructions.
A cubic-converging formula (or formulas?) definitely does exist although I don't have a link at hand,
but as I recall it just takes too many steps to be worthwhile.  The whole thing would have to beat the
5 bit good formula, which takes 2 cycles to get the initial inverse estimate, and then another 6 cycles
using Dumas' method [http://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html] to get to 20 good
bits - that's only 8 cycles total.  Since the 6 bit formula takes 3 cycles, the cubic refinement to 18
good bits would need to complete in 4 cycles or less.  This was a half year ago but I recall the cubic
method wasn't fast enough as presented - IIRC each cubic refinement took something like 6 or 7 cycles.


On Oct 8, 2020, at 1:33 PM, Ernst Mayer wrote:

That's interesting re. the 6-bits-good inverse formula - I had similarly considered writing a search
program for useful formulae beyond 5 bits but decided against it, since my sense was that the complexity
by the time one got to 8-bits would be more than simply using a LUT. But nice to see you found something
neat for 6-bits, at least.

One alternate way of going to 8 bits is via "splicing". Here a table of the high 3 bits of 8-bits-good
inverses for q = 1,3,...,255, where some fiddling shows useful patterning in a 16-column view of those data:

	0,5,6,5,1,5,6,7,7,0,1,5,1,0,1,6,
	7,4,5,4,0,4,5,6,6,7,0,4,0,7,0,5,
	6,3,4,3,7,3,4,5,5,6,7,3,7,6,7,4,
	5,2,3,2,6,2,3,4,4,5,6,2,6,5,6,3,
	4,1,2,1,5,1,2,3,3,4,5,1,5,4,5,2,
	3,0,1,0,4,0,1,2,2,3,4,0,4,3,4,1,
	2,7,0,7,3,7,0,1,1,2,3,7,3,2,3,0,
	1,6,7,6,2,6,7,0,0,1,2,6,2,1,2,7,

We see that each column contains a circular permutation of [7,6,5,4,3,2,1,0]; the left-shift counts of said
"base permutation" by column are:

	7 2 1 2 6 2 1 0 0 7 6 2 6 7 6 1

We can simply concatenate those 16 perm-values into a 48-bit int, but find it more convenient to use a
uint64, 4 bits per L-perm count, the corresponding hexits run in reverse order of the above:

	lperm = 0x1676267001262127

Then, given an odd int q whose mod-R inverse is sought:

0. Bits 1-4 and 5-7 of q correspond to the column and row indices of the above 8-row, 16-col table:

	qinv8 = (q&0xff)>>1, col = qinv8&0xf, row = qinv8>>4;

1. (3*q)^2 gives low-5-bits-good inverse, whose higher bits we mask off: qinv8 = ((q+q+q)^2)&0x1f;

2. Leftward-shift count needed for circular permutation of [7,6,5,4,3,2,1,0] is ls = (lperm >> (col<<2))&0xf;

3. High 3 bits of 8-bit inverse are hi3 = (7 - ls - row)&0x7;

4. Splice together to get 8-bits-good inverse: qinv8 += (hi3<<5);

That needs 13 1-cycle ops beyond what is needed to compute (3*q)^2.
As you noted, that much work beyond the simple 5-bits-good formula is probably an academic exercise, since
a 128-entry LUT of precomputed 8-bit inverses is going to be faster.
*/

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

/* bit-reversal, RNG, Set-X-equal-to-Y, Set-X-equal-to-scalar-A, Set-X-equal-to-0, set selected bit: */
DEV void	mi64_brev				(uint64 x[], uint32 n);
DEV void	mi64_rand				(uint64 x[], uint32 len);
DEV void	mi64_set_eq				(uint64 x[], const uint64 y[], uint32 len);
DEV void	mi64_set_eq_scalar		(uint64 x[], const uint64   a, uint32 len);
DEV void	mi64_clear				(uint64 x[], uint32 len);

/* bitwise shift y = (x <<|>> nbits); '..._short()' are specialized less-than-64-bit left and rightward shifts: */
DEV void	mi64_shlc				(const uint64 x[], uint64 y[], uint32 nbits, uint32 nshift, uint32 len, uint32 sign_flip);
DEV void	mi64_shrc				(const uint64 x[], uint64 y[], uint32 nbits, uint32 nshift, uint32 len, uint32 sign_flip);
DEV uint32	mi64_shlc_bits_align	(const uint64 x[], uint64 y[], uint32 nbits);	// This actually a hybrid of circular-shift and compare
DEV uint32	mi64_shlc_bits_limb0	(const uint64 x0, const uint64 y[], uint32 nbits);
DEV uint64	mi64_shl 				(const uint64 x[], uint64 y[], uint32 nshift, uint32 len);
DEV uint64	mi64_shrl				(const uint64 x[], uint64 y[], uint32 nshift, uint32 len, uint32 output_len);
DEV uint64	mi64_shl_short_ref		(const uint64 x[], uint64 y[], uint32 nshift, uint32 len);
DEV uint64	mi64_shl_short			(const uint64 x[], uint64 y[], uint32 nshift, uint32 len);
DEV uint64	mi64_shrl_short_ref		(const uint64 x[], uint64 y[], uint32 nshift, uint32 len);
DEV uint64	mi64_shrl_short			(const uint64 x[], uint64 y[], uint32 nshift, uint32 len);

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
DEV int		mi64_test_bit			(const uint64 x[], uint32 bit);
DEV void	mi64_set_bit			(      uint64 x[], uint32 bit, uint32 len, uint32 val);
DEV void	mi64_flip_bit			(      uint64 x[], uint32 bit, uint32 len);
DEV int		mi64_ith_set_bit		(const uint64 x[], uint32 bit, uint32 len);
DEV uint32	mi64_trailz				(const uint64 x[], uint32 len);
DEV uint32	mi64_leadz				(const uint64 x[], uint32 len);
// v20: MD5 hash. Declare x[] non-const to enable message-preprocessing, but restore x to its input value prior to return:
DEV void	mi64_md5				(uint64 x[], uint32 len, uint64 md5[], char*const md5_str);

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
DEV void	mi64_mul_vector_hi_trunc(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);	// Rhombus-truncated multiword MULH
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
	void	mi64_modpow_lr	(const uint64 a[], const uint64 b[], const uint64 n[], uint32 len, uint64 c[]);
	void	mi64_scalar_modpow_lr( uint64 a  , const uint64 b[], const uint64 n[], uint32 len, uint64 c[]);
	uint32	mi64_init_mers_or_ferm_modulus(uint64 exp, int modtype, uint64 mvec[]);
	uint32	mi64_pprimeF	(const uint64 p[], uint64 z, uint32 len);

// Hand-rolled vector double/quad conversions which work on KNL, i.e. need just AVX512F,CD:
DEV	void	mi64_vcvtuqq2pd(const uint64 a[], double b[]);
DEV	void	mi64_vcvtpd2uqq(const double b[], uint64 a[]);
// Batched 16-way[avx2] or 64-way[avx512] 53-bit vector-DP-math/FMA-based modmul:
DEV	void	mi64_modmul53_batch(const double a[], const double b[], const double m[], double r[]);
// Generic 64-bit modmul using x86 FPU for quotient:
DEV uint64	mi64_modmul64(const uint64 a, const uint64 b, const uint64 m);

/* Fast division based on Montgomery multiply: */
	uint64	radix_power64(const uint64 q, const uint64 qinv, uint32 n);
DEV int		mi64_div				(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint64 r[]);	// Wrapper for the next 2 routines
	int		mi64_div_mont			(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint64 r[]);
DEV uint64	mi64_div_by_scalar64	(const uint64 x[], uint64 a, uint32 lenX, uint64 q[]);
	// x declared non-const in folded versions to permit 0-padding:
DEV uint64	mi64_div_by_scalar64_u2	(uint64 x[], uint64 a, uint32 lenX, uint64 q[]);	// 2-way interleaved-|| loop
DEV uint64	mi64_div_by_scalar64_u4	(uint64 x[], uint64 a, uint32 lenX, uint64 q[]);	// 4-way interleaved-|| loop
/* Slow bit-at-a-time division to obtain quotient q = x/y and/or remainder r = x%y: */
	int		mi64_div_binary			(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint32*lenQ, uint64 r[]);

/* Fast is-divisible-by-scalar, div-by-scalar y]: */
DEV int		mi64_is_div_by_scalar32 (const uint32 x[], uint32 a, uint32 len);
DEV uint32	mi64_is_div_by_scalar32_x4(const uint32 x[], uint32 tryq0, uint32 tryq1, uint32 tryq2, uint32 tryq3, uint32 len);
DEV uint32	mi64_is_div_by_scalar32_x8(const uint32 x[], uint32 tryq0, uint32 tryq1, uint32 tryq2, uint32 tryq3, uint32 tryq4, uint32 tryq5, uint32 tryq6, uint32 tryq7, uint32 len);
DEV int		mi64_is_div_by_scalar64	(const uint64 x[], uint64 a, uint32 len);
DEV int		mi64_is_div_by_scalar64_x4(const uint64 x[], uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint32 len);
DEV int		mi64_is_div_by_scalar64_u2	(const uint64 x[], uint64 a, uint32 len);	// 2-way interleaved-|| loop
DEV int		mi64_is_div_by_scalar64_u4	(const uint64 x[], uint64 a, uint32 len);	// 4-way interleaved-|| loop
DEV uint32	mi64_div_y32			(const uint64 x[], uint32 y, uint64 q[], uint32 len);

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
		uint64 l_lo,l_hi;					\
		MUL_LOHI64(__x,__y,l_lo,l_hi);		\
		MULL64(__qinv,l_lo,l_lo);			\
		MULH64(__q,l_lo,l_lo);				\
		__z = l_hi - l_lo + ((-(int64)(l_hi < l_lo)) & __q);	/* did we have a borrow from (l_hi-l_lo)? */\
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

