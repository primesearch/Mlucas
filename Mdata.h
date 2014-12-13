/*******************************************************************************
*                                                                              *
*   (C) 1997-2013 by Ernst W. Mayer.                                           *
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
#ifndef Mdata_h_included
#define Mdata_h_included

#include "masterdefs.h"
#include "types.h"
#include "mi64.h"
#include "imul_macro.h"
#include "util.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
!...Module for stashing global parameters used by various Mlucas routines.
*/

/* X87 FPU control: */
extern unsigned short FPU_64RND, FPU_64CHOP;

/* Default length of character arrays used for I/O: */
#define STR_MAX_LEN	1024

/*...program version with patch suffix - value set in Mlucas.c: */
extern const char VERSION[];
extern const int CHAROFFSET;
extern int len_a;

/* These must match the smallest and largest values in the switch() in get_fft_radices(): */
#define MIN_FFT_LENGTH_IN_K			2
#define MAX_FFT_LENGTH_IN_K			245760
/* This next one should be set to log2(MAX_FFT_LENGTH_IN_K) + 10 + 4,
i.e. max. 16 bits per digit of the transform vector: */
#define MAX_PRIMALITY_TEST_BITS		63

#define MAX_SELFTEST_ITERS			1000000
extern int ITERS_BETWEEN_CHECKPOINTS;	/* number of iterations between checkpoints */

#undef	FACTOR_PASS_MAX
#define	FACTOR_PASS_MAX	16

/* Max. # of bits allowed for trial factoring is set by the largest # of
64-bit words for which we have a MUL macro in the included imul_macro.h file .

All but the "else" stuff below is specific to factor.c built in standalone mode:
*/
#if(defined(NWORD))	/* N-word P: */
	#ifndef FACTOR_STANDALONE
		#error NWORD flag requires FACTOR_STANDALONE flag, only allowed for factor.c built in standalone mode!
	#endif
	#define MAX_BITS_P	64000000
	#ifdef USE_FLOAT
		#error	USE_FLOAT not currently supported for exponents of 2 or more words!
	#else
		#define MAX_BITS_Q	64000000
	#endif
#elif(defined(P4WORD))	/* 4-word P: */
	#ifndef FACTOR_STANDALONE
		#error P4WORD flag requires FACTOR_STANDALONE flag, only allowed for factor.c built in standalone mode!
	#endif
	#ifdef USE_FLOAT
		#error	USE_FLOAT not currently supported for exponents of 4 words!
	//	#define MAX_BITS_P	128
	//	#define MAX_BITS_Q	192	// This will eventually be raised to ~210 bits
	#else
		#define MAX_BITS_P	242
		#define MAX_BITS_Q	256
	#endif
#elif(defined(P3WORD))	/* 3-word P: */
	#ifndef FACTOR_STANDALONE
		#error P3WORD flag requires FACTOR_STANDALONE flag, only allowed for factor.c built in standalone mode!
	#endif
	#ifdef USE_FLOAT
		#define MAX_BITS_P	178
		#define MAX_BITS_Q	192	// This will eventually be raised to ~210 bits
	#else
		#define MAX_BITS_P	178
		#define MAX_BITS_Q	192
	#endif
#elif(defined(P2WORD))	/* 2-word P: */
	#ifndef FACTOR_STANDALONE
		#error P2WORD flag requires FACTOR_STANDALONE flag, only allowed for factor.c built in standalone mode!
	#endif
	#define MAX_BITS_P	114
	#ifdef USE_FLOAT
		#error	USE_FLOAT not currently supported for exponents of 2 words!
	#else
		#define MAX_BITS_Q	128
	#endif
#else	/* 1-word P limit is set by #of bits p can have so 120^2*p < 2^64: */
	#define MAX_BITS_P	50
	#ifdef USE_FLOAT
		#define MAX_BITS_Q	78
	#else
		#define MAX_BITS_Q	96
	#endif
#endif

#undef  MAX_EXPO_BITS
#define MAX_EXPO_BITS	MAX_BITS_P

#undef  MAX_FACT_BITS
#define MAX_FACT_BITS	MAX_BITS_Q

// ernstMain() return-value Error codes: these must be <= 255:
// ***MAKE SURE*** to also update Mlucas.c:printMlucasErrCode() whenever make any changes to this table!
#define ERR_INCORRECT_RES64			1
#define ERR_RADIX0_UNAVAILABLE		2
#define ERR_RADIXSET_UNAVAILABLE	3
#define ERR_TESTTYPE_UNSUPPORTED	4
#define ERR_EXPONENT_ILLEGAL		5
#define ERR_FFTLENGTH_ILLEGAL		6
#define ERR_ECHECK_NOTBOOL			7
#define ERR_TESTITERS_OUTOFRANGE	8
#define ERR_ROUNDOFF				9
#define ERR_CARRY					10
#define ERR_RUN_SELFTEST_FORLENGTH	11	/* If this is the value of the lowest byte of the return value, upper 3 bytes
										assumed to contain FFT length (in K) for which we need to run a timing self-test */
#define ERR_ASSERT					12	// Assert-type fail but when we need execution to continue (e.g. self-test mode)
#define ERR_UNKNOWN_FATAL			13	// Halt execution, do not proceed to next worktodo.ini entry (e.g. data corruption suspected)
#define ERR_MAX	ERR_UNKNOWN_FATAL

/***********************************************************************************************/
/* Globals. Unless specified otherwise, these are declared in Mdata.h and defined in Mlucas.c: */
/***********************************************************************************************/

/* ESTRING stores the exponent E of the number in question, whose representation is indicated by the value
of the MODULUS_TYPE global declared below. We assume the number is of the form N = A^E + B,
with A and E positive ints and B an int of either sign). PSTRING stores some conventional HRF
representation of N, specifically for the following currently-supported modulus types:

	MODULUS_TYPE_MERSENNE:	A = 2, E is a prime, B = -1, pstring = "M{E}"	({} denote string-appending)
	MODULUS_TYPE_MERSMERS:	A = 2, E = 2^p-1 is itself a Mersenne-prime, B = -1, pstring = "M(M(p))"
	MODULUS_TYPE_FERMAT:	A = 2, E = 2^findex, B = +1, pstring = "F{findex}"
*/
extern char ESTRING[STR_MAX_LEN];	/* Exponent in string form */
extern char PSTRING[STR_MAX_LEN];	/* Number being tested in string form, typically estring concatenated with several other descriptors, e.g. strcat("M",estring) for Mersennes */

#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
	extern const uint32 mask02,
		br8[8],		// length-8 index-scramble array for mapping from scalar-complex to AVX (re,re,re,re,im,im,im,im)
		brinv8[8];	// length-8 index-unscramble array: br[brinv[i]] = brinv[br[i]] = i .
#endif
#ifdef USE_SSE2
	extern const uint32 mask01,
		br4[4];	// length-4 index-scramble array for mapping from scalar-complex (re,im,re,im) to SSE2 (re,re,im,im)
				// For length-4 this is its own inverse.
#endif

extern const int hex_chars[16];
extern char cbuf[STR_MAX_LEN];
extern char in_line[STR_MAX_LEN];
extern char *char_addr;
extern int char_offset;
extern FILE *fp, *fq;

/* Character constants and string global related to file access mode.
For simplicity's sake, we choose to only allow update mode ('+') in ascii, which allows
file access mode to be a 3-character string, with the last of these being a mandatory null terminator.
*/
extern char FILE_ACCESS_MODE[3];
#define FILE_ACCESS_READ	'r'
#define FILE_ACCESS_WRITE	'w'
#define FILE_ACCESS_APPEND	'a'
#define FILE_ACCESS_UPDATE	'+'
#define FILE_FORMAT_ASCII	'\0'
#define FILE_FORMAT_BINARY	'b'

extern uint32 N2,NRT,NRT_BITS,NRTM1;
extern uint32 SW_DIV_N;	/* Needed for the subset of radix* carry routines which support Fermat-mod. Defined in fermat_mod_square. */

/* Compare these at runtime against the analogous qfloat consts: */
#define ISRT2 (const double)0.7071067811865476	//extern const double ISRT2;
#define SQRT2 (const double)1.41421356237309504880

extern double AME,MME;	/* Avg and Max per-iteration fractional error for a given iteration interval */

/* Iteration # at which to start collecting RO Err data for AME & MME computation: */
#define AME_ITER_START (30u)

extern uint32 PMIN;	/* minimum exponent allowed */
extern uint32 PMAX;	/* maximum exponent allowed */

/* These must be <= 255 to be compatible with bytewise savefile format! */
extern uint32 TEST_TYPE;
#define TEST_TYPE_PRIMALITY			1
#define TEST_TYPE_TRIALFACTORING	2
#define TEST_TYPE_PMINUS1			3
#define TEST_TYPE_ECM				4
#define TEST_TYPE_MAX				2
#define TEST_TYPE_DIM				5	/* Set = to largest *declared* test_type value + 1 */

/* These must be <= 255 to be compatible with bytewise savefile format! */
extern uint32 MODULUS_TYPE;
#define MODULUS_TYPE_MERSENNE		1	// First [MODULUS_TYPE_MAX] defines here must all be software-supported
#define MODULUS_TYPE_MERSMERS		2
#define MODULUS_TYPE_FERMAT			3
#define MODULUS_TYPE_MAX			3	/* Set = to largest *supported* modulus type value */
#define MODULUS_TYPE_GENMERSENNE	4
#define MODULUS_TYPE_GENFERMAT		5
#define MODULUS_TYPE_PROTH			6
#define MODULUS_TYPE_EISENSTEIN		7
#define MODULUS_TYPE_DIM			8	/* Set 1 larger than largest *declared* modulus_type value */

extern uint32 TRANSFORM_TYPE;
#define REAL_WRAPPER		1
#define RIGHT_ANGLE			2
#define TRANSFORM_TYPE_MAX	2

extern const char OFILE[], RANGEFILE[];
extern const char LOCAL_INI_FILE[];
extern char CONFIGFILE[];
extern char STATFILE[];
extern char RESTARTFILE[];

extern int INTERACT;

#undef LOG2
#define LOG2	(double)0.6931471805599453

/* Allocated/defined in util.c: */
extern double RND_A, RND_B;	/* "Magic" numbers for fast floating-double NINT */

/* Various useful precomputed powers of 2 in floating-double form: */
extern double             TWO13FLINV;	/* (double)2^13 inverse */
extern double TWO25FLOAT, TWO25FLINV;	/* (double)2^25 and inverse */
extern double TWO26FLOAT, TWO26FLINV;	/* (double)2^26 and inverse */
extern double TWO32FLOAT, TWO32FLINV;	/* (double)2^32 and inverse */
extern double TWO48FLOAT, TWO48FLINV;	/* (double)2^48 and inverse */
extern double TWO50FLOAT, TWO50FLINV;	/* (double)2^50 and inverse */
extern double TWO51FLOAT, TWO51FLINV;	/* (double)2^51 and inverse */
extern double TWO52FLOAT, TWO52FLINV;	/* (double)2^52 and inverse */
extern double TWO53FLOAT, TWO53FLINV;	/* (double)2^53 and inverse */
extern double TWO54FLOAT;	/* (double)2^54 */
extern double TWO63FLOAT;	/* (double)2^63 */
extern double TWO64FLOAT, TWO64FLINV;	/* (double)2^64 and inverse */

/* Use 2^16 as fixed base for generic FFT-based mul for now: */
#define FFT_MUL_BITS	16	/* FFT_MUL_BASE = 2^(FFT_MUL_BITS) */

extern double FFT_MUL_BASE, FFT_MUL_BASE_INV;	/* Base for generic FFT-based mul */

/* Allocated/defined in FFTmul.c: */

/*...Quasi-generic cache-and-bank-conflict-avoiding array padding parameters are here. These MUST be a power of two to be
!   compatible with the simple padded-index calculation scheme we use, so we define them in terms of their base-2 logarithm.
!   Unctuous TV announcer: "Are you suffering from cache flow problems? Bank foreclosures? Then try my new array padding!"
*/
/* Make these signed so can use value < 0 as indicating uninitialized. */
extern int32 DAT_BITS, PAD_BITS;

#define DAT_BITS_DEF (10u)	/* Number of 8-byte array data in each contiguous-data block = 2^datbits.
			!* This should be chosen so a complete data block, along with a roughly equal
			!* number of FFT sincos or DWT weights data, fits in the L1 cache. 512 8-byte
			!* floats seems a good choice for an 8KB L1 cache.
			!* NOTES: - the blocklength must divide the (unpadded) vector length.
			!*        - dat_bits must be at least 5 to be compatible with radix16_wrapper_square.
			!*        - dat_bits must be at least 6 to be compatible with radix32_wrapper_square.
			!*        - dat_bits must be at least 7 to be compatible with radix64_wrapper_square.
			*/
// Number of doubles per padding insert = (1 << PAD_BITS_DEF) .
// In SIMD mode, number of padding elements between data blocks must be a multiple of SIMD vector-double count:
#ifdef USE_AVX512	// AVX512 both uses 512-bit registers [8 doubles]
	#define PAD_BITS_DEF ( 3u)
#else	// AVX uses 256-bit registers [4 doubles]; SSE2 uses 128-bit [2 doubles], but make 4 pad-doubles the minimum:
	#define PAD_BITS_DEF ( 2u)
#endif

#if !defined(PAD_BITS_DEF) || !(PAD_BITS_DEF > 0)
	#error PAD_BITS_DEF not properly defined in Mdata.h!
#endif

// In SIMD mode, number of padding elements between data blocks must be a multiple of SIMD vector-double count:
//	AVX512 both uses 512-bit registers [8 doubles]
//	AVX and AVX2 both use 256-bit registers [4 doubles]
//	SSE2 uses 128-bit [2 doubles]
#ifdef USE_SSE2
	#if ((1 << PAD_BITS_DEF) % RE_IM_STRIDE)
		#error (PAD_BITS_DEF != RE_IM_STRIDE) in Mdata.h
	#endif
#endif

/* Constant used to effect speedy NINT(x) = (x + RND) - RND.
!* We define two separate versions of the same constant to keep the compiler from
!* optimizing away the very add/subtract sequence we need to perform the rounding.
!* The general form of RND is .75*2^{number of mantissa bits in FP representation},
!* e.g. on the x86 family with its 64-bit-mantissa register float we'd use .75*2.0**64.
!* For this NINT to work, must prevent overaggressive compiler optimizations
!* which eliminate sequences like (x + y) - y, e.g. under DEC Unix
!* we invoke the < -assume accuracy_sensitive > compiler flag.
*/

#if defined(CPU_IS_X86)
/* ...x86 versions need an extra 2^11, due to 64-bit register mantissas. */

/*typedef long double reg_double;	* This defines a register-length double, which we can use for
				   local floating constants if there is no significant speed penalty for doing so. */
/* Since different compilers have different ways of identifying the x86, it's better
   to set these dynamically at runtime, in the util.c function check_nbits_in_types(). */
/*#define RND_A ((long double) 3.0*0x4000000*0x2000000*0x800)*/
/*#define RND_B ((long double)12.0*0x2000000*0x1000000*0x800)*/

#else
/* These assume IEEE64-compliant double-precision hardware arithmetic. */

/*typedef double reg_double;	* Don't want long double on non-x86, since e.g. real*16 won't fit in a floating register */
/*#define RND_A ((     double) 3.0*0x4000000*0x2000000)*/
/*#define RND_B ((     double)12.0*0x2000000*0x1000000)*/

#endif

/*...Parameters used frequently in the complex floating-point transform. */
/*
#define twopi (long double)(6.283185307179586476925286766559L)
#define ISRT2 (     double)(0.707106781186547524400844362105L)
*/
/*...Parameters used frequently in the modular transform over the ring of Gaussian integers GF(q^2)... */
/*
const uint64 q     = 2305843009213693951ull;
const uint64 qhalf = 1152921504606846976ull;
const uint64 q2    = 4611686018427387902ull;
const uint64 q4    = 9223372036854775804ull;
const long double inv61   = (long double)1.0/61;
const long double inv_m61 = (long double)1.0/2305843009213693951ull;
*/

extern int NRADICES, RADIX_VEC[10];	/* RADIX[] stores sequence of complex FFT radices used.	*/

extern int ROE_ITER;	// Iteration of any dangerously high ROE encountered during the current iteration interval. This must be > 0, but make signed to allow sign-flip encoding of retry-fail.
extern double ROE_VAL;	// Value (must be in (0, 0.5)) of dangerously high ROE encountered during the current iteration interval

extern FILE *dbg_file;
extern double*ADDR0;	// Allows for easy debug on address-read-or-write than setting a watchpoint
#ifdef MULTITHREAD
	// Alas must do one-threa-at-a-time here (and then assemble the resulting data dumps) to prevent overlapping file writes:
	#define	TARG_THREAD_ID	0
#else
	#define	TARG_THREAD_ID	-1
#endif

#ifdef __cplusplus
}
#endif

#endif	/* Mdata_h_included */

