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

/********************************************************************************************************/
/* Declare various Globals (externs). Unless specified otherwise, these are defined/inited in Mlucas.c: */
/********************************************************************************************************/
#ifndef INCLUDE_HWLOC
	#define INCLUDE_HWLOC 0	// v21: Make INCLUDE_HWLOC = 0 (FALSE) the default; set = 1 in compile-args to enable hwloc,
							// = 2 to further enable hwloc-related debug-printing (e.g. detailed core assignments in util.c)
#else
	#if(INCLUDE_HWLOC != 0 && INCLUDE_HWLOC != 1 && INCLUDE_HWLOC != 2)
		#error INCLUDE_HWLOC flag, if invoked, must be set = 0, 1 or 2! Aborting.
	#endif
#endif
#if INCLUDE_HWLOC
	//#warning Including HWLOC library.
	#include <hwloc.h>
	// Make sure it's at least hwloc v1:
	#if HWLOC_API_VERSION < 0x00010000
		#error Application requires installed hwloc version to be >= 1.0
	#endif

	extern hwloc_topology_t hw_topology;
	extern int HWLOC_AFFINITY;	// Is per-thread LPU-binding (affinity) supported?
#else
	//#warning NOT including HWLOC library.
#endif
// System-related globals:
extern uint32 SYSTEM_RAM, MAX_RAM_USE;	// Total usable main memory size, and max. amount of that to use per instance, in MB

// Used to force local-data-tables-reinits in cases of suspected table-data corruption:
extern int REINIT_LOCAL_DATA_TABLES;
// Normally = True; set = False on quit-signal-received to allow desired code sections to and take appropriate action:
extern int MLUCAS_KEEP_RUNNING;
typedef void sigfunc(int);
sigfunc *signal(int, sigfunc*);

// v18: Enable access to cmd-args outside main():
extern char **global_argv;

/* X87 FPU control: */
extern unsigned short FPU_64RND, FPU_64CHOP;

/* Default length of character arrays used for I/O: */
#define STR_MAX_LEN	1024

extern const char HOMEPAGE[];
/*...program version with patch suffix - value set in Mlucas.c: */
extern const char VERSION[];
extern const int CHAROFFSET;
extern int len_a;

/* These must match the smallest and largest values in the switch() in get_fft_radices(): */
#define MIN_FFT_LENGTH_IN_K			1
#define MAX_FFT_LENGTH_IN_K			524288
/* This next one should be set to log2(MAX_FFT_LENGTH_IN_K) + 10 + 4,
i.e. max. 16 bits per digit of the transform vector: */
#define MAX_PRIMALITY_TEST_BITS		63

#define MAX_SELFTEST_ITERS			1000000
extern int ITERS_BETWEEN_CHECKPOINTS;	/* number of iterations between checkpoints */
extern int DO_GCHECK;	// Mersenne/PRP or Fermat/Pepin case
extern int ITERS_BETWEEN_GCHECK_UPDATES;	// #iterations between Gerbicz-checksum updates
extern int ITERS_BETWEEN_GCHECKS;			// #iterations between Gerbicz-checksum residue-integrity checks
extern uint32 NERR_GCHECK;	// v20: Add counter for Gerbicz-check errors encountered during test

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
#else	/* 1-word P limit is set by #of bits p can have so 120*p < 2^64: */
	#define MAX_BITS_P	57
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
#define ERR_SKIP_RADIX_SET			14	// In context of self-testing, not fatal for run overall but skip the current set of FFT radices
#define ERR_INTERRUPT				15	// On one of several interrupt SIGs, exit iteration loop prematurely, write savefiles and exit
#define ERR_GERBICZ_CHECK			16
#define ERR_MAX		ERR_GERBICZ_CHECK

/***********************************************************************************************/
/* Globals. Unless specified otherwise, these are declared in Mdata.h and defined in Mlucas.c: */
/***********************************************************************************************/

extern const char *err_code[ERR_MAX];

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

// The index following 'mask' here = log2(#doubles in SIMD register) = log2(#bits in SIMD register) - 6 :
#ifdef USE_AVX512
	extern const uint32 mask03,
		br16[16],	// length-16 index-scramble array for mapping from scalar-complex to AVX512 (8 x re,8 x im)
		brinv16[16];// length-16 index-unscramble array: br[brinv[i]] = brinv[br[i]] = i .
#endif
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
extern char cbuf[STR_MAX_LEN*2], cstr[STR_MAX_LEN];
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

extern double AME,MME;			/* Avg and Max per-iteration fractional error for a given iteration interval */
extern uint32 AME_ITER_START;	/* Iteration # at which to start collecting RO Err data for AME & MME computation: */

extern uint64 PMIN;	/* minimum exponent allowed */
extern uint64 PMAX;	/* maximum exponent allowed */

// Apr 2018 - Support shift count for rotated-residue computations:
extern uint64 RES_SHIFT, GCHECK_SHIFT;
// Feb 2020: added a uint32 to keep track of the shifted-residue sign, needed for rotated residue Fermat-mod arithmetic:
extern uint32 RES_SIGN;
extern uint64 *BIGWORD_BITMAP;	/* Needed for fast how-many-residue-bits-up-to-this-array-word lookup, which the carry routines
								use to figure out where to inject the -2 in a rotated-residue LL test using Crandall/Fagin IBDWT.
								A 1-bit means a bigword (p/n+1 bits); 0 means a smallword (p/n bits), where p/n is integer-div. */
extern uint32 *BIGWORD_NBITS;	/* Array in which the (k)th element stores total number of ones bits in elts 0,...,k-1 - i.e. in all
								elts up to but not including the (k)th one - of BIGWORD_BITMAP. Init at start of LL test, then use
								to quickly compute the double-array word of a specified hift count. Cf Mlucas.c:shift_word() for use.

								For an n-word main-array, BIGWORD_BITMAP and BIGWORD_NBITS have (n/64) elts each,
								thus need 1/64 + 1/128 the total storage of the main-array. */

// For PRP tests, the base. For Pépin tests via the "FermatTest" worktype, the base defaults to 3; to do a
// Pépin test to another base, the more-general PRP-worktype must be specified with appropriate parameters.
extern uint32 PRP_BASE;
extern uint64 *BASE_MULTIPLIER_BITS;
// Nov 2020: p-1 stuff:
extern uint64 *PM1_S1_PRODUCT, PM1_S1_PROD_RES64;	// Vector to hold Stage 1 prime-powers product product in
					// most-significant-bit-deleted-and-result-bit-reversed form, and (mod 2^64) checksum on same.
extern uint32 PM1_S1_PROD_B1, PM1_S1_PROD_BITS;	// Stage 1 bound to which the current value of PM1_S1_PRODUCT corresponds, and #bits in the MSBDARBR result
extern uint32 PM1_S2_NBUF;	// # of floating-double residue-length memblocks available for Stage 2
// Allow Stage 2 bounds to be > 2^32; B2_start defaults to B1, but can be set > B1 to allow for arbitrary Stage 2 prime intervals:
extern uint32 B1;
extern uint64 B2,B2_start;
// Bit-depth of TF done on a given exponent. This is currently only used for auto-setting p-1 bounds:
extern uint32 TF_BITS;

// Work types. These must be <= 255 to be compatible with bytewise savefile format!
// Make sure to bunch all the currently-supported test types at top of list, with no numerical gaps:
extern uint32 TEST_TYPE;
#define TEST_TYPE_PRIMALITY			1
// Dec 2017: Add PRP option. Note in FermatTest mode, the default Pépin test is identical to a 3-PRP; to do a
// Pépin test to another base, the more-general PRP-worktype must be specified with appropriate parameters.
/*
July 2019: e-mail from George W. re. the different kinds of PRP residues supported by the Primenet server:

Prime95 uses a few final GMP steps to produce the same residue as if you had done a standard PRP step.
If I decipher the code below properly, if you do not do this you would be calculating residue type 3.
I've tried to convince everyone to return a standard Fermat PRP residue, gpuowl returned type 4 for quite a while.

From the code, here is the description of the 5 residue types --- quite a mess.  Ignore the "mul interim" phase as that has to do with printing interim residues.

	Below is some pseudo-code to show how various input numbers are handled for each PRP residue type (rt=residue type,
	a=PRP base, E is number for the binary exponentiation code, KF=known factors).  We can do Gerbicz error checking if b=2 and
	there are a long string of squarings -- which also greatly reduces the number of mul-by-small-consts when c<0.

	   k * 2^n + c
	if rt=1,5 E=k*2^n, gerbicz after a^k, mul interim by a^-1 if c<0, mul final by a^(c-1), compare to 1
	if rt=2   E=k*2^(n-1), gerbicz after a^k, mul interim by a^-1 if c<0, mul final by a^((c-1)/2), compare to +/-1
	if rt=3   E=k*2^n, gerbicz after a^k, mul interim by a^-1 if c<0, compare to a^-(c-1)
	if rt=4   E=k*2^(n-1), gerbicz after a^k, mul interim by a^-1 if c<0, compare to +/-a^-((c-1)/2)

	   (k * 2^n + c) / KF
	if rt=1-4 go to general case
	if rt=5   E=k*2^n, gerbicz after a^k, mul interim by a^-1 if c<0, mul final by a^(c-1), compare to a^(KF-1) mod (N/KF)

	   (k * b^n + c) / KF
	if rt=1   E=(k*b^n+c)/KF-1, compare to 1
	if rt=2   E=((k*b^n+c)/KF-1)/2, compare to +/-1
	if rt=3   E=(k*b^n+c)/KF+1, compare to a^2	 (pointless case, make rt=1)
	if rt=4   E=((k*b^n+c)/KF+1)/2, compare to +/-a	 (pointless case, make rt=2)
	if rt=5   E=k*b^n+c-1, compare mod (N/KF) to a^(KF-1)

EWM: For the v19 release, Mlucas will support only PRP type 1|5 (these 2 only differ when accompanied by an ensuing cofactor-PRP test).
*/
#define TEST_TYPE_PRP				2	// PRP-CF (Mersenne-cofactor-PRP tests) have same functionality as PRP, but with one
										// or more known factors also supplied for an ensuing Suyama-style cofactor-PRP test.
										// We differentiate the 2 test types by PRP-CF having > 0 entries in the KNOWN_FACTORS array.
#define TEST_TYPE_PM1				3
#define TEST_TYPE_MAX				3	// Set = to largest currently-supported test_type
// Remaining work types not currently supported:
#define TEST_TYPE_TF				4
#define TEST_TYPE_ECM				5
#define TEST_TYPE_DIM				(TEST_TYPE_ECM + 1)	// Set = to largest *declared* test_type value + 1

/* These must be <= 255 to be compatible with bytewise savefile format! */
extern uint32 MODULUS_TYPE;
#define MODULUS_TYPE_MERSENNE		1	// First [MODULUS_TYPE_MAX] defines here must all be software-supported
#define MODULUS_TYPE_MERSMERS		2
#define MODULUS_TYPE_FERMAT			3
#define MODULUS_TYPE_GENFFTMUL		4	// Generic-FFT-mul using [re,im]-paired 0-padded input vectors
#define MODULUS_TYPE_MAX			4	/* Set = to largest *supported* modulus type value */
// Remaining modulus types not currently supported:
#define MODULUS_TYPE_GENMERSENNE	5
#define MODULUS_TYPE_GENFERMAT		6
#define MODULUS_TYPE_PROTH			7
#define MODULUS_TYPE_EISENSTEIN		8
#define MODULUS_TYPE_DIM			(MODULUS_TYPE_EISENSTEIN + 1)	// Set = to largest *declared* modulus_type value + 1

extern uint32 TRANSFORM_TYPE;
#define REAL_WRAPPER		1
#define RIGHT_ANGLE			2
#define TRANSFORM_TYPE_MAX	2

extern const char OFILE[], WORKFILE[];
extern const char MLUCAS_INI_FILE[];
extern char CONFIGFILE[];
extern char STATFILE[];
extern char RESTARTFILE[];
extern uint64 KNOWN_FACTORS[40];	// Known prime-factors input to p-1 runs ... for now limit to 10 factors, each < 2^256

extern int INTERACT;

#undef LOG2
#define LOG2	(double)0.6931471805599453094172321215

#undef ILG2
#define ILG2	(double)1.442695040888963407359924681

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
#ifdef USE_AVX512	// AVX512 uses 512-bit registers [8 doubles]
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
#define MAX_RADIX	4096	// Largest leading-radix currently supported
extern uint32 NRADICES, RADIX_VEC[10];	/* NRADICES, RADIX_VEC[] store number & set of complex FFT radices used.	*/

extern int ROE_ITER;	// Iteration of any dangerously high ROE encountered during the current iteration interval.
						// This must be > 0, but make signed to allow sign-flip encoding of retry-fail.
extern double ROE_VAL;	// Value (must be in (0, 0.5)) of dangerously high ROE encountered during the current iteration interval
extern uint32 NERR_ROE;	// v20: Add counter for dangerously high ROEs encountered during test

extern int USE_SHORT_CY_CHAIN;
#define USE_SHORT_CY_CHAIN_MAX	3	// This is a fixed const built into the carry-step implementation via the incr_long|med|short values
									// v20: Add support for one more, hiacc, replacing the earlier compile-time HIACC flag in CY routines
extern FILE *dbg_file;
extern double*ADDR0;	// Allows for easy debug on address-read-or-write than setting a watchpoint

#ifdef MULTITHREAD
	#define MAX_CORES	1024				// Must be > 0 and a multiple of 64
	extern uint64 CORE_SET[MAX_CORES>>6];	// Bitmap for user-controlled affinity setting, as specified via the -cpu flag

	// Alas must do one-thread-at-a-time here (and then assemble the resulting data dumps)
	// to prevent overlapping file writes:
	#define	TARG_THREAD_ID	0
#else
	#define	TARG_THREAD_ID	-1
#endif

#ifdef __cplusplus
}
#endif

#endif	/* Mdata_h_included */

