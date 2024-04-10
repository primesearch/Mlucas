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

/******************************************************************************************
*   COMPILING AND RUNNING THE PROGRAM: see https://www.mersenneforum.org/mayer/README.html *
******************************************************************************************/
#include "Mlucas.h"
#ifndef imul_macro_h_included
	#error imul_macro.h file not included in build!
#endif

// Oct 2021: Fixed non-threadsafe bug in v20-added FFT-length reversion code, but add preprocessor
// flag to allow builders to disable it at compile time if they encounter any further bugs or performance issues:
#ifndef USE_FFTLEN_REVERSION
	#define USE_FFTLEN_REVERSION 1	// Builder would invoke -DUSE_FFTLEN_REVERSION=0 at compile time to override the default
#endif
/* Make sure none of the factoring-module-only flags are active: */
#if(defined(P4WORD) || defined(P3WORD) || defined(P2WORD))
	#error multiword exponents only allowed for factor.c built in standalone mode!
#endif
#ifdef FACTOR_STANDALONE
	#error FACTOR_STANDALONE flag only allowed for factor.c built in standalone mode!
#endif

/******************************************************************************************************/
/* Allocate storage for Globals (externs). Unless specified otherwise, these are declared in Mdata.h: */
/******************************************************************************************************/
#if INCLUDE_HWLOC
	hwloc_topology_t hw_topology;
	int HWLOC_AFFINITY = 0;	// Is per-thread LPU-binding (affinity) supported?
#endif
// System-related globals:
uint32 SYSTEM_RAM = 0;	// Total usable main memory size in MB, and max. % of that to use per instance.
#ifdef OS_TYPE_LINUX			// Linux: this is based on the value of the sysinfo "freeram" field; default is 90%.
	uint32 MAX_RAM_USE = 90;
#else							// MacOS: I've not found a reliable way to obtain free-RAM numbers, so default is 50%
	uint32 MAX_RAM_USE = 50;	// of available RAM, based on the sysctl "hw.memsize" value.
#endif

// Used to force local-data-tables-reinits in cases of suspected table-data corruption:
int REINIT_LOCAL_DATA_TABLES = 0;
// Normally = True; set = False on quit-signal-received to allow desired code sections to and take appropriate action:
int MLUCAS_KEEP_RUNNING = 1;
// v18: Enable savefile-on-interrupt-signal, access to argc/argv outside main():
char **global_argv;

FILE *dbg_file = 0x0;
double*ADDR0 = 0x0;	// Allows for easy debug on address-read-or-write than setting a watchpoint

// Define FFT-related globals (declared in Mdata.h):
uint32 N2,NRT,NRT_BITS,NRTM1;
int PFETCH_BLOCK_IDX[MAX_RADIX];// Need this for prefetch-block-index arrays
uint32 NRADICES, RADIX_VEC[10];	// NRADICES, RADIX_VEC[] store number & set of complex FFT radices used.
#ifdef MULTITHREAD
	uint64 CORE_SET[MAX_CORES>>6];	// Bitmap for user-controlled affinity setting, as specified via the -cpu flag
#endif
int ROE_ITER = 0;		// Iteration of any dangerously high ROE encountered during the current iteration interval.
uint32 NERR_ROE = 0;	// v20: Add counter for dangerously high ROEs encountered during test
						// This must be > 0, but make signed to allow sign-flip encoding of retry-fail.
double ROE_VAL = 0.0;	// Value (must be in (0, 0.5)) of dangerously high ROE encountered during the current iteration interval

int USE_SHORT_CY_CHAIN = 0;

int ITERS_BETWEEN_CHECKPOINTS;	/* number of iterations between checkpoints */
int DO_GCHECK = FALSE;	// If Mersenne/PRP or Fermat/Peoin test, Toggle to TRUE at runtime
uint32 NERR_GCHECK = 0;	// v20: Add counter for Gerbicz-check errors encountered during test
int ITERS_BETWEEN_GCHECK_UPDATES = 1000;	// iterations between Gerbicz-checkproduct updates
int ITERS_BETWEEN_GCHECKS     = 1000000;	// #iterations between Gerbicz-checksum residue-integrity checks

char ESTRING[STR_MAX_LEN];	// Mersenne exponent or Fermat-number index in string form - for M(p) this == p, for F(m) this == m
char BIN_EXP[STR_MAX_LEN];	// Binary exponent in string form - for M(p) this == p, for F(m) this == 2^m
char PSTRING[STR_MAX_LEN];	// Modulus being used in string form, e.g. "M110916929" and "F33".

// The index following 'mask' here = log2(#doubles in SIMD register) = log2(#bits in SIMD register) - 6.
// The corrsponding mask masks off one my bit than this, since we are dealing with complex FFT data which
// occupy pairs of such registers:
#ifdef USE_AVX512
	const uint32 mask03 = 0xfffffff0,
		br16[16]    = {0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15},	// length-16 index-scramble array for mapping from scalar-complex to AVX512 (8 x re,8 x im)
		brinv16[16] = {0,2,4,6,8,10,12,14,1,3,5,7,9,11,13,15};	// length-16 index-unscramble array: br[brinv[i]] = brinv[br[i]] = i .
											// EZ way to enumerate (i)th element: 'in which slot of br16[] is i?'
#endif
#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
	const uint32 mask02 = 0xfffffff8,
		br8[8]    = {0,4,1,5,2,6,3,7},	// length-8 index-scramble array for mapping from scalar-complex to AVX (re,re,re,re,im,im,im,im)
		brinv8[8] = {0,2,4,6,1,3,5,7};	// length-8 index-unscramble array: br[brinv[i]] = brinv[br[i]] = i .
#endif
#ifdef USE_SSE2
	const uint32 mask01 = 0xfffffffc,
		br4[4]    = {0,2,1,3};	// length-4 index-scramble array for mapping from scalar-complex (re,im,re,im) to SSE2 (re,re,im,im)
								// For length-4 this is its own inverse.
#endif

const int hex_chars[16] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
char cbuf[STR_MAX_LEN*2],cstr[STR_MAX_LEN];
char in_line[STR_MAX_LEN];
char *char_addr;

FILE *fp = 0x0, *fq = 0x0;	// Our convention is to always set these = 0x0 on fclose, so nullity correlates with "no open file attached"

/* File access mode is a 3-character string, with the last of these being a mandatory null terminator: */
char FILE_ACCESS_MODE[3] = {'x','y','\0'};

/* Matrix of supported moduli/test-types - init'ed in Mlucas_init */
int ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_DIM][TEST_TYPE_DIM];

/* These need to be kept updated to match the #defines in Mdata.h: */
const char *err_code[ERR_MAX] = {
	"ERR_INCORRECT_RES64",
	"ERR_RADIX0_UNAVAILABLE",
	"ERR_RADIXSET_UNAVAILABLE",
	"ERR_TESTTYPE_UNSUPPORTED",
	"ERR_EXPONENT_ILLEGAL",
	"ERR_FFTLENGTH_ILLEGAL",
	"ERR_ECHECK_NOTBOOL",
	"ERR_TESTITERS_OUTOFRANGE",
	"ERR_ROUNDOFF",
	"ERR_CARRY",
	"ERR_RUN_SELFTEST_FORLENGTH",
	"ERR_ASSERT",
	"ERR_UNKNOWN_FATAL",
	"ERR_SKIP_RADIX_SET",
	"ERR_INTERRUPT",
	"ERR_GERBICZ_CHECK"
};

// Shift count and auxiliary arrays used to support rotated-residue computations:
uint64 RES_SHIFT = 0xFFFFFFFFFFFFFFFFull;	// 0 is a valid value here, so init to UINT64_MAX, which value is treated as "uninited"
uint64 GCHECK_SHIFT = 0ull;
uint32 RES_SIGN = 0;	// Feb 2020: uint32 to keep track of shifted-residue sign flips, needed for rotated residue Fermat-mod arithmetic.
uint64 *BIGWORD_BITMAP = 0x0;
uint32 *BIGWORD_NBITS = 0x0;

// For PRP tests, the base. For Pépin tests via the "FermatTest" worktype, the base defaults to 3; to do a
// Pépin test to another base, the more-general PRP-worktype must be specified with appropriate parameters.
uint32 PRP_BASE = 0;
uint64 *BASE_MULTIPLIER_BITS = 0x0;	// Runtime-allocated bitwise multiply-by-base array
// Nov 2020: p-1 stuff:
uint64 *PM1_S1_PRODUCT = 0x0, PM1_S1_PROD_RES64 = 0ull;	// Vector to hold Stage 1 prime-powers product, and (mod 2^64) checksum on same
uint32 PM1_S1_PROD_B1 = 0, PM1_S1_PROD_BITS = 0;	// Stage 1 bound to which the current value of PM1_S1_PRODUCT corresponds, and its #bits
uint32 PM1_S2_NBUF = 0;	// # of floating-double residue-length memblocks available for Stage 2
// Allow Stage 2 bounds to be > 2^32; B2_start defaults to B1, but can be set > B1 to allow for arbitrary Stage 2 prime intervals:
uint32 B1 = 0;
uint64 B2 = 0ull, B2_start = 0ull;
// Bit-depth of TF done on a given exponent. This is currently only used for auto-setting p-1 bounds:
uint32 TF_BITS = 0;

/* These should all be set to a valid (nonzero) value at the time the appropriate test is begun */
uint32 TEST_TYPE		= 0;
uint32 MODULUS_TYPE		= 0;
uint32 TRANSFORM_TYPE	= 0;

const char HOMEPAGE  [] = "https://www.mersenneforum.org/mayer/README.html";

/* Program version with patch suffix:

For the foreseeable future, version numbers will be of form x.yz, with x a 1-digit integer,
y an int having 3 digits or less, and z an optional alphabetic patch suffix. We choose these such that
retiming is only necessary between versions that differ in x or in the leading digit of y,
i.e. we'd need to regenerate .cfg files in going from version 2.9 to 3.0 or from 3.0 to 3.141,
but not between versions 3.0x and 3.0g or between 3.1 and 3.141.

This reduces the "retime?" decision in the cfgNeedsUpdating() function to a simple comparison
of the leading 3 characters of the two version strings in question.
*/
/* Dec 2014: Implement version numbering according to the scheme:
	Major index = year - 2000
	Minor index = release # of that year, zero-indexed.
A version suffix of x, y, or z following the above numeric index indicates an [alpha,beta,gamma] (experimental,unstable) code.
A third index following release # indicates a patch number relative to that release. No 3rd index can be read as "patch number 0".
*/
const char VERSION   [] = "21.0.1";

const char OFILE     [] = "results.txt";	/* ASCII logfile containing FINAL RESULT ONLY for each
											assignment - detailed intermediate results for each assignment
											are written to the exponent-specific STATFILE (see below). */
const char WORKFILE [] = "worktodo.txt";	/* File containing exponents to be tested. New exponents
											may be appended at will, while the program is running. */

const char MLUCAS_INI_FILE[] = "mlucas.ini";	/* File containing user-customizable configuration settings [currently unused] */

char CONFIGFILE[15];						/* Configuration File: contains allowed FFT lengths
											and allows user to control (at runtime, and in
											a modifiable way) which of a predefined set
											of FFT radices gets used for each runlength.
											This file is re-read after each checkpoint.
											File is named mlucas.cfg or fermat.cfg, depending
											on whether Mersenne or Fermat number test is being done. */

/* These are set at runtime, based on the exponent being processed. */
char STATFILE   [STR_MAX_LEN];	/* ASCII logfile for the current exponent */
char RESTARTFILE[STR_MAX_LEN];	/* Restart file name(s) */
uint64 KNOWN_FACTORS[40];	// Known prime-factors input to p-1 runs ... for now limit to 10 factors, each < 2^256
int INTERACT;
double AME,MME;			/* Avg and Max per-iteration fractional error for a given iteration interval */
uint32 AME_ITER_START;	/* Iteration # at which to start collecting RO Err data for AME & MME computation: */

// These externs declared in platform.h:
int MAX_THREADS = 0;		/* Max. allowable No. of threads. */
int NTHREADS = 0;			/* actual No. of threads. If multithreading disabled, set = 1. */

uint64 PMIN;		/* minimum exponent allowed */
uint64 PMAX;		/* maximum exponent allowed depends on max. FFT length allowed
					   and will be determined at runtime, via call to given_N_get_maxP(). */

/****** END(Allocate storage for Globals (externs)). ******/

#ifndef NO_USE_SIGNALS
	void sig_handler(int signo)
	{
		if (signo == SIGINT) {
			fprintf(stderr,"received SIGINT signal.\n");	sprintf(cbuf,"received SIGINT signal.\n");
		} else if(signo == SIGTERM) {
			fprintf(stderr,"received SIGTERM signal.\n");	sprintf(cbuf,"received SIGTERM signal.\n");
	#ifndef __MINGW32__
		} else if(signo == SIGHUP) {
			fprintf(stderr,"received SIGHUP signal.\n");	sprintf(cbuf,"received SIGHUP signal.\n");
		} else if(signo == SIGALRM) {
			fprintf(stderr,"received SIGALRM signal.\n");	sprintf(cbuf,"received SIGALRM signal.\n");
		} else if(signo == SIGUSR1) {
			fprintf(stderr,"received SIGUSR1 signal.\n");	sprintf(cbuf,"received SIGUSR1 signal.\n");
		} else if(signo == SIGUSR2) {
			fprintf(stderr,"received SIGUSR2 signal.\n");	sprintf(cbuf,"received SIGUSR2 signal.\n");
	#endif
		}
	// Dec 2021: Until resolve run-to-run inconsistencies in signal handling, kill it with fire:
	exit(1);
		// Toggle a global to allow desired code sections to detect signal-received and take appropriate action:
		MLUCAS_KEEP_RUNNING = 0;
	}
#endif

/****************/

/*
!...Code to test primality of Mersenne numbers, using arbitrary-precision (array-integer) arithmetic.
!   Author: Ernst W. Mayer.
!
!   Accomplishes the Mersenne-mod squaring via the weighted discrete Fourier transform technique
!   of Crandall and Fagin (Mathematics of Computation 62 (205), pp.305-324, January 1994; also
!   available online at http://www.faginfamily.net/barry/Papers/Discrete%20Weighted%20Transforms.pdf ).
!
!   For each binary exponent p, generates first p-2 terms of the Lucas-Lehmer sequence,
!
!       s(n+1) = s(n)**2 - 2, s(0) = 4
!
!   modulo N = 2^p - 1. If N divides s(p-2) (specifically, the (p-2)th residue == 0 modulo N), N is prime.
!
!   See https://www.mersenneforum.org/mayer/README.html for build instructions, recent revision history
!   of the code and a summary of what has been changed/added/deleted in the latest release.
!
!***TO DO LIST:
!
!   Performance-related:
!   (0) Continue to improve cache/TLB behavior of code. This should include the following enhancements at a minimum:
!
!    - Extend block-FFT strategy to include more than just the initial FFT pass, thus allowing the blocks
!      to be smaller than n/radix(pass 1). This should help at very large runlengths and on systems with
!      small L1 caches.
!
!   (1) Continue to optimize SIMD assembler; add support for next-gen AVX512 vector instructions.
!
!   (2) Implement hybrid complex floating-point FFT with a modular complex (i.e. Gaussian integer) transform
!       over a suitable Galois field GF(Mp^2), where Mp = 2^31 - 1. Based on vector integer-math capabilities of
!       current bleeding-edge x86 processors, p = 31 is the prime candidate (pun intended).
!       For background on the mathematics, see Richard Crandall's preprint, "Integer convolution via split-radix fast
!       Galois transform", freely available at: http://academic.reed.edu/physics/faculty/crandall/papers/confgt.pdf .
!
!   Functionality-related:
!   (*) P-1 factoring module (together with hand-rolled subquadratic, memory-efficient GCD)
!   (*) An (optional) GUI would be nice...
!
!
!***ACKNOWLEDGEMENTS: special thanks are due to the following people:
!
!  * Richard Crandall and Barry Fagin - for the discrete weighted transform.
!
!  * Peter Montgomery - for help with theoretical and implementation questions,
!    loan of his fast factoring code, and for testing early versions of the program on his SGI
!    (and uncovering several compiler bugs along the way).
!
!  * Luke Welsh - for historical perspective, und besonders fu"r das Max und Moritz Weissbier Glas.
!
!  * George Woltman - for organizing the Great Internet Mersenne Prime Search,
!    and for freely sharing his time and keen computational insights.
!
!  * Scott Kurowski - for the GIMPS PrimeNet server.
!
!  * Jason Papadopoulos - for his valuable perspectives regarding FFT machine
!    implementations and algorithmic optimization for various machine architectures.
!
!  * Guillermo Ballester Valor - for many interesting discussions on optimization
!       and useful suggestions for improving the code.
!
!  * John Pierce - for hosting the ftp archive on his hogranch.com server for many years.
!
!  * David Stanfill - for generous access to high-end Xeon server machinery and for hosting the "GIMPS KNL".
!
!    ...as well as the numerous people who have been kind enough to try the code out
!    and send timings, compiler options for various architectures and suggestions
!    for improvement.
!
!***For UPDATES, PATCHES and USAGE INSTRUCTIONS, see
!
!    https://www.mersenneforum.org/mayer/README.html
!
!***Send QUESTIONS, COMMENTS or SUGGESTIONS to me at <ewmayer@aol.com>.
!
!    Alright, let's go hunt for some primes!
!
*/

// Mod-2^64 sum of elements of uint64 array a[] ... in the context of the G-checkproduct integrity
// checking we send it a double-float array treated as a uint64 one via type-punning cast of a[]:
uint64 sum64(uint64 a[], uint32 n) {
	uint32 i;
	uint64 sum = 0ull;
	for(i = 0; i < n; i++)
		sum += a[i];
	return sum;
}

// Simple majority-vote consensus:
uint64 consensus_checksum(uint64 s1, uint64 s2, uint64 s3) {
	if(s1 == s2) return s1;
	if(s1 == s3) return s1;
	if(s2 == s3) return s2;
	return 0;
}

/* The return value here is considered as consisting of bytewise subfields. The lowest byte contains
one of the error codes declared in Mdata.h (or 0 if successful exit); the upper bytes should only be nonzero
for certain specific value of the lowest-order byte, as indicated in Mdata.h:
*/
uint32	ernstMain
(
	int		mod_type,
	int		test_type,
	uint64	exponent,
	uint32	fft_length,
	int		radix_set,
	uint32	maxFFT,
	uint32	iterations,	/* Use to store log2[max factor depth] in TF mode; ignore in p-1 mode */
	uint64	*sh0,	/* Reference/Computed mod-(2^64) residue */
	uint64	*sh1,	/* Reference/Computed Selfridge/Hurwitz Residue mod 2^35-1 */
	uint64	*sh2,	/* Reference/Computed Selfridge/Hurwitz Residue mod 2^36-1 */
	int		scrnFlag,
	double	*runtime
)
{
	int resFlag = 0;
/*...Various general run parameters.. */
	uint32 timing_test_iters = 0;
/*...scalars and fixed-size arrays...	*/
	uint32 i,j,k = 0;
	/* TODO: some of these need to become 64-bit: */
	uint32 dum = 0,findex = 0,ierr = 0,ilo = 0,ihi = 0,iseed,isprime,kblocks = 0,maxiter = 0,n = 0,npad = 0;
	uint64 itmp64,cy, s1 = 0ull,s2 = 0ull,s3 = 0ull;	// s1,2,3: Triply-redundant whole-array checksum on b,c-arrays used in the G-check
	uint32 mode_flag = 0, first_sub, last_sub;
	/* Exponent of number to be tested - note that for trial-factoring, we represent p
	strictly in string[STR_MAX_LEN] form in this module, only converting it to numeric
	form in the factoring module. For all other types of assignments uint64 should suffice: */
	uint64 p = 0, i1,i2,i3, rmodb,mmodb;
	uint32 nbits_in_p = 0, nfac;
	/* Res64 and Selfridge-Hurwitz residues: */
	uint64 Res64, Res35m1, Res36m1;
/*...Known Mersenne prime exponents. This array must be null-terminated.	*/
	// Dec 2018: Including M51, there are (31, 19) p = 1,3 (mod 4), resp., vs (25.8, 24.5) predicted
	// (for p < 10^8) by the Lenstra/Wagstaff heuristic (cf. est_num_mp_in_interval() in util.c):
	const uint32 knowns[] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701	// M#1-25
		,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011		// M#26-40
		,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,0x0};	// M#41-51

/*...What a bunch of characters...	*/
	char *cptr = 0x0, *endp, gcd_str[STR_MAX_LEN], aid[33] = "\0";	// 32-hexit Primenet assignment id needs 33rd char for \0
/*...initialize logicals and factoring parameters...	*/
	int restart = FALSE, use_lowmem = 0, check_interval = 0;

#if INCLUDE_TF
	uint32 bit_depth_todo = 0;
	uint64 factor_k_start = 0;
	uint32 factor_pass_start = 0, factor_pass_hi = 0;
	double log2_min_factor = 0, log2_max_factor = 0;
#endif
	double tests_saved = 0.0;	// v21: make this a dfloat to allow fractional-parts
	uint32 pm1_done = FALSE, split_curr_assignment = FALSE, s2_continuation = FALSE, s2_partial = FALSE;
	uint32 pm1_bigstep = 0, pm1_stage2_mem_multiple = 0, psmall = 0;
/*...allocatable data arrays and associated params: */
	static uint64 nbytes = 0, nalloc = 0, arrtmp_alloc = 0, s1p_alloc = 0;
	static double *a_ptmp = 0x0, *a = 0x0, *b = 0x0, *c = 0x0, *d = 0x0, *e = 0x0;
	// uint64 scratch array and 4 pointers used to store cast-to-(uint64*) version of above b,c,d,e-pointers
	static uint64 *arrtmp = 0x0, *b_uint64_ptr = 0x0, *c_uint64_ptr = 0x0, *d_uint64_ptr = 0x0, *e_uint64_ptr = 0x0;
	double final_res_offset;

/*...time-related stuff. clock_t is typically an int (signed 32-bit)
	and a typical value of CLOCKS_PER_SEC (e.g. as defined in <machine/machtime.h>
	on Alpha TruUnix) is 1000000, so we should accumulate clock_t-stored time
	differences at least roughly every half hour (preferably every second or so)
	to avoid weirdness due to flipping of the sign bit or integer overflow.
*/
	/*clock_t clock1, clock2;	Moved these to [mers|fermat]_mod_square.c */
	double tdiff,tdif2;
  #define SIZE 256
	time_t calendar_time;
	struct tm *local_time, *gm_time;
	char timebuffer[SIZE];

/*...entry point for one or more Lucas-Lehmer tests is here.	*/
	MODULUS_TYPE = mod_type;
	TEST_TYPE = test_type;
	INTERACT = FALSE;

RANGE_BEG:

	p = 0ull; ierr = 0;
	USE_SHORT_CY_CHAIN = 0;		// v19: Reset carry-chain length fiddler to default (faster/lower-accuracy) at start of each run:
	ROE_ITER = 0; ROE_VAL = 0.0;
	NERR_GCHECK = NERR_ROE = 0;	// v20: Add counters for Gerbicz-check errors and dangerously high ROEs encountered
								// during test - if a restart, will re-read actual cumulative values from checkpoint file.
	// Clear out any FFT-radix or known-factor data that might remain from a just-completed run:
	for(i = 0; i < 10; i++) { RADIX_VEC[i] = 0; }
	nfac = 0; mi64_clear(KNOWN_FACTORS,40);
	NRADICES = 0;
	RESTARTFILE[0] = STATFILE[0] = '\0';
	restart = FALSE;
	B1 = 0; B2 = B2_start = 0ull; gcd_str[0] = '\0'; split_curr_assignment = s2_continuation = s2_partial = FALSE;
	pm1_bigstep = pm1_stage2_mem_multiple = psmall = 0;

	// Check for user-set value of various flags. Failure-to-find-or-parse results in isNaN(dtmp) = TRUE, print nothing in that case:
	double dtmp = mlucas_getOptVal(MLUCAS_INI_FILE,"LowMem");
	if(dtmp != 0) {
		if(dtmp != dtmp) {	// isNaN is C99, want something that also works on pre-C99 platforms
			sprintf(cbuf,"User did not set LowMem in %s ... allowing all test types.\n",MLUCAS_INI_FILE);
		} else if(dtmp == 1) {
			sprintf(cbuf,"User set LowMem = 1 in %s ... this allows PRP-testing but excludes p-1 stage 2.\n",MLUCAS_INI_FILE);	use_lowmem = 1;
		} else if(dtmp == 2) {
			sprintf(cbuf,"User set LowMem = 2 in %s ... this excludes both PRP-testing and p-1 stage 2.\n",MLUCAS_INI_FILE);	use_lowmem = 2;
		} else {
			sprintf(cbuf,"User set unsupported value LowMem = %f in %s ... ignoring.\n",dtmp,MLUCAS_INI_FILE);
		}
		mlucas_fprint(cbuf,1);
	}
	dtmp = mlucas_getOptVal(MLUCAS_INI_FILE,"CheckInterval");
	if(dtmp != 0) {
		if(dtmp != dtmp) {
			sprintf(cbuf,"User did not set CheckInterval in %s ... using default.\n",MLUCAS_INI_FILE);
		} else if(dtmp < 1000 || dtmp > 1000000) {
			sprintf(cbuf,"User set CheckInterval = %f in %s ... values < 10^3 or > 10^6 are not supported, ignoring.\n",dtmp,MLUCAS_INI_FILE);
		} else if(DNINT(dtmp) != dtmp) {
			sprintf(cbuf,"User set non-whole-number CheckInterval = %f in %s ... ignoring.\n",dtmp,MLUCAS_INI_FILE);
		} else {
			sprintf(cbuf,"User set CheckInterval = %d in %s.\n",(int)dtmp,MLUCAS_INI_FILE);	check_interval = (int)dtmp;
		}
		mlucas_fprint(cbuf,1);
	}

/*  ...If multithreading enabled, set max. # of threads based on # of available (logical) processors,
with the default #threads = 1 and affinity set to logical core 0, unless user overrides those via -nthread or -cpu:
*/
#ifdef MULTITHREAD

  #ifdef USE_OMP
	// OpenMP not currently supported (attempting to build with this #define enabled barfs in
	// preprocessing via #error in platform.h), this is merely placeholder for possible future use:
	ASSERT(MAX_THREADS = omp_get_num_procs(), "Illegal #Cores value stored in MAX_THREADS");
  #elif(defined(USE_PTHREAD))
	ASSERT(MAX_THREADS =     get_num_cores(), "Illegal #Cores value stored in MAX_THREADS");
  #else
	#error Unrecognized multithreading model!
  #endif
	// MAX_THREADS based on number of processing cores will most often be a power of 2, but don't assume that.
	ASSERT(MAX_THREADS > 0,"MAX_THREADS must be > 0");
	ASSERT(MAX_THREADS <= MAX_CORES,"MAX_THREADS exceeds the MAX_CORES setting in Mdata.h .");

	if(!NTHREADS) {
		NTHREADS = 1;
		fprintf(stderr,"No CPU set or threadcount specified ... running single-threaded.\n");
		// Use the same affinity-setting code here as for the -cpu option, but simply for cores [0:NTHREADS-1]:
		sprintf(cbuf,"0:%d",NTHREADS-1);
		parseAffinityString(cbuf);
	} else if(NTHREADS > MAX_CORES) {
		sprintf(cbuf,"ERROR: NTHREADS = %d exceeds the MAX_CORES setting in Mdata.h = %d\n", NTHREADS, MAX_CORES);
		ASSERT(0, cbuf);
	} else {	// In timing-test mode, allow #threads > #cores
		if(NTHREADS > MAX_THREADS) {
			fprintf(stderr,"WARN: NTHREADS = %d exceeds number of cores = %d\n", NTHREADS, MAX_THREADS);
		}
		fprintf(stderr,"NTHREADS = %d\n", NTHREADS);
	}

#else

		MAX_THREADS = NTHREADS = 1;

#endif	// #ifdef MULTITHREAD ?

	/* Make number of iterations between checkpoints dependent on #threads -
	don't want excessively frequent savefile writes, at most 1 or 2 an hour is needed:
	*/
// Oct 2021: Allow this to be set via CheckInterval option in mlucas.ini and read at program-start time,
// albeit for PRP-tests subject to constraints related to Gerbicz checking:
	if(!check_interval) {
		if(NTHREADS > 4)
			ITERS_BETWEEN_CHECKPOINTS = 100000;
		else
			ITERS_BETWEEN_CHECKPOINTS =  10000;
	} else if(check_interval < 1000) {
		ASSERT(0,"User-set value of check_interval must >= 1000.");
	} else
		ITERS_BETWEEN_CHECKPOINTS = check_interval;

	fprintf(stderr,"Setting ITERS_BETWEEN_CHECKPOINTS = %u.\n",ITERS_BETWEEN_CHECKPOINTS);

	i = ITERS_BETWEEN_GCHECKS;
	j = ITERS_BETWEEN_GCHECK_UPDATES;
	ASSERT(i == j*j, "#iterations between Gerbicz-checksum updates must = sqrt(#iterations between residue-integrity checks)");
	// v19: If PRP test, make sure Gerbicz-checkproduct interval divides checkpoint-writing one.
	// If not true, merely warn here because user may be doing LL/DC/p-1 and not PRP-tests:
	k = ITERS_BETWEEN_CHECKPOINTS;
	if(i%k != 0 || k%j != 0)
		fprintf(stderr,"WARN: G-checkproduct update interval must divide savefile-update one, which must divide the G-check interval ... this will trigger an assertion-exit if a PRP-test is attempted.");

	// Alloc bitwise multiply-by-base array, needed to support P-1 factoring and PRP testing:
	if(!BASE_MULTIPLIER_BITS) {
		j = ((ITERS_BETWEEN_CHECKPOINTS+63) >> 6) + 1;	// Add 1 pad element in case compiler does not 64-bit align
		BASE_MULTIPLIER_BITS = ALLOC_UINT64(BASE_MULTIPLIER_BITS, j);	if(!BASE_MULTIPLIER_BITS){ sprintf(cbuf, "ERROR: unable to allocate BASE_MULTIPLIER_BITS array in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		BASE_MULTIPLIER_BITS = ALIGN_UINT64(BASE_MULTIPLIER_BITS);	ASSERT(((intptr_t)BASE_MULTIPLIER_BITS & 63) == 0x0,"BASE_MULTIPLIER_BITS[] not aligned on 64-byte boundary!");
		for(i = 0; i < j; i++) { BASE_MULTIPLIER_BITS[i] = 0ull; }	// v20: Init = 0 here, in case we jump directly into p-1 stage 2 on restart
	}

	/* Look for work file...	*/
 fp = 0x0;
	if (!exponent || (exponent!=0 && fft_length!=0 && iterations==0)) {
		fprintf(stderr," looking for %s file...\n",WORKFILE);
		fp = mlucas_fopen(WORKFILE, "r");
	}

	/***********************************************************************************/
	/* Automated exponent-dispatch mode - Note that in the assignment-syntax comments, */
	/* <> and [] mean required and optional arguments, respectively:                   */
	/***********************************************************************************/
	if (fp) {
		fprintf(stderr," %s file found...reading next assignment...\n",WORKFILE);

	  read_next_assignment:	// Read first line of worktodo.ini file into 1K character array:
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			fprintf(stderr,"Hit EOF while attempting to read next line of worktodo ... quitting.\n");
			exit(0);
		}
		fprintf(stderr," %s entry: %s\n",WORKFILE,in_line);
		/* Skip any whitespace at beginning of the line: */
		char_addr = in_line;
		while(isspace(*char_addr)) {
			++char_addr;
		}
		// v20.1.1: Parse all lines whose 1st non-WS char is alphabetic; print "Ignoring (copy of workfile line)" for all entries not so.
		// NB: Discontinue support for numeric leading char, i.e. Mersenne-exponent-only legacy format:
		if(!isalpha(*char_addr)) {
			fprintf(stderr," Leading non-WS char of %s entry is not alphabetic ... skipping to next entry.\n",WORKFILE);
			goto read_next_assignment;
		}

		// Otherwise assume Prime95-style ini file format, with a possible modulus-specific leading keyword;
		// Default "Test=" means Mersenne, unless "Test" preceded by an explicit modulus-type string:
		MODULUS_TYPE = MODULUS_TYPE_MERSENNE;
		/* Re. the recently-added-to-Primenet PRP assignment type, On Dec 19, 2017, at 5:07 PM, George Woltman wrote:

		In "PRP=[aid],1,2,75869377,-1,75,0,3,4"		([aid] stands for an optional 32-hexit assignment ID)
			The first four (numeric) values are k,b,n,c as in modulus = k*b^n + c
			75 is how far factored
			0 is the number of PRP tests that will be saved if P-1 is done and finds a factor
			3 is the PRP base		(used for the first-time test)
			4 is the residue type	(used for the first-time test)

		You should only need to support residue types 1 and 5, which are identical if there are no cofactors
		(==> in JSON-result, return residue type 1).

		Example:	PRP=C42540C352E54E906108D48FA5D89488,1,2,80340397,-1,75,0,3,1
		means a PRP test of 1.2^80340397-1 = M80340397, TFed to 2^75 (ignore), 0 saved tests (ignore), PRP base 3 and type 1.

		EWM: The above are PRP-DCs ... first-time tests lack the last 2 numeric args above,
		and default to base = 3, residue type 1. Further, there is also a PRPDC= format for such.
		PRP-CF (Mersenne-cofactor-PRP tests) have same format as PRP, but append one or more known factors:

			PRP=[aid],1,2,3215747,-1,99,0,3,1,"4457025343,185822885311153245017"
		[EWM: F25 cofactor-test:
			PRP=n/a,1,2,33554432,+1,99,0,3,5,"2170072644496392193"
		]
		Ref. for residue type definitions: https://www.mersenneforum.org/showthread.php?p=468378#post468378
		George: "There are (at least) 5 PRP residue types:
																		EWM: Needed LR-binary-modpow modmul sequence, N = M(p) = 2^p-1 = [p binary ones],
																		initial seed x = a, and example final result for a = 3 and M(23) (not prime), M(31) (prime):
																													p = 23	p = 31	bc code:
			1: Fermat PRP. N is PRP if a^(N-1) = 1 mod N				p-1 (x = a.x^2), followed by one (x = x^2)	5884965	+1	m = 2^p-1; x = a = 3; for(i = 0; i < (p-2); i++) { x = a*x^2 % m; }; x = x^2 % m; x
			2: SPRP variant, N is PRP if a^((N-1)/2) = +/-1 mod N		p-1 (x = a.x^2)								7511964	-1	m = 2^p-1; x = a = 3; for(i = 0; i < (p-2); i++) { x = a*x^2 % m; }; x
			3: Fermat PRP variant. N is PRP if a^(N+1) = a^2 mod N		p   mod-squarings (x = x^2)					2633043	 9	m = 2^p-1; x = a = 3; for(i = 0; i < p; i++) { x = x^2 % m; }; x
			4: SPRP variant. N is PRP if a^((N+1)/2) = +/-a mod N		p-1 mod-squarings (x = x^2)					5758678	-3	m = 2^p-1; x = a = 3; for(i = 0; i < (p-1); i++) { x = x^2 % m; }; x
			5: Fermat PRP cofactor variant, N/d is PRP if a^(N-1) = a^d mod N/d

		I encourage programs to return type 1 residues as that has been the standard for prime95, PFGW, LLR for many years."

		EWM: Note for PRP w/Gerbicz-check we need a mod-squaring-only sequence, so we *compute* type 3, but then do a
		mod-div-by-scalar to get a^(N-1) mod N from a^2 mod N and *report* the resulting type 1 residue to the server.
		Example:
			Compute type 3 (Fermat-PRP) residue for N = M(23) using p = 23 mod-squarings of initial seed a = 3, yielding R3 = 2633043.
		Now we need to mod-divide by a^2 = 9 to get R1, which amounts to finding R1 such that a^2*R1 == R3 (mod N). 2 options here:

		1: Compute the inverse of a^2 (mod N) = 1864135, R1 = (a^-2)*R3 (mod N) = 5884965 .

		2: Find a small multiplier k such that R3 + k*N is divisible by a^2 = 9. Using trial and error starting from k = 0:
			k	R3 + k*N (mod 9)
			--	-------
			0	3
			1	7
			2	2
			3	6
			4	1
			5	5
			6	0, so R1 = (R3 + 6*N)/9 = 5884965 .
		More cleverly, we compute R3 == 3 (mod 9) and N == 4 (mod 9), so we want to solve 3 + 4*k == 0 (mod 9) for the index k.
		That can be done by simply rearranging to 4*k == -3 (mod 9), then finding inverse of 4 (mod 9) and multiplying both sides
		by that (mod 9) to get k. Note the difference between this mod-inverse computation and the one in option [1]: This one is
		not the mod-inverse w.r.to N, it's one w.r.to a^2, which is tiny compared to N.
			Calling that with args x = 4, n = 9 gives our k = -3*modinv(4,9) = -3*-2 = 6, same as the above trial-and-error approach.
		*/
		if((char_addr = strstr(in_line, "PRP")) != 0)	// This also handles the PRPDC= format
		{
			TEST_TYPE = TEST_TYPE_PRP;
			char_addr += 3;
			// Check [k,b,n,c] portion of in_line:
			cptr = check_kbnc(char_addr, &p);
			ASSERT(cptr != 0x0, "[k,b,n,c] portion of in_line fails to parse correctly!");
			// Next 2 entries in in_line are how-far-factored and "# of PRP tests that will be saved if P-1 is done and finds a factor":
			TF_BITS = 0xffffffff; tests_saved = 0.0;
			if((char_addr = strstr(cptr, ",")) != 0x0) {
				cptr++;
				// Only check if there's an appropriate TF_BITS entry in the input line
				TF_BITS = strtoul(++char_addr, &endp, 10);
				ASSERT((char_addr = strstr(cptr, ",")) != 0x0,"Expected ',' not found after TF_BITS field in assignment-specifying line!");	cptr++;
				tests_saved = strtod(++char_addr, &endp);
				if(tests_saved < 0 || tests_saved > 2) {
					sprintf(cbuf, "ERROR: the specified tests_saved field [%10.5f] should be in the range [0,2]!\n",tests_saved);	ASSERT(0,cbuf);
				}
				// char_addr now points to leftmost char of tests_saved field, which we will overwrite with 0;
				// endp points to to-be-appended leftover portion
			}
			pm1_done = (tests_saved == 0);
			// If there is still factoring remaining to be done, modify the assignment type appropriately.
			if(pm1_done) {	// pm1_done == TRUE is more or less a no-op, translating to "proceed with primality test"
				cptr = char_addr;	// ...but we do need to advance cptr past the ,TF_BITS,tests_saved char-block
			} else {
				// Create p-1 assignment, then edit original assignment line appropriately
				TEST_TYPE = TEST_TYPE_PM1;
				kblocks = get_default_fft_length(p);
				ASSERT(pm1_set_bounds(p, kblocks<<10, TF_BITS, tests_saved), "Failed to set p-1 bounds!");
				// Format the p-1 assignment into cbuf - use cptr here, as need to preserve value of char_addr:
				cptr = strstr(in_line, "=");	ASSERT(cptr != 0x0,"Malformed assignment!");
				cptr++;	while(isspace(*cptr)) { ++cptr; }	// Skip any whitespace following the equals sign
				if(is_hex_string(cptr, 32)) {
					strncpy(aid,cptr,32);	sprintf(cbuf,"Pminus1=%s,1,2,%llu,-1,%u,%llu\n",aid,p,B1,B2);	// If we get here, it's a M(p), not F(m)
				} else
					sprintf(cbuf,"Pminus1=1,2,%llu,-1,%u,%llu\n",p,B1,B2);
				// Copy up to the final (tests_saved) char of the assignment into cstr and append tests_saved = 0;
				// A properly formatted tests_saved field is 1 char wide and begins at the current value of char_addr:
				i = char_addr - in_line; strncpy(cstr,in_line, i); cstr[i] = '0'; cstr[i+1] = '\0';
				// Append the rest of the original assignment. If original lacked a linefeed, add one to the edited copy:
				strcat(cstr,endp);
				if(cstr[strlen(cstr)-1] != '\n') {
					strcat(cstr,"\n");
				}
				split_curr_assignment = TRUE;	// This will trigger the corresponding code following the goto:
				goto GET_NEXT_ASSIGNMENT;
			}	// First-time PRP test ... !cptr check is for assignments ending with [k,b,n,c] like "PRP=1,2,93018301,-1":
			if(!cptr || (char_addr = strstr(cptr, ",")) == 0x0) {
				PRP_BASE = 3;
				TEST_TYPE = TEST_TYPE_PRP;
			} else {	// PRP double-check:
				// NB: Hit a gcc compiler bug (which left i = 0 for e.g. char_addr = ", 3 ,...") using -O0 here ... clang compiled correctly, as did gcc -O1:
				i = (int)strtol(char_addr+1, &cptr, 10); // PRP bases other than 3 allowed; see https://github.com/primesearch/Mlucas/issues/18 //	ASSERT(i == 3,"PRP-test base must be 3!");
				PRP_BASE = i;
				ASSERT((char_addr = strstr(cptr, ",")) != 0x0,"Expected ',' not found in assignment-specifying line!");
				i = (int)strtol(char_addr+1, &cptr, 10); ASSERT(i == 1 || i == 5,"Only PRP-tests of type 1 (PRP-only) and type 5 (PRP and subsequent cofactor-PRP check) supported!");
				// Read in known prime-factors, if any supplied - resulting factors end up in KNOWN_FACTORS[]:
				if(*cptr == ',')						//vv--- Pass in unused file-ptr fq here in case function emits any messages:
					nfac = extract_known_factors(p,cptr+1);
				// Use 0-or-not-ness of KNOWN_FACTORS[0] to differentiate between PRP-only and PRP-CF:
				if(KNOWN_FACTORS[0] != 0ull) {
					ASSERT(i == 5,"Only PRP-CF tests of type 5 supported!");
					if (MODULUS_TYPE == MODULUS_TYPE_FERMAT) ASSERT(PRP_BASE == 3, "PRP-CF test base for Fermat numbers must be 3!");
				}
			}
			goto GET_EXPO;
		}
		else if((char_addr = strstr(in_line, "Fermat")) != 0)
		{
			char_addr += 6;
			/* Look for comma following the modulus keyword and position next-keyword search right after it: */
			if(!STREQN(char_addr,",",1))
				ASSERT(0,"Expected ',' not found in input following modulus type specifier!");
			else
				char_addr++;

			MODULUS_TYPE = MODULUS_TYPE_FERMAT;
			PRP_BASE = 2;	// v20: Pépin test doesn't use this as the initial seed (that defaults to 3), but rather for the random-shift
							// offsets used to prevent the shift count from modding to 0 as a result of repeated doublings (mod 2^m)
		}
		/* "Mersenne" is the default and hence not required, but allow it: */
		else if((char_addr = strstr(in_line, "Mersenne")) != 0)
		{
			char_addr += 8;
			/* Look for comma following the modulus keyword and position next-keyword search right after it: */
			if(!STREQN(char_addr,",",1))
				ASSERT(0,"Expected ',' not found in input following modulus type specifier!");
			else
				char_addr++;
		}

		// Pépin tests are assigned via "Fermat,Test=<Fermat number index>", so catch this clause by starting new if/else()
		if((char_addr = strstr(in_line, "Test")) != 0)
		{
			TEST_TYPE = TEST_TYPE_PRIMALITY;
			char_addr +=  4;
		}
		else if((char_addr = strstr(in_line, "DoubleCheck")) != 0)
		{
			TEST_TYPE = TEST_TYPE_PRIMALITY;
			char_addr += 11;
		}
	#if INCLUDE_TF
		else if((char_addr = strstr(in_line, "Factor")) != 0)
		{
			TEST_TYPE = TEST_TYPE_TF;
			char_addr +=  6;
		}
	#endif
		/* 11/23/2020: George W. re. p-1 assignment formats:
		"There is no documentation on the worktodo lines (except the source code):

			Pminus1=[aid,]k,b,n,c,B1,B2[,TF_BITS][,B2_start][,known_factors]	(,-separated list of known factors bookended with "")
			Pfactor=[aid,]k,b,n,c,TF_BITS,ll_tests_saved_if_factor_found

		(***EWM: Mlucas v20 converts PRP-with-(p-1)-needed assignments to paired Pminus1|PRP
			assignments to support the fused (p-1)-S1|First-part-of-PRP-test algorithm. ***)

		A tests_saved value of 0.0 will bypass any P-1 factoring. The PRP residue type is defined in primenet.h .
		AFAIK, the server issues Pfactor= lines not Pminus1= lines."
		*/
		// Use case-insensitive analog of strstr, stristr():
		else if((char_addr = stristr(in_line, "pminus1")) != 0)
		{
			TEST_TYPE = TEST_TYPE_PM1;
			char_addr += 7;
			// Check [k,b,n,c] portion of in_line:
			cptr = check_kbnc(char_addr, &p);
			ASSERT(cptr != 0x0, "[k,b,n,c] portion of in_line fails to parse correctly!");
			ASSERT((char_addr = strstr(cptr, ",")) != 0x0 ,"Expected ',' not found in assignment-specifying line!");
			B1 = (uint32)strtoul (char_addr+1, &cptr, 10);
			ASSERT((char_addr = strstr(cptr, ",")) != 0x0 ,"Expected ',' not found in assignment-specifying line!");
			/* The C11 standard re. strtoull: "On success the function returns the converted integer as unsigned long long int type
			and sets endPtr to point to the first character after the input number. On failure it returns 0 and sets endPtr to
			point to NULL. It handles integer overflows efficiently and return ULONG_LONG_MAX on overflow."
			However, in gdb-under-Mac testing of a decimal-digit input > 2^64, I found it returned ULONG_LONG_MAX, but
			also set endPtr to point to the first character after the input, which leaves some ambiguity - what if the
			input was in fact == ULONG_LONG_MAX? We assume here that nobody will use a p-1 stage bound so large:
			*/
			B2 = (uint64)strtoull(char_addr+1, &cptr, 10);	ASSERT(B2 != -1ull, "strtoull() overflow detected.");
			// Remaining args optional, with the 2 numerics presumed in-order, e.g. we only look for ',B2_start' field if ',TF_BITS' was present:
			if((char_addr = strstr(cptr, ",")) != 0x0) {
				TF_BITS = (int)strtoul(char_addr+1, &cptr, 10);	ASSERT(TF_BITS < 100 ,"TF_BITS value read from assignment is out of range.");
				if((char_addr = strstr(cptr, ",")) != 0x0) {
					B2_start = (uint64)strtoull(char_addr+1, &cptr, 10);	ASSERT(B2_start != -1ull, "strtoull() overflow detected.");
					if(B2_start > B1)	// It's a stage 2 continuation run
						s2_continuation = TRUE;
					// Read in known prime-factors, if any supplied - resulting factors end up in KNOWN_FACTORS[]:
					if(*cptr == ',') nfac = extract_known_factors(p,cptr+1);
				} else if((char_addr = strstr(cptr, "\"")) != 0x0) {	// Known-factors list need not be preceded by TF_BITS or B2_start
					nfac = extract_known_factors(p,cptr);	// cptr, not cptr+1 here, since need to preserve leading " bracketing factors-list
				}
			}
		}
		else if((char_addr = stristr(in_line, "pfactor")) != 0)	// Caseless substring-match as with pminus 1
		{
			TEST_TYPE = TEST_TYPE_PM1;
			PRP_BASE = 3;
			char_addr += 7;
			// Check [k,b,n,c] portion of in_line:
			cptr = check_kbnc(char_addr, &p);
			ASSERT(cptr != 0x0, "[k,b,n,c] portion of in_line fails to parse correctly!");
			ASSERT((char_addr = strstr(cptr, ",")) != 0x0 ,"Expected ',' not found in assignment-specifying line!");
			TF_BITS = (int)strtoul(char_addr+1, &cptr, 10);
			ASSERT((char_addr = strstr(cptr, ",")) != 0x0 ,"Expected ',' not found in assignment-specifying line!");
			tests_saved = strtod(++char_addr, &endp);
			if(tests_saved < 0 || tests_saved > 2) {
				sprintf(cbuf, "ERROR: the specified tests_saved field [%10.5f] should be in the range [0,2]!\n",tests_saved);	ASSERT(0,cbuf);
			}
			ASSERT(pm1_set_bounds(p, get_default_fft_length(p)<<10, TF_BITS, tests_saved), "Failed to set p-1 bounds!");
		}
	#if INCLUDE_ECM
		else if(strstr(char_addr, "ECM"))
		{
			TEST_TYPE = TEST_TYPE_ECM;
			char_addr +=  3;
		}
	#endif
		else
		{
			snprintf(cbuf,STR_MAX_LEN*2,"WARN: Unrecognized/Unsupported option or empty assignment line. The ini file entry was %s\n",in_line);
			fprintf(stderr,"%s",cbuf);
			goto read_next_assignment;
		}

		if(!p) {	// For legacy assignment types, set p here
			ASSERT((char_addr = strstr(char_addr, "=")) != 0x0,"Expected '=' not found in assignment-specifying line!");
			char_addr++;
			/* Skip any whitespace following the equals sign:*/
			while(isspace(*char_addr)) { ++char_addr; }
			/* Check for a 32-hex-digit PrimeNet v5 assignment ID preceding the exponent: */
			if(is_hex_string(char_addr, 32))
				char_addr += 33;
			else if(STREQN_NOCASE(char_addr,"n/a",3))
				char_addr = strstr(char_addr, ",") + 1;

			p = strtoull(char_addr, &cptr, 10);	ASSERT(p != -1ull, "strtoull() overflow detected.");
		}

	GET_EXPO:
		// Need to init this for savefile-naming code
		ASSERT(p != 0ull, "Exponent has not been set!");
		sprintf(ESTRING,"%llu",p);

		// In PRP-test case, have already read the exponent from the worktodo line
		/* Special case of user forcing a non-default FFT length for an exponent in the worktodo file: */
		if(exponent && (p != exponent)) {	// || (MODULUS_TYPE != MODULUS_TYPE_MERSENNE))	15. Oct 2012: Need same flexibility for Fermat numbers (e.g. F27 @ 7168k) as for Mersennes, so disable modulus-type part of conditional
			sprintf(cbuf,"User-supplied exponent and FFT-length for full-length test requires an exponent-matching 'Test=<exponent>' or 'DoubleCheck=<exponent>' %s entry!",WORKFILE);
			ASSERT(0,cbuf);
		}

		/* Check #bits in the Mersenne exponent vs. the allowed maximum: */
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			nbits_in_p = 64 - leadz64(p);
		}
		/* If it's a Fermat number, need to check size of 2^ESTRING: */
		else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
			findex = (uint32)p;
			if(findex <= MAX_PRIMALITY_TEST_BITS)
				p = (uint64)1 << findex;
			else
				ASSERT(0,"nbits_in_p <= MAX_PRIMALITY_TEST_BITS");
			// For purposes of the bits-in-p limit, treat 2^findex as having (findex) rather than (findex+1) bits:
			nbits_in_p = findex;
		}
		else
			ASSERT(0,"MODULUS_TYPE unknown!");

		ASSERT(nbits_in_p <= MAX_EXPO_BITS,"Require nbits_in_p <= MAX_EXPO_BITS");

	#if INCLUDE_TF

	  INIT_TF:

		/* If nbits_in_p > MAX_PRIMALITY_TEST_BITS, it better be a TF run: */
		if(TEST_TYPE == TEST_TYPE_TF) {
			/* Currently TF only supported for Mersennes: */
			if(MODULUS_TYPE != MODULUS_TYPE_MERSENNE) {
				sprintf(cbuf, "ERROR: Trial-factoring Currently only supported for Mersenne numbers. The ini file entry was %s\n",in_line);
				fprintf(stderr,"%s",cbuf);
				goto GET_NEXT_ASSIGNMENT;
			}

			/* For now, always start at k = 1: */
			log2_min_factor = 0.0;
			log2_max_factor = get_default_factoring_depth(p);
			ASSERT(log2_max_factor <= MAX_FACT_BITS, "log2_max_factor > MAX_FACT_BITS!");

			/* Field following the exponent is the already-factored-to depth: if none found, use defaults. */
			char_addr = strstr(char_addr, ",");
			if(char_addr++) {
				/* Convert the ensuing numeric digits to ulong: */
				TF_BITS = strtoul(char_addr, &endp, 10);
				/* Specified already-factored-to depth is larger than default factor-to depth - no more factoring to be done. */
				if(TF_BITS > log2_max_factor) {
					sprintf(cbuf, "INFO: the specified already-factored-to depth of %u bits exceeds the default %10.4f bits - no more factoring to be done.\n", TF_BITS, log2_max_factor);
					fprintf(stderr,"%s",cbuf);
					goto GET_NEXT_ASSIGNMENT;
				}
			}
			/* Mode to allow user to override normal automated-test TF default limits:
			If a second ,{#bits} argument is present in the input line, then this user-set
			desired factoring depth overrides the normal default for the exponent in question.
			In this mode, warn if TF-to depth greater than automated-mode default, but allow: */
			if(char_addr)
				char_addr = strstr(char_addr, ",");
			if(char_addr++) {
				bit_depth_todo = strtoul(char_addr, &endp, 10);
				if(bit_depth_todo > MAX_FACT_BITS) {
					sprintf(cbuf, "ERROR: factor-to bit_depth of %u > max. allowed of %u. The ini file entry was %s\n",TF_BITS,MAX_FACT_BITS,in_line);
					fprintf(stderr,"%s",cbuf);
					goto GET_NEXT_ASSIGNMENT;
				}
				else if(bit_depth_todo <= TF_BITS) {
					sprintf(cbuf, "ERROR: factor-to bit_depth of %u < already-done depth of %u. The ini file entry was %s\n",TF_BITS,TF_BITS,in_line);
					fprintf(stderr,"%s",cbuf);
					goto GET_NEXT_ASSIGNMENT;
				}
				else if(bit_depth_todo > log2_max_factor) {
					sprintf(cbuf, "WARN: the specified factor-to depth of %u bits exceeds the default %10.4f bits - I hope you know what you're doing.\n",TF_BITS,log2_max_factor);
					log2_max_factor = bit_depth_todo;
					fprintf(stderr,"%s",cbuf);
				}
				log2_max_factor = bit_depth_todo;
			}
		}
		else if(nbits_in_p > MAX_PRIMALITY_TEST_BITS)
		{
			sprintf(cbuf, "ERROR: Inputs this large only permitted for trial-factoring. The ini file entry was %s\n",in_line);
			fprintf(stderr,"%s",cbuf);
			goto GET_NEXT_ASSIGNMENT;
		}
	#endif 	// INCLUDE_TF

		/* If "Test..." or "DoubleCheck", check for TF_BITS and pm1_done fields following the = sign:
		if present and there is still factoring remaining to be done, modify the assignment type appropriately.
		Assignment format:
			<Test|DoubleCheck>=[aid],p,TF_BITS,pm1_done?
		*/
		if(TEST_TYPE == TEST_TYPE_PRIMALITY) {
			/* TF_BITS: */
			TF_BITS = 0xffffffff;	/* Only check if there's an appropriate entry in the input line */
			char_addr = strstr(char_addr, ",");
			if(char_addr++) {
				/* Convert the ensuing numeric digits to ulong: */
				TF_BITS = strtoul(char_addr, &endp, 10);
			#if INCLUDE_TF
				if(TF_BITS > MAX_FACT_BITS) {
					snprintf(cbuf,STR_MAX_LEN*2,"ERROR: TF_BITS of %u > max. allowed of %u. The ini file entry was %s\n", TF_BITS, MAX_FACT_BITS, in_line);
					fprintf(stderr,"%s",cbuf);
					goto GET_NEXT_ASSIGNMENT;
				}
				/* If TF_BITS less than default TF depth for this exponent
				and this platform is approved for trial factoring, switch task type to TF:
				*/
				log2_max_factor = get_default_factoring_depth(p);
				if(TF_BITS < log2_max_factor) {
					TEST_TYPE = TEST_TYPE_TF;
					// For now, always start at k = 1:
					log2_min_factor = 0.0;
				}
				else if(TF_BITS > log2_max_factor) {
					sprintf(cbuf, "WARN: the specified already-factored-to depth of %u bits exceeds the default %10.4f bits - no more factoring to be done.\n", TF_BITS, log2_max_factor);
					fprintf(stderr,"%s",cbuf);
				}
			#endif
				// p-1 factoring already done?
				char_addr = strstr(char_addr, ",");
				if(char_addr++) {
					/* Convert the ensuing numeric digits to ulong: */
					pm1_done = strtoul(char_addr, &endp, 10);
					if(pm1_done > 1) {
						sprintf(cbuf, "ERROR: the specified pm1_done field [%u] should be 0 or 1!\n",pm1_done);
						ASSERT(0,cbuf);
					}
					if(!pm1_done) {	// pm1_done == TRUE is a no-op, translating to "proceed with primality test"
						// Don't actually use this in pm1_set_bounds(), due to the rise of the single-shot PRP-with-proof paradigm, but for form's sake:
						tests_saved = 1;
						// Create p-1 assignment, then edit original assignment line appropriately
						TEST_TYPE = TEST_TYPE_PM1;
						kblocks = get_default_fft_length(p);
						ASSERT(pm1_set_bounds(p, kblocks<<10, TF_BITS, tests_saved), "Failed to set p-1 bounds!");
						// Format the p-1 assignment into cbuf:
						char_addr = strstr(in_line, "=");	ASSERT(char_addr != 0x0,"Malformed assignment!");
						char_addr++;	while(isspace(*char_addr)) { ++char_addr; }	// Skip any whitespace following the equals sign
						if(is_hex_string(char_addr, 32)) {
							strncpy(aid,char_addr,32);	sprintf(cbuf,"Pminus1=%s,1,2,%llu,-1,%u,%llu\n",aid,p,B1,B2);	// If we get here, it's a M(p), not F(m)
						} else
							sprintf(cbuf,"Pminus1=1,2,%llu,-1,%u,%llu\n",p,B1,B2);

						// Copy all but the final (pm1_done) char of the assignment into cstr and append pm1_done = 1. If in_line ends with newline, first --j:
						j = strlen(in_line) - 1;	j -= (in_line[j] == '\n');
						strncpy(cstr,in_line,j);	cstr[j] = '\0';	strcat(cstr,"1\n");
						split_curr_assignment = TRUE;	// This will trigger the corresponding code following the goto:
						goto GET_NEXT_ASSIGNMENT;
					}
				}
			}
		}	/* if(TEST_TYPE == TEST_TYPE_PRIMALITY) */
		if(fp) {	// Close workfile, if still open
			fclose(fp); fp = 0x0;
		}
	}	// endif(fp)
	/****************************************/
	/* Self-test or User-supplied exponent: */
	/****************************************/
	else if(exponent != 0)	/* elseif((found WORKFILE) == FALSE) */
	{
		p = exponent;
		fprintf(stderr," %s file not found...using user-supplied command-line exponent p = %llu\n",WORKFILE,p);
		/* This takes care of the number-to-char conversion and leading-whitespace-removal
		in one step - use PSTRING for temporary storage here: */
		strcpy(ESTRING, &PSTRING[convert_uint64_base10_char(PSTRING, p)]);

		/* If it's a Fermat number, get the real exponent: */
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
			findex = (uint32)p;
			if(findex <= MAX_PRIMALITY_TEST_BITS)
				p = (uint64)1 << findex;
			else
				ASSERT(0,"nbits_in_p <= MAX_PRIMALITY_TEST_BITS");

			/* For purposes of the bits-in-p limit, treat 2^findex as having
			(findex) rather than (findex+1) bits: */
			nbits_in_p = findex;
		} else
			nbits_in_p = 64 - leadz64(p);

		INTERACT=TRUE;

		ASSERT(TEST_TYPE,"TEST_TYPE not set!");
		ASSERT(TEST_TYPE <= TEST_TYPE_MAX,"TEST_TYPE out of range!");

		/* If nbits_in_p > MAX_PRIMALITY_TEST_BITS, it better be a TF run: */
		if(TEST_TYPE == TEST_TYPE_TF)
		{
		#if INCLUDE_TF
			/* Currently TF only supported for Mersennes: */
			ASSERT((MODULUS_TYPE == MODULUS_TYPE_MERSENNE), "Trial-factoring Currently only supported for Mersenne numbers");
			/* For now, always start at k = 1: */
			log2_min_factor = 0.0;
			if(iterations) {
				log2_max_factor = iterations;
				iterations = 0;
			}
			else
				log2_max_factor = get_default_factoring_depth(p);

			ASSERT(log2_max_factor >=             0, "log2_max_factor must be positive!");
			ASSERT(log2_max_factor <= MAX_FACT_BITS, "log2_max_factor exceeds MAX_FACT_BITS!");
		#else
			ASSERT(0, "Trial-factoring not supported for this build/platform.");
		#endif
		}
		else if(TEST_TYPE == TEST_TYPE_PM1)	/* P-1 factoring attempt */
		{
			ASSERT(nbits_in_p <= MAX_PRIMALITY_TEST_BITS, "Inputs this large only permitted for trial-factoring.");
			pm1_check_bounds();
			// Proper setting of timing_test_iters in this case needs us to compute the stage 1 prime-powers product:
			// Compute stage 1 prime-powers product, store in PM1_S1_PRODUCT and store #bits of same in PM1_S1_PROD_BITS:
			s1p_alloc = compute_pm1_s1_product(p);
			RES_SHIFT = 0ull;	// Must set = 0 here to make sure BASE_MULTIPLIER_BITS array gets set = 0 below
			iterations = PM1_S1_PROD_BITS;
			if(iterations > MAX_SELFTEST_ITERS) {
				fprintf(stderr, " Stage 1: %u iterations resulting from bound b1 = %u exceeds self-test limit of %u.\n",iterations,B1,MAX_SELFTEST_ITERS);
				return ERR_TESTITERS_OUTOFRANGE;
			}
			timing_test_iters = iterations;
		}
		else	/* Primality or PRP test */
		{
		/*	fprintf(stderr, "P = %u, nbits_in_p = %d\n",p,nbits_in_p);	*/
			ASSERT(nbits_in_p <= MAX_PRIMALITY_TEST_BITS, "Inputs this large only permitted for trial-factoring.");
			ASSERT(iterations != 0,"Timing test with User-supplied exponent requires number of iterations to be specified via the -iters flag!");
			if(iterations <= 0) {
				fprintf(stderr, " Specified %u self-test iterations : must be > 0.\n", iterations);
				return ERR_TESTITERS_OUTOFRANGE;
			} else if(iterations > MAX_SELFTEST_ITERS) {
				fprintf(stderr, " Specified %u iterations exceeds self-test limit of %u.\n", iterations, MAX_SELFTEST_ITERS);
				return ERR_TESTITERS_OUTOFRANGE;
			}
			timing_test_iters = iterations;
		}
	} else {
		fprintf(stderr,"No %s file not found, nor user-supplied command-line exponent.\n",WORKFILE);
		print_help();
		ASSERT(0, "Unsupported combination of command-line args. Note that if you are trying to\nrun a single-FFT-length self-test, you *must* explicitly specify the iteration\ncount, e.g. './Mlucas -fft 7168 <-iters [+int]> [-cpu <args>]'");
	}	// endif(found WORKFILE?)

	// If production run (not self-test), echo assignment to per-exponent logfile:
	if(!INTERACT) {
		snprintf(cbuf,STR_MAX_LEN*2," %s entry: %s\n",WORKFILE,in_line);
		mlucas_fprint(cbuf,0);
	}

/*************************************************************************************/

	/*...If a valid assignment, look for a matching restart file.
	The first character of the restart file name indicates the
	assignment type: 'p' for primality test and p-1 factoring, (with 'q' for
	the backup restart file in these cases), 't' for trial-factoring.
	*/
	/* gcc with optimizations turned on wasn't initing all elements of restart file names = \0, so insert one manually after the p and q, before calling strcat() */
	RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f'); RESTARTFILE[1] = '\0';
	strcat(RESTARTFILE, ESTRING);
	/* The statfile for a given exponent is 'p{exponent}.stat'irrespective of assignment type: */
	strcpy(STATFILE, RESTARTFILE);
	strcat(STATFILE, ".stat");
	/*fprintf(stderr, "STATFILE = %s\n",STATFILE);	*/
	ASSERT(TEST_TYPE,"TEST_TYPE not set!");
	ASSERT(TEST_TYPE <= TEST_TYPE_MAX,"TEST_TYPE out of range!");

	/* Fom this point onward the first character of restart filenames is context-dependent: */
#if INCLUDE_TF
	if(TEST_TYPE == TEST_TYPE_TF)
	{
		/* Do any needed Trial-Factoring and then update the worktodo file appropriately: */
		/* Let the factoring module handle the restart-file processing in this case. */
		if(!ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE][TEST_TYPE_TF])
		{
			sprintf(cbuf, "TEST_TYPE_TF with MODULUS_TYPE = %u not supported!\n", MODULUS_TYPE);
			ASSERT(0, cbuf);
		}

		factor(ESTRING, log2_min_factor, log2_max_factor);
		goto GET_NEXT_ASSIGNMENT;
	} else
#endif
	if(TEST_TYPE > TEST_TYPE_MAX)
	{
		ASSERT(0,"ERROR: Unrecognized assignment type in savefile processing.\n");
	}
	/* endif(TEST_TYPE == ...) */

/********************* P-1, primality, or PRP Test: ***********************************************/

	if(p < PMIN) {
		fprintf(stderr, " p must be at least %llu.\n",PMIN);
		return ERR_EXPONENT_ILLEGAL;
	} else if(p > PMAX) {
		fprintf(stderr, " p must be no greater than %llu.\n",PMAX);
		return ERR_EXPONENT_ILLEGAL;
	}

	// Only check production-run exponents here - self-tests already have an iteration limit < 2^32:
	if(TEST_TYPE == TEST_TYPE_PRIMALITY || TEST_TYPE == TEST_TYPE_PRP) {
		if(p > (1ull << 32) && !INTERACT) {	// No point emitting the warning for self-tests since those are guaranteed to have (maxiter < MAX_SELFTEST_ITERS)
			fprintf(stderr,"Full-length LL/PRP/Pepin tests on exponents this large not supported; will exit at self-test limit of %u.\n",MAX_SELFTEST_ITERS);
			maxiter = MAX_SELFTEST_ITERS;
		} else
			maxiter = (uint32)p;
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		/* TODO: if/when TF code is hooked in, make its small-primes table a global and
		init it at the start of execution, then use it to trial-divide p up to sqrt(p):
		*/
		/* Allow non-prime-exponents for p-1 test, but not for LL: */
		if(!isPRP64(p) && ((TEST_TYPE == TEST_TYPE_PRIMALITY) || (TEST_TYPE == TEST_TYPE_PRP))) {
			fprintf(stderr, "Mersenne Primality/PRP-tests require a prime exponent!\n");
			return ERR_EXPONENT_ILLEGAL;
		}
		// PRP implies Mersenne-mod, base must be > 1 and not a power of 2:
		if(TEST_TYPE == TEST_TYPE_PRP && (PRP_BASE < 2 || isPow2(PRP_BASE))) {
			fprintf(stderr, "PRP-test requires base > 1 and not a power of 2!\n");
			return ERR_EXPONENT_ILLEGAL;
		}

		TRANSFORM_TYPE = REAL_WRAPPER;
		snprintf(PSTRING,STR_MAX_LEN, "M%s", ESTRING);
		/* v19:
		Unlike standard mod-M(p) Fermat-PRP test, x0^(N-1) ?== 1 (mod N) which for N = M(p) gives N-1 = 2^p-2
		= 0b111[p-1 binary 1s]1110 and thus requires [p-2 (x := x^2*base) steps followed by 1 final squaring], the
		Gerbicz check requires a pure sequence of squarings. RG: simply replace the standard Fermat-PRP test with a
		modified one where we add 2 to the computed power and check whether the result == x0^2 (mod n).
		Thus instead of the standard base-x0 Fermat PRP test, [x0^(N-1) = x0^(2^p-2) == 1 (mod N)], we check whether
			x0^(n+1) = x0^(2^p) == x0^2 (mod n), which needs p mod-squarings.
		For a Pepin test of N = F(m) = 2^(2^m) + 1 we already have a (2^m-1)-pure-squaring sequence, no changes needed:
		*/
		if(TEST_TYPE == TEST_TYPE_PRIMALITY) {	// LL does p-2 iterations; Fermat-PRP modified as described in above commentary does p
			maxiter -= 2;
		} else if(TEST_TYPE == TEST_TYPE_PRP) {
			RES_SHIFT = 0ull;	// Must set = 0 here to make sure BASE_MULTIPLIER_BITS array gets set = 0 below
		} else if(TEST_TYPE == TEST_TYPE_PM1) {
			// Compute stage 1 prime-powers product, store in PM1_S1_PRODUCT, store #bits of same in PM1_S1_PROD_BITS:
			pm1_check_bounds();
			s1p_alloc = compute_pm1_s1_product(p);
			maxiter = PM1_S1_PROD_BITS;	// NOTE: In this case we don't want to override the PRP_BASE = 3 value set in compute_pm1_s1_product()
			ASSERT(B1 > 0 && maxiter > B1, "P-1 b1 and/or maxiter unset!");
			RES_SHIFT = 0ull;	// Must set = 0 here to make sure BASE_MULTIPLIER_BITS array gets set = 0 below
		} else
			ASSERT(0,"Unsupported test type! (Neither LL,PRP nor P-1)");
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
	#ifdef USE_ARM_V8_SIMD
		ASSERT(0, "ARMv8 SIMD builds do not support Fermat-number testing!");
	#endif
		ASSERT(findex >= 13 && findex < 64, "Fermat number index must be in range [13,63]!\n");
		// This takes care of the number-to-char conversion and leading-whitespace-removal
		// in one step - use PSTRING for temporary storage here:
		strcpy(ESTRING, &PSTRING[convert_uint64_base10_char(PSTRING, (uint64)findex)]);
		ASSERT((p >> findex) == 1,"Require (p >> findex) == 1");
		sprintf(BIN_EXP,"%llu",p);	// May need this for workfile postprocessing if assignment is in KBNC format
		TRANSFORM_TYPE = RIGHT_ANGLE;
		sprintf(PSTRING, "F%u", findex);
		if(TEST_TYPE == TEST_TYPE_PRIMALITY) {
			maxiter -= 1;	// Pepin-test does p-1 iterations
			PRP_BASE = 2;	// v20: Pépin test doesn't use this as the initial seed (that defaults to 3), but rather for the random-shift
							// offsets used to prevent the shift count from modding to 0 as a result of repeated doublings (mod 2^m)
		} else if(TEST_TYPE == TEST_TYPE_PRP) {
			ASSERT(KNOWN_FACTORS[0] != 0, "Fermat-mod PRP test implies a PRP-CF run, but no known-factors provided!");
			RES_SHIFT = 0ull;	// Must set = 0 here to make sure BASE_MULTIPLIER_BITS array gets set = 0 below
		} else if(TEST_TYPE == TEST_TYPE_PM1) {
			// Compute stage 1 prime-powers product, store in PM1_S1_PRODUCT, store #bits of same in PM1_S1_PROD_BITS:
			pm1_check_bounds();
			s1p_alloc = compute_pm1_s1_product(p);
			maxiter = PM1_S1_PROD_BITS;	// NOTE: In this case we don't want to override the PRP_BASE = 3 value set in compute_pm1_s1_product()
			ASSERT(B1 > 0 && maxiter > B1, "P-1 b1 and/or maxiter unset!");
			RES_SHIFT = 0ull;	// Must set = 0 here to make sure BASE_MULTIPLIER_BITS array gets set = 0 below
		} else
			ASSERT(0,"Unsupported test type! (Neither Pepin-primality nor P-1)");

		j = ((ITERS_BETWEEN_CHECKPOINTS+63) >> 6);
		if(RES_SHIFT == 0ull) {
			for(i = 0; i < j; i++) { BASE_MULTIPLIER_BITS[i] = 0ull; }
		} else {
			// In order to avoid shift = 0 after [findex] repeated doublings (mod 2^findex), try using a randomized
			// bitmap, 1 bit read each squaring, if bit = 1 we multiply the residue by 2 and increment the shift by 1:
			// "Prime' the RNG by doing [24-bit random] calls:
		#if 0	// This is problematic in terms of nixing the kind of repeatable-deterministic-pseudorandomness needed for debug
			k = (int)time(NULL) & 0x00ffffff;	itmp64 = 0ull;
			printf("Shift = %u: Priming RNG with %u rand() calls...\n",RES_SHIFT,k);
			for(i = 0; i < k; i++) { itmp64 += rng_isaac_rand(); }	// Use rand() outout to increment itmp64, to eliminate "expression result unused" warns
		#endif
			for(i = 0; i < j; i++) { BASE_MULTIPLIER_BITS[i] = rng_isaac_rand(); }; i = ITERS_BETWEEN_CHECKPOINTS&63;	// i = #low bits used in high limb
			printf("random-bit-array popcount = %u\n",mi64_popcount(BASE_MULTIPLIER_BITS,j) - popcount64(BASE_MULTIPLIER_BITS[j-1] >> i));
		}
	}
	else {
		ASSERT(0,"Unknown Self-Test Modulus Type!");
	}	/* endif(MODULUS_TYPE) */

	// mi64_shlc currently limited to 32-bit shift counts - for technical reasons described in comments at top of that function,
	// the largest exponent testable-with-shift must satisfy condition below, which yields largest M(p) with p = 4294967231 = 2^32-65:
	if(RES_SHIFT && (p+63) > 0xFFFFFFFFull) {
		sprintf(cbuf,"ERROR: Exponents this large do not support residue shift! Please run with '-shift 0'.\n");
		ASSERT(0,cbuf);
	}
	/* In production-run (INTERACT = False) mode, allow command-line-forced FFT lengths which are at most
	"one size too large" relative to the default length for the exponent in question. Supported lengths
	are of the form [8,9,10,11,12,13,14,15]*2^k, so base our acceptance threshold on the largest adjacent-
	FFT-lengths ratio permitted under this schema, 9/8:
	*/
	kblocks = get_default_fft_length(p);
	if(!fft_length || (!INTERACT && MODULUS_TYPE == MODULUS_TYPE_MERSENNE && 8*fft_length > 9*kblocks)) {
		if(!kblocks) {
			fprintf(stderr,"ERROR detected in get_default_fft_length for p = %llu.\n",p);
			return ERR_FFTLENGTH_ILLEGAL;
		}
	} else {
		kblocks = fft_length;
		if(kblocks <= 0) {
			fprintf(stderr, "ERROR: %d K must be >= 0.\n",kblocks);
			return ERR_FFTLENGTH_ILLEGAL;
		}
	}
	/*...calculate the unpadded runlength...	*/
	n = kblocks << 10;

	/* Make sure the FFT length and radix set are supported: */
	if(INTERACT)
	{
		if(radix_set < 0) {
			fprintf(stderr, " Specified radix set %u for self-test : must be >= 0.\n", radix_set);
			return ERR_RADIXSET_UNAVAILABLE;
		}

		/* Now that we have the preferred radix set index, get the corresponding FFT radices: */
		i = get_fft_radices(kblocks, radix_set, &NRADICES, RADIX_VEC, 10);
		if(i == ERR_FFTLENGTH_ILLEGAL)
		{
			fprintf(stderr, "ERROR: length %d = %d K not available.\n",n,kblocks);
			return ERR_FFTLENGTH_ILLEGAL;
		}
		if(i == ERR_RADIXSET_UNAVAILABLE)
		{
			fprintf(stderr, " Specified radix set %u for self-test unavailable.\n", radix_set);
			return ERR_RADIXSET_UNAVAILABLE;
		}

		if(timing_test_iters > maxiter) {
			fprintf(stderr, " This exceeds the primality-test limit; will perform %u iterations for timing test.\n",maxiter);
			timing_test_iters = maxiter;
		}
		else if(maxiter > timing_test_iters)
			maxiter = timing_test_iters;
	}
	else
	{
	SETUP_FFT:
		/* Look for a best-FFT-radix-set entry in the .cfg file: */
		dum = get_preferred_fft_radix(kblocks);
		if(!dum) {	// Need to run a timing self-test at this FFT length before proceeding:
			sprintf(cbuf, "INFO: FFT length %d = %d K not found in the '%s' file.\n", n, kblocks, CONFIGFILE);
			fprintf(stderr, "%s", cbuf); // Extra information on the default FFT selected. The following line allows the FFT to be overridden for Fermat exponents; see https://github.com/primesearch/Mlucas/pull/11
			if (!fft_length || MODULUS_TYPE == MODULUS_TYPE_MERSENNE) return ERR_RUN_SELFTEST_FORLENGTH + (kblocks << 8);
		}
		else if(dum != kblocks)
		{
			/* If return value != kblocks, extract the FFT length it encodes: */
			i = extractFFTlengthFrom32Bit(dum);
			/* Only allow lengths that are <= 2x default */
			if( !(i >= kblocks && i <= (kblocks<<1) ) )
			{
				sprintf(cbuf,"Call to get_preferred_fft_radix returns out-of-range FFT length: asked for %u, returned %u, packed value= 0x%8X\n", kblocks, i, dum);
				ASSERT(0, cbuf);
			}
			else	/* If length acceptable, extract the FFT-radix data encoded and populate the NRADICES and RADIX_VEC[] globals */
			{
				extractFFTradicesFrom32Bit(dum);
				kblocks = i;
				/* Make sure the FFT length is supported: */
				if(get_fft_radices(kblocks, 0, 0x0, 0x0, 0) != 0)
				{
					ASSERT(get_fft_radices(kblocks, 0, 0x0, 0x0, 0) == ERR_FFTLENGTH_ILLEGAL, "Unexpected return value for get_fft_radices()");
					sprintf(cbuf, "ERROR: length %d = %d K not available.\n",n,kblocks);
					fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
				}
			}
		}
	}

	/*...calculate the unpadded runlength...	*/
	n = kblocks << 10;
	if(kblocks != (n>>10)) {
		fprintf(stderr, "ERROR: length %d K overflows N = 1024*K.\n",kblocks);
		return ERR_FFTLENGTH_ILLEGAL;
	}

	/* If specified FFT length smaller than default for this exponent [only an issue for user-specified FFT length],
	print a warning if the p/pmax ratio > 1 to an acceptably small degree; error out if the ratio is unreasonably > 1:
	*/
	uint64 pmax_rec = given_N_get_maxP(n);	double exp_ratio =  (double)p/pmax_rec;
	fprintf(stderr, "INFO: Maximum recommended exponent for FFT length (%u Kdbl) = %llu; p[ = %llu]/pmax_rec = %12.10f.\n",kblocks,pmax_rec,p,exp_ratio);
	// Set initial value of USE_SHORT_CY_CHAIN based on how close p/pmax is to 1.0, but only if current chain length is longer
	// (e.g. if ROE-retry logic has led to a shorter-than-default chain length, don't revert to default):
	if(exp_ratio > 0.99 && USE_SHORT_CY_CHAIN < 3)
		USE_SHORT_CY_CHAIN = 3;
	else if(exp_ratio > 0.98 && USE_SHORT_CY_CHAIN < 2)
		USE_SHORT_CY_CHAIN = 2;
	else if(exp_ratio > 0.97 && USE_SHORT_CY_CHAIN < 1)
		USE_SHORT_CY_CHAIN = 1;
	const char*arr_sml[] = {"long","medium","short","hiacc"};
	fprintf(stderr,"Initial DWT-multipliers chain length = [%s] in carry step.\n",arr_sml[USE_SHORT_CY_CHAIN]);
	// v20: If exp_ratio > 0.98, set ITERS_BETWEEN_CHECKPOINTS = 10000 irrespective of #threads or user-forced greater value;
	if(USE_SHORT_CY_CHAIN)
		ITERS_BETWEEN_CHECKPOINTS = MIN(ITERS_BETWEEN_CHECKPOINTS,10000);

	if(kblocks < (i = get_default_fft_length(p))) {
		/* If it's at least close, allow it but print a warning; otherwise error out: */
		if((i << 10) == get_nextlarger_fft_length(n) && exp_ratio < 1.05)
			fprintf(stderr, "INFO: specified FFT length %d K is less than recommended %d K for this p.\n",kblocks,i);
		else {
			sprintf(cbuf, "ERROR: specified FFT length %d K is much too small: Recommended length for this p = %d K ... quitting.\n",kblocks,i);
			ASSERT(0, cbuf);
		}
	}

	// Set the array padding parameters - only use array padding elements for runlengths > 32K:
	if(kblocks > 32) {
		DAT_BITS = DAT_BITS_DEF;
		PAD_BITS = PAD_BITS_DEF;
	} else {	// This causes the padding to go away:
		DAT_BITS =31;
		PAD_BITS = 0;
	}
	/*...If array padding turned on, check that the blocklength divides the unpadded runlength...	*/
	if((DAT_BITS < 31) && ((n >> DAT_BITS) << DAT_BITS) != n)
		ASSERT(0,"ERROR: blocklength does not divide runlength!");

	/*...Find padded array length...	*/
	npad = n + ( (n >> DAT_BITS) << PAD_BITS );	/* length of padded data array.	*/
	/* If the residue and other modulus-size-dependent data arrays too small for the new assignment, deallocate them: */
	if(nalloc > 0 && npad > nalloc)
	{
		ASSERT(a_ptmp != 0x0 && a != 0x0 && b != 0x0 && c != 0x0 && d != 0x0,"Require (a_ptmp,a,b,c,d) != 0x0");
		free((void *)a_ptmp); a_ptmp = a = b = c = d = e = 0x0; b_uint64_ptr = c_uint64_ptr = d_uint64_ptr = e_uint64_ptr = 0x0;
		free((void *)arrtmp); arrtmp=0x0;
		free((void *)BIGWORD_BITMAP);	BIGWORD_BITMAP = 0x0;
		free((void *)BIGWORD_NBITS);	BIGWORD_NBITS = 0x0;
		nbytes = nalloc = 0;
	}

	if(nalloc == 0) {
		// If it's a multi-exponent self-test, alloc to the maximum FFT length which will be used:
		if(maxFFT > kblocks) {
			i = (maxFFT << 10);	npad = i + ( (i >> DAT_BITS) << PAD_BITS );
		}	/*** MUST PRESERVE I THROUGH ALLOC OF ARRTMP[] BELOW ***/
		j = 0;
		if(npad & 7)
			j = 8 - (npad & 7);
		nalloc = npad + j;	ASSERT((nalloc & 7) == 0,"nalloc must be a multiple of 8!");	// This is so b,c,d enjoy same 64-byte alignment as a[]
		nbytes = nalloc<<3;
		ASSERT(a_ptmp == 0x0 && a == 0x0 && b == 0x0 && c == 0x0 && d == 0x0 && e == 0x0 && arrtmp == 0x0,"Require (a_ptmp,b,c,d,e,arrtmp) == 0x0");
		if(use_lowmem == 2) {	// Handy for huge-FFT self-tests on low-mem systems
			sprintf(cbuf,"WARN: Low-memory[%u] run mode disallows PRP-testing|Gerbicz-check and p-1 stage 2.\n",use_lowmem);
			mlucas_fprint(cbuf,1);
			j = 1;
		} else {
			j = 5;
		}
		a_ptmp = ALLOC_DOUBLE(a_ptmp, j*nalloc);	if(!a_ptmp){ sprintf(cbuf, "ERROR: unable to allocate array A in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		a      = ALIGN_DOUBLE(a_ptmp);
		ASSERT(((intptr_t)a & 63) == 0x0,"a[] not aligned on 64-byte boundary!");
		if(((intptr_t)a & 127) != 0x0)
			fprintf(stderr, "WARN: a[] = 0x%08lX not aligned on 128-byte boundary!\n", (intptr_t)a);
		// v19: Add three more full-residue arrays to support 2-input FFT-modmul needed for Gerbicz check (and later, p-1 support):
		if(use_lowmem < 2) {
			b = a + nalloc;	c = b + nalloc;	d = c + nalloc, e = d + nalloc;
			b_uint64_ptr = (uint64*)b; c_uint64_ptr = (uint64*)c; d_uint64_ptr = (uint64*)d; e_uint64_ptr = (uint64*)e;
		}

		// For this residue (and scratch) byte-array, conservatively figure at least 4 bits per float-double residue word.
		// For multi-FFT-length self-tests, conservatively figure as many as 20 bits (2.5 bytes) per float-double residue word:
		// v20: for largest currently supported FFT of 512Mdoubles, i still -barely - fits in a uint32, but 2.5*i does not:
		arrtmp_alloc = i; arrtmp_alloc = MAX((p+63)>>2, (uint64)(arrtmp_alloc*2.5)) >> 3;	// #limb needed to store p bits = (p+63)>>6, so alloc at least 2x this
		arrtmp = ALLOC_UINT64(arrtmp, arrtmp_alloc);if(!arrtmp ){ sprintf(cbuf, "ERROR: unable to allocate array ARRTMP with %u bytes in main.\n",i); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

		// For an n-word main-array, BIGWORD_BITMAP and BIGWORD_NBITS have (n/64) elts each, thus need 1/64 + 1/32 the total
		// storage of the main-array. Use uint64 alloc-macro for both, so halve the num-elts arg for the BIGWORD_NBITS alloc.
		// As with above arrays, for multi-length self-test, alloc based on max. FFT length used (i) rather than current length (n).
		// Don't need any array padding on these bitmap arrays, but since nalloc includes padding, no harm in using it:
		BIGWORD_BITMAP =           ALLOC_UINT64(BIGWORD_BITMAP, nalloc>>6);	if(!BIGWORD_BITMAP){ sprintf(cbuf, "ERROR: unable to allocate array BIGWORD_BITMAP in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		BIGWORD_NBITS  = (uint32 *)ALLOC_UINT64(BIGWORD_NBITS , nalloc>>7);	if(!BIGWORD_NBITS ){ sprintf(cbuf, "ERROR: unable to allocate array BIGWORD_NBITS in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
	}
// Multithreaded-code debug: Set address to watch:
#ifdef MULTITHREAD
	ADDR0 = a;
#endif

	/* Make sure we start with primary restart file: */
	RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
	strcpy(cstr, RESTARTFILE);
	// G-check applies to primality/Fermat and PRP/Mersenne (but just to this full-PRP phase, not to any ensuing PRP-CF step):
	DO_GCHECK = ( (TEST_TYPE == TEST_TYPE_PRIMALITY) && (MODULUS_TYPE == MODULUS_TYPE_FERMAT) && (use_lowmem < 2) )
				|| ( (TEST_TYPE == TEST_TYPE_PRP) );
	// v19: If PRP test, make sure Gerbicz-checkproduct interval divides checkpoint-writing one:
	if(DO_GCHECK) {
		i = ITERS_BETWEEN_GCHECKS;
		j = ITERS_BETWEEN_GCHECK_UPDATES;
		k = ITERS_BETWEEN_CHECKPOINTS;
		ASSERT(i == j*j, "#iterations between Gerbicz-checksum updates must = sqrt(#iterations between residue-integrity checks)");
		ASSERT(i%k == 0 && k%j == 0, "G-checkproduct update interval must divide savefile-update one, which must divide the G-check interval");
	}

	// PRP-test: Init bitwise multiply-by-base array - cf. comment re. modified Fermat-PRP needed by Gerbicz check
	// above ==> all bits = 0 for Mersenne-PRP-test, rather than all-ones-with-least-significant-bit-0 as for the
	// unmodified Fermat-PRP test of an M(p). Thus no need to call mi64_brev to put bitmap in BRed form used by LR modpow:
	if(TEST_TYPE == TEST_TYPE_PRP) {	// No need to replace with if(DO_GCHECK) here, since Pepin test is all-squarings
		j = (ITERS_BETWEEN_CHECKPOINTS+63) >> 6;
		mi64_clear(BASE_MULTIPLIER_BITS,j);
	}

READ_RESTART_FILE:

	if(!INTERACT)	// Only do the restart-file stuff if it's not a self-test
	{	// 27 Nov 2021: If hit successive G-check errors, on try 2 cstr already has the .G extension, end up with doubled .G.G and
		if(ierr == ERR_GERBICZ_CHECK) {	// "file not found" assertion-exit. So add a re-init of cstr == RESTARTFILE before strcat()
			strcpy(cstr, RESTARTFILE); strcat(cstr, ".G");
		} else if(s2_continuation) {
			strcpy(cstr, RESTARTFILE); strcat(cstr, ".s1");
		}
		/* See if there's a restart file: */
		fp = mlucas_fopen(cstr, "rb");
		/* If so, read the savefile: */
		if(fp) {
			if(TEST_TYPE == TEST_TYPE_PRP) {
				dum = PRP_BASE;
				ASSERT(use_lowmem < 2, "PRP-test mode not available in Low-memory[2] run mode!");
			}
			i = read_ppm1_savefiles(cstr, p, &j, fp, &itmp64,
												(uint8*)arrtmp      , &Res64,&Res35m1,&Res36m1,	// Primality-test residue
												(uint8*)e_uint64_ptr, &i1   ,&i2     ,&i3     );// v19: G-check residue
			fclose(fp); fp = 0x0;
			ilo = itmp64;	// v20: E.g. distributed deep p-1 S2 may use B2 >= 2^32, so made nsquares field in savefiles r/w a uint64
			if(!i) {
				/* First print any error message that may have been issued during the above function call: */
				if(strstr(cbuf, "read_ppm1_savefiles"))
					mlucas_fprint(cbuf,1);
				/* And now for the official spokesmessage: */
				snprintf(cbuf,STR_MAX_LEN*2, "ERROR: read_ppm1_savefiles Failed on savefile %s!\n",cstr);
				mlucas_fprint(cbuf,1);

				if(ierr == ERR_GERBICZ_CHECK) {
					sprintf(cbuf,"Failed to correctly read last-good-Gerbicz-check data savefile!");
					mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
				} else if(cstr[0] != 'q') {
					cstr[0] = 'q';	goto READ_RESTART_FILE;
				} else {
					sprintf(cbuf,"Failed to correctly read both primary or secondary savefile!");
					mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
				}
			}
			// If user attempts to restart run with different PRP base than it was started with, ignore the new value and continue with the initial one:
			if(TEST_TYPE == TEST_TYPE_PRP && dum != PRP_BASE) {
				fprintf(stderr,"INFO: User-specified PRP-test base %u differs from value of %u read from savefile %s ... using the latter.\n",dum,PRP_BASE,cstr);
			}
			// If FFT-length-in-K field returns nonzero and is greater than kblocks (e.g. if ROEs caused switch to larger FFT length),
			// ***and user did not override via cmd-line*** (i.e. fft_length arg not set) reset FFT length and radix set accordingly.
			// v20: If FFTlen was upped due to ROE but err-freq sufficiently low, allow reversion of FFTlen, by skipping the 'if' body:
			//*** Oct 2021: to work safely, this needs to trigger reinit of the default-FFT-length-and-radix-set carry-thread data ***
		  #if USE_FFTLEN_REVERSION
			/* Clause #3: Use larger FFT length in savefile if ROE-frequecy is insufficiently low, defined as
			ROE-frequecy > [one larger-FFT-triggering ROE every 64 checkpoint intervals, on average], thus
			[NERR_ROE > ilo/(64*ITERS_BETWEEN_CHECKPOINTS)], or sans divide, (NERR_ROE<<6)*ITERS_BETWEEN_CHECKPOINTS > ilo]:
			*/
			if((j && j > kblocks) && !fft_length && (NERR_ROE<<6)*ITERS_BETWEEN_CHECKPOINTS > ilo) {
		  #else
			#warning INFO: Disallowing FFT-length reversion for borderline-ROE cases.
			if((j && j > kblocks) && !fft_length) { //^^^^^^^^^^^^^^^^^^^^^^
		  #endif
				kblocks = j;
				// Clear out current FFT-radix data, since get_preferred_fft_radix() expects that:
				for(i = 0; i < NRADICES; i++) { RADIX_VEC[i] = 0; }
				NRADICES = 0;
				goto SETUP_FFT;
			}
			/* On gcheck-error restart, if p near max for the given FFT length, there is a small chance the GEC-failure was caused
			by ROE aliasing, e.g. a floating-point convolution output X.4375 or Y.5625 which gets rounded to X or Y but where the
			exact integer output is really X+1 or Y-1, respectively. In such an case, restart from the last GEC checkpoint data with
			the same circular residue shift will simply reproduce the failure, i.e. report a ROE = |NINT(output) - output| = 0.4375.
			To avoid this infinite-loop-of-GEC-failures scenario, pseudorandomize the FFT inputs by fiddling the residue shift in such
			cases - easiest way is to simple double (mod p) the shift count, then do as normal, i.e. read the shift-removed packed-bit
			residue in last-good savefile and circular-shift it by the modified shift count before calling the convert_res_bytewise_FP()
			function to convert it to balanced-digit floating-point form. There is still a chance that the thus-fiddled FFT inputs will
			again lead to a fatal ROE in the resulting convolution outputs, but the odds of such an ROE again being of the "silent but
			deadly" aliased-ROE type are negligibly small. Even should such an improbability occur, if it does the program will once
			more pseudorandomize the FFT inputs by again mod-doubling the shift count, i.e. we'll never get stuck:
			*/
			if(ierr == ERR_GERBICZ_CHECK) {
				MOD_ADD64(RES_SHIFT,RES_SHIFT,p,RES_SHIFT);
				snprintf(cbuf,STR_MAX_LEN*2, "Gerbicz-check-error restart: Mod-doubling residue shift to avoid repeating any possible fractional-error aliasing in retry, new shift = %llu\n",RES_SHIFT);
				mlucas_fprint(cbuf,1);
			}
			/* Allocate floating-point residue array and convert savefile bytewise residue to floating-point form, after
			first applying required circular shift read into the global RES_SHIFT during the above bytewise-savefile read.
			*/
			if(!convert_res_bytewise_FP((uint8*)arrtmp, a, n, p)) {
				snprintf(cbuf,STR_MAX_LEN*2, "ERROR: convert_res_bytewise_FP Failed on primality-test residue read from savefile %s!\n",cstr);
				mlucas_fprint(cbuf,0);
				if(cstr[0] != 'q' && !(ierr == ERR_GERBICZ_CHECK)) {	// Secondary savefile only exists for regular checkpoint files
					cstr[0] = 'q';
					goto READ_RESTART_FILE;
				} else {
					ASSERT(0,cbuf);
				}
			}
			// v19: G-check residue - we only create savefile for PRP-phase of any PRP-CF run, i.e. always expect a G-check residue:
		  if(DO_GCHECK) {
			if(!convert_res_bytewise_FP((uint8*)e_uint64_ptr, b, n, p)) {
				snprintf(cbuf,STR_MAX_LEN*2, "ERROR: convert_res_bytewise_FP Failed on Gerbicz-check residue read from savefile %s!\n",cstr);
				mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
			} else {
				ierr = 0;
				s1 = sum64(b_uint64_ptr, n); s2 = s3 = s1;	// Init triply-redundant checksum of G-checkproduct
			}
		  }
			ASSERT(ilo > 0,"Require ilo > 0!");
			ihi = ilo+ITERS_BETWEEN_CHECKPOINTS;
			/* If for some reason last checkpoint was at a non-multiple of ITERS_BETWEEN_CHECKPOINTS, round down: */
			ihi-= ihi%ITERS_BETWEEN_CHECKPOINTS;
			// Will check below: if ihi > maxiter, will reset == maxiter
			restart = TRUE;
			/*** already closed restart file ***/
		}
		else /* if(!fp) */
		{
			/* If we're on the primary restart file, set up for secondary: */
			if(ierr == ERR_GERBICZ_CHECK || s2_continuation) {	// Secondary savefile only exists for regular checkpoint files
				snprintf(cbuf,STR_MAX_LEN*2, "INFO: Needed restart file %s not found...moving on to next assignment in %s.\n",cstr,WORKFILE);
				mlucas_fprint(cbuf,1);
				goto GET_NEXT_ASSIGNMENT;
			} else if(cstr[0] != 'q') {
				snprintf(cbuf,STR_MAX_LEN*2, "INFO: primary restart file %s not found...looking for secondary...\n",cstr);
				mlucas_fprint(cbuf,1);
				cstr[0] = 'q';
				goto READ_RESTART_FILE;
			} else {
				sprintf(cbuf, "INFO: no restart file found...starting run from scratch.\n");
				mlucas_fprint(cbuf,1);
				if(ierr == ERR_GERBICZ_CHECK) {
					ierr = 0; restart = FALSE;
				}
			}
		}	/* endif(fp) */
	}	/* endif(!INTERACT)	*/

	if(!restart) {
		/*...set initial iteration loop parameters...	*/
		ilo = 0;
		if(INTERACT)
			ihi = timing_test_iters;
		else
			ihi = ITERS_BETWEEN_CHECKPOINTS;
	}

	ASSERT(MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	/* If at the start of a p-1 or primality test, set the initial seed for the run: */
	ASSERT(TEST_TYPE <= TEST_TYPE_MAX,"Given TEST_TYPE not supported!");
	if(ilo == 0)
	{
		memset(a, 0, npad*sizeof(double));
		if(b) memset(b, 0, npad*sizeof(double));
		if(c) memset(c, 0, npad*sizeof(double));
		if(d) memset(d, 0, npad*sizeof(double));

		/* Always use 3 as the p-1 and Pepin-test seed, and 4 for the LL-test seed. For PRP-test, use seed set in worktodo assignment line: */
		if(TEST_TYPE == TEST_TYPE_PM1) {
			iseed = PRP_BASE;
		} else if(TEST_TYPE == TEST_TYPE_PRP) {	// v21: Enable G-check also for Fermat Pepin test, but use separate clause below due to the PRP_BASE-used-for-something-else issue
			iseed = b[0] = d[0] = PRP_BASE;	// Init the Gerbicz residue-product accumulator b[] and its redundant copy d[].
											// On restart b[] inited via full-bytewise-array read from savefile.
			s1 = s2 = s3 = PRP_BASE;	// Init triply-redundant checksum of G-checkproduct
		} else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
			iseed = 4;
		else	// Pepin test: If using G-check, init these = 3 rather than = PRP_BASE, since for Fermat moduli, PRP_BASE
			iseed = b[0] = d[0] = s1 = s2 = s3 = 3;	// instead is used to implement the random-residue-shift offset scheme

		// In theory could allow residue-shift during P-1, at least in stage 1, but in practice need the BASE_MULTIPLIER_BITS array
		// to hold the part of the stage 1 prime-powers product needed for the current iteration interval of the stage 1 powering:
		if(TEST_TYPE == TEST_TYPE_PM1) {
			ASSERT(RES_SHIFT == 0ull, "Shifted residues unsupported for p-1!\n");
			RES_SHIFT = 0ull; a[0] = iseed;
		} else {
			// Apply initial-residue shift - if user has not set one via cmd-line or current value >= p, randomly choose a value in [0,p).
			// [Note that the RNG is inited as part of the standard program-start-sequence, via function host_init().]
			if(RES_SHIFT == -1ull || RES_SHIFT >= p) {	// Since 0 a valid value of RES_SHIFT, treat -1 as "uninited"
				itmp64 = 0ull;
				while(!itmp64) {	// Eliminate the small probability of the roulette wheel coming up 0
					itmp64 = rng_isaac_rand() % p;
				}
				RES_SHIFT = itmp64;
			} else if(RES_SHIFT == 0) {	// This is needed for multi-exponent runs in which the first assignment is a partially-complete v17 run:
				itmp64 = parse_cmd_args_get_shift_value();
				if(itmp64 == -1ull)	// -1 indicates "-shift not specified on command line", in which case we init-to-random since this is a v18 run-from-scratch
					RES_SHIFT = rng_isaac_rand() % p;
				else	// Otherwise set shift to whatever value was specified on command line
					RES_SHIFT = itmp64 % p;
			}
			// Since residue is otherwise 0, use shifted-carryin function on double-precision padded-array residue:
			itmp64 = shift_word(a, n, p, RES_SHIFT, (double)iseed);	// Note return value (specifically high 7 bytes thereof) is an unpadded index
			ASSERT((itmp64 >>  8) < n                , "Return value of shift_word(): unpadded-array-index out of range!");
			ASSERT((itmp64 & 255) < ceil((double)p/n), "Return value of shift_word(): bit-in-array-word value out of range!");
		}
	} else if(DO_GCHECK) {
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT && TEST_TYPE == TEST_TYPE_PRIMALITY && !INTERACT) {	// Allow shift in timing-test mode
			ASSERT(RES_SHIFT == 0ull, "Shifted residues unsupported for Pépin test with Gerbicz check!\n");
		}
		memcpy(d, b, nbytes);	// If doing a PRP test, init redundant copy d[] Gerbicz residue-product accumulator b[].
	}

	if(restart) {
		if (MODULUS_TYPE == MODULUS_TYPE_FERMAT) snprintf(cbuf,STR_MAX_LEN*2, "Restarting %s at iteration = %u, residue shift count = %llu.\nRes64,Res35m1,Res36m1: %016llX,%llu,%llu\n",PSTRING,ilo,RES_SHIFT,Res64,Res35m1,Res36m1);
		else snprintf(cbuf,STR_MAX_LEN*2, "Restarting %s at iteration = %u. Res64: %016llX, residue shift count = %llu\n",PSTRING,ilo,Res64,RES_SHIFT);
		mlucas_fprint(cbuf,0);
	}

	/*...Restart and FFT info.	*/
	snprintf(cbuf,STR_MAX_LEN*2,"%s: using FFT length %uK = %u 8-byte floats, initial residue shift count = %llu\n",PSTRING,kblocks,n,RES_SHIFT);
	sprintf(cstr,"This gives an average %20.15f bits per digit\n",1.0*p/n);	strcat(cbuf,cstr);
	if(TEST_TYPE == TEST_TYPE_PRP) {
		sprintf(cstr,"The test will be done in form of a %u-PRP test.\n",PRP_BASE);	strcat(cbuf,cstr);
	}
	// If self-test (INTERACT = True), echo only to stderr (uint32 flag > 1);
	// otherwise echo only to logfile (flag = 0); use -INTERACT to set flag appropriately
	mlucas_fprint(cbuf,-INTERACT);

/*...main loop...	*/
/******************* AVX debug stuff: *******************/
#if 0
	int jj,ipad;
	// Use RNG to populate data array:
	rng_isaac_init(TRUE);
	double pow2_dmult = 1024.0*128.0;	// Restrict inputs to 18 bits, which in balanced-digit representation
										// means restricting multiplier of random-inputs-in-[-1,+1] below to 2^17
	for(jj = 0; jj < n; jj += 16) {
		ipad = jj + ( (jj >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		// All the inits are w.r.to an un-SIMD-rearranged ...,re,im,re,im,... pattern:
	#ifdef USE_AVX512
		a[ipad+br16[ 0]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
		a[ipad+br16[ 1]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
		a[ipad+br16[ 2]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
		a[ipad+br16[ 3]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
		a[ipad+br16[ 4]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
		a[ipad+br16[ 5]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
		a[ipad+br16[ 6]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
		a[ipad+br16[ 7]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
		a[ipad+br16[ 8]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
		a[ipad+br16[ 9]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
		a[ipad+br16[10]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
		a[ipad+br16[11]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
		a[ipad+br16[12]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
		a[ipad+br16[13]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
		a[ipad+br16[14]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
		a[ipad+br16[15]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
	#elif defined(USE_AVX)
		a[ipad+br8[0]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
		a[ipad+br8[1]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
		a[ipad+br8[2]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
		a[ipad+br8[3]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
		a[ipad+br8[4]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
		a[ipad+br8[5]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
		a[ipad+br8[6]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
		a[ipad+br8[7]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
		a[ipad+br8[0]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
		a[ipad+br8[1]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
		a[ipad+br8[2]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
		a[ipad+br8[3]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
		a[ipad+br8[4]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
		a[ipad+br8[5]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
		a[ipad+br8[6]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
		a[ipad+br8[7]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
	#elif defined(USE_SSE2)
		a[ipad+br4[0]   ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
		a[ipad+br4[1]   ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
		a[ipad+br4[2]   ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
		a[ipad+br4[3]   ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
		a[ipad+br4[0]+4 ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
		a[ipad+br4[1]+4 ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
		a[ipad+br4[2]+4 ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
		a[ipad+br4[3]+4 ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
		a[ipad+br4[0]+8 ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
		a[ipad+br4[1]+8 ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
		a[ipad+br4[2]+8 ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
		a[ipad+br4[3]+8 ]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
		a[ipad+br4[0]+12]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
		a[ipad+br4[1]+12]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
		a[ipad+br4[2]+12]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
		a[ipad+br4[3]+12]= DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
	#else
		a[ipad+ 0]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
		a[ipad+ 1]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
		a[ipad+ 2]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
		a[ipad+ 3]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
		a[ipad+ 4]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
		a[ipad+ 5]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
		a[ipad+ 6]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
		a[ipad+ 7]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
		a[ipad+ 8]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
		a[ipad+ 9]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
		a[ipad+10]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
		a[ipad+11]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
		a[ipad+12]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
		a[ipad+13]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
		a[ipad+14]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
		a[ipad+15]       = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
	#endif
	}
#endif
/********************************************************/

	// Set function pointer to point to [mers|fermat]_mod_square based on modulus type:
	int	(*func_mod_square)(double [], int [], int, int, int, uint64, uint64, int, double *, int, double *)
						= ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? mers_mod_square : fermat_mod_square);
	int update_shift = (RES_SHIFT != 0ull);	// If shift = 0 at outset, don't update (only need for Fermat-mod, due to the random-bit aspect there)

	if(TEST_TYPE == TEST_TYPE_PM1 && ilo >= maxiter) {
		ASSERT(ilo == maxiter && ilo == PM1_S1_PROD_BITS,"For completed S1 expect ilo == maxiter == PM1_S1_PROD_BITS!");
		snprintf(cbuf,STR_MAX_LEN*2, "%s: p-1 stage 1 to b1 = %u already done -- proceeding to stage 2.\n",PSTRING,B1);
		fprintf(stderr,"%s",cbuf);
		ilo = ihi;		// Need this to differentiate between just-completed S1 and S1 residue read from restart file,
		goto PM1_STAGE2;// in terms of whether we need to do a GCD before proceeding to S2
	} else if(KNOWN_FACTORS[0] != 0ull) {	// PRP-CF - but if ilo < (p-1) it's in the PRP-phase, handle like regular PRP run until that completes
		ASSERT(TEST_TYPE == TEST_TYPE_PRP,"One or more known-factors in workfile entry requires a PRP= assignment type!");
		if( ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) && (ilo >= p))
		 || ((MODULUS_TYPE == MODULUS_TYPE_FERMAT) && (ilo >= p-1)) )
			goto PM1_STAGE2;	// The CF-handling is a clause of the if/else beginning at this label
	}

	for(;;)
	{
		ASSERT(maxiter > 0,"Require (uint32)maxiter > 0");
		if(ihi > maxiter)
			ihi = maxiter;
		// If p-1: start of each iteration cycle, copy bits ilo:ihi-1 of PM1_S1_PRODUCT into low bits of BASE_MULTIPLIER_BITS vector:
		if(TEST_TYPE == TEST_TYPE_PM1) {
			k = (ihi+63)>>6;// #limbs of PM1_S1_PRODUCT needed, INCL. BITS 0:ILO-1 WHICH WILL BE GETTING OFF-SHIFTED - only extract ilo:ihi-1
			j = ITERS_BETWEEN_CHECKPOINTS;	// Need full iter-interval's worth of bits here because for interval ilo:ihi-1, cy-routines
			// read bits starting with bit (ilo % ITERS_BETWEEN_CHECKPOINTS). If user stops run, rebuilds with higher value of the latter
			// and restarts, cy-routine will read bits (ilo % ITERS_BETWEEN_CHECKPOINTS):(ihi % ITERS_BETWEEN_CHECKPOINTS - 1) of the
			// BASE_MULTIPLIER_BITS vector. That further means my original mi64_shrl(PM1_S1_PRODUCT,BASE_MULTIPLIER_BITS, ilo, k,i) call
			// below needs ilo to be replaced by ilo - (ilo % ITERS_BETWEEN_CHECKPOINTS), which simply == ilo in the absence of such
			// within-run ITERS_BETWEEN_CHECKPOINTS-fiddling. 1st-release of v20 had this bug, which hosed my 1st run of F33 s1.
			i = (j+63)>>6; j &= 63;			// i = #limbs needed to hold current bit ilo:ihi-1 window; j = #low bits set in high uint64 of same
			// Copy the needed limbs from arrtmp into BASE_MULTIPLIER_BITS...
			mi64_shrl(PM1_S1_PRODUCT,BASE_MULTIPLIER_BITS, ilo - (ilo % ITERS_BETWEEN_CHECKPOINTS), k,i);
			itmp64 = ~(-1ull << j); BASE_MULTIPLIER_BITS[i-1] &= itmp64;// ...and zero any excess bits at the high end.
			for(i = 0, itmp64 = 0ull; i < s1p_alloc; i++) { itmp64 += PM1_S1_PRODUCT[i]; }
			if(itmp64 != PM1_S1_PROD_RES64) {
				snprintf(cbuf,STR_MAX_LEN*2,"PM1_S1_PRODUCT (mod 2^64_ checksum mismatch! (Current[%llu] != Reference[%llu]). Aborting due to suspected data corruption.\n",itmp64,PM1_S1_PROD_RES64);
				mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
			}
		}
		/* Here's the big one - (ITERS_BETWEEN_CHECKPOINTS) squaring steps.
		XYZ_mod_square returns 0 if no errors detected during this iteration cycle.
		If fatal error was encountered, skip to next assignment in worktodo file
		but keep restart files for current assignment around (can finish using hiacc code.)
		*/
		AME = MME = 0.0;	/* Init Avg. & Max. RO Error */
		AME_ITER_START = 30;/* Start collecting AME after allowing residue to "fill up" in initial few tens of iters */
	#ifdef USE_FGT61
		ierr = func_mod_square  (a,c, (int*)arrtmp, n, ilo, ihi, 0ull, p, scrnFlag, &tdiff, update_shift);
	#else
	  // v19 Gerbicz-checkproduct updates: break iteration loop into subsegments of length ITERS_BETWEEN_GCHECK_UPDATES.
	  if(DO_GCHECK) {
		// Only do check-related stuff for full production-run iteration intervals:
		if(MLUCAS_KEEP_RUNNING && (ihi-ilo) >= ITERS_BETWEEN_GCHECK_UPDATES) {
			i = ilo;	tdiff = 0.0;	// Need 2 timers here - tdif2 for the individual func_mod_square calls, accumulate in tdiff
			while(!ierr && MLUCAS_KEEP_RUNNING && i < ihi) {
				// See G-check code for why this logfile-print of initial-G-check-update residue shift value is needed in Fermat-mod case:
				if(i == ITERS_BETWEEN_GCHECK_UPDATES) { sprintf(cbuf,"At iter ITERS_BETWEEN_GCHECK_UPDATES = %u: RES_SHIFT = %llu\n",i,RES_SHIFT); mlucas_fprint(cbuf,1); }
				/* If restart-after-interrupt and thus ilo neither a non-multiple of ITERS_BETWEEN_CHECKPOINTS nor of
				ITERS_BETWEEN_GCHECK_UPDATES, round first i-update > ilo to nearest multiple of ITERS_BETWEEN_GCHECK_UPDATES:
				*/
				uint32 itodo = ITERS_BETWEEN_GCHECK_UPDATES - i%ITERS_BETWEEN_GCHECK_UPDATES;
				// And this handles case of final odd-length iteration interval:
				if(i+itodo > maxiter)
					itodo = maxiter-i;
				/* iteration 0: Both [a],[b] in pure-integer form; do fwt-weight and initial-fwd-FFT-pass of [a], followed by
				ITERS_BETWEEN_GCHECK_UPDATES FFT mod-squarings.
				mode_flag = 0 for all this since we are FFT-auto-squaring, not FFT-based 2-input modmul, thus fwd_fft_only = 0:
				*/
				first_sub = (i == ilo); last_sub = (i+itodo == ihi);
			// First subinterval: [a] needs fwd-weighting and initial-fwd-FFT-pass done on entry, !undone on exit: mode_flag = 10_2
			// Last  subinterval: [a] !need fwd-weighting and initial-fwd-FFT-pass done on entry but done on exit: mode_flag = 01_2
			// Intermediate subs: [a] !need fwd-weighting and initial-fwd-FFT-pass done on entry, !undone on exit: mode_flag = 11_2
				mode_flag = 3 - first_sub - (last_sub<<1);
				ierr = func_mod_square  (a, (int*)arrtmp, n, i,i+itodo, (uint64)mode_flag, p, scrnFlag, &tdif2, update_shift, 0x0);	tdiff += tdif2;
				if(ierr) {
					fprintf(stderr,"At iteration %d: mod_square returned with error code[%u] = %s\n",ROE_ITER,ierr,returnMlucasErrCode(ierr));
					/* If interrupt *and* we're past the first subinterval, need to undo initial-fwd-FFT-pass and DWT-weighting on b[],
					whose value will reflect the last multiple-of-ITERS_BETWEEN_GCHECK_UPDATES iteration - prior to writing it,
					along with the current PRP residue, to savefile: */
					if(ierr == ERR_INTERRUPT && !first_sub)
						ierr = func_mod_square  (b, (int*)arrtmp, n, i,i+1, 8ull, p, scrnFlag, &tdif2, FALSE, 0x0);
					break;
				}
				/* At end of each subinterval, do a single modmul of current residue a[] with Gerbicz-checkproduct to update the latter:
				Make a copy of a[], e.g. c[] = a[], undo the fwd-FFT-radix pass via a call to radix*_dit_pass1(), then fwd-transform that, then
				compute b *= c (mod n). But again wasteful, would rather just complete the fwd-FFT of c[] instead of undoing and then immediately
				redoing the initial fwd-FFT-radix pass of it. But, we do this infrequently enough that we don't care!
				[EWM: Update: Added the FFT mode_flag control mechanism to eliminate this redundant work] */
				memcpy(c, a, nbytes);	// Copy a into c and do fwd-FFT-only of c:
				// If going to proceed to actual Gerbicz-check, need to save an un-updated copy of the checkproduct:
				i += itodo;
			//	fprintf(stderr,"At Iter %u: copy a[] -> c[]\n",i);
				if(i % ITERS_BETWEEN_GCHECKS == 0) {
					memcpy(d, b, nbytes);
				}
			// All-but-Last  sub: [c] !need fwd-weighting and initial-fwd-FFT-pass done on entry, exit moot since fwd-FFT-only: mode_flag = 01_2:
			// Last  subinterval: [c] needs fwd-weighting and initial-fwd-FFT-pass done on entry, exit moot since fwd-FFT-only: mode_flag = 00_2:
			/* Note: Interrupt during this step should not be a problem, except in the sense that the resulting MLUCAS_KEEP_RUNNING = False
			would prevent the ensuing func_mod_square call to compute FFT(b)*FFT(c) from doing so ... so instead break out of loop: */
				mode_flag = 1 - last_sub;
				ierr = func_mod_square  (c, (int*)arrtmp, n, i,i+1,      4ull + (uint64)mode_flag, p, scrnFlag, &tdif2, FALSE, 0x0);
				if(ierr) {
					if(ierr == ERR_INTERRUPT) {
						fprintf(stderr,"Caught interrupt in fFFT(c) step.\n");
						break;
					} else {
						snprintf(cbuf,STR_MAX_LEN*2,"Unhandled Error of type[%u] = %s in fFFT(c) step - please send e-mail to ewmayer@aol.com with copy of the p*.stat file attached. Proceeding to next assignment...\n",ierr,returnMlucasErrCode(ierr));
						mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
					//	goto GET_NEXT_ASSIGNMENT;
					}
				}
			//	fprintf(stderr,"fFFT(c), mode = %u\n",mode_flag);

				// prior to each b[]-update, check integrity of array data:
				if(!mi64_cmp_eq(b_uint64_ptr,d_uint64_ptr,n)) {	// Houston, we have a problem
					s1 = consensus_checksum(s1,s2,s3);
					if(s1 == sum64(b_uint64_ptr, n)) {	/* b-data good; no-op */
					} else if(s1 == sum64(d_uint64_ptr, n)) {	// c-data good, copy back into b
						memcpy(d, b, nbytes);
					} else	// Catastrophic data corruption
						ASSERT(0, "Catastrophic data corruption detected in G-checkproduct integrity validation ... rolling back to last good G-check. ");
				}

			// First subinterval: [b] needs fwd-weighting and initial-fwd-FFT-pass done on entry, !undone on exit: mode_flag = 10_2
			// Intermediate subs: [b] !need fwd-weighting and initial-fwd-FFT-pass done on entry, !undone on exit: mode_flag = 11_2
			// Last  subinterval: [b] !need fwd-weighting and initial-fwd-FFT-pass done on entry,  undone on exit: mode_flag = 01_2
			/* Note: Interrupt during this step should not be a problem, the handling code in func_mod_square will complete the FFT-mul
			step and force the undo-initial-FFT-pass-and-DWT-weighting step, leaving a pure-int G-check residue ready for savefile-writing: */
				mode_flag = 3 - first_sub - (last_sub<<1);
			//	printf("Iter %u: FFT(b)*FFT(c) step.\n",i);
				ierr = func_mod_square  (b, (int*)arrtmp, n, i,i+1, (uint64)c + (uint64)mode_flag, p, scrnFlag, &tdif2, FALSE, 0x0);
				if(ierr) {
					if(ierr == ERR_INTERRUPT) {
						fprintf(stderr,"Caught interrupt in FFT(b)*FFT(c) step.\n");
						break;
					} else {
						snprintf(cbuf,STR_MAX_LEN*2,"Unhandled Error of type[%u] = %s in FFT(b)*FFT(c) step - please send e-mail to ewmayer@aol.com with copy of the p*.stat file attached. Proceeding to next assignment...\n",ierr,returnMlucasErrCode(ierr));
						mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
					//	goto GET_NEXT_ASSIGNMENT;
					}
				}
			//	fprintf(stderr,"FFT(b)*FFT(c), mode = %u\n",mode_flag);
				// If not going to proceed to actual Gerbicz-check, save a post-updated copy of the checkproduct,
				// in order to guard against single-bit or other data corruption in b[] (h/t George Woltman):
				if(i % ITERS_BETWEEN_GCHECKS != 0) {
					memcpy(d, b, nbytes);
					s1 = sum64(b_uint64_ptr, n); s2 = s3 = s1;	// Init triply-redundant checksum of G-checkproduct
				}
				/**************************************************************************************************************
				Here is some simple *nix bc code illustrating the G-check algorithm extended to circularly-shifted
				PRP initial seeds, in which the b-checkproduct is updated every 100 iters (ITERS_BETWEEN_GCHECK_UPDATES = 1000
				in Mlucas) and the d-product is checked every 1000 iters (10^6 iters in Mlucas). Note to run this code
				one needs to convert the C++ //-comments to C-style ones:

					define gerbicz(p,shift) {
						auto i,j, n,x,b,d,s;
						n = 2^p-1;
						i = 0; s = shift; x = 3*2^s % n; b = x;
						print "Initial residue u0 = 3 shifted s = ",s," places, u*2^s (mod n) = ",x,"\n";
						while(i < p) {
							i = i+1;
							s = 2*s % p;
							x = x^2 % n;
							if(i%100 == 0) {
								// Every so often, supplement Gerbicz checksum-update (*= current residue (mod n)) [2] with check [3]:
								if(i%1000 == 0) {
									d = b; j = 0;
									while(j < 100) {
										j = j+1;
										d = d^2 % n;
									}
									print "Iteration ",i,": Doing Gerbicz check, current shift count s = ",s,"\n";
									// d needs initial-shift applied prior to final scalar multiply:
									d = d*2^shift % n;
									d = 3*d % n;	// Update needed for Gerbicz check [3]
									b = b*x % n;	// Update Gerbicz checksum of form [2]
									if(b == d) {
										print "i = ",i,": Gerbicz checksums match.\n";
									} else {
										print "i = ",i,": Gerbicz checksums mismatch!\n";
										if(shift != 0) {
											j = 0;
											while(j < p) {
												j = j+1;
												if((d*2^j % n) == b) {
													print "Checksums match for b == d*2^",j," (mod n), i.e. d == b*2^",p-j," (mod n)\n";
												}
											}
										}
									}
								} else {	// Every L iters, do Gerbicz checksum-update of form [2]:
									b = b*x % n;
									print "Iteration ",i,": Updated Gerbicz checksum, current shift count s = ",s,"\n";
								}
							}
						}
						print "Final residue = ",x,"\n";
					}
				**************************************************************************************************************/
			}	// end while(!ierr && MLUCAS_KEEP_RUNNING && i < ihi)
		} else if(MLUCAS_KEEP_RUNNING) {	// Final partial-length interval skips G-check
			ierr = func_mod_square  (a, (int*)arrtmp, n, ilo,ihi, 0ull, p, scrnFlag, &tdiff, update_shift, 0x0);
		}
	  } else {
			// For straight LL-test there is (at least at this writing) no known analog of the Gerbicz check:
			ierr = func_mod_square  (a, (int*)arrtmp, n, ilo,ihi, 0ull, p, scrnFlag, &tdiff, update_shift, 0x0);
	  }	// endif(DO_GCHECK)
	#endif	// #ifdef USE_FGT61 [strictly experimental]

		// Roundoff-retry scheme is detailed in comments of the above func_mod_square() function:
		if(ierr) {
			// v19: For "nonzero exit carry" errors due to residue data corruption (not uncommon on e.g. cellphones)
			// added handling - if at least one previous savefile write, try restart from that:
			if((ierr == ERR_CARRY) && !INTERACT && ihi > ITERS_BETWEEN_CHECKPOINTS) {
				sprintf(cbuf,"Retrying from most-recent savefile, in case the issue was due to residue-data corruption.\n");
				mlucas_fprint(cbuf,1);
				goto READ_RESTART_FILE;
			}
			// v20: Simplify the logic here - skip previous interval-retry step:
			if((ierr == ERR_ROUNDOFF) && !INTERACT) {
				ASSERT(ROE_ITER > 0, "ERR_ROUNDOFF returned but ROE_ITER <= 0!");
				n = get_nextlarger_fft_length(n);	kblocks = (n >> 10);
				sprintf(cbuf," Switching to next-larger available FFT length %uK and restarting from last checkpoint file.\n",kblocks);
				mlucas_fprint(cbuf,1);
				NERR_ROE++;
				USE_SHORT_CY_CHAIN = 0;
				ROE_ITER = 0;
				ierr = 0;	// v19: Need to explicitly clear ierr flag here, otherwise get oo retry loop in PRP-test mode
				// Clear out current FFT-radix data, since get_preferred_fft_radix() expects that:
				for(i = 0; i < NRADICES; i++) { RADIX_VEC[i] = 0; }
				NRADICES = 0;
				goto SETUP_FFT;
			}
			else if(ierr == ERR_UNKNOWN_FATAL)
				return(EXIT_FAILURE);
			if(ierr != ERR_INTERRUPT)	// Unknown error ... try simply going to next assignment
				goto GET_NEXT_ASSIGNMENT;
		}

	/*...Every (ITERS_BETWEEN_CHECKPOINTS)th iteration, print timings to stdout or STATFILE.
	If it's a regular (i.e. non-timing) test, also write current residue to restart files.
	*/
		/* v18: Moved term-output below convert_res_FP_bytewise call since in the presence of cshifted residues
		need convert_res_FP_bytewise call to remove cshift and compute SH residues: */
		// Zero high uint64s of target arrays, since double-to-int residue conversion is bytewise & may leave >=1 MSBs in high word untouched:
		// Fermat-mod residue formally needs an extra bit, though said bit should == 1
		// only in the highly unlikely case of a prime-Fermat Pepin-test result:
		j = (p+63+(MODULUS_TYPE == MODULUS_TYPE_FERMAT))>>6;	arrtmp[j-1] = 0ull;
		convert_res_FP_bytewise(	a, (uint8*)      arrtmp, n, p, &Res64, &Res35m1, &Res36m1);	// LL/PRP-test/[p-1 stage 1] residue
		// G-check residue...must not touch i1,i2,i3 again until ensuing write_ppm1_savefiles call!
		if(DO_GCHECK) {
			e_uint64_ptr[j-1] = 0ull;
			convert_res_FP_bytewise(b, (uint8*)e_uint64_ptr, n, p, &i1,&i2,&i3);
		}

		// In interactive-timing-test (e.g. self-tests) mode, do immediate-exit-sans-savefile-write on signal:
		if(INTERACT && (ierr == ERR_INTERRUPT))
			exit(0);

		// In non-interactive (production-run) mode, write savefiles and exit gracefully on signal:
		if(ierr == ERR_INTERRUPT) {
			// First print the signal-handler-generated message:
			mlucas_fprint(cbuf,1);
			ihi = ROE_ITER;	// Last-iteration-completed-before-interrupt saved here
		/*** Nov 2021: interrupt-handling still not stable ... runs that have been underway for a day or more refuse to quit. Just clean-exit w/o savefile write for now: ***/
		//	sprintf(cbuf,"Iter = %u: Writing savefiles and exiting.\n",ihi);
			sprintf(cbuf,"Exiting at Iter = %u.\n",ihi); mlucas_fprint(cbuf,1);
			exit(1);
		}

		/*...Done?	*/
		if(!INTERACT) {
			AME /= (ihi - ilo);	// Don't /= ITERS_BETWEEN_CHECKPOINTS here since final interval is a partial one
			/*...get a quick timestamp...	*/
			calendar_time = time(NULL); local_time = localtime(&calendar_time);
			strftime(timebuffer,SIZE,"%Y-%m-%d %H:%M:%S",local_time);
			const char*iter_or_stage[] = {"Iter#","S1 bit"};	// Tag indicates Primality/PRP-test or p-1 S1 iteration
			/*...print [date in hh:mm:ss | p | iter-count-or-stage progress | %-complete | time | per-iter time | Res64 | max ROE | residue-shift] */
			snprintf(cbuf,STR_MAX_LEN*2, "[%s] %s %s = %u [%5.2f%% complete] clocks =%s [%8.4f msec/iter] Res64: %016llX. AvgMaxErr = %10.9f. MaxErr = %10.9f. Residue shift count = %llu.\n"
				, timebuffer, PSTRING, iter_or_stage[TEST_TYPE == TEST_TYPE_PM1], ihi, (float)ihi / (float)maxiter * 100,get_time_str(tdiff)
				, 1000*get_time(tdiff)/(ihi - ilo), Res64, AME, MME, RES_SHIFT);
			mlucas_fprint(cbuf,scrnFlag);
		}

		// Do not save a final residue unless p-1 (if not, still leave penultimate residue file intact).
		// We don't save "final residue" in cofactor-PRP mode, since in the (mod M(p)) case this is for p+1 squarings (G-check needs this),
		// i.e. needs a mod-div-by-base^2 postprocessing step to put in form of the p-1 squarings of the standard Fermat-PRP test:
		// Comment by Catherine Cowie, 2024: for Pepin tests we actually do need a final residue saved, as the worktodo format for a Pepin test does not include possibility of cofactor testing. See https://github.com/primesearch/Mlucas/pull/11 for more information.
		if ((ihi == maxiter) && (INTERACT || (TEST_TYPE != TEST_TYPE_PM1 && !(TEST_TYPE == TEST_TYPE_PRIMALITY && MODULUS_TYPE == MODULUS_TYPE_FERMAT))))
				break;

		/* If Mersenne/PRP or Fermat test, do Gerbicz-check every million squarings. Make sure current run has done at least
		one intermediate G-checkproduct update. How do we know that current iter-interval (ilo,ihi] contains such an update?
		(That is, a multiple of ITERS_BETWEEN_GCHECK_UPDATES): Set j =  ITERS_BETWEEN_GCHECK_UPDATES, test if ihi/j > ilo/j.
		Examples, all using ITERS_BETWEEN_GCHECK_UPDATES = 1000:
		o ilo = 0, ihi = 1000: integer divide gives ilo/1000 = 0, ihi/1000 = 1 difference > 0, thus contains an update
		o ilo = 1000, ihi = 1507: ilo/1000 = 1, ihi/1000 = 1 difference = 0, thus !contains an update
		*/
		i = ITERS_BETWEEN_GCHECKS; j = ITERS_BETWEEN_GCHECK_UPDATES;
		if(MLUCAS_KEEP_RUNNING && DO_GCHECK && (ilo/j < ihi/j) && (ihi % i) == 0) {
			// Un-updated copy of the checkproduct saved in [d]; square that ITERS_BETWEEN_GCHECK_UPDATES times...
			/*
			Mar 2022: User hit assertion-exit below ... had set CheckInterval = 1000 = ITERS_BETWEEN_GCHECK_UPDATES,
			smallest allowable value, in his mlucas.ini ... that means the iter-interval between 999-1000k, leading
			up to first G-check, had both first_sub - last_sub = 1.
			My original version of the mode_flag setting below implicitly assumed first_sub = 0, last_sub = 1.
			*/
			// Last subinterval always TRUE at this point, i.e. [d] needs fwd-weighting and initial-fwd-FFT-pass undone on exit: mode_flag = 0*_2
			// If First subinterval also TRUE, [d] needs fwd-weighting and initial-fwd-FFT-pass done on entry: mode_flag = 00_2.
			/* Note: Interrupt during this step should not be a problem, the handling code in func_mod_square will complete the FFT-mul
			step and force the undo-initial-FFT-pass-and-DWT-weighting step, leaving a pure-int G-check residue ready for savefile-writing: */
			mode_flag = 1 - first_sub;
			ierr = func_mod_square  (d,0x0, n, ihi,ihi+ITERS_BETWEEN_GCHECK_UPDATES, (uint64)mode_flag, p, scrnFlag, &tdiff, FALSE, 0x0);
			if(ierr) {
				if(ierr == ERR_INTERRUPT) {
					fprintf(stderr,"Caught interrupt in Gerbicz-checkproduct mod-squaring update ... skipping G-check and savefile-update and performing immediate-exit.\n");
					exit(1);
				} else {
					snprintf(cbuf,STR_MAX_LEN*2,"Unhandled Error of type[%u] = %s in Gerbicz-checkproduct mod-squaring update - please send e-mail to ewmayer@aol.com with copy of the p*.stat file attached. Proceeding to next assignment...\n",ierr,returnMlucasErrCode(ierr));
					mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
				//	goto GET_NEXT_ASSIGNMENT;
				}
			}
			// ...then multiply squaring result by initial PRP seed and compare result to updated copy of the checkproduct saved in b[].
			// First zero the high uint64s of the target array, since the double-to-int residue conversion is bytewise, i.e. may leave
			// 1 or more MSBs in high word untouched:
			j = (p+63)>>6;	/*** Jun 2021: cf. convert_res_FP_bytewise() for why we don't include the extra Fermat-modulus bit here ***/
			c_uint64_ptr[j-1] = 0ull;
			// [1] Convert b[],d[] to bytewise form, former assumed already in e[] doubles-array, latter into currently-unused c[] doubles-array:
			convert_res_FP_bytewise(d, (uint8*)c_uint64_ptr, n, p, 0x0,0x0,0x0);
			// Only need to compute this for initial interval - after that the needed adjustment-shift remains constant
			if(ihi == ITERS_BETWEEN_GCHECKS && RES_SHIFT) {
				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
					/* d[] needs initial-shift applied prior to final scalar multiply, but don't explicitly,
					store the initial shift, so need to recompute it from a current value s at iteration i:
						s = s0.2^i (mod p)
					Repeatedly div-by-2 (mod p) of y = 2.x (mod p):
						y even: y = y>>1
						y odd : y = (y+p)>>1
					*/
					itmp64 = RES_SHIFT;
					for(i = ITERS_BETWEEN_GCHECK_UPDATES; i < ITERS_BETWEEN_GCHECKS; i++) {	// Recover shift at initial ITERS_BETWEEN_GCHECK_UPDATES-iteration subinterval from that at initial savefile-checkpoint
						if(itmp64 & 1)	// y odd
							itmp64 = (itmp64+p)>>1;
						else			// y even
							itmp64 >>= 1;
					}
				} else {
					// In Fermat-mod case, with its random-bit shift offset, simple repeated-mod-halving does not work,
					// need to also account for the said per-iter offset ... but that's no good either, since only store
					// the random-offset bits for latest ITERS_BETWEEN_CHECKPOINTS iters in BASE_MULTIPLIER_BITS. So instead
					// must write RES_SHIFT value at iter = ITERS_BETWEEN_GCHECKS to logfile and read back here:
				#if 1
					if(filegrep(STATFILE,"ITERS_BETWEEN_GCHECK_UPDATES",cbuf,0)) {
						char_addr = strstr(cbuf,"RES_SHIFT = ") + 12;	// Skip ahead by length of search-substring
						itmp64 = strtoull(char_addr, &cptr, 10);
						ASSERT(itmp64 != -1ull, "strtoull() overflow detected.");
					}
				#else
					itmp64 = RES_SHIFT;
					// Unlike Mers-mod case, need to run this loop in reverse in order to duplicate
					// the actual iteration counts and their corr. random-bit shift offsets:
					for(i = ITERS_BETWEEN_GCHECKS; i >= ITERS_BETWEEN_GCHECK_UPDATES; i--) {	// Recover shift at initial ITERS_BETWEEN_GCHECK_UPDATES-iteration subinterval from that at initial savefile-checkpoint
						uint32 nhalvings,curr_bit = ((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1);	// No mod needed on this add, since result of pvs line even and < p, which is itself even in the Fermat-mod case (p = 2^m)
						// If current random-offset bit = 1, do 2 mod-halvings; otherwise do just one:
						for(nhalvings = 0; nhalvings <= curr_bit; nhalvings++) {
							if(itmp64 & 1)	// y odd
								itmp64 = (itmp64+p)>>1;
							else			// y even
								itmp64 >>= 1;
						}
					}
				#endif
				}
				fprintf(stderr,"Recovered initial shift %llu\n",itmp64);
				ASSERT((itmp64>>32) == 0ull,"Shift must be < 2^32!");
				GCHECK_SHIFT = itmp64;
			}
			mi64_shlc(c_uint64_ptr, c_uint64_ptr, (uint32)p, (uint32)GCHECK_SHIFT, j, (MODULUS_TYPE == MODULUS_TYPE_FERMAT));
			/*** Now that have undone shift, include extra modulus bit for Fermat-mod case ***/
			if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) { c_uint64_ptr[j++] = 0ull; }
			// Use mi64 routines to compute d[]*PRP_BASE and do ensuing equality check:
			itmp64 = ((MODULUS_TYPE == MODULUS_TYPE_FERMAT) ? 3ull : (uint64)PRP_BASE);	// Fermat-mod uses PRP_BASE to store 2 for random-shift-offset scheme
			c_uint64_ptr[j] = mi64_mul_scalar(c_uint64_ptr, itmp64, c_uint64_ptr, j);
			ASSERT(c_uint64_ptr[j] == 0ull, "d[]*PRP_BASE result has unexpected carryout!");
			// Need to (mod N) ... store modulus N in d[] doubles-array, which is freed up by above convert_res_FP_bytewise(d,...) call:
			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
				// Loop rather than call to mi64_set_eq_scalar here, since need to set all elts = -1:
				for(i = 0; i < j; i++) { d_uint64_ptr[i] = -1ull; }
				d_uint64_ptr[j-1] >>= 64-(p&63);	// Leading word needs >> to leave just low p%64 bits set
			} else {
				// j = uint64 vector length; init sans the leading '1' word, then increment prior to mi64_div
				mi64_clear(d_uint64_ptr,j);
				d_uint64_ptr[j-1] = d_uint64_ptr[0] = 1ull;
			}
			for(i = 1; i < itmp64; i++) {
				if(!c_uint64_ptr[j] && mi64_cmpult(c_uint64_ptr,d_uint64_ptr,j)) break;
				cy = mi64_sub(c_uint64_ptr,d_uint64_ptr,c_uint64_ptr,j);	// c -= d, with d = 2^p-1
				c_uint64_ptr[j] -= cy;	//ASSERT(cy == 0ull, "mi64_sub result has unexpected borrow!");
			}
			ASSERT(mi64_cmpult(c_uint64_ptr,d_uint64_ptr,j), "Gerbicz checkproduct reduction (mod 2^p-1) failed!");
			if(mi64_cmp_eq(e_uint64_ptr,c_uint64_ptr,j)) {
				sprintf(cbuf,"At iteration %u, shift = %llu: Gerbicz check passed.\n",ihi,RES_SHIFT);
				mlucas_fprint(cbuf,0);
				// In G-check case we need b[] for that, thus skipped the d = b redundancy-copy ... do that now:
				memcpy(d, b, nbytes);
				s1 = sum64(b_uint64_ptr, n); s2 = s3 = s1;	// Init triply-redundant checksum of G-checkproduct
			} else {
				i = mi64_shlc_bits_align(e_uint64_ptr,c_uint64_ptr,p);
				if(i != -1) {
					sprintf(cbuf,"Gerbicz check passes if D *= 2^%u (mod 2^p-1)\n",i);
					mlucas_fprint(cbuf,0);
					// In G-check case we need b[] for that, thus skipped the d = b redundancy-copy ... do that now:
					memcpy(d, b, nbytes);
					s1 = sum64(b_uint64_ptr, n); s2 = s3 = s1;	// Init triply-redundant checksum of G-checkproduct
				} else {
					if(ihi == ITERS_BETWEEN_GCHECKS)
						sprintf(cbuf,"Gerbicz check iteration %u failed! Restarting from scratch.\n",ihi);
					else
						sprintf(cbuf,"Gerbicz check iteration %u failed! Restarting from last-good-Gerbicz-check data.\n",ihi);
					mlucas_fprint(cbuf,0);
					ierr = ERR_GERBICZ_CHECK;
					NERR_GCHECK++;
					goto READ_RESTART_FILE;
				}
			}
		}

		/* Make sure we start with primary restart file: */
		RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');

		/* Oct 2014: Add every-10-million-iter file-checkpointing: deposit a unique-named restart file
		             p[exponent].xM every 10 million iterations, on top of the usual checkpointing.
		To avoid having to write an extra copy of the p-savefile, wait for the *next* checkpoint -
		i.e. ihi = (x million + ITERS_BETWEEN_CHECKPOINTS), or more simply, ilo = (x million) -
		then simply rename the (not yet updated) p-savefile to add the -M extension, and open a
		new version of the p-savefile on the ensuing checkpointing:
		*/
		if((ilo > 0) && (ilo%10000000 == 0)) {
			sprintf(cbuf, ".%dM", ilo/1000000);
			strcpy(cstr, RESTARTFILE);
			strcat(cstr, cbuf);
			if(rename(RESTARTFILE, cstr)) {
				snprintf(cbuf,STR_MAX_LEN*2,"ERROR: unable to rename %s restart file ==> %s ... skipping every-10M-iteration restart file archiving\n",WORKFILE,cstr);
				fprintf(stderr,"%s",cbuf);
			}
		}	// ilo a multiple of 10 million?

	WRITE_RESTART_FILE:

		itmp64 = ihi;
		// If Pepin test is at final iteration, change PRP base to 3 for final write to file (cf. earlier assignment at line 1214). More info: https://github.com/primesearch/Mlucas/pull/11
		if (ihi == maxiter && TEST_TYPE == TEST_TYPE_PRIMALITY && MODULUS_TYPE == MODULUS_TYPE_FERMAT) PRP_BASE = 3;
		fp = mlucas_fopen(RESTARTFILE, "wb");
		if(fp) {		// In the non-PRP-test case, write_ppm1_savefiles() treats the latter 4 args as null:
			write_ppm1_savefiles(RESTARTFILE,p,n,fp, itmp64, (uint8*)arrtmp,Res64,Res35m1,Res36m1, (uint8*)e_uint64_ptr,i1,i2,i3);
			fclose(fp); fp = 0x0;
			/* If we're on the primary restart file, set up for secondary: */
			if(RESTARTFILE[0] != 'q') {
				RESTARTFILE[0] = 'q';	goto WRITE_RESTART_FILE;
			} else {
				RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
			}
		} else {
			snprintf(cbuf,STR_MAX_LEN*2, "ERROR: unable to open restart file %s for write of checkpoint data.\n",RESTARTFILE);
			mlucas_fprint(cbuf,1);
			/*
			Don't want to assert here - asllow processing to continue, in case this is a transient failure-to-open.
			e.g. due to a backup utility having temporarily locked the savefile, or an out-of-disk-space problem.
			*/
		}

		// v19: For PRP tests, add every-ITERS_BETWEEN_GCHECKS-iter last-good-Gerbicz-check savefile,
		// whose content gets updated that frequently but whose name remains fixed:
		RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
		if(DO_GCHECK) {
			if(ihi%ITERS_BETWEEN_GCHECKS == 0) {
				strcpy(cstr, RESTARTFILE);
				strcat(cstr, ".G");
				fp = mlucas_fopen(cstr, "wb");
				if(fp) {
					write_ppm1_savefiles(cstr,p,n,fp, itmp64, (uint8*)arrtmp,Res64,Res35m1,Res36m1, (uint8*)e_uint64_ptr,i1,i2,i3);
					fclose(fp); fp = 0x0;
				} else {
					snprintf(cbuf,STR_MAX_LEN*2, "ERROR: unable to open Gerbicz-check savefile %s for write of checkpoint data.\n",cstr);
					mlucas_fprint(cbuf,1);
				}
			}	// ihi a multiple of ITERS_BETWEEN_GCHECKS?
		}

		if(ierr == ERR_INTERRUPT) exit(0);

		// For Fermats and cofactor-PRP tests of either modulus type, exit only after writing final-residue checkpoint file:
		if( ihi == maxiter ) break;

		if((TEST_TYPE == TEST_TYPE_PRIMALITY || TEST_TYPE == TEST_TYPE_PRP) && p > (1ull << 32) && ihi >= MAX_SELFTEST_ITERS) {
			fprintf(stderr,"Full-length LL/PRP tests on exponents this large not supported; exiting at self-test limit of %u.\n",MAX_SELFTEST_ITERS);
			exit(0);
		}

	  #if USE_FFTLEN_REVERSION
		/* If FFTlen was pvsly auto-increased due to ROE but latest interval clean and err-freq sufficiently low, revert
		FFTlen to default. 2nd clause in the if() guards against doing this if the user has forced the larger FFT.
		Our criterion for such FFT-length reversion is
			(64 * NERR_ROE)*ITERS_BETWEEN_CHECKPOINTS <= ihi ==> NERR_ROE <= ihi/ITERS_BETWEEN_CHECKPOINTS/64 ,
		or in words "if cumulative count of worrisome ROEs <= 1 every 64th iteration interval, on average, revert".
		*/
		if(kblocks > get_default_fft_length(p) && kblocks > fft_length && (NERR_ROE<<6)*ITERS_BETWEEN_CHECKPOINTS <= ihi) {
			kblocks = get_default_fft_length(p);	// Default FFT length in Kdoubles for this exponent
			if(n > (kblocks << 10))		// We are already running at a larger-than-default FFT length
				n = kblocks << 10;
			USE_SHORT_CY_CHAIN = USE_SHORT_CY_CHAIN_MAX;
			// Clear out current FFT-radix data, since get_preferred_fft_radix() expects that:
			for(i = 0; i < NRADICES; i++) { RADIX_VEC[i] = 0; }
			NRADICES = 0;
			goto SETUP_FFT;
		}
	  #endif
		/*...reset loop parameters and begin next iteration cycle...	*/
		ilo = ihi;
		ihi = ilo+ITERS_BETWEEN_CHECKPOINTS;
		// If Fermat-mod-with-residue-shift, re-init the random-offset-bit array:
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT && update_shift) {
			j = ((ITERS_BETWEEN_CHECKPOINTS+63) >> 6);
			for(i = 0; i < j; i++) { BASE_MULTIPLIER_BITS[i] = rng_isaac_rand(); }; i = ITERS_BETWEEN_CHECKPOINTS&63;	// i = #low bits used in high limb
		//	fprintf(stderr,"Iter %u: random-bit-array popcount = %u\n",ilo,mi64_popcount(BASE_MULTIPLIER_BITS,j) - popcount64(BASE_MULTIPLIER_BITS[j-1] >> i));
		}
	}

	/*...For timing tests, print timing and 64-bit residue and return timing via arglist runtime ptr...	*/
	*runtime = tdiff;

	/* If Selftest mode... */
	if(INTERACT) {
		fprintf(stderr, "%u iterations of %s with FFT length %u = %u K, final residue shift count = %llu\n",timing_test_iters,PSTRING,n,kblocks,RES_SHIFT);
		// If TEST_TYPE non-default (e.g. PRP for Mersennes), add text indicating that:
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRP)
			sprintf(cbuf,"PRP-%u ",PRP_BASE);
		else
			cbuf[0] = '\0';

		/* If Fermat number, make sure exponent a power of 2: */
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
			ASSERT((p >> findex) == 1,"Require (p >> findex) == 1");

		if(timing_test_iters > AME_ITER_START) {
			AME /= (timing_test_iters - AME_ITER_START);
			fprintf(stderr, "%sRes64: %016llX. AvgMaxErr = %10.9f. MaxErr = %10.9f. Program: E%s\n", cbuf, Res64, AME, MME, VERSION);
		} else {
			fprintf(stderr, "%sRes64: %016llX. AvgMaxErr N/A. MaxErr = %10.9f. Program: E%s\n", cbuf, Res64, MME, VERSION);
		}
		/* MSVC/.NET incorrectly output these when using uint64 and %20llu format, so cast to double and print: */
		fprintf(stderr, "Res mod 2^35 - 1 = %20.0f\n",(double)Res35m1);
		fprintf(stderr, "Res mod 2^36 - 1 = %20.0f\n",(double)Res36m1);
		/* If they are provided, check the Selfridge-Hurwitz residues: */
		if(sh0) {
			if(*sh0 != 0) {
				if (Res64 == *sh0)
					resFlag = 0;
				else {
					resFlag = 1;	/* False */
					fprintf(stderr, "  ***   Res64 Error   ***\n");
					fprintf(stderr, " current   = %20.0f\n", (double)Res64);
					fprintf(stderr, " should be = %20.0f\n", (double) *sh0);
				}
			} else
				*sh0 = Res64;
		}
		if(sh1) {
			if(*sh1 != 0) {
				if (Res35m1 == *sh1)
					resFlag = 0;
				else {
					resFlag = 1;	/* False */
					fprintf(stderr, "  ***   Res35m1 Error   ***\n");
					fprintf(stderr, " current   = %20.0f\n", (double)Res35m1);
					fprintf(stderr, " should be = %20.0f\n", (double) *sh1);
				}
			} else
				*sh1 = Res35m1;
		}
		if(sh2) {
			if(*sh2 != 0) {
				if (Res36m1 == *sh2)
					resFlag = 0;
				else {
					resFlag = 1;	/* False */
					fprintf(stderr, "  ***   Res36m1 Error   ***\n");
					fprintf(stderr, " current   = %20.0f\n", (double)Res36m1);
					fprintf(stderr, " should be = %20.0f\n", (double) *sh2);
				}
			} else
				*sh2 = Res36m1;
		}
		/*...print runtime in hh:mm:ss format.	*/
		fprintf(stderr, "Clocks =%s\n",get_time_str(tdiff) );
		/*exit(EXIT_SUCCESS);*/
 		return(resFlag);
	}	/* endif(INTERACT) */

PM1_STAGE2:	// Stage 2 invocation is several hundred lines below, but this needs to go at start of big if/else block containing that.
	/*...Check (probable) primality:	*/
	if(TEST_TYPE == TEST_TYPE_PRIMALITY || TEST_TYPE == TEST_TYPE_PRP)
	{
		/* v19: for Mersenne-mod PRP, we use a Gerbicz-check modified-PRP test with 2 more squarings than standard Fermat-mod PRP,
		at the end of which PRP-ness is indicated by the final residue == PRP_BASE^2. Since PRP_BASE is a uint32, the latter is
		as large as 64 bits, thus simply check whether Res64 == PRP_BASE^2 and infer all-higher-words-0-ness by additionally
		checking whether Res64 % 2^35-1 == Res35m1 and Res64 % 2^36-1 == Res36m1:
		*/
		if(TEST_TYPE == TEST_TYPE_PRP && MODULUS_TYPE != MODULUS_TYPE_FERMAT)	// Applies only to mod-M(p) case,
		{																		// Pepin-test and LL are handled in next clause.
			ASSERT(ihi == p, "Gerbicz-check-modified PRP-test requires p mod-squarings!");
			/* Final PRP-residue which is *reported*, OTOH, is the standard Fermat-style (p-2)-squaring one.
			That requires us to do 2 mod-divs of the 2-squares-too-many prp-residue r by the PRP-test base b.
			If b divides r, we're good. Otherwise, need to find multiple of modulus m = 2^p-1 which needs to
			be added to r to make result divisible by b, i.e. k such that (r + k.m) == 0 (mod b); here's how:

			Compute r' := r (mod b) and m' := m (mod b). Then, since b small (typically O(1)) can simply start
			with r' and loop, each pass adding m' until sum == 0 (mod b), needing at most b-1 loop execs.
			Or better, use mod-inverses to solve
									r' + k.m' == 0 (mod b)	[*]
			Compute mi' such that m'.mi' == 1 (mod b). Then multiplying [*] by mi' we have
									r'.mi' + k.m'.mi' == r'.mi' + k == 0 (mod b),
			which we solve for k:
									k = -r'.mi' (mod b) .
			The best part of the mod-inverse approach is that since b is 32-bit, we can handle both mod-div-by-b steps
			in one go, by instead working the above computations (mod b^2):
				1. Compute r" := r (mod b^2) and m" := m (mod b^2);
				2. Compute mi" such that m".mi" == 1 (mod b^2)
				3. Solve for k = -r".mi" (mod b^2) .
			*/
			j = (p+63)>>6;	itmp64 = PRP_BASE*(uint64)PRP_BASE;
			// Init vector to hold m = 2^p-1:
			d_uint64_ptr[j-1] = (1ull<<(p&63)) - 1ull;	for(i = 0; i < (j-1); i++) { d_uint64_ptr[i] = -1ull; }
		  // [1]:
			mi64_div(arrtmp, &itmp64, j,1, 0x0,&rmodb);	// Omit quotient computation; remainder r" in rmodb
			mmodb = twopmmodq64(p,itmp64);				// m"
			// In the most common case PRP_BASE = 3, use that 2^6 == 1 (mod 9), thus 2^p == 2^(p mod 6) (mod 9)
			if(PRP_BASE == 3)
				ASSERT(mmodb == (1ull<<(p % 6)) % 9,"2^p == 2^(p mod 6) (mod 9) fails!");
			// mmodb = (2^p-1) % base ... for reasons unknown, the macro MOD_SUB64 was not inlined properly under gdb
			if(mmodb)
				mmodb--;
			else
				mmodb = itmp64-1ull;	//MOD_SUB64(mmodb,1ull,itmp64, mmodb);
		  // [2]: mi" in i1:
			i1 = modinv64(mmodb,itmp64);
		  // [3]: k ends up in i2, and may need reduction (mod b^2):
		  #ifdef MUL_LOHI64_SUBROUTINE
			MUL_LOHI64(rmodb,i1,&i2,&i3);
		  #else
			MUL_LOHI64(rmodb,i1, i2, i3);
		  #endif
			i2 %= itmp64;	ASSERT(i3 == 0ull, "K-multiplier needs 64-bit reduction (mod b^2)!");
			if(i2) i2 = itmp64 - i2;	// if(k) k = -r".mi" (mod b^2) = b^2 - r".mi" .
			// i2 contains the needed multiplier k. Since ensuing quotient computation needs separate arrays
			// for dividend and quotient, stash output of mi64_mul_scalar_add_vec2 in c[] and ensuing quotient back in arrtmp[]:
			c_uint64_ptr[j] = mi64_mul_scalar_add_vec2(d_uint64_ptr,i2,arrtmp, c_uint64_ptr, j);
			// Now short-div - allowing for the possibility of a carryout from above mi64_mul_scalar_add_vec2() call -
			// by base and check that remainder 0. Note that we do want the quotient now, as that is our reside/base:
			mi64_div(c_uint64_ptr, &itmp64, j+1,1, arrtmp,&rmodb);	ASSERT(rmodb == 0ull,"After short-div, R != 0 (mod B)");
			// And recompute the S-H residues:
			res_SH(arrtmp,j,&Res64,&Res35m1,&Res36m1);
			// Now that residue is standard Fermat-PRP-test one, check if == 1:
			isprime = (arrtmp[0] == 1ull);
			if(isprime) { // Check the arrtmp[] array for non-zero elements; see https://github.com/primesearch/Mlucas/issues/15
				for(i = 1; i < j; i++) {
					if(arrtmp[i] != 0x0) { isprime = 0; break; }
				}
			}
		} else {	// older impl. of LL-test isprime parsed the entire double-float residue array:
			isprime = 1;
			/* For Fermat numbers, in balanced-digit form, it's prime if the lowest-order digit = -1, all others 0: */
			final_res_offset = (MODULUS_TYPE == MODULUS_TYPE_FERMAT);
			a[0] += final_res_offset;
			for(i = 0; i < n; i++) {
				j = i + ( (i >> DAT_BITS) << PAD_BITS );
				if(a[j] != 0.0) { isprime = 0; break; }
			}
			a[0] -= final_res_offset;
		}

	/************************************************************************************************************************/
	// At this point, for LL-test-of-M(p)/Pepin-test-of-F(m) we have redundant uint64/double residues in arrtmp[] and a[].
	// For PRP-test of M(p), arrtmp[] = uint64 base-3 Fermat-PRP residue gotten via above a[]/9 (mod M(p)) postprocessing step.
	// Thus, for PRP-CF code below, for M(p) we're good to go; for F(m) we need 1 more FFT-square of a[] and convert-to-int,
	// which logic is handled in the Suyama_CF_PRP() routine.
	// At present the Primenet server handles PRP-CF reporting by way of the Res64 for the preceding M(p) PRP-test,
	// the list of known factors used, and the basic "is PRP?" result for the corresponding cofactor.
	/************************************************************************************************************************/

		char Res2048[513];
		// Must save Res2048 before PRP cofactor test: https://github.com/primesearch/Mlucas/issues/25
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			for (int i = 31; i >= 0; i--) sprintf(Res2048+496-i*16, "%016llX", arrtmp[i]);

		// v21: PRP-CF: Cofactor-PRP test applies to primality/Fermat (which we follow by 1 additional mod-squaring
		// to convert the base^((N-1)/2) Pepin/Euler-PRP residue to a base^(N-1) Fermat-PRP one) and PRP/Mersenne residues:
		if(KNOWN_FACTORS[0])	// This is automatically false for LL-test
			isprime = Suyama_CF_PRP(p, &Res64, nfac, a,b,arrtmp, ilo, func_mod_square, n, scrnFlag, &tdiff, gcd_str);

		// JSON-formatted result report for LL/PRP/cofactor-PRP tests of M(p):
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			// UTC is for the result as submitted to the server:
			calendar_time = time(NULL);
			gm_time = gmtime(&calendar_time);
			if(!gm_time)	// If UTC not available for some reason, just substitute the local time:
				gm_time = localtime(&calendar_time);
			// Want 'UTC' instead of 'GMT', so include that in lieu of the %Z format specifier (but also see https://www.mersenneforum.org/showthread.php?p=653672#post653672 as UTC is actually redundant)
			strftime(timebuffer,SIZE,"%Y-%m-%d %H:%M:%S",gm_time);
			// Trio of p-1 fields all 0; cstr holds the formatted output line here. Need to differentiate
			// between this PRP-CF result and the preceding PRP; set the otherwise-unused s2_partial flag:
			generate_JSON_report(isprime,p,n,Res64,Res2048,timebuffer, 0,0ull,0x0,s2_partial, cstr);
		}
		// For Non-cofactor LL/PRP runs, print summary status to logfile:
		if(!KNOWN_FACTORS[0]) {
			// Unbelievable - I must be hallucinating:
			if(isprime) {
				if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
					/*... this gets written both to file and to stdout, the latter irrespective of whether the run is in interactive mode...	*/
					snprintf(cbuf,STR_MAX_LEN*2, "%s is a new FERMAT PRIME!!!\nPlease send e-mail to ewmayer@aol.com.\n",PSTRING);
					mlucas_fprint(cbuf,1);
				}
				else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					for(i=0; knowns[i] != 0; i++) {
						if(p == knowns[i])
							break;
					}
					if(knowns[i] != 0) {
						snprintf(cbuf,STR_MAX_LEN*2, "%s is a known MERSENNE PRIME.\n",PSTRING);
						mlucas_fprint(cbuf,(INTERACT || scrnFlag));	// Latter clause == "Echo output to stderr?"
					} else {
						// This gets written both to file and to stderr, the latter irrespective of whether the run is in interactive mode:
						snprintf(cbuf,STR_MAX_LEN*2, "%s is a (probable) new MERSENNE PRIME!!!\nPlease send e-mail to ewmayer@aol.com and woltman@alum.mit.edu.\n",PSTRING);
						mlucas_fprint(cbuf,1);
					}
				}
				else
					ASSERT(0, "Unsupported modulus type!");
			}
			/*
			The more likely scenario - it's not prime, so we form a 64-bit residue and write that.
			If residue has < 64 bits, print a warning.
			*/
			else {
				// Otherwise, write the 64-bit hex residue. As of v19, we write the old-style HRF-formatted result
				// just to the exponent-specific logfile, and the server-expected JSON-formatted result to the results file:
				// Note that Fermat primality tests are not submitted to server, so accordingly we slightly modify the output. More info: https://github.com/primesearch/Mlucas/pull/11
				snprintf(cbuf,STR_MAX_LEN*2, "%s is not prime. Program: E%s. Final residue shift count = %llu.\n",PSTRING,VERSION,RES_SHIFT);
				mlucas_fprint(cbuf,1);
				if (MODULUS_TYPE == MODULUS_TYPE_FERMAT) snprintf(cbuf,STR_MAX_LEN*2, "Selfridge-Hurwitz residues Res64,Res35m1,Res36m1 = %016llX,%11llu,%11llu.\n",Res64,Res35m1,Res36m1);
				else {
					snprintf(cbuf,STR_MAX_LEN*2, "If using the manual results submission form at mersenne.org, paste the following JSON-formatted results line:\n%s\n",cstr);
					// v19: Finish with the JSON-formatted result line:
					fp = mlucas_fopen(OFILE,"a");
					if(fp) {
						fprintf(fp,"\n%s",cstr); fclose(fp); fp = 0x0;
					}
				}
				mlucas_fprint(cbuf,1);
			}
		} else if (MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {	// Cofactor-PRP run:
			snprintf(cbuf,STR_MAX_LEN*2,"If using the manual results submission form at mersenne.org, paste the following JSON-formatted results line:\n%s\n",cstr);
			mlucas_fprint(cbuf,1);
			// Write JSON-formatted result line to results file:
			fp = mlucas_fopen(OFILE,"a");
			if(fp) {
				fprintf(fp,"\n%s",cstr); fclose(fp); fp = 0x0;
			}
		}
	} else if(TEST_TYPE == TEST_TYPE_PM1) {
		// If just completed S1, do a GCD. (ihi == maxiter) is true of both just-completed S1 and completed-S1 residue read from savefile,
		// but in the latter case set ilo == ihi to differentiate between the two. ***6/22/21: BUT! If run halted mid-GCD, on restart
		// will have ilo == ihi ... supplement with what amounts to 'grep GCD [STATFILE]', if found, then GCD completed:
		if( ilo < ihi || !filegrep(STATFILE,"GCD",cbuf,0))
		{	// j = #limbs; clear high limb before filling arrtmp[0:j-1] with bytewise residue just to be sure:
			j = (p+63+(MODULUS_TYPE == MODULUS_TYPE_FERMAT))>>6; arrtmp[j-1] = 0ull;
			convert_res_FP_bytewise(a,(uint8*)arrtmp,n,p,0x0,0x0,0x0);
			arrtmp[0] -= 1;	// S1 GCD needs residue-1
			i = gcd(1,p,arrtmp,0x0,j,gcd_str);	// 1st arg = stage just completed
			// If factor found, gcd() will have done needed status-file-writes:
			if(i || B2 <= B1) {	// Need to also account for the possibility of no-stage-2, in which case B2 <= B1
				// Write JSON output and go to next assignment:
				B2_start = B2 = 0ull;	// Need JSON output to reflect that no S2 was run
				calendar_time = time(NULL);
				gm_time = gmtime(&calendar_time);
				if(!gm_time)	// If UTC not available for some reason, just substitute the local time:
					gm_time = localtime(&calendar_time);
				strftime(timebuffer,SIZE,"%Y-%m-%d %H:%M:%S UTC",gm_time);
				generate_JSON_report(0,p,n,0ull,NULL,timebuffer, B1,B2,gcd_str,s2_partial, cstr);	// cstr holds JSONified output
				snprintf(cbuf,STR_MAX_LEN*2, "If using the manual results submission form at mersenne.org, paste the following JSON-formatted results line:\n%s\n",cstr);
				mlucas_fprint(cbuf,0);
				fp = mlucas_fopen(OFILE,"a");
				if(fp) {
					fprintf(fp,"\n%s",cstr); fclose(fp); fp = 0x0;
				}
			}
		}
		// Huge-FFT runs on low-mem systems only allow Stage 1:
		if(!use_lowmem) {
			// Stage 2: set #buffers based on available system RAM
			if(B2_start < B2) {
				// Must clear BASE_MULTIPLIER_BITS in prep for S2:
				j = ((ITERS_BETWEEN_CHECKPOINTS+63) >> 6);
				for(i = 0; i < j; i++) {
					BASE_MULTIPLIER_BITS[i] = 0ull;
				}
				// If user did not preset via -pm1_s2_nbuf, set PM1_S2_NBUF based on available RAM.
				// n doubles need n/2^17 MB; use KB and double-math for the intermediates here:
				j = SYSTEM_RAM*(double)MAX_RAM_USE*0.01;
				if(!PM1_S2_NBUF) {
				#ifdef OS_TYPE_LINUX
					printf("INFO: %u MB of free system RAM detected; will use up to %u%% = %u MB of that for p-1 stage 2.\n",SYSTEM_RAM,MAX_RAM_USE,j);
				#else
					printf("INFO: %u MB of available system RAM detected; will use up to %u%% = %u MB of that for p-1 stage 2.\n",SYSTEM_RAM,MAX_RAM_USE,j);
				#endif
					// Assume 5 of those n-double chunks are already used for Stage 1 residue and other stuff
					PM1_S2_NBUF = (uint32)(j*1024./(n>>7)) - 5;
					fprintf(stderr,"Available memory allows up to %u Stage 2 residue-sized memory buffers.\n",PM1_S2_NBUF);
				} else {
					fprintf(stderr,"User specified a maximum of %u Stage 2 residue-sized memory buffers.\n",PM1_S2_NBUF);
					if( PM1_S2_NBUF > ((uint32)(j*1024./(n>>7)) - 5) )
						fprintf(stderr,"WARNING: User-specified maximum number of Stage 2 buffers may exceed %u MB of available RAM.\n",j);
				}
				ASSERT(PM1_S2_NBUF >= 24,"p-1 Stage 2 requires at least 24 residue-sized memory buffers!\n");
				// See if S2 restart file exists:
				strcpy(cstr,RESTARTFILE); cstr[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f'); strcat(cstr, ".s2");
				// If a regular (non-continuation, i.e. B2_start = B1) stage 2 and S2 restart file exists, read
				// bounds and bigstep D from logfile to infer relocation-prime psmall from D and seed the call to
				// pm1_bigstep_size() with that to ensure our restart-run bigstep shares the same relocation-prime:
				if(!s2_continuation) {
					fp = mlucas_fopen(cstr,"rb");
					if(fp) {
						// This snip reads the relocation-prime from the high byte of the nsquares field, byte 10 of the S2 savefile:
						i = fgetc(fp);
						if(!test_types_compatible(i, TEST_TYPE)) {
							snprintf(cbuf,STR_MAX_LEN*2, "%s: TEST_TYPE != fgetc(fp)\n",cstr); ASSERT(0,cbuf);
						}
						if((i = fgetc(fp)) != MODULUS_TYPE) {
							snprintf(cbuf,STR_MAX_LEN*2, "ERROR: %s: MODULUS_TYPE != fgetc(fp)\n",cstr); ASSERT(0,cbuf);
						}
						itmp64 = 0ull; 	for(j = 0; j < 8; j++) { i = fgetc(fp);	itmp64 += (uint64)i << (8*j); }
						fclose(fp); fp = 0x0;
						if(i != EOF)	// Needed to handle case where .s2 file was touched but ended up empty or < 10 bytes long
							psmall = i;
						itmp64 &= 0x00FFFFFFFFFFFFFFull;	// Mask off psmall to get stage 2 q of checkpoint data
						fprintf(stderr,"Read iter = %llu and relocation-prime psmall = %u from savefile %s.\n",itmp64,psmall,cstr);
						// Now parse logfile to get proper B2 and validate corresponding B2_start vs B2/[psmall from .s2 file].
						// Logfiles can be messy and include one or more aborted-restarts; we want the last B2_start-containing
						// entry followed by a savefile-write entry, as inferred from presence of a "% complete" substring:
						j = filegrep(STATFILE,"% complete",cbuf,-1);	// Trailing -1 means return last such match, if any, in cbuf
						filegrep(STATFILE,"B2_start = ",cbuf,j);	// Trailing j-arg means return last such match before line j, if any, in cbuf
						// If match "B2_start =", read bigstep D from match-line and infer relocation-prime psmall from D:
						if(strlen(cbuf)) {
							char_addr = strstr(cbuf,"B2_start = ");
							B2_start = (uint64)strtoull(char_addr+11, &cptr, 10);	ASSERT(B2_start != -1ull, "strtoull() overflow detected.");
							char_addr = strstr(cbuf,"B2 = ");
							B2 = (uint64)strtoull(char_addr+5, &cptr, 10);	ASSERT(B2 != -1ull, "strtoull() overflow detected.");
							char_addr = strstr(cbuf,"Bigstep = ");
							if(char_addr) {
								i = strtoul(char_addr+10, &endp, 10);
								if((i%210) == 0)
									i = 11;
								else if((i%330) == 0)
									i = 7;
								else {
									fprintf(stderr,"WARNING: Bigstep value %u read from logfile %s unsupported ... ignoring.\n",i,cstr);
									i = 0;
								}
							}
							// Now compare the params from the restartfile vs those captured in the log:
							if(psmall)
								ASSERT(psmall == i && B2_start == B2/psmall, "Stage 2 params mismatch those captured in the .stat logfile!");
							else
								psmall = i;
							// If stage 2 q of checkpoint >= B2, proceed directly to GCD:
							if(itmp64 >= B2)
								goto PM1_STAGE2_GCD;
							// Lastly, must reset B2_start = B1, since stage 2 code expects that to properly (re)init relocation-params:
							B2_start = B1;
						}
					}	// endif( S2 restart file exists? )
				}
				// Now feed any relocation-prime value (psmall, stored in i here) to the pm1_bigstep_size() call, which resets
				// PM1_S2_NBUF to the largest S2 buffer count <= (input PM1_S2_NBUF value) which is also compatible with psmall:
				pm1_bigstep_size(&PM1_S2_NBUF, &pm1_bigstep, &pm1_stage2_mem_multiple, psmall);
				fprintf(stderr,"Using %u Stage 2 memory buffers, D = %u, M = %u, psmall = %u.\n",PM1_S2_NBUF,pm1_bigstep,pm1_stage2_mem_multiple,psmall);
				double*mult[4] = {b,c,d,e};	// Pointers to scratch storage in form of 4 residue-length subarrays
				ierr = pm1_stage2(p, pm1_bigstep, pm1_stage2_mem_multiple,
					a,mult,arrtmp,	// Pointers stage 1 residue in a[], double** scratch storage 4-vector mult[], arrtmp
					func_mod_square, n, scrnFlag, &tdiff, gcd_str);
				if(ierr && (ierr % ERR_ROUNDOFF) == 0) {	// Workaround for pthreaded-run issue where a single ROE sometimes shows up as an integer multiple of ERR_ROUNDOFF
					n = get_nextlarger_fft_length(n);	kblocks = (n >> 10);
					// Clear out current FFT-radix data, since get_preferred_fft_radix() expects that:
					for(i = 0; i < NRADICES; i++) { RADIX_VEC[i] = 0; }		NRADICES = 0;
					goto SETUP_FFT;
				} else if(ierr == ERR_INTERRUPT) {
					// First print the signal-handler-generated message:
					mlucas_fprint(cbuf,1);
					exit(1);
				} else if(ierr) {
					sprintf(cbuf,"p-1 stage 2 hit an unhandled error of type[%u] = %s! Aborting.",ierr,returnMlucasErrCode(ierr));
					ASSERT(0,cbuf);
				}
				// If gcd_str non-empty on return, it means one of the intermediate S2 GCDs turned up a factor,
				// prompting an early-return, In this case the S2 code will have reset B2 to reflect the actual interval run.
				// Otherwise do end-of-scheduled-S2 GCD - S2 residue returned in arrtmp, no need to call convert_res_FP_bytewise():
			PM1_STAGE2_GCD:
				if(strlen(gcd_str)) {
					s2_partial = TRUE;	// Clue the JSON-generating function to partial-s2-ness
				} else {
					j = (p+63+(MODULUS_TYPE == MODULUS_TYPE_FERMAT))>>6;
					i = gcd(2,p,arrtmp,0x0,j,gcd_str);	// 1st arg = stage just completed
				}
				// If factor found, gcd() will have done needed status-file-writes; unlike S1, leave B2 as-is, need it for the JSON output.
				// Write JSON output and go to next assignment:
				calendar_time = time(NULL);
				gm_time = gmtime(&calendar_time);
				if(!gm_time)	// If UTC not available for some reason, just substitute the local time:
					gm_time = localtime(&calendar_time);
				strftime(timebuffer,SIZE,"%Y-%m-%d %H:%M:%S UTC",gm_time);
				generate_JSON_report(0,p,n,0ull,NULL,timebuffer, B1,B2,gcd_str,s2_partial, cstr);	// cstr holds JSONified output
				snprintf(cbuf,STR_MAX_LEN*2, "If using the manual results submission form at mersenne.org, paste the following JSON-formatted results line:\n%s\n",cstr);
				mlucas_fprint(cbuf,0);
				fp = mlucas_fopen(OFILE,"a");
				if(fp){
					fprintf(fp,"\n%s",cstr); fclose(fp); fp = 0x0;
				}
			}
		} else {
			fprintf(stderr,"User specified low-mem run mode ... no stage 2.\n");
		}
	} else {
		ASSERT(0, "Unrecognized test type!");
	}	/* endif(TEST_TYPE == TEST_TYPE_PRIMALITY) */

	/*...If successful completion, delete the secondary restart files...save the primary in case it's a prime,
	so we can re-run the last (p % ITERS_BETWEEN_CHECKPOINTS) iterations by way of quick partial verification: */
	/* v20.1: If not an assignment-to-split and just-completed test was a regular (non-stage-2-continuation) p-1
	run, rename the primary savefile add a ".s1" to the name so as to not collide with any ensuing LL/PRP test.
	***BUT:*** If our restart ended up being from the secondary q-savefile due to a missing/corrupted p-savefile,
	we want to rename *that* one to [p|f][expo].s1.
		And in fact we have a similar scenario for LL/PRP tests - we want to preserve 1 copy of final savefile
	in case we need to rerun-last-few-thousand iters for an alleged prime discovery, so fiddle logic for both cases.
	*/
	// Start with the primary savefile, loop over both, if needed
	RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
	if(TEST_TYPE == TEST_TYPE_PM1) {
		/*
		Dec 2021: If just-completed p-1 run, delete any .s2 savefile created during stage 2, in case user wants to do
		one or more stage 2 continuation runs. Using the s2 residue instead of the s1 should be OK in such cases, but
		it offers no advantage in terms of covering the s2 primes, and if the s2 residue was somehow silently corrupted,
		re-using it would kibosh subsequent stage 2 continuation runs. Safer to start with s1 residue for those:
		*/
		strcpy(cstr, RESTARTFILE); strcat(cstr, ".s2");
		if(remove(cstr)) {
			snprintf(cbuf,STR_MAX_LEN*2,"INFO: Unable to remove stage 2 savefile %s.\n",cstr);
			mlucas_fprint(cbuf,1);
		}
		if(!s2_continuation) {
			strcpy(cstr, RESTARTFILE); strcat(cstr, ".s1");	// cstr = [p|f][expo].s1
		}
	} else if(TEST_TYPE == TEST_TYPE_PRIMALITY || TEST_TYPE == TEST_TYPE_PRP) {
		strcpy(cstr, RESTARTFILE); cstr[0] = 'q';		// cstr = q[expo]
	}
	for(ierr = 0; ; RESTARTFILE[0] = 'q') {	// Start with the p-savefile, inrement to q-savefile on looping
		if(restart_file_valid(RESTARTFILE, p, (uint8*)arrtmp, (uint8*)e_uint64_ptr)) {
			// If end of a regular (non-s2-continuation) p-1 run and primary good, rename it from [p|f][expo] ==> [p|f][expo].s1;
			// if primary missing/corrupt, rename secondary q[expo] ==> [p|f][expo].s1:
			if(TEST_TYPE == TEST_TYPE_PM1 && !s2_continuation) {
				if(rename(RESTARTFILE, cstr)) {
					snprintf(cbuf,STR_MAX_LEN*2,"ERROR: unable to rename the p-1 stage 1 savefile %s ==> %s ... any ensuing LL/PRP test will overwrite.\n",RESTARTFILE,cstr);
					mlucas_fprint(cbuf,1);
				}
			} else if(TEST_TYPE == TEST_TYPE_PRIMALITY || TEST_TYPE == TEST_TYPE_PRP) {
				// If primary was missing/corrupt, i.e. we're on the secondary, rename secondary to primary-name
				if(RESTARTFILE[0] == 'q') {
					RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
					if(rename(cstr, RESTARTFILE)) {
						snprintf(cbuf,STR_MAX_LEN*2,"ERROR: Primary savefile missing/corrupt, but unable to rename the secondary %s ==> %s ... any ensuing LL/PRP test will overwrite.\n",RESTARTFILE,cstr);
						mlucas_fprint(cbuf,1);
					}
				} else if(remove(cstr))	// ...otherwise delete the secondary
					fprintf(stderr, "Unable to delete secondary restart file %s.\n",cstr);
			}
		} else
			ierr = 1;
		// If no errors on primary, break; otherwise loop
		if(!ierr || RESTARTFILE[0] == 'q')
			break;
	}
	// If completion of LL/PRP run, or secondary was not used as a backup for p-1 stage 1 savefile rename, delete it now:
	RESTARTFILE[0] = 'q';
	if(remove(RESTARTFILE))
		fprintf(stderr, "Unable to delete secondary restart file %s\n",RESTARTFILE);

	RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');

	/*...If in non-interactive (production) mode, delete the just-completed assignment from the workfile.
	Do this by copying all the currently assigned exponents (except the first) from the workfile to a temporary file,
	then deleting the workfile and replacing it with the temp-file.
	v20: Added logic to take Primality or PRP-test assignment with p-1 still needing to to be done and split into:
		1. A Pminus1 assignment - Should this find a factor, both assignments are subsequently deleted;
		2. A Primality or PRP-test assignment with p-1 already done.
	*/
GET_NEXT_ASSIGNMENT:

	// If workfile is still open, close it:
	if(fp) { fclose(fp); fp = 0x0; }

	if(!INTERACT) {
		//*** IN THIS CASE MUST MAKE SURE CBUF,CSTR ONLY GET OVERWRITTEN ON ERROR ERROR, SINCE THEY CONTAIN THE SPLIT ASSIGNMENT! ***
		if(split_curr_assignment) {
			sprintf(ESTRING,"%llu",p);	// Set ESTRING here, as this bypasses the normal route for getting to GET_NEXT_ASSIGNMENT
			ASSERT(TEST_TYPE == TEST_TYPE_PM1,"GET_NEXT_ASSIGNMENT: split_curr_assignment = TRUE, but TEST_TYPE != PM1.");
		}

		fp = mlucas_fopen(WORKFILE,"r");
		if(!fp) {
			sprintf(cbuf,"ERROR: unable to open %s file for reading.\n",WORKFILE);
			ASSERT(0,cbuf);
		}
		/* Remove any WINI.TMP file that may be present: */
		remove("WINI.TMP");
		fq = mlucas_fopen("WINI.TMP", "w");
		if(!fq) {
			sprintf(cbuf, "Unable to open WINI.TMP file for writing.\n");
			ASSERT(0,cbuf);
		}

	GET_NEXT:
		/* Delete or suitably modify current-assignment line (line 1) of worktodo file: */
		i = 0;	// This counter tells how many *additional* assignments exist in worktodo
		if(!fgets(in_line, STR_MAX_LEN, fp)) {
			sprintf(cbuf, "ERROR: %s file not found at end of current-assignment processing\n", WORKFILE);
			ASSERT(0,cbuf);
		}
		// v20.1.1: Parse all lines whose 1st non-WS char is alphabetic;
		char_addr = in_line;	j = 0;
		while(isspace(in_line[j])) { ++j; }
		char_addr += j;
		if(!isalpha(in_line[j]))
			goto GET_NEXT;

		// Look for m in first eligible assignment; for F[m], need to also look for 2^m in case assignment is in KBNC format:
		if(!strstr(in_line, ESTRING) && !(MODULUS_TYPE == MODULUS_TYPE_FERMAT && strstr(in_line, BIN_EXP)) ) {
			snprintf(cbuf,STR_MAX_LEN*2, "ERROR: Current exponent %s not found in line 1 of %s file - quitting.\n", ESTRING, WORKFILE);
			ASSERT(0,cbuf);
		} else {
			/* If we just finished the TF or p-1 preprocessing step of an LL or PRP test,
			update the current-assignment line to reflect that and write it out: */
			if(strstr(in_line, "PRP") || strstr(in_line, "Test") || strstr(in_line, "DoubleCheck")) {
			#if INCLUDE_TF
				if(TEST_TYPE == TEST_TYPE_TF) {
					/* Factor depth assumed to follow the first comma in in_line: */
					char_addr = strstr(char_addr, ",");
					ASSERT(char_addr != 0x0,"Null char_addr");
					sprintf(++char_addr, "%u", TF_BITS);
					fputs(in_line, fq);
				}
			#endif
				// This imples TEST_TYPE == TEST_TYPE_PM1; note that this flag gets cleared on cycling back to RANGE_BEG:
				if(split_curr_assignment) {
					/*0/1 flag indicating whether P-1 has been done assumed to follow second comma in in_line: */
					fputs(cbuf, fq);	// The Pminus1 assignment
					fputs(cstr, fq);	// The PRP or LL assignment, with trailing 0 indicating p-1 done (true by time we get to it)
					i = 2;	// And reset remaining-assignments counter
				}
			} else if(stristr(in_line, "pminus1")) {
				// If current p-1 assignment found a factor and resulted from splitting of a PRP/LL assignment -
				// note that split_curr_assignment == TRUE only at time of the initial splitting - delete them both:
				ASSERT(TEST_TYPE == TEST_TYPE_PM1,"GET_NEXT_ASSIGNMENT: current assignment is Pminus1=, but TEST_TYPE != PM1.");
				if(strlen(gcd_str) != 0) {	// Found a factor?
					char_addr = strstr(in_line, "=");	ASSERT(char_addr != 0x0,"Malformed assignment!");
					char_addr++;
					if(is_hex_string(char_addr, 32)) {
						strncpy(aid,char_addr,32);
					} else if(STREQN_NOCASE(char_addr,"n/a",3)) {
						strncpy(aid,char_addr, 3);
					} else {
						snprintf(cbuf,STR_MAX_LEN*2,"INFO: Assignment \"%s\" lacks a valid assignment ID ... proceeding anyway.\n",in_line);
						mlucas_fprint(cbuf,1);
						aid[0] = '\0';	// This guarantees that the strstr(in_line,aid) part of on the next-assignment search below succeeds.
					}
					// If next assignment exists, is an LL/PRP, and has same exponent and AID, delete it (by doing nothing), otherwise
					// copy it to the fq-file. We use both the AID and binary exponent here, since the former could be some made-up
					// value (e.g. all-0s) and the latter has a nonzero chance of appearing in the AID of an unrelated next-assignment:
					in_line[0] = '\0';	// Careful here - if next fgets(in_line) fails, are left with old in_line, must zero it first.
					if(fgets(in_line, STR_MAX_LEN, fp)) {
						if( (strstr(in_line, "PRP") || strstr(in_line, "Test") || strstr(in_line, "DoubleCheck"))
						 && strstr(in_line,ESTRING) && strstr(in_line,aid) )
						{
							/* Lose the assignment (by way of no-op) */
						} else {
							fputs(in_line, fq); i = 1;	// Copy PRP/LL assignment and reset remaining-assignments counter
						}
					}
				}
			}
			/* Otherwise lose the current line (by way of no-op) */
		}
		/* Copy the remaining ones; */
		while(fgets(in_line, STR_MAX_LEN, fp))
		{
			fputs(in_line, fq);	++i;
		}
		fclose(fp); fp = 0x0;
		fclose(fq); fq = 0x0;

		/* Now blow away the old worktodo file and rename WINI.TMP ==> worktodo.txt...	*/
		remove(WORKFILE);
		if(rename("WINI.TMP", WORKFILE))
		{
			sprintf(cbuf,"ERROR: unable to rename WINI.TMP file ==> %s ... attempting line-by-line copy instead.\n",WORKFILE);
			fprintf(stderr,"%s",cbuf);

			/* If attempting to simply rename the TMP file fails, do it the hard way: */
			fp = mlucas_fopen(WORKFILE,"w");
			if(!fp) {
				sprintf(cbuf,"ERROR: unable to open %s file for writing.\n", WORKFILE);
				ASSERT(0,cbuf);
			}

			fq = mlucas_fopen("WINI.TMP", "r");
			if(!fq) {
				sprintf(cbuf,"Unable to open WINI.TMP file for reading.\n");
				ASSERT(0,cbuf);
			}
			while(fgets(in_line, STR_MAX_LEN, fq)) {
				fputs(in_line, fp);
			}
			fclose(fp); fp = 0x0;
			fclose(fq); fq = 0x0;

			/*...Then remove the WINI.TMP file:	*/

			remove("WINI.TMP");
		}
		/* if one or more exponents left in rangefile, go back for more; otherwise exit. */
		if (i > 0) {
			// Reset any globals which are optionally user-forced here; e.g. PM1_S2_NBUF can be user or auto-set,
			// no easy way to determine which it was, so only respect for first p-1 stage 2 done since program start:
			PM1_S2_NBUF = 0;
			goto RANGE_BEG;   /* CYCLE RANGE */
		} else
			exit(EXIT_SUCCESS);
	}

	/* If it was an interactive run, return: */
	return(ierr);
}

/***************/

/* Do some basic self-tests: */
void 	Mlucas_init(void)
{
	uint32 i,j,n;

	host_init();

	/* Set approved test types in ASSIGNMENT_TYPE_MATRIX matrix: */
	for(i = 0; i < MODULUS_TYPE_DIM; i++)
	{
		for(j = 0; j < TEST_TYPE_DIM; j++)
		{
			ASSIGNMENT_TYPE_MATRIX[i][j] = 0;
		}
	}
#if INCLUDE_TF
	/* TF currently only supported for Mersennes: */
	ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_MERSENNE   ][TEST_TYPE_TF] = TRUE;
#endif
	ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_MERSENNE   ][TEST_TYPE_PRIMALITY     ] = TRUE;
	ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_FERMAT     ][TEST_TYPE_PRIMALITY     ] = TRUE;

	/* Simple self-tester for get_fft_radices(): */
	fprintf(stderr, "INFO: testing FFT radix tables...\n");
	test_fft_radixtables();

	/* Set min. exponent (in terms of power of 2) that can be tested: */
	/* Check that the purported min. FFT length is actually supported: */
	ASSERT(get_fft_radices(MIN_FFT_LENGTH_IN_K, 0, 0x0, 0x0, 0) == 0,"Require get_fft_radices(MIN_FFT_LENGTH_IN_K, 0) == 0");
	n = (MIN_FFT_LENGTH_IN_K << 10);
	/* Make sure N didn't overflow */
	ASSERT((n >> 10) == MIN_FFT_LENGTH_IN_K,"Require (n >> 10) == MIN_FFT_LENGTH_IN_K");
	PMIN = 2*n;	/* 2 bits per input is about the smallest we can test without getting nonzero-carry errors */

	/* Set max. exponent (in terms of power of 2) that can be tested: */
	/* Check that the purported max. FFT length is actually supported: */
	ASSERT(get_fft_radices(MAX_FFT_LENGTH_IN_K, 0, 0x0, 0x0, 0) == 0,"Require get_fft_radices(MAX_FFT_LENGTH_IN_K, 0) == 0");
	n = (MAX_FFT_LENGTH_IN_K << 10);
	/* Make sure N didn't overflow */
	ASSERT((n >> 10) == MAX_FFT_LENGTH_IN_K,"Require (n >> 10) == MAX_FFT_LENGTH_IN_K");
	PMAX = 1.05*given_N_get_maxP(n);	// Allow same wiggle room here as in ernstMain

	ASSERT(PMAX > PMIN,"Require PMAX > PMIN");

#if INCLUDE_TF
	/* Simple self-tester for sieve factoring routines: */
	fprintf(stderr, "INFO: testing trial-factoring routines...\n");
	if(test_fac() != 0)
	{
		sprintf(cbuf, "Mlucas_init : Trial-factoring self-test failed.\n");
		fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
	}
#endif
#if 0	// v20: Use GMP GCD, own-rolled O(n*(log n)^2) one simply not in the cards.
	/* Simple self-tester for GCD routines in gcd_lehmer.c: */
	fprintf(stderr, "INFO: testing GCD routines...\n");
	if(test_gcd() != 0)
	{
		sprintf(cbuf, "Mlucas_init : GCD test failed.\n");
		fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
	}
#endif
}

/********************/

/* Needed to support circularly-shifted residues - given a bitshift value, finds the padded-array
index of the word containing the given bit and computes the corresponding carry-in (which will be either
the initial seed at start of the run for LL and PRP tests, or the -2 added to the convolution-
based mod-squaring output on each ensuing LL iteration) needing to be injected into said word
at the appropriate within-word bit offset.

If first_entry == TRUE, inits the BIGWORD_BITMAP and BIGWORD_NBITS arrays and directly injects
the carry into the residue array a[], a valid pointer to which is assumed to be supplied in the arglist.
In this case we have cy_in != 0, and we also forward-weight the resulting shifted-carryin
in preparation for the initial forward FFT. (This single-word init-and-fwd-weighting saves us from having
to run through the whole-array version of same in mers_mod_square() used for fully-populated interim
LL residues read from a savefile .)

If first_entry == FALSE, uses the previously-inited BIGWORD_BITMAP and BIGWORD_NBITS arrays and
does a fast lookup based on those. In this case whether the carry is injected into the residue array
is up to the caller. Most commonly these subsequent calls will be made from one of the DFT/carry functions,
in which case we do not have time for a leisurely initial-residue-array-forward-weighting as on the initial
iteration, and thus the caller should set cy_in = 0.0. Rather, the carry routine is expected to take the
index+shift data and use them to inject the carry into the carry-macro call sequence in appropriate fashion,
e.g. by on-the-fly computing the corresponding forward DWT weight for the array element into which the carry
is injected, multiplying cy_in*(n/2)*wt_fwd to match the scaling of convolution outputs prior to the carry-
macro step, adding that scaled carry to the appropriated convolution output and then doing the carry-macro
step as usual.

OTOH during e.g. self-tests we will also be calling this function at run-start for each of the various
sets of FFT radices being tried at a given FFT length, in which case first_entry = false but we still
want to init a single main-array element based on the properly shifted initial LL residue. In this case
we have cy_in != 0, and we use both the fast lookup scheme and forward-weight the resulting shifted-carryin
in preparation for the initial forward FFT.

Return value:
o High 7 bytes store *unpadded* index of double-array word containing the bit of the specified shift;
o Low byte stores bit-within-that-double-array word of the specified shift.
*/
uint64 	shift_word(double a[], int n, const uint64 p, const uint64 shift, const double cy_in)
{
	static int first_entry=TRUE, nsave = 0;
	static uint64 psave = 0ull;
	static double bits_per_word = 0.0, words_per_bit = 0.0;
	int pow2_fft,bimodn,curr_bit64,curr_wd64,w64,curr_wd_bits,mod64,findex,i,j,j1,j2;
	static uint32 bw = 0,sw = 0,nwt = 0,sw_div_n,bits[2];
	uint32 sw_idx_modn,ii;	// Aug 2021: ii needs to be unsigned for the (shift < ii) compare when shift >= 2^31
	double cy,theta,wt_fwd,wt_cos,wt_sin;
	uint64 nbits, itmp64;
	 int64 retval = -1;	// Make this signed to ease "not yet set?" check
#ifdef USE_FGT61
	ASSERT(0,"shift_word() needs to be modified to support FGT!");
#endif
	if(n != nsave || p != psave) {
		first_entry = TRUE;	for(j = 0; j < (n>>6); j++) { BIGWORD_BITMAP[j] = 0ull; }	// Need to clear bitmap in case of multi-FFT-length run
		bw = p%n; sw = n-bw;
		/* If Fermat number, make sure exponent a power of 2: */
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
			ASSERT(TRANSFORM_TYPE == RIGHT_ANGLE,"Require TRANSFORM_TYPE == RIGHT_ANGLE");
			findex = trailz64(p);
			ASSERT((p >> findex) == 1,"Require (p >> findex) == 1");
			/* For Fermat-mod, only need IBDWT weights table if it's a non-power-of-2-length transform, in which
			case the table has {nwt = odd part of N} distinct elements. Avoid if() logic related to power-of-2-or-not
			by initing a single DWT weight = 1.0 in the power-of-2 case and = 2^((j%nwt)/n) otherwise:
			*/
			nwt = (n >> trailz32(n));
			sw_div_n = sw*nwt/n;
		}
		else
			ASSERT(TRANSFORM_TYPE == REAL_WRAPPER,"Require TRANSFORM_TYPE == REAL_WRAPPER");

		/* Vector length a power of 2? */
		pow2_fft = (n >> trailz32(n)) == 1;

		bits[0] = p/n;
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT && pow2_fft == TRUE)
			bits[1] =     bits[0];
		else
			bits[1] = 1 + bits[0];
	}
	// May 2018: Check that BIGWORD_BITMAP and BIGWORD_NBITS arrays have been alloc'ed and use fast lookup based on those.
	// Timing loop on my Core2 macbook indicates this fast-lookup needs ~160 cycles @1Mdouble FFT, not horrible but slower
	// than I'd like, likely due to cache impacts of doing random-word lookups in the resulting 128kB and 64kB BIGWORD* arrays.
	// Also had the "adjusting..." printfs enabled during the timing tests, 0 such adjustments needed for 10^9 random-shifts:
	if(!first_entry) {
	//	ASSERT(BIGWORD_BITMAP != 0x0 && BIGWORD_NBITS != 0x0, "BIGWORD_BITMAP and BIGWORD_NBITS arrays not alloc'ed!");
		// Divide [shift] by the average bits per word to get a quick estimate of which word contains the corresponding bit:
		j = shift*words_per_bit;	w64 = j>>6; mod64 = j&63;
		// Then exactly compute the bitcount at the resulting word, by adding the BIGWORD_NBITS-array-stored exact
		// total bitcount at the next-lower index-multiple-of-64 to the number of bits in the next (mod64) words,
		// not including the current word, hence (64-mod64) rather than (63-mod64) as the BIGWORD_BITMAP[w64] shift count.
		// NOTE: That (64 - ...) shiftcount also necessitates an added AND-mask to zero the result in the mod64 == 0 case:
		itmp64 = BIGWORD_BITMAP[w64];	i = (int)bits_per_word;
		ii = BIGWORD_NBITS[w64] + i*mod64 + popcount64( (itmp64<<(64-mod64)) & -(mod64 != 0) );
		// Loop up or down (should be by at most 1 word in either direction) if the resulting #bits is such that RES_SHIFT maps to a <> word:
		curr_wd_bits = i + ( (int64)(itmp64<<(63-mod64)) < 0 );
		// Can gain a few % speed by commenting out this correction-step code, but even though I've encountered
		// no cases where it's used in my (admittedly quite limited) testing, better safe than sorry:
		if(shift < ii) {
		//	printf("shift[%llu] < ii [%u] ... adjusting downward.\n",shift,ii);
			while(shift < ii) {
				if(--j < 0) {	// Note j is signed
					j += 64;	w64 = j>>6; mod64 = j&63;	// Go to next-lower word of BIGWORD_BITMAP
				} else {
					--mod64;
				}
				curr_wd_bits = i + ( (int64)(BIGWORD_BITMAP[w64]<<(63-mod64)) < 0 );
				ii -= curr_wd_bits;
			}
		} else if(shift >= (ii + curr_wd_bits) ) {
		//	printf("shift[%llu] >= (ii + curr_wd_bits) [%u] ... adjusting upward.\n",shift,(ii + curr_wd_bits));
			while(shift >= (ii + curr_wd_bits) ) {
				if(++j > 63) {
					j -= 64;	w64 = j>>6; mod64 = j&63;	// Go to next-higher word of BIGWORD_BITMAP
				} else
					++mod64;
				curr_wd_bits = i + ( (int64)(BIGWORD_BITMAP[w64]<<(63-mod64)) < 0 );
				ii += curr_wd_bits;
			}
		}
		// Must account for the right-angle-transform data layout here:
		if(TRANSFORM_TYPE == RIGHT_ANGLE) {
			j <<= 1;
			if(j >= n) {
				j -= n; i = 1;
			} else {
				i = 0;
			}
		}
	#ifdef USE_AVX512
		j1 = (j & mask03) + br16[j&15];
	#elif defined(USE_AVX)
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		// Carry-in needs to be properly left-shifted w.r.to the residue word into which it goes;
		// don't attempt to propagate any carry which would result from a proper base-normalization of the shifted carryin:
		// Since the stdlib pow() function is expensive (~100-200 cycles), wrap this bit in an if(cy_in):
		if(cy_in) {
			if(TRANSFORM_TYPE == REAL_WRAPPER) {
				// Compute wt = 2^(j*sw % n)/n:
				sw_idx_modn = ((uint64)j*sw) % n;	// N is 32-bit, so only use 64-bit to hold intermediate product
				wt_fwd = pow(2.0, (double)sw_idx_modn/n);
				a[j1] = cy_in * (double)(1ull << (shift - ii)) * wt_fwd;
			} else {
				j2 = j1 + RE_IM_STRIDE;
				// Compute both the DWT weight (if n = 2^k this = 1) and the acyclic-twisting root of unity:
				sw_idx_modn = ((uint64)(j>>1)*sw_div_n) % nwt;	// Analogous loop in fermat_mod_square increments ii += SW_DIV_N (mod nwt) for every j += 2, so need to halve j here
				wt_fwd = pow(2.0, (double)sw_idx_modn/nwt);
				wt_fwd *= (double)(1ull << (shift - ii));
				theta  = j * qfdbl(QPIHALF) / n;
				wt_cos = wt_fwd*cos(theta);	wt_sin = wt_fwd*sin(theta);
				if(!i) {	// [Re,0] * exp((j/2)*Pi/n) = [ Re*cos, Re*sin]:
					a[j1] =  cy_in*wt_cos;	a[j2] = cy_in*wt_sin;
				} else {	// [0,Im] * exp((j/2)*Pi/n) = [-Im*sin, Im*cos]:
					a[j1] = -cy_in*wt_sin;	a[j2] = cy_in*wt_cos;
				}
			}
		}
		retval = ((uint64)j<<8) + (shift - ii);
		goto DONE;
	}

	first_entry = FALSE;
	psave = p; nsave = n; bits_per_word = (double)p/n; words_per_bit = 1.0/bits_per_word;
	ASSERT(MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");
	ASSERT(TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");
	ASSERT(p > shift,"Specified shift count out of range!");

	nbits = 0;	/* Total bits accumulated so far in the residue words processed	*/

	if(TRANSFORM_TYPE == REAL_WRAPPER)
	{
		bimodn = 0;
		ii = 0;		/* Index storing whether current word is a bigword (ii = 1) or a littleword (ii = 0). */
		/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
		if(bw > 0)
			ii = 1;
		curr_wd64 = -1; curr_bit64 = 0;
		for(j = 0; j < n; j++)
		{
			// Every 64th word create an entry in how-many-residue-bits-in-all-elts-below-this-one array, incr curr_wd64 and reset curr_bit64:
			if(!(curr_bit64 & 63)) {
				BIGWORD_NBITS[++curr_wd64] = nbits;
				curr_bit64 = 0;
			}
			nbits += bits[ii];
			BIGWORD_BITMAP[curr_wd64] |= ((uint64)ii) << curr_bit64++;
			// Inject the carry-in and set the return value:
			if(nbits > shift) {
				if(retval < 0) {	// retval has not yet been set
					curr_wd_bits = shift - (nbits - bits[ii]);	retval = ((uint64)j<<8) + curr_wd_bits;
					cy = cy_in;
				//	printf("Hit target bit %llu in a[%u] (=> BIGWORD_BITMAP[%u]), bit %u of <0:%u>, bitmap-word bit = %u\n",shift,j,curr_wd64,curr_wd_bits,bits[ii]-1,curr_bit64-1);	ASSERT(curr_wd_bits <= bits[ii]-1,"GAH!");
				}
			#ifdef USE_AVX512
				j1 = (j & mask03) + br16[j&15];
			#elif defined(USE_AVX)
				j1 = (j & mask02) + br8[j&7];
			#elif defined(USE_SSE2)
				j1 = (j & mask01) + br4[j&3];
			#else
				j1 = j;
			#endif
				j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				// Carry-in needs to be properly left-shifted w.r.to the residue word into which it goes;
				// don't attempt to propagate any carry which would result from a proper base-normalization of the shifted carryin:
				if(cy) {	// Compute wt = 2^(j*sw % n)/n:
					sw_idx_modn = ((uint64)j*sw) % n;	// N is 32-bit, so only use 64-bit to hold intermediate product
					wt_fwd = pow(2.0, (double)sw_idx_modn/n);
					a[j1] = cy * (double)(1ull << (shift - (nbits - bits[ii]))) * wt_fwd;
					cy = 0.0;
				}
			}
			bimodn += bw;
			if(bimodn >= n) bimodn -= n;
			ii = (uint32)(sw - bimodn) >> 31;
		}
	}
	else
	{
		ASSERT(TRANSFORM_TYPE == RIGHT_ANGLE, "Invalid or uninited TRANSFORM_TYPE!");
		curr_wd64 = -1; curr_bit64 = 0;
	  for(i = 0; i < 2; i++)	// Two stride-2 loops to cover even and odd-indexed array elements, respectively:
	  {
		bimodn = n;
		for(j = 0; j < n; j += 2)
		{
			ii = (bimodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */
			bimodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */
			bimodn += ( ((int)ii-1) & n);		/*       add 0 if a bigword,   N if a smallword */
			// Every 64th word create an entry in how-many-residue-bits-in-all-elts-below-this-one array, incr curr_wd64 and reset curr_bit64:
			if(!(curr_bit64 & 63)) {
				BIGWORD_NBITS[++curr_wd64] = nbits;
				curr_bit64 = 0;
			}
			nbits += bits[ii];
			BIGWORD_BITMAP[curr_wd64] |= ((uint64)ii) << curr_bit64++;
			// Inject the carry-in and set the return value:
			if(nbits > shift) {
				if(retval < 0) {	// retval has not yet been set
					curr_wd_bits = shift - (nbits - bits[ii]);	retval = ((uint64)j<<8) + curr_wd_bits;
					cy = cy_in;
				//	printf("Hit target bit %llu in a[%u] (=> BIGWORD_BITMAP[%u]), bit %u of <0:%u>, bitmap-word bit = %u\n",shift,j,curr_wd64,curr_wd_bits,bits[ii]-1,curr_bit64-1);	ASSERT(curr_wd_bits <= bits[ii]-1,"GAH!");
				}
			#ifdef USE_AVX512
				j1 = (j & mask03) + br16[j&15];
			#elif defined(USE_AVX)
				j1 = (j & mask02) + br8[j&7];
			#elif defined(USE_SSE2)
				j1 = (j & mask01) + br4[j&3];
			#else
				j1 = j;
			#endif
				j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1 + RE_IM_STRIDE;
				// Carry-in needs to be properly left-shifted w.r.to the residue word into which it goes;
				// don't attempt to propagate any carry which would result from a proper base-normalization of the shifted carryin:
				if(cy) {	// Compute both the DWT weight (if n = 2^k this = 1) and the acyclic-twisting root of unity:
					sw_idx_modn = ((uint64)(j>>1)*sw_div_n) % nwt;	// Analogous loop in fermat_mod_square increments ii += SW_DIV_N (mod nwt) for every j += 2, so need to halve j here
					wt_fwd = pow(2.0, (double)sw_idx_modn/nwt);
					cy *= (double)(1ull << (shift - (nbits - bits[ii])));
					cy *= wt_fwd;
					theta  = j * qfdbl(QPIHALF) / n;
					wt_cos = cos(theta);	wt_sin = sin(theta);
					if(!i) {	// [Re,0] * exp((j/2)*Pi/n) = [ Re*cos, Re*sin]:
						a[j1] =  cy*wt_cos;	a[j2] = cy*wt_sin;
					} else {	// [0,Im] * exp((j/2)*Pi/n) = [-Im*sin, Im*cos]:
						a[j1] = -cy*wt_sin;	a[j2] = cy*wt_cos;
					}
				//	printf("shift_word(): set a[%u] = %20.12f, a[%u] = %20.12f\n",j1,a[j1],j2,a[j2]);
					cy = 0.0;
				}
			}
		}
	  }
	}
DONE:
	return retval;
}

/*
Suyama cofactor-PRP test for Fermat/Mersenne modulus N, using base-3 PRP-test to illustrate:
	Let P = 3^[(N-1)/2] (mod N) be the base-3 Euler (Pepin-test) residue for N;
	Let F = product of known small-prime factors of N, C the corresponding cofactor, i.e. N = F*C;
	Let A = 3^(N-1) = P^2 (mod N) be the base-3 Fermat-PRP residue for N;
	Let B = 3^(F-1) (mod N).
Then C is composite if R := (A - B) (mod C) = 0; otherwise if R = 0, C is a Fermat probable prime to base 3^F.
Proof: We compute A - B = 3^(F*C-1) - 3^(F-1) = 3^(F-1) * [ 3^(F*(C-1)-1) - 1 ] (mod N) and then reduce (mod C).
Then we note that R == (A - B) == 0 (mod C)
	iff 3^(F*(C-1)-1) == 1 (mod C), which is just a base-3 Fermat-PRP test to base 3^F. QED
Alternatively, we can get rid of some of the (-1)s by multiplying by the base:
	C is a PRP iff A' - B' := 3^(F*C) - 3^F == 0 (mod C), i.e. if R' := 3^(F*C) - 3^F == 0 (mod N) is divisible bY C.
Example: N = M(109) = 2^109-1 = 745988807.870035986098720987332873; let F = 745988807. Then
	A = 3^(N-1) ==  91603550890398533251179409971537 (mod N)
	B = 3^(F-1) == 610082383855388688949555345767473 (mod N), and we verify that (A - B) == 0 (mod C), thus C is a PRP.
Similarly, if we let A' = 3^N == 3.A (mod N) and B' = 3^F == 3.B (mod N), that also satisfies (A' - B') == 0 (mod C).
*/
uint32 Suyama_CF_PRP(uint64 p, uint64*Res64, uint32 nfac, double a[], double b[], uint64 ci[], uint32 ilo,
	int	(*func_mod_square)(double [], int [], int, int, int, uint64, uint64, int, double *, int, double *),
	int n, int scrnFlag, double *tdiff, char*const gcd_str)
{
	uint64 *ai = (uint64 *)a, *bi = (uint64 *)b;	// Handy 'precast' pointers
	uint32 i,j,k, isprime, ierr = 0, ihi, fbits,lenf;
	uint32 kblocks = n>>10, npad = n + ( (n >> DAT_BITS) << PAD_BITS );	// npad = length of padded data array
	uint64 itmp64, Res35m1, Res36m1;	// Res64 from original PRP passed in via pointer; these are locally-def'd
	cbuf[0] = '\0';
	snprintf(cbuf,STR_MAX_LEN*2,"Suyama-PRP on cofactors of %s: using FFT length %uK = %u 8-byte floats.\n",PSTRING,kblocks,n);//	strcat(cbuf,cstr);
	// sprintf(cstr, " this gives an average %20.15f bits per digit\n",1.0*p/n);	strcat(cbuf,cstr);
	mlucas_fprint(cbuf,1);
	// Pepin-test output = P, vs Mersenne-PRP (type 1) residue = A; thus only need an initial mod-squaring for:
	// the former. Compute Fermat-PRP residue [A] from Euler-PRP (= Pepin-test) residue via a single mod-squaring:
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
		ASSERT(ilo == p-1, "Fermat-mod cofactor-PRP test requires p-1 mod-squarings!");
		snprintf(cbuf,STR_MAX_LEN*2,"Doing one mod-%s squaring of iteration-%u residue [Res64 = %016llX] to get Fermat-PRP residue\n",PSTRING,ilo,*Res64);
		mlucas_fprint(cbuf,1);
		ilo = 0;	ihi = ilo+1;	// Have checked that savefile residue is for a complete PRP test, so reset iteration counter
		BASE_MULTIPLIER_BITS[0] = 0ull;
/*A*/	ierr = func_mod_square(a, (int*)ci, n, ilo,ihi, 0ull, p, scrnFlag, tdiff, TRUE, 0x0);
		convert_res_FP_bytewise(a, (uint8*)ci, n, p, Res64, &Res35m1, &Res36m1);	// Overwrite passed-in Pepin-Res64 with Fermat-PRP one
		snprintf(cbuf,STR_MAX_LEN,"MaxErr = %10.9f\n",MME); mlucas_fprint(cbuf,1);
	} else if (MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {	// Mersenne PRP-CF doesn't have the Res35m1 or Res36m1 values passed in,
		res_SH(ci,n,&itmp64,&Res35m1,&Res36m1);			// so we refresh these; see https://github.com/primesearch/Mlucas/issues/27
	}
	if(ierr) {
		snprintf(cbuf,STR_MAX_LEN*2,"Error of type[%u] = %s in mod-squaring ... aborting\n",ierr,returnMlucasErrCode(ierr));
		mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
	}
	sprintf(cbuf, "Fermat-PRP residue (A)     = 0x%016llX,%11llu,%11llu\n",*Res64,Res35m1,Res36m1);
	mlucas_fprint(cbuf,1);
	j = (p+63)>>6;	// j = uint64 vector length; Omit leading '1' bit in Fermat case since PRP-residue only has that set if a Fermat prime
	mi64_set_eq(ai,ci,j);	// Copy packed-bit result back into low ceiling(p/64) bytes of A-vec (treated as a uint64 array)
	// Compute "prime-factor product residue" [B] from Euler-PRP (= Pepin-test) residue ... first init bitwise mul-by-base array = F, i.e. storing product of known small-prime factors:
	if(!nfac) {
		sprintf(cbuf, "Cofactor-PRP test requires one or more known factors!");
		mlucas_fprint(cbuf,0); ASSERT(0, cbuf);
	}
	BASE_MULTIPLIER_BITS[0] = 1ull;	lenf = 1;
	// Multiply each known-factor with current partial product of factors.
	// Use BASE_MULTIPLIER_BITS to store factor product here, but need curr_fac[] for intermediate partial products:
	uint64 curr_fac[20];
	for(i = 0; KNOWN_FACTORS[i] != 0ull; i += 4) {
		k = mi64_getlen(KNOWN_FACTORS+i,4);	// k = number of nonzero limbs in curr_fac (alloc 4 limbs per in KNOWN_FACTORS[])
		// Multiply factor into current partial product of factors; use curr_fac[] array to store product to work around none-of-3-input-pointers-may-coincide restriction in mi64_mul_vector:
		mi64_mul_vector(BASE_MULTIPLIER_BITS,lenf, KNOWN_FACTORS+i,k, curr_fac,&lenf);
		mi64_set_eq(BASE_MULTIPLIER_BITS,curr_fac,lenf);
	}
	ASSERT((i>>2) == nfac, "Number of known-factors mismatch!");
	ASSERT(lenf <= 20, "Product of known-factors too large to fit into curr_fac[]!");
	for(i = 0; i < lenf; i++) { curr_fac[i] = 0ull; }	// Re-zero the elts of curr_fac[] used as tmps in above loop
	fbits = (lenf<<6) - mi64_leadz(BASE_MULTIPLIER_BITS, lenf);
	// Now that have F stored in BASE_MULTIPLIER_BITS array, do powmod to get B = base^(F-1) (mod N):
	BASE_MULTIPLIER_BITS[0] -= 1ull;	// F-1; no chance of a borrow here
	for(i = 0; i < npad; i++) { b[i] = 0; }	// Zero the elements of the floating-point array b[]
	/****** Note: For Fermat *cofactor* PRP check we use a PRP assignment (not Pepin-test, though we need that residue as our
	input), meaning that PRP_BASE = 3, not the speecial value 2 it has for residue-shift purposes in Pepin test mode: ******/
	b[0] = PRP_BASE;	ASSERT(PRP_BASE < (1 << (uint32)ceil(1.0*p/n)), "PRP_BASE out of range!");
	ilo = 0;	ihi = fbits-1;	// LR modpow; init b[0] = PRP_BASE takes cares of leftmots bit
	RES_SHIFT = 0ull;	// Zero the residue-shift so as to not have to play games with where-to-inject-the-initial-seed
	mi64_brev(BASE_MULTIPLIER_BITS,ihi);	// bit-reverse low [ihi] bits of BASE_MULTIPLIER_BITS:
/*B*/	ierr = func_mod_square(b, (int*)ci, n, ilo,ihi, 0ull, p, scrnFlag, tdiff, TRUE, 0x0);
	if(ierr) {
		snprintf(cbuf,STR_MAX_LEN*2,"Error of type[%u] = %s on iteration %u of mod-squaring chain ... aborting\n",ierr,returnMlucasErrCode(ierr),ROE_ITER);
		mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
	}
	sprintf(cbuf,"Processed %u bits in binary modpow; MaxErr = %10.9f\n",ihi,MME);
	convert_res_FP_bytewise(b, (uint8*)ci, n, p, &itmp64, &Res35m1, &Res36m1);	// Res64 reserved for Fermat-PRP result; use itmp64 here
	sprintf(cstr, "%u^(F-1) residue (B)        = 0x%016llX,%11llu,%11llu\n",PRP_BASE,itmp64,Res35m1,Res36m1);
	strcat(cbuf,cstr);	mlucas_fprint(cbuf,1);
	ASSERT(j = (p+63)>>6,"uint64 vector length got clobbered!");
	mi64_set_eq(bi,ci,j);	// Copy packed-bit result into low j limbs of B-vec (treated as a uint64 array)
	itmp64 = mi64_sub(ai,bi, ai,j);
	// If result < 0, need to add Modulus - for N = Fm,Mp this means +-1 in LSW, respectively.
	// For Fermat case, the borrow out of the high limb in the preceding vector-sub is canceled by the
	// leading binary '1' in F[m]; in the Mersenne case, need to explicitly add 2^(p%64) to high limb:
	if(itmp64) {
		ASSERT(itmp64 == 1ull,"Carryout = 1 expected!");
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			itmp64 = mi64_sub_scalar(ai,1ull, ai,j);
			ai[j-1] += 1ull << (p&63);
		} else {
			itmp64 = mi64_add_scalar(ai,1ull, ai,j);
		}	ASSERT(itmp64 == 0ull,"Carryout = 0 expected!");
	}
	// B-array again free, re-use in uint64-cast form to compute C = Fm/F and (A-B) mod C:
	// Compute Modulus ... note mi64-vecs have no cache-oriented element padding:
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		itmp64 = -1ull;
		// Loop rather than call to mi64_set_eq_scalar here, since need to set all elts = -1:
		for(i = 0; i < j; i++) { bi[i] = itmp64; }
		bi[j-1] >>= 64-(p&63);	// Leading word needs >> to leave just low p%64 bits set
	} else {
		// j = uint64 vector length; init sans the leading '1' word, then increment prior to mi64_div
		mi64_clear(bi,j);
		bi[j++] = bi[0] = 1ull;	// Post-increment here
	}
// (A - B) in ai[], F in BASE_MULTIPLIER_BITS[], C = N/F returned in ci[]:
	mi64_brev(BASE_MULTIPLIER_BITS,ihi);// 2nd BR of low [ihi] bits of BASE_MULTIPLIER_BITS to recover the factored part F-1, sans leftmost bit...
	BASE_MULTIPLIER_BITS[lenf-1] += 1ull << (fbits-1);	// Restore leftmost bit ...
	BASE_MULTIPLIER_BITS[     0] += 1ull;	// ... and add 1 to recover F; no chance of a carryout here
	// Since F << N, use Mont-mul-div for C - quotient overwrites N, no rem-vec needed, just verify that F is in fact a divisor:
	ASSERT(1 == mi64_div(bi,BASE_MULTIPLIER_BITS, j,lenf, ci,0x0), "C = N/F should have 0 remainder!");	// C in ci[]
	j -= (MODULUS_TYPE == MODULUS_TYPE_FERMAT);	// In Fermat case, undo the above j++ used to insert the leading bit in F[m]
	i = j;	j = mi64_getlen(ci, j);	// *** Apr 2022 bug: don't add extra limb for Fermat-case to i here since (A-B) < N ***
// R = (A - B) mod C in B-array (bi[]); store Q = (A - B)/C in curr_fac[] in case want to remultiply and verify Q*C + R = (A - B):
	sprintf(cbuf,"(A - B) Res64 = 0x%016llX, C Res64 = 0x%016llX\n",ai[0],ci[0]);
	mlucas_fprint(cbuf,1);
	mi64_div_binary(ai,ci, i,j, curr_fac,(uint32 *)&k, bi);	// On return, k has quotient length; curr_fac[] = quo, bi[] = rem
	snprintf(cbuf,STR_MAX_LEN*2,"(A - B)/C: Quotient = %s, Remainder Res64 = 0x%016llX\n",&cstr[convert_mi64_base10_char(cstr,curr_fac,k,0)],bi[0]);
	mlucas_fprint(cbuf,1);
	// For 1-word quotient q, double-check binary-div result by computing (q*denominator + r) and comparing vs numerator:
  #if 0	/*** May 2022: This overwrites ci[], which hoses the is-cofactor-a-prime-power GCD() below ***/
	if(k == 1) {
		ASSERT(0 == mi64_mul_scalar_add_vec2(ci, curr_fac[0], bi, ci, i), "Unexpected carryout!");
		ASSERT(1 == mi64_cmp_eq(ai,ci,i), "Q*C + R = (A - B) check fails!");
	}
  #endif
	snprintf(cbuf,STR_MAX_LEN*2,"Suyama Cofactor-PRP test of %s",PSTRING);
	// Base-2 log of cofactor = lg(Fm/F) = lg(Fm) - lg(F) ~= 2^m - lg(F). 2^m stored in p, sub lg(F) in loop below:
	double lg_cof = p,lg_fac,log10_2 = 0.30102999566398119521;	// Use lg_fac to store log2 of each factor as we recompute it
	for(i = 0; KNOWN_FACTORS[i] != 0ull; i += 4) {
		k = mi64_getlen(KNOWN_FACTORS+i,4);	// k = number of nonzero limbs in curr_fac (alloc 4 limbs per in KNOWN_FACTORS[])
		strcat( cbuf, " / " );
		strcat( cbuf, &cstr[convert_mi64_base10_char(cstr, KNOWN_FACTORS+i, k, 0)] );
		lg_fac  = (double)mi64_extract_lead64(KNOWN_FACTORS+i, k, &itmp64) - 64;
		lg_fac += log((double)itmp64)*ILG2;
		lg_cof -= lg_fac;
	}
	i = ceil(lg_cof*log10_2);	// #decimal digits of cofactor
	j = mi64_getlen(bi,j);	// Returns 0 iff all limbs of remainder == 0
	isprime = !j;
	if(!j) {
		sprintf(cbuf,"This cofactor is PROBABLE PRIME [PRP%u].\n",i);	mlucas_fprint(cbuf,1);
	} else {
		res_SH(bi,j,&itmp64,&Res35m1,&Res36m1);	// Res64 reserved for Fermat-PRP result; use itmp64 here
		sprintf(cstr," with FFT length %u = %u K:\n\t(A - B) mod C has Res64,35m1,36m1: 0x%016llX,%11llu,%11llu.\n",n,kblocks,itmp64,Res35m1,Res36m1);
		strcat(cbuf,cstr);	mlucas_fprint(cbuf,1);
		/* Compute gcd(A - B,C) [cf. Phil Moore post: https://mersenneforum.org/showpost.php?p=210599&postcount=67]
			"Take the GCD of the difference of these two residues (A - B) with C. If the GCD is equal to 1,
			C cannot be a prime power. (If it is not equal to 1, we have discovered a new factor of C.)"
		*/
		sprintf(cbuf,"This cofactor is COMPOSITE [C%u]. Checking prime-power-ness via GCD(A - B,C) ... \n",i); mlucas_fprint(cbuf,1);
		i = gcd(0,0ull,ai,ci,j,gcd_str);	// 1st arg = stage of (p-1 or ecm) just completed, does not apply here
		if(i)
			sprintf(cbuf,"Cofactor is a prime power! GCD(A - B,C) = %s.\n",gcd_str);
		else
			sprintf(cbuf,"Cofactor is not a prime power.\n");
		mlucas_fprint(cbuf,1);
	}
	return isprime;
}

/***********************************************************************************************/
/****** Thanks to Tom Cage (RIP) for the initial version of the tuning-self-test harness: ******/
/***********************************************************************************************/
#ifdef macintosh
	#include <console.h>	/* Macintosh */
#endif

/* Number of distinct FFT lengths supported for self-tests: */
#define numTest				136	// = sum of all the subranges below
/* Number of FFT lengths in the various subranges of the full self-test suite: */
#define numTeensy			15	// v21: added 'Teensy', moved 8 smallest 'Tiny' into it,
#define numTiny 			32	// changed counts from [-,32,32,16,24,16,9,0,0] to [15,32,24,20,20,16,9,0,0]
#define numSmall			24
#define numMedium			20
#define numLarge			20
#define numHuge				16
/* Adding larger FFT lengths to test vectors requires supporting changes to Mdata.h:MAX_FFT_LENGTH_IN_K and get_fft_radices.c */
#define numEgregious		 9
#define numBrobdingnagian	 0
#define numGodzillian		 0
#if(numTeensy + numTiny + numSmall + numMedium + numLarge + numHuge + numEgregious + numBrobdingnagian + numGodzillian != numTest)
	#error Sum(numTeensy + ...) != numTest in main!
#endif

struct res_triplet{
	uint64 sh0;	/* Res64 (assuming initial LL test seed = 4) for this exponent. */
	uint64 sh1;	/* Residue mod 2^35-1 */
	uint64 sh2;	/* Residue mod 2^36-1 */
};

/* Array of distinct test cases for Mersenne self-tests. Add one extra slot to vector for user-specified self-test exponents;
We use p's given by given_N_get_maxP(), so maximum RO errors should be consistent and around the target 0.25 value, a little
bit higher in SSE2 mode.

The first letter of the descriptor for each size range serves as a menmonic for the -[*] option which runs the self-tests
for that range, e.g. -s runs the [Small] range, -e the [Egregious], etc.

HERE IS THE PROCEDURE FOR ADDING A BEW ENTRY TO THE EXPONENT/RESIDUE TABLE BELOW:

1. Increment the num**** entry above corr. to the self-test subset which is being augmented; also ++numTest;

2. Go to get_fft_radices.c and use the PARI code in the commentary immediately below given_N_get_maxP()
	to determine maxP for the new runlength;

3. Use PARI isprime() to find the largest prime <= maxP;

4. Run 100 and 1000-iteration self-tests at the next-higher runlength already appearing in the self-test table;

5. Use the results - specifically the hexadecimal Res64 and the mod 2^35-1 and mod 2^36-1 SH residues - for the
	2 runs to create a table entry for the new runlength;

6. Rebuild this file, re-link and repeat the 2 self-tests on the exponent in the new table row, at the new runlength.
*/
struct testMers{
	int fftLength;		/* FFT length in K (i.e. 4 means an array of 4K doubles) */
	uint64 exponent;	/* Test exponent for this FFT length. */
	struct res_triplet	res_t[3];	/* 100,1000 and 10000-iteration SH-residue triplets */
};

// Reference LL-test residues:
struct testMers MersVec[numTest+1] =
{
/*                                         100-iteration residues:	                               1000-iteration residues:                */
/*	  FFTlen(K)     p              Res64           mod 2^35-1      mod 2^36-1               Res64           mod 2^35-1      mod 2^36-1     */
/*	    -----    --------     ----------------     -----------     -----------         ----------------     -----------     -----------    */
	/* Teensy:                                     [%34359738367  ][%68719476735  ]                         [%34359738367  ][%68719476735  ] */
	{     1,     22679ull, { {0x27F79B322B31FDA9ull, 23663495240ull,  2817376218ull}, {0xF602AA618A821226ull,  1646512433ull, 39303046116ull}, {0x31E33974FE4AB5F8ull, 26644986014ull, 46355958559ull} } },
	{     2,     44657ull, { {0xC15789922E67B055ull,  1083638452ull, 32161054300ull}, {0xAFA7CCAFEEE4A0AFull, 20332793592ull, 55476105191ull}, {0xE1061D8AAAA1E4E4ull, 17639694995ull, 28777282852ull} } },
	{     3,     66431ull, { {0xC17597154ED43DF7ull, 18774471399ull, 20870587619ull}, {0x1EC9A4159C4067BCull, 16780639455ull,  8129191272ull}, {0x6075773CD190248Dull, 20838041038ull, 20210614176ull} } },
	{     4,     88019ull, { {0x7AFD2F95BF1D2361ull, 11148380131ull, 43485107986ull}, {0x87FE95F54CC56993ull,   802983452ull, 49338968618ull}, {0x14F089BA4C45AFECull, 10078158375ull, 43611889204ull} } },
	{     5,    109481ull, { {0x43530D2B6398F339ull, 34170186731ull, 26358746864ull}, {0x9327D1EB998E1F9Full, 33838644122ull, 34880982988ull}, {0x5A2E80997F65C230ull,  9603989727ull, 12982729662ull} } },
	{     6,    130873ull, { {0x1ECD6D4A5257DF87ull, 30276540713ull, 35782617858ull}, {0x21DDB2C246D2712Cull, 31677605153ull, 20356415217ull}, {0x506CB5100C615C81ull,  3178546981ull, 32034984726ull} } },
	{     7,    152197ull, { {0xCD9D42BD1FC8898Aull,  6862900694ull, 32580164643ull}, {0x61AF4697185D79CFull,  9022613594ull, 37933406056ull}, {0xA4D18D8B27B0CC74ull, 15251090953ull, 41762157108ull} } },
	{     8,    173431ull, { {0x85301536E4CA9B11ull,  3707224323ull, 36851834664ull}, {0x2FD5120BEC41F449ull, 28734955954ull, 23103392631ull}, {0x139D1D396F173696ull, 12716541843ull, 58117202214ull} } },
	{     9,    194609ull, { {0xC711AF1008612BC6ull,  1574019740ull, 37260026270ull}, {0x5153F6E040CD1BE6ull, 15446410924ull,  3291404673ull}, {0x33E19077F35070A3ull, 34231248547ull, 24411805292ull} } },
	{    10,    215767ull, { {0x4428783BC62760F0ull,  7466284975ull, 53916123655ull}, {0xED46A8C001908815ull,   739143119ull, 36950829937ull}, {0xCBE0AD544E96FDB9ull,  7625169722ull, 52027104104ull} } },
	{    11,    236813ull, { {0x592D849AF4D1336Full, 29025996994ull, 48905971124ull}, {0xB4EEB63BB656F424ull,  5361680901ull, 31850818767ull}, {0x5AEBAE493B085903ull, 20240457506ull, 42866015746ull} } },
	{    12,    257903ull, { {0x1D8121DE28B60996ull, 22402402961ull, 65583959369ull}, {0x54F2BE961A674CB1ull, 25601315553ull, 54298899520ull}, {0x50DFBA28D0D1A8C3ull, 17306049864ull, 68068537809ull} } },
	{    13,    278917ull, { {0xE3BC90B0E652C7C0ull, 21244206101ull, 51449948145ull}, {0x93AF8994F95F2E50ull, 16427368469ull, 10707190710ull}, {0x1674AAA04F7BD61Aull, 12079507298ull, 56593045102ull} } },
	{    14,    299903ull, { {0xDB8E39C67F8CCA0Aull, 20506717562ull, 44874927985ull}, {0x4E7CCB446371C470ull, 34135369163ull, 61575700812ull}, {0x04ACC83FFE9CEAD4ull, 26179715264ull, 65445483729ull} } },
	{    15,    320851ull, { {0xB3C5A1C03E26BB17ull, 22101045153ull,  4420560161ull}, {0x923A9870D65BC73Dull, 29411268414ull, 30739991617ull}, {0xB3F1ACF3A26C4D72ull, 32179253815ull, 68615042306ull} } },
	/* Tiny: */
	{    16,    341749ull, { {0x8223DF939E46A0FFull, 32377771756ull, 38218252095ull}, {0xC6A5D4B6034A34B8ull, 31917858141ull, 59888258577ull}, {0x93EF44581866E318ull, 18805111197ull,  8333640393ull} } },
	{    18,    383521ull, { {0xBF30D4AF5ADF87C8ull, 15059093425ull, 52618040649ull}, {0x9F453732B3FE3C04ull,  4385160151ull, 47987324636ull}, {0x0DBF50D7F2142148ull,  1608653720ull, 52016825449ull} } },
	{    20,    425149ull, { {0x6951388C3B99EEC0ull,  4401287495ull, 19242775142ull}, {0x501CEC2CB2080627ull, 21816565170ull, 41043945930ull}, {0x5A9A9BF4608090A2ull, 27025233803ull, 68581005187ull} } },
	{    22,    466733ull, { {0xD95F8EC0F32B4756ull, 19305723506ull, 26588871256ull}, {0xB1F58184918D94B6ull,  8443388060ull, 11738516313ull}, {0xAC4B1F499BF2C2DAull,  7322105347ull, 15747709958ull} } },
	{    24,    508223ull, { {0xDA46E41316F8BCCAull, 25471180026ull,  1635203275ull}, {0x27A5B285281466B9ull, 11438869313ull,  7226774009ull}, {0x4ABED2868B800F7Dull,  7783533092ull, 66921500486ull} } },
	{    26,    549623ull, { {0x6649D9D6CD4E0CE1ull, 25445908581ull, 26118212198ull}, {0x1A4F280627A15B3Cull, 13286323782ull, 31550278005ull}, {0x86404E236E99B3C4ull, 17401894517ull, 40934891751ull} } },
	{    28,    590963ull, { {0x4ADDB6C4A76465AFull,  6532108269ull, 54921134131ull}, {0x3063D08A7BABD7B8ull,  4777711548ull, 39733274344ull}, {0xBE2ABBB09336F32Eull, 30656127523ull, 50296089656ull} } },
	{    30,    632251ull, { {0x0811FAA40601EB1Dull, 16369365746ull,  6888026123ull}, {0xF324E4DEC564AF91ull, 10236920023ull, 34068699974ull}, {0xAA622CF2A48F6085ull, 22315931502ull,  1049914969ull} } },
	{    32,    673469ull, { {0x1A4EF8A0D172FBAAull, 32667536946ull, 11393278588ull}, {0xA4DFD62B928F68A4ull, 11900420802ull, 66610946021ull}, {0xFA3993AC9CE7BEEDull,   117685830ull, 39361570554ull} } },
	{    36,    755737ull, { {0x13B13C61298088DCull, 34092803628ull,  7584858890ull}, {0x33A2A43DE8782CCCull,  2953124985ull, 62434716987ull}, {0xD0DF76911349551Bull, 28919846011ull, 30432969648ull} } },
	{    40,    837817ull, { {0x88555D9AAD3FF2DDull,  8573348747ull, 67896670216ull}, {0xEAC1676D914878C0ull, 34312095136ull, 45077378164ull}, {0x89E1C4D06BB0F9F3ull,  6272358557ull, 24712951618ull} } },
	{    44,    919729ull, { {0x6ACC03213A37BA5Bull,  3870201113ull, 48145739792ull}, {0xDA98B49CC83C60CBull, 15886769401ull, 62221100895ull}, {0xF1DA20E7D5A89638ull,  4633752262ull, 20941274692ull} } },
	{    48,   1001467ull, { {0x6B1C76AB5431FDA4ull,  6795121241ull, 65308927583ull}, {0xBD99FD21F4136BFCull, 26386988063ull, 61607603549ull}, {0x9136554E4718BFA9ull, 18451197932ull, 37688798842ull} } },
	{    52,   1083077ull, { {0xA591637EC8CF3FE4ull,  4769775755ull, 65114872367ull}, {0xE59C08B13B00E6FFull,  1383963096ull, 26100699764ull}, {0x48CCA2242A1F9352ull, 30318043361ull, 12067176371ull} } },
	{    56,   1164533ull, { {0xEC4F2579E4533584ull,  5456769127ull, 59922459736ull}, {0xF7D2BF94C2767D36ull, 30727892629ull, 48141128220ull}, {0xE332E3891AE98AD7ull,  7024598607ull, 65691841143ull} } },
	{    60,   1245877ull, { {0xC91002E1A4EE7E07ull,  6217476228ull, 40164514288ull}, {0xEABE9E1A31DF5877ull,   831216169ull, 29591771932ull}, {0xCB85101F6857519Dull, 30425930108ull,  2198194326ull} } },
	{    64,   1327099ull, { {0xAC070112281229E0ull, 14226353454ull,  1640524016ull}, {0xF25AA54053C5BB64ull, 32455038659ull, 53547160776ull}, {0x3854D019CE12CC9Aull, 29589836279ull,  2174826233ull} } },
	{    72,   1489223ull, { {0x6674518EA19B3D6Aull, 32383400562ull, 53234746310ull}, {0xEB312091097F6C3Bull,  3980687908ull,  8568698675ull}, {0xC225FAF24A093590ull, 22407813999ull, 30924932017ull} } },
	{    80,   1650959ull, { {0xE5326E754F3202A8ull,  5593233426ull, 33337128557ull}, {0xFC3E8CDA60AF5CF8ull, 11466296968ull, 12651602524ull}, {0x4D68554B73674A60ull,  2253999911ull, 55045374456ull} } },
	{    88,   1812347ull, { {0x81BDD3AC63DF3F73ull, 19957199947ull, 61002681293ull}, {0x3D3E429D7427C4EAull, 25342898119ull, 34322985438ull}, {0x5769D7B47C49436Full,  5234049262ull, 26872574292ull} } },
	{    96,   1973431ull, { {0x901C8305DA9FF95Aull, 32611811878ull, 55986702160ull}, {0x0790CA11ADAA47E3ull, 17075140845ull, 12883521448ull}, {0xF18EADC267DE6FC1ull, 24308841307ull, 31678116890ull} } },
	{   104,   2134201ull, { {0x59BDA0D80F3279EDull, 17901153436ull,  3927067335ull}, {0x2F81B21BC680C861ull, 18443771511ull, 45465079919ull}, {0x439245FA16A38116ull, 20996570088ull,   489289103ull} } },
	{   112,   2294731ull, { {0xC44ACC96D268625Full, 10331638988ull,  2292055445ull}, {0xED20577E16E128DEull, 32248607028ull, 14903460370ull}, {0xCB862A1B42B230A2ull, 23316229090ull, 23891565685ull} } },
	{   120,   2455003ull, { {0xC5F7DB23F174A67Dull, 32991574397ull, 31642856976ull}, {0x401670254012E5ABull, 33626385418ull, 66465546971ull}, {0x20AB396E327C09C1ull, 13309965383ull, 60492105240ull} } },
	{   128,   2614999ull, { {0x040918890E98F8DAull, 14867710211ull, 47602627318ull}, {0x1A184504D2DE2D3Cull,  5934292942ull,  4090378120ull}, {0xE7126F512D3FD742ull, 17101849610ull, 66501661438ull} } },
	{   144,   2934479ull, { {0x1B90A27301980A3Aull,  7043479338ull, 38327130996ull}, {0x8C3045C6534867C6ull, 12456621644ull, 52801948293ull}, {0xF17F4A594A281B94ull,  5970782987ull, 68371435254ull} } },
	{   160,   3253153ull, { {0x9AFD3618C164D1B4ull, 16551334620ull, 55616214582ull}, {0x1493A70897A8D058ull, 34082962858ull, 60773088284ull}, {0x57D3F1A090E78729ull, 26902546905ull, 49396480035ull} } },
	{   176,   3571153ull, { {0xA016F25779902477ull, 21500047857ull,  9150810891ull}, {0x8E6F248EC96445FFull, 22443629034ull, 16625023722ull}, {0xFFD8B840C06A5EACull, 32183737848ull, 42577197579ull} } },
	{   192,   3888509ull, { {0x71E61322CCFB396Cull, 29259839105ull, 50741070790ull}, {0x3CEDB241702D2907ull,  6177258458ull, 21951191321ull}, {0xAF407D11B2D74C3Cull, 31653650180ull, 27299459944ull} } },
	{   208,   4205303ull, { {0xC08562DA75132764ull,  7099101614ull, 36784779697ull}, {0xAD381B4FE91D46FDull,  7173420823ull, 51721175527ull}, {0xC70061EF9537C4E1ull,  9945894076ull,  2301956793ull} } },
	{   224,   4521557ull, { {0xE68210464F96D6A6ull, 20442129364ull, 11338970081ull}, {0x3B06B74F5D4C0E35ull,  7526060994ull, 28782225212ull}, {0xB720ACD1D69A7ECFull, 28103212586ull, 10983125296ull} } },
	{   240,   4837331ull, { {0xB0D0E72B7C87C174ull, 15439682274ull, 46315054895ull}, {0x3AA14E0E90D16317ull,  5730133308ull, 50944347816ull}, {0x12CFBF6001E59FF7ull, 26877054587ull, 60322521357ull} } },
	/* Small: */
	{   256,   5152643ull, { {0x074879D86679CB5Bull,  1208548964ull, 48525653083ull}, {0x98AF5E14C824A252ull,   783196824ull,  6594302302ull}, {0x7DA0D3B9EFEA4931ull, 32608975347ull, 43428286760ull} } },
	{   288,   5782013ull, { {0x9869BE81D9AB1564ull, 15509103769ull, 49640026911ull}, {0x7C998719C6001318ull, 23749848147ull, 19853218689ull}, {0xE2E246D9094EBFD7ull, 26657044660ull,  7091330955ull} } },
	{   320,   6409849ull, { {0x20739E43A693A937ull, 27970131351ull, 15334307151ull}, {0xE20A76DCEB6774A6ull, 14260757089ull, 68560882840ull}, {0xCEC786F8883D8D1Full,  5597853948ull, 57984323163ull} } },
	{   352,   7036339ull, { {0xD6A226BAB5E14D62ull,  9444972171ull, 28639488018ull}, {0x0579D28296F29D92ull, 18964853245ull, 30111685201ull}, {0xB9CF9FE489BD34CFull, 11297283566ull, 16782498229ull} } },
	{   384,   7661567ull, { {0x3A929F577AC9725Full, 23890492835ull, 64821764767ull}, {0x2ECBA785576E6D58ull, 26446200615ull, 60269798452ull}, {0x22AA4A0A7A3676A7ull, 23262227373ull, 26736591016ull} } },
	{   416,   8285659ull, { {0xDCA138D55C36E40Cull, 10452583294ull,  4308453248ull}, {0x1FEE7F79E32229A6ull, 11936075329ull, 16061515794ull}, {0x9B5324C99CCE498Dull, 14697684311ull, 16439873760ull} } },
	{   448,   8908723ull, { {0x436494C7EA194FA1ull,  2976149044ull, 21645125251ull}, {0xD2FCEDE29E818A26ull, 15260150529ull, 11678366985ull}, {0x952F59F4972F830Aull, 17286796157ull, 11216464341ull} } },
	{   480,   9530803ull, { {0x6D52BCCB796A46E9ull, 28296985362ull, 66294636306ull}, {0x786CB610A809B762ull, 24654197494ull, 57943258783ull}, {0x40564770A8540610ull, 27197548273ull, 32971949076ull} } },
	{   512,  10151971ull, { {0x7BA3DD9C38878B83ull, 24395728676ull, 12632037498ull}, {0x8E9FA3093ACD81C1ull,  6345070464ull, 65203567231ull}, {0xBBD3BD10BA9983B5ull,  7822185024ull, 14668688365ull} } },
	{   576,  11391823ull, { {0x78D8769B9F75FB2Bull,  7286184017ull, 17306427758ull}, {0x5834BDA7558DD43Cull, 11189092321ull, 23026923116ull}, {0x793D548FCB76B28Cull, 11310068490ull, 24716315416ull} } },
	{   640,  12628613ull, { {0xF951B697F5C5032Aull,  9262385660ull, 57463106915ull}, {0x93B526040205BA27ull, 30538079080ull, 32317022014ull}, {0x25544FC69ADB28F9ull, 28253742327ull,  1823182110ull} } },
	{   704,  13862759ull, { {0xBB2F69275D79A9EEull, 12230183669ull, 68684647134ull}, {0x7343ECC160AA00D5ull, 24655585170ull, 51102704879ull}, {0x1402EEF49394CDC7ull,  5500341204ull, 59999916295ull} } },
	{   768,  15094403ull, { {0xF6895EB66EADE9C5ull,  8490184692ull, 23393343807ull}, {0xF673A8D6413923A9ull, 20026779905ull, 67516766223ull}, {0xD63752CA13598971ull, 24773095342ull, 29303310893ull} } },
	{   832,  16323773ull, { {0xEB8890F379392B2Full, 27289972116ull, 63975275393ull}, {0xD681EDD3A1EC3780ull, 12515962698ull, 40155157152ull}, {0xEB9C9477368BF584ull, 13378242091ull,  9365072054ull} } },
	{   896,  17551099ull, { {0xAB1180428ED65EE0ull,  3105108668ull, 66518734167ull}, {0x31813367849BBF49ull,  9516734777ull, 18271834608ull}, {0x95C2E1F201FCE598ull,  6264820675ull, 49312303312ull} } },
	{   960,  18776473ull, { {0xCA7D81B22AE24935ull, 24317941565ull, 67706175547ull}, {0x02EB980A49E7B60Full,  5730644436ull, 48386545950ull}, {0xBC6503AA5C062308ull, 29760131532ull, 31603724687ull} } },
	{  1024,  19800083ull, { {0x95AFD7A5269F14F6ull, 26677613826ull, 37493068952ull}, {0x18A53602BAA0E197ull,   624371978ull, 15180896714ull}, {0xF79CD0274644183Dull, 21538070258ull, 26190157173ull} } },
	{  1152,  22217791ull, { {0x9EFE3D8C08D89E48ull,  8307313500ull, 40673806995ull}, {0xE1F94CD14457EDADull, 23322294478ull, 42677160325ull}, {0xF60CFBDEA4ADF55Full,  5469936596ull, 35790203222ull} } },
	{  1280,  24629621ull, { {0x85D92483D90E5029ull,  6310387878ull, 18231127032ull}, {0xBEE63CF182681345ull,  7247112769ull, 24502530130ull}, {0x89977592EE68F853ull, 31511858957ull, 44154237710ull} } },
	{  1408,  27036157ull, { {0x3BAC320B542307B0ull, 28733834194ull, 35177623244ull}, {0xBAE66098A85CD960ull, 25462246119ull, 52444387524ull}, {0x0569BF7834267F10ull, 27834108530ull, 44136002879ull} } },
	{  1536,  29437799ull, { {0x9E279D6E3752A61Cull,  6467723895ull, 50774708769ull}, {0x5A2D06C05D222BF2ull, 22780146874ull, 60521558933ull}, {0xC6D352199D949DCDull, 30008450771ull, 37092077102ull} } },
	{  1664,  31835017ull, { {0xE10EFEAEDCF46110ull, 21648491345ull, 41207849648ull}, {0x8717F387BB55E1ECull, 17021700047ull, 59165499154ull}, {0x5CDA924B80871209ull, 22299652758ull, 41171353038ull} } },
	{  1792,  34228133ull, { {0x9FC8394655E0334Eull,  1603847275ull, 51947401644ull}, {0x1128E4676929F8C8ull, 33342766910ull, 55912489998ull}, {0xE697B18729853BEEull, 25335158858ull, 63783469524ull} } },
	{  1920,  36617407ull, { {0xB7FA68741ABA807Aull,  2831791183ull, 44591522415ull}, {0x0F635E0B2DC95F81ull, 24912661453ull, 45929081959ull}, {0x2B5BCFC6BA94177Aull,  7824653065ull, 62951001255ull} } },
	/* Medium: */
	{  2048,  39003229ull, { {0xEC810981F56D5EC7ull, 29671814311ull, 35851198865ull}, {0xC7979AEE894F6DDEull,  6017514065ull, 41670805527ull}, {0x4669B1DD352C46BCull,  4806501770ull, 33957162576ull} } },
	{  2304,  43765019ull, { {0x672821643AEE9552ull, 19240824825ull,  8017358641ull}, {0x90CF100A46BE5B64ull, 24299840772ull, 29209909852ull}, {0xA280EC055C5ADAE2ull, 28830689436ull, 48353443940ull} } },
	{  2560,  48515021ull, { {0x015D3A5DC74024D1ull,  4691430475ull, 36435156228ull}, {0xFCFC9219B830DA28ull, 11946340838ull, 45970544027ull}, {0x88EBED279C6A1DE3ull,  9620966362ull,  3964242016ull} } },
	{  2816,  53254447ull, { {0x75F9FC3C84FF5B40ull,   745331612ull, 37955073618ull}, {0xBECD92EA2134D5FEull, 18833421072ull, 15730516149ull}, {0x57B3BAD33AFBC633ull, 26171947880ull, 67540452305ull} } },
	{  3072,  57984131ull, { {0xC7632FA366FCA848ull, 19810840296ull,  5919463661ull}, {0x678F955DDBA52AA4ull, 11995136918ull, 51932258375ull}, {0xE24A37C7B3884058ull, 17264313237ull,  9074782273ull} } },
	{  3328,  62705077ull, { {0xE2CC67B70119ECA9ull,  7229868717ull, 32316011766ull}, {0xD2F94003F0B47BBCull, 19281845997ull, 52908700180ull}, {0xAF81565673A4D22Bull, 11464285328ull, 30793044388ull} } },
	{  3584,  67417873ull, { {0x536142387279C6B7ull, 16480189705ull, 28663441842ull}, {0x4398ADDF46CE1F19ull, 20678941264ull, 64584283338ull}, {0xB88E95C5FA3D4C5Bull,  3311195148ull,  4541514528ull} } },
	{  3840,  72123137ull, { {0x99D05C8A80D2C4ABull, 23223478983ull,  7287108796ull}, {0x48F0329C29D5464Bull,  7693484624ull, 34460970218ull}, {0xA1EA526351C5A64Eull, 10198942223ull, 29783198404ull} } },
	{  4096,  76821337ull, { {0x18F108A0AEC92DA2ull, 25631216955ull, 19217786538ull}, {0x91B09D853C38FC2Bull,  1172724107ull, 62760945815ull}, {0xFB663A86824AAFC5ull,  9621633845ull, 68132687365ull} } },
	{  4608,  86198291ull, { {0x238B15CD7C7C8936ull, 34119649670ull, 50469517589ull}, {0xF1A033E0DEC1330Eull, 16165796253ull, 41865067374ull}, {0x2F1992A097857555ull,  7327358291ull, 32105027991ull} } },
	{  5120,  95551873ull, { {0x4AA89191484E7B56ull, 23537039943ull, 42089192670ull}, {0x811090D60FA9522Cull,  9010676320ull, 44606915786ull}, {0x7BE331521D4C32C2ull, 29429879592ull, 10276666456ull} } },
	{  5632, 104884309ull, { {0x1B2A391D6B7404CDull,  7884768330ull, 55887740575ull}, {0x77CA1E7005AFA4E0ull, 26306026902ull, 21589856845ull}, {0x0F08D303912B27B4ull,   630913159ull, 40892070508ull} } },
	{  6144, 114197579ull, { {0x6903363805E74E29ull,  6950779993ull, 63247095474ull}, {0xEB9D9802854B00A5ull, 24150112753ull, 67039249183ull}, {0xB64714689FA05366ull, 28265310306ull, 41895437602ull} } },
	{  6656, 123493333ull, { {0x651AA3D64C84FA54ull,  9429656635ull, 48354702979ull}, {0x967C90E697CCE6D9ull,  6779392734ull, 18484099736ull}, {0x6C9511DD6EA12528ull, 27647756232ull, 21104526614ull} } },
	{  7168, 132772789ull, { {0xDD02AEFE839F92D5ull,  7411321303ull, 16339659737ull}, {0xED6E26868AC2833Eull, 14154101692ull, 46327957293ull}, {0x80610E8FC3EB92E2ull, 19290762572ull, 46994666267ull} } },
	{  7680, 142037359ull, { {0x9CD0C494D16CB432ull, 23235865558ull, 14066262122ull}, {0xFD6240B21A370394ull, 16216979592ull, 44514519060ull}, {0x097C240EA1436743ull,  5457504643ull, 58797441684ull} } },
	{  8192, 152816047ull, { {0xB58E6FA510DC5049ull, 27285140530ull, 16378703918ull}, {0x75F2841AEBE29216ull, 13527336804ull,   503424366ull}, {0x99F8960CD890E06Aull,  8967321988ull, 43646415661ull} } },
	{  9216, 171465013ull, { {0x60FE24EF89D6140Eull, 25324379967ull,  3841674711ull}, {0x6753411471AD8945ull, 17806860702ull,  3977771754ull}, {0xED3635BF88F37FEFull,  7478721112ull, 47452797377ull} } },
	{ 10240, 190066777ull, { {0x65CF47927C02AC8Eull, 33635344843ull, 67530958158ull}, {0xBADA7FD24D959D21ull, 12777066809ull, 67273129313ull}, {0x82F65495D24A985Full, 22254800275ull, 49183722280ull} } },
	{ 11264, 208626181ull, { {0x6FC0151B81E5173Full, 29164620640ull, 19126254587ull}, {0xD74AA66757A5345Eull, 17524190590ull, 14029371481ull}, {0xDCF9ED39C7EB15B8ull, 34266921309ull, 65896285387ull} } },
	/* Large: */
	{ 12288, 227147083ull, { {0xE01AE9C859ADB03Aull,  7273133358ull,   681418986ull}, {0x303F142E1E88D5B4ull, 28479237457ull, 42044197589ull}, {0x3102781BC131D263ull, 24437355640ull, 48518577431ull} } },
	{ 13312, 245632679ull, { {0x0A6ACB405ADC0354ull,    39452330ull, 38999048555ull}, {0xB38B02A4F195762Full,  3280152282ull, 30314100936ull}, {0xF020F5041AE2CABEull, 24388185991ull, 16285954298ull} } },
	{ 14336, 264085733ull, { {0x5ACE4CCE3B925A81ull,  4584210608ull, 36618317213ull}, {0x02F5EC0CBB1C2032ull, 27165893636ull,   687123146ull}, {0xC6D65BD8A6087F08ull, 15586314376ull, 54717373852ull} } },
	{ 15360, 282508657ull, { {0xE7B08ED3A92EC6ECull,   875689313ull, 41754616020ull}, {0xD08FBAFF5CA5096Full, 30398073011ull, 62088094181ull}, {0xD6B7357DF761AA51ull, 28631146088ull, 26883666300ull} } },
	{ 16384, 300903377ull, { {0xA23E8D2F532F05E6ull, 17871262795ull, 53388776441ull}, {0x14F20059083BF452ull, 16549596802ull, 56184170215ull}, {0x76B8A857EC9B3042ull, 14094306048ull, 61845793513ull} } },
	{ 18432, 337615277ull, { {0xAEB976D153A4176Bull, 15040345558ull, 14542578090ull}, {0x503B443CB1E0CD2Dull, 29149628739ull,  5785599363ull}, {0x2D3047CEFF2F5A6Dull,  6100949709ull, 36303747216ull} } },
	{ 20480, 374233309ull, { {0x6D95C0E62C8F9606ull,  3426866174ull, 39787406588ull}, {0xD08FB9031D460B7Eull, 30083048700ull, 30636357797ull}, {0x3A58018C387FBB68ull, 26771468430ull,  7763681227ull} } },
	{ 22528, 410766953ull, { {0x254572CAB2014E6Cull, 17794160487ull, 13396317974ull}, {0x6202B11AA8602531ull,  9896689730ull, 13179771096ull}, {0x57919ED698DB2058ull, 12202175362ull, 54676834262ull} } },
	{ 24576, 447223969ull, { {0x547DF462C7DAD1F6ull, 21523243373ull, 60663902864ull}, {0x06AB4D8E6FD9D69Bull, 23751736171ull, 10720018557ull}, {0x2B0DE10480A159B7ull,  8420637297ull, 56684836957ull} } },
	{ 26624, 483610763ull, { {0xED3E248A29A1C6A8ull,  3863520901ull, 56560420765ull}, {0xC29358F8206746D6ull, 28829535828ull,  8160695393ull}, {0x56442E62439686ABull, 16425477242ull, 62275447263ull} } },
	{ 28672, 519932827ull, { {0xCA7B3A76819D67F7ull,  7078504016ull, 32836389262ull}, {0x5799C8BE8E02B56Full, 22269194969ull, 11617462155ull}, {0xBA8FC230F7B5EA7Dull, 27830917747ull, 28996845727ull} } },
	{ 30720, 556194803ull, { {0x56CDF80EFF67C1C1ull, 15611746359ull, 45837150207ull}, {0xCBC88456B4B47AC0ull,  3535843471ull, 18652930008ull}, {0xB0B6B164FC22EA57ull, 23758922431ull, 63785512003ull} } },
	{ 32768, 592400713ull, { {0xD8C829884B234EB4ull, 13278383896ull, 54835531046ull}, {0x994DD6B24F452451ull, 28289805997ull, 11462134131ull}, {0x8EEAD14850955F52ull,  3931092242ull,  2483613485ull} } },
	{ 36864, 664658101ull, { {0x2A809B0C735BAC4Bull, 10513414423ull, 54347266174ull}, {0xAB2147D9BAA22BB4ull, 12259954326ull, 67125404781ull}, {0xB4CFF625B3E3FD79ull, 11392807713ull, 32757679957ull} } },
	{ 40960, 736728527ull, { {0xB9AC3EC848FF60A5ull,  7352037613ull,  7261166574ull}, {0x3D623A79D0F14EFFull, 31246254654ull, 49195074754ull}, {0x08F1C4771CDDC601ull, 26432814693ull, 42011833744ull} } },
	{ 45056, 808631029ull, { {0x9D543D67F48AF766ull,  4288733345ull, 27338399854ull}, {0x62A4DF80612E897Bull, 32232536296ull, 47296762118ull}, {0xC42E0660237238BFull,  2574282651ull, 59023604757ull} } },
	{ 49152, 880380937ull, { {0x82ED59E22D585BF6ull, 34028441473ull, 54282489944ull}, {0x7F110FD687DB7CB5ull, 14440703450ull, 57885720403ull}, {0xFB7D2A3EA41594FEull,  2563979873ull, 33101322377ull} } },
	{ 53248, 951990961ull, { {0xE98B9D69E5F350D1ull, 14692034938ull,  2987127112ull}, {0x634157BCC596287Dull, 31799414844ull, 64590776355ull}, {0x2D82590AC2DCB435ull, 25060823443ull, 13978506048ull} } },
	{ 57344,1023472049ull, { {0x960CE3D7029BDB70ull,  2075307892ull, 59259408155ull}, {0x9D98FF9A3FD21D1Bull, 10507186734ull,  3891581073ull}, {0x1CAE6ACE4720BCE8ull, 34199730716ull, 12202402908ull} } },
	{ 61440,1094833457ull, { {0x9A2592E96BC1C827ull, 18678252759ull,   949397216ull}, {0x79BF2E2F5AE7985Bull,  8407199527ull, 64114889744ull}, {0x911DB4B5D8EFD861ull,  1735549534ull, 50988019040ull} } },
	/* Huge: */
	{ 65536,1154422469ull, { {0x192E6EAFFD43E9FAull, 30072967508ull, 30146559703ull}, {0xF2F00B7E84C2E1BFull, 19943239289ull, 24420286066ull}, {0xAB8AC7C62F338C6Bull, 12524447086ull, 59476256925ull} } },
	{ 73728,1295192531ull, { {0xB27DA804065CE0F0ull, 32028474870ull,  4673589134ull}, {0xC27FFBC8C26ADECAull, 14821731491ull, 12633158186ull}, {0xB59F4E88ED5BD206ull,  8793702325ull, 23847832170ull} } },
	{ 81920,1435594063ull, { {0xEAA4B1860B83B83Aull,  5360760716ull, 30032711420ull}, {0x8ED1F6797052808Bull, 12465651506ull, 15167335169ull}, {0x11D81C6F0C26E1DDull, 29918698499ull, 14822236568ull} } },
	{ 90112,1575664327ull, { {0xB1208BF90F8A5F2Bull, 14482434352ull, 53904130189ull}, {0x91DCC94B91975220ull, 11918171509ull, 24914196890ull}, {0x14634D2F1AE28E71ull, 27426869580ull, 55674782001ull} } },
	{ 98304,1715433593ull, { {0x888C5077BF7B35C6ull, 11014256284ull,  7058435784ull}, {0xE2900883E7942E82ull,  3249945856ull, 52796723415ull}, {0xF508DA785F13900Aull, 12331159768ull,  8263352568ull} } },
	{106496,1854927187ull, { {0xE214E9D9B5DF88C6ull, 31246037001ull, 41196573397ull}, {0x97B90915FCC84F3Dull,  4066557343ull,  8643832976ull}, {0xC925A7C98B071EA3ull,  4584603566ull,  5494202371ull} } },
	{114688,1994166553ull, { {0x97CB9BB7B10ACDB6ull, 29731218900ull, 41079146164ull}, {0x731E1DEA0C44C6F0ull, 17898260503ull, 50239810521ull}, {0xE3B5A66B1E7BD939ull, 18336439338ull,  8747708778ull} } },
	{122880,2133169847ull, { {0xE47D8B5358EF54C3ull, 23695933555ull,  9095062555ull}, {0x67AB11594B0259E2ull, 28996300980ull, 42317029883ull}, {0xF1BE5C86FC344316ull, 15007099960ull, 48305102359ull} } },
	{131072,2271952979ull, { {0xF8137FB0E00BCC55ull, 26051697479ull, 27441427672ull}, {0x62A2E1A9A055358Bull,  4627749088ull, 28710079083ull}, {0x0ull, 0ull, 0ull} } },
	{147456,2548912547ull, { {0x5E213FBAFF82685Aull, 22451221320ull, 39039734900ull}, {0xC389436CC479E2E1ull,  7838917905ull, 66548008864ull}, {0x0ull, 0ull, 0ull} } },
	{163840,2825137853ull, { {0xAD3B129CD55821CFull,  2579580241ull, 65442845493ull}, {0x925F66B7986320BEull,  5355291557ull, 43930449714ull}, {0x0ull, 0ull, 0ull} } },
	{180224,3100703087ull, { {0x0B9AADC5599AA013ull, 17824815619ull, 32339793121ull}, {0x15A52A9E4B638559ull, 18530057795ull, 63757289449ull}, {0x0ull, 0ull, 0ull} } },
	{196608,3375668707ull, { {0x487B04D336385B62ull, 22477489910ull, 66079895176ull}, {0xC01F6617CD154EE8ull,  2883760509ull,  3966976333ull}, {0x0ull, 0ull, 0ull} } },
	{212992,3650084989ull, { {0x4F061E9732F6D86Eull, 25292467156ull, 29551251146ull}, {0x4BA7E5B4BF1051BEull, 15573066947ull, 25395261391ull}, {0x0ull, 0ull, 0ull} } },
	{229376,3923994593ull, { {0x1BCAD347E137ABECull, 19940191823ull, 48110255963ull}, {0x32677A46B243C568ull,   602580754ull, 43290561406ull}, {0x0ull, 0ull, 0ull} } },
	{245760,4197433843ull, { {0xC198853698CCAC7Full,  4937317521ull, 39532729153ull}, {0x74BA0019FC712DA6ull,  4461023520ull,  2965110143ull}, {0x0ull, 0ull, 0ull} } },
/* Larger require -shift 0: */
	/* Egregious: */
	{262144,4515590323ull, { {0x3B720039AA646317ull, 27060729660ull,  9603836312ull}, {0x9F8DD1E56E1063D9ull, 18665901522ull, 17899962874ull}, {0x0ull, 0ull, 0ull} } },
	{294912,5065885219ull, { {0x48B9C434ADA5C932ull, 16019066807ull, 17346586962ull}, {0x7F65539D123420A8ull,  2842988794ull, 23960120686ull}, {0x0ull, 0ull, 0ull} } },
	{327680,5614702259ull, { {0xE900D20947E96181ull, 12677175309ull, 46263325668ull}, {0xC0D077CD8E03FE07ull, 33414013825ull, 32326844283ull}, {0x0ull, 0ull, 0ull} } },
	{360448,6162190477ull, { {0x9E271959008B80D3ull, 31437330830ull, 36751682930ull}, {0xAE9E9BE256655FE1ull, 30096752173ull, 10472606829ull}, {0x0ull, 0ull, 0ull} } },
	{393216,6708471481ull, { {0x9DC838BAC4CCF77Cull, 26649031745ull,   200515496ull}, {0x59077693710A8DF3ull, 32893924229ull, 51130131345ull}, {0x0ull, 0ull, 0ull} } },
	{425984,7253646773ull, { {0xF04B3315B0250C1Aull, 25225037379ull, 54925972151ull}, {0x2227BF75276EDDE1ull,  6501658866ull, 45908319782ull}, {0x0ull, 0ull, 0ull} } },
	{458752,7797801821ull, { {0xDFB6E25E4DAE64D5ull,  2368273704ull, 20371186304ull}, {0x1D0AC271BAFB7B3Full, 19243408286ull, 36957506781ull}, {0x0ull, 0ull, 0ull} } },
	{491520,8341009997ull, { {0x15ED525006050B39ull, 32048847842ull, 54898451004ull}, {0x15EAF415588F56D8ull, 14874774582ull, 22782105046ull}, {0x0ull, 0ull, 0ull} } },
	{524288,8883334793ull, { {0x61776F2413CD7E79ull, 32946827752ull, 24059669414ull}, {0xE496C31B2534F520ull, 11663456039ull, 21852180932ull}, {0x0ull, 0ull, 0ull} } },
	/* Brobdingnagian: */
	/* Godzillian: */
	{     0,         0ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } }
};

// Reference Mersenne base-3 PRP residues - as with above MersVec[], use a 0-padding slot to store any user-set-exponent data:
struct testMers MvecPRP[numTest+1] =
{
/*                                         100-iteration residues:	                               1000-iteration residues:                */
/*	  FFTlen(K)     p              Res64           mod 2^35-1      mod 2^36-1               Res64           mod 2^35-1      mod 2^36-1     */
/*	    -----    --------     ----------------     -----------     -----------         ----------------     -----------     -----------    */
	/* Teensy:                                     [%34359738367  ][%68719476735  ]                         [%34359738367  ][%68719476735  ] */
	{     1,     22679ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{     2,     44657ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{     3,     66431ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{     4,     88019ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{     5,    109481ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{     6,    130873ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{     7,    152197ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{     8,    173431ull, { {0x5626DEB6B07622F3ull, 32131772281ull, 45561861712ull}, {0x2CFA9D066FA0AC59ull, 19159609871ull,  3629111748ull}, {0x91968AE0B907BC06ull, 18780821129ull, 60175052951ull} } },
	{     9,    194609ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    10,    215767ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    11,    236813ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    12,    257903ull, { {0x582A66DD4D8799EEull, 29535927244ull, 68660141246ull}, {0xEA5F22C54501846Full, 28855133772ull, 33619980259ull}, {0x5869D6DAA84D8F95ull, 13664256258ull,  9507781795ull} } },
	{    13,    278917ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    14,    299903ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    15,    320851ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	/* Tiny: */
	{    16,    341749ull, { {0xC736F1C2D213F1C1ull, 12211349626ull, 66509411860ull}, {0x8640A90521C2F7CCull, 26292223176ull, 67588668714ull}, {0x9EBE2EF30FB7D464ull, 12905240924ull, 64076380848ull} } },
	{    18,    383521ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    20,    425149ull, { {0xD3C192FF131CFC3Bull, 24519426320ull, 28595715402ull}, {0xCF0C17092AA78E04ull,  2271546708ull, 64056281496ull}, {0xE1F87989962DF48Eull, 26999099038ull, 52441645398ull} } },
	{    22,    466733ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    24,    508223ull, { {0x104F46B5B27DA109ull, 30634543418ull, 10636054863ull}, {0xFF40B3BCD8A23424ull, 30297179510ull, 62834095876ull}, {0x760A533DA9DE9552ull, 21594588451ull, 56782339492ull} } },
	{    26,    549623ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    28,    590963ull, { {0xECC45D40806A111Cull,  5787783130ull, 59666588997ull}, {0x4B8E8B5D4FA21D2Dull, 18930182921ull, 43534886672ull}, {0x876053199A486A8Full,  6807528947ull, 12093315482ull} } },
	{    30,    632251ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{    32,    673469ull, { {0xBD9646EE9AB18A7Aull, 26182028770ull, 59419684231ull}, {0x739B932AC594AA0Aull, 27600316683ull, 23680060816ull}, {0x0090EA7A76098D32ull, 10354482855ull, 28199252596ull} } },
	{    36,    755737ull, { {0xD3D7FCAAEC0D1CFAull, 16915465692ull, 12970156038ull}, {0x2ACB2B184E935C33ull,  2577059552ull, 25431200361ull}, {0xAD4F3FEAD246F78Dull, 11105939650ull, 27677803777ull} } },
	{    40,    837817ull, { {0x4BB1FE516DA56352ull, 14130534986ull, 57171707410ull}, {0xEDE384F60E7930DEull,  6782691233ull, 37080068160ull}, {0x985E59A168A2C374ull, 19638026482ull, 42674666120ull} } },
	{    44,    919729ull, { {0xBB36618A59F029BEull, 34052740788ull, 15139731592ull}, {0xAD419B6CE6F8A9B3ull, 25643929440ull,  3277180760ull}, {0xF1CCE9D2493ED0E7ull, 22901033162ull, 64906607173ull} } },
	{    48,   1001467ull, { {0xA750B4E8AEB33412ull, 14158218491ull, 67000020730ull}, {0x037EAC5E83798D97ull, 18831036867ull, 45448571503ull}, {0x8500430AFA74A003ull, 32712078146ull, 26108052951ull} } },
	{    52,   1083077ull, { {0x114B2898842D031Dull, 29992192737ull, 56153605976ull}, {0xC94D7EB86E333093ull, 25051877532ull, 36547212533ull}, {0x7B3B4B2B1E4AB4F1ull, 15083485543ull, 45190699824ull} } },
	{    56,   1164533ull, { {0xF3845BFD27BE17C0ull,  6122435706ull, 16830581570ull}, {0xEEA2A73E1B5B9F5Eull, 14009513529ull, 10409904152ull}, {0x8A1B59E0A1013082ull, 17983199887ull, 45703273488ull} } },
	{    60,   1245877ull, { {0xD3563072BBF5635Full, 26075147849ull, 59833238241ull}, {0x8228F9186BCCB825ull, 23336209236ull, 33753236087ull}, {0x8F50A4E950A208ECull, 10931878951ull, 43004788142ull} } },
	{    64,   1327099ull, { {0xA27D809BEA4B39EAull,  9126841940ull, 61797063350ull}, {0x5A458E090D0A11FBull, 12180510644ull, 38241257645ull}, {0xBAC955A84544EA74ull, 19644122492ull, 31588337102ull} } },
	{    72,   1489223ull, { {0x8115A8E974FE98BEull, 11713639178ull, 10466124455ull}, {0x853E9866ADE224FEull, 33773857851ull, 32975791386ull}, {0x30582D7F92A6D1FCull, 33748010249ull, 27388102096ull} } },
	{    80,   1650959ull, { {0x1D2DD2DF1933C0FCull,  6345481587ull, 37786169874ull}, {0x2D535BBF29F7CCD2ull,  4940281875ull, 42673535389ull}, {0xAF0A3BA384876D7Aull, 11399245247ull, 30073539441ull} } },
	{    88,   1812347ull, { {0x7FC9EA6769384320ull,  7891088336ull, 25044784554ull}, {0x0D1A4A0CD22AEB9Aull,  9596030110ull, 59365967842ull}, {0xA68AE2A59E471BE5ull, 15180142327ull, 52575114331ull} } },
	{    96,   1973431ull, { {0xB4CCA30966A4CE5Dull, 10270866356ull,   231429614ull}, {0xFCFC64781863674Eull, 23326776387ull, 52335987311ull}, {0x24EA8DEDBC00B415ull,  7671704365ull, 48042022252ull} } },
	{   104,   2134201ull, { {0x95DC857850180714ull, 11171528275ull, 64711382522ull}, {0x9FDC4093379BDB26ull, 18930952224ull,  1611454415ull}, {0xA141EA8967903F5Bull, 29876787665ull, 33823887927ull} } },
	{   112,   2294731ull, { {0x061CBC59DFB9CA06ull, 30425834940ull, 58336162616ull}, {0x580DF81D69FF6AD5ull,  8559569627ull, 10292910153ull}, {0x11BBBAD23995BB7Eull, 24990081881ull, 54500653284ull} } },
	{   120,   2455003ull, { {0x6580BD648772EDCEull, 13614650659ull, 27233253889ull}, {0x123237E17FF20692ull,  2520412521ull, 58497023977ull}, {0xF534CA2A0954221Bull, 16035573032ull, 55020965205ull} } },
	{   128,   2614999ull, { {0x61B3FF2749C421B2ull,  2127565925ull, 16135835978ull}, {0x5BF36C6CEA901199ull, 10377428760ull, 51455129205ull}, {0x44AE5B04991410FAull, 21353616283ull, 33587198015ull} } },
	{   144,   2934479ull, { {0xA322CF3BF0E43C86ull, 11064021322ull, 40809116502ull}, {0x556A5ADD65F0FB19ull, 11232143236ull, 16154650197ull}, {0x371C7434AB5079E1ull,  8183503318ull, 20346531393ull} } },
	{   160,   3253153ull, { {0xBE895574A3B5D29Dull,  4066068489ull, 27925792309ull}, {0x01972AEA521B5CDBull,  3673290729ull, 52399520094ull}, {0x37626B30822C90C2ull, 13476423895ull, 31936880647ull} } },
	{   176,   3571153ull, { {0xEB2F80B8B07B0C2Eull, 25525358118ull, 23794747344ull}, {0x669E15DC5325F9F3ull, 19996304981ull, 57830333278ull}, {0x20084A9830401673ull,  2415150132ull, 51845829949ull} } },
	{   192,   3888509ull, { {0x79FA4A9B6D077D31ull, 11817088595ull, 25793400384ull}, {0x986ED026065B3A3Dull, 28358614794ull, 61747377645ull}, {0xCD9CF17AA323A403ull, 25995389972ull, 16101059938ull} } },
	{   208,   4205303ull, { {0x20026FC75F1DACF2ull,   552886193ull, 35476382470ull}, {0x242488772DD4950Bull, 30174013897ull,  1712275816ull}, {0x370F8FD2CCB44245ull,  1046366655ull, 34353654023ull} } },
	{   224,   4521557ull, { {0xA5A3CA4D85A4DC25ull,  9322238657ull, 22846601560ull}, {0xD3046AC342DF5023ull,  1943060627ull, 46395479024ull}, {0xDCD16AAA39454FE8ull,  3893524834ull,  2191423213ull} } },
	{   240,   4837331ull, { {0x5E3C81457025A604ull, 12304143407ull, 29892018271ull}, {0x5E3C01108CBDE306ull, 28432606765ull, 35060019602ull}, {0xB82D80D8A8E662E0ull,  5053875088ull,  3674489982ull} } },
	/* Small: */
	{   256,   5152643ull, { {0x0A0ECD28374A3886ull, 10006324455ull, 31836113189ull}, {0xFAB1B0725803FC5Bull, 23503205738ull, 22955043669ull}, {0x25EB5A52ED298638ull, 21533519724ull, 56771934948ull} } },
	{   288,   5782013ull, { {0x8337D5DC89B9D7D8ull,  2750208646ull, 43710371927ull}, {0x81E025B2C36AF97Aull, 24587374922ull, 13476375540ull}, {0x7136299FEAE393A1ull,   873504035ull, 29264996801ull} } },
	{   320,   6409849ull, { {0xC1F07CFE2F8A2C90ull, 10288499049ull, 20419861489ull}, {0x962E7949F9A7E628ull, 15432666362ull, 30085592572ull}, {0x61E84019BC81B682ull, 17697124115ull, 46942817024ull} } },
	{   352,   7036339ull, { {0x336D9F3394B6F164ull, 17579907539ull, 13739842033ull}, {0x73AB57CC1CE38E7Cull, 30295988192ull, 23080433211ull}, {0x484554266231F919ull, 12854439638ull, 51687722348ull} } },
	{   384,   7661567ull, { {0x8CD5DBE23CC7D207ull, 27522372639ull,  4855853542ull}, {0x46453E1E483825B3ull, 26570117855ull,  8838577209ull}, {0x90838723D9A646F0ull, 19365339919ull, 54968311494ull} } },
	{   416,   8285659ull, { {0x2937D7BB53DB3A2Eull, 11180491559ull, 25143810934ull}, {0x453BB746F323DEAEull, 33021626783ull, 27443129753ull}, {0x4BAF1BA1324437A9ull,  7860885520ull, 59647517322ull} } },
	{   448,   8908723ull, { {0x547ADE1396F9C514ull, 18521157282ull, 59551009514ull}, {0x50092676FDE714A2ull, 16701877626ull, 43153619853ull}, {0x9E1092E438D8BCC8ull, 22087345475ull,  5454828580ull} } },
	{   480,   9530803ull, { {0x78D7C3BE6CEDAC01ull,  9241094902ull, 30808621041ull}, {0xC779524E4064A4B0ull, 18207111510ull, 33501852597ull}, {0x3EA6480E92D43A4Eull, 30247093504ull, 48234064574ull} } },
	{   512,  10151971ull, { {0x7145494EBE359D55ull, 23248741638ull, 48651777267ull}, {0x8BDB481C8392BFD5ull, 26285780620ull,  7866437445ull}, {0x08F52C235D33F535ull, 10426106456ull, 60796101403ull} } },
	{   576,  11391823ull, { {0xF9B279685D09CFCBull, 33118971289ull,  8380201274ull}, {0x631163D4E003D666ull, 15205345326ull, 56430348449ull}, {0xF001086F1B17E323ull, 14146155388ull, 43443184054ull} } },
	{   640,  12628613ull, { {0x19F38A1A24F47141ull, 28053480100ull, 63839732215ull}, {0x46CC84D00A6C3798ull,  7876783548ull, 52638210432ull}, {0x52ACFB11F83FB959ull, 28102345675ull, 48183960941ull} } },
	{   704,  13862759ull, { {0x7841B14AA240BC68ull,  9071998639ull,  8214384913ull}, {0x0518A8C7C3FB7B6Full,  9386753257ull,  4271506407ull}, {0xAB7F59730EAF000Cull, 25492803550ull, 39530850106ull} } },
	{   768,  15094403ull, { {0xD2BAF6916DE41B4Aull, 33907731243ull, 60082525243ull}, {0xE4780B94895718A8ull,  9102718496ull, 57395342555ull}, {0xA8AFFF234B1E3C32ull, 16763634198ull, 28601819201ull} } },
	{   832,  16323773ull, { {0xA87D786264D0EC8Aull, 12495933155ull, 37472327802ull}, {0xAE66B1040E98A5E8ull, 16343234311ull, 38946126334ull}, {0xCB0B755797CF128Bull, 10365202833ull, 20742591729ull} } },
	{   896,  17551099ull, { {0x59DAAE490AE3000Bull, 28683449829ull, 28452420814ull}, {0x766B2B663B520EF7ull,  2975453465ull, 42346868910ull}, {0x1CF5E66D07CD953Dull, 25387533318ull, 26555619006ull} } },
	{   960,  18776473ull, { {0xCDE871BF25E32E3Bull, 31673267310ull, 52228110774ull}, {0x5366B5BE07F98107ull, 23248985466ull, 26763145241ull}, {0x58DB393E9B69FE09ull,  1673827116ull, 65412337112ull} } },
	{  1024,  20000047ull, { {0xF804962A99E06D0Cull,  1860502464ull, 53057883331ull}, {0x89356E1C8F6F1E48ull, 22320943580ull, 35169637716ull}, {0x94C9EE9C64F737E4ull, 31567228626ull, 56393904459ull} } },
	{  1152,  22442237ull, { {0x0D9595BF361781F9ull, 20413575841ull, 37106981458ull}, {0xE5BFCB9E43A3E9F3ull, 20960264993ull, 64719249630ull}, {0x44515095B41E5606ull, 24782945059ull, 39032256339ull} } },
	{  1280,  24878401ull, { {0x2CBF037107A1847Dull, 22690755724ull, 56236025943ull}, {0x59F656025983667Bull,  5142037023ull, 11990415177ull}, {0x3A343F6424076489ull, 18282185034ull, 25811451519ull} } },
	{  1408,  27309229ull, { {0xBF16EFC49A336692ull, 31238875220ull, 58423375260ull}, {0xB3B79B9B95B76297ull, 26402154389ull, 62028968784ull}, {0x2E0C3A1000F7B0F9ull, 32941800112ull,   747986439ull} } },
	{  1536,  29735137ull, { {0x1C67BD358F25FE5Aull, 18399031066ull, 37515566273ull}, {0xF246C8AB58A1C74Eull, 21054258805ull, 47020117222ull}, {0x207FE2BADA1C43C5ull,  3284146171ull,  2281976028ull} } },
	{  1664,  32156581ull, { {0x98A7C4D249C14A4Full, 20144152830ull, 56678356879ull}, {0xE5AEA4A688D5928Cull, 27069838966ull, 37265104256ull}, {0xE63FD6C1CF930EB2ull, 23334182391ull, 19780994645ull} } },
	{  1792,  34573867ull, { {0xD42C4E48467FFFB7ull, 31980303049ull, 62848444803ull}, {0xE5C6F0D767DCBCB0ull, 26394307125ull, 19768489106ull}, {0xC8F429962FC0C412ull, 31086854861ull, 48294508690ull} } },
	{  1920,  36987271ull, { {0xFF48564F7201B487ull, 25019473293ull, 10999816574ull}, {0x432CACFAC7536F31ull, 23972386709ull, 39242248517ull}, {0xA5FDDBF4742D26EDull,  7858251857ull, 16067544118ull} } },
	/* Medium: */
	{  2048,  39397201ull, { {0x926CA6804251B3EEull, 28060637476ull, 49790249351ull}, {0x6489E97CBC5580AAull,  9625584894ull, 62826847344ull}, {0x053BC856B8BE4270ull,  4691188379ull, 60624374837ull} } },
	{  2304,  44207087ull, { {0x02BB72063670F51Bull, 12694994014ull, 31300999610ull}, {0xDB77D5DBFD8297EDull, 20254637099ull, 52711100645ull}, {0x391A0B2307926600ull, 33289407481ull, 33078655462ull} } },
	{  2560,  49005071ull, { {0xA58C6E638188B85Full, 30546989250ull, 52706842241ull}, {0xE2126248B4846CD6ull,   685992662ull, 35816889167ull}, {0xCA778F9AB26281B8ull, 11413212083ull,  1212979690ull} } },
	{  2816,  53792327ull, { {0xE5B6FD041C5C5B08ull, 16802279919ull, 45744609376ull}, {0x905C64E6B30CBBABull, 10202620340ull, 13224782271ull}, {0x056DDDAFC35968E5ull,  6286892553ull, 37514012293ull} } },
	{  3072,  58569809ull, { {0x7F39D5ED6103EBABull, 16756172503ull, 66487365574ull}, {0x8830F2A0B6737836ull,  2942497587ull, 46149388407ull}, {0x0CFE2DB14EEFBF07ull, 26575645314ull, 46479029624ull} } },
	{  3328,  63338459ull, { {0x984077364B7DC067ull, 31200355828ull, 53223232350ull}, {0xC460BCC10E7C4633ull, 28854750816ull, 17126512478ull}, {0x4E012FD2746144CAull, 28091592352ull,  9463381637ull} } },
	{  3584,  68098843ull, { {0x2629F428E8AAF210ull, 31208335732ull, 28724687211ull}, {0x27C14BD10BBCC373ull, 25428252468ull, 24335967243ull}, {0xFA54F0C508E93DFDull, 10745895462ull,  6297999612ull} } },
	{  3840,  72851621ull, { {0xC9790A9AE318D8A0ull,  7965585773ull, 38904230464ull}, {0xD40DA752D00A58B8ull,  3810739518ull,  7705774890ull}, {0x0ED03C062FF2A606ull,  2639280055ull, 37383026335ull} } },
	{  4096,  77597293ull, { {0x8A37C0FF09636A7Full,   830455861ull, 40725642853ull}, {0x525B3401B12A2878ull,  6322052778ull, 50560800498ull}, {0x7E130AB5C1CF9FC3ull,   821663107ull, 49875216353ull} } },
	{  4608,  87068977ull, { {0xAFD146065C8A2F8Aull,  3265972347ull, 27188046914ull}, {0xC3E83DF5C838746Dull,  9185368986ull, 65232442937ull}, {0x4273C31CB66A685Bull, 22667296906ull, 46963850727ull} } },
	{  5120,  96517019ull, { {0x3BBBD10EC4BD191Eull, 13411649980ull, 56736468345ull}, {0xB06AF4101E3CC028ull, 14290899010ull, 58971130200ull}, {0xDD87E6A5A509FFC2ull,  2239081373ull, 53398491455ull} } },
	{  5632, 105943723ull, { {0xA0204A10ADA1D777ull,   325219564ull, 48693854932ull}, {0x61D45F7DAFD716FCull,  2255330486ull, 33206674222ull}, {0x0818EBD02C01A08Eull, 31797592248ull, 31752563546ull} } },
	{  6144, 115351063ull, { {0x35C7DE3E3F3F8D94ull, 31618972943ull, 55520653728ull}, {0x956FB2D53BEC843Cull, 11148946878ull, 33245408201ull}, {0x4CAA4EBEB4EB4951ull, 11777136074ull, 25148258060ull} } },
	{  6656, 124740697ull, { {0x4681EC451768F31Aull, 12959366213ull, 18617894816ull}, {0x05DFD257F3AE1F02ull, 12757038996ull,  2260817322ull}, {0x1CC290E749E46688ull,  6693583950ull, 43197953245ull} } },
	{  7168, 134113933ull, { {0xE6A04B7C79535282ull, 31728172652ull, 58237681894ull}, {0x3D1CF94DCA031EA6ull,  6612981966ull, 36523212505ull}, {0x5AAB8CCA656621B2ull,  7673677458ull, 48483227415ull} } },
	{  7680, 143472073ull, { {0x6EAAFFBBB2B27868ull, 12622370397ull,  9274290901ull}, {0x513410FA3142A1ECull, 27159455539ull, 38629970432ull}, {0x1D3D8F16FA48CC4Cull, 27219308286ull, 15713263754ull} } },
	{  8192, 152816047ull, { {0x55D755491E9A5BE7ull, 17408991436ull, 55204850562ull}, {0xA160152F199821D0ull, 23171738795ull, 32522593027ull}, {0x10447B5958C7153Dull, 25649892134ull, 63692031319ull} } },
	{  9216, 171465013ull, { {0x7FC9A3D6580A67DBull,  2493822844ull, 32653389776ull}, {0xEDBDBED649AC1C07ull, 31826896222ull, 51396833159ull}, {0x1259EFF4D4B88CE2ull, 20123348812ull, 61034401374ull} } },
	{ 10240, 190066777ull, { {0xC875BAE2D9D23F8Eull, 10787379418ull, 62215501884ull}, {0xD0534B8C3FD4FEBDull, 22561335508ull, 19377764663ull}, {0xD7B93BF968F0F50Full, 24193017740ull,  1698596575ull} } },
	{ 11264, 208626181ull, { {0x121330EF1C9C65D4ull, 15421852243ull, 16454197259ull}, {0x9342E60165C9515Bull, 21714755947ull,  9528514349ull}, {0xB73DA3A3DCFC2715ull, 18785773681ull,  5941932830ull} } },
	/* Large: */
	{ 12288, 227147083ull, { {0x939900344B3A9CF5ull, 10791587738ull, 51376518661ull}, {0x9B6A405153110744ull, 33911048215ull, 10244698095ull}, {0x372C128DE18F44A8ull, 22661022105ull, 28319883481ull} } },
	{ 13312, 245632679ull, { {0x77BAD267187DB572ull,  6978161239ull, 31566453664ull}, {0x843B2C80EFD985D4ull, 31638655555ull, 34886706969ull}, {0xFEA3E15FF92C0B9Eull, 25087739150ull, 65596702716ull} } },
	{ 14336, 264085733ull, { {0x3B6C5DD137A06F3Cull,  2041442582ull, 41068697072ull}, {0x1D4E6F465FB90F6Eull, 30561782032ull,  1263429588ull}, {0xAA93E30434811C8Dull,  2228362920ull,  8335956471ull} } },
	{ 15360, 282508657ull, { {0x006B5DC0A65002D1ull, 31937249556ull, 54892782386ull}, {0xCC239E0E7FCBC7E2ull, 22387326720ull, 24840066078ull}, {0x804B979EB89922BAull, 26535016499ull, 63049971720ull} } },
	{ 16384, 300903377ull, { {0xE1FFB7FA51A666BDull, 16509353676ull, 33290659747ull}, {0x729F90E2F3D1C751ull, 11110715400ull, 57675528900ull}, {0xE584214AFA269421ull, 33543551870ull, 17530147104ull} } },
	{ 18432, 337615277ull, { {0x4EA0D844C4D4A158ull, 15144962431ull, 63334776550ull}, {0xB23CA173BE980CD5ull, 34244619245ull, 21356443084ull}, {0xBCE24CE8A48EF4C8ull, 18408960714ull, 40047346011ull} } },
	{ 20480, 374233309ull, { {0xAC8A22C76BFE8A7Dull, 30655495979ull, 20979915815ull}, {0xB31E512AD3D23426ull, 19514358995ull, 11856185197ull}, {0x38B9508197A2F880ull, 17914520610ull, 40675602543ull} } },
	{ 22528, 410766953ull, { {0xA5A013AD452E9CA1ull,  6158071328ull, 30004265076ull}, {0x5863D58C9141BAADull,  5535389685ull,  4594364305ull}, {0x6D0DA537655ADB8Full, 28534262546ull, 56406949849ull} } },
	{ 24576, 447223969ull, { {0xEC329D5F8EA6B575ull,  8127808506ull, 48282079888ull}, {0xBAD0E2BA110A8248ull,   551404710ull, 48950059099ull}, {0x294D8DA9F0957D69ull, 18477018631ull,  4404820879ull} } },
	{ 26624, 483610763ull, { {0x51652EAE126ACC75ull, 20340182567ull, 55729433233ull}, {0x184D46BA91B4286Dull, 31213394420ull, 54386141290ull}, {0x7AD715347087A95Dull, 29627983156ull, 55271838273ull} } },
	{ 28672, 519932827ull, { {0xE2C5CA7A92E9D9AFull, 10027203195ull, 55451827923ull}, {0x6154C17BA65CDE41ull,  1952185223ull,  7287620568ull}, {0x87268C55E4BCB257ull, 23433413895ull, 47801213689ull} } },
	{ 30720, 556194803ull, { {0xEA0DD84F89EC9475ull,  2640397683ull, 10711381600ull}, {0x9366059D474B5372ull, 26599065025ull, 33329252442ull}, {0x3D17BC75E5E3DDEFull, 32874051574ull, 57034690579ull} } },
	{ 32768, 592400713ull, { {0xA4869F5FE9F96DBDull, 19624228853ull, 26960586543ull}, {0x092EC7F7DA96A330ull, 20551131024ull, 16345984908ull}, {0x38C8047BB23D3146ull, 18285343891ull, 62090151768ull} } },
	{ 36864, 664658101ull, { {0x6DB632AA8386525Bull,  9946879119ull, 35779904199ull}, {0xCFA85C70B49A5B3Full, 14532066374ull, 10923711144ull}, {0x680471EBEF0500BCull, 29423921028ull, 21769143331ull} } },
	{ 40960, 736728527ull, { {0x999B28F81B1B37C0ull, 15281289408ull, 12055047400ull}, {0x851D858E64D2D7DEull, 15017284283ull,  6105199500ull}, {0x2022E24AC582ABDFull, 26611160419ull, 12727750935ull} } },
	{ 45056, 808631029ull, { {0xC6E167CAB66F6C66ull, 24733824729ull, 55248890737ull}, {0x7FE8B4A8C80520CDull, 26074748676ull, 17574876945ull}, {0x1BD35D3BB502A327ull,  7929585559ull, 13807068803ull} } },
	{ 49152, 880380937ull, { {0xBE4EDE9750058245ull,  1628233733ull, 63365088770ull}, {0xF855B49E2A5653E2ull,  7420489469ull, 44867201251ull}, {0xB98AA627050B671Aull, 10196185952ull,  1461018327ull} } },
	{ 53248, 951990961ull, { {0xCDCEB6E9A3CF5747ull, 13166317856ull, 46957323254ull}, {0x1EC13F17D33FC1F0ull, 30008301556ull, 19194697819ull}, {0x91827F71C70B72E2ull,  6915899554ull, 51031186257ull} } },
	{ 57344,1023472049ull, { {0x428763AC63EEBB06ull,  4046036533ull, 65887115399ull}, {0x9819D4EECD27A260ull,  2725007656ull, 27008383182ull}, {0x0B9FF561BB7E6A8Eull, 27229975800ull, 68594865215ull} } },
	{ 61440,1094833457ull, { {0xD8419FF6C0F1BBE0ull,  1663876259ull, 34786957257ull}, {0x2A79788D7D8D853Bull, 20881114859ull, 33590133948ull}, {0x2D3034F9BEFD96D5ull,  2175887418ull, 16970800813ull} } },
	/* Huge: */
	{ 65536,1154422469ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 73728,1295192531ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 81920,1435594063ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 90112,1575664327ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 98304,1715433593ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{106496,1854927187ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{114688,1994166553ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{122880,2133169847ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{131072,2271952979ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{147456,2548912547ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{163840,2825137853ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{180224,3100703087ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{196608,3375668707ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{212992,3650084989ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{229376,3923994593ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{245760,4197433843ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
/* Larger require -shift 0: */
	/* Egregious: */
	{262144,4515590323ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{294912,5065885219ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{327680,5614702259ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{360448,6162190477ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{393216,6708471481ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{425984,7253646773ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{458752,7797801821ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{491520,8341009997ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{524288,8883334793ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	/* Brobdingnagian: */
	/* Godzillian: */
	{     0,         0ull, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } }
};

/* PARI/GP script for generating more FFT-length/maxP entries in the above table uses condensed code from given_N_get_maxP():

Bmant = 53;
AsympConst = 0.6;
ln2inv = 1.0/log(2.0);
N = [some starting power of 2, e.g. 128M --> m = 2^27]
currMult = 8; // <-- cycle through the supported FFT lengths starting with m
N *= currMult/(currMult-1); N/1024	// <-- Get Kdoubles form of N
ln_N=log(1.0*N);lnln_N=log(ln_N);l2_N=ln2inv*ln_N;lnl2_N=log(l2_N);l2l2_N=ln2inv*lnl2_N;lnlnln_N=log(lnln_N);l2lnln_N=ln2inv*lnlnln_N;Wbits=0.5*(Bmant-AsympConst-0.5*(l2_N+l2l2_N)-1.5*(l2lnln_N));maxExp2=Wbits*N;i=round(maxExp2-0.5)
[Now increment i if even, then iterate until TRUE:]
i-=2;isprime(i)
*/

struct testFerm{
	int fftLength;		/* FFT length in K (i.e. 4 means an array of 4K doubles) */
	uint32 Fidx;			/* Fermat number index */
	struct res_triplet	res_t[3];	/* 100,1000 and 10000-iteration SH-residue triplets */
};

/* Array of 100-iteration reference values for Fermat self-tests. Only allow Fidx >= 14: */
#define FermArrayIdxOffset		14	/* Amount to subtract from m to get the prestored residues for F_m */
#define numFerm		20
struct testFerm FermVec[numFerm+1] =
{
/*                                   100-iteration residues:                                  1000-iteration residues:              */
/* FFTlen(K) Fidx           Res64           mod 2^35-1      mod 2^36-1               Res64           mod 2^35-1      mod 2^36-1     */
/*   ------- ----      ----------------     -----------     -----------         ----------------     -----------     -----------    */
	/*                                    [%34359738367  ][%68719476735  ]                         [%34359738367  ][%68719476735  ] */
//	{     2,  13u, { {0xC8FC67EA3A1AC788ull, 29592689237ull, 35156594447ull}, {0xBE489C2CF00D582Aull,  3108135315ull, 47125592449ull}, {0x0ull, 0ull, 0ull} } },
	{     1,  14u, { {0xDB9AC520C403CB21ull,   342168579ull, 59244817440ull}, {0xF111F12732CCCB0Full, 24848612524ull, 66609820796ull}, {0x78738D068D641C2Cull, 12664273769ull, 29297626750ull} } },
	{     2,  15u, { {0x3B21A6E55ED13454ull, 28379302213ull, 15546218647ull}, {0x4784657F2A36BE74ull,   617376037ull, 44891093359ull}, {0x589BFE53458FFC14ull, 12200390606ull, 46971422957ull} } },
	{     4,  16u, { {0xAAE76C15C2B37465ull, 20013824731ull,  2261076122ull}, {0x42CC2CBE97C728E6ull, 30814966349ull, 44505312792ull}, {0xEED00D8AE6886440ull, 19057922301ull, 53800020279ull} } },
	{     8,  17u, { {0xFFA16CDC8C87483Cull, 20917337408ull, 26110327818ull}, {0x43CAB295FFB2661Full, 18197605796ull,  9842643677ull}, {0x5C7B0049549D4174ull,  8923959253ull, 40303785249ull} } },
	{    16,  18u, { {0x7C6B681485EB86DBull,  5745147782ull, 50521157289ull}, {0x8193BD41931E9DE8ull, 19662968587ull, 51102742548ull}, {0x9D6A242467D28700ull, 12912307491ull, 23293425575ull} } },
	{    32,  19u, { {0x529E54642A813995ull, 17797950508ull, 32039741221ull}, {0xE24EAE4B153EE86Bull, 11155350666ull, 49866866361ull}, {0x787FABD98DD5FEC0ull, 10784283812ull, 50254721650ull} } },
	{    64,  20u, { {0x64629CED6E218018ull,  8485981669ull, 53977437340ull}, {0xA380121F6FD26B2Aull, 15876203498ull, 36314727556ull}, {0x352CB92DABA82A8Bull,  6625203514ull, 20044302250ull} } },
	{   128,  21u, { {0x1DE0171591038250ull, 33758422990ull,  8269940507ull}, {0x1B331FBB41AF33D7ull, 17971032338ull,  2929392342ull}, {0x2E7F9D30EAFC7D47ull, 28599205675ull, 44913594527ull} } },
	{   256,  22u, { {0x9201143390F3828Dull,  1749100092ull, 46602233256ull}, {0x1B331FBB41AF33D7ull, 17971032338ull,  2929392342ull}, {0x4380A384A80C079Eull, 30799637462ull, 32805936529ull} } },
	{   512,  23u, { {0x9C3F8E29B397B32Bull,  1094055486ull, 13316822657ull}, {0xBD642EA0479D8FF0ull, 31625967305ull, 57187857233ull}, {0xC0F18B226C20AE50ull, 31325111828ull, 10029614864ull} } },
	{  1024,  24u, { {0xDB9F01963ED9DC8Bull, 27887793041ull, 13169874547ull}, {0x40F2DECE9C351236ull,  9074286032ull, 38590877049ull}, {0x84FEE61771F20003ull, 19739002855ull, 42109658313ull} } },
	{  2048,  25u, { {0x376C33921E5F675Full, 13022327996ull, 46818697393ull}, {0xA51F8577A407CB75ull,  9865976783ull, 35171498411ull}, {0x2B1D829E05C85CA3ull, 27914095914ull, 58025822027ull} } },
	{  4096,  26u, { {0xA42BECD80DAEC4CBull, 10087739060ull, 25252768685ull}, {0xECC9408A7295401Dull,  5904751941ull, 58967745948ull}, {0x0C054D0376BD8E9Eull, 24761357761ull, 23215069147ull} } },
	{  8192,  27u, { {0xFB69E377519D8CE6ull, 15449775614ull, 51221672039ull}, {0x24898E3BEB59DCE6ull, 24957168001ull,  2072452827ull}, {0x443DE289690070EDull, 30674247970ull, 29955017069ull} } },
	{ 16384,  28u, { {0xA4FF6F8C3CB38B85ull, 18933356966ull, 30899345457ull}, {0x8B451AF25E8CC50Eull,   674652743ull, 39963850167ull}, {0x0F38D26B35C6794Cull,  7003106751ull, 60469235270ull} } },
	{ 32768,  29u, { {0xAFBF110B593E26F6ull, 32666279868ull, 18995112582ull}, {0xA6B643FF24C6ADC1ull, 15753158767ull, 13965270144ull}, {0xB194096866F68C59ull, 19667273394ull,  3552165634ull} } },
	{ 65536,  30u, { {0x68B1BDA5D6BAE04Bull,  3347054148ull, 47892955488ull}, {0x361F30024AF9FE26ull, 26693502373ull, 67933083515ull}, {0x3A70D98DA3ED809Aull, 10877397803ull, 41746600776ull} } },
	{131072,  31u, { {0x00C36B4AB38FC326ull, 12487821860ull, 64847210796ull}, {0x243D46185B4FAAC8ull, 14946978930ull, 42870390813ull}, {0xDD3651B25892F424ull,  6237393494ull, 22880395670ull} } },
	{262144,  32u, { {0xFB5A3146BE2CA886ull, 33629944997ull, 19559359478ull}, {0x0479E68514EAD529ull,    59701496ull, 58397447084ull}, {0x78CD466FDDA4F156ull,  5752278548ull, 18141301956ull} } },
	{524288,  33u, { {0x25D9822BD1555EB9ull, 10052071893ull, 20645266428ull}, {0xACF521F811129C9Eull, 25612341094ull, 61523429244ull}, {0x269405A9B18B3310ull,  4910011908ull,  2895484666ull} } },
	{     0,   0u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } }
};

/***************************************************************************************
Main program: User-argument parsing, self-testing and production-run dispatch.
By way of reference, here is the complete list of currently supported assignment types
and related params as per the PrimeNet v5 API - Mlucas supports only a subset of these:

	struct work_unit {        // One line from the worktodo file
		int    work_type;    // Type of work to do
		char    assignment_uid[33]; // Primenet assignment ID
		char    extension[9];    // Optional save file extension
		double    k;        // K in k*b^n+c
		unsigned long b;    // B in k*b^n+c
		unsigned long n;    // N in k*b^n+c
		signed long c;        // C in k*b^n+c
		unsigned long minimum_fftlen;// Minimum FFT length to use.  Zero means default fftlen.
		double    sieve_depth;    // How far it has been trial factored
		double    factor_to;    // How far we should trial factor to
		int    pminus1ed;    // TRUE if has been P-1 factored
		double    B1;        // ECM, P-1, P+1 - Stage 1 bound
		double    B2;        // ECM, P-1, P+1 - Stage 2 bound
		double    B2_start;    // P-1 - Stage 2 start
		int    nth_run;    // P+1 - 1 for start 2/7, 2 for start 6/5, 3+ for random start
		unsigned int curves_to_do; // ECM - curves to try
		double    curve;        // ECM - Specific curve to test (debug tool)
		double    tests_saved;    // Pfactor - primality tests saved if a factor is found
		unsigned int prp_base;    // PRP base to use
		int    prp_residue_type; // PRP residue to output -- see primenet.h
		int    prp_dblchk;    // True if this is a doublecheck of a previous PRP
		int    cert_squarings; // Number of squarings required for PRP proof certification
		char    *known_factors;    // ECM, P-1, P+1, PRP - list of known factors
		char    *comment;    // Comment line in worktodo.txt
	}
***************************************************************************************/
int 	main(int argc, char *argv[])
{
	int		retVal=0;
	uint64	Res64, Res35m1, Res36m1;
	char	stFlag[STR_MAX_LEN];
	uint64	i64arg, expo = 0, lo,hi;
	uint32	iarg = 0, iters = 0, k = 0, maxFFT, findex = 0;
	double	darg;
	int		new_cfg = FALSE;
	int		i,j, idum, nargs, scrnFlag, maxAllocSet = FALSE, nbufSet = FALSE;
	int		start = -1, finish = -1, modType = 0, testType = 0, selfTest = 0, userSetExponent = 0, xNum = 0;
#ifdef MULTITHREAD
	// Vars for mgmt of mutually exclusive arg sets; 'core' is specifically for hwloc-including builds:
	int		nthread = 0, cpu = 0, core = 0;
#endif
	char *cptr = 0x0;
	int		quick_self_test = 0, fftlen = 0, radset = -1;
	uint32 numrad = 0, rad_prod = 0, rvec[10], rvec2[10];	/* Temporary storage for FFT radices */
	double	runtime,wruntime, runtime_best,wruntime_best, tdiff;	// v20: w-prefixed are weighted by associated ROEs
	double	roerr_avg = 0, roerr_max = 0;
	int		radix_set, radix_best, nradix_set_succeed;

	uint32 mvec_res_t_idx = 0;	/* Lookup index into the res_triplet table */
	uint32 new_data;
	struct res_triplet new_res = {0ull,0ull,0ull};
	struct testMers*MvecPtr = MersVec;	// Set this to point at either MersVec (the default) or MvecPRP, depending on test type

/* Enable this and set upper loop bound appropriately to regenerate a quick list of primes
just below the upper limit for each FFT lengh in some subrange of the self-tests:
*/
#if 0
	// Note this is 32-bit only!
	for(i = 68; i < numTest; i++) {
		j = given_N_get_maxP(MersVec[i].fftLength << 10);
		if((j&1) == 0) --j; 	// make sure it's odd
		for(;;) {
			if(is_prime(j)) {
				fprintf(stderr, "%10u\n",j);
				break;
			}
			j -= 2;
		}
	}
	exit(0);
#endif

	// In 32-bit Linux, may need to up the stacklimit from the defaults
	// to avoid SIGSEGV faults during the alloc-heavy self-test suites.
	set_stacklimit_restart(argv);
	fprintf(stderr, "\n    Mlucas %s\n", VERSION);
	fprintf(stderr, "\n    %s\n\n", HOMEPAGE);

#ifdef macintosh
	argc = ccommand(&argv);			/* Macintosh CW */
#endif
	// v18: Enable access to argc/argv outside main():
	global_argv = argv;

	ASSERT((MersVec[numTest-1].fftLength != 0) &&  (MersVec[numTest].fftLength == 0), "numTest != MersVec allocated size!");
	ASSERT((MvecPRP[numTest-1].fftLength != 0) &&  (MvecPRP[numTest].fftLength == 0), "numTest != MvecPRP allocated size!");
	ASSERT((FermVec[numFerm-1].fftLength != 0) &&  (FermVec[numFerm].fftLength == 0), "numFerm != FermVec allocated size!");

	/*...check that various data types are of the assumed length
	and do some other basic sanity checks:
	*/
	Mlucas_init();

	/* Uncomment this and the print statement inside the given_N_get_maxP
	function to see the full list of maxP's for all supported FFT lengths: */
  #if 0
	get_default_fft_length(PMAX);
  #endif

	scrnFlag = 0;	/* Do Not echo output to stddev */
	// Currently only primality-test and PRP modes supported - init to former by default, switch to PRP if -prp flag set:
	testType = TEST_TYPE_PRIMALITY;

	/******** command-line-argument processing while() loop: ********/
	nargs = 1;
	while(argv[nargs])
	{
		strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
		if(nargs > argc) {	// == no longer applies since e.g. -prp requires no numeric arg and can come last:
			fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
			print_help();
		}

		if(stFlag[0] != '-') {
			fprintf(stderr, "*** ERROR: Illegal command-line option '%s'\n", stFlag);
			print_help();
		}

		if(STREQ(stFlag, "-h")) { print_help(); }

		/* Mersenne self-test: requires a user-set exponent, FFT length or one of the supported -s arguments below: */
		if(STREQ(stFlag, "-s"))
		{
			selfTest = TRUE;
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			for(;;) {
				if(STREQ(stFlag, "a") || STREQ(stFlag, "all")) {	/* all, which really means all the non-Huge-and-larger sets */
					start = 0; finish = numTeensy + numTiny + numSmall + numMedium + numLarge;
					break;
				}

				finish = 0;

				start = finish; finish += numTeensy;
				if(STREQ(stFlag, "tt") || STREQ(stFlag, "teensy"))	/* teensy   */
					break;
				start = finish; finish += numTiny;
				if(STREQ(stFlag, "t") || STREQ(stFlag, "tiny"))		/* tiny   */
					break;
				start = finish; finish += numSmall;
				if(STREQ(stFlag, "s") || STREQ(stFlag, "small"))	/* small  */
					break;
				start = finish; finish += numMedium;
				if(STREQ(stFlag, "m") || STREQ(stFlag, "medium"))	/* medium */
					break;
				start = finish; finish += numLarge;
				if(STREQ(stFlag, "l") || STREQ(stFlag, "large"))	/* large  */
					break;
				start = finish; finish += numHuge;
				if(STREQ(stFlag, "h") || STREQ(stFlag, "huge") || STREQ(stFlag, "xl"))		/* huge   */
					break;
				start = finish; finish += numEgregious;
				if(STREQ(stFlag, "e") || STREQ(stFlag, "egregious") || STREQ(stFlag, "xxl"))
					break;
				start = finish; finish += numBrobdingnagian;
				if(STREQ(stFlag, "b") || STREQ(stFlag, "brobdingnagian") || STREQ(stFlag, "xxxl"))
					break;
				start = finish; finish += numGodzillian;
				if(STREQ(stFlag, "g") || STREQ(stFlag, "godzillian") || STREQ(stFlag, "gojira"))
					break;

				fprintf(stderr, "*** ERROR: Illegal argument '%s' to -s flag\n", stFlag);
				print_help();
			}
			modType  = MODULUS_TYPE_MERSENNE;
		}

		else if(STREQ(stFlag, "-maxalloc"))	// maxalloc arg is max %-of-available-mem to use
		{
			ASSERT(nbufSet == FALSE, "Only one of -maxalloc and -pm1_s2_buf flags may be used!");
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			darg = strtod(stFlag,&cptr);
			// Must be > 0:
			ASSERT((darg > 0), "maxalloc (%%-of-available-mem to use) argument must be > 0 ... halting.");
			// Max-%-of-RAM-to-use currently stored in MAX_RAM_USE ... later will multiply by (available system RAM in MB):
			MAX_RAM_USE = darg;
			maxAllocSet = TRUE;
			printf("INFO: User specified -maxalloc: will use up to %u%% of available system RAM.\n",MAX_RAM_USE);
		}

		else if(STREQ(stFlag, "-pm1_s2_nbuf"))	// pm1_s2_nbuf arg is max %-of-available-mem to use
		{
			ASSERT(maxAllocSet == FALSE, "Only one of -maxalloc and -pm1_s2_buf flags may be used!");
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			darg = strtod(stFlag,&cptr);
			// Must be > 0:
			ASSERT((darg > 0), "pm1_s2_nbuf argument must be integer ... halting.");
			// Max-%-of-RAM-to-use currently stored in MAX_RAM_USE ... later will convert to floating-fraction and multiply by (available system RAM in MB):
			PM1_S2_NBUF = darg;
			nbufSet = TRUE;
			printf("INFO: User specified -pm1_s2_nbuf = %u.\n",PM1_S2_NBUF);
		}

		else if(STREQ(stFlag, "-iters"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			i64arg = atol(stFlag);
			// Must be < 2^32:
			ASSERT(!(i64arg>>32), "#iters argument must be < 2^32 ... halting.");
			iters = (uint32)i64arg;
		}

		else if(STREQ(stFlag, "-fft") || STREQ(stFlag, "-fftlen"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			// v20: default is still integer-FFT-length in Kdoubles, but add support for [float]M,
			// where floating-point arg must be exactly representable, such that [float]*2^10 is integer:
			i64arg = -1ull;
			darg = strtod(stFlag,&cptr);
			if(strlen(cptr)) {
				if(STREQ(cptr,"M"))
					i64arg = darg*1024;
				else if(STREQ(cptr,"K"))
					i64arg = darg;
				else {
					ASSERT(0, "The only non-numeric suffixes allowed for the argument to -fft are K and M");
				}
			} else
				i64arg = darg;
			// Must be in range [MIN_FFT_LENGTH_IN_K,MAX_FFT_LENGTH_IN_K], def'd in Mdata.h:
			if(i64arg < MIN_FFT_LENGTH_IN_K || i64arg > MAX_FFT_LENGTH_IN_K) {
				sprintf(cbuf  , "ERROR: FFT-length argument = %llu, must be in range [%u,%u]K\n",i64arg,MIN_FFT_LENGTH_IN_K,MAX_FFT_LENGTH_IN_K);
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}
			fftlen = (uint32)i64arg;	// Note this is the REAL-vector FFT length
			if((i = get_fft_radices(fftlen, 0, 0x0, 0x0, 0)) != 0) {
				sprintf(cbuf  , "ERROR: FFT length %d K not available.\n",fftlen);
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}
			// If user has supplied a set of complex-FFT radices, their product must equal half the real-FFT length:
			if(rad_prod) { ASSERT((rad_prod>>9) == fftlen,"Product of user-supplied set of complex-FFT radices must equal half the real-FFT length!"); }
		}

		/* v19.1: Enhance the -radset flag to take either an index into the big table in get_fft_radices(),
		or an actual set of comma-separated FFT radices? If the expected -radset[whitespace]numeric arg-pair
		is immediately followed by a comma, assume it's a set of radices, read those in, check whether said
		set is supported and if so, set radset to the corresponding table-index numeric value: */
		else if(STREQ(stFlag, "-radset"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);

			// Check if it's a comma-separated actual set of complex-FFT radices:
			char_addr = stFlag;
			cptr = strchr(char_addr,',');
			if(!cptr) {	// It's a radix-set index
				i64arg = atol(stFlag);
				// Must be < 2^32:
				ASSERT(i64arg < 20, "radset-index argument must be < 2^32 ... halting.");
				radset = (uint32)i64arg;
			} else {	// It's a set of complex-FFT radices
				numrad = 0;
				while(0x0 != (cptr = strchr(char_addr,','))) {
					// Copy substring into cbuf and null-terminate:
					strncpy(cbuf,char_addr,(cptr-char_addr));	cbuf[cptr-char_addr] = '\0';
					// Convert current radix to long and sanity-check:
					i64arg = atol(cbuf);	ASSERT(!(i64arg>>12), "user-supplied radices must be < 2^12 ... halting.");
					rvec[numrad++] = (uint32)i64arg;
					char_addr = cptr+1;
				}
				// A properly formatted radix-set arg will end with ',[numeric]', with the numeric in char_addr:
				i64arg = atol(char_addr);	ASSERT(!(i64arg>>12), "user-supplied radices must be < 2^12 ... halting.");
				rvec[numrad++] = (uint32)i64arg;
				rvec[numrad] = 0;	// Null-terminate the vector just for aesthetics
				// Compute the radix product and make sure it's < 2^30, constraint due to the (fftlen < 2^31) one:
				rad_prod = 1; i64arg = 1ull;
				for(i = 0; i < numrad; i++) {
					i64arg *= rvec[i];	ASSERT(!(i64arg>>30), "Product of complex-FFT radices supplied via -radset argument must be < 2^32 ... halting.");
				}
				rad_prod = (uint32)i64arg;
				// If user has supplied a real-FFT length (in Kdoubles) via -fftlen, product of the complex-FFT radices must equal half that value:
				if(fftlen) {
					ASSERT((rad_prod>>9) == fftlen,"Product of user-supplied set of complex-FFT radices must equal half the real-FFT length!");
				} else {
					fftlen = rad_prod>>9;	// If user supplies fftlen via cmd-line arg after -radset, that's OK,
								// we'll overwrite fftlen with user-supplied value and repeat the above check then
				}
				// Now loop over the available radix sets for this FFT length and see if we find a match,
				// in which case j holds #radices and rvec2[] holds the corresponding radices on return.
				// Note that get_fft_radices takes real-FFT length in terms of Kdoubles:
				i = 0;
				while(!get_fft_radices(fftlen, i++, (uint32*)&j, rvec2, 10)) {	// 0-return means radset index is in range
					if(j != numrad) continue;
					// #radices matches, see if actual complex radices do
					for(j = 0; j < numrad; j++) {
						if(rvec[j] != rvec2[j]) break;
					}
					if(j == numrad) {	// i -= 1 here to undo post-increment in above get_fft_radices() call
						radset = i-1; break;
					}
				}
				// The init-value of radset -1 getting overwritten with something >= 0 means success:
				ASSERT(radset >= 0, "User-supplied set of complex-FFT radices not supported.");
			}
		}

		else if(STREQ(stFlag, "-shift"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			i64arg = atol(stFlag);
			// Must be < 2^32, though store in a uint64 for later bignum-upgrades:
			ASSERT(!(i64arg>>32), "shift argument must be < 2^32 ... halting.");
			RES_SHIFT = i64arg;
		}

		// v20: Add p-1 support:
		else if(STREQ(stFlag, "-b1"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			i64arg = atol(stFlag);
			// Must be < 2^32:
			ASSERT(!(i64arg>>32), "P-1 Stage 1 bound must be < 2^32 ... halting.");
			B1 = (uint32)i64arg;
			ASSERT(testType != TEST_TYPE_PRP, "b1-argument implies P-1 factoring; that and PRP-test types not simultaneously specifiable.");
			testType = TEST_TYPE_PM1;
		}
		else if(STREQ(stFlag, "-b2")) {
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			// Allow Stage 2 bounds to be > 2^32:
			B2 = atol(stFlag);
			ASSERT(testType != TEST_TYPE_PRP, "b2-argument implies P-1 factoring; that and PRP-test types not simultaneously specifiable.");
			testType = TEST_TYPE_PM1;
		}
		else if(STREQ(stFlag, "-b2_start")) {
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			// Allow Stage 2 bounds to be > 2^32:
			B2_start = atol(stFlag);
			ASSERT(testType != TEST_TYPE_PRP, "b2_start-argument implies P-1 factoring; that and PRP-test types not simultaneously specifiable.");
			testType = TEST_TYPE_PM1;
		}

		else if(STREQ(stFlag, "-nthread"))
		{
		#ifndef MULTITHREAD
			ASSERT(0,"Multithreading must be enabled in build to permit -nthread argument!");
		#else
			ASSERT(cpu == FALSE && core == FALSE,"Only one of -nthread, -cpu and -core flags permitted!");
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			i64arg = atol(stFlag);
			// Must be < 2^32:
			ASSERT(!(i64arg>>32), "nthread argument must be < 2^32 ... halting.");
			NTHREADS = (uint32)i64arg;
			nthread = TRUE;
			// Use the same affinity-setting code here as for the -cpu option, but simply for cores [0:NTHREADS-1]:
			sprintf(cbuf,"0:%d",NTHREADS-1);
			parseAffinityString(cbuf);
		#endif
		}

		else if(STREQ(stFlag, "-cpu"))
		{
		#ifndef MULTITHREAD
			ASSERT(0,"Multithreading must be enabled in build to permit -cpu argument!");
		#else
			ASSERT(nthread == FALSE && core == FALSE,"Only one of -nthread, -cpu and -core flags permitted!");
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			parseAffinityString(stFlag);
			cpu = TRUE;
		#endif
		}

	#if INCLUDE_HWLOC
		else if(STREQ(stFlag, "-core"))
		{
		#ifndef MULTITHREAD
			ASSERT(0,"Multithreading must be enabled in build to permit -core argument!");
		#else
			ASSERT(cpu == FALSE && nthread == FALSE,"Only one of -nthread, -cpu and -core flags permitted!");
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			NTHREADS = parseAffinityTriplet(stFlag,TRUE);	// 2nd-arg = TRUE: Use hwloc-generated topology, via '-core lo:hi[:threads_per_core]'
			if(NTHREADS > MAX_THREADS) {
				fprintf(stderr,"ERROR: NTHREADS [ = %d] must not exceed those of available logical cores = 0-%d!\n",NTHREADS,MAX_THREADS-1);
				exit(EXIT_FAILURE);
			}
			core = TRUE;
		#endif
		}
	#endif // endif(INCLUDE_HWLOC?)

		/******************************************************************************************/
		/* MAKE SURE ALL OTHER FLAG-PROCESSING SECTIONS SET userSetExponent TO A NONZERO VALUE!!! */
		/******************************************************************************************/

		else if(STREQ(stFlag, "-m") || STREQ(stFlag, "-mersenne"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			expo = atol(stFlag);
			userSetExponent = 1;
			// Use 0-pad slot in MvecPtr[] to store user-set-exponent data - that can point to either MersVec
			// or MvecPRP. if -prp is invoked after the -m {+int} flag, at that point copy any exponent set
			// here from MersVec[numTest].exponent to MvecPRP[numTest].exponent :
			MvecPtr[numTest].exponent = expo;
			start = numTest; finish = start+1;
			modType = MODULUS_TYPE_MERSENNE;
		}

		else if(STREQ(stFlag, "-prp"))	// This flag optionally takes a numeric base arg, and trips us into PRP-test mode
		{
			if(nargs < argc) {
				strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
				if(isdigit(stFlag[0])) {
					PRP_BASE = atol(stFlag);
					if(PRP_BASE+1 == 0) {
						snprintf(cbuf,STR_MAX_LEN*2, "*** ERROR: Numeric arg to -prp flag, '%s', overflows uint32 field.\n", stFlag);
						ASSERT(0,cbuf);
					}
				}
				else
					--nargs;
			}
			// Use 0-pad slot in MvecPRP[] to store user-set-exponent data:
			ASSERT(MvecPtr == MersVec,"-prp flag invoked, but MvecPtr does not reflect the default MersVec init-value!");
			MvecPtr = MvecPRP;
			if(MersVec[numTest].exponent) {
				MvecPtr[numTest].exponent = MersVec[numTest].exponent;
				MersVec[numTest].exponent = 0ull;
			}
			if(!PRP_BASE)	// If user has not (yet) specified base, default = 3; if -base comes later in cmd-line, it will override
				PRP_BASE = 3;
			modType = MODULUS_TYPE_MERSENNE;
			testType = TEST_TYPE_PRP;
		}

		else if(STREQ(stFlag, "-base"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			i64arg = atol(stFlag);
			PRP_BASE = (uint32)i64arg;
		}

		else if(STREQ(stFlag, "-f") || STREQ(stFlag, "-fermat"))
		{
			strncpy(stFlag, argv[nargs++], STR_MAX_LEN-1);
			i64arg = atol(stFlag);
			// Must be < 2^32:
			ASSERT(!(i64arg>>32), "Fermat-number-index argument must be < 2^32 ... halting.");
			findex = (uint32)i64arg;
			/* Make sure the Fermat number index is in range: */
			if(findex < 13 || findex > 63) {
				fprintf(stderr, " Fermat number index must be in the range [13,63].\n");
				return ERR_EXPONENT_ILLEGAL;
			}
			userSetExponent = 1;
			FermVec[numFerm].Fidx = findex;
			start = numFerm; finish = start+1;
			modType = MODULUS_TYPE_FERMAT;
		}

		else
		{
			fprintf(stderr, "*** ERROR: Unrecognized flag %s.\n", stFlag);
			print_help();
		}
	}	/* end of command-line-argument processing while() loop */

	// Nov 2020: Sanity-check any p-1 bounds:
	if(testType == TEST_TYPE_PM1) {
		ASSERT((modType == MODULUS_TYPE_MERSENNE || modType == MODULUS_TYPE_FERMAT) && userSetExponent, "P-1 in command-line mode requires a Mersenne or Fermat-number modulus to be specified via '-m [int]' or '-f [int]'.");
		pm1_check_bounds();
	}

	if(!modType)
		modType = MODULUS_TYPE_MERSENNE;

	// Now that have determined the modType, copy any user-set FFT length into the appropriate field:
	if(fftlen) {
		/* Don't set userSetExponent here, since -fftlen can be invoked without an explicit exponent */
		if(modType == MODULUS_TYPE_FERMAT) {
			FermVec[numFerm].fftLength = fftlen;
			start = numFerm; finish = start+1;
		} else {
			MvecPtr[numTest].fftLength = fftlen;
			start = numTest; finish = start+1;
		}
	}
	// If user has specified a radix set, make sure an FFT length has also specified:
	if(radset != -1) {
		if(modType == MODULUS_TYPE_FERMAT)
			iarg = FermVec[numFerm].fftLength;
		else
			iarg = MvecPtr[numTest].fftLength;

		if(iarg == 0) {
			sprintf(cbuf  , "*** ERROR: Must specify a valid FFT length on command line before -radset argument!\n");
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
		}

		/* Make sure it's a valid radix set index for this FFT length: */
		if((i = get_fft_radices(iarg, radset, (uint32*)&idum, 0x0, 0)) != 0) {
			if     (i == ERR_FFTLENGTH_ILLEGAL)
				sprintf(cbuf  , "ERROR: FFT length %d K illegal!\n", iarg);
			else if(i == ERR_RADIXSET_UNAVAILABLE)
				sprintf(cbuf  , "ERROR: radix set index %d for FFT length %d K exceeds maximum allowable of %d.\n",radset, iarg, idum-1);
			else
				sprintf(cbuf  , "ERROR: Unknown error-code value %d from get_fft_radices(), called with radix set index %d, FFT length %d K\n",i,radset, iarg);

			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
		}

	}

	// Use selfTest == TRUE or -iters (in the single-FFT-length-timing case) to differentiate between timing-test and production runs.
	// In fact it eases the logic to explicitly set selfTest = TRUE whenever iters is set (but not nec. the converse), so do that here:
	if(iters) selfTest = TRUE;

	if(modType == MODULUS_TYPE_MERSENNE && !selfTest)
	{
		if(userSetExponent) {
			ASSERT(start > 0, "userSetExponent = TRUE but self-test starting-index unset!");
			sprintf(cbuf, "ERROR: Production-run-mode [-iters not invoked] does not allow command-line\nsetting of exponent - that must be read from the %s file.\n",WORKFILE);
			ASSERT(0,cbuf);
		} else if(start == -1) {
			start = numTest; finish = start+1;
		}
		if(radset != -1) {
			sprintf(cbuf, "ERROR: Production-run-mode [-iters not invoked] allows command-line setting of\nFFT length, but not the radix set - that must be read from the mlucas.cfg file.\n");
			ASSERT(0,cbuf);
		}
	ERNST_MAIN:
		if((retVal = ernstMain(modType,testType,0,MvecPtr[start].fftLength,0,0,0,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime)) != 0)
		{
			printMlucasErrCode(retVal);

			/* If need to run a timing self-test at a particular FFT length, do that and then try again... */
			if((retVal & 0xff) == ERR_RUN_SELFTEST_FORLENGTH) {
				quick_self_test = TRUE;
				selfTest = TRUE;
				k = (uint32)(retVal >> 8);
				if((i = get_fft_radices(k, 0, 0x0, 0x0, 0)) != 0) {
					sprintf(cbuf, "ERROR: FFT length %d K not available.\n",k);
					fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
				}

			/**** IF POSSIBLE, USE ONE OF THE STANDARD TEST EXPONENTS HERE, SO CAN CHECK RES64s!!! ****/
				for(i = 0; i < numTest; i++) {
					if(MvecPtr[i].fftLength == k) {
						userSetExponent = 0;
						start = i; finish = start+1;
						break;
					}
				}
				if(i == numTest) {
					userSetExponent = 1;
					MvecPtr[numTest].exponent = convert_base10_char_uint64(ESTRING);
					MvecPtr[numTest].fftLength = k;
					start = numTest; finish = start+1;
				}

				modType = MODULUS_TYPE_MERSENNE;
				goto TIMING_TEST_LOOP;
			}
			/* ...Otherwise barf. */
			else {
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}
		} else {
			fprintf(stderr, "\n  Done ...\n\n");
			exit(EXIT_SUCCESS);
		}
	}
	/* If user specified iters and FFT length but no exponent, get default Mersenne exponent for that FFT length: */
	else if(modType == MODULUS_TYPE_MERSENNE)
	{
		if(MvecPtr[start].exponent == 0)
		{
			i = MvecPtr[start].fftLength;
			ASSERT(i > 0                  ,"Require i > 0                  ");
			ASSERT(i <=MAX_FFT_LENGTH_IN_K,"Require i <=MAX_FFT_LENGTH_IN_K");

			// If FFT length is not represented in reference-residue array, find nearest prime <= 0.99*given_N_get_maxP(FFT length):
			for(j = 0; j < numTest; j++) {
				if(i == MvecPtr[j].fftLength) break;
			}
			if(i != MvecPtr[j].fftLength) {
				hi = (99*given_N_get_maxP(i<<10)/100) | 0x1;	// Make sure starting value is odd. v21: Cut to 99% of pmax_rec
				lo = hi - 1000;	if(lo < PMIN) lo = PMIN;
				for(expo = hi; expo >=lo; expo -= 2) {
					if(isPRP64(expo)) {
						MvecPtr[numTest].exponent = expo;
						break;
					}
				}
				if(expo < lo || lo >= hi) {
					fprintf(stderr, "ERROR: unable to find a prime in the interval %llu <= x <= %llu.\n", lo, hi);
					ASSERT(0,"0");
				}
			} else {	/* Use the corresponding entry of MvecPtr: */
				start = j; finish = start+1;
			}
		}
		/* If user specified exponent but no FFT length, get default FFT length for that exponent: */
		else if(MvecPtr[numTest].exponent && (MvecPtr[numTest].fftLength == 0))
			MvecPtr[numTest].fftLength = get_default_fft_length((uint64)(MvecPtr[numTest].exponent));
	}
	else if(modType == MODULUS_TYPE_FERMAT)
	{
		if(FermVec[start].Fidx == 0) {
			i = FermVec[start].fftLength;
			ASSERT(i > 0                  ,"Require i > 0                  ");
			ASSERT(i <=MAX_FFT_LENGTH_IN_K,"Require i <=MAX_FFT_LENGTH_IN_K");

			if(i > FermVec[numFerm-1].fftLength)	/* Computing a new-largest entry? */
				FermVec[numFerm].Fidx = (i << 4);
			else {	/* Find the corresponding entry of FermVec: */
				for(lo = 0; lo < numFerm; lo++) {
					if(FermVec[lo].fftLength >= i) {
						start = lo; finish = start+1;	/* Using >= here allows for non-power-of-2 FFT lengths */
						break;
					}
				}
				if(lo >= numFerm) {
					fprintf(stderr, "ERROR: unable to find FFT length %d K in the Reference Residue table.\n", i);
					ASSERT(0,"0");
				}
			}
		}
		/* If user specified exponent but no FFT length, get default FFT length for that exponent: */
		else if(findex && (FermVec[numFerm].fftLength == 0))
			FermVec[numFerm].fftLength = get_default_fft_length((uint64)1 << findex);
	} else{
		ASSERT(0,"modType not recognized!");
	}

TIMING_TEST_LOOP:

	if(selfTest) {
		fprintf(stderr, "\n           Mlucas selftest running.....\n\n");
		/* We have precomputed 100, 1000 and 10000-iteration residues for the predefined self-test exponents: */
		if( userSetExponent && (modType == MODULUS_TYPE_MERSENNE) ) {
			fprintf(stderr, "\n********** Non-default exponent - you will need to manually verify that the residue **********");
			fprintf(stderr, "\n********** triplets output for this self-test match for all FFT radix combinations! **********\n\n");
		} else if( iters && (iters != 100) && (iters != 1000) && (iters != 10000) ) {
			fprintf(stderr, "\n********** #Iters not one of [100,1000,10000] - you will need to manually verify that  **********");
			fprintf(stderr, "\n********** residue triplets output for self-test match for all FFT radix combinations! **********\n\n");
		} else if( testType == TEST_TYPE_PRP && PRP_BASE != 3 ) {
			fprintf(stderr, "\n********** Non-default PRP base - you will need to manually verify that the residue **********");
			fprintf(stderr, "\n********** triplets output for this self-test match for all FFT radix combinations! **********\n\n");
		}
		fprintf(stderr, "/****************************************************************************/\n\n");
	}
	/* FFT-radix configuration file is named mlucas.cfg or fermat.cfg,
	depending on whether Mersenne or Fermat number test is being done.

	The .cfg file stores data about the optimal FFT parameters at each runlength
	on the host platform. This file is generated by running one of the automated self-tests.

	The first line in the file is assumed to contain the VERSION datum of the version program
	which was used to generate the timing data contained therein. For pre-v3.0 (for which the
	format was different than this) this is of course not true, but that's OK, since in such
	a case (or one in which the .cfg file is unreadable or nonexistent) a new set of timing
	runs is triggered, by way of the cfgNeedsUpdating() function returning a nonzero value.

	Each subsequent line in the file is assumed to start with a pair of integers, followed by an optional
	comment or supplemental data (we may choose to make some such data mandatory or at least meaningful in
	terms of their influencing execution if present, at some future date). The first integer is the FFT length
	(in units of 1K), i.e. the number of doubles used to store the residues for the number being tested.
	(Note that the program does a complex-data FFT, i.e. treats this vector as a complex vector of half the length.)
	The second integer on each line is the index of the optimal set of FFT radices with respect to the supported
	FFT radix sets for the particular FFT length. (The available FFT radices can be found in the large table in
	the switch(n){} structure of the function get_fft_radices() in the get_fft_radices.c source file.)
	The optimal set of complex radices and the per-iteration timing it yields on the host platform are written
	to the comment field right of each FFT-length/optimal-radix-set pair.

	NOTE that the ordering of the {FFT length | Index of best radix set} lines is not considered to be
	important, e.g. by the get_preferred_fft_radix() routine, thus if we find that we need timing data for
	a particular FFT length, we can append the results to any pre-existing ones in the relevant .cfg file.
	That saves us the headache of having to contemplate doing on-the-fly re-sorting of the file contents.
	*/
	if(modType == MODULUS_TYPE_FERMAT)
		strcpy(CONFIGFILE,"fermat.cfg");
	else
		strcpy(CONFIGFILE,"mlucas.cfg");

	/* Determine mode in which to open the file. There are several cases to consider here:

	0) There is no valid (i.e. readable, containing valid data) .cfg file in the run directory.

	ACTION: In this case we open the .cfg file in write mode (i.e. blow away any previous file contents
	or create a new file if none existed before) and write fresh timing entries
	for the FFT lengths of the particular self-test being performed.

	1) There is a valid .cfg file, but with a version string corresponding to a code version
	which is sufficiently older that new timings are needed.

	ACTION: Same as (0) - generate new .cfg file data.

	2) There is a valid .cfg file with a compatible (in the sense of not requiring retiming)
	program-version entry on line 1.

	ACTION: 2 possibilities here:

		(a) If a timing test: Open the file in append mode and generate/write timing data for FFT lengths
		in the current set of self-tests. If some FFT lengths of the old and new self-test timing data overlap,
		we will wind up with mutliple entries for the repeated FFT lengths, but that's OK, since in such cases
		the get_preferred_fft_radix() function will simply grab the entry with the best timing from among those.

		(b) If a regular LL test, look for an entry for the requisite FFT length. If multiple entries for that
		length, get_preferred_fft_radix() function will simply grab the entry with the best timing from among those.
		If a valid entry is found for that length but an entry for a larger FFT length is also found having a better
		timing, the larger FFT length will be used for the run. If no valid entry is found for the default FFT length,
		a quick set of timing tests at the length will be run in order to generate one, then go back to square (b).
	*/
	FILE_ACCESS_MODE[0] = FILE_ACCESS_READ;
	FILE_ACCESS_MODE[1] = FILE_ACCESS_UPDATE;	/* Include this in mlucas_fopen, so can
													(a) check whether it exists and is writeable;
													(b) if so, read the first line;
												in one fell swoop. */
	fp = mlucas_fopen(CONFIGFILE,FILE_ACCESS_MODE);
	if(!fp) {
		sprintf(cbuf, "INFO: Unable to find/open %s file in %s mode ... creating from scratch.\n", CONFIGFILE, FILE_ACCESS_MODE);
		fprintf(stderr,"%s",cbuf);
		new_cfg = TRUE;
	} else {
		/* Program version string assumed to be stored in line 1 of .cfg file -
		if it's not, the .cfg file is by definition outdated:
		*/
		if(fgets(in_line, STR_MAX_LEN, fp))
			new_cfg = cfgNeedsUpdating(in_line);
		else
			new_cfg = TRUE;

		fclose(fp); fp = 0x0;
	}

	/* If existing .cfg file is for current program version, file-append mode is appropriate;
	otherwise set up to re-open in write mode, which will clear any previous file contents:
	*/
	FILE_ACCESS_MODE[1] = FILE_FORMAT_ASCII;
	if(new_cfg)
		FILE_ACCESS_MODE[0] = FILE_ACCESS_WRITE;
	else
		FILE_ACCESS_MODE[0] = FILE_ACCESS_APPEND;

	/* What's the max. FFT length (in K) for the set of self-tests? */
	maxFFT = MvecPtr[finish-1].fftLength;

	for (xNum = start; xNum < finish; xNum++)    /* Step through the exponents */
	{
		new_data = FALSE;	Res64 = Res36m1 = Res35m1 = 0ull;

		/* If it's a self-test [i.e. timing test] and user hasn't specified #iters, set to default: */
		if(selfTest && !iters) {
			if(NTHREADS > 4)
				iters = 1000;
			else
				iters = 100;
		}

		if(iters == 100 || iters == 1000 || iters == 10000) {
			mvec_res_t_idx = NINT( log((double)iters)/log(10.) ) - 2;	/* log10(iters) - 2, use slower NINT rather than DNINT here since latter needs correct rounding mode */
			ASSERT(mvec_res_t_idx < 3,"main: mvec_res_t_idx out of range!");
			// Use empty-data-slot at top of MersVec[] or MvecPRP[], respectively, for primality & prp single-case tests:
			if( (modType == MODULUS_TYPE_MERSENNE && MvecPtr[xNum].res_t[mvec_res_t_idx].sh0 == 0)
			 || (modType == MODULUS_TYPE_FERMAT   && FermVec[xNum].res_t[mvec_res_t_idx].sh0 == 0) )
			{	// New self-test residue being computed:
				new_data = TRUE;
				new_res.sh0 = new_res.sh1 = new_res.sh2 = 0ull;
			}
		}

		/* Init best-radix-set, best-runtime and #radix-sets-which-succeeded for this FFT length: */
		runtime_best = wruntime_best = 0.0;
		radix_best = -1;
		nradix_set_succeed = 0;

		/* If user-specified radix set, do only that one: */
		if(radset >= 0)
			radix_set = radset;
		else
			radix_set = 0;

		if(modType == MODULUS_TYPE_MERSENNE)
			iarg = MvecPtr[xNum].fftLength;
		else if(modType == MODULUS_TYPE_FERMAT)
			iarg = FermVec[xNum].fftLength;

		// For SIMD builds only use FFT lengths which are at least a multiple of 4 Kdoubles:
	  #if 0//def USE_SSE2	*** v21: Removed array padding for FFT < 32K, disabled this "skip"
		if(iarg & 3) continue;
	  #endif

		while(get_fft_radices(iarg, radix_set, &NRADICES, RADIX_VEC, 10) == 0)	/* Try all the radix sets available for this FFT length. */
		{
			if(modType == MODULUS_TYPE_FERMAT)
			{
				Res64   = FermVec[xNum].res_t[mvec_res_t_idx].sh0;
				Res35m1 = FermVec[xNum].res_t[mvec_res_t_idx].sh1;
				Res36m1 = FermVec[xNum].res_t[mvec_res_t_idx].sh2;
				retVal = ernstMain(modType,testType,(uint64)FermVec[xNum].Fidx    ,iarg,radix_set,maxFFT,iters,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime);
			}
			else if(modType == MODULUS_TYPE_MERSENNE)
			{
				Res64   = MvecPtr[xNum].res_t[mvec_res_t_idx].sh0;
				Res35m1 = MvecPtr[xNum].res_t[mvec_res_t_idx].sh1;
				Res36m1 = MvecPtr[xNum].res_t[mvec_res_t_idx].sh2;
				retVal = ernstMain(modType,testType,(uint64)MvecPtr[xNum].exponent,iarg,radix_set,maxFFT,iters,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime);
			}
			else if(testType == TEST_TYPE_PM1) {
				retVal = ernstMain(modType,testType,(uint64)MvecPtr[xNum].exponent,iarg,radix_set,maxFFT,iters,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime);
			}
			else
				ASSERT(0,"Unsupported modulus and/or test type!");

			// (retVal != 0) relates to dangerously high ROEs, use maxErr to decide whether to accept radix set.
			/*** (to-do: factor in #occurrences) ***/
			// AME,MME contain avgMaxErr and maxMaxErr for iter-interval, ROE_ITER holds #maxErr > 0.40625;	ROE_VAL recapitulates MME
			// (meaning we can use it to store some additional ROE-related measure should that become desirable in the future):
			if(selfTest && ( !userSetExponent && ((iters == 100) || (iters == 1000) || (iters == 10000)) )
			&&	 ( (iters ==   100 && MME > 0.40625)	// Even a single ROE unacceptable for such a short run
				|| (iters ==  1000 && MME > 0.42   )
				|| (iters == 10000 && MME >=0.4375 ) ) )
			{
				fprintf(stderr, "***** Excessive level of roundoff error detected - this radix set will not be used. *****\n");
				if(radset >= 0)	// If user-specified radix set, do only that one:
					goto DONE;
				runtime = 0.0; ++radix_set; continue;
			}
			else if(retVal)	// Bzzzzzzzzzzzt!!! That answer is incorrect. The penalty is death:
			{
				printMlucasErrCode(retVal);
				if( !userSetExponent && ((iters == 100) || (iters == 1000) || (iters == 10000)) )
					fprintf(stderr, "Error detected - this radix set will not be used.\n");
				if(radset >= 0)	// If user-specified radix set, do only that one:
					goto DONE;
				runtime = 0.0; ++radix_set; continue;
			}
			else if(radset >= 0)	// If user-specified radix set, do only that one:
			{
				goto DONE;
			}
			else if(new_data)	// New self-test residues being computed - write to .cfg file if get a consensus value:
			{
				if(!new_res.sh0)	// First of the available radix sets:
				{
					new_res.sh0 = Res64  ;
					new_res.sh1 = Res35m1;
					new_res.sh2 = Res36m1;
					nradix_set_succeed++;	// This assumes when doing such a new-data-selftest the initial radset gives the correct result,
											// which may not be justified, but we can add fancier consensus-weighing logic later, if needed.
				}
				else if			// All subsequent radix sets must match to produce a consensus result:
				(
					new_res.sh0 == Res64
				 &&	new_res.sh1 == Res35m1
				 &&	new_res.sh2 == Res36m1
				)
				{
					nradix_set_succeed++;
				} else {
					runtime = 0.0; ++radix_set; continue;
				}
			} else {	// If not a new-data self-tests (i.e. it's a regular -s one), getting here means the current radset succeeded:
				nradix_set_succeed++;
			}

			/* 16 Dec 2007: Added the (runtime != 0) here to workaround the valid-timing-test-but-runtime = 0
			issue on systems with round-to-nearest-second granularity of the clock() function: */
		#if 1
			if( runtime_best == 0.0 || ((runtime != 0) && (runtime < runtime_best)) ) {
				runtime_best = runtime;
		#else	//*** This ROE-weighted scheme is having some unexpected suboptimal results - revisit later: ***
			// v20: Multiply runtime in sec by 1 + [weighted sum of ROE data] to prefer
			// more-accurate radix-combos among ones with similar runtimes:
			wruntime = runtime*(1 + 4*AME + MME);	// Give AME more weight as it is less noisy
			if( runtime_best == 0.0 || ((runtime != 0) && (wruntime < wruntime_best)) ) {
				// Use ROE-weighted runtimes in above comparison but save actual runtime for later cfg-file writing:
				runtime_best = runtime;
				wruntime_best = wruntime;
		#endif
				radix_best = radix_set;
				// Dec 2014: Changed from reporting [min,max]MME for all radsets @ the given FFT length - which is
				// confusing since what one really cares about is the ROE levels of the best-timed radset, not the ones
				// which were discarded as being slower - to reporting [AME,MME] for the best-timed radset:
				// May 2015: It would help if I actually moved this code inside the (new best runtime?) if()...
				roerr_avg = AME;
				roerr_max = MME;
			}

			fprintf(stderr, "\n");
			++radix_set;
		}

		// If get no successful reference-Res64-matching results, or less than half of results @this FFT length match, skip it:
		if(radix_best < 0 || runtime_best == 0.0 || nradix_set_succeed < (radix_set+1)/2)
		{
			sprintf(cbuf, "WARNING: %d of %d radix-sets at FFT length %u K passed - skipping it. PLEASE CHECK YOUR BUILD OPTIONS.\n",nradix_set_succeed,radix_set,iarg);
			fprintf(stderr,"%s", cbuf);
		}
		/* If get a nonzero best-runtime, write the corresponding radix set index to the .cfg file: */
		else
		{
			sprintf(cbuf, "INFO: %d of %d radix-sets at FFT length %u K passed - writing cfg-file entry.\n",nradix_set_succeed,radix_set,iarg);
			fprintf(stderr,"%s", cbuf);

			/* Divide by the number of iterations done in the self-test: */
		#ifdef MULTITHREAD	// In || mode the mod_square routines use getRealTime() to accumulate wall-clock time, thus CLOCKS_PER_SEC not needed
			tdiff = runtime_best/iters;
		#else
			tdiff = runtime_best/((double)iters*CLOCKS_PER_SEC);
		#endif
			fp = mlucas_fopen(CONFIGFILE,FILE_ACCESS_MODE);
			if(!fp) {
				sprintf(cbuf  , "INFO: Unable to open %s file in %s mode ... \n", CONFIGFILE, FILE_ACCESS_MODE);
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}

			/* Put code version on line 1.
			Only want to do this once; subsequent mlucas_fopen/fprintf are in append mode:
			*/
			if(new_cfg && FILE_ACCESS_MODE[0] == FILE_ACCESS_WRITE) {
				fprintf(fp, "%s\n", VERSION);
				FILE_ACCESS_MODE[0]=FILE_ACCESS_APPEND;
			}

			/* Print per-iteration time in sec and best FFT radix set into comment field.
			   Use fixed-width padded formats so #bytes per line is constant (needed to allow update mode).

			   5/01/2007: eliminate the print of radix_best, since that makes things dependent on the tables
			   in get_fft_radices() never changing - we only care about the actual *radices*, not the index into
			   the get_fft_radices case{} selector which yields them.
			*/
			fprintf(fp, "%10d  msec/iter = %7.2f  ROE[avg,max] = [%10.9f, %10.9f]  radices = ", iarg, 1000*tdiff, roerr_avg, roerr_max);
			/*
				Double-check existence of, and print info about,
				best radix set found for this FFT length:
			*/
			if(get_fft_radices(iarg, radix_best, &NRADICES, RADIX_VEC, 10) != 0) {
				sprintf(cbuf  , "ERROR: alleged best-radix-set index %u is unsupported.\n",radix_best);
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}
			/* Zero-pad the radices-printing to the full length of the RADIX_VEC array
			so each line has same length (needed to allow update mode):
			*/
			for(i = 0; i < 10; i++){ fprintf(fp,"%3u",RADIX_VEC[i]); };

			/* If it's a new self-test residue being computed, add the SH residues to the .cfg file line */
			if(new_data)
				fprintf(fp, "\tp = %s: %d-iter Res mod 2^64, 2^35-1, 2^36-1 = %016llX, %11.0f, %11.0f",ESTRING,iters,new_res.sh0,(double)new_res.sh1,(double)new_res.sh2);

			fprintf(fp,"\n");
			fclose(fp); fp = 0x0;

			/* if just adding entry for a single FFT length needed for current exponent, return to here: */
			if (quick_self_test) {
				quick_self_test = selfTest = 0; start = numTest; finish = start+1;
				goto ERNST_MAIN;
			}
		}

		fprintf(stderr, "/ **************************************************************************** /\n\n");
	}

DONE:
	fprintf(stderr, "\n  Done ...\n\n");
	return(0);
}

/******************/

/* v18: Add reparse-command-line-args-to-extract-desired-flag-value function for residue shift.
Assumes the basic program-start command-parsing has been done, thus dispenses with the former's bad-arg checking.
Can expand as needed to cover other options as the need arises.
*/
uint64	parse_cmd_args_get_shift_value(void)
{
	uint64	i64arg = -1ull;	// 0 is a valid value, so init = -1 (= 0xfff...fff in unsigned form)
	char	stFlag[STR_MAX_LEN];
	/******** command-line-argument processing while() loop: ********/
	int i, nargs = 1;
	while(global_argv[nargs])
	{
		strncpy(stFlag, global_argv[nargs++], STR_MAX_LEN-1);
		if(STREQ(stFlag, "-shift"))
		{
			strncpy(stFlag, global_argv[nargs++], STR_MAX_LEN-1);
			/* Convert the shift argument to a uint64: */
			i64arg = 0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++) {
				if(isdigit(stFlag[i])) {
					i64arg = 10*i64arg + (stFlag[i]-CHAROFFSET);
					/* Check for overflow: */
					if(i64arg % (uint64)10 != (uint64)(stFlag[i]-CHAROFFSET))
					{
						snprintf(cbuf,STR_MAX_LEN*2, "*** ERROR: -shift argument %s overflows uint64 field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
					}
				} else {
					snprintf(cbuf,STR_MAX_LEN*2, "*** ERROR: Non-numeric character encountered in -shift argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
				}
			}
		}
	}
	return i64arg;
}

/******************/

void print_help(void)
{
	fprintf(stderr, "Please refer to the help.txt file for the full list of command line options.\n");
	exit(EXIT_SUCCESS);
}

/******************/

/* Given an input line (presumed to be the first line of the mlucas.cfg or fermat.cfg file,
which is presumed to contain the program version that was used to generate said .cfg file),
returns nonzero (which should be taken as a proxy for TRUE) if the .cfg file needs to be updated
via a new round of timing tests, FALSE otherwise.
*/
int	cfgNeedsUpdating(char*in_line)
{
	/* For the foreseeable future, version numbers will be of form x.yz, with x a 2-digit integer,
	y an int having 3 digits or less, and z an optional alphabetic suffix. We choose these such that
	retiming is only necessary between versions that differ in x or in the leading digit of y,
	i.e. we'd need to regenerate .cfg files in going from version 16.9 to 17.0 or from 17.0 to 17.141,
	but not between versions 17.0x and 17.0g or between 17.1 and 17.141.

	This reduces the "retime?" decision to a simple comparison of the leading 4 characters of the version string.
	*/
	return STRNEQN(in_line, VERSION, 4);
}

/******************/

const char*returnMlucasErrCode(uint32 ierr)
{
	ASSERT(ierr < ERR_MAX, "Error code out of range!");
	return err_code[ierr-1];
}

void	printMlucasErrCode(uint32 ierr)
{
	/* Error type indicated by lowest byte: */
	uint32 i = (ierr & 0xff);

	if(i == 0)
		fprintf(stderr, "\n Return with error code 0 - no errors.\n");
	else if(i < ERR_MAX)
		fprintf(stderr, "\n Return with code %s\n", err_code[i-1]);
	else
		fprintf(stderr, "\n Return with unknown error error code %u - suggest running under a debugger.\n\n",(uint32)ierr);

	/* High bytes should only be nonzero if low byte == ERR_RUN_SELFTEST_FORLENGTH: */
	if((ierr>>8) != 0)
	{
		ASSERT(i==ERR_RUN_SELFTEST_FORLENGTH, "High bytes should only be nonzero if low byte == ERR_RUN_SELFTEST_FORLENGTH!");
	}
}

/******************/

/* Set of functions to Read/Write full-length residue data in bytewise format
from a primality-test or p-1 factoring run.

The format is an extension (to multiple moduli and assignment types) of the following -
EWM's notes added in {}:

	From: Richard Crandall <crandall>
	Date: Tue, 15 Dec 98 10:35:03 -0800
	To: Jason Stratos Papadopoulos <jasonp@glue.umd.edu>
	Subject: FORMAT definition for Pepin residue files
	Cc: F24 Millenium Group -- George Woltman <woltman@magicnet.magicnet.net>, Joe
	Buhler <jpb@reed.edu>, Richard Crandall <crandall@cascade.reed.edu>

	PEPIN RESIDUE FILE DEFINITION	{EWM - we use for any kind of full-modulus-length
									residue, be it from a primality test, p-1, or ECM run}

	(F_24 Millenium Group, advent T'giving 1998,
	first residues stored 5 Dec 1998 via REC's
	integer convolution checker.)

	ENTITIES IN FILE

	The (endian-platform-independent) residue file contains:

	{t}: EWM - add a leading byte whose value indicates the assignment type,
		i.e. stores the corresponding TEST_TYPE value as defined in Mdata.h
		(Even though those #defines get stored in long ints, we assume a max.
		value of 255.)

	{m}: EWM - add a leading byte whose value indicates the modulus type,
		i.e. stores the corresponding MODULUS_TYPE value as defined in Mdata.h
		(as with TEST_TYPE field, assume a max.	value of 255.)

	s:  The number of squarings of 3 {EWM: or LL iterations, or whatever; expanded to 8 bytes}
		that have been performed to generate said file, stored as four bytes,
		the low byte first.

		{EWM: extend s to 8 bytes, and store the following, depending on the case:

		NONFACTORIAL PRIMALITY/COMPOSITENESS TEST: number of iterations (typically
			modular squarings or some slight variant thereof) performed so far;

		P-1 FACTORIZATION:	If stage 1 using LR binary powering, #of bits processed
							so far (i.e. again = modular squarings performed).

							If stage 1 using prime-by-prime powering (e.g. if we are
							extending the depth of a previous p-1 run), most-recent
							stage 1 prime (or power thereof) processed. Assumes stage
							1 primes processed in strictly ascending order.

							If stage 2, most-recent stage 2 prime processed.

		ECM FACTORIZATION:	{TODO: Need to fill in details here}

		NOTE: For p-1 and ECM, the stage 1 and 2 bounds can be found in
		character form following the residue data (see below.)
		}

	r: For modulus N, exactly ceiling(log2(N)/8) residue bytes,
	   starting with low byte first in the file, ending with high byte. If balanced-digit
	   representation is used during the computation, this must be adjusted at file-write
	   time so that the whole residue is a standard, least non-negative
	   one (i.e. all bytes deemed positive).

	{c0}: EWM - an 8-byte unsigned integer checksum containing residue modulo 2^64, i.e.
		the least-significant 64 residue bits. This is all-but-useless as a residue integrity
		check, but is included for arcane code-historical reasons. This field may be repurposed
		at some later point.

	{c1}: EWM - a  5-byte unsigned integer checksum containing residue modulo 2^35 - 1.

	{c2}: EWM - a  5-byte unsigned integer checksum containing residue modulo 2^36 - 1.

	{a}: EWM - an arbitrary number (but of reasonably limited length, say <= 2^20 bytes)
		of ascii character bytes storing any other data needed to fully characterize the
		computation, e.g. stage bounds, possibly even a condensed run history, if more
		than one person has contributed pieces of the full computation.

EWM - added in v18:
	{kblocks}: 4 bytes to store the FFT length (in Kdoubles) actually
		used at time of savefile write. This is needed to support the code auto-switching
		to a larger-than-default FFT length during the run as a result of excessive ROE.

	{shift}: 8 bytes to store the circular bitwise residue-shift value at the time of the
		savefile write. Note that the unshifted residue is what is actually written to the file.

EWM - added in v19:
	For PRP-tests, second bytewise residue array and associated [c0,c1,c2] checksum triplet.
	The 8-byte {c0}-checksum field is repurposed here to store a circular bitwise residue-shift
	value which to apply to the G-check residues after reading it.

EWM - added in v20:
	{e1} A 4-byte field storing the cumulative number of occurrences of convolution-output fractional
	roundoff errors (ROE) >= 0.4375 (>= for LL, > for PRP)for the test in question.

	{e2} A 4-byte field storing the number of occurrences of Gerbicz-check errors for the test in question.
*/

/* Dec 2017: For Fermat case the Pepin primality test is indistinguishable from an Euler-PRP test and
e.g. cofactor-PRP tests start from Pepin residue. Thus for this modulus type these two
disparate-seeming test types are acceptable; add this little utility function to enforce that:
*/
int test_types_compatible(uint32 t1, uint32 t2)
{
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		return (t1 <= TEST_TYPE_MAX) && (t2 <= TEST_TYPE_MAX);
	else	// For Mersennes cofactor-PRP test must start from a PRP, not LL, residue, so continue to enforce strict equality of test types:
		return t1 == t2;
}

/*** READ: Assumes a valid file pointer has been gotten via a call of the form
fp = mlucas_fopen(RESTARTFILE,"rb");
***/
// Reads an [nbytes]-long LL|Pepin|PRP|p-1 residue and validates the associated mod(2^64,2^35-1,2^36-1) checksums:
int read_ppm1_residue(const uint32 nbytes, FILE*fp, uint8 arr_tmp[], uint64*Res64, uint64*Res35m1, uint64*Res36m1)
{
	const char func[] = "read_ppm1_residue";
	uint32 i,j;
	uint64 itmp64;
	const uint64 two35m1 = (uint64)0x7FFFFFFFFull, two36m1 = (uint64)0xFFFFFFFFFull;	/* 2^35,36-1 */

	*Res64 = *Res35m1 = *Res36m1 = 0ull;	// 0 value on return indicates failure of some kind

	i = fread(arr_tmp, sizeof(char), nbytes, fp);		/* Read bytewise residue...	*/
	if(i != nbytes)	{ sprintf(cbuf, "%s: Error reading bytewise residue array.\n",func)										; return 0; }
	if(ferror(fp))	{ sprintf(cbuf, "%s: Unknown Error reading bytewise residue array.\n",func)								; return 0; }
	if(feof(fp))	{ sprintf(cbuf, "%s: End-of-file encountered while attempting to read bytewise residue array.\n",func)	; return 0; }

	/* 8 bytes for Res64: */
	for(j = 0; j < 8; j++) {
		i = fgetc(fp);	*Res64 += (uint64)i << (8*j);
	}
	/* 5 bytes for Res35m1: */
	for(j = 0; j < 5; j++) {
		i = fgetc(fp);	*Res35m1 += (uint64)i << (8*j);
	}
	if(*Res35m1 > two35m1) { sprintf(cbuf, "%s: *Res35m1 should be <= 2^35",func); return 0; }
	/* 5 bytes for Res36m1: */
	for(j = 0; j < 5; j++) {
		i = fgetc(fp);	*Res36m1 += (uint64)i << (8*j);
	}
	if(*Res36m1 > two36m1) { sprintf(cbuf, "%s: *Res36m1 should be <= 2^36",func); return 0; }

	// Checksums: Check Res64,Res35m1,Res36m1 checksum vs the ones extracted from the residue byte-array read from the file:
	/* Since arr_tmp may hold previous-exponent data in its high bits (if previous p > current one), need to zero
	those out before calling any mi64 functions which treat arr_tmp as a 64-bit-limb array rather than a byte-array:
	*/
	j = 8 - (nbytes&7);	// nbytes%8 = #significant bytes in high limb; j = #bytes needing zeroing at the high end of the limb
	for(i = nbytes; i < nbytes+j; i++) arr_tmp[i] = 0;
	itmp64 = ((uint64*)arr_tmp)[0];
	if(*Res64 != itmp64) {
		sprintf(cbuf, "%s: On restart: Res64 checksum error! Got %llX, expected %llX\n"  ,func,itmp64,*Res64); return 0;
	}
	// For big-endian CPUs, casting byte-array to uint64* gives byte-reversed limbs, so use a direct bitwise mod:
  #ifdef USE_BIG_ENDIAN
	/*
	My original approach here was wrong - here was the reasoning, illustrated using the (mod 2^35-1) of the S-H mod pair:
		For low 64-bit limb x0, x0 % two35m1 = x0.lo35 + x0.hi29. Now say our current 64-bit limb starts at bit B, i.e.
		represents itmp64*2^B. That means we need to lcshift itmp64 by (B % 35) bits before doing the above 29|35 folding:
		x = b63...b0, 64-bit, we seek x*2^64 (mod 2^35-1).
	Here's why that doesn't work - consider the next-lowest 64-bit limb in our residue, call it x1:
		2^64  = 2^29 (mod 2^35-1), thus we seek x*2^29 (mod 2^35-1), which in bitwise terms is
		b63...b0[29 binary 0s] (mod 2^35-1), i.e. bits 35 thru 63+29 = 92 of b63...b0[29 binary 0s] get folded back by adding to low bits:
									34.33.32.31.30.29.28.27.26.25.24.23.22.21.20.19.18.17.16.15.14.13.12.11.10.09.08.07.06.05.04.03.02.01.00
			b5...b0[29 binary 0s]	 5. 4. 3. 2. 1. 0.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
		+	b40...b6				40.39.38.37.36.35.34.33.32.31.30.29.28.27.26.25.24.23.22.21.20.19.18.17.16.15.14.13.12.11.10.09.08.07.06
		+	b63...b41				                                       63.62.61.60.59.58.57.56.55.54.43.52.51.49.48.47.46.45.44.43.42.41
		But the circular-shift scheme instead gives
		[b34...b0][b63...b35], folding that on the 35-bit boundary gives
			[b5...b0][b63...b35]
		+	[b34...b6], i.e. bits 35-40 and 41-63 end up misaligned.
		What we really want is to start with the current value of bmod35, 29, and do like so:
		bmod35 = 29: (x << 29) & two35m1
					+(x >> (35-29)) & two35m1
					+(x >> (70-29)) & two35m1
	*/
	int bmod35 = 0, bmod36 = 0, rshift;
	uint64 rmod35 = 0ull, rmod36 = 0ull, cmod64;
	for(i = 0; i < nbytes; i += 8) {
		// Assemble 64-bit limb from byte octet:
		itmp64 = ((uint64)arr_tmp[i]<<56)+((uint64)arr_tmp[i+1]<<48)+((uint64)arr_tmp[i+2]<<40)+((uint64)arr_tmp[i+3]<<32)
				+((uint64)arr_tmp[i+4]<<24)+((uint64)arr_tmp[i+5]<<16)+((uint64)arr_tmp[i+6]<<8)+(uint64)arr_tmp[i+7];
		cmod64 = (itmp64 << bmod35) & two35m1; rshift = 35 - bmod35;
		while(rshift < 64) { cmod64 += (itmp64 >> rshift) & two35m1; rshift += 35; }
		rmod35 += cmod64;
		cmod64 = (itmp64 << bmod36) & two36m1; rshift = 36 - bmod36;
		while(rshift < 64) { cmod64 += (itmp64 >> rshift) & two36m1; rshift += 36; }
		rmod36 += cmod64;
		MOD_ADD64(bmod35,29,35,bmod35); MOD_ADD64(bmod36,28,36,bmod36);	// bmod35|36 += 29|28 (mod 35|36)
	}
	rmod35 = (rmod35 & two35m1) + (rmod35 >> 35); rmod36 = (rmod36 & two36m1) + (rmod36 >> 36);	// And do a final pair of folds to get mods
	if(*Res35m1 != rmod35)	{ sprintf(cbuf, "%s: On restart: Res35m1 checksum error! Got %llX, expected %llX\n",func,rmod35,*Res35m1); return 0; }
	if(*Res36m1 != rmod36)	{ sprintf(cbuf, "%s: On restart: Res36m1 checksum error! Got %llX, expected %llX\n",func,rmod36,*Res36m1); return 0; }
  #else
	i = (nbytes+7)>>3;	// # of 64-bit limbs
	itmp64 = mi64_div_by_scalar64((uint64*)arr_tmp,two35m1,i,0x0);
	if(*Res35m1 != itmp64) {
		sprintf(cbuf, "%s: On restart: Res35m1 checksum error! Got %llX, expected %llX\n",func,itmp64,*Res35m1); return 0;
	}
	itmp64 = mi64_div_by_scalar64((uint64*)arr_tmp,two36m1,i,0x0);
	if(*Res36m1 != itmp64) {
		sprintf(cbuf, "%s: On restart: Res36m1 checksum error! Got %llX, expected %llX\n",func,itmp64,*Res36m1); return 0;
	}
  #endif
	return 1;
}

// Returns 1 on successful read, 0 otherwise:
// v19: For PRP-tests, also write a second Gerbicz-check residue array [arr2] and associated S-H checksum triplet [i1,i2,i3]
// v20: Distributed deep p-1 S2 may use B2 >= 2^32, so make ilo a uint64-ptr; add filename arg since S2 appends '.s2' to RESTARTFILE:
int read_ppm1_savefiles(const char*fname, uint64 p, uint32*kblocks, FILE*fp, uint64*ilo,
	uint8 arr1[], uint64*Res64, uint64*Res35m1, uint64*Res36m1,
	uint8 arr2[], uint64*i1   , uint64*i2     , uint64*i3     )
{
	const char func[] = "read_ppm1_savefiles";
	uint32 i,j,k,len,nbytes = 0,nerr;
	uint64 itmp64, nsquares = 0ull, *avec = (uint64*)arr1, *bvec = (uint64*)arr2, exp[4],pow[4],rem[4];
	uint128 ui128,vi128; uint192 ui192,vi192; uint256 ui256,vi256;	// Fixed-length 2/3/4-word ints for stashing results of multiword modexp.
	*Res64 = 0ull;	// 0 value on return indicates failure of some kind
	mi64_clear(pow,4); mi64_clear(rem,4);
	ASSERT(arr1 != 0x0, "Null arr1 pointer!");
	if(!file_valid(fp)) {
		sprintf(cbuf, "%s: File pointer invalid for read!\n",func);	ASSERT(0, cbuf);
	}
	fprintf(stderr, " INFO: restart file %s found...reading...\n",fname);
	/* t: */
	i = fgetc(fp);
	if(!test_types_compatible(i, TEST_TYPE)) {
		sprintf(cbuf, "%s: TEST_TYPE != fgetc(fp)\n",func);
		return 0;
	}
	/* m: */
	if((i = fgetc(fp)) != MODULUS_TYPE) {
		// For some reason, this fubared in my rerun-final-F25-iterations-from-33.55m (fgetc = 176, MODULUS_TYPE = 3)
		// but residue OK, so emit error msg but allow execution past it:
		sprintf(cbuf, "ERROR: %s: MODULUS_TYPE != fgetc(fp)\n",func);
		return 0;
	}
	/* s: */
 	for(j = 0; j < 8; j++) {
		i = fgetc(fp);	nsquares += (uint64)i << (8*j);
	}
	// v20: E.g. distributed deep p-1 S2 may use B2 >= 2^32: Only allow nsquares >= 2^32 if it's an S2 restart:
	if(TEST_TYPE == TEST_TYPE_PM1) {
		if(strstr(fname, ".s2")) {
			if(nsquares > 0xFFFFFFFFull)
				ASSERT(B2_start <= nsquares, "P-1 stage 2 restart requires (B2_start in worktodo assignment) <= (savefile nsquares field)!");
		} else {	// It's a stage 1 restart:
			ASSERT(nsquares <= 0xFFFFFFFFull && nsquares < 1.5*(double)B1, "P-1 stage 1 restart: savefile nsquares value out of bounds!");
		}
		// If S2 restart and (nsquares > B2_start), read the ensuing S2 interim residue; if (nsquares == B2_start)
		// it means S2 started but was aborted for some reason before writing an interim S2 residue. That will set
		// *ilo = B2_start below and trigger a 'return 0' following the ensuing read_ppm1_residue() call, which
		// the S2 code interprets as "start stage 2 from B2_start."
	} else {	// For primality-tests, make sure nsquares < 2^32 and copy to ilo:
		if(nsquares > p) {	// v21: change from >= p to > p, since Mersenne-PRP restart-to-check-CF will have nsquares == p:
			sprintf(cbuf,"%s: nsquares = %llu out of range, should be < p = %llu\n",func, nsquares, p);
			return 0;
		} else if(nsquares > 0xFFFFFFFFull) {
			sprintf(cbuf,"%s: nsquares = %llu out of range, current limit = 2^32-1.\n",func, nsquares);
			return 0;
		}
	}
	*ilo = nsquares;

	/* Set the expected number of residue bytes, depending on the modulus: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		nbytes = (p + 7)/8;
		TRANSFORM_TYPE = REAL_WRAPPER;
	} else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
		sprintf(cbuf, "%s: MODULUS_TYPE_FERMAT but (p mod 8) != 0",func);
		ASSERT((p & 7) == 0,cbuf);
		nbytes = (p>>3) + 1;
		TRANSFORM_TYPE = RIGHT_ANGLE;
	}

	i = read_ppm1_residue(nbytes, fp, arr1, Res64,Res35m1,Res36m1);
	if(!i) return 0;

	/* For PRP-tests of N with known-factors (e.g. prelude to a cofactor-PRP postprocessing step),
	R. Gerbicz [post #79 of https://mersenneforum.org/showthread.php?t=18748&page=8] suggests, using
	the specific case of the (m)th Fermat number F_m to illustrate:

		Have you checked that:
			(3^(2^(2^m)) mod F_m) mod q = 3^(2^(2^m)) mod q
		is true for the stored interim residues?

		Since you have 3^(2^(2^m)) mod F_n then to get the left side is trivial and you can get quickly
		the right side, because q is small , the check is cheap. Though it is not that very-very strong
		check, but the two known factors could find a mismatch with probability = 1-(1e-10) which is
		still quite impressive. Note that we don't get the strength of 1/(q-1) or 1/znorder(Mod(3,q))
		because the order is divisible by a "large" power of two, and that will be lost in the iterated
		squarings of a wrong residue.

	(For Pepin-test residue, we'd actually check (3^(2^(2^m - 1)) mod F_m) mod q = 3^(2^(2^m - 1)) mod q.)
	More generally, for the powermod residue R = a^pow (mod N) where N has prime factor q, we check whether

		a^pow (mod N) = a^pow (mod q).

	If pow > q, the RHS can be simplified using Fermat's little theorem, a^(q-1) == 1 (mod q), thus
	a^pow == a^pow' (mod q), where pow' = pow (mod q-1).

	Ex: For the final Pépin-test residue R for F30 (N = 2^2^30 + 1, pow = (N-1)/2 = 2^(2^30 - 1)), with
	2 known small prime factors p = 640126220763137 and q = 1095981164658689, we get pow' =

		p: pow' = 2^(2^30 - 1) (mod p-1) = 433448099512320, R == 367500396407933 (mod p)
			[bc code: n = 2^30-1; pow = modpow_lr2(n,p-1); modpow_lr(3,pow,p)]
		q: pow' = 2^(2^30 - 1) (mod q-1) =  96576634617856, R ==  95971534699540 (mod q).

	Note: This FLT-shortcut does not work for products of prime factors, since even though products of
	Fermat-number factors are base-2 Fermat pseudoprimes, they are not base-3 such. For example, for
	F25, the product of 3 known factors: f = 11528509095030405010053836714337362416244581466113; n = 2^25-1; pow = modpow_lr2(n,f-1); modpow_lr(3,pow,f)
	gives:
			pow = 8069349577210892231936784610857835297469527752704
		3^pow%f =  820422345454042236417736966133820291898601368899 ,
	but my mi64_div routine gives the correct remainder 3398198221728141144755902061331087657421092620571, which clearly differs from 3^pow%f result.
	(Also added CRT routine to nt_utils.txt, which takes the 3 remainders mod the known prime factors and confirms the DIV result.)
	*/
	if(TEST_TYPE == TEST_TYPE_PRP) {
		len = (nbytes+7)>>3; j = p&63; itmp64 = avec[len-1];	ASSERT((itmp64 >> j) == 0ull, "High limb of residue array1 does not have upper bits cleared!");
		for(i = 0; KNOWN_FACTORS[i] != 0ull; i += 4) {
			j = mi64_getlen(KNOWN_FACTORS+i,4);	// j = number of nonzero limbs in curr_fac (alloc 4 limbs per in KNOWN_FACTORS[])
			sprintf(cstr,"Computing %llu-squaring residue R (mod known prime q = %s)\n",nsquares,&cbuf[convert_mi64_base10_char(cbuf, KNOWN_FACTORS+i, j, 0)] ); mlucas_fprint(cstr,1);
			mi64_div(avec,KNOWN_FACTORS+i, len,j, 0x0,rem);	// R (mod p) returned in rem[]
			k = mi64_getlen(rem,4);	// j = number of nonzero limbs in remainder
			sprintf(cstr,"\tA: R == %s (mod q)\n",&cbuf[convert_mi64_base10_char(cbuf, rem, k, 0)] ); mlucas_fprint(cstr,1);
			if(j == 1) {
				exp[0] = twopmmodq64(nsquares,KNOWN_FACTORS[i]-1);	// pow' = 2^nsquares (mod p-1)
			} else if(j == 2) {
				vi128.d0 = nsquares; vi128.d1 = 0ull;
				ui128.d0 = KNOWN_FACTORS[i]-1; ui128.d1 = KNOWN_FACTORS[i+1];
				ui128 = twopmmodq128(vi128,ui128);	// pow' = 2^nsquares (mod p-1)
				exp[0] = ui128.d0; exp[1] = ui128.d1;
			} else if(j == 3) {
				vi192.d0 = nsquares; vi192.d1 = vi192.d2 = 0ull;
				ui192.d0 = KNOWN_FACTORS[i]-1; ui192.d1 = KNOWN_FACTORS[i+1]; ui192.d2 = KNOWN_FACTORS[i+2];
				ui192 = twopmmodq192(vi192,ui192);	// pow' = 2^nsquares (mod p-1)
				exp[0] = ui192.d0; exp[1] = ui192.d1; exp[2] = ui192.d2;
			} else if(j == 4) {
				vi256.d0 = nsquares; vi256.d1 = vi256.d2 = vi256.d3 = 0ull;
				ui256.d0 = KNOWN_FACTORS[i]-1; ui256.d1 = KNOWN_FACTORS[i+1]; ui256.d2 = KNOWN_FACTORS[i+2]; ui256.d3 = KNOWN_FACTORS[i+3];
				ui256 = twopmmodq256(vi256,ui256);	// pow' = 2^nsquares (mod p-1)
				exp[0] = ui256.d0; exp[1] = ui256.d1; exp[2] = ui256.d2; exp[3] = ui256.d3;
			} else
				ASSERT(0, "Only known-factors < 2^256 supported!");
			// Raise PRP base (usually but not always 3) to the just-computed power; result in 4-limb local-array pow[]:
			mi64_scalar_modpow_lr(PRP_BASE, exp, KNOWN_FACTORS+i, j, pow);
			sprintf(cstr,"\tB: R == %s (mod q)\n",&cbuf[convert_mi64_base10_char(cbuf, pow, j, 0)] ); mlucas_fprint(cstr,1);
			if (mi64_getlen(pow,4) != k || !mi64_cmp_eq(pow,rem,k)) {
				snprintf(cbuf,STR_MAX_LEN,"Full-residue == %u^nsquares (mod q) check fails!", PRP_BASE); mlucas_fprint(cbuf,0);
				ASSERT(0, cbuf);
			}
		}
	}
#if 0
/* Above pow' = 2^(2^30 - 1) (mod p-1) = 433448099512320 = 2^35.3.5.29^2 example:
My 2^-p (mod q) code needs q odd ... we note q = (p-1) = 2^32.103.1447 = 2^32.149041, so first try positive-power version:
q = 2^32.qodd, so compute r' = 2^(p-32) (mod qodd = 149041) = 100920 = 2^3.3.5.29^2, whence r = 2^32.r' = 433448099512320
The inverse of r' (mod qodd) = 38455, and that is just the output of the negative-powering r'' = 2^(-(p-32)) (mod 149041).
Thus r'.r'' = 100920.38455 == 1 (mod 149041), BUT (2^32.r').(2^32.r'') = 433448099512320.165162967367680 !== 1 (mod p-1).
Thus if we use a negative-power algo, to recover 2^p (mod q = 2^k.qodd):
1. Compute negative-power r'' = 2^(-(p-k)) (mod qodd);
2. Compute mod-inverse r' = (r'')^-1 (mod qodd);
3. Then 2^k.r' = 2^p (mod q = 2^k.qodd). */
#endif

	// FFT length in K (3 bytes) - first added this in v18:
	i = *kblocks = 0;
	for(j = 0; j < 3 && i != EOF; j++) {
		i = fgetc(fp);	*kblocks += (uint64)i << (8*j);
	}
	if(i == EOF) {
		*kblocks = 0;
		sprintf(cbuf,"%s: Hit EOF in read of FFT-kblocks ... assuming a pre-v18 savefile.\n",func); fprintf(stderr,"%s", cbuf); return 1;
	}
	/* May 2018: 8 bytes for circular-shift to apply to the (unshifted) residue read from the file: */
	i = 0; RES_SHIFT = 0ull;
	for(j = 0; j < 8 && i != EOF; j++) {
		i = fgetc(fp);	RES_SHIFT += (uint64)i << (8*j);
	}
	if(i == EOF) {
		RES_SHIFT = 0ull;
		sprintf(cbuf,"%s: Hit EOF in read of FFT-kblocks ... assuming a pre-v18 savefile.\n",func); fprintf(stderr,"%s", cbuf); return 1;
	}

  // v19: For PRP-tests, also read a second Gerbicz-check residue array [arr2] and associated S-H checksum triplet [i1,i2,i3]:
  if(DO_GCHECK) {	// v21: Change to key off DO_GCHECK, to allow Fermat-mod Pepin-tests to use the Gerbicz check, too
	ASSERT(arr2 != 0x0, "Null arr2 pointer!");
	PRP_BASE = 0ull;
	for(j = 0; j < 4; j++) {
		i = fgetc(fp);	PRP_BASE += i << (8*j);
	}
	i = read_ppm1_residue(nbytes, fp, arr2, i1,i2,i3);
	if(!i) return 0;
	// G-check residues all need to be clshifted by residue-shift count at the ITERS_BETWEEN_GCHECK_UPDATESth PRP-test iteration:
	GCHECK_SHIFT = 0ull;
	for(j = 0; j < 8; j++) {
		i = fgetc(fp);	GCHECK_SHIFT += (uint64)i << (8*j);
	}
  }

  // v20: Read cumulative #errs for ROE >= 0.4375 (>= for LL, > for PRP) and Gerbicz-check for the test in question.
  // If it's a restart from a v19 run, emit a warning and simply start these from their main()-init values of 0:
	nerr = 0ull;
	for(j = 0; j < 4; j++) {
		i = fgetc(fp);
		if(i == EOF) {
			if(!j) {
				sprintf(cbuf, "%s: Restart from v19 savefile - will start tracking #errors encountered at this point.\n",func);
				fprintf(stderr,"%s", cbuf);
				return 1;
			} else {	// If at least the first of the 3 bytes exists, all 3 had better be there:
				sprintf(cbuf, "%s: Expected 4 nerr bytes!",func);
				fprintf(stderr,"%s", cbuf);
				return 0;
			}
		}
		nerr += i << (8*j);
	}
	NERR_ROE = MAX(nerr,NERR_ROE);	// If restart-from-savefile as result of hitting an ROE, preserve the runtime-incremented value of NERR_ROE:
	// Similar handling for G-check error count:
	nerr = 0ull;
	for(j = 0; j < 4; j++) {
		i = fgetc(fp);	nerr += i << (8*j);
	}
	NERR_GCHECK = MAX(nerr,NERR_GCHECK);
	/* Don't deallocate arr1 here, since we'll need it later for savefile writes. */
	return 1;
}

/* WRITE: Assumes the following:
[1] a valid file pointer has been gotten via a call of the form fp = mlucas_fopen(RESTARTFILE,"wb");
[2] any circular shift stored in the global RES_SHIFT has been removed in the preceding call to convert_res_FP_bytewise.
v19: For PRP-tests, also write a second Gerbicz-check residue array [arr2] and associated S-H checksum triplet [i1,i2,i3]
*/
void write_ppm1_residue(const uint32 nbytes, FILE*fp, const uint8 arr_tmp[], const uint64 Res64, const uint64 Res35m1, const uint64 Res36m1)
{
	const char func[] = "write_ppm1_residue";
	uint32 i;
	/* Write bytewise residue r... */
	i = fwrite(arr_tmp, sizeof(char), nbytes, fp);
	if(i != nbytes) {
		fclose(fp); fp = 0x0;
		snprintf(cbuf,STR_MAX_LEN*2,"%s: Error writing residue to restart file.\n",func);
		mlucas_fprint(cbuf,0);	ASSERT(0,cbuf);
	}
	/* ...and checksums:	*/
	/* Res64: */
	for(i = 0; i < 64; i+=8)
		fputc((int)(Res64 >> i) & 0xff, fp);
	/* Res35m1: */
	for(i = 0; i < 40; i+=8)
		fputc((int)(Res35m1 >> i) & 0xff, fp);
	/* Res36m1: */
	for(i = 0; i < 40; i+=8)
		fputc((int)(Res36m1 >> i) & 0xff, fp);
}

// v20: E.g. distributed deep p-1 S2 may use B2 >= 2^32, so make ihi a uint64; add filename arg since S2 appends '.s2' to RESTARTFILE:
void write_ppm1_savefiles(const char*fname, uint64 p, int n, FILE*fp, uint64 ihi,
	uint8 arr1[], uint64 Res64, uint64 Res35m1, uint64 Res36m1,
	uint8 arr2[], uint64 i1   , uint64 i2     , uint64 i3     )
{
	uint32 i,kblocks,nbytes = 0;
	ASSERT(file_valid(fp),"write_ppm1_savefiles: File pointer invalid for write!");
	// Make sure n is a proper (unpadded) FFT-length, i.e. is a multiple of 1K:
	kblocks = (n >> 10);
	ASSERT(n == (kblocks << 10),"Not a proper unpadded FFT length");

	/* See the function read_ppm1_savefiles() for the file format here: */
	/* t: */
	fputc(TEST_TYPE, fp);
	/* m: */
	fputc(MODULUS_TYPE, fp);
	/* s: */
	for(i = 0; i < 64; i+=8)
		fputc((ihi >> i) & 0xff, fp);

	/* Set the expected number of residue bytes, depending on the modulus: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		nbytes = (p + 7)/8;
		TRANSFORM_TYPE = REAL_WRAPPER;
	} else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
		ASSERT(p % 8 == 0,"write_ppm1_savefiles: p % 8 == 0");
		nbytes = (p/8) + 1;	// We don't expect > p bits except in the highly unlikely case of a prime-Fermat Pepin-test result
		TRANSFORM_TYPE = RIGHT_ANGLE;
	}

	write_ppm1_residue(nbytes, fp, arr1, Res64,Res35m1,Res36m1);

	// v18: FFT length in K (3 bytes):
	for(i = 0; i < 24; i+=8)
		fputc((kblocks >> i) & 0xff, fp);
	// v18: circular-shift to apply to the (unshifted) residue read from the file (8 bytes):
	for(i = 0; i < 64; i+=8)
		fputc((RES_SHIFT >> i) & 0xff, fp);

  // v19: For PRP-tests, also write a second Gerbicz-check residue array [arr2] and associated S-H checksum triplet [i1,i2,i3]:
  if(DO_GCHECK) {	// v21: Change to key off DO_GCHECK, to allow Fermat-mod Pepin-tests to use the Gerbicz check, too
	for(i = 0; i < 32; i+=8)
		fputc((PRP_BASE >> i) & 0xff, fp);
	write_ppm1_residue(nbytes, fp, arr2, i1,i2,i3);
	// G-check residues all need to be clshifted by residue-shift count at the ITERS_BETWEEN_GCHECK_UPDATESth PRP-test iteration:
	for(i = 0; i < 64; i+=8)
		fputc((GCHECK_SHIFT >> i) & 0xff, fp);
  }
	// v20: Write cumulative #errs for ROE >= 0.4375 (>= for LL, > for PRP) and Gerbicz-check for the test in question:
	for(i = 0; i < 32; i+=8)
		fputc((NERR_ROE >> i) & 0xff, fp);
	for(i = 0; i < 32; i+=8)
		fputc((NERR_GCHECK >> i) & 0xff, fp);
}

/*********************/

/*
Function to take an n-digit Mersenne or Fermat-mod residue in bytewise savefile form,
apply the required circular shift read into the global RES_SHIFT during the preceding
bytewise-savefile read and convert it to balanced-digit floating-point form.

In the Mersenne-mod case the residue digits are stored consecutively in the a[] array.

In the Fermat-mod case the digits are arranged in (j,j+n/2) (i.e. right-angle transform) order.
*/
int 	convert_res_bytewise_FP(const uint8 ui64_arr_in[], double a[], int n, const uint64 p)
{
	uint32 nbytes;
	uint64 nbits;
	int bimodn,curr_char,findex,ii,j = 0,j1 = 0,k,pass,rbits;
	int bw,sw,bits[2];
	uint64 base[2];	/* Assume base may be > 2^32 (e.g. for mixed FFT/FGT)
						but < 2^53, i.e. fits in a double */
	int64 itmp,cy;
	uint64 curr_word, curr_wd64;
	int pow2_fft;

	ASSERT(MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	ASSERT(TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");

	/* Set the number of residue bytes, which is the same for Mersenne (2^p-1) and Fermat-mod (2^p+1, with p = 2^findex)
	despite the fact the latter can formally be as large as 2^p, since only ever hit that if it`s the last residue of
	a Pepin test and the number hqppens to be prime. (We would love for that exception to break some other ASSERTion in the code): */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		ASSERT(TRANSFORM_TYPE == REAL_WRAPPER,"convert_res_bytewise_FP: TRANSFORM_TYPE == REAL_WRAPPER");
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(TRANSFORM_TYPE == RIGHT_ANGLE,"convert_res_bytewise_FP: TRANSFORM_TYPE == RIGHT_ANGLE");
		/* If Fermat number, make sure exponent a power of 2: */
		findex = trailz64(p);
		ASSERT((p >> findex) == 1,"convert_res_bytewise_FP: (p >> findex) == 1");

		ASSERT(p % 8 == 0,"convert_res_bytewise_FP: p % 8 == 0");
	}
	nbytes = (p + 7)/8;
	// Apply the circular shift:
	uint32 is_ferm = (MODULUS_TYPE == MODULUS_TYPE_FERMAT);
	mi64_shlc((uint64*)ui64_arr_in, (uint64*)ui64_arr_in, p, RES_SHIFT, (p+63+is_ferm)>>6, is_ferm);

	/* Vector length a power of 2? */
	pow2_fft = (n >> trailz32(n)) == 1;

	bits[0] = p/n;		ASSERT(bits[0] > 1,"convert_res_bytewise_FP: bits[0] > 1");
	base[0] = 1 << bits[0];

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT && pow2_fft == TRUE)
		bits[1] =     bits[0];
	else
		bits[1] = 1 + bits[0];
	base[1] = 1 << bits[1];

	bw = p%n;	/* cf. mers_mod_square.c	*/
	sw = n - bw;

	/*...Now form the SH residues, converting to positive-digit form along the way...	*/

	curr_char = 0;	/* Current byte to be read from the ui64_arr_in array */
	nbits = 0;		/* Total bits accumulated so far in the residue	*/
	rbits = 0;		/* # of bits left in our 64-bit buffer after processing of previous word	*/
	curr_wd64 = 0;	/*      bits left in our 64-bit buffer after processing of previous word	*/
	cy = 0;

	/* For right-angle transform, need to process odd and even-index elements
	separately, via a pair of stride-2 passes through the array:
	*/
	if(TRANSFORM_TYPE == REAL_WRAPPER)
	{
		bimodn = 0;
		ii = 0;		/* Index into the BASE array. */
		/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
		if(bw > 0)
			ii = 1;

		for(j = 0; j < n; j++)
		{
		#ifdef USE_AVX512
			j1 = (j & mask03) + br16[j&15];
		#elif defined(USE_AVX)
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			/* If rbits (# of bits left in current 64-bit window of the bytewise residue)
			is less than bits[ii] (number we need to fill the current digit of a[]), write
			the bits we have to curr_word, grab the next 8 bytes of ui64_arr_in[] and finish
			the write of the remaining (bits[ii] - rbits) of curr_word:
			*/
			if(rbits < bits[ii])
			{
				itmp = curr_wd64;
				ASSERT(itmp < (1ull<<rbits),"convert_res_bytewise_FP: itmp >= 2^rbits!");

				/* Now grab the next 64 bits of the bytewise residue... */
				curr_wd64 = 0;
				for(k = 0; k < 8; k++)
				{
					curr_wd64 += (uint64)ui64_arr_in[curr_char++] << (k<<3);	/* left-shift current residue byte k*8 bits and accumulate */
					if(curr_char == nbytes)
						break;
				}
				/* ...and use the LS (bits[ii] - rbits) of the just-grabbed 64
				to fill in the high part of the bits[ii] bits of itmp: */
				curr_word = (curr_wd64 << (64 - (bits[ii] - rbits)));	/* Off-shift everything but the bits we need... */
				curr_word = (curr_word >> (64 -  bits[ii]));	/* Now right-shift so LSB of the bits we need goes in <rbits> slot... */
				itmp +=  curr_word;
				curr_wd64 >>= (bits[ii] - rbits);
				rbits = 64 - (bits[ii] - rbits);	/* # of bits now left in curr_wd64 */
			}
			else	/* Grab the LS bits[ii] bits of curr_wd64 and stick them in itmp */
			{
				curr_word = (curr_wd64 << (64 - bits[ii]));			/* Off-shift everything but the bits we need... */
				itmp = (curr_word >> (64 - bits[ii]));
				curr_wd64 >>= bits[ii];
				rbits -= bits[ii];					/* # of bits now left in curr_wd64 */
			}

			itmp += cy;		/* Add in any carry from the previous digit
							and normalize result so |current digit| <= base[ii]/2. */
			if(itmp > 0) {
				cy = itmp >> bits[ii];
				itmp -= (cy << bits[ii]);
				if(itmp > (base[ii]>>1))
				{
					itmp -= base[ii];
					cy += 1;
				}
			} else {
				cy = 0;
			}
			a[j1]= (double)itmp;
			nbits += bits[ii];
			bimodn += bw;
			if(bimodn >= n) bimodn -= n;
			ii = (uint32)(sw - bimodn) >> 31;
		}
	}
	else
	{
	  for(pass = 0; pass <=1; pass++)
	  {
		bimodn = n;

		for(j = pass; j < n; j += 2)
		{
		#ifdef USE_AVX512
			j1 = (j & mask03) + br16[j&15];
		#elif defined(USE_AVX)
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			ii = (bimodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */
			bimodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */
			bimodn += ( ((int)ii-1) & n);		/*       add 0 if a bigword,   N if a smallword */

			/* If rbits (# of bits left in current 64-bit window of the bytewise residue)
			is less than bits[ii] (number we need to fill the current digit of a[]), write
			the bits we have to curr_word, grab the next 8 bytes of ui64_arr_in[] and finish
			the write of the remaining (bits[ii] - rbits) of curr_word:
			*/
			if(rbits < bits[ii])
			{
				itmp = curr_wd64;
				ASSERT(itmp < (1<<rbits),"convert_res_bytewise_FP: itmp >= 2^rbits!");

				/* Now grab the next 64 bits of the bytewise residue... */
				curr_wd64 = 0;
				for(k = 0; k < 8; k++)
				{
					curr_wd64 += (uint64)ui64_arr_in[curr_char++] << (k<<3);	/* left-shift current residue byte k*8 bits and accumulate */
				}
				/* ...and use the LS (bits[ii] - rbits) of the just-grabbed 64
				to fill in the high part of the bits[ii] bits of itmp: */
				curr_word = (curr_wd64 << (64 - (bits[ii] - rbits)));	/* Off-shift everything but the bits we need... */
				curr_word = (curr_word >> (64 -  bits[ii]));	/* Now right-shift so LSB of the bits we need goes in <rbits> slot... */
				itmp +=  curr_word;
				curr_wd64 >>= (bits[ii] - rbits);
				rbits = 64 - (bits[ii] - rbits);	/* # of bits now left in curr_wd64 */
			}
			else	/* Grab the LS bits[ii] bits of curr_wd64 and stick them in itmp */
			{
				curr_word = (curr_wd64 << (64 - bits[ii]));			/* Off-shift everything but the bits we need... */
				itmp = (curr_word >> (64 - bits[ii]));
				curr_wd64 >>= bits[ii];
				rbits -= bits[ii];					/* # of bits now left in curr_wd64 */
			}

			itmp += cy;		/* Add in any carry from the previous digit
							and normalize result so |current digit| <= base[ii]/2. */
			if(itmp > 0) {
				cy = itmp >> bits[ii];
				itmp -= (cy << bits[ii]);
				if(itmp > (base[ii]>>1))
				{
					itmp -= base[ii];
					cy += 1;
				}
			} else {
				cy = 0;
			}
			a[j1]= (double)itmp;
			nbits += bits[ii];
		}
	  }
	}

	ASSERT(curr_char == nbytes, "convert_res_bytewise_FP: curr_char == (p+7)/8");
	ASSERT(nbits == p    ,"convert_res_bytewise_FP: nbits == p    ");
	ASSERT(curr_wd64 == 0,"convert_res_bytewise_FP: curr_word == 0");

	/*
	Fold any carryout from the conversion to balanced-representation form
	into the LS word of the residue (e.g. if cy = 1, this amounts to subtracting
	the modulus from the positive-digit form to get the balanced-digit form):
	*/
	/* Should have carryout of +1 Iff MS word < 0; otherwise expect 0 carry: */
	if(cy && (a[j1] >= 0 || cy != +1))
	{
		sprintf(cbuf, "convert_res_bytewise_FP: Illegal combination of nonzero carry = %lld, most sig. word = %20.4f\n", cy, a[j]);
		ASSERT(0, cbuf);
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		a[0] += cy;
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		a[0] -= cy;
	else
		ASSERT(0,"Illegal modulus type!");
	return TRUE;
}

/*********************/

/*
Function to take an n-digit Mersenne or Fermat-mod residue in balanced-digit
floating-point form, convert it to bytewise form, un-apply circular shift
stored in the global RES_SHIFT, and (optionally) generate Selfridge-Hurwitz checksums from the result.

In the Mersenne-mod case the residue digits are assumed to be stored consecutively in the a[] array.

In the Fermat-mod case digits assumed arranged in (j,j+n/2), i.e. right-angle transform, order.

***NOTE:*** To convert a generic multiword int (i.e. not w.r.to a particular
modulus) from/to balanced-digit fixed-base floating-point form and uint64
form, use the mi64_cvt_double_uint64() and mi64_cvt_uint64_double() functions in mi64.c .
*/
void	convert_res_FP_bytewise(const double a[], uint8 ui64_arr_out[], int n, const uint64 p, uint64*Res64, uint64*Res35m1, uint64*Res36m1)
{
	const char func[] = "convert_res_FP_bytewise";
	int bimodn,curr_bits,curr_char,cy,findex,ii,j,j1,k,pass,rbits,msw_lt0,bw,sw,bits[2],pow2_fft;
	uint64 nbits,curr_wd64,base[2];	/* Assume base may be > 2^32 (e.g. for mixed FFT/FGT) but < 2^53, i.e. fits in a double */
	double atmp;
	int64 itmp;
	const uint64 two35m1 = (uint64)0x00000007FFFFFFFFull, two36m1 = (uint64)0x0000000FFFFFFFFFull;	/* 2^35,36-1 */
	uint64*u64_ptr = (uint64*)ui64_arr_out;

	ASSERT(MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");
	ASSERT(TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");

	/* If Fermat number, make sure exponent a power of 2: */
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(TRANSFORM_TYPE == RIGHT_ANGLE,"convert_res_FP_bytewise: TRANSFORM_TYPE == RIGHT_ANGLE");
		findex = trailz64(p);
		ASSERT((p >> findex) == 1,"convert_res_FP_bytewise: (p >> findex) == 1");
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		ASSERT(TRANSFORM_TYPE == REAL_WRAPPER,"convert_res_FP_bytewise: TRANSFORM_TYPE == REAL_WRAPPER");
	else
		ASSERT(0,"Illegal modulus type!");

	/* Vector length a power of 2? */
	pow2_fft = (n >> trailz32(n)) == 1;

	bits[0] = p/n;
	base[0] = 1 << bits[0];

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT && pow2_fft == TRUE)
		bits[1] =     bits[0];
	else
		bits[1] = 1 + bits[0];
	base[1] = 1 << bits[1];

	bw = p%n;	/* cf. mers_mod_square.c	*/
	sw = n - bw;

	/*
	If most-significant digit in the balanced-representation form is < 0, add the modulus to the residue.
	For Mersenne (2^p-1) and Fermat (2^p+1) moduli, combine this with the normalize-to-nonnegative-digit
	step (which we do in any event) by simply initializing the carry into the latter to -1 or +1, respectively.
	In this case we expect the carryout of the normalization loop to = -1, indicating that the MS word has
	been accordingly normalized during the loop - add an assertion check to that effect.
	*/
	cy=0;		/* init carry.	*/
	msw_lt0 = 0;
	atmp = 0.0;	// To make sure the find-most-significant-nonzero-element check below actually finds a nonzero element
	for(j = n-1; j >= 0; j -= TRANSFORM_TYPE)
	{
	#ifdef USE_AVX512
		j1 = (j & mask03) + br16[j&15];
	#elif defined(USE_AVX)
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		atmp = a[j1];
		if(atmp != 0.0) {
		//	fprintf(stderr,"convert_res_FP_bytewise: found MSW[%u] = %20.15f, LSW[%u] = %20.15f\n",j,atmp,0,a[0]);
			if(atmp < 0.0) {
				msw_lt0 = 1;	/* MS word was < 0 prior to normalization */
				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
					cy = -1;
				else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
					cy = +1;
			}
			break;
		}
	}
	// v18: In Fermat-mod case, if above loop over the odd-indexed elements (the upper residue half
	// in the right-angle transform data layout) found no nonzero element, try the even-indexed ones:
	if(atmp == 0.0 && TRANSFORM_TYPE == RIGHT_ANGLE) {
		for(j = n-2; j >= 0; j -= TRANSFORM_TYPE)
		{
		#ifdef USE_AVX512
			j1 = (j & mask03) + br16[j&15];
		#elif defined(USE_AVX)
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
			atmp = a[j1];
			if(atmp != 0.0) {
			//	fprintf(stderr,"convert_res_FP_bytewise: found MSW[%u] = %20.15f, LSW[%u] = %20.15f\n",j,atmp,0,a[0]);
				if(atmp < 0.0) {
					msw_lt0 = 1;	/* MS word was < 0 prior to normalization */
					if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
						cy = -1;
					else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
						cy = +1;
				}
				break;
			}
		}
	}

	/*...Now form the SH residues, converting to positive-digit form along the way...	*/

	curr_char = 0;	/* Current byte to be written to in the ui64_arr_out array */
	nbits = 0;		/* Total bits accumulated so far in the residue	*/
	rbits = 0;		/* # of Upper bits left over from processing of previous word	*/
	curr_wd64 = 0;	/*      Upper bits left over from processing of previous word	*/

	if(TRANSFORM_TYPE == REAL_WRAPPER)
	{
		bimodn = 0;
		ii = 0;		/* Index into the BASE array. */
		/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
		if(bw > 0)
			ii = 1;

		for(j = 0; j < n; j++)
		{
		#ifdef USE_AVX512
			j1 = (j & mask03) + br16[j&15];
		#elif defined(USE_AVX)
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
			atmp = a[j1];
			if(atmp != NINT(atmp)) {
				sprintf(cbuf,"%s: Input float-residue elements must have 0 fractional part! A[%u (of %u)] = %20.10f",func,j,n,atmp);
				ASSERT(0, cbuf);
			}
			itmp = (int64)(atmp+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0) {			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
				itmp += (base[ii]);
				cy = -1;
			} else {
				cy = 0;
			}
			ASSERT(itmp >= 0,"convert_res_FP_bytewise: itmp >= 0");

		/* Update 8-byte residue buffer last, since this one modifies itmp: */
			ASSERT(rbits < 8,"convert_res_FP_bytewise: rbits < 8");
			ASSERT(curr_wd64 < (1<<rbits),"convert_res_FP_bytewise: curr_wd64 >= 2^rbits!");

			itmp = (itmp << rbits) + curr_wd64;
			curr_bits = bits[ii] + rbits;

			rbits = curr_bits;
			for(k = 0; k < curr_bits/8; k++) {
				ui64_arr_out[curr_char++] = itmp & 255;
				itmp = (uint64)itmp>>8;
				rbits -= 8;
			}
			curr_wd64 = (int)itmp;

			nbits += bits[ii];

			bimodn += bw;
			if(bimodn >= n) bimodn -= n;
			ii = (uint32)(sw - bimodn) >> 31;
		}
	}
	else	/* Complex right-angle transform */
	{
	  /* For right-angle transform, need to process odd and even-index elements
	  separately, via a pair of stride-2 passes through the array:
	  */
	  for(pass = 0; pass <=1; pass++)
	  {
		bimodn = n;

		for(j = pass; j < n; j += 2)
		{
			ii = (bimodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */
			bimodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */
			bimodn += ( ((int)ii-1) & n);		/*       add 0 if a bigword,   N if a smallword */
		#ifdef USE_AVX512
			j1 = (j & mask03) + br16[j&15];
		#elif defined(USE_AVX)
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
			atmp = a[j1];
			if(atmp != NINT(atmp)) {
				sprintf(cbuf,"%s: Input float-residue elements must have 0 fractional part! A[%u (of %u) = %20.10f] = ",func,j,n,atmp);
				ASSERT(0, cbuf);
			}
			itmp = (int64)(atmp+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0) {			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
				itmp += (base[ii]);
				cy = -1;
			} else {
				cy = 0;
			}
			ASSERT(itmp >= 0,"convert_res_FP_bytewise: itmp >= 0");

		/* Update 8-byte residue buffer last, since this one modifies itmp: */
			ASSERT(rbits < 8,"convert_res_FP_bytewise: rbits < 8");
			ASSERT(curr_wd64 < (1<<rbits),"convert_res_FP_bytewise: curr_wd64 >= 2^rbits!");

			itmp = (itmp << rbits) + curr_wd64;
			curr_bits = bits[ii] + rbits;

			rbits = curr_bits;
			for(k = 0; k < curr_bits/8; k++) {
				ui64_arr_out[curr_char++] = itmp & 255;
				itmp = (uint64)itmp>>8;
				rbits -= 8;
			}
			curr_wd64 = (int)itmp;
			nbits += bits[ii];
		}
	  }
	}
	/* Should have carryout of -1 Iff MS word < 0; otherwise expect 0 carry: */
	if(cy && (!msw_lt0 || cy != -1))
	{
		sprintf(cbuf, "convert_res_FP_bytewise: Illegal combination of nonzero carry = %d, msw_lt0 = %d\n", cy, msw_lt0);
		ASSERT(0, cbuf);
	}
	/* Residue should contain ceiling(p/8) bytes: */
	ASSERT(rbits < 8, "rbits >= 8");
	if(rbits) {
		ASSERT(curr_wd64 < (1<<rbits),"convert_res_FP_bytewise: curr_wd64 >= 2^rbits!");
		ui64_arr_out[curr_char++] = curr_wd64 & 255;
		curr_wd64 >>= 8;
	}
	ASSERT(curr_char == (p+7)/8,"convert_res_FP_bytewise: curr_char == (p+7)/8");
	ASSERT(nbits == p          ,"convert_res_FP_bytewise: nbits == p          ");
	ASSERT(curr_wd64 == 0      ,"convert_res_FP_bytewise: curr_wd64 == 0      ");

	// Remove the circular shift ... have no mi64_shrc function, so use that b-bit rightward cshift equivalent to (p-b)-bit left-cshift.
	// (But must guard against RES_SHIFT = 0, since in that case the left-shift count == p and mi64_shlc requires shift count strictly < p):
	/*** Jun 2022: For Fermat-mod case With nonzero shift, clause below kicks in, leading to j 1-too-large in mi64_shlc call if we use
	this formal count of #limbs-in-modulus. Since high residue bit only = 1 in the case of the final Pepin-test residue of a Fermat prime,
	must omit said high limb in residue-shift-and-sign-flip below, hence no p+63+(MODULUS_TYPE == MODULUS_TYPE_FERMAT):
	***/
	j = (p+63)>>6;	// # of 64-bit limbs
	if(RES_SHIFT) {
	//	fprintf(stderr,"convert_res_FP_bytewise: removing shift = %llu\n",RES_SHIFT);
		uint32 sign_flip = (MODULUS_TYPE == MODULUS_TYPE_FERMAT);
		mi64_shlc(u64_ptr, u64_ptr, p, p-RES_SHIFT,j,sign_flip);
		// If current residue R needed a sign-flip - again, this can only happen in the Fermat-mod case -
		// our shrc-done-as-shlc already took care of it. If not, need explicit negation. Rather than doing an
		// explicit Fm - R, can simply do a bitwise-complement of the residue vector and further += 2 of limb 0:
		if(sign_flip && !RES_SIGN) {	// sign_flip needed here since only do for Fermat case
		//	fprintf(stderr,"%s: Flipping sign of residue...\n",func);
			for(ii = 0; ii < j; ii++) { u64_ptr[ii] = ~u64_ptr[ii]; }	u64_ptr[0] += 2;
		}
	}
	/* Checksums: */
	if(Res64  ) *Res64 = ((uint64*)ui64_arr_out)[0];
	if(Res35m1) *Res35m1 = mi64_div_by_scalar64((uint64*)ui64_arr_out,two35m1,j,0x0);
	if(Res36m1) *Res36m1 = mi64_div_by_scalar64((uint64*)ui64_arr_out,two36m1,j,0x0);
//	fprintf(stderr,"Res35m1,Res36m1: %llu,%llu\n",*Res35m1,*Res36m1);
}

/*********************/
// Takes len-limb residue vector a[] and returns Selfridge-Hurwitz residues via the 3 arglist pointers:
void res_SH(uint64 a[], uint32 len, uint64*Res64, uint64*Res35m1, uint64*Res36m1)
{
	const uint64 two35m1 = (uint64)0x7FFFFFFFFull, two36m1 = (uint64)0xFFFFFFFFFull;	/* 2^35-1,36-1 */
	*Res64  = a[0];
	*Res35m1 = mi64_div_by_scalar64(a,two35m1,len,0x0);
	*Res36m1 = mi64_div_by_scalar64(a,two36m1,len,0x0);
}

/*********************/
uint32 get_default_factoring_depth(uint64 p)
{
/* Sample: here's how to set things to factor to a constant k-depth: */
#if 1
	const uint32 KMAX_BITS = 40;
	return (uint32) ceil(log(1.0*p)/log(2.0)) + 1 + KMAX_BITS;

#else

	uint32 qbitsmax;
/* These default depths are designed to match those of Prime95 v24, as described here:
	https://www.mersenneforum.org/showthread.php?t=4213
**** To-do: reverse order, add cases for > pmax, < pmin *****/
	     if(p <=23390000)
		qbitsmax = 66;
	else if(p > 23390000)	/* k ~= 40.5 bits */
		qbitsmax = 66;
	else if(p > 29690000)
		qbitsmax = 67;
	else if(p > 37800000)
		qbitsmax = 68;
	else if(p > 47450000)
		qbitsmax = 69;
	else if(p > 58520000)
		qbitsmax = 70;
	else if(p > 75670000)
		qbitsmax = 71;
	else if(p > 96830000)	/* k ~= 44.5 bits at the starting value of p*/
		qbitsmax = 72;

	/* If it's hardware on which factoring (at least my implementation) is slow, allow 64 bits max: */
  #if defined(INTEGER_MUL_32)
	if(qbitsmax > 64)
		qbitsmax = 64;
  #endif

	return qbitsmax;
#endif


}

/*********************/

int	is_hex_string(char*s, int len)
{
	int i;
	ASSERT(s != 0x0, "Null ptr to is_hex_string()");
	for(i = 0; i < len; ++i)
	{
		if( !isxdigit(s[i]) )
			return FALSE;
	}
	return TRUE;
}

/*********************/

/* In assignment lines in which the first four [decimal-numeric] values following the '=' are expected to contain ints
k,b,n,c of types [uint32,uint32,uint64, int32] which specify a modulus M = k*b^n+c, parse and sanity-check same.

Assumes:
	1. The arglist char* points to the first character following the test-type specifier (e.g. PRP, Test, Pminus1, etc)
		in the assignment line in question;
	2. Since my code at this time only supports Mersenne-number (M = 2^p - 1) and Fermat-number (M = 2^2^m + 1)
		moduli, assumes that k = 1, b = 2, c = (-1 or +1) and (n prime if c = -1, n  =  2^m if c = +1).

Returns two values:
	1. Return value of the function: For legal [k,b,n,c]-quartet specifiers, a char-pointer to the next character in
		the input string after the numerical value for the c-element of the quartet;
	2. The Mersenne-number exponent (p in 2^p - 1) or Fermat-number index (m in  2^2^m + 1).

Inputs in_str violating the above assumptions cause both return value to be set = 0.

Side effect:
	For legal assignment lines of the specified form, sets the value of the global MODULUS_TYPE,
	to MODULUS_TYPE_MERSENNE for c = -1 and MODULUS_TYPE_FERMAT for c = +1;
*/
char*check_kbnc(char*in_str, uint64*p) {
	const char func[] = "check_kbnc";
	char*char_addr = in_str, *cptr = 0x0;
	int i = 0;
	while(1) {
		if((char_addr = strstr(char_addr, "=")) == 0x0) {
			fprintf(stderr,"Expected '=' not found in assignment-specifying line!"); break;
		}
		char_addr++;
		while(isspace(*char_addr)) { ++char_addr; }	// Skip any whitespace following the equals sign
		if(is_hex_string(char_addr, 32)) {
			cptr = char_addr + 32;
			if((char_addr = strstr(cptr, ",")) == 0x0) {
				fprintf(stderr,"%s: Expected ',' not found in assignment-specifying line!\n",func); break;
			} else
				++char_addr;
		} else if(STREQN_NOCASE(char_addr,"n/a",3)) {
			cptr = char_addr + 3;
			if((char_addr = strstr(cptr, ",")) == 0x0) {
				fprintf(stderr,"%s: Expected ',' not found in assignment-specifying line!\n",func); break;
			} else
				++char_addr;
		}
		/*
		Since my code only supports Mersenne and Fermat-number moduli, check to ensure
		that k = 1, b = 2, c = (+1 or -1) and (n prime if c = -1, n  =  2^m if c = +1):
		*/
		i = (int)strtol(char_addr, &cptr, 10);
		if(i != 1) {
			fprintf(stderr,"%s: In modulus expression m = k*b^n+c, only k = 1 currently supported!\n",func); break;
		}
		if((char_addr = strstr(cptr, ",")) == 0x0) {
			fprintf(stderr,"%s: Expected ',' not found in assignment-specifying line!\n",func); break;
		}
		i = (int)strtol(char_addr+1, &cptr, 10);
		if(i != 2) {
			fprintf(stderr,"%s: In modulus expression m = k*b^n+c, only b = 2 currently supported!\n",func); break;
		}
		if((char_addr = strstr(cptr, ",")) == 0x0) {
			fprintf(stderr,"%s: Expected ',' not found in assignment-specifying line!\n",func); break;
		}
		*p = strtoull(char_addr+1, &cptr, 10);	ASSERT(*p != -1ull, "strtoull() overflow detected.");
		if(*p > PMAX) {
			fprintf(stderr,"%s: Exponent n in modulus expression m = k*b^n+c exceeds limit! (Suggest checking for unsigned overflow.)\n",func); break;
		}
		if((char_addr = strstr(cptr, ",")) == 0x0) {
			fprintf(stderr,"%s: Expected ',' not found in assignment-specifying line!\n",func); break;
		}
		i = (int)strtol(char_addr+1, &cptr, 10);
		if(ABS(i) != 1) {
			fprintf(stderr,"%s: In modulus expression m = k*b^n+c, only c = +-1 currently supported!\n",func); break;
		}
		if(i == -1) {
			uint32 phi32 = (*p >> 32);
			if((phi32 && !isPRP64(*p)) || (!phi32 && !is_prime((uint32)*p))) {
				fprintf(stderr,"%s: Mersenne exponent must be prime!\n",func); break;
			}
			MODULUS_TYPE = MODULUS_TYPE_MERSENNE;
		} else if(i == 1) {
			if(!isPow2_64(*p)) {
				fprintf(stderr,"%s: Fermat exponent must be a power of 2!\n",func); break;
			}
			*p = trailz64(*p);	// Code requires p contain Fermat-number *index* (not exponent), just as in a FermatTest assignment
			MODULUS_TYPE = MODULUS_TYPE_FERMAT;
		}
		return cptr;
	}
	// If we get here, return null *p and function values, indicating 'invalid assignment format':
	*p = 0ull;
	return 0x0;
}

/*********************/

/* v19: Add basic JSON-formatted result report. Examples:
	{"status":"C", "exponent":86749043, "worktype":"LL", "res64":"9EEA7CAD97A07648", "fft-length":4718592, "shift-count":6030412, "error-code":"00000000", "program":{"name":"Mlucas", "version":"18.0"}, "timestamp":"2019-11-11 01:23:45", "user":"madpoo", "computer":"mediaboy", "aid":"ABCDEF012456789 ABCDEF012456789"}
	{"status":"C", "exponent":110527, "worktype":"PRP-3", "res64":"E95075F756DD7BEB", "residue-type":1, "res2048":"E5E2DB84978E0355041AE377E588931B54FC75DCAD705044F21F17D0C8D5F524E98C535101C6DA9799F1433934FBAC2090761B1F4D8EA1F91AD63D03D477312E42F1CE7666C5E776A49A5BBDA146543C3CB1D74E0400CF6E81DF35173741289C76E69DB909726E50ECEE697F69A92E6BDF27A6AC6C9591EF97753F4555BABBFB26A385F78497ACAA4F7738A1E3C01564975DBDD3306C89FE7946B1523698BA334FAF2F53D74060BFDDDAF9D643E116D50DA7FDF2EB96CBC5D074602FDEDBD88E2706E0DED9324CEC6AD702016547E300748D6E5685C123CBB93744176B340B7DA6C1E478A685C774D554AAE5335C1FEBADD07A382A33BE1CE95075F756DD7BEB", "fft-length":6144, "shift-count":97673, "error-code":"00000000", "security-code":"F580BEFA", "program":{"name":"Prime95", "version":"29.8", "build":4, "port":10}, "timestamp":"2019-10-24 19:44:06", "errors":{"gerbicz":0}, "user":"gw_2", "computer":"Macbook_Pro"}
*/
/* v20: Added support for p-1 factoring results. Examples (these were actually found using gpuowl on GPU):
Factor found in stage 1:
	{"exponent":"102973951", "worktype":"PM1", "status":"F", "program":{"name":"Mlucas", "version":"20.0"}, "timestamp":"2020-02-25 03:35:38 UTC", "computer":"gfx906+sram-ecc-0", "aid":"341B9C13D982425C9C24DD5AA3AF96E7", "fft-length":5767168, "B1":1000000, "factors":["470377562071431809697977"]}
Factor found in stage 2:
	{"status":"F", "exponent":"107373143", "worktype":"PM1", "B1":"5500000", "fft-length":"6291456", "factors":["262356824147950958931679"], "program":{"name":"gpuowl", "version":"v6.11-520-g28dbf88"}, "computer":"gfx906+sram-ecc-0", "aid":"0A7F328173A3434EF690D1E7FF785484", "timestamp":"2021-02-25 14:12:05 UTC"}
No factor found:
	{"exponent":"103984877", "worktype":"PM1", "status":"NF", "program":{"name":"gpuowl", "version":"v6.11-142-gf54af2e"}, "timestamp":"2020-02-04 22:44:16 UTC", "user":"ewmayer", "computer":"gfx906+sram-ecc-0", "aid":"2927B92CB20B120F1353AC2F6FFB1C88", "fft-length":5767168, "B1":1000000, "B2":30000000}
2. B2 <= B1 means no stage 2 was run.
*/
/* v21:
Need to differentiate between PRP-CF results and the preceding PRP; set the otherwise-unused s2_partial flag to indicate PRP-CF.
*/
void generate_JSON_report(
	const uint32 isprime, const uint64 p, const uint32 n, const uint64 Res64, const char* Res2048, const char*timebuffer,
	const uint32 B1, const uint64 B2, const char*factor, const uint32 s2_partial,	// Quartet of p-1 fields
	char*cstr)	// cstr, takes the formatted output line; the preceding const-ones are inputs for that:
{
	int i,j,k;
	char ttype[11] = "\0", aid[33] = "\0";	// [test-type needs 11th | aid needs 33rd] char for \0
	const char prp_status[2] = {'C','P'};
	const char*pm1_status[2] = {"NF","F"};
	const char*false_or_true[2] = {"false","true"};
	// Attempt to read 32-hex-char Primenet assignment ID for current assignment (first line of WORKFILE):
	ASSERT((fp = mlucas_fopen(WORKFILE, "r")) != 0x0,"Workfile not found!");
	// v20.1.1: Parse first line whose leading non-WS char is alphabetic:
	char_addr = 0x0;
	while(fgets(in_line, STR_MAX_LEN, fp) != 0x0) {
		char_addr = in_line; while(isspace(*char_addr)) { ++char_addr; }
		if(isalpha(*char_addr)) break;
	}
	fclose(fp); fp = 0x0;
	ASSERT(strlen(char_addr) != 0 && isalpha(*char_addr),"Eligible assignment (leading non-WS char alphabetic) not found in workfile!");
	if(!strstr(in_line, ESTRING) && !(MODULUS_TYPE == MODULUS_TYPE_FERMAT && strstr(in_line, BIN_EXP)) ) {
		snprintf(cbuf,STR_MAX_LEN*2, "ERROR: Current exponent %s not found in %s file!\n",ESTRING,WORKFILE);
		ASSERT(0,cbuf);
	}
	// Is there a Primenet-server 32-hexit assignment ID in the assignment line? If so, include it in the JSON output:
	char_addr = strstr(in_line, "=");
	if(char_addr) {
		char_addr++;
		while(isspace(*char_addr)) { ++char_addr; }	// Skip any whitespace following the equals sign
		if(is_hex_string(char_addr, 32) && STRNEQN(char_addr,"00000000000000000000000000000000",32))
			strncpy(aid,char_addr,32);
	}
	// Write the result line. The 2 nested conditionals here are LL-or-PRP and has-AID-or-not:
	if(TEST_TYPE == TEST_TYPE_PRIMALITY) {
		snprintf(ttype,10,"LL");
		if(*aid) {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"worktype\":\"%s\", \"res64\":\"%016llX\", \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\", \"aid\":\"%s\"}\n",prp_status[isprime],p,ttype,Res64,n,RES_SHIFT,VERSION,timebuffer,aid);
		} else {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"worktype\":\"%s\", \"res64\":\"%016llX\", \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\"}\n",prp_status[isprime],p,ttype,Res64,n,RES_SHIFT,VERSION,timebuffer);
		}
	} else if(TEST_TYPE == TEST_TYPE_PRP && KNOWN_FACTORS[0]) {	// PRP-CF result
		// Print list of known factors used for CF test. Unlike the Primenet assignment formtting on the input side,
		// where all known factors are wrapped in a single bookending "" pair, the JSON needs ["factor1","factor2",..."].
		// We print that into cbuf, then include in full JSON result in cstr below:
		strcpy( cbuf, "[");	// Use cbuf as accumulator for loop below
		for(i = 0; KNOWN_FACTORS[i] != 0ull; i += 4) {
			k = mi64_getlen(KNOWN_FACTORS+i,4);	// k = number of nonzero limbs in curr_fac (alloc 4 limbs per in KNOWN_FACTORS[])
			// j = index of leading nonzero digit of decimal-string-printed factor:
			strcat( cbuf, "\"" );
			j = convert_mi64_base10_char(cstr, KNOWN_FACTORS+i, k, 0); strcat( cbuf, cstr+j );
			strcat( cbuf, "\"" );
			if(i < 36 && KNOWN_FACTORS[i+4])	// KNOWN_FACTORS alloc'ed for at most 10 factors of at most 4 limbs each
				strcat( cbuf, "," );
		}
		strcat( cbuf, "]");
		snprintf(ttype,10,"PRP-%u",PRP_BASE);
		if(*aid) {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"known-factors\":%s, \"worktype\":\"%s\", \"res64\":\"%016llX\", \"residue-type\":5, \"res2048\":\"%s\", \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\", \"aid\":\"%s\"}\n",prp_status[isprime],p,cbuf,ttype,Res64,Res2048,n,RES_SHIFT,VERSION,timebuffer,aid);
		} else {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"known-factors\":%s, \"worktype\":\"%s\", \"res64\":\"%016llX\", \"residue-type\":5, \"res2048\":\"%s\", \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\"}\n",prp_status[isprime],p,cbuf,ttype,Res64,Res2048,n,RES_SHIFT,VERSION,timebuffer);
		}
	} else if(TEST_TYPE == TEST_TYPE_PRP) {	// Only support type-1 PRP tests, so hardcode that subfield:
		snprintf(ttype,10,"PRP-%u",PRP_BASE);
		if(*aid) {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"worktype\":\"%s\", \"res64\":\"%016llX\", \"residue-type\":1, \"res2048\":\"%s\", \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\", \"aid\":\"%s\"}\n",prp_status[isprime],p,ttype,Res64,Res2048,n,RES_SHIFT,VERSION,timebuffer,aid);
		} else {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"worktype\":\"%s\", \"res64\":\"%016llX\", \"residue-type\":1, \"res2048\":\"%s\", \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\"}\n",prp_status[isprime],p,ttype,Res64,Res2048,n,RES_SHIFT,VERSION,timebuffer);
		}
	} else if(TEST_TYPE == TEST_TYPE_PM1) {	// For p-1 assume there was an AID in the assignment, even if an all-0s one:
		snprintf(ttype,10,"PM1");
		if(!strlen(factor)) {	// No factor was found:
		  if(*aid) {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%s\", \"exponent\":%llu, \"worktype\":\"%s\", \"fft-length\":%u, \"B1\":%u, \"B2\":%llu, \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\", \"aid\":\"%s\"}\n",pm1_status[0],p,ttype,n,B1,B2,VERSION,timebuffer,aid);
		  } else {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%s\", \"exponent\":%llu, \"worktype\":\"%s\", \"fft-length\":%u, \"B1\":%u, \"B2\":%llu, \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\"}\n",pm1_status[0],p,ttype,n,B1,B2,VERSION,timebuffer);
		  }
		} else {	// The factor in the eponymous arglist field was found:
		  if(B2 <= B1) {	// No stage 2 was run
		   if(*aid) {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%s\", \"exponent\":%llu, \"worktype\":\"%s\", \"fft-length\":%u, \"B1\":%u, \"factors\":[\"%s\"], \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\", \"aid\":\"%s\"}\n",pm1_status[1],p,ttype,n,B1,factor,VERSION,timebuffer,aid);
		   } else {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%s\", \"exponent\":%llu, \"worktype\":\"%s\", \"fft-length\":%u, \"B1\":%u, \"factors\":[\"%s\"], \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\"}\n",pm1_status[1],p,ttype,n,B1,factor,VERSION,timebuffer);
		   }
		  } else {	// Include B2 and flag indicating whether the s2 interval was completely covered or not. Factor must be in "" due to possibility of > 64-bit, which overflows a JSON int:
		   if(*aid) {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%s\", \"exponent\":%llu, \"worktype\":\"%s\", \"fft-length\":%u, \"B1\":%u, \"B2\":%llu, \"partial-stage-2\":%s, \"factors\":[\"%s\"], \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\", \"aid\":\"%s\"}\n",pm1_status[1],p,ttype,n,B1,B2,false_or_true[s2_partial],factor,VERSION,timebuffer,aid);
		   } else {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%s\", \"exponent\":%llu, \"worktype\":\"%s\", \"fft-length\":%u, \"B1\":%u, \"B2\":%llu, \"partial-stage-2\":%s, \"factors\":[\"%s\"], \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\"}\n",pm1_status[1],p,ttype,n,B1,B2,false_or_true[s2_partial],factor,VERSION,timebuffer);
		   }
		  }
	}
	} else
		ASSERT(0, "Unsupported test type!");
}

/*********************/

// Sets function pointers for DIF|DIT pass1 based on value of radix0:
void dif1_dit1_func_name(
	const int radix0,
	void (**func_dif_pass1)(double [], int),
	void (**func_dit_pass1)(double [], int)
) {
  switch(radix0) {
	case  5:	*func_dif_pass1 =  radix5_dif_pass1;	*func_dit_pass1 =  radix5_dit_pass1;	break;
	case  6:	*func_dif_pass1 =  radix6_dif_pass1;	*func_dit_pass1 =  radix6_dit_pass1;	break;
	case  7:	*func_dif_pass1 =  radix7_dif_pass1;	*func_dit_pass1 =  radix7_dit_pass1;	break;
	case  8:	*func_dif_pass1 =  radix8_dif_pass1;	*func_dit_pass1 =  radix8_dit_pass1;	break;
	case  9:	*func_dif_pass1 =  radix9_dif_pass1;	*func_dit_pass1 =  radix9_dit_pass1;	break;
	case 10:	*func_dif_pass1 = radix10_dif_pass1;	*func_dit_pass1 = radix10_dit_pass1;	break;
	case 11:	*func_dif_pass1 = radix11_dif_pass1;	*func_dit_pass1 = radix11_dit_pass1;	break;
	case 12:	*func_dif_pass1 = radix12_dif_pass1;	*func_dit_pass1 = radix12_dit_pass1;	break;
	case 13:	*func_dif_pass1 = radix13_dif_pass1;	*func_dit_pass1 = radix13_dit_pass1;	break;
	case 14:	*func_dif_pass1 = radix14_dif_pass1;	*func_dit_pass1 = radix14_dit_pass1;	break;
	case 15:	*func_dif_pass1 = radix15_dif_pass1;	*func_dit_pass1 = radix15_dit_pass1;	break;
	case 16:	*func_dif_pass1 = radix16_dif_pass1;	*func_dit_pass1 = radix16_dit_pass1;	break;
	case 18:	*func_dif_pass1 = radix18_dif_pass1;	*func_dit_pass1 = radix18_dit_pass1;	break;
	case 20:	*func_dif_pass1 = radix20_dif_pass1;	*func_dit_pass1 = radix20_dit_pass1;	break;
	case 22:	*func_dif_pass1 = radix22_dif_pass1;	*func_dit_pass1 = radix22_dit_pass1;	break;
	case 24:	*func_dif_pass1 = radix24_dif_pass1;	*func_dit_pass1 = radix24_dit_pass1;	break;
	case 26:	*func_dif_pass1 = radix26_dif_pass1;	*func_dit_pass1 = radix26_dit_pass1;	break;
	case 28:	*func_dif_pass1 = radix28_dif_pass1;	*func_dit_pass1 = radix28_dit_pass1;	break;
	case 30:	*func_dif_pass1 = radix30_dif_pass1;	*func_dit_pass1 = radix30_dit_pass1;	break;
	case 32:	*func_dif_pass1 = radix32_dif_pass1;	*func_dit_pass1 = radix32_dit_pass1;	break;
	case 36:	*func_dif_pass1 = radix36_dif_pass1;	*func_dit_pass1 = radix36_dit_pass1;	break;
	case 40:	*func_dif_pass1 = radix40_dif_pass1;	*func_dit_pass1 = radix40_dit_pass1;	break;
	case 44:	*func_dif_pass1 = radix44_dif_pass1;	*func_dit_pass1 = radix44_dit_pass1;	break;
	case 48:	*func_dif_pass1 = radix48_dif_pass1;	*func_dit_pass1 = radix48_dit_pass1;	break;
	case 52:	*func_dif_pass1 = radix52_dif_pass1;	*func_dit_pass1 = radix52_dit_pass1;	break;
	case 56:	*func_dif_pass1 = radix56_dif_pass1;	*func_dit_pass1 = radix56_dit_pass1;	break;
	case 60:	*func_dif_pass1 = radix60_dif_pass1;	*func_dit_pass1 = radix60_dit_pass1;	break;
	case 63:	*func_dif_pass1 = radix63_dif_pass1;	*func_dit_pass1 = radix63_dit_pass1;	break;
	case 64:	*func_dif_pass1 = radix64_dif_pass1;	*func_dit_pass1 = radix64_dit_pass1;	break;
//	case 112:	*func_dif_pass1 = radix112_dif_pass1;	*func_dit_pass1 = radix112_dit_pass1;	break;
//	case 120:	*func_dif_pass1 = radix120_dif_pass1;	*func_dit_pass1 = radix120_dit_pass1;	break;
	case 128:	*func_dif_pass1 = radix128_dif_pass1;	*func_dit_pass1 = radix128_dit_pass1;	break;
	case 144:	*func_dif_pass1 = radix144_dif_pass1;	*func_dit_pass1 = radix144_dit_pass1;	break;
	case 160:	*func_dif_pass1 = radix160_dif_pass1;	*func_dit_pass1 = radix160_dit_pass1;	break;
	case 176:	*func_dif_pass1 = radix176_dif_pass1;	*func_dit_pass1 = radix176_dit_pass1;	break;
	case 192:	*func_dif_pass1 = radix192_dif_pass1;	*func_dit_pass1 = radix192_dit_pass1;	break;
	case 208:	*func_dif_pass1 = radix208_dif_pass1;	*func_dit_pass1 = radix208_dit_pass1;	break;
	case 224:	*func_dif_pass1 = radix224_dif_pass1;	*func_dit_pass1 = radix224_dit_pass1;	break;
	case 240:	*func_dif_pass1 = radix240_dif_pass1;	*func_dit_pass1 = radix240_dit_pass1;	break;
	case 256:	*func_dif_pass1 = radix256_dif_pass1;	*func_dit_pass1 = radix256_dit_pass1;	break;
	case 288:	*func_dif_pass1 = radix288_dif_pass1;	*func_dit_pass1 = radix288_dit_pass1;	break;
	case 320:	*func_dif_pass1 = radix320_dif_pass1;	*func_dit_pass1 = radix320_dit_pass1;	break;
	case 352:	*func_dif_pass1 = radix352_dif_pass1;	*func_dit_pass1 = radix352_dit_pass1;	break;
	case 384:	*func_dif_pass1 = radix384_dif_pass1;	*func_dit_pass1 = radix384_dit_pass1;	break;
	case 512:	*func_dif_pass1 = radix512_dif_pass1;	*func_dit_pass1 = radix512_dit_pass1;	break;
	case 768:	*func_dif_pass1 = radix768_dif_pass1;	*func_dit_pass1 = radix768_dit_pass1;	break;
	case 960:	*func_dif_pass1 = radix960_dif_pass1;	*func_dit_pass1 = radix960_dit_pass1;	break;
	case 992:	*func_dif_pass1 = radix992_dif_pass1;	*func_dit_pass1 = radix992_dit_pass1;	break;
	case 1008:	*func_dif_pass1 = radix1008_dif_pass1;	*func_dit_pass1 = radix1008_dit_pass1;	break;
	case 1024:	*func_dif_pass1 = radix1024_dif_pass1;	*func_dit_pass1 = radix1024_dit_pass1;	break;
	case 4032:	*func_dif_pass1 = radix4032_dif_pass1;	*func_dit_pass1 = radix4032_dit_pass1;	break;
//	case 4096:	*func_dif_pass1 = radix4096_dif_pass1;	*func_dit_pass1 = radix4096_dit_pass1;	break;
	default:
		sprintf(cbuf,"ERROR: radix %d not available for [dif,dit] pass1. Halting...\n",radix0); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
  }
}

/*********************/

/* Read known factors starting from a char-ptr corr. to a line of worktodo file or substring of such, in
double-quote-and comma-delimited "factor1,factor2,...,factorN" format. Each alleged known-factor is checked for
correctness (in both the is-base-2-Fermat-PRP and divides-the-modulus senses), with factors limited to < 2^256;
the factors are stored in the global KNOWN_FACTORS[], with a limit of 10 known factors per household.
Assumes: Modulus type and binary exponent have been set prior to function call;
Returns: The number of factors read in.
*/
uint32 extract_known_factors(uint64 p, char*fac_start) {
	uint32 i, findex, fbits, lenf, nchar, nfac = 0;
	uint64 *fac = 0x0, twop[4], quo[4],rem[4];	// fac = ptr to each mi64-converted factor input string;
	uint256 p256,q256,res256;
	char*cptr = fac_start+1;
	ASSERT(fac_start[0] == '\"',"Known-factors line of worktodo must consist of a comma-separated list of such enclosed in double-quotes!");
	/* If it's a Fermat number, need to check size of 2^ESTRING: */
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
		findex = (uint32)p;
		if(findex <= MAX_PRIMALITY_TEST_BITS)
			p = (uint64)1 << findex;
		else
			ASSERT(0,"nbits_in_p <= MAX_PRIMALITY_TEST_BITS");
	}
	// Factors separated by commas (first clause of while()); list terminated with " (2nd clause):
	while((char_addr = strstr(cptr,",")) != 0x0 || (char_addr = strstr(cptr,"\"")) != 0x0) {
		nchar = char_addr - cptr;
		strncpy(cbuf,cptr,nchar);	cbuf[nchar] = '\0';	// Extract current-factor-as-string into cbuf
		// Convert stringified factor f to mi64 form:
		lenf = 0; fac = convert_base10_char_mi64(cbuf, &lenf);	// This does the mem-alloc for us
		ASSERT(lenf > 0, "Error converting known-factor string!");
		ASSERT(lenf < 5, "known-factor out of range, must be < 2^256!");
		fbits = (lenf<<6) - mi64_leadz(fac, lenf);
		// Make sure the alleged factor is of the proper form:
		// For Mersenne M(p), q = 2.k.p + 1, with p prime; For Fermat F_n = 2^2^n+1, q = k.2^(n+2) + 1
		// and we store the binary exponent 2^n in p, and 2^(n+2) in twop (yes, a misnomer in this case):
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			mi64_set_eq_scalar(twop,p<<1,lenf);	ASSERT(p < 0x8000000000000000ull, "Mersenne exponent limited to 63 bits!");
		} else {
			mi64_set_eq_scalar(twop,p<<2,lenf);	ASSERT(p < 0x4000000000000000ull, "Fermat-number index must be < 62!");
		}
		mi64_div(fac,twop, lenf,lenf, quo,rem);
		i = mi64_cmp_eq_scalar(rem,1ull,lenf);	ASSERT(i,"Factor not of required form!");
		// Alloc 4 limbs per factor in KNOWN_FACTORS; if current factor needs just 1 there's no uninited
		// problem with the high limbs since KNOWN_FACTORS is zeroed at start of each new assignment:
		ASSERT(nfac < 10, "Limit of 10 known factors!");
		mi64_set_eq(KNOWN_FACTORS + 4*nfac++,fac,lenf);
		// Verify that F is a base-3 Fermat-PRP via binary modpow, 3^(q-1) == 1 (mod q):
		ASSERT(mi64_pprimeF(fac,3ull,lenf),"Factor-is-base-3-PRP check fails!");
		// Verify that it's a factor via binary modpow:
		p256.d0 = p; p256.d1 = p256.d2 = p256.d3 = 0ull;
		q256.d0 = KNOWN_FACTORS[0];	q256.d1 = KNOWN_FACTORS[1];	q256.d2 = KNOWN_FACTORS[2];	q256.d3 = KNOWN_FACTORS[3];
		res256 = twopmmodq256(p256,q256);
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			ASSERT(CMPEQ256(res256,ONE256),"Factor-divides-modulus check fails!");
		} else {
			res256.d0 += 1ull;	// Fermat case: check that 2^p == -1 == q - 1 (mod q):
			ASSERT(CMPEQ256(res256,q256),"Factor-divides-modulus check fails!");
		}
		// If find any duplicate-entries in input list, warn & remove:
		if(nfac > 1) {
			for(i = 0; i < nfac-1; i++) {
				// Incremented nfac already, so just-added entry in slot (nfac-1):
				if(mi64_cmp_eq(KNOWN_FACTORS + 4*i, KNOWN_FACTORS + 4*(nfac-1), 4)) {
					mi64_clear(KNOWN_FACTORS + 4*(--nfac), 4);
					// Using cbuf as both string-arg and target string is problematic, so use 2nd string-global cstr as target:
					snprintf(cstr,STR_MAX_LEN, "WARNING: p = %llu, known-factor list entry %s is a duplicate ... removing.\n",p,cbuf);
					fprintf(stderr,"%s",cstr);
				}
			}
		}
		cptr = char_addr+1;	// Advance 1-char past the current , or "
	}
	if(char_addr != 0x0) {
		sprintf(cbuf,"%s: Unrecognized token sequence in parsing known-factors portion of assignment: \"%s\".",WORKFILE,fac_start);	ASSERT(0,cbuf);
	}
	ASSERT(nfac != 0,"Must specify at least one known factor!");
// A bit of just-for-fun code: For smaller moduli N, use mi64 utils to see if cofactor C is a base-3 PRP:
#if 0
	const char mod_type[2] = {'-','+'}, *is_prp[] = {"is not","is"}, exclam[2] = {'.','!'};
	int j,k;
	uint64 *mvec = 0x0, *qvec = 0x0, itmp64;
	j = (p+63+(MODULUS_TYPE == MODULUS_TYPE_FERMAT))>>6;	// # of 64-bit limbs
	mvec = ALLOC_UINT64(mvec,j);	// Modulus N
	qvec = ALLOC_UINT64(qvec,j);	// Quotient stores cofactor C = N/F
	if(!mvec || !qvec) {
		sprintf(cbuf, "ERROR: unable to allocate arrays mvec,qvec in extract_known_factors.\n"); fprintf(stderr,"%s", cbuf);
		ASSERT(0,cbuf);
	}

	// Compute Modulus N ... note mi64-vecs have no cache-oriented element padding:
	mvec[j-1] = 0ull;
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		itmp64 = -1ull;
		// Loop rather than call to mi64_set_eq_scalar here, since need to set all elts = -1:
		for(i = 0; i < j; i++) { mvec[i] = itmp64; }
		mvec[j-1] >>= 64-(p&63);	// Leading word needs >> to leave just low p%64 bits set
	} else {
		// j = uint64 vector length; init sans the leading '1' word, then increment prior to mi64_div
		mi64_clear(mvec,j);
		mvec[j++] = mvec[0] = 1ull;	// Post-increment here
	}

	// Multiply each known-factor with current partial product of factors.
	// Use BASE_MULTIPLIER_BITS to store factor product here, but need curr_fac[] for intermediate partial products:
	BASE_MULTIPLIER_BITS[0] = 1ull;	lenf = 1;
	uint64 curr_fac[20];
	for(i = 0; KNOWN_FACTORS[i] != 0ull; i += 4) {
		k = mi64_getlen(KNOWN_FACTORS+i,4);	// k = number of nonzero limbs in curr_fac (alloc 4 limbs per in KNOWN_FACTORS[])
		// Multiply factor into current partial product of factors; use curr_fac[] array to store product to work around none-of-3-input-pointers-may-coincide restriction in mi64_mul_vector:
		mi64_mul_vector(BASE_MULTIPLIER_BITS,lenf, KNOWN_FACTORS+i,k, curr_fac,&lenf);
		mi64_set_eq(BASE_MULTIPLIER_BITS,curr_fac,lenf);
	}
	ASSERT(lenf <= 20, "Product of factors too large to fit into curr_fac[]!");

	// Since F << N, use Mont-mul-div for C - quotient overwrites N, no rem-vec needed, just verify that F is in fact a divisor:
	ASSERT(1 == mi64_div(mvec,BASE_MULTIPLIER_BITS, j,lenf, qvec,0x0), "C = N/F should have 0 remainder!");
	k = mi64_getlen(qvec,j);	// j = number of nonzero limbs in cofactor C
	i = mi64_pprimeF(qvec,3,k);
	printf("2^%llu %c 1 %s a base-3 Fermat-PRP%c\n",p,mod_type[MODULUS_TYPE == MODULUS_TYPE_FERMAT],is_prp[i],exclam[i]);
	free((void *)mvec);	mvec = 0x0;
	exit(0);
#endif
	return nfac;
}

/*********************/

/*
If p != 0 (this requires vec2 == 0x0):
	For (MODULUS_TYPE == MODULUS_TYPE_[FERMAT|MERSENNE]), take GCD of 2^p[+|-]1 and nlimb resarr[],
	with nlimb := ceiling(p/64). A nonzero return value indicates a nontrivial GCD, after dividing
	out any known factors stored in the extern KNOWN_FACTORS array.
If p == 0 (this requires vec2 != 0x0):
	Take GCD of the 2 nlimb-inputs vec1[] and vec2[].
	A nonzero return value indicates a nontrivial GCD.
The decimal value of the GCD is returned in gcd_str, presumed to be dimensioned >= 1024 chars:
*/
uint32 gcd(uint32 stage, uint64 p, uint64*vec1, uint64*vec2, uint32 nlimb, char*const gcd_str) {
#if !INCLUDE_GMP
	#warning INCLUDE_GMP defined == 0 at compile time ... No GCDs will be done on p-1 outputs.
	snprintf(cbuf,STR_MAX_LEN*2,"INCLUDE_GMP defined == 0 at compile time ... No GCD will be done.\n");
	mlucas_fprint(cbuf,1);
	return 0;	// If user turns off p-1 support, keep the decl of gcd() to allow pm1.c to build
#else
	// Unlike standard types and Mlucas internal structs, GMP objects must be declared before any expressions,
	// else GCC emits "error: a label can only be part of a statement and a declaration is not a statement":
	mpz_t gmp_arr1, gmp_arr2, gmp_one, gmp_d, gmp_r, gmp_q;
	mp_bitcnt_t gmp_exp;
	size_t gmp_size, sz1,sz2;
	uint32 i, retval = 0;
	double tdiff = 0.0, clock1, clock2;
	clock1 = getRealTime();
	ASSERT(vec1 != 0x0, "Null-pointer vec1 input to GCD()!");
	ASSERT(!(p && vec2), "One and only one of p and vec2 args to GCD() must be non-null!");
	mpz_init(gmp_arr1); mpz_init(gmp_arr2);
	// Init divisor, remainder, quotient, in case of nontrivial raw GCD and >= 1 known factors:
	mpz_init(gmp_d); mpz_init(gmp_r); mpz_init(gmp_q);
	mpz_init_set_ui(gmp_one,1ull); gmp_exp = p;
	// Import vec1 into GMP array1, least-sign. element first, host byte order within each word, at 64-bit width:
	mpz_import(gmp_arr1, nlimb, -1, sizeof(uint64), 0, 0, vec1);
	if(p != 0) {
		ASSERT(nlimb == (p + 63 + (MODULUS_TYPE == MODULUS_TYPE_FERMAT))>>6, "Bad inputs to GCD()!");
		mpz_mul_2exp(gmp_arr2, gmp_one,gmp_exp);
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)	// 2^p-1:
			mpz_sub(gmp_arr2, gmp_arr2,gmp_one);
		else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)// F(m): p holds 2^m, so F(m) = 2^p+1:
			mpz_add(gmp_arr2, gmp_arr2,gmp_one);
	} else {
		mpz_import(gmp_arr2, nlimb, -1, sizeof(uint64), 0, 0, vec2);
	}
	sz1 = mpz_sizeinbase(gmp_arr1,2);// gmp_printf("Input1 has %llu bits\n",sz1);
	sz2 = mpz_sizeinbase(gmp_arr2,2);// gmp_printf("Input2 has %llu bits\n",sz2);
	// Take gcd and return in gmp_arr1:
	mpz_gcd(gmp_arr1, gmp_arr1,gmp_arr2);
	gmp_size = mpz_sizeinbase(gmp_arr1,2);
	if(gmp_size < 2) {
		goto gcd_return;	// GCD = 0 or 1
	} else {
		if(KNOWN_FACTORS[0]) fprintf(stderr,"Raw GCD has %llu bits ... dividing out any known factors...\n",(uint64)gmp_size);
		for(i = 0; i < 40; i += 4) {	// Current limit = 10 factors, each stored in a 4-limb field, i.e. < 2^256
			if(!KNOWN_FACTORS[i])
				break;
			mpz_import(gmp_d, 4, -1, sizeof(uint64), 0, 0, KNOWN_FACTORS+i);
			mpz_tdiv_qr(gmp_q,gmp_r, gmp_arr1,gmp_d);
			if(!mpz_size(gmp_r))	// This known factor divides the GCD; replace the latter with the quotient:
				mpz_set(gmp_arr1, gmp_q);
		}
	}
	// Recompute bitlength of GCD
	gmp_size = mpz_sizeinbase(gmp_arr1,2);
	if(gmp_size < 2)
		goto gcd_return;	// GCD = 0 or 1
	// Now the base-10 digit count:
	gmp_size = mpz_sizeinbase(gmp_arr1,10);
	// Anything >= 900 digits (~90% the value of our STR_MAX_LEN dimensioning of I/O strings) treated as suspect:
	if(gmp_size >= 900) {
		snprintf(cbuf,STR_MAX_LEN*2, "GCD has %u digits -- possible data corruption, aborting.\n",(uint32)gmp_size);
		mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
	}
	retval = 1;
gcd_return:
	if(!p) {
		gmp_snprintf(gcd_str,STR_MAX_LEN,"%Zd",gmp_arr1);
		gmp_snprintf(cbuf,STR_MAX_LEN*2,"GCD(A[%llu bits], B[%llu bits]) = %s\n",sz1,sz2,gcd_str);
	} else if(retval) {
		gmp_snprintf(gcd_str,STR_MAX_LEN,"%Zd",gmp_arr1);
		gmp_snprintf(cbuf,STR_MAX_LEN*2,"Found %u-digit factor in Stage %u: %s\n",gmp_size,stage,gcd_str);
	} else {	// Caller can use either return value or empty gcd_str as proxy for "no factor found"
		gcd_str[0] = '\0';
		gmp_snprintf(cbuf,STR_MAX_LEN*2,"Stage %u: No factor found.\n",stage);
	}
	mlucas_fprint(cbuf,1);
	clock2 = getRealTime(); tdiff = clock2 - clock1;
	snprintf(cbuf,STR_MAX_LEN*2,"Time for GCD =%s\n",get_time_str(tdiff));
	mlucas_fprint(cbuf,1);
	// Done with the GMP arrays:
	mpz_clear(gmp_arr1); mpz_clear(gmp_arr2); mpz_clear(gmp_d); mpz_clear(gmp_r); mpz_clear(gmp_q);
	return retval;
#endif	// INCLUDE_GMP ?
}

/*********************/

/*
For (MODULUS_TYPE == MODULUS_TYPE_[FERMAT|MERSENNE]), computes inverse (mod 2^p[+|-]1)
of vec1[] and returns it in vec2[].
GMP mod-inverse; arglist as for mpz_gcd but also returns int:

		int mpz_invert (mpz t rop, const mpz t op1, const mpz t op2)

	Compute the inverse of op1 modulo op2 and put the result in rop. If the inverse exists, the
	return value is non-zero and rop will satisfy 0 ≤ rop < |op2| (with rop = 0 possible only when
	|op2| = 1, i.e., in the somewhat degenerate zero ring). If an inverse doesn’t exist the return
	value is zero and rop is undefined. The behaviour of this function is undefined when op2 = 0.
*/
void modinv(uint64 p, uint64*vec1, uint64*vec2, uint32 nlimb) {
#if !INCLUDE_GMP
	#warning INCLUDE_GMP defined == 0 at compile time ... No MODINV support.
	// If user turns off p-1 support, keep the decl of modinv() to allow pm1.c to build
	mi64_clear(vec2,nlimb);	// Return 0
#else
	// Unlike standard types and Mlucas internal structs, GMP objects must be declared before any expressions,
	// else GCC emits "error: a label can only be part of a statement and a declaration is not a statement":
	mpz_t gmp_arr1, gmp_arr2, gmp_one;
	mp_bitcnt_t gmp_exp;
	size_t gmp_size;
	uint32 i, retval = 0;
	size_t inv_limbs;
	uint64 *export_result_addr;
	double tdiff = 0.0, clock1, clock2;
	clock1 = getRealTime();
	ASSERT(vec1 != 0x0 && vec2 != 0x0, "Null-pointer input to MODINV()!");
	mpz_init(gmp_arr1); mpz_init(gmp_arr2);
	mpz_init_set_ui(gmp_one,1ull); gmp_exp = p;
	// Import vec1 into GMP array1, least-sign. element first, host byte order within each word, at 64-bit width:
	// void mpz_import (mpz_t rop, size_t count, int order, size_t size, int	[Function] endian, size_t nails, const void *op)
	mpz_import(gmp_arr1, nlimb, -1, sizeof(uint64), 0, 0, vec1);
	ASSERT((p != 0) && (nlimb == (p + 63 + (MODULUS_TYPE == MODULUS_TYPE_FERMAT))>>6), "Bad inputs to MODINV()!");
	mpz_mul_2exp(gmp_arr2, gmp_one,gmp_exp);
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)	// 2^p-1:
		mpz_sub(gmp_arr2, gmp_arr2,gmp_one);
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)// F(m): p holds 2^m, so F(m) = 2^p+1:
		mpz_add(gmp_arr2, gmp_arr2,gmp_one);
//	gmp_printf("Input1 has %llu bits\n",mpz_sizeinbase(gmp_arr1,2));
//	gmp_printf("Input2 has %llu bits\n",mpz_sizeinbase(gmp_arr2,2));
	/*
	GMP mod-inverse; arglist as for mpz_gcd but also returns int:

		int mpz_invert (mpz t rop, const mpz t op1, const mpz t op2)

	"Compute the inverse of op1 modulo op2 and put the result in rop. If the inverse exists, the
	return value is non-zero and rop will satisfy 0 ≤ rop < |op2| (with rop = 0 possible only when
	|op2| = 1, i.e., in the somewhat degenerate zero ring). If an inverse doesn’t exist the return
	value is zero and rop is undefined. The behaviour of this function is undefined when op2 = 0."
	*/
	// Return inverse in gmp_arr1:
	retval = mpz_invert(gmp_arr1, gmp_arr1,gmp_arr2);
	gmp_size = mpz_sizeinbase(gmp_arr1,2);
	if(!retval) {
		snprintf(cbuf,STR_MAX_LEN*2,"MODINV: Fatal error: inverse does not exist.\n");
		mlucas_fprint(cbuf,0); ASSERT(0,cbuf);
	}
	// Export the result from gmp_arr1 to destination array vec2:
	// void * mpz_export (void *rop, size_t *countp, int order, size_t size, int	[Function] endian, size_t nails, const mpz_t op)
	export_result_addr = mpz_export(vec2, &inv_limbs, -1, sizeof(uint64), 0, 0, gmp_arr1);
	ASSERT(inv_limbs <= nlimb && export_result_addr == vec2, "GMP was unable to export result to the specified target array!");
	// Explicitly zero any excess limbs left at top of vec2:
	for(i = inv_limbs; i < nlimb; i++) {
		vec2[i] = 0ull;
	}
	// Done with the GMP arrays:
	mpz_clear(gmp_arr1); mpz_clear(gmp_arr2);
#endif	// INCLUDE_GMP ?
}

/*********************/
// Tries to read restartfile with name [fname], allegedly containing restart data for Mersenne or Fermat
// modulus with binary exponent [expo]. Returns 1 on successful read and checksum validation, 0 otherwise.
// Assumes the globals TEST_TYPE and MODULUS_TYPE have been properly set in main.
// Requires pointers arr1 (and arr2 if a PRP-test residue) to scratch arrays with sufficient allocation
// to hold the bytewise residue(s) stored in the file to be passed in.
int restart_file_valid(const char*fname, const uint64 p, uint8*arr1, uint8*arr2)
{
	int retval = 0;
	uint32 j;
	uint64 Res64,Res35m1,Res36m1, i1,i2,i3, itmp64;
	FILE*fptr = mlucas_fopen(fname,"r");
	if(fptr) {
		retval = read_ppm1_savefiles(fname, p, &j, fptr, &itmp64,
										arr1, &Res64,&Res35m1,&Res36m1,	// Primality-test residue
										arr2, &i1   ,&i2     ,&i3     );// [PRP-test] G-check residues
		fclose(fptr);
	}
	if(!retval && strstr(cbuf, "read_ppm1_savefiles"))
		mlucas_fprint(cbuf,1);
	return retval;
}

/*********************/
// Returns line number of substring find_string found in file named fname. Assumes file exists in run directory - if not, hard-ASSERT.
// Depending on whether the associated integer arg find_before_line_number = 0 or not, the mandatory cstr arg will contain a copy of
// either the first or last-line-before-line-number matching find_str, respectively. If no matches, 0 is returned and cstr = '\0'.
// To find the line number and content of the last occurrence of find_str, period, set find_before_line_number = -1:
uint32 filegrep(const char*fname, const char*find_str, char*cstr, uint32 find_before_line_number)
{
	uint32 curr_line = 0, found_line = 0;
	ASSERT(cstr != 0x0, "filegrep(): cstr pointer argument must be non-null!");
	cstr[0] = '\0';
	if(strlen(find_str) == 0)	// Nothing to find
		return 0;
	FILE*fptr = mlucas_fopen(fname,"r");
	if(fptr){
		while(fgets(in_line, STR_MAX_LEN, fptr)) {
			++curr_line;
			if(find_before_line_number && curr_line >= find_before_line_number)
				break;
			if(strstr(in_line,find_str) != 0x0) {
				found_line = curr_line;
				strcpy(cstr,in_line);
				if(!find_before_line_number)	// if find_before_line_number == 0, return first occurrence:
					break;
			}
		}
		fclose(fptr);
	} else {
		sprintf(cbuf,"filegrep error: file %s not found.\n",fname);
		ASSERT(0, cbuf);
	}
	if(strlen(cstr) != 0)
		return found_line;
	else return 0;
}

/*********************/

void write_fft_debug_data(double a[], int jlo, int jhi)
{
	int j,j1;
	const char dbg_fname[] = "FFT_DEBUG.txt";
	ASSERT(dbg_file == 0x0, "dbg_file != 0x0 prior to mlucas_fopen");
	dbg_file = mlucas_fopen(dbg_fname, "a");
	ASSERT(dbg_file != 0x0, "Unable to open dbg_file!");
	fprintf(dbg_file, "RE_IM_STRIDE = %d\n", RE_IM_STRIDE);
	fprintf(dbg_file, "%s\n", cbuf);

    for(j=jlo; j < jhi; j += 2)
    {
	#ifdef USE_AVX512
		j1 = (j & mask03) + br16[j&15];
	#elif defined(USE_AVX)
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

		fprintf(dbg_file, "j = %8u : %20.5f  %20.5f\n", j, a[j1], a[j1+RE_IM_STRIDE]);
	}

	fclose(dbg_file);	dbg_file = 0x0;
}

/*********************/

/*********************/
