/*******************************************************************************
*                                                                              *
*   (C) 1997-2020 by Ernst W. Mayer.                                           *
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
*   COMPILING AND RUNNING THE PROGRAM: see http://www.mersenneforum.org/mayer/README.html *
******************************************************************************************/
#include "Mlucas.h"
#ifndef imul_macro_h_included
	#error imul_macro.h file not included in build!
#endif

/* 7 Jun 2007: Better to force the builder to invoke these (or not) via compile flag:
#define INCLUDE_TF	// Set to activate trial-factoring code inclusion
#define INCLUDE_PM1	// Set to activate (p-1)-factoring code inclusion
*/

/* Make sure none of the factoring-module-only flags are active: */
#if(defined(P4WORD) || defined(P3WORD) || defined(P2WORD))
	#error multiword exponents only allowed for factor.c built in standalone mode!
#endif
#ifdef FACTOR_STANDALONE
	#error FACTOR_STANDALONE flag only allowed for factor.c built in standalone mode!
#endif

/*********************************************************************************/
/* Globals. Unless specified otherwise, these are declared in Mdata.h:           */
/*********************************************************************************/

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
int NRADICES, RADIX_VEC[10];	// RADIX_VEC[] stores sequence of complex FFT radices used.
#ifdef MULTITHREAD
	uint64 CORE_SET[MAX_CORES>>6];	// Bitmap for user-controlled affinity setting, as specified via the -cpu flag
#endif
int ROE_ITER = 0;		// Iteration of any dangerously high ROE encountered during the current iteration interval.
						// This must be > 0, but make signed to allow sign-flip encoding of retry-fail.
double ROE_VAL = 0.0;	// Value (must be in (0, 0.5)) of dangerously high ROE encountered during the current iteration interval

int USE_SHORT_CY_CHAIN = 0;

int ITERS_BETWEEN_CHECKPOINTS;	/* number of iterations between checkpoints */
int ITERS_BETWEEN_GCHECK_UPDATES = 1000;	// iterations between Gerbicz-checkproduct updates
int ITERS_BETWEEN_GCHECKS     = 1000000;	// #iterations between Gerbicz-checksum residue-integrity checks

char ESTRING[STR_MAX_LEN];	/* Exponent in string form */
char PSTRING[STR_MAX_LEN];	/* Number being tested in string form, typically estring concatenated with several other descriptors, e.g. strcat("M",estring) */

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
char cbuf[STR_MAX_LEN],cstr[STR_MAX_LEN];
char in_line[STR_MAX_LEN];
char *char_addr;
int char_offset;
FILE *fp, *fq;

/* File access mode is a 3-character string, with the last of these being a mandatory null terminator: */
char FILE_ACCESS_MODE[3] = {'x','y','\0'};

/* Matrix of supported moduli/test-types - init'ed in Mlucas_init */
int ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_DIM][TEST_TYPE_DIM];

/****** Allocate storage for Globals (externs): ******/

// Shift count and auxiliary arrays used to support rotated-residue computations:
uint64 RES_SHIFT = 0xFFFFFFFFFFFFFFFFull;	// 0 is a valid value here, so init to UINT64_MAX, which value is treated as "uninited"
uint64 GCHECK_SHIFT = 0ull;
uint64 *BIGWORD_BITMAP = 0x0;
uint32 *BIGWORD_NBITS = 0x0;

// For PRP tests, the base. For Pépin tests via the "FermatTest" worktype, the base defaults to 3; to do a
// Pépin test to another base, the more-general PRP-worktype must be specified with appropriate parameters.
uint32 PRP_BASE = 0;
uint64 *BASE_MULTIPLIER_BITS = 0x0;	// Runtime-allocated bitwise multiply-by-base array

/* These should all be set to a valid (nonzero) value at the time the appropriate test is begun */
uint32 TEST_TYPE		= 0;
uint32 MODULUS_TYPE		= 0;
uint32 TRANSFORM_TYPE	= 0;

const char HOMEPAGE  [] = "http://www.mersenneforum.org/mayer/README.html";

/* Program version with patch suffix:

For the foreseeable future, version numbers will be of form x.yz, with x a 1-digit integer,
y an int having 3 digits or less, and z an optional alphabetic patch suffix. We choose these such that
retiming is only necessary between versions that differ in x or in the leading digit of y,
i.e. we'd need to regenerate .cfg files in going from version 2.9 to 3.0 or from 3.0 to 3.141,
but not between versions 3.0x and 3.0g or between 3.1 and 3.141.

This reduces the "retime?" decision in the cfgNeedsUpdating() function to a simple comparison
of the leading 3 characters of the two version strings in question.
*/
// Dec 2014: After many years frozen at 3.0x (since x86 code wildly uncompetitive with Prime95),
//           now that I'm < 2x slower (on Haswell+) resume version numbering according to the scheme:
//				Major index = year - 2000
//				Minor index = release # of that year, zero-indexed.
// As before, a patch suffix of x, y, or z following the above numeric index indicates an [alpha,beta,gamma] (experimental,unstable) code.
const char VERSION   [] = "19.0";			// e.g. Dec 2014 was 2nd release of that year, thus 14.1 = [20]14.[--2]

const char OFILE     [] = "results.txt";	/* ASCII logfile containing FINAL RESULT ONLY for each
											assignment - detailed intermediate results for each assignment
											are written to the exponent-specific STATFILE (see below). */
const char RANGEFILE [] = "worktodo.ini";	/* File containing exponents to be tested. New exponents
											may be appended at will, while the program is running. */

const char LOCAL_INI_FILE[] = "local.ini";	/* File containing user-customizable configuration settings [currently unused] */

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

int INTERACT;

double AME,MME;	/* Avg and Max per-iteration fractional error for a given iteration interval */

// These externs declared in platform.h:
int MAX_THREADS = 0;		/* Max. allowable No. of threads. */
int NTHREADS = 0;			/* actual No. of threads. If multithreading disabled, set = 1. */

uint64 PMIN;		/* minimum exponent allowed */
uint64 PMAX;		/* maximum exponent allowed depends on max. FFT length allowed
					   and will be determined at runtime, via call to given_N_get_maxP(). */

void sig_handler(int signo)
{
	if (signo == SIGINT) {
		fprintf(stderr,"received SIGINT signal.\n");	sprintf(cbuf,"received SIGINT signal.\n");
	} else if(signo == SIGTERM) {
		fprintf(stderr,"received SIGTERM signal.\n");	sprintf(cbuf,"received SIGTERM signal.\n");
	} else if(signo == SIGHUP) {
		fprintf(stderr,"received SIGHUP signal.\n");	sprintf(cbuf,"received SIGHUP signal.\n");
	} else if(signo == SIGALRM) {
		fprintf(stderr,"received SIGALRM signal.\n");	sprintf(cbuf,"received SIGALRM signal.\n");
	} else if(signo == SIGUSR1) {
		fprintf(stderr,"received SIGUSR1 signal.\n");	sprintf(cbuf,"received SIGUSR1 signal.\n");
	} else if(signo == SIGUSR2) {
		fprintf(stderr,"received SIGUSR2 signal.\n");	sprintf(cbuf,"received SIGUSR2 signal.\n");
	}
	// Toggle a global to allow desired code sections to detect signal-received and take appropriate action:
	MLUCAS_KEEP_RUNNING = 0;
}

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
!   See http://www.mersenneforum.org/mayer/README.html for build instructions, recent revision history
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
!    http://www.mersenneforum.org/mayer/README.html
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
	uint32	iterations,	/* Use to store log2[max factor depth] in TF mode */
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
	int i,j,k = 0;
	/* TODO: some of these need to become 64-bit: */
	uint32 dum,findex = 0,ierr = 0,ilo = 0,ihi = 0,iseed,isprime,kblocks = 0,maxiter = 0,n = 0,npad = 0;
	uint64 itmp64, fwd_fft_only = 0ull, s1,s2,s3;	// s1,2,3: Triply-redundant whole-array checksum on b,c-arrays used in the G-check
	uint32 mode_flag = 0;
	/* Exponent of number to be tested - note that for trial-factoring, we represent p
	strictly in string[STR_MAX_LEN] form in this module, only converting it to numeric
	form in the factoring module. For all other types of assignments uint64 should suffice: */
	uint64 p = 0, i1,i2,i3, rmodb,mmodb;
	uint32 nbits_in_p = 0;
	/* Res64 and Selfridge-Hurwitz residues: */
	const uint64 two35m1 = (uint64)0x7FFFFFFFFull, two36m1 = (uint64)0xFFFFFFFFFull;	/* 2^35,36-1 */
	uint64 Res64, Res35m1, Res36m1;
	uint64 factor_k[16], *fac, twop, quo64,rem64;	// factor_k = Array storing factor k's (assumed < 2^64) of any known-factors for PRP-cofactor
	// testing; fac = ptr to each mi64-converted factor input string; twop = 2*[binary exponent of modulus]; qu064,rem64 = 64-bit quotient & remainder
	uint32 fbits, lenf, nfac, nchar;
	uint128 p128,q128,res128;
/*...Known Mersenne prime exponents. This array must be null-terminated.	*/
	// Dec 2018: Including M51, there are (31, 19) p = 1,3 (mod 4), resp., vs (25.8, 24.5) predicted
	// (for p < 10^8) by the Lenstra/Wagstaff heuristic (cf. est_num_mp_in_interval() in util.c):
	const uint32 knowns[] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941
		,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593
		,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,0x0};

/*...What a bunch of characters...	*/
	char *cptr = 0x0;
/*...initialize logicals and factoring parameters...	*/
	int restart = FALSE, start_run = TRUE;
	uint32 bit_depth_done = 0;

#ifdef INCLUDE_TF
	uint32 bit_depth_todo = 0;
	uint64 factor_k_start = 0;
	uint32 factor_pass_start = 0, factor_pass_hi = 0;
	double log2_min_factor = 0, log2_max_factor = 0;
#endif
#ifdef INCLUDE_PM1
	pm1_done = 0;
#endif

/*...allocatable data arrays...	*/
	static int32 nbytes = 0, nalloc = 0;
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
	/*clock_t clock1, clock2;	Moved these to mers_mod_square.c */
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

	// v19: Reset carry-chain length fiddler to default (faster/lower-accuracy) at start of each run:
	USE_SHORT_CY_CHAIN = 0;
	/* Clear out any FFT-radix data that might remain from a just-completed run: */
	for(i = 0; i < 10; i++) {
		RADIX_VEC[i] = 0;
	}
	NRADICES = 0;

	RESTARTFILE[0] = STATFILE[0] = '\0';
	restart=FALSE;

/*  ...If multithreading enabled, set max. # of threads based on # of available (logical) processors,
with the default #threads = 1 and affinity set to logical core 0, unless user overrides those via -nthread or -cpu:
*/
#ifdef MULTITHREAD

  #ifdef USE_OMP
	// OpenMP not currently supported (attempting to build with this #define enabled barfs in
	// preprocessing via #error in platform.h), this is merely placeholder for possible future use:
	ASSERT(HERE, MAX_THREADS = omp_get_num_procs(), "Illegal #Cores value stored in MAX_THREADS");
  #elif(defined(USE_PTHREAD))
	ASSERT(HERE, MAX_THREADS =     get_num_cores(), "Illegal #Cores value stored in MAX_THREADS");
  #else
	#error Unrecognized multithreading model!
  #endif
	// MAX_THREADS based on number of processing cores will most often be a power of 2, but don't assume that.
	ASSERT(HERE, MAX_THREADS > 0,"MAX_THREADS must be > 0");
	ASSERT(HERE, MAX_THREADS <= MAX_CORES,"MAX_THREADS exceeds the MAX_CORES setting in Mdata.h .");

	if(!NTHREADS) {
		NTHREADS = 1;
		fprintf(stderr,"No CPU set or threadcount specified ... running single-threaded.\n");
		// Use the same affinity-setting code here as for the -cpu option, but simply for cores [0:NTHREADS-1]:
		sprintf(cbuf,"0:%d",NTHREADS-1);
		parseAffinityString(cbuf);
	} else if(NTHREADS > MAX_CORES) {
		sprintf(cbuf,"FATAL: NTHREADS = %d exceeds the MAX_CORES setting in Mdata.h = %d\n", NTHREADS, MAX_CORES);
		ASSERT(HERE, 0, cbuf);
	} else {	// In timing-test mode, allow #threads > #cores
		if(NTHREADS > MAX_THREADS) {
			fprintf(stderr,"WARN: NTHREADS = %d exceeds number of cores = %d\n", NTHREADS, MAX_THREADS);
		}
		if(start_run) fprintf(stderr,"NTHREADS = %d\n", NTHREADS);
	}

  #if 0//defined(USE_PTHREAD) && defined(OS_TYPE_MACOSX)

	thread_t thr = mach_thread_self();
	thread_extended_policy_data_t epolicy;
	epolicy.timeshare = FALSE;
	kern_return_t ret = thread_policy_set(
		thr, THREAD_EXTENDED_POLICY,
		(thread_policy_t) &epolicy, THREAD_EXTENDED_POLICY_COUNT);
	if (ret != KERN_SUCCESS) {
		printf("thread_policy_set returned %d", ret);
		exit(-1);
	}

	thread_affinity_policy_data_t apolicy;
	int cpui = MAX_THREADS - 1;	// get cpu mask using sequential thread ID modulo #available cores
	apolicy.affinity_tag = cpui; // set affinity tag

	printf("Setting CPU = %d affinity of main thread, mach_id = %u\n", cpui, thr);

	ret = thread_policy_set(
		thr, THREAD_EXTENDED_POLICY,
		(thread_policy_t) &apolicy, THREAD_EXTENDED_POLICY_COUNT);
	if (ret != KERN_SUCCESS) {
		printf("thread_policy_set returned %d", ret);
		exit(-1);
	}

  #endif

#else

		MAX_THREADS = NTHREADS = 1;

#endif	// #ifdef MULTITHREAD ?

	/* Make number of iterations between checkpoints dependent on #threads -
	don't want excessively frequent savefile writes, at most 1 or 2 an hour is needed:
	*/
	if(NTHREADS > 4)
		ITERS_BETWEEN_CHECKPOINTS = 100000;
	else
		ITERS_BETWEEN_CHECKPOINTS =  10000;

	// v19: If PRP test, make sure Gerbicz-checkproduct interval divides checkpoint-writing one:
	i = ITERS_BETWEEN_GCHECKS;
	j = ITERS_BETWEEN_GCHECK_UPDATES;
	k = ITERS_BETWEEN_CHECKPOINTS;
	ASSERT(HERE, i == j*j, "#iterations between Gerbicz-checksum updates must = sqrt(#iterations between residue-integrity checks)");
	ASSERT(HERE, i%k == 0 && k%j == 0, "G-checkproduct update interval must divide savefile-update one, which must divide the G-check interval");

	// Alloc bitwise multiply-by-base array, needed to support PRP testing:
	if(!BASE_MULTIPLIER_BITS) {
		i = ((ITERS_BETWEEN_CHECKPOINTS+63) >> 6) + 1;	// Add 1 pad element in case compiler does not 64-bit align
		BASE_MULTIPLIER_BITS = ALLOC_UINT64(BASE_MULTIPLIER_BITS, i);	if(!BASE_MULTIPLIER_BITS){ sprintf(cbuf, "FATAL: unable to allocate BASE_MULTIPLIER_BITS array in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		BASE_MULTIPLIER_BITS = ALIGN_UINT64(BASE_MULTIPLIER_BITS);
		ASSERT(HERE, ((long)BASE_MULTIPLIER_BITS & 63) == 0x0,"BASE_MULTIPLIER_BITS[] not aligned on 64-byte boundary!");
	}

/*  ...look for worktodo.ini file...	*/
	fp = 0x0;
	if (start_run && (!exponent || (exponent!=0 && fft_length!=0 && iterations==0)))
	{
		fprintf(stderr," looking for %s file...\n", RANGEFILE);
		fp = mlucas_fopen(RANGEFILE, "r");
	}

	/****************************************/
	/*...Automated exponent-dispatch mode...*/
	/****************************************/
	if (fp) {
		if (start_run) fprintf(stderr," worktodo.ini file found...checking next exponent in range...\n");

		/* 7/06/04: need to handle both exponent-only (i.e. pre-networked Mlucas) and Prime95 ini file formats:
		fscanf(fp,"%u", &p);
		fclose(fp);	fp = 0x0;
		*/
		/* Read first line of worktodo.ini file into 1K character array... */
		fgets(in_line, STR_MAX_LEN, fp);

		/* Skip any whitespace at beginning of the line: */
		char_addr = in_line;
		i = 0;
		while(isspace(in_line[i])) { ++i; }

		/* If leading char is numeric, assume it's Mersenne-exponent-only (legacy) format: */
		if(isdigit(in_line[i]))
		{
			char_addr += i;
			goto GET_EXPO;
		}

		// Otherwise assume Prime95-style ini file format, with a possible modulus-specific leading keyword;
		// Default "Test=" means Mersenne, unless "Test" preceded by an explicit modulus-type string:
		MODULUS_TYPE = MODULUS_TYPE_MERSENNE;
		/*
		Re. the recently-added-to-Primenet PRP assignment type, On Dec 19, 2017, at 5:07 PM, George Woltman wrote:

			In "PRP= xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx,1,2,75869377,-1,75,0,3,4"

			The first four [numeric] values are k,b,n,c as in k*b^n+c
			75 is how far factored
			0 is the number of PRP tests that will be saved if P-1 is done and finds a factor
			3 is the PRP base		[used for the first-time test]
			4 is the residue type	[used for the first-time test]

			You should only need to support type-1 and type-5, which are identical if there are no cofactors (in JSON return residue type 1).

		Example:	PRP=C42540C352E54E906108D48FA5D89488,1,2,80340397,-1,75,0,3,1
		means a PRP test of 1.2^80340397-1 = M80340397, TFed to 2^75 (ignore), 0 saved tests (ignore), PRP base 3 and type 1.
		EWM: The above are PRP-DCs ... first-time tests lack the last 2 numeric args above, and default to base = 3, residue type 1.
		*/
		if((char_addr = strstr(in_line, "PRP")) != 0)
		{
			TEST_TYPE = TEST_TYPE_PRP;
			/* No dicking around with whitespace allowed here: */
			if(char_addr != in_line)
				ASSERT(HERE, 0,"Assignment type specifier must occur at beginning of worktodo.ini file entry, with no whitespace!");
			else
				char_addr += 3;

			ASSERT(HERE, (char_addr = strstr(char_addr, "=")) != 0x0,"Expected '=' not found in assignment-specifying line!");
			char_addr++;
			while(isspace(*char_addr)) { ++char_addr; }	// Skip any whitespace following the equals sign
			ASSERT(HERE, is_hex_string(char_addr, 32), "Expect a 32-hex-digit PrimeNet v5 assignment ID following the work type specifier!");
			cptr = char_addr + 32;
			/*
			Since my code only supports Mersenne and Fermat-number moduli, check to ensure
			that k = 1, b = 2, c = (+1 or -1) and (n prime if c = -1, n  =  2^m if c = +1):
			*/
			ASSERT(HERE, (char_addr = strstr(cptr, ",")) != 0x0,"Expected ',' not found in assignment-specifying line!");
			i = (int)strtol(char_addr+1, &cptr, 10);	ASSERT(HERE, i == 1,"In modulus expression m = k*b^n+c, only k = 1 currently supported!");
			ASSERT(HERE, (char_addr = strstr(cptr, ",")) != 0x0,"Expected ',' not found in assignment-specifying line!");
			i = (int)strtol(char_addr+1, &cptr, 10);	ASSERT(HERE, i == 2,"In modulus expression m = k*b^n+c, only b = 2 currently supported!");
			ASSERT(HERE, (char_addr = strstr(cptr, ",")) != 0x0,"Expected ',' not found in assignment-specifying line!");
			p =     strtoul(char_addr+1, &cptr, 10);	ASSERT(HERE, p < PMAX,"Exponent n in modulus expression m = k*b^n+c exceeds limit! (Suggest checking for unsigned overflow.)");
			ASSERT(HERE, (char_addr = strstr(cptr, ",")) != 0x0,"Expected ',' not found in assignment-specifying line!");
			i = (int)strtol(char_addr+1, &cptr, 10);	ASSERT(HERE, ABS(i) == 1,"In modulus expression m = k*b^n+c, only c = +-1 currently supported!");
			if(i == -1) {
				MODULUS_TYPE = MODULUS_TYPE_MERSENNE;
				ASSERT(HERE, !(p >> 32), "p must be 32-bit or less!");	// Future versions will need to loosen this p < 2^32 restriction
				ASSERT(HERE, isPRP((uint32)p),"Mersenne exponent must be prime!");
			} else if(i == 1) {
				MODULUS_TYPE = MODULUS_TYPE_FERMAT;
				ASSERT(HERE, isPow2_64(p),"Fermat exponent must be a power of 2!");
				p = trailz64(p);	// Code below requires p contain Fermat-number *index* (not exponent), just as in a FermatTest assignment
			}
			/* Already asserted that there is no 'else' needed */
			sprintf(ESTRING,"%llu",p);	// Need to init this for filenaming code
			// Skip next 2 entries in in-string, how-far-factored and "# of PRP tests that will be saved if P-1 is done and finds a factor":
			ASSERT(HERE, (cptr      = strstr(cptr, ",")) != 0x0,"Expected ',' not found in assignment-specifying line!");	cptr++;
			ASSERT(HERE, (cptr      = strstr(cptr, ",")) != 0x0,"Expected ',' not found in assignment-specifying line!");	cptr++;
			if((char_addr = strstr(cptr, ",")) == 0x0) {	// First-time PRP test:
				PRP_BASE = 3;
				TEST_TYPE = TEST_TYPE_PRP;
			} else {	// PRP double-check:
				// NB: Hit a gcc compiler bug (which left i = 0 for e.g. char_addr = ", 3 ,...") using -O0 here ... clang compiled correctly, as did gcc -O1:
				i = (int)strtol(char_addr+1, &cptr, 10);	ASSERT(HERE, i > 1 && IS_ODD(i),"PRP-test base must be > 1 and odd!");
				PRP_BASE = i;
				ASSERT(HERE, (char_addr = strstr(cptr, ",")) != 0x0,"Expected ',' not found in assignment-specifying line!");
				i = (int)strtol(char_addr+1, &cptr, 10);	ASSERT(HERE, i == 1 || i == 5,"Only PRP-tests of type 1 (PRP-only) and type 5 (PRP and subsequent cofactor-PRP check) supported!");
				if(i == 1 && i == 5)
					TEST_TYPE = TEST_TYPE_PRP;
			}
		/*
			else if(i == 5) {	// To-Do!
				TEST_TYPE = TEST_TYPE_PRPCOFACTOR;
				// Read known factors from next line of worktodo, in double-quote-delimited "factor1,factor2,..." format:
				ASSERT(HERE, fgets(in_line, STR_MAX_LEN, fp) != 0x0, "Hit EOF while attempting to read known factors from next line of worktodo!");
				ASSERT(HERE, in_line[0] == '\"',"Known-factors line of worktodo must consist of a comma-separated list of such enclosed in double-quotes!");
				cptr = in_line+1;	nfac = 0;	fac = 0x0;
				while((char_addr = strstr(cptr, "\"")) != 0x0) {
					nchar = char_addr - cptr;
					strncpy(cbuf,cptr,nchar);	cbuf[nchar] = '\0';	// Extract current-factor-as-string into cbuf
					// Convert stringified factor f to mi64 form:
					fac = convert_base10_char_mi64(cbuf, &lenf);	// This does the mem-alloc for us
					ASSERT(HERE, lenf > 0, "Error converting known-factor string!");
					ASSERT(HERE, lenf < 3, "known-factor out of range, must be < 2^128!");
					fbits = (lenf<<6) - mi64_leadz(fac, lenf);
					// Extract the factor k and make sure k < 2^64:
					if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
						twop = p<<1;
						mi64_div(fac,&twop,lenf,1,&quo64,&rem64);
					} else {
						twop = 1ull<<(p+1);
						mi64_div(fac,&twop,lenf,1,&quo64,&rem64);
					}
					ASSERT(HERE, rem64 == 1ull,"Factor not of required form!");
					factor_k[nfac++] = quo64;
					// Verify that it is a factor via binary modpow - the latter routine auto-determines the modulus type based on the exponent:
					if(lenf == 1)
						ASSERT(HERE, twopmodq64(twop>>1,fac[0]) == 1ull,"Factor-check fails!");
					else {
						p128.d0 = twop>>1;	p128.d1 = 0ull;
						q128.d0 = fac[0];	q128.d1 = fac[1];
						res128 = twopmodq128(p128,q128);
						ASSERT(HERE, CMPEQ128(res128,ONE128),"Factor-check fails!");
					}
					if(char_addr[0] == '\"' && char_addr[1] == ',' && char_addr[2] == '\"')
						cptr = char_addr+3;
					else if(char_addr[0] == '\"' && char_addr[1] != ',')
						break;
					else
						ASSERT(HERE,0,"Unrecognized token sequence in parsing known-factors line of worktodo file!");
				}
				ASSERT(HERE, nfac != 0,"Must specify at least one known factor!");
			}
		*/
			goto PRP;
		}
		else if((char_addr = strstr(in_line, "Fermat")) != 0)
		{
			if(char_addr != in_line)
				ASSERT(HERE, 0,"Modulus type specifier must occur at beginning of worktodo.ini file entry, with no whitespace!");
			else
				char_addr += 6;
			/* Look for comma following the modulus keyword and position next-keyword search right after it: */
			if(!STREQN(char_addr,",",1))
				ASSERT(HERE, 0,"Expected ',' not found in input following modulus type specifier!");
			else
				char_addr++;

			MODULUS_TYPE = MODULUS_TYPE_FERMAT;
			PRP_BASE = 2;	// v18: If use residue shift in context of Pépin test, need this whenever the 'shift = 2*shift + random[0,1]' update gets a 1-bit in the random slot
		}
		/* "Mersenne" is the default and hence not required, but allow it: */
		else if((char_addr = strstr(in_line, "Mersenne")) != 0)
		{
			if(char_addr != in_line)
				ASSERT(HERE, 0,"Modulus type specifier must occur at beginning of worktodo.ini file entry, with no whitespace!");
			else
				char_addr += 8;
			/* Look for comma following the modulus keyword and position next-keyword search right after it: */
			if(!STREQN(char_addr,",",1))
				ASSERT(HERE, 0,"Expected ',' not found in input following modulus type specifier!");
			else
				char_addr++;
		}
		else
			char_addr = in_line;

		if(strstr(char_addr, "Test"))
		{
			TEST_TYPE = TEST_TYPE_PRIMALITY;
			char_offset =  4;
		}
		else if(strstr(char_addr, "DoubleCheck"))
		{
			TEST_TYPE = TEST_TYPE_PRIMALITY;
			char_offset = 11;
		}
	#ifdef INCLUDE_TF
		else if(strstr(char_addr, "Factor"))
		{
			TEST_TYPE = TEST_TYPE_TRIALFACTORING;
			char_offset =  6;
		}
	#endif
	#ifdef INCLUDE_PM1
		else if(strstr(char_addr, "Pminus1") || strstr(char_addr, "PMinus1"))
		{
			TEST_TYPE = TEST_TYPE_PMINUS1;
			char_offset =  7;
		}
	#endif
	#ifdef INCLUDE_ECM
		else if(strstr(char_addr, "ECM"))
		{
			TEST_TYPE = TEST_TYPE_ECM;
			char_offset =  3;
		}
	#endif
		else
		{
			sprintf(cbuf, "ERROR: Unrecognized/Unsupported option. The ini file entry was %s\n", char_addr);
			fprintf(stderr,"%s",cbuf);
			fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf); fclose(fp);	fp = 0x0; }
			ASSERT(HERE, 0, "");
		}

		char_addr += char_offset;

		ASSERT(HERE, (char_addr = strstr(char_addr, "=")) != 0x0,"Expected '=' not found in assignment-specifying line!");
		char_addr++;
		/* Skip any whitespace following the equals sign:*/
		while(isspace(*char_addr))
		{
			++char_addr;
		}
		/* Check for a 32-hex-digit PrimeNet v5 assignment ID preceding the exponent: */
		if(is_hex_string(char_addr, 32))
		{
			char_addr += 33;
		}

	GET_EXPO:
		/* This takes care of the number-to-char conversion and leading-whitespace-removal
		in one step - use PSTRING for temporary storage here:
		*/

	/* Copy the exponent to ESTRING and null-terminate: */
		ASSERT(HERE, isdigit(*char_addr),"False result for isdigit(*char_addr)");
		i = 0;
		while(isdigit(*char_addr))
		{
			ESTRING[i++] = *char_addr;
			++char_addr;
		}
		ESTRING[i++] = '\0';
		p = convert_base10_char_uint64(ESTRING);

	PRP:	// In PRP-test case, have already read the exponent from the worktodo line
		/* Special case of user forcing a non-default FFT length for an exponent in the worktodo.ini file: */
		if(exponent)
		{
  			if((p != exponent))// || (MODULUS_TYPE != MODULUS_TYPE_MERSENNE))	15. Oct 2012: Need same flexibility for Fermat numbers (e.g. F27 @ 7168k) as for Mersennes, so disable modulus-type part of conditional
				ASSERT(HERE, 0,"User-supplied exponent and FFT-length for full-length test requires an exponent-matching 'Test=[exponent]' or 'DoubleCheck=[exponent]' worktodo.ini entry!");
		}

		/* Check #bits in the power-of-2 exponent vs. the allowed maximum: */
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			nbits_in_p = 64 - leadz64(p);
		}
		/* If it's a Fermat number, need to check size of 2^ESTRING: */
		else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			findex = (uint32)p;
			if(findex <= MAX_PRIMALITY_TEST_BITS)
				p = (uint64)1 << findex;
			else
				ASSERT(HERE, 0,"nbits_in_p <= MAX_PRIMALITY_TEST_BITS");

			/* For purposes of the bits-in-p limit, treat 2^findex as having
			(findex) rather than (findex+1) bits: */
			nbits_in_p = findex;
		}
		else
			ASSERT(HERE, 0,"MODULUS_TYPE unknown!");

		ASSERT(HERE, nbits_in_p <= MAX_EXPO_BITS,"Require nbits_in_p <= MAX_EXPO_BITS");

	#ifdef INCLUDE_TF

	  INIT_TF:

		/* If nbits_in_p > MAX_PRIMALITY_TEST_BITS, it better be a TF run: */
		if(TEST_TYPE == TEST_TYPE_TRIALFACTORING)
		{
			/* Currently TF only supported for Mersennes: */
			if(MODULUS_TYPE != MODULUS_TYPE_MERSENNE)
			{
				sprintf(cbuf, "ERROR: Trial-factoring Currently only supported for Mersenne numbers. The ini file entry was %s\n", in_line);
									           fprintf(stderr,"%s",cbuf);
				fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				goto GET_NEXT_ASSIGNMENT;
			}

			/* For now, always start at k = 1: */
			log2_min_factor = 0.0;
			log2_max_factor = get_default_factoring_depth(p);
			ASSERT(HERE, log2_max_factor <= MAX_FACT_BITS, "log2_max_factor > MAX_FACT_BITS!");

			/* Field following the exponent is the already-factored-to depth: if none found, use defaults. */
			char_addr = strstr(char_addr, ",");
			if(char_addr)
			{
				char_addr++;
				/* Convert the ensuing numeric digits to ulong: */
				bit_depth_done = strtoul(char_addr, (char**)NULL, 10);

				/* Specified already-factored-to depth is larger than default factor-to depth - no more factoring to be done. */
				if(bit_depth_done > log2_max_factor)
				{
					sprintf(cbuf, "INFO: the specified already-factored-to depth of %u bits exceeds the default %10.4f bits - no more factoring to be done.\n", bit_depth_done, log2_max_factor);
										           fprintf(stderr,"%s",cbuf);
					fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					goto GET_NEXT_ASSIGNMENT;
				}
			}
			/* Mode to allow user to override normal automated-test TF default limits:
			If a second ,{#bits} argument is present in the input line, then this user-set
			desired factoring depth overrides the normal default for the exponent in question.
			In this mode, warn if TF-to depth greater than automated-mode default, but allow: */
			if(char_addr)char_addr = strstr(char_addr, ",");
			if(char_addr)
			{
				char_addr++;
				bit_depth_todo = strtoul(char_addr, (char**)NULL, 10);
				if(bit_depth_todo > MAX_FACT_BITS)
				{
					sprintf(cbuf, "ERROR: factor-to bit_depth of %u > max. allowed of %u. The ini file entry was %s\n", bit_depth_done, MAX_FACT_BITS, in_line);
										           fprintf(stderr,"%s",cbuf);
					fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					goto GET_NEXT_ASSIGNMENT;
				}
				else if(bit_depth_todo <= bit_depth_done)
				{
					sprintf(cbuf, "ERROR: factor-to bit_depth of %u < already-done depth of %u. The ini file entry was %s\n", bit_depth_done, bit_depth_done, in_line);
										           fprintf(stderr,"%s",cbuf);
					fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					goto GET_NEXT_ASSIGNMENT;
				}
				else if(bit_depth_todo > log2_max_factor)
				{
					sprintf(cbuf, "WARN: the specified factor-to depth of %u bits exceeds the default %10.4f bits - I hope you know what you're doing.\n", bit_depth_done, log2_max_factor);
					log2_max_factor = bit_depth_todo;
										           fprintf(stderr,"%s",cbuf);
					fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				}
				log2_max_factor = bit_depth_todo;
			}
		}
		else if(nbits_in_p > MAX_PRIMALITY_TEST_BITS)
		{
			sprintf(cbuf, "ERROR: Inputs this large only permitted for trial-factoring. The ini file entry was %s\n", in_line);
								           fprintf(stderr,"%s",cbuf);
			fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
			goto GET_NEXT_ASSIGNMENT;
		}
	#endif 	// INCLUDE_TF

		/* If "Test..." or "DoubleCheck", check for bit_depth_done and pm1_done fields following the = sign:
		if present and there is still factoring remaining to be done, modify the assignment type appropriately: */
		if(TEST_TYPE == TEST_TYPE_PRIMALITY)
		{
			/* bit_depth_done: */
			bit_depth_done = 0xffffffff;	/* Only check if there's an appropriate entry in the input line */
			char_addr = strstr(char_addr, ",");
			if(char_addr)
			{
				char_addr++;
				/* FFT length specified? (E.g. if user wants non-default) */
				if(strstr(char_addr, "fftlen"))
				{
					char_addr = strstr(char_addr+6, "=");
					if(char_addr)
					{
						char_addr++;
						/* Convert the ensuing numeric digits to ulong: */
						fft_length = strtoul(char_addr, (char**)NULL, 10);
						kblocks = get_default_fft_length(p);	/* Default FFT length for this exponent */
						/* Check that user-specified FFT length is >= default, or that p <= 1.01*(max exponent for given length): */
						if( (kblocks > fft_length) && (p > 1.01*given_N_get_maxP(fft_length<<10)) )
						{
							snprintf(cbuf,STR_MAX_LEN,"ERROR: Illegal 'fftlen = ' argument - suggested FFT length for this p = %u. The ini file entry was %s\n", kblocks, in_line);
														   fprintf(stderr,"%s",cbuf);
							fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
							goto GET_NEXT_ASSIGNMENT;
						}
					}
				}
				else
				{
					/* Convert the ensuing numeric digits to ulong: */
					bit_depth_done = strtoul(char_addr, (char**)NULL, 10);
					if(bit_depth_done > MAX_FACT_BITS)
					{
						snprintf(cbuf,STR_MAX_LEN,"ERROR: bit_depth_done of %u > max. allowed of %u. The ini file entry was %s\n", bit_depth_done, MAX_FACT_BITS, in_line);
													   fprintf(stderr,"%s",cbuf);
						fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
						goto GET_NEXT_ASSIGNMENT;
					}

					/* If bit_depth_done less than default TF depth for this exponent
					and this platform is approved for trial factoring, switch task type to TF:
					*/
				#ifdef INCLUDE_TF
					log2_max_factor = get_default_factoring_depth(p);
					if(bit_depth_done < log2_max_factor)
					{
						TEST_TYPE = TEST_TYPE_TRIALFACTORING;
						/* For now, always start at k = 1: */
						log2_min_factor = 0.0;
					}
					else if(bit_depth_done > log2_max_factor)
					{
						sprintf(cbuf, "WARN: the specified already-factored-to depth of %u bits exceeds the default %10.4f bits - no more factoring to be done.\n", bit_depth_done, log2_max_factor);
													   fprintf(stderr,"%s",cbuf);
						fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					}
				  #ifdef INCLUDE_PM1
					else
					{
						/* pm1_done: */
						char_addr = strstr(char_addr, ",");
						if(char_addr)
						{
							char_addr++;
							/* Convert the ensuing numeric digits to ulong: */
							pm1_done = strtoul(char_addr, (char**)NULL, 10);
							/* This flag should be 0 or 1: */
							if(pm1_done == 1)
							{
								/* no-op */
							}
							else if(pm1_done == 0)
							{
								TEST_TYPE = TEST_TYPE_PMINUS1;
							}
							else
							{
								sprintf(cbuf, "ERROR: Illegal p-1 status flag value of %lu (should be 0 or 1). The ini file entry was %s\n", pm1_done, in_line);
															   fprintf(stderr,"%s",cbuf);
								fp = mlucas_fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
								goto GET_NEXT_ASSIGNMENT;
							}
						}
					}
				  #endif
				#endif
				}
			}	/* if(char_addr) */
		}	/* if(TEST_TYPE == TEST_TYPE_PRIMALITY) */
	}
	/****************************************/
	/* Self-test or User-supplied exponent: */
	/****************************************/
	else if(exponent != 0)	/* elseif((found RANGEFILE) == FALSE) */
	{
		p = exponent;

		/* This takes care of the number-to-char conversion and leading-whitespace-removal
		in one step - use PSTRING for temporary storage here:
		*/
		strcpy(ESTRING, &PSTRING[convert_uint64_base10_char(PSTRING, p)]);

		/* If it's a Fermat number, get the real exponent: */
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			findex = (uint32)p;
			if(findex <= MAX_PRIMALITY_TEST_BITS)
				p = (uint64)1 << findex;
			else
				ASSERT(HERE, 0,"nbits_in_p <= MAX_PRIMALITY_TEST_BITS");

			/* For purposes of the bits-in-p limit, treat 2^findex as having
			(findex) rather than (findex+1) bits: */
			nbits_in_p = findex;
		}
		else
			nbits_in_p = 64 - leadz64(p);

		INTERACT=TRUE;

		ASSERT(HERE,TEST_TYPE,"TEST_TYPE not set!");
		ASSERT(HERE,TEST_TYPE <= TEST_TYPE_MAX,"TEST_TYPE out of range!");

		/* If nbits_in_p > MAX_PRIMALITY_TEST_BITS, it better be a TF run: */
		if(TEST_TYPE == TEST_TYPE_TRIALFACTORING)
		{
		#ifdef INCLUDE_TF
			/* Currently TF only supported for Mersennes: */
			ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE), "Trial-factoring Currently only supported for Mersenne numbers");

			/* For now, always start at k = 1: */
			log2_min_factor = 0.0;

			if(iterations)
			{
				log2_max_factor = iterations;
				iterations = 0;
			}
			else
				log2_max_factor = get_default_factoring_depth(p);

			ASSERT(HERE, log2_max_factor >=             0, "log2_max_factor must be positive!");
			ASSERT(HERE, log2_max_factor <= MAX_FACT_BITS, "log2_max_factor exceeds MAX_FACT_BITS!");
		#else
			ASSERT(HERE, 0, "Trial-factoring not supported for this build/platform.");
		#endif
		}
		else	/* Primality or PRP test */
		{
		/*	fprintf(stderr, "P = %u, nbits_in_p = %d\n",p,nbits_in_p);	*/
			ASSERT(HERE, nbits_in_p <= MAX_PRIMALITY_TEST_BITS, "Inputs this large only permitted for trial-factoring.");

			ASSERT(HERE,iterations != 0,"Timing test with User-supplied exponent requires number of iterations to be specified via the -iters flag!");

			if((int)iterations <= 0)
			{
				fprintf(stderr, " Specified %u self-test iterations : must be > 0.\n", iterations);
				return ERR_TESTITERS_OUTOFRANGE;
			}

			if(iterations > MAX_SELFTEST_ITERS)
			{
				fprintf(stderr, " Specified %u iterations exceeds self-test limit of %u.\n", iterations, MAX_SELFTEST_ITERS);
				return ERR_TESTITERS_OUTOFRANGE;
			}

			timing_test_iters = iterations;
		}
	}
	else
	{
		ASSERT(HERE, 0, "Illegal combination of command args - please run with -h to see help menu. Note that if you are trying to run a single-FFT-length self-test, you *must* explicitly specify the iteration count, e.g. './Mlucas -fftlen 7168 -iters [100|1000|10000] [-cpu [args]]'");
	}
	/* endif(found RANGEFILE?)	*/

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

	ASSERT(HERE,TEST_TYPE,"TEST_TYPE not set!");
	ASSERT(HERE,TEST_TYPE <= TEST_TYPE_MAX,"TEST_TYPE out of range!");

	/* Fom this point onward the first character of restart filenames is context-dependent: */
#ifdef INCLUDE_TF
	if(TEST_TYPE == TEST_TYPE_TRIALFACTORING)
	{
		/* Do any needed Trial-Factoring and then update the worktodo.ini file appropriately: */
		/* Let the factoring module handle the restart-file processing in this case. */
		if(!ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE][TEST_TYPE_TRIALFACTORING])
		{
			sprintf(cbuf, "TEST_TYPE_TRIALFACTORING with MODULUS_TYPE = %u not supported!\n", MODULUS_TYPE);
			ASSERT(HERE, 0, cbuf);
		}

		factor(ESTRING, log2_min_factor, log2_max_factor);
		goto GET_NEXT_ASSIGNMENT;
	} else
#endif
	if(TEST_TYPE > TEST_TYPE_MAX)
	{
		ASSERT(HERE, 0,"FATAL: Unrecognized assignment type in savefile processing.\n");
	}
	/* endif(TEST_TYPE == ...) */

/********************* Primality Test: ***********************************************/

	if(p < PMIN) {
		fprintf(stderr, " p must be at least %llu.\n",PMIN);
		return ERR_EXPONENT_ILLEGAL;
	} else if(p > PMAX) {
		fprintf(stderr, " p must be no greater than %llu.\n",PMAX);
		return ERR_EXPONENT_ILLEGAL;
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		ASSERT(HERE, !(p >> 32), "p must be 32-bit or less!");	/* Future versions will need to loosen this p < 2^32 restriction: */
		/*  ...make sure p is prime...	*/
		/* TODO: once factoring code is hooked in, make its small-primes table
		a global and init it at the start of execution, then use it to trial-divide
		p up to sqrt(p):
		*/
		/* Allow non-prime-exponents for p-1 test, but not for LL: */
		if(!isPRP((uint32)p))
		{
			fprintf(stderr, " p is not prime!\n");
//			return ERR_EXPONENT_ILLEGAL;
		}

		TRANSFORM_TYPE = REAL_WRAPPER;
		snprintf(PSTRING,STR_MAX_LEN, "M%s", ESTRING);
		/*
		v19: Unlike standard Fermat-PRP test with its [p-2 squarings-and-mul-by-base, followed by 1 final pure-squaring],
		Gerbicz check requires a pure sequence of squarings. RG: simply replace the standard Fermat-PRP test with a modified
		one where we add 2 to the computed power and check whether the result == x0^2 (mod n).
		E.g. for n = 2^p-1, instead of the standard base-x0 Fermat PRP test,
			x0^(n-1) = x0^(2^p-2) == 1 (mod n)
		we check whether
			x0^(n+1) = x0^(2^p) == x0^2 (mod n).
		That means start with initial seed x0 and do p mod-squarings:
		*/
		if(TEST_TYPE == TEST_TYPE_PRIMALITY)
			maxiter = (uint32)p-2;	// LL-test
		else if(TEST_TYPE == TEST_TYPE_PRP)
			maxiter = (uint32)p;	// Fermat-PRP test modified as described in above commentary
		else
			ASSERT(HERE,0,"Unsupported test type! (Neither LL nor PRP)");
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
	#ifdef USE_ARM_V8_SIMD
		ASSERT(HERE, 0, "ARMv8 SIMD builds do not support Fermat-number testing!");
	#endif
		ASSERT(HERE, findex >=14, "Fermat number index must be at least 14!\n");
		ASSERT(HERE, findex < 64, "Fermat number index must be < 64!\n"       );

		// This takes care of the number-to-char conversion and leading-whitespace-removal
		// in one step - use PSTRING for temporary storage here:
		strcpy(ESTRING, &PSTRING[convert_uint64_base10_char(PSTRING, (uint64)findex)]);
		ASSERT(HERE, (p >> findex) == 1,"Require (p >> findex) == 1");
		TRANSFORM_TYPE = RIGHT_ANGLE;
		sprintf(PSTRING, "F%u", findex);
		maxiter = (uint32)p-1;
		// v18: If nonzero (but 'inited', i.e. != -1) residue shift, init BASE_MULTIPLIER_BITS to random bits:
		if((RES_SHIFT != -1ull) && !PRP_BASE) {
			PRP_BASE = 2;	ASSERT(HERE, RES_SHIFT < p,"User-set initial shift must be less than binary exponent!");
		} else if(RES_SHIFT == -1ull) {
			PRP_BASE = 2;	RES_SHIFT = rng_isaac_rand() % p;
		}
		j = ((ITERS_BETWEEN_CHECKPOINTS+63) >> 6);
		itmp64 = -(RES_SHIFT != 0ull);
		// "Prime' the RNG by doing [RES_SHIFT] calls:
		for(i = 0; i < RES_SHIFT; i++) { rng_isaac_rand(); }
		for(i = 0; i < j; i++) { BASE_MULTIPLIER_BITS[i] = itmp64 & rng_isaac_rand(); }
	//	printf("random-bit-array[0] = 0x%016llX\n",BASE_MULTIPLIER_BITS[0]);
	}
	else
	{
		ASSERT(HERE, 0,"Unknown Self-Test Modulus Type!");
	}	/* endif(MODULUS_TYPE) */

	/* In production-run (INTERACT = False) mode, allow command-line-forced FFT lengths which are at most
	"one size too large" relative to the default length for the exponent in question. Supported lengths
	are of the form [8,9,10,11,12,13,14,15]*2^k, so base our acceptance threshold on the largest adjacent-
	FFT-lengths ratio permitted under this schema, 9/8:
	*/
	kblocks = get_default_fft_length(p);
	if(!fft_length || (!INTERACT && 8*fft_length > 9*kblocks)) {
		if(!kblocks) {
			fprintf(stderr,"ERROR detected in get_default_fft_length for p = %u.\n", (uint32)p);
			return ERR_FFTLENGTH_ILLEGAL;
		}
	} else {
		kblocks = fft_length;
		if((int)kblocks <= 0) {
			fprintf(stderr, "ERROR: %d K must be >= 0.\n",kblocks);
			return ERR_FFTLENGTH_ILLEGAL;
		}
	}

	/*...calculate the unpadded runlength...	*/
	n = kblocks << 10;

	/* Make sure the FFT length and radix set are supported: */
	if(INTERACT)
	{
		if(radix_set < 0)
		{
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

		if(timing_test_iters > maxiter)
		{
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
		if(!dum)
		{
			/* Need to run a timing self-test at this FFT length before proceeding: */
			return ERR_RUN_SELFTEST_FORLENGTH + (kblocks << 8);
		}
		else if(dum != kblocks)
		{
			/* If return value != kblocks, extract the FFT length it encodes: */
			i = extractFFTlengthFrom32Bit(dum);

			/* Only allow lengths that are <= 2x default */
			if( !(i >= kblocks && i <= (kblocks<<1) ) )
			{
				sprintf(cbuf,"Call to get_preferred_fft_radix returns out-of-range FFT length: asked for %u, returned %u, packed value= 0x%8X\n", kblocks, i, dum);
				ASSERT(HERE, 0, cbuf);
			}
			else	/* If length acceptable, extract the FFT-radix data encoded and populate the NRADICES and RADIX_VEC[] globals */
			{
				extractFFTradicesFrom32Bit(dum);
				kblocks = i;
				/* Make sure the FFT length is supported: */
				if(get_fft_radices(kblocks, 0, 0x0, 0x0, 0) != 0)
				{
					ASSERT(HERE, get_fft_radices(kblocks, 0, 0x0, 0x0, 0) == ERR_FFTLENGTH_ILLEGAL, "Unexpected return value for get_fft_radices()");
					sprintf(cbuf, "ERROR: length %d = %d K not available.\n",n,kblocks);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}
		}
	}

	/*...calculate the unpadded runlength...	*/
	n = kblocks << 10;

	if(kblocks != (n>>10))
	{
		fprintf(stderr, "ERROR: length %d K overflows N = 1024*K.\n",kblocks);
		return ERR_FFTLENGTH_ILLEGAL;
	}

	/* If specified FFT length smaller than default for this exponent [only an issue for user-specified FFT length],
	print a warning if the p/pmax ratio > 1 to an acceptably small degree; error out if the ratio is unreasonably > 1:
	*/
	uint64 pmax_rec = given_N_get_maxP(n);	double exp_ratio =  (double)p/pmax_rec;
	fprintf(stderr, "INFO: Maximum recommended exponent for this runlength = %llu; p[ = %llu]/pmax_rec = %12.10f.\n",given_N_get_maxP(n),p,exp_ratio);
	// Set initial value of USE_SHORT_CY_CHAIN based on how close p/pmax is to 1.0, but only if current chain length is longer
	// (e.g. if ROE-retry logic has led to a shorter-than-default chain length, don't revert to default):
	if(exp_ratio > 0.99 && USE_SHORT_CY_CHAIN < 2)
		USE_SHORT_CY_CHAIN = 2;
	else if(exp_ratio > 0.98 && USE_SHORT_CY_CHAIN < 1)
		USE_SHORT_CY_CHAIN = 1;
	const char*arr_sml[] = {"long","medium","short"};
	fprintf(stderr,"Initial DWT-multipliers chain length = [%s] in carry step.\n",arr_sml[USE_SHORT_CY_CHAIN]);

	if(kblocks < (i = get_default_fft_length(p)))
	{
		/* If it's at least close, allow it but print a warning; otherwise error out: */
		if((i << 10) == get_nextlarger_fft_length(n) && exp_ratio < 1.05)
		{
			fprintf(stderr, "specified FFT length %d K is less than recommended %d K for this p.\n",kblocks,i);
		}
		else
		{
			sprintf(cbuf, "specified FFT length %d K is much too small: Recommended length for this p = %d K.\n",kblocks,i);
			ASSERT(HERE, 0, cbuf);
		}
	}

	/*...If array padding turned on, check that the blocklength divides the unpadded runlength...	*/
	if((DAT_BITS < 31) && ((n >> DAT_BITS) << DAT_BITS) != n)
	{
		ASSERT(HERE, 0,"ERROR: blocklength does not divide runlength!");
	}

	/*...Find padded array length...	*/
	npad = n + ( (n >> DAT_BITS) << PAD_BITS );	/* length of padded data array.	*/
	/* If the residue and other modulus-size-dependent data arrays too small for the new assignment, deallocate them: */
	if(nalloc > 0 && npad > nalloc)
	{
		ASSERT(HERE, a_ptmp != 0x0 && a != 0x0 && b != 0x0 && c != 0x0 && d != 0x0,"Require (a_ptmp,a,b,c,d) != 0x0");
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
		}	// Must preserve i through alloc of arrtmp[] below
		j = 0;
		if(npad & 7)
			j = 8 - (npad & 7);
		nalloc = npad + j;	ASSERT(HERE, (nalloc & 7) == 0,"nalloc must be a multiple of 8!");	// This is so b,c,d enjoy same 64-byte alignment as a[]
		nbytes = nalloc<<3;
		ASSERT(HERE, a_ptmp == 0x0 && a == 0x0 && b == 0x0 && c == 0x0 && d == 0x0 && e == 0x0 && arrtmp == 0x0,"Require (a_ptmp,b,c,d,e,arrtmp) == 0x0");
		a_ptmp = ALLOC_DOUBLE(a_ptmp, 5*nalloc);	if(!a_ptmp){ sprintf(cbuf, "FATAL: unable to allocate array A in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		a      = ALIGN_DOUBLE(a_ptmp);
		ASSERT(HERE, ((long)a & 63) == 0x0,"a[] not aligned on 64-byte boundary!");
		if(((long)a & 127) != 0x0)
			fprintf(stderr, "WARN: a[] = 0x%08lX not aligned on 128-byte boundary!\n", (long)a);
		// v19: Add three more full-residue arrays to support 2-input FFT-modmul needed for Gerbicz check (and later, p-1 support):
		b = a + nalloc;	c = b + nalloc;	d = c + nalloc, e = d + nalloc;
		b_uint64_ptr = (uint64*)b; c_uint64_ptr = (uint64*)c; d_uint64_ptr = (uint64*)d; e_uint64_ptr = (uint64*)e;
		// For this residue (and scratch) byte-array, conservatively figure at least 4 bits per float-double residue word.
		// For multi-FFT-length self-tests, conservatively figure as many as 20 bits (2.5 bytes) per float-double residue word:
		i = MAX(p>>2, i*2.5);
		arrtmp = ALLOC_UINT64(arrtmp, i>>3);if(!arrtmp ){ sprintf(cbuf, "FATAL: unable to allocate array ARRTMP  with %u bytes in main.\n",i); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		// For an n-word main-array, BIGWORD_BITMAP and BIGWORD_NBITS have (n/64) elts each, thus need 1/64 + 1/32 the total
		// storage of the main-array. Use uint64 alloc-macro for both, so halve the num-elts arg for the BIGWORD_NBITS alloc.
		// As with above arrays, for multi-length self-test, alloc based on max. FFT length used (i) rather than current length (n).
		// Don't need any array padding on these bitmap arrays, but since nalloc includes padding, no harm in using it:
		BIGWORD_BITMAP =           ALLOC_UINT64(BIGWORD_BITMAP, nalloc>>6);	if(!BIGWORD_BITMAP){ sprintf(cbuf, "FATAL: unable to allocate array BIGWORD_BITMAP in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		BIGWORD_NBITS  = (uint32 *)ALLOC_UINT64(BIGWORD_NBITS , nalloc>>7);	if(!BIGWORD_NBITS ){ sprintf(cbuf, "FATAL: unable to allocate array BIGWORD_NBITS in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
	}

// Multithreaded-code debug: Set address to watch:
#ifdef MULTITHREAD
	ADDR0 = a;
#endif

	/* Make sure we start with primary restart file: */
	RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
	strcpy(cstr, RESTARTFILE);

READ_RESTART_FILE:

	if(!INTERACT)	/* Only do the restart-file stuff if it's not a self-test */
	{
		if(ierr == ERR_GERBICZ_CHECK)
			strcat(cstr, ".G");
		/* See if there's a restart file: */
		fp = mlucas_fopen(cstr, "rb");
		/* If so, read the savefile: */
		if(fp) {
			if(TEST_TYPE == TEST_TYPE_PRP) dum = PRP_BASE;
			i = read_ppm1_savefiles (p, &j, fp, &ilo, (uint8*)arrtmp, &Res64,&Res35m1,&Res36m1,	// Primality-test residue
												(uint8*)e_uint64_ptr, &i1   ,&i2     ,&i3     );// v19: G-check residue
			fclose(fp); fp = 0x0;
			if(!i) {
				fp = mlucas_fopen(   OFILE,"a");
				fq = mlucas_fopen(STATFILE,"a");
				/* First print any error message that may have been issued during the above function call: */
				if(strstr(cbuf, "read_ppm1_savefiles")) {
							fprintf(stderr,"%s",cbuf);
					if(fp){ fprintf(	fp,"%s",cbuf); }
					if(fq){ fprintf(	fq,"%s",cbuf); }
				}
				/* And now for the official spokesmessage: */
				snprintf(cbuf,STR_MAX_LEN, "ERROR: read_ppm1_savefiles Failed on savefile %s!\n",cstr);
						fprintf(stderr,"%s",cbuf);
				if(fp){ fprintf(	fp,"%s",cbuf); fclose(fp); fp = 0x0; }
				if(fq){ fprintf(	fq,"%s",cbuf); fclose(fq); fq = 0x0; }

				if(ierr == ERR_GERBICZ_CHECK)
					ASSERT(HERE, 0,"Failed to correctly read last-good-Gerbicz-check data savefile!");
				else if(cstr[0] != 'q') {
					cstr[0] = 'q';
					goto READ_RESTART_FILE;
				} else
					ASSERT(HERE, 0,"Failed to correctly read both primary or secondary savefile!");
			}
			// If user attempts to restart run with different PRP base than it was started with, ignore the new value and continue with the initial one:
			if(TEST_TYPE == TEST_TYPE_PRP && dum != PRP_BASE) {
				fprintf(stderr,"INFO: User-specified PRP-test base %u differs from value of %u read from savefile %s ... using the latter.\n",dum,PRP_BASE,cstr);
			}
			// If FFT-length-in-K field returns nonzero and is greater than kblocks (e.g. if ROEs caused switch to larger FFT length),
			// ***and user did not override via cmd-line*** (i.e. fft_length arg not set) reset FFT length and radix set accordingly:
			if((j && j > kblocks) && !fft_length) {
				kblocks = j;
				// Clear out current FFT-radix data, since get_preferred_fft_radix() expects that:
				for(i = 0; i < NRADICES; i++) { RADIX_VEC[i] = 0; }
				NRADICES = 0;
				goto SETUP_FFT;	// Do this for both ROE_ITER < 0 and > 0; in the latter case with unchanged FFT params
			}
			/* Allocate floating-point residue array and convert savefile bytewise residue to floating-point form, after
			first applying required circular shift read into the global RES_SHIFT during the above bytewise-savefile read.
			*/
			if(!convert_res_bytewise_FP((uint8*)arrtmp, a, n, p)) {
				snprintf(cbuf,STR_MAX_LEN, "ERROR: convert_res_bytewise_FP Failed on primality-test residue read from savefile %s!\n",cstr);
															fprintf(stderr,"%s",cbuf);
				fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				if(cstr[0] != 'q' && !(ierr == ERR_GERBICZ_CHECK)) {	// Secondary savefile only exists for regular checkpoint files
					cstr[0] = 'q';
					goto READ_RESTART_FILE;
				} else {
					ASSERT(HERE, 0,"0");
				}
			}
			// v19: G-check residue
		  if(TEST_TYPE == TEST_TYPE_PRP) {
			if(!convert_res_bytewise_FP((uint8*)e_uint64_ptr, b, n, p)) {
				snprintf(cbuf,STR_MAX_LEN, "ERROR: convert_res_bytewise_FP Failed on Gerbicz-check residue read from savefile %s!\n",cstr);
															fprintf(stderr,"%s",cbuf);
				fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
			} else {
				ierr = 0;
				s1 = sum64(b_uint64_ptr, n); s2 = s3 = s1;	// Init triply-redundant checksum of G-checkproduct
			}
		  }
			ASSERT(HERE, ilo > 0,"Require ilo > 0!");
			if(TEST_TYPE != TEST_TYPE_PRPCOFACTOR) {
				ASSERT(HERE, ilo < maxiter,"Require ilo < maxiter!");
				ihi = ilo+ITERS_BETWEEN_CHECKPOINTS;
				/* If for some reason last checkpoint was at a non-multiple of ITERS_BETWEEN_CHECKPOINTS, round down: */
				ihi-= ihi%ITERS_BETWEEN_CHECKPOINTS;
			} else {	// Init scratch array b[] needed for cofactor-PRP test:
INTERACT = 1;
				uint64 *ai = (uint64 *)a, *bi = (uint64 *)b, *atmp = (uint64 *)arrtmp;	// Handy 'precast' pointers
				/* Suyama cofactor-PRP test for Fermat/Mersenne modulus N, using base-3 PRP-test to illustrate:
					Let P = 3^[(N-1)/2] (mod N) be the base-3 Euler (Pepin-test) residue for N;
					Let F = product of known small-prime factors of N, C the corresponding cofactor, i.e. N = F*C;
					Let A = 3^(N-1) = P^2 (mod N) be the base-3 Fermat-PRP residue for N;
					Let B = 3^(F-1) (mod N).
				Then C is composite if R := (A - B) (mod C) = 0; otherwise if R = 0, C is a Fermat probable prime to base 3^F.
				Proof: A - B = 3^(F*C-1) - 3^(F-1) = 3^(F-1) * [ 3^(F*(C-1)-1) - 1 ] == 0 (mod N) and thus also == 0 (mod C)
				iff 3^(F*(C-1)-1) == 1 (mod C), which is just a base-3 Fermat-PRP test to base 3^F. QED
				*/
				fprintf(stderr,"Doing a cofactor-PRP test: maxiter = %u, savefile data are for iter = %u\n",maxiter,ilo);
				ilo = 0;	ihi = ilo+1;	// Have checked that savefile residue is for a complete PRP test, so reset iteration counter
				fprintf(stderr, "Doing one mod-%s squaring of iteration-%u residue [Res64 = %016llX] to get Fermat-PRP residue\n",PSTRING,ilo,Res64);
				fprintf(stderr, "%s: using FFT length %uK = %u 8-byte floats.\n",PSTRING,kblocks,n);
				fprintf(stderr, " this gives an average %20.15f bits per digit\n",1.0*p/n);
				// Compute Fermat-PRP residue [A] from Euler-PRP (= Pepin-test) residue via a single mod-squaring:
				BASE_MULTIPLIER_BITS[0] = 0ull;
		/*A*/	ierr = fermat_mod_square(a, (int*)arrtmp, n, ilo, ihi, p, scrnFlag, &tdiff);
				convert_res_FP_bytewise(a, (uint8*)arrtmp, n, p, &Res64, &Res35m1, &Res36m1);
				fprintf(stderr, "Fermat-PRP residue: Res64 = %llX\n",Res64);
				mi64_set_eq(ai,atmp,(p+63)>>6);	// Copy packed-bit result back into low ceiling(p/64) bytes of A-vec (treated as a uint64 array)
				// Compute "prime-factor product residue" [B] from Euler-PRP (= Pepin-test) residue ... first init bitwise mul-by-base array = F, i.e. storing product of known small-prime factors:
				ASSERT(HERE, nfac > 0,"Cofactor-PRP test requires one or more known factors!");
				BASE_MULTIPLIER_BITS[0] = 1ull;	lenf = 1;
				// Must recompute each factor from its k-value before multiplying with current partial product of factors:
				uint64 curr_fac[10];
				for(i = 0; i < nfac; i++) {
					curr_fac[0] = factor_k[i];	curr_fac[1] = 0ull;	j = 1;	// j = number of nonzero limbs in curr_fac
					j += mi64_mul_scalar(curr_fac, twop, curr_fac, j);	curr_fac[0] += 1ull;
					// Multiply factor into current partial product of factors; use b-array to store product to work around none-of-3-input-pointers-may-coincide restriction in mi64_mul_vector:
					mi64_mul_vector(BASE_MULTIPLIER_BITS,lenf, curr_fac,j, bi,&lenf);
					mi64_set_eq(BASE_MULTIPLIER_BITS,bi,lenf);
				}
				// Only use curr_fac to store one individual factor at a time here, but need it for all of F later on
				ASSERT(HERE, lenf <= 10, "Product of factors too large to fit into curr_fac[]!");
				for(i = 0; i < lenf; i++) { b[i] = 0.0; }	// Re-zero the elts of b used as tmps in above loop
				fbits = (lenf<<6) - mi64_leadz(BASE_MULTIPLIER_BITS, lenf);
				// Now that have F stored in BASE_MULTIPLIER_BITS array, do powmod to get B = base^(F-1) (mod N):
				BASE_MULTIPLIER_BITS[0] -= 1ull;	// F-1; no chance of a borrow here
				b[0] = PRP_BASE;	ASSERT(HERE, PRP_BASE < (1 << (uint32)ceil(1.0*p/n)), "PRP_BASE out of range!");
				ilo = 0;	ihi = fbits-1;	// LR modpow; init b[0] = PRP_BASE takes cares of leftmots bit
/************** For mers-mod support, do we need to replace p>>6 with (p+63)>>6 in mi64 calls below? ***************/
				mi64_brev(BASE_MULTIPLIER_BITS,ihi);	// bit-reverse low [ihi] bits of BASE_MULTIPLIER_BITS:
		/*B*/	ierr = fermat_mod_square(b, (int*)arrtmp, n, ilo, ihi, p, scrnFlag, &tdiff);
				convert_res_FP_bytewise(b, (uint8*)arrtmp, n, p, &Res64, &Res35m1, &Res36m1);	// Res64 = 0x25f5ab0ffc728c87
				fprintf(stderr, "Processed %u bits in binary modpow. %u^(F-1) residue: Res64 = %llX\n",ihi,PRP_BASE,Res64);
				mi64_set_eq(bi,atmp,p>>6);	// Copy packed-bit result back into low ceiling(p/8) bytes of A-vec (treated as a unit64 array)
				itmp64 = mi64_sub(ai,bi, ai,p>>6);
				// If result < 0, need to add Modulus - for N = Fm,Mp this means +-1 in LSW, respectively:
				if(itmp64) {	ASSERT(HERE, itmp64 == 1ull,"Carryout = 1 expected!");
					if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
						itmp64 = mi64_sub_scalar(ai,1ull, ai,p>>6);
					} else {
						itmp64 = mi64_add_scalar(ai,1ull, ai,p>>6);
					}	ASSERT(HERE, itmp64 == 0ull,"Carryout = 0 expected!");
				}
				// B-array again free, re-use in uint64-cast form to compute C = Fm/F and (A-B) mod C:
				// Compute Modulus ... note mi64-vecs have no cache-oriented element padding:
				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
					itmp64 =-1ull;	j = (p >> 6);	// j = uint64 vector length; init sans the leading word, which needs special handling
					mi64_set_eq_scalar(bi,itmp64,j);
					itmp64 >>= (p&63);	b[++j] = *((double *)&itmp64);	// Pre-increment here
				} else {
					itmp64 = 0ull;	j = (p >> 6);	// j = uint64 vector length; init sans the leading '1' word, then increment prior to mi64_div
					mi64_set_eq_scalar(bi,itmp64,j);
					itmp64 = 1ull;	b[j++] = b[0] = *((double *)&itmp64);	// Post-increment here
				}
		// C = N/F:
				mi64_brev(BASE_MULTIPLIER_BITS,ihi);// 2nd BR of low [ihi] bits of BASE_MULTIPLIER_BITS to recover the factored part F-1, sans leftmost bit...
				BASE_MULTIPLIER_BITS[lenf-1] += 1ull << (fbits-1);	// Restore leftmost bit ...
				BASE_MULTIPLIER_BITS[     0] += 1ull;	// ... and add 1 to recover F; no chance of a carryout here
				// Since F << N, use Mont-mul-div for C - quotient overwrites N, no rem-vec needed, just verify that F is in fact a divisor:
				ASSERT(HERE, 1 == mi64_div(bi,BASE_MULTIPLIER_BITS, j,lenf, atmp,0x0), "C = N/F should have 0 remainder!");
				i = (p+63)>>6;	j = mi64_getlen(atmp, j);
				// R = (A - B) mod C in B-array; store Q = (A - B)/C in curr_fac[] in case want to remultiply and verify Q*C + R = (A - B):
				mi64_div_binary(ai,atmp, i,j, curr_fac,(uint32 *)&j, bi);	// On return, j has quotient length
				// For 1-word quotient q, double-check binary-div result by computing (q*denominator + r) and comparing vs numerator:
				if(j == 1) {
					ASSERT(HERE, 0 == mi64_mul_scalar_add_vec2(atmp, curr_fac[0], bi, atmp, i), "Unexpected carryout!");
					ASSERT(HERE, 1 == mi64_cmp_eq(ai,atmp,i), "Q*C + R = (A - B) check fails!");
				}
				/* Compute S-H residues of B-array here or just use Res64 = ((uint64*)b)[0]? */
				printf("Suyama Cofactor-PRP test of %s",PSTRING);
				// Base-2 log of cofactor = lg(Fm/F) = lg(Fm) - lg(F) ~= 2^m - lg(F). 2^m stored in p, sub lg(F) in loop below:
				double lg_cof = p,lg_fac,log10_2 = 0.30102999566398119521;	// Use lg_fac to store log2 of each factor as we recompute it
				for(i = 0; i < nfac; i++) {
					curr_fac[0] = factor_k[i];	curr_fac[1] = 0ull;	j = 1;	// j = number of nonzero limbs in curr_fac
					j += mi64_mul_scalar(curr_fac, twop, curr_fac, j);	curr_fac[0] += 1ull;
					printf(" / %s",&cbuf[convert_mi64_base10_char(cbuf, curr_fac, j, 0)] );
					lg_fac  = (double)mi64_extract_lead64(curr_fac, j, &itmp64) - 64;
					lg_fac += log((double)itmp64)*ILG2;
					lg_cof -= lg_fac;
				}
				i = (p+63)>>6;	j = mi64_getlen(bi, i);	// Returns 0 iff all limbs of remainder == 0
				printf(" with FFT length %u = %u K: Res64: %16llX.\n",n,kblocks,bi[0]);
				i = ceil(lg_cof*log10_2);
				if(!j)
					printf("This cofactor is PROBABLE PRIME [PRP%u].\n",i);
				else
					printf("This cofactor is COMPOSITE [C%u].\n",i);
exit(0);
			}
			restart = TRUE;
		}
		else /* if(!fp) */
		{
			/* If we're on the primary restart file, set up for secondary: */
			if(cstr[0] != 'q' && !(ierr == ERR_GERBICZ_CHECK)) {	// Secondary savefile only exists for regular checkpoint files
				snprintf(cbuf,STR_MAX_LEN, "INFO: primary restart file %s not found...looking for secondary...\n",cstr);
													fprintf(stderr,"%s",cbuf);
				fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				cstr[0] = 'q';
				goto READ_RESTART_FILE;
			} else {
				sprintf(cbuf, "INFO: no restart file found...starting run from scratch.\n");	fprintf(stderr,"%s",cbuf);
				fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				if(ierr == ERR_GERBICZ_CHECK) { ierr = 0; restart = FALSE; }
			}
		}	/* endif(fp) */
	}	/* endif(!INTERACT)	*/

	if(!restart) {
		/*...set initial iteration loop parameters...	*/
		ilo=0;
		if(INTERACT)
			ihi=timing_test_iters;
		else
			ihi=ITERS_BETWEEN_CHECKPOINTS;
	}

	// PRP-test: Init bitwise multiply-by-base array - cf. comment re. modified Fermat-PRP needed by Gerbicz check
	// ~line 1070 above ==> all bits = 0 for Mersenne-PRP-test, not all-1s-with-LS-bit-0 as for the unmodfied Fermat-PRP
	// test of a Mersenne. Thus also no need to call mi64_brev to put bitmap in BRed form used by LR modpow:
	if(TEST_TYPE == TEST_TYPE_PRP) {
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			itmp64 = 0ull;//0xFFFFFFFFFFFFFFFFull;
		else
			itmp64 = 0ull;
		for(i = 0; i < ((ITERS_BETWEEN_CHECKPOINTS+63) >> 6); i++) {
			BASE_MULTIPLIER_BITS[i] = itmp64;
		}
	}
	ASSERT(HERE,MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(HERE,MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	/* If at the start of a p-1 or primality test, set the initial seed for the run: */
	ASSERT(HERE, TEST_TYPE <= TEST_TYPE_MAX,"Given TEST_TYPE not supported!");
	if(ilo == 0)
	{
		memset(a, 0, npad*sizeof(double));
		memset(b, 0, npad*sizeof(double));
		memset(c, 0, npad*sizeof(double));
		memset(d, 0, npad*sizeof(double));

		/* Always use 3 as the p-1 and Pepin-test seed, and 4 for the LL-test seed. For PRP-test, use seed set in worktodo assignment line: */
		if(TEST_TYPE == TEST_TYPE_PRP || TEST_TYPE == TEST_TYPE_PRPCOFACTOR) {
			iseed = b[0] = d[0] = PRP_BASE;	// If doing a PRP test, init the Gerbicz residue-product accumulator b[] and its redundant copy d[]. On restart b[] inited via full-bytewise-array read from savefile.
			s1 = s2 = s3 = PRP_BASE;	// Init triply-redundant checksum of G-checkproduct
		} else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
			iseed = 4;
		else	// Pepin test:
			iseed = 3;

	#ifndef USE_SSE2	// Only support residue-shift for SIMD builds:
		RES_SHIFT = 0ull; a[0] = iseed;
	#else
		if(RES_SHIFT && MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
			fprintf(stderr,"Shifted residues currently unsupported for Fermat-mod! Setting shift = 0 for run.\n");
			RES_SHIFT = 0ull; a[0] = iseed;
		} else {
			// Apply initial-residue shift - if user has not set one via cmd-line or current value >= p, randomly choose a value in [0,p).
			// [Note that the RNG is inited as part of the standard program-start-sequence, via function host_init().]
			if(RES_SHIFT == -1ull || RES_SHIFT >= p)	// Since 0 a valid value of RES_SHIFT, treat -1 as "uninited"
				RES_SHIFT = rng_isaac_rand() % p;
			else if(RES_SHIFT == 0) {	// This is needed for multi-exponent runs in which the first assignment is a partially-complete v17 run:
				itmp64 = parse_cmd_args_get_shift_value();
				if(itmp64 == -1ull)	// -1 indicates "-shift not specified on command line", in which case we init-to-random since this is a v18 run-from-scratch
					RES_SHIFT = rng_isaac_rand() % p;
				else	// Otherwise set shift to whatever value was specified on command line
					RES_SHIFT = itmp64 % p;
			}
			// Since residue is otherwise 0, use shifted-carryin function on double-precision padded-array residue:
			itmp64 = shift_word(a, n, p, RES_SHIFT, (double)iseed);	// Note return value (specifically high 7 bytes thereof) is an unpadded index
			ASSERT(HERE, (itmp64 >>  8) < n                , "Return value of shift_word(): unpadded-array-index out of range!");
			ASSERT(HERE, (itmp64 & 255) < ceil((double)p/n), "Return value of shift_word(): bit-in-array-word value out of range!");
		}
	#endif
	} else if(TEST_TYPE == TEST_TYPE_PRP || TEST_TYPE == TEST_TYPE_PRPCOFACTOR)
		memcpy(d, b, nbytes);	// If doing a PRP test, init redundant copy d[] Gerbicz residue-product accumulator b[].

	if(restart)
	{
		snprintf(cbuf,STR_MAX_LEN, "Restarting %s at iteration = %u. Res64: %016llX, residue shift count = %llu\n",PSTRING,ilo,Res64,RES_SHIFT);
	}

	/*...Restart and FFT info.	*/
	if(INTERACT)
	{
		fprintf(stderr, "%s: using FFT length %uK = %u 8-byte floats, initial residue shift count = %llu\n",PSTRING,kblocks,n,RES_SHIFT);
		fprintf(stderr, " this gives an average %20.15f bits per digit\n",1.0*p/n);
		if(TEST_TYPE == TEST_TYPE_PRP || TEST_TYPE == TEST_TYPE_PRPCOFACTOR)
			fprintf(stderr, "The test will be done in form of a %u-PRP test.\n",PRP_BASE);
	}
	else
	{
		fp = mlucas_fopen(STATFILE,"a");
		if(fp)
		{
			fprintf(fp,"%s",cbuf);
			fprintf(fp,"%s: using FFT length %uK = %u 8-byte floats, initial residue shift count = %llu\n",PSTRING,kblocks,n,RES_SHIFT);
			fprintf(fp," this gives an average %20.15f bits per digit\n",1.0*p/n);
			if(TEST_TYPE == TEST_TYPE_PRP || TEST_TYPE == TEST_TYPE_PRPCOFACTOR)
				fprintf(fp, "The test will be done in form of a %u-PRP test.\n",PRP_BASE);
			fclose(fp);	fp = 0x0;
		}
	}

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

	for(;;)
	{
		ASSERT(HERE, (int)maxiter > 0,"Require (int)maxiter > 0");
		if(ihi > maxiter)
			ihi = maxiter;
		/* Here's the big one - (ITERS_BETWEEN_CHECKPOINTS) squaring steps.
		XYZ_mod_square returns 0 if no errors detected during this iteration cycle.

		If fatal error was encountered, skip to next assignment in worktodo.ini file
		but keep restart files for current assignment around (can finish using hiacc code.)
		*/
		AME = MME = 0.0;	/* Init Avg. & Max. RO Error */
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_FGT61
			fwd_fft_only = 0ull;
			ierr = mers_mod_square  (a,c, (int*)arrtmp, n, ilo, ihi, fwd_fft_only, p, scrnFlag, &tdiff, TRUE);
		#elif USE_2STEP	/****************** test out 2-input modmul ******************/
		  fwd_fft_only = (uint64)b;	// Copy a into b and do fwd-FFT-only of b:
		  uint32 iii;
		  for(iii = ilo; iii < ihi; iii++) {
			memcpy(b, a, nbytes);	// Copy a into b and do fwd-FFT-only of b:
			ierr = mers_mod_square  (b, (int*)arrtmp, n, iii,iii+1, 1ull        , p, scrnFlag, &tdiff, TRUE);
			ierr = mers_mod_square  (a, (int*)arrtmp, n, iii,iii+1, fwd_fft_only, p, scrnFlag, &tdiff, TRUE);
		  }
		#else
		  // v19 Gerbicz-checkproduct updates: break iteration loop into subsegments of length ITERS_BETWEEN_GCHECK_UPDATES.
		  if(TEST_TYPE == TEST_TYPE_PRP) {
			// Only do check-related stuff for full production-run iteration intervals:
			if(MLUCAS_KEEP_RUNNING && (ihi-ilo) >= ITERS_BETWEEN_GCHECK_UPDATES) {
				i = ilo;	tdiff = 0.0;	// Need 2 timers here - tdif2 for the individual mers_mod_square calls, accumulate in tdiff
				while(!ierr && MLUCAS_KEEP_RUNNING && i < ihi) {
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
					uint32 first_sub = (i == ilo), last_sub = (i+itodo == ihi);
				// First subinterval: [a] needs fwd-weighting and initial-fwd-FFT-pass done on entry, !undone on exit: mode_flag = 10_2
				// Last  subinterval: [a] !need fwd-weighting and initial-fwd-FFT-pass done on entry but done on exit: mode_flag = 01_2
				// Intermediate subs: [a] !need fwd-weighting and initial-fwd-FFT-pass done on entry, !undone on exit: mode_flag = 11_2
					mode_flag = 3 - first_sub - (last_sub<<1);
					ierr = mers_mod_square  (a, (int*)arrtmp, n, i,i+itodo, (uint64)mode_flag, p, scrnFlag, &tdif2, TRUE);	tdiff += tdif2;
					if(ierr) {
						fprintf(stderr,"At iteration %u: mers_mod_square returned with err code %u\n",ROE_ITER,ierr);
						/* If interrupt *and* we're past the first subinterval, need to undo initial-fwd-FFT-pass and DWT-weighting on b[],
						whose value will reflect the last multiple-of-ITERS_BETWEEN_GCHECK_UPDATES iteration - prior to writing it,
						along with the current PRP residue, to savefile: */
						if(ierr == ERR_INTERRUPT && !first_sub)
							ierr = mers_mod_square  (b, (int*)arrtmp, n, i,i+1, 8ull, p, scrnFlag, &tdif2, FALSE);
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
				would prevent the ensuing mers_mod_square call to compute FFT(b)*FFT(c) from doing so ... so instead break out of loop: */
					mode_flag = 1 - last_sub;
					ierr = mers_mod_square  (c, (int*)arrtmp, n, i,i+1,      4ull + (uint64)mode_flag, p, scrnFlag, &tdif2, FALSE);
					if(ierr) {
						if(ierr == ERR_INTERRUPT) {
							fprintf(stderr,"Caught interrupt in fFFT(c) step.\n");
							break;
						} else {
							returnMlucasErrCode(ierr,cstr);
							sprintf(cbuf,"Unhandled error of type %s in fFFT(c) step - please send e-mail to ewmayer@aol.com with copy of the p*.stat file attached. Proceeding to next assignment...\n",cstr);
							fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
							fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
							ASSERT(HERE,0,cbuf);
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
						ASSERT(HERE, 0, "Catastrophic data corruption detected in G-checkproduct integrity validation ... rolling back to last good G-check. ");
				}

				// First subinterval: [b] needs fwd-weighting and initial-fwd-FFT-pass done on entry, !undone on exit: mode_flag = 10_2
				// Intermediate subs: [b] !need fwd-weighting and initial-fwd-FFT-pass done on entry, !undone on exit: mode_flag = 11_2
				// Last  subinterval: [b] !need fwd-weighting and initial-fwd-FFT-pass done on entry,  undone on exit: mode_flag = 01_2
				/* Note: Interrupt during this step should not be a problem, the handling code in mers_mod_square will complete the FFT-mul
				step and force the undo-initial-FFT-pass-and-DWT-weighting step, leaving a pure-int G-check residue ready for savefile-writing: */
					mode_flag = 3 - first_sub - (last_sub<<1);
				//	printf("Iter %u: FFT(b)*FFT(c) step.\n",i);
					ierr = mers_mod_square  (b, (int*)arrtmp, n, i,i+1, (uint64)c + (uint64)mode_flag, p, scrnFlag, &tdif2, FALSE);
					if(ierr) {
						if(ierr == ERR_INTERRUPT) {
							fprintf(stderr,"Caught interrupt in FFT(b)*FFT(c) step.\n");
							break;
						} else {
							returnMlucasErrCode(ierr,cstr);
							sprintf(cbuf,"Unhandled error of type %s in FFT(b)*FFT(c) step - please send e-mail to ewmayer@aol.com with copy of the p*.stat file attached. Proceeding to next assignment...\n",cstr);
							fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
							fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
							ASSERT(HERE,0,cbuf);
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
				}
			} else if(MLUCAS_KEEP_RUNNING) {	// Final partial-length interval skips G-check
				ierr = mers_mod_square  (a, (int*)arrtmp, n, ilo,ihi, 0ull, p, scrnFlag, &tdiff, TRUE);
			}
		  } else {
				// For straight LL-test there is (at least at this writing) no known analog of the Gerbicz check:
				ierr = mers_mod_square  (a, (int*)arrtmp, n, ilo,ihi, 0ull, p, scrnFlag, &tdiff, TRUE);
		  }
		#endif
		}
		else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
		#if 0//def USE_FGT61	**** First get things working for LL case ****
			ierr = fermat_mod_square(a,c, (int*)arrtmp, n, ilo, ihi, p, scrnFlag, &tdiff);
		#else
			ierr = fermat_mod_square(a,   (int*)arrtmp, n, ilo, ihi, p, scrnFlag, &tdiff);
		#endif
		}

		// Roundoff-retry scheme is detailed in comments of the above fermat_mod_square() function:
		if(ierr) {
			// v19: For "nonzero exit carry" errors due to residue data corruption (not uncommon on e.g. cellphones)
			// added handling - if at least one previous savefile write, try restart from that:
			if((ierr == ERR_CARRY) && !INTERACT && ihi > ITERS_BETWEEN_CHECKPOINTS) {
				sprintf(cbuf,"Retrying from most-recent savefile, in case the issue was due to residue-data corruption.\n");
				fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				goto READ_RESTART_FILE;
			}
			if((ierr == ERR_ROUNDOFF) && !INTERACT) {
				if(!ROE_ITER) {
					ASSERT(HERE, 0, "ERR_ROUNDOFF returned but ROE_ITER not set!");
				} else if(ROE_ITER > 0) {
					// [1] [condition] (ROE_ITER > 0) :
					//     [action] Rerun the interval with the same FFT params (via goto below)
				} else if(ROE_ITER < 0) {
					// [2b] [condition] (ROE_ITER < 0) The error is reproducible, i.e. recurs on same iter, with same fracmax.
					//     [action] Rerun interval needs (and the rest of the run) with the next-larger available FFT length:
					if(ROE_VAL == 0.0) {	// Case [2b]
						// v19: If already tried reducing DWT-multipliers chain length in carry step,
						// revert on switch to larger FFT length; otherwise try reduced chain length:
						fp = mlucas_fopen(   OFILE,"a");
						fq = mlucas_fopen(STATFILE,"a");
						if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {
							sprintf(cbuf," Reducing DWT-multipliers chain length from [%s] to [%s] in carry step to see if this prevents further danger-level fractional parts.\n",arr_sml[USE_SHORT_CY_CHAIN],arr_sml[USE_SHORT_CY_CHAIN+1]);
							USE_SHORT_CY_CHAIN++;
						} else {
							sprintf(cbuf," Switching to next-larger available FFT length %uK and retrying.\n",(n >> 10));
							n = get_nextlarger_fft_length(n);
							USE_SHORT_CY_CHAIN = 0;
						}
						fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
						fprintf(fq,"%s",cbuf);	fclose(fq);	fq = 0x0;
					} else {				// Case [2c]
						kblocks = get_default_fft_length(p);	// Default FFT length in Kdoubles for this exponent
						if(n > (kblocks << 10))		// We are already running at a larger-than-default FFT length
							n = kblocks << 10;
						else {
							n = get_nextlarger_fft_length(n);
							// v19: If already tried reducing DWT-multipliers chain length in carry step, revert on switch to larger FFT length:
							USE_SHORT_CY_CHAIN = 0;
						}
					}
					ROE_ITER = 0;	kblocks = (n >> 10);
				}
				ierr = 0;	// v19: Need to explicitly clear ierr flag here, otherwise get oo retry loop in PRP-test mode
				// Clear out current FFT-radix data, since get_preferred_fft_radix() expects that:
				for(i = 0; i < NRADICES; i++) { RADIX_VEC[i] = 0; }
				NRADICES = 0;
				goto SETUP_FFT;	// Do this for both ROE_ITER < 0 and > 0; in the latter case with unchanged FFT params
			}
			else if(ierr == ERR_UNKNOWN_FATAL)
				return(EXIT_FAILURE);
			if(ierr != ERR_INTERRUPT)	// Unknown error ... try simply going to next assignment
				goto GET_NEXT_ASSIGNMENT;
		}

	/*...Every (ITERS_BETWEEN_CHECKPOINTS)th iteration, print timings to stdout or STATFILE.
	If it's a regular (i.e. non-timing) test, also write current residue to restart files.
	*/
	//	if(RES_SHIFT) for(j = 0; j < n; j++) a[j] = -a[j];	// *************** DEBUG (attempted) of Fermat-mod-with-shift *******************

		/* v18: Moved term-output below convert_res_FP_bytewise call since in the presence of cshifted residues
		need convert_res_FP_bytewise call to remove cshift and compute SH residues: */
		// Zero high uint64s of target arrays, since double-to-int residue conversion is bytewise & may leave >=1 MSBs in high word untouched:
		j = (p+63)>>6;	arrtmp[j-1] = e_uint64_ptr[j-1] = 0ull;
		convert_res_FP_bytewise(	a, (uint8*)      arrtmp, n, p, &Res64, &Res35m1, &Res36m1);	// LL/PRP-test residue
		if(TEST_TYPE == TEST_TYPE_PRP) {
			convert_res_FP_bytewise(b, (uint8*)e_uint64_ptr, n, p, &i1,&i2,&i3);	// G-check residue...must not touch i1,i2,i3 again until ensuing write_ppm1_savefiles call!
		}
		// In interactive-timing-test (e.g. self-tests) mode, do immediate-exit-sans-savefile-write on signal:
		if(INTERACT && (ierr == ERR_INTERRUPT))
			exit(0);

		// In non-interactive (production-run) mode, write savefiles and exit gracefully on signal:
		if(ierr == ERR_INTERRUPT) {
			fp = mlucas_fopen(STATFILE,"a");
			if(fp){ fprintf(fp,"%s",cbuf); }	// First print the signal-handler-generated message
			ihi = ROE_ITER;	// Last-iteration-completed-before-interrupt saved here
			fprintf(stderr,"Iter = %u: Writing savefiles and exiting.\n",ihi);
			sprintf(cbuf  ,"Iter = %u: Writing savefiles and exiting.\n",ihi);
			if(fp){	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;	}
		}

		/*...Done?	*/
		if(!INTERACT) {
			AME /= (ihi - ilo);	// Don't /= ITERS_BETWEEN_CHECKPOINTS here since final interval is a partial one
			/*...get a quick timestamp...	*/
			calendar_time = time(NULL);
			local_time = localtime(&calendar_time);
			strftime(timebuffer,SIZE,"%Y-%m-%d %H:%M:%S",local_time);

			/*...print runtime in hh:mm:ss format.	*/
			snprintf(cbuf,STR_MAX_LEN, "[%s] %s Iter# = %u [%5.2f%% complete] clocks =%s [%8.4f msec/iter] Res64: %016llX. AvgMaxErr = %10.9f. MaxErr = %10.9f. Residue shift count = %llu.\n"
				, timebuffer, PSTRING, ihi, (float)ihi / (float)p * 100,get_time_str(tdiff)
				, 1000*get_time(tdiff)/(ihi - ilo), Res64, AME, MME, RES_SHIFT);

			fp = mlucas_fopen(STATFILE,"a");	  if(fp){	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;	}
			if (scrnFlag)
				fprintf(stderr,"%s",cbuf);	/* Echo output to stddev */
		}

		// Do not save a final residue unless cofactor-PRP mode (but this leaves the penultimate residue file intact):
		if( (ihi == maxiter) && (TEST_TYPE != TEST_TYPE_PRPCOFACTOR)) break;

		/* v19: If PRP test, do Gerbicz-check every million squarings. Must make sure current run has done at least one
		intermediate G-checkproduct update. How do we know that current iter-interval (ilo,ihi] contains such an update?
		(That is, a multiple of ITERS_BETWEEN_GCHECK_UPDATES): Test if (ihi % ITERS_BETWEEN_GCHECK_UPDATES > ilo):
		Examples, all using ITERS_BETWEEN_GCHECK_UPDATES = 1000:
		o ilo = 0, ihi = 1000: integer divide gives ilo/1000 = 0, ihi/1000 = 1 difference > 0, thus contains an update
		o ilo = 1000, ihi = 1507: ilo/1000 = 1, ihi/1000 = 1 difference = 0, thus !contains an update
		*/
		i = ITERS_BETWEEN_GCHECKS; j = ITERS_BETWEEN_GCHECK_UPDATES;
		if(MLUCAS_KEEP_RUNNING && TEST_TYPE == TEST_TYPE_PRP && (ilo/j < ihi/j) && (ihi % i) == 0)
		{
			// Un-updated copy of the checkproduct saved in [d]; square that ITERS_BETWEEN_GCHECK_UPDATES times...
			// [d] !need fwd-weighting and initial-fwd-FFT-pass done on entry but undone on exit: mode_flag = 01_2:
			mode_flag = 1;
			ierr = mers_mod_square  (d,0x0, n, ihi,ihi+ITERS_BETWEEN_GCHECK_UPDATES, (uint64)mode_flag, p, scrnFlag, &tdiff, FALSE);
			if(ierr) {
				if(ierr == ERR_INTERRUPT) {
					fprintf(stderr,"Caught interrupt in Gerbicz-checkproduct mod-squaring update ... skipping G-check and savefile-update and performing immediate-exit.\n");
					exit(0);
				} else {
					returnMlucasErrCode(ierr,cstr);
					sprintf(cbuf,"Unhandled error of type %s in Gerbicz-checkproduct mod-squaring update - please send e-mail to ewmayer@aol.com with copy of the p*.stat file attached. Proceeding to next assignment...\n",cstr);
					fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					ASSERT(HERE,0,cbuf);
				//	goto GET_NEXT_ASSIGNMENT;
				}
			}
			// ...then multiply squaring result by initial PRP seed and compare result to updated copy of the checkproduct saved in b[].
			// First zero the high uint64s of the target array, since the double-to-int residue conversion is bytewise, i.e. may leave 1 or more MSBs in high word untouched:
			j = (p+63)>>6;	c_uint64_ptr[j-1] = 0ull;
			// [1] Convert b[],d[] to bytewise form, former assumed already in e[] doubles-array, latter into currently-unused c[] doubles-array:
			convert_res_FP_bytewise(d, (uint8*)c_uint64_ptr, n, p, 0x0,0x0,0x0);
			// Only need to compute this for initial interval - after that the needed adjustment-shift remains constant
			if(ihi == ITERS_BETWEEN_GCHECKS && RES_SHIFT) {
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
				fprintf(stderr,"Recovered initial shift %llu\n",itmp64);
				ASSERT(HERE, (itmp64>>32) == 0ull,"Shift must be < 2^32!");
				GCHECK_SHIFT = itmp64;
			}
			mi64_shlc(c_uint64_ptr, c_uint64_ptr, (uint32)p, (uint32)GCHECK_SHIFT, j);
			// Use mi64 routines to compute d[]*PRP_BASE and do ensuing equality check:
			c_uint64_ptr[j] = mi64_mul_scalar(c_uint64_ptr, (uint64)PRP_BASE, c_uint64_ptr, j);
		//	ASSERT(HERE, c_uint64_ptr[j] == 0ull, "d[]*PRP_BASE result has unexpected carryout!");
			// Need to (mod 2^p-1) ... store 2^p-1 in d[] doubles-array, which is freed up by above convert_res_FP_bytewise(d,...) call:
			d_uint64_ptr[j-1] = (1ull<<(p&63)) - 1ull;
			for(i = 0; i < (j-1); i++) { d_uint64_ptr[i] = -1ull; }
			for(i = 1; i < PRP_BASE; i++) {
				if(!c_uint64_ptr[j] && mi64_cmpult(c_uint64_ptr,d_uint64_ptr,j)) break;
				itmp64 = mi64_sub(c_uint64_ptr,d_uint64_ptr,c_uint64_ptr,j);	// c -= d, with d = 2^p-1
				c_uint64_ptr[j] -= itmp64;	//ASSERT(HERE, itmp64 == 0ull, "mi64_sub result has unexpected borrow!");
			}
			ASSERT(HERE, mi64_cmpult(c_uint64_ptr,d_uint64_ptr,j), "Gerbicz checkproduct reduction (mod 2^p-1) failed!");
			if(mi64_cmp_eq(e_uint64_ptr,c_uint64_ptr,j)) {
				sprintf(cbuf,"At iteration %u, shift = %llu: Gerbicz check passed.\n",ihi,RES_SHIFT);
				fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				// In G-check case we need b[] for that, thus skipped the d = b redundancy-copy ... do that now:
				memcpy(d, b, nbytes);
				s1 = sum64(b_uint64_ptr, n); s2 = s3 = s1;	// Init triply-redundant checksum of G-checkproduct
			} else {
				i = mi64_shlc_bits_align(e_uint64_ptr,c_uint64_ptr,p);
				if(i != -1) {
					sprintf(cbuf,"Gerbicz check passes if D *= 2^%u (mod 2^p-1)\n",i);
					fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					// In G-check case we need b[] for that, thus skipped the d = b redundancy-copy ... do that now:
					memcpy(d, b, nbytes);
					s1 = sum64(b_uint64_ptr, n); s2 = s3 = s1;	// Init triply-redundant checksum of G-checkproduct
				} else {
					sprintf(cbuf,"Gerbicz check failed! Restarting from last-good-Gerbicz-check data, or from scratch if iteration < %u\n",ITERS_BETWEEN_GCHECKS);
					fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					ierr = ERR_GERBICZ_CHECK;
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
				snprintf(cbuf,STR_MAX_LEN,"ERROR: unable to rename %s restart file ==> %s ... skipping every-10M-iteration restart file archiving\n",RANGEFILE,cstr);
				fprintf(stderr,"%s",cbuf);
			}
		}	// ilo a multiple of 10 million?

	WRITE_RESTART_FILE:

		fp = mlucas_fopen(RESTARTFILE, "wb");
		if(fp) {
			write_ppm1_savefiles(p,n,fp, ihi, (uint8*)arrtmp,Res64,Res35m1,Res36m1, (uint8*)e_uint64_ptr,i1,i2,i3);
			fclose(fp);	fp = 0x0;
			/* If we're on the primary restart file, set up for secondary: */
			if(RESTARTFILE[0] != 'q') {
				RESTARTFILE[0] = 'q';
				goto WRITE_RESTART_FILE;
			} else {
				RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
			}
		} else {
			snprintf(cbuf,STR_MAX_LEN, "ERROR: unable to open restart file %s for write of checkpoint data.\n",RESTARTFILE);
										        fprintf(stderr,"%s",cbuf);
			fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
			fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
			/*
			Don't want to assert here - asllow processing to continue, in case this is a transient failure-to-open.
			e.g. due to a backup utility having temporarily locked the savefile, or an out-of-disk-space problem.
			*/
		}

		// v19: For PRP tests, add every-ITERS_BETWEEN_GCHECKS-iter last-good-Gerbicz-check savefile,
		// whose content gets updated that frequently but whose name remains fixed:
		RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
		if(TEST_TYPE == TEST_TYPE_PRP) {
			if(ihi%ITERS_BETWEEN_GCHECKS == 0) {
				strcpy(cstr, RESTARTFILE);
				strcat(cstr, ".G");
				fp = mlucas_fopen(cstr, "wb");
				if(fp) {
					write_ppm1_savefiles(p,n,fp, ihi, (uint8*)arrtmp,Res64,Res35m1,Res36m1, (uint8*)e_uint64_ptr,i1,i2,i3);
					fclose(fp);	fp = 0x0;
				} else {
					snprintf(cbuf,STR_MAX_LEN, "ERROR: unable to open Gerbicz-check savefile %s for write of checkpoint data.\n",cstr);
														fprintf(stderr,"%s",cbuf);
					fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				}
			}	// ihi a multiple of ITERS_BETWEEN_GCHECKS?
		}

		if(ierr == ERR_INTERRUPT) exit(0);

		// For Fermats and cofactor-PRP tests of either modulus type, exit only after writing final-residue checkpoint file:
		if( ihi == maxiter ) break;

		/*...reset loop parameters and begin next iteration cycle...	*/
		ilo = ihi;
		ihi = ilo+ITERS_BETWEEN_CHECKPOINTS;
	}

	/*...For timing tests, print timing and 64-bit residue and return timing via arglist runtime ptr...	*/
	*runtime = tdiff;

	/* If Selftest mode... */
	if(INTERACT) {
		fprintf(stderr, "%u iterations of %s with FFT length %u = %u K, final residue shift count = %llu\n",timing_test_iters,PSTRING,n,kblocks,RES_SHIFT);
		// If testType non-default (e.g. PRP for Mersennes), add text indicating that:
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRP)
			sprintf(cbuf,"PRP-%u ",PRP_BASE);
		else
			sprintf(cbuf,"");

		/* If Fermat number, make sure exponent a power of 2: */
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
			ASSERT(HERE, (p >> findex) == 1,"Require (p >> findex) == 1");

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

	/*...Check (probable) primality:	*/
	if(TEST_TYPE == TEST_TYPE_PRIMALITY || TEST_TYPE == TEST_TYPE_PRP)
	{
		// v19: for Mersenne-mod PRP, we use a Gerbicz-check modified-PRP test with 2 more squarings than standard Fermat-mod PRP,
		// at the end of which PRP-ness is indicated by the final residue == PRP_BASE^2. Since PRP_BASE is a uint32, the latter is
		// as large as 64 bits, thus simply check whether Res64 == PRP_BASE^2 and infer all-higher-words-0-ness by additionally
		// checking whether Res64 % 2^35-1 == Res35m1 and Res64 % 2^36-1 == Res36m1:
		if(TEST_TYPE == TEST_TYPE_PRP) {
			ASSERT(HERE, ihi == p, "Gerbicz-check-modified PRP-test requires p mod-squarings!");
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
			i2 %= itmp64;	ASSERT(HERE, i3 == 0ull, "K-multiplier needs 64-bit reduction (mod b^2)!");
			if(i2) i2 = itmp64 - i2;	// if(k) k = -r".mi" (mod b^2) = b^2 - r".mi" .
			// i2 contains the needed multiplier k. Since ensuing quotient computation needs separate arrays
			// for dividend and quotient, stash output of mi64_mul_scalar_add_vec2 in c[] and ensuing quotient back in arrtmp[]:
			c_uint64_ptr[j] = mi64_mul_scalar_add_vec2(d_uint64_ptr,i2,arrtmp, c_uint64_ptr, j);
			// Now short-div - allowing for the possibility of a carryout from above mi64_mul_scalar_add_vec2() call -
			// by base and check that remainder 0. Note that we do want the quotient now, as that is our reside/base:
			mi64_div(c_uint64_ptr, &itmp64, j+1,1, arrtmp,&rmodb);	ASSERT(HERE, rmodb == 0ull,"After short-div, R != 0 (mod B)");
			// And recompute the S-H residues:
			Res64 = arrtmp[0];
			Res35m1 = mi64_div_by_scalar64((uint64*)arrtmp,two35m1,i,0x0);
			Res36m1 = mi64_div_by_scalar64((uint64*)arrtmp,two36m1,i,0x0);
			// Now that residue is standard Fermat-PRP-test one, check if == 1:
			isprime = (arrtmp[0] == 1ull);
			if(isprime) {
				for(i = 1; i < n; i++)
				{
					j = i + ( (i >> DAT_BITS) << PAD_BITS );
					if(a[j] != 0.0) { isprime = 0; break; }
				}
			}
		} else {	// older impl. of LL-test isprime parsed the entire double-float residue array:
			isprime = 1;
			/* For Fermat numbers, in balanced-digit form, it's prime if the lowest-order digit = -1, all others 0: */
			if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
				final_res_offset = 1;
			else
				final_res_offset = 0;
			a[0] += final_res_offset;
			for(i = 0; i < n; i++)
			{
				j = i + ( (i >> DAT_BITS) << PAD_BITS );
				if(a[j] != 0.0) { isprime = 0; break; }
			}
			a[0] -= final_res_offset;
		}

		// v19: Add basic JSON-formatted result report for M-number tests:
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			// UTC is for the result as submitted to the server:
			gm_time = gmtime(&calendar_time);
			strftime(timebuffer,SIZE,"%Y-%m-%d %H:%M:%S %Z",gm_time);
			generate_JSON_report(isprime,p,n,Res64,timebuffer, cstr);	// cstr holds the formatted output line here
		}

		/*...Unbelievable - I must be hallucinating.	*/
		if(isprime)
		{
			if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
			{
				/*... this gets written both to file and to stdout, the latter irrespective of whether the run is in interactive mode...	*/
				snprintf(cbuf,STR_MAX_LEN, "%s is a new FERMAT PRIME!!!\nPlease send e-mail to ewmayer@aol.com.\n",PSTRING);
											        fprintf(stderr,"%s",cbuf);
				fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0; }
				fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0; }
			}
			else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				for(i=0; knowns[i] != 0; i++)
				{
					if(p == knowns[i]) { break; }
				}

				if(knowns[i] != 0)
				{
					snprintf(cbuf,STR_MAX_LEN, "%s is a known MERSENNE PRIME.\n",PSTRING);

					if(INTERACT || scrnFlag)	/* Echo output to stddev */
						fprintf(stderr,"%s",cbuf);
					else
					{
						fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(fp,"%s",cbuf);	fprintf(fp,"\n%s",cstr);	fclose(fp);	fp = 0x0; }
						fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(fp,"%s",cbuf);	fprintf(fp,"\n%s",cstr);	fclose(fp);	fp = 0x0; }
					}
				}
				else
				{
					/*... this gets written both to file and to stdout, the latter irrespective of whether the run is in interactive mode...	*/
					snprintf(cbuf,STR_MAX_LEN, "%s is a (probable) new MERSENNE PRIME!!!\nPlease send e-mail to ewmayer@aol.com and woltman@alum.mit.edu.\n",PSTRING);
												fprintf(stderr,"%s",cbuf);
					fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(fp,"%s",cbuf);	fprintf(fp,"\n%s",cstr);	fclose(fp);	fp = 0x0; }
					fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(fp,"%s",cbuf);	fprintf(fp,"\n%s",cstr);	fclose(fp);	fp = 0x0; }
				}
			}
			else
				ASSERT(HERE, 0, "Unsupported modulus type!");
		}
		/*
		The more likely scenario - it's not prime, so we form a 64-bit residue and write that.
		If residue has < 64 bits, print a warning.
		*/
		else {
			// Otherwise, write the 64-bit hex residue. As of v19, we write the old-style HRF-formatted result
			// just to the exponent-specific logfile, and the server-expected JSON-formatted result to the results file:
			fp = mlucas_fopen(STATFILE,"a");
			if(fp) {
				snprintf(cbuf,STR_MAX_LEN, "%s is not prime. Res64: %016llX. Program: E%s. Final residue shift count = %llu\n",PSTRING,Res64,VERSION,RES_SHIFT);
				fprintf(fp,"%s",cbuf);
				// Also write the 3 supplemental SH residues to the .stat file:
				if(!INTERACT) {
					snprintf(cbuf,STR_MAX_LEN, "%s mod 2^35 - 1 = %20.0f\n",PSTRING,(double)Res35m1);
					if(fp){ fprintf(fp,"%s",cbuf); }
					snprintf(cbuf,STR_MAX_LEN, "%s mod 2^36 - 1 = %20.0f\n",PSTRING,(double)Res36m1);
					if(fp){ fprintf(fp,"%s",cbuf); }
				}
				fclose(fp);	fp = 0x0;
			}
			// v19: Finish with the JSON-formatted result line:
			fp = mlucas_fopen(OFILE,"a");
			if(fp){ fprintf(fp,"\n%s",cstr);	fclose(fp);	fp = 0x0; }
		}
	} else {
		ASSERT(HERE, 0, "Unrecognized test type!");
	}	/* endif(TEST_TYPE == TEST_TYPE_PRIMALITY) */

	/*...If successful completion, delete the secondary restart files...save the primary in case it's a prime,
	so we can re-run the last p%2000 iterations by way of quick partial verification:
	*/
	RESTARTFILE[0] = 'q';
	if(remove(RESTARTFILE))
		fprintf(stderr, "Unable to delete secondary restart file %s\n",RESTARTFILE);
	RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');

	/*...If in non-interactive (range test) mode, delete the just-completed exponent from the rangefile,
	or (in the case of Mersenne primality pre-testing), munge the assignment line to reflect that trial
	factoring or p-1 factoring has completed.	*/
	/*
		Copy all the currently assigned exponents (except the first) from the worktodo.ini file to a temporary file...
	*/
GET_NEXT_ASSIGNMENT:

	if(!INTERACT)
	{
		fp = mlucas_fopen(RANGEFILE,"r");
		if(!fp)
		{
			sprintf(cbuf,"ERROR: unable to open %s file for reading.\n",RANGEFILE);
			fp = mlucas_fopen(OFILE,"a");
			if(fp)
			{
				fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
			}
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		/* Remove any WINI.TMP file that may be present: */
		remove("WINI.TMP");
		fq = mlucas_fopen("WINI.TMP", "w");
		if(!fq)
		{
			fprintf(stderr, "Unable to open WINI.TMP file for writing.\n");
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		/* Delete or suitably modify current-assignment line (line 1) of worktodo.ini file: */
		if(!fgets(in_line, STR_MAX_LEN, fp))
		{
			sprintf(cbuf, "ERROR: %s file not found at end of current-assignment processing\n", RANGEFILE);
			fprintf(stderr,"%s",cbuf);
			fp = mlucas_fopen(OFILE, "a");	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
			ASSERT(HERE,0,"");
		}
		else if(!strstr(in_line, ESTRING))
		{
			snprintf(cbuf,STR_MAX_LEN, "WARNING: Current exponent %s not found in line 1 of %s file - skipping line deletion\n", ESTRING, RANGEFILE);
			fprintf(stderr,"%s",cbuf);
			fp = mlucas_fopen(OFILE, "a");	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
			/* In this case we copy line 1 to WINI.TMP, rather than deleting it: */
			fputs(in_line, fq);
		}
		else
		{
			/* If we just finished the TF or p-1 preprocessing step of an LL test,
			update the current-assignment line to reflect that and write it out: */
			if(strstr(in_line, "Test") || strstr(in_line, "DoubleCheck"))
			{
				if(TEST_TYPE == TEST_TYPE_TRIALFACTORING)
				{
					/* Factor depth assumed to follow the first comma in in_line: */
					char_addr = strstr(char_addr, ",");
					ASSERT(HERE, char_addr != 0x0,"Null char_addr");
					sprintf(++char_addr, "%u", bit_depth_done);
					fputs(in_line, fq);
				}
				else if(TEST_TYPE == TEST_TYPE_PMINUS1)
				{
					/*0/1 flag indicating whether P-1 has been done assumed to follow second comma in in_line: */
					char_addr = strstr(char_addr, ",");
					ASSERT(HERE, char_addr != 0x0,"Null char_addr");
					char_addr++;
					char_addr = strstr(char_addr, ",");
					ASSERT(HERE, char_addr != 0x0,"Null char_addr");
					sprintf(++char_addr, "1");
					fputs(in_line, fq);
				}
			}
			/* Otherwise lose the current line (by way of no-op) */
		}
		/* Copy the remaining ones; */
		i = 0;
		while(fgets(in_line, STR_MAX_LEN, fp))
		{
			fputs(in_line, fq);	++i;
		}
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;

		/* Now blow away the old worktodo.ini file and rename WINI.TMP ==> worktodo.ini...	*/
		remove(RANGEFILE);
		if(rename("WINI.TMP", RANGEFILE))
		{
			sprintf(cbuf,"ERROR: unable to rename WINI.TMP file ==> %s ... attempting line-by-line copy instead.\n",RANGEFILE);
			fprintf(stderr,"%s",cbuf);
			fp = mlucas_fopen(OFILE,"a");	if(fp){	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;	}

			/* If attempting to simply rename the TMP file fails, do it the hard way: */
			fp = mlucas_fopen(RANGEFILE,"w");
			if(!fp)
			{
				sprintf(cbuf,"ERROR: unable to open %s file for writing.\n", RANGEFILE);
				fp = mlucas_fopen(OFILE,"a");
				if(fp)
				{
					fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
				}
				if(scrnFlag)		/* Echo output to stddev */
				{
					fprintf(stderr,"%s",cbuf);
				}
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}

			fq = mlucas_fopen("WINI.TMP", "r");
			if(!fq)
			{
				sprintf(cbuf,"Unable to open WINI.TMP file for reading.\n");
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			while(fgets(in_line, STR_MAX_LEN, fq))
			{
				fputs(in_line, fp);
			}
			fclose(fp);	fp = 0x0;
			fclose(fq); fq = 0x0;

			/*...Then remove the WINI.TMP file:	*/

			remove("WINI.TMP");
		}
		/* if one or more exponents left in rangefile, go back for more; otherwise exit. */
		if (i > 0)
			goto RANGE_BEG;   /* CYCLE RANGE */
		else
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
#ifdef INCLUDE_TF
	/* TF currently only supported for Mersennes: */
	ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_MERSENNE   ][TEST_TYPE_TRIALFACTORING] = TRUE;
#endif
	ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_MERSENNE   ][TEST_TYPE_PRIMALITY     ] = TRUE;
	ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_FERMAT     ][TEST_TYPE_PRIMALITY     ] = TRUE;

	/* Simple self-tester for get_fft_radices(): */
	fprintf(stderr, "INFO: testing FFT radix tables...\n");
	test_fft_radixtables();

	/* Set min. exponent (in terms of power of 2) that can be tested: */
	/* Check that the purported min. FFT length is actually supported: */
	ASSERT(HERE, get_fft_radices(MIN_FFT_LENGTH_IN_K, 0, 0x0, 0x0, 0) == 0,"Require get_fft_radices(MIN_FFT_LENGTH_IN_K, 0) == 0");
	n = (MIN_FFT_LENGTH_IN_K << 10);
	/* Make sure N didn't overflow */
	ASSERT(HERE, (n >> 10) == MIN_FFT_LENGTH_IN_K,"Require (n >> 10) == MIN_FFT_LENGTH_IN_K");
	PMIN = 2*n;	/* 2 bits per input is about the smallest we can test without getting nonzero-carry errors */

	/* Set max. exponent (in terms of power of 2) that can be tested: */
	/* Check that the purported max. FFT length is actually supported: */
	ASSERT(HERE, get_fft_radices(MAX_FFT_LENGTH_IN_K, 0, 0x0, 0x0, 0) == 0,"Require get_fft_radices(MAX_FFT_LENGTH_IN_K, 0) == 0");
	n = (MAX_FFT_LENGTH_IN_K << 10);
	/* Make sure N didn't overflow */
	ASSERT(HERE, (n >> 10) == MAX_FFT_LENGTH_IN_K,"Require (n >> 10) == MAX_FFT_LENGTH_IN_K");
	PMAX = given_N_get_maxP(n);

	ASSERT(HERE, PMAX > PMIN,"Require PMAX > PMIN");

#ifdef INCLUDE_TF
	/* Simple self-tester for sieve factoring routines: */
	fprintf(stderr, "INFO: testing trial-factoring routines...\n");
	if(test_fac() != 0)
	{
		sprintf(cbuf, "Mlucas_init : Trial-factoring self-test failed.\n");
		fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}
#endif
#ifdef INCLUDE_PM1
	/* Simple self-tester for GCD routines in gcd_lehmer.c: */
	fprintf(stderr, "INFO: testing GCD routines...\n");
	if(test_gcd() != 0)
	{
		sprintf(cbuf, "Mlucas_init : GCD test failed.\n");
		fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
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
	int pow2_fft,bimodn,curr_bit64,curr_wd64,w64,curr_wd_bits,mod64,findex,i,ii,j,j1,j2;
	uint32 bw = p%n,sw = n-bw,nwt,sw_div_n,base[2],baseinv[2],bits[2],sw_idx_modn;
	double cy,theta,wt_fwd,wt_cos,wt_sin;
	uint64 nbits, itmp64;
	 int64 retval = -1;	// Make this signed to ease "not yet set?" check
#ifdef USE_FGT61
	ASSERT(HERE,0,"shift_word() needs to be modified to support FGT!");
#endif
	if(n != nsave || p != psave) {
		first_entry = TRUE;	for(j = 0; j < (n>>6); j++) { BIGWORD_BITMAP[j] = 0ull; }	// Need to clear bitmap in case of multi-FFT-length run
	}
	// May 2018: Check that BIGWORD_BITMAP and BIGWORD_NBITS arrays have been alloc'ed and use fast lookup based on those.
	// Timing loop on my Core2 macbook indicates this fast-lookup needs ~160 cycles @1Mdouble FFT, not horrible but slower
	// than I'd like, likely due to cache impacts of doing random-word lookups in the resulting 128kB and 64kB BIGWORD* arrays.
	// Also had the "adjusting..." printfs enabled during the timing tests, 0 such adjustments needed for 10^9 random-shifts:
	if(!first_entry) {
	//	ASSERT(HERE, BIGWORD_BITMAP != 0x0 && BIGWORD_NBITS != 0x0, "BIGWORD_BITMAP and BIGWORD_NBITS arrays not alloc'ed!");
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
		// Any carry-in needs to be properly left-shifted w.r.to the residue word into which it goes.
		// Since the stdlib pow() function is expensive (~100-200 cycles), wrap this bit in an if(cy_in):
		if(cy_in) {
			// Compute wt = 2^(j*sw % n)/n:
			sw_idx_modn = ((uint64)j*sw) % n;	// N is 32-bit, so only use 64-bit to hold intermediate product
			wt_fwd = pow(2.0, (double)sw_idx_modn/n);
			a[j1] = cy_in * (double)(1ull << (shift - ii)) * wt_fwd;
		}
		retval = ((uint64)j<<8) + (shift - ii);
		goto DONE;
	}

	first_entry = FALSE;
	psave = p; nsave = n; bits_per_word = (double)p/n; words_per_bit = 1.0/bits_per_word;
	ASSERT(HERE,MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(HERE,MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");
	ASSERT(HERE,TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(HERE,TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");
	ASSERT(HERE,p > shift,"Specified shift count out of range!");

	/* If Fermat number, make sure exponent a power of 2: */
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
		ASSERT(HERE, TRANSFORM_TYPE == RIGHT_ANGLE,"Require TRANSFORM_TYPE == RIGHT_ANGLE");
		findex = trailz64(p);
		ASSERT(HERE, (p >> findex) == 1,"Require (p >> findex) == 1");
		/* For Fermat-mod, only need IBDWT weights table if it's a non-power-of-2-length transform, in which
		case the table has {nwt = odd part of N} distinct elements. Avoid if() logic related to power-of-2-or-not
		by initing a single DWT weight = 1.0 in the power-of-2 case and = 2^((j%nwt)/n) otherwise:
		*/
		nwt = (n >> trailz32(n));
		sw_div_n = sw*nwt/n;
	}
	else
		ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER,"Require TRANSFORM_TYPE == REAL_WRAPPER");

	/* Vector length a power of 2? */
	pow2_fft = (n >> trailz32(n)) == 1;

	bits[0] = p/n;
	base[0] = 1 << bits[0];	baseinv[0] = 1.0/base[0];

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT && pow2_fft == TRUE)
		bits[1] =     bits[0];
	else
		bits[1] = 1 + bits[0];

	base[1] = 1 << bits[1];	baseinv[1] = 1.0/base[1];

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
				//	printf("Hit target bit %llu in a[%u] (=> BIGWORD_BITMAP[%u]), bit %u of <0:%u>, bitmap-word bit = %u\n",shift,j,curr_wd64,curr_wd_bits,bits[ii]-1,curr_bit64-1);	ASSERT(HERE, curr_wd_bits <= bits[ii]-1,"GAH!");
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
		ASSERT(HERE, TRANSFORM_TYPE == RIGHT_ANGLE, "Invalid or uninited TRANSFORM_TYPE!");
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
				//	printf("Hit target bit %llu in a[%u] (=> BIGWORD_BITMAP[%u]), bit %u of <0:%u>, bitmap-word bit = %u\n",shift,j,curr_wd64,curr_wd_bits,bits[ii]-1,curr_bit64-1);	ASSERT(HERE, curr_wd_bits <= bits[ii]-1,"GAH!");
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


/******************Thanks to Tom Cage for the initial version of this: ************************/

#ifdef macintosh
	#include <console.h>	/* Macintosh */
#endif

/* Number of distinct FFT lengths supported for self-tests: */
#define numTest				120	// = sum of all the subranges below
/* Number of FFT lengths in the various subranges of the full self-test suite: */
#define numTiny 			32
#define numSmall			32
#define numMedium			16	// v18: Moved 8 from medium to small to increase starting FFTlen of -m selftests from 1024K to 2048K
#define numLarge			24
#define numHuge				16
/* Adding larger FFT lengths to test vectors requires supporting changes to Mdata.h:MAX_FFT_LENGTH_IN_K and get_fft_radices.c */
#define numEgregious		 0
#define numBrobdingnagian	 0
#define numGodzillian		 0
#if(numTiny+numSmall+numMedium+numLarge+numHuge+numEgregious+numBrobdingnagian+numGodzillian != numTest)
	#error Sum(numTiny+...) != numTest in main!
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
	uint32 exponent;	/* Test exponent for this FFT length. */
	struct res_triplet	res_t[3];	/* 100,1000 and 10000-iteration SH-residue triplets */
};

// Reference LL-test residues:
struct testMers MersVec[numTest+1] =
{
/*                                         100-iteration residues:	                               1000-iteration residues:                */
/*	  FFTlen(K)     p              Res64           mod 2^35-1      mod 2^36-1               Res64           mod 2^35-1      mod 2^36-1     */
/*	    -----    --------     ----------------     -----------     -----------         ----------------     -----------     -----------    */
	/* Tiny:                                     [%34359738367  ][%68719476735  ]                         [%34359738367  ][%68719476735  ] */
	{     8,    173431u, { {0x85301536E4CA9B11ull,  3707224323ull, 36851834664ull}, {0x2FD5120BEC41F449ull, 28734955954ull, 23103392631ull}, {0x139D1D396F173696ull, 12716541843ull, 58117202214ull} } },
	{     9,    194609u, { {0xC711AF1008612BC6ull,  1574019740ull, 37260026270ull}, {0x5153F6E040CD1BE6ull, 15446410924ull,  3291404673ull}, {0x33E19077F35070A3ull, 34231248547ull, 24411805292ull} } },
	{    10,    215767u, { {0x4428783BC62760F0ull,  7466284975ull, 53916123655ull}, {0xED46A8C001908815ull,   739143119ull, 36950829937ull}, {0xCBE0AD544E96FDB9ull,  7625169722ull, 52027104104ull} } },
	{    11,    236813u, { {0x592D849AF4D1336Full, 29025996994ull, 48905971124ull}, {0xB4EEB63BB656F424ull,  5361680901ull, 31850818767ull}, {0x5AEBAE493B085903ull, 20240457506ull, 42866015746ull} } },
	{    12,    257903u, { {0x1D8121DE28B60996ull, 22402402961ull, 65583959369ull}, {0x54F2BE961A674CB1ull, 25601315553ull, 54298899520ull}, {0x50DFBA28D0D1A8C3ull, 17306049864ull, 68068537809ull} } },
	{    13,    278917u, { {0xE3BC90B0E652C7C0ull, 21244206101ull, 51449948145ull}, {0x93AF8994F95F2E50ull, 16427368469ull, 10707190710ull}, {0x1674AAA04F7BD61Aull, 12079507298ull, 56593045102ull} } },
	{    14,    299903u, { {0xDB8E39C67F8CCA0Aull, 20506717562ull, 44874927985ull}, {0x4E7CCB446371C470ull, 34135369163ull, 61575700812ull}, {0x04ACC83FFE9CEAD4ull, 26179715264ull, 65445483729ull} } },
	{    15,    320851u, { {0xB3C5A1C03E26BB17ull, 22101045153ull,  4420560161ull}, {0x923A9870D65BC73Dull, 29411268414ull, 30739991617ull}, {0xB3F1ACF3A26C4D72ull, 32179253815ull, 68615042306ull} } },
	{    16,    341749u, { {0x8223DF939E46A0FFull, 32377771756ull, 38218252095ull}, {0xC6A5D4B6034A34B8ull, 31917858141ull, 59888258577ull}, {0x93EF44581866E318ull, 18805111197ull,  8333640393ull} } },
	{    18,    383521u, { {0xBF30D4AF5ADF87C8ull, 15059093425ull, 52618040649ull}, {0x9F453732B3FE3C04ull,  4385160151ull, 47987324636ull}, {0x0DBF50D7F2142148ull,  1608653720ull, 52016825449ull} } },
	{    20,    425149u, { {0x6951388C3B99EEC0ull,  4401287495ull, 19242775142ull}, {0x501CEC2CB2080627ull, 21816565170ull, 41043945930ull}, {0x5A9A9BF4608090A2ull, 27025233803ull, 68581005187ull} } },
	{    22,    466733u, { {0xD95F8EC0F32B4756ull, 19305723506ull, 26588871256ull}, {0xB1F58184918D94B6ull,  8443388060ull, 11738516313ull}, {0xAC4B1F499BF2C2DAull,  7322105347ull, 15747709958ull} } },
	{    24,    508223u, { {0xDA46E41316F8BCCAull, 25471180026ull,  1635203275ull}, {0x27A5B285281466B9ull, 11438869313ull,  7226774009ull}, {0x4ABED2868B800F7Dull,  7783533092ull, 66921500486ull} } },
	{    26,    549623u, { {0x6649D9D6CD4E0CE1ull, 25445908581ull, 26118212198ull}, {0x1A4F280627A15B3Cull, 13286323782ull, 31550278005ull}, {0x86404E236E99B3C4ull, 17401894517ull, 40934891751ull} } },
	{    28,    590963u, { {0x4ADDB6C4A76465AFull,  6532108269ull, 54921134131ull}, {0x3063D08A7BABD7B8ull,  4777711548ull, 39733274344ull}, {0xBE2ABBB09336F32Eull, 30656127523ull, 50296089656ull} } },
	{    30,    632251u, { {0x0811FAA40601EB1Dull, 16369365746ull,  6888026123ull}, {0xF324E4DEC564AF91ull, 10236920023ull, 34068699974ull}, {0xAA622CF2A48F6085ull, 22315931502ull,  1049914969ull} } },
	{    32,    673469u, { {0x1A4EF8A0D172FBAAull, 32667536946ull, 11393278588ull}, {0xA4DFD62B928F68A4ull, 11900420802ull, 66610946021ull}, {0xFA3993AC9CE7BEEDull,   117685830ull, 39361570554ull} } },
	{    36,    755737u, { {0x13B13C61298088DCull, 34092803628ull,  7584858890ull}, {0x33A2A43DE8782CCCull,  2953124985ull, 62434716987ull}, {0xD0DF76911349551Bull, 28919846011ull, 30432969648ull} } },
	{    40,    837817u, { {0x88555D9AAD3FF2DDull,  8573348747ull, 67896670216ull}, {0xEAC1676D914878C0ull, 34312095136ull, 45077378164ull}, {0x89E1C4D06BB0F9F3ull,  6272358557ull, 24712951618ull} } },
	{    44,    919729u, { {0x6ACC03213A37BA5Bull,  3870201113ull, 48145739792ull}, {0xDA98B49CC83C60CBull, 15886769401ull, 62221100895ull}, {0xF1DA20E7D5A89638ull,  4633752262ull, 20941274692ull} } },
	{    48,   1001467u, { {0x6B1C76AB5431FDA4ull,  6795121241ull, 65308927583ull}, {0xBD99FD21F4136BFCull, 26386988063ull, 61607603549ull}, {0x9136554E4718BFA9ull, 18451197932ull, 37688798842ull} } },
	{    52,   1083077u, { {0xA591637EC8CF3FE4ull,  4769775755ull, 65114872367ull}, {0xE59C08B13B00E6FFull,  1383963096ull, 26100699764ull}, {0x48CCA2242A1F9352ull, 30318043361ull, 12067176371ull} } },
	{    56,   1164533u, { {0xEC4F2579E4533584ull,  5456769127ull, 59922459736ull}, {0xF7D2BF94C2767D36ull, 30727892629ull, 48141128220ull}, {0xE332E3891AE98AD7ull,  7024598607ull, 65691841143ull} } },
	{    60,   1245877u, { {0xC91002E1A4EE7E07ull,  6217476228ull, 40164514288ull}, {0xEABE9E1A31DF5877ull,   831216169ull, 29591771932ull}, {0xCB85101F6857519Dull, 30425930108ull,  2198194326ull} } },
	{    64,   1327099u, { {0xAC070112281229E0ull, 14226353454ull,  1640524016ull}, {0xF25AA54053C5BB64ull, 32455038659ull, 53547160776ull}, {0x3854D019CE12CC9Aull, 29589836279ull,  2174826233ull} } },
	{    72,   1489223u, { {0x6674518EA19B3D6Aull, 32383400562ull, 53234746310ull}, {0xEB312091097F6C3Bull,  3980687908ull,  8568698675ull}, {0xC225FAF24A093590ull, 22407813999ull, 30924932017ull} } },
	{    80,   1650959u, { {0xE5326E754F3202A8ull,  5593233426ull, 33337128557ull}, {0xFC3E8CDA60AF5CF8ull, 11466296968ull, 12651602524ull}, {0x4D68554B73674A60ull,  2253999911ull, 55045374456ull} } },
	{    88,   1812347u, { {0x81BDD3AC63DF3F73ull, 19957199947ull, 61002681293ull}, {0x3D3E429D7427C4EAull, 25342898119ull, 34322985438ull}, {0x5769D7B47C49436Full,  5234049262ull, 26872574292ull} } },
	{    96,   1973431u, { {0x901C8305DA9FF95Aull, 32611811878ull, 55986702160ull}, {0x0790CA11ADAA47E3ull, 17075140845ull, 12883521448ull}, {0xF18EADC267DE6FC1ull, 24308841307ull, 31678116890ull} } },
	{   104,   2134201u, { {0x59BDA0D80F3279EDull, 17901153436ull,  3927067335ull}, {0x2F81B21BC680C861ull, 18443771511ull, 45465079919ull}, {0x439245FA16A38116ull, 20996570088ull,   489289103ull} } },
	{   112,   2294731u, { {0xC44ACC96D268625Full, 10331638988ull,  2292055445ull}, {0xED20577E16E128DEull, 32248607028ull, 14903460370ull}, {0xCB862A1B42B230A2ull, 23316229090ull, 23891565685ull} } },
	{   120,   2455003u, { {0xC5F7DB23F174A67Dull, 32991574397ull, 31642856976ull}, {0x401670254012E5ABull, 33626385418ull, 66465546971ull}, {0x20AB396E327C09C1ull, 13309965383ull, 60492105240ull} } },
	/* Small: */
	{   128,   2614999u, { {0x040918890E98F8DAull, 14867710211ull, 47602627318ull}, {0x1A184504D2DE2D3Cull,  5934292942ull,  4090378120ull}, {0xE7126F512D3FD742ull, 17101849610ull, 66501661438ull} } },
	{   144,   2934479u, { {0x1B90A27301980A3Aull,  7043479338ull, 38327130996ull}, {0x8C3045C6534867C6ull, 12456621644ull, 52801948293ull}, {0xF17F4A594A281B94ull,  5970782987ull, 68371435254ull} } },
	{   160,   3253153u, { {0x9AFD3618C164D1B4ull, 16551334620ull, 55616214582ull}, {0x1493A70897A8D058ull, 34082962858ull, 60773088284ull}, {0x57D3F1A090E78729ull, 26902546905ull, 49396480035ull} } },
	{   176,   3571153u, { {0xA016F25779902477ull, 21500047857ull,  9150810891ull}, {0x8E6F248EC96445FFull, 22443629034ull, 16625023722ull}, {0xFFD8B840C06A5EACull, 32183737848ull, 42577197579ull} } },
	{   192,   3888509u, { {0x71E61322CCFB396Cull, 29259839105ull, 50741070790ull}, {0x3CEDB241702D2907ull,  6177258458ull, 21951191321ull}, {0xAF407D11B2D74C3Cull, 31653650180ull, 27299459944ull} } },
	{   208,   4205303u, { {0xC08562DA75132764ull,  7099101614ull, 36784779697ull}, {0xAD381B4FE91D46FDull,  7173420823ull, 51721175527ull}, {0xC70061EF9537C4E1ull,  9945894076ull,  2301956793ull} } },
	{   224,   4521557u, { {0xE68210464F96D6A6ull, 20442129364ull, 11338970081ull}, {0x3B06B74F5D4C0E35ull,  7526060994ull, 28782225212ull}, {0xB720ACD1D69A7ECFull, 28103212586ull, 10983125296ull} } },
	{   240,   4837331u, { {0xB0D0E72B7C87C174ull, 15439682274ull, 46315054895ull}, {0x3AA14E0E90D16317ull,  5730133308ull, 50944347816ull}, {0x12CFBF6001E59FF7ull, 26877054587ull, 60322521357ull} } },
	{   256,   5152643u, { {0x074879D86679CB5Bull,  1208548964ull, 48525653083ull}, {0x98AF5E14C824A252ull,   783196824ull,  6594302302ull}, {0x7DA0D3B9EFEA4931ull, 32608975347ull, 43428286760ull} } },
	{   288,   5782013u, { {0x9869BE81D9AB1564ull, 15509103769ull, 49640026911ull}, {0x7C998719C6001318ull, 23749848147ull, 19853218689ull}, {0xE2E246D9094EBFD7ull, 26657044660ull,  7091330955ull} } },
	{   320,   6409849u, { {0x20739E43A693A937ull, 27970131351ull, 15334307151ull}, {0xE20A76DCEB6774A6ull, 14260757089ull, 68560882840ull}, {0xCEC786F8883D8D1Full,  5597853948ull, 57984323163ull} } },
	{   352,   7036339u, { {0xD6A226BAB5E14D62ull,  9444972171ull, 28639488018ull}, {0x0579D28296F29D92ull, 18964853245ull, 30111685201ull}, {0xB9CF9FE489BD34CFull, 11297283566ull, 16782498229ull} } },
	{   384,   7661567u, { {0x3A929F577AC9725Full, 23890492835ull, 64821764767ull}, {0x2ECBA785576E6D58ull, 26446200615ull, 60269798452ull}, {0x22AA4A0A7A3676A7ull, 23262227373ull, 26736591016ull} } },
	{   416,   8285659u, { {0xDCA138D55C36E40Cull, 10452583294ull,  4308453248ull}, {0x1FEE7F79E32229A6ull, 11936075329ull, 16061515794ull}, {0x9B5324C99CCE498Dull, 14697684311ull, 16439873760ull} } },
	{   448,   8908723u, { {0x436494C7EA194FA1ull,  2976149044ull, 21645125251ull}, {0xD2FCEDE29E818A26ull, 15260150529ull, 11678366985ull}, {0x952F59F4972F830Aull, 17286796157ull, 11216464341ull} } },
	{   480,   9530803u, { {0x6D52BCCB796A46E9ull, 28296985362ull, 66294636306ull}, {0x786CB610A809B762ull, 24654197494ull, 57943258783ull}, {0x40564770A8540610ull, 27197548273ull, 32971949076ull} } },
	{   512,  10151971u, { {0x7BA3DD9C38878B83ull, 24395728676ull, 12632037498ull}, {0x8E9FA3093ACD81C1ull,  6345070464ull, 65203567231ull}, {0xBBD3BD10BA9983B5ull,  7822185024ull, 14668688365ull} } },
	{   576,  11391823u, { {0x78D8769B9F75FB2Bull,  7286184017ull, 17306427758ull}, {0x5834BDA7558DD43Cull, 11189092321ull, 23026923116ull}, {0x793D548FCB76B28Cull, 11310068490ull, 24716315416ull} } },
	{   640,  12628613u, { {0xF951B697F5C5032Aull,  9262385660ull, 57463106915ull}, {0x93B526040205BA27ull, 30538079080ull, 32317022014ull}, {0x25544FC69ADB28F9ull, 28253742327ull,  1823182110ull} } },
	{   704,  13862759u, { {0xBB2F69275D79A9EEull, 12230183669ull, 68684647134ull}, {0x7343ECC160AA00D5ull, 24655585170ull, 51102704879ull}, {0x1402EEF49394CDC7ull,  5500341204ull, 59999916295ull} } },
	{   768,  15094403u, { {0xF6895EB66EADE9C5ull,  8490184692ull, 23393343807ull}, {0xF673A8D6413923A9ull, 20026779905ull, 67516766223ull}, {0xD63752CA13598971ull, 24773095342ull, 29303310893ull} } },
	{   832,  16323773u, { {0xEB8890F379392B2Full, 27289972116ull, 63975275393ull}, {0xD681EDD3A1EC3780ull, 12515962698ull, 40155157152ull}, {0xEB9C9477368BF584ull, 13378242091ull,  9365072054ull} } },
	{   896,  17551099u, { {0xAB1180428ED65EE0ull,  3105108668ull, 66518734167ull}, {0x31813367849BBF49ull,  9516734777ull, 18271834608ull}, {0x95C2E1F201FCE598ull,  6264820675ull, 49312303312ull} } },
	{   960,  18776473u, { {0xCA7D81B22AE24935ull, 24317941565ull, 67706175547ull}, {0x02EB980A49E7B60Full,  5730644436ull, 48386545950ull}, {0xBC6503AA5C062308ull, 29760131532ull, 31603724687ull} } },
	{  1024,  19800083u, { {0x95AFD7A5269F14F6ull, 26677613826ull, 37493068952ull}, {0x18A53602BAA0E197ull,   624371978ull, 15180896714ull}, {0xF79CD0274644183Dull, 21538070258ull, 26190157173ull} } },
	{  1152,  22217791u, { {0x9EFE3D8C08D89E48ull,  8307313500ull, 40673806995ull}, {0xE1F94CD14457EDADull, 23322294478ull, 42677160325ull}, {0xF60CFBDEA4ADF55Full,  5469936596ull, 35790203222ull} } },
	{  1280,  24629621u, { {0x85D92483D90E5029ull,  6310387878ull, 18231127032ull}, {0xBEE63CF182681345ull,  7247112769ull, 24502530130ull}, {0x89977592EE68F853ull, 31511858957ull, 44154237710ull} } },
	{  1408,  27036157u, { {0x3BAC320B542307B0ull, 28733834194ull, 35177623244ull}, {0xBAE66098A85CD960ull, 25462246119ull, 52444387524ull}, {0x0569BF7834267F10ull, 27834108530ull, 44136002879ull} } },
	{  1536,  29437799u, { {0x9E279D6E3752A61Cull,  6467723895ull, 50774708769ull}, {0x5A2D06C05D222BF2ull, 22780146874ull, 60521558933ull}, {0xC6D352199D949DCDull, 30008450771ull, 37092077102ull} } },
	{  1664,  31835017u, { {0xE10EFEAEDCF46110ull, 21648491345ull, 41207849648ull}, {0x8717F387BB55E1ECull, 17021700047ull, 59165499154ull}, {0x5CDA924B80871209ull, 22299652758ull, 41171353038ull} } },
	{  1792,  34573867u, { {0x9FC8394655E0334Eull,  1603847275ull, 51947401644ull}, {0x1128E4676929F8C8ull, 33342766910ull, 55912489998ull}, {0xE697B18729853BEEull, 25335158858ull, 63783469524ull} } },
	{  1920,  36617407u, { {0xB7FA68741ABA807Aull,  2831791183ull, 44591522415ull}, {0x0F635E0B2DC95F81ull, 24912661453ull, 45929081959ull}, {0x2B5BCFC6BA94177Aull,  7824653065ull, 62951001255ull} } },
	/* Medium: */
	{  2048,  39003229u, { {0xEC810981F56D5EC7ull, 29671814311ull, 35851198865ull}, {0xC7979AEE894F6DDEull,  6017514065ull, 41670805527ull}, {0x4669B1DD352C46BCull,  4806501770ull, 33957162576ull} } },
	{  2304,  43765019u, { {0x672821643AEE9552ull, 19240824825ull,  8017358641ull}, {0x90CF100A46BE5B64ull, 24299840772ull, 29209909852ull}, {0xA280EC055C5ADAE2ull, 28830689436ull, 48353443940ull} } },
	{  2560,  48515021u, { {0x015D3A5DC74024D1ull,  4691430475ull, 36435156228ull}, {0xFCFC9219B830DA28ull, 11946340838ull, 45970544027ull}, {0x88EBED279C6A1DE3ull,  9620966362ull,  3964242016ull} } },
	{  2816,  53254447u, { {0x75F9FC3C84FF5B40ull,   745331612ull, 37955073618ull}, {0xBECD92EA2134D5FEull, 18833421072ull, 15730516149ull}, {0x57B3BAD33AFBC633ull, 26171947880ull, 67540452305ull} } },
	{  3072,  57984131u, { {0xC7632FA366FCA848ull, 19810840296ull,  5919463661ull}, {0x678F955DDBA52AA4ull, 11995136918ull, 51932258375ull}, {0xE24A37C7B3884058ull, 17264313237ull,  9074782273ull} } },
	{  3328,  62705077u, { {0xE2CC67B70119ECA9ull,  7229868717ull, 32316011766ull}, {0xD2F94003F0B47BBCull, 19281845997ull, 52908700180ull}, {0xAF81565673A4D22Bull, 11464285328ull, 30793044388ull} } },
	{  3584,  67417873u, { {0x536142387279C6B7ull, 16480189705ull, 28663441842ull}, {0x4398ADDF46CE1F19ull, 20678941264ull, 64584283338ull}, {0xB88E95C5FA3D4C5Bull,  3311195148ull,  4541514528ull} } },
	{  3840,  72123137u, { {0x99D05C8A80D2C4ABull, 23223478983ull,  7287108796ull}, {0x48F0329C29D5464Bull,  7693484624ull, 34460970218ull}, {0xA1EA526351C5A64Eull, 10198942223ull, 29783198404ull} } },
	{  4096,  76821337u, { {0x18F108A0AEC92DA2ull, 25631216955ull, 19217786538ull}, {0x91B09D853C38FC2Bull,  1172724107ull, 62760945815ull}, {0xFB663A86824AAFC5ull,  9621633845ull, 68132687365ull} } },
	{  4608,  86198291u, { {0x238B15CD7C7C8936ull, 34119649670ull, 50469517589ull}, {0xF1A033E0DEC1330Eull, 16165796253ull, 41865067374ull}, {0x2F1992A097857555ull,  7327358291ull, 32105027991ull} } },
	{  5120,  95551873u, { {0x4AA89191484E7B56ull, 23537039943ull, 42089192670ull}, {0x811090D60FA9522Cull,  9010676320ull, 44606915786ull}, {0x7BE331521D4C32C2ull, 29429879592ull, 10276666456ull} } },
	{  5632, 104884309u, { {0x1B2A391D6B7404CDull,  7884768330ull, 55887740575ull}, {0x77CA1E7005AFA4E0ull, 26306026902ull, 21589856845ull}, {0x0F08D303912B27B4ull,   630913159ull, 40892070508ull} } },
	{  6144, 114197579u, { {0x6903363805E74E29ull,  6950779993ull, 63247095474ull}, {0xEB9D9802854B00A5ull, 24150112753ull, 67039249183ull}, {0xB64714689FA05366ull, 28265310306ull, 41895437602ull} } },
	{  6656, 123493333u, { {0x651AA3D64C84FA54ull,  9429656635ull, 48354702979ull}, {0x967C90E697CCE6D9ull,  6779392734ull, 18484099736ull}, {0x6C9511DD6EA12528ull, 27647756232ull, 21104526614ull} } },
	{  7168, 132772789u, { {0xDD02AEFE839F92D5ull,  7411321303ull, 16339659737ull}, {0xED6E26868AC2833Eull, 14154101692ull, 46327957293ull}, {0x80610E8FC3EB92E2ull, 19290762572ull, 46994666267ull} } },
	{  7680, 142037359u, { {0x9CD0C494D16CB432ull, 23235865558ull, 14066262122ull}, {0xFD6240B21A370394ull, 16216979592ull, 44514519060ull}, {0x097C240EA1436743ull,  5457504643ull, 58797441684ull} } },
	/* Large: */
	{  8192, 152816047u, { {0xB58E6FA510DC5049ull, 27285140530ull, 16378703918ull}, {0x75F2841AEBE29216ull, 13527336804ull,   503424366ull}, {0x99F8960CD890E06Aull,  8967321988ull, 43646415661ull} } },
	{  9216, 171465013u, { {0x60FE24EF89D6140Eull, 25324379967ull,  3841674711ull}, {0x6753411471AD8945ull, 17806860702ull,  3977771754ull}, {0xED3635BF88F37FEFull,  7478721112ull, 47452797377ull} } },
	{ 10240, 190066777u, { {0x65CF47927C02AC8Eull, 33635344843ull, 67530958158ull}, {0xBADA7FD24D959D21ull, 12777066809ull, 67273129313ull}, {0x82F65495D24A985Full, 22254800275ull, 49183722280ull} } },
	{ 11264, 208626181u, { {0x6FC0151B81E5173Full, 29164620640ull, 19126254587ull}, {0xD74AA66757A5345Eull, 17524190590ull, 14029371481ull}, {0xDCF9ED39C7EB15B8ull, 34266921309ull, 65896285387ull} } },
	{ 12288, 227147083u, { {0xE01AE9C859ADB03Aull,  7273133358ull,   681418986ull}, {0x303F142E1E88D5B4ull, 28479237457ull, 42044197589ull}, {0x3102781BC131D263ull, 24437355640ull, 48518577431ull} } },
	{ 13312, 245632679u, { {0x0A6ACB405ADC0354ull,    39452330ull, 38999048555ull}, {0xB38B02A4F195762Full,  3280152282ull, 30314100936ull}, {0xF020F5041AE2CABEull, 24388185991ull, 16285954298ull} } },
	{ 14336, 264085733u, { {0x5ACE4CCE3B925A81ull,  4584210608ull, 36618317213ull}, {0x02F5EC0CBB1C2032ull, 27165893636ull,   687123146ull}, {0xC6D65BD8A6087F08ull, 15586314376ull, 54717373852ull} } },
	{ 15360, 282508657u, { {0xE7B08ED3A92EC6ECull,   875689313ull, 41754616020ull}, {0xD08FBAFF5CA5096Full, 30398073011ull, 62088094181ull}, {0xD6B7357DF761AA51ull, 28631146088ull, 26883666300ull} } },
	{ 16384, 300903377u, { {0xA23E8D2F532F05E6ull, 17871262795ull, 53388776441ull}, {0x14F20059083BF452ull, 16549596802ull, 56184170215ull}, {0x76B8A857EC9B3042ull, 14094306048ull, 61845793513ull} } },
	{ 18432, 337615277u, { {0xAEB976D153A4176Bull, 15040345558ull, 14542578090ull}, {0x503B443CB1E0CD2Dull, 29149628739ull,  5785599363ull}, {0x2D3047CEFF2F5A6Dull,  6100949709ull, 36303747216ull} } },
	{ 20480, 374233309u, { {0x6D95C0E62C8F9606ull,  3426866174ull, 39787406588ull}, {0xD08FB9031D460B7Eull, 30083048700ull, 30636357797ull}, {0x3A58018C387FBB68ull, 26771468430ull,  7763681227ull} } },
	{ 22528, 410766953u, { {0x254572CAB2014E6Cull, 17794160487ull, 13396317974ull}, {0x6202B11AA8602531ull,  9896689730ull, 13179771096ull}, {0x57919ED698DB2058ull, 12202175362ull, 54676834262ull} } },
	{ 24576, 447223969u, { {0x547DF462C7DAD1F6ull, 21523243373ull, 60663902864ull}, {0x06AB4D8E6FD9D69Bull, 23751736171ull, 10720018557ull}, {0x2B0DE10480A159B7ull,  8420637297ull, 56684836957ull} } },
	{ 26624, 483610763u, { {0xED3E248A29A1C6A8ull,  3863520901ull, 56560420765ull}, {0xC29358F8206746D6ull, 28829535828ull,  8160695393ull}, {0x56442E62439686ABull, 16425477242ull, 62275447263ull} } },
	{ 28672, 519932827u, { {0xCA7B3A76819D67F7ull,  7078504016ull, 32836389262ull}, {0x5799C8BE8E02B56Full, 22269194969ull, 11617462155ull}, {0xBA8FC230F7B5EA7Dull, 27830917747ull, 28996845727ull} } },
	{ 30720, 556194803u, { {0x56CDF80EFF67C1C1ull, 15611746359ull, 45837150207ull}, {0xCBC88456B4B47AC0ull,  3535843471ull, 18652930008ull}, {0xB0B6B164FC22EA57ull, 23758922431ull, 63785512003ull} } },
	{ 32768, 592400713u, { {0xD8C829884B234EB4ull, 13278383896ull, 54835531046ull}, {0x994DD6B24F452451ull, 28289805997ull, 11462134131ull}, {0x8EEAD14850955F52ull,  3931092242ull,  2483613485ull} } },
	{ 36864, 664658101u, { {0x2A809B0C735BAC4Bull, 10513414423ull, 54347266174ull}, {0xAB2147D9BAA22BB4ull, 12259954326ull, 67125404781ull}, {0xB4CFF625B3E3FD79ull, 11392807713ull, 32757679957ull} } },
	{ 40960, 736728527u, { {0xB9AC3EC848FF60A5ull,  7352037613ull,  7261166574ull}, {0x3D623A79D0F14EFFull, 31246254654ull, 49195074754ull}, {0x08F1C4771CDDC601ull, 26432814693ull, 42011833744ull} } },
	{ 45056, 808631029u, { {0x9D543D67F48AF766ull,  4288733345ull, 27338399854ull}, {0x62A4DF80612E897Bull, 32232536296ull, 47296762118ull}, {0xC42E0660237238BFull,  2574282651ull, 59023604757ull} } },
	{ 49152, 880380937u, { {0x82ED59E22D585BF6ull, 34028441473ull, 54282489944ull}, {0x7F110FD687DB7CB5ull, 14440703450ull, 57885720403ull}, {0xFB7D2A3EA41594FEull,  2563979873ull, 33101322377ull} } },
	{ 53248, 951990961u, { {0xE98B9D69E5F350D1ull, 14692034938ull,  2987127112ull}, {0x634157BCC596287Dull, 31799414844ull, 64590776355ull}, {0x2D82590AC2DCB435ull, 25060823443ull, 13978506048ull} } },
	{ 57344,1023472049u, { {0x960CE3D7029BDB70ull,  2075307892ull, 59259408155ull}, {0x9D98FF9A3FD21D1Bull, 10507186734ull,  3891581073ull}, {0x1CAE6ACE4720BCE8ull, 34199730716ull, 12202402908ull} } },
	{ 61440,1094833457u, { {0x9A2592E96BC1C827ull, 18678252759ull,   949397216ull}, {0x79BF2E2F5AE7985Bull,  8407199527ull, 64114889744ull}, {0x911DB4B5D8EFD861ull,  1735549534ull, 50988019040ull} } },
	/* Huge: */
	{ 65536,1166083297u, { {0xA0FE3066C834E360ull, 29128713056ull,  7270151463ull}, {0x4CC11F1C55F95C9Bull, 17910403663ull, 19262812656ull}, {                 0ull,           0ull,           0ull} } },
	{ 73728,1308275261u, { {0xE97D94366B318338ull, 27568210744ull, 15277294114ull}, {0x0913DB564B5FB2C9ull, 22759704201ull,  4713378634ull}, {                 0ull,           0ull,           0ull} } },
	{ 81920,1450094993u, { {0x0F0945ACF6C46D21ull, 55679675681ull,  4898080803ull}, {0x12FB0E6628819A5Aull, 26252443199ull, 65921528136ull}, {                 0ull,           0ull,           0ull} } },
	{ 90112,1591580099u, { {0xD8B487D8F5AFDF1Aull, 38481682202ull, 32203366427ull}, {0xAA5F7E4D5EC75AD3ull, 16596905178ull,  8253283767ull}, {                 0ull,           0ull,           0ull} } },
	{ 98304,1732761197u, { {0xB3ABD7E277F6EFB4ull, 10602606516ull, 18939646037ull}, {0x2170E14E6B46F009ull, 25931343295ull, 29383122206ull}, {                 0ull,           0ull,           0ull} } },
	{106496,1873663819u, { {0x938367CFDA83235Cull, 68090536796ull,   536356188ull}, {0xE4F8B81BD2E42112ull, 25813505565ull, 59441464991ull}, {                 0ull,           0ull,           0ull} } },
	{114688,2014309639u, { {0x84C8DED05749A21Full,  1464443423ull,  2648601266ull}, {0xDE2C03AA5068CC6Aull, 17222310589ull, 55630799101ull}, {                 0ull,           0ull,           0ull} } },
	{122880,2154717007u, { {0xE3056185380C64CFull, 22415172815ull, 18918434052ull}, {0x8DE279575CD1154Eull, 32123477176ull, 52211425638ull}, {                 0ull,           0ull,           0ull} } },
	{131072,2294901977u, { {0xB8CC9297E6E7512Cull, 33938690348ull, 32834437586ull}, {0x77FB54BF318C2DBAull, 20797246137ull, 29632517743ull}, {                 0ull,           0ull,           0ull} } },
	{147456,2574659081u, { {0x715208AF863EBA5Full, 66676767327ull, 17913296158ull}, {0xAA2E017A840AE336ull, 17798910323ull, 59557563397ull}, {                 0ull,           0ull,           0ull} } },
	{163840,2853674573u, { {0x2DF960079E128064ull, 32716783716ull,   830562867ull}, {0x341EB1C8403CE6AEull, 11934818175ull, 64851277709ull}, {                 0ull,           0ull,           0ull} } },
	{180224,3132023311u, { {0xD320ADB4B1483D7Eull, 20154170750ull,  9281699703ull}, {0x5008B456D49F551Full,  9130349955ull, 45637727974ull}, {                 0ull,           0ull,           0ull} } },
	{196608,3409766351u, { {0x1760B9AA1C6C9DFDull, 43426553341ull, 28002093744ull}, {0xE831AB643E67738Aull, 17588493473ull, 38449686407ull}, {                 0ull,           0ull,           0ull} } },
	{212992,3686954519u, { {0x9A3231FF531C8617ull, 10531609552ull, 42475524136ull}, {0x0CD00ADEE0D4D0DBull,  1995928646ull, 59987670753ull}, {                 0ull,           0ull,           0ull} } },
	{229376,3963630893u, { {0x427266E0A64E77F8ull,   135586632ull,  3161584476ull}, {0xBA1F2CE3D59018E7ull,  3171773156ull, 31951413199ull}, {                 0ull,           0ull,           0ull} } },
	{245760,4239832153u, { {0x05FF04CA8AEA8E54ull,   718789070ull, 67021029350ull}, {0xDCA27889AB5F4D92ull,    55229805ull, 33882881510ull}, {                 0ull,           0ull,           0ull} } },
/* Larger require 64-bit exponent support: */
	/* Egregious: */
	/* Brobdingnagian: */
	/* Godzillian: */
	{     0,         0u, { {                 0ull,           0ull,           0ull}, {                 0ull,           0ull,           0ull} } }
};

// Reference Mersenne base-3 PRP residues - we use the 0-padding slot in the above MersVec[]
// to store user-set-exponent data irrespective of whether LL-test or PRP-test, so no 0-pad here:
struct testMers MersVecPRP[numTest] =
{
/*                                         100-iteration residues:	                               1000-iteration residues:                */
/*	  FFTlen(K)     p              Res64           mod 2^35-1      mod 2^36-1               Res64           mod 2^35-1      mod 2^36-1     */
/*	    -----    --------     ----------------     -----------     -----------         ----------------     -----------     -----------    */
	/* Tiny:                                     [%34359738367  ][%68719476735  ]                         [%34359738367  ][%68719476735  ] */
	{     8,    173431u, { {0x59CB251F7D619429ull,  8164951424ull,  9846474576ull}, {0x918DD39863CB308Full,  7494791413ull, 15913265901ull}, {0x0ull, 0ull, 0ull} } },
	{     9,    194609u, { {0x6D47B1862C6AD488ull,  7469936757ull, 15426710219ull}, {0xE27427193254DABAull, 23325429059ull, 42526935522ull}, {0x0ull, 0ull, 0ull} } },
	{    10,    215767u, { {0xB7AA48826FA5622Eull, 18605212920ull, 51958758083ull}, {0x9F4A039DD6E9FAEEull, 19792956637ull, 29606652933ull}, {0x0ull, 0ull, 0ull} } },
	{    11,    236813u, { {0xE521AA8C68D3EA8Dull, 29439531719ull, 51299632013ull}, {0x2F13F1F265F2F2D3ull, 25527516254ull, 56556340711ull}, {0x0ull, 0ull, 0ull} } },
	{    12,    257903u, { {0x78D18DAA0211746Cull, 24776275699ull, 23960399868ull}, {0x5BA635C1FA4FC4C1ull, 16504851407ull, 15421897878ull}, {0x0ull, 0ull, 0ull} } },
	{    13,    278917u, { {0x466481AA8F523774ull,  4879692182ull, 35261457136ull}, {0x703FA3FE37EFFAA3ull, 17262836411ull, 57590389198ull}, {0x0ull, 0ull, 0ull} } },
	{    14,    299903u, { {0x2FFDF13628DD9271ull, 27842284148ull, 39770512868ull}, {0x8D10C61FEA78088Cull, 30967510385ull, 17203398344ull}, {0x0ull, 0ull, 0ull} } },
	{    15,    320851u, { {0xBFC67AB3CD20A5D1ull, 19065610760ull, 43681250762ull}, {0x3BC6C1D7281A842Full,  9947758605ull, 28952573372ull}, {0x0ull, 0ull, 0ull} } },
	{    16,    341749u, { {0x35EF82F656B51F1Full,  7029818683ull, 28919393248ull}, {0x4523FBE344DB54E1ull, 10285169590ull, 13266208357ull}, {0x0ull, 0ull, 0ull} } },
	{    18,    383521u, { {0x27490630F66623F0ull, 18175042063ull, 46462537136ull}, {0xED1B1867004B2043ull, 25094109590ull, 17095709313ull}, {0x0ull, 0ull, 0ull} } },
	{    20,    425149u, { {0x8D76E39499E2084Cull, 24854409897ull, 22759143673ull}, {0x93458CC621904526ull, 21512676920ull, 16851627129ull}, {0x0ull, 0ull, 0ull} } },
	{    22,    466733u, { {0x879E8E552CD9DA78ull, 23261845992ull, 68280331649ull}, {0x6236D9EC41B6B5B6ull,  4903056247ull, 48171921146ull}, {0x0ull, 0ull, 0ull} } },
	{    24,    508223u, { {0x60A1C833EFA3DB08ull, 20330461575ull, 33681508814ull}, {0x5478EE8AECEF292Dull, 14773583446ull, 53242362800ull}, {0x0ull, 0ull, 0ull} } },
	{    26,    549623u, { {0x7712E54BC57A6453ull, 22830350274ull,   975086245ull}, {0xCFCCD199B73D359Eull, 27986465261ull,  9872179127ull}, {0x0ull, 0ull, 0ull} } },
	{    28,    590963u, { {0x854A6C0BF1C0CA1Eull, 28465583444ull, 36207128387ull}, {0xD5494A4A30A6FF6Cull, 31443909484ull,  5696651641ull}, {0x0ull, 0ull, 0ull} } },
	{    30,    632251u, { {0xB9CE4D9CA145E37Eull,  2871789007ull, 25988330025ull}, {0x4FC369BD2786569Aull,  6701883912ull, 61126990820ull}, {0x0ull, 0ull, 0ull} } },
	{    32,    673469u, { {0x6FBFD30EA9AFA1DFull,  8567510574ull, 68209730811ull}, {0xDCAE0C41C5226178ull,  9502071634ull, 54856537405ull}, {0x0ull, 0ull, 0ull} } },
	{    36,    755737u, { {0x192B125399E579F8ull,   686716901ull, 14476488268ull}, {0x70C7BD4ECBFBBE50ull,  8165021045ull, 35292725948ull}, {0x0ull, 0ull, 0ull} } },
	{    40,    837817u, { {0x926C5D9C8E7DB701ull, 28915417288ull, 67081161725ull}, {0x90CF0AB873151579ull,  5754527985ull, 19567107501ull}, {0x0ull, 0ull, 0ull} } },
	{    44,    919729u, { {0x7E6995C2CCF4EF68ull, 28772226367ull, 18123291009ull}, {0x82AD751476F1686Aull, 10164437921ull,   713706080ull}, {0x0ull, 0ull, 0ull} } },
	{    48,   1001467u, { {0x42377D7C601947B7ull,  9525663359ull, 23726952031ull}, {0x8FC9D307DECC9D27ull,  9068970885ull, 10905178413ull}, {0x0ull, 0ull, 0ull} } },
	{    52,   1083077u, { {0x421844FFA21114FDull, 31935019641ull, 23327236236ull}, {0xEB311A383414F206ull, 32908672229ull, 27953087816ull}, {0x0ull, 0ull, 0ull} } },
	{    56,   1164533u, { {0xA53DFEE326A1047Dull, 16620085594ull,  8781848624ull}, {0xCADC99271CD41A6Full,  6810587643ull, 31020450647ull}, {0x0ull, 0ull, 0ull} } },
	{    60,   1245877u, { {0x375B40C9457701FCull,  2166978503ull, 39855448230ull}, {0x215D54D0A8FFFF60ull, 28855568754ull, 66990790203ull}, {0x0ull, 0ull, 0ull} } },
	{    64,   1327099u, { {0xE103AADBBB1C0647ull,  5046837446ull, 54383774262ull}, {0x2ECF95D374E31C8Bull, 12248705176ull, 62538924593ull}, {0x0ull, 0ull, 0ull} } },
	{    72,   1489223u, { {0xF6529A6352653F58ull, 13581202031ull,   526519996ull}, {0x889AB3D63246A31Eull, 13004174458ull, 26135310351ull}, {0x0ull, 0ull, 0ull} } },
	{    80,   1650959u, { {0x5FED95FC59F7565Cull, 15548760878ull, 61249993123ull}, {0xE4E8534D42053642ull,  7337120128ull, 20789742952ull}, {0x0ull, 0ull, 0ull} } },
	{    88,   1812347u, { {0x9AD66AD3088328A4ull, 12597406758ull, 63502374600ull}, {0xBF286CDFAD530AC0ull, 10289356709ull, 46915628424ull}, {0x0ull, 0ull, 0ull} } },
	{    96,   1973431u, { {0xAEE0C19E8742E3F4ull, 13726807213ull, 35724802848ull}, {0xBC213307E70191C5ull, 22512982343ull, 57049970336ull}, {0x0ull, 0ull, 0ull} } },
	{   104,   2134201u, { {0x098CAA3C0D107D6Bull, 22532136240ull, 26338802769ull}, {0x323774250D9B9516ull, 23273779605ull, 62209779235ull}, {0x0ull, 0ull, 0ull} } },
	{   112,   2294731u, { {0x47BC4AAA2FCCB1B9ull, 12907399344ull, 18863682101ull}, {0x91C50C372C12C703ull, 20908991540ull,  3451873052ull}, {0x0ull, 0ull, 0ull} } },
	{   120,   2455003u, { {0x63B512C24CB1233Bull,  5184435613ull, 27026274870ull}, {0x32E590C747238DA0ull, 21275876272ull, 57125434840ull}, {0x0ull, 0ull, 0ull} } },	/* Small: */
	{   128,   2614999u, { {0xBE4B6CE5EE76C4C0ull, 18066416827ull, 29451205104ull}, {0xB891A4CBE3EB0176ull,  2435134785ull, 53572085501ull}, {0x0ull, 0ull, 0ull} } },
	{   144,   2934479u, { {0x88A4403B9664A68Cull, 23316820590ull,  9544550461ull}, {0xFAD3E1AF502B1890ull, 32036343077ull, 67901574685ull}, {0x0ull, 0ull, 0ull} } },
	{   160,   3253153u, { {0x169BE706A0A0E185ull, 14483224776ull, 24483883970ull}, {0x5AED9CABD303E0A2ull, 31855596295ull, 44842948697ull}, {0x0ull, 0ull, 0ull} } },
	{   176,   3571153u, { {0x9C8A2EE1911F431Full,  1629566665ull, 66119740380ull}, {0x9C3CBB337EFFEC88ull,  5934946818ull,  8480402165ull}, {0x0ull, 0ull, 0ull} } },
	{   192,   3888509u, { {0x6CBAD54404C0AD1Aull, 12134100075ull, 18473223759ull}, {0x7A932AEA2AD84232ull,  6089810731ull, 44138905118ull}, {0x0ull, 0ull, 0ull} } },
	{   208,   4205303u, { {0x7D26A1F82679CECEull, 19082741636ull, 13798992673ull}, {0xEED98C7E1CC9ECC3ull, 13732641502ull, 49458162824ull}, {0x0ull, 0ull, 0ull} } },
	{   224,   4521557u, { {0x209349496454CEE6ull,  9940761503ull, 40783122082ull}, {0xE714F3E1C1F76042ull, 29521914431ull, 66307263755ull}, {0x0ull, 0ull, 0ull} } },
	{   240,   4837331u, { {0x4921BBB418D6C421ull, 32777132087ull, 40599659577ull}, {0xAF3009961984E8AFull,   901839485ull, 27501073879ull}, {0x0ull, 0ull, 0ull} } },
	{   256,   5152643u, { {0x2C67CB0282445407ull, 23339732910ull, 66204717452ull}, {0xC98A287287F9A1D7ull, 27688432919ull,  5877304726ull}, {0x0ull, 0ull, 0ull} } },
	{   288,   5782013u, { {0x985729C933D5E2FCull, 18320412728ull, 63235668703ull}, {0xAD12AB708012B28Eull, 24649900938ull, 62536412254ull}, {0x0ull, 0ull, 0ull} } },
	{   320,   6409849u, { {0xC58B1A6B98D84124ull, 18260049968ull,  6810024189ull}, {0x9B92CB21587AD6A6ull,  2110758259ull, 35105008257ull}, {0x0ull, 0ull, 0ull} } },
	{   352,   7036339u, { {0xC4364C1187F97A3Dull, 32150872197ull, 22781304616ull}, {0xA4D21311394D4FA7ull,  6324578481ull, 63084024685ull}, {0x0ull, 0ull, 0ull} } },
	{   384,   7661567u, { {0x1BBD63B803CE21D3ull, 16526828977ull, 42457373301ull}, {0xF4925EF17FC3E524ull, 20165927952ull,  9629593178ull}, {0x0ull, 0ull, 0ull} } },
	{   416,   8285659u, { {0xCA23151CCBE0740Full, 32337110525ull, 42837664157ull}, {0xA26DE72D0F7C6FD5ull, 11507652174ull, 49006121983ull}, {0x0ull, 0ull, 0ull} } },
	{   448,   8908723u, { {0xA95ECF765A7DE485ull, 26599751692ull, 19767035912ull}, {0x49222AF41653543Full, 25829109676ull, 51898058114ull}, {0x0ull, 0ull, 0ull} } },
	{   480,   9530803u, { {0x347EFC46B7C24927ull,  6753764868ull, 35076313672ull}, {0x4F47C7D80A2A809Eull, 15357957620ull, 19681757341ull}, {0x0ull, 0ull, 0ull} } },
	{   512,  10151971u, { {0xB3ED60F4AFBFB3A9ull, 22642977085ull,  1231176073ull}, {0x0EC12794332A0E0Eull,  6716458209ull, 51444232175ull}, {0x0ull, 0ull, 0ull} } },
	{   576,  11391823u, { {0x03ED6E37A156364Full,  1386523443ull, 30193596545ull}, {0xA4F7D77F194536D4ull, 33981920699ull, 60276616145ull}, {0x0ull, 0ull, 0ull} } },
	{   640,  12628613u, { {0xEAB53ECA446C7805ull,  3501890600ull, 56340973024ull}, {0x43BE1E7645CCBF2Cull, 32711506328ull, 41455841852ull}, {0x0ull, 0ull, 0ull} } },
	{   704,  13862759u, { {0x8D00941F37977416ull, 29194902133ull, 53354630705ull}, {0x9D080C44C635F531ull, 28299345740ull, 14583867538ull}, {0x0ull, 0ull, 0ull} } },
	{   768,  15094403u, { {0xC06D680B4315439Eull,   572672330ull, 17897298649ull}, {0xF5E8A6E2EB14969Dull, 26424320220ull, 33311538128ull}, {0x0ull, 0ull, 0ull} } },
	{   832,  16323773u, { {0x8A384DAF318BA11Cull, 29227782395ull, 49747663115ull}, {0x3A9A2BC7040F8F9Bull, 13221935219ull, 45353311168ull}, {0x0ull, 0ull, 0ull} } },
	{   896,  17551099u, { {0x236222CD1EF76404ull, 21554725722ull,  3020230570ull}, {0x15E08DB95124C451ull,  6571842111ull,   415617468ull}, {0x0ull, 0ull, 0ull} } },
	{   960,  18776473u, { {0x64DB852FB1D962E3ull, 30312455711ull,  3637871395ull}, {0xA654C7C6210CD8F2ull, 18496147373ull, 43979502788ull}, {0x0ull, 0ull, 0ull} } },
	{  1024,  20000047u, { {0x9D7E8C7336639EC4ull, 26990419379ull, 31343788296ull}, {0xEC57D79E86049AB4ull, 14020695804ull, 57562607914ull}, {0x5CC79106F08DC6D6ull, 14144638043ull,   651426155ull} } },
	{  1152,  22442237u, { {0xED82318E04A2A345ull,  8753176558ull, 32504299898ull}, {0x3B646DBA5F3653A2ull, 16747135485ull, 54376147696ull}, {0x2A8DCFFF3810D237ull,   597977664ull,  4018709428ull} } },
	{  1280,  24878401u, { {0xCE58A90EAC0139F1ull, 21686423079ull, 52735040981ull}, {0xE41022DE4EBE5A0Eull, 31073614876ull, 50104107291ull}, {0x71DAC0A6112FF806ull,  1019654400ull, 54335318694ull} } },
	{  1408,  27309229u, { {0xA070B7876AA20DA7ull, 26032112892ull, 32165824910ull}, {0x5B7DB0F2C5CA7EE0ull, 33557517402ull, 49136716314ull}, {0xC5285004A4A8738Eull, 23326684308ull, 15608345791ull} } },
	{  1536,  29735137u, { {0x70EF0D4A10AF820Cull, 31410231016ull, 19708126433ull}, {0xBA715ACC32497841ull, 30597754543ull, 12259094014ull}, {0x5505A14DFD239B10ull,  4793841111ull, 10366011218ull} } },
	{  1664,  32156581u, { {0x0ABAAF28B6995BF5ull, 25429439627ull, 29692330080ull}, {0xA3A513F3864F2486ull,  7487969602ull, 28002054176ull}, {0xBEAB8973C2FAAB17ull, 19372157564ull,  4721013073ull} } },
	{  1792,  34573867u, { {0x93B9E2F8EB75BFCDull, 28229043348ull, 23811859510ull}, {0xAF09EB79E9FD2368ull, 14328542178ull, 22451054390ull}, {0xDCC79070AEFDF160ull, 28503147261ull,  9196986298ull} } },
	{  1920,  36987271u, { {0x1FD9118DC82BEE78ull, 30917851138ull, 46406380675ull}, {0x33D7AD617976CEC1ull, 12975636226ull, 29187764572ull}, {0x5E856539865A176Cull,  9594612383ull, 46445499555ull} } },
	/* Medium: */
	{  2048,  39397201u, { {0x6D2DA75A04F538E6ull, 27183858872ull, 28453479981ull}, {0xF71AB763CBA394DDull, 20849861319ull, 65329085551ull}, {0xDCFC108BB85E2FA7ull, 24908241672ull, 53956544053ull} } },
	{  2304,  44207087u, { {0x17A51F8EBE222971ull, 20700958815ull, 19557281705ull}, {0x9603404F386381D1ull, 13688982479ull, 15674920724ull}, {0xC7B242D05373D3CAull,  2707265651ull,  5882282433ull} } },
	{  2560,  49005071u, { {0x778BE03E5722FE96ull, 25483813566ull, 14727845285ull}, {0x798CFC88734572AEull, 30521003775ull, 33673503490ull}, {0x9737799C41D18480ull, 18819745044ull, 37090058475ull} } },
	{  2816,  53792327u, { {0x7BE008D91E80DB49ull,  1791708263ull,  1480510204ull}, {0xDD858F590924665Full, 26857094761ull, 17153737658ull}, {0x8139A3391FC82B64ull, 14978416401ull, 23884206625ull} } },
	{  3072,  58569809u, { {0xECE8A3CA1F66D83Aull,  9805362291ull, 43153620706ull}, {0xFBA8C97760C5FE86ull, 11145859641ull, 16209757579ull}, {0x29FF471C9A28D07Full, 32759739090ull,  5331109071ull} } },
	{  3328,  63338459u, { {0x24A0F7BE95910DFCull,  6775505897ull, 39005906148ull}, {0x2B63FF56792E239Cull, 15844211278ull, 37865572416ull}, {0xDAEF98E50ECF3AABull, 28953421873ull, 49571783147ull} } },
	{  3584,  68098843u, { {0xAE2AED0851C344FFull, 29322480531ull,  7619686742ull}, {0x796289EFF35C1625ull, 20303644792ull, 16113781047ull}, {0xDEF7187C62A2FE6Eull,  1552555257ull, 23682222360ull} } },
	{  3840,  72851621u, { {0x41BB879F5BF89744ull, 15213631662ull, 47597325409ull}, {0x13AD29D142D2DD68ull,  2896274647ull, 28918042313ull}, {0x08BFA29C23275DD5ull, 15547700048ull, 36389342409ull} } },
	{  4096,  77597293u, { {0x5CA5EB4B2C403B40ull, 16067213695ull, 16060551561ull}, {0x8A91734601A3CA83ull, 25368547633ull, 15106685038ull}, {0x1B1628C39003FC08ull, 31104368473ull, 30825982076ull} } },
	{  4608,  87068977u, { {0x790629521C99A54Full,  6004740145ull, 58203809974ull}, {0x8C70DA87F106F6F6ull,  5407783662ull, 22056008955ull}, {0x45DC86154DE680CDull, 23780020949ull, 25229002925ull} } },
	{  5120,  96517019u, { {0x50BE76BDB9D43685ull, 33893245159ull, 32100820412ull}, {0x4D05925CD78F753Cull, 13802941362ull, 22609697002ull}, {0x954387C749FFA9F3ull, 14415933176ull, 19405933610ull} } },
	{  5632, 105943723u, { {0x1236DCECCCD91F1Bull, 11308642381ull, 28977161847ull}, {0x79A6575BCA014534ull, 30115274802ull,   828334544ull}, {0x4B8E27CD1CA2ED7Bull, 14367849395ull, 43143949281ull} } },
	{  6144, 115351063u, { {0x952594D45AC0CD66ull, 33154662725ull, 46573658072ull}, {0x71ABDF14CCB6F83Cull,  2372123324ull, 52765668028ull}, {0x82D8CD4BA692CB52ull, 14510733169ull, 67876507716ull} } },
	{  6656, 124740697u, { {0xF8E8C1410A022709ull, 13190590977ull, 39845130919ull}, {0x36F1A6E76E58E9C0ull, 27234521341ull, 16590439770ull}, {0xEE5448CFEDDBDEFEull,  2841484794ull, 12463515449ull} } },
	{  7168, 134113933u, { {0x1A946511EF44970Full, 21181711241ull, 30706184956ull}, {0x685E4EDB1A55944Aull,  6228216981ull, 24476904464ull}, {0xA96D930D64AA9158ull,  3314680740ull, 44392410575ull} } },
	{  7680, 143472073u, { {0xBD65A37723E2747Dull, 22515080503ull,  9322142584ull}, {0xB97E419D5BB91AE7ull,  4358792147ull, 10943432917ull}, {0xC16E2E327A370783ull, 31343185282ull, 21808110777ull} } },	/* Large: */
	{  8192, 152816047u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{  9216, 171465013u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 10240, 190066777u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 11264, 208626181u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 12288, 227147083u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 13312, 245632679u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 14336, 264085733u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 15360, 282508657u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 16384, 300903377u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 18432, 337615277u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 20480, 374233309u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 22528, 410766953u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 24576, 447223969u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 26624, 483610763u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 28672, 519932827u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 30720, 556194803u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 32768, 592400713u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 36864, 664658101u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 40960, 736728527u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 45056, 808631029u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 49152, 880380937u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 53248, 951990961u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 57344,1023472049u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	{ 61440,1094833457u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull} } },
	/* Huge: */
	{ 65536,1166083297u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{ 73728,1308275261u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{ 81920,1450094993u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{ 90112,1591580099u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{ 98304,1732761197u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{106496,1873663819u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{114688,2014309639u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{122880,2154717007u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{131072,2294901977u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{147456,2574659081u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{163840,2853674573u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{180224,3132023311u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{196608,3409766351u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{212992,3686954519u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{229376,3963630893u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } },
	{245760,4239832153u, { {0x0ull, 0ull, 0ull}, {0x0ull, 0ull, 0ull}, {                 0ull,           0ull,           0ull} } }
/* Larger require 64-bit exponent support: */
	/* Egregious: */
	/* Brobdingnagian: */
	/* Godzillian: */
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
	int Fidx;			/* Fermat number index */
	struct res_triplet	res_t[3];	/* 100,1000 and 10000-iteration SH-residue triplets */
};

/* Array of 100-iteration reference values for Fermat self-tests. Only allow Fidx >= 14: */
#define FermArrayIdxOffset		14	/* Amount to subtract from m to get the prestored residues for F_m */
#define numFerm		20
struct testMers FermVec[numFerm+1] =
{
/*                                   100-iteration residues:                                  1000-iteration residues:              */
/* FFTlen(K) Fidx           Res64           mod 2^35-1      mod 2^36-1               Res64           mod 2^35-1      mod 2^36-1     */
/*   ------- ----      ----------------     -----------     -----------         ----------------     -----------     -----------    */
	/*                                    [%34359738367  ][%68719476735  ]                         [%34359738367  ][%68719476735  ] */
	{     1,  14, { {0xDB9AC520C403CB21ull,   342168579ull, 59244817440ull}, {0xF111F12732CCCB0Full, 24848612524ull, 66609820796ull}, {0x78738D068D641C2Cull, 12664273769ull, 29297626750ull} } },
	{     2,  15, { {0x3B21A6E55ED13454ull, 28379302213ull, 15546218647ull}, {0x4784657F2A36BE74ull,   617376037ull, 44891093359ull}, {0x589BFE53458FFC14ull, 12200390606ull, 46971422957ull} } },
	{     4,  16, { {0xAAE76C15C2B37465ull, 20013824731ull,  2261076122ull}, {0x42CC2CBE97C728E6ull, 30814966349ull, 44505312792ull}, {0xEED00D8AE6886440ull, 19057922301ull, 53800020279ull} } },
	{     8,  17, { {0xFFA16CDC8C87483Cull, 20917337408ull, 26110327818ull}, {0x43CAB295FFB2661Full, 18197605796ull,  9842643677ull}, {0x5C7B0049549D4174ull,  8923959253ull, 40303785249ull} } },
	{    16,  18, { {0x7C6B681485EB86DBull,  5745147782ull, 50521157289ull}, {0x8193BD41931E9DE8ull, 19662968587ull, 51102742548ull}, {0x9D6A242467D28700ull, 12912307491ull, 23293425575ull} } },
	{    32,  19, { {0x529E54642A813995ull, 17797950508ull, 32039741221ull}, {0xE24EAE4B153EE86Bull, 11155350666ull, 49866866361ull}, {0x787FABD98DD5FEC0ull, 10784283812ull, 50254721650ull} } },
	{    64,  20, { {0x64629CED6E218018ull,  8485981669ull, 53977437340ull}, {0xA380121F6FD26B2Aull, 15876203498ull, 36314727556ull}, {0x352CB92DABA82A8Bull,  6625203514ull, 20044302250ull} } },
	{   128,  21, { {0x1DE0171591038250ull, 33758422990ull,  8269940507ull}, {0x1B331FBB41AF33D7ull, 17971032338ull,  2929392342ull}, {0x2E7F9D30EAFC7D47ull, 28599205675ull, 44913594527ull} } },
	{   256,  22, { {0x9201143390F3828Dull,  1749100092ull, 46602233256ull}, {0x1B331FBB41AF33D7ull, 17971032338ull,  2929392342ull}, {0x4380A384A80C079Eull, 30799637462ull, 32805936529ull} } },
	{   512,  23, { {0x9C3F8E29B397B32Bull,  1094055486ull, 13316822657ull}, {0xBD642EA0479D8FF0ull, 31625967305ull, 57187857233ull}, {0xC0F18B226C20AE50ull, 31325111828ull, 10029614864ull} } },
	{  1024,  24, { {0xDB9F01963ED9DC8Bull, 27887793041ull, 13169874547ull}, {0x40F2DECE9C351236ull,  9074286032ull, 38590877049ull}, {0x84FEE61771F20003ull, 19739002855ull, 42109658313ull} } },
	{  2048,  25, { {0x376C33921E5F675Full, 13022327996ull, 46818697393ull}, {0xA51F8577A407CB75ull,  9865976783ull, 35171498411ull}, {0x2B1D829E05C85CA3ull, 27914095914ull, 58025822027ull} } },
	{  4096,  26, { {0xA42BECD80DAEC4CBull, 10087739060ull, 25252768685ull}, {0xECC9408A7295401Dull,  5904751941ull, 58967745948ull}, {0x0C054D0376BD8E9Eull, 24761357761ull, 23215069147ull} } },
	{  8192,  27, { {0xFB69E377519D8CE6ull, 15449775614ull, 51221672039ull}, {0x24898E3BEB59DCE6ull, 24957168001ull,  2072452827ull}, {0x443DE289690070EDull, 30674247970ull, 29955017069ull} } },
	{ 16384,  28, { {0xA4FF6F8C3CB38B85ull, 18933356966ull, 30899345457ull}, {0x8B451AF25E8CC50Eull,   674652743ull, 39963850167ull}, {0x0F38D26B35C6794Cull,  7003106751ull, 60469235270ull} } },
	{ 32768,  29, { {0xAFBF110B593E26F6ull, 32666279868ull, 18995112582ull}, {0xA6B643FF24C6ADC1ull, 15753158767ull, 13965270144ull}, {0xB194096866F68C59ull, 19667273394ull,  3552165634ull} } },
	{ 65536,  30, { {0x68B1BDA5D6BAE04Bull,  3347054148ull, 47892955488ull}, {0x361F30024AF9FE26ull, 26693502373ull, 67933083515ull}, {0x3A70D98DA3ED809Aull, 10877397803ull, 41746600776ull} } },
	{131072,  31, { {0x00C36B4AB38FC326ull, 12487821860ull, 64847210796ull}, {0x243D46185B4FAAC8ull, 14946978930ull, 42870390813ull}, {0xDD3651B25892F424ull,  6237393494ull, 22880395670ull} } },
	{262144,  32, { {0xFB5A3146BE2CA886ull, 33629944997ull, 19559359478ull}, {0x0479E68514EAD529ull,    59701496ull, 58397447084ull}, {0x78CD466FDDA4F156ull,  5752278548ull, 18141301956ull} } },
	{524288,  33, { {0x25D9822BD1555EB9ull, 10052071893ull, 20645266428ull}, {0xACF521F811129C9Eull, 25612341094ull, 61523429244ull}, {0x269405A9B18B3310ull,  4910011908ull,  2895484666ull} } },
	{     0,   0, { {                 0ull,           0ull,           0ull}, {                 0ull,           0ull,           0ull}, {                 0ull,           0ull,           0ull} } }
};

/***************************************************************************************/

int 	main(int argc, char *argv[])
{
	int		retVal=0;
	uint64	Res64, Res35m1, Res36m1;
	char	stFlag[STR_MAX_LEN];
	uint64	i64arg;
	uint32	iarg = 0, iters = 0, k = 0, maxFFT, expo = 0, findex = 0;
	double	darg;
	int		new_cfg = FALSE;
	int		i,j, idum, lo, hi, nargs, scrnFlag;
	int		start = -1, finish = -1, modType = 0, testType = 0, selfTest = 0, userSetExponent = 0, xNum = 0;
#ifdef MULTITHREAD
	// Vars for mgmt of mutually exclusive arg sets:
	int		nthread = 0, cpu = 0;
#endif
	int		quick_self_test = 0, fftlen = 0, radset = -1;
	double	runtime, runtime_best, tdiff;
	double	roerr_avg, roerr_max;
	int		radix_set, radix_best, nradix_set_succeed;

	uint32 mvec_res_t_idx = 0;	/* Lookup index into the res_triplet table */
	uint32 new_data;
	struct res_triplet new_res = {0ull,0ull,0ull};

/* Enable this and set upper loop bound appropriately to regenerate a quick list of PRPs
just below the upper limit for each FFT lengh in some subrange of the self-tests:
*/
#if 0
	for(i = 68; i < numTest; i++)
	{
		j = given_N_get_maxP(MersVec[i].fftLength << 10);
		if((j&1) == 0) --j; 	// make sure it's odd
		for(;;)
		{
			if(isPRP(j))
			{
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

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		ASSERT(HERE, (MersVec[numTest-1].fftLength != 0) &&  (MersVec[numTest].fftLength == 0), "numTest != MersVec allocated size!");
	} else {
		ASSERT(HERE, (FermVec[numFerm-1].fftLength != 0) &&  (FermVec[numFerm].fftLength == 0), "numFerm != FermVec allocated size!");
	}

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
		strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

		if(stFlag[0] != '-')
		{
			fprintf(stderr, "*** ERROR: Illegal command-line option %s\n", stFlag);
			fprintf(stderr, "*** All command-line options must be of form -{flag} [argument]\n\n");
			print_help("");
		}

		if(STREQ(stFlag, "-h"))
		{
			print_help("");
		}

		if(STREQ(stFlag, "-topic"))
		{
			print_help("topic");
		}

		/* Mersenne self-test: requires a user-set exponent, FFT length or one of the supported -s arguments below: */
		if(STREQ(stFlag, "-s"))
		{
			selfTest = TRUE;

			if(nargs >= argc) {
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("self_test");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			for(;;)
			{
				if(STREQ(stFlag, "a") || STREQ(stFlag, "all"))		/* all, which really means all the non-Huge-and-larger sets */
				{
					start = 0; finish = numTiny + numSmall + numMedium + numLarge;
					break;
				}

				finish = 0;

				start = finish; finish += numTiny;
				if(STREQ(stFlag, "t") || STREQ(stFlag, "tiny"))		/* tiny   */
				{
					break;
				}
				start = finish; finish += numSmall;
				if(STREQ(stFlag, "s") || STREQ(stFlag, "small"))	/* small  */
				{
					break;
				}
				start = finish; finish += numMedium;
				if(STREQ(stFlag, "m") || STREQ(stFlag, "medium"))	/* medium */
				{
					break;
				}
				start = finish; finish += numLarge;
				if(STREQ(stFlag, "l") || STREQ(stFlag, "large"))	/* large  */
				{
					break;
				}
				start = finish; finish += numHuge;
				if(STREQ(stFlag, "h") || STREQ(stFlag, "xl"))		/* huge   */
				{
					break;
				}
				start = finish; finish += numEgregious;
				if(STREQ(stFlag, "e") || STREQ(stFlag, "xxl"))
				{
					break;
				}
				start = finish; finish += numBrobdingnagian;
				if(STREQ(stFlag, "b") || STREQ(stFlag, "xxxl"))
				{
					break;
				}
				start = finish; finish += numGodzillian;
				if(STREQ(stFlag, "g") || STREQ(stFlag, "gojira"))
				{
					break;
				}

				fprintf(stderr, "*** ERROR: Illegal argument %s to -s flag\n", stFlag);
				print_help("self_test");
			}

			modType  = MODULUS_TYPE_MERSENNE;
		}

		else if(STREQ(stFlag, "-iters"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("iters");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			/* Convert the #iters argument to an int: */
			iarg = 0;
			darg = 0.0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					darg = 10.0*darg + (stFlag[i]-CHAROFFSET);
					iarg = (uint32)darg;
					/* Check for overflow: */
					if((double)iarg != darg)
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -iters argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** ERROR: Non-numeric character encountered in -iters argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}

			iters = (uint32)iarg;
		}

		else if(STREQ(stFlag, "-fftlen"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("fftlen");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			/* Convert the length argument to an int: */
			iarg = 0;
			darg = 0.0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					darg = 10.0*darg + (stFlag[i]-CHAROFFSET);
					iarg = (uint32)darg;
					/* Check for overflow: */
					if((double)iarg != darg)
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -fftlen argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** ERROR: Non-numeric character encountered in -fftlen argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}

			fftlen = (uint32)iarg;
			if((i = get_fft_radices(fftlen, 0, 0x0, 0x0, 0)) != 0)
			{
				sprintf(cbuf  , "ERROR: FFT length %d K not available.\n",fftlen);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
		}

		else if(STREQ(stFlag, "-radset"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("radset");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			/* Convert the length argument to an int: */
			iarg = 0;
			darg = 0.0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					darg = 10.0*darg + (stFlag[i]-CHAROFFSET);
					iarg = (uint32)darg;
					/* Check for overflow: */
					if((double)iarg != darg)
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -radset argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** ERROR: Non-numeric character encountered in -radset argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}
			radset = iarg;
		}

		else if(STREQ(stFlag, "-shift"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("shift");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			/* Convert the shift argument to a uint64: */
			i64arg = 0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					i64arg = 10*i64arg + (stFlag[i]-CHAROFFSET);
					/* Check for overflow: */
					if(i64arg % (uint64)10 != (uint64)(stFlag[i]-CHAROFFSET))
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -shift argument %s overflows uint64 field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** ERROR: Non-numeric character encountered in -shift argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}
			RES_SHIFT = i64arg;
		#ifndef USE_SSE2	// Only support residue-shift for SIMD builds:
			if(RES_SHIFT) {
				fprintf(stderr, "*** INFO: Only support residue-shift for SIMD builds ... ignoring user-set value.\n");
				fprintf(stderr,"%s", cbuf);
				RES_SHIFT = 0;
			}
		#endif	// USE_SSE2	?
		}

		else if(STREQ(stFlag, "-prp"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("prp");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			/* Convert the initial-seed argument to an int: */
			iarg = 0;
			darg = 0.0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					darg = 10.0*darg + (stFlag[i]-CHAROFFSET);
					iarg = (uint32)darg;
					/* Check for overflow: */
					if((double)iarg != darg)
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -prp argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** ERROR: Non-numeric character encountered in -prp argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}
			PRP_BASE = iarg;
			testType = TEST_TYPE_PRP;
		}

		else if(STREQ(stFlag, "-nthread"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("nthread");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			/* Convert the length argument to an int: */
			iarg = 0;
			darg = 0.0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					darg = 10.0*darg + (stFlag[i]-CHAROFFSET);
					iarg = (uint32)darg;
					/* Check for overflow: */
					if((double)iarg != darg)
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -nthread argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** ERROR: Non-numeric character encountered in -nthread argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}

		#ifndef MULTITHREAD
			ASSERT(HERE,0,"Multithreading must be enabled in build to permit -nthread argument!");
		#else
			NTHREADS = iarg;
			ASSERT(HERE,cpu == FALSE,"Only one of -nthread and -cpu flags permitted!");
			nthread = TRUE;
			// Use the same affinity-setting code here as for the -cpu option, but simply for cores [0:NTHREADS-1]:
			sprintf(cbuf,"0:%d",NTHREADS-1);
			parseAffinityString(cbuf);
		#endif
		}

		else if(STREQ(stFlag, "-cpu"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("nthread");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

		#ifndef MULTITHREAD
			ASSERT(HERE,0,"Multithreading must be enabled in build to permit -cpu argument!");
		#else
			parseAffinityString(stFlag);
			ASSERT(HERE,nthread == FALSE,"Only one of -nthread and -cpu flags permitted!");
			cpu = TRUE;
		#endif
		}

		/******************************************************************************************/
		/* MAKE SURE ALL OTHER FLAG-PROCESSING SECTIONS SET userSetExponent TO A NONZERO VALUE!!! */
		/******************************************************************************************/

		else if(STREQ(stFlag, "-m") || STREQ(stFlag, "-mersenne"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("mersenne");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			/* Convert the exponent argument to an int: */
			iarg = 0;
			darg = 0.0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					darg = 10.0*darg + (stFlag[i]-CHAROFFSET);
					iarg = (uint32)darg;
					/* Check for overflow: */
					if((double)iarg != darg)
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -m argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** INFO: Non-numeric character encountered in -m argument %s ... using default exponent for whatever FFT length is supplied.\n", stFlag);
					fprintf(stderr,"%s", cbuf);
					--nargs;
					goto SET_MERS;
				}
			}

			expo = (uint32)iarg;
			userSetExponent = 1;
			//*** use 0-pad slot in MersVec[] to store user-set-exponent data irrespective of whether LL-test or PRP-test: ***
			MersVec[numTest].exponent = expo;
			start = numTest; finish = start+1;
		SET_MERS:
			modType = MODULUS_TYPE_MERSENNE;
		}

		else if(STREQ(stFlag, "-f") || STREQ(stFlag, "-fermat"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				print_help("fermat");
			}

			strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

			/* Convert the exponent argument to an int: */
			iarg = 0;
			darg = 0.0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					darg = 10.0*darg + (stFlag[i]-CHAROFFSET);
					iarg = (uint32)darg;
					/* Check for overflow: */
					if((double)iarg != darg)
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -f argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** INFO: Non-numeric character encountered in -f argument %s ... using default exponent for whatever FFT length is supplied.\n", stFlag);
					fprintf(stderr,"%s", cbuf);
					--nargs;
					goto SET_FERM;
				}
			}

			/* Make sure the Fermat number index is in range: */
			if(iarg < 14)
			{
				fprintf(stderr, " Fermat number index must be at least 14.\n");
				return ERR_EXPONENT_ILLEGAL;
			}
			if(iarg > 63)
			{
				fprintf(stderr, " Fermat number index must be < 64.\n");
				return ERR_EXPONENT_ILLEGAL;
			}

			findex = iarg;
			userSetExponent = 1;
			FermVec[numFerm].exponent = findex;
			start = numFerm; finish = start+1;
		SET_FERM:
			modType = MODULUS_TYPE_FERMAT;
		}

		else
		{
			fprintf(stderr, "*** ERROR: Unrecognized flag %s.\n", stFlag);
			print_help("");
		}
	}	/* end of command-line-argument processing while() loop */

	if(!modType)
	{
		modType = MODULUS_TYPE_MERSENNE;
	}
	// Now that have determined the modType, copy any user-set FFT length into the appropriate field:
	if(fftlen) {
		/* Don't set userSetExponent here, since -fftlen can be invoked without an explicit exponent */
		if(modType == MODULUS_TYPE_FERMAT)
		{
			FermVec[numFerm].fftLength = fftlen;
			start = numFerm; finish = start+1;
		}
		else
		{
			MersVec[numTest].fftLength = fftlen;
			start = numTest; finish = start+1;
		}
	}
	// If user has specified a radix set, make sure an FFT length has also specified:
	if(radset != -1) {
		if(modType == MODULUS_TYPE_FERMAT)
		{
			iarg = FermVec[numFerm].fftLength;
		}
		else
		{
			iarg = MersVec[numTest].fftLength;
		}
		if(iarg == 0)
		{
			sprintf(cbuf  , "*** ERROR: Must specify a valid FFT length on command line before -radset argument!\n");
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		/* Make sure it's a valid radix set index for this FFT length: */
		if((i = get_fft_radices(iarg, radset, &idum, 0x0, 0)) != 0)
		{
			if     (i == ERR_FFTLENGTH_ILLEGAL)
				sprintf(cbuf  , "ERROR: FFT length %d K illegal!\n", iarg);
			else if(i == ERR_RADIXSET_UNAVAILABLE)
				sprintf(cbuf  , "ERROR: radix set index %d for FFT length %d K exceeds maximum allowable of %d.\n",radset, iarg, idum-1);
			else
				sprintf(cbuf  , "ERROR: Unknown error-code value %d from get_fft_radices(), called with radix set index %d, FFT length %d K\n",i,radset, iarg);

			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

	}

	// Use selfTest == TRUE or -iters (in the single-FFT-length-timing case) to differentiate between timing-test and production runs.
	// In fact it eases the logic to explicitly set selfTest = TRUE whenever iters is set (but not nec. the converse), so do that here:
	if(iters) selfTest = TRUE;

	if(modType == MODULUS_TYPE_MERSENNE && !selfTest)
	{
		if(userSetExponent) {
			ASSERT(HERE, start > 0, "userSetExponent = TRUE but self-test starting-index unset!");
			sprintf(cbuf, "ERROR: Production-run-mode [-iters not invoked] does not allow command-line\nsetting of exponent - that must be read from the worktodo.ini file.\n");
			ASSERT(HERE, 0,cbuf);
		} else if(start == -1) {
			start = numTest; finish = start+1;
		}
		if(radset != -1) {
			sprintf(cbuf, "ERROR: Production-run-mode [-iters not invoked] allows command-line setting of\nFFT length, but not the radix set - that must be read from the mlucas.cfg file.\n");
			ASSERT(HERE, 0,cbuf);
		}
	ERNST_MAIN:
		if((retVal = ernstMain(modType,testType,0,MersVec[start].fftLength,0,0,0,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime)) != 0)
		{
			printMlucasErrCode(retVal);

			/* If need to run a timing self-test at a particular FFT length, do that and then try again... */
			if((retVal & 0xff) == ERR_RUN_SELFTEST_FORLENGTH)
			{
				quick_self_test = TRUE;
				selfTest = TRUE;
				k = (uint32)(retVal >> 8);
				if((i = get_fft_radices(k, 0, 0x0, 0x0, 0)) != 0)
				{
					sprintf(cbuf, "ERROR: FFT length %d K not available.\n",k);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}

			/**** IF POSSIBLE, USE ONE OF THE STANDARD TEST EXPONENTS HERE, SO CAN CHECK RES64s!!! ****/
				for(i = 0; i < numTest; i++)
				{
					if(MersVec[i].fftLength == k)
					{
						userSetExponent = 0;
						start = i; finish = start+1;
						break;
					}
				}
				if(i == numTest)
				{
					userSetExponent = 1;
					MersVec[numTest].exponent = convert_base10_char_uint64(ESTRING);
					MersVec[numTest].fftLength = k;
					start = numTest; finish = start+1;
				}

				modType = MODULUS_TYPE_MERSENNE;
				goto TIMING_TEST_LOOP;
			}
			/* ...Otherwise barf. */
			else
			{
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
		}
		else
		{
			fprintf(stderr, "\n  Done ...\n\n");
			exit(EXIT_SUCCESS);
		}
	}
	/* If user specified iters and FFT length but no exponent, get default Mersenne exponent for that FFT length: */
	else if(modType == MODULUS_TYPE_MERSENNE)
	{
		if(MersVec[start].exponent == 0)
		{
			i = MersVec[start].fftLength;
			ASSERT(HERE, i > 0                  ,"Require i > 0                  ");
			ASSERT(HERE, i <=MAX_FFT_LENGTH_IN_K,"Require i <=MAX_FFT_LENGTH_IN_K");

			/* If the FFT length is not represented in MersVec[],
			find the nearest prime <= given_N_get_maxP(FFT length):
			*/
			for(j = 0; j < numTest; j++) {
				if(i == MersVec[j].fftLength) break;
			}
			if(i != MersVec[j].fftLength)
			{
				hi = given_N_get_maxP(i<<10) | 0x1;	/* Make sure starting value is odd */
				lo = hi - 1000;	if(lo < PMIN) lo = PMIN;
				for(i = hi; i >=lo; i -= 2)
				{
					if(isPRP(i))
					{
						MersVec[numTest].exponent = i;
						break;
					}
				}
				if(i < lo || lo >= hi)
				{
					fprintf(stderr, "ERROR: unable to find a prime in the interval %u <= x <= %u.\n", lo, hi);
					ASSERT(HERE, 0,"0");
				}
			}
			else	/* Use the corresponding entry of MersVec: */
			{
				start = j; finish = start+1;
			}
		}
		/* If user specified exponent but no FFT length, get default FFT length for that exponent: */
		else if(MersVec[numTest].exponent && (MersVec[numTest].fftLength == 0))
		{
			MersVec[numTest].fftLength = get_default_fft_length((uint64)(MersVec[numTest].exponent));
		}
	}
	else if(modType == MODULUS_TYPE_FERMAT)
	{
		if(FermVec[start].exponent == 0)
		{
			i = FermVec[start].fftLength;
			ASSERT(HERE, i > 0                  ,"Require i > 0                  ");
			ASSERT(HERE, i <=MAX_FFT_LENGTH_IN_K,"Require i <=MAX_FFT_LENGTH_IN_K");

			if(i > FermVec[numFerm-1].fftLength)	/* Computing a new-largest entry? */
			{
				FermVec[numFerm].exponent = (i << 4);
			}
			else	/* Find the corresponding entry of FermVec: */
			{
				for(lo = 0; lo < numFerm; lo++)
				{
					if(FermVec[lo].fftLength >= i)
					{
						start = lo; finish = start+1;	/* Using >= here allows for non-power-of-2 FFT lengths */
						break;
					}
				}
				if(lo >= numFerm)
				{
					fprintf(stderr, "ERROR: unable to find FFT length %d K in the Reference Residue table.\n", i);
					ASSERT(HERE, 0,"0");
				}
			}
		}
		/* If user specified exponent but no FFT length, get default FFT length for that exponent: */
		else if(findex && (FermVec[numFerm].fftLength == 0))
		{
			FermVec[numFerm].fftLength = get_default_fft_length((uint64)1 << findex);
		}
	}
	else
	{
		ASSERT(HERE, 0,"modType not recognized!");
	}

TIMING_TEST_LOOP:

	if(selfTest) {
		fprintf(stderr, "\n           Mlucas selftest running.....\n\n");
		/* We have precomputed 100, 1000 and 10000-iteration residues for the predefined self-test exponents: */
		if( (userSetExponent && (modType == MODULUS_TYPE_MERSENNE)) || (iters && (iters != 100) && (iters != 1000) && (iters != 10000)) )
		{
			fprintf(stderr, "\n********** Non-default exponent or #iters not one of [100,1000,10000] - You will have to manually  **********");
			fprintf(stderr, "\n********** verify that the Res64s output for this self-test match for all FFT radix combinations!! **********\n\n\n");
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
	if(!fp)
	{
		sprintf(cbuf, "INFO: Unable to find/open %s file in %s mode ... creating from scratch.\n", CONFIGFILE, FILE_ACCESS_MODE);
		fprintf(stderr,"%s",cbuf);
		new_cfg = TRUE;
	}
	else
	{
		/* Program version string assumed to be stored in line 1 of .cfg file -
		if it's not, the .cfg file is by definition outdated:
		*/
		if(fgets(in_line, STR_MAX_LEN, fp))
			new_cfg = cfgNeedsUpdating(in_line);
		else
			new_cfg = TRUE;

		fclose(fp);	fp = 0x0;
	}

	/* If existing .cfg file is for current program version, file-append mode is appropriate;
	otherwise set up to re-open in write mode, which will clear any previous file contents:
	*/
	if(new_cfg)
	{
		FILE_ACCESS_MODE[0] = FILE_ACCESS_WRITE;
	}
	else
	{
		FILE_ACCESS_MODE[0] = FILE_ACCESS_APPEND;
	}
	FILE_ACCESS_MODE[1] = FILE_FORMAT_ASCII;

	/* What's the max. FFT length (in K) for the set of self-tests? */
	maxFFT = MersVec[finish-1].fftLength;

	for (xNum = start; xNum < finish; xNum++)    /* Step through the exponents */
	{
		new_data = FALSE;	Res64 = Res36m1 = Res35m1 = 0ull;

		/* If it's a self-test [i.e. timing test] and user hasn't specified #iters, set to default: */
		if(selfTest && !iters) {
			if(NTHREADS >= 4)
				iters = 1000;
			else
				iters = 100;
		}

		if(iters == 100 || iters == 1000 || iters == 10000) {
			mvec_res_t_idx = NINT( log((double)iters)/log(10.) ) - 2;	/* log10(iters) - 2, use slower NINT rather than DNINT here since latter needs correct rounding mode */
			ASSERT(HERE, mvec_res_t_idx < 3,"main: mvec_res_t_idx out of range!");
			// Since use empty-data-slot at top of MersVec[] for both primality & prp single-case tests that,
			// rather than MersVecPRP[], appears here in the PRP-part of the conditional:
			if( (modType == MODULUS_TYPE_MERSENNE && testType == TEST_TYPE_PRIMALITY && MersVec[xNum].res_t[mvec_res_t_idx].sh0 == 0)
			 || (modType == MODULUS_TYPE_MERSENNE && testType == TEST_TYPE_PRP       && MersVec[xNum].res_t[mvec_res_t_idx].sh0 == 0)
			 || (modType == MODULUS_TYPE_FERMAT   && FermVec[xNum].res_t[mvec_res_t_idx].sh0 == 0) )	/* New self-test residue being computed */
			{
				new_data = TRUE;
				new_res.sh0 = new_res.sh1 = new_res.sh2 = 0ull;
			}
		}

		/* Init best-radix-set, best-runtime and #radix-sets-which-succeeded for this FFT length: */
		runtime_best = 0.0;
		radix_best = -1;
		nradix_set_succeed = 0;

		/* If user-specified radix set, do only that one: */
		if(radset >= 0)
			radix_set = radset;
		else
			radix_set = 0;

		if(modType == MODULUS_TYPE_MERSENNE)
			iarg = MersVec[xNum].fftLength;
		else if(modType == MODULUS_TYPE_FERMAT)
			iarg = FermVec[xNum].fftLength;

		while(get_fft_radices(iarg, radix_set, &NRADICES, RADIX_VEC, 10) == 0)	/* Try all the radix sets available for this FFT length. */
		{
			if(modType == MODULUS_TYPE_FERMAT)
			{
				Res64   = FermVec[xNum].res_t[mvec_res_t_idx].sh0;
				Res35m1 = FermVec[xNum].res_t[mvec_res_t_idx].sh1;
				Res36m1 = FermVec[xNum].res_t[mvec_res_t_idx].sh2;
				retVal = ernstMain(modType,testType,(uint64)FermVec[xNum].exponent,iarg,radix_set,maxFFT,iters,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime);
			}
			else if(modType == MODULUS_TYPE_MERSENNE && testType == TEST_TYPE_PRIMALITY)
			{
				Res64   = MersVec[xNum].res_t[mvec_res_t_idx].sh0;
				Res35m1 = MersVec[xNum].res_t[mvec_res_t_idx].sh1;
				Res36m1 = MersVec[xNum].res_t[mvec_res_t_idx].sh2;
				retVal = ernstMain(modType,testType,(uint64)MersVec[xNum].exponent,iarg,radix_set,maxFFT,iters,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime);
			}
			else if(modType == MODULUS_TYPE_MERSENNE && testType == TEST_TYPE_PRP)
			{
				Res64   = MersVecPRP[xNum].res_t[mvec_res_t_idx].sh0;
				Res35m1 = MersVecPRP[xNum].res_t[mvec_res_t_idx].sh1;
				Res36m1 = MersVecPRP[xNum].res_t[mvec_res_t_idx].sh2;
				retVal = ernstMain(modType,testType,(uint64)MersVec[xNum].exponent,iarg,radix_set,maxFFT,iters,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime);
			}
			else
				ASSERT(HERE, 0,"Unsupported modulus and/or test type!");

			// If retVal != 0 relates to dangerously high ROEs, use (to-do: factor in #occurrences) maxErr to decide whether to accept the radix set.
			// AME,MME contain avgMaxErr and maxMaxErr for iter-interval, ROE_ITER holds #maxErr > 0.40625;	ROE_VAL recapitulates MME
			// (meaning we can use it to store some additional ROE-related measure should that become desirable in the future):
			if(selfTest && ( !userSetExponent && ((iters == 100) || (iters == 1000) || (iters == 10000)) )
			&&	 ( (iters ==   100 && MME > 0.40625)	// Even a single ROE unacceptable for such a short run
				|| (iters ==  1000 && MME > 0.42   )
				|| (iters == 10000 && MME >=0.4375 ) ) )
			{
				fprintf(stderr, "***** Excessive level of roundoff error detected - this radix set will not be used. *****\n");

				/* If user-specified radix set, do only that one: */
				if(radset >= 0) goto DONE;

				runtime = 0.0;
				++radix_set;
				continue;
			}
			else if(retVal)	/* Bzzzzzzzzzzzt!!! That answer is incorrect. The penalty is death... */
			{
				printMlucasErrCode(retVal);
				if( !userSetExponent && ((iters == 100) || (iters == 1000) || (iters == 10000)) )
				{
					fprintf(stderr, "Error detected - this radix set will not be used.\n");
				}
				/* If user-specified radix set, do only that one: */
				if(radset >= 0) goto DONE;

				runtime = 0.0;
				++radix_set;
				continue;
			}
			else if(radset >= 0)	/* If user-specified radix set, do only that one: */
			{
				goto DONE;
			}
			else if(new_data)		/* New self-test residue being computed - write detailed SH residues to .cfg file if get a consensus value */
			{
				if(!new_res.sh0)	/* First of the available radix sets... */
				{
					new_res.sh0 = Res64  ;
					new_res.sh1 = Res35m1;
					new_res.sh2 = Res36m1;
					nradix_set_succeed++;	// This assumes when doing such a new-data-selftest the initial radset gives the correct result,
											// which may not be justified, but we can add fancier consensus-weighing logic later, if needed.
				}
				else if				/* All subsequent radix sets must match to produce a consensus result: */
				(
					new_res.sh0 == Res64
				 &&	new_res.sh1 == Res35m1
				 &&	new_res.sh2 == Res36m1
				)
				{
					nradix_set_succeed++;
				}
				else
				{
					runtime = 0.0;
					++radix_set;
					continue;
				}
			} else {	// If not a new-data self-tests (i.e. it's a regular -s one), getting here means the current radset succeeded:
				nradix_set_succeed++;
			}

			/* 16 Dec 2007: Added the (runtime != 0) here to workaround the valid-timing-test-but-runtime = 0
			issue on systems with round-to-nearest-second granularity of the clock() function:
			*/
			if( runtime_best == 0.0 || ((runtime != 0) && (runtime < runtime_best)) )
			{
				runtime_best = runtime;
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
			if(!fp)
			{
				sprintf(cbuf  , "INFO: Unable to open %s file in %s mode ... \n", CONFIGFILE, FILE_ACCESS_MODE);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}

			/* Put code version on line 1.
			Only want to do this once; subsequent mlucas_fopen/fprintf are in append mode:
			*/
			if(new_cfg && FILE_ACCESS_MODE[0]==FILE_ACCESS_WRITE)
			{
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
			if(get_fft_radices(iarg, radix_best, &NRADICES, RADIX_VEC, 10) != 0)
			{
				sprintf(cbuf  , "ERROR: alleged best-radix-set index %u is unsupported.\n",radix_best);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			/* Zero-pad the radices-printing to the full length of the RADIX_VEC array
			so each line has same length (needed to allow update mode):
			*/
			for(i = 0; i < 10; i++){ fprintf(fp,"%3u",RADIX_VEC[i]); };

			/* If it's a new self-test residue being computed, add the SH residues to the .cfg file line */
			if(new_data)
			{
				fprintf(fp, "\t%d-iteration Res mod 2^64, 2^35-1, 2^36-1 = %16llX, %11.0f, %11.0f",iters,new_res.sh0,(double)new_res.sh1,(double)new_res.sh2);
			}

			fprintf(fp,"\n");
			fclose(fp);	fp = 0x0;

			/* if just adding entry for a single FFT length needed for current exponent, return to here: */
			if (quick_self_test)
			{
				quick_self_test = 0;
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
		strncpy(stFlag, global_argv[nargs++], STR_MAX_LEN);
		if(STREQ(stFlag, "-shift"))
		{
		#ifndef USE_SSE2	// Only support residue-shift for SIMD builds:
			fprintf(stderr, "*** INFO: Only support residue-shift for SIMD builds ... ignoring user-set value.\n");
			fprintf(stderr,"%s", cbuf);
			nargs++;
		#else
			strncpy(stFlag, global_argv[nargs++], STR_MAX_LEN);
			/* Convert the shift argument to a uint64: */
			i64arg = 0;
			for(i = 0; i < STR_MAX_LEN && stFlag[i] != '\0'; i++)
			{
				if(isdigit(stFlag[i]))
				{
					i64arg = 10*i64arg + (stFlag[i]-CHAROFFSET);
					/* Check for overflow: */
					if(i64arg % (uint64)10 != (uint64)(stFlag[i]-CHAROFFSET))
					{
						snprintf(cbuf,STR_MAX_LEN, "*** ERROR: -shift argument %s overflows uint64 field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					snprintf(cbuf,STR_MAX_LEN, "*** ERROR: Non-numeric character encountered in -shift argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}
		#endif	// USE_SSE2	?
		}
	}
	return i64arg;
}

/******************/

void print_help(char*option)
{
	int lo, hi, printall = 0;
	if(STREQ(option,"")) printall = TRUE;

	fprintf(stderr, "For the full list of command line options, run the program with the -h flag.\n");
	fprintf(stderr, "For a list of command-line options grouped by type, run the program with the -topic flag.\n\n");

  if(printall)
	fprintf(stderr, "Mlucas command line options:\n\
	\n\
	 Symbol and abbreviation key:\n\
	       <CR> :  carriage return\n\
	        |   :  separator for one-of-the-following multiple-choice menus\n\
	       []   :  encloses optional arguments\n\
	       {}   :  denotes user-supplied numerical arguments of the type noted.\n\
	              ({int} means nonnegative integer, {+int} = positive int, {float} = float.)\n\
	  -argument :  Vertical stacking indicates argument short 'nickname' options,\n\
	  -arg      :  e.g. in this example '-arg' can be used in place of '-argument'.\n\
	\n\
	 Supported arguments:\n\
	\n\
	 <CR>        Default mode: looks for a %s file in the local\n\
	             directory; if none found, prompts for manual keyboard entry\n\
	\n", RANGEFILE);

  if(printall || STREQ(option,"topic")) {
	fprintf(stderr, "Help submenus by topic. No additional arguments may follow the displayed ones:\n");
	fprintf(stderr, " -s            Post-build self-testing for various FFT-length rnages.\n");
	fprintf(stderr, " -fftlen       FFT-length setting.\n");
	fprintf(stderr, " -radset       FFT radix-set specification.\n");
	fprintf(stderr, " -m[ersenne]   Mersenne-number primality testing.\n");
	fprintf(stderr, " -f[ermat]     Fermat-number primality testing.\n");
	fprintf(stderr, " -shift        ***SIMD builds only*** Number of bits by which to shift the initial seed (= iteration-0 residue).\n");
	fprintf(stderr, " -prp          Probable-primality testing mode.\n");
	fprintf(stderr, " -iters        Iteration-number setting.\n");
	fprintf(stderr, " -nthread|cpu  Setting threadcount and CPU core affinity.\n");
	fprintf(stderr, "\n");
  }

  if(printall || STREQ(option,"self_test")) {
	fprintf(stderr, " *** NOTE: *** The following self-test options will cause an mlucas.cfg file containing\n\
     the optimal FFT radix set for the runlength(s) tested to be created (if one did not\n\
     exist previously) or appended (if one did) with new timing data. Such a file-write is\n\
     triggered by each complete set of FFT radices available at a given FFT length being\n\
     tested, i.e. by a self-test without a user-specified -radset argument.\n\
     (A user-specific Mersenne exponent may be supplied via the -m flag; if none is specified,\n\
     the program will use the largest permissible exponent for the given FFT length, based on\n\
     its internal length-setting algorithm). The user must specify the number of iterations for\n\
     the self-test via the -iters flag; while it is not required, it is strongly recommended to\n\
     stick to one of the standard timing-test values of -iters = [100,1000,10000], with the larger\n\
     values being preferred for multithreaded timing tests, in order to assure a decently large\n\
     slice of CPU time. Similarly, it is recommended to not use the -m flag for such tests, unless\n\
     roundoff error levels on a given compute platform are such that the default exponent at one or\n\
     more FFT lengths of interest prevents a reasonable sampling of available radix sets at same.\n\
        If the user lets the program set the exponent and uses one of the aforementioned standard\n\
     self-test iteration counts, the resulting best-timing FFT radix set will only be written to the\n\
     resulting mlucas.cfg file if the timing-test result matches the internally- stored precomputed\n\
     one for the given default exponent at the iteration count in question, with eligible radix sets\n\
     consisting of those for which the roundoff error remains below an acceptable threshold.\n\
     If the user instead specifies the exponent (only allowed for a single-FFT-length timing test)****************\n\
     and/or a non-default iteration number, the resulting best-timing FFT radix set will only be\n\
     written to the resulting mlucas.cfg file if the timing-test results match each other? ********* check logic here *******\n");
	fprintf(stderr, "     This is important for tuning code parameters to your particular platform.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "   FOR BEST RESULTS, RUN ANY SELF-TESTS UNDER ZERO- OR CONSTANT-LOAD CONDITIONS\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -s {...}    Self-test, user must also supply exponent [via -m or -f] and/or FFT length to use.\n");
	fprintf(stderr, "\n");
	lo = 0; hi = numTiny;
	fprintf(stderr, " -s tiny     Runs 100-iteration self-tests on set of %3d Mersenne exponents, ranging from %d to %d\n", numTiny, MersVec[lo].exponent, MersVec[hi-1].exponent);
	fprintf(stderr, " -s t        This will take around 1 minute on a fast CPU..\n");
	fprintf(stderr, "\n");
	lo = 0; hi = numSmall;
	fprintf(stderr, " -s small    Runs 100-iteration self-tests on set of %3d Mersenne exponents, ranging from %d to %d\n", numSmall, MersVec[lo].exponent, MersVec[hi-1].exponent);
	fprintf(stderr, " -s s        This will take around 10 minutes on a fast CPU..\n");
	fprintf(stderr, "\n");
	lo =hi; hi+= numMedium;
	fprintf(stderr, "**** THIS IS THE ONLY SELF-TEST ORDINARY USERS ARE RECOMMENDED TO DO: ******\n");
	fprintf(stderr, "*                                                                          *\n");
	fprintf(stderr, "* -s medium   Runs set of %3d Mersenne exponents, ranging from %d to %d\n", numMedium,MersVec[lo].exponent, MersVec[hi-1].exponent);
	fprintf(stderr, "* -s m        This will take around an hour on a fast CPU.                 *\n");
	fprintf(stderr, "*                                                                          *\n");
	fprintf(stderr, "****************************************************************************\n");
	fprintf(stderr, "\n");
	lo =hi; hi+= numLarge;
	fprintf(stderr, " -s large    Runs set of %3d Mersenne exponents, ranging from %d to %d\n", numLarge, MersVec[lo].exponent, MersVec[hi-1].exponent);
	fprintf(stderr, " -s l        This will take around an hour on a fast CPU.\n");
	fprintf(stderr, "\n");
	lo =hi; hi+= numHuge;
	fprintf(stderr, " -s huge     Runs set of %3d Mersenne exponents, ranging from %d to %d\n", numHuge,  MersVec[lo].exponent, MersVec[hi-1].exponent);
	fprintf(stderr, " -s h        This will take a couple of hours on a fast CPU.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -s all      Runs 100-iteration self-tests of all test Mersenne exponents and all FFT radix sets.\n");
	fprintf(stderr, " -s a        This will take several hours on a fast CPU.\n");
	fprintf(stderr, "\n");
  }
  if(printall || STREQ(option,"fftlen")) {
	fprintf(stderr, " -fftlen {+int}   If {+int} is one of the available FFT lengths (in Kilodoubles), runs all\n");
	fprintf(stderr, "             all available FFT radices available at that length, unless the -radset flag is\n");
	fprintf(stderr, "             invoked (see below for details). If -fftlen is invoked without the -iters flag,\n");
	fprintf(stderr, "             it is assumed the user wishes to do a production run with a non-default FFT length,\n");
	fprintf(stderr, "             In this case the program requires a valid %s-file entry with exponent\n",RANGEFILE);
	fprintf(stderr, "             not more than 5%% larger than the default maximum for that FFT length.\n");
	fprintf(stderr, "                  If -fftlen is invoked with a user-supplied value of -iters but without a\n");
	fprintf(stderr, "             user-supplied exponent, the program will do the specified number of iterations\n");
	fprintf(stderr, "             using the default self-test Mersenne or Fermat exponent for that FFT length.\n");
	fprintf(stderr, "                  If -fftlen is invoked with a user-supplied value of -iters and either the\n");
	fprintf(stderr, "             -m or -f flag and a user-supplied exponent, the program will do the specified\n");
	fprintf(stderr, "             number of iterations of either the Lucas-Lehmer test with starting value 4 (-m)\n");
	fprintf(stderr, "             or the Pe'pin test with starting value 3 (-f) on the user-specified modulus.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "             In either of the latter 2 cases, the program will produce a cfg-file entry based\n");
	fprintf(stderr, "             on the timing results, assuming at least one radix set ran the specified #iters\n");
	fprintf(stderr, "             to completion without suffering a fatal error of some kind.\n");
	fprintf(stderr, "             Use this to find the optimal radix set for a single FFT length on your hardware.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "             NOTE: IF YOU USE OTHER THAN THE DEFAULT MODULUS OR #ITERS FOR SUCH A SINGLE-FFT-\n");
	fprintf(stderr, "             LENGTH TIMING TEST, IT IS UP TO YOU TO MANUALLY VERIFY THAT THE RESIDUES OUTPUT\n");
	fprintf(stderr, "             MATCH FOR ALL FFT RADIX COMBINATIONS AND THE ROUNDOFF ERRORS ARE REASONABLE!\n");
	fprintf(stderr, "\n");
  }
  if(printall || STREQ(option,"radset")) {
	fprintf(stderr, " -radset {int}    Specific index of a set of complex FFT radices to use, based on the big\n");
	fprintf(stderr, "             select table in the function get_fft_radices(). Requires a supported value of\n");
	fprintf(stderr, "             -fftlen to also be specified, as well as a value of -iters for the timing test.\n");
	fprintf(stderr, "\n");
  }
  if(printall || STREQ(option,"mersenne")) {
	fprintf(stderr, " -m [{+int}] Performs a Lucas-Lehmer primality test of the Mersenne number M(int) = 2^int - 1,\n");
	fprintf(stderr, "             where int must be an odd prime. If -iters is also invoked, this indicates a timing test.\n");
	fprintf(stderr, "             and requires suitable added arguments (-fftlen and, optionally, -radset) to be supplied.\n");
	fprintf(stderr, "                If the -fftlen option (and optionally -radset) is also invoked but -iters is not, the\n");
	fprintf(stderr, "             program first checks the first line of the worktodo.ini file to see if the assignment\n");
	fprintf(stderr, "             specified there is a Lucas-Lehmer test with the same exponent as specified via the -m\n");
	fprintf(stderr, "             argument. If so, the -fftlen argument is treated as a user override of the default FFT\n");
	fprintf(stderr, "             length for the exponent. If -radset is also invoked, this is similarly treated as a user-\n");
	fprintf(stderr, "             specified radix set for the user-set FFT length; otherwise the program will use the cfg file\n");
	fprintf(stderr, "             to select the radix set to be used for the user-forced FFT length.\n");
	fprintf(stderr, "                If the worktodo.ini file entry does not match the -m value, a set of timing self-tests is\n");
	fprintf(stderr, "             run on the user-specified Mersenne number using all sets of FFT radices available at the\n");
	fprintf(stderr, "             specified FFT length.\n");
	fprintf(stderr, "                If the -fftlen option is not invoked, the self-tests use all sets of\n");
	fprintf(stderr, "             FFT radices available at that exponent's default FFT length.\n");
	fprintf(stderr, "                Use this to find the optimal radix set for a single given Mersenne number\n");
	fprintf(stderr, "             exponent on your hardware, similarly to the -fftlen option.\n");
	fprintf(stderr, "                Performs as many iterations as specified via the -iters flag [required].\n");
	fprintf(stderr, "\n");
  }
  if(printall || STREQ(option,"fermat")) {
	fprintf(stderr, " -f {int}    Performs a base-3 Pe'pin test on the Fermat number F(num) = 2^(2^num) + 1.\n");
	fprintf(stderr, "                If desired this can be invoked together with the -fftlen option.\n");
	fprintf(stderr, "             as for the Mersenne-number self-tests (see notes about the -m flag;\n");
	fprintf(stderr, "             note that not all FFT lengths supported for -m are available for -f).\n");
	fprintf(stderr, "             Optimal radix sets and timings are written to a fermat.cfg file.\n");
	fprintf(stderr, "                Performs as many iterations as specified via the -iters flag [required].\n");
	fprintf(stderr, "\n");
  }
  if(printall || STREQ(option,"shift")) {
	fprintf(stderr, " -shift         ***SIMD builds only*** Bits by which to circular-left-shift the initial seed.\n");
	fprintf(stderr, "             This shift count is doubled (modulo the number of bits of the modulus being tested)\n");
	fprintf(stderr, "             each iteration. Savefile residues are rightward-shifted by the current shift count\n");
	fprintf(stderr, "             before being written to the file; thus savefiles contain the unshifted residue, and\n");
	fprintf(stderr, "             separately the current shift count, which the program uses to leftward-shift the\n");
	fprintf(stderr, "             savefile residue when the program is restarted from interrupt.\n");
	fprintf(stderr, "                The shift count is a 64-bit unsigned int (e.g. to accommodate Fermat numbers > F32).\n");
	fprintf(stderr, "\n");
  }
  if(printall || STREQ(option,"prp")) {
	fprintf(stderr, " -prp {int}     Instead of running the rigorous primality test defined for the modulus type\n");
	fprintf(stderr, "             in question (Lucas-Lehmer test for Mersenne numbers, Pe'pin test for Fermat numbers\n");
	fprintf(stderr, "             do a probably-primality test to the specified integer base b = {int}.\n");
	fprintf(stderr, "                For a Mersenne number M(p), starting with initial seed x = b (which must not = 2\n");
	fprintf(stderr, "             or a power of 2), this means do a Fermat-PRP test, consisting of (p-2) iterations of\n");
	fprintf(stderr, "             form x = b*x^2 (mod M(p)) plus a final mod-squaring x = x^2 (mod M(p)), with M(p) being\n");
	fprintf(stderr, "             a probable-prime to base b if the result == 1.\n");
	fprintf(stderr, "                For a Fermat number F(m), starting with initial seed x = b (which must not = 2\n");
	fprintf(stderr, "             or a power of 2), this means do an Euler-PRP test (referred to as a Pe'pin test for these\n");
	fprintf(stderr, "             moduli), i.e. do 2^m-1 iterations of form x = b*x^2 (mod M(p)), with M(p) being not merely\n");
	fprintf(stderr, "             a probable prime but in fact deterministically a prime if the result == -1. The reason we\n");
	fprintf(stderr, "             still use the -prp flag in the Fermat case is for legacy-code compatibility: All pre-v18\n");
	fprintf(stderr, "             Mlucas versions supported only Pe'pin testing to base b = 3; now the user can use the -prp\n");
	fprintf(stderr, "             flag with a suitable base-value to override this default choice of base.\n");
	fprintf(stderr, "\n");
  }
  if(printall || STREQ(option,"iters")) {
	fprintf(stderr, " -iters {int}   Do {int} self-test iterations of the type determined by the\n");
	fprintf(stderr, "             modulus-related options (-s/-m = Lucas-Lehmer test iterations with\n");
	fprintf(stderr, "             initial seed 4, -f = Pe'pin-test squarings with initial seed 3.\n");
	fprintf(stderr, "\n");
  }
 #ifdef MULTITHREAD
  if(printall || STREQ(option,"nthread")) {
	fprintf(stderr, " -nthread {int}   For multithread-enabled builds, run with this many threads.\n\
     If the user does not specify a thread count, the default is to run single-threaded\n\
     with that thread's affinity set to logical core 0.\n\
     \n\
     AFFINITY: The code will attempt to set the affinity of the resulting threads\n\
     0:n-1 to the same-indexed processor cores - whether this means distinct physical\n\
     cores is entirely up to the CPU vendor - E.g. Intel uses such a numbering scheme\n\
     but AMD does not. For this reason as of v17 this option is deprecated in favor of\n\
     the -cpu flag, whose usage is detailed below, with the online README page providing\n\
     guidance for the core-numbering schemes of popular CPU vendors.\n\
     \n\
     If n exceeds the available number of logical processor cores (call it #cpu), the\n\
     program will halt with an error message.\n\
     \n\
     For greater control over affinity setting, use the -cpu option, which supports two\n\
     distinct core-specification syntaxes (which may be mixed together), as follows:\n\
     \n\
     -cpu {lo[:hi[:incr]]}   (All args {int} here) Set thread/CPU affinity.\n\
     NOTE: This flag and -nthread are mutually exclusive: If -cpu is used, the threadcount\n\
     is inferred from the numeric-argument-triplet which follows. If only the 'lo' argument\n\
     of the triplet is supplied, this means 'run single-threaded with affinity to CPU {lo}.'\n\
     If the increment (third) argument of the triplet is omitted, it is taken as incr = 1.\n\
     The CPU set encoded by the integer-triplet argument to -cpu corresponds to the\n\
     values of the integer loop index i in the C-loop for(i = lo; i <= hi; i += incr),\n\
     excluding the loop-exit value of i. Thus '-cpu 0:3' and '-cpu 0:3:1' are both\n\
     exactly equivalent to '-nthread 4', whereas '-cpu 0:6:2' and '-cpu 0:7:2' both\n\
     specify affinity setting to cores 0,2,4,6, assuming said cores exist.\n\
     Lastly, note that no whitespace is permitted within the colon-separated numeric field.\n\
     \n\
     -cpu {triplet0[,triplet1,...]}   This is simply an extended version of the above affinity-\n\
     setting syntax in which each of the comma-separated 'triplet' subfields is in the above\n\
     form and, analogously to the one-triplet-only version, no whitespace is permitted within\n\
     the colon-and-comma-separated numeric field. Thus '-cpu 0:3,8:11' and '-cpu 0:3:1,8:11:1'\n\
     both specify an 8-threaded run with affinity set to the core quartets 0-3 and 8-11,\n\
     whereas '-cpu 0:3:2,8:11:2' means run 4-threaded on cores 0,2,8,10. As described for the\n\
     -nthread option, it is an error for any core index to exceed the available number of logical\n\
     processor cores.\n");
  }
 #endif
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

/* These need to be kept updated to match the #defines in Mdata.h: */
const char *err_code[] = {
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
	"ERR_INTERRUPT"
};

void	returnMlucasErrCode(uint32 ierr, char*s)
{
	ASSERT(HERE, s != 0x0, "Null string pointer!");
	ASSERT(HERE, ierr < ERR_MAX, "Error code out of range!");
	if(ierr == 0)
		s[0] = '\0';
	else
		strcpy(s, err_code[ierr-1]);
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
		ASSERT(HERE, i==ERR_RUN_SELFTEST_FORLENGTH, "High bytes should only be nonzero if low byte == ERR_RUN_SELFTEST_FORLENGTH!");
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

The input arguments Res64, Res35m1, Res36m1 are 3 checksums (residue modulo 2^64, 2^35-1, 2^36-1)
intended to be compared to the same gotten via subsequent processing (e.g. conversion to floating-
point form) of the residue array itself. (These are essentially the Selfridge-Hurwitz residues
commonly used in Fermat-number primality testing, but with the mod-2^36 component widened to 64 bits).
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
	itmp64 = ((uint64*)arr_tmp)[0];	if(*Res64 != itmp64) { sprintf(cbuf, "%s: On restart: Res64 checksum error! Got %llX, expected %llX\n"  ,func,itmp64,*Res64); return 0; }
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
	itmp64 = mi64_div_by_scalar64((uint64*)arr_tmp,two35m1,i,0x0);	if(*Res35m1 != itmp64)	{ sprintf(cbuf, "%s: On restart: Res35m1 checksum error! Got %llX, expected %llX\n",func,itmp64,*Res35m1); return 0; }
	itmp64 = mi64_div_by_scalar64((uint64*)arr_tmp,two36m1,i,0x0);	if(*Res36m1 != itmp64)	{ sprintf(cbuf, "%s: On restart: Res36m1 checksum error! Got %llX, expected %llX\n",func,itmp64,*Res36m1); return 0; }
  #endif
	return 1;
}

// Returns 1 on successful read, 0 otherwise:
// v19: For PRP-tests, also write a second Gerbicz-check residue array [arr2] and associated S-H checksum triplet [i1,i2,i3]
int read_ppm1_savefiles(uint64 p, int*kblocks, FILE*fp, uint32*ilo,
	uint8 arr1[], uint64*Res64, uint64*Res35m1, uint64*Res36m1,
	uint8 arr2[], uint64*i1   , uint64*i2     , uint64*i3     )
{
	const char func[] = "read_ppm1_savefiles";
	uint32 i,j,nbytes;
	uint64 itmp64, nsquares = 0ull;

	*Res64 = 0ull;	// 0 value on return indicates failure of some kind

	sprintf(cbuf, "%s: p must be 32-bit or less!",func);
	ASSERT(HERE, !(p >> 32), cbuf);	/* Future versions will need to loosen this p < 2^32 restriction: */

	if(!file_valid(fp)) {
		sprintf(cbuf, "%s: File pointer invalid for read!\n",func);	ASSERT(HERE, 0, cbuf);
	}
	fprintf(stderr, " INFO: restart file %s found...reading...\n",RESTARTFILE);
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
	/* For now, just make sure nsquares < 2^32 and copy to ilo: */
	if(nsquares >= p) {
		sprintf(cbuf,"%s: nsquares = %llu out of range, should be < p = %llu\n",func, nsquares, p);
		return 0;
	} else if(nsquares > 0xFFFFFFFFull) {
		sprintf(cbuf,"%s: nsquares = %llu out of range, current limit = 2^32-1.\n",func, nsquares);
		return 0;
	}
	*ilo = nsquares;

	/* Set the expected number of residue bytes, depending on the modulus: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		nbytes = (p + 7)/8;
		TRANSFORM_TYPE = REAL_WRAPPER;
	} else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
		sprintf(cbuf, "%s: MODULUS_TYPE_FERMAT but (p mod 8) != 0",func);
		ASSERT(HERE, (p & 7) == 0,cbuf);
		nbytes = (p>>3) + 1;
		TRANSFORM_TYPE = RIGHT_ANGLE;
	}

	i = read_ppm1_residue(nbytes, fp, arr1, Res64,Res35m1,Res36m1);
	if(!i) return 0;

	// FFT length in K (3 bytes) - first added this in v18:
	*kblocks = 0;
	for(j = 0; j < 3; j++) {
		i = fgetc(fp);	*kblocks += (uint64)i << (8*j);
	}
	/* May 2018: 8 bytes for circular-shift to apply to the (unshifted) residue read from the file: */
	RES_SHIFT = 0ull;
	for(j = 0; j < 8; j++) {
		i = fgetc(fp);	RES_SHIFT += (uint64)i << (8*j);
	}

  // v19: For PRP-tests, also read a second Gerbicz-check residue array [arr2] and associated S-H checksum triplet [i1,i2,i3]:
  if(TEST_TYPE == TEST_TYPE_PRP) {
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
		fclose(fp);
		snprintf(cbuf,STR_MAX_LEN,"%s: Error writing residue to restart file %s.\n",func,RESTARTFILE);
									        fprintf(stderr,"%s",cbuf);
		fp = mlucas_fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
		fp = mlucas_fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
		ASSERT(HERE, 0,"");
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

void write_ppm1_savefiles(uint64 p, int n, FILE*fp, uint32 ihi,
	uint8 arr1[], uint64 Res64, uint64 Res35m1, uint64 Res36m1,
	uint8 arr2[], uint64 i1   , uint64 i2     , uint64 i3     )
{
	uint32 i,kblocks,nbytes = 0;
	uint64 itmp64;

	ASSERT(HERE,file_valid(fp),"write_ppm1_savefiles: File pointer invalid for write!");
	// Make sure n is a proper (unpadded) FFT-length, i.e. is a multiple of 1K:
	kblocks = (n >> 10);
	ASSERT(HERE,n == (kblocks << 10),"Not a proper unpadded FFT length");

	/* See the function read_ppm1_savefiles() for the file format here: */
	/* t: */
	fputc(TEST_TYPE, fp);
	/* m: */
	fputc(MODULUS_TYPE, fp);
	/* s: */
	itmp64 = ihi;
	for(i = 0; i < 64; i+=8)
		fputc((itmp64 >> i) & 0xff, fp);

	/* Set the expected number of residue bytes, depending on the modulus: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		nbytes = (p + 7)/8;
		TRANSFORM_TYPE = REAL_WRAPPER;
	} else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT) {
		ASSERT(HERE, p % 8 == 0,"write_ppm1_savefiles: p % 8 == 0");
		nbytes = (p/8) + 1;
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
  if(TEST_TYPE == TEST_TYPE_PRP) {
	for(i = 0; i < 32; i+=8)
		fputc((PRP_BASE >> i) & 0xff, fp);
	write_ppm1_residue(nbytes, fp, arr2, i1,i2,i3);
	// G-check residues all need to be clshifted by residue-shift count at the ITERS_BETWEEN_GCHECK_UPDATESth PRP-test iteration:
	for(i = 0; i < 64; i+=8)
		fputc((GCHECK_SHIFT >> i) & 0xff, fp);
  }
}

/*********************/

/*
Function to take an n-digit Mersenne or Fermat-mod residue in bytewise savefile form,
apply the required circular shift read into the global RES_SHIFT during the preceding
bytewise-savefile read and convert it to balanced-digit floating-point form.

In the Mersenne-mod case the residue digits are stored
consecutively in the a[] array.

In the Fermat-mod case the digits are arranged in (j,j+n/2)
(i.e. right-angle transform) order.
*/
int 	convert_res_bytewise_FP(const uint8 arr_tmp[], double a[], int n, const uint64 p)
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

	ASSERT(HERE,MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(HERE,MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	ASSERT(HERE,TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(HERE,TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");

	ASSERT(HERE, !(p >> 32), "p must be 32-bit or less!");	/* Future versions will need to loosen this p < 2^32 restriction: */

	/* Set the number of residue bytes, which is the same for Mersenne (2^p-1) and Fermat-mod (2^p+1, with p = 2^findex)
	despite the fact the latter can formally be as large as 2^p, since only ever hit that if it`s the last residue of
	a Pepin test and the number hqppens to be prime. (We would love for that exception to break some other ASSERTion in the code): */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER,"convert_res_bytewise_FP: TRANSFORM_TYPE == REAL_WRAPPER");
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, TRANSFORM_TYPE == RIGHT_ANGLE,"convert_res_bytewise_FP: TRANSFORM_TYPE == RIGHT_ANGLE");
		/* If Fermat number, make sure exponent a power of 2: */
		findex = trailz64(p);
		ASSERT(HERE, (p >> findex) == 1,"convert_res_bytewise_FP: (p >> findex) == 1");

		ASSERT(HERE, p % 8 == 0,"convert_res_bytewise_FP: p % 8 == 0");
	}
	nbytes = (p + 7)/8;

	// Apply the circular shift:
	mi64_shlc((uint64*)arr_tmp, (uint64*)arr_tmp, p, RES_SHIFT, (p+63)>>6);

	/* Vector length a power of 2? */
	pow2_fft = (n >> trailz32(n)) == 1;

	bits[0] = p/n;		ASSERT(HERE, bits[0] > 1,"convert_res_bytewise_FP: bits[0] > 1");
	base[0] = 1 << bits[0];

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT && pow2_fft == TRUE)
	{
		bits[1] =     bits[0];
	}
	else
	{
		bits[1] = 1 + bits[0];
	}
	base[1] = 1 << bits[1];

	bw = p%n;	/* cf. mers_mod_square.c	*/
	sw = n - bw;

	/*...Now form the SH residues, converting to positive-digit form along the way...	*/

	curr_char = 0;	/* Current byte to be read from the arr_tmp array */
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
			the bits we have to curr_word, grab the next 8 bytes of arr_tmp[] and finish
			the write of the remaining (bits[ii] - rbits) of curr_word:
			*/
			if(rbits < bits[ii])
			{
				itmp = curr_wd64;
				ASSERT(HERE, itmp < (1ull<<rbits),"convert_res_bytewise_FP: itmp >= 2^rbits!");

				/* Now grab the next 64 bits of the bytewise residue... */
				curr_wd64 = 0;
				for(k = 0; k < 8; k++)
				{
					curr_wd64 += (uint64)arr_tmp[curr_char++] << (k<<3);	/* left-shift current residue byte k*8 bits and accumulate */
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
			the bits we have to curr_word, grab the next 8 bytes of arr_tmp[] and finish
			the write of the remaining (bits[ii] - rbits) of curr_word:
			*/
			if(rbits < bits[ii])
			{
				itmp = curr_wd64;
				ASSERT(HERE, itmp < (1<<rbits),"convert_res_bytewise_FP: itmp >= 2^rbits!");

				/* Now grab the next 64 bits of the bytewise residue... */
				curr_wd64 = 0;
				for(k = 0; k < 8; k++)
				{
					curr_wd64 += (uint64)arr_tmp[curr_char++] << (k<<3);	/* left-shift current residue byte k*8 bits and accumulate */
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

	ASSERT(HERE, curr_char == nbytes, "convert_res_bytewise_FP: curr_char == (p+7)/8");
	ASSERT(HERE, nbits == p    ,"convert_res_bytewise_FP: nbits == p    ");
	ASSERT(HERE, curr_wd64 == 0,"convert_res_bytewise_FP: curr_word == 0");

	/*
	Fold any carryout from the conversion to balanced-representation form
	into the LS word of the residue (e.g. if cy = 1, this amounts to subtracting
	the modulus from the positive-digit form to get the balanced-digit form):
	*/
	/* Should have carryout of +1 Iff MS word < 0; otherwise expect 0 carry: */
	if(cy && (a[j1]>= 0 || cy != +1))
	{
		sprintf(cbuf, "convert_res_bytewise_FP: Illegal combination of nonzero carry = %lld, most sig. word = %20.4f\n", cy, a[j]);
		ASSERT(HERE, 0, cbuf);
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		a[0] += cy;
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		a[0] -= cy;
	else
		ASSERT(HERE, 0,"Illegal modulus type!");
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
void	convert_res_FP_bytewise(const double a[], uint8 arr_tmp[], int n, const uint64 p, uint64*Res64, uint64*Res35m1, uint64*Res36m1)
{
	int bimodn,curr_bits,curr_char,cy,findex,ii,j,j1,k,pass,rbits,msw_lt0,bw,sw,bits[2],pow2_fft;
	uint64 nbits,curr_wd64,base[2];	/* Assume base may be > 2^32 (e.g. for mixed FFT/FGT) but < 2^53, i.e. fits in a double */
	double atmp;
	int64 itmp;
	const uint64 two35m1 = (uint64)0x00000007FFFFFFFFull, two36m1 = (uint64)0x0000000FFFFFFFFFull;	/* 2^35,36-1 */

	ASSERT(HERE,MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(HERE,MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	ASSERT(HERE,TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(HERE,TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");

	/* If Fermat number, make sure exponent a power of 2: */
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, TRANSFORM_TYPE == RIGHT_ANGLE,"convert_res_FP_bytewise: TRANSFORM_TYPE == RIGHT_ANGLE");
		findex = trailz64(p);
		ASSERT(HERE, (p >> findex) == 1,"convert_res_FP_bytewise: (p >> findex) == 1");
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER,"convert_res_FP_bytewise: TRANSFORM_TYPE == REAL_WRAPPER");
	else
		ASSERT(HERE, 0,"Illegal modulus type!");

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
	for(j=n-1; j >= 0; j -= TRANSFORM_TYPE)
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
		if(atmp != 0.0)
		{
			if(atmp < 0.0)
			{
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
			if(atmp != 0.0)
			{
				if(atmp < 0.0)
				{
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

	curr_char = 0;	/* Current byte to be written to in the arr_tmp array */
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
			ASSERT(HERE, atmp == NINT(atmp), "Input float-residue elements must have 0 fractional part!");
			itmp = (int64)(atmp+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0) {			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
				itmp += (base[ii]);
				cy = -1;
			} else {
				cy = 0;
			}
			ASSERT(HERE, itmp >= 0,"convert_res_FP_bytewise: itmp >= 0");

		/* Update 8-byte residue buffer last, since this one modifies itmp: */
			ASSERT(HERE, rbits < 8,"convert_res_FP_bytewise: rbits < 8");
			ASSERT(HERE, curr_wd64 < (1<<rbits),"convert_res_FP_bytewise: curr_wd64 >= 2^rbits!");

			itmp = (itmp << rbits) + curr_wd64;
			curr_bits = bits[ii] + rbits;

			rbits = curr_bits;
			for(k = 0; k < curr_bits/8; k++)
			{
				arr_tmp[curr_char++] = itmp & 255;
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
			ASSERT(HERE, atmp == NINT(atmp), "Input float-residue elements must have 0 fractional part!");
			itmp = (int64)(atmp+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0) {			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
				itmp += (base[ii]);
				cy = -1;
			} else {
				cy = 0;
			}
			ASSERT(HERE, itmp >= 0,"convert_res_FP_bytewise: itmp >= 0");

		/* Update 8-byte residue buffer last, since this one modifies itmp: */
			ASSERT(HERE, rbits < 8,"convert_res_FP_bytewise: rbits < 8");
			ASSERT(HERE, curr_wd64 < (1<<rbits),"convert_res_FP_bytewise: curr_wd64 >= 2^rbits!");

			itmp = (itmp << rbits) + curr_wd64;
			curr_bits = bits[ii] + rbits;

			rbits = curr_bits;
			for(k = 0; k < curr_bits/8; k++)
			{
				arr_tmp[curr_char++] = itmp & 255;
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
		ASSERT(HERE, 0, cbuf);
	}

	/* Residue should contain ceiling(p/8) bytes: */
	ASSERT(HERE, rbits < 8, "rbits >= 8");
	if(rbits)
	{
		ASSERT(HERE, curr_wd64 < (1<<rbits),"convert_res_FP_bytewise: curr_wd64 >= 2^rbits!");
		arr_tmp[curr_char++] = curr_wd64 & 255;
		curr_wd64 >>= 8;
	}
	ASSERT(HERE, curr_char == (p+7)/8,"convert_res_FP_bytewise: curr_char == (p+7)/8");
	ASSERT(HERE, nbits == p          ,"convert_res_FP_bytewise: nbits == p          ");
	ASSERT(HERE, curr_wd64 == 0      ,"convert_res_FP_bytewise: curr_wd64 == 0      ");

	// Remove the circular shift ... have no mi64_rshc function, so use that b-bit rightward cshift equivalent to (p-b)-bit left-cshift.
	// (But must guard against RES_SHIFT = 0, since in that case the left-shift count == p and mi64_shlc requires shift count strictly < p):
	j = (p+63)>>6;	// # of 64-bit limbs
	if(RES_SHIFT) {
	//	fprintf(stderr,"convert_res_FP_bytewise: removing shift = %llu\n",RES_SHIFT);
		mi64_shlc((uint64*)arr_tmp, (uint64*)arr_tmp, p, p-RES_SHIFT,j);
	}
	/* Checksums: */
	if(Res64  ) *Res64 = ((uint64*)arr_tmp)[0];
	if(Res35m1) *Res35m1 = mi64_div_by_scalar64((uint64*)arr_tmp,two35m1,j,0x0);
	if(Res36m1) *Res36m1 = mi64_div_by_scalar64((uint64*)arr_tmp,two36m1,j,0x0);
}

/*********************/

uint32 	get_default_factoring_depth(uint64 p)
{
/* Sample: here's how to set things to factor to a constant k-depth: */
#if 1
	const uint32 KMAX_BITS = 40;
	return (uint32) ceil(log(1.0*p)/log(2.0)) + 1 + KMAX_BITS;

#else

	uint32 qbitsmax;

/* These default depths are designed to match those of Prime95 v24, as described here:

	http://www.mersenneforum.org/showthread.php?t=4213
*/
****reverse order, add cases for > pmax, < pmin ***
	else if(p <=23390000)
	{
		qbitsmax = 66;
	}
	else if(p > 23390000)	/* k ~= 40.5 bits */
	{
		qbitsmax = 66;
	}
	else if(p > 29690000)
	{
		qbitsmax = 67;
	}
	else if(p > 37800000)
	{
		qbitsmax = 68;
	}
	else if(p > 47450000)
	{
		qbitsmax = 69;
	}
	else if(p > 58520000)
	{
		qbitsmax = 70;
	}
	else if(p > 75670000)
	{
		qbitsmax = 71;
	}
	else if(p > 96830000)	/* k ~= 44.5 bits at the starting value of p*/
	{
		qbitsmax = 72;
	}

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
	ASSERT(HERE, s != 0x0, "Null ptr to is_hex_string()");
	for(i = 0; i < len; ++i)
	{
		if( !isxdigit(s[i]) )
			return FALSE;
	}
	return TRUE;
}

/*********************/

/* v19: Add basic JSON-formatted result report. Examples:
{"status":"C", "exponent":86749043, "worktype":"LL", "res64":"9EEA7CAD97A07648", "fft-length":4718592, "shift-count":6030412, "error-code":"00000000", "program":{"name":"Mlucas", "version":"18.0"}, "timestamp":"2019-11-11 01:23:45", "user":"madpoo", "computer":"mediaboy", "aid":"ABCDEF012456789 ABCDEF012456789"}
{"status":"C", "exponent":110527, "worktype":"PRP-3", "res64":"E95075F756DD7BEB", "residue-type":1, "res2048":"E5E2DB84978E0355041AE377E588931B54FC75DCAD705044F21F17D0C8D5F524E98C535101C6DA9799F1433934FBAC2090761B1F4D8EA1F91AD63D03D477312E42F1CE7666C5E776A49A5BBDA146543C3CB1D74E0400CF6E81DF35173741289C76E69DB909726E50ECEE697F69A92E6BDF27A6AC6C9591EF97753F4555BABBFB26A385F78497ACAA4F7738A1E3C01564975DBDD3306C89FE7946B1523698BA334FAF2F53D74060BFDDDAF9D643E116D50DA7FDF2EB96CBC5D074602FDEDBD88E2706E0DED9324CEC6AD702016547E300748D6E5685C123CBB93744176B340B7DA6C1E478A685C774D554AAE5335C1FEBADD07A382A33BE1CE95075F756DD7BEB", "fft-length":6144, "shift-count":97673, "error-code":"00000000", "security-code":"F580BEFA", "program":{"name":"Prime95", "version":"29.8", "build":4, "port":10}, "timestamp":"2019-10-24 19:44:06", "errors":{"gerbicz":0}, "user":"gw_2", "computer":"Macbook_Pro"}
*/
void generate_JSON_report(const uint32 isprime, const uint64 p, const uint32 n, const uint64 Res64, const char*timebuffer, char*cstr)
{
	char ttype[11] = "\0", aid[33] = "\0";	// aid needs 33rd chard for \0
	const char status[2] = {'C','P'};
	// Attempt to read 32-hex-char Primenet assignment ID for current assignment (first line of RANGEFILE):
	fp = mlucas_fopen(RANGEFILE, "r");
	ASSERT(HERE,fp != 0x0,"Rangefile not found!");
	fgets(in_line, STR_MAX_LEN, fp);
	fclose(fp);	fp = 0x0;
	// Is there a Primenet-style 32-hexit assignment ID in the assignment line? If so, include it in the JSON output:
	char_addr = strstr(in_line, "=");
	if(char_addr) {
		char_addr++;
		while(isspace(*char_addr)) { ++char_addr; }	// Skip any whitespace following the equals sign
		ASSERT(HERE, is_hex_string(char_addr, 32), "Expect a 32-hex-digit PrimeNet v5 assignment ID following the work type specifier!");
		strncpy(aid,char_addr,32);
	}
	// Write the result line. The 2 nested conditionals here are LL-or-PRP and AID-found-or-not:
	if(TEST_TYPE == TEST_TYPE_PRIMALITY) {
		snprintf(ttype,10,"LL");
		if(char_addr) {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"worktype\":\"%s\", \"res64\":\"%16llX\", \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\", \"aid\":\"%s\"}\n",status[isprime],p,ttype,Res64,n,RES_SHIFT,VERSION,timebuffer,aid);
		} else {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"worktype\":\"%s\", \"res64\":\"%16llX\", \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\"}\n",status[isprime],p,ttype,Res64,n,RES_SHIFT,VERSION,timebuffer);
		}
	} else if(TEST_TYPE == TEST_TYPE_PRP) {	// Only support type-1 PRP tests, so hardcode that subfield:
		snprintf(ttype,10,"PRP-%u",PRP_BASE);
		if(char_addr) {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"worktype\":\"%s\", \"res64\":\"%16llX\", \"residue-type\":1, \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\", \"aid\":\"%s\"}\n",status[isprime],p,ttype,Res64,n,RES_SHIFT,VERSION,timebuffer,aid);
		} else {
			snprintf(cstr,STR_MAX_LEN,"{\"status\":\"%c\", \"exponent\":%llu, \"worktype\":\"%s\", \"res64\":\"%16llX\", \"residue-type\":1, \"fft-length\":%u, \"shift-count\":%llu, \"error-code\":\"00000000\", \"program\":{\"name\":\"Mlucas\", \"version\":\"%s\"}, \"timestamp\":\"%s\"}\n",status[isprime],p,ttype,Res64,n,RES_SHIFT,VERSION,timebuffer);
		}
	} else
		ASSERT(HERE, 0, "Unsupported test type!");
}

/*********************/


/*********************/


/*********************/

void write_fft_debug_data(double a[], int jlo, int jhi)
{
	int j,j1;
	const char dbg_fname[] = "FFT_DEBUG.txt";
	ASSERT(HERE, dbg_file == 0x0, "dbg_file != 0x0 prior to mlucas_fopen");
	dbg_file = mlucas_fopen(dbg_fname, "a");
	ASSERT(HERE, dbg_file != 0x0, "Unable to open dbg_file!");
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

