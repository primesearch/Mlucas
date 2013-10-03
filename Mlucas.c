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

/*********************************************************************************
*   COMPILING AND RUNNING THE PROGRAM: see http://hogranch.com/mayer/README.html *
*********************************************************************************/
#include "Mlucas.h"

int	is_hex_string(char*s, int len);

/*** DEBUG ***/
#define 	OLD_FORMAT	0
#if OLD_FORMAT
void	write_ppm1_savefiles(uint64 p, FILE*fp, uint32 ihi, uint32 n, int32 arr_tmp[], double a[]);
#else
void	write_ppm1_savefiles(uint64 p, FILE*fp, uint32 ihi, uint8 arr_tmp[], uint64 Res64, uint64 Res35m1, uint64 Res36m1);
#endif

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

FILE *dbg_file = 0x0;
double*ADDR0 = 0x0;	// Allows for easy debug on address-read-or-write than setting a watchpoint

/* Define FFT-related globals (declared in Mdata.h) */
uint32 N2,NRT,NRT_BITS,NRTM1;
int NRADICES, RADIX_VEC[10];	/* RADIX_VEC[] stores sequence of complex FFT radices used.	*/

int ITERS_BETWEEN_CHECKPOINTS;	/* number of iterations between checkpoints */

char ESTRING[STR_MAX_LEN];	/* Exponent in string form */
char PSTRING[STR_MAX_LEN];	/* Number being tested in string form, typically estring concatenated with several other descriptors, e.g. strcat("M",estring) */

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
char cbuf[STR_MAX_LEN];
char in_line[STR_MAX_LEN];
char *char_addr;
int char_offset;
FILE *fp, *fq;

/* File access mode is a 3-character string, with the last of these being a mandatory null terminator: */
char FILE_ACCESS_MODE[3] = {'x','y','\0'};

/* Matrix of supported moduli/test-types - init'ed in Mlucas_init */
int ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE_DIM][TEST_TYPE_DIM];

/* Allocate storage for Globals */

/* These should all be set to a valid (nonzero) value at the time the appropriate test is begun */
uint32 TEST_TYPE		= 0;
uint32 MODULUS_TYPE		= 0;
uint32 TRANSFORM_TYPE	= 0;

/* Program version with patch suffix:

For the foreseeable future, version numbers will be of form x.yz, with x a 1-digit integer,
y an int having 3 digits or less, and z an optional alphabetic patch suffix. We choose these such that
retiming is only necessary between versions that differ in x or in the leading digit of y,
i.e. we'd need to regenerate .cfg files in going from version 2.9 to 3.0 or from 3.0 to 3.141,
but not between versions 3.0x and 3.0g or between 3.1 and 3.141.

This reduces the "retime?" decision in the cfgNeedsUpdating() function to a simple comparison
of the leading 3 characters of the two version strings in question.
*/
const char VERSION   [] = "3.0x";			/* <--- patch suffix of x, y, or z indicates a beta code. */

const char OFILE     [] = "results.txt";	/* ASCII logfile containing FINAL RESULT ONLY for each
											assignment - detailed intermediate results for each assignment
											are written to the exponent-specific STATFILE (see below). */
const char RANGEFILE [] = "worktodo.ini";	/* File containing exponents to be tested. New exponents
											may be appended at will, while the program is running. */

const char LOCAL_INI_FILE[] = "local.ini";	/* File containing user-customizable configuration settings */

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

#ifdef MULTITHREAD
	const char ntfile[] = "nthreads.ini";	/* Allows us to control the no. of threads dynamically at runtime. */
#endif

int INTERACT;

double AME,MME;	/* Avg and Max per-iteration fractional error for a given iteration interval */

int MAX_THREADS = 0;		/* Max. allowable No. of threads. */
int NTHREADS = 0;			/* actual No. of threads. If multithreading disabled, set = 1. */

uint32 PMIN;		/* minimum exponent allowed */
uint32 PMAX;		/* maximum exponent allowed depends on max. FFT length allowed
					   and will be determined at runtime, via call to given_N_get_maxP(). */

/****************/

/*
!...Code to test primality of Mersenne numbers, using arbitrary-precision (array-integer) arithmetic.
!   Author: Ernst W. Mayer.
!
!   This version (v3.0x) dated 06 Dec 2012.
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
!   See the Whatsnew.txt file for a brief overview of the various module dependencies and revision history of the code
!   along with a summary of what has been changed/added/deleted in the latest release.
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
!   (1) Continue to optimize SSE2/3/4/etc assembler; add support for next-gen 256-bit vector instructions.
!
!   (2) Implement hybrid complex floating-point FFT with a modular complex (i.e. Gaussian integer) transform
!       over the Galois field GF(M61^2), where M61 = 2^61 - 1. See Richard Crandall's preprint, "Integer convolution
!       via split-radix fast Galois transform" for details on the fast Galois transform (FGT), freely available at:
!       http://www.perfsci.com/free/techpapers/confgt.ps. On hardware with good 64-bit integer capabilities, this
!       will really rock.
!
!
!   Functionality-related:
!   (*) Multithreading support, especially for SSE2 builds.
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
!  * John Pierce - for hosting the ftp archive on his hogranch.com server.
!
!
!    ...as well as the numerous people who have been kind enough to try the code out
!    and send timings, compiler options for various architectures and suggestions
!    for improvement.
!
!***For UPDATES, PATCHES and USAGE INSTRUCTIONS, see
!
!    http://hogranch.com/mayer/README.html
!
!***Send QUESTIONS, COMMENTS or SUGGESTIONS to me at <ewmayer@aol.com>.
!
!    Alright, let's go hunt for some primes!
!
*/

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
	int		error_checking,
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
	int i,j = 0,nbits;
	/* TODO: some of these need to become 64-bit: */
	uint32 dum,err_iter = 0,findex = 0,ierr = 0,ilo = 0,ihi = 0,iseed,isprime,kblocks = 0,maxiter = 0,n = 0,npad = 0;
	uint64 itmp64;

	/* Exponent of number to be tested - note that for trial-factoring, we represent p
	strictly in string[STR_MAX_LEN] form in this module, only converting it to numeric
	form in the factoring module. For all other types of assignments uint64 should suffice:
	*/
	uint64 p = 0;
	uint32 nbits_in_p = 0;

	/* Res64 and Selfridge-Hurwitz residues: */
	uint64 Res64, Res35m1, Res36m1;
	uint64 sum0, sum1, sum2;

/*...Known Mersenne prime exponents. This array must be null-terminated.	*/

	const uint32 knowns[] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941
		,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593
		,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,0x0};

/*...What a bunch of characters...	*/
	char hex_res[17];

/*...initialize logicals and factoring parameters...	*/
	int restart = FALSE, start_run=TRUE, relax_err = 0;
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
	static int32 nalloc = 0, *arrtmp = 0x0;
	static double *a = 0x0, *a_ptmp = 0x0;
	double final_res_offset;

/*...time-related stuff. clock_t is typically an int (signed 32-bit)
	and a typical value of CLOCKS_PER_SEC (e.g. as defined in <machine/machtime.h>
	on Alpha TruUnix) is 1000000, so we should accumulate clock_t-stored time
	differences at least roughly every half hour (preferably every second or so)
	to avoid weirdness due to flipping of the sign bit or integer overflow.
*/
	/*clock_t clock1, clock2;	Moved these to mers_mod_square.c */
	double tdiff;

	#define SIZE 256
	time_t calendar_time;
	struct tm *local_time;
	char timebuffer[SIZE];

/*...entry point for one or more Lucas-Lehmer tests is here.	*/

	MODULUS_TYPE = mod_type;
	TEST_TYPE = test_type;
	INTERACT = FALSE;

RANGE_BEG:

	/* Clear out any FFT-radix data that might remain from a just-completed run: */
	for(i = 0; i < 10; i++)
	{
		RADIX_VEC[i] = 0;
	}
	NRADICES = 0;

	RESTARTFILE[0] = STATFILE[0] = '\0';
	restart=FALSE;
	relax_err=TRUE;	/* This variable determines whether roundoff error checking may be disabled after	*/
					/* a trial number of error-checked iterations completes with no error warnings.		*/

/*  ...If multithreading enabled, set max. # of threads based on # of available processors, and use
that as the actual # of threads, unless there's a different (and smaller) user-set value of #threads
in an nthreads.ini file : */

#ifdef MULTITHREAD

  #ifdef USE_OMP 
	MAX_THREADS = omp_get_num_procs();
  #elif(defined(USE_PTHREAD))
	MAX_THREADS = get_num_cores();
	ASSERT(HERE, MAX_THREADS > 0, "Illegal #Cores value stored in MAX_THREADS");
  #else
	#error Unrecognized multithreading model!
  #endif
	// MAX_THREADS based on number of processing cores will most often be a power of 2, but don't assume that.
	ASSERT(HERE, MAX_THREADS > 0,"Mlucas.c: MAX_THREADS > 0");

	if(!NTHREADS)	/* User may have already set via -nthread argument, in which case we skip this stuff: */
	{
		if (start_run) fprintf(stderr," looking for number of threads to use in %s file...\n", ntfile);

		fp = fopen(ntfile, "r");
		if ((fp))
		{
			fgets(in_line, STR_MAX_LEN, fp);
			sprintf(cbuf, "Input conversion problem in reading %s file\n", ntfile);
			ASSERT(HERE, sscanf(in_line,"%d", &NTHREADS) == 1, cbuf);
			fclose(fp);	fp = 0x0;

			if(NTHREADS > MAX_THREADS)
			{
				fprintf(stderr,"NTHREADS = %d specified in %s file exceeds number of cores - reducing to %d\n", NTHREADS, ntfile, MAX_THREADS);
				NTHREADS = MAX_THREADS;
			}
			if(start_run) fprintf(stderr,"NTHREADS = %d\n", NTHREADS);
		}
	}
	else {	// In timing-test mode, allow #threads > #cores
		if(NTHREADS > MAX_THREADS)
		{
			fprintf(stderr,"WARN: NTHREADS = %d exceeds number of cores = %d\n", NTHREADS, MAX_THREADS);
		}
		if(start_run) fprintf(stderr,"NTHREADS = %d\n", NTHREADS);
	}

	if(!NTHREADS)
	{
		NTHREADS = MAX_THREADS;
		fprintf(stderr,"Using NTHREADS = #CPUs = %d.\n", NTHREADS);
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

#endif

	/* Make number of iterations between checkpoints dependent on #threads -
	don't want excessively frequent savefile writes, at most 1 or 2 an hour is needed:
	*/
	if(NTHREADS > 4)
		ITERS_BETWEEN_CHECKPOINTS = 100000;
	else
		ITERS_BETWEEN_CHECKPOINTS =  10000;

/*  ...look for worktodo.ini file...	*/

	fp = 0x0;
	if (start_run && (!exponent || (exponent!=0 && fft_length!=0 && iterations==0)))
	{
		fprintf(stderr," looking for %s file...\n", RANGEFILE);
		fp = fopen(RANGEFILE, "r");
	}

	/****************************************/
	/*...Automated exponent-dispatch mode...*/
	/****************************************/
	if (fp)
	{
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
		while(isspace(in_line[i]))
		{
			++i;
		}

		/* If leading char is numeric, assume it's Mersenne-exponent-only (legacy) format: */
		if(isdigit(in_line[i]))
		{
			char_addr += i;
			goto GET_EXPO;
		}

		/* Otherwise assume Prime95-style ini file format, with a possible modulus-specific leading keyword: */
		/* Default is Mersenne, unless overriden by an explicit modulus-type string: */
		MODULUS_TYPE = MODULUS_TYPE_MERSENNE;

		if((char_addr = strstr(in_line, "Fermat")) != 0)
		{
			/* No dicking around with whitespace allowed here: */
			if(char_addr != in_line)
				ASSERT(HERE, 0,"Mlucas.c: modulus type specifier must occur at beginning of worktodo.ini file entry, with no whitespace!");
			else
				char_addr += 6;
			/* Look for comma following the modulus keyword and position next-keyword search right after it: */
			if(!STREQN(char_addr,",",1))
				ASSERT(HERE, 0,"Mlucas.c: expected ',' not found in input following modulus type specifier!");
			else
				char_addr++;

			MODULUS_TYPE = MODULUS_TYPE_FERMAT;
		}
		/* "Mersenne" is the default and hence not required, but allow it: */
		else if((char_addr = strstr(in_line, "Mersenne")) != 0)
		{
			if(char_addr != in_line)
				ASSERT(HERE, 0,"Mlucas.c: modulus type specifier must occur at beginning of worktodo.ini file entry, with no whitespace!");
			else
				char_addr += 8;
			/* Look for comma following the modulus keyword and position next-keyword search right after it: */
			if(!STREQN(char_addr,",",1))
				ASSERT(HERE, 0,"Mlucas.c: expected ',' not found in input following modulus type specifier!");
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
			fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf); fclose(fp);	fp = 0x0; }
			ASSERT(HERE, 0, "");
		}

		char_addr += char_offset;

		ASSERT(HERE, (char_addr = strstr(char_addr, "=")) != 0x0,"Mlucas.c: expected '=' not found in assignment-specifying line!");
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
		ASSERT(HERE, isdigit(*char_addr),"Mlucas.c: isdigit(*char_addr)");
		i = 0;
		while(isdigit(*char_addr))
		{
			ESTRING[i++] = *char_addr;
			++char_addr;
		}
		ESTRING[i++] = '\0';

		p = convert_base10_char_uint64(ESTRING);
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

		ASSERT(HERE, nbits_in_p <= MAX_EXPO_BITS,"Mlucas.c: p <= MAX_EXPO_BITS");

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
				fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
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
					fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
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
					fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					goto GET_NEXT_ASSIGNMENT;
				}
				else if(bit_depth_todo <= bit_depth_done)
				{
					sprintf(cbuf, "ERROR: factor-to bit_depth of %u < already-done depth of %u. The ini file entry was %s\n", bit_depth_done, bit_depth_done, in_line);
										           fprintf(stderr,"%s",cbuf);
					fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
					goto GET_NEXT_ASSIGNMENT;
				}
				else if(bit_depth_todo > log2_max_factor)
				{
					sprintf(cbuf, "WARN: the specified factor-to depth of %u bits exceeds the default %10.4f bits - I hope you know what you're doing.\n", bit_depth_done, log2_max_factor);
					log2_max_factor = bit_depth_todo;
										           fprintf(stderr,"%s",cbuf);
					fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				}
				log2_max_factor = bit_depth_todo;
			}
		}
		else if(nbits_in_p > MAX_PRIMALITY_TEST_BITS)
		{
			sprintf(cbuf, "ERROR: Inputs this large only permitted for trial-factoring. The ini file entry was %s\n", in_line);
								           fprintf(stderr,"%s",cbuf);
			fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
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
						kblocks = get_default_fft_length((uint32)p);	/* Default FFT length for this exponent */
						/* Check that user-specified FFT length is >= default, or that p <= 1.01*(max exponent for given length): */
						if( (kblocks > fft_length) && (p > 1.01*given_N_get_maxP(fft_length<<10)) )
						{
							sprintf(cbuf, "ERROR: Illegal 'fftlen = ' argument - suggested FFT length for this p = %u. The ini file entry was %s\n", kblocks, in_line);
														   fprintf(stderr,"%s",cbuf);
							fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
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
						sprintf(cbuf,"ERROR: bit_depth_done of %u > max. allowed of %u. The ini file entry was %s\n", bit_depth_done, MAX_FACT_BITS, in_line);
													   fprintf(stderr,"%s",cbuf);
						fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
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
						fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
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
								fp = fopen(OFILE,"a"); if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
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
		else	/* Primality test */
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
		ASSERT(HERE, 0, "Illegal combination of command args - please run with -h to see help menu.");
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
	if(TEST_TYPE == TEST_TYPE_TRIALFACTORING)
	{
		/* Do any needed Trial-Factoring and then update the worktodo.ini file appropriately: */
		/* Let the factoring module handle the restart-file processing in this case. */
		if(!ASSIGNMENT_TYPE_MATRIX[MODULUS_TYPE][TEST_TYPE_TRIALFACTORING])
		{
			sprintf(cbuf, "Mlucas.c: TEST_TYPE_TRIALFACTORING with MODULUS_TYPE = %u not supported!\n", MODULUS_TYPE);
			ASSERT(HERE, 0, cbuf);
		}

	#ifndef INCLUDE_TF
		sprintf(cbuf, "Mlucas.c: TEST_TYPE_TRIALFACTORING requires INCLUDE_TF flag to be defined at build time!\n");
		ASSERT(HERE, 0, cbuf);
	#else
		factor(ESTRING, log2_min_factor, log2_max_factor);
		goto GET_NEXT_ASSIGNMENT;
	#endif
	}
	else if(TEST_TYPE != TEST_TYPE_PRIMALITY)
	{
		ASSERT(HERE, 0,"FATAL: Unrecognized assignment type in savefile processing.\n");
	}
	/* endif(TEST_TYPE == ...) */

/********************* Primality Test: ***********************************************/

	if(p < PMIN)
	{
		fprintf(stderr, " p must be at least %u.\n",PMIN);
		return ERR_EXPONENT_ILLEGAL;
	}
	if(p > PMAX)
	{
		fprintf(stderr, " p must be no greater than %u.\n",PMAX);
		return ERR_EXPONENT_ILLEGAL;
	}

	ASSERT(HERE, !(p >> 32), "p must be 32-bit or less!");	/* Future versions will need to loosen this p < 2^32 restriction: */

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		/*  ...make sure p is prime...	*/
		/* TODO: once factoring code is hooked in, make its small-primes table
		a global and init it at the start of execution, then use it to trial-divide
		p up to sqrt(p):
		*/
		/* Allow non-prime-exponents for p-1 test, but not for LL: */
		if(!isPRP((uint32)p))
		{
			fprintf(stderr, " p is not prime!\n");
			return ERR_EXPONENT_ILLEGAL;
		}

		TRANSFORM_TYPE = REAL_WRAPPER;
		sprintf(PSTRING, "M%s", ESTRING);

		maxiter = (uint32)p-2;
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, findex >=14, "Fermat number index must be at least 14!\n");
		ASSERT(HERE, findex < 64, "Fermat number index must be < 64!\n"       );

		/* This takes care of the number-to-char conversion and leading-whitespace-removal
		in one step - use PSTRING for temporary storage here:
		*/
		strcpy(ESTRING, &PSTRING[convert_uint64_base10_char(PSTRING, (uint64)findex)]);
		ASSERT(HERE, (p >> findex) == 1,"Mlucas.c: (p >> findex) == 1");

		TRANSFORM_TYPE = RIGHT_ANGLE;
		sprintf(PSTRING, "F%u", findex);

		maxiter = (uint32)p-1;
	}
	else
	{
		ASSERT(HERE, 0,"Unknown Self-Test Modulus Type!");
	}	/* endif(MODULUS_TYPE) */

	if(!fft_length)
	{
		kblocks = get_default_fft_length((uint32)p);
		if(!kblocks)
		{
			fprintf(stderr,"ERROR detected in get_default_fft_length for p = %u.\n", (uint32)p);
			return ERR_FFTLENGTH_ILLEGAL;
		}
	}
	else
	{
		kblocks = fft_length;
		if((int)kblocks <= 0)
		{
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

		if(error_checking)
			err_iter = maxiter;
		else
			err_iter = 0;
	}
	else
	{
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

		/* See if user has set a particular value of err_iter via the ROCheckIter flag in local.ini: */
		fp = fopen(LOCAL_INI_FILE,"r");
		if(fp)
		{
			while(fgets(in_line, STR_MAX_LEN, fp))		/* Skip over the first block of 3 comment lines, which begin with #	*/
			{
				if(STREQN(in_line, "ROCheckIter",11))
				{
					ASSERT(HERE, (char_addr = strstr(in_line, "=")) != 0x0,"Mlucas.c: expected '=' not found in -specifying line!");
					char_addr++;
					/* Skip any whitespace following the equals sign:*/
					while(isspace(*char_addr))
					{
						++char_addr;
					}
					/* Copy the exponent to cbuf and null-terminate: */
					ASSERT(HERE, isdigit(*char_addr),"Mlucas.c: isdigit(*char_addr)");
					i = 0;
					while(isdigit(*char_addr))
					{
						cbuf[i++] = *char_addr;
						++char_addr;
					}
					cbuf[i++] = '\0';

					/* Check #bits in the power-of-2 exponent vs. the allowed maximum: */
					err_iter = convert_base10_char_uint64(cbuf);
					break;
				}
				else
					continue;
			}
			fclose(fp);	fp = 0x0;
		}
		else
		{
			fprintf(stderr, "INFO: local.ini file %s not found.\n",LOCAL_INI_FILE);
		}
		/* endif(found LOCAL_INI_FILE) */
	}

	/* If specified FFT length smaller than default for this exponent
	[only an issue for user-specified FFT length], print a warning:
	*/
	if(kblocks < (i = get_default_fft_length((uint32)p)))
	{
		fprintf(stderr, "INFO:    maximum recommended exponent for this runlength = %u.\n",given_N_get_maxP(n));
		/* If it's at least close, allow it but print a warning; otherwise error out: */
		if((i << 10) == get_nextlarger_fft_length(n))
		{
			fprintf(stderr, "specified FFT length %d K is less than recommended %d K for this p.\n",kblocks,i);
		}
		else
		{
			sprintf(cbuf, "specified FFT length %d K is much too small: Recommended length for this p = %d K.\n",kblocks,i);
			ASSERT(HERE, 0, cbuf);
		}
	}

	/*...calculate the unpadded runlength...	*/
	n = kblocks << 10;

	if(kblocks != (n>>10))
	{
		fprintf(stderr, "ERROR: length %d K overflows N = 1024*K.\n",kblocks);
		return ERR_FFTLENGTH_ILLEGAL;
	}

	/*...Set the array padding parameters - only use array padding elements for runlengths > 32K. */
	if(kblocks > 32)
	{
		DAT_BITS = DAT_BITS_DEF;
		PAD_BITS = PAD_BITS_DEF;
		/*...If array padding turned on, check that the blocklength divides the unpadded runlength...	*/
		if((DAT_BITS < 31) && ((n >> DAT_BITS) << DAT_BITS) != n)
		{
			ASSERT(HERE, 0,"ERROR: blocklength does not divide runlength!");
		}
	}
	else
	{
		DAT_BITS =31;	/* This causes the padding to go away */
		PAD_BITS = 0;
	}

	/*...Find padded array length...	*/
	npad = n + ( (n >> DAT_BITS) << PAD_BITS );	/* length of padded data array.	*/
	/* If the residue data array is too small for the new assignment, deallocate it: */
	if(nalloc > 0 && npad > nalloc)
	{
		ASSERT(HERE, a != 0x0 && a_ptmp != 0x0,"Mlucas.c: a != 0x0 && a_ptmp != 0x0");
		free((void *)a_ptmp); a_ptmp=0x0; a=0x0;
		free((void *)arrtmp); arrtmp=0x0;
		nalloc = 0;
	}

	if(nalloc == 0)
	{
		/* If it's a multi-exponent self-test, alloc to the maximum FFT length
		that will be used: */
		if(maxFFT > kblocks)
		{
			i = (maxFFT << 10);
			npad = i + ( (i >> DAT_BITS) << PAD_BITS );
		}
		nalloc = npad;
		ASSERT(HERE, a == 0x0 && a_ptmp == 0x0 && arrtmp == 0x0,"Mlucas.c: a == 0x0 && a_ptmp == 0x0 && arrtmp == 0x0");
		a_ptmp = ALLOC_DOUBLE(a_ptmp, nalloc);	if(!a_ptmp){ sprintf(cbuf, "FATAL: unable to allocate array A in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		a      = ALIGN_DOUBLE(a_ptmp);
		ASSERT(HERE, ((long)((void *)a) & 63) == 0x0,"Mlucas.c: a[] not aligned on 64-byte boundary!");
		if(((long)((void *)a) & 127) != 0x0)fprintf(stderr, "WARN: Mlucas.c: a[] = 0x%08X not aligned on 128-byte boundary!\n", (uint32)((void *)a));

		/* Alloc this to the old-style nalloc*int size, so it'll be large
		enough to hold either an old-style integer residue or a (more compact)
		bytewise residue:
		*/
		arrtmp = (int*)malloc(nalloc* sizeof(int));	if(!arrtmp){ sprintf(cbuf, "FATAL: unable to allocate array ARRTMP in main.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
	}

// Multithreaded-code debug: Set address to watch:
#ifdef MULTITHREAD
	ADDR0 = a;
#endif

	/* Make sure we start with primary restart file: */
	RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');

READ_RESTART_FILE:

	if(!INTERACT)	/* Skip the restart-file stuff if it's a self-test */
	{
		/* See if there's a restart file: */
		fp = fopen(RESTARTFILE, "rb");
		/* If so, read the savefile: */
		if(fp)
		{
			i = read_ppm1_savefiles (p, fp, &ilo, (uint8*)arrtmp, &Res64, &Res35m1, &Res36m1);
			fclose(fp); fp = 0x0;

			if(!i)	/*	&& !convert_LL_savefiles(p, fp, &ilo, n, (int32*)arrtmp, a))	*/
			{
				fp = fopen(   OFILE,"a");
				fq = fopen(STATFILE,"a");

				/* First print any error message that may have been issued during the above function call: */
				if(strstr(cbuf, "read_ppm1_savefiles"))
				{
						fprintf(stderr,"%s",cbuf);
				if(fp){ fprintf(	fp,"%s",cbuf); }
				if(fq){ fprintf(	fq,"%s",cbuf); }
				}
				/* And now for the official spokesmessage: */
				sprintf(cbuf, "ERROR: read_ppm1_savefiles Failed on savefile %s!\n",RESTARTFILE);
						fprintf(stderr,"%s",cbuf);
				if(fp){ fprintf(	fp,"%s",cbuf); fclose(fp); fp = 0x0; }
				if(fq){ fprintf(	fq,"%s",cbuf); fclose(fq); fq = 0x0; }

				if(RESTARTFILE[0] != 'q')
				{
					RESTARTFILE[0] = 'q';
					goto READ_RESTART_FILE;
				}
				else
					ASSERT(HERE, 0,"Mlucas.c: Failed to correctly read both primary or secondary savefile!");
			}
			ASSERT(HERE, ilo > 0,"Mlucas.c: ilo <= 0!");
			ASSERT(HERE, ilo < maxiter,"Mlucas.c: ilo >= maxiter!");

			/* Allocate the floating-point residue array and convert the
			savefile bytewise residue to the appropriate floating-point form.
			Checksum validation is performed within the convert_res_bytewise_FP function:
			*/
			if(!convert_res_bytewise_FP((uint8*)arrtmp, a, n, p, Res64, Res35m1, Res36m1))
			{
				sprintf(cbuf, "ERROR: convert_res_bytewise_FP Failed on savefile %s!\n",RESTARTFILE);
													fprintf(stderr,"%s",cbuf);
				fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				fp = fopen(STATFILE,"a");	if(fq){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				if(RESTARTFILE[0] != 'q')
				{
					RESTARTFILE[0] = 'q';
					goto READ_RESTART_FILE;
				}
				else
					ASSERT(HERE, 0,"0");
			}

			restart = TRUE;
			ihi = ilo+ITERS_BETWEEN_CHECKPOINTS;
			/* If for some reason last checkpoint was at a non-multiple of ITERS_BETWEEN_CHECKPOINTS, round down: */
			ihi-= ihi%ITERS_BETWEEN_CHECKPOINTS;
		}
		else /* if(!fp) */
		{
			/* If we're on the primary restart file, set up for secondary: */
			if(RESTARTFILE[0] != 'q')
			{
				sprintf(cbuf, "INFO: primary restart file %s not found...looking for secondary...\n",RESTARTFILE);
													fprintf(stderr,"%s",cbuf);
				fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				fp = fopen(STATFILE,"a");	if(fq){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }

				RESTARTFILE[0] = 'q';
				goto READ_RESTART_FILE;
			}
			else
			{
				sprintf(cbuf, "INFO: no restart file found...starting run from scratch.\n");
													fprintf(stderr,"%s",cbuf);
				fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
				fp = fopen(STATFILE,"a");	if(fq){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }

				RESTARTFILE[0] = 'p';
			}
		}	/* endif(fp) */
	}	/* endif(!INTERACT)	*/

	/* If a restart, parse the .stat file to see if any RO warnings were previously issued.	*/
	/* If so, will require RO error checking to be performed on each iteration, even if iter > err_iter.	*/
	if(restart)
	{
		fp = fopen(STATFILE,"r");
		if(fp)
		{
			for(;;)
			{
				if(fgets(in_line, STR_MAX_LEN, fp))
				{
					if(strstr(in_line,"Roundoff"))	/* If any RO warnings were issued during initial 200K iterations, leave error checking on... */
					{
						if(!strstr(in_line,"0.40625"))	/* ...but only if the RO errors were > 0.40625 . */
						{
							fprintf(stderr, "INFO: RO warnings > 0.40625 detected in first 200K iterations;\n");
							fprintf(stderr, "Leaving error checking on for rest of the LL test of this exponent.\n");
							relax_err = FALSE;
							break;
						}
					}
				}
				else
				{
					break;
				}
			}
			fclose(fp);	fp = 0x0;
		}
	}
	else
	{
		/*...set initial iteration loop parameters...	*/
		ilo=0;
		if(INTERACT)
			ihi=timing_test_iters;
		else
			ihi=ITERS_BETWEEN_CHECKPOINTS;
	}

	/*...If a restart and RO warnings were previously issued to the .stat file, require error checking
	to be done on each iteration, even if iter > err_iter. We implement this simply by making err_iter > maxiter.	*/
	if(!relax_err)
	{
		err_iter = maxiter + 1;
	}

	/* If it's not a timing test, make sure that no matter what, we do error checking for at least the first few 200k iter */
	if(!INTERACT && err_iter < 200000)
	{
		err_iter = 200000;
		fprintf(stderr, "INFO: Using default value of ROCheckIter = %u.\n", err_iter);
	}

	ASSERT(HERE,MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(HERE,MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	/* If at the start of a p-1 or primality test, set the initial seed for the run: */
	ASSERT(HERE, TEST_TYPE == TEST_TYPE_PMINUS1 || TEST_TYPE == TEST_TYPE_PRIMALITY,"Mlucas.c: TEST_TYPE == TEST_TYPE_PMINUS1 || TEST_TYPE == TEST_TYPE_PRIMALITY");
	if(ilo == 0)
	{
		/* Always use 3 as the p-1 and Fermat-test seed, and 4 for the LL-test seed: */
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
			iseed = 4;
		else
			iseed = 3;

		memset(a, 0, npad*sizeof(double));
		a[0] = iseed;
	}

	if(restart)
	{
		itmp64 = res64(a, n, p,&nbits,hex_res);
		ASSERT(HERE,Res64 == itmp64,"On restart: Res64 != itmp64");
		sprintf(cbuf, "Restarting %s at iteration = %u. Res%2d: %s\n",PSTRING,ilo,nbits,hex_res);
	}

	/*...Restart and FFT info.	*/
	if(INTERACT)
	{
		fprintf(stderr, "%s: using FFT length %uK = %u 8-byte floats.\n",PSTRING,kblocks,n);
		fprintf(stderr, " this gives an average %20.15f bits per digit\n",1.0*p/n);
	}
	else
	{
		fp = fopen(STATFILE,"a");
		if(fp)
		{
			fprintf(fp,"%s",cbuf);
			fprintf(fp,"%s: using FFT length %uK = %u 8-byte floats.\n",PSTRING,kblocks,n);
			fprintf(fp," this gives an average %20.15f bits per digit\n",1.0*p/n);
			fclose(fp);	fp = 0x0;
		}
	}

/*...main loop...	*/

	for(;;)
	{
		ASSERT(HERE, (int)maxiter > 0,"Mlucas.c: (int)maxiter > 0");
		if(ihi > maxiter) ihi=maxiter;

		/* Here's the big one - (ITERS_BETWEEN_CHECKPOINTS) squaring steps.
		XYZ_mod_square returns 0 if no errors detected during this iteration cycle.

		If fatal error was encountered, skip to next assignment in worktodo.ini file
		but keep restart files for current assignment around (can finish using hiacc code.)
		*/
		AME = MME = 0.0;	/* Init Avg. & Max. RO Error */

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			ierr = mers_mod_square  (a, arrtmp, n, ilo, ihi, p, &err_iter, scrnFlag, &tdiff);
			if(ierr)
				goto GET_NEXT_ASSIGNMENT;
		}
		else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			ierr = fermat_mod_square(a, arrtmp, n, ilo, ihi, p, &err_iter, scrnFlag, &tdiff);
			if(ierr)
				goto GET_NEXT_ASSIGNMENT;
		}

		/*...Done?	*/
		if(INTERACT) break;	/* C version uses zero-offset iteration counter. */

		AME /= ITERS_BETWEEN_CHECKPOINTS;

		/*...Every (ITERS_BETWEEN_CHECKPOINTS)th iteration, print timings to stdout or STATFILE.
		If it's a regular (i.e. non-timing) test, also write current residue to restart files.
		*/

		/*...get a quick mod-2^64 residue...	*/
		Res64 = res64(a, n, p,&nbits,hex_res);
				resSH(a, n, p, &Res35m1, &Res36m1);

		/*...get a quick timestamp...	*/
		calendar_time = time(NULL);
		local_time = localtime(&calendar_time);
		strftime(timebuffer,SIZE,"%b %d %H:%M:%S",local_time);

		/*...print runtime in hh:mm:ss format.	*/
		sprintf(cbuf, "[%s] %s Iter# = %u clocks =%s [%8.4f sec/iter] Res%2d: %s. AvgMaxErr = %10.9f. MaxErr = %10.9f\n"
			,timebuffer,PSTRING,ihi,get_time_str(tdiff)
		,tdiff/(ihi - ilo),nbits,hex_res, AME, MME);

		if(INTERACT)
			fprintf(stderr,"%s",cbuf);
		else
		{
			fp = fopen(STATFILE,"a");	  if(fp){	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;	}
			if (scrnFlag)	fprintf(stderr,"%s",cbuf);	/* Echo output to stddev */
		}
		// For Mersennes we do not save a final residue (but this leaves the penultimate residue file intact)
		if( (MODULUS_TYPE == MODULUS_TYPE_MERSENNE) && (ihi == maxiter) ) break;

		/* We've already calculated the SH residues here, but recalculating them during floating-to-bytewise
		residue array conversion is sufficiently cheap that the redundant calls to res64() and resSH() above don't
		really matter - use as an extra check on convert_res_FP_bytewise():
		*/
		sum0 = Res64; sum1 = Res35m1; sum2 = Res36m1;

		convert_res_FP_bytewise(a, (uint8*)arrtmp, n, p, &Res64, &Res35m1, &Res36m1);
		ASSERT(HERE, sum0 == Res64  , "Res64   returned by convert_res_FP_bytewise differs from res64()!");
		ASSERT(HERE, sum1 == Res35m1, "Res35m1 returned by convert_res_FP_bytewise differs from resSH()!");
		ASSERT(HERE, sum2 == Res36m1, "Res36m1 returned by convert_res_FP_bytewise differs from resSH()!");

		/* Make sure we start with primary restart file: */
		RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');
	WRITE_RESTART_FILE:
		fp = fopen(RESTARTFILE, "wb");
		if(fp)
		{
		#if OLD_FORMAT
			write_ppm1_savefiles(p, fp, ihi, n, (int32*)arrtmp, a);
		#else
			write_ppm1_savefiles(p, fp, ihi, (uint8*)arrtmp, Res64, Res35m1, Res36m1);
		#endif
			fclose(fp);	fp = 0x0;

			/* If we're on the primary restart file, set up for secondary: */
			if(RESTARTFILE[0] != 'q')
			{
				RESTARTFILE[0] = 'q';
				goto WRITE_RESTART_FILE;
			}
			else
				RESTARTFILE[0] = ((MODULUS_TYPE == MODULUS_TYPE_MERSENNE) ? 'p' : 'f');

		}
		else
		{
			sprintf(cbuf, "ERROR: unable to open restart file %s for write of checkpoint data.\n",RESTARTFILE);
										        fprintf(stderr,"%s",cbuf);
			fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
			fp = fopen(STATFILE,"a");	if(fq){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
			/*
			Don't want to assert here - asllow processing to continue, in case this is a transient failure-to-open.
			e.g. due to a backup utility having temporarily locked the savefile, or an out-of-disk-space problem.
			*/
		}	/* endif(fp) */

		// For Fermats, exit only after writing final-residue checkpoint file:
		if( (MODULUS_TYPE == MODULUS_TYPE_FERMAT) && (ihi == maxiter) ) break;

		/*...reset loop parameters and begin next iteration cycle...	*/
		ilo=ihi;
		ihi=ilo+ITERS_BETWEEN_CHECKPOINTS;
	}

	/*...For timing tests, print timing and 64-bit residue and return timing via arglist runtime ptr...	*/

	*runtime = tdiff;
	Res64 = res64(a, n, p,&nbits,hex_res);
			resSH(a, n, p, &Res35m1, &Res36m1);

	/* If Selftest mode... */
	if(INTERACT)
	{
		fprintf(stderr, "%u iterations of %s with FFT length %u = %u K\n",timing_test_iters,PSTRING,n,kblocks);

		/* If Fermat number, make sure exponent a power of 2: */
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			ASSERT(HERE, (p >> findex) == 1,"Mlucas.c: (p >> findex) == 1");
		}

		if(nbits < 64) fprintf(stderr, "WARNING: this residue contains only %u bits\n",nbits);

		if(timing_test_iters > AME_ITER_START)
		{
			AME /= (timing_test_iters - AME_ITER_START);
			fprintf(stderr, "Res%2d: %s. AvgMaxErr = %10.9f. MaxErr = %10.9f. Program: E%s\n", nbits, hex_res, AME, MME, VERSION);
		}
		else
		{
			fprintf(stderr, "Res%2d: %s. AvgMaxErr N/A. MaxErr = %10.9f. Program: E%s\n", nbits, hex_res, MME, VERSION);
		}
		/* MSVC/.NET incorrectly output these when using uint64 and %20llu format, so cast to double and print: */
		fprintf(stderr, "Res mod 2^36     = %20.0f\n",(double)(Res64 & (uint64)0x0000000FFFFFFFFFull));
		fprintf(stderr, "Res mod 2^35 - 1 = %20.0f\n",(double)Res35m1);
		fprintf(stderr, "Res mod 2^36 - 1 = %20.0f\n",(double)Res36m1);

		/* If they are provided, check the Selfridge-Hurwitz residues: */
		if(sh0)
		{
			if(*sh0 != 0)
			{
				if (Res64 == *sh0)
					resFlag = 0;
				else
				{
					resFlag = 1;	/* False */
					fprintf(stderr, "  ***   Res64 Error   ***\n");
					fprintf(stderr, " current   = %20.0f\n", (double)Res64);
					fprintf(stderr, " should be = %20.0f\n", (double) *sh0);
				}
			}
			else
				*sh0 = Res64;
		}
		if(sh1)
		{
			if(*sh1 != 0)
			{
				if (Res35m1 == *sh1)
					resFlag = 0;
				else
				{
					resFlag = 1;	/* False */
					fprintf(stderr, "  ***   Res35m1 Error   ***\n");
					fprintf(stderr, " current   = %20.0f\n", (double)Res35m1);
					fprintf(stderr, " should be = %20.0f\n", (double) *sh1);
				}
			}
			else
				*sh1 = Res35m1;
		}
		if(sh2)
		{
			if(*sh2 != 0)
			{
				if (Res36m1 == *sh2)
					resFlag = 0;
				else
				{
					resFlag = 1;	/* False */
					fprintf(stderr, "  ***   Res36m1 Error   ***\n");
					fprintf(stderr, " current   = %20.0f\n", (double)Res36m1);
					fprintf(stderr, " should be = %20.0f\n", (double) *sh2);
				}
			}
			else
				*sh2 = Res36m1;
		}

		/*...print runtime in hh:mm:ss format.	*/
		fprintf(stderr, "Clocks =%s\n",get_time_str(tdiff) );

		/*exit(EXIT_SUCCESS);*/
 		return(resFlag);

	}	/* endif(INTERACT) */

	/*...Check primality:	*/
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		isprime = 1;

		/* For Fermat numbers, in balanced-digit form, it's prime if the lowest-order digit = -1, all others 0: */
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
			final_res_offset = 1;
		else
			final_res_offset = 0;

		a[0] += final_res_offset;
		for(i=0; i < n; i++)
		{
			j = i + ( (i >> DAT_BITS) << PAD_BITS );
			if(a[j] != 0.0) { isprime = 0; break; }
		}
		a[0] -= final_res_offset;

		/*...Unbelievable - I must be hallucinating.	*/
		if(isprime)
		{
			if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
			{
				/*... this gets written both to file and to stdout, the latter irrespective of whether the run is in interactive mode...	*/
				sprintf(cbuf, "%s is a new FERMAT PRIME!!!\nPlease send e-mail to ewmayer@aol.com and crandall@reed.edu.\n",PSTRING);
											        fprintf(stderr,"%s",cbuf);
				fp = fopen(STATFILE,"a");	if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp);	fp = 0x0; }
				fp = fopen(   OFILE,"a");	if(fq){ fprintf(    fp,"%s",cbuf);	fclose(fp);	fp = 0x0; }
			}
			else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				for(i=0; knowns[i] != 0; i++)
				{
					if(p == knowns[i]) { break; }
				}

				if(knowns[i] != 0)
				{
					sprintf(cbuf, "%s is a known MERSENNE PRIME.\n",PSTRING);

					if(INTERACT || scrnFlag)	/* Echo output to stddev */
						fprintf(stderr,"%s",cbuf);
					else
					{
						fp = fopen(STATFILE,"a");	if(fp){ fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0; }
						fp = fopen(   OFILE,"a");	if(fq){ fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0; }
					}
				}
				else
				{
					/*... this gets written both to file and to stdout, the latter irrespective of whether the run is in interactive mode...	*/
					sprintf(cbuf, "%s is a new MERSENNE PRIME!!!\nPlease send e-mail to ewmayer@aol.com and woltman@alum.mit.edu.\n",PSTRING);
												fprintf(stderr,"%s",cbuf);
					fp = fopen(STATFILE,"a");	if(fp){ fprintf(    fp,"%s",cbuf);	fclose(fp);	fp = 0x0; }
					fp = fopen(   OFILE,"a");	if(fq){ fprintf(    fp,"%s",cbuf);	fclose(fp);	fp = 0x0; }
				}
			}
			else
				ASSERT(HERE, 0, "Unsupported modulus type!");
		}
		/*
		The more likely scenario - it's not prime, so we form a 64-bit residue and write that.
		If residue has < 64 bits, print a warning.
		*/
		else
		{
			/*...otherwise, write the 64-bit hex residue.	*/
			fp = fopen(STATFILE,"a");
			fq = fopen(   OFILE,"a");
			sprintf(cbuf, "%s is not prime. Res%2d: %s. Program: E%s\n",PSTRING,nbits,hex_res,VERSION);
			if(fp){ fprintf(fp,"%s",cbuf); }
			if(fq){ fprintf(fq,"%s",cbuf); }
			if(INTERACT)
			{
				if(nbits < 64) sprintf(cbuf, "WARNING: this residue contains only %u bits\n",nbits);
				fprintf(stderr,"%s",cbuf);
			}
			else
			{
				if(nbits < 64)
				{
					if (scrnFlag)		/* Echo output to stddev */
					{
						fprintf(stderr,"%s",cbuf);
					}
					if(fp){ fprintf(fp,"WARNING: this residue contains only %u bits\n",nbits); }
					if(fq){ fprintf(fq,"WARNING: this residue contains only %u bits\n",nbits); }
				}

				sprintf(cbuf, "%s mod 2^36     = %20.0f\n",PSTRING,(double)(Res64 & (uint64)0x0000000FFFFFFFFFull));
				if(fp){ fprintf(fp,"%s",cbuf); }
				if(fq){ fprintf(fq,"%s",cbuf); }
				sprintf(cbuf, "%s mod 2^35 - 1 = %20.0f\n",PSTRING,(double)Res35m1);
				if(fp){ fprintf(fp,"%s",cbuf); }
				if(fq){ fprintf(fq,"%s",cbuf); }
				sprintf(cbuf, "%s mod 2^36 - 1 = %20.0f\n",PSTRING,(double)Res36m1);
				if(fp){ fprintf(fp,"%s",cbuf); }
				if(fq){ fprintf(fq,"%s",cbuf); }

			}
			if(fp){ fclose(fp);	fp = 0x0; }
			if(fq){ fclose(fq);	fq = 0x0; }
		}
	}	/* endif(TEST_TYPE == TEST_TYPE_PRIMALITY) */

	/*...If successful completion, delete the secondary restart files...save the primary in case it's a prime,
	so we can re-run the last p%2000 iterations by way of quick partial verification:
	*/
#if 0
	if(remove(RESTARTFILE)) fprintf(stderr, "Unable to delete primary restart file %s\n",RESTARTFILE);
#endif
	RESTARTFILE[0] = 'q';
	if(remove(RESTARTFILE)) fprintf(stderr, "Unable to delete secondary restart file %s\n",RESTARTFILE);
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
		fp = fopen(RANGEFILE,"r");
		if(!fp)
		{
			sprintf(cbuf,"ERROR: unable to open %s file for reading.\n",RANGEFILE);
			fp = fopen(OFILE,"a");
			if(fp)
			{
				fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
			}
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		/* Remove any WINI.TMP file that may be present: */
		remove("WINI.TMP");
		fq = fopen("WINI.TMP", "w");
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
			fp = fopen(OFILE, "a");	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
			ASSERT(HERE,0,"");
		}
		else if(!strstr(in_line, ESTRING))
		{
			sprintf(cbuf, "WARNING: Current exponent %s not found in line 1 of %s file - skipping line deletion\n", ESTRING, RANGEFILE);
			fprintf(stderr,"%s",cbuf);
			fp = fopen(OFILE, "a");	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
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
					ASSERT(HERE, char_addr != 0x0,"Mlucas.c: char_addr");
					sprintf(++char_addr, "%u", bit_depth_done);
					fputs(in_line, fq);
				}
				else if(TEST_TYPE == TEST_TYPE_PMINUS1)
				{
					/*0/1 flag indicating whether P-1 has been done assumed to follow second comma in in_line: */
					char_addr = strstr(char_addr, ",");
					ASSERT(HERE, char_addr != 0x0,"Mlucas.c: char_addr");
					char_addr++;
					char_addr = strstr(char_addr, ",");
					ASSERT(HERE, char_addr != 0x0,"Mlucas.c: char_addr");
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
			fp = fopen(OFILE,"a");	if(fp){	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;	}

			/* If attempting to simply rename the TMP file fails, do it the hard way: */
			fp = fopen(RANGEFILE,"w");
			if(!fp)
			{
				sprintf(cbuf,"ERROR: unable to open %s file for writing.\n", RANGEFILE);
				fp = fopen(OFILE,"a");
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

			fq = fopen("WINI.TMP", "r");
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

		/* if one or more exponents left in rangefile, go back for more... */

		if (i > 0) goto RANGE_BEG;   /* CYCLE RANGE */
	}

/* If it was an interactive run or the last exponent in the rangefile, we're done. */

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
	ASSERT(HERE, get_fft_radices(MIN_FFT_LENGTH_IN_K, 0, 0x0, 0x0, 0) == 0,"Mlucas.c: get_fft_radices(MIN_FFT_LENGTH_IN_K, 0) == 0");
	n = (MIN_FFT_LENGTH_IN_K << 10);
	/* Make sure N didn't overflow */
	ASSERT(HERE, (n >> 10) == MIN_FFT_LENGTH_IN_K,"Mlucas.c: (n >> 10) == MIN_FFT_LENGTH_IN_K");
	PMIN = 2*n;	/* 2 bits per input is about the smallest we can test without getting nonzero-carry errors */

	/* Set max. exponent (in terms of power of 2) that can be tested: */
	/* Check that the purported max. FFT length is actually supported: */
	ASSERT(HERE, get_fft_radices(MAX_FFT_LENGTH_IN_K, 0, 0x0, 0x0, 0) == 0,"Mlucas.c: get_fft_radices(MAX_FFT_LENGTH_IN_K, 0) == 0");
	n = (MAX_FFT_LENGTH_IN_K << 10);
	/* Make sure N didn't overflow */
	ASSERT(HERE, (n >> 10) == MAX_FFT_LENGTH_IN_K,"Mlucas.c: (n >> 10) == MAX_FFT_LENGTH_IN_K");
	PMAX = given_N_get_maxP(n);

	ASSERT(HERE, PMAX > PMIN,"Mlucas.c: PMAX > PMIN");

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

/*
Function to take an n-digit Mersenne or Fermat-mod residue, convert it to
positive-digit form, and pretty-print the lowest-order 16 hex digits.
*/
uint64 	res64(double a[], int n, const uint64 p, int *nbits, char *hex_res)
{
	int bimodn,cy,findex,ii,j,j1;
	int bw,sw,base[2],bits[2];
	int64  itmp;
	uint64 retval = 0, res64;
	int pow2_fft;

	ASSERT(HERE,MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(HERE,MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	ASSERT(HERE,TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(HERE,TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");

	/* If Fermat number, make sure exponent a power of 2: */
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, TRANSFORM_TYPE == RIGHT_ANGLE,"Mlucas.c: TRANSFORM_TYPE == RIGHT_ANGLE");
		findex = trailz64(p);
		ASSERT(HERE, (p >> findex) == 1,"Mlucas.c: (p >> findex) == 1");
	}
	else
		ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER,"Mlucas.c: TRANSFORM_TYPE == REAL_WRAPPER");

	/* Vector length a power of 2? */
	if((n >> trailz32(n)) == 1)
	{
		pow2_fft = TRUE;
	}
	else
	{
		pow2_fft = FALSE;
	}

	hex_res[0]='\0';
	hex_res[16]='\0';
	bits[0] = p/n;
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
	/*
	If most-significant digit in the balanced-representation form is < 0, add the modulus to the residue.
	For Mersenne (2^p-1) and Fermat (2^p+1) moduli, can combine this with the normalize-to-nonnegative-digit
	step (which we do in any event) by simply initializing the carry into the latter to -1 or +1, respectively:
	*/
	cy=0;		/* init carry.	*/

	for(j=n-1; j >= 0; j -= TRANSFORM_TYPE)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

		if(a[j1] != 0.0)
		{
			if(a[j1] < 0.0)
			{
				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					cy = -1;
				}
				else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
				{
					cy = +1;
				}
			}
			break;
		}
	}

	/*...Now form the hexadecimal residue, converting digits to positive-digit form along the way...	*/

	*nbits = 0;	/* Total bits accumulated so far in the hex residue	*/

	if(TRANSFORM_TYPE == REAL_WRAPPER)
	{
		bimodn = 0;
		ii = 0;		/* Index into the BASE array. */
		/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
		if(bw > 0)
			ii = 1;

		for(j=0; j < n; j++)
		{
		#ifdef USE_AVX
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			itmp = (int64)(a[j1] + cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0)			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
			{
				itmp += (base[ii]);
				cy = -1;
			}
			else
			{
				cy = 0;
			}

			retval += (itmp << *nbits);
			*nbits += bits[ii];
			if(*nbits >= 64){ *nbits = 64; break; }

			bimodn += bw;
			if(bimodn >= n) bimodn -= n;
			ii = (uint32)(sw - bimodn) >> 31;
		}
	}
	else
	{
		bimodn = n;

		for(j=0; j < n; j += 2)
		{
		#ifdef USE_AVX
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

			itmp = (int64)(a[j1]+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0)			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
			{
				itmp += (base[ii]);
				cy = -1;
			}
			else
			{
				cy = 0;
			}

			retval += (itmp << *nbits);
			*nbits += bits[ii];
			if(*nbits >= 64){ *nbits = 64; break; }
		}
	}

	res64 = retval;
	for(ii=15; ii >= 16-(*nbits>>2); ii--)
	{
		hex_res[ii] = hex_chars[res64 & 15];
		res64 >>= 4;
	}
	for(ii=15-(*nbits>>2); ii >= 0; ii--)
	{
		hex_res[ii] = 'x';
	}

	return retval;
}

/********************/

/*
Function to take an n-digit Mersenne or Fermat-mod residue, convert it to
positive-digit form, and return the Selfridge-Hurwitz residues,
i.e. the final residue mod 2^35-1 and 2^36-1.

In the Mersenne-mod case it is assumed the residue digits are stored
consecutively in the a[] array.

In the Fermat-mod case we assume the digits are arranged in (j,j+n/2)
(i.e. right-angle transform) order.
*/

void 	resSH(double a[], int n, const uint64 p, uint64 *Res35m1, uint64 *Res36m1)
{
	uint64 nbits;
	int bimodn,cy,findex,ii,j,j1,pass,shift;
	int bw,sw,bits[2];
	uint64 base[2];	/* Assume base may be > 2^32 (e.g. for mixed FFT/FGT) but < 2^53, i.e. fits in a double */
int bs_count[2];
	int64 itmp;
	uint64 curr_word = 0, mod1=0, mod2=0;
	const uint64 two35m1 = (uint64)0x00000007FFFFFFFFull, two36m1 = (uint64)0x0000000FFFFFFFFFull;	/* 2^35,36-1 */
	int pow2_fft;

	ASSERT(HERE,MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(HERE,MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	ASSERT(HERE,TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(HERE,TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");

	/* If Fermat number, make sure exponent a power of 2: */
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, TRANSFORM_TYPE == RIGHT_ANGLE,"Mlucas.c: TRANSFORM_TYPE == RIGHT_ANGLE");
		findex = trailz64(p);
		ASSERT(HERE, (p >> findex) == 1,"Mlucas.c: (p >> findex) == 1");
	}
	else
		ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER,"Mlucas.c: TRANSFORM_TYPE == REAL_WRAPPER");

	/* Vector length a power of 2? */
	if((n >> trailz32(n)) == 1)
	{
		pow2_fft = TRUE;
	}
	else
	{
		pow2_fft = FALSE;
	}

	bits[0] = p/n;
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
	/*
	If most-significant digit in the balanced-representation form is < 0, add the modulus to the residue.
	For Mersenne (2^p-1) and Fermat (2^p+1) moduli, can combine this with the normalize-to-nonnegative-digit
	step (which we do in any event) by simply initializing the carry into the latter to -1 or +1, respectively:
	*/
	cy=0;		/* init carry.	*/

	for(j=n-1; j >= 0; j -= TRANSFORM_TYPE)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

		if(a[j1]!= 0.0)
		{
			if(a[j1]< 0.0)
			{
				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					cy = -1;
				}
				else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
				{
					cy = +1;
				}
			}
			break;
		}
	}

	/*...Now form the SH residues, converting to positive-digit form along the way...	*/

	nbits = 0;		/* Total bits accumulated so far in the residue	*/

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
		#ifdef USE_AVX
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			itmp = (int64)(a[j1]+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0)			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
			{
				itmp += (base[ii]);
				cy = -1;
			}
			else
			{
				cy = 0;
			}

			ASSERT(HERE, itmp >= 0,"Mlucas.c: itmp >= 0");

		/* Mod-(2^35-1) residue: */
			curr_word = (uint64)itmp;
			/* Current word must be left-shifted by nbits%35 before folding into residue: */
			shift = (nbits%35);
			mod1 += (curr_word << shift) & two35m1;
			mod1 = (mod1 >> 35) + (mod1 & two35m1);
			curr_word >>= (35-shift);
			while(curr_word)
			{
				mod1 += curr_word & two35m1;
				mod1 = (mod1 >> 35) + (mod1 & two35m1);
				curr_word >>= 35;
			}
		/* Mod-(2^36-1) residue: */
			curr_word = (uint64)itmp;
			/* Current word must be left-shifted by nbits%36 before folding into residue: */
			shift = (nbits%36);
			mod2 += (curr_word << shift) & two36m1;
			mod2 = (mod2 >> 36) + (mod2 & two36m1);
			curr_word >>= (36-shift);
			while(curr_word)
			{
				mod2 += curr_word & two36m1;
				mod2 = (mod2 >> 36) + (mod2 & two36m1);
				curr_word >>= 36;
			}

			nbits += bits[ii];

			bimodn += bw;
			if(bimodn >= n) bimodn -= n;
			ii = (uint32)(sw - bimodn) >> 31;
		}
	}
	else
	{
bs_count[0] = bs_count[1] = 0;
	  for(pass = 0; pass <=1; pass++)
	  {
		bimodn = n;
		for(j = pass; j < n; j += 2)
		{
		#ifdef USE_AVX
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			ii = (bimodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */
++bs_count[ii];
			bimodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */
			bimodn += ( ((int)ii-1) & n);		/*       add 0 if a bigword,   N if a smallword */

			itmp = (int64)(a[j1]+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0)			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
			{
				itmp += (base[ii]);
				cy = -1;
			}
			else
			{
				cy = 0;
			}
			if(itmp < 0)
			{
				fprintf(stderr, "Warning: itmp < 0 detected: value = %lld\n", itmp);
				ASSERT(HERE, itmp >= 0,"Mlucas.c: itmp >= 0");
			}
		/* Mod-(2^35-1) residue: */
			curr_word = (uint64)itmp;
			/* Current word must be left-shifted by nbits%35 before folding into residue: */
			shift = (nbits%35);
			mod1 += (curr_word << shift) & two35m1;
			mod1 = (mod1 >> 35) + (mod1 & two35m1);
			curr_word >>= (35-shift);
			while(curr_word)
			{
				mod1 += curr_word & two35m1;
				mod1 = (mod1 >> 35) + (mod1 & two35m1);
				curr_word >>= 35;
			}
		/* Mod-(2^36-1) residue: */
			curr_word = (uint64)itmp;
			/* Current word must be left-shifted by nbits%36 before folding into residue: */
			shift = (nbits%36);
			mod2 += (curr_word << shift) & two36m1;
			mod2 = (mod2 >> 36) + (mod2 & two36m1);
			curr_word >>= (36-shift);
			while(curr_word)
			{
				mod2 += curr_word & two36m1;
				mod2 = (mod2 >> 36) + (mod2 & two36m1);
				curr_word >>= 36;
			}

			nbits += bits[ii];
		}
	  }
	}

	if(nbits != p)
	{
		sprintf(cbuf, "resSH : nbits [%llu] != p [%llu]\n",nbits,p);
		fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}
	ASSERT(HERE, nbits == p    ,"Mlucas.c: nbits == p    ");
	ASSERT(HERE, curr_word == 0,"Mlucas.c: curr_word == 0");

	*Res35m1 = mod1;
	*Res36m1 = mod2;
	return;
}

/********************/

/*
Function to take an n-digit Mersenne or Fermat-mod residue and (if the file pointer
in the argument list is non-null) pretty-print it in hex form to a file.
*/
void 	hex_res_printtofile(double a[], int n, const uint64 p, int timing_test_iters, FILE *fp)
{
	static char *hex_res;
	static
	char c;
	uint64 nbits;
	int bimodn,curr_bits,curr_hexd,cy,findex,ii,j,j1,k,pass,rbits;
	uint64 curr_word;
	int bw,sw,base[2],bits[2];
	int64 itmp;
	int pow2_fft;
	static int first_entry=TRUE;

	ASSERT(HERE, fp != 0x0,"Mlucas.c: fp != 0x0");

	if(first_entry)
	{
		first_entry = FALSE;
		hex_res = (char *)calloc((p+3)/4, sizeof(char));
	}

	ASSERT(HERE,MODULUS_TYPE,"MODULUS_TYPE not set!");
	ASSERT(HERE,MODULUS_TYPE <= MODULUS_TYPE_MAX,"MODULUS_TYPE out of range!");

	ASSERT(HERE,TRANSFORM_TYPE,"TRANSFORM_TYPE not set!");
	ASSERT(HERE,TRANSFORM_TYPE <= TRANSFORM_TYPE_MAX,"TRANSFORM_TYPE out of range!");

	/* If Fermat number, make sure exponent a power of 2: */
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, TRANSFORM_TYPE == RIGHT_ANGLE,"Mlucas.c: TRANSFORM_TYPE == RIGHT_ANGLE");
		findex = trailz64(p);
		ASSERT(HERE, (p >> findex) == 1,"Mlucas.c: (p >> findex) == 1");
		fprintf(fp, "%u iterations of F%u with FFT length %u\n",timing_test_iters,findex,n);
	}
	else
	{
		ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER,"Mlucas.c: TRANSFORM_TYPE == REAL_WRAPPER");
		fprintf(fp, "%u iterations of M%u with FFT length %u\n",timing_test_iters,(uint32)p,n);
	}

	/* Vector length a power of 2? */
	if((n >> trailz32(n)) == 1)
	{
		pow2_fft = TRUE;
	}
	else
	{
		pow2_fft = FALSE;
	}

	bits[0] = p/n;
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
	/*
	If most-significant digit in the balanced-representation form is < 0, add the modulus to the residue.
	For Mersenne (2^p-1) and Fermat (2^p+1) moduli, can combine this with the normalize-to-nonnegative-digit
	step (which we do in any event) by simply initializing the carry into the latter to -1 or +1, respectively:
	*/
	cy=0;		/* init carry.	*/

	for(j=n-1; j >= 0; j -= TRANSFORM_TYPE)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

		if(a[j1]!= 0.0)
		{
			if(a[j1]< 0.0)
			{
				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					cy = -1;
				}
				else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
				{
					cy = +1;
				}
			}
			break;
		}
	}

/*...Now form the hexadecimal residue, converting digits to positive-digit form along the way...	*/

	curr_hexd = 0;	/* Current digit to be written to in the hex_res array */
	nbits = 0;		/* Total bits accumulated so far in the hex residue	*/
	rbits = 0;		/* # of Upper bits left over from processing of previous word	*/
	curr_word = 0;	/*      Upper bits left over from processing of previous word	*/


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
		#ifdef USE_AVX
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			itmp = (int64)(a[j1]+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0)			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
			{
				itmp += (base[ii]);
				cy = -1;
			}
			else
			{
				cy = 0;
			}

			ASSERT(HERE, itmp >= 0,"Mlucas.c: itmp >= 0");
			ASSERT(HERE, rbits < 4,"Mlucas.c: rbits < 4");
			ASSERT(HERE, curr_word < (1<<rbits),"Mlucas.c: curr_word >= 2^rbits!");

			itmp = (itmp << rbits) + curr_word;
			curr_bits = bits[ii] + rbits;

			rbits = curr_bits;
			for(k = 0; k < curr_bits/4; k++)
			{
				hex_res[curr_hexd++] = hex_chars[itmp & 15];
				itmp >>= 4;
				rbits -= 4;
			}
			curr_word = (int)itmp;

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
		#ifdef USE_AVX
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

			itmp = (int64)(a[j1]+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0)			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
			{
				itmp += (base[ii]);
				cy = -1;
			}
			else
			{
				cy = 0;
			}

			ASSERT(HERE, itmp >= 0,"Mlucas.c: itmp >= 0");
			ASSERT(HERE, rbits < 4,"Mlucas.c: rbits < 4");
			ASSERT(HERE, curr_word < (1<<rbits),"Mlucas.c: curr_word >= 2^rbits!");

			itmp = (itmp << rbits) + curr_word;
			nbits += bits[ii];
			curr_bits = bits[ii] + rbits;

			rbits = curr_bits;
			for(k = 0; k < curr_bits/4; k++)
			{
				hex_res[curr_hexd++] = hex_chars[itmp & 15];
				itmp >>= 4;
				rbits -= 4;
			}
			curr_word = (int)itmp;
		}
	  }
	}

	/* Residue should contain ceiling(p/8) bytes: */
	ASSERT(HERE, rbits < 4, "rbits >= 4");
	if(rbits)
	{
		ASSERT(HERE, curr_word < (1<<rbits),"hex_res_printtofile: curr_word >= 2^rbits!");
		hex_res[curr_hexd++] = hex_chars[curr_word & 15];
		curr_word >>= 4;
	}
	ASSERT(HERE, curr_hexd == (p+3)/4,"Mlucas.c: curr_hexd == (p+3)/4");
	ASSERT(HERE, curr_word == 0      ,"Mlucas.c: curr_word == 0      ");
	/* Reverse the string prior to writing to the file: */
	for(ii = 0; ii < curr_hexd/2; ii++)
	{
		c = hex_res[curr_hexd-1-ii]; hex_res[curr_hexd-1-ii] = hex_res[ii]; hex_res[ii] = c;
	}
	fwrite(hex_res, sizeof(char), curr_hexd, fp);
	fprintf(fp, "\n");
	fclose(fp);	fp = 0x0;
	free((void *)hex_res);
	return;
}


/******************Thanks to Tom Cage for the initial version of this: ************************/

#ifdef macintosh
	#include <console.h>	/* Macintosh */
#endif

int 	main(int argc, char *argv[])
{
	int		retVal=0;
	uint64	Res64, Res35m1, Res36m1;
	char	stFlag[STR_MAX_LEN];
	uint32	iters = 0, k = 0, maxFFT, expo = 0, findex = 0;
	double	darg;
	int		new_cfg = FALSE;
	int		i, iarg = 0, idum, lo, hi, start = -1, finish = -1, nargs, scrnFlag, modType = 0, testType = 0, selfTest = 0, userSetExponent = 0, xNum = 0;
	// Force errCheck-always-on in || mode, since we only multithreaded the err-checking versions of the carry routines:
#ifdef MULTITHREAD
	int errCheck = 1;
#else
	int errCheck = 0;
#endif
	int		quick_self_test = 0;
	int		radset = -1;
	double	runtime, runtime_best, tdiff;
	double	roerr_min, roerr_max;
	int		radix_set, radix_best;

	/* Number of distinct FFT lengths supported for self-tests: */
	#define numTest				120	// = sum of all the subranges below
	/* Number of FFT lengths in the various subranges of the full self-test suite: */
	#define numTiny 			32
	#define numSmall			22
	#define numMedium			18
	#define numLarge			 9
	#define numHuge				 9
	#define numEgregious		15
	#define numBrobdingnagian	15
	#define numGodzillian		 0	/* Adding larger FFT lengths to test vectors requires supporting changes to Mdata.h:MAX_FFT_LENGTH_IN_K and get_fft_radices.c */

	#if(numTiny+numSmall+numMedium+numLarge+numHuge+numEgregious+numBrobdingnagian+numGodzillian != numTest)
		#error Sum(numTiny+...) != numTest in main!
	#endif

	struct res_triplet{
		uint64 sh0;	/* Res64 (assuming initial LL test seed = 4) for this exponent. */
		uint64 sh1;	/* Residue mod 2^35-1 */
		uint64 sh2;	/* Residue mod 2^36-1 */
	};

	struct testMers{
		int fftLength;		/* FFT length in K (i.e. 4 means an array of 4K doubles) */
		uint32 exponent;	/* Test exponent for this FFT length. */
		struct res_triplet	res_t[2];	/* 100 and 1000-iteration residue triplet */
	};
	uint32 mvec_res_t_idx = 0;	/* Lookup index into the res_triplet table */
	uint32 new_data;
	struct res_triplet new_res = {0ull,0ull,0ull};

	/* Array of distinct test cases for Mersenne self-tests. Add one extra slot to vector for user-specified self-test exponents;
	We use p's given by given_N_get_maxP(), so maximum RO errors should be consistent and around the target 0.25 value, a little
	bit higher in SSE2 mode.

	The first letter of the descriptoor for each size range serves as a menmonic for the -[*] option which runs the self-tests
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
	struct testMers MersVec[numTest+1] =
	{
	/*                                         100-iteration residues:	                               1000-iteration residues:                */
	/*	  FFTlen(K)     p              Res64           mod 2^35-1      mod 2^36-1               Res64           mod 2^35-1      mod 2^36-1     */
	/*	    -----    --------     ----------------     -----------     -----------         ----------------     -----------     -----------    */
		/* Tiny:                                     [%34359738367  ][%68719476735  ]                         [%34359738367  ][%68719476735  ] */
/*
./Mlucas_new -fftlen  8 -radset 0 -iters 1000
./Mlucas_new -fftlen  9 -radset 0 -iters 1000
./Mlucas_new -fftlen 10 -radset 0 -iters 1000*
./Mlucas_new -fftlen 11 -radset 0 -iters 1000
./Mlucas_new -fftlen 12 -radset 0 -iters 1000
./Mlucas_new -fftlen 13 -radset 0 -iters 1000
./Mlucas_new -fftlen 14 -radset 0 -iters 1000
./Mlucas_new -fftlen 15 -radset 0 -iters 1000
./Mlucas_new -fftlen 16 -radset 0 -iters 1000
./Mlucas_new -fftlen 18 -radset 0 -iters 1000
./Mlucas_new -fftlen 20 -radset 0 -iters 1000
./Mlucas_new -fftlen 22 -radset 0 -iters 1000
./Mlucas_new -fftlen 24 -radset 0 -iters 1000
./Mlucas_new -fftlen 26 -radset 0 -iters 1000
./Mlucas_new -fftlen 28 -radset 0 -iters 1000
./Mlucas_new -fftlen 30 -radset 2 -iters 1000
./Mlucas_new -fftlen 32 -radset 0 -iters 1000
./Mlucas_new -fftlen 36 -radset 0 -iters 1000
./Mlucas_new -fftlen 40 -radset 0 -iters 1000
./Mlucas_new -fftlen 44 -radset 0 -iters 1000
./Mlucas_new -fftlen 48 -radset 0 -iters 1000
./Mlucas_new -fftlen 52 -radset 0 -iters 1000
./Mlucas_new -fftlen 56 -radset 0 -iters 1000
./Mlucas_new -fftlen 60 -radset 3 -iters 1000
./Mlucas_new -fftlen 64 -radset 0 -iters 1000
./Mlucas_new -fftlen 72 -radset 0 -iters 1000
./Mlucas_new -fftlen 80 -radset 0 -iters 1000
./Mlucas_new -fftlen 88 -radset 0 -iters 1000
./Mlucas_new -fftlen 96 -radset 0 -iters 1000
./Mlucas_new -fftlen 104 -radset 0 -iters 1000
./Mlucas_new -fftlen 112 -radset 0 -iters 1000
./Mlucas_new -fftlen 120 -radset 1 -iters 1000

Radix-15 routines need debug:

M2455003: using FFT length 120K = 122880 8-byte floats.
 this gives an average   19.978865559895834 bits per digit
Using complex FFT radices        15        16        16        16
Output a[8194] =        -568516.00000 out of range: base[1] =         1048576
ERROR: at line 1990 of file mers_mod_square.c
Assertion failed: Output out of range!
*/
		{     8,    173431u, { {0x85301536E4CA9B11ull,  3707224323ull, 36851834664ull}, {0x2FD5120BEC41F449ull, 28734955954ull, 23103392631ull} } },
		{     9,    194609u, { {0xC711AF1008612BC6ull,  1574019740ull, 37260026270ull}, {0x5153F6E040CD1BE6ull, 15446410924ull,  3291404673ull} } },
		{    10,    215767u, { {0x4428783BC62760F0ull,  7466284975ull, 53916123655ull}, {0xED46A8C001908815ull,   739143119ull, 36950829937ull} } },
		{    11,    236813u, { {0x592D849AF4D1336Full, 29025996994ull, 48905971124ull}, {0xB4EEB63BB656F424ull,  5361680901ull, 31850818767ull} } },
		{    12,    257903u, { {0x1D8121DE28B60996ull, 22402402961ull, 65583959369ull}, {0x54F2BE961A674CB1ull, 25601315553ull, 54298899520ull} } },
		{    13,    278917u, { {0xE3BC90B0E652C7C0ull, 21244206101ull, 51449948145ull}, {0x93AF8994F95F2E50ull, 16427368469ull, 10707190710ull} } },
		{    14,    299903u, { {0xDB8E39C67F8CCA0Aull, 20506717562ull, 44874927985ull}, {0x4E7CCB446371C470ull, 34135369163ull, 61575700812ull} } },
		{    15,    320851u, { {0xB3C5A1C03E26BB17ull, 22101045153ull,  4420560161ull}, {0x923A9870D65BC73Dull, 29411268414ull, 30739991617ull} } },
		{    16,    341749u, { {0x8223DF939E46A0FFull, 32377771756ull, 38218252095ull}, {0xC6A5D4B6034A34B8ull, 31917858141ull, 59888258577ull} } },
		{    18,    383521u, { {0xBF30D4AF5ADF87C8ull, 15059093425ull, 52618040649ull}, {0x9F453732B3FE3C04ull,  4385160151ull, 47987324636ull} } },
		{    20,    425149u, { {0x6951388C3B99EEC0ull,  4401287495ull, 19242775142ull}, {0x501CEC2CB2080627ull, 21816565170ull, 41043945930ull} } },
		{    22,    466733u, { {0xD95F8EC0F32B4756ull, 19305723506ull, 26588871256ull}, {0xB1F58184918D94B6ull,  8443388060ull, 11738516313ull} } },
		{    24,    508223u, { {0xDA46E41316F8BCCAull, 25471180026ull,  1635203275ull}, {0x27A5B285281466B9ull, 11438869313ull,  7226774009ull} } },
		{    26,    549623u, { {0x6649D9D6CD4E0CE1ull, 25445908581ull, 26118212198ull}, {0x1A4F280627A15B3Cull, 13286323782ull, 31550278005ull} } },
		{    28,    590963u, { {0x4ADDB6C4A76465AFull,  6532108269ull, 54921134131ull}, {0x3063D08A7BABD7B8ull,  4777711548ull, 39733274344ull} } },
		{    30,    632251u, { {0x0811FAA40601EB1Dull, 16369365746ull,  6888026123ull}, {0xF324E4DEC564AF91ull, 10236920023ull, 34068699974ull} } },
		{    32,    673469u, { {0x1A4EF8A0D172FBAAull, 32667536946ull, 11393278588ull}, {0xA4DFD62B928F68A4ull, 11900420802ull, 66610946021ull} } },
		{    36,    755737u, { {0x13B13C61298088DCull, 34092803628ull,  7584858890ull}, {0x33A2A43DE8782CCCull,  2953124985ull, 62434716987ull} } },
		{    40,    837817u, { {0x88555D9AAD3FF2DDull,  8573348747ull, 67896670216ull}, {0xEAC1676D914878C0ull, 34312095136ull, 45077378164ull} } },
		{    44,    919729u, { {0x6ACC03213A37BA5Bull,  3870201113ull, 48145739792ull}, {0xDA98B49CC83C60CBull, 15886769401ull, 62221100895ull} } },
		{    48,   1001467u, { {0x6B1C76AB5431FDA4ull,  6795121241ull, 65308927583ull}, {0xBD99FD21F4136BFCull, 26386988063ull, 61607603549ull} } },
		{    52,   1083077u, { {0xA591637EC8CF3FE4ull,  4769775755ull, 65114872367ull}, {0xE59C08B13B00E6FFull,  1383963096ull, 26100699764ull} } },
		{    56,   1164533u, { {0xEC4F2579E4533584ull,  5456769127ull, 59922459736ull}, {0xF7D2BF94C2767D36ull, 30727892629ull, 48141128220ull} } },
		{    60,   1245877u, { {0xC91002E1A4EE7E07ull,  6217476228ull, 40164514288ull}, {0xEABE9E1A31DF5877ull,   831216169ull, 29591771932ull} } },
		{    64,   1327099u, { {0xAC070112281229E0ull, 14226353454ull,  1640524016ull}, {0xF25AA54053C5BB64ull, 32455038659ull, 53547160776ull} } },
		{    72,   1489223u, { {0x6674518EA19B3D6Aull, 32383400562ull, 53234746310ull}, {0xEB312091097F6C3Bull,  3980687908ull,  8568698675ull} } },
		{    80,   1650959u, { {0xE5326E754F3202A8ull,  5593233426ull, 33337128557ull}, {0xFC3E8CDA60AF5CF8ull, 11466296968ull, 12651602524ull} } },
		{    88,   1812347u, { {0x81BDD3AC63DF3F73ull, 19957199947ull, 61002681293ull}, {0x3D3E429D7427C4EAull, 25342898119ull, 34322985438ull} } },
		{    96,   1973431u, { {0x901C8305DA9FF95Aull, 32611811878ull, 55986702160ull}, {0x0790CA11ADAA47E3ull, 17075140845ull, 12883521448ull} } },
		{   104,   2134201u, { {0x59BDA0D80F3279EDull, 17901153436ull,  3927067335ull}, {0x2F81B21BC680C861ull, 18443771511ull, 45465079919ull} } },
		{   112,   2294731u, { {0xC44ACC96D268625Full, 10331638988ull,  2292055445ull}, {0xED20577E16E128DEull, 32248607028ull, 14903460370ull} } },
		{   120,   2455003u, { {0xC5F7DB23F174A67Dull, 32991574397ull, 31642856976ull}, {0x401670254012E5ABull, 33626385418ull, 66465546971ull} } },
		/* Small: */
		{   128,   2614999u, { {0x040918890E98F8DAull, 14867710211ull, 47602627318ull}, {0x1A184504D2DE2D3Cull,  5934292942ull,  4090378120ull} } },
		{   144,   2934479u, { {0x1B90A27301980A3Aull,  7043479338ull, 38327130996ull}, {0x8C3045C6534867C6ull, 12456621644ull, 52801948293ull} } },
		{   160,   3253153u, { {0x9AFD3618C164D1B4ull, 16551334620ull, 55616214582ull}, {0x1493A70897A8D058ull, 34082962858ull, 60773088284ull} } },
		{   176,   3571153u, { {0xA016F25779902477ull, 21500047857ull,  9150810891ull}, {0x8E6F248EC96445FFull, 22443629034ull, 16625023722ull} } },
		{   192,   3888509u, { {0x71E61322CCFB396Cull, 29259839105ull, 50741070790ull}, {0x3CEDB241702D2907ull,  6177258458ull, 21951191321ull} } },
		{   208,   4205303u, { {0xC08562DA75132764ull,  7099101614ull, 36784779697ull}, {0xAD381B4FE91D46FDull,  7173420823ull, 51721175527ull} } },
		{   240,   4837331u, { {0xB0D0E72B7C87C174ull, 15439682274ull, 46315054895ull}, {0x3AA14E0E90D16317ull,  5730133308ull, 50944347816ull} } },
		{   224,   4521557u, { {0xE68210464F96D6A6ull, 20442129364ull, 11338970081ull}, {0x3B06B74F5D4C0E35ull,  7526060994ull, 28782225212ull} } },
		{   256,   5152643u, { {0x074879D86679CB5Bull,  1208548964ull, 48525653083ull}, {0x98AF5E14C824A252ull,   783196824ull,  6594302302ull} } },
		{   288,   5782013u, { {0x9869BE81D9AB1564ull, 15509103769ull, 49640026911ull}, {0x7C998719C6001318ull, 23749848147ull, 19853218689ull} } },
		{   320,   6409849u, { {0x20739E43A693A937ull, 27970131351ull, 15334307151ull}, {0xE20A76DCEB6774A6ull, 14260757089ull, 68560882840ull} } },
		{   352,   7036339u, { {0xD6A226BAB5E14D62ull,  9444972171ull, 28639488018ull}, {0x0579D28296F29D92ull, 18964853245ull, 30111685201ull} } },
		{   384,   7661567u, { {0x3A929F577AC9725Full, 23890492835ull, 64821764767ull}, {0x2ECBA785576E6D58ull, 26446200615ull, 60269798452ull} } },
		{   416,   8285659u, { {0xDCA138D55C36E40Cull, 10452583294ull,  4308453248ull}, {0x1FEE7F79E32229A6ull, 11936075329ull, 16061515794ull} } },
		{   448,   8908723u, { {0x436494C7EA194FA1ull,  2976149044ull, 21645125251ull}, {0xD2FCEDE29E818A26ull, 15260150529ull, 11678366985ull} } },
		{   480,   9530803u, { {0x6D52BCCB796A46E9ull, 28296985362ull, 66294636306ull}, {0x786CB610A809B762ull, 24654197494ull, 57943258783ull} } },
		{   512,  10151971u, { {0x7BA3DD9C38878B83ull, 24395728676ull, 12632037498ull}, {0x8E9FA3093ACD81C1ull,  6345070464ull, 65203567231ull} } },
		{   576,  11391823u, { {0x78D8769B9F75FB2Bull,  7286184017ull, 17306427758ull}, {0x5834BDA7558DD43Cull, 11189092321ull, 23026923116ull} } },
		{   640,  12628613u, { {0xF951B697F5C5032Aull,  9262385660ull, 57463106915ull}, {0x93B526040205BA27ull, 30538079080ull, 32317022014ull} } },
		{   704,  13862759u, { {0xBB2F69275D79A9EEull, 12230183669ull, 68684647134ull}, {0x7343ECC160AA00D5ull, 24655585170ull, 51102704879ull} } },
		{   768,  15094403u, { {0xF6895EB66EADE9C5ull,  8490184692ull, 23393343807ull}, {0xF673A8D6413923A9ull, 20026779905ull, 67516766223ull} } },
		{   832,  16323773u, { {0xEB8890F379392B2Full, 27289972116ull, 63975275393ull}, {0xD681EDD3A1EC3780ull, 12515962698ull, 40155157152ull} } },
		{   896,  17551099u, { {0xAB1180428ED65EE0ull,  3105108668ull, 66518734167ull}, {0x31813367849BBF49ull,  9516734777ull, 18271834608ull} } },
		{   960,  18776473u, { {0xCA7D81B22AE24935ull, 24317941565ull, 67706175547ull}, {0x02EB980A49E7B60Full,  5730644436ull, 48386545950ull} } },
		/* Medium: */
		{  1024,  20000047u, { {0xDD61B3E031F1E0BAull,  6238131189ull, 41735145962ull}, {0x6F0542D58BE1D854ull, 33550477926ull,  4560911335ull} } },
		{  1152,  22442237u, { {0x62C479B03F3E9DD9ull, 16186007845ull, 66602070649ull}, {0xF8214922C0726AA9ull,  2103619488ull, 63860960236ull} } },
		{  1280,  24878401u, { {0x8A8644FC94CB0A8Bull, 31904286697ull, 27818942965ull}, {0xAB247C597DCD678Aull, 24768590666ull, 25504625244ull} } },
		{  1408,  27309229u, { {0xCCE2DF04E61DC922ull, 20296295440ull, 23445145000ull}, {0xBC6E94B454AF005Dull, 20163791999ull, 57504620236ull} } },
		{  1536,  29735137u, { {0x2D26046FFAAEBC2Bull,  2381591353ull, 48693163035ull}, {0x00E993C62E8A75D6ull, 18680287655ull, 52664049258ull} } },
		{  1664,  32156581u, { {0x77E274E6C29C203Eull,  8887595531ull,  5248686295ull}, {0xAF77EA82895CC6A8ull,  5759391333ull, 18391720566ull} } },
		{  1792,  34573867u, { {0x87F3FA16F713EF7Full,   679541321ull, 62692450676ull}, {0x10EA1A34F05813C6ull,  9097156999ull, 67076927840ull} } },
		{  1920,  36987271u, { {0x5EFE558EF126B1C8ull, 22078146080ull,   989319562ull}, {0x4865E07EB07F7CE5ull, 10419976906ull,  2724830069ull} } },
		{  2048,  39397201u, { {0x6179CD26EC3B3274ull,  8060072069ull, 29249383388ull}, {0x81AEAC0C7E6089BBull, 25132671466ull, 41950605021ull} } },
		{  2304,  44207087u, { {0x6EE085510D8C8F39ull, 24377809749ull, 55912147302ull}, {0x821F029CE53B2FDAull, 30673953109ull, 44311760422ull} } },
		{  2560,  49005071u, { {0x07EFE3EF1F78E763ull, 22407816581ull, 54111649274ull}, {0xB8D8C7CFFCF4125Cull,  8330528102ull,  6182508157ull} } },
		{  2816,  53792327u, { {0x6C658F13E0F4A102ull,  1998570299ull, 59196879534ull}, {0x95A2D7F362429603ull, 33703733096ull, 60474218222ull} } },
		{  3072,  58569809u, { {0x3CBFF8D5BC5CDDB6ull,  1737566556ull, 47007808842ull}, {0xDC08E05E46E39ECFull,  5122766298ull, 63677326652ull} } },
		{  3328,  63338459u, { {0xF6D5A31C4CE769D0ull, 32776678711ull, 12897379614ull}, {0x59BA4BBB9F240C23ull, 17218869802ull, 56649053764ull} } },
		{  3584,  68098843u, { {0x46675D1C9B599EFEull,  8451437941ull,  2732573183ull}, {0x0AF9DBE03FE3939Full, 28618044864ull, 16348034528ull} } },
		{  3840,  72851621u, { {0x5A313D90DECAEF25ull, 20007549674ull, 66588163113ull}, {0xE54B85058DB80689ull,  6488756762ull, 45980449388ull} } },
		/* Large: */
		{  4096,  77597293u, { {0x8CC30E314BF3E556ull, 22305398329ull, 64001568053ull}, {0x5F87421FA9DD8F1Full, 26302807323ull, 54919604018ull} } },
		{  4608,  87068977u, { {0xE6FDFBFC6600B9D8ull, 15739469668ull, 29270706514ull}, {0x1CE06D2ADF8238DCull, 32487573228ull, 43042840938ull} } },
		{  5120,  96517019u, { {0x49D44A10AC9EA1D7ull,  1040777858ull, 40875102228ull}, {0xC528F96FE4E06C17ull, 12554196331ull, 58274243273ull} } },
		{  5632, 105943723u, { {0x031835135A45506Cull, 25880831521ull, 24662569542ull}, {0x52811E5306BD680Aull, 11950909209ull, 30005399045ull} } },
		{  6144, 115351063u, { {0x32246D04038A48A5ull,  3869276047ull, 58028049428ull}, {0xC2568F0B9908ACBEull, 20470377290ull, 53793586344ull} } },
		{  6656, 124740697u, { {0x0362F83AAE178D1Aull, 11483908585ull, 30666543873ull}, {0xAD454FE846C84018ull, 17578866610ull, 41948771933ull} } },
		{  7168, 134113933u, { {0x1C97C7AAEC5CB8C1ull,  2508417425ull, 49620821143ull}, {0xED969839710C0872ull, 27772820714ull, 45853381796ull} } },
		{  7680, 143472073u, { {0x8FBCF315FA8C8BDAull, 15967283992ull, 62343341476ull}, {0xA0A5A19324FB17DFull,  4384978845ull, 65655725783ull} } },
		{  8192, 152816047u, { {0xB58E6FA510DC5049ull, 27285140530ull, 16378703918ull}, {0x75F2841AEBE29216ull, 13527336804ull,   503424366ull} } },
		/* Huge: */
		{  9216, 171465013u, { {0x60FE24EF89D6140Eull, 25324379967ull,  3841674711ull}, {0x6753411471AD8945ull, 17806860702ull,  3977771754ull} } },
		{ 10240, 190066777u, { {0x65CF47927C02AC8Eull, 33635344843ull, 67530958158ull}, {0xBADA7FD24D959D21ull, 12777066809ull, 67273129313ull} } },
		{ 11264, 208626181u, { {0x6FC0151B81E5173Full, 29164620640ull, 19126254587ull}, {0xD74AA66757A5345Eull, 17524190590ull, 14029371481ull} } },
		{ 12288, 227147083u, { {0xE01AE9C859ADB03Aull,  7273133358ull,   681418986ull}, {0x303F142E1E88D5B4ull, 28479237457ull, 42044197589ull} } },
		{ 13312, 245632679u, { {0x0A6ACB405ADC0354ull,    39452330ull, 38999048555ull}, {0xB38B02A4F195762Full,  3280152282ull, 30314100936ull} } },
		{ 14336, 264085733u, { {0x5ACE4CCE3B925A81ull,  4584210608ull, 36618317213ull}, {0x02F5EC0CBB1C2032ull, 27165893636ull,   687123146ull} } },
		{ 15360, 282508657u, { {0xE7B08ED3A92EC6ECull,   875689313ull, 41754616020ull}, {0xD08FBAFF5CA5096Full, 30398073011ull, 62088094181ull} } },
		{ 16384, 300903377u, { {0xA23E8D2F532F05E6ull, 17871262795ull, 53388776441ull}, {0x14F20059083BF452ull, 16549596802ull, 56184170215ull} } },
		{ 18432, 337615277u, { {0xAEB976D153A4176Bull, 15040345558ull, 14542578090ull}, {0x503B443CB1E0CD2Dull, 29149628739ull,  5785599363ull} } },
		/* Egregious: */
		{ 20480, 374233309u, { {0x6D95C0E62C8F9606ull,  3426866174ull, 39787406588ull}, {0xD08FB9031D460B7Eull, 30083048700ull, 30636357797ull} } },
		{ 22528, 410766953u, { {0x254572CAB2014E6Cull, 17794160487ull, 13396317974ull}, {0x6202B11AA8602531ull,  9896689730ull, 13179771096ull} } },
		{ 24576, 447223969u, { {0x547DF462C7DAD1F6ull, 21523243373ull, 60663902864ull}, {0x06AB4D8E6FD9D69Bull, 23751736171ull, 10720018557ull} } },
		{ 26624, 483610763u, { {0xED3E248A29A1C6A8ull,  3863520901ull, 56560420765ull}, {0xC29358F8206746D6ull, 28829535828ull,  8160695393ull} } },
		{ 28672, 519932827u, { {0xCA7B3A76819D67F7ull,  7078504016ull, 32836389262ull}, {0x5799C8BE8E02B56Full, 22269194969ull, 11617462155ull} } },
		{ 30720, 556194803u, { {0x56CDF80EFF67C1C1ull, 15611746359ull, 45837150207ull}, {0xCBC88456B4B47AC0ull,  3535843471ull, 18652930008ull} } },
		{ 32768, 592400713u, { {0xD8C829884B234EB4ull, 13278383896ull, 54835531046ull}, {0x994DD6B24F452451ull, 28289805997ull, 11462134131ull} } },
		{ 36864, 664658101u, { {0x2A809B0C735BAC4Bull, 10513414423ull, 54347266174ull}, {0xAB2147D9BAA22BB4ull, 12259954326ull, 67125404781ull} } },
		{ 40960, 736728527u, { {0xB9AC3EC848FF60A5ull,  7352037613ull,  7261166574ull}, {0x3D623A79D0F14EFFull, 31246254654ull, 49195074754ull} } },
		{ 45056, 808631029u, { {0x9D543D67F48AF766ull,  4288733345ull, 27338399854ull}, {0x62A4DF80612E897Bull, 32232536296ull, 47296762118ull} } },
		{ 49152, 880380937u, { {0x82ED59E22D585BF6ull, 34028441473ull, 54282489944ull}, {0x7F110FD687DB7CB5ull, 14440703450ull, 57885720403ull} } },
		{ 53248, 951990961u, { {0xE98B9D69E5F350D1ull, 14692034938ull,  2987127112ull}, {0x634157BCC596287Dull, 31799414844ull, 64590776355ull} } },
		{ 57344,1023472049u, { {0x960CE3D7029BDB70ull, 30108539760ull,  2075307892ull}, {0x9D98FF9A3FD21D1Bull, 10507186734ull,  3891581073ull} } },
		{ 61440,1094833457u, { {0x9A2592E96BC1C827ull, 40462567463ull, 18678252759ull}, {0x79BF2E2F5AE7985Bull,  8407199527ull, 64114889744ull} } },
		{ 65536,1166083297u, { {0xA0FE3066C834E360ull, 29128713056ull,  7270151463ull}, {0x4CC11F1C55F95C9Bull, 17910403663ull, 19262812656ull} } },
		/* Brobdingnagian: */
		{ 73728,1308275261u, { {0xE97D94366B318338ull, 27568210744ull, 15277294114ull}, {0x0913DB564B5FB2C9ull, 22759704201ull,  4713378634ull} } },
		{ 81920,1450094993u, { {0x0F0945ACF6C46D21ull, 55679675681ull,  4898080803ull}, {0x12FB0E6628819A5Aull, 26252443199ull, 65921528136ull} } },
		{ 90112,1591580099u, { {0xD8B487D8F5AFDF1Aull, 38481682202ull, 32203366427ull}, {0xAA5F7E4D5EC75AD3ull, 16596905178ull,  8253283767ull} } },
		{ 98304,1732761197u, { {0xB3ABD7E277F6EFB4ull, 10602606516ull, 18939646037ull}, {0x2170E14E6B46F009ull, 25931343295ull, 29383122206ull} } },
		{106496,1873663819u, { {0x938367CFDA83235Cull, 68090536796ull,   536356188ull}, {0xE4F8B81BD2E42112ull, 25813505565ull, 59441464991ull} } },
		{114688,2014309639u, { {0x84C8DED05749A21Full,  1464443423ull,  2648601266ull}, {0xDE2C03AA5068CC6Aull, 17222310589ull, 55630799101ull} } },
		{122880,2154717007u, { {0xE3056185380C64CFull, 22415172815ull, 18918434052ull}, {0x8DE279575CD1154Eull, 32123477176ull, 52211425638ull} } },
		{131072,2294901977u, { {0xB8CC9297E6E7512Cull, 33938690348ull, 32834437586ull}, {0x77FB54BF318C2DBAull, 20797246137ull, 29632517743ull} } },
		{147456,2574659081u, { {0x715208AF863EBA5Full, 66676767327ull, 17913296158ull}, {0xAA2E017A840AE336ull, 17798910323ull, 59557563397ull} } },
		{163840,2853674573u, { {0x2DF960079E128064ull, 32716783716ull,   830562867ull}, {0x341EB1C8403CE6AEull, 11934818175ull, 64851277709ull} } },
		{180224,3132023311u, { {0xD320ADB4B1483D7Eull, 20154170750ull,  9281699703ull}, {0x5008B456D49F551Full,  9130349955ull, 45637727974ull} } },
		{196608,3409766351u, { {0x1760B9AA1C6C9DFDull, 43426553341ull, 28002093744ull}, {0x0000000000000000ull,           0ull,           0ull} } },
		{212992,3686954519u, { {0x9A3231FF531C8617ull, 10531609552ull, 42475524136ull}, {0x0CD00ADEE0D4D0DBull,  1995928646ull, 59987670753ull} } },
		{229376,3963630893u, { {0x427266E0A64E77F8ull,   135586632ull,  3161584476ull}, {0xBA1F2CE3D59018E7ull,  3171773156ull, 31951413199ull} } },
		{245760,4239832153u, { {0x05FF04CA8AEA8E54ull,   718789070ull, 67021029350ull}, {0xDCA27889AB5F4D92ull,    55229805ull, 33882881510ull} } },
		/* Godzillian requires 64-bit exponent support: */
		{     0,         0u, { {                 0ull,           0ull,           0ull}, {                 0ull,           0ull,           0ull} } }
	};
/*
Using complex FFT radices        18         8         8         8        16        16        16
1000 iterations of M1308275261 with FFT length 75497472 = 73728 K
Res64: . AvgMaxErr = 0.226149827. MaxErr = 0.298042297. Program: E3.0x
Clocks = 03:13:02.209

Using complex FFT radices        11         8         8        16        16        16        16
1000 iterations of M1591580099 with FFT length 92274688 = 90112 K
Res64: . AvgMaxErr = 0.245961766. MaxErr = 0.309280396. Program: E3.0x
Clocks = 04:33:31.880

Using complex FFT radices        13         8         8        16        16        16        16
1000 iterations of M1873663819 with FFT length 109051904 = 106496 K
Res64: . AvgMaxErr = 0.265383741. MaxErr = 0.322921753. Program: E3.0x
Clocks = 05:10:53.619

Using complex FFT radices        14         8         8        16        16        16        16
1000 iterations of M2014309639 with FFT length 117440512 = 114688 K
Res64: . AvgMaxErr = 0.235015343. MaxErr = 0.285949707. Program: E3.0x
Clocks = 05:42:18.229	<*** 4/22/09: SSE2-enabled radix-28+32,16,16,16,16 needs 4h50m***

Using complex FFT radices        15         8         8        16        16        16        16
1000 iterations of M2154717007 with FFT length 125829120 = 122880 K
Res64: . AvgMaxErr = 0.244326867. MaxErr = 0.320281982. Program: E3.0x
Clocks = 05:57:22.259

Using complex FFT radices        16         8         8        16        16        16        16
1000 iterations of M2294901977 with FFT length 134217728 = 131072 K
Res64: . AvgMaxErr = 0.252516914. MaxErr = 0.312500000. Program: E3.0x
Clocks = 05:39:16.959	<*** 4/22/09: SSE2-enabled radix-32+32,16,16,16,16 needs 5h34m***

Using complex FFT radices        18         8         8        16        16        16        16
1000 iterations of M2574659081 with FFT length 150994944 = 147456 K
Res64: . AvgMaxErr = 0.226303838. MaxErr = 0.272964478. Program: E3.0x
Clocks = 07:06:10.860

Using complex FFT radices        20        32        16        16        32        16
1000 iterations of M2853674573 with FFT length 167772160 = 163840 K
Res64: . AvgMaxErr = 0.241956951. MaxErr = 0.296875000. Program: E3.0x
Clocks = 06:48:30.500

Using complex FFT radices        11         8        16        16        16        16        16
1000 iterations of M3132023311 with FFT length 184549376 = 180224 K
Res64: . AvgMaxErr = 0.239856285. MaxErr = 0.295349121. Program: E3.0x
Clocks = 09:14:00.330

M3409766351: using FFT length 196608K = 201326592 8-byte floats.
 this gives an average   16.936492676536243 bits per digit
Using complex FFT radices        24        32        16        16        32        16
FATAL: iter =          2; nonzero exit carry in radix24_ditN_cy_dif1 - input wordsize may be too small.

Using complex FFT radices        13         8        16        16        16        16        16
1000 iterations of M3686954519 with FFT length 218103808 = 212992 K
Res64: . AvgMaxErr = 0.259983692. MaxErr = 0.312774658. Program: E3.0x
Clocks = 10:29:35.589

Using complex FFT radices        14         8        16        16        16        16        16
1000 iterations of M3963630893 with FFT length 234881024 = 229376 K
Res64: . AvgMaxErr = 0.233737114. MaxErr = 0.301208496. Program: E3.0x
Clocks = 11:27:59.089

Using complex FFT radices        15         8        16        16        16        16        16
1000 iterations of M4239832153 with FFT length 251658240 = 245760 K
Res64: . AvgMaxErr = 0.242665772. MaxErr = 0.300384521. Program: E3.0x
Clocks = 12:03:41.449

*/

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
		struct res_triplet	res_t[2];	/* 100 and 1000-iteration residue triplet */
	};

	/* Array of 100-iteration reference values for Fermat self-tests. Only allow Fidx >= 14: */
	#define FermArrayIdxOffset		14	/* Amount to subtract from m to get the prestored residues for F_m */
	#define numFerm		17
	struct testMers FermVec[numFerm+1] =
	{
	/*                                   100-iteration residues:                                  1000-iteration residues:              */
	/* FFTlen(K) Fidx           Res64           mod 2^35-1      mod 2^36-1               Res64           mod 2^35-1      mod 2^36-1     */
	/*   ------- ----      ----------------     -----------     -----------         ----------------     -----------     -----------    */
		/*                                    [%34359738367  ][%68719476735  ]                         [%34359738367  ][%68719476735  ] */
		{     1,  14, { {0xDB9AC520C403CB21ull,   342168579ull, 59244817440ull}, {0xF111F12732CCCB0Full, 24848612524ull, 66609820796ull} } },
		{     2,  15, { {0x3B21A6E55ED13454ull, 28379302213ull, 15546218647ull}, {0x4784657F2A36BE74ull,   617376037ull, 44891093359ull} } },
		{     4,  16, { {0xAAE76C15C2B37465ull, 20013824731ull,  2261076122ull}, {0x42CC2CBE97C728E6ull, 30814966349ull, 44505312792ull} } },
		{     8,  17, { {0xFFA16CDC8C87483Cull, 20917337408ull, 26110327818ull}, {0x43CAB295FFB2661Full, 18197605796ull,  9842643677ull} } },
		{    16,  18, { {0x7C6B681485EB86DBull,  5745147782ull, 50521157289ull}, {0x8193BD41931E9DE8ull, 19662968587ull, 51102742548ull} } },
		{    32,  19, { {0x529E54642A813995ull, 17797950508ull, 32039741221ull}, {0xE24EAE4B153EE86Bull, 11155350666ull, 49866866361ull} } },
		{    64,  20, { {0x64629CED6E218018ull,  8485981669ull, 53977437340ull}, {0xA380121F6FD26B2Aull, 15876203498ull, 36314727556ull} } },
		{   128,  21, { {0x1DE0171591038250ull, 33758422990ull,  8269940507ull}, {0x1B331FBB41AF33D7ull, 17971032338ull,  2929392342ull} } },
		{   256,  22, { {0x9201143390F3828Dull,  1749100092ull, 46602233256ull}, {0x1B331FBB41AF33D7ull, 17971032338ull,  2929392342ull} } },
		{   512,  23, { {0x9C3F8E29B397B32Bull,  1094055486ull, 13316822657ull}, {0xBD642EA0479D8FF0ull, 31625967305ull, 57187857233ull} } },
		{  1024,  24, { {0xDB9F01963ED9DC8Bull, 27887793041ull, 13169874547ull}, {0x40F2DECE9C351236ull,  9074286032ull, 38590877049ull} } },
		{  2048,  25, { {0x376C33921E5F675Full, 13022327996ull, 46818697393ull}, {0xA51F8577A407CB75ull,  9865976783ull, 35171498411ull} } },
		{  4096,  26, { {0xA42BECD80DAEC4CBull, 10087739060ull, 25252768685ull}, {0xECC9408A7295401Dull,  5904751941ull, 58967745948ull} } },
		{  8192,  27, { {0xFB69E377519D8CE6ull, 15449775614ull, 51221672039ull}, {0x24898E3BEB59DCE6ull, 24957168001ull,  2072452827ull} } },
		{ 16384,  28, { {0xA4FF6F8C3CB38B85ull, 18933356966ull, 30899345457ull}, {0x8B451AF25E8CC50Eull,   674652743ull, 39963850167ull} } },
		{ 32768,  29, { {0xAFBF110B593E26F6ull, 32666279868ull, 18995112582ull}, {0xA6B643FF24C6ADC1ull, 15753158767ull, 13965270144ull} } },
		{ 65536,  30, { {0x68B1BDA5D6BAE04Bull,  3347054148ull, 47892955488ull}, {0x0000000000000000ull,           0ull,           0ull} } },
		{     0,   0, { {                 0ull,           0ull,           0ull}, {0x0000000000000000ull,           0ull,           0ull} } }
	};

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


	fprintf(stderr, "\n    Mlucas %s\n", VERSION);
	fprintf(stderr, "\n    http://hogranch.com/mayer/README.html\n\n");

#ifdef macintosh
	argc = ccommand(&argv);			/* Macintosh CW */
#endif

if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
{
	ASSERT(HERE, (MersVec[numTest-1].fftLength != 0) &&  (MersVec[numTest].fftLength == 0), "numTest != MersVec allocated size!");
}
else
{
	ASSERT(HERE, (FermVec[numFerm-1].fftLength != 0) &&  (MersVec[numTest].fftLength == 0), "numTest != FermVec allocated size!");
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

	/* Currently only primality-test mode is supported: */
	testType = TEST_TYPE_PRIMALITY;

	if (argc == 1)	/* Default Mode */
	{
	ERNST_MAIN:
		if((retVal = ernstMain(MODULUS_TYPE_MERSENNE,TEST_TYPE_PRIMALITY,0,0,0,0,0,FALSE,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime)) != 0)
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

	/******** command-line-argument processing while() loop: ********/
	nargs = 1;
	while(argv[nargs])
	{
		strncpy(stFlag, argv[nargs++], STR_MAX_LEN);

		if(stFlag[0] != '-')
		{
			fprintf(stderr, "*** ERROR: Illegal command-line option %s\n", stFlag);
			fprintf(stderr, "*** All command-line options must be of form -{flag} [argument]\n\n");
			goto MLUCAS_HELP;
		}

		if(STREQ(stFlag, "-h"))
		{
			goto MLUCAS_HELP;
		}

		/* Mersenne self-test: requires a user-set exponent, FFT length or one of the supported -s arguments below: */
		if(STREQ(stFlag, "-s"))
		{
			selfTest = TRUE;

		/* -s [arg] cannot be combined with any other flag except -rocheck [on|off]:
			if(userSetExponent)
			{
				sprintf(cbuf  , "*** ERROR: The -s [arg] flag cannot be combined with any other except -rocheck [on|off].\n");
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
		*/

			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				goto MLUCAS_HELP;
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
				if(STREQ(stFlag, "b") || STREQ(stFlag, "gojira"))
				{
					break;
				}

				fprintf(stderr, "*** ERROR: Illegal argument %s to -s flag\n", stFlag);
				goto MLUCAS_HELP;
			}

			modType  = MODULUS_TYPE_MERSENNE;
			/* Don't break yet, since we need to see if -rocheck was invoked by the caller. */
		}

		else if(STREQ(stFlag, "-iters"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				goto MLUCAS_HELP;
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
						sprintf(cbuf  , "*** ERROR: -iters argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					sprintf(cbuf  , "*** ERROR: Non-numeric character encountered in -iters argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}

			iters = (uint32)iarg;
		}

	// Force errCheck-always-on in || mode, since we only multithreaded the err-checking versions of the carry routines:
	#ifndef MULTITHREAD
		else if(STREQ(stFlag, "-rocheck"))
		{
			errCheck = 1;
		}
	#endif

		else if(STREQ(stFlag, "-fftlen"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				goto MLUCAS_HELP;
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
						sprintf(cbuf  , "*** ERROR: -fftlen argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					sprintf(cbuf  , "*** ERROR: Non-numeric character encountered in -fftlen argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}

			k = (uint32)iarg;
			if((i = get_fft_radices(k, 0, 0x0, 0x0, 0)) != 0)
			{
				sprintf(cbuf  , "ERROR: FFT length %d K not available.\n",k);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}

			/* Don't set userSetExponent here, since -fftlen can be invoked without an explicit exponent */
			if(modType == MODULUS_TYPE_FERMAT)
			{
				FermVec[numFerm].fftLength = k;	/* For this to work reqire any explicit mod-type invocation to occur before fftlen in input arglist */
				start = numFerm; finish = start+1;
			}
			else
			{
				MersVec[numTest].fftLength = k;
				start = numTest; finish = start+1;
			}
		}

		else if(STREQ(stFlag, "-radset"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				goto MLUCAS_HELP;
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
						sprintf(cbuf  , "*** ERROR: -radset argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					sprintf(cbuf  , "*** ERROR: Non-numeric character encountered in -radset argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}
			radset = iarg;

			/* Make sure user has specified an FFT length: */
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


		else if(STREQ(stFlag, "-nthread"))
		{
			if(nargs >= argc)
			{
				fprintf(stderr, "*** ERROR: Unterminated command-line option or malformed argument.\n");
				goto MLUCAS_HELP;
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
						sprintf(cbuf  , "*** ERROR: -nthread argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					sprintf(cbuf  , "*** ERROR: Non-numeric character encountered in -nthread argument %s.\n", stFlag);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}

		#ifndef MULTITHREAD
			fprintf(stderr,"Multithreading not enabled; ignoring -nthread argument.");
		#else
			NTHREADS = iarg;
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
				goto MLUCAS_HELP;
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
						sprintf(cbuf  , "*** ERROR: -m argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					sprintf(cbuf  , "*** INFO: Non-numeric character encountered in -m argument %s ... using default exponent for whatever FFT length is supplied.\n", stFlag);
					fprintf(stderr,"%s", cbuf);
					--nargs;
					goto SET_MERS;
				}
			}

			expo = (uint32)iarg;
			userSetExponent = 1;
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
				goto MLUCAS_HELP;
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
						sprintf(cbuf  , "*** ERROR: -f argument %s overflows integer field.\n", stFlag);
						fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
					}
				}
				else
				{
					sprintf(cbuf  , "*** INFO: Non-numeric character encountered in -f argument %s ... using default exponent for whatever FFT length is supplied.\n", stFlag);
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
			goto MLUCAS_HELP;
		}
	}	/* end of command-line-argument processing while() loop */

	if(!modType)
	{
		modType = MODULUS_TYPE_MERSENNE;
	}

	/* If user specified FFT length but no exponent, get default Mersenne exponent for that FFT length: */
	if(modType == MODULUS_TYPE_MERSENNE)
	{
		/* Special case of user forcing a non-default FFT length for an exponent in the worktodo.ini file -
		this assumes the arglist exponent matches the first entry in worktodo.ini and the specified FFT length
		already has a best-radix-set entry in the mlucas.cfg file:
		*/
		if (argc == 5 && MersVec[start].exponent && MersVec[start].fftLength )
		{
			if((retVal = ernstMain(MODULUS_TYPE_MERSENNE,TEST_TYPE_PRIMALITY,MersVec[start].exponent,MersVec[start].fftLength,0,0,0,FALSE,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime)) != 0)
			{
				printMlucasErrCode(retVal);
			}
		}
		else if(MersVec[start].exponent == 0)
		{
			i = MersVec[start].fftLength;
			ASSERT(HERE, i > 0                  ,"Mlucas.c: i > 0                  ");
			ASSERT(HERE, i <=MAX_FFT_LENGTH_IN_K,"Mlucas.c: i <=MAX_FFT_LENGTH_IN_K");

			/* If the FFT length is not represented in MersVec[],
			find the nearest prime <= given_N_get_maxP(FFT length):
			*/
			if(i < MersVec[0].fftLength || i > MersVec[numTest-1].fftLength)
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
			else	/* Find the corresponding entry of MersVec: */
			{
				for(lo = 0; lo < numTest; lo++)
				{
					if(MersVec[lo].fftLength == i)
					{
						start = lo; finish = start+1;
						break;
					}
				}
			}
		}
		/* If user specified exponent but no FFT length, get default FFT length for that exponent: */
		else if(MersVec[numTest].exponent && (MersVec[numTest].fftLength == 0))
		{
			MersVec[numTest].fftLength = get_default_fft_length(MersVec[numTest].exponent);
		}
	}
	else if(modType == MODULUS_TYPE_FERMAT)
	{
		if(FermVec[start].exponent == 0)
		{
			i = FermVec[start].fftLength;
			ASSERT(HERE, i > 0                  ,"Mlucas.c: i > 0                  ");
			ASSERT(HERE, i <=MAX_FFT_LENGTH_IN_K,"Mlucas.c: i <=MAX_FFT_LENGTH_IN_K");

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
		/* If user specified exponent but no FFT length, get default power-of-2 FFT length for that exponent: */
		else if(findex && (FermVec[numFerm].fftLength == 0))
		{
			FermVec[numFerm].fftLength = (uint32)1 << (findex - FermArrayIdxOffset);	/* Default is 16 bits per word, #kblocks = (f/16)/1024 */
		/*	FermVec[numFerm].fftLength = get_default_fft_length((uint32)1 << findex);	// Default is 16 bits per word */

			if(findex <= FermVec[numFerm-1].exponent)	/* Find the corresponding entry of FermVec: */
			{
				start = (findex - FermArrayIdxOffset); finish = start+1;
			}
		}
		/* User specified both exponent and FFT length: */
		else if(0)
		{
			if(findex <= FermVec[numFerm-1].exponent)	/* Find the corresponding entry of FermVec: */
			{
				start = (findex - FermArrayIdxOffset); finish = start+1;
			}
		}
	}
	else
	{
		ASSERT(HERE, 0,"modType not recognized!");
	}


TIMING_TEST_LOOP:

	fprintf(stderr, "\n           Mlucas selftest running.....\n\n");
	/* We have precomputed 100 and 1000-iteration residues for the predefined self-test exponents: */
	if( (userSetExponent && (modType == MODULUS_TYPE_MERSENNE)) || ((iters != 100) && (iters != 1000)) )
	{
		fprintf(stderr, "\n**** You will need to manually verify that the Res64s output");
		fprintf(stderr, "\n**** for this user-set p match for all FFT radix combinations!!\n\n\n");
	}
	fprintf(stderr, "/****************************************************************************/\n\n");

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
	{
		strcpy(CONFIGFILE,"fermat.cfg");
	}
	else	/* For now, anything else gets done using the Mersenne-mod .cfg file */
	{
		strcpy(CONFIGFILE,"mlucas.cfg");
	}

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
	FILE_ACCESS_MODE[1] = FILE_ACCESS_UPDATE;	/* Include this in fopen, so can
													(a) check whether it exists and is writeable;
													(b) if so, read the first line;
												in one fell swoop. */
	fp = fopen(CONFIGFILE,FILE_ACCESS_MODE);
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
		new_data = FALSE;
		Res64   = 0ull;
		Res36m1 = 0ull;
		Res35m1 = 0ull;

		/* If it's a self-test [i.e. timing test] and user hasn't specified #iters, set to default: */
		if(selfTest && !iters)
		{
		if(NTHREADS > 4)
			iters = 1000;
		else
			iters = 100;
		}

		if(iters == 100 || iters == 1000)
		{
			mvec_res_t_idx = NINT( log((double)iters)/log(10.) ) - 2;	/* log10(iters) - 2, use slower NINT rather than DNINT here since latter needs correct rounding mode */
			ASSERT(HERE, mvec_res_t_idx < 2,"main: mvec_res_t_idx out of range!");

			if( (modType == MODULUS_TYPE_MERSENNE && MersVec[  xNum].res_t[mvec_res_t_idx].sh0 == 0)
			 || (modType == MODULUS_TYPE_FERMAT   && FermVec[  xNum].res_t[mvec_res_t_idx].sh0 == 0) )	/* New self-test residue being computed */
			{
				new_data = TRUE;
				new_res.sh0 = Res64  ;
				new_res.sh1 = Res35m1;
				new_res.sh2 = Res36m1;
			}
		}

		/* Init best-radix-set and best-runtime for this FFT length: */
		runtime_best = 0.0;
		radix_best = -1;

		/* If user-specified radix set, do only that one: */
		if(radset >= 0)
			radix_set = radset;
		else
			radix_set = 0;

		roerr_min = 1.0;
		roerr_max = 0.0;

		if(modType == MODULUS_TYPE_MERSENNE)
		{
			iarg = MersVec[xNum].fftLength;
		}
		else if(modType == MODULUS_TYPE_FERMAT)
		{
			iarg = FermVec[xNum].fftLength;
		}

		while(get_fft_radices(iarg, radix_set, &NRADICES, RADIX_VEC, 10) == 0)	/* Try all the radix sets available for this FFT length. */
		{

			/* errCheck = 0/1: Error Checking is OFF/ON */
			if     (modType == MODULUS_TYPE_FERMAT)
			{
				Res64   = FermVec[  xNum].res_t[mvec_res_t_idx].sh0;
				Res35m1 = FermVec[  xNum].res_t[mvec_res_t_idx].sh1;
				Res36m1 = FermVec[  xNum].res_t[mvec_res_t_idx].sh2;
				retVal = ernstMain(modType,testType,(uint64)FermVec[xNum].exponent,iarg,radix_set,maxFFT,iters,errCheck,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime);
			}
			else if(modType == MODULUS_TYPE_MERSENNE)
			{
				Res64   = MersVec[  xNum].res_t[mvec_res_t_idx].sh0;
				Res35m1 = MersVec[  xNum].res_t[mvec_res_t_idx].sh1;
				Res36m1 = MersVec[  xNum].res_t[mvec_res_t_idx].sh2;
				retVal = ernstMain(modType,testType,(uint64)MersVec[xNum].exponent,iarg,radix_set,maxFFT,iters,errCheck,&Res64,&Res35m1,&Res36m1,scrnFlag,&runtime);
			}
			else
				ASSERT(HERE, 0,"0");

			/* Bzzzzzzzzzzzt!!! That answer is incorrect. The penalty is death... */
			if(retVal)
			{
				printMlucasErrCode(retVal);
				if( !userSetExponent && ((iters == 100) || (iters == 1000)) )
				{
					fprintf(stderr, "Error detected - this radix set will not be used.\n");
				}
				/* If user-specified radix set, do only that one: */
				if(radset >= 0)
					goto DONE;

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
				}
				else if				/* All subsequent radix sets must match to produce a consensus result: */
				(
					new_res.sh0 != Res64
				 ||	new_res.sh1 != Res35m1
				 ||	new_res.sh2 != Res36m1
				)
				{
					radix_best = -1;
					break;
				}
			}

			/* 16 Dec 2007: Added the (runtime != 0) here to workaround the valid-timing-test-but-runtime = 0
			issue on systems with round-to-nearest-second granularity of the clock() function:
			*/
			if( runtime_best == 0.0 || ((runtime != 0) && (runtime < runtime_best)) )
			{
				runtime_best = runtime;
				radix_best = radix_set;
			}

			if(MME < roerr_min)
				roerr_min = MME;

			if(MME > roerr_max)
				roerr_max = MME;

			fprintf(stderr, "\n");
			++radix_set;
		}

		/* If get no successful reference-Res64-matching results, inform user: */
		if(radix_best < 0 || runtime_best == 0.0)
		{
			sprintf(cbuf  , "WARNING: No valid best-radix-set found for FFT length %u K.\n",iarg << 10);
			fprintf(stderr,"%s", cbuf);
		}
		/* If get a nonzero best-runtime, write the corresponding radix set index to the .cfg file: */
		else
		{
			/* Divide by the number of iterations done in the self-test: */
			tdiff = runtime_best/iters;

			fp = fopen(CONFIGFILE,FILE_ACCESS_MODE);
			if(!fp)
			{
				sprintf(cbuf  , "INFO: Unable to open %s file in %s mode ... \n", CONFIGFILE, FILE_ACCESS_MODE);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}

			/* Put code version on line 1.
			Only want to do this once; subsequent fopen/fprintf are in append mode:
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
			fprintf(fp, "%10d  sec/iter = %8.3f  ROE[min,max] = [%10.9f, %10.9f]  radices = ", iarg, tdiff, roerr_min, roerr_max);
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


MLUCAS_HELP:
	fprintf(stderr, " Mlucas command line options:\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " Symbol and abbreviation key:\n");
	fprintf(stderr, "       <CR> :  carriage return\n");
	fprintf(stderr, "       []   :  encloses optional arguments\n");
	fprintf(stderr, "       {}   :  denotes user-supplied numerical arguments of the type noted.\n");
	fprintf(stderr, "              ({int} means nonnegative integer, {+int} = positive int, {float} = float.)\n");
	fprintf(stderr, "  -argument :  Vertical stacking indicates argument short 'nickname' options,\n");
	fprintf(stderr, "  -arg      :  e.g. in this example '-arg' can be used in place of '-argument'.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " Supported arguments:\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " <CR>        Default mode: looks for a %s file in the local;\n", RANGEFILE);
	fprintf(stderr, "             directory; if none found, prompts for manual keyboard entry\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -h          Prints this help file and exits\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " *** NOTE: *** The following options may cause an mlucas.cfg file containing\n");
	fprintf(stderr, "     the optimal FFT radix set for the runlength(s) tested, to be created (if one\n");
	fprintf(stderr, "     did not exist previously) and/or appended with new timing data.\n");
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
	fprintf(stderr, " -fftlen {+int}   If {+int} is one of the available FFT lengths (in Kilodoubles), runs all\n");
	fprintf(stderr, "             all available FFT radices available at that length, unless the -radset\n");
	fprintf(stderr, "             flag is also invoked (see below for details). If -fftlen is invoked\n");
	fprintf(stderr, "             with either the -m or -f flag, the self-tests will do either 100 iterations\n");
	fprintf(stderr, "             of a Lucas-Lehmer test (-m) or Pe'pin test (-f) on the user-specified Mersenne\n");
	fprintf(stderr, "             or Fermat number. If no user-set exponent is invoked, does 100 LL-test iterations\n");
	fprintf(stderr, "             using the default self-test Mersenne or Fermat exponent for that FFT length\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "             Use this to find the optimal radix set for a single FFT length on your hardware.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -radset {int}    Specific index of a set of complex FFT radices to use, based on\n");
	fprintf(stderr, "             the big select table in the function get_fft_radices(). Requires a\n");
	fprintf(stderr, "             supported value of -fftlen to be specified.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -m [{+int}] Performs a Lucas-Lehmer primality test of the Mersenne number M(int) = 2^int - 1,\n");
	fprintf(stderr, "             where int must be an odd prime. If -iters is also invoked, this indicates a timing test.\n");
	fprintf(stderr, "             and requires suitable added arguments (-fftlen, -iters, optionally -radset) to be supplied.\n");
	fprintf(stderr, "\n        If the -fftlen option (and optionally -radset) is also invoked but -iters is not,");
	fprintf(stderr, "             the program first checks the first line of the worktodo.ini file voked, the tests use all sets of\n");
	fprintf(stderr, "             FFT radices available at the specified FFT length.\n");
	fprintf(stderr, "                If the -fftlen option is not invoked, the tests use all sets of\n");
	fprintf(stderr, "             FFT radices available at that exponent's default FFT length.\n");
	fprintf(stderr, "                Use this to find the optimal radix set for a single given Mersenne\n");
	fprintf(stderr, "             exponent on your hardware, similarly to the -fftlen option.\n");
	fprintf(stderr, "                Performs 100 iterations, or as many as specified via the -iters flag.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -f {int}    Performs a base-3 Pe'pin test on the Fermat number F(num) = 2^(2^num) + 1.\n");
	fprintf(stderr, "                If desired this can be invoked together with the -fftlen option.\n");
	fprintf(stderr, "             as for the Mersenne-number self-tests (see notes about the -m flag;\n");
	fprintf(stderr, "             note that not all FFT lengths supported for -m are available for -f).\n");
	fprintf(stderr, "             Optimal radix sets and timings are written to a fermat.cfg file.\n");
	fprintf(stderr, "                Performs 100 iterations, or as many as specified via the -iters flag.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -rocheck    Turn per-iteration roundoff error checking on/off.\n");
	fprintf(stderr, "             Program will typically run 3-10%% slower with checking enabled.\n");
	fprintf(stderr, "             Default is OFF, i.e. run with roundoff error checking disabled.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -iters {int}   Do {int} self-test iterations of the type determined by the\n");
	fprintf(stderr, "             modulus-related options (-s/-m = Lucas-Lehmer test iterations with\n");
	fprintf(stderr, "             initial seed 4, -f = Pe'pin-test squarings with initial seed 3.\n");
	fprintf(stderr, "             Default is 100 iterations.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " -nthread {int}   For multithread-enabled builds, run with this many threads.\n");
	return(0);
}

/******************/

/* Given an input line (presumed to be the first line of the mlucas.cfg or fermat.cfg file,
which is presumed to contain the program version that was used to generate said .cfg file),
returns nonzero (which should be taken as a proxy for TRUE) if the .cfg file needs to be updated
via a new round of timing tests, FALSE otherwise.
*/
int	cfgNeedsUpdating(char*in_line)
{
	/* For the foreseeable future, version numbers will be of form x.yz, with x a 1-digit integer,
	y an int having 3 digits or less, and z an optional alphabetic suffix. We choose these such that
	retiming is only necessary between versions that differ in x or in the leading digit of y,
	i.e. we'd need to regenerate .cfg files in going from version 2.9 to 3.0 or from 3.0 to 3.141,
	but not between versions 3.0x and 3.0g or between 3.1 and 3.141.

	This reduces the "retime?" decision to a simple comparison of the leading 3 characters of the version string.
	*/
	return STRNEQN(in_line, VERSION, 3);
}

/******************/

void	printMlucasErrCode(int retVal)
{
	/* These need to be kept updated to match the #defines in Mdata.h: */
	const char *err_code[] = {
		"ERR_INCORRECT_RES64		",
		"ERR_RADIX0_UNAVAILABLE		",
		"ERR_RADIXSET_UNAVAILABLE	",
		"ERR_TESTTYPE_UNSUPPORTED	",
		"ERR_EXPONENT_ILLEGAL		",
		"ERR_FFTLENGTH_ILLEGAL		",
		"ERR_ECHECK_NOTBOOL			",
		"ERR_TESTITERS_OUTOFRANGE	",
		"ERR_ROUNDOFF				",
		"ERR_CARRY					",
		"ERR_RUN_SELFTEST_FORLENGTH	"
	};

	/* Error type indicated by lowest byte: */
	uint32 i = (retVal & 0xff);

	if(i == 0)
		fprintf(stderr, "\n Return with error code 0 - no errors.\n");
	else if(i <= ERR_MAX)
		fprintf(stderr, "\n Return with code %s\n", err_code[i-1]);
	else
		fprintf(stderr, "\n Return with unknown error error code %u - suggest running under a debugger.\n\n",(uint32)retVal);

	/* High bytes should only be nonzero if low byte == ERR_RUN_SELFTEST_FORLENGTH: */
	if((retVal>>8) != 0)
	{
		ASSERT(HERE, i==ERR_RUN_SELFTEST_FORLENGTH, "High bytes should only be nonzero if low byte == ERR_RUN_SELFTEST_FORLENGTH!");
	}
}

/******************/

/* Read the data contained in an old-style (pre-V3.0) Mlucas LL test savefile
(as this is not in a platform-independent form, requires the savefile to have
been generated on the the same type of system as the binary was built for)
and then rename the old savefile file p{exponent}.bak .
*/
#if 0
/*
	EWM: For some reason the old-style savefile write combined with this gave weird extra bytes
	in front of the exponent - simply wasn't worth any more time trying to debug, just switch to bytewise
	savefile read/write and ask users to finish any pre-v3 assignments using their v2.8 build.
*/
	#error obsolete code - remove from build!

int 	convert_LL_savefiles(uint64 psave, FILE*fp, uint32*ilo, uint32 ndim, int32 arr_tmp[], double a[])
{
	uint32 p,i,j,j1,n;
	double sum1, sum2;
	char BACKUPFILE[STR_MAX_LEN];

	if(fp)
	{
		fprintf(stderr, " INFO: restart file %s found...reading...\n",RESTARTFILE);

		i = fread(&p, sizeof(int), 1, fp);
		if(i != 1)		{ fprintf(stderr, " ERROR: Error reading p.\n")										; fclose(fp);	fp = 0x0; return FALSE; }
		/* Check that the restart-file exponent matches the passed one: */
		if(p != psave)	{ fprintf(stderr, " ERROR: savefile p does not match passed p!\n")					; fclose(fp);	fp = 0x0; return FALSE; }
		if(feof(fp))	{ fprintf(stderr, " ERROR: End-of-file encountered while attempting to read p.\n")	; fclose(fp);	fp = 0x0; return FALSE; }

		i = fread(ilo, sizeof(int), 1, fp);
		if(i != 1)		{ fprintf(stderr, " ERROR: Error reading ilo.\n")									; fclose(fp);	fp = 0x0; return FALSE; }
		if(feof(fp))	{ fprintf(stderr, " ERROR: End-of-file encountered while attempting to read ilo.\n")	; fclose(fp);	fp = 0x0; return FALSE; }

		i = fread(&n, sizeof(int), 1, fp);
		if(i != 1)		{ fprintf(stderr, " ERROR: Error reading n.\n")										; fclose(fp);	fp = 0x0; return FALSE; }
		if(feof(fp))	{ fprintf(stderr, " ERROR: End-of-file encountered while attempting to read n.\n")	; fclose(fp);	fp = 0x0; return FALSE; }
		if(n != ndim)	{ fprintf(stderr, " ERROR: n != ndim in convert_LL_savefiles\n")					; fclose(fp);	fp = 0x0; return FALSE; }

		i = fread(arr_tmp, sizeof(int), n, fp);		/* Read integer residues...	*/
		if(i != n)		{ fprintf(stderr, "ERROR: Error reading integer residue array.\n")										; fclose(fp);	fp = 0x0; return FALSE; }
		if(feof(fp))	{ fprintf(stderr, "ERROR: End-of-file encountered while attempting to read integer residue array.\n")	; fclose(fp);	fp = 0x0; return FALSE; }

		i = fread(&sum1, sizeof(double), 1, fp);		/* ...and checksum...	*/
		if(i != 1)		{ fprintf(stderr, "ERROR: Error reading checksum1.\n")													; fclose(fp);	fp = 0x0; return FALSE; }
		if(feof(fp))	{ fprintf(stderr, "ERROR: End-of-file encountered while attempting to read checksum1.\n")				; fclose(fp);	fp = 0x0; return FALSE; }

		fclose(fp);	fp = 0x0;

		sum2=0.0;
		for(j=0; j < n; j++)	/* ...then move elements of residue into their normal doubles slots.	*/
		{
		#ifdef USE_AVX
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			a[j1]= arr_tmp[i];
			sum2 += a[j1];
		}
		/* Don't deallocate arr_tmp here, since we'll need it later for savefile writes. */
		/*...we can re-use I/O unit number 1 here since if program execution gets to this point, the savefile previously opened on unit 1 will have been closed.	*/

		if(sum1 != sum2)
		{
										        fprintf(stderr,"ERROR: Checksum error for restart file%s\n", RESTARTFILE);
			fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Checksum error for restart file%s\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
			fp = fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Checksum error for restart file%s\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
			return FALSE;
		}
	}
	else
	{
									        fprintf(stderr,"ERROR: Unable to open restart file %s in convert_LL_savefiles\n", RESTARTFILE);
		fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Unable to open restart file %s in convert_LL_savefiles\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
		fp = fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Unable to open restart file %s in convert_LL_savefiles\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
		return FALSE;
	}

	/* now rename the old savefile by appending a '.bak' to the name: */
	strcpy(BACKUPFILE, RESTARTFILE);
	strcat(BACKUPFILE,".bak");
	if(!rename(RESTARTFILE, BACKUPFILE))
	{
		sprintf(cbuf  ,"FATAL: Unable to rename the savefile %s ==> %s\n",RESTARTFILE, BACKUPFILE);
		fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf); }
		fp = fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf); }
		fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

	return TRUE;
}
#endif

/******************/

/* Pair of functions to Read/Write full-length residue data in bytewise format
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

	s:  The number of squarings of 3 {EWM: or LL iterations, or whatever}
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

The input arguments Res64, Res35m1, Res36m1 are 3 checksums (residue modulo 2^64, 2^35-1, 2^36-1)
intended to be compared to the same gotten via subsequent processing (e.g. conversion to floating-
point form) of the residue array itself. (These are essentially the Selfridge-Hurwitz residues
commonly used in Fermat-number primality testing, but with the mod-2^36 component widened to 64 bits).
*/
/*** READ: Assumes a valid file pointer has been gotten via a call of the form
fp = fopen(RESTARTFILE,"rb");
***/
int	read_ppm1_savefiles(uint64 p, FILE*fp, uint32*ilo, uint8 arr_tmp[], uint64*Res64, uint64*Res35m1, uint64*Res36m1)
{
	uint32 i;
	uint32 nbytes;
	uint64 nsquares= 0;

	ASSERT(HERE, !(p >> 32), "read_ppm1_savefiles: p must be 32-bit or less!");	/* Future versions will need to loosen this p < 2^32 restriction: */

	if(!file_valid(fp))
	{
		sprintf(cbuf, "read_ppm1_savefiles: File pointer invalid for read!\n");
		return FALSE;
	}
	fprintf(stderr, " INFO: restart file %s found...reading...\n",RESTARTFILE);
	/* t: */
	if((i = fgetc(fp)) != TEST_TYPE)
	{
		sprintf(cbuf, "read_ppm1_savefiles: TEST_TYPE != fgetc(fp)\n");
	//	return FALSE;
	}
	/* m: */
	if((i = fgetc(fp)) != MODULUS_TYPE)
	{
		// For some reason, this fubared in my rerun-final-F25-iterations-from-33.55m (fgetc = 176, MODULUS_TYPE = 3)
		// but residue OK, so emit error msg but allow execution past it:
		sprintf(cbuf, "ERROR: read_ppm1_savefiles: MODULUS_TYPE != fgetc(fp)\n");
	//	return FALSE;
	}
	/* s: */
	for(nbytes = 0; nbytes < 8; nbytes++)
	{
		i = fgetc(fp);
		nsquares += (uint64)i << (8*nbytes);
	}
	/* For now, just make sure nsquares < 2^32 and copy to ilo: */
	if(nsquares >= p)
	{
		sprintf(cbuf,"read_ppm1_savefiles: nsquares = %llu out of range, should be < p = %llu\n", nsquares, p);
	//	return FALSE;
	}
	*ilo = nsquares;

	/* r: */
	/* Set the expected number of residue bytes, depending on the modulus: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		nbytes = (p + 7)/8;
		TRANSFORM_TYPE = REAL_WRAPPER;
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, p % 8 == 0,"read_ppm1_savefiles: p % 8 == 0");
		nbytes = (p/8) + 1;
		TRANSFORM_TYPE = RIGHT_ANGLE;
	}

	i = fread(arr_tmp, sizeof(char), nbytes, fp);		/* Read bytewise residue...	*/
	if(i != nbytes)	{ sprintf(cbuf, "read_ppm1_savefiles: Error reading bytewise residue array.\n")										; return FALSE; }
	if(ferror(fp))	{ sprintf(cbuf, "read_ppm1_savefiles: Unknown Error reading bytewise residue array.\n")								; return FALSE; }
	if(feof(fp))	{ sprintf(cbuf, "read_ppm1_savefiles: End-of-file encountered while attempting to read bytewise residue array.\n")	; return FALSE; }

	/* 8 bytes for Res64: */
	*Res64 = 0;
	for(nbytes = 0; nbytes < 8; nbytes++)
	{
		i = fgetc(fp);
		*Res64 += (uint64)i << (8*nbytes);
	}
	/* 5 bytes for Res35m1: */
	*Res35m1 = 0;
	for(nbytes = 0; nbytes < 5; nbytes++)
	{
		i = fgetc(fp);
		*Res35m1 += (uint64)i << (8*nbytes);
	}
	ASSERT(HERE, *Res35m1 <= 0x00000007FFFFFFFFull,"read_ppm1_savefiles: *Res35m1 <= 0x00000007ffffffff");
	/* 5 bytes for Res36m1: */
	*Res36m1 = 0;
	for(nbytes = 0; nbytes < 5; nbytes++)
	{
		i = fgetc(fp);
		*Res36m1 += (uint64)i << (8*nbytes);
	}
	ASSERT(HERE, *Res36m1 <= 0x0000000FFFFFFFFFull,"read_ppm1_savefiles: *Res36m1 <= 0x0000000fffffffff");
	/* Don't deallocate arr_tmp here, since we'll need it later for savefile writes. */
	return TRUE;
}

/*** WRITE:  Assumes a valid file pointer has been gotten via a call of the form
fp = fopen(RESTARTFILE,"wb");
***/
#if OLD_FORMAT
void write_ppm1_savefiles(uint64 p, FILE*fp, uint32 ihi, uint32 n, int32 arr_tmp[], double a[])
{
	double sum1 = 0.0;
	uint32 i,j;
	uint32 p32 = (uint32)p;
#else
void write_ppm1_savefiles(uint64 p, FILE*fp, uint32 ihi, uint8 arr_tmp[], uint64 Res64, uint64 Res35m1, uint64 Res36m1)
{
	uint32 i;
	uint32 nbytes = 0;
	uint64 itmp64;
#endif
	if(!file_valid(fp))
	{
		ASSERT(HERE,0,"write_ppm1_savefiles: File pointer invalid for write!");
	}

	/* See the function read_ppm1_savefiles() for the file format here: */
#if OLD_FORMAT
	#error This segment of code is obsolete!
	/* Obsolete (platform-dependent) savefile format - invoke this code only to test convert_LL_savefiles functionality: */
	for(i = 0; i < n; i++)
	{
		j = i + ( (i >> DAT_BITS) << PAD_BITS );
		arr_tmp[i]=(int)a[j];
		sum1 += a[j];
	}
/********************************************
#error This fails on fread of the same file!
********************************************/
	i = fwrite(&p32, sizeof(int), 1, fp);
	ASSERT(HERE, i == 1, "write_ppm1_savefiles: Error writing p32.");

	i = fwrite(&ihi, sizeof(int), 1, fp);
	ASSERT(HERE, i == 1, "write_ppm1_savefiles: Error writing ihi.");

	i = fwrite(&n, sizeof(int), 1, fp);
	ASSERT(HERE, i == 1, "write_ppm1_savefiles: Error writing n.");

	i = fwrite(arr_tmp, sizeof(int), n, fp);		/* Write integer residue...	*/
	ASSERT(HERE, i == n, "write_ppm1_savefiles: Error writing residue.");

	i = fwrite(&sum1, sizeof(double), 1, fp);		/* ...and checksum.	*/
	ASSERT(HERE, i == 1, "write_ppm1_savefiles: Error writing checksum1.");
#else
	/* Bytewise (platform-independent) savefile format: */
	/* t: */
	fputc(TEST_TYPE, fp);
	/* m: */
	fputc(MODULUS_TYPE, fp);
	/* s: */
	itmp64 = ihi;
	for(i = 0; i < 64; i+=8)
	{
		fputc((itmp64 >> i) & 0xff, fp);
	}

	/* Set the expected number of residue bytes, depending on the modulus: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		nbytes = (p + 7)/8;
		TRANSFORM_TYPE = REAL_WRAPPER;
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, p % 8 == 0,"write_ppm1_savefiles: p % 8 == 0");
		nbytes = (p/8) + 1;
		TRANSFORM_TYPE = RIGHT_ANGLE;
	}

	/* Write bytewise residue r... */
	i = fwrite(arr_tmp, sizeof(char), nbytes, fp);
	if(i != nbytes)
	{
		fclose(fp);
		sprintf(cbuf, "Error writing residue to restart file %s.\n",RESTARTFILE);
									        fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
		fp = fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"%s",cbuf);	fclose(fp); fp = 0x0; }
		ASSERT(HERE, 0,"");
	}
	/* ...and checksums:	*/
	/* Res64: */
	for(i = 0; i < 64; i+=8)
	{
		fputc((int)(Res64 >> i) & 0xff, fp);
	}
	/* Res35m1: */
	for(i = 0; i < 40; i+=8)
	{
		fputc((int)(Res35m1 >> i) & 0xff, fp);
	}
	/* Res36m1: */
	for(i = 0; i < 40; i+=8)
	{
		fputc((int)(Res36m1 >> i) & 0xff, fp);
	}
#endif
}

/*********************/

/*
Function to take an n-digit Mersenne or Fermat-mod residue in bytewise savefile form,
convert it to balanced-digit floating-point form, calculate the 3 Selfridge-Hurwitz residues
and perform checksum validation - it is assumed that the 3 residues will have been separately
read from e.g. a savefile, and these reference values are passed as const ints.

In the Mersenne-mod case the residue digits are stored
consecutively in the a[] array.

In the Fermat-mod case the digits are arranged in (j,j+n/2)
(i.e. right-angle transform) order.
*/
int 	convert_res_bytewise_FP(const uint8 arr_tmp[], double a[], int n, const uint64 p, const uint64 Res64, const uint64 Res35m1, const uint64 Res36m1)
{
	char hex_res[17];
	uint32 nbytes;
	uint64 nbits,sum0,sum1,sum2;
	int bimodn,curr_char,findex,ii,j = 0,j1 = 0,k,pass,rbits;
	int bw,sw,bits[2];
	uint64 base[2];	/* Assume base may be > 2^32 (e.g. for mixed FFT/FGT)
						but < 2^53, i.e. fits in a double */
	int64 itmp,cy;
	uint64 curr_word, curr_wd64;
	int pow2_fft;
	FILE *fp;

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

	/* Vector length a power of 2? */
	if((n >> trailz32(n)) == 1)
	{
		pow2_fft = TRUE;
	}
	else
	{
		pow2_fft = FALSE;
	}

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
		#ifdef USE_AVX
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
				ASSERT(HERE, itmp < (1<<rbits),"convert_res_bytewise_FP: itmp >= 2^rbits!");

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
			if(itmp > 0)
			{
				cy = itmp >> bits[ii];
				itmp -= (cy << bits[ii]);
				if(itmp > (base[ii]>>1))
				{
					itmp -= base[ii];
					cy += 1;
				}
			}
			else
			{
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
		#ifdef USE_AVX
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
			if(itmp > 0)
			{
				cy = itmp >> bits[ii];
				itmp -= (cy << bits[ii]);
				if(itmp > (base[ii]>>1))
				{
					itmp -= base[ii];
					cy += 1;
				}
			}
			else
			{
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
	{
		a[0] += cy;
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		a[0] -= cy;
	}
	else
	{
		ASSERT(HERE, 0,"Illegal modulus type!");
	}

	/* Checksum validation: */
	sum0 = res64(a,n,p,&j,hex_res);
	if(Res64 != sum0)
	{
									        fprintf(stderr,"ERROR: Res64 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE);
		fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Res64 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
		fp = fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Res64 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
		return FALSE;
	}

	resSH(a,n,p,&sum1,&sum2);
	if(Res35m1 != sum1)
	{
									        fprintf(stderr,"ERROR: Res35m1 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE);
		fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Res35m1 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
		fp = fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Res35m1 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
		return FALSE;
	}
	if(Res36m1 != sum2)
	{
									        fprintf(stderr,"ERROR: Res36m1 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE);
		fp = fopen(   OFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Res36m1 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
		fp = fopen(STATFILE,"a");	if(fp){ fprintf(	fp,"ERROR: Res36m1 Checksum error for restart file %s in convert_res_bytewise_FP\n", RESTARTFILE); fclose(fp);	fp = 0x0; }
		return FALSE;
	}

	return TRUE;
}

/*********************/

/*
Function to take an n-digit Mersenne or Fermat-mod residue in balanced-digit
floating-point form, convert it to bytewise form, and generate Selfridge-Hurwitz checksums.

In the Mersenne-mod case the residue digits are assumed to be stored
consecutively in the a[] array.

In the Fermat-mod case the digits are assumed to be arranged in (j,j+n/2)
(i.e. right-angle transform) order.

***NOTE:*** To convert a generic multiword int (i.e. not w.r.to a particular
modulus) from/to balanced-digit fixed-base floating-point form and uint64
form, use the mi64_cvt_double_uint64() and mi64_cvt_uint64_double() functions in mi64.c .
*/
void	convert_res_FP_bytewise(const double a[], uint8 arr_tmp[], int n, const uint64 p, uint64*Res64, uint64*Res35m1, uint64*Res36m1)
{
	uint64 nbits;
	int bimodn,curr_bits,curr_char,cy,findex,ii,j,j1,k,pass,shift,rbits,msw_lt0;
	int bw,sw,bits[2];
	uint64 base[2];	/* Assume base may be > 2^32 (e.g. for mixed FFT/FGT)
						but < 2^53, i.e. fits in a double */
	int64 itmp;
	uint64 curr_word, curr_wd64, mod1=0, mod2=0;
	const uint64 two35m1 = (uint64)0x00000007FFFFFFFFull, two36m1 = (uint64)0x0000000FFFFFFFFFull;	/* 2^35,36-1 */
	int pow2_fft;

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
	{
		ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER,"convert_res_FP_bytewise: TRANSFORM_TYPE == REAL_WRAPPER");
	}
	else
	{
		ASSERT(HERE, 0,"Illegal modulus type!");
	}

	/* Vector length a power of 2? */
	if((n >> trailz32(n)) == 1)
	{
		pow2_fft = TRUE;
	}
	else
	{
		pow2_fft = FALSE;
	}

	bits[0] = p/n;
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
	/*
	If most-significant digit in the balanced-representation form is < 0, add the modulus to the residue.
	For Mersenne (2^p-1) and Fermat (2^p+1) moduli, combine this with the normalize-to-nonnegative-digit
	step (which we do in any event) by simply initializing the carry into the latter to -1 or +1, respectively.
	In this case we expect the carryout of the normalization loop to = -1, indicating that the MS word has
	been accordingly normalized during the loop - add an assertion check to that effect.
	*/
	cy=0;		/* init carry.	*/
	msw_lt0 = 0;

	for(j=n-1; j >= 0; j -= TRANSFORM_TYPE)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

		if(a[j1]!= 0.0)
		{
			if(a[j1]< 0.0)
			{
				msw_lt0 = 1;	/* MS word was < 0 prior to normalization */

				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					cy = -1;
				}
				else if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
				{
					cy = +1;
				}
			}
			break;
		}
	}

	/*...Now form the SH residues, converting to positive-digit form along the way...	*/

	curr_char = 0;	/* Current byte to be written to in the arr_tmp array */
	nbits = 0;		/* Total bits accumulated so far in the residue	*/
	rbits = 0;		/* # of Upper bits left over from processing of previous word	*/
	curr_wd64 = 0;	/*      Upper bits left over from processing of previous word	*/

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
		#ifdef USE_AVX
			j1 = (j & mask02) + br8[j&7];
		#elif defined(USE_SSE2)
			j1 = (j & mask01) + br4[j&3];
		#else
			j1 = j;
		#endif
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			itmp = (int64)(a[j1]+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0)			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
			{
				itmp += (base[ii]);
				cy = -1;
			}
			else
			{
				cy = 0;
			}

			ASSERT(HERE, itmp >= 0,"convert_res_FP_bytewise: itmp >= 0");

		/* Mod-(2^35-1) residue: */
			curr_word = (uint64)itmp;
			/* Current word must be left-shifted by nbits%35 before folding into residue: */
			shift = (nbits%35);
			mod1 += (curr_word << shift) & two35m1;
			mod1 = (mod1 >> 35) + (mod1 & two35m1);
			curr_word >>= (35-shift);
			while(curr_word)
			{
				mod1 += curr_word & two35m1;
				mod1 = (mod1 >> 35) + (mod1 & two35m1);
				curr_word >>= 35;
			}
		/* Mod-(2^36-1) residue: */
			curr_word = (uint64)itmp;
			/* Current word must be left-shifted by nbits%36 before folding into residue: */
			shift = (nbits%36);
			mod2 += (curr_word << shift) & two36m1;
			mod2 = (mod2 >> 36) + (mod2 & two36m1);
			curr_word >>= (36-shift);
			while(curr_word)
			{
				mod2 += curr_word & two36m1;
				mod2 = (mod2 >> 36) + (mod2 & two36m1);
				curr_word >>= 36;
			}

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
	  for(pass = 0; pass <=1; pass++)
	  {
		bimodn = n;

		for(j = pass; j < n; j += 2)
		{
		#ifdef USE_AVX
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

			itmp = (int64)(a[j1]+ cy);	/* current digit in int64 form, subtracting any borrow from the previous digit.	*/
			if(itmp < 0)			/* If current digit < 0, add the current base and set carry into next-higher digit = -1	*/
			{
				itmp += (base[ii]);
				cy = -1;
			}
			else
			{
				cy = 0;
			}

			ASSERT(HERE, itmp >= 0,"convert_res_FP_bytewise: itmp >= 0");

		/* Mod-(2^35-1) residue: */
			curr_word = (uint64)itmp;
			/* Current word must be left-shifted by nbits%35 before folding into residue: */
			shift = (nbits%35);
			mod1 += (curr_word << shift) & two35m1;
			mod1 = (mod1 >> 35) + (mod1 & two35m1);
			curr_word >>= (35-shift);
			while(curr_word)
			{
				mod1 += curr_word & two35m1;
				mod1 = (mod1 >> 35) + (mod1 & two35m1);
				curr_word >>= 35;
			}
		/* Mod-(2^36-1) residue: */
			curr_word = (uint64)itmp;
			/* Current word must be left-shifted by nbits%36 before folding into residue: */
			shift = (nbits%36);
			mod2 += (curr_word << shift) & two36m1;
			mod2 = (mod2 >> 36) + (mod2 & two36m1);
			curr_word >>= (36-shift);
			while(curr_word)
			{
				mod2 += curr_word & two36m1;
				mod2 = (mod2 >> 36) + (mod2 & two36m1);
				curr_word >>= 36;
			}

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

	/* Checksums: */
	*Res64 = 0;
	for(ii= 0, j = 0; ii< 64;ii+=8, j++)
	{
		*Res64 += (uint64)arr_tmp[j] << ii;
	}
	*Res35m1 = mod1;
	*Res36m1 = mod2;
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


/*********************/


/*********************/


/*********************/

void write_fft_debug_data(double a[], int jlo, int jhi)
{
	int j,j1;
	const char dbg_fname[] = "FFT_DEBUG.txt";
	ASSERT(HERE, dbg_file == 0x0, "dbg_file != 0x0 prior to fopen");
	dbg_file = fopen(dbg_fname, "a");
	ASSERT(HERE, dbg_file != 0x0, "Unable to open dbg_file!");
	fprintf(dbg_file, "RE_IM_STRIDE = %d\n", RE_IM_STRIDE);
	fprintf(dbg_file, "%s\n", cbuf);

    for(j=jlo; j < jhi; j += 2)
    {
	#ifdef USE_AVX
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

