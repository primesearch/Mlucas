/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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

#include "Mlucas.h"

#define RADIX 36	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#ifndef PFETCH_DIST
  #ifdef USE_AVX
	#define PFETCH_DIST	32	// This seems to work best on my Haswell, even though 64 bytes seems more logical in AVX mode
  #else
	#define PFETCH_DIST	32
  #endif
#endif

#ifdef MULTITHREAD
	#ifndef USE_PTHREAD
		#error Pthreads is only thread model currently supported!
	#endif
#endif

#ifdef USE_SSE2

	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */
	#if defined(USE_AVX) && !defined(ERR_CHECK_ALL)
		#error ERR_CHECK_ALL *required* for AVX-mode builds!
	#endif

	#define EPS 1e-10

  // For Mersenne-mod we need (16 [SSE2] or 64 [AVX]) + 4 added slots for the half_arr lookup tables.
  // Add relevant number (half_arr_offset36 + RADIX) to get required value of radix36_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset36 = 0xa5;	// + RADIX = 0xc9; Used for thread local-storage-integrity checking
	const int radix36_creals_in_local_store = 0x110;	// += 68 (= 0x44) and round up to nearest multiple of 4
	#if OS_BITS == 64
	  #define GCC_ASM_FULL_INLINE	0	// This *must* be #defined = 1 in AVX mode!
	#endif
  #else
	const int half_arr_offset36 = 0xae;	// + RADIX = 0xd2; Used for thread local-storage-integrity checking
	const int radix36_creals_in_local_store = 0xe8;	// += 20 (= 0x14) and round up to nearest multiple of 4
	#if OS_BITS == 64
	  // #define to either (if left undefined) use small-macro form below, or (if defined) to inline the fused macros as single big blob of asm (64-bit only):
	  #define GCC_ASM_FULL_INLINE	0	// sse2/core 2 timings better for small-macro form due to smaller obj-code size.
	#endif
  #endif

	// Only have a fully-fused macro for the DIT here, not enough timing gain to justify more work for a DIF:
	#if GCC_ASM_FULL_INLINE
		#if OS_BITS == 32
			#include "radix36_ditN_cy_dif1_gcc32.h"
		#else
			#include "radix36_ditN_cy_dif1_gcc64.h"
		#endif
	#else
		#include "sse2_macro.h"
		#include "radix09_sse_macro.h"
	#endif

#endif	// SSE2

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
		int iter;
		int tid;
		int ndivr;

		int khi;
		int i;
		int jstart;
		int jhi;
		int col;
		int co2;
		int co3;
		int sw;
		int nwt;

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
	#ifdef USE_SSE2
		vec_dbl *r00;
		vec_dbl *half_arr;
	#else
		double *r00;
		double *half_arr;
	#endif
		uint32 bjmodnini;
		int bjmodn0;
	// For large radix0 use thread-local arrays for DWT indices/carries - only caveat is these must be SIMD-aligned:
	#if GCC_EVER_GETS_ITS_ACT_TOGETHER_HERE
	/* Jan 2014: Bloody hell - turns out GCC uses __BIGGEST_ALIGNMENT__ = 16 on x86, which is too small to be useful for avx data!
		int bjmodn[RADIX] __attribute__ ((aligned (32)));
		double cy[RADIX] __attribute__ ((aligned (32)));
	*/
	#else
	// Thus, we are forced to resort to fugly hackage - add pad slots to a garbage-named struct-internal array along with
	// a pointer-to-be-inited-at-runtime, when we set ptr to the lowest-index array element having the desired alginment:
		double *cy;
		double cy_dat[RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	#endif
	};

#endif

/****************/

int radix36_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-36 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-36 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix36_ditN_cy_dif1";
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	const int sz_vd = sizeof(vec_dbl), sz_vd_m1 = sz_vd-1;
	// lg(sizeof(vec_dbl)):
  #ifdef USE_AVX
	const int l2_sz_vd = 5;
  #else
	const int l2_sz_vd = 4;
  #endif
#else
	const int sz_vd = sizeof(double), sz_vd_m1 = sz_vd-1;
#endif

	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32;
	static int poff[RADIX>>2];	// Store mults of p04 offset for loop control
#if defined(USE_SSE2) || !defined(MULTITHREAD)
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
#endif
	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX], *tptr;
	int *itmp;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #ifdef USE_AVX2
	// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
	vec_dbl *rad9_iptr[9], *rad9_optr[9];
   #endif
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #endif

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7, *add8;	/* Addresses into array sections */
  #endif

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	vec_dbl *tmp,*tm1,*tm2;	// Non-static utility ptrs
	static vec_dbl *two,*one, *cc1, *ss1, *cc2, *ss2, *cc3m1, *ss3, *cc4, *ss4, *max_err, *sse2_rnd, *half_arr
		,*r00,*r02,*r04,*r06,*r08,*r0a,*r0c,*r0e,*r0g
		,*r10,*r12,*r14,*r16,*r18,*r1a,*r1c,*r1e,*r1g
		,*r20,*r22,*r24,*r26,*r28,*r2a,*r2c,*r2e,*r2g
		,*r30,*r32,*r34,*r36,*r38,*r3a,*r3c,*r3e,*r3g
		,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r
		,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r
		,*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx

#endif

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy36_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double *addr;
	int bjmodn[RADIX];
	double rt,it,temp,frac,cy[RADIX], re;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodnini = 0x0,*_bjmodn[RADIX];
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,*_cy[RADIX];
	if(!_maxerr) {
		_cy[0] = 0x0;	// First of these used as an "already inited consts?" sentinel, must init = 0x0 at same time do so for non-array static ptrs
	}

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/RADIX in %s.\n",iter,func);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry: */

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the smallest
		power of 2 that is >= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else
		{
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> (i-1);
		}

		if(CY_THREADS > MAX_THREADS)
		{
		//	CY_THREADS = MAX_THREADS;
			fprintf(stderr,"WARN: CY_THREADS = %d exceeds number of cores = %d\n", CY_THREADS, MAX_THREADS);
		}
		if(CY_THREADS < NTHREADS)	{ WARN(HERE, "CY_THREADS < NTHREADS", "", 1); return(ERR_ASSERT); }
		if(!isPow2(CY_THREADS))		{ WARN(HERE, "CY_THREADS not a power of 2!", "", 1); return(ERR_ASSERT); }
		if(CY_THREADS > 1)
		{
			if(NDIVR    %CY_THREADS != 0) { WARN(HERE, "NDIVR    %CY_THREADS != 0", "", 1); return(ERR_ASSERT); }
			if(n_div_nwt%CY_THREADS != 0) { WARN(HERE, "n_div_nwt%CY_THREADS != 0", "", 1); return(ERR_ASSERT); }
		}

	  #ifdef USE_PTHREAD

		j = (uint32)sizeof(struct cy_thread_data_t);
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, j);

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#if 0//def OS_TYPE_MACOSX

			if(CY_THREADS > 1) {
				main_work_units = CY_THREADS/2;
				pool_work_units = CY_THREADS - main_work_units;
				ASSERT(HERE, 0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
				printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
			} else {
				main_work_units = 1;
				printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
			}

		#else

			pool_work_units = CY_THREADS;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

		#endif

		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);

	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		// int data:
			tdat[ithread].tid = ithread;
			tdat[ithread].ndivr = NDIVR;

			tdat[ithread].sw  = sw;
			tdat[ithread].nwt = nwt;

		// pointer data:
			tdat[ithread].arrdat = a;			/* Main data array */
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
			tdat[ithread].si  = si;

		// This array pointer must be set based on vec_dbl-sized alignment at runtime for each thread:
			for(l = 0; l < 4; l++) {
				if( ((uint32)&tdat[ithread].cy_dat[l] & sz_vd_m1) == 0 ) {
					tdat[ithread].cy = &tdat[ithread].cy_dat[l];
				//	fprintf(stderr,"%d-byte-align cy_dat array at element[%d]\n",sz_vd,l);
					break;
				}
			}
			ASSERT(HERE, l < 4, "Failed to align cy_dat array!");
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix36_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix36_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 192 16-byte slots of sc_arr for r-and-s temporaries, next 7 for the nontrivial complex 16th roots,
	next 36 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the   table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
		tmp = sc_ptr;			tm2 = tmp + 0x48;
		r00	= tmp + 0x00;		s1p00r = tm2 + 0x00;
		r02	= tmp + 0x02;		s1p01r = tm2 + 0x02;
		r04	= tmp + 0x04;		s1p02r = tm2 + 0x04;
		r06	= tmp + 0x06;		s1p03r = tm2 + 0x06;
		r08	= tmp + 0x08;		s1p04r = tm2 + 0x08;
		r0a	= tmp + 0x0a;		s1p05r = tm2 + 0x0a;
		r0c	= tmp + 0x0c;		s1p06r = tm2 + 0x0c;
		r0e	= tmp + 0x0e;		s1p07r = tm2 + 0x0e;
		r0g	= tmp + 0x10;		s1p08r = tm2 + 0x10;
		r10	= tmp + 0x12;		s1p09r = tm2 + 0x12;
		r12	= tmp + 0x14;		s1p10r = tm2 + 0x14;
		r14	= tmp + 0x16;		s1p11r = tm2 + 0x16;
		r16	= tmp + 0x18;		s1p12r = tm2 + 0x18;
		r18	= tmp + 0x1a;		s1p13r = tm2 + 0x1a;
		r1a	= tmp + 0x1c;		s1p14r = tm2 + 0x1c;
		r1c	= tmp + 0x1e;		s1p15r = tm2 + 0x1e;
		r1e	= tmp + 0x20;		s1p16r = tm2 + 0x20;
		r1g	= tmp + 0x22;		s1p17r = tm2 + 0x22;
		r20	= tmp + 0x24;		s1p18r = tm2 + 0x24;
		r22	= tmp + 0x26;		s1p19r = tm2 + 0x26;
		r24	= tmp + 0x28;		s1p20r = tm2 + 0x28;
		r26	= tmp + 0x2a;		s1p21r = tm2 + 0x2a;
		r28	= tmp + 0x2c;		s1p22r = tm2 + 0x2c;
		r2a	= tmp + 0x2e;		s1p23r = tm2 + 0x2e;
		r2c	= tmp + 0x30;		s1p24r = tm2 + 0x30;
		r2e	= tmp + 0x32;		s1p25r = tm2 + 0x32;
		r2g	= tmp + 0x34;		s1p26r = tm2 + 0x34;
		r30	= tmp + 0x36;		s1p27r = tm2 + 0x36;
		r32	= tmp + 0x38;		s1p28r = tm2 + 0x38;
		r34	= tmp + 0x3a;		s1p29r = tm2 + 0x3a;
		r36	= tmp + 0x3c;		s1p30r = tm2 + 0x3c;
		r38	= tmp + 0x3e;		s1p31r = tm2 + 0x3e;
		r3a	= tmp + 0x40;		s1p32r = tm2 + 0x40;
		r3c	= tmp + 0x42;		s1p33r = tm2 + 0x42;
		r3e	= tmp + 0x44;		s1p34r = tm2 + 0x44;
		r3g	= tmp + 0x46;		s1p35r = tm2 + 0x46;
		tmp	+= 0x92;	// Extra 2 slots here for two,one below - added those late, too lazy to rejigger all the existing offsets following
		two    = tmp - 2;	// AVX+ versions of Radix-32 DFT macros assume consts 2.0,1.0,sqrt2,isrt2 laid out thusly
		one    = tmp - 1;
		cc1    = tmp + 0;
		ss1    = tmp + 1;
		cc2    = tmp + 2;
		ss2    = tmp + 3;
		cc3m1  = tmp + 4;
		ss3    = tmp + 5;
		cc4    = tmp + 6;
		ss4    = tmp + 7;
		tmp += 0x8;	// sc_ptr += 0x9a
	#ifdef USE_AVX
		cy = tmp;		tmp += 9;
	#else
		cy = tmp;		tmp += 18;
	#endif
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 2 = 0xa5 [avx] or 0xae [sse2]; This is where the value of half_arr_offset36 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */

		ASSERT(HERE, (radix36_creals_in_local_store << l2_sz_vd) >= ((long)half_arr - (long)r00) + (20 << l2_sz_vd), "radix36_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0 );		VEC_DBL_INIT(one, 1.0 );
		VEC_DBL_INIT(cc1  , c	);		VEC_DBL_INIT(ss1, s );
		VEC_DBL_INIT(cc2  , c2  );		VEC_DBL_INIT(ss2, s2);
		VEC_DBL_INIT(cc3m1, c3m1);		VEC_DBL_INIT(ss3, s3);
		VEC_DBL_INIT(cc4  , c4  );		VEC_DBL_INIT(ss4, s4);
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		VEC_DBL_INIT(sse2_rnd, crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)cy - (int)two;	// #bytes in above sincos block of data
		tmp = two;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}
		nbytes = sz_vd;	// sse2_rnd is a solo (in the SIMD-vector) datum
		tmp = sse2_rnd;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		/* SSE2 version of the one_half array - we have a 2-bit lookup, low bit is from the low word of the carry pair,
		high bit from the high, i.e. based on this lookup index [listed with LSB at right], we have:

			index	half_lo	half_hi
			00		1.0		1.0
			01		.50		1.0
			10		1.0		.50
			11		.50		.50

		The inverse-weights computation uses a similar table, but with all entries multiplied by .50:

			index2	half_lo	half_hi
			00		.50		.50
			01		.25		.50
			10		.50		.25
			11		.25		.25

		We do similarly for the base[] and baseinv[] table lookups - each of these get 4 further slots in half_arr.
		We also allocate a further 4 16-byte slots [uninitialized] for storage of the wtl,wtn,wtlp1,wtnm1 locals.
		*/
		tmp = half_arr;

	#ifdef USE_AVX
		/* Forward-weight multipliers: */
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		/* Forward-base[] multipliers: */
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;

		nbytes = 64 << l2_sz_vd;

	#elif defined(USE_SSE2)

		ctmp = (struct complex *)tmp;
		/* Forward-weight multipliers: */
		ctmp->re = 1.0;	ctmp->im = 1.0;	++ctmp;
		ctmp->re = .50;	ctmp->im = 1.0;	++ctmp;
		ctmp->re = 1.0;	ctmp->im = .50;	++ctmp;
		ctmp->re = .50;	ctmp->im = .50;	++ctmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		ctmp->re = .50;	ctmp->im = .50;	++ctmp;
		ctmp->re = .25;	ctmp->im = .50;	++ctmp;
		ctmp->re = .50;	ctmp->im = .25;	++ctmp;
		ctmp->re = .25;	ctmp->im = .25;	++ctmp;
		/* Forward-base[] multipliers: */
		ctmp->re = base   [0];	ctmp->im = base   [0];	++ctmp;
		ctmp->re = base   [1];	ctmp->im = base   [0];	++ctmp;
		ctmp->re = base   [0];	ctmp->im = base   [1];	++ctmp;
		ctmp->re = base   [1];	ctmp->im = base   [1];	++ctmp;
		/* Inverse-base[] multipliers: */
		ctmp->re = baseinv[0];	ctmp->im = baseinv[0];	++ctmp;
		ctmp->re = baseinv[1];	ctmp->im = baseinv[0];	++ctmp;
		ctmp->re = baseinv[0];	ctmp->im = baseinv[1];	++ctmp;
		ctmp->re = baseinv[1];	ctmp->im = baseinv[1];	++ctmp;

		nbytes = 16 << l2_sz_vd;

	#endif

		// Propagate the above consts to the remaining threads:
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		/* Floating-point sign mask used for FABS on packed doubles: */
		sign_mask = sm_ptr;
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sign_mask+i) = (uint64)0x7FFFFFFFFFFFFFFFull;
		}

		// Set up the SIMD-tupled-32-bit-int SSE constants used by the carry macros:
		sse_bw  = sm_ptr + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		tmp64 = (uint64)bw;
		tmp64 = tmp64 + (tmp64 << 32);
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sse_bw+i) = tmp64;
		}

		sse_sw  = sse_bw + RE_IM_STRIDE;
		tmp64 = (uint64)sw;
		tmp64 = tmp64 + (tmp64 << 32);
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sse_sw+i) = tmp64;
		}

		sse_n   = sse_sw + RE_IM_STRIDE;
		tmp64 = (uint64)n;
		tmp64 = tmp64 + (tmp64 << 32);
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sse_n +i) = tmp64;
		}

		nbytes = 4 << l2_sz_vd;

	#ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;
		nbytes += 64;;
	#endif

		// Propagate the above consts to the remaining threads:
		tmp = (vec_dbl *)sm_ptr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	// For large radices, array-access to bjmodn means only init base-ptr here:
	#ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	#else
		bjmodn = (int*)(sse_n   + RE_IM_STRIDE);
	#endif

	#endif	// USE_SSE2

		pini = NDIVR/CY_THREADS;
		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );

		poff[0] =   0; poff[1] = p04    ; poff[2] = p08; poff[3] = p04+p08;
		poff[4] = p16; poff[5] = p04+p16; poff[6] = p24; poff[7] = p04+p24;
		poff[8] = p32;

		if(_cy[0])	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;
			for(i = 0; i < RADIX; i++) {
				free((void *)_bjmodn[i]); _bjmodn[i] = 0x0;
				free((void *)    _cy[i]);     _cy[i] = 0x0;
			}
			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;
			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		j = CY_THREADS*sizeof(int);
		_i       	= (int *)malloc(j);	ptr_prod += (uint32)(_i== 0x0);
		for(i = 0; i < RADIX; i++) {
			_bjmodn[i]	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn[i]== 0x0);
		}
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		for(i = 0; i < RADIX; i++) {
			_cy[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy[i]== 0x0);
		}
		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/RADIX-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		jhi = NDIVR/CY_THREADS;

		for(j=0; j < jhi; j++)
		{
			_bjmodnini[1] -= sw; _bjmodnini[1] = _bjmodnini[1] + ( (-(int)((uint32)_bjmodnini[1] >> 31)) & n);
		}

		if(CY_THREADS > 1)
		{
			for(ithread = 2; ithread <= CY_THREADS; ithread++)
			{
				_bjmodnini[ithread] = _bjmodnini[ithread-1] + _bjmodnini[1] - n; _bjmodnini[ithread] = _bjmodnini[ithread] + ( (-(int)((uint32)_bjmodnini[ithread] >> 31)) & n);
			}
		}
		/* Check upper element against scalar value, as precomputed in single-thread mode: */
		bjmodnini=0;
		for(j=0; j < jhi*CY_THREADS; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
		#ifdef USE_SSE2
			tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].r00 + ((long)half_arr - (long)r00);
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r00      = (double *)base;
			tdat[ithread].half_arr = (double *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

	}	/* endif(first_entry) */

/*...The radix-36 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy[i][ithread] = 0;
		}
	}
	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy[0][      0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = 0.0;
	}

for(outer=0; outer <= 1; outer++)
{
	_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (_i[0] = 1).	*/

	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	/*
	Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
	then simply overwrite it with 1 prior to starting the k-loop.
	*/
	khi = n_div_nwt/CY_THREADS;
	j = _bjmodnini[CY_THREADS];
	// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn[0][ithread] = _bjmodnini[ithread];
		for(i = 1; i < RADIX; i++) {
			MOD_ADD32(_bjmodn[i-1][ithread], j, n, _bjmodn[i][ithread]);
		}
		_jstart[ithread] = ithread*NDIVR/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
	}

#if defined(USE_SSE2) && defined(USE_PTHREAD)

	tmp = max_err;	VEC_DBL_INIT(tmp, 0.0);
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, sz_vd);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

#endif	// USE_PTHREAD

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	}

#ifdef USE_PTHREAD
	/* Populate the thread-specific data structs - use the invariant terms as memchecks: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		tdat[ithread].iter = iter;
	// int data:
		ASSERT(HERE, tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");

		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(HERE, tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].nwt == nwt, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = _maxerr[ithread];
		tdat[ithread].scale = scale;

	// pointer data:
		ASSERT(HERE, tdat[ithread].arrdat == a, "thread-local memcheck fail!");			/* Main data array */
		ASSERT(HERE, tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].si  == si, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(HERE, tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	#endif
		/* init carries: */
		for(i = 0; i < RADIX; i++) {
			tdat[ithread].cy[i] = _cy[i][ithread];
		}
	}
#endif

#ifdef USE_PTHREAD

	// If also using main thread to do work units, that task-dispatch occurs after all the threadpool-task launches:
	for(ithread = 0; ithread < pool_work_units; ithread++)
	{
		task_control.data = (void*)(&tdat[ithread]);
		threadpool_add_task(tpool, &task_control, task_is_blocking);

#else

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		/***** DEC/HP CC doesn't properly copy init value of maxerr = 0 into threads,
		so need to set once again explicitly for each: *****/
		maxerr = 0.0;
	#ifdef USE_SSE2
	//	VEC_DBL_INIT(max_err, 0.0);	*** must do this in conjunction with thread-local-data-copy
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		col = _col[ithread];
		co2 = _co2[ithread];
		co3 = _co3[ithread];

		for(l = 0; l < RADIX; l++) {
			bjmodn[l] = _bjmodn[l][ithread];
		}
		/* init carries	*/
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			tmp->d0 = _cy[l  ][ithread];
			tmp->d1 = _cy[l+1][ithread];
			tmp->d2 = _cy[l+2][ithread];
			tmp->d3 = _cy[l+3][ithread];
		}
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			tmp->d0 = _cy[l  ][ithread];
			tmp->d1 = _cy[l+1][ithread];
		}
	#else
		for(l = 0; l < RADIX; l++) {
			cy[l] = _cy[l][ithread];
		}
	#endif

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix36_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
			_cy[l+2][ithread] = tmp->d2;
			_cy[l+3][ithread] = tmp->d3;
		}
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
		}
		maxerr = MAX(max_err->d0,max_err->d1);
	#else
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = cy[l];
		}
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #if 0//def OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(HERE, 0x0 == cy36_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix40_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = tdat[ithread].cy[l];
		}
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here:

	The cleanup carries from the end of each length-N/RADIX set of contiguous data into the begining of the next
	can all be neatly processed as follows:

	(1) Invert the forward DIF FFT of the first block of RADIX complex elements in A and unweight;
	(2) Propagate cleanup carries among the real and imaginary parts of the RADIX outputs of (1);
	(3) Reweight and perform a forward DIF FFT on the result of (2);
	(4) If any of the exit carries from (2) are nonzero, advance to the next RADIX elements and repeat (1-4).
	*/
	for(l = 0; l < RADIX; l++) {
		t[l].re = _cy[l][CY_THREADS - 1];
	}
	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = _cy[l][ithread-1];
		}
	}
	_cy[0][0] =+t[RADIX-1].re;	/* ...The wraparound carry is here: */
	for(l = 1; l < RADIX; l++) {
		_cy[l][0] = t[l-1].re;
	}

	full_pass = 0;
	scale = 1;
	j_jhi = 7;

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			// Generate padded version of j, since prepadding pini is thread-count unsafe:
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
			for(l = 0; l < RADIX>>2; l++) {
				jt = j1 + poff[l ];	// poff[] = p04,p08,...,p56
				a[jt    ] *= radix_inv;
				a[jt+p01] *= radix_inv;
				a[jt+p02] *= radix_inv;
				a[jt+p03] *= radix_inv;
			}
		}
	}
}	/* endfor(outer) */

	dtmp = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(l = 0; l < RADIX; l++) {
			dtmp += fabs(_cy[l][ithread]);
		}
		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}
	if(dtmp != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	return(0);
}

/****************/

int	radix36_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-36 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-36 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int NDIVR, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17
		,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35
		,i,j,jt,jp,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17
		,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	static double radix_inv, n2inv;
	double rt,it,re
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i
		,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23,cy24,cy25,cy26,cy27,cy28,cy29,cy30,cy31,cy32,cy33,cy34,cy35
		,temp,scale;

#if PFETCH
	double *addr, *addp;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/36;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/36 in radix36_ditN_cy_dif1_nochk.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)36));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;
		p28 = p27 + p01;
		p29 = p28 + p01;
		p30 = p29 + p01;
		p31 = p30 + p01;
		p32 = p31 + p01;
		p33 = p32 + p01;
		p34 = p33 + p01;
		p35 = p34 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p29 = p29 + ( (p29 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p31 = p31 + ( (p31 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p33 = p33 + ( (p33 >> DAT_BITS) << PAD_BITS );
		p34 = p34 + ( (p34 >> DAT_BITS) << PAD_BITS );
		p35 = p35 + ( (p35 >> DAT_BITS) << PAD_BITS );

		bjmodnini=0;
		for(j=0; j < NDIVR; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-36 final DIT pass is here.	*/

	cy00= 0;		/* init carries	*/
	cy01= 0;
	cy02= 0;
	cy03= 0;
	cy04= 0;
	cy05= 0;
	cy06= 0;
	cy07= 0;
	cy08= 0;
	cy09= 0;
	cy10= 0;
	cy11= 0;
	cy12= 0;
	cy13= 0;
	cy14= 0;
	cy15= 0;
	cy16= 0;
	cy17= 0;
	cy18= 0;
	cy19= 0;
	cy20= 0;
	cy21= 0;
	cy22= 0;
	cy23= 0;
	cy24= 0;
	cy25= 0;
	cy26= 0;
	cy27= 0;
	cy28= 0;
	cy29= 0;
	cy30= 0;
	cy31= 0;
	cy32= 0;
	cy33= 0;
	cy34= 0;
	cy35= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy00= -2;
	}
	else
	{
		ASSERT(HERE,0,"Radix-36 currently only supports LL test mode!");
	}

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	jstart = 0;
	jhi = jstart+nwt-1;
	khi = n_div_nwt;

for(outer=0; outer <= 1; outer++)
{
	i = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (i = 1).	*/

	bjmodn00= 0;
	bjmodn01= bjmodnini;
	bjmodn02= bjmodn01+bjmodnini-n; bjmodn02= bjmodn02+ ( (-(int)((uint32)bjmodn02>> 31)) & n);
	bjmodn03= bjmodn02+bjmodnini-n; bjmodn03= bjmodn03+ ( (-(int)((uint32)bjmodn03>> 31)) & n);
	bjmodn04= bjmodn03+bjmodnini-n; bjmodn04= bjmodn04+ ( (-(int)((uint32)bjmodn04>> 31)) & n);
	bjmodn05= bjmodn04+bjmodnini-n; bjmodn05= bjmodn05+ ( (-(int)((uint32)bjmodn05>> 31)) & n);
	bjmodn06= bjmodn05+bjmodnini-n; bjmodn06= bjmodn06+ ( (-(int)((uint32)bjmodn06>> 31)) & n);
	bjmodn07= bjmodn06+bjmodnini-n; bjmodn07= bjmodn07+ ( (-(int)((uint32)bjmodn07>> 31)) & n);
	bjmodn08= bjmodn07+bjmodnini-n; bjmodn08= bjmodn08+ ( (-(int)((uint32)bjmodn08>> 31)) & n);
	bjmodn09= bjmodn08+bjmodnini-n; bjmodn09= bjmodn09+ ( (-(int)((uint32)bjmodn09>> 31)) & n);
	bjmodn10= bjmodn09+bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);
	bjmodn12= bjmodn11+bjmodnini-n; bjmodn12= bjmodn12+ ( (-(int)((uint32)bjmodn12>> 31)) & n);
	bjmodn13= bjmodn12+bjmodnini-n; bjmodn13= bjmodn13+ ( (-(int)((uint32)bjmodn13>> 31)) & n);
	bjmodn14= bjmodn13+bjmodnini-n; bjmodn14= bjmodn14+ ( (-(int)((uint32)bjmodn14>> 31)) & n);
	bjmodn15= bjmodn14+bjmodnini-n; bjmodn15= bjmodn15+ ( (-(int)((uint32)bjmodn15>> 31)) & n);
	bjmodn16= bjmodn15+bjmodnini-n; bjmodn16= bjmodn16+ ( (-(int)((uint32)bjmodn16>> 31)) & n);
	bjmodn17= bjmodn16+bjmodnini-n; bjmodn17= bjmodn17+ ( (-(int)((uint32)bjmodn17>> 31)) & n);
	bjmodn18= bjmodn17+bjmodnini-n; bjmodn18= bjmodn18+ ( (-(int)((uint32)bjmodn18>> 31)) & n);
	bjmodn19= bjmodn18+bjmodnini-n; bjmodn19= bjmodn19+ ( (-(int)((uint32)bjmodn19>> 31)) & n);
	bjmodn20= bjmodn19+bjmodnini-n; bjmodn20= bjmodn20+ ( (-(int)((uint32)bjmodn20>> 31)) & n);
	bjmodn21= bjmodn20+bjmodnini-n; bjmodn21= bjmodn21+ ( (-(int)((uint32)bjmodn21>> 31)) & n);
	bjmodn22= bjmodn21+bjmodnini-n; bjmodn22= bjmodn22+ ( (-(int)((uint32)bjmodn22>> 31)) & n);
	bjmodn23= bjmodn22+bjmodnini-n; bjmodn23= bjmodn23+ ( (-(int)((uint32)bjmodn23>> 31)) & n);
	bjmodn24= bjmodn23+bjmodnini-n; bjmodn24= bjmodn24+ ( (-(int)((uint32)bjmodn24>> 31)) & n);
	bjmodn25= bjmodn24+bjmodnini-n; bjmodn25= bjmodn25+ ( (-(int)((uint32)bjmodn25>> 31)) & n);
	bjmodn26= bjmodn25+bjmodnini-n; bjmodn26= bjmodn26+ ( (-(int)((uint32)bjmodn26>> 31)) & n);
	bjmodn27= bjmodn26+bjmodnini-n; bjmodn27= bjmodn27+ ( (-(int)((uint32)bjmodn27>> 31)) & n);
	bjmodn28= bjmodn27+bjmodnini-n; bjmodn28= bjmodn28+ ( (-(int)((uint32)bjmodn28>> 31)) & n);
	bjmodn29= bjmodn28+bjmodnini-n; bjmodn29= bjmodn29+ ( (-(int)((uint32)bjmodn29>> 31)) & n);
	bjmodn30= bjmodn29+bjmodnini-n; bjmodn30= bjmodn30+ ( (-(int)((uint32)bjmodn30>> 31)) & n);
	bjmodn31= bjmodn30+bjmodnini-n; bjmodn31= bjmodn31+ ( (-(int)((uint32)bjmodn31>> 31)) & n);
	bjmodn32= bjmodn31+bjmodnini-n; bjmodn32= bjmodn32+ ( (-(int)((uint32)bjmodn32>> 31)) & n);
	bjmodn33= bjmodn32+bjmodnini-n; bjmodn33= bjmodn33+ ( (-(int)((uint32)bjmodn33>> 31)) & n);
	bjmodn34= bjmodn33+bjmodnini-n; bjmodn34= bjmodn34+ ( (-(int)((uint32)bjmodn34>> 31)) & n);
	bjmodn35= bjmodn34+bjmodnini-n; bjmodn35= bjmodn35+ ( (-(int)((uint32)bjmodn35>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+36;
	co3=co2-36;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

	for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
	{
		for(j=jstart; j<jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
		{
		#ifdef USE_SSE2
			j1 = (j & mask01) + br4[j&3];
			j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
		#else
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		#endif
			j2 = j1+RE_IM_STRIDE;
		/*
		!...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do a radix-36 DIT transform...
		*/
	/*
	Twiddleless version requires us to swap inputs as follows:
	indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35
		  -> 0,32,28,24,20,16,12, 8, 4,27,23,19,15,11, 7, 3,35,31,18,14,10, 6, 2,34,30,26,22, 9, 5, 1,33,29,25,21,17,13

	I.e. start out with first quartet of indices {0,9,18,27}, permute those according to
	  {0,9,18,27}*35%36 = {0,27,18,9}, then each is head of a length-9 list of indices with decrement 4 in the radix-9 DFTs.

	Remember, inputs to DIT are bit-reversed, so
	a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35] contain
	x[0,18, 9,27, 3,21,12,30, 6,24,15,33, 1,19,10,28, 4,22,13,31, 7,25,16,34, 2,20,11,29, 5,23,14,32, 8,26,17,35], which get swapped to
	x[0,18,27, 9,24, 6,15,33,12,30, 3,21,32,14,23, 5,20, 2,11,29, 8,26,35,17,28,10,19, 1,16,34, 7,25, 4,22,31,13], which means the a-indices get swapped as
	a[0, 1, 3, 2| 9, 8,10,11| 6, 7, 4, 5|31,30,29,28|25,24,26,27|32,33,35,34|15,14,13,12|22,23,20,21|16,17,19,18]. These are the 9 quartets going into the radix-4 DFTs.
	*/
	#if 0

	/*...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do 9 radix-4 transforms...*/
						 /*          inputs           */ /*                                      outputs                                      */
		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00,t01,t10,t11,t20,t21,t30,t31,rt,it);
		RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],t02,t03,t12,t13,t22,t23,t32,t33,rt,it);
		RADIX_04_DIT(a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],t04,t05,t14,t15,t24,t25,t34,t35,rt,it);
		RADIX_04_DIT(a[j1+p31],a[j2+p31],a[j1+p30],a[j2+p30],a[j1+p29],a[j2+p29],a[j1+p28],a[j2+p28],t06,t07,t16,t17,t26,t27,t36,t37,rt,it);
		RADIX_04_DIT(a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p26],a[j2+p26],a[j1+p27],a[j2+p27],t08,t09,t18,t19,t28,t29,t38,t39,rt,it);
		RADIX_04_DIT(a[j1+p32],a[j2+p32],a[j1+p33],a[j2+p33],a[j1+p35],a[j2+p35],a[j1+p34],a[j2+p34],t0a,t0b,t1a,t1b,t2a,t2b,t3a,t3b,rt,it);
		RADIX_04_DIT(a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],t0c,t0d,t1c,t1d,t2c,t2d,t3c,t3d,rt,it);
		RADIX_04_DIT(a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],t0e,t0f,t1e,t1f,t2e,t2f,t3e,t3f,rt,it);
		RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],t0g,t0h,t1g,t1h,t2g,t2h,t3g,t3h,rt,it);

	/*...and now do 4 radix-5 transforms...*/
					 /*                            inputs                                 */ /*                 outputs                   */
		RADIX_09_DIT(t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,a1p00r,a1p00i,a1p32r,a1p32i,a1p28r,a1p28i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,rt,it,re);
		RADIX_09_DIT(t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p35r,a1p35i,a1p31r,a1p31i,rt,it,re);
		RADIX_09_DIT(t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,a1p18r,a1p18i,a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p34r,a1p34i,a1p30r,a1p30i,a1p26r,a1p26i,a1p22r,a1p22i,rt,it,re);
		RADIX_09_DIT(t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p33r,a1p33i,a1p29r,a1p29i,a1p25r,a1p25i,a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,rt,it,re);

	#else

		/*...Block 1:	*/
			jt = j1;	jp = j2;
			t00=a[jt    ];	t01=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t10=t00-rt;	t00=t00+rt;
			t11=t01-it;	t01=t01+it;

			t20=a[jt+p03];	t21=a[jp+p03];
			rt =a[jt+p02];	it =a[jp+p02];
			t30=t20-rt;	t20=t20+rt;
			t31=t21-it;	t21=t21+it;

			rt =t20;	t20=t00-rt;	t00=t00+rt;
			it =t21;	t21=t01-it;	t01=t01+it;

			rt =t30;	t30=t10-t31;	t10=t10+t31;
					t31=t11+rt;	t11=t11-rt;

		/*...Block 2:	*/
			jt = j1+p08;	jp = j2 + p08;
			t02=a[jt+p01];	t03=a[jp+p01];
			rt =a[jt    ];	it =a[jp    ];
			t12=t02-rt;	t02=t02+rt;
			t13=t03-it;	t03=t03+it;

			t22=a[jt+p02];	t23=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t32=t22-rt;	t22=t22+rt;
			t33=t23-it;	t23=t23+it;

			rt =t22;	t22=t02-rt;	t02=t02+rt;
			it =t23;	t23=t03-it;	t03=t03+it;

			rt =t32;	t32=t12-t33;	t12=t12+t33;
					t33=t13+rt;	t13=t13-rt;

		/*...Block 3:	*/
			jt = j1+p04;	jp = j2 + p04;
			t04=a[jt+p02];	t05=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t14=t04-rt;	t04=t04+rt;
			t15=t05-it;	t05=t05+it;

			t24=a[jt    ];	t25=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t34=t24-rt;	t24=t24+rt;
			t35=t25-it;	t25=t25+it;

			rt =t24;	t24=t04-rt;	t04=t04+rt;
			it =t25;	t25=t05-it;	t05=t05+it;

			rt =t34;	t34=t14-t35;	t14=t14+t35;
					t35=t15+rt;	t15=t15-rt;

		/*...Block 4:	*/
			jt = j1+p28;	jp = j2 + p28;
			t06=a[jt+p03];	t07=a[jp+p03];
			rt =a[jt+p02];	it =a[jp+p02];
			t16=t06-rt;	t06=t06+rt;
			t17=t07-it;	t07=t07+it;

			t26=a[jt+p01];	t27=a[jp+p01];
			rt =a[jt    ];	it =a[jp    ];
			t36=t26-rt;	t26=t26+rt;
			t37=t27-it;	t27=t27+it;

			rt =t26;	t26=t06-rt;	t06=t06+rt;
			it =t27;	t27=t07-it;	t07=t07+it;

			rt =t36;	t36=t16-t37;	t16=t16+t37;
					t37=t17+rt;	t17=t17-rt;

		/*...Block 5:	*/
			jt = j1+p24;	jp = j2 + p24;
			t08=a[jt+p01];	t09=a[jp+p01];
			rt =a[jt    ];	it =a[jp    ];
			t18=t08-rt;	t08=t08+rt;
			t19=t09-it;	t09=t09+it;

			t28=a[jt+p02];	t29=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t38=t28-rt;	t28=t28+rt;
			t39=t29-it;	t29=t29+it;

			rt =t28;	t28=t08-rt;	t08=t08+rt;
			it =t29;	t29=t09-it;	t09=t09+it;

			rt =t38;	t38=t18-t39;	t18=t18+t39;
					t39=t19+rt;	t19=t19-rt;

		/*...Block 6:	*/
			jt = j1+p32;	jp = j2 + p32;
			t0a=a[jt    ];	t0b=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t1a=t0a-rt;	t0a=t0a+rt;
			t1b=t0b-it;	t0b=t0b+it;

			t2a=a[jt+p03];	t2b=a[jp+p03];
			rt =a[jt+p02];	it =a[jp+p02];
			t3a=t2a-rt;	t2a=t2a+rt;
			t3b=t2b-it;	t2b=t2b+it;

			rt =t2a;	t2a=t0a-rt;	t0a=t0a+rt;
			it =t2b;	t2b=t0b-it;	t0b=t0b+it;

			rt =t3a;	t3a=t1a-t3b;	t1a=t1a+t3b;
					t3b=t1b+rt;	t1b=t1b-rt;

		/*...Block 7:	*/
			jt = j1+p12;	jp = j2 + p12;
			t0c=a[jt+p03];	t0d=a[jp+p03];
			rt =a[jt+p02];	it =a[jp+p02];
			t1c=t0c-rt;	t0c=t0c+rt;
			t1d=t0d-it;	t0d=t0d+it;

			t2c=a[jt+p01];	t2d=a[jp+p01];
			rt =a[jt    ];	it =a[jp    ];
			t3c=t2c-rt;	t2c=t2c+rt;
			t3d=t2d-it;	t2d=t2d+it;

			rt =t2c;	t2c=t0c-rt;	t0c=t0c+rt;
			it =t2d;	t2d=t0d-it;	t0d=t0d+it;

			rt =t3c;	t3c=t1c-t3d;	t1c=t1c+t3d;
					t3d=t1d+rt;	t1d=t1d-rt;

		/*...Block 8:	*/
			jt = j1+p20;	jp = j2 + p20;
			t0e=a[jt+p02];	t0f=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t1e=t0e-rt;	t0e=t0e+rt;
			t1f=t0f-it;	t0f=t0f+it;

			t2e=a[jt    ];	t2f=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t3e=t2e-rt;	t2e=t2e+rt;
			t3f=t2f-it;	t2f=t2f+it;

			rt =t2e;	t2e=t0e-rt;	t0e=t0e+rt;
			it =t2f;	t2f=t0f-it;	t0f=t0f+it;

			rt =t3e;	t3e=t1e-t3f;	t1e=t1e+t3f;
					t3f=t1f+rt;	t1f=t1f-rt;

		/*...Block 9:	*/
			jt = j1+p16;	jp = j2 + p16;
			t0g=a[jt    ];	t0h=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t1g=t0g-rt;	t0g=t0g+rt;
			t1h=t0h-it;	t0h=t0h+it;

			t2g=a[jt+p03];	t2h=a[jp+p03];
			rt =a[jt+p02];	it =a[jp+p02];
			t3g=t2g-rt;	t2g=t2g+rt;
			t3h=t2h-it;	t2h=t2h+it;

			rt =t2g;	t2g=t0g-rt;	t0g=t0g+rt;
			it =t2h;	t2h=t0h-it;	t0h=t0h+it;

			rt =t3g;	t3g=t1g-t3h;	t1g=t1g+t3h;
					t3h=t1h+rt;	t1h=t1h-rt;

	/*       ...and now do four radix-9 transforms.	*/

		/*...First radix-9 transform:	*/
			rt  =t02;			it  =t03;
			t02 =rt+t04;		t03 =it+t05;
			t04 =rt-t04;		t05 =it-t05;
			t00 =t00+t02;		t01 =t01+t03;
			t02 =t00+c3m1*t02;	t03 =t01+c3m1*t03;
			rt  =s3*t04;		it  =s3*t05;
			t04 =t02-it;		t05 =t03+rt;
			t02 =t02+it;		t03 =t03-rt;

			rt  =t08;			it  =t09;
			t08 =rt+t0a;		t09 =it+t0b;
			t0a =rt-t0a;		t0b =it-t0b;
			t06 =t06+t08;		t07 =t07+t09;
			t08 =t06+c3m1*t08;	t09 =t07+c3m1*t09;
			rt  =s3*t0a;		it  =s3*t0b;
			t0a =t08-it;		t0b =t09+rt;
			t08 =t08+it;		t09 =t09-rt;

			rt  =t0e;			it  =t0f;
			t0e =rt+t0g;		t0f =it+t0h;
			t0g =rt-t0g;		t0h =it-t0h;
			t0c =t0c+t0e;		t0d =t0d+t0f;
			t0e =t0c+c3m1*t0e;	t0f =t0d+c3m1*t0f;
			rt  =s3*t0g;		it  =s3*t0h;
			t0g =t0e-it;		t0h =t0f+rt;
			t0e =t0e+it;		t0f =t0f-rt;
		/* Twiddles: */
			rt  =t06;			it  =t07;
			t06 =rt+t0c;		t07 =it+t0d;
			t0c =rt-t0c;		t0d =it-t0d;
			t00 =t00+t06;		t01 =t01+t07;
			a1p00r=t00;		a1p00i=t01;
			t06 =t00+c3m1*t06;	t07 =t01+c3m1*t07;
			rt  =s3*t0c;		it  =s3*t0d;
			a1p24r=t06+it;		a1p24i=t07-rt;
			a1p12r=t06-it;		a1p12i=t07+rt;

			rt  =t08*c +t09*s;	it  =t09*c -t08*s;
			re  =t0e*c2+t0f*s2;	t0f =t0f*c2-t0e*s2;	t0e=re;
			t08 =rt+t0e;		t09 =it+t0f;
			t0e =rt-t0e;		t0f =it-t0f;
			t02 =t02+t08;		t03 =t03+t09;
			a1p08r=t02;		a1p08i=t03;
			t08 =t02+c3m1*t08;	t09=t03+c3m1*t09;
			rt  =s3*t0e;		it  =s3*t0f;
			a1p32r=t08+it;		a1p32i=t09-rt;
			a1p20r=t08-it;		a1p20i=t09+rt;

			rt  =t0a*c2+t0b*s2;	it  =t0b*c2-t0a*s2;
			re  =t0g*c4+t0h*s4;	t0h =t0h*c4-t0g*s4;	t0g=re;
			t0a =rt+t0g;		t0b =it+t0h;
			t0g =rt-t0g;		t0h =it-t0h;
			t04 =t04+t0a;		t05 =t05+t0b;
			a1p16r=t04;		a1p16i=t05;
			t0a =t04+c3m1*t0a;	t0b =t05+c3m1*t0b;
			rt  =s3*t0g;		it  =s3*t0h;
			a1p04r=t0a+it;		a1p04i=t0b-rt;
			a1p28r=t0a-it;		a1p28i=t0b+rt;

		/*...Second radix-9 transform:	*/
			rt  =t12;			it  =t13;
			t12 =rt+t14;		t13 =it+t15;
			t14 =rt-t14;		t15 =it-t15;
			t10 =t10+t12;		t11 =t11+t13;
			t12 =t10+c3m1*t12;	t13 =t11+c3m1*t13;
			rt  =s3*t14;		it  =s3*t15;
			t14 =t12-it;		t15 =t13+rt;
			t12 =t12+it;		t13 =t13-rt;

			rt  =t18;			it  =t19;
			t18 =rt+t1a;		t19 =it+t1b;
			t1a =rt-t1a;		t1b =it-t1b;
			t16 =t16+t18;		t17 =t17+t19;
			t18 =t16+c3m1*t18;	t19 =t17+c3m1*t19;
			rt  =s3*t1a;		it  =s3*t1b;
			t1a =t18-it;		t1b =t19+rt;
			t18 =t18+it;		t19 =t19-rt;

			rt  =t1e;			it  =t1f;
			t1e =rt+t1g;		t1f =it+t1h;
			t1g =rt-t1g;		t1h =it-t1h;
			t1c =t1c+t1e;		t1d =t1d+t1f;
			t1e =t1c+c3m1*t1e;	t1f =t1d+c3m1*t1f;
			rt  =s3*t1g;		it  =s3*t1h;
			t1g =t1e-it;		t1h =t1f+rt;
			t1e =t1e+it;		t1f =t1f-rt;
		/* Twiddles: */
			rt  =t16;			it  =t17;
			t16 =rt+t1c;		t17 =it+t1d;
			t1c =rt-t1c;		t1d =it-t1d;
			t10 =t10+t16;		t11 =t11+t17;
			a1p27r=t10;		a1p27i=t11;
			t16 =t10+c3m1*t16;	t17 =t11+c3m1*t17;
			rt  =s3*t1c;		it  =s3*t1d;
			a1p15r=t16+it;		a1p15i=t17-rt;
			a1p03r=t16-it;		a1p03i=t17+rt;

			rt  =t18*c +t19*s;	it  =t19*c -t18*s;
			re  =t1e*c2+t1f*s2;	t1f =t1f*c2-t1e*s2;	t1e=re;
			t18 =rt+t1e;		t19 =it+t1f;
			t1e =rt-t1e;		t1f =it-t1f;
			t12 =t12+t18;		t13 =t13+t19;
			a1p35r=t12;		a1p35i=t13;
			t18 =t12+c3m1*t18;	t19=t13+c3m1*t19;
			rt  =s3*t1e;		it  =s3*t1f;
			a1p23r=t18+it;		a1p23i=t19-rt;
			a1p11r=t18-it;		a1p11i=t19+rt;

			rt  =t1a*c2+t1b*s2;	it  =t1b*c2-t1a*s2;
			re  =t1g*c4+t1h*s4;	t1h =t1h*c4-t1g*s4;	t1g=re;
			t1a =rt+t1g;		t1b =it+t1h;
			t1g =rt-t1g;		t1h =it-t1h;
			t14 =t14+t1a;		t15 =t15+t1b;
			a1p07r=t14;		a1p07i=t15;
			t1a =t14+c3m1*t1a;	t1b =t15+c3m1*t1b;
			rt  =s3*t1g;		it  =s3*t1h;
			a1p31r=t1a+it;		a1p31i=t1b-rt;
			a1p19r=t1a-it;		a1p19i=t1b+rt;

		/*...Third radix-9 transform:	*/
			rt  =t22;			it  =t23;
			t22 =rt+t24;		t23 =it+t25;
			t24 =rt-t24;		t25 =it-t25;
			t20 =t20+t22;		t21 =t21+t23;
			t22 =t20+c3m1*t22;	t23 =t21+c3m1*t23;
			rt  =s3*t24;		it  =s3*t25;
			t24 =t22-it;		t25 =t23+rt;
			t22 =t22+it;		t23 =t23-rt;

			rt  =t28;			it  =t29;
			t28 =rt+t2a;		t29 =it+t2b;
			t2a =rt-t2a;		t2b =it-t2b;
			t26 =t26+t28;		t27 =t27+t29;
			t28 =t26+c3m1*t28;	t29 =t27+c3m1*t29;
			rt  =s3*t2a;		it  =s3*t2b;
			t2a =t28-it;		t2b =t29+rt;
			t28 =t28+it;		t29 =t29-rt;

			rt  =t2e;			it  =t2f;
			t2e =rt+t2g;		t2f =it+t2h;
			t2g =rt-t2g;		t2h =it-t2h;
			t2c =t2c+t2e;		t2d =t2d+t2f;
			t2e =t2c+c3m1*t2e;	t2f =t2d+c3m1*t2f;
			rt  =s3*t2g;		it  =s3*t2h;
			t2g =t2e-it;		t2h =t2f+rt;
			t2e =t2e+it;		t2f =t2f-rt;
		/* Twiddles: */
			rt  =t26;			it  =t27;
			t26 =rt+t2c;		t27 =it+t2d;
			t2c =rt-t2c;		t2d =it-t2d;
			t20 =t20+t26;		t21 =t21+t27;
			a1p18r=t20;		a1p18i=t21;
			t26 =t20+c3m1*t26;	t27 =t21+c3m1*t27;
			rt  =s3*t2c;		it  =s3*t2d;
			a1p06r=t26+it;		a1p06i=t27-rt;
			a1p30r=t26-it;		a1p30i=t27+rt;

			rt  =t28*c +t29*s;	it  =t29*c -t28*s;
			re  =t2e*c2+t2f*s2;	t2f =t2f*c2-t2e*s2;	t2e=re;
			t28 =rt+t2e;		t29 =it+t2f;
			t2e =rt-t2e;		t2f =it-t2f;
			t22 =t22+t28;		t23 =t23+t29;
			a1p26r=t22;		a1p26i=t23;
			t28 =t22+c3m1*t28;	t29=t23+c3m1*t29;
			rt  =s3*t2e;		it  =s3*t2f;
			a1p14r=t28+it;		a1p14i=t29-rt;
			a1p02r=t28-it;		a1p02i=t29+rt;

			rt  =t2a*c2+t2b*s2;	it  =t2b*c2-t2a*s2;
			re  =t2g*c4+t2h*s4;	t2h =t2h*c4-t2g*s4;	t2g=re;
			t2a =rt+t2g;		t2b =it+t2h;
			t2g =rt-t2g;		t2h =it-t2h;
			t24 =t24+t2a;		t25 =t25+t2b;
			a1p34r=t24;		a1p34i=t25;
			t2a =t24+c3m1*t2a;	t2b =t25+c3m1*t2b;
			rt  =s3*t2g;		it  =s3*t2h;
			a1p22r=t2a+it;		a1p22i=t2b-rt;
			a1p10r=t2a-it;		a1p10i=t2b+rt;

		/*...Fourth radix-9 transform:	*/
			rt  =t32;			it  =t33;
			t32 =rt+t34;		t33 =it+t35;
			t34 =rt-t34;		t35 =it-t35;
			t30 =t30+t32;		t31 =t31+t33;
			t32 =t30+c3m1*t32;	t33 =t31+c3m1*t33;
			rt  =s3*t34;		it  =s3*t35;
			t34 =t32-it;		t35 =t33+rt;
			t32 =t32+it;		t33 =t33-rt;

			rt  =t38;			it  =t39;
			t38 =rt+t3a;		t39 =it+t3b;
			t3a =rt-t3a;		t3b =it-t3b;
			t36 =t36+t38;		t37 =t37+t39;
			t38 =t36+c3m1*t38;	t39 =t37+c3m1*t39;
			rt  =s3*t3a;		it  =s3*t3b;
			t3a =t38-it;		t3b =t39+rt;
			t38 =t38+it;		t39 =t39-rt;

			rt  =t3e;			it  =t3f;
			t3e =rt+t3g;		t3f =it+t3h;
			t3g =rt-t3g;		t3h =it-t3h;
			t3c =t3c+t3e;		t3d =t3d+t3f;
			t3e =t3c+c3m1*t3e;	t3f =t3d+c3m1*t3f;
			rt  =s3*t3g;		it  =s3*t3h;
			t3g =t3e-it;		t3h =t3f+rt;
			t3e =t3e+it;		t3f =t3f-rt;
		/* Twiddles: */
			rt  =t36;			it  =t37;
			t36 =rt+t3c;		t37 =it+t3d;
			t3c =rt-t3c;		t3d =it-t3d;
			t30 =t30+t36;		t31 =t31+t37;
			a1p09r=t30;		a1p09i=t31;
			t36 =t30+c3m1*t36;	t37 =t31+c3m1*t37;
			rt  =s3*t3c;		it  =s3*t3d;
			a1p33r=t36+it;		a1p33i=t37-rt;
			a1p21r=t36-it;		a1p21i=t37+rt;

			rt  =t38*c +t39*s;	it  =t39*c -t38*s;
			re  =t3e*c2+t3f*s2;	t3f =t3f*c2-t3e*s2;	t3e=re;
			t38 =rt+t3e;		t39 =it+t3f;
			t3e =rt-t3e;		t3f =it-t3f;
			t32 =t32+t38;		t33 =t33+t39;
			a1p17r=t32;		a1p17i=t33;
			t38 =t32+c3m1*t38;	t39=t33+c3m1*t39;
			rt  =s3*t3e;		it  =s3*t3f;
			a1p05r=t38+it;		a1p05i=t39-rt;
			a1p29r=t38-it;		a1p29i=t39+rt;

			rt  =t3a*c2+t3b*s2;	it  =t3b*c2-t3a*s2;
			re  =t3g*c4+t3h*s4;	t3h =t3h*c4-t3g*s4;	t3g=re;
			t3a =rt+t3g;		t3b =it+t3h;
			t3g =rt-t3g;		t3h =it-t3h;
			t34 =t34+t3a;		t35 =t35+t3b;
			a1p25r=t34;		a1p25i=t35;
			t3a =t34+c3m1*t3a;	t3b =t35+c3m1*t3b;
			rt  =s3*t3g;		it  =s3*t3h;
			a1p13r=t3a+it;		a1p13i=t3b-rt;
			a1p01r=t3a-it;		a1p01i=t3b+rt;

	#endif

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 36 separate blocks of the A-array, we need 36 separate carries.	*/

			l= j & (nwt-1);
			n_minus_sil   = n-si[l  ];
			n_minus_silp1 = n-si[l+1];
			sinwt   = si[nwt-l  ];
			sinwtm1 = si[nwt-l-1];

			wtl     =wt0[    l  ];
			wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
			wtlp1   =wt0[    l+1];
			wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

			/*...set0 is slightly different from others:	*/
		   cmplx_carry_norm_nocheck0(a1p00r,a1p00i,cy00,bjmodn00   );
			cmplx_carry_norm_nocheck(a1p01r,a1p01i,cy01,bjmodn01,1 );
			cmplx_carry_norm_nocheck(a1p02r,a1p02i,cy02,bjmodn02,2 );
			cmplx_carry_norm_nocheck(a1p03r,a1p03i,cy03,bjmodn03,3 );
			cmplx_carry_norm_nocheck(a1p04r,a1p04i,cy04,bjmodn04,4 );
			cmplx_carry_norm_nocheck(a1p05r,a1p05i,cy05,bjmodn05,5 );
			cmplx_carry_norm_nocheck(a1p06r,a1p06i,cy06,bjmodn06,6 );
			cmplx_carry_norm_nocheck(a1p07r,a1p07i,cy07,bjmodn07,7 );
			cmplx_carry_norm_nocheck(a1p08r,a1p08i,cy08,bjmodn08,8 );
			cmplx_carry_norm_nocheck(a1p09r,a1p09i,cy09,bjmodn09,9 );
			cmplx_carry_norm_nocheck(a1p10r,a1p10i,cy10,bjmodn10,10);
			cmplx_carry_norm_nocheck(a1p11r,a1p11i,cy11,bjmodn11,11);
			cmplx_carry_norm_nocheck(a1p12r,a1p12i,cy12,bjmodn12,12);
			cmplx_carry_norm_nocheck(a1p13r,a1p13i,cy13,bjmodn13,13);
			cmplx_carry_norm_nocheck(a1p14r,a1p14i,cy14,bjmodn14,14);
			cmplx_carry_norm_nocheck(a1p15r,a1p15i,cy15,bjmodn15,15);
			cmplx_carry_norm_nocheck(a1p16r,a1p16i,cy16,bjmodn16,16);
			cmplx_carry_norm_nocheck(a1p17r,a1p17i,cy17,bjmodn17,17);
			cmplx_carry_norm_nocheck(a1p18r,a1p18i,cy18,bjmodn18,18);
			cmplx_carry_norm_nocheck(a1p19r,a1p19i,cy19,bjmodn19,19);
			cmplx_carry_norm_nocheck(a1p20r,a1p20i,cy20,bjmodn20,20);
			cmplx_carry_norm_nocheck(a1p21r,a1p21i,cy21,bjmodn21,21);
			cmplx_carry_norm_nocheck(a1p22r,a1p22i,cy22,bjmodn22,22);
			cmplx_carry_norm_nocheck(a1p23r,a1p23i,cy23,bjmodn23,23);
			cmplx_carry_norm_nocheck(a1p24r,a1p24i,cy24,bjmodn24,24);
			cmplx_carry_norm_nocheck(a1p25r,a1p25i,cy25,bjmodn25,25);
			cmplx_carry_norm_nocheck(a1p26r,a1p26i,cy26,bjmodn26,26);
			cmplx_carry_norm_nocheck(a1p27r,a1p27i,cy27,bjmodn27,27);
			cmplx_carry_norm_nocheck(a1p28r,a1p28i,cy28,bjmodn28,28);
			cmplx_carry_norm_nocheck(a1p29r,a1p29i,cy29,bjmodn29,29);
			cmplx_carry_norm_nocheck(a1p30r,a1p30i,cy30,bjmodn30,30);
			cmplx_carry_norm_nocheck(a1p31r,a1p31i,cy31,bjmodn31,31);
			cmplx_carry_norm_nocheck(a1p32r,a1p32i,cy32,bjmodn32,32);
			cmplx_carry_norm_nocheck(a1p33r,a1p33i,cy33,bjmodn33,33);
			cmplx_carry_norm_nocheck(a1p34r,a1p34i,cy34,bjmodn34,34);
			cmplx_carry_norm_nocheck(a1p35r,a1p35i,cy35,bjmodn35,35);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		/*...The radix-36 DIF pass is here:	*/
	#if 0
		/*...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do 4 radix-9 transforms...*/
							 /*                                                inputs                                                   */ /*                 outputs                   */
		RADIX_09_DIF(a1p00r,a1p00i,a1p32r,a1p32i,a1p28r,a1p28i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,rt,it,re);
		RADIX_09_DIF(a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p35r,a1p35i,a1p31r,a1p31i,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,rt,it,re);
		RADIX_09_DIF(a1p18r,a1p18i,a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p34r,a1p34i,a1p30r,a1p30i,a1p26r,a1p26i,a1p22r,a1p22i,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,rt,it,re);
		RADIX_09_DIF(a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p33r,a1p33i,a1p29r,a1p29i,a1p25r,a1p25i,a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,rt,it,re);

		/*...and now do 9 radix-4 transforms...*/
					/*          inputs           */ /*                                      outputs                                      */
		RADIX_04_DIF(t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it);
		RADIX_04_DIF(t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p32],a[j2+p32],a[j1+p33],a[j2+p33],a[j1+p35],a[j2+p35],a[j1+p34],a[j2+p34],rt,it);
		RADIX_04_DIF(t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],rt,it);
		RADIX_04_DIF(t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it);
		RADIX_04_DIF(t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p31],a[j2+p31],a[j1+p30],a[j2+p30],a[j1+p29],a[j2+p29],a[j1+p28],a[j2+p28],rt,it);
		RADIX_04_DIF(t0a,t0b,t1a,t1b,t2a,t2b,t3a,t3b,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it);
		RADIX_04_DIF(t0c,t0d,t1c,t1d,t2c,t2d,t3c,t3d,a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],rt,it);
		RADIX_04_DIF(t0e,t0f,t1e,t1f,t2e,t2f,t3e,t3f,a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p26],a[j2+p26],a[j1+p27],a[j2+p27],rt,it);
		RADIX_04_DIF(t0g,t0h,t1g,t1h,t2g,t2h,t3g,t3h,a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],rt,it);

	#else

		#if PFETCH
		addr = &a[j1];
		prefetch_p_doubles(addr);
		#endif
	/*...First radix-9 transform:	*/
		/* permuted version of jt + p00,p12,p24 */
			t00 =a1p00r;			t01 =a1p00i;
			t02 =a1p24r+a1p12r;		t03 =a1p24i+a1p12i;
			t04 =a1p24r-a1p12r;		t05 =a1p24i-a1p12i;
			t00 =t00+t02;			t01 =t01+t03;
			t02 =t00+c3m1*t02;		t03 =t01+c3m1*t03;
			rt  =s3*t04;			it  =s3*t05;
			t04 =t02+it;			t05 =t03-rt;
			t02 =t02-it;			t03 =t03+rt;
		#if PFETCH
		addp = addr+p01;
		prefetch_p_doubles(addp);
		#endif
		/* permuted version of jt + p04,p16,p28 */
			t06 =a1p32r;			t07 =a1p32i;
			t08 =a1p20r+a1p08r;		t09 =a1p20i+a1p08i;
			t0a =a1p20r-a1p08r;		t0b =a1p20i-a1p08i;
			t06 =t06+t08;			t07 =t07+t09;
			t08 =t06+c3m1*t08;		t09 =t07+c3m1*t09;
			rt  =s3*t0a;			it  =s3*t0b;
			t0a =t08+it;			t0b =t09-rt;
			t08 =t08-it;			t09 =t09+rt;
		#if PFETCH
		addp = addr+p02;
		prefetch_p_doubles(addp);
		#endif
		/* permuted version of jt + p08,p20,p32 */
			t0c =a1p28r;			t0d =a1p28i;
			t0e =a1p16r+a1p04r;		t0f =a1p16i+a1p04i;
			t0g =a1p16r-a1p04r;		t0h =a1p16i-a1p04i;
			t0c =t0c+t0e;			t0d =t0d+t0f;
			t0e =t0c+c3m1*t0e;		t0f =t0d+c3m1*t0f;
			rt  =s3*t0g;			it  =s3*t0h;
			t0g =t0e+it;			t0h =t0f-rt;
			t0e =t0e-it;			t0f =t0f+rt;
		#if PFETCH
		addp = addr+p03;
		prefetch_p_doubles(addp);
		#endif
		/* Twiddles: */
			rt  =t06;			it  =t07;
			t06 =rt+t0c;			t07 =it+t0d;
			t0c =rt-t0c;			t0d =it-t0d;
			t00 =t00+t06;			t01 =t01+t07;
			t06 =t00+c3m1*t06;		t07 =t01+c3m1*t07;
			rt  =s3*t0c;			it  =s3*t0d;
			t0c =t06+it;			t0d =t07-rt;
			t06 =t06-it;			t07 =t07+rt;
		#if PFETCH
		addp = addr+p04;
		prefetch_p_doubles(addp);
		#endif
			rt  =t08*c -t09*s;		it  =t08*s +t09*c;
			re  =t0e*c2-t0f*s2;		t0f =t0e*s2+t0f*c2;	t0e=re;
			t08 =rt+t0e;			t09 =it+t0f;
			t0e =rt-t0e;			t0f =it-t0f;
			t02 =t02+t08;			t03 =t03+t09;
			t08 =t02+c3m1*t08;		t09 =t03+c3m1*t09;
			rt  =s3*t0e;			it  =s3*t0f;
			t0e =t08+it;			t0f =t09-rt;
			t08 =t08-it;			t09 =t09+rt;
		#if PFETCH
		addp = addr+p05;
		prefetch_p_doubles(addp);
		#endif
			rt  =t0a*c2-t0b*s2;		it  =t0a*s2+t0b*c2;
			re  =t0g*c4-t0h*s4;		t0h =t0g*s4+t0h*c4;	t0g=re;
			t0a =rt+t0g;			t0b =it+t0h;
			t0g =rt-t0g;			t0h =it-t0h;
			t04 =t04+t0a;			t05 =t05+t0b;
			t0a =t04+c3m1*t0a;		t0b =t05+c3m1*t0b;
			rt  =s3*t0g;			it  =s3*t0h;
			t0g =t0a+it;			t0h =t0b-rt;
			t0a =t0a-it;			t0b =t0b+rt;
		#if PFETCH
		addp = addr+p06;
		prefetch_p_doubles(addp);
		#endif

	/*...Second radix-9 transform:	*/
		/* permuted version of jt + p01,p13,p25 */
			t10 =a1p27r;			t11 =a1p27i;
			t12 =a1p15r+a1p03r;		t13 =a1p15i+a1p03i;
			t14 =a1p15r-a1p03r;		t15 =a1p15i-a1p03i;
			t10 =t10+t12;			t11 =t11+t13;
			t12 =t10+c3m1*t12;		t13 =t11+c3m1*t13;
			rt  =s3*t14;			it  =s3*t15;
			t14 =t12+it;			t15 =t13-rt;
			t12 =t12-it;			t13 =t13+rt;
		#if PFETCH
		addp = addr+p07;
		prefetch_p_doubles(addp);
		#endif
		/* permuted version of jt + p05,p17,p29 */
			t16 =a1p23r;			t17 =a1p23i;
			t18 =a1p11r+a1p35r;		t19 =a1p11i+a1p35i;
			t1a =a1p11r-a1p35r;		t1b =a1p11i-a1p35i;
			t16 =t16+t18;			t17 =t17+t19;
			t18 =t16+c3m1*t18;		t19 =t17+c3m1*t19;
			rt  =s3*t1a;			it  =s3*t1b;
			t1a =t18+it;			t1b =t19-rt;
			t18 =t18-it;			t19 =t19+rt;
		#if PFETCH
		addp = addr+p08;
		prefetch_p_doubles(addp);
		#endif
		/* permuted version of jt + p09,p21,p33 */
			t1c =a1p19r;			t1d =a1p19i;
			t1e =a1p07r+a1p31r;		t1f =a1p07i+a1p31i;
			t1g =a1p07r-a1p31r;		t1h =a1p07i-a1p31i;
			t1c =t1c+t1e;			t1d =t1d+t1f;
			t1e =t1c+c3m1*t1e;		t1f =t1d+c3m1*t1f;
			rt  =s3*t1g;			it  =s3*t1h;
			t1g =t1e+it;			t1h =t1f-rt;
			t1e =t1e-it;			t1f =t1f+rt;
		#if PFETCH
		addp = addr+p09;
		prefetch_p_doubles(addp);
		#endif
		/* Twiddles: */
			rt  =t16;			it  =t17;
			t16 =rt+t1c;			t17 =it+t1d;
			t1c =rt-t1c;			t1d =it-t1d;
			t10 =t10+t16;			t11 =t11+t17;
			t16 =t10+c3m1*t16;		t17 =t11+c3m1*t17;
			rt  =s3*t1c;			it  =s3*t1d;
			t1c =t16+it;			t1d =t17-rt;
			t16 =t16-it;			t17 =t17+rt;
		#if PFETCH
		addp = addr+p10;
		prefetch_p_doubles(addp);
		#endif
			rt  =t18*c -t19*s;		it  =t18*s +t19*c;
			re  =t1e*c2-t1f*s2;		t1f =t1e*s2+t1f*c2;	t1e=re;
			t18 =rt+t1e;			t19 =it+t1f;
			t1e =rt-t1e;			t1f =it-t1f;
			t12 =t12+t18;			t13 =t13+t19;
			t18 =t12+c3m1*t18;		t19 =t13+c3m1*t19;
			rt  =s3*t1e;			it  =s3*t1f;
			t1e =t18+it;			t1f =t19-rt;
			t18 =t18-it;			t19 =t19+rt;
		#if PFETCH
		addp = addr+p11;
		prefetch_p_doubles(addp);
		#endif
			rt  =t1a*c2-t1b*s2;		it  =t1a*s2+t1b*c2;
			re  =t1g*c4-t1h*s4;		t1h =t1g*s4+t1h*c4;	t1g=re;
			t1a =rt+t1g;			t1b =it+t1h;
			t1g =rt-t1g;			t1h =it-t1h;
			t14 =t14+t1a;			t15 =t15+t1b;
			t1a =t14+c3m1*t1a;		t1b =t15+c3m1*t1b;
			rt  =s3*t1g;			it  =s3*t1h;
			t1g =t1a+it;			t1h =t1b-rt;
			t1a =t1a-it;			t1b =t1b+rt;
		#if PFETCH
		addp = addr+p12;
		prefetch_p_doubles(addp);
		#endif

	/*...Third radix-9 transform:	*/
		/* permuted version of jt + p02,p14,p26 */
			t20 =a1p18r;			t21 =a1p18i;
			t22 =a1p06r+a1p30r;		t23 =a1p06i+a1p30i;
			t24 =a1p06r-a1p30r;		t25 =a1p06i-a1p30i;
			t20 =t20+t22;			t21 =t21+t23;
			t22 =t20+c3m1*t22;		t23 =t21+c3m1*t23;
			rt  =s3*t24;			it  =s3*t25;
			t24 =t22+it;			t25 =t23-rt;
			t22 =t22-it;			t23 =t23+rt;
		#if PFETCH
		addp = addr+p13;
		prefetch_p_doubles(addp);
		#endif
		/* permuted version of jt + p06,p18,p30 */
			t26 =a1p14r;			t27 =a1p14i;
			t28 =a1p02r+a1p26r;		t29 =a1p02i+a1p26i;
			t2a =a1p02r-a1p26r;		t2b =a1p02i-a1p26i;
			t26 =t26+t28;			t27 =t27+t29;
			t28 =t26+c3m1*t28;		t29 =t27+c3m1*t29;
			rt  =s3*t2a;			it  =s3*t2b;
			t2a =t28+it;			t2b =t29-rt;
			t28 =t28-it;			t29 =t29+rt;
		#if PFETCH
		addp = addr+p14;
		prefetch_p_doubles(addp);
		#endif
		/* permuted version of jt + p10,p22,p34 */
			t2c =a1p10r;			t2d =a1p10i;
			t2e =a1p34r+a1p22r;		t2f =a1p34i+a1p22i;
			t2g =a1p34r-a1p22r;		t2h =a1p34i-a1p22i;
			t2c =t2c+t2e;			t2d =t2d+t2f;
			t2e =t2c+c3m1*t2e;		t2f =t2d+c3m1*t2f;
			rt  =s3*t2g;			it  =s3*t2h;
			t2g =t2e+it;			t2h =t2f-rt;
			t2e =t2e-it;			t2f =t2f+rt;
		#if PFETCH
		addp = addr+p15;
		prefetch_p_doubles(addp);
		#endif
		/* Twiddles: */
			rt  =t26;			it  =t27;
			t26 =rt+t2c;			t27 =it+t2d;
			t2c =rt-t2c;			t2d =it-t2d;
			t20 =t20+t26;			t21 =t21+t27;
			t26 =t20+c3m1*t26;		t27 =t21+c3m1*t27;
			rt  =s3*t2c;			it  =s3*t2d;
			t2c =t26+it;			t2d =t27-rt;
			t26 =t26-it;			t27 =t27+rt;
		#if PFETCH
		addp = addr+p16;
		prefetch_p_doubles(addp);
		#endif
			rt  =t28*c -t29*s;		it  =t28*s +t29*c;
			re  =t2e*c2-t2f*s2;		t2f =t2e*s2+t2f*c2;	t2e=re;
			t28 =rt+t2e;			t29 =it+t2f;
			t2e =rt-t2e;			t2f =it-t2f;
			t22 =t22+t28;			t23 =t23+t29;
			t28 =t22+c3m1*t28;		t29 =t23+c3m1*t29;
			rt  =s3*t2e;			it  =s3*t2f;
			t2e =t28+it;			t2f =t29-rt;
			t28 =t28-it;			t29 =t29+rt;
		#if PFETCH
		addp = addr+p17;
		prefetch_p_doubles(addp);
		#endif
			rt  =t2a*c2-t2b*s2;		it  =t2a*s2+t2b*c2;
			re  =t2g*c4-t2h*s4;		t2h =t2g*s4+t2h*c4;	t2g=re;
			t2a =rt+t2g;			t2b =it+t2h;
			t2g =rt-t2g;			t2h =it-t2h;
			t24 =t24+t2a;			t25 =t25+t2b;
			t2a =t24+c3m1*t2a;		t2b =t25+c3m1*t2b;
			rt  =s3*t2g;			it  =s3*t2h;
			t2g =t2a+it;			t2h =t2b-rt;
			t2a =t2a-it;			t2b =t2b+rt;
		#if PFETCH
		addp = addr+p18;
		prefetch_p_doubles(addp);
		#endif

	/*...Fourth radix-9 transform:	*/
		/* permuted version of jt + p03,p15,p27 */
			t30 =a1p09r;			t31 =a1p09i;
			t32 =a1p33r+a1p21r;		t33 =a1p33i+a1p21i;
			t34 =a1p33r-a1p21r;		t35 =a1p33i-a1p21i;
			t30 =t30+t32;			t31 =t31+t33;
			t32 =t30+c3m1*t32;		t33 =t31+c3m1*t33;
			rt  =s3*t34;			it  =s3*t35;
			t34 =t32+it;			t35 =t33-rt;
			t32 =t32-it;			t33 =t33+rt;
		#if PFETCH
		addp = addr+p19;
		prefetch_p_doubles(addp);
		#endif
		/* permuted version of jt + p07,p19,p31 */
			t36 =a1p05r;			t37 =a1p05i;
			t38 =a1p29r+a1p17r;		t39 =a1p29i+a1p17i;
			t3a =a1p29r-a1p17r;		t3b =a1p29i-a1p17i;
			t36 =t36+t38;			t37 =t37+t39;
			t38 =t36+c3m1*t38;		t39 =t37+c3m1*t39;
			rt  =s3*t3a;			it  =s3*t3b;
			t3a =t38+it;			t3b =t39-rt;
			t38 =t38-it;			t39 =t39+rt;
		#if PFETCH
		addp = addr+p20;
		prefetch_p_doubles(addp);
		#endif
		/* permuted version of jt + p11,p23,p35 */
			t3c =a1p01r;			t3d =a1p01i;
			t3e =a1p25r+a1p13r;		t3f =a1p25i+a1p13i;
			t3g =a1p25r-a1p13r;		t3h =a1p25i-a1p13i;
			t3c =t3c+t3e;			t3d =t3d+t3f;
			t3e =t3c+c3m1*t3e;		t3f =t3d+c3m1*t3f;
			rt  =s3*t3g;			it  =s3*t3h;
			t3g =t3e+it;			t3h =t3f-rt;
			t3e =t3e-it;			t3f =t3f+rt;
		#if PFETCH
		addp = addr+p21;
		prefetch_p_doubles(addp);
		#endif
		/* Twiddles: */
			rt  =t36;			it  =t37;
			t36 =rt+t3c;			t37 =it+t3d;
			t3c =rt-t3c;			t3d =it-t3d;
			t30 =t30+t36;			t31 =t31+t37;
			t36 =t30+c3m1*t36;		t37 =t31+c3m1*t37;
			rt  =s3*t3c;			it  =s3*t3d;
			t3c =t36+it;			t3d =t37-rt;
			t36 =t36-it;			t37 =t37+rt;
		#if PFETCH
		addp = addr+p22;
		prefetch_p_doubles(addp);
		#endif
			rt  =t38*c -t39*s;		it  =t38*s +t39*c;
			re  =t3e*c2-t3f*s2;		t3f =t3e*s2+t3f*c2;	t3e=re;
			t38 =rt+t3e;			t39 =it+t3f;
			t3e =rt-t3e;			t3f =it-t3f;
			t32 =t32+t38;			t33 =t33+t39;
			t38 =t32+c3m1*t38;		t39 =t33+c3m1*t39;
			rt  =s3*t3e;			it  =s3*t3f;
			t3e =t38+it;			t3f =t39-rt;
			t38 =t38-it;			t39 =t39+rt;
		#if PFETCH
		addp = addr+p23;
		prefetch_p_doubles(addp);
		#endif
			rt  =t3a*c2-t3b*s2;		it  =t3a*s2+t3b*c2;
			re  =t3g*c4-t3h*s4;		t3h =t3g*s4+t3h*c4;	t3g=re;
			t3a =rt+t3g;			t3b =it+t3h;
			t3g =rt-t3g;			t3h =it-t3h;
			t34 =t34+t3a;			t35 =t35+t3b;
			t3a =t34+c3m1*t3a;		t3b =t35+c3m1*t3b;
			rt  =s3*t3g;			it  =s3*t3h;
			t3g =t3a+it;			t3h =t3b-rt;
			t3a =t3a-it;			t3b =t3b+rt;
		#if PFETCH
		addp = addr+p24;
		prefetch_p_doubles(addp);
		#endif

	/*...and now do nine radix-4 transforms:	*/

	/*...Block 1: */
	jt = j1;	jp = j2;
			rt =t20;	t20=t00-rt;	t00=t00+rt;
			it =t21;	t21=t01-it;	t01=t01+it;

			rt =t30;	t30=t10-rt;	t10=t10+rt;
			it =t31;	t31=t11-it;	t11=t11+it;

			a[jt    ]=t00+t10;	a[jp    ]=t01+t11;
			a[jt+p01]=t00-t10;	a[jp+p01]=t01-t11;

			a[jt+p03]=t20-t31;	a[jp+p03]=t21+t30;	/* mpy by I is inlined here...	*/
			a[jt+p02]=t20+t31;	a[jp+p02]=t21-t30;
		#if PFETCH
		addp = addr+p25;
		prefetch_p_doubles(addp);
		#endif
	/*...Block 2: */
	jt = j1+p32;	jp = j2 + p32;
			rt =t22;	t22=t02-rt;	t02=t02+rt;
			it =t23;	t23=t03-it;	t03=t03+it;

			rt =t32;	t32=t12-rt;	t12=t12+rt;
			it =t33;	t33=t13-it;	t13=t13+it;

			a[jt    ]=t02+t12;	a[jp    ]=t03+t13;
			a[jt+p01]=t02-t12;	a[jp+p01]=t03-t13;

			a[jt+p03]=t22-t33;	a[jp+p03]=t23+t32;	/* mpy by I is inlined here...	*/
			a[jt+p02]=t22+t33;	a[jp+p02]=t23-t32;
		#if PFETCH
		addp = addr+p26;
		prefetch_p_doubles(addp);
		#endif
	/*...Block 3: */
	jt = j1+p20;	jp = j2 + p20;
			rt =t24;	t24=t04-rt;	t04=t04+rt;
			it =t25;	t25=t05-it;	t05=t05+it;

			rt =t34;	t34=t14-rt;	t14=t14+rt;
			it =t35;	t35=t15-it;	t15=t15+it;

			a[jt+p02]=t04+t14;	a[jp+p02]=t05+t15;
			a[jt+p03]=t04-t14;	a[jp+p03]=t05-t15;

			a[jt    ]=t24-t35;	a[jp    ]=t25+t34;	/* mpy by I is inlined here...	*/
			a[jt+p01]=t24+t35;	a[jp+p01]=t25-t34;
		#if PFETCH
		addp = addr+p27;
		prefetch_p_doubles(addp);
		#endif
	/*...Block 4: */
	jt = j1+p08;	jp = j2 + p08;
			rt =t26;	t26=t06-rt;	t06=t06+rt;
			it =t27;	t27=t07-it;	t07=t07+it;

			rt =t36;	t36=t16-rt;	t16=t16+rt;
			it =t37;	t37=t17-it;	t17=t17+it;

			a[jt+p01]=t06+t16;	a[jp+p01]=t07+t17;
			a[jt    ]=t06-t16;	a[jp    ]=t07-t17;

			a[jt+p02]=t26-t37;	a[jp+p02]=t27+t36;	/* mpy by I is inlined here...	*/
			a[jt+p03]=t26+t37;	a[jp+p03]=t27-t36;
		#if PFETCH
		addp = addr+p28;
		prefetch_p_doubles(addp);
		#endif
	/*...Block 5: */
	jt = j1+p28;	jp = j2 + p28;
			rt =t28;	t28=t08-rt;	t08=t08+rt;
			it =t29;	t29=t09-it;	t09=t09+it;

			rt =t38;	t38=t18-rt;	t18=t18+rt;
			it =t39;	t39=t19-it;	t19=t19+it;

			a[jt+p03]=t08+t18;	a[jp+p03]=t09+t19;
			a[jt+p02]=t08-t18;	a[jp+p02]=t09-t19;

			a[jt+p01]=t28-t39;	a[jp+p01]=t29+t38;	/* mpy by I is inlined here...	*/
			a[jt    ]=t28+t39;	a[jp    ]=t29-t38;
		#if PFETCH
		addp = addr+p29;
		prefetch_p_doubles(addp);
		#endif
	/*...Block 6: */
	jt = j1+p16;	jp = j2 + p16;
			rt =t2a;	t2a=t0a-rt;	t0a=t0a+rt;
			it =t2b;	t2b=t0b-it;	t0b=t0b+it;

			rt =t3a;	t3a=t1a-rt;	t1a=t1a+rt;
			it =t3b;	t3b=t1b-it;	t1b=t1b+it;

			a[jt    ]=t0a+t1a;	a[jp    ]=t0b+t1b;
			a[jt+p01]=t0a-t1a;	a[jp+p01]=t0b-t1b;

			a[jt+p03]=t2a-t3b;	a[jp+p03]=t2b+t3a;	/* mpy by I is inlined here...	*/
			a[jt+p02]=t2a+t3b;	a[jp+p02]=t2b-t3a;
		#if PFETCH
		addp = addr+p30;
		prefetch_p_doubles(addp);
		#endif
	/*...Block 7: */
	jt = j1+p04;	jp = j2 + p04;
			rt =t2c;	t2c=t0c-rt;	t0c=t0c+rt;
			it =t2d;	t2d=t0d-it;	t0d=t0d+it;

			rt =t3c;	t3c=t1c-rt;	t1c=t1c+rt;
			it =t3d;	t3d=t1d-it;	t1d=t1d+it;

			a[jt+p02]=t0c+t1c;	a[jp+p02]=t0d+t1d;
			a[jt+p03]=t0c-t1c;	a[jp+p03]=t0d-t1d;

			a[jt    ]=t2c-t3d;	a[jp    ]=t2d+t3c;	/* mpy by I is inlined here...	*/
			a[jt+p01]=t2c+t3d;	a[jp+p01]=t2d-t3c;
		#if PFETCH
		addp = addr+p31;
		prefetch_p_doubles(addp);
		#endif
	/*...Block 8: */
	jt = j1+p24;	jp = j2 + p24;
			rt =t2e;	t2e=t0e-rt;	t0e=t0e+rt;
			it =t2f;	t2f=t0f-it;	t0f=t0f+it;

			rt =t3e;	t3e=t1e-rt;	t1e=t1e+rt;
			it =t3f;	t3f=t1f-it;	t1f=t1f+it;

			a[jt+p01]=t0e+t1e;	a[jp+p01]=t0f+t1f;
			a[jt    ]=t0e-t1e;	a[jp    ]=t0f-t1f;

			a[jt+p02]=t2e-t3f;	a[jp+p02]=t2f+t3e;	/* mpy by I is inlined here...	*/
			a[jt+p03]=t2e+t3f;	a[jp+p03]=t2f-t3e;
		#if PFETCH
		addp = addr+p32;
		prefetch_p_doubles(addp);
		#endif
	/*...Block 9: */
	jt = j1+p12;	jp = j2 + p12;
			rt =t2g;	t2g=t0g-rt;	t0g=t0g+rt;
			it =t2h;	t2h=t0h-it;	t0h=t0h+it;

			rt =t3g;	t3g=t1g-rt;	t1g=t1g+rt;
			it =t3h;	t3h=t1h-it;	t1h=t1h+it;

			a[jt+p03]=t0g+t1g;	a[jp+p03]=t0h+t1h;
			a[jt+p02]=t0g-t1g;	a[jp+p02]=t0h-t1h;

			a[jt+p01]=t2g-t3h;	a[jp+p01]=t2h+t3g;	/* mpy by I is inlined here...	*/
			a[jt    ]=t2g+t3h;	a[jp    ]=t2h-t3g;
		#if PFETCH
		addp = addr+p33;
		prefetch_p_doubles(addp);
		#endif
		#if PFETCH
		addp = addr+p34;
		prefetch_p_doubles(addp);
		#endif
		#if PFETCH
		addp = addr+p35;
		prefetch_p_doubles(addp);
		#endif

	#endif	/* #if 0|1 */

			iroot += root_incr;		/* increment sincos index.	*/

		}

		jstart += nwt;
		jhi    += nwt;
		col += 36;
		co3 -= 36;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-36 forward DIF FFT of the first block of 36 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 36 outputs of (1);
!   (3) Reweight and perform a radix-36 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 36 elements and repeat (1-4).
*/
	t00 = cy35;
	cy35= cy34;
	cy34= cy33;
	cy33= cy32;
	cy32= cy31;
	cy31= cy30;
	cy30= cy29;
	cy29= cy28;
	cy28= cy27;
	cy27= cy26;
	cy26= cy25;
	cy25= cy24;
	cy24= cy23;
	cy23= cy22;
	cy22= cy21;
	cy21= cy20;
	cy20= cy19;
	cy19= cy18;
	cy18= cy17;
	cy17= cy16;
	cy16= cy15;
	cy15= cy14;
	cy14= cy13;
	cy13= cy12;
	cy12= cy11;
	cy11= cy10;
	cy10= cy09;
	cy09= cy08;
	cy08= cy07;
	cy07= cy06;
	cy06= cy05;
	cy05= cy04;
	cy04= cy03;
	cy03= cy02;
	cy02= cy01;
	cy01= cy00;
	cy00= t00 ;

	iroot = 0;
	root_incr = 0;
	scale = 1;

	jstart = 0;
	jhi = 7;
	khi = 1;

	for(j=0; j<=jhi; j++)
	{
		a[j    ] *= radix_inv;
		a[j+p01] *= radix_inv;
		a[j+p02] *= radix_inv;
		a[j+p03] *= radix_inv;
		a[j+p04] *= radix_inv;
		a[j+p05] *= radix_inv;
		a[j+p06] *= radix_inv;
		a[j+p07] *= radix_inv;
		a[j+p08] *= radix_inv;
		a[j+p09] *= radix_inv;
		a[j+p10] *= radix_inv;
		a[j+p11] *= radix_inv;
		a[j+p12] *= radix_inv;
		a[j+p13] *= radix_inv;
		a[j+p14] *= radix_inv;
		a[j+p15] *= radix_inv;
		a[j+p16] *= radix_inv;
		a[j+p17] *= radix_inv;
		a[j+p18] *= radix_inv;
		a[j+p19] *= radix_inv;
		a[j+p20] *= radix_inv;
		a[j+p21] *= radix_inv;
		a[j+p22] *= radix_inv;
		a[j+p23] *= radix_inv;
		a[j+p24] *= radix_inv;
		a[j+p25] *= radix_inv;
		a[j+p26] *= radix_inv;
		a[j+p27] *= radix_inv;
		a[j+p28] *= radix_inv;
		a[j+p29] *= radix_inv;
		a[j+p30] *= radix_inv;
		a[j+p31] *= radix_inv;
		a[j+p32] *= radix_inv;
		a[j+p33] *= radix_inv;
		a[j+p34] *= radix_inv;
		a[j+p35] *= radix_inv;
	}
}

	if(fabs(cy00)+fabs(cy01)+fabs(cy02)+fabs(cy03)+fabs(cy04)+fabs(cy05)+fabs(cy06)+fabs(cy07)+fabs(cy08)+fabs(cy09)+fabs(cy10)+fabs(cy11)+fabs(cy12)+fabs(cy13)+fabs(cy14)+fabs(cy15)+fabs(cy16)+fabs(cy17)+fabs(cy18)+fabs(cy19)+fabs(cy20)+fabs(cy21)+fabs(cy22)+fabs(cy23)+fabs(cy24)+fabs(cy25)+fabs(cy26)+fabs(cy27)+fabs(cy28)+fabs(cy29)+fabs(cy30)+fabs(cy31)+fabs(cy32)+fabs(cy33)+fabs(cy34)+fabs(cy35) != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix36_ditN_cy_dif1_nochk - input wordsize may be too small.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	return(0);
}

/***************/

void radix36_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-36 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,jt,jp,j1,j2;
	static int NDIVR,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32;
	static double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double rt,it,re;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/RADIX;

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-36 pass is here.	*/

	for(j=0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;
	/*
	Twiddleless version arranges 4 sets of radix-9 DFT inputs as follows: 0 in upper left corner, decrement 4 horizontally and 9 vertically:

		RADIX_09_DFT(00,32,28,24,20,16,12,08,04)
		RADIX_09_DFT(27,23,19,15,11,07,03,35,31)
		RADIX_09_DFT(18,14,10,06,02,34,30,26,22)
		RADIX_09_DFT(09,05,01,33,29,25,21,17,13)

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-4 DFT outputs.
	*/
		/*...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do 4 radix-9 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_09_DIF(a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, rt,it,re);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIF(a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, rt,it,re);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIF(a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, rt,it,re);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIF(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, rt,it,re);
		/*...and now do 9 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

		/* Totals: 4*radix09 + 9*radix04 = 4*(68 FADD, 40 FMUL)	+ 9*(16 FADD, 0 FMUL) = 416 FADD, 160 FMUL	*/
	}
}

/***************/

void radix36_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-36 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix36_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int j,jt,jp,j1,j2;
	static int NDIVR,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32;
	static double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double rt,it,re;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/RADIX;

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-36 pass is here.	*/

	for(j=0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;
	/*
	Twiddleless version uses same linear-index-vector-form permutation as in DIF:

		[0,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]
		00,32,28,24,20,16,12,08,04,27,23,19,15,11,07,03,35,31,18,14,10,06,02,34,30,26,22,09,05,01,33,29,25,21,17,13].	(*)

	Remember, inputs to DIT are bit-reversed, so

		a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35] contain
		x[0,18, 9,27, 3,21,12,30, 6,24,15,33, 1,19,10,28, 4,22,13,31, 7,25,16,34, 2,20,11,29, 5,23,14,32, 8,26,17,35], which get swapped [using the permutation (*) on the index *values*] to
		x[0,18,27, 9,24, 6,15,33,12,30, 3,21,32,14,23, 5,20, 2,11,29, 8,26,35,17,28,10,19, 1,16,34, 7,25, 4,22,31,13], which means the a-indices get swapped as
		a[0, 1, 3, 2| 9, 8,10,11| 6, 7, 4, 5|31,30,29,28|25,24,26,27|32,33,35,34|15,14,13,12|22,23,20,21|16,17,19,18]. These are the 9 quartets going into the radix-4 DFTs.

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF] to properly permute the radix-5 DFT outputs.
	*/
		/*...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do 9 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);
		/*...and now do 4 radix-9 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],rt,it,re);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],rt,it,re);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],rt,it,re);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],rt,it,re);

		/* Totals: 4*radix09 + 9*radix04 = 4*(68 FADD, 40 FMUL)	+ 9*(16 FADD, 0 FMUL) = 416 FADD, 160 FMUL	*/
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy36_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32;
		int poff[RADIX>>2];	// Store mults of p04 offset for loop control
		int j,j1,j2,jt,jp,k,l,ntmp;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	  #ifdef USE_AVX2
		// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
		vec_dbl *rad9_iptr[9], *rad9_optr[9];
	  #endif
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		int *itmp;	// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *two,*one, *cc1, *ss1, *cc2, *ss2, *cc3m1, *ss3, *cc4, *ss4, *max_err, *sse2_rnd, *half_arr
			,*r00,*r02,*r04,*r06,*r08,*r0a,*r0c,*r0e,*r0g
			,*r10,*r12,*r14,*r16,*r18,*r1a,*r1c,*r1e,*r1g
			,*r20,*r22,*r24,*r26,*r28,*r2a,*r2c,*r2e,*r2g
			,*r30,*r32,*r34,*r36,*r38,*r3a,*r3c,*r3e,*r3g
			,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r
			,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r
			,*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
		vec_dbl *tmp,*tm1,*tm2;	// Non-static utility ptrs
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
						s   =  0.64278760968653932631,	/* sin(2*pi/9) */
						c2  =  0.17364817766693034887,	/* cos(2*u) */
						s2  =  0.98480775301220805936,	/* sin(2*u) */
						c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
						s3  =  0.86602540378443864677,	/* sin(3*u) */
						c4  = -0.93969262078590838404,	/* cos(4*u) */
						s4  =  0.34202014332566873307;	/* sin(4*u) */
		double *base, *baseinv;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy = thread_arg->cy, rt,it,re, temp,frac;
		struct complex t[RADIX], *tptr;	// Bizarre: If we try using [RADIX] here, gcc gives "error: storage size of ‘r’ isn’t constant"
		int *itmp;	// Pointer into the bjmodn array

	#endif

	// int data:
		int iter = thread_arg->iter;
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int khi    = thread_arg->khi;
		int i      = thread_arg->i;	/* Pointer to the BASE and BASEINV arrays.	*/
		int jstart = thread_arg->jstart;
		int jhi    = thread_arg->jhi;
		int col = thread_arg->col;
		int co2 = thread_arg->co2;
		int co3 = thread_arg->co3;
		int sw  = thread_arg->sw;
		int nwt = thread_arg->nwt;

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;	int full_pass = scale < 0.5;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );

		poff[0] =   0; poff[1] = p04    ; poff[2] = p08; poff[3] = p04+p08;
		poff[4] = p16; poff[5] = p04+p16; poff[6] = p24; poff[7] = p04+p24;
		poff[8] = p32;

	#ifdef USE_SSE2
		tmp = thread_arg->r00;	tm2 = tmp + 0x48;
		r00	= tmp + 0x00;		s1p00r = tm2 + 0x00;
		r02	= tmp + 0x02;		s1p01r = tm2 + 0x02;
		r04	= tmp + 0x04;		s1p02r = tm2 + 0x04;
		r06	= tmp + 0x06;		s1p03r = tm2 + 0x06;
		r08	= tmp + 0x08;		s1p04r = tm2 + 0x08;
		r0a	= tmp + 0x0a;		s1p05r = tm2 + 0x0a;
		r0c	= tmp + 0x0c;		s1p06r = tm2 + 0x0c;
		r0e	= tmp + 0x0e;		s1p07r = tm2 + 0x0e;
		r0g	= tmp + 0x10;		s1p08r = tm2 + 0x10;
		r10	= tmp + 0x12;		s1p09r = tm2 + 0x12;
		r12	= tmp + 0x14;		s1p10r = tm2 + 0x14;
		r14	= tmp + 0x16;		s1p11r = tm2 + 0x16;
		r16	= tmp + 0x18;		s1p12r = tm2 + 0x18;
		r18	= tmp + 0x1a;		s1p13r = tm2 + 0x1a;
		r1a	= tmp + 0x1c;		s1p14r = tm2 + 0x1c;
		r1c	= tmp + 0x1e;		s1p15r = tm2 + 0x1e;
		r1e	= tmp + 0x20;		s1p16r = tm2 + 0x20;
		r1g	= tmp + 0x22;		s1p17r = tm2 + 0x22;
		r20	= tmp + 0x24;		s1p18r = tm2 + 0x24;
		r22	= tmp + 0x26;		s1p19r = tm2 + 0x26;
		r24	= tmp + 0x28;		s1p20r = tm2 + 0x28;
		r26	= tmp + 0x2a;		s1p21r = tm2 + 0x2a;
		r28	= tmp + 0x2c;		s1p22r = tm2 + 0x2c;
		r2a	= tmp + 0x2e;		s1p23r = tm2 + 0x2e;
		r2c	= tmp + 0x30;		s1p24r = tm2 + 0x30;
		r2e	= tmp + 0x32;		s1p25r = tm2 + 0x32;
		r2g	= tmp + 0x34;		s1p26r = tm2 + 0x34;
		r30	= tmp + 0x36;		s1p27r = tm2 + 0x36;
		r32	= tmp + 0x38;		s1p28r = tm2 + 0x38;
		r34	= tmp + 0x3a;		s1p29r = tm2 + 0x3a;
		r36	= tmp + 0x3c;		s1p30r = tm2 + 0x3c;
		r38	= tmp + 0x3e;		s1p31r = tm2 + 0x3e;
		r3a	= tmp + 0x40;		s1p32r = tm2 + 0x40;
		r3c	= tmp + 0x42;		s1p33r = tm2 + 0x42;
		r3e	= tmp + 0x44;		s1p34r = tm2 + 0x44;
		r3g	= tmp + 0x46;		s1p35r = tm2 + 0x46;
		tmp	+= 0x92;	// Extra 2 slots here for two,one below - added those late, too lazy to rejigger all the existing offsets following
		two    = tmp - 2;	// AVX+ versions of Radix-32 DFT macros assume consts 2.0,1.0,sqrt2,isrt2 laid out thusly
		one    = tmp - 1;
		cc1    = tmp + 0;
		ss1    = tmp + 1;
		cc2    = tmp + 2;
		ss2    = tmp + 3;
		cc3m1  = tmp + 4;
		ss3    = tmp + 5;
		cc4    = tmp + 6;
		ss4    = tmp + 7;
		tmp += 0x8;
	  #ifdef USE_AVX
		cy = tmp;		tmp += 9;
	  #else
		cy = tmp;		tmp += 18;
	  #endif
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */

		ASSERT(HERE, (r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
		tmp = half_arr;
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix36_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;

		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n + RE_IM_STRIDE);
	  #endif

	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->r00     ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		/* Init DWT-indices: */
		uint32 bjmodnini = thread_arg->bjmodnini;
		bjmodn[0] = thread_arg->bjmodn0;
		for(l = 1; l < RADIX; l++) {	// must use e.g. l for loop idx here as i is used for dwt indexing
			MOD_ADD32(bjmodn[l-1], bjmodnini, n, bjmodn[l]);
		}

		/* init carries	*/
		addr = thread_arg->cy;
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			tmp->d0 = *(addr+l  );
			tmp->d1 = *(addr+l+1);
			tmp->d2 = *(addr+l+2);
			tmp->d3 = *(addr+l+3);
		}
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			tmp->d0 = *(addr+l  );
			tmp->d1 = *(addr+l+1);
		}
	#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
		for(l = 0; l < RADIX; l++) {
			cy[l] = *(addr+l);
		}
	#endif

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix36_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		addr = thread_arg->cy;
	#ifdef USE_AVX
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			*(addr+l  ) = tmp->d0;
			*(addr+l+1) = tmp->d1;
			*(addr+l+2) = tmp->d2;
			*(addr+l+3) = tmp->d3;
		}
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			*(addr+l  ) = tmp->d0;
			*(addr+l+1) = tmp->d1;
		}
		maxerr = MAX(max_err->d0,max_err->d1);
	#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
		for(l = 0; l < RADIX; l++) {
			*(addr+l) = cy[l];
		}
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}

		return 0x0;
	}
#endif

#undef RADIX
#undef PFETCH_DIST
