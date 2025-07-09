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

#include "Mlucas.h"

#define RADIX 36	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

#ifndef COMPACT_OBJ	// Toggle for parametrized-loop-DFT compact-object code scheme
  #ifdef USE_AVX512
	#define COMPACT_OBJ	1	// AVX-512 builds [Intel KNL, specifically] only ones to show clear gain from this for radices < 64
  #endif
#endif

#ifndef PFETCH_DIST
  #ifdef USE_AVX512
	#define PFETCH_DIST	64	// Feb 2017: Test on KNL point to this as best
  #elif defined(USE_AVX)
	#define PFETCH_DIST	32	// This seems to work best on my Haswell, even though 64 bytes seems more logical in AVX mode
  #else
	#define PFETCH_DIST	32
  #endif
#endif

#ifdef USE_SSE2

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 96 [AVX]), add to (half_arr_offset36 + RADIX) to get value of radix36_creals_in_local_store.
  #ifdef USE_AVX512	// RADIX/8 = 4 fewer carry slots than AVX:
	const int half_arr_offset36 = 0xa1;	// + RADIX = 0xc5; Used for thread local-storage-integrity checking
	const int radix36_creals_in_local_store = 0x14c;	// += 132 (= 0x84) and round up to nearest multiple of 4
  #elif defined(USE_AVX)
	const int half_arr_offset36 = 0xa5;	// + RADIX = 0xc9; Used for thread local-storage-integrity checking
	const int radix36_creals_in_local_store = 0x150;	// += 132 (= 0x84) and round up to nearest multiple of 4
  #else
	const int half_arr_offset36 = 0xae;	// + RADIX = 0xd2; Used for thread local-storage-integrity checking
	const int radix36_creals_in_local_store = 0xfc;	// += 40 (= 0x28) and round up to nearest multiple of 4
  #endif

	#include "sse2_macro_gcc64.h"
	#include "radix09_sse_macro.h"

#endif	// SSE2

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data:
		int iter;
		int tid;
		int ndivr;
		int target_idx, target_set;	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
		double target_cy;

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
		double prp_mult;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		double *wts_mult, *inv_mult;
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
	  #ifdef USE_AVX512
		double cy_dat[RADIX+8] __attribute__ ((__aligned__(8)));
	  #else
		double cy_dat[RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	  #endif
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
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
  #endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
  #if !defined(MULTITHREAD) && (!defined(USE_SSE2) || (COMPACT_OBJ != 0))
	static uint32 p0123[4];
  #endif
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
  #endif
  #ifdef USE_AVX512
	const int jhi_wrap = 15;
  #else
	const int jhi_wrap =  7;
  #endif
	int NDIVR,i,j,j1,jt,jhi,full_pass,khi,l,outer;
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	int j2,jp;
  #endif
  #if !defined(MULTITHREAD) && (!defined(USE_SSE2) || defined(USE_AVX512))
	int ntmp;
  #endif
  #ifndef MULTITHREAD
	int jstart;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	int incr;
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = --|--|9 for avx512,avx,sse, respectively:
	const int incr_long = 9,incr_med = 3, incr_short = 3;
  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, so use "goes to 11" in LOACC mode via an incr_hiacc = 1:
  #if defined(USE_AVX512) || defined(USE_ARM_V8_SIMD)
	const int incr_hiacc = 1;
  #else
	const int incr_hiacc = 0;
  #endif
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(USE_SHORT_CY_CHAIN == 0)
		incr = incr_long;
	else if(USE_SHORT_CY_CHAIN == 1)
		incr = incr_med;
	else if(USE_SHORT_CY_CHAIN == 2)
		incr = incr_short;
	else
		incr = incr_hiacc;
  #endif

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32, nsave = 0;
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
  #ifndef MULTITHREAD
	double *addr;
  #endif
	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	struct complex *tptr;
  #endif
  #ifndef MULTITHREAD
	int *itmp;	// Pointer into the bjmodn array
  #endif
  #if !defined(MULTITHREAD) && defined(USE_AVX) && !defined(USE_AVX512)
	int *itm2;	// Pointer into the bjmodn array
  #endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
  #ifndef MULTITHREAD
	int col,co2,co3;
   #ifdef USE_AVX512
	double t0,t1,t2,t3;
	static struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #elif defined(USE_AVX)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
   #endif
   #ifdef USE_AVX2
	// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
	vec_dbl *rad9_iptr[9], *rad9_optr[9];
   #endif
  #endif

#ifdef USE_SSE2

  #if !defined(MULTITHREAD) && defined(USE_AVX512)
	int i0,i1,i2,i3;
  #endif
	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
  #endif

  #ifndef MULTITHREAD
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
  #endif
  #ifndef USE_AVX512
	const double crnd = 3.0*0x4000000*0x2000000;
  #endif
  #ifndef USE_AVX
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
  #endif
	vec_dbl *tmp,*tm2;	// Non-static utility ptrs
  #ifndef MULTITHREAD
	vec_dbl *tm1;	// Non-static utility ptrs
  #endif
	static vec_dbl *two,*one, *cc1, *ss1, *cc2, *ss2, *cc3m1, *ss3, *cc4, *ss4, *max_err, *sse2_rnd, *half_arr
		,*r00
	  #if !defined(MULTITHREAD)
	   #if !COMPACT_OBJ
	    ,*r02,*r04,*r06,*r08,*r0a,*r0c,*r0e,*r0g
		,*r10,*r12,*r14,*r16,*r18,*r1a,*r1c,*r1e,*r1g
		,*r20,*r22,*r24,*r26,*r28,*r2a,*r2c,*r2e,*r2g
		,*r30,*r32,*r34,*r36,*r38,*r3a,*r3c,*r3e,*r3g
	   #endif
		,*s1p00
	   #if !COMPACT_OBJ
		,*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p09,*s1p10,*s1p11,*s1p12,*s1p13,*s1p14,*s1p15,*s1p16,*s1p17
		,*s1p18,*s1p19,*s1p20,*s1p21,*s1p22,*s1p23,*s1p24,*s1p25,*s1p26,*s1p27,*s1p28,*s1p29,*s1p30,*s1p31,*s1p32,*s1p33,*s1p34,*s1p35
	   #endif
	  #endif
		,*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx

#endif

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
  #if 0//def OS_TYPE_MACOSX
	static int main_work_units = 0;
  #endif
	static int pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy36_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int bjmodn[RADIX];
	double rt,it,temp,frac,cy[RADIX], re;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodnini = 0x0,*_bjmodn[RADIX];
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_cy[RADIX];
	// In unthreaded-build mode, the _cy[] array having dimension == 4 (mod 8) causes the
	// AVX512 8-way carry-init loops to run beyond the end of the _cy[] array:
  #if defined(USE_AVX512) && !defined(USE_THREADS)
	WARN(HERE, "radix36_ditN_cy_dif1: AVX-512 support requires multithreaded build (-DUSE_THREADS); Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif

	if(!_jhi) {
		_cy[0] = 0x0;	// First of these used as an "already inited consts?" sentinel, must init = 0x0 at same time do so for non-array static ptrs
	}

  #ifdef USE_IMCI512
	WARN(HERE, "radix36_ditN_cy_dif1: No k10m/IMCI-512 support; Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(0, "Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

  #ifndef MULTITHREAD
	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;
  #endif
	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	if((TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;	ndivr_inv = (double)RADIX/n;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		fprintf(stderr,"ERROR: iter = %10d; NWT_BITS does not divide N/RADIX in %s.\n",iter,func);
		err = ERR_SKIP_RADIX_SET;
		return(err);
	}

	if(p != psave || n != nsave
	#ifdef USE_PTHREAD	// Oct 2021: cf. radix176_ditN_cy_dif1.c for why I added this
		|| (tdat != 0x0 && tdat[0].wt1 != wt1)
	#endif
	) {	/* Exponent or array length change triggers re-init */
		first_entry=TRUE;
		/* To-do: Support #thread change here! */
	}

/*...initialize things upon first entry: */

	if(first_entry)
	{
		psave = p;	nsave = n;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	  #ifdef USE_AVX	// AVX LOACC: Make CARRY_8_WAY the default for this mode
		i = 8;
	  #elif defined(USE_SSE2)	// AVX and SSE2 modes use 4-way carry macros
		i = 4;
	  #else	// Scalar-double mode:
		i = 1;
	  #endif

		// For n a power of 2 don't need to worry about 32-bit integer overflow in the sw*NDIVR term,
		// but for non-power-of-2 n we must cast-to-uint64 to avoid such overflows fubaring the result:
		struct qfloat qt,qn;
		qt = i64_to_q(i*(uint64)sw*NDIVR % n);
		qn = i64_to_q((int64) n);
		qt = qfdiv(qt, qn);		// x = (sw*NDIVR (mod n))/n
		qt = qfmul(qt, QLN2);	// x*ln(2)...
		qt = qfexp(qt);			// ...and get 2^x via exp[x*ln(2)].
		wts_mult[0] = qfdbl(qt);		// a = 2^(x/n), with x = sw
		inv_mult[0] = qfdbl(qfinv(qt));	// Double-based inversion (1.0 / wts_mult_a[0]) often gets LSB wrong
		ASSERT(fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		//curr have w, 2/w, separate-mul-by-1-or-0.5 gives [w,w/2] and [1/w,2/w] for i = 0,1, resp:
		wts_mult[1] = 0.5*wts_mult[0];
		inv_mult[1] = 2.0*inv_mult[0];
		ASSERT(fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the largest
		power of 2 that is <= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else
		{
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> i;
		}

		if(CY_THREADS > MAX_THREADS)
		{
		//	CY_THREADS = MAX_THREADS;
			fprintf(stderr,"WARN: CY_THREADS = %d exceeds number of cores = %d\n", CY_THREADS, MAX_THREADS);
		}
		if(!isPow2(CY_THREADS))		{ WARN(HERE, "CY_THREADS not a power of 2!", "", 1); return(ERR_ASSERT); }
		if(CY_THREADS > 1)
		{
			if(NDIVR    %CY_THREADS != 0) { WARN(HERE, "NDIVR    %CY_THREADS != 0 ... likely more threads than this leading radix can handle.", "", 1); return(ERR_ASSERT); }
			if(n_div_nwt%CY_THREADS != 0) { WARN(HERE, "n_div_nwt%CY_THREADS != 0 ... likely more threads than this leading radix can handle.", "", 1); return(ERR_ASSERT); }
		}

	  #ifdef USE_PTHREAD
		if(tdat == 0x0) {
			j = (uint32)sizeof(struct cy_thread_data_t);
			tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

			// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
			// so on that platform try to be clever and interleave main-thread and threadpool-work processing
			#if 0//def OS_TYPE_MACOSX

				if(CY_THREADS > 1) {
					main_work_units = CY_THREADS/2;
					pool_work_units = CY_THREADS - main_work_units;
					ASSERT(0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
					printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
				} else {
					main_work_units = 1;
					printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
				}

			#else

				pool_work_units = CY_THREADS;
				ASSERT(0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

			#endif

			fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
		}
	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].iter = iter;
		// int data:
			tdat[ithread].tid = ithread;
			tdat[ithread].ndivr = NDIVR;

			tdat[ithread].sw  = sw;
			tdat[ithread].nwt = nwt;

		// pointer data:
		//	tdat[ithread].arrdat = a;			/* Main data array */
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
			tdat[ithread].wts_mult = wts_mult;
			tdat[ithread].inv_mult = inv_mult;
			tdat[ithread].si  = si;

		// This array pointer must be set based on vec_dbl-sized alignment at runtime for each thread:
			for(l = 0; l < RE_IM_STRIDE; l++) {
				if( ((intptr_t)&tdat[ithread].cy_dat[l] & SZ_VDM1) == 0 ) {
					tdat[ithread].cy = &tdat[ithread].cy_dat[l];
				//	fprintf(stderr,"%d-byte-align cy_dat array at element[%d]\n",SZ_VD,l);
					break;
				}
			}
			ASSERT(l < RE_IM_STRIDE, "Failed to align cy_dat array!");
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(((intptr_t)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(((intptr_t)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix36_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix36_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 192 16-byte slots of sc_arr for r-and-s temporaries, next 7 for the nontrivial complex 16th roots,
	next 36 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the   table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
		tmp = sc_ptr;			tm2 = tmp + 0x48;
		r00	= tmp + 0x00;
	#ifndef MULTITHREAD
								s1p00 = tm2 + 0x00;
	  #if !COMPACT_OBJ
		r02	= tmp + 0x02;		s1p01 = tm2 + 0x02;
		r04	= tmp + 0x04;		s1p02 = tm2 + 0x04;
		r06	= tmp + 0x06;		s1p03 = tm2 + 0x06;
		r08	= tmp + 0x08;		s1p04 = tm2 + 0x08;
		r0a	= tmp + 0x0a;		s1p05 = tm2 + 0x0a;
		r0c	= tmp + 0x0c;		s1p06 = tm2 + 0x0c;
		r0e	= tmp + 0x0e;		s1p07 = tm2 + 0x0e;
		r0g	= tmp + 0x10;		s1p08 = tm2 + 0x10;
		r10	= tmp + 0x12;		s1p09 = tm2 + 0x12;
		r12	= tmp + 0x14;		s1p10 = tm2 + 0x14;
		r14	= tmp + 0x16;		s1p11 = tm2 + 0x16;
		r16	= tmp + 0x18;		s1p12 = tm2 + 0x18;
		r18	= tmp + 0x1a;		s1p13 = tm2 + 0x1a;
		r1a	= tmp + 0x1c;		s1p14 = tm2 + 0x1c;
		r1c	= tmp + 0x1e;		s1p15 = tm2 + 0x1e;
		r1e	= tmp + 0x20;		s1p16 = tm2 + 0x20;
		r1g	= tmp + 0x22;		s1p17 = tm2 + 0x22;
		r20	= tmp + 0x24;		s1p18 = tm2 + 0x24;
		r22	= tmp + 0x26;		s1p19 = tm2 + 0x26;
		r24	= tmp + 0x28;		s1p20 = tm2 + 0x28;
		r26	= tmp + 0x2a;		s1p21 = tm2 + 0x2a;
		r28	= tmp + 0x2c;		s1p22 = tm2 + 0x2c;
		r2a	= tmp + 0x2e;		s1p23 = tm2 + 0x2e;
		r2c	= tmp + 0x30;		s1p24 = tm2 + 0x30;
		r2e	= tmp + 0x32;		s1p25 = tm2 + 0x32;
		r2g	= tmp + 0x34;		s1p26 = tm2 + 0x34;
		r30	= tmp + 0x36;		s1p27 = tm2 + 0x36;
		r32	= tmp + 0x38;		s1p28 = tm2 + 0x38;
		r34	= tmp + 0x3a;		s1p29 = tm2 + 0x3a;
		r36	= tmp + 0x3c;		s1p30 = tm2 + 0x3c;
		r38	= tmp + 0x3e;		s1p31 = tm2 + 0x3e;
		r3a	= tmp + 0x40;		s1p32 = tm2 + 0x40;
		r3c	= tmp + 0x42;		s1p33 = tm2 + 0x42;
		r3e	= tmp + 0x44;		s1p34 = tm2 + 0x44;
		r3g	= tmp + 0x46;		s1p35 = tm2 + 0x46;
	  #endif
	#endif
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
	#ifdef USE_AVX512
		cy = tmp;		tmp += 5;		// RADIX/8 and round up
	#elif defined(USE_AVX)
		cy = tmp;		tmp += 9;
	#else
		cy = tmp;		tmp += 18;
	#endif
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 2 = 0xa5 [avx] or 0xae [sse2]; This is where the value of half_arr_offset36 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */

		ASSERT((radix36_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix36_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0 );		VEC_DBL_INIT(one, 1.0 );
		VEC_DBL_INIT(cc1  , c	);		VEC_DBL_INIT(ss1, s );
		VEC_DBL_INIT(cc2  , c2  );		VEC_DBL_INIT(ss2, s2);
		VEC_DBL_INIT(cc3m1, c3m1);		VEC_DBL_INIT(ss3, s3);
		VEC_DBL_INIT(cc4  , c4  );		VEC_DBL_INIT(ss4, s4);
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy - (intptr_t)two;	// #bytes in above sincos block of data
		tmp = two;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}
		nbytes = SZ_VD;	// sse2_rnd is a solo (in the SIMD-vector) datum
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

		In 4-way SIMD (AVX) mode, we expand this from 2^2 2-vector table entries to 2^4 4-vector entries.
		*/
		tmp = half_arr;

	#ifdef USE_AVX512
		// Each lookup-category in the 'mini-tables' used in AVX mode balloons from 16x32-bytes to 64x64-bytes,
		// so switch to an opmask-based scheme which starts with e.g. a broadcast constant and onditional doubling.
		// Here are the needed consts and opmasks:
		// [1] Fwd-wt multipliers: Init = 0.50 x 8, anytime AVX-style lookup into 1st table below would have bit = 0, double the corr. datum
		// [2] Inv-wt multipliers: Init = 0.25 x 8, anytime AVX-style lookup into 2nd table below would have bit = 0, double the corr. datum
		// [3] Fwd-base mults: Init = base[0] x 8, anytime AVX-style lookup into 3rd table below would have bit = 1, double the corr. datum
		// [4] Inv-base mults: Init = binv[1] x 8, anytime AVX-style lookup into 4th table below would have bit = 0, double the corr. datum
		// [5] [LOACC] Init = wts_mult[1] x 8, anytime AVX-style lookup into 5th table below would have bit = 0, double the corr. datum
		// [6] [LOACC] Init = inv_mult[0] x 8, anytime AVX-style lookup into 6th table below would have bit = 1, double the corr. datum
		nbytes = 0;
	#elif defined(USE_AVX)
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
		// In LOACC mode, put wts_mult and their inverses in the first 32 slots below in place of the 1/2-stuff:
		/* wts_mult:*/
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		/* inv_mult: */
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		nbytes = 96 << L2_SZ_VD;

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
		// In LOACC mode, put wts_mult and their inverses in the first 8 slots below in place of the 1/2-stuff:
		/* wts_mult:*/
		ctmp->re = wts_mult[0];	ctmp->im = wts_mult[0];	++ctmp;
		ctmp->re = wts_mult[1];	ctmp->im = wts_mult[0];	++ctmp;
		ctmp->re = wts_mult[0];	ctmp->im = wts_mult[1];	++ctmp;
		ctmp->re = wts_mult[1];	ctmp->im = wts_mult[1];	++ctmp;
		/* inv_mult:*/
		ctmp->re = inv_mult[0];	ctmp->im = inv_mult[0];	++ctmp;
		ctmp->re = inv_mult[1];	ctmp->im = inv_mult[0];	++ctmp;
		ctmp->re = inv_mult[0];	ctmp->im = inv_mult[1];	++ctmp;
		ctmp->re = inv_mult[1];	ctmp->im = inv_mult[1];	++ctmp;
		nbytes = 24 << L2_SZ_VD;

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

		nbytes = 4 << L2_SZ_VD;

	#ifdef USE_AVX512
	  #ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
	  #endif
		nbytes += 128;
	#elif defined(USE_AVX)
	  #ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;
	  #endif
		nbytes += 64;
	#endif

		// Propagate the above consts to the remaining threads:
		tmp = (vec_dbl *)sm_ptr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	// For large radices, array-access to bjmodn means only init base-ptr here:
	#ifndef MULTITHREAD
	  #ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n   + RE_IM_STRIDE);
	  #endif
	#endif

	#endif	// USE_SSE2

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
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
	#if !defined(MULTITHREAD) && (!defined(USE_SSE2) || (COMPACT_OBJ != 0))
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif
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

		ASSERT(ptr_prod == 0, "ERROR: unable to allocate one or more auxiliary arrays.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/RADIX-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodnini in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
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
		ASSERT(_bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
		#ifdef USE_SSE2
			tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].r00 + ((intptr_t)half_arr - (intptr_t)r00));
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r00      = (double *)base;
			tdat[ithread].half_arr = (double *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

		first_entry=FALSE;
	}	/* endif(first_entry) */

	// Jun 2018: If LL test and shift applied, compute target index for data-processing loop.
	// Note that only 1 thread of the carry-processing set will hit the target, but all need the same logic to check for a hit:
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY) {
		if(RES_SHIFT) {
			itmp64 = shift_word(a, n, p, RES_SHIFT, 0.0);	// Note return value (specifically high 7 bytes thereof) is an unpadded index
			target_idx = (int)(itmp64 >>  8);	// This still needs to be (mod NDIVR)'ed, but first use unmodded form to compute needed DWT weights
			// Compute wt = 2^(target_idx*sw % n)/n and its reciprocal:
			uint32 sw_idx_modn = ((uint64)target_idx*sw) % n;	// N is 32-bit, so only use 64-bit to hold intermediate product
			double target_wtfwd = pow(2.0, sw_idx_modn*0.5*n2inv);	// 0.5*n2inv = 0.5/(n/2) = 1.0/n
			target_set = target_idx*ndivr_inv;	// Which of the [RADIX] independent sub-carry-chains contains the target index?
			target_idx -= target_set*NDIVR;		// Fast computation of target_idx = (target_idx % NDIVR)
			// Now compute the doubles-pointer offset of the target double w.r.to the SIMD s1p00-... data layout:
			tidx_mod_stride = target_idx & (stride-1);	// Stride a power of 2, so can use AND-minus-1 for mod
			target_idx -= tidx_mod_stride;
		//	printf("Iter %d: cy_shift = %d, target_idx,tidx_mod_stride,target_set = %d,%d,%d\n",iter,(itmp64 & 255),target_idx,tidx_mod_stride,target_set);
		#ifdef USE_AVX512
			tidx_mod_stride = br16[tidx_mod_stride];
		#elif defined(USE_AVX)
			tidx_mod_stride = br8[tidx_mod_stride];
		#elif defined(USE_SSE2)
			tidx_mod_stride = br4[tidx_mod_stride];
		#endif
			target_set = (target_set<<(L2_SZ_VD-2)) + tidx_mod_stride;
			target_cy  = target_wtfwd * (-(int)(2u << (itmp64 & 255)));
		} else {
			target_idx = target_set = 0;
			target_cy = -2.0;
		}
	}

/*...The radix-36 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy[i][ithread] = 0;
		}
	}
  #if 0	//ndef USE_SSE2	*** v20: Non-SIMD builds now also support shifted-residue
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy[0][0] = -2;
	}
  #endif
	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

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
			_jhi[ithread] = _jstart[ithread] + jhi_wrap;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
	}

#ifdef USE_PTHREAD
	for(ithread = 0; ithread < CY_THREADS; ++ithread) { tdat[ithread].iter = iter; }
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		// Carry-injection location for the shifted-residue -2 addend is only needed for full pass:
		if(full_pass) {
			tdat[0].target_idx = target_idx;
			tdat[0].target_set = target_set;
			tdat[0].target_cy  = target_cy;
		} else {
			tdat[0].target_idx = -1;
			tdat[0].target_set = 0;
			tdat[0].target_cy  = 0;
		}
		// Copy to the remaining threads:
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			tdat[ithread].target_idx = tdat[0].target_idx;
			tdat[ithread].target_set = tdat[0].target_set;
			tdat[ithread].target_cy  = tdat[0].target_cy;
		}
	}
#endif
#ifdef USE_SSE2

	tmp = max_err;	VEC_DBL_INIT(tmp, 0.0);
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, SZ_VD);
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
		ASSERT(tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");

		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].nwt == nwt, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = 0.0;
		tdat[ithread].scale = scale;
		tdat[ithread].prp_mult = prp_mult;

	// pointer data:
		tdat[ithread].arrdat = a;			/* Main data array */
		ASSERT(tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].si  == si, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		ASSERT(((tmp-1)->d0 == base[0] && (tmp-1)->d1 == baseinv[1] && (tmp-1)->d2 == wts_mult[1] && (tmp-1)->d3 == inv_mult[0]), "thread-local memcheck failed!");
	  #else
		ASSERT(((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	  #endif
	  #ifdef USE_AVX512
			/* No-Op */
	  #elif defined(USE_AVX)
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	#endif
		/* init carries: */
		for(i = 0; i < RADIX; i++) {
			tdat[ithread].cy[i] = _cy[i][ithread];
		}
	}
#endif

/******************* AVX debug stuff: *******************/
#if 0
	int ipad;
	ASSERT(p01 >= 16, "Smallest array-stride must be large enough to hold an AVX-512 vec_cmplx!");
	// Use RNG to populate data array:
	rng_isaac_init(TRUE);
	double dtmp = 128.0*1024.0*1024.0*1024.0*1024.0;	// 2^47
	for(i = 0; i < n; i += 16) {
		ipad = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		// All the inits are w.r.to an un-SIMD-rearranged ...,re,im,re,im,... pattern:
	#ifdef USE_AVX512
		a[ipad+br16[ 0]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re0
		a[ipad+br16[ 1]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im0
		a[ipad+br16[ 2]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re1
		a[ipad+br16[ 3]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im1
		a[ipad+br16[ 4]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re2
		a[ipad+br16[ 5]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im2
		a[ipad+br16[ 6]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re3
		a[ipad+br16[ 7]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im3
		a[ipad+br16[ 8]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re4
		a[ipad+br16[ 9]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im4
		a[ipad+br16[10]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re5
		a[ipad+br16[11]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im5
		a[ipad+br16[12]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re6
		a[ipad+br16[13]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im6
		a[ipad+br16[14]] = dtmp*rng_isaac_rand_double_norm_pm1();	// re7
		a[ipad+br16[15]] = dtmp*rng_isaac_rand_double_norm_pm1();	// im7
	#elif defined(USE_AVX)
		a[ipad+br8[0]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re0
		a[ipad+br8[1]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im0
		a[ipad+br8[2]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re1
		a[ipad+br8[3]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im1
		a[ipad+br8[4]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re2
		a[ipad+br8[5]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im2
		a[ipad+br8[6]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re3
		a[ipad+br8[7]  ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im3
		a[ipad+br8[0]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re4
		a[ipad+br8[1]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im4
		a[ipad+br8[2]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re5
		a[ipad+br8[3]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im5
		a[ipad+br8[4]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re6
		a[ipad+br8[5]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im6
		a[ipad+br8[6]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re7
		a[ipad+br8[7]+8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im7
	#elif defined(USE_SSE2)
		a[ipad+br4[0]   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re0
		a[ipad+br4[1]   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im0
		a[ipad+br4[2]   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re1
		a[ipad+br4[3]   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im1
		a[ipad+br4[0]+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// re2
		a[ipad+br4[1]+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// im2
		a[ipad+br4[2]+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// re3
		a[ipad+br4[3]+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// im3
		a[ipad+br4[0]+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re4
		a[ipad+br4[1]+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im4
		a[ipad+br4[2]+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re5
		a[ipad+br4[3]+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im5
		a[ipad+br4[0]+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// re6
		a[ipad+br4[1]+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// im6
		a[ipad+br4[2]+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// re7
		a[ipad+br4[3]+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// im7
	#else
		a[ipad+0   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re0
		a[ipad+1   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im0
		a[ipad+2   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// re1
		a[ipad+3   ] = dtmp*rng_isaac_rand_double_norm_pm1();	// im1
		a[ipad+0+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// re2
		a[ipad+1+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// im2
		a[ipad+2+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// re3
		a[ipad+3+ 4] = dtmp*rng_isaac_rand_double_norm_pm1();	// im3
		a[ipad+0+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re4
		a[ipad+1+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im4
		a[ipad+2+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// re5
		a[ipad+3+ 8] = dtmp*rng_isaac_rand_double_norm_pm1();	// im5
		a[ipad+0+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// re6
		a[ipad+1+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// im6
		a[ipad+2+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// re7
		a[ipad+3+12] = dtmp*rng_isaac_rand_double_norm_pm1();	// im7
	#endif
  #if 0	// print DFT inputs in linear array fashion, 1st w.r.to non-SIMD {re,im,re,im,...} layout, then w.r.to actual SIMD layout:
	#ifdef USE_AVX512
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 0, a[ipad+br16[ 0]],ipad+ 0, a[ipad+ 0]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 1, a[ipad+br16[ 1]],ipad+ 1, a[ipad+ 1]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 2, a[ipad+br16[ 2]],ipad+ 2, a[ipad+ 2]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 3, a[ipad+br16[ 3]],ipad+ 3, a[ipad+ 3]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 4, a[ipad+br16[ 4]],ipad+ 4, a[ipad+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 5, a[ipad+br16[ 5]],ipad+ 5, a[ipad+ 5]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 6, a[ipad+br16[ 6]],ipad+ 6, a[ipad+ 6]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 7, a[ipad+br16[ 7]],ipad+ 7, a[ipad+ 7]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 8, a[ipad+br16[ 8]],ipad+ 8, a[ipad+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i, 9, a[ipad+br16[ 9]],ipad+ 9, a[ipad+ 9]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,10, a[ipad+br16[10]],ipad+10, a[ipad+10]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,11, a[ipad+br16[11]],ipad+11, a[ipad+11]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,12, a[ipad+br16[12]],ipad+12, a[ipad+12]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,13, a[ipad+br16[13]],ipad+13, a[ipad+13]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,14, a[ipad+br16[14]],ipad+14, a[ipad+14]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,15, a[ipad+br16[15]],ipad+15, a[ipad+15]);
	#elif defined(USE_AVX)
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0  ,a[ipad+br8[0]  ],ipad+0  ,a[ipad+0  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1  ,a[ipad+br8[1]  ],ipad+1  ,a[ipad+1  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2  ,a[ipad+br8[2]  ],ipad+2  ,a[ipad+2  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3  ,a[ipad+br8[3]  ],ipad+3  ,a[ipad+3  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,4  ,a[ipad+br8[4]  ],ipad+4  ,a[ipad+4  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,5  ,a[ipad+br8[5]  ],ipad+5  ,a[ipad+5  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,6  ,a[ipad+br8[6]  ],ipad+6  ,a[ipad+6  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,7  ,a[ipad+br8[7]  ],ipad+7  ,a[ipad+7  ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0+8,a[ipad+br8[0]+8],ipad+0+8,a[ipad+0+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1+8,a[ipad+br8[1]+8],ipad+1+8,a[ipad+1+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2+8,a[ipad+br8[2]+8],ipad+2+8,a[ipad+2+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3+8,a[ipad+br8[3]+8],ipad+3+8,a[ipad+3+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,4+8,a[ipad+br8[4]+8],ipad+4+8,a[ipad+4+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,5+8,a[ipad+br8[5]+8],ipad+5+8,a[ipad+5+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,6+8,a[ipad+br8[6]+8],ipad+6+8,a[ipad+6+8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,7+8,a[ipad+br8[7]+8],ipad+7+8,a[ipad+7+8]);
	#elif defined(USE_SSE2)
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0   ,a[ipad+br4[0]   ],ipad+0   ,a[ipad+0   ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1   ,a[ipad+br4[1]   ],ipad+1   ,a[ipad+1   ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2   ,a[ipad+br4[2]   ],ipad+2   ,a[ipad+2   ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3   ,a[ipad+br4[3]   ],ipad+3   ,a[ipad+3   ]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0+ 4,a[ipad+br4[0]+ 4],ipad+0+ 4,a[ipad+0+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1+ 4,a[ipad+br4[1]+ 4],ipad+1+ 4,a[ipad+1+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2+ 4,a[ipad+br4[2]+ 4],ipad+2+ 4,a[ipad+2+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3+ 4,a[ipad+br4[3]+ 4],ipad+3+ 4,a[ipad+3+ 4]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0+ 8,a[ipad+br4[0]+ 8],ipad+0+ 8,a[ipad+0+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1+ 8,a[ipad+br4[1]+ 8],ipad+1+ 8,a[ipad+1+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2+ 8,a[ipad+br4[2]+ 8],ipad+2+ 8,a[ipad+2+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3+ 8,a[ipad+br4[3]+ 8],ipad+3+ 8,a[ipad+3+ 8]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,0+12,a[ipad+br4[0]+12],ipad+0+12,a[ipad+0+12]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,1+12,a[ipad+br4[1]+12],ipad+1+12,a[ipad+1+12]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,2+12,a[ipad+br4[2]+12],ipad+2+12,a[ipad+2+12]);
		printf("A[%3d][%2d] = %20.10e; SIMD: A[%2d] = %20.10e\n",i,3+12,a[ipad+br4[3]+12],ipad+3+12,a[ipad+3+12]);
	#else
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,0   ,a[ipad+0   ],ipad+0   );
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,1   ,a[ipad+1   ],ipad+1   );
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,2   ,a[ipad+2   ],ipad+2   );
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,3   ,a[ipad+3   ],ipad+3   );
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,0+ 4,a[ipad+0+ 4],ipad+0+ 4);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,1+ 4,a[ipad+1+ 4],ipad+1+ 4);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,2+ 4,a[ipad+2+ 4],ipad+2+ 4);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,3+ 4,a[ipad+3+ 4],ipad+3+ 4);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,0+ 8,a[ipad+0+ 8],ipad+0+ 8);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,1+ 8,a[ipad+1+ 8],ipad+1+ 8);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,2+ 8,a[ipad+2+ 8],ipad+2+ 8);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,3+ 8,a[ipad+3+ 8],ipad+3+ 8);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,0+12,a[ipad+0+12],ipad+0+12);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,1+12,a[ipad+1+12],ipad+1+12);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,2+12,a[ipad+2+12],ipad+2+12);
		printf("A[%3d][%2d] = %20.10e; ipad: %2d\n",i,3+12,a[ipad+3+12],ipad+3+12);
	#endif
	if(i+16 >= n) exit(0);	// If printing the above inputs, exit immediately.
  #endif
	}
#endif
/********************************************************/

#ifdef USE_PTHREAD

	// If also using main thread to do work units, that task-dispatch occurs after all the threadpool-task launches:
	for(ithread = 0; ithread < pool_work_units; ithread++)
	{
		task_control.data = (void*)(&tdat[ithread]);
		threadpool_add_task(tpool, &task_control, task_is_blocking);

#else

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		if(full_pass) maxerr = 0.0;
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
	#ifdef USE_AVX512
		tmp = cy;
		for(l = 0; l < RADIX; l += 8, ++tmp) {
			tmp->d0 = _cy[l  ][ithread];
			tmp->d1 = _cy[l+1][ithread];
			tmp->d2 = _cy[l+2][ithread];
			tmp->d3 = _cy[l+3][ithread];
			tmp->d4 = _cy[l+4][ithread];
			tmp->d5 = _cy[l+5][ithread];
			tmp->d6 = _cy[l+6][ithread];
			tmp->d7 = _cy[l+7][ithread];
		}
	#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
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
	#ifdef USE_AVX512
		tmp = cy;
		for(l = 0; l < RADIX; l += 8, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
			_cy[l+2][ithread] = tmp->d2;
			_cy[l+3][ithread] = tmp->d3;
			_cy[l+4][ithread] = tmp->d4;
			_cy[l+5][ithread] = tmp->d5;
			_cy[l+6][ithread] = tmp->d6;
			_cy[l+7][ithread] = tmp->d7;
		}
		if(full_pass) {
			t0 = MAX(max_err->d0,max_err->d1);
			t1 = MAX(max_err->d2,max_err->d3);
			t2 = MAX(max_err->d4,max_err->d5);
			t3 = MAX(max_err->d6,max_err->d7);
			maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
		}
	#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
			_cy[l+2][ithread] = tmp->d2;
			_cy[l+3][ithread] = tmp->d3;
		}
		if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
		}
		if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
	#else
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = cy[l];
		}
	#endif

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #if 0//def OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(0x0 == cy36_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}
//	printf("radix40_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		if(maxerr < tdat[ithread].maxerr) {
			maxerr = tdat[ithread].maxerr;
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
	scale = prp_mult = 1;
	j_jhi = jhi_wrap;

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
		*fracmax = maxerr;
	}
	if(dtmp != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
		mlucas_fprint(cbuf,INTERACT);
		err = ERR_CARRY;
		return(err);
	}

	return(0);
}

/****************/

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
		NDIVR = n/RADIX;

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

	for(j = 0; j < NDIVR; j += 2)
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
	  #if 1
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
	  #else	// _FMA version of DIF uses different output ordering, thus yields different operm for radix-36 DIT:
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_09_DIF_FMA(a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIF_FMA(a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIF_FMA(a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIF_FMA(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, c,s,c2,s2,c3,s3,c4,s4);
		/*...and now do 9 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
	  #endif
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
					c3  = -0.50000000000000000000,	/* cos(3*u) */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double rt,it;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;

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

	for(j = 0; j < NDIVR; j += 2)
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
	  #if 0
		jt = j1    ; jp = j2    ;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],rt,it,re);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],rt,it,re);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],rt,it,re);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],rt,it,re);
	  #else	// _FMA version of DIT uses different output ordering, thus yields different operm for radix-36 DIT:
		jt = j1    ; jp = j2    ;	RADIX_09_DIT_FMA(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p32],a[jp+p32],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28], c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIT_FMA(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p24],a[jp+p24],a[jt+p32],a[jp+p32],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16], c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIT_FMA(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p32],a[jp+p32],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08], c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIT_FMA(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p32],a[jp+p32],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt    ],a[jp    ], c,s,c2,s2,c3,s3,c4,s4);
	  #endif
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
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32;
		int poff[RADIX>>2];	// Store mults of p04 offset for loop control
		int j,j1,l;
	#ifndef USE_SSE2
		int j2,jt,jp;
	#endif
	#if !defined(USE_SSE2) || defined(USE_AVX512)
		int ntmp;
	#endif
	#ifdef USE_SSE2
		int incr;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = --|--|9 for avx512,avx,sse, respectively:
		const int incr_long = 9,incr_med = 3, incr_short = 3;
	  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, so use "goes to 11" in LOACC mode via an incr_hiacc = 1:
	  #if defined(USE_AVX512) || defined(USE_ARM_V8_SIMD)
		const int incr_hiacc = 1;
	  #else
		const int incr_hiacc = 0;
	  #endif
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(USE_SHORT_CY_CHAIN == 0)
			incr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			incr = incr_med;
		else if(USE_SHORT_CY_CHAIN == 2)
			incr = incr_short;
		else
			incr = incr_hiacc;
	#endif

	#ifndef USE_AVX
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#endif
	#ifdef USE_AVX512
		double t0,t1,t2,t3;
		struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#elif defined(USE_AVX)
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
	#ifdef USE_AVX2
		// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
		vec_dbl *rad9_iptr[9], *rad9_optr[9];
	#endif
	#ifndef USE_SSE2
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
	#if !defined(USE_SSE2) || (COMPACT_OBJ != 0)
		uint32 p0123[4];
	#endif
	#ifdef USE_SSE2
	  #ifdef USE_AVX512
		int i0,i1,i2,i3;
	  #else
		const double crnd = 3.0*0x4000000*0x2000000;
	  #endif
		int *itmp;	// Pointer into the bjmodn array
	  #if defined(USE_AVX) && !defined(USE_AVX512)
		int *itm2;	// Pointer into the bjmodn array
	  #endif
	  #ifndef USE_AVX
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	  #endif
		double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl
		  #ifdef USE_AVX2
			*two,
		  #endif
		  /* *one, */ *cc1, /* *ss1, *cc2, *ss2, *cc3m1, *ss3, *cc4, *ss4, */ *max_err, *half_arr
		  #ifndef USE_AVX512
			,*sse2_rnd
		  #endif
			,*r00
		  #if !COMPACT_OBJ
			,*r02,*r04,*r06,*r08,*r0a,*r0c,*r0e,*r0g
			,*r10,*r12,*r14,*r16,*r18,*r1a,*r1c,*r1e,*r1g
			,*r20,*r22,*r24,*r26,*r28,*r2a,*r2c,*r2e,*r2g
			,*r30,*r32,*r34,*r36,*r38,*r3a,*r3c,*r3e,*r3g
		  #endif
			,*s1p00
		  #if !COMPACT_OBJ
			,*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p09,*s1p10,*s1p11,*s1p12,*s1p13,*s1p14,*s1p15,*s1p16,*s1p17
			,*s1p18,*s1p19,*s1p20,*s1p21,*s1p22,*s1p23,*s1p24,*s1p25,*s1p26,*s1p27,*s1p28,*s1p29,*s1p30,*s1p31,*s1p32,*s1p33,*s1p34,*s1p35
		  #endif
			,*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
		vec_dbl *tmp,*tm1,*tm2;	// Non-static utility ptrs
	  #ifndef USE_AVX512
		double dtmp;
	  #endif
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
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy = thread_arg->cy, rt,it,re, temp,frac;
		struct complex t[RADIX], *tptr;	// Bizarre: If we try using [RADIX] here, gcc gives "error: storage size of ‘r’ isn’t constant"
		int *itmp;	// Pointer into the bjmodn array

	#endif	// SIMD or scalar?

	// int data:
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int target_idx = thread_arg->target_idx;
		int target_set = thread_arg->target_set;
		double target_cy  = thread_arg->target_cy;
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
		double scale = thread_arg->scale;
		double prp_mult = thread_arg->prp_mult;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		double *wts_mult = thread_arg->wts_mult;	// Const Intra-block wts-multiplier...
		double *inv_mult = thread_arg->inv_mult;	// ...and 2*(its multiplicative inverse).
		ASSERT(fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		ASSERT(fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");
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
	#if !defined(USE_SSE2) || (COMPACT_OBJ != 0)
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif
		poff[0] =   0; poff[1] = p04    ; poff[2] = p08; poff[3] = p04+p08;
		poff[4] = p16; poff[5] = p04+p16; poff[6] = p24; poff[7] = p04+p24;
		poff[8] = p32;

	#ifdef USE_SSE2
		tmp = thread_arg->r00;	tm2 = tmp + 0x48;
		r00	= tmp + 0x00;		s1p00 = tm2 + 0x00;
	  #if !COMPACT_OBJ
		r02	= tmp + 0x02;		s1p01 = tm2 + 0x02;
		r04	= tmp + 0x04;		s1p02 = tm2 + 0x04;
		r06	= tmp + 0x06;		s1p03 = tm2 + 0x06;
		r08	= tmp + 0x08;		s1p04 = tm2 + 0x08;
		r0a	= tmp + 0x0a;		s1p05 = tm2 + 0x0a;
		r0c	= tmp + 0x0c;		s1p06 = tm2 + 0x0c;
		r0e	= tmp + 0x0e;		s1p07 = tm2 + 0x0e;
		r0g	= tmp + 0x10;		s1p08 = tm2 + 0x10;
		r10	= tmp + 0x12;		s1p09 = tm2 + 0x12;
		r12	= tmp + 0x14;		s1p10 = tm2 + 0x14;
		r14	= tmp + 0x16;		s1p11 = tm2 + 0x16;
		r16	= tmp + 0x18;		s1p12 = tm2 + 0x18;
		r18	= tmp + 0x1a;		s1p13 = tm2 + 0x1a;
		r1a	= tmp + 0x1c;		s1p14 = tm2 + 0x1c;
		r1c	= tmp + 0x1e;		s1p15 = tm2 + 0x1e;
		r1e	= tmp + 0x20;		s1p16 = tm2 + 0x20;
		r1g	= tmp + 0x22;		s1p17 = tm2 + 0x22;
		r20	= tmp + 0x24;		s1p18 = tm2 + 0x24;
		r22	= tmp + 0x26;		s1p19 = tm2 + 0x26;
		r24	= tmp + 0x28;		s1p20 = tm2 + 0x28;
		r26	= tmp + 0x2a;		s1p21 = tm2 + 0x2a;
		r28	= tmp + 0x2c;		s1p22 = tm2 + 0x2c;
		r2a	= tmp + 0x2e;		s1p23 = tm2 + 0x2e;
		r2c	= tmp + 0x30;		s1p24 = tm2 + 0x30;
		r2e	= tmp + 0x32;		s1p25 = tm2 + 0x32;
		r2g	= tmp + 0x34;		s1p26 = tm2 + 0x34;
		r30	= tmp + 0x36;		s1p27 = tm2 + 0x36;
		r32	= tmp + 0x38;		s1p28 = tm2 + 0x38;
		r34	= tmp + 0x3a;		s1p29 = tm2 + 0x3a;
		r36	= tmp + 0x3c;		s1p30 = tm2 + 0x3c;
		r38	= tmp + 0x3e;		s1p31 = tm2 + 0x3e;
		r3a	= tmp + 0x40;		s1p32 = tm2 + 0x40;
		r3c	= tmp + 0x42;		s1p33 = tm2 + 0x42;
		r3e	= tmp + 0x44;		s1p34 = tm2 + 0x44;
		r3g	= tmp + 0x46;		s1p35 = tm2 + 0x46;
	  #endif
		tmp	+= 0x92;	// Extra 2 slots here for two,one below - added those late, too lazy to rejigger all the existing offsets following
	   #ifdef USE_AVX2
		two    = tmp - 2;	// AVX+ versions of Radix-32 DFT macros assume consts 2.0,1.0,sqrt2,isrt2 laid out thusly
	   #endif
		//one    = tmp - 1;
		cc1    = tmp + 0;
		//ss1    = tmp + 1;
		//cc2    = tmp + 2;
		//ss2    = tmp + 3;
		//cc3m1  = tmp + 4;
		//ss3    = tmp + 5;
		//cc4    = tmp + 6;
		//ss4    = tmp + 7;
		tmp += 0x8;
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 5;		// RADIX/8 and round up
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 9;
	  #else
		cy = tmp;		tmp += 18;
	  #endif
		max_err = tmp + 0x00;
	  #ifndef USE_AVX512
		sse2_rnd= tmp + 0x01;
	  #endif
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */

		ASSERT((r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
	  #ifndef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts:
		ASSERT((sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
	  #endif
		tmp = half_arr;
	  #ifdef USE_AVX512
		/* No-Op */
	  #elif defined(USE_AVX)
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix36_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX512
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
	  #elif defined(USE_AVX)
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;
	  #endif
	  #ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n + RE_IM_STRIDE);
	  #endif
	  #if defined(USE_AVX512) && defined(X4_ZMM)
		#warning Fused AVX_cmplx_carry_fast_errcheck_X4_ZMM macro not yet ready for primetime!
		// For fused X4_ZMM version, place a copy of bjmodn[0:3] in bjmodn[4:7] slots and increment by 7*bw (mod n):
		int bw_x7_modn = (7 * (int)sse_bw[0]) % n;
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
	#ifdef USE_AVX512
		tmp = cy;
		for(l = 0; l < RADIX; l += 8, ++tmp) {
			tmp->d0 = *(addr+l  );
			tmp->d1 = *(addr+l+1);
			tmp->d2 = *(addr+l+2);
			tmp->d3 = *(addr+l+3);
			tmp->d4 = *(addr+l+4);
			tmp->d5 = *(addr+l+5);
			tmp->d6 = *(addr+l+6);
			tmp->d7 = *(addr+l+7);
		}
	#elif defined(USE_AVX)
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
	#ifdef USE_AVX512
		tmp = cy;
		for(l = 0; l < RADIX; l += 8, ++tmp) {
			*(addr+l  ) = tmp->d0;
			*(addr+l+1) = tmp->d1;
			*(addr+l+2) = tmp->d2;
			*(addr+l+3) = tmp->d3;
			*(addr+l+4) = tmp->d4;
			*(addr+l+5) = tmp->d5;
			*(addr+l+6) = tmp->d6;
			*(addr+l+7) = tmp->d7;
		}
		t0 = MAX(max_err->d0,max_err->d1);
		t1 = MAX(max_err->d2,max_err->d3);
		t2 = MAX(max_err->d4,max_err->d5);
		t3 = MAX(max_err->d6,max_err->d7);
		maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
	#elif defined(USE_AVX)
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
