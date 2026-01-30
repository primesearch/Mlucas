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

#define RADIX 12	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

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

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (24 [SSE2] or 68 [AVX])[LOACC] added slots for half_arr lookup tables.
  // Max = (40 = 0x28 [SSE2]; 132 = 0x84 [AVX]), add to (half_arr_offset12 + RADIX) to get value of radix12_creals_in_local_store.
  #ifdef USE_AVX
	const int half_arr_offset12 = 0x38;	// + RADIX = 0x44; Used for thread local-storage-integrity checking
	const int radix12_creals_in_local_store = 0xbc;	// (half_arr_offset12 + RADIX) + 132 and round up to nearest multiple of 4
  #else
	const int half_arr_offset12 = 0x3b;	// + RADIX = 0x47; Used for thread local-storage-integrity checking
	const int radix12_creals_in_local_store = 0x64;	// (half_arr_offset12 + RADIX) + 40 and round up to nearest multiple of 4
  #endif

  #ifdef COMPILER_TYPE_MSVC
	#error Mlucas build requires GCC-compatible compiler!
  #endif

	#include "sse2_macro_gcc64.h"

#endif	/* USE_SSE2 */

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
		int _pad0;	// Pads to make sizeof this struct a multiple of 16 bytes

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

		int bjmodn00;
		int bjmodn01;
		int bjmodn02;
		int bjmodn03;
		int bjmodn04;
		int bjmodn05;
		int bjmodn06;
		int bjmodn07;
		int bjmodn08;
		int bjmodn09;
		int bjmodn10;
		int bjmodn11;
		/* carries: */
		double cy00;
		double cy01;
		double cy02;
		double cy03;
		double cy04;
		double cy05;
		double cy06;
		double cy07;
		double cy08;
		double cy09;
		double cy10;
		double cy11;
	};

#endif

/**************/

int radix12_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-12 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-12 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix12_ditN_cy_dif1";
	static	int NDIVR;
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
  #endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int i,j,j2,jhi,full_pass,khi,outer;
  #ifndef MULTITHREAD
	int j1,jstart,l;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
  #endif
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p8, nsave = 0;
	#if !defined(MULTITHREAD) && (!defined(USE_SSE2) || !defined(USE_AVX))
	static int poff[RADIX>>2];	// Store [RADIX/4] mults of p4 offset for loop control
  #endif
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	static int p0123[4];	// Store [RADIX/4] mults of p4 offset for loop control
  #endif
  #if !defined(MULTITHREAD) || defined(USE_SSE2)
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
  #endif
	static double radix_inv, n2inv;
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	double rt,it;
  #endif
	double scale
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11;
  #if defined(MULTITHREAD) && defined(USE_SSE2)
	double dtmp;
  #endif
	double maxerr = 0.0;
  #if !defined(MULTITHREAD) && defined(USE_SSE2) && !defined(USE_AVX)
	int *itmp;	// Pointer into the bjmodn array
  #endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
  #ifndef MULTITHREAD
	int col,co2,co3;
   #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
   #endif
  #endif

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error x86 SIMD build requires GCC-compatible compiler!
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *addr, *add0,*add1,*add2,*add3;	/* Addresses into array sections */
  #endif

	const double crnd = 3.0*0x4000000*0x2000000;
  #ifndef USE_AVX
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
  #endif
	static vec_dbl *cc0, *cc3, *max_err, *sse2_rnd, *half_arr, *tmp,*tm2, *two
		,*r00
  #ifndef MULTITHREAD
		,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r10,*r11
		,*s1p00,*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p09,*s1p10,*s1p11
  #endif
		;
  #if !defined(MULTITHREAD) && !defined(USE_AVX)
	static vec_dbl *tm1;
  #endif

  #ifndef MULTITHREAD
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11;
	static vec_dbl *cy00,*cy04,*cy08;
   #ifndef USE_AVX
	static vec_dbl *cy02,*cy06,*cy10;
   #endif
  #endif

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
	static task_control_t   task_control = {NULL, (void*)cy12_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
  #if PFETCH && defined(MULTITHREAD)
	double *addr, *addp;
  #endif
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double wt,wtinv,wtA,wtB,wtC;
	int jt,jp,m,m2;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11;
	double t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23;
	double temp,frac
		,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11;
#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double
	*_cy00 = 0x0,*_cy01 = 0x0,*_cy02 = 0x0,*_cy03 = 0x0,*_cy04 = 0x0,*_cy05 = 0x0,*_cy06 = 0x0,*_cy07 = 0x0,*_cy08 = 0x0,*_cy09 = 0x0,*_cy10 = 0x0,*_cy11 = 0x0;

  #ifdef USE_AVX512
	WARN(HERE, "radix12_ditN_cy_dif1: No AVX-512 support; Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(0, "radix12_ditN_cy_dif1: Fermat-mod only available for radices 7,8,9,15 and their multiples!");
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
	}

/*...initialize things upon first entry: */

	if(first_entry)
	{
		psave = p;	nsave = n;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	  #ifdef USE_SSE2	// For this small radix, all vector modes use 4-way carry macros
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
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(((intptr_t)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(((intptr_t)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix12_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix12_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 48 16-byte slots of sc_arr for temporaries, next 2 for the nontrivial complex roots,
	next 6 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
	#ifdef USE_AVX
		r00	= sc_ptr + 0x00;
	  #ifndef MULTITHREAD
		r01	= sc_ptr + 0x02;
		r02	= sc_ptr + 0x04;
		r03	= sc_ptr + 0x06;
		r04	= sc_ptr + 0x08;
		r05	= sc_ptr + 0x0a;
		r06	= sc_ptr + 0x0c;
		r07	= sc_ptr + 0x0e;
		r08	= sc_ptr + 0x10;
		r09	= sc_ptr + 0x12;
		r10	= sc_ptr + 0x14;
		r11	= sc_ptr + 0x16;
	  #endif
		tmp	= sc_ptr + 0x18;
	  #ifndef MULTITHREAD
		s1p00    = tmp + 0x00;
		s1p01    = tmp + 0x02;
		s1p02    = tmp + 0x04;
		s1p03    = tmp + 0x06;
		s1p04    = tmp + 0x08;
		s1p05    = tmp + 0x0a;
		s1p06    = tmp + 0x0c;
		s1p07    = tmp + 0x0e;
		s1p08    = tmp + 0x10;
		s1p09    = tmp + 0x12;
		s1p10    = tmp + 0x14;
		s1p11    = tmp + 0x16;
	  #endif
		two      = tmp + 0x18;
		cc3      = tmp + 0x19;
		cc0      = tmp + 0x1a;
	  #ifndef MULTITHREAD
		cy00     = tmp + 0x1b;
		cy04     = tmp + 0x1c;
		cy08     = tmp + 0x1d;
	  #endif
		max_err  = tmp + 0x1e;
		sse2_rnd = tmp + 0x1f;
		half_arr = tmp + 0x20;	/* This table needs 20x16 bytes */
								// half_ar = sc_ptr + 0x38; This is where the value of half_arr_offset12 comes from
	#else
		r00	= sc_ptr + 0x00;
	  #ifndef MULTITHREAD
		r01	= sc_ptr + 0x02;
		r02	= sc_ptr + 0x04;
		r03	= sc_ptr + 0x06;
		r04	= sc_ptr + 0x08;
		r05	= sc_ptr + 0x0a;
		r06	= sc_ptr + 0x0c;
		r07	= sc_ptr + 0x0e;
		r08	= sc_ptr + 0x10;
		r09	= sc_ptr + 0x12;
		r10	= sc_ptr + 0x14;
		r11	= sc_ptr + 0x16;
	  #endif
		tmp	= sc_ptr + 0x18;
	  #ifndef MULTITHREAD
		s1p00    = tmp + 0x00;
		s1p01    = tmp + 0x02;
		s1p02    = tmp + 0x04;
		s1p03    = tmp + 0x06;
		s1p04    = tmp + 0x08;
		s1p05    = tmp + 0x0a;
		s1p06    = tmp + 0x0c;
		s1p07    = tmp + 0x0e;
		s1p08    = tmp + 0x10;
		s1p09    = tmp + 0x12;
		s1p10    = tmp + 0x14;
		s1p11    = tmp + 0x16;
	  #endif
		two      = tmp + 0x18;
		cc3      = tmp + 0x19;
		cc0      = tmp + 0x1a;
	  #ifndef MULTITHREAD
		cy00     = tmp + 0x1b;
		cy02     = tmp + 0x1c;
		cy04     = tmp + 0x1d;
		cy06     = tmp + 0x1e;
		cy08     = tmp + 0x1f;
		cy10     = tmp + 0x20;
	  #endif
		max_err  = tmp + 0x21;
		sse2_rnd = tmp + 0x22;
		half_arr = tmp + 0x23;	/* This table needs 20x16 bytes */
								// half_arr = sc_ptr + 0x3b; This is where the value of half_arr_offset12 comes from
	#endif
		/* These remain fixed: */
		VEC_DBL_INIT(two, 2.0);
		VEC_DBL_INIT(cc3  , c3m1 );
		VEC_DBL_INIT(cc0  , s    );
		VEC_DBL_INIT(sse2_rnd,crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cc0 - (intptr_t)two + SZ_VD;	// #bytes in above sincos block of data
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

	#ifdef USE_AVX
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

	#ifndef MULTITHREAD
	  #ifdef USE_AVX
		bjmodn00 = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn00 = (int*)(sse_n   + RE_IM_STRIDE);
	  #endif
		bjmodn01 = bjmodn00 + 1;
		bjmodn02 = bjmodn01 + 1;
		bjmodn03 = bjmodn02 + 1;
		bjmodn04 = bjmodn03 + 1;
		bjmodn05 = bjmodn04 + 1;
		bjmodn06 = bjmodn05 + 1;
		bjmodn07 = bjmodn06 + 1;
		bjmodn08 = bjmodn07 + 1;
		bjmodn09 = bjmodn08 + 1;
		bjmodn10 = bjmodn09 + 1;
		bjmodn11 = bjmodn10 + 1;
	#endif

	#endif	// USE_SSE2

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].r00 + ((intptr_t)half_arr - (intptr_t)r00));
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r00      = (double *)base;
			tdat[ithread].half_arr = (double *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

	/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		p1 = NDIVR;
		p2 = p1 + p1;
		p3 = p2 + p1;
		p4 = p3 + p1;
		p8 = p4 + p4;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );

	#if !defined(MULTITHREAD) && !defined(USE_SSE2)
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
	#if !defined(MULTITHREAD) && (!defined(USE_SSE2) || !defined(USE_AVX))
		poff[0] =   0; poff[1] = p4; poff[2] = p8;
	#endif

		if(_cy00)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn00); _bjmodn00 = 0x0;
			free((void *)_bjmodn01); _bjmodn01 = 0x0;
			free((void *)_bjmodn02); _bjmodn02 = 0x0;
			free((void *)_bjmodn03); _bjmodn03 = 0x0;
			free((void *)_bjmodn04); _bjmodn04 = 0x0;
			free((void *)_bjmodn05); _bjmodn05 = 0x0;
			free((void *)_bjmodn06); _bjmodn06 = 0x0;
			free((void *)_bjmodn07); _bjmodn07 = 0x0;
			free((void *)_bjmodn08); _bjmodn08 = 0x0;
			free((void *)_bjmodn09); _bjmodn09 = 0x0;
			free((void *)_bjmodn10); _bjmodn10 = 0x0;
			free((void *)_bjmodn11); _bjmodn11 = 0x0;

			free((void *)_cy00); _cy00 = 0x0;
			free((void *)_cy01); _cy01 = 0x0;
			free((void *)_cy02); _cy02 = 0x0;
			free((void *)_cy03); _cy03 = 0x0;
			free((void *)_cy04); _cy04 = 0x0;
			free((void *)_cy05); _cy05 = 0x0;
			free((void *)_cy06); _cy06 = 0x0;
			free((void *)_cy07); _cy07 = 0x0;
			free((void *)_cy08); _cy08 = 0x0;
			free((void *)_cy09); _cy09 = 0x0;
			free((void *)_cy10); _cy10 = 0x0;
			free((void *)_cy11); _cy11 = 0x0;

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
		_bjmodn00	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn00== 0x0);
		_bjmodn01	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn01== 0x0);
		_bjmodn02	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn02== 0x0);
		_bjmodn03	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn03== 0x0);
		_bjmodn04	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn04== 0x0);
		_bjmodn05	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn05== 0x0);
		_bjmodn06	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn06== 0x0);
		_bjmodn07	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn07== 0x0);
		_bjmodn08	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn08== 0x0);
		_bjmodn09	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn09== 0x0);
		_bjmodn10	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn10== 0x0);
		_bjmodn11	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn11== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy00== 0x0);
		_cy01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy01== 0x0);
		_cy02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy02== 0x0);
		_cy03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy03== 0x0);
		_cy04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy04== 0x0);
		_cy05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy05== 0x0);
		_cy06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy06== 0x0);
		_cy07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy07== 0x0);
		_cy08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy08== 0x0);
		_cy09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy09== 0x0);
		_cy10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy10== 0x0);
		_cy11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy11== 0x0);

		ASSERT(ptr_prod == 0, "ERROR: unable to allocate one or more auxiliary arrays in radix12_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/20-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodnini in radix12_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
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

/*...The radix-12 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy00[ithread] = 0;
		_cy01[ithread] = 0;
		_cy02[ithread] = 0;
		_cy03[ithread] = 0;
		_cy04[ithread] = 0;
		_cy05[ithread] = 0;
		_cy06[ithread] = 0;
		_cy07[ithread] = 0;
		_cy08[ithread] = 0;
		_cy09[ithread] = 0;
		_cy10[ithread] = 0;
		_cy11[ithread] = 0;
	}
  #if 0	//ndef USE_SSE2	*** v20: Non-SIMD builds now also support shifted-residue
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy00[0] = -2;
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
		_bjmodn00[ithread] = _bjmodnini[ithread];
		MOD_ADD32(_bjmodn00[ithread], j, n, _bjmodn01[ithread]);
		MOD_ADD32(_bjmodn01[ithread], j, n, _bjmodn02[ithread]);
		MOD_ADD32(_bjmodn02[ithread], j, n, _bjmodn03[ithread]);
		MOD_ADD32(_bjmodn03[ithread], j, n, _bjmodn04[ithread]);
		MOD_ADD32(_bjmodn04[ithread], j, n, _bjmodn05[ithread]);
		MOD_ADD32(_bjmodn05[ithread], j, n, _bjmodn06[ithread]);
		MOD_ADD32(_bjmodn06[ithread], j, n, _bjmodn07[ithread]);
		MOD_ADD32(_bjmodn07[ithread], j, n, _bjmodn08[ithread]);
		MOD_ADD32(_bjmodn08[ithread], j, n, _bjmodn09[ithread]);
		MOD_ADD32(_bjmodn09[ithread], j, n, _bjmodn10[ithread]);
		MOD_ADD32(_bjmodn10[ithread], j, n, _bjmodn11[ithread]);

		_jstart[ithread] = ithread*NDIVR/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
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
		ASSERT(((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	#endif

		tdat[ithread].bjmodn00 = _bjmodn00[ithread];
		tdat[ithread].bjmodn01 = _bjmodn01[ithread];
		tdat[ithread].bjmodn02 = _bjmodn02[ithread];
		tdat[ithread].bjmodn03 = _bjmodn03[ithread];
		tdat[ithread].bjmodn04 = _bjmodn04[ithread];
		tdat[ithread].bjmodn05 = _bjmodn05[ithread];
		tdat[ithread].bjmodn06 = _bjmodn06[ithread];
		tdat[ithread].bjmodn07 = _bjmodn07[ithread];
		tdat[ithread].bjmodn08 = _bjmodn08[ithread];
		tdat[ithread].bjmodn09 = _bjmodn09[ithread];
		tdat[ithread].bjmodn10 = _bjmodn10[ithread];
		tdat[ithread].bjmodn11 = _bjmodn11[ithread];
		/* init carries	*/
		tdat[ithread].cy00 = _cy00[ithread];
		tdat[ithread].cy01 = _cy01[ithread];
		tdat[ithread].cy02 = _cy02[ithread];
		tdat[ithread].cy03 = _cy03[ithread];
		tdat[ithread].cy04 = _cy04[ithread];
		tdat[ithread].cy05 = _cy05[ithread];
		tdat[ithread].cy06 = _cy06[ithread];
		tdat[ithread].cy07 = _cy07[ithread];
		tdat[ithread].cy08 = _cy08[ithread];
		tdat[ithread].cy09 = _cy09[ithread];
		tdat[ithread].cy10 = _cy10[ithread];
		tdat[ithread].cy11 = _cy11[ithread];
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

	#ifdef USE_SSE2
		*bjmodn00 = _bjmodn00[ithread];
		*bjmodn01 = _bjmodn01[ithread];
		*bjmodn02 = _bjmodn02[ithread];
		*bjmodn03 = _bjmodn03[ithread];
		*bjmodn04 = _bjmodn04[ithread];
		*bjmodn05 = _bjmodn05[ithread];
		*bjmodn06 = _bjmodn06[ithread];
		*bjmodn07 = _bjmodn07[ithread];
		*bjmodn08 = _bjmodn08[ithread];
		*bjmodn09 = _bjmodn09[ithread];
		*bjmodn10 = _bjmodn10[ithread];
		*bjmodn11 = _bjmodn11[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];
		bjmodn01 = _bjmodn01[ithread];
		bjmodn02 = _bjmodn02[ithread];
		bjmodn03 = _bjmodn03[ithread];
		bjmodn04 = _bjmodn04[ithread];
		bjmodn05 = _bjmodn05[ithread];
		bjmodn06 = _bjmodn06[ithread];
		bjmodn07 = _bjmodn07[ithread];
		bjmodn08 = _bjmodn08[ithread];
		bjmodn09 = _bjmodn09[ithread];
		bjmodn10 = _bjmodn10[ithread];
		bjmodn11 = _bjmodn11[ithread];
	#endif
		/* init carries	*/
	#ifdef USE_AVX//def USE_AVX	// AVX and AVX2 both use 256-bit registers
		cy00->d0 = _cy00[ithread];	cy00->d1 = _cy01[ithread];	cy00->d2 = _cy02[ithread];	cy00->d3 = _cy03[ithread];
		cy04->d0 = _cy04[ithread];	cy04->d1 = _cy05[ithread];	cy04->d2 = _cy06[ithread];	cy04->d3 = _cy07[ithread];
		cy08->d0 = _cy08[ithread];	cy08->d1 = _cy09[ithread];	cy08->d2 = _cy10[ithread];	cy08->d3 = _cy11[ithread];
	#elif defined(USE_SSE2)
		cy00->d0 = _cy00[ithread];	cy00->d1 = _cy01[ithread];
		cy02->d0 = _cy02[ithread];	cy02->d1 = _cy03[ithread];
		cy04->d0 = _cy04[ithread];	cy04->d1 = _cy05[ithread];
		cy06->d0 = _cy06[ithread];	cy06->d1 = _cy07[ithread];
		cy08->d0 = _cy08[ithread];	cy08->d1 = _cy09[ithread];
		cy10->d0 = _cy10[ithread];	cy10->d1 = _cy11[ithread];
	#else
		cy00 = _cy00[ithread];
		cy01 = _cy01[ithread];
		cy02 = _cy02[ithread];
		cy03 = _cy03[ithread];
		cy04 = _cy04[ithread];
		cy05 = _cy05[ithread];
		cy06 = _cy06[ithread];
		cy07 = _cy07[ithread];
		cy08 = _cy08[ithread];
		cy09 = _cy09[ithread];
		cy10 = _cy10[ithread];
		cy11 = _cy11[ithread];
	#endif

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix12_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_AVX
		_cy00[ithread] = cy00->d0;	_cy01[ithread] = cy00->d1;	_cy02[ithread] = cy00->d2;	_cy03[ithread] = cy00->d3;
		_cy04[ithread] = cy04->d0;	_cy05[ithread] = cy04->d1;	_cy06[ithread] = cy04->d2;	_cy07[ithread] = cy04->d3;
		_cy08[ithread] = cy08->d0;	_cy09[ithread] = cy08->d1;	_cy10[ithread] = cy08->d2;	_cy11[ithread] = cy08->d3;
		if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		_cy00[ithread] = cy00->d0;	_cy01[ithread] = cy00->d1;
		_cy02[ithread] = cy02->d0;	_cy03[ithread] = cy02->d1;
		_cy04[ithread] = cy04->d0;	_cy05[ithread] = cy04->d1;
		_cy06[ithread] = cy06->d0;	_cy07[ithread] = cy06->d1;
		_cy08[ithread] = cy08->d0;	_cy09[ithread] = cy08->d1;
		_cy10[ithread] = cy10->d0;	_cy11[ithread] = cy10->d1;
		if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
	#else
		_cy00[ithread] = cy00;
		_cy01[ithread] = cy01;
		_cy02[ithread] = cy02;
		_cy03[ithread] = cy03;
		_cy04[ithread] = cy04;
		_cy05[ithread] = cy05;
		_cy06[ithread] = cy06;
		_cy07[ithread] = cy07;
		_cy08[ithread] = cy08;
		_cy09[ithread] = cy09;
		_cy10[ithread] = cy10;
		_cy11[ithread] = cy11;
	#endif

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #if 0//def OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(0x0 == cy12_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		if(maxerr < tdat[ithread].maxerr) {
			maxerr = tdat[ithread].maxerr;
		}

		_cy00[ithread] = tdat[ithread].cy00;
		_cy01[ithread] = tdat[ithread].cy01;
		_cy02[ithread] = tdat[ithread].cy02;
		_cy03[ithread] = tdat[ithread].cy03;
		_cy04[ithread] = tdat[ithread].cy04;
		_cy05[ithread] = tdat[ithread].cy05;
		_cy06[ithread] = tdat[ithread].cy06;
		_cy07[ithread] = tdat[ithread].cy07;
		_cy08[ithread] = tdat[ithread].cy08;
		_cy09[ithread] = tdat[ithread].cy09;
		_cy10[ithread] = tdat[ithread].cy10;
		_cy11[ithread] = tdat[ithread].cy11;
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-12 forward DIF FFT of the first block of 12 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 12 outputs of (1);
	!   (3) Reweight and perform a radix-12 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 12 elements and repeat (1-4).
	*/
	t00= _cy00[CY_THREADS - 1];
	t01= _cy01[CY_THREADS - 1];
	t02= _cy02[CY_THREADS - 1];
	t03= _cy03[CY_THREADS - 1];
	t04= _cy04[CY_THREADS - 1];
	t05= _cy05[CY_THREADS - 1];
	t06= _cy06[CY_THREADS - 1];
	t07= _cy07[CY_THREADS - 1];
	t08= _cy08[CY_THREADS - 1];
	t09= _cy09[CY_THREADS - 1];
	t10= _cy10[CY_THREADS - 1];
	t11= _cy11[CY_THREADS - 1];

	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		ASSERT(CY_THREADS > 1,"radix20_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
		_cy00[ithread] = _cy00[ithread-1];
		_cy01[ithread] = _cy01[ithread-1];
		_cy02[ithread] = _cy02[ithread-1];
		_cy03[ithread] = _cy03[ithread-1];
		_cy04[ithread] = _cy04[ithread-1];
		_cy05[ithread] = _cy05[ithread-1];
		_cy06[ithread] = _cy06[ithread-1];
		_cy07[ithread] = _cy07[ithread-1];
		_cy08[ithread] = _cy08[ithread-1];
		_cy09[ithread] = _cy09[ithread-1];
		_cy10[ithread] = _cy10[ithread-1];
		_cy11[ithread] = _cy11[ithread-1];
	}

	_cy00[0] =+t11;	/* ...The wraparound carry is here: */
	_cy01[0] = t00;
	_cy02[0] = t01;
	_cy03[0] = t02;
	_cy04[0] = t03;
	_cy05[0] = t04;
	_cy06[0] = t05;
	_cy07[0] = t06;
	_cy08[0] = t07;
	_cy09[0] = t08;
	_cy10[0] = t09;
	_cy11[0] = t10;

	full_pass = 0;
	scale = prp_mult = 1;

	jhi = 7;

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + jhi; j++)
		{
			// Generate padded version of j, since prepadding pini is thread-count unsafe:
			j2 = j + ( (j >> DAT_BITS) << PAD_BITS );
			a[j2      ] *= radix_inv;
			a[j2+p1   ] *= radix_inv;
			a[j2+p2   ] *= radix_inv;
			a[j2+p3   ] *= radix_inv;
			a[j2+p4   ] *= radix_inv;
			a[j2+p1+p4] *= radix_inv;
			a[j2+p2+p4] *= radix_inv;
			a[j2+p3+p4] *= radix_inv;
			a[j2+p8   ] *= radix_inv;
			a[j2+p1+p8] *= radix_inv;
			a[j2+p2+p8] *= radix_inv;
			a[j2+p3+p8] *= radix_inv;
		}
	}
}	/* endfor(outer) */

	t00 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00 += fabs(_cy00[0])+fabs(_cy01[0])+fabs(_cy02[0])+fabs(_cy03[0])+fabs(_cy04[0])+fabs(_cy05[0])+fabs(_cy06[0])+fabs(_cy07[0])+fabs(_cy08[0])+fabs(_cy09[0])+fabs(_cy10[0])+fabs(_cy11[0]);
		*fracmax = maxerr;
	}

	if(t00 != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
		mlucas_fprint(cbuf,INTERACT);
		err = ERR_CARRY;
		return(err);
	}

	return(0);
}

/***************/

void radix12_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-12 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Uses an optimized [radix-3] x [radix-4] transform.
*/
	int j,j1,j2;
	static int n12,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, first_entry=TRUE;
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
	double rt,it
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23;

	if(!first_entry && (n/12) != n12)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n12=n/12;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n12;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;
		p11= p10+p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-12 pass is here.	*/

      for(j=0; j < n12; j += 2)
      {
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j; // LD: double-check this additional init is what is expected (comes from radix24)
	#endif
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

/*...gather the needed data (12 64-bit complex, i.e. 24 64-bit reals) and do four radix-3 transforms...	*/

/* Twiddleless version requires us to swap inputs x4 <-> x8, x2 <-> x6, x7 <-> x11 and x1 <-> x9...
indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11
      -> 0,-3,-6,-9, 8, 5, 2,-1, 4, 1,-2,-5
      == 0, 9, 6, 3, 8, 5, 2,11, 4, 1,10, 7 modulo 12.
I.e. start out with first triplet of indices {0,4,8}, permute those according to
{0,4,8}*11%12 = {0,8,4}, then each is head of a length-4 list of indices with decrement 3.
*/
#if 0
		t1 =a[j1    ];				t2 =a[j2    ];
		t3 =a[j1+p8 ]+a[j1+p4 ];	t4 =a[j2+p8 ]+a[j2+p4 ];
		t5 =a[j1+p8 ]-a[j1+p4 ];	t6 =a[j2+p8 ]-a[j2+p4 ];
		t1 =t1+t3;					t2 =t2+t4;
		t3 =t1+c3m1*t3;				t4 =t2+c3m1*t4;
		rt =s*t5;					it =s*t6;
		t5 =t3+it;					t6 =t4-rt;
		t3 =t3-it;					t4 =t4+rt;

		t7 =a[j1+p9 ];				t8 =a[j2+p9 ];
		t9 =a[j1+p5 ]+a[j1+p1 ];	t10=a[j2+p5 ]+a[j2+p1 ];
		t11=a[j1+p5 ]-a[j1+p1 ];	t12=a[j2+p5 ]-a[j2+p1 ];
		t7 =t7+t9;					t8 =t8+t10;
		t9 =t7+c3m1*t9;				t10=t8+c3m1*t10;
		rt =s*t11;					it =s*t12;
		t11=t9+it;					t12=t10-rt;
		t9 =t9-it;					t10=t10+rt;

		t13=a[j1+p6 ];				t14=a[j2+p6 ];
		t15=a[j1+p2 ]+a[j1+p10];	t16=a[j2+p2 ]+a[j2+p10];
		t17=a[j1+p2 ]-a[j1+p10];	t18=a[j2+p2 ]-a[j2+p10];
		t13=t13+t15;				t14=t14+t16;
		t15=t13+c3m1*t15;			t16=t14+c3m1*t16;
		rt =s*t17;					it =s*t18;
		t17=t15+it;					t18=t16-rt;
		t15=t15-it;					t16=t16+rt;

		t19=a[j1+p3 ];				t20=a[j2+p3 ];
		t21=a[j1+p11]+a[j1+p7 ];	t22=a[j2+p11]+a[j2+p7 ];
		t23=a[j1+p11]-a[j1+p7 ];	t24=a[j2+p11]-a[j2+p7 ];
		t19=t19+t21;				t20=t20+t22;
		t21=t19+c3m1*t21;			t22=t20+c3m1*t22;
		rt =s*t23;					it =s*t24;
		t23=t21+it;					t24=t22-rt;
		t21=t21-it;					t22=t22+rt;
/*
!...and now do three radix-4 transforms:
*/
/*...Block 1: t1,7,13,19	*/
		rt =t13;	t13=t1 -rt;		t1 =t1 +rt;
		it =t14 ;	t14=t2 -it;		t2 =t2 +it;

		rt =t19;	t19=t7 -rt;		t7 =t7 +rt;
		it =t20;	t20=t8 -it;		t8 =t8 +it;

		a[j1    ]=t1+t7;			a[j2    ]=t2+t8;
		a[j1+p1 ]=t1-t7;			a[j2+p1 ]=t2-t8;

		a[j1+p2 ]=t13-t20;			a[j2+p2 ]=t14+t19;	/* mpy by I is inlined here...	*/
		a[j1+p3 ]=t13+t20;			a[j2+p3 ]=t14-t19;

/*...Block 2: t3,9,15,21	*/
		rt =t15;	t15=t3 -rt;		t3 =t3 +rt;
		it =t16;	t16=t4 -it;		t4 =t4 +it;

		rt =t21;	t21=t9 -rt;		t9 =t9 +rt;
		it =t22;	t22=t10-it;		t10=t10+it;

		a[j1+p9 ]=t3+t9;			a[j2+p9 ]=t4+t10;
		a[j1+p8 ]=t3-t9;			a[j2+p8 ]=t4-t10;

		a[j1+p11]=t15-t22;			a[j2+p11]=t16+t21;	/* mpy by I is inlined here...	*/
		a[j1+p10]=t15+t22;			a[j2+p10]=t16-t21;

/*...Block 3: t5,11,17,23	*/
		rt =t17;	t17=t5 -rt;		t5 =t5 +rt;
		it =t18;	t18=t6 -it;		t6 =t6 +it;

		rt =t23;	t23=t11-rt;		t11=t11+rt;
		it =t24;	t24=t12-it;		t12=t12+it;

		a[j1+p6 ]=t5+t11;			a[j2+p6 ]=t6+t12;
		a[j1+p7 ]=t5-t11;			a[j2+p7 ]=t6-t12;

		a[j1+p5 ]=t17-t24;			a[j2+p5 ]=t18+t23;	/* mpy by I is inlined here...	*/
		a[j1+p4 ]=t17+t24;			a[j2+p4 ]=t18-t23;
								/* Totals: 96 FADD, 16 FMUL	*/
#elif 0
/*
Alternatively, can use the 3-decrement version:
indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11
      -> 0, 3, 6, 9, 8,11,14,17, 4, 7,10,13
      == 0, 3, 6, 9, 8,11, 2, 5, 4, 7,10, 1 modulo 12.
I.e. start out with first triplet of indices {0,4,8}, permute those according to
{0,4,8}*11%12 = {0,8,4}, then each is head of a length-4 list of indices with decrement 3.

This requires us to swap output pairs 4/5 and 10/11 w.r.to the increment version.
*/
		t1 =a[j1    ];				t2 =a[j2    ];
		t3 =a[j1+p8 ]+a[j1+p4 ];	t4 =a[j2+p8 ]+a[j2+p4 ];
		t5 =a[j1+p8 ]-a[j1+p4 ];	t6 =a[j2+p8 ]-a[j2+p4 ];
		t1 =t1+t3;					t2 =t2+t4;
		t3 =t1+c3m1*t3;				t4 =t2+c3m1*t4;
		rt =s*t5;					it =s*t6;
		t5 =t3+it;					t6 =t4-rt;
		t3 =t3-it;					t4 =t4+rt;

		t7 =a[j1+p3 ];				t8 =a[j2+p3 ];
		t9 =a[j1+p11]+a[j1+p7 ];	t10=a[j2+p11]+a[j2+p7 ];
		t11=a[j1+p11]-a[j1+p7 ];	t12=a[j2+p11]-a[j2+p7 ];
		t7 =t7+t9;					t8 =t8+t10;
		t9 =t7+c3m1*t9;				t10=t8+c3m1*t10;
		rt =s*t11;					it =s*t12;
		t11=t9+it;					t12=t10-rt;
		t9 =t9-it;					t10=t10+rt;

		t13=a[j1+p6 ];				t14=a[j2+p6 ];
		t15=a[j1+p2 ]+a[j1+p10];	t16=a[j2+p2 ]+a[j2+p10];
		t17=a[j1+p2 ]-a[j1+p10];	t18=a[j2+p2 ]-a[j2+p10];
		t13=t13+t15;				t14=t14+t16;
		t15=t13+c3m1*t15;			t16=t14+c3m1*t16;
		rt =s*t17;					it =s*t18;
		t17=t15+it;					t18=t16-rt;
		t15=t15-it;					t16=t16+rt;

		t19=a[j1+p9 ];				t20=a[j2+p9 ];
		t21=a[j1+p5 ]+a[j1+p1 ];	t22=a[j2+p5 ]+a[j2+p1 ];
		t23=a[j1+p5 ]-a[j1+p1 ];	t24=a[j2+p5 ]-a[j2+p1 ];
		t19=t19+t21;				t20=t20+t22;
		t21=t19+c3m1*t21;			t22=t20+c3m1*t22;
		rt =s*t23;					it =s*t24;
		t23=t21+it;					t24=t22-rt;
		t21=t21-it;					t22=t22+rt;
/*
!...and now do three radix-4 transforms:
*/
/*...Block 1: t1,7,13,19	*/
		rt =t13;	t13=t1 -rt;	t1 =t1 +rt;
		it =t14 ;	t14=t2 -it;	t2 =t2 +it;

		rt =t19;	t19=t7 -rt;	t7 =t7 +rt;
		it =t20;	t20=t8 -it;	t8 =t8 +it;

		a[j1    ]=t1+t7;		a[j2    ]=t2+t8;
		a[j1+p1 ]=t1-t7;		a[j2+p1 ]=t2-t8;

		a[j1+p2 ]=t13-t20;		a[j2+p2 ]=t14+t19;	/* mpy by I is inlined here...	*/
		a[j1+p3 ]=t13+t20;		a[j2+p3 ]=t14-t19;

/*...Block 2: t3,9,15,21	*/
		rt =t15;	t15=t3 -rt;	t3 =t3 +rt;
		it =t16;	t16=t4 -it;	t4 =t4 +it;

		rt =t21;	t21=t9 -rt;	t9 =t9 +rt;
		it =t22;	t22=t10-it;	t10=t10+it;

		a[j1+p9 ]=t3+t9;		a[j2+p9 ]=t4+t10;
		a[j1+p8 ]=t3-t9;		a[j2+p8 ]=t4-t10;

		a[j1+p10]=t15-t22;		a[j2+p10]=t16+t21;	/* mpy by I is inlined here...	*/
		a[j1+p11]=t15+t22;		a[j2+p11]=t16-t21;

/*...Block 3: t5,11,17,23	*/
		rt =t17;	t17=t5 -rt;	t5 =t5 +rt;
		it =t18;	t18=t6 -it;	t6 =t6 +it;

		rt =t23;	t23=t11-rt;	t11=t11+rt;
		it =t24;	t24=t12-it;	t12=t12+it;

		a[j1+p6 ]=t5+t11;		a[j2+p6 ]=t6+t12;
		a[j1+p7 ]=t5-t11;		a[j2+p7 ]=t6-t12;

		a[j1+p4 ]=t17-t24;		a[j2+p4 ]=t18+t23;	/* mpy by I is inlined here...	*/
		a[j1+p5 ]=t17+t24;		a[j2+p5 ]=t18-t23;
								/* Totals: 96 FADD, 16 FMUL	*/
#else
		int jt,jp;
		double u1,u2,u3,u4,u5,u6
			,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09;
		// Twiddleless version requires us to swap inputs [0,4,8|1,5,9|2,6,10|3,7,11] => [0,8,4|9,5,1|6,2,10|3,11,7]:
		RADIX_03_DFT(s,c3m1,a[j1   ],a[j2   ],a[j1+p8 ],a[j2+p8 ],a[j1+p4 ],a[j2+p4 ],u1,u2,u3,u4,u5,u6,t00,t01,t02,t03,t04,t05);	jt = j1+p1; jp = j2+p1;
		RADIX_03_DFT(s,c3m1,a[jt+p8],a[jp+p8],a[jt+p4 ],a[jp+p4 ],a[jt    ],a[jp    ],u1,u2,u3,u4,u5,u6,t06,t07,t08,t09,t10,t11);	jt = j1+p2; jp = j2+p2;
		RADIX_03_DFT(s,c3m1,a[jt+p4],a[jp+p4],a[jt    ],a[jp    ],a[jt+p8 ],a[jp+p8 ],u1,u2,u3,u4,u5,u6,t12,t13,t14,t15,t16,t17);	jt = j1+p3; jp = j2+p3;
		RADIX_03_DFT(s,c3m1,a[jt   ],a[jp   ],a[jt+p8 ],a[jp+p8 ],a[jt+p4 ],a[jp+p4 ],u1,u2,u3,u4,u5,u6,t18,t19,t20,t21,t22,t23);
		// And do a trio of 4-DFTs; Required output index permutation = [0,1,2,3|9,8,b,a|6,7,5,4]:
		RADIX_04_DIF(t00,t01,t06,t07,t12,t13,t18,t19,a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],rt,it);	jt = j1+p8; jp = j2+p8;
		RADIX_04_DIF(t02,t03,t08,t09,t14,t15,t20,t21,a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],rt,it);	jt = j1+p4; jp = j2+p4;
		RADIX_04_DIF(t04,t05,t10,t11,t16,t17,t22,t23,a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],rt,it);
#endif
	}
}

/***************/

void radix12_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-12 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Uses an optimized [radix-3] x [radix-4] transform.
*/
	int j,j1,j2;
	static int n12,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, first_entry=TRUE;
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
	double rt,it
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23;

	if(!first_entry && (n/12) != n12)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n12=n/12;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n12;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;
		p11= p10+p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-12 pass is here.	*/

      for(j=0; j < n12; j += 2)
      {
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j; // LD: double-check this additional init is what is expected (comes from radix24)
	#endif
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

/*...gather the needed data (12 64-bit complex, i.e. 24 64-bit reals) and do four radix-3 transforms...	*/

	/*
	Twiddleless version requires us to swap inputs x3 <-> x9, x1 <-> x4, x7 <-> x10, x2 <-> x8.
	indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11
		  -> 0, 4, 8, 9, 1, 5, 6,10, 2, 3, 7,11 modulo 12.
	I.e. start out with first quartet of indices {0,3,6,9}, permute those according to
	  {0,3,6,9}*11%12 = {0,9,6,3}, then each is head of a length-3 list of indices with increment 4, so it seems
	that we can use either decrement 4 (see version C of this routine) or increment 4, whichever minimizes the number of index swaps.

	Remember, inputs to DIT are bit-reversed, so
	a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11] contain
	x[0, 6, 3, 9, 1, 7, 4,10, 2, 8, 5,11], which get swapped to
	x[0, 6, 9, 3, 4,10, 1, 7, 8, 2, 5,11], which means the a-indices get swapped as
	a[0, 1, 3, 2, 6, 7, 4, 5, 9, 8,10,11],
	i.e.
	swapping x3 and x9  means swapping a[2] and a[3];
	swapping x1 and x4  means swapping a[4] and a[6];
	swapping x7 and x10 means swapping a[5] and a[7];
	swapping x2 and x8  means swapping a[8] and a[9].
	*/
#if 0
	/*...Block 1:	*/
		t1 =a[j1    ];		t2 =a[j2    ];
		rt =a[j1+p1 ];		it =a[j2+p1 ];
		t3 =t1 -rt;  		t1 =t1 +rt;
		t4 =t2 -it;			t2 =t2 +it;

		t5 =a[j1+p3 ];		t6 =a[j2+p3 ];
		rt =a[j1+p2 ];		it =a[j2+p2 ];
		t7 =t5 -rt;  		t5 =t5 +rt;
		t8 =t6 -it;  		t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
				t8 =t4 +rt;	t4 =t4 -rt;

	/*...Block 2:	*/
		t9 =a[j1+p6 ];		t10=a[j2+p6 ];
		rt =a[j1+p7 ];		it =a[j2+p7 ];
		t11=t9 -rt;  		t9 =t9 +rt;
		t12=t10-it;			t10=t10+it;

		t13=a[j1+p4 ];		t14=a[j2+p4 ];
		rt =a[j1+p5 ];		it =a[j2+p5 ];
		t15=t13-rt;  		t13=t13+rt;
		t16=t14-it;			t14=t14+it;

		rt =t13;t13=t9 -rt;	t9 =t9 +rt;
		it =t14;t14=t10-it;	t10=t10+it;

		rt =t15;t15=t11-t16;t11=t11+t16;
				t16=t12+rt;	t12=t12-rt;

	/*...Block 3:	*/
		t17=a[j1+p9 ];		t18=a[j2+p9 ];
		rt =a[j1+p8 ];		it =a[j2+p8 ];
		t19=t17-rt;  		t17=t17+rt;
		t20=t18-it;			t18=t18+it;

		t21=a[j1+p10];		t22=a[j2+p10];
		rt =a[j1+p11];		it =a[j2+p11];
		t23=t21-rt;  		t21=t21+rt;
		t24=t22-it;			t22=t22+it;

		rt =t21;t21=t17-rt;	t17=t17+rt;
		it =t22;t22=t18-it;	t18=t18+it;

		rt =t23;t23=t19-t24;t19=t19+t24;
				t24=t20+rt;	t20=t20-rt;
	/*
	!...and now do four radix-3 transforms.
	*/
	/*...Block 1: t1,9,17	*/
		rt =t9;				it =t10;
		t9 =rt+t17;			t10=it+t18;
		t17=rt-t17;			t18=it-t18;
		t1 =t1+t9;			t2 =t2+t10;
		a[j1    ]=t1;		a[j2    ]=t2;
		t9 =t1+c3m1*t9;		t10=t2+c3m1*t10;
		rt =s*t17;			it =s*t18;
		a[j1+p4 ]=t9+it;	a[j2+p4 ]=t10-rt;
		a[j1+p8 ]=t9-it;	a[j2+p8 ]=t10+rt;
	/*...Block 2: t3,11,19	*/
		rt =t11;			it =t12;
		t11=rt+t19;			t12=it+t20;
		t19=rt-t19;			t20=it-t20;
		t3 =t3+t11;			t4 =t4+t12;
		a[j1+p3 ]=t3;		a[j2+p3 ]=t4;
		t11=t3+c3m1*t11;	t12=t4+c3m1*t12;
		rt =s*t19;			it =s*t20;
		a[j1+p7 ]=t11+it;	a[j2+p7 ]=t12-rt;
		a[j1+p11]=t11-it;	a[j2+p11]=t12+rt;
	/*...Block 3: t5,13,21	*/
		rt =t13;			it =t14;
		t13=rt+t21;			t14=it+t22;
		t21=rt-t21;			t22=it-t22;
		t5 =t5+t13;			t6 =t6+t14;
		a[j1+p6 ]=t5;		a[j2+p6 ]=t6;
		t13=t5+c3m1*t13;	t14=t6+c3m1*t14;
		rt =s*t21;			it =s*t22;
		a[j1+p10]=t13+it;	a[j2+p10]=t14-rt;
		a[j1+p2 ]=t13-it;	a[j2+p2 ]=t14+rt;
	/*...Block 3: t7,15,23	*/
		rt =t15;			it =t16;
		t15=rt+t23;			t16=it+t24;
		t23=rt-t23;			t24=it-t24;
		t7 =t7+t15;			t8 =t8+t16;
		a[j1+p9 ]=t7;		a[j2+p9 ]=t8;
		t15=t7+c3m1*t15;	t16=t8+c3m1*t16;
		rt =s*t23;			it =s*t24;
		a[j1+p1 ]=t15+it;	a[j2+p1 ]=t16-rt;
		a[j1+p5 ]=t15-it;	a[j2+p5 ]=t16+rt;
						/* Totals: 96 FADD, 16 FMUL.	*/
#else
		int jt,jp;
		double u1,u2,u3,u4,u5,u6
			,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09;
		// Do a trio of 4-DFTs; combined DIT input-scramble array = [0,1,3,2|9,8,a,b|6,7,4,5]:
		RADIX_04_DIT(a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p3],a[j2+p3],a[j1+p2],a[j2+p2],t00,t01,t06,t07,t12,t13,t18,t19,rt,it);	jt = j1+p8; jp = j2+p8;
		RADIX_04_DIT(a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],t02,t03,t08,t09,t14,t15,t20,t21,rt,it);	jt = j1+p4; jp = j2+p4;
		RADIX_04_DIT(a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],t04,t05,t10,t11,t16,t17,t22,t23,rt,it);
		/* And now do a quarter of 3-DFTs. Required output index permutation: [0,1,2,3,4,5,6,7,8,9,a,b] => [0,3,6,9,4,7,a,1,8,b,2,5],
		i.e. in terms of our trial-output ordering into trios, [0,4,8|1,5,9|2,6,a|3,7,b] => [0,4,8|3,7,b|6,a,2|9,1,5]:
		*/
		RADIX_03_DFT(s,c3m1,t00,t01,t02,t03,t04,t05,u1,u2,u3,u4,u5,u6,a[j1   ],a[j2   ],a[j1+p4],a[j2+p4],a[j1+p8],a[j2+p8]);	jt = j1+p3; jp = j2+p3;
		RADIX_03_DFT(s,c3m1,t06,t07,t08,t09,t10,t11,u1,u2,u3,u4,u5,u6,a[jt   ],a[jp   ],a[jt+p4],a[jp+p4],a[jt+p8],a[jp+p8]);	jt = j1+p2; jp = j2+p2;
		RADIX_03_DFT(s,c3m1,t12,t13,t14,t15,t16,t17,u1,u2,u3,u4,u5,u6,a[jt+p4],a[jp+p4],a[jt+p8],a[jp+p8],a[jt   ],a[jp   ]);	jt = j1+p1; jp = j2+p1;
		RADIX_03_DFT(s,c3m1,t18,t19,t20,t21,t22,t23,u1,u2,u3,u4,u5,u6,a[jt+p8],a[jp+p8],a[jt   ],a[jp   ],a[jt+p4],a[jp+p4]);
#endif
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy12_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p8;
	#if !defined(USE_SSE2) || !defined(USE_AVX)
		int poff[RADIX>>2];	// Store [RADIX/4] mults of p4 offset for loop control
	#endif
	#ifndef USE_SSE2
		int p0123[4];
	#endif
		int j,j1,l;
	#ifndef USE_SSE2
		int j2;
	#endif
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
	#ifndef USE_SSE2
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif

	#ifdef USE_SSE2

		double *addr, *add0,*add1,*add2,*add3;
		const double crnd = 3.0*0x4000000*0x2000000;
	  #ifndef USE_AVX
		int *itmp;	// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	  #endif
		vec_dbl /* *cc0, */ *cc3, *max_err, *sse2_rnd, *half_arr, *tmp/* , *two */
			,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r10,*r11
			,*s1p00,*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p09,*s1p10,*s1p11;
	  #ifndef USE_AVX
		vec_dbl *tm1,*tm2;
	  #endif
		vec_dbl *cy00,*cy04,*cy08;
	  #ifndef USE_AVX
		vec_dbl *cy02,*cy06,*cy10;
	  #endif
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11;
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
						s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
		double *base, *baseinv;
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		double rt,it,temp,frac
			,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23
			,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11;
		int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11;

	#endif

		struct cy_thread_data_t* thread_arg = targ;
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
		p1 = NDIVR;
		p2 = p1 + p1;
		p3 = p2 + p1;
		p4 = p3 + p1;
		p8 = p4 + p4;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );

	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif

	#if !defined(USE_SSE2) || !defined(USE_AVX)
		poff[0] =   0; poff[1] = p4; poff[2] = p8;
	#endif

	#ifdef USE_SSE2
		r00	= thread_arg->r00;
	#ifdef USE_AVX
								tmp	= r00 + 0x18;
								s1p00 = tmp + 0x00;		//two     = tmp + 0x18;
		r01	= r00 + 0x02;		s1p01 = tmp + 0x02;		cc3     = tmp + 0x19;
		r02	= r00 + 0x04;		s1p02 = tmp + 0x04;		//cc0     = tmp + 0x1a;
		r03	= r00 + 0x06;		s1p03 = tmp + 0x06;		cy00    = tmp + 0x1b;
		r04	= r00 + 0x08;		s1p04 = tmp + 0x08;		cy04    = tmp + 0x1c;
		r05	= r00 + 0x0a;		s1p05 = tmp + 0x0a;		cy08    = tmp + 0x1d;
		r06	= r00 + 0x0c;		s1p06 = tmp + 0x0c;
		r07	= r00 + 0x0e;		s1p07 = tmp + 0x0e;
		r08	= r00 + 0x10;		s1p08 = tmp + 0x10;
		r09	= r00 + 0x12;		s1p09 = tmp + 0x12;		max_err = tmp + 0x1e;
		r10	= r00 + 0x14;		s1p10 = tmp + 0x14;		sse2_rnd= tmp + 0x1f;
		r11	= r00 + 0x16;		s1p11 = tmp + 0x16;		half_arr= tmp + 0x20;	/* This table needs 20x16 bytes */
																// half_ar = r00 + 0x38; This is where the value of half_arr_offset12 comes from
	#else
								tmp	= r00 + 0x18;
								s1p00 = tmp + 0x00;		//two     = tmp + 0x18;
		r01	= r00 + 0x02;		s1p01 = tmp + 0x02;		cc3     = tmp + 0x19;
		r02	= r00 + 0x04;		s1p02 = tmp + 0x04;		//cc0     = tmp + 0x1a;
		r03	= r00 + 0x06;		s1p03 = tmp + 0x06;		cy00    = tmp + 0x1b;
		r04	= r00 + 0x08;		s1p04 = tmp + 0x08;		cy02    = tmp + 0x1c;
		r05	= r00 + 0x0a;		s1p05 = tmp + 0x0a;		cy04    = tmp + 0x1d;
		r06	= r00 + 0x0c;		s1p06 = tmp + 0x0c;		cy06    = tmp + 0x1e;
		r07	= r00 + 0x0e;		s1p07 = tmp + 0x0e;		cy08    = tmp + 0x1f;
		r08	= r00 + 0x10;		s1p08 = tmp + 0x10;		cy10    = tmp + 0x20;
		r09	= r00 + 0x12;		s1p09 = tmp + 0x12;		max_err = tmp + 0x21;
		r10	= r00 + 0x14;		s1p10 = tmp + 0x14;		sse2_rnd= tmp + 0x22;
		r11	= r00 + 0x16;		s1p11 = tmp + 0x16;		half_arr= tmp + 0x23;	/* This table needs 20x16 bytes */
																// half_arr = r00 + 0x3b; This is where the value of half_arr_offset12 comes from
	#endif
		ASSERT((r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT((sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
		tmp = half_arr;
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix12_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;

		bjmodn00 = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn00 = (int*)(sse_n + RE_IM_STRIDE);
	  #endif
		bjmodn01 = bjmodn00 +  1;
		bjmodn02 = bjmodn00 +  2;
		bjmodn03 = bjmodn00 +  3;
		bjmodn04 = bjmodn00 +  4;
		bjmodn05 = bjmodn00 +  5;
		bjmodn06 = bjmodn00 +  6;
		bjmodn07 = bjmodn00 +  7;
		bjmodn08 = bjmodn00 +  8;
		bjmodn09 = bjmodn00 +  9;
		bjmodn10 = bjmodn00 + 10;
		bjmodn11 = bjmodn00 + 11;
	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->r00     ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		/* Init DWT-indices: */	/* init carries	*/
	#ifdef USE_AVX
		*bjmodn00 = thread_arg->bjmodn00;		cy00->d0 = thread_arg->cy00;
		*bjmodn01 = thread_arg->bjmodn01;		cy00->d1 = thread_arg->cy01;
		*bjmodn02 = thread_arg->bjmodn02;		cy00->d2 = thread_arg->cy02;
		*bjmodn03 = thread_arg->bjmodn03;		cy00->d3 = thread_arg->cy03;
		*bjmodn04 = thread_arg->bjmodn04;		cy04->d0 = thread_arg->cy04;
		*bjmodn05 = thread_arg->bjmodn05;		cy04->d1 = thread_arg->cy05;
		*bjmodn06 = thread_arg->bjmodn06;		cy04->d2 = thread_arg->cy06;
		*bjmodn07 = thread_arg->bjmodn07;		cy04->d3 = thread_arg->cy07;
		*bjmodn08 = thread_arg->bjmodn08;		cy08->d0 = thread_arg->cy08;
		*bjmodn09 = thread_arg->bjmodn09;		cy08->d1 = thread_arg->cy09;
		*bjmodn10 = thread_arg->bjmodn10;		cy08->d2 = thread_arg->cy10;
		*bjmodn11 = thread_arg->bjmodn11;		cy08->d3 = thread_arg->cy11;
	#elif defined(USE_SSE2)
		*bjmodn00 = thread_arg->bjmodn00;		cy00->d0 = thread_arg->cy00;
		*bjmodn01 = thread_arg->bjmodn01;		cy00->d1 = thread_arg->cy01;
		*bjmodn02 = thread_arg->bjmodn02;		cy02->d0 = thread_arg->cy02;
		*bjmodn03 = thread_arg->bjmodn03;		cy02->d1 = thread_arg->cy03;
		*bjmodn04 = thread_arg->bjmodn04;		cy04->d0 = thread_arg->cy04;
		*bjmodn05 = thread_arg->bjmodn05;		cy04->d1 = thread_arg->cy05;
		*bjmodn06 = thread_arg->bjmodn06;		cy06->d0 = thread_arg->cy06;
		*bjmodn07 = thread_arg->bjmodn07;		cy06->d1 = thread_arg->cy07;
		*bjmodn08 = thread_arg->bjmodn08;		cy08->d0 = thread_arg->cy08;
		*bjmodn09 = thread_arg->bjmodn09;		cy08->d1 = thread_arg->cy09;
		*bjmodn10 = thread_arg->bjmodn10;		cy10->d0 = thread_arg->cy10;
		*bjmodn11 = thread_arg->bjmodn11;		cy10->d1 = thread_arg->cy11;
	#else
		bjmodn00 = thread_arg->bjmodn00;		cy00 = thread_arg->cy00;
		bjmodn01 = thread_arg->bjmodn01;		cy01 = thread_arg->cy01;
		bjmodn02 = thread_arg->bjmodn02;		cy02 = thread_arg->cy02;
		bjmodn03 = thread_arg->bjmodn03;		cy03 = thread_arg->cy03;
		bjmodn04 = thread_arg->bjmodn04;		cy04 = thread_arg->cy04;
		bjmodn05 = thread_arg->bjmodn05;		cy05 = thread_arg->cy05;
		bjmodn06 = thread_arg->bjmodn06;		cy06 = thread_arg->cy06;
		bjmodn07 = thread_arg->bjmodn07;		cy07 = thread_arg->cy07;
		bjmodn08 = thread_arg->bjmodn08;		cy08 = thread_arg->cy08;
		bjmodn09 = thread_arg->bjmodn09;		cy09 = thread_arg->cy09;
		bjmodn10 = thread_arg->bjmodn10;		cy10 = thread_arg->cy10;
		bjmodn11 = thread_arg->bjmodn11;		cy11 = thread_arg->cy11;
	#endif	// SSE2 or AVX?

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix12_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_AVX
		thread_arg->cy00 = cy00->d0;
		thread_arg->cy01 = cy00->d1;
		thread_arg->cy02 = cy00->d2;
		thread_arg->cy03 = cy00->d3;
		thread_arg->cy04 = cy04->d0;
		thread_arg->cy05 = cy04->d1;
		thread_arg->cy06 = cy04->d2;
		thread_arg->cy07 = cy04->d3;
		thread_arg->cy08 = cy08->d0;
		thread_arg->cy09 = cy08->d1;
		thread_arg->cy10 = cy08->d2;
		thread_arg->cy11 = cy08->d3;
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		thread_arg->cy00 = cy00->d0;
		thread_arg->cy01 = cy00->d1;
		thread_arg->cy02 = cy02->d0;
		thread_arg->cy03 = cy02->d1;
		thread_arg->cy04 = cy04->d0;
		thread_arg->cy05 = cy04->d1;
		thread_arg->cy06 = cy06->d0;
		thread_arg->cy07 = cy06->d1;
		thread_arg->cy08 = cy08->d0;
		thread_arg->cy09 = cy08->d1;
		thread_arg->cy10 = cy10->d0;
		thread_arg->cy11 = cy10->d1;
		maxerr = MAX(max_err->d0,max_err->d1);
	#else
		thread_arg->cy00 = cy00;
		thread_arg->cy01 = cy01;
		thread_arg->cy02 = cy02;
		thread_arg->cy03 = cy03;
		thread_arg->cy04 = cy04;
		thread_arg->cy05 = cy05;
		thread_arg->cy06 = cy06;
		thread_arg->cy07 = cy07;
		thread_arg->cy08 = cy08;
		thread_arg->cy09 = cy09;
		thread_arg->cy10 = cy10;
		thread_arg->cy11 = cy11;
	#endif	// SSE2 or AVX?

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
