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

#define RADIX 24	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

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

  // For Mersenne-mod we need (16 [SSE2] or 64 [AVX]) + 4 added slots for the half_arr lookup tables.
  // Add relevant number (half_arr_offset24 + RADIX) to get required value of radix24_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset24 = 91;	// + RADIX = 115; Used for thread local-storage-integrity checking
	const int radix24_creals_in_local_store = 184;	// (half_arr_offset24 + RADIX) + 68 and round up to nearest multiple of 4
  #else
	const int half_arr_offset24 = 97;	// + RADIX = 121; Used for thread local-storage-integrity checking
	const int radix24_creals_in_local_store = 144;	// (half_arr_offset24 + RADIX) + 20 and round up to nearest multiple of 4
  #endif

	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */
	#if defined(USE_AVX) && !defined(ERR_CHECK_ALL)
		#error ERR_CHECK_ALL *required* for AVX-mode builds!
	#endif

	#define EPS 1e-10

	#ifdef COMPILER_TYPE_MSVC
		#include "sse2_macro.h"
	#endif

	#if defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix24_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix24_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

#endif

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data:
		int iter;
		int tid;
		int ndivr;
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

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
	#ifdef USE_SSE2
		vec_dbl *s1p00r;
		vec_dbl *half_arr;
	#else
		double *s1p00r;
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
		int bjmodn12;
		int bjmodn13;
		int bjmodn14;
		int bjmodn15;
		int bjmodn16;
		int bjmodn17;
		int bjmodn18;
		int bjmodn19;
		int bjmodn20;
		int bjmodn21;
		int bjmodn22;
		int bjmodn23;
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
		double cy12;
		double cy13;
		double cy14;
		double cy15;
		double cy16;
		double cy17;
		double cy18;
		double cy19;
		double cy20;
		double cy21;
		double cy22;
		double cy23;
	};

#endif

/**************/

int radix24_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-24 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-24 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	const int sz_vd = sizeof(vec_dbl);
	// lg(sizeof(vec_dbl)):
  #ifdef USE_AVX
	const int l2_sz_vd = 5;
  #else
	const int l2_sz_vd = 4;
  #endif
#endif
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p16;
	static double radix_inv, n2inv;
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
	double scale,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23;
	double dtmp, maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
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
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
   #ifndef COMPILER_TYPE_GCC
		double *add4, *add5, *add6, *add7;
   #endif
  #endif

	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	static vec_dbl *isrt2, *cc0, *cc3, *max_err, *sse2_rnd, *half_arr, *tmp,*tm2;
	static vec_dbl *s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r;
	/* Only explicitly reference the even-indexed carries in SSE2 mode: */
	static vec_dbl
		*cy00,*cy04,*cy08,*cy12,*cy16,*cy20;
  #ifndef USE_AVX
	static vec_dbl
		*cy02,*cy06,*cy10,*cy14,*cy18,*cy22;
  #endif
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23;

#endif

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy24_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
  #if PFETCH
	double *add0, *addr;
  #endif
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double wt,wtinv,wtA,wtB,wtC;
	int m,m2;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23;
	double rt,it,temp,frac
		,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r
		,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i
		,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0, *_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy_00 = 0x0,*_cy_01 = 0x0,*_cy_02 = 0x0,*_cy_03 = 0x0,*_cy_04 = 0x0,*_cy_05 = 0x0,*_cy_06 = 0x0,*_cy_07 = 0x0,*_cy_08 = 0x0,*_cy_09 = 0x0,*_cy_10 = 0x0,*_cy_11 = 0x0,*_cy_12 = 0x0,*_cy_13 = 0x0,*_cy_14 = 0x0,*_cy_15 = 0x0,*_cy_16 = 0x0,*_cy_17 = 0x0,*_cy_18 = 0x0,*_cy_19 = 0x0,*_cy_20 = 0x0,*_cy_21 = 0x0,*_cy_22 = 0x0,*_cy_23 = 0x0;

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "radix24_ditN_cy_dif1: Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/RADIX in radix24_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
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
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

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
			tdat[ithread].iter = iter;
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
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix24_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix24_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 48 16-byte slots of sc_arr for temporaries, next 2 for the doubled cos and c3m1 terms,
	next 12 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
	#ifdef USE_AVX
		s1p00r = sc_ptr + 0x00;		isrt2	= sc_ptr + 0x30;
		s1p01r = sc_ptr + 0x02;		cc3		= sc_ptr + 0x31;
		s1p02r = sc_ptr + 0x04;		cc0  	= sc_ptr + 0x32;
		s1p03r = sc_ptr + 0x06;		cy00	= sc_ptr + 0x33;
		s1p04r = sc_ptr + 0x08;		cy04	= sc_ptr + 0x34;
		s1p05r = sc_ptr + 0x0a;		cy08	= sc_ptr + 0x35;
		s1p06r = sc_ptr + 0x0c;		cy12	= sc_ptr + 0x36;
		s1p07r = sc_ptr + 0x0e;		cy16	= sc_ptr + 0x37;
		s1p08r = sc_ptr + 0x10;		cy20	= sc_ptr + 0x38;
		s1p09r = sc_ptr + 0x12;		max_err = sc_ptr + 0x39;
		s1p10r = sc_ptr + 0x14;		sse2_rnd= sc_ptr + 0x3a;
		s1p11r = sc_ptr + 0x16;		half_arr= sc_ptr + 0x3b;	/* This table needs 20x16 bytes */
		s1p12r = sc_ptr + 0x18;		// half_arr = sc_ptr + 0x3b; This is where the value of half_arr_offset24 comes from
		s1p13r = sc_ptr + 0x1a;
		s1p14r = sc_ptr + 0x1c;
		s1p15r = sc_ptr + 0x1e;
		s1p16r = sc_ptr + 0x20;
		s1p17r = sc_ptr + 0x22;
		s1p18r = sc_ptr + 0x24;
		s1p19r = sc_ptr + 0x26;
		s1p20r = sc_ptr + 0x28;
		s1p21r = sc_ptr + 0x2a;
		s1p22r = sc_ptr + 0x2c;
		s1p23r = sc_ptr + 0x2e;
	#else
		s1p00r = sc_ptr + 0x00;		isrt2	= sc_ptr + 0x30;
		s1p01r = sc_ptr + 0x02;		cc3		= sc_ptr + 0x31;
		s1p02r = sc_ptr + 0x04;		cc0  	= sc_ptr + 0x32;
		s1p03r = sc_ptr + 0x06;		cy00	= sc_ptr + 0x33;
		s1p04r = sc_ptr + 0x08;		cy02	= sc_ptr + 0x34;
		s1p05r = sc_ptr + 0x0a;		cy04	= sc_ptr + 0x35;
		s1p06r = sc_ptr + 0x0c;		cy06	= sc_ptr + 0x36;
		s1p07r = sc_ptr + 0x0e;		cy08	= sc_ptr + 0x37;
		s1p08r = sc_ptr + 0x10;		cy10	= sc_ptr + 0x38;
		s1p09r = sc_ptr + 0x12;		cy12	= sc_ptr + 0x39;
		s1p10r = sc_ptr + 0x14;		cy14	= sc_ptr + 0x3a;
		s1p11r = sc_ptr + 0x16;		cy16	= sc_ptr + 0x3b;
		s1p12r = sc_ptr + 0x18;		cy18	= sc_ptr + 0x3c;
		s1p13r = sc_ptr + 0x1a;		cy20	= sc_ptr + 0x3d;
		s1p14r = sc_ptr + 0x1c;		cy22	= sc_ptr + 0x3e;
		s1p15r = sc_ptr + 0x1e;		max_err = sc_ptr + 0x3f;
		s1p16r = sc_ptr + 0x20;		sse2_rnd= sc_ptr + 0x40;
		s1p17r = sc_ptr + 0x22;		half_arr= sc_ptr + 0x41;	/* This table needs 20x16 bytes */
		s1p18r = sc_ptr + 0x24;		// half_arr = sc_ptr + 0x41; This is where the value of half_arr_offset24 comes from
		s1p19r = sc_ptr + 0x26;
		s1p20r = sc_ptr + 0x28;
		s1p21r = sc_ptr + 0x2a;
		s1p22r = sc_ptr + 0x2c;
		s1p23r = sc_ptr + 0x2e;
	#endif
		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);
		VEC_DBL_INIT(cc3  , c3m1 );
		VEC_DBL_INIT(cc0  , s    );
		VEC_DBL_INIT(sse2_rnd,crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)cc0 - (int)isrt2 + sz_vd;	// #bytes in above sincos block of data
		tmp = isrt2;
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
		bjmodn12 = bjmodn11 + 1;
		bjmodn13 = bjmodn12 + 1;
		bjmodn14 = bjmodn13 + 1;
		bjmodn15 = bjmodn14 + 1;
		bjmodn16 = bjmodn15 + 1;
		bjmodn17 = bjmodn16 + 1;
		bjmodn18 = bjmodn17 + 1;
		bjmodn19 = bjmodn18 + 1;
		bjmodn20 = bjmodn19 + 1;
		bjmodn21 = bjmodn20 + 1;
		bjmodn22 = bjmodn21 + 1;
		bjmodn23 = bjmodn22 + 1;

	  // Now use targeted copy of just the per-iteration-need-to-re-init data in the local block;
	  // This is the older copy-whole-local-block approach:
	  #if 0//def USE_PTHREAD
		s1p00r = __r0 + cslots_in_local_store;
		/* Init thread 1-CY_THREADS's local stores and pointers: */
		for(i = 1; i < CY_THREADS; ++i) {
			/* Only care about the constants for each thread here, but easier to just copy the entire thread0 local store: */
			memcpy(s1p00r, __r0, cslots_in_local_store << l2_sz_vd);	// bytewise copy treats complex and uint64 subdata the same
			s1p00r += cslots_in_local_store;
		}
	  #endif

	#endif	// USE_SSE2

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].s1p00r = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].s1p00r + ((long)half_arr - (long)s1p00r);
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].s1p00r   = (double *)base;
			tdat[ithread].half_arr = (double *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;
		ASSERT(HERE, p16 == p08+p08, "p16 != p08+p08; radix24 ASM macro requires this!");
		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );

		if(_cy_00)	/* If it's a new exponent of a range test, need to deallocate these. */
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
			free((void *)_bjmodn12); _bjmodn12 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;
			free((void *)_bjmodn20); _bjmodn20 = 0x0;
			free((void *)_bjmodn21); _bjmodn21 = 0x0;
			free((void *)_bjmodn22); _bjmodn22 = 0x0;
			free((void *)_bjmodn23); _bjmodn23 = 0x0;

			free((void *)_cy_00); _cy_00 = 0x0;
			free((void *)_cy_01); _cy_01 = 0x0;
			free((void *)_cy_02); _cy_02 = 0x0;
			free((void *)_cy_03); _cy_03 = 0x0;
			free((void *)_cy_04); _cy_04 = 0x0;
			free((void *)_cy_05); _cy_05 = 0x0;
			free((void *)_cy_06); _cy_06 = 0x0;
			free((void *)_cy_07); _cy_07 = 0x0;
			free((void *)_cy_08); _cy_08 = 0x0;
			free((void *)_cy_09); _cy_09 = 0x0;
			free((void *)_cy_10); _cy_10 = 0x0;
			free((void *)_cy_11); _cy_11 = 0x0;
			free((void *)_cy_12); _cy_12 = 0x0;
			free((void *)_cy_13); _cy_13 = 0x0;
			free((void *)_cy_14); _cy_14 = 0x0;
			free((void *)_cy_15); _cy_15 = 0x0;
			free((void *)_cy_16); _cy_16 = 0x0;
			free((void *)_cy_17); _cy_17 = 0x0;
			free((void *)_cy_18); _cy_18 = 0x0;
			free((void *)_cy_19); _cy_19 = 0x0;
			free((void *)_cy_20); _cy_20 = 0x0;
			free((void *)_cy_21); _cy_21 = 0x0;
			free((void *)_cy_22); _cy_22 = 0x0;
			free((void *)_cy_23); _cy_23 = 0x0;

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
		_bjmodn12	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn12== 0x0);
		_bjmodn13	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn13== 0x0);
		_bjmodn14	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn14== 0x0);
		_bjmodn15	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn15== 0x0);
		_bjmodn16	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn16== 0x0);
		_bjmodn17	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn17== 0x0);
		_bjmodn18	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn18== 0x0);
		_bjmodn19	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn19== 0x0);
		_bjmodn20	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy_00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_00== 0x0);
		_cy_01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_01== 0x0);
		_cy_02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_02== 0x0);
		_cy_03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_03== 0x0);
		_cy_04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_04== 0x0);
		_cy_05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_05== 0x0);
		_cy_06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_06== 0x0);
		_cy_07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_07== 0x0);
		_cy_08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_08== 0x0);
		_cy_09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_09== 0x0);
		_cy_10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_10== 0x0);
		_cy_11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_11== 0x0);
		_cy_12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_12== 0x0);
		_cy_13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_13== 0x0);
		_cy_14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_14== 0x0);
		_cy_15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_15== 0x0);
		_cy_16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_16== 0x0);
		_cy_17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_17== 0x0);
		_cy_18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_18== 0x0);
		_cy_19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_19== 0x0);
		_cy_20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_20== 0x0);
		_cy_21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_21== 0x0);
		_cy_22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_22== 0x0);
		_cy_23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_23== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix24_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/24-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix24_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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

	}	/* endif(first_entry) */

/*...The radix-24 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_00[ithread] = 0;
		_cy_01[ithread] = 0;
		_cy_02[ithread] = 0;
		_cy_03[ithread] = 0;
		_cy_04[ithread] = 0;
		_cy_05[ithread] = 0;
		_cy_06[ithread] = 0;
		_cy_07[ithread] = 0;
		_cy_08[ithread] = 0;
		_cy_09[ithread] = 0;
		_cy_10[ithread] = 0;
		_cy_11[ithread] = 0;
		_cy_12[ithread] = 0;
		_cy_13[ithread] = 0;
		_cy_14[ithread] = 0;
		_cy_15[ithread] = 0;
		_cy_16[ithread] = 0;
		_cy_17[ithread] = 0;
		_cy_18[ithread] = 0;
		_cy_19[ithread] = 0;
		_cy_20[ithread] = 0;
		_cy_21[ithread] = 0;
		_cy_22[ithread] = 0;
		_cy_23[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_00[      0] = -2;
	}
	else
	{
		ASSERT(HERE,0,"Radix-24 currently only supports LL test mode!");
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
		MOD_ADD32(_bjmodn11[ithread], j, n, _bjmodn12[ithread]);
		MOD_ADD32(_bjmodn12[ithread], j, n, _bjmodn13[ithread]);
		MOD_ADD32(_bjmodn13[ithread], j, n, _bjmodn14[ithread]);
		MOD_ADD32(_bjmodn14[ithread], j, n, _bjmodn15[ithread]);
		MOD_ADD32(_bjmodn15[ithread], j, n, _bjmodn16[ithread]);
		MOD_ADD32(_bjmodn16[ithread], j, n, _bjmodn17[ithread]);
		MOD_ADD32(_bjmodn17[ithread], j, n, _bjmodn18[ithread]);
		MOD_ADD32(_bjmodn18[ithread], j, n, _bjmodn19[ithread]);
		MOD_ADD32(_bjmodn19[ithread], j, n, _bjmodn20[ithread]);
		MOD_ADD32(_bjmodn20[ithread], j, n, _bjmodn21[ithread]);
		MOD_ADD32(_bjmodn21[ithread], j, n, _bjmodn22[ithread]);
		MOD_ADD32(_bjmodn22[ithread], j, n, _bjmodn23[ithread]);

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
		ASSERT(HERE, tdat[ithread].s1p00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
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
		tdat[ithread].bjmodn12 = _bjmodn12[ithread];
		tdat[ithread].bjmodn13 = _bjmodn13[ithread];
		tdat[ithread].bjmodn14 = _bjmodn14[ithread];
		tdat[ithread].bjmodn15 = _bjmodn15[ithread];
		tdat[ithread].bjmodn16 = _bjmodn16[ithread];
		tdat[ithread].bjmodn17 = _bjmodn17[ithread];
		tdat[ithread].bjmodn18 = _bjmodn18[ithread];
		tdat[ithread].bjmodn19 = _bjmodn19[ithread];
		tdat[ithread].bjmodn20 = _bjmodn20[ithread];
		tdat[ithread].bjmodn21 = _bjmodn21[ithread];
		tdat[ithread].bjmodn22 = _bjmodn22[ithread];
		tdat[ithread].bjmodn23 = _bjmodn23[ithread];
		/* init carries	*/
		tdat[ithread].cy00 = _cy_00[ithread];
		tdat[ithread].cy01 = _cy_01[ithread];
		tdat[ithread].cy02 = _cy_02[ithread];
		tdat[ithread].cy03 = _cy_03[ithread];
		tdat[ithread].cy04 = _cy_04[ithread];
		tdat[ithread].cy05 = _cy_05[ithread];
		tdat[ithread].cy06 = _cy_06[ithread];
		tdat[ithread].cy07 = _cy_07[ithread];
		tdat[ithread].cy08 = _cy_08[ithread];
		tdat[ithread].cy09 = _cy_09[ithread];
		tdat[ithread].cy10 = _cy_10[ithread];
		tdat[ithread].cy11 = _cy_11[ithread];
		tdat[ithread].cy12 = _cy_12[ithread];
		tdat[ithread].cy13 = _cy_13[ithread];
		tdat[ithread].cy14 = _cy_14[ithread];
		tdat[ithread].cy15 = _cy_15[ithread];
		tdat[ithread].cy16 = _cy_16[ithread];
		tdat[ithread].cy17 = _cy_17[ithread];
		tdat[ithread].cy18 = _cy_18[ithread];
		tdat[ithread].cy19 = _cy_19[ithread];
		tdat[ithread].cy20 = _cy_20[ithread];
		tdat[ithread].cy21 = _cy_21[ithread];
		tdat[ithread].cy22 = _cy_22[ithread];
		tdat[ithread].cy23 = _cy_23[ithread];
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
		*bjmodn12 = _bjmodn12[ithread];
		*bjmodn13 = _bjmodn13[ithread];
		*bjmodn14 = _bjmodn14[ithread];
		*bjmodn15 = _bjmodn15[ithread];
		*bjmodn16 = _bjmodn16[ithread];
		*bjmodn17 = _bjmodn17[ithread];
		*bjmodn18 = _bjmodn18[ithread];
		*bjmodn19 = _bjmodn19[ithread];
		*bjmodn20 = _bjmodn20[ithread];
		*bjmodn21 = _bjmodn21[ithread];
		*bjmodn22 = _bjmodn22[ithread];
		*bjmodn23 = _bjmodn23[ithread];
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
		bjmodn12 = _bjmodn12[ithread];
		bjmodn13 = _bjmodn13[ithread];
		bjmodn14 = _bjmodn14[ithread];
		bjmodn15 = _bjmodn15[ithread];
		bjmodn16 = _bjmodn16[ithread];
		bjmodn17 = _bjmodn17[ithread];
		bjmodn18 = _bjmodn18[ithread];
		bjmodn19 = _bjmodn19[ithread];
		bjmodn20 = _bjmodn20[ithread];
		bjmodn21 = _bjmodn21[ithread];
		bjmodn22 = _bjmodn22[ithread];
		bjmodn23 = _bjmodn23[ithread];
	#endif
		/* init carries	*/
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		cy00->d0 = _cy_00[ithread];		cy00->d1 = _cy_01[ithread];		cy00->d2 = _cy_02[ithread];		cy00->d3 = _cy_03[ithread];
		cy04->d0 = _cy_04[ithread];		cy04->d1 = _cy_05[ithread];		cy04->d2 = _cy_06[ithread];		cy04->d3 = _cy_07[ithread];
		cy08->d0 = _cy_08[ithread];		cy08->d1 = _cy_09[ithread];		cy08->d2 = _cy_10[ithread];		cy08->d3 = _cy_11[ithread];
		cy12->d0 = _cy_12[ithread];		cy12->d1 = _cy_13[ithread];		cy12->d2 = _cy_14[ithread];		cy12->d3 = _cy_15[ithread];
		cy16->d0 = _cy_16[ithread];		cy16->d1 = _cy_17[ithread];		cy16->d2 = _cy_18[ithread];		cy16->d3 = _cy_19[ithread];
		cy20->d0 = _cy_20[ithread];		cy20->d1 = _cy_21[ithread];		cy20->d2 = _cy_22[ithread];		cy20->d3 = _cy_23[ithread];
	#elif defined(USE_SSE2)
		cy00->d0 = _cy_00[ithread];		cy00->d1 = _cy_01[ithread];
		cy02->d0 = _cy_02[ithread];		cy02->d1 = _cy_03[ithread];
		cy04->d0 = _cy_04[ithread];		cy04->d1 = _cy_05[ithread];
		cy06->d0 = _cy_06[ithread];		cy06->d1 = _cy_07[ithread];
		cy08->d0 = _cy_08[ithread];		cy08->d1 = _cy_09[ithread];
		cy10->d0 = _cy_10[ithread];		cy10->d1 = _cy_11[ithread];
		cy12->d0 = _cy_12[ithread];		cy12->d1 = _cy_13[ithread];
		cy14->d0 = _cy_14[ithread];		cy14->d1 = _cy_15[ithread];
		cy16->d0 = _cy_16[ithread];		cy16->d1 = _cy_17[ithread];
		cy18->d0 = _cy_18[ithread];		cy18->d1 = _cy_19[ithread];
		cy20->d0 = _cy_20[ithread];		cy20->d1 = _cy_21[ithread];
		cy22->d0 = _cy_22[ithread];		cy22->d1 = _cy_23[ithread];
	#else
		cy00 = _cy_00[ithread];
		cy01 = _cy_01[ithread];
		cy02 = _cy_02[ithread];
		cy03 = _cy_03[ithread];
		cy04 = _cy_04[ithread];
		cy05 = _cy_05[ithread];
		cy06 = _cy_06[ithread];
		cy07 = _cy_07[ithread];
		cy08 = _cy_08[ithread];
		cy09 = _cy_09[ithread];
		cy10 = _cy_10[ithread];
		cy11 = _cy_11[ithread];
		cy12 = _cy_12[ithread];
		cy13 = _cy_13[ithread];
		cy14 = _cy_14[ithread];
		cy15 = _cy_15[ithread];
		cy16 = _cy_16[ithread];
		cy17 = _cy_17[ithread];
		cy18 = _cy_18[ithread];
		cy19 = _cy_19[ithread];
		cy20 = _cy_20[ithread];
		cy21 = _cy_21[ithread];
		cy22 = _cy_22[ithread];
		cy23 = _cy_23[ithread];
	#endif

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix24_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_AVX
		_cy_00[ithread] = cy00->d0;	_cy_01[ithread] = cy00->d1;	_cy_02[ithread] = cy00->d2;	_cy_03[ithread] = cy00->d3;
		_cy_04[ithread] = cy04->d0;	_cy_05[ithread] = cy04->d1;	_cy_06[ithread] = cy04->d2;	_cy_07[ithread] = cy04->d3;
		_cy_08[ithread] = cy08->d0;	_cy_09[ithread] = cy08->d1;	_cy_10[ithread] = cy08->d2;	_cy_11[ithread] = cy08->d3;
		_cy_12[ithread] = cy12->d0;	_cy_13[ithread] = cy12->d1;	_cy_14[ithread] = cy12->d2;	_cy_15[ithread] = cy12->d3;
		_cy_16[ithread] = cy16->d0;	_cy_17[ithread] = cy16->d1;	_cy_18[ithread] = cy16->d2;	_cy_19[ithread] = cy16->d3;
		_cy_20[ithread] = cy20->d0;	_cy_21[ithread] = cy20->d1;	_cy_22[ithread] = cy20->d2;	_cy_23[ithread] = cy20->d3;
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		_cy_00[ithread] = cy00->d0;	_cy_01[ithread] = cy00->d1;
		_cy_02[ithread] = cy02->d0;	_cy_03[ithread] = cy02->d1;
		_cy_04[ithread] = cy04->d0;	_cy_05[ithread] = cy04->d1;
		_cy_06[ithread] = cy06->d0;	_cy_07[ithread] = cy06->d1;
		_cy_08[ithread] = cy08->d0;	_cy_09[ithread] = cy08->d1;
		_cy_10[ithread] = cy10->d0;	_cy_11[ithread] = cy10->d1;
		_cy_12[ithread] = cy12->d0;	_cy_13[ithread] = cy12->d1;
		_cy_14[ithread] = cy14->d0;	_cy_15[ithread] = cy14->d1;
		_cy_16[ithread] = cy16->d0;	_cy_17[ithread] = cy16->d1;
		_cy_18[ithread] = cy18->d0;	_cy_19[ithread] = cy18->d1;
		_cy_20[ithread] = cy20->d0;	_cy_21[ithread] = cy20->d1;
		_cy_22[ithread] = cy22->d0;	_cy_23[ithread] = cy22->d1;
		maxerr = MAX(max_err->d0,max_err->d1);
	#else
		_cy_00[ithread] = cy00;
		_cy_01[ithread] = cy01;
		_cy_02[ithread] = cy02;
		_cy_03[ithread] = cy03;
		_cy_04[ithread] = cy04;
		_cy_05[ithread] = cy05;
		_cy_06[ithread] = cy06;
		_cy_07[ithread] = cy07;
		_cy_08[ithread] = cy08;
		_cy_09[ithread] = cy09;
		_cy_10[ithread] = cy10;
		_cy_11[ithread] = cy11;
		_cy_12[ithread] = cy12;
		_cy_13[ithread] = cy13;
		_cy_14[ithread] = cy14;
		_cy_15[ithread] = cy15;
		_cy_16[ithread] = cy16;
		_cy_17[ithread] = cy17;
		_cy_18[ithread] = cy18;
		_cy_19[ithread] = cy19;
		_cy_20[ithread] = cy20;
		_cy_21[ithread] = cy21;
		_cy_22[ithread] = cy22;
		_cy_23[ithread] = cy23;
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
		ASSERT(HERE, 0x0 == cy24_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix32_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}

		_cy_00[ithread] = tdat[ithread].cy00;
		_cy_01[ithread] = tdat[ithread].cy01;
		_cy_02[ithread] = tdat[ithread].cy02;
		_cy_03[ithread] = tdat[ithread].cy03;
		_cy_04[ithread] = tdat[ithread].cy04;
		_cy_05[ithread] = tdat[ithread].cy05;
		_cy_06[ithread] = tdat[ithread].cy06;
		_cy_07[ithread] = tdat[ithread].cy07;
		_cy_08[ithread] = tdat[ithread].cy08;
		_cy_09[ithread] = tdat[ithread].cy09;
		_cy_10[ithread] = tdat[ithread].cy10;
		_cy_11[ithread] = tdat[ithread].cy11;
		_cy_12[ithread] = tdat[ithread].cy12;
		_cy_13[ithread] = tdat[ithread].cy13;
		_cy_14[ithread] = tdat[ithread].cy14;
		_cy_15[ithread] = tdat[ithread].cy15;
		_cy_16[ithread] = tdat[ithread].cy16;
		_cy_17[ithread] = tdat[ithread].cy17;
		_cy_18[ithread] = tdat[ithread].cy18;
		_cy_19[ithread] = tdat[ithread].cy19;
		_cy_20[ithread] = tdat[ithread].cy20;
		_cy_21[ithread] = tdat[ithread].cy21;
		_cy_22[ithread] = tdat[ithread].cy22;
		_cy_23[ithread] = tdat[ithread].cy23;
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-24 forward DIF FFT of the first block of 24 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 24 outputs of (1);
	!   (3) Reweight and perform a radix-24 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 24 elements and repeat (1-4).
	*/
	t00= _cy_00[CY_THREADS - 1];
	t01= _cy_01[CY_THREADS - 1];
	t02= _cy_02[CY_THREADS - 1];
	t03= _cy_03[CY_THREADS - 1];
	t04= _cy_04[CY_THREADS - 1];
	t05= _cy_05[CY_THREADS - 1];
	t06= _cy_06[CY_THREADS - 1];
	t07= _cy_07[CY_THREADS - 1];
	t08= _cy_08[CY_THREADS - 1];
	t09= _cy_09[CY_THREADS - 1];
	t10= _cy_10[CY_THREADS - 1];
	t11= _cy_11[CY_THREADS - 1];
	t12= _cy_12[CY_THREADS - 1];
	t13= _cy_13[CY_THREADS - 1];
	t14= _cy_14[CY_THREADS - 1];
	t15= _cy_15[CY_THREADS - 1];
	t16= _cy_16[CY_THREADS - 1];
	t17= _cy_17[CY_THREADS - 1];
	t18= _cy_18[CY_THREADS - 1];
	t19= _cy_19[CY_THREADS - 1];
	t20= _cy_20[CY_THREADS - 1];
	t21= _cy_21[CY_THREADS - 1];
	t22= _cy_22[CY_THREADS - 1];
	t23= _cy_23[CY_THREADS - 1];

	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		ASSERT(HERE, CY_THREADS > 1,"radix24_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
		_cy_00[ithread] = _cy_00[ithread-1];
		_cy_01[ithread] = _cy_01[ithread-1];
		_cy_02[ithread] = _cy_02[ithread-1];
		_cy_03[ithread] = _cy_03[ithread-1];
		_cy_04[ithread] = _cy_04[ithread-1];
		_cy_05[ithread] = _cy_05[ithread-1];
		_cy_06[ithread] = _cy_06[ithread-1];
		_cy_07[ithread] = _cy_07[ithread-1];
		_cy_08[ithread] = _cy_08[ithread-1];
		_cy_09[ithread] = _cy_09[ithread-1];
		_cy_10[ithread] = _cy_10[ithread-1];
		_cy_11[ithread] = _cy_11[ithread-1];
		_cy_12[ithread] = _cy_12[ithread-1];
		_cy_13[ithread] = _cy_13[ithread-1];
		_cy_14[ithread] = _cy_14[ithread-1];
		_cy_15[ithread] = _cy_15[ithread-1];
		_cy_16[ithread] = _cy_16[ithread-1];
		_cy_17[ithread] = _cy_17[ithread-1];
		_cy_18[ithread] = _cy_18[ithread-1];
		_cy_19[ithread] = _cy_19[ithread-1];
		_cy_20[ithread] = _cy_20[ithread-1];
		_cy_21[ithread] = _cy_21[ithread-1];
		_cy_22[ithread] = _cy_22[ithread-1];
		_cy_23[ithread] = _cy_23[ithread-1];
	}

	_cy_00[0] =+t23;	/* ...The wraparound carry is here: */
	_cy_01[0] = t00;
	_cy_02[0] = t01;
	_cy_03[0] = t02;
	_cy_04[0] = t03;
	_cy_05[0] = t04;
	_cy_06[0] = t05;
	_cy_07[0] = t06;
	_cy_08[0] = t07;
	_cy_09[0] = t08;
	_cy_10[0] = t09;
	_cy_11[0] = t10;
	_cy_12[0] = t11;
	_cy_13[0] = t12;
	_cy_14[0] = t13;
	_cy_15[0] = t14;
	_cy_16[0] = t15;
	_cy_17[0] = t16;
	_cy_18[0] = t17;
	_cy_19[0] = t18;
	_cy_20[0] = t19;
	_cy_21[0] = t20;
	_cy_22[0] = t21;
	_cy_23[0] = t22;

	full_pass = 0;
	scale = 1;

	jhi = 7;

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + jhi; j++)
		{
			// Generate padded version of j, since prepadding pini is thread-count unsafe:
			j2 = j + ( (j >> DAT_BITS) << PAD_BITS );
			j1 = j2;
			a[j1    ] *= radix_inv;
			a[j1+p01] *= radix_inv;
			a[j1+p02] *= radix_inv;
			a[j1+p03] *= radix_inv;
			a[j1+p04] *= radix_inv;
			a[j1+p05] *= radix_inv;
			a[j1+p06] *= radix_inv;
			a[j1+p07] *= radix_inv;
			j1 = j2 + p08;
			a[j1    ] *= radix_inv;
			a[j1+p01] *= radix_inv;
			a[j1+p02] *= radix_inv;
			a[j1+p03] *= radix_inv;
			a[j1+p04] *= radix_inv;
			a[j1+p05] *= radix_inv;
			a[j1+p06] *= radix_inv;
			a[j1+p07] *= radix_inv;
			j1 = j2 + p16;
			a[j1    ] *= radix_inv;
			a[j1+p01] *= radix_inv;
			a[j1+p02] *= radix_inv;
			a[j1+p03] *= radix_inv;
			a[j1+p04] *= radix_inv;
			a[j1+p05] *= radix_inv;
			a[j1+p06] *= radix_inv;
			a[j1+p07] *= radix_inv;
		}
	}
}	/* endfor(outer) */

	t00 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00 += fabs(_cy_00[0])+fabs(_cy_01[0])+fabs(_cy_02[0])+fabs(_cy_03[0])+fabs(_cy_04[0])+fabs(_cy_05[0])+fabs(_cy_06[0])+fabs(_cy_07[0])+fabs(_cy_08[0])+fabs(_cy_09[0])+fabs(_cy_10[0])+fabs(_cy_11[0])+fabs(_cy_12[0])+fabs(_cy_13[0])+fabs(_cy_14[0])+fabs(_cy_15[0])+fabs(_cy_16[0])+fabs(_cy_17[0])+fabs(_cy_18[0])+fabs(_cy_19[0]+_cy_20[0])+fabs(_cy_21[0])+fabs(_cy_22[0])+fabs(_cy_23[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}

	if(t00 != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix24_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

int radix24_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter                 , uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-24 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-24 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int NDIVR, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double radix_inv, n2inv;
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i
	,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23,temp,scale;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "radix24_ditN_cy_dif1_nochk: Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/24;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/24 in radix24_ditN_cy_dif1.\n",iter);
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
	  radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)24));
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

	  bjmodnini=0;
	  for(j=0; j < NDIVR; j++)
	  {
	    bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
	  }
	}

/*...The radix-24 final DIT pass is here.	*/

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

	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy00= -2;
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

	col=0;
	co2=(n >> nwt_bits)-1+24;
	co3=co2-24;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

/*...The radix-24 DIT pass is here:	*/
/*...gather the needed data (24 64-bit complex, i.e. 48 64-bit reals) and do 3 radix-8 transforms,	*/
	                 /*                                                                                 inputs                                                                                    */ /*                 intermediates                    */ /*                                                         outputs                                                           */
	RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,rt,it);
	RADIX_08_DIT(a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,rt,it);
	RADIX_08_DIT(a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],a[j1+p23],a[j2+p23],a[j1+p22],a[j2+p22],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,rt,it);

/*...and now do 8 in-place radix-3 transforms.	*/
	                 /*                  inputs                   */ /*             intermediates                 */ /*                       outputs                               */
	RADIX_03_DFT(s,c3m1,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,t1,t2,t3,t4,t5,t6,a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i);
	RADIX_03_DFT(s,c3m1,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,t1,t2,t3,t4,t5,t6,a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i);
	RADIX_03_DFT(s,c3m1,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,t1,t2,t3,t4,t5,t6,a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i);
	RADIX_03_DFT(s,c3m1,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,t1,t2,t3,t4,t5,t6,a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i);
	RADIX_03_DFT(s,c3m1,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,t1,t2,t3,t4,t5,t6,a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i);
	RADIX_03_DFT(s,c3m1,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,t1,t2,t3,t4,t5,t6,a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i);
	RADIX_03_DFT(s,c3m1,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,t1,t2,t3,t4,t5,t6,a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i);
	RADIX_03_DFT(s,c3m1,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,t1,t2,t3,t4,t5,t6,a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i);

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 24 separate blocks of the A-array, we need 24 separate carries.	*/

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

	    i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
	    co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			   and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif

/*...The radix-24 DIF pass is here:	*/

/* EWM: 10/21/04: We swap the following outputs of the radix-3 transforms: {1,9,17}<=>{5,13,21}, {3,11,19}<=>{7,15,23}, so that the indexing
	              winds up being in-place. This allows us to properly re-use the ajp1 variables in the carry-pass version of this routine.
*/

/*...gather the needed data (24 64-bit complex, i.e. 48 64-bit reals) and do 8 in-place radix-3 transforms...*/
	                 /*                        inputs                               */ /*             intermediates                 */ /*                 outputs                   */
#if PFETCH
	RADIX_03_DFT_PFETCH(s,c3m1,a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,t1,t2,t3,t4,t5,t6,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,p01);
	RADIX_03_DFT_PFETCH(s,c3m1,a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,t1,t2,t3,t4,t5,t6,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,p02);
	RADIX_03_DFT_PFETCH(s,c3m1,a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,t1,t2,t3,t4,t5,t6,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,p03);
	RADIX_03_DFT_PFETCH(s,c3m1,a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,t1,t2,t3,t4,t5,t6,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,p04);
	RADIX_03_DFT_PFETCH(s,c3m1,a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,t1,t2,t3,t4,t5,t6,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,p05);
	RADIX_03_DFT_PFETCH(s,c3m1,a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,t1,t2,t3,t4,t5,t6,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,p06);
	RADIX_03_DFT_PFETCH(s,c3m1,a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,t1,t2,t3,t4,t5,t6,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,p07);
	RADIX_03_DFT_PFETCH(s,c3m1,a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,p08);
#else
	RADIX_03_DFT       (s,c3m1,a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,t1,t2,t3,t4,t5,t6,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i);
	RADIX_03_DFT       (s,c3m1,a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,t1,t2,t3,t4,t5,t6,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i);
	RADIX_03_DFT       (s,c3m1,a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,t1,t2,t3,t4,t5,t6,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i);
	RADIX_03_DFT       (s,c3m1,a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,t1,t2,t3,t4,t5,t6,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i);
	RADIX_03_DFT       (s,c3m1,a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,t1,t2,t3,t4,t5,t6,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i);
	RADIX_03_DFT       (s,c3m1,a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,t1,t2,t3,t4,t5,t6,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i);
	RADIX_03_DFT       (s,c3m1,a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,t1,t2,t3,t4,t5,t6,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i);
	RADIX_03_DFT       (s,c3m1,a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i);
#endif

/*...and now do 3 radix-8 transforms:	*/
	                 /*                                                          inputs                                                           */ /*                 intermediates                    */ /*                 outputs                   */
#if PFETCH
	RADIX_08_DIF_PFETCH(a1p00r,a1p00i,a1p05r,a1p05i,a1p02r,a1p02i,a1p07r,a1p07i,a1p04r,a1p04i,a1p01r,a1p01i,a1p06r,a1p06i,a1p03r,a1p03i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],rt,it,p09,p10,p11,p12,p13);
	RADIX_08_DIF_PFETCH(a1p08r,a1p08i,a1p13r,a1p13i,a1p10r,a1p10i,a1p15r,a1p15i,a1p12r,a1p12i,a1p09r,a1p09i,a1p14r,a1p14i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],a[j1+p17],a[j2+p17],a[j1+p16],a[j2+p16],a[j1+p23],a[j2+p23],a[j1+p22],a[j2+p22],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],rt,it,p14,p15,p16,p17,p18);
	RADIX_08_DIF_PFETCH(a1p16r,a1p16i,a1p21r,a1p21i,a1p18r,a1p18i,a1p23r,a1p23i,a1p20r,a1p20i,a1p17r,a1p17i,a1p22r,a1p22i,a1p19r,a1p19i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],rt,it,p19,p20,p21,p22,p23);
#else
	RADIX_08_DIF       (a1p00r,a1p00i,a1p05r,a1p05i,a1p02r,a1p02i,a1p07r,a1p07i,a1p04r,a1p04i,a1p01r,a1p01i,a1p06r,a1p06i,a1p03r,a1p03i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],rt,it);
	RADIX_08_DIF       (a1p08r,a1p08i,a1p13r,a1p13i,a1p10r,a1p10i,a1p15r,a1p15i,a1p12r,a1p12i,a1p09r,a1p09i,a1p14r,a1p14i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],a[j1+p17],a[j2+p17],a[j1+p16],a[j2+p16],a[j1+p23],a[j2+p23],a[j1+p22],a[j2+p22],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],rt,it);
	RADIX_08_DIF       (a1p16r,a1p16i,a1p21r,a1p21i,a1p18r,a1p18i,a1p23r,a1p23i,a1p20r,a1p20i,a1p17r,a1p17i,a1p22r,a1p22i,a1p19r,a1p19i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],rt,it);
#endif
	  iroot += root_incr;		/* increment sincos index.	*/

	  }

	  jstart += nwt;
	  jhi    += nwt;
	  col += 24;
	  co3 -= 24;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-24 forward DIF FFT of the first block of 24 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 24 outputs of (1);
!   (3) Reweight and perform a radix-24 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 24 elements and repeat (1-4).
*/
	t1  = cy23;
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
	cy00= t1;

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
	}
}

	if(fabs(cy00)+fabs(cy01)+fabs(cy02)+fabs(cy03)+fabs(cy04)+fabs(cy05)+fabs(cy06)+fabs(cy07)+fabs(cy08)+fabs(cy09)+fabs(cy10)+fabs(cy11)+fabs(cy12)+fabs(cy13)+fabs(cy14)+fabs(cy15)+fabs(cy16)+fabs(cy17)+fabs(cy18)+fabs(cy19)+fabs(cy20)+fabs(cy21)+fabs(cy22)+fabs(cy23) != 0.0)
	{
	    sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix24_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

void radix24_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-24 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix3,8_dif_pass for details on the radix-3,8 subtransforms.
*/
	int j,j1,j2,jt,jp;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p16, first_entry=TRUE;
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
	double *add0, *addr;
	double rt,it
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i;

	if(!first_entry && (n/24) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/24;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-24 pass is here.	*/

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

		add0 = &a[j1];
	/*
	Twiddleless version arranges 8 sets of radix-3 DFT inputs as follows: 0 in upper left corner, decrement 8 horizontally and 3 vertically:

		RADIX_03_DFT(00,16,08)
		RADIX_03_DFT(21,13,05)
		RADIX_03_DFT(18,10,02)
		RADIX_03_DFT(15,07,23)
		RADIX_03_DFT(12,04,20)
		RADIX_03_DFT(09,01,17)
		RADIX_03_DFT(06,22,14)
		RADIX_03_DFT(03,19,11)

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-8 DFT outputs.
	*/

	/* EWM: 10/21/04: We swap the following outputs of the radix-3 transforms: {1,9,17}<=>{5,13,21}, {3,11,19}<=>{7,15,23}, so that the indexing
					  winds up being in-place. This allows us to properly re-use the ajp1 variables in the carry-pass version of this routine.
	*/
	/*...gather the needed data (24 64-bit complex) and do 8 radix-3 transforms...*/
					 /*                        inputs                         */ /*             intermediates                 */ /*                 outputs                   */
		RADIX_03_DFT(s,c3m1,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08],t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i);	jt = j1+p05; jp = j2+p05;
		RADIX_03_DFT(s,c3m1,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],t01,t02,t03,t04,t05,t06,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i);	jt = j1+p02; jp = j2+p02;
		RADIX_03_DFT(s,c3m1,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],t01,t02,t03,t04,t05,t06,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i);	jt = j1+p07; jp = j2+p07;
		RADIX_03_DFT(s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i);	jt = j1+p04; jp = j2+p04;
		RADIX_03_DFT(s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i);	jt = j1+p01; jp = j2+p01;
		RADIX_03_DFT(s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i);	jt = j1+p06; jp = j2+p06;
		RADIX_03_DFT(s,c3m1,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i);	jt = j1+p03; jp = j2+p03;
		RADIX_03_DFT(s,c3m1,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i);

	/*...and now do 3 radix-8 transforms:	*/
					 /*                                                          inputs                                                           */ /*                 intermediates                    */ /*                 outputs                   */
		RADIX_08_DIF_PFETCH(a1p00r,a1p00i,a1p05r,a1p05i,a1p02r,a1p02i,a1p07r,a1p07i,a1p04r,a1p04i,a1p01r,a1p01i,a1p06r,a1p06i,a1p03r,a1p03i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],rt,it,p08+p01,p08+p02,p08+p03,p08+p04,p08+p05);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIF_PFETCH(a1p08r,a1p08i,a1p13r,a1p13i,a1p10r,a1p10i,a1p15r,a1p15i,a1p12r,a1p12i,a1p09r,a1p09i,a1p14r,a1p14i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it,p08+p06,p08+p07,p16    ,p16+p01,p16+p02);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF_PFETCH(a1p16r,a1p16i,a1p21r,a1p21i,a1p18r,a1p18i,a1p23r,a1p23i,a1p20r,a1p20i,a1p17r,a1p17i,a1p22r,a1p22i,a1p19r,a1p19i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,p16+p03,p16+p04,p16+p05,p16+p06,p16+p07);
	}
}

/***************/

void radix24_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Time
!
!...Subroutine to perform an initial radix-24 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2,jt,jp;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p16, first_entry=TRUE;
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;	/* sin(twopi/3)	... Same as for DIF, since we use index permutations to effect the desired inverse-DFT */
	double rt,it
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i;

	if(!first_entry && (n/24) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/24;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-24 pass is here.	*/

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

		00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23

		00,16,08,21,13,05,18,10,02,15,07,23,12,04,20,09,01,17,06,22,14,03,19,11.	(*)

	Remember, inputs to DIT are bit-reversed, so

		a[00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23] contain [cf. bit-reversal index array output by test_fft_radix()]
		x[00,12,06,18,03,15,09,21,01,13,07,19,04,16,10,22,02,14,08,20,05,17,11,23], which get swapped [using the permutation (*)] to
		x[00,12,18,06,21,09,15,03,16,04,10,22,13,01,07,19,08,20,02,14,05,17,23,11], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04,13,12,14,15,09,08,10,11,18,19,16,17,20,21,23,22]. These are the 3 octets going into the radix-8 DFTs.

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF] to properly permute the radix-3 DFT outputs.
	*/
	/*...gather the needed data (24 64-bit complex) and do 3 radix-8 transforms,	*/

		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,rt,it);

	/*...and now do 8 radix-3 transforms.	*/

		RADIX_03_DFT(s,c3m1,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,t01,t02,t03,t04,t05,t06,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08]);	jt = j1+p01; jp = j2+p01;
		RADIX_03_DFT(s,c3m1,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,t01,t02,t03,t04,t05,t06,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);	jt = j1+p02; jp = j2+p02;
		RADIX_03_DFT(s,c3m1,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,t01,t02,t03,t04,t05,t06,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ]);	jt = j1+p03; jp = j2+p03;
		RADIX_03_DFT(s,c3m1,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,t01,t02,t03,t04,t05,t06,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08]);	jt = j1+p04; jp = j2+p04;
		RADIX_03_DFT(s,c3m1,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,t01,t02,t03,t04,t05,t06,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);	jt = j1+p05; jp = j2+p05;
		RADIX_03_DFT(s,c3m1,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,t01,t02,t03,t04,t05,t06,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ]);	jt = j1+p06; jp = j2+p06;
		RADIX_03_DFT(s,c3m1,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,t01,t02,t03,t04,t05,t06,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08]);	jt = j1+p07; jp = j2+p07;
		RADIX_03_DFT(s,c3m1,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,t01,t02,t03,t04,t05,t06,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy24_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p05,p06,p07,p08,p16;
		int j,j1,j2,k,l;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif

	#ifdef USE_SSE2

		double *add0, *add1, *add2, *add3;
		const double crnd = 3.0*0x4000000*0x2000000;
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *isrt2, *cc0, *cc3, *max_err, *sse2_rnd, *half_arr, *tmp, *s1p00r,*s1p04r,*s1p08r,*s1p12r,*s1p16r,*s1p20r;
		vec_dbl *cy00,*cy04,*cy08,*cy12,*cy16,*cy20;
	  #ifndef USE_AVX
		vec_dbl *cy02,*cy06,*cy10,*cy14,*cy18,*cy22;
	  #endif
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23;
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
						s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
		double *base, *baseinv;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		double rt,it,temp,frac
			,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15
			,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r
			,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i
			,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23;
		int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23;

	#endif

		struct cy_thread_data_t* thread_arg = targ;
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
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );

	#ifdef USE_SSE2
		s1p00r	= thread_arg->s1p00r;
		s1p04r	= s1p00r + 0x08;
		s1p08r	= s1p00r + 0x10;
		s1p12r	= s1p00r + 0x18;
		s1p16r	= s1p00r + 0x20;
		s1p20r	= s1p00r + 0x28;
		isrt2	= s1p00r + 0x30;
		cc3		= s1p00r + 0x31;
		cc0  	= s1p00r + 0x32;
	  #ifdef USE_AVX
		isrt2	= s1p00r + 0x30;
		cc3		= s1p00r + 0x31;
		cc0  	= s1p00r + 0x32;
		cy00	= s1p00r + 0x33;
		cy04	= s1p00r + 0x34;
		cy08	= s1p00r + 0x35;
		cy12	= s1p00r + 0x36;
		cy16	= s1p00r + 0x37;
		cy20	= s1p00r + 0x38;
		max_err = s1p00r + 0x39;
		sse2_rnd= s1p00r + 0x3a;
		half_arr= s1p00r + 0x3b;
	  #else
		cy00	= s1p00r + 0x33;
		cy02	= s1p00r + 0x34;
		cy04	= s1p00r + 0x35;
		cy06	= s1p00r + 0x36;
		cy08	= s1p00r + 0x37;
		cy10	= s1p00r + 0x38;
		cy12	= s1p00r + 0x39;
		cy14	= s1p00r + 0x3a;
		cy16	= s1p00r + 0x3b;
		cy18	= s1p00r + 0x3c;
		cy20	= s1p00r + 0x3d;
		cy22	= s1p00r + 0x3e;
		max_err = s1p00r + 0x3f;
		sse2_rnd= s1p00r + 0x40;
		half_arr= s1p00r + 0x41;	/* This table needs 20x16 bytes */
	  #endif
		ASSERT(HERE, (s1p00r == thread_arg->s1p00r), "thread-local memcheck failed!");
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

		sign_mask = (uint64*)(s1p00r + radix24_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
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
		bjmodn12 = bjmodn00 + 12;
		bjmodn13 = bjmodn00 + 13;
		bjmodn14 = bjmodn00 + 14;
		bjmodn15 = bjmodn00 + 15;
		bjmodn16 = bjmodn00 + 16;
		bjmodn17 = bjmodn00 + 17;
		bjmodn18 = bjmodn00 + 18;
		bjmodn19 = bjmodn00 + 19;
		bjmodn20 = bjmodn00 + 20;
		bjmodn21 = bjmodn00 + 21;
		bjmodn22 = bjmodn00 + 22;
		bjmodn23 = bjmodn00 + 23;
	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->s1p00r  ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		/* Init DWT-indices: */	/* init carries	*/
	#ifdef USE_AVX
		cy00->d0 = thread_arg->cy00;		*bjmodn00 = thread_arg->bjmodn00;
		cy00->d1 = thread_arg->cy01;		*bjmodn01 = thread_arg->bjmodn01;
		cy00->d2 = thread_arg->cy02;		*bjmodn02 = thread_arg->bjmodn02;
		cy00->d3 = thread_arg->cy03;		*bjmodn03 = thread_arg->bjmodn03;
		cy04->d0 = thread_arg->cy04;		*bjmodn04 = thread_arg->bjmodn04;
		cy04->d1 = thread_arg->cy05;		*bjmodn05 = thread_arg->bjmodn05;
		cy04->d2 = thread_arg->cy06;		*bjmodn06 = thread_arg->bjmodn06;
		cy04->d3 = thread_arg->cy07;		*bjmodn07 = thread_arg->bjmodn07;
		cy08->d0 = thread_arg->cy08;		*bjmodn08 = thread_arg->bjmodn08;
		cy08->d1 = thread_arg->cy09;		*bjmodn09 = thread_arg->bjmodn09;
		cy08->d2 = thread_arg->cy10;		*bjmodn10 = thread_arg->bjmodn10;
		cy08->d3 = thread_arg->cy11;		*bjmodn11 = thread_arg->bjmodn11;
		cy12->d0 = thread_arg->cy12;		*bjmodn12 = thread_arg->bjmodn12;
		cy12->d1 = thread_arg->cy13;		*bjmodn13 = thread_arg->bjmodn13;
		cy12->d2 = thread_arg->cy14;		*bjmodn14 = thread_arg->bjmodn14;
		cy12->d3 = thread_arg->cy15;		*bjmodn15 = thread_arg->bjmodn15;
		cy16->d0 = thread_arg->cy16;		*bjmodn16 = thread_arg->bjmodn16;
		cy16->d1 = thread_arg->cy17;		*bjmodn17 = thread_arg->bjmodn17;
		cy16->d2 = thread_arg->cy18;		*bjmodn18 = thread_arg->bjmodn18;
		cy16->d3 = thread_arg->cy19;		*bjmodn19 = thread_arg->bjmodn19;
		cy20->d0 = thread_arg->cy20;		*bjmodn20 = thread_arg->bjmodn20;
		cy20->d1 = thread_arg->cy21;		*bjmodn21 = thread_arg->bjmodn21;
		cy20->d2 = thread_arg->cy22;		*bjmodn22 = thread_arg->bjmodn22;
		cy20->d3 = thread_arg->cy23;		*bjmodn23 = thread_arg->bjmodn23;

	#elif defined(USE_SSE2)

		cy00->d0 = thread_arg->cy00;		*bjmodn00 = thread_arg->bjmodn00;
		cy00->d1 = thread_arg->cy01;		*bjmodn01 = thread_arg->bjmodn01;
		cy02->d0 = thread_arg->cy02;		*bjmodn02 = thread_arg->bjmodn02;
		cy02->d1 = thread_arg->cy03;		*bjmodn03 = thread_arg->bjmodn03;
		cy04->d0 = thread_arg->cy04;		*bjmodn04 = thread_arg->bjmodn04;
		cy04->d1 = thread_arg->cy05;		*bjmodn05 = thread_arg->bjmodn05;
		cy06->d0 = thread_arg->cy06;		*bjmodn06 = thread_arg->bjmodn06;
		cy06->d1 = thread_arg->cy07;		*bjmodn07 = thread_arg->bjmodn07;
		cy08->d0 = thread_arg->cy08;		*bjmodn08 = thread_arg->bjmodn08;
		cy08->d1 = thread_arg->cy09;		*bjmodn09 = thread_arg->bjmodn09;
		cy10->d0 = thread_arg->cy10;		*bjmodn10 = thread_arg->bjmodn10;
		cy10->d1 = thread_arg->cy11;		*bjmodn11 = thread_arg->bjmodn11;
		cy12->d0 = thread_arg->cy12;		*bjmodn12 = thread_arg->bjmodn12;
		cy12->d1 = thread_arg->cy13;		*bjmodn13 = thread_arg->bjmodn13;
		cy14->d0 = thread_arg->cy14;		*bjmodn14 = thread_arg->bjmodn14;
		cy14->d1 = thread_arg->cy15;		*bjmodn15 = thread_arg->bjmodn15;
		cy16->d0 = thread_arg->cy16;		*bjmodn16 = thread_arg->bjmodn16;
		cy16->d1 = thread_arg->cy17;		*bjmodn17 = thread_arg->bjmodn17;
		cy18->d0 = thread_arg->cy18;		*bjmodn18 = thread_arg->bjmodn18;
		cy18->d1 = thread_arg->cy19;		*bjmodn19 = thread_arg->bjmodn19;
		cy20->d0 = thread_arg->cy20;		*bjmodn20 = thread_arg->bjmodn20;
		cy20->d1 = thread_arg->cy21;		*bjmodn21 = thread_arg->bjmodn21;
		cy22->d0 = thread_arg->cy22;		*bjmodn22 = thread_arg->bjmodn22;
		cy22->d1 = thread_arg->cy23;		*bjmodn23 = thread_arg->bjmodn23;

	#else

		cy00 = thread_arg->cy00;		bjmodn00 = thread_arg->bjmodn00;
		cy01 = thread_arg->cy01;		bjmodn01 = thread_arg->bjmodn01;
		cy02 = thread_arg->cy02;		bjmodn02 = thread_arg->bjmodn02;
		cy03 = thread_arg->cy03;		bjmodn03 = thread_arg->bjmodn03;
		cy04 = thread_arg->cy04;		bjmodn04 = thread_arg->bjmodn04;
		cy05 = thread_arg->cy05;		bjmodn05 = thread_arg->bjmodn05;
		cy06 = thread_arg->cy06;		bjmodn06 = thread_arg->bjmodn06;
		cy07 = thread_arg->cy07;		bjmodn07 = thread_arg->bjmodn07;
		cy08 = thread_arg->cy08;		bjmodn08 = thread_arg->bjmodn08;
		cy09 = thread_arg->cy09;		bjmodn09 = thread_arg->bjmodn09;
		cy10 = thread_arg->cy10;		bjmodn10 = thread_arg->bjmodn10;
		cy11 = thread_arg->cy11;		bjmodn11 = thread_arg->bjmodn11;
		cy12 = thread_arg->cy12;		bjmodn12 = thread_arg->bjmodn12;
		cy13 = thread_arg->cy13;		bjmodn13 = thread_arg->bjmodn13;
		cy14 = thread_arg->cy14;		bjmodn14 = thread_arg->bjmodn14;
		cy15 = thread_arg->cy15;		bjmodn15 = thread_arg->bjmodn15;
		cy16 = thread_arg->cy16;		bjmodn16 = thread_arg->bjmodn16;
		cy17 = thread_arg->cy17;		bjmodn17 = thread_arg->bjmodn17;
		cy18 = thread_arg->cy18;		bjmodn18 = thread_arg->bjmodn18;
		cy19 = thread_arg->cy19;		bjmodn19 = thread_arg->bjmodn19;
		cy20 = thread_arg->cy20;		bjmodn20 = thread_arg->bjmodn20;
		cy21 = thread_arg->cy21;		bjmodn21 = thread_arg->bjmodn21;
		cy22 = thread_arg->cy22;		bjmodn22 = thread_arg->bjmodn22;
		cy23 = thread_arg->cy23;		bjmodn23 = thread_arg->bjmodn23;

	#endif	// SSE2 or AVX?

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix24_main_carry_loop.h"

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
		thread_arg->cy12 = cy12->d0;
		thread_arg->cy13 = cy12->d1;
		thread_arg->cy14 = cy12->d2;
		thread_arg->cy15 = cy12->d3;
		thread_arg->cy16 = cy16->d0;
		thread_arg->cy17 = cy16->d1;
		thread_arg->cy18 = cy16->d2;
		thread_arg->cy19 = cy16->d3;
		thread_arg->cy20 = cy20->d0;
		thread_arg->cy21 = cy20->d1;
		thread_arg->cy22 = cy20->d2;
		thread_arg->cy23 = cy20->d3;
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
		thread_arg->cy12 = cy12->d0;
		thread_arg->cy13 = cy12->d1;
		thread_arg->cy14 = cy14->d0;
		thread_arg->cy15 = cy14->d1;
		thread_arg->cy16 = cy16->d0;
		thread_arg->cy17 = cy16->d1;
		thread_arg->cy18 = cy18->d0;
		thread_arg->cy19 = cy18->d1;
		thread_arg->cy20 = cy20->d0;
		thread_arg->cy21 = cy20->d1;
		thread_arg->cy22 = cy22->d0;
		thread_arg->cy23 = cy22->d1;
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
		thread_arg->cy12 = cy12;
		thread_arg->cy13 = cy13;
		thread_arg->cy14 = cy14;
		thread_arg->cy15 = cy15;
		thread_arg->cy16 = cy16;
		thread_arg->cy17 = cy17;
		thread_arg->cy18 = cy18;
		thread_arg->cy19 = cy19;
		thread_arg->cy20 = cy20;
		thread_arg->cy21 = cy21;
		thread_arg->cy22 = cy22;
		thread_arg->cy23 = cy23;

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
