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

#include "Mlucas.h"

#define USE_SCALAR_DFT_MACRO	1

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
  // For Fermat-mod we use RADIX*4 = 64 [note there is no LOACC option for this power-of-2 DFT] more
  // slots in AVX mode for the compact negacyclic-roots chained-multiply scheme. Add larger of the 2 numbers -
  // 68 for AVX, 20 for SSE2 - to (half_arr_offset16 + RADIX) to get AVX value of radix16_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset16 = 46;	// + RADIX = 62; Used for thread local-storage-integrity checking
	const int radix16_creals_in_local_store = 132;	// (half_arr_offset16 + RADIX) + 68 and round up to nearest multiple of 4
  #else
	const int half_arr_offset16 = 54;	// + RADIX = 70; Used for thread local-storage-integrity checking
	const int radix16_creals_in_local_store = 92;	// (half_arr_offset16 + RADIX) + 20 and round up to nearest multiple of 4
  #endif

	#include "sse2_macro.h"

	#ifdef COMPILER_TYPE_MSVC
		/*  */
	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix16_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix16_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

#endif	/* USE_SSE2 */

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
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
		int _pad_;

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *rn0;
		struct complex *rn1;
		vec_dbl *r00;
		vec_dbl *half_arr;

		int bjmodn0;
		int bjmodn1;
		int bjmodn2;
		int bjmodn3;
		int bjmodn4;
		int bjmodn5;
		int bjmodn6;
		int bjmodn7;
		int bjmodn8;
		int bjmodn9;
		int bjmodnA;
		int bjmodnB;
		int bjmodnC;
		int bjmodnD;
		int bjmodnE;
		int bjmodnF;
		/* carries: */
		double cy_r0;
		double cy_r1;
		double cy_r2;
		double cy_r3;
		double cy_r4;
		double cy_r5;
		double cy_r6;
		double cy_r7;
		double cy_r8;
		double cy_r9;
		double cy_rA;
		double cy_rB;
		double cy_rC;
		double cy_rD;
		double cy_rE;
		double cy_rF;

		double cy_i0;
		double cy_i1;
		double cy_i2;
		double cy_i3;
		double cy_i4;
		double cy_i5;
		double cy_i6;
		double cy_i7;
		double cy_i8;
		double cy_i9;
		double cy_iA;
		double cy_iB;
		double cy_iC;
		double cy_iD;
		double cy_iE;
		double cy_iF;
	};

#endif

/***************/

/* If using the FFT routines for a standalone build of the GCD code,
don't need the special-number carry routines:
*/
#ifdef GCD_STANDALONE

int radix16_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
	ASSERT(HERE, 0,"radix16_ditN_cy_dif1 should not be called if GCD_STANDALONE is set!");
	return 0;
}

int radix16_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter                 , uint64 p)
{
	ASSERT(HERE, 0,"radix16_ditN_cy_dif1_nochk should not be called if GCD_STANDALONE is set!");
	return 0;
}

#else

/***************/

int radix16_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-16 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-16 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix16_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const uint32 RADIX = 16;
	static int NDIVR;
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
	int i,j,j1,j2,jstart,jhi,full_pass,k,khi,l,outer,nbytes;
	int col,co2,co3;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */
	static uint32 bjmodnini;
	static uint64 psave=0;
	static uint32 bw,sw,nm1,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15;
	const  double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;
	static double radix_inv,n2inv,scale;	/* Need scale to be static since in the MSVC-only pure-ASM vesrion of the carry step save the address in a static pointer below */
	double dtmp,maxerr = 0.0;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;

#ifdef USE_SSE2

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
//	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	static int idx_offset, idx_incr;
//	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
  #endif

	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	static int *bjmodn0,*bjmodn1,*bjmodn2,*bjmodn3,*bjmodn4,*bjmodn5,*bjmodn6,*bjmodn7,*bjmodn8,*bjmodn9,*bjmodnA,*bjmodnB,*bjmodnC,*bjmodnD,*bjmodnE,*bjmodnF;
	static int *si_ptr;
	static double *wt0_ptr, *wt1_ptr, *scale_ptr = &scale;
	static vec_dbl *cc0, *ss0, *isrt2, *max_err, *sse2_rnd, *half_arr
		,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0A,*r0B,*r0C,*r0D,*r0E,*r0F
		,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1A,*r1B,*r1C,*r1D,*r1E,*r1F
		,*cy_r0,*cy_i0,*cy_r4,*cy_i4,*cy_r8,*cy_i8,*cy_rC,*cy_iC;
  #ifndef USE_AVX
	static vec_dbl
		 *cy_r2,*cy_i2,*cy_r6,*cy_i6,*cy_rA,*cy_iA,*cy_rE,*cy_iE;
  #else
	static vec_dbl *base_negacyclic_root;
	int k1,k2;
  #endif

	vec_dbl *tmp,*tm2;	// Non-static utility ptrs

#endif

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy16_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int jt,jp,k1,k2,m,m2,ntmp;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #if PFETCH
	double *addr, *addp;
  #endif
	int bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF;
	double temp,frac
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi
		,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF
		,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 nt_save = 0xffffffff, CY_THREADS = 0,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn0 = 0x0,*_bjmodn1 = 0x0,*_bjmodn2 = 0x0,*_bjmodn3 = 0x0,*_bjmodn4 = 0x0,*_bjmodn5 = 0x0,*_bjmodn6 = 0x0,*_bjmodn7 = 0x0,*_bjmodn8 = 0x0,*_bjmodn9 = 0x0,*_bjmodnA = 0x0,*_bjmodnB = 0x0,*_bjmodnC = 0x0,*_bjmodnD = 0x0,*_bjmodnE = 0x0,*_bjmodnF = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy_r0 = 0x0,*_cy_r1 = 0x0,*_cy_r2 = 0x0,*_cy_r3 = 0x0,*_cy_r4 = 0x0,*_cy_r5 = 0x0,*_cy_r6 = 0x0,*_cy_r7 = 0x0,*_cy_r8 = 0x0,*_cy_r9 = 0x0,*_cy_rA = 0x0,*_cy_rB = 0x0,*_cy_rC = 0x0,*_cy_rD = 0x0,*_cy_rE = 0x0,*_cy_rF = 0x0,
	*_cy_i0 = 0x0,*_cy_i1 = 0x0,*_cy_i2 = 0x0,*_cy_i3 = 0x0,*_cy_i4 = 0x0,*_cy_i5 = 0x0,*_cy_i6 = 0x0,*_cy_i7 = 0x0,*_cy_i8 = 0x0,*_cy_i9 = 0x0,*_cy_iA = 0x0,*_cy_iB = 0x0,*_cy_iC = 0x0,*_cy_iD = 0x0,*_cy_iE = 0x0,*_cy_iF = 0x0;

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/16;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/16.\n",iter);
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

	if(p != psave)	/* Exponent or #thread change triggers re-init */
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)16));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		nm1   = n-1;

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
		ASSERT(HERE, CY_THREADS >= NTHREADS,"CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, NDIVR    %CY_THREADS == 0,"NDIVR    %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"n_div_nwt%CY_THREADS != 0");
		}

	  #ifdef USE_PTHREAD

		j = (uint32)sizeof(struct cy_thread_data_t);
		if(0 != (j & 0xf)) {
			printf("sizeof(cy_thread_data_t) = %x\n",j);
			ASSERT(HERE, 0, "struct cy_thread_data_t not 16-byte size multiple!");
		}
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#ifdef OS_TYPE_MACOSX

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

			main_work_units = 0;
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
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of 128 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix16_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix16_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	next 16 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
		r00 = sc_ptr + 0x00;		r01 = sc_ptr + 0x01;
		r02 = sc_ptr + 0x02;		r03 = sc_ptr + 0x03;
		r04 = sc_ptr + 0x04;		r05 = sc_ptr + 0x05;
		r06 = sc_ptr + 0x06;		r07 = sc_ptr + 0x07;
		r08 = sc_ptr + 0x08;		r09 = sc_ptr + 0x09;
		r0A = sc_ptr + 0x0a;		r0B = sc_ptr + 0x0b;
		r0C = sc_ptr + 0x0c;		r0D = sc_ptr + 0x0d;
		r0E = sc_ptr + 0x0e;		r0F = sc_ptr + 0x0f;
		r10 = sc_ptr + 0x10;		r11 = sc_ptr + 0x11;
		r12 = sc_ptr + 0x12;		r13 = sc_ptr + 0x13;
		r14 = sc_ptr + 0x14;		r15 = sc_ptr + 0x15;
		r16 = sc_ptr + 0x16;		r17 = sc_ptr + 0x17;
		r18 = sc_ptr + 0x18;		r19 = sc_ptr + 0x19;
		r1A = sc_ptr + 0x1a;		r1B = sc_ptr + 0x1b;
		r1C = sc_ptr + 0x1c;		r1D = sc_ptr + 0x1d;
		r1E = sc_ptr + 0x1e;		r1F = sc_ptr + 0x1f;
	  isrt2 = sc_ptr + 0x20;
		cc0 = sc_ptr + 0x21;
		ss0 = sc_ptr + 0x22;	// +34 vec_dbl slots
		tmp = ss0 + 0x2;	// +36; Only need incr = 1 but prefer an even offset
	#ifdef USE_AVX
		 cy_r0 = tmp + 0x00;
		 cy_r4 = tmp + 0x01;
		 cy_r8 = tmp + 0x02;
		 cy_rC = tmp + 0x03;
		 cy_i0 = tmp + 0x04;
		 cy_i4 = tmp + 0x05;
		 cy_i8 = tmp + 0x06;
		 cy_iC = tmp + 0x07;
		max_err = tmp + 0x08;
		sse2_rnd= tmp + 0x09;
		half_arr= tmp + 0x0a;	// 36 + 10 = 46 = half_arr_offset16; This table needs 20x16 bytes
	#else
		 cy_r0 = tmp + 0x00;		cy_r2 = tmp + 0x01;
		 cy_r4 = tmp + 0x02;		cy_r6 = tmp + 0x03;
		 cy_r8 = tmp + 0x04;		cy_rA = tmp + 0x05;
		 cy_rC = tmp + 0x06;		cy_rE = tmp + 0x07;
		 cy_i0 = tmp + 0x08;		cy_i2 = tmp + 0x09;
		 cy_i4 = tmp + 0x0a;		cy_i6 = tmp + 0x0b;
		 cy_i8 = tmp + 0x0c;		cy_iA = tmp + 0x0d;
		 cy_iC = tmp + 0x0e;		cy_iE = tmp + 0x0f;
		max_err = tmp + 0x10;
		sse2_rnd= tmp + 0x11;
		half_arr= tmp + 0x12;	// 36 + 18 = 54 = half_arr_offset16; This table needs 20x16 bytes
	#endif
		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);
		VEC_DBL_INIT(cc0  , c	);
		VEC_DBL_INIT(ss0  , s	);
		VEC_DBL_INIT(sse2_rnd, crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)ss0 - (int)isrt2 + sz_vd;	// #bytes in 1st of above block of consts
		tmp = isrt2;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		nbytes = sz_vd;	// #bytes in 2nd block of consts, which contains just sse2_rnd
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

	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		/* In Fermat-mod mode, use first 2 128-bit slots for base and 1/base: */
		VEC_DBL_INIT(tmp, base   [0]);	++tmp;
		VEC_DBL_INIT(tmp, baseinv[0]);	++tmp;
		/* [+2] slot is for [scale,scale] */

		// Propagate the above consts to the remaining threads:
		nbytes = 2*sz_vd;
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#if 0 
	/* Simple qfloat-based loop to crank out the roots needed by th AVX negacyclic-mults scheme below: */
		struct qfloat qc,qs,qx,qy,qt,qn,qmul;
		qt = qfdiv(QPIHALF, i64_to_q((int64)RADIX));	// Pi/2/RADIX
		qc = qfcos(qt);	qs = qfsin(qt);
		qx = QONE;		qy = QZRO;
		for(j = 1; j < RADIX; j++) {
			// Up-multiply the complex exponential:
			qn = qfmul(qx, qc); qt = qfmul(qy, qs); qmul = qfsub(qn, qt);	// Store qxnew in qmul for now.
			qn = qfmul(qx, qs); qt = qfmul(qy, qc); qy   = qfadd(qn, qt); qx = qmul;
			printf("j = %3u: cos[j*Pi/2] = 0x%16llX, sin[j*Pi/2] = 0x%16llX\n",j,qfdbl_as_uint64(qx),qfdbl_as_uint64(qy));
		}
		exit(0);
	#endif

	#ifdef USE_AVX

		base_negacyclic_root = half_arr + RADIX;

		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of RADIX DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
		(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

		See the radix28 version of this routine for additional details.
		*/
		tmp = base_negacyclic_root + RADIX*2;	// First 120 = 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
		tm2 = tmp + RADIX/2 - 1;
										tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = 0x3FEFD88DA3D12526ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(01*I*Pi/32) = sin(31*I*Pi/32) */
		tmp64 = 0x3FEF6297CFF75CB0ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(02*I*Pi/32) = sin(30*I*Pi/32) */
		tmp64 = 0x3FEE9F4156C62DDAull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(03*I*Pi/32) = sin(29*I*Pi/32) */	tmp += 2;
		tmp64 = 0x3FED906BCF328D46ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(04*I*Pi/32) = sin(28*I*Pi/32) */	tm2 -= 2;
		tmp64 = 0x3FEC38B2F180BDB1ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(05*I*Pi/32) = sin(27*I*Pi/32) */
		tmp64 = 0x3FEA9B66290EA1A3ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(06*I*Pi/32) = sin(26*I*Pi/32) */
		tmp64 = 0x3FE8BC806B151741ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(07*I*Pi/32) = sin(25*I*Pi/32) */	tmp += 2;
		tmp64 = 0x3FE6A09E667F3BCDull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(08*I*Pi/32) = sin(24*I*Pi/32) */	tm2 -= 2;
		tmp64 = 0x3FE44CF325091DD6ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(09*I*Pi/32) = sin(23*I*Pi/32) */
		tmp64 = 0x3FE1C73B39AE68C8ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(10*I*Pi/32) = sin(22*I*Pi/32) */
		tmp64 = 0x3FDE2B5D3806F63Bull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(11*I*Pi/32) = sin(21*I*Pi/32) */	tmp += 2;
		tmp64 = 0x3FD87DE2A6AEA963ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(12*I*Pi/32) = sin(20*I*Pi/32) */	tm2 -= 2;
		tmp64 = 0x3FD294062ED59F06ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(13*I*Pi/32) = sin(19*I*Pi/32) */
		tmp64 = 0x3FC8F8B83C69A60Bull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(14*I*Pi/32) = sin(18*I*Pi/32) */
		tmp64 = 0x3FB917A6BC29B42Cull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(15*I*Pi/32) = sin(17*I*Pi/32) */	tmp += 2;

		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX << (l2_sz_vd-1);	// RADIX*sz_vd/2; 7 AVX-register-sized complex data

		// Propagate the above consts to the remaining threads:
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#endif	// AVX?
	}
	else
	{
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

	}	// TRANSFORM_TYPE toggle

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

	sse_nm1 = sse_sw + RE_IM_STRIDE;
	tmp64 = (uint64)nm1;
	tmp64 = tmp64 + (tmp64 << 32);
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sse_nm1 +i) = tmp64;
	}

	nbytes = 4 << l2_sz_vd;

#ifdef USE_AVX
	n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
	n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
	sinwt         = (struct uint32x4 *)sse_nm1 + 3;
	sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;
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
		bjmodn0 = (int*)(sinwtm1 + RE_IM_STRIDE);
	#else
		bjmodn0 = (int*)(sse_nm1 + RE_IM_STRIDE);
	#endif
		bjmodn1 = bjmodn0 + 1;
		bjmodn2 = bjmodn1 + 1;
		bjmodn3 = bjmodn2 + 1;
		bjmodn4 = bjmodn3 + 1;
		bjmodn5 = bjmodn4 + 1;
		bjmodn6 = bjmodn5 + 1;
		bjmodn7 = bjmodn6 + 1;
		bjmodn8 = bjmodn7 + 1;
		bjmodn9 = bjmodn8 + 1;
		bjmodnA = bjmodn9 + 1;
		bjmodnB = bjmodnA + 1;
		bjmodnC = bjmodnB + 1;
		bjmodnD = bjmodnC + 1;
		bjmodnE = bjmodnD + 1;
		bjmodnF = bjmodnE + 1;

	#endif	// USE_SSE2

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].r00 + ((long)half_arr - (long)r00);
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r00      = (vec_dbl *)base;
			tdat[ithread].half_arr = (vec_dbl *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

		/*   constant index offsets for load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
		p1  = NDIVR;
		p2  = p1  + p1;
		p3  = p2  + p1;
		p4  = p3  + p1;
	#ifndef USE_SSE2
		p5  = p4  + p1;
		p6  = p5  + p1;
		p7  = p6  + p1;
		p8  = p7  + p1;
		p9  = p8  + p1;
		p10 = p9  + p1;
		p11 = p10 + p1;
		p12 = p11 + p1;
		p13 = p12 + p1;
		p14 = p13 + p1;
		p15 = p14 + p1;
	#endif	// USE_SSE2
		p1  = p1  + ( (p1  >> DAT_BITS) << PAD_BITS );
		p2  = p2  + ( (p2  >> DAT_BITS) << PAD_BITS );
		p3  = p3  + ( (p3  >> DAT_BITS) << PAD_BITS );
		p4  = p4  + ( (p4  >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p5  = p5  + ( (p5  >> DAT_BITS) << PAD_BITS );
		p6  = p6  + ( (p6  >> DAT_BITS) << PAD_BITS );
		p7  = p7  + ( (p7  >> DAT_BITS) << PAD_BITS );
		p8  = p8  + ( (p8  >> DAT_BITS) << PAD_BITS );
		p9  = p9  + ( (p9  >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
	#endif	// USE_SSE2

		if(_cy_r0)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn0); _bjmodn0 = 0x0;
			free((void *)_bjmodn1); _bjmodn1 = 0x0;
			free((void *)_bjmodn2); _bjmodn2 = 0x0;
			free((void *)_bjmodn3); _bjmodn3 = 0x0;
			free((void *)_bjmodn4); _bjmodn4 = 0x0;
			free((void *)_bjmodn5); _bjmodn5 = 0x0;
			free((void *)_bjmodn6); _bjmodn6 = 0x0;
			free((void *)_bjmodn7); _bjmodn7 = 0x0;
			free((void *)_bjmodn8); _bjmodn8 = 0x0;
			free((void *)_bjmodn9); _bjmodn9 = 0x0;
			free((void *)_bjmodnA); _bjmodnA = 0x0;
			free((void *)_bjmodnB); _bjmodnB = 0x0;
			free((void *)_bjmodnC); _bjmodnC = 0x0;
			free((void *)_bjmodnD); _bjmodnD = 0x0;
			free((void *)_bjmodnE); _bjmodnE = 0x0;
			free((void *)_bjmodnF); _bjmodnF = 0x0;

			free((void *)_cy_r0); _cy_r0 = 0x0;	free((void *)_cy_i0); _cy_i0 = 0x0;
			free((void *)_cy_r1); _cy_r1 = 0x0;	free((void *)_cy_i1); _cy_i1 = 0x0;
			free((void *)_cy_r2); _cy_r2 = 0x0;	free((void *)_cy_i2); _cy_i2 = 0x0;
			free((void *)_cy_r3); _cy_r3 = 0x0;	free((void *)_cy_i3); _cy_i3 = 0x0;
			free((void *)_cy_r4); _cy_r4 = 0x0;	free((void *)_cy_i4); _cy_i4 = 0x0;
			free((void *)_cy_r5); _cy_r5 = 0x0;	free((void *)_cy_i5); _cy_i5 = 0x0;
			free((void *)_cy_r6); _cy_r6 = 0x0;	free((void *)_cy_i6); _cy_i6 = 0x0;
			free((void *)_cy_r7); _cy_r7 = 0x0;	free((void *)_cy_i7); _cy_i7 = 0x0;
			free((void *)_cy_r8); _cy_r8 = 0x0;	free((void *)_cy_i8); _cy_i8 = 0x0;
			free((void *)_cy_r9); _cy_r9 = 0x0;	free((void *)_cy_i9); _cy_i9 = 0x0;
			free((void *)_cy_rA); _cy_rA = 0x0;	free((void *)_cy_iA); _cy_iA = 0x0;
			free((void *)_cy_rB); _cy_rB = 0x0;	free((void *)_cy_iB); _cy_iB = 0x0;
			free((void *)_cy_rC); _cy_rC = 0x0;	free((void *)_cy_iC); _cy_iC = 0x0;
			free((void *)_cy_rD); _cy_rD = 0x0;	free((void *)_cy_iD); _cy_iD = 0x0;
			free((void *)_cy_rE); _cy_rE = 0x0;	free((void *)_cy_iE); _cy_iE = 0x0;
			free((void *)_cy_rF); _cy_rF = 0x0;	free((void *)_cy_iF); _cy_iF = 0x0;

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
		_bjmodn0	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn0== 0x0);
		_bjmodn1	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn1== 0x0);
		_bjmodn2	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn2== 0x0);
		_bjmodn3	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn3== 0x0);
		_bjmodn4	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn4== 0x0);
		_bjmodn5	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn5== 0x0);
		_bjmodn6	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn6== 0x0);
		_bjmodn7	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn7== 0x0);
		_bjmodn8	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn8== 0x0);
		_bjmodn9	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn9== 0x0);
		_bjmodnA	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnA== 0x0);
		_bjmodnB	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnB== 0x0);
		_bjmodnC	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnC== 0x0);
		_bjmodnD	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnD== 0x0);
		_bjmodnE	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnE== 0x0);
		_bjmodnF	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnF== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy_r0	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r0== 0x0);
		_cy_r1	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r1== 0x0);
		_cy_r2	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r2== 0x0);
		_cy_r3	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r3== 0x0);
		_cy_r4	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r4== 0x0);
		_cy_r5	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r5== 0x0);
		_cy_r6	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r6== 0x0);
		_cy_r7	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r7== 0x0);
		_cy_r8	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r8== 0x0);
		_cy_r9	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r9== 0x0);
		_cy_rA	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rA== 0x0);
		_cy_rB	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rB== 0x0);
		_cy_rC	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rC== 0x0);
		_cy_rD	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rD== 0x0);
		_cy_rE	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rE== 0x0);
		_cy_rF	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rF== 0x0);

		_cy_i0	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i0== 0x0);
		_cy_i1	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i1== 0x0);
		_cy_i2	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i2== 0x0);
		_cy_i3	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i3== 0x0);
		_cy_i4	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i4== 0x0);
		_cy_i5	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i5== 0x0);
		_cy_i6	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i6== 0x0);
		_cy_i7	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i7== 0x0);
		_cy_i8	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i8== 0x0);
		_cy_i9	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i9== 0x0);
		_cy_iA	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iA== 0x0);
		_cy_iB	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iB== 0x0);
		_cy_iC	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iC== 0x0);
		_cy_iD	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iD== 0x0);
		_cy_iE	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iE== 0x0);
		_cy_iF	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iF== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays.");

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
			i.e. the one that n2/16-separated FFT outputs need:
			*/
			_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
			_bjmodnini[0] = 0;
			_bjmodnini[1] = 0;
			for(j=0; j < NDIVR/CY_THREADS; j++)
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
			for(j=0; j < NDIVR; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
			ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");
		}
	}	/* endif(first_entry) */

/*...The radix-16 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r0[ithread] = 0;	_cy_i0[ithread] = 0;
		_cy_r1[ithread] = 0;	_cy_i1[ithread] = 0;
		_cy_r2[ithread] = 0;	_cy_i2[ithread] = 0;
		_cy_r3[ithread] = 0;	_cy_i3[ithread] = 0;
		_cy_r4[ithread] = 0;	_cy_i4[ithread] = 0;
		_cy_r5[ithread] = 0;	_cy_i5[ithread] = 0;
		_cy_r6[ithread] = 0;	_cy_i6[ithread] = 0;
		_cy_r7[ithread] = 0;	_cy_i7[ithread] = 0;
		_cy_r8[ithread] = 0;	_cy_i8[ithread] = 0;
		_cy_r9[ithread] = 0;	_cy_i9[ithread] = 0;
		_cy_rA[ithread] = 0;	_cy_iA[ithread] = 0;
		_cy_rB[ithread] = 0;	_cy_iB[ithread] = 0;
		_cy_rC[ithread] = 0;	_cy_iC[ithread] = 0;
		_cy_rD[ithread] = 0;	_cy_iD[ithread] = 0;
		_cy_rE[ithread] = 0;	_cy_iE[ithread] = 0;
		_cy_rF[ithread] = 0;	_cy_iF[ithread] = 0;

		_maxerr[ithread] = 0.0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r0[0] = -2;
	}

	*fracmax = 0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

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

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_bjmodn0[ithread] = _bjmodnini[ithread];	j = _bjmodnini[CY_THREADS];
			MOD_ADD32(_bjmodn0[ithread], j, n, _bjmodn1[ithread]);
			MOD_ADD32(_bjmodn1[ithread], j, n, _bjmodn2[ithread]);
			MOD_ADD32(_bjmodn2[ithread], j, n, _bjmodn3[ithread]);
			MOD_ADD32(_bjmodn3[ithread], j, n, _bjmodn4[ithread]);
			MOD_ADD32(_bjmodn4[ithread], j, n, _bjmodn5[ithread]);
			MOD_ADD32(_bjmodn5[ithread], j, n, _bjmodn6[ithread]);
			MOD_ADD32(_bjmodn6[ithread], j, n, _bjmodn7[ithread]);
			MOD_ADD32(_bjmodn7[ithread], j, n, _bjmodn8[ithread]);
			MOD_ADD32(_bjmodn8[ithread], j, n, _bjmodn9[ithread]);
			MOD_ADD32(_bjmodn9[ithread], j, n, _bjmodnA[ithread]);
			MOD_ADD32(_bjmodnA[ithread], j, n, _bjmodnB[ithread]);
			MOD_ADD32(_bjmodnB[ithread], j, n, _bjmodnC[ithread]);
			MOD_ADD32(_bjmodnC[ithread], j, n, _bjmodnD[ithread]);
			MOD_ADD32(_bjmodnD[ithread], j, n, _bjmodnE[ithread]);
			MOD_ADD32(_bjmodnE[ithread], j, n, _bjmodnF[ithread]);

			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		khi = 1;

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 15;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}
	}

#if defined(USE_SSE2) && defined(USE_PTHREAD)

	/* Init thread 1-CY_THREADS's local stores and pointers: */
   #if 1	// Targeted copy of just the per-iteration-need-to-re-init data in the local block:

	tmp = max_err;	VEC_DBL_INIT(tmp, 0.0);
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, sz_vd);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

   #else	// This is the copy-whole-local-block approach:

	VEC_DBL_INIT(max_err, 0.0);
	/* Init thread 1-CY_THREADS's local stores and pointers: */
	tmp = r00;	VEC_DBL_INIT(tmp, 0.0);
	nbytes = cslots_in_local_store << l2_sz_vd;
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, nbytes);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

   #endif

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
		ASSERT(HERE, tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(HERE, tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].r00;
		ASSERT(HERE, ((tmp + 0x20)->d0 == ISRT2 && (tmp + 0x20)->d1 == ISRT2), "thread-local memcheck failed!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	#endif

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			tdat[ithread].bjmodn0 = _bjmodn0[ithread];
			tdat[ithread].bjmodn1 = _bjmodn1[ithread];
			tdat[ithread].bjmodn2 = _bjmodn2[ithread];
			tdat[ithread].bjmodn3 = _bjmodn3[ithread];
			tdat[ithread].bjmodn4 = _bjmodn4[ithread];
			tdat[ithread].bjmodn5 = _bjmodn5[ithread];
			tdat[ithread].bjmodn6 = _bjmodn6[ithread];
			tdat[ithread].bjmodn7 = _bjmodn7[ithread];
			tdat[ithread].bjmodn8 = _bjmodn8[ithread];
			tdat[ithread].bjmodn9 = _bjmodn9[ithread];
			tdat[ithread].bjmodnA = _bjmodnA[ithread];
			tdat[ithread].bjmodnB = _bjmodnB[ithread];
			tdat[ithread].bjmodnC = _bjmodnC[ithread];
			tdat[ithread].bjmodnD = _bjmodnD[ithread];
			tdat[ithread].bjmodnE = _bjmodnE[ithread];
			tdat[ithread].bjmodnF = _bjmodnF[ithread];
			/* init carries	*/
			tdat[ithread].cy_r0 = _cy_r0[ithread];
			tdat[ithread].cy_r1 = _cy_r1[ithread];
			tdat[ithread].cy_r2 = _cy_r2[ithread];
			tdat[ithread].cy_r3 = _cy_r3[ithread];
			tdat[ithread].cy_r4 = _cy_r4[ithread];
			tdat[ithread].cy_r5 = _cy_r5[ithread];
			tdat[ithread].cy_r6 = _cy_r6[ithread];
			tdat[ithread].cy_r7 = _cy_r7[ithread];
			tdat[ithread].cy_r8 = _cy_r8[ithread];
			tdat[ithread].cy_r9 = _cy_r9[ithread];
			tdat[ithread].cy_rA = _cy_rA[ithread];
			tdat[ithread].cy_rB = _cy_rB[ithread];
			tdat[ithread].cy_rC = _cy_rC[ithread];
			tdat[ithread].cy_rD = _cy_rD[ithread];
			tdat[ithread].cy_rE = _cy_rE[ithread];
			tdat[ithread].cy_rF = _cy_rF[ithread];
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_SSE2
			// This is slightly different for power-of-2 DFTs: Here, scale is in the +2 slot, base & baseinv remain fixed in 0,+1 slots:
			dtmp = tmp->d0 * (tmp+1)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = tmp->d1 * (tmp+1)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			// scale gets set immediately prior to calling carry macro, hence no use checking it here.
		#endif
			/* init carries	*/
			tdat[ithread].cy_r0 = _cy_r0[ithread];	tdat[ithread].cy_i0 = _cy_i0[ithread];
			tdat[ithread].cy_r1 = _cy_r1[ithread];	tdat[ithread].cy_i1 = _cy_i1[ithread];
			tdat[ithread].cy_r2 = _cy_r2[ithread];	tdat[ithread].cy_i2 = _cy_i2[ithread];
			tdat[ithread].cy_r3 = _cy_r3[ithread];	tdat[ithread].cy_i3 = _cy_i3[ithread];
			tdat[ithread].cy_r4 = _cy_r4[ithread];	tdat[ithread].cy_i4 = _cy_i4[ithread];
			tdat[ithread].cy_r5 = _cy_r5[ithread];	tdat[ithread].cy_i5 = _cy_i5[ithread];
			tdat[ithread].cy_r6 = _cy_r6[ithread];	tdat[ithread].cy_i6 = _cy_i6[ithread];
			tdat[ithread].cy_r7 = _cy_r7[ithread];	tdat[ithread].cy_i7 = _cy_i7[ithread];
			tdat[ithread].cy_r8 = _cy_r8[ithread];	tdat[ithread].cy_i8 = _cy_i8[ithread];
			tdat[ithread].cy_r9 = _cy_r9[ithread];	tdat[ithread].cy_i9 = _cy_i9[ithread];
			tdat[ithread].cy_rA = _cy_rA[ithread];	tdat[ithread].cy_iA = _cy_iA[ithread];
			tdat[ithread].cy_rB = _cy_rB[ithread];	tdat[ithread].cy_iB = _cy_iB[ithread];
			tdat[ithread].cy_rC = _cy_rC[ithread];	tdat[ithread].cy_iC = _cy_iC[ithread];
			tdat[ithread].cy_rD = _cy_rD[ithread];	tdat[ithread].cy_iD = _cy_iD[ithread];
			tdat[ithread].cy_rE = _cy_rE[ithread];	tdat[ithread].cy_iE = _cy_iE[ithread];
			tdat[ithread].cy_rF = _cy_rF[ithread];	tdat[ithread].cy_iF = _cy_iF[ithread];
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

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_SSE2
			*bjmodn0 = _bjmodn0[ithread];
			*bjmodn1 = _bjmodn1[ithread];
			*bjmodn2 = _bjmodn2[ithread];
			*bjmodn3 = _bjmodn3[ithread];
			*bjmodn4 = _bjmodn4[ithread];
			*bjmodn5 = _bjmodn5[ithread];
			*bjmodn6 = _bjmodn6[ithread];
			*bjmodn7 = _bjmodn7[ithread];
			*bjmodn8 = _bjmodn8[ithread];
			*bjmodn9 = _bjmodn9[ithread];
			*bjmodnA = _bjmodnA[ithread];
			*bjmodnB = _bjmodnB[ithread];
			*bjmodnC = _bjmodnC[ithread];
			*bjmodnD = _bjmodnD[ithread];
			*bjmodnE = _bjmodnE[ithread];
			*bjmodnF = _bjmodnF[ithread];
		#else
			bjmodn0 = _bjmodn0[ithread];
			bjmodn1 = _bjmodn1[ithread];
			bjmodn2 = _bjmodn2[ithread];
			bjmodn3 = _bjmodn3[ithread];
			bjmodn4 = _bjmodn4[ithread];
			bjmodn5 = _bjmodn5[ithread];
			bjmodn6 = _bjmodn6[ithread];
			bjmodn7 = _bjmodn7[ithread];
			bjmodn8 = _bjmodn8[ithread];
			bjmodn9 = _bjmodn9[ithread];
			bjmodnA = _bjmodnA[ithread];
			bjmodnB = _bjmodnB[ithread];
			bjmodnC = _bjmodnC[ithread];
			bjmodnD = _bjmodnD[ithread];
			bjmodnE = _bjmodnE[ithread];
			bjmodnF = _bjmodnF[ithread];
		#endif
			/* init carries	*/
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_r1[ithread];	cy_r0->d2 = _cy_r2[ithread];	cy_r0->d3 = _cy_r3[ithread];
			cy_r4->d0 = _cy_r4[ithread];	cy_r4->d1 = _cy_r5[ithread];	cy_r4->d2 = _cy_r6[ithread];	cy_r4->d3 = _cy_r7[ithread];
			cy_r8->d0 = _cy_r8[ithread];	cy_r8->d1 = _cy_r9[ithread];	cy_r8->d2 = _cy_rA[ithread];	cy_r8->d3 = _cy_rB[ithread];
			cy_rC->d0 = _cy_rC[ithread];	cy_rC->d1 = _cy_rD[ithread];	cy_rC->d2 = _cy_rE[ithread];	cy_rC->d3 = _cy_rF[ithread];
		#elif defined(USE_SSE2)
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_r1[ithread];
			cy_r2->d0 = _cy_r2[ithread];	cy_r2->d1 = _cy_r3[ithread];
			cy_r4->d0 = _cy_r4[ithread];	cy_r4->d1 = _cy_r5[ithread];
			cy_r6->d0 = _cy_r6[ithread];	cy_r6->d1 = _cy_r7[ithread];
			cy_r8->d0 = _cy_r8[ithread];	cy_r8->d1 = _cy_r9[ithread];
			cy_rA->d0 = _cy_rA[ithread];	cy_rA->d1 = _cy_rB[ithread];
			cy_rC->d0 = _cy_rC[ithread];	cy_rC->d1 = _cy_rD[ithread];
			cy_rE->d0 = _cy_rE[ithread];	cy_rE->d1 = _cy_rF[ithread];
		#else
			cy_r0 = _cy_r0[ithread];
			cy_r1 = _cy_r1[ithread];
			cy_r2 = _cy_r2[ithread];
			cy_r3 = _cy_r3[ithread];
			cy_r4 = _cy_r4[ithread];
			cy_r5 = _cy_r5[ithread];
			cy_r6 = _cy_r6[ithread];
			cy_r7 = _cy_r7[ithread];
			cy_r8 = _cy_r8[ithread];
			cy_r9 = _cy_r9[ithread];
			cy_rA = _cy_rA[ithread];
			cy_rB = _cy_rB[ithread];
			cy_rC = _cy_rC[ithread];
			cy_rD = _cy_rD[ithread];
			cy_rE = _cy_rE[ithread];
			cy_rF = _cy_rF[ithread];
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#ifdef USE_AVX
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_r1[ithread];	cy_r0->d2 = _cy_r2[ithread];	cy_r0->d3 = _cy_r3[ithread];
			cy_r4->d0 = _cy_r4[ithread];	cy_r4->d1 = _cy_r5[ithread];	cy_r4->d2 = _cy_r6[ithread];	cy_r4->d3 = _cy_r7[ithread];
			cy_r8->d0 = _cy_r8[ithread];	cy_r8->d1 = _cy_r9[ithread];	cy_r8->d2 = _cy_rA[ithread];	cy_r8->d3 = _cy_rB[ithread];
			cy_rC->d0 = _cy_rC[ithread];	cy_rC->d1 = _cy_rD[ithread];	cy_rC->d2 = _cy_rE[ithread];	cy_rC->d3 = _cy_rF[ithread];

			cy_i0->d0 = _cy_i0[ithread];	cy_i0->d1 = _cy_i1[ithread];	cy_i0->d2 = _cy_i2[ithread];	cy_i0->d3 = _cy_i3[ithread];
			cy_i4->d0 = _cy_i4[ithread];	cy_i4->d1 = _cy_i5[ithread];	cy_i4->d2 = _cy_i6[ithread];	cy_i4->d3 = _cy_i7[ithread];
			cy_i8->d0 = _cy_i8[ithread];	cy_i8->d1 = _cy_i9[ithread];	cy_i8->d2 = _cy_iA[ithread];	cy_i8->d3 = _cy_iB[ithread];
			cy_iC->d0 = _cy_iC[ithread];	cy_iC->d1 = _cy_iD[ithread];	cy_iC->d2 = _cy_iE[ithread];	cy_iC->d3 = _cy_iF[ithread];
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_i0[ithread];
			cy_r2->d0 = _cy_r1[ithread];	cy_r2->d1 = _cy_i1[ithread];
			cy_r4->d0 = _cy_r2[ithread];	cy_r4->d1 = _cy_i2[ithread];
			cy_r6->d0 = _cy_r3[ithread];	cy_r6->d1 = _cy_i3[ithread];
			cy_r8->d0 = _cy_r4[ithread];	cy_r8->d1 = _cy_i4[ithread];
			cy_rA->d0 = _cy_r5[ithread];	cy_rA->d1 = _cy_i5[ithread];
			cy_rC->d0 = _cy_r6[ithread];	cy_rC->d1 = _cy_i6[ithread];
			cy_rE->d0 = _cy_r7[ithread];	cy_rE->d1 = _cy_i7[ithread];
			cy_i0->d0 = _cy_r8[ithread];	cy_i0->d1 = _cy_i8[ithread];
			cy_i2->d0 = _cy_r9[ithread];	cy_i2->d1 = _cy_i9[ithread];
			cy_i4->d0 = _cy_rA[ithread];	cy_i4->d1 = _cy_iA[ithread];
			cy_i6->d0 = _cy_rB[ithread];	cy_i6->d1 = _cy_iB[ithread];
			cy_i8->d0 = _cy_rC[ithread];	cy_i8->d1 = _cy_iC[ithread];
			cy_iA->d0 = _cy_rD[ithread];	cy_iA->d1 = _cy_iD[ithread];
			cy_iC->d0 = _cy_rE[ithread];	cy_iC->d1 = _cy_iE[ithread];
			cy_iE->d0 = _cy_rF[ithread];	cy_iE->d1 = _cy_iF[ithread];
		#else
			cy_r0 = _cy_r0[ithread];
			cy_r1 = _cy_r1[ithread];
			cy_r2 = _cy_r2[ithread];
			cy_r3 = _cy_r3[ithread];
			cy_r4 = _cy_r4[ithread];
			cy_r5 = _cy_r5[ithread];
			cy_r6 = _cy_r6[ithread];
			cy_r7 = _cy_r7[ithread];
			cy_r8 = _cy_r8[ithread];
			cy_r9 = _cy_r9[ithread];
			cy_rA = _cy_rA[ithread];
			cy_rB = _cy_rB[ithread];
			cy_rC = _cy_rC[ithread];
			cy_rD = _cy_rD[ithread];
			cy_rE = _cy_rE[ithread];
			cy_rF = _cy_rF[ithread];

			cy_i0 = _cy_i0[ithread];
			cy_i1 = _cy_i1[ithread];
			cy_i2 = _cy_i2[ithread];
			cy_i3 = _cy_i3[ithread];
			cy_i4 = _cy_i4[ithread];
			cy_i5 = _cy_i5[ithread];
			cy_i6 = _cy_i6[ithread];
			cy_i7 = _cy_i7[ithread];
			cy_i8 = _cy_i8[ithread];
			cy_i9 = _cy_i9[ithread];
			cy_iA = _cy_iA[ithread];
			cy_iB = _cy_iB[ithread];
			cy_iC = _cy_iC[ithread];
			cy_iD = _cy_iD[ithread];
			cy_iE = _cy_iE[ithread];
			cy_iF = _cy_iF[ithread];
		#endif
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
			(But only ever need to explicitly do this in debug mode).
			*/
			for(j = jstart; j < jhi; j += stride)
			{
				j1 =  j;
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1 + RE_IM_STRIDE;

		/*...The radix-32 DIT pass is here:	*/

		#ifdef USE_SSE2

		  #if defined(COMPILER_TYPE_MSVC)

			/*...Block 1: */
			#if 1
				add0 = &a[j1];
				__asm	mov eax, add0
				__asm	mov ebx, p1
				__asm	mov ecx, p2
				__asm	mov edx, p3
				__asm	mov	edi, p4		/* edi will store copy of p4 throughout */
				__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
				__asm	shl	ecx, 3
				__asm	shl	edx, 3
				__asm	shl	edi, 3
				__asm	add ebx, eax
				__asm	add ecx, eax
				__asm	add edx, eax
				SSE2_RADIX4_DIT_0TWIDDLE_B(r00)
			#else
				add0 = &a[j1];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;
				/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r00)
			#endif

			/*...Block 2: */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
				SSE2_RADIX4_DIT_0TWIDDLE_B(r08)
			#else
				add0 = &a[j1+p4];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;
				/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r08)
			#endif

			/*...Block 3: */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
				SSE2_RADIX4_DIT_0TWIDDLE_B(r10)
			#else
				add0 = &a[j1+p8];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;
				/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r10)
			#endif

			/*...Block 4: */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
				SSE2_RADIX4_DIT_0TWIDDLE_B(r18)
			#else
				add0 = &a[j1+p12];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;
				/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r18)
			#endif

			/****************************************************************************************************
			!...and now do four more radix-4 transforms, including the internal [no external]twiddle factors:   !
			****************************************************************************************************/

			/*...Block 1: Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				__asm	mov	eax, r00

				__asm	movaps	xmm0,[eax      ]	/* t1  */
				__asm	movaps	xmm1,[eax+0x010]	/* t2  */
				__asm	movaps	xmm2,[eax+0x080]	/* t9  */
				__asm	movaps	xmm3,[eax+0x090]	/* t10 */

				__asm	subpd	xmm0,[eax+0x080]	/*~t9 =t1 -t9 */
				__asm	subpd	xmm1,[eax+0x090]	/*~t10=t2 -t10*/
				__asm	addpd	xmm2,[eax      ]	/*~t1 =t9 +t1 */
				__asm	addpd	xmm3,[eax+0x010]	/*~t2 =t10+t2 */

				__asm	movaps	xmm4,[eax+0x100]	/* t17 */
				__asm	movaps	xmm5,[eax+0x110]	/* t18 */
				__asm	movaps	xmm6,[eax+0x180]	/* t25 */
				__asm	movaps	xmm7,[eax+0x190]	/* t26 */

				__asm	subpd	xmm4,[eax+0x180]	/*~t25=t17-t25*/
				__asm	subpd	xmm5,[eax+0x190]	/*~t26=t18-t26*/
				__asm	addpd	xmm6,[eax+0x100]	/*~t17=t25+t17*/
				__asm	addpd	xmm7,[eax+0x110]	/*~t18=t26+t18*/

			/*
				t1       =t1+t17;				t2       =t2+t18;
				t17     *=     2;				t18     *=     2;
				t17      =t1-t17;				t18      =t2-t18;
			*/
				__asm	subpd	xmm2,xmm6		/* t1 -t17 */
				__asm	subpd	xmm3,xmm7		/* t2 -t18 */
				__asm	movaps	[eax+0x100],xmm2	/* a[jt+p8 ]/t17 */
				__asm	movaps	[eax+0x110],xmm3	/* a[jp+p8 ]/t18 */
				__asm	addpd	xmm6,xmm6		/*   2*t17 */
				__asm	addpd	xmm7,xmm7		/*   2*t18 */
				__asm	addpd	xmm6,xmm2		/*~t17 <- t1 +t17 */
				__asm	addpd	xmm7,xmm3		/*~t18 <- t2 +t18 */
				__asm	movaps	[eax      ],xmm6	/* a[jt+p0 ]/t1  */
				__asm	movaps	[eax+0x010],xmm7	/* a[jp+p0 ]/t2  */

			/*
				t9       =t9 +t26;				t10      =t10-t25;	// mpy by E^-4 = -I is inlined here...
				t26     *=     2;				t25     *=     2;
				t26      =t9 -t26;				t25      =t10+t25;
			*/
				__asm	subpd	xmm0,xmm5		/* t9 -t26 */
				__asm	subpd	xmm1,xmm4		/* t10-t25 */
				__asm	movaps	[eax+0x180],xmm0	/* a[jt+p12]/t26 */
				__asm	movaps	[eax+0x090],xmm1	/* a[jp+p4 ]/t10 */
				__asm	addpd	xmm5,xmm5		/*   2*t26 */
				__asm	addpd	xmm4,xmm4		/*   2*t25 */
				__asm	addpd	xmm5,xmm0		/* t9 +t26 */
				__asm	addpd	xmm4,xmm1		/* t10+t25 */
				__asm	movaps	[eax+0x080],xmm5	/* a[jt+p4 ]/t9  */
				__asm	movaps	[eax+0x190],xmm4	/* a[jp+p12]/t25 */

			/*...Block 2: Cost: 19 MOVapd, 24 ADD/SUBpd,  4 MULpd */
				__asm	mov	eax, r04
				__asm	mov	ebx, isrt2
				__asm	movaps	xmm2,[ebx]	/* isrt2 */

				__asm	movaps	xmm4,[eax+0x100]	/* t21 */
				__asm	movaps	xmm5,[eax+0x110]	/* t22 */
				__asm	movaps	xmm0,[eax+0x180]	/* t29 */
				__asm	movaps	xmm1,[eax+0x190]	/* t30 */

				__asm	addpd	xmm4,[eax+0x110]	/*~t21=t21+t22*/
				__asm	subpd	xmm5,[eax+0x100]	/*~t22=t22-t21*/
				__asm	subpd	xmm0,[eax+0x190]	/* rt =t29-t30*/
				__asm	addpd	xmm1,[eax+0x180]	/* it =t30+t29*/
				__asm	mulpd	xmm4,xmm2
				__asm	mulpd	xmm5,xmm2
				__asm	mulpd	xmm0,xmm2
				__asm	mulpd	xmm1,xmm2
				__asm	movaps	xmm6,xmm4			/* t21 copy */
				__asm	movaps	xmm7,xmm5			/* t22 copy */

				__asm	subpd	xmm4,xmm0			/*~t21=t21-rt */
				__asm	subpd	xmm5,xmm1			/*~t22=t22-it */
				__asm	addpd	xmm6,xmm0			/*~t29=t21+rt */
				__asm	addpd	xmm7,xmm1			/*~t30=t22+it */

				__asm	movaps	xmm0,[eax      ]	/* t5  */
				__asm	movaps	xmm1,[eax+0x010]	/* t6  */
				__asm	movaps	xmm2,[eax+0x080]	/* t13 */
				__asm	movaps	xmm3,[eax+0x090]	/* t14 */

				__asm	subpd	xmm0,[eax+0x090]	/*~t13=t5 -t14*/
				__asm	subpd	xmm1,[eax+0x080]	/*~t6 =t6 -t13*/
				__asm	addpd	xmm3,[eax      ]	/*~t5 =t14+t5 */
				__asm	addpd	xmm2,[eax+0x010]	/*~t14=t13+t6 */

			/*
				t5       =t5 +t21;			t6       =t6 +t22;
				t21     *=     2;			t22     *=     2;
				t21      =t5 -t21;			t22      =t6 -t22;
			*/
				__asm	subpd	xmm3,xmm4		/*~t21 <- t5 -t21 */
				__asm	subpd	xmm1,xmm5		/*~t22 <- t6 -t22 */
				__asm	movaps	[eax+0x100],xmm3	/* a[jt+p8 ]/t21 */
				__asm	movaps	[eax+0x110],xmm1	/* a[jp+p8 ]/t22 */
				__asm	addpd	xmm4,xmm4		/*          2*t21 */
				__asm	addpd	xmm5,xmm5		/*          2*t22 */
				__asm	addpd	xmm4,xmm3		/*~t5  <- t5 +t21 */
				__asm	addpd	xmm5,xmm1		/*~t6  <- t6 +t22 */
				__asm	movaps	[eax      ],xmm4	/* a[jt+p0 ]/t5  */
				__asm	movaps	[eax+0x010],xmm5	/* a[jp+p0 ]/t6  */
			/*
				t13      =t13+t30;			t14      =t14-t29;	// mpy by E^-4 = -I is inlined here...
				t30     *=     2;			t29     *=     2;
				t30      =t13-t30;			t29      =t14+t29;
			*/
				__asm	subpd	xmm0,xmm7		/*~t30 <- t13-t30 */
				__asm	subpd	xmm2,xmm6		/*~t14 <- t14-t29 */
				__asm	movaps	[eax+0x180],xmm0	/* a[jt+p12]/t30 */
				__asm	movaps	[eax+0x090],xmm2	/* a[jp+p4 ]/t14 */
				__asm	addpd	xmm7,xmm7		/*          2*t30 */
				__asm	addpd	xmm6,xmm6		/*          2*t29 */
				__asm	addpd	xmm7,xmm0		/*~t13 <- t13+t30 */
				__asm	addpd	xmm6,xmm2		/*~t29 <- t14+t29 */
				__asm	movaps	[eax+0x080],xmm7	/* a[jt+p4 ]/t13 */
				__asm	movaps	[eax+0x190],xmm6	/* a[jp+p12]/t29 */

			/*...Block 3: Cost: 22 MOVapd, 28 ADD/SUBpd, 10 MULpd */
				__asm	mov	eax, r02
				__asm	mov	ebx, isrt2
				__asm	mov	ecx, cc0

				__asm	movaps	xmm4,[eax+0x100]	/* t19 */				__asm	movaps	xmm0,[eax+0x180]	/* t27 */
				__asm	movaps	xmm5,[eax+0x110]	/* t20 */				__asm	movaps	xmm1,[eax+0x190]	/* t28 */
				__asm	movaps	xmm6,[eax+0x100]	/* xmm2 <- cpy t19 */	__asm	movaps	xmm2,[eax+0x180]	/* xmm6 <- cpy t27 */
				__asm	movaps	xmm7,[eax+0x110]	/* xmm3 <- cpy t20 */	__asm	movaps	xmm3,[eax+0x190]	/* xmm7 <- cpy t28 */

				__asm	mulpd	xmm4,[ecx     ]	/* t19*c */					__asm	mulpd	xmm0,[ecx+0x10]	/* t27*s */
				__asm	mulpd	xmm5,[ecx     ]	/* t20*c */					__asm	mulpd	xmm1,[ecx+0x10]	/* t28*s */
				__asm	mulpd	xmm6,[ecx+0x10]	/* t19*s */					__asm	mulpd	xmm2,[ecx     ]	/* t27*c */
				__asm	mulpd	xmm7,[ecx+0x10]	/* t20*s */					__asm	mulpd	xmm3,[ecx     ]	/* t28*c */
				__asm	subpd	xmm5,xmm6	/* xmm1 <-~t20*/				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
				__asm	addpd	xmm4,xmm7	/* xmm0 <-~t19*/				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
				__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t20*/
				__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t19*/

				__asm	addpd	xmm4,xmm0	/* ~t19 <- t19+rt */
				__asm	addpd	xmm5,xmm1	/* ~t20 <- t20+it */
				__asm	subpd	xmm6,xmm0	/* ~t27 <- t19-rt */
				__asm	subpd	xmm7,xmm1	/* ~t28 <- t20-it */

				__asm	movaps	xmm2,[eax+0x080]	/* t11 */
				__asm	movaps	xmm3,[eax+0x090]	/* t12 */
				__asm	movaps	xmm0,[eax      ]	/* t3  */
				__asm	movaps	xmm1,[eax+0x010]	/* t4  */
				__asm	addpd	xmm2,[eax+0x090]	/*~t11=t11+t12*/
				__asm	subpd	xmm3,[eax+0x080]	/*~t12=t12-t11*/
				__asm	mulpd	xmm2,[ebx]	/* rt */
				__asm	mulpd	xmm3,[ebx]	/* it */

				__asm	subpd	xmm0,xmm2	/*~t11 <- t3 - rt */
				__asm	subpd	xmm1,xmm3	/*~t12 <- t4 - it */
				__asm	addpd	xmm2,xmm2	/*          2* rt */
				__asm	addpd	xmm3,xmm3	/*          2* it */
				__asm	addpd	xmm2,xmm0	/*~t3  <- t3 + rt */
				__asm	addpd	xmm3,xmm1	/*~t4  <- t4 + it */

			/*
				t3       =t3 +t19;			t4       =t4 +t20;
				t19     *=     2;			t20     *=     2;
				t19      =t3 -t19;			t20      =t4 -t20;
			*/
				__asm	subpd	xmm2,xmm4		/*~t19 <- t3 -t19 */
				__asm	subpd	xmm3,xmm5		/*~t20 <- t4 -t20 */
				__asm	movaps	[eax+0x100],xmm2	/* a[jt+p8 ]/t19 */
				__asm	movaps	[eax+0x110],xmm3	/* a[jp+p8 ]/t20 */
				__asm	addpd	xmm4,xmm4		/*          2*t19 */
				__asm	addpd	xmm5,xmm5		/*          2*t20 */
				__asm	addpd	xmm4,xmm2		/* rt  <- t3 +t19 */
				__asm	addpd	xmm5,xmm3		/* it  <- t4 +t20 */
				__asm	movaps	[eax      ],xmm4	/* a[jt    ]/t3  */
				__asm	movaps	[eax+0x010],xmm5	/* a[jp    ]/t4  */

			/*
				t11      =t11+t28;			t12      =t12-t27;	// mpy by E^-4 = -I is inlined here...
				t28     *=     2;			t27     *=     2;
				t28      =t11-t28;			t27      =t12+t27;
			*/
				__asm	subpd	xmm0,xmm7		/*~t28 <- t11-t28 */
				__asm	subpd	xmm1,xmm6		/*~t12 <- t12-t27 */
				__asm	movaps	[eax+0x180],xmm0	/* a[jt+p12]/t28 */
				__asm	movaps	[eax+0x090],xmm1	/* a[jp+p4 ]/t12 */
				__asm	addpd	xmm7,xmm7		/*          2*t28 */
				__asm	addpd	xmm6,xmm6		/*          2*t27 */
				__asm	addpd	xmm7,xmm0		/*~t11 <- t11+t28 */
				__asm	addpd	xmm6,xmm1		/*~t27 <- t12+t27 */
				__asm	movaps	[eax+0x080],xmm7	/* a[jt+p4 ]/t11 */
				__asm	movaps	[eax+0x190],xmm6	/* a[jp+p12]/t27 */

			/*...Block 4: Cost: 22 MOVapd, 28 ADD/SUBpd, 10 MULpd */
				__asm	mov	eax, r06
				__asm	mov	ebx, isrt2
				__asm	mov	ecx, cc0

				__asm	movaps	xmm4,[eax+0x100]	/* t23 */				__asm	movaps	xmm0,[eax+0x180]	/* t31 */
				__asm	movaps	xmm5,[eax+0x110]	/* t24 */				__asm	movaps	xmm1,[eax+0x190]	/* t32 */
				__asm	movaps	xmm6,[eax+0x100]	/* xmm2 <- cpy t23 */	__asm	movaps	xmm2,[eax+0x180]	/* xmm6 <- cpy t31 */
				__asm	movaps	xmm7,[eax+0x110]	/* xmm3 <- cpy t24 */	__asm	movaps	xmm3,[eax+0x190]	/* xmm7 <- cpy t32 */

				__asm	mulpd	xmm4,[ecx+0x10]	/* t23*s */					__asm	mulpd	xmm0,[ecx     ]	/* t31*c */
				__asm	mulpd	xmm5,[ecx+0x10]	/* t24*s */					__asm	mulpd	xmm1,[ecx     ]	/* t32*c */
				__asm	mulpd	xmm6,[ecx     ]	/* t23*c */					__asm	mulpd	xmm2,[ecx+0x10]	/* t31*s */
				__asm	mulpd	xmm7,[ecx     ]	/* t24*c */					__asm	mulpd	xmm3,[ecx+0x10]	/* t32*s */
				__asm	subpd	xmm5,xmm6	/* xmm1 <-~t24*/				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
				__asm	addpd	xmm4,xmm7	/* xmm0 <-~t23*/				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
				__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t24*/
				__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t23*/

				__asm	addpd	xmm4,xmm0	/* ~t31 <- t23+rt */
				__asm	addpd	xmm5,xmm1	/* ~t32 <- t24+it */
				__asm	subpd	xmm6,xmm0	/* ~t23 <- t23-rt */
				__asm	subpd	xmm7,xmm1	/* ~t24 <- t24-it */

				__asm	movaps	xmm2,[eax+0x080]	/* t15 */
				__asm	movaps	xmm3,[eax+0x090]	/* t16 */
				__asm	movaps	xmm0,[eax      ]	/* t7  */
				__asm	movaps	xmm1,[eax+0x010]	/* t8  */
				__asm	subpd	xmm2,[eax+0x090]	/*~t15=t15-t16*/
				__asm	addpd	xmm3,[eax+0x080]	/*~t16=t16+t15*/
				__asm	mulpd	xmm2,[ebx]	/* rt */
				__asm	mulpd	xmm3,[ebx]	/* it */

				__asm	subpd	xmm0,xmm2	/*~t7  <- t7 - rt */
				__asm	subpd	xmm1,xmm3	/*~t8  <- t8 - it */
				__asm	addpd	xmm2,xmm2	/*          2* rt */
				__asm	addpd	xmm3,xmm3	/*          2* it */
				__asm	addpd	xmm2,xmm0	/*~t15 <- t7 + rt */
				__asm	addpd	xmm3,xmm1	/*~t16 <- t8 + it */

			/*
				t7       =t7 +t23;			t8       =t8 +t24;
				t23     *=     2;			t24     *=     2;
				t23      =t7 -t23;			t24      =t8 -t24;
			*/
				__asm	subpd	xmm0,xmm6		/*~t23 <- t7 -t23 */
				__asm	subpd	xmm1,xmm7		/*~t24 <- t8 -t24 */
				__asm	movaps	[eax+0x100],xmm0	/* a[jt+p8 ]/t23 */
				__asm	movaps	[eax+0x110],xmm1	/* a[jp+p8 ]/t24 */
				__asm	addpd	xmm6,xmm6		/*          2*t23 */
				__asm	addpd	xmm7,xmm7		/*          2*t24 */
				__asm	addpd	xmm6,xmm0		/* rt  <- t7 +t23 */
				__asm	addpd	xmm7,xmm1		/* it  <- t8 +t24 */
				__asm	movaps	[eax      ],xmm6	/* a[jt+p0 ]/t7  */
				__asm	movaps	[eax+0x010],xmm7	/* a[jp+p0 ]/t8  */

			/*
				t15      =t15+t32;			t16      =t16-t31;	// mpy by E^-4 = -I is inlined here...
				t32     *=     2;			t31     *=     2;
				t32      =t15-t32;			t31      =t16+t31;
			*/
				__asm	subpd	xmm2,xmm5		/*~t32 <- t15-t32 */
				__asm	subpd	xmm3,xmm4		/*~t16 <- t16-t31 */
				__asm	movaps	[eax+0x180],xmm2	/* a[jt+p12]/t19 */
				__asm	movaps	[eax+0x090],xmm3	/* a[jp+p4 ]/t20 */
				__asm	addpd	xmm5,xmm5		/*          2*t32 */
				__asm	addpd	xmm4,xmm4		/*          2*t31 */
				__asm	addpd	xmm5,xmm2		/*~t15 <- t15+t32 */
				__asm	addpd	xmm4,xmm3		/*~t31 <- t16+t31 */
				__asm	movaps	[eax+0x080],xmm5	/* a[jt+p4 ]/t19 */
				__asm	movaps	[eax+0x190],xmm4	/* a[jp+p12]/t20 */

				/***************************************************/
				/* DIT Totals: 143 MOVapd, 180 ADD/SUBpd, 24 MULpd */
				/***************************************************/

		  #else	/* GCC-style inline ASM: */

				add0 = &a[j1];
				SSE2_RADIX16_DIT_NOTWIDDLE(add0,p1,p2,p3,p4,r00,r02,r04,r06,r08,r0a,r10,r18,isrt2,cc0);

		  #endif

		#else	/* !USE_SSE2 */

			#if USE_SCALAR_DFT_MACRO
				RADIX_16_DIT(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
							,a1p0r,a1p0i,a1p1r,a1p1i,a1p2r,a1p2i,a1p3r,a1p3i,a1p4r,a1p4i,a1p5r,a1p5i,a1p6r,a1p6i,a1p7r,a1p7i,a1p8r,a1p8i,a1p9r,a1p9i,a1pAr,a1pAi,a1pBr,a1pBi,a1pCr,a1pCi,a1pDr,a1pDi,a1pEr,a1pEi,a1pFr,a1pFi
							,c,s)
			#else
				/*...Block 1:	*/
				t1 =a[j1    ];	t2 =a[j2    ];
				rt =a[j1+p1 ];	it =a[j2+p1 ];
				t3 =t1 -rt;		t1 =t1 +rt;
				t4 =t2 -it;		t2 =t2 +it;

				t5 =a[j1+p2 ];	t6 =a[j2+p2 ];
				rt =a[j1+p3 ];	it =a[j2+p3 ];
				t7 =t5 -rt;  	t5 =t5 +rt;
				t8 =t6 -it;  	t6 =t6 +it;

				rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
				it =t6;	t6 =t2 -it;	t2 =t2 +it;

				rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
				t8 =t4 +rt;	t4 =t4 -rt;

				/*...Block 2:	*/
				t9 =a[j1+p4 ];	t10=a[j2+p4 ];
				rt =a[j1+p5 ];	it =a[j2+p5 ];
				t11=t9 -rt;		t9 =t9 +rt;
				t12=t10-it;		t10=t10+it;

				t13=a[j1+p6 ];	t14=a[j2+p6 ];
				rt =a[j1+p7 ];	it =a[j2+p7 ];
				t15=t13-rt;  	t13=t13+rt;
				t16=t14-it;		t14=t14+it;

				rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
				it =t14;	t14=t10-it;	t10=t10+it;

				rt =t15;	t15=t11-t16;	t11=t11+t16;
				t16=t12+rt;	t12=t12-rt;

				/*...Block 3:	*/
				t17=a[j1+p8 ];	t18=a[j2+p8 ];
				rt =a[j1+p9 ];	it =a[j2+p9 ];
				t19=t17-rt;		t17=t17+rt;
				t20=t18-it;		t18=t18+it;

				t21=a[j1+p10];	t22=a[j2+p10];
				rt =a[j1+p11];	it =a[j2+p11];
				t23=t21-rt;  	t21=t21+rt;
				t24=t22-it;		t22=t22+it;

				rt =t21;	t21=t17-rt;	t17=t17+rt;
				it =t22;	t22=t18-it;	t18=t18+it;

				rt =t23;	t23=t19-t24;	t19=t19+t24;
				t24=t20+rt;	t20=t20-rt;

				/*...Block 4:	*/
				t25=a[j1+p12];	t26=a[j2+p12];
				rt =a[j1+p13];	it =a[j2+p13];
				t27=t25-rt;		t25=t25+rt;
				t28=t26-it;		t26=t26+it;

				t29=a[j1+p14];	t30=a[j2+p14];
				rt =a[j1+p15];	it =a[j2+p15];
				t31=t29-rt;  	t29=t29+rt;
				t32=t30-it;		t30=t30+it;

				rt =t29;	t29=t25-rt;	t25=t25+rt;
				it =t30;	t30=t26-it;	t26=t26+it;

				rt =t31;	t31=t27-t32;	t27=t27+t32;
				t32=t28+rt;	t28=t28-rt;

				/*...and now do four more radix-4 transforms, including the internal twiddle factors:
				1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
				1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
				1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
				(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
				I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
											 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
				and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

				/*...Block 1: t1,9,17,25	*/
				rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
				it =t10;	t10=t2 -it;	t2 =t2 +it;

				rt =t25;	t25=t17-rt;	t17=t17+rt;
				it =t26;	t26=t18-it;	t18=t18+it;

				a1p0r =t1+t17;	a1p0i =t2+t18;
				a1p8r =t1-t17;	a1p8i =t2-t18;

				a1p4r =t9 +t26;	a1p4i =t10-t25;	/* mpy by E^-4 = -I is inlined here...	*/
				a1pCr=t9 -t26;	a1pCi=t10+t25;

				/*...Block 3: t5,13,21,29	*/
				rt =t13;	t13=t5 -t14;	t5 =t5 +t14;		/* twiddle mpy by E^4 =-I	*/
				t14=t6 +rt;	t6 =t6 -rt;

				rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;	/* twiddle mpy by E^-2	*/
				rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
				t29=t21+rt;		t21=t21-rt;			/* ...and get E^-6 by flipping signs here.	*/
				t30=t22+it;		t22=t22-it;

				a1p2r =t5+t21;	a1p2i =t6+t22;
				a1pAr=t5-t21;	a1pAi=t6-t22;

				a1p6r =t13+t30;	a1p6i =t14-t29;	/* mpy by E^-4 = -I is inlined here...	*/
				a1pEr=t13-t30;	a1pEi=t14+t29;

				/*...Block 2: t3,11,19,27	*/
				rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;		/* twiddle mpy by E^-2	*/
				t11=t3 -rt;		t3 =t3 +rt;
				t12=t4 -it;		t4 =t4 +it;

				rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;	/* twiddle mpy by E^-1	*/
				rt =t27*s + t28*c;	it =t28*s - t27*c;		/* twiddle mpy by E^-3	*/
				t27=t19-rt;		t19=t19+rt;
				t28=t20-it;		t20=t20+it;

				a1p1r =t3+t19;	a1p1i =t4+t20;
				a1p9r =t3-t19;	a1p9i =t4-t20;

				a1p5r =t11+t28;	a1p5i =t12-t27;	/* mpy by E^-4 = -I is inlined here...	*/
				a1pDr=t11-t28;	a1pDi=t12+t27;

				/*...Block 4: t7,15,23,31	*/
				rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
				t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
				t16=t8 +it;		t8 =t8 -it;

				rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;	/* twiddle mpy by E^-3	*/
				rt =t31*c + t32*s;	it =t32*c - t31*s;		/* twiddle mpy by E^-1 = -E^-9...	*/
				t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
				t32=t24+it;		t24=t24-it;

				a1p3r =t7+t23;	a1p3i =t8+t24;
				a1pBr=t7-t23;	a1pBi=t8-t24;

				a1p7r =t15+t32;	a1p7i =t16-t31;	/* mpy by E^-4 = -I is inlined here...	*/
				a1pFr=t15-t32;	a1pFi=t16+t31;
			#endif	// USE_SCALAR_DFT_MACRO ?

		#endif	/* USE_SSE2 */

			/*...and combine those to complete the radix-16 transform and do the carries. Since the outputs would
			normally be getting dispatched to 16 separate blocks of the A-array, we need 16 separate carries.	*/

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
			#ifdef USE_AVX

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

				l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
				n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				AVX_cmplx_carry_norm_pow2_errcheck0_X4(r00,add1,add2,add3,cy_r0,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(r08,add1,add2,add3,cy_r4,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(r10,add1,add2,add3,cy_r8,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(r18,add1,add2,add3,cy_rC,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

				i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/

			#elif defined(USE_SSE2)

				/* In SSE2 mode, the data are arranged in memory like so, where we view things in 16-byte chunks:

					R0:	a0.re,b0.re		I0:	a0.im,b0.im
					R1:	a1.re,b1.re		I1:	a1.im,b1.im
					R2:	a2.re,b2.re		I2:	a2.im,b2.im
					R3:	a3.re,b3.re		I3:	a3.im,b3.im
					R4:	a4.re,b4.re		I4:	a4.im,b4.im
					R5:	a5.re,b5.re		I5:	a5.im,b5.im
					R6:	a6.re,b6.re		I6:	a6.im,b6.im
					R7:	a7.re,b7.re		I7:	a7.im,b7.im
					R8:	a8.re,b8.re		I8:	a8.im,b8.im
					R9:	a9.re,b9.re		I9:	a9.im,b9.im
					Ra:	aA.re,bA.re		Ia:	aA.im,bA.im
					Rb:	aB.re,bB.re		Ib:	aB.im,bB.im
					Rc:	aC.re,bC.re		Ic:	aC.im,bC.im
					Rd:	aD.re,bD.re		Id:	aD.im,bD.im
					Re:	aE.re,bE.re		Ie:	aE.im,bE.im
					Rf:	aF.re,bF.re		If:	aF.im,bF.im

				Where the R's and I's map to the local temps as follows: R0:f ==> r00:1E:2, I0:f ==> r01:1F:2 , and the
				a's and b's of each pair represent doubles which in non-SSE2 mode would be getting processed in the same relative
				position on subsequent passes through the main-array for-loop, i.e. among which the carry propagation proceeds as

					a0.re -> a0.im -> b0.re -> b0.im .

				Because of the undesirable intra-xmm-register data dependency this leads to, we instead require data arranged as

					R0/r00:	a0.re,a1.re		I0/r01:	a0.im,a1.im
					R1/r02:	b0.re,b1.re		I1/r03:	b0.im,b1.im

				We need to interleave these pairwise so as to swap the high word of each even-indexed R-and-I-pair
				with the low word of the subsequent odd-indexed pair, e.g. for R0/r00 and R1/r02:

						low		high	low		high
					R0	[a0.re,b0.re]	[a1.re,b1.re]	R1
						   |      \       /      |
						   |        \   /        |
						   |          x          |
						   |        /   \        |
						   V      /       \      V
					R0~	[a0.re,a1.re]	[b0.re,b1.re]	R1~, and analogously for I0/r01 and I1/r03.

				This is the same butterfly swap pattern as is used in the wrapper_square routines. The other nice things about this:

					1) Even though e.g. a0 and a1 appear adjacent, they are actually n/16 memory locations apart, i.e. there
					   is no carry propagation between them;

					2) Processing a[j] and a[j+1] together means we access the following elements of the wt1[] array paiwise in the carry step:

						xmm.lo:			xmm.hi:
						wt1[col+j]		wt1[col+(j+1)]
						wt1[co2-j]		wt1[co2-(j+1)]
						wt1[co3-j]		wt1[co3-(j+1)]

					Thus these wt-array elements are also adjacent in memory and can be loaded pairwise into an XMM register
					[With an unaligned movupd load and a shufpd-based lo/hi-word swap needed on the latter two.]
				*/
				/* These indices remain constant throughout the carry block below - only change when loop index j does: */

				l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			  #if defined(COMPILER_TYPE_MSVC)

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC);
			   #else
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC);
			   #endif

			  #else	/* GCC-style inline ASM: */

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #endif

			  #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			  #if defined(COMPILER_TYPE_MSVC)

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r00,add1,add2,cy_r0,cy_r2,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r08,add1,add2,cy_r4,cy_r6,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r10,add1,add2,cy_r8,cy_rA,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r18,add1,add2,cy_rC,cy_rE,bjmodnC);
			   #else
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r00,add1,add2,cy_r0,cy_r2,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r08,add1,add2,cy_r4,cy_r6,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r10,add1,add2,cy_r8,cy_rA,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r18,add1,add2,cy_rC,cy_rE,bjmodnC);
			   #endif

			  #else	/* GCC-style inline ASM: */

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r00,add1,add2,cy_r0,cy_r2,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r08,add1,add2,cy_r4,cy_r6,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r10,add1,add2,cy_r8,cy_rA,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r18,add1,add2,cy_rC,cy_rE,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r00,add1,add2,cy_r0,cy_r2,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r08,add1,add2,cy_r4,cy_r6,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r10,add1,add2,cy_r8,cy_rA,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r18,add1,add2,cy_rC,cy_rE,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #endif

			  #endif

				i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/

			#else	// Scalar-double mode:

				l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

			/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_pow2_errcheck0(a1p0r,a1p0i,cy_r0,bjmodn0    );
				cmplx_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,bjmodn1,0x1);
				cmplx_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,bjmodn2,0x2);
				cmplx_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,bjmodn3,0x3);
				cmplx_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,bjmodn4,0x4);
				cmplx_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,bjmodn5,0x5);
				cmplx_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,bjmodn6,0x6);
				cmplx_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,bjmodn7,0x7);
				cmplx_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,bjmodn8,0x8);
				cmplx_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,bjmodn9,0x9);
				cmplx_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,bjmodnA,0xA);
				cmplx_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,bjmodnB,0xB);
				cmplx_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,bjmodnC,0xC);
				cmplx_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,bjmodnD,0xD);
				cmplx_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,bjmodnE,0xE);
				cmplx_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,bjmodnF,0xF);

				i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			#endif	// USE_AVX?

			}
			else	/* MODULUS_TYPE_FERMAT */
			{

			#ifdef USE_AVX

				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.
				tmp = half_arr+2;
				VEC_DBL_INIT(tmp, scale);
				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				tmp = base_negacyclic_root;	tm2 = tmp+1;

				// Hi-accuracy version needs 4 copies of each base root:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);

				// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
				// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:

				// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
				// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
				// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
				// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
				tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r00,tmp,0x400,cy_r0,cy_i0,half_arr,sign_mask);
				tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r08,tmp,0x340,cy_r4,cy_i4,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r10,tmp,0x280,cy_r8,cy_i8,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r18,tmp,0x1c0,cy_rC,cy_iC,half_arr,sign_mask);

			#elif defined(USE_SSE2)

				/* In SSE2 mode, carry propagation proceeds as

					a0.re -> b0.re;		a0.im -> b0.im, where these imaginary parts really represent elements
					                                    a0.im = a[n/2] and b0.im = a[n/2+1] of the right-angle transform.

				This data layout is ideal for the negacyclic unweighting/reweighting step bracketing the carry step, but in the latter,
				because of the undesirable intra-xmm-register data dependency this leads to, we instead require data arranged as

					R0/r00:	a0.re,a0.im		I0/r01:	b0.re,b0.im, i.e. the non-SSE2 data layout works best in the carry step!

				We need to interleave these pairwise so as to swap the high word of each R-element
				with the low word of the corresponding I-element, e.g. for R0/r00 and I0/r01:

						low		high	low		high
					R0	[a0.re,b0.re]	[a0.im,b0.im]	I0
						   |      \       /      |
						   |        \   /        |
						   |          x          |
						   |        /   \        |
						   V      /       \      V
					R0~	[a0.re,a0.im]	[b0.re,b0.im]	I0~.

				Note that even though e.g. a0 and a1 appear adjacent in terms of their a-subscripts, they are actually
				n/16 memory locations apart, i.e. there is no carry propagation between them.
				*/

				tmp = half_arr+2;
				VEC_DBL_INIT(tmp, scale);
				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

			  #if defined(COMPILER_TYPE_MSVC)
				/* The cy_[r|i]_idx[A|B] names here are not meaningful, each simple stores one [re,im] carry pair,
				e.g. cy_r0 stores the carries our of [a0.re,a0.im], cy_r2 stores the carries our of [a1.re,a1.im], etc.
				Here is the actual mapping between these SSE2-mode 2-vector carry pairs and the scalar carries:
													  2-vector                               Scalar
													 ----------                            ----------- */
				SSE2_fermat_carry_norm_pow2_errcheck(r00,cy_r0,idx_offset,idx_incr);	/* cy_r0,cy_i0 */
				SSE2_fermat_carry_norm_pow2_errcheck(r02,cy_r2,idx_offset,idx_incr);	/* cy_r1,cy_i1 */
				SSE2_fermat_carry_norm_pow2_errcheck(r04,cy_r4,idx_offset,idx_incr);	/* cy_r2,cy_i2 */
				SSE2_fermat_carry_norm_pow2_errcheck(r06,cy_r6,idx_offset,idx_incr);	/* cy_r3,cy_i3 */
				SSE2_fermat_carry_norm_pow2_errcheck(r08,cy_r8,idx_offset,idx_incr);	/* cy_r4,cy_i4 */
				SSE2_fermat_carry_norm_pow2_errcheck(r0A,cy_rA,idx_offset,idx_incr);	/* cy_r5,cy_i5 */
				SSE2_fermat_carry_norm_pow2_errcheck(r0C,cy_rC,idx_offset,idx_incr);	/* cy_r6,cy_i6 */
				SSE2_fermat_carry_norm_pow2_errcheck(r0E,cy_rE,idx_offset,idx_incr);	/* cy_r7,cy_i7 */
				SSE2_fermat_carry_norm_pow2_errcheck(r10,cy_i0,idx_offset,idx_incr);	/* cy_r8,cy_i8 */
				SSE2_fermat_carry_norm_pow2_errcheck(r12,cy_i2,idx_offset,idx_incr);	/* cy_r9,cy_i9 */
				SSE2_fermat_carry_norm_pow2_errcheck(r14,cy_i4,idx_offset,idx_incr);	/* cy_rA,cy_iA */
				SSE2_fermat_carry_norm_pow2_errcheck(r16,cy_i6,idx_offset,idx_incr);	/* cy_rB,cy_iB */
				SSE2_fermat_carry_norm_pow2_errcheck(r18,cy_i8,idx_offset,idx_incr);	/* cy_rC,cy_iC */
				SSE2_fermat_carry_norm_pow2_errcheck(r1A,cy_iA,idx_offset,idx_incr);	/* cy_rD,cy_iD */
				SSE2_fermat_carry_norm_pow2_errcheck(r1C,cy_iC,idx_offset,idx_incr);	/* cy_rE,cy_iE */
				SSE2_fermat_carry_norm_pow2_errcheck(r1E,cy_iE,idx_offset,idx_incr);	/* cy_rF,cy_iF */

			  #elif (OS_BITS == 32)

				SSE2_fermat_carry_norm_pow2_errcheck(r00,cy_r0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r02,cy_r2,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r04,cy_r4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r06,cy_r6,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r08,cy_r8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r0A,cy_rA,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r0C,cy_rC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r0E,cy_rE,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r10,cy_i0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r12,cy_i2,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r14,cy_i4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r16,cy_i6,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r18,cy_i8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r1A,cy_iA,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r1C,cy_iC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck(r1E,cy_iE,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);

			  #else	// 64-bit SSE2

				SSE2_fermat_carry_norm_pow2_errcheck_X2(r00,cy_r0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r04,cy_r4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r08,cy_r8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r0C,cy_rC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r10,cy_i0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r14,cy_i4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r18,cy_i8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r1C,cy_iC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);

			  #endif

			#else	// Scalar-double mode:

				ntmp = 0;
				fermat_carry_norm_pow2_errcheck(a1p0r,a1p0i,cy_r0,cy_i0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,cy_i1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,cy_i2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,cy_i3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,cy_i4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,cy_i5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,cy_i6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,cy_i7,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,cy_i8,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,cy_i9,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,cy_iA,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,cy_iB,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,cy_iC,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,cy_iD,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,cy_iE,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,cy_iF,ntmp,NRTM1,NRT_BITS);

			#endif	/* #ifdef USE_SSE2 */

			}	/* if(MODULUS_TYPE == ...) */

		/*...The radix-16 DIF pass is here:	*/

		/* Four DIF radix-4 subconvolution, sans twiddles.	Cost each: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */

		#ifdef USE_SSE2

		  #if defined(COMPILER_TYPE_MSVC)

				SSE2_RADIX4_DIF_IN_PLACE(r00, r10, r08, r18)
				SSE2_RADIX4_DIF_IN_PLACE(r04, r14, r0C, r1C)
				SSE2_RADIX4_DIF_IN_PLACE(r02, r12, r0A, r1A)
				SSE2_RADIX4_DIF_IN_PLACE(r06, r16, r0E, r1E)

			/****************************************************************************************
			!...and now do four more radix-4 transforms, including the internal twiddle factors.	!
			!																						!
			!	This is identical to latter half of radix16 DIF, except for the r-vector indexing,	!
			!	which permutes as follows:															!
			!																						!
			!			t1	t3	t5	t7	t9	t11	t13	t15	t17	t19	t21	t23	t25	t27	t29	t31				!
			!		==>	r00	r08	r10	r18	r04	r0C	r14	r1C	r02	r0A	r12	r1A	r06	r0E	r16	r1E				!
			!																						!
			****************************************************************************************/

			/*...Block 1: t1,9,17,25 ==> r00,04,02,06	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
			#if 1	// if(1) - test out pure-asm version
			//	add0 = &a[j1];
				__asm	mov eax, add0
				__asm	mov ebx, p1
				__asm	mov ecx, p2
				__asm	mov edx, p3
				__asm	mov	edi, p4		/* edi will store copy of p4 throughout */
				__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
				__asm	shl	ecx, 3
				__asm	shl	edx, 3
				__asm	shl	edi, 3
				__asm	add ebx, eax
				__asm	add ecx, eax
				__asm	add edx, eax
			#else
				add0 = &a[j1];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;

				__asm	mov	eax, add0	/* restore main-array indices */
				__asm	mov	ebx, add1
				__asm	mov	ecx, add2
				__asm	mov	edx, add3
			#endif
				__asm	mov	esi, r00

				__asm	movaps	xmm0,[esi      ]	/* t1  */
				__asm	movaps	xmm1,[esi+0x010]	/* t2  */
				__asm	movaps	xmm2,[esi+0x040]	/* t9  */
				__asm	movaps	xmm3,[esi+0x050]	/* t10 */

				__asm	subpd	xmm0,[esi+0x040]	/* t9 =t1 -rt */
				__asm	subpd	xmm1,[esi+0x050]	/* t10=t2 -it */
				__asm	addpd	xmm2,[esi      ]	/* t1 =t1 +rt */
				__asm	addpd	xmm3,[esi+0x010]	/* t2 =t2 +it */

				__asm	movaps	xmm4,[esi+0x020]	/* t17 */
				__asm	movaps	xmm5,[esi+0x030]	/* t18 */
				__asm	movaps	xmm6,[esi+0x060]	/* t25 */
				__asm	movaps	xmm7,[esi+0x070]	/* t26 */

				__asm	subpd	xmm4,[esi+0x060]	/* t25=t17-rt */
				__asm	subpd	xmm5,[esi+0x070]	/* t26=t18-it */
				__asm	addpd	xmm6,[esi+0x020]	/* t17=t17+rt */
				__asm	addpd	xmm7,[esi+0x030]	/* t18=t18+it */

				__asm	subpd	xmm2,xmm6		/* t1  <- t1 -t17 */
				__asm	subpd	xmm3,xmm7		/* t2  <- t2 -t18 */
				__asm	addpd	xmm6,xmm6		/*          2*t17 */
				__asm	addpd	xmm7,xmm7		/*          2*t18 */
				__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */
				__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */
				__asm	addpd	xmm6,xmm2		/* t17 <- t1 +t17 */
				__asm	addpd	xmm7,xmm3		/* t18 <- t2 +t18 */
				__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */
				__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */

				__asm	subpd	xmm0,xmm5		/* t9  <- t9 -t26 */
				__asm	subpd	xmm1,xmm4		/* t10 <- t10-t25 */
				__asm	addpd	xmm5,xmm5		/*          2*t26 */
				__asm	addpd	xmm4,xmm4		/*          2*t25 */
				__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
				__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
				__asm	addpd	xmm5,xmm0		/* t26 <- t9 +t26 */
				__asm	addpd	xmm4,xmm1		/* t25 <- t10+t25 */
				__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
				__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

			/*...Block 3: t5,13,21,29 ==> r10,14,12,16	Cost: 16 MOVapd, 26 ADD/SUBpd,  4 MULpd */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
			#else
				add0 = &a[j1+p4];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;

				__asm	mov	eax, add0	/* restore main-array indices */
				__asm	mov	ebx, add1
				__asm	mov	ecx, add2
				__asm	mov	edx, add3
			#endif
				__asm	mov	esi, r10

				__asm	movaps	xmm0,[esi      ]	/* t5  */
				__asm	movaps	xmm1,[esi+0x010]	/* t6  */
				__asm	movaps	xmm2,[esi+0x040]	/* t13 */
				__asm	movaps	xmm3,[esi+0x050]	/* t14 */

				__asm	subpd	xmm0,[esi+0x050]	/* t5 =t5 -t14*/
				__asm	subpd	xmm1,[esi+0x040]	/* t14=t6 -t13*/
				__asm	addpd	xmm2,[esi+0x010]	/* t6 =t13+t6 */
				__asm	addpd	xmm3,[esi      ]	/* t13=t14+t5 */

				__asm	movaps	xmm4,[esi+0x020]	/* t21 */
				__asm	movaps	xmm5,[esi+0x030]	/* t22 */
				__asm	movaps	xmm6,[esi+0x060]	/* t29 */
				__asm	movaps	xmm7,[esi+0x070]	/* t30 */

				__asm	subpd	xmm4,[esi+0x030]	/* t21-t22 */
				__asm	addpd	xmm5,[esi+0x020]	/* t22+t21 */
				__asm	addpd	xmm6,[esi+0x070]	/* t29+t30 */
				__asm	subpd	xmm7,[esi+0x060]	/* t30-t29 */

				__asm	mov	esi, isrt2
				__asm	mulpd	xmm4,[esi]	/* t21 = (t21-t22)*ISRT2 */
				__asm	mulpd	xmm5,[esi]	/* t22 = (t22+t21)*ISRT2 */
				__asm	mulpd	xmm6,[esi]	/*  rt = (t29+t30)*ISRT2 */
				__asm	mulpd	xmm7,[esi]	/*  it = (t30-t29)*ISRT2 */

				__asm	subpd	xmm4,xmm6		/* t21=t21-rt */
				__asm	subpd	xmm5,xmm7		/* t22=t22-it */
				__asm	addpd	xmm6,xmm6		/*      2* rt */
				__asm	addpd	xmm7,xmm7		/*      2* it */
				__asm	addpd	xmm6,xmm4		/* t29=t21+rt */
				__asm	addpd	xmm7,xmm5		/* t30=t22+it */

				__asm	subpd	xmm0,xmm4		/* t5 -t21 */
				__asm	subpd	xmm2,xmm5		/* t6 -t22 */
				__asm	addpd	xmm4,xmm4		/*   2*t21 */
				__asm	addpd	xmm5,xmm5		/*   2*t22 */

				__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */
				__asm	movaps	[ebx+0x10],xmm2	/* a[jp+p1 ] */
				__asm	addpd	xmm4,xmm0		/* t5 +t21 */
				__asm	addpd	xmm5,xmm2		/* t6 +t22 */
				__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */
				__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */

				__asm	subpd	xmm3,xmm7		/* t13-t30 */
				__asm	subpd	xmm1,xmm6		/* t14-t29 */
				__asm	addpd	xmm7,xmm7		/*   2*t30 */
				__asm	addpd	xmm6,xmm6		/*   2*t29 */
				__asm	movaps	[ecx     ],xmm3	/* a[jt+p2 ] */
				__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
				__asm	addpd	xmm7,xmm3		/* t13+t30 */
				__asm	addpd	xmm6,xmm1		/* t14+t29 */
				__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
				__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

			/*...Block 2: t3,11,19,27 ==> r08,0C,0A,0E	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
			#else
				add0 = &a[j1+p8];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;

				__asm	mov	eax, add0	/* restore main-array indices */
				__asm	mov	ebx, add1
				__asm	mov	ecx, add2
				__asm	mov	edx, add3
			#endif
				__asm	mov	esi, r08
				__asm	mov	edi, cc0
				__asm	movaps	xmm4,[esi+0x020]	/* t19 */		__asm	movaps	xmm6,[esi+0x060]	/* t27 */
				__asm	movaps	xmm5,[esi+0x030]	/* t20 */		__asm	movaps	xmm7,[esi+0x070]	/* t28 */
				__asm	movaps	xmm0,xmm4			/* copy t19 */	__asm	movaps	xmm2,xmm6			/* copy t27 */
				__asm	movaps	xmm1,xmm5			/* copy t20 */	__asm	movaps	xmm3,xmm7			/* copy t28 */

				__asm	mulpd	xmm4,[edi     ]	/* t19*c */			__asm	mulpd	xmm6,[edi+0x10]	/* t27*s */
				__asm	mulpd	xmm1,[edi+0x10]	/* t20*s */			__asm	mulpd	xmm3,[edi     ]	/* t28*c */
				__asm	mulpd	xmm5,[edi     ]	/* t20*c */			__asm	mulpd	xmm7,[edi+0x10]	/* t28*s */
				__asm	mulpd	xmm0,[edi+0x10]	/* t19*s */			__asm	mulpd	xmm2,[edi     ]	/* t27*c */
				__asm	subpd	xmm4,xmm1	/* ~t19 */				__asm	subpd	xmm6,xmm3	/* rt */
				__asm	addpd	xmm5,xmm0	/* ~t20 */				__asm	addpd	xmm7,xmm2	/* it */

				__asm	subpd	xmm4,xmm6		/*~t27=t19-rt */
				__asm	subpd	xmm5,xmm7		/*~t28=t20-it */
				__asm	addpd	xmm6,xmm6		/*      2* rt */
				__asm	addpd	xmm7,xmm7		/*      2* it */
				__asm	addpd	xmm6,xmm4		/*~t19=t19+rt */
				__asm	addpd	xmm7,xmm5		/*~t20=t20+it */

				__asm	mov	edi, isrt2
				__asm	movaps	xmm2,[esi+0x040]	/* t11 */
				__asm	movaps	xmm3,[esi+0x050]	/* t12 */
				__asm	subpd	xmm2,[esi+0x050]	/* t11-t12 */
				__asm	addpd	xmm3,[esi+0x040]	/* t12+t11 */
				__asm	mulpd	xmm2,[edi]	/* rt = (t11-t12)*ISRT2 */
				__asm	mulpd	xmm3,[edi]	/* it = (t12+t11)*ISRT2 */

				__asm	movaps	xmm0,[esi      ]	/* t3  */
				__asm	movaps	xmm1,[esi+0x010]	/* t4  */

				__asm	mov	edi, p4		/* restore p4-based value of edi */
				__asm	shl	edi, 3

				__asm	subpd	xmm0,xmm2			/*~t11=t3 -rt */
				__asm	subpd	xmm1,xmm3			/*~t12=t4 -it */
				__asm	addpd	xmm2,[esi      ]	/*~t3 =rt +t3 */
				__asm	addpd	xmm3,[esi+0x010]	/*~t4 =it +t4 */

				__asm	subpd	xmm2,xmm6		/* t3 -t19 */
				__asm	subpd	xmm3,xmm7		/* t4 -t20 */
				__asm	addpd	xmm6,xmm6		/*   2*t19 */
				__asm	addpd	xmm7,xmm7		/*   2*t20 */
				__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */
				__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */
				__asm	addpd	xmm6,xmm2		/* t3 +t19 */
				__asm	addpd	xmm7,xmm3		/* t4 +t20 */
				__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */
				__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */

				__asm	subpd	xmm0,xmm5		/* t11-t28 */
				__asm	subpd	xmm1,xmm4		/* t12-t27 */
				__asm	addpd	xmm5,xmm5		/*          2*t28 */
				__asm	addpd	xmm4,xmm4		/*          2*t27 */
				__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
				__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
				__asm	addpd	xmm5,xmm0		/* t11+t28 */
				__asm	addpd	xmm4,xmm1		/* t12+t27 */
				__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
				__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

			/*...Block 4: t7,15,23,31 ==> r18,1C,1A,1E	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
			#else
				add0 = &a[j1+p12];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;

				__asm	mov	eax, add0	/* restore main-array indices */
				__asm	mov	ebx, add1
				__asm	mov	ecx, add2
				__asm	mov	edx, add3
			#endif

				__asm	mov	esi, r18
				__asm	mov	edi, cc0
				__asm	movaps	xmm4,[esi+0x020]	/* t23 */		__asm	movaps	xmm6,[esi+0x060]	/* t31 */
				__asm	movaps	xmm5,[esi+0x030]	/* t24 */		__asm	movaps	xmm7,[esi+0x070]	/* t32 */
				__asm	movaps	xmm0,xmm4			/* copy t23 */	__asm	movaps	xmm2,xmm6			/* copy t31 */
				__asm	movaps	xmm1,xmm5			/* copy t24 */	__asm	movaps	xmm3,xmm7			/* copy t32 */

				__asm	mulpd	xmm4,[edi+0x10]	/* t23*s */			__asm	mulpd	xmm6,[edi     ]	/* t31*c */
				__asm	mulpd	xmm1,[edi     ]	/* t24*c */			__asm	mulpd	xmm3,[edi+0x10]	/* t32*s */
				__asm	mulpd	xmm5,[edi+0x10]	/* t24*s */			__asm	mulpd	xmm7,[edi     ]	/* t32*c */
				__asm	mulpd	xmm0,[edi     ]	/* t23*c */			__asm	mulpd	xmm2,[edi+0x10]	/* t31*s */
				__asm	subpd	xmm4,xmm1	/* ~t23 */				__asm	subpd	xmm6,xmm3	/* rt */
				__asm	addpd	xmm5,xmm0	/* ~t24 */				__asm	addpd	xmm7,xmm2	/* it */

				__asm	subpd	xmm4,xmm6		/*~t23=t23-rt */
				__asm	subpd	xmm5,xmm7		/*~t24=t24-it */
				__asm	addpd	xmm6,xmm6		/*      2* rt */
				__asm	addpd	xmm7,xmm7		/*      2* it */
				__asm	addpd	xmm6,xmm4		/*~t31=t23+rt */
				__asm	addpd	xmm7,xmm5		/*~t32=t24+it */

				__asm	mov	edi, isrt2
				__asm	movaps	xmm2,[esi+0x040]	/* t15 */
				__asm	movaps	xmm3,[esi+0x050]	/* t16 */
				__asm	addpd	xmm2,[esi+0x050]	/* t15+t16 */
				__asm	subpd	xmm3,[esi+0x040]	/* t16-t15 */
				__asm	mulpd	xmm2,[edi]	/* rt = (t15+t16)*ISRT2 */
				__asm	mulpd	xmm3,[edi]	/* it = (t16-t15)*ISRT2 */

				__asm	movaps	xmm0,[esi      ]	/* t7  */
				__asm	movaps	xmm1,[esi+0x010]	/* t8  */

				__asm	subpd	xmm0,xmm2			/*~t7 =t7 -rt */
				__asm	subpd	xmm1,xmm3			/*~t8 =t8 -it */
				__asm	addpd	xmm2,[esi      ]	/*~t15=rt +t7 */
				__asm	addpd	xmm3,[esi+0x010]	/*~t16=it +t8 */

				__asm	subpd	xmm0,xmm4		/* t7 -t23 */
				__asm	subpd	xmm1,xmm5		/* t8 -t24 */
				__asm	addpd	xmm4,xmm4		/*   2*t23 */
				__asm	addpd	xmm5,xmm5		/*   2*t24 */
				__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */
				__asm	movaps	[ebx+0x10],xmm1	/* a[jp+p1 ] */
				__asm	addpd	xmm4,xmm0		/* t7 +t23 */
				__asm	addpd	xmm5,xmm1		/* t8 +t24 */
				__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */
				__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */

				__asm	subpd	xmm2,xmm7		/* t15-t32 */
				__asm	subpd	xmm3,xmm6		/* t16-t31 */
				__asm	addpd	xmm7,xmm7		/*   2*t32 */
				__asm	addpd	xmm6,xmm6		/*   2*t31 */
				__asm	movaps	[ecx     ],xmm2	/* a[jt+p2 ] */
				__asm	movaps	[edx+0x10],xmm3	/* a[jp+p3 ] */
				__asm	addpd	xmm7,xmm2		/* t15+t32 */
				__asm	addpd	xmm6,xmm3		/* t16+t31 */
				__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
				__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

				/***************************************************/
				/* DIF Totals: 132 MOVapd, 182 ADD/SUBpd, 24 MULpd */
				/***************************************************/

		  #else	/* GCC-style inline ASM: */

			SSE2_RADIX16_DIF_NOTWIDDLE(add0,p1,p2,p3,p4,r00,r02,r04,r06,r08,r0A,r0C,r0E,r10,r12,r14,r16,r18,r1A,r1C,r1E,isrt2,cc0);

		  #endif

		#else	/* !USE_SSE2 */

			#if USE_SCALAR_DFT_MACRO
				RADIX_16_DIF(a1p0r,a1p0i,a1p1r,a1p1i,a1p2r,a1p2i,a1p3r,a1p3i,a1p4r,a1p4i,a1p5r,a1p5i,a1p6r,a1p6i,a1p7r,a1p7i,a1p8r,a1p8i,a1p9r,a1p9i,a1pAr,a1pAi,a1pBr,a1pBi,a1pCr,a1pCi,a1pDr,a1pDi,a1pEr,a1pEi,a1pFr,a1pFi
							,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
							,c,s)
			#else
				/*       gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms.	*/
			  #if PFETCH
				addr = &a[j1];
				prefetch_p_doubles(addr);
			  #endif
				/*...Block 1:	*/
				t3 =a1p0r -a1p8r;  	t1 =a1p0r +a1p8r;
				t4 =a1p0i -a1p8i;	t2 =a1p0i +a1p8i;

				t7 =a1p4r -a1pCr;	t5 =a1p4r +a1pCr;
				t8 =a1p4i -a1pCi;	t6 =a1p4i +a1pCi;
			  #if PFETCH
				addp = addr+p1;
				prefetch_p_doubles(addp);
			  #endif
				rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
				it =t6;	t6 =t2 -it;	t2 =t2 +it;

				rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
						t8 =t4 -rt;	t4 =t4 +rt;
			  #if PFETCH
				addp = addr+p2;
				prefetch_p_doubles(addp);
			  #endif
				/*...Block 2:	*/
				t11=a1p2r -a1pAr;	t9 =a1p2r +a1pAr;
				t12=a1p2i -a1pAi;	t10=a1p2i +a1pAi;

				t15=a1p6r -a1pEr;	t13=a1p6r +a1pEr;
				t16=a1p6i -a1pEi;	t14=a1p6i +a1pEi;
			  #if PFETCH
				addp = addr+p3;
				prefetch_p_doubles(addp);
			  #endif
				rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
				it =t14;	t14=t10-it;	t10=t10+it;

				rt =t15;	t15=t11+t16;t11=t11-t16;
							t16=t12-rt;	t12=t12+rt;
			  #if PFETCH
				addp = addr+p4;
				prefetch_p_doubles(addp);
			  #endif
				/*...Block 3:	*/
				t19=a1p1r -a1p9r;  	t17=a1p1r +a1p9r;
				t20=a1p1i -a1p9i;	t18=a1p1i +a1p9i;

				t23=a1p5r -a1pDr;	t21=a1p5r +a1pDr;
				t24=a1p5i -a1pDi;	t22=a1p5i +a1pDi;
			  #if PFETCH
				addp = addr+p5;
				prefetch_p_doubles(addp);
			  #endif
				rt =t21;	t21=t17-rt;	t17=t17+rt;
				it =t22;	t22=t18-it;	t18=t18+it;

				rt =t23;	t23=t19+t24;	t19=t19-t24;
							t24=t20-rt;	t20=t20+rt;
			  #if PFETCH
				addp = addr+p6;
				prefetch_p_doubles(addp);
			  #endif
				/*...Block 4:	*/
				t27=a1p3r -a1pBr;	t25=a1p3r +a1pBr;
				t28=a1p3i -a1pBi;	t26=a1p3i +a1pBi;

				t31=a1p7r -a1pFr;	t29=a1p7r +a1pFr;
				t32=a1p7i -a1pFi;	t30=a1p7i +a1pFi;
			  #if PFETCH
				addp = addr+p7;
				prefetch_p_doubles(addp);
			  #endif
				rt =t29;	t29=t25-rt;	t25=t25+rt;
				it =t30;	t30=t26-it;	t26=t26+it;

				rt =t31;	t31=t27+t32;t27=t27-t32;
							t32=t28-rt;	t28=t28+rt;
			  #if PFETCH
				addp = addr+p8;
				prefetch_p_doubles(addp);
			  #endif

				/*...and now do four more radix-4 transforms, including the internal twiddle factors:
				1, exp(i* 1*twopi/16) =       ( c, s), exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 3*twopi/16) =       ( s, c) (for inputs to transform block 2)
				1, exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 4*twopi/16) =       ( 0, 1), exp(i* 6*twopi/16) = ISRT2*(-1, 1) (for inputs to transform block 3)
				1, exp(i* 3*twopi/16) =       ( s, c), exp(i* 6*twopi/16) = ISRT2*(-1, 1), exp(i* 9*twopi/16) =       (-c,-s) (for inputs to transform block 4).
				(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
				I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
											 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
				and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

				/*...Block 1: t1,9,17,25	*/
				rt =t9 ;	t9 =t1 -rt;	t1 =t1 +rt;
				it =t10;	t10=t2 -it;	t2 =t2 +it;

				rt =t25;	t25=t17-rt;	t17=t17+rt;
				it =t26;	t26=t18-it;	t18=t18+it;
			  #if PFETCH
				addp = addr+p9;
				prefetch_p_doubles(addp);
			  #endif
				a[j1    ]=t1+t17;	a[j2    ]=t2+t18;
				a[j1+p1 ]=t1-t17;	a[j2+p1 ]=t2-t18;

				a[j1+p2 ]=t9 -t26;	a[j2+p2 ]=t10+t25;	/* mpy by E^4=i is inlined here...	*/
				a[j1+p3 ]=t9 +t26;	a[j2+p3 ]=t10-t25;
			  #if PFETCH
				addp = addr+p10;
				prefetch_p_doubles(addp);
			  #endif
				/*...Block 3: t5,13,21,29	*/
				rt =t13;	t13=t5 +t14;	t5 =t5 -t14;		/* twiddle mpy by E^4 = I	*/
				t14=t6 -rt;	t6 =t6 +rt;

				rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30	*/
				rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
				t29=t21+rt;		t21=t21-rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
				t30=t22+it;		t22=t22-it;
			  #if PFETCH
				addp = addr+p11;
				prefetch_p_doubles(addp);
			  #endif
				a[j1+p4 ]=t5+t21;	a[j2+p4 ]=t6+t22;
				a[j1+p5 ]=t5-t21;	a[j2+p5 ]=t6-t22;

				a[j1+p6 ]=t13-t30;	a[j2+p6 ]=t14+t29;	/* mpy by E^4=i is inlined here...	*/
				a[j1+p7 ]=t13+t30;	a[j2+p7 ]=t14-t29;
			  #if PFETCH
				addp = addr+p12;
				prefetch_p_doubles(addp);
			  #endif
				/*...Block 2: t3,11,19,27	*/
				rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;		/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12	*/
				t11=t3 -rt;		t3 =t3 +rt;
				t12=t4 -it;		t4 =t4 +it;

				rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;	/* twiddle mpy by E^1	*/
				rt =t27*s - t28*c;	it =t28*s + t27*c;		/* twiddle mpy by E^3	*/
				t27=t19-rt;		t19=t19+rt;
				t28=t20-it;		t20=t20+it;
			  #if PFETCH
				addp = addr+p13;
				prefetch_p_doubles(addp);
			  #endif
				a[j1+p8 ]=t3+t19;	a[j2+p8 ]=t4+t20;
				a[j1+p9 ]=t3-t19;	a[j2+p9 ]=t4-t20;

				a[j1+p10]=t11-t28;	a[j2+p10]=t12+t27;	/* mpy by E^4=i is inlined here...	*/
				a[j1+p11]=t11+t28;	a[j2+p11]=t12-t27;
			  #if PFETCH
				addp = addr+p14;
				prefetch_p_doubles(addp);
			  #endif
				/*...Block 4: t7,15,23,31	*/
				rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
				t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
				t16=t8 +it;		t8 =t8 -it;

				rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;	/* twiddle mpy by E^3	*/
				rt =t31*c - t32*s;	it =t32*c + t31*s;		/* twiddle mpy by E^1 = -E^9...	*/
				t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
				t32=t24+it;		t24=t24-it;			/* Note: t23+rt = t23*(s+1)	*/
			  #if PFETCH
				addp = addr+p15;
				prefetch_p_doubles(addp);
			  #endif
				a[j1+p12]=t7+t23;	a[j2+p12]=t8+t24;
				a[j1+p13]=t7-t23;	a[j2+p13]=t8-t24;

				a[j1+p14]=t15-t32;	a[j2+p14]=t16+t31;	/* mpy by E^4=i is inlined here...	*/
				a[j1+p15]=t15+t32;	a[j2+p15]=t16-t31;
			#endif	// USE_SCALAR_DFT_MACRO ?

		#endif	/* USE_SSE2 */

			}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += 16;
				co3 -= 16;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			_cy_r0[ithread] = cy_r0->d0;	_cy_r1[ithread] = cy_r0->d1;	_cy_r2[ithread] = cy_r0->d2;	_cy_r3[ithread] = cy_r0->d3;
			_cy_r4[ithread] = cy_r4->d0;	_cy_r5[ithread] = cy_r4->d1;	_cy_r6[ithread] = cy_r4->d2;	_cy_r7[ithread] = cy_r4->d3;
			_cy_r8[ithread] = cy_r8->d0;	_cy_r9[ithread] = cy_r8->d1;	_cy_rA[ithread] = cy_r8->d2;	_cy_rB[ithread] = cy_r8->d3;
			_cy_rC[ithread] = cy_rC->d0;	_cy_rD[ithread] = cy_rC->d1;	_cy_rE[ithread] = cy_rC->d2;	_cy_rF[ithread] = cy_rC->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			_cy_r0[ithread] = cy_r0->d0;	_cy_r1[ithread] = cy_r0->d1;
			_cy_r2[ithread] = cy_r2->d0;	_cy_r3[ithread] = cy_r2->d1;
			_cy_r4[ithread] = cy_r4->d0;	_cy_r5[ithread] = cy_r4->d1;
			_cy_r6[ithread] = cy_r6->d0;	_cy_r7[ithread] = cy_r6->d1;
			_cy_r8[ithread] = cy_r8->d0;	_cy_r9[ithread] = cy_r8->d1;
			_cy_rA[ithread] = cy_rA->d0;	_cy_rB[ithread] = cy_rA->d1;
			_cy_rC[ithread] = cy_rC->d0;	_cy_rD[ithread] = cy_rC->d1;
			_cy_rE[ithread] = cy_rE->d0;	_cy_rF[ithread] = cy_rE->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			_cy_r0[ithread] = cy_r0;
			_cy_r1[ithread] = cy_r1;
			_cy_r2[ithread] = cy_r2;
			_cy_r3[ithread] = cy_r3;
			_cy_r4[ithread] = cy_r4;
			_cy_r5[ithread] = cy_r5;
			_cy_r6[ithread] = cy_r6;
			_cy_r7[ithread] = cy_r7;
			_cy_r8[ithread] = cy_r8;
			_cy_r9[ithread] = cy_r9;
			_cy_rA[ithread] = cy_rA;
			_cy_rB[ithread] = cy_rB;
			_cy_rC[ithread] = cy_rC;
			_cy_rD[ithread] = cy_rD;
			_cy_rE[ithread] = cy_rE;
			_cy_rF[ithread] = cy_rF;
		#endif
		}
		else
		{
		#ifdef USE_AVX
			_cy_r0[ithread] = cy_r0->d0;	_cy_r1[ithread] = cy_r0->d1;	_cy_r2[ithread] = cy_r0->d2;	_cy_r3[ithread] = cy_r0->d3;
			_cy_r4[ithread] = cy_r4->d0;	_cy_r5[ithread] = cy_r4->d1;	_cy_r6[ithread] = cy_r4->d2;	_cy_r7[ithread] = cy_r4->d3;
			_cy_r8[ithread] = cy_r8->d0;	_cy_r9[ithread] = cy_r8->d1;	_cy_rA[ithread] = cy_r8->d2;	_cy_rB[ithread] = cy_r8->d3;
			_cy_rC[ithread] = cy_rC->d0;	_cy_rD[ithread] = cy_rC->d1;	_cy_rE[ithread] = cy_rC->d2;	_cy_rF[ithread] = cy_rC->d3;

			_cy_i0[ithread] = cy_i0->d0;	_cy_i1[ithread] = cy_i0->d1;	_cy_i2[ithread] = cy_i0->d2;	_cy_i3[ithread] = cy_i0->d3;
			_cy_i4[ithread] = cy_i4->d0;	_cy_i5[ithread] = cy_i4->d1;	_cy_i6[ithread] = cy_i4->d2;	_cy_i7[ithread] = cy_i4->d3;
			_cy_i8[ithread] = cy_i8->d0;	_cy_i9[ithread] = cy_i8->d1;	_cy_iA[ithread] = cy_i8->d2;	_cy_iB[ithread] = cy_i8->d3;
			_cy_iC[ithread] = cy_iC->d0;	_cy_iD[ithread] = cy_iC->d1;	_cy_iE[ithread] = cy_iC->d2;	_cy_iF[ithread] = cy_iC->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			_cy_r0[ithread] = cy_r0->d0;	_cy_i0[ithread] = cy_r0->d1;
			_cy_r1[ithread] = cy_r2->d0;	_cy_i1[ithread] = cy_r2->d1;
			_cy_r2[ithread] = cy_r4->d0;	_cy_i2[ithread] = cy_r4->d1;
			_cy_r3[ithread] = cy_r6->d0;	_cy_i3[ithread] = cy_r6->d1;
			_cy_r4[ithread] = cy_r8->d0;	_cy_i4[ithread] = cy_r8->d1;
			_cy_r5[ithread] = cy_rA->d0;	_cy_i5[ithread] = cy_rA->d1;
			_cy_r6[ithread] = cy_rC->d0;	_cy_i6[ithread] = cy_rC->d1;
			_cy_r7[ithread] = cy_rE->d0;	_cy_i7[ithread] = cy_rE->d1;
			_cy_r8[ithread] = cy_i0->d0;	_cy_i8[ithread] = cy_i0->d1;
			_cy_r9[ithread] = cy_i2->d0;	_cy_i9[ithread] = cy_i2->d1;
			_cy_rA[ithread] = cy_i4->d0;	_cy_iA[ithread] = cy_i4->d1;
			_cy_rB[ithread] = cy_i6->d0;	_cy_iB[ithread] = cy_i6->d1;
			_cy_rC[ithread] = cy_i8->d0;	_cy_iC[ithread] = cy_i8->d1;
			_cy_rD[ithread] = cy_iA->d0;	_cy_iD[ithread] = cy_iA->d1;
			_cy_rE[ithread] = cy_iC->d0;	_cy_iE[ithread] = cy_iC->d1;
			_cy_rF[ithread] = cy_iE->d0;	_cy_iF[ithread] = cy_iE->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			_cy_r0[ithread] = cy_r0;
			_cy_r1[ithread] = cy_r1;
			_cy_r2[ithread] = cy_r2;
			_cy_r3[ithread] = cy_r3;
			_cy_r4[ithread] = cy_r4;
			_cy_r5[ithread] = cy_r5;
			_cy_r6[ithread] = cy_r6;
			_cy_r7[ithread] = cy_r7;
			_cy_r8[ithread] = cy_r8;
			_cy_r9[ithread] = cy_r9;
			_cy_rA[ithread] = cy_rA;
			_cy_rB[ithread] = cy_rB;
			_cy_rC[ithread] = cy_rC;
			_cy_rD[ithread] = cy_rD;
			_cy_rE[ithread] = cy_rE;
			_cy_rF[ithread] = cy_rF;

			_cy_i0[ithread] = cy_i0;
			_cy_i1[ithread] = cy_i1;
			_cy_i2[ithread] = cy_i2;
			_cy_i3[ithread] = cy_i3;
			_cy_i4[ithread] = cy_i4;
			_cy_i5[ithread] = cy_i5;
			_cy_i6[ithread] = cy_i6;
			_cy_i7[ithread] = cy_i7;
			_cy_i8[ithread] = cy_i8;
			_cy_i9[ithread] = cy_i9;
			_cy_iA[ithread] = cy_iA;
			_cy_iB[ithread] = cy_iB;
			_cy_iC[ithread] = cy_iC;
			_cy_iD[ithread] = cy_iD;
			_cy_iE[ithread] = cy_iE;
			_cy_iF[ithread] = cy_iF;
		#endif
		}

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #ifdef OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(HERE, 0x0 == cy16_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix16_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			_cy_r0[ithread] = tdat[ithread].cy_r0;
			_cy_r1[ithread] = tdat[ithread].cy_r1;
			_cy_r2[ithread] = tdat[ithread].cy_r2;
			_cy_r3[ithread] = tdat[ithread].cy_r3;
			_cy_r4[ithread] = tdat[ithread].cy_r4;
			_cy_r5[ithread] = tdat[ithread].cy_r5;
			_cy_r6[ithread] = tdat[ithread].cy_r6;
			_cy_r7[ithread] = tdat[ithread].cy_r7;
			_cy_r8[ithread] = tdat[ithread].cy_r8;
			_cy_r9[ithread] = tdat[ithread].cy_r9;
			_cy_rA[ithread] = tdat[ithread].cy_rA;
			_cy_rB[ithread] = tdat[ithread].cy_rB;
			_cy_rC[ithread] = tdat[ithread].cy_rC;
			_cy_rD[ithread] = tdat[ithread].cy_rD;
			_cy_rE[ithread] = tdat[ithread].cy_rE;
			_cy_rF[ithread] = tdat[ithread].cy_rF;
		}
		else
		{
			_cy_r0[ithread] = tdat[ithread].cy_r0;	_cy_i0[ithread] = tdat[ithread].cy_i0;
			_cy_r1[ithread] = tdat[ithread].cy_r1;	_cy_i1[ithread] = tdat[ithread].cy_i1;
			_cy_r2[ithread] = tdat[ithread].cy_r2;	_cy_i2[ithread] = tdat[ithread].cy_i2;
			_cy_r3[ithread] = tdat[ithread].cy_r3;	_cy_i3[ithread] = tdat[ithread].cy_i3;
			_cy_r4[ithread] = tdat[ithread].cy_r4;	_cy_i4[ithread] = tdat[ithread].cy_i4;
			_cy_r5[ithread] = tdat[ithread].cy_r5;	_cy_i5[ithread] = tdat[ithread].cy_i5;
			_cy_r6[ithread] = tdat[ithread].cy_r6;	_cy_i6[ithread] = tdat[ithread].cy_i6;
			_cy_r7[ithread] = tdat[ithread].cy_r7;	_cy_i7[ithread] = tdat[ithread].cy_i7;
			_cy_r8[ithread] = tdat[ithread].cy_r8;	_cy_i8[ithread] = tdat[ithread].cy_i8;
			_cy_r9[ithread] = tdat[ithread].cy_r9;	_cy_i9[ithread] = tdat[ithread].cy_i9;
			_cy_rA[ithread] = tdat[ithread].cy_rA;	_cy_iA[ithread] = tdat[ithread].cy_iA;
			_cy_rB[ithread] = tdat[ithread].cy_rB;	_cy_iB[ithread] = tdat[ithread].cy_iB;
			_cy_rC[ithread] = tdat[ithread].cy_rC;	_cy_iC[ithread] = tdat[ithread].cy_iC;
			_cy_rD[ithread] = tdat[ithread].cy_rD;	_cy_iD[ithread] = tdat[ithread].cy_iD;
			_cy_rE[ithread] = tdat[ithread].cy_rE;	_cy_iE[ithread] = tdat[ithread].cy_iE;
			_cy_rF[ithread] = tdat[ithread].cy_rF;	_cy_iF[ithread] = tdat[ithread].cy_iF;
		}
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here:

	The cleanup carries from the end of each length-N/16 block into the begining of the next
	can all be neatly processed as follows:

	(1) Invert the radix-16 forward DIF FFT of the first block of 16 complex elements in A and unweight;
	(2) Propagate cleanup carries among the real and imaginary parts of the 16 outputs of (1);
	(3) Reweight and perform a radix-16 forward DIF FFT on the result of (2);
	(4) If any of the exit carries from (2) are nonzero, advance to the next 16 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1 = _cy_r0[CY_THREADS - 1];
		t3 = _cy_r1[CY_THREADS - 1];
		t5 = _cy_r2[CY_THREADS - 1];
		t7 = _cy_r3[CY_THREADS - 1];
		t9 = _cy_r4[CY_THREADS - 1];
		t11= _cy_r5[CY_THREADS - 1];
		t13= _cy_r6[CY_THREADS - 1];
		t15= _cy_r7[CY_THREADS - 1];
		t17= _cy_r8[CY_THREADS - 1];
		t19= _cy_r9[CY_THREADS - 1];
		t21= _cy_rA[CY_THREADS - 1];
		t23= _cy_rB[CY_THREADS - 1];
		t25= _cy_rC[CY_THREADS - 1];
		t27= _cy_rD[CY_THREADS - 1];
		t29= _cy_rE[CY_THREADS - 1];
		t31= _cy_rF[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"");	/* Make sure loop only gets executed if multiple threads */
			_cy_r0[ithread] = _cy_r0[ithread-1];
			_cy_r1[ithread] = _cy_r1[ithread-1];
			_cy_r2[ithread] = _cy_r2[ithread-1];
			_cy_r3[ithread] = _cy_r3[ithread-1];
			_cy_r4[ithread] = _cy_r4[ithread-1];
			_cy_r5[ithread] = _cy_r5[ithread-1];
			_cy_r6[ithread] = _cy_r6[ithread-1];
			_cy_r7[ithread] = _cy_r7[ithread-1];
			_cy_r8[ithread] = _cy_r8[ithread-1];
			_cy_r9[ithread] = _cy_r9[ithread-1];
			_cy_rA[ithread] = _cy_rA[ithread-1];
			_cy_rB[ithread] = _cy_rB[ithread-1];
			_cy_rC[ithread] = _cy_rC[ithread-1];
			_cy_rD[ithread] = _cy_rD[ithread-1];
			_cy_rE[ithread] = _cy_rE[ithread-1];
			_cy_rF[ithread] = _cy_rF[ithread-1];
		}

		_cy_r0[0] =+t31;	/* ...The wraparound carry is here: */
		_cy_r1[0] = t1 ;
		_cy_r2[0] = t3 ;
		_cy_r3[0] = t5 ;
		_cy_r4[0] = t7 ;
		_cy_r5[0] = t9 ;
		_cy_r6[0] = t11;
		_cy_r7[0] = t13;
		_cy_r8[0] = t15;
		_cy_r9[0] = t17;
		_cy_rA[0] = t19;
		_cy_rB[0] = t21;
		_cy_rC[0] = t23;
		_cy_rD[0] = t25;
		_cy_rE[0] = t27;
		_cy_rF[0] = t29;
	}
	else
	{
		t1 = _cy_r0[CY_THREADS - 1];	t2 = _cy_i0[CY_THREADS - 1];
		t3 = _cy_r1[CY_THREADS - 1];	t4 = _cy_i1[CY_THREADS - 1];
		t5 = _cy_r2[CY_THREADS - 1];	t6 = _cy_i2[CY_THREADS - 1];
		t7 = _cy_r3[CY_THREADS - 1];	t8 = _cy_i3[CY_THREADS - 1];
		t9 = _cy_r4[CY_THREADS - 1];	t10= _cy_i4[CY_THREADS - 1];
		t11= _cy_r5[CY_THREADS - 1];	t12= _cy_i5[CY_THREADS - 1];
		t13= _cy_r6[CY_THREADS - 1];	t14= _cy_i6[CY_THREADS - 1];
		t15= _cy_r7[CY_THREADS - 1];	t16= _cy_i7[CY_THREADS - 1];
		t17= _cy_r8[CY_THREADS - 1];	t18= _cy_i8[CY_THREADS - 1];
		t19= _cy_r9[CY_THREADS - 1];	t20= _cy_i9[CY_THREADS - 1];
		t21= _cy_rA[CY_THREADS - 1];	t22= _cy_iA[CY_THREADS - 1];
		t23= _cy_rB[CY_THREADS - 1];	t24= _cy_iB[CY_THREADS - 1];
		t25= _cy_rC[CY_THREADS - 1];	t26= _cy_iC[CY_THREADS - 1];
		t27= _cy_rD[CY_THREADS - 1];	t28= _cy_iD[CY_THREADS - 1];
		t29= _cy_rE[CY_THREADS - 1];	t30= _cy_iE[CY_THREADS - 1];
		t31= _cy_rF[CY_THREADS - 1];	t32= _cy_iF[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"");	/* Make sure loop only gets executed if multiple threads */
			_cy_r0[ithread] = _cy_r0[ithread-1];	_cy_i0[ithread] = _cy_i0[ithread-1];
			_cy_r1[ithread] = _cy_r1[ithread-1];	_cy_i1[ithread] = _cy_i1[ithread-1];
			_cy_r2[ithread] = _cy_r2[ithread-1];	_cy_i2[ithread] = _cy_i2[ithread-1];
			_cy_r3[ithread] = _cy_r3[ithread-1];	_cy_i3[ithread] = _cy_i3[ithread-1];
			_cy_r4[ithread] = _cy_r4[ithread-1];	_cy_i4[ithread] = _cy_i4[ithread-1];
			_cy_r5[ithread] = _cy_r5[ithread-1];	_cy_i5[ithread] = _cy_i5[ithread-1];
			_cy_r6[ithread] = _cy_r6[ithread-1];	_cy_i6[ithread] = _cy_i6[ithread-1];
			_cy_r7[ithread] = _cy_r7[ithread-1];	_cy_i7[ithread] = _cy_i7[ithread-1];
			_cy_r8[ithread] = _cy_r8[ithread-1];	_cy_i8[ithread] = _cy_i8[ithread-1];
			_cy_r9[ithread] = _cy_r9[ithread-1];	_cy_i9[ithread] = _cy_i9[ithread-1];
			_cy_rA[ithread] = _cy_rA[ithread-1];	_cy_iA[ithread] = _cy_iA[ithread-1];
			_cy_rB[ithread] = _cy_rB[ithread-1];	_cy_iB[ithread] = _cy_iB[ithread-1];
			_cy_rC[ithread] = _cy_rC[ithread-1];	_cy_iC[ithread] = _cy_iC[ithread-1];
			_cy_rD[ithread] = _cy_rD[ithread-1];	_cy_iD[ithread] = _cy_iD[ithread-1];
			_cy_rE[ithread] = _cy_rE[ithread-1];	_cy_iE[ithread] = _cy_iE[ithread-1];
			_cy_rF[ithread] = _cy_rF[ithread-1];	_cy_iF[ithread] = _cy_iF[ithread-1];
		}

		_cy_r0[0] =-t32;	_cy_i0[0] =+t31;	/* ...The 2 Mo"bius carries are here: */
		_cy_r1[0] = t1 ;	_cy_i1[0] = t2 ;
		_cy_r2[0] = t3 ;	_cy_i2[0] = t4 ;
		_cy_r3[0] = t5 ;	_cy_i3[0] = t6 ;
		_cy_r4[0] = t7 ;	_cy_i4[0] = t8 ;
		_cy_r5[0] = t9 ;	_cy_i5[0] = t10;
		_cy_r6[0] = t11;	_cy_i6[0] = t12;
		_cy_r7[0] = t13;	_cy_i7[0] = t14;
		_cy_r8[0] = t15;	_cy_i8[0] = t16;
		_cy_r9[0] = t17;	_cy_i9[0] = t18;
		_cy_rA[0] = t19;	_cy_iA[0] = t20;
		_cy_rB[0] = t21;	_cy_iB[0] = t22;
		_cy_rC[0] = t23;	_cy_iC[0] = t24;
		_cy_rD[0] = t25;	_cy_iD[0] = t26;
		_cy_rE[0] = t27;	_cy_iE[0] = t28;
		_cy_rF[0] = t29;	_cy_iF[0] = t30;
	}

	full_pass = 0;
	scale = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		j_jhi =15;
	}
	else
	{
		j_jhi = 7;
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			a[j    ] *= radix_inv;
			a[j+p1 ] *= radix_inv;
			a[j+p2 ] *= radix_inv;
			a[j+p3 ] *= radix_inv;
			a[j+p4 ] *= radix_inv;
			a[j+p5 ] *= radix_inv;
			a[j+p6 ] *= radix_inv;
			a[j+p7 ] *= radix_inv;
			a[j+p8 ] *= radix_inv;
			a[j+p9 ] *= radix_inv;
			a[j+p10] *= radix_inv;
			a[j+p11] *= radix_inv;
			a[j+p12] *= radix_inv;
			a[j+p13] *= radix_inv;
			a[j+p14] *= radix_inv;
			a[j+p15] *= radix_inv;
		}
	}
}	/* endfor(outer) */

	t1 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t1 += fabs(_cy_r0[0])+fabs(_cy_r1[0])+fabs(_cy_r2[0])+fabs(_cy_r3[0])+fabs(_cy_r4[0])+fabs(_cy_r5[0])+fabs(_cy_r6[0])+fabs(_cy_r7[0])+fabs(_cy_r8[0])+fabs(_cy_r9[0])+fabs(_cy_rA[0])+fabs(_cy_rB[0])+fabs(_cy_rC[0])+fabs(_cy_rD[0])+fabs(_cy_rE[0])+fabs(_cy_rF[0]);
		t1 += fabs(_cy_i0[0])+fabs(_cy_i1[0])+fabs(_cy_i2[0])+fabs(_cy_i3[0])+fabs(_cy_i4[0])+fabs(_cy_i5[0])+fabs(_cy_i6[0])+fabs(_cy_i7[0])+fabs(_cy_i8[0])+fabs(_cy_i9[0])+fabs(_cy_iA[0])+fabs(_cy_iB[0])+fabs(_cy_iC[0])+fabs(_cy_iD[0])+fabs(_cy_iE[0])+fabs(_cy_iF[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}
//fprintf(stderr, "radix16_carry: A[0-3] = %20.5f %20.5f %20.5f %20.5f\n",a[0],a[1],a[2],a[3]);
	if(t1 != 0.0)
	{
	    sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix16_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

/**************/

int radix16_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter                 , uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-16 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-16 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix16_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	static int n16;
	int i,j,j1,j2,jstart,jhi,full_pass,k1,k2,k,khi,l,ntmp,outer;
	int bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF;
	static uint32 bjmodnini;
	static uint64 psave=0;
	static uint32 bw,sw,nm1,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599, radix_inv,n2inv;
	double rt,it;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	double temp,scale
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi
		,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF
		,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF;

#if PFETCH
	double *addr, *addp;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double wt_re,wt_im;									/* Fermat-mod weights stuff */

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 nt_save = 0xffffffff, CY_THREADS = 0,pini;
	int ithread,j_jhi;
	static int **_bjmodn = 0x0, *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double **_cy_r = 0x0, **_cy_i = 0x0;
	/* Temporaries, used for storing each row-vector address during 2-D array allocs */
	int *iptr;
	double *dptr;

	// Mersenne-mod version of this routine broken-due-to-cause unknown (perhaps related to my carry-vector experiment here):
	ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE), "This function is known to have a bug ... please use radix16_ditN_cy_dif1-with-error-checking until fixed!");

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	bjmodn0=bjmodn1=bjmodn2=bjmodn3=bjmodn4=bjmodn5=bjmodn6=bjmodn7=bjmodn8=bjmodn9=bjmodnA=bjmodnB=bjmodnC=bjmodnD=bjmodnE=bjmodnF=-1;

/*...change n16 and n_div_wt to non-static to work around a gcc compiler bug. */
	n16   = n/16;
	n_div_nwt = n16 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n16)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/16.\n",iter);
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

	if(p != psave || nt_save != CY_THREADS)	/* Exponent or #thread change triggers re-init */
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)16));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		nm1   = n-1;

		/*   constant index offsets for load/stores are here.	*/
		p1 = n16 + ( (n16 >> DAT_BITS) << PAD_BITS );
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
		p12= p11+p1;
		p13= p12+p1;
		p14= p13+p1;
		p15= p14+p1;

		if(_cy_r)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			for(i=0; i < nt_save; i++)	/* Use the old thread count for the dealloc here */
			{
				free((void *)&_bjmodn[i][0]);
				free((void *)&_cy_r  [i][0]);
				free((void *)&_cy_i  [i][0]);
			}
			free((void *)_bjmodn); _bjmodn = 0x0;
			free((void *)_cy_r  ); _cy_r   = 0x0;
			free((void *)_cy_i  ); _cy_i   = 0x0;

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

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
			CY_THREADS = MAX_THREADS;

		ASSERT(HERE, CY_THREADS >= NTHREADS,"CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, n16      %CY_THREADS == 0,"n16      %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"n_div_nwt%CY_THREADS != 0");
		}

		nt_save = CY_THREADS;

	#ifdef MULTITHREAD
		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
	#endif

		pini = n16/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );

		_i       = (int *)malloc(CY_THREADS*sizeof(int));
		_jstart  = (int *)malloc(CY_THREADS*sizeof(int));
		_jhi     = (int *)malloc(CY_THREADS*sizeof(int));
		_col     = (int *)malloc(CY_THREADS*sizeof(int));
		_co2     = (int *)malloc(CY_THREADS*sizeof(int));
		_co3     = (int *)malloc(CY_THREADS*sizeof(int));

		_bjmodn  = (int    **)malloc(CY_THREADS*sizeof(void *));
		_cy_r    = (double **)malloc(CY_THREADS*sizeof(void *));
		_cy_i    = (double **)malloc(CY_THREADS*sizeof(void *));

		for(i=0; i < CY_THREADS; i++)
		{
			iptr = (   int *)malloc(16*sizeof(   int));
			_bjmodn[i] = iptr;
			dptr = (double *)malloc(16*sizeof(double));
			_cy_r  [i] = dptr;
			dptr = (double *)malloc(16*sizeof(double));
			_cy_i  [i] = dptr;
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
			i.e. the one that n2/16-separated FFT outputs need:
			*/
			_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int)); if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
			_bjmodnini[0] = 0;
			_bjmodnini[1] = 0;
			for(j=0; j < n16/CY_THREADS; j++)
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
			for(j=0; j < n16; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
			ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");
		}
	}	/* endif(first_entry) */

/*...The radix-16 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i=0; i < 16; i++)
		{
			_cy_r[ithread][i] = 0;	_cy_i[ithread][i] = 0;
		}
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
			_cy_r[      0][0] = -2;
	}

	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/

	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

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

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_bjmodn[ithread][0] = _bjmodnini[ithread];

			for(i=1; i < 16; i++)
			{
				_bjmodn[ithread][i] = _bjmodn[ithread][i-1] + _bjmodnini[CY_THREADS] - n; _bjmodn[ithread][i-1] = _bjmodn[ithread][i-1] + ( (-(int)((uint32)_bjmodn[ithread][i-1] >> 31)) & n);
			}

			_jstart[ithread] = ithread*n16/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*16);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+16 - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-16;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		khi = 1;

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*n16/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 15;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}
	}

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		col = _col[ithread];
		co2 = _co2[ithread];
		co3 = _co3[ithread];

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			bjmodn0 = _bjmodn[ithread][0x0];
			bjmodn1 = _bjmodn[ithread][0x1];
			bjmodn2 = _bjmodn[ithread][0x2];
			bjmodn3 = _bjmodn[ithread][0x3];
			bjmodn4 = _bjmodn[ithread][0x4];
			bjmodn5 = _bjmodn[ithread][0x5];
			bjmodn6 = _bjmodn[ithread][0x6];
			bjmodn7 = _bjmodn[ithread][0x7];
			bjmodn8 = _bjmodn[ithread][0x8];
			bjmodn9 = _bjmodn[ithread][0x9];
			bjmodnA = _bjmodn[ithread][0xA];
			bjmodnB = _bjmodn[ithread][0xB];
			bjmodnC = _bjmodn[ithread][0xC];
			bjmodnD = _bjmodn[ithread][0xD];
			bjmodnE = _bjmodn[ithread][0xE];
			bjmodnF = _bjmodn[ithread][0xF];
		}

		if(1)//full_pass)	/* Only do this on the main pass, not the cleanup-tails mini-pass: */
		{
			/* init carries	*/
			cy_r0 = _cy_r[ithread][0x0];	cy_i0 = _cy_i[ithread][0x0];
			cy_r1 = _cy_r[ithread][0x1];	cy_i1 = _cy_i[ithread][0x1];
			cy_r2 = _cy_r[ithread][0x2];	cy_i2 = _cy_i[ithread][0x2];
			cy_r3 = _cy_r[ithread][0x3];	cy_i3 = _cy_i[ithread][0x3];
			cy_r4 = _cy_r[ithread][0x4];	cy_i4 = _cy_i[ithread][0x4];
			cy_r5 = _cy_r[ithread][0x5];	cy_i5 = _cy_i[ithread][0x5];
			cy_r6 = _cy_r[ithread][0x6];	cy_i6 = _cy_i[ithread][0x6];
			cy_r7 = _cy_r[ithread][0x7];	cy_i7 = _cy_i[ithread][0x7];
			cy_r8 = _cy_r[ithread][0x8];	cy_i8 = _cy_i[ithread][0x8];
			cy_r9 = _cy_r[ithread][0x9];	cy_i9 = _cy_i[ithread][0x9];
			cy_rA = _cy_r[ithread][0xA];	cy_iA = _cy_i[ithread][0xA];
			cy_rB = _cy_r[ithread][0xB];	cy_iB = _cy_i[ithread][0xB];
			cy_rC = _cy_r[ithread][0xC];	cy_iC = _cy_i[ithread][0xC];
			cy_rD = _cy_r[ithread][0xD];	cy_iD = _cy_i[ithread][0xD];
			cy_rE = _cy_r[ithread][0xE];	cy_iE = _cy_i[ithread][0xE];
			cy_rF = _cy_r[ithread][0xF];	cy_iF = _cy_i[ithread][0xF];
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
		#ifdef USE_SSE2
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
				j1 = (j & mask01) + br4[j&3];
		#else
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
				j1 =  j;
		#endif
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

				/*...Block 1:	*/
				t1 =a[j1    ];	t2 =a[j2    ];
				rt =a[j1+p1 ];	it =a[j2+p1 ];
				t3 =t1 -rt;		t1 =t1 +rt;
				t4 =t2 -it;		t2 =t2 +it;

				t5 =a[j1+p2 ];	t6 =a[j2+p2 ];
				rt =a[j1+p3 ];	it =a[j2+p3 ];
				t7 =t5 -rt;  	t5 =t5 +rt;
				t8 =t6 -it;  	t6 =t6 +it;

				rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
				it =t6;	t6 =t2 -it;	t2 =t2 +it;

				rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
				t8 =t4 +rt;	t4 =t4 -rt;

				/*...Block 2:	*/
				t9 =a[j1+p4 ];	t10=a[j2+p4 ];
				rt =a[j1+p5 ];	it =a[j2+p5 ];
				t11=t9 -rt;		t9 =t9 +rt;
				t12=t10-it;		t10=t10+it;

				t13=a[j1+p6 ];	t14=a[j2+p6 ];
				rt =a[j1+p7 ];	it =a[j2+p7 ];
				t15=t13-rt;  	t13=t13+rt;
				t16=t14-it;		t14=t14+it;

				rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
				it =t14;	t14=t10-it;	t10=t10+it;

				rt =t15;	t15=t11-t16;	t11=t11+t16;
				t16=t12+rt;	t12=t12-rt;

				/*...Block 3:	*/
				t17=a[j1+p8 ];	t18=a[j2+p8 ];
				rt =a[j1+p9 ];	it =a[j2+p9 ];
				t19=t17-rt;		t17=t17+rt;
				t20=t18-it;		t18=t18+it;

				t21=a[j1+p10];	t22=a[j2+p10];
				rt =a[j1+p11];	it =a[j2+p11];
				t23=t21-rt;  	t21=t21+rt;
				t24=t22-it;		t22=t22+it;

				rt =t21;	t21=t17-rt;	t17=t17+rt;
				it =t22;	t22=t18-it;	t18=t18+it;

				rt =t23;	t23=t19-t24;	t19=t19+t24;
				t24=t20+rt;	t20=t20-rt;

				/*...Block 4:	*/
				t25=a[j1+p12];	t26=a[j2+p12];
				rt =a[j1+p13];	it =a[j2+p13];
				t27=t25-rt;		t25=t25+rt;
				t28=t26-it;		t26=t26+it;

				t29=a[j1+p14];	t30=a[j2+p14];
				rt =a[j1+p15];	it =a[j2+p15];
				t31=t29-rt;  	t29=t29+rt;
				t32=t30-it;		t30=t30+it;

				rt =t29;	t29=t25-rt;	t25=t25+rt;
				it =t30;	t30=t26-it;	t26=t26+it;

				rt =t31;	t31=t27-t32;	t27=t27+t32;
				t32=t28+rt;	t28=t28-rt;

				/*...and now do four more radix-4 transforms, including the internal twiddle factors:
				1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
				1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
				1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
				(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
				I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
											 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
				and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

				/*...Block 1: t1,9,17,25	*/
				rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
				it =t10;	t10=t2 -it;	t2 =t2 +it;

				rt =t25;	t25=t17-rt;	t17=t17+rt;
				it =t26;	t26=t18-it;	t18=t18+it;

				a1p0r =t1+t17;	a1p0i =t2+t18;
				a1p8r =t1-t17;	a1p8i =t2-t18;

				a1p4r =t9 +t26;	a1p4i =t10-t25;	/* mpy by E^-4 = -I is inlined here...	*/
				a1pCr=t9 -t26;	a1pCi=t10+t25;

				/*...Block 3: t5,13,21,29	*/
				rt =t13;	t13=t5 -t14;	t5 =t5 +t14;		/* twiddle mpy by E^4 =-I	*/
				t14=t6 +rt;	t6 =t6 -rt;

				rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;	/* twiddle mpy by E^-2	*/
				rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
				t29=t21+rt;		t21=t21-rt;			/* ...and get E^-6 by flipping signs here.	*/
				t30=t22+it;		t22=t22-it;

				a1p2r =t5+t21;	a1p2i =t6+t22;
				a1pAr=t5-t21;	a1pAi=t6-t22;

				a1p6r =t13+t30;	a1p6i =t14-t29;	/* mpy by E^-4 = -I is inlined here...	*/
				a1pEr=t13-t30;	a1pEi=t14+t29;

				/*...Block 2: t3,11,19,27	*/
				rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;		/* twiddle mpy by E^-2	*/
				t11=t3 -rt;		t3 =t3 +rt;
				t12=t4 -it;		t4 =t4 +it;

				rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;	/* twiddle mpy by E^-1	*/
				rt =t27*s + t28*c;	it =t28*s - t27*c;		/* twiddle mpy by E^-3	*/
				t27=t19-rt;		t19=t19+rt;
				t28=t20-it;		t20=t20+it;

				a1p1r =t3+t19;	a1p1i =t4+t20;
				a1p9r =t3-t19;	a1p9i =t4-t20;

				a1p5r =t11+t28;	a1p5i =t12-t27;	/* mpy by E^-4 = -I is inlined here...	*/
				a1pDr=t11-t28;	a1pDi=t12+t27;

				/*...Block 4: t7,15,23,31	*/
				rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
				t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
				t16=t8 +it;		t8 =t8 -it;

				rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;	/* twiddle mpy by E^-3	*/
				rt =t31*c + t32*s;	it =t32*c - t31*s;		/* twiddle mpy by E^-1 = -E^-9...	*/
				t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
				t32=t24+it;		t24=t24-it;

				a1p3r =t7+t23;	a1p3i =t8+t24;
				a1pBr=t7-t23;	a1pBi=t8-t24;

				a1p7r =t15+t32;	a1p7i =t16-t31;	/* mpy by E^-4 = -I is inlined here...	*/
				a1pFr=t15-t32;	a1pFi=t16+t31;

				/*...and combine those to complete the radix-16 transform and do the carries. Since the outputs would
				normally be getting dispatched to 16 separate blocks of the A-array, we need 16 separate carries.	*/

				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
					n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
					n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
					sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
					sinwtm1 = si[nwt-l-1];

					wtl     =wt0[    l  ];
					wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
					wtlp1   =wt0[    l+1];
					wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

					/*...set0 is slightly different from others:	*/
				   cmplx_carry_norm_pow2_nocheck0(a1p0r,a1p0i,cy_r0,bjmodn0    );
					cmplx_carry_norm_pow2_nocheck(a1p1r,a1p1i,cy_r1,bjmodn1,0x1);
					cmplx_carry_norm_pow2_nocheck(a1p2r,a1p2i,cy_r2,bjmodn2,0x2);
					cmplx_carry_norm_pow2_nocheck(a1p3r,a1p3i,cy_r3,bjmodn3,0x3);
					cmplx_carry_norm_pow2_nocheck(a1p4r,a1p4i,cy_r4,bjmodn4,0x4);
					cmplx_carry_norm_pow2_nocheck(a1p5r,a1p5i,cy_r5,bjmodn5,0x5);
					cmplx_carry_norm_pow2_nocheck(a1p6r,a1p6i,cy_r6,bjmodn6,0x6);
					cmplx_carry_norm_pow2_nocheck(a1p7r,a1p7i,cy_r7,bjmodn7,0x7);
					cmplx_carry_norm_pow2_nocheck(a1p8r,a1p8i,cy_r8,bjmodn8,0x8);
					cmplx_carry_norm_pow2_nocheck(a1p9r,a1p9i,cy_r9,bjmodn9,0x9);
					cmplx_carry_norm_pow2_nocheck(a1pAr,a1pAi,cy_rA,bjmodnA,0xA);
					cmplx_carry_norm_pow2_nocheck(a1pBr,a1pBi,cy_rB,bjmodnB,0xB);
					cmplx_carry_norm_pow2_nocheck(a1pCr,a1pCi,cy_rC,bjmodnC,0xC);
					cmplx_carry_norm_pow2_nocheck(a1pDr,a1pDi,cy_rD,bjmodnD,0xD);
					cmplx_carry_norm_pow2_nocheck(a1pEr,a1pEi,cy_rE,bjmodnE,0xE);
					cmplx_carry_norm_pow2_nocheck(a1pFr,a1pFi,cy_rF,bjmodnF,0xF);

					i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
					co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
				}
				else
				{
					ntmp = 0;
					fermat_carry_norm_pow2_nocheck(a1p0r,a1p0i,cy_r0,cy_i0,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p1r,a1p1i,cy_r1,cy_i1,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p2r,a1p2i,cy_r2,cy_i2,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p3r,a1p3i,cy_r3,cy_i3,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p4r,a1p4i,cy_r4,cy_i4,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p5r,a1p5i,cy_r5,cy_i5,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p6r,a1p6i,cy_r6,cy_i6,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p7r,a1p7i,cy_r7,cy_i7,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p8r,a1p8i,cy_r8,cy_i8,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p9r,a1p9i,cy_r9,cy_i9,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pAr,a1pAi,cy_rA,cy_iA,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pBr,a1pBi,cy_rB,cy_iB,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pCr,a1pCi,cy_rC,cy_iC,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pDr,a1pDi,cy_rD,cy_iD,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pEr,a1pEi,cy_rE,cy_iE,ntmp,NRTM1,NRT_BITS);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pFr,a1pFi,cy_rF,cy_iF,ntmp,NRTM1,NRT_BITS);
				}

				/*...The radix-16 DIF pass is here:	*/
				/*       gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms.	*/
			#if PFETCH
				addr = &a[j1];
				prefetch_p_doubles(addr);
			#endif
				/*...Block 1:	*/
				t3 =a1p0r -a1p8r;  	t1 =a1p0r +a1p8r;
				t4 =a1p0i -a1p8i;	t2 =a1p0i +a1p8i;

				t7 =a1p4r -a1pCr;	t5 =a1p4r +a1pCr;
				t8 =a1p4i -a1pCi;	t6 =a1p4i +a1pCi;
			#if PFETCH
				addp = addr+p1;
				prefetch_p_doubles(addp);
			#endif
				rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
				it =t6;	t6 =t2 -it;	t2 =t2 +it;

				rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
						t8 =t4 -rt;	t4 =t4 +rt;
			#if PFETCH
				addp = addr+p2;
				prefetch_p_doubles(addp);
			#endif
				/*...Block 2:	*/
				t11=a1p2r -a1pAr;	t9 =a1p2r +a1pAr;
				t12=a1p2i -a1pAi;	t10=a1p2i +a1pAi;

				t15=a1p6r -a1pEr;	t13=a1p6r +a1pEr;
				t16=a1p6i -a1pEi;	t14=a1p6i +a1pEi;
			#if PFETCH
				addp = addr+p3;
				prefetch_p_doubles(addp);
			#endif
				rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
				it =t14;	t14=t10-it;	t10=t10+it;

				rt =t15;	t15=t11+t16;t11=t11-t16;
							t16=t12-rt;	t12=t12+rt;
			#if PFETCH
				addp = addr+p4;
				prefetch_p_doubles(addp);
			#endif
				/*...Block 3:	*/
				t19=a1p1r -a1p9r;  	t17=a1p1r +a1p9r;
				t20=a1p1i -a1p9i;	t18=a1p1i +a1p9i;

				t23=a1p5r -a1pDr;	t21=a1p5r +a1pDr;
				t24=a1p5i -a1pDi;	t22=a1p5i +a1pDi;
			#if PFETCH
				addp = addr+p5;
				prefetch_p_doubles(addp);
			#endif
				rt =t21;	t21=t17-rt;	t17=t17+rt;
				it =t22;	t22=t18-it;	t18=t18+it;

				rt =t23;	t23=t19+t24;	t19=t19-t24;
							t24=t20-rt;	t20=t20+rt;
			#if PFETCH
				addp = addr+p6;
				prefetch_p_doubles(addp);
			#endif
				/*...Block 4:	*/
				t27=a1p3r -a1pBr;	t25=a1p3r +a1pBr;
				t28=a1p3i -a1pBi;	t26=a1p3i +a1pBi;

				t31=a1p7r -a1pFr;	t29=a1p7r +a1pFr;
				t32=a1p7i -a1pFi;	t30=a1p7i +a1pFi;
			#if PFETCH
				addp = addr+p7;
				prefetch_p_doubles(addp);
			#endif
				rt =t29;	t29=t25-rt;	t25=t25+rt;
				it =t30;	t30=t26-it;	t26=t26+it;

				rt =t31;	t31=t27+t32;t27=t27-t32;
							t32=t28-rt;	t28=t28+rt;
			#if PFETCH
				addp = addr+p8;
				prefetch_p_doubles(addp);
			#endif

				/*...and now do four more radix-4 transforms, including the internal twiddle factors:
				1, exp(i* 1*twopi/16) =       ( c, s), exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 3*twopi/16) =       ( s, c) (for inputs to transform block 2)
				1, exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 4*twopi/16) =       ( 0, 1), exp(i* 6*twopi/16) = ISRT2*(-1, 1) (for inputs to transform block 3)
				1, exp(i* 3*twopi/16) =       ( s, c), exp(i* 6*twopi/16) = ISRT2*(-1, 1), exp(i* 9*twopi/16) =       (-c,-s) (for inputs to transform block 4).
				(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
				I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
											 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
				and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

				/*...Block 1: t1,9,17,25	*/
				rt =t9 ;	t9 =t1 -rt;	t1 =t1 +rt;
				it =t10;	t10=t2 -it;	t2 =t2 +it;

				rt =t25;	t25=t17-rt;	t17=t17+rt;
				it =t26;	t26=t18-it;	t18=t18+it;
			#if PFETCH
				addp = addr+p9;
				prefetch_p_doubles(addp);
			#endif
				a[j1    ]=t1+t17;	a[j2    ]=t2+t18;
				a[j1+p1 ]=t1-t17;	a[j2+p1 ]=t2-t18;

				a[j1+p2 ]=t9 -t26;	a[j2+p2 ]=t10+t25;	/* mpy by E^4=i is inlined here...	*/
				a[j1+p3 ]=t9 +t26;	a[j2+p3 ]=t10-t25;
			#if PFETCH
				addp = addr+p10;
				prefetch_p_doubles(addp);
			#endif
				/*...Block 3: t5,13,21,29	*/
				rt =t13;	t13=t5 +t14;	t5 =t5 -t14;		/* twiddle mpy by E^4 = I	*/
				t14=t6 -rt;	t6 =t6 +rt;

				rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30	*/
				rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
				t29=t21+rt;		t21=t21-rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
				t30=t22+it;		t22=t22-it;
			#if PFETCH
				addp = addr+p11;
				prefetch_p_doubles(addp);
			#endif
				a[j1+p4 ]=t5+t21;	a[j2+p4 ]=t6+t22;
				a[j1+p5 ]=t5-t21;	a[j2+p5 ]=t6-t22;

				a[j1+p6 ]=t13-t30;	a[j2+p6 ]=t14+t29;	/* mpy by E^4=i is inlined here...	*/
				a[j1+p7 ]=t13+t30;	a[j2+p7 ]=t14-t29;
			#if PFETCH
				addp = addr+p12;
				prefetch_p_doubles(addp);
			#endif
				/*...Block 2: t3,11,19,27	*/
				rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;		/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12	*/
				t11=t3 -rt;		t3 =t3 +rt;
				t12=t4 -it;		t4 =t4 +it;

				rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;	/* twiddle mpy by E^1	*/
				rt =t27*s - t28*c;	it =t28*s + t27*c;		/* twiddle mpy by E^3	*/
				t27=t19-rt;		t19=t19+rt;
				t28=t20-it;		t20=t20+it;
			#if PFETCH
				addp = addr+p13;
				prefetch_p_doubles(addp);
			#endif
				a[j1+p8 ]=t3+t19;	a[j2+p8 ]=t4+t20;
				a[j1+p9 ]=t3-t19;	a[j2+p9 ]=t4-t20;

				a[j1+p10]=t11-t28;	a[j2+p10]=t12+t27;	/* mpy by E^4=i is inlined here...	*/
				a[j1+p11]=t11+t28;	a[j2+p11]=t12-t27;
			#if PFETCH
				addp = addr+p14;
				prefetch_p_doubles(addp);
			#endif
				/*...Block 4: t7,15,23,31	*/
				rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
				t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
				t16=t8 +it;		t8 =t8 -it;

				rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;	/* twiddle mpy by E^3	*/
				rt =t31*c - t32*s;	it =t32*c + t31*s;		/* twiddle mpy by E^1 = -E^9...	*/
				t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
				t32=t24+it;		t24=t24-it;			/* Note: t23+rt = t23*(s+1)	*/
			#if PFETCH
				addp = addr+p15;
				prefetch_p_doubles(addp);
			#endif
				a[j1+p12]=t7+t23;	a[j2+p12]=t8+t24;
				a[j1+p13]=t7-t23;	a[j2+p13]=t8-t24;

				a[j1+p14]=t15-t32;	a[j2+p14]=t16+t31;	/* mpy by E^4=i is inlined here...	*/
				a[j1+p15]=t15+t32;	a[j2+p15]=t16-t31;

			}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += 16;
				co3 -= 16;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		_cy_r[ithread][0x0] = cy_r0;	_cy_i[ithread][0x0] = cy_i0;
		_cy_r[ithread][0x1] = cy_r1;	_cy_i[ithread][0x1] = cy_i1;
		_cy_r[ithread][0x2] = cy_r2;	_cy_i[ithread][0x2] = cy_i2;
		_cy_r[ithread][0x3] = cy_r3;	_cy_i[ithread][0x3] = cy_i3;
		_cy_r[ithread][0x4] = cy_r4;	_cy_i[ithread][0x4] = cy_i4;
		_cy_r[ithread][0x5] = cy_r5;	_cy_i[ithread][0x5] = cy_i5;
		_cy_r[ithread][0x6] = cy_r6;	_cy_i[ithread][0x6] = cy_i6;
		_cy_r[ithread][0x7] = cy_r7;	_cy_i[ithread][0x7] = cy_i7;
		_cy_r[ithread][0x8] = cy_r8;	_cy_i[ithread][0x8] = cy_i8;
		_cy_r[ithread][0x9] = cy_r9;	_cy_i[ithread][0x9] = cy_i9;
		_cy_r[ithread][0xA] = cy_rA;	_cy_i[ithread][0xA] = cy_iA;
		_cy_r[ithread][0xB] = cy_rB;	_cy_i[ithread][0xB] = cy_iB;
		_cy_r[ithread][0xC] = cy_rC;	_cy_i[ithread][0xC] = cy_iC;
		_cy_r[ithread][0xD] = cy_rD;	_cy_i[ithread][0xD] = cy_iD;
		_cy_r[ithread][0xE] = cy_rE;	_cy_i[ithread][0xE] = cy_iE;
		_cy_r[ithread][0xF] = cy_rF;	_cy_i[ithread][0xF] = cy_iF;
	}	/******* END OF PARALLEL FOR-LOOP ********/

	if(!full_pass)break;

	/*   Wraparound carry cleanup loop is here:

	The cleanup carries from the end of each length-N/16 block into the begining of the next
	can all be neatly processed as follows:

	(1) Invert the radix-16 forward DIF FFT of the first block of 16 complex elements in A and unweight;
	(2) Propagate cleanup carries among the real and imaginary parts of the 16 outputs of (1);
	(3) Reweight and perform a radix-16 forward DIF FFT on the result of (2);
	(4) If any of the exit carries from (2) are nonzero, advance to the next 16 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1 = _cy_r[CY_THREADS - 1][0x0];
		t3 = _cy_r[CY_THREADS - 1][0x1];
		t5 = _cy_r[CY_THREADS - 1][0x2];
		t7 = _cy_r[CY_THREADS - 1][0x3];
		t9 = _cy_r[CY_THREADS - 1][0x4];
		t11= _cy_r[CY_THREADS - 1][0x5];
		t13= _cy_r[CY_THREADS - 1][0x6];
		t15= _cy_r[CY_THREADS - 1][0x7];
		t17= _cy_r[CY_THREADS - 1][0x8];
		t19= _cy_r[CY_THREADS - 1][0x9];
		t21= _cy_r[CY_THREADS - 1][0xA];
		t23= _cy_r[CY_THREADS - 1][0xB];
		t25= _cy_r[CY_THREADS - 1][0xC];
		t27= _cy_r[CY_THREADS - 1][0xD];
		t29= _cy_r[CY_THREADS - 1][0xE];
		t31= _cy_r[CY_THREADS - 1][0xF];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"");	/* Make sure loop only gets executed if multiple threads */
			for(i=0; i < 16; i++)
			{
				_cy_r[ithread][i] = _cy_r[ithread-1][i];
			}
		}

		_cy_r[0][0x0] =+t31;	/* ...The wraparound carry is here: */
		_cy_r[0][0x1] = t1 ;
		_cy_r[0][0x2] = t3 ;
		_cy_r[0][0x3] = t5 ;
		_cy_r[0][0x4] = t7 ;
		_cy_r[0][0x5] = t9 ;
		_cy_r[0][0x6] = t11;
		_cy_r[0][0x7] = t13;
		_cy_r[0][0x8] = t15;
		_cy_r[0][0x9] = t17;
		_cy_r[0][0xA] = t19;
		_cy_r[0][0xB] = t21;
		_cy_r[0][0xC] = t23;
		_cy_r[0][0xD] = t25;
		_cy_r[0][0xE] = t27;
		_cy_r[0][0xF] = t29;
	}
	else
	{
		t1 = _cy_r[CY_THREADS - 1][0x0];	t2 = _cy_i[CY_THREADS - 1][0x0];
		t3 = _cy_r[CY_THREADS - 1][0x1];	t4 = _cy_i[CY_THREADS - 1][0x1];
		t5 = _cy_r[CY_THREADS - 1][0x2];	t6 = _cy_i[CY_THREADS - 1][0x2];
		t7 = _cy_r[CY_THREADS - 1][0x3];	t8 = _cy_i[CY_THREADS - 1][0x3];
		t9 = _cy_r[CY_THREADS - 1][0x4];	t10= _cy_i[CY_THREADS - 1][0x4];
		t11= _cy_r[CY_THREADS - 1][0x5];	t12= _cy_i[CY_THREADS - 1][0x5];
		t13= _cy_r[CY_THREADS - 1][0x6];	t14= _cy_i[CY_THREADS - 1][0x6];
		t15= _cy_r[CY_THREADS - 1][0x7];	t16= _cy_i[CY_THREADS - 1][0x7];
		t17= _cy_r[CY_THREADS - 1][0x8];	t18= _cy_i[CY_THREADS - 1][0x8];
		t19= _cy_r[CY_THREADS - 1][0x9];	t20= _cy_i[CY_THREADS - 1][0x9];
		t21= _cy_r[CY_THREADS - 1][0xA];	t22= _cy_i[CY_THREADS - 1][0xA];
		t23= _cy_r[CY_THREADS - 1][0xB];	t24= _cy_i[CY_THREADS - 1][0xB];
		t25= _cy_r[CY_THREADS - 1][0xC];	t26= _cy_i[CY_THREADS - 1][0xC];
		t27= _cy_r[CY_THREADS - 1][0xD];	t28= _cy_i[CY_THREADS - 1][0xD];
		t29= _cy_r[CY_THREADS - 1][0xE];	t30= _cy_i[CY_THREADS - 1][0xE];
		t31= _cy_r[CY_THREADS - 1][0xF];	t32= _cy_i[CY_THREADS - 1][0xF];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"");	/* Make sure loop only gets executed if multiple threads */
			for(i=0; i < 16; i++)
			{
				_cy_r[ithread][i] = _cy_r[ithread-1][i-1];	_cy_i[ithread][i] = _cy_i[ithread-1][i-1];
			}
		}

		_cy_r[0][0x0] =-t32;	_cy_i[0][0x0] =+t31;	/* ...The 2 Mo"bius carries are here: */
		_cy_r[0][0x1] = t1 ;	_cy_i[0][0x1] = t2 ;
		_cy_r[0][0x2] = t3 ;	_cy_i[0][0x2] = t4 ;
		_cy_r[0][0x3] = t5 ;	_cy_i[0][0x3] = t6 ;
		_cy_r[0][0x4] = t7 ;	_cy_i[0][0x4] = t8 ;
		_cy_r[0][0x5] = t9 ;	_cy_i[0][0x5] = t10;
		_cy_r[0][0x6] = t11;	_cy_i[0][0x6] = t12;
		_cy_r[0][0x7] = t13;	_cy_i[0][0x7] = t14;
		_cy_r[0][0x8] = t15;	_cy_i[0][0x8] = t16;
		_cy_r[0][0x9] = t17;	_cy_i[0][0x9] = t18;
		_cy_r[0][0xA] = t19;	_cy_i[0][0xA] = t20;
		_cy_r[0][0xB] = t21;	_cy_i[0][0xB] = t22;
		_cy_r[0][0xC] = t23;	_cy_i[0][0xC] = t24;
		_cy_r[0][0xD] = t25;	_cy_i[0][0xD] = t26;
		_cy_r[0][0xE] = t27;	_cy_i[0][0xE] = t28;
		_cy_r[0][0xF] = t29;	_cy_i[0][0xF] = t30;
	}

	full_pass = 0;
	scale = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		j_jhi =15;
	}
	else
	{
		j_jhi = 7;
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			a[j    ] *= radix_inv;
			a[j+p1 ] *= radix_inv;
			a[j+p2 ] *= radix_inv;
			a[j+p3 ] *= radix_inv;
			a[j+p4 ] *= radix_inv;
			a[j+p5 ] *= radix_inv;
			a[j+p6 ] *= radix_inv;
			a[j+p7 ] *= radix_inv;
			a[j+p8 ] *= radix_inv;
			a[j+p9 ] *= radix_inv;
			a[j+p10] *= radix_inv;
			a[j+p11] *= radix_inv;
			a[j+p12] *= radix_inv;
			a[j+p13] *= radix_inv;
			a[j+p14] *= radix_inv;
			a[j+p15] *= radix_inv;
		}
	}
}	/* endfor(outer) */

	t1 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i=0; i < 16; i++)
		{
			t1 += fabs(_cy_r[ithread][i]);
			t1 += fabs(_cy_i[ithread][i]);
		}
	}

	if(t1 != 0.0)
	{
	    sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix16_ditN_cy_dif1_nochk - input wordsize may be too small.\n",iter);
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

#endif	/* #endifdef(GCD_STANDALONE) */

/***************/

void radix16_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-16 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n16,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, first_entry=TRUE;
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)] */
  #if !USE_SCALAR_DFT_MACRO
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
  #endif
	if(!first_entry && (n >> 4) != n16)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n16=n/16;

	  p1 = n16 + ( (n16 >> DAT_BITS) << PAD_BITS );
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
	  p12= p11+p1;
	  p13= p12+p1;
	  p14= p13+p1;
	  p15= p14+p1;
	}

/*...The radix-16 pass is here.	*/

	for(j = 0; j < n16; j += 2)
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

	#if USE_SCALAR_DFT_MACRO
		RADIX_16_DIF(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,c,s)
		// Test FMA-based DIF macro, with twiddles set = unity:
//		RADIX_16_DIF_FMA(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
//						,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
//						,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,c,s)
	#else
		/*...Block 1:	*/
	    t1 =a[j1    ];	t2 =a[j2    ];
	    rt =a[j1+p8 ];	it =a[j2+p8 ];
	    t3 =t1 -rt;  	t1 =t1 +rt;
	    t4 =t2 -it;		t2 =t2 +it;

	    t5 =a[j1+p4 ];	t6 =a[j2+p4 ];
	    rt =a[j1+p12];	it =a[j2+p12];
	    t7 =t5 -rt;  	t5 =t5 +rt;
	    t8 =t6 -it;  	t6 =t6 +it;

	    rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
	    it =t6;	t6 =t2 -it;	t2 =t2 +it;

	    rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
			t8 =t4 -rt;	t4 =t4 +rt;

		/*...Block 2:	*/
	    t9 =a[j1+p2 ];	t10=a[j2+p2 ];
	    rt =a[j1+p10];	it =a[j2+p10];
	    t11=t9 -rt;  	t9 =t9 +rt;
	    t12=t10-it;		t10=t10+it;

	    t13=a[j1+p6 ];	t14=a[j2+p6 ];
	    rt =a[j1+p14];	it =a[j2+p14];
	    t15=t13-rt;  	t13=t13+rt;
	    t16=t14-it;		t14=t14+it;

	    rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
	    it =t14;	t14=t10-it;	t10=t10+it;

	    rt =t15;	t15=t11+t16;	t11=t11-t16;
			t16=t12-rt;	t12=t12+rt;

		/*...Block 3:	*/
	    t17=a[j1+p1 ];	t18=a[j2+p1 ];
	    rt =a[j1+p9 ];	it =a[j2+p9 ];
	    t19=t17-rt;  	t17=t17+rt;
	    t20=t18-it;		t18=t18+it;

	    t21=a[j1+p5 ];	t22=a[j2+p5 ];
	    rt =a[j1+p13];	it =a[j2+p13];
	    t23=t21-rt;  	t21=t21+rt;
	    t24=t22-it;		t22=t22+it;

	    rt =t21;	t21=t17-rt;	t17=t17+rt;
	    it =t22;	t22=t18-it;	t18=t18+it;

	    rt =t23;	t23=t19+t24;	t19=t19-t24;
			t24=t20-rt;	t20=t20+rt;

		/*...Block 4:	*/
	    t25=a[j1+p3 ];	t26=a[j2+p3 ];
	    rt =a[j1+p11];	it =a[j2+p11];
	    t27=t25-rt;  	t25=t25+rt;
	    t28=t26-it;		t26=t26+it;

	    t29=a[j1+p7 ];	t30=a[j2+p7 ];
	    rt =a[j1+p15];	it =a[j2+p15];
	    t31=t29-rt;  	t29=t29+rt;
	    t32=t30-it;		t30=t30+it;

	    rt =t29;	t29=t25-rt;	t25=t25+rt;
	    it =t30;	t30=t26-it;	t26=t26+it;

	    rt =t31;	t31=t27+t32;	t27=t27-t32;
			t32=t28-rt;	t28=t28+rt;

		/*...and now do four more radix-4 transforms, including the internal twiddle factors:
			1, exp(i* 1*twopi/16) =       ( c, s), exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 3*twopi/16) =       ( s, c) (for inputs to transform block 2)
			1, exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 4*twopi/16) =       ( 0, 1), exp(i* 6*twopi/16) = ISRT2*(-1, 1) (for inputs to transform block 3)
			1, exp(i* 3*twopi/16) =       ( s, c), exp(i* 6*twopi/16) = ISRT2*(-1, 1), exp(i* 9*twopi/16) =       (-c,-s) (for inputs to transform block 4).
			(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
			 I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
													 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
			 and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

		/*...Block 1: t1,9,17,25	*/
	    rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
	    it =t10;	t10=t2 -it;	t2 =t2 +it;

	    rt =t25;	t25=t17-rt;	t17=t17+rt;
	    it =t26;	t26=t18-it;	t18=t18+it;

	    a[j1    ]=t1+t17;	a[j2    ]=t2+t18;
	    a[j1+p1 ]=t1-t17;	a[j2+p1 ]=t2-t18;

	    a[j1+p2 ]=t9 -t26;	a[j2+p2 ]=t10+t25;	/* mpy by E^4=i is inlined here...	*/
	    a[j1+p3 ]=t9 +t26;	a[j2+p3 ]=t10-t25;

		/*...Block 3: t5,13,21,29	*/
	    rt =t13;	t13=t5 +t14;	t5 =t5 -t14;		/* twiddle mpy by E^4 = I	*/
			t14=t6 -rt;	t6 =t6 +rt;

	    rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30	*/
	    rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
	    t29=t21+rt;		t21=t21-rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
	    t30=t22+it;		t22=t22-it;

	    a[j1+p4 ]=t5+t21;	a[j2+p4 ]=t6+t22;
	    a[j1+p5 ]=t5-t21;	a[j2+p5 ]=t6-t22;

	    a[j1+p6 ]=t13-t30;	a[j2+p6 ]=t14+t29;	/* mpy by E^4=i is inlined here...	*/
	    a[j1+p7 ]=t13+t30;	a[j2+p7 ]=t14-t29;

		/*...Block 2: t3,11,19,27	*/
	    rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;		/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12	*/
	    t11=t3 -rt;		t3 =t3 +rt;
	    t12=t4 -it;		t4 =t4 +it;

	    rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;	/* twiddle mpy by E^1	*/
	    rt =t27*s - t28*c;	it =t28*s + t27*c;		/* twiddle mpy by E^3	*/
	    t27=t19-rt;		t19=t19+rt;
	    t28=t20-it;		t20=t20+it;

	    a[j1+p8 ]=t3+t19;	a[j2+p8 ]=t4+t20;
	    a[j1+p9 ]=t3-t19;	a[j2+p9 ]=t4-t20;

	    a[j1+p10]=t11-t28;	a[j2+p10]=t12+t27;	/* mpy by E^4=i is inlined here...	*/
	    a[j1+p11]=t11+t28;	a[j2+p11]=t12-t27;

		/*...Block 4: t7,15,23,31	*/
	    rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
	    t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
	    t16=t8 +it;		t8 =t8 -it;

	    rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;	/* twiddle mpy by E^3	*/
	    rt =t31*c - t32*s;	it =t32*c + t31*s;		/* twiddle mpy by E^1 = -E^9...	*/
	    t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
	    t32=t24+it;		t24=t24-it;

		/* Note: t23+rt = t23*(s+1)	*/

	    a[j1+p12]=t7+t23;	a[j2+p12]=t8+t24;
	    a[j1+p13]=t7-t23;	a[j2+p13]=t8-t24;

	    a[j1+p14]=t15-t32;	a[j2+p14]=t16+t31;	/* mpy by E^4=i is inlined here...	*/
	    a[j1+p15]=t15+t32;	a[j2+p15]=t16-t31;
	#endif
	}
}

/**************/

void radix16_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-16 complex inverse DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details.
*/
	int j,j1,j2;
	static int n16,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, first_entry=TRUE;
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]	*/
  #if !USE_SCALAR_DFT_MACRO
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
  #endif
	if(!first_entry && (n >> 4) != n16)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n16=n/16;

	  p1 = n16 + ( (n16 >> DAT_BITS) << PAD_BITS );
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
	  p12= p11+p1;
	  p13= p12+p1;
	  p14= p13+p1;
	  p15= p14+p1;
	}

/*...The radix-16 pass is here.	*/

	for(j = 0; j < n16; j += 2)
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

	#if USE_SCALAR_DFT_MACRO
		RADIX_16_DIT(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,c,s)
	#else
		/*...Block 1:	*/
	    t1 =a[j1    ];	t2 =a[j2    ];
	    rt =a[j1+p1 ];	it =a[j2+p1 ];
	    t3 =t1 -rt;  	t1 =t1 +rt;
	    t4 =t2 -it;		t2 =t2 +it;

	    t5 =a[j1+p2 ];	t6 =a[j2+p2 ];
	    rt =a[j1+p3 ];	it =a[j2+p3 ];
	    t7 =t5 -rt;  	t5 =t5 +rt;
	    t8 =t6 -it;  	t6 =t6 +it;

	    rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
	    it =t6;	t6 =t2 -it;	t2 =t2 +it;

	    rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
			t8 =t4 +rt;	t4 =t4 -rt;

		/*...Block 2:	*/
	    t9 =a[j1+p4 ];	t10=a[j2+p4 ];
	    rt =a[j1+p5 ];	it =a[j2+p5 ];
	    t11=t9 -rt;  	t9 =t9 +rt;
	    t12=t10-it;		t10=t10+it;

	    t13=a[j1+p6 ];	t14=a[j2+p6 ];
	    rt =a[j1+p7 ];	it =a[j2+p7 ];
	    t15=t13-rt;  	t13=t13+rt;
	    t16=t14-it;		t14=t14+it;

	    rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
	    it =t14;	t14=t10-it;	t10=t10+it;

	    rt =t15;	t15=t11-t16;	t11=t11+t16;
			t16=t12+rt;	t12=t12-rt;

		/*...Block 3:	*/
	    t17=a[j1+p8 ];	t18=a[j2+p8 ];
	    rt =a[j1+p9 ];	it =a[j2+p9 ];
	    t19=t17-rt;  	t17=t17+rt;
	    t20=t18-it;		t18=t18+it;

	    t21=a[j1+p10];	t22=a[j2+p10];
	    rt =a[j1+p11];	it =a[j2+p11];
	    t23=t21-rt;  	t21=t21+rt;
	    t24=t22-it;		t22=t22+it;

	    rt =t21;	t21=t17-rt;	t17=t17+rt;
	    it =t22;	t22=t18-it;	t18=t18+it;

	    rt =t23;	t23=t19-t24;	t19=t19+t24;
			t24=t20+rt;	t20=t20-rt;

		/*...Block 4:	*/
	    t25=a[j1+p12];	t26=a[j2+p12];
	    rt =a[j1+p13];	it =a[j2+p13];
	    t27=t25-rt;  	t25=t25+rt;
	    t28=t26-it;		t26=t26+it;

	    t29=a[j1+p14];	t30=a[j2+p14];
	    rt =a[j1+p15];	it =a[j2+p15];
	    t31=t29-rt;  	t29=t29+rt;
	    t32=t30-it;		t30=t30+it;

	    rt =t29;	t29=t25-rt;	t25=t25+rt;
	    it =t30;	t30=t26-it;	t26=t26+it;

	    rt =t31;	t31=t27-t32;	t27=t27+t32;
			t32=t28+rt;	t28=t28-rt;

		/*...and now do four more radix-4 transforms, including the internal twiddle factors:
			1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
			1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
			1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
			(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
			 I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
													 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
			 and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

		/*...Block 1: t1,9,17,25	*/
	    rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
	    it =t10;	t10=t2 -it;	t2 =t2 +it;

	    rt =t25;	t25=t17-rt;	t17=t17+rt;
	    it =t26;	t26=t18-it;	t18=t18+it;

	    a[j1    ]=t1+t17;	a[j2    ]=t2+t18;
	    a[j1+p8 ]=t1-t17;	a[j2+p8 ]=t2-t18;

	    a[j1+p4 ]=t9 +t26;	a[j2+p4 ]=t10-t25;	/* mpy by E^-4 = -I is inlined here...	*/
	    a[j1+p12]=t9 -t26;	a[j2+p12]=t10+t25;

		/*...Block 3: t5,13,21,29	*/
	    rt =t13;	t13=t5 -t14;	t5 =t5 +t14;		/* twiddle mpy by E^4 =-I	*/
			t14=t6 +rt;	t6 =t6 -rt;

	    rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;	/* twiddle mpy by E^-2	*/
	    rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
	    t29=t21+rt;		t21=t21-rt;			/* ...and get E^-6 by flipping signs here.	*/
	    t30=t22+it;		t22=t22-it;

	    a[j1+p2 ]=t5+t21;	a[j2+p2 ]=t6+t22;
	    a[j1+p10]=t5-t21;	a[j2+p10]=t6-t22;

	    a[j1+p6 ]=t13+t30;	a[j2+p6 ]=t14-t29;	/* mpy by E^-4 = -I is inlined here...	*/
	    a[j1+p14]=t13-t30;	a[j2+p14]=t14+t29;

		/*...Block 2: t3,11,19,27	*/
	    rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;		/* twiddle mpy by E^-2	*/
	    t11=t3 -rt;		t3 =t3 +rt;
	    t12=t4 -it;		t4 =t4 +it;

	    rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;	/* twiddle mpy by E^-1	*/
	    rt =t27*s + t28*c;	it =t28*s - t27*c;		/* twiddle mpy by E^-3	*/
	    t27=t19-rt;		t19=t19+rt;
	    t28=t20-it;		t20=t20+it;

	    a[j1+p1 ]=t3+t19;	a[j2+p1 ]=t4+t20;
	    a[j1+p9 ]=t3-t19;	a[j2+p9 ]=t4-t20;

	    a[j1+p5 ]=t11+t28;	a[j2+p5 ]=t12-t27;	/* mpy by E^-4 = -I is inlined here...	*/
	    a[j1+p13]=t11-t28;	a[j2+p13]=t12+t27;

		/*...Block 4: t7,15,23,31	*/
	    rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
	    t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
	    t16=t8 +it;		t8 =t8 -it;

	    rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;	/* twiddle mpy by E^-3	*/
	    rt =t31*c + t32*s;	it =t32*c - t31*s;		/* twiddle mpy by E^-1 = -E^-9...	*/
	    t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
	    t32=t24+it;		t24=t24-it;

	    a[j1+p3 ]=t7+t23;	a[j2+p3 ]=t8+t24;
	    a[j1+p11]=t7-t23;	a[j2+p11]=t8-t24;

	    a[j1+p7 ]=t15+t32;	a[j2+p7 ]=t16-t31;	/* mpy by E^-4 = -I is inlined here...	*/
	    a[j1+p15]=t15-t32;	a[j2+p15]=t16+t31;
	#endif
	}
}

/******************** Multithreaded function body: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy16_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 16;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		int idx_offset,idx_incr;
		int j,j1,j2,k,l;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
		double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */

	#ifdef USE_SSE2

		uint32 p1,p2,p3,p4;
		double *add0, *add1, *add2, *add3;
		const double crnd = 3.0*0x4000000*0x2000000;
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *cc0, *ss0, *isrt2, *max_err, *sse2_rnd, *half_arr, *tmp,*tm2
			,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0A,*r0B,*r0C,*r0D,*r0E,*r0F
			,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1A,*r1B,*r1C,*r1D,*r1E,*r1F
			,*cy_r0,*cy_i0,*cy_r4,*cy_i4,*cy_r8,*cy_i8,*cy_rC,*cy_iC;
	  #ifndef USE_AVX
		vec_dbl
			 *cy_r2,*cy_i2,*cy_r6,*cy_i6,*cy_rA,*cy_iA,*cy_rE,*cy_iE;
	  #else
		int k1,k2;
		vec_dbl *base_negacyclic_root;
	  #endif
		int *bjmodn0,*bjmodn1,*bjmodn2,*bjmodn3,*bjmodn4,*bjmodn5,*bjmodn6,*bjmodn7,*bjmodn8,*bjmodn9,*bjmodnA,*bjmodnB,*bjmodnC,*bjmodnD,*bjmodnE,*bjmodnF;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
		double dtmp;

	#else

		const  double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;
		double *base, *baseinv;
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15;
		// Vars needed in scalar mode only:
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,k1,k2,m,m2,ntmp;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		double temp,frac
			,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr
			,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi
			,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF
			,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF;
		int bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF;

	#endif

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
		int ithread = thread_arg->tid;	/* unique thread index (use for debug) */

		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX, nm1 = n-1;
		int khi    = thread_arg->khi;
		int i      = thread_arg->i;	/* Pointer to the BASE and BASEINV arrays.	*/
		int jstart = thread_arg->jstart;
		int jhi    = thread_arg->jhi;
		int col = thread_arg->col;
		int co2 = thread_arg->co2;
		int co3 = thread_arg->co3;
		int sw  = thread_arg->sw, bw = n - sw;
		int nwt = thread_arg->nwt;

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

		p1  = NDIVR;
		p2  = p1  + p1;
		p3  = p2  + p1;
		p4  = p3  + p1;
	#ifndef USE_SSE2
		p5  = p4  + p1;
		p6  = p5  + p1;
		p7  = p6  + p1;
		p8  = p7  + p1;
		p9  = p8  + p1;
		p10 = p9  + p1;
		p11 = p10 + p1;
		p12 = p11 + p1;
		p13 = p12 + p1;
		p14 = p13 + p1;
		p15 = p14 + p1;
	#endif	// USE_SSE2
		p1  = p1  + ( (p1  >> DAT_BITS) << PAD_BITS );
		p2  = p2  + ( (p2  >> DAT_BITS) << PAD_BITS );
		p3  = p3  + ( (p3  >> DAT_BITS) << PAD_BITS );
		p4  = p4  + ( (p4  >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p5  = p5  + ( (p5  >> DAT_BITS) << PAD_BITS );
		p6  = p6  + ( (p6  >> DAT_BITS) << PAD_BITS );
		p7  = p7  + ( (p7  >> DAT_BITS) << PAD_BITS );
		p8  = p8  + ( (p8  >> DAT_BITS) << PAD_BITS );
		p9  = p9  + ( (p9  >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
	#endif	// USE_SSE2

	#ifdef USE_SSE2
		r00 = thread_arg->r00;	r01 = r00 + 0x01;
		r02 = r00 + 0x02;		r03 = r00 + 0x03;
		r04 = r00 + 0x04;		r05 = r00 + 0x05;
		r06 = r00 + 0x06;		r07 = r00 + 0x07;
		r08 = r00 + 0x08;		r09 = r00 + 0x09;
		r0A = r00 + 0x0a;		r0B = r00 + 0x0b;
		r0C = r00 + 0x0c;		r0D = r00 + 0x0d;
		r0E = r00 + 0x0e;		r0F = r00 + 0x0f;
		r10 = r00 + 0x10;		r11 = r00 + 0x11;
		r12 = r00 + 0x12;		r13 = r00 + 0x13;
		r14 = r00 + 0x14;		r15 = r00 + 0x15;
		r16 = r00 + 0x16;		r17 = r00 + 0x17;
		r18 = r00 + 0x18;		r19 = r00 + 0x19;
		r1A = r00 + 0x1a;		r1B = r00 + 0x1b;
		r1C = r00 + 0x1c;		r1D = r00 + 0x1d;
		r1E = r00 + 0x1e;		r1F = r00 + 0x1f;
	  isrt2 = r00 + 0x20;
		cc0 = r00 + 0x21;
		ss0 = r00 + 0x22;
		tmp = ss0 + 0x2;	// Only need +1 but prefer an even offset
	  #ifdef USE_AVX
		 cy_r0 = tmp + 0x00;
		 cy_r4 = tmp + 0x01;
		 cy_r8 = tmp + 0x02;
		 cy_rC = tmp + 0x03;
		 cy_i0 = tmp + 0x04;
		 cy_i4 = tmp + 0x05;
		 cy_i8 = tmp + 0x06;
		 cy_iC = tmp + 0x07;
		max_err = tmp + 0x08;
		sse2_rnd= tmp + 0x09;
		half_arr= tmp + 0x0a;	/* This table needs 20x16 bytes */
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
		 cy_r0 = tmp + 0x00;
		 cy_r2 = tmp + 0x01;
		 cy_r4 = tmp + 0x02;
		 cy_r6 = tmp + 0x03;
		 cy_r8 = tmp + 0x04;
		 cy_rA = tmp + 0x05;
		 cy_rC = tmp + 0x06;
		 cy_rE = tmp + 0x07;
		 cy_i0 = tmp + 0x08;
		 cy_i2 = tmp + 0x09;
		 cy_i4 = tmp + 0x0a;
		 cy_i6 = tmp + 0x0b;
		 cy_i8 = tmp + 0x0c;
		 cy_iA = tmp + 0x0d;
		 cy_iC = tmp + 0x0e;
		 cy_iE = tmp + 0x0f;
		max_err = tmp + 0x10;
		sse2_rnd= tmp + 0x11;
		half_arr= tmp + 0x12;	/* This table needs 20x16 bytes */
	  #endif

		ASSERT(HERE, (isrt2->d0 == ISRT2 && isrt2->d1 == ISRT2), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			tmp = half_arr;
		#ifdef USE_AVX
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#else	// SSE2:
			dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
		} else {
			dtmp = (half_arr)->d0 * (half_arr+1)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (half_arr)->d1 * (half_arr+1)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix16_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_nm1 = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
		sinwt         = (struct uint32x4 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;

		bjmodn0 = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn0 = (int*)(sse_nm1 + RE_IM_STRIDE);
	  #endif
		bjmodn1 = bjmodn0 + 0x01;
		bjmodn2 = bjmodn0 + 0x02;
		bjmodn3 = bjmodn0 + 0x03;
		bjmodn4 = bjmodn0 + 0x04;
		bjmodn5 = bjmodn0 + 0x05;
		bjmodn6 = bjmodn0 + 0x06;
		bjmodn7 = bjmodn0 + 0x07;
		bjmodn8 = bjmodn0 + 0x08;
		bjmodn9 = bjmodn0 + 0x09;
		bjmodnA = bjmodn0 + 0x0A;
		bjmodnB = bjmodn0 + 0x0B;
		bjmodnC = bjmodn0 + 0x0C;
		bjmodnD = bjmodn0 + 0x0D;
		bjmodnE = bjmodn0 + 0x0E;
		bjmodnF = bjmodn0 + 0x0F;

	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->r00     ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */	/* init carries	*/
		#ifdef USE_AVX
			*bjmodn0 = thread_arg->bjmodn0;		cy_r0->d0 = thread_arg->cy_r0;
			*bjmodn1 = thread_arg->bjmodn1;		cy_r0->d1 = thread_arg->cy_r1;
			*bjmodn2 = thread_arg->bjmodn2;		cy_r0->d2 = thread_arg->cy_r2;
			*bjmodn3 = thread_arg->bjmodn3;		cy_r0->d3 = thread_arg->cy_r3;
			*bjmodn4 = thread_arg->bjmodn4;		cy_r4->d0 = thread_arg->cy_r4;
			*bjmodn5 = thread_arg->bjmodn5;		cy_r4->d1 = thread_arg->cy_r5;
			*bjmodn6 = thread_arg->bjmodn6;		cy_r4->d2 = thread_arg->cy_r6;
			*bjmodn7 = thread_arg->bjmodn7;		cy_r4->d3 = thread_arg->cy_r7;
			*bjmodn8 = thread_arg->bjmodn8;		cy_r8->d0 = thread_arg->cy_r8;
			*bjmodn9 = thread_arg->bjmodn9;		cy_r8->d1 = thread_arg->cy_r9;
			*bjmodnA = thread_arg->bjmodnA;		cy_r8->d2 = thread_arg->cy_rA;
			*bjmodnB = thread_arg->bjmodnB;		cy_r8->d3 = thread_arg->cy_rB;
			*bjmodnC = thread_arg->bjmodnC;		cy_rC->d0 = thread_arg->cy_rC;
			*bjmodnD = thread_arg->bjmodnD;		cy_rC->d1 = thread_arg->cy_rD;
			*bjmodnE = thread_arg->bjmodnE;		cy_rC->d2 = thread_arg->cy_rE;
			*bjmodnF = thread_arg->bjmodnF;		cy_rC->d3 = thread_arg->cy_rF;
		#elif defined(USE_SSE2)
			*bjmodn0 = thread_arg->bjmodn0;		cy_r0->d0 = thread_arg->cy_r0;
			*bjmodn1 = thread_arg->bjmodn1;		cy_r0->d1 = thread_arg->cy_r1;
			*bjmodn2 = thread_arg->bjmodn2;		cy_r2->d0 = thread_arg->cy_r2;
			*bjmodn3 = thread_arg->bjmodn3;		cy_r2->d1 = thread_arg->cy_r3;
			*bjmodn4 = thread_arg->bjmodn4;		cy_r4->d0 = thread_arg->cy_r4;
			*bjmodn5 = thread_arg->bjmodn5;		cy_r4->d1 = thread_arg->cy_r5;
			*bjmodn6 = thread_arg->bjmodn6;		cy_r6->d0 = thread_arg->cy_r6;
			*bjmodn7 = thread_arg->bjmodn7;		cy_r6->d1 = thread_arg->cy_r7;
			*bjmodn8 = thread_arg->bjmodn8;		cy_r8->d0 = thread_arg->cy_r8;
			*bjmodn9 = thread_arg->bjmodn9;		cy_r8->d1 = thread_arg->cy_r9;
			*bjmodnA = thread_arg->bjmodnA;		cy_rA->d0 = thread_arg->cy_rA;
			*bjmodnB = thread_arg->bjmodnB;		cy_rA->d1 = thread_arg->cy_rB;
			*bjmodnC = thread_arg->bjmodnC;		cy_rC->d0 = thread_arg->cy_rC;
			*bjmodnD = thread_arg->bjmodnD;		cy_rC->d1 = thread_arg->cy_rD;
			*bjmodnE = thread_arg->bjmodnE;		cy_rE->d0 = thread_arg->cy_rE;
			*bjmodnF = thread_arg->bjmodnF;		cy_rE->d1 = thread_arg->cy_rF;
		#else
			bjmodn0 = thread_arg->bjmodn0;		cy_r0 = thread_arg->cy_r0;
			bjmodn1 = thread_arg->bjmodn1;		cy_r1 = thread_arg->cy_r1;
			bjmodn2 = thread_arg->bjmodn2;		cy_r2 = thread_arg->cy_r2;
			bjmodn3 = thread_arg->bjmodn3;		cy_r3 = thread_arg->cy_r3;
			bjmodn4 = thread_arg->bjmodn4;		cy_r4 = thread_arg->cy_r4;
			bjmodn5 = thread_arg->bjmodn5;		cy_r5 = thread_arg->cy_r5;
			bjmodn6 = thread_arg->bjmodn6;		cy_r6 = thread_arg->cy_r6;
			bjmodn7 = thread_arg->bjmodn7;		cy_r7 = thread_arg->cy_r7;
			bjmodn8 = thread_arg->bjmodn8;		cy_r8 = thread_arg->cy_r8;
			bjmodn9 = thread_arg->bjmodn9;		cy_r9 = thread_arg->cy_r9;
			bjmodnA = thread_arg->bjmodnA;		cy_rA = thread_arg->cy_rA;
			bjmodnB = thread_arg->bjmodnB;		cy_rB = thread_arg->cy_rB;
			bjmodnC = thread_arg->bjmodnC;		cy_rC = thread_arg->cy_rC;
			bjmodnD = thread_arg->bjmodnD;		cy_rD = thread_arg->cy_rD;
			bjmodnE = thread_arg->bjmodnE;		cy_rE = thread_arg->cy_rE;
			bjmodnF = thread_arg->bjmodnF;		cy_rF = thread_arg->cy_rF;
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#ifdef USE_AVX
			cy_r0->d0 = thread_arg->cy_r0;		cy_i0->d0 = thread_arg->cy_i0;
			cy_r0->d1 = thread_arg->cy_r1;		cy_i0->d1 = thread_arg->cy_i1;
			cy_r0->d2 = thread_arg->cy_r2;		cy_i0->d2 = thread_arg->cy_i2;
			cy_r0->d3 = thread_arg->cy_r3;		cy_i0->d3 = thread_arg->cy_i3;
			cy_r4->d0 = thread_arg->cy_r4;		cy_i4->d0 = thread_arg->cy_i4;
			cy_r4->d1 = thread_arg->cy_r5;		cy_i4->d1 = thread_arg->cy_i5;
			cy_r4->d2 = thread_arg->cy_r6;		cy_i4->d2 = thread_arg->cy_i6;
			cy_r4->d3 = thread_arg->cy_r7;		cy_i4->d3 = thread_arg->cy_i7;
			cy_r8->d0 = thread_arg->cy_r8;		cy_i8->d0 = thread_arg->cy_i8;
			cy_r8->d1 = thread_arg->cy_r9;		cy_i8->d1 = thread_arg->cy_i9;
			cy_r8->d2 = thread_arg->cy_rA;		cy_i8->d2 = thread_arg->cy_iA;
			cy_r8->d3 = thread_arg->cy_rB;		cy_i8->d3 = thread_arg->cy_iB;
			cy_rC->d0 = thread_arg->cy_rC;		cy_iC->d0 = thread_arg->cy_iC;
			cy_rC->d1 = thread_arg->cy_rD;		cy_iC->d1 = thread_arg->cy_iD;
			cy_rC->d2 = thread_arg->cy_rE;		cy_iC->d2 = thread_arg->cy_iE;
			cy_rC->d3 = thread_arg->cy_rF;		cy_iC->d3 = thread_arg->cy_iF;
		#elif defined(USE_SSE2)
			cy_r0->d0 = thread_arg->cy_r0;		cy_r0->d1 = thread_arg->cy_i0;
			cy_r2->d0 = thread_arg->cy_r1;		cy_r2->d1 = thread_arg->cy_i1;
			cy_r4->d0 = thread_arg->cy_r2;		cy_r4->d1 = thread_arg->cy_i2;
			cy_r6->d0 = thread_arg->cy_r3;		cy_r6->d1 = thread_arg->cy_i3;
			cy_r8->d0 = thread_arg->cy_r4;		cy_r8->d1 = thread_arg->cy_i4;
			cy_rA->d0 = thread_arg->cy_r5;		cy_rA->d1 = thread_arg->cy_i5;
			cy_rC->d0 = thread_arg->cy_r6;		cy_rC->d1 = thread_arg->cy_i6;
			cy_rE->d0 = thread_arg->cy_r7;		cy_rE->d1 = thread_arg->cy_i7;
			cy_i0->d0 = thread_arg->cy_r8;		cy_i0->d1 = thread_arg->cy_i8;
			cy_i2->d0 = thread_arg->cy_r9;		cy_i2->d1 = thread_arg->cy_i9;
			cy_i4->d0 = thread_arg->cy_rA;		cy_i4->d1 = thread_arg->cy_iA;
			cy_i6->d0 = thread_arg->cy_rB;		cy_i6->d1 = thread_arg->cy_iB;
			cy_i8->d0 = thread_arg->cy_rC;		cy_i8->d1 = thread_arg->cy_iC;
			cy_iA->d0 = thread_arg->cy_rD;		cy_iA->d1 = thread_arg->cy_iD;
			cy_iC->d0 = thread_arg->cy_rE;		cy_iC->d1 = thread_arg->cy_iE;
			cy_iE->d0 = thread_arg->cy_rF;		cy_iE->d1 = thread_arg->cy_iF;
		#else
			cy_r0 = thread_arg->cy_r0;			cy_i0 = thread_arg->cy_i0;
			cy_r1 = thread_arg->cy_r1;			cy_i1 = thread_arg->cy_i1;
			cy_r2 = thread_arg->cy_r2;			cy_i2 = thread_arg->cy_i2;
			cy_r3 = thread_arg->cy_r3;			cy_i3 = thread_arg->cy_i3;
			cy_r4 = thread_arg->cy_r4;			cy_i4 = thread_arg->cy_i4;
			cy_r5 = thread_arg->cy_r5;			cy_i5 = thread_arg->cy_i5;
			cy_r6 = thread_arg->cy_r6;			cy_i6 = thread_arg->cy_i6;
			cy_r7 = thread_arg->cy_r7;			cy_i7 = thread_arg->cy_i7;
			cy_r8 = thread_arg->cy_r8;			cy_i8 = thread_arg->cy_i8;
			cy_r9 = thread_arg->cy_r9;			cy_i9 = thread_arg->cy_i9;
			cy_rA = thread_arg->cy_rA;			cy_iA = thread_arg->cy_iA;
			cy_rB = thread_arg->cy_rB;			cy_iB = thread_arg->cy_iB;
			cy_rC = thread_arg->cy_rC;			cy_iC = thread_arg->cy_iC;
			cy_rD = thread_arg->cy_rD;			cy_iD = thread_arg->cy_iD;
			cy_rE = thread_arg->cy_rE;			cy_iE = thread_arg->cy_iE;
			cy_rF = thread_arg->cy_rF;			cy_iF = thread_arg->cy_iF;
		#endif
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += stride)
			{
				j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

			#ifdef USE_SSE2
				add0 = &a[j1];
				SSE2_RADIX16_DIT_NOTWIDDLE(add0,p1,p2,p3,p4,r00,r02,r04,r06,r08,r0a,r10,r18,isrt2,cc0);
			#else
				RADIX_16_DIT(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
							,a1p0r,a1p0i,a1p1r,a1p1i,a1p2r,a1p2i,a1p3r,a1p3i,a1p4r,a1p4i,a1p5r,a1p5i,a1p6r,a1p6i,a1p7r,a1p7i,a1p8r,a1p8i,a1p9r,a1p9i,a1pAr,a1pAi,a1pBr,a1pBi,a1pCr,a1pCi,a1pDr,a1pDi,a1pEr,a1pEi,a1pFr,a1pFi
							,c,s)
			#endif

			/*...Now do the carries. Since the outputs would
			normally be getting dispatched to 16 separate blocks of the A-array, we need 16 separate carries.	*/

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
			#ifdef USE_AVX

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

				l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
				n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				AVX_cmplx_carry_norm_pow2_errcheck0_X4(r00,add1,add2,add3,cy_r0,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(r08,add1,add2,add3,cy_r4,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(r10,add1,add2,add3,cy_r8,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(r18,add1,add2,add3,cy_rC,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

				i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/

			#elif defined(USE_SSE2)

				l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r00,add1,add2,cy_r0,cy_r2,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r08,add1,add2,cy_r4,cy_r6,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r10,add1,add2,cy_r8,cy_rA,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r18,add1,add2,cy_rC,cy_rE,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r00,add1,add2,cy_r0,cy_r2,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r08,add1,add2,cy_r4,cy_r6,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r10,add1,add2,cy_r8,cy_rA,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r18,add1,add2,cy_rC,cy_rE,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #endif

				i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/

			#else	// Scalar-double mode:

				l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

			/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_pow2_errcheck0(a1p0r,a1p0i,cy_r0,bjmodn0    );
				cmplx_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,bjmodn1,0x1);
				cmplx_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,bjmodn2,0x2);
				cmplx_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,bjmodn3,0x3);
				cmplx_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,bjmodn4,0x4);
				cmplx_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,bjmodn5,0x5);
				cmplx_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,bjmodn6,0x6);
				cmplx_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,bjmodn7,0x7);
				cmplx_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,bjmodn8,0x8);
				cmplx_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,bjmodn9,0x9);
				cmplx_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,bjmodnA,0xA);
				cmplx_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,bjmodnB,0xB);
				cmplx_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,bjmodnC,0xC);
				cmplx_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,bjmodnD,0xD);
				cmplx_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,bjmodnE,0xE);
				cmplx_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,bjmodnF,0xF);

				i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			#endif	// USE_AVX?

			}
			else	/* MODULUS_TYPE_FERMAT */
			{

			#ifdef USE_AVX

				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.
				tmp = half_arr+2;
				VEC_DBL_INIT(tmp, scale);
				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				tmp = base_negacyclic_root;	tm2 = tmp+1;

				// Hi-accuracy version needs 4 copies of each base root:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);

				tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r00,tmp,0x400,cy_r0,cy_i0,half_arr,sign_mask);
				tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r08,tmp,0x340,cy_r4,cy_i4,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r10,tmp,0x280,cy_r8,cy_i8,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r18,tmp,0x1c0,cy_rC,cy_iC,half_arr,sign_mask);

			#elif defined(USE_SSE2)

				tmp = half_arr+2;
				VEC_DBL_INIT(tmp, scale);
				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				SSE2_fermat_carry_norm_pow2_errcheck_X2(r00,cy_r0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r04,cy_r4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r08,cy_r8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r0C,cy_rC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r10,cy_i0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r14,cy_i4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r18,cy_i8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(r1C,cy_iC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);

			#else	// Scalar-double mode:

				ntmp = 0;
				fermat_carry_norm_pow2_errcheck(a1p0r,a1p0i,cy_r0,cy_i0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,cy_i1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,cy_i2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,cy_i3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,cy_i4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,cy_i5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,cy_i6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,cy_i7,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,cy_i8,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,cy_i9,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,cy_iA,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,cy_iB,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,cy_iC,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,cy_iD,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,cy_iE,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,cy_iF,ntmp,NRTM1,NRT_BITS);

			#endif	/* #ifdef USE_SSE2 */

			}	/* if(MODULUS_TYPE == ...) */

			/*...The radix-16 DIF pass is here:	*/
			#ifdef USE_SSE2
				SSE2_RADIX16_DIF_NOTWIDDLE(add0,p1,p2,p3,p4,r00,r02,r04,r06,r08,r0A,r0C,r0E,r10,r12,r14,r16,r18,r1A,r1C,r1E,isrt2,cc0);
			#else
				RADIX_16_DIF(a1p0r,a1p0i,a1p1r,a1p1i,a1p2r,a1p2i,a1p3r,a1p3i,a1p4r,a1p4i,a1p5r,a1p5i,a1p6r,a1p6i,a1p7r,a1p7i,a1p8r,a1p8i,a1p9r,a1p9i,a1pAr,a1pAi,a1pBr,a1pBi,a1pCr,a1pCi,a1pDr,a1pDi,a1pEr,a1pEi,a1pFr,a1pFi
							,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
							,c,s)
			#endif

			}	/* end for(j=_jstart; j < _jhi; j += 2) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += RADIX;
				co3 -= RADIX;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			thread_arg->cy_r0 = cy_r0->d0;
			thread_arg->cy_r1 = cy_r0->d1;
			thread_arg->cy_r2 = cy_r0->d2;
			thread_arg->cy_r3 = cy_r0->d3;
			thread_arg->cy_r4 = cy_r4->d0;
			thread_arg->cy_r5 = cy_r4->d1;
			thread_arg->cy_r6 = cy_r4->d2;
			thread_arg->cy_r7 = cy_r4->d3;
			thread_arg->cy_r8 = cy_r8->d0;
			thread_arg->cy_r9 = cy_r8->d1;
			thread_arg->cy_rA = cy_r8->d2;
			thread_arg->cy_rB = cy_r8->d3;
			thread_arg->cy_rC = cy_rC->d0;
			thread_arg->cy_rD = cy_rC->d1;
			thread_arg->cy_rE = cy_rC->d2;
			thread_arg->cy_rF = cy_rC->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			thread_arg->cy_r0 = cy_r0->d0;
			thread_arg->cy_r1 = cy_r0->d1;
			thread_arg->cy_r2 = cy_r2->d0;
			thread_arg->cy_r3 = cy_r2->d1;
			thread_arg->cy_r4 = cy_r4->d0;
			thread_arg->cy_r5 = cy_r4->d1;
			thread_arg->cy_r6 = cy_r6->d0;
			thread_arg->cy_r7 = cy_r6->d1;
			thread_arg->cy_r8 = cy_r8->d0;
			thread_arg->cy_r9 = cy_r8->d1;
			thread_arg->cy_rA = cy_rA->d0;
			thread_arg->cy_rB = cy_rA->d1;
			thread_arg->cy_rC = cy_rC->d0;
			thread_arg->cy_rD = cy_rC->d1;
			thread_arg->cy_rE = cy_rE->d0;
			thread_arg->cy_rF = cy_rE->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			thread_arg->cy_r0 = cy_r0;
			thread_arg->cy_r1 = cy_r1;
			thread_arg->cy_r2 = cy_r2;
			thread_arg->cy_r3 = cy_r3;
			thread_arg->cy_r4 = cy_r4;
			thread_arg->cy_r5 = cy_r5;
			thread_arg->cy_r6 = cy_r6;
			thread_arg->cy_r7 = cy_r7;
			thread_arg->cy_r8 = cy_r8;
			thread_arg->cy_r9 = cy_r9;
			thread_arg->cy_rA = cy_rA;
			thread_arg->cy_rB = cy_rB;
			thread_arg->cy_rC = cy_rC;
			thread_arg->cy_rD = cy_rD;
			thread_arg->cy_rE = cy_rE;
			thread_arg->cy_rF = cy_rF;
		#endif	// SSE2 or AVX?
		}
		else
		{
		#ifdef USE_AVX
			thread_arg->cy_r0 = cy_r0->d0;		thread_arg->cy_i0 = cy_i0->d0;
			thread_arg->cy_r1 = cy_r0->d1;		thread_arg->cy_i1 = cy_i0->d1;
			thread_arg->cy_r2 = cy_r0->d2;		thread_arg->cy_i2 = cy_i0->d2;
			thread_arg->cy_r3 = cy_r0->d3;		thread_arg->cy_i3 = cy_i0->d3;
			thread_arg->cy_r4 = cy_r4->d0;		thread_arg->cy_i4 = cy_i4->d0;
			thread_arg->cy_r5 = cy_r4->d1;		thread_arg->cy_i5 = cy_i4->d1;
			thread_arg->cy_r6 = cy_r4->d2;		thread_arg->cy_i6 = cy_i4->d2;
			thread_arg->cy_r7 = cy_r4->d3;		thread_arg->cy_i7 = cy_i4->d3;
			thread_arg->cy_r8 = cy_r8->d0;		thread_arg->cy_i8 = cy_i8->d0;
			thread_arg->cy_r9 = cy_r8->d1;		thread_arg->cy_i9 = cy_i8->d1;
			thread_arg->cy_rA = cy_r8->d2;		thread_arg->cy_iA = cy_i8->d2;
			thread_arg->cy_rB = cy_r8->d3;		thread_arg->cy_iB = cy_i8->d3;
			thread_arg->cy_rC = cy_rC->d0;		thread_arg->cy_iC = cy_iC->d0;
			thread_arg->cy_rD = cy_rC->d1;		thread_arg->cy_iD = cy_iC->d1;
			thread_arg->cy_rE = cy_rC->d2;		thread_arg->cy_iE = cy_iC->d2;
			thread_arg->cy_rF = cy_rC->d3;		thread_arg->cy_iF = cy_iC->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			thread_arg->cy_r0 = cy_r0->d0;		thread_arg->cy_i0 = cy_r0->d1;
			thread_arg->cy_r1 = cy_r2->d0;		thread_arg->cy_i1 = cy_r2->d1;
			thread_arg->cy_r2 = cy_r4->d0;		thread_arg->cy_i2 = cy_r4->d1;
			thread_arg->cy_r3 = cy_r6->d0;		thread_arg->cy_i3 = cy_r6->d1;
			thread_arg->cy_r4 = cy_r8->d0;		thread_arg->cy_i4 = cy_r8->d1;
			thread_arg->cy_r5 = cy_rA->d0;		thread_arg->cy_i5 = cy_rA->d1;
			thread_arg->cy_r6 = cy_rC->d0;		thread_arg->cy_i6 = cy_rC->d1;
			thread_arg->cy_r7 = cy_rE->d0;		thread_arg->cy_i7 = cy_rE->d1;
			thread_arg->cy_r8 = cy_i0->d0;		thread_arg->cy_i8 = cy_i0->d1;
			thread_arg->cy_r9 = cy_i2->d0;		thread_arg->cy_i9 = cy_i2->d1;
			thread_arg->cy_rA = cy_i4->d0;		thread_arg->cy_iA = cy_i4->d1;
			thread_arg->cy_rB = cy_i6->d0;		thread_arg->cy_iB = cy_i6->d1;
			thread_arg->cy_rC = cy_i8->d0;		thread_arg->cy_iC = cy_i8->d1;
			thread_arg->cy_rD = cy_iA->d0;		thread_arg->cy_iD = cy_iA->d1;
			thread_arg->cy_rE = cy_iC->d0;		thread_arg->cy_iE = cy_iC->d1;
			thread_arg->cy_rF = cy_iE->d0;		thread_arg->cy_iF = cy_iE->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			thread_arg->cy_r0 = cy_r0;			thread_arg->cy_i0 = cy_i0;
			thread_arg->cy_r1 = cy_r1;			thread_arg->cy_i1 = cy_i1;
			thread_arg->cy_r2 = cy_r2;			thread_arg->cy_i2 = cy_i2;
			thread_arg->cy_r3 = cy_r3;			thread_arg->cy_i3 = cy_i3;
			thread_arg->cy_r4 = cy_r4;			thread_arg->cy_i4 = cy_i4;
			thread_arg->cy_r5 = cy_r5;			thread_arg->cy_i5 = cy_i5;
			thread_arg->cy_r6 = cy_r6;			thread_arg->cy_i6 = cy_i6;
			thread_arg->cy_r7 = cy_r7;			thread_arg->cy_i7 = cy_i7;
			thread_arg->cy_r8 = cy_r8;			thread_arg->cy_i8 = cy_i8;
			thread_arg->cy_r9 = cy_r9;			thread_arg->cy_i9 = cy_i9;
			thread_arg->cy_rA = cy_rA;			thread_arg->cy_iA = cy_iA;
			thread_arg->cy_rB = cy_rB;			thread_arg->cy_iB = cy_iB;
			thread_arg->cy_rC = cy_rC;			thread_arg->cy_iC = cy_iC;
			thread_arg->cy_rD = cy_rD;			thread_arg->cy_iD = cy_iD;
			thread_arg->cy_rE = cy_rE;			thread_arg->cy_iE = cy_iE;
			thread_arg->cy_rF = cy_rF;			thread_arg->cy_iF = cy_iF;
		#endif	// SSE2 or AVX?
		}

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}

		return 0x0;
	}
#endif

