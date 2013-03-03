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

#undef FFT_DEBUG
#define FFT_DEBUG	0

#if FFT_DEBUG
	#include "carry_dbg.h"
#endif

#ifdef USE_SSE2

	const int radix16_creals_in_local_store = 80;

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

	#include "sse2_macro.h"

//	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */

	#ifdef COMPILER_TYPE_MSVC
		/*  */
	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix16_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix16_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

  #ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
		int tid;
		int ndivr;
		int _pad0;	// Pads to make sizeof this struct a multiple of 16 bytes
		int _pad1;
	
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
		struct complex *r1;

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

#endif	/* USE_SSE2 */

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
	const double crnd = 3.0*0x4000000*0x2000000;
	int i,j,j1,j2,jstart,jhi,full_pass,k,khi,outer;
	int l,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	static uint32 bjmodnini;
	static uint64 psave=0;
	static uint32 bw,sw,nm1,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15;
	const  double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;
	static double radix_inv,n2inv,scale;	/* Need scale to be static since in the MSVC-only pure-ASM vesrion of the carry step save the address in a static pointer below */
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	double maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;

#ifdef USE_SSE2

	static int cslots_in_local_store;
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
	uint64 tmp64;

  #ifdef MULTITHREAD

	#ifdef USE_PTHREAD
		static struct complex *__r0;	// Base address for discrete per-thread local stores
		static struct cy_thread_data_t *tdat = 0x0;
		// Threadpool-based dispatch stuff:
		static int main_work_units = 0, pool_work_units = 0;
		static struct threadpool *tpool = 0x0;
		static int task_is_blocking = TRUE;
		static thread_control_t thread_control = {0,0,0};
		// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
		static task_control_t   task_control = {NULL, (void*)cy16_process_chunk, NULL, 0x0};
	#endif

  #else

//	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	static int idx_offset, idx_incr;
//	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
/*...stuff for the reduced-length DWT weights array is here:	*/
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */

  #endif

	static int *si_ptr;
	static double *wt0_ptr, *wt1_ptr, *scale_ptr = &scale;
	static struct complex *cc0, *ss0, *isrt2, *max_err, *sse2_rnd, *half_arr, *tmp;
	static struct complex *r1,*r2,*r3,*r4,*r5,*r6,*r7,*r8,*r9,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r30,*r31,*r32;
	static struct complex *cy_r01,*cy_r23,*cy_r45,*cy_r67,*cy_r89,*cy_rAB,*cy_rCD,*cy_rEF,*cy_i01,*cy_i23,*cy_i45,*cy_i67,*cy_i89,*cy_iAB,*cy_iCD,*cy_iEF;
	static int *bjmodn0,*bjmodn1,*bjmodn2,*bjmodn3,*bjmodn4,*bjmodn5,*bjmodn6,*bjmodn7,*bjmodn8,*bjmodn9,*bjmodnA,*bjmodnB,*bjmodnC,*bjmodnD,*bjmodnE,*bjmodnF;

#else

	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double rt,it;
	int jt,jp,k1,k2,m,m2,ntmp;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	static double wt_re,wt_im;									/* Fermat-mod weights stuff */
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
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/16 in radix16_ditN_cy_dif1.\n",iter);
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

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 128 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix16_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_COMPLEX(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
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
		r1  = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;
		r2  = sc_ptr + 0x01;		cc0 = sc_ptr + 0x21;
		r3  = sc_ptr + 0x02;		ss0 = sc_ptr + 0x22;
		r4  = sc_ptr + 0x03;	 cy_r01 = sc_ptr + 0x23;
		r5  = sc_ptr + 0x04;	 cy_r23 = sc_ptr + 0x24;
		r6  = sc_ptr + 0x05;	 cy_r45 = sc_ptr + 0x25;
		r7  = sc_ptr + 0x06;	 cy_r67 = sc_ptr + 0x26;
		r8  = sc_ptr + 0x07;	 cy_r89 = sc_ptr + 0x27;
		r9  = sc_ptr + 0x08;	 cy_rAB = sc_ptr + 0x28;
		r10 = sc_ptr + 0x09;	 cy_rCD = sc_ptr + 0x29;
		r11 = sc_ptr + 0x0a;	 cy_rEF = sc_ptr + 0x2a;
		r12 = sc_ptr + 0x0b;	 cy_i01 = sc_ptr + 0x2b;
		r13 = sc_ptr + 0x0c;	 cy_i23 = sc_ptr + 0x2c;
		r14 = sc_ptr + 0x0d;	 cy_i45 = sc_ptr + 0x2d;
		r15 = sc_ptr + 0x0e;	 cy_i67 = sc_ptr + 0x2e;
		r16 = sc_ptr + 0x0f;	 cy_i89 = sc_ptr + 0x2f;
		r17 = sc_ptr + 0x10;	 cy_iAB = sc_ptr + 0x30;
		r18 = sc_ptr + 0x11;	 cy_iCD = sc_ptr + 0x31;
		r19 = sc_ptr + 0x12;	 cy_iEF = sc_ptr + 0x32;
		r20 = sc_ptr + 0x13;	max_err = sc_ptr + 0x33;
		r21 = sc_ptr + 0x14;	sse2_rnd= sc_ptr + 0x34;
		r22 = sc_ptr + 0x15;	half_arr= sc_ptr + 0x35;	/* This table needs 20x16 bytes */
		r23 = sc_ptr + 0x16;
		r24 = sc_ptr + 0x17;
		r25 = sc_ptr + 0x18;
		r26 = sc_ptr + 0x19;
		r27 = sc_ptr + 0x1a;
		r28 = sc_ptr + 0x1b;
		r29 = sc_ptr + 0x1c;
		r30 = sc_ptr + 0x1d;
		r31 = sc_ptr + 0x1e;
		r32 = sc_ptr + 0x1f;

		/* These remain fixed: */
		isrt2->re = isrt2->im = ISRT2;
		cc0  ->re = cc0  ->im = c	;
		ss0  ->re = ss0  ->im = s	;
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		sse2_rnd->re = sse2_rnd->im = crnd;

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
		tmp->re = base   [0];	tmp->im = base   [0];	++tmp;
		tmp->re = baseinv[0];	tmp->im = baseinv[0];	++tmp;
		/* [+2] slot is for [scale,scale] */
	}
	else
	{
		/* These static pointers are for [experimental] pure-ASM version of Mersenne-mod carry step: */
		si_ptr = &si[0];
		wt0_ptr = &wt0[0];	wt1_ptr = &wt1[0];
		/* Forward-weight multipliers: */
		tmp->re = 1.0;	tmp->im = 1.0;	++tmp;
		tmp->re = .50;	tmp->im = 1.0;	++tmp;
		tmp->re = 1.0;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .50;	++tmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		tmp->re = .50;	tmp->im = .50;	++tmp;
		tmp->re = .25;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .25;	++tmp;
		tmp->re = .25;	tmp->im = .25;	++tmp;
		/* Forward-base[] multipliers: */
		tmp->re = base   [0];	tmp->im = base   [0];	++tmp;
		tmp->re = base   [1];	tmp->im = base   [0];	++tmp;
		tmp->re = base   [0];	tmp->im = base   [1];	++tmp;
		tmp->re = base   [1];	tmp->im = base   [1];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->re = baseinv[0];	tmp->im = baseinv[0];	++tmp;
		tmp->re = baseinv[1];	tmp->im = baseinv[0];	++tmp;
		tmp->re = baseinv[0];	tmp->im = baseinv[1];	++tmp;
		tmp->re = baseinv[1];	tmp->im = baseinv[1];	++tmp;
	}

		/* Floating-point sign mask used for FABS on packed doubles: */
		sign_mask = sm_ptr;
		*sign_mask++ = (uint64)0x7FFFFFFFFFFFFFFFull;
		*sign_mask-- = (uint64)0x7FFFFFFFFFFFFFFFull;

		/* Test it: */
		tmp  = half_arr + 16;	/* ptr to local SSE2-floating-point storage */
		tmp->re = -ISRT2;	tmp->im = -ISRT2;

	  #if 0
		__asm	mov	eax, tmp
		__asm	mov	ebx, sign_mask
		__asm	movaps	xmm0,[eax]
		__asm	andpd	xmm0,[ebx]
		__asm	movaps	[eax],xmm0

		ASSERT(HERE, tmp->re == ISRT2, "sign_mask0");
		ASSERT(HERE, tmp->im == ISRT2, "sign_mask1");

		// Set up the quadrupled-32-bit-int SSE constants used by the carry macros:
		sse_bw  = sm_ptr + 2;
		__asm	mov	eax, bw
		__asm	mov	ebx, sse_bw
		__asm	movd	xmm0,eax	/* Move actual *value* of reg eax into low 32 bits of xmm0 */
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0

		sse_sw  = sm_ptr + 4;
		__asm	lea	eax, sw
		__asm	mov	ebx, sse_sw
		__asm	movd	xmm0,[eax]	/* Variant 2: Move contents of address pointed to by reg eax into low 32 bits of xmm0 */
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0

		sse_nm1 = sm_ptr + 6;
		__asm	lea	eax, nm1
		__asm	mov	ebx, sse_nm1
		__asm	movd	xmm0,[eax]
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0
	  #else
		sse_bw  = sm_ptr + 2;
		tmp64 = (uint64)bw;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_bw++ = tmp64;
		*sse_bw-- = tmp64;

		sse_sw  = sm_ptr + 4;
		tmp64 = (uint64)sw;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_sw++ = tmp64;
		*sse_sw-- = tmp64;

		sse_nm1 = sm_ptr + 6;
		tmp64 = (uint64)nm1;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_nm1++ = tmp64;
		*sse_nm1-- = tmp64;
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
		tdat[ithread].r1 = __r0 + ithread*cslots_in_local_store;
	}
#endif

		bjmodn0 = (uint32*)(sm_ptr + 8);
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

	#ifdef USE_PTHREAD
		r1 = __r0 + cslots_in_local_store;
		/* Init thread 1-CY_THREADS's local stores and pointers: */
		for(i = 1; i < CY_THREADS; ++i) {
			/* Only care about the constants for each thread here, but easier to just copy the entire thread0 local store: */
			memcpy(r1, __r0, cslots_in_local_store<<4);	// bytewise copy treats complex and uint64 subdata the same
			r1 += cslots_in_local_store;
		}
	#endif

	#endif	// USE_SSE2

		/*   constant index offsets for load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
		p1 = NDIVR + ( (NDIVR >> DAT_BITS) << PAD_BITS );
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

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix16_ditN_cy_dif1.");

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
			i.e. the one that n2/16-separated FFT outputs need:
			*/
			_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix16_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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

    /* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
    if(!full_pass)
    {
		khi = 1;
	#if FFT_DEBUG
		fprintf(stderr, "radix16_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars addr & addp for this to compile properly: */
#ifdef USE_OMP
	omp_set_num_threads(CY_THREADS);
/*#undef PFETCH	*/
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,k1,k2,l,col,co2,co3,m,m2,\
		n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,\
		rt,it,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,\
		a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr,\
		a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi,\
		bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF,\
		cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF,\
		cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF\
	) default(shared) schedule(static)
#endif

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
		ASSERT(HERE, tdat[ithread].r1 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].r1;
		ASSERT(HERE, ((tmp + 0x20)->re == ISRT2 && (tmp + 0x20)->im == ISRT2), "thread-local memcheck failed!");
		ASSERT(HERE, ((tmp + 0x34)->re == crnd && (tmp + 0x34)->im == crnd), "thread-local memcheck failed!");
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		ASSERT(HERE, (tmp + 0x35+10)->re * (tmp + 0x35+14)->re == 1.0 && (tmp + 0x35+10)->im * (tmp + 0x35+14)->im == 1.0, "thread-local memcheck failed!");
	} else {
		ASSERT(HERE, (tmp + 0x35)->re * (tmp + 0x35+1)->re == 1.0 && (tmp + 0x35)->im * (tmp + 0x35+1)->im == 1.0, "thread-local memcheck failed!");
	}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
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
		max_err->re = 0.0;	max_err->im = 0.0;
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
		#ifdef USE_SSE2
			cy_r01->re = _cy_r0[ithread];
			cy_r01->im = _cy_r1[ithread];
			cy_r23->re = _cy_r2[ithread];
			cy_r23->im = _cy_r3[ithread];
			cy_r45->re = _cy_r4[ithread];
			cy_r45->im = _cy_r5[ithread];
			cy_r67->re = _cy_r6[ithread];
			cy_r67->im = _cy_r7[ithread];
			cy_r89->re = _cy_r8[ithread];
			cy_r89->im = _cy_r9[ithread];
			cy_rAB->re = _cy_rA[ithread];
			cy_rAB->im = _cy_rB[ithread];
			cy_rCD->re = _cy_rC[ithread];
			cy_rCD->im = _cy_rD[ithread];
			cy_rEF->re = _cy_rE[ithread];
			cy_rEF->im = _cy_rF[ithread];
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
		#ifdef USE_SSE2
			cy_r01->re = _cy_r0[ithread];	cy_r01->im = _cy_i0[ithread];
			cy_r23->re = _cy_r1[ithread];	cy_r23->im = _cy_i1[ithread];
			cy_r45->re = _cy_r2[ithread];	cy_r45->im = _cy_i2[ithread];
			cy_r67->re = _cy_r3[ithread];	cy_r67->im = _cy_i3[ithread];
			cy_r89->re = _cy_r4[ithread];	cy_r89->im = _cy_i4[ithread];
			cy_rAB->re = _cy_r5[ithread];	cy_rAB->im = _cy_i5[ithread];
			cy_rCD->re = _cy_r6[ithread];	cy_rCD->im = _cy_i6[ithread];
			cy_rEF->re = _cy_r7[ithread];	cy_rEF->im = _cy_i7[ithread];
			cy_i01->re = _cy_r8[ithread];	cy_i01->im = _cy_i8[ithread];
			cy_i23->re = _cy_r9[ithread];	cy_i23->im = _cy_i9[ithread];
			cy_i45->re = _cy_rA[ithread];	cy_i45->im = _cy_iA[ithread];
			cy_i67->re = _cy_rB[ithread];	cy_i67->im = _cy_iB[ithread];
			cy_i89->re = _cy_rC[ithread];	cy_i89->im = _cy_iC[ithread];
			cy_iAB->re = _cy_rD[ithread];	cy_iAB->im = _cy_iD[ithread];
			cy_iCD->re = _cy_rE[ithread];	cy_iCD->im = _cy_iE[ithread];
			cy_iEF->re = _cy_rF[ithread];	cy_iEF->im = _cy_iF[ithread];
		#else
			cy_r0 = _cy_r0[ithread];	cy_i0 = _cy_i0[ithread];
			cy_r1 = _cy_r1[ithread];	cy_i1 = _cy_i1[ithread];
			cy_r2 = _cy_r2[ithread];	cy_i2 = _cy_i2[ithread];
			cy_r3 = _cy_r3[ithread];	cy_i3 = _cy_i3[ithread];
			cy_r4 = _cy_r4[ithread];	cy_i4 = _cy_i4[ithread];
			cy_r5 = _cy_r5[ithread];	cy_i5 = _cy_i5[ithread];
			cy_r6 = _cy_r6[ithread];	cy_i6 = _cy_i6[ithread];
			cy_r7 = _cy_r7[ithread];	cy_i7 = _cy_i7[ithread];
			cy_r8 = _cy_r8[ithread];	cy_i8 = _cy_i8[ithread];
			cy_r9 = _cy_r9[ithread];	cy_i9 = _cy_i9[ithread];
			cy_rA = _cy_rA[ithread];	cy_iA = _cy_iA[ithread];
			cy_rB = _cy_rB[ithread];	cy_iB = _cy_iB[ithread];
			cy_rC = _cy_rC[ithread];	cy_iC = _cy_iC[ithread];
			cy_rD = _cy_rD[ithread];	cy_iD = _cy_iD[ithread];
			cy_rE = _cy_rE[ithread];	cy_iE = _cy_iE[ithread];
			cy_rF = _cy_rF[ithread];	cy_iF = _cy_iF[ithread];
		#endif
		}
#if FFT_DEBUG
	fprintf(stderr, "\nIter = %d, full_pass = %d:\n",iter,full_pass);
#endif
		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
		#ifdef USE_SSE2
			for(j = jstart; j < jhi; j += 4)
			{
			/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
			*/
				j1 = (j & mask01) + br4[j&3];
		#elif defined(USE_SSE2)	/* This allows us to use #if 0 above and disable sse2-based *computation*, while still using sse2-style data layout */
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

		#ifdef DEBUG_SSE2
			rng_isaac_init(TRUE);
			jt = j1;		jp = j2;
			rt = 1024.0*1024.0*1024.0*1024.0;
			a[jt   ] = rt*rng_isaac_rand_double_norm_pm1();	a[jp   ] = rt*rng_isaac_rand_double_norm_pm1();	a[jt   +1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp   +1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p1] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p1+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p1+1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p2] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p2] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p2+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p2+1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p3] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p3] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p3+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p3+1] = rt*rng_isaac_rand_double_norm_pm1();	jt = j1 + p4;	jp = j2 + p4;
			a[jt   ] = rt*rng_isaac_rand_double_norm_pm1();	a[jp   ] = rt*rng_isaac_rand_double_norm_pm1();	a[jt   +1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp   +1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p1] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p1+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p1+1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p2] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p2] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p2+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p2+1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p3] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p3] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p3+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p3+1] = rt*rng_isaac_rand_double_norm_pm1();	jt = j1 + p8;	jp = j2 + p8;
			a[jt   ] = rt*rng_isaac_rand_double_norm_pm1();	a[jp   ] = rt*rng_isaac_rand_double_norm_pm1();	a[jt   +1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp   +1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p1] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p1+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p1+1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p2] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p2] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p2+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p2+1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p3] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p3] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p3+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p3+1] = rt*rng_isaac_rand_double_norm_pm1();	jt = j1 + p12;	jp = j2 + p12;
			a[jt   ] = rt*rng_isaac_rand_double_norm_pm1();	a[jp   ] = rt*rng_isaac_rand_double_norm_pm1();	a[jt   +1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp   +1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p1] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p1+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p1+1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p2] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p2] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p2+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p2+1] = rt*rng_isaac_rand_double_norm_pm1();
			a[jt+p3] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p3] = rt*rng_isaac_rand_double_norm_pm1();	a[jt+p3+1] = rt*rng_isaac_rand_double_norm_pm1();	a[jp+p3+1] = rt*rng_isaac_rand_double_norm_pm1();
		#endif
		#if 0//def DEBUG_SSE2
			if(j < 4 && full_pass && iter <= 1)
			{
			jt = j1;		jp = j2;
			fprintf(stderr, "radix16_carry: A_in[00] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt   ],a[jt   +1],a[jp   ],a[jp   +1]);
			fprintf(stderr, "radix16_carry: A_in[01] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p1],a[jt+p1+1],a[jp+p1],a[jp+p1+1]);
			fprintf(stderr, "radix16_carry: A_in[02] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p2],a[jt+p2+1],a[jp+p2],a[jp+p2+1]);
			fprintf(stderr, "radix16_carry: A_in[03] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p3],a[jt+p3+1],a[jp+p3],a[jp+p3+1]);	jt = j1 + p4;	jp = j2 + p4;
			fprintf(stderr, "radix16_carry: A_in[04] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt   ],a[jt   +1],a[jp   ],a[jp   +1]);
			fprintf(stderr, "radix16_carry: A_in[05] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p1],a[jt+p1+1],a[jp+p1],a[jp+p1+1]);
			fprintf(stderr, "radix16_carry: A_in[06] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p2],a[jt+p2+1],a[jp+p2],a[jp+p2+1]);
			fprintf(stderr, "radix16_carry: A_in[07] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p3],a[jt+p3+1],a[jp+p3],a[jp+p3+1]);	jt = j1 + p8;	jp = j2 + p8;
			fprintf(stderr, "radix16_carry: A_in[08] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt   ],a[jt   +1],a[jp   ],a[jp   +1]);
			fprintf(stderr, "radix16_carry: A_in[09] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p1],a[jt+p1+1],a[jp+p1],a[jp+p1+1]);
			fprintf(stderr, "radix16_carry: A_in[10] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p2],a[jt+p2+1],a[jp+p2],a[jp+p2+1]);
			fprintf(stderr, "radix16_carry: A_in[11] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p3],a[jt+p3+1],a[jp+p3],a[jp+p3+1]);	jt = j1 + p12;	jp = j2 + p12;
			fprintf(stderr, "radix16_carry: A_in[12] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt   ],a[jt   +1],a[jp   ],a[jp   +1]);
			fprintf(stderr, "radix16_carry: A_in[13] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p1],a[jt+p1+1],a[jp+p1],a[jp+p1+1]);
			fprintf(stderr, "radix16_carry: A_in[14] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p2],a[jt+p2+1],a[jp+p2],a[jp+p2+1]);
			fprintf(stderr, "radix16_carry: A_in[15] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p3],a[jt+p3+1],a[jp+p3],a[jp+p3+1]);
			fprintf(stderr, "\n");
			}
		#endif

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
				SSE2_RADIX4_DIT_0TWIDDLE_B(r1)
			#else
				add0 = &a[j1];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;
				/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r1)
			#endif

			/*...Block 2: */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
				SSE2_RADIX4_DIT_0TWIDDLE_B(r9)
			#else
				add0 = &a[j1+p4];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;
				/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r9)
			#endif

			/*...Block 3: */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
				SSE2_RADIX4_DIT_0TWIDDLE_B(r17)
			#else
				add0 = &a[j1+p8];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;
				/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r17)
			#endif

			/*...Block 4: */
			#if 1
				__asm	add eax, edi
				__asm	add ebx, edi
				__asm	add ecx, edi
				__asm	add edx, edi
				SSE2_RADIX4_DIT_0TWIDDLE_B(r25)
			#else
				add0 = &a[j1+p12];
				add1 = add0+p1;
				add2 = add0+p2;
				add3 = add0+p3;
				/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r25)
			#endif

			/****************************************************************************************************
			!...and now do four more radix-4 transforms, including the internal [no external]twiddle factors:   !
			****************************************************************************************************/

			/*...Block 1: Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
				__asm	mov	eax, r1

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
				__asm	mov	eax, r5
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
				__asm	mov	eax, r3
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
				__asm	mov	eax, r7
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

				SSE2_RADIX16_DIT_NOTWIDDLE(add0,p1,p2,p3,p4,r1,r3,r5,r7,r9,r11,r17,r25,isrt2,cc0);

		  #endif

		#else	/* !USE_SSE2 */

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

		#endif	/* USE_SSE2 */

				/*...and combine those to complete the radix-16 transform and do the carries. Since the outputs would
				normally be getting dispatched to 16 separate blocks of the A-array, we need 16 separate carries.	*/

			#ifndef USE_SSE2

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
				}
				else	/* MODULUS_TYPE_FERMAT */
				{
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
				}

			#elif defined(USE_SSE2)

				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
				/* In SSE2 mode, the data are arranged in memory like so, where we view things in 16-byte chunks:

					R0/r1 :	a0.re,b0.re		I0/r2 :	a0.im,b0.im
					R1/r3 :	a1.re,b1.re		I1/r4 :	a1.im,b1.im
					R2/r5 :	a2.re,b2.re		I2/r6 :	a2.im,b2.im
					R3/r7 :	a3.re,b3.re		I3/r8 :	a3.im,b3.im
					R4/r9 :	a4.re,b4.re		I4/r10:	a4.im,b4.im
					R5/r11:	a5.re,b5.re		I5/r12:	a5.im,b5.im
					R6/r13:	a6.re,b6.re		I6/r14:	a6.im,b6.im
					R7/r15:	a7.re,b7.re		I7/r16:	a7.im,b7.im
					R8/r17:	a8.re,b8.re		I8/r18:	a8.im,b8.im
					R9/r19:	a9.re,b9.re		I9/r20:	a9.im,b9.im
					Ra/r21:	aA.re,bA.re		Ia/r22:	aA.im,bA.im
					Rb/r23:	aB.re,bB.re		Ib/r24:	aB.im,bB.im
					Rc/r25:	aC.re,bC.re		Ic/r26:	aC.im,bC.im
					Rd/r27:	aD.re,bD.re		Id/r28:	aD.im,bD.im
					Re/r29:	aE.re,bE.re		Ie/r30:	aE.im,bE.im
					Rf/r31:	aF.re,bF.re		If/r32:	aF.im,bF.im

				Where the R's and I's map to the local temps as follows: R0:f ==> r1:31:2, I0:f ==> r2:32:2 , and the
				a's and b's of each pair represent doubles which in non-SSE2 mode would be getting processed in the same relative
				position on subsequent passes through the main-array for-loop, i.e. among which the carry propagation proceeds as

					a0.re -> a0.im -> b0.re -> b0.im .

				Because of the indesirable intra-xmm-register data dependency this leads to, we instead require data arranged as

					R0/r1 :	a0.re,a1.re		I0/r2 :	a0.im,a1.im
					R1/r3 :	b0.re,b1.re		I1/r4 :	b0.im,b1.im

				We need to interleave these pairwise so as to swap the high word of each even-indexed R-and-I-pair
				with the low word of the subsequent odd-indexed pair, e.g. for R0/r1 and R1/r3:

						low		high	low		high
					R0	[a0.re,b0.re]	[a1.re,b1.re]	R1
						   |      \       /      |
						   |        \   /        |
						   |          x          |
						   |        /   \        |
						   V      /       \      V
					R0~	[a0.re,a1.re]	[b0.re,b1.re]	R1~, and analogously for I0/r2 and I1/r4.

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
				#if 0

					l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
					n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
					n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
					sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
					sinwtm1 = si[nwt-l-1];

					wtl     =wt0[    l  ];
					wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
					wtlp1   =wt0[    l+1];
					wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

					tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
					tmp->re = wtl;		tmp->im = wtl;	++tmp;
					tmp->re = wtn;		tmp->im = wtn;	++tmp;
					tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
					tmp->re = wtnm1;	tmp->im = wtnm1;

					add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
					add2 = &wt1[co2-1];
					add3 = &wt1[co3-1];

				#elif defined(COMPILER_TYPE_MSVC)

					__asm	mov		eax, j
					__asm	mov		ebx, nwt
					__asm	mov		ecx, n
					__asm	sub		ebx, 1
					__asm	and		eax, ebx	// eax = l
					__asm	add		ebx, 1
					__asm	sub		ebx, eax	// ebx = nwt-l
					//asm	mov		  l, eax
					__asm	shl		eax, 2		// 4 bytes for array-of-ints
					__asm	shl		ebx, 2		// 4 bytes for array-of-ints
					__asm	mov		esi, si_ptr	// Master copy of si_ptr
					__asm	mov		edi, esi
					__asm	add		edi, eax	// &si[l]
					__asm	mov		edx,[edi]	//  si[l]
					__asm	sub		ecx, edx	// ecx = n-si[l]
					__asm	mov		n_minus_sil  , ecx
					__asm	add		ecx, edx	// ecx = n
					__asm	add		edi, 0x4	// &si[l+1]
					__asm	mov		edx,[edi]	//  si[l+1]
					__asm	sub		ecx, edx	// ecx = n-si[l]
					__asm	mov		n_minus_silp1, ecx
					__asm	mov		edi, esi
					__asm	add		edi, ebx	// &si[nwt-l]
					__asm	mov		edx,[edi]	//  si[nwt-l]
					__asm	mov		sinwt  , edx
					__asm	sub		edi, 0x4	// &si[nwt-l-1]
					__asm	mov		edx,[edi]	//  si[nwt-l-1]
					__asm	mov		sinwtm1, edx

					__asm	shl		eax, 1			// 8 bytes for array-of-doubles
					__asm	shl		ebx, 1			// 8 bytes for array-of-doubles
					__asm	mov		esi, wt0_ptr	// Master copy of wt0_ptr
					__asm	mov		edi, esi
					__asm	add		edi, eax		// &wt0[l]
					__asm	movlpd	xmm0,[edi    ]	//  wtl		NOTE: movhpd/movlpd preferable to movupd/shufpd
					__asm	movhpd	xmm0,[edi    ]
					__asm	movlpd	xmm1,[edi+0x8]	//  wtlp1
					__asm	movhpd	xmm1,[edi+0x8]
					__asm	mov		edi, esi
					__asm	add		edi, ebx		// &wt0[nwt-l]
					__asm	movlpd	xmm2,[edi    ]	//  wtn		NOTE: movhpd/movlpd preferable to movupd/shufpd
					__asm	movhpd	xmm2,[edi    ]
					__asm	movlpd	xmm3,[edi-0x8]	//  wtnm1
					__asm	movhpd	xmm3,[edi-0x8]
					__asm	mov		ecx, scale_ptr
					__asm	movlpd	xmm4,[ecx]
					__asm	movhpd	xmm4,[ecx]
					__asm	mulpd	xmm2,xmm4
					__asm	mulpd	xmm3,xmm4
					__asm	mov		edx, half_arr
					__asm	add		edx, 0x100		// ptr to local SSE2-floating-point storage
					__asm	movaps	[edx     ],xmm0	// wtl
					__asm	movaps	[edx+0x10],xmm2	// wtn
					__asm	movaps	[edx+0x20],xmm1	// wtlp1`
					__asm	movaps	[edx+0x30],xmm3	// wtnm1

					__asm	mov		eax, col
					__asm	mov		ebx, co2
					__asm	mov		ecx, co3
					__asm	mov		esi, wt1_ptr	// Master copy of wt0_ptr
					__asm	mov		edi, esi
					__asm	mov		edx, esi
					__asm	shl		eax, 3			// 8 bytes for array-of-doubles
					__asm	shl		ebx, 3			// 8 bytes for array-of-doubles
					__asm	shl		ecx, 3			// 8 bytes for array-of-doubles
					__asm	sub		ebx, 0x8
					__asm	sub		ecx, 0x8
					__asm	add		esi, eax		// add1 = &wt1[col  ];
					__asm	add		edi, ebx		// add2 = &wt1[co2-1];
					__asm	add		edx, ecx		// add3 = &wt1[co3-1];
					__asm	mov		add1,esi
					__asm	mov		add2,edi
					__asm	mov		add3,edx

				#else	/* GCC-style inline ASM: */

					#if OS_BITS == 32

					__asm__ volatile (\
						"movl	%[__j],%%eax			\n\t"\
						"movl	%[__nwt],%%ebx			\n\t"\
						"movl	%[__n],%%ecx			\n\t"\
						"subl	$1,%%ebx			\n\t"\
						"andl	%%ebx,%%eax				\n\t"\
						"addl	$1,%%ebx			\n\t"\
						"subl	%%eax,%%ebx				\n\t"\
						"shll	$2,%%eax				\n\t"\
						"shll	$2,%%ebx				\n\t"\
						"movl	%[__si_ptr],%%esi		\n\t"\
						"movl	%%esi,%%edi				\n\t"\
						"addl	%%eax,%%edi				\n\t"\
						"movl	(%%edi),%%edx			\n\t"\
						"subl	%%edx,%%ecx				\n\t"\
						"movl	%%ecx,%[__n_minus_sil]	\n\t"\
						"addl	%%edx,%%ecx				\n\t"\
						"addl	$0x4,%%edi				\n\t"\
						"movl	(%%edi),%%edx			\n\t"\
						"subl	%%edx,%%ecx				\n\t"\
						"movl	%%ecx,%[__n_minus_silp1]\n\t"\
						"movl	%%esi,%%edi				\n\t"\
						"addl	%%ebx,%%edi				\n\t"\
						"movl	(%%edi),%%edx			\n\t"\
						"movl	%%edx,%[__sinwt]		\n\t"\
						"subl	$0x4,%%edi				\n\t"\
						"movl	(%%edi),%%edx			\n\t"\
						"movl	%%edx,%[__sinwtm1]		\n\t"\
						"shll	$1,%%eax				\n\t"\
						"shll	$1,%%ebx				\n\t"\
						"movl	%[__wt0_ptr],%%esi		\n\t"\
						"movl	%%esi,%%edi				\n\t"\
						"addl	%%eax,%%edi				\n\t"\
						"movlpd	   (%%edi),%%xmm0		\n\t"\
						"movhpd	   (%%edi),%%xmm0		\n\t"\
						"movlpd	0x8(%%edi),%%xmm1		\n\t"\
						"movhpd	0x8(%%edi),%%xmm1		\n\t"\
						"movl	%%esi,%%edi				\n\t"\
						"addl	%%ebx,%%edi				\n\t"\
						"movlpd	    (%%edi),%%xmm2		\n\t"\
						"movhpd	    (%%edi),%%xmm2		\n\t"\
						"movlpd	-0x8(%%edi),%%xmm3		\n\t"\
						"movhpd	-0x8(%%edi),%%xmm3		\n\t"\
						"movl	%[__scale_ptr],%%ecx	\n\t"\
						"movlpd	(%%ecx),%%xmm4			\n\t"\
						"movhpd	(%%ecx),%%xmm4			\n\t"\
						"mulpd	%%xmm4,%%xmm2			\n\t"\
						"mulpd	%%xmm4,%%xmm3			\n\t"\
						"movl	%[__half_arr],%%edx		\n\t"\
						"addl	$0x100,%%edx			\n\t"\
						"movaps	%%xmm0,    (%%edx)		\n\t"\
						"movaps	%%xmm2,0x10(%%edx)		\n\t"\
						"movaps	%%xmm1,0x20(%%edx)		\n\t"\
						"movaps	%%xmm3,0x30(%%edx)		\n\t"\
						"movl	%[__col],%%eax			\n\t"\
						"movl	%[__co2],%%ebx			\n\t"\
						"movl	%[__co3],%%ecx			\n\t"\
						"movl	%[__wt1_ptr],%%esi		\n\t"\
						"movl	%%esi,%%edi				\n\t"\
						"movl	%%esi,%%edx				\n\t"\
						"shll	$3,%%eax				\n\t"\
						"shll	$3,%%ebx				\n\t"\
						"shll	$3,%%ecx				\n\t"\
						"subl	$0x8,%%ebx				\n\t"\
						"subl	$0x8,%%ecx				\n\t"\
						"addl	%%eax,%%esi				\n\t"\
						"addl	%%ebx,%%edi				\n\t"\
						"addl	%%ecx,%%edx				\n\t"\
						"movl	%%esi,%[__add1]			\n\t"\
						"movl	%%edi,%[__add2]			\n\t"\
						"movl	%%edx,%[__add3]			\n\t"\
					:						/* outputs: none */\
					:	[__j]				"m" (j)			/* All inputs from memory addresses here */\
					,	[__nwt]				"m" (nwt)\
					,	[__n]				"m" (n)\
					,	[__si_ptr]			"m" (si_ptr)\
					,	[__n_minus_sil]		"m" (n_minus_sil)\
					,	[__n_minus_silp1]	"m" (n_minus_silp1)\
					,	[__sinwt]			"m" (sinwt)\
					,	[__sinwtm1]			"m" (sinwtm1)\
					,	[__wt0_ptr]			"m" (wt0_ptr)\
					,	[__wt1_ptr]			"m" (wt1_ptr)\
					,	[__scale_ptr]		"m" (scale_ptr)\
					,	[__half_arr]		"m" (half_arr)\
					,	[__col]				"m" (col)\
					,	[__co2]				"m" (co2)\
					,	[__co3]				"m" (co3)\
					,	[__add1]			"m" (add1)\
					,	[__add2]			"m" (add2)\
					,	[__add3]			"m" (add3)\
						: "eax","ebx","ecx","edx","esi","edi"	/* Clobbered registers, excluding xmm* in 32-bit mode */\
					);

					#else

					__asm__ volatile (\
						"movslq	%[__j],%%rax			\n\t"\
						"movslq	%[__nwt],%%rbx			\n\t"\
						"movslq	%[__n],%%rcx			\n\t"\
						"subq	$1,%%rbx				\n\t"\
						"andq	%%rbx,%%rax				\n\t"\
						"addq	$1,%%rbx				\n\t"\
						"subq	%%rax,%%rbx				\n\t"\
						"shlq	$2,%%rax				\n\t"\
						"shlq	$2,%%rbx				\n\t"\
						"movq	%[__si_ptr],%%rsi		\n\t"\
						"movq	%%rsi,%%rdi				\n\t"\
						"addq	%%rax,%%rdi				\n\t"\
						"movslq	(%%rdi),%%rdx			\n\t"\
						"subq	%%rdx,%%rcx				\n\t"\
						"movl	%%ecx,%[__n_minus_sil]	\n\t"\
						"addq	%%rdx,%%rcx				\n\t"\
						"addq	$0x4,%%rdi				\n\t"\
						"movslq	(%%rdi),%%rdx			\n\t"\
						"subq	%%rdx,%%rcx				\n\t"\
						"movl	%%ecx,%[__n_minus_silp1]\n\t"\
						"movq	%%rsi,%%rdi				\n\t"\
						"addq	%%rbx,%%rdi				\n\t"\
						"movslq	(%%rdi),%%rdx			\n\t"\
						"movl	%%edx,%[__sinwt]		\n\t"\
						"subq	$0x4,%%rdi				\n\t"\
						"movq	(%%rdi),%%rdx			\n\t"\
						"movl	%%edx,%[__sinwtm1]		\n\t"\
						"shlq	$1,%%rax				\n\t"\
						"shlq	$1,%%rbx				\n\t"\
						"movq	%[__wt0_ptr],%%rsi		\n\t"\
						"movq	%%rsi,%%rdi				\n\t"\
						"addq	%%rax,%%rdi				\n\t"\
						"movlpd	   (%%rdi),%%xmm0		\n\t"\
						"movhpd	   (%%rdi),%%xmm0		\n\t"\
						"movlpd	0x8(%%rdi),%%xmm1		\n\t"\
						"movhpd	0x8(%%rdi),%%xmm1		\n\t"\
						"movq	%%rsi,%%rdi				\n\t"\
						"addq	%%rbx,%%rdi				\n\t"\
						"movlpd	    (%%rdi),%%xmm2		\n\t"\
						"movhpd	    (%%rdi),%%xmm2		\n\t"\
						"movlpd	-0x8(%%rdi),%%xmm3		\n\t"\
						"movhpd	-0x8(%%rdi),%%xmm3		\n\t"\
						"movq	%[__scale_ptr],%%rcx	\n\t"\
						"movlpd	(%%rcx),%%xmm4			\n\t"\
						"movhpd	(%%rcx),%%xmm4			\n\t"\
						"mulpd	%%xmm4,%%xmm2			\n\t"\
						"mulpd	%%xmm4,%%xmm3			\n\t"\
						"movq	%[__half_arr],%%rdx		\n\t"\
						"addq	$0x100,%%rdx			\n\t"\
						"movaps	%%xmm0,    (%%rdx)		\n\t"\
						"movaps	%%xmm2,0x10(%%rdx)		\n\t"\
						"movaps	%%xmm1,0x20(%%rdx)		\n\t"\
						"movaps	%%xmm3,0x30(%%rdx)		\n\t"\
						"movslq	%[__col],%%rax			\n\t"\
						"movslq	%[__co2],%%rbx			\n\t"\
						"movslq	%[__co3],%%rcx			\n\t"\
						"movq	%[__wt1_ptr],%%rsi		\n\t"\
						"movq	%%rsi,%%rdi				\n\t"\
						"movq	%%rsi,%%rdx				\n\t"\
						"shlq	$3,%%rax				\n\t"\
						"shlq	$3,%%rbx				\n\t"\
						"shlq	$3,%%rcx				\n\t"\
						"subq	$0x8,%%rbx				\n\t"\
						"subq	$0x8,%%rcx				\n\t"\
						"addq	%%rax,%%rsi				\n\t"\
						"addq	%%rbx,%%rdi				\n\t"\
						"addq	%%rcx,%%rdx				\n\t"\
						"movq	%%rsi,%[__add1]			\n\t"\
						"movq	%%rdi,%[__add2]			\n\t"\
						"movq	%%rdx,%[__add3]			\n\t"\
					:						/* outputs: none */\
					:	[__j]				"m" (j)			/* All inputs from memory addresses here */\
					,	[__nwt]				"m" (nwt)\
					,	[__n]				"m" (n)\
					,	[__si_ptr]			"m" (si_ptr)\
					,	[__n_minus_sil]		"m" (n_minus_sil)\
					,	[__n_minus_silp1]	"m" (n_minus_silp1)\
					,	[__sinwt]			"m" (sinwt)\
					,	[__sinwtm1]			"m" (sinwtm1)\
					,	[__wt0_ptr]			"m" (wt0_ptr)\
					,	[__wt1_ptr]			"m" (wt1_ptr)\
					,	[__scale_ptr]		"m" (scale_ptr)\
					,	[__half_arr]		"m" (half_arr)\
					,	[__col]				"m" (col)\
					,	[__co2]				"m" (co2)\
					,	[__co3]				"m" (co3)\
					,	[__add1]			"m" (add1)\
					,	[__add2]			"m" (add2)\
					,	[__add3]			"m" (add3)\
						: "rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4"	/* Clobbered registers */\
					);

					#endif	/* 32/64-bit */

				#endif	/* #if 0 */
/*
				SSE2_cmplx_carry_norm_pow2_errcheck0(r1 ,add1,add2,add3,cy_r01,bjmodn0,bjmodn1);//add1 += 2;	add2 -= 2;	add3 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck1(r5 ,add1,add2,add3,cy_r23,bjmodn2,bjmodn3);//add1 += 2;	add2 -= 2;	add3 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck1(r9 ,add1,add2,add3,cy_r45,bjmodn4,bjmodn5);//add1 += 2;	add2 -= 2;	add3 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck1(r13,add1,add2,add3,cy_r67,bjmodn6,bjmodn7);//add1 += 2;	add2 -= 2;	add3 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck1(r17,add1,add2,add3,cy_r89,bjmodn8,bjmodn9);//add1 += 2;	add2 -= 2;	add3 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck1(r21,add1,add2,add3,cy_rAB,bjmodnA,bjmodnB);//add1 += 2;	add2 -= 2;	add3 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck1(r25,add1,add2,add3,cy_rCD,bjmodnC,bjmodnD);//add1 += 2;	add2 -= 2;	add3 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck1(r29,add1,add2,add3,cy_rEF,bjmodnE,bjmodnF);
//
				SSE2_cmplx_carry_norm_pow2_errcheck0_2x(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0,bjmodn1,bjmodn2,bjmodn3);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2x(r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4,bjmodn5,bjmodn6,bjmodn7);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2x(r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8,bjmodn9,bjmodnA,bjmodnB);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2x(r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC,bjmodnD,bjmodnE,bjmodnF);
*/
			#if defined(COMPILER_TYPE_MSVC)

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC);
			  #else
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC);
			  #endif

			#else	/* GCC-style inline ASM: */

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			  #endif

				/* Bizarre - when I disabled the diagnostic prints above and below, the resulting GCC build immediately gave
					fatal roundoff errors starting on iteration #5 - so insert the bogus [never taken] if() here as a workaround.
					Equally bizarre, inserting the bogus if() *before* the 4 carry-macro calls above gave the correct result as well,
					but ran fully 10% slower. Good old GCC...
				*/
				if(j < 0)
				{
					fprintf(stderr, "Iter %3d\n",iter);
				}
			#endif

				#if 0

					l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
					n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
					n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
					sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
					sinwtm1 = si[nwt-l-1];

					wtl     =wt0[    l  ];
					wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
					wtlp1   =wt0[    l+1];
					wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

					tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
					tmp->re = wtl;		tmp->im = wtl;	++tmp;
					tmp->re = wtn;		tmp->im = wtn;	++tmp;
					tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
					tmp->re = wtnm1;	tmp->im = wtnm1;

				/*	i =((uint32)(sw - bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

					co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

					add1 = &wt1[col  ];
					add2 = &wt1[co2-1];

				#elif defined(COMPILER_TYPE_MSVC)

					__asm	mov		eax, j
					__asm	mov		ebx, nwt
					__asm	mov		ecx, n
					__asm	add		eax, 2
					__asm	sub		ebx, 1
					__asm	and		eax, ebx	// eax = l
					__asm	add		ebx, 1
					__asm	sub		ebx, eax	// ebx = nwt-l
					//asm	mov		  l, eax
					__asm	shl		eax, 2		// 4 bytes for array-of-ints
					__asm	shl		ebx, 2		// 4 bytes for array-of-ints
					__asm	mov		esi, si_ptr	// Master copy of si_ptr
					__asm	mov		edi, esi
					__asm	add		edi, eax	// &si[l]
					__asm	mov		edx,[edi]	//  si[l]
					__asm	sub		ecx, edx	// ecx = n-si[l]
					__asm	mov		n_minus_sil  , ecx
					__asm	add		ecx, edx	// ecx = n
					__asm	add		edi, 0x4	// &si[l+1]
					__asm	mov		edx,[edi]	//  si[l+1]
					__asm	sub		ecx, edx	// ecx = n-si[l]
					__asm	mov		n_minus_silp1, ecx
					__asm	mov		edi, esi
					__asm	add		edi, ebx	// &si[nwt-l]
					__asm	mov		edx,[edi]	//  si[nwt-l]
					__asm	mov		sinwt  , edx
					__asm	sub		edi, 0x4	// &si[nwt-l-1]
					__asm	mov		edx,[edi]	//  si[nwt-l-1]
					__asm	mov		sinwtm1, edx

					__asm	shl		eax, 1			// 8 bytes for array-of-doubles
					__asm	shl		ebx, 1			// 8 bytes for array-of-doubles
					__asm	mov		esi, wt0_ptr	// Master copy of wt0_ptr
					__asm	mov		edi, esi
					__asm	add		edi, eax		// &wt0[l]
					__asm	movlpd	xmm0,[edi    ]	//  wtl		NOTE: movhpd/movlpd preferable to movupd/shufpd
					__asm	movhpd	xmm0,[edi    ]
					__asm	movlpd	xmm1,[edi+0x8]	//  wtlp1
					__asm	movhpd	xmm1,[edi+0x8]
					__asm	mov		edi, esi
					__asm	add		edi, ebx		// &wt0[nwt-l]
					__asm	movlpd	xmm2,[edi    ]	//  wtn		NOTE: movhpd/movlpd preferable to movupd/shufpd
					__asm	movhpd	xmm2,[edi    ]
					__asm	movlpd	xmm3,[edi-0x8]	//  wtnm1
					__asm	movhpd	xmm3,[edi-0x8]
					__asm	mov		ecx, scale_ptr
					__asm	movlpd	xmm4,[ecx]
					__asm	movhpd	xmm4,[ecx]
					__asm	mulpd	xmm2,xmm4
					__asm	mulpd	xmm3,xmm4
					__asm	mov		edx, half_arr
					__asm	add		edx, 0x100		// ptr to local SSE2-floating-point storage
					__asm	movaps	[edx     ],xmm0	// wtl
					__asm	movaps	[edx+0x10],xmm2	// wtn
					__asm	movaps	[edx+0x20],xmm1	// wtlp1`
					__asm	movaps	[edx+0x30],xmm3	// wtnm1

					__asm	mov		eax, col
					__asm	mov		ecx, co3
					__asm	mov		esi, wt1_ptr	// Master copy of wt0_ptr
					__asm	mov		edi, esi
					__asm	mov		co2, ecx		// co2 = co3
					__asm	shl		eax, 3			// 8 bytes for array-of-doubles
					__asm	shl		ecx, 3			// 8 bytes for array-of-doubles
					__asm	sub		ecx, 0x8
					__asm	add		esi, eax		// add1 = &wt1[col  ];
					__asm	add		edi, ecx		// add2 = &wt1[co2-1];
					__asm	mov		add1,esi
					__asm	mov		add2,edi

				#else	/* GCC-style inline ASM: */

				  #if OS_BITS == 32

					__asm__ volatile (\
					  "movl	%[__j],%%eax			\n\t"\
					  "movl	%[__nwt],%%ebx			\n\t"\
					  "movl	%[__n],%%ecx			\n\t"\
					  "addl	$2,%%eax				\n\t"\
					  "subl	$1,%%ebx			\n\t"\
					  "andl	%%ebx,%%eax				\n\t"\
					  "addl	$1,%%ebx			\n\t"\
					  "subl	%%eax,%%ebx				\n\t"\
					  "shll	$2,%%eax				\n\t"\
					  "shll	$2,%%ebx				\n\t"\
					  "movl	%[__si_ptr],%%esi		\n\t"\
					  "movl	%%esi,%%edi				\n\t"\
					  "addl	%%eax,%%edi				\n\t"\
					  "movl	(%%edi),%%edx			\n\t"\
					  "subl	%%edx,%%ecx				\n\t"\
					  "movl	%%ecx,%[__n_minus_sil]	\n\t"\
					  "addl	%%edx,%%ecx				\n\t"\
					  "addl	$0x4,%%edi				\n\t"\
					  "movl	(%%edi),%%edx			\n\t"\
					  "subl	%%edx,%%ecx				\n\t"\
					  "movl	%%ecx,%[__n_minus_silp1]\n\t"\
					  "movl	%%esi,%%edi				\n\t"\
					  "addl	%%ebx,%%edi				\n\t"\
					  "movl	(%%edi),%%edx			\n\t"\
					  "movl	%%edx,%[__sinwt]		\n\t"\
					  "subl	$0x4,%%edi				\n\t"\
					  "movl	(%%edi),%%edx			\n\t"\
					  "movl	%%edx,%[__sinwtm1]		\n\t"\
					  "shll	$1,%%eax				\n\t"\
					  "shll	$1,%%ebx				\n\t"\
					  "movl	%[__wt0_ptr],%%esi		\n\t"\
					  "movl	%%esi,%%edi				\n\t"\
					  "addl	%%eax,%%edi				\n\t"\
					  "movlpd	   (%%edi),%%xmm0		\n\t"\
					  "movhpd	   (%%edi),%%xmm0		\n\t"\
					  "movlpd	0x8(%%edi),%%xmm1		\n\t"\
					  "movhpd	0x8(%%edi),%%xmm1		\n\t"\
					  "movl	%%esi,%%edi				\n\t"\
					  "addl	%%ebx,%%edi				\n\t"\
					  "movlpd	    (%%edi),%%xmm2		\n\t"\
					  "movhpd	    (%%edi),%%xmm2		\n\t"\
					  "movlpd	-0x8(%%edi),%%xmm3		\n\t"\
					  "movhpd	-0x8(%%edi),%%xmm3		\n\t"\
					  "movl	%[__scale_ptr],%%ecx	\n\t"\
					  "movlpd	(%%ecx),%%xmm4			\n\t"\
					  "movhpd	(%%ecx),%%xmm4			\n\t"\
					  "mulpd	%%xmm4,%%xmm2			\n\t"\
					  "mulpd	%%xmm4,%%xmm3			\n\t"\
					  "movl	%[__half_arr],%%edx		\n\t"\
					  "addl	$0x100,%%edx			\n\t"\
					  "movaps	%%xmm0,    (%%edx)		\n\t"\
					  "movaps	%%xmm2,0x10(%%edx)		\n\t"\
					  "movaps	%%xmm1,0x20(%%edx)		\n\t"\
					  "movaps	%%xmm3,0x30(%%edx)		\n\t"\
					  "movl	%[__col],%%eax			\n\t"\
					  "movl	%[__co3],%%ecx			\n\t"\
					  "movl	%[__wt1_ptr],%%esi		\n\t"\
					  "movl	%%esi,%%edi				\n\t"\
					  "movl	%%ecx,%[__co2]			\n\t"\
					  "shll	$3,%%eax				\n\t"\
					  "shll	$3,%%ecx				\n\t"\
					  "subl	$0x8,%%ecx				\n\t"\
					  "addl	%%eax,%%esi				\n\t"\
					  "addl	%%ecx,%%edi				\n\t"\
					  "movl	%%esi,%[__add1]			\n\t"\
					  "movl	%%edi,%[__add2]			\n\t"\
					  :						/* outputs: none */\
					  :	[__j]				"m" (j)			/* All inputs from memory addresses here */\
					  ,	[__nwt]				"m" (nwt)\
					  ,	[__n]				"m" (n)\
					  ,	[__si_ptr]			"m" (si_ptr)\
					  ,	[__n_minus_sil]		"m" (n_minus_sil)\
					  ,	[__n_minus_silp1]	"m" (n_minus_silp1)\
					  ,	[__sinwt]			"m" (sinwt)\
					  ,	[__sinwtm1]			"m" (sinwtm1)\
					  ,	[__wt0_ptr]			"m" (wt0_ptr)\
					  ,	[__wt1_ptr]			"m" (wt1_ptr)\
					  ,	[__scale_ptr]		"m" (scale_ptr)\
					  ,	[__half_arr]		"m" (half_arr)\
					  ,	[__col]				"m" (col)\
					  ,	[__co2]				"m" (co2)\
					  ,	[__co3]				"m" (co3)\
					  ,	[__add1]			"m" (add1)\
					  ,	[__add2]			"m" (add2)\
					  , [__add3]			"m" (add3)\
					  : "eax","ebx","ecx","edx","esi","edi"	/* Clobbered registers, excluding xmm* in 32-bit mode */\
					  );

				  #else

					__asm__ volatile (\
					  "movslq	%[__j],%%rax			\n\t"\
					  "movslq	%[__nwt],%%rbx			\n\t"\
					  "movslq	%[__n],%%rcx			\n\t"\
					  "addq	$2,%%rax				\n\t"\
					  "subq	$1,%%rbx				\n\t"\
					  "andq	%%rbx,%%rax				\n\t"\
					  "addq	$1,%%rbx				\n\t"\
					  "subq	%%rax,%%rbx				\n\t"\
					  "shlq	$2,%%rax				\n\t"\
					  "shlq	$2,%%rbx				\n\t"\
					  "movq	%[__si_ptr],%%rsi		\n\t"\
					  "movq	%%rsi,%%rdi				\n\t"\
					  "addq	%%rax,%%rdi				\n\t"\
					  "movslq	(%%rdi),%%rdx			\n\t"\
					  "subq	%%rdx,%%rcx				\n\t"\
					  "movl	%%ecx,%[__n_minus_sil]	\n\t"\
					  "addq	%%rdx,%%rcx				\n\t"\
					  "addq	$0x4,%%rdi				\n\t"\
					  "movslq	(%%rdi),%%rdx			\n\t"\
					  "subq	%%rdx,%%rcx				\n\t"\
					  "movl	%%ecx,%[__n_minus_silp1]\n\t"\
					  "movq	%%rsi,%%rdi				\n\t"\
					  "addq	%%rbx,%%rdi				\n\t"\
					  "movslq	(%%rdi),%%rdx			\n\t"\
					  "movl	%%edx,%[__sinwt]		\n\t"\
					  "subq	$0x4,%%rdi				\n\t"\
					  "movq	(%%rdi),%%rdx			\n\t"\
					  "movl	%%edx,%[__sinwtm1]		\n\t"\
					  "shlq	$1,%%rax				\n\t"\
					  "shlq	$1,%%rbx				\n\t"\
					  "movq	%[__wt0_ptr],%%rsi		\n\t"\
					  "movq	%%rsi,%%rdi				\n\t"\
					  "addq	%%rax,%%rdi				\n\t"\
					  "movlpd	   (%%rdi),%%xmm0		\n\t"\
					  "movhpd	   (%%rdi),%%xmm0		\n\t"\
					  "movlpd	0x8(%%rdi),%%xmm1		\n\t"\
					  "movhpd	0x8(%%rdi),%%xmm1		\n\t"\
					  "movq	%%rsi,%%rdi				\n\t"\
					  "addq	%%rbx,%%rdi				\n\t"\
					  "movlpd	    (%%rdi),%%xmm2		\n\t"\
					  "movhpd	    (%%rdi),%%xmm2		\n\t"\
					  "movlpd	-0x8(%%rdi),%%xmm3		\n\t"\
					  "movhpd	-0x8(%%rdi),%%xmm3		\n\t"\
					  "movq	%[__scale_ptr],%%rcx	\n\t"\
					  "movlpd	(%%rcx),%%xmm4			\n\t"\
					  "movhpd	(%%rcx),%%xmm4			\n\t"\
					  "mulpd	%%xmm4,%%xmm2			\n\t"\
					  "mulpd	%%xmm4,%%xmm3			\n\t"\
					  "movq	%[__half_arr],%%rdx		\n\t"\
					  "addq	$0x100,%%rdx			\n\t"\
					  "movaps	%%xmm0,    (%%rdx)		\n\t"\
					  "movaps	%%xmm2,0x10(%%rdx)		\n\t"\
					  "movaps	%%xmm1,0x20(%%rdx)		\n\t"\
					  "movaps	%%xmm3,0x30(%%rdx)		\n\t"\
					  "movslq	%[__col],%%rax			\n\t"\
					  "movslq	%[__co3],%%rcx			\n\t"\
					  "movq	%[__wt1_ptr],%%rsi		\n\t"\
					  "movq	%%rsi,%%rdi				\n\t"\
					  "movl	%%ecx,%[__co2]			\n\t"\
					  "shlq	$3,%%rax				\n\t"\
					  "shlq	$3,%%rcx				\n\t"\
					  "subq	$0x8,%%rcx				\n\t"\
					  "addq	%%rax,%%rsi				\n\t"\
					  "addq	%%rcx,%%rdi				\n\t"\
					  "movq	%%rsi,%[__add1]			\n\t"\
					  "movq	%%rdi,%[__add2]			\n\t"\
					  :						/* outputs: none */\
					  :	[__j]				"m" (j)			/* All inputs from memory addresses here */\
					  ,	[__nwt]				"m" (nwt)\
					  ,	[__n]				"m" (n)\
					  ,	[__si_ptr]			"m" (si_ptr)\
					  ,	[__n_minus_sil]		"m" (n_minus_sil)\
					  ,	[__n_minus_silp1]	"m" (n_minus_silp1)\
					  ,	[__sinwt]			"m" (sinwt)\
					  ,	[__sinwtm1]			"m" (sinwtm1)\
					  ,	[__wt0_ptr]			"m" (wt0_ptr)\
					  ,	[__wt1_ptr]			"m" (wt1_ptr)\
					  ,	[__scale_ptr]		"m" (scale_ptr)\
					  ,	[__half_arr]		"m" (half_arr)\
					  ,	[__col]				"m" (col)\
					  ,	[__co2]				"m" (co2)\
					  ,	[__co3]				"m" (co3)\
					  ,	[__add1]			"m" (add1)\
					  ,	[__add2]			"m" (add2)\
					  ,	[__add3]			"m" (add3)\
					  /*: "rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4"	/ Clobbered registers */\
					  );

				  #endif	/* if OS_BITS == 32 */

				#endif	/* if COMPILER_TYPE == ... */
			/*
				SSE2_cmplx_carry_norm_pow2_errcheck2(r1 ,add1,add2     ,cy_r01,bjmodn0,bjmodn1);//add1 += 2;	add2 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck2(r5 ,add1,add2     ,cy_r23,bjmodn2,bjmodn3);//add1 += 2;	add2 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck2(r9 ,add1,add2     ,cy_r45,bjmodn4,bjmodn5);//add1 += 2;	add2 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck2(r13,add1,add2     ,cy_r67,bjmodn6,bjmodn7);//add1 += 2;	add2 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck2(r17,add1,add2     ,cy_r89,bjmodn8,bjmodn9);//add1 += 2;	add2 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck2(r21,add1,add2     ,cy_rAB,bjmodnA,bjmodnB);//add1 += 2;	add2 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck2(r25,add1,add2     ,cy_rCD,bjmodnC,bjmodnD);//add1 += 2;	add2 -= 2;
				SSE2_cmplx_carry_norm_pow2_errcheck2(r29,add1,add2     ,cy_rEF,bjmodnE,bjmodnF);
//
				SSE2_cmplx_carry_norm_pow2_errcheck2_2x(r1 ,add1,add2,cy_r01,cy_r23,bjmodn0,bjmodn1,bjmodn2,bjmodn3);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2x(r9 ,add1,add2,cy_r45,cy_r67,bjmodn4,bjmodn5,bjmodn6,bjmodn7);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2x(r17,add1,add2,cy_r89,cy_rAB,bjmodn8,bjmodn9,bjmodnA,bjmodnB);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2x(r25,add1,add2,cy_rCD,cy_rEF,bjmodnC,bjmodnD,bjmodnE,bjmodnF);
*/
			#if defined(COMPILER_TYPE_MSVC)

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r1 ,add1,add2,cy_r01,cy_r23,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r9 ,add1,add2,cy_r45,cy_r67,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r17,add1,add2,cy_r89,cy_rAB,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r25,add1,add2,cy_rCD,cy_rEF,bjmodnC);
			  #else
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r1 ,add1,add2,cy_r01,cy_r23,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r9 ,add1,add2,cy_r45,cy_r67,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r17,add1,add2,cy_r89,cy_rAB,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r25,add1,add2,cy_rCD,cy_rEF,bjmodnC);
			  #endif

			#else	/* GCC-style inline ASM: */

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r1 ,add1,add2,cy_r01,cy_r23,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r9 ,add1,add2,cy_r45,cy_r67,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r17,add1,add2,cy_r89,cy_rAB,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r25,add1,add2,cy_rCD,cy_rEF,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r1 ,add1,add2,cy_r01,cy_r23,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r9 ,add1,add2,cy_r45,cy_r67,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r17,add1,add2,cy_r89,cy_rAB,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r25,add1,add2,cy_rCD,cy_rEF,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			  #endif

			#endif

				#if 1
					i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/
				#else
					__asm	mov		edi, bjmodn0/* bjmodn0 (pointer to int) */
					__asm	mov		ecx, [edi]	/* dereference the pointer */
					__asm	mov		esi, sw
					__asm	mov		edx, i
					__asm	sub		esi, ecx	/* sw - *bjmodn0 */
					__asm	shr		esi, 31		/*	    ((uint32)(sw - *bjmodn0) >> 31);	*/
					__asm	mov		i, esi		/*	i = ((uint32)(sw - *bjmodn0) >> 31);	*/
				#endif

				}
				else	/* Fermat-mod carry in SSE2 mode */
				{
				/* In SSE2 mode, the data are arranged in memory like so, where we view things in 16-byte chunks:

					R0/r1 :	a0.re,b0.re		I0/r2 :	a0.im,b0.im
					R1/r3 :	a1.re,b1.re		I1/r4 :	a1.im,b1.im
					R2/r5 :	a2.re,b2.re		I2/r6 :	a2.im,b2.im
					R3/r7 :	a3.re,b3.re		I3/r8 :	a3.im,b3.im
					R4/r9 :	a4.re,b4.re		I4/r10:	a4.im,b4.im
					R5/r11:	a5.re,b5.re		I5/r12:	a5.im,b5.im
					R6/r13:	a6.re,b6.re		I6/r14:	a6.im,b6.im
					R7/r15:	a7.re,b7.re		I7/r16:	a7.im,b7.im
					R8/r17:	a8.re,b8.re		I8/r18:	a8.im,b8.im
					R9/r19:	a9.re,b9.re		I9/r20:	a9.im,b9.im
					Ra/r21:	aA.re,bA.re		Ia/r22:	aA.im,bA.im
					Rb/r23:	aB.re,bB.re		Ib/r24:	aB.im,bB.im
					Rc/r25:	aC.re,bC.re		Ic/r26:	aC.im,bC.im
					Rd/r27:	aD.re,bD.re		Id/r28:	aD.im,bD.im
					Re/r29:	aE.re,bE.re		Ie/r30:	aE.im,bE.im
					Rf/r31:	aF.re,bF.re		If/r32:	aF.im,bF.im

				Where the R's and I's map to the local temps as follows: R0:f ==> r1:31:2, I0:f ==> r2:32:2 , and the
				a's and b's of each pair represent doubles which in non-SSE2 mode would be getting processed in the same relative
				position on subsequent passes through the main-array for-loop, i.e. among which the carry propagation proceeds as

					a0.re -> b0.re;		a0.im -> b0.im, where these imaginary parts really represent elements
					                                    a0.im = a[n/2] and b0.im = a[n/2+1] of the right-angle transform.

				This data layout is ideal for the negacyclic unweighting/reweighting step bracketing the carry step, but in the latter,
				because of the indesirable intra-xmm-register data dependency this leads to, we instead require data arranged as

					R0/r1 :	a0.re,a0.im		I0/r2 :	b0.re,b0.im, i.e. the non-SSE2 data layout works best in the carry step!

				We need to interleave these pairwise so as to swap the high word of each R-element
				with the low word of the corresponding I-element, e.g. for R0/r1 and I0/r2:

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
					tmp->re = tmp->im = scale;

					/* Get the needed Nth root of -1: */
					add1 = (double *)&rn0[0];
					add2 = (double *)&rn1[0];

					idx_offset = j;
					idx_incr = NDIVR;

		#if 0//def DEBUG_SSE2
		// *** Change stderr ==> dbg_file to dump-to-file ***
			if(j < 4 && full_pass && iter <= 1)
			{
				fprintf(stderr, "Iter = %d, J = %d, Before carries:\n",iter,j);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP00 = %20.2f, %20.2f, %20.2f, %20.2f, CY00 = %20.5f, %20.5f\n",r1 ->re,r1 ->im,r2 ->re,r2 ->im,cy_r01->re,cy_r01->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP01 = %20.2f, %20.2f, %20.2f, %20.2f, CY01 = %20.5f, %20.5f\n",r3 ->re,r3 ->im,r4 ->re,r4 ->im,cy_r23->re,cy_r23->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP02 = %20.2f, %20.2f, %20.2f, %20.2f, CY02 = %20.5f, %20.5f\n",r5 ->re,r5 ->im,r6 ->re,r6 ->im,cy_r45->re,cy_r45->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP03 = %20.2f, %20.2f, %20.2f, %20.2f, CY03 = %20.5f, %20.5f\n",r7 ->re,r7 ->im,r8 ->re,r8 ->im,cy_r67->re,cy_r67->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP04 = %20.2f, %20.2f, %20.2f, %20.2f, CY04 = %20.5f, %20.5f\n",r9 ->re,r9 ->im,r10->re,r10->im,cy_r89->re,cy_r89->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP05 = %20.2f, %20.2f, %20.2f, %20.2f, CY05 = %20.5f, %20.5f\n",r11->re,r11->im,r12->re,r12->im,cy_rAB->re,cy_rAB->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP06 = %20.2f, %20.2f, %20.2f, %20.2f, CY06 = %20.5f, %20.5f\n",r13->re,r13->im,r14->re,r14->im,cy_rCD->re,cy_rCD->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP07 = %20.2f, %20.2f, %20.2f, %20.2f, CY07 = %20.5f, %20.5f\n",r15->re,r15->im,r16->re,r16->im,cy_rEF->re,cy_rEF->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP08 = %20.2f, %20.2f, %20.2f, %20.2f, CY08 = %20.5f, %20.5f\n",r17->re,r17->im,r18->re,r18->im,cy_i01->re,cy_i01->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP09 = %20.2f, %20.2f, %20.2f, %20.2f, CY09 = %20.5f, %20.5f\n",r19->re,r19->im,r20->re,r20->im,cy_i23->re,cy_i23->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP10 = %20.2f, %20.2f, %20.2f, %20.2f, CY10 = %20.5f, %20.5f\n",r21->re,r21->im,r22->re,r22->im,cy_i45->re,cy_i45->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP11 = %20.2f, %20.2f, %20.2f, %20.2f, CY11 = %20.5f, %20.5f\n",r23->re,r23->im,r24->re,r24->im,cy_i67->re,cy_i67->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP12 = %20.2f, %20.2f, %20.2f, %20.2f, CY12 = %20.5f, %20.5f\n",r25->re,r25->im,r26->re,r26->im,cy_i89->re,cy_i89->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP13 = %20.2f, %20.2f, %20.2f, %20.2f, CY13 = %20.5f, %20.5f\n",r27->re,r27->im,r28->re,r28->im,cy_iAB->re,cy_iAB->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP14 = %20.2f, %20.2f, %20.2f, %20.2f, CY14 = %20.5f, %20.5f\n",r29->re,r29->im,r30->re,r30->im,cy_iCD->re,cy_iCD->im);
				fprintf(stderr, "radix16_ditN_cy_dif1: TMP15 = %20.2f, %20.2f, %20.2f, %20.2f, CY15 = %20.5f, %20.5f\n",r31->re,r31->im,r32->re,r32->im,cy_iEF->re,cy_iEF->im);
				fprintf(stderr, "\n");
			}
		#endif

				#if defined(COMPILER_TYPE_MSVC)
				/* The cy_[r|i]_idx[A|B] names here are not meaningful, each simple stores one [re,im] carry pair,
					e.g. cy_r01 stores the carries our of [a0.re,a0.im], cy_r23 stores the carries our of [a1.re,a1.im], etc.
					Here is the actual mapping between these SSE2-mode 2-vector carry pairs and the scalar carries:
					                                      2-vector                               Scalar
					                                     ----------                            ----------- */
					SSE2_fermat_carry_norm_pow2_errcheck(r1 ,cy_r01,idx_offset,idx_incr);	/* cy_r0,cy_i0 */
					SSE2_fermat_carry_norm_pow2_errcheck(r3 ,cy_r23,idx_offset,idx_incr);	/* cy_r1,cy_i1 */
					SSE2_fermat_carry_norm_pow2_errcheck(r5 ,cy_r45,idx_offset,idx_incr);	/* cy_r2,cy_i2 */
					SSE2_fermat_carry_norm_pow2_errcheck(r7 ,cy_r67,idx_offset,idx_incr);	/* cy_r3,cy_i3 */
					SSE2_fermat_carry_norm_pow2_errcheck(r9 ,cy_r89,idx_offset,idx_incr);	/* cy_r4,cy_i4 */
					SSE2_fermat_carry_norm_pow2_errcheck(r11,cy_rAB,idx_offset,idx_incr);	/* cy_r5,cy_i5 */
					SSE2_fermat_carry_norm_pow2_errcheck(r13,cy_rCD,idx_offset,idx_incr);	/* cy_r6,cy_i6 */
					SSE2_fermat_carry_norm_pow2_errcheck(r15,cy_rEF,idx_offset,idx_incr);	/* cy_r7,cy_i7 */
					SSE2_fermat_carry_norm_pow2_errcheck(r17,cy_i01,idx_offset,idx_incr);	/* cy_r8,cy_i8 */
					SSE2_fermat_carry_norm_pow2_errcheck(r19,cy_i23,idx_offset,idx_incr);	/* cy_r9,cy_i9 */
					SSE2_fermat_carry_norm_pow2_errcheck(r21,cy_i45,idx_offset,idx_incr);	/* cy_rA,cy_iA */
					SSE2_fermat_carry_norm_pow2_errcheck(r23,cy_i67,idx_offset,idx_incr);	/* cy_rB,cy_iB */
					SSE2_fermat_carry_norm_pow2_errcheck(r25,cy_i89,idx_offset,idx_incr);	/* cy_rC,cy_iC */
					SSE2_fermat_carry_norm_pow2_errcheck(r27,cy_iAB,idx_offset,idx_incr);	/* cy_rD,cy_iD */
					SSE2_fermat_carry_norm_pow2_errcheck(r29,cy_iCD,idx_offset,idx_incr);	/* cy_rE,cy_iE */
					SSE2_fermat_carry_norm_pow2_errcheck(r31,cy_iEF,idx_offset,idx_incr);	/* cy_rF,cy_iF */

				#else

				  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit 
					SSE2_fermat_carry_norm_pow2_errcheck(r1 ,cy_r01,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r3 ,cy_r23,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r5 ,cy_r45,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r7 ,cy_r67,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r9 ,cy_r89,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r11,cy_rAB,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r13,cy_rCD,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r15,cy_rEF,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r17,cy_i01,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r19,cy_i23,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r21,cy_i45,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r23,cy_i67,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r25,cy_i89,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r27,cy_iAB,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r29,cy_iCD,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r31,cy_iEF,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				  #else
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r1 ,cy_r01,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r5 ,cy_r45,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r9 ,cy_r89,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r13,cy_rCD,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r17,cy_i01,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r21,cy_i45,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r25,cy_i89,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r29,cy_iCD,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				  #endif

				#endif

				}	/* if(MODULUS_TYPE == ...) */

			#endif	/* USE_SSE2 */

		/*...The radix-16 DIF pass is here:	*/

		/* Four DIF radix-4 subconvolution, sans twiddles.	Cost each: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */

		#ifdef USE_SSE2

		  #if defined(COMPILER_TYPE_MSVC)

				SSE2_RADIX4_DIF_IN_PLACE(r1 , r17, r9 , r25)
				SSE2_RADIX4_DIF_IN_PLACE(r5 , r21, r13, r29)
				SSE2_RADIX4_DIF_IN_PLACE(r3 , r19, r11, r27)
				SSE2_RADIX4_DIF_IN_PLACE(r7 , r23, r15, r31)

			/****************************************************************************************
			!...and now do four more radix-4 transforms, including the internal twiddle factors.	!
			!																						!
			!	This is identical to latter half of radix16 DIF, except for the r-vector indexing,	!
			!	which permutes as follows:															!
			!																						!
			!			t1	t3	t5	t7	t9	t11	t13	t15	t17	t19	t21	t23	t25	t27	t29	t31				!
			!		==>	r1	r9	r17	r25	r5	r13	r21	r29	r3	r11	r19	r27	r7	r15	r23	r31				!
			!																						!
			****************************************************************************************/

			/*...Block 1: t1,9,17,25 ==> r1,5,3,7	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
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
				__asm	mov	esi, r1

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

			/*...Block 3: t5,13,21,29 ==> r17,21,19,23	Cost: 16 MOVapd, 26 ADD/SUBpd,  4 MULpd */
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
				__asm	mov	esi, r17

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

			/*...Block 2: t3,11,19,27 ==> r9,13,11,15	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
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
				__asm	mov	esi, r9
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

			/*...Block 4: t7,15,23,31 ==> r25,29,27,31	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
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

				__asm	mov	esi, r25
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

			SSE2_RADIX16_DIF_NOTWIDDLE(add0,p1,p2,p3,p4,r1,r3,r5,r7,r9,r11,r13,r15,r17,r19,r21,r23,r25,r27,r29,r31,isrt2,cc0);

		  #endif

		  #if 0//def DEBUG_SSE2
			if(j < 4 && full_pass && iter <= 100)
			{
				jt = j1;		jp = j2;
				fprintf(stderr, "Iter = %d, full_pass = %d, J = %d, After carries:\n",iter,full_pass,j);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[00] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt   ],a[jp   ],a[jt   +1],a[jp   +1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[01] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p1],a[jp+p1],a[jt+p1+1],a[jp+p1+1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[02] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p2],a[jp+p2],a[jt+p2+1],a[jp+p2+1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[03] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p3],a[jp+p3],a[jt+p3+1],a[jp+p3+1]);			jt = j1+p4;	jp = j2+p4;
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[04] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt   ],a[jp   ],a[jt   +1],a[jp   +1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[05] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p1],a[jp+p1],a[jt+p1+1],a[jp+p1+1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[06] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p2],a[jp+p2],a[jt+p2+1],a[jp+p2+1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[07] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p3],a[jp+p3],a[jt+p3+1],a[jp+p3+1]);			jt = j1+p8;	jp = j2+p8;
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[08] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt   ],a[jp   ],a[jt   +1],a[jp   +1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[09] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p1],a[jp+p1],a[jt+p1+1],a[jp+p1+1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[10] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p2],a[jp+p2],a[jt+p2+1],a[jp+p2+1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[11] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p3],a[jp+p3],a[jt+p3+1],a[jp+p3+1]);			jt = j1+p12;	jp = j2+p12;
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[12] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt   ],a[jp   ],a[jt   +1],a[jp   +1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[13] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p1],a[jp+p1],a[jt+p1+1],a[jp+p1+1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[14] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p2],a[jp+p2],a[jt+p2+1],a[jp+p2+1]);
				fprintf(stderr, "radix16_ditN_cy_dif1: A_out[15] = %20.5f, %20.5f %20.5f, %20.5f\n",a[jt+p3],a[jp+p3],a[jt+p3+1],a[jp+p3+1]);
				fprintf(stderr, "\n");
			}
			if(j == 0 && iter > 100)
			{
				exit(0);
			}
		  #endif

		#else	/* !USE_SSE2 */
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
	#ifdef USE_SSE2
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			_cy_r0[ithread] = cy_r01->re;
			_cy_r1[ithread] = cy_r01->im;
			_cy_r2[ithread] = cy_r23->re;
			_cy_r3[ithread] = cy_r23->im;
			_cy_r4[ithread] = cy_r45->re;
			_cy_r5[ithread] = cy_r45->im;
			_cy_r6[ithread] = cy_r67->re;
			_cy_r7[ithread] = cy_r67->im;
			_cy_r8[ithread] = cy_r89->re;
			_cy_r9[ithread] = cy_r89->im;
			_cy_rA[ithread] = cy_rAB->re;
			_cy_rB[ithread] = cy_rAB->im;
			_cy_rC[ithread] = cy_rCD->re;
			_cy_rD[ithread] = cy_rCD->im;
			_cy_rE[ithread] = cy_rEF->re;
			_cy_rF[ithread] = cy_rEF->im;
		}
		else
		{
			_cy_r0[ithread] = cy_r01->re;	_cy_i0[ithread] = cy_r01->im;
			_cy_r1[ithread] = cy_r23->re;	_cy_i1[ithread] = cy_r23->im;
			_cy_r2[ithread] = cy_r45->re;	_cy_i2[ithread] = cy_r45->im;
			_cy_r3[ithread] = cy_r67->re;	_cy_i3[ithread] = cy_r67->im;
			_cy_r4[ithread] = cy_r89->re;	_cy_i4[ithread] = cy_r89->im;
			_cy_r5[ithread] = cy_rAB->re;	_cy_i5[ithread] = cy_rAB->im;
			_cy_r6[ithread] = cy_rCD->re;	_cy_i6[ithread] = cy_rCD->im;
			_cy_r7[ithread] = cy_rEF->re;	_cy_i7[ithread] = cy_rEF->im;
			_cy_r8[ithread] = cy_i01->re;	_cy_i8[ithread] = cy_i01->im;
			_cy_r9[ithread] = cy_i23->re;	_cy_i9[ithread] = cy_i23->im;
			_cy_rA[ithread] = cy_i45->re;	_cy_iA[ithread] = cy_i45->im;
			_cy_rB[ithread] = cy_i67->re;	_cy_iB[ithread] = cy_i67->im;
			_cy_rC[ithread] = cy_i89->re;	_cy_iC[ithread] = cy_i89->im;
			_cy_rD[ithread] = cy_iAB->re;	_cy_iD[ithread] = cy_iAB->im;
			_cy_rE[ithread] = cy_iCD->re;	_cy_iE[ithread] = cy_iCD->im;
			_cy_rF[ithread] = cy_iEF->re;	_cy_iF[ithread] = cy_iEF->im;
		}
		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
	#else
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
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
		}
		else
		{
			_cy_r0[ithread] = cy_r0;		_cy_i0[ithread] = cy_i0;
			_cy_r1[ithread] = cy_r1;		_cy_i1[ithread] = cy_i1;
			_cy_r2[ithread] = cy_r2;		_cy_i2[ithread] = cy_i2;
			_cy_r3[ithread] = cy_r3;		_cy_i3[ithread] = cy_i3;
			_cy_r4[ithread] = cy_r4;		_cy_i4[ithread] = cy_i4;
			_cy_r5[ithread] = cy_r5;		_cy_i5[ithread] = cy_i5;
			_cy_r6[ithread] = cy_r6;		_cy_i6[ithread] = cy_i6;
			_cy_r7[ithread] = cy_r7;		_cy_i7[ithread] = cy_i7;
			_cy_r8[ithread] = cy_r8;		_cy_i8[ithread] = cy_i8;
			_cy_r9[ithread] = cy_r9;		_cy_i9[ithread] = cy_i9;
			_cy_rA[ithread] = cy_rA;		_cy_iA[ithread] = cy_iA;
			_cy_rB[ithread] = cy_rB;		_cy_iB[ithread] = cy_iB;
			_cy_rC[ithread] = cy_rC;		_cy_iC[ithread] = cy_iC;
			_cy_rD[ithread] = cy_rD;		_cy_iD[ithread] = cy_iD;
			_cy_rE[ithread] = cy_rE;		_cy_iE[ithread] = cy_iE;
			_cy_rF[ithread] = cy_rF;		_cy_iF[ithread] = cy_iF;
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

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	bjmodn0=bjmodn1=bjmodn2=bjmodn3=bjmodn4=bjmodn5=bjmodn6=bjmodn7=bjmodn8=bjmodn9=bjmodnA=bjmodnB=bjmodnC=bjmodnD=bjmodnE=bjmodnF=-1;

/*...change n16 and n_div_wt to non-static to work around a gcc compiler bug. */
	n16   = n/16;
	n_div_nwt = n16 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n16)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/16 in radix16_ditN_cy_dif1_nochk.\n",iter);
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

		ASSERT(HERE, CY_THREADS >= NTHREADS,"radix16_ditN_cy_dif1_nochk: CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"radix16_ditN_cy_dif1_nochk: CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, n16      %CY_THREADS == 0,"radix16_ditN_cy_dif1_nochk: n16      %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"radix16_ditN_cy_dif1_nochk: n_div_nwt%CY_THREADS != 0");
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
			_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int)); if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix16_ditN_cy_dif1_nochk.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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

/* Needed to remove the prefetch-address vars addr & addp for this to compile properly: */
#ifdef USE_OMP
	omp_set_num_threads(CY_THREADS);
  #undef PFETCH
	#pragma omp parallel for private(         i,j,j1,jstart,jhi,k,k1,k2,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,rt,it,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi,temp,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF) default(shared) schedule(static)
#endif

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
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]	*/
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;

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

      for(j=0; j < n16; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		/*       gather the needed data (16 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.	*/
	#if 1
		RADIX_16_DIF(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
					,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,rt,it,c,s)
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
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;

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

      for(j=0; j < n16; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		/*       gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms.	*/
	#if 1
		RADIX_16_DIT(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
					,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,rt,it,c,s)
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

	#ifndef USE_SSE2
		#error pthreaded carry code requires SSE2-enabled build!
	#endif
	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy16_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 16;
		const double crnd = 3.0*0x4000000*0x2000000;
		int idx_offset,idx_incr;
		int j,j1,j2,k;
		int l,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		uint32 p1,p2,p3,p4;
		double *add0, *add1, *add2, *add3;
		struct complex *cc0, *ss0, *isrt2, *max_err, *sse2_rnd, *half_arr, *tmp;
		struct complex *r1,*r2,*r3,*r4,*r5,*r6,*r7,*r8,*r9,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r30,*r31,*r32;
		struct complex *cy_r01,*cy_r23,*cy_r45,*cy_r67,*cy_r89,*cy_rAB,*cy_rCD,*cy_rEF,*cy_i01,*cy_i23,*cy_i45,*cy_i67,*cy_i89,*cy_iAB,*cy_iCD,*cy_iEF;
		int *bjmodn0,*bjmodn1,*bjmodn2,*bjmodn3,*bjmodn4,*bjmodn5,*bjmodn6,*bjmodn7,*bjmodn8,*bjmodn9,*bjmodnA,*bjmodnB,*bjmodnC,*bjmodnD,*bjmodnE,*bjmodnF;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_nm1;

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
		int ithread = thread_arg->tid;	/* unique thread index (use for debug) */

//	printf("cy16_process_chunk: thread %d\n", ithread);

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
		double scale = thread_arg->scale;

	// pointer data:
		double *a = thread_arg->arrdat;	
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;	

		p1 = NDIVR;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );

		r1 = thread_arg->r1;  isrt2 = r1 + 0x20;
		r2  = r1 + 0x01;		cc0 = r1 + 0x21;
		r3  = r1 + 0x02;		ss0 = r1 + 0x22;
		r4  = r1 + 0x03;	 cy_r01 = r1 + 0x23;
		r5  = r1 + 0x04;	 cy_r23 = r1 + 0x24;
		r6  = r1 + 0x05;	 cy_r45 = r1 + 0x25;
		r7  = r1 + 0x06;	 cy_r67 = r1 + 0x26;
		r8  = r1 + 0x07;	 cy_r89 = r1 + 0x27;
		r9  = r1 + 0x08;	 cy_rAB = r1 + 0x28;
		r10 = r1 + 0x09;	 cy_rCD = r1 + 0x29;
		r11 = r1 + 0x0a;	 cy_rEF = r1 + 0x2a;
		r12 = r1 + 0x0b;	 cy_i01 = r1 + 0x2b;
		r13 = r1 + 0x0c;	 cy_i23 = r1 + 0x2c;
		r14 = r1 + 0x0d;	 cy_i45 = r1 + 0x2d;
		r15 = r1 + 0x0e;	 cy_i67 = r1 + 0x2e;
		r16 = r1 + 0x0f;	 cy_i89 = r1 + 0x2f;
		r17 = r1 + 0x10;	 cy_iAB = r1 + 0x30;
		r18 = r1 + 0x11;	 cy_iCD = r1 + 0x31;
		r19 = r1 + 0x12;	 cy_iEF = r1 + 0x32;
		r20 = r1 + 0x13;	max_err = r1 + 0x33;
		r21 = r1 + 0x14;	sse2_rnd= r1 + 0x34;
		r22 = r1 + 0x15;	half_arr= r1 + 0x35;	/* This table needs 20x16 bytes */
		r23 = r1 + 0x16;
		r24 = r1 + 0x17;
		r25 = r1 + 0x18;
		r26 = r1 + 0x19;
		r27 = r1 + 0x1a;
		r28 = r1 + 0x1b;
		r29 = r1 + 0x1c;
		r30 = r1 + 0x1d;
		r31 = r1 + 0x1e;
		r32 = r1 + 0x1f;

		ASSERT(HERE, (isrt2->re == ISRT2 && isrt2->im == ISRT2), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->re == crnd && sse2_rnd->im == crnd), "thread-local memcheck failed!");
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		ASSERT(HERE, (half_arr+10)->re * (half_arr+14)->re == 1.0 && (half_arr+10)->im * (half_arr+14)->im == 1.0, "thread-local memcheck failed!");
	} else {
		ASSERT(HERE, half_arr->re * (half_arr+1)->re == 1.0 && half_arr->im * (half_arr+1)->im == 1.0, "thread-local memcheck failed!");
	}

		max_err->re = 0.0;	max_err->im = 0.0;

		sign_mask = (uint64*)(r1 + radix16_creals_in_local_store);
		sse_bw  = sign_mask + 2;
		sse_sw  = sign_mask + 4;
		sse_nm1 = sign_mask + 6;
		bjmodn0 = (int*)(sign_mask + 8);
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

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */	/* init carries	*/
			*bjmodn0 = thread_arg->bjmodn0;		cy_r01->re = thread_arg->cy_r0;
			*bjmodn1 = thread_arg->bjmodn1;		cy_r01->im = thread_arg->cy_r1;
			*bjmodn2 = thread_arg->bjmodn2;		cy_r23->re = thread_arg->cy_r2;
			*bjmodn3 = thread_arg->bjmodn3;		cy_r23->im = thread_arg->cy_r3;
			*bjmodn4 = thread_arg->bjmodn4;		cy_r45->re = thread_arg->cy_r4;
			*bjmodn5 = thread_arg->bjmodn5;		cy_r45->im = thread_arg->cy_r5;
			*bjmodn6 = thread_arg->bjmodn6;		cy_r67->re = thread_arg->cy_r6;
			*bjmodn7 = thread_arg->bjmodn7;		cy_r67->im = thread_arg->cy_r7;
			*bjmodn8 = thread_arg->bjmodn8;		cy_r89->re = thread_arg->cy_r8;
			*bjmodn9 = thread_arg->bjmodn9;		cy_r89->im = thread_arg->cy_r9;
			*bjmodnA = thread_arg->bjmodnA;		cy_rAB->re = thread_arg->cy_rA;
			*bjmodnB = thread_arg->bjmodnB;		cy_rAB->im = thread_arg->cy_rB;
			*bjmodnC = thread_arg->bjmodnC;		cy_rCD->re = thread_arg->cy_rC;
			*bjmodnD = thread_arg->bjmodnD;		cy_rCD->im = thread_arg->cy_rD;
			*bjmodnE = thread_arg->bjmodnE;		cy_rEF->re = thread_arg->cy_rE;
			*bjmodnF = thread_arg->bjmodnF;		cy_rEF->im = thread_arg->cy_rF;
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
			cy_r01->re = thread_arg->cy_r0;		cy_r01->im = thread_arg->cy_i0;
			cy_r23->re = thread_arg->cy_r1;		cy_r23->im = thread_arg->cy_i1;
			cy_r45->re = thread_arg->cy_r2;		cy_r45->im = thread_arg->cy_i2;
			cy_r67->re = thread_arg->cy_r3;		cy_r67->im = thread_arg->cy_i3;
			cy_r89->re = thread_arg->cy_r4;		cy_r89->im = thread_arg->cy_i4;
			cy_rAB->re = thread_arg->cy_r5;		cy_rAB->im = thread_arg->cy_i5;
			cy_rCD->re = thread_arg->cy_r6;		cy_rCD->im = thread_arg->cy_i6;
			cy_rEF->re = thread_arg->cy_r7;		cy_rEF->im = thread_arg->cy_i7;
			cy_i01->re = thread_arg->cy_r8;		cy_i01->im = thread_arg->cy_i8;
			cy_i23->re = thread_arg->cy_r9;		cy_i23->im = thread_arg->cy_i9;
			cy_i45->re = thread_arg->cy_rA;		cy_i45->im = thread_arg->cy_iA;
			cy_i67->re = thread_arg->cy_rB;		cy_i67->im = thread_arg->cy_iB;
			cy_i89->re = thread_arg->cy_rC;		cy_i89->im = thread_arg->cy_iC;
			cy_iAB->re = thread_arg->cy_rD;		cy_iAB->im = thread_arg->cy_iD;
			cy_iCD->re = thread_arg->cy_rE;		cy_iCD->im = thread_arg->cy_iE;
			cy_iEF->re = thread_arg->cy_rF;		cy_iEF->im = thread_arg->cy_iF;
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += 4)
			{
			/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
			*/
				j1 = (j & mask01) + br4[j&3];
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

				add0 = &a[j1];
				SSE2_RADIX16_DIT_NOTWIDDLE(add0,p1,p2,p3,p4,r1,r3,r5,r7,r9,r11,r17,r25,isrt2,cc0);

			/*...Now do the carries. Since the outputs would
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

					tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
					tmp->re = wtl;		tmp->im = wtl;	++tmp;
					tmp->re = wtn;		tmp->im = wtn;	++tmp;
					tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
					tmp->re = wtnm1;	tmp->im = wtnm1;

					add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
					add2 = &wt1[co2-1];
					add3 = &wt1[co3-1];

				#ifdef ERR_CHECK_ALL
					SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				#else
					SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
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

					tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
					tmp->re = wtl;		tmp->im = wtl;	++tmp;
					tmp->re = wtn;		tmp->im = wtn;	++tmp;
					tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
					tmp->re = wtnm1;	tmp->im = wtnm1;

				/*	i =((uint32)(sw - bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

					co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

					add1 = &wt1[col  ];
					add2 = &wt1[co2-1];

				#ifdef ERR_CHECK_ALL
					SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r1 ,add1,add2,cy_r01,cy_r23,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r9 ,add1,add2,cy_r45,cy_r67,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r17,add1,add2,cy_r89,cy_rAB,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r25,add1,add2,cy_rCD,cy_rEF,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				#else
					SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r1 ,add1,add2,cy_r01,cy_r23,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r9 ,add1,add2,cy_r45,cy_r67,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r17,add1,add2,cy_r89,cy_rAB,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
					SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r25,add1,add2,cy_rCD,cy_rEF,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				#endif

					i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/
				}
				else	/* Fermat-mod carry in SSE2 mode */
				{
					tmp = half_arr+2;
					tmp->re = tmp->im = scale;

					/* Get the needed Nth root of -1: */
					add1 = (double *)&rn0[0];
					add2 = (double *)&rn1[0];

					idx_offset = j;
					idx_incr = NDIVR;

				  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit 
					SSE2_fermat_carry_norm_pow2_errcheck(r1 ,cy_r01,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r3 ,cy_r23,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r5 ,cy_r45,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r7 ,cy_r67,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r9 ,cy_r89,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r11,cy_rAB,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r13,cy_rCD,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r15,cy_rEF,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r17,cy_i01,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r19,cy_i23,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r21,cy_i45,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r23,cy_i67,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r25,cy_i89,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r27,cy_iAB,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r29,cy_iCD,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck(r31,cy_iEF,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				  #else
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r1 ,cy_r01,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r5 ,cy_r45,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r9 ,cy_r89,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r13,cy_rCD,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r17,cy_i01,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r21,cy_i45,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r25,cy_i89,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
					SSE2_fermat_carry_norm_pow2_errcheck_X2(r29,cy_iCD,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				  #endif
				}	/* if(MODULUS_TYPE == ...) */

			/*...The radix-16 DIF pass is here:	*/
	
				SSE2_RADIX16_DIF_NOTWIDDLE(add0,p1,p2,p3,p4,r1,r3,r5,r7,r9,r11,r13,r15,r17,r19,r21,r23,r25,r27,r29,r31,isrt2,cc0);

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
			thread_arg->cy_r0 = cy_r01->re;
			thread_arg->cy_r1 = cy_r01->im;
			thread_arg->cy_r2 = cy_r23->re;
			thread_arg->cy_r3 = cy_r23->im;
			thread_arg->cy_r4 = cy_r45->re;
			thread_arg->cy_r5 = cy_r45->im;
			thread_arg->cy_r6 = cy_r67->re;
			thread_arg->cy_r7 = cy_r67->im;
			thread_arg->cy_r8 = cy_r89->re;
			thread_arg->cy_r9 = cy_r89->im;
			thread_arg->cy_rA = cy_rAB->re;
			thread_arg->cy_rB = cy_rAB->im;
			thread_arg->cy_rC = cy_rCD->re;
			thread_arg->cy_rD = cy_rCD->im;
			thread_arg->cy_rE = cy_rEF->re;
			thread_arg->cy_rF = cy_rEF->im;
		}
		else
		{
			thread_arg->cy_r0 = cy_r01->re;		thread_arg->cy_i0 = cy_r01->im;
			thread_arg->cy_r1 = cy_r23->re;		thread_arg->cy_i1 = cy_r23->im;
			thread_arg->cy_r2 = cy_r45->re;		thread_arg->cy_i2 = cy_r45->im;
			thread_arg->cy_r3 = cy_r67->re;		thread_arg->cy_i3 = cy_r67->im;
			thread_arg->cy_r4 = cy_r89->re;		thread_arg->cy_i4 = cy_r89->im;
			thread_arg->cy_r5 = cy_rAB->re;		thread_arg->cy_i5 = cy_rAB->im;
			thread_arg->cy_r6 = cy_rCD->re;		thread_arg->cy_i6 = cy_rCD->im;
			thread_arg->cy_r7 = cy_rEF->re;		thread_arg->cy_i7 = cy_rEF->im;
			thread_arg->cy_r8 = cy_i01->re;		thread_arg->cy_i8 = cy_i01->im;
			thread_arg->cy_r9 = cy_i23->re;		thread_arg->cy_i9 = cy_i23->im;
			thread_arg->cy_rA = cy_i45->re;		thread_arg->cy_iA = cy_i45->im;
			thread_arg->cy_rB = cy_i67->re;		thread_arg->cy_iB = cy_i67->im;
			thread_arg->cy_rC = cy_i89->re;		thread_arg->cy_iC = cy_i89->im;
			thread_arg->cy_rD = cy_iAB->re;		thread_arg->cy_iD = cy_iAB->im;
			thread_arg->cy_rE = cy_iCD->re;		thread_arg->cy_iE = cy_iCD->im;
			thread_arg->cy_rF = cy_iEF->re;		thread_arg->cy_iF = cy_iEF->im;
		}
		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}

		return 0x0;
	}
#endif

