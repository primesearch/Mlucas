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

#ifdef CTIME	// define at compile time to enable internal timing diagnostics
	double dt_fwd, dt_inv, dt_cy, dt_tot;
	clock_t clock1, clock2, clock3;
#endif

#ifdef USE_SSE2

	const int radix36_creals_in_local_store = 216;

  #ifdef USE_PTHREAD

	#ifndef USE_SSE2
		#error Pthreading only available in SSE2 mode!
	#endif

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

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *r00;
		struct complex *half_arr;

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
		int bjmodn24;
		int bjmodn25;
		int bjmodn26;
		int bjmodn27;
		int bjmodn28;
		int bjmodn29;
		int bjmodn30;
		int bjmodn31;
		int bjmodn32;
		int bjmodn33;
		int bjmodn34;
		int bjmodn35;
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
		double cy24;
		double cy25;
		double cy26;
		double cy27;
		double cy28;
		double cy29;
		double cy30;
		double cy31;
		double cy32;
		double cy33;
		double cy34;
		double cy35;
	};

  #endif

//	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */

	#if defined(COMPILER_TYPE_MSVC)

		#include "sse2_macro.h"

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix36_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix36_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

#endif

/**************/

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
	const uint32 RADIX = 36;
	const double crnd = 3.0*0x4000000*0x2000000;
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35;
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
		,temp,frac,scale;
	double maxerr = 0.0;

	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;
	double wt_re,wt_im;									/* Fermat-mod weights stuff */
	int ii00,ii01,ii02,ii03,ii04,ii05,ii06,ii07,ii08,ii09,ii10,ii11,ii12,ii13,ii14,ii15,ii16,ii17,ii18,ii19,ii20,ii21,ii22,ii23,ii24,ii25,ii26,ii27,ii28,ii29,ii30,ii31,ii32,ii33,ii34,ii35;	/* indices into weights arrays (mod NWT) */

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static int cslots_in_local_store;
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
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
		static task_control_t   task_control = {NULL, (void*)cy36_process_chunk, NULL, 0x0};
	#endif

  #else
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7, *add8;	/* Addresses into array sections */
  #endif

	static struct complex *isrt2, *cc1, *ss1, *cc2, *ss2, *cc3m1, *ss3, *cc4, *ss4, *max_err, *sse2_rnd, *half_arr, *tmp;
	static struct complex
	 *r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0a,*r0b,*r0c,*r0d,*r0e,*r0f,*r0g,*r0h
	,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1a,*r1b,*r1c,*r1d,*r1e,*r1f,*r1g,*r1h
	,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2a,*r2b,*r2c,*r2d,*r2e,*r2f,*r2g,*r2h
	,*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3a,*r3b,*r3c,*r3d,*r3e,*r3f,*r3g,*r3h
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r
	,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i,*s1p20i,*s1p21i,*s1p22i,*s1p23i,*s1p24i,*s1p25i,*s1p26i,*s1p27i,*s1p28i,*s1p29i,*s1p30i,*s1p31i,*s1p32i,*s1p33i,*s1p34i,*s1p35i;
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35;
	static struct complex *cy_r00,*cy_r02,*cy_r04,*cy_r06,*cy_r08,*cy_r10,*cy_r12,*cy_r14,*cy_r16,*cy_r18,*cy_r20,*cy_r22,*cy_r24,*cy_r26,*cy_r28,*cy_r30,*cy_r32,*cy_r34;
	static struct complex *cy_i00,*cy_i02,*cy_i04,*cy_i06,*cy_i08,*cy_i10,*cy_i12,*cy_i14,*cy_i16,*cy_i18,*cy_i20,*cy_i22,*cy_i24,*cy_i26,*cy_i28,*cy_i30,*cy_i32,*cy_i34;

#else

	double wt,wtinv,wtA,wtB,wtC;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35;
	double
	 a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i
	,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35
	,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,cy_i28,cy_i29,cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0,*_bjmodn28 = 0x0,*_bjmodn29 = 0x0,*_bjmodn30 = 0x0,*_bjmodn31 = 0x0,*_bjmodn32 = 0x0,*_bjmodn33 = 0x0,*_bjmodn34 = 0x0,*_bjmodn35 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy_r00 = 0x0,*_cy_r01 = 0x0,*_cy_r02 = 0x0,*_cy_r03 = 0x0,*_cy_r04 = 0x0,*_cy_r05 = 0x0,*_cy_r06 = 0x0,*_cy_r07 = 0x0,*_cy_r08 = 0x0,*_cy_r09 = 0x0,*_cy_r10 = 0x0,*_cy_r11 = 0x0,*_cy_r12 = 0x0,*_cy_r13 = 0x0,*_cy_r14 = 0x0,*_cy_r15 = 0x0,*_cy_r16 = 0x0,*_cy_r17 = 0x0,*_cy_r18 = 0x0,*_cy_r19 = 0x0,*_cy_r20 = 0x0,*_cy_r21 = 0x0,*_cy_r22 = 0x0,*_cy_r23 = 0x0,*_cy_r24 = 0x0,*_cy_r25 = 0x0,*_cy_r26 = 0x0,*_cy_r27 = 0x0,*_cy_r28 = 0x0,*_cy_r29 = 0x0,*_cy_r30 = 0x0,*_cy_r31 = 0x0,*_cy_r32 = 0x0,*_cy_r33 = 0x0,*_cy_r34 = 0x0,*_cy_r35 = 0x0,
	*_cy_i00 = 0x0,*_cy_i01 = 0x0,*_cy_i02 = 0x0,*_cy_i03 = 0x0,*_cy_i04 = 0x0,*_cy_i05 = 0x0,*_cy_i06 = 0x0,*_cy_i07 = 0x0,*_cy_i08 = 0x0,*_cy_i09 = 0x0,*_cy_i10 = 0x0,*_cy_i11 = 0x0,*_cy_i12 = 0x0,*_cy_i13 = 0x0,*_cy_i14 = 0x0,*_cy_i15 = 0x0,*_cy_i16 = 0x0,*_cy_i17 = 0x0,*_cy_i18 = 0x0,*_cy_i19 = 0x0,*_cy_i20 = 0x0,*_cy_i21 = 0x0,*_cy_i22 = 0x0,*_cy_i23 = 0x0,*_cy_i24 = 0x0,*_cy_i25 = 0x0,*_cy_i26 = 0x0,*_cy_i27 = 0x0,*_cy_i28 = 0x0,*_cy_i29 = 0x0,*_cy_i30 = 0x0,*_cy_i31 = 0x0,*_cy_i32 = 0x0,*_cy_i33 = 0x0,*_cy_i34 = 0x0,*_cy_i35 = 0x0;

#ifdef CTIME
	const double ICPS = 1.0/CLOCKS_PER_SEC;
	clock1 = clock();
	dt_fwd = dt_inv = dt_cy = dt_tot = 0.0;
#endif

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

			pool_work_units = CY_THREADS;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

		#endif

		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);

	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE), "SSE2 currently only supports Mersenne-mod!");
		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix36_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_COMPLEX(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix36_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 192 16-byte slots of sc_arr for r-and-s temporaries, next 7 for the nontrivial complex 16th roots,
	next 36 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
									tmp	= sc_ptr + 0x48;
		r00	= sc_ptr + 0x00;		s1p00r = tmp + 0x00;	cc1		= tmp + 0x48;
		r01	= sc_ptr + 0x01;		s1p00i = tmp + 0x01;	ss1		= tmp + 0x49;
		r02	= sc_ptr + 0x02;		s1p01r = tmp + 0x02;	cc2		= tmp + 0x4a;
		r03	= sc_ptr + 0x03;		s1p01i = tmp + 0x03;	ss2		= tmp + 0x4b;
		r04	= sc_ptr + 0x04;		s1p02r = tmp + 0x04;	cc3m1	= tmp + 0x4c;
		r05	= sc_ptr + 0x05;		s1p02i = tmp + 0x05;	ss3		= tmp + 0x4d;
		r06	= sc_ptr + 0x06;		s1p03r = tmp + 0x06;	cc4		= tmp + 0x4e;
		r07	= sc_ptr + 0x07;		s1p03i = tmp + 0x07;	ss4		= tmp + 0x4f;
		r08	= sc_ptr + 0x08;		s1p04r = tmp + 0x08;	cy_r00	= tmp + 0x50;
		r09	= sc_ptr + 0x09;		s1p04i = tmp + 0x09;	cy_r02	= tmp + 0x51;
		r0a	= sc_ptr + 0x0a;		s1p05r = tmp + 0x0a;	cy_r04	= tmp + 0x52;
		r0b	= sc_ptr + 0x0b;		s1p05i = tmp + 0x0b;	cy_r06	= tmp + 0x53;
		r0c	= sc_ptr + 0x0c;		s1p06r = tmp + 0x0c;	cy_r08	= tmp + 0x54;
		r0d	= sc_ptr + 0x0d;		s1p06i = tmp + 0x0d;	cy_r10	= tmp + 0x55;
		r0e	= sc_ptr + 0x0e;		s1p07r = tmp + 0x0e;	cy_r12	= tmp + 0x56;
		r0f	= sc_ptr + 0x0f;		s1p07i = tmp + 0x0f;	cy_r14	= tmp + 0x57;
		r0g	= sc_ptr + 0x10;		s1p08r = tmp + 0x10;	cy_r16	= tmp + 0x58;
		r0h	= sc_ptr + 0x11;		s1p08i = tmp + 0x11;	cy_r18	= tmp + 0x59;
		r10	= sc_ptr + 0x12;		s1p09r = tmp + 0x12;	cy_r20	= tmp + 0x5a;
		r11	= sc_ptr + 0x13;		s1p09i = tmp + 0x13;	cy_r22	= tmp + 0x5b;
		r12	= sc_ptr + 0x14;		s1p10r = tmp + 0x14;	cy_r24	= tmp + 0x5c;
		r13	= sc_ptr + 0x15;		s1p10i = tmp + 0x15;	cy_r26	= tmp + 0x5d;
		r14	= sc_ptr + 0x16;		s1p11r = tmp + 0x16;	cy_r28	= tmp + 0x5e;
		r15	= sc_ptr + 0x17;		s1p11i = tmp + 0x17;	cy_r30	= tmp + 0x5f;
		r16	= sc_ptr + 0x18;		s1p12r = tmp + 0x18;	cy_r32	= tmp + 0x60;
		r17	= sc_ptr + 0x19;		s1p12i = tmp + 0x19;	cy_r34	= tmp + 0x61;
		r18	= sc_ptr + 0x1a;		s1p13r = tmp + 0x1a;	cy_i00	= tmp + 0x62;
		r19	= sc_ptr + 0x1b;		s1p13i = tmp + 0x1b;	cy_i02	= tmp + 0x63;
		r1a	= sc_ptr + 0x1c;		s1p14r = tmp + 0x1c;	cy_i04	= tmp + 0x64;
		r1b	= sc_ptr + 0x1d;		s1p14i = tmp + 0x1d;	cy_i06	= tmp + 0x65;
		r1c	= sc_ptr + 0x1e;		s1p15r = tmp + 0x1e;	cy_i08	= tmp + 0x66;
		r1d	= sc_ptr + 0x1f;		s1p15i = tmp + 0x1f;	cy_i10	= tmp + 0x67;
		r1e	= sc_ptr + 0x20;		s1p16r = tmp + 0x20;	cy_i12	= tmp + 0x68;
		r1f	= sc_ptr + 0x21;		s1p16i = tmp + 0x21;	cy_i14	= tmp + 0x69;
		r1g	= sc_ptr + 0x22;		s1p17r = tmp + 0x22;	cy_i16	= tmp + 0x6a;
		r1h	= sc_ptr + 0x23;		s1p17i = tmp + 0x23;	cy_i18	= tmp + 0x6b;
		r20	= sc_ptr + 0x24;		s1p18r = tmp + 0x24;	cy_i20	= tmp + 0x6c;
		r21	= sc_ptr + 0x25;		s1p18i = tmp + 0x25;	cy_i22	= tmp + 0x6d;
		r22	= sc_ptr + 0x26;		s1p19r = tmp + 0x26;	cy_i24	= tmp + 0x6e;
		r23	= sc_ptr + 0x27;		s1p19i = tmp + 0x27;	cy_i26	= tmp + 0x6f;
		r24	= sc_ptr + 0x28;		s1p20r = tmp + 0x28;	cy_i28	= tmp + 0x70;
		r25	= sc_ptr + 0x29;		s1p20i = tmp + 0x29;	cy_i30	= tmp + 0x71;
		r26	= sc_ptr + 0x2a;		s1p21r = tmp + 0x2a;	cy_i32	= tmp + 0x72;
		r27	= sc_ptr + 0x2b;		s1p21i = tmp + 0x2b;	cy_i34	= tmp + 0x73;
		r28	= sc_ptr + 0x2c;		s1p22r = tmp + 0x2c;	max_err = tmp + 0x74;
		r29	= sc_ptr + 0x2d;		s1p22i = tmp + 0x2d;	sse2_rnd= tmp + 0x75;
		r2a	= sc_ptr + 0x2e;		s1p23r = tmp + 0x2e;	half_arr= tmp + 0x76;	/* This table needs 20x16 bytes */
		r2b	= sc_ptr + 0x2f;		s1p23i = tmp + 0x2f;
		r2c	= sc_ptr + 0x30;		s1p24r = tmp + 0x30;
		r2d	= sc_ptr + 0x31;		s1p24i = tmp + 0x31;
		r2e	= sc_ptr + 0x32;		s1p25r = tmp + 0x32;
		r2f	= sc_ptr + 0x33;		s1p25i = tmp + 0x33;
		r2g	= sc_ptr + 0x34;		s1p26r = tmp + 0x34;
		r2h	= sc_ptr + 0x35;		s1p26i = tmp + 0x35;
		r30	= sc_ptr + 0x36;		s1p27r = tmp + 0x36;
		r31	= sc_ptr + 0x37;		s1p27i = tmp + 0x37;
		r32	= sc_ptr + 0x38;		s1p28r = tmp + 0x38;
		r33	= sc_ptr + 0x39;		s1p28i = tmp + 0x39;
		r34	= sc_ptr + 0x3a;		s1p29r = tmp + 0x3a;
		r35	= sc_ptr + 0x3b;		s1p29i = tmp + 0x3b;
		r36	= sc_ptr + 0x3c;		s1p30r = tmp + 0x3c;
		r37	= sc_ptr + 0x3d;		s1p30i = tmp + 0x3d;
		r38	= sc_ptr + 0x3e;		s1p31r = tmp + 0x3e;
		r39	= sc_ptr + 0x3f;		s1p31i = tmp + 0x3f;
		r3a	= sc_ptr + 0x40;		s1p32r = tmp + 0x40;
		r3b	= sc_ptr + 0x41;		s1p32i = tmp + 0x41;
		r3c	= sc_ptr + 0x42;		s1p33r = tmp + 0x42;
		r3d	= sc_ptr + 0x43;		s1p33i = tmp + 0x43;
		r3e	= sc_ptr + 0x44;		s1p34r = tmp + 0x44;
		r3f	= sc_ptr + 0x45;		s1p34i = tmp + 0x45;
		r3g	= sc_ptr + 0x46;		s1p35r = tmp + 0x46;
		r3h	= sc_ptr + 0x47;		s1p35i = tmp + 0x47;
		ASSERT(HERE, (radix36_creals_in_local_store << 4) >= ((long)half_arr - (long)r00) + (20 << 4), "radix36_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		cc1  ->re = cc1  ->im = c	;		ss1->re = ss1->im = s ;
		cc2  ->re = cc2  ->im = c2  ;		ss2->re = ss2->im = s2;
		cc3m1->re = cc3m1->im = c3m1;		ss3->re = ss3->im = s3;
		cc4  ->re = cc4  ->im = c4  ;		ss4->re = ss4->im = s4;

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
		/* Forward-weight multipliers: */
		tmp->re = 1.0;	tmp->im = 1.0;	++tmp;
		tmp->re = .50;	tmp->im = 1.0;	++tmp;
		tmp->re = 1.0;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .50;	++tmp;
		/* Inverse-weight multipliers: */
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

		/* Floating-point sign mask used for FABS on packed doubles: */
		sign_mask = sm_ptr;
		*sign_mask++ = (uint64)0x7FFFFFFFFFFFFFFFull;
		*sign_mask-- = (uint64)0x7FFFFFFFFFFFFFFFull;

	#if 0
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

		sse_n   = sm_ptr + 6;
		__asm	lea	eax, n
		__asm	mov	ebx, sse_n
		__asm	movd	xmm0,[eax]
		__asm	pshufd	xmm0,xmm0,0	// Broadcast low 32 bits of xmm0 to all 4 slots of xmm0
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

		sse_n   = sm_ptr + 6;
		tmp64 = (uint64)n;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_n++ = tmp64;
		*sse_n-- = tmp64;
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
		tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
		tdat[ithread].half_arr = (long)tdat[ithread].r00 + ((long)half_arr - (long)r00);
	}
#endif

		bjmodn00 = (uint32*)(sm_ptr + 8);
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
		bjmodn24 = bjmodn23 + 1;
		bjmodn25 = bjmodn24 + 1;
		bjmodn26 = bjmodn25 + 1;
		bjmodn27 = bjmodn26 + 1;
		bjmodn28 = bjmodn27 + 1;
		bjmodn29 = bjmodn28 + 1;
		bjmodn30 = bjmodn29 + 1;
		bjmodn31 = bjmodn30 + 1;
		bjmodn32 = bjmodn31 + 1;
		bjmodn33 = bjmodn32 + 1;
		bjmodn34 = bjmodn33 + 1;
		bjmodn35 = bjmodn34 + 1;

	  #ifdef USE_PTHREAD
		tmp = __r0 + cslots_in_local_store;
		/* Init thread 1-CY_THREADS's local stores and pointers: */
		for(i = 1; i < CY_THREADS; ++i) {
			/* Only care about the constants for each thread here, but easier to just copy the entire thread0 local store: */
			memcpy(tmp, __r0, cslots_in_local_store<<4);	// bytewise copy treats complex and uint64 subdata the same
			tmp += cslots_in_local_store;
		}
	  #endif

	#endif	// USE_SSE2

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
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

		ASSERT(HERE, p01+p01 == p02, "p01+p01 != p02");
		ASSERT(HERE, p02+p02 == p04, "p02+p02 != p04");
		ASSERT(HERE, p04+p04 == p08, "p04+p04 != p08");
		ASSERT(HERE, p08+p04 == p12, "p08+p04 != p12");
		ASSERT(HERE, p12+p04 == p16, "p12+p04 != p16");
		ASSERT(HERE, p16+p04 == p20, "p16+p04 != p20");
		ASSERT(HERE, p20+p04 == p24, "p20+p04 != p24");
		ASSERT(HERE, p24+p04 == p28, "p24+p04 != p28");
		ASSERT(HERE, p28+p04 == p32, "p28+p04 != p32");

		if(_cy_r00)	/* If it's a new exponent of a range test, need to deallocate these. */
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
			free((void *)_bjmodn24); _bjmodn24 = 0x0;
			free((void *)_bjmodn25); _bjmodn25 = 0x0;
			free((void *)_bjmodn26); _bjmodn26 = 0x0;
			free((void *)_bjmodn27); _bjmodn27 = 0x0;
			free((void *)_bjmodn28); _bjmodn28 = 0x0;
			free((void *)_bjmodn29); _bjmodn29 = 0x0;
			free((void *)_bjmodn30); _bjmodn30 = 0x0;
			free((void *)_bjmodn31); _bjmodn31 = 0x0;
			free((void *)_bjmodn32); _bjmodn32 = 0x0;
			free((void *)_bjmodn33); _bjmodn33 = 0x0;
			free((void *)_bjmodn34); _bjmodn34 = 0x0;
			free((void *)_bjmodn35); _bjmodn35 = 0x0;

			free((void *)_cy_r00); _cy_r00 = 0x0;		free((void *)_cy_i00); _cy_i00 = 0x0;
			free((void *)_cy_r01); _cy_r01 = 0x0;		free((void *)_cy_i01); _cy_i01 = 0x0;
			free((void *)_cy_r02); _cy_r02 = 0x0;		free((void *)_cy_i02); _cy_i02 = 0x0;
			free((void *)_cy_r03); _cy_r03 = 0x0;		free((void *)_cy_i03); _cy_i03 = 0x0;
			free((void *)_cy_r04); _cy_r04 = 0x0;		free((void *)_cy_i04); _cy_i04 = 0x0;
			free((void *)_cy_r05); _cy_r05 = 0x0;		free((void *)_cy_i05); _cy_i05 = 0x0;
			free((void *)_cy_r06); _cy_r06 = 0x0;		free((void *)_cy_i06); _cy_i06 = 0x0;
			free((void *)_cy_r07); _cy_r07 = 0x0;		free((void *)_cy_i07); _cy_i07 = 0x0;
			free((void *)_cy_r08); _cy_r08 = 0x0;		free((void *)_cy_i08); _cy_i08 = 0x0;
			free((void *)_cy_r09); _cy_r09 = 0x0;		free((void *)_cy_i09); _cy_i09 = 0x0;
			free((void *)_cy_r10); _cy_r10 = 0x0;		free((void *)_cy_i10); _cy_i10 = 0x0;
			free((void *)_cy_r11); _cy_r11 = 0x0;		free((void *)_cy_i11); _cy_i11 = 0x0;
			free((void *)_cy_r12); _cy_r12 = 0x0;		free((void *)_cy_i12); _cy_i12 = 0x0;
			free((void *)_cy_r13); _cy_r13 = 0x0;		free((void *)_cy_i13); _cy_i13 = 0x0;
			free((void *)_cy_r14); _cy_r14 = 0x0;		free((void *)_cy_i14); _cy_i14 = 0x0;
			free((void *)_cy_r15); _cy_r15 = 0x0;		free((void *)_cy_i15); _cy_i15 = 0x0;
			free((void *)_cy_r16); _cy_r16 = 0x0;		free((void *)_cy_i16); _cy_i16 = 0x0;
			free((void *)_cy_r17); _cy_r17 = 0x0;		free((void *)_cy_i17); _cy_i17 = 0x0;
			free((void *)_cy_r18); _cy_r18 = 0x0;		free((void *)_cy_i18); _cy_i18 = 0x0;
			free((void *)_cy_r19); _cy_r19 = 0x0;		free((void *)_cy_i19); _cy_i19 = 0x0;
			free((void *)_cy_r20); _cy_r20 = 0x0;		free((void *)_cy_i20); _cy_i20 = 0x0;
			free((void *)_cy_r21); _cy_r21 = 0x0;		free((void *)_cy_i21); _cy_i21 = 0x0;
			free((void *)_cy_r22); _cy_r22 = 0x0;		free((void *)_cy_i22); _cy_i22 = 0x0;
			free((void *)_cy_r23); _cy_r23 = 0x0;		free((void *)_cy_i23); _cy_i23 = 0x0;
			free((void *)_cy_r24); _cy_r24 = 0x0;		free((void *)_cy_i24); _cy_i24 = 0x0;
			free((void *)_cy_r25); _cy_r25 = 0x0;		free((void *)_cy_i25); _cy_i25 = 0x0;
			free((void *)_cy_r26); _cy_r26 = 0x0;		free((void *)_cy_i26); _cy_i26 = 0x0;
			free((void *)_cy_r27); _cy_r27 = 0x0;		free((void *)_cy_i27); _cy_i27 = 0x0;
			free((void *)_cy_r28); _cy_r28 = 0x0;		free((void *)_cy_i28); _cy_i28 = 0x0;
			free((void *)_cy_r29); _cy_r29 = 0x0;		free((void *)_cy_i29); _cy_i29 = 0x0;
			free((void *)_cy_r30); _cy_r30 = 0x0;		free((void *)_cy_i30); _cy_i30 = 0x0;
			free((void *)_cy_r31); _cy_r31 = 0x0;		free((void *)_cy_i31); _cy_i31 = 0x0;
			free((void *)_cy_r32); _cy_r32 = 0x0;		free((void *)_cy_i32); _cy_i32 = 0x0;
			free((void *)_cy_r33); _cy_r33 = 0x0;		free((void *)_cy_i33); _cy_i33 = 0x0;
			free((void *)_cy_r34); _cy_r34 = 0x0;		free((void *)_cy_i34); _cy_i34 = 0x0;
			free((void *)_cy_r35); _cy_r35 = 0x0;		free((void *)_cy_i35); _cy_i35 = 0x0;

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
		_bjmodn24	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_bjmodn28	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn28== 0x0);
		_bjmodn29	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn29== 0x0);
		_bjmodn30	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn30== 0x0);
		_bjmodn31	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn31== 0x0);
		_bjmodn32	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn32== 0x0);
		_bjmodn33	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn33== 0x0);
		_bjmodn34	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn34== 0x0);
		_bjmodn35	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn35== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy_r00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r00== 0x0);
		_cy_r01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r01== 0x0);
		_cy_r02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r02== 0x0);
		_cy_r03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r03== 0x0);
		_cy_r04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r04== 0x0);
		_cy_r05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r05== 0x0);
		_cy_r06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r06== 0x0);
		_cy_r07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r07== 0x0);
		_cy_r08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r08== 0x0);
		_cy_r09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r09== 0x0);
		_cy_r10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r10== 0x0);
		_cy_r11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r11== 0x0);
		_cy_r12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r12== 0x0);
		_cy_r13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r13== 0x0);
		_cy_r14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r14== 0x0);
		_cy_r15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r15== 0x0);
		_cy_r16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r16== 0x0);
		_cy_r17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r17== 0x0);
		_cy_r18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r18== 0x0);
		_cy_r19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r19== 0x0);
		_cy_r20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r20== 0x0);
		_cy_r21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r21== 0x0);
		_cy_r22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r22== 0x0);
		_cy_r23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r23== 0x0);
		_cy_r24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r24== 0x0);
		_cy_r25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r25== 0x0);
		_cy_r26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r26== 0x0);
		_cy_r27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r27== 0x0);
		_cy_r28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r28== 0x0);
		_cy_r29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r29== 0x0);
		_cy_r30	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r30== 0x0);
		_cy_r31	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r31== 0x0);
		_cy_r32	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r32== 0x0);
		_cy_r33	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r33== 0x0);
		_cy_r34	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r34== 0x0);
		_cy_r35	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r35== 0x0);

		_cy_i00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i00== 0x0);
		_cy_i01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i01== 0x0);
		_cy_i02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i02== 0x0);
		_cy_i03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i03== 0x0);
		_cy_i04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i04== 0x0);
		_cy_i05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i05== 0x0);
		_cy_i06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i06== 0x0);
		_cy_i07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i07== 0x0);
		_cy_i08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i08== 0x0);
		_cy_i09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i09== 0x0);
		_cy_i10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i10== 0x0);
		_cy_i11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i11== 0x0);
		_cy_i12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i12== 0x0);
		_cy_i13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i13== 0x0);
		_cy_i14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i14== 0x0);
		_cy_i15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i15== 0x0);
		_cy_i16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i16== 0x0);
		_cy_i17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i17== 0x0);
		_cy_i18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i18== 0x0);
		_cy_i19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i19== 0x0);
		_cy_i20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i20== 0x0);
		_cy_i21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i21== 0x0);
		_cy_i22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i22== 0x0);
		_cy_i23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i23== 0x0);
		_cy_i24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i24== 0x0);
		_cy_i25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i25== 0x0);
		_cy_i26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i26== 0x0);
		_cy_i27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i27== 0x0);
		_cy_i28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i28== 0x0);
		_cy_i29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i29== 0x0);
		_cy_i30	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i30== 0x0);
		_cy_i31	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i31== 0x0);
		_cy_i32	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i32== 0x0);
		_cy_i33	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i33== 0x0);
		_cy_i34	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i34== 0x0);
		_cy_i35	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i35== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix36_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/36-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix36_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jhi = NDIVR/CY_THREADS;
		}
		else
		{
			jhi = NDIVR/CY_THREADS/2;
		}

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

/*...The radix-36 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r00[ithread] = 0;	_cy_i00[ithread] = 0;
		_cy_r01[ithread] = 0;	_cy_i01[ithread] = 0;
		_cy_r02[ithread] = 0;	_cy_i02[ithread] = 0;
		_cy_r03[ithread] = 0;	_cy_i03[ithread] = 0;
		_cy_r04[ithread] = 0;	_cy_i04[ithread] = 0;
		_cy_r05[ithread] = 0;	_cy_i05[ithread] = 0;
		_cy_r06[ithread] = 0;	_cy_i06[ithread] = 0;
		_cy_r07[ithread] = 0;	_cy_i07[ithread] = 0;
		_cy_r08[ithread] = 0;	_cy_i08[ithread] = 0;
		_cy_r09[ithread] = 0;	_cy_i09[ithread] = 0;
		_cy_r10[ithread] = 0;	_cy_i10[ithread] = 0;
		_cy_r11[ithread] = 0;	_cy_i11[ithread] = 0;
		_cy_r12[ithread] = 0;	_cy_i12[ithread] = 0;
		_cy_r13[ithread] = 0;	_cy_i13[ithread] = 0;
		_cy_r14[ithread] = 0;	_cy_i14[ithread] = 0;
		_cy_r15[ithread] = 0;	_cy_i15[ithread] = 0;
		_cy_r16[ithread] = 0;	_cy_i16[ithread] = 0;
		_cy_r17[ithread] = 0;	_cy_i17[ithread] = 0;
		_cy_r18[ithread] = 0;	_cy_i18[ithread] = 0;
		_cy_r19[ithread] = 0;	_cy_i19[ithread] = 0;
		_cy_r20[ithread] = 0;	_cy_i20[ithread] = 0;
		_cy_r21[ithread] = 0;	_cy_i21[ithread] = 0;
		_cy_r22[ithread] = 0;	_cy_i22[ithread] = 0;
		_cy_r23[ithread] = 0;	_cy_i23[ithread] = 0;
		_cy_r24[ithread] = 0;	_cy_i24[ithread] = 0;
		_cy_r25[ithread] = 0;	_cy_i25[ithread] = 0;
		_cy_r26[ithread] = 0;	_cy_i26[ithread] = 0;
		_cy_r27[ithread] = 0;	_cy_i27[ithread] = 0;
		_cy_r28[ithread] = 0;	_cy_i28[ithread] = 0;
		_cy_r29[ithread] = 0;	_cy_i29[ithread] = 0;
		_cy_r30[ithread] = 0;	_cy_i30[ithread] = 0;
		_cy_r31[ithread] = 0;	_cy_i31[ithread] = 0;
		_cy_r32[ithread] = 0;	_cy_i32[ithread] = 0;
		_cy_r33[ithread] = 0;	_cy_i33[ithread] = 0;
		_cy_r34[ithread] = 0;	_cy_i34[ithread] = 0;
		_cy_r35[ithread] = 0;	_cy_i35[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r00[      0] = -2;
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
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

		khi = n_div_nwt/CY_THREADS;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
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
		_i[0] = 0;		/* Pointer to the BASE and BASEINV arrays. If n divides p, lowest-order digit is always a smallword (_i[0] = 0).	*/

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

		/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
		so for even radix0 values still only need [radix0 >> trailz(radix0)] bjmodn and ii's:
		*/
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii00= 0;
		ii01= (SW_DIV_N*NDIVR/2) % nwt;
		MOD_ADD32(ii01,ii01,nwt,ii02);
		MOD_ADD32(ii02,ii01,nwt,ii03);
		MOD_ADD32(ii03,ii01,nwt,ii04);
		MOD_ADD32(ii04,ii01,nwt,ii05);
		MOD_ADD32(ii05,ii01,nwt,ii06);
		MOD_ADD32(ii06,ii01,nwt,ii07);
		MOD_ADD32(ii07,ii01,nwt,ii08);
		MOD_ADD32(ii08,ii01,nwt,ii09);
		MOD_ADD32(ii09,ii01,nwt,ii10);
		MOD_ADD32(ii10,ii01,nwt,ii11);
		MOD_ADD32(ii11,ii01,nwt,ii12);
		MOD_ADD32(ii12,ii01,nwt,ii13);
		MOD_ADD32(ii13,ii01,nwt,ii14);
		MOD_ADD32(ii14,ii01,nwt,ii15);
		MOD_ADD32(ii15,ii01,nwt,ii16);
		MOD_ADD32(ii16,ii01,nwt,ii17);
		MOD_ADD32(ii17,ii01,nwt,ii18);
		MOD_ADD32(ii18,ii01,nwt,ii19);
		MOD_ADD32(ii19,ii01,nwt,ii20);
		MOD_ADD32(ii20,ii01,nwt,ii21);
		MOD_ADD32(ii21,ii01,nwt,ii22);
		MOD_ADD32(ii22,ii01,nwt,ii23);
		MOD_ADD32(ii23,ii01,nwt,ii24);
		MOD_ADD32(ii24,ii01,nwt,ii25);
		MOD_ADD32(ii25,ii01,nwt,ii26);
		MOD_ADD32(ii26,ii01,nwt,ii27);
		MOD_ADD32(ii27,ii01,nwt,ii28);
		MOD_ADD32(ii28,ii01,nwt,ii29);
		MOD_ADD32(ii29,ii01,nwt,ii30);
		MOD_ADD32(ii30,ii01,nwt,ii31);
		MOD_ADD32(ii31,ii01,nwt,ii32);
		MOD_ADD32(ii32,ii01,nwt,ii33);
		MOD_ADD32(ii33,ii01,nwt,ii34);
		MOD_ADD32(ii34,ii01,nwt,ii35);
	}

	// In non-power-of-2-runlength case, both Mersenne and Fermat-mod share these next 2 loops:
	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
	j = _bjmodnini[CY_THREADS];
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
		MOD_ADD32(_bjmodn23[ithread], j, n, _bjmodn24[ithread]);
		MOD_ADD32(_bjmodn24[ithread], j, n, _bjmodn25[ithread]);
		MOD_ADD32(_bjmodn25[ithread], j, n, _bjmodn26[ithread]);
		MOD_ADD32(_bjmodn26[ithread], j, n, _bjmodn27[ithread]);
		MOD_ADD32(_bjmodn27[ithread], j, n, _bjmodn28[ithread]);
		MOD_ADD32(_bjmodn28[ithread], j, n, _bjmodn29[ithread]);
		MOD_ADD32(_bjmodn29[ithread], j, n, _bjmodn30[ithread]);
		MOD_ADD32(_bjmodn30[ithread], j, n, _bjmodn31[ithread]);
		MOD_ADD32(_bjmodn31[ithread], j, n, _bjmodn32[ithread]);
		MOD_ADD32(_bjmodn32[ithread], j, n, _bjmodn33[ithread]);
		MOD_ADD32(_bjmodn33[ithread], j, n, _bjmodn34[ithread]);
		MOD_ADD32(_bjmodn34[ithread], j, n, _bjmodn35[ithread]);

		// Every (odd_radix)th bjmodn initializer needs to be forced-to-bigword in fermat-mod DWT case:
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
			fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
			*/
			_bjmodn00[ithread] = n;
			_bjmodn09[ithread] = n;
			_bjmodn18[ithread] = n;
			_bjmodn27[ithread] = n;

		}
	}

    /* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
    if(!full_pass)
    {
		khi = 1;
	#if FFT_DEBUG
		fprintf(dbg_file, "radix36_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifdef USE_OMP
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,k1,k2,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,rt,it,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,cy_i28,cy_i29,cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35) default(shared) schedule(static)
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
		ASSERT(HERE, tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->re == crnd && (tmp-1)->im == crnd), "thread-local memcheck failed!");
		ASSERT(HERE, (tmp+10)->re * (tmp+14)->re == 1.0 && (tmp+10)->im * (tmp+14)->im == 1.0, "thread-local memcheck failed!");

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
		tdat[ithread].bjmodn24 = _bjmodn24[ithread];
		tdat[ithread].bjmodn25 = _bjmodn25[ithread];
		tdat[ithread].bjmodn26 = _bjmodn26[ithread];
		tdat[ithread].bjmodn27 = _bjmodn27[ithread];
		tdat[ithread].bjmodn28 = _bjmodn28[ithread];
		tdat[ithread].bjmodn29 = _bjmodn29[ithread];
		tdat[ithread].bjmodn30 = _bjmodn30[ithread];
		tdat[ithread].bjmodn31 = _bjmodn31[ithread];
		tdat[ithread].bjmodn32 = _bjmodn32[ithread];
		tdat[ithread].bjmodn33 = _bjmodn33[ithread];
		tdat[ithread].bjmodn34 = _bjmodn34[ithread];
		tdat[ithread].bjmodn35 = _bjmodn35[ithread];
		/* init carries	*/
		tdat[ithread].cy00 = _cy_r00[ithread];
		tdat[ithread].cy01 = _cy_r01[ithread];
		tdat[ithread].cy02 = _cy_r02[ithread];
		tdat[ithread].cy03 = _cy_r03[ithread];
		tdat[ithread].cy04 = _cy_r04[ithread];
		tdat[ithread].cy05 = _cy_r05[ithread];
		tdat[ithread].cy06 = _cy_r06[ithread];
		tdat[ithread].cy07 = _cy_r07[ithread];
		tdat[ithread].cy08 = _cy_r08[ithread];
		tdat[ithread].cy09 = _cy_r09[ithread];
		tdat[ithread].cy10 = _cy_r10[ithread];
		tdat[ithread].cy11 = _cy_r11[ithread];
		tdat[ithread].cy12 = _cy_r12[ithread];
		tdat[ithread].cy13 = _cy_r13[ithread];
		tdat[ithread].cy14 = _cy_r14[ithread];
		tdat[ithread].cy15 = _cy_r15[ithread];
		tdat[ithread].cy16 = _cy_r16[ithread];
		tdat[ithread].cy17 = _cy_r17[ithread];
		tdat[ithread].cy18 = _cy_r18[ithread];
		tdat[ithread].cy19 = _cy_r19[ithread];
		tdat[ithread].cy20 = _cy_r20[ithread];
		tdat[ithread].cy21 = _cy_r21[ithread];
		tdat[ithread].cy22 = _cy_r22[ithread];
		tdat[ithread].cy23 = _cy_r23[ithread];
		tdat[ithread].cy24 = _cy_r24[ithread];
		tdat[ithread].cy25 = _cy_r25[ithread];
		tdat[ithread].cy26 = _cy_r26[ithread];
		tdat[ithread].cy27 = _cy_r27[ithread];
		tdat[ithread].cy28 = _cy_r28[ithread];
		tdat[ithread].cy29 = _cy_r29[ithread];
		tdat[ithread].cy30 = _cy_r30[ithread];
		tdat[ithread].cy31 = _cy_r31[ithread];
		tdat[ithread].cy32 = _cy_r32[ithread];
		tdat[ithread].cy33 = _cy_r33[ithread];
		tdat[ithread].cy34 = _cy_r34[ithread];
		tdat[ithread].cy35 = _cy_r35[ithread];
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
			*bjmodn24 = _bjmodn24[ithread];
			*bjmodn25 = _bjmodn25[ithread];
			*bjmodn26 = _bjmodn26[ithread];
			*bjmodn27 = _bjmodn27[ithread];
			*bjmodn28 = _bjmodn28[ithread];
			*bjmodn29 = _bjmodn29[ithread];
			*bjmodn30 = _bjmodn30[ithread];
			*bjmodn31 = _bjmodn31[ithread];
			*bjmodn32 = _bjmodn32[ithread];
			*bjmodn33 = _bjmodn33[ithread];
			*bjmodn34 = _bjmodn34[ithread];
			*bjmodn35 = _bjmodn35[ithread];
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
			bjmodn24 = _bjmodn24[ithread];
			bjmodn25 = _bjmodn25[ithread];
			bjmodn26 = _bjmodn26[ithread];
			bjmodn27 = _bjmodn27[ithread];
			bjmodn28 = _bjmodn28[ithread];
			bjmodn29 = _bjmodn29[ithread];
			bjmodn30 = _bjmodn30[ithread];
			bjmodn31 = _bjmodn31[ithread];
			bjmodn32 = _bjmodn32[ithread];
			bjmodn33 = _bjmodn33[ithread];
			bjmodn34 = _bjmodn34[ithread];
			bjmodn35 = _bjmodn35[ithread];
		#endif
		}

	#ifdef USE_SSE2
		/* init carries	*/
		cy_r00->re = _cy_r00[ithread];	cy_i00->re = _cy_i00[ithread];
		cy_r00->im = _cy_r01[ithread];	cy_i00->im = _cy_i01[ithread];
		cy_r02->re = _cy_r02[ithread];	cy_i02->re = _cy_i02[ithread];
		cy_r02->im = _cy_r03[ithread];	cy_i02->im = _cy_i03[ithread];
		cy_r04->re = _cy_r04[ithread];	cy_i04->re = _cy_i04[ithread];
		cy_r04->im = _cy_r05[ithread];	cy_i04->im = _cy_i05[ithread];
		cy_r06->re = _cy_r06[ithread];	cy_i06->re = _cy_i06[ithread];
		cy_r06->im = _cy_r07[ithread];	cy_i06->im = _cy_i07[ithread];
		cy_r08->re = _cy_r08[ithread];	cy_i08->re = _cy_i08[ithread];
		cy_r08->im = _cy_r09[ithread];	cy_i08->im = _cy_i09[ithread];
		cy_r10->re = _cy_r10[ithread];	cy_i10->re = _cy_i10[ithread];
		cy_r10->im = _cy_r11[ithread];	cy_i10->im = _cy_i11[ithread];
		cy_r12->re = _cy_r12[ithread];	cy_i12->re = _cy_i12[ithread];
		cy_r12->im = _cy_r13[ithread];	cy_i12->im = _cy_i13[ithread];
		cy_r14->re = _cy_r14[ithread];	cy_i14->re = _cy_i14[ithread];
		cy_r14->im = _cy_r15[ithread];	cy_i14->im = _cy_i15[ithread];
		cy_r16->re = _cy_r16[ithread];	cy_i16->re = _cy_i16[ithread];
		cy_r16->im = _cy_r17[ithread];	cy_i16->im = _cy_i17[ithread];
		cy_r18->re = _cy_r18[ithread];	cy_i18->re = _cy_i18[ithread];
		cy_r18->im = _cy_r19[ithread];	cy_i18->im = _cy_i19[ithread];
		cy_r20->re = _cy_r20[ithread];	cy_i20->re = _cy_i20[ithread];
		cy_r20->im = _cy_r21[ithread];	cy_i20->im = _cy_i21[ithread];
		cy_r22->re = _cy_r22[ithread];	cy_i22->re = _cy_i22[ithread];
		cy_r22->im = _cy_r23[ithread];	cy_i22->im = _cy_i23[ithread];
		cy_r24->re = _cy_r24[ithread];	cy_i24->re = _cy_i24[ithread];
		cy_r24->im = _cy_r25[ithread];	cy_i24->im = _cy_i25[ithread];
		cy_r26->re = _cy_r26[ithread];	cy_i26->re = _cy_i26[ithread];
		cy_r26->im = _cy_r27[ithread];	cy_i26->im = _cy_i27[ithread];
		cy_r28->re = _cy_r28[ithread];	cy_i28->re = _cy_i28[ithread];
		cy_r28->im = _cy_r29[ithread];	cy_i28->im = _cy_i29[ithread];
		cy_r30->re = _cy_r30[ithread];	cy_i30->re = _cy_i30[ithread];
		cy_r30->im = _cy_r31[ithread];	cy_i30->im = _cy_i31[ithread];
		cy_r32->re = _cy_r32[ithread];	cy_i32->re = _cy_i32[ithread];
		cy_r32->im = _cy_r33[ithread];	cy_i32->im = _cy_i33[ithread];
		cy_r34->re = _cy_r34[ithread];	cy_i34->re = _cy_i34[ithread];
		cy_r34->im = _cy_r35[ithread];	cy_i34->im = _cy_i35[ithread];
	#else
		/* init carries	*/
		cy_r00 = _cy_r00[ithread];	cy_i00 = _cy_i00[ithread];
		cy_r01 = _cy_r01[ithread];	cy_i01 = _cy_i01[ithread];
		cy_r02 = _cy_r02[ithread];	cy_i02 = _cy_i02[ithread];
		cy_r03 = _cy_r03[ithread];	cy_i03 = _cy_i03[ithread];
		cy_r04 = _cy_r04[ithread];	cy_i04 = _cy_i04[ithread];
		cy_r05 = _cy_r05[ithread];	cy_i05 = _cy_i05[ithread];
		cy_r06 = _cy_r06[ithread];	cy_i06 = _cy_i06[ithread];
		cy_r07 = _cy_r07[ithread];	cy_i07 = _cy_i07[ithread];
		cy_r08 = _cy_r08[ithread];	cy_i08 = _cy_i08[ithread];
		cy_r09 = _cy_r09[ithread];	cy_i09 = _cy_i09[ithread];
		cy_r10 = _cy_r10[ithread];	cy_i10 = _cy_i10[ithread];
		cy_r11 = _cy_r11[ithread];	cy_i11 = _cy_i11[ithread];
		cy_r12 = _cy_r12[ithread];	cy_i12 = _cy_i12[ithread];
		cy_r13 = _cy_r13[ithread];	cy_i13 = _cy_i13[ithread];
		cy_r14 = _cy_r14[ithread];	cy_i14 = _cy_i14[ithread];
		cy_r15 = _cy_r15[ithread];	cy_i15 = _cy_i15[ithread];
		cy_r16 = _cy_r16[ithread];	cy_i16 = _cy_i16[ithread];
		cy_r17 = _cy_r17[ithread];	cy_i17 = _cy_i17[ithread];
		cy_r18 = _cy_r18[ithread];	cy_i18 = _cy_i18[ithread];
		cy_r19 = _cy_r19[ithread];	cy_i19 = _cy_i19[ithread];
		cy_r20 = _cy_r20[ithread];	cy_i20 = _cy_i20[ithread];
		cy_r21 = _cy_r21[ithread];	cy_i21 = _cy_i21[ithread];
		cy_r22 = _cy_r22[ithread];	cy_i22 = _cy_i22[ithread];
		cy_r23 = _cy_r23[ithread];	cy_i23 = _cy_i23[ithread];
		cy_r24 = _cy_r24[ithread];	cy_i24 = _cy_i24[ithread];
		cy_r25 = _cy_r25[ithread];	cy_i25 = _cy_i25[ithread];
		cy_r26 = _cy_r26[ithread];	cy_i26 = _cy_i26[ithread];
		cy_r27 = _cy_r27[ithread];	cy_i27 = _cy_i27[ithread];
		cy_r28 = _cy_r28[ithread];	cy_i28 = _cy_i28[ithread];
		cy_r29 = _cy_r29[ithread];	cy_i29 = _cy_i29[ithread];
		cy_r30 = _cy_r30[ithread];	cy_i30 = _cy_i30[ithread];
		cy_r31 = _cy_r31[ithread];	cy_i31 = _cy_i31[ithread];
		cy_r32 = _cy_r32[ithread];	cy_i32 = _cy_i32[ithread];
		cy_r33 = _cy_r33[ithread];	cy_i33 = _cy_i33[ithread];
		cy_r34 = _cy_r34[ithread];	cy_i34 = _cy_i34[ithread];
		cy_r35 = _cy_r35[ithread];	cy_i35 = _cy_i35[ithread];
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
		/*
		!...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do a radix-36 DIT transform...
		*/
#ifdef CTIME
	clock2 = clock();
#endif
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
		#ifdef USE_SSE2

		  #ifdef COMPILER_TYPE_MSVC

			/* Outputs in SSE2 modes are temps 2*9*16 = 18*16 = 0x120 bytes apart: */
		   #if 0
			add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x120)
			add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x120)
			add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x120)
			add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x120)
			add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x120)
			add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x120)
			add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x120)
			add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x120)
			add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0g, 0x120)
		   #else
		  /* 03/15/2008: Decided to try to do all the array-indexing by custom asm, as well, to see if the register-state-saves
		  incurred by mixing asm and hll code represent an appreciable timing penalty. Doing JUST THIS ONE BLOCK this way cut runtime by 3% -
		  so we should avoid such code-mixing, if reasonably possible.
		  */
			add0 = &a[j1    ];
		//	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;
			__asm	mov	eax, add0	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_B */
			__asm	mov	edx, add0
			__asm	mov	esi, p01	/* esi will store power-of-2 multiples of p01 throughout */
			__asm	shl	esi, 3		/* Pointer offset for floating doubles */
			__asm	add edx, esi
			__asm	mov ebx, edx	/* add1 = add0+p01 */
			__asm	add edx, esi
			__asm	mov ecx, edx	/* add3 = add0+p02 */
			__asm	add edx, esi	/* add2 = add0+p03 */
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r00)

		//	add1,0,2,3 = &a[j1+p08]+p0,1,2,3
			__asm	mov	esi, p08
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p08]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x120, 0x240, r02)

		//	add2,3,0,1 = &a[j1+p04]+p0,1,2,3
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p04]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x120, 0x240, r04)

		//	add3,2,1,0 = &a[j1+p28]+p0,1,2,3
			__asm	mov	esi, p24
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p28]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x120, 0x240, r06)

		//	add1,0,2,3 = &a[j1+p24]+p0,1,2,3
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p24]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x120, 0x240, r08)

		//	add0,1,3,2 = &a[j1+p32]+p0,1,2,3
			__asm	mov	esi, p08
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p32]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r0a)

		//	add3,2,1,0 = &a[j1+p12]+p0,1,2,3
			__asm	mov	esi, p20
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p12]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x120, 0x240, r0c)

		//	add2,3,0,1 = &a[j1+p20]+p0,1,2,3
			__asm	mov	esi, p08
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p20]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x120, 0x240, r0e)

		//	add0,1,3,2 = &a[j1+p16]+p0,1,2,3
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p16]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r0g)
		   #endif	// (0)

			/* Radix-9 DFT uses adjacent temps, i.e. stride = 2*16 bytes: */
			SSE2_RADIX_09_DIT_0TWIDDLE(r00,cc1,cc2,cc3m1,cc4,s1p00r,s1p32r,s1p28r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r)
			SSE2_RADIX_09_DIT_0TWIDDLE(r10,cc1,cc2,cc3m1,cc4,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,s1p07r,s1p03r,s1p35r,s1p31r)
			SSE2_RADIX_09_DIT_0TWIDDLE(r20,cc1,cc2,cc3m1,cc4,s1p18r,s1p14r,s1p10r,s1p06r,s1p02r,s1p34r,s1p30r,s1p26r,s1p22r)
			SSE2_RADIX_09_DIT_0TWIDDLE(r30,cc1,cc2,cc3m1,cc4,s1p09r,s1p05r,s1p01r,s1p33r,s1p29r,s1p25r,s1p21r,s1p17r,s1p13r)

		  #else	/* GCC-style inline ASM: */

			add0 = &a[j1    ];
			SSE2_RADIX36_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,r00,s1p00r,cc1);

		  #endif // COMPILER_TYPE_MSVC?

		#else	/* !USE_SSE2 */

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

		/*...and now do 4 radix-9 transforms...*/
						 /*                            inputs                                 */ /*                 outputs                   */
			RADIX_09_DIT(t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,a1p00r,a1p00i,a1p32r,a1p32i,a1p28r,a1p28i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,rt,it,re);
			RADIX_09_DIT(t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p35r,a1p35i,a1p31r,a1p31i,rt,it,re);
			RADIX_09_DIT(t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,a1p18r,a1p18i,a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p34r,a1p34i,a1p30r,a1p30i,a1p26r,a1p26i,a1p22r,a1p22i,rt,it,re);
			RADIX_09_DIT(t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p33r,a1p33i,a1p29r,a1p29i,a1p25r,a1p25i,a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,rt,it,re);

		#endif

#ifdef CTIME
	clock3 = clock();
	dt_fwd += (double)(clock3 - clock2);
	clock2 = clock3;
#endif
	/*...Now do the carries. Since the outputs would
		normally be getting dispatched to 36 separate blocks of the A-array, we need 36 separate carries.	*/

//		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
//		{
#if 1
				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
		#ifdef USE_SSE2

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

		  #if defined(COMPILER_TYPE_MSVC)

			#ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32);
			#else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32);
			#endif

		  #else	/* GCC-style inline ASM: */

			#if ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			#else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			#endif

				/* Bizarre - when I disabled the diagnostic prints above and below, the resulting GCC build immediately gave
					fatal roundoff errors starting on iteration #5 - so insert the bogus [never taken] if() here as a workaround.
					Equally bizarre, inserting the bogus if() *before* the 4 carry-macro calls above gave the correct result as well,
					but ran fully 10% slower. Good old GCC...
				Dec 2011: Suspect this was a side effect of my gcc asm macros not including cc/memory in the clobber list, because
				the code now runs correctly without this hack ... but the code runs sign. faster with iy left in. So still "bizarre" but in a new way.
				*/

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

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

		  #if defined(COMPILER_TYPE_MSVC)

			#if ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32);
			#else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32);
			#endif

		  #else	/* GCC-style inline ASM: */
			#if ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			#else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			#endif
		  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

		#else	/* #ifdef USE_SSE2 */

				/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);
				cmplx_carry_norm_errcheck(a1p28r,a1p28i,cy_r28,bjmodn28,28);
				cmplx_carry_norm_errcheck(a1p29r,a1p29i,cy_r29,bjmodn29,29);
				cmplx_carry_norm_errcheck(a1p30r,a1p30i,cy_r30,bjmodn30,30);
				cmplx_carry_norm_errcheck(a1p31r,a1p31i,cy_r31,bjmodn31,31);
				cmplx_carry_norm_errcheck(a1p32r,a1p32i,cy_r32,bjmodn32,32);
				cmplx_carry_norm_errcheck(a1p33r,a1p33i,cy_r33,bjmodn33,33);
				cmplx_carry_norm_errcheck(a1p34r,a1p34i,cy_r34,bjmodn34,34);
				cmplx_carry_norm_errcheck(a1p35r,a1p35i,cy_r35,bjmodn35,35);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	/* #ifdef USE_SSE2 */
#elif 0
			}
			else
			{
				ASSERT(HERE, 0, "Fermat-mod carries not yet supported for SSE2!");
				fermat_carry_norm_errcheck(a1p00r,a1p00i,cy_r00,cy_i00,ii00,bjmodn00,0 *NDIVR);
				fermat_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,cy_i01,ii01,bjmodn01,1 *NDIVR);
				fermat_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,cy_i02,ii02,bjmodn02,2 *NDIVR);
				fermat_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,cy_i03,ii03,bjmodn03,3 *NDIVR);
				fermat_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,cy_i04,ii04,bjmodn04,4 *NDIVR);
				fermat_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,cy_i05,ii05,bjmodn05,5 *NDIVR);
				fermat_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,cy_i06,ii06,bjmodn06,6 *NDIVR);
				fermat_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,cy_i07,ii07,bjmodn07,7 *NDIVR);
				fermat_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,cy_i08,ii08,bjmodn08,8 *NDIVR);
				fermat_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,cy_i09,ii09,bjmodn09,9 *NDIVR);
				fermat_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,cy_i10,ii10,bjmodn10,10*NDIVR);
				fermat_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,cy_i11,ii11,bjmodn11,11*NDIVR);
				fermat_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,cy_i12,ii12,bjmodn12,12*NDIVR);
				fermat_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,cy_i13,ii13,bjmodn13,13*NDIVR);
				fermat_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,cy_i14,ii14,bjmodn14,14*NDIVR);
				fermat_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,cy_i15,ii15,bjmodn15,15*NDIVR);
				fermat_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,cy_i16,ii16,bjmodn16,16*NDIVR);
				fermat_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,cy_i17,ii17,bjmodn17,17*NDIVR);
				fermat_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,cy_i18,ii18,bjmodn18,18*NDIVR);
				fermat_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,cy_i19,ii19,bjmodn19,19*NDIVR);
				fermat_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,cy_i20,ii20,bjmodn20,20*NDIVR);
				fermat_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,cy_i21,ii21,bjmodn21,21*NDIVR);
				fermat_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,cy_i22,ii22,bjmodn22,22*NDIVR);
				fermat_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,cy_i23,ii23,bjmodn23,23*NDIVR);
				fermat_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,cy_i24,ii24,bjmodn24,24*NDIVR);
				fermat_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,cy_i25,ii25,bjmodn25,25*NDIVR);
				fermat_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,cy_i26,ii26,bjmodn26,26*NDIVR);
				fermat_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,cy_i27,ii27,bjmodn27,27*NDIVR);
				fermat_carry_norm_errcheck(a1p28r,a1p28i,cy_r28,cy_i28,ii28,bjmodn28,28*NDIVR);
				fermat_carry_norm_errcheck(a1p29r,a1p29i,cy_r29,cy_i29,ii29,bjmodn29,29*NDIVR);
				fermat_carry_norm_errcheck(a1p30r,a1p30i,cy_r30,cy_i30,ii30,bjmodn30,30*NDIVR);
				fermat_carry_norm_errcheck(a1p31r,a1p31i,cy_r31,cy_i31,ii31,bjmodn31,31*NDIVR);
				fermat_carry_norm_errcheck(a1p32r,a1p32i,cy_r32,cy_i32,ii32,bjmodn32,32*NDIVR);
				fermat_carry_norm_errcheck(a1p33r,a1p33i,cy_r33,cy_i33,ii33,bjmodn33,33*NDIVR);
				fermat_carry_norm_errcheck(a1p34r,a1p34i,cy_r34,cy_i34,ii34,bjmodn34,34*NDIVR);
				fermat_carry_norm_errcheck(a1p35r,a1p35i,cy_r35,cy_i35,ii35,bjmodn35,35*NDIVR);
			}
#endif

		/*...The radix-36 DIF pass is here:	*/
#ifdef CTIME
	clock3 = clock();
	dt_cy += (double)(clock3 - clock2);
	clock2 = clock3;
#endif

		#ifdef USE_SSE2

		  #ifdef COMPILER_TYPE_MSVC

			/* Radix-9 DFT uses adjacent temps, i.e. stride = 2*16 bytes: */
			SSE2_RADIX_09_DIF_0TWIDDLE(r00,cc1,cc2,cc3m1,cc4,s1p00r,s1p32r,s1p28r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r)
			SSE2_RADIX_09_DIF_0TWIDDLE(r10,cc1,cc2,cc3m1,cc4,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,s1p07r,s1p03r,s1p35r,s1p31r)
			SSE2_RADIX_09_DIF_0TWIDDLE(r20,cc1,cc2,cc3m1,cc4,s1p18r,s1p14r,s1p10r,s1p06r,s1p02r,s1p34r,s1p30r,s1p26r,s1p22r)
			SSE2_RADIX_09_DIF_0TWIDDLE(r30,cc1,cc2,cc3m1,cc4,s1p09r,s1p05r,s1p01r,s1p33r,s1p29r,s1p25r,s1p21r,s1p17r,s1p13r)

			/* Outputs in SSE2 modes are temps 2*9*16 = 18*16 = 0x120 bytes apart: */
		   #if 0
			add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x120)
			add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x120)
			add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x120)
			add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x120)
			add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x120)
			add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x120)
			add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x120)
			add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x120)
			add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0g, 0x120)
		   #else
		  /* In the DIF case, the _B version of the SSE2_RADIX4_DIF_0TWIDDLE_STRIDE macro similarly leaves eax [in fact all 4 of e*x] unchanged,
		  but if using that version of the macro, need to restore the pointer offset in esi after each macro call. We can do this by subtracting
		  one of e*x from one of ebx,ecx,edx, in order to leave eax untouched.
		  */

		//	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;
			__asm	mov	eax, add0	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */
			__asm	mov	edx, add0
			__asm	mov	esi, p01	/* esi will store power-of-2 multiples of p01 throughout */
			__asm	shl	esi, 3		/* Pointer offset for floating doubles */
			__asm	add edx, esi
			__asm	mov ebx, edx	/* add0+p01 */
			__asm	add edx, esi
			__asm	mov ecx, edx	/* add0+p02 */
			__asm	add edx, esi	/* add0+p03 */
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x120, 0x240, eax,ebx,edx,ecx)

		//	add0,1,3,2 = &a[j1+p32]+p0,1,2,3
			__asm	mov	esi, p32
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p32]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, 0x120, 0x240, eax,ebx,edx,ecx)

		//	add2,3,0,1 = &a[j1+p20]+p0,1,2,3
			__asm	mov	esi, p12
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p20]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, 0x120, 0x240, ecx,edx,eax,ebx)

		//	add1,0,2,3 = &a[j1+p08]+p0,1,2,3
			__asm	mov	esi, p12
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p08]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, 0x120, 0x240, ebx,eax,ecx,edx)

		//	add3,2,1,0 = &a[j1+p28]+p0,1,2,3
			__asm	mov	esi, p20
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p28]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, 0x120, 0x240, edx,ecx,ebx,eax)

		//	add0,1,3,2 = &a[j1+p16]+p0,1,2,3
			__asm	mov	esi, p12
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p16]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0a, 0x120, 0x240, eax,ebx,edx,ecx)

		//	add2,3,0,1 = &a[j1+p04]+p0,1,2,3
			__asm	mov	esi, p12
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p04]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0c, 0x120, 0x240, ecx,edx,eax,ebx)

		//	add1,0,2,3 = &a[j1+p24]+p0,1,2,3
			__asm	mov	esi, p20
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p24]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0e, 0x120, 0x240, ebx,eax,ecx,edx)

		//	add,2,1,03 = &a[j1+p12]+p0,1,2,3
			__asm	mov	esi, p12
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p12]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0g, 0x120, 0x240, edx,ecx,ebx,eax)

		   #endif	// #if(0)

		  #else	/* GCC-style inline ASM: */

			add0 = &a[j1    ];
			SSE2_RADIX36_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,r00,s1p00r,cc1);

		  #endif

		#else	/* !USE_SSE2 */

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

		#endif

#ifdef CTIME
	clock3 = clock();
	dt_inv += (double)(clock3 - clock2);
	clock2 = clock3;
#endif
			}

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
	#ifdef USE_SSE2
		_cy_r00[ithread] = cy_r00->re;	_cy_i00[ithread] = cy_i00->re;
		_cy_r01[ithread] = cy_r00->im;	_cy_i01[ithread] = cy_i00->im;
		_cy_r02[ithread] = cy_r02->re;	_cy_i02[ithread] = cy_i02->re;
		_cy_r03[ithread] = cy_r02->im;	_cy_i03[ithread] = cy_i02->im;
		_cy_r04[ithread] = cy_r04->re;	_cy_i04[ithread] = cy_i04->re;
		_cy_r05[ithread] = cy_r04->im;	_cy_i05[ithread] = cy_i04->im;
		_cy_r06[ithread] = cy_r06->re;	_cy_i06[ithread] = cy_i06->re;
		_cy_r07[ithread] = cy_r06->im;	_cy_i07[ithread] = cy_i06->im;
		_cy_r08[ithread] = cy_r08->re;	_cy_i08[ithread] = cy_i08->re;
		_cy_r09[ithread] = cy_r08->im;	_cy_i09[ithread] = cy_i08->im;
		_cy_r10[ithread] = cy_r10->re;	_cy_i10[ithread] = cy_i10->re;
		_cy_r11[ithread] = cy_r10->im;	_cy_i11[ithread] = cy_i10->im;
		_cy_r12[ithread] = cy_r12->re;	_cy_i12[ithread] = cy_i12->re;
		_cy_r13[ithread] = cy_r12->im;	_cy_i13[ithread] = cy_i12->im;
		_cy_r14[ithread] = cy_r14->re;	_cy_i14[ithread] = cy_i14->re;
		_cy_r15[ithread] = cy_r14->im;	_cy_i15[ithread] = cy_i14->im;
		_cy_r16[ithread] = cy_r16->re;	_cy_i16[ithread] = cy_i16->re;
		_cy_r17[ithread] = cy_r16->im;	_cy_i17[ithread] = cy_i16->im;
		_cy_r18[ithread] = cy_r18->re;	_cy_i18[ithread] = cy_i18->re;
		_cy_r19[ithread] = cy_r18->im;	_cy_i19[ithread] = cy_i18->im;
		_cy_r20[ithread] = cy_r20->re;	_cy_i20[ithread] = cy_i20->re;
		_cy_r21[ithread] = cy_r20->im;	_cy_i21[ithread] = cy_i20->im;
		_cy_r22[ithread] = cy_r22->re;	_cy_i22[ithread] = cy_i22->re;
		_cy_r23[ithread] = cy_r22->im;	_cy_i23[ithread] = cy_i22->im;
		_cy_r24[ithread] = cy_r24->re;	_cy_i24[ithread] = cy_i24->re;
		_cy_r25[ithread] = cy_r24->im;	_cy_i25[ithread] = cy_i24->im;
		_cy_r26[ithread] = cy_r26->re;	_cy_i26[ithread] = cy_i26->re;
		_cy_r27[ithread] = cy_r26->im;	_cy_i27[ithread] = cy_i26->im;
		_cy_r28[ithread] = cy_r28->re;	_cy_i28[ithread] = cy_i28->re;
		_cy_r29[ithread] = cy_r28->im;	_cy_i29[ithread] = cy_i28->im;
		_cy_r30[ithread] = cy_r30->re;	_cy_i30[ithread] = cy_i30->re;
		_cy_r31[ithread] = cy_r30->im;	_cy_i31[ithread] = cy_i30->im;
		_cy_r32[ithread] = cy_r32->re;	_cy_i32[ithread] = cy_i32->re;
		_cy_r33[ithread] = cy_r32->im;	_cy_i33[ithread] = cy_i32->im;
		_cy_r34[ithread] = cy_r34->re;	_cy_i34[ithread] = cy_i34->re;
		_cy_r35[ithread] = cy_r34->im;	_cy_i35[ithread] = cy_i34->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
	#else
		_cy_r00[ithread] = cy_r00;	_cy_i00[ithread] = cy_i00;
		_cy_r01[ithread] = cy_r01;	_cy_i01[ithread] = cy_i01;
		_cy_r02[ithread] = cy_r02;	_cy_i02[ithread] = cy_i02;
		_cy_r03[ithread] = cy_r03;	_cy_i03[ithread] = cy_i03;
		_cy_r04[ithread] = cy_r04;	_cy_i04[ithread] = cy_i04;
		_cy_r05[ithread] = cy_r05;	_cy_i05[ithread] = cy_i05;
		_cy_r06[ithread] = cy_r06;	_cy_i06[ithread] = cy_i06;
		_cy_r07[ithread] = cy_r07;	_cy_i07[ithread] = cy_i07;
		_cy_r08[ithread] = cy_r08;	_cy_i08[ithread] = cy_i08;
		_cy_r09[ithread] = cy_r09;	_cy_i09[ithread] = cy_i09;
		_cy_r10[ithread] = cy_r10;	_cy_i10[ithread] = cy_i10;
		_cy_r11[ithread] = cy_r11;	_cy_i11[ithread] = cy_i11;
		_cy_r12[ithread] = cy_r12;	_cy_i12[ithread] = cy_i12;
		_cy_r13[ithread] = cy_r13;	_cy_i13[ithread] = cy_i13;
		_cy_r14[ithread] = cy_r14;	_cy_i14[ithread] = cy_i14;
		_cy_r15[ithread] = cy_r15;	_cy_i15[ithread] = cy_i15;
		_cy_r16[ithread] = cy_r16;	_cy_i16[ithread] = cy_i16;
		_cy_r17[ithread] = cy_r17;	_cy_i17[ithread] = cy_i17;
		_cy_r18[ithread] = cy_r18;	_cy_i18[ithread] = cy_i18;
		_cy_r19[ithread] = cy_r19;	_cy_i19[ithread] = cy_i19;
		_cy_r20[ithread] = cy_r20;	_cy_i20[ithread] = cy_i20;
		_cy_r21[ithread] = cy_r21;	_cy_i21[ithread] = cy_i21;
		_cy_r22[ithread] = cy_r22;	_cy_i22[ithread] = cy_i22;
		_cy_r23[ithread] = cy_r23;	_cy_i23[ithread] = cy_i23;
		_cy_r24[ithread] = cy_r24;	_cy_i24[ithread] = cy_i24;
		_cy_r25[ithread] = cy_r25;	_cy_i25[ithread] = cy_i25;
		_cy_r26[ithread] = cy_r26;	_cy_i26[ithread] = cy_i26;
		_cy_r27[ithread] = cy_r27;	_cy_i27[ithread] = cy_i27;
		_cy_r28[ithread] = cy_r28;	_cy_i28[ithread] = cy_i28;
		_cy_r29[ithread] = cy_r29;	_cy_i29[ithread] = cy_i29;
		_cy_r30[ithread] = cy_r30;	_cy_i30[ithread] = cy_i30;
		_cy_r31[ithread] = cy_r31;	_cy_i31[ithread] = cy_i31;
		_cy_r32[ithread] = cy_r32;	_cy_i32[ithread] = cy_i32;
		_cy_r33[ithread] = cy_r33;	_cy_i33[ithread] = cy_i33;
		_cy_r34[ithread] = cy_r34;	_cy_i34[ithread] = cy_i34;
		_cy_r35[ithread] = cy_r35;	_cy_i35[ithread] = cy_i35;
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
		ASSERT(HERE, 0x0 == cy36_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that

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

		_cy_r00[ithread] = tdat[ithread].cy00;
		_cy_r01[ithread] = tdat[ithread].cy01;
		_cy_r02[ithread] = tdat[ithread].cy02;
		_cy_r03[ithread] = tdat[ithread].cy03;
		_cy_r04[ithread] = tdat[ithread].cy04;
		_cy_r05[ithread] = tdat[ithread].cy05;
		_cy_r06[ithread] = tdat[ithread].cy06;
		_cy_r07[ithread] = tdat[ithread].cy07;
		_cy_r08[ithread] = tdat[ithread].cy08;
		_cy_r09[ithread] = tdat[ithread].cy09;
		_cy_r10[ithread] = tdat[ithread].cy10;
		_cy_r11[ithread] = tdat[ithread].cy11;
		_cy_r12[ithread] = tdat[ithread].cy12;
		_cy_r13[ithread] = tdat[ithread].cy13;
		_cy_r14[ithread] = tdat[ithread].cy14;
		_cy_r15[ithread] = tdat[ithread].cy15;
		_cy_r16[ithread] = tdat[ithread].cy16;
		_cy_r17[ithread] = tdat[ithread].cy17;
		_cy_r18[ithread] = tdat[ithread].cy18;
		_cy_r19[ithread] = tdat[ithread].cy19;
		_cy_r20[ithread] = tdat[ithread].cy20;
		_cy_r21[ithread] = tdat[ithread].cy21;
		_cy_r22[ithread] = tdat[ithread].cy22;
		_cy_r23[ithread] = tdat[ithread].cy23;
		_cy_r24[ithread] = tdat[ithread].cy24;
		_cy_r25[ithread] = tdat[ithread].cy25;
		_cy_r26[ithread] = tdat[ithread].cy26;
		_cy_r27[ithread] = tdat[ithread].cy27;
		_cy_r28[ithread] = tdat[ithread].cy28;
		_cy_r29[ithread] = tdat[ithread].cy29;
		_cy_r30[ithread] = tdat[ithread].cy30;
		_cy_r31[ithread] = tdat[ithread].cy31;
		_cy_r32[ithread] = tdat[ithread].cy32;
		_cy_r33[ithread] = tdat[ithread].cy33;
		_cy_r34[ithread] = tdat[ithread].cy34;
		_cy_r35[ithread] = tdat[ithread].cy35;
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-36 forward DIF FFT of the first block of 36 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 36 outputs of (1);
	!   (3) Reweight and perform a radix-36 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 36 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t00= _cy_r00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];
		t0a= _cy_r05[CY_THREADS - 1];
		t0c= _cy_r06[CY_THREADS - 1];
		t0e= _cy_r07[CY_THREADS - 1];
		t0g= _cy_r08[CY_THREADS - 1];
		t10= _cy_r09[CY_THREADS - 1];
		t12= _cy_r10[CY_THREADS - 1];
		t14= _cy_r11[CY_THREADS - 1];
		t16= _cy_r12[CY_THREADS - 1];
		t18= _cy_r13[CY_THREADS - 1];
		t1a= _cy_r14[CY_THREADS - 1];
		t1c= _cy_r15[CY_THREADS - 1];
		t1e= _cy_r16[CY_THREADS - 1];
		t1g= _cy_r17[CY_THREADS - 1];
		t20= _cy_r18[CY_THREADS - 1];
		t22= _cy_r19[CY_THREADS - 1];
		t24= _cy_r20[CY_THREADS - 1];
		t26= _cy_r21[CY_THREADS - 1];
		t28= _cy_r22[CY_THREADS - 1];
		t2a= _cy_r23[CY_THREADS - 1];
		t2c= _cy_r24[CY_THREADS - 1];
		t2e= _cy_r25[CY_THREADS - 1];
		t2g= _cy_r26[CY_THREADS - 1];
		t30= _cy_r27[CY_THREADS - 1];
		t32= _cy_r28[CY_THREADS - 1];
		t34= _cy_r29[CY_THREADS - 1];
		t36= _cy_r30[CY_THREADS - 1];
		t38= _cy_r31[CY_THREADS - 1];
		t3a= _cy_r32[CY_THREADS - 1];
		t3c= _cy_r33[CY_THREADS - 1];
		t3e= _cy_r34[CY_THREADS - 1];
		t3g= _cy_r35[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"radix36_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];
			_cy_r28[ithread] = _cy_r28[ithread-1];
			_cy_r29[ithread] = _cy_r29[ithread-1];
			_cy_r30[ithread] = _cy_r30[ithread-1];
			_cy_r31[ithread] = _cy_r31[ithread-1];
			_cy_r32[ithread] = _cy_r32[ithread-1];
			_cy_r33[ithread] = _cy_r33[ithread-1];
			_cy_r34[ithread] = _cy_r34[ithread-1];
			_cy_r35[ithread] = _cy_r35[ithread-1];
		}

		_cy_r00[0] =+t3g;	/* ...The wraparound carry is here: */
		_cy_r01[0] = t00;
		_cy_r02[0] = t02;
		_cy_r03[0] = t04;
		_cy_r04[0] = t06;
		_cy_r05[0] = t08;
		_cy_r06[0] = t0a;
		_cy_r07[0] = t0c;
		_cy_r08[0] = t0e;
		_cy_r09[0] = t0g;
		_cy_r10[0] = t10;
		_cy_r11[0] = t12;
		_cy_r12[0] = t14;
		_cy_r13[0] = t16;
		_cy_r14[0] = t18;
		_cy_r15[0] = t1a;
		_cy_r16[0] = t1c;
		_cy_r17[0] = t1e;
		_cy_r18[0] = t1g;
		_cy_r19[0] = t20;
		_cy_r20[0] = t22;
		_cy_r21[0] = t24;
		_cy_r22[0] = t26;
		_cy_r23[0] = t28;
		_cy_r24[0] = t2a;
		_cy_r25[0] = t2c;
		_cy_r26[0] = t2e;
		_cy_r27[0] = t2g;
		_cy_r28[0] = t30;
		_cy_r29[0] = t32;
		_cy_r30[0] = t34;
		_cy_r31[0] = t36;
		_cy_r32[0] = t38;
		_cy_r33[0] = t3a;
		_cy_r34[0] = t3c;
		_cy_r35[0] = t3e;
	}
	else
	{
		t00= _cy_r00[CY_THREADS - 1];	t01= _cy_i00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];	t03= _cy_i01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];	t05= _cy_i02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];	t07= _cy_i03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];	t09= _cy_i04[CY_THREADS - 1];
		t0a= _cy_r05[CY_THREADS - 1];	t0b= _cy_i05[CY_THREADS - 1];
		t0c= _cy_r06[CY_THREADS - 1];	t0d= _cy_i06[CY_THREADS - 1];
		t0e= _cy_r07[CY_THREADS - 1];	t0f= _cy_i07[CY_THREADS - 1];
		t0g= _cy_r08[CY_THREADS - 1];	t0h= _cy_i08[CY_THREADS - 1];
		t10= _cy_r09[CY_THREADS - 1];	t11= _cy_i09[CY_THREADS - 1];
		t12= _cy_r10[CY_THREADS - 1];	t13= _cy_i10[CY_THREADS - 1];
		t14= _cy_r11[CY_THREADS - 1];	t15= _cy_i11[CY_THREADS - 1];
		t16= _cy_r12[CY_THREADS - 1];	t17= _cy_i12[CY_THREADS - 1];
		t18= _cy_r13[CY_THREADS - 1];	t19= _cy_i13[CY_THREADS - 1];
		t1a= _cy_r14[CY_THREADS - 1];	t1b= _cy_i14[CY_THREADS - 1];
		t1c= _cy_r15[CY_THREADS - 1];	t1d= _cy_i15[CY_THREADS - 1];
		t1e= _cy_r16[CY_THREADS - 1];	t1f= _cy_i16[CY_THREADS - 1];
		t1g= _cy_r17[CY_THREADS - 1];	t1h= _cy_i17[CY_THREADS - 1];
		t20= _cy_r18[CY_THREADS - 1];	t21= _cy_i18[CY_THREADS - 1];
		t22= _cy_r19[CY_THREADS - 1];	t23= _cy_i19[CY_THREADS - 1];
		t24= _cy_r20[CY_THREADS - 1];	t25= _cy_i20[CY_THREADS - 1];
		t26= _cy_r21[CY_THREADS - 1];	t27= _cy_i21[CY_THREADS - 1];
		t28= _cy_r22[CY_THREADS - 1];	t29= _cy_i22[CY_THREADS - 1];
		t2a= _cy_r23[CY_THREADS - 1];	t2b= _cy_i23[CY_THREADS - 1];
		t2c= _cy_r24[CY_THREADS - 1];	t2d= _cy_i24[CY_THREADS - 1];
		t2e= _cy_r25[CY_THREADS - 1];	t2f= _cy_i25[CY_THREADS - 1];
		t2g= _cy_r26[CY_THREADS - 1];	t2h= _cy_i26[CY_THREADS - 1];
		t30= _cy_r27[CY_THREADS - 1];	t31= _cy_i27[CY_THREADS - 1];
		t32= _cy_r28[CY_THREADS - 1];	t33= _cy_i28[CY_THREADS - 1];
		t34= _cy_r29[CY_THREADS - 1];	t35= _cy_i29[CY_THREADS - 1];
		t36= _cy_r30[CY_THREADS - 1];	t37= _cy_i30[CY_THREADS - 1];
		t38= _cy_r31[CY_THREADS - 1];	t39= _cy_i31[CY_THREADS - 1];
		t3a= _cy_r32[CY_THREADS - 1];	t3b= _cy_i32[CY_THREADS - 1];
		t3c= _cy_r33[CY_THREADS - 1];	t3d= _cy_i33[CY_THREADS - 1];
		t3e= _cy_r34[CY_THREADS - 1];	t3f= _cy_i34[CY_THREADS - 1];
		t3g= _cy_r35[CY_THREADS - 1];	t3h= _cy_i35[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"radix36_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];		_cy_i00[ithread] = _cy_i00[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];		_cy_i01[ithread] = _cy_i01[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];		_cy_i02[ithread] = _cy_i02[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];		_cy_i03[ithread] = _cy_i03[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];		_cy_i04[ithread] = _cy_i04[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];		_cy_i05[ithread] = _cy_i05[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];		_cy_i06[ithread] = _cy_i06[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];		_cy_i07[ithread] = _cy_i07[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];		_cy_i08[ithread] = _cy_i08[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];		_cy_i09[ithread] = _cy_i09[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];		_cy_i10[ithread] = _cy_i10[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];		_cy_i11[ithread] = _cy_i11[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];		_cy_i12[ithread] = _cy_i12[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];		_cy_i13[ithread] = _cy_i13[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];		_cy_i14[ithread] = _cy_i14[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];		_cy_i15[ithread] = _cy_i15[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];		_cy_i16[ithread] = _cy_i16[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];		_cy_i17[ithread] = _cy_i17[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];		_cy_i18[ithread] = _cy_i18[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];		_cy_i19[ithread] = _cy_i19[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];		_cy_i20[ithread] = _cy_i20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];		_cy_i21[ithread] = _cy_i21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];		_cy_i22[ithread] = _cy_i22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];		_cy_i23[ithread] = _cy_i23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];		_cy_i24[ithread] = _cy_i24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];		_cy_i25[ithread] = _cy_i25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];		_cy_i26[ithread] = _cy_i26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];		_cy_i27[ithread] = _cy_i27[ithread-1];
			_cy_r28[ithread] = _cy_r28[ithread-1];		_cy_i28[ithread] = _cy_i28[ithread-1];
			_cy_r29[ithread] = _cy_r29[ithread-1];		_cy_i29[ithread] = _cy_i29[ithread-1];
			_cy_r30[ithread] = _cy_r30[ithread-1];		_cy_i30[ithread] = _cy_i30[ithread-1];
			_cy_r31[ithread] = _cy_r31[ithread-1];		_cy_i31[ithread] = _cy_i31[ithread-1];
			_cy_r32[ithread] = _cy_r32[ithread-1];		_cy_i32[ithread] = _cy_i32[ithread-1];
			_cy_r33[ithread] = _cy_r33[ithread-1];		_cy_i33[ithread] = _cy_i33[ithread-1];
			_cy_r34[ithread] = _cy_r34[ithread-1];		_cy_i34[ithread] = _cy_i34[ithread-1];
			_cy_r35[ithread] = _cy_r35[ithread-1];		_cy_i35[ithread] = _cy_i35[ithread-1];
		}

		_cy_r00[0] =-t3g;	_cy_i00[0] =+t3h;	/* ...The 2 Mo"bius carries are here: */
		_cy_r01[0] = t00;	_cy_i01[0] = t01;
		_cy_r02[0] = t02;	_cy_i02[0] = t03;
		_cy_r03[0] = t04;	_cy_i03[0] = t05;
		_cy_r04[0] = t06;	_cy_i04[0] = t07;
		_cy_r05[0] = t08;	_cy_i05[0] = t09;
		_cy_r06[0] = t0a;	_cy_i06[0] = t0b;
		_cy_r07[0] = t0c;	_cy_i07[0] = t0d;
		_cy_r08[0] = t0e;	_cy_i08[0] = t0f;
		_cy_r09[0] = t0g;	_cy_i09[0] = t0h;
		_cy_r10[0] = t10;	_cy_i10[0] = t11;
		_cy_r11[0] = t12;	_cy_i11[0] = t13;
		_cy_r12[0] = t14;	_cy_i12[0] = t15;
		_cy_r13[0] = t16;	_cy_i13[0] = t17;
		_cy_r14[0] = t18;	_cy_i14[0] = t19;
		_cy_r15[0] = t1a;	_cy_i15[0] = t1b;
		_cy_r16[0] = t1c;	_cy_i16[0] = t1d;
		_cy_r17[0] = t1e;	_cy_i17[0] = t1f;
		_cy_r18[0] = t1g;	_cy_i18[0] = t1h;
		_cy_r19[0] = t20;	_cy_i19[0] = t21;
		_cy_r20[0] = t22;	_cy_i20[0] = t23;
		_cy_r21[0] = t24;	_cy_i21[0] = t25;
		_cy_r22[0] = t26;	_cy_i22[0] = t27;
		_cy_r23[0] = t28;	_cy_i23[0] = t29;
		_cy_r24[0] = t2a;	_cy_i24[0] = t2b;
		_cy_r25[0] = t2c;	_cy_i25[0] = t2d;
		_cy_r26[0] = t2e;	_cy_i26[0] = t2f;
		_cy_r27[0] = t2g;	_cy_i27[0] = t2h;
		_cy_r28[0] = t30;	_cy_i28[0] = t31;
		_cy_r29[0] = t32;	_cy_i29[0] = t33;
		_cy_r30[0] = t34;	_cy_i30[0] = t35;
		_cy_r31[0] = t36;	_cy_i31[0] = t37;
		_cy_r32[0] = t38;	_cy_i32[0] = t39;
		_cy_r33[0] = t3a;	_cy_i33[0] = t3b;
		_cy_r34[0] = t3c;	_cy_i34[0] = t3d;
		_cy_r35[0] = t3e;	_cy_i35[0] = t3f;
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
}	/* endfor(outer) */

#ifdef CTIME
	clock2 = clock();
	dt_tot = (double)(clock2 - clock1);
	printf("radix36_carry cycle times: total = %10.5f, fwd = %10.5f, inv = %10.5f, cy = %10.5f\n", dt_tot*ICPS, dt_fwd*ICPS, dt_inv*ICPS, dt_cy*ICPS);
#endif

    t00 = 0;
    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		t00 += fabs(_cy_r00[0])+fabs(_cy_r01[0])+fabs(_cy_r02[0])+fabs(_cy_r03[0])+fabs(_cy_r04[0])+fabs(_cy_r05[0])+fabs(_cy_r06[0])+fabs(_cy_r07[0])+fabs(_cy_r08[0])+fabs(_cy_r09[0])+fabs(_cy_r10[0])+fabs(_cy_r11[0])+fabs(_cy_r12[0])+fabs(_cy_r13[0])+fabs(_cy_r14[0])+fabs(_cy_r15[0])+fabs(_cy_r16[0])+fabs(_cy_r17[0])+fabs(_cy_r18[0])+fabs(_cy_r19[0])+fabs(_cy_r20[0])+fabs(_cy_r21[0])+fabs(_cy_r22[0])+fabs(_cy_r23[0])+fabs(_cy_r24[0])+fabs(_cy_r25[0])+fabs(_cy_r26[0])+fabs(_cy_r27[0])+fabs(_cy_r28[0])+fabs(_cy_r29[0])+fabs(_cy_r30[0])+fabs(_cy_r31[0])+fabs(_cy_r32[0])+fabs(_cy_r33[0])+fabs(_cy_r34[0])+fabs(_cy_r35[0]);
		t00 += fabs(_cy_i00[0])+fabs(_cy_i01[0])+fabs(_cy_i02[0])+fabs(_cy_i03[0])+fabs(_cy_i04[0])+fabs(_cy_i05[0])+fabs(_cy_i06[0])+fabs(_cy_i07[0])+fabs(_cy_i08[0])+fabs(_cy_i09[0])+fabs(_cy_i10[0])+fabs(_cy_i11[0])+fabs(_cy_i12[0])+fabs(_cy_i13[0])+fabs(_cy_i14[0])+fabs(_cy_i15[0])+fabs(_cy_i16[0])+fabs(_cy_i17[0])+fabs(_cy_i18[0])+fabs(_cy_i19[0])+fabs(_cy_i20[0])+fabs(_cy_i21[0])+fabs(_cy_i22[0])+fabs(_cy_i23[0])+fabs(_cy_i24[0])+fabs(_cy_i25[0])+fabs(_cy_i26[0])+fabs(_cy_i27[0])+fabs(_cy_i28[0])+fabs(_cy_i29[0])+fabs(_cy_i30[0])+fabs(_cy_i31[0])+fabs(_cy_i32[0])+fabs(_cy_i33[0])+fabs(_cy_i34[0])+fabs(_cy_i35[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
    }

	if(t00 != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix36_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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
	double rt,it,re
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h;

	if(!first_entry && (n/36) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/36;

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
	}

/*...The radix-36 pass is here.	*/

    for(j=0; j < NDIVR; j += 2)
    {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		/*...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do 4 radix-9 transforms...*/
	/*
	Twiddleless version arranges 4 sets of radix-9 DFT inputs as follows: 0 in upper left corner, decrement 4 horizontally and 9 vertically:

		RADIX_09_DFT(00,32,28,24,20,16,12,08,04)
		RADIX_09_DFT(27,23,19,15,11,07,03,35,31)
		RADIX_09_DFT(18,14,10,06,02,34,30,26,22)
		RADIX_09_DFT(09,05,01,33,29,25,21,17,13)

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-4 DFT outputs.
	*/
		RADIX_09_DIF(a[j1    ],a[j2    ],a[j1+p32],a[j2+p32],a[j1+p28],a[j2+p28],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,rt,it,re);	jt = j1+p03; jp = j2+p03;
		RADIX_09_DIF(a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,rt,it,re);	jt = j1+p02; jp = j2+p02;
		RADIX_09_DIF(a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,rt,it,re);	jt = j1+p01; jp = j2+p01;
		RADIX_09_DIF(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,rt,it,re);

		/*...and now do 9 radix-4 transforms...*/
					/*          inputs           */ /*                                      outputs                                      */
		RADIX_04_DIF(t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it);	jt = j1+p32; jp = j2+p32;
		RADIX_04_DIF(t02,t03,t12,t13,t22,t23,t32,t33,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_04_DIF(t04,t05,t14,t15,t24,t25,t34,t35,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_04_DIF(t06,t07,t16,t17,t26,t27,t36,t37,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_04_DIF(t08,t09,t18,t19,t28,t29,t38,t39,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_04_DIF(t0a,t0b,t1a,t1b,t2a,t2b,t3a,t3b,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_04_DIF(t0c,t0d,t1c,t1d,t2c,t2d,t3c,t3d,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p24; jp = j2+p24;
		RADIX_04_DIF(t0e,t0f,t1e,t1f,t2e,t2f,t3e,t3f,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p12; jp = j2+p12;
		RADIX_04_DIF(t0g,t0h,t1g,t1h,t2g,t2h,t3g,t3h,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

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
	double rt,it,re
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h;

	if(!first_entry && (n/36) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/36;

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
	}

/*...The radix-36 pass is here.	*/

      for(j=0; j < NDIVR; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
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
						 /*          inputs           */ /*                                      outputs                                      */
		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00,t01,t10,t11,t20,t21,t30,t31,rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t02,t03,t12,t13,t22,t23,t32,t33,rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t04,t05,t14,t15,t24,t25,t34,t35,rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t06,t07,t16,t17,t26,t27,t36,t37,rt,it);	jt = j1+p24; jp = j2+p24;
		RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t08,t09,t18,t19,t28,t29,t38,t39,rt,it);	jt = j1+p32; jp = j2+p32;
		RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0a,t0b,t1a,t1b,t2a,t2b,t3a,t3b,rt,it);	jt = j1+p12; jp = j2+p12;
		RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t0c,t0d,t1c,t1d,t2c,t2d,t3c,t3d,rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0e,t0f,t1e,t1f,t2e,t2f,t3e,t3f,rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0g,t0h,t1g,t1h,t2g,t2h,t3g,t3h,rt,it);

		/*...and now do 4 radix-9 transforms...*/
					 /*                            inputs                                 */ /*                 outputs                   */
		RADIX_09_DIT(t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,a[j1    ],a[j2    ],a[j1+p32],a[j2+p32],a[j1+p28],a[j2+p28],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04],rt,it,re);	jt = j1+p03; jp = j2+p03;
		RADIX_09_DIT(t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],rt,it,re);	jt = j1+p02; jp = j2+p02;
		RADIX_09_DIT(t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],rt,it,re);	jt = j1+p01; jp = j2+p01;
		RADIX_09_DIT(t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],rt,it,re);

		/* Totals: 4*radix09 + 9*radix04 = 4*(68 FADD, 40 FMUL)	+ 9*(16 FADD, 0 FMUL) = 416 FADD, 160 FMUL	*/
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
	cy36_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 36;
		const double crnd = 3.0*0x4000000*0x2000000;
		int j,j1,j2,k;
		int l,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32;
		double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09
			,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19
			,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29
			,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35;
		struct complex *cc1, *max_err, *sse2_rnd, *half_arr, *tmp, *r00
			,*s1p00r,*s1p04r,*s1p08r,*s1p12r,*s1p16r,*s1p20r,*s1p24r,*s1p28r,*s1p32r;
		struct complex *cy00,*cy02,*cy04,*cy06,*cy08,*cy10,*cy12,*cy14,*cy16,*cy18,*cy20,*cy22,*cy24,*cy26,*cy28,*cy30,*cy32,*cy34;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
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

		r00	= thread_arg->r00;
		tmp	= r00 + 0x48;
		s1p00r = tmp + 0x00;
		s1p04r = tmp + 0x08;
		s1p08r = tmp + 0x10;
		s1p12r = tmp + 0x18;
		s1p16r = tmp + 0x20;
		s1p20r = tmp + 0x28;
		s1p24r = tmp + 0x30;
		s1p28r = tmp + 0x38;
		s1p32r = tmp + 0x40;
		cc1		= tmp + 0x48;
		cy00	= tmp + 0x50;
		cy02	= tmp + 0x51;
		cy04	= tmp + 0x52;
		cy06	= tmp + 0x53;
		cy08	= tmp + 0x54;
		cy10	= tmp + 0x55;
		cy12	= tmp + 0x56;
		cy14	= tmp + 0x57;
		cy16	= tmp + 0x58;
		cy18	= tmp + 0x59;
		cy20	= tmp + 0x5a;
		cy22	= tmp + 0x5b;
		cy24	= tmp + 0x5c;
		cy26	= tmp + 0x5d;
		cy28	= tmp + 0x5e;
		cy30	= tmp + 0x5f;
		cy32	= tmp + 0x60;
		cy34	= tmp + 0x61;
	/* For future Fermat-mod option:
		cy_i00	= tmp + 0x62;
		cy_i02	= tmp + 0x63;
		cy_i04	= tmp + 0x64;
		cy_i06	= tmp + 0x65;
		cy_i08	= tmp + 0x66;
		cy_i10	= tmp + 0x67;
		cy_i12	= tmp + 0x68;
		cy_i14	= tmp + 0x69;
		cy_i16	= tmp + 0x6a;
		cy_i18	= tmp + 0x6b;
		cy_i20	= tmp + 0x6c;
		cy_i22	= tmp + 0x6d;
		cy_i24	= tmp + 0x6e;
		cy_i26	= tmp + 0x6f;
		cy_i28	= tmp + 0x70;
		cy_i30	= tmp + 0x71;
		cy_i32	= tmp + 0x72;
		cy_i34	= tmp + 0x73;
	*/
		max_err = tmp + 0x74;
		sse2_rnd= tmp + 0x75;
		half_arr= tmp + 0x76;	/* This table needs 20x16 bytes */

		ASSERT(HERE, (r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->re == crnd && sse2_rnd->im == crnd), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr+10)->re * (half_arr+14)->re == 1.0 && (half_arr+10)->im * (half_arr+14)->im == 1.0, "thread-local memcheck failed!");

		max_err->re = 0.0;	max_err->im = 0.0;

		sign_mask = (uint64*)(r00 + radix36_creals_in_local_store);
		sse_bw  = sign_mask + 2;
		sse_sw  = sign_mask + 4;
		sse_n   = sign_mask + 6;
		bjmodn00 = (int*)(sign_mask + 8);
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
		bjmodn24 = bjmodn00 + 24;
		bjmodn25 = bjmodn00 + 25;
		bjmodn26 = bjmodn00 + 26;
		bjmodn27 = bjmodn00 + 27;
		bjmodn28 = bjmodn00 + 28;
		bjmodn29 = bjmodn00 + 29;
		bjmodn30 = bjmodn00 + 30;
		bjmodn31 = bjmodn00 + 31;
		bjmodn32 = bjmodn00 + 32;
		bjmodn33 = bjmodn00 + 33;
		bjmodn34 = bjmodn00 + 34;
		bjmodn35 = bjmodn00 + 35;

		/* Init DWT-indices: */	/* init carries	*/
		*bjmodn00 = thread_arg->bjmodn00;	cy00->re = thread_arg->cy00;
		*bjmodn01 = thread_arg->bjmodn01;	cy00->im = thread_arg->cy01;
		*bjmodn02 = thread_arg->bjmodn02;	cy02->re = thread_arg->cy02;
		*bjmodn03 = thread_arg->bjmodn03;	cy02->im = thread_arg->cy03;
		*bjmodn04 = thread_arg->bjmodn04;	cy04->re = thread_arg->cy04;
		*bjmodn05 = thread_arg->bjmodn05;	cy04->im = thread_arg->cy05;
		*bjmodn06 = thread_arg->bjmodn06;	cy06->re = thread_arg->cy06;
		*bjmodn07 = thread_arg->bjmodn07;	cy06->im = thread_arg->cy07;
		*bjmodn08 = thread_arg->bjmodn08;	cy08->re = thread_arg->cy08;
		*bjmodn09 = thread_arg->bjmodn09;	cy08->im = thread_arg->cy09;
		*bjmodn10 = thread_arg->bjmodn10;	cy10->re = thread_arg->cy10;
		*bjmodn11 = thread_arg->bjmodn11;	cy10->im = thread_arg->cy11;
		*bjmodn12 = thread_arg->bjmodn12;	cy12->re = thread_arg->cy12;
		*bjmodn13 = thread_arg->bjmodn13;	cy12->im = thread_arg->cy13;
		*bjmodn14 = thread_arg->bjmodn14;	cy14->re = thread_arg->cy14;
		*bjmodn15 = thread_arg->bjmodn15;	cy14->im = thread_arg->cy15;
		*bjmodn16 = thread_arg->bjmodn16;	cy16->re = thread_arg->cy16;
		*bjmodn17 = thread_arg->bjmodn17;	cy16->im = thread_arg->cy17;
		*bjmodn18 = thread_arg->bjmodn18;	cy18->re = thread_arg->cy18;
		*bjmodn19 = thread_arg->bjmodn19;	cy18->im = thread_arg->cy19;
		*bjmodn20 = thread_arg->bjmodn20;	cy20->re = thread_arg->cy20;
		*bjmodn21 = thread_arg->bjmodn21;	cy20->im = thread_arg->cy21;
		*bjmodn22 = thread_arg->bjmodn22;	cy22->re = thread_arg->cy22;
		*bjmodn23 = thread_arg->bjmodn23;	cy22->im = thread_arg->cy23;
		*bjmodn24 = thread_arg->bjmodn24;	cy24->re = thread_arg->cy24;
		*bjmodn25 = thread_arg->bjmodn25;	cy24->im = thread_arg->cy25;
		*bjmodn26 = thread_arg->bjmodn26;	cy26->re = thread_arg->cy26;
		*bjmodn27 = thread_arg->bjmodn27;	cy26->im = thread_arg->cy27;
		*bjmodn28 = thread_arg->bjmodn28;	cy28->re = thread_arg->cy28;
		*bjmodn29 = thread_arg->bjmodn29;	cy28->im = thread_arg->cy29;
		*bjmodn30 = thread_arg->bjmodn30;	cy30->re = thread_arg->cy30;
		*bjmodn31 = thread_arg->bjmodn31;	cy30->im = thread_arg->cy31;
		*bjmodn32 = thread_arg->bjmodn32;	cy32->re = thread_arg->cy32;
		*bjmodn33 = thread_arg->bjmodn33;	cy32->im = thread_arg->cy33;
		*bjmodn34 = thread_arg->bjmodn34;	cy34->re = thread_arg->cy34;
		*bjmodn35 = thread_arg->bjmodn35;	cy34->im = thread_arg->cy35;

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

				/* GCC-style fully-inlined ASM (64-bit only): */

				add0 = &a[j1    ];
				SSE2_RADIX36_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,r00,s1p00r,cc1);

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
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
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p28r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p32r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				l= (j+2) & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
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

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p28r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p32r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p28r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p32r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

				/*...The radix-36 DIF pass is here:	*/

				SSE2_RADIX36_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,r00,s1p00r,cc1);
			
			}	/* end for(j=_jstart; j < _jhi; j += 2) */

			jstart += nwt;
			jhi    += nwt;

			col += RADIX;
			co3 -= RADIX;
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		thread_arg->cy00 = cy00->re;
		thread_arg->cy01 = cy00->im;
		thread_arg->cy02 = cy02->re;
		thread_arg->cy03 = cy02->im;
		thread_arg->cy04 = cy04->re;
		thread_arg->cy05 = cy04->im;
		thread_arg->cy06 = cy06->re;
		thread_arg->cy07 = cy06->im;
		thread_arg->cy08 = cy08->re;
		thread_arg->cy09 = cy08->im;
		thread_arg->cy10 = cy10->re;
		thread_arg->cy11 = cy10->im;
		thread_arg->cy12 = cy12->re;
		thread_arg->cy13 = cy12->im;
		thread_arg->cy14 = cy14->re;
		thread_arg->cy15 = cy14->im;
		thread_arg->cy16 = cy16->re;
		thread_arg->cy17 = cy16->im;
		thread_arg->cy18 = cy18->re;
		thread_arg->cy19 = cy18->im;
		thread_arg->cy20 = cy20->re;
		thread_arg->cy21 = cy20->im;
		thread_arg->cy22 = cy22->re;
		thread_arg->cy23 = cy22->im;
		thread_arg->cy24 = cy24->re;
		thread_arg->cy25 = cy24->im;
		thread_arg->cy26 = cy26->re;
		thread_arg->cy27 = cy26->im;
		thread_arg->cy28 = cy28->re;
		thread_arg->cy29 = cy28->im;
		thread_arg->cy30 = cy30->re;
		thread_arg->cy31 = cy30->im;
		thread_arg->cy32 = cy32->re;
		thread_arg->cy33 = cy32->im;
		thread_arg->cy34 = cy34->re;
		thread_arg->cy35 = cy34->im;

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

