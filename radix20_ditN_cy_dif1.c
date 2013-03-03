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

#ifdef USE_SSE2

	const int radix20_creals_in_local_store = 124;

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

//	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */

	#ifdef COMPILER_TYPE_MSVC
		#include "sse2_macro.h"
	#endif

	#if defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix20_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix20_ditN_cy_dif1_gcc64.h"

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

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *r00;

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
	};

  #endif

#endif	/* USE_SSE2 */

/**************/

int radix20_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-20 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-20 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const uint32 RADIX = 20;
	const double crnd = 3.0*0x4000000*0x2000000;
	int NDIVR,i,j,j1,j2,jstart,jhi,full_pass,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19;
	static double	uc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					uc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					us1 =  0.95105651629515357211,	/*  sin(u) */
					us2 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					us3 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	static double radix_inv, n2inv;
	double scale
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19;
	double maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;

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
		static task_control_t   task_control = {NULL, (void*)cy20_process_chunk, NULL, 0x0};
	#endif

  #else
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
  #endif

	static struct complex *cc1, *ss1, *cc2, *ss2, *ss3, *max_err, *sse2_rnd, *half_arr, *tmp, *two, *r00,*r10,*r20,*r30
  #ifndef COMPILER_TYPE_GCC
		,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09
		,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19
		,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29
		,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39
	,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i
  #endif
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r;

	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19;
	static struct complex *cy00,*cy02,*cy04,*cy06,*cy08,*cy10,*cy12,*cy14,*cy16,*cy18;

  #ifdef DEBUG_SSE2
	int jt,jp;
  #endif

#else

  #if PFETCH
	double *addr, *addp;
  #endif
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double rt,it,temp,frac;
	double wt,wtinv,wtA,wtB,wtC;
	int m,m2;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19;
	double t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39;
	double
	 a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i
	,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy00 = 0x0,*_cy01 = 0x0,*_cy02 = 0x0,*_cy03 = 0x0,*_cy04 = 0x0,*_cy05 = 0x0,*_cy06 = 0x0,*_cy07 = 0x0,*_cy08 = 0x0,*_cy09 = 0x0,*_cy10 = 0x0,*_cy11 = 0x0,*_cy12 = 0x0,*_cy13 = 0x0,*_cy14 = 0x0,*_cy15 = 0x0,*_cy16 = 0x0,*_cy17 = 0x0,*_cy18 = 0x0,*_cy19 = 0x0;

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "radix20_ditN_cy_dif1: Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/20 in radix20_ditN_cy_dif1.\n",iter);
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

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix20_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_COMPLEX(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix20_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 80 16-byte slots of sc_arr for temporaries, next 5 for the nontrivial complex 16th roots,
	next 10 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
	  #ifdef COMPILER_TYPE_GCC
		r00	= sc_ptr + 0x00;		tmp	= sc_ptr + 0x28;
									s1p00r = tmp + 0x00;		two     = tmp + 0x28;
									s1p01r = tmp + 0x02;		cc1     = tmp + 0x29;
									s1p02r = tmp + 0x04;		cc2     = tmp + 0x2a;
									s1p03r = tmp + 0x06;		ss1     = tmp + 0x2b;
									s1p04r = tmp + 0x08;		ss2     = tmp + 0x2c;
		r10	= sc_ptr + 0x0a;		s1p05r = tmp + 0x0a;		ss3     = tmp + 0x2d;
									s1p06r = tmp + 0x0c;		cy00    = tmp + 0x2e;
									s1p07r = tmp + 0x0e;		cy02    = tmp + 0x2f;
									s1p08r = tmp + 0x10;		cy04    = tmp + 0x30;
									s1p09r = tmp + 0x12;		cy06    = tmp + 0x31;
		r20	= sc_ptr + 0x14;		s1p10r = tmp + 0x14;		cy08    = tmp + 0x32;
									s1p11r = tmp + 0x16;		cy10    = tmp + 0x33;
									s1p12r = tmp + 0x18;		cy12    = tmp + 0x34;
									s1p13r = tmp + 0x1a;		cy14    = tmp + 0x35;
									s1p14r = tmp + 0x1c;		cy16    = tmp + 0x36;
		r30	= sc_ptr + 0x1e;		s1p15r = tmp + 0x1e;		cy18    = tmp + 0x37;
									s1p16r = tmp + 0x20;		max_err = tmp + 0x38;
									s1p17r = tmp + 0x22;		sse2_rnd= tmp + 0x39;
									s1p18r = tmp + 0x24;		half_arr= tmp + 0x3a;	/* This table needs 20x16 bytes */
									s1p19r = tmp + 0x26;
	#else
									tmp	= sc_ptr + 0x28;
		r00	= sc_ptr + 0x00;		s1p00r = tmp + 0x00;		two     = tmp + 0x28;
		r01	= sc_ptr + 0x01;		s1p00i = tmp + 0x01;		cc1     = tmp + 0x29;
		r02	= sc_ptr + 0x02;		s1p01r = tmp + 0x02;		cc2     = tmp + 0x2a;
		r03	= sc_ptr + 0x03;		s1p01i = tmp + 0x03;		ss1     = tmp + 0x2b;
		r04	= sc_ptr + 0x04;		s1p02r = tmp + 0x04;		ss2     = tmp + 0x2c;
		r05	= sc_ptr + 0x05;		s1p02i = tmp + 0x05;		ss3     = tmp + 0x2d;
		r06	= sc_ptr + 0x06;		s1p03r = tmp + 0x06;		cy00    = tmp + 0x2e;
		r07	= sc_ptr + 0x07;		s1p03i = tmp + 0x07;		cy02    = tmp + 0x2f;
		r08	= sc_ptr + 0x08;		s1p04r = tmp + 0x08;		cy04    = tmp + 0x30;
		r09	= sc_ptr + 0x09;		s1p04i = tmp + 0x09;		cy06    = tmp + 0x31;
		r10	= sc_ptr + 0x0a;		s1p05r = tmp + 0x0a;		cy08    = tmp + 0x32;
		r11	= sc_ptr + 0x0b;		s1p05i = tmp + 0x0b;		cy10    = tmp + 0x33;
		r12	= sc_ptr + 0x0c;		s1p06r = tmp + 0x0c;		cy12    = tmp + 0x34;
		r13	= sc_ptr + 0x0d;		s1p06i = tmp + 0x0d;		cy14    = tmp + 0x35;
		r14	= sc_ptr + 0x0e;		s1p07r = tmp + 0x0e;		cy16    = tmp + 0x36;
		r15	= sc_ptr + 0x0f;		s1p07i = tmp + 0x0f;		cy18    = tmp + 0x37;
		r16	= sc_ptr + 0x10;		s1p08r = tmp + 0x10;		max_err = tmp + 0x38;
		r17	= sc_ptr + 0x11;		s1p08i = tmp + 0x11;		sse2_rnd= tmp + 0x39;
		r18	= sc_ptr + 0x12;		s1p09r = tmp + 0x12;		half_arr= tmp + 0x3a;	/* This table needs 20x16 bytes */
		r19	= sc_ptr + 0x13;		s1p09i = tmp + 0x13;
		r20	= sc_ptr + 0x14;		s1p10r = tmp + 0x14;
		r21	= sc_ptr + 0x15;		s1p10i = tmp + 0x15;
		r22	= sc_ptr + 0x16;		s1p11r = tmp + 0x16;
		r23	= sc_ptr + 0x17;		s1p11i = tmp + 0x17;
		r24	= sc_ptr + 0x18;		s1p12r = tmp + 0x18;
		r25	= sc_ptr + 0x19;		s1p12i = tmp + 0x19;
		r26	= sc_ptr + 0x1a;		s1p13r = tmp + 0x1a;
		r27	= sc_ptr + 0x1b;		s1p13i = tmp + 0x1b;
		r28	= sc_ptr + 0x1c;		s1p14r = tmp + 0x1c;
		r29	= sc_ptr + 0x1d;		s1p14i = tmp + 0x1d;
		r30	= sc_ptr + 0x1e;		s1p15r = tmp + 0x1e;
		r31	= sc_ptr + 0x1f;		s1p15i = tmp + 0x1f;
		r32	= sc_ptr + 0x20;		s1p16r = tmp + 0x20;
		r33	= sc_ptr + 0x21;		s1p16i = tmp + 0x21;
		r34	= sc_ptr + 0x22;		s1p17r = tmp + 0x22;
		r35	= sc_ptr + 0x23;		s1p17i = tmp + 0x23;
		r36	= sc_ptr + 0x24;		s1p18r = tmp + 0x24;
		r37	= sc_ptr + 0x25;		s1p18i = tmp + 0x25;
		r38	= sc_ptr + 0x26;		s1p19r = tmp + 0x26;
		r39	= sc_ptr + 0x27;		s1p19i = tmp + 0x27;
	#endif
		/* These remain fixed: */
		two->re = two->im = 2.0;
		cc1->re = cc1->im = uc1;
		cc2->re = cc2->im = uc2;
		ss1->re = ss1->im = us1;
		ss2->re = ss2->im = us2;
		ss3->re = ss3->im = us3;

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		sse2_rnd->re = crnd;
		sse2_rnd->im = crnd;

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

	  #ifdef USE_PTHREAD
		r00 = __r0 + cslots_in_local_store;
		/* Init thread 1-CY_THREADS's local stores and pointers: */
		for(i = 1; i < CY_THREADS; ++i) {
			/* Only care about the constants for each thread here, but easier to just copy the entire thread0 local store: */
			memcpy(r00, __r0, cslots_in_local_store<<4);	// bytewise copy treats complex and uint64 subdata the same
			r00 += cslots_in_local_store;
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
			free((void *)_bjmodn12); _bjmodn12 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;

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
			free((void *)_cy12); _cy12 = 0x0;
			free((void *)_cy13); _cy13 = 0x0;
			free((void *)_cy14); _cy14 = 0x0;
			free((void *)_cy15); _cy15 = 0x0;
			free((void *)_cy16); _cy16 = 0x0;
			free((void *)_cy17); _cy17 = 0x0;
			free((void *)_cy18); _cy18 = 0x0;
			free((void *)_cy19); _cy19 = 0x0;

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
		_cy12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy12== 0x0);
		_cy13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy13== 0x0);
		_cy14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy14== 0x0);
		_cy15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy15== 0x0);
		_cy16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy16== 0x0);
		_cy17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy17== 0x0);
		_cy18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy18== 0x0);
		_cy19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy19== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix20_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/20-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix20_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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

/*...The radix-20 final DIT pass is here.	*/

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
		_cy12[ithread] = 0;
		_cy13[ithread] = 0;
		_cy14[ithread] = 0;
		_cy15[ithread] = 0;
		_cy16[ithread] = 0;
		_cy17[ithread] = 0;
		_cy18[ithread] = 0;
		_cy19[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy00[      0] = -2;
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

		_jstart[ithread] = ithread*NDIVR/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
	}

    /* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
    if(!full_pass)
    {
		khi = 1;
	#if FFT_DEBUG
		fprintf(dbg_file, "radix20_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifdef USE_OMP
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,l,col,co2,co3,m,m2,\
		n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,rt,it,\
		t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,\
		a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,\
		a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,\
		bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,\
		cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19\
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
		ASSERT(HERE, tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].r00 + 0x28 + 0x39;	// sse2_rnd;
		ASSERT(HERE, ((tmp + 0)->re == crnd && (tmp + 0)->im == crnd), "thread-local memcheck failed!");
		ASSERT(HERE, (tmp + 1+10)->re * (tmp + 1+14)->re == 1.0 && (tmp + 1+10)->im * (tmp + 1+14)->im == 1.0, "thread-local memcheck failed!");

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
		tdat[ithread].cy12 = _cy12[ithread];
		tdat[ithread].cy13 = _cy13[ithread];
		tdat[ithread].cy14 = _cy14[ithread];
		tdat[ithread].cy15 = _cy15[ithread];
		tdat[ithread].cy16 = _cy16[ithread];
		tdat[ithread].cy17 = _cy17[ithread];
		tdat[ithread].cy18 = _cy18[ithread];
		tdat[ithread].cy19 = _cy19[ithread];
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

	#ifdef USE_SSE2									/* init carries	*/
		*bjmodn00 = _bjmodn00[ithread];				cy00->re = _cy00[ithread];
		*bjmodn01 = _bjmodn01[ithread];				cy00->im = _cy01[ithread];
		*bjmodn02 = _bjmodn02[ithread];				cy02->re = _cy02[ithread];
		*bjmodn03 = _bjmodn03[ithread];				cy02->im = _cy03[ithread];
		*bjmodn04 = _bjmodn04[ithread];				cy04->re = _cy04[ithread];
		*bjmodn05 = _bjmodn05[ithread];				cy04->im = _cy05[ithread];
		*bjmodn06 = _bjmodn06[ithread];				cy06->re = _cy06[ithread];
		*bjmodn07 = _bjmodn07[ithread];				cy06->im = _cy07[ithread];
		*bjmodn08 = _bjmodn08[ithread];				cy08->re = _cy08[ithread];
		*bjmodn09 = _bjmodn09[ithread];				cy08->im = _cy09[ithread];
		*bjmodn10 = _bjmodn10[ithread];				cy10->re = _cy10[ithread];
		*bjmodn11 = _bjmodn11[ithread];				cy10->im = _cy11[ithread];
		*bjmodn12 = _bjmodn12[ithread];				cy12->re = _cy12[ithread];
		*bjmodn13 = _bjmodn13[ithread];				cy12->im = _cy13[ithread];
		*bjmodn14 = _bjmodn14[ithread];				cy14->re = _cy14[ithread];
		*bjmodn15 = _bjmodn15[ithread];				cy14->im = _cy15[ithread];
		*bjmodn16 = _bjmodn16[ithread];				cy16->re = _cy16[ithread];
		*bjmodn17 = _bjmodn17[ithread];				cy16->im = _cy17[ithread];
		*bjmodn18 = _bjmodn18[ithread];				cy18->re = _cy18[ithread];
		*bjmodn19 = _bjmodn19[ithread];				cy18->im = _cy19[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];				cy00 = _cy00[ithread];
		bjmodn01 = _bjmodn01[ithread];				cy01 = _cy01[ithread];
		bjmodn02 = _bjmodn02[ithread];				cy02 = _cy02[ithread];
		bjmodn03 = _bjmodn03[ithread];				cy03 = _cy03[ithread];
		bjmodn04 = _bjmodn04[ithread];				cy04 = _cy04[ithread];
		bjmodn05 = _bjmodn05[ithread];				cy05 = _cy05[ithread];
		bjmodn06 = _bjmodn06[ithread];				cy06 = _cy06[ithread];
		bjmodn07 = _bjmodn07[ithread];				cy07 = _cy07[ithread];
		bjmodn08 = _bjmodn08[ithread];				cy08 = _cy08[ithread];
		bjmodn09 = _bjmodn09[ithread];				cy09 = _cy09[ithread];
		bjmodn10 = _bjmodn10[ithread];				cy10 = _cy10[ithread];
		bjmodn11 = _bjmodn11[ithread];				cy11 = _cy11[ithread];
		bjmodn12 = _bjmodn12[ithread];				cy12 = _cy12[ithread];
		bjmodn13 = _bjmodn13[ithread];				cy13 = _cy13[ithread];
		bjmodn14 = _bjmodn14[ithread];				cy14 = _cy14[ithread];
		bjmodn15 = _bjmodn15[ithread];				cy15 = _cy15[ithread];
		bjmodn16 = _bjmodn16[ithread];				cy16 = _cy16[ithread];
		bjmodn17 = _bjmodn17[ithread];				cy17 = _cy17[ithread];
		bjmodn18 = _bjmodn18[ithread];				cy18 = _cy18[ithread];
		bjmodn19 = _bjmodn19[ithread];				cy19 = _cy19[ithread];
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
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p04] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p04] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p05] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p05] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p06] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p06] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p07] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p07] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p08] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p08] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p09] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p09] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p10] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p10] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p11] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p11] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p12] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p12] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p13] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p13] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p14] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p14] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p15] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p15] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p16] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p16] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p17] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p17] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p18] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p18] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p19] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p19] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix20_ditN_cy_dif1 inputs:\n");
			fprintf(stderr, "A_in[00] = %20.5f, %20.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_in[01] = %20.5f, %20.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_in[02] = %20.5f, %20.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_in[03] = %20.5f, %20.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_in[04] = %20.5f, %20.5f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_in[05] = %20.5f, %20.5f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_in[06] = %20.5f, %20.5f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_in[07] = %20.5f, %20.5f\n",a[jt+p07],a[jp+p07]);
			fprintf(stderr, "A_in[08] = %20.5f, %20.5f\n",a[jt+p08],a[jp+p08]);
			fprintf(stderr, "A_in[09] = %20.5f, %20.5f\n",a[jt+p09],a[jp+p09]);
			fprintf(stderr, "A_in[10] = %20.5f, %20.5f\n",a[jt+p10],a[jp+p10]);
			fprintf(stderr, "A_in[11] = %20.5f, %20.5f\n",a[jt+p11],a[jp+p11]);
			fprintf(stderr, "A_in[12] = %20.5f, %20.5f\n",a[jt+p12],a[jp+p12]);
			fprintf(stderr, "A_in[13] = %20.5f, %20.5f\n",a[jt+p13],a[jp+p13]);
			fprintf(stderr, "A_in[14] = %20.5f, %20.5f\n",a[jt+p14],a[jp+p14]);
			fprintf(stderr, "A_in[15] = %20.5f, %20.5f\n",a[jt+p15],a[jp+p15]);
			fprintf(stderr, "A_in[16] = %20.5f, %20.5f\n",a[jt+p16],a[jp+p16]);
			fprintf(stderr, "A_in[17] = %20.5f, %20.5f\n",a[jt+p17],a[jp+p17]);
			fprintf(stderr, "A_in[18] = %20.5f, %20.5f\n",a[jt+p18],a[jp+p18]);
			fprintf(stderr, "A_in[19] = %20.5f, %20.5f\n",a[jt+p19],a[jp+p19]);
			fprintf(stderr, "\n");
		#endif

			/*
			!...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do a radix-20 DIT transform...
			*/
	#ifdef USE_SSE2

			add0 = &a[j1    ];
		//	add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p03;	add3 = add0+p02;

		#if defined(COMPILER_TYPE_MSVC)

			/* Outputs in SSE2 modes are temps 2*5*16 = 10*16 = 0x0a0 bytes apart: */
			__asm	mov	eax, add0	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C */
			__asm	mov	edx, add0
			__asm	mov	esi, p01	/* esi will store power-of-2 multiples of p01 throughout */
			__asm	shl	esi, 3		/* Pointer offset for floating doubles */
			__asm	add edx, esi
			__asm	mov ebx, edx	/* add1 = add0+p01 */
			__asm	add edx, esi
			__asm	mov ecx, edx	/* add3 = add0+p02 */
			__asm	add edx, esi	/* add2 = add0+p03 */
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x0a0, 0x140, r00)

		/*	add0,1,2,3 = &a[j1+p04]+p3,2,1,0 */
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	add	eax, esi	/* &a[j1+p04] */
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x0a0, 0x140, r02)

		/*	add0,1,2,3 = &a[j1+p08]+p1,0,2,3 */
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	add	eax, esi	/* &a[j1+p08] */
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x0a0, 0x140, r04)

		/*	add0,1,2,3 = &a[j1+p12]+p2,3,0,1 */
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	add	eax, esi	/* &a[j1+p12] */
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x0a0, 0x140, r06)

		/*	add0,1,2,3 = &a[j1+p16]+p0,1,3,2 */
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	add	eax, esi	/* &a[j1+p16] */
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x0a0, 0x140, r08)

			/* Radix-5 DFT uses adjacent temps, i.e. stride = 2*16 bytes: */
			SSE2_RADIX_05_DFT_0TWIDDLE(r00,r02,r04,r06,r08,cc1,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r)
			SSE2_RADIX_05_DFT_0TWIDDLE(r10,r12,r14,r16,r18,cc1,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r)
			SSE2_RADIX_05_DFT_0TWIDDLE(r20,r22,r24,r26,r28,cc1,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r)
			SSE2_RADIX_05_DFT_0TWIDDLE(r30,r32,r34,r36,r38,cc1,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r)

		#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

			SSE2_RADIX20_DIT_NOTWIDDLE(add0,p01,p04,r00,r10,r20,r30,cc1,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r);
#if 0
	if(j < 2 && !full_pass && iter <= 10)
	{
  #ifdef USE_SSE2
		fprintf(stderr, "Iter %3d: err_re = %10.8f, err_im = %10.8f, maxerr = %10.8f\n", iter, max_err->re,max_err->im, MAX(max_err->re,max_err->im));
		fprintf(stderr, "sp00-sp01 = %20.5f %20.5f %20.5f %20.5f\n",s1p00r->re,s1p00r->im,s1p01r->re,s1p01r->im);
		fprintf(stderr, "sp02-sp03 = %20.5f %20.5f %20.5f %20.5f\n",s1p02r->re,s1p02r->im,s1p03r->re,s1p03r->im);
		fprintf(stderr, "sp04-sp05 = %20.5f %20.5f %20.5f %20.5f\n",s1p04r->re,s1p04r->im,s1p05r->re,s1p05r->im);
		fprintf(stderr, "sp06-sp07 = %20.5f %20.5f %20.5f %20.5f\n",s1p06r->re,s1p06r->im,s1p07r->re,s1p07r->im);
		fprintf(stderr, "sp08-sp09 = %20.5f %20.5f %20.5f %20.5f\n",s1p08r->re,s1p08r->im,s1p09r->re,s1p09r->im);
		fprintf(stderr, "sp10-sp11 = %20.5f %20.5f %20.5f %20.5f\n",s1p10r->re,s1p10r->im,s1p11r->re,s1p11r->im);
		fprintf(stderr, "sp12-sp13 = %20.5f %20.5f %20.5f %20.5f\n",s1p12r->re,s1p12r->im,s1p13r->re,s1p13r->im);
		fprintf(stderr, "sp14-sp15 = %20.5f %20.5f %20.5f %20.5f\n",s1p14r->re,s1p14r->im,s1p15r->re,s1p15r->im);
		fprintf(stderr, "sp16-sp17 = %20.5f %20.5f %20.5f %20.5f\n",s1p16r->re,s1p16r->im,s1p17r->re,s1p17r->im);
		fprintf(stderr, "sp18-sp19 = %20.5f %20.5f %20.5f %20.5f\n",s1p18r->re,s1p18r->im,s1p19r->re,s1p19r->im);
	//	fprintf(stderr, "Iter %3d: a0-3_in = %20.5f %20.5f %20.5f %20.5f, cy0-3 = %20.5f %20.5f %20.5f %20.5f, maxerr = %10.8f\n"
	//	,iter,s1p00r->re,s1p01r->re,s1p00r->im,s1p01r->im,cy00->re,cy00->im,cy02->re,cy02->im,MAX(max_err->re,max_err->im));
  #else
		fprintf(stderr, "Iter %3d: a0-3_in = %20.5f %20.5f %20.5f %20.5f, cy0-3 = %20.5f %20.5f %20.5f %20.5f, maxerr = %10.8f\n"
		,iter,a1p00r,a1p00i,a1p1r,a1p1i,cy00,cy01,cy02,cy03,maxerr);
  #endif
		fflush(stderr);
	}
#endif

		#endif

	#else
							 /*          inputs           */ /*                                      outputs                                      */
			RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00,t01,t10,t11,t20,t21,t30,t31,rt,it);
			RADIX_04_DIT(a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t02,t03,t12,t13,t22,t23,t32,t33,rt,it);
			RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],t04,t05,t14,t15,t24,t25,t34,t35,rt,it);
			RADIX_04_DIT(a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],t06,t07,t16,t17,t26,t27,t36,t37,rt,it);
			RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],t08,t09,t18,t19,t28,t29,t38,t39,rt,it);

		/*...and now do 4 radix-5 transforms...*/
							 /*                                                inputs                                                   */ /*                 outputs                   */
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,rt,it);

	#endif

	#ifdef DEBUG_SSE2
		fprintf(stderr, "radix20_dit: R00= %20.5f, %20.5f\n",r00->re,r01->re);
		fprintf(stderr, "radix20_dit: R02= %20.5f, %20.5f\n",r02->re,r03->re);
		fprintf(stderr, "radix20_dit: R04= %20.5f, %20.5f\n",r04->re,r05->re);
		fprintf(stderr, "radix20_dit: R06= %20.5f, %20.5f\n",r06->re,r07->re);
		fprintf(stderr, "radix20_dit: R08= %20.5f, %20.5f\n",r08->re,r09->re);
		fprintf(stderr, "radix20_dit: R10= %20.5f, %20.5f\n",r10->re,r11->re);
		fprintf(stderr, "radix20_dit: R12= %20.5f, %20.5f\n",r12->re,r13->re);
		fprintf(stderr, "radix20_dit: R14= %20.5f, %20.5f\n",r14->re,r15->re);
		fprintf(stderr, "radix20_dit: R16= %20.5f, %20.5f\n",r16->re,r17->re);
		fprintf(stderr, "radix20_dit: R18= %20.5f, %20.5f\n",r18->re,r19->re);
		fprintf(stderr, "radix20_dit: R20= %20.5f, %20.5f\n",r20->re,r21->re);
		fprintf(stderr, "radix20_dit: R22= %20.5f, %20.5f\n",r22->re,r23->re);
		fprintf(stderr, "radix20_dit: R24= %20.5f, %20.5f\n",r24->re,r25->re);
		fprintf(stderr, "radix20_dit: R26= %20.5f, %20.5f\n",r26->re,r27->re);
		fprintf(stderr, "radix20_dit: R28= %20.5f, %20.5f\n",r28->re,r29->re);
		fprintf(stderr, "radix20_dit: R30= %20.5f, %20.5f\n",r30->re,r31->re);
		fprintf(stderr, "radix20_dit: R32= %20.5f, %20.5f\n",r32->re,r33->re);
		fprintf(stderr, "radix20_dit: R34= %20.5f, %20.5f\n",r34->re,r35->re);
		fprintf(stderr, "radix20_dit: R36= %20.5f, %20.5f\n",r36->re,r37->re);
		fprintf(stderr, "radix20_dit: R38= %20.5f, %20.5f\n",r38->re,r39->re);
		fprintf(stderr, "\n");
		fprintf(stderr, "radix20_dit: s00= %20.5f, %20.5f\n",s1p00r->re,s1p00i->re);
		fprintf(stderr, "radix20_dit: s01= %20.5f, %20.5f\n",s1p01r->re,s1p01i->re);
		fprintf(stderr, "radix20_dit: s02= %20.5f, %20.5f\n",s1p02r->re,s1p02i->re);
		fprintf(stderr, "radix20_dit: s03= %20.5f, %20.5f\n",s1p03r->re,s1p03i->re);
		fprintf(stderr, "radix20_dit: s04= %20.5f, %20.5f\n",s1p04r->re,s1p04i->re);
		fprintf(stderr, "radix20_dit: s05= %20.5f, %20.5f\n",s1p05r->re,s1p05i->re);
		fprintf(stderr, "radix20_dit: s06= %20.5f, %20.5f\n",s1p06r->re,s1p06i->re);
		fprintf(stderr, "radix20_dit: s07= %20.5f, %20.5f\n",s1p07r->re,s1p07i->re);
		fprintf(stderr, "radix20_dit: s08= %20.5f, %20.5f\n",s1p08r->re,s1p08i->re);
		fprintf(stderr, "radix20_dit: s09= %20.5f, %20.5f\n",s1p09r->re,s1p09i->re);
		fprintf(stderr, "radix20_dit: s10= %20.5f, %20.5f\n",s1p10r->re,s1p10i->re);
		fprintf(stderr, "radix20_dit: s11= %20.5f, %20.5f\n",s1p11r->re,s1p11i->re);
		fprintf(stderr, "radix20_dit: s12= %20.5f, %20.5f\n",s1p12r->re,s1p12i->re);
		fprintf(stderr, "radix20_dit: s13= %20.5f, %20.5f\n",s1p13r->re,s1p13i->re);
		fprintf(stderr, "radix20_dit: s14= %20.5f, %20.5f\n",s1p14r->re,s1p14i->re);
		fprintf(stderr, "radix20_dit: s15= %20.5f, %20.5f\n",s1p15r->re,s1p15i->re);
		fprintf(stderr, "radix20_dit: s16= %20.5f, %20.5f\n",s1p16r->re,s1p16i->re);
		fprintf(stderr, "radix20_dit: s17= %20.5f, %20.5f\n",s1p17r->re,s1p17i->re);
		fprintf(stderr, "radix20_dit: s18= %20.5f, %20.5f\n",s1p18r->re,s1p18i->re);
		fprintf(stderr, "radix20_dit: s19= %20.5f, %20.5f\n",s1p19r->re,s1p19i->re);
		exit(0);
	#endif

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 20 separate blocks of the A-array, we need 20 separate carries.	*/

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
			SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
			SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
			SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
			SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
			SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
		  #else
			SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
		  #endif

		  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		  #ifdef ERR_CHECK_ALL
			SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		  #else
			SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
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

		  #ifdef ERR_CHECK_ALL
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00);
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04);
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08);
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12);
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16);
		  #else
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00);
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04);
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08);
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12);
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16);
		  #endif

		  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		  #ifdef ERR_CHECK_ALL
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		  #else
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		  #endif

		  #endif

			i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

		#else	/* USE_SSE2 */

		/*...set0 is slightly different from others:	*/
			cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00   );
			cmplx_carry_norm_errcheck (a1p01r,a1p01i,cy01,bjmodn01,1 );
			cmplx_carry_norm_errcheck (a1p02r,a1p02i,cy02,bjmodn02,2 );
			cmplx_carry_norm_errcheck (a1p03r,a1p03i,cy03,bjmodn03,3 );
			cmplx_carry_norm_errcheck (a1p04r,a1p04i,cy04,bjmodn04,4 );
			cmplx_carry_norm_errcheck (a1p05r,a1p05i,cy05,bjmodn05,5 );
			cmplx_carry_norm_errcheck (a1p06r,a1p06i,cy06,bjmodn06,6 );
			cmplx_carry_norm_errcheck (a1p07r,a1p07i,cy07,bjmodn07,7 );
			cmplx_carry_norm_errcheck (a1p08r,a1p08i,cy08,bjmodn08,8 );
			cmplx_carry_norm_errcheck (a1p09r,a1p09i,cy09,bjmodn09,9 );
			cmplx_carry_norm_errcheck (a1p10r,a1p10i,cy10,bjmodn10,10);
			cmplx_carry_norm_errcheck (a1p11r,a1p11i,cy11,bjmodn11,11);
			cmplx_carry_norm_errcheck (a1p12r,a1p12i,cy12,bjmodn12,12);
			cmplx_carry_norm_errcheck (a1p13r,a1p13i,cy13,bjmodn13,13);
			cmplx_carry_norm_errcheck (a1p14r,a1p14i,cy14,bjmodn14,14);
			cmplx_carry_norm_errcheck (a1p15r,a1p15i,cy15,bjmodn15,15);
			cmplx_carry_norm_errcheck (a1p16r,a1p16i,cy16,bjmodn16,16);
			cmplx_carry_norm_errcheck (a1p17r,a1p17i,cy17,bjmodn17,17);
			cmplx_carry_norm_errcheck (a1p18r,a1p18i,cy18,bjmodn18,18);
			cmplx_carry_norm_errcheck (a1p19r,a1p19i,cy19,bjmodn19,19);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	/* USE_SSE2 */

	#ifdef USE_SSE2

	#ifdef DEBUG_SSE2
		fprintf(stderr, "radix20_carry_out: s00= %20.5f, %20.5f\n",s1p00r->re,s1p00i-> re);
		fprintf(stderr, "radix20_carry_out: s01= %20.5f, %20.5f\n",s1p01r->re,s1p01i-> re);
		fprintf(stderr, "radix20_carry_out: s02= %20.5f, %20.5f\n",s1p02r->re,s1p02i-> re);
		fprintf(stderr, "radix20_carry_out: s03= %20.5f, %20.5f\n",s1p03r->re,s1p03i-> re);
		fprintf(stderr, "radix20_carry_out: s04= %20.5f, %20.5f\n",s1p04r->re,s1p04i-> re);
		fprintf(stderr, "radix20_carry_out: s05= %20.5f, %20.5f\n",s1p05r->re,s1p05i-> re);
		fprintf(stderr, "radix20_carry_out: s06= %20.5f, %20.5f\n",s1p06r->re,s1p06i-> re);
		fprintf(stderr, "radix20_carry_out: s07= %20.5f, %20.5f\n",s1p07r->re,s1p07i-> re);
		fprintf(stderr, "radix20_carry_out: s08= %20.5f, %20.5f\n",s1p08r->re,s1p08i-> re);
		fprintf(stderr, "radix20_carry_out: s09= %20.5f, %20.5f\n",s1p09r->re,s1p09i-> re);
		fprintf(stderr, "radix20_carry_out: s10= %20.5f, %20.5f\n",s1p10r->re,s1p10i-> re);
		fprintf(stderr, "radix20_carry_out: s11= %20.5f, %20.5f\n",s1p11r->re,s1p11i-> re);
		fprintf(stderr, "radix20_carry_out: s12= %20.5f, %20.5f\n",s1p12r->re,s1p12i-> re);
		fprintf(stderr, "radix20_carry_out: s13= %20.5f, %20.5f\n",s1p13r->re,s1p13i-> re);
		fprintf(stderr, "radix20_carry_out: s14= %20.5f, %20.5f\n",s1p14r->re,s1p14i-> re);
		fprintf(stderr, "radix20_carry_out: s15= %20.5f, %20.5f\n",s1p15r->re,s1p15i-> re);
		fprintf(stderr, "radix20_carry_out: s16= %20.5f, %20.5f\n",s1p16r->re,s1p16i-> re);
		fprintf(stderr, "radix20_carry_out: s17= %20.5f, %20.5f\n",s1p17r->re,s1p17i-> re);
		fprintf(stderr, "radix20_carry_out: s18= %20.5f, %20.5f\n",s1p18r->re,s1p18i-> re);
		fprintf(stderr, "radix20_carry_out: s19= %20.5f, %20.5f\n",s1p19r->re,s1p19i-> re);
		fprintf(stderr, "\n");
	#endif

		#if defined(COMPILER_TYPE_MSVC)

			/* Index patterns of s1p-terms here are the same as in DIT DFT: */
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,cc1,r00,r02,r04,r08,r06)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,cc1,r10,r12,r14,r18,r16)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,cc1,r20,r22,r24,r28,r26)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p05r,s1p01r,s1p17r,s1p13r,s1p09r,cc1,r30,r32,r34,r38,r36)

		//	add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p03;	add3 = add0+p02;
			__asm	mov	eax, add0	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */
			__asm	mov	edx, add0
			__asm	mov	esi, p01	/* esi will store power-of-2 multiples of p01 throughout */
			__asm	shl	esi, 3		/* Pointer offset for floating doubles */
			__asm	add edx, esi
			__asm	mov ebx, edx	/* add0+p01 */
			__asm	add edx, esi
			__asm	mov ecx, edx	/* add0+p02 */
			__asm	add edx, esi	/* add0+p03 */
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x0a0, 0x140, eax,ebx,edx,ecx)

		//	add0,1,2,3 = &a[j1+p16]+p0,1,3,2
			__asm	mov	esi, p16
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p16]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, 0x0a0, 0x140, eax,ebx,edx,ecx)

		//	add0,1,2,3 = &a[j1+p12]+p2,3,0,1
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p12]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, 0x0a0, 0x140, ecx,edx,eax,ebx)

		//	add0,1,2,3 = &a[j1+p04]+p3,2,1,0
			__asm	mov	esi, p08
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p04]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, 0x0a0, 0x140, edx,ecx,ebx,eax)

		//	add0,1,2,3 = &a[j1+p08]+p1,0,2,3
			__asm	mov	esi, p04
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p08]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, 0x0a0, 0x140, ebx,eax,ecx,edx)

		#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

			SSE2_RADIX20_DIF_NOTWIDDLE(add0,p01,p04,p08,p16,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r,cc1,r00,r10,r20,r30)

		#endif

		#ifdef DEBUG_SSE2
			fprintf(stderr, "radix20: A_out[00] = %20.5f, %20.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix20: A_out[01] = %20.5f, %20.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix20: A_out[02] = %20.5f, %20.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix20: A_out[03] = %20.5f, %20.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "radix20: A_out[04] = %20.5f, %20.5f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "radix20: A_out[05] = %20.5f, %20.5f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "radix20: A_out[06] = %20.5f, %20.5f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "radix20: A_out[07] = %20.5f, %20.5f\n",a[jt+p07],a[jp+p07]);
			fprintf(stderr, "radix20: A_out[08] = %20.5f, %20.5f\n",a[jt+p08],a[jp+p08]);
			fprintf(stderr, "radix20: A_out[09] = %20.5f, %20.5f\n",a[jt+p09],a[jp+p09]);
			fprintf(stderr, "radix20: A_out[10] = %20.5f, %20.5f\n",a[jt+p10],a[jp+p10]);
			fprintf(stderr, "radix20: A_out[11] = %20.5f, %20.5f\n",a[jt+p11],a[jp+p11]);
			fprintf(stderr, "radix20: A_out[12] = %20.5f, %20.5f\n",a[jt+p12],a[jp+p12]);
			fprintf(stderr, "radix20: A_out[13] = %20.5f, %20.5f\n",a[jt+p13],a[jp+p13]);
			fprintf(stderr, "radix20: A_out[14] = %20.5f, %20.5f\n",a[jt+p14],a[jp+p14]);
			fprintf(stderr, "radix20: A_out[15] = %20.5f, %20.5f\n",a[jt+p15],a[jp+p15]);
			fprintf(stderr, "radix20: A_out[16] = %20.5f, %20.5f\n",a[jt+p16],a[jp+p16]);
			fprintf(stderr, "radix20: A_out[17] = %20.5f, %20.5f\n",a[jt+p17],a[jp+p17]);
			fprintf(stderr, "radix20: A_out[18] = %20.5f, %20.5f\n",a[jt+p18],a[jp+p18]);
			fprintf(stderr, "radix20: A_out[19] = %20.5f, %20.5f\n",a[jt+p19],a[jp+p19]);
			exit(0);
		#endif

	#else

		/*...The radix-20 DIF pass is here:	*/
		#if PFETCH
			addr = &a[j1];
			prefetch_p_doubles(addr);
		#endif

		/*...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do 4 radix-5 transforms...*/
							 /*                                                inputs                                                   */ /*                 outputs                   */
		#if PFETCH																																/*[--y3-] [--y4-] <<<<< swap last 2 outputs to undo swap of these in macro */
			RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it,addr,addp,p01,p02);
			RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it,addr,addp,p03,p04);
			RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it,addr,addp,p05,p06);
			RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it,addr,addp,p07,p08);
		#else
			RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it);
			RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it);
			RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it);
			RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it);
		#endif

		/*...and now do 5 radix-4 transforms...*/
							 /*          inputs           */ /*                                      outputs                                      */
		#if PFETCH
			addp = addr+p09;
			prefetch_p_doubles(addp);

			RADIX_04_DIF_PFETCH(t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it,addr,addp,p10,p11);
			RADIX_04_DIF_PFETCH(t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it,addr,addp,p12,p13);
			RADIX_04_DIF_PFETCH(t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it,addr,addp,p14,p15);
			RADIX_04_DIF_PFETCH(t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it,addr,addp,p16,p17);
			RADIX_04_DIF_PFETCH(t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it,addr,addp,p18,p19);
		#else
			RADIX_04_DIF       (t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it);
			RADIX_04_DIF       (t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it);
			RADIX_04_DIF       (t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
			RADIX_04_DIF       (t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
			RADIX_04_DIF       (t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it);
		#endif

	#endif	/* USE_SSE2 */
			}

			jstart += nwt;
			jhi    += nwt;
			col += RADIX;
			co3 -= RADIX;

		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_SSE2
		_cy00[ithread] = cy00->re;
		_cy01[ithread] = cy00->im;
		_cy02[ithread] = cy02->re;
		_cy03[ithread] = cy02->im;
		_cy04[ithread] = cy04->re;
		_cy05[ithread] = cy04->im;
		_cy06[ithread] = cy06->re;
		_cy07[ithread] = cy06->im;
		_cy08[ithread] = cy08->re;
		_cy09[ithread] = cy08->im;
		_cy10[ithread] = cy10->re;
		_cy11[ithread] = cy10->im;
		_cy12[ithread] = cy12->re;
		_cy13[ithread] = cy12->im;
		_cy14[ithread] = cy14->re;
		_cy15[ithread] = cy14->im;
		_cy16[ithread] = cy16->re;
		_cy17[ithread] = cy16->im;
		_cy18[ithread] = cy18->re;
		_cy19[ithread] = cy18->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
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
		_cy12[ithread] = cy12;
		_cy13[ithread] = cy13;
		_cy14[ithread] = cy14;
		_cy15[ithread] = cy15;
		_cy16[ithread] = cy16;
		_cy17[ithread] = cy17;
		_cy18[ithread] = cy18;
		_cy19[ithread] = cy19;
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
		ASSERT(HERE, 0x0 == cy20_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that

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
		_cy12[ithread] = tdat[ithread].cy12;
		_cy13[ithread] = tdat[ithread].cy13;
		_cy14[ithread] = tdat[ithread].cy14;
		_cy15[ithread] = tdat[ithread].cy15;
		_cy16[ithread] = tdat[ithread].cy16;
		_cy17[ithread] = tdat[ithread].cy17;
		_cy18[ithread] = tdat[ithread].cy18;
		_cy19[ithread] = tdat[ithread].cy19;
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-20 forward DIF FFT of the first block of 20 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 20 outputs of (1);
	!   (3) Reweight and perform a radix-20 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 20 elements and repeat (1-4).
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
	t12= _cy12[CY_THREADS - 1];
	t13= _cy13[CY_THREADS - 1];
	t14= _cy14[CY_THREADS - 1];
	t15= _cy15[CY_THREADS - 1];
	t16= _cy16[CY_THREADS - 1];
	t17= _cy17[CY_THREADS - 1];
	t18= _cy18[CY_THREADS - 1];
	t19= _cy19[CY_THREADS - 1];

	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		ASSERT(HERE, CY_THREADS > 1,"radix20_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
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
		_cy12[ithread] = _cy12[ithread-1];
		_cy13[ithread] = _cy13[ithread-1];
		_cy14[ithread] = _cy14[ithread-1];
		_cy15[ithread] = _cy15[ithread-1];
		_cy16[ithread] = _cy16[ithread-1];
		_cy17[ithread] = _cy17[ithread-1];
		_cy18[ithread] = _cy18[ithread-1];
		_cy19[ithread] = _cy19[ithread-1];
	}

	_cy00[0] =+t19;	/* ...The wraparound carry is here: */
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
	_cy12[0] = t11;
	_cy13[0] = t12;
	_cy14[0] = t13;
	_cy15[0] = t14;
	_cy16[0] = t15;
	_cy17[0] = t16;
	_cy18[0] = t17;
	_cy19[0] = t18;

	full_pass = 0;
	scale = 1;

	jhi = 7;

    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		for(j = ithread*pini; j <= ithread*pini + jhi; j++)
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
		}
	}
}	/* endfor(outer) */

    t00 = 0;
    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		t00 += fabs(_cy00[0])+fabs(_cy01[0])+fabs(_cy02[0])+fabs(_cy03[0])+fabs(_cy04[0])+fabs(_cy05[0])+fabs(_cy06[0])+fabs(_cy07[0])+fabs(_cy08[0])+fabs(_cy09[0])+fabs(_cy10[0])+fabs(_cy11[0])+fabs(_cy12[0])+fabs(_cy13[0])+fabs(_cy14[0])+fabs(_cy15[0])+fabs(_cy16[0])+fabs(_cy17[0])+fabs(_cy18[0])+fabs(_cy19[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
    }

	if(t00 != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix20_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

int radix20_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter                 , uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-20 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-20 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n20, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09
		,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double uc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					uc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					us1 =  0.95105651629515357211,	/*  sin(u) */
					us2 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					us3 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	static double radix_inv, n2inv;
	double rt,it
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i
		,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,temp,scale;
#if PFETCH
	double *addr, *addp;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "radix20_ditN_cy_dif1_nochk: Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

/*...change n20 and n_div_wt to non-static to work around a gcc compiler bug. */
	n20   = n/20;
	n_div_nwt = n20 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n20)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/20 in radix20_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)20));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p01 = n20;
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

		bjmodnini=0;
		for(j=0; j < n20; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-20 final DIT pass is here.	*/

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

	col=0;
	co2=(n >> nwt_bits)-1+20;
	co3=co2-20;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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
	!...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do a radix-20 DIT transform...
	*/
	#if 1
		/*...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do 5 radix-4 transforms...*/
							 /*          inputs           */ /*                                      outputs                                      */
			RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00,t01,t10,t11,t20,t21,t30,t31,rt,it);
			RADIX_04_DIT(a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t02,t03,t12,t13,t22,t23,t32,t33,rt,it);
			RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],t04,t05,t14,t15,t24,t25,t34,t35,rt,it);
			RADIX_04_DIT(a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],t06,t07,t16,t17,t26,t27,t36,t37,rt,it);
			RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],t08,t09,t18,t19,t28,t29,t38,t39,rt,it);

		/*...and now do 4 radix-5 transforms...*/
							 /*                                                inputs                                                   */ /*                 outputs                   */
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,rt,it);
		#else

		/*...Block 1:	*/
			t00=a[j1    ];	t01=a[j2    ];
			rt =a[j1+p01];	it =a[j2+p01];
			t10=t00-rt;  	t00=t00+rt;
			t11=t01-it;	t01=t01+it;

			t20=a[j1+p03];	t21=a[j2+p03];
			rt =a[j1+p02];	it =a[j2+p02];
			t30=t20-rt;  	t20=t20+rt;
			t31=t21-it;  	t21=t21+it;

			rt =t20;	t20=t00-rt;	t00=t00+rt;
			it =t21;	t21=t01-it;	t01=t01+it;

			rt =t30;	t30=t10-t31;	t10=t10+t31;
					t31=t11+rt;	t11=t11-rt;

		/*...Block 2:	*/
			t02=a[j1+p07];	t03=a[j2+p07];
			rt =a[j1+p06];	it =a[j2+p06];
			t12=t02-rt;  	t02=t02+rt;
			t13=t03-it;	t03=t03+it;

			t22=a[j1+p05];	t23=a[j2+p05];
			rt =a[j1+p04];	it =a[j2+p04];
			t32=t22-rt;  	t22=t22+rt;
			t33=t23-it;  	t23=t23+it;

			rt =t22;	t22=t02-rt;	t02=t02+rt;
			it =t23;	t23=t03-it;	t03=t03+it;

			rt =t32;	t32=t12-t33;	t12=t12+t33;
					t33=t13+rt;	t13=t13-rt;

		/*...Block 3:	*/
			t04=a[j1+p09];	t05=a[j2+p09];
			rt =a[j1+p08];	it =a[j2+p08];
			t14=t04-rt;  	t04=t04+rt;
			t15=t05-it;	t05=t05+it;

			t24=a[j1+p10];	t25=a[j2+p10];
			rt =a[j1+p11];	it =a[j2+p11];
			t34=t24-rt;  	t24=t24+rt;
			t35=t25-it;  	t25=t25+it;

			rt =t24;	t24=t04-rt;	t04=t04+rt;
			it =t25;	t25=t05-it;	t05=t05+it;

			rt =t34;	t34=t14-t35;	t14=t14+t35;
					t35=t15+rt;	t15=t15-rt;

		/*...Block 4:	*/
			t06=a[j1+p14];	t07=a[j2+p14];
			rt =a[j1+p15];	it =a[j2+p15];
			t16=t06-rt;  	t06=t06+rt;
			t17=t07-it;	t07=t07+it;

			t26=a[j1+p12];	t27=a[j2+p12];
			rt =a[j1+p13];	it =a[j2+p13];
			t36=t26-rt;  	t26=t26+rt;
			t37=t27-it;  	t27=t27+it;

			rt =t26;	t26=t06-rt;	t06=t06+rt;
			it =t27;	t27=t07-it;	t07=t07+it;

			rt =t36;	t36=t16-t37;	t16=t16+t37;
					t37=t17+rt;	t17=t17-rt;

		/*...Block 5:	*/
			t08=a[j1+p16];	t09=a[j2+p16];
			rt =a[j1+p17];	it =a[j2+p17];
			t18=t08-rt;  	t08=t08+rt;
			t19=t09-it;	t09=t09+it;

			t28=a[j1+p19];	t29=a[j2+p19];
			rt =a[j1+p18];	it =a[j2+p18];
			t38=t28-rt;  	t28=t28+rt;
			t39=t29-it;  	t29=t29+it;

			rt =t28;	t28=t08-rt;	t08=t08+rt;
			it =t29;	t29=t09-it;	t09=t09+it;

			rt =t38;	t38=t18-t39;	t18=t18+t39;
					t39=t19+rt;	t19=t19-rt;

	/*       ...and now do four radix-5 transforms.	*/
	/*...First radix-5 transform:	*/
			rt  =t08;			it  =t09;
			t08 =t02 -rt;		t09 =t03 -it;
			t02 =t02 +rt;		t03 =t03 +it;
			rt  =t06;			it  =t07;
			t06 =t04 -rt;		t07 =t05 -it;
			t04 =t04 +rt;		t05 =t05 +it;
			rt  = t02+t04;		it  = t03+t05;
			t00 = t00+rt;		t01 = t01+it;
			rt  = t00+uc1*rt;		it  = t01+uc1*it;
			t04 = uc2*(t02-t04);	t05 = uc2*(t03-t05);
			t02 = rt+t04;		t03 = it+t05;
			t04 = rt-t04;		t05 = it-t05;
			rt  = us1*(t06-t08);	it  = us1*(t07-t09);
			t06 = us2* t06;		t07 = us2* t07;
			t08 = us3* t08;		t09 = us3* t09;
			t06 = rt-t06;		t07 = it-t07;
			t08 = rt+t08;		t09 = it+t09;
			a1p00r=t00;		a1p00i=t01;
			a1p04r=t02-t07;		a1p04i=t03+t06;
			a1p08r=t04-t09;		a1p08i=t05+t08;
			a1p12r=t04+t09;		a1p12i=t05-t08;
			a1p16r=t02+t07;		a1p16i=t03-t06;

	/*...Second radix-5 transform:	*/
			rt  =t18;			it  =t19;
			t18 =t12 -rt;		t19 =t13 -it;
			t12 =t12 +rt;		t13 =t13 +it;
			rt  =t16;			it  =t17;
			t16 =t14 -rt;		t17 =t15 -it;
			t14 =t14 +rt;		t15 =t15 +it;
			rt  = t12+t14;		it  = t13+t15;
			t10 = t10+rt;		t11 = t11+it;
			rt  = t10+uc1*rt;		it  = t11+uc1*it;
			t14 = uc2*(t12-t14);	t15 = uc2*(t13-t15);
			t12 = rt+t14;		t13 = it+t15;
			t14 = rt-t14;		t15 = it-t15;
			rt  = us1*(t16-t18);	it  = us1*(t17-t19);
			t16 = us2* t16;		t17 = us2* t17;
			t18 = us3* t18;		t19 = us3* t19;
			t16 = rt-t16;		t17 = it-t17;
			t18 = rt+t18;		t19 = it+t19;
			a1p15r=t10;		a1p15i=t11;
			a1p19r=t12-t17;		a1p19i=t13+t16;
			a1p03r=t14-t19;		a1p03i=t15+t18;
			a1p07r=t14+t19;		a1p07i=t15-t18;
			a1p11r=t12+t17;		a1p11i=t13-t16;

	/*...Third radix-5 transform:	*/
			rt  =t28;			it  =t29;
			t28 =t22 -rt;		t29 =t23 -it;
			t22 =t22 +rt;		t23 =t23 +it;
			rt  =t26;			it  =t27;
			t26 =t24 -rt;		t27 =t25 -it;
			t24 =t24 +rt;		t25 =t25 +it;
			rt  = t22+t24;		it  = t23+t25;
			t20 = t20+rt;		t21 = t21+it;
			rt  = t20+uc1*rt;		it  = t21+uc1*it;
			t24 = uc2*(t22-t24);	t25 = uc2*(t23-t25);
			t22 = rt+t24;		t23 = it+t25;
			t24 = rt-t24;		t25 = it-t25;
			rt  = us1*(t26-t28);	it  = us1*(t27-t29);
			t26 = us2* t26;		t27 = us2* t27;
			t28 = us3* t28;		t29 = us3* t29;
			t26 = rt-t26;		t27 = it-t27;
			t28 = rt+t28;		t29 = it+t29;
			a1p10r=t20;		a1p10i=t21;
			a1p14r=t22-t27;		a1p14i=t23+t26;
			a1p18r=t24-t29;		a1p18i=t25+t28;
			a1p02r=t24+t29;		a1p02i=t25-t28;
			a1p06r=t22+t27;		a1p06i=t23-t26;

	/*...Fourth radix-5 transform:	*/
			rt  =t38;			it  =t39;
			t38 =t32 -rt;		t39 =t33 -it;
			t32 =t32 +rt;		t33 =t33 +it;
			rt  =t36;			it  =t37;
			t36 =t34 -rt;		t37 =t35 -it;
			t34 =t34 +rt;		t35 =t35 +it;
			rt  = t32+t34;		it  = t33+t35;
			t30 = t30+rt;		t31 = t31+it;
			rt  = t30+uc1*rt;		it  = t31+uc1*it;
			t34 = uc2*(t32-t34);	t35 = uc2*(t33-t35);
			t32 = rt+t34;		t33 = it+t35;
			t34 = rt-t34;		t35 = it-t35;
			rt  = us1*(t36-t38);	it  = us1*(t37-t39);
			t36 = us2* t36;		t37 = us2* t37;
			t38 = us3* t38;		t39 = us3* t39;
			t36 = rt-t36;		t37 = it-t37;
			t38 = rt+t38;		t39 = it+t39;
			a1p05r=t30;		a1p05i=t31;
			a1p09r=t32-t37;		a1p09i=t33+t36;
			a1p13r=t34-t39;		a1p13i=t35+t38;
			a1p17r=t34+t39;		a1p17i=t35-t38;
			a1p01r=t32+t37;		a1p01i=t33-t36;
#endif

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 20 separate blocks of the A-array, we need 20 separate carries.	*/

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
			cmplx_carry_norm_nocheck (a1p01r,a1p01i,cy01,bjmodn01,1 );
			cmplx_carry_norm_nocheck (a1p02r,a1p02i,cy02,bjmodn02,2 );
			cmplx_carry_norm_nocheck (a1p03r,a1p03i,cy03,bjmodn03,3 );
			cmplx_carry_norm_nocheck (a1p04r,a1p04i,cy04,bjmodn04,4 );
			cmplx_carry_norm_nocheck (a1p05r,a1p05i,cy05,bjmodn05,5 );
			cmplx_carry_norm_nocheck (a1p06r,a1p06i,cy06,bjmodn06,6 );
			cmplx_carry_norm_nocheck (a1p07r,a1p07i,cy07,bjmodn07,7 );
			cmplx_carry_norm_nocheck (a1p08r,a1p08i,cy08,bjmodn08,8 );
			cmplx_carry_norm_nocheck (a1p09r,a1p09i,cy09,bjmodn09,9 );
			cmplx_carry_norm_nocheck (a1p10r,a1p10i,cy10,bjmodn10,10);
			cmplx_carry_norm_nocheck (a1p11r,a1p11i,cy11,bjmodn11,11);
			cmplx_carry_norm_nocheck (a1p12r,a1p12i,cy12,bjmodn12,12);
			cmplx_carry_norm_nocheck (a1p13r,a1p13i,cy13,bjmodn13,13);
			cmplx_carry_norm_nocheck (a1p14r,a1p14i,cy14,bjmodn14,14);
			cmplx_carry_norm_nocheck (a1p15r,a1p15i,cy15,bjmodn15,15);
			cmplx_carry_norm_nocheck (a1p16r,a1p16i,cy16,bjmodn16,16);
			cmplx_carry_norm_nocheck (a1p17r,a1p17i,cy17,bjmodn17,17);
			cmplx_carry_norm_nocheck (a1p18r,a1p18i,cy18,bjmodn18,18);
			cmplx_carry_norm_nocheck (a1p19r,a1p19i,cy19,bjmodn19,19);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	/*...The radix-20 DIF pass is here:	*/
	#if PFETCH
		addr = &a[j1];
		prefetch_p_doubles(addr);
	#endif

	#if 1
		/*...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do 4 radix-5 transforms...*/
							 /*                                                inputs                                                   */ /*                 outputs                   */
		#if PFETCH																																/*[--y3-] [--y4-] <<<<< swap last 2 outputs to undo swap of these in macro */
			RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it,addr,addp,p01,p02);
			RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it,addr,addp,p03,p04);
			RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it,addr,addp,p05,p06);
			RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it,addr,addp,p07,p08);
		#else
			RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it);
			RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it);
			RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it);
			RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it);
		#endif
		/*...and now do 5 radix-4 transforms...*/
							 /*          inputs           */ /*                                      outputs                                      */
		#if PFETCH
			addp = addr+p09;
			prefetch_p_doubles(addp);

			RADIX_04_DIF_PFETCH(t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it,addr,addp,p10,p11);
			RADIX_04_DIF_PFETCH(t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it,addr,addp,p12,p13);
			RADIX_04_DIF_PFETCH(t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it,addr,addp,p14,p15);
			RADIX_04_DIF_PFETCH(t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it,addr,addp,p16,p17);
			RADIX_04_DIF_PFETCH(t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it,addr,addp,p18,p19);
		#else
			RADIX_04_DIF       (t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it);
			RADIX_04_DIF       (t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it);
			RADIX_04_DIF       (t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
			RADIX_04_DIF       (t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
			RADIX_04_DIF       (t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it);
		#endif

	#else

	/*...First radix-5 transform:	*/
		/* Swap 4<->16, 8<->12 below to get permuted version of j1 + p00,04,08,12,16: */
			t00=a1p00r;			t01=a1p00i;
			t02=a1p16r+a1p04r;		t03=a1p16i+a1p04i;
			t04=a1p12r+a1p08r;		t05=a1p12i+a1p08i;
			t06=a1p12r-a1p08r;		t07=a1p12i-a1p08i;
			t08=a1p16r-a1p04r;		t09=a1p16i-a1p04i;
			rt = t02+t04;			it = t03+t05;
			t00= t00+rt;			t01= t01+it;		/* y0	*/
			rt = t00+uc1*rt;		it = t01+uc1*it;
			t04= uc2*(t02-t04);		t05= uc2*(t03-t05);
			t02= rt+t04;			t03= it+t05;
			t04= rt-t04;			t05= it-t05;
		  #if PFETCH
			addp = addr+p01;
			prefetch_p_doubles(addp);
		  #endif
			rt = us1*(t08-t06);		it = us1*(t09-t07);
			t06= us2* t06;		t07= us2* t07;
			t08= us3* t08;		t09= us3* t09;
			t06= rt+t06;			t07= it+t07;
			t08= rt-t08;			t09= it-t09;
			rt = t06;			it = t07;
			t06= t02+it;			t07= t03-rt;		/* y4 - note swap with y3 below!	*/
			t02= t02-it;			t03= t03+rt;		/* y1	*/
			rt = t08;			it = t09;
			t08= t04+it;			t09= t05-rt;		/* y3	*/
			t04= t04-it;			t05= t05+rt;		/* y2	*/
		  #if PFETCH
			addp = addr+p02;
			prefetch_p_doubles(addp);
		  #endif

	/*...Second radix-5 transform:	*/
		/* Swap 01->15, 05->11, 09->07, 13->03, 17->19 below to get permuted version of j1 + p01,05,09,13,17: */
			t10=a1p15r;			t11=a1p15i;
			t12=a1p11r+a1p19r;		t13=a1p11i+a1p19i;
			t14=a1p07r+a1p03r;		t15=a1p07i+a1p03i;
			t16=a1p07r-a1p03r;		t17=a1p07i-a1p03i;
			t18=a1p11r-a1p19r;		t19=a1p11i-a1p19i;
			rt = t12+t14;			it = t13+t15;
			t10= t10+rt;			t11= t11+it;		/* y0	*/
			rt = t10+uc1*rt;		it = t11+uc1*it;
			t14= uc2*(t12-t14);		t15= uc2*(t13-t15);
			t12= rt+t14;			t13= it+t15;
			t14= rt-t14;			t15= it-t15;
		  #if PFETCH
			addp = addr+p03;
			prefetch_p_doubles(addp);
		  #endif
			rt = us1*(t18-t16);		it = us1*(t19-t17);
			t16= us2* t16;		t17= us2* t17;
			t18= us3* t18;		t19= us3* t19;
			t16= rt+t16;			t17= it+t17;
			t18= rt-t18;			t19= it-t19;
			rt = t16;			it = t17;
			t16= t12+it;			t17= t13-rt;		/* y4 - note swap with y3 below!	*/
			t12= t12-it;			t13= t13+rt;		/* y1	*/
			rt = t18;			it = t19;
			t18= t14+it;			t19= t15-rt;		/* y3	*/
			t14= t14-it;			t15= t15+rt;		/* y2	*/
		  #if PFETCH
			addp = addr+p04;
			prefetch_p_doubles(addp);
		  #endif

	/*...Third radix-5 transform:	*/
		/* Swap 02<->10, 14<->18 below to get permuted version of j1 + p02,06,10,14,18: */
			t20=a1p10r;			t21=a1p10i;
			t22=a1p06r+a1p14r;		t23=a1p06i+a1p14i;
			t24=a1p02r+a1p18r;		t25=a1p02i+a1p18i;
			t26=a1p02r-a1p18r;		t27=a1p02i-a1p18i;
			t28=a1p06r-a1p14r;		t29=a1p06i-a1p14i;
			rt = t22+t24;			it = t23+t25;
			t20= t20+rt;			t21= t21+it;		/* y0	*/
			rt = t20+uc1*rt;		it = t21+uc1*it;
			t24= uc2*(t22-t24);		t25= uc2*(t23-t25);
			t22= rt+t24;			t23= it+t25;
			t24= rt-t24;			t25= it-t25;
		  #if PFETCH
			addp = addr+p05;
			prefetch_p_doubles(addp);
		  #endif
			rt = us1*(t28-t26);		it = us1*(t29-t27);
			t26= us2* t26;		t27= us2* t27;
			t28= us3* t28;		t29= us3* t29;
			t26= rt+t26;			t27= it+t27;
			t28= rt-t28;			t29= it-t29;
			rt = t26;			it = t27;
			t26= t22+it;			t27= t23-rt;		/* y4 - note swap with y3 below!	*/
			t22= t22-it;			t23= t23+rt;		/* y1	*/
			rt = t28;			it = t29;
			t28= t24+it;			t29= t25-rt;		/* y3	*/
			t24= t24-it;			t25= t25+rt;		/* y2	*/
		  #if PFETCH
			addp = addr+p06;
			prefetch_p_doubles(addp);
		  #endif

	/*...Fourth radix-5 transform:	*/
		/* Swap 03->05, 07->01, 11->17, 15->13, 19->09 below to get permuted version of j1 + p03,07,11,15,19: */
			t30=a1p05r;			t31=a1p05i;
			t32=a1p01r+a1p09r;		t33=a1p01i+a1p09i;
			t34=a1p17r+a1p13r;		t35=a1p17i+a1p13i;
			t36=a1p17r-a1p13r;		t37=a1p17i-a1p13i;
			t38=a1p01r-a1p09r;		t39=a1p01i-a1p09i;
			rt = t32+t34;			it = t33+t35;
			t30= t30+rt;			t31= t31+it;		/* y0	*/
			rt = t30+uc1*rt;		it = t31+uc1*it;
			t34= uc2*(t32-t34);		t35= uc2*(t33-t35);
			t32= rt+t34;			t33= it+t35;
			t34= rt-t34;			t35= it-t35;
		  #if PFETCH
			addp = addr+p07;
			prefetch_p_doubles(addp);
		  #endif
			rt = us1*(t38-t36);		it = us1*(t39-t37);
			t36= us2* t36;		t37= us2* t37;
			t38= us3* t38;		t39= us3* t39;
			t36= rt+t36;			t37= it+t37;
			t38= rt-t38;			t39= it-t39;
			rt = t36;			it = t37;
			t36= t32+it;			t37= t33-rt;		/* y4 - note swap with y3 below!	*/
			t32= t32-it;			t33= t33+rt;		/* y1	*/
			rt = t38;			it = t39;
			t38= t34+it;			t39= t35-rt;		/* y3	*/
			t34= t34-it;			t35= t35+rt;		/* y2	*/
		  #if PFETCH
			addp = addr+p08;
			prefetch_p_doubles(addp);
		  #endif

	/*...and now do five radix-4 transforms:	*/
			rt =t20;	t20=t00-rt;	t00=t00+rt;
			it =t21;	t21=t01-it;	t01=t01+it;

			rt =t30;	t30=t10-rt;	t10=t10+rt;
			it =t31;	t31=t11-it;	t11=t11+it;
		#if PFETCH
			addp = addr+p09;
			prefetch_p_doubles(addp);
		#endif
			a[j1    ]=t00+t10;	a[j2    ]=t01+t11;
			a[j1+p01]=t00-t10;	a[j2+p01]=t01-t11;
			a[j1+p03]=t20-t31;	a[j2+p03]=t21+t30;	/* mpy by I is inlined here...	*/
			a[j1+p02]=t20+t31;	a[j2+p02]=t21-t30;
		#if PFETCH
			addp = addr+p10;
			prefetch_p_doubles(addp);
		#endif
	/*...Block 2: */
			rt =t22;	t22=t02-rt;	t02=t02+rt;
			it =t23;	t23=t03-it;	t03=t03+it;

			rt =t32;	t32=t12-rt;	t12=t12+rt;
			it =t33;	t33=t13-it;	t13=t13+it;
		#if PFETCH
			addp = addr+p11;
			prefetch_p_doubles(addp);
		#endif
			a[j1+p16]=t02+t12;	a[j2+p16]=t03+t13;
			a[j1+p17]=t02-t12;	a[j2+p17]=t03-t13;
			a[j1+p19]=t22-t33;	a[j2+p19]=t23+t32;	/* mpy by I is inlined here...	*/
			a[j1+p18]=t22+t33;	a[j2+p18]=t23-t32;
		#if PFETCH
			addp = addr+p12;
			prefetch_p_doubles(addp);
		#endif
	/*...Block 3: */
			rt =t24;	t24=t04-rt;	t04=t04+rt;
			it =t25;	t25=t05-it;	t05=t05+it;

			rt =t34;	t34=t14-rt;	t14=t14+rt;
			it =t35;	t35=t15-it;	t15=t15+it;
		#if PFETCH
			addp = addr+p13;
			prefetch_p_doubles(addp);
		#endif
			a[j1+p14]=t04+t14;	a[j2+p14]=t05+t15;
			a[j1+p15]=t04-t14;	a[j2+p15]=t05-t15;
			a[j1+p12]=t24-t35;	a[j2+p12]=t25+t34;	/* mpy by I is inlined here...	*/
			a[j1+p13]=t24+t35;	a[j2+p13]=t25-t34;
		#if PFETCH
			addp = addr+p14;
			prefetch_p_doubles(addp);
		#endif
	/*...Block 4: */
			rt =t26;	t26=t06-rt;	t06=t06+rt;
			it =t27;	t27=t07-it;	t07=t07+it;

			rt =t36;	t36=t16-rt;	t16=t16+rt;
			it =t37;	t37=t17-it;	t17=t17+it;
		#if PFETCH
			addp = addr+p15;
			prefetch_p_doubles(addp);
		#endif
			a[j1+p07]=t06+t16;	a[j2+p07]=t07+t17;
			a[j1+p06]=t06-t16;	a[j2+p06]=t07-t17;
			a[j1+p05]=t26-t37;	a[j2+p05]=t27+t36;	/* mpy by I is inlined here...	*/
			a[j1+p04]=t26+t37;	a[j2+p04]=t27-t36;
		#if PFETCH
			addp = addr+p16;
			prefetch_p_doubles(addp);
		#endif
	/*...Block 5: */
			rt =t28;	t28=t08-rt;	t08=t08+rt;
			it =t29;	t29=t09-it;	t09=t09+it;

			rt =t38;	t38=t18-rt;	t18=t18+rt;
			it =t39;	t39=t19-it;	t19=t19+it;
		#if PFETCH
			addp = addr+p17;
			prefetch_p_doubles(addp);
		#endif
			a[j1+p09]=t08+t18;	a[j2+p09]=t09+t19;
			a[j1+p08]=t08-t18;	a[j2+p08]=t09-t19;
			a[j1+p10]=t28-t39;	a[j2+p10]=t29+t38;	/* mpy by I is inlined here...	*/
			a[j1+p11]=t28+t39;	a[j2+p11]=t29-t38;
		#if PFETCH
			addp = addr+p18;
			prefetch_p_doubles(addp);
		#endif
		#if PFETCH
			addp = addr+p19;
			prefetch_p_doubles(addp);
		#endif

	#endif	/* endif(1) */

			iroot += root_incr;		/* increment sincos index.	*/

		}

		jstart += nwt;
		jhi    += nwt;
		col += 20;
		co3 -= 20;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-36 forward DIF FFT of the first block of 36 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 36 outputs of (1);
!   (3) Reweight and perform a radix-36 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 36 elements and repeat (1-4).
*/
	t00 = cy19;
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
	}
}

	if(fabs(cy00)+fabs(cy01)+fabs(cy02)+fabs(cy03)+fabs(cy04)+fabs(cy05)+fabs(cy06)+fabs(cy07)+fabs(cy08)+fabs(cy09)+fabs(cy10)+fabs(cy11)+fabs(cy12)+fabs(cy13)+fabs(cy14)+fabs(cy15)+fabs(cy16)+fabs(cy17)+fabs(cy18)+fabs(cy19) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix20_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

void radix20_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-20 complex DIF FFT pass on the data in the length-N real vector A.
!
!   Uses an optimized radix-5 transform a la Nussbaumer (2nd ed., p.146).
*/
	int j,j1,j2;
	static int n20, first_entry=TRUE
				,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19;
	static double uc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					uc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					us1 =  0.95105651629515357211,	/*  sin(u) */
					us2 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					us3 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39;

	if(!first_entry && (n/20) != n20)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n20=n/20;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n20;
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
	}

/*...The radix-20 pass is here.	*/

	for(j=0; j < n20; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
	/*
	Twiddleless version arranges 4 sets of radix-9 DFT inputs as follows: 0 in upper left corner, decrement 4 horizontally and 5 vertically:

		RADIX_05_DFT(00,16,12,08,04)
		RADIX_05_DFT(15,11,07,03,19)
		RADIX_05_DFT(10,06,02,18,14)
		RADIX_05_DFT(05,01,17,13,09)

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-4 DFT outputs.
	*/

	/*...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do 4 radix-5 transforms...*/
																																									/*[--y3-] [--y4-] <<<<< swap last 2 outputs to undo swap of these in macro */
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04],t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[j1+p15],a[j2+p15],a[j1+p11],a[j2+p11],a[j1+p07],a[j2+p07],a[j1+p03],a[j2+p03],a[j1+p19],a[j2+p19],t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[j1+p10],a[j2+p10],a[j1+p06],a[j2+p06],a[j1+p02],a[j2+p02],a[j1+p18],a[j2+p18],a[j1+p14],a[j2+p14],t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,a[j1+p05],a[j2+p05],a[j1+p01],a[j2+p01],a[j1+p17],a[j2+p17],a[j1+p13],a[j2+p13],a[j1+p09],a[j2+p09],t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it);

	/*...and now do 5 radix-4 transforms...*/

		RADIX_04_DIF(t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it);
		RADIX_04_DIF(t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it);
		RADIX_04_DIF(t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
		RADIX_04_DIF(t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
		RADIX_04_DIF(t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it);
	}
}

/***************/

void radix20_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-20 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix20_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int j,j1,j2;
	static int n20, first_entry=TRUE
				,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19;
	static double uc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					uc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					us1 =  0.95105651629515357211,	/*  sin(u) */
					us2 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					us3 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39;

	if(!first_entry && (n/20) != n20)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n20=n/20;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n20;
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
	}

/*...The radix-20 pass is here.	*/

	for(j=0; j < n20; j += 2)
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

		[0,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19]
		00,16,12,08,04,15,11,07,03,19,10,06,02,18,14,05,01,17,13,09].	(*)

		Remember, inputs to DIT are bit-reversed, so
		a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19] contain
		x[0,10, 5,15, 1,11, 6,16, 2,12, 7,17, 3,13, 8,18, 4,14, 9,19], which get swapped [using the permutation (*) on the index *values*] to
		x[0,10,15, 5,16, 6,11, 1,12, 2, 7,17, 8,18, 3,13, 4,14,19, 9], which means the a-indices get swapped as
		a[0, 1, 3, 2| 7, 6, 5, 4| 9, 8,10,11|14,15,12,13|16,17,19,18]. These are the 5 quartets going into the radix-4 DFTs.

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF] to properly permute the radix-4 DFT outputs.
	*/
	/*...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do 5 radix-4 transforms...*/

		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00,t01,t10,t11,t20,t21,t30,t31,rt,it);
		RADIX_04_DIT(a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t02,t03,t12,t13,t22,t23,t32,t33,rt,it);
		RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],t04,t05,t14,t15,t24,t25,t34,t35,rt,it);
		RADIX_04_DIT(a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],t06,t07,t16,t17,t26,t27,t36,t37,rt,it);
		RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],t08,t09,t18,t19,t28,t29,t38,t39,rt,it);

	/*...and now do 4 radix-5 transforms...*/

		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04],rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,a[j1+p15],a[j2+p15],a[j1+p11],a[j2+p11],a[j1+p07],a[j2+p07],a[j1+p03],a[j2+p03],a[j1+p19],a[j2+p19],rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,a[j1+p10],a[j2+p10],a[j1+p06],a[j2+p06],a[j1+p02],a[j2+p02],a[j1+p18],a[j2+p18],a[j1+p14],a[j2+p14],rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,a[j1+p05],a[j2+p05],a[j1+p01],a[j2+p01],a[j1+p17],a[j2+p17],a[j1+p13],a[j2+p13],a[j1+p09],a[j2+p09],rt,it);
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
	cy20_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 20;
		const double crnd = 3.0*0x4000000*0x2000000;
		int j,j1,j2,k;
		int l,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		uint32 p01,p04,p08,p16;
		double *add0, *add1, *add2, *add3;
		struct complex *cc1, *max_err, *sse2_rnd, *half_arr, *tmp, *r00,*r10,*r20,*r30,\
			*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r;
		struct complex *cy00,*cy02,*cy04,*cy06,*cy08,*cy10,*cy12,*cy14,*cy16,*cy18;
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
#if FFT_DEBUG
	int ithread = thread_arg->tid;	/* unique thread index (use for debug) */
	fprintf(dbg_file,"cy20_process_chunk: thread %d, NDIVR = %d, NWT = %d, &rn0,1 = %llx %llx\n"\
		, ithread, thread_arg->ndivr, thread_arg->nwt, (uint64)thread_arg->rn0, (uint64)thread_arg->rn1);
#endif
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
		p04 = p01<<2;
		p08 = p01<<3;
		p16 = p01<<4;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );

		r00	= thread_arg->r00;		tmp	= r00 + 0x28;
									s1p00r = tmp + 0x00;
									s1p01r = tmp + 0x02;		cc1     = tmp + 0x29;
									s1p02r = tmp + 0x04;
									s1p03r = tmp + 0x06;
									s1p04r = tmp + 0x08;
		r10	= r00 + 0x0a;			s1p05r = tmp + 0x0a;
									s1p06r = tmp + 0x0c;		cy00    = tmp + 0x2e;
									s1p07r = tmp + 0x0e;		cy02    = tmp + 0x2f;
									s1p08r = tmp + 0x10;		cy04    = tmp + 0x30;
									s1p09r = tmp + 0x12;		cy06    = tmp + 0x31;
		r20	= r00 + 0x14;			s1p10r = tmp + 0x14;		cy08    = tmp + 0x32;
									s1p11r = tmp + 0x16;		cy10    = tmp + 0x33;
									s1p12r = tmp + 0x18;		cy12    = tmp + 0x34;
									s1p13r = tmp + 0x1a;		cy14    = tmp + 0x35;
									s1p14r = tmp + 0x1c;		cy16    = tmp + 0x36;
		r30	= r00 + 0x1e;			s1p15r = tmp + 0x1e;		cy18    = tmp + 0x37;
									s1p16r = tmp + 0x20;		max_err = tmp + 0x38;
									s1p17r = tmp + 0x22;		sse2_rnd= tmp + 0x39;
									s1p18r = tmp + 0x24;		half_arr= tmp + 0x3a;	/* This table needs 20x16 bytes */
									s1p19r = tmp + 0x26;

		tmp = sse2_rnd;
		ASSERT(HERE, ((tmp + 0)->re == crnd && (tmp + 0)->im == crnd), "thread-local memcheck failed!");
		ASSERT(HERE, (tmp + 1+10)->re * (tmp + 1+14)->re == 1.0 && (tmp + 1+10)->im * (tmp + 1+14)->im == 1.0, "thread-local memcheck failed!");

		max_err->re = 0.0;	max_err->im = 0.0;

		sign_mask = (uint64*)(r00 + radix20_creals_in_local_store);
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
				SSE2_RADIX20_DIT_NOTWIDDLE(add0,p01,p04,r00,r10,r20,r30,cc1,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r);

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
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
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
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

				SSE2_RADIX20_DIF_NOTWIDDLE(add0,p01,p04,p08,p16,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r,cc1,r00,r10,r20,r30)

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

