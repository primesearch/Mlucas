/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

#ifdef USE_SSE2

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)	/* GCC-style inline ASM: */

		#if OS_BITS == 64

		  #define GCC_ASM_FULL_INLINE	1	// 0 to use small-macro form below, 1 to inline the fused macros as a single big blob of assembler (64-bit only)
		  #if GCC_ASM_FULL_INLINE
			#include "radix40_ditN_cy_dif1_gcc64.h"
		  #endif

		#endif

	#endif

#endif	/* USE_SSE2 */

/**************/

int radix40_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-40 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-40 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n40,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32;
	static double	uc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					uc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					us1 =  0.95105651629515357211,	/*  sin(u) */
					us2 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					us3 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	static double radix_inv, n2inv;
	double
		 t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15
		,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31
		,t32,t33,t34,t35,t36,t37,t38,t39
		,scale;
	double maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;

#ifdef USE_SSE2

	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */

  #if defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC)

	static struct complex *sc_arr = 0x0, *sc_ptr;
	static struct complex *isrt2, *cc1, *ss1, *cc2, *ss2, *ss3, *max_err, *sse2_rnd, *half_arr, *tmp
		,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09
		,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19
		,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29
		,*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39
		,*r40,*r41,*r42,*r43,*r44,*r45,*r46,*r47,*r48,*r49
		,*r50,*r51,*r52,*r53,*r54,*r55,*r56,*r57,*r58,*r59
		,*r60,*r61,*r62,*r63,*r64,*r65,*r66,*r67,*r68,*r69
		,*r70,*r71,*r72,*r73,*r74,*r75,*r76,*r77,*r78,*r79
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r
	,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i
	,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r
	,*s1p20i,*s1p21i,*s1p22i,*s1p23i,*s1p24i,*s1p25i,*s1p26i,*s1p27i,*s1p28i,*s1p29i,*s1p30i,*s1p31i,*s1p32i,*s1p33i,*s1p34i,*s1p35i,*s1p36i,*s1p37i,*s1p38i,*s1p39i;

	static uint64 *sm_arr = 0x0, *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19
	          ,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39;
	static struct complex *cy_r00,*cy_r02,*cy_r04,*cy_r06,*cy_r08,*cy_r10,*cy_r12,*cy_r14,*cy_r16,*cy_r18,*cy_r20,*cy_r22,*cy_r24,*cy_r26,*cy_r28,*cy_r30,*cy_r32,*cy_r34,*cy_r36,*cy_r38;

  #else

	#error SSE2 code not supported for this compiler!

  #endif

  #ifdef DEBUG_SSE2
	double
	 a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p36r,a1p37r,a1p38r,a1p39r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,a1p36i,a1p37i,a1p38i,a1p39i;
  #endif

#else

	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double t40,t41,t42,t43,t44,t45,t46,t47
		,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63
		,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
		,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71;
	double rt,it,temp,frac;
	double wt,wtinv,wtA,wtB,wtC;
	int m,m2;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39;
	double
	 a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p36r,a1p37r,a1p38r,a1p39r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,a1p36i,a1p37i,a1p38i,a1p39i
	,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_r36,cy_r37,cy_r38,cy_r39;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0;
	static int *_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0,*_bjmodn28 = 0x0,*_bjmodn29 = 0x0,*_bjmodn30 = 0x0,*_bjmodn31 = 0x0,*_bjmodn32 = 0x0,*_bjmodn33 = 0x0,*_bjmodn34 = 0x0,*_bjmodn35 = 0x0,*_bjmodn36 = 0x0,*_bjmodn37 = 0x0,*_bjmodn38 = 0x0,*_bjmodn39 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy_r00 = 0x0,*_cy_r01 = 0x0,*_cy_r02 = 0x0,*_cy_r03 = 0x0,*_cy_r04 = 0x0,*_cy_r05 = 0x0,*_cy_r06 = 0x0,*_cy_r07 = 0x0,*_cy_r08 = 0x0,*_cy_r09 = 0x0,*_cy_r10 = 0x0,*_cy_r11 = 0x0,*_cy_r12 = 0x0,*_cy_r13 = 0x0,*_cy_r14 = 0x0,*_cy_r15 = 0x0,*_cy_r16 = 0x0,*_cy_r17 = 0x0,*_cy_r18 = 0x0,*_cy_r19 = 0x0,
	*_cy_r20 = 0x0,*_cy_r21 = 0x0,*_cy_r22 = 0x0,*_cy_r23 = 0x0,*_cy_r24 = 0x0,*_cy_r25 = 0x0,*_cy_r26 = 0x0,*_cy_r27 = 0x0,*_cy_r28 = 0x0,*_cy_r29 = 0x0,*_cy_r30 = 0x0,*_cy_r31 = 0x0,*_cy_r32 = 0x0,*_cy_r33 = 0x0,*_cy_r34 = 0x0,*_cy_r35 = 0x0,*_cy_r36 = 0x0,*_cy_r37 = 0x0,*_cy_r38 = 0x0,*_cy_r39 = 0x0;

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "radix40_ditN_cy_dif1: Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

/*...change n40 and n_div_wt to non-static to work around a gcc compiler bug. */
	n40   = n/40;
	n_div_nwt = n40 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n40)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/40 in radix40_ditN_cy_dif1.\n",iter);
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

		ASSERT(HERE, CY_THREADS >= NTHREADS,"radix40_ditN_cy_dif1.c: CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"radix40_ditN_cy_dif1.c: CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, n40      %CY_THREADS == 0,"radix40_ditN_cy_dif1.c: n40      %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"radix40_ditN_cy_dif1.c: n_div_nwt%CY_THREADS != 0");
		}

	#ifdef MULTITHREAD
		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
	#endif

		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)40));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef USE_SSE2
		ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE), "SSE2 currently only supports Mersenne-mod!");
		sc_arr = ALLOC_COMPLEX(sc_arr, 212);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		/* Size here is [8 + radix/2 + 4] 8-byte elements */
		sm_arr = ALLOC_UINT64(sm_arr, 32);	if(!sm_arr){ sprintf(cbuf, "FATAL: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sm_ptr = ALIGN_UINT64(sm_arr);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 80 16-byte slots of sc_arr for temporaries, next 5 for the nontrivial complex 16th roots,
	next 80 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/								tmp	= sc_ptr + 0x50;		isrt2	= tmp + 0x50;
		r00	= sc_ptr + 0x00;		s1p00r = tmp + 0x00;		cc1		= tmp + 0x51;
		r01	= sc_ptr + 0x01;		s1p00i = tmp + 0x01;		cc2  	= tmp + 0x52;
		r02	= sc_ptr + 0x02;		s1p01r = tmp + 0x02;		ss1		= tmp + 0x53;
		r03	= sc_ptr + 0x03;		s1p01i = tmp + 0x03;		ss2		= tmp + 0x54;
		r04	= sc_ptr + 0x04;		s1p02r = tmp + 0x04;		ss3		= tmp + 0x55;
		r05	= sc_ptr + 0x05;		s1p02i = tmp + 0x05;		cy_r00	= tmp + 0x56;
		r06	= sc_ptr + 0x06;		s1p03r = tmp + 0x06;		cy_r02	= tmp + 0x57;
		r07	= sc_ptr + 0x07;		s1p03i = tmp + 0x07;		cy_r04	= tmp + 0x58;
		r08	= sc_ptr + 0x08;		s1p04r = tmp + 0x08;		cy_r06	= tmp + 0x59;
		r09	= sc_ptr + 0x09;		s1p04i = tmp + 0x09;		cy_r08	= tmp + 0x5a;
		r10	= sc_ptr + 0x0a;		s1p05r = tmp + 0x0a;		cy_r10	= tmp + 0x5b;
		r11	= sc_ptr + 0x0b;		s1p05i = tmp + 0x0b;		cy_r12	= tmp + 0x5c;
		r12	= sc_ptr + 0x0c;		s1p06r = tmp + 0x0c;		cy_r14	= tmp + 0x5d;
		r13	= sc_ptr + 0x0d;		s1p06i = tmp + 0x0d;		cy_r16	= tmp + 0x5e;
		r14	= sc_ptr + 0x0e;		s1p07r = tmp + 0x0e;		cy_r18	= tmp + 0x5f;
		r15	= sc_ptr + 0x0f;		s1p07i = tmp + 0x0f;		cy_r20	= tmp + 0x60;
		r16	= sc_ptr + 0x10;		s1p08r = tmp + 0x10;		cy_r22	= tmp + 0x61;
		r17	= sc_ptr + 0x11;		s1p08i = tmp + 0x11;		cy_r24	= tmp + 0x62;
		r18	= sc_ptr + 0x12;		s1p09r = tmp + 0x12;		cy_r26	= tmp + 0x63;
		r19	= sc_ptr + 0x13;		s1p09i = tmp + 0x13;		cy_r28	= tmp + 0x64;
		r20	= sc_ptr + 0x14;		s1p10r = tmp + 0x14;		cy_r30	= tmp + 0x65;
		r21	= sc_ptr + 0x15;		s1p10i = tmp + 0x15;		cy_r32	= tmp + 0x66;
		r22	= sc_ptr + 0x16;		s1p11r = tmp + 0x16;		cy_r34	= tmp + 0x67;
		r23	= sc_ptr + 0x17;		s1p11i = tmp + 0x17;		cy_r36	= tmp + 0x68;
		r24	= sc_ptr + 0x18;		s1p12r = tmp + 0x18;		cy_r38	= tmp + 0x69;
		r25	= sc_ptr + 0x19;		s1p12i = tmp + 0x19;		max_err = tmp + 0x6a;
		r26	= sc_ptr + 0x1a;		s1p13r = tmp + 0x1a;		sse2_rnd= tmp + 0x6b;
		r27	= sc_ptr + 0x1b;		s1p13i = tmp + 0x1b;		half_arr= tmp + 0x6c;	/* This table needs 20x16 bytes */
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
		r40	= sc_ptr + 0x28;		s1p20r = tmp + 0x28;
		r41	= sc_ptr + 0x29;		s1p20i = tmp + 0x29;
		r42	= sc_ptr + 0x2a;		s1p21r = tmp + 0x2a;
		r43	= sc_ptr + 0x2b;		s1p21i = tmp + 0x2b;
		r44	= sc_ptr + 0x2c;		s1p22r = tmp + 0x2c;
		r45	= sc_ptr + 0x2d;		s1p22i = tmp + 0x2d;
		r46	= sc_ptr + 0x2e;		s1p23r = tmp + 0x2e;
		r47	= sc_ptr + 0x2f;		s1p23i = tmp + 0x2f;
		r48	= sc_ptr + 0x30;		s1p24r = tmp + 0x30;
		r49	= sc_ptr + 0x31;		s1p24i = tmp + 0x31;
		r50	= sc_ptr + 0x32;		s1p25r = tmp + 0x32;
		r51	= sc_ptr + 0x33;		s1p25i = tmp + 0x33;
		r52	= sc_ptr + 0x34;		s1p26r = tmp + 0x34;
		r53	= sc_ptr + 0x35;		s1p26i = tmp + 0x35;
		r54	= sc_ptr + 0x36;		s1p27r = tmp + 0x36;
		r55	= sc_ptr + 0x37;		s1p27i = tmp + 0x37;
		r56	= sc_ptr + 0x38;		s1p28r = tmp + 0x38;
		r57	= sc_ptr + 0x39;		s1p28i = tmp + 0x39;
		r58	= sc_ptr + 0x3a;		s1p29r = tmp + 0x3a;
		r59	= sc_ptr + 0x3b;		s1p29i = tmp + 0x3b;
		r60	= sc_ptr + 0x3c;		s1p30r = tmp + 0x3c;
		r61	= sc_ptr + 0x3d;		s1p30i = tmp + 0x3d;
		r62	= sc_ptr + 0x3e;		s1p31r = tmp + 0x3e;
		r63	= sc_ptr + 0x3f;		s1p31i = tmp + 0x3f;
		r64	= sc_ptr + 0x40;		s1p32r = tmp + 0x40;
		r65	= sc_ptr + 0x41;		s1p32i = tmp + 0x41;
		r66	= sc_ptr + 0x42;		s1p33r = tmp + 0x42;
		r67	= sc_ptr + 0x43;		s1p33i = tmp + 0x43;
		r68	= sc_ptr + 0x44;		s1p34r = tmp + 0x44;
		r69	= sc_ptr + 0x45;		s1p34i = tmp + 0x45;
		r70	= sc_ptr + 0x46;		s1p35r = tmp + 0x46;
		r71	= sc_ptr + 0x47;		s1p35i = tmp + 0x47;
		r72	= sc_ptr + 0x48;		s1p36r = tmp + 0x48;
		r73	= sc_ptr + 0x49;		s1p36i = tmp + 0x49;
		r74	= sc_ptr + 0x4a;		s1p37r = tmp + 0x4a;
		r75	= sc_ptr + 0x4b;		s1p37i = tmp + 0x4b;
		r76	= sc_ptr + 0x4c;		s1p38r = tmp + 0x4c;
		r77	= sc_ptr + 0x4d;		s1p38i = tmp + 0x4d;
		r78	= sc_ptr + 0x4e;		s1p39r = tmp + 0x4e;
		r79	= sc_ptr + 0x4f;		s1p39i = tmp + 0x4f;

		/* These remain fixed: */
		isrt2->re = ISRT2;	isrt2->im = ISRT2;
		cc1->re = uc1;	cc1->im = uc1;
		cc2->re = uc2;	cc2->im = uc2;
		ss1->re = us1;	ss1->im = us1;
		ss2->re = us2;	ss2->im = us2;
		ss3->re = us3;	ss3->im = us3;

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		sse2_rnd->re = 3.0*0x4000000*0x2000000;
		sse2_rnd->im = 3.0*0x4000000*0x2000000;

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
		bjmodn36 = bjmodn35 + 1;
		bjmodn37 = bjmodn36 + 1;
		bjmodn38 = bjmodn37 + 1;
		bjmodn39 = bjmodn38 + 1;
	#endif

/*   constant index offsets for array load/stores are here.	*/

		p01 = n40;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );

		if(_cy_r00)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn00); _bjmodn00 = 0x0;			free((void *)_cy_r00); _cy_r00 = 0x0;
			free((void *)_bjmodn01); _bjmodn01 = 0x0;			free((void *)_cy_r01); _cy_r01 = 0x0;
			free((void *)_bjmodn02); _bjmodn02 = 0x0;			free((void *)_cy_r02); _cy_r02 = 0x0;
			free((void *)_bjmodn03); _bjmodn03 = 0x0;			free((void *)_cy_r03); _cy_r03 = 0x0;
			free((void *)_bjmodn04); _bjmodn04 = 0x0;			free((void *)_cy_r04); _cy_r04 = 0x0;
			free((void *)_bjmodn05); _bjmodn05 = 0x0;			free((void *)_cy_r05); _cy_r05 = 0x0;
			free((void *)_bjmodn06); _bjmodn06 = 0x0;			free((void *)_cy_r06); _cy_r06 = 0x0;
			free((void *)_bjmodn07); _bjmodn07 = 0x0;			free((void *)_cy_r07); _cy_r07 = 0x0;
			free((void *)_bjmodn08); _bjmodn08 = 0x0;			free((void *)_cy_r08); _cy_r08 = 0x0;
			free((void *)_bjmodn09); _bjmodn09 = 0x0;			free((void *)_cy_r09); _cy_r09 = 0x0;
			free((void *)_bjmodn10); _bjmodn10 = 0x0;			free((void *)_cy_r10); _cy_r10 = 0x0;
			free((void *)_bjmodn11); _bjmodn11 = 0x0;			free((void *)_cy_r11); _cy_r11 = 0x0;
			free((void *)_bjmodn12); _bjmodn12 = 0x0;			free((void *)_cy_r12); _cy_r12 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;			free((void *)_cy_r13); _cy_r13 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;			free((void *)_cy_r14); _cy_r14 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;			free((void *)_cy_r15); _cy_r15 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;			free((void *)_cy_r16); _cy_r16 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;			free((void *)_cy_r17); _cy_r17 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;			free((void *)_cy_r18); _cy_r18 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;			free((void *)_cy_r19); _cy_r19 = 0x0;
			free((void *)_bjmodn20); _bjmodn20 = 0x0;			free((void *)_cy_r20); _cy_r20 = 0x0;
			free((void *)_bjmodn21); _bjmodn21 = 0x0;			free((void *)_cy_r21); _cy_r21 = 0x0;
			free((void *)_bjmodn22); _bjmodn22 = 0x0;			free((void *)_cy_r22); _cy_r22 = 0x0;
			free((void *)_bjmodn23); _bjmodn23 = 0x0;			free((void *)_cy_r23); _cy_r23 = 0x0;
			free((void *)_bjmodn24); _bjmodn24 = 0x0;			free((void *)_cy_r24); _cy_r24 = 0x0;
			free((void *)_bjmodn25); _bjmodn25 = 0x0;			free((void *)_cy_r25); _cy_r25 = 0x0;
			free((void *)_bjmodn26); _bjmodn26 = 0x0;			free((void *)_cy_r26); _cy_r26 = 0x0;
			free((void *)_bjmodn27); _bjmodn27 = 0x0;			free((void *)_cy_r27); _cy_r27 = 0x0;
			free((void *)_bjmodn28); _bjmodn28 = 0x0;			free((void *)_cy_r28); _cy_r28 = 0x0;
			free((void *)_bjmodn29); _bjmodn29 = 0x0;			free((void *)_cy_r29); _cy_r29 = 0x0;
			free((void *)_bjmodn30); _bjmodn30 = 0x0;			free((void *)_cy_r30); _cy_r30 = 0x0;
			free((void *)_bjmodn31); _bjmodn31 = 0x0;			free((void *)_cy_r31); _cy_r31 = 0x0;
			free((void *)_bjmodn32); _bjmodn32 = 0x0;			free((void *)_cy_r32); _cy_r32 = 0x0;
			free((void *)_bjmodn33); _bjmodn33 = 0x0;			free((void *)_cy_r33); _cy_r33 = 0x0;
			free((void *)_bjmodn34); _bjmodn34 = 0x0;			free((void *)_cy_r34); _cy_r34 = 0x0;
			free((void *)_bjmodn35); _bjmodn35 = 0x0;			free((void *)_cy_r35); _cy_r35 = 0x0;
			free((void *)_bjmodn36); _bjmodn36 = 0x0;			free((void *)_cy_r36); _cy_r36 = 0x0;
			free((void *)_bjmodn37); _bjmodn37 = 0x0;			free((void *)_cy_r37); _cy_r37 = 0x0;
			free((void *)_bjmodn38); _bjmodn38 = 0x0;			free((void *)_cy_r38); _cy_r38 = 0x0;
			free((void *)_bjmodn39); _bjmodn39 = 0x0;			free((void *)_cy_r39); _cy_r39 = 0x0;

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		_i       	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_i== 0x0);
		_bjmodn00	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn00== 0x0);
		_bjmodn01	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn01== 0x0);
		_bjmodn02	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn02== 0x0);
		_bjmodn03	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn03== 0x0);
		_bjmodn04	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn04== 0x0);
		_bjmodn05	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn05== 0x0);
		_bjmodn06	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn06== 0x0);
		_bjmodn07	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn07== 0x0);
		_bjmodn08	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn08== 0x0);
		_bjmodn09	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn09== 0x0);
		_bjmodn10	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn10== 0x0);
		_bjmodn11	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn11== 0x0);
		_bjmodn12	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn12== 0x0);
		_bjmodn13	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn13== 0x0);
		_bjmodn14	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn14== 0x0);
		_bjmodn15	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn15== 0x0);
		_bjmodn16	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn16== 0x0);
		_bjmodn17	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn17== 0x0);
		_bjmodn18	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn18== 0x0);
		_bjmodn19	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn19== 0x0);
		_bjmodn20	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_bjmodn24	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_bjmodn28	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn28== 0x0);
		_bjmodn29	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn29== 0x0);
		_bjmodn30	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn30== 0x0);
		_bjmodn31	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn31== 0x0);
		_bjmodn32	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn32== 0x0);
		_bjmodn33	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn33== 0x0);
		_bjmodn34	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn34== 0x0);
		_bjmodn35	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn35== 0x0);
		_bjmodn36	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn36== 0x0);
		_bjmodn37	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn37== 0x0);
		_bjmodn38	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn38== 0x0);
		_bjmodn39	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn39== 0x0);
		_jstart  	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co3     == 0x0);

		_cy_r00	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r00== 0x0);
		_cy_r01	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r01== 0x0);
		_cy_r02	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r02== 0x0);
		_cy_r03	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r03== 0x0);
		_cy_r04	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r04== 0x0);
		_cy_r05	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r05== 0x0);
		_cy_r06	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r06== 0x0);
		_cy_r07	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r07== 0x0);
		_cy_r08	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r08== 0x0);
		_cy_r09	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r09== 0x0);
		_cy_r10	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r10== 0x0);
		_cy_r11	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r11== 0x0);
		_cy_r12	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r12== 0x0);
		_cy_r13	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r13== 0x0);
		_cy_r14	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r14== 0x0);
		_cy_r15	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r15== 0x0);
		_cy_r16	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r16== 0x0);
		_cy_r17	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r17== 0x0);
		_cy_r18	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r18== 0x0);
		_cy_r19	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r19== 0x0);
		_cy_r20	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r20== 0x0);
		_cy_r21	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r21== 0x0);
		_cy_r22	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r22== 0x0);
		_cy_r23	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r23== 0x0);
		_cy_r24	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r24== 0x0);
		_cy_r25	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r25== 0x0);
		_cy_r26	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r26== 0x0);
		_cy_r27	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r27== 0x0);
		_cy_r28	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r28== 0x0);
		_cy_r29	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r29== 0x0);
		_cy_r30	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r30== 0x0);
		_cy_r31	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r31== 0x0);
		_cy_r32	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r32== 0x0);
		_cy_r33	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r33== 0x0);
		_cy_r34	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r34== 0x0);
		_cy_r35	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r35== 0x0);
		_cy_r36	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r36== 0x0);
		_cy_r37	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r37== 0x0);
		_cy_r38	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r38== 0x0);
		_cy_r39	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r39== 0x0);

		_maxerr	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix40_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/40-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix40_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		jhi = n40/CY_THREADS;

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
		for(j=0; j < n40; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	}	/* endif(first_entry) */

/*...The radix-40 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r00[ithread] = 0;
		_cy_r01[ithread] = 0;
		_cy_r02[ithread] = 0;
		_cy_r03[ithread] = 0;
		_cy_r04[ithread] = 0;
		_cy_r05[ithread] = 0;
		_cy_r06[ithread] = 0;
		_cy_r07[ithread] = 0;
		_cy_r08[ithread] = 0;
		_cy_r09[ithread] = 0;
		_cy_r10[ithread] = 0;
		_cy_r11[ithread] = 0;
		_cy_r12[ithread] = 0;
		_cy_r13[ithread] = 0;
		_cy_r14[ithread] = 0;
		_cy_r15[ithread] = 0;
		_cy_r16[ithread] = 0;
		_cy_r17[ithread] = 0;
		_cy_r18[ithread] = 0;
		_cy_r19[ithread] = 0;
		_cy_r20[ithread] = 0;
		_cy_r21[ithread] = 0;
		_cy_r22[ithread] = 0;
		_cy_r23[ithread] = 0;
		_cy_r24[ithread] = 0;
		_cy_r25[ithread] = 0;
		_cy_r26[ithread] = 0;
		_cy_r27[ithread] = 0;
		_cy_r28[ithread] = 0;
		_cy_r29[ithread] = 0;
		_cy_r30[ithread] = 0;
		_cy_r31[ithread] = 0;
		_cy_r32[ithread] = 0;
		_cy_r33[ithread] = 0;
		_cy_r34[ithread] = 0;
		_cy_r35[ithread] = 0;
		_cy_r36[ithread] = 0;
		_cy_r37[ithread] = 0;
		_cy_r38[ithread] = 0;
		_cy_r39[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
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

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_jstart[ithread] = ithread*n40/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*40);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+40 - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-40;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn00[ithread] = _bjmodnini[ithread];
		_bjmodn01[ithread] = _bjmodn00[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn01[ithread] = _bjmodn01[ithread] + ( (-(int)((uint32)_bjmodn01[ithread] >> 31)) & n);
		_bjmodn02[ithread] = _bjmodn01[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn02[ithread] = _bjmodn02[ithread] + ( (-(int)((uint32)_bjmodn02[ithread] >> 31)) & n);
		_bjmodn03[ithread] = _bjmodn02[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn03[ithread] = _bjmodn03[ithread] + ( (-(int)((uint32)_bjmodn03[ithread] >> 31)) & n);
		_bjmodn04[ithread] = _bjmodn03[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn04[ithread] = _bjmodn04[ithread] + ( (-(int)((uint32)_bjmodn04[ithread] >> 31)) & n);
		_bjmodn05[ithread] = _bjmodn04[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn05[ithread] = _bjmodn05[ithread] + ( (-(int)((uint32)_bjmodn05[ithread] >> 31)) & n);
		_bjmodn06[ithread] = _bjmodn05[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn06[ithread] = _bjmodn06[ithread] + ( (-(int)((uint32)_bjmodn06[ithread] >> 31)) & n);
		_bjmodn07[ithread] = _bjmodn06[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn07[ithread] = _bjmodn07[ithread] + ( (-(int)((uint32)_bjmodn07[ithread] >> 31)) & n);
		_bjmodn08[ithread] = _bjmodn07[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn08[ithread] = _bjmodn08[ithread] + ( (-(int)((uint32)_bjmodn08[ithread] >> 31)) & n);
		_bjmodn09[ithread] = _bjmodn08[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn09[ithread] = _bjmodn09[ithread] + ( (-(int)((uint32)_bjmodn09[ithread] >> 31)) & n);
		_bjmodn10[ithread] = _bjmodn09[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn10[ithread] = _bjmodn10[ithread] + ( (-(int)((uint32)_bjmodn10[ithread] >> 31)) & n);
		_bjmodn11[ithread] = _bjmodn10[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn11[ithread] = _bjmodn11[ithread] + ( (-(int)((uint32)_bjmodn11[ithread] >> 31)) & n);
		_bjmodn12[ithread] = _bjmodn11[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn12[ithread] = _bjmodn12[ithread] + ( (-(int)((uint32)_bjmodn12[ithread] >> 31)) & n);
		_bjmodn13[ithread] = _bjmodn12[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn13[ithread] = _bjmodn13[ithread] + ( (-(int)((uint32)_bjmodn13[ithread] >> 31)) & n);
		_bjmodn14[ithread] = _bjmodn13[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn14[ithread] = _bjmodn14[ithread] + ( (-(int)((uint32)_bjmodn14[ithread] >> 31)) & n);
		_bjmodn15[ithread] = _bjmodn14[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn15[ithread] = _bjmodn15[ithread] + ( (-(int)((uint32)_bjmodn15[ithread] >> 31)) & n);
		_bjmodn16[ithread] = _bjmodn15[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn16[ithread] = _bjmodn16[ithread] + ( (-(int)((uint32)_bjmodn16[ithread] >> 31)) & n);
		_bjmodn17[ithread] = _bjmodn16[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn17[ithread] = _bjmodn17[ithread] + ( (-(int)((uint32)_bjmodn17[ithread] >> 31)) & n);
		_bjmodn18[ithread] = _bjmodn17[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn18[ithread] = _bjmodn18[ithread] + ( (-(int)((uint32)_bjmodn18[ithread] >> 31)) & n);
		_bjmodn19[ithread] = _bjmodn18[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn19[ithread] = _bjmodn19[ithread] + ( (-(int)((uint32)_bjmodn19[ithread] >> 31)) & n);
		_bjmodn20[ithread] = _bjmodn19[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn20[ithread] = _bjmodn20[ithread] + ( (-(int)((uint32)_bjmodn20[ithread] >> 31)) & n);
		_bjmodn21[ithread] = _bjmodn20[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn21[ithread] = _bjmodn21[ithread] + ( (-(int)((uint32)_bjmodn21[ithread] >> 31)) & n);
		_bjmodn22[ithread] = _bjmodn21[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn22[ithread] = _bjmodn22[ithread] + ( (-(int)((uint32)_bjmodn22[ithread] >> 31)) & n);
		_bjmodn23[ithread] = _bjmodn22[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn23[ithread] = _bjmodn23[ithread] + ( (-(int)((uint32)_bjmodn23[ithread] >> 31)) & n);
		_bjmodn24[ithread] = _bjmodn23[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn24[ithread] = _bjmodn24[ithread] + ( (-(int)((uint32)_bjmodn24[ithread] >> 31)) & n);
		_bjmodn25[ithread] = _bjmodn24[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn25[ithread] = _bjmodn25[ithread] + ( (-(int)((uint32)_bjmodn25[ithread] >> 31)) & n);
		_bjmodn26[ithread] = _bjmodn25[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn26[ithread] = _bjmodn26[ithread] + ( (-(int)((uint32)_bjmodn26[ithread] >> 31)) & n);
		_bjmodn27[ithread] = _bjmodn26[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn27[ithread] = _bjmodn27[ithread] + ( (-(int)((uint32)_bjmodn27[ithread] >> 31)) & n);
		_bjmodn28[ithread] = _bjmodn27[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn28[ithread] = _bjmodn28[ithread] + ( (-(int)((uint32)_bjmodn28[ithread] >> 31)) & n);
		_bjmodn29[ithread] = _bjmodn28[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn29[ithread] = _bjmodn29[ithread] + ( (-(int)((uint32)_bjmodn29[ithread] >> 31)) & n);
		_bjmodn30[ithread] = _bjmodn29[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn30[ithread] = _bjmodn30[ithread] + ( (-(int)((uint32)_bjmodn30[ithread] >> 31)) & n);
		_bjmodn31[ithread] = _bjmodn30[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn31[ithread] = _bjmodn31[ithread] + ( (-(int)((uint32)_bjmodn31[ithread] >> 31)) & n);
		_bjmodn32[ithread] = _bjmodn31[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn32[ithread] = _bjmodn32[ithread] + ( (-(int)((uint32)_bjmodn32[ithread] >> 31)) & n);
		_bjmodn33[ithread] = _bjmodn32[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn33[ithread] = _bjmodn33[ithread] + ( (-(int)((uint32)_bjmodn33[ithread] >> 31)) & n);
		_bjmodn34[ithread] = _bjmodn33[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn34[ithread] = _bjmodn34[ithread] + ( (-(int)((uint32)_bjmodn34[ithread] >> 31)) & n);
		_bjmodn35[ithread] = _bjmodn34[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn35[ithread] = _bjmodn35[ithread] + ( (-(int)((uint32)_bjmodn35[ithread] >> 31)) & n);
		_bjmodn36[ithread] = _bjmodn35[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn36[ithread] = _bjmodn36[ithread] + ( (-(int)((uint32)_bjmodn36[ithread] >> 31)) & n);
		_bjmodn37[ithread] = _bjmodn36[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn37[ithread] = _bjmodn37[ithread] + ( (-(int)((uint32)_bjmodn37[ithread] >> 31)) & n);
		_bjmodn38[ithread] = _bjmodn37[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn38[ithread] = _bjmodn38[ithread] + ( (-(int)((uint32)_bjmodn38[ithread] >> 31)) & n);
		_bjmodn39[ithread] = _bjmodn38[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn39[ithread] = _bjmodn39[ithread] + ( (-(int)((uint32)_bjmodn39[ithread] >> 31)) & n);
	}

    /* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
    if(!full_pass)
    {
		khi = 1;
	#if FFT_DEBUG
		fprintf(dbg_file, "radix40_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifdef MULTITHREAD
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,rt,it,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p36r,a1p37r,a1p38r,a1p39r,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,a1p36i,a1p37i,a1p38i,a1p39i,bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_r36,cy_r37,cy_r38,cy_r39) default(shared) schedule(static)
#endif

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
		*bjmodn00 = _bjmodn00[ithread];				cy_r00->re = _cy_r00[ithread];
		*bjmodn01 = _bjmodn01[ithread];				cy_r00->im = _cy_r01[ithread];
		*bjmodn02 = _bjmodn02[ithread];				cy_r02->re = _cy_r02[ithread];
		*bjmodn03 = _bjmodn03[ithread];				cy_r02->im = _cy_r03[ithread];
		*bjmodn04 = _bjmodn04[ithread];				cy_r04->re = _cy_r04[ithread];
		*bjmodn05 = _bjmodn05[ithread];				cy_r04->im = _cy_r05[ithread];
		*bjmodn06 = _bjmodn06[ithread];				cy_r06->re = _cy_r06[ithread];
		*bjmodn07 = _bjmodn07[ithread];				cy_r06->im = _cy_r07[ithread];
		*bjmodn08 = _bjmodn08[ithread];				cy_r08->re = _cy_r08[ithread];
		*bjmodn09 = _bjmodn09[ithread];				cy_r08->im = _cy_r09[ithread];
		*bjmodn10 = _bjmodn10[ithread];				cy_r10->re = _cy_r10[ithread];
		*bjmodn11 = _bjmodn11[ithread];				cy_r10->im = _cy_r11[ithread];
		*bjmodn12 = _bjmodn12[ithread];				cy_r12->re = _cy_r12[ithread];
		*bjmodn13 = _bjmodn13[ithread];				cy_r12->im = _cy_r13[ithread];
		*bjmodn14 = _bjmodn14[ithread];				cy_r14->re = _cy_r14[ithread];
		*bjmodn15 = _bjmodn15[ithread];				cy_r14->im = _cy_r15[ithread];
		*bjmodn16 = _bjmodn16[ithread];				cy_r16->re = _cy_r16[ithread];
		*bjmodn17 = _bjmodn17[ithread];				cy_r16->im = _cy_r17[ithread];
		*bjmodn18 = _bjmodn18[ithread];				cy_r18->re = _cy_r18[ithread];
		*bjmodn19 = _bjmodn19[ithread];				cy_r18->im = _cy_r19[ithread];
		*bjmodn20 = _bjmodn20[ithread];				cy_r20->re = _cy_r20[ithread];
		*bjmodn21 = _bjmodn21[ithread];				cy_r20->im = _cy_r21[ithread];
		*bjmodn22 = _bjmodn22[ithread];				cy_r22->re = _cy_r22[ithread];
		*bjmodn23 = _bjmodn23[ithread];				cy_r22->im = _cy_r23[ithread];
		*bjmodn24 = _bjmodn24[ithread];				cy_r24->re = _cy_r24[ithread];
		*bjmodn25 = _bjmodn25[ithread];				cy_r24->im = _cy_r25[ithread];
		*bjmodn26 = _bjmodn26[ithread];				cy_r26->re = _cy_r26[ithread];
		*bjmodn27 = _bjmodn27[ithread];				cy_r26->im = _cy_r27[ithread];
		*bjmodn28 = _bjmodn28[ithread];				cy_r28->re = _cy_r28[ithread];
		*bjmodn29 = _bjmodn29[ithread];				cy_r28->im = _cy_r29[ithread];
		*bjmodn30 = _bjmodn30[ithread];				cy_r30->re = _cy_r30[ithread];
		*bjmodn31 = _bjmodn31[ithread];				cy_r30->im = _cy_r31[ithread];
		*bjmodn32 = _bjmodn32[ithread];				cy_r32->re = _cy_r32[ithread];
		*bjmodn33 = _bjmodn33[ithread];				cy_r32->im = _cy_r33[ithread];
		*bjmodn34 = _bjmodn34[ithread];				cy_r34->re = _cy_r34[ithread];
		*bjmodn35 = _bjmodn35[ithread];				cy_r34->im = _cy_r35[ithread];
		*bjmodn36 = _bjmodn36[ithread];				cy_r36->re = _cy_r36[ithread];
		*bjmodn37 = _bjmodn37[ithread];				cy_r36->im = _cy_r37[ithread];
		*bjmodn38 = _bjmodn38[ithread];				cy_r38->re = _cy_r38[ithread];
		*bjmodn39 = _bjmodn39[ithread];				cy_r38->im = _cy_r39[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];				cy_r00 = _cy_r00[ithread];
		bjmodn01 = _bjmodn01[ithread];				cy_r01 = _cy_r01[ithread];
		bjmodn02 = _bjmodn02[ithread];				cy_r02 = _cy_r02[ithread];
		bjmodn03 = _bjmodn03[ithread];				cy_r03 = _cy_r03[ithread];
		bjmodn04 = _bjmodn04[ithread];				cy_r04 = _cy_r04[ithread];
		bjmodn05 = _bjmodn05[ithread];				cy_r05 = _cy_r05[ithread];
		bjmodn06 = _bjmodn06[ithread];				cy_r06 = _cy_r06[ithread];
		bjmodn07 = _bjmodn07[ithread];				cy_r07 = _cy_r07[ithread];
		bjmodn08 = _bjmodn08[ithread];				cy_r08 = _cy_r08[ithread];
		bjmodn09 = _bjmodn09[ithread];				cy_r09 = _cy_r09[ithread];
		bjmodn10 = _bjmodn10[ithread];				cy_r10 = _cy_r10[ithread];
		bjmodn11 = _bjmodn11[ithread];				cy_r11 = _cy_r11[ithread];
		bjmodn12 = _bjmodn12[ithread];				cy_r12 = _cy_r12[ithread];
		bjmodn13 = _bjmodn13[ithread];				cy_r13 = _cy_r13[ithread];
		bjmodn14 = _bjmodn14[ithread];				cy_r14 = _cy_r14[ithread];
		bjmodn15 = _bjmodn15[ithread];				cy_r15 = _cy_r15[ithread];
		bjmodn16 = _bjmodn16[ithread];				cy_r16 = _cy_r16[ithread];
		bjmodn17 = _bjmodn17[ithread];				cy_r17 = _cy_r17[ithread];
		bjmodn18 = _bjmodn18[ithread];				cy_r18 = _cy_r18[ithread];
		bjmodn19 = _bjmodn19[ithread];				cy_r19 = _cy_r19[ithread];
		bjmodn20 = _bjmodn20[ithread];				cy_r20 = _cy_r20[ithread];
		bjmodn21 = _bjmodn21[ithread];				cy_r21 = _cy_r21[ithread];
		bjmodn22 = _bjmodn22[ithread];				cy_r22 = _cy_r22[ithread];
		bjmodn23 = _bjmodn23[ithread];				cy_r23 = _cy_r23[ithread];
		bjmodn24 = _bjmodn24[ithread];				cy_r24 = _cy_r24[ithread];
		bjmodn25 = _bjmodn25[ithread];				cy_r25 = _cy_r25[ithread];
		bjmodn26 = _bjmodn26[ithread];				cy_r26 = _cy_r26[ithread];
		bjmodn27 = _bjmodn27[ithread];				cy_r27 = _cy_r27[ithread];
		bjmodn28 = _bjmodn28[ithread];				cy_r28 = _cy_r28[ithread];
		bjmodn29 = _bjmodn29[ithread];				cy_r29 = _cy_r29[ithread];
		bjmodn30 = _bjmodn30[ithread];				cy_r30 = _cy_r30[ithread];
		bjmodn31 = _bjmodn31[ithread];				cy_r31 = _cy_r31[ithread];
		bjmodn32 = _bjmodn32[ithread];				cy_r32 = _cy_r32[ithread];
		bjmodn33 = _bjmodn33[ithread];				cy_r33 = _cy_r33[ithread];
		bjmodn34 = _bjmodn34[ithread];				cy_r34 = _cy_r34[ithread];
		bjmodn35 = _bjmodn35[ithread];				cy_r35 = _cy_r35[ithread];
		bjmodn36 = _bjmodn36[ithread];				cy_r36 = _cy_r36[ithread];
		bjmodn37 = _bjmodn37[ithread];				cy_r37 = _cy_r37[ithread];
		bjmodn38 = _bjmodn38[ithread];				cy_r38 = _cy_r38[ithread];
		bjmodn39 = _bjmodn39[ithread];				cy_r39 = _cy_r39[ithread];
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
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1+p08;	jp = j2+p08;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1+p16;	jp = j2+p16;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1+p24;	jp = j2+p24;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1+p32;	jp = j2+p32;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p04] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p05] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p06] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p07] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			jt = j1;		jp = j2;
			fprintf(stderr, "radix40_ditN_cy_dif1 inputs:\n");
			fprintf(stderr, "A_in[00]: %20.2f, %20.2f\n",a[j1    ],a[j2    ]);
			fprintf(stderr, "A_in[01]: %20.2f, %20.2f\n",a[j1+p01],a[j2+p01]);
			fprintf(stderr, "A_in[02]: %20.2f, %20.2f\n",a[j1+p02],a[j2+p02]);
			fprintf(stderr, "A_in[03]: %20.2f, %20.2f\n",a[j1+p03],a[j2+p03]);
			fprintf(stderr, "A_in[04]: %20.2f, %20.2f\n",a[j1+p04],a[j2+p04]);
			fprintf(stderr, "A_in[05]: %20.2f, %20.2f\n",a[j1+p05],a[j2+p05]);
			fprintf(stderr, "A_in[06]: %20.2f, %20.2f\n",a[j1+p06],a[j2+p06]);
			fprintf(stderr, "A_in[07]: %20.2f, %20.2f\n",a[j1+p07],a[j2+p07]);	jt = j1+p08;	jp = j2+p08;
			fprintf(stderr, "A_in[08]: %20.2f, %20.2f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_in[09]: %20.2f, %20.2f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_in[10]: %20.2f, %20.2f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_in[11]: %20.2f, %20.2f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_in[12]: %20.2f, %20.2f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_in[13]: %20.2f, %20.2f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_in[14]: %20.2f, %20.2f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_in[15]: %20.2f, %20.2f\n",a[jt+p07],a[jp+p07]);	jt = j1+p16;	jp = j2+p16;
			fprintf(stderr, "A_in[16]: %20.2f, %20.2f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_in[17]: %20.2f, %20.2f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_in[18]: %20.2f, %20.2f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_in[19]: %20.2f, %20.2f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_in[20]: %20.2f, %20.2f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_in[21]: %20.2f, %20.2f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_in[22]: %20.2f, %20.2f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_in[23]: %20.2f, %20.2f\n",a[jt+p07],a[jp+p07]);	jt = j1+p24;	jp = j2+p24;
			fprintf(stderr, "A_in[24]: %20.2f, %20.2f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_in[25]: %20.2f, %20.2f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_in[26]: %20.2f, %20.2f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_in[27]: %20.2f, %20.2f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_in[28]: %20.2f, %20.2f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_in[29]: %20.2f, %20.2f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_in[30]: %20.2f, %20.2f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_in[31]: %20.2f, %20.2f\n",a[jt+p07],a[jp+p07]);	jt = j1+p32;	jp = j2+p32;
			fprintf(stderr, "A_in[32]: %20.2f, %20.2f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_in[33]: %20.2f, %20.2f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_in[34]: %20.2f, %20.2f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_in[35]: %20.2f, %20.2f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_in[36]: %20.2f, %20.2f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_in[37]: %20.2f, %20.2f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_in[38]: %20.2f, %20.2f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_in[39]: %20.2f, %20.2f\n",a[jt+p07],a[jp+p07]);
			fprintf(stderr, "\n");
		#endif	/* DEBUG_SSE2 */

		/*
		!...gather the needed data (40 64-bit complex, i.e. 80 64-bit reals) and do a radix-40 DIT transform...
		*/
	#ifdef USE_SSE2

			add0 = &a[j1    ];
		//	add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p03;	add3 = add0+p02;

		#if defined(COMPILER_TYPE_MSVC) || !GCC_ASM_FULL_INLINE

				add0 = &a[j1    ];
				add1 = add0+p01;
				add2 = add0+p03;
				add3 = add0+p02;
				add4 = add0+p07;
				add5 = add0+p06;
				add6 = add0+p05;
				add7 = add0+p04;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00, isrt2)

				add3 = &a[j1+p16];
				add0 = add3+p03;
				add1 = add3+p02;
				add2 = add3+p01;
				add4 = add3+p05;
				add5 = add3+p04;
				add6 = add3+p06;
				add7 = add3+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r16, isrt2)

				add1 = &a[j1+p32];
				add0 = add1+p01;
				add2 = add1+p02;
				add3 = add1+p03;
				add4 = add1+p06;
				add5 = add1+p07;
				add6 = add1+p04;
				add7 = add1+p05;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r32, isrt2)

				add6 = &a[j1+p08];
				add0 = add6+p06;
				add1 = add6+p07;
				add2 = add6+p04;
				add3 = add6+p05;
				add4 = add6+p02;
				add5 = add6+p03;
				add7 = add6+p01;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r48, isrt2)

				add4 = &a[j1+p24];
				add0 = add4+p04;
				add1 = add4+p05;
				add2 = add4+p07;
				add3 = add4+p06;
				add5 = add4+p01;
				add6 = add4+p03;
				add7 = add4+p02;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r64, isrt2)

		  #ifdef DEBUG_SSE2
			fprintf(stderr, "radix40 DIT08x5 outputs A:\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r00->re,r01->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r02->re,r03->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r04->re,r05->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r06->re,r07->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r08->re,r09->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r10->re,r11->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r12->re,r13->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r14->re,r15->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r16->re,r17->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r18->re,r19->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r20->re,r21->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r22->re,r23->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r24->re,r25->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r26->re,r27->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r28->re,r29->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r30->re,r31->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r32->re,r33->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r34->re,r35->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r36->re,r37->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r38->re,r39->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r40->re,r41->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r42->re,r43->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r44->re,r45->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r46->re,r47->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r48->re,r49->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r50->re,r51->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r52->re,r53->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r54->re,r55->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r56->re,r57->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r58->re,r59->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r60->re,r61->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r62->re,r63->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r64->re,r65->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r66->re,r67->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r68->re,r69->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r70->re,r71->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r72->re,r73->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r74->re,r75->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r76->re,r77->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r78->re,r79->re);
			fprintf(stderr, "\n");
		  #endif

				SSE2_RADIX_05_DFT_0TWIDDLE(r00,r16,r32,r48,r64,cc1,s1p00r,s1p16r,s1p32r,s1p08r,s1p24r)
				SSE2_RADIX_05_DFT_0TWIDDLE(r02,r18,r34,r50,r66,cc1,s1p25r,s1p01r,s1p17r,s1p33r,s1p09r)
				SSE2_RADIX_05_DFT_0TWIDDLE(r04,r20,r36,r52,r68,cc1,s1p10r,s1p26r,s1p02r,s1p18r,s1p34r)
				SSE2_RADIX_05_DFT_0TWIDDLE(r06,r22,r38,r54,r70,cc1,s1p35r,s1p11r,s1p27r,s1p03r,s1p19r)
				SSE2_RADIX_05_DFT_0TWIDDLE(r08,r24,r40,r56,r72,cc1,s1p20r,s1p36r,s1p12r,s1p28r,s1p04r)
				SSE2_RADIX_05_DFT_0TWIDDLE(r10,r26,r42,r58,r74,cc1,s1p05r,s1p21r,s1p37r,s1p13r,s1p29r)
				SSE2_RADIX_05_DFT_0TWIDDLE(r12,r28,r44,r60,r76,cc1,s1p30r,s1p06r,s1p22r,s1p38r,s1p14r)
				SSE2_RADIX_05_DFT_0TWIDDLE(r14,r30,r46,r62,r78,cc1,s1p15r,s1p31r,s1p07r,s1p23r,s1p39r)

		#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

			SSE2_RADIX40_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p16,p24,p32,r00,cc1,s1p00r,s1p08r,s1p16r,s1p24r,s1p32r);

		#endif

		  #ifdef DEBUG_SSE2
			fprintf(stderr, "radix40 DIT outputs:\n");
			fprintf(stderr, "s1p00 = %20.2f, %20.2f\n",s1p00r->re,s1p00i->re);
			fprintf(stderr, "s1p08 = %20.2f, %20.2f\n",s1p08r->re,s1p08i->re);
			fprintf(stderr, "s1p16 = %20.2f, %20.2f\n",s1p16r->re,s1p16i->re);
			fprintf(stderr, "s1p24 = %20.2f, %20.2f\n",s1p24r->re,s1p24i->re);
			fprintf(stderr, "s1p32 = %20.2f, %20.2f\n",s1p32r->re,s1p32i->re);		fprintf(stderr, "\n");
			fprintf(stderr, "s1p01 = %20.2f, %20.2f\n",s1p01r->re,s1p01i->re);
			fprintf(stderr, "s1p09 = %20.2f, %20.2f\n",s1p09r->re,s1p09i->re);
			fprintf(stderr, "s1p17 = %20.2f, %20.2f\n",s1p17r->re,s1p17i->re);
			fprintf(stderr, "s1p25 = %20.2f, %20.2f\n",s1p25r->re,s1p25i->re);
			fprintf(stderr, "s1p33 = %20.2f, %20.2f\n",s1p33r->re,s1p33i->re);		fprintf(stderr, "\n");
			fprintf(stderr, "s1p02 = %20.2f, %20.2f\n",s1p02r->re,s1p02i->re);
			fprintf(stderr, "s1p10 = %20.2f, %20.2f\n",s1p10r->re,s1p10i->re);
			fprintf(stderr, "s1p18 = %20.2f, %20.2f\n",s1p18r->re,s1p18i->re);
			fprintf(stderr, "s1p26 = %20.2f, %20.2f\n",s1p26r->re,s1p26i->re);
			fprintf(stderr, "s1p34 = %20.2f, %20.2f\n",s1p34r->re,s1p34i->re);		fprintf(stderr, "\n");
			fprintf(stderr, "s1p03 = %20.2f, %20.2f\n",s1p03r->re,s1p03i->re);
			fprintf(stderr, "s1p11 = %20.2f, %20.2f\n",s1p11r->re,s1p11i->re);
			fprintf(stderr, "s1p19 = %20.2f, %20.2f\n",s1p19r->re,s1p19i->re);
			fprintf(stderr, "s1p27 = %20.2f, %20.2f\n",s1p27r->re,s1p27i->re);
			fprintf(stderr, "s1p35 = %20.2f, %20.2f\n",s1p35r->re,s1p35i->re);		fprintf(stderr, "\n");
			fprintf(stderr, "s1p04 = %20.2f, %20.2f\n",s1p04r->re,s1p04i->re);
			fprintf(stderr, "s1p12 = %20.2f, %20.2f\n",s1p12r->re,s1p12i->re);
			fprintf(stderr, "s1p20 = %20.2f, %20.2f\n",s1p20r->re,s1p20i->re);
			fprintf(stderr, "s1p28 = %20.2f, %20.2f\n",s1p28r->re,s1p28i->re);
			fprintf(stderr, "s1p36 = %20.2f, %20.2f\n",s1p36r->re,s1p36i->re);		fprintf(stderr, "\n");
			fprintf(stderr, "s1p05 = %20.2f, %20.2f\n",s1p05r->re,s1p05i->re);
			fprintf(stderr, "s1p13 = %20.2f, %20.2f\n",s1p13r->re,s1p13i->re);
			fprintf(stderr, "s1p21 = %20.2f, %20.2f\n",s1p21r->re,s1p21i->re);
			fprintf(stderr, "s1p29 = %20.2f, %20.2f\n",s1p29r->re,s1p29i->re);
			fprintf(stderr, "s1p37 = %20.2f, %20.2f\n",s1p37r->re,s1p37i->re);		fprintf(stderr, "\n");
			fprintf(stderr, "s1p06 = %20.2f, %20.2f\n",s1p06r->re,s1p06i->re);
			fprintf(stderr, "s1p14 = %20.2f, %20.2f\n",s1p14r->re,s1p14i->re);
			fprintf(stderr, "s1p22 = %20.2f, %20.2f\n",s1p22r->re,s1p22i->re);
			fprintf(stderr, "s1p30 = %20.2f, %20.2f\n",s1p30r->re,s1p30i->re);
			fprintf(stderr, "s1p38 = %20.2f, %20.2f\n",s1p38r->re,s1p38i->re);		fprintf(stderr, "\n");
			fprintf(stderr, "s1p07 = %20.2f, %20.2f\n",s1p07r->re,s1p07i->re);
			fprintf(stderr, "s1p15 = %20.2f, %20.2f\n",s1p15r->re,s1p15i->re);
			fprintf(stderr, "s1p23 = %20.2f, %20.2f\n",s1p23r->re,s1p23i->re);
			fprintf(stderr, "s1p31 = %20.2f, %20.2f\n",s1p31r->re,s1p31i->re);
			fprintf(stderr, "s1p39 = %20.2f, %20.2f\n",s1p39r->re,s1p39i->re);
		  #endif

	#else

		/*...gather the needed data (40 64-bit complex) and do 5 radix-8 transforms,	*/
						 /*                                                                                 inputs                                                                    */ /*                 intermediates                             */ /*                        outputs                            */
			RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,rt,it);	jt = j1+p16; jp = j2+p16;
			RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,rt,it);	jt = j1+p32; jp = j2+p32;
			RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,rt,it);	jt = j1+p08; jp = j2+p08;
			RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,rt,it);	jt = j1+p24; jp = j2+p24;
			RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,rt,it);
		  #ifdef DEBUG_SSE2
			fprintf(stderr, "radix40 DIT08x5 outputs A:\n");
			fprintf(stderr, "%20.2f, %20.2f\n",t00,t01);
			fprintf(stderr, "%20.2f, %20.2f\n",t02,t03);
			fprintf(stderr, "%20.2f, %20.2f\n",t04,t05);
			fprintf(stderr, "%20.2f, %20.2f\n",t06,t07);
			fprintf(stderr, "%20.2f, %20.2f\n",t08,t09);
			fprintf(stderr, "%20.2f, %20.2f\n",t10,t11);
			fprintf(stderr, "%20.2f, %20.2f\n",t12,t13);
			fprintf(stderr, "%20.2f, %20.2f\n",t14,t15);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",t16,t17);
			fprintf(stderr, "%20.2f, %20.2f\n",t18,t19);
			fprintf(stderr, "%20.2f, %20.2f\n",t20,t21);
			fprintf(stderr, "%20.2f, %20.2f\n",t22,t23);
			fprintf(stderr, "%20.2f, %20.2f\n",t24,t25);
			fprintf(stderr, "%20.2f, %20.2f\n",t26,t27);
			fprintf(stderr, "%20.2f, %20.2f\n",t28,t29);
			fprintf(stderr, "%20.2f, %20.2f\n",t30,t31);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",t32,t33);
			fprintf(stderr, "%20.2f, %20.2f\n",t34,t35);
			fprintf(stderr, "%20.2f, %20.2f\n",t36,t37);
			fprintf(stderr, "%20.2f, %20.2f\n",t38,t39);
			fprintf(stderr, "%20.2f, %20.2f\n",t40,t41);
			fprintf(stderr, "%20.2f, %20.2f\n",t42,t43);
			fprintf(stderr, "%20.2f, %20.2f\n",t44,t45);
			fprintf(stderr, "%20.2f, %20.2f\n",t46,t47);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",t48,t49);
			fprintf(stderr, "%20.2f, %20.2f\n",t50,t51);
			fprintf(stderr, "%20.2f, %20.2f\n",t52,t53);
			fprintf(stderr, "%20.2f, %20.2f\n",t54,t55);
			fprintf(stderr, "%20.2f, %20.2f\n",t56,t57);
			fprintf(stderr, "%20.2f, %20.2f\n",t58,t59);
			fprintf(stderr, "%20.2f, %20.2f\n",t60,t61);
			fprintf(stderr, "%20.2f, %20.2f\n",t62,t63);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",t64,t65);
			fprintf(stderr, "%20.2f, %20.2f\n",t66,t67);
			fprintf(stderr, "%20.2f, %20.2f\n",t68,t69);
			fprintf(stderr, "%20.2f, %20.2f\n",t70,t71);
			fprintf(stderr, "%20.2f, %20.2f\n",t72,t73);
			fprintf(stderr, "%20.2f, %20.2f\n",t74,t75);
			fprintf(stderr, "%20.2f, %20.2f\n",t76,t77);
			fprintf(stderr, "%20.2f, %20.2f\n",t78,t79);
			fprintf(stderr, "\n");
		  #endif

		/*...and now do 8 radix-5 transforms.	*/
						 /*   sincos     */ /*             inputs                 */ /*                       outputs                                   */
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t00,t01,t16,t17,t32,t33,t48,t49,t64,t65,a1p00r,a1p00i,a1p16r,a1p16i,a1p32r,a1p32i,a1p08r,a1p08i,a1p24r,a1p24i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t02,t03,t18,t19,t34,t35,t50,t51,t66,t67,a1p25r,a1p25i,a1p01r,a1p01i,a1p17r,a1p17i,a1p33r,a1p33i,a1p09r,a1p09i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t04,t05,t20,t21,t36,t37,t52,t53,t68,t69,a1p10r,a1p10i,a1p26r,a1p26i,a1p02r,a1p02i,a1p18r,a1p18i,a1p34r,a1p34i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t06,t07,t22,t23,t38,t39,t54,t55,t70,t71,a1p35r,a1p35i,a1p11r,a1p11i,a1p27r,a1p27i,a1p03r,a1p03i,a1p19r,a1p19i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t08,t09,t24,t25,t40,t41,t56,t57,t72,t73,a1p20r,a1p20i,a1p36r,a1p36i,a1p12r,a1p12i,a1p28r,a1p28i,a1p04r,a1p04i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t10,t11,t26,t27,t42,t43,t58,t59,t74,t75,a1p05r,a1p05i,a1p21r,a1p21i,a1p37r,a1p37i,a1p13r,a1p13i,a1p29r,a1p29i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t12,t13,t28,t29,t44,t45,t60,t61,t76,t77,a1p30r,a1p30i,a1p06r,a1p06i,a1p22r,a1p22i,a1p38r,a1p38i,a1p14r,a1p14i,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,t14,t15,t30,t31,t46,t47,t62,t63,t78,t79,a1p15r,a1p15i,a1p31r,a1p31i,a1p07r,a1p07i,a1p23r,a1p23i,a1p39r,a1p39i,rt,it)
		  #ifdef DEBUG_SSE2
			fprintf(stderr, "radix40 DIT05x8 DIT outputs A:\n");
			fprintf(stderr, "s1p00 = %20.2f, %20.2f\n",a1p00r,a1p00i);
			fprintf(stderr, "s1p08 = %20.2f, %20.2f\n",a1p08r,a1p08i);
			fprintf(stderr, "s1p16 = %20.2f, %20.2f\n",a1p16r,a1p16i);
			fprintf(stderr, "s1p24 = %20.2f, %20.2f\n",a1p24r,a1p24i);
			fprintf(stderr, "s1p32 = %20.2f, %20.2f\n",a1p32r,a1p32i);			fprintf(stderr, "\n");
			fprintf(stderr, "s1p01 = %20.2f, %20.2f\n",a1p01r,a1p01i);
			fprintf(stderr, "s1p09 = %20.2f, %20.2f\n",a1p09r,a1p09i);
			fprintf(stderr, "s1p17 = %20.2f, %20.2f\n",a1p17r,a1p17i);
			fprintf(stderr, "s1p25 = %20.2f, %20.2f\n",a1p25r,a1p25i);
			fprintf(stderr, "s1p33 = %20.2f, %20.2f\n",a1p33r,a1p33i);			fprintf(stderr, "\n");
			fprintf(stderr, "s1p02 = %20.2f, %20.2f\n",a1p02r,a1p02i);
			fprintf(stderr, "s1p10 = %20.2f, %20.2f\n",a1p10r,a1p10i);
			fprintf(stderr, "s1p18 = %20.2f, %20.2f\n",a1p18r,a1p18i);
			fprintf(stderr, "s1p26 = %20.2f, %20.2f\n",a1p26r,a1p26i);
			fprintf(stderr, "s1p34 = %20.2f, %20.2f\n",a1p34r,a1p34i);			fprintf(stderr, "\n");
			fprintf(stderr, "s1p03 = %20.2f, %20.2f\n",a1p03r,a1p03i);
			fprintf(stderr, "s1p11 = %20.2f, %20.2f\n",a1p11r,a1p11i);
			fprintf(stderr, "s1p19 = %20.2f, %20.2f\n",a1p19r,a1p19i);
			fprintf(stderr, "s1p27 = %20.2f, %20.2f\n",a1p27r,a1p27i);
			fprintf(stderr, "s1p35 = %20.2f, %20.2f\n",a1p35r,a1p35i);			fprintf(stderr, "\n");
			fprintf(stderr, "s1p04 = %20.2f, %20.2f\n",a1p04r,a1p04i);
			fprintf(stderr, "s1p12 = %20.2f, %20.2f\n",a1p12r,a1p12i);
			fprintf(stderr, "s1p20 = %20.2f, %20.2f\n",a1p20r,a1p20i);
			fprintf(stderr, "s1p28 = %20.2f, %20.2f\n",a1p28r,a1p28i);
			fprintf(stderr, "s1p36 = %20.2f, %20.2f\n",a1p36r,a1p36i);			fprintf(stderr, "\n");
			fprintf(stderr, "s1p05 = %20.2f, %20.2f\n",a1p05r,a1p05i);
			fprintf(stderr, "s1p13 = %20.2f, %20.2f\n",a1p13r,a1p13i);
			fprintf(stderr, "s1p21 = %20.2f, %20.2f\n",a1p21r,a1p21i);
			fprintf(stderr, "s1p29 = %20.2f, %20.2f\n",a1p29r,a1p29i);
			fprintf(stderr, "s1p37 = %20.2f, %20.2f\n",a1p37r,a1p37i);			fprintf(stderr, "\n");
			fprintf(stderr, "s1p06 = %20.2f, %20.2f\n",a1p06r,a1p06i);
			fprintf(stderr, "s1p14 = %20.2f, %20.2f\n",a1p14r,a1p14i);
			fprintf(stderr, "s1p22 = %20.2f, %20.2f\n",a1p22r,a1p22i);
			fprintf(stderr, "s1p30 = %20.2f, %20.2f\n",a1p30r,a1p30i);
			fprintf(stderr, "s1p38 = %20.2f, %20.2f\n",a1p38r,a1p38i);			fprintf(stderr, "\n");
			fprintf(stderr, "s1p07 = %20.2f, %20.2f\n",a1p07r,a1p07i);
			fprintf(stderr, "s1p15 = %20.2f, %20.2f\n",a1p15r,a1p15i);
			fprintf(stderr, "s1p23 = %20.2f, %20.2f\n",a1p23r,a1p23i);
			fprintf(stderr, "s1p31 = %20.2f, %20.2f\n",a1p31r,a1p31i);
			fprintf(stderr, "s1p39 = %20.2f, %20.2f\n",a1p39r,a1p39i);
		  #endif

	#endif

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 40 separate blocks of the A-array, we need 40 separate carries.	*/

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

			SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy_r36,cy_r38,bjmodn36);

		#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

			SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

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

			SSE2_cmplx_carry_norm_nocheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p36r,add1,add2,     cy_r36,cy_r38,bjmodn36);

		  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

			SSE2_cmplx_carry_norm_nocheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			SSE2_cmplx_carry_norm_nocheck2_2B(s1p36r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

		  #endif

			i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

		  #ifdef DEBUG_SSE2
			a1p00r = s1p00r->re;	a1p00i = s1p00i->re;	a1p01r = s1p01r->re;	a1p01i = s1p01i->re;	a1p02r = s1p02r->re;	a1p02i = s1p02i->re;	a1p03r = s1p03r->re;	a1p03i = s1p03i->re;	a1p04r = s1p04r->re;	a1p04i = s1p04i->re;	a1p05r = s1p05r->re;	a1p05i = s1p05i->re;	a1p06r = s1p06r->re;	a1p06i = s1p06i->re;	a1p07r = s1p07r->re;	a1p07i = s1p07i->re;	a1p08r = s1p08r->re;	a1p08i = s1p08i->re;	a1p09r = s1p09r->re;	a1p09i = s1p09i->re;	a1p10r = s1p10r->re;	a1p10i = s1p10i->re;	a1p11r = s1p11r->re;	a1p11i = s1p11i->re;	a1p12r = s1p12r->re;	a1p12i = s1p12i->re;	a1p13r = s1p13r->re;	a1p13i = s1p13i->re;	a1p14r = s1p14r->re;	a1p14i = s1p14i->re;	a1p15r = s1p15r->re;	a1p15i = s1p15i->re;	a1p16r = s1p16r->re;	a1p16i = s1p16i->re;	a1p17r = s1p17r->re;	a1p17i = s1p17i->re;	a1p18r = s1p18r->re;	a1p18i = s1p18i->re;	a1p19r = s1p19r->re;	a1p19i = s1p19i->re;	a1p20r = s1p20r->re;	a1p20i = s1p20i->re;	a1p21r = s1p21r->re;	a1p21i = s1p21i->re;	a1p22r = s1p22r->re;	a1p22i = s1p22i->re;	a1p23r = s1p23r->re;	a1p23i = s1p23i->re;	a1p24r = s1p24r->re;	a1p24i = s1p24i->re;	a1p25r = s1p25r->re;	a1p25i = s1p25i->re;	a1p26r = s1p26r->re;	a1p26i = s1p26i->re;	a1p27r = s1p27r->re;	a1p27i = s1p27i->re;	a1p28r = s1p28r->re;	a1p28i = s1p28i->re;	a1p29r = s1p29r->re;	a1p29i = s1p29i->re;	a1p30r = s1p30r->re;	a1p30i = s1p30i->re;	a1p31r = s1p31r->re;	a1p31i = s1p31i->re;	a1p32r = s1p32r->re;	a1p32i = s1p32i->re;	a1p33r = s1p33r->re;	a1p33i = s1p33i->re;	a1p34r = s1p34r->re;	a1p34i = s1p34i->re;	a1p35r = s1p35r->re;	a1p35i = s1p35i->re;	a1p36r = s1p36r->re;	a1p36i = s1p36i->re;	a1p37r = s1p37r->re;	a1p37i = s1p37i->re;	a1p38r = s1p38r->re;	a1p38i = s1p38i->re;	a1p39r = s1p39r->re;	a1p39i = s1p39i->re;
		  #endif

		#else	/* #ifdef USE_SSE2 */

/*...set0 is slightly different from others:	*/
			cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
			cmplx_carry_norm_errcheck (a1p01r,a1p01i,cy_r01,bjmodn01,1 );
			cmplx_carry_norm_errcheck (a1p02r,a1p02i,cy_r02,bjmodn02,2 );
			cmplx_carry_norm_errcheck (a1p03r,a1p03i,cy_r03,bjmodn03,3 );
			cmplx_carry_norm_errcheck (a1p04r,a1p04i,cy_r04,bjmodn04,4 );
			cmplx_carry_norm_errcheck (a1p05r,a1p05i,cy_r05,bjmodn05,5 );
			cmplx_carry_norm_errcheck (a1p06r,a1p06i,cy_r06,bjmodn06,6 );
			cmplx_carry_norm_errcheck (a1p07r,a1p07i,cy_r07,bjmodn07,7 );
			cmplx_carry_norm_errcheck (a1p08r,a1p08i,cy_r08,bjmodn08,8 );
			cmplx_carry_norm_errcheck (a1p09r,a1p09i,cy_r09,bjmodn09,9 );
			cmplx_carry_norm_errcheck (a1p10r,a1p10i,cy_r10,bjmodn10,10);
			cmplx_carry_norm_errcheck (a1p11r,a1p11i,cy_r11,bjmodn11,11);
			cmplx_carry_norm_errcheck (a1p12r,a1p12i,cy_r12,bjmodn12,12);
			cmplx_carry_norm_errcheck (a1p13r,a1p13i,cy_r13,bjmodn13,13);
			cmplx_carry_norm_errcheck (a1p14r,a1p14i,cy_r14,bjmodn14,14);
			cmplx_carry_norm_errcheck (a1p15r,a1p15i,cy_r15,bjmodn15,15);
			cmplx_carry_norm_errcheck (a1p16r,a1p16i,cy_r16,bjmodn16,16);
			cmplx_carry_norm_errcheck (a1p17r,a1p17i,cy_r17,bjmodn17,17);
			cmplx_carry_norm_errcheck (a1p18r,a1p18i,cy_r18,bjmodn18,18);
			cmplx_carry_norm_errcheck (a1p19r,a1p19i,cy_r19,bjmodn19,19);
			cmplx_carry_norm_errcheck (a1p20r,a1p20i,cy_r20,bjmodn20,20);
			cmplx_carry_norm_errcheck (a1p21r,a1p21i,cy_r21,bjmodn21,21);
			cmplx_carry_norm_errcheck (a1p22r,a1p22i,cy_r22,bjmodn22,22);
			cmplx_carry_norm_errcheck (a1p23r,a1p23i,cy_r23,bjmodn23,23);
			cmplx_carry_norm_errcheck (a1p24r,a1p24i,cy_r24,bjmodn24,24);
			cmplx_carry_norm_errcheck (a1p25r,a1p25i,cy_r25,bjmodn25,25);
			cmplx_carry_norm_errcheck (a1p26r,a1p26i,cy_r26,bjmodn26,26);
			cmplx_carry_norm_errcheck (a1p27r,a1p27i,cy_r27,bjmodn27,27);
			cmplx_carry_norm_errcheck (a1p28r,a1p28i,cy_r28,bjmodn28,28);
			cmplx_carry_norm_errcheck (a1p29r,a1p29i,cy_r29,bjmodn29,29);
			cmplx_carry_norm_errcheck (a1p30r,a1p30i,cy_r30,bjmodn30,30);
			cmplx_carry_norm_errcheck (a1p31r,a1p31i,cy_r31,bjmodn31,31);
			cmplx_carry_norm_errcheck (a1p32r,a1p32i,cy_r32,bjmodn32,32);
			cmplx_carry_norm_errcheck (a1p33r,a1p33i,cy_r33,bjmodn33,33);
			cmplx_carry_norm_errcheck (a1p34r,a1p34i,cy_r34,bjmodn34,34);
			cmplx_carry_norm_errcheck (a1p35r,a1p35i,cy_r35,bjmodn35,35);
			cmplx_carry_norm_errcheck (a1p36r,a1p36i,cy_r36,bjmodn36,36);
			cmplx_carry_norm_errcheck (a1p37r,a1p37i,cy_r37,bjmodn37,37);
			cmplx_carry_norm_errcheck (a1p38r,a1p38i,cy_r38,bjmodn38,38);
			cmplx_carry_norm_errcheck (a1p39r,a1p39i,cy_r39,bjmodn39,39);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	/* #ifdef USE_SSE2 */

		#ifdef DEBUG_SSE2
			fprintf(stderr, "radix40 DIF 40 inputs:\n");
			fprintf(stderr, "s1p00 = %20.2f, %20.2f\n",a1p00r,a1p00i);
			fprintf(stderr, "s1p01 = %20.2f, %20.2f\n",a1p01r,a1p01i);
			fprintf(stderr, "s1p02 = %20.2f, %20.2f\n",a1p02r,a1p02i);
			fprintf(stderr, "s1p03 = %20.2f, %20.2f\n",a1p03r,a1p03i);
			fprintf(stderr, "s1p04 = %20.2f, %20.2f\n",a1p04r,a1p04i);
			fprintf(stderr, "s1p05 = %20.2f, %20.2f\n",a1p05r,a1p05i);
			fprintf(stderr, "s1p06 = %20.2f, %20.2f\n",a1p06r,a1p06i);
			fprintf(stderr, "s1p07 = %20.2f, %20.2f\n",a1p07r,a1p07i);
			fprintf(stderr, "s1p08 = %20.2f, %20.2f\n",a1p08r,a1p08i);
			fprintf(stderr, "s1p09 = %20.2f, %20.2f\n",a1p09r,a1p09i);
			fprintf(stderr, "s1p10 = %20.2f, %20.2f\n",a1p10r,a1p10i);
			fprintf(stderr, "s1p11 = %20.2f, %20.2f\n",a1p11r,a1p11i);
			fprintf(stderr, "s1p12 = %20.2f, %20.2f\n",a1p12r,a1p12i);
			fprintf(stderr, "s1p13 = %20.2f, %20.2f\n",a1p13r,a1p13i);
			fprintf(stderr, "s1p14 = %20.2f, %20.2f\n",a1p14r,a1p14i);
			fprintf(stderr, "s1p15 = %20.2f, %20.2f\n",a1p15r,a1p15i);
			fprintf(stderr, "s1p16 = %20.2f, %20.2f\n",a1p16r,a1p16i);
			fprintf(stderr, "s1p17 = %20.2f, %20.2f\n",a1p17r,a1p17i);
			fprintf(stderr, "s1p18 = %20.2f, %20.2f\n",a1p18r,a1p18i);
			fprintf(stderr, "s1p19 = %20.2f, %20.2f\n",a1p19r,a1p19i);
			fprintf(stderr, "s1p20 = %20.2f, %20.2f\n",a1p20r,a1p20i);
			fprintf(stderr, "s1p21 = %20.2f, %20.2f\n",a1p21r,a1p21i);
			fprintf(stderr, "s1p22 = %20.2f, %20.2f\n",a1p22r,a1p22i);
			fprintf(stderr, "s1p23 = %20.2f, %20.2f\n",a1p23r,a1p23i);
			fprintf(stderr, "s1p24 = %20.2f, %20.2f\n",a1p24r,a1p24i);
			fprintf(stderr, "s1p25 = %20.2f, %20.2f\n",a1p25r,a1p25i);
			fprintf(stderr, "s1p26 = %20.2f, %20.2f\n",a1p26r,a1p26i);
			fprintf(stderr, "s1p27 = %20.2f, %20.2f\n",a1p27r,a1p27i);
			fprintf(stderr, "s1p28 = %20.2f, %20.2f\n",a1p28r,a1p28i);
			fprintf(stderr, "s1p29 = %20.2f, %20.2f\n",a1p29r,a1p29i);
			fprintf(stderr, "s1p30 = %20.2f, %20.2f\n",a1p30r,a1p30i);
			fprintf(stderr, "s1p31 = %20.2f, %20.2f\n",a1p31r,a1p31i);
			fprintf(stderr, "s1p32 = %20.2f, %20.2f\n",a1p32r,a1p32i);
			fprintf(stderr, "s1p33 = %20.2f, %20.2f\n",a1p33r,a1p33i);
			fprintf(stderr, "s1p34 = %20.2f, %20.2f\n",a1p34r,a1p34i);
			fprintf(stderr, "s1p35 = %20.2f, %20.2f\n",a1p35r,a1p35i);
			fprintf(stderr, "s1p36 = %20.2f, %20.2f\n",a1p36r,a1p36i);
			fprintf(stderr, "s1p37 = %20.2f, %20.2f\n",a1p37r,a1p37i);
			fprintf(stderr, "s1p38 = %20.2f, %20.2f\n",a1p38r,a1p38i);
			fprintf(stderr, "s1p39 = %20.2f, %20.2f\n",a1p39r,a1p39i);
		#endif

		#ifndef USE_SSE2

		/*...gather the needed data (40 64-bit complex) and do 8 radix-5 transforms...*/
						 /*   sincos     */ /*                              inputs                              */ /*          outputs                  */
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p32r,a1p32i,a1p24r,a1p24i,a1p16r,a1p16i,a1p08r,a1p08i,t00,t01,t16,t17,t32,t33,t48,t49,t64,t65,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,a1p35r,a1p35i,a1p27r,a1p27i,a1p19r,a1p19i,a1p11r,a1p11i,a1p03r,a1p03i,t02,t03,t18,t19,t34,t35,t50,t51,t66,t67,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,a1p30r,a1p30i,a1p22r,a1p22i,a1p14r,a1p14i,a1p06r,a1p06i,a1p38r,a1p38i,t04,t05,t20,t21,t36,t37,t52,t53,t68,t69,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,a1p25r,a1p25i,a1p17r,a1p17i,a1p09r,a1p09i,a1p01r,a1p01i,a1p33r,a1p33i,t06,t07,t22,t23,t38,t39,t54,t55,t70,t71,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,a1p20r,a1p20i,a1p12r,a1p12i,a1p04r,a1p04i,a1p36r,a1p36i,a1p28r,a1p28i,t08,t09,t24,t25,t40,t41,t56,t57,t72,t73,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p07r,a1p07i,a1p39r,a1p39i,a1p31r,a1p31i,a1p23r,a1p23i,t10,t11,t26,t27,t42,t43,t58,t59,t74,t75,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p02r,a1p02i,a1p34r,a1p34i,a1p26r,a1p26i,a1p18r,a1p18i,t12,t13,t28,t29,t44,t45,t60,t61,t76,t77,rt,it);
			RADIX_05_DFT(uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p37r,a1p37i,a1p29r,a1p29i,a1p21r,a1p21i,a1p13r,a1p13i,t14,t15,t30,t31,t46,t47,t62,t63,t78,t79,rt,it)
		  #ifdef DEBUG_SSE2
			fprintf(stderr, "radix40 DIF05x8 outputs A:\n");
			fprintf(stderr, "%20.2f, %20.2f\n",t00,t01);
			fprintf(stderr, "%20.2f, %20.2f\n",t02,t03);
			fprintf(stderr, "%20.2f, %20.2f\n",t04,t05);
			fprintf(stderr, "%20.2f, %20.2f\n",t06,t07);
			fprintf(stderr, "%20.2f, %20.2f\n",t08,t09);
			fprintf(stderr, "%20.2f, %20.2f\n",t10,t11);
			fprintf(stderr, "%20.2f, %20.2f\n",t12,t13);
			fprintf(stderr, "%20.2f, %20.2f\n",t14,t15);
			fprintf(stderr, "%20.2f, %20.2f\n",t16,t17);
			fprintf(stderr, "%20.2f, %20.2f\n",t18,t19);
			fprintf(stderr, "%20.2f, %20.2f\n",t20,t21);
			fprintf(stderr, "%20.2f, %20.2f\n",t22,t23);
			fprintf(stderr, "%20.2f, %20.2f\n",t24,t25);
			fprintf(stderr, "%20.2f, %20.2f\n",t26,t27);
			fprintf(stderr, "%20.2f, %20.2f\n",t28,t29);
			fprintf(stderr, "%20.2f, %20.2f\n",t30,t31);
			fprintf(stderr, "%20.2f, %20.2f\n",t32,t33);
			fprintf(stderr, "%20.2f, %20.2f\n",t34,t35);
			fprintf(stderr, "%20.2f, %20.2f\n",t36,t37);
			fprintf(stderr, "%20.2f, %20.2f\n",t38,t39);
			fprintf(stderr, "%20.2f, %20.2f\n",t40,t41);
			fprintf(stderr, "%20.2f, %20.2f\n",t42,t43);
			fprintf(stderr, "%20.2f, %20.2f\n",t44,t45);
			fprintf(stderr, "%20.2f, %20.2f\n",t46,t47);
			fprintf(stderr, "%20.2f, %20.2f\n",t48,t49);
			fprintf(stderr, "%20.2f, %20.2f\n",t50,t51);
			fprintf(stderr, "%20.2f, %20.2f\n",t52,t53);
			fprintf(stderr, "%20.2f, %20.2f\n",t54,t55);
			fprintf(stderr, "%20.2f, %20.2f\n",t56,t57);
			fprintf(stderr, "%20.2f, %20.2f\n",t58,t59);
			fprintf(stderr, "%20.2f, %20.2f\n",t60,t61);
			fprintf(stderr, "%20.2f, %20.2f\n",t62,t63);
			fprintf(stderr, "%20.2f, %20.2f\n",t64,t65);
			fprintf(stderr, "%20.2f, %20.2f\n",t66,t67);
			fprintf(stderr, "%20.2f, %20.2f\n",t68,t69);
			fprintf(stderr, "%20.2f, %20.2f\n",t70,t71);
			fprintf(stderr, "%20.2f, %20.2f\n",t72,t73);
			fprintf(stderr, "%20.2f, %20.2f\n",t74,t75);
			fprintf(stderr, "%20.2f, %20.2f\n",t76,t77);
			fprintf(stderr, "%20.2f, %20.2f\n",t78,t79);
			fprintf(stderr, "\n");
		  #endif

		/*...and now do 5 radix-8 transforms, swapping the t[48+i] <--> t[64+i] pairs to undo the last-2-outputs-swap in the RADIX_05_DFT macro:	*/
						 /*                            inputs                         */ /*                      intermediates                        */ /*                 outputs                   */
			RADIX_08_DIF(t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],rt,it);	jt = j1+p32; jp = j2+p32;
			RADIX_08_DIF(t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],rt,it);	jt = j1+p24; jp = j2+p24;
			RADIX_08_DIF(t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p08; jp = j2+p08;
			RADIX_08_DIF(t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p16; jp = j2+p16;
			RADIX_08_DIF(t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],rt,it);

		#else

		  #if defined(COMPILER_TYPE_MSVC) || !GCC_ASM_FULL_INLINE

		/*...gather the needed data (40 64-bit complex) and do 8 radix-5 transforms...*/
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p00r,s1p32r,s1p24r,s1p16r,s1p08r,cc1,r00,r16,r32,r48,r64)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p35r,s1p27r,s1p19r,s1p11r,s1p03r,cc1,r02,r18,r34,r50,r66)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p30r,s1p22r,s1p14r,s1p06r,s1p38r,cc1,r04,r20,r36,r52,r68)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p25r,s1p17r,s1p09r,s1p01r,s1p33r,cc1,r06,r22,r38,r54,r70)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p20r,s1p12r,s1p04r,s1p36r,s1p28r,cc1,r08,r24,r40,r56,r72)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p15r,s1p07r,s1p39r,s1p31r,s1p23r,cc1,r10,r26,r42,r58,r74)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p10r,s1p02r,s1p34r,s1p26r,s1p18r,cc1,r12,r28,r44,r60,r76)
			SSE2_RADIX_05_DFT_0TWIDDLE(s1p05r,s1p37r,s1p29r,s1p21r,s1p13r,cc1,r14,r30,r46,r62,r78)

		/*...and now do 5 radix-8 transforms, swapping the t[48+i] <--> t[64+i] pairs to undo the last-2-outputs-swap in the RADIX_05_DFT macro:	*/

			add0 = &a[j1    ];
			add1 = add0+p01;
			add2 = add0+p03;
			add3 = add0+p02;
			add4 = add0+p06;
			add5 = add0+p07;
			add6 = add0+p04;
			add7 = add0+p05;
			SSE2_RADIX8_DIF_0TWIDDLE(r00,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);

			add1 = &a[j1+p32];
			add0 = add1+p01;
			add2 = add1+p02;
			add3 = add1+p03;
			add4 = add1+p07;
			add5 = add1+p06;
			add6 = add1+p05;
			add7 = add1+p04;
			SSE2_RADIX8_DIF_0TWIDDLE(r16,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);

			add5 = &a[j1+p24];
			add0 = add5+p04;
			add1 = add5+p05;
			add2 = add5+p07;
			add3 = add5+p06;
			add4 = add5+p01;
			add6 = add5+p02;
			add7 = add5+p03;
			SSE2_RADIX8_DIF_0TWIDDLE(r32,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);

			add3 = &a[j1+p16];
			add0 = add3+p03;
			add1 = add3+p02;
			add2 = add3+p01;
			add4 = add3+p04;
			add5 = add3+p05;
			add6 = add3+p07;
			add7 = add3+p06;
			SSE2_RADIX8_DIF_0TWIDDLE(r48,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);

			add7 = &a[j1+p08];
			add0 = add7+p06;
			add1 = add7+p07;
			add2 = add7+p04;
			add3 = add7+p05;
			add4 = add7+p03;
			add5 = add7+p02;
			add6 = add7+p01;
			SSE2_RADIX8_DIF_0TWIDDLE(r64,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);

		  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

			add0 = &a[j1    ];
			SSE2_RADIX40_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32,r00,cc1,s1p00r,s1p08r,s1p16r,s1p24r,s1p32r);

		  #endif

		  #if 0//DEBUG_SSE2
			fprintf(stderr, "radix40 DIF05x8 outputs B:\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r00->re,r01->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r16->re,r17->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r32->re,r33->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r48->re,r49->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r64->re,r65->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r02->re,r03->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r18->re,r19->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r34->re,r35->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r50->re,r51->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r66->re,r67->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r04->re,r05->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r20->re,r21->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r36->re,r37->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r52->re,r53->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r68->re,r69->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r06->re,r07->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r22->re,r23->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r38->re,r39->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r54->re,r55->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r70->re,r71->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r08->re,r09->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r24->re,r25->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r40->re,r41->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r56->re,r57->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r72->re,r73->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r10->re,r11->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r26->re,r27->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r42->re,r43->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r58->re,r59->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r74->re,r75->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r12->re,r13->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r28->re,r29->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r44->re,r45->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r60->re,r61->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r76->re,r77->re);
			fprintf(stderr, "\n");
			fprintf(stderr, "%20.2f, %20.2f\n",r14->re,r15->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r30->re,r31->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r46->re,r47->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r62->re,r63->re);
			fprintf(stderr, "%20.2f, %20.2f\n",r78->re,r79->re);
			fprintf(stderr, "\n");
		  #endif

	#endif	/* #ifdef USE_SSE2 */

		#ifdef DEBUG_SSE2
			jt = j1;		jp = j2;
			fprintf(stderr, "radix40_ditN_cy_dif1_outputs:\n");
			fprintf(stderr, "A_out[00]: %20.2f, %20.2f\n",a[j1    ],a[j2    ]);
			fprintf(stderr, "A_out[01]: %20.2f, %20.2f\n",a[j1+p01],a[j2+p01]);
			fprintf(stderr, "A_out[02]: %20.2f, %20.2f\n",a[j1+p02],a[j2+p02]);
			fprintf(stderr, "A_out[03]: %20.2f, %20.2f\n",a[j1+p03],a[j2+p03]);
			fprintf(stderr, "A_out[04]: %20.2f, %20.2f\n",a[j1+p04],a[j2+p04]);
			fprintf(stderr, "A_out[05]: %20.2f, %20.2f\n",a[j1+p05],a[j2+p05]);
			fprintf(stderr, "A_out[06]: %20.2f, %20.2f\n",a[j1+p06],a[j2+p06]);
			fprintf(stderr, "A_out[07]: %20.2f, %20.2f\n",a[j1+p07],a[j2+p07]);	jt = j1+p08;	jp = j2+p08;
			fprintf(stderr, "A_out[08]: %20.2f, %20.2f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_out[09]: %20.2f, %20.2f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_out[10]: %20.2f, %20.2f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_out[11]: %20.2f, %20.2f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_out[12]: %20.2f, %20.2f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_out[13]: %20.2f, %20.2f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_out[14]: %20.2f, %20.2f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_out[15]: %20.2f, %20.2f\n",a[jt+p07],a[jp+p07]);	jt = j1+p16;	jp = j2+p16;
			fprintf(stderr, "A_out[16]: %20.2f, %20.2f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_out[17]: %20.2f, %20.2f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_out[18]: %20.2f, %20.2f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_out[19]: %20.2f, %20.2f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_out[20]: %20.2f, %20.2f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_out[21]: %20.2f, %20.2f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_out[22]: %20.2f, %20.2f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_out[23]: %20.2f, %20.2f\n",a[jt+p07],a[jp+p07]);	jt = j1+p24;	jp = j2+p24;
			fprintf(stderr, "A_out[24]: %20.2f, %20.2f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_out[25]: %20.2f, %20.2f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_out[26]: %20.2f, %20.2f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_out[27]: %20.2f, %20.2f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_out[28]: %20.2f, %20.2f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_out[29]: %20.2f, %20.2f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_out[30]: %20.2f, %20.2f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_out[31]: %20.2f, %20.2f\n",a[jt+p07],a[jp+p07]);	jt = j1+p32;	jp = j2+p32;
			fprintf(stderr, "A_out[32]: %20.2f, %20.2f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "A_out[33]: %20.2f, %20.2f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "A_out[34]: %20.2f, %20.2f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "A_out[35]: %20.2f, %20.2f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "A_out[36]: %20.2f, %20.2f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "A_out[37]: %20.2f, %20.2f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "A_out[38]: %20.2f, %20.2f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "A_out[39]: %20.2f, %20.2f\n",a[jt+p07],a[jp+p07]);
			exit(0);
		#endif	/* DEBUG_SSE2 */

			}

			jstart += nwt;
			jhi    += nwt;
			col += 40;
			co3 -= 40;

		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_SSE2
		_cy_r00[ithread] = cy_r00->re;		_cy_r20[ithread] = cy_r20->re;
		_cy_r01[ithread] = cy_r00->im;		_cy_r21[ithread] = cy_r20->im;
		_cy_r02[ithread] = cy_r02->re;		_cy_r22[ithread] = cy_r22->re;
		_cy_r03[ithread] = cy_r02->im;		_cy_r23[ithread] = cy_r22->im;
		_cy_r04[ithread] = cy_r04->re;		_cy_r24[ithread] = cy_r24->re;
		_cy_r05[ithread] = cy_r04->im;		_cy_r25[ithread] = cy_r24->im;
		_cy_r06[ithread] = cy_r06->re;		_cy_r26[ithread] = cy_r26->re;
		_cy_r07[ithread] = cy_r06->im;		_cy_r27[ithread] = cy_r26->im;
		_cy_r08[ithread] = cy_r08->re;		_cy_r28[ithread] = cy_r28->re;
		_cy_r09[ithread] = cy_r08->im;		_cy_r29[ithread] = cy_r28->im;
		_cy_r10[ithread] = cy_r10->re;		_cy_r30[ithread] = cy_r30->re;
		_cy_r11[ithread] = cy_r10->im;		_cy_r31[ithread] = cy_r30->im;
		_cy_r12[ithread] = cy_r12->re;		_cy_r32[ithread] = cy_r32->re;
		_cy_r13[ithread] = cy_r12->im;		_cy_r33[ithread] = cy_r32->im;
		_cy_r14[ithread] = cy_r14->re;		_cy_r34[ithread] = cy_r34->re;
		_cy_r15[ithread] = cy_r14->im;		_cy_r35[ithread] = cy_r34->im;
		_cy_r16[ithread] = cy_r16->re;		_cy_r36[ithread] = cy_r36->re;
		_cy_r17[ithread] = cy_r16->im;		_cy_r37[ithread] = cy_r36->im;
		_cy_r18[ithread] = cy_r18->re;		_cy_r38[ithread] = cy_r38->re;
		_cy_r19[ithread] = cy_r18->im;		_cy_r39[ithread] = cy_r38->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
	#else
		_cy_r00[ithread] = cy_r00;			_cy_r20[ithread] = cy_r20;
		_cy_r01[ithread] = cy_r01;			_cy_r21[ithread] = cy_r21;
		_cy_r02[ithread] = cy_r02;			_cy_r22[ithread] = cy_r22;
		_cy_r03[ithread] = cy_r03;			_cy_r23[ithread] = cy_r23;
		_cy_r04[ithread] = cy_r04;			_cy_r24[ithread] = cy_r24;
		_cy_r05[ithread] = cy_r05;			_cy_r25[ithread] = cy_r25;
		_cy_r06[ithread] = cy_r06;			_cy_r26[ithread] = cy_r26;
		_cy_r07[ithread] = cy_r07;			_cy_r27[ithread] = cy_r27;
		_cy_r08[ithread] = cy_r08;			_cy_r28[ithread] = cy_r28;
		_cy_r09[ithread] = cy_r09;			_cy_r29[ithread] = cy_r29;
		_cy_r10[ithread] = cy_r10;			_cy_r30[ithread] = cy_r30;
		_cy_r11[ithread] = cy_r11;			_cy_r31[ithread] = cy_r31;
		_cy_r12[ithread] = cy_r12;			_cy_r32[ithread] = cy_r32;
		_cy_r13[ithread] = cy_r13;			_cy_r33[ithread] = cy_r33;
		_cy_r14[ithread] = cy_r14;			_cy_r34[ithread] = cy_r34;
		_cy_r15[ithread] = cy_r15;			_cy_r35[ithread] = cy_r35;
		_cy_r16[ithread] = cy_r16;			_cy_r36[ithread] = cy_r36;
		_cy_r17[ithread] = cy_r17;			_cy_r37[ithread] = cy_r37;
		_cy_r18[ithread] = cy_r18;			_cy_r38[ithread] = cy_r38;
		_cy_r19[ithread] = cy_r19;			_cy_r39[ithread] = cy_r39;
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

	}	/******* END OF PARALLEL FOR-LOOP ********/

	if(!full_pass)break;

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-40 forward DIF FFT of the first block of 40 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 40 outputs of (1);
	!   (3) Reweight and perform a radix-40 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 40 elements and repeat (1-4).
	*/
	t00= _cy_r00[CY_THREADS - 1];		t20= _cy_r20[CY_THREADS - 1];
	t01= _cy_r01[CY_THREADS - 1];		t21= _cy_r21[CY_THREADS - 1];
	t02= _cy_r02[CY_THREADS - 1];		t22= _cy_r22[CY_THREADS - 1];
	t03= _cy_r03[CY_THREADS - 1];		t23= _cy_r23[CY_THREADS - 1];
	t04= _cy_r04[CY_THREADS - 1];		t24= _cy_r24[CY_THREADS - 1];
	t05= _cy_r05[CY_THREADS - 1];		t25= _cy_r25[CY_THREADS - 1];
	t06= _cy_r06[CY_THREADS - 1];		t26= _cy_r26[CY_THREADS - 1];
	t07= _cy_r07[CY_THREADS - 1];		t27= _cy_r27[CY_THREADS - 1];
	t08= _cy_r08[CY_THREADS - 1];		t28= _cy_r28[CY_THREADS - 1];
	t09= _cy_r09[CY_THREADS - 1];		t29= _cy_r29[CY_THREADS - 1];
	t10= _cy_r10[CY_THREADS - 1];		t30= _cy_r30[CY_THREADS - 1];
	t11= _cy_r11[CY_THREADS - 1];		t31= _cy_r31[CY_THREADS - 1];
	t12= _cy_r12[CY_THREADS - 1];		t32= _cy_r32[CY_THREADS - 1];
	t13= _cy_r13[CY_THREADS - 1];		t33= _cy_r33[CY_THREADS - 1];
	t14= _cy_r14[CY_THREADS - 1];		t34= _cy_r34[CY_THREADS - 1];
	t15= _cy_r15[CY_THREADS - 1];		t35= _cy_r35[CY_THREADS - 1];
	t16= _cy_r16[CY_THREADS - 1];		t36= _cy_r36[CY_THREADS - 1];
	t17= _cy_r17[CY_THREADS - 1];		t37= _cy_r37[CY_THREADS - 1];
	t18= _cy_r18[CY_THREADS - 1];		t38= _cy_r38[CY_THREADS - 1];
	t19= _cy_r19[CY_THREADS - 1];		t39= _cy_r39[CY_THREADS - 1];

	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		ASSERT(HERE, CY_THREADS > 1,"radix40_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
		_cy_r00[ithread] = _cy_r00[ithread-1];		_cy_r20[ithread] = _cy_r20[ithread-1];
		_cy_r01[ithread] = _cy_r01[ithread-1];		_cy_r21[ithread] = _cy_r21[ithread-1];
		_cy_r02[ithread] = _cy_r02[ithread-1];		_cy_r22[ithread] = _cy_r22[ithread-1];
		_cy_r03[ithread] = _cy_r03[ithread-1];		_cy_r23[ithread] = _cy_r23[ithread-1];
		_cy_r04[ithread] = _cy_r04[ithread-1];		_cy_r24[ithread] = _cy_r24[ithread-1];
		_cy_r05[ithread] = _cy_r05[ithread-1];		_cy_r25[ithread] = _cy_r25[ithread-1];
		_cy_r06[ithread] = _cy_r06[ithread-1];		_cy_r26[ithread] = _cy_r26[ithread-1];
		_cy_r07[ithread] = _cy_r07[ithread-1];		_cy_r27[ithread] = _cy_r27[ithread-1];
		_cy_r08[ithread] = _cy_r08[ithread-1];		_cy_r28[ithread] = _cy_r28[ithread-1];
		_cy_r09[ithread] = _cy_r09[ithread-1];		_cy_r29[ithread] = _cy_r29[ithread-1];
		_cy_r10[ithread] = _cy_r10[ithread-1];		_cy_r30[ithread] = _cy_r30[ithread-1];
		_cy_r11[ithread] = _cy_r11[ithread-1];		_cy_r31[ithread] = _cy_r31[ithread-1];
		_cy_r12[ithread] = _cy_r12[ithread-1];		_cy_r32[ithread] = _cy_r32[ithread-1];
		_cy_r13[ithread] = _cy_r13[ithread-1];		_cy_r33[ithread] = _cy_r33[ithread-1];
		_cy_r14[ithread] = _cy_r14[ithread-1];		_cy_r34[ithread] = _cy_r34[ithread-1];
		_cy_r15[ithread] = _cy_r15[ithread-1];		_cy_r35[ithread] = _cy_r35[ithread-1];
		_cy_r16[ithread] = _cy_r16[ithread-1];		_cy_r36[ithread] = _cy_r36[ithread-1];
		_cy_r17[ithread] = _cy_r17[ithread-1];		_cy_r37[ithread] = _cy_r37[ithread-1];
		_cy_r18[ithread] = _cy_r18[ithread-1];		_cy_r38[ithread] = _cy_r38[ithread-1];
		_cy_r19[ithread] = _cy_r19[ithread-1];		_cy_r39[ithread] = _cy_r39[ithread-1];
	}

	_cy_r00[0] =+t39;	/* ...The wraparound carry is here: */
	_cy_r01[0] = t00;
	_cy_r02[0] = t01;
	_cy_r03[0] = t02;
	_cy_r04[0] = t03;
	_cy_r05[0] = t04;
	_cy_r06[0] = t05;
	_cy_r07[0] = t06;
	_cy_r08[0] = t07;
	_cy_r09[0] = t08;
	_cy_r10[0] = t09;
	_cy_r11[0] = t10;
	_cy_r12[0] = t11;
	_cy_r13[0] = t12;
	_cy_r14[0] = t13;
	_cy_r15[0] = t14;
	_cy_r16[0] = t15;
	_cy_r17[0] = t16;
	_cy_r18[0] = t17;
	_cy_r19[0] = t18;
	_cy_r20[0] = t19;
	_cy_r21[0] = t20;
	_cy_r22[0] = t21;
	_cy_r23[0] = t22;
	_cy_r24[0] = t23;
	_cy_r25[0] = t24;
	_cy_r26[0] = t25;
	_cy_r27[0] = t26;
	_cy_r28[0] = t27;
	_cy_r29[0] = t28;
	_cy_r30[0] = t29;
	_cy_r31[0] = t30;
	_cy_r32[0] = t31;
	_cy_r33[0] = t32;
	_cy_r34[0] = t33;
	_cy_r35[0] = t34;
	_cy_r36[0] = t35;
	_cy_r37[0] = t36;
	_cy_r38[0] = t37;
	_cy_r39[0] = t38;

	full_pass = 0;
	scale = 1;

	jhi = 7;

    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		for(j = ithread*pini; j <= ithread*pini + jhi; j++)
		{
			jt = j;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p08;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p16;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p24;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p32;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
		}
	}
}	/* endfor(outer) */

    t00 = 0;
    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		t00 += fabs(_cy_r00[0])+fabs(_cy_r01[0])+fabs(_cy_r02[0])+fabs(_cy_r03[0])+fabs(_cy_r04[0])+fabs(_cy_r05[0])+fabs(_cy_r06[0])+fabs(_cy_r07[0])+fabs(_cy_r08[0])+fabs(_cy_r09[0])+fabs(_cy_r10[0])+fabs(_cy_r11[0])+fabs(_cy_r12[0])+fabs(_cy_r13[0])+fabs(_cy_r14[0])+fabs(_cy_r15[0])+fabs(_cy_r16[0])+fabs(_cy_r17[0])+fabs(_cy_r18[0])+fabs(_cy_r19[0]);
		t00 += fabs(_cy_r20[0])+fabs(_cy_r21[0])+fabs(_cy_r22[0])+fabs(_cy_r23[0])+fabs(_cy_r24[0])+fabs(_cy_r25[0])+fabs(_cy_r26[0])+fabs(_cy_r27[0])+fabs(_cy_r28[0])+fabs(_cy_r29[0])+fabs(_cy_r30[0])+fabs(_cy_r31[0])+fabs(_cy_r32[0])+fabs(_cy_r33[0])+fabs(_cy_r34[0])+fabs(_cy_r35[0])+fabs(_cy_r36[0])+fabs(_cy_r37[0])+fabs(_cy_r38[0])+fabs(_cy_r39[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
    }

	if(t00 != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix40_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

void radix40_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-40 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in dft_macro.h for details on the radix-5,8 subtransforms.
*/
	int j,j1,j2,jp,jt;
	static int n40,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32, first_entry=TRUE;
	static double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					s2  =  0.95105651629515357211,	/*  sin(u) */
					ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09
	,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
	,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
	,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
	,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
	,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69
	,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
	,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71;

	if(!first_entry && (n/40) != n40)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n40=n/40;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n40;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-40 pass is here.	*/

	for(j=0; j < n40; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
	/*
	Twiddleless version arranges 8 sets of radix-5 DFT inputs as follows: 0 in upper left corner, decrement 8 horizontally and 5 vertically:

		RADIX_05_DFT(00,32,24,16,08)
		RADIX_05_DFT(35,27,19,11,03)
		RADIX_05_DFT(30,22,14,06,38)
		RADIX_05_DFT(25,17,09,01,33)
		RADIX_05_DFT(20,12,04,36,28)
		RADIX_05_DFT(15,07,39,31,23)
		RADIX_05_DFT(10,02,34,26,18)
		RADIX_05_DFT(05,37,29,21,13)

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-8 DFT outputs.
	*/
	/*...gather the needed data (40 64-bit complex) and do 8 radix-5 transforms...*/
            	     /*   sincos     */ /*                                            inputs                                             */ /*          outputs                  */
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[j1    ],a[j2    ],a[j1+p32],a[j2+p32],a[j1+p24],a[j2+p24],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,rt,it);	jt = j1+p03; jp = j2+p03;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,rt,it);	jt = j1+p06; jp = j2+p06;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,rt,it);	jt = j1+p01; jp = j2+p01;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,rt,it);	jt = j1+p07; jp = j2+p07;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,rt,it);	jt = j1+p02; jp = j2+p02;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,rt,it);	jt = j1+p05; jp = j2+p05;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,rt,it)

	/*...and now do 5 radix-8 transforms, swapping the [t*6,t*7] <--> [t*8,t*9] pairs to undo the last-2-outputs-swap in the RADIX_05_DFT macro:	*/
					 /*                            inputs                         */ /*                      intermediates                        */ /*                 outputs                   */
		RADIX_08_DIF(t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],rt,it);	jt = j1+p32; jp = j2+p32;
		RADIX_08_DIF(t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],rt,it);	jt = j1+p24; jp = j2+p24;
		RADIX_08_DIF(t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF(t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIF(t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],rt,it);
	}
}

/***************/

void radix40_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-40 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in dft_macro.h for details on the radix-5,8 subtransforms.
*/
	int j,j1,j2,jp,jt;
	static int n40,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32, first_entry=TRUE;
	static double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					s2  =  0.95105651629515357211,	/*  sin(u) */
					ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09
	,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
	,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
	,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
	,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
	,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69
	,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79
	,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71;

	if(!first_entry && (n/40) != n40)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n40=n/40;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n40;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-40 pass is here.	*/

	for(j=0; j < n40; j += 2)
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

		00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39

		00,32,24,16,08,35,27,19,11,03,30,22,14,06,38,25,17,09,01,33,20,12,04,36,28,15,07,39,31,23,10,02,34,26,18,05,37,29,21,13.	(*)

	Remember, inputs to DIT are bit-reversed, so

		a[00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39] contain [cf. bit-reversal index array output by test_fft_radix()]
		x[00,20,10,30,05,25,15,35,01,21,11,31,06,26,16,36,02,22,12,32,07,27,17,37,03,23,13,33,08,28,18,38,04,24,14,34,09,29,19,39], which get swapped [using the permutation (*)] to
		x[00,20,30,10,35,15,25,05,32,12,22,02,27,07,17,37,24,04,14,34,19,39,09,29,16,36,06,26,11,31,01,21,08,28,38,18,03,23,33,13], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04|19,18,17,16,21,20,22,23|33,32,34,35,38,39,36,37|14,15,12,13,10,11,08,09|28,29,31,30,24,25,27,26]. These are the 5 octets going into the radix-8 DFTs.

		RADIX_08_DFT(00,01,03,02,07,06,05,04)
		RADIX_08_DFT(19,18,17,16,21,20,22,23)
		RADIX_08_DFT(33,32,34,35,38,39,36,37)
		RADIX_08_DFT(14,15,12,13,10,11,08,09)
		RADIX_08_DFT(28,29,31,30,24,25,27,26)

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF] to properly permute the radix-5 DFT outputs.
	*/
	/*...gather the needed data (40 64-bit complex) and do 5 radix-8 transforms,	*/
					 /*                                                                        inputs                                                                             */ /*                 intermediates                             */ /*                   outputs                                 */
		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,rt,it);	jt = j1+p32; jp = j2+p32;
		RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,rt,it);	jt = j1+p24; jp = j2+p24;
		RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,rt,it);

	/*...and now do 8 radix-5 transforms.	*/
					 /*   sincos     */ /*           inputs                  */ /*                                            outputs                                            */
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p32],a[j2+p32],a[j1+p08],a[j2+p08],a[j1+p24],a[j2+p24],rt,it);	jt = j1+p01; jp = j2+p01;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],rt,it);	jt = j1+p02; jp = j2+p02;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it);	jt = j1+p03; jp = j2+p03;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],rt,it);	jt = j1+p05; jp = j2+p05;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],rt,it);	jt = j1+p06; jp = j2+p06;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],rt,it);	jt = j1+p07; jp = j2+p07;
		RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it)
	}
}

