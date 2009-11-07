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

#ifdef USE_SSE2

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

	#ifdef COMPILER_TYPE_MSVC

		#include "sse2_macro.h"

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix24_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix24_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

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
	int n24,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p16;
	static double radix_inv, n2inv;
	static double c = .86602540378443864676, c3m1 = -1.5;
	double scale,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23;
	double maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;

#ifdef USE_SSE2

	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static struct complex *isrt2, *cc0, *cc3, *max_err, *sse2_rnd, *half_arr, *tmp
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r
	,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i,*s1p20i,*s1p21i,*s1p22i,*s1p23i;
	/* Only explicitly reference the odd-indexed carries in SSE2 mode: */
	static struct complex *cy00,*cy02,*cy04,*cy06,*cy08,*cy10,*cy12,*cy14,*cy16,*cy18,*cy20,*cy22;
	static uint64 *sm_arr = 0x0, *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23;

#else

  #if PFETCH
	double *add0, *addr;
  #endif
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double rt,it,temp,frac;
	double wt,wtinv,wtA,wtB,wtC;
	int m,m2;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23;
	double
	 a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r
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

/*...change n24 and n_div_wt to non-static to work around a gcc compiler bug. */
	n24   = n/24;
	n_div_nwt = n24 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n24)
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

		ASSERT(HERE, CY_THREADS >= NTHREADS,"radix24_ditN_cy_dif1.c: CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"radix24_ditN_cy_dif1.c: CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, n24      %CY_THREADS == 0,"radix24_ditN_cy_dif1.c: n24      %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"radix24_ditN_cy_dif1.c: n_div_nwt%CY_THREADS != 0");
		}

	#ifdef MULTITHREAD
		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
	#endif

		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)24));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef USE_SSE2

		sc_arr = ALLOC_COMPLEX(sc_arr, 88);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		/* Size here is [8 + radix/2 + 4] 8-byte elements */
		sm_arr = ALLOC_UINT64(sm_arr, 24);	if(!sm_arr){ sprintf(cbuf, "FATAL: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sm_ptr = ALIGN_UINT64(sm_arr);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 48 16-byte slots of sc_arr for temporaries, next 2 for the doubled cos and c3m1 terms,
	next 12 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
		s1p00r = sc_ptr + 0x00;		isrt2	= sc_ptr + 0x30;
		s1p00i = sc_ptr + 0x01;		cc3		= sc_ptr + 0x31;
		s1p01r = sc_ptr + 0x02;		cc0  	= sc_ptr + 0x32;
		s1p01i = sc_ptr + 0x03;		cy00	= sc_ptr + 0x33;
		s1p02r = sc_ptr + 0x04;		cy02	= sc_ptr + 0x34;
		s1p02i = sc_ptr + 0x05;		cy04	= sc_ptr + 0x35;
		s1p03r = sc_ptr + 0x06;		cy06	= sc_ptr + 0x36;
		s1p03i = sc_ptr + 0x07;		cy08	= sc_ptr + 0x37;
		s1p04r = sc_ptr + 0x08;		cy10	= sc_ptr + 0x38;
		s1p04i = sc_ptr + 0x09;		cy12	= sc_ptr + 0x39;
		s1p05r = sc_ptr + 0x0a;		cy14	= sc_ptr + 0x3a;
		s1p05i = sc_ptr + 0x0b;		cy16	= sc_ptr + 0x3b;
		s1p06r = sc_ptr + 0x0c;		cy18	= sc_ptr + 0x3c;
		s1p06i = sc_ptr + 0x0d;		cy20	= sc_ptr + 0x3d;
		s1p07r = sc_ptr + 0x0e;		cy22	= sc_ptr + 0x3e;
		s1p07i = sc_ptr + 0x0f;		max_err = sc_ptr + 0x3f;
		s1p08r = sc_ptr + 0x10;		sse2_rnd= sc_ptr + 0x40;
		s1p08i = sc_ptr + 0x11;		half_arr= sc_ptr + 0x41;	/* This table needs 20x16 bytes */
		s1p09r = sc_ptr + 0x12;
		s1p09i = sc_ptr + 0x13;
		s1p10r = sc_ptr + 0x14;
		s1p10i = sc_ptr + 0x15;
		s1p11r = sc_ptr + 0x16;
		s1p11i = sc_ptr + 0x17;
		s1p12r = sc_ptr + 0x18;
		s1p12i = sc_ptr + 0x19;
		s1p13r = sc_ptr + 0x1a;
		s1p13i = sc_ptr + 0x1b;
		s1p14r = sc_ptr + 0x1c;
		s1p14i = sc_ptr + 0x1d;
		s1p15r = sc_ptr + 0x1e;
		s1p15i = sc_ptr + 0x1f;
		s1p16r = sc_ptr + 0x20;
		s1p16i = sc_ptr + 0x21;
		s1p17r = sc_ptr + 0x22;
		s1p17i = sc_ptr + 0x23;
		s1p18r = sc_ptr + 0x24;
		s1p18i = sc_ptr + 0x25;
		s1p19r = sc_ptr + 0x26;
		s1p19i = sc_ptr + 0x27;
		s1p20r = sc_ptr + 0x28;
		s1p20i = sc_ptr + 0x29;
		s1p21r = sc_ptr + 0x2a;
		s1p21i = sc_ptr + 0x2b;
		s1p22r = sc_ptr + 0x2c;
		s1p22i = sc_ptr + 0x2d;
		s1p23r = sc_ptr + 0x2e;
		s1p23i = sc_ptr + 0x2f;

		/* These remain fixed: */
		isrt2->re = ISRT2;	isrt2->im = ISRT2;
		cc3  ->re = c3m1 ;	cc3  ->im = c3m1 ;
		cc0  ->re = c    ;	cc0  ->im = c    ;

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

	#endif

	/*   constant index offsets for array load/stores are here.	*/
		p01 = n24;
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
		_jstart  	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co3     == 0x0);

		_cy_00	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_00== 0x0);
		_cy_01	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_01== 0x0);
		_cy_02	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_02== 0x0);
		_cy_03	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_03== 0x0);
		_cy_04	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_04== 0x0);
		_cy_05	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_05== 0x0);
		_cy_06	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_06== 0x0);
		_cy_07	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_07== 0x0);
		_cy_08	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_08== 0x0);
		_cy_09	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_09== 0x0);
		_cy_10	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_10== 0x0);
		_cy_11	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_11== 0x0);
		_cy_12	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_12== 0x0);
		_cy_13	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_13== 0x0);
		_cy_14	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_14== 0x0);
		_cy_15	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_15== 0x0);
		_cy_16	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_16== 0x0);
		_cy_17	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_17== 0x0);
		_cy_18	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_18== 0x0);
		_cy_19	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_19== 0x0);
		_cy_20	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_20== 0x0);
		_cy_21	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_21== 0x0);
		_cy_22	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_22== 0x0);
		_cy_23	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_23== 0x0);

		_maxerr	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix24_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/24-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix24_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		jhi = n24/CY_THREADS;

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
		for(j=0; j < n24; j++)
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

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_jstart[ithread] = ithread*n24/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*24);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+24 - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-24;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
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
	}

    /* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
    if(!full_pass)
    {
		khi = 1;
	#if FFT_DEBUG
		fprintf(dbg_file, "radix24_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifdef MULTITHREAD
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,rt,it,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23) default(shared) schedule(static)
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

	#ifdef USE_SSE2
		*bjmodn00 = _bjmodn00[ithread];				cy00->re = _cy_00[ithread];
		*bjmodn01 = _bjmodn01[ithread];				cy00->im = _cy_01[ithread];
		*bjmodn02 = _bjmodn02[ithread];				cy02->re = _cy_02[ithread];
		*bjmodn03 = _bjmodn03[ithread];				cy02->im = _cy_03[ithread];
		*bjmodn04 = _bjmodn04[ithread];				cy04->re = _cy_04[ithread];
		*bjmodn05 = _bjmodn05[ithread];				cy04->im = _cy_05[ithread];
		*bjmodn06 = _bjmodn06[ithread];				cy06->re = _cy_06[ithread];
		*bjmodn07 = _bjmodn07[ithread];				cy06->im = _cy_07[ithread];
		*bjmodn08 = _bjmodn08[ithread];				cy08->re = _cy_08[ithread];
		*bjmodn09 = _bjmodn09[ithread];				cy08->im = _cy_09[ithread];
		*bjmodn10 = _bjmodn10[ithread];				cy10->re = _cy_10[ithread];
		*bjmodn11 = _bjmodn11[ithread];				cy10->im = _cy_11[ithread];
		*bjmodn12 = _bjmodn12[ithread];				cy12->re = _cy_12[ithread];
		*bjmodn13 = _bjmodn13[ithread];				cy12->im = _cy_13[ithread];
		*bjmodn14 = _bjmodn14[ithread];				cy14->re = _cy_14[ithread];
		*bjmodn15 = _bjmodn15[ithread];				cy14->im = _cy_15[ithread];
		*bjmodn16 = _bjmodn16[ithread];				cy16->re = _cy_16[ithread];
		*bjmodn17 = _bjmodn17[ithread];				cy16->im = _cy_17[ithread];
		*bjmodn18 = _bjmodn18[ithread];				cy18->re = _cy_18[ithread];
		*bjmodn19 = _bjmodn19[ithread];				cy18->im = _cy_19[ithread];
		*bjmodn20 = _bjmodn20[ithread];				cy20->re = _cy_20[ithread];
		*bjmodn21 = _bjmodn21[ithread];				cy20->im = _cy_21[ithread];
		*bjmodn22 = _bjmodn22[ithread];				cy22->re = _cy_22[ithread];
		*bjmodn23 = _bjmodn23[ithread];				cy22->im = _cy_23[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];				cy00 = _cy_00[ithread];
		bjmodn01 = _bjmodn01[ithread];				cy01 = _cy_01[ithread];
		bjmodn02 = _bjmodn02[ithread];				cy02 = _cy_02[ithread];
		bjmodn03 = _bjmodn03[ithread];				cy03 = _cy_03[ithread];
		bjmodn04 = _bjmodn04[ithread];				cy04 = _cy_04[ithread];
		bjmodn05 = _bjmodn05[ithread];				cy05 = _cy_05[ithread];
		bjmodn06 = _bjmodn06[ithread];				cy06 = _cy_06[ithread];
		bjmodn07 = _bjmodn07[ithread];				cy07 = _cy_07[ithread];
		bjmodn08 = _bjmodn08[ithread];				cy08 = _cy_08[ithread];
		bjmodn09 = _bjmodn09[ithread];				cy09 = _cy_09[ithread];
		bjmodn10 = _bjmodn10[ithread];				cy10 = _cy_10[ithread];
		bjmodn11 = _bjmodn11[ithread];				cy11 = _cy_11[ithread];
		bjmodn12 = _bjmodn12[ithread];				cy12 = _cy_12[ithread];
		bjmodn13 = _bjmodn13[ithread];				cy13 = _cy_13[ithread];
		bjmodn14 = _bjmodn14[ithread];				cy14 = _cy_14[ithread];
		bjmodn15 = _bjmodn15[ithread];				cy15 = _cy_15[ithread];
		bjmodn16 = _bjmodn16[ithread];				cy16 = _cy_16[ithread];
		bjmodn17 = _bjmodn17[ithread];				cy17 = _cy_17[ithread];
		bjmodn18 = _bjmodn18[ithread];				cy18 = _cy_18[ithread];
		bjmodn19 = _bjmodn19[ithread];				cy19 = _cy_19[ithread];
		bjmodn20 = _bjmodn20[ithread];				cy20 = _cy_20[ithread];
		bjmodn21 = _bjmodn21[ithread];				cy21 = _cy_21[ithread];
		bjmodn22 = _bjmodn22[ithread];				cy22 = _cy_22[ithread];
		bjmodn23 = _bjmodn23[ithread];				cy23 = _cy_23[ithread];
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
			fprintf(stderr, "radix24_wrapper: A_in[00] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix24_wrapper: A_in[01] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix24_wrapper: A_in[02] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix24_wrapper: A_in[03] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "radix24_wrapper: A_in[04] = %24.5f, %24.5f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "radix24_wrapper: A_in[05] = %24.5f, %24.5f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "radix24_wrapper: A_in[06] = %24.5f, %24.5f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "radix24_wrapper: A_in[07] = %24.5f, %24.5f\n",a[jt+p07],a[jp+p07]);
			jt = j1+p08;	jp = j2+p08;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p04] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p04] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p05] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p05] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p06] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p06] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p07] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p07] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix24_wrapper: A_in[08] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix24_wrapper: A_in[09] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix24_wrapper: A_in[10] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix24_wrapper: A_in[11] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "radix24_wrapper: A_in[12] = %24.5f, %24.5f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "radix24_wrapper: A_in[13] = %24.5f, %24.5f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "radix24_wrapper: A_in[14] = %24.5f, %24.5f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "radix24_wrapper: A_in[15] = %24.5f, %24.5f\n",a[jt+p07],a[jp+p07]);
			jt = j1+p16;	jp = j2+p16;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p04] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p04] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p05] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p05] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p06] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p06] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p07] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p07] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix24_wrapper: A_in[16] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix24_wrapper: A_in[17] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix24_wrapper: A_in[18] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix24_wrapper: A_in[29] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "radix24_wrapper: A_in[20] = %24.5f, %24.5f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "radix24_wrapper: A_in[21] = %24.5f, %24.5f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "radix24_wrapper: A_in[22] = %24.5f, %24.5f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "radix24_wrapper: A_in[23] = %24.5f, %24.5f\n",a[jt+p07],a[jp+p07]);
			fprintf(stderr, "\n");
		#endif

			/*...The radix-24 DIT pass is here:	*/
			/*...gather the needed data (24 64-bit complex, i.e. 48 64-bit reals) and do 3 radix-8 transforms,	*/

	#ifdef USE_SSE2

		#if defined(COMPILER_TYPE_MSVC)

				add0 = &a[j1    ];
				add1 = add0+p01;
				add2 = add0+p03;
				add3 = add0+p02;
				add4 = add0+p07;
				add5 = add0+p06;
				add6 = add0+p05;
				add7 = add0+p04;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p00r, isrt2)

				add5 = &a[j1+p08];
				add0 = add5+p05;
				add1 = add5+p04;
				add2 = add5+p06;
				add3 = add5+p07;
				add4 = add5+p01;
				add6 = add5+p02;
				add7 = add5+p03;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p08r, isrt2)

				add2 = &a[j1+p16];
				add0 = add2+p02;
				add1 = add2+p03;
				add3 = add2+p01;
				add4 = add2+p04;
				add5 = add2+p05;
				add6 = add2+p07;
				add7 = add2+p06;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p16r, isrt2)

				SSE2_RADIX_03_DFT(s1p00r,s1p08r,s1p16r,cc3,s1p00r,s1p16r,s1p08r)
				SSE2_RADIX_03_DFT(s1p01r,s1p09r,s1p17r,cc3,s1p09r,s1p01r,s1p17r)
				SSE2_RADIX_03_DFT(s1p02r,s1p10r,s1p18r,cc3,s1p18r,s1p10r,s1p02r)
				SSE2_RADIX_03_DFT(s1p03r,s1p11r,s1p19r,cc3,s1p03r,s1p19r,s1p11r)
				SSE2_RADIX_03_DFT(s1p04r,s1p12r,s1p20r,cc3,s1p12r,s1p04r,s1p20r)
				SSE2_RADIX_03_DFT(s1p05r,s1p13r,s1p21r,cc3,s1p21r,s1p13r,s1p05r)
				SSE2_RADIX_03_DFT(s1p06r,s1p14r,s1p22r,cc3,s1p06r,s1p22r,s1p14r)
				SSE2_RADIX_03_DFT(s1p07r,s1p15r,s1p23r,cc3,s1p15r,s1p07r,s1p23r)

		#else	/* GCC-style inline ASM: */

			add0 = &a[j1    ];
			SSE2_RADIX24_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p16,s1p00r,isrt2,cc3);

		#endif

	#else

				RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,rt,it);	jt = j1+p08; jp = j2+p08;
				RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,rt,it);	jt = j1+p16; jp = j2+p16;
				RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,rt,it);

			/*...and now do 8 in-place radix-3 transforms.	*/

				RADIX_03_DFT(a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,rt,it);
				RADIX_03_DFT(a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,t01,t02,t03,t04,t05,t06,a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,rt,it);
				RADIX_03_DFT(a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,t01,t02,t03,t04,t05,t06,a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,rt,it);
				RADIX_03_DFT(a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,rt,it);
				RADIX_03_DFT(a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,t01,t02,t03,t04,t05,t06,a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,rt,it);
				RADIX_03_DFT(a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,t01,t02,t03,t04,t05,t06,a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,rt,it);
				RADIX_03_DFT(a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,rt,it);
				RADIX_03_DFT(a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,t01,t02,t03,t04,t05,t06,a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,rt,it);

	#endif	/* USE_SSE2 */

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

				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy20,cy22,bjmodn20);

			#else	/* GCC-style inline ASM: */

				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

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

				SSE2_cmplx_carry_norm_nocheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p20r,add1,add2,     cy20,cy22,bjmodn20);

			#else	/* GCC-style inline ASM: */

				SSE2_cmplx_carry_norm_nocheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

			#endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

		#else	/* USE_SSE2 */

		/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy01,bjmodn01,1 );
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy02,bjmodn02,2 );
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy03,bjmodn03,3 );
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy04,bjmodn04,4 );
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy05,bjmodn05,5 );
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy06,bjmodn06,6 );
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy07,bjmodn07,7 );
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy08,bjmodn08,8 );
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy09,bjmodn09,9 );
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy23,bjmodn23,23);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					   and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	/* USE_SSE2 */

		#ifdef USE_SSE2

			#if defined(COMPILER_TYPE_MSVC)

				SSE2_RADIX_03_DFT(s1p00r,s1p16r,s1p08r,cc3,s1p00r,s1p08r,s1p16r)
				SSE2_RADIX_03_DFT(s1p21r,s1p13r,s1p05r,cc3,s1p05r,s1p13r,s1p21r)
				SSE2_RADIX_03_DFT(s1p18r,s1p10r,s1p02r,cc3,s1p02r,s1p10r,s1p18r)
				SSE2_RADIX_03_DFT(s1p15r,s1p07r,s1p23r,cc3,s1p07r,s1p15r,s1p23r)
				SSE2_RADIX_03_DFT(s1p12r,s1p04r,s1p20r,cc3,s1p04r,s1p12r,s1p20r)
				SSE2_RADIX_03_DFT(s1p09r,s1p01r,s1p17r,cc3,s1p01r,s1p09r,s1p17r)
				SSE2_RADIX_03_DFT(s1p06r,s1p22r,s1p14r,cc3,s1p06r,s1p14r,s1p22r)
				SSE2_RADIX_03_DFT(s1p03r,s1p19r,s1p11r,cc3,s1p03r,s1p11r,s1p19r)

				add0 = &a[j1    ];
				add1 = add0+p01;
				add2 = add0+p02;
				add3 = add0+p03;
				add4 = add0+p05;
				add5 = add0+p04;
				add6 = add0+p07;
				add7 = add0+p06;
				SSE2_RADIX8_DIF_0TWIDDLE(s1p00r,0xa0,0x40,0xe0,0x80,0x20,0xc0,0x60, add0,add1,add2,add3,add4,add5,add6,add7, isrt2)

				add3 = &a[j1+p16];
				add0 = add3+p02;
				add1 = add3+p03;
				add2 = add3+p01;
				add4 = add3+p07;
				add5 = add3+p06;
				add6 = add3+p04;
				add7 = add3+p05;
				SSE2_RADIX8_DIF_0TWIDDLE(s1p08r,0xa0,0x40,0xe0,0x80,0x20,0xc0,0x60, add0,add1,add2,add3,add4,add5,add6,add7, isrt2)

				add7 = &a[j1+p08];
				add0 = add7+p05;
				add1 = add7+p04;
				add2 = add7+p07;
				add3 = add7+p06;
				add4 = add7+p02;
				add5 = add7+p03;
				add6 = add7+p01;
				SSE2_RADIX8_DIF_0TWIDDLE(s1p16r,0xa0,0x40,0xe0,0x80,0x20,0xc0,0x60, add0,add1,add2,add3,add4,add5,add6,add7, isrt2)

			#else	/* GCC-style inline ASM: */

				add0 = &a[j1    ];
				SSE2_RADIX24_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p16,s1p00r,isrt2,cc3);

			#endif

	#ifdef DEBUG_SSE2
			jt = j1;		jp = j2;
			fprintf(stderr, "radix24_wrapper: A_out[00] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix24_wrapper: A_out[01] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix24_wrapper: A_out[02] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix24_wrapper: A_out[03] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "radix24_wrapper: A_out[04] = %24.5f, %24.5f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "radix24_wrapper: A_out[05] = %24.5f, %24.5f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "radix24_wrapper: A_out[06] = %24.5f, %24.5f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "radix24_wrapper: A_out[07] = %24.5f, %24.5f\n",a[jt+p07],a[jp+p07]);
			jt = j1+p08;	jp = j2+p08;
			fprintf(stderr, "radix24_wrapper: A_out[08] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix24_wrapper: A_out[09] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix24_wrapper: A_out[10] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix24_wrapper: A_out[11] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "radix24_wrapper: A_out[12] = %24.5f, %24.5f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "radix24_wrapper: A_out[13] = %24.5f, %24.5f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "radix24_wrapper: A_out[14] = %24.5f, %24.5f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "radix24_wrapper: A_out[15] = %24.5f, %24.5f\n",a[jt+p07],a[jp+p07]);
			jt = j1+p16;	jp = j2+p16;
			fprintf(stderr, "radix24_wrapper: A_out[16] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix24_wrapper: A_out[17] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix24_wrapper: A_out[18] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix24_wrapper: A_out[29] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "radix24_wrapper: A_out[20] = %24.5f, %24.5f\n",a[jt+p04],a[jp+p04]);
			fprintf(stderr, "radix24_wrapper: A_out[21] = %24.5f, %24.5f\n",a[jt+p05],a[jp+p05]);
			fprintf(stderr, "radix24_wrapper: A_out[22] = %24.5f, %24.5f\n",a[jt+p06],a[jp+p06]);
			fprintf(stderr, "radix24_wrapper: A_out[23] = %24.5f, %24.5f\n",a[jt+p07],a[jp+p07]);
			fprintf(stderr, "\n");
		exit(0);
	#endif

		#else

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
				RADIX_03_DFT_PFETCH(a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,rt,it,p01);
				RADIX_03_DFT_PFETCH(a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,t01,t02,t03,t04,t05,t06,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,rt,it,p02);
				RADIX_03_DFT_PFETCH(a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,t01,t02,t03,t04,t05,t06,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,rt,it,p03);
				RADIX_03_DFT_PFETCH(a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,t01,t02,t03,t04,t05,t06,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,rt,it,p04);
				RADIX_03_DFT_PFETCH(a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,t01,t02,t03,t04,t05,t06,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,rt,it,p05);
				RADIX_03_DFT_PFETCH(a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,t01,t02,t03,t04,t05,t06,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,rt,it,p06);
				RADIX_03_DFT_PFETCH(a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,rt,it,p07);
				RADIX_03_DFT_PFETCH(a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,rt,it,p08);
			#else
				RADIX_03_DFT       (a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,rt,it);
				RADIX_03_DFT       (a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,t01,t02,t03,t04,t05,t06,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,rt,it);
				RADIX_03_DFT       (a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,t01,t02,t03,t04,t05,t06,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,rt,it);
				RADIX_03_DFT       (a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,t01,t02,t03,t04,t05,t06,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,rt,it);
				RADIX_03_DFT       (a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,t01,t02,t03,t04,t05,t06,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,rt,it);
				RADIX_03_DFT       (a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,t01,t02,t03,t04,t05,t06,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,rt,it);
				RADIX_03_DFT       (a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,rt,it);
				RADIX_03_DFT       (a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,rt,it);
			#endif

			/*...and now do 3 radix-8 transforms:	*/
								 /*                                                          inputs                                                           */ /*                 intermediates                    */ /*                 outputs                   */
			#if PFETCH
				RADIX_08_DIF_PFETCH(a1p00r,a1p00i,a1p05r,a1p05i,a1p02r,a1p02i,a1p07r,a1p07i,a1p04r,a1p04i,a1p01r,a1p01i,a1p06r,a1p06i,a1p03r,a1p03i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],rt,it,p08+p01,p08+p02,p08+p03,p08+p04,p08+p05);	jt = j1+p16; jp = j2+p16;
				RADIX_08_DIF_PFETCH(a1p08r,a1p08i,a1p13r,a1p13i,a1p10r,a1p10i,a1p15r,a1p15i,a1p12r,a1p12i,a1p09r,a1p09i,a1p14r,a1p14i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it,p08+p06,p08+p07,p16    ,p16+p01,p16+p02);	jt = j1+p08; jp = j2+p08;
				RADIX_08_DIF_PFETCH(a1p16r,a1p16i,a1p21r,a1p21i,a1p18r,a1p18i,a1p23r,a1p23i,a1p20r,a1p20i,a1p17r,a1p17i,a1p22r,a1p22i,a1p19r,a1p19i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,p16+p03,p16+p04,p16+p05,p16+p06,p16+p07);
			#else
				RADIX_08_DIF       (a1p00r,a1p00i,a1p05r,a1p05i,a1p02r,a1p02i,a1p07r,a1p07i,a1p04r,a1p04i,a1p01r,a1p01i,a1p06r,a1p06i,a1p03r,a1p03i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],rt,it);	jt = j1+p16; jp = j2+p16;
				RADIX_08_DIF       (a1p08r,a1p08i,a1p13r,a1p13i,a1p10r,a1p10i,a1p15r,a1p15i,a1p12r,a1p12i,a1p09r,a1p09i,a1p14r,a1p14i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it);	jt = j1+p08; jp = j2+p08;
				RADIX_08_DIF       (a1p16r,a1p16i,a1p21r,a1p21i,a1p18r,a1p18i,a1p23r,a1p23i,a1p20r,a1p20i,a1p17r,a1p17i,a1p22r,a1p22i,a1p19r,a1p19i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
			#endif

		#endif	/* USE_SSE2 */
			}

			jstart += nwt;
			jhi    += nwt;
			col += 24;
			co3 -= 24;

		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_SSE2
		_cy_00[ithread] = cy00->re;
		_cy_01[ithread] = cy00->im;
		_cy_02[ithread] = cy02->re;
		_cy_03[ithread] = cy02->im;
		_cy_04[ithread] = cy04->re;
		_cy_05[ithread] = cy04->im;
		_cy_06[ithread] = cy06->re;
		_cy_07[ithread] = cy06->im;
		_cy_08[ithread] = cy08->re;
		_cy_09[ithread] = cy08->im;
		_cy_10[ithread] = cy10->re;
		_cy_11[ithread] = cy10->im;
		_cy_12[ithread] = cy12->re;
		_cy_13[ithread] = cy12->im;
		_cy_14[ithread] = cy14->re;
		_cy_15[ithread] = cy14->im;
		_cy_16[ithread] = cy16->re;
		_cy_17[ithread] = cy16->im;
		_cy_18[ithread] = cy18->re;
		_cy_19[ithread] = cy18->im;
		_cy_20[ithread] = cy20->re;
		_cy_21[ithread] = cy20->im;
		_cy_22[ithread] = cy22->re;
		_cy_23[ithread] = cy22->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
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

	}	/******* END OF PARALLEL FOR-LOOP ********/

	if(!full_pass)break;

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
		{								jt=j+p08;					jp=j+p16;
			a[j     ] *= radix_inv;		a[jt    ] *= radix_inv;		a[jp    ] *= radix_inv;
			a[j +p01] *= radix_inv;		a[jt+p01] *= radix_inv;		a[jp+p01] *= radix_inv;
			a[j +p02] *= radix_inv;		a[jt+p02] *= radix_inv;		a[jp+p02] *= radix_inv;
			a[j +p03] *= radix_inv;		a[jt+p03] *= radix_inv;		a[jp+p03] *= radix_inv;
			a[j +p04] *= radix_inv;		a[jt+p04] *= radix_inv;		a[jp+p04] *= radix_inv;
			a[j +p05] *= radix_inv;		a[jt+p05] *= radix_inv;		a[jp+p05] *= radix_inv;
			a[j +p06] *= radix_inv;		a[jt+p06] *= radix_inv;		a[jp+p06] *= radix_inv;
			a[j +p07] *= radix_inv;		a[jt+p07] *= radix_inv;		a[jp+p07] *= radix_inv;
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
	int n24, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double radix_inv, n2inv;
	static double c = .86602540378443864676, c3m1 = -1.5;
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

/*...change n24 and n_div_wt to non-static to work around a gcc compiler bug. */
	n24   = n/24;
	n_div_nwt = n24 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n24)
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

	  p01 = n24;
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
	  for(j=0; j < n24; j++)
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
	RADIX_03_DFT(a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,t1,t2,t3,t4,t5,t6,a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,rt,it);
	RADIX_03_DFT(a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,t1,t2,t3,t4,t5,t6,a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,rt,it);
	RADIX_03_DFT(a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,t1,t2,t3,t4,t5,t6,a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,rt,it);
	RADIX_03_DFT(a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,t1,t2,t3,t4,t5,t6,a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,rt,it);
	RADIX_03_DFT(a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,t1,t2,t3,t4,t5,t6,a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,rt,it);
	RADIX_03_DFT(a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,t1,t2,t3,t4,t5,t6,a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,rt,it);
	RADIX_03_DFT(a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,t1,t2,t3,t4,t5,t6,a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,rt,it);
	RADIX_03_DFT(a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,t1,t2,t3,t4,t5,t6,a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,rt,it);

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
	RADIX_03_DFT_PFETCH(a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,t1,t2,t3,t4,t5,t6,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,rt,it,p01);
	RADIX_03_DFT_PFETCH(a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,t1,t2,t3,t4,t5,t6,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,rt,it,p02);
	RADIX_03_DFT_PFETCH(a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,t1,t2,t3,t4,t5,t6,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,rt,it,p03);
	RADIX_03_DFT_PFETCH(a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,t1,t2,t3,t4,t5,t6,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,rt,it,p04);
	RADIX_03_DFT_PFETCH(a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,t1,t2,t3,t4,t5,t6,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,rt,it,p05);
	RADIX_03_DFT_PFETCH(a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,t1,t2,t3,t4,t5,t6,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,rt,it,p06);
	RADIX_03_DFT_PFETCH(a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,t1,t2,t3,t4,t5,t6,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,rt,it,p07);
	RADIX_03_DFT_PFETCH(a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,rt,it,p08);
#else
	RADIX_03_DFT       (a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,t1,t2,t3,t4,t5,t6,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,rt,it);
	RADIX_03_DFT       (a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,t1,t2,t3,t4,t5,t6,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,rt,it);
	RADIX_03_DFT       (a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,t1,t2,t3,t4,t5,t6,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,rt,it);
	RADIX_03_DFT       (a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,t1,t2,t3,t4,t5,t6,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,rt,it);
	RADIX_03_DFT       (a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,t1,t2,t3,t4,t5,t6,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,rt,it);
	RADIX_03_DFT       (a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,t1,t2,t3,t4,t5,t6,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,rt,it);
	RADIX_03_DFT       (a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,t1,t2,t3,t4,t5,t6,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,rt,it);
	RADIX_03_DFT       (a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,rt,it);
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
	static int n24,p01,p02,p03,p04,p05,p06,p07,p08,p16, first_entry=TRUE;
	static double c = .86602540378443864676, c3m1 = -1.5;
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i;

	if(!first_entry && (n/24) != n24)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n24=n/24;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = n24;
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

	for(j=0; j < n24; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
	/*
	Twiddleless version arranges 8 sets of radix-5 DFT inputs as follows: 0 in upper left corner, decrement 8 horizontally and 3 vertically:

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
		RADIX_03_DFT(a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08],a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,rt,it);	jt = j1+p05; jp = j2+p05;
		RADIX_03_DFT(a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,rt,it);	jt = j1+p02; jp = j2+p02;
		RADIX_03_DFT(a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,rt,it);	jt = j1+p07; jp = j2+p07;
		RADIX_03_DFT(a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_03_DFT(a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,rt,it);	jt = j1+p01; jp = j2+p01;
		RADIX_03_DFT(a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,rt,it);	jt = j1+p06; jp = j2+p06;
		RADIX_03_DFT(a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,rt,it);	jt = j1+p03; jp = j2+p03;
		RADIX_03_DFT(a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,rt,it);

	/*...and now do 3 radix-8 transforms:	*/
					 /*                                                          inputs                                                           */ /*                 intermediates                    */ /*                 outputs                   */
		RADIX_08_DIF(a1p00r,a1p00i,a1p05r,a1p05i,a1p02r,a1p02i,a1p07r,a1p07i,a1p04r,a1p04i,a1p01r,a1p01i,a1p06r,a1p06i,a1p03r,a1p03i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIF(a1p08r,a1p08i,a1p13r,a1p13i,a1p10r,a1p10i,a1p15r,a1p15i,a1p12r,a1p12i,a1p09r,a1p09i,a1p14r,a1p14i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF(a1p16r,a1p16i,a1p21r,a1p21i,a1p18r,a1p18i,a1p23r,a1p23i,a1p20r,a1p20i,a1p17r,a1p17i,a1p22r,a1p22i,a1p19r,a1p19i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
	}
}

/***************/

void radix24_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-24 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2,jt,jp;
	static int n24,p01,p02,p03,p04,p05,p06,p07,p08,p16, first_entry=TRUE;
	static double c = .86602540378443864676, c3m1 = -1.5;
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i;

	if(!first_entry && (n/24) != n24)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n24=n/24;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = n24;
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

	for(j=0; j < n24; j += 2)
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

		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,rt,it);

	/*...and now do 8 radix-3 transforms.	*/

		RADIX_03_DFT(a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08],rt,it);	jt = j1+p01; jp = j2+p01;
		RADIX_03_DFT(a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);	jt = j1+p02; jp = j2+p02;
		RADIX_03_DFT(a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],rt,it);	jt = j1+p03; jp = j2+p03;
		RADIX_03_DFT(a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_03_DFT(a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);	jt = j1+p05; jp = j2+p05;
		RADIX_03_DFT(a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],rt,it);	jt = j1+p06; jp = j2+p06;
		RADIX_03_DFT(a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],rt,it);	jt = j1+p07; jp = j2+p07;
		RADIX_03_DFT(a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);
	}
}

