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

#undef FFT_DEBUG
#define FFT_DEBUG	0

#if FFT_DEBUG
	#include "carry_dbg.h"
#endif

#ifdef USE_SSE2

	#include "sse2_macro.h"

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

//	#define	USE_SCALAR_CARRY	// Uncomment if want to use original non-SSE carry macros

	#if 0//COMPILER_TYPE_MSVC

		/*
		In the Fermat-mod negacyclic-DWT carry scheme, real & imaginary parts
		are carried separately due to the right-angle transform:
		*/
		#define SSE2_fermat_carry_norm_pow2_errcheck(__data,__cy,__iscale,__half_arr,........idx_offset)\
		{\
				__asm	mov	eax, __data\
				__asm	mov	ebx, __half_arr\
				/* Data enter in R0 = [a0.re,b0.re], I0 = [a0.im,b0.im]: */\
				__asm	movaps	xmm0,[eax     ]	/* R0 */\
				__asm	movaps	xmm1,[eax+0x10]	/* I0 */\
				__asm	mulpd	xmm0,[ebx+__iscale]	/* R0 *= scale */\
				__asm	mulpd	xmm1,[ebx+__iscale]	/* I0 *= scale */\
			/* Get the needed Nth root of -1: */\
				l = ((j + idx_offset) >> 1);\
				k1=(l & NRTM1);\
				k2=(l >> NRT_BITS);\
				temp=rn0[k1].re;		wt_im=rn0[k1].im;\
				rt  =rn1[k2].re;		it   =rn1[k2].im;\
				wt_re =temp*rt-wt_im*it;wt_im =temp*it+wt_im*rt;\
				__asm	mov	edi, __wt_re_ptr\
				__asm	mov	esi, __wt_im_ptr\
				__asm	movlpd	xmm2,[edi]\
				__asm	movlpd	xmm3,[esi]\
				l += 2;\
				k1=(l & NRTM1);\
				k2=(l >> NRT_BITS);\
				temp=rn0[k1].re;		wt_im=rn0[k1].im;\
				rt  =rn1[k2].re;		it   =rn1[k2].im;\
				wt_re =temp*rt-wt_im*it;wt_im =temp*it+wt_im*rt;\
				__asm	mov	edi, __wt_re_ptr\
				__asm	mov	esi, __wt_im_ptr\
				__asm	movhpd	xmm4,[edi]\
				__asm	movhpd	xmm5,[esi]\
				/* Init carries: (scalar) carries cx,cy will get added into low halves, high half just gets unweighted: */\
				__asm	mov	edi, &__cx\
				__asm	mov	esi, &__cy\
				__asm	xorpd	xmm6,xmm6	/* zero */\
				__asm	xorpd	xmm7,xmm7	/* zero */\
				__asm	movlpd	xmm6,[edi]	/* cx in low half */\
				__asm	movlpd	xmm7,[esi]	/* cy in low half */\
				\
			/* Inverse weight is (wt_re, -wt_im): */\
				__asm	movaps	xmm2,xmm0	/* cpy R0 */\
				__asm	movaps	xmm3,xmm1	/* cpy I0 */\
				__asm	mulpd	xmm0,xmm4	/* rt*wt_re */\
				__asm	mulpd	xmm1,xmm4	/* it*wt_re */\
				__asm	mulpd	xmm3,xmm5	/* it*wt_im */\
				__asm	mulpd	xmm2,xmm5	/* rt*wt_im */\
				__asm	addpd	xmm0,xmm3	/* x = rt*wt_re - it*wt_im */\
				__asm	subpd	xmm1,xmm2	/* y = rt*wt_im + it*wt_re */\
				__asm	movaps	xmm2,[ebx-0x10]	/* rnd_const */\
				__asm	addpd	xmm0,xmm6	/* rt = x*wt_re + y*wt_im + cx */\
				__asm	addpd	xmm1,xmm7	/* it = y*wt_re - x*wt_im + cy */\
				__asm	mulpd	xmm0,[ebx+0x40]	/* R0 *= [baseinv,1] */\
				__asm	mulpd	xmm1,[ebx+0x40]	/* I0 *= [baseinv,1]; don't need to round to compute carry into high halves */\
				__asm	haddpd	xmm0,xmm1	/* horizontal-add of carry from low into high halves */\
				__asm	mulpd	xmm0,[ebx+0x50]	/* R0 *= [1,baseinv] */\
				__asm	mulpd	xmm1,[ebx+0x50]	/* I0 *= [1,baseinv]; don't need to round to compute carry into high halves */\
				__asm	movaps	xmm0,xmm6	/* cpy R0 */\
				__asm	movaps	xmm1,xmm7	/* cpy I0 */\
				__asm	addpd	xmm6,xmm2\
				__asm	addpd	xmm7,xmm2\
				__asm	subpd	xmm6,xmm2	/* cx */\
				__asm	subpd	xmm7,xmm2	/* cy */\
				__asm	movhpd	[edi],xmm6	/* cx in high half */\
				__asm	movhpd	[esi],xmm7	/* cy in high half */\
				__asm	mulpd	xmm6,[ebx+0x30]	/* cx *= [base,base] */\
				__asm	mulpd	xmm7,[ebx+0x30]	/* cy *= [base,base] */\
				__asm	subpd	xmm0,xmm6	/* rt = temp-cx*base[0] */\
				__asm	subpd	xmm1,xmm7	/* it = temp-cy*base[0] */\
			/* Forward weight is (wt_re, +wt_im): */\
				__asm	movaps	xmm2,xmm0	/* cpy R0 */\
				__asm	movaps	xmm3,xmm1	/* cpy I0 */\
				__asm	mulpd	xmm0,xmm4	/* rt*wt_re */\
				__asm	mulpd	xmm1,xmm4	/* it*wt_re */\
				__asm	mulpd	xmm3,xmm5	/* it*wt_im */\
				__asm	mulpd	xmm2,xmm5	/* rt*wt_im */\
				__asm	subpd	xmm0,xmm3	/* x = rt*wt_re - it*wt_im */\
				__asm	addpd	xmm1,xmm2	/* y = rt*wt_im + it*wt_re */\
		}

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix16_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix16_ditN_cy_dif1_gcc64.h"

		#endif

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
	static int n16;
	int i,j,j1,j2,jstart,jhi,full_pass,k1,k2,k,khi,l,ntmp,outer;
	static uint32 bjmodnini;
	static uint64 psave=0;
	static uint32 bw,sw,nm1,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599, radix_inv,n2inv;
	double rt,it;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	double maxerr = 0.0;
#if PFETCH
	double *addr, *addp;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	static double wt_re,wt_im;									/* Fermat-mod weights stuff */
	static double *wt_re_ptr = &wt_re, *wt_im_ptr = &wt_im;

#ifdef USE_SSE2

	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	double scale, *add0, *add1, *add2, *add3;	/* Addresses into array sections */
#if 1
	int *si_ptr = &si[0];
	double *wt0_ptr = &wt0[0], *wt1_ptr = &wt1[0], *scale_ptr = &scale;
#endif

	static struct complex *sc_arr = 0x0,*sc_ptr;
	static struct complex *cc0, *ss0, *isrt2, *max_err, *sse2_rnd, *half_arr, *tmp;
	static struct complex *r1,*r2,*r3,*r4,*r5,*r6,*r7,*r8,*r9,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r30,*r31,*r32;
	static struct complex *cy_r01,*cy_r23,*cy_r45,*cy_r67,*cy_r89,*cy_rAB,*cy_rCD,*cy_rEF,*cy_i01,*cy_i23,*cy_i45,*cy_i67,*cy_i89,*cy_iAB,*cy_iCD,*cy_iEF;
	static uint64 *sm_arr = 0x0, *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
	uint64 tmp64;
	static int *bjmodn0,*bjmodn1,*bjmodn2,*bjmodn3,*bjmodn4,*bjmodn5,*bjmodn6,*bjmodn7,*bjmodn8,*bjmodn9,*bjmodnA,*bjmodnB,*bjmodnC,*bjmodnD,*bjmodnE,*bjmodnF;

  #ifdef DEBUG_SSE2
	int jt,jp;
  #endif

  #ifdef USE_SCALAR_CARRY
	double temp,frac
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi
		,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF
		,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF;
  #endif

#else

	int bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF;
	double temp,frac,scale
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi
		,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF
		,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS = 0,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn0 = 0x0,*_bjmodn1 = 0x0,*_bjmodn2 = 0x0,*_bjmodn3 = 0x0,*_bjmodn4 = 0x0,*_bjmodn5 = 0x0,*_bjmodn6 = 0x0,*_bjmodn7 = 0x0,*_bjmodn8 = 0x0,*_bjmodn9 = 0x0,*_bjmodnA = 0x0,*_bjmodnB = 0x0,*_bjmodnC = 0x0,*_bjmodnD = 0x0,*_bjmodnE = 0x0,*_bjmodnF = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy_r0 = 0x0,*_cy_r1 = 0x0,*_cy_r2 = 0x0,*_cy_r3 = 0x0,*_cy_r4 = 0x0,*_cy_r5 = 0x0,*_cy_r6 = 0x0,*_cy_r7 = 0x0,*_cy_r8 = 0x0,*_cy_r9 = 0x0,*_cy_rA = 0x0,*_cy_rB = 0x0,*_cy_rC = 0x0,*_cy_rD = 0x0,*_cy_rE = 0x0,*_cy_rF = 0x0,
	*_cy_i0 = 0x0,*_cy_i1 = 0x0,*_cy_i2 = 0x0,*_cy_i3 = 0x0,*_cy_i4 = 0x0,*_cy_i5 = 0x0,*_cy_i6 = 0x0,*_cy_i7 = 0x0,*_cy_i8 = 0x0,*_cy_i9 = 0x0,*_cy_iA = 0x0,*_cy_iB = 0x0,*_cy_iC = 0x0,*_cy_iD = 0x0,*_cy_iE = 0x0,*_cy_iF = 0x0;

#if FFT_DEBUG
	int len_a;
	FILE *dbg_file;
	const char dbg_fname[] = "CY16_MT_DEBUG.txt";
	ASSERT(HERE, n < 10000, "n too large!");
	dbg_file = fopen(dbg_fname, "a");
	fprintf(dbg_file, "radix16_ditN_cy_dif1: Multithreaded, n = %d\n", n);
#endif

/*...change n16 and n_div_wt to non-static to work around a gcc compiler bug. */
	n16   = n/16;
	n_div_nwt = n16 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n16)
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

		ASSERT(HERE, CY_THREADS >= NTHREADS,"radix16_ditN_cy_dif1.c: CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"radix16_ditN_cy_dif1.c: CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, n16      %CY_THREADS == 0,"radix16_ditN_cy_dif1.c: n32      %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"radix16_ditN_cy_dif1.c: n_div_nwt%CY_THREADS != 0");
		}

	#ifdef MULTITHREAD
		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
	#endif

		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)16));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		nm1   = n-1;

	#ifdef USE_SSE2
	  #ifndef USE_SCALAR_CARRY
		ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE), "SSE2 currently only supports Mersenne-mod!");
	  #endif
		sc_arr = ALLOC_COMPLEX(sc_arr, 80);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		sm_arr = ALLOC_UINT64(sm_arr, 12);	if(!sm_arr){ sprintf(cbuf, "FATAL: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sm_ptr = ALIGN_UINT64(sm_arr);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	next 16 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
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
		isrt2->re = ISRT2;	isrt2->im = ISRT2;
		cc0  ->re = c	;	cc0  ->im = c	;
		ss0  ->re = s	;	ss0  ->im = s	;
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
	if(TRANSFORM_TYPE == RIGHT_ANGLE)	/* In Fermat-mod mode, use first 2 128-bit slots for the needed combos of the 2 scaling factors 1/n2 and 1: */
	{
		tmp->re = n2inv;	tmp->im = n2inv;	++tmp;
		/* In wraparound-carry step, scale = 1: */
		tmp->re =   1.0;	tmp->im =   1.0;	++tmp;
		/* Forward-base[] multipliers: */
/*		tmp->re = base   [0];	tmp->im =        1.0;	++tmp;*/
		tmp->re =        1.0;	tmp->im = base   [0];	++tmp;
		tmp->re = base   [0];	tmp->im = base   [0];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->re = baseinv[0];	tmp->im =        1.0;	++tmp;
		tmp->re =        1.0;	tmp->im = baseinv[0];	++tmp;
	}
	else
	{
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

	#endif

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
		_i       	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_i== 0x0);
		_bjmodn0	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn0== 0x0);
		_bjmodn1	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn1== 0x0);
		_bjmodn2	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn2== 0x0);
		_bjmodn3	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn3== 0x0);
		_bjmodn4	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn4== 0x0);
		_bjmodn5	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn5== 0x0);
		_bjmodn6	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn6== 0x0);
		_bjmodn7	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn7== 0x0);
		_bjmodn8	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn8== 0x0);
		_bjmodn9	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn9== 0x0);
		_bjmodnA	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodnA== 0x0);
		_bjmodnB	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodnB== 0x0);
		_bjmodnC	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodnC== 0x0);
		_bjmodnD	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodnD== 0x0);
		_bjmodnE	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodnE== 0x0);
		_bjmodnF	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodnF== 0x0);
		_jstart  	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co3     == 0x0);

		_cy_r0	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r0== 0x0);
		_cy_r1	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r1== 0x0);
		_cy_r2	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r2== 0x0);
		_cy_r3	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r3== 0x0);
		_cy_r4	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r4== 0x0);
		_cy_r5	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r5== 0x0);
		_cy_r6	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r6== 0x0);
		_cy_r7	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r7== 0x0);
		_cy_r8	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r8== 0x0);
		_cy_r9	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r9== 0x0);
		_cy_rA	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_rA== 0x0);
		_cy_rB	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_rB== 0x0);
		_cy_rC	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_rC== 0x0);
		_cy_rD	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_rD== 0x0);
		_cy_rE	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_rE== 0x0);
		_cy_rF	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_rF== 0x0);

		_cy_i0	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i0== 0x0);
		_cy_i1	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i1== 0x0);
		_cy_i2	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i2== 0x0);
		_cy_i3	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i3== 0x0);
		_cy_i4	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i4== 0x0);
		_cy_i5	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i5== 0x0);
		_cy_i6	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i6== 0x0);
		_cy_i7	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i7== 0x0);
		_cy_i8	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i8== 0x0);
		_cy_i9	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i9== 0x0);
		_cy_iA	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_iA== 0x0);
		_cy_iB	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_iB== 0x0);
		_cy_iC	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_iC== 0x0);
		_cy_iD	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_iD== 0x0);
		_cy_iE	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_iE== 0x0);
		_cy_iF	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_iF== 0x0);

		_maxerr	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix16_ditN_cy_dif1.");

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
			i.e. the one that n2/16-separated FFT outputs need:
			*/
			_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix16_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r0[      0] = -2;
	}

	*fracmax = 0;	/* init max. fractional error	*/
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
			_bjmodn0[ithread] = _bjmodnini[ithread];
			_bjmodn1[ithread] = _bjmodn0[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn1[ithread] = _bjmodn1[ithread] + ( (-(int)((uint32)_bjmodn1[ithread] >> 31)) & n);
			_bjmodn2[ithread] = _bjmodn1[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn2[ithread] = _bjmodn2[ithread] + ( (-(int)((uint32)_bjmodn2[ithread] >> 31)) & n);
			_bjmodn3[ithread] = _bjmodn2[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn3[ithread] = _bjmodn3[ithread] + ( (-(int)((uint32)_bjmodn3[ithread] >> 31)) & n);
			_bjmodn4[ithread] = _bjmodn3[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn4[ithread] = _bjmodn4[ithread] + ( (-(int)((uint32)_bjmodn4[ithread] >> 31)) & n);
			_bjmodn5[ithread] = _bjmodn4[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn5[ithread] = _bjmodn5[ithread] + ( (-(int)((uint32)_bjmodn5[ithread] >> 31)) & n);
			_bjmodn6[ithread] = _bjmodn5[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn6[ithread] = _bjmodn6[ithread] + ( (-(int)((uint32)_bjmodn6[ithread] >> 31)) & n);
			_bjmodn7[ithread] = _bjmodn6[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn7[ithread] = _bjmodn7[ithread] + ( (-(int)((uint32)_bjmodn7[ithread] >> 31)) & n);
			_bjmodn8[ithread] = _bjmodn7[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn8[ithread] = _bjmodn8[ithread] + ( (-(int)((uint32)_bjmodn8[ithread] >> 31)) & n);
			_bjmodn9[ithread] = _bjmodn8[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn9[ithread] = _bjmodn9[ithread] + ( (-(int)((uint32)_bjmodn9[ithread] >> 31)) & n);
			_bjmodnA[ithread] = _bjmodn9[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodnA[ithread] = _bjmodnA[ithread] + ( (-(int)((uint32)_bjmodnA[ithread] >> 31)) & n);
			_bjmodnB[ithread] = _bjmodnA[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodnB[ithread] = _bjmodnB[ithread] + ( (-(int)((uint32)_bjmodnB[ithread] >> 31)) & n);
			_bjmodnC[ithread] = _bjmodnB[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodnC[ithread] = _bjmodnC[ithread] + ( (-(int)((uint32)_bjmodnC[ithread] >> 31)) & n);
			_bjmodnD[ithread] = _bjmodnC[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodnD[ithread] = _bjmodnD[ithread] + ( (-(int)((uint32)_bjmodnD[ithread] >> 31)) & n);
			_bjmodnE[ithread] = _bjmodnD[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodnE[ithread] = _bjmodnE[ithread] + ( (-(int)((uint32)_bjmodnE[ithread] >> 31)) & n);
			_bjmodnF[ithread] = _bjmodnE[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodnF[ithread] = _bjmodnF[ithread] + ( (-(int)((uint32)_bjmodnF[ithread] >> 31)) & n);

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
	#if FFT_DEBUG
		fprintf(dbg_file, "radix16_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars addr & addp for this to compile properly: */
#ifdef MULTITHREAD
	omp_set_num_threads(CY_THREADS);
/*#undef PFETCH	*/
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,k1,k2,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,rt,it,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF) default(shared) schedule(static)
#else
/*                 empty_function(temp,frac,maxerr,i,j,j1,jstart,jhi,k,k1,k2,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,rt,it,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF);	*/
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
		}

		/* init carries	*/
	#if defined(USE_SSE2) && !defined(USE_SCALAR_CARRY)
		cy_r01->re = _cy_r0[ithread];	cy_i01->re = _cy_i0[ithread];
		cy_r01->im = _cy_r1[ithread];	cy_i01->im = _cy_i1[ithread];
		cy_r23->re = _cy_r2[ithread];	cy_i23->re = _cy_i2[ithread];
		cy_r23->im = _cy_r3[ithread];	cy_i23->im = _cy_i3[ithread];
		cy_r45->re = _cy_r4[ithread];	cy_i45->re = _cy_i4[ithread];
		cy_r45->im = _cy_r5[ithread];	cy_i45->im = _cy_i5[ithread];
		cy_r67->re = _cy_r6[ithread];	cy_i67->re = _cy_i6[ithread];
		cy_r67->im = _cy_r7[ithread];	cy_i67->im = _cy_i7[ithread];
		cy_r89->re = _cy_r8[ithread];	cy_i89->re = _cy_i8[ithread];
		cy_r89->im = _cy_r9[ithread];	cy_i89->im = _cy_i9[ithread];
		cy_rAB->re = _cy_rA[ithread];	cy_iAB->re = _cy_iA[ithread];
		cy_rAB->im = _cy_rB[ithread];	cy_iAB->im = _cy_iB[ithread];
		cy_rCD->re = _cy_rC[ithread];	cy_iCD->re = _cy_iC[ithread];
		cy_rCD->im = _cy_rD[ithread];	cy_iCD->im = _cy_iD[ithread];
		cy_rEF->re = _cy_rE[ithread];	cy_iEF->re = _cy_iE[ithread];
		cy_rEF->im = _cy_rF[ithread];	cy_iEF->im = _cy_iF[ithread];
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
			a[jt   ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp   ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p4;	jp = j2 + p4;
			a[jt   ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp   ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p8;	jp = j2 + p8;
			a[jt   ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp   ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p12;	jp = j2 + p12;
			a[jt   ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp   ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			jt = j1;		jp = j2;
			fprintf(stderr, "radix16_carry: A_in[0] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
			fprintf(stderr, "radix16_carry: A_in[1] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
			fprintf(stderr, "radix16_carry: A_in[2] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
			fprintf(stderr, "radix16_carry: A_in[3] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p4;	jp = j2 + p4;
			fprintf(stderr, "radix16_carry: A_in[4] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
			fprintf(stderr, "radix16_carry: A_in[5] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
			fprintf(stderr, "radix16_carry: A_in[6] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
			fprintf(stderr, "radix16_carry: A_in[7] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p8;	jp = j2 + p8;
			fprintf(stderr, "radix16_carry: A_in[8] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
			fprintf(stderr, "radix16_carry: A_in[9] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
			fprintf(stderr, "radix16_carry: A_in[A] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
			fprintf(stderr, "radix16_carry: A_in[B] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p12;	jp = j2 + p12;
			fprintf(stderr, "radix16_carry: A_in[C] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
			fprintf(stderr, "radix16_carry: A_in[D] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
			fprintf(stderr, "radix16_carry: A_in[E] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
			fprintf(stderr, "radix16_carry: A_in[F] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);
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

	#ifdef DEBUG_SSE2
		fprintf(stderr, "radix16_carry_ins: R1:= %20.5f, %20.5f\n",r1 ->re,r2 -> re);
		fprintf(stderr, "radix16_carry_ins: R3:= %20.5f, %20.5f\n",r3 ->re,r4 -> re);
		fprintf(stderr, "radix16_carry_ins: R5:= %20.5f, %20.5f\n",r5 ->re,r6 -> re);
		fprintf(stderr, "radix16_carry_ins: R7:= %20.5f, %20.5f\n",r7 ->re,r8 -> re);
		fprintf(stderr, "radix16_carry_ins: R9:= %20.5f, %20.5f\n",r9 ->re,r10-> re);
		fprintf(stderr, "radix16_carry_ins: R11= %20.5f, %20.5f\n",r11->re,r12-> re);
		fprintf(stderr, "radix16_carry_ins: R13= %20.5f, %20.5f\n",r13->re,r14-> re);
		fprintf(stderr, "radix16_carry_ins: R15= %20.5f, %20.5f\n",r15->re,r16-> re);
		fprintf(stderr, "radix16_carry_ins: R17= %20.5f, %20.5f\n",r17->re,r18-> re);
		fprintf(stderr, "radix16_carry_ins: R19= %20.5f, %20.5f\n",r19->re,r20-> re);
		fprintf(stderr, "radix16_carry_ins: R21= %20.5f, %20.5f\n",r21->re,r22-> re);
		fprintf(stderr, "radix16_carry_ins: R23= %20.5f, %20.5f\n",r23->re,r24-> re);
		fprintf(stderr, "radix16_carry_ins: R25= %20.5f, %20.5f\n",r25->re,r26-> re);
		fprintf(stderr, "radix16_carry_ins: R27= %20.5f, %20.5f\n",r27->re,r28-> re);
		fprintf(stderr, "radix16_carry_ins: R29= %20.5f, %20.5f\n",r29->re,r30-> re);
		fprintf(stderr, "radix16_carry_ins: R31= %20.5f, %20.5f\n",r31->re,r32-> re);
	#endif

			/* DEBUG: To do carries using the scalar [non-SSE2] carry macros, move appropriate half of the xmm temps into a1p's: */
			#ifdef USE_SCALAR_CARRY

			SSE_LOOP:

				if((j&3) == 0)
				{
					a1p0r=r1 ->re;	a1p0i=r2 ->re;
					a1p1r=r3 ->re;	a1p1i=r4 ->re;
					a1p2r=r5 ->re;	a1p2i=r6 ->re;
					a1p3r=r7 ->re;	a1p3i=r8 ->re;
					a1p4r=r9 ->re;	a1p4i=r10->re;
					a1p5r=r11->re;	a1p5i=r12->re;
					a1p6r=r13->re;	a1p6i=r14->re;
					a1p7r=r15->re;	a1p7i=r16->re;
					a1p8r=r17->re;	a1p8i=r18->re;
					a1p9r=r19->re;	a1p9i=r20->re;
					a1pAr=r21->re;	a1pAi=r22->re;
					a1pBr=r23->re;	a1pBi=r24->re;
					a1pCr=r25->re;	a1pCi=r26->re;
					a1pDr=r27->re;	a1pDi=r28->re;
					a1pEr=r29->re;	a1pEi=r30->re;
					a1pFr=r31->re;	a1pFi=r32->re;
				}
				else
				{
					a1p0r=r1 ->im;	a1p0i=r2 ->im;
					a1p1r=r3 ->im;	a1p1i=r4 ->im;
					a1p2r=r5 ->im;	a1p2i=r6 ->im;
					a1p3r=r7 ->im;	a1p3i=r8 ->im;
					a1p4r=r9 ->im;	a1p4i=r10->im;
					a1p5r=r11->im;	a1p5i=r12->im;
					a1p6r=r13->im;	a1p6i=r14->im;
					a1p7r=r15->im;	a1p7i=r16->im;
					a1p8r=r17->im;	a1p8i=r18->im;
					a1p9r=r19->im;	a1p9i=r20->im;
					a1pAr=r21->im;	a1pAi=r22->im;
					a1pBr=r23->im;	a1pBi=r24->im;
					a1pCr=r25->im;	a1pCi=r26->im;
					a1pDr=r27->im;	a1pDi=r28->im;
					a1pEr=r29->im;	a1pEi=r30->im;
					a1pFr=r31->im;	a1pFi=r32->im;
				}
			#endif

		#else	/* USE_SSE2 */

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

			#if defined(USE_SCALAR_CARRY) || !defined(USE_SSE2)

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

				#if defined(USE_SCALAR_CARRY)

				/*...set0 is slightly different from others:	*/
				   cmplx_carry_norm_pow2_errcheck0(a1p0r,a1p0i,cy_r0,*bjmodn0    );
					cmplx_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,*bjmodn1,0x1);
					cmplx_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,*bjmodn2,0x2);
					cmplx_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,*bjmodn3,0x3);
					cmplx_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,*bjmodn4,0x4);
					cmplx_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,*bjmodn5,0x5);
					cmplx_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,*bjmodn6,0x6);
					cmplx_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,*bjmodn7,0x7);
					cmplx_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,*bjmodn8,0x8);
					cmplx_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,*bjmodn9,0x9);
					cmplx_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,*bjmodnA,0xA);
					cmplx_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,*bjmodnB,0xB);
					cmplx_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,*bjmodnC,0xC);
					cmplx_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,*bjmodnD,0xD);
					cmplx_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,*bjmodnE,0xE);
					cmplx_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,*bjmodnF,0xF);

					i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/

				#else	/* SSE2 mode uses pointers for the bjmodn's, non-SSE2 uses scalars: */

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

				#endif

					co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
				}
				else
				{
					ntmp = 0;
				#if FFT_DEBUG
					DBG_fermat_carry_norm_pow2_errcheck(a1p0r,a1p0i,cy_r0,cy_i0,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,cy_i1,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,cy_i2,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,cy_i3,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,cy_i4,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,cy_i5,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,cy_i6,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,cy_i7,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,cy_i8,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,cy_i9,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,cy_iA,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,cy_iB,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,cy_iC,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,cy_iD,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,cy_iE,ntmp);	ntmp += n16;
					DBG_fermat_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,cy_iF,ntmp);
				#else
						fermat_carry_norm_pow2_errcheck(a1p0r,a1p0i,cy_r0,cy_i0,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,cy_i1,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,cy_i2,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,cy_i3,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,cy_i4,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,cy_i5,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,cy_i6,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,cy_i7,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,cy_i8,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,cy_i9,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,cy_iA,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,cy_iB,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,cy_iC,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,cy_iD,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,cy_iE,ntmp);	ntmp += n16;
						fermat_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,cy_iF,ntmp);
				#endif
				}

			#elif defined(USE_SSE2)	/* !defined(USE_SCALAR_CARRY) */

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
					__asm	mov		ecx, scale_ptr	// Master copy of wt0_ptr
					__asm	movlpd	xmm4,[ecx]
					__asm	movhpd	xmm4,[ecx]
					__asm	mulpd	xmm2,xmm4
					__asm	mulpd	xmm3,xmm4
					__asm	mov		edx, half_arr	// Master copy of wt0_ptr
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
						[__add3]			"m" (add3)\
						: "eax","ebx","ecx","edx","esi","edi","xmm0","xmm1","xmm2","xmm3","xmm4"	/* Clobbered registers */\
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

					#endif

				#endif
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
//
if(j < 2 && full_pass && iter <= 10)
{
	fprintf(stderr, "Iter %3d: a0-3_in = %20.5f %20.5f %20.5f %20.5f, cy0-3 = %20.5f %20.5f %20.5f %20.5f, maxerr = %10.8f\n",iter,r1 ->re,r2 ->re,r1 ->im,r2 ->im,cy_r01->re,cy_r01->im,cy_r23->re,cy_r23->im,maxerr);
}
*/
			#if defined(COMPILER_TYPE_MSVC)

				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B (r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B (r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B (r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC);

			#else	/* GCC-style inline ASM: */

				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r1 ,add1,add2,add3,cy_r01,cy_r23,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r9 ,add1,add2,add3,cy_r45,cy_r67,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r17,add1,add2,add3,cy_r89,cy_rAB,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r25,add1,add2,add3,cy_rCD,cy_rEF,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);

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
					__asm	mov		ecx, scale_ptr	// Master copy of wt0_ptr
					__asm	movlpd	xmm4,[ecx]
					__asm	movhpd	xmm4,[ecx]
					__asm	mulpd	xmm2,xmm4
					__asm	mulpd	xmm3,xmm4
					__asm	mov		edx, half_arr	// Master copy of wt0_ptr
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
					  [__add3]			"m" (add3)\
					  : "eax","ebx","ecx","edx","esi","edi","xmm0","xmm1","xmm2","xmm3","xmm4"	/* Clobbered registers */\
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

				SSE2_cmplx_carry_norm_pow2_errcheck2_2B (r1 ,add1,add2,cy_r01,cy_r23,bjmodn0);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B (r9 ,add1,add2,cy_r45,cy_r67,bjmodn4);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B (r17,add1,add2,cy_r89,cy_rAB,bjmodn8);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B (r25,add1,add2,cy_rCD,cy_rEF,bjmodnC);

			#else	/* GCC-style inline ASM: */

				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r1 ,add1,add2,cy_r01,cy_r23,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r9 ,add1,add2,cy_r45,cy_r67,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r17,add1,add2,cy_r89,cy_rAB,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r25,add1,add2,cy_rCD,cy_rEF,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);

			#endif

	#ifdef DEBUG_SSE2
		fprintf(stderr, "radix16_carry_out: R1:= %20.5f, %20.5f\n",r1 ->re,r2 -> re);
		fprintf(stderr, "radix16_carry_out: R3:= %20.5f, %20.5f\n",r3 ->re,r4 -> re);
		fprintf(stderr, "radix16_carry_out: R5:= %20.5f, %20.5f\n",r5 ->re,r6 -> re);
		fprintf(stderr, "radix16_carry_out: R7:= %20.5f, %20.5f\n",r7 ->re,r8 -> re);
		fprintf(stderr, "radix16_carry_out: R9:= %20.5f, %20.5f\n",r9 ->re,r10-> re);
		fprintf(stderr, "radix16_carry_out: R11= %20.5f, %20.5f\n",r11->re,r12-> re);
		fprintf(stderr, "radix16_carry_out: R13= %20.5f, %20.5f\n",r13->re,r14-> re);
		fprintf(stderr, "radix16_carry_out: R15= %20.5f, %20.5f\n",r15->re,r16-> re);
		fprintf(stderr, "radix16_carry_out: R17= %20.5f, %20.5f\n",r17->re,r18-> re);
		fprintf(stderr, "radix16_carry_out: R19= %20.5f, %20.5f\n",r19->re,r20-> re);
		fprintf(stderr, "radix16_carry_out: R21= %20.5f, %20.5f\n",r21->re,r22-> re);
		fprintf(stderr, "radix16_carry_out: R23= %20.5f, %20.5f\n",r23->re,r24-> re);
		fprintf(stderr, "radix16_carry_out: R25= %20.5f, %20.5f\n",r25->re,r26-> re);
		fprintf(stderr, "radix16_carry_out: R27= %20.5f, %20.5f\n",r27->re,r28-> re);
		fprintf(stderr, "radix16_carry_out: R29= %20.5f, %20.5f\n",r29->re,r30-> re);
		fprintf(stderr, "radix16_carry_out: R31= %20.5f, %20.5f\n",r31->re,r32-> re);
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
					                                    a[n/2] and a[n/2+1] of the right-angle transform.

				Because of the indesirable intra-xmm-register data dependency this leads to, we instead require data arranged as

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

				Note that even though e.g. a0 and a1 appear adjacent, they are actually n/16 memory locations apart, i.e. there
				is no carry propagation between them.
				*/
				}	/* if(MODULUS_TYPE == ...) */

			#endif	/* USE_SSE2 */

				/*...The radix-16 DIF pass is here:	*/

	#ifdef USE_SSE2

		#ifdef USE_SCALAR_CARRY

			if((j&3) == 0)
			{
				r1 ->re=a1p0r;	r2 ->re=a1p0i;
				r17->re=a1p8r;	r18->re=a1p8i;

				r9 ->re=a1p4r;	r10->re=a1p4i;
				r25->re=a1pCr;	r26->re=a1pCi;

				r5 ->re=a1p2r;	r6 ->re=a1p2i;
				r21->re=a1pAr;	r22->re=a1pAi;

				r13->re=a1p6r;	r14->re=a1p6i;
				r29->re=a1pEr;	r30->re=a1pEi;

				r3 ->re=a1p1r;	r4 ->re=a1p1i;
				r19->re=a1p9r;	r20->re=a1p9i;

				r11->re=a1p5r;	r12->re=a1p5i;
				r27->re=a1pDr;	r28->re=a1pDi;

				r7 ->re=a1p3r;	r8 ->re=a1p3i;
				r23->re=a1pBr;	r24->re=a1pBi;

				r15->re=a1p7r;	r16->re=a1p7i;
				r31->re=a1pFr;	r32->re=a1pFi;

				j += 2;
				goto SSE_LOOP;	/* Go back and do 2nd set of carries */
			}
			else
			{
				r1 ->im=a1p0r;	r2 ->im=a1p0i;
				r17->im=a1p8r;	r18->im=a1p8i;

				r9 ->im=a1p4r;	r10->im=a1p4i;
				r25->im=a1pCr;	r26->im=a1pCi;

				r5 ->im=a1p2r;	r6 ->im=a1p2i;
				r21->im=a1pAr;	r22->im=a1pAi;

				r13->im=a1p6r;	r14->im=a1p6i;
				r29->im=a1pEr;	r30->im=a1pEi;

				r3 ->im=a1p1r;	r4 ->im=a1p1i;
				r19->im=a1p9r;	r20->im=a1p9i;

				r11->im=a1p5r;	r12->im=a1p5i;
				r27->im=a1pDr;	r28->im=a1pDi;

				r7 ->im=a1p3r;	r8 ->im=a1p3i;
				r23->im=a1pBr;	r24->im=a1pBi;

				r15->im=a1p7r;	r16->im=a1p7i;
				r31->im=a1pFr;	r32->im=a1pFi;

				j -= 2;
			}

		#endif

		/* Four DIF radix-4 subconvolution, sans twiddles.	Cost each: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */

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

	#ifdef DEBUG_SSE2
		jt = j1;		jp = j2;
		fprintf(stderr, "radix16_carry: Aout[0] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
		fprintf(stderr, "radix16_carry: Aout[1] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
		fprintf(stderr, "radix16_carry: Aout[2] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
		fprintf(stderr, "radix16_carry: Aout[3] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p4;	jp = j2 + p4;
		fprintf(stderr, "radix16_carry: Aout[4] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
		fprintf(stderr, "radix16_carry: Aout[5] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
		fprintf(stderr, "radix16_carry: Aout[6] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
		fprintf(stderr, "radix16_carry: Aout[7] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p8;	jp = j2 + p8;
		fprintf(stderr, "radix16_carry: Aout[8] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
		fprintf(stderr, "radix16_carry: Aout[9] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
		fprintf(stderr, "radix16_carry: Aout[A] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
		fprintf(stderr, "radix16_carry: Aout[B] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p12;	jp = j2 + p12;
		fprintf(stderr, "radix16_carry: Aout[C] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
		fprintf(stderr, "radix16_carry: Aout[D] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
		fprintf(stderr, "radix16_carry: Aout[E] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
		fprintf(stderr, "radix16_carry: Aout[F] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);
		exit(0);
	#endif

	#else	/* USE_SSE2 */
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
	#if defined(USE_SSE2) && !defined(USE_SCALAR_CARRY)
		_cy_r0[ithread] = cy_r01->re;	_cy_i0[ithread] = cy_i01->re;
		_cy_r1[ithread] = cy_r01->im;	_cy_i1[ithread] = cy_i01->im;
		_cy_r2[ithread] = cy_r23->re;	_cy_i2[ithread] = cy_i23->re;
		_cy_r3[ithread] = cy_r23->im;	_cy_i3[ithread] = cy_i23->im;
		_cy_r4[ithread] = cy_r45->re;	_cy_i4[ithread] = cy_i45->re;
		_cy_r5[ithread] = cy_r45->im;	_cy_i5[ithread] = cy_i45->im;
		_cy_r6[ithread] = cy_r67->re;	_cy_i6[ithread] = cy_i67->re;
		_cy_r7[ithread] = cy_r67->im;	_cy_i7[ithread] = cy_i67->im;
		_cy_r8[ithread] = cy_r89->re;	_cy_i8[ithread] = cy_i89->re;
		_cy_r9[ithread] = cy_r89->im;	_cy_i9[ithread] = cy_i89->im;
		_cy_rA[ithread] = cy_rAB->re;	_cy_iA[ithread] = cy_iAB->re;
		_cy_rB[ithread] = cy_rAB->im;	_cy_iB[ithread] = cy_iAB->im;
		_cy_rC[ithread] = cy_rCD->re;	_cy_iC[ithread] = cy_iCD->re;
		_cy_rD[ithread] = cy_rCD->im;	_cy_iD[ithread] = cy_iCD->im;
		_cy_rE[ithread] = cy_rEF->re;	_cy_iE[ithread] = cy_iEF->re;
		_cy_rF[ithread] = cy_rEF->im;	_cy_iF[ithread] = cy_iEF->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
	#else
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
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

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
			ASSERT(HERE, CY_THREADS > 1,"radix16_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
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
			ASSERT(HERE, CY_THREADS > 1,"radix16_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
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

#if FFT_DEBUG
	fclose(dbg_file);
	dbg_file = 0x0;
	sprintf(cbuf, "Wrote debug file %s", dbg_fname);
	if(iter == 20)
		ASSERT(HERE, 0, cbuf);
#endif

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

#if FFT_DEBUG
	int len_a;
	FILE *dbg_file;
	const char dbg_fname[] = "CY16_MT_DEBUG.txt";
	ASSERT(HERE, n < 10000, "n too large!");
	dbg_file = fopen(dbg_fname, "a");
	fprintf(dbg_file, "radix16_ditN_cy_dif1_nochk: Multithreaded, n = %d\n", n);
#endif

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
	#if FFT_DEBUG
		fprintf(dbg_file, "radix16_ditn_cy_dif1_nochk: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars addr & addp for this to compile properly: */
#ifdef MULTITHREAD
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
					fermat_carry_norm_pow2_nocheck(a1p0r,a1p0i,cy_r0,cy_i0,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p1r,a1p1i,cy_r1,cy_i1,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p2r,a1p2i,cy_r2,cy_i2,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p3r,a1p3i,cy_r3,cy_i3,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p4r,a1p4i,cy_r4,cy_i4,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p5r,a1p5i,cy_r5,cy_i5,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p6r,a1p6i,cy_r6,cy_i6,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p7r,a1p7i,cy_r7,cy_i7,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p8r,a1p8i,cy_r8,cy_i8,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1p9r,a1p9i,cy_r9,cy_i9,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pAr,a1pAi,cy_rA,cy_iA,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pBr,a1pBi,cy_rB,cy_iB,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pCr,a1pCi,cy_rC,cy_iC,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pDr,a1pDi,cy_rD,cy_iD,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pEr,a1pEi,cy_rE,cy_iE,ntmp);	ntmp += n16;
					fermat_carry_norm_pow2_nocheck(a1pFr,a1pFi,cy_rF,cy_iF,ntmp);
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
			ASSERT(HERE, CY_THREADS > 1,"radix16_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
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
			ASSERT(HERE, CY_THREADS > 1,"radix16_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
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

