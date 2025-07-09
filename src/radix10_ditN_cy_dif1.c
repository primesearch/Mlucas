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

#undef MULTITHREAD	/* EWM - no thread support here for now */

/***************/

int radix10_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-10 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-10 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n10,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9
		,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double	cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					sn2 =  0.95105651629515357211,	/*  sin(u) */
					ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					ss2 =  0.36327126400268044292,	/* [sin(u)-sin(2u)] */
					radix_inv, n2inv;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,temp,frac,scale;
	double maxerr = 0.0;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;
	double wt_re,wt_im, wi_re,wi_im;					/* Fermat/LOACC weights stuff */

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;

	if(RES_SHIFT) {
		WARN(HERE, "CY routines with radix < 16 do not support shifted residues!", "", 1);
		return(ERR_ASSERT);
	}

	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	if((TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change n10 and n_div_wt to non-static to work around a gcc compiler bug. */
	n10   = n/10;
	n_div_nwt = n10 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n10)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/10 in radix10_ditN_cy_dif1.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave || n != nsave) {	/* Exponent or array length change triggers re-init */
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;	nsave = n;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)10));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p1 = n10;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );

		bjmodnini=0;
		for(j=0; j < n10; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-10 final DIT pass is here.	*/

	cy0 = 0;		/* init carries	*/
	cy1 = 0;
	cy2 = 0;
	cy3 = 0;
	cy4 = 0;
	cy5 = 0;
	cy6 = 0;
	cy7 = 0;
	cy8 = 0;
	cy9 = 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy0 = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/

	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

	jstart = 0;
	jhi = jstart+nwt-1;
	khi = n_div_nwt;

for(outer=0; outer <= 1; outer++)
{
	i = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (i = 1).	*/

	bjmodn0 = 0;
	bjmodn1 = bjmodnini;
	bjmodn2 = bjmodn1 +bjmodnini-n; bjmodn2 = bjmodn2 + ( (-(int)((uint32)bjmodn2 >> 31)) & n);
	bjmodn3 = bjmodn2 +bjmodnini-n; bjmodn3 = bjmodn3 + ( (-(int)((uint32)bjmodn3 >> 31)) & n);
	bjmodn4 = bjmodn3 +bjmodnini-n; bjmodn4 = bjmodn4 + ( (-(int)((uint32)bjmodn4 >> 31)) & n);
	bjmodn5 = bjmodn4 +bjmodnini-n; bjmodn5 = bjmodn5 + ( (-(int)((uint32)bjmodn5 >> 31)) & n);
	bjmodn6 = bjmodn5 +bjmodnini-n; bjmodn6 = bjmodn6 + ( (-(int)((uint32)bjmodn6 >> 31)) & n);
	bjmodn7 = bjmodn6 +bjmodnini-n; bjmodn7 = bjmodn7 + ( (-(int)((uint32)bjmodn7 >> 31)) & n);
	bjmodn8 = bjmodn7 +bjmodnini-n; bjmodn8 = bjmodn8 + ( (-(int)((uint32)bjmodn8 >> 31)) & n);
	bjmodn9 = bjmodn8 +bjmodnini-n; bjmodn9 = bjmodn9 + ( (-(int)((uint32)bjmodn9 >> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+10;
	co3=co2-10;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

			t1 =a[j1   ];		t2 =a[j2   ];
			rt =a[j1+p1];		it =a[j2+p1];
			t3 =t1 -rt;			t4 =t2 -it;
			t1 =t1 +rt;			t2 =t2 +it;

			t5 =a[j1+p7];		t6 =a[j2+p7];
			rt =a[j1+p6];		it =a[j2+p6];
			t7 =t5 -rt;			t8 =t6 -it;
			t5 =t5 +rt;			t6 =t6 +it;

			t9 =a[j1+p3];		t10=a[j2+p3];
			rt =a[j1+p2];		it =a[j2+p2];
			t11=t9 -rt;			t12=t10-it;
			t9 =t9 +rt;			t10=t10+it;

			t13=a[j1+p8];		t14=a[j2+p8];
			rt =a[j1+p9];		it =a[j2+p9];
			t15=t13-rt;			t16=t14-it;
			t13=t13+rt;			t14=t14+it;

			t17=a[j1+p4];		t18=a[j2+p4];
			rt =a[j1+p5];		it =a[j2+p5];
			t19=t17-rt;			t20=t18-it;
			t17=t17+rt;			t18=t18+it;

			rt =t17;			it =t18;
			t17=t5 -rt;			t18=t6 -it;
			t5 =t5 +rt;			t6 =t6 +it;
			rt =t13;			it =t14;
			t13=t9 -rt;			t14=t10 -it;
			t9 =t9 +rt;			t10=t10 +it;
			rt = t5+t9;			it = t6+t10;
			t1 = t1+rt;			t2 = t2+it;
			rt = t1+cc1*rt;		it = t2+cc1*it;
			t9 = cc2*(t5-t9);	t10= cc2*(t6-t10);
			t5 = rt+t9;			t6 = it+t10;
			t9 = rt-t9;			t10= it-t10;
			rt = sn2*(t13-t17);	it = sn2*(t14-t18);
			t13= ss1* t13;		t14= ss1* t14;
			t17= ss2* t17;		t18= ss2* t18;
			t13= rt-t13;		t14= it-t14;
			t17= rt+t17;		t18= it+t18;
			a1p0r  =t1;			a1p0i  =t2;
			a1p4r =t5-t14;		a1p4i =t6 +t13;
			a1p8r =t9-t18;		a1p8i =t10+t17;
			a1p2r =t9+t18;		a1p2i =t10-t17;
			a1p6r =t5+t14;		a1p6i =t6 -t13;

			rt =t19;			it =t20;
			t19=t7 -rt;			t20=t8 -it;
			t7 =t7 +rt;			t8 =t8 +it;
			rt =t15;			it =t16;
			t15=t11 -rt;		t16=t12 -it;
			t11=t11 +rt;		t12=t12 +it;
			rt = t7+t11;		it = t8+t12;
			t3 = t3+rt;			t4 = t4+it;
			rt = t3+cc1*rt;		it = t4+cc1*it;
			t11= cc2*(t7-t11);	t12= cc2*(t8-t12);
			t7 = rt+t11;		t8 = it+t12;
			t11= rt-t11;		t12= it-t12;
			rt = sn2*(t15-t19);	it = sn2*(t16-t20);
			t15= ss1* t15;		t16= ss1* t16;
			t19= ss2* t19;		t20= ss2* t20;
			t15= rt-t15;		t16= it-t16;
			t19= rt+t19;		t20= it+t20;
			a1p5r =t3;			a1p5i =t4;
			a1p9r =t7 -t16;		a1p9i =t8 +t15;
			a1p3r =t11-t20;		a1p3i =t12+t19;
			a1p7r =t11+t20;		a1p7i =t12-t19;
			a1p1r =t7 +t16;		a1p1i =t8 -t15;

			/*...Now do the carries. Since the outputs would
			normally be getting dispatched to 10 separate blocks of the A-array, we need 10 separate carries.
			*/
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
			cmplx_carry_norm_errcheck0(a1p0r ,a1p0i ,cy0 ,bjmodn0 ,0 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p1r ,a1p1i ,cy1 ,bjmodn1 ,1 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p2r ,a1p2i ,cy2 ,bjmodn2 ,2 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p3r ,a1p3i ,cy3 ,bjmodn3 ,3 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p4r ,a1p4i ,cy4 ,bjmodn4 ,4 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p5r ,a1p5i ,cy5 ,bjmodn5 ,5 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p6r ,a1p6i ,cy6 ,bjmodn6 ,6 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p7r ,a1p7i ,cy7 ,bjmodn7 ,7 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p8r ,a1p8i ,cy8 ,bjmodn8 ,8 ,prp_mult);
			cmplx_carry_norm_errcheck (a1p9r ,a1p9i ,cy9 ,bjmodn9 ,9 ,prp_mult);

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 		and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			/*...The radix-10 DIF pass is here:	*/
			/*
			Twiddleless version requires us to swap inputs x2 <-> x8, x4 <-> x6, x1 <-> x5 and x7 <-> x9...
			indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9
				  -> 0,-5, 8, 3, 6, 1, 4,-1, 2,-3
				  == 0, 5, 8, 3, 6, 1, 4, 9, 2, 7 modulo 10.
			I.e. start out with first quintet of indices {0,2,4,6,8}, permute those according to
			{0,2,4,6,8}*9%10 = {0,8,6,4,2}, then each is head of a length-2 list of indices with decrement 5.
			Outputs are permuted as with dif_pass1.
			*/
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif
			t1 =a1p0r;				t2 =a1p0i;
			t3 =a1p8r+a1p2r;		t4 =a1p8i+a1p2i;
			t5 =a1p6r+a1p4r;		t6 =a1p6i+a1p4i;
			t7 =a1p6r-a1p4r;		t8 =a1p6i-a1p4i;
			t9 =a1p8r-a1p2r;		t10=a1p8i-a1p2i;
			rt = t3+t5;				it = t4+t6;
			t1 = t1+rt;				t2 = t2+it;		/* y0	*/

			rt = t1+cc1*rt;			it = t2+cc1*it;
			t5 = cc2*(t3-t5);		t6 = cc2*(t4-t6);
			t3 = rt+t5;				t4 = it+t6;
			t5 = rt-t5;				t6 = it-t6;
			rt = sn2*(t9-t7);		it = sn2*(t10-t8);
			t7 = ss1* t7;			t8 = ss1* t8;
			t9 = ss2* t9;			t10= ss2* t10;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
			t7 = rt+t7;				t8 = it+t8;
			t9 = rt-t9;				t10= it-t10;
			rt = t7;				it = t8;
			t7 = t3+it;				t8 = t4-rt;		/* y4	*/
			t3 = t3-it;				t4 = t4+rt;		/* y1	*/
			rt = t9;				it = t10;
			t9 = t5+it;				t10= t6-rt;		/* y3	*/
			t5 = t5-it;				t6 = t6+rt;		/* y2	*/
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif

			t11 =a1p5r;				t12 =a1p5i;
			t13 =a1p3r+a1p7r;		t14 =a1p3i+a1p7i;
			t15 =a1p1r+a1p9r;		t16 =a1p1i+a1p9i;
			t17 =a1p1r-a1p9r;		t18 =a1p1i-a1p9i;
			t19 =a1p3r-a1p7r;		t20 =a1p3i-a1p7i;
			rt  = t13+t15;			it  = t14+t16;
			t11 = t11+rt;			t12 = t12+it;		/* z0	*/
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
			rt  = t11+cc1*rt;		it  = t12+cc1*it;
			t15 = cc2*(t13-t15);		t16 = cc2*(t14-t16);
			t13 = rt+t15;			t14 = it+t16;
			t15 = rt-t15;			t16 = it-t16;
			rt  = sn2*(t19-t17);		it  = sn2*(t20-t18);
			t17 = ss1* t17;			t18 = ss1* t18;
			t19 = ss2* t19;			t20 = ss2* t20;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
			t17 = rt+t17;			t18 = it+t18;
			t19 = rt-t19;			t20 = it-t20;
			rt  = t17;				it  = t18;
			t17 = t13+it;			t18 = t14-rt;		/* z4	*/
			t13 = t13-it;			t14 = t14+rt;		/* z1	*/
			rt  = t19;				it  = t20;
			t19 = t15+it;			t20 = t16-rt;		/* z3	*/
			t15 = t15-it;			t16 = t16+rt;		/* z2	*/

			/* ...and now do five radix-2 transforms.	*/
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
			a[j1   ]=t1+t11;		a[j2   ]=t2+t12;
			a[j1+p1]=t1-t11;		a[j2+p1]=t2-t12;
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif
			a[j1+p8]=t3 +t13;		a[j2+p8]=t4 +t14;
			a[j1+p9]=t3 -t13;		a[j2+p9]=t4 -t14;
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
			a[j1+p7]=t5 +t15;		a[j2+p7]=t6 +t16;
			a[j1+p6]=t5 -t15;		a[j2+p6]=t6 -t16;
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif
			a[j1+p4]=t9 +t19;		a[j2+p4]=t10+t20;
			a[j1+p5]=t9 -t19;		a[j2+p5]=t10-t20;
#if PFETCH
addr = add0+p9;
prefetch_p_doubles(addr);
#endif
			a[j1+p3]=t7 +t17;		a[j2+p3]=t8 +t18;
			a[j1+p2]=t7 -t17;		a[j2+p2]=t8 -t18;
		}

		jstart += nwt;
		jhi    += nwt;
		col += 10;
		co3 -= 10;

	}

	if(root_incr==0)break;

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-10 forward DIF FFT of the first block of 10 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 10 outputs of (1);
	!   (3) Reweight and perform a radix-10 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 10 elements and repeat (1-4).
	*/
	/*
	printf("iter = %10d; maxerr = %20.15f\n",iter,maxerr);
	printf("carries = %10d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",(int)cy0,(int)cy1,(int)cy2,(int)cy3,(int)cy4,(int)cy5,(int)cy6,(int)cy7,(int)cy8,(int)cy9);
	*/
	t1  = cy9;
	cy9 = cy8;
	cy8 = cy7;
	cy7 = cy6;
	cy6 = cy5;
	cy5 = cy4;
	cy4 = cy3;
	cy3 = cy2;
	cy2 = cy1;
	cy1 = cy0;
	cy0 = t1;

	root_incr = 0;
	scale = prp_mult = 1;

	jstart = 0;
	jhi = 7;
	khi = 1;

	for(j=0; j<=jhi; j++)
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
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9) != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix10_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix10_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-10 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Uses an optimized radix-5 transform a la Nussbaumer (2nd ed., p.146).
*/
	int j,j1,j2;
	static int n10,p1,p2,p3,p4,p5,p6,p7,p8,p9, first_entry=TRUE;
	static double	cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					s2  =  0.95105651629515357211,	/*  sin(u) */
					ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20;

	if(!first_entry && (n/10) != n10)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n10=n/10;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n10;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-10 pass is here.	*/

	for(j=0; j < n10; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		/*       gather the needed data (10 64-bit complex, i.e. 20 64-bit reals) and do two radix-5 transforms...	*/

		/*
		Twiddleless version requires us to swap inputs x2 <-> x8, x4 <-> x6, x1 <-> x5 and x7 <-> x9...
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9
			  -> 0,-5, 8, 3, 6, 1, 4,-1, 2,-3
			  == 0, 5, 8, 3, 6, 1, 4, 9, 2, 7 modulo 10.
		I.e. start out with first quintet of indices {0,2,4,6,8}, permute those according to
		{0,2,4,6,8}*9%10 = {0,8,6,4,2}, then each is head of a length-2 list of indices with decrement 5.
		*/
		t1 =a[j1   ];			t2 =a[j2   ];
		t3 =a[j1+p8]+a[j1+p2];	t4 =a[j2+p8]+a[j2+p2];
		t5 =a[j1+p6]+a[j1+p4];	t6 =a[j2+p6]+a[j2+p4];
		t7 =a[j1+p6]-a[j1+p4];	t8 =a[j2+p6]-a[j2+p4];
		t9 =a[j1+p8]-a[j1+p2];	t10=a[j2+p8]-a[j2+p2];
		rt = t3+t5;				it = t4+t6;
		t1 = t1+rt;				t2 = t2+it;		/* y0	*/
		rt = t1+cc1*rt;			it = t2+cc1*it;
		t5 = cc2*(t3-t5);		t6 = cc2*(t4-t6);
		t3 = rt+t5;				t4 = it+t6;
		t5 = rt-t5;				t6 = it-t6;
		rt = s2 *(t9-t7);		it = s2 *(t10-t8);
		t7 = ss1* t7;			t8 = ss1* t8;
		t9 = ss2* t9;			t10= ss2* t10;
		t7 = rt+t7;				t8 = it+t8;
		t9 = rt-t9;				t10= it-t10;
		rt = t7;				it = t8;
		t7 = t3+it;				t8 = t4-rt;		/* y4 - note swap with y3 below!	*/
		t3 = t3-it;				t4 = t4+rt;		/* y1	*/
		rt = t9;				it = t10;
		t9 = t5+it;				t10= t6-rt;		/* y3	*/
		t5 = t5-it;				t6 = t6+rt;		/* y2	*/

		t11 =a[j1+p5];			t12 =a[j2+p5];
		t13 =a[j1+p3]+a[j1+p7];	t14 =a[j2+p3]+a[j2+p7];
		t15 =a[j1+p1]+a[j1+p9];	t16 =a[j2+p1]+a[j2+p9];
		t17 =a[j1+p1]-a[j1+p9];	t18 =a[j2+p1]-a[j2+p9];
		t19 =a[j1+p3]-a[j1+p7];	t20 =a[j2+p3]-a[j2+p7];
		rt  = t13+t15;			it  = t14+t16;
		t11 = t11+rt;			t12 = t12+it;	/* z0	*/
		rt  = t11+cc1*rt;		it  = t12+cc1*it;
		t15 = cc2*(t13-t15);	t16 = cc2*(t14-t16);
		t13 = rt+t15;			t14 = it+t16;
		t15 = rt-t15;			t16 = it-t16;
		rt  = s2 *(t19-t17);	it  = s2 *(t20-t18);
		t17 = ss1* t17;			t18 = ss1* t18;
		t19 = ss2* t19;			t20 = ss2* t20;
		t17 = rt+t17;			t18 = it+t18;
		t19 = rt-t19;			t20 = it-t20;
		rt  = t17;				it  = t18;
		t17 = t13+it;			t18 = t14-rt;		/* z4 - note swap with z3 below!	*/
		t13 = t13-it;			t14 = t14+rt;		/* z1	*/
		rt  = t19;				it  = t20;
		t19 = t15+it;			t20 = t16-rt;		/* z3	*/
		t15 = t15-it;			t16 = t16+rt;		/* z2	*/

/*       ...and now do five radix-2 transforms: */
		a[j1   ]=t1 +t11;		a[j2   ]=t2 +t12;
		a[j1+p1]=t1 -t11;		a[j2+p1]=t2 -t12;

		a[j1+p8]=t3 +t13;		a[j2+p8]=t4 +t14;
		a[j1+p9]=t3 -t13;		a[j2+p9]=t4 -t14;

		a[j1+p7]=t5 +t15;		a[j2+p7]=t6 +t16;
		a[j1+p6]=t5 -t15;		a[j2+p6]=t6 -t16;

		a[j1+p4]=t9 +t19;		a[j2+p4]=t10+t20;
		a[j1+p5]=t9 -t19;		a[j2+p5]=t10-t20;

		a[j1+p3]=t7 +t17;		a[j2+p3]=t8 +t18;
		a[j1+p2]=t7 -t17;		a[j2+p2]=t8 -t18;
								/* Totals: 17*4+5*4 = 88 FADD, 4*5 = 20 FMUL.	*/
	}
}

/***************/

void radix10_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-5 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix10_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
!
!   See the documentation in radix16_dif_pass for further details on the array indexing,
!   and radix10_dif_pass1 for details on the algorithm.
!
!   See the documentation in radix5_dif_pass1 for notes on the algorithm.
*/
	int j,j1,j2;
	static int n10,p1,p2,p3,p4,p5,p6,p7,p8,p9, first_entry=TRUE;
	static double	cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					s2  =  0.95105651629515357211,	/*  sin(u) */
					ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20;

	if(!first_entry && (n/10) != n10)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n10=n/10;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n10;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-10 pass is here.	*/

	for(j=0; j < n10; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		/*       gather the needed data (5 64-bit complex, i.e. 10 64-bit reals) and do five radix-2 transforms.
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9
			  -> 0,-2,-4,-6,-8, 5, 3, 1,-1,-3
			  == 0, 8, 6, 4, 2, 5, 3, 1, 9, 7 modulo 10.
		I.e. start out with first pair of indices {0,5}, permute those according to
		{0,5}*9%10 = {0,5}, then each is head of a length-5 list of indices with decrement 2.

		Remember, inputs to DIT are bit-reversed, so a[0,1,2,3,4,5,6,7,8,9]
											 contain x[0,5,1,6,2,7,3,8,4,9], so
		x1->x8	means a2 -> a7
		x2->x6	means a4 -> a3
		x3->x4	means a6 -> a8
		x4->x2	means a8 -> a4
		x6->x3	means a3 -> a6
		x7->x1	means a5 -> a2
		x8->x9	means a7 -> a9
		x9->x7	means a9 -> a5
		*/
		t1 =a[j1   ];		t2 =a[j2   ];
		rt =a[j1+p1];		it =a[j2+p1];
		t3 =t1 -rt;			t4 =t2 -it;
		t1 =t1 +rt;			t2 =t2 +it;

		t5 =a[j1+p7];		t6 =a[j2+p7];
		rt =a[j1+p6];		it =a[j2+p6];
		t7 =t5 -rt;			t8 =t6 -it;
		t5 =t5 +rt;			t6 =t6 +it;

		t9 =a[j1+p3];		t10=a[j2+p3];
		rt =a[j1+p2];		it =a[j2+p2];
		t11=t9 -rt;			t12=t10-it;
		t9 =t9 +rt;			t10=t10+it;

		t13=a[j1+p8];		t14=a[j2+p8];
		rt =a[j1+p9];		it =a[j2+p9];
		t15=t13-rt;			t16=t14-it;
		t13=t13+rt;			t14=t14+it;

		t17=a[j1+p4];		t18=a[j2+p4];
		rt =a[j1+p5];		it =a[j2+p5];
		t19=t17-rt;			t20=t18-it;
		t17=t17+rt;			t18=t18+it;

		/* ...and now do two radix-5 transforms.	*/

		rt =t17;			it =t18;
		t17=t5 -rt;			t18=t6 -it;
		t5 =t5 +rt;			t6 =t6 +it;
		rt =t13;			it =t14;
		t13=t9 -rt;			t14=t10 -it;
		t9 =t9 +rt;			t10=t10 +it;
		rt = t5+t9;			it = t6+t10;
		t1 = t1+rt;			t2 = t2+it;
		rt = t1+cc1*rt;		it = t2+cc1*it;
		t9 = cc2*(t5-t9);	t10= cc2*(t6-t10);
		t5 = rt+t9;			t6 = it+t10;
		t9 = rt-t9;			t10= it-t10;
		rt = s2 *(t13-t17);	it = s2 *(t14-t18);
		t13= ss1* t13;		t14= ss1* t14;
		t17= ss2* t17;		t18= ss2* t18;
		t13= rt-t13;		t14= it-t14;
		t17= rt+t17;		t18= it+t18;
		a[j1   ]=t1;		a[j2   ]=t2;
		a[j1+p4]=t5-t14;	a[j2+p4]=t6 +t13;
		a[j1+p8]=t9-t18;	a[j2+p8]=t10+t17;
		a[j1+p2]=t9+t18;	a[j2+p2]=t10-t17;
		a[j1+p6]=t5+t14;	a[j2+p6]=t6 -t13;

		rt =t19;			it =t20;
		t19=t7 -rt;			t20=t8 -it;
		t7 =t7 +rt;			t8 =t8 +it;
		rt =t15;			it =t16;
		t15=t11 -rt;		t16=t12 -it;
		t11=t11 +rt;		t12=t12 +it;
		rt = t7+t11;		it = t8+t12;
		t3 = t3+rt;			t4 = t4+it;
		rt = t3+cc1*rt;		it = t4+cc1*it;
		t11= cc2*(t7-t11);	t12= cc2*(t8-t12);
		t7 = rt+t11;		t8 = it+t12;
		t11= rt-t11;		t12= it-t12;
		rt = s2 *(t15-t19);	it = s2 *(t16-t20);
		t15= ss1* t15;		t16= ss1* t16;
		t19= ss2* t19;		t20= ss2* t20;
		t15= rt-t15;		t16= it-t16;
		t19= rt+t19;		t20= it+t20;
		a[j1+p5]=t3;		a[j2+p5]=t4;
		a[j1+p9]=t7 -t16;	a[j2+p9]=t8 +t15;
		a[j1+p3]=t11-t20;	a[j2+p3]=t12+t19;
		a[j1+p7]=t11+t20;	a[j2+p7]=t12-t19;
		a[j1+p1]=t7 +t16;	a[j2+p1]=t8 -t15;
	}
}

