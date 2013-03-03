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
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
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

/*...change n10 and n_div_wt to non-static to work around a gcc compiler bug. */
	n10   = n/10;
	n_div_nwt = n10 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n10)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/10 in radix10_ditN_cy_dif1.\n",iter);
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

	cy0 =-2;		/* init carries	*/
	cy1 = 0;
	cy2 = 0;
	cy3 = 0;
	cy4 = 0;
	cy5 = 0;
	cy6 = 0;
	cy7 = 0;
	cy8 = 0;
	cy9 = 0;

	*fracmax=0;	/* init max. fractional error	*/

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

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
			cmplx_carry_norm_errcheck0(a1p0r ,a1p0i ,cy0 ,bjmodn0    );
			cmplx_carry_norm_errcheck (a1p1r ,a1p1i ,cy1 ,bjmodn1 ,1 );
			cmplx_carry_norm_errcheck (a1p2r ,a1p2i ,cy2 ,bjmodn2 ,2 );
			cmplx_carry_norm_errcheck (a1p3r ,a1p3i ,cy3 ,bjmodn3 ,3 );
			cmplx_carry_norm_errcheck (a1p4r ,a1p4i ,cy4 ,bjmodn4 ,4 );
			cmplx_carry_norm_errcheck (a1p5r ,a1p5i ,cy5 ,bjmodn5 ,5 );
			cmplx_carry_norm_errcheck (a1p6r ,a1p6i ,cy6 ,bjmodn6 ,6 );
			cmplx_carry_norm_errcheck (a1p7r ,a1p7i ,cy7 ,bjmodn7 ,7 );
			cmplx_carry_norm_errcheck (a1p8r ,a1p8i ,cy8 ,bjmodn8 ,8 );
			cmplx_carry_norm_errcheck (a1p9r ,a1p9i ,cy9 ,bjmodn9 ,9 );

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

			iroot += root_incr;		/* increment sincos index.	*/
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

	iroot = 0;
	root_incr = 0;
	scale = 1;

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
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix10_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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
	*fracmax = maxerr;
	return(0);
}

/***************/

int radix10_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter                 , uint64 p)
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

/**************** MULTI-THREADED VERSION ********************/

#ifdef MULTITHREAD

	int n10,ithread,i,j,j1,j2,iroot,root_incr,k1,k2,k,khi,l,outer;
	int nt_sz_int;
	int nt_sz_dbl;
	static int *_i, *bjmodn0 = 0x0, *bjmodn1 = 0x0, *bjmodn2 = 0x0, *bjmodn3 = 0x0, *bjmodn4 = 0x0, *bjmodn5 = 0x0, *bjmodn6 = 0x0, *bjmodn7 = 0x0, *bjmodn8 = 0x0, *bjmodn9 = 0x0
	, *bjmodnini = 0x0, *jstart = 0x0, *jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static uint64 psave=0;
	static uint32 bw,sw,nm1,pini,p1,p2,p3,p4,p5,p6,p7,p8,p9;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double	cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					sn2 =  0.95105651629515357211,	/*  sin(u) */
					ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					ss2 =  0.36327126400268044292,	/* [sin(u)-sin(2u)] */
					radix_inv, n2inv
	, *cy0 = 0x0, *cy1 = 0x0, *cy2 = 0x0, *cy3 = 0x0, *cy4 = 0x0, *cy5 = 0x0, *cy6 = 0x0, *cy7 = 0x0, *cy8 = 0x0, *cy9 = 0x0;
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
	,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r
	,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i
	,temp,scale;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;
	static uint32 CY_THREADS;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;

/*...change n10 and n_div_wt to non-static to work around a gcc compiler bug. */
	n10   = n/10;
	n_div_nwt = n10 >> nwt_bits;
	if(NTHREADS > 1)
	{
		ASSERT(HERE, n10      %NTHREADS == 0,"radix10_ditN_cy_dif1.c: n10      %NTHREADS == 0");
		ASSERT(HERE, n_div_nwt%NTHREADS == 0,"radix10_ditN_cy_dif1.c: n_div_nwt%NTHREADS == 0");
	}

	if((n_div_nwt << nwt_bits) != n10)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/10 in radix10_ditN_cy_dif1.\n",iter);
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
		/* #threads used in carry step must be a power of 2, so use the largest
		power of 2 that is <= the value of the global NTHREADS:
		*/
		CY_THREADS = (((uint32)NTHREADS << leadz32(NTHREADS)) & 0x80000000) >> leadz32(NTHREADS);

		ASSERT(HERE,  CY_THREADS <= NTHREADS                ,"radix16_ditN_cy_dif1.c: ");
		ASSERT(HERE, (CY_THREADS >> trailz32(CY_THREADS)) == 1,"radix16_ditN_cy_dif1.c: ");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, CY_THREADS <= NTHREADS   ,"radix16_ditN_cy_dif1.c: ");
			ASSERT(HERE, n16      %CY_THREADS == 0,"radix16_ditN_cy_dif1.c: ");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"radix16_ditN_cy_dif1.c: ");
		}

	#if 1	/* CY_DBG */CY_THREADS
		fprintf(stderr,"Using %d threads in carry step\n", );
	#endif

		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)10));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		/* constant index offsets for array load/stores are here.	*/

		pini = n10/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
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

		if(cy0)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)bjmodn0); bjmodn0 = 0x0;
			free((void *)bjmodn1); bjmodn1 = 0x0;
			free((void *)bjmodn2); bjmodn2 = 0x0;
			free((void *)bjmodn3); bjmodn3 = 0x0;
			free((void *)bjmodn4); bjmodn4 = 0x0;
			free((void *)bjmodn5); bjmodn5 = 0x0;
			free((void *)bjmodn6); bjmodn6 = 0x0;
			free((void *)bjmodn7); bjmodn7 = 0x0;
			free((void *)bjmodn8); bjmodn8 = 0x0;
			free((void *)bjmodn9); bjmodn9 = 0x0;

			free((void *)cy0    ); cy0     = 0x0;
			free((void *)cy1    ); cy1     = 0x0;
			free((void *)cy2    ); cy2     = 0x0;
			free((void *)cy3    ); cy3     = 0x0;
			free((void *)cy4    ); cy4     = 0x0;
			free((void *)cy5    ); cy5     = 0x0;
			free((void *)cy6    ); cy6     = 0x0;
			free((void *)cy7    ); cy7     = 0x0;
			free((void *)cy8    ); cy8     = 0x0;
			free((void *)cy9    ); cy9     = 0x0;

			free((void *)jstart ); jstart  = 0x0;
			free((void *)jhi    ); jhi     = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)bjmodnini); bjmodnini = 0x0;
		}

		nt_sz_int = NTHREADS*sizeof(int);
		nt_sz_dbl = NTHREADS*sizeof(double);

		_i      = (int *)malloc(nt_sz_int); if(!_i     ){ sprintf(cbuf,"FATAL: unable to allocate array _i      in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn0 = (int *)malloc(nt_sz_int); if(!bjmodn0){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn0 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn1 = (int *)malloc(nt_sz_int); if(!bjmodn1){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn1 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn2 = (int *)malloc(nt_sz_int); if(!bjmodn2){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn2 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn3 = (int *)malloc(nt_sz_int); if(!bjmodn3){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn3 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn4 = (int *)malloc(nt_sz_int); if(!bjmodn4){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn4 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn5 = (int *)malloc(nt_sz_int); if(!bjmodn5){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn5 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn6 = (int *)malloc(nt_sz_int); if(!bjmodn6){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn6 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn7 = (int *)malloc(nt_sz_int); if(!bjmodn7){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn7 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn8 = (int *)malloc(nt_sz_int); if(!bjmodn8){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn8 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodn9 = (int *)malloc(nt_sz_int); if(!bjmodn9){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn9 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		jstart  = (int *)malloc(nt_sz_int); if(!jstart ){ sprintf(cbuf,"FATAL: unable to allocate array jstart  in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		jhi     = (int *)malloc(nt_sz_int); if(!jhi    ){ sprintf(cbuf,"FATAL: unable to allocate array jhi     in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_col    = (int *)malloc(nt_sz_int); if(!_col   ){ sprintf(cbuf,"FATAL: unable to allocate array _col    in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_co2    = (int *)malloc(nt_sz_int); if(!_co2   ){ sprintf(cbuf,"FATAL: unable to allocate array _co2    in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_co3    = (int *)malloc(nt_sz_int); if(!_co3   ){ sprintf(cbuf,"FATAL: unable to allocate array _co3    in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }

		cy0     = (double *)malloc(nt_sz_dbl); if(!bjmodn0){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn0 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy1     = (double *)malloc(nt_sz_dbl); if(!bjmodn1){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn1 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy2     = (double *)malloc(nt_sz_dbl); if(!bjmodn2){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn2 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy3     = (double *)malloc(nt_sz_dbl); if(!bjmodn3){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn3 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy4     = (double *)malloc(nt_sz_dbl); if(!bjmodn4){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn4 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy5     = (double *)malloc(nt_sz_dbl); if(!bjmodn5){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn5 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy6     = (double *)malloc(nt_sz_dbl); if(!bjmodn6){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn6 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy7     = (double *)malloc(nt_sz_dbl); if(!bjmodn7){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn7 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy8     = (double *)malloc(nt_sz_dbl); if(!bjmodn8){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn8 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		cy9     = (double *)malloc(nt_sz_dbl); if(!bjmodn9){ sprintf(cbuf,"FATAL: unable to allocate array bjmodn9 in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }

		bjmodnini = (int *)malloc((NTHREADS + 1)*sizeof(int)); if(!bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array bjmodnini in radix10_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		bjmodnini[0] = 0;
		bjmodnini[1] = 0;
		for(j=0; j < n10/NTHREADS; j++)
		{
			bjmodnini[1] -= sw; bjmodnini[1] = bjmodnini[1] + ( (-(int)((uint32)bjmodnini[1] >> 31)) & n);
		}

		if(NTHREADS > 1)
		{
			for(ithread = 2; ithread <= NTHREADS; ithread++)
			{
				bjmodnini[ithread] = bjmodnini[ithread-1] + bjmodnini[1] - n; bjmodnini[ithread] = bjmodnini[ithread] + ( (-(int)((uint32)bjmodnini[ithread] >> 31)) & n);
            }
		}
	}

/*...The radix-10 final DIT pass is here.	*/

    for(ithread = 0; ithread < NTHREADS; ithread++)
    {
		/* init carries	*/
		cy0[ithread] = 0;
		cy1[ithread] = 0;
		cy2[ithread] = 0;
		cy3[ithread] = 0;
		cy4[ithread] = 0;
		cy5[ithread] = 0;
		cy6[ithread] = 0;
		cy7[ithread] = 0;
		cy8[ithread] = 0;
		cy9[ithread] = 0;
    }
    cy0[0] = -2;

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

for(outer=0; outer <= 1; outer++)
{
	_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/
	if(NTHREADS > 1)
	{
		for(ithread = 1; ithread < NTHREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - bjmodnini[ithread]) >> 31);
		}
	}

	/*
	Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
	then simply overwrite it with 1 prior to starting the k-loop.
	*/
	khi = n_div_nwt/NTHREADS;

    for(ithread = 0; ithread < NTHREADS; ithread++)
    {
		bjmodn0[ithread] = bjmodnini[ithread];
		bjmodn1[ithread] = bjmodn0[ithread] + bjmodnini[NTHREADS] - n; bjmodn1[ithread] = bjmodn1[ithread] + ( (-(int)((uint32)bjmodn1[ithread] >> 31)) & n);
		bjmodn2[ithread] = bjmodn1[ithread] + bjmodnini[NTHREADS] - n; bjmodn2[ithread] = bjmodn2[ithread] + ( (-(int)((uint32)bjmodn2[ithread] >> 31)) & n);
		bjmodn3[ithread] = bjmodn2[ithread] + bjmodnini[NTHREADS] - n; bjmodn3[ithread] = bjmodn3[ithread] + ( (-(int)((uint32)bjmodn3[ithread] >> 31)) & n);
		bjmodn4[ithread] = bjmodn3[ithread] + bjmodnini[NTHREADS] - n; bjmodn4[ithread] = bjmodn4[ithread] + ( (-(int)((uint32)bjmodn4[ithread] >> 31)) & n);
		bjmodn5[ithread] = bjmodn4[ithread] + bjmodnini[NTHREADS] - n; bjmodn5[ithread] = bjmodn5[ithread] + ( (-(int)((uint32)bjmodn5[ithread] >> 31)) & n);
		bjmodn6[ithread] = bjmodn5[ithread] + bjmodnini[NTHREADS] - n; bjmodn6[ithread] = bjmodn6[ithread] + ( (-(int)((uint32)bjmodn6[ithread] >> 31)) & n);
		bjmodn7[ithread] = bjmodn6[ithread] + bjmodnini[NTHREADS] - n; bjmodn7[ithread] = bjmodn7[ithread] + ( (-(int)((uint32)bjmodn7[ithread] >> 31)) & n);
		bjmodn8[ithread] = bjmodn7[ithread] + bjmodnini[NTHREADS] - n; bjmodn8[ithread] = bjmodn8[ithread] + ( (-(int)((uint32)bjmodn8[ithread] >> 31)) & n);
		bjmodn9[ithread] = bjmodn8[ithread] + bjmodnini[NTHREADS] - n; bjmodn9[ithread] = bjmodn9[ithread] + ( (-(int)((uint32)bjmodn9[ithread] >> 31)) & n);

		jstart[ithread] = ithread*n10/NTHREADS;
		if(root_incr==0)
			jhi   [ithread] = jstart[ithread] + 3;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			jhi   [ithread] = jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*10);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+10 - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-10;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
    }
    /* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
    if(root_incr==0) khi = 1;

	/*printf("Using %d threads in carry step\n", NTHREADS);*/

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifndef USE_OMP
	#error Current code is OMP-only!
#endif
omp_set_num_threads(NTHREADS);
#pragma omp parallel for private(i,j,j1,k,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,rt,it,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,temp) default(shared) schedule(dynamic)

	for(ithread = 0; ithread < NTHREADS; ithread++)
	{
		i = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/

	/*printf("KHI = %d\n",khi);*/
		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
	/*printf("JLO = %d, JHI = %d\n",jstart[ithread],jhi[ithread]);*/
			for(j=jstart[ithread]; j < jhi[ithread]; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
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
		normally be getting dispatched to 10 separate blocks of the A-array, we need 10 separate carries.	*/

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
				col = _col[ithread];
				co2 = _co2[ithread];
				co3 = _co3[ithread];

				cmplx_carry_norm_nocheck0(a1p0r ,a1p0i ,cy0[ithread] ,bjmodn0[ithread]    );
				cmplx_carry_norm_nocheck (a1p1r ,a1p1i ,cy1[ithread] ,bjmodn1[ithread] ,1 );
				cmplx_carry_norm_nocheck (a1p2r ,a1p2i ,cy2[ithread] ,bjmodn2[ithread] ,2 );
				cmplx_carry_norm_nocheck (a1p3r ,a1p3i ,cy3[ithread] ,bjmodn3[ithread] ,3 );
				cmplx_carry_norm_nocheck (a1p4r ,a1p4i ,cy4[ithread] ,bjmodn4[ithread] ,4 );
				cmplx_carry_norm_nocheck (a1p5r ,a1p5i ,cy5[ithread] ,bjmodn5[ithread] ,5 );
				cmplx_carry_norm_nocheck (a1p6r ,a1p6i ,cy6[ithread] ,bjmodn6[ithread] ,6 );
				cmplx_carry_norm_nocheck (a1p7r ,a1p7i ,cy7[ithread] ,bjmodn7[ithread] ,7 );
				cmplx_carry_norm_nocheck (a1p8r ,a1p8i ,cy8[ithread] ,bjmodn8[ithread] ,8 );
				cmplx_carry_norm_nocheck (a1p9r ,a1p9i ,cy9[ithread] ,bjmodn9[ithread] ,9 );

				i =((uint32)(sw - bjmodn0[ithread]) >> 31);	/* get ready for the next set...	*/
				_co2[ithread] = _co3[ithread];	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

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
				t1 =a1p0r;			t2 =a1p0i;
				t3 =a1p8r+a1p2r;	t4 =a1p8i+a1p2i;
				t5 =a1p6r+a1p4r;	t6 =a1p6i+a1p4i;
				t7 =a1p6r-a1p4r;	t8 =a1p6i-a1p4i;
				t9 =a1p8r-a1p2r;	t10=a1p8i-a1p2i;
				rt = t3+t5;			it = t4+t6;
				t1 = t1+rt;			t2 = t2+it;		/* y0	*/

				rt = t1+cc1*rt;		it = t2+cc1*it;
				t5 = cc2*(t3-t5);	t6 = cc2*(t4-t6);
				t3 = rt+t5;			t4 = it+t6;
				t5 = rt-t5;			t6 = it-t6;
				rt = sn2*(t9-t7);	it = sn2*(t10-t8);
				t7 = ss1* t7;		t8 = ss1* t8;
				t9 = ss2* t9;		t10= ss2* t10;
	#if PFETCH
	addr = add0+p1;
	prefetch_p_doubles(addr);
	#endif
				t7 = rt+t7;			t8 = it+t8;
				t9 = rt-t9;			t10= it-t10;
				rt = t7;			it = t8;
				t7 = t3+it;			t8 = t4-rt;		/* y4	*/
				t3 = t3-it;			t4 = t4+rt;		/* y1	*/
				rt = t9;			it = t10;
				t9 = t5+it;			t10= t6-rt;		/* y3	*/
				t5 = t5-it;			t6 = t6+rt;		/* y2	*/
	#if PFETCH
	addr = add0+p2;
	prefetch_p_doubles(addr);
	#endif
				t11 =a1p5r;			t12 =a1p5i;
				t13 =a1p3r+a1p7r;	t14 =a1p3i+a1p7i;
				t15 =a1p1r+a1p9r;	t16 =a1p1i+a1p9i;
				t17 =a1p1r-a1p9r;	t18 =a1p1i-a1p9i;
				t19 =a1p3r-a1p7r;	t20 =a1p3i-a1p7i;
				rt  = t13+t15;		it  = t14+t16;
				t11 = t11+rt;		t12 = t12+it;		/* z0	*/
	#if PFETCH
	addr = add0+p3;
	prefetch_p_doubles(addr);
	#endif
				rt  = t11+cc1*rt;	it  = t12+cc1*it;
				t15 = cc2*(t13-t15);t16 = cc2*(t14-t16);
				t13 = rt+t15;		t14 = it+t16;
				t15 = rt-t15;		t16 = it-t16;
				rt  = sn2*(t19-t17);it  = sn2*(t20-t18);
				t17 = ss1* t17;		t18 = ss1* t18;
				t19 = ss2* t19;		t20 = ss2* t20;
	#if PFETCH
	addr = add0+p4;
	prefetch_p_doubles(addr);
	#endif
				t17 = rt+t17;		t18 = it+t18;
				t19 = rt-t19;		t20 = it-t20;
				rt  = t17;			it  = t18;
				t17 = t13+it;		t18 = t14-rt;		/* z4	*/
				t13 = t13-it;		t14 = t14+rt;		/* z1	*/
				rt  = t19;			it  = t20;
				t19 = t15+it;		t20 = t16-rt;		/* z3	*/
				t15 = t15-it;		t16 = t16+rt;		/* z2	*/

				/* ...and now do five radix-2 transforms.	*/
	#if PFETCH
	addr = add0+p5;
	prefetch_p_doubles(addr);
	#endif
				a[j1   ]=t1+t11;	a[j2   ]=t2+t12;
				a[j1+p1]=t1-t11;	a[j2+p1]=t2-t12;
	#if PFETCH
	addr = add0+p6;
	prefetch_p_doubles(addr);
	#endif
				a[j1+p8]=t3 +t13;	a[j2+p8]=t4 +t14;
				a[j1+p9]=t3 -t13;	a[j2+p9]=t4 -t14;
	#if PFETCH
	addr = add0+p7;
	prefetch_p_doubles(addr);
	#endif
				a[j1+p7]=t5 +t15;	a[j2+p7]=t6 +t16;
				a[j1+p6]=t5 -t15;	a[j2+p6]=t6 -t16;
	#if PFETCH
	addr = add0+p8;
	prefetch_p_doubles(addr);
	#endif
				a[j1+p4]=t9 +t19;	a[j2+p4]=t10+t20;
				a[j1+p5]=t9 -t19;	a[j2+p5]=t10-t20;
	#if PFETCH
	addr = add0+p9;
	prefetch_p_doubles(addr);
	#endif
				a[j1+p3]=t7 +t17;	a[j2+p3]=t8 +t18;
				a[j1+p2]=t7 -t17;	a[j2+p2]=t8 -t18;

				iroot += root_incr;		/* increment sincos index.	*/
			}	/* endfor(j) */

			jstart[ithread] += nwt;
			jhi   [ithread] += nwt;
			_col[ithread] += 10;
			_co3[ithread] -= 10;
		}	/* endfor(k) */
	}	/* endfor(ithread) */
	/******* END OF PARALLEL FOR-LOOP ********/

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-10 forward DIF FFT of the first block of 10 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 10 outputs of (1);
!   (3) Reweight and perform a radix-10 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 10 elements and repeat (1-4).
*/
	t1  = cy9[NTHREADS - 1];

    if(NTHREADS > 1)
    {
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy9[ithread] = cy9[ithread-1]; }; cy9[0] = cy8[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy8[ithread] = cy8[ithread-1]; }; cy8[0] = cy7[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy7[ithread] = cy7[ithread-1]; }; cy7[0] = cy6[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy6[ithread] = cy6[ithread-1]; }; cy6[0] = cy5[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy5[ithread] = cy5[ithread-1]; }; cy5[0] = cy4[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy4[ithread] = cy4[ithread-1]; }; cy4[0] = cy3[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy3[ithread] = cy3[ithread-1]; }; cy3[0] = cy2[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy2[ithread] = cy2[ithread-1]; }; cy2[0] = cy1[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy1[ithread] = cy1[ithread-1]; }; cy1[0] = cy0[NTHREADS - 1];
		for(ithread = NTHREADS - 1; ithread > 0; ithread--){ cy0[ithread] = cy0[ithread-1]; }; cy0[0] = t1;
    }
    else
    {
		cy9[0]= cy8[0];
		cy8[0]= cy7[0];
		cy7[0]= cy6[0];
		cy6[0]= cy5[0];
		cy5[0]= cy4[0];
		cy4[0]= cy3[0];
		cy3[0]= cy2[0];
		cy2[0]= cy1[0];
		cy1[0]= cy0[0];
		cy0[0]= t1;
    }

	iroot = 0;
	root_incr = 0;
	scale = 1;

    for(ithread = 0; ithread < NTHREADS; ithread++)
    {
		for(j = ithread*pini; j <= ithread*pini + 3; j++)
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
}	/* endfor(outer) */

    t1 = 0;
    for(ithread = 0; ithread < NTHREADS; ithread++)
    {
		t1 += fabs(cy0[ithread])+fabs(cy1[ithread])+fabs(cy2[ithread])+fabs(cy3[ithread])+fabs(cy4[ithread])+fabs(cy5[ithread])+fabs(cy6[ithread])+fabs(cy7[ithread])+fabs(cy8[ithread])+fabs(cy9[ithread]);
    }

	if(t1 != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix10_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

#else	/**************** SINGLE-THREADED VERSION ********************/

	int n10,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					sn2 =  0.95105651629515357211,	/*  sin(u) */
					ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					ss2 =  0.36327126400268044292,	/* [sin(u)-sin(2u)] */
					radix_inv, n2inv;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,temp,scale;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;

/*...change n10 and n_div_wt to non-static to work around a gcc compiler bug. */
	n10   = n/10;
	n_div_nwt = n10 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n10)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/10 in radix10_ditN_cy_dif1.\n",iter);
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

	cy0 =-2;		/* init carries	*/
	cy1 = 0;
	cy2 = 0;
	cy3 = 0;
	cy4 = 0;
	cy5 = 0;
	cy6 = 0;
	cy7 = 0;
	cy8 = 0;
	cy9 = 0;

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

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
			cmplx_carry_norm_nocheck0(a1p0r ,a1p0i ,cy0 ,bjmodn0    );
			cmplx_carry_norm_nocheck (a1p1r ,a1p1i ,cy1 ,bjmodn1 ,1 );
			cmplx_carry_norm_nocheck (a1p2r ,a1p2i ,cy2 ,bjmodn2 ,2 );
			cmplx_carry_norm_nocheck (a1p3r ,a1p3i ,cy3 ,bjmodn3 ,3 );
			cmplx_carry_norm_nocheck (a1p4r ,a1p4i ,cy4 ,bjmodn4 ,4 );
			cmplx_carry_norm_nocheck (a1p5r ,a1p5i ,cy5 ,bjmodn5 ,5 );
			cmplx_carry_norm_nocheck (a1p6r ,a1p6i ,cy6 ,bjmodn6 ,6 );
			cmplx_carry_norm_nocheck (a1p7r ,a1p7i ,cy7 ,bjmodn7 ,7 );
			cmplx_carry_norm_nocheck (a1p8r ,a1p8i ,cy8 ,bjmodn8 ,8 );
			cmplx_carry_norm_nocheck (a1p9r ,a1p9i ,cy9 ,bjmodn9 ,9 );

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
			t1 =a1p0r;			t2 =a1p0i;
			t3 =a1p8r+a1p2r;	t4 =a1p8i+a1p2i;
			t5 =a1p6r+a1p4r;	t6 =a1p6i+a1p4i;
			t7 =a1p6r-a1p4r;	t8 =a1p6i-a1p4i;
			t9 =a1p8r-a1p2r;	t10=a1p8i-a1p2i;
			rt = t3+t5;			it = t4+t6;
			t1 = t1+rt;			t2 = t2+it;		/* y0	*/

			rt = t1+cc1*rt;		it = t2+cc1*it;
			t5 = cc2*(t3-t5);	t6 = cc2*(t4-t6);
			t3 = rt+t5;			t4 = it+t6;
			t5 = rt-t5;			t6 = it-t6;
			rt = sn2*(t9-t7);	it = sn2*(t10-t8);
			t7 = ss1* t7;		t8 = ss1* t8;
			t9 = ss2* t9;		t10= ss2* t10;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
			t7 = rt+t7;			t8 = it+t8;
			t9 = rt-t9;			t10= it-t10;
			rt = t7;			it = t8;
			t7 = t3+it;			t8 = t4-rt;		/* y4	*/
			t3 = t3-it;			t4 = t4+rt;		/* y1	*/
			rt = t9;			it = t10;
			t9 = t5+it;			t10= t6-rt;		/* y3	*/
			t5 = t5-it;			t6 = t6+rt;		/* y2	*/
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif
			t11 =a1p5r;			t12 =a1p5i;
			t13 =a1p3r+a1p7r;	t14 =a1p3i+a1p7i;
			t15 =a1p1r+a1p9r;	t16 =a1p1i+a1p9i;
			t17 =a1p1r-a1p9r;	t18 =a1p1i-a1p9i;
			t19 =a1p3r-a1p7r;	t20 =a1p3i-a1p7i;
			rt  = t13+t15;		it  = t14+t16;
			t11 = t11+rt;		t12 = t12+it;		/* z0	*/
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
			rt  = t11+cc1*rt;	it  = t12+cc1*it;
			t15 = cc2*(t13-t15);t16 = cc2*(t14-t16);
			t13 = rt+t15;		t14 = it+t16;
			t15 = rt-t15;		t16 = it-t16;
			rt  = sn2*(t19-t17);it  = sn2*(t20-t18);
			t17 = ss1* t17;		t18 = ss1* t18;
			t19 = ss2* t19;		t20 = ss2* t20;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
			t17 = rt+t17;		t18 = it+t18;
			t19 = rt-t19;		t20 = it-t20;
			rt  = t17;			it  = t18;
			t17 = t13+it;		t18 = t14-rt;		/* z4	*/
			t13 = t13-it;		t14 = t14+rt;		/* z1	*/
			rt  = t19;			it  = t20;
			t19 = t15+it;		t20 = t16-rt;		/* z3	*/
			t15 = t15-it;		t16 = t16+rt;		/* z2	*/

		/* ...and now do five radix-2 transforms.	*/
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
			a[j1   ]=t1+t11;	a[j2   ]=t2+t12;
			a[j1+p1]=t1-t11;	a[j2+p1]=t2-t12;
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif
			a[j1+p8]=t3 +t13;	a[j2+p8]=t4 +t14;
			a[j1+p9]=t3 -t13;	a[j2+p9]=t4 -t14;
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
			a[j1+p7]=t5 +t15;	a[j2+p7]=t6 +t16;
			a[j1+p6]=t5 -t15;	a[j2+p6]=t6 -t16;
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif
			a[j1+p4]=t9 +t19;	a[j2+p4]=t10+t20;
			a[j1+p5]=t9 -t19;	a[j2+p5]=t10-t20;
#if PFETCH
addr = add0+p9;
prefetch_p_doubles(addr);
#endif
			a[j1+p3]=t7 +t17;	a[j2+p3]=t8 +t18;
			a[j1+p2]=t7 -t17;	a[j2+p2]=t8 -t18;

			iroot += root_incr;		/* increment sincos index.	*/
		}	/* endfor(j) */

		jstart += nwt;
		jhi    += nwt;
		col += 10;
		co3 -= 10;

	}	/* endfor(k) */

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-10 forward DIF FFT of the first block of 10 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 10 outputs of (1);
!   (3) Reweight and perform a radix-10 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 10 elements and repeat (1-4).
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

	iroot = 0;
	root_incr = 0;
	scale = 1;

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
}	/* endfor(outer) */

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix10_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

#endif	/**************** END MULTITHREAD #IFDEF ********************/
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

