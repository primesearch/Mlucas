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

/***************/

int radix12_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-12 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-14 complex DIF pass on the data in
!   the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n12,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10
		,bjmodn11,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double c = .86602540378443864676, c4m1 = -1.5, radix_inv, n2inv;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,cy10,cy11,temp,frac,scale;
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

/*...change n12 and n_div_wt to non-static to work around a gcc compiler bug. */
	n12   = n/12;
	n_div_nwt = n12 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n12)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/12 in radix12_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)12));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p1 = n12;
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

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );

		bjmodnini=0;
		for(j=0; j < n12; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-12 final DIT pass is here.	*/

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
	cy10= 0;
	cy11= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy0 = -2;
	}
	else
	{
		ASSERT(HERE,0,"Radix-12 currently only supports LL test mode!");
	}

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
	bjmodn10= bjmodn9 +bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+12;
	co3=co2-12;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

	/*...Block 1:	*/
			t1 =a[j1    ];			t2 =a[j2    ];
			rt =a[j1+p1 ];			it =a[j2+p1 ];
			t3 =t1 -rt;				t1 =t1 +rt;
			t4 =t2 -it;				t2 =t2 +it;

			t5 =a[j1+p3 ];			t6 =a[j2+p3 ];
			rt =a[j1+p2 ];			it =a[j2+p2 ];
			t7 =t5 -rt;  			t5 =t5 +rt;
			t8 =t6 -it;  			t6 =t6 +it;

			rt =t5;	t5 =t1 -rt;		t1 =t1 +rt;
			it =t6;	t6 =t2 -it;		t2 =t2 +it;

			rt =t7;	t7 =t3 -t8;		t3 =t3 +t8;
					t8 =t4 +rt;		t4 =t4 -rt;
	/*...Block 2:	*/
			t9 =a[j1+p6 ];			t10=a[j2+p6 ];
			rt =a[j1+p7 ];			it =a[j2+p7 ];
			t11=t9 -rt;				t9 =t9 +rt;
			t12=t10-it;				t10=t10+it;

			t13=a[j1+p4 ];			t14=a[j2+p4 ];
			rt =a[j1+p5 ];			it =a[j2+p5 ];
			t15=t13-rt;  			t13=t13+rt;
			t16=t14-it;				t14=t14+it;

			rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
			it =t14;	t14=t10-it;	t10=t10+it;

			rt =t15;	t15=t11-t16;t11=t11+t16;
						t16=t12+rt;	t12=t12-rt;
	/*...Block 3:	*/
			t17=a[j1+p9 ];			t18=a[j2+p9 ];
			rt =a[j1+p8 ];			it =a[j2+p8 ];
			t19=t17-rt;  			t17=t17+rt;
			t20=t18-it;				t18=t18+it;

			t21=a[j1+p10];			t22=a[j2+p10];
			rt =a[j1+p11];			it =a[j2+p11];
			t23=t21-rt;  			t21=t21+rt;
			t24=t22-it;				t22=t22+it;

			rt =t21;	t21=t17-rt;	t17=t17+rt;
			it =t22;	t22=t18-it;	t18=t18+it;

			rt =t23;	t23=t19-t24;t19=t19+t24;
			t24=t20+rt;	t20=t20-rt;
	/*
	!...and now do four radix-3 transforms:
	*/
	/*...Block 1: t1,9,17	*/
			rt =t9;					it =t10;
			t9 =rt+t17;				t10=it+t18;
			t17=rt-t17;				t18=it-t18;
			t1 =t1+t9;				t2 =t2+t10;
			aj1p0r =t1;				aj1p0i =t2;
			t9 =t1+c4m1*t9;			t10=t2+c4m1*t10;
			rt =c*t17;				it =c*t18;
			aj1p4r =t9+it;			aj1p4i =t10-rt;
			aj1p8r =t9-it;			aj1p8i =t10+rt;
	/*...Block 2: t3,11,19	*/
			rt =t11;				it =t12;
			t11=rt+t19;				t12=it+t20;
			t19=rt-t19;				t20=it-t20;
			t3 =t3+t11;				t4 =t4+t12;
			aj1p3r =t3;				aj1p3i =t4;
			t11=t3+c4m1*t11;		t12=t4+c4m1*t12;
			rt =c*t19;				it =c*t20;
			aj1p7r =t11+it;			aj1p7i =t12-rt;
			aj1p11r=t11-it;			aj1p11i=t12+rt;
	/*...Block 3: t5,13,21	*/
			rt =t13;				it =t14;
			t13=rt+t21;				t14=it+t22;
			t21=rt-t21;				t22=it-t22;
			t5 =t5+t13;				t6 =t6+t14;
			aj1p6r =t5;				aj1p6i =t6;
			t13=t5+c4m1*t13;		t14=t6+c4m1*t14;
			rt =c*t21;				it =c*t22;
			aj1p10r=t13+it;			aj1p10i=t14-rt;
			aj1p2r =t13-it;			aj1p2i =t14+rt;
	/*...Block 3: t7,15,23	*/
			rt =t15;				it =t16;
			t15=rt+t23;				t16=it+t24;
			t23=rt-t23;				t24=it-t24;
			t7 =t7+t15;				t8 =t8+t16;
			aj1p9r =t7;				aj1p9i =t8;
			t15=t7+c4m1*t15;		t16=t8+c4m1*t16;
			rt =c*t23;				it =c*t24;
			aj1p1r =t15+it;			aj1p1i =t16-rt;
			aj1p5r =t15-it;			aj1p5i =t16+rt;

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 12 separate blocks of the A-array, we need 12 separate carries.	*/

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
		 cmplx_carry_norm_errcheck0(aj1p0r ,aj1p0i ,cy0 ,bjmodn0 );
			cmplx_carry_norm_errcheck(aj1p1r ,aj1p1i ,cy1 ,bjmodn1 ,1 );
			cmplx_carry_norm_errcheck(aj1p2r ,aj1p2i ,cy2 ,bjmodn2 ,2 );
			cmplx_carry_norm_errcheck(aj1p3r ,aj1p3i ,cy3 ,bjmodn3 ,3 );
			cmplx_carry_norm_errcheck(aj1p4r ,aj1p4i ,cy4 ,bjmodn4 ,4 );
			cmplx_carry_norm_errcheck(aj1p5r ,aj1p5i ,cy5 ,bjmodn5 ,5 );
			cmplx_carry_norm_errcheck(aj1p6r ,aj1p6i ,cy6 ,bjmodn6 ,6 );
			cmplx_carry_norm_errcheck(aj1p7r ,aj1p7i ,cy7 ,bjmodn7 ,7 );
			cmplx_carry_norm_errcheck(aj1p8r ,aj1p8i ,cy8 ,bjmodn8 ,8 );
			cmplx_carry_norm_errcheck(aj1p9r ,aj1p9i ,cy9 ,bjmodn9 ,9 );
			cmplx_carry_norm_errcheck(aj1p10r,aj1p10i,cy10,bjmodn10,10);
			cmplx_carry_norm_errcheck(aj1p11r,aj1p11i,cy11,bjmodn11,11);

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-12 DIF pass is here:	*/

/*...gather the needed data (12 64-bit complex, i.e. 24 64-bit reals) and do four radix-3 transforms...	*/
/* Twiddleless version requires us to swap inputs x4 <-> x8, x2 <-> x6, x7 <-> x11, x1 <-> x9... */
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif
			t1 =aj1p0r;				t2 =aj1p0i;
			t3 =aj1p8r +aj1p4r;		t4 =aj1p8i +aj1p4i;
			t5 =aj1p8r -aj1p4r;		t6 =aj1p8i -aj1p4i;
			t1 =t1+t3;				t2 =t2+t4;
			t3 =t1+c4m1*t3;			t4 =t2+c4m1*t4;
			rt =c*t5;				it =c*t6;
			t5 =t3+it;				t6 =t4-rt;
			t3 =t3-it;				t4 =t4+rt;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif

			t7 =aj1p9r;				t8 =aj1p9i;
			t9 =aj1p5r +aj1p1r;		t10=aj1p5i +aj1p1i;
			t11=aj1p5r -aj1p1r;		t12=aj1p5i -aj1p1i;
			t7 =t7+t9;				t8 =t8+t10;
			t9 =t7+c4m1*t9;			t10=t8+c4m1*t10;
			rt =c*t11;				it =c*t12;
			t11=t9+it;				t12=t10-rt;
			t9 =t9-it;				t10=t10+rt;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif

			t13=aj1p6r;				t14=aj1p6i;
			t15=aj1p2r +aj1p10r;	t16=aj1p2i +aj1p10i;
			t17=aj1p2r -aj1p10r;	t18=aj1p2i -aj1p10i;
			t13=t13+t15;			t14=t14+t16;
			t15=t13+c4m1*t15;		t16=t14+c4m1*t16;
			rt =c*t17;				it =c*t18;
			t17=t15+it;				t18=t16-rt;
			t15=t15-it;				t16=t16+rt;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif

			t19=aj1p3r;				t20=aj1p3i;
			t21=aj1p11r+aj1p7r;		t22=aj1p11i+aj1p7i;
			t23=aj1p11r-aj1p7r;		t24=aj1p11i-aj1p7i;
			t19=t19+t21;			t20=t20+t22;
			t21=t19+c4m1*t21;		t22=t20+c4m1*t22;
			rt =c*t23;				it =c*t24;
			t23=t21+it;				t24=t22-rt;
			t21=t21-it;				t22=t22+rt;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
/*
!...and now do three radix-4 transforms:
*/
/*...Block 1: t1,7,13,19	*/
			rt =t13;	t13=t1 -rt;	t1 =t1 +rt;	/* t 1, 2 = b0+b6	*/
			it =t14 ;	t14=t2 -it;	t2 =t2 +it;	/* t13,14 = b0-b6	*/

			rt =t19;	t19=t7 -rt;	t7 =t7 +rt;	/* t 7, 8 = b3+b9	*/
			it =t20;	t20=t8 -it;	t8 =t8 +it;	/* t19,20 = b3-b9	*/
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
			a[j1    ]=t1+t7;		a[j2    ]=t2+t8;	/* b0+b3+b6+b9	*/
			a[j1+p1 ]=t1-t7;		a[j2+p1 ]=t2-t8;	/* b0-b3+b6-b9	*/

			a[j1+p2 ]=t13-t20;		a[j2+p2 ]=t14+t19;	/* mpy by I is inlined here...	*/
			a[j1+p3 ]=t13+t20;		a[j2+p3 ]=t14-t19;
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif

/*...Block 2: t3,9,15,21	*/
			rt =t15;	t15=t3 -rt;	t3 =t3 +rt;
			it =t16;	t16=t4 -it;	t4 =t4 +it;

			rt =t21;	t21=t9 -rt;	t9 =t9 +rt;
			it =t22;	t22=t10-it;	t10=t10+it;
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
			a[j1+p9 ]=t3+t9;		a[j2+p9 ]=t4+t10;
			a[j1+p8 ]=t3-t9;		a[j2+p8 ]=t4-t10;

			a[j1+p11]=t15-t22;		a[j2+p11]=t16+t21;	/* mpy by I is inlined here...	*/
			a[j1+p10]=t15+t22;		a[j2+p10]=t16-t21;
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif

/*...Block 3: t5,11,17,23	*/
			rt =t17;	t17=t5 -rt;	t5 =t5 +rt;
			it =t18;	t18=t6 -it;	t6 =t6 +it;

			rt =t23;	t23=t11-rt;	t11=t11+rt;
			it =t24;	t24=t12-it;	t12=t12+it;
#if PFETCH
addr = add0+p9;
prefetch_p_doubles(addr);
#endif
			a[j1+p6 ]=t5+t11;		a[j2+p6 ]=t6+t12;
			a[j1+p7 ]=t5-t11;		a[j2+p7 ]=t6-t12;

			a[j1+p5 ]=t17-t24;		a[j2+p5 ]=t18+t23;	/* mpy by I is inlined here...	*/
			a[j1+p4 ]=t17+t24;		a[j2+p4 ]=t18-t23;
#if PFETCH
addr = add0+p10;
prefetch_p_doubles(addr);
#endif

#if PFETCH
addr = add0+p11;
prefetch_p_doubles(addr);
#endif
		iroot += root_incr;		/* increment sincos index.	*/

		}

		jstart += nwt;
		jhi    += nwt;
		col += 12;
		co3 -= 12;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-12 forward DIF FFT of the first block of 12 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 12 outputs of (1);
!   (3) Reweight and perform a radix-12 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 12 elements and repeat (1-4).
*/

	t1  = cy11;
	cy11= cy10;
	cy10= cy9;
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
		a[j+p10] *= radix_inv;
		a[j+p11] *= radix_inv;
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9)+fabs(cy10)+fabs(cy11) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix12_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

int radix12_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter                 , uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-12 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-14 complex DIF pass on the data in
!   the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n12,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10
		,bjmodn11,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double c = .86602540378443864676, c4m1 = -1.5, radix_inv, n2inv;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,cy10,cy11,temp,scale;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;

/*...change n12 and n_div_wt to non-static to work around a gcc compiler bug. */
	n12   = n/12;
	n_div_nwt = n12 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n12)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/12 in radix12_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)12));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p1 = n12;
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

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );

		bjmodnini=0;
		for(j=0; j < n12; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-12 final DIT pass is here.	*/

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
	cy10= 0;
	cy11= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy0 = -2;
	}
	else
	{
		ASSERT(HERE,0,"Radix-12 currently only supports LL test mode!");
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
	bjmodn10= bjmodn9 +bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+12;
	co3=co2-12;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

/*...Block 1:	*/
			t1 =a[j1    ];	t2 =a[j2    ];
			rt =a[j1+p1 ];	it =a[j2+p1 ];
			t3 =t1 -rt;		t1 =t1 +rt;
			t4 =t2 -it;		t2 =t2 +it;

			t5 =a[j1+p3 ];	t6 =a[j2+p3 ];
			rt =a[j1+p2 ];	it =a[j2+p2 ];
			t7 =t5 -rt;  	t5 =t5 +rt;
			t8 =t6 -it;  	t6 =t6 +it;

			rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
			it =t6;	t6 =t2 -it;	t2 =t2 +it;

			rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
			t8 =t4 +rt;	t4 =t4 -rt;
/*...Block 2:	*/
			t9 =a[j1+p6 ];	t10=a[j2+p6 ];
			rt =a[j1+p7 ];	it =a[j2+p7 ];
			t11=t9 -rt;		t9 =t9 +rt;
			t12=t10-it;		t10=t10+it;

			t13=a[j1+p4 ];	t14=a[j2+p4 ];
			rt =a[j1+p5 ];	it =a[j2+p5 ];
			t15=t13-rt;  	t13=t13+rt;
			t16=t14-it;		t14=t14+it;

			rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
			it =t14;	t14=t10-it;	t10=t10+it;

			rt =t15;	t15=t11-t16;	t11=t11+t16;
			t16=t12+rt;	t12=t12-rt;
/*...Block 3:	*/
			t17=a[j1+p9 ];	t18=a[j2+p9 ];
			rt =a[j1+p8 ];	it =a[j2+p8 ];
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
/*
!...and now do four radix-3 transforms:
*/
/*...Block 1: t1,9,17	*/
		rt =t9;			it =t10;
		t9 =rt+t17;			t10=it+t18;
		t17=rt-t17;			t18=it-t18;
		t1 =t1+t9;			t2 =t2+t10;
		aj1p0r =t1;			aj1p0i =t2;
		t9 =t1+c4m1*t9;		t10=t2+c4m1*t10;
		rt =c*t17;			it =c*t18;
		aj1p4r =t9+it;		aj1p4i =t10-rt;
		aj1p8r =t9-it;		aj1p8i =t10+rt;
/*...Block 2: t3,11,19	*/
		rt =t11;			it =t12;
		t11=rt+t19;			t12=it+t20;
		t19=rt-t19;			t20=it-t20;
		t3 =t3+t11;			t4 =t4+t12;
		aj1p3r =t3;			aj1p3i =t4;
		t11=t3+c4m1*t11;		t12=t4+c4m1*t12;
		rt =c*t19;			it =c*t20;
		aj1p7r =t11+it;		aj1p7i =t12-rt;
		aj1p11r=t11-it;		aj1p11i=t12+rt;
/*...Block 3: t5,13,21	*/
		rt =t13;			it =t14;
		t13=rt+t21;			t14=it+t22;
		t21=rt-t21;			t22=it-t22;
		t5 =t5+t13;			t6 =t6+t14;
		aj1p6r =t5;			aj1p6i =t6;
		t13=t5+c4m1*t13;		t14=t6+c4m1*t14;
		rt =c*t21;			it =c*t22;
		aj1p10r=t13+it;		aj1p10i=t14-rt;
		aj1p2r =t13-it;		aj1p2i =t14+rt;
/*...Block 3: t7,15,23	*/
		rt =t15;			it =t16;
		t15=rt+t23;			t16=it+t24;
		t23=rt-t23;			t24=it-t24;
		t7 =t7+t15;			t8 =t8+t16;
		aj1p9r =t7;			aj1p9i =t8;
		t15=t7+c4m1*t15;		t16=t8+c4m1*t16;
		rt =c*t23;			it =c*t24;
		aj1p1r =t15+it;		aj1p1i =t16-rt;
		aj1p5r =t15-it;		aj1p5i =t16+rt;

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 12 separate blocks of the A-array, we need 12 separate carries.	*/

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
		 cmplx_carry_norm_nocheck0(aj1p0r ,aj1p0i ,cy0 ,bjmodn0 );
			cmplx_carry_norm_nocheck(aj1p1r ,aj1p1i ,cy1 ,bjmodn1 ,1 );
			cmplx_carry_norm_nocheck(aj1p2r ,aj1p2i ,cy2 ,bjmodn2 ,2 );
			cmplx_carry_norm_nocheck(aj1p3r ,aj1p3i ,cy3 ,bjmodn3 ,3 );
			cmplx_carry_norm_nocheck(aj1p4r ,aj1p4i ,cy4 ,bjmodn4 ,4 );
			cmplx_carry_norm_nocheck(aj1p5r ,aj1p5i ,cy5 ,bjmodn5 ,5 );
			cmplx_carry_norm_nocheck(aj1p6r ,aj1p6i ,cy6 ,bjmodn6 ,6 );
			cmplx_carry_norm_nocheck(aj1p7r ,aj1p7i ,cy7 ,bjmodn7 ,7 );
			cmplx_carry_norm_nocheck(aj1p8r ,aj1p8i ,cy8 ,bjmodn8 ,8 );
			cmplx_carry_norm_nocheck(aj1p9r ,aj1p9i ,cy9 ,bjmodn9 ,9 );
			cmplx_carry_norm_nocheck(aj1p10r,aj1p10i,cy10,bjmodn10,10);
			cmplx_carry_norm_nocheck(aj1p11r,aj1p11i,cy11,bjmodn11,11);

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-12 DIF pass is here:	*/

/*...gather the needed data (12 64-bit complex, i.e. 24 64-bit reals) and do four radix-3 transforms...	*/
/* Twiddleless version requires us to swap inputs x4 <-> x8, x2 <-> x6, x7 <-> x11, x1 <-> x9... */
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif
		t1 =aj1p0r;			t2 =aj1p0i;
		t3 =aj1p8r +aj1p4r;		t4 =aj1p8i +aj1p4i;
		t5 =aj1p8r -aj1p4r;		t6 =aj1p8i -aj1p4i;
		t1 =t1+t3;			t2 =t2+t4;
		t3 =t1+c4m1*t3;		t4 =t2+c4m1*t4;
		rt =c*t5;			it =c*t6;
		t5 =t3+it;			t6 =t4-rt;
		t3 =t3-it;			t4 =t4+rt;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif

		t7 =aj1p9r;			t8 =aj1p9i;
		t9 =aj1p5r +aj1p1r;		t10=aj1p5i +aj1p1i;
		t11=aj1p5r -aj1p1r;		t12=aj1p5i -aj1p1i;
		t7 =t7+t9;			t8 =t8+t10;
		t9 =t7+c4m1*t9;		t10=t8+c4m1*t10;
		rt =c*t11;			it =c*t12;
		t11=t9+it;			t12=t10-rt;
		t9 =t9-it;			t10=t10+rt;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif

		t13=aj1p6r;			t14=aj1p6i;
		t15=aj1p2r +aj1p10r;		t16=aj1p2i +aj1p10i;
		t17=aj1p2r -aj1p10r;		t18=aj1p2i -aj1p10i;
		t13=t13+t15;			t14=t14+t16;
		t15=t13+c4m1*t15;		t16=t14+c4m1*t16;
		rt =c*t17;			it =c*t18;
		t17=t15+it;			t18=t16-rt;
		t15=t15-it;			t16=t16+rt;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif

		t19=aj1p3r;			t20=aj1p3i;
		t21=aj1p11r+aj1p7r;		t22=aj1p11i+aj1p7i;
		t23=aj1p11r-aj1p7r;		t24=aj1p11i-aj1p7i;
		t19=t19+t21;			t20=t20+t22;
		t21=t19+c4m1*t21;		t22=t20+c4m1*t22;
		rt =c*t23;			it =c*t24;
		t23=t21+it;			t24=t22-rt;
		t21=t21-it;			t22=t22+rt;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
/*
!...and now do three radix-4 transforms:
*/
/*...Block 1: t1,7,13,19	*/
			rt =t13;	t13=t1 -rt;	t1 =t1 +rt;	/* t 1, 2 = b0+b6	*/
			it =t14 ;	t14=t2 -it;	t2 =t2 +it;	/* t13,14 = b0-b6	*/

			rt =t19;	t19=t7 -rt;	t7 =t7 +rt;	/* t 7, 8 = b3+b9	*/
			it =t20;	t20=t8 -it;	t8 =t8 +it;	/* t19,20 = b3-b9	*/
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
			a[j1    ]=t1+t7;	a[j2    ]=t2+t8;	/* b0+b3+b6+b9	*/
			a[j1+p1 ]=t1-t7;	a[j2+p1 ]=t2-t8;	/* b0-b3+b6-b9	*/

			a[j1+p2 ]=t13-t20;	a[j2+p2 ]=t14+t19;	/* mpy by I is inlined here...	*/
			a[j1+p3 ]=t13+t20;	a[j2+p3 ]=t14-t19;
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif

/*...Block 2: t3,9,15,21	*/
			rt =t15;	t15=t3 -rt;	t3 =t3 +rt;
			it =t16;	t16=t4 -it;	t4 =t4 +it;

			rt =t21;	t21=t9 -rt;	t9 =t9 +rt;
			it =t22;	t22=t10-it;	t10=t10+it;
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
			a[j1+p9 ]=t3+t9;	a[j2+p9 ]=t4+t10;
			a[j1+p8 ]=t3-t9;	a[j2+p8 ]=t4-t10;

			a[j1+p11]=t15-t22;	a[j2+p11]=t16+t21;	/* mpy by I is inlined here...	*/
			a[j1+p10]=t15+t22;	a[j2+p10]=t16-t21;
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif

/*...Block 3: t5,11,17,23	*/
			rt =t17;	t17=t5 -rt;	t5 =t5 +rt;
			it =t18;	t18=t6 -it;	t6 =t6 +it;

			rt =t23;	t23=t11-rt;	t11=t11+rt;
			it =t24;	t24=t12-it;	t12=t12+it;
#if PFETCH
addr = add0+p9;
prefetch_p_doubles(addr);
#endif
			a[j1+p6 ]=t5+t11;	a[j2+p6 ]=t6+t12;
			a[j1+p7 ]=t5-t11;	a[j2+p7 ]=t6-t12;

			a[j1+p5 ]=t17-t24;	a[j2+p5 ]=t18+t23;	/* mpy by I is inlined here...	*/
			a[j1+p4 ]=t17+t24;	a[j2+p4 ]=t18-t23;
#if PFETCH
addr = add0+p10;
prefetch_p_doubles(addr);
#endif

#if PFETCH
addr = add0+p11;
prefetch_p_doubles(addr);
#endif
		iroot += root_incr;		/* increment sincos index.	*/

		}

		jstart += nwt;
		jhi    += nwt;
		col += 12;
		co3 -= 12;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-12 forward DIF FFT of the first block of 12 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 12 outputs of (1);
!   (3) Reweight and perform a radix-12 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 12 elements and repeat (1-4).
*/

	t1  = cy11;
	cy11= cy10;
	cy10= cy9;
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
		a[j+p10] *= radix_inv;
		a[j+p11] *= radix_inv;
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9)+fabs(cy10)+fabs(cy11) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix12_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

void radix12_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-12 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Uses an optimized [radix-3] x [radix-4] transform.
*/
	int j,j1,j2;
	static int n12,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, first_entry=TRUE;
	static double c = .86602540378443864676, c4m1 = -1.5;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24;

	if(!first_entry && (n/12) != n12)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n12=n/12;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n12;
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

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-12 pass is here.	*/

      for(j=0; j < n12; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*...gather the needed data (12 64-bit complex, i.e. 24 64-bit reals) and do four radix-3 transforms...	*/

/* Twiddleless version requires us to swap inputs x4 <-> x8, x2 <-> x6, x7 <-> x11 and x1 <-> x9...
indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11
      -> 0,-3,-6,-9, 8, 5, 2,-1, 4, 1,-2,-5
      == 0, 9, 6, 3, 8, 5, 2,11, 4, 1,10, 7 modulo 12.
I.e. start out with first triplet of indices {0,4,8}, permute those according to
{0,4,8}*11%12 = {0,8,4}, then each is head of a length-4 list of indices with decrement 3.
*/
#if 0
		t1 =a[j1    ];				t2 =a[j2    ];
		t3 =a[j1+p8 ]+a[j1+p4 ];	t4 =a[j2+p8 ]+a[j2+p4 ];
		t5 =a[j1+p8 ]-a[j1+p4 ];	t6 =a[j2+p8 ]-a[j2+p4 ];
		t1 =t1+t3;					t2 =t2+t4;
		t3 =t1+c4m1*t3;				t4 =t2+c4m1*t4;
		rt =c*t5;					it =c*t6;
		t5 =t3+it;					t6 =t4-rt;
		t3 =t3-it;					t4 =t4+rt;

		t7 =a[j1+p9 ];				t8 =a[j2+p9 ];
		t9 =a[j1+p5 ]+a[j1+p1 ];	t10=a[j2+p5 ]+a[j2+p1 ];
		t11=a[j1+p5 ]-a[j1+p1 ];	t12=a[j2+p5 ]-a[j2+p1 ];
		t7 =t7+t9;					t8 =t8+t10;
		t9 =t7+c4m1*t9;				t10=t8+c4m1*t10;
		rt =c*t11;					it =c*t12;
		t11=t9+it;					t12=t10-rt;
		t9 =t9-it;					t10=t10+rt;

		t13=a[j1+p6 ];				t14=a[j2+p6 ];
		t15=a[j1+p2 ]+a[j1+p10];	t16=a[j2+p2 ]+a[j2+p10];
		t17=a[j1+p2 ]-a[j1+p10];	t18=a[j2+p2 ]-a[j2+p10];
		t13=t13+t15;				t14=t14+t16;
		t15=t13+c4m1*t15;			t16=t14+c4m1*t16;
		rt =c*t17;					it =c*t18;
		t17=t15+it;					t18=t16-rt;
		t15=t15-it;					t16=t16+rt;

		t19=a[j1+p3 ];				t20=a[j2+p3 ];
		t21=a[j1+p11]+a[j1+p7 ];	t22=a[j2+p11]+a[j2+p7 ];
		t23=a[j1+p11]-a[j1+p7 ];	t24=a[j2+p11]-a[j2+p7 ];
		t19=t19+t21;				t20=t20+t22;
		t21=t19+c4m1*t21;			t22=t20+c4m1*t22;
		rt =c*t23;					it =c*t24;
		t23=t21+it;					t24=t22-rt;
		t21=t21-it;					t22=t22+rt;
/*
!...and now do three radix-4 transforms:
*/
/*...Block 1: t1,7,13,19	*/
		rt =t13;	t13=t1 -rt;		t1 =t1 +rt;
		it =t14 ;	t14=t2 -it;		t2 =t2 +it;

		rt =t19;	t19=t7 -rt;		t7 =t7 +rt;
		it =t20;	t20=t8 -it;		t8 =t8 +it;

		a[j1    ]=t1+t7;			a[j2    ]=t2+t8;
		a[j1+p1 ]=t1-t7;			a[j2+p1 ]=t2-t8;

		a[j1+p2 ]=t13-t20;			a[j2+p2 ]=t14+t19;	/* mpy by I is inlined here...	*/
		a[j1+p3 ]=t13+t20;			a[j2+p3 ]=t14-t19;

/*...Block 2: t3,9,15,21	*/
		rt =t15;	t15=t3 -rt;		t3 =t3 +rt;
		it =t16;	t16=t4 -it;		t4 =t4 +it;

		rt =t21;	t21=t9 -rt;		t9 =t9 +rt;
		it =t22;	t22=t10-it;		t10=t10+it;

		a[j1+p9 ]=t3+t9;			a[j2+p9 ]=t4+t10;
		a[j1+p8 ]=t3-t9;			a[j2+p8 ]=t4-t10;

		a[j1+p11]=t15-t22;			a[j2+p11]=t16+t21;	/* mpy by I is inlined here...	*/
		a[j1+p10]=t15+t22;			a[j2+p10]=t16-t21;

/*...Block 3: t5,11,17,23	*/
		rt =t17;	t17=t5 -rt;		t5 =t5 +rt;
		it =t18;	t18=t6 -it;		t6 =t6 +it;

		rt =t23;	t23=t11-rt;		t11=t11+rt;
		it =t24;	t24=t12-it;		t12=t12+it;

		a[j1+p6 ]=t5+t11;			a[j2+p6 ]=t6+t12;
		a[j1+p7 ]=t5-t11;			a[j2+p7 ]=t6-t12;

		a[j1+p5 ]=t17-t24;			a[j2+p5 ]=t18+t23;	/* mpy by I is inlined here...	*/
		a[j1+p4 ]=t17+t24;			a[j2+p4 ]=t18-t23;
								/* Totals: 96 FADD, 16 FMUL	*/
#else
/*
Alternatively, can use the 3-decrement version:
indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11
      -> 0, 3, 6, 9, 8,11,14,17, 4, 7,10,13
      == 0, 3, 6, 9, 8,11, 2, 5, 4, 7,10, 1 modulo 12.
I.e. start out with first triplet of indices {0,4,8}, permute those according to
{0,4,8}*11%12 = {0,8,4}, then each is head of a length-4 list of indices with decrement 3.

This requires us to swap output pairs 4/5 and 10/11 w.r.to the increment version.
*/
		t1 =a[j1    ];				t2 =a[j2    ];
		t3 =a[j1+p8 ]+a[j1+p4 ];	t4 =a[j2+p8 ]+a[j2+p4 ];
		t5 =a[j1+p8 ]-a[j1+p4 ];	t6 =a[j2+p8 ]-a[j2+p4 ];
		t1 =t1+t3;					t2 =t2+t4;
		t3 =t1+c4m1*t3;				t4 =t2+c4m1*t4;
		rt =c*t5;					it =c*t6;
		t5 =t3+it;					t6 =t4-rt;
		t3 =t3-it;					t4 =t4+rt;

		t7 =a[j1+p3 ];				t8 =a[j2+p3 ];
		t9 =a[j1+p11]+a[j1+p7 ];	t10=a[j2+p11]+a[j2+p7 ];
		t11=a[j1+p11]-a[j1+p7 ];	t12=a[j2+p11]-a[j2+p7 ];
		t7 =t7+t9;					t8 =t8+t10;
		t9 =t7+c4m1*t9;				t10=t8+c4m1*t10;
		rt =c*t11;					it =c*t12;
		t11=t9+it;					t12=t10-rt;
		t9 =t9-it;					t10=t10+rt;

		t13=a[j1+p6 ];				t14=a[j2+p6 ];
		t15=a[j1+p2 ]+a[j1+p10];	t16=a[j2+p2 ]+a[j2+p10];
		t17=a[j1+p2 ]-a[j1+p10];	t18=a[j2+p2 ]-a[j2+p10];
		t13=t13+t15;				t14=t14+t16;
		t15=t13+c4m1*t15;			t16=t14+c4m1*t16;
		rt =c*t17;					it =c*t18;
		t17=t15+it;					t18=t16-rt;
		t15=t15-it;					t16=t16+rt;

		t19=a[j1+p9 ];				t20=a[j2+p9 ];
		t21=a[j1+p5 ]+a[j1+p1 ];	t22=a[j2+p5 ]+a[j2+p1 ];
		t23=a[j1+p5 ]-a[j1+p1 ];	t24=a[j2+p5 ]-a[j2+p1 ];
		t19=t19+t21;				t20=t20+t22;
		t21=t19+c4m1*t21;			t22=t20+c4m1*t22;
		rt =c*t23;					it =c*t24;
		t23=t21+it;					t24=t22-rt;
		t21=t21-it;					t22=t22+rt;
/*
!...and now do three radix-4 transforms:
*/
/*...Block 1: t1,7,13,19	*/
			rt =t13;	t13=t1 -rt;	t1 =t1 +rt;
			it =t14 ;	t14=t2 -it;	t2 =t2 +it;

			rt =t19;	t19=t7 -rt;	t7 =t7 +rt;
			it =t20;	t20=t8 -it;	t8 =t8 +it;

			a[j1    ]=t1+t7;		a[j2    ]=t2+t8;
			a[j1+p1 ]=t1-t7;		a[j2+p1 ]=t2-t8;

			a[j1+p2 ]=t13-t20;		a[j2+p2 ]=t14+t19;	/* mpy by I is inlined here...	*/
			a[j1+p3 ]=t13+t20;		a[j2+p3 ]=t14-t19;

/*...Block 2: t3,9,15,21	*/
			rt =t15;	t15=t3 -rt;	t3 =t3 +rt;
			it =t16;	t16=t4 -it;	t4 =t4 +it;

			rt =t21;	t21=t9 -rt;	t9 =t9 +rt;
			it =t22;	t22=t10-it;	t10=t10+it;

			a[j1+p9 ]=t3+t9;		a[j2+p9 ]=t4+t10;
			a[j1+p8 ]=t3-t9;		a[j2+p8 ]=t4-t10;

			a[j1+p10]=t15-t22;		a[j2+p10]=t16+t21;	/* mpy by I is inlined here...	*/
			a[j1+p11]=t15+t22;		a[j2+p11]=t16-t21;

/*...Block 3: t5,11,17,23	*/
			rt =t17;	t17=t5 -rt;	t5 =t5 +rt;
			it =t18;	t18=t6 -it;	t6 =t6 +it;

			rt =t23;	t23=t11-rt;	t11=t11+rt;
			it =t24;	t24=t12-it;	t12=t12+it;

			a[j1+p6 ]=t5+t11;		a[j2+p6 ]=t6+t12;
			a[j1+p7 ]=t5-t11;		a[j2+p7 ]=t6-t12;

			a[j1+p4 ]=t17-t24;		a[j2+p4 ]=t18+t23;	/* mpy by I is inlined here...	*/
			a[j1+p5 ]=t17+t24;		a[j2+p5 ]=t18-t23;
								/* Totals: 96 FADD, 16 FMUL	*/
#endif
	}
}

/***************/

void radix12_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-12 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Uses an optimized [radix-3] x [radix-4] transform.
*/
	int j,j1,j2;
	static int n12,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, first_entry=TRUE;
	static double c = .86602540378443864676, c4m1 = -1.5;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24;

	if(!first_entry && (n/12) != n12)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n12=n/12;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n12;
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

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-12 pass is here.	*/

      for(j=0; j < n12; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*...gather the needed data (12 64-bit complex, i.e. 24 64-bit reals) and do four radix-3 transforms...	*/

	/*
	Twiddleless version requires us to swap inputs x3 <-> x9, x1 <-> x4, x7 <-> x10, x2 <-> x8.
	indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11
		  -> 0, 4, 8, 9, 1, 5, 6,10, 2, 3, 7,11 modulo 12.
	I.e. start out with first quartet of indices {0,3,6,9}, permute those according to
	  {0,3,6,9}*11%12 = {0,9,6,3}, then each is head of a length-3 list of indices with increment 4, so it seems
	that we can use either decrement 4 (see version C of this routine) or increment 4, whichever minimizes the number of index swaps.

	Remember, inputs to DIT are bit-reversed, so
	a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11] contain
	x[0, 6, 3, 9, 1, 7, 4,10, 2, 8, 5,11], which get swapped to
	x[0, 6, 9, 3, 4,10, 1, 7, 8, 2, 5,11], which means the a-indices get swapped as
	a[0, 1, 3, 2, 6, 7, 4, 5, 9, 8,10,11],
	i.e.
	swapping x3 and x9  means swapping a[2] and a[3];
	swapping x1 and x4  means swapping a[4] and a[6];
	swapping x7 and x10 means swapping a[5] and a[7];
	swapping x2 and x8  means swapping a[8] and a[9].
	*/
	/*...Block 1:	*/
			t1 =a[j1    ];		t2 =a[j2    ];
			rt =a[j1+p1 ];		it =a[j2+p1 ];
			t3 =t1 -rt;  		t1 =t1 +rt;
			t4 =t2 -it;			t2 =t2 +it;

			t5 =a[j1+p3 ];		t6 =a[j2+p3 ];
			rt =a[j1+p2 ];		it =a[j2+p2 ];
			t7 =t5 -rt;  		t5 =t5 +rt;
			t8 =t6 -it;  		t6 =t6 +it;

			rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
			it =t6;	t6 =t2 -it;	t2 =t2 +it;

			rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
					t8 =t4 +rt;	t4 =t4 -rt;

	/*...Block 2:	*/
			t9 =a[j1+p6 ];		t10=a[j2+p6 ];
			rt =a[j1+p7 ];		it =a[j2+p7 ];
			t11=t9 -rt;  		t9 =t9 +rt;
			t12=t10-it;			t10=t10+it;

			t13=a[j1+p4 ];		t14=a[j2+p4 ];
			rt =a[j1+p5 ];		it =a[j2+p5 ];
			t15=t13-rt;  		t13=t13+rt;
			t16=t14-it;			t14=t14+it;

			rt =t13;t13=t9 -rt;	t9 =t9 +rt;
			it =t14;t14=t10-it;	t10=t10+it;

			rt =t15;t15=t11-t16;t11=t11+t16;
					t16=t12+rt;	t12=t12-rt;

	/*...Block 3:	*/
			t17=a[j1+p9 ];		t18=a[j2+p9 ];
			rt =a[j1+p8 ];		it =a[j2+p8 ];
			t19=t17-rt;  		t17=t17+rt;
			t20=t18-it;			t18=t18+it;

			t21=a[j1+p10];		t22=a[j2+p10];
			rt =a[j1+p11];		it =a[j2+p11];
			t23=t21-rt;  		t21=t21+rt;
			t24=t22-it;			t22=t22+it;

			rt =t21;t21=t17-rt;	t17=t17+rt;
			it =t22;t22=t18-it;	t18=t18+it;

			rt =t23;t23=t19-t24;t19=t19+t24;
					t24=t20+rt;	t20=t20-rt;
	/*
	!...and now do four radix-3 transforms.
	*/
	/*...Block 1: t1,9,17	*/
			rt =t9;				it =t10;
			t9 =rt+t17;			t10=it+t18;
			t17=rt-t17;			t18=it-t18;
			t1 =t1+t9;			t2 =t2+t10;
			a[j1    ]=t1;		a[j2    ]=t2;
			t9 =t1+c4m1*t9;		t10=t2+c4m1*t10;
			rt =c*t17;			it =c*t18;
			a[j1+p4 ]=t9+it;	a[j2+p4 ]=t10-rt;
			a[j1+p8 ]=t9-it;	a[j2+p8 ]=t10+rt;
	/*...Block 2: t3,11,19	*/
			rt =t11;			it =t12;
			t11=rt+t19;			t12=it+t20;
			t19=rt-t19;			t20=it-t20;
			t3 =t3+t11;			t4 =t4+t12;
			a[j1+p3 ]=t3;		a[j2+p3 ]=t4;
			t11=t3+c4m1*t11;	t12=t4+c4m1*t12;
			rt =c*t19;			it =c*t20;
			a[j1+p7 ]=t11+it;	a[j2+p7 ]=t12-rt;
			a[j1+p11]=t11-it;	a[j2+p11]=t12+rt;
	/*...Block 3: t5,13,21	*/
			rt =t13;			it =t14;
			t13=rt+t21;			t14=it+t22;
			t21=rt-t21;			t22=it-t22;
			t5 =t5+t13;			t6 =t6+t14;
			a[j1+p6 ]=t5;		a[j2+p6 ]=t6;
			t13=t5+c4m1*t13;	t14=t6+c4m1*t14;
			rt =c*t21;			it =c*t22;
			a[j1+p10]=t13+it;	a[j2+p10]=t14-rt;
			a[j1+p2 ]=t13-it;	a[j2+p2 ]=t14+rt;
	/*...Block 3: t7,15,23	*/
			rt =t15;			it =t16;
			t15=rt+t23;			t16=it+t24;
			t23=rt-t23;			t24=it-t24;
			t7 =t7+t15;			t8 =t8+t16;
			a[j1+p9 ]=t7;		a[j2+p9 ]=t8;
			t15=t7+c4m1*t15;	t16=t8+c4m1*t16;
			rt =c*t23;			it =c*t24;
			a[j1+p1 ]=t15+it;	a[j2+p1 ]=t16-rt;
			a[j1+p5 ]=t15-it;	a[j2+p5 ]=t16+rt;
							/* Totals: 96 FADD, 16 FMUL.	*/
	}
}

