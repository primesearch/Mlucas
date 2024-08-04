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

/***************/

int radix11_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-11 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-14 complex DIF pass on the data in
!   the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n11,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10
		,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
#if LO_ADD
	const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
#else
/*...Fast length-5 cyclic convolution scheme needs the following: */
	const double a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -0.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5	*/

			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
	/* Here's how to use bc -l to get the above:
        scale=20
        p=8*a(1)
        a = c(  p/11)
        b = c(2*p/11)
        c = c(3*p/11)
        d = c(4*p/11)
        e = c(5*p/11)
        (   a      -  d+  c-  e)
        (         b-  d+  c-  e)
        (-2*a-2*b+3*d-2*c+3*e)/5
        (-  a+  b-  d+  c      )
        (-  a+  b-  d      +  e)
        ( 3*a-2*b+3*d-2*c-2*e)/5
        (            -  d+  c      )
        (         b-  d            )
        (-  a-  b+4*d-  c-  e)/5
        (   a+  b+  d+  c+  e)/5
        a = s(  p/11)
        b = s(2*p/11)
        c = s(3*p/11)
        d = s(4*p/11)
        e = s(5*p/11)
	(   a      -  d+  c-  e)
	(      -  b-  d+  c-  e)
	(-2*a+2*b+3*d-2*c+3*e)/5
	(-  a-  b-  d+  c      )
	(-  a-  b-  d      +  e)
	( 3*a+2*b+3*d-2*c-2*e)/5
	(            -  d+  c      )
	(      -  b-  d            )
	(-  a+  b+4*d-  c-  e)/5
	(   a-  b+  d+  c+  e)/5
	*/
	double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10;
#endif
	static double radix_inv, n2inv;
	double rt,it
		,cr1,cr2,cr3,cr4,cr5,ci1,ci2,ci3,ci4,ci5
		,sr1,sr2,sr3,sr4,sr5,si1,si2,si3,si4,si5
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,cy10,temp,frac,scale;
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

/*...change n11 and n_div_wt to non-static to work around a gcc compiler bug. */
	n11   = n/11;
	n_div_nwt = n11 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n11)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/11 in radix11_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)11));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p1 = n11;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;

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

		bjmodnini=0;
		for(j=0; j < n11; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-11 final DIT pass is here.	*/

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
	bjmodn10= bjmodn9 +bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+11;
	co3=co2-11;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

		t1 =a[j1    ];			t2 =a[j2    ];

		t3 =a[j1+p1 ];			t4 =a[j2+p1 ];
		rt =a[j1+p10];			it =a[j2+p10];
		t21=t3 -rt;				t22=t4 -it;
		t3 =t3 +rt;				t4 =t4 +it;

		t5 =a[j1+p2 ];			t6 =a[j2+p2 ];
		rt =a[j1+p9 ];			it =a[j2+p9 ];
		t19=t5 -rt;				t20=t6 -it;
		t5 =t5 +rt;				t6 =t6 +it;

		t7 =a[j1+p3 ];			t8 =a[j2+p3 ];
		rt =a[j1+p8 ];			it =a[j2+p8 ];
		t17=t7 -rt;				t18=t8 -it;
		t7 =t7 +rt;				t8 =t8 +it;

		t9 =a[j1+p4 ];			t10=a[j2+p4 ];
		rt =a[j1+p7 ];			it =a[j2+p7 ];
		t15=t9 -rt;				t16=t10-it;
		t9 =t9 +rt;				t10=t10+it;

		t11=a[j1+p5 ];			t12=a[j2+p5 ];
		rt =a[j1+p6 ];			it =a[j2+p6 ];
		t13=t11-rt;				t14=t12-it;
		t11=t11+rt;				t12=t12+it;

#if LO_ADD
		aj1p0r = t1+t3+t5+t7+t9+t11;		aj1p0i = t2+t4+t6+t8+t10+t12;	/* X0	*/

		cr1= t1+cc1*t3+cc2*t5+cc3*t7+cc4*t9+cc5*t11;	ci1= t2+cc1*t4+cc2*t6+cc3*t8+cc4*t10+cc5*t12;	/* C1	*/
		cr2= t1+cc2*t3+cc4*t5+cc5*t7+cc3*t9+cc1*t11;	ci2= t2+cc2*t4+cc4*t6+cc5*t8+cc3*t10+cc1*t12;	/* C2	*/
		cr3= t1+cc3*t3+cc5*t5+cc2*t7+cc1*t9+cc4*t11;	ci3= t2+cc3*t4+cc5*t6+cc2*t8+cc1*t10+cc4*t12;	/* C3	*/
		cr4= t1+cc4*t3+cc3*t5+cc1*t7+cc5*t9+cc2*t11;	ci4= t2+cc4*t4+cc3*t6+cc1*t8+cc5*t10+cc2*t12;	/* C4	*/
		cr5= t1+cc5*t3+cc1*t5+cc4*t7+cc2*t9+cc3*t11;	ci5= t2+cc5*t4+cc1*t6+cc4*t8+cc2*t10+cc3*t12;	/* C5	*/

		sr1= ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1= ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
		sr2= ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2= ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
		sr3= ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3= ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
		sr4= ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4= ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
		sr5= ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5= ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/
#else
/*  Here are the 5 cosine terms:
    let y0 = (x1+x10) = t3,4, y4 = (x2+x9) = t5,6, y2 = (x3+x8) = t7,8, y3 = (x4+x7) = t9,10, y1 = (x5+x6) = t11,12, then form	*/
		c1 = t3-t5;				s1 = t4-t6;
		c2 = t11-t5;			s2 = t12-t6;
		c3 = c1+c2;				s3 = s1+s2;
		c4 = t7-t5;				s4 = t8-t6;
		c5 = t9-t5;				s5 = t10-t6;
		c6 = c4+c5;				s6 = s4+s5;
		c7 = c1-c4;				s7 = s1-s4;
		c8 = c2-c5;				s8 = s2-s5;
		c9 = c3-c6;				s9 = s3-s6;
		c10= t3+t5+t7+t9+t11;	s10= t4+t6+t8+t10+t12;

		aj1p0r = t1+c10;		aj1p0i = t2+s10;	/* X0	*/

		c3 = a2*c3;				s3 = a2*s3;
		c6 = a5*c6;				s6 = a5*s6;
		c9 = a8*c9;				s9 = a8*s9;
		c10= a9*c10+t1;			s10= a9*s10+t2;

		c1 = a0*c1+c3;			s1 = a0*s1+s3;
		c2 = a1*c2+c3;			s2 = a1*s2+s3;
		c3 = a3*c4+c6;			s3 = a3*s4+s6;
		c4 = a4*c5+c6;			s4 = a4*s5+s6;
		c5 = a6*c7+c9;			s5 = a6*s7+s9;
		c6 = a7*c8+c9;			s6 = a7*s8+s9;

		cr1 = c10+c1-c5;		ci1 = s10+s1-s5;
		cr2 = c10-c1-c2-c3-c4;	ci2 = s10-s1-s2-s3-s4;
		cr3 = c10+c3+c5;		ci3 = s10+s3+s5;
		cr4 = c10+c4+c6;		ci4 = s10+s4+s6;
		cr5 = c10+c2-c6;		ci5 = s10+s2-s6;
/*  Here are the 5 sine terms:
    let y0 = (x1-x10) = t21,22, y4 = (x9-x2) = -t19,20, y2 = (x3-x8) = t17,18, y3 = (x4-x7) = t15,16, y1 = (x5-x6) = t13,t14, then form	*/
		c1 = t21+t19;			s1 = t22+t20;
		c2 = t13+t19;			s2 = t14+t20;
		c3 = c1+c2;				s3 = s1+s2;
		c4 = t17+t19;			s4 = t18+t20;
		c5 = t15+t19;			s5 = t16+t20;
		c6 = c4+c5;				s6 = s4+s5;
		c7 = c1-c4;				s7 = s1-s4;
		c8 = c2-c5;				s8 = s2-s5;
		c9 = c3-c6;				s9 = s3-s6;
		c10= t21-t19+t17+t15+t13;	s10= t22-t20+t18+t16+t14;

		c3 = b2*c3;				s3 = b2*s3;
		c6 = b5*c6;				s6 = b5*s6;
		c9 = b8*c9;				s9 = b8*s9;
		c10= b9*c10;			s10= b9*s10;

		c1 = b0*c1+c3;			s1 = b0*s1+s3;
		c2 = b1*c2+c3;			s2 = b1*s2+s3;
		c3 = b3*c4+c6;			s3 = b3*s4+s6;
		c4 = b4*c5+c6;			s4 = b4*s5+s6;
		c5 = b6*c7+c9;			s5 = b6*s7+s9;
		c6 = b7*c8+c9;			s6 = b7*s8+s9;

		sr1 = c10+c1-c5;		si1 = s10+s1-s5;
		sr2 = c1+c2+c3+c4-c10;	si2 = s1+s2+s3+s4-s10;
		sr3 = c10+c3+c5;		si3 = s10+s3+s5;
		sr4 = c10+c4+c6;		si4 = s10+s4+s6;
		sr5 = c10+c2-c6;		si5 = s10+s2-s6;
#endif
/*...Inline multiply of sine parts by +-I into finishing phase...	*/

		aj1p1r =cr1+si1;		aj1p1i =ci1-sr1;
		aj1p2r =cr2+si2;		aj1p2i =ci2-sr2;
		aj1p3r =cr3+si3;		aj1p3i =ci3-sr3;
		aj1p4r =cr4+si4;		aj1p4i =ci4-sr4;
		aj1p5r =cr5+si5;		aj1p5i =ci5-sr5;
		aj1p6r =cr5-si5;		aj1p6i =ci5+sr5;
		aj1p7r =cr4-si4;		aj1p7i =ci4+sr4;
		aj1p8r =cr3-si3;		aj1p8i =ci3+sr3;
		aj1p9r =cr2-si2;		aj1p9i =ci2+sr2;
		aj1p10r=cr1-si1;		aj1p10i=ci1+sr1;

/*..Now do the carries. Since the outputs would
    normally be getting dispatched to 11 separate blocks of the A-array, we need 11 separate carries.	*/

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
	   cmplx_carry_norm_errcheck0(aj1p0r ,aj1p0i ,cy0 ,bjmodn0 ,0 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p1r ,aj1p1i ,cy1 ,bjmodn1 ,1 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p2r ,aj1p2i ,cy2 ,bjmodn2 ,2 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p3r ,aj1p3i ,cy3 ,bjmodn3 ,3 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p4r ,aj1p4i ,cy4 ,bjmodn4 ,4 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p5r ,aj1p5i ,cy5 ,bjmodn5 ,5 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p6r ,aj1p6i ,cy6 ,bjmodn6 ,6 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p7r ,aj1p7i ,cy7 ,bjmodn7 ,7 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p8r ,aj1p8i ,cy8 ,bjmodn8 ,8 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p9r ,aj1p9i ,cy9 ,bjmodn9 ,9 ,prp_mult);
		cmplx_carry_norm_errcheck(aj1p10r,aj1p10i,cy10,bjmodn10,10,prp_mult);

		i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
		co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-11 DIF pass is here:	*/

/*...gather the needed data (11 64-bit complex, i.e. 22 64-bit reals)...	*/

		t1 =aj1p0r;				t2 =aj1p0i;
		t3 =aj1p1r+aj1p10r;		t4 =aj1p1i+aj1p10i;
		t5 =aj1p2r+aj1p9r;		t6 =aj1p2i+aj1p9i;
		t7 =aj1p3r+aj1p8r;		t8 =aj1p3i+aj1p8i;
		t9 =aj1p4r+aj1p7r;		t10=aj1p4i+aj1p7i;
		t11=aj1p5r+aj1p6r;		t12=aj1p5i+aj1p6i;
		t13=aj1p5r-aj1p6r;		t14=aj1p5i-aj1p6i;
		t15=aj1p4r-aj1p7r;		t16=aj1p4i-aj1p7i;
		t17=aj1p3r-aj1p8r;		t18=aj1p3i-aj1p8i;
		t19=aj1p2r-aj1p9r;		t20=aj1p2i-aj1p9i;
		t21=aj1p1r-aj1p10r;		t22=aj1p1i-aj1p10i;
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif

#if LO_ADD
		a[j1   ] = t1+t3+t5+t7+t9+t11;				a[j2] = t2+t4+t6+t8+t10+t12;	/* X0	*/

		cr1= t1+cc1*t3+cc2*t5+cc3*t7+cc4*t9+cc5*t11;	ci1= t2+cc1*t4+cc2*t6+cc3*t8+cc4*t10+cc5*t12;	/* C1	*/
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
		cr2= t1+cc2*t3+cc4*t5+cc5*t7+cc3*t9+cc1*t11;	ci2= t2+cc2*t4+cc4*t6+cc5*t8+cc3*t10+cc1*t12;	/* C2	*/
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif
		cr3= t1+cc3*t3+cc5*t5+cc2*t7+cc1*t9+cc4*t11;	ci3= t2+cc3*t4+cc5*t6+cc2*t8+cc1*t10+cc4*t12;	/* C3	*/
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
		cr4= t1+cc4*t3+cc3*t5+cc1*t7+cc5*t9+cc2*t11;	ci4= t2+cc4*t4+cc3*t6+cc1*t8+cc5*t10+cc2*t12;	/* C4	*/
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
		cr5= t1+cc5*t3+cc1*t5+cc4*t7+cc2*t9+cc3*t11;	ci5= t2+cc5*t4+cc1*t6+cc4*t8+cc2*t10+cc3*t12;	/* C5	*/

		sr1= ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1= ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
		sr2= ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2= ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif
		sr3= ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3= ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
		sr4= ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4= ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif
		sr5= ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5= ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/
#else
/*  Here are the 5 cosine terms:
    let y0 = (x1+x10) = t3,4, y4 = (x2+x9) = t5,6, y2 = (x3+x8) = t7,8, y3 = (x4+x7) = t9,10, y1 = (x5+x6) = t11,12, then form	*/
		c1 = t3-t5;			s1 = t4-t6;
		c2 = t11-t5;			s2 = t12-t6;
		c3 = c1+c2;			s3 = s1+s2;
		c4 = t7-t5;			s4 = t8-t6;
		c5 = t9-t5;			s5 = t10-t6;
		c6 = c4+c5;			s6 = s4+s5;
		c7 = c1-c4;			s7 = s1-s4;
		c8 = c2-c5;			s8 = s2-s5;
		c9 = c3-c6;			s9 = s3-s6;
		c10= t3+t5+t7+t9+t11;		s10= t4+t6+t8+t10+t12;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
		a[j1  ] = t1+c10;		a[j2] = t2+s10;	/* X0	*/

		c3 = a2*c3;			s3 = a2*s3;
		c6 = a5*c6;			s6 = a5*s6;
		c9 = a8*c9;			s9 = a8*s9;
		c10= a9*c10+t1;		s10= a9*s10+t2;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif
		c1 = a0*c1+c3;		s1 = a0*s1+s3;
		c2 = a1*c2+c3;		s2 = a1*s2+s3;
		c3 = a3*c4+c6;		s3 = a3*s4+s6;
		c4 = a4*c5+c6;		s4 = a4*s5+s6;
		c5 = a6*c7+c9;		s5 = a6*s7+s9;
		c6 = a7*c8+c9;		s6 = a7*s8+s9;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
		cr1 = c10+c1-c5;		ci1 = s10+s1-s5;
		cr2 = c10-c1-c2-c3-c4;	ci2 = s10-s1-s2-s3-s4;
		cr3 = c10+c3+c5;		ci3 = s10+s3+s5;
		cr4 = c10+c4+c6;		ci4 = s10+s4+s6;
		cr5 = c10+c2-c6;		ci5 = s10+s2-s6;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif

/*  Here are the 5 sine terms:
    let y0 = (x1-x10) = t21,22, y4 = (x9-x2) = -t19,20, y2 = (x3-x8) = t17,18, y3 = (x4-x7) = t15,16, y1 = (x5-x6) = t13,t14, then form	*/

		c1 = t21+t19;			s1 = t22+t20;
		c2 = t13+t19;			s2 = t14+t20;
		c3 = c1+c2;			s3 = s1+s2;
		c4 = t17+t19;			s4 = t18+t20;
		c5 = t15+t19;			s5 = t16+t20;
		c6 = c4+c5;			s6 = s4+s5;
		c7 = c1-c4;			s7 = s1-s4;
		c8 = c2-c5;			s8 = s2-s5;
		c9 = c3-c6;			s9 = s3-s6;
		c10= t21-t19+t17+t15+t13;	s10= t22-t20+t18+t16+t14;
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
		c3 = b2*c3;			s3 = b2*s3;
		c6 = b5*c6;			s6 = b5*s6;
		c9 = b8*c9;			s9 = b8*s9;
		c10= b9*c10;			s10= b9*s10;
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif
		c1 = b0*c1+c3;		s1 = b0*s1+s3;
		c2 = b1*c2+c3;		s2 = b1*s2+s3;
		c3 = b3*c4+c6;		s3 = b3*s4+s6;
		c4 = b4*c5+c6;		s4 = b4*s5+s6;
		c5 = b6*c7+c9;		s5 = b6*s7+s9;
		c6 = b7*c8+c9;		s6 = b7*s8+s9;
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
		sr1 = c10+c1-c5;		si1 = s10+s1-s5;
		sr2 = c1+c2+c3+c4-c10;	si2 = s1+s2+s3+s4-s10;
		sr3 = c10+c3+c5;		si3 = s10+s3+s5;
		sr4 = c10+c4+c6;		si4 = s10+s4+s6;
		sr5 = c10+c2-c6;		si5 = s10+s2-s6;
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif

#endif

/*...Inline multiply of sine parts by +-I into finishing phase...	*/

		a[j1+p1 ]=cr1-si1;		a[j2+p1 ]=ci1+sr1;
		a[j1+p2 ]=cr2-si2;		a[j2+p2 ]=ci2+sr2;
		a[j1+p3 ]=cr3-si3;		a[j2+p3 ]=ci3+sr3;
		a[j1+p4 ]=cr4-si4;		a[j2+p4 ]=ci4+sr4;
		a[j1+p5 ]=cr5-si5;		a[j2+p5 ]=ci5+sr5;
#if PFETCH
addr = add0+p9;
prefetch_p_doubles(addr);
#endif
		a[j1+p6 ]=cr5+si5;		a[j2+p6 ]=ci5-sr5;
		a[j1+p7 ]=cr4+si4;		a[j2+p7 ]=ci4-sr4;
		a[j1+p8 ]=cr3+si3;		a[j2+p8 ]=ci3-sr3;
		a[j1+p9 ]=cr2+si2;		a[j2+p9 ]=ci2-sr2;
		a[j1+p10]=cr1+si1;		a[j2+p10]=ci1-sr1;
#if PFETCH
addr = add0+p10;
prefetch_p_doubles(addr);
#endif
		}

		jstart += nwt;
		jhi    += nwt;
		col += 11;
		co3 -= 11;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-11 forward DIF FFT of the first block of 11 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 11 outputs of (1);
!   (3) Reweight and perform a radix-11 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 11 elements and repeat (1-4).
*/
	t1  = cy10;
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
		a[j+p10] *= radix_inv;
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9)+fabs(cy10) != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix11_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix11_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-11 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Given complex inputs (x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10), we need the following outputs
!   (here cJ = cos(2*J*pi/11), sJ = sin(2*J*pi/11)):
!
!	X0 = C0,          where C0 = x0+   (x1+x10)+   (x2+x9)+   (x3+x8)+   (x4+x7)+   (x5+x6),
!
!	X1 = C1 + I*S1          C1 = x0+c1*(x1+x10)+c2*(x2+x9)+c3*(x3+x8)+c4*(x4+x7)+c5*(x5+x6)
!	X2 = C2 + I*S2          C2 = x0+c2*(x1+x10)+c4*(x2+x9)+c5*(x3+x8)+c3*(x4+x7)+c1*(x5+x6)
!	X3 = C3 + I*S3          C3 = x0+c3*(x1+x10)+c5*(x2+x9)+c2*(x3+x8)+c1*(x4+x7)+c4*(x5+x6)
!	X4 = C4 + I*S4          C4 = x0+c4*(x1+x10)+c3*(x2+x9)+c1*(x3+x8)+c5*(x4+x7)+c2*(x5+x6)
!	X5 = C5 + I*S5          C5 = x0+c5*(x1+x10)+c1*(x2+x9)+c4*(x3+x8)+c2*(x4+x7)+c3*(x5+x6),
!												To make all signs positive:
!	X6 = C5 - I*S5          S1 =    s1*(x1-x10)+s2*(x2-x9)+s3*(x3-x8)+s4*(x4-x7)+s5*(x5-x6) I=2 = 01000 (lsb first):
!	X7 = C4 - I*S4          S2 =    s2*(x1-x10)+s4*(x2-x9)-s5*(x3-x8)-s3*(x4-x7)-s1*(x5-x6) -change sign of s2, (x2-x9)
!	X8 = C3 - I*S3          S3 =    s3*(x1-x10)-s5*(x2-x9)-s2*(x3-x8)+s1*(x4-x7)+s4*(x5-x6) -change sign of entire row 2
!	X9 = C2 - I*S2          S4 =    s4*(x1-x10)-s3*(x2-x9)+s1*(x3-x8)+s5*(x4-x7)-s2*(x5-x6) I=29= 10111 (lsb first):
!	X10= C1 - I*S1          S5 =    s5*(x1-x10)-s1*(x2-x9)+s4*(x3-x8)-s2*(x4-x7)+s3*(x5-x6).-simply the 31-comp. of I=2.
!
!   We refer to the terms C1,2,3,4,5 (which do not explicitly involving the imaginary constant I)
!   as the "cosine part" of the output, and S1,2,3,4,5 (those multiplied by I) as the "sine part."
!											opcount for general odd-prime radix R:
!   Form      (x1+-x10),  (x2+-x9),  (x3+-x8),  (x4+-x7),  (x5+-x6) :  0 FMUL, 20 FADD		0 fmul	2*(R-1)       fadd
!   Form X0                                                         :  0 FMUL, 10 FADD			(R-1)         fadd
!   Form x0+c1*(x1+x10)+c2*(x2+x9)+c3*(x3+x8)+c4*(x4+x7)+c5*(x5+x6) : 10 FMUL, 10 FADD	(R-1)^2/2 fmul	(R-1)*(R-1)/2 fadd
!   Form x0+c2*(x1+x10)+c4*(x2+x9)+c5*(x3+x8)+c3*(x4+x7)+c1*(x5+x6) : 10 FMUL, 10 FADD
!   Form x0+c3*(x1+x10)+c5*(x2+x9)+c2*(x3+x8)+c1*(x4+x7)+c4*(x5+x6) : 10 FMUL, 10 FADD
!   Form x0+c4*(x1+x10)+c3*(x2+x9)+c1*(x3+x8)+c5*(x4+x7)+c2*(x5+x6) : 10 FMUL, 10 FADD
!   Form x0+c5*(x1+x10)+c1*(x2+x9)+c4*(x3+x8)+c2*(x4+x7)+c3*(x5+x6) : 10 FMUL, 10 FADD
!
!   Form    s1*(x1-x10)+s2*(x2-x9)+s3*(x3-x8)+s4*(x4-x7)+s5*(x5-x6) : 10 FMUL,  8 FADD	(R-1)^2/2 fmul	(R-3)*(R-1)/2 fadd
!   Form    s2*(x1-x10)+s4*(x2-x9)-s5*(x3-x8)-s3*(x4-x7)-s1*(x5-x6) : 10 FMUL,  8 FADD
!   Form    s3*(x1-x10)-s5*(x2-x9)-s2*(x3-x8)+s1*(x4-x7)+s4*(x5-x6) : 10 FMUL,  8 FADD
!   Form    s4*(x1-x10)-s3*(x2-x9)+s1*(x3-x8)+s5*(x4-x7)-s2*(x5-x6) : 10 FMUL,  8 FADD
!   Form    s5*(x1-x10)-s1*(x2-x9)+s4*(x3-x8)-s2*(x4-x7)+s3*(x5-x6) : 10 FMUL,  8 FADD
!   Form X1,2,3,4,5,6,7,8,9,10                                      :  0 FMUL, 20 FADD		0 fmul	2*(R-1)       fadd
!
!   Totals :                                                        100 FMUL, 140 FADD,		(R-1)^2 fmul	(R+3)*(R-1) fadd
!                                                        compared to 16 FMUL,  96 FADD for radix-12. (Ouch!)
!
!   Relative cost := #FADD/(radix*lg2(radix)) = 3.679 .
!
!   This routine gets called only infrequently, so no optimizations beyond exploitation of the obvious symmetries are done.
!   On the other hand, the radix11_ditN_cy_dif1 subroutine is called once every iteration, so there we use an optimized
!   version which needs just 40 FMUL and 168 FADD.
*/
	int j,j1,j2;
	static int n11,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, first_entry=TRUE;
#if LO_ADD
	const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
#else
/*...Fast length-5 cyclic convolution scheme needs the following: */
	const double a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1; This helps save a register over the "old" scheme (discovered while mapping out the SSE2 version of the radix-11 DFT) */

			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
	/* Here's how to use bc -l to get the above:
        scale=20
        p=8*a(1)
        a = c(  p/11)
        b = c(2*p/11)
        c = c(3*p/11)
        d = c(4*p/11)
        e = c(5*p/11)
        (   a      -  d+  c-  e)
        (         b-  d+  c-  e)
        (-2*a-2*b+3*d-2*c+3*e)/5
        (-  a+  b-  d+  c      )
        (-  a+  b-  d      +  e)
        ( 3*a-2*b+3*d-2*c-2*e)/5
        (            -  d+  c      )
        (         b-  d            )
        (-  a-  b+4*d-  c-  e)/5
        (   a+  b+  d+  c+  e)/5
        a = s(  p/11)
        b = s(2*p/11)
        c = s(3*p/11)
        d = s(4*p/11)
        e = s(5*p/11)
	(   a      -  d+  c-  e)
	(      -  b-  d+  c-  e)
	(-2*a+2*b+3*d-2*c+3*e)/5
	(-  a-  b-  d+  c      )
	(-  a-  b-  d      +  e)
	( 3*a+2*b+3*d-2*c-2*e)/5
	(            -  d+  c      )
	(      -  b-  d            )
	(-  a+  b+4*d-  c-  e)/5
	(   a-  b+  d+  c+  e)/5
	*/
#endif

	if(!first_entry && (n/11) != n11)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n11=n/11;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n11;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;

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
	}

/*...The radix-11 pass is here.	*/

	for(j=0; j < n11; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
			j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

#if LO_ADD
		RADIX_11_DFT_BASIC(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10]
					,a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);
#else
		RADIX_11_DFT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10]
					,a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);
#endif
	}
}

/***************/

void radix11_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform a final radix-11 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n11,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, first_entry=TRUE;
#if LO_ADD
	const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
#else
/*...Fast length-5 cyclic convolution scheme needs the following: */
	const double a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1; This helps save a register over the "old" scheme (discovered while mapping out the SSE2 version of the radix-11 DFT) */

			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
#endif

	if(!first_entry && (n/11) != n11)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n11=n/11;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n11;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;

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
	}

/*...The radix-11 pass is here.	*/

	for(j=0; j < n11; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
		/* Call same radix-11 DFT macro as for DIF, but replace indices [0,1,2,3,4,5,6,7,8,9,10] with j*10%11, j = 0, ..., 10: */
#if LO_ADD
		RADIX_11_DFT_BASIC(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10]
					,a[j1],a[j2],a[j1+p10],a[j2+p10],a[j1+p9],a[j2+p9],a[j1+p8],a[j2+p8],a[j1+p7],a[j2+p7],a[j1+p6],a[j2+p6],a[j1+p5],a[j2+p5],a[j1+p4],a[j2+p4],a[j1+p3],a[j2+p3],a[j1+p2],a[j2+p2],a[j1+p1],a[j2+p1], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);
#else
		RADIX_11_DFT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10]
					,a[j1],a[j2],a[j1+p10],a[j2+p10],a[j1+p9],a[j2+p9],a[j1+p8],a[j2+p8],a[j1+p7],a[j2+p7],a[j1+p6],a[j2+p6],a[j1+p5],a[j2+p5],a[j1+p4],a[j2+p4],a[j1+p3],a[j2+p3],a[j1+p2],a[j2+p2],a[j1+p1],a[j2+p1], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);
#endif
	}
}

