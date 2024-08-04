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

int radix18_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-18 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-18 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n18,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17
		,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double c   =  0.76604444311897803520,	/* cos(2*pi/9) */
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
	,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r,aj1p16r,aj1p17r
	,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i,aj1p16i,aj1p17i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,temp,frac,scale;
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
	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	if((TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change n18 and n_div_wt to non-static to work around a gcc compiler bug. */
	n18   = n/18;
	n_div_nwt = n18 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n18)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/18 in radix18_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)18));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p1 = n18;
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
		p16= p15+p1;
		p17= p16+p1;

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
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
		p15= p15+ ( (p15>> DAT_BITS) << PAD_BITS );
		p16= p16+ ( (p16>> DAT_BITS) << PAD_BITS );
		p17= p17+ ( (p17>> DAT_BITS) << PAD_BITS );

		bjmodnini=0;
		for(j=0; j < n18; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-18 final DIT pass is here.	*/

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
	cy12= 0;
	cy13= 0;
	cy14= 0;
	cy15= 0;
	cy16= 0;
	cy17= 0;

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
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);
	bjmodn12= bjmodn11+bjmodnini-n; bjmodn12= bjmodn12+ ( (-(int)((uint32)bjmodn12>> 31)) & n);
	bjmodn13= bjmodn12+bjmodnini-n; bjmodn13= bjmodn13+ ( (-(int)((uint32)bjmodn13>> 31)) & n);
	bjmodn14= bjmodn13+bjmodnini-n; bjmodn14= bjmodn14+ ( (-(int)((uint32)bjmodn14>> 31)) & n);
	bjmodn15= bjmodn14+bjmodnini-n; bjmodn15= bjmodn15+ ( (-(int)((uint32)bjmodn15>> 31)) & n);
	bjmodn16= bjmodn15+bjmodnini-n; bjmodn16= bjmodn16+ ( (-(int)((uint32)bjmodn16>> 31)) & n);
	bjmodn17= bjmodn16+bjmodnini-n; bjmodn17= bjmodn17+ ( (-(int)((uint32)bjmodn17>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+18;
	co3=co2-18;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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
	!...gather the needed data (18 64-bit complex, i.e. 36 64-bit reals) and do a radix-18 DIT transform...
	*/
			t00=a[j1    ];		t01=a[j2    ];
			rt =a[j1+p1 ];		it =a[j2+p1 ];
			t10=t00-rt;			t11=t01-it;
			t00=t00+rt;			t01=t01+it;

			t02=a[j1+p3 ];		t03=a[j2+p3 ];
			rt =a[j1+p2 ];		it =a[j2+p2 ];
			t12=t02-rt;			t13=t03-it; ;
			t02=t02+rt;			t03=t03+it; ;

			t04=a[j1+p4 ];		t05=a[j2+p4 ];
			rt =a[j1+p5 ];		it =a[j2+p5 ];
			t14=t04-rt;			t15=t05-it; ;
			t04=t04+rt;			t05=t05+it; ;

			t06=a[j1+p11];		t07=a[j2+p11];
			rt =a[j1+p10];		it =a[j2+p10];
			t16=t06-rt;			t17=t07-it; ;
			t06=t06+rt;			t07=t07+it; ;

			t08=a[j1+p7 ];		t09=a[j2+p7 ];
			rt =a[j1+p6 ];		it =a[j2+p6 ];
			t18=t08-rt;			t19=t09-it; ;
			t08=t08+rt;			t09=t09+it; ;

			t0a=a[j1+p8 ];		t0b=a[j2+p8 ];
			rt =a[j1+p9 ];		it =a[j2+p9 ];
			t1a=t0a-rt;			t1b=t0b-it; ;
			t0a=t0a+rt;			t0b=t0b+it; ;

			t0c=a[j1+p15];		t0d=a[j2+p15];
			rt =a[j1+p14];		it =a[j2+p14];
			t1c=t0c-rt;			t1d=t0d-it; ;
			t0c=t0c+rt;			t0d=t0d+it; ;

			t0e=a[j1+p16];		t0f=a[j2+p16];
			rt =a[j1+p17];		it =a[j2+p17];
			t1e=t0e-rt;			t1f=t0f-it; ;
			t0e=t0e+rt;			t0f=t0f+it; ;

			t0g=a[j1+p12];		t0h=a[j2+p12];
			rt =a[j1+p13];		it =a[j2+p13];
			t1g=t0g-rt;			t1h=t0h-it;
			t0g=t0g+rt;			t0h=t0h+it;

	/*       ...and now do two radix-9 transforms.	*/
	/*...t0[0:16:2]r use t[00:0g:2]; t0[1:17:2]i use t[10:1g:2]	*/

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
			aj1p0r =t00;		aj1p0i =t01;
			t06 =t00+c3m1*t06;	t07 =t01+c3m1*t07;
			rt  =s3*t0c;		it  =s3*t0d;
			aj1p6r =t06+it;		aj1p6i =t07-rt;
			aj1p12r=t06-it;		aj1p12i=t07+rt;

			rt  =t08*c +t09*s;	it  =t09*c -t08*s;
			re  =t0e*c2+t0f*s2;	t0f =t0f*c2-t0e*s2;	t0e=re;
			t08 =rt+t0e;		t09 =it+t0f;
			t0e =rt-t0e;		t0f =it-t0f;
			t02 =t02+t08;		t03 =t03+t09;
			aj1p8r =t02;		aj1p8i =t03;
			t08 =t02+c3m1*t08;	t09=t03+c3m1*t09;
			rt  =s3*t0e;		it  =s3*t0f;
			aj1p14r=t08+it;		aj1p14i=t09-rt;
			aj1p2r =t08-it;		aj1p2i =t09+rt;

			rt  =t0a*c2+t0b*s2;	it  =t0b*c2-t0a*s2;
			re  =t0g*c4+t0h*s4;	t0h =t0h*c4-t0g*s4;	t0g=re;
			t0a =rt+t0g;		t0b =it+t0h;
			t0g =rt-t0g;		t0h =it-t0h;
			t04 =t04+t0a;		t05 =t05+t0b;
			aj1p16r=t04;		aj1p16i=t05;
			t0a =t04+c3m1*t0a;	t0b =t05+c3m1*t0b;
			rt  =s3*t0g;		it  =s3*t0h;
			aj1p4r =t0a+it;		aj1p4i =t0b-rt;
			aj1p10r=t0a-it;		aj1p10i=t0b+rt;

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
			aj1p9r =t10;		aj1p9i =t11;
			t16 =t10+c3m1*t16;	t17 =t11+c3m1*t17;
			rt  =s3*t1c;		it  =s3*t1d;
			aj1p15r=t16+it;		aj1p15i=t17-rt;
			aj1p3r =t16-it;		aj1p3i =t17+rt;

			rt  =t18*c +t19*s;	it  =t19*c -t18 *s;
			re  =t1e*c2+t1f*s2;	t1f =t1f*c2-t1e*s2;	t1e=re;
			t18 =rt+t1e;		t19 =it+t1f;
			t1e =rt-t1e;		t1f =it-t1f;
			t12 =t12+t18;		t13 =t13+t19;
			aj1p17r=t12;		aj1p17i=t13;
			t18 =t12+c3m1*t18;	t19=t13+c3m1*t19;
			rt  =s3*t1e;		it  =s3*t1f;
			aj1p5r =t18+it;		aj1p5i =t19-rt;
			aj1p11r=t18-it;		aj1p11i=t19+rt;

			rt  =t1a*c2+t1b*s2;	it  =t1b*c2-t1a*s2;
			re  =t1g*c4+t1h*s4;	t1h =t1h*c4-t1g*s4;	t1g=re;
			t1a =rt+t1g;		t1b =it+t1h;
			t1g =rt-t1g;		t1h =it-t1h;
			t14 =t14+t1a;		t15 =t15+t1b;
			aj1p7r =t14;		aj1p7i =t15;
			t1a =t14+c3m1*t1a;	t1b =t15+c3m1*t1b;
			rt  =s3*t1g;		it  =s3*t1h;
			aj1p13r=t1a+it;		aj1p13i=t1b-rt;
			aj1p1r =t1a-it;		aj1p1i =t1b+rt;

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 18 separate blocks of the A-array, we need 18 separate carries.	*/

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
			cmplx_carry_norm_errcheck(aj1p11r,aj1p11i,cy11,bjmodn11,11,prp_mult);
			cmplx_carry_norm_errcheck(aj1p12r,aj1p12i,cy12,bjmodn12,12,prp_mult);
			cmplx_carry_norm_errcheck(aj1p13r,aj1p13i,cy13,bjmodn13,13,prp_mult);
			cmplx_carry_norm_errcheck(aj1p14r,aj1p14i,cy14,bjmodn14,14,prp_mult);
			cmplx_carry_norm_errcheck(aj1p15r,aj1p15i,cy15,bjmodn15,15,prp_mult);
			cmplx_carry_norm_errcheck(aj1p16r,aj1p16i,cy16,bjmodn16,16,prp_mult);
			cmplx_carry_norm_errcheck(aj1p17r,aj1p17i,cy17,bjmodn17,17,prp_mult);

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-18 DIF pass is here:	*/
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif

/*...First radix-9 transform:	*/
		t00 =aj1p0r;			t01 =aj1p0i;
		t02 =aj1p12r+aj1p6r;	t03 =aj1p12i+aj1p6i ;
		t04 =aj1p12r-aj1p6r;	t05 =aj1p12i-aj1p6i ;
		t00 =t00+t02;			t01 =t01+t03;
		t02 =t00+c3m1*t02;		t03 =t01+c3m1*t03;
		rt  =s3*t04;			it  =s3*t05;
		t04 =t02+it;			t05 =t03-rt;
		t02 =t02-it;			t03 =t03+rt;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
		t06 =aj1p16r;			t07 =aj1p16i;
		t08 =aj1p10r+aj1p4r;	t09 =aj1p10i+aj1p4i ;
		t0a =aj1p10r-aj1p4r;	t0b =aj1p10i-aj1p4i ;
		t06 =t06+t08;			t07 =t07+t09;
		t08 =t06+c3m1*t08;		t09 =t07+c3m1*t09;
		rt  =s3*t0a;			it  =s3*t0b;
		t0a =t08+it;			t0b =t09-rt;
		t08 =t08-it;			t09 =t09+rt;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif
		t0c =aj1p14r;			t0d =aj1p14i;
		t0e =aj1p8r +aj1p2r;	t0f =aj1p8i +aj1p2i ;
		t0g =aj1p8r -aj1p2r;	t0h =aj1p8i -aj1p2i ;
		t0c =t0c+t0e;			t0d =t0d+t0f;
		t0e =t0c+c3m1*t0e;		t0f =t0d+c3m1*t0f;
		rt  =s3*t0g;			it  =s3*t0h;
		t0g =t0e+it;			t0h =t0f-rt;
		t0e =t0e-it;			t0f =t0f+rt;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif	/* Twiddles: */
		rt  =t06;				it  =t07;
		t06 =rt+t0c;			t07 =it+t0d;
		t0c =rt-t0c;			t0d =it-t0d;
		t00 =t00+t06;			t01 =t01+t07;
		t06 =t00+c3m1*t06;		t07 =t01+c3m1*t07;
		rt  =s3*t0c;			it  =s3*t0d;
		t0c =t06+it;			t0d =t07-rt;
		t06 =t06-it;			t07 =t07+rt;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
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
addr = add0+p5;
prefetch_p_doubles(addr);
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
addr = add0+p6;
prefetch_p_doubles(addr);
#endif
/*...Second radix-9 transform:	*/
		t10 =aj1p9r;			t11 =aj1p9i;
		t12 =aj1p3r +aj1p15r;	t13 =aj1p3i +aj1p15i;
		t14 =aj1p3r -aj1p15r;	t15 =aj1p3i -aj1p15i;
		t10 =t10+t12;			t11 =t11+t13;
		t12 =t10+c3m1*t12;		t13 =t11+c3m1*t13;
		rt  =s3*t14;			it  =s3*t15;
		t14 =t12+it;			t15 =t13-rt;
		t12 =t12-it;			t13 =t13+rt;
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
		t16 =aj1p7r;			t17 =aj1p7i;
		t18 =aj1p1r +aj1p13r;	t19 =aj1p1i +aj1p13i;
		t1a =aj1p1r -aj1p13r;	t1b =aj1p1i -aj1p13i;
		t16 =t16+t18;			t17 =t17+t19;
		t18 =t16+c3m1*t18;		t19 =t17+c3m1*t19;
		rt  =s3*t1a;			it  =s3*t1b;
		t1a =t18+it;			t1b =t19-rt;
		t18 =t18-it;			t19 =t19+rt;
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif
		t1c =aj1p5r;			t1d =aj1p5i;
		t1e =aj1p17r+aj1p11r;	t1f =aj1p17i+aj1p11i;
		t1g =aj1p17r-aj1p11r;	t1h =aj1p17i-aj1p11i;
		t1c =t1c+t1e;			t1d =t1d+t1f;
		t1e =t1c+c3m1*t1e;		t1f =t1d+c3m1*t1f;
		rt  =s3*t1g;			it  =s3*t1h;
		t1g =t1e+it;			t1h =t1f-rt;
		t1e =t1e-it;			t1f =t1f+rt;
#if PFETCH
addr = add0+p9;
prefetch_p_doubles(addr);
#endif	/* Twiddles: */
		rt  =t16;				it  =t17;
		t16 =rt+t1c;			t17 =it+t1d;
		t1c =rt-t1c;			t1d =it-t1d;
		t10 =t10+t16;			t11 =t11+t17;
		t16 =t10+c3m1*t16;		t17 =t11+c3m1*t17;
		rt  =s3*t1c;			it  =s3*t1d;
		t1c =t16+it;			t1d =t17-rt;
		t16 =t16-it;			t17 =t17+rt;
#if PFETCH
addr = add0+p10;
prefetch_p_doubles(addr);
#endif
		rt  =t18 *c -t19*s;		it  =t18 *s +t19*c;
		re  =t1e*c2-t1f*s2;		t1f =t1e*s2+t1f*c2;	t1e=re;
		t18 =rt+t1e;			t19 =it+t1f;
		t1e =rt-t1e;			t1f =it-t1f;
		t12 =t12+t18;			t13 =t13+t19;
		t18 =t12+c3m1*t18;		t19 =t13+c3m1*t19;
		rt  =s3*t1e;			it  =s3*t1f;
		t1e =t18+it;			t1f =t19-rt;
		t18 =t18-it;			t19 =t19+rt;
#if PFETCH
addr = add0+p11;
prefetch_p_doubles(addr);
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
addr = add0+p12;
prefetch_p_doubles(addr);
#endif
/*...and now do nine radix-2 transforms:	*/

		a[j1    ]=t00+t10;		a[j2    ]=t01+t11;
		a[j1+p1 ]=t00-t10;		a[j2+p1 ]=t01-t11;

		a[j1+p4 ]=t06+t16;		a[j2+p4 ]=t07+t17;
		a[j1+p5 ]=t06-t16;		a[j2+p5 ]=t07-t17;
#if PFETCH
addr = add0+p13;
prefetch_p_doubles(addr);
#endif
		a[j1+p3 ]=t0c+t1c;		a[j2+p3 ]=t0d+t1d;
		a[j1+p2 ]=t0c-t1c;		a[j2+p2 ]=t0d-t1d;

		a[j1+p16]=t02+t12;		a[j2+p16]=t03+t13;
		a[j1+p17]=t02-t12;		a[j2+p17]=t03-t13;
#if PFETCH
addr = add0+p14;
prefetch_p_doubles(addr);
#endif
		a[j1+p15]=t08+t18;		a[j2+p15]=t09+t19;
		a[j1+p14]=t08-t18;		a[j2+p14]=t09-t19;

		a[j1+p12]=t0e+t1e;		a[j2+p12]=t0f+t1f;
		a[j1+p13]=t0e-t1e;		a[j2+p13]=t0f-t1f;
#if PFETCH
addr = add0+p15;
prefetch_p_doubles(addr);
#endif
		a[j1+p11]=t04+t14;		a[j2+p11]=t05+t15;
		a[j1+p10]=t04-t14;		a[j2+p10]=t05-t15;

		a[j1+p8 ]=t0a+t1a;		a[j2+p8 ]=t0b+t1b;
		a[j1+p9 ]=t0a-t1a;		a[j2+p9 ]=t0b-t1b;
#if PFETCH
addr = add0+p16;
prefetch_p_doubles(addr);
#endif
		a[j1+p7 ]=t0g+t1g;		a[j2+p7 ]=t0h+t1h;
		a[j1+p6 ]=t0g-t1g;		a[j2+p6 ]=t0h-t1h;

#if PFETCH
addr = add0+p17;
prefetch_p_doubles(addr);
#endif
		}

		jstart += nwt;
		jhi    += nwt;
		col += 18;
		co3 -= 18;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-18 forward DIF FFT of the first block of 18 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 18 outputs of (1);
!   (3) Reweight and perform a radix-18 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 18 elements and repeat (1-4).
*/
	t00 = cy17;
	cy17= cy16;
	cy16= cy15;
	cy15= cy14;
	cy14= cy13;
	cy13= cy12;
	cy12= cy11;
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
	cy0 = t00;
/*
     if(cy0 != 0)printf("cy18!= 0\n");
else if(cy17!= 0)printf("cy17!= 0\n");
else if(cy16!= 0)printf("cy16!= 0\n");
else if(cy15!= 0)printf("cy15!= 0\n");
else if(cy14!= 0)printf("cy14!= 0\n");
else if(cy13!= 0)printf("cy13!= 0\n");
else if(cy12!= 0)printf("cy12!= 0\n");
else if(cy11!= 0)printf("cy11!= 0\n");
else if(cy10!= 0)printf("cy10!= 0\n");
else if(cy9 != 0)printf("cy9 != 0\n");
else if(cy8 != 0)printf("cy8 != 0\n");
else if(cy7 != 0)printf("cy7 != 0\n");
else if(cy6 != 0)printf("cy6 != 0\n");
else if(cy5 != 0)printf("cy5 != 0\n");
else if(cy4 != 0)printf("cy4 != 0\n");
else if(cy3 != 0)printf("cy3 != 0\n");
else if(cy2 != 0)printf("cy2 != 0\n");
else if(cy1 != 0)printf("cy1 != 0\n");
*/
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
		a[j+p11] *= radix_inv;
		a[j+p12] *= radix_inv;
		a[j+p13] *= radix_inv;
		a[j+p14] *= radix_inv;
		a[j+p15] *= radix_inv;
		a[j+p16] *= radix_inv;
		a[j+p17] *= radix_inv;
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9)
		+fabs(cy10)+fabs(cy11)+fabs(cy12)+fabs(cy13)+fabs(cy14)+fabs(cy15)+fabs(cy16)+fabs(cy17) != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix18_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix18_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-18 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,j1,j2;
	static int n18,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17, first_entry=TRUE;
	static double c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double rt,it,re
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h;

	if(!first_entry && (n/18) != n18)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n18=n/18;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n18;
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
		p16= p15+p1;
		p17= p16+p1;

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
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
		p15= p15+ ( (p15>> DAT_BITS) << PAD_BITS );
		p16= p16+ ( (p16>> DAT_BITS) << PAD_BITS );
		p17= p17+ ( (p17>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-18 pass is here.	*/

	for(j=0; j < n18; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
/*
Twiddleless version requires us to swap inputs as follows:
indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17
      -> 0,-9,16, 7,14, 5,12, 3,10, 1, 8,-1, 6,-3, 4,-5, 2,-7
      == 0, 9,16, 7,14, 5,12, 3,10, 1, 8,17, 6,15, 4,13, 2,11 modulo 18.
I.e. start out with first nonet of indices {0,2,4,6,8,10,12,14,16}, permute those according to
{0,2,4,6,8,10,12,14,16}*17%18 = {0,16,14,12,10,8,6,4,2}, then each is head of a length-2 list of indices with decrement 9.
*/

/*...First radix-9 transform:	*/
		t00 =a[j1    ];				t01 =a[j2    ];
		t02 =a[j1+p12]+a[j1+p6 ];	t03 =a[j2+p12]+a[j2+p6 ];
		t04 =a[j1+p12]-a[j1+p6 ];	t05 =a[j2+p12]-a[j2+p6 ];
		t00 =t00+t02;				t01 =t01+t03;
		t02 =t00+c3m1*t02;			t03 =t01+c3m1*t03;
		rt  =s3*t04;				it  =s3*t05;
		t04 =t02+it;				t05 =t03-rt;
		t02 =t02-it;				t03 =t03+rt;

		t06 =a[j1+p16];				t07 =a[j2+p16];
		t08 =a[j1+p10]+a[j1+p4 ];	t09 =a[j2+p10]+a[j2+p4 ];
		t0a =a[j1+p10]-a[j1+p4 ];	t0b =a[j2+p10]-a[j2+p4 ];
		t06 =t06+t08;				t07 =t07+t09;
		t08 =t06+c3m1*t08;			t09 =t07+c3m1*t09;
		rt  =s3*t0a;				it  =s3*t0b;
		t0a =t08+it;				t0b =t09-rt;
		t08 =t08-it;				t09 =t09+rt;

		t0c =a[j1+p14];				t0d =a[j2+p14];
		t0e =a[j1+p8 ]+a[j1+p2 ];	t0f =a[j2+p8 ]+a[j2+p2 ];
		t0g =a[j1+p8 ]-a[j1+p2 ];	t0h =a[j2+p8 ]-a[j2+p2 ];
		t0c =t0c+t0e;				t0d =t0d+t0f;
		t0e =t0c+c3m1*t0e;			t0f =t0d+c3m1*t0f;
		rt  =s3*t0g;				it  =s3*t0h;
		t0g =t0e+it;				t0h =t0f-rt;
		t0e =t0e-it;				t0f =t0f+rt;
	/* Twiddles: */
		rt  =t06;					it  =t07;
		t06 =rt+t0c;				t07 =it+t0d;
		t0c =rt-t0c;				t0d =it-t0d;
		t00 =t00+t06;				t01 =t01+t07;
		t06 =t00+c3m1*t06;			t07 =t01+c3m1*t07;
		rt  =s3*t0c;				it  =s3*t0d;
		t0c =t06+it;				t0d =t07-rt;
		t06 =t06-it;				t07 =t07+rt;

		rt  =t08*c -t09*s;			it  =t08*s +t09*c;
		re  =t0e*c2-t0f*s2;			t0f =t0e*s2+t0f*c2;	t0e=re;
		t08 =rt+t0e;				t09 =it+t0f;
		t0e =rt-t0e;				t0f =it-t0f;
		t02 =t02+t08;				t03 =t03+t09;
		t08 =t02+c3m1*t08;			t09 =t03+c3m1*t09;
		rt  =s3*t0e;				it  =s3*t0f;
		t0e =t08+it;				t0f =t09-rt;
		t08 =t08-it;				t09 =t09+rt;

		rt  =t0a*c2-t0b*s2;			it  =t0a*s2+t0b*c2;
		re  =t0g*c4-t0h*s4;			t0h =t0g*s4+t0h*c4;	t0g=re;
		t0a =rt+t0g;				t0b =it+t0h;
		t0g =rt-t0g;				t0h =it-t0h;
		t04 =t04+t0a;				t05 =t05+t0b;
		t0a =t04+c3m1*t0a;			t0b =t05+c3m1*t0b;
		rt  =s3*t0g;				it  =s3*t0h;
		t0g =t0a+it;				t0h =t0b-rt;
		t0a =t0a-it;				t0b =t0b+rt;

/*...Second radix-9 transform:	*/
		t10 =a[j1+p9 ];				t11 =a[j2+p9 ];
		t12 =a[j1+p3 ]+a[j1+p15];	t13 =a[j2+p3 ]+a[j2+p15];
		t14 =a[j1+p3 ]-a[j1+p15];	t15 =a[j2+p3 ]-a[j2+p15];
		t10 =t10+t12;				t11 =t11+t13;
		t12 =t10+c3m1*t12;			t13 =t11+c3m1*t13;
		rt  =s3*t14;				it  =s3*t15;
		t14 =t12+it;				t15 =t13-rt;
		t12 =t12-it;				t13 =t13+rt;

		t16 =a[j1+p7 ];				t17 =a[j2+p7 ];
		t18 =a[j1+p1 ]+a[j1+p13];	t19 =a[j2+p1 ]+a[j2+p13];
		t1a =a[j1+p1 ]-a[j1+p13];	t1b =a[j2+p1 ]-a[j2+p13];
		t16 =t16+t18;				t17 =t17+t19;
		t18 =t16+c3m1*t18;			t19 =t17+c3m1*t19;
		rt  =s3*t1a;				it  =s3*t1b;
		t1a =t18+it;				t1b =t19-rt;
		t18 =t18-it;				t19 =t19+rt;

		t1c =a[j1+p5 ];				t1d =a[j2+p5 ];
		t1e =a[j1+p17]+a[j1+p11];	t1f =a[j2+p17]+a[j2+p11];
		t1g =a[j1+p17]-a[j1+p11];	t1h =a[j2+p17]-a[j2+p11];
		t1c =t1c+t1e;				t1d =t1d+t1f;
		t1e =t1c+c3m1*t1e;			t1f =t1d+c3m1*t1f;
		rt  =s3*t1g;				it  =s3*t1h;
		t1g =t1e+it;				t1h =t1f-rt;
		t1e =t1e-it;				t1f =t1f+rt;
	/* Twiddles: */
		rt  =t16;					it  =t17;
		t16 =rt+t1c;				t17 =it+t1d;
		t1c =rt-t1c;				t1d =it-t1d;
		t10 =t10+t16;				t11 =t11+t17;
		t16 =t10+c3m1*t16;			t17 =t11+c3m1*t17;
		rt  =s3*t1c;				it  =s3*t1d;
		t1c =t16+it;				t1d =t17-rt;
		t16 =t16-it;				t17 =t17+rt;

		rt  =t18 *c -t19*s;			it  =t18 *s +t19*c;
		re  =t1e*c2-t1f*s2;			t1f =t1e*s2+t1f*c2;	t1e=re;
		t18 =rt+t1e;				t19 =it+t1f;
		t1e =rt-t1e;				t1f =it-t1f;
		t12 =t12+t18;				t13 =t13+t19;
		t18 =t12+c3m1*t18;			t19 =t13+c3m1*t19;
		rt  =s3*t1e;				it  =s3*t1f;
		t1e =t18+it;				t1f =t19-rt;
		t18 =t18-it;				t19 =t19+rt;

		rt  =t1a*c2-t1b*s2;			it  =t1a*s2+t1b*c2;
		re  =t1g*c4-t1h*s4;			t1h =t1g*s4+t1h*c4;	t1g=re;
		t1a =rt+t1g;				t1b =it+t1h;
		t1g =rt-t1g;				t1h =it-t1h;
		t14 =t14+t1a;				t15 =t15+t1b;
		t1a =t14+c3m1*t1a;			t1b =t15+c3m1*t1b;
		rt  =s3*t1g;				it  =s3*t1h;
		t1g =t1a+it;				t1h =t1b-rt;
		t1a =t1a-it;				t1b =t1b+rt;

/*...and now do nine radix-2 transforms:	*/

		a[j1    ]=t00+t10;			a[j2    ]=t01+t11;
		a[j1+p1 ]=t00-t10;			a[j2+p1 ]=t01-t11;

		a[j1+p4 ]=t06+t16;			a[j2+p4 ]=t07+t17;
		a[j1+p5 ]=t06-t16;			a[j2+p5 ]=t07-t17;

		a[j1+p3 ]=t0c+t1c;			a[j2+p3 ]=t0d+t1d;
		a[j1+p2 ]=t0c-t1c;			a[j2+p2 ]=t0d-t1d;

		a[j1+p16]=t02+t12;			a[j2+p16]=t03+t13;
		a[j1+p17]=t02-t12;			a[j2+p17]=t03-t13;

		a[j1+p15]=t08+t18;			a[j2+p15]=t09+t19;
		a[j1+p14]=t08-t18;			a[j2+p14]=t09-t19;

		a[j1+p12]=t0e+t1e;			a[j2+p12]=t0f+t1f;
		a[j1+p13]=t0e-t1e;			a[j2+p13]=t0f-t1f;

		a[j1+p11]=t04+t14;			a[j2+p11]=t05+t15;
		a[j1+p10]=t04-t14;			a[j2+p10]=t05-t15;

		a[j1+p8 ]=t0a+t1a;			a[j2+p8 ]=t0b+t1b;
		a[j1+p9 ]=t0a-t1a;			a[j2+p9 ]=t0b-t1b;

		a[j1+p7 ]=t0g+t1g;			a[j2+p7 ]=t0h+t1h;
		a[j1+p6 ]=t0g-t1g;			a[j2+p6 ]=t0h-t1h;
								/* Totals: 196 FADD, 80 FMUL.	*/
	}
}

/***************/

void radix18_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-18 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix18_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int j,j1,j2;
	static int n18,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17, first_entry=TRUE;
	static double c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double rt,it,re
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h;

	if(!first_entry && (n/18) != n18)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n18=n/18;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n18;
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
		p16= p15+p1;
		p17= p16+p1;

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
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
		p15= p15+ ( (p15>> DAT_BITS) << PAD_BITS );
		p16= p16+ ( (p16>> DAT_BITS) << PAD_BITS );
		p17= p17+ ( (p17>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-18 pass is here.	*/

	for(j=0; j < n18; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (18 64-bit complex, i.e. 36 64-bit reals) and do nine radix-2 transforms.	*/

	/*
	Twiddleless version requires us to swap inputs as follows:
	indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17
		  -> 0,-2,-4,-6,-8,-A,-C,-D,-G, 9, 7, 5, 3, 1,-1,-3,-5,-7
		  == 0,16,14,12,10, 8, 6, 4, 2, 9, 7, 5, 3, 1,17,15,13,11 modulo 18.
	I.e. start out with first pair of indices {0,9}, permute those according to
	{0,9}*17%18 = {0,9}, then each is head of a length-9 list of indices with decrement 2.

	Remember, inputs to DIT are bit-reversed, so
	a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17] contain
	x[0, 9, 3,12, 6,15, 1,10, 4,13, 7,16, 2,11, 5,14, 8,17], which get swapped to
	x[0, 9,12, 3, 6,15,16, 7,10, 1, 4,13,14, 5, 8,17, 2,11], which means the a-indices get swapped as
	a[0, 1, 3, 2, 4, 5,11,10, 7, 6, 8, 9,15,14,16,17,12,13].
	*/
		t00=a[j1    ];			t01=a[j2    ];
		rt =a[j1+p1 ];			it =a[j2+p1 ];
		t10=t00-rt;				t11=t01-it;
		t00=t00+rt;				t01=t01+it;

		t02=a[j1+p3 ];			t03=a[j2+p3 ];
		rt =a[j1+p2 ];			it =a[j2+p2 ];
		t12=t02-rt;				t13=t03-it; ;
		t02=t02+rt;				t03=t03+it; ;

		t04=a[j1+p4 ];			t05=a[j2+p4 ];
		rt =a[j1+p5 ];			it =a[j2+p5 ];
		t14=t04-rt;				t15=t05-it; ;
		t04=t04+rt;				t05=t05+it; ;

		t06=a[j1+p11];			t07=a[j2+p11];
		rt =a[j1+p10];			it =a[j2+p10];
		t16=t06-rt;				t17=t07-it; ;
		t06=t06+rt;				t07=t07+it; ;

		t08=a[j1+p7 ];			t09=a[j2+p7 ];
		rt =a[j1+p6 ];			it =a[j2+p6 ];
		t18=t08-rt;				t19=t09-it; ;
		t08=t08+rt;				t09=t09+it; ;

		t0a=a[j1+p8 ];			t0b=a[j2+p8 ];
		rt =a[j1+p9 ];			it =a[j2+p9 ];
		t1a=t0a-rt;				t1b=t0b-it; ;
		t0a=t0a+rt;				t0b=t0b+it; ;

		t0c=a[j1+p15];			t0d=a[j2+p15];
		rt =a[j1+p14];			it =a[j2+p14];
		t1c=t0c-rt;				t1d=t0d-it; ;
		t0c=t0c+rt;				t0d=t0d+it; ;

		t0e=a[j1+p16];			t0f=a[j2+p16];
		rt =a[j1+p17];			it =a[j2+p17];
		t1e=t0e-rt;				t1f=t0f-it; ;
		t0e=t0e+rt;				t0f=t0f+it; ;

		t0g=a[j1+p12];			t0h=a[j2+p12];
		rt =a[j1+p13];			it =a[j2+p13];
		t1g=t0g-rt;				t1h=t0h-it;
		t0g=t0g+rt;				t0h=t0h+it;

/*       ...and now do two radix-9 transforms.	*/
/*...t0[0:16:2]r use t[00:0g:2]; t0[1:17:2]i use t[10:1g:2]	*/

/*...First radix-9 transform:	*/
		rt  =t02;				it  =t03;
		t02 =rt+t04;			t03 =it+t05;
		t04 =rt-t04;			t05 =it-t05;
		t00 =t00+t02;			t01 =t01+t03;
		t02 =t00+c3m1*t02;		t03 =t01+c3m1*t03;
		rt  =s3*t04;			it  =s3*t05;
		t04 =t02-it;			t05 =t03+rt;
		t02 =t02+it;			t03 =t03-rt;

		rt  =t08;				it  =t09;
		t08 =rt+t0a;			t09 =it+t0b;
		t0a =rt-t0a;			t0b =it-t0b;
		t06 =t06+t08;			t07 =t07+t09;
		t08 =t06+c3m1*t08;		t09 =t07+c3m1*t09;
		rt  =s3*t0a;			it  =s3*t0b;
		t0a =t08-it;			t0b =t09+rt;
		t08 =t08+it;			t09 =t09-rt;

		rt  =t0e;				it  =t0f;
		t0e =rt+t0g;			t0f =it+t0h;
		t0g =rt-t0g;			t0h =it-t0h;
		t0c =t0c+t0e;			t0d =t0d+t0f;
		t0e =t0c+c3m1*t0e;		t0f =t0d+c3m1*t0f;
		rt  =s3*t0g;			it  =s3*t0h;
		t0g =t0e-it;			t0h =t0f+rt;
		t0e =t0e+it;			t0f =t0f-rt;
	/* Twiddles: */
		rt  =t06;				it  =t07;
		t06 =rt+t0c;			t07 =it+t0d;
		t0c =rt-t0c;			t0d =it-t0d;
		t00 =t00+t06;			t01 =t01+t07;
		a[j1    ]=t00;			a[j2    ]=t01;
		t06 =t00+c3m1*t06;		t07 =t01+c3m1*t07;
		rt  =s3*t0c;			it  =s3*t0d;
		a[j1+p6 ]=t06+it;		a[j2+p6 ]=t07-rt;
		a[j1+p12]=t06-it;		a[j2+p12]=t07+rt;

		rt  =t08*c +t09*s;		it  =t09*c -t08*s;
		re  =t0e*c2+t0f*s2;		t0f =t0f*c2-t0e*s2;	t0e=re;
		t08 =rt+t0e;			t09 =it+t0f;
		t0e =rt-t0e;			t0f =it-t0f;
		t02 =t02+t08;			t03 =t03+t09;
		a[j1+p8 ]=t02;			a[j2+p8 ]=t03;
		t08 =t02+c3m1*t08;		t09=t03+c3m1*t09;
		rt  =s3*t0e;			it  =s3*t0f;
		a[j1+p14]=t08+it;		a[j2+p14]=t09-rt;
		a[j1+p2 ]=t08-it;		a[j2+p2 ]=t09+rt;

		rt  =t0a*c2+t0b*s2;		it  =t0b*c2-t0a*s2;
		re  =t0g*c4+t0h*s4;		t0h =t0h*c4-t0g*s4;	t0g=re;
		t0a =rt+t0g;			t0b =it+t0h;
		t0g =rt-t0g;			t0h =it-t0h;
		t04 =t04+t0a;			t05 =t05+t0b;
		a[j1+p16]=t04;			a[j2+p16]=t05;
		t0a =t04+c3m1*t0a;		t0b =t05+c3m1*t0b;
		rt  =s3*t0g;			it  =s3*t0h;
		a[j1+p4 ]=t0a+it;		a[j2+p4 ]=t0b-rt;
		a[j1+p10]=t0a-it;		a[j2+p10]=t0b+rt;

/*...Second radix-9 transform:	*/
		rt  =t12;				it  =t13;
		t12 =rt+t14;			t13 =it+t15;
		t14 =rt-t14;			t15 =it-t15;
		t10 =t10+t12;			t11 =t11+t13;
		t12 =t10+c3m1*t12;		t13 =t11+c3m1*t13;
		rt  =s3*t14;			it  =s3*t15;
		t14 =t12-it;			t15 =t13+rt;
		t12 =t12+it;			t13 =t13-rt;

		rt  =t18;				it  =t19;
		t18 =rt+t1a;			t19 =it+t1b;
		t1a =rt-t1a;			t1b =it-t1b;
		t16 =t16+t18;			t17 =t17+t19;
		t18 =t16+c3m1*t18;		t19 =t17+c3m1*t19;
		rt  =s3*t1a;			it  =s3*t1b;
		t1a =t18-it;			t1b =t19+rt;
		t18 =t18+it;			t19 =t19-rt;

		rt  =t1e;				it  =t1f;
		t1e =rt+t1g;			t1f =it+t1h;
		t1g =rt-t1g;			t1h =it-t1h;
		t1c =t1c+t1e;			t1d =t1d+t1f;
		t1e =t1c+c3m1*t1e;		t1f =t1d+c3m1*t1f;
		rt  =s3*t1g;			it  =s3*t1h;
		t1g =t1e-it;			t1h =t1f+rt;
		t1e =t1e+it;			t1f =t1f-rt;
	/* Twiddles: */
		rt  =t16;				it  =t17;
		t16 =rt+t1c;			t17 =it+t1d;
		t1c =rt-t1c;			t1d =it-t1d;
		t10 =t10+t16;			t11 =t11+t17;
		a[j1+p9 ]=t10;			a[j2+p9 ]=t11;
		t16 =t10+c3m1*t16;		t17 =t11+c3m1*t17;
		rt  =s3*t1c;			it  =s3*t1d;
		a[j1+p15]=t16+it;		a[j2+p15]=t17-rt;
		a[j1+p3 ]=t16-it;		a[j2+p3 ]=t17+rt;

		rt  =t18*c +t19*s;		it  =t19*c -t18 *s;
		re  =t1e*c2+t1f*s2;		t1f =t1f*c2-t1e*s2;	t1e=re;
		t18 =rt+t1e;			t19 =it+t1f;
		t1e =rt-t1e;			t1f =it-t1f;
		t12 =t12+t18;			t13 =t13+t19;
		a[j1+p17]=t12;			a[j2+p17]=t13;
		t18 =t12+c3m1*t18;		t19=t13+c3m1*t19;
		rt  =s3*t1e;			it  =s3*t1f;
		a[j1+p5 ]=t18+it;		a[j2+p5 ]=t19-rt;
		a[j1+p11]=t18-it;		a[j2+p11]=t19+rt;

		rt  =t1a*c2+t1b*s2;		it  =t1b*c2-t1a*s2;
		re  =t1g*c4+t1h*s4;		t1h =t1h*c4-t1g*s4;	t1g=re;
		t1a =rt+t1g;			t1b =it+t1h;
		t1g =rt-t1g;			t1h =it-t1h;
		t14 =t14+t1a;			t15 =t15+t1b;
		a[j1+p7 ]=t14;			a[j2+p7 ]=t15;
		t1a =t14+c3m1*t1a;		t1b =t15+c3m1*t1b;
		rt  =s3*t1g;			it  =s3*t1h;
		a[j1+p13]=t1a+it;		a[j2+p13]=t1b-rt;
		a[j1+p1 ]=t1a-it;		a[j2+p1 ]=t1b+rt;

	}
}

