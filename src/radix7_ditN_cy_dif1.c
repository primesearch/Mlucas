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

int radix7_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[],double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-7 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-7 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n7,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6
		,i,j,j1,j2,jstart,jhi,root_incr,k1,k2,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
#if LO_ADD
	static double cc1 =  0.62348980185873353053,	/* cos(2*pi/7) */
		      ss1 =  0.78183148246802980870,	/* sin(2*pi/7) */
		      cc2 = -0.22252093395631440427,	/* cos(2*u) */
		      ss2 =  0.97492791218182360702,	/* sin(2*u) */
		      cc3 = -0.90096886790241912622,	/* cos(3*u) */
		      ss3 =  0.43388373911755812050;	/* sin(3*u) */
	double re,im;
#else
	static double cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3	*/
		      cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
		      cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
		      cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
		/* Switch the sign of ss3 in these: */
		      sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
		      sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
		      sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
		      sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#endif
	static double radix_inv, n2inv;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i
		,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6
		,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6
		,temp,frac,scale;
	double maxerr = 0.0;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double wt_re,wt_im, wi_re,wi_im;					/* Fermat/LOACC weights stuff */
	int ii0,ii1,ii2,ii3,ii4,ii5,ii6;					/* indices into weights arrays (mod NWT) */

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii0=ii1=ii2=ii3=ii4=ii5=ii6=0;

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

/*...change n7 and n_div_wt to non-static to work around a gcc compiler bug. */
	n7   = n/7;
	n_div_nwt = n7 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n7)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/7 in radix7_ditN_cy_dif1.\n",iter);
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
		ASSERT(LO_ADD,"radix7_ditN_cy_dif1.c: LO_ADD");
		psave = p;	nsave = n;
		first_entry=FALSE;

		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)7));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		/*   constant index offsets for array load/stores are here.	*/

		p1 = n7;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			bjmodnini=0;
			for(j=0; j < n7; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
		else
		{
			bjmodnini=0;
			for(j=0; j < n7/2; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
	}

/*...The radix-7 final DIT pass is here.	*/

	/* init carries	*/
	cy_r0 = 0;	cy_i0 = 0;
	cy_r1 = 0;	cy_i1 = 0;
	cy_r2 = 0;	cy_i2 = 0;
	cy_r3 = 0;	cy_i3 = 0;
	cy_r4 = 0;	cy_i4 = 0;
	cy_r5 = 0;	cy_i5 = 0;
	cy_r6 = 0;	cy_i6 = 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy_r0 = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/

	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/

	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart = 0;
		jhi = jstart+nwt-1;
		khi = n_div_nwt;
	}
	else
	{
		jstart = 0;
		jhi = n_div_nwt;
		khi = 1;
	}

for(outer=0; outer <= 1; outer++)
{
	i = 0;		/* Index into the BASE and BASEINV arrays. */
	/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
	if(bw > 0)
		i = 1;

	bjmodn0 = 0;
	bjmodn1 = bjmodnini;
	bjmodn2 = bjmodn1 +bjmodnini-n; bjmodn2 = bjmodn2 + ( (-(int)((uint32)bjmodn2 >> 31)) & n);
	bjmodn3 = bjmodn2 +bjmodnini-n; bjmodn3 = bjmodn3 + ( (-(int)((uint32)bjmodn3 >> 31)) & n);
	bjmodn4 = bjmodn3 +bjmodnini-n; bjmodn4 = bjmodn4 + ( (-(int)((uint32)bjmodn4 >> 31)) & n);
	bjmodn5 = bjmodn4 +bjmodnini-n; bjmodn5 = bjmodn5 + ( (-(int)((uint32)bjmodn5 >> 31)) & n);
	bjmodn6 = bjmodn5 +bjmodnini-n; bjmodn6 = bjmodn6 + ( (-(int)((uint32)bjmodn6 >> 31)) & n);

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		col=0;
		co2=(n >> nwt_bits)-1+7;
		co3=co2-7;		/* At the start of each new j-loop, co3=co2-radix(1)	*/
	}
	else
	{
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii0 = 0;
		ii1 = (SW_DIV_N*n7/2) % nwt;
		ii2 = (ii1 + ii1) % nwt;
		ii3 = (ii2 + ii1) % nwt;
		ii4 = (ii3 + ii1) % nwt;
		ii5 = (ii4 + ii1) % nwt;
		ii6 = (ii5 + ii1) % nwt;

		/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
		fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
		*/
		bjmodn0 = n;
	}

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

#if 0
	// We have just a single radix-7 DFT macro, which matches the DIF - for DIT must flip order of the last 6 outputs:
	RADIX_07_DFT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p0r,a1p0i,a1p6r,a1p6i,a1p5r,a1p5i,a1p4r,a1p4i,a1p3r,a1p3i,a1p2r,a1p2i,a1p1r,a1p1i,cc1,ss1,cc2,ss2,cc3,ss3,rt,it,re,im)

//	RADIX_07_DFT_NUSS(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p0r,a1p0i,a1p6r,a1p6i,a1p5r,a1p5i,a1p4r,a1p4i,a1p3r,a1p3i,a1p2r,a1p2i,a1p1r,a1p1i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it)

#else

	    t1 =a[j1   ];	t2 =a[j2   ];

	    t3 =a[j1+p1];	t4 =a[j2+p1];
	    rt =a[j1+p6];	it =a[j2+p6];
	    t13=t3 -rt;		t14=t4 -it;
	    t3 =t3 +rt;		t4 =t4 +it;

	    t5 =a[j1+p2];	t6 =a[j2+p2];
	    rt =a[j1+p5];	it =a[j2+p5];
	    t11=t5 -rt;		t12=t6 -it;
	    t5 =t5 +rt;		t6 =t6 +it;

	    t7 =a[j1+p3];	t8 =a[j2+p3];
	    rt =a[j1+p4];	it =a[j2+p4];
	    t9 =t7 -rt;		t10=t8 -it;
	    t7 =t7 +rt;		t8 =t8 +it;

	/* Version 1: direct computation of 3x3 matrix-vector product. Cost: 60 ADD, 36 MUL */
	#if LO_ADD
	    a1p0r = t1+t3+t5+t7;			a1p0i = t2+t4+t6+t8	;
	    rt = t1+cc1*t3+cc2*t5+cc3*t7;	it = t2+cc1*t4+cc2*t6+cc3*t8;
	    re = t1+cc2*t3+cc3*t5+cc1*t7;	im = t2+cc2*t4+cc3*t6+cc1*t8;
	    t1 = t1+cc3*t3+cc1*t5+cc2*t7;	t2 = t2+cc3*t4+cc1*t6+cc2*t8;

	    t3 = ss1*t13+ss2*t11+ss3*t9;	t4 = ss1*t14+ss2*t12+ss3*t10;
	    t5 = ss2*t13-ss3*t11-ss1*t9;	t6 = ss2*t14-ss3*t12-ss1*t10;
	    t7 = ss3*t13-ss1*t11+ss2*t9;	t8 = ss3*t14-ss1*t12+ss2*t10;

	    a1p1r=rt+t4;					a1p1i=it-t3;
	    a1p2r=re+t6;					a1p2i=im-t5;
	    a1p3r=t1+t8;					a1p3i=t2-t7;
	    a1p4r=t1-t8;					a1p4i=t2+t7;
	    a1p5r=re-t6;					a1p5i=im+t5;
	    a1p6r=rt-t4;					a1p6i=it+t3;

	/* Version 2: get cyclic convolution of (cc1,cc2,cc3) with (t3,t7,t5) and with (t4,t8,t6), (note switch of
	last 2 t-elements in each case), add t1 and t2 to respective outputs. Cost: 72 ADD, 16 MUL */
	#else
		rt = t3+t7+t5; a1p0r=rt+t1;	it = t4+t8+t6; a1p0i=it+t2;
		t1 = rt*cx0+t1;					t2 = it*cx0+t2;
		t3 = t3-t5;						t4 = t4-t6;
		t5 = t7-t5;						t6 = t8-t6;
		t7 =(t3+t5)*cx3;				t8 =(t4+t6)*cx3;
		t3 = t3*cx1;					t4 = t4*cx1;
		t5 = t5*cx2;					t6 = t6*cx2;
		rt = t3-t7;						it = t4-t8;
		t5 = t5-t7;						t6 = t6-t8;

		t3 = t1-rt-t5;					t4 = t2-it-t6;
		t5 = t1+t5;						t6 = t2+t6;
		t1 = t1+rt;						t2 = t2+it;

		t7 =(t13-t9 +t11)*sx0;			t8 =(t14-t10+t12)*sx0;
		t13= t13-t11;					t14= t14-t12;
		t11= t9 +t11;					t12= t10+t12;
		t9 =(t11-t13)*sx3;				t10=(t12-t14)*sx3;
		t13= t13*sx1;					t14= t14*sx1;
		t11= t11*sx2;					t12= t12*sx2;
		t13= t9+t13;					t14= t10+t14;
		t11= t9-t11;					t12= t10-t12;

		t9 = t7-t13-t11;				t10= t8-t14-t12;
		t11= t7+t11;					t12= t8+t12;
		t7 = t7+t13;					t8 = t8+t14;

		a1p1r =t1+t8;					a1p1i =t2-t7;
		a1p2r =t3+t10;					a1p2i =t4-t9;
		a1p3r =t5-t12;					a1p3i =t6+t11;
		a1p4r =t5+t12;					a1p4i =t6-t11;
		a1p5r =t3-t10;					a1p5i =t4+t9;
		a1p6r =t1-t8;					a1p6i =t2+t7;
	#endif
#endif

/*...and combine those to complete the radix-7 transform and do the carries. Since the outputs would
    normally be getting dispatched to 7 separate blocks of the A-array, we need 7 separate carries.	*/

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
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
			cmplx_carry_norm_errcheck0(a1p0r,a1p0i,cy_r0,bjmodn0,0,prp_mult);
			cmplx_carry_norm_errcheck (a1p1r,a1p1i,cy_r1,bjmodn1,1,prp_mult);
			cmplx_carry_norm_errcheck (a1p2r,a1p2i,cy_r2,bjmodn2,2,prp_mult);
			cmplx_carry_norm_errcheck (a1p3r,a1p3i,cy_r3,bjmodn3,3,prp_mult);
			cmplx_carry_norm_errcheck (a1p4r,a1p4i,cy_r4,bjmodn4,4,prp_mult);
			cmplx_carry_norm_errcheck (a1p5r,a1p5i,cy_r5,bjmodn5,5,prp_mult);
			cmplx_carry_norm_errcheck (a1p6r,a1p6i,cy_r6,bjmodn6,6,prp_mult);

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		}
		else
		{
			fermat_carry_norm_errcheck (a1p0r,a1p0i,cy_r0,cy_i0,ii0,bjmodn0,0x0*n7,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck (a1p1r,a1p1i,cy_r1,cy_i1,ii1,bjmodn1,0x1*n7,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck (a1p2r,a1p2i,cy_r2,cy_i2,ii2,bjmodn2,0x2*n7,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck (a1p3r,a1p3i,cy_r3,cy_i3,ii3,bjmodn3,0x3*n7,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck (a1p4r,a1p4i,cy_r4,cy_i4,ii4,bjmodn4,0x4*n7,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck (a1p5r,a1p5i,cy_r5,cy_i5,ii5,bjmodn5,0x5*n7,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck (a1p6r,a1p6i,cy_r6,cy_i6,ii6,bjmodn6,0x6*n7,NRTM1,NRT_BITS,prp_mult);
		}

/*...The radix-7 DIF pass is here:	*/
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif
/*...gather the needed data (7 64-bit complex, i.e. 14 64-bit reals) and begin the transform...	*/
	  t1 =a1p0r;					t2 =a1p0i;
	  t3 =a1p1r+a1p6r;			t4 =a1p1i+a1p6i;
	  t5 =a1p2r+a1p5r;			t6 =a1p2i+a1p5i;
	  t7 =a1p3r+a1p4r;			t8 =a1p3i+a1p4i;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
	  t9 =a1p3r-a1p4r;			t10=a1p3i-a1p4i;
	  t11=a1p2r-a1p5r;			t12=a1p2i-a1p5i;
	  t13=a1p1r-a1p6r;			t14=a1p1i-a1p6i;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif

/* Version 1: */
#if LO_ADD
	  a[j1   ] = t1+t3+t5+t7;		a[j2   ] = t2+t4+t6+t8	;
	  rt = t1+cc1*t3+cc2*t5+cc3*t7;	it = t2+cc1*t4+cc2*t6+cc3*t8;
	  re = t1+cc2*t3+cc3*t5+cc1*t7;	im = t2+cc2*t4+cc3*t6+cc1*t8;
	  t1 = t1+cc3*t3+cc1*t5+cc2*t7;	t2 = t2+cc3*t4+cc1*t6+cc2*t8;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
	  t3 = ss1*t13+ss2*t11+ss3*t9;	t4 = ss1*t14+ss2*t12+ss3*t10;
	  t5 = ss2*t13-ss3*t11-ss1*t9;	t6 = ss2*t14-ss3*t12-ss1*t10;
	  t7 = ss3*t13-ss1*t11+ss2*t9;	t8 = ss3*t14-ss1*t12+ss2*t10;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
	  a[j1+p1]=rt-t4;				a[j2+p1]=it+t3;
	  a[j1+p2]=re-t6;				a[j2+p2]=im+t5;
	  a[j1+p3]=t1-t8;				a[j2+p3]=t2+t7;
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
	  a[j1+p4]=t1+t8;				a[j2+p4]=t2-t7;
	  a[j1+p5]=re+t6;				a[j2+p5]=im-t5;
	  a[j1+p6]=rt+t4;				a[j2+p6]=it-t3;
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif

/* Version 2: */
#else
	  rt = t3+t7+t5; a[j1]=rt+t1;	it = t4+t8+t6; a[j2]=it+t2;
	  t1 = rt*cx0+t1;				t2 = it*cx0+t2;
	  t3 = t3-t5;					t4 = t4-t6;
	  t5 = t7-t5;					t6 = t8-t6;
	  t7 =(t3+t5)*cx3;				t8 =(t4+t6)*cx3;
	  t3 = t3*cx1;					t4 = t4*cx1;
	  t5 = t5*cx2;					t6 = t6*cx2;
	  rt = t3-t7;					it = t4-t8;
	  t5 = t5-t7;					t6 = t6-t8;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
	  t3 = t1-rt-t5;				t4 = t2-it-t6;
	  t5 = t1+t5;					t6 = t2+t6;
	  t1 = t1+rt;					t2 = t2+it;

#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
	  t7 =(t13-t9 +t11)*sx0;		t8 =(t14-t10+t12)*sx0;
	  t13= t13-t11;					t14= t14-t12;
	  t11= t9 +t11;					t12= t10+t12;
	  t9 =(t11-t13)*sx3;			t10=(t12-t14)*sx3;
	  t13= t13*sx1;					t14= t14*sx1;
	  t11= t11*sx2;					t12= t12*sx2;
	  t13= t9+t13;					t14= t10+t14;
	  t11= t9-t11;					t12= t10-t12;
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
	  t9 = t7-t13-t11;				t10= t8-t14-t12;
	  t11= t7+t11;					t12= t8+t12;
	  t7 = t7+t13;					t8 = t8+t14;
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif
	  a[j1+p1] =t1-t8;				a[j1+p1+1] =t2+t7;
	  a[j1+p2] =t3-t10;				a[j1+p2+1] =t4+t9;
	  a[j1+p3] =t5+t12;				a[j1+p3+1] =t6-t11;
	  a[j1+p4] =t5-t12;				a[j1+p4+1] =t6+t11;
	  a[j1+p5] =t3+t10;				a[j1+p5+1] =t4-t9;
	  a[j1+p6] =t1+t8;				a[j1+p6+1] =t2-t7;
#endif	/* endif #LO_ADD */
	  }

	  if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	  {
	    jstart += nwt;
	    jhi    += nwt;

	    col += 7;
	    co3 -= 7;
	  }
	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-7 forward DIF FFT of the first block of 7 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 7 outputs of (1);
!   (3) Reweight and perform a radix-7 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 7 elements and repeat (1-4).
*/
	/*
printf("iter = %10d; maxerr = %20.15f\n",iter,maxerr);
printf("carries = %10d %10d %10d %10d %10d %10d %10d\n",(int)cy_r0,(int)cy_r1,(int)cy_r2,(int)cy_r3,(int)cy_r4,(int)cy_r5,(int)cy_r6);
*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1    = cy_r6;
		cy_r6 = cy_r5;
		cy_r5 = cy_r4;
		cy_r4 = cy_r3;
		cy_r3 = cy_r2;
		cy_r2 = cy_r1;
		cy_r1 = cy_r0;
		cy_r0 =    t1;
	}
	else
	{
		/* ...The 2 Mo"bius carries are here: */
		t1    = cy_r6;	t2    = cy_i6;
		cy_r6 = cy_r5;	cy_i6 = cy_i5;
		cy_r5 = cy_r4;	cy_i5 = cy_i4;
		cy_r4 = cy_r3;	cy_i4 = cy_i3;
		cy_r3 = cy_r2;	cy_i3 = cy_i2;
		cy_r2 = cy_r1;	cy_i2 = cy_i1;
		cy_r1 = cy_r0;	cy_i1 = cy_i0;
		cy_r0 =   -t2;	cy_i0 =   +t1;
	}

	root_incr = 0;
	scale = prp_mult = 1;

	jstart = 0;
	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		jhi =15;
	}
	else
	{
		jhi = 7;
	}
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
	}
}

	if(fabs(cy_r0)+fabs(cy_r1)+fabs(cy_r2)+fabs(cy_r3)+fabs(cy_r4)+fabs(cy_r5)+fabs(cy_r6)
	  +fabs(cy_i0)+fabs(cy_i1)+fabs(cy_i2)+fabs(cy_i3)+fabs(cy_i4)+fabs(cy_i5)+fabs(cy_i6) != 0.0)
	{
	    sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix7_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
	    mlucas_fprint(cbuf,INTERACT);
	    err=ERR_CARRY;
	    return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix7_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-7 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Given complex inputs (x0,x1,x2,x3,x4,x5,x6), we need the following outputs
!   (here cJ = cos(2*J*pi/7), sJ = sin(2*J*pi/7)):
!
!	X0 = x0+x1+x2+x3+x4+x5+x6
!	X1 = x0+c1*(x1+x6)+c2*(x2+x5)+c3*(x3+x4) + I*[ s1*(x1-x6)+s2*(x2-x5)+s3*(x3-x4)] = C1 + I*S1
!	X2 = x0+c2*(x1+x6)+c3*(x2+x5)+c1*(x3+x4) + I*[ s2*(x1-x6)-s3*(x2-x5)-s1*(x3-x4)] = C2 + I*S2
!	X3 = x0+c3*(x1+x6)+c1*(x2+x5)+c2*(x3+x4) + I*[ s3*(x1-x6)-s1*(x2-x5)+s2*(x3-x4)] = C3 + I*S3
!	X4 = x0+c3*(x1+x6)+c1*(x2+x5)+c2*(x3+x4) - I*[ s3*(x1-x6)+s1*(x2-x5)+s2*(x3-x4)] = C3 - I*S3
!	X5 = x0+c2*(x1+x6)+c3*(x2+x5)+c1*(x3+x4) - I*[ s2*(x1-x6)-s3*(x2-x5)-s1*(x3-x4)] = C2 - I*S2
!	X6 = x0+c1*(x1+x6)+c2*(x2+x5)+c3*(x3+x4) - I*[ s1*(x1-x6)-s2*(x2-x5)+s3*(x3-x4)] = C1 - I*S1

!    We refer to the terms C1,2,3 (which do not explicitly involving the imaginary constant I)
!    as the "cosine part" of the output, and S1,2,3 (those multiplied by I) as the "sine part."
!
!    Form      (x1+-x6), (x2+-x5), (x3+-x4)   :  0 FMUL, 12 FADD
!    Form X0                                  :  0 FMUL,  6 FADD
!    Form x0+c1*(x1+x6)+c2*(x2+x5)+c3*(x3+x4) :  6 FMUL,  6 FADD
!    Form x0+c2*(x1+x6)+c3*(x2+x5)+c1*(x3+x4) :  6 FMUL,  6 FADD
!    Form x0+c3*(x1+x6)+c1*(x2+x5)+c2*(x3+x4) :  6 FMUL,  6 FADD
!    Form    s1*(x1-x6)+s2*(x2-x5)+s3*(x3-x4) :  6 FMUL,  4 FADD
!    Form    s2*(x1-x6)-s3*(x2-x5)-s1*(x3-x4) :  6 FMUL,  4 FADD
!    Form    s3*(x1-x6)-s1*(x2-x5)+s2*(x3-x4) :  6 FMUL,  4 FADD
!    Form X1,2,3,4,5,6                        :  0 FMUL, 12 FADD
!
!    Totals :                                   36 FMUL, 60 FADD.
!
!   Alternatively, we can get the X0 and the 3 x0+cosine terms for the cost of 2 length-3 real cyclic convolutions plus 4 extra adds.
!   We can write the 3 sine terms as the outputs of another pair of length-3 real cyclic convolutions by changing the
!   signs of s3 and (x3-x4) and changing the sign of the third resulting convolution output. Thus, we replace the
!   middle 7 "Form..." steps above with 4 length-3 real cyclic convolutions plus 4 extra adds, which costs just 4*4 = 16 FMUL
!   and 4*11+4 = 48 FADD, for a total of   16 FMUL, 72 FADD.

!   The low-multiply scheme may actually require a few more clock cycles to complete
!   on hardware which can complete both an FADD and an FMUL every clock cycle, but will have better
!   floating-point accuracy and is much more desirable for a modular implementation, where multiplies are expensive.
*/
	int j,j1,j2;
	static int n7,p1,p2,p3,p4,p5,p6, first_entry=TRUE;
	static double cc1 =  0.62348980185873353053,	/* cos(2*pi/7) */
		      ss1 =  0.78183148246802980870,	/* sin(2*pi/7) */
		      cc2 = -0.22252093395631440427,	/* cos(2*u) */
		      ss2 =  0.97492791218182360702,	/* sin(2*u) */
		      cc3 = -0.90096886790241912622,	/* cos(3*u) */
		      ss3 =  0.43388373911755812050;	/* sin(3*u) */
	double rt,it,re,im
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14;

	if(!first_entry && (n/7) != n7)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n7=n/7;

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n7;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;
	  p6 = p5 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	  p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	  p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-7 pass is here.	*/

      for(j=0; j < n7; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*...gather the needed data (7 64-bit complex, i.e. 14 64-bit reals) and begin the transform...	*/

	  t1 =a[j1   ];					t2 =a[j2   ];			/* x0	*/
	  t3 =a[j1+p1]+a[j1+p6];		t4 =a[j2+p1]+a[j2+p6];	/* x1 + x6	*/
	  t5 =a[j1+p2]+a[j1+p5];		t6 =a[j2+p2]+a[j2+p5];	/* x2 + x5	*/
	  t7 =a[j1+p3]+a[j1+p4];		t8 =a[j2+p3]+a[j2+p4];	/* x3 + x4	*/
	  t9 =a[j1+p3]-a[j1+p4];		t10=a[j2+p3]-a[j2+p4];	/* x3 - x4	*/
	  t11=a[j1+p2]-a[j1+p5];		t12=a[j2+p2]-a[j2+p5];	/* x2 - x5	*/
	  t13=a[j1+p1]-a[j1+p6];		t14=a[j2+p1]-a[j2+p6];	/* x1 - x6	*/

/*      ...now complete the first complex datum of the result.	*/

	  a[j1   ] = t1+t3+t5+t7;		a[j2   ] = t2+t4+t6+t8;			/* X0	*/
	  rt = t1+cc1*t3+cc2*t5+cc3*t7;	it = t2+cc1*t4+cc2*t6+cc3*t8;	/* C1	*/
	  re = t1+cc2*t3+cc3*t5+cc1*t7;	im = t2+cc2*t4+cc3*t6+cc1*t8;	/* C2	*/
	  t1 = t1+cc3*t3+cc1*t5+cc2*t7;	t2 = t2+cc3*t4+cc1*t6+cc2*t8;	/* C3	*/

	  t3 = ss1*t13+ss2*t11+ss3*t9;	t4 = ss1*t14+ss2*t12+ss3*t10;	/* S1	*/
	  t5 = ss2*t13-ss3*t11-ss1*t9;	t6 = ss2*t14-ss3*t12-ss1*t10;	/* S2	*/
	  t7 = ss3*t13-ss1*t11+ss2*t9;	t8 = ss3*t14-ss1*t12+ss2*t10;	/* S3	*/

/*...Inline multiply of sine parts by +-I into finishing phase...	*/

	  a[j1+p1]=rt-t4;				a[j2+p1]=it+t3;	/* X1 = C1 + I*S1	*/
	  a[j1+p2]=re-t6;				a[j2+p2]=im+t5;	/* X2 = C2 + I*S2	*/
	  a[j1+p3]=t1-t8;				a[j2+p3]=t2+t7;	/* X3 = C3 + I*S3	*/
	  a[j1+p4]=t1+t8;				a[j2+p4]=t2-t7;	/* X4 =	C3 - I*S3	*/
	  a[j1+p5]=re+t6;				a[j2+p5]=im-t5;	/* X5 =	C2 - I*S2	*/
	  a[j1+p6]=rt+t4;				a[j2+p6]=it-t3;	/* X6 =	C1 - I*S1	*/
                                                                /* Totals: 60 FADD, 36 FMUL.	*/
	}
}

/***************/

void radix7_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-7 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix7_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
!
!   See the documentation in radix16_dif_pass for further details on the array indexing,
!   and radix7_dif_pass1 for details on the algorithm.
*/
	int j,j1,j2;
	static int n7,p1,p2,p3,p4,p5,p6, first_entry=TRUE;
	static double cc1 =  0.62348980185873353053,	/* cos(2*pi/7) */
		      ss1 =  0.78183148246802980870,	/* sin(2*pi/7) */
		      cc2 = -0.22252093395631440427,	/* cos(2*u) */
		      ss2 =  0.97492791218182360702,	/* sin(2*u) */
		      cc3 = -0.90096886790241912622,	/* cos(3*u) */
		      ss3 =  0.43388373911755812050;	/* sin(3*u) */
	double rt,it,re,im
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14;

	if(!first_entry && (n/7) != n7)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n7=n/7;

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n7;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;
	  p6 = p5 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	  p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	  p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-7 pass is here.	*/

      for(j=0; j < n7; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*...gather the needed data (7 64-bit complex, i.e. 14 64-bit reals) and begin the transform...	*/

	  t1 =a[j1   ];					t2 =a[j2   ]		;	/* x0	*/
	  t3 =a[j1+p1]+a[j1+p6];		t4 =a[j2+p1]+a[j2+p6];	/* x1 + x6	*/
	  t5 =a[j1+p2]+a[j1+p5];		t6 =a[j2+p2]+a[j2+p5];	/* x2 + x5	*/
	  t7 =a[j1+p3]+a[j1+p4];		t8 =a[j2+p3]+a[j2+p4];	/* x3 + x4	*/
	  t9 =a[j1+p3]-a[j1+p4];		t10=a[j2+p3]-a[j2+p4];	/* x3 - x4	*/
	  t11=a[j1+p2]-a[j1+p5];		t12=a[j2+p2]-a[j2+p5];	/* x2 - x5	*/
	  t13=a[j1+p1]-a[j1+p6];		t14=a[j2+p1]-a[j2+p6];	/* x1 - x6	*/

/*      ...now complete the first complex datum of the result.	*/

	  a[j1   ] = t1+t3+t5+t7;		a[j2   ] = t2+t4+t6+t8	;	/* X0	*/
	  rt = t1+cc1*t3+cc2*t5+cc3*t7;	it = t2+cc1*t4+cc2*t6+cc3*t8;	/* C1	*/
	  re = t1+cc2*t3+cc3*t5+cc1*t7;	im = t2+cc2*t4+cc3*t6+cc1*t8;	/* C2	*/
	  t1 = t1+cc3*t3+cc1*t5+cc2*t7;	t2 = t2+cc3*t4+cc1*t6+cc2*t8;	/* C3	*/

	  t3 = ss1*t13+ss2*t11+ss3*t9;	t4 = ss1*t14+ss2*t12+ss3*t10;	/* S1	*/
	  t5 = ss2*t13-ss3*t11-ss1*t9;	t6 = ss2*t14-ss3*t12-ss1*t10;	/* S2	*/
	  t7 = ss3*t13-ss1*t11+ss2*t9;	t8 = ss3*t14-ss1*t12+ss2*t10;	/* S3	*/

/*...Inline multiply of sine parts by +-I into finishing phase...	*/

	  a[j1+p1]=rt+t4;				a[j2+p1]=it-t3;	/* X1 = C1 - I*S1	*/
	  a[j1+p2]=re+t6;				a[j2+p2]=im-t5;	/* X2 = C2 - I*S2	*/
	  a[j1+p3]=t1+t8;				a[j2+p3]=t2-t7;	/* X3 = C3 - I*S3	*/
	  a[j1+p4]=t1-t8;				a[j2+p4]=t2+t7;	/* X4 =	C3 + I*S3	*/
	  a[j1+p5]=re-t6;				a[j2+p5]=im+t5;	/* X5 =	C2 + I*S2	*/
	  a[j1+p6]=rt-t4;				a[j2+p6]=it+t3;	/* X6 =	C1 + I*S1	*/
                                                                /* Totals: 60 FADD, 36 FMUL.	*/

	}
}

