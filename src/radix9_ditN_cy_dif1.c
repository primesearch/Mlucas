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

int radix9_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-9 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-9 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n9,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8
		,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
#if LO_ADD
	static double c   =  0.76604444311897803520,	/* cos(2*pi/9) */
		      s   =  0.64278760968653932631,	/* sin(2*pi/9) */
		      c2  =  0.17364817766693034887,	/* cos(2*u) */
		      s2  =  0.98480775301220805936,	/* sin(2*u) */
		      c3m1= -1.50000000000000000000,	/* cos(3*u)-1] */
		      s3  =  0.86602540378443864677,	/* sin(3*u)] */
		      c4  = -0.93969262078590838404,	/* cos(4*u) */
		      s4  =  0.34202014332566873307;	/* sin(4*u) */
#else
	static double cx1 = 1.70573706390488641924, 	/*  cc1-cc4		*/
		      cx2 = 1.11334079845283873291, 	/*  cc2-cc4		*/
		      cx3 = 0.93969262078590838405,	/* (cc1+cc2-2*cc4)/3	*/
		/* Switch the sign of ss2 in these: */
		      sx1 = 0.30076746636087059324, 	/*  ss1-ss4		*/
		      sx2 =-1.32682789633787679243, 	/* -ss2-ss4		*/
		      ss3 = 0.86602540378443864677,	/*  ss3			*/
		      sx3 =-0.34202014332566873306;	/* (ss1-ss2-2*ss4)/3	*/
	double im,r2,i2;
#endif
	static double radix_inv, n2inv;
	double rt,it,re
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,temp,frac,scale;
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

/*...change n9 and n_div_wt to non-static to work around a gcc compiler bug. */
	n9   = n/9;
	n_div_nwt = n9 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n9)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/9 in radix9_ditN_cy_dif1.\n",iter);
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
	  radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)9));
	  n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

	  bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
	  sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n9;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;
	  p6 = p5 +p1;
	  p7 = p6 +p1;
	  p8 = p7 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	  p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	  p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
	  p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
	  p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );

	  bjmodnini=0;
	  for(j=0; j < n9; j++)
	  {
	    bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
	  }
	}

/*...The radix-9 final DIT pass is here.	*/

	cy0 = cy1 = cy2 = cy3 = cy4 = cy5 = cy6 = cy7 = cy8 = 0;	/* init carries	*/

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

	col=0;
	co2=(n >> nwt_bits)-1+9;
	co3=co2-9;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

#if LO_ADD
/*
!...gather the needed data (9 64-bit complex, i.e. 18 64-bit reals) and do three radix-3 transforms...
!	including the external twiddles on the inputs...
*/
	  t1 =a[j1   ];			t2 =a[j2   ];
	  t3 =a[j1+p1];			t4 =a[j2+p1];
	  rt =a[j1+p2];			it =a[j2+p2];
	  t5 =t3 -rt;			t6 =t4 -it;
	  t3 =t3 +rt;			t4 =t4 +it;
	  t1 =t1+t3;			t2 =t2+t4;
	  t3 =t1+c3m1*t3;		t4 =t2+c3m1*t4;
	  rt =s3*t5;			it =s3*t6;
	  t5 =t3-it;			t6 =t4+rt;
	  t3 =t3+it;			t4 =t4-rt;

	  t7 =a[j1+p3];			t8 =a[j2+p3];
	  t9 =a[j1+p4];			t10=a[j2+p4];
	  rt =a[j1+p5];			it =a[j2+p5];
	  t11=t9 -rt;			t12=t10-it;
	  t9 =t9 +rt;			t10=t10+it;
	  t7 =t7+t9;			t8 =t8+t10;
	  t9 =t7+c3m1*t9;		t10=t8+c3m1*t10;
	  rt =s3*t11;			it =s3*t12;
	  t11=t9-it;			t12=t10+rt;
	  t9 =t9+it;			t10=t10-rt;

	  t13=a[j1+p6];			t14=a[j2+p6];
	  t15=a[j1+p7];			t16=a[j2+p7];
	  rt =a[j1+p8];			it =a[j2+p8];
	  t17=t15-rt;			t18=t16-it;
	  t15=t15+rt;			t16=t16+it;
	  t13=t13+t15;			t14=t14+t16;
	  t15=t13+c3m1*t15;		t16=t14+c3m1*t16;
	  rt =s3*t17;			it =s3*t18;
	  t17=t15-it;			t18=t16+rt;
	  t15=t15+it;			t16=t16-rt;
/*
!...and now do three more radix-3 transforms, including the twiddle factors:
!                                           1, exp(-i*2*pi/9), exp(-i*4*pi/9) (for inputs to transform block 2)
!                                           1, exp(-i*4*pi/9), exp(-i*8*pi/9) (for inputs to transform block 3).
!   I.e. do similar as above, except inputs a(j1  +p0:8:1) are replaced by t1:17:2,
!                                           a(j2+p0:8:1) are replaced by t2:18:2, and v.v. for outputs,
!   and only the last 2 inputs to radix-3 transforms 2 and 3 are multiplied by non-unity twiddles.
*/
	  rt =t7;				it =t8;
	  t7 =rt+t13;			t8 =it+t14;
	  t13=rt-t13;			t14=it-t14;
	  t1 =t1+t7;			t2 =t2+t8;
	  aj1p0r=t1;			aj1p0i=t2;
	  t7 =t1+c3m1*t7;		t8 =t2+c3m1*t8;
	  rt =s3*t13;			it =s3*t14;
	  aj1p3r=t7+it;			aj1p3i=t8-rt;
	  aj1p6r=t7-it;			aj1p6i=t8+rt;

	  rt =t9 *c +t10*s;		it =t10*c -t9 *s;		/* twiddle mpy by E	*/
	  re =t15*c2+t16*s2;	t16=t16*c2-t15*s2;	t15=re;	/* twiddle mpy by E^2	*/
	  t9 =rt+t15;			t10=it+t16;
	  t15=rt-t15;			t16=it-t16;
	  t3 =t3+t9;			t4 =t4+t10;
	  aj1p1r=t3;			aj1p1i=t4;
	  t9 =t3+c3m1*t9;		t10=t4+c3m1*t10;
	  rt =s3*t15;			it =s3*t16;
	  aj1p4r=t9+it;			aj1p4i=t10-rt;
	  aj1p7r=t9-it;			aj1p7i=t10+rt;

	  rt =t11*c2+t12*s2;	it =t12*c2-t11*s2;		/* twiddle mpy by E^2	*/
	  re =t17*c4+t18*s4;	t18=t18*c4-t17*s4;	t17=re;	/* twiddle mpy by E^4	*/
	  t11=rt+t17;			t12=it+t18;
	  t17=rt-t17;			t18=it-t18;
	  t5 =t5+t11;			t6 =t6+t12;
	  aj1p2r=t5;			aj1p2i=t6;
	  t11=t5+c3m1*t11;		t12=t6+c3m1*t12;
	  rt =s3*t17;			it =s3*t18;
	  aj1p5r=t11+it;		aj1p5i=t12-rt;
	  aj1p8r=t11-it;		aj1p8i=t12+rt;
#else
/*...gather the needed data and do a subconvolution-based transform.
     We expect an input order [0,3,6,1,4,7,2,5,8], so scramble array indices accordingly. */
	  t1 =a[j1   ];			t2 =a[j2   ];	/* x0		*/

	  t3 =a[j1+p3];			t4 =a[j2+p3];
	  rt =a[j1+p8];			it =a[j2+p8];
	  t17=t3-rt;			t18=t4-it;	/* x1 - x8	*/
	  t3 =t3+rt;			t4 =t4+it;	/* x1 + x8	*/

	  t5 =a[j1+p6];			t6 =a[j2+p6];
	  rt =a[j1+p5];			it =a[j2+p5];
	  t15=t5-rt;			t16=t6-it;	/* x2 - x7	*/
	  t5 =t5+rt;			t6 =t6+it;	/* x2 + x7	*/

	  t7 =a[j1+p1];			t8 =a[j2+p1];
	  rt =a[j1+p2];			it =a[j2+p2];
	  t13=t7-rt;			t14=t8-it;	/* x3 - x6	*/
	  t7 =t7+rt;			t8 =t8+it;	/* x3 + x6	*/

	  t9 =a[j1+p4];			t10=a[j2+p4];
	  rt =a[j1+p7];			it =a[j2+p7];
	  t11=t9-rt;			t12=t10-it;	/* x4 - x5	*/
	  t9 =t9+rt;			t10=t10+it;	/* x4 + x5	*/

	  rt = t1+t7;			it = t2+t8;
	  r2 = t3+t5+t9;		i2 = t4+t6+t10;
	  aj1p0r = rt+r2;		aj1p0i = it+i2;	/* X0	*/

	  re = t1-0.5*t7;		im = t2-0.5*t8;
	  t7 = rt-0.5*r2;		t8 = it-0.5*i2;			/* C3 calculated separately. */
	  t3 = t3-t5;			t4 = t4 -t6;
	  t5 = t9-t5;			t6 = t10-t6;
	  t9 =(t3+t5)*cx3;		t10=(t4+t6)*cx3;
	  t3 = t3*cx1;			t4 = t4*cx1;
	  t5 = t5*cx2;			t6 = t6*cx2;
	  rt = t3-t9;			it = t4-t10;
	  t5 = t5-t9;			t6 = t6-t10;

	  t3 = re-rt-t5;		t4 = im-it-t6;	/* C2 */
	  t5 = re+t5;			t6 = im+t6;	/* C4 */
	  t1 = re+rt;			t2 = im+it;	/* C1 */
	/* For sine convo, replace t3,t5,t9 --> t17,-t15,t11 on input */
	  re = t13*ss3;			im = t14*ss3;
	  t13= (t17-t15+t11)*ss3;t14= (t18-t16+t12)*ss3;		/* S3 calculated separately. */
	  t17= t17+t15;			t18= t18+t16;
	  t15= t11+t15;			t16= t12+t16;
	  t11=(t17+t15)*sx3;	t12=(t18+t16)*sx3;
	  t17= t17*sx1;			t18= t18*sx1;
	  t15= t15*sx2;			t16= t16*sx2;
	  rt = t17-t11;			it = t18-t12;
	  t15= t15-t11;			t16= t16-t12;

	  t17= re-rt-t15;		t18= im-it-t16;	/*-S2 */
	  t15= re+t15;			t16= im+t16;	/* S4 */
	  t9 = re+rt;			t10= im+it;	/* S1 */

/*...Inline multiply of sine parts by +-I into finishing phase.	Outputs are ordered. */

	  aj1p1r=t1+t10;		aj1p1i=t2-t9 ;	/* X1 = C1 - I*S1	*/
	  aj1p2r=t3-t18;		aj1p2i=t4+t17;	/* X2 = C2 - I*S2	*/
	  aj1p3r=t7+t14;		aj1p3i=t8-t13;	/* X3 = C3 - I*S3	*/
	  aj1p4r=t5+t16;		aj1p4i=t6-t15;	/* X4 = C4 - I*S4	*/
	  aj1p5r=t5-t16;		aj1p5i=t6+t15;	/* X5 =	C4 + I*S4	*/
	  aj1p6r=t7-t14;		aj1p6i=t8+t13;	/* X6 =	C3 + I*S3	*/
	  aj1p7r=t3+t18;		aj1p7i=t4-t17;	/* X7 =	C2 + I*S2	*/
	  aj1p8r=t1-t10;		aj1p8i=t2+t9 ;	/* X8 =	C1 + I*S1	*/
#endif

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 9 separate blocks of the A-array, we need 9 separate carries.	*/

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

	    i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
	    co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			   and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	/*...The radix-9 DIF pass is here:	*/
	#if PFETCH
	add0 = &a[j1];
	prefetch_p_doubles(add0);
	#endif

#if LO_ADD
	  t1 =aj1p0r;			t2 =aj1p0i;
	  t3 =aj1p3r+aj1p6r;	t4 =aj1p3i+aj1p6i;
	  t5 =aj1p3r-aj1p6r;	t6 =aj1p3i-aj1p6i;
	  t1 =t1+t3;			t2 =t2+t4;
	  t3 =t1+c3m1*t3;		t4 =t2+c3m1*t4;
	  rt =s3*t5;			it =s3*t6;
	  t5 =t3+it;			t6 =t4-rt;
	  t3 =t3-it;			t4 =t4+rt;
	#if PFETCH
	addr = add0+p1;
	prefetch_p_doubles(addr);
	#endif
	  t7 =aj1p1r;			t8 =aj1p1i;
	  t9 =aj1p4r+aj1p7r;	t10=aj1p4i+aj1p7i;
	  t11=aj1p4r-aj1p7r;	t12=aj1p4i-aj1p7i;
	  t7 =t7+t9;			t8 =t8+t10;
	  t9 =t7+c3m1*t9;		t10=t8+c3m1*t10;
	  rt =s3*t11;			it =s3*t12;
	  t11=t9+it;			t12=t10-rt;
	  t9 =t9-it;			t10=t10+rt;
	#if PFETCH
	addr = add0+p2;
	prefetch_p_doubles(addr);
	#endif
	  t13=aj1p2r;			t14=aj1p2i;
	  t15=aj1p5r+aj1p8r;	t16=aj1p5i+aj1p8i;
	  t17=aj1p5r-aj1p8r;	t18=aj1p5i-aj1p8i;
	  t13=t13+t15;			t14=t14+t16;
	  t15=t13+c3m1*t15;		t16=t14+c3m1*t16;
	  rt =s3*t17;			it =s3*t18;
	  t17=t15+it;			t18=t16-rt;
	  t15=t15-it;			t16=t16+rt;
	#if PFETCH
	addr = add0+p3;
	prefetch_p_doubles(addr);
	#endif
/*
!...and now do three more radix-3 transforms, including the twiddle factors:
!                                           1, exp(i*2*pi/9), exp(i*4*pi/9) (for inputs to transform block 2)
!                                           1, exp(i*4*pi/9), exp(i*8*pi/9) (for inputs to transform block 3).
!   I.e. do similar as above, except inputs a(j1  +p0:8:1) are replaced by t1:17:2,
!                                           a(j2+p0:8:1) are replaced by t2:18:2, and v.v. for outputs,
!   and only the last 2 inputs to radix-3 transforms 2 and 3 are multiplied by non-unity twiddles.
*/
	  rt =t7;				it =t8;
	  t7 =rt+t13;			t8 =it+t14;
	  t13=rt-t13;			t14=it-t14;
	  t1 =t1+t7;			t2 =t2+t8;
	  a[j1   ]=t1;			a[j2   ]=t2;
	  t7 =t1+c3m1*t7;		t8 =t2+c3m1*t8;
	  rt =s3*t13;			it =s3*t14;
	  a[j1+p1]=t7-it;		a[j2+p1]=t8+rt;
	  a[j1+p2]=t7+it;		a[j2+p2]=t8-rt;
	#if PFETCH
	addr = add0+p4;
	prefetch_p_doubles(addr);
	#endif
	  rt =t9 *c -t10*s;		it =t9 *s +t10*c;		/* twiddle mpy by E^-1	*/
	  re =t15*c2-t16*s2;	t16=t15*s2+t16*c2;	t15=re;	/* twiddle mpy by E^-2	*/
	  t9 =rt+t15;			t10=it+t16;
	  t15=rt-t15;			t16=it-t16;
	  t3 =t3+t9;			t4 =t4+t10;
	#if PFETCH
	addr = add0+p5;
	prefetch_p_doubles(addr);
	#endif
	  a[j1+p3]=t3;			a[j2+p3]=t4;
	  t9 =t3+c3m1*t9;		t10=t4+c3m1*t10;
	  rt =s3*t15;			it =s3*t16;
	  a[j1+p4]=t9-it;		a[j2+p4]=t10+rt;
	  a[j1+p5]=t9+it;		a[j2+p5]=t10-rt;
	#if PFETCH
	addr = add0+p6;
	prefetch_p_doubles(addr);
	#endif
	  rt =t11*c2-t12*s2;	it =t11*s2+t12*c2;		/* twiddle mpy by E^-2	*/
	  re =t17*c4-t18*s4;	t18=t17*s4+t18*c4;	t17=re;	/* twiddle mpy by E^-4	*/
	  t11=rt+t17;			t12=it+t18;
	  t17=rt-t17;			t18=it-t18;
	  t5 =t5+t11;			t6 =t6+t12;
	#if PFETCH
	addr = add0+p7;
	prefetch_p_doubles(addr);
	#endif
	  a[j1+p6]=t5;			a[j2+p6]=t6;
	  t11=t5+c3m1*t11;		t12=t6+c3m1*t12;
	  rt =s3*t17;			it =s3*t18;
	  a[j1+p7]=t11-it;		a[j2+p7]=t12+rt;
	  a[j1+p8]=t11+it;		a[j2+p8]=t12-rt;
	#if PFETCH
	addr = add0+p8;
	prefetch_p_doubles(addr);
	#endif

#else
/*...gather the needed data and do a subconvolution-based transform. */

	  t1 =aj1p0r;			t2 =aj1p0i;		/* x0		*/
	  t3 =aj1p1r+aj1p8r;	t4 =aj1p1i+aj1p8i;	/* x1 + x8	*/
	  t5 =aj1p2r+aj1p7r;	t6 =aj1p2i+aj1p7i;	/* x2 + x7	*/
	  t7 =aj1p3r+aj1p6r;	t8 =aj1p3i+aj1p6i;	/* x3 + x6	*/
	  t9 =aj1p4r+aj1p5r;	t10=aj1p4i+aj1p5i;	/* x4 + x5	*/
	#if PFETCH
	addr = add0+p1;
	prefetch_p_doubles(addr);
	#endif
	  t11=aj1p4r-aj1p5r;	t12=aj1p4i-aj1p5i;	/* x4 - x5	*/
	  t13=aj1p3r-aj1p6r;	t14=aj1p3i-aj1p6i;	/* x3 - x6	*/
	  t15=aj1p2r-aj1p7r;	t16=aj1p2i-aj1p7i;	/* x2 - x7	*/
	  t17=aj1p1r-aj1p8r;	t18=aj1p1i-aj1p8i;	/* x1 - x8	*/
	#if PFETCH
	addr = add0+p2;
	prefetch_p_doubles(addr);
	#endif
	  rt = t1+t7;			it = t2+t8;
	  r2 = t3+t5+t9;		i2 = t4+t6+t10;
	  a[j1   ] = rt+r2;		a[j1+ 1] = it+i2;		/* X0	*/

	  re = t1-0.5*t7;		im = t2-0.5*t8;
	  t7 = rt-0.5*r2;		t8 = it-0.5*i2;			/* C3 calculated separately. */
	  t3 = t3-t5;			t4 = t4 -t6;
	  t5 = t9-t5;			t6 = t10-t6;
	#if PFETCH
	addr = add0+p3;
	prefetch_p_doubles(addr);
	#endif
	  t9 =(t3+t5)*cx3;		t10=(t4+t6)*cx3;
	  t3 = t3*cx1;			t4 = t4*cx1;
	  t5 = t5*cx2;			t6 = t6*cx2;
	  rt = t3-t9;			it = t4-t10;
	  t5 = t5-t9;			t6 = t6-t10;

	  t3 = re-rt-t5;		t4 = im-it-t6;	/* C2 */
	  t5 = re+t5;			t6 = im+t6;	/* C4 */
	  t1 = re+rt;			t2 = im+it;	/* C1 */
	#if PFETCH
	addr = add0+p4;
	prefetch_p_doubles(addr);
	#endif
	/* For sine convo, replace t3,t5,t9 --> t17,-t15,t11 on input */
	  re = t13*ss3;			im = t14*ss3;
	  t13= (t17-t15+t11)*ss3;t14= (t18-t16+t12)*ss3;		/* S3 calculated separately. */
	  t17= t17+t15;			t18= t18+t16;
	  t15= t11+t15;			t16= t12+t16;
	#if PFETCH
	addr = add0+p5;
	prefetch_p_doubles(addr);
	#endif
	  t11=(t17+t15)*sx3;	t12=(t18+t16)*sx3;
	  t17= t17*sx1;			t18= t18*sx1;
	  t15= t15*sx2;			t16= t16*sx2;
	  rt = t17-t11;			it = t18-t12;
	  t15= t15-t11;			t16= t16-t12;

	  t17= re-rt-t15;		t18= im-it-t16;	/*-S2 */
	  t15= re+t15;			t16= im+t16;	/* S4 */
	  t9 = re+rt;			t10= im+it;	/* S1 */
	#if PFETCH
	addr = add0+p6;
	prefetch_p_doubles(addr);
	#endif
/*...Inline multiply of sine parts by +-I into finishing phase.
     We desire an output order [0,3,6,1,4,7,2,5,8], so scramble array indices accordingly. */

	  a[j1+p3]=t1-t10;		a[j2+p3]=t2+t9 ;	/* X1 = C1 + I*S1	*/
	  a[j1+p6]=t3+t18;		a[j2+p6]=t4-t17;	/* X2 = C2 + I*S2	*/
	  a[j1+p1]=t7-t14;		a[j2+p1]=t8+t13;	/* X3 = C3 + I*S3	*/
	  a[j1+p4]=t5-t16;		a[j2+p4]=t6+t15;	/* X4 = C4 + I*S4	*/
	#if PFETCH
	addr = add0+p7;
	prefetch_p_doubles(addr);
	#endif
	  a[j1+p7]=t5+t16;		a[j2+p7]=t6-t15;	/* X5 =	C4 - I*S4	*/
	  a[j1+p2]=t7+t14;		a[j2+p2]=t8-t13;	/* X6 =	C3 - I*S3	*/
	  a[j1+p5]=t3-t18;		a[j2+p5]=t4+t17;	/* X7 =	C2 - I*S2	*/
	  a[j1+p8]=t1+t10;		a[j2+p8]=t2-t9 ;	/* X8 =	C1 - I*S1	*/
	#if PFETCH
	addr = add0+p8;
	prefetch_p_doubles(addr);
	#endif

#endif
	  }

	  jstart += nwt;
	  jhi    += nwt;
	  col += 9;
	  co3 -= 9;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-9 forward DIF FFT of the first block of 9 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 9 outputs of (1);
!   (3) Reweight and perform a radix-9 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 9 elements and repeat (1-4).
*/
	t1  = cy8;
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
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8) != 0.0)
	{
	    sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix9_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
	    mlucas_fprint(cbuf,INTERACT);
	    err=ERR_CARRY;
	    return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix9_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-9 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Uses an optimized [radix-3]^2 transform with twiddles (if LO_ADD is set),
!   or a subconvolution-based version similar to that used for prime radix 7 (if not).
!
!   DFT multiplier matrix is (with E := exp(i*2*pi/9), using that E^j = *E^(9-j) and E^k = E^(k mod 9))
!
!	1	 1	 1	 1	 1	 1	 1	 1	 1
!	1	 E	 E^2	 E^3	 E^4	*E^4	*E^3	*E^2	*E
!	1	 E^2	 E^4	*E^3	*E	 E	 E^3	*E^4	*E^2
!	1	 E^3	*E^3	 1	 E^3	*E^3	 1	 E^3	*E^3
!	1	 E^4	*E	 E^3	*E^2	 E^2	*E^3	 E	*E^4
!	1	*E^4	 E	*E^3	 E^2	*E^2	 E^3	*E	 E^4
!	1	*E^3	 E^3	 1	*E^3	 E^3	 1	*E^3	 E^3
!	1	*E^2	*E^4	 E^3	 E	*E	*E^3	 E^4	 E^2
!	1	*E	*E^2	*E^3	*E^4	 E^4	 E^3	 E^2	 E
!
!   Totals : LO_ADD = 0: 20 FMUL, 84 FADD (same #mul, 4 less add than Nussbaumer.)
!            LO_ADD = 1: 40 FMUL, 80 FADD.
!
*/
	int j,j1,j2;
	static int n9,p1,p2,p3,p4,p5,p6,p7,p8, first_entry=TRUE;
#if LO_ADD
	static double c   =  0.76604444311897803520,	/* cos(2*pi/9) */
		      s   =  0.64278760968653932631,	/* sin(2*pi/9) */
		      c2  =  0.17364817766693034887,	/* cos(2*u) */
		      s2  =  0.98480775301220805936,	/* sin(2*u) */
		      c3m1= -1.50000000000000000000,	/* cos(3*u)-1] */
		      s3  =  0.86602540378443864677,	/* sin(3*u)] */
		      c4  = -0.93969262078590838404,	/* cos(4*u) */
		      s4  =  0.34202014332566873307;	/* sin(4*u) */
#else
	static double /*cx0 = 0.00000000000000000000,	 (cc1+cc2+cc4)/3 = 0	*/
		      cx1 = 1.70573706390488641924, 	/*  cc1-cc4		*/
		      cx2 = 1.11334079845283873291, 	/*  cc2-cc4		*/
		      cx3 = 0.93969262078590838405,	/* (cc1+cc2-2*cc4)/3 = -cc4	*/
			/*
			Cosine subconvo is of form

			|c1 c2 c3 c4| |x0|
			|c2 c4 c3 c1|*|x1|
			|c3 c3 1  c3| |x2|
			|c4 c1 c3 c2| |x3|

			We'll calculate the convolution-wrecking c2/x2 terms separately.
			Deleting row 3 and col 3 of the matrix gives

			|c1 c2 c4| |x0|
			|c2 c4 c1|*|x1|
			|c4 c1 c2| |x3|

			Which is in the form a length-3 cyclic convo if we re-index (d0,d1,d2) = (c1,c2,c4)
			and swap the x1 and x3 inputs (i.e. swap cols 2 and 3 of the 3x3 multiplier matrix.)
			*/

		/* Switch the sign of ss2 in these: */
		      /*sx0 = 0.00000000000000000000,	 (ss1-ss2+ss4)/3	*/
		      sx1 = 0.30076746636087059328, 	/*  ss1-ss4		*/
		      sx2 =-1.32682789633787679241, 	/* -ss2-ss4		*/
		      s3  = 0.86602540378443864676,	/*  ss3			*/
		      sx3 =-0.34202014332566873304;	/* (ss1-ss2-2*ss4)/3	*/
			/*
			Sine subconvo is of form

			|s1  s2  s3  s4| |y0|
			|s2  s4 -s3 -s1|*|y1|
			|s3 -s3  0   s3| |y2|
			|s4 -s1  s3 -s2| |y3|

			We'll calculate the convolution-wrecking s2/y2 terms separately.
			Deleting row 3 and col 3 of the matrix gives

			|s1  s2  s4| |y0|
			|s2  s4 -s1|*|y1|
			|s4 -s1 -s2| |y3|

			Which is in the form a length-3 cyclic convo if we set (e0,e1,e2) = (s1,-s2,s4)
			and (z0,z1,z2) = (y0,y3,-y1) and negate the sign of the second convo output:

			| e0  e2  e1| |z0|
			|-e1 -e0 -e2|*|z1|
			| e2  e1  e0| |z2|
			*/
	double im,r2,i2;
#endif
	double rt,it,re
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;

	if(!first_entry && (n/9) != n9)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n9=n/9;

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n9;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;
	  p6 = p5 +p1;
	  p7 = p6 +p1;
	  p8 = p7 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	  p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	  p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
	  p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
	  p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-9 pass is here.	*/

      for(j=0; j < n9; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

#if LO_ADD
    #if 1
		/*...gather the needed data (9 64-bit complex, i.e. 18 64-bit reals) and do three radix-3 transforms...	*/
		t1 =a[j1   ];			t2 =a[j2   ];
		t3 =a[j1+p3]+a[j1+p6];t4 =a[j2+p3]+a[j2+p6];
		t5 =a[j1+p3]-a[j1+p6];t6 =a[j2+p3]-a[j2+p6];
		t1 =t1+t3;			t2 =t2+t4;
		t3 =t1+c3m1*t3;		t4 =t2+c3m1*t4;
		rt =s3*t5;			it =s3*t6;
		t5 =t3+it;			t6 =t4-rt;
		t3 =t3-it;			t4 =t4+rt;

		t7 =a[j1+p1];			t8 =a[j2+p1];
		t9 =a[j1+p4]+a[j1+p7];t10=a[j2+p4]+a[j2+p7];
		t11=a[j1+p4]-a[j1+p7];t12=a[j2+p4]-a[j2+p7];
		t7 =t7+t9;			t8 =t8+t10;
		t9 =t7+c3m1*t9;		t10=t8+c3m1*t10;
		rt =s3*t11;			it =s3*t12;
		t11=t9+it;			t12=t10-rt;
		t9 =t9-it;			t10=t10+rt;

		t13=a[j1+p2];			t14=a[j2+p2];
		t15=a[j1+p5]+a[j1+p8];t16=a[j2+p5]+a[j2+p8];
		t17=a[j1+p5]-a[j1+p8];t18=a[j2+p5]-a[j2+p8];
		t13=t13+t15;			t14=t14+t16;
		t15=t13+c3m1*t15;		t16=t14+c3m1*t16;
		rt =s3*t17;			it =s3*t18;
		t17=t15+it;			t18=t16-rt;
		t15=t15-it;			t16=t16+rt;
		/*
		!...and now do three more radix-3 transforms, including the twiddle factors:
		!                                           1, exp(i*2*pi/9), exp(i*4*pi/9) (for inputs to transform block 2)
		!                                           1, exp(i*4*pi/9), exp(i*8*pi/9) (for inputs to transform block 3).
		!   I.e. do similar as above, except inputs a[j1  +p0:8:1) are replaced by t1:17:2,
		!                                           a[j2+p0:8:1) are replaced by t2:18:2, and v.v. for outputs,
		!   and only the last 2 inputs to radix-3 transforms 2 and 3 are multiplied by non-unity twiddles.
		!   t1,2/7,8/13,14:
		*/
		rt =t7;				it =t8;
		t7 =rt+t13;			t8 =it+t14;
		t13=rt-t13;			t14=it-t14;
		t1 =t1+t7;			t2 =t2+t8;
		a[j1   ]=t1;			a[j2   ]=t2;
		t7 =t1+c3m1*t7;		t8 =t2+c3m1*t8;
		rt =s3*t13;			it =s3*t14;
		a[j1+p1]=t7-it;		a[j2+p1]=t8+rt;
		a[j1+p2]=t7+it;		a[j2+p2]=t8-rt;
		/*  t3,4/9,10/15,16:	*/
		rt =t9 *c -t10*s;		it =t9 *s +t10*c;		/* twiddle mpy by E	*/
		re =t15*c2-t16*s2;	t16=t15*s2+t16*c2;	t15=re;	/* twiddle mpy by E^2	*/
		t9 =rt+t15;			t10=it+t16;
		t15=rt-t15;			t16=it-t16;
		t3 =t3+t9;			t4 =t4+t10;
		a[j1+p3]=t3;			a[j2+p3]=t4;
		t9 =t3+c3m1*t9;		t10=t4+c3m1*t10;
		rt =s3*t15;			it =s3*t16;
		a[j1+p4]=t9-it;		a[j2+p4]=t10+rt;
		a[j1+p5]=t9+it;		a[j2+p5]=t10-rt;
		/*   5,6/11,12/17,18:	*/
		rt =t11*c2-t12*s2;	it =t11*s2+t12*c2;		/* twiddle mpy by E^2	*/
		re =t17*c4-t18*s4;	t18=t17*s4+t18*c4;	t17=re;	/* twiddle mpy by E^4	*/
		t11=rt+t17;			t12=it+t18;
		t17=rt-t17;			t18=it-t18;
		t5 =t5+t11;			t6 =t6+t12;
		a[j1+p6]=t5;			a[j2+p6]=t6;
		t11=t5+c3m1*t11;		t12=t6+c3m1*t12;
		rt =s3*t17;			it =s3*t18;
		a[j1+p7]=t11-it;		a[j2+p7]=t12+rt;
		a[j1+p8]=t11+it;		a[j2+p8]=t12-rt;
    #else
//*** The radix-9 non-FMA macros were designed as a basis for radix-36 and higher ... in standalone mode they return permuted outputs!
		// The non-FMA version of the RADIX_09_DIF macro does not support in-place-ness, so copy inputs into temps first:
		t1  = a[j1   ];	t2  = a[j2   ];
		t3  = a[j1+p1];	t4  = a[j2+p1];
		t5  = a[j1+p2];	t6  = a[j2+p2];
		t7  = a[j1+p3];	t8  = a[j2+p3];
		t9  = a[j1+p4];	t10 = a[j2+p4];
		t11 = a[j1+p5];	t12 = a[j2+p5];
		t13 = a[j1+p6];	t14 = a[j2+p6];
		t15 = a[j1+p7];	t16 = a[j2+p7];
		t17 = a[j1+p8];	t18 = a[j2+p8];
		RADIX_09_DIF(
			t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],
            rt,it,re);
	/*
		// FMA version is OK with in-place-ness:
		RADIX_09_DIF_FMA(
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8]
			);
	*/
    #endif
#else
/*...gather the needed data and do a subconvolution-based transform. */

	  t1 =a[j1   ];				t2 =a[j2   ];			/* x0		*/
	  t3 =a[j1+p1]+a[j1+p8];	t4 =a[j2+p1]+a[j2+p8];	/* x1 + x8	*/
	  t5 =a[j1+p2]+a[j1+p7];	t6 =a[j2+p2]+a[j2+p7];	/* x2 + x7	*/
	  t7 =a[j1+p3]+a[j1+p6];	t8 =a[j2+p3]+a[j2+p6];	/* x3 + x6	*/
	  t9 =a[j1+p4]+a[j1+p5];	t10=a[j2+p4]+a[j2+p5];	/* x4 + x5	*/
	  t11=a[j1+p4]-a[j1+p5];	t12=a[j2+p4]-a[j2+p5];	/* x4 - x5	*/
	  t13=a[j1+p3]-a[j1+p6];	t14=a[j2+p3]-a[j2+p6];	/* x3 - x6	*/
	  t15=a[j1+p2]-a[j1+p7];	t16=a[j2+p2]-a[j2+p7];	/* x2 - x7	*/
	  t17=a[j1+p1]-a[j1+p8];	t18=a[j2+p1]-a[j2+p8];	/* x1 - x8	*/
/*
printf("Method A Cosine terms:\n");
printf("%20.15f  %20.15f\n", t1+c *t3+c2*t5-.5*t7+c4*t9, t2+c *t4+c2*t6-.5*t8+c4*t10);
printf("%20.15f  %20.15f\n", t1+c2*t3+c4*t5-.5*t7+c *t9, t2+c2*t4+c4*t6-.5*t8+c *t10);
printf("%20.15f  %20.15f\n", t1-.5*t3-.5*t5+   t7-.5*t9, t2-.5*t4-.5*t6+   t8-.5*t10);
printf("%20.15f  %20.15f\n", t1+c4*t3+c *t5-.5*t7+c2*t9, t2+c4*t4+c *t6-.5*t8+c2*t10);
printf("Method A Sine terms:\n");
printf("%20.15f  %20.15f\n", s *t17+s2*t15+s3*t13+s4*t11, s *t18+s2*t16+s3*t14+s4*t12);
printf("%20.15f  %20.15f\n", s2*t17+s4*t15-s3*t13-s *t11, s2*t18+s4*t16-s3*t14-s *t12);
printf("%20.15f  %20.15f\n", s3*t17-s3*t15       +s3*t11, s3*t18-s3*t16       +s3*t12);
printf("%20.15f  %20.15f\n", s4*t17-s *t15+s3*t13-s2*t11, s4*t18-s *t16+s3*t14-s2*t12);
*/
	  rt = t1+t7;				it = t2+t8;
	  r2 = t3+t5+t9;			i2 = t4+t6+t10;
	  a[j1   ] = rt+r2;			a[j2   ] = it+i2;		/* X0	*/

	  re = t1-0.5*t7;			im = t2-0.5*t8;
	  t7 = rt-0.5*r2;			t8 = it-0.5*i2;			/* C3 calculated separately. */
	  t3 = t3-t5;				t4 = t4 -t6;
	  t5 = t9-t5;				t6 = t10-t6;
	  t9 =(t3+t5)*cx3;			t10=(t4+t6)*cx3;
	  t3 = t3*cx1;				t4 = t4*cx1;
	  t5 = t5*cx2;				t6 = t6*cx2;
	  rt = t3-t9;				it = t4-t10;
	  t5 = t5-t9;				t6 = t6-t10;

	  t3 = re-rt-t5;			t4 = im-it-t6;	/* C2 */
	  t5 = re+t5;				t6 = im+t6;	/* C4 */
	  t1 = re+rt;				t2 = im+it;	/* C1 */
/*
printf("Method B Cosine terms:\n");
printf("%20.15f  %20.15f\n", t1 ,t2 );
printf("%20.15f  %20.15f\n", t3 ,t4 );
printf("%20.15f  %20.15f\n", t7 ,t8 );
printf("%20.15f  %20.15f\n", t5 ,t6 );
*/
	/* For sine convo, replace t3,t5,t9 --> t17,-t15,t11 on input */
	  re = t13*s3;				im = t14*s3;
	  t13= (t17-t15+t11)*s3;	t14= (t18-t16+t12)*s3;		/* S3 calculated separately. */
	  t17= t17+t15;				t18= t18+t16;
	  t15= t11+t15;				t16= t12+t16;
	  t11=(t17+t15)*sx3;		t12=(t18+t16)*sx3;
	  t17= t17*sx1;				t18= t18*sx1;
	  t15= t15*sx2;				t16= t16*sx2;
	  rt = t17-t11;				it = t18-t12;
	  t15= t15-t11;				t16= t16-t12;

	  t17= re-rt-t15;			t18= im-it-t16;	/*-S2 */
	  t15= re+t15;				t16= im+t16;	/* S4 */
	  t9 = re+rt;				t10= im+it;	/* S1 */
/*
printf("Method B Sine terms:\n");
printf("%20.15f  %20.15f\n", t9 , t10);
printf("%20.15f  %20.15f\n", t17, t18);
printf("%20.15f  %20.15f\n", t13, t14);
printf("%20.15f  %20.15f\n", t15, t16);
printf("");
*/
/*...Inline multiply of sine parts by +-I into finishing phase.
     We desire an output order [0,3,6,1,4,7,2,5,8], so scramble array indices accordingly. */

	  a[j1+p3]=t1-t10;			a[j2+p3]=t2+t9 ;	/* X1 = C1 + I*S1	*/
	  a[j1+p6]=t3+t18;			a[j2+p6]=t4-t17;	/* X2 = C2 + I*S2	*/
	  a[j1+p1]=t7-t14;			a[j2+p1]=t8+t13;	/* X3 = C3 + I*S3	*/
	  a[j1+p4]=t5-t16;			a[j2+p4]=t6+t15;	/* X4 = C4 + I*S4	*/
	  a[j1+p7]=t5+t16;			a[j2+p7]=t6-t15;	/* X5 =	C4 - I*S4	*/
	  a[j1+p2]=t7+t14;			a[j2+p2]=t8-t13;	/* X6 =	C3 - I*S3	*/
	  a[j1+p5]=t3-t18;			a[j2+p5]=t4+t17;	/* X7 =	C2 - I*S2	*/
	  a[j1+p8]=t1+t10;			a[j2+p8]=t2-t9 ;	/* X8 =	C1 - I*S1	*/
#endif
	}
}

/***************/

void radix9_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-9 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix9_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int j,j1,j2;
	static int n9,p1,p2,p3,p4,p5,p6,p7,p8, first_entry=TRUE;
#if LO_ADD
	static double c   =  0.76604444311897803520,	/* cos(2*pi/9) */
		      s   =  0.64278760968653932631,	/* sin(2*pi/9) */
		      c2  =  0.17364817766693034887,	/* cos(2*u) */
		      s2  =  0.98480775301220805936,	/* sin(2*u) */
		      c3m1= -1.50000000000000000000,	/* cos(3*u)-1] */
		      s3  =  0.86602540378443864677,	/* sin(3*u)] */
		      c4  = -0.93969262078590838404,	/* cos(4*u) */
		      s4  =  0.34202014332566873307;	/* sin(4*u) */
#else
	static double cx1 = 1.70573706390488641924, 	/*  cc1-cc4		*/
		      cx2 = 1.11334079845283873291, 	/*  cc2-cc4		*/
		      cx3 = 0.93969262078590838405,	/* (cc1+cc2-2*cc4)/3	*/
		/* Switch the sign of ss2 in these: */
		      sx1 = 0.30076746636087059324, 	/*  ss1-ss4		*/
		      sx2 =-1.32682789633787679243, 	/* -ss2-ss4		*/
		      s3  = 0.86602540378443864677,	/*  ss3			*/
		      sx3 =-0.34202014332566873306;	/* (ss1-ss2-2*ss4)/3	*/
	double im,r2,i2;
#endif
	double rt,it,re
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;

	if(!first_entry && (n/9) != n9)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n9=n/9;

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n9;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;
	  p6 = p5 +p1;
	  p7 = p6 +p1;
	  p8 = p7 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	  p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	  p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
	  p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
	  p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-9 pass is here.	*/

      for(j=0; j < n9; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

#if LO_ADD
    #if 1
/*...gather the needed data (9 64-bit complex, i.e. 18 64-bit reals) and do three radix-3 transforms...	*/
	  t1 =a[j1   ];				t2 =a[j2   ];
	  t3 =a[j1+p1]+a[j1+p2];	t4 =a[j2+p1]+a[j2+p2];
	  t5 =a[j1+p1]-a[j1+p2];	t6 =a[j2+p1]-a[j2+p2];
	  t1 =t1+t3;				t2 =t2+t4;
	  t3 =t1+c3m1*t3;			t4 =t2+c3m1*t4;
	  rt =s3*t5;				it =s3*t6;
	  t5 =t3-it;				t6 =t4+rt;
	  t3 =t3+it;				t4 =t4-rt;

	  t7 =a[j1+p3];				t8 =a[j2+p3];
	  t9 =a[j1+p4]+a[j1+p5];	t10=a[j2+p4]+a[j2+p5];
	  t11=a[j1+p4]-a[j1+p5];	t12=a[j2+p4]-a[j2+p5];
	  t7 =t7+t9;				t8 =t8+t10;
	  t9 =t7+c3m1*t9;			t10=t8+c3m1*t10;
	  rt =s3*t11;				it =s3*t12;
	  t11=t9-it;				t12=t10+rt;
	  t9 =t9+it;				t10=t10-rt;

	  t13=a[j1+p6];				t14=a[j2+p6];
	  t15=a[j1+p7]+a[j1+p8];	t16=a[j2+p7]+a[j2+p8];
	  t17=a[j1+p7]-a[j1+p8];	t18=a[j2+p7]-a[j2+p8];
	  t13=t13+t15;				t14=t14+t16;
	  t15=t13+c3m1*t15;			t16=t14+c3m1*t16;
	  rt =s3*t17;				it =s3*t18;
	  t17=t15-it;				t18=t16+rt;
	  t15=t15+it;				t16=t16-rt;
/*
!...and now do three more radix-3 transforms, including the twiddle factors:
!                                           1, exp(-i*2*pi/9), exp(-i*4*pi/9) (for inputs to transform block 2)
!                                           1, exp(-i*4*pi/9), exp(-i*8*pi/9) (for inputs to transform block 3).
!   I.e. do similar as above, except inputs a[j1  +p0:8:1) are replaced by t1:17:2,
!                                           a[j2+p0:8:1) are replaced by t2:18:2, and v.v. for outputs,
!   and only the last 2 inputs to radix-3 transforms 2 and 3 are multiplied by non-unity twiddles.
!   t1,2/7,8/13,14:
*/
	  rt =t7;					it =t8;
	  t7 =rt+t13;				t8 =it+t14;
	  t13=rt-t13;				t14=it-t14;
	  t1 =t1+t7;				t2 =t2+t8;
	  a[j1   ]=t1;				a[j2   ]=t2;
	  t7 =t1+c3m1*t7;			t8 =t2+c3m1*t8;
	  rt =s3*t13;				it =s3*t14;
	  a[j1+p3]=t7+it;			a[j2+p3]=t8-rt;
	  a[j1+p6]=t7-it;			a[j2+p6]=t8+rt;
/*  t3,4/9,10/15,16:	*/
	  rt =t9 *c +t10*s;			it =t10*c -t9 *s;		/* twiddle mpy by E^-1	*/
	  re =t15*c2+t16*s2;		t16=t16*c2-t15*s2;	t15=re;	/* twiddle mpy by E^-2	*/
	  t9 =rt+t15;				t10=it+t16;
	  t15=rt-t15;				t16=it-t16;
	  t3 =t3+t9;				t4 =t4+t10;
	  a[j1+p1]=t3;				a[j2+p1]=t4;
	  t9 =t3+c3m1*t9;			t10=t4+c3m1*t10;
	  rt =s3*t15;				it =s3*t16;
	  a[j1+p4]=t9+it;			a[j2+p4]=t10-rt;
	  a[j1+p7]=t9-it;			a[j2+p7]=t10+rt;
/*  5,6/11,12/17,18:	*/
	  rt =t11*c2+t12*s2;		it =t12*c2-t11*s2;		/* twiddle mpy by E^-2	*/
	  re =t17*c4+t18*s4;		t18=t18*c4-t17*s4;	t17=re;	/* twiddle mpy by E^-4	*/
	  t11=rt+t17;				t12=it+t18;
	  t17=rt-t17;				t18=it-t18;
	  t5 =t5+t11;				t6 =t6+t12;
	  a[j1+p2]=t5;				a[j2+p2]=t6;
	  t11=t5+c3m1*t11;			t12=t6+c3m1*t12;
	  rt =s3*t17;				it =s3*t18;
	  a[j1+p5]=t11+it;			a[j2+p5]=t12-rt;
	  a[j1+p8]=t11-it;			a[j2+p8]=t12+rt;
    #else
//*** The radix-9 non-FMA macros were designed as a basis for radix-36 and higher ... in standalone mode they return permuted outputs!
		// The non-FMA version of the RADIX_09_DIF macro does not support in-place-ness, so copy inputs into temps first:
		t1  = a[j1   ];	t2  = a[j2   ];
		t3  = a[j1+p1];	t4  = a[j2+p1];
		t5  = a[j1+p2];	t6  = a[j2+p2];
		t7  = a[j1+p3];	t8  = a[j2+p3];
		t9  = a[j1+p4];	t10 = a[j2+p4];
		t11 = a[j1+p5];	t12 = a[j2+p5];
		t13 = a[j1+p6];	t14 = a[j2+p6];
		t15 = a[j1+p7];	t16 = a[j2+p7];
		t17 = a[j1+p8];	t18 = a[j2+p8];
		RADIX_09_DIT(
			t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],
            rt,it,re);
	/*
		// FMA version is OK with in-place-ness:
		RADIX_09_DIT_FMA(
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8]
			);
	*/
    #endif
#else
/*...gather the needed data (11 64-bit complex, i.e. 22 64-bit reals).
     We expect an input order [0,3,6,1,4,7,2,5,8], so scramble array indices accordingly. */
	  t1 =a[j1   ];				t2 =a[j2   ];			/* x0		*/
	  t3 =a[j1+p3]+a[j1+p8];	t4 =a[j2+p3]+a[j2+p8];	/* x1 + x8	*/
	  t5 =a[j1+p6]+a[j1+p5];	t6 =a[j2+p6]+a[j2+p5];	/* x2 + x7	*/
	  t7 =a[j1+p1]+a[j1+p2];	t8 =a[j2+p1]+a[j2+p2];	/* x3 + x6	*/
	  t9 =a[j1+p4]+a[j1+p7];	t10=a[j2+p4]+a[j2+p7];	/* x4 + x5	*/
	  t11=a[j1+p4]-a[j1+p7];	t12=a[j2+p4]-a[j2+p7];	/* x4 - x5	*/
	  t13=a[j1+p1]-a[j1+p2];	t14=a[j2+p1]-a[j2+p2];	/* x3 - x6	*/
	  t15=a[j1+p6]-a[j1+p5];	t16=a[j2+p6]-a[j2+p5];	/* x2 - x7	*/
	  t17=a[j1+p3]-a[j1+p8];	t18=a[j2+p3]-a[j2+p8];	/* x1 - x8	*/

	  rt = t1+t7;				it = t2+t8;
	  r2 = t3+t5+t9;			i2 = t4+t6+t10;
	  a[j1   ] = rt+r2;			a[j1+ 1] = it+i2;		/* X0	*/

	  re = t1-0.5*t7;			im = t2-0.5*t8;
	  t7 = rt-0.5*r2;			t8 = it-0.5*i2;			/* C3 calculated separately. */
	  t3 = t3-t5;				t4 = t4 -t6;
	  t5 = t9-t5;				t6 = t10-t6;
	  t9 =(t3+t5)*cx3;			t10=(t4+t6)*cx3;
	  t3 = t3*cx1;				t4 = t4*cx1;
	  t5 = t5*cx2;				t6 = t6*cx2;
	  rt = t3-t9;				it = t4-t10;
	  t5 = t5-t9;				t6 = t6-t10;

	  t3 = re-rt-t5;			t4 = im-it-t6;	/* C2 */
	  t5 = re+t5;				t6 = im+t6;		/* C4 */
	  t1 = re+rt;				t2 = im+it;		/* C1 */
	/* For sine convo, replace t3,t5,t9 --> t17,-t15,t11 on input */
	  re = t13*s3;				im = t14*s3;
	  t13= (t17-t15+t11)*s3;	t14= (t18-t16+t12)*s3;		/* S3 calculated separately. */
	  t17= t17+t15;				t18= t18+t16;
	  t15= t11+t15;				t16= t12+t16;
	  t11=(t17+t15)*sx3;		t12=(t18+t16)*sx3;
	  t17= t17*sx1;				t18= t18*sx1;
	  t15= t15*sx2;				t16= t16*sx2;
	  rt = t17-t11;				it = t18-t12;
	  t15= t15-t11;				t16= t16-t12;

	  t17= re-rt-t15;			t18= im-it-t16;	/*-S2 */
	  t15= re+t15;				t16= im+t16;	/* S4 */
	  t9 = re+rt;				t10= im+it;		/* S1 */

/*...Inline multiply of sine parts by +-I into finishing phase.	Outputs are ordered. */

	  a[j1+p1]=t1+t10;			a[j2+p1]=t2-t9 ;	/* X1 = C1 - I*S1	*/
	  a[j1+p2]=t3-t18;			a[j2+p2]=t4+t17;	/* X2 = C2 - I*S2	*/
	  a[j1+p3]=t7+t14;			a[j2+p3]=t8-t13;	/* X3 = C3 - I*S3	*/
	  a[j1+p4]=t5+t16;			a[j2+p4]=t6-t15;	/* X4 = C4 - I*S4	*/
	  a[j1+p5]=t5-t16;			a[j2+p5]=t6+t15;	/* X5 =	C4 + I*S4	*/
	  a[j1+p6]=t7-t14;			a[j2+p6]=t8+t13;	/* X6 =	C3 + I*S3	*/
	  a[j1+p7]=t3+t18;			a[j2+p7]=t4-t17;	/* X7 =	C2 + I*S2	*/
	  a[j1+p8]=t1-t10;			a[j2+p8]=t2+t9 ;	/* X8 =	C1 + I*S1	*/
#endif
	}
}

