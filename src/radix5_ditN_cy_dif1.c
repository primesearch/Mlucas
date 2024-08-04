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

int radix5_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-5 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-5 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n5,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
		      cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
		      s   =  0.95105651629515357211,	/*  sin(u) */
		      ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
		      ss2 =  0.36327126400268044292,	/* [sin(u)-sin(2u)] */
		      radix_inv, n2inv;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i
		,cy0,cy1,cy2,cy3,cy4,temp,frac,scale;
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

/*...change n5 and n_div_wt to non-static to work around a gcc compiler bug. */
	n5   = n/5;
	n_div_nwt = n5 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n5)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/5 in radix5_ditN_cy_dif1.\n",iter);
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
	  radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)5));
	  n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

	  bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
	  sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n5;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );

	  bjmodnini=0;
	  for(j=0; j < n5; j++)
	  {
	    bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
	  }
	}

/*...The radix-5 final DIT pass is here.	*/

	/* init carries	*/
	cy0 = 0;
	cy1 = 0;
	cy2 = 0;
	cy3 = 0;
	cy4 = 0;

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

	col=0;
	co2=(n >> nwt_bits)-1+5;
	co3=co2-5;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

	    t3 =a[j1+p1];		t4 =a[j2+p1];
	    rt =a[j1+p4];		it =a[j2+p4];
	    t9 =t3 -rt;			t10=t4 -it;
	    t3 =t3 +rt;			t4 =t4 +it;

	    t5 =a[j1+p2];		t6 =a[j2+p2];
	    rt =a[j1+p3];		it =a[j2+p3];
	    t7 =t5 -rt;			t8 =t6 -it;
	    t5 =t5 +rt;			t6 =t6 +it;

/*       ...now complete the result.	*/

	    rt = t3+t5;			it = t4+t6;
	    t1 = t1+rt;			t2 = t2+it;
	    aj1p0r = t1;		aj1p0i = t2;

	    rt = t1+cc1*rt;		it = t2+cc1*it;
	    t5 = cc2*(t3-t5);	t6 = cc2*(t4-t6);

	    t3 = rt+t5;			t4 = it+t6;
	    t5 = rt-t5;			t6 = it-t6;

/*	...only difference from the forward transform is that in this next part, we inline minus signs in front of s,ss1,ss2.	*/

	    rt = s  *(t7-t9);	it = s  *(t8-t10);
	    t7 = ss1* t7;		t8 = ss1* t8;
	    t9 = ss2* t9;		t10= ss2* t10;

	    t7 = rt-t7;			t8 = it-t8;
	    t9 = rt+t9;			t10= it+t10;

	    aj1p1r=t3-t8;		aj1p1i=t4+t7;
	    aj1p2r=t5-t10;		aj1p2i=t6+t9;
	    aj1p3r=t5+t10;		aj1p3i=t6-t9;
	    aj1p4r=t3+t8;		aj1p4i=t4-t7;

/*...and combine those to complete the radix-5 transform and do the carries. Since the outputs would
     normally be getting dispatched to 5 separate blocks of the A-array, we need 5 separate carries.	*/

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

	    i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
	    co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			   and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-5 DIF pass is here:	*/

/*       gather the needed data (5 64-bit complex, i.e. 10 64-bit reals) and begin the transform...	*/
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif
	  t1 =aj1p0r;			t2 =aj1p0i;
	  t3 =aj1p1r+aj1p4r;	t4 =aj1p1i+aj1p4i;
	  t5 =aj1p2r+aj1p3r;	t6 =aj1p2i+aj1p3i;
	  t7 =aj1p2r-aj1p3r;	t8 =aj1p2i-aj1p3i;
	  t9 =aj1p1r-aj1p4r;	t10=aj1p1i-aj1p4i;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
/*       ...now complete the first complex datum of the result.	*/

	  rt = t3+t5;			it = t4+t6;
	  t1 = t1+rt;			t2 = t2+it;

	  rt = t1+cc1*rt;		it = t2+cc1*it;
	  t5 = cc2*(t3-t5);		t6 = cc2*(t4-t6);

	  t3 = rt+t5;			t4 = it+t6;
	  t5 = rt-t5;			t6 = it-t6;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif
	  rt = s  *(t9-t7);		it = s  *(t10-t8);
	  t7 = ss1* t7;			t8 = ss1* t8;
	  t9 = ss2* t9;			t10= ss2* t10;

	  t7 = rt+t7;			t8 = it+t8;
	  t9 = rt-t9;			t10= it-t10;

/*...Inline multiply of sine part by I into finishing phase...	*/
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
	  a[j1   ]=t1;			a[j2   ]=t2;
	  a[j1+p1]=t3-t8;		a[j2+p1]=t4+t7;
	  a[j1+p2]=t5-t10;		a[j2+p2]=t6+t9;
	  a[j1+p3]=t5+t10;		a[j2+p3]=t6-t9;
	  a[j1+p4]=t3+t8;		a[j2+p4]=t4-t7;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
	  }

	  jstart += nwt;
	  jhi    += nwt;
	  col += 5;
	  co3 -= 5;

	}

	if(root_incr==0)break;

/*  Wraparound carry cleanup loop is here: ***
!
!   The cleanup carries from the end of each length-N/5 block into the begining of the next
!   can all be neatly processed as follows:
!
!   (1) Invert the radix-5 forward DIF FFT of the first block of 5 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 5 outputs of (1);
!   (3) Reweight and perform a radix-5 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 5 elements and repeat (1-4).
*/
	t1 =cy4;
	cy4=cy3;
	cy3=cy2;
	cy2=cy1;
	cy1=cy0;
	cy0=t1;

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
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4) != 0.0)
	{
	    sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix5_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
	    mlucas_fprint(cbuf,INTERACT);
	    err=ERR_CARRY;
	    return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix5_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-5 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Uses an optimized transform a la Nussbaumer (2nd ed., p.146).
!
!   Given complex input data (x0,x1,x2,x3,x4), we need the following combinations on output
!   (here c1 = cos(2*pi/5), s1 = sin(2*pi/5),  c2 = cos(4*pi/5), s2 = sin(4*pi/5)):
!
!	X0 = x0+x1+x2+x3+x4
!	X1 = x0+c1*(x1+x4)+c2*(x2+x3) + I*[ s1*(x1-x4)+s2*(x2-x3)]
!	X2 = x0+c2*(x1+x4)+c1*(x2+x3) + I*[ s2*(x1-x4)-s1*(x2-x3)]
!	X3 = x0+c2*(x1+x4)+c1*(x2+x3) + I*[-s2*(x1-x4)+s1*(x2-x3)]
!	X4 = x0+c1*(x1+x4)+c2*(x2+x3) + I*[-s1*(x1-x4)-s2*(x2-x3)] .
!
!   We refer to the terms not explicitly involving the imaginary constant I
!   as the "cosine part" of the output, and those multiplied by I as the "sine part."
!
!   We seek to form both cosine and sine parts via a minimal set of algebraic operations,
!   with as few multiplies as possible. (Multiplies are as cheap as adds for floating data,
!   but are more expensive for modular data.)
!
!   To form the cosine part we seek a final operation pair of the form
!
!	A+B = 0+c1*(1+4)+c2*(2+3)  (cosine part of X1,X4)
!	A-B = 0+c2*(1+4)+c1*(2+3)  (cosine part of X2,X3).
!
! Thus	A   = 0+(1+2+3+4)*(c1+c2)/2
! and	B   =   (1-2-3+4)*(c1-c2)/2, where we precompute the trigonometric constants (c1+-c2)/2.
!
!   To form the sine part we similarly seek a final operation pair of the form
!
!	C+D = s1*(1-4)+s2*(2-3)  (sine part of X1,-X4)
!	C-D = s2*(1-4)-s1*(2-3)  (sine part of X2,-X3).
!
! Thus	C   = (1-4)*(s1+s2)/2 - (2-3)*(s1-s2)/2
! and	D   = (1-4)*(s1-s2)/2 + (2-3)*(s1+s2)/2,
!
!   but these are not in a desirable form, since they need four multiplies.
!   To reduce the number of multiplies, let
!
!	C' = (1-2+3-4)* s1
!	D' = (  2-3  )*[s1+s2]
!	E' = (1    -4)*[s1-s2],
!
!   where we precompute [s1+-s2] and s1 (i.e. we need just three multiplies for this),
!
! then	C'+D' = s1*(1-4)+s2*(2-3)
! and	C'-E' = s2*(1-4)-s1*(2-3) .
!
!   Totals :  5 nontrivial multiplications (all by real constants, i.e. amounting to 10 FMUL)
!            17 complex additions (amounting to 34 FADD).
*/
	int j,j1,j2;
	static int n5,p1,p2,p3,p4, first_entry=TRUE;
	static double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
		      cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
		      s   =  0.95105651629515357211,	/*  sin(u) */
		      ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
		      ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;

	if(!first_entry && (n/5) != n5)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n5=n/5;

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n5;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-5 pass is here.	*/

      for(j=0; j < n5; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*       gather the needed data (5 64-bit complex, i.e. 10 64-bit reals) and begin the transform...	*/
									/** Quantity				= Nussbaumer's	*/
	  t1 =a[j1    ];			t2 =a[j2    ];			/** x0					x0	*/
	  t3 =a[j1+p1 ]+a[j1+p4 ];	t4 =a[j2+p1 ]+a[j2+p4 ];	/** x1 + x4				t1	*/
	  t5 =a[j1+p2 ]+a[j1+p3 ];	t6 =a[j2+p2 ]+a[j2+p3 ];	/** x2 + x3				t2	*/
	  t7 =a[j1+p2 ]-a[j1+p3 ];	t8 =a[j2+p2 ]-a[j2+p3 ];	/** x2 - x3				-t4	*/
	  t9 =a[j1+p1 ]-a[j1+p4 ];	t10=a[j2+p1 ]-a[j2+p4 ];	/** x1 - x4				t3	*/

/*       ...now complete the first complex datum of the result.	*/

	  rt = t3+t5;				it = t4+t6;			/**    x1+x2+x3+x4			t1+t2=t5	*/
	  t1 = t1+rt;				t2 = t2+it;			/** x0+x1+x2+x3+x4			x0+t5=m0	*/

	  rt = t1+cc1*rt;			it = t2+cc1*it;			/** A   = x0+(x1+x2+x3+x4)*[c1+c2]/2	m0+m1=s1	*/
	  t5 = cc2*(t3-t5);			t6 = cc2*(t4-t6);		/** B   =    (x1-x2-x3+x4)*[c1-c2]/2	cc2*(t1-t2)=m2	*/

	  t3 = rt+t5;				t4 = it+t6;			/** A+B = x0+c1*(x1+x4)+c2*(x2+x3)	s1+m2=s2	*/
	  t5 = rt-t5;				t6 = it-t6;			/** A-B = x0+c2*(x1+x4)+c1*(x2+x3)	s1-m2=s4	*/

	  rt = s  *(t9-t7);			it = s  *(t10-t8);		/** C' = (x1-x2+x3-x4)* s1		s*(t3+t4)=-Im(m3)	*/
	  t7 = ss1* t7;				t8 = ss1* t8;			/** D' = (   x2-x3   )*[s1+s2]		ss1*(-t4)= Im(m4)	*/
	  t9 = ss2* t9;				t10= ss2* t10;			/** E' = (x1      -x4)*[s1-s2],		ss2*( t3)= Im(m5)	*/

	  t7 = rt+t7;				t8 = it+t8;			/** C'+D'   = s1*(x1-x4)+s2*(x2-x3)	-Im(m3)+Im(m4)=-Im(s3)	*/
	  t9 = rt-t9;				t10= it-t10;			/** C'-E'   = s2*(x1-x4)-s1*(x2-x3)	-Im(m3)-Im(m4)=-Im(s5)	*/

/*...Inline multiply of sine part by I into finishing phase...	*/

	  a[j1   ]=t1;				a[j2   ]=t2;			/** X0, first output datum		m0	*/
	  a[j1+p1]=t3-t8;			a[j2+p1]=t4+t7;		/** X1	=				s2+s3	*/
	  a[j1+p2]=t5-t10;			a[j2+p2]=t6+t9;		/** X2	=				s4+s5	*/
	  a[j1+p3]=t5+t10;			a[j2+p3]=t6-t9;		/** X3	=				s4-s5	*/
	  a[j1+p4]=t3+t8;			a[j2+p4]=t4-t7;		/** X4	=				s2-s3	*/
	}
}

/***************/

void radix5_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-5 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix5_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
!
!   See the documentation in radix16_dif_pass for further details on the array indexing.
!
!   See the documentation in radix5_dif_pass1 for notes on the algorithm.
*/
	int j,j1,j2;
	static int n5,p1,p2,p3,p4, first_entry=TRUE;
	static double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
		      cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
		      s   =  0.95105651629515357211,	/*  sin(u) */
		      ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
		      ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;

	if(!first_entry && (n/5) != n5)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n5=n/5;

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n5;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-5 pass is here.	*/

      for(j=0; j < n5; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*       gather the needed data (5 64-bit complex, i.e. 10 64-bit reals) and begin the transform...	*/

	    t1 =a[j1   ];		t2 =a[j2   ];

	    t3 =a[j1+p1];		t4 =a[j2+p1];
	    rt =a[j1+p4];		it =a[j2+p4];
	    t9 =t3 -rt;			t10=t4 -it;
	    t3 =t3 +rt;			t4 =t4 +it;

	    t5 =a[j1+p2];		t6 =a[j2+p2];
	    rt =a[j1+p3];		it =a[j2+p3];
	    t7 =t5 -rt;			t8 =t6 -it;
	    t5 =t5 +rt;			t6 =t6 +it;

/*       ...now complete the result.	*/

	    rt = t3+t5;			it = t4+t6;
	    t1 = t1+rt;			t2 = t2+it;

	    rt = t1+cc1*rt;		it = t2+cc1*it;
	    t5 = cc2*(t3-t5);	t6 = cc2*(t4-t6);

	    t3 = rt+t5;			t4 = it+t6;
	    t5 = rt-t5;			t6 = it-t6;

/*	...only difference from the forward transform is that in this next part, we inline minus signs in front of s,ss1,ss2.	*/

	    rt = s  *(t7-t9);		it = s  *(t8-t10);
	    t7 = ss1* t7;		t8 = ss1* t8;
	    t9 = ss2* t9;		t10= ss2* t10;

	    t7 = rt-t7;			t8 = it-t8;
	    t9 = rt+t9;			t10= it+t10;

	    a[j1    ]=t1;		a[j2    ]=t2;
	    a[j1+p1 ]=t3-t8;	a[j2+p1 ]=t4+t7;
	    a[j1+p2 ]=t5-t10;	a[j2+p2 ]=t6+t9;
	    a[j1+p3 ]=t5+t10;	a[j2+p3 ]=t6-t9;
	    a[j1+p4 ]=t3+t8;	a[j2+p4 ]=t4-t7;
	}
}

