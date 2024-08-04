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

int radix6_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-6 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-6 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n6,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double s = 0.86602540378443864675, c3m1 = -1.5, radix_inv,n2inv;	/* exp[i*(twopi/6)], cos(twopi/3)-1	*/
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i
		,cy0,cy1,cy2,cy3,cy4,cy5,temp,frac,scale;
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

/*...change n6 and n_div_wt to non-static to work around a gcc compiler bug. */
	n6   = n/6;
	n_div_nwt = n6 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n6)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/6 in radix6_ditN_cy_dif1.\n",iter);
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
	  radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)6));
	  n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

	  bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
	  sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n6;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	  p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );

	  bjmodnini=0;
	  for(j=0; j < n6; j++)
	  {
	    bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
	  }
	}

/*...The radix-6 final DIT pass is here.	*/

	cy0 = 0;		/* init carries	*/
	cy1 = 0;
	cy2 = 0;
	cy3 = 0;
	cy4 = 0;
	cy5 = 0;

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

	col=0;
	co2=(n >> nwt_bits)-1+6;
	co3=co2-6;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

	    t1 =a[j1   ];			t2 =a[j2   ];
	    rt =a[j1+p1];			it =a[j2+p1];
	    t3 =t1 -rt;				t4 =t2 -it;
	    t1 =t1 +rt;				t2 =t2 +it;

	    t5 =a[j1+p3];			t6 =a[j2+p3];
	    rt =a[j1+p2];			it =a[j2+p2];
	    t7 =t5 -rt;				t8 =t6 -it;
	    t5 =t5 +rt;				t6 =t6 +it;

	    t9 =a[j1+p4];			t10=a[j2+p4];
	    rt =a[j1+p5];			it =a[j2+p5];
	    t11=t9 -rt;				t12=t10-it;
	    t9 =t9 +rt;				t10=t10+it;

	    rt =t9;					it =t10;
	    t9 =t5 -rt;				t10=t6 -it;
	    t5 =t5 +rt;				t6 =t6 +it;
	    t1 =t1 +t5;				t2 =t2 +t6;
	    aj1p0r =t1;				aj1p0i =t2;
	    t5 =t1+c3m1*t5;			t6 =t2+c3m1*t6;
	    t9 =      s*t9;			t10=      s*t10;
	    aj1p2r=t5+t10;			aj1p2i=t6-t9;
	    aj1p4r=t5-t10;			aj1p4i=t6+t9;
/* ...and to order these three outputs as 3,5,1 rather than 1,3,5. */
	    rt =t11;				it =t12;
	    t11=t7 -rt;				t12=t8 -it;
	    t7 =t7 +rt;				t8 =t8 +it;
	    t3 =t3 +t7;				t4 =t4 +t8;
	    aj1p3r =t3;				aj1p3i =t4;
	    t7 =t3+c3m1*t7;			t8 =t4+c3m1*t8;
	    t11=      s*t11;			t12=      s*t12;
	    aj1p5r=t7+t12;			aj1p5i=t8-t11;
	    aj1p1r=t7-t12;			aj1p1i=t8+t11;

/*...and combine those to complete the radix-6 transform and do the carries. Since the outputs would
    normally be getting dispatched to 6 separate blocks of the A-array, we need 6 separate carries.	*/

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

	    i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
	    co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			   and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-6 DIF pass is here:	*/

/*       gather the needed data (6 64-bit complex, i.e. 8 64-bit reals) and do two radix-3 transforms...	*/
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif
/* Twiddleless version requires us to swap inputs x2 <-> x4 and x1 <-> x3... */
	  t1 =aj1p0r;			t2 =aj1p0i;
	  t3 =aj1p4r+aj1p2r;	t4 =aj1p4i+aj1p2i;
	  t5 =aj1p4r-aj1p2r;	t6 =aj1p4i-aj1p2i;
	  t1 =t1+t3;			t2 =t2+t4;
	  t3=t1+c3m1*t3;		t4=t2+c3m1*t4;
	  t5=      s*t5;		t6=      s*t6;
	  rt=t5;				it=t6;
	  t5=t3+it;				t6=t4-rt;
	  t3=t3-it;				t4=t4+rt;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
	  t7 =aj1p3r;			t8 =aj1p3i;
	  t9 =aj1p1r+aj1p5r;	t10=aj1p1i+aj1p5i;
	  t11=aj1p1r-aj1p5r;	t12=aj1p1i-aj1p5i;
	  t7 =t7+t9;			t8 =t8+t10;
	  t9 =t7+c3m1*t9;		t10=t8+c3m1*t10;
	  t11=      s*t11;		t12=      s*t12;
	  rt=t11;				it=t12;
	  t11=t9+it;			t12=t10-rt;
	  t9 =t9-it;			t10=t10+rt;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif
/* ...and to order the outputs as 0,1,4,5,3,2, rtaher than 0,1,2,3,4,5. */
	  a[j1   ]=t1+t7;		a[j2   ]=t2+t8;
	  a[j1+p1]=t1-t7;		a[j2+p1]=t2-t8;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
	  a[j1+p4]=t3 +t9;		a[j2+p4]=t4+t10;
	  a[j1+p5]=t3 -t9;		a[j2+p5]=t4-t10;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
	  a[j1+p3]=t5 +t11;		a[j2+p3]=t6+t12;
	  a[j1+p2]=t5 -t11;		a[j2+p2]=t6-t12;
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
	  }

	  jstart += nwt;
	  jhi    += nwt;
	  col += 6;
	  co3 -= 6;

	}

	if(root_incr==0)break;

/**** Wraparound carry cleanup loop is here: ***

!   The cleanup carries from the end of each length-N/6 block into the begining of the next
!   can all be neatly processed as follows:
!
!   (1) Invert the radix-6 forward DIF FFT of the first block of 6 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 6 outputs of (1);
!   (3) Reweight and perform a radix-6 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next three elements and repeat (1-4).
*/
	t1  = cy5;
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
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5) != 0.0)
	{
	    sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix6_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
	    mlucas_fprint(cbuf,INTERACT);
	    err=ERR_CARRY;
	    return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix6_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-6 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   Uses an optimized radix-3 transform a la Nussbaumer (2nd ed., p.145).
*/
	int j,j1,j2;
	static int n6,p1,p2,p3,p4,p5, first_entry=TRUE;
	static double s = 0.86602540378443864675, c3m1 = -1.5;	/* exp[i*(twopi/6)], cos(twopi/3)-1	*/
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;

	if(!first_entry && (n/6) != n6)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n6=n/6;

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n6;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	  p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-6 pass is here.	*/

      for(j=0; j < n6; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*       gather the needed data (6 64-bit complex, i.e. 8 64-bit reals) and do two radix-3 transforms...	*/

/* Twiddleless version requires us to swap inputs x2 <-> x4 and x1 <-> x3...
indices  0, 1, 2, 3, 4, 5
      -> 0,-3, 4, 1, 2,-1
      == 0, 3, 4, 1, 2, 5 modulo 6.
I.e. start out with first triplet of indices {0,2,4}, permute those according to
{0,2,4}*5%6 = {0,4,2}, then each is head of a length-2 list of indices with decrement 3.
*/
	  t1 =a[j1   ];				t2 =a[j2   ];
	  t3 =a[j1+p4]+a[j1+p2 ];	t4 =a[j2+p4]+a[j2+p2];
	  t5 =a[j1+p4]-a[j1+p2 ];	t6 =a[j2+p4]-a[j2+p2];
	  t1 =t1+t3;				t2 =t2+t4;
	  t3=t1+c3m1*t3;			t4=t2+c3m1*t4;
	  t5=      s*t5;			t6=      s*t6;
	  rt=t5;					it=t6;
	  t5=t3+it;					t6=t4-rt;
	  t3=t3-it;					t4=t4+rt;

	  t7 =a[j1+p3];				t8 =a[j2+p3];
	  t9 =a[j1+p1]+a[j1+p5 ];	t10=a[j2+p1]+a[j2+p5];
	  t11=a[j1+p1]-a[j1+p5 ];	t12=a[j2+p1]-a[j2+p5];
	  t7 =t7+t9;				t8 =t8+t10;
	  t9 =t7+c3m1*t9;			t10=t8+c3m1*t10;
	  t11=      s*t11;			t12=      s*t12;
	  rt=t11;					it=t12;
	  t11=t9+it;				t12=t10-rt;
	  t9 =t9-it;				t10=t10+rt;

/* ...and to order the outputs as 0,1,4,5,3,2, rather than 0,1,2,3,4,5. */

	  a[j1   ]=t1+t7;			a[j2   ]=t2+t8;
	  a[j1+p1]=t1-t7;			a[j2+p1]=t2-t8;

	  a[j1+p4]=t3 +t9;			a[j2+p4]=t4+t10;
	  a[j1+p5]=t3 -t9;			a[j2+p5]=t4-t10;

	  a[j1+p3]=t5 +t11;			a[j2+p3]=t6+t12;
	  a[j1+p2]=t5 -t11;			a[j2+p2]=t6-t12;
								/* Totals: 2*12+3*4 = 36 FADD, 2*4+3*0 = 8 FMUL.	*/
	}
}

/**************/

void radix6_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-6 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix6_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
!
!   See the documentation in radix16_dif_pass for further details on the array indexing,
!   and radix6_dif_pass1 for details on the algorithm.
*/
	int j,j1,j2;
	static int n6,p1,p2,p3,p4,p5, first_entry=TRUE;
	static double s = 0.86602540378443864675, c3m1 = -1.5;	/* exp[i*(twopi/6)], cos(twopi/3)-1	*/
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;

	if(!first_entry && (n/6) != n6)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n6=n/6;

/*   constant index offsets for array load/stores are here.	*/

	  p1 = n6;
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;

	  p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	  p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	  p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	  p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	  p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-6 pass is here.	*/

      for(j=0; j < n6; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*       gather the needed data (6 64-bit complex, i.e. 8 64-bit reals) and do two radix-3 transforms...	*/

/* Twiddleless version requires us to swap inputs x1 <-> x4:
indices  0, 1, 2, 3, 4, 5
      -> 0,-2,-4, 3, 1,-1
      == 0, 4, 2, 3, 1, 5 modulo 6.
I.e. start out with first pair of indices {0,3}, permute those according to
{0,3}*5%6 = {0,3}, then each is head of a length-3 list of indices with decrement 2.

Remember, inputs to DIT are bit-reversed, so a[0,1,2,3,4,5] contain x[0,3,1,4,2,5], so swapping x1 and x4 means swapping a[2] and a[3].
*/
	    t1 =a[j1   ];			t2 =a[j2   ];
	    rt =a[j1+p1];			it =a[j2+p1];
	    t3 =t1 -rt;				t4 =t2 -it;
	    t1 =t1 +rt;				t2 =t2 +it;

	    t5 =a[j1+p3];			t6 =a[j2+p3];
	    rt =a[j1+p2];			it =a[j2+p2];
	    t7 =t5 -rt;				t8 =t6 -it;
	    t5 =t5 +rt;				t6 =t6 +it;

	    t9 =a[j1+p4];			t10=a[j2+p4];
	    rt =a[j1+p5];			it =a[j2+p5];
	    t11=t9 -rt;				t12=t10-it;
	    t9 =t9 +rt;				t10=t10+it;

/* ...and to order the sencond set of three outputs as 3,5,1 rather than 1,3,5. */

	    rt =t9;					it =t10;
	    t9 =t5 -rt;				t10=t6 -it;
	    t5 =t5 +rt;				t6 =t6 +it;
	    t1 =t1+t5;				t2 =t2+t6;
	    t5 =t1+c3m1*t5;			t6 =t2+c3m1*t6;
	    t9 =      s*t9;			t10=      s*t10;

	    a[j1   ]=t1;			a[j2   ]=t2;
	    a[j1+p2]=t5+t10;		a[j2+p2]=t6-t9;
	    a[j1+p4]=t5-t10;		a[j2+p4]=t6+t9;

	    rt =t11;				it =t12;
	    t11=t7 -rt;				t12=t8 -it;
	    t7 =t7 +rt;				t8 =t8 +it;
	    t3 =t3 +t7;				t4 =t4 +t8;
	    t7 =t3+c3m1*t7;			t8 =t4+c3m1*t8;
	    t11=      s*t11;		t12=      s*t12;

	    a[j1+p3]=t3;			a[j2+p3]=t4;
	    a[j1+p5]=t7+t12;		a[j2+p5]=t8-t11;
	    a[j1+p1]=t7-t12;		a[j2+p1]=t8+t11;
	}
}

