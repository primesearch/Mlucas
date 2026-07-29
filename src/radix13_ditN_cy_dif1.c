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
#include "radix13.h"

/***************/

int radix13_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-13 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-14 complex DIF pass on the data in
!   the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n13,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10,bjmodn11,bjmodn12
		,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double radix_inv, n2inv;
	double cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,cy10,cy11,cy12,temp,frac,scale
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i;
	double maxerr = 0.0;
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

/*...change n13 and n_div_wt to non-static to work around a gcc compiler bug. */
	n13   = n/13;
	n_div_nwt = n13 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n13)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/13 in radix13_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)13));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
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

		bjmodnini=0;
		for(j=0; j < n13; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-13 final DIT pass is here.	*/

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

	col=0;
	co2=(n >> nwt_bits)-1+13;
	co3=co2-13;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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

		/* Do the final radix-13 DIT */
		RADIX_13_DFT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12]
					,aj1p0r,aj1p0i,aj1p12r,aj1p12i,aj1p11r,aj1p11i,aj1p10r,aj1p10i,aj1p9r,aj1p9i,aj1p8r,aj1p8i,aj1p7r,aj1p7i,aj1p6r,aj1p6i,aj1p5r,aj1p5i,aj1p4r,aj1p4i,aj1p3r,aj1p3i,aj1p2r,aj1p2i,aj1p1r,aj1p1i);

		/*..Now do the carries. Since the outputs would
		normally be getting dispatched to 13 separate blocks of the A-array, we need 13 separate carries.	*/

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

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		/*...The radix-13 DIF pass is here:	*/
		RADIX_13_DFT(aj1p0r,aj1p0i,aj1p1r,aj1p1i,aj1p2r,aj1p2i,aj1p3r,aj1p3i,aj1p4r,aj1p4i,aj1p5r,aj1p5i,aj1p6r,aj1p6i,aj1p7r,aj1p7i,aj1p8r,aj1p8i,aj1p9r,aj1p9i,aj1p10r,aj1p10i,aj1p11r,aj1p11i,aj1p12r,aj1p12i
					,a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12]);
		}

		jstart += nwt;
		jhi    += nwt;
		col += 13;
		co3 -= 13;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-13 forward DIF FFT of the first block of 13 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 13 outputs of (1);
!   (3) Reweight and perform a radix-13 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 13 elements and repeat (1-4).
*/

	temp= cy12;
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
	cy0 = temp;

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
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9)+fabs(cy10)+fabs(cy11)+fabs(cy12) != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix13_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}
	*fracmax = maxerr;
	return(0);
}


/***************/

void radix13_dif_pass1(double a[], int n)
{
/*
!
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-13 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n13,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, first_entry=TRUE;
#if 0	/* Keep these variables and the original code sequence around just for reference: */
	double xr1, xi1, xr2, xi2, xr3, xi3, xr4, xi4, xr5, xi5, xr6, xi6, xr7, xi7;
	double yr1, yi1, yr2, yi2, yr3, yi3, yr4, yi4, yr5, yi5, yr6, yi6;
	double pr1, pi1, pr2, pi2, pr3, pi3, pr4, pi4, pr5, pi5, pr6, pi6, pr7, pi7;
	double ur1, ui1, ur2, ui2, ur3, ui3, ur4, ui4, ur5, ui5, ur6, ui6, ur7, ui7;
	double sr1, si1, sr2, si2, sr3, si3, sr4, si4, sr5, si5, sr6, si6, sr7, si7;
	double qr1, qi1, qr2, qi2, qr3, qi3, qr4, qi4, qr5, qi5, qr6, qi6, qr34, qi34, qr56, qi56;
	double vr1, vi1, vr2, vi2, vr3, vi3, vr4, vi4, vr5, vi5, vr6, vi6, vr34, vi34, vr56, vi56;
	double tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4, tr5, ti5, tr6, ti6;
	const double dc1 = -0.66277335223429382656331008444690382,
		dc2  =  0.73124599097534822519618254560377760,
		dc3  =  1.0070740657275332544937477077369340,
		dc4  = -0.30816846519175820059367219184688543,
		dc5  =  0.81698338691215549726750306085822509,
		dc6  =  0.22376403323791637458136993583748652,
		ds1  =  0.57514072947400312136838554745545335,
		ds2  = -0.17413860115213590500566079492926474,
		ds3  = -0.33582506518644535421963182119524142,
		ds4  =  4.5240494294812713569277280991401412E-0002,
		ds5  =  1.1543953381323634420147226757584967,
		ds6  = -8.7981928766792081008399945211312915E-0002,
		ds7  =  0.90655220171271016880349079977456841,
		ds8  =  1.1971367726043428094538453399784083,
		ds9  = -0.24784313641965327321123187598392850,
		ds10 = -0.86131170741789745523421351878316690,
		ds11 = -4.2741434471979367439122664219911502E-0002;
#endif
	if(!first_entry && (n/13) != n13)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n13=n/13;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
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
	}

/*...The radix-13 pass is here.	*/

	for(j=0; j < n13; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

#if 0
	/* gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		xr1 = a[j1	 ];					xi1 = a[j2	 ];
		xr5 = a[j1+p8] + a[j1+p5 ];		xi5 = a[j2+p8] + a[j2+p5 ];
		xr7 = a[j1+p6] + a[j1+p7 ];		xi7 = a[j2+p6] + a[j2+p7 ];
		xr4 = a[j1+p4] + a[j1+p9 ];		xi4 = a[j2+p4] + a[j2+p9 ];
		xr6 = a[j1+p3] + a[j1+p10];		xi6 = a[j2+p3] + a[j2+p10];
		xr3 = a[j1+p2] + a[j1+p11];		xi3 = a[j2+p2] + a[j2+p11];
		xr2 = a[j1+p1] + a[j1+p12];		xi2 = a[j2+p1] + a[j2+p12];
		yr1 = a[j1+p1] - a[j1+p12];		yi1 = a[j2+p1] - a[j2+p12];
		yr2 = a[j1+p2] - a[j1+p11];		yi2 = a[j2+p2] - a[j2+p11];
		yr5 = a[j1+p3] - a[j1+p10];		yi5 = a[j2+p3] - a[j2+p10];
		yr3 = a[j1+p4] - a[j1+p9 ];		yi3 = a[j2+p4] - a[j2+p9 ];
		yr6 = a[j1+p6] - a[j1+p7 ];		yi6 = a[j2+p6] - a[j2+p7 ];
		yr4 = a[j1+p8] - a[j1+p5 ];		yi4 = a[j2+p8] - a[j2+p5 ];

		pr2 = xr2+xr5;					pi2 = xi2+xi5;
		pr5 = xr2-xr5;					pi5 = xi2-xi5;
		pr3 = xr3+xr6;					pi3 = xi3+xi6;
		pr6 = xr3-xr6;					pi6 = xi3-xi6;
		pr4 = xr4+xr7;					pi4 = xi4+xi7;
		pr7 = xr4-xr7;					pi7 = xi4-xi7;
		pr1 = pr2+pr3+pr4;				pi1 = pi2+pi3+pi4;
		sr1 = xr1+pr1;					si1 = xi1+pi1;
		ur1 = xr1+dc1*pr1;				ui1 = xi1+dc1*pi1;
		ur2 = ur1+dc2*pr3+dc3*pr4;		ui2 = ui1+dc2*pi3+dc3*pi4;
		ur3 = ur1+dc2*pr4+dc3*pr2;		ui3 = ui1+dc2*pi4+dc3*pi2;
		ur4 = ur1+dc2*pr2+dc3*pr3;		ui4 = ui1+dc2*pi2+dc3*pi3;
		ur5 = dc4*pr5+dc5*pr6+dc6*pr7;	ui5 = dc4*pi5+dc5*pi6+dc6*pi7;
		ur6 = dc4*pr6+dc5*pr7-dc6*pr5;	ui6 = dc4*pi6+dc5*pi7-dc6*pi5;
		ur7 = dc4*pr7-dc5*pr5-dc6*pr6;	ui7 = dc4*pi7-dc5*pi5-dc6*pi6;
		sr2 = ur2+ur5;					si2 = ui2+ui5;
		sr5 = ur2-ur5;					si5 = ui2-ui5;
		sr3 = ur3+ur6;					si3 = ui3+ui6;
		sr6 = ur3-ur6;					si6 = ui3-ui6;
		sr4 = ur4+ur7;					si4 = ui4+ui7;
		sr7 = ur4-ur7;					si7 = ui4-ui7;

		qr1 = yr1-yr3+yr5;				qi1 = yi1-yi3+yi5;
		qr3 = yr1+yr3;					qi3 = yi1+yi3;
		qr5 = yr3+yr5;		 			qi5 = yi3+yi5;
		qr2 = yr2-yr4+yr6;				qi2 = yi2-yi4+yi6;
		qr4 = yr2+yr4;					qi4 = yi2+yi4;
		qr6 = yr4+yr6;		 			qi6 = yi4+yi6;
		vr1 = ds1*qr1+ds2*qr2;			vi1 = ds1*qi1+ds2*qi2;
		vr2 = ds1*qr2-ds2*qr1;			vi2 = ds1*qi2-ds2*qi1;
		qr34= qr3+qr4;		 			qi34= qi3+qi4;
		qr56= qr5+qr6;					qi56= qi5+qi6;
		vr34= ds3*qr34-ds6*qr56;		vi34= ds3*qi34-ds6*qi56;
		vr3 = vr34+ds4*qr4-ds7*qr6;		vi3 = vi34+ds4*qi4-ds7*qi6;
		vr4 = vr34+ds5*qr3-ds8*qr5;		vi4 = vi34+ds5*qi3-ds8*qi5;
		vr56= ds9*qr34-ds3*qr56;		vi56= ds9*qi34-ds3*qi56;
		vr5 = vr56+ds10*qr4-ds4*qr6;	vi5 = vi56+ds10*qi4-ds4*qi6;
		vr6 = vr56+ds11*qr3-ds5*qr5;	vi6 = vi56+ds11*qi3-ds5*qi5;
		tr1 = vr1+vr3;					ti1 = vi1+vi3;
		tr3 = vr3-vr1-vr5;				ti3 = vi3-vi1-vi5;
		tr5 = vr1-vr5;					ti5 = vi1-vi5;
		tr2 = vr2+vr4;					ti2 = vi2+vi4;
		tr4 = vr4-vr2-vr6;				ti4 = vi4-vi2-vi6;
		tr6 = vr2-vr6;					ti6 = vi2-vi6;

		a[j1	] = sr1;				a[j2	] = si1;
		a[j1+p1 ] = sr7-ti6;	 		a[j2+p1 ] = si7+tr6;
		a[j1+p2 ] = sr6-ti5;	 		a[j2+p2 ] = si6+tr5;
		a[j1+p3 ] = sr3-ti2;	 		a[j2+p3 ] = si3+tr2;
		a[j1+p4 ] = sr5-ti4;	 		a[j2+p4 ] = si5+tr4;
		a[j1+p5 ] = sr4+ti3;			a[j2+p5 ] = si4-tr3;
		a[j1+p6 ] = sr2-ti1;			a[j2+p6 ] = si2+tr1;
		a[j1+p7 ] = sr2+ti1;			a[j2+p7 ] = si2-tr1;
		a[j1+p8 ] = sr4-ti3;			a[j2+p8 ] = si4+tr3;
		a[j1+p9 ] = sr5+ti4;			a[j2+p9 ] = si5-tr4;
		a[j1+p10] = sr3+ti2;			a[j2+p10] = si3-tr2;
		a[j1+p11] = sr6+ti5;			a[j2+p11] = si6-tr5;
		a[j1+p12] = sr7+ti6;			a[j2+p12] = si7-tr6;
	/* Cost: 164 ADD, 64 MUL (compare to 188 Add, 40 MUL for Johnson/Burrus [http://cnx.org/content/m17413/latest/?collection=col10569/latest] */

#elif 0
		RADIX_13_DFT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12]
					,a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12]);
#else
		RADIX_13_XYZ(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12]
					,a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12]);
#endif
	}
}

/* DIT - Since the radix is prime - differs from DIF only in having the + and - of the final output-combinations swapped.
	We also use the DIT to test out an in-place scheme prototyping the kind of memory-location reuse we need in assembler,
	where we pay the price of some extra arithmetic in order to do things in-register, e.g. we replace the sequence
		x = a+b, y = a-b
	with
		a = a+b, b *= 2, b = a-b, which costs an extra MUL but is in-place.
*/
void radix13_dit_pass1(double a[], int n)
{
/*
!
!...Acronym: DIF = Decimation In Time
!
!...Subroutine to perform an initial radix-13 complex DFI FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n13,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, first_entry=TRUE;

	if(!first_entry && (n/13) != n13)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n13=n/13;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
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
	}

/*...The radix-13 pass is here.	*/

	for(j=0; j < n13; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

#if 0
	/* gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		xr0 = a[j1	 ];					xi0 = a[j2	 ];
		xr7 = a[j1+p6] + a[j1+p7 ];		xi7 = a[j2+p6] + a[j2+p7 ];
		xr5 = a[j1+p5] + a[j1+p8 ];		xi5 = a[j2+p5] + a[j2+p8 ];
		xr4 = a[j1+p4] + a[j1+p9 ];		xi4 = a[j2+p4] + a[j2+p9 ];
		xr6 = a[j1+p3] + a[j1+p10];		xi6 = a[j2+p3] + a[j2+p10];
		xr3 = a[j1+p2] + a[j1+p11];		xi3 = a[j2+p2] + a[j2+p11];
		xr1 = a[j1+p1] + a[j1+p12];		xi1 = a[j2+p1] + a[j2+p12];

		yr1 = a[j1+p1] - a[j1+p12];		yi1 = a[j2+p1] - a[j2+p12];
		yr2 = a[j1+p2] - a[j1+p11];		yi2 = a[j2+p2] - a[j2+p11];
		yr5 = a[j1+p3] - a[j1+p10];		yi5 = a[j2+p3] - a[j2+p10];
		yr3 = a[j1+p4] - a[j1+p9 ];		yi3 = a[j2+p4] - a[j2+p9 ];
		yr4 = a[j1+p5] - a[j1+p8 ];		yi4 = a[j2+p5] - a[j2+p8 ];
		yr6 = a[j1+p6] - a[j1+p7 ];		yi6 = a[j2+p6] - a[j2+p7 ];

// x-terms need 8 registers for each side:
		xr2 = xr1+xr5;					xi2 = xi1+xi5;
		xr5 = xr1-xr5;					xi5 = xi1-xi5;
		xr3 = xr3+xr6;	xr6 *= 2.0;		xi3 = xi3+xi6;	xi6 *= 2.0;
		xr6 = xr3-xr6;					xi6 = xi3-xi6;
		xr4 = xr4+xr7;	xr7 *= 2.0;		xi4 = xi4+xi7;	xi7 *= 2.0;
		xr7 = xr4-xr7;					xi7 = xi4-xi7;
		xr1 = xr2+xr3+xr4;				xi1 = xi2+xi3+xi4;
		xr0 = xr0+xr1;					xi0 = xi0+xi1;
		a[j1	] = xr0;				a[j2	] = xi0;
		xr0 +=    DC1*xr1;				xi0 +=    DC1*xi1;
		xr1 = xr0+DC2*xr2+DC3*xr3;		xi1 = xi0+DC2*xi2+DC3*xi3;
		xr2 = xr0+DC2*xr4+DC3*xr2;		xi2 = xi0+DC2*xi4+DC3*xi2;
		xr3 = xr0+DC2*xr3+DC3*xr4;		xi3 = xi0+DC2*xi3+DC3*xi4;
		xr0 = DC4*xr5+DC5*xr6+DC6*xr7;	xi0 = DC4*xi5+DC5*xi6+DC6*xi7;
		xr4 = DC4*xr6+DC5*xr7-DC6*xr5;	xi4 = DC4*xi6+DC5*xi7-DC6*xi5;
		xr7 = DC4*xr7-DC5*xr5-DC6*xr6;	xi7 = DC4*xi7-DC5*xi5-DC6*xi6;
		xr5 = xr3+xr0;					xi5 = xi3+xi0;
		xr3 = xr3-xr0;					xi3 = xi3-xi0;
		xr6 = xr2+xr4;					xi6 = xi2+xi4;
		xr2 = xr2-xr4;					xi2 = xi2-xi4;
		xr4 = xr1+xr7;					xi4 = xi1+xi7;
		xr7 = xr1-xr7;					xi7 = xi1-xi7;
// x0,1 free...now do y-terms:
		xr0 = yr1-yr3+yr5;				xi0 = yi1-yi3+yi5;
		yr1 = yr1+yr3;					yi1 = yi1+yi3;
		yr5 = yr3+yr5;		 			yi5 = yi3+yi5;
		yr3 = yr2+yr4+yr6;				yi3 = yi2+yi4+yi6;
		yr2 = yr2-yr4;					yi2 = yi2-yi4;
		yr4 = yr6-yr4;		 			yi4 = yi6-yi4;
		xr1 = DS1*xr0+DS2*yr3;			xi1 = DS1*xi0+DS2*yi3;
		ttr = DS1*yr3-DS2*xr0;			tti = DS1*yi3-DS2*xi0;
		yr3 = yr1+yr2;		 			yi3 = yi1+yi2;
		yr6 = yr5+yr4;					yi6 = yi5+yi4;
		xr0 = DS3*yr3-DS6*yr6;			xi0 = DS3*yi3-DS6*yi6;
		yr3 = DS9*yr3-DS3*yr6;			yi3 = DS9*yi3-DS3*yi6;
		yr6 = xr0+DS4*yr2-DS7*yr4;		yi6 = xi0+DS4*yi2-DS7*yi4;
		yr2 = yr3+DSa*yr2-DS4*yr4;		yi2 = yi3+DSa*yi2-DS4*yi4;
		yr4 = xr0+DS5*yr1-DS8*yr5;		yi4 = xi0+DS5*yi1-DS8*yi5;
		yr3 = yr3+DSb*yr1-DS5*yr5;		yi3 = yi3+DSb*yi1-DS5*yi5;
		yr5 = xr1+yr6;					yi5 = xi1+yi6;
		yr6 = yr6-xr1-yr2;				yi6 = yi6-xi1-yi2;
		yr2 = xr1-yr2;					yi2 = xi1-yi2;
		yr1 = ttr+yr4;					yi1 = tti+yi4;
		yr4 = yr4-ttr-yr3;				yi4 = yi4-tti-yi3;
		yr3 = ttr-yr3;					yi3 = tti-yi3;
// In ASM, do xr and yi-terms and combine to get real parts of outputs,
// then do xi and yr-terms and combine to get imaginary parts:
		a[j1+p1 ] = xr7+yi3;	 		a[j2+p1 ] = xi7-yr3;
		a[j1+p12] = xr7-yi3;			a[j2+p12] = xi7+yr3;
		a[j1+p2 ] = xr2+yi2;	 		a[j2+p2 ] = xi2-yr2;
		a[j1+p11] = xr2-yi2;			a[j2+p11] = xi2+yr2;
		a[j1+p3 ] = xr6+yi1;	 		a[j2+p3 ] = xi6-yr1;
		a[j1+p10] = xr6-yi1;			a[j2+p10] = xi6+yr1;
		a[j1+p4 ] = xr3+yi4;	 		a[j2+p4 ] = xi3-yr4;
		a[j1+p9 ] = xr3-yi4;			a[j2+p9 ] = xi3+yr4;
		a[j1+p5 ] = xr4-yi6;			a[j2+p5 ] = xi4+yr6;
		a[j1+p8 ] = xr4+yi6;			a[j2+p8 ] = xi4-yr6;
		a[j1+p6 ] = xr5+yi5;			a[j2+p6 ] = xi5-yr5;
		a[j1+p7 ] = xr5-yi5;			a[j2+p7 ] = xi5+yr5;
	/* Cost: 164 ADD, 64 MUL (compare to 188 Add, 40 MUL for Johnson/Burrus [http://cnx.org/content/m17413/latest/?collection=col10569/latest] */

#elif 0
		RADIX_13_DFT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12]
					,a[j1],a[j2],a[j1+p12],a[j2+p12],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],a[j1+p9],a[j2+p9],a[j1+p8],a[j2+p8],a[j1+p7],a[j2+p7],a[j1+p6],a[j2+p6],a[j1+p5],a[j2+p5],a[j1+p4],a[j2+p4],a[j1+p3],a[j2+p3],a[j1+p2],a[j2+p2],a[j1+p1],a[j2+p1]);
#else
		RADIX_13_XYZ(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12]
					,a[j1],a[j2],a[j1+p12],a[j2+p12],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],a[j1+p9],a[j2+p9],a[j1+p8],a[j2+p8],a[j1+p7],a[j2+p7],a[j1+p6],a[j2+p6],a[j1+p5],a[j2+p5],a[j1+p4],a[j2+p4],a[j1+p3],a[j2+p3],a[j1+p2],a[j2+p2],a[j1+p1],a[j2+p1]);
#endif
	}
}

