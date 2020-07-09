/*******************************************************************************
*                                                                              *
*   (C) 1997-2019 by Ernst W. Mayer.                                           *
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

// This main loop is same for un-and-multithreaded, so stick into a header file
// (can't use a macro because of the #if-enclosed stuff).

for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
{
	for(j = jstart; j < jhi; j += stride)
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-63 DIT pass is here:	*/

	//...gather the needed data (63 64-bit complex, i.e. 126 64-bit reals) and do 7 radix-9 transforms:
		tptr = t; iptr = dit_iperm;
		for(l = 0; l < 7; l++) {
			k0 = p[*iptr]; k1 = p[*(iptr+1)]; k2 = p[*(iptr+2)]; k3 = p[*(iptr+3)]; k4 = p[*(iptr+4)]; k5 = p[*(iptr+5)]; k6 = p[*(iptr+6)]; k7 = p[*(iptr+7)]; k8 = p[*(iptr+8)];
			RADIX_09_DIT(
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],a[j1+k7],a[j2+k7],a[j1+k8],a[j2+k8],
				tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im,
				rt,it,re
			);	tptr += 9; iptr += 9;
		}
	//...and now do 9 radix-7 transforms:
		tptr = t; iptr = dit_operm;
		for(l = 0; l < 9; l++) {
			k0 = p[*iptr]; k1 = p[*(iptr+1)]; k2 = p[*(iptr+2)]; k3 = p[*(iptr+3)]; k4 = p[*(iptr+4)]; k5 = p[*(iptr+5)]; k6 = p[*(iptr+6)];
			RADIX_07_DFT(
				tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im,
				t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],
				uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
			);	tptr++; iptr += 7;
		}

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 63 separate blocks of the A-array, we need 63 separate carries.	*/

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

		/*...set0 is slightly different from others; divide work into blocks of RADIX/4 macro calls, 1st set of which gets pulled out of loop: */
// Apr 2014: Fermat-mod works fine, but mers-mod barfs immediately with what looks like a bad a0 value,
// div-by-n/2 should give 16, but instead see
//	iter 1, full = 1, a0in =   15.492078993055555
//	iter 1, full = 1, a0out =   13.000000000000000
//	Iter = 1, maxerr =    0.492078993055555
//if(!j)printf("iter %d, full = %d, a0in = %20.15f\n",iter,full_pass,a[0]/(n>>1));
		l = 0; addr = cy_r; itmp = bjmodn;
		jt = j1; jp = j2;
	   cmplx_carry_norm_errcheck0(a[j1   ],a[j2   ],*addr,*itmp,0,prp_mult); ++l; ++addr; ++itmp;
		// Next 15 quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < 16; ntmp++) {
			cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			jt = j1 + p[ntmp<<2]; jp = j2 + p[ntmp<<2];
			cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		}
		// Cleanup of final 2 sets of carries:
		cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
//if(!j)printf("iter %d, full = %d, a0out = %20.15f\n",iter,full_pass,a[0]);
		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
	}
	else	/* MODULUS_TYPE_FERMAT */
	{
		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		ntmp = 0; addr = cy_r; addi = cy_i; ic = 0;	// ic = idx into icycle mini-array, gets incremented (mod ODD_RADIX) between macro calls
		jt = j1; jp = j2;
		fermat_carry_norm_errcheckB(a[jt   ],a[jp   ],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; ++ic;
		fermat_carry_norm_errcheckB(a[jt+p1],a[jp+p1],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; ++ic;
		fermat_carry_norm_errcheckB(a[jt+p2],a[jp+p2],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; ++ic;
		for(m = 1; m < 16; m++) {
			fermat_carry_norm_errcheckB(a[jt+p3],a[jp+p3],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; ++ic;
			jt = j1 + p[m<<2]; jp = j2 + p[m<<2];
			fermat_carry_norm_errcheckB(a[jt   ],a[jp   ],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; ++ic;
			fermat_carry_norm_errcheckB(a[jt+p1],a[jp+p1],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; ++ic;
			fermat_carry_norm_errcheckB(a[jt+p2],a[jp+p2],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; ++ic;
		}
		for(l = 0; l < ODD_RADIX; l++) {
			icycle[l] += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
			icycle[l] += ( (-(int)((uint32)icycle[l] >> 31)) & nwt);
		}
	}	/* if(MODULUS_TYPE == ...) */

	/*...The radix-63 DIF pass is here:	*/

	//...gather the needed data (63 64-bit complex, i.e. 126 64-bit reals) and do 9 radix-7 transforms:
		tptr = t; iptr = dif_iperm;
		for(l = 0; l < 9; l++) {
			k0 = p[*iptr]; k1 = p[*(iptr+1)]; k2 = p[*(iptr+2)]; k3 = p[*(iptr+3)]; k4 = p[*(iptr+4)]; k5 = p[*(iptr+5)]; k6 = p[*(iptr+6)];
			RADIX_07_DFT(
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],
				t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
				tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im,
				uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
			);	tptr++; iptr += 7;
		}
	//...and now do 7 radix-9 transforms:
		tptr = t; iptr = dif_operm;
		for(l = 0; l < 7; l++) {
			k0 = p[*iptr]; k1 = p[*(iptr+1)]; k2 = p[*(iptr+2)]; k3 = p[*(iptr+3)]; k4 = p[*(iptr+4)]; k5 = p[*(iptr+5)]; k6 = p[*(iptr+6)]; k7 = p[*(iptr+7)]; k8 = p[*(iptr+8)];
			RADIX_09_DIF(
				tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im,
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],a[j1+k7],a[j2+k7],a[j1+k8],a[j2+k8],
				rt,it,re
			);	tptr += 9; iptr += 9;
		}
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(k=1; k <= khi; k++) */

