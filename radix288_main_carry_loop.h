/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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

	/*...The radix-288 DIT pass is here:	*/
	#ifdef USE_SSE2

	//...gather the needed data (288 64-bit complex, i.e. 576 64-bit reals) and do 9 radix-32 transforms...
		tmp = r00;
		for(l = 0; l < ODD_RADIX; l++) {
			add0 = &a[j1+dit_phi[l]];	itmp = (int *)dit_offsets+(l<<5);
			SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	tmp += 64;
		}
	//...and now do 32 radix-9 transforms, with the columns of t*[r,i] output pairs in the above 9x radix-32 set now acting as input rows.
		tmp = r00;
		for(l = 0; l < 32; l++) {
			int kk = dit_p20_lo_offset[l];
			// Extract index (in [0-8]) into circ-shift array used for high parts of p-mults. The [0-8] value is
			// in low 4 bits of kk; the "which length-17 half of the dit_p20_cperms array?" selector is via (kk < 0):
			int ic = ((-(kk < 0)) & 17)	// +/- sign on kk puts us into lower/upper half of the cshift array (base index 0/17)
						+ (kk & 0xf);	// ...and low 4 bits give the element index w.r.to the array-half in question.
			int k0 = dit_p20_cperms[ic], k1 = dit_p20_cperms[ic+1], k2 = dit_p20_cperms[ic+2], k3 = dit_p20_cperms[ic+3], k4 = dit_p20_cperms[ic+4], k5 = dit_p20_cperms[ic+5], k6 = dit_p20_cperms[ic+6], k7 = dit_p20_cperms[ic+7], k8 = dit_p20_cperms[ic+8];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs:
			kk = (kk & 0x7fffffff) >> 4;
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;
			va1 = tmp + 0x40;
			va2 = tmp + 0x80;
			va3 = tmp + 0xc0;
			va4 = tmp + 0x100;
			va5 = tmp + 0x140;
			va6 = tmp + 0x180;
			va7 = tmp + 0x1c0;
			va8 = tmp + 0x200;
			tm1 = s1p00+kk; vb0=tm1+k0; vb1=tm1+k1; vb2=tm1+k2; vb3=tm1+k3; vb4=tm1+k4; vb5=tm1+k5; vb6=tm1+k6; vb7=tm1+k7; vb8=tm1+k8;
			SSE2_RADIX_09_DIT(
				va0,va1,va2,va3,va4,va5,va6,va7,va8,
				ycc1,
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8
			);	tmp += 2;
		}

	#else	// USE_SSE2 = False:

	//...gather the needed data (288 64-bit complex) and do 9 radix-32 transforms:
		tptr = t;
		for(l = 0; l < ODD_RADIX; l++) {
			jt = j1+dit_phi[l]; RADIX_32_DIT((a+jt),dit_offsets+(l<<5),RE_IM_STRIDE, (double *)tptr,t_offsets,1);	tptr += 32;
		}
	//...and now do 32 radix-9 transforms, with the columns of t*[r,i] output pairs in the above 9x radix-32 set now acting as input rows.
		tptr = t;
		for(l = 0; l < 32; l++) {
			int kk = dit_p20_lo_offset[l];
			// Extract index (in [0-8]) into circ-shift array used for high parts of p-mults. The [0-8] value is
			// in low 4 bits of kk; the "which length-17 half of the dit_p20_cperms array?" selector is via (kk < 0):
			int ic = ((-(kk < 0)) & 17)	// +/- sign on kk puts us into lower/upper half of the cshift array (base index 0/17)
						+ (kk & 0xf);	// ...and low 4 bits give the element index w.r.to the array-half in question.
			int k0 = dit_p20_cperms[ic], k1 = dit_p20_cperms[ic+1], k2 = dit_p20_cperms[ic+2], k3 = dit_p20_cperms[ic+3], k4 = dit_p20_cperms[ic+4], k5 = dit_p20_cperms[ic+5], k6 = dit_p20_cperms[ic+6], k7 = dit_p20_cperms[ic+7], k8 = dit_p20_cperms[ic+8];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs:
			kk = (kk & 0x7fffffff) >> 4;
			jt = j1+kk; jp = j2+kk;
			RADIX_09_DIT(
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x100)->re,(tptr+0x100)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],
				rt,it,re
			);	tptr++;
		}

	#endif	/* USE_SSE2 */

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need 28 separate carries.	*/

	#ifdef USE_AVX

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
		n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

	/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		tm1 = s1p00; tmp = cy; itmp = bjmodn;
		AVX_cmplx_carry_norm_errcheck0_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		tm1 += 8; tmp += 1; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			AVX_cmplx_carry_norm_errcheck1_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			tm1 += 8; tmp += 1; itmp += 4;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

		l= j & (nwt-1);
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwt-l  ];
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		SSE2_cmplx_carry_norm_errcheck0_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

		l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

	/*	i =((uint32)(sw - bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#else	// Scalar-double mode:

		l= j & (nwt-1);
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwt-l  ];
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */		
		l = 0; addr = cy; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1   ],a[j2   ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p1],a[j2+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p2],a[j2+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p3],a[j2+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p4,p8,...
			cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	/* #ifdef USE_SSE2 */

/*...The radix-288 DIF pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (288 64-bit complex) and do 32 radix-9 transforms:
		tmp = r00;
		for(l = 0; l < 32; l++) {
			// Input pointers are into s1p** memblock:
			int kk = dif_p20_lo_offset[l];
			// Extract index (in [0-8]) into circ-shift array used for high parts of p-mults. The [0-8] value is
			// in low 4 bits of kk; the "which length-17 half of the dif_p20_cperms array?" selector is via (kk < 0):
			int ic = ((-(kk < 0)) & 17)	// +/- sign on ll puts us into lower/upper half of the cshift array (base index 0/17)
						+ (kk & 0xf);	// ...and low 4 bits give the element index w.r.to the array-half in question.
			int k0 = dif_p20_cperms[ic], k1 = dif_p20_cperms[ic+1], k2 = dif_p20_cperms[ic+2], k3 = dif_p20_cperms[ic+3], k4 = dif_p20_cperms[ic+4], k5 = dif_p20_cperms[ic+5], k6 = dif_p20_cperms[ic+6], k7 = dif_p20_cperms[ic+7], k8 = dif_p20_cperms[ic+8];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs:
			kk = (kk & 0x7fffffff) >> 4;
			tm1 = s1p00+kk; vb0=tm1+k0; vb1=tm1+k1; vb2=tm1+k2; vb3=tm1+k3; vb4=tm1+k4; vb5=tm1+k5; vb6=tm1+k6; vb7=tm1+k7; vb8=tm1+k8;
			// Output-ptrs [va/vb swap roles here vs DIT] are regular-stride offsets of r00:
			va0 = tmp;
			va1 = tmp + 0x40;
			va2 = tmp + 0x80;
			va3 = tmp + 0xc0;
			va4 = tmp + 0x100;
			va5 = tmp + 0x140;
			va6 = tmp + 0x180;
			va7 = tmp + 0x1c0;
			va8 = tmp + 0x200;
			SSE2_RADIX_09_DIF(
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,
				ycc1,
				va0,va1,va2,va3,va4,va5,va6,va7,va8
			);	tmp += 2;
		}
	//...and now do 9 radix-32 transforms:
		tmp = r00;
		for(l = 0; l < ODD_RADIX; l++) {
			add0 = &a[j1+dif_phi[l]];	itmp = (int *)dif_offsets+(l<<5);
			SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);	tmp += 64;
		}

	#else	// USE_SSE2 = False:

	//...gather the needed data (288 64-bit complex) and do 32 radix-9 transforms:
		tptr = t;
		for(l = 0; l < 32; l++) {
			int kk = dif_p20_lo_offset[l];
			// Extract index (in [0-8]) into circ-shift array used for high parts of p-mults. The [0-8] value is
			// in low 4 bits of kk; the "which length-17 half of the dif_p20_cperms array?" selector is via (kk < 0):
			int ic = ((-(kk < 0)) & 17)	// +/- sign on ll puts us into lower/upper half of the cshift array (base index 0/17)
						+ (kk & 0xf);	// ...and low 4 bits give the element index w.r.to the array-half in question.
			int k0 = dif_p20_cperms[ic], k1 = dif_p20_cperms[ic+1], k2 = dif_p20_cperms[ic+2], k3 = dif_p20_cperms[ic+3], k4 = dif_p20_cperms[ic+4], k5 = dif_p20_cperms[ic+5], k6 = dif_p20_cperms[ic+6], k7 = dif_p20_cperms[ic+7], k8 = dif_p20_cperms[ic+8];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs:
			kk = (kk & 0x7fffffff) >> 4;
			jt = j1+kk; jp = j2+kk;
			RADIX_09_DIF(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x100)->re,(tptr+0x100)->im,
				rt,it,re
			);	tptr++;
		}
	//...and now do 9 radix-32 transforms:
		tptr = t;
		for(l = 0; l < ODD_RADIX; l++) {
			jt = j1+dif_phi[l]; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets+(l<<5),RE_IM_STRIDE);	tptr += 32;
		}

	#endif	/* !USE_SSE2 */

	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(k=1; k <= khi; k++) */
