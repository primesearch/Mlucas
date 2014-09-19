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

	/*...The radix-176 DIT pass is here:	*/
	#ifdef USE_SSE2

	//...gather the needed data (176 64-bit complex) and do 11 radix-16 transforms:
		tmp = r00;
		kk = 0;
		for(l = 0; l < 11; l++) {
			i64 = dit16_iidx_lo[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			addr = &a[j1] + phi[kk];	// offset = p10*[0,6,1,7,2,8,3,9,4,a,5], start idx = 0 and decr 5 (mod 11) each loop
			add0=addr+k0; add1=addr+k1; add2=addr+k2; add3=addr+k3; add4=addr+k4; add5=addr+k5; add6=addr+k6; add7=addr+k7; add8=addr+k8; add9=addr+k9; adda=addr+ka; addb=addr+kb; addc=addr+kc; addd=addr+kd; adde=addr+ke; addf=addr+kf;
			SSE2_RADIX16_DIT_0TWIDDLE(
				add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf,
				isrt2,
				tmp,OFF1,OFF2,OFF3,OFF4
			);	tmp += 0x20;
			kk -= 5; kk += ((-(kk < 0)) & 11);
		}
	//...and now do 16 radix-11 transforms:
		tmp = r00;
		for(l = 0; l < 16; l++) {
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;
			va1 = tmp + 0x20;
			va2 = tmp + 0x40;
			va3 = tmp + 0x60;
			va4 = tmp + 0x80;
			va5 = tmp + 0xa0;
			va6 = tmp + 0xc0;
			va7 = tmp + 0xe0;
			va8 = tmp + 0x100;
			va9 = tmp + 0x120;
			vaa = tmp + 0x140;
			tm1 = s1p00 + (((16 - l)&0xf)<<1);	// Low-part offset = p0,f,e,...,2,1
			// Hi-part of p-offset indices:
			// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
			// replacing phi[*] with (*)<<5 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
			iptr = dit_pcshft + dit_ncshft[l];
			k0 = *iptr;			vb0 = tm1 + (k0<<5);
			k1 = *(iptr+0x1);	vb1 = tm1 + (k1<<5);
			k2 = *(iptr+0x2);	vb2 = tm1 + (k2<<5);
			k3 = *(iptr+0x3);	vb3 = tm1 + (k3<<5);
			k4 = *(iptr+0x4);	vb4 = tm1 + (k4<<5);
			k5 = *(iptr+0x5);	vb5 = tm1 + (k5<<5);
			k6 = *(iptr+0x6);	vb6 = tm1 + (k6<<5);
			k7 = *(iptr+0x7);	vb7 = tm1 + (k7<<5);
			k8 = *(iptr+0x8);	vb8 = tm1 + (k8<<5);
			k9 = *(iptr+0x9);	vb9 = tm1 + (k9<<5);
			ka = *(iptr+0xa);	vba = tm1 + (ka<<5);
			SSE2_RADIX_11_DFT(
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,
				ua0,
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba
			);	tmp += 2;
		}

	#else	// USE_SSE2 = False:

	//...gather the needed data (176 64-bit complex) and do 11 radix-16 transforms:
		tptr = t;
		kk = 0;
		for(l = 0; l < 11; l++) {
			i64 = dit16_iidx_lo[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			jp = phi[kk];	// = p10*[0,6,1,7,2,8,3,9,4,a,5], start idx = 0 and decr 5 (mod 11) each loop
			jt = j1 + jp; jp += j2;
			RADIX_16_DIT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
				c16,s16
			);	tptr += 0x10;
			kk -= 5; kk += ((-(kk < 0)) & 11);
		}
	//...and now do 16 radix-11 transforms:
		tptr = t;
		for(l = 0; l < 16; l++) {
			iptr = dit_pcshft + dit_ncshft[l];
			// Hi-part of p-offset indices:
			k0 = phi[*iptr];
			k1 = phi[*(iptr+0x1)];
			k2 = phi[*(iptr+0x2)];
			k3 = phi[*(iptr+0x3)];
			k4 = phi[*(iptr+0x4)];
			k5 = phi[*(iptr+0x5)];
			k6 = phi[*(iptr+0x6)];
			k7 = phi[*(iptr+0x7)];
			k8 = phi[*(iptr+0x8)];
			k9 = phi[*(iptr+0x9)];
			ka = phi[*(iptr+0xa)];
			jp = plo[(16 - l)&0xf];	// Low-part offset = p0,f,e,...,2,1
			jt = j1 + jp; jp += j2;
			RADIX_11_DFT(
				tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],
				a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9
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

/*...The radix-176 DIF pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (176 64-bit complex) and do 16 radix-11 transforms:
		tmp = r00;
		for(l = 0; l < 16; l++) {
			tm1 = s1p00 + ((((l<<2)+l) & 0xf)<<1);	// Low-part offset = p[5*l (mod 16)]
			iptr = dif_pcshft + dif_ncshft[l];
			// Hi-part of p-offset indices:
			// Since SIMD code stores DIF-ins into contig-local mem rather than back into large-strided main-array locs,
			// replacing phi[*] with (*)<<5 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
			k0 = *iptr;			vb0 = tm1 + (k0<<5);
			k1 = *(iptr+0x1);	vb1 = tm1 + (k1<<5);
			k2 = *(iptr+0x2);	vb2 = tm1 + (k2<<5);
			k3 = *(iptr+0x3);	vb3 = tm1 + (k3<<5);
			k4 = *(iptr+0x4);	vb4 = tm1 + (k4<<5);
			k5 = *(iptr+0x5);	vb5 = tm1 + (k5<<5);
			k6 = *(iptr+0x6);	vb6 = tm1 + (k6<<5);
			k7 = *(iptr+0x7);	vb7 = tm1 + (k7<<5);
			k8 = *(iptr+0x8);	vb8 = tm1 + (k8<<5);
			k9 = *(iptr+0x9);	vb9 = tm1 + (k9<<5);
			ka = *(iptr+0xa);	vba = tm1 + (ka<<5);
			// Output-ptrs [va/vb swap roles here vs DIT] are regular-stride offsets of r00:
			va0 = tmp;
			va1 = tmp + 0x20;
			va2 = tmp + 0x40;
			va3 = tmp + 0x60;
			va4 = tmp + 0x80;
			va5 = tmp + 0xa0;
			va6 = tmp + 0xc0;
			va7 = tmp + 0xe0;
			va8 = tmp + 0x100;
			va9 = tmp + 0x120;
			vaa = tmp + 0x140;
			SSE2_RADIX_11_DFT(
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,
				ua0,
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa
			);	tmp += 2;
		}
	//...and now do 11 radix-16 transforms:
		tmp = r00;
		for(l = 0; l < 11; l++) {
			i64 = dif16_oidx_lo[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			addr = &a[j1] + phi[dif_pcshft[l]];	// offset = p0,pa0,p90,...,p10
			add0=addr+k0; add1=addr+k1; add2=addr+k2; add3=addr+k3; add4=addr+k4; add5=addr+k5; add6=addr+k6; add7=addr+k7; add8=addr+k8; add9=addr+k9; adda=addr+ka; addb=addr+kb; addc=addr+kc; addd=addr+kd; adde=addr+ke; addf=addr+kf;
			SSE2_RADIX16_DIF_0TWIDDLE(
				tmp,OFF1,OFF2,OFF3,OFF4,
				isrt2,
				add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf
			);	tmp += 0x20;
		}

	#else	// USE_SSE2 = False:

	//...gather the needed data (176 64-bit complex) and do 16 radix-11 transforms:
		tptr = t;
		for(l = 0; l < 16; l++) {
			iptr = dif_pcshft + dif_ncshft[l];
			// Hi-part of p-offset indices:
			k0 = phi[*iptr];
			k1 = phi[*(iptr+0x1)];
			k2 = phi[*(iptr+0x2)];
			k3 = phi[*(iptr+0x3)];
			k4 = phi[*(iptr+0x4)];
			k5 = phi[*(iptr+0x5)];
			k6 = phi[*(iptr+0x6)];
			k7 = phi[*(iptr+0x7)];
			k8 = phi[*(iptr+0x8)];
			k9 = phi[*(iptr+0x9)];
			ka = phi[*(iptr+0xa)];
			jp = plo[((l<<2)+l) & 0xf];	// Low-part offset = p[5*l (mod 16)] ...
			jt = j1 + jp; jp += j2;		// ... = p0,5,a,f,4,9...
			RADIX_11_DFT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],
				tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,
				a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9
			);	tptr++;
		}
	//...and now do 11 radix-16 transforms:
		tptr = t;
		for(l = 0; l < 11; l++) {
			i64 = dif16_oidx_lo[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			jp = phi[dif_pcshft[l]];	// = p0,pa0,p90,...,p10
			jt = j1 + jp; jp += j2;
			RADIX_16_DIF(
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				c16,s16);	tptr += 0x10;
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
