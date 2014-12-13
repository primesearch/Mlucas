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

	/*...The radix-768 DIT pass is here:	*/

	#ifdef USE_SSE2

		// Use s1p000-0ff for scratch for these 3 DFTs:
		SSE2_RADIX256_DIT(
			(a+j1     ),(int *)(dit_offsets_lo+0x00), dit_idx1, dit_offsets_hi1, 
			s1p000, isrt2,two,twid0,
			r000
		);	// Outputs in r0[ 00- ff]
		SSE2_RADIX256_DIT(
			(a+j1+p200),(int *)(dit_offsets_lo+0x20), dit_idx2, dit_offsets_hi2, 
			s1p000, isrt2,two,twid0,
			r100
		);	// Outputs in r1[100-1ff]
		SSE2_RADIX256_DIT(
			(a+j1+p100),(int *)(dit_offsets_lo+0x40), dit_idx3, dit_offsets_hi3, 
			s1p000, isrt2,two,twid0,
			r200
		);	// Outputs in r2[200-2ff]

	/*...and now do 256 radix-3 transforms: */

		vec_dbl *vd0,*vd1,*vd2;
		// Loop-based compact-obj-code impl exploits above index pattern to group interior [sandwiched
		// between single leading and 15 trailing] macro calls into sets of 16 with neatly cutoff index groupings:
		l = 2; tmp = r000+l; tm1 = tmp+512; tm2 = tmp+1024;	// Skip 0-term, which gets saved for wraparound
		for(l1 = 0; l1 < 48; l1 += 3) {
			k0 = dit_triplets[l1]; k1 = dit_triplets[l1+1]; k2 = dit_triplets[l1+2];
			l2 = 0x1e;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x1c;	vd0 = s1p000+l2; vd1 = vd0+k2; vd2 = vd0+k0; vd0 += k1; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x1a;	vd0 = s1p000+l2; vd1 = vd0+k0; vd2 = vd0+k1; vd0 += k2; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x18;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x16;	vd0 = s1p000+l2; vd1 = vd0+k2; vd2 = vd0+k0; vd0 += k1; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x14;	vd0 = s1p000+l2; vd1 = vd0+k0; vd2 = vd0+k1; vd0 += k2; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x12;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x10;	vd0 = s1p000+l2; vd1 = vd0+k2; vd2 = vd0+k0; vd0 += k1; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x0e;	vd0 = s1p000+l2; vd1 = vd0+k0; vd2 = vd0+k1; vd0 += k2; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x0c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x0a;	vd0 = s1p000+l2; vd1 = vd0+k2; vd2 = vd0+k0; vd0 += k1; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x08;	vd0 = s1p000+l2; vd1 = vd0+k0; vd2 = vd0+k1; vd0 += k2; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x06;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x04;	vd0 = s1p000+l2; vd1 = vd0+k2; vd2 = vd0+k0; vd0 += k1; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x02;	vd0 = s1p000+l2; vd1 = vd0+k0; vd2 = vd0+k1; vd0 += k2; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; l &= 0x1ff; tmp = r000+l; tm1 = tmp+512; tm2 = tmp+1024;	/* needed for final-loop-pass wraparound of 0-term */
						vd0 = s1p000+k0; vd1 = s1p000+k1; vd2 = s1p000+k2;      SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd1,vd2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
		}

	#else	/* !USE_SSE2 */

	/*...gather the needed data (768 64-bit complex) and do 3 radix-256 transforms,	*/

		RADIX_256_DIT((a+j1     ),RE_IM_STRIDE,(int *)(dit_offsets_lo+0x00),dit_idx1,dit_offsets_hi1, (double *)(t+0x000),1,t_offsets_lo,t_offsets_hi);	/* Outputs in t[ 00- ff] */
		RADIX_256_DIT((a+j1+p200),RE_IM_STRIDE,(int *)(dit_offsets_lo+0x20),dit_idx2,dit_offsets_hi2, (double *)(t+0x100),1,t_offsets_lo,t_offsets_hi);	/* Outputs in t[100-1ff] */
		RADIX_256_DIT((a+j1+p100),RE_IM_STRIDE,(int *)(dit_offsets_lo+0x40),dit_idx3,dit_offsets_hi3, (double *)(t+0x200),1,t_offsets_lo,t_offsets_hi);	/* Outputs in t[200-2ff] */

	/*...and now do 256 radix-3 transforms: */

		// Loop-based compact-obj-code impl exploits above index pattern to group interior [sandwiched between single leading
		// and 15 trailing, which also get fused into a final group] macro calls into sets of 16 with neatly cutoff index groupings:
		l = 1; l1 = l+256; l2 = l+512;	// Skip 0-term, which gets saved for wraparound
		for(m = 0; m < 48; m += 3) {
			k0 = dit_triplets[m]; k1 = dit_triplets[m+1]; k2 = dit_triplets[m+2];
			jt = j1 + pf; jp = j2 + pf;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pe; jp = j2 + pe;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + pd; jp = j2 + pd;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + pc; jp = j2 + pc;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pb; jp = j2 + pb;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + pa; jp = j2 + pa;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p9; jp = j2 + p9;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p8; jp = j2 + p8;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p7; jp = j2 + p7;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p6; jp = j2 + p6;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p5; jp = j2 + p5;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p4; jp = j2 + p4;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p2; jp = j2 + p2;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p1; jp = j2 + p1;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1]); ++l; l &= 0xff; l1 = l+256; l2 = l+512;	// <*** needed for final-loop-pass wraparound of 0-term
										RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2]); ++l; ++l1; ++l2;
		}

	#endif

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 48 separate blocks of the A-array, we need 48 separate carries.	*/

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
		tm1 = s1p000; tmp = cy; itmp = bjmodn;
		// Each AVX carry macro call also processes 4 prefetches of main-array data
		tm2 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_norm_errcheck0_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, tm2,p1,p2,p3);
		tm1 += 8; tmp += 1; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = a + j1 + pfetch_dist + poff[l];
			AVX_cmplx_carry_norm_errcheck1_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, tm2,p1,p2,p3);
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

	/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p000; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
		add0 = a + j1 + pfetch_dist;
		SSE2_cmplx_carry_norm_errcheck0_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p000; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p2,p3);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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

		/*...set0 is slightly different from others:	*/
		l = 0; addr = cy; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1   ],a[j2   ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p1],a[j2+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p2],a[j2+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p3],a[j2+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(l1 = 1; l1 < RADIX>>2; l1++, l2 -= 3) {
			jt = j1 + poff[l1]; jp = j2 + poff[l1];
			cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}
	//	ASSERT(HERE, k0 == p2f0, "loop-exit value of k0 bad!");
		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

	/*...The radix-768 DIF pass is here:	*/

	#ifdef USE_SSE2

	/*...gather the needed data (768 64-bit complex) and do 256 radix-3 transforms - We want unit-strides in the radix768-DFT macro, so use large output strides here: */

		l = 2; tmp = r000+l; tm1 = tmp+512; tm2 = tmp+1024;	// Skip 0-term, which gets saved for wraparound
		for(l1 = 0; l1 < 144; l1 += 3) {
			k0 = dif_triplets[l1]; k1 = dif_triplets[l1+1]; k2 = dif_triplets[l1+2]; l1 += 3;
			l2 = 0x1a;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x14;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x0e;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x08;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x02;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
						k0 = dif_triplets[l1]; k1 = dif_triplets[l1+1]; k2 = dif_triplets[l1+2]; l1 += 3;
			l2 = 0x1c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x16;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x10;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x0a;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x04;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
						k0 = dif_triplets[l1]; k1 = dif_triplets[l1+1]; k2 = dif_triplets[l1+2];	// Skip k-incr here since loop control handles this one
			l2 = 0x1e;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x18;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x12;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x0c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
			l2 = 0x06;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; l &= 0x1ff; tmp = r000+l; tm1 = tmp+512; tm2 = tmp+1024;	/* needed for final-loop-pass wraparound of 0-term */
						vd0 = s1p000+k0; vd1 = s1p000+k1; vd2 = s1p000+k2;      SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp = r000+l; tm1+=2; tm2+=2;
		}

	/*...and now do 3 radix-256 transforms: */

		SSE2_RADIX256_DIF(	// Inputs in r0[ 00- ff]
			r000,
			s1p000, isrt2,two,twid0,
			(a+j1     ),dif_offsets_lo, dif_idx1, dif_offsets_hi1
		);
		SSE2_RADIX256_DIF(	// Inputs in r1[100-1ff]
			r100,
			s1p000, isrt2,two,twid0,
			(a+j1+p200),dif_offsets_lo, dif_idx2, dif_offsets_hi2
		);
		SSE2_RADIX256_DIF(	// Inputs in r2[200-2ff]
			r200,
			s1p000, isrt2,two,twid0,
			(a+j1+p100),dif_offsets_lo, dif_idx3, dif_offsets_hi3
		);

	#else	/* !USE_SSE2 */

	/*...gather the needed data (768 64-bit complex) and do 256 radix-3 transforms - We want unit-strides in the radix768-DFT macro, so use large output strides here: */

		// Loop-based compact-obj-code impl exploits above index pattern to group interior [sandwiched between single leading
		// and 15 trailing, which also get fused into a final group] macro calls into sets of 16 with neatly cutoff index groupings:
		l = 1; l1 = l+256; l2 = l+512;	// Skip 0-term, which gets saved for wraparound
		for(m = 0; m < 144; m += 3) {
			k0 = dif_triplets[m]; k1 = dif_triplets[m+1]; k2 = dif_triplets[m+2]; m += 3;
			jt = j1 + pd; jp = j2 + pd;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pa; jp = j2 + pa;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p7; jp = j2 + p7;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p4; jp = j2 + p4;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p1; jp = j2 + p1;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			k0 = dif_triplets[m]; k1 = dif_triplets[m+1]; k2 = dif_triplets[m+2]; m += 3;
			jt = j1 + pe; jp = j2 + pe;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pb; jp = j2 + pb;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p8; jp = j2 + p8;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p5; jp = j2 + p5;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p2; jp = j2 + p2;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			k0 = dif_triplets[m]; k1 = dif_triplets[m+1]; k2 = dif_triplets[m+2];	// Skip k-incr here since loop control handles this one
			jt = j1 + pf; jp = j2 + pf;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pc; jp = j2 + pc;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p9; jp = j2 + p9;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p6; jp = j2 + p6;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; l &= 0xff; l1 = l+256; l2 = l+512;	// <*** needed for final-loop-pass wraparound of 0-term
										RADIX_03_DFT(s,c3m1, a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
		}

	/*...and now do 3 radix-256 transforms: */

		RADIX_256_DIF((double *)(t+0x000),1,t_offsets_lo,t_offsets_hi, (a+j1     ),RE_IM_STRIDE,dif_offsets_lo,dif_idx1,dif_offsets_hi1);	/* Inputs in t[ 00- ff] */
		RADIX_256_DIF((double *)(t+0x100),1,t_offsets_lo,t_offsets_hi, (a+j1+p200),RE_IM_STRIDE,dif_offsets_lo,dif_idx2,dif_offsets_hi2);	/* Inputs in t[100-1ff] */
		RADIX_256_DIF((double *)(t+0x200),1,t_offsets_lo,t_offsets_hi, (a+j1+p100),RE_IM_STRIDE,dif_offsets_lo,dif_idx3,dif_offsets_hi3);	/* Inputs in t[200-2ff] */

	#endif
	}

	jstart += nwt;
	jhi    += nwt;

	col += RADIX;
	co3 -= RADIX;
}	/* end for(k=1; k <= khi; k++) */

