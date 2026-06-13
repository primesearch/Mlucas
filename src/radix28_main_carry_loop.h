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

// This main loop is same for un-and-multithreaded, so stick into a header file
// (can't use a macro because of the #if-enclosed stuff).

for(int k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
{
	/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
	Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
	but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
	(But only ever need to explicitly do this in debug mode).
	*/
	for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

/*...The radix-28 DIT pass is here:	*/
/* EWM: 10/18/04: We swap the odd-index outputs of each of the radix-4 DIT transforms (1<=>3, 5<=>7, etc.) so that the indexing
				  of all the radix-7 transforms (really just the 2nd and 4th of these) winds up being in-place. This allows us
				  to properly re-use the ajp1 variables in the carry-pass version of this routine.
*/
#ifndef USE_SSE2

  /*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 7 radix-4 transforms...*/
					 /*                                      outputs                                      */ /*                          inputs                           */
	RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a1p00r,a1p00i,a1p03r,a1p03i,a1p02r,a1p02i,a1p01r,a1p01i,rt,it);	jt = j1+p12; jp = j2+p12;
	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a1p04r,a1p04i,a1p07r,a1p07i,a1p06r,a1p06i,a1p05r,a1p05i,rt,it);	jt = j1+p24; jp = j2+p24;
	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a1p08r,a1p08i,a1p11r,a1p11i,a1p10r,a1p10i,a1p09r,a1p09i,rt,it);	jt = j1+p08; jp = j2+p08;
	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a1p12r,a1p12i,a1p15r,a1p15i,a1p14r,a1p14i,a1p13r,a1p13i,rt,it);	jt = j1+p20; jp = j2+p20;
	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a1p16r,a1p16i,a1p19r,a1p19i,a1p18r,a1p18i,a1p17r,a1p17i,rt,it);	jt = j1+p04; jp = j2+p04;
	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a1p20r,a1p20i,a1p23r,a1p23i,a1p22r,a1p22i,a1p21r,a1p21i,rt,it);	jt = j1+p16; jp = j2+p16;
	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a1p24r,a1p24i,a1p27r,a1p27i,a1p26r,a1p26i,a1p25r,a1p25i,rt,it);

  /*...and now do 4 radix-7 transforms...*/
						/*                                     inputs                                                */  /*                 intermediates                     */  /*                                        outputs                                              */  /*   sincos consts   */  /* tmps */
  #if LO_ADD
	RADIX_07_DFT_FMA (a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1    ],a[j2    ],a[j1+p08],a[j2+p08],a[j1+p16],a[j2+p16],a[j1+p24],a[j2+p24],a[j1+p04],a[j2+p04],a[j1+p12],a[j2+p12],a[j1+p20],a[j2+p20], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	jt = j1+p03; jp = j2+p03;
	RADIX_07_DFT_FMA (a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	jt = j1+p02; jp = j2+p02;
	RADIX_07_DFT_FMA (a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p04],a[jp+p04], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	jt = j1+p01; jp = j2+p01;
	RADIX_07_DFT_FMA (a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);
  #else
	RADIX_07_DFT_NUSS(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1    ],a[j2    ],a[j1+p08],a[j2+p08],a[j1+p16],a[j2+p16],a[j1+p24],a[j2+p24],a[j1+p04],a[j2+p04],a[j1+p12],a[j2+p12],a[j1+p20],a[j2+p20], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	jt = j1+p03; jp = j2+p03;
	RADIX_07_DFT_NUSS(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	jt = j1+p02; jp = j2+p02;
	RADIX_07_DFT_NUSS(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p04],a[jp+p04], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	jt = j1+p01; jp = j2+p01;
	RADIX_07_DFT_NUSS(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);
  #endif

#else	// USE_SSE2 = True:

  #if !GCC_ASM_FULL_INLINE

	/* Since doing radix-7 in-place here, outputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */

	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p00r,s1p03r,s1p02r,s1p01r)
	add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p04r,s1p07r,s1p06r,s1p05r)
	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p08r,s1p11r,s1p10r,s1p09r)
	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p12r,s1p15r,s1p14r,s1p13r)
	add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p16r,s1p19r,s1p18r,s1p17r)
	add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p20r,s1p23r,s1p22r,s1p21r)
	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p24r,s1p27r,s1p26r,s1p25r)

	/*...and now do 4 radix-7 transforms...*/

   #ifdef USE_AVX2
	SSE2_RADIX_07_DFT(s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r, cc0,two, s1p00r,s1p08r,s1p16r,s1p24r,s1p04r,s1p12r,s1p20r);
	SSE2_RADIX_07_DFT(s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r, cc0,two, s1p07r,s1p15r,s1p23r,s1p03r,s1p11r,s1p19r,s1p27r);
	SSE2_RADIX_07_DFT(s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r, cc0,two, s1p14r,s1p22r,s1p02r,s1p10r,s1p18r,s1p26r,s1p06r);
	SSE2_RADIX_07_DFT(s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r, cc0,two, s1p21r,s1p01r,s1p09r,s1p17r,s1p25r,s1p05r,s1p13r);
   #else
	SSE2_RADIX_07_DFT(s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r, cc0,     s1p00r,s1p08r,s1p16r,s1p24r,s1p04r,s1p12r,s1p20r);
	SSE2_RADIX_07_DFT(s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r, cc0,     s1p07r,s1p15r,s1p23r,s1p03r,s1p11r,s1p19r,s1p27r);
	SSE2_RADIX_07_DFT(s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r, cc0,     s1p14r,s1p22r,s1p02r,s1p10r,s1p18r,s1p26r,s1p06r);
	SSE2_RADIX_07_DFT(s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r, cc0,     s1p21r,s1p01r,s1p09r,s1p17r,s1p25r,s1p05r,s1p13r);
   #endif

  #else	/* GCC-style inline ASM: */

	add0 = &a[j1    ];
	SSE2_RADIX28_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,s1p00r,cc0);

  #endif

#endif	// USE_SSE2 ?

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 28 separate blocks of the A-array, we need 28 separate carries.	*/

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		// Check if current index-interval contains the target index for rotated-residue carry injection.
		// In data-init we set target_idx = -1 on wraparound-carry mini-pass, so if() only taken on full pass:
		if(target_idx == j) {
		#ifdef USE_SSE2
			addr = (double *)s1p00r + target_set;
			*addr += target_cy*(n>>1);	// target_cy = [-2 << within-word-shift]*[DWT weight]*n/2, i.e. includes fwd DWT weight and n/2 factor
		#else
			// target_set in [0,2*RADIX); tidx_mod_stride [even|odd] means shifted-carry goes into [Re|Im] part of the complex FFT datum:
			l = target_set&1;	target_set >>= 1;
			a[j1+poff[target_set>>2]+p0123[target_set&3]+l] += target_cy*(n>>1);
		#endif
			target_idx = -1;
		}

	#ifdef USE_AVX

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		l= j & (nwt-1);						tmp = half_arr + 128;	/* ptr to local storage for the doubled wtl,wtn terms: */
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

	 #ifndef USE_AVX512	// Radix-12,20,24,28 carry routines are special, require avx-512 builds to be done in HIACC-only mode
	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining
														// (Note Radix-28 forbidden for AVX-512, so HIACC handling for that moot here).
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
		i = (!j);
		addr = &prp_mult;
		tmp = s1p00r; tm1 = cy_r; tm2 = cy_r+1; itmp = bjmodn; itm2 = bjmodn+4;
		for(l = 0; l < RADIX>>3; l++) {
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l+l];
			AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
			tmp += 16; tm1 += 2; tm2 += 2; itmp += 8; itm2 += 8; i = 0;
		}
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1,     itmp,      half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);

	  } else {	// HiACC:
	 #endif	// #ifndef USE_AVX512
		/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		tm1 = s1p00r; tmp = cy_r; itmp = bjmodn;
		i = (!j);
		addr = &prp_mult;
		for(l = 0; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[l]);	// poff[] = p0,4,8,...
			AVX_cmplx_carry_norm_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, tm2,p01,p02,p03, addr);
			tm1 += 8; tmp += 1; itmp += 4; i = 0;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
	 #ifndef USE_AVX512
	  }	// USE_AVX: (8-way or 4-way LOACC) or (4-way HIACC) ?
	 #endif

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	  if(inc_arr[0]) {	// LOACC with tunable DWT-weights chaining, incr != 0 check in lieu of (USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) is for ARMv8

		uint32 ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		tm1 = s1p00r; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every few passes:
		incr = inc_arr;
		for(loop = 0; loop < nloop; loop += *incr++)
		{
			ii = loop << 2;	// Reflects 4 independent carry chains being done in each SSE2_cmplx_carry_fast_pow2_errcheck call
			/*** wt_re,wi_re,wt_im,wi_im inits. Cf. radix16_main_carry_loop.h for scalar-macro prototyping of this: ***/
			l = j & (nwt-1);	nwtml = nwt-l;
			n_minus_sil   = n-si[l  ];
			n_minus_silp1 = n-si[l+1];
			sinwt   = si[nwtml  ];
			sinwtm1 = si[nwtml-1];
			wtl     = wt0[    l  ];
			wtn     = wt0[nwtml  ]*scale;
			wtlp1   = wt0[    l+1];
			wtnm1   = wt0[nwtml-1]*scale;

			co2 = co2save;	// Need this for all wts-inits beynd the initial set, due to the co2 = co3 preceding the (j+2) data
			ctmp = (struct complex *)half_arr + 24;	// ptr to local storage for the doubled wtl,wtn terms:
			// (j)-data occupy the 8 xmm-sized slots above the 16 used by fixed auxiliary-data, and overwrite these inits:
			ctmp->re = ctmp->im = wtl;		ctmp += 2;
			ctmp->re = ctmp->im = wtn;		ctmp += 2;
			ctmp->re = ctmp->im = wtlp1;	ctmp += 2;
			ctmp->re = ctmp->im = wtnm1;

			l = (j+2) & (nwt-1);	nwtml = nwt-l;;
			uint32 k0 = n-si[l  ];
			uint32 k1 = n-si[l+1];
			uint32 k2 = si[nwtml  ];
			uint32 k3 = si[nwtml-1];
			wtn     = wt0[nwtml  ]*scale;
			wtl     = wt0[    l  ];
			wtlp1   = wt0[    l+1];
			wtnm1   = wt0[nwtml-1]*scale;

			ctmp = (struct complex *)half_arr + 32;	// (j+2) data start at ctmp + 8
			ctmp->re = ctmp->im = wtl;		ctmp += 2;
			ctmp->re = ctmp->im = wtn;		ctmp += 2;
			ctmp->re = ctmp->im = wtlp1;	ctmp += 2;
			ctmp->re = ctmp->im = wtnm1;

			add1 = &wt1[col  +ii];	/* Don't use add0 here, to avoid need to reload main-array address */
			add2 = &wt1[co2-1-ii];
			add3 = &wt1[co3-1-ii];

			// Since use wt1-array in the wtsinit macro, need to fiddle this here:
			co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
			// *But*: since the init macro does an on-the-fly version of this between j,j+2 portions, external code co2=co3 must come *after* both ctmp-data octets are inited.
			add0 = (double*)(bjmodn+ii);
			SSE2_cmplx_carry_fast_wtsinit(add1,add2,add3, add0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, k0,k1,k2,k3, sse_bw,sse_n)

			addr = &prp_mult;
			for(l = loop; l < loop+*incr; l++) {
				// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data:
				add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
				SSE2_cmplx_carry_fast_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
				tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
			}
		}

	  } else {	// HiACC:

		l= j & (nwt-1);
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwt-l  ];
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		addr = &prp_mult;
		tm1 = s1p00r; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01, addr);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
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

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

	/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00r; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p02,p03, addr);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

	  }	// LOACC or HIACC?

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

	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining

		//...set0 is slightly different from others; divide work into blocks of 4 macro calls:
		l = 0; addr = cy_r; itmp = bjmodn;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			// Re-init weights every 4th macro invocation to keep errors under control:
			cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		}

	  } else {	// HiACC:

		//...set0 is slightly different from others; divide work into blocks of 4 macro calls:
		l = 0; addr = cy_r; itmp = bjmodn;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p01],a[jp+p01],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p02],a[jp+p02],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p03],a[jp+p03],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		}

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

	}		/************************************************************************/
	else	/*                MODULUS_TYPE_FERMAT:                                 */
	{		/************************************************************************/
		addr = &prp_mult;

	#ifdef USE_AVX512

		/* In AVX512 mode, the data are arranged in memory like so, where we view things in 64-byte chunks with R0
		being short for s1p00r, I0 for s1p00i, etc:

			R0 :	a00.re,b00.re,c00.re, ... ,h00.re		I0 :	a00.im,b00.im,c00.im, ... ,h00.im
			R1 :	a01.re,b01.re,c01.re, ... ,h01.re		I1 :	a01.im,b01.im,c01.im, ... ,h01.im
			R2 :	a02.re,b02.re,c02.re, ... ,h02.re		I2 :	a02.im,b02.im,c02.im, ... ,h02.im
			R3 :	a03.re,b03.re,c03.re, ... ,h03.re		I3 :	a03.im,b03.im,c03.im, ... ,h03.im
			R4 :	a04.re,b04.re,c04.re, ... ,h04.re		I4 :	a04.im,b04.im,c04.im, ... ,h04.im
			R5 :	a05.re,b05.re,c05.re, ... ,h05.re		I5 :	a05.im,b05.im,c05.im, ... ,h05.im
			R6 :	a06.re,b06.re,c06.re, ... ,h06.re		I6 :	a06.im,b06.im,c06.im, ... ,h06.im
			R7 :	a07.re,b07.re,c07.re, ... ,h07.re		I7 :	a07.im,b07.im,c07.im, ... ,h07.im
			R8 :	a08.re,b08.re,c08.re, ... ,h08.re		I8 :	a08.im,b08.im,c08.im, ... ,h08.im
			R9 :	a09.re,b09.re,c09.re, ... ,h09.re		I9 :	a09.im,b09.im,c09.im, ... ,h09.im
			R10:	a10.re,b10.re,c10.re, ... ,h10.re		I10:	a10.im,b10.im,c10.im, ... ,h10.im
			R11:	a11.re,b11.re,c11.re, ... ,h11.re		I11:	a11.im,b11.im,c11.im, ... ,h11.im
			R12:	a12.re,b12.re,c12.re, ... ,h12.re		I12:	a12.im,b12.im,c12.im, ... ,h12.im
			R13:	a13.re,b13.re,c13.re, ... ,h13.re		I13:	a13.im,b13.im,c13.im, ... ,h13.im
			R14:	a14.re,b14.re,c14.re, ... ,h14.re		I14:	a14.im,b14.im,c14.im, ... ,h14.im
			R15:	a15.re,b15.re,c15.re, ... ,h15.re		I15:	a15.im,b15.im,c15.im, ... ,h15.im
			R16:	a16.re,b16.re,c16.re, ... ,h16.re		I16:	a16.im,b16.im,c16.im, ... ,h16.im
			R17:	a17.re,b17.re,c17.re, ... ,h17.re		I17:	a17.im,b17.im,c17.im, ... ,h17.im
			R18:	a18.re,b18.re,c18.re, ... ,h18.re		I18:	a18.im,b18.im,c18.im, ... ,h18.im
			R19:	a19.re,b19.re,c19.re, ... ,h19.re		I19:	a19.im,b19.im,c19.im, ... ,h19.im
			R20:	a20.re,b20.re,c20.re, ... ,h20.re		I20:	a20.im,b20.im,c20.im, ... ,h20.im
			R21:	a21.re,b21.re,c21.re, ... ,h21.re		I21:	a21.im,b21.im,c21.im, ... ,h21.im
			R22:	a22.re,b22.re,c22.re, ... ,h22.re		I22:	a22.im,b22.im,c22.im, ... ,h22.im
			R23:	a23.re,b23.re,c23.re, ... ,h23.re		I23:	a23.im,b23.im,c23.im, ... ,h23.im
			R24:	a24.re,b24.re,c24.re, ... ,h24.re		I24:	a24.im,b24.im,c24.im, ... ,h24.im
			R25:	a25.re,b25.re,c25.re, ... ,h25.re		I25:	a25.im,b25.im,c25.im, ... ,h25.im
			R26:	a26.re,b26.re,c26.re, ... ,h26.re		I26:	a26.im,b26.im,c26.im, ... ,h26.im
			R27:	a27.re,b27.re,c27.re, ... ,h27.re		I27:	a27.im,b27.im,c27.im, ... ,h27.im

		The a-h's of each quartet represent doubles which in non-SIMD mode would be getting processed in the same relative
		position on subsequent passes through the main-array for-loop, i.e. among which the carry propagation proceeds as

			a00.re -> b00.re -> ... -> h00.re;		a00.im -> b00.im -> ... -> h00.im ,

		where the imaginary parts really represent elements [a-h]00.im = a[n/2,n/2+1,n/2+2,...,n/2+7] of the right-angle transform.

		Now a crucial thing to note about the above data is that for any value of the main loop index j,
		the length-[odd radix] cyclical nature of the IBDWT weights as a function of j means that all of the a-terms
		in the above get the same IBDWT weights (corr. to index j), all the b-terms get the same (j+2) weights,
		c-terms all get the same (j+4) weights, ... , and h-terms all get the same (j+14) weights.

		The SIMD layout is ideal for the negacyclic unweighting/reweighting step bracketing the carry step, but in the latter,
		because of the undesirable intra-SIMD-register data dependency this leads to, we instead shuffle the above data using
		the same kind of 8 x 8 real-double-array transposition as used around the dyadic-square step. After the shuffle, R0-7
		and I0-7, respectively, contain their original 64 doubles, transposed like so:

										ON INPUT TO CARRY STEP:

			R0 :	a00.re,b00.re,c00.re, ... ,h00.re		I0 :	a00.im,b00.im,c00.im, ... ,h00.im
			R1 :	a01.re,b01.re,c01.re, ... ,h01.re		I1 :	a01.im,b01.im,c01.im, ... ,h01.im
			R2 :	a02.re,b02.re,c02.re, ... ,h02.re		I2 :	a02.im,b02.im,c02.im, ... ,h02.im
														...
			R7 :	a03.re,b03.re,c03.re, ... ,h03.re		I7 :	a03.im,b03.im,c03.im, ... ,h03.im

										SHUFFLE STEP GIVES:

			R0 :	a00.re,a01.re,a02.re, ... ,a07.re		I0 :	a00.im,a01.im,a02.im, ... ,a07.im
			R1 :	b00.re,b01.re,b02.re, ... ,b07.re		I1 :	b00.im,b01.im,b02.im, ... ,b07.im
			R2 :	c00.re,c01.re,c02.re, ... ,c07.re		I2 :	c00.im,c01.im,c02.im, ... ,c07.im
														...
			R3 :	h00.re,h01.re,h02.re, ... ,h07.re		I3 :	h00.im,h01.im,h02.im, ... ,h07.im

		NOTE: Even though e.g. a00 and a01 appear adjacent in terms of their a-subscripts, they are actually
		n/28 memory locations apart, i.e. there is no carry propagation between them.
		*/
		#warning No avx512 version of AVX_cmplx_carry_fast_pow2_wtsinit_X8 available yet!

	#elif defined(USE_AVX)

		/* For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.
		In AVX mode, the data are arranged in memory like so, where we view things in 32-byte chunks with R0
		being short for s1p00r, I0 for s1p00i, etc:

			R0 :	a00.re,b00.re,c00.re,d00.re		I0 :	a00.im,b00.im,c00.im,d00.im
			R1 :	a01.re,b01.re,c01.re,d01.re		I1 :	a01.im,b01.im,c01.im,d01.im
			R2 :	a02.re,b02.re,c02.re,d02.re		I2 :	a02.im,b02.im,c02.im,d02.im
			R3 :	a03.re,b03.re,c03.re,d03.re		I3 :	a03.im,b03.im,c03.im,d03.im
			R4 :	a04.re,b04.re,c04.re,d04.re		I4 :	a04.im,b04.im,c04.im,d04.im
			R5 :	a05.re,b05.re,c05.re,d05.re		I5 :	a05.im,b05.im,c05.im,d05.im
			R6 :	a06.re,b06.re,c06.re,d06.re		I6 :	a06.im,b06.im,c06.im,d06.im
			R7 :	a07.re,b07.re,c07.re,d07.re		I7 :	a07.im,b07.im,c07.im,d07.im
			R8 :	a08.re,b08.re,c08.re,d08.re		I8 :	a08.im,b08.im,c08.im,d08.im
			R9 :	a09.re,b09.re,c09.re,d09.re		I9 :	a09.im,b09.im,c09.im,d09.im
			R10:	a10.re,b10.re,c10.re,d10.re		I10:	a10.im,b10.im,c10.im,d10.im
			R11:	a11.re,b11.re,c11.re,d11.re		I11:	a11.im,b11.im,c11.im,d11.im
			R12:	a12.re,b12.re,c12.re,d12.re		I12:	a12.im,b12.im,c12.im,d12.im
			R13:	a13.re,b13.re,c13.re,d13.re		I13:	a13.im,b13.im,c13.im,d13.im
			R14:	a14.re,b14.re,c14.re,d14.re		I14:	a14.im,b14.im,c14.im,d14.im
			R15:	a15.re,b15.re,c15.re,d15.re		I15:	a15.im,b15.im,c15.im,d15.im
			R16:	a16.re,b16.re,c16.re,d16.re		I16:	a16.im,b16.im,c16.im,d16.im
			R17:	a17.re,b17.re,c17.re,d17.re		I17:	a17.im,b17.im,c17.im,d17.im
			R18:	a18.re,b18.re,c18.re,d18.re		I18:	a18.im,b18.im,c18.im,d18.im
			R19:	a19.re,b19.re,c19.re,d19.re		I19:	a19.im,b19.im,c19.im,d19.im
			R20:	a20.re,b20.re,c20.re,d20.re		I20:	a20.im,b20.im,c20.im,d20.im
			R21:	a21.re,b21.re,c21.re,d21.re		I21:	a21.im,b21.im,c21.im,d21.im
			R22:	a22.re,b22.re,c22.re,d22.re		I22:	a22.im,b22.im,c22.im,d22.im
			R23:	a23.re,b23.re,c23.re,d23.re		I23:	a23.im,b23.im,c23.im,d23.im
			R24:	a24.re,b24.re,c24.re,d24.re		I24:	a24.im,b24.im,c24.im,d24.im
			R25:	a25.re,b25.re,c25.re,d25.re		I25:	a25.im,b25.im,c25.im,d25.im
			R26:	a26.re,b26.re,c26.re,d26.re		I26:	a26.im,b26.im,c26.im,d26.im
			R27:	a27.re,b27.re,c27.re,d27.re		I27:	a27.im,b27.im,c27.im,d27.im

		The a-d's of each quartet represent doubles which in non-SIMD mode would be getting processed in the same relative
		position on subsequent passes through the main-array for-loop, i.e. among which the carry propagation proceeds as

			a00.re -> b00.re -> c00.re -> d00.re;		a00.im -> b00.im -> c00.im -> d00.im ,

		where the imaginary parts really represent elements [a-d]00.im = a[n/2,n/2+1,n/2+2,n/2+3] of the right-angle transform.

		Now a crucial thing to note about the above data is that for any value of the main loop index j,
		the length-[odd radix] cyclical nature of the IBDWT weights as a function of j means that all of the a-terms
		in the above get the same IBDWT weights (corr. to index j), all the b-terms get the same (j+2) weights,
		c-terms all get the same (j+4) weights and d-terms all get the same (j+6) weights.

		The SIMD layout is ideal for the negacyclic unweighting/reweighting step bracketing the carry step, but in the latter,
		because of the undesirable intra-SIMD-register data dependency this leads to, we instead shuffle the above data using
		the same kind of 4 x 4 real-double-array transposition as used around the dyadic-square step. After the shuffle, R0-3
		and I0-3, respectively, contain their original 16 doubles, transposed like so:

										ON INPUT TO CARRY STEP:

			R0 :	a00.re,b00.re,c00.re,d00.re		I0 :	a00.im,b00.im,c00.im,d00.im
			R1 :	a01.re,b01.re,c01.re,d01.re		I1 :	a01.im,b01.im,c01.im,d01.im
			R2 :	a02.re,b02.re,c02.re,d02.re		I2 :	a02.im,b02.im,c02.im,d02.im
			R3 :	a03.re,b03.re,c03.re,d03.re		I3 :	a03.im,b03.im,c03.im,d03.im

										SHUFFLE STEP GIVES:

			R0 :	a00.re,a01.re,a02.re,a03.re		I0 :	a00.im,a01.im,a02.im,a03.im
			R1 :	b00.re,b01.re,b02.re,b03.re		I1 :	b00.im,b01.im,b02.im,b03.im
			R2 :	c00.re,c01.re,c02.re,c03.re		I2 :	c00.im,c01.im,c02.im,c03.im
			R3 :	d00.re,d01.re,d02.re,d03.re		I3 :	d00.im,d01.im,d02.im,d03.im

		NOTE #1: This is slightly different for AVX than SSE2 - in SSE2 we further interleaved re/im data.

		NOTE #2: Even though e.g. a00 and a01 appear adjacent in terms of their a-subscripts, they are actually
		n/28 memory locations apart, i.e. there is no carry propagation between them.
		*/

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		tmp = base_negacyclic_root;	tm2 = tmp+1;

	  #ifdef HIACC
		uint32 k1, k2;
		// Hi-accuracy version needs RADIX/4 copies of each base root:
		l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 8) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 8) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 8) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 8) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}

		// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indices j,j+2,j+4,j+6:

		/* The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
		which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls: */

		tm0 = s1p00r; tmp = base_negacyclic_root; l = 0x700;
		tm1 = cy_r; // *cycle[] indices increment by +4 (mod ODD_RADIX) between macro calls
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic_idx = 0; jc_idx = 1; kc_idx = 2; lc_idx = 3;
		while(tm0 < s1p27r)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{	// NB: (int)(tmp-cy_r) < RADIX (as used for SSE2 build) no good here, since just 1 vec_dbl increment
			// per 4 Re+Im-carries; but (int)(tmp-cy_r) < (RADIX>>1) would work
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			uint32 k1 = icycle[ic_idx], k5 = jcycle[ic_idx], k6 = kcycle[ic_idx], k7 = lcycle[ic_idx];
			uint32 k2 = icycle[jc_idx];
			uint32 k3 = icycle[kc_idx];
			uint32 k4 = icycle[lc_idx];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0xe0, 0xe0,0x1c0,0x2a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p01,p02,p03, addr);
			tm0 += 8; tm1++; tmp += 8; l -= 0xc0;
			MOD_ADD32(ic_idx, 4, ODD_RADIX, ic_idx);
			MOD_ADD32(jc_idx, 4, ODD_RADIX, jc_idx);
			MOD_ADD32(kc_idx, 4, ODD_RADIX, kc_idx);
			MOD_ADD32(lc_idx, 4, ODD_RADIX, lc_idx);
		}

	  #else	// HIACC = false:

		// Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
		// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
		l = (j >> 1);
		for(i = 0; i < RE_IM_STRIDE; i++) {
			uint32 k1=(l & NRTM1), k2=(l >> NRT_BITS);
			dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
			rt  =rn1[k2].re;			it   =rn1[k2].im;
			wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
			VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;
			l += 1;
		}

		// The above need some inits to prepare for the AVX version of the Fermat-mod carry macro:
		SSE2_fermat_carry_init_loacc(base_negacyclic_root);

		tm0 = s1p00r; tmp = base_negacyclic_root;	// tmp *not* incremented between macro calls in loacc version
		tm1 = cy_r;
		ic_idx = 0; jc_idx = 1; kc_idx = 2; lc_idx = 3;

		for(l = 0; l < RADIX>>2; l++) {	// RADIX/4 loop passes
			uint32 k1 = icycle[ic_idx];
			uint32 k2 = icycle[jc_idx], k5 = jcycle[ic_idx];
			uint32 k3 = icycle[kc_idx], k6 = kcycle[ic_idx];
			uint32 k4 = icycle[lc_idx], k7 = lcycle[ic_idx];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
												/* (cy_i_cy_r) --vvvvv  vvvvvvvvvvvvvvvvv--[1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0xe0, 0xe0,0x1c0,0x2a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p01,p02,p03, addr);
			tm0 += 8; tm1++;
			MOD_ADD32(ic_idx, 4, ODD_RADIX, ic_idx);
			MOD_ADD32(jc_idx, 4, ODD_RADIX, jc_idx);
			MOD_ADD32(kc_idx, 4, ODD_RADIX, kc_idx);
			MOD_ADD32(lc_idx, 4, ODD_RADIX, lc_idx);
		}

	  #endif	/* HIACC? */

	#elif defined(USE_SSE2)

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],jcycle[ic],icycle[ic+1],jcycle[ic+1] of the non-looped version with icycle[ic],jcycle[ic],icycle[jc],jcycle[jc]:
		ic_idx = 0; jc_idx = 1;
		tm1 = s1p00r; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while((int)(tmp-cy_r) < RADIX) {
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			uint32 k1 = icycle[ic_idx];
			uint32 k2 = jcycle[ic_idx];
			uint32 k3 = icycle[jc_idx];
			uint32 k4 = jcycle[jc_idx];
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tmp-cy_r)>>2]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
			tm2 += (-((int)((tmp-cy_r)>>1)&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_fermat_carry_norm_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2,k3,k4, tm2,p01, addr);
			tm1 += 4; tmp += 2;
			MOD_ADD32(ic_idx, 2, ODD_RADIX, ic_idx);
			MOD_ADD32(jc_idx, 2, ODD_RADIX, jc_idx);
		}

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		uint32 k1, k2;
		ntmp = 0; addr = cy_r; addi = cy_i; ic_idx = 0;	// ic_idx = idx into icycle mini-array, gets incremented (mod ODD_RADIX) between macro calls
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];
			fermat_carry_norm_errcheckB(a[jt    ],a[jp    ],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p01],a[jp+p01],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p02],a[jp+p02],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p03],a[jp+p03],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
		}
	//	if(!jt)printf("Iter = %d, full = %u, a[0]_out = %20.15f\n",iter,full_pass,a[0]);
		for(ntmp = 0; ntmp < ODD_RADIX; ntmp++) {
			icycle[ntmp] += wts_idx_incr;	// Inside the loop use this, as it is faster than general-mod '% nwt'
			icycle[ntmp] += ( (-(int)((uint32)icycle[ntmp] >> 31)) & nwt);
		}

	#endif	/* #ifdef USE_SSE2 */

	// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
	#ifdef USE_SSE2
		for(ntmp = 0; ntmp < ODD_RADIX; ntmp++)
		{
			icycle[ntmp] += wts_idx_inc2;		icycle[ntmp] += ( (-(icycle[ntmp] < 0)) & nwt16);
			jcycle[ntmp] += wts_idx_inc2;		jcycle[ntmp] += ( (-(jcycle[ntmp] < 0)) & nwt16);
		#ifdef USE_AVX
			kcycle[ntmp] += wts_idx_inc2;		kcycle[ntmp] += ( (-(kcycle[ntmp] < 0)) & nwt16);
			lcycle[ntmp] += wts_idx_inc2;		lcycle[ntmp] += ( (-(lcycle[ntmp] < 0)) & nwt16);
		#endif
		}
	#endif

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-28 DIF pass is here:	*/

#ifndef USE_SSE2

  #if PFETCH
	addr = &a[j1];
	prefetch_p_doubles(addr);
  #endif

/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 4 radix-7 transforms...*/
						/*                                     inputs                                                  */  /*                 intermediates                     */  /*                                        outputs                                              */  /*   sincos consts   */  /* tmps, prefetch addresses */
  #if LO_ADD
   #if PFETCH
	RADIX_07_DFT_PFETCH(a[j1    ],a[j2    ],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im, addr,addp,p01,p02,p03);			jt = j1+p01; jp = j2+p01;
	RADIX_07_DFT_PFETCH(a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im, addr,addp,p04,p04+p01,p04+p02);	jt = j1+p02; jp = j2+p02;
	RADIX_07_DFT_PFETCH(a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im, addr,addp,p04+p01,p08,p08+p01);	jt = j1+p03; jp = j2+p03;
	RADIX_07_DFT_PFETCH(a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im, addr,addp,p08+p02,p08+p03,p12);
   #else
	RADIX_07_DFT       (a[j1    ],a[j2    ],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	jt = j1+p01; jp = j2+p01;
	RADIX_07_DFT       (a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	jt = j1+p02; jp = j2+p02;
	RADIX_07_DFT       (a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	jt = j1+p03; jp = j2+p03;
	RADIX_07_DFT       (a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);
   #endif
  #else
	RADIX_07_DFT_NUSS  (a[j1    ],a[j2    ],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	jt = j1+p01; jp = j2+p01;
	RADIX_07_DFT_NUSS  (a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	jt = j1+p02; jp = j2+p02;
	RADIX_07_DFT_NUSS  (a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	jt = j1+p03; jp = j2+p03;
	RADIX_07_DFT_NUSS  (a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);
  #endif

/*...and now do 7 radix-4 transforms...*/
					 /*                          inputs                           */ /*                                      outputs                                      */
  #if PFETCH
	addp = addr+p12+p01;
	prefetch_p_doubles(addp);

	RADIX_04_DIF_PFETCH(a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it,addr,addp,p12+p02,p12+p03);	jt = j1+p24; jp = j2+p24;
	RADIX_04_DIF_PFETCH(a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it,addr,addp,p16    ,p16+p01);	jt = j1+p20; jp = j2+p20;
	RADIX_04_DIF_PFETCH(a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,addr,addp,p16+p02,p16+p03);	jt = j1+p16; jp = j2+p16;
	RADIX_04_DIF_PFETCH(a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it,addr,addp,p20    ,p20+p01);	jt = j1+p12; jp = j2+p12;
	RADIX_04_DIF_PFETCH(a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it,addr,addp,p20+p02,p20+p03);	jt = j1+p08; jp = j2+p08;
	RADIX_04_DIF_PFETCH(a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it,addr,addp,p24    ,p20+p01);	jt = j1+p04; jp = j2+p04;
	RADIX_04_DIF_PFETCH(a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,addr,addp,p20+p02,p20+p03);
  #else
	RADIX_04_DIF       (a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);	jt = j1+p24; jp = j2+p24;
	RADIX_04_DIF       (a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p20; jp = j2+p20;
	RADIX_04_DIF       (a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p16; jp = j2+p16;
	RADIX_04_DIF       (a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p12; jp = j2+p12;
	RADIX_04_DIF       (a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p08; jp = j2+p08;
	RADIX_04_DIF       (a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p04; jp = j2+p04;
	RADIX_04_DIF       (a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
  #endif

#else	// USE_SSE2 = True:

  #if !GCC_ASM_FULL_INLINE

   #ifdef USE_AVX2
	SSE2_RADIX_07_DFT(s1p00r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r, cc0,two, s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r);
	SSE2_RADIX_07_DFT(s1p21r,s1p17r,s1p13r,s1p09r,s1p05r,s1p01r,s1p25r, cc0,two, s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r);
	SSE2_RADIX_07_DFT(s1p14r,s1p10r,s1p06r,s1p02r,s1p26r,s1p22r,s1p18r, cc0,two, s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r);
	SSE2_RADIX_07_DFT(s1p07r,s1p03r,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r, cc0,two, s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r);
   #else
	SSE2_RADIX_07_DFT(s1p00r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r, cc0,     s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r);
	SSE2_RADIX_07_DFT(s1p21r,s1p17r,s1p13r,s1p09r,s1p05r,s1p01r,s1p25r, cc0,     s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r);
	SSE2_RADIX_07_DFT(s1p14r,s1p10r,s1p06r,s1p02r,s1p26r,s1p22r,s1p18r, cc0,     s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r);
	SSE2_RADIX_07_DFT(s1p07r,s1p03r,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r, cc0,     s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r);
   #endif

	/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */

	add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00r,s1p01r,s1p02r,s1p03r, add0,add1,add2,add3)
	add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04r,s1p05r,s1p06r,s1p07r, add0,add1,add2,add3)
	add3 = &a[j1+p20];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08r,s1p09r,s1p10r,s1p11r, add0,add1,add2,add3)
	add0 = &a[j1+p16];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12r,s1p13r,s1p14r,s1p15r, add0,add1,add2,add3)
	add2 = &a[j1+p12];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16r,s1p17r,s1p18r,s1p19r, add0,add1,add2,add3)
	add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20r,s1p21r,s1p22r,s1p23r, add0,add1,add2,add3)
	add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24r,s1p25r,s1p26r,s1p27r, add0,add1,add2,add3)

  #else

	add0 = &a[j1    ];
	SSE2_RADIX28_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,s1p00r,cc0);

  #endif

#endif	// USE_SSE2 ?

	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(int k=1; k <= khi; k++) */
