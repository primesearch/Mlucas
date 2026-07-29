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
	for(j = jstart; j < jhi; j += stride)
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	/*...The radix-24 DIT pass is here:	*/
	/*...gather the needed data (24 64-bit complex, i.e. 48 64-bit reals) and do 3 radix-8 transforms,	*/

	#ifdef USE_SSE2

	  #if USE_SMALL_MACROS

		// SSE2_RADIX8_DIT_0TWIDDLE( add0     +p[0,1,3,2,7,6,5,4], s1p00)
		add0 = &a[j1];
		add1 = add0+p01; add2 = add0+p03; add3 = add0+p02; add4 = add0+p07; add5 = add0+p06; add6 = add0+p05; add7 = add0+p04;
	#ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p00, isrt2,two);
	#else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p00, isrt2);
	#endif
		// SSE2_RADIX8_DIT_0TWIDDLE([add0+p08]+p[5,4,6,7,1,0,2,3], s1p08)
		add5 = &a[j1] + p08;
		add0 = add5+p05; add1 = add5+p04; add2 = add5+p06; add3 = add5+p07; add4 = add5+p01; add6 = add5+p02; add7 = add5+p03;
	#ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p08, isrt2,two);
	#else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p08, isrt2);
	#endif
		// SSE2_RADIX8_DIT_0TWIDDLE([add0+p16]+p[2,3,0,1,4,5,7,6], s1p16)
		add2 = &a[j1] + p16;
		add0 = add2+p02; add1 = add2+p03; add3 = add2+p01; add4 = add2+p04; add5 = add2+p05; add6 = add2+p07; add7 = add2+p06;
	#ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p16, isrt2,two);
	#else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p16, isrt2);
	#endif
	   #if 0	// doubled-data _X2 version slower on x86 SSE2 but slightly faster in ARMv8 SIMD, latter is
	   			// our main target for this assemble-24-DFTs-from-small-macros code, so make _X2 the default:
		SSE2_RADIX_03_DFT(s1p00,s1p08,s1p16,cc3,s1p00,s1p16,s1p08);
		SSE2_RADIX_03_DFT(s1p01,s1p09,s1p17,cc3,s1p09,s1p01,s1p17);
		SSE2_RADIX_03_DFT(s1p02,s1p10,s1p18,cc3,s1p18,s1p10,s1p02);
		SSE2_RADIX_03_DFT(s1p03,s1p11,s1p19,cc3,s1p03,s1p19,s1p11);
		SSE2_RADIX_03_DFT(s1p04,s1p12,s1p20,cc3,s1p12,s1p04,s1p20);
		SSE2_RADIX_03_DFT(s1p05,s1p13,s1p21,cc3,s1p21,s1p13,s1p05);
		SSE2_RADIX_03_DFT(s1p06,s1p14,s1p22,cc3,s1p06,s1p22,s1p14);
		SSE2_RADIX_03_DFT(s1p07,s1p15,s1p23,cc3,s1p15,s1p07,s1p23);
	   #else
		SSE2_RADIX_03_DFT_X2(cc3, s1p00,s1p08,s1p16, s1p00,s1p16,s1p08,  s1p01,s1p09,s1p17, s1p09,s1p01,s1p17);
		SSE2_RADIX_03_DFT_X2(cc3, s1p02,s1p10,s1p18, s1p18,s1p10,s1p02,  s1p03,s1p11,s1p19, s1p03,s1p19,s1p11);
		SSE2_RADIX_03_DFT_X2(cc3, s1p04,s1p12,s1p20, s1p12,s1p04,s1p20,  s1p05,s1p13,s1p21, s1p21,s1p13,s1p05);
		SSE2_RADIX_03_DFT_X2(cc3, s1p06,s1p14,s1p22, s1p06,s1p22,s1p14,  s1p07,s1p15,s1p23, s1p15,s1p07,s1p23);
	   #endif

	  #else

		add0 = &a[j1    ];
		SSE2_RADIX24_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p16,s1p00,isrt2,cc3);

	  #endif

	#else

		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,rt,it);

	/*...and now do 8 in-place radix-3 transforms.	*/

		RADIX_03_DFT(s,c3m1,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,t01,t02,t03,t04,t05,t06,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08]);	jt = j1+p01; jp = j2+p01;
		RADIX_03_DFT(s,c3m1,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,t01,t02,t03,t04,t05,t06,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);	jt = j1+p02; jp = j2+p02;
		RADIX_03_DFT(s,c3m1,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,t01,t02,t03,t04,t05,t06,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ]);	jt = j1+p03; jp = j2+p03;
		RADIX_03_DFT(s,c3m1,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,t01,t02,t03,t04,t05,t06,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08]);	jt = j1+p04; jp = j2+p04;
		RADIX_03_DFT(s,c3m1,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,t01,t02,t03,t04,t05,t06,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);	jt = j1+p05; jp = j2+p05;
		RADIX_03_DFT(s,c3m1,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,t01,t02,t03,t04,t05,t06,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ]);	jt = j1+p06; jp = j2+p06;
		RADIX_03_DFT(s,c3m1,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,t01,t02,t03,t04,t05,t06,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08]);	jt = j1+p07; jp = j2+p07;
		RADIX_03_DFT(s,c3m1,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,t01,t02,t03,t04,t05,t06,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);

	#endif	/* USE_SSE2 */

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 24 separate blocks of the A-array, we need 24 separate carries.	*/

		// Check if current index-interval contains the target index for rotated-residue carry injection.
		// In data-init we set target_idx = -1 on wraparound-carry mini-pass, so if() only taken on full pass:
		if(target_idx == j) {
		#ifdef USE_SSE2
			addr = (double *)s1p00 + target_set;
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

		/* ptr to local storage for the doubled wtl,wtn terms: */
	  #ifdef USE_AVX512
		tmp = half_arr +  64;	// No lookup-tables used in avx-512; instead use opmasked conditional-doubling; 1st 32 slots hold outputs of wtsinit call
	  #else
		tmp = half_arr + 128;	// 1st 64 slots are basic-4 LUTs, next 32 are the additional 2 LOACC LUTs, next 32 hold outputs of wtsinit call
	  #endif
		l= j & (nwt-1);						// These rcol wts-terms are for individual-double-broadcast-to-full-vector-width,
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
	  #ifdef USE_AVX512
		l= (j+8) & (nwt-1);					tmp -= 3;	// Reset to same tmp-startval as above, now copy data into d4-7 slots of vec_dbl
		n_minus_sil  ->d4 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d4 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d4 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d4 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+10) & (nwt-1);				++tmp;
		n_minus_sil  ->d5 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d5 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d5 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d5 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+12) & (nwt-1);				++tmp;
		n_minus_sil  ->d6 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d6 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d6 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d6 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+14) & (nwt-1);				++tmp;
		n_minus_sil  ->d7 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d7 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d7 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d7 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;
	  #endif

	 #ifndef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512, so just use LOACC irrespective of carry-chain setting
	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining
	 #endif	// #ifndef USE_AVX512
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn00, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

		i = (!j);
		addr = &prp_mult;
		// Each carry macro call also processes 8 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
	  #ifdef USE_AVX512	// In AVX-512 mode, the 4 doubles base[0],baseinv[1],wts_mult[1],inv_mult[0] are in the d0-3 slots of the otherwise-unused sse2_rnd vec_dbl:
		AVX_cmplx_carry_fast_errcheck_X8(s1p00, cy00     , bjmodn00         , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr); i = 0;
		add0 = a + j1 + pfetch_dist + p08;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_fast_errcheck_X8(s1p08, cy08     , bjmodn08         , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
		add0 = a + j1 + pfetch_dist + p16;
		AVX_cmplx_carry_fast_errcheck_X8(s1p16, cy16     , bjmodn16         , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
	  #else	// USE_AVX:
		AVX_cmplx_carry_fast_errcheck_X8(s1p00, cy00,cy04, bjmodn00,bjmodn04, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr); i = 0;
		add0 = a + j1 + pfetch_dist + p08;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_fast_errcheck_X8(s1p08, cy08,cy12, bjmodn08,bjmodn12, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
		add0 = a + j1 + pfetch_dist + p16;
		AVX_cmplx_carry_fast_errcheck_X8(s1p16, cy16,cy20, bjmodn16,bjmodn20, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
	  #endif

	 #ifndef USE_AVX512	// Radix-12,20,24,28 carry routines are special, require avx-512 builds to be done in HIACC-only mode
	  } else {	// HiACC:

		// Each AVX carry macro call also processes 4 prefetches of main-array data
		i = (!j);
		addr = &prp_mult;
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_norm_errcheck_X4(s1p00,add1,add2,add3,cy00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr); i = 0;
		add0 += p04;
		AVX_cmplx_carry_norm_errcheck_X4(s1p04,add1,add2,add3,cy04,bjmodn04,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
		add0 = a + j1 + pfetch_dist + p08;
		AVX_cmplx_carry_norm_errcheck_X4(s1p08,add1,add2,add3,cy08,bjmodn08,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
		add0 += p04;
		AVX_cmplx_carry_norm_errcheck_X4(s1p12,add1,add2,add3,cy12,bjmodn12,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
		add0 = a + j1 + pfetch_dist + p16;
		AVX_cmplx_carry_norm_errcheck_X4(s1p16,add1,add2,add3,cy16,bjmodn16,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
		add0 += p04;
		AVX_cmplx_carry_norm_errcheck_X4(s1p20,add1,add2,add3,cy20,bjmodn20,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  }	// USE_AVX: (8-way or 4-way LOACC) or (4-way HIACC) ?
	 #endif
		i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	 #ifdef USE_ARM_V8_SIMD
	  if(1) {	// No HIACC mode for ARMv8
	 #else
	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining
	 #endif
		/*** wt_re,wi_re,wt_im,wi_im inits. Cf. radix16_main_carry_loop.h for scalar-macro prototyping of this: ***/
		uint32 k0,k1,k2,k3, nwtml;
		l = j & (nwt-1);	nwtml = nwt-l;
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwtml  ];
		sinwtm1 = si[nwtml-1];
		wtl     = wt0[    l  ];
		wtn     = wt0[nwtml  ]*scale;
		wtlp1   = wt0[    l+1];
		wtnm1   = wt0[nwtml-1]*scale;

		ctmp = (struct complex *)half_arr + 24;	// ptr to local storage for the doubled wtl,wtn terms:
		// (j)-data occupy the 8 xmm-sized slots above the 16 used by fixed auxiliary-data, and overwrite these inits:
		ctmp->re = ctmp->im = wtl;		ctmp += 2;
		ctmp->re = ctmp->im = wtn;		ctmp += 2;
		ctmp->re = ctmp->im = wtlp1;	ctmp += 2;
		ctmp->re = ctmp->im = wtnm1;

		l = (j+2) & (nwt-1);	nwtml = nwt-l;
		k0 = n-si[l  ];
		k1 = n-si[l+1];
		k2 = si[nwtml  ];
		k3 = si[nwtml-1];
		wtl     = wt0[    l  ];
		wtn     = wt0[nwtml  ]*scale;
		wtlp1   = wt0[    l+1];
		wtnm1   = wt0[nwtml-1]*scale;

		ctmp = (struct complex *)half_arr + 32;	// (j+2) data start at ctmp + 8
		ctmp->re = ctmp->im = wtl;		ctmp += 2;
		ctmp->re = ctmp->im = wtn;		ctmp += 2;
		ctmp->re = ctmp->im = wtlp1;	ctmp += 2;
		ctmp->re = ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
		// *But*: since the init macro does an on-the-fly version of this between j,j+2 portions, external code co2=co3 must come *after* both ctmp-data octets are inited.
		SSE2_cmplx_carry_fast_wtsinit(add1,add2,add3, bjmodn00, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, k0,k1,k2,k3, sse_bw,sse_n)

		i = (!j);
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy00; tm2 = cy00 + 1; itmp = bjmodn00;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			SSE2_cmplx_carry_fast_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
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

	/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
		i = (!j);
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy00; tm2 = cy00 + 1; itmp = bjmodn00;
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

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy00; tm2 = cy00 + 1; itmp = bjmodn00;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p02,p03, addr);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

	  }	// LOACC or HIACC?

		i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

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

		// Insert a HIACC-version of the cy macro every now and again to reseed weights:
		jt = j1; jp = j2;
		cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],cy00,bjmodn00,0 ,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],cy01,bjmodn01,1 ,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],cy02,bjmodn02,2 ,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],cy03,bjmodn03,3 ,prp_mult);
		jt = j1+p04; jp = j2+p04;
		cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],cy04,bjmodn04,4 ,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],cy05,bjmodn05,5 ,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],cy06,bjmodn06,6 ,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],cy07,bjmodn07,7 ,prp_mult);
		jt = j1+p08; jp = j2+p08;
		cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],cy08,bjmodn08,8 ,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],cy09,bjmodn09,9 ,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],cy10,bjmodn10,10,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],cy11,bjmodn11,11,prp_mult);
		jt = j1+p08+p04; jp = j2+p08+p04;
		cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],cy12,bjmodn12,12,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],cy13,bjmodn13,13,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],cy14,bjmodn14,14,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],cy15,bjmodn15,15,prp_mult);
		jt = j1+p16; jp = j2+p16;
		cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],cy16,bjmodn16,16,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],cy17,bjmodn17,17,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],cy18,bjmodn18,18,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],cy19,bjmodn19,19,prp_mult);
		jt = j1+p16+p04; jp = j2+p16+p04;
		cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],cy20,bjmodn20,20,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],cy21,bjmodn21,21,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],cy22,bjmodn22,22,prp_mult);
		cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],cy23,bjmodn23,23,prp_mult);

	  } else {	// HiACC:

		/*...set0 is slightly different from others:	*/
		jt = j1; jp = j2;
		cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],cy00,bjmodn00,0 ,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p01],a[jp+p01],cy01,bjmodn01,1 ,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p02],a[jp+p02],cy02,bjmodn02,2 ,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p03],a[jp+p03],cy03,bjmodn03,3 ,prp_mult);
		jt = j1+p04; jp = j2+p04;
		cmplx_carry_norm_errcheck (a[jt    ],a[jp    ],cy04,bjmodn04,4 ,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p01],a[jp+p01],cy05,bjmodn05,5 ,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p02],a[jp+p02],cy06,bjmodn06,6 ,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p03],a[jp+p03],cy07,bjmodn07,7 ,prp_mult);
		jt = j1+p08; jp = j2+p08;
		cmplx_carry_norm_errcheck (a[jt    ],a[jp    ],cy08,bjmodn08,8 ,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p01],a[jp+p01],cy09,bjmodn09,9 ,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p02],a[jp+p02],cy10,bjmodn10,10,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p03],a[jp+p03],cy11,bjmodn11,11,prp_mult);
		jt = j1+p08+p04; jp = j2+p08+p04;
		cmplx_carry_norm_errcheck (a[jt    ],a[jp    ],cy12,bjmodn12,12,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p01],a[jp+p01],cy13,bjmodn13,13,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p02],a[jp+p02],cy14,bjmodn14,14,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p03],a[jp+p03],cy15,bjmodn15,15,prp_mult);
		jt = j1+p16; jp = j2+p16;
		cmplx_carry_norm_errcheck (a[jt    ],a[jp    ],cy16,bjmodn16,16,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p01],a[jp+p01],cy17,bjmodn17,17,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p02],a[jp+p02],cy18,bjmodn18,18,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p03],a[jp+p03],cy19,bjmodn19,19,prp_mult);
		jt = j1+p16+p04; jp = j2+p16+p04;
		cmplx_carry_norm_errcheck (a[jt    ],a[jp    ],cy20,bjmodn20,20,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p01],a[jp+p01],cy21,bjmodn21,21,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p02],a[jp+p02],cy22,bjmodn22,22,prp_mult);
		cmplx_carry_norm_errcheck (a[jt+p03],a[jp+p03],cy23,bjmodn23,23,prp_mult);

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
		co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			   and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

/*...The radix-24 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #if USE_SMALL_MACROS

	   #ifdef USE_ARM_V8_SIMD
		const uint32 OFF1 = 0xa0;
		const uint32 OFF2 = 0x40;
		const uint32 OFF3 = 0xe0;
		const uint32 OFF4 = 0x80;
		const uint32 OFF5 = 0x20;
		const uint32 OFF6 = 0xc0;
		const uint32 OFF7 = 0x60;
	   #elif defined(USE_AVX)
		#define OFF1	2*0xa0
		#define OFF2	2*0x40
		#define OFF3	2*0xe0
		#define OFF4	2*0x80
		#define OFF5	2*0x20
		#define OFF6	2*0xc0
		#define OFF7	2*0x60
	   #else
		#define OFF1	0xa0
		#define OFF2	0x40
		#define OFF3	0xe0
		#define OFF4	0x80
		#define OFF5	0x20
		#define OFF6	0xc0
		#define OFF7	0x60
	   #endif

	   #if 0
		SSE2_RADIX_03_DFT(s1p00,s1p16,s1p08,cc3,s1p00,s1p08,s1p16);
		SSE2_RADIX_03_DFT(s1p09,s1p01,s1p17,cc3,s1p01,s1p09,s1p17);
		SSE2_RADIX_03_DFT(s1p18,s1p10,s1p02,cc3,s1p02,s1p10,s1p18);
		SSE2_RADIX_03_DFT(s1p03,s1p19,s1p11,cc3,s1p03,s1p11,s1p19);
		SSE2_RADIX_03_DFT(s1p12,s1p04,s1p20,cc3,s1p04,s1p12,s1p20);
		SSE2_RADIX_03_DFT(s1p21,s1p13,s1p05,cc3,s1p05,s1p13,s1p21);
		SSE2_RADIX_03_DFT(s1p06,s1p22,s1p14,cc3,s1p06,s1p14,s1p22);
		SSE2_RADIX_03_DFT(s1p15,s1p07,s1p23,cc3,s1p07,s1p15,s1p23);
	   #else
		SSE2_RADIX_03_DFT_X2(cc3, s1p00,s1p16,s1p08, s1p00,s1p08,s1p16,  s1p09,s1p01,s1p17, s1p01,s1p09,s1p17);
		SSE2_RADIX_03_DFT_X2(cc3, s1p18,s1p10,s1p02, s1p02,s1p10,s1p18,  s1p03,s1p19,s1p11, s1p03,s1p11,s1p19);
		SSE2_RADIX_03_DFT_X2(cc3, s1p12,s1p04,s1p20, s1p04,s1p12,s1p20,  s1p21,s1p13,s1p05, s1p05,s1p13,s1p21);
		SSE2_RADIX_03_DFT_X2(cc3, s1p06,s1p22,s1p14, s1p06,s1p14,s1p22,  s1p15,s1p07,s1p23, s1p07,s1p15,s1p23);
	   #endif

		// SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p00 + 0x[0a4e82c6]0, o[0-7] = add0 + p[01235476])
		add0 = &a[j1];
		add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p05; add5 = add0+p04; add6 = add0+p07; add7 = add0+p06;
	#ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(s1p00,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	#else
		SSE2_RADIX8_DIF_0TWIDDLE(s1p00,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	#endif
		// SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p16 + 0x[0a4e82c6]0, o[0-7] = add8 + p[54762310])
		add7 = &a[j1] + p08;
		add0 = add7+p05; add1 = add7+p04; add2 = add7+p07; add3 = add7+p06; add4 = add7+p02; add5 = add7+p03; add6 = add7+p01;
	#ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(s1p16,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	#else
		SSE2_RADIX8_DIF_0TWIDDLE(s1p16,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	#endif
		// SSE2_RADIX8_DIF_0TWIDDLE( i[0-7] = s1p08 + 0x[0a4e82c6]0, o[0-7] = add16+ p[23107645])
		add3 = &a[j1] + p16;
		add0 = add3+p02; add1 = add3+p03; add2 = add3+p01; add4 = add3+p07; add5 = add3+p06; add6 = add3+p04; add7 = add3+p05;
	#ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(s1p08,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	#else
		SSE2_RADIX8_DIF_0TWIDDLE(s1p08,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	#endif

	   #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
		#undef OFF5
		#undef OFF6
		#undef OFF7
	   #endif

	  #else	// USE_SMALL_MACROS = False:

		add0 = &a[j1    ];
		SSE2_RADIX24_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p16,s1p00,isrt2,cc3);

	  #endif

	#else

	  #if PFETCH
		add0 = &a[j1];
		prefetch_p_doubles(add0);
	  #endif

	/*...The radix-24 DIF pass is here:	*/

	/* EWM: 10/21/04: We swap the following outputs of the radix-3 transforms: {1,9,17}<=>{5,13,21}, {3,11,19}<=>{7,15,23}, so that the indexing
					  winds up being in-place. This allows us to properly re-use the ajp1 variables in the carry-pass version of this routine.
	*/

	/*...gather the needed data (24 64-bit complex, i.e. 48 64-bit reals) and do 8 in-place radix-3 transforms...*/
						 /*                        inputs                               */ /*             intermediates                 */ /*                 outputs                   */
	  #if PFETCH
		RADIX_03_DFT_PFETCH(s,c3m1,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08],t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,p01);	jt = j1+p05; jp = j2+p05;
		RADIX_03_DFT_PFETCH(s,c3m1,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],t01,t02,t03,t04,t05,t06,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,p02);	jt = j1+p02; jp = j2+p02;
		RADIX_03_DFT_PFETCH(s,c3m1,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],t01,t02,t03,t04,t05,t06,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,p03);	jt = j1+p07; jp = j2+p07;
		RADIX_03_DFT_PFETCH(s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,p04);	jt = j1+p04; jp = j2+p04;
		RADIX_03_DFT_PFETCH(s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,p05);	jt = j1+p01; jp = j2+p01;
		RADIX_03_DFT_PFETCH(s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,p06);	jt = j1+p06; jp = j2+p06;
		RADIX_03_DFT_PFETCH(s,c3m1,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,p07);	jt = j1+p03; jp = j2+p03;
		RADIX_03_DFT_PFETCH(s,c3m1,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,p08);
	  #else
		RADIX_03_DFT       (s,c3m1,a[j1    ],a[j2    ],a[j1+p16],a[j2+p16],a[j1+p08],a[j2+p08],t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i);	jt = j1+p05; jp = j2+p05;
		RADIX_03_DFT       (s,c3m1,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],t01,t02,t03,t04,t05,t06,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i);	jt = j1+p02; jp = j2+p02;
		RADIX_03_DFT       (s,c3m1,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],t01,t02,t03,t04,t05,t06,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i);	jt = j1+p07; jp = j2+p07;
		RADIX_03_DFT       (s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i);	jt = j1+p04; jp = j2+p04;
		RADIX_03_DFT       (s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i);	jt = j1+p01; jp = j2+p01;
		RADIX_03_DFT       (s,c3m1,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],t01,t02,t03,t04,t05,t06,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i);	jt = j1+p06; jp = j2+p06;
		RADIX_03_DFT       (s,c3m1,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i);	jt = j1+p03; jp = j2+p03;
		RADIX_03_DFT       (s,c3m1,a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i);
	  #endif

	/*...and now do 3 radix-8 transforms:	*/
						 /*                                                          inputs                                                           */ /*                 intermediates                    */ /*                 outputs                   */
	  #if PFETCH
		RADIX_08_DIF_PFETCH(a1p00r,a1p00i,a1p05r,a1p05i,a1p02r,a1p02i,a1p07r,a1p07i,a1p04r,a1p04i,a1p01r,a1p01i,a1p06r,a1p06i,a1p03r,a1p03i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],rt,it,p08+p01,p08+p02,p08+p03,p08+p04,p08+p05);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIF_PFETCH(a1p08r,a1p08i,a1p13r,a1p13i,a1p10r,a1p10i,a1p15r,a1p15i,a1p12r,a1p12i,a1p09r,a1p09i,a1p14r,a1p14i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it,p08+p06,p08+p07,p16    ,p16+p01,p16+p02);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF_PFETCH(a1p16r,a1p16i,a1p21r,a1p21i,a1p18r,a1p18i,a1p23r,a1p23i,a1p20r,a1p20i,a1p17r,a1p17i,a1p22r,a1p22i,a1p19r,a1p19i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,p16+p03,p16+p04,p16+p05,p16+p06,p16+p07);
	  #else
		RADIX_08_DIF       (a1p00r,a1p00i,a1p05r,a1p05i,a1p02r,a1p02i,a1p07r,a1p07i,a1p04r,a1p04i,a1p01r,a1p01i,a1p06r,a1p06i,a1p03r,a1p03i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIF       (a1p08r,a1p08i,a1p13r,a1p13i,a1p10r,a1p10i,a1p15r,a1p15i,a1p12r,a1p12i,a1p09r,a1p09i,a1p14r,a1p14i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF       (a1p16r,a1p16i,a1p21r,a1p21i,a1p18r,a1p18i,a1p23r,a1p23i,a1p20r,a1p20i,a1p17r,a1p17i,a1p22r,a1p22i,a1p19r,a1p19i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
	  #endif

	#endif	/* USE_SSE2 */
	}

	jstart += nwt;
	jhi    += nwt;
	col += RADIX;
	co3 -= RADIX;

}	/* end for(int k=1; k <= khi; k++) */
