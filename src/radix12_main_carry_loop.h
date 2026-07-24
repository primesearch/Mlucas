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

	/*...The radix-12 DIT pass is here:	*/
	#ifdef USE_SSE2

	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF = 0x60;
	  #elif defined(USE_AVX512)
		#define OFF	0x60*4
	  #elif defined(USE_AVX)
		#define OFF	0x60*2
	  #else
		#define OFF	0x60
	  #endif
		add0 = &a[j1   ];	add1 = add0+p1; add2 = add0+p2; add3 = add0+p3;	tmp = r00;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0,add1,add3,add2, tmp, OFF);
		add0 = &a[j1+p8];	add1 = add0+p1; add2 = add0+p2; add3 = add0+p3;	tmp += 2;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add1,add0,add2,add3, tmp, OFF);
		add0 = &a[j1+p4];	add1 = add0+p1; add2 = add0+p2; add3 = add0+p3;	tmp += 2;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add2,add3,add0,add1, tmp, OFF);
		// ^^^^^^^^^^^^ add0 += p4-p8 on above line gives wrong ptr-value
		// Wide-stride outouts for above 4-DFTs means contig-stride inputs for the ensuing 3-DFTs:
		SSE2_RADIX_03_DFT(r00,r01,r02,cc3,s1p00,s1p04,s1p08);
		SSE2_RADIX_03_DFT(r03,r04,r05,cc3,s1p03,s1p07,s1p11);
		SSE2_RADIX_03_DFT(r06,r07,r08,cc3,s1p06,s1p10,s1p02);
		SSE2_RADIX_03_DFT(r09,r10,r11,cc3,s1p09,s1p01,s1p05);

	#else

		double u1,u2,u3,u4,u5,u6;
		// Do a trio of 4-DFTs; combined DIT input-scramble array = [0,1,3,2|9,8,a,b|6,7,4,5]:
		RADIX_04_DIT(a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p3],a[j2+p3],a[j1+p2],a[j2+p2],t00,t01,t06,t07,t12,t13,t18,t19,rt,it);	jt = j1+p8; jp = j2+p8;
		RADIX_04_DIT(a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],t02,t03,t08,t09,t14,t15,t20,t21,rt,it);	jt = j1+p4; jp = j2+p4;
		RADIX_04_DIT(a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],t04,t05,t10,t11,t16,t17,t22,t23,rt,it);
		/* And now do a quarter of 3-DFTs. Required output index permutation: [0,1,2,3,4,5,6,7,8,9,a,b] => [0,3,6,9,4,7,a,1,8,b,2,5],
		i.e. in terms of our trial-output ordering into trios, [0,4,8|1,5,9|2,6,a|3,7,b] => [0,4,8|3,7,b|6,a,2|9,1,5]:
		*/
		RADIX_03_DFT(s,c3m1,t00,t01,t02,t03,t04,t05,u1,u2,u3,u4,u5,u6,a[j1   ],a[j2   ],a[j1+p4],a[j2+p4],a[j1+p8],a[j2+p8]);	jt = j1+p3; jp = j2+p3;
		RADIX_03_DFT(s,c3m1,t06,t07,t08,t09,t10,t11,u1,u2,u3,u4,u5,u6,a[jt   ],a[jp   ],a[jt+p4],a[jp+p4],a[jt+p8],a[jp+p8]);	jt = j1+p2; jp = j2+p2;
		RADIX_03_DFT(s,c3m1,t12,t13,t14,t15,t16,t17,u1,u2,u3,u4,u5,u6,a[jt+p4],a[jp+p4],a[jt+p8],a[jp+p8],a[jt   ],a[jp   ]);	jt = j1+p1; jp = j2+p1;
		RADIX_03_DFT(s,c3m1,t18,t19,t20,t21,t22,t23,u1,u2,u3,u4,u5,u6,a[jt+p8],a[jp+p8],a[jt   ],a[jp   ],a[jt+p4],a[jp+p4]);


	#endif	/* USE_SSE2 */

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 12 separate blocks of the A-array, we need 12 separate carries.	*/

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

		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
		AVX_cmplx_carry_fast_wtsinit_X4(add1,add2,add3, bjmodn00, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

		i = (!j);
		addr = &prp_mult;
		// Each carry macro call also processes 4 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_fast_errcheck_X4(s1p00, cy00, bjmodn00, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr); i = 0;
		add0 += p4;
		AVX_cmplx_carry_fast_errcheck_X4(s1p04, cy04, bjmodn04, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr);
		add0 = a + j1 + pfetch_dist + p8;
		AVX_cmplx_carry_fast_errcheck_X4(s1p08, cy08, bjmodn08, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr);

	  } else {	// HIACC:
	 #endif
		// Each AVX carry macro call also processes 4 prefetches of main-array data
		i = (!j);
		addr = &prp_mult;
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_norm_errcheck_X4(s1p00,add1,add2,add3,cy00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr); i = 0;
		add0 += p4;
		AVX_cmplx_carry_norm_errcheck_X4(s1p04,add1,add2,add3,cy04,bjmodn04,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr);
		add0 = a + j1 + pfetch_dist + p8;
		AVX_cmplx_carry_norm_errcheck_X4(s1p08,add1,add2,add3,cy08,bjmodn08,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr);

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
	 #ifndef USE_AVX512
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
			SSE2_cmplx_carry_fast_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr);
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
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1, addr);
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
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p2,p3, addr);
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

		/*...set0 is slightly different from others:	*/
		jt = j1; jp = j2;
	   cmplx_carry_norm_errcheck0(a[jt   ],a[jp   ],cy00,bjmodn00,0 ,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],cy01,bjmodn01,1 ,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],cy02,bjmodn02,2 ,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],cy03,bjmodn03,3 ,prp_mult);
		jt = j1+p4; jp = j2+p4;
		cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],cy04,bjmodn04,4 ,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],cy05,bjmodn05,5 ,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],cy06,bjmodn06,6 ,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],cy07,bjmodn07,7 ,prp_mult);
		jt = j1+p8; jp = j2+p8;
		cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],cy08,bjmodn08,8 ,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],cy09,bjmodn09,9 ,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],cy10,bjmodn10,10,prp_mult);
		cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],cy11,bjmodn11,11,prp_mult);

		i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
		co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			   and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

/*...The radix-12 DIF pass is here:	*/

	#ifdef USE_SSE2

		SSE2_RADIX_03_DFT(s1p00,s1p08,s1p04,cc3,r00,r01,r02);
		SSE2_RADIX_03_DFT(s1p09,s1p05,s1p01,cc3,r03,r04,r05);
		SSE2_RADIX_03_DFT(s1p06,s1p02,s1p10,cc3,r06,r07,r08);
		SSE2_RADIX_03_DFT(s1p03,s1p11,s1p07,cc3,r09,r10,r11);
		// Contig-stride outouts for above 3-DFTs means wide-stride inputs for the ensuing 4-DFTs:
		tmp = r00;	add0 = &a[j1   ]; add1 = add0+p1; add2 = add0+p2; add3 = add0+p3;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0,add1,add2,add3, tmp, OFF)
		tmp += 2;	add0 = &a[j1+p8]; add1 = add0+p1; add2 = add0+p2; add3 = add0+p3;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add1,add0,add3,add2, tmp, OFF)
		tmp += 2;	add0 = &a[j1+p4]; add1 = add0+p1; add2 = add0+p2; add3 = add0+p3;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add2,add3,add1,add0, tmp, OFF)
	  #ifndef USE_ARM_V8_SIMD
		#undef OFF
	  #endif

	#else
		// Twiddleless version requires us to swap inputs [0,4,8|1,5,9|2,6,10|3,7,11] => [0,8,4|9,5,1|6,2,10|3,11,7]:
		RADIX_03_DFT(s,c3m1,a[j1   ],a[j2   ],a[j1+p8],a[j2+p8],a[j1+p4],a[j2+p4],u1,u2,u3,u4,u5,u6,t00,t01,t02,t03,t04,t05);	jt = j1+p1; jp = j2+p1;
		RADIX_03_DFT(s,c3m1,a[jt+p8],a[jp+p8],a[jt+p4],a[jp+p4],a[jt   ],a[jp   ],u1,u2,u3,u4,u5,u6,t06,t07,t08,t09,t10,t11);	jt = j1+p2; jp = j2+p2;
		RADIX_03_DFT(s,c3m1,a[jt+p4],a[jp+p4],a[jt   ],a[jp   ],a[jt+p8],a[jp+p8],u1,u2,u3,u4,u5,u6,t12,t13,t14,t15,t16,t17);	jt = j1+p3; jp = j2+p3;
		RADIX_03_DFT(s,c3m1,a[jt   ],a[jp   ],a[jt+p8],a[jp+p8],a[jt+p4],a[jp+p4],u1,u2,u3,u4,u5,u6,t18,t19,t20,t21,t22,t23);
		// And do a trio of 4-DFTs; Required output index permutation = [0,1,2,3|9,8,b,a|6,7,5,4]:
		RADIX_04_DIF(t00,t01,t06,t07,t12,t13,t18,t19,a[j1   ],a[j2   ],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],rt,it);	jt = j1+p8; jp = j2+p8;
		RADIX_04_DIF(t02,t03,t08,t09,t14,t15,t20,t21,a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],rt,it);	jt = j1+p4; jp = j2+p4;
		RADIX_04_DIF(t04,t05,t10,t11,t16,t17,t22,t23,a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],rt,it);

	#endif	/* USE_SSE2 */
	}

	jstart += nwt;
	jhi    += nwt;
	col += RADIX;
	co3 -= RADIX;

}	/* end for(int k=1; k <= khi; k++) */

