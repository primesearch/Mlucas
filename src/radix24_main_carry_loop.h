/*******************************************************************************
*                                                                              *
*   (C) 1997-2016 by Ernst W. Mayer.                                           *
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
	/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
	Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
	but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
	(But only ever need to explicitly do this in debug mode).
	*/
	for(j = jstart; j < jhi; j += stride)
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-24 DIT pass is here:	*/
	/*...gather the needed data (24 64-bit complex, i.e. 48 64-bit reals) and do 3 radix-8 transforms,	*/

	#ifdef USE_SSE2

		add0 = &a[j1    ];
		SSE2_RADIX24_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p16,s1p00r,isrt2,cc3);

	#else

		RADIX_08_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_08_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,rt,it);

	/*...and now do 8 in-place radix-3 transforms.	*/

		RADIX_03_DFT(s,c3m1,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i);
		RADIX_03_DFT(s,c3m1,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,t01,t02,t03,t04,t05,t06,a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i);
		RADIX_03_DFT(s,c3m1,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,t01,t02,t03,t04,t05,t06,a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i);
		RADIX_03_DFT(s,c3m1,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i);
		RADIX_03_DFT(s,c3m1,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,t01,t02,t03,t04,t05,t06,a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i);
		RADIX_03_DFT(s,c3m1,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,t01,t02,t03,t04,t05,t06,a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i);
		RADIX_03_DFT(s,c3m1,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i);
		RADIX_03_DFT(s,c3m1,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,t01,t02,t03,t04,t05,t06,a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i);

	#endif	/* USE_SSE2 */

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 24 separate blocks of the A-array, we need 24 separate carries.	*/

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

	 #ifdef LOACC

		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn00, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

		i = (!j);

		// Each carry macro call also processes 8 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_fast_errcheck_X8(s1p00r, cy00,cy04, bjmodn00,bjmodn04, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04); i = 0;
		add0 = a + j1 + pfetch_dist + p08;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_fast_errcheck_X8(s1p08r, cy08,cy12, bjmodn08,bjmodn12, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04);
		add0 = a + j1 + pfetch_dist + p16;
		AVX_cmplx_carry_fast_errcheck_X8(s1p16r, cy16,cy20, bjmodn16,bjmodn20, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04);

	 #else	// USE_AVX: Hi-accuracy 4-way carry is the default:

		// Each AVX carry macro call also processes 4 prefetches of main-array data
		i = (!j);
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_norm_errcheck_X4(s1p00r,add1,add2,add3,cy00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03); i = 0;
		add0 += p04;
		AVX_cmplx_carry_norm_errcheck_X4(s1p04r,add1,add2,add3,cy04,bjmodn04,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		add0 = a + j1 + pfetch_dist + p08;
		AVX_cmplx_carry_norm_errcheck_X4(s1p08r,add1,add2,add3,cy08,bjmodn08,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		add0 += p04;
		AVX_cmplx_carry_norm_errcheck_X4(s1p12r,add1,add2,add3,cy12,bjmodn12,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		add0 = a + j1 + pfetch_dist + p16;
		AVX_cmplx_carry_norm_errcheck_X4(s1p16r,add1,add2,add3,cy16,bjmodn16,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		add0 += p04;
		AVX_cmplx_carry_norm_errcheck_X4(s1p20r,add1,add2,add3,cy20,bjmodn20,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	 #endif	// USE_AVX: (8-way or 4-way LOACC) or (4-way HIACC) ?

		i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	  #ifdef LOACC

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
		tm1 = s1p00r; tmp = cy00; tm2 = cy00 + 1; itmp = bjmodn00;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			SSE2_cmplx_carry_fast_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
		}

	  #else	// Hi-accuracy is the default:

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

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
		i = (!j);
		tm1 = s1p00r; tmp = cy00; tm2 = cy00 + 1; itmp = bjmodn00;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01);
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

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00r; tmp = cy00; tm2 = cy00 + 1; itmp = bjmodn00;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p02,p03);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

	  #endif	// LOACC or HIACC?

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
	   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00,0 );
		cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy01,bjmodn01,1 );
		cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy02,bjmodn02,2 );
		cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy03,bjmodn03,3 );
		cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy04,bjmodn04,4 );
		cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy05,bjmodn05,5 );
		cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy06,bjmodn06,6 );
		cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy07,bjmodn07,7 );
		cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy08,bjmodn08,8 );
		cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy09,bjmodn09,9 );
		cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy10,bjmodn10,10);
		cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy11,bjmodn11,11);
		cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy12,bjmodn12,12);
		cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy13,bjmodn13,13);
		cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy14,bjmodn14,14);
		cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy15,bjmodn15,15);
		cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy16,bjmodn16,16);
		cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy17,bjmodn17,17);
		cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy18,bjmodn18,18);
		cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy19,bjmodn19,19);
		cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy20,bjmodn20,20);
		cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy21,bjmodn21,21);
		cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy22,bjmodn22,22);
		cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy23,bjmodn23,23);

		i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
		co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			   and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

/*...The radix-24 DIF pass is here:	*/

	#ifdef USE_SSE2

		add0 = &a[j1    ];
		SSE2_RADIX24_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p16,s1p00r,isrt2,cc3);

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
		RADIX_03_DFT_PFETCH(s,c3m1,a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,p01);
		RADIX_03_DFT_PFETCH(s,c3m1,a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,t01,t02,t03,t04,t05,t06,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i,p02);
		RADIX_03_DFT_PFETCH(s,c3m1,a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,t01,t02,t03,t04,t05,t06,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,p03);
		RADIX_03_DFT_PFETCH(s,c3m1,a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,t01,t02,t03,t04,t05,t06,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,p04);
		RADIX_03_DFT_PFETCH(s,c3m1,a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,t01,t02,t03,t04,t05,t06,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,p05);
		RADIX_03_DFT_PFETCH(s,c3m1,a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,t01,t02,t03,t04,t05,t06,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,p06);
		RADIX_03_DFT_PFETCH(s,c3m1,a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i,p07);
		RADIX_03_DFT_PFETCH(s,c3m1,a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,p08);
	  #else
		RADIX_03_DFT       (s,c3m1,a1p00r,a1p00i,a1p16r,a1p16i,a1p08r,a1p08i,t01,t02,t03,t04,t05,t06,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i);
		RADIX_03_DFT       (s,c3m1,a1p21r,a1p21i,a1p13r,a1p13i,a1p05r,a1p05i,t01,t02,t03,t04,t05,t06,a1p05r,a1p05i,a1p13r,a1p13i,a1p21r,a1p21i);
		RADIX_03_DFT       (s,c3m1,a1p18r,a1p18i,a1p10r,a1p10i,a1p02r,a1p02i,t01,t02,t03,t04,t05,t06,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i);
		RADIX_03_DFT       (s,c3m1,a1p15r,a1p15i,a1p07r,a1p07i,a1p23r,a1p23i,t01,t02,t03,t04,t05,t06,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i);
		RADIX_03_DFT       (s,c3m1,a1p12r,a1p12i,a1p04r,a1p04i,a1p20r,a1p20i,t01,t02,t03,t04,t05,t06,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i);
		RADIX_03_DFT       (s,c3m1,a1p09r,a1p09i,a1p01r,a1p01i,a1p17r,a1p17i,t01,t02,t03,t04,t05,t06,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i);
		RADIX_03_DFT       (s,c3m1,a1p06r,a1p06i,a1p22r,a1p22i,a1p14r,a1p14i,t01,t02,t03,t04,t05,t06,a1p06r,a1p06i,a1p14r,a1p14i,a1p22r,a1p22i);
		RADIX_03_DFT       (s,c3m1,a1p03r,a1p03i,a1p19r,a1p19i,a1p11r,a1p11i,t01,t02,t03,t04,t05,t06,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i);
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

}	/* end for(k=1; k <= khi; k++) */

