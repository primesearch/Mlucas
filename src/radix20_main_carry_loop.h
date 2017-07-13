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

		/*
		!...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do a radix-20 DIT transform...
		*/
#ifdef USE_SSE2

		add0 = &a[j1    ];
		SSE2_RADIX20_DIT_NOTWIDDLE(add0,p01,p04,r00,r10,r20,r30,cc1,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r);

#else
						 /*          inputs           */ /*                                      outputs                                      */
		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00,t01,t10,t11,t20,t21,t30,t31,rt,it);
		RADIX_04_DIT(a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t02,t03,t12,t13,t22,t23,t32,t33,rt,it);
		RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],t04,t05,t14,t15,t24,t25,t34,t35,rt,it);
		RADIX_04_DIT(a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],t06,t07,t16,t17,t26,t27,t36,t37,rt,it);
		RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],t08,t09,t18,t19,t28,t29,t38,t39,rt,it);

	/*...and now do 4 radix-5 transforms...*/
						 /*                                                inputs                                                   */ /*                 outputs                   */
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,rt,it);

#endif

/*...Now do the carries. Since the outputs would
normally be getting dispatched to 20 separate blocks of the A-array, we need 20 separate carries.	*/

	#ifdef USE_AVX	// AVX: can select between carry macros processing 4 and 8 independent carries in LOACC mode:

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

	  #ifdef CARRY_8_WAY
		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn00, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
	  #else
		AVX_cmplx_carry_fast_wtsinit_X4(add1,add2,add3, bjmodn00, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
	  #endif

		i = (!j);

	  #ifdef CARRY_8_WAY

		// Each carry macro call also processes 8 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_fast_errcheck_X8(s1p00r, cy00,cy04, bjmodn00,bjmodn04, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04); i = 0;
		add0 = a + j1 + pfetch_dist + p08;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_fast_errcheck_X8(s1p08r, cy08,cy12, bjmodn08,bjmodn12, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04);
		add0 = a + j1 + pfetch_dist + p16;
		AVX_cmplx_carry_fast_errcheck_X4(s1p16r, cy16,      bjmodn16,          half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);

	  #else	// USE_AVX, LOACC, 4-way carry:
	  	//*** 4-way is Lower accuracy than 8-way, thus default is to use 4-way only for cleanup, e.g. 20 = 2*8 + 4 here ***:

		// Each carry macro call also processes 4 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_fast_errcheck_X4(s1p00r, cy00, bjmodn00, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03); i = 0;
		add0 += p04;
		AVX_cmplx_carry_fast_errcheck_X4(s1p04r, cy04, bjmodn04, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		add0 = a + j1 + pfetch_dist + p08;
		AVX_cmplx_carry_fast_errcheck_X4(s1p08r, cy08, bjmodn08, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		add0 += p04;
		AVX_cmplx_carry_fast_errcheck_X4(s1p12r, cy12, bjmodn12, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		add0 = a + j1 + pfetch_dist + p16;
		AVX_cmplx_carry_fast_errcheck_X4(s1p16r, cy16, bjmodn16, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);

	  #endif	// USE_AVX, LOACC, (8-way or 4-way carry) ?

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

		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
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

	#ifdef LOACC
		// Insert a HIACC-version of the cy macro every now and again to reseed weights:
		cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00,0 );
		cmplx_carry_fast_errcheck (a1p01r,a1p01i,cy01,bjmodn01,1 );
		cmplx_carry_fast_errcheck (a1p02r,a1p02i,cy02,bjmodn02,2 );
		cmplx_carry_fast_errcheck (a1p03r,a1p03i,cy03,bjmodn03,3 );
		cmplx_carry_norm_errcheck0(a1p04r,a1p04i,cy04,bjmodn04,4 );
		cmplx_carry_fast_errcheck (a1p05r,a1p05i,cy05,bjmodn05,5 );
		cmplx_carry_fast_errcheck (a1p06r,a1p06i,cy06,bjmodn06,6 );
		cmplx_carry_fast_errcheck (a1p07r,a1p07i,cy07,bjmodn07,7 );
		cmplx_carry_norm_errcheck0(a1p08r,a1p08i,cy08,bjmodn08,8 );
		cmplx_carry_fast_errcheck (a1p09r,a1p09i,cy09,bjmodn09,9 );
		cmplx_carry_fast_errcheck (a1p10r,a1p10i,cy10,bjmodn10,10);
		cmplx_carry_fast_errcheck (a1p11r,a1p11i,cy11,bjmodn11,11);
		cmplx_carry_norm_errcheck0(a1p12r,a1p12i,cy12,bjmodn12,12);
		cmplx_carry_fast_errcheck (a1p13r,a1p13i,cy13,bjmodn13,13);
		cmplx_carry_fast_errcheck (a1p14r,a1p14i,cy14,bjmodn14,14);
		cmplx_carry_fast_errcheck (a1p15r,a1p15i,cy15,bjmodn15,15);
		cmplx_carry_norm_errcheck0(a1p16r,a1p16i,cy16,bjmodn16,16);
		cmplx_carry_fast_errcheck (a1p17r,a1p17i,cy17,bjmodn17,17);
		cmplx_carry_fast_errcheck (a1p18r,a1p18i,cy18,bjmodn18,18);
		cmplx_carry_fast_errcheck (a1p19r,a1p19i,cy19,bjmodn19,19);
	#else
		/*...set0 is slightly different from others:	*/
		cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00,0 );
		cmplx_carry_norm_errcheck (a1p01r,a1p01i,cy01,bjmodn01,1 );
		cmplx_carry_norm_errcheck (a1p02r,a1p02i,cy02,bjmodn02,2 );
		cmplx_carry_norm_errcheck (a1p03r,a1p03i,cy03,bjmodn03,3 );
		cmplx_carry_norm_errcheck (a1p04r,a1p04i,cy04,bjmodn04,4 );
		cmplx_carry_norm_errcheck (a1p05r,a1p05i,cy05,bjmodn05,5 );
		cmplx_carry_norm_errcheck (a1p06r,a1p06i,cy06,bjmodn06,6 );
		cmplx_carry_norm_errcheck (a1p07r,a1p07i,cy07,bjmodn07,7 );
		cmplx_carry_norm_errcheck (a1p08r,a1p08i,cy08,bjmodn08,8 );
		cmplx_carry_norm_errcheck (a1p09r,a1p09i,cy09,bjmodn09,9 );
		cmplx_carry_norm_errcheck (a1p10r,a1p10i,cy10,bjmodn10,10);
		cmplx_carry_norm_errcheck (a1p11r,a1p11i,cy11,bjmodn11,11);
		cmplx_carry_norm_errcheck (a1p12r,a1p12i,cy12,bjmodn12,12);
		cmplx_carry_norm_errcheck (a1p13r,a1p13i,cy13,bjmodn13,13);
		cmplx_carry_norm_errcheck (a1p14r,a1p14i,cy14,bjmodn14,14);
		cmplx_carry_norm_errcheck (a1p15r,a1p15i,cy15,bjmodn15,15);
		cmplx_carry_norm_errcheck (a1p16r,a1p16i,cy16,bjmodn16,16);
		cmplx_carry_norm_errcheck (a1p17r,a1p17i,cy17,bjmodn17,17);
		cmplx_carry_norm_errcheck (a1p18r,a1p18i,cy18,bjmodn18,18);
		cmplx_carry_norm_errcheck (a1p19r,a1p19i,cy19,bjmodn19,19);
	#endif

		i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
		co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

/*...The radix-20 DIF pass is here:	*/

	#ifdef USE_SSE2

		add0 = &a[j1    ];	// re-init this, because ptr used as a prefetch address in carry step above
		SSE2_RADIX20_DIF_NOTWIDDLE(add0,p01,p04,p08,p16,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r,cc1,r00,r10,r20,r30)

	#else	// Scalar-double mode:

	/*...The radix-20 DIF pass is here:	*/
	  #if PFETCH
		addr = &a[j1];
		prefetch_p_doubles(addr);
	  #endif

	/*...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do 4 radix-5 transforms...*/
						 /*                                                inputs                                                   */ /*                 outputs                   */
	  #if PFETCH																																/*[--y3-] [--y4-] <<<<< swap last 2 outputs to undo swap of these in macro */
		RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it,addr,addp,p01,p02);
		RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it,addr,addp,p03,p04);
		RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it,addr,addp,p05,p06);
		RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it,addr,addp,p07,p08);
	  #else
		RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it);
		RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it);
		RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it);
		RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it);
	  #endif

	/*...and now do 5 radix-4 transforms...*/
						 /*          inputs           */ /*                                      outputs                                      */
	  #if PFETCH
		addp = addr+p09;
		prefetch_p_doubles(addp);

		RADIX_04_DIF_PFETCH(t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it,addr,addp,p10,p11);
		RADIX_04_DIF_PFETCH(t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it,addr,addp,p12,p13);
		RADIX_04_DIF_PFETCH(t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it,addr,addp,p14,p15);
		RADIX_04_DIF_PFETCH(t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it,addr,addp,p16,p17);
		RADIX_04_DIF_PFETCH(t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it,addr,addp,p18,p19);
	  #else
		RADIX_04_DIF       (t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it);
		RADIX_04_DIF       (t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it);
		RADIX_04_DIF       (t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
		RADIX_04_DIF       (t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
		RADIX_04_DIF       (t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it);
	  #endif

	#endif	/* USE_SSE2 */
	}

	jstart += nwt;
	jhi    += nwt;
	col += RADIX;
	co3 -= RADIX;

}	/* end for(k=1; k <= khi; k++) */

