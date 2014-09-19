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

	  #if defined(COMPILER_TYPE_MSVC)

		add0 = &a[j1    ];
		add1 = add0+p01;
		add2 = add0+p03;
		add3 = add0+p02;
		add4 = add0+p07;
		add5 = add0+p06;
		add6 = add0+p05;
		add7 = add0+p04;
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p00r, isrt2)

		add5 = &a[j1+p08];
		add0 = add5+p05;
		add1 = add5+p04;
		add2 = add5+p06;
		add3 = add5+p07;
		add4 = add5+p01;
		add6 = add5+p02;
		add7 = add5+p03;
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p08r, isrt2)

		add2 = &a[j1+p16];
		add0 = add2+p02;
		add1 = add2+p03;
		add3 = add2+p01;
		add4 = add2+p04;
		add5 = add2+p05;
		add6 = add2+p07;
		add7 = add2+p06;
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, s1p16r, isrt2)

		SSE2_RADIX_03_DFT(s1p00r,s1p08r,s1p16r,cc3,s1p00r,s1p16r,s1p08r)
		SSE2_RADIX_03_DFT(s1p01r,s1p09r,s1p17r,cc3,s1p09r,s1p01r,s1p17r)
		SSE2_RADIX_03_DFT(s1p02r,s1p10r,s1p18r,cc3,s1p18r,s1p10r,s1p02r)
		SSE2_RADIX_03_DFT(s1p03r,s1p11r,s1p19r,cc3,s1p03r,s1p19r,s1p11r)
		SSE2_RADIX_03_DFT(s1p04r,s1p12r,s1p20r,cc3,s1p12r,s1p04r,s1p20r)
		SSE2_RADIX_03_DFT(s1p05r,s1p13r,s1p21r,cc3,s1p21r,s1p13r,s1p05r)
		SSE2_RADIX_03_DFT(s1p06r,s1p14r,s1p22r,cc3,s1p06r,s1p22r,s1p14r)
		SSE2_RADIX_03_DFT(s1p07r,s1p15r,s1p23r,cc3,s1p15r,s1p07r,s1p23r)

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

			add0 = &a[j1    ];
			SSE2_RADIX24_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p16,s1p00r,isrt2,cc3);

	  #endif

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

		AVX_cmplx_carry_norm_errcheck0_X4(s1p00r,add1,add2,add3,cy00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p04r,add1,add2,add3,cy04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p08r,add1,add2,add3,cy08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p12r,add1,add2,add3,cy12,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p16r,add1,add2,add3,cy16,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p20r,add1,add2,add3,cy20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

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

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

	  #if defined(COMPILER_TYPE_MSVC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy20,cy22,bjmodn20);
	   #else
		SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy20,cy22,bjmodn20);
	   #endif

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
	   #else
		SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
	   #endif

	  #endif

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

	  #if defined(COMPILER_TYPE_MSVC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy20,cy22,bjmodn20);
	   #else
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy20,cy22,bjmodn20);
	   #endif

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
	   #else
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
	   #endif

	  #endif

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
	   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00   );
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

	  #if defined(COMPILER_TYPE_MSVC)

		SSE2_RADIX_03_DFT(s1p00r,s1p16r,s1p08r,cc3,s1p00r,s1p08r,s1p16r)
		SSE2_RADIX_03_DFT(s1p21r,s1p13r,s1p05r,cc3,s1p05r,s1p13r,s1p21r)
		SSE2_RADIX_03_DFT(s1p18r,s1p10r,s1p02r,cc3,s1p02r,s1p10r,s1p18r)
		SSE2_RADIX_03_DFT(s1p15r,s1p07r,s1p23r,cc3,s1p07r,s1p15r,s1p23r)
		SSE2_RADIX_03_DFT(s1p12r,s1p04r,s1p20r,cc3,s1p04r,s1p12r,s1p20r)
		SSE2_RADIX_03_DFT(s1p09r,s1p01r,s1p17r,cc3,s1p01r,s1p09r,s1p17r)
		SSE2_RADIX_03_DFT(s1p06r,s1p22r,s1p14r,cc3,s1p06r,s1p14r,s1p22r)
		SSE2_RADIX_03_DFT(s1p03r,s1p19r,s1p11r,cc3,s1p03r,s1p11r,s1p19r)

		add0 = &a[j1    ];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
		add4 = add0+p05;
		add5 = add0+p04;
		add6 = add0+p07;
		add7 = add0+p06;
		SSE2_RADIX8_DIF_0TWIDDLE(s1p00r,0xa0,0x40,0xe0,0x80,0x20,0xc0,0x60, add0,add1,add2,add3,add4,add5,add6,add7, isrt2)

		add3 = &a[j1+p16];
		add0 = add3+p02;
		add1 = add3+p03;
		add2 = add3+p01;
		add4 = add3+p07;
		add5 = add3+p06;
		add6 = add3+p04;
		add7 = add3+p05;
		SSE2_RADIX8_DIF_0TWIDDLE(s1p08r,0xa0,0x40,0xe0,0x80,0x20,0xc0,0x60, add0,add1,add2,add3,add4,add5,add6,add7, isrt2)

		add7 = &a[j1+p08];
		add0 = add7+p05;
		add1 = add7+p04;
		add2 = add7+p07;
		add3 = add7+p06;
		add4 = add7+p02;
		add5 = add7+p03;
		add6 = add7+p01;
		SSE2_RADIX8_DIF_0TWIDDLE(s1p16r,0xa0,0x40,0xe0,0x80,0x20,0xc0,0x60, add0,add1,add2,add3,add4,add5,add6,add7, isrt2)

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		add0 = &a[j1    ];
		SSE2_RADIX24_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p16,s1p00r,isrt2,cc3);

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

