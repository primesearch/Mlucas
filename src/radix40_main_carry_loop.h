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
	for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;
	/*
	!...gather the needed data (40 64-bit complex, i.e. 80 64-bit reals) and do a radix-40 DIT transform...
	*/
	#ifdef USE_SSE2

	  #if defined(COMPILER_TYPE_MSVC) || !GCC_ASM_FULL_INLINE

		add0 = &a[j1    ]; add1 = add0+p01; add2 = add0+p03; add3 = add0+p02; add4 = add0+p07; add5 = add0+p06; add6 = add0+p05; add7 = add0+p04;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00, isrt2)
	   #endif

		add3 = &a[j1+p16]; add0 = add3+p03; add1 = add3+p02; add2 = add3+p01; add4 = add3+p05; add5 = add3+p04; add6 = add3+p06; add7 = add3+p07;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r16, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r16, isrt2)
	   #endif

		add1 = &a[j1+p32]; add0 = add1+p01; add2 = add1+p02; add3 = add1+p03; add4 = add1+p06; add5 = add1+p07; add6 = add1+p04; add7 = add1+p05;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r32, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r32, isrt2)
	   #endif

		add6 = &a[j1+p08]; add0 = add6+p06; add1 = add6+p07; add2 = add6+p04; add3 = add6+p05; add4 = add6+p02; add5 = add6+p03; add7 = add6+p01;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r48, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r48, isrt2)
	   #endif

		add4 = &a[j1+p24]; add0 = add4+p04; add1 = add4+p05; add2 = add4+p07; add3 = add4+p06; add5 = add4+p01; add6 = add4+p03; add7 = add4+p02;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r64, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r64, isrt2)
	   #endif

	   #if OS_BITS == 64
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, r00,r16,r32,r48,r64, s1p00r,s1p16r,s1p32r,s1p08r,s1p24r,\
										r02,r18,r34,r50,r66, s1p25r,s1p01r,s1p17r,s1p33r,s1p09r)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, r04,r20,r36,r52,r68, s1p10r,s1p26r,s1p02r,s1p18r,s1p34r,\
										r06,r22,r38,r54,r70, s1p35r,s1p11r,s1p27r,s1p03r,s1p19r)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, r08,r24,r40,r56,r72, s1p20r,s1p36r,s1p12r,s1p28r,s1p04r,\
										r10,r26,r42,r58,r74, s1p05r,s1p21r,s1p37r,s1p13r,s1p29r)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, r12,r28,r44,r60,r76, s1p30r,s1p06r,s1p22r,s1p38r,s1p14r,\
										r14,r30,r46,r62,r78, s1p15r,s1p31r,s1p07r,s1p23r,s1p39r)
	   #else
		SSE2_RADIX_05_DFT_0TWIDDLE(r00,r16,r32,r48,r64,xcc1,s1p00r,s1p16r,s1p32r,s1p08r,s1p24r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r02,r18,r34,r50,r66,xcc1,s1p25r,s1p01r,s1p17r,s1p33r,s1p09r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r04,r20,r36,r52,r68,xcc1,s1p10r,s1p26r,s1p02r,s1p18r,s1p34r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r06,r22,r38,r54,r70,xcc1,s1p35r,s1p11r,s1p27r,s1p03r,s1p19r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r08,r24,r40,r56,r72,xcc1,s1p20r,s1p36r,s1p12r,s1p28r,s1p04r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r10,r26,r42,r58,r74,xcc1,s1p05r,s1p21r,s1p37r,s1p13r,s1p29r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r12,r28,r44,r60,r76,xcc1,s1p30r,s1p06r,s1p22r,s1p38r,s1p14r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r14,r30,r46,r62,r78,xcc1,s1p15r,s1p31r,s1p07r,s1p23r,s1p39r)
	   #endif

	  #else	/* GCC-style fully-inlined ASM (64-bit only); */

		add0 = &a[j1    ];
		SSE2_RADIX40_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p16,p24,p32,r00, xcc1,isrt2, s1p00r,s1p08r,s1p16r,s1p24r,s1p32r);

	  #endif

	#else

	/*...gather the needed data (40 64-bit complex) and do 5 radix-8 transforms,	*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_08_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it)
	/*...and now do 8 radix-5 transforms.	*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],rt,it);	tptr += 5;
		jt = j1+p01; jp = j2+p01;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],rt,it);	tptr += 5;
		jt = j1+p02; jp = j2+p02;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it);	tptr += 5;
		jt = j1+p03; jp = j2+p03;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);	tptr += 5;
		jt = j1+p04; jp = j2+p04;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],rt,it);	tptr += 5;
		jt = j1+p05; jp = j2+p05;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],rt,it);	tptr += 5;
		jt = j1+p06; jp = j2+p06;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],rt,it);	tptr += 5;
		jt = j1+p07; jp = j2+p07;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it)

	#endif

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need RADIX separate carries.	*/

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
		tm1 = s1p00r; tmp = cy; itmp = bjmodn;
		// Each AVX carry macro call also processes 4 prefetches of main-array data
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_norm_errcheck0_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		tm1 += 8; tmp += 1; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];
			AVX_cmplx_carry_norm_errcheck1_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
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

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00r; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
		add0 = a + j1 + pfetch_dist;
		SSE2_cmplx_carry_norm_errcheck0_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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

		tm1 = s1p00r; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p02,p03);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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
	   cmplx_carry_norm_errcheck0(a[j1    ],a[j2    ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p01],a[j2+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p02],a[j2+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p03],a[j2+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,p08,...
			cmplx_carry_norm_errcheck(a[jt    ],a[jp    ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p01],a[jp+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p02],a[jp+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p03],a[jp+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

/*...The radix-40 DIF pass is here:	*/

	#ifndef USE_SSE2

	/*...gather the needed data (40 64-bit complex) and do 8 radix-5 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p03; jp = j2+p03;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p06; jp = j2+p06;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p01; jp = j2+p01;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p04; jp = j2+p04;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p07; jp = j2+p07;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p02; jp = j2+p02;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p05; jp = j2+p05;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);

	/*...and now do 5 radix-8 transforms, ***SWAPPING THE ORDER OF THE FINAL 2*** to undo the last-2-outputs-swap in the RADIX_05_DFT macro:	*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

	#else

	  #if defined(COMPILER_TYPE_MSVC) || !GCC_ASM_FULL_INLINE

	/*...gather the needed data (40 64-bit complex) and do 8 radix-5 transforms...*/
	   #if OS_BITS == 64
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, s1p00r,s1p32r,s1p24r,s1p16r,s1p08r, r00,r16,r32,r48,r64,\
										s1p35r,s1p27r,s1p19r,s1p11r,s1p03r, r02,r18,r34,r50,r66)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, s1p30r,s1p22r,s1p14r,s1p06r,s1p38r, r04,r20,r36,r52,r68,\
										s1p25r,s1p17r,s1p09r,s1p01r,s1p33r, r06,r22,r38,r54,r70)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, s1p20r,s1p12r,s1p04r,s1p36r,s1p28r, r08,r24,r40,r56,r72,\
										s1p15r,s1p07r,s1p39r,s1p31r,s1p23r, r10,r26,r42,r58,r74)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, s1p10r,s1p02r,s1p34r,s1p26r,s1p18r, r12,r28,r44,r60,r76,\
										s1p05r,s1p37r,s1p29r,s1p21r,s1p13r, r14,r30,r46,r62,r78)
	   #else
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p00r,s1p32r,s1p24r,s1p16r,s1p08r,xcc1,r00,r16,r32,r48,r64)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p35r,s1p27r,s1p19r,s1p11r,s1p03r,xcc1,r02,r18,r34,r50,r66)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p30r,s1p22r,s1p14r,s1p06r,s1p38r,xcc1,r04,r20,r36,r52,r68)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p25r,s1p17r,s1p09r,s1p01r,s1p33r,xcc1,r06,r22,r38,r54,r70)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p20r,s1p12r,s1p04r,s1p36r,s1p28r,xcc1,r08,r24,r40,r56,r72)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p15r,s1p07r,s1p39r,s1p31r,s1p23r,xcc1,r10,r26,r42,r58,r74)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p10r,s1p02r,s1p34r,s1p26r,s1p18r,xcc1,r12,r28,r44,r60,r76)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p05r,s1p37r,s1p29r,s1p21r,s1p13r,xcc1,r14,r30,r46,r62,r78)
	   #endif

	/*...and now do 5 radix-8 transforms, swapping the t[48+i] <--> t[64+i] pairs to undo the last-2-outputs-swap in the RADIX_05_DFT macro:	*/
	   #ifdef USE_AVX
		#define OFF1	0x40
		#define OFF2	0x80
		#define OFF3	0xc0
		#define OFF4	0x100
		#define OFF5	0x140
		#define OFF6	0x180
		#define OFF7	0x1c0
	   #else
		#define OFF1	0x20
		#define OFF2	0x40
		#define OFF3	0x60
		#define OFF4	0x80
		#define OFF5	0xa0
		#define OFF6	0xc0
		#define OFF7	0xe0
	   #endif

		add0 = &a[j1    ]; add1 = add0+p01; add2 = add0+p03; add3 = add0+p02; add4 = add0+p06; add5 = add0+p07; add6 = add0+p04; add7 = add0+p05;
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r00,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	  #else
		SSE2_RADIX8_DIF_0TWIDDLE(r00,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	  #endif

		add1 = &a[j1+p32]; add0 = add1+p01; add2 = add1+p02; add3 = add1+p03; add4 = add1+p07; add5 = add1+p06; add6 = add1+p05; add7 = add1+p04;
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r16,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	  #else
		SSE2_RADIX8_DIF_0TWIDDLE(r16,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	  #endif

		add5 = &a[j1+p24]; add0 = add5+p04; add1 = add5+p05; add2 = add5+p07; add3 = add5+p06; add4 = add5+p01; add6 = add5+p02; add7 = add5+p03;
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r32,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	  #else
		SSE2_RADIX8_DIF_0TWIDDLE(r32,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	  #endif

		add3 = &a[j1+p16]; add0 = add3+p03; add1 = add3+p02; add2 = add3+p01; add4 = add3+p04; add5 = add3+p05; add6 = add3+p07; add7 = add3+p06;
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r48,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	  #else
		SSE2_RADIX8_DIF_0TWIDDLE(r48,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	  #endif

		add7 = &a[j1+p08]; add0 = add7+p06; add1 = add7+p07; add2 = add7+p04; add3 = add7+p05; add4 = add7+p03; add5 = add7+p02; add6 = add7+p01;
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r64,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	  #else
		SSE2_RADIX8_DIF_0TWIDDLE(r64,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	  #endif

		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
		#undef OFF5
		#undef OFF6
		#undef OFF7

	  #else	/* GCC-style inline ASM: */

		add0 = &a[j1    ];
		SSE2_RADIX40_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32,r00, xcc1,isrt2, s1p00r,s1p08r,s1p16r,s1p24r,s1p32r);

	  #endif

	#endif	/* #ifdef USE_SSE2 */
	}

	jstart += nwt;
	jhi    += nwt;

	col += RADIX;
	co3 -= RADIX;
}	/* end for(k=1; k <= khi; k++) */
