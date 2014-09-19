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

	/*...The radix-48 DIT pass is here:	*/

	#ifdef USE_SSE2

	// Macros use literal ostrides which are [1,2,3,4]-multiples of 1 vector-complex = 0x20 bytes for sse2, 0x40 bytes for avx
	// Offsets: 	00,01,03,02,07,06,05,04,15,14,13,12,11,10,09,08:
		add0 = &a[j1    ];	addf = add0 + p08;
		add1 = add0 + p01;	add8 = addf + p07;
		add2 = add0 + p03;	add9 = addf + p06;
		add3 = add0 + p02;	adda = addf + p05;
		add4 = add0 + p07;	addb = addf + p04;
		add5 = add0 + p06;	addc = addf + p03;
		add6 = add0 + p05;	addd = addf + p02;
		add7 = add0 + p04;	adde = addf + p01;
	  #ifdef USE_AVX
		SSE2_RADIX16_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, isrt2, r00r,0x40,0x80,0xc0,0x100)
	  #else
		SSE2_RADIX16_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, isrt2, r00r,0x20,0x40,0x60,0x80)
	  #endif

	// Offsets: 32+	05,04,06,07,01,00,02,03,09,08,10,11,14,15,12,13:
		add5 = &a[j1+p32];	add9 = add5 + p08;
		add0 = add5 + p05;	add8 = add9 + p01;
		add1 = add5 + p04;	adda = add9 + p02;
		add2 = add5 + p06;	addb = add9 + p03;
		add3 = add5 + p07;	addc = add9 + p06;
		add4 = add5 + p01;	addd = add9 + p07;
		add6 = add5 + p02;	adde = add9 + p04;
		add7 = add5 + p03;	addf = add9 + p05;
	  #ifdef USE_AVX
		SSE2_RADIX16_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, isrt2, r16r,0x40,0x80,0xc0,0x100)
	  #else
		SSE2_RADIX16_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, isrt2, r16r,0x20,0x40,0x60,0x80)
	  #endif

	// Offsets: 16+ 10,11,08,09,12,13,15,14,02,03,00,01,04,05,07,06:
		adda = &a[j1+p16];	add2 = adda + p08;
		add0 = add2 + p02;	add8 = adda + p02;
		add1 = add2 + p03;	add9 = adda + p03;
		add3 = add2 + p01;	addb = adda + p01;
		add4 = add2 + p04;	addc = adda + p04;
		add5 = add2 + p05;	addd = adda + p05;
		add6 = add2 + p07;	adde = adda + p07;
		add7 = add2 + p06;	addf = adda + p06;
	  #ifdef USE_AVX
		SSE2_RADIX16_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, isrt2, r32r,0x40,0x80,0xc0,0x100)
	  #else
		SSE2_RADIX16_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, isrt2, r32r,0x20,0x40,0x60,0x80)
	  #endif

	  #if OS_BITS == 64
		SSE2_RADIX_03_DFT_X2(cc1, r01r,r17r,r33r, s1p15r,s1p31r,s1p47r, r00r,r16r,r32r, s1p00r,s1p16r,s1p32r)
		SSE2_RADIX_03_DFT_X2(cc1, r03r,r19r,r35r, s1p45r,s1p13r,s1p29r, r02r,r18r,r34r, s1p30r,s1p46r,s1p14r)
		SSE2_RADIX_03_DFT_X2(cc1, r05r,r21r,r37r, s1p27r,s1p43r,s1p11r, r04r,r20r,r36r, s1p12r,s1p28r,s1p44r)
		SSE2_RADIX_03_DFT_X2(cc1, r07r,r23r,r39r, s1p09r,s1p25r,s1p41r, r06r,r22r,r38r, s1p42r,s1p10r,s1p26r)
		SSE2_RADIX_03_DFT_X2(cc1, r09r,r25r,r41r, s1p39r,s1p07r,s1p23r, r08r,r24r,r40r, s1p24r,s1p40r,s1p08r)
		SSE2_RADIX_03_DFT_X2(cc1, r11r,r27r,r43r, s1p21r,s1p37r,s1p05r, r10r,r26r,r42r, s1p06r,s1p22r,s1p38r)
		SSE2_RADIX_03_DFT_X2(cc1, r13r,r29r,r45r, s1p03r,s1p19r,s1p35r, r12r,r28r,r44r, s1p36r,s1p04r,s1p20r)
		SSE2_RADIX_03_DFT_X2(cc1, r15r,r31r,r47r, s1p33r,s1p01r,s1p17r, r14r,r30r,r46r, s1p18r,s1p34r,s1p02r)
	  #else
		SSE2_RADIX_03_DFT(r00r,r16r,r32r, cc1, s1p00r,s1p16r,s1p32r)
		SSE2_RADIX_03_DFT(r01r,r17r,r33r, cc1, s1p15r,s1p31r,s1p47r)
		SSE2_RADIX_03_DFT(r02r,r18r,r34r, cc1, s1p30r,s1p46r,s1p14r)
		SSE2_RADIX_03_DFT(r03r,r19r,r35r, cc1, s1p45r,s1p13r,s1p29r)
		SSE2_RADIX_03_DFT(r04r,r20r,r36r, cc1, s1p12r,s1p28r,s1p44r)
		SSE2_RADIX_03_DFT(r05r,r21r,r37r, cc1, s1p27r,s1p43r,s1p11r)
		SSE2_RADIX_03_DFT(r06r,r22r,r38r, cc1, s1p42r,s1p10r,s1p26r)
		SSE2_RADIX_03_DFT(r07r,r23r,r39r, cc1, s1p09r,s1p25r,s1p41r)
		SSE2_RADIX_03_DFT(r08r,r24r,r40r, cc1, s1p24r,s1p40r,s1p08r)
		SSE2_RADIX_03_DFT(r09r,r25r,r41r, cc1, s1p39r,s1p07r,s1p23r)
		SSE2_RADIX_03_DFT(r10r,r26r,r42r, cc1, s1p06r,s1p22r,s1p38r)
		SSE2_RADIX_03_DFT(r11r,r27r,r43r, cc1, s1p21r,s1p37r,s1p05r)
		SSE2_RADIX_03_DFT(r12r,r28r,r44r, cc1, s1p36r,s1p04r,s1p20r)
		SSE2_RADIX_03_DFT(r13r,r29r,r45r, cc1, s1p03r,s1p19r,s1p35r)
		SSE2_RADIX_03_DFT(r14r,r30r,r46r, cc1, s1p18r,s1p34r,s1p02r)
		SSE2_RADIX_03_DFT(r15r,r31r,r47r, cc1, s1p33r,s1p01r,s1p17r)
	  #endif

	#else	/* !USE_SSE2 */

	/*...gather the needed data (48 64-bit complex) and do 3 radix-16 transforms,	*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_16_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p08+p07],a[jp+p08+p07],a[jt+p08+p06],a[jp+p08+p06],a[jt+p08+p05],a[jp+p08+p05],a[jt+p08+p04],a[jp+p08+p04],a[jt+p08+p03],a[jp+p08+p03],a[jt+p08+p02],a[jp+p08+p02],a[jt+p08+p01],a[jp+p08+p01],a[jt+p08    ],a[jp+p08    ], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im, c16,s16);	tptr += 0x10;
		jt = j1+p32; jp = j2+p32;	RADIX_16_DIT(a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p08+p01],a[jp+p08+p01],a[jt+p08    ],a[jp+p08    ],a[jt+p08+p02],a[jp+p08+p02],a[jt+p08+p03],a[jp+p08+p03],a[jt+p08+p06],a[jp+p08+p06],a[jt+p08+p07],a[jp+p08+p07],a[jt+p08+p04],a[jp+p08+p04],a[jt+p08+p05],a[jp+p08+p05], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im, c16,s16);	tptr += 0x10;
		jt = j1+p16; jp = j2+p16;	RADIX_16_DIT(a[jt+p08+p02],a[jp+p08+p02],a[jt+p08+p03],a[jp+p08+p03],a[jt+p08    ],a[jp+p08    ],a[jt+p08+p01],a[jp+p08+p01],a[jt+p08+p04],a[jp+p08+p04],a[jt+p08+p05],a[jp+p08+p05],a[jt+p08+p07],a[jp+p08+p07],a[jt+p08+p06],a[jp+p08+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im, c16,s16);

	/*...and now do 16 radix-3 transforms.	*/
		tptr = t;
		jt = j1        ; jp = j2        ;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32]);	tptr++;
		jt = j1+p08+p07; jp = j2+p08+p07;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32]);	tptr++;
		jt = j1+p08+p06; jp = j2+p08+p06;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ]);	tptr++;
		jt = j1+p08+p05; jp = j2+p08+p05;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);	tptr++;
		jt = j1+p08+p04; jp = j2+p08+p04;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32]);	tptr++;
		jt = j1+p08+p03; jp = j2+p08+p03;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ]);	tptr++;
		jt = j1+p08+p02; jp = j2+p08+p02;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);	tptr++;
		jt = j1+p08+p01; jp = j2+p08+p01;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32]);	tptr++;
		jt = j1+p08    ; jp = j2+p08    ;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ]);	tptr++;
		jt = j1    +p07; jp = j2    +p07;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);	tptr++;
		jt = j1    +p06; jp = j2    +p06;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32]);	tptr++;
		jt = j1    +p05; jp = j2    +p05;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ]);	tptr++;
		jt = j1    +p04; jp = j2    +p04;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);	tptr++;
		jt = j1    +p03; jp = j2    +p03;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32]);	tptr++;
		jt = j1    +p02; jp = j2    +p02;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ]);	tptr++;
		jt = j1    +p01; jp = j2    +p01;	RADIX_03_DFT(s,c3m1, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im, t00,t01,t02,t03,t04,t05, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16]);

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
		tm1 = s1p00r; tmp = cy; itmp = bjmodn;
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

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00r; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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

	/*...The radix-48 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #if OS_BITS == 64
		SSE2_RADIX_03_DFT_X2(cc1, s1p00r,s1p32r,s1p16r, r00r,r01r,r02r, s1p45r,s1p29r,s1p13r, r03r,r04r,r05r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p42r,s1p26r,s1p10r, r06r,r07r,r08r, s1p39r,s1p23r,s1p07r, r09r,r10r,r11r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p36r,s1p20r,s1p04r, r12r,r13r,r14r, s1p33r,s1p17r,s1p01r, r15r,r16r,r17r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p30r,s1p14r,s1p46r, r18r,r19r,r20r, s1p27r,s1p11r,s1p43r, r21r,r22r,r23r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p24r,s1p08r,s1p40r, r24r,r25r,r26r, s1p21r,s1p05r,s1p37r, r27r,r28r,r29r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p18r,s1p02r,s1p34r, r30r,r31r,r32r, s1p15r,s1p47r,s1p31r, r33r,r34r,r35r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p12r,s1p44r,s1p28r, r36r,r37r,r38r, s1p09r,s1p41r,s1p25r, r39r,r40r,r41r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p06r,s1p38r,s1p22r, r42r,r43r,r44r, s1p03r,s1p35r,s1p19r, r45r,r46r,r47r)
	  #else
		SSE2_RADIX_03_DFT(s1p00r,s1p32r,s1p16r, cc1, r00r,r01r,r02r)
		SSE2_RADIX_03_DFT(s1p45r,s1p29r,s1p13r, cc1, r03r,r04r,r05r)
		SSE2_RADIX_03_DFT(s1p42r,s1p26r,s1p10r, cc1, r06r,r07r,r08r)
		SSE2_RADIX_03_DFT(s1p39r,s1p23r,s1p07r, cc1, r09r,r10r,r11r)
		SSE2_RADIX_03_DFT(s1p36r,s1p20r,s1p04r, cc1, r12r,r13r,r14r)
		SSE2_RADIX_03_DFT(s1p33r,s1p17r,s1p01r, cc1, r15r,r16r,r17r)
		SSE2_RADIX_03_DFT(s1p30r,s1p14r,s1p46r, cc1, r18r,r19r,r20r)
		SSE2_RADIX_03_DFT(s1p27r,s1p11r,s1p43r, cc1, r21r,r22r,r23r)
		SSE2_RADIX_03_DFT(s1p24r,s1p08r,s1p40r, cc1, r24r,r25r,r26r)
		SSE2_RADIX_03_DFT(s1p21r,s1p05r,s1p37r, cc1, r27r,r28r,r29r)
		SSE2_RADIX_03_DFT(s1p18r,s1p02r,s1p34r, cc1, r30r,r31r,r32r)
		SSE2_RADIX_03_DFT(s1p15r,s1p47r,s1p31r, cc1, r33r,r34r,r35r)
		SSE2_RADIX_03_DFT(s1p12r,s1p44r,s1p28r, cc1, r36r,r37r,r38r)
		SSE2_RADIX_03_DFT(s1p09r,s1p41r,s1p25r, cc1, r39r,r40r,r41r)
		SSE2_RADIX_03_DFT(s1p06r,s1p38r,s1p22r, cc1, r42r,r43r,r44r)
		SSE2_RADIX_03_DFT(s1p03r,s1p35r,s1p19r, cc1, r45r,r46r,r47r)
	  #endif

	// istride of [3 vector-complex]*[1,2,3,4] = [1,2,3,4]*0x60 bytes for sse2, [1,2,3,4]*0xc0 bytes for avx
	// Offsets: 	00,01,02,03,05,04,07,06,10,11,09,08,15,14,12,13
		add0 = &a[j1    ];	addb = add0 + p08;
		add1 = add0 + p01;	add8 = addb + p02;
		add2 = add0 + p02;	add9 = addb + p03;
		add3 = add0 + p03;	adda = addb + p01;
		add4 = add0 + p05;	addc = addb + p07;
		add5 = add0 + p04;	addd = addb + p06;
		add6 = add0 + p07;	adde = addb + p04;
		add7 = add0 + p06;	addf = addb + p05;
	  #ifdef USE_AVX
		SSE2_RADIX16_DIF_0TWIDDLE(r00r,0x0c0,0x180,0x240,0x300, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)
	  #else
		SSE2_RADIX16_DIF_0TWIDDLE(r00r,0x060,0x0c0,0x120,0x180, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)
	  #endif

	// Offsets: 32+	05,04,07,06,02,03,01,00,15,14,12,13,09,08,11,10
		add7 = &a[j1+p32];	addd = add7 + p08;
		add0 = add7 + p05;	add8 = addd + p07;
		add1 = add7 + p04;	add9 = addd + p06;
		add2 = add7 + p07;	adda = addd + p04;
		add3 = add7 + p06;	addb = addd + p05;
		add4 = add7 + p02;	addc = addd + p01;
		add5 = add7 + p03;	adde = addd + p03;
		add6 = add7 + p01;	addf = addd + p02;
	  #ifdef USE_AVX
		SSE2_RADIX16_DIF_0TWIDDLE(r01r,0x0c0,0x180,0x240,0x300, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)
	  #else
		SSE2_RADIX16_DIF_0TWIDDLE(r01r,0x060,0x0c0,0x120,0x180, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)
	  #endif

	// Offsets: 16+ 10,11,09,08,15,14,12,13,05,04,07,06,02,03,01,00:
		addf = &a[j1+p16];	add3 = addf + p08;
		add0 = add3 + p02;	add8 = addf + p05;
		add1 = add3 + p03;	add9 = addf + p04;
		add2 = add3 + p01;	adda = addf + p07;
		add4 = add3 + p07;	addb = addf + p06;
		add5 = add3 + p06;	addc = addf + p02;
		add6 = add3 + p04;	addd = addf + p03;
		add7 = add3 + p05;	adde = addf + p01;
	  #ifdef USE_AVX
		SSE2_RADIX16_DIF_0TWIDDLE(r02r,0x0c0,0x180,0x240,0x300, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)
	  #else
		SSE2_RADIX16_DIF_0TWIDDLE(r02r,0x060,0x0c0,0x120,0x180, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)
	  #endif

	#else	/* !USE_SSE2 */

	/*...gather the needed data (48 64-bit complex) and do 16 radix-3 transforms...*/
		tptr = t;
		jt = j1        ; jp = j2        ;	RADIX_03_DFT(s,c3m1,a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1+p08+p05; jp = j2+p08+p05;	RADIX_03_DFT(s,c3m1,a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1+p08+p02; jp = j2+p08+p02;	RADIX_03_DFT(s,c3m1,a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1    +p07; jp = j2    +p07;	RADIX_03_DFT(s,c3m1,a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1    +p04; jp = j2    +p04;	RADIX_03_DFT(s,c3m1,a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1    +p01; jp = j2    +p01;	RADIX_03_DFT(s,c3m1,a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1+p08+p06; jp = j2+p08+p06;	RADIX_03_DFT(s,c3m1,a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1+p08+p03; jp = j2+p08+p03;	RADIX_03_DFT(s,c3m1,a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1    +p08; jp = j2    +p08;	RADIX_03_DFT(s,c3m1,a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1    +p05; jp = j2    +p05;	RADIX_03_DFT(s,c3m1,a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1    +p02; jp = j2    +p02;	RADIX_03_DFT(s,c3m1,a[jt+p16],a[jp+p16], a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1+p08+p07; jp = j2+p08+p07;	RADIX_03_DFT(s,c3m1,a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1+p08+p04; jp = j2+p08+p04;	RADIX_03_DFT(s,c3m1,a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1+p08+p01; jp = j2+p08+p01;	RADIX_03_DFT(s,c3m1,a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1    +p06; jp = j2    +p06;	RADIX_03_DFT(s,c3m1,a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);	tptr++;
		jt = j1    +p03; jp = j2    +p03;	RADIX_03_DFT(s,c3m1,a[jt    ],a[jp    ], a[jt+p32],a[jp+p32], a[jt+p16],a[jp+p16], t00,t01,t02,t03,t04,t05, tptr->re,tptr->im,(tptr+16)->re,(tptr+16)->im,(tptr+32)->re,(tptr+32)->im);

	/*...and now do 3 radix-16 transforms:	*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_16_DIF(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p08+p02],a[jp+p08+p02],a[jt+p08+p03],a[jp+p08+p03],a[jt+p08+p01],a[jp+p08+p01],a[jt+p08    ],a[jp+p08    ],a[jt+p08+p07],a[jp+p08+p07],a[jt+p08+p06],a[jp+p08+p06],a[jt+p08+p04],a[jp+p08+p04],a[jt+p08+p05],a[jp+p08+p05], c16,s16);	tptr += 0x10;
		jt = j1+p32; jp = j2+p32;	RADIX_16_DIF(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im, a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p08+p07],a[jp+p08+p07],a[jt+p08+p06],a[jp+p08+p06],a[jt+p08+p04],a[jp+p08+p04],a[jt+p08+p05],a[jp+p08+p05],a[jt+p08+p01],a[jp+p08+p01],a[jt+p08    ],a[jp+p08    ],a[jt+p08+p03],a[jp+p08+p03],a[jt+p08+p02],a[jp+p08+p02], c16,s16);	tptr += 0x10;
		jt = j1+p16; jp = j2+p16;	RADIX_16_DIF(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im, a[jt+p08+p02],a[jp+p08+p02],a[jt+p08+p03],a[jp+p08+p03],a[jt+p08+p01],a[jp+p08+p01],a[jt+p08    ],a[jp+p08    ],a[jt+p08+p07],a[jp+p08+p07],a[jt+p08+p06],a[jp+p08+p06],a[jt+p08+p04],a[jp+p08+p04],a[jt+p08+p05],a[jp+p08+p05],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], c16,s16);

	#endif	/* !USE_SSE2 */

	}

	jstart += nwt;
	jhi    += nwt;

	col += RADIX;
	co3 -= RADIX;
}	/* end for(k=1; k <= khi; k++) */
