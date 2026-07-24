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
	for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	/*...The radix-56 DIT pass is here:	*/

	#ifdef USE_SSE2

		/* Outputs of radix-8 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between r00r and r01r: */

	  #if COMPACT_OBJ

		// Indices into above 8-elt table; each DFT-8 needs eight 3-bit indices, thus needs low 3 bytes of a uint32.
		// Ex: 1st DFT-8 has add0-7 p01-multiple offsets 01327654; reverse that to get p_idx[0] = 45672310:
		const uint32 p_id1[7] = {045672310,076450123,032017645,054763201,010325476,067541032,023106754};	// Leading 0s = octal consts
		const uint8  p_od1[7] = {0,12,10,8,6,4,2};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00r; l < 7; l++, tmp+=16) {
			i7 = p_id1[l]; i0 = i7&7; i1 = (i7>>3)&7; i2 = (i7>>6)&7; i3 = (i7>>9)&7; i4 = (i7>>12)&7; i5 = (i7>>15)&7; i6 = (i7>>18)&7; i7 = (i7>>21);
			addr = &a[j1+poff[p_od1[l]]];
			add0 = addr+pp07[i0]; add1 = addr+pp07[i1]; add2 = addr+pp07[i2]; add3 = addr+pp07[i3]; add4 = addr+pp07[i4]; add5 = addr+pp07[i5]; add6 = addr+pp07[i6]; add7 = addr+pp07[i7];
		#ifdef USE_AVX2
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, tmp, isrt2,two)
		#else
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, tmp, isrt2)
		#endif
		}
	/*...and now do 8 radix-7 transforms, with the columns of r*[r,i] output pairs in the above 7x radix-8 set now acting as input rows: */
		// Bytewise array of output-pointer offsets w.r.to the s1p00 base-pointer:
		const uint8 optr_off[RADIX] = {
			 0,16,32,48,64,80,96,
			 98,2,18,34,50,66,82,
			 84,100,4,20,36,52,68,
			 70,86,102,6,22,38,54,
			 56,72,88,104,8,24,40,
			 42,58,74,90,106,10,26,
			 28,44,60,76,92,108,12,
			 14,30,46,62,78,94,110};
		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,	// I-ptrs
		*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6;	// O-ptrs
		tmp = r00r;
		for(l = 0, ntmp = 0; l < 8; l++, ntmp += 7) {
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;			vb0 = s1p00r + optr_off[ntmp  ];
			va1 = tmp + 0x10;	vb1 = s1p00r + optr_off[ntmp+1];
			va2 = tmp + 0x20;	vb2 = s1p00r + optr_off[ntmp+2];
			va3 = tmp + 0x30;	vb3 = s1p00r + optr_off[ntmp+3];
			va4 = tmp + 0x40;	vb4 = s1p00r + optr_off[ntmp+4];
			va5 = tmp + 0x50;	vb5 = s1p00r + optr_off[ntmp+5];
			va6 = tmp + 0x60;	vb6 = s1p00r + optr_off[ntmp+6];
							/*            inputs        */         /*        outputs          */
		  #ifdef USE_AVX2
			SSE2_RADIX_07_DFT(va0,va1,va2,va3,va4,va5,va6, cc0,two, vb0,vb1,vb2,vb3,vb4,vb5,vb6);
		  #else
			SSE2_RADIX_07_DFT(va0,va1,va2,va3,va4,va5,va6, cc0,     vb0,vb1,vb2,vb3,vb4,vb5,vb6);
		  #endif
			tmp += 2;
		}

	  #else

		add0 = &a[j1    ]; add1 = add0+p01; add2 = add0+p03; add3 = add0+p02; add4 = add0+p07; add5 = add0+p06; add6 = add0+p05; add7 = add0+p04;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00r, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00r, isrt2)
	   #endif

		add3 = &a[j1+p30]; add0 = add3+p03; add1 = add3+p02; add2 = add3+p01; add4 = add3+p05; add5 = add3+p04; add6 = add3+p06; add7 = add3+p07;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r08r, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r08r, isrt2)
	   #endif

		add5 = &a[j1+p28]; add0 = add5+p05; add1 = add5+p04; add2 = add5+p06; add3 = add5+p07; add4 = add5+p01; add6 = add5+p02; add7 = add5+p03;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r16r, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r16r, isrt2)
	   #endif

		add1 = &a[j1+p20]; add0 = add1+p01; add2 = add1+p02; add3 = add1+p03; add4 = add1+p06; add5 = add1+p07; add6 = add1+p04; add7 = add1+p05;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r24r, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r24r, isrt2)
	   #endif

		add6 = &a[j1+p18]; add0 = add6+p06; add1 = add6+p07; add2 = add6+p04; add3 = add6+p05; add4 = add6+p02; add5 = add6+p03; add7 = add6+p01;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r32r, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r32r, isrt2)
	   #endif

		add2 = &a[j1+p10]; add0 = add2+p02; add1 = add2+p03; add3 = add2+p01; add4 = add2+p04; add5 = add2+p05; add6 = add2+p07; add7 = add2+p06;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r40r, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r40r, isrt2)
	   #endif

		add4 = &a[j1+p08]; add0 = add4+p04; add1 = add4+p05; add2 = add4+p07; add3 = add4+p06; add5 = add4+p01; add6 = add4+p03; add7 = add4+p02;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r48r, isrt2,two)
	   #else
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r48r, isrt2)
	   #endif

	/*...and now do 8 radix-7 transforms, with the columns of r*[r,i] output pairs in the above 7x radix-8 set now acting as input rows: */
						/*            inputs        */  /* sincos ptr */ /*            outputs                   */
	   #ifdef USE_AVX2
		SSE2_RADIX_07_DFT(r00r,r08r,r16r,r24r,r32r,r40r,r48r, cc0,two, s1p00r,s1p08r,s1p16r,s1p24r,s1p32r,s1p40r,s1p48r);
		SSE2_RADIX_07_DFT(r01r,r09r,r17r,r25r,r33r,r41r,r49r, cc0,two, s1p49r,s1p01r,s1p09r,s1p17r,s1p25r,s1p33r,s1p41r);
		SSE2_RADIX_07_DFT(r02r,r10r,r18r,r26r,r34r,r42r,r50r, cc0,two, s1p42r,s1p50r,s1p02r,s1p10r,s1p18r,s1p26r,s1p34r);
		SSE2_RADIX_07_DFT(r03r,r11r,r19r,r27r,r35r,r43r,r51r, cc0,two, s1p35r,s1p43r,s1p51r,s1p03r,s1p11r,s1p19r,s1p27r);
		SSE2_RADIX_07_DFT(r04r,r12r,r20r,r28r,r36r,r44r,r52r, cc0,two, s1p28r,s1p36r,s1p44r,s1p52r,s1p04r,s1p12r,s1p20r);
		SSE2_RADIX_07_DFT(r05r,r13r,r21r,r29r,r37r,r45r,r53r, cc0,two, s1p21r,s1p29r,s1p37r,s1p45r,s1p53r,s1p05r,s1p13r);
		SSE2_RADIX_07_DFT(r06r,r14r,r22r,r30r,r38r,r46r,r54r, cc0,two, s1p14r,s1p22r,s1p30r,s1p38r,s1p46r,s1p54r,s1p06r);
		SSE2_RADIX_07_DFT(r07r,r15r,r23r,r31r,r39r,r47r,r55r, cc0,two, s1p07r,s1p15r,s1p23r,s1p31r,s1p39r,s1p47r,s1p55r);
	   #else
		SSE2_RADIX_07_DFT(r00r,r08r,r16r,r24r,r32r,r40r,r48r, cc0,     s1p00r,s1p08r,s1p16r,s1p24r,s1p32r,s1p40r,s1p48r);
		SSE2_RADIX_07_DFT(r01r,r09r,r17r,r25r,r33r,r41r,r49r, cc0,     s1p49r,s1p01r,s1p09r,s1p17r,s1p25r,s1p33r,s1p41r);
		SSE2_RADIX_07_DFT(r02r,r10r,r18r,r26r,r34r,r42r,r50r, cc0,     s1p42r,s1p50r,s1p02r,s1p10r,s1p18r,s1p26r,s1p34r);
		SSE2_RADIX_07_DFT(r03r,r11r,r19r,r27r,r35r,r43r,r51r, cc0,     s1p35r,s1p43r,s1p51r,s1p03r,s1p11r,s1p19r,s1p27r);
		SSE2_RADIX_07_DFT(r04r,r12r,r20r,r28r,r36r,r44r,r52r, cc0,     s1p28r,s1p36r,s1p44r,s1p52r,s1p04r,s1p12r,s1p20r);
		SSE2_RADIX_07_DFT(r05r,r13r,r21r,r29r,r37r,r45r,r53r, cc0,     s1p21r,s1p29r,s1p37r,s1p45r,s1p53r,s1p05r,s1p13r);
		SSE2_RADIX_07_DFT(r06r,r14r,r22r,r30r,r38r,r46r,r54r, cc0,     s1p14r,s1p22r,s1p30r,s1p38r,s1p46r,s1p54r,s1p06r);
		SSE2_RADIX_07_DFT(r07r,r15r,r23r,r31r,r39r,r47r,r55r, cc0,     s1p07r,s1p15r,s1p23r,s1p31r,s1p39r,s1p47r,s1p55r);
	   #endif

	  #endif	// COMPACT_OBJ?

	#else	// USE_SSE2 = False:

	/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 7 radix-8 transforms...*/
		tptr = t;								/*                                    inputs                                                                                                                         */  /*                  intermediates                            */  /*                         outputs                   */
		jt = j1    ; jp = j2    ;	RADIX_08_DIT(a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, rt,it);	tptr += 8;
		jt = j1+p30; jp = j2+p30;	RADIX_08_DIT(a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, rt,it);	tptr += 8;
		jt = j1+p28; jp = j2+p28;	RADIX_08_DIT(a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, rt,it);	tptr += 8;
		jt = j1+p20; jp = j2+p20;	RADIX_08_DIT(a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, rt,it);	tptr += 8;
		jt = j1+p18; jp = j2+p18;	RADIX_08_DIT(a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, rt,it);	tptr += 8;
		jt = j1+p10; jp = j2+p10;	RADIX_08_DIT(a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, rt,it);	tptr += 8;
		jt = j1+p08; jp = j2+p08;	RADIX_08_DIT(a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, rt,it);

	/*...and now do 8 radix-7 transforms, with the columns of r*[r,i] output pairs in the above 7x radix-8 set now acting as input rows: */
		tptr = t;
	  #if LO_ADD										/*                                                                           inputs                                                                                                */   /*                         intermediates             */  /*                                                               outputs                                                                 */  /*      sincos + misc temps      */
		jt = j1    ; jp = j2    ;	RADIX_07_DFT     (tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p01; jp = j2+p01;	RADIX_07_DFT     (tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p02; jp = j2+p02;	RADIX_07_DFT     (tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p03; jp = j2+p03;	RADIX_07_DFT     (tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_07_DFT     (tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p05; jp = j2+p05;	RADIX_07_DFT     (tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p06; jp = j2+p06;	RADIX_07_DFT     (tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p07; jp = j2+p07;	RADIX_07_DFT     (tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	  #else
		jt = j1    ; jp = j2    ;	RADIX_07_DFT_NUSS(tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p01; jp = j2+p01;	RADIX_07_DFT_NUSS(tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p02; jp = j2+p02;	RADIX_07_DFT_NUSS(tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p03; jp = j2+p03;	RADIX_07_DFT_NUSS(tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_07_DFT_NUSS(tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p05; jp = j2+p05;	RADIX_07_DFT_NUSS(tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p06; jp = j2+p06;	RADIX_07_DFT_NUSS(tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt    ],a[jp    ], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p07; jp = j2+p07;	RADIX_07_DFT_NUSS(tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30], cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);
	  #endif

	#endif	// SIMD or not?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need 28 separate carries.	*/

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
		/* ptr to local storage for the doubled wtl,wtn terms: */
	  #ifdef USE_AVX512
		tmp = half_arr +  64;	// No lookup-tables used in avx-512; instead use opmasked conditional-doubling;
								// 1st 64 slots hold outputs of wtsinit call. Only half of said slots used in 8-way-init mode.
	  #else
		tmp = half_arr + 128;	// 1st 64 slots are basic-4 LUTs, next 32 are the additional 2 LOACC LUTs, next 32 hold outputs of wtsinit call
	  #endif
		l= j & (nwt-1);						// These rcol wts-terms are for individual-double-broadcast-to-full-vector-width,
		n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];	// hence the mixing of fwd/inv wts, which is normally taboo.
		n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+4) & (nwt-1);					++tmp;
		n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+6) & (nwt-1);					++tmp;
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

	  if(inc_arr[0]) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of the low
						// inc_arr element to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
		i = (!j);
		addr = &prp_mult;
		tmp = s1p00r; tm1 = cy_r; tm2 = cy_r+1; itmp = bjmodn;
	  #ifndef USE_AVX512
		itm2 = bjmodn+4;
	  #endif
		for(l = 0; l < RADIX>>3; l++) {
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l+l];
		  #ifdef USE_AVX512	// In AVX-512 mode, the 4 doubles base[0],baseinv[1],wts_mult[1],inv_mult[0] are in the d0-3 slots of the otherwise-unused sse2_rnd vec_dbl:
			AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
			tmp += 16; tm1 += 1;           itmp += 8;            i = 0;	// CY-ptr only advances 1 in AVX-512 mode, since all 8 dbl-carries fit in a single vec_dbl
		  #else	// USE_AVX:
			AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
			tmp += 16; tm1 += 2; tm2 += 2; itmp += 8; itm2 += 8; i = 0;
		  #endif
		}

	  } else {	// HiACC:

		/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		tm1 = s1p00r; tmp = cy_r; itmp = bjmodn;
		addr = &prp_mult;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[l]);	// poff[] = p0,4,8,...
			AVX_cmplx_carry_norm_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, tm2,p01,p02,p03, addr);
			tm1 += 8; tmp += 1; itmp += 4; i = 0;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	  if(inc_arr[0]) {	// Have no specialized HIACC carry macro in ARM_V8_SIMD, so use 0-or-not-ness of incr
						// in lieu of (USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) is for ARMv8
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		uint32 k0,k1,k2,k3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;

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
			k0 = n-si[l  ];
			k1 = n-si[l+1];
			k2 = si[nwtml  ];
			k3 = si[nwtml-1];
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

		add1 = &wt1[col  ];
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

	/*	i =((uint32)(sw - bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

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

		// AVX-custom 4-way carry macro - each macro call contains 4 of the RADIX stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indices j,j+2,j+4,j+6.
		// For non-power-of-2 FFT lengths we have 2 versions of the AVX carry sequence, tradong off speed (3-5%) vs accuracy:
	#ifdef USE_AVX
		int k3,k4,k5,k6,k7;
		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		tmp = base_negacyclic_root;	tm2 = tmp+1;

	  #ifdef HIACC
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

		tm0 = s1p00r; tmp = base_negacyclic_root; l = 0xe00;
		tm1 = cy_r; // *cycle[] indices increment by +4 (mod ODD_RADIX) between macro calls
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic_idx = 0; jc_idx = 1; kc_idx = 2; lc_idx = 3;
		while(tm0 < s1p55r)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{	// NB: (int)(tmp-cy_r) < RADIX (as used for SSE2 build) no good here, since just 1 vec_dbl increment
			// per 4 Re+Im-carries; but (int)(tmp-cy_r) < (RADIX>>1) would work
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic_idx];	k5 = jcycle[ic_idx];	k6 = kcycle[ic_idx];	k7 = lcycle[ic_idx];
			k2 = icycle[jc_idx];
			k3 = icycle[kc_idx];
			k4 = icycle[lc_idx];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x1c0, 0xe0,0x1c0,0x2a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p01,p02,p03, addr);
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
			k1=(l & NRTM1);		k2=(l >> NRT_BITS);
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
	   #ifdef USE_AVX512
		mc_idx = 4; nc_idx = 5; oc_idx = 6; pc_idx = 7;
	   #endif

	   #ifdef USE_AVX512
		int k8,k9,ka,kb,kc,kd,ke,kf;

		for(l = 0; l < RADIX>>3; l++) {	// RADIX/8 loop passes
			k1 = icycle[ic_idx];
			k2 = icycle[jc_idx];	k9 = jcycle[ic_idx];
			k3 = icycle[kc_idx];	ka = kcycle[ic_idx];
			k4 = icycle[lc_idx];	kb = lcycle[ic_idx];
			k5 = icycle[mc_idx];	kc = mcycle[ic_idx];
			k6 = icycle[nc_idx];	kd = ncycle[ic_idx];
			k7 = icycle[oc_idx];	ke = ocycle[ic_idx];
			k8 = icycle[pc_idx];	kf = pcycle[ic_idx];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
												/* (cy_i_cy_r) --vvvvv  vvvvvvvvvvvvvvvvv--[1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X8_loacc(tm0,tmp,tm1,0x1c0, 0x1c0,0x380,0x540, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, tm2,p01,p02,p03,p04, addr);
			tm0 += 16; tm1++;
			// *** BUG: *** Aug 2021: Needed to reduce the constant addend 8 mod-ODD_RADIX, i.e. 8%7 = 1:
			MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			MOD_ADD32(jc_idx, 1, ODD_RADIX, jc_idx);
			MOD_ADD32(kc_idx, 1, ODD_RADIX, kc_idx);
			MOD_ADD32(lc_idx, 1, ODD_RADIX, lc_idx);
			MOD_ADD32(mc_idx, 1, ODD_RADIX, mc_idx);
			MOD_ADD32(nc_idx, 1, ODD_RADIX, nc_idx);
			MOD_ADD32(oc_idx, 1, ODD_RADIX, oc_idx);
			MOD_ADD32(pc_idx, 1, ODD_RADIX, pc_idx);
		}

	   #else	// AVX / AVX2

		for(l = 0; l < RADIX>>2; l++) {	// RADIX/4 loop passes
			k1 = icycle[ic_idx];
			k2 = icycle[jc_idx];	k5 = jcycle[ic_idx];
			k3 = icycle[kc_idx];	k6 = kcycle[ic_idx];
			k4 = icycle[lc_idx];	k7 = lcycle[ic_idx];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
												/* (cy_i_cy_r) --vvvvv  vvvvvvvvvvvvvvvvv--[1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x1c0, 0xe0,0x1c0,0x2a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p01,p02,p03, addr);
			tm0 += 8; tm1++;
			MOD_ADD32(ic_idx, 4, ODD_RADIX, ic_idx);
			MOD_ADD32(jc_idx, 4, ODD_RADIX, jc_idx);
			MOD_ADD32(kc_idx, 4, ODD_RADIX, kc_idx);
			MOD_ADD32(lc_idx, 4, ODD_RADIX, lc_idx);
		}

	   #endif

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
			k1 = icycle[ic_idx];
			k2 = jcycle[ic_idx];
			int k3 = icycle[jc_idx];
			int k4 = jcycle[jc_idx];
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
		ntmp = 0; addr = cy_r; addi = cy_i; ic_idx = 0;	// ic_idx = idx into icycle mini-array, gets incremented (mod ODD_RADIX) between macro calls
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];
			fermat_carry_norm_errcheckB(a[jt    ],a[jp    ],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p01],a[jp+p01],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p02],a[jp+p02],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p03],a[jp+p03],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
		}
		//if(!jt)printf("Iter = %d, full = %u, a[0]_out = %20.15f\n",iter,full_pass,a[0]);
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
		#ifdef USE_AVX512
			mcycle[ntmp] += wts_idx_inc2;		mcycle[ntmp] += ( (-(mcycle[ntmp] < 0)) & nwt16);
			ncycle[ntmp] += wts_idx_inc2;		ncycle[ntmp] += ( (-(ncycle[ntmp] < 0)) & nwt16);
			ocycle[ntmp] += wts_idx_inc2;		ocycle[ntmp] += ( (-(ocycle[ntmp] < 0)) & nwt16);
			pcycle[ntmp] += wts_idx_inc2;		pcycle[ntmp] += ( (-(pcycle[ntmp] < 0)) & nwt16);
		#endif
		}
	#endif

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-56 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7;
		OFF1 = 0x20;
		OFF2 = 0x40;
		OFF3 = 0x60;
		OFF4 = 0x80;
		OFF5 = 0xa0;
		OFF6 = 0xc0;
		OFF7 = 0xe0;
	  #elif defined(USE_AVX512)
		#define OFF1	0x20*4
		#define OFF2	0x40*4
		#define OFF3	0x60*4
		#define OFF4	0x80*4
		#define OFF5	0xa0*4
		#define OFF6	0xc0*4
		#define OFF7	0xe0*4
	   #elif defined(USE_AVX)
		#define OFF1	0x20*2
		#define OFF2	0x40*2
		#define OFF3	0x60*2
		#define OFF4	0x80*2
		#define OFF5	0xa0*2
		#define OFF6	0xc0*2
		#define OFF7	0xe0*2
	   #else
		#define OFF1	0x20
		#define OFF2	0x40
		#define OFF3	0x60
		#define OFF4	0x80
		#define OFF5	0xa0
		#define OFF6	0xc0
		#define OFF7	0xe0
	   #endif

	  #if COMPACT_OBJ

	/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 8 radix-7 transforms...*/
		// Bytewise array of input-pointer offsets w.r.to the s1p00 base-pointer:
		const uint8 iptr_off[RADIX] = {
			 0,96,80,64,48,32,16,
			98,82,66,50,34,18, 2,
			84,68,52,36,20,4,100,
			70,54,38,22,6,102,86,
			56,40,24,8,104,88,72,
			42,26,10,106,90,74,58,
			28,12,108,92,76,60,44,
			14,110,94,78,62,46,30};
		tmp = r00r;
		for(l = 0, ntmp = 0; l < 8; l++, ntmp += 7) {
			// Output-ptrs [va* here] are regular-stride offsets of r00:
			va0 = tmp;			vb0 = s1p00r + iptr_off[ntmp  ];
			va1 = tmp + 0x10;	vb1 = s1p00r + iptr_off[ntmp+1];
			va2 = tmp + 0x20;	vb2 = s1p00r + iptr_off[ntmp+2];
			va3 = tmp + 0x30;	vb3 = s1p00r + iptr_off[ntmp+3];
			va4 = tmp + 0x40;	vb4 = s1p00r + iptr_off[ntmp+4];
			va5 = tmp + 0x50;	vb5 = s1p00r + iptr_off[ntmp+5];
			va6 = tmp + 0x60;	vb6 = s1p00r + iptr_off[ntmp+6];
							/*            inputs        */         /*        outputs          */
		  #ifdef USE_AVX2
			SSE2_RADIX_07_DFT(vb0,vb1,vb2,vb3,vb4,vb5,vb6, cc0,two, va0,va1,va2,va3,va4,va5,va6);
		  #else
			SSE2_RADIX_07_DFT(vb0,vb1,vb2,vb3,vb4,vb5,vb6, cc0,     va0,va1,va2,va3,va4,va5,va6);
		  #endif
			tmp += 2;
		}
	/*...and now do 7 radix-8 transforms: */
		const uint32 p_id2[7] = {076543210,054671023,010236745,067452301,023014576,045760132,001327654};	// Leading 0s = octal consts
		for(l = 0, tmp = r00r; l < 7; l++, tmp+=16) {
			i7 = p_id2[l]; i0 = i7&7; i1 = (i7>>3)&7; i2 = (i7>>6)&7; i3 = (i7>>9)&7; i4 = (i7>>12)&7; i5 = (i7>>15)&7; i6 = (i7>>18)&7; i7 = (i7>>21);
			addr = &a[j1+poff[p_od1[l]]];
			add0 = addr+pp07[i0]; add1 = addr+pp07[i1]; add2 = addr+pp07[i2]; add3 = addr+pp07[i3]; add4 = addr+pp07[i4]; add5 = addr+pp07[i5]; add6 = addr+pp07[i6]; add7 = addr+pp07[i7];
		#ifdef USE_AVX2
			SSE2_RADIX8_DIF_0TWIDDLE(tmp,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
		#else
			SSE2_RADIX8_DIF_0TWIDDLE(tmp,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
		#endif
		}

	  #else

	/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 8 radix-7 transforms...*/
					 /*                   inputs                    */ /* sincos */ /*         outputs           */
	   #ifdef USE_AVX2
		SSE2_RADIX_07_DFT(s1p00r,s1p48r,s1p40r,s1p32r,s1p24r,s1p16r,s1p08r, cc0,two, r00r,r08r,r16r,r24r,r32r,r40r,r48r);
		SSE2_RADIX_07_DFT(s1p49r,s1p41r,s1p33r,s1p25r,s1p17r,s1p09r,s1p01r, cc0,two, r01r,r09r,r17r,r25r,r33r,r41r,r49r);
		SSE2_RADIX_07_DFT(s1p42r,s1p34r,s1p26r,s1p18r,s1p10r,s1p02r,s1p50r, cc0,two, r02r,r10r,r18r,r26r,r34r,r42r,r50r);
		SSE2_RADIX_07_DFT(s1p35r,s1p27r,s1p19r,s1p11r,s1p03r,s1p51r,s1p43r, cc0,two, r03r,r11r,r19r,r27r,r35r,r43r,r51r);
		SSE2_RADIX_07_DFT(s1p28r,s1p20r,s1p12r,s1p04r,s1p52r,s1p44r,s1p36r, cc0,two, r04r,r12r,r20r,r28r,r36r,r44r,r52r);
		SSE2_RADIX_07_DFT(s1p21r,s1p13r,s1p05r,s1p53r,s1p45r,s1p37r,s1p29r, cc0,two, r05r,r13r,r21r,r29r,r37r,r45r,r53r);
		SSE2_RADIX_07_DFT(s1p14r,s1p06r,s1p54r,s1p46r,s1p38r,s1p30r,s1p22r, cc0,two, r06r,r14r,r22r,r30r,r38r,r46r,r54r);
		SSE2_RADIX_07_DFT(s1p07r,s1p55r,s1p47r,s1p39r,s1p31r,s1p23r,s1p15r, cc0,two, r07r,r15r,r23r,r31r,r39r,r47r,r55r);
	   #else
		SSE2_RADIX_07_DFT(s1p00r,s1p48r,s1p40r,s1p32r,s1p24r,s1p16r,s1p08r, cc0,     r00r,r08r,r16r,r24r,r32r,r40r,r48r);
		SSE2_RADIX_07_DFT(s1p49r,s1p41r,s1p33r,s1p25r,s1p17r,s1p09r,s1p01r, cc0,     r01r,r09r,r17r,r25r,r33r,r41r,r49r);
		SSE2_RADIX_07_DFT(s1p42r,s1p34r,s1p26r,s1p18r,s1p10r,s1p02r,s1p50r, cc0,     r02r,r10r,r18r,r26r,r34r,r42r,r50r);
		SSE2_RADIX_07_DFT(s1p35r,s1p27r,s1p19r,s1p11r,s1p03r,s1p51r,s1p43r, cc0,     r03r,r11r,r19r,r27r,r35r,r43r,r51r);
		SSE2_RADIX_07_DFT(s1p28r,s1p20r,s1p12r,s1p04r,s1p52r,s1p44r,s1p36r, cc0,     r04r,r12r,r20r,r28r,r36r,r44r,r52r);
		SSE2_RADIX_07_DFT(s1p21r,s1p13r,s1p05r,s1p53r,s1p45r,s1p37r,s1p29r, cc0,     r05r,r13r,r21r,r29r,r37r,r45r,r53r);
		SSE2_RADIX_07_DFT(s1p14r,s1p06r,s1p54r,s1p46r,s1p38r,s1p30r,s1p22r, cc0,     r06r,r14r,r22r,r30r,r38r,r46r,r54r);
		SSE2_RADIX_07_DFT(s1p07r,s1p55r,s1p47r,s1p39r,s1p31r,s1p23r,s1p15r, cc0,     r07r,r15r,r23r,r31r,r39r,r47r,r55r);
	   #endif
	/*...and now do 7 radix-8 transforms: */
					 /*                                   inputs                                  */ /*                       intermediates                       */ /*                 outputs                   */
		add0 = &a[j1    ]; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r00r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	   #else
		SSE2_RADIX8_DIF_0TWIDDLE(r00r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	   #endif

		add2 = &a[j1+p30]; add0 = add2+p03; add1 = add2+p02; add3 = add2+p01; add4 = add2+p07; add5 = add2+p06; add6 = add2+p04; add7 = add2+p05;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r08r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	   #else
		SSE2_RADIX8_DIF_0TWIDDLE(r08r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	   #endif

		add6 = &a[j1+p28]; add0 = add6+p05; add1 = add6+p04; add2 = add6+p07; add3 = add6+p06; add4 = add6+p03; add5 = add6+p02; add7 = add6+p01;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r16r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	   #else
		SSE2_RADIX8_DIF_0TWIDDLE(r16r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	   #endif

		add1 = &a[j1+p20]; add0 = add1+p01; add2 = add1+p03; add3 = add1+p02; add4 = add1+p05; add5 = add1+p04; add6 = add1+p07; add7 = add1+p06;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r24r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	   #else
		SSE2_RADIX8_DIF_0TWIDDLE(r24r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	   #endif

		add5 = &a[j1+p18]; add0 = add5+p06; add1 = add5+p07; add2 = add5+p05; add3 = add5+p04; add4 = add5+p01; add6 = add5+p03; add7 = add5+p02;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r32r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	   #else
		SSE2_RADIX8_DIF_0TWIDDLE(r32r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	   #endif

		add3 = &a[j1+p10]; add0 = add3+p02; add1 = add3+p03; add2 = add3+p01; add4 = add3+p06; add5 = add3+p07; add6 = add3+p05; add7 = add3+p04;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r40r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	   #else
		SSE2_RADIX8_DIF_0TWIDDLE(r40r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	   #endif

		add7 = &a[j1+p08]; add0 = add7+p04; add1 = add7+p05; add2 = add7+p06; add3 = add7+p07; add4 = add7+p02; add5 = add7+p03; add6 = add7+p01;
	   #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(r48r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
	   #else
		SSE2_RADIX8_DIF_0TWIDDLE(r48r,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
	   #endif

	  #endif	// COMPACT_OBJ?

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
		#undef OFF5
		#undef OFF6
		#undef OFF7
	  #endif

	#else	// USE_SSE2 = False:

	/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 8 radix-7 transforms...*/
		tptr = t;
														/*                                                                  inputs                                                                   */  /*                 intermediates                     */  /*                                        outputs                                                                                                                                    */  /*   sincos consts   */  /* tmps */
	  #if LO_ADD
		jt = j1    ; jp = j2    ;	RADIX_07_DFT     (a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p01; jp = j2+p01;	RADIX_07_DFT     (a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p02; jp = j2+p02;	RADIX_07_DFT     (a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p03; jp = j2+p03;	RADIX_07_DFT     (a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_07_DFT     (a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p05; jp = j2+p05;	RADIX_07_DFT     (a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p06; jp = j2+p06;	RADIX_07_DFT     (a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p07; jp = j2+p07;	RADIX_07_DFT     (a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);
	  #else
		jt = j1    ; jp = j2    ;	RADIX_07_DFT_NUSS(a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p01; jp = j2+p01;	RADIX_07_DFT_NUSS(a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p02; jp = j2+p02;	RADIX_07_DFT_NUSS(a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p03; jp = j2+p03;	RADIX_07_DFT_NUSS(a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_07_DFT_NUSS(a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p05; jp = j2+p05;	RADIX_07_DFT_NUSS(a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p06; jp = j2+p06;	RADIX_07_DFT_NUSS(a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);	tptr++;
		jt = j1+p07; jp = j2+p07;	RADIX_07_DFT_NUSS(a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+8)->re,(tptr+8)->im,(tptr+16)->re,(tptr+16)->im,(tptr+24)->re,(tptr+24)->im,(tptr+32)->re,(tptr+32)->im,(tptr+40)->re,(tptr+40)->im,(tptr+48)->re,(tptr+48)->im, cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3, rt,it);
	  #endif

	/*...and now do 7 radix-8 transforms: */
		tptr = t;								 /*                                                          inputs                                                                                                                                   */  /*                       intermediates                       */ /*                                                              outputs                                                                                      */  /* temps */
		jt = j1    ; jp = j2    ;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07], rt,it);	tptr += 8;
		jt = j1+p30; jp = j2+p30;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05], rt,it);	tptr += 8;
		jt = j1+p28; jp = j2+p28;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], rt,it);	tptr += 8;
		jt = j1+p20; jp = j2+p20;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06], rt,it);	tptr += 8;
		jt = j1+p18; jp = j2+p18;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], rt,it);	tptr += 8;
		jt = j1+p10; jp = j2+p10;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04], rt,it);	tptr += 8;
		jt = j1+p08; jp = j2+p08;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], rt,it);

	#endif	/* !USE_SSE2 */
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(int k=1; k <= khi; k++) */
