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
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-56 DIT pass is here:	*/

	#ifdef USE_SSE2

		/* Outputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between r00r and r01r: */

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

	#endif	// USE_SSE2?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need 28 separate carries.	*/

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
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
		tm1 = s1p00r; tmp = cy_r; itmp = bjmodn;
		// Each AVX carry macro call also processes 4 prefetches of main-array data
		tm2 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_norm_errcheck0_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, tm2,p01,p02,p03);
		tm1 += 8; tmp += 1; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			AVX_cmplx_carry_norm_errcheck1_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, tm2,p01,p02,p03);
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

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00r; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
		add0 = a + j1 + pfetch_dist;
		SSE2_cmplx_carry_norm_errcheck0_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01);
		tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p02,p03);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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
		l = 0; addr = cy_r; itmp = bjmodn;
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

	}		/************************************************************************/
	else	/*                MODULUS_TYPE_FERMAT:                                 */
	{		/************************************************************************/

		// AVX-custom 4-way carry macro - each macro call contains 4 of the RADIX stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indices j,j+2,j+4,j+6.
		// For non-power-of-2 FFT lengths we have 2 versions of the AVX carry sequence, tradong off speed (3-5%) vs accuracy:
	#ifdef USE_AVX
		int k3,k4,k5,k6,k7;
		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

		tmp = base_negacyclic_root;	tm2 = tmp+1;

	  #if HIACC
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

	  #else	// HIACC = false:

		// Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
		// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
		l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

		// The above need some inits to prepare for the AVX version of the Fermat-mod carry macro:
		SSE2_fermat_carry_init_loacc(base_negacyclic_root);

	  #endif

		// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indices j,j+2,j+4,j+6:
	  #if HIACC

		/* The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
		which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls: */

		tm0 = s1p00r; tmp = base_negacyclic_root; l = 0xe00;
		tm1 = cy_r; // *cycle[] indices increment by +4 (mod ODD_RADIX) between macro calls
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic = 0; jc = 1; kc = 2; lc = 3;
		while(tm0 < s1p55r)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[jc];
			k3 = icycle[kc];
			k4 = icycle[lc];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)];	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x1c0, 0xe0,0x1c0,0x2a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p01,p02,p03);
			tm0 += 8; tm1++; tmp += 8; l -= 0xc0;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #else	/* HIACC = false: */

		tm0 = s1p00r; tmp = base_negacyclic_root;	// tmp *not* incremented between macro calls in loacc version
		tm1 = cy_r;
		ic = 0; jc = 1; kc = 2; lc = 3;
		for(l = 0; l < RADIX>>2; l++) {	// RADIX/4 loop passes
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[jc];
			k3 = icycle[kc];
			k4 = icycle[lc];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)];	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x1c0, 0xe0,0x1c0,0x2a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p01,p02,p03);
			tm0 += 8; tm1++;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #endif	/* HIACC? */

	#elif defined(USE_SSE2)

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

	  #if (OS_BITS == 64)

		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],jcycle[ic],icycle[ic+1],jcycle[ic+1] of the non-looped version with icycle[ic],jcycle[ic],icycle[jc],jcycle[jc]:
		ic = 0; jc = 1;
		tm1 = s1p00r; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while(tm1 < s1p55r) {
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];
			k2 = jcycle[ic];
			int k3 = icycle[jc];
			int k4 = jcycle[jc];
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			tm2 = a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)];	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
			tm2 += (-((int)(tm1-cy_r)&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_fermat_carry_norm_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2,k3,k4, tm2,p01);
			tm1 += 4; tmp += 2;
			MOD_ADD32(ic, 2, ODD_RADIX, ic);
			MOD_ADD32(jc, 2, ODD_RADIX, jc);
		}

	  #else // Mar 2014: Worked around the out-of-regs compiler issues with the _X2 version of this macro (the
			// code in carry_gcc64.h has details), but keep non-X2 version in case hit out-of-regs again at some point

		ic = 0;	// ic = idx into [i|j]cycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls
		tm1 = s1p00r; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		l = ODD_RADIX << 4;	// 32-bit version needs preshifted << 4 input value
		while(tm1 <= s1p55r) {
			//Sep 2014: Even with reduced-register version of the 32-bit Fermat-mod carry macro,
			// GCC runs out of registers on this one, without some playing-around-with-alternate code-sequences ...
			// Pulling the array-refs out of the carry-macro call like so solves the problem:
			k1 = icycle[ic];
			k2 = jcycle[ic];
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			tm2 = a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)];	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
			tm2 += p01*((int)(tm1-cy_r)&0x3);	// Added offset cycles among p0,1,2,3
			SSE2_fermat_carry_norm_errcheck(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2, tm2);
			tm1 += 2; tmp++;
			MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}

	  #endif

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		ntmp = 0; addr = cy_r; addi = cy_i; ic = 0;	// ic = idx into icycle mini-array, gets incremented (mod ODD_RADIX) between macro calls
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];	// poff[] = p04,p08,...,p56
			fermat_carry_norm_errcheckB(a[jt    ],a[jp    ],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p01],a[jp+p01],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p02],a[jp+p02],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p03],a[jp+p03],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}

		icycle[ 0] += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
		icycle[ 1] += wts_idx_incr;
		icycle[ 2] += wts_idx_incr;
		icycle[ 3] += wts_idx_incr;
		icycle[ 4] += wts_idx_incr;
		icycle[ 5] += wts_idx_incr;
		icycle[ 6] += wts_idx_incr;
		icycle[ 0] += ( (-(int)((uint32)icycle[ 0] >> 31)) & nwt);
		icycle[ 1] += ( (-(int)((uint32)icycle[ 1] >> 31)) & nwt);
		icycle[ 2] += ( (-(int)((uint32)icycle[ 2] >> 31)) & nwt);
		icycle[ 3] += ( (-(int)((uint32)icycle[ 3] >> 31)) & nwt);
		icycle[ 4] += ( (-(int)((uint32)icycle[ 4] >> 31)) & nwt);
		icycle[ 5] += ( (-(int)((uint32)icycle[ 5] >> 31)) & nwt);
		icycle[ 6] += ( (-(int)((uint32)icycle[ 6] >> 31)) & nwt);

	#endif	/* #ifdef USE_SSE2 */

	// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
	#ifdef USE_SSE2

		icycle[ 0] += wts_idx_inc2;		icycle[ 0] += ( (-(icycle[ 0] < 0)) & nwt16);
		icycle[ 1] += wts_idx_inc2;		icycle[ 1] += ( (-(icycle[ 1] < 0)) & nwt16);
		icycle[ 2] += wts_idx_inc2;		icycle[ 2] += ( (-(icycle[ 2] < 0)) & nwt16);
		icycle[ 3] += wts_idx_inc2;		icycle[ 3] += ( (-(icycle[ 3] < 0)) & nwt16);
		icycle[ 4] += wts_idx_inc2;		icycle[ 4] += ( (-(icycle[ 4] < 0)) & nwt16);
		icycle[ 5] += wts_idx_inc2;		icycle[ 5] += ( (-(icycle[ 5] < 0)) & nwt16);
		icycle[ 6] += wts_idx_inc2;		icycle[ 6] += ( (-(icycle[ 6] < 0)) & nwt16);

		jcycle[ 0] += wts_idx_inc2;		jcycle[ 0] += ( (-(jcycle[ 0] < 0)) & nwt16);
		jcycle[ 1] += wts_idx_inc2;		jcycle[ 1] += ( (-(jcycle[ 1] < 0)) & nwt16);
		jcycle[ 2] += wts_idx_inc2;		jcycle[ 2] += ( (-(jcycle[ 2] < 0)) & nwt16);
		jcycle[ 3] += wts_idx_inc2;		jcycle[ 3] += ( (-(jcycle[ 3] < 0)) & nwt16);
		jcycle[ 4] += wts_idx_inc2;		jcycle[ 4] += ( (-(jcycle[ 4] < 0)) & nwt16);
		jcycle[ 5] += wts_idx_inc2;		jcycle[ 5] += ( (-(jcycle[ 5] < 0)) & nwt16);
		jcycle[ 6] += wts_idx_inc2;		jcycle[ 6] += ( (-(jcycle[ 6] < 0)) & nwt16);

	  #ifdef USE_AVX
		kcycle[ 0] += wts_idx_inc2;		kcycle[ 0] += ( (-(kcycle[ 0] < 0)) & nwt16);
		kcycle[ 1] += wts_idx_inc2;		kcycle[ 1] += ( (-(kcycle[ 1] < 0)) & nwt16);
		kcycle[ 2] += wts_idx_inc2;		kcycle[ 2] += ( (-(kcycle[ 2] < 0)) & nwt16);
		kcycle[ 3] += wts_idx_inc2;		kcycle[ 3] += ( (-(kcycle[ 3] < 0)) & nwt16);
		kcycle[ 4] += wts_idx_inc2;		kcycle[ 4] += ( (-(kcycle[ 4] < 0)) & nwt16);
		kcycle[ 5] += wts_idx_inc2;		kcycle[ 5] += ( (-(kcycle[ 5] < 0)) & nwt16);
		kcycle[ 6] += wts_idx_inc2;		kcycle[ 6] += ( (-(kcycle[ 6] < 0)) & nwt16);

		lcycle[ 0] += wts_idx_inc2;		lcycle[ 0] += ( (-(lcycle[ 0] < 0)) & nwt16);
		lcycle[ 1] += wts_idx_inc2;		lcycle[ 1] += ( (-(lcycle[ 1] < 0)) & nwt16);
		lcycle[ 2] += wts_idx_inc2;		lcycle[ 2] += ( (-(lcycle[ 2] < 0)) & nwt16);
		lcycle[ 3] += wts_idx_inc2;		lcycle[ 3] += ( (-(lcycle[ 3] < 0)) & nwt16);
		lcycle[ 4] += wts_idx_inc2;		lcycle[ 4] += ( (-(lcycle[ 4] < 0)) & nwt16);
		lcycle[ 5] += wts_idx_inc2;		lcycle[ 5] += ( (-(lcycle[ 5] < 0)) & nwt16);
		lcycle[ 6] += wts_idx_inc2;		lcycle[ 6] += ( (-(lcycle[ 6] < 0)) & nwt16);
	  #endif
	#endif

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-56 DIF pass is here:	*/

	#ifdef USE_SSE2

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

		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
		#undef OFF5
		#undef OFF6
		#undef OFF7

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
}	/* end for(k=1; k <= khi; k++) */
