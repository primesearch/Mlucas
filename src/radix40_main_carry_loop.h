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
	/*
	!...gather the needed data (40 64-bit complex, i.e. 80 64-bit reals) and do a radix-40 DIT transform...
	*/
	#ifdef USE_SSE2

	  #if COMPACT_OBJ

		// Indices into above 8-elt table; each DFT-8 needs eight 3-bit indices, thus needs low 3 bytes of a uint32.
		// Ex: 1st DFT-8 has add0-7 p01-multiple offsets 01327654; reverse that to get p_idx[0] = 45672310:
		const uint32 p_id1[5] = {045672310,076450123,054763201,010325476,023106754};	// Leading 0s = octal consts
		const uint8  p_od1[5] = {0,4,8,2,6};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00; l < 5; l++, tmp+=16) {
			i7 = p_id1[l]; i0 = i7&7; i1 = (i7>>3)&7; i2 = (i7>>6)&7; i3 = (i7>>9)&7; i4 = (i7>>12)&7; i5 = (i7>>15)&7; i6 = (i7>>18)&7; i7 = (i7>>21);
			addr = &a[j1+poff[p_od1[l]]];
			add0 = addr+pp07[i0]; add1 = addr+pp07[i1]; add2 = addr+pp07[i2]; add3 = addr+pp07[i3]; add4 = addr+pp07[i4]; add5 = addr+pp07[i5]; add6 = addr+pp07[i6]; add7 = addr+pp07[i7];
		#ifdef USE_AVX2
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, tmp, isrt2,two)
		#else
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, tmp, isrt2)
		#endif
		}
	/*...and now do 8 radix-5 transforms: */
		// Bytewise array of output-pointer offsets w.r.to the s1p00 base-pointer.
		// Radix-5 DFT outputs are (cyclic) with vec_dbl-pointer += 32 (mod 80) between successive outputs:
		const uint8 optr_off[RADIX] = {
			 0,32,64,16,48,
			 50,2,34,66,18,
			 20,52,4,36,68,
			 70,22,54,6,38,
			 40,72,24,56,8,
			 10,42,74,26,58,
			 60,12,44,76,28,
			 30,62,14,46,78};
		vec_dbl
		*va0,*va1,*va2,*va3,*va4,	// I-ptrs
		*vb0,*vb1,*vb2,*vb3,*vb4;	// O-ptrs
		vec_dbl
		*vc0,*vc1,*vc2,*vc3,*vc4,	// I-ptrs
		*vd0,*vd1,*vd2,*vd3,*vd4;	// O-ptrs
		for(l = 0, tmp = r00, ntmp = 0; l < 4; l++, ntmp += 10)
		{
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;			vb0 = s1p00r + optr_off[ntmp  ];
			va1 = tmp + 0x10;	vb1 = s1p00r + optr_off[ntmp+1];
			va2 = tmp + 0x20;	vb2 = s1p00r + optr_off[ntmp+2];
			va3 = tmp + 0x30;	vb3 = s1p00r + optr_off[ntmp+3];
			va4 = tmp + 0x40;	vb4 = s1p00r + optr_off[ntmp+4];

			vc0 = tmp + 0x02;	vd0 = s1p00r + optr_off[ntmp+5];
			vc1 = tmp + 0x12;	vd1 = s1p00r + optr_off[ntmp+6];
			vc2 = tmp + 0x22;	vd2 = s1p00r + optr_off[ntmp+7];
			vc3 = tmp + 0x32;	vd3 = s1p00r + optr_off[ntmp+8];
			vc4 = tmp + 0x42;	vd4 = s1p00r + optr_off[ntmp+9];
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two,
				va0,va1,va2,va3,va4, vb0,vb1,vb2,vb3,vb4,
				vc0,vc1,vc2,vc3,vc4, vd0,vd1,vd2,vd3,vd4)
			tmp += 4;
		}

	  #else

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

		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, r00,r16,r32,r48,r64, s1p00r,s1p16r,s1p32r,s1p08r,s1p24r,\
										r02,r18,r34,r50,r66, s1p25r,s1p01r,s1p17r,s1p33r,s1p09r)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, r04,r20,r36,r52,r68, s1p10r,s1p26r,s1p02r,s1p18r,s1p34r,\
										r06,r22,r38,r54,r70, s1p35r,s1p11r,s1p27r,s1p03r,s1p19r)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, r08,r24,r40,r56,r72, s1p20r,s1p36r,s1p12r,s1p28r,s1p04r,\
										r10,r26,r42,r58,r74, s1p05r,s1p21r,s1p37r,s1p13r,s1p29r)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, r12,r28,r44,r60,r76, s1p30r,s1p06r,s1p22r,s1p38r,s1p14r,\
										r14,r30,r46,r62,r78, s1p15r,s1p31r,s1p07r,s1p23r,s1p39r)

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
		tmp = half_arr +  64;	// No lookup-tables used in avx-512; instead use opmasked conditional-doubling; 1st 32 slots hold outputs of wtsinit call
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

	  if(incr) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of incr
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

		i = (!j);
		addr = &prp_mult;
		tmp = s1p00r; tm1 = cy;
	  #ifndef USE_AVX512
		tm2 = cy+1;
	  #endif
		itmp = bjmodn;
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
		tm1 = s1p00r; tmp = cy; itmp = bjmodn;
		addr = &prp_mult;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];
			AVX_cmplx_carry_norm_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
			tm1 += 8; tmp += 1; itmp += 4; i = 0;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  }	// USE_AVX: (8-way or 4-way LOACC) or (4-way HIACC) ?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	  if(incr) {	// Have no specialized HIACC carry macro in ARM_V8_SIMD, so use 0-or-not-ness of incr
						// in lieu of (USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) is for ARMv8
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
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
		SSE2_cmplx_carry_fast_wtsinit(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, k0,k1,k2,k3, sse_bw,sse_n)

		i = (!j);
		addr = &prp_mult;
		tm1 = s1p00r; tmp = cy; tm2 = cy + 1; itmp = bjmodn;
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

		i = (!j);
		addr = &prp_mult;
		tm1 = s1p00r; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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

		tm1 = s1p00r; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; addr = cy; itmp = bjmodn;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			// Re-init weights every 4th macro invocation to keep errors under control:
			cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		}

	  } else {	// HiACC:

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; addr = cy; itmp = bjmodn;
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

/*...The radix-40 DIF pass is here:	*/

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

	/*...gather the needed data (40 64-bit complex, i.e. 80 64-bit reals) and do 8 radix-5 transforms...*/
		// Bytewise array of input-pointer offsets w.r.to the s1p00 base-pointer.
		// Radix-5 DFT inputs are (cyclic) with vec_dbl-pointer -= 16 (mod 80) between successive outputs:
		const uint8 iptr_off[RADIX] = {
			 0,64,48,32,16,
			 70,54,38,22,6,
			 60,44,28,12,76,
			 50,34,18,2,66,
			 40,24,8,72,56,
			 30,14,78,62,46,
			 20,4,68,52,36,
			 10,74,58,42,26};
		for(l = 0, tmp = r00, ntmp = 0; l < 4; l++, ntmp += 10)
		{
			// Output-ptrs [va* here] are regular-stride offsets of r00:
			va0 = tmp;			vb0 = s1p00r + iptr_off[ntmp  ];
			va1 = tmp + 0x10;	vb1 = s1p00r + iptr_off[ntmp+1];
			va2 = tmp + 0x20;	vb2 = s1p00r + iptr_off[ntmp+2];
			va3 = tmp + 0x30;	vb3 = s1p00r + iptr_off[ntmp+3];
			va4 = tmp + 0x40;	vb4 = s1p00r + iptr_off[ntmp+4];

			vc0 = tmp + 0x02;	vd0 = s1p00r + iptr_off[ntmp+5];
			vc1 = tmp + 0x12;	vd1 = s1p00r + iptr_off[ntmp+6];
			vc2 = tmp + 0x22;	vd2 = s1p00r + iptr_off[ntmp+7];
			vc3 = tmp + 0x32;	vd3 = s1p00r + iptr_off[ntmp+8];
			vc4 = tmp + 0x42;	vd4 = s1p00r + iptr_off[ntmp+9];
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two,
				vb0,vb1,vb2,vb3,vb4, va0,va1,va2,va3,va4,
				vd0,vd1,vd2,vd3,vd4, vc0,vc1,vc2,vc3,vc4)
			tmp += 4;
		}
	/*...and now do 5 radix-8 transforms: */
		const uint32 p_id2[5] = {054762310,045673201,032016754,067540123,001235476};	// Leading 0s = octal consts
		const uint8  p_od2[5] = {0,8,6,4,2};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00; l < 5; l++, tmp+=16) {
			i7 = p_id2[l]; i0 = i7&7; i1 = (i7>>3)&7; i2 = (i7>>6)&7; i3 = (i7>>9)&7; i4 = (i7>>12)&7; i5 = (i7>>15)&7; i6 = (i7>>18)&7; i7 = (i7>>21);
			addr = &a[j1+poff[p_od2[l]]];
			add0 = addr+pp07[i0]; add1 = addr+pp07[i1]; add2 = addr+pp07[i2]; add3 = addr+pp07[i3]; add4 = addr+pp07[i4]; add5 = addr+pp07[i5]; add6 = addr+pp07[i6]; add7 = addr+pp07[i7];
		  #ifdef USE_AVX2
			SSE2_RADIX8_DIF_0TWIDDLE(tmp,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two);
		  #else
			SSE2_RADIX8_DIF_0TWIDDLE(tmp,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
		  #endif
		}

	  #else

	/*...gather the needed data (40 64-bit complex) and do 8 radix-5 transforms...*/

		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, s1p00r,s1p32r,s1p24r,s1p16r,s1p08r, r00,r16,r32,r48,r64,\
										s1p35r,s1p27r,s1p19r,s1p11r,s1p03r, r02,r18,r34,r50,r66)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, s1p30r,s1p22r,s1p14r,s1p06r,s1p38r, r04,r20,r36,r52,r68,\
										s1p25r,s1p17r,s1p09r,s1p01r,s1p33r, r06,r22,r38,r54,r70)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, s1p20r,s1p12r,s1p04r,s1p36r,s1p28r, r08,r24,r40,r56,r72,\
										s1p15r,s1p07r,s1p39r,s1p31r,s1p23r, r10,r26,r42,r58,r74)
		SSE2_RADIX_05_DFT_0TWIDDLE_X2(xcc1,two, s1p10r,s1p02r,s1p34r,s1p26r,s1p18r, r12,r28,r44,r60,r76,\
										s1p05r,s1p37r,s1p29r,s1p21r,s1p13r, r14,r30,r46,r62,r78)

	/*...and now do 5 radix-8 transforms, swapping the t[48+i] <--> t[64+i] pairs to undo the last-2-outputs-swap in the RADIX_05_DFT macro:	*/
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

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
		#undef OFF5
		#undef OFF6
		#undef OFF7
	  #endif

	  #endif

	#else

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

	#endif	/* #ifdef USE_SSE2 */
	}

	jstart += nwt;
	jhi    += nwt;

	col += RADIX;
	co3 -= RADIX;
}	/* end for(int k=1; k <= khi; k++) */
