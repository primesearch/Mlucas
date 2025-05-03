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

	/*...The radix-48 DIT pass is here:	*/
	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4;
		OFF1 = 0x20;
		OFF2 = 0x40;
		OFF3 = 0x60;
		OFF4 = 0x80;
	  #elif defined(USE_AVX512)
		#define OFF1	0x20*4
		#define OFF2	0x40*4
		#define OFF3	0x60*4
		#define OFF4	0x80*4
	  #elif defined(USE_AVX)
		#define OFF1	0x20*2
		#define OFF2	0x40*2
		#define OFF3	0x60*2
		#define OFF4	0x80*2
	  #else
		#define OFF1	0x20
		#define OFF2	0x40
		#define OFF3	0x60
		#define OFF4	0x80
	  #endif

	#ifdef USE_SSE2

	  #if COMPACT_OBJ

		// Indices into 16-elt p01*[0-f] table; each DFT-16 needs sixteen 4-bit indices, thus needs a uint64.
		// Offsets: 	00,01,03,02,07,06,05,04,15,14,13,12,11,10,09,08, reverse order ==> 0x89abcdef45672310
		// Offsets: 32+	05,04,06,07,01,00,02,03,09,08,10,11,14,15,12,13, reverse order ==> 0xdcfeba8932017645
		// Offsets: 16+ 10,11,08,09,12,13,15,14,02,03,00,01,04,05,07,06, reverse order ==> 0x67541032efdc98ba
		const uint64 p_id1[3] = {0x89ABCDEF45672310ull,0xDCFEBA8932017645ull,0x67541032EFDC98BAull};	uint64 tmp64;
		const uint8  p_od1[3] = {0,8,4};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00r; l < 3; l++, tmp+=32) {
			tmp64 = p_id1[l];
			k0 =  tmp64     &15; k1 = (tmp64>> 4)&15; k2 = (tmp64>> 8)&15; k3 = (tmp64>>12)&15; k4 = (tmp64>>16)&15; k5 = (tmp64>>20)&15; k6 = (tmp64>>24)&15; k7 = (tmp64>>28)&15;
			k8 = (tmp64>>32)&15; k9 = (tmp64>>36)&15; ka = (tmp64>>40)&15; kb = (tmp64>>44)&15; kc = (tmp64>>48)&15; kd = (tmp64>>52)&15; ke = (tmp64>>56)&15; kf = (tmp64>>60);
			addr = &a[j1+poff[p_od1[l]]];
			po_kperm[0x0] = plo[k0];
			po_kperm[0x1] = plo[k1];
			po_kperm[0x2] = plo[k2];
			po_kperm[0x3] = plo[k3];
			po_kperm[0x4] = plo[k4];
			po_kperm[0x5] = plo[k5];
			po_kperm[0x6] = plo[k6];
			po_kperm[0x7] = plo[k7];
			po_kperm[0x8] = plo[k8];
			po_kperm[0x9] = plo[k9];
			po_kperm[0xa] = plo[ka];
			po_kperm[0xb] = plo[kb];
			po_kperm[0xc] = plo[kc];
			po_kperm[0xd] = plo[kd];
			po_kperm[0xe] = plo[ke];
			po_kperm[0xf] = plo[kf];
			SSE2_RADIX16_DIT_0TWIDDLE(
				addr,po_ptr,
				isrt2,two,
				tmp, OFF1,OFF2,OFF3,OFF4
			)
		}

	/*...and now do 16 radix-3 transforms. */
		// Bytewise array of output-pointer offsets w.r.to the s1p00 base-pointer:
		const uint8 optr_off[RADIX] = {
			 0,32,64,
			30,62,94,
			60,92,28,
			90,26,58,
			24,56,88,
			54,86,22,
			84,20,52,
			18,50,82,
			48,80,16,
			78,14,46,
			12,44,76,
			42,74,10,
			72, 8,40,
			 6,38,70,
			36,68, 4,
			66, 2,34};
		vec_dbl
		*va0,*va1,*va2,	// I-ptrs
		*vb0,*vb1,*vb2;	// O-ptrs
		vec_dbl
		*vc0,*vc1,*vc2,	// I-ptrs
		*vd0,*vd1,*vd2;	// O-ptrs
		for(l = 0, tmp = r00r, ntmp = 0; l <  8; l++, ntmp += 6)
		{
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;			vb0 = s1p00r + optr_off[ntmp  ];
			va1 = tmp + 0x20;	vb1 = s1p00r + optr_off[ntmp+1];
			va2 = tmp + 0x40;	vb2 = s1p00r + optr_off[ntmp+2];

			vc0 = tmp + 0x02;	vd0 = s1p00r + optr_off[ntmp+3];
			vc1 = tmp + 0x22;	vd1 = s1p00r + optr_off[ntmp+4];
			vc2 = tmp + 0x42;	vd2 = s1p00r + optr_off[ntmp+5];
			SSE2_RADIX_03_DFT_X2(cc1, va0,va1,va2, vb0,vb1,vb2, vc0,vc1,vc2, vd0,vd1,vd2)
			tmp += 4;
		}

	  #else

	// Macros use literal ostrides which are [1,2,3,4]-multiples of 1 vector-complex = 0x20 bytes for sse2, 0x40 bytes for avx
	// Offsets: 	00,01,03,02,07,06,05,04,15,14,13,12,11,10,09,08:
		addr = &a[j1    ];
		po_kperm[0x0] =   0;
		po_kperm[0x1] = p01;
		po_kperm[0x2] = p03;
		po_kperm[0x3] = p02;
		po_kperm[0x4] = p07;
		po_kperm[0x5] = p06;
		po_kperm[0x6] = p05;
		po_kperm[0x7] = p04;
		po_kperm[0x8] = p07 + p08;
		po_kperm[0x9] = p06 + p08;
		po_kperm[0xa] = p05 + p08;
		po_kperm[0xb] = p04 + p08;
		po_kperm[0xc] = p03 + p08;
		po_kperm[0xd] = p02 + p08;
		po_kperm[0xe] = p01 + p08;
		po_kperm[0xf] =       p08;
		SSE2_RADIX16_DIT_0TWIDDLE(addr,po_ptr, isrt2,two, r00r, OFF1,OFF2,OFF3,OFF4)

	// Offsets: 32+	05,04,06,07,01,00,02,03,09,08,10,11,14,15,12,13:
		addr = &a[j1+p32];
		po_kperm[0x0] = p05;
		po_kperm[0x1] = p04;
		po_kperm[0x2] = p06;
		po_kperm[0x3] = p07;
		po_kperm[0x4] = p01;
		po_kperm[0x5] =   0;
		po_kperm[0x6] = p02;
		po_kperm[0x7] = p03;
		po_kperm[0x8] = p01 + p08;
		po_kperm[0x9] =       p08;
		po_kperm[0xa] = p02 + p08;
		po_kperm[0xb] = p03 + p08;
		po_kperm[0xc] = p06 + p08;
		po_kperm[0xd] = p07 + p08;
		po_kperm[0xe] = p04 + p08;
		po_kperm[0xf] = p05 + p08;
		SSE2_RADIX16_DIT_0TWIDDLE(addr,po_ptr, isrt2,two, r16r, OFF1,OFF2,OFF3,OFF4)

	// Offsets: 16+ 10,11,08,09,12,13,15,14,02,03,00,01,04,05,07,06:
		addr = &a[j1+p16];
		po_kperm[0x0] = p02 + p08;
		po_kperm[0x1] = p03 + p08;
		po_kperm[0x2] =       p08;
		po_kperm[0x3] = p01 + p08;
		po_kperm[0x4] = p04 + p08;
		po_kperm[0x5] = p05 + p08;
		po_kperm[0x6] = p07 + p08;
		po_kperm[0x7] = p06 + p08;
		po_kperm[0x8] = p02;
		po_kperm[0x9] = p03;
		po_kperm[0xa] =   0;
		po_kperm[0xb] = p01;
		po_kperm[0xc] = p04;
		po_kperm[0xd] = p05;
		po_kperm[0xe] = p07;
		po_kperm[0xf] = p06;
		SSE2_RADIX16_DIT_0TWIDDLE(addr,po_ptr, isrt2,two, r32r, OFF1,OFF2,OFF3,OFF4)

		SSE2_RADIX_03_DFT_X2(cc1, r01r,r17r,r33r, s1p15r,s1p31r,s1p47r, r00r,r16r,r32r, s1p00r,s1p16r,s1p32r)
		SSE2_RADIX_03_DFT_X2(cc1, r03r,r19r,r35r, s1p45r,s1p13r,s1p29r, r02r,r18r,r34r, s1p30r,s1p46r,s1p14r)
		SSE2_RADIX_03_DFT_X2(cc1, r05r,r21r,r37r, s1p27r,s1p43r,s1p11r, r04r,r20r,r36r, s1p12r,s1p28r,s1p44r)
		SSE2_RADIX_03_DFT_X2(cc1, r07r,r23r,r39r, s1p09r,s1p25r,s1p41r, r06r,r22r,r38r, s1p42r,s1p10r,s1p26r)
		SSE2_RADIX_03_DFT_X2(cc1, r09r,r25r,r41r, s1p39r,s1p07r,s1p23r, r08r,r24r,r40r, s1p24r,s1p40r,s1p08r)
		SSE2_RADIX_03_DFT_X2(cc1, r11r,r27r,r43r, s1p21r,s1p37r,s1p05r, r10r,r26r,r42r, s1p06r,s1p22r,s1p38r)
		SSE2_RADIX_03_DFT_X2(cc1, r13r,r29r,r45r, s1p03r,s1p19r,s1p35r, r12r,r28r,r44r, s1p36r,s1p04r,s1p20r)
		SSE2_RADIX_03_DFT_X2(cc1, r15r,r31r,r47r, s1p33r,s1p01r,s1p17r, r14r,r30r,r46r, s1p18r,s1p34r,s1p02r)

	  #endif	// COMPACT_OBJ ?

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

	  if(incr) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of the low
						// inc_arr element to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
		AVX_cmplx_carry_fast_wtsinit_X16(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
	  #else
		const uint32 nloop = RADIX>>3;
		AVX_cmplx_carry_fast_wtsinit_X8 (add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
	  #endif
		i = (!j);
		addr = &prp_mult;
		tmp = s1p00r; tm1 = cy; itmp = bjmodn;
	  #ifndef USE_AVX512
		tm2 = cy+1; itm2 = bjmodn+4;
	  #endif
		for(l = 0; l < nloop; l++) {
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l+l];
		  #ifdef USE_AVX512
		   #ifdef CARRY_16_WAY
			AVX_cmplx_carry_fast_errcheck_X16(tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
			tmp += 32; tm1 += 2;           itmp += 16;           i = 0;
		   #else
			AVX_cmplx_carry_fast_errcheck_X8 (tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
			tmp += 16; tm1 += 1;           itmp +=  8;           i = 0;	// CY-ptr only advances 1 in AVX-512/CARRY_8_WAY mode, since all 8 dbl-carries fit in a single vec_dbl
		   #endif
		  #else	// USE_AVX:
			AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
			tmp += 16; tm1 += 2; tm2 += 2; itmp += 8; itm2 += 8; i = 0;
		  #endif
		}

	  } else {	// HiACC:

		/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		i = (!j);
		addr = &prp_mult;
		tm1 = s1p00r; tmp = cy; itmp = bjmodn;
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
		uint32 k0,k1,k2,k3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00r; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every few passes:
		for(loop = 0; loop < nloop; loop += incr)
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
			wtl     = wt0[    l  ];
			wtn     = wt0[nwtml  ]*scale;
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

			for(l = loop; l < loop+incr; l++) {
				// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data
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

	/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00r; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		addr = &prp_mult;
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

	/*...The radix-48 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif
	  #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x060;
		OFF2 = 0x0c0;
		OFF3 = 0x120;
		OFF4 = 0x180;
	  #elif defined(USE_AVX512)
		#define OFF1	0x060*4
		#define OFF2	0x0c0*4
		#define OFF3	0x120*4
		#define OFF4	0x180*4
	  #elif defined(USE_AVX)
		#define OFF1	0x060*2
		#define OFF2	0x0c0*2
		#define OFF3	0x120*2
		#define OFF4	0x180*2
	  #else
		#define OFF1	0x060
		#define OFF2	0x0c0
		#define OFF3	0x120
		#define OFF4	0x180
	  #endif

	  #if COMPACT_OBJ

		// Bytewise array of input-pointer offsets w.r.to the s1p00 base-pointer -
		// offsets decr -= 32 (mod 96) left-to-right within each row, and -= 6 (mod 96) downward in each column:
		const uint8 iptr_off[RADIX] = {
			 0,64,32,
			90,58,26,
			84,52,20,
			78,46,14,
			72,40, 8,
			66,34, 2,
			60,28,92,
			54,22,86,
			48,16,80,
			42,10,74,
			36, 4,68,
			30,94,62,
			24,88,56,
			18,82,50,
			12,76,44,
			 6,70,38};
		for(l = 0, tmp = r00r, ntmp = 0; l <  8; l++, ntmp += 6)
		{
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;			vb0 = s1p00r + iptr_off[ntmp  ];
			va1 = tmp + 0x02;	vb1 = s1p00r + iptr_off[ntmp+1];
			va2 = tmp + 0x04;	vb2 = s1p00r + iptr_off[ntmp+2];

			vc0 = tmp + 0x06;	vd0 = s1p00r + iptr_off[ntmp+3];
			vc1 = tmp + 0x08;	vd1 = s1p00r + iptr_off[ntmp+4];
			vc2 = tmp + 0x0a;	vd2 = s1p00r + iptr_off[ntmp+5];
			SSE2_RADIX_03_DFT_X2(cc1, vb0,vb1,vb2, va0,va1,va2, vd0,vd1,vd2, vc0,vc1,vc2)
			tmp += 12;
		}

		// Indices into 16-elt p01*[0-f] table; each DFT-16 needs sixteen 4-bit indices, thus needs a uint64.
		// Offsets: 	00,01,02,03,05,04,07,06,10,11,09,08,15,14,12,13, reverse order ==> 0xDCEF89BA67453210
		// Offsets: 32+	05,04,07,06,02,03,01,00,15,14,12,13,09,08,11,10, reverse order ==> 0xAB89DCEF01326745
		// Offsets: 16+ 10,11,09,08,15,14,12,13,05,04,07,06,02,03,01,00, reverse order ==> 0x01326745DCEF89BA
		const uint64 p_id2[3] = {0xDCEF89BA67453210ull,0xAB89DCEF01326745ull,0x01326745DCEF89BAull};
		for(l = 0, tmp = r00r; l < 3; l++) {
			tmp64 = p_id2[l];
			po_kperm[0x0] = plo[ tmp64     &15];
			po_kperm[0x1] = plo[(tmp64>> 4)&15];
			po_kperm[0x2] = plo[(tmp64>> 8)&15];
			po_kperm[0x3] = plo[(tmp64>>12)&15];
			po_kperm[0x4] = plo[(tmp64>>16)&15];
			po_kperm[0x5] = plo[(tmp64>>20)&15];
			po_kperm[0x6] = plo[(tmp64>>24)&15];
			po_kperm[0x7] = plo[(tmp64>>28)&15];
			po_kperm[0x8] = plo[(tmp64>>32)&15];
			po_kperm[0x9] = plo[(tmp64>>36)&15];
			po_kperm[0xa] = plo[(tmp64>>40)&15];
			po_kperm[0xb] = plo[(tmp64>>44)&15];
			po_kperm[0xc] = plo[(tmp64>>48)&15];
			po_kperm[0xd] = plo[(tmp64>>52)&15];
			po_kperm[0xe] = plo[(tmp64>>56)&15];
			po_kperm[0xf] = plo[(tmp64>>60)];
			addr = &a[j1+poff[p_od1[l]]];	// p_od2[] would be same, no need for a 2nd such array
			SSE2_RADIX16_DIF_0TWIDDLE(
				tmp,OFF1,OFF2,OFF3,OFF4,
				isrt2,two,
				addr,po_ptr
			);	tmp += 0x2;
		}

	  #else

		SSE2_RADIX_03_DFT_X2(cc1, s1p00r,s1p32r,s1p16r, r00r,r01r,r02r, s1p45r,s1p29r,s1p13r, r03r,r04r,r05r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p42r,s1p26r,s1p10r, r06r,r07r,r08r, s1p39r,s1p23r,s1p07r, r09r,r10r,r11r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p36r,s1p20r,s1p04r, r12r,r13r,r14r, s1p33r,s1p17r,s1p01r, r15r,r16r,r17r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p30r,s1p14r,s1p46r, r18r,r19r,r20r, s1p27r,s1p11r,s1p43r, r21r,r22r,r23r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p24r,s1p08r,s1p40r, r24r,r25r,r26r, s1p21r,s1p05r,s1p37r, r27r,r28r,r29r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p18r,s1p02r,s1p34r, r30r,r31r,r32r, s1p15r,s1p47r,s1p31r, r33r,r34r,r35r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p12r,s1p44r,s1p28r, r36r,r37r,r38r, s1p09r,s1p41r,s1p25r, r39r,r40r,r41r)
		SSE2_RADIX_03_DFT_X2(cc1, s1p06r,s1p38r,s1p22r, r42r,r43r,r44r, s1p03r,s1p35r,s1p19r, r45r,r46r,r47r)

	// istride of [3 vector-complex]*[1,2,3,4] = [1,2,3,4]*0x60 bytes for sse2, [1,2,3,4]*0xc0 bytes for avx
	// Offsets: 	00,01,02,03,05,04,07,06,10,11,09,08,15,14,12,13
		addr = &a[j1    ];
		po_kperm[0x0] =   0;
		po_kperm[0x1] = p01;
		po_kperm[0x2] = p02;
		po_kperm[0x3] = p03;
		po_kperm[0x4] = p05;
		po_kperm[0x5] = p04;
		po_kperm[0x6] = p07;
		po_kperm[0x7] = p06;
		po_kperm[0x8] = p02 + p08;
		po_kperm[0x9] = p03 + p08;
		po_kperm[0xa] = p01 + p08;
		po_kperm[0xb] =       p08;
		po_kperm[0xc] = p07 + p08;
		po_kperm[0xd] = p06 + p08;
		po_kperm[0xe] = p04 + p08;
		po_kperm[0xf] = p05 + p08;
		SSE2_RADIX16_DIF_0TWIDDLE(r00r,OFF1,OFF2,OFF3,OFF4, isrt2,two, addr,po_ptr)

	// Offsets: 32+	05,04,07,06,02,03,01,00,15,14,12,13,09,08,11,10
		addr = &a[j1+p32];
		po_kperm[0x0] = p05;
		po_kperm[0x1] = p04;
		po_kperm[0x2] = p07;
		po_kperm[0x3] = p06;
		po_kperm[0x4] = p02;
		po_kperm[0x5] = p03;
		po_kperm[0x6] = p01;
		po_kperm[0x7] =   0;
		po_kperm[0x8] = p07 + p08;
		po_kperm[0x9] = p06 + p08;
		po_kperm[0xa] = p04 + p08;
		po_kperm[0xb] = p05 + p08;
		po_kperm[0xc] = p01 + p08;
		po_kperm[0xd] =       p08;
		po_kperm[0xe] = p03 + p08;
		po_kperm[0xf] = p02 + p08;
		SSE2_RADIX16_DIF_0TWIDDLE(r01r,OFF1,OFF2,OFF3,OFF4, isrt2,two, addr,po_ptr)

	// Offsets: 16+ 10,11,09,08,15,14,12,13,05,04,07,06,02,03,01,00:
		addr = &a[j1+p16];
		po_kperm[0x0] = p02 + p08;
		po_kperm[0x1] = p03 + p08;
		po_kperm[0x2] = p01 + p08;
		po_kperm[0x3] =       p08;
		po_kperm[0x4] = p07 + p08;
		po_kperm[0x5] = p06 + p08;
		po_kperm[0x6] = p04 + p08;
		po_kperm[0x7] = p05 + p08;
		po_kperm[0x8] = p05;
		po_kperm[0x9] = p04;
		po_kperm[0xa] = p07;
		po_kperm[0xb] = p06;
		po_kperm[0xc] = p02;
		po_kperm[0xd] = p03;
		po_kperm[0xe] = p01;
		po_kperm[0xf] =   0;
		SSE2_RADIX16_DIF_0TWIDDLE(r02r,OFF1,OFF2,OFF3,OFF4, isrt2,two, addr,po_ptr)

	  #endif	// COMPACT_OBJ ?

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
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
}	/* end for(int k=1; k <= khi; k++) */
