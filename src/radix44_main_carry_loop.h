/*******************************************************************************
*                                                                              *
*   (C) 1997-2017 by Ernst W. Mayer.                                           *
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

	/*...The radix-44 DIT pass is here:	*/

	#ifdef USE_SSE2

		// Radix-4 Outputs in SSE2 mode are temps 2*11*sz_vd bytes apart:
	  #ifdef USE_AVX512
		#define OFF	0x160*4
	  #elif defined(USE_AVX)
		#define OFF	0x160*2
	  #else
		#define OFF	0x160
	  #endif

	  #if COMPACT_OBJ

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 11 radix-4 transforms...*/
		/* Outputs in SSE2 modes are temps 2*11*16 = 22*16 = 0x160 bytes apart: */
		// Indices into above 4-elt table; each DFT-4 needs four 2-bit indices, thus gets 1 byte
		// Ex: 1st DFT-4 has add0-3 p01-multiple offsets 0,1,3,2; bit-reverse that to get p_idx[0] = 2310_4 = 0xb4:
		const uint8 p_id1[11] = {0xb4,0x1b,0x1b,0xe1,0xe1,0xe1,0x4e,0x4e,0x4e,0xb4,0xb4};
		const uint8 p_od1[11] = {0,7,3,10,6,2,9,5,1,8,4};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00; l < 11; l++, tmp+=2) {
			i3 = p_id1[l]; i0 = i3&3; i1 = (i3>>2)&3; i2 = (i3>>4)&3; i3 = (i3>>6);
			addr = &a[j1+poff[p_od1[l]]];
			add0 = addr+p0123[i0];	add1 = addr+p0123[i1];	add2 = addr+p0123[i2];	add3 = addr+p0123[i3];
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0,add1,add2,add3, tmp, OFF);
		}
		/*...and now do 4 radix-11 transforms: */
		// Radix-11 DFT outputs are (cyclic) with vec_dbl-pointer += 24 (mod 88) between successive outputs:
		const uint8 optr_off[RADIX] = {
			 0,24,48,72,8,32,56,80,16,40,64,
			22,46,70,6,30,54,78,14,38,62,86,
			44,68,4,28,52,76,12,36,60,84,20,
			66,2,26,50,74,10,34,58,82,18,42};
		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,*va9,*vaa,	// I-ptrs
		*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8,*vb9,*vba;	// O-ptrs
		for(l = 0, tmp = r00, ntmp = 0; l < 4; l++, ntmp += 11) {
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;		vb0 = s1p00 + optr_off[ntmp  ];
			va1 = tmp +  2;	vb1 = s1p00 + optr_off[ntmp+1];
			va2 = tmp +  4;	vb2 = s1p00 + optr_off[ntmp+2];
			va3 = tmp +  6;	vb3 = s1p00 + optr_off[ntmp+3];
			va4 = tmp +  8;	vb4 = s1p00 + optr_off[ntmp+4];
			va5 = tmp + 10;	vb5 = s1p00 + optr_off[ntmp+5];
			va6 = tmp + 12;	vb6 = s1p00 + optr_off[ntmp+6];
			va7 = tmp + 14;	vb7 = s1p00 + optr_off[ntmp+7];
			va8 = tmp + 16;	vb8 = s1p00 + optr_off[ntmp+8];
			va9 = tmp + 18;	vb9 = s1p00 + optr_off[ntmp+9];
			vaa = tmp + 20;	vba = s1p00 + optr_off[ntmp+10];
			SSE2_RADIX_11_DFT(
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,	/* inputs */
				ua0,	/* auxiliary-consts */
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba		/* outputs */
			);
			tmp += 22;
		}

	  #else

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 11 radix-4 transforms...*/
		/* Outputs in SSE2 modes are temps 2*11*16 = 22*16 = 0x160 bytes apart: */
		add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, OFF);
		add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, OFF);
		add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, OFF);
		add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, OFF);
		add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, OFF);
		add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, OFF);
		add2 = &a[j1+p36];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, OFF);
		add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, OFF);
		add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, OFF);
		add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, OFF);
		add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, OFF);

		/*...and now do 4 radix-11 transforms...*/
		SSE2_RADIX_11_DFT(r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a, ua0, s1p00,s1p0c,s1p18,s1p24,s1p04,s1p10,s1p1c,s1p28,s1p08,s1p14,s1p20)
		SSE2_RADIX_11_DFT(r0b,r0c,r0d,r0e,r0f,r10,r11,r12,r13,r14,r15, ua0, s1p0b,s1p17,s1p23,s1p03,s1p0f,s1p1b,s1p27,s1p07,s1p13,s1p1f,s1p2b)
		SSE2_RADIX_11_DFT(r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e,r1f,r20, ua0, s1p16,s1p22,s1p02,s1p0e,s1p1a,s1p26,s1p06,s1p12,s1p1e,s1p2a,s1p0a)
		SSE2_RADIX_11_DFT(r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b, ua0, s1p21,s1p01,s1p0d,s1p19,s1p25,s1p05,s1p11,s1p1d,s1p29,s1p09,s1p15)

	  #endif	// AVX or SSE2?

	#else	/* !USE_SSE2 */

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 11 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);
		/*...and now do 4 radix-11 transforms...*/
	  #ifdef LO_ADD	// cf. platform.h; this gives higher MUL count but improved accuracy
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_11_DFT_BASIC(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT_BASIC(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT_BASIC(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT_BASIC(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);
	  #else
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_11_DFT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);
	  #endif	// LO_ADD ?

	#endif

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need RADIX separate carries.	*/

	#ifdef USE_AVX

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];
		/* ptr to local storage for the doubled wtl,wtn terms: */
	  #ifdef USE_AVX512
		tmp = half_arr +  32;	// No lookup-tables used in avx-512; instead use opmasked conditional-doubling; 1st 32 slots hold outputs of wtsinit call
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

	  #ifdef LOACC

		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

		i = (!j);
		tmp = s1p00; tm1 = cy; tm2 = cy+1; itmp = bjmodn; itm2 = bjmodn+4;
		for(l = 0; l < RADIX>>3; l++) {
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l+l];
		  #ifdef USE_AVX512	// In AVX-512 mode, the 4 doubles base[0],baseinv[1],wts_mult[1],inv_mult[0] are in the d0-3 slots of the otherwie-unused sse2_rnd vec_dbl:
			AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04);
			tmp += 16; tm1 += 1;           itmp += 8;            i = 0;	// CY-ptr only advances 1 in AVX-512 mode, since all 8 dbl-carries fit in a single vec_dbl
		  #else	// USE_AVX:
			AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04);
			tmp += 16; tm1 += 2; tm2 += 2; itmp += 8; itm2 += 8; i = 0;
		  #endif
		}
		// And do a final 4-fold pass for 32-35:
		add0 = a + j1 + pfetch_dist + poff[l+l];
	   #ifdef USE_AVX512
		// AVX-512 mode calls this macro twice, with Call 2ptr-offsets fiddled as described in comments to that version of the macro:
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,0x000, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
		add0 += p04;				// prefetch of a + [prefetch offset] + p4,5,6,7
		tmp  = (double *)tmp +  4;	// Call 2 will handle the .d4-7 doubles of our 4 input zmm register-sized vector data
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,0x400, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);	// Call 2 wts-data pointers += 0x400
	   #else
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,i,     sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
	   #endif

	  #else	// USE_AVX: Hi-accuracy 4-way carry is the default:

		/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		tm1 = s1p00; tmp = cy; itmp = bjmodn;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];
			AVX_cmplx_carry_norm_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
			tm1 += 8; tmp += 1; itmp += 4; i = 0;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  #endif	// USE_AVX: (8-way or 4-way LOACC) or (4-way HIACC) ?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	  #ifdef LOACC

		uint32 k0,k1,k2,k3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;
		const uint32 *incr,inc_arr[] = {4,3,4};	// incr = 6+5 too large here, so divide 11 into 3 nearly-equal parts
		i = (!j);	// Need this to force 0-wod to be bigword
		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every few passes:
		incr = inc_arr;	// When RADIX/4 prime, fiddle the loop to allow e.g. 4+5+4 macro calls separated by reseeding steps.
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

			for(l = loop; l < loop+*incr; l++) {
				// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data:
				add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
				SSE2_cmplx_carry_fast_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03);
				tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
			}
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

		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		i = (!j);
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

		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p02,p03);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

	  #endif	// LOACC or HIACC?

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

	  #ifdef LOACC

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; addr = cy; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1    ],a[j2    ],*addr,*itmp,0); ++l; ++addr; ++itmp;
		cmplx_carry_fast_errcheck(a[j1+p01],a[j2+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_fast_errcheck(a[j1+p02],a[j2+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_fast_errcheck(a[j1+p03],a[j2+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			// Re-init weights every 4th macro invocatin to keep errors under control:
			cmplx_carry_norm_errcheck0(a[jt    ],a[jp    ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p01],a[jp+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p02],a[jp+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p03],a[jp+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}

	  #else	// Hi-accuracy is the default:

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; addr = cy; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1    ],a[j2    ],*addr,*itmp,0); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p01],a[j2+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p02],a[j2+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p03],a[j2+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			cmplx_carry_norm_errcheck (a[jt    ],a[jp    ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p01],a[jp+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p02],a[jp+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p03],a[jp+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}

	  #endif

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous)ghis assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

		/*...The radix-44 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #if COMPACT_OBJ

		/* Do 4 radix-11 transforms: */
		// Radix-11 DFT inputs are (cyclic) with vec_dbl-pointer -= 8 (mod 88) between successive outputs
		const uint8 iptr_off[RADIX] = {
			 0,80,72,64,56,48,40,32,24,16,8,
			66,58,50,42,34,26,18,10,2,82,74,
			44,36,28,20,12,4,84,76,68,60,52,
			22,14,6,86,78,70,62,54,46,38,30};
		for(l = 0, tmp = r00, ntmp = 0; l < 4; l++, ntmp += 11) {
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;		vb0 = s1p00 + iptr_off[ntmp  ];
			va1 = tmp +  2;	vb1 = s1p00 + iptr_off[ntmp+1];
			va2 = tmp +  4;	vb2 = s1p00 + iptr_off[ntmp+2];
			va3 = tmp +  6;	vb3 = s1p00 + iptr_off[ntmp+3];
			va4 = tmp +  8;	vb4 = s1p00 + iptr_off[ntmp+4];
			va5 = tmp + 10;	vb5 = s1p00 + iptr_off[ntmp+5];
			va6 = tmp + 12;	vb6 = s1p00 + iptr_off[ntmp+6];
			va7 = tmp + 14;	vb7 = s1p00 + iptr_off[ntmp+7];
			va8 = tmp + 16;	vb8 = s1p00 + iptr_off[ntmp+8];
			va9 = tmp + 18;	vb9 = s1p00 + iptr_off[ntmp+9];
			vaa = tmp + 20;	vba = s1p00 + iptr_off[ntmp+10];
			SSE2_RADIX_11_DFT(
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,	/* inputs */
				ua0,	/* auxiliary-consts */
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa		/* outputs */
			);
			tmp += 22;
		}
		// Indices into above 4-elt table; each DFT-4 needs four 2-bit indices, thus gets 1 byte
		// Ex: 1st DFT-4 has add0-3 p01-multiple offsets 0,1,2,3; bit-reverse that to get p_idx[0] = 3210_4 = 0xe4:
		const uint8 p_id2[11] = {0xe4,0xb1,0x1e,0xe4,0x4b,0xb1,0x1e,0xe4,0x4b,0xb1,0x1e};
		const uint8 p_od2[11] = {0,10,9,8,7,6,5,4,3,2,1};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00; l < 11; l++, tmp+=2) {
			i3 = p_id2[l]; i0 = i3&3; i1 = (i3>>2)&3; i2 = (i3>>4)&3; i3 = (i3>>6);
			addr = &a[j1+poff[p_od2[l]]];
			add0 = addr+p0123[i0];	add1 = addr+p0123[i1];	add2 = addr+p0123[i2];	add3 = addr+p0123[i3];
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0,add1,add2,add3, tmp, OFF);
		}

	  #else

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-11 transforms...*/
		SSE2_RADIX_11_DFT(s1p00,s1p28,s1p24,s1p20,s1p1c,s1p18,s1p14,s1p10,s1p0c,s1p08,s1p04, ua0, r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a)
		SSE2_RADIX_11_DFT(s1p21,s1p1d,s1p19,s1p15,s1p11,s1p0d,s1p09,s1p05,s1p01,s1p29,s1p25, ua0, r0b,r0c,r0d,r0e,r0f,r10,r11,r12,r13,r14,r15)
		SSE2_RADIX_11_DFT(s1p16,s1p12,s1p0e,s1p0a,s1p06,s1p02,s1p2a,s1p26,s1p22,s1p1e,s1p1a, ua0, r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e,r1f,r20)
		SSE2_RADIX_11_DFT(s1p0b,s1p07,s1p03,s1p2b,s1p27,s1p23,s1p1f,s1p1b,s1p17,s1p13,s1p0f, ua0, r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b)

		/*...and now do 11 radix-4 transforms...*/
		// Radix-4 Inputs in SIMD mode are temps 2*11*sz_vd bytes apart. Notice how the add* indices after the first row repeat with period 4:OFF:
		add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, OFF);
		add1 = &a[j1+p40];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, OFF);
		add3 = &a[j1+p36];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, OFF);
		add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, OFF);
		add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, OFF);
		add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, OFF);
		add3 = &a[j1+p20];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, OFF);
		add0 = &a[j1+p16];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, OFF);
		add2 = &a[j1+p12];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, OFF);
		add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, OFF);
		add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, OFF);

	  #endif

	  #undef OFF

	#else	/* !USE_SSE2 */

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-11 transforms...*/
	  #ifdef LO_ADD	// cf. platform.h; this gives higher MUL count but improved accuracy
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_11_DFT_BASIC(a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT_BASIC(a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT_BASIC(a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT_BASIC(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);
	  #else
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_11_DFT(a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT(a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT(a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);
	  #endif	// LO_ADD ?
		/*...and now do 11 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

	#endif	/* !USE_SSE2 */

	}	/* end for(j=_jstart; j < _jhi; j += 2) */

	jstart += nwt;
	jhi    += nwt;

	col += RADIX;
	co3 -= RADIX;
}	/* end for(k=1; k <= khi; k++) */