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
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	/*...The radix-36 DIT pass is here:	*/

	#ifdef USE_SSE2

		/*...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do 9 radix-4 transforms...*/
		/* Outputs in SSE2 modes are temps 2*9*16 = 18*16 = 0x120 bytes apart: */
	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF = 0x120;
	  #elif defined(USE_AVX512)
		#define OFF	0x120*4
	   #elif defined(USE_AVX)
		#define OFF	0x120*2
	  #else
		#define OFF	0x120
	  #endif

	  #if COMPACT_OBJ

		// Indices into above 4-elt table; each DFT-4 needs four 2-bit indices, thus gets 1 byte
		// Ex: 1st DFT-4 has add0-3 p01-multiple offsets 0,1,3,2; bit-reverse that to get p_idx[0] = 2310_4 = 0xb4:
		const uint8 p_id1[9] = {0xb4,0xe1,0x4e,0x1b,0xe1,0xb4,0x1b,0x4e,0xb4};
		const uint8 p_od1[9] = {0,2,1,7,6,8,3,5,4};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00; l < 9; l++, tmp+=2) {
			i3 = p_id1[l]; i0 = i3&3; i1 = (i3>>2)&3; i2 = (i3>>4)&3; i3 = (i3>>6);
			addr = &a[j1+poff[p_od1[l]]];
			add0 = addr+p0123[i0];	add1 = addr+p0123[i1];	add2 = addr+p0123[i2];	add3 = addr+p0123[i3];
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0,add1,add2,add3, tmp, OFF);
		}
		/*...and now do 4 radix-9 transforms: */
		// Radix-9 DFT outputs are (cyclic) with vec_dbl-pointer -= 8 (mod 72) between successive outputs:
		const uint8 optr_off[RADIX] = {
			 0,64,56,48,40,32,24,16,8,
			54,46,38,30,22,14,6,70,62,
			36,28,20,12,4,68,60,52,44,
			18,10,2,66,58,50,42,34,26};
		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,	// I-ptrs
		*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8;	// O-ptrs
	   #ifdef USE_AVX2
		// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
		tm1 = (vec_dbl *)rad9_iptr;	// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity
		tm2 = (vec_dbl *)rad9_optr;
		for(l = 0, tmp = r00, ntmp = 0; l < 2; l++, ntmp += 18) {
	   #else
		for(l = 0, tmp = r00, ntmp = 0; l < 4; l++, ntmp += 9) {
	   #endif
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
		   #ifdef USE_AVX2
			// Pointer patterns here same as for DIF, just need to swap I/O by reversing order of tm1,tm2 --> tm2,tm1 in macro arglists:
			rad9_optr[0] = tmp + 18;	rad9_iptr[0] = s1p00 + optr_off[ntmp+ 9];
			rad9_optr[1] = tmp + 20;	rad9_iptr[1] = s1p00 + optr_off[ntmp+10];
			rad9_optr[2] = tmp + 22;	rad9_iptr[2] = s1p00 + optr_off[ntmp+11];
			rad9_optr[3] = tmp + 24;	rad9_iptr[3] = s1p00 + optr_off[ntmp+12];
			rad9_optr[4] = tmp + 26;	rad9_iptr[4] = s1p00 + optr_off[ntmp+13];
			rad9_optr[5] = tmp + 28;	rad9_iptr[5] = s1p00 + optr_off[ntmp+14];
			rad9_optr[6] = tmp + 30;	rad9_iptr[6] = s1p00 + optr_off[ntmp+15];
			rad9_optr[7] = tmp + 32;	rad9_iptr[7] = s1p00 + optr_off[ntmp+16];
			rad9_optr[8] = tmp + 34;	rad9_iptr[8] = s1p00 + optr_off[ntmp+17];
			SSE2_RADIX_09_DIT_X2(
				va0,va1,va2,va3,va4,va5,va6,va7,va8,	/* inputs  1 */
				cc1,two,	/* auxiliary-consts */
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,	/* outputs 1 */
				tm2,tm1									/* in/outputs 2 */
			);
			tmp += 36;
		   #else
			SSE2_RADIX_09_DIT(
				va0,va1,va2,va3,va4,va5,va6,va7,va8,	/* inputs  1 */
				cc1,	/* auxiliary-consts */
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8		/* outputs 1 */
			);
			tmp += 18;
		   #endif
		}

	  #else

		add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, OFF)
		add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, OFF)
		add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, OFF)
		add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, OFF)
		add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, OFF)
		add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, OFF)
		add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, OFF)
		add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, OFF)
		add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0g, OFF)

		/* Radix-9 DFT uses adjacent temps, i.e. stride = 2*16 bytes: */
	   #ifdef USE_AVX2
		// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
		tm1 = (vec_dbl *)rad9_iptr;	// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity
		tm2 = (vec_dbl *)rad9_optr;
		// Pointer patterns here same as for DIF, just need to swap I/O by reversing order of tm1,tm2 --> tm2,tm1 in macro arglists:
		rad9_iptr[0] = s1p27; rad9_iptr[1] = s1p23; rad9_iptr[2] = s1p19; rad9_iptr[3] = s1p15; rad9_iptr[4] = s1p11; rad9_iptr[5] = s1p07; rad9_iptr[6] = s1p03; rad9_iptr[7] = s1p35; rad9_iptr[8] = s1p31;
		rad9_optr[0] = r10; rad9_optr[1] = r12; rad9_optr[2] = r14; rad9_optr[3] = r16; rad9_optr[4] = r18; rad9_optr[5] = r1a; rad9_optr[6] = r1c; rad9_optr[7] = r1e; rad9_optr[8] = r1g;
		SSE2_RADIX_09_DIT_X2(r00,r02,r04,r06,r08,r0a,r0c,r0e,r0g, cc1,two, s1p00,s1p32,s1p28,s1p24,s1p20,s1p16,s1p12,s1p08,s1p04,
				tm2,tm1)  // r10,r12,r14,r16,r18,r1a,r1c,r1e,r1g, cc1,     s1p27,s1p23,s1p19,s1p15,s1p11,s1p07,s1p03,s1p35,s1p31r)

		rad9_iptr[0] = s1p09; rad9_iptr[1] = s1p05; rad9_iptr[2] = s1p01; rad9_iptr[3] = s1p33; rad9_iptr[4] = s1p29; rad9_iptr[5] = s1p25; rad9_iptr[6] = s1p21; rad9_iptr[7] = s1p17; rad9_iptr[8] = s1p13;
		rad9_optr[0] = r30; rad9_optr[1] = r32; rad9_optr[2] = r34; rad9_optr[3] = r36; rad9_optr[4] = r38; rad9_optr[5] = r3a; rad9_optr[6] = r3c; rad9_optr[7] = r3e; rad9_optr[8] = r3g;
		SSE2_RADIX_09_DIT_X2(r20,r22,r24,r26,r28,r2a,r2c,r2e,r2g, cc1,two, s1p18,s1p14,s1p10,s1p06,s1p02,s1p34,s1p30,s1p26,s1p22,
				tm2,tm1)  // r30,r32,r34,r36,r38,r3a,r3c,r3e,r3g, cc1,     s1p09,s1p05,s1p01,s1p33,s1p29,s1p25,s1p21,s1p17,s1p13r)
	   #else
		SSE2_RADIX_09_DIT(r00,r02,r04,r06,r08,r0a,r0c,r0e,r0g, cc1,     s1p00,s1p32,s1p28,s1p24,s1p20,s1p16,s1p12,s1p08,s1p04)
		SSE2_RADIX_09_DIT(r10,r12,r14,r16,r18,r1a,r1c,r1e,r1g, cc1,     s1p27,s1p23,s1p19,s1p15,s1p11,s1p07,s1p03,s1p35,s1p31)
		SSE2_RADIX_09_DIT(r20,r22,r24,r26,r28,r2a,r2c,r2e,r2g, cc1,     s1p18,s1p14,s1p10,s1p06,s1p02,s1p34,s1p30,s1p26,s1p22)
		SSE2_RADIX_09_DIT(r30,r32,r34,r36,r38,r3a,r3c,r3e,r3g, cc1,     s1p09,s1p05,s1p01,s1p33,s1p29,s1p25,s1p21,s1p17,s1p13)
	   #endif

	  #endif	// COMPACT_OBJ ?

	#else	// Scalar-double mode:

		/*...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do 9 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, rt,it);
		/*...and now do 4 radix-9 transforms...*/
		tptr = t;
	  #if 1	// 2-level radix-3-triplet radix-9 DFT less accurate but faster than _FMA-named version (used strictly for opcount-costing-out an FMA-ized option) below
		jt = j1    ; jp = j2    ;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],rt,it,re);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],rt,it,re);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],rt,it,re);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],rt,it,re);
	  #else	// _FMA version of DIT uses different output ordering, thus yields different operm for radix-36 DIT:
		jt = j1    ; jp = j2    ;	RADIX_09_DIT_FMA(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p32],a[jp+p32],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28], c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIT_FMA(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p24],a[jp+p24],a[jt+p32],a[jp+p32],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16], c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIT_FMA(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p32],a[jp+p32],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt    ],a[jp    ],a[jt+p08],a[jp+p08], c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIT_FMA(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, a[jt+p08],a[jp+p08],a[jt+p16],a[jp+p16],a[jt+p24],a[jp+p24],a[jt+p32],a[jp+p32],a[jt+p04],a[jp+p04],a[jt+p12],a[jp+p12],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt    ],a[jp    ], c,s,c2,s2,c3,s3,c4,s4);
	  #endif

	#endif	// SIMD or not?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need RADIX separate carries.	*/

		// Check if current index-interval contains the target index for rotated-residue carry injection.
		// In data-init we set target_idx = -1 on wraparound-carry mini-pass, so if() only taken on full pass:
		if(target_idx == j) {
		#ifdef USE_SSE2
			addr = (double *)s1p00 + target_set;
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

	  if(incr) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of incr
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

		i = (!j);
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy; tm2 = cy+1; itmp = bjmodn;
	   #ifndef USE_AVX512
		itm2 = bjmodn+4;
	   #endif
		for(l = 0; l < RADIX>>3; l++) {
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l+l];
		  #ifdef USE_AVX512	// In AVX-512 mode, the 4 doubles base[0],baseinv[1],wts_mult[1],inv_mult[0] are in the d0-3 slots of the otherwie-unused sse2_rnd vec_dbl:
			AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
			tmp += 16; tm1 += 1;           itmp += 8;            i = 0;	// CY-ptr only advances 1 in AVX-512 mode, since all 8 dbl-carries fit in a single vec_dbl
		  #else	// USE_AVX:
			AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03,p04, addr);
			tmp += 16; tm1 += 2; tm2 += 2; itmp += 8; itm2 += 8; i = 0;
		  #endif
		}
		// And do a final 4-fold pass for 32-35:
		add0 = a + j1 + pfetch_dist + poff[l+l];
	   #ifdef USE_AVX512
	   #ifdef X4_ZMM
		// For fused X4x2 version, place a copy of bjmodn[0:3] in bjmodn[4:7] slots and increment by 7*bw (mod n):
		MOD_ADD32(((struct uint32x8*)itmp)->d0,bw_x7_modn,n, ((struct uint32x8*)itmp)->d4);
		MOD_ADD32(((struct uint32x8*)itmp)->d1,bw_x7_modn,n, ((struct uint32x8*)itmp)->d5);
		MOD_ADD32(((struct uint32x8*)itmp)->d2,bw_x7_modn,n, ((struct uint32x8*)itmp)->d6);
		MOD_ADD32(((struct uint32x8*)itmp)->d3,bw_x7_modn,n, ((struct uint32x8*)itmp)->d7);
		AVX_cmplx_carry_fast_errcheck_X4_ZMM(tmp, tm1, itmp, half_arr,0x000, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
	   #else
		// AVX-512 mode calls this macro twice, with Call 2ptr-offsets fiddled as described in comments to that version of the macro:
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,0x000, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
/*
=== 4-WAY ===
Input:
*tmp = {-516770855357715.75, 190672422644522.5, -211956860295766.72, -86904961926692.984, -55584894075895.883, -71126886553534.016, -206117597153738.62, 1097001158675664.9}
*tm1 = 0 x 8
*(struct uint32x8*)itmp = {d0 = 262144, d1 = 196608, d2 = 131072, d3 = 65536, d4 = 0, d5 = 0, d6 = 0, d7 = 0}
Call 1 Output:
*tmp = {-218730.49378017406, -48056.699550186735, -93280.269136510979, 52231.276634893475, [d4-7 same as above]}
*tm1 = {d0 = -78, d1 = 1402, d2 = -650, d3 = -149, d4 = 0, d5 = 0, d6 = 0, d7 = 0};
*itmp = {d0 = 1206488, d1 = 1140952, d2 = 1075416, d3 = 1009880, [d4-7 = 0]}
Call 2 Output:
*tmp = {-218730.49378017406, -48056.699550186735, -93280.269136510979, 52231.276634893475, -39592.837455274152, 102333.6974532048, 86549.392108867847, -268614.23201552214}
*tm1 = {d0 = -365, d1 = 1336, d2 = -256, d3 = 314, d4 = 0, d5 = 0, d6 = 0, d7 = 0}
*itmp = {d0 = 2150832, d1 = 2085296, d2 = 2019760, d3 = 1954224, d4 = 0, d5 = 0, d6 = 0, d7 = 0}

=== 8-WAY ===
tmp[0] = {-516770855357715.75, 190672422644522.5, -211956860295766.72, -86904961926692.984, -55584894075895.883, -71126886553534.016, -206117597153738.62, 1097001158675664.9}
tmp[7] = {815558989568825.25, 166479510996781, 602333754119585.88, -184329755184470.81, 118987226527178.23, -581635903644156.12, -75851190430049.625, -103403585708466.94}
	Set bit in k1 if sw < bjmodn[0:3]: sw = n - (p%n) = 2359296 - (43765019 % 2359296) = 1061605, which is < all 4 imtp[4-7] values below,
	thus
*tm1 = 0 x 8
*(struct uint32x8*)itmp = {d0 = 262144, d1 = 196608, d2 = 131072, d3 = 65536, d4 = 2268093, d5 = 2202557, d6 = 2137021, d7 = 2071485}
After rcol-carry-estimate step:
$zmm3.v8_int32 = {262144, 196608, 131072, 65536, 1206488, 1140952, 1075416, 1009880} d4-7 now match above Call 1 Output ones
$zmm1.v8_double = {0, 0, 0, 0, 1284, 257, 912, -274}, d4-7 nowhere near the correct {-78,1402,-650,-149} values!
Output:
*/
		add0 += p04;				// prefetch of a + [prefetch offset] + p4,5,6,7
		tmp  = (vec_dbl *)((double *)tmp +  4);	// Call 2 will handle the .d4-7 doubles of our 4 input zmm register-sized vector data
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,0x800, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);	// Call 2 wts-data pointers += 0x400
	   #endif
/*
IDEA: In avx-512 build, fuse the above 2 _X4 calls into 1, in which we use full-width zmm-registers rather than half-width ymm.
Added "AVX debug" rng-init to radix36*c, examined data:

Args: AVX_cmplx_carry_fast_errcheck_4(data,cy,bjmod_0, half_arr,doff, sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, prp_mult)\
vec_dbl*ptr := ((void *)half_arr + doff) = ((vec_dbl*)half_arr + [0,32]) is start of wt/wtinv block for current dataset
For the 8 data-quartets (A.re,A.im,B.re,B.im,C.re,C.im,D.re,D.im):
 	Residue-data in data + 0-7
	wt_fwd in ptr + 0,4,8,...
	wt_inv in ptr + 2,6,10,... [why 2-strides between wt_fwd/inv ... perhaps to allow for 16-way carry macros?]
	maxerr in (vec_dbl*)half_arr - 2
	cy, bjmodn data in the associated ptrs
Call 1:
  In:
	data: {d0 = -4037272307482.1543, d1 = 1489628301910.332, d2 = -1655912971060.6775, d3 = -678945015052.28894, d4 = -434256984967.93658, d5 = -555678801199.4845, d6 = -1610293727763.583, d7 = 8570321552153.6318}
	cy	: {d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0, d6 = 0, d7 = 0}
	bjmod {d0 = 262144, d1 = 196608, d2 = 131072, d3 = 65536
	p *(struct uint32x8*)itmp = {d0 = 262144, d1 = 196608, d2 = 131072, d3 = 65536, d4 = 0, d5 = 0, d6 = 0, d7 = 0}
	half_arr + ...:
		0	{d0 = 1.8517494245745807, d1 = 1.8877486253633868, d2 = 1.9244476737882898, d3 = 1.9618601753378295, d4 = 1.9999999999999996, d5 = 1.0194406437021448, d6 = 1.0392592260318434, d7 = 1.0594630943592951}
		4	{d0 = 1.2647530720435787, d1 = 1.2893406858883709, d2 = 1.3144062987734058, d3 = 1.3399592033077143, d4 = 1.3660088727546293, d5 = 1.3925649645438212, d6 = 1.4196373238516071, d7 = 1.4472359872508724}
		...
		28	{d0 = 1.0271571120694201, d1 = 1.0471257075112859, d2 = 1.0674825053023689, d3 = 1.0882350523462252, d4 = 1.1093910422630731, d5 = 1.1309583182420606, d6 = 1.1529448759489815, d7 = 1.1753588664905192}
		32	{d0 = 1.4031057287998494, d1 = 1.4303830073498853, d2 = 1.4581905737533771, d3 = 1.4865387371475425, d4 = 1.515438007085864, d5 = 1.5448990974343091, d6 = 1.5749329303432944, d7 = 1.6055506402968731}
		2	{d0 = 4.5778899251823675e-07, d1 = 4.4905899656477821e-07, d2 = 4.4049548086880288e-07, d3 = 4.3209527066639561e-07, d4 = 4.2385525173611109e-07, d5 = 8.3154473848886652e-07, d6 = 8.1568725322650911e-07, d7 = 8.0013216881789595e-07}
		6	{d0 = 6.7025771449797521e-07, d1 = 6.5747595864326529e-07, d2 = 6.4493794975214226e-07, d3 = 6.3263903959137917e-07, d4 = 6.2057466856914965e-07, d5 = 6.0874036404464391e-07, d6 = 5.9713173866991963e-07, d7 = 5.8574448876337593e-07}
		...
		30	{d0 = 8.2529779866327764e-07, d1 = 8.0955944199573126e-07, d2 = 7.9412121441007049e-07, d3 = 7.7897739247101577e-07, d4 = 7.641223618886962e-07, d5 = 7.4955061543725714e-07, d6 = 7.3525675091315763e-07, d7 = 7.2123546913240532e-07}
		34	{d0 = 6.041672313584761e-07, d1 = 5.9264581522315561e-07, d2 = 5.8134411148346563e-07, d3 = 5.7025793024328346e-07, d4 = 5.5938316150743785e-07, d5 = 5.4871577365800598e-07, d6 = 5.3825181195966429e-07, d7 = 5.2798739709354578e-07}
  Out:
	data: {d0 = -24463.461648054785, d1 = -95907.80709410111, d2 = -136076.03802495691, d3 = -181309.34733751958, d4 = -434256984967.93658, d5 = -555678801199.4845, d6 = -1610293727763.583, d7 = 8570321552153.6318}
	cy	: {d0 = -1, d1 = 11, d2 = -5, d3 = -1, d4 = 0, d5 = 0, d6 = 0, d7 = 0}
	bjmod {d0 = 1206488, d1 = 1140952, d2 = 1075416, d3 = 1009880, d4 = 0, d5 = 0, d6 = 0, d7 = 0}
	p *(struct uint32x8*)itmp = {d0 = 1206488, d1 = 1140952, d2 = 1075416, d3 = 1009880, d4 = 0, d5 = 0, d6 = 0, d7 = 0}
Call 2:
  In:
	data: [above data.d4-7]
	cy	: [above cy  .d0-3]
	bjmod [above bjmod.d0-3]
  Out:
	data:
	cy	:
	bjmod
***No Point Going Further*** -
Forgot that the two succeeding _X4 calls have data dependency, bjmodn/carries output by Call 1 are inputs to Call 2!
To make this work would need to quick-estimate carries out of Call 1, feed those as input-carries to the Call 2 data,
then at end of fused _X4x2 macro fold the differences between the actual Call 1 o-carries (i.e. the ones resulting from
the full-length Call 1 carry-propagation sequence) and the estimated ones into the starting data (A.re) for Call 2,
for which we could omit the usual fwd-weighting following the carry step until after said folding-in has been done.

Will also be doing lots of 256-bit half-register merges ... data quartets for Call 1|2 in low|high half of zmm.
Options -- assume Call 1|2 data in ymm0,ymm1, will concatenate into zmm0:
															Latency/thruput:
	Instruction sequence:					KNL											Skylake-X	Comment:
vpslldq 32,zmm0,zmm0						2/1												1/1		n/a on k1om; use vpermf32x4 78,zmm0,zmm0
valignd  8,zmm1,zmm0,zmm0					3-6/1											3/1		k1om has just valignd (not q), so use for both architectures

vshuff64x2 0b01000100,zmm1,zmm0,zmm0		4-7/2											3/1

vinsertf64x4 1,ymm1,zmm0,zmm0	 3-6/1 for y,z,z, 7/1 for m256,z,z			3/1 for y,z,z, 5/0.5 for m256,z,z	<=== BEST ===
*/
	   #else
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,i,     sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
	   #endif

	  } else {	// HiACC:

		/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		i = (!j);
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy; itmp = bjmodn;
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
		uint32 k0,k1,k2,k3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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
		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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

		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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

	/*...The radix-36 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #if COMPACT_OBJ

		/* Do 4 radix-9 transforms: */
		// Radix-9 DFT inputs can use same optr_off[] perm-index array as DIT:
	   #ifdef USE_AVX2
		// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
		tm1 = (vec_dbl *)rad9_iptr;	// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity
		tm2 = (vec_dbl *)rad9_optr;
		for(l = 0, tmp = r00, ntmp = 0; l < 2; l++, ntmp += 18) {
	   #else
		for(l = 0, tmp = r00, ntmp = 0; l < 4; l++, ntmp += 9) {
	   #endif
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
		   #ifdef USE_AVX2
			// Pointer patterns here same as for DIF, just need to swap I/O by reversing order of tm1,tm2 --> tm2,tm1 in macro arglists:
			rad9_optr[0] = tmp + 18;	rad9_iptr[0] = s1p00 + optr_off[ntmp+ 9];
			rad9_optr[1] = tmp + 20;	rad9_iptr[1] = s1p00 + optr_off[ntmp+10];
			rad9_optr[2] = tmp + 22;	rad9_iptr[2] = s1p00 + optr_off[ntmp+11];
			rad9_optr[3] = tmp + 24;	rad9_iptr[3] = s1p00 + optr_off[ntmp+12];
			rad9_optr[4] = tmp + 26;	rad9_iptr[4] = s1p00 + optr_off[ntmp+13];
			rad9_optr[5] = tmp + 28;	rad9_iptr[5] = s1p00 + optr_off[ntmp+14];
			rad9_optr[6] = tmp + 30;	rad9_iptr[6] = s1p00 + optr_off[ntmp+15];
			rad9_optr[7] = tmp + 32;	rad9_iptr[7] = s1p00 + optr_off[ntmp+16];
			rad9_optr[8] = tmp + 34;	rad9_iptr[8] = s1p00 + optr_off[ntmp+17];
			SSE2_RADIX_09_DIF_X2(
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,	/* inputs  1 */
				cc1,two,	/* auxiliary-consts */
				va0,va1,va2,va3,va4,va5,va6,va7,va8,	/* outputs 1 */
				tm1,tm2									/* in/outputs 2 */
			);
			tmp += 36;
		   #else
			SSE2_RADIX_09_DIF(
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,	/* inputs  1 */
				cc1,	/* auxiliary-consts */
				va0,va1,va2,va3,va4,va5,va6,va7,va8		/* outputs 1 */
			);
			tmp += 18;
		   #endif
		}
		// Indices into above 4-elt table; each DFT-4 needs four 2-bit indices, thus gets 1 byte
		// Ex: 1st DFT-4 has add0-3 p01-multiple offsets 0,1,3,2; bit-reverse that to get p_idx[0] = 2310_4 = 0xb4:
		const uint8 p_id2[9] = {0xb4,0xb4,0x4e,0xe1,0x1b,0xb4,0x4e,0xe1,0x1b};
		const uint8 p_od2[9] = {0,8,5,2,7,4,1,6,3};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00; l < 9; l++, tmp+=2) {
			i3 = p_id2[l]; i0 = i3&3; i1 = (i3>>2)&3; i2 = (i3>>4)&3; i3 = (i3>>6);
			addr = &a[j1+poff[p_od2[l]]];
			add0 = addr+p0123[i0];	add1 = addr+p0123[i1];	add2 = addr+p0123[i2];	add3 = addr+p0123[i3];
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0,add1,add2,add3, tmp, OFF);
		}

	  #else

		/* Radix-9 DFT uses adjacent temps, i.e. stride = 2*16 bytes: */
	   #ifdef USE_AVX2
		// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
		tm1 = (vec_dbl *)rad9_iptr;	// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity
		tm2 = (vec_dbl *)rad9_optr;
		rad9_iptr[0] = s1p27; rad9_iptr[1] = s1p23; rad9_iptr[2] = s1p19; rad9_iptr[3] = s1p15; rad9_iptr[4] = s1p11; rad9_iptr[5] = s1p07; rad9_iptr[6] = s1p03; rad9_iptr[7] = s1p35; rad9_iptr[8] = s1p31;
		rad9_optr[0] = r10; rad9_optr[1] = r12; rad9_optr[2] = r14; rad9_optr[3] = r16; rad9_optr[4] = r18; rad9_optr[5] = r1a; rad9_optr[6] = r1c; rad9_optr[7] = r1e; rad9_optr[8] = r1g;
		SSE2_RADIX_09_DIF_X2(s1p00,s1p32,s1p28,s1p24,s1p20,s1p16,s1p12,s1p08,s1p04, cc1,two, r00,r02,r04,r06,r08,r0a,r0c,r0e,r0g,
				tm1,tm2)  // s1p27,s1p23,s1p19,s1p15,s1p11,s1p07,s1p03,s1p35,s1p31, cc1, r10,r12,r14,r16,r18,r1a,r1c,r1e,r1g)

		rad9_iptr[0] = s1p09; rad9_iptr[1] = s1p05; rad9_iptr[2] = s1p01; rad9_iptr[3] = s1p33; rad9_iptr[4] = s1p29; rad9_iptr[5] = s1p25; rad9_iptr[6] = s1p21; rad9_iptr[7] = s1p17; rad9_iptr[8] = s1p13;
		rad9_optr[0] = r30; rad9_optr[1] = r32; rad9_optr[2] = r34; rad9_optr[3] = r36; rad9_optr[4] = r38; rad9_optr[5] = r3a; rad9_optr[6] = r3c; rad9_optr[7] = r3e; rad9_optr[8] = r3g;
		SSE2_RADIX_09_DIF_X2(s1p18,s1p14,s1p10,s1p06,s1p02,s1p34,s1p30,s1p26,s1p22, cc1,two, r20,r22,r24,r26,r28,r2a,r2c,r2e,r2g,
				tm1,tm2)  // s1p09,s1p05,s1p01,s1p33,s1p29,s1p25,s1p21,s1p17,s1p13, cc1, r30,r32,r34,r36,r38,r3a,r3c,r3e,r3g)
	   #else
		SSE2_RADIX_09_DIF(s1p00,s1p32,s1p28,s1p24,s1p20,s1p16,s1p12,s1p08,s1p04, cc1, r00,r02,r04,r06,r08,r0a,r0c,r0e,r0g)
		SSE2_RADIX_09_DIF(s1p27,s1p23,s1p19,s1p15,s1p11,s1p07,s1p03,s1p35,s1p31, cc1, r10,r12,r14,r16,r18,r1a,r1c,r1e,r1g)
		SSE2_RADIX_09_DIF(s1p18,s1p14,s1p10,s1p06,s1p02,s1p34,s1p30,s1p26,s1p22, cc1, r20,r22,r24,r26,r28,r2a,r2c,r2e,r2g)
		SSE2_RADIX_09_DIF(s1p09,s1p05,s1p01,s1p33,s1p29,s1p25,s1p21,s1p17,s1p13, cc1, r30,r32,r34,r36,r38,r3a,r3c,r3e,r3g)
	   #endif

		/* Outputs in SSE2 modes are temps 2*9*16 = 18*16 = 0X120 bytes apart: */
		add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, OFF)
		add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, OFF)
		add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, OFF)
		add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, OFF)
		add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, OFF)
		add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, OFF)
		add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, OFF)
		add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, OFF)
		add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0g, OFF)

	  #endif	// COMPACT_OBJ ?

	#else	/* !USE_SSE2 */

		/*...gather the needed data (36 64-bit complex, i.e. 72 64-bit reals) and do 4 radix-9 transforms...*/
	  #if 1	// 2-level radix-3-triplet radix-9 DFT less accurate but faster than _FMA-named version (used strictly for opcount-costing-out an FMA-ized option) below
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_09_DIF(a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, rt,it,re);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIF(a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, rt,it,re);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIF(a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, rt,it,re);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIF(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, rt,it,re);
		/*...and now do 9 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
	  #else	// _FMA version of DIF uses different output ordering, thus yields different operm for radix-36 DIT:
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_09_DIF_FMA(a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p03; jp = j2+p03;	RADIX_09_DIF_FMA(a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p02; jp = j2+p02;	RADIX_09_DIF_FMA(a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, c,s,c2,s2,c3,s3,c4,s4);	tptr += 9;
		jt = j1+p01; jp = j2+p01;	RADIX_09_DIF_FMA(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im, c,s,c2,s2,c3,s3,c4,s4);
		/*...and now do 9 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
	  #endif

	#endif	/* !USE_SSE2 */

	}	/* end for(j=_jstart; j < _jhi; j += 2) */

	jstart += nwt;
	jhi    += nwt;

	col += RADIX;
	co3 -= RADIX;
}	/* end for(int k=1; k <= khi; k++) */

#ifndef USE_ARM_V8_SIMD
  #undef OFF
#endif
