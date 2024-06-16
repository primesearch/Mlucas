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

for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
{
	/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
	Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
	but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
	(But only ever need to explicitly do this in debug mode).
	*/
	for(j = jstart; j < jhi; j += stride)
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-60 DIT pass is here:	*/

	  #ifdef USE_ARM_V8_SIMD
		const uint32 OFF = 0x1e0;
	  #elif defined(USE_AVX512)
		#define OFF	0x1e0*4
	  #elif defined(USE_AVX)
		#define OFF	0x1e0*2
	  #else
		#define OFF	0x1e0
	  #endif

	#ifdef USE_SSE2

		/*...gather the needed data and do 15 radix-4 transforms...*/
	  #if COMPACT_OBJ
		#warning using COMPACT_OBJ code.
		/* Outputs in SSE2 modes are temps 2*15*16 = 22*16 = 0x160 bytes apart: */
		// Indices into above 4-elt table; each DFT-4 needs four 2-bit indices, thus gets 1 byte
		// Ex: 1st DFT-4 has add0-3 p01-multiple offsets 0,1,3,2; bit-reverse that to get p_idx[0] = 2310_4 = 0xb4:
		const uint8 p_id1[15] = {0xb4,0xe1,0x4e,0x1b,0xe1,0x4e,0x1b,0xe1,0xb4,0x1b,0x4e,0xb4,0xe1,0x4e,0xb4};
		const uint8 p_od1[15] = {0,2,1,5,4,3,7,6,8,9,11,10,14,13,12};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00; l < 15; l++, tmp+=2) {
			i3 = p_id1[l]; i0 = i3&3; i1 = (i3>>2)&3; i2 = (i3>>4)&3; i3 = (i3>>6);
			addr = &a[j1+poff[p_od1[l]]];
			add0 = addr+p0123[i0];	add1 = addr+p0123[i1];	add2 = addr+p0123[i2];	add3 = addr+p0123[i3];
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0,add1,add2,add3, tmp, OFF);
		}
		/*...and now do 4 radix-15 transforms. */
		// Radix-15 DFT outputs are (cyclic) with vec_dbl-pointer -= 32 (mod 120) between successive outputs
		const uint8 optr_off[RADIX] = {
			 0, 88,56,24,112,80,48,16,104,72,40,8,96,64,32,
			30,118,86,54,22,110,78,46,14,102,70,38,6,94,62,
			60,28,116,84,52,20,108,76,44,12,100,68,36,4,92,
			90,58,26,114,82,50,18,106,74,42,10,98,66,34, 2};
		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,*va9,*vaa,*vab,*vac,*vad,*vae,	// I-ptrs
		*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8,*vb9,*vba,*vbb,*vbc,*vbd,*vbe;	// O-ptrs
		vec_dbl
		*vc0,*vc1,*vc2,*vc3,*vc4,*vc5,*vc6,*vc7,*vc8,*vc9,*vca,*vcb,*vcc,*vcd,*vce,	// I-ptrs
		*vd0,*vd1,*vd2,*vd3,*vd4,*vd5,*vd6,*vd7,*vd8,*vd9,*vda,*vdb,*vdc,*vdd,*vde;	// O-ptrs
		for(l = 0, tmp = r00, ntmp = 0; l < 2; l++, ntmp += 30)
		{
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
			vab = tmp + 22;	vbb = s1p00 + optr_off[ntmp+11];
			vac = tmp + 24;	vbc = s1p00 + optr_off[ntmp+12];
			vad = tmp + 26;	vbd = s1p00 + optr_off[ntmp+13];
			vae = tmp + 28;	vbe = s1p00 + optr_off[ntmp+14];

			vc0 = tmp + 30;	vd0 = s1p00 + optr_off[ntmp+15];
			vc1 = tmp + 32;	vd1 = s1p00 + optr_off[ntmp+16];
			vc2 = tmp + 34;	vd2 = s1p00 + optr_off[ntmp+17];
			vc3 = tmp + 36;	vd3 = s1p00 + optr_off[ntmp+18];
			vc4 = tmp + 38;	vd4 = s1p00 + optr_off[ntmp+19];
			vc5 = tmp + 40;	vd5 = s1p00 + optr_off[ntmp+20];
			vc6 = tmp + 42;	vd6 = s1p00 + optr_off[ntmp+21];
			vc7 = tmp + 44;	vd7 = s1p00 + optr_off[ntmp+22];
			vc8 = tmp + 46;	vd8 = s1p00 + optr_off[ntmp+23];
			vc9 = tmp + 48;	vd9 = s1p00 + optr_off[ntmp+24];
			vca = tmp + 50;	vda = s1p00 + optr_off[ntmp+25];
			vcb = tmp + 52;	vdb = s1p00 + optr_off[ntmp+26];
			vcc = tmp + 54;	vdc = s1p00 + optr_off[ntmp+27];
			vcd = tmp + 56;	vdd = s1p00 + optr_off[ntmp+28];
			vce = tmp + 58;	vde = s1p00 + optr_off[ntmp+29];
			SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,two,
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,	/* inputs  1 */
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,	/* scratch 1 */
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,vbb,vbc,vbd,vbe,	/* outputs 1 */
				vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,	/* inputs  2 */
				y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,	/* scratch 2 */
				vd0,vd1,vd2,vd3,vd4,vd5,vd6,vd7,vd8,vd9,vda,vdb,vdc,vdd,vde		/* outputs 2 */
			);
			tmp += 60;
		}

	  #else

		/*...gather the needed data and do 15 radix-4 transforms...*/
	/* Outputs in SSE2 mode are temps 2*15*16 = 30*16 = 0x1e0 bytes apart: */
	/* Block 01: */	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, OFF)
	/* Block 02: */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, OFF)
	/* Block 03: */	add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, OFF)
	/* Block 04: */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, OFF)
	/* Block 05: */	add1 = &a[j1+p16];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, OFF)
	/* Block 06: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, OFF)
	/* Block 07: */	add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, OFF)
	/* Block 08: */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, OFF)
	/* Block 09: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, OFF)
	/* Block 10: */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, OFF)
	/* Block 11: */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, OFF)
	/* Block 12: */	add0 = &a[j1+p40];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, OFF)
	/* Block 13: */	add1 = &a[j1+p56];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, OFF)
	/* Block 14: */	add2 = &a[j1+p52];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, OFF)
	/* Block 15: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, OFF)

		/*...and now do 4 radix-15 transforms. */
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,two,
			r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1p00,s1p2e,s1p1d,s1p0c,s1p3b,s1p2a,s1p19,s1p08,s1p37,s1p26,s1p15,s1p04,s1p33,s1p22,s1p11,
			r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1p10,s1p3e,s1p2d,s1p1c,s1p0b,s1p3a,s1p29,s1p18,s1p07,s1p36,s1p25,s1p14,s1p03,s1p32,s1p21);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,two,
			r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1p20,s1p0e,s1p3d,s1p2c,s1p1b,s1p0a,s1p39,s1p28,s1p17,s1p06,s1p35,s1p24,s1p13,s1p02,s1p31,
			r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1p30,s1p1e,s1p0d,s1p3c,s1p2b,s1p1a,s1p09,s1p38,s1p27,s1p16,s1p05,s1p34,s1p23,s1p12,s1p01);

	  #endif	// COMPACT_OBJ?

	#else	// Scalar-double mode:

		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p44; jp = j2+p44;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p56; jp = j2+p56;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p52; jp = j2+p52;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);	tptr++;
		jt = j1+p48; jp = j2+p48;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,rt,it);

		/*...and now do 4 radix-15 transforms.
		The required output permutation is [in terms of decimal indices...we use base-15 indexing in the actual code below]

			00,44,28,12,56,40,24,08,52,36,20,04,48,32,16   00,44,28,12,56,40,24,08,52,36,20,04,48,32,16 + p0
			15,59,43,27,11,55,39,23,07,51,35,19,03,47,31 = 12,56,40,24,08,52,36,20,04,48,32,16,00,44,28 + p3
			30,14,58,42,26,10,54,38,22,06,50,34,18,02,46   28,12,56,40,24,08,52,36,20,04,48,32,16,00,44 + p2
			45,29,13,57,41,25,09,53,37,21,05,49,33,17,01   44,28,12,56,40,24,08,52,36,20,04,48,32,16,00 + p1
		*/
		tptr = t;	jt = j1    ; jp = j2    ;
		RADIX_15_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im
				,a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16])
		tptr += 15;	jt = j1+p03; jp = j2+p03;
		RADIX_15_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im
				,a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28])
		tptr += 15;	jt = j1+p02; jp = j2+p02;
		RADIX_15_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im
				,a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44])
		tptr += 15;	jt = j1+p01; jp = j2+p01;
		RADIX_15_DIT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im
				,a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ])

	#endif	// SIMD or not?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need 28 separate carries.	*/

/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
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

	  if(inc_arr[0]) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of the low
						// inc_arr element to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

		i = (!j);
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy_r; tm2 = cy_r+1; itmp = bjmodn; itm2 = bjmodn+4;
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
		// And do a final 4-fold pass for 56-59:
		add0 = a + j1 + pfetch_dist + poff[l+l];
	   #ifdef USE_AVX512
		// AVX-512 mode calls this macro twice, with Call 2ptr-offsets fiddled as described in comments to that version of the macro:
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,0x000, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
		add0 += p04;				// prefetch of a + [prefetch offset] + p4,5,6,7
		tmp  = (double *)tmp +  4;	// Call 2 will handle the .d4-7 doubles of our 4 input zmm register-sized vector data
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,0x800, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);	// Call 2 wts-data pointers += 0x400
	   #else
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,i,     sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
	   #endif

	  } else {	// HiACC:

		/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		tm1 = s1p00; tmp = cy_r; itmp = bjmodn;
		i = (!j);
		addr = &prp_mult;
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
		uint32 ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
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
			k1 = n-si[l  ];
			k2 = n-si[l+1];
			k3 = si[nwtml  ];
			k4 = si[nwtml-1];
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
			SSE2_cmplx_carry_fast_wtsinit(add1,add2,add3, add0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, k1,k2,k3,k4, sse_bw,sse_n)

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

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
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

	/*	i =((uint32)(sw - bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
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

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
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
		int k5,k6,k7;
		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

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

		/* The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
		which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls: */

		tm0 = s1p00; tmp = base_negacyclic_root; l = 0xf00;
		tm1 = cy_r; // *cycle[] indices increment by +4 (mod ODD_RADIX) between macro calls
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic_idx = 0; jc_idx = 1; kc_idx = 2; lc_idx = 3;
		while(tm0 < s1p3e)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic_idx];	k5 = jcycle[ic_idx];	k6 = kcycle[ic_idx];	k7 = lcycle[ic_idx];
			k2 = icycle[jc_idx];
			k3 = icycle[kc_idx];
			k4 = icycle[lc_idx];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x1e0, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p01,p02,p03, addr);
			tm0 += 8; tm1++; tmp += 8; l -= 0xc0;
			MOD_ADD32(ic_idx, 4, ODD_RADIX, ic_idx);
			MOD_ADD32(jc_idx, 4, ODD_RADIX, jc_idx);
			MOD_ADD32(kc_idx, 4, ODD_RADIX, kc_idx);
			MOD_ADD32(lc_idx, 4, ODD_RADIX, lc_idx);
		}

	  #else	// HIACC = false:

		// Oct 2014: Try getting most of the LOACC speedup with better accuracy by breaking the complex-roots-of-(-1)
		// chaining into 2 or more equal-sized subchains, each starting with 'fresh' (unchained) complex roots:
		tm0 = s1p00; tm1 = cy_r; // tm2 = cy_i;	*** replace with literal-byte-offset in macro call to save a reg
		ic_idx = 0; jc_idx = 1; kc_idx = 2; lc_idx = 3;
	  #ifdef USE_AVX512
		mc_idx = 4; nc_idx = 5; oc_idx = 6; pc_idx = 7;
	  #endif

		uint32 naccum = 0;	// Stores sum of [0-ntmp]th elements of inc_arr[]
		for(ntmp = 0; ntmp < (1 << nfold); ++ntmp)
		{
			// E.g.: nfold = 1 (==> 2^nfold = 2-subchains) means L takes its value
			// from (j) at start of 1st inner-loop exec, and from (j + n/2) at start of 2nd:
		//	l = (j + ntmp*(n>>nfold)) >> 1;	*** Only works if RADIX divisible by 2^(lg(RE_IM_STRIDE)+nfold)
			l = (j + naccum*NDIVR*RE_IM_STRIDE) >> 1;	naccum += inc_arr[ntmp];

		// Get the needed quartet (octet if AVX512) of Nth roots of -1: This is the same code as in the scalar
		// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
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

			// The other ptrs need to carry over from pvs loop, but this one needs resetting due to above 'multipliers refresh'
			tmp = base_negacyclic_root;	// tmp *not* incremented between macro calls in loacc version

		#ifdef USE_AVX512
			// will never hit this since have same assert in preprocessing code - just a placeholder/reminder:
			ASSERT(0, "radix60_ditN_cy_dif1: No AVX-512 support for Fermat-mod; Skipping this leading radix.");

		#else	// AVX / AVX2

			for(l = 0; l < inc_arr[ntmp]; l++) {
				k1 = icycle[ic_idx];
				k2 = icycle[jc_idx];	k5 = jcycle[ic_idx];
				k3 = icycle[kc_idx];	k6 = kcycle[ic_idx];
				k4 = icycle[lc_idx];	k7 = lcycle[ic_idx];
				// Each AVX carry macro call also processes 4 prefetches of main-array data
				tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
													/* (cy_i_cy_r) --vvvvv  vvvvvvvvvvvvvvvvv--[1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
				SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x1e0, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p01,p02,p03, addr);
				tm0 += 8; tm1++;
				MOD_ADD32(ic_idx, 4, ODD_RADIX, ic_idx);
				MOD_ADD32(jc_idx, 4, ODD_RADIX, jc_idx);
				MOD_ADD32(kc_idx, 4, ODD_RADIX, kc_idx);
				MOD_ADD32(lc_idx, 4, ODD_RADIX, lc_idx);
			}

		#endif
		}	// Outer (ntmp-indexed) loop

	  #endif	/* HIACC? */

	#elif defined(USE_SSE2)

		// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],jcycle[ic],icycle[ic+1],jcycle[ic+1] of the non-looped version with icycle[ic],jcycle[ic],icycle[jc],jcycle[jc]:
		ic_idx = 0; jc_idx = 1;
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while(tm1 < s1p3e) {
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic_idx];
			k2 = jcycle[ic_idx];
			k3 = icycle[jc_idx];
			k4 = jcycle[jc_idx];
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
		for(l = 0; l < ODD_RADIX; l++) {
			icycle[l] += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
			icycle[l] += ( (-(int)((uint32)icycle[l] >> 31)) & nwt);
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

	/*...The radix-60 DIF pass is here:	*/

	#ifdef USE_SSE2

	/* General indexing for radix-15 done as 3 radix-5 followed by 5 radix-3 is
		RADIX_15_DIF(00,01,02,03,04,05,06,07,08,09,0A,0B,0C,0D,0E) ==>

		RADIX_05_DFT(i0,iC,i9,i6,i3, t0,t1,t2,t3,t4)
		RADIX_05_DFT(iA,i7,i4,i1,iD, t5,t6,t7,t8,t9)
		RADIX_05_DFT(i5,i2,iE,iB,i8, tA,tB,tC,tD,tE)
			RADIX_03_DFT(t0,t5,tA, o0,o1,o2,)
			RADIX_03_DFT(t1,t6,tB, oD,oE,oB,)
			RADIX_03_DFT(t2,t7,tC, o9,oA,oB,)
			RADIX_03_DFT(t3,t8,tD, o8,o6,o7,)
			RADIX_03_DFT(t4,t9,tE, o4,o5,o3,)
		*/

	  #if COMPACT_OBJ

		/* Do 4 radix-15 transforms: */
		// Radix-15 DFT inputs are (cyclic) with vec_dbl-pointer -= 8 (mod 120) between successive outputs
		const uint8 iptr_off[RADIX] = {
			 0,112,104,96,88,80,72,64,56,48,40,32,24,16,8,
			90,82,74,66,58,50,42,34,26,18,10,2,114,106,98,
			60,52,44,36,28,20,12,4,116,108,100,92,84,76,68,
			30,22,14,6,118,110,102,94,86,78,70,62,54,46,38};
		for(l = 0, tmp = r00, ntmp = 0; l < 2; l++, ntmp += 30)
		{
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
			vab = tmp + 22;	vbb = s1p00 + iptr_off[ntmp+11];
			vac = tmp + 24;	vbc = s1p00 + iptr_off[ntmp+12];
			vad = tmp + 26;	vbd = s1p00 + iptr_off[ntmp+13];
			vae = tmp + 28;	vbe = s1p00 + iptr_off[ntmp+14];

			vc0 = tmp + 30;	vd0 = s1p00 + iptr_off[ntmp+15];
			vc1 = tmp + 32;	vd1 = s1p00 + iptr_off[ntmp+16];
			vc2 = tmp + 34;	vd2 = s1p00 + iptr_off[ntmp+17];
			vc3 = tmp + 36;	vd3 = s1p00 + iptr_off[ntmp+18];
			vc4 = tmp + 38;	vd4 = s1p00 + iptr_off[ntmp+19];
			vc5 = tmp + 40;	vd5 = s1p00 + iptr_off[ntmp+20];
			vc6 = tmp + 42;	vd6 = s1p00 + iptr_off[ntmp+21];
			vc7 = tmp + 44;	vd7 = s1p00 + iptr_off[ntmp+22];
			vc8 = tmp + 46;	vd8 = s1p00 + iptr_off[ntmp+23];
			vc9 = tmp + 48;	vd9 = s1p00 + iptr_off[ntmp+24];
			vca = tmp + 50;	vda = s1p00 + iptr_off[ntmp+25];
			vcb = tmp + 52;	vdb = s1p00 + iptr_off[ntmp+26];
			vcc = tmp + 54;	vdc = s1p00 + iptr_off[ntmp+27];
			vcd = tmp + 56;	vdd = s1p00 + iptr_off[ntmp+28];
			vce = tmp + 58;	vde = s1p00 + iptr_off[ntmp+29];
			SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,two,
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,vbb,vbc,vbd,vbe,	/* inputs  1 */
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,	/* scratch 1 */
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,	/* outputs 1 */
				vd0,vd1,vd2,vd3,vd4,vd5,vd6,vd7,vd8,vd9,vda,vdb,vdc,vdd,vde,	/* inputs  2 */
				y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,	/* scratch 2 */
				vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce		/* outputs 2 */
			);
			tmp += 60;
		}
		/*...and now do 15 radix-4 transforms...*/
		// Indices into above 4-elt table; each DFT-4 needs four 2-bit indices, thus gets 1 byte
		// Ex: 1st DFT-4 has add0-3 p01-multiple offsets 0,1,2,3; bit-reverse that to get p_idx[0] = 3210_4 = 0xe4:
		const uint8 p_id2[15] = {0xe4,0xb1,0x1e,0xb1,0x1e,0xe4,0x1e,0xe4,0x4b,0xe4,0x4b,0xb1,0x4b,0xb1,0x1e};
		const uint8 p_od2[15] = {0,2,1,14,13,12,11,10,9,8,7,6,5,4,3};	// Indices into poff[]-array of p04-multiples
		for(l = 0, tmp = r00; l < 15; l++, tmp+=2) {
			i3 = p_id2[l]; i0 = i3&3; i1 = (i3>>2)&3; i2 = (i3>>4)&3; i3 = (i3>>6);
		//	printf("Compact: a[]-offsets = %d + %d,%d,%d,%d\n",poff[p_od2[l]]/p01, p0123[i0]/p01,p0123[i1]/p01,p0123[i2]/p01,p0123[i3]/p01);
			addr = &a[j1+poff[p_od2[l]]];
			add0 = addr+p0123[i0];	add1 = addr+p0123[i1];	add2 = addr+p0123[i2];	add3 = addr+p0123[i3];
			SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0,add1,add2,add3, tmp, OFF);
		}

	  #else

		/* Do 4 radix-15 transforms: */
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,two, s1p00,s1p3b,s1p37,s1p33,s1p2e,s1p2a,s1p26,s1p22,s1p1d,s1p19,s1p15,s1p11,s1p0c,s1p08,s1p04, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e, s1p30,s1p2b,s1p27,s1p23,s1p1e,s1p1a,s1p16,s1p12,s1p0d,s1p09,s1p05,s1p01,s1p3c,s1p38,s1p34, y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e, r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,two, s1p20,s1p1b,s1p17,s1p13,s1p0e,s1p0a,s1p06,s1p02,s1p3d,s1p39,s1p35,s1p31,s1p2c,s1p28,s1p24, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e, s1p10,s1p0b,s1p07,s1p03,s1p3e,s1p3a,s1p36,s1p32,s1p2d,s1p29,s1p25,s1p21,s1p1c,s1p18,s1p14, y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e, r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e);

		/*...and now do 15 radix-4 transforms...*/
	/* Block 01: */	add0 = &a[j1    ];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, OFF)
	/* Block 02: */	add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, OFF)
	/* Block 03: */	add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, OFF)
	/* Block 04: */	add1 = &a[j1+p56];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, OFF)
	/* Block 05: */	add3 = &a[j1+p52];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, OFF)
	/* Block 06: */	add0 = &a[j1+p48];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, OFF)
	/* Block 07: */	add3 = &a[j1+p44];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, OFF)
	/* Block 08: */	add0 = &a[j1+p40];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, OFF)
	/* Block 09: */	add2 = &a[j1+p36];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, OFF)
	/* Block 10: */	add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, OFF)
	/* Block 11: */	add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, OFF)
	/* Block 12: */	add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, OFF)
	/* Block 13: */	add2 = &a[j1+p20];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, OFF)
	/* Block 14: */	add1 = &a[j1+p16];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, OFF)
	/* Block 15: */	add3 = &a[j1+p12];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, OFF)

	  #endif	// COMPACT_OBJ?

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF
	  #endif

	#else	/* !USE_SSE2 */

		tptr = t;	jt = j1    ; jp = j2    ;	RADIX_15_DIF(a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],
															tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im);
		tptr += 15;	jt = j1+p01; jp = j2+p01;	RADIX_15_DIF(a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],
															tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im);
		tptr += 15;	jt = j1+p02; jp = j2+p02;	RADIX_15_DIF(a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],
															tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im);
		tptr += 15;	jt = j1+p03; jp = j2+p03;	RADIX_15_DIF(a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],
															tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im);
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p56; jp = j2+p56;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p52; jp = j2+p52;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p48; jp = j2+p48;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p44; jp = j2+p44;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

	#endif	// SSE2 or AVX?

	}	/* end for(j=_jstart; j < _jhi; j += 2) */

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(k=1; k <= khi; k++) */

