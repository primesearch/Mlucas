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

// Aug 2014: Need arbitrary-pointer-offsets to support I/O permutations needed by
// larger-radix DFTs of length 15 * 2^n:
#ifdef USE_LITERAL_BYTE_OFFSETS
	#error USE_LITERAL_BYTE_OFFSETS no longer supported for radix-15 DFTs!
#endif

// This main loop is same for un-and-multithreaded, so stick into a header file
// (can't use a macro because of the #if-enclosed stuff).

for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
{
	for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

/*...The radix-60 DIT pass is here:	*/

	/* The 15 radix-4 subblocks have several different SIMD flow modes: */
	#ifdef USE_AVX
	/* Outputs in AVX mode are temps 2*15*32 = 30*32 = 0x3c0 bytes apart: */
	// Reorder blocks to yield sequentially increasing a-array offsets:
		add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x3c0)
		add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x3c0)
		add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x3c0)
		add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x3c0)
		add1 = &a[j1+p16];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x3c0)
		add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x3c0)
		add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x3c0)
		add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x3c0)
		add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x3c0)
		add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x3c0)
		add0 = &a[j1+p40];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x3c0)
		add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x3c0)
		add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x3c0)
		add2 = &a[j1+p52];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x3c0)
		add1 = &a[j1+p56];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x3c0)

	#elif defined(USE_SSE2)

	/* Outputs in SSE2 mode are temps 2*15*16 = 30*16 = 0x1e0 bytes apart: */
	// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 01: */	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x1e0)
/* Block 03: */	add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x1e0)
/* Block 02: */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x1e0)
/* Block 06: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x1e0)
/* Block 05: */	add1 = &a[j1+p16];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x1e0)
/* Block 04: */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x1e0)
/* Block 08: */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x1e0)
/* Block 07: */	add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x1e0)
/* Block 09: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x1e0)
/* Block 10: */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x1e0)
/* Block 12: */	add0 = &a[j1+p40];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x1e0)
/* Block 11: */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x1e0)
/* Block 15: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x1e0)
/* Block 14: */	add2 = &a[j1+p52];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x1e0)
/* Block 13: */	add1 = &a[j1+p56];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x1e0)

	#endif	// AVX or SSE2?

	#ifdef USE_SSE2	// USE_AVX under same umbrella here

	  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.
		/* Radix-15 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 12 (24*16 bytes = 0x180) between successive outputs: */
	   #ifndef USE_LITERAL_BYTE_OFFSETS
		SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r00,  r01,  r02,  r03,  r04,  r05,  r06,  r07,  r08,  r09,  r0a,  r0b,  r0c,  r0d,  r0e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p00,s1p2e,s1p1d,s1p0c,s1p3b,s1p2a,s1p19,s1p08,s1p37,s1p26,s1p15,s1p04,s1p33,s1p22,s1p11);
		SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r10,  r11,  r12,  r13,  r14,  r15,  r16,  r17,  r18,  r19,  r1a,  r1b,  r1c,  r1d,  r1e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p10,s1p3e,s1p2d,s1p1c,s1p0b,s1p3a,s1p29,s1p18,s1p07,s1p36,s1p25,s1p14,s1p03,s1p32,s1p21);
		SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r20,  r21,  r22,  r23,  r24,  r25,  r26,  r27,  r28,  r29,  r2a,  r2b,  r2c,  r2d,  r2e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p20,s1p0e,s1p3d,s1p2c,s1p1b,s1p0a,s1p39,s1p28,s1p17,s1p06,s1p35,s1p24,s1p13,s1p02,s1p31);
		SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r30,  r31,  r32,  r33,  r34,  r35,  r36,  r37,  r38,  r39,  r3a,  r3b,  r3c,  r3d,  r3e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p30,s1p1e,s1p0d,s1p3c,s1p2b,s1p1a,s1p09,s1p38,s1p27,s1p16,s1p05,s1p34,s1p23,s1p12,s1p01);
	   #else	// Versions using base-address-plus-literal-byte-offsets:                                                                                       s1p**r complex-element offsets w.r.to s1p00 (all negative offsets get added to +120):     0    +88    +56    +24    -8     -40    -72    -104   -16    -48    -80    -112   -24    -56    -88 ... when offset < -(row-start offset below), reset by adding 120*16 = 0x780 bytes
		SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00, 0x000, 0x580, 0x380, 0x180, 0x700, 0x500, 0x300, 0x100, 0x680, 0x480, 0x280, 0x080, 0x600, 0x400, 0x200);	// end offset of each s1p-sequence = start-offset + 0x200
		SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r10,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00, 0x1e0, 0x760, 0x560, 0x360, 0x160, 0x6e0, 0x4e0, 0x2e0, 0x0e0, 0x660, 0x460, 0x260, 0x060, 0x5e0, 0x3e0);	// s1p10r = +30 complex = +0x1e0 bytes
		SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r20,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00, 0x3c0, 0x1c0, 0x740, 0x540, 0x340, 0x140, 0x6c0, 0x4c0, 0x2c0, 0x0c0, 0x640, 0x440, 0x240, 0x040, 0x5c0);	// s1p20r = +60 complex = +0x3c0 bytes
		SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r30,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00, 0x5a0, 0x3a0, 0x1a0, 0x720, 0x520, 0x320, 0x120, 0x6a0, 0x4a0, 0x2a0, 0x0a0, 0x620, 0x420, 0x220, 0x020);	// s1p30r = +90 complex = +0x5a0 bytes
	   #endif
	  #else
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,\
			r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e,\
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,\
			s1p00,s1p2e,s1p1d,s1p0c,s1p3b,s1p2a,s1p19,s1p08,s1p37,s1p26,s1p15,s1p04,s1p33,s1p22,s1p11,\
			r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e,\
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,\
			s1p10,s1p3e,s1p2d,s1p1c,s1p0b,s1p3a,s1p29,s1p18,s1p07,s1p36,s1p25,s1p14,s1p03,s1p32,s1p21);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,\
			r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e,\
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,\
			s1p20,s1p0e,s1p3d,s1p2c,s1p1b,s1p0a,s1p39,s1p28,s1p17,s1p06,s1p35,s1p24,s1p13,s1p02,s1p31,\
			r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e,\
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,\
			s1p30,s1p1e,s1p0d,s1p3c,s1p2b,s1p1a,s1p09,s1p38,s1p27,s1p16,s1p05,s1p34,s1p23,s1p12,s1p01);
	  #endif

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
		tm1 = s1p00; tmp = cy_r; itmp = bjmodn;
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

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
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

	/*	i =((uint32)(sw - bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
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
		l = 0; addr = cy_r; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1    ],a[j2    ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p01],a[j2+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p02],a[j2+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p03],a[j2+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,p08,...,p56
			cmplx_carry_norm_errcheck(a[jt    ],a[jp    ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p01],a[jp+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p02],a[jp+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p03],a[jp+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

	}
	else	/* MODULUS_TYPE_FERMAT */
	{

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

		tm0 = s1p00; tmp = base_negacyclic_root; l = 0xf00;
		tm1 = cy_r; // *cycle[] indices increment by +4 (mod ODD_RADIX) between macro calls
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic = 0; jc = 1; kc = 2; lc = 3;
		while(tm0 < s1p3e)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x1e0, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
			tm0 += 8; tm1++; tmp += 8; l -= 0xc0;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #else	/* HIACC = false: */

		tm0 = s1p00; tmp = base_negacyclic_root;	// tmp *not* incremented between macro calls in loacc version
		tm1 = cy_r;
		ic = 0; jc = 1; kc = 2; lc = 3;
		for(l = 0; l < RADIX>>2; l++) {	// RADIX/4 loop passes
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x1e0, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
			tm0 += 8; tm1++;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #endif	/* HIACC? */

	#elif defined(USE_SSE2)

		// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

	  #if (OS_BITS == 64)	// Run out of registers here in serial-build mode, so use (threaded or not?) to toggle carry-macro version selection here:

		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],jcycle[ic],icycle[ic+1],jcycle[ic+1] of the non-looped version with icycle[ic],jcycle[ic],icycle[jc],jcycle[jc]:
		ic = 0; jc = 1;
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while(tm1 < s1p3e) {
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];
			k2 = jcycle[ic];
			int k3 = icycle[jc];
			int k4 = jcycle[jc];
			SSE2_fermat_carry_norm_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2,k3,k4);
			tm1 += 4; tmp += 2;
			MOD_ADD32(ic, 2, ODD_RADIX, ic);
			MOD_ADD32(jc, 2, ODD_RADIX, jc);
		}

	  #else // Mar 2014: Worked around the out-of-regs compiler issues with the _X2 version of this macro (the
			// code in carry_gcc64.h has details), but keep non-X2 version in case hit out-of-regs again at some point

		ic = 0;	// ic = idx into [i|j]cycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		l = ODD_RADIX << 4;	// 32-bit version needs preshifted << 4 input value
		while(tm1 <= s1p3e) {
			//Sep 2014: Even with reduced-register version of the 32-bit Fermat-mod carry macro,
			// GCC runs out of registers on this one, without some playing-around-with-alternate code-sequences ...
			// Pulling the array-refs out of the carry-macro call like so solves the problem:
			k1 = icycle[ic];
			k2 = jcycle[ic];
			SSE2_fermat_carry_norm_errcheck(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2);
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
		for(l = 0; l < ODD_RADIX; l++) {
			icycle[l] += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
			icycle[l] += ( (-(int)((uint32)icycle[l] >> 31)) & nwt);
		}

	#endif	/* #ifdef USE_SSE2 */

	// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
	#ifdef USE_SSE2

		for(l = 0; l < ODD_RADIX; l++) {
			icycle[l] += wts_idx_inc2;		icycle[l] += ( (-(icycle[l] < 0)) & nwt16);
			jcycle[l] += wts_idx_inc2;		jcycle[l] += ( (-(jcycle[l] < 0)) & nwt16);
		  #ifdef USE_AVX
			kcycle[l] += wts_idx_inc2;		kcycle[l] += ( (-(kcycle[l] < 0)) & nwt16);
			lcycle[l] += wts_idx_inc2;		lcycle[l] += ( (-(lcycle[l] < 0)) & nwt16);
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
	/*
	RADIX_05_DFT(00,0c,19,26,33, r00,r01,r02,r03,r04);
	RADIX_05_DFT(15,22,2e,3b,08, r05,r06,r07,r08,r09);
	RADIX_05_DFT(2a,37,04,11,1d, r0a,r0b,r0c,r0d,r0e);

	RADIX_03_DFT(r00,r05,r0a,,00,01,02,);
	RADIX_03_DFT(r01,r06,r0b,,13,14,12,);
	RADIX_03_DFT(r02,r07,r0c,,09,10,11,);
	RADIX_03_DFT(r03,r08,r0d,,08,06,07,);
	RADIX_03_DFT(r04,r09,r0e,,04,05,03,);
	*/
	/*
	RADIX_05_DFT(30,3c,09,16,23, r10,r11,r12,r13,r14);
	RADIX_05_DFT(05,12,1e,2b,38, r15,r16,r17,r18,r19);
	RADIX_05_DFT(1a,27,34,01,0d, r1a,r1b,r1c,r1d,r1e);
	*/
	/*
	RADIX_05_DFT(20,2c,39,06,13, r20,r21,r22,r23,r24);
	RADIX_05_DFT(35,02,0e,1b,28, r25,r26,r27,r28,r29);
	RADIX_05_DFT(0a,17,24,31,3d, r2a,r2b,r2c,r2d,r2e);
	*/
	/*
	RADIX_05_DFT(10,1c,29,36,03, r30,r31,r32,r33,r34);
	RADIX_05_DFT(25,32,3e,0b,18, r35,r36,r37,r38,r39);
	RADIX_05_DFT(3a,07,14,21,2d, r3a,r3b,r3c,r3d,r3e);
	*/
	  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.
		/* NOTE that to permit us to re-use the s1p-data as both inputs and outputs of the radix-15 DIF DFT, must swap s1p[0123] -> s1p[0321] in output rows: */
	   #ifndef USE_LITERAL_BYTE_OFFSETS
		SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p00,s1p3b,s1p37,s1p33,s1p2e,s1p2a,s1p26,s1p22,s1p1d,s1p19,s1p15,s1p11,s1p0c,s1p08,s1p04, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e);
		SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p30,s1p2b,s1p27,s1p23,s1p1e,s1p1a,s1p16,s1p12,s1p0d,s1p09,s1p05,s1p01,s1p3c,s1p38,s1p34, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e);
		SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p20,s1p1b,s1p17,s1p13,s1p0e,s1p0a,s1p06,s1p02,s1p3d,s1p39,s1p35,s1p31,s1p2c,s1p28,s1p24, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e);
		SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p10,s1p0b,s1p07,s1p03,s1p3e,s1p3a,s1p36,s1p32,s1p2d,s1p29,s1p25,s1p21,s1p1c,s1p18,s1p14, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e);
	   #else	// Versions using base-address-plus-literal-byte-offsets:
		SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00, 0x000, 0x700, 0x680, 0x600, 0x580, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
		SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00, 0x5a0, 0x520, 0x4a0, 0x420, 0x3a0, 0x320, 0x2a0, 0x220, 0x1a0, 0x120, 0x0a0, 0x020, 0x720, 0x6a0, 0x620, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r10,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
		SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00, 0x3c0, 0x340, 0x2c0, 0x240, 0x1c0, 0x140, 0x0c0, 0x040, 0x740, 0x6c0, 0x640, 0x5c0, 0x540, 0x4c0, 0x440, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r20,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
		SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00, 0x1e0, 0x160, 0x0e0, 0x060, 0x760, 0x6e0, 0x660, 0x5e0, 0x560, 0x4e0, 0x460, 0x3e0, 0x360, 0x2e0, 0x260, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r30,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
	   #endif
	  #else
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1, s1p00,s1p3b,s1p37,s1p33,s1p2e,s1p2a,s1p26,s1p22,s1p1d,s1p19,s1p15,s1p11,s1p0c,s1p08,s1p04, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e, s1p30,s1p2b,s1p27,s1p23,s1p1e,s1p1a,s1p16,s1p12,s1p0d,s1p09,s1p05,s1p01,s1p3c,s1p38,s1p34, y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e, r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1, s1p20,s1p1b,s1p17,s1p13,s1p0e,s1p0a,s1p06,s1p02,s1p3d,s1p39,s1p35,s1p31,s1p2c,s1p28,s1p24, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e, s1p10,s1p0b,s1p07,s1p03,s1p3e,s1p3a,s1p36,s1p32,s1p2d,s1p29,s1p25,s1p21,s1p1c,s1p18,s1p14, y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e, r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e);
	  #endif

	#ifdef USE_AVX
	// Reorder blocks to yield sequentially increasing a-array offsets:
		add0 = &a[j1    ];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x3c0)
		add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x3c0)
		add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x3c0)
		add3 = &a[j1+p12];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x3c0)
		add1 = &a[j1+p16];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x3c0)
		add2 = &a[j1+p20];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x3c0)
		add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x3c0)
		add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x3c0)
		add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x3c0)
		add2 = &a[j1+p36];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x3c0)
		add0 = &a[j1+p40];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x3c0)
		add3 = &a[j1+p44];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x3c0)
		add0 = &a[j1+p48];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x3c0)
		add3 = &a[j1+p52];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x3c0)
		add1 = &a[j1+p56];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x3c0)

	#else

	// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 01: */	add0 = &a[j1    ];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x1e0)
/* Block 02: */	add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x1e0)
/* Block 03: */	add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x1e0)
/* Block 04: */	add1 = &a[j1+p56];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x1e0)
/* Block 05: */	add3 = &a[j1+p52];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x1e0)
/* Block 06: */	add0 = &a[j1+p48];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x1e0)
/* Block 07: */	add3 = &a[j1+p44];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x1e0)
/* Block 08: */	add0 = &a[j1+p40];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x1e0)
/* Block 09: */	add2 = &a[j1+p36];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x1e0)
/* Block 10: */	add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x1e0)
/* Block 11: */	add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x1e0)
/* Block 12: */	add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x1e0)
/* Block 13: */	add2 = &a[j1+p20];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x1e0)
/* Block 14: */	add1 = &a[j1+p16];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x1e0)
/* Block 15: */	add3 = &a[j1+p12];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x1e0)

	  #endif	// AVX or SSE2?

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

