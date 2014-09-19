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

#ifndef USE_COMPACT_OBJ_CODE
	#error USE_COMPACT_OBJ_CODE expected, but not found!
#endif

// This main loop is same for un-and-multithreaded, so stick into a header file
// (can't use a macro because of the #if-enclosed stuff).

for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
{
	for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

	/*...The radix-240 DIT pass is here:	*/

	#ifdef USE_SSE2

	  #if USE_COMPACT_OBJ_CODE

		tmp = r00;
		for(l = 0; l < 15; l++) {
			i64 = dit_perm16[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			addr = &a[j1] + p_in_hi[l];	// p_in_hi[] = p0,p20,...,p30
			add0=addr+k0; add1=addr+k1; add2=addr+k2; add3=addr+k3; add4=addr+k4; add5=addr+k5; add6=addr+k6; add7=addr+k7; add8=addr+k8; add9=addr+k9; adda=addr+ka; addb=addr+kb; addc=addr+kc; addd=addr+kd; adde=addr+ke; addf=addr+kf;
			SSE2_RADIX16_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, isrt2, tmp,OFF1,OFF2,OFF3,OFF4)
			tmp += 2;
		}

	  #else
		#error need to restore original code here!
	  #endif

		/*...and now do 16 radix-15 transforms. See radix240_dit_pass1() for details on the oindex patterning.
		In 64-bit mode we use the 16-register/64-bit-mode doubled-radix-15-DFT macros as detailed in the radix60 carry routine.
		*/
	  #if (OS_BITS == 32)

		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,*va9,*vaa,*vab,*vac,*vad,*vae,
		*vc0,*vc1,*vc2,*vc3,*vc4,*vc5,*vc6,*vc7,*vc8,*vc9,*vca,*vcb,*vcc,*vcd,*vce;

		tmp = r00;
		tm2 = s1p00;
		for(l = 0; l < 16; l++) {
			// NB: All offsets end up being shifted "one too many" bits leftward since
			// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:

			// Low bits of the output-ptr offsets: 2*[0,f,e,...,2,1]
			jt = ((16 - l)&0xf) << 1;	// 2*[0,f,e,...,2,1]
			k0 = (l-1)&(-(l>0)); 		// 0,0,1,2...e
												// Now get the resulting shifted offsets, adding the low bits to each:
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = (k0 << 5) + jt;
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = (k1 << 5) + jt;
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = (k2 << 5) + jt;
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = (k3 << 5) + jt;
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = (k4 << 5) + jt;
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = (k5 << 5) + jt;
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = (k6 << 5) + jt;
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = (k7 << 5) + jt;
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = (k8 << 5) + jt;
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = (k9 << 5) + jt;
			kb = ka-1; kb += (-(kb < 0))&15;		ka = (ka << 5) + jt;
			kc = kb-1; kc += (-(kc < 0))&15;		kb = (kb << 5) + jt;
			kd = kc-1; kd += (-(kd < 0))&15;		kc = (kc << 5) + jt;
			ke = kd-1; ke += (-(ke < 0))&15;		kd = (kd << 5) + jt;
													ke = (ke << 5) + jt;
			// Output ptrs:
			vc0 = tm2+k0;
			vc1 = tm2+k1;
			vc2 = tm2+k2;
			vc3 = tm2+k3;
			vc4 = tm2+k4;
			vc5 = tm2+k5;
			vc6 = tm2+k6;
			vc7 = tm2+k7;
			vc8 = tm2+k8;
			vc9 = tm2+k9;
			vca = tm2+ka;
			vcb = tm2+kb;
			vcc = tm2+kc;
			vcd = tm2+kd;
			vce = tm2+ke;
			
			// Input ptrs:		
			va0 = tmp     ;
			va1 = tmp+0x02;
			va2 = tmp+0x04;
			va3 = tmp+0x06;
			va4 = tmp+0x08;
			va5 = tmp+0x0a;
			va6 = tmp+0x0c;
			va7 = tmp+0x0e;
			va8 = tmp+0x10;
			va9 = tmp+0x12;
			vaa = tmp+0x14;
			vab = tmp+0x16;
			vac = tmp+0x18;
			vad = tmp+0x1a;
			vae = tmp+0x1c;

			SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1,
/* r00 +... */	va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
/* s1p00 +... */	vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce
			);
			tmp += 0x1e;	// advance ptr by 15 vec_cmplx [= 30 vec_dbl] elements
		}

	  #else	// (OS_BITS == 64):

	   #if USE_COMPACT_OBJ_CODE
		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,*va9,*vaa,*vab,*vac,*vad,*vae,
		*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8,*vb9,*vba,*vbb,*vbc,*vbd,*vbe,
		*vc0,*vc1,*vc2,*vc3,*vc4,*vc5,*vc6,*vc7,*vc8,*vc9,*vca,*vcb,*vcc,*vcd,*vce,
		*vd0,*vd1,*vd2,*vd3,*vd4,*vd5,*vd6,*vd7,*vd8,*vd9,*vda,*vdb,*vdc,*vdd,*vde;

		tmp = r00;
		tm2 = s1p00;
		for(l = 0; l < 16; l += 2) {
			// NB: All offsets end up being shifted "one too many" bits leftward since
			// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:

			// Low bits of the output-ptr offsets: 2*[0,f,e,...,2,1]
			jt = ((16 - l)&0xf) << 1;
			k0 = (l-1)&(-(l>0)); /* even-order terms of 0,0,1,2...e, i.e. 0,1,3,... */
												// Now get the resulting shifted offsets, adding the low bits to each:
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = (k0 << 5) + jt;
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = (k1 << 5) + jt;
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = (k2 << 5) + jt;
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = (k3 << 5) + jt;
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = (k4 << 5) + jt;
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = (k5 << 5) + jt;
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = (k6 << 5) + jt;
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = (k7 << 5) + jt;
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = (k8 << 5) + jt;
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = (k9 << 5) + jt;
			kb = ka-1; kb += (-(kb < 0))&15;		ka = (ka << 5) + jt;
			kc = kb-1; kc += (-(kc < 0))&15;		kb = (kb << 5) + jt;
			kd = kc-1; kd += (-(kd < 0))&15;		kc = (kc << 5) + jt;
			ke = kd-1; ke += (-(ke < 0))&15;		kd = (kd << 5) + jt;
													ke = (ke << 5) + jt;
			// 1st set of Output ptrs:
			vc0 = tm2+k0;
			vc1 = tm2+k1;
			vc2 = tm2+k2;
			vc3 = tm2+k3;
			vc4 = tm2+k4;
			vc5 = tm2+k5;
			vc6 = tm2+k6;
			vc7 = tm2+k7;
			vc8 = tm2+k8;
			vc9 = tm2+k9;
			vca = tm2+ka;
			vcb = tm2+kb;
			vcc = tm2+kc;
			vcd = tm2+kd;
			vce = tm2+ke;
			
			// Now emulate the elided odd-index loop pass (replace l with l+1 in the above sequence) to get next set of ks:
			jt = (15 - l) << 1;	// Subtract can't yield < 0 here, so no need for mask of result
			k0 = l; /* odd-order terms of 0,0,1,2...e, i.e. 0,2,4,... */
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = (k0 << 5) + jt;
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = (k1 << 5) + jt;
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = (k2 << 5) + jt;
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = (k3 << 5) + jt;
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = (k4 << 5) + jt;
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = (k5 << 5) + jt;
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = (k6 << 5) + jt;
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = (k7 << 5) + jt;
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = (k8 << 5) + jt;
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = (k9 << 5) + jt;
			kb = ka-1; kb += (-(kb < 0))&15;		ka = (ka << 5) + jt;
			kc = kb-1; kc += (-(kc < 0))&15;		kb = (kb << 5) + jt;
			kd = kc-1; kd += (-(kd < 0))&15;		kc = (kc << 5) + jt;
			ke = kd-1; ke += (-(ke < 0))&15;		kd = (kd << 5) + jt;
													ke = (ke << 5) + jt;
			// 2nd set of Output ptrs:
			vd0 = tm2+k0;
			vd1 = tm2+k1;
			vd2 = tm2+k2;
			vd3 = tm2+k3;
			vd4 = tm2+k4;
			vd5 = tm2+k5;
			vd6 = tm2+k6;
			vd7 = tm2+k7;
			vd8 = tm2+k8;
			vd9 = tm2+k9;
			vda = tm2+ka;
			vdb = tm2+kb;
			vdc = tm2+kc;
			vdd = tm2+kd;
			vde = tm2+ke;

			// 2 sets of Input ptrs:		
							tm1 = tmp + 0x1e;
			va0 = tmp     ;	vb0 = tm1     ;
			va1 = tmp+0x02;	vb1 = tm1+0x02;
			va2 = tmp+0x04;	vb2 = tm1+0x04;
			va3 = tmp+0x06;	vb3 = tm1+0x06;
			va4 = tmp+0x08;	vb4 = tm1+0x08;
			va5 = tmp+0x0a;	vb5 = tm1+0x0a;
			va6 = tmp+0x0c;	vb6 = tm1+0x0c;
			va7 = tmp+0x0e;	vb7 = tm1+0x0e;
			va8 = tmp+0x10;	vb8 = tm1+0x10;
			va9 = tmp+0x12;	vb9 = tm1+0x12;
			vaa = tmp+0x14;	vba = tm1+0x14;
			vab = tmp+0x16;	vbb = tm1+0x16;
			vac = tmp+0x18;	vbc = tm1+0x18;
			vad = tmp+0x1a;	vbd = tm1+0x1a;
			vae = tmp+0x1c;	vbe = tm1+0x1c;

			SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
/* r00 +... */	va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
/* s1p00 +... */	vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,
/* r0f +... */	vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,vbb,vbc,vbd,vbe,
				y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
/* s1p00+jt +... */	vd0,vd1,vd2,vd3,vd4,vd5,vd6,vd7,vd8,vd9,vda,vdb,vdc,vdd,vde
			);
			tmp += 0x3c;	// advance ptr by 30 vec_cmplx [= 60 vec_dbl] elements
		}

	   #else

		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
			r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1p00,s1pe0,s1pd0,s1pc0,s1pb0,s1pa0,s1p90,s1p80,s1p70,s1p60,s1p50,s1p40,s1p30,s1p20,s1p10,
			r0f,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1p0f,s1pef,s1pdf,s1pcf,s1pbf,s1paf,s1p9f,s1p8f,s1p7f,s1p6f,s1p5f,s1p4f,s1p3f,s1p2f,s1p1f);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
			r1e,r1f,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1p1e,s1p0e,s1pee,s1pde,s1pce,s1pbe,s1pae,s1p9e,s1p8e,s1p7e,s1p6e,s1p5e,s1p4e,s1p3e,s1p2e,
			r2d,r2e,r2f,r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1p2d,s1p1d,s1p0d,s1ped,s1pdd,s1pcd,s1pbd,s1pad,s1p9d,s1p8d,s1p7d,s1p6d,s1p5d,s1p4d,s1p3d);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
			r3c,r3d,r3e,r3f,r40,r41,r42,r43,r44,r45,r46,r47,r48,r49,r4a,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1p3c,s1p2c,s1p1c,s1p0c,s1pec,s1pdc,s1pcc,s1pbc,s1pac,s1p9c,s1p8c,s1p7c,s1p6c,s1p5c,s1p4c,
			r4b,r4c,r4d,r4e,r4f,r50,r51,r52,r53,r54,r55,r56,r57,r58,r59,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1p4b,s1p3b,s1p2b,s1p1b,s1p0b,s1peb,s1pdb,s1pcb,s1pbb,s1pab,s1p9b,s1p8b,s1p7b,s1p6b,s1p5b);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
			r5a,r5b,r5c,r5d,r5e,r5f,r60,r61,r62,r63,r64,r65,r66,r67,r68,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1p5a,s1p4a,s1p3a,s1p2a,s1p1a,s1p0a,s1pea,s1pda,s1pca,s1pba,s1paa,s1p9a,s1p8a,s1p7a,s1p6a,
			r69,r6a,r6b,r6c,r6d,r6e,r6f,r70,r71,r72,r73,r74,r75,r76,r77,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1p69,s1p59,s1p49,s1p39,s1p29,s1p19,s1p09,s1pe9,s1pd9,s1pc9,s1pb9,s1pa9,s1p99,s1p89,s1p79);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
			r78,r79,r7a,r7b,r7c,r7d,r7e,r7f,r80,r81,r82,r83,r84,r85,r86,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1p78,s1p68,s1p58,s1p48,s1p38,s1p28,s1p18,s1p08,s1pe8,s1pd8,s1pc8,s1pb8,s1pa8,s1p98,s1p88,
			r87,r88,r89,r8a,r8b,r8c,r8d,r8e,r8f,r90,r91,r92,r93,r94,r95,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1p87,s1p77,s1p67,s1p57,s1p47,s1p37,s1p27,s1p17,s1p07,s1pe7,s1pd7,s1pc7,s1pb7,s1pa7,s1p97);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
			r96,r97,r98,r99,r9a,r9b,r9c,r9d,r9e,r9f,ra0,ra1,ra2,ra3,ra4,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1p96,s1p86,s1p76,s1p66,s1p56,s1p46,s1p36,s1p26,s1p16,s1p06,s1pe6,s1pd6,s1pc6,s1pb6,s1pa6,
			ra5,ra6,ra7,ra8,ra9,raa,rab,rac,rad,rae,raf,rb0,rb1,rb2,rb3,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1pa5,s1p95,s1p85,s1p75,s1p65,s1p55,s1p45,s1p35,s1p25,s1p15,s1p05,s1pe5,s1pd5,s1pc5,s1pb5);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
			rb4,rb5,rb6,rb7,rb8,rb9,rba,rbb,rbc,rbd,rbe,rbf,rc0,rc1,rc2,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1pb4,s1pa4,s1p94,s1p84,s1p74,s1p64,s1p54,s1p44,s1p34,s1p24,s1p14,s1p04,s1pe4,s1pd4,s1pc4,
			rc3,rc4,rc5,rc6,rc7,rc8,rc9,rca,rcb,rcc,rcd,rce,rcf,rd0,rd1,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1pc3,s1pb3,s1pa3,s1p93,s1p83,s1p73,s1p63,s1p53,s1p43,s1p33,s1p23,s1p13,s1p03,s1pe3,s1pd3);
		SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
			rd2,rd3,rd4,rd5,rd6,rd7,rd8,rd9,rda,rdb,rdc,rdd,rde,rdf,re0,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			s1pd2,s1pc2,s1pb2,s1pa2,s1p92,s1p82,s1p72,s1p62,s1p52,s1p42,s1p32,s1p22,s1p12,s1p02,s1pe2,
			re1,re2,re3,re4,re5,re6,re7,re8,re9,rea,reb,rec,red,ree,ref,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			s1pe1,s1pd1,s1pc1,s1pb1,s1pa1,s1p91,s1p81,s1p71,s1p61,s1p51,s1p41,s1p31,s1p21,s1p11,s1p01);
	   #endif

	 #endif	// (OS_BITS == 32)

	#else	/* !USE_SSE2 */

	// Gather the needed data (240 64-bit complex) and do 15 radix-16 transforms.
	// The loop-based compact-object-code version of the radix-240 DIT below saves ~250kB off the *.o file
	// [similar for the corresponding DIF]. Compact-obj also saves similarly on build and - best of all - runtimes:
		tptr = t;
		for(ntmp = 0; ntmp < 15; ntmp++) {
			i64 = dit_perm16[ntmp];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			jt = j1 + p_in_hi[ntmp]; jp = j2 + p_in_hi[ntmp];	// p_in_hi[] = p0,p20,...,p30
			RADIX_16_DIT(a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
				c16,s16);	tptr++;
		}
		//...and now do 16 radix-15 transforms; cf. radix240_dit_pass1 for indexing details:
		tptr = t;
		for(ntmp = 0; ntmp < 16; ntmp++) {
			k0 = plo[(16 - ntmp)&0xf];	// p0,f,e,...,2,1
			jt = j1 + k0; jp = j2 + k0;
			k0 = (ntmp-1)&(-(ntmp>0)); /* 0,0,1,2...e */	// Now get the resulting p* offsets:
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = p_out_hi[k0];
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = p_out_hi[k1];
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = p_out_hi[k2];
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = p_out_hi[k3];
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = p_out_hi[k4];
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = p_out_hi[k5];
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = p_out_hi[k6];
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = p_out_hi[k7];
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = p_out_hi[k8];
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = p_out_hi[k9];
			kb = ka-1; kb += (-(kb < 0))&15;		ka = p_out_hi[ka];
			kc = kb-1; kc += (-(kc < 0))&15;		kb = p_out_hi[kb];
			kd = kc-1; kd += (-(kd < 0))&15;		kc = p_out_hi[kc];
			ke = kd-1; ke += (-(ke < 0))&15;		kd = p_out_hi[kd];
													ke = p_out_hi[ke];
			RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke]
			);	tptr += 0xf;
		}

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

	/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

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

		/*...set0 is slightly different from others; divide work into blocks of RADIX/4 macro calls, 1st set of which gets pulled out of loop: */		
		l = 0; addr = cy_r; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1   ],a[j2   ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p1],a[j2+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p2],a[j2+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p3],a[j2+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];
			cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
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

		tm0 = s1p00; tmp = base_negacyclic_root; l = 0x3c00;
		tm1 = cy_r; // tm2 = cy_i;	*** replace with literal-byte-offset in macro call to save a reg
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic = 0; jc = 1; kc = 2; lc = 3;
		while(tm0 < s1pef)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x780, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
			tm0 += 8; tm1++; tmp += 8; l -= 0xc0;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #else	/* HIACC = false: */

		tm0 = s1p00; tmp = base_negacyclic_root;	// tmp *not* incremented between macro calls in loacc version
		tm1 = cy_r; // tm2 = cy_i;	*** replace with literal-byte-offset in macro call to save a reg
		ic = 0; jc = 1; kc = 2; lc = 3;
		for(l = 0; l < RADIX>>2; l++) {	// RADIX/4 loop passes
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x780, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
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
		while(tm1 < s1pef) {
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
		while(tm1 <= s1pef) {
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
			jt = j1 + poff[m]; jp = j2 + poff[m];
			fermat_carry_norm_errcheckB(a[jt   ],a[jp   ],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p1],a[jp+p1],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p2],a[jp+p2],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p3],a[jp+p3],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}

		icycle[ 0] += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
		icycle[ 1] += wts_idx_incr;
		icycle[ 2] += wts_idx_incr;
		icycle[ 3] += wts_idx_incr;
		icycle[ 4] += wts_idx_incr;
		icycle[ 5] += wts_idx_incr;
		icycle[ 6] += wts_idx_incr;
		icycle[ 7] += wts_idx_incr;
		icycle[ 8] += wts_idx_incr;
		icycle[ 9] += wts_idx_incr;
		icycle[10] += wts_idx_incr;
		icycle[11] += wts_idx_incr;
		icycle[12] += wts_idx_incr;
		icycle[13] += wts_idx_incr;
		icycle[14] += wts_idx_incr;
		icycle[ 0] += ( (-(int)((uint32)icycle[ 0] >> 31)) & nwt);
		icycle[ 1] += ( (-(int)((uint32)icycle[ 1] >> 31)) & nwt);
		icycle[ 2] += ( (-(int)((uint32)icycle[ 2] >> 31)) & nwt);
		icycle[ 3] += ( (-(int)((uint32)icycle[ 3] >> 31)) & nwt);
		icycle[ 4] += ( (-(int)((uint32)icycle[ 4] >> 31)) & nwt);
		icycle[ 5] += ( (-(int)((uint32)icycle[ 5] >> 31)) & nwt);
		icycle[ 6] += ( (-(int)((uint32)icycle[ 6] >> 31)) & nwt);
		icycle[ 7] += ( (-(int)((uint32)icycle[ 7] >> 31)) & nwt);
		icycle[ 8] += ( (-(int)((uint32)icycle[ 8] >> 31)) & nwt);
		icycle[ 9] += ( (-(int)((uint32)icycle[ 9] >> 31)) & nwt);
		icycle[10] += ( (-(int)((uint32)icycle[10] >> 31)) & nwt);
		icycle[11] += ( (-(int)((uint32)icycle[11] >> 31)) & nwt);
		icycle[12] += ( (-(int)((uint32)icycle[12] >> 31)) & nwt);
		icycle[13] += ( (-(int)((uint32)icycle[13] >> 31)) & nwt);
		icycle[14] += ( (-(int)((uint32)icycle[14] >> 31)) & nwt);

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
		icycle[ 7] += wts_idx_inc2;		icycle[ 7] += ( (-(icycle[ 7] < 0)) & nwt16);
		icycle[ 8] += wts_idx_inc2;		icycle[ 8] += ( (-(icycle[ 8] < 0)) & nwt16);
		icycle[ 9] += wts_idx_inc2;		icycle[ 9] += ( (-(icycle[ 9] < 0)) & nwt16);
		icycle[10] += wts_idx_inc2;		icycle[10] += ( (-(icycle[10] < 0)) & nwt16);
		icycle[11] += wts_idx_inc2;		icycle[11] += ( (-(icycle[11] < 0)) & nwt16);
		icycle[12] += wts_idx_inc2;		icycle[12] += ( (-(icycle[12] < 0)) & nwt16);
		icycle[13] += wts_idx_inc2;		icycle[13] += ( (-(icycle[13] < 0)) & nwt16);
		icycle[14] += wts_idx_inc2;		icycle[14] += ( (-(icycle[14] < 0)) & nwt16);

		jcycle[ 0] += wts_idx_inc2;		jcycle[ 0] += ( (-(jcycle[ 0] < 0)) & nwt16);
		jcycle[ 1] += wts_idx_inc2;		jcycle[ 1] += ( (-(jcycle[ 1] < 0)) & nwt16);
		jcycle[ 2] += wts_idx_inc2;		jcycle[ 2] += ( (-(jcycle[ 2] < 0)) & nwt16);
		jcycle[ 3] += wts_idx_inc2;		jcycle[ 3] += ( (-(jcycle[ 3] < 0)) & nwt16);
		jcycle[ 4] += wts_idx_inc2;		jcycle[ 4] += ( (-(jcycle[ 4] < 0)) & nwt16);
		jcycle[ 5] += wts_idx_inc2;		jcycle[ 5] += ( (-(jcycle[ 5] < 0)) & nwt16);
		jcycle[ 6] += wts_idx_inc2;		jcycle[ 6] += ( (-(jcycle[ 6] < 0)) & nwt16);
		jcycle[ 7] += wts_idx_inc2;		jcycle[ 7] += ( (-(jcycle[ 7] < 0)) & nwt16);
		jcycle[ 8] += wts_idx_inc2;		jcycle[ 8] += ( (-(jcycle[ 8] < 0)) & nwt16);
		jcycle[ 9] += wts_idx_inc2;		jcycle[ 9] += ( (-(jcycle[ 9] < 0)) & nwt16);
		jcycle[10] += wts_idx_inc2;		jcycle[10] += ( (-(jcycle[10] < 0)) & nwt16);
		jcycle[11] += wts_idx_inc2;		jcycle[11] += ( (-(jcycle[11] < 0)) & nwt16);
		jcycle[12] += wts_idx_inc2;		jcycle[12] += ( (-(jcycle[12] < 0)) & nwt16);
		jcycle[13] += wts_idx_inc2;		jcycle[13] += ( (-(jcycle[13] < 0)) & nwt16);
		jcycle[14] += wts_idx_inc2;		jcycle[14] += ( (-(jcycle[14] < 0)) & nwt16);

	  #ifdef USE_AVX
		kcycle[ 0] += wts_idx_inc2;		kcycle[ 0] += ( (-(kcycle[ 0] < 0)) & nwt16);
		kcycle[ 1] += wts_idx_inc2;		kcycle[ 1] += ( (-(kcycle[ 1] < 0)) & nwt16);
		kcycle[ 2] += wts_idx_inc2;		kcycle[ 2] += ( (-(kcycle[ 2] < 0)) & nwt16);
		kcycle[ 3] += wts_idx_inc2;		kcycle[ 3] += ( (-(kcycle[ 3] < 0)) & nwt16);
		kcycle[ 4] += wts_idx_inc2;		kcycle[ 4] += ( (-(kcycle[ 4] < 0)) & nwt16);
		kcycle[ 5] += wts_idx_inc2;		kcycle[ 5] += ( (-(kcycle[ 5] < 0)) & nwt16);
		kcycle[ 6] += wts_idx_inc2;		kcycle[ 6] += ( (-(kcycle[ 6] < 0)) & nwt16);
		kcycle[ 7] += wts_idx_inc2;		kcycle[ 7] += ( (-(kcycle[ 7] < 0)) & nwt16);
		kcycle[ 8] += wts_idx_inc2;		kcycle[ 8] += ( (-(kcycle[ 8] < 0)) & nwt16);
		kcycle[ 9] += wts_idx_inc2;		kcycle[ 9] += ( (-(kcycle[ 9] < 0)) & nwt16);
		kcycle[10] += wts_idx_inc2;		kcycle[10] += ( (-(kcycle[10] < 0)) & nwt16);
		kcycle[11] += wts_idx_inc2;		kcycle[11] += ( (-(kcycle[11] < 0)) & nwt16);
		kcycle[12] += wts_idx_inc2;		kcycle[12] += ( (-(kcycle[12] < 0)) & nwt16);
		kcycle[13] += wts_idx_inc2;		kcycle[13] += ( (-(kcycle[13] < 0)) & nwt16);
		kcycle[14] += wts_idx_inc2;		kcycle[14] += ( (-(kcycle[14] < 0)) & nwt16);

		lcycle[ 0] += wts_idx_inc2;		lcycle[ 0] += ( (-(lcycle[ 0] < 0)) & nwt16);
		lcycle[ 1] += wts_idx_inc2;		lcycle[ 1] += ( (-(lcycle[ 1] < 0)) & nwt16);
		lcycle[ 2] += wts_idx_inc2;		lcycle[ 2] += ( (-(lcycle[ 2] < 0)) & nwt16);
		lcycle[ 3] += wts_idx_inc2;		lcycle[ 3] += ( (-(lcycle[ 3] < 0)) & nwt16);
		lcycle[ 4] += wts_idx_inc2;		lcycle[ 4] += ( (-(lcycle[ 4] < 0)) & nwt16);
		lcycle[ 5] += wts_idx_inc2;		lcycle[ 5] += ( (-(lcycle[ 5] < 0)) & nwt16);
		lcycle[ 6] += wts_idx_inc2;		lcycle[ 6] += ( (-(lcycle[ 6] < 0)) & nwt16);
		lcycle[ 7] += wts_idx_inc2;		lcycle[ 7] += ( (-(lcycle[ 7] < 0)) & nwt16);
		lcycle[ 8] += wts_idx_inc2;		lcycle[ 8] += ( (-(lcycle[ 8] < 0)) & nwt16);
		lcycle[ 9] += wts_idx_inc2;		lcycle[ 9] += ( (-(lcycle[ 9] < 0)) & nwt16);
		lcycle[10] += wts_idx_inc2;		lcycle[10] += ( (-(lcycle[10] < 0)) & nwt16);
		lcycle[11] += wts_idx_inc2;		lcycle[11] += ( (-(lcycle[11] < 0)) & nwt16);
		lcycle[12] += wts_idx_inc2;		lcycle[12] += ( (-(lcycle[12] < 0)) & nwt16);
		lcycle[13] += wts_idx_inc2;		lcycle[13] += ( (-(lcycle[13] < 0)) & nwt16);
		lcycle[14] += wts_idx_inc2;		lcycle[14] += ( (-(lcycle[14] < 0)) & nwt16);
	  #endif
	#endif

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-240 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #if (OS_BITS == 32)
	
		tmp = r00;
		tm2 = s1p00;
		for(l = 0; l < 16; l++) {
			// NB: All offsets end up being shifted "one too many" bits leftward since
			// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:

			// Low bits of the input-ptr offsets: 2*[0,1,2,...,e,f]
			jt = l << 1;
			k0 = (15-l)&(-(l>0)); // 0,e,d,...,0
												// Now get the resulting shifted offsets, adding the low bits to each:
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = (k0 << 5) + jt;
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = (k1 << 5) + jt;
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = (k2 << 5) + jt;
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = (k3 << 5) + jt;
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = (k4 << 5) + jt;
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = (k5 << 5) + jt;
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = (k6 << 5) + jt;
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = (k7 << 5) + jt;
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = (k8 << 5) + jt;
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = (k9 << 5) + jt;
			kb = ka-1; kb += (-(kb < 0))&15;		ka = (ka << 5) + jt;
			kc = kb-1; kc += (-(kc < 0))&15;		kb = (kb << 5) + jt;
			kd = kc-1; kd += (-(kd < 0))&15;		kc = (kc << 5) + jt;
			ke = kd-1; ke += (-(ke < 0))&15;		kd = (kd << 5) + jt;
													ke = (ke << 5) + jt;
			// 1st set of input ptrs:
			vc0 = tm2+k0;
			vc1 = tm2+k1;
			vc2 = tm2+k2;
			vc3 = tm2+k3;
			vc4 = tm2+k4;
			vc5 = tm2+k5;
			vc6 = tm2+k6;
			vc7 = tm2+k7;
			vc8 = tm2+k8;
			vc9 = tm2+k9;
			vca = tm2+ka;
			vcb = tm2+kb;
			vcc = tm2+kc;
			vcd = tm2+kd;
			vce = tm2+ke;
			
			// Output ptrs:		
			va0 = tmp     ;
			va1 = tmp+0x02;
			va2 = tmp+0x04;
			va3 = tmp+0x06;
			va4 = tmp+0x08;
			va5 = tmp+0x0a;
			va6 = tmp+0x0c;
			va7 = tmp+0x0e;
			va8 = tmp+0x10;
			va9 = tmp+0x12;
			vaa = tmp+0x14;
			vab = tmp+0x16;
			vac = tmp+0x18;
			vad = tmp+0x1a;
			vae = tmp+0x1c;

			SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,
/* s1p00 +... */	vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
/* r00 +... */	va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae
			);
			tmp += 0x1e;	// advance ptr by 15 vec_cmplx [= 30 vec_dbl] elements
		}
	
	  #else	// (OS_BITS == 64):
	
	   #if USE_COMPACT_OBJ_CODE

		tmp = r00;
		tm2 = s1p00;
		for(l = 0; l < 16; l += 2) {
			// NB: All offsets end up being shifted "one too many" bits leftward since
			// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:

			// Low bits of the input-ptr offsets: 2*[0,1,2,...,e,f]
			jt = l << 1;
			k0 = (15-l)&(-(l>0)); /* even-order terms of 0,e,d,...,0, i.e. 0,d,b,9,...,1 */
												// Now get the resulting shifted offsets, adding the low bits to each:
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = (k0 << 5) + jt;
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = (k1 << 5) + jt;
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = (k2 << 5) + jt;
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = (k3 << 5) + jt;
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = (k4 << 5) + jt;
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = (k5 << 5) + jt;
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = (k6 << 5) + jt;
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = (k7 << 5) + jt;
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = (k8 << 5) + jt;
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = (k9 << 5) + jt;
			kb = ka-1; kb += (-(kb < 0))&15;		ka = (ka << 5) + jt;
			kc = kb-1; kc += (-(kc < 0))&15;		kb = (kb << 5) + jt;
			kd = kc-1; kd += (-(kd < 0))&15;		kc = (kc << 5) + jt;
			ke = kd-1; ke += (-(ke < 0))&15;		kd = (kd << 5) + jt;
													ke = (ke << 5) + jt;
			// 1st set of input ptrs:
			vc0 = tm2+k0;
			vc1 = tm2+k1;
			vc2 = tm2+k2;
			vc3 = tm2+k3;
			vc4 = tm2+k4;
			vc5 = tm2+k5;
			vc6 = tm2+k6;
			vc7 = tm2+k7;
			vc8 = tm2+k8;
			vc9 = tm2+k9;
			vca = tm2+ka;
			vcb = tm2+kb;
			vcc = tm2+kc;
			vcd = tm2+kd;
			vce = tm2+ke;
			
			// Now emulate the elided odd-index loop pass (replace l with l+1 in the above sequence) to get next set of ks:
			jt += 2;
			k0 = (14-l); /* odd-order terms of 0,e,d,...,0, i.e. e,c,a,8,...,0
					(Subtract can't yield < 0 here, so no need for mask of result) */
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = (k0 << 5) + jt;
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = (k1 << 5) + jt;
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = (k2 << 5) + jt;
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = (k3 << 5) + jt;
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = (k4 << 5) + jt;
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = (k5 << 5) + jt;
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = (k6 << 5) + jt;
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = (k7 << 5) + jt;
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = (k8 << 5) + jt;
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = (k9 << 5) + jt;
			kb = ka-1; kb += (-(kb < 0))&15;		ka = (ka << 5) + jt;
			kc = kb-1; kc += (-(kc < 0))&15;		kb = (kb << 5) + jt;
			kd = kc-1; kd += (-(kd < 0))&15;		kc = (kc << 5) + jt;
			ke = kd-1; ke += (-(ke < 0))&15;		kd = (kd << 5) + jt;
													ke = (ke << 5) + jt;
			// 2nd set of input ptrs:
			vd0 = tm2+k0;
			vd1 = tm2+k1;
			vd2 = tm2+k2;
			vd3 = tm2+k3;
			vd4 = tm2+k4;
			vd5 = tm2+k5;
			vd6 = tm2+k6;
			vd7 = tm2+k7;
			vd8 = tm2+k8;
			vd9 = tm2+k9;
			vda = tm2+ka;
			vdb = tm2+kb;
			vdc = tm2+kc;
			vdd = tm2+kd;
			vde = tm2+ke;

			// 2 sets of output ptrs:		
							tm1 = tmp + 0x1e;
			va0 = tmp     ;	vb0 = tm1     ;
			va1 = tmp+0x02;	vb1 = tm1+0x02;
			va2 = tmp+0x04;	vb2 = tm1+0x04;
			va3 = tmp+0x06;	vb3 = tm1+0x06;
			va4 = tmp+0x08;	vb4 = tm1+0x08;
			va5 = tmp+0x0a;	vb5 = tm1+0x0a;
			va6 = tmp+0x0c;	vb6 = tm1+0x0c;
			va7 = tmp+0x0e;	vb7 = tm1+0x0e;
			va8 = tmp+0x10;	vb8 = tm1+0x10;
			va9 = tmp+0x12;	vb9 = tm1+0x12;
			vaa = tmp+0x14;	vba = tm1+0x14;
			vab = tmp+0x16;	vbb = tm1+0x16;
			vac = tmp+0x18;	vbc = tm1+0x18;
			vad = tmp+0x1a;	vbd = tm1+0x1a;
			vae = tmp+0x1c;	vbe = tm1+0x1c;

			SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
/* s1p00 +... */	vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
/* r00 +... */	va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,
/* s1p00+jt +... */	vd0,vd1,vd2,vd3,vd4,vd5,vd6,vd7,vd8,vd9,vda,vdb,vdc,vdd,vde,
				y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
/* r0f +... */	vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,vbb,vbc,vbd,vbe
			);
			tmp += 0x3c;	// advance ptr by 30 vec_cmplx [= 60 vec_dbl] elements
		}

	  #else

		/* We use the 16-register/64-bit-mode doubled-radix-15-DFT macros as detailed in the radix60 carry routine: */
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
			s1p00,s1pe0,s1pd0,s1pc0,s1pb0,s1pa0,s1p90,s1p80,s1p70,s1p60,s1p50,s1p40,s1p30,s1p20,s1p10,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e,
			s1pe1,s1pd1,s1pc1,s1pb1,s1pa1,s1p91,s1p81,s1p71,s1p61,s1p51,s1p41,s1p31,s1p21,s1p11,s1p01,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			r0f,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
			s1pd2,s1pc2,s1pb2,s1pa2,s1p92,s1p82,s1p72,s1p62,s1p52,s1p42,s1p32,s1p22,s1p12,s1p02,s1pe2,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			r1e,r1f,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,
			s1pc3,s1pb3,s1pa3,s1p93,s1p83,s1p73,s1p63,s1p53,s1p43,s1p33,s1p23,s1p13,s1p03,s1pe3,s1pd3,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			r2d,r2e,r2f,r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
			s1pb4,s1pa4,s1p94,s1p84,s1p74,s1p64,s1p54,s1p44,s1p34,s1p24,s1p14,s1p04,s1pe4,s1pd4,s1pc4,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			r3c,r3d,r3e,r3f,r40,r41,r42,r43,r44,r45,r46,r47,r48,r49,r4a,
			s1pa5,s1p95,s1p85,s1p75,s1p65,s1p55,s1p45,s1p35,s1p25,s1p15,s1p05,s1pe5,s1pd5,s1pc5,s1pb5,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			r4b,r4c,r4d,r4e,r4f,r50,r51,r52,r53,r54,r55,r56,r57,r58,r59);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
			s1p96,s1p86,s1p76,s1p66,s1p56,s1p46,s1p36,s1p26,s1p16,s1p06,s1pe6,s1pd6,s1pc6,s1pb6,s1pa6,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			r5a,r5b,r5c,r5d,r5e,r5f,r60,r61,r62,r63,r64,r65,r66,r67,r68,
			s1p87,s1p77,s1p67,s1p57,s1p47,s1p37,s1p27,s1p17,s1p07,s1pe7,s1pd7,s1pc7,s1pb7,s1pa7,s1p97,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			r69,r6a,r6b,r6c,r6d,r6e,r6f,r70,r71,r72,r73,r74,r75,r76,r77);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
			s1p78,s1p68,s1p58,s1p48,s1p38,s1p28,s1p18,s1p08,s1pe8,s1pd8,s1pc8,s1pb8,s1pa8,s1p98,s1p88,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			r78,r79,r7a,r7b,r7c,r7d,r7e,r7f,r80,r81,r82,r83,r84,r85,r86,
			s1p69,s1p59,s1p49,s1p39,s1p29,s1p19,s1p09,s1pe9,s1pd9,s1pc9,s1pb9,s1pa9,s1p99,s1p89,s1p79,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			r87,r88,r89,r8a,r8b,r8c,r8d,r8e,r8f,r90,r91,r92,r93,r94,r95);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
			s1p5a,s1p4a,s1p3a,s1p2a,s1p1a,s1p0a,s1pea,s1pda,s1pca,s1pba,s1paa,s1p9a,s1p8a,s1p7a,s1p6a,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			r96,r97,r98,r99,r9a,r9b,r9c,r9d,r9e,r9f,ra0,ra1,ra2,ra3,ra4,
			s1p4b,s1p3b,s1p2b,s1p1b,s1p0b,s1peb,s1pdb,s1pcb,s1pbb,s1pab,s1p9b,s1p8b,s1p7b,s1p6b,s1p5b,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			ra5,ra6,ra7,ra8,ra9,raa,rab,rac,rad,rae,raf,rb0,rb1,rb2,rb3);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
			s1p3c,s1p2c,s1p1c,s1p0c,s1pec,s1pdc,s1pcc,s1pbc,s1pac,s1p9c,s1p8c,s1p7c,s1p6c,s1p5c,s1p4c,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			rb4,rb5,rb6,rb7,rb8,rb9,rba,rbb,rbc,rbd,rbe,rbf,rc0,rc1,rc2,
			s1p2d,s1p1d,s1p0d,s1ped,s1pdd,s1pcd,s1pbd,s1pad,s1p9d,s1p8d,s1p7d,s1p6d,s1p5d,s1p4d,s1p3d,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			rc3,rc4,rc5,rc6,rc7,rc8,rc9,rca,rcb,rcc,rcd,rce,rcf,rd0,rd1);
		SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
			s1p1e,s1p0e,s1pee,s1pde,s1pce,s1pbe,s1pae,s1p9e,s1p8e,s1p7e,s1p6e,s1p5e,s1p4e,s1p3e,s1p2e,
			x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
			rd2,rd3,rd4,rd5,rd6,rd7,rd8,rd9,rda,rdb,rdc,rdd,rde,rdf,re0,
			s1p0f,s1pef,s1pdf,s1pcf,s1pbf,s1paf,s1p9f,s1p8f,s1p7f,s1p6f,s1p5f,s1p4f,s1p3f,s1p2f,s1p1f,
			y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
			re1,re2,re3,re4,re5,re6,re7,re8,re9,rea,reb,rec,red,ree,ref);

	  #endif	// USE_COMPACT_OBJ_CODE?

	 #endif	// (OS_BITS == 32)?
	
		/*...and now do 15 radix-16 transforms, with index-perms as described in radix240_dif_pass1().
		inputs of SSE2_RADIX16_DIT_0TWIDDLE from 30*vec_dbl - separated memlocs, same offsets as already set for DIT: */

	  #if USE_COMPACT_OBJ_CODE

		tmp = r00;
		for(l = 0; l < 15; l++) {
			i64 = dif_perm16[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			addr = &a[j1] + p_in_hi[l];	// p_in_hi[] = p0,p20,...,p30
			add0=addr+k0; add1=addr+k1; add2=addr+k2; add3=addr+k3; add4=addr+k4; add5=addr+k5; add6=addr+k6; add7=addr+k7; add8=addr+k8; add9=addr+k9; adda=addr+ka; addb=addr+kb; addc=addr+kc; addd=addr+kd; adde=addr+ke; addf=addr+kf;
			SSE2_RADIX16_DIF_0TWIDDLE(tmp,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)
			tmp += 2;
		}

	  #else

		addr = &a[j1]    ; add0=addr   ; add1=addr+p1; add2=addr+p2; add3=addr+p3; add4=addr+p4; add5=addr+p5; add6=addr+p6; add7=addr+p7; add8=addr+p8; add9=addr+p9; adda=addr+pa; addb=addr+pb; addc=addr+pc; addd=addr+pd; adde=addr+pe; addf=addr+pf;
		SSE2_RADIX16_DIF_0TWIDDLE(r00,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p20; add0=addr+p5; add1=addr+p4; add2=addr+p7; add3=addr+p6; add4=addr+p3; add5=addr+p2; add6=addr   ; add7=addr+p1; add8=addr+pd; add9=addr+pc; adda=addr+pf; addb=addr+pe; addc=addr+pb; addd=addr+pa; adde=addr+p8; addf=addr+p9;
		SSE2_RADIX16_DIF_0TWIDDLE(r01,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p10; add0=addr+pa; add1=addr+pb; add2=addr+p9; add3=addr+p8; add4=addr+pe; add5=addr+pf; add6=addr+pd; add7=addr+pc; add8=addr+p6; add9=addr+p7; adda=addr+p5; addb=addr+p4; addc=addr+p1; addd=addr   ; adde=addr+p3; addf=addr+p2;
		SSE2_RADIX16_DIF_0TWIDDLE(r02,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+pe0; add0=addr+p7; add1=addr+p6; add2=addr+p4; add3=addr+p5; add4=addr   ; add5=addr+p1; add6=addr+p2; add7=addr+p3; add8=addr+pf; add9=addr+pe; adda=addr+pc; addb=addr+pd; addc=addr+p8; addd=addr+p9; adde=addr+pa; addf=addr+pb;
		SSE2_RADIX16_DIF_0TWIDDLE(r03,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+pd0; add0=addr+p9; add1=addr+p8; add2=addr+pb; add3=addr+pa; add4=addr+pd; add5=addr+pc; add6=addr+pf; add7=addr+pe; add8=addr+p5; add9=addr+p4; adda=addr+p7; addb=addr+p6; addc=addr+p3; addd=addr+p2; adde=addr   ; addf=addr+p1;
		SSE2_RADIX16_DIF_0TWIDDLE(r04,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+pc0; add0=addr+p2; add1=addr+p3; add2=addr+p1; add3=addr   ; add4=addr+p6; add5=addr+p7; add6=addr+p5; add7=addr+p4; add8=addr+pa; add9=addr+pb; adda=addr+p9; addb=addr+p8; addc=addr+pe; addd=addr+pf; adde=addr+pd; addf=addr+pc;
		SSE2_RADIX16_DIF_0TWIDDLE(r05,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+pb0; add0=addr+pb; add1=addr+pa; add2=addr+p8; add3=addr+p9; add4=addr+pf; add5=addr+pe; add6=addr+pc; add7=addr+pd; add8=addr+p7; add9=addr+p6; adda=addr+p4; addb=addr+p5; addc=addr   ; addd=addr+p1; adde=addr+p2; addf=addr+p3;
		SSE2_RADIX16_DIF_0TWIDDLE(r06,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+pa0; add0=addr+p1; add1=addr   ; add2=addr+p3; add3=addr+p2; add4=addr+p5; add5=addr+p4; add6=addr+p7; add7=addr+p6; add8=addr+p9; add9=addr+p8; adda=addr+pb; addb=addr+pa; addc=addr+pd; addd=addr+pc; adde=addr+pf; addf=addr+pe;
		SSE2_RADIX16_DIF_0TWIDDLE(r07,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p90; add0=addr+pc; add1=addr+pd; add2=addr+pe; add3=addr+pf; add4=addr+pa; add5=addr+pb; add6=addr+p9; add7=addr+p8; add8=addr+p2; add9=addr+p3; adda=addr+p1; addb=addr   ; addc=addr+p6; addd=addr+p7; adde=addr+p5; addf=addr+p4;
		SSE2_RADIX16_DIF_0TWIDDLE(r08,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p80; add0=addr+p3; add1=addr+p2; add2=addr   ; add3=addr+p1; add4=addr+p7; add5=addr+p6; add6=addr+p4; add7=addr+p5; add8=addr+pb; add9=addr+pa; adda=addr+p8; addb=addr+p9; addc=addr+pf; addd=addr+pe; adde=addr+pc; addf=addr+pd;
		SSE2_RADIX16_DIF_0TWIDDLE(r09,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p70; add0=addr+pe; add1=addr+pf; add2=addr+pd; add3=addr+pc; add4=addr+p9; add5=addr+p8; add6=addr+pb; add7=addr+pa; add8=addr+p1; add9=addr   ; adda=addr+p3; addb=addr+p2; addc=addr+p5; addd=addr+p4; adde=addr+p7; addf=addr+p6;
		SSE2_RADIX16_DIF_0TWIDDLE(r0a,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p60; add0=addr+p4; add1=addr+p5; add2=addr+p6; add3=addr+p7; add4=addr+p2; add5=addr+p3; add6=addr+p1; add7=addr   ; add8=addr+pc; add9=addr+pd; adda=addr+pe; addb=addr+pf; addc=addr+pa; addd=addr+pb; adde=addr+p9; addf=addr+p8;
		SSE2_RADIX16_DIF_0TWIDDLE(r0b,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p50; add0=addr+pd; add1=addr+pc; add2=addr+pf; add3=addr+pe; add4=addr+pb; add5=addr+pa; add6=addr+p8; add7=addr+p9; add8=addr+p3; add9=addr+p2; adda=addr   ; addb=addr+p1; addc=addr+p7; addd=addr+p6; adde=addr+p4; addf=addr+p5;
		SSE2_RADIX16_DIF_0TWIDDLE(r0c,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p40; add0=addr+p6; add1=addr+p7; add2=addr+p5; add3=addr+p4; add4=addr+p1; add5=addr   ; add6=addr+p3; add7=addr+p2; add8=addr+pe; add9=addr+pf; adda=addr+pd; addb=addr+pc; addc=addr+p9; addd=addr+p8; adde=addr+pb; addf=addr+pa;
		SSE2_RADIX16_DIF_0TWIDDLE(r0d,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

		addr = &a[j1]+p30; add0=addr+p8; add1=addr+p9; add2=addr+pa; add3=addr+pb; add4=addr+pc; add5=addr+pd; add6=addr+pe; add7=addr+pf; add8=addr+p4; add9=addr+p5; adda=addr+p6; addb=addr+p7; addc=addr+p2; addd=addr+p3; adde=addr+p1; addf=addr   ;
		SSE2_RADIX16_DIF_0TWIDDLE(r0e,OFF1,OFF2,OFF3,OFF4, isrt2, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf)

	   #endif	// USE_COMPACT_OBJ_CODE

	#else	/* !USE_SSE2 */

		tptr = t;
		for(ntmp = 0; ntmp < 16; ntmp++) {
		// [1] here:
			k0 = plo[ntmp];	// p0,..,f
			jt = j1 + k0; jp = j2 + k0;
		// [2] here:
			k0 = (15-ntmp)&(-(ntmp>0)); /* 0,e,d,...,0 */	// Now get the resulting p* offsets:
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = p_out_hi[k0];
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = p_out_hi[k1];
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = p_out_hi[k2];
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = p_out_hi[k3];
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = p_out_hi[k4];
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = p_out_hi[k5];
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = p_out_hi[k6];
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = p_out_hi[k7];
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = p_out_hi[k8];
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = p_out_hi[k9];
			kb = ka-1; kb += (-(kb < 0))&15;		ka = p_out_hi[ka];
			kc = kb-1; kc += (-(kc < 0))&15;		kb = p_out_hi[kb];
			kd = kc-1; kd += (-(kd < 0))&15;		kc = p_out_hi[kc];
			ke = kd-1; ke += (-(ke < 0))&15;		kd = p_out_hi[kd];
													ke = p_out_hi[ke];
			RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im
			);	tptr += 0xf;
		}
		//...and now do 15 radix-16 transforms:
		tptr = t;
		for(ntmp = 0; ntmp < 15; ntmp++) {
			i64 = dif_perm16[ntmp];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			jt = j1 + p_in_hi[ntmp]; jp = j2 + p_in_hi[ntmp];	// p_in_hi[] = p0,p20,...,p30
			RADIX_16_DIF(
				tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				c16,s16);	tptr++;
		}

	#endif	/* if(USE_SSE2) */
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(k=1; k <= khi; k++) */

