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

	/*...The radix-240 DIT pass is here:	*/

	#ifdef USE_SSE2

	  #ifdef USE_ARM_V8_SIMD
		const uint32 OFF1 = 0x1e0;
		const uint32 OFF2 = 0x3c0;
		const uint32 OFF3 = 0x5a0;
		const uint32 OFF4 = 0x780;
	  #elif defined(USE_AVX512)
		#define OFF1	0x1e0*4
		#define OFF2	0x3c0*4
		#define OFF3	0x5a0*4
		#define OFF4	0x780*4
	  #elif defined(USE_AVX)
		#define OFF1	0x1e0*2
		#define OFF2	0x3c0*2
		#define OFF3	0x5a0*2
		#define OFF4	0x780*2
	  #else
		#define OFF1	0x1e0
		#define OFF2	0x3c0
		#define OFF3	0x5a0
		#define OFF4	0x780
	  #endif

		tmp = r00;
		for(l = 0; l < 15; l++) {
			i64 = dit_perm16[l];
			addr = &a[j1] + p_in_hi[l];	// p_in_hi[] = p0,p20,...,p30
			// p-offset indices encoded in little-endian hex-char fashion:
			po_kperm[0x0] = plo[(i64 >> 60)&0xf];
			po_kperm[0x1] = plo[(i64 >> 56)&0xf];
			po_kperm[0x2] = plo[(i64 >> 52)&0xf];
			po_kperm[0x3] = plo[(i64 >> 48)&0xf];
			po_kperm[0x4] = plo[(i64 >> 44)&0xf];
			po_kperm[0x5] = plo[(i64 >> 40)&0xf];
			po_kperm[0x6] = plo[(i64 >> 36)&0xf];
			po_kperm[0x7] = plo[(i64 >> 32)&0xf];
			po_kperm[0x8] = plo[(i64 >> 28)&0xf];
			po_kperm[0x9] = plo[(i64 >> 24)&0xf];
			po_kperm[0xa] = plo[(i64 >> 20)&0xf];
			po_kperm[0xb] = plo[(i64 >> 16)&0xf];
			po_kperm[0xc] = plo[(i64 >> 12)&0xf];
			po_kperm[0xd] = plo[(i64 >>  8)&0xf];
			po_kperm[0xe] = plo[(i64 >>  4)&0xf];
			po_kperm[0xf] = plo[(i64      )&0xf];
			SSE2_RADIX16_DIT_0TWIDDLE(
				addr,po_ptr,
				isrt2,two,
				tmp,OFF1,OFF2,OFF3,OFF4
			);	tmp += 2;
		}

		/*...and now do 16 radix-15 transforms. See radix240_dit_pass1() for details on the oindex patterning.
		In 64-bit mode we use the 16-register/64-bit-mode doubled-radix-15-DFT macros as detailed in the radix60 carry routine.
		*/
	  #if !RAD_15_2FOLD
		#warning Using 1-fold 15-DFT
		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,*va9,*vaa,*vab,*vac,*vad,*vae,
		*vc0,*vc1,*vc2,*vc3,*vc4,*vc5,*vc6,*vc7,*vc8,*vc9,*vca,*vcb,*vcc,*vcd,*vce;

		tmp = r00;
		tm2 = s1p00;
		for(l = 0; l < 16; l++) {
			// NB: All offsets end up being shifted "one too many" bits leftward since
			// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:

		#if 1	// HLL/ASM toggle

			// Low bits of the output-ptr offsets: 2*[0,f,e,...,2,1]
			jt = ((16 - l)&0xf) << 1;	// 2*[0,f,e,...,2,1]
			k0 = (l-1)&(-(l>0)); 		// 0,0,1,2...e
												// Now get the resulting shifted offsets, adding the low bits to each:
			k1 = k0-1; k1 += (-(k1 < 0))&15;	k0 = (k0 << 5) + jt;
			k2 = k1-1; k2 += (-(k2 < 0))&15;	k1 = (k1 << 5) + jt;
			k3 = k2-1; k3 += (-(k3 < 0))&15;	k2 = (k2 << 5) + jt;
			k4 = k3-1; k4 += (-(k4 < 0))&15;	k3 = (k3 << 5) + jt;
			k5 = k4-1; k5 += (-(k5 < 0))&15;	k4 = (k4 << 5) + jt;
			k6 = k5-1; k6 += (-(k6 < 0))&15;	k5 = (k5 << 5) + jt;
			k7 = k6-1; k7 += (-(k7 < 0))&15;	k6 = (k6 << 5) + jt;
			k8 = k7-1; k8 += (-(k8 < 0))&15;	k7 = (k7 << 5) + jt;
			k9 = k8-1; k9 += (-(k9 < 0))&15;	k8 = (k8 << 5) + jt;
			ka = k9-1; ka += (-(ka < 0))&15;	k9 = (k9 << 5) + jt;
			kb = ka-1; kb += (-(kb < 0))&15;	ka = (ka << 5) + jt;
			kc = kb-1; kc += (-(kc < 0))&15;	kb = (kb << 5) + jt;
			kd = kc-1; kd += (-(kd < 0))&15;	kc = (kc << 5) + jt;
			ke = kd-1; ke += (-(ke < 0))&15;	kd = (kd << 5) + jt;
												ke = (ke << 5) + jt;
		//	printf("15-DFT #%2u: [k0-E]/2 = %u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u\n",l,k0/2,k1/2,k2/2,k3/2,k4/2,k5/2,k6/2,k7/2,k8/2,k9/2,ka/2,kb/2,kc/2,kd/2,ke/2);
		//	printf("0x0%2X%2X%2X%2X%2X%2X%2X,%#2X%2X%2X%2X%2X%2X%2X%2X\n",ke/2,kd/2,kc/2,kb/2,ka/2,k9/2,k8/2,k7/2,k6/2,k5/2,k4/2,k3/2,k2/2,k1/2,k0/2);
			// Input ptrs:		// Output ptrs:
			va0 = tmp     ;		vc0 = tm2 + k0;
			va1 = tmp+0x02;		vc1 = tm2 + k1;
			va2 = tmp+0x04;		vc2 = tm2 + k2;
			va3 = tmp+0x06;		vc3 = tm2 + k3;
			va4 = tmp+0x08;		vc4 = tm2 + k4;
			va5 = tmp+0x0a;		vc5 = tm2 + k5;
			va6 = tmp+0x0c;		vc6 = tm2 + k6;
			va7 = tmp+0x0e;		vc7 = tm2 + k7;
			va8 = tmp+0x10;		vc8 = tm2 + k8;
			va9 = tmp+0x12;		vc9 = tm2 + k9;
			vaa = tmp+0x14;		vca = tm2 + ka;
			vab = tmp+0x16;		vcb = tm2 + kb;
			vac = tmp+0x18;		vcc = tm2 + kc;
			vad = tmp+0x1a;		vcd = tm2 + kd;
			vae = tmp+0x1c;		vce = tm2 + ke;

		#else

		  #ifdef FOO//USE_AVX512
			#warning This experimental code barfs with fatal ROE after several iters for reasons unknown, and needs debug!
			// Same table as below, but in uint64 form, and the 16 bytes within the two elements
			// of each row corresponding to 0x[k7,k6,k5,k4,k3,k2,k1,k0], 0x[0,ke,kd,kc,kb,ka,k8]:
			const uint64 dit15_offs[32] = {
				0x8090A0B0C0D0E000ull,0x0010203040506070ull,
				0x8F9FAFBFCFDFEF0Full,0x001F2F3F4F5F6F7Full,
				0x9EAEBECEDEEE0E1Eull,0x002E3E4E5E6E7E8Eull,
				0xADBDCDDDED0D1D2Dull,0x003D4D5D6D7D8D9Dull,
				0xBCCCDCEC0C1C2C3Cull,0x004C5C6C7C8C9CACull,
				0xCBDBEB0B1B2B3B4Bull,0x005B6B7B8B9BABBBull,
				0xDAEA0A1A2A3A4A5Aull,0x006A7A8A9AAABACAull,
				0xE909192939495969ull,0x00798999A9B9C9D9ull,
				0x0818283848586878ull,0x008898A8B8C8D8E8ull,
				0x1727374757677787ull,0x0097A7B7C7D7E707ull,
				0x2636465666768696ull,0x00A6B6C6D6E60616ull,
				0x35455565758595A5ull,0x00B5C5D5E5051525ull,
				0x445464748494A4B4ull,0x00C4D4E404142434ull,
				0x5363738393A3B3C3ull,0x00D3E30313233343ull,
				0x62728292A2B2C2D2ull,0x00E2021222324252ull,
				0x718191A1B1C1D1E1ull,0x0001112131415161ull
			};
			// Use spare slots of local-alloc aligned strorage for ptr64_vecs, so MOVs in ASM below can be aligned.
			// Each of these arrays-of-pointers needs two 512-bit chunks, hence the 1/3 increments w.r.to sinwtm1:
			vec_dbl **pvec1 = (vec_dbl*)sinwtm1+1, **pvec2 = (vec_dbl*)sinwtm1+3;
			uint64 *ptr64 = dit15_offs + (l<<1);
			__asm__ volatile (\
				"movq	%[ptr64],%%rax		\n\t"\
			"movq $0x0E0C0A0806040200,%%rbx	\n\t"/* 64-bit register w/byte offsets  0:14:2, bytes ordered left-to-right in decreasing significance */\
			"movq $0x001C1A1816141210,%%rcx	\n\t"/* 64-bit register w/byte offsets 16:28:2, high byte = 0 */\
				"vmovq	   (%%rax),%%xmm0	\n\t"/* Load word1 into low qword (64 bits) of zmm8 [NB: avx-512 only supports MOVQ to/from 128-bit vector regs] */\
				"vmovq	0x8(%%rax),%%xmm1	\n\t"/* word2 */\
				"vmovq		%%rbx ,%%xmm2	\n\t"\
				"vmovq		%%rcx ,%%xmm3	\n\t"\
				"vpmovzxbq	%%xmm0,%%zmm0	\n\t"/* each word1 byte ==> low byte of qword */\
				"vpmovzxbq	%%xmm1,%%zmm1	\n\t"/* each word2 byte ==> low byte of qword */\
				"vpmovzxbq	%%xmm2,%%zmm2	\n\t"\
				"vpmovzxbq	%%xmm3,%%zmm3	\n\t"\
				"vpslld	$6,%%zmm0,%%zmm0	\n\t"/* << 6 to turn 0:28:2 index-offsets into proper 64-byte vec_dbl-pointer offsets */\
				"vpslld	$6,%%zmm1,%%zmm1	\n\t"\
				"vpslld	$6,%%zmm2,%%zmm2	\n\t"\
				"vpslld	$6,%%zmm3,%%zmm3	\n\t"\
				"movq	%[tmp],%%rbx		\n\t	vpbroadcastq  %%rbx,%%zmm4	\n\t"\
				"movq	%[tm2],%%rcx		\n\t	vpbroadcastq  %%rcx,%%zmm5	\n\t"\
			"vpaddq	%%zmm0,%%zmm0,%%zmm0	\n\t"/* 2*k0-7 */\
			"vpaddq	%%zmm1,%%zmm1,%%zmm1	\n\t"/* 2*k8-e */\
			/* Note: load-with-broadcast-from-mem64, e.g. "vpaddq (%%rbx)%{1to8%},%%zmm0,%%zmm0", not what
			we want for the tmp/tm2-base-pointers here, since in the present context tmp/tm2 are not being
			treated as memory locations but rather as 64-bit ints which we want to broadcast and vector-add to: */\
			"vpaddq %%zmm5,%%zmm0,%%zmm0	\n\t"/* vc0-7 */\
			"vpaddq %%zmm5,%%zmm1,%%zmm1	\n\t"/* vc8-e */\
			"vpaddq %%zmm4,%%zmm2,%%zmm2	\n\t"/* va0-7 */\
			"vpaddq %%zmm4,%%zmm3,%%zmm3	\n\t"/* va8-e */\
				/* write to ptr64-vecs: */
				"movq	%[pvec2],%%rax	\n\t"/* vc-ptrs are outputs for 15-DIT */\
				"movq	%[pvec1],%%rbx	\n\t"\
				"vmovaps %%zmm0,    (%%rax)	\n\t"/* pvec2[0:7] = vc0-7 */\
				"vmovaps %%zmm1,0x40(%%rax)	\n\t"/* pvec2[8:e] = vc8-e */\
				"vmovaps %%zmm2,    (%%rbx)	\n\t"/* pvec1[0:7] = va0-7 */\
				"vmovaps %%zmm3,0x40(%%rbx)	\n\t"/* pvec1[8:e] = va8-e */\
				:					// outputs: none
				: [tmp] "m" (tmp)	// All inputs from memory addresses here
				 ,[tm2] "m" (tm2)
				 ,[ptr64] "m" (ptr64)
				 ,[pvec1] "m" (pvec1)
				 ,[pvec2] "m" (pvec2)
				: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3"		/* Clobbered registers */\
			);
		//	printf("ASM I-ptr offs[%2u] = %u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u\n",l,pvec1[0x0]-r00,pvec1[0x1]-r00,pvec1[0x2]-r00,pvec1[0x3]-r00,pvec1[0x4]-r00,pvec1[0x5]-r00,pvec1[0x6]-r00,pvec1[0x7]-r00,pvec1[0x8]-r00,pvec1[0x9]-r00,pvec1[0xa]-r00,pvec1[0xb]-r00,pvec1[0xc]-r00,pvec1[0xd]-r00,pvec1[0xe]-r00);
		//	printf("ASM O-ptr offs[%2u] = %u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u\n",l,pvec2[0x0]-tm2,pvec2[0x1]-tm2,pvec2[0x2]-tm2,pvec2[0x3]-tm2,pvec2[0x4]-tm2,pvec2[0x5]-tm2,pvec2[0x6]-tm2,pvec2[0x7]-tm2,pvec2[0x8]-tm2,pvec2[0x9]-tm2,pvec2[0xa]-tm2,pvec2[0xb]-tm2,pvec2[0xc]-tm2,pvec2[0xd]-tm2,pvec2[0xe]-tm2);
			va0 = pvec1[0x0];va1 = pvec1[0x1];va2 = pvec1[0x2];va3 = pvec1[0x3];va4 = pvec1[0x4];va5 = pvec1[0x5];va6 = pvec1[0x6];va7 = pvec1[0x7];va8 = pvec1[0x8];va9 = pvec1[0x9];vaa = pvec1[0xa];vab = pvec1[0xb];vac = pvec1[0xc];vad = pvec1[0xd];vae = pvec1[0xe];
			vc0 = pvec2[0x0];vc1 = pvec2[0x1];vc2 = pvec2[0x2];vc3 = pvec2[0x3];vc4 = pvec2[0x4];vc5 = pvec2[0x5];vc6 = pvec2[0x6];vc7 = pvec2[0x7];vc8 = pvec2[0x8];vc9 = pvec2[0x9];vca = pvec2[0xa];vcb = pvec2[0xb];vcc = pvec2[0xc];vcd = pvec2[0xd];vce = pvec2[0xe];

		  #else

			const uint8 dit15_offs[RADIX] = {
				0,224,208,192,176,160,144,128,112,96,80,64,48,32,16,
				15,239,223,207,191,175,159,143,127,111,95,79,63,47,31,
				30,14,238,222,206,190,174,158,142,126,110,94,78,62,46,
				45,29,13,237,221,205,189,173,157,141,125,109,93,77,61,
				60,44,28,12,236,220,204,188,172,156,140,124,108,92,76,
				75,59,43,27,11,235,219,203,187,171,155,139,123,107,91,
				90,74,58,42,26,10,234,218,202,186,170,154,138,122,106,
				105,89,73,57,41,25,9,233,217,201,185,169,153,137,121,
				120,104,88,72,56,40,24,8,232,216,200,184,168,152,136,
				135,119,103,87,71,55,39,23,7,231,215,199,183,167,151,
				150,134,118,102,86,70,54,38,22,6,230,214,198,182,166,
				165,149,133,117,101,85,69,53,37,21,5,229,213,197,181,
				180,164,148,132,116,100,84,68,52,36,20,4,228,212,196,
				195,179,163,147,131,115,99,83,67,51,35,19,3,227,211,
				210,194,178,162,146,130,114,98,82,66,50,34,18,2,226,
				225,209,193,177,161,145,129,113,97,81,65,49,33,17,1
			};
			i = 15*l;
			k0 = dit15_offs[i+ 0];k1 = dit15_offs[i+ 1];k2 = dit15_offs[i+ 2];k3 = dit15_offs[i+ 3];k4 = dit15_offs[i+ 4];
			k5 = dit15_offs[i+ 5];k6 = dit15_offs[i+ 6];k7 = dit15_offs[i+ 7];k8 = dit15_offs[i+ 8];k9 = dit15_offs[i+ 9];
			ka = dit15_offs[i+10];kb = dit15_offs[i+11];kc = dit15_offs[i+12];kd = dit15_offs[i+13];ke = dit15_offs[i+14];
			// Input ptrs:		// Output ptrs:
			va0 = tmp     ;		vc0 = tm2 + (k0<<1);
			va1 = tmp+0x02;		vc1 = tm2 + (k1<<1);
			va2 = tmp+0x04;		vc2 = tm2 + (k2<<1);
			va3 = tmp+0x06;		vc3 = tm2 + (k3<<1);
			va4 = tmp+0x08;		vc4 = tm2 + (k4<<1);
			va5 = tmp+0x0a;		vc5 = tm2 + (k5<<1);
			va6 = tmp+0x0c;		vc6 = tm2 + (k6<<1);
			va7 = tmp+0x0e;		vc7 = tm2 + (k7<<1);
			va8 = tmp+0x10;		vc8 = tm2 + (k8<<1);
			va9 = tmp+0x12;		vc9 = tm2 + (k9<<1);
			vaa = tmp+0x14;		vca = tm2 + (ka<<1);
			vab = tmp+0x16;		vcb = tm2 + (kb<<1);
			vac = tmp+0x18;		vcc = tm2 + (kc<<1);
			vad = tmp+0x1a;		vcd = tm2 + (kd<<1);
			vae = tmp+0x1c;		vce = tm2 + (ke<<1);
		//	printf("HLL I-ptr offs[%2u] = %u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u\n",l,va0-r00,va1-r00,va2-r00,va3-r00,va4-r00,va5-r00,va6-r00,va7-r00,va8-r00,va9-r00,vaa-r00,vab-r00,vac-r00,vad-r00,vae-r00);
		//	printf("HLL O-ptr offs[%2u] = %u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u\n",l,vc0-tm2,vc1-tm2,vc2-tm2,vc3-tm2,vc4-tm2,vc5-tm2,vc6-tm2,vc7-tm2,vc8-tm2,vc9-tm2,vca-tm2,vcb-tm2,vcc-tm2,vcd-tm2,vce-tm2);

		  #endif	// USE_AVX512 ?

		#endif	// HLL/ASM toggle

			SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1,
/* r00 +... */	va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
/* s1p00 +... */	vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce
			);
			tmp += 0x1e;	// advance ptr by 15 vec_cmplx [= 30 vec_dbl] elements
		}

	  #else	// RAD_15_2FOLD == True:
		#warning Using 2-fold 15-DFT
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
												// Shifted offsets:	// 1st set of Output ptrs:
			k1 = k0-1; k1 += (-(k1 < 0))&15;	k0 = (k0 << 5) + jt;	vc0 = tm2 + k0;
			k2 = k1-1; k2 += (-(k2 < 0))&15;	k1 = (k1 << 5) + jt;	vc1 = tm2 + k1;
			k3 = k2-1; k3 += (-(k3 < 0))&15;	k2 = (k2 << 5) + jt;	vc2 = tm2 + k2;
			k4 = k3-1; k4 += (-(k4 < 0))&15;	k3 = (k3 << 5) + jt;	vc3 = tm2 + k3;
			k5 = k4-1; k5 += (-(k5 < 0))&15;	k4 = (k4 << 5) + jt;	vc4 = tm2 + k4;
			k6 = k5-1; k6 += (-(k6 < 0))&15;	k5 = (k5 << 5) + jt;	vc5 = tm2 + k5;
			k7 = k6-1; k7 += (-(k7 < 0))&15;	k6 = (k6 << 5) + jt;	vc6 = tm2 + k6;
			k8 = k7-1; k8 += (-(k8 < 0))&15;	k7 = (k7 << 5) + jt;	vc7 = tm2 + k7;
			k9 = k8-1; k9 += (-(k9 < 0))&15;	k8 = (k8 << 5) + jt;	vc8 = tm2 + k8;
			ka = k9-1; ka += (-(ka < 0))&15;	k9 = (k9 << 5) + jt;	vc9 = tm2 + k9;
			kb = ka-1; kb += (-(kb < 0))&15;	ka = (ka << 5) + jt;	vca = tm2 + ka;
			kc = kb-1; kc += (-(kc < 0))&15;	kb = (kb << 5) + jt;	vcb = tm2 + kb;
			kd = kc-1; kd += (-(kd < 0))&15;	kc = (kc << 5) + jt;	vcc = tm2 + kc;
			ke = kd-1; ke += (-(ke < 0))&15;	kd = (kd << 5) + jt;	vcd = tm2 + kd;
												ke = (ke << 5) + jt;	vce = tm2 + ke;

			// Now emulate the elided odd-index loop pass (replace l with l+1 in the above sequence) to get next set of ks:
			jt = (15 - l) << 1;	// Subtract can't yield < 0 here, so no need for mask of result
			k0 = l; // odd-order terms of 0,0,1,2...e, i.e. 0,2,4,...	// 2nd set of Output ptrs:
			k1 = k0-1; k1 += (-(k1 < 0))&15;	k0 = (k0 << 5) + jt;	vd0 = tm2 + k0;
			k2 = k1-1; k2 += (-(k2 < 0))&15;	k1 = (k1 << 5) + jt;	vd1 = tm2 + k1;
			k3 = k2-1; k3 += (-(k3 < 0))&15;	k2 = (k2 << 5) + jt;	vd2 = tm2 + k2;
			k4 = k3-1; k4 += (-(k4 < 0))&15;	k3 = (k3 << 5) + jt;	vd3 = tm2 + k3;
			k5 = k4-1; k5 += (-(k5 < 0))&15;	k4 = (k4 << 5) + jt;	vd4 = tm2 + k4;
			k6 = k5-1; k6 += (-(k6 < 0))&15;	k5 = (k5 << 5) + jt;	vd5 = tm2 + k5;
			k7 = k6-1; k7 += (-(k7 < 0))&15;	k6 = (k6 << 5) + jt;	vd6 = tm2 + k6;
			k8 = k7-1; k8 += (-(k8 < 0))&15;	k7 = (k7 << 5) + jt;	vd7 = tm2 + k7;
			k9 = k8-1; k9 += (-(k9 < 0))&15;	k8 = (k8 << 5) + jt;	vd8 = tm2 + k8;
			ka = k9-1; ka += (-(ka < 0))&15;	k9 = (k9 << 5) + jt;	vd9 = tm2 + k9;
			kb = ka-1; kb += (-(kb < 0))&15;	ka = (ka << 5) + jt;	vda = tm2 + ka;
			kc = kb-1; kc += (-(kc < 0))&15;	kb = (kb << 5) + jt;	vdb = tm2 + kb;
			kd = kc-1; kd += (-(kd < 0))&15;	kc = (kc << 5) + jt;	vdc = tm2 + kc;
			ke = kd-1; ke += (-(ke < 0))&15;	kd = (kd << 5) + jt;	vdd = tm2 + kd;
												ke = (ke << 5) + jt;	vde = tm2 + ke;
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

			SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,two,
/* r00 +... */	va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
/* s1p00 +... */	vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,
/* r0f +... */	vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,vbb,vbc,vbd,vbe,
				y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
/* s1p00+jt +... */	vd0,vd1,vd2,vd3,vd4,vd5,vd6,vd7,vd8,vd9,vda,vdb,vdc,vdd,vde
			);
			tmp += 0x3c;	// advance ptr by 30 vec_cmplx [= 60 vec_dbl] elements
		}

	 #endif	// RAD_15_2FOLD ?

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

	  if(incr) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of incr = incr_hiacc
					// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		uint32 ii,loop, co2save = co2;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every 8th pass or better:
		// incr must divide nloop [RADIX/8 = 30 or RADIX/16 = 15, depending on whether we use 8-or-16-way carry macros]!
	  #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
	  #else
		const uint32 nloop = RADIX>>3;
	  #endif
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy_r; itmp = bjmodn;
	  #ifndef USE_AVX512
		tm2 = cy_r+1; itm2 = bjmodn+4;	// tm2,itm2 not used in AVX-512 mode
	  #endif
		for(loop = 0; loop < nloop; loop += incr)
		{
			co2 = co2save;	// Need this for all wts-inits beynd the initial set, due to the co2 = co3 preceding the (j+2) data
		  #ifdef CARRY_16_WAY
			ii = loop << 4;	// Reflects 16 independent carry chains being done in each AVX_cmplx_carry_fast_errcheck_X8 call
		  #else
			ii = loop << 3;	// Reflects  8 independent carry chains being done in each AVX_cmplx_carry_fast_errcheck_X8 call
		  #endif
			add1 = &wt1[col  +ii];	/* Don't use add0 here, to avoid need to reload main-array address */
			add2 = &wt1[co2-1-ii];
			add3 = &wt1[co3-1-ii];

			// Since use wt1-array in the wtsinit macro, need to fiddle this here:
			co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
			// *But*: since the init macro does an on-the-fly version of this between j,j+2 portions, external code co2=co3 must come *after* both ctmp-data octets are inited.
		  #ifdef CARRY_16_WAY
			AVX_cmplx_carry_fast_wtsinit_X16(add1,add2,add3, itmp, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
		  #else
			AVX_cmplx_carry_fast_wtsinit_X8 (add1,add2,add3, itmp, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
		  #endif
			for(l = loop; l < loop+incr; l++) {
				// Each AVX carry macro call also processes 8 prefetches of main-array data
				add0 = a + j1 + pfetch_dist + poff[l+l];
			  // In AVX-512 mode, the 4 doubles base[0],baseinv[1],wts_mult[1],inv_mult[0] are in the d0-3 slots of the otherwise-unused sse2_rnd vec_dbl:
			  #ifdef USE_AVX512
			   #ifdef CARRY_16_WAY
				AVX_cmplx_carry_fast_errcheck_X16(tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 32; tm1 += 2;           itmp += 16;           i = 0;
			   #else
				AVX_cmplx_carry_fast_errcheck_X8 (tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 16; tm1 += 1;           itmp +=  8;           i = 0;	// CY-ptr only advances 1 in AVX-512/CARRY_8_WAY mode, since all 8 dbl-carries fit in a single vec_dbl
			   #endif
			  #else	// USE_AVX:
				AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 16; tm1 += 2; tm2 += 2; itmp += 8; itm2 += 8; i = 0;
			  #endif
			}
		}

	  } else {	// HiACC:

		/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		tm1 = s1p00; tmp = cy_r; itmp = bjmodn;
		addr = &prp_mult;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[l]);	// poff[] = p0,4,8,...
			AVX_cmplx_carry_norm_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, tm2,p1,p2,p3, addr);
			tm1 += 8; tmp += 1; itmp += 4; i = 0;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	  if(incr) {	// Have no specialized HIACC carry macro in ARM_V8_SIMD, so use 0-or-not-ness of incr
						// in lieu of (USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) is for ARMv8
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		uint32 i0,i1,i2,i3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
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
			i0 = n-si[l  ];
			i1 = n-si[l+1];
			i2 = si[nwtml  ];
			i3 = si[nwtml-1];
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
			SSE2_cmplx_carry_fast_wtsinit(add1,add2,add3, add0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, i0,i1,i2,i3, sse_bw,sse_n)

			for(l = loop; l < loop+incr; l++) {
				// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data:
				add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
				SSE2_cmplx_carry_fast_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr);
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

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		addr = &prp_mult;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1, addr);
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

	/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p2,p3, addr);
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
			cmplx_carry_norm_errcheck0(a[jt   ],a[jp   ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p1],a[jp+p1],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p2],a[jp+p2],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p3],a[jp+p3],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		}

	  } else {	// HiACC:

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; addr = cy_r; itmp = bjmodn;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			cmplx_carry_norm_errcheck0(a[jt   ],a[jp   ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p1],a[jp+p1],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p2],a[jp+p2],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p3],a[jp+p3],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
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

		/* The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
		which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls: */

		tm0 = s1p00; tmp = base_negacyclic_root; l = 0x3c00;
		tm1 = cy_r; // tm2 = cy_i;	*** replace with literal-byte-offset in macro call to save a reg
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic_idx = 0; jc_idx = 1; kc_idx = 2; lc_idx = 3;
		for(ntmp = 0; ntmp < RADIX; ntmp += 4)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{
			k1 = icycle[ic_idx];	k5 = jcycle[ic_idx];	k6 = kcycle[ic_idx];	k7 = lcycle[ic_idx];
			k2 = icycle[jc_idx];
			k3 = icycle[kc_idx];
			k4 = icycle[lc_idx];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x780, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p1,p2,p3, addr);
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

			for(l = 0; l < inc_arr[ntmp]; l++) {
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
				SSE2_fermat_carry_norm_errcheck_X8_loacc(tm0,tmp,tm1,0x780, 0x3c0,0x780,0xb40, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, tm2,p1,p2,p3,p4, addr);
				tm0 += 16; tm1++;
				MOD_ADD32(ic_idx, 8, ODD_RADIX, ic_idx);
				MOD_ADD32(jc_idx, 8, ODD_RADIX, jc_idx);
				MOD_ADD32(kc_idx, 8, ODD_RADIX, kc_idx);
				MOD_ADD32(lc_idx, 8, ODD_RADIX, lc_idx);
				MOD_ADD32(mc_idx, 8, ODD_RADIX, mc_idx);
				MOD_ADD32(nc_idx, 8, ODD_RADIX, nc_idx);
				MOD_ADD32(oc_idx, 8, ODD_RADIX, oc_idx);
				MOD_ADD32(pc_idx, 8, ODD_RADIX, pc_idx);
			}

		#else	// AVX / AVX2

			for(l = 0; l < inc_arr[ntmp]; l++) {
				k1 = icycle[ic_idx];
				k2 = icycle[jc_idx];	k5 = jcycle[ic_idx];
				k3 = icycle[kc_idx];	k6 = kcycle[ic_idx];
				k4 = icycle[lc_idx];	k7 = lcycle[ic_idx];
				// Each AVX carry macro call also processes 4 prefetches of main-array data
				tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
													/* (cy_i_cy_r) --vvvvv  vvvvvvvvvvvvvvvvv--[1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
				SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x780, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p1,p2,p3, addr);
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
		while((int)(tmp-cy_r) < RADIX) {
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic_idx];
			k2 = jcycle[ic_idx];
			k3 = icycle[jc_idx];
			k4 = jcycle[jc_idx];
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tmp-cy_r)>>2]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
			tm2 += (-((int)((tmp-cy_r)>>1)&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_fermat_carry_norm_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2,k3,k4, tm2,p1, addr);
			tm1 += 4; tmp += 2;
			MOD_ADD32(ic_idx, 2, ODD_RADIX, ic_idx);
			MOD_ADD32(jc_idx, 2, ODD_RADIX, jc_idx);
		}

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		ntmp = 0; addr = cy_r; addi = cy_i; ic_idx = 0;	// ic_idx = idx into icycle mini-array, gets incremented (mod ODD_RADIX) between macro calls
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];
			fermat_carry_norm_errcheckB(a[jt   ],a[jp   ],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p1],a[jp+p1],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p2],a[jp+p2],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
			fermat_carry_norm_errcheckB(a[jt+p3],a[jp+p3],*addr,*addi,icycle[ic_idx],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic_idx, 1, ODD_RADIX, ic_idx);
		}
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

/*...The radix-240 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #if !RAD_15_2FOLD

		tmp = r00;
		tm2 = s1p00;
		for(l = 0; l < 16; l++) {
			// NB: All offsets end up being shifted "one too many" bits leftward since
			// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:

		#if 1	// HLL/ASM toggle

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
		//	printf("15-DFT #%2u: [k0-E]/2 = %u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u\n",l,k0/2,k1/2,k2/2,k3/2,k4/2,k5/2,k6/2,k7/2,k8/2,k9/2,ka/2,kb/2,kc/2,kd/2,ke/2);
			// Input ptrs:		// Output ptrs:
			vc0 = tm2 + k0;		va0 = tmp     ;
			vc1 = tm2 + k1;		va1 = tmp+0x02;
			vc2 = tm2 + k2;		va2 = tmp+0x04;
			vc3 = tm2 + k3;		va3 = tmp+0x06;
			vc4 = tm2 + k4;		va4 = tmp+0x08;
			vc5 = tm2 + k5;		va5 = tmp+0x0a;
			vc6 = tm2 + k6;		va6 = tmp+0x0c;
			vc7 = tm2 + k7;		va7 = tmp+0x0e;
			vc8 = tm2 + k8;		va8 = tmp+0x10;
			vc9 = tm2 + k9;		va9 = tmp+0x12;
			vca = tm2 + ka;		vaa = tmp+0x14;
			vcb = tm2 + kb;		vab = tmp+0x16;
			vcc = tm2 + kc;		vac = tmp+0x18;
			vcd = tm2 + kd;		vad = tmp+0x1a;
			vce = tm2 + ke;		vae = tmp+0x1c;

		#else

		  #if 0	//def USE_AVX512***************************

			// Same table as below, but in uint64 form, and the 16 bytes within the two elements
			// of each row corresponding to 0x[k7,k6,k5,k4,k3,k2,k1,k0], 0x[0,ke,kd,kc,kb,ka,k8]:
			const uint64 dif15_offs[32] = {
				0x8090A0B0C0D0E000ull,0x0010203040506070ull,
				0x718191A1B1C1D1E1ull,0x0001112131415161ull,
				0x62728292A2B2C2D2ull,0x00E2021222324252ull,
				0x5363738393A3B3C3ull,0x00D3E30313233343ull,
				0x445464748494A4B4ull,0x00C4D4E404142434ull,
				0x35455565758595A5ull,0x00B5C5D5E5051525ull,
				0x2636465666768696ull,0x00A6B6C6D6E60616ull,
				0x1727374757677787ull,0x0097A7B7C7D7E707ull,
				0x0818283848586878ull,0x008898A8B8C8D8E8ull,
				0xE909192939495969ull,0x00798999A9B9C9D9ull,
				0xDAEA0A1A2A3A4A5Aull,0x006A7A8A9AAABACAull,
				0xCBDBEB0B1B2B3B4Bull,0x005B6B7B8B9BABBBull,
				0xBCCCDCEC0C1C2C3Cull,0x004C5C6C7C8C9CACull,
				0xADBDCDDDED0D1D2Dull,0x003D4D5D6D7D8D9Dull,
				0x9EAEBECEDEEE0E1Eull,0x002E3E4E5E6E7E8Eull,
				0x8F9FAFBFCFDFEF0Full,0x001F2F3F4F5F6F7Full,
			};

		  #else

			// Same as dit15_offs, but with rows 2-16 reverse-ordered:
			const uint8 dif15_offs[RADIX] = {
				0,224,208,192,176,160,144,128,112,96,80,64,48,32,16,
				225,209,193,177,161,145,129,113,97,81,65,49,33,17,1,
				210,194,178,162,146,130,114,98,82,66,50,34,18,2,226,
				195,179,163,147,131,115,99,83,67,51,35,19,3,227,211,
				180,164,148,132,116,100,84,68,52,36,20,4,228,212,196,
				165,149,133,117,101,85,69,53,37,21,5,229,213,197,181,
				150,134,118,102,86,70,54,38,22,6,230,214,198,182,166,
				135,119,103,87,71,55,39,23,7,231,215,199,183,167,151,
				120,104,88,72,56,40,24,8,232,216,200,184,168,152,136,
				105,89,73,57,41,25,9,233,217,201,185,169,153,137,121,
				90,74,58,42,26,10,234,218,202,186,170,154,138,122,106,
				75,59,43,27,11,235,219,203,187,171,155,139,123,107,91,
				60,44,28,12,236,220,204,188,172,156,140,124,108,92,76,
				45,29,13,237,221,205,189,173,157,141,125,109,93,77,61,
				30,14,238,222,206,190,174,158,142,126,110,94,78,62,46,
				15,239,223,207,191,175,159,143,127,111,95,79,63,47,31
			};
			i = 15*l;
			k0 = dif15_offs[i+ 0];k1 = dif15_offs[i+ 1];k2 = dif15_offs[i+ 2];k3 = dif15_offs[i+ 3];k4 = dif15_offs[i+ 4];
			k5 = dif15_offs[i+ 5];k6 = dif15_offs[i+ 6];k7 = dif15_offs[i+ 7];k8 = dif15_offs[i+ 8];k9 = dif15_offs[i+ 9];
			ka = dif15_offs[i+10];kb = dif15_offs[i+11];kc = dif15_offs[i+12];kd = dif15_offs[i+13];ke = dif15_offs[i+14];
			// Input ptrs:			// Output ptrs:
			vc0 = tm2 + (k0<<1);	va0 = tmp     ;
			vc1 = tm2 + (k1<<1);	va1 = tmp+0x02;
			vc2 = tm2 + (k2<<1);	va2 = tmp+0x04;
			vc3 = tm2 + (k3<<1);	va3 = tmp+0x06;
			vc4 = tm2 + (k4<<1);	va4 = tmp+0x08;
			vc5 = tm2 + (k5<<1);	va5 = tmp+0x0a;
			vc6 = tm2 + (k6<<1);	va6 = tmp+0x0c;
			vc7 = tm2 + (k7<<1);	va7 = tmp+0x0e;
			vc8 = tm2 + (k8<<1);	va8 = tmp+0x10;
			vc9 = tm2 + (k9<<1);	va9 = tmp+0x12;
			vca = tm2 + (ka<<1);	vaa = tmp+0x14;
			vcb = tm2 + (kb<<1);	vab = tmp+0x16;
			vcc = tm2 + (kc<<1);	vac = tmp+0x18;
			vcd = tm2 + (kd<<1);	vad = tmp+0x1a;
			vce = tm2 + (ke<<1);	vae = tmp+0x1c;

		  #endif	// USE_AVX512 ?

		#endif	// HLL/ASM toggle

			SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,
/* s1p00 +... */	vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
/* r00 +... */	va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae
			);
			tmp += 0x1e;	// advance ptr by 15 vec_cmplx [= 30 vec_dbl] elements
		}

	  #else	// RAD_15_2FOLD == True:

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

			SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,two,
/* s1p00 +... */	vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,
/* r00 +... */	va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,
/* s1p00+jt +... */	vd0,vd1,vd2,vd3,vd4,vd5,vd6,vd7,vd8,vd9,vda,vdb,vdc,vdd,vde,
				y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,
/* r0f +... */	vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,vbb,vbc,vbd,vbe
			);
			tmp += 0x3c;	// advance ptr by 30 vec_cmplx [= 60 vec_dbl] elements
		}

	 #endif	// RAD_15_2FOLD ?

		/*...and now do 15 radix-16 transforms, with index-perms as described in radix240_dif_pass1().
		inputs of SSE2_RADIX16_DIT_0TWIDDLE from 30*vec_dbl - separated memlocs, same offsets as already set for DIT: */

		tmp = r00;
		for(l = 0; l < 15; l++) {
			i64 = dif_perm16[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			po_kperm[0x0] = plo[(i64 >> 60)&0xf];
			po_kperm[0x1] = plo[(i64 >> 56)&0xf];
			po_kperm[0x2] = plo[(i64 >> 52)&0xf];
			po_kperm[0x3] = plo[(i64 >> 48)&0xf];
			po_kperm[0x4] = plo[(i64 >> 44)&0xf];
			po_kperm[0x5] = plo[(i64 >> 40)&0xf];
			po_kperm[0x6] = plo[(i64 >> 36)&0xf];
			po_kperm[0x7] = plo[(i64 >> 32)&0xf];
			po_kperm[0x8] = plo[(i64 >> 28)&0xf];
			po_kperm[0x9] = plo[(i64 >> 24)&0xf];
			po_kperm[0xa] = plo[(i64 >> 20)&0xf];
			po_kperm[0xb] = plo[(i64 >> 16)&0xf];
			po_kperm[0xc] = plo[(i64 >> 12)&0xf];
			po_kperm[0xd] = plo[(i64 >>  8)&0xf];
			po_kperm[0xe] = plo[(i64 >>  4)&0xf];
			po_kperm[0xf] = plo[(i64      )&0xf];
			addr = &a[j1] + p_in_hi[l];	// p_in_hi[] = p0,p20,...,p30
			SSE2_RADIX16_DIF_0TWIDDLE(
				tmp,OFF1,OFF2,OFF3,OFF4,
				isrt2,two,
				addr,po_ptr
			);	tmp += 0x2;
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

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

	#endif	// USE_SSE2 ?
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(int k=1; k <= khi; k++) */
