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

	/*...The radix-44 DIT pass is here:	*/

	#ifdef USE_SSE2

		// Radix-4 Outputs in SSE2 mode are temps 2*11*sz_vd bytes apart:
	  #ifdef USE_ARM_V8_SIMD
		const uint32 OFF = 0x160;
	  #elif defined(USE_AVX512)
		#define OFF	0x160*4
	  #elif defined(USE_AVX)
		#define OFF	0x160*2
	  #else
		#define OFF	0x160
	  #endif

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
		const uint32 optr_off[RADIX] = {
			 0<<L2_SZ_VD,24<<L2_SZ_VD,48<<L2_SZ_VD,72<<L2_SZ_VD,8<<L2_SZ_VD,32<<L2_SZ_VD,56<<L2_SZ_VD,80<<L2_SZ_VD,16<<L2_SZ_VD,40<<L2_SZ_VD,64<<L2_SZ_VD,
			22<<L2_SZ_VD,46<<L2_SZ_VD,70<<L2_SZ_VD,6<<L2_SZ_VD,30<<L2_SZ_VD,54<<L2_SZ_VD,78<<L2_SZ_VD,14<<L2_SZ_VD,38<<L2_SZ_VD,62<<L2_SZ_VD,86<<L2_SZ_VD,
			44<<L2_SZ_VD,68<<L2_SZ_VD,4<<L2_SZ_VD,28<<L2_SZ_VD,52<<L2_SZ_VD,76<<L2_SZ_VD,12<<L2_SZ_VD,36<<L2_SZ_VD,60<<L2_SZ_VD,84<<L2_SZ_VD,20<<L2_SZ_VD,
			66<<L2_SZ_VD,2<<L2_SZ_VD,26<<L2_SZ_VD,50<<L2_SZ_VD,74<<L2_SZ_VD,10<<L2_SZ_VD,34<<L2_SZ_VD,58<<L2_SZ_VD,82<<L2_SZ_VD,18<<L2_SZ_VD,42<<L2_SZ_VD
		};
		for(l = 0, tmp = r00, ntmp = 0; l < 4; l++, ntmp += 11) {
			// I-ptrs are regular-stride offsets of r00; O-ptrs are offset w.r.to s1p00;
			// the needed pointer-arithmetic shift has been incorporated into both sets of offsets,
			// so cast both base-pointers to (uint64) to avoid need for add-with-one-shifted-addend:
			// In the DIT-context 11-DFT macro invocation, I-offsets are constant-stride and O-offsets permuted:
			ui32_ptr = &(optr_off[ntmp]);
			SSE2_RADIX_11_DFT(
				tmp,dft11_offptr,
				ua0,
				s1p00,ui32_ptr
			);	tmp += 22;
		}

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
	  #if LO_ADD	// cf. masterdefs.h; this gives higher MUL count but improved accuracy
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

	  if(inc_arr[0]) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of incr
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

		i = (!j);
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy; itmp = bjmodn;
	  #ifndef USE_AVX512
		tm2 = cy+1;
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

		// And do a final 4-fold pass for 40-43:
		add0 = a + j1 + pfetch_dist + poff[l+l];
	   #ifdef USE_AVX512
		// AVX-512 mode calls this macro twice, with Call 2ptr-offsets fiddled as described in comments to that version of the macro:
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,0x000, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);
		add0 += p04;				// prefetch of a + [prefetch offset] + p4,5,6,7
		tmp  = (vec_dbl *)((double *)tmp +  4);	// Call 2 will handle the .d4-7 doubles of our 4 input zmm register-sized vector data
		AVX_cmplx_carry_fast_errcheck_X4(tmp, tm1, itmp, half_arr,0x800, sign_mask,sse_bw,sse_n,sse_sw, add0,p01,p02,p03, addr);	// Call 2 wts-data pointers += 0x400
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

	  if(inc_arr[0]) {	// Have no specialized HIACC carry macro in ARM_V8_SIMD, so use 0-or-not-ness of incr
						// in lieu of (USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) is for ARMv8
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		uint32 k0,k1,k2,k3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every few passes:
		incr = inc_arr;	// When RADIX/4 prime, fiddle the loop to allow e.g. 4+3+4 macro calls separated by reseeding steps.
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
					(and only then: for all subsequent blocks it's superfluous)ghis assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

		/*...The radix-44 DIF pass is here:	*/

	#ifdef USE_SSE2

		/* Do 4 radix-11 transforms: */
		// Radix-11 DFT inputs are (cyclic) with vec_dbl-pointer -= 8 (mod 88) between successive outputs
		const uint32 iptr_off[RADIX] = {
			 0<<L2_SZ_VD,80<<L2_SZ_VD,72<<L2_SZ_VD,64<<L2_SZ_VD,56<<L2_SZ_VD,48<<L2_SZ_VD,40<<L2_SZ_VD,32<<L2_SZ_VD,24<<L2_SZ_VD,16<<L2_SZ_VD,8<<L2_SZ_VD,
			66<<L2_SZ_VD,58<<L2_SZ_VD,50<<L2_SZ_VD,42<<L2_SZ_VD,34<<L2_SZ_VD,26<<L2_SZ_VD,18<<L2_SZ_VD,10<<L2_SZ_VD,2<<L2_SZ_VD,82<<L2_SZ_VD,74<<L2_SZ_VD,
			44<<L2_SZ_VD,36<<L2_SZ_VD,28<<L2_SZ_VD,20<<L2_SZ_VD,12<<L2_SZ_VD,4<<L2_SZ_VD,84<<L2_SZ_VD,76<<L2_SZ_VD,68<<L2_SZ_VD,60<<L2_SZ_VD,52<<L2_SZ_VD,
			22<<L2_SZ_VD,14<<L2_SZ_VD,6<<L2_SZ_VD,86<<L2_SZ_VD,78<<L2_SZ_VD,70<<L2_SZ_VD,62<<L2_SZ_VD,54<<L2_SZ_VD,46<<L2_SZ_VD,38<<L2_SZ_VD,30<<L2_SZ_VD
		};
		for(l = 0, tmp = r00, ntmp = 0; l < 4; l++, ntmp += 11) {
			// O-ptrs are regular-stride offsets of r00; I-ptrs are offset w.r.to s1p00;
			// the needed pointer-arithmetic shift has been incorporated into both sets of offsets,
			// so cast both base-pointers to (uint64) to avoid need for add-with-one-shifted-addend:
			ui32_ptr = &(iptr_off[ntmp]);
			// In the DIF-context 11-DFT macro invocation, I-offsets are permuted and O-offsets constant-stride:
			SSE2_RADIX_11_DFT(
				s1p00,ui32_ptr,
				ua0,
				tmp,dft11_offptr
			);	tmp += 22;
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

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF
	  #endif

	#else	/* !USE_SSE2 */

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-11 transforms...*/
	  #if LO_ADD
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
}	/* end for(int k=1; k <= khi; k++) */
