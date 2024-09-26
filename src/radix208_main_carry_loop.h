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

	/*...The radix-208 DIT pass is here:	*/
	#ifdef USE_SSE2

	  #ifdef USE_ARM_V8_SIMD
		const uint32 OFF1 = 0x20;
		const uint32 OFF2 = 0x40;
		const uint32 OFF3 = 0x60;
		const uint32 OFF4 = 0x80;
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
		tmp = r00;
		kk = 0;
		for(l = 0; l < 13; l++) {
			i64 = dit16_iidx_lo[l];
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
			addr = &a[j1] + phi[kk];	// offset = p10*[0,a,7,4,1,b,8,5,2,c,9,6,3], start idx = 0 and decr 3 (mod 13) each loop
			SSE2_RADIX16_DIT_0TWIDDLE(
				addr,po_ptr,
				isrt2,two,
				tmp,OFF1,OFF2,OFF3,OFF4
			);	tmp += 0x20;
			kk -= 3; kk += ((-(kk < 0)) & 13);
		}
	//...and now do 16 radix-13 transforms:
		tmp = r00;
		for(l = 0; l < 16; l++) {
			tm1 = s1p00 + (((16 - l)&0xf)<<1);	// Low-part offset = p0,f,e,...,2,1
			iptr = dit_pcshft2 + dit_ncshft[l];	// dit_ncshft[l] is an index of the starting element of dit_pcshft[] to use for the current 11-DFT's O-addressing
			// I-ptrs are regular-stride offsets of r00; O-ptrs are offset w.r.to s1p00;
			// the needed pointer-arithmetic shift has been incorporated into both sets of offsets,
			// so cast both base-pointers to (uint64) to avoid need for add-with-one-shifted-addend:
			// In the DIT-context 11-DFT macro invocation, I-offsets are constant-stride and O-offsets permuted:
			SSE2_RADIX_13_DFT(
				tmp,dft13_offptr,
				rad13_const,
				tm1,iptr
			);	tmp += 2;
		}

	#else	// USE_SSE2 = False:

	//...gather the needed data (208 64-bit complex) and do 13 radix-16 transforms:
		tptr = t;
		kk = 0;
		for(l = 0; l < 13; l++) {
			i64 = dit16_iidx_lo[l];
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
			jp = phi[kk];	// = p10*[0,a,7,4,1,b,8,5,2,c,9,6,3], start idx = 0 and decr 3 (mod 13) each loop
			jt = j1 + jp; jp += j2;	// phi[] = p0,pc0,pb0,...,p10
			RADIX_16_DIT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
				c16,s16
			);	tptr += 0x10;
			kk -= 3; kk += ((-(kk < 0)) & 13);
		}
	//...and now do 16 radix-13 transforms:
		tptr = t;
		for(l = 0; l < 16; l++) {
			iptr = dit_pcshft + dit_ncshft[l];
			// Hi-part of p-offset indices:
			k0 = phi[*iptr];
			k1 = phi[*(iptr+0x1)];
			k2 = phi[*(iptr+0x2)];
			k3 = phi[*(iptr+0x3)];
			k4 = phi[*(iptr+0x4)];
			k5 = phi[*(iptr+0x5)];
			k6 = phi[*(iptr+0x6)];
			k7 = phi[*(iptr+0x7)];
			k8 = phi[*(iptr+0x8)];
			k9 = phi[*(iptr+0x9)];
			ka = phi[*(iptr+0xa)];
			kb = phi[*(iptr+0xb)];
			kc = phi[*(iptr+0xc)];
			jp = plo[(16 - l)&0xf];	// Low-part offset = p0,f,e,...,2,1
			jt = j1 + jp; jp += j2;
			RADIX_13_DFT(
				tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc]
			);	tptr++;
		}

	#endif	// SIMD or not?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need 28 separate carries.	*/

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
		uint32 ii,loop, co2save = co2;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every 8th pass or better:
	  #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
	  #else
		const uint32 nloop = RADIX>>3;
	  #endif
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy; itmp = bjmodn;	// tm2,itm2 not used in AVX-512 mode
	  #ifndef USE_AVX512
		tm2 = cy+1; itm2 = bjmodn+4;
	  #endif
		incr = inc_arr;
		for(loop = 0; loop < nloop; loop += *incr++)
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
			for(l = loop; l < loop+*incr; l++) {
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
		tm1 = s1p00; tmp = cy; itmp = bjmodn;
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

	  if(inc_arr[0]) {	// Have no specialized HIACC carry macro in ARM_V8_SIMD, so use 0-or-not-ness of incr
						// in lieu of (USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) is for ARMv8
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		// DIF/DIT loops use k1,k2 and there those mustbbe signed, so use i0-3 here:
		uint32 i0,i1,i2,i3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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

			for(l = loop; l < loop+*incr; l++) {
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

		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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

	/*	i =((uint32)(sw - bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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
		l = 0; addr = cy; itmp = bjmodn;
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
		l = 0; addr = cy; itmp = bjmodn;
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
	#endif	/* #ifdef USE_SSE2 */

/*...The radix-208 DIF pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (208 64-bit complex) and do 16 radix-13 transforms:
		tmp = r00;
		for(l = 0; l < 16; l++) {
			tm1 = s1p00 + ((((l<<1)+l) & 0xf)<<1);	// Low-part offset = p[3*l (mod 16)] = p0,3,6,9,c,f,2,5,...
			iptr = dif_pcshft2 + dif_ncshft[l];
			// O-ptrs are regular-stride offsets of r00; I-ptrs are offset w.r.to s1p00;
			// the needed pointer-arithmetic shift has been incorporated into both sets of offsets,
			// so cast both base-pointers to (uint64) to avoid need for add-with-one-shifted-addend:
			// In the DIF-context 11-DFT macro invocation, I-offsets are permuted and O-offsets constant-stride:
			SSE2_RADIX_13_DFT(
				tm1,iptr,
				rad13_const,
				tmp,dft13_offptr
			);	tmp += 2;
		}
	//...and now do 13 radix-16 transforms:
		tmp = r00;
		for(l = 0; l < 13; l++) {
			i64 = dif16_oidx_lo[l];
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
			addr = &a[j1] + phi[dif_pcshft[l]];	// offset = p0,pc0,pb0,...,p10
			SSE2_RADIX16_DIF_0TWIDDLE(
				tmp,OFF1,OFF2,OFF3,OFF4,
				isrt2,two,
				addr,po_ptr
			);	tmp += 0x20;
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

	#else	// USE_SSE2 = False:

	//...gather the needed data (208 64-bit complex) and do 16 radix-13 transforms:
		tptr = t;
		for(l = 0; l < 16; l++) {
			iptr = dif_pcshft + dif_ncshft[l];
			// Hi-part of p-offset indices:
			k0 = phi[*iptr];
			k1 = phi[*(iptr+0x1)];
			k2 = phi[*(iptr+0x2)];
			k3 = phi[*(iptr+0x3)];
			k4 = phi[*(iptr+0x4)];
			k5 = phi[*(iptr+0x5)];
			k6 = phi[*(iptr+0x6)];
			k7 = phi[*(iptr+0x7)];
			k8 = phi[*(iptr+0x8)];
			k9 = phi[*(iptr+0x9)];
			ka = phi[*(iptr+0xa)];
			kb = phi[*(iptr+0xb)];
			kc = phi[*(iptr+0xc)];
			jp = plo[((l<<1)+l) & 0xf];	// Low-part offset = p[3*l (mod 16)] ...
			jt = j1 + jp; jp += j2;		// ... = p0,3,6,9,c,f,2,5,...
			RADIX_13_DFT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],
				tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im
			);	tptr++;
		}
	//...and now do 13 radix-16 transforms:
		tptr = t;
		for(l = 0; l < 13; l++) {
			i64 = dif16_oidx_lo[l];
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
			jp = phi[dif_pcshft[l]];	// = p0,pc0,pb0,...,p10
			jt = j1 + jp; jp += j2;	// dif_phi[] = p0,pc0,pb0,...,p10
			RADIX_16_DIF(
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				c16,s16);	tptr += 0x10;
		}

	#endif	/* !USE_SSE2 */

	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(int k=1; k <= khi; k++) */
