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
	for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-1008 DIT pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (1008 64-bit complex, i.e. 2016 64-bit reals) and do 63 radix-16 transforms...
	  // Radix-16 DFT local-array basic strides OFF1-4 = [1-4] * (2*sizeof(vec_dbl)) * 63 = [1-4] * [sse2:0x20, avx:0x40] * 0x3f:
	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4;
		OFF1 = 0x07e0;
		OFF2 = 0x0fc0;
		OFF3 = 0x17a0;
		OFF4 = 0x1f80;
	  #elif defined(USE_AVX512)
		#define OFF1	0x07e0*4
		#define OFF2	0x0fc0*4
		#define OFF3	0x17a0*4
		#define OFF4	0x1f80*4
	  #elif defined(USE_AVX)
		#define OFF1	0x07e0*2
		#define OFF2	0x0fc0*2
		#define OFF3	0x17a0*2
		#define OFF4	0x1f80*2
	  #else
		#define OFF1	0x07e0
		#define OFF2	0x0fc0
		#define OFF3	0x17a0
		#define OFF4	0x1f80
	  #endif

		tmp = r00;
		for(l = 0; l < ODD_RADIX; ++l) {
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
			addr = &a[j1] + phi[dit_phi[l]];
			SSE2_RADIX16_DIT_0TWIDDLE(
				addr,po_ptr,
				isrt2,two,
				tmp,OFF1,OFF2,OFF3,OFF4	/*** Outs 1,2,6,7,9,14,15 bad! ***/
			);	tmp += 2;
		}

	//...and now do 16 radix-63 transforms:
		tmp = r00;
		for(kk = 0; kk < 16; ++kk)
		{
			k0 = (16-kk) & (-(kk > 0));	// Common index into plo[] and circ-shift count
			iptr = dit_p10_cperms + k0;
			for(l = 0; l < ODD_RADIX; l++)
			{
				ioff[l] = (int)iptr[l]<<5;	// SIMD: Local-mem analog of phi[l] ==> l*16, cast-to-vec-dbl needs another 2x
			}
			SSE2_RADIX_63_DIT( FALSE, thr_id,
				tmp, toff,
				s1p00 + (k0<<1), ioff
			);
			tmp += (ODD_RADIX<<1);
		}

	#else	/* !USE_SSE2 */

	//...gather the needed data (1008 64-bit complex, i.e. 2016 64-bit reals) and do 63 radix-16 transforms...
		tptr = t;
		for(l = 0; l < ODD_RADIX; ++l) {
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
			jp = phi[dit_phi[l]];
			jt = j1 + jp; jp += j2;
			RADIX_16_DIT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				tptr->re,tptr->im,(tptr+0x03f)->re,(tptr+0x03f)->im,(tptr+0x07e)->re,(tptr+0x07e)->im,(tptr+0x0bd)->re,(tptr+0x0bd)->im,(tptr+0x0fc)->re,(tptr+0x0fc)->im,(tptr+0x13b)->re,(tptr+0x13b)->im,(tptr+0x17a)->re,(tptr+0x17a)->im,(tptr+0x1b9)->re,(tptr+0x1b9)->im,(tptr+0x1f8)->re,(tptr+0x1f8)->im,(tptr+0x237)->re,(tptr+0x237)->im,(tptr+0x276)->re,(tptr+0x276)->im,(tptr+0x2b5)->re,(tptr+0x2b5)->im,(tptr+0x2f4)->re,(tptr+0x2f4)->im,(tptr+0x333)->re,(tptr+0x333)->im,(tptr+0x372)->re,(tptr+0x372)->im,(tptr+0x3b1)->re,(tptr+0x3b1)->im,
				c16,s16);	tptr++;
		}
		//...and now do 16 radix-63 transforms:
		tptr = t;
		for(kk = 0; kk < 16; ++kk)
		{
			k0 = (16-kk) & (-(kk > 0));	// Common index into plo[] and circ-shift count
			iptr = dit_p10_cperms + k0;
			for(l = 0; l < ODD_RADIX; l++)
			{
				ioff[l] = phi[iptr[l]];
			}
			RADIX_63_DIT(
				(double *)tptr, toff, 1,
				a+j1+plo[k0], ioff, RE_IM_STRIDE
			);
			tptr += ODD_RADIX;
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

	#ifdef USE_AVX512

		#warning No avx-512 mers-mod carry support yet!

	#elif defined(USE_AVX)

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		l= j & (nwt-1);						tmp = half_arr + 128;	/* ptr to local storage for the doubled wtl,wtn terms: */
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

	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining

		uint32 ii,incr,loop,nloop = RADIX>>3, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy_r; tm2 = cy_r+1; itmp = bjmodn; itm2 = bjmodn+4;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every 8th pass:
		incr = 4;	// incr must divide RADIX/8!
		for(loop = 0; loop < nloop; loop += incr)
		{
			co2 = co2save;	// Need this for all wts-inits beynd the initial set, due to the co2 = co3 preceding the (j+2) data
			ii = loop << 3;	// Reflects 8 independent carry chains being done in each AVX_cmplx_carry_fast_pow2_errcheck_X8 call
			add1 = &wt1[col  +ii];	/* Don't use add0 here, to avoid need to reload main-array address */
			add2 = &wt1[co2-1-ii];
			add3 = &wt1[co3-1-ii];

			// Since use wt1-array in the wtsinit macro, need to fiddle this here:
			co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
			// *But*: since the init macro does an on-the-fly version of this between j,j+2 portions, external code co2=co3 must come *after* both ctmp-data octets are inited.
		  #ifdef USE_AVX512
			ASSERT(0, "AVX-512 version of AVX_cmplx_carry_fast_wtsinit_X8 not yet ported!");
		  #endif
			AVX_cmplx_carry_fast_wtsinit_X8(add1,add2,add3, itmp, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)

			for(l = loop; l < loop+incr; l++) {
				// Each AVX carry macro call also processes 8 prefetches of main-array data
				add0 = a + j1 + pfetch_dist + poff[l+l];
				AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 16; tm1 += 2; tm2 += 2; itmp += 8; itm2 += 8; i = 0;
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

	 #ifdef USE_ARM_V8_SIMD
	  if(1) {	// No HIACC mode for ARMv8
	 #else
	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining
	 #endif
		uint32 i0,i1,i2,i3, ii,incr,nwtml, loop,nloop = RADIX>>2, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every few passes:
		incr = 4;	// incr must divide radix/4!
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

		/* The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
		which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls: */

		tm0 = s1p00; tmp = base_negacyclic_root; l = 0xfc00;
		tm1 = cy_r; // tm2 = cy_i;	*** replace with literal-byte-offset in macro call to save a reg
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic_idx = 0; jc_idx = 1; kc_idx = 2; lc_idx = 3;
		while(tm0 < two)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic_idx];	k5 = jcycle[ic_idx];	k6 = kcycle[ic_idx];	k7 = lcycle[ic_idx];
			k2 = icycle[jc_idx];
			k3 = icycle[kc_idx];
			k4 = icycle[lc_idx];
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x1f80, 0x7e0,0xfc0,0x17a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p1,p2,p3, addr);
			tm0 += 8; tm1++; tmp += 8; l -= 0xc0;
			MOD_ADD32(ic_idx, 4, ODD_RADIX, ic_idx);
			MOD_ADD32(jc_idx, 4, ODD_RADIX, jc_idx);
			MOD_ADD32(kc_idx, 4, ODD_RADIX, kc_idx);
			MOD_ADD32(lc_idx, 4, ODD_RADIX, lc_idx);
		}

	  #else	// HIACC = false:

		// Oct 2014: Try getting most of the LOACC speedup with better accuracy by breaking the complex-roots-of-(-1)
		// chaining into 2 or more equal-sized subchains, each starting with 'fresh' (unchained) complex roots:
		#if (LOACC == 0)
			#warning LOACC = 0
			#define NFOLD 0
		#elif (LOACC == 1)
			#warning LOACC = 1
			#define NFOLD 1
		#elif (LOACC == 2)
			#warning LOACC = 2
			#define NFOLD 2
		#elif (LOACC == 3)
			#warning LOACC = 3
			#define NFOLD 3
		#elif (LOACC == 4)
			#warning LOACC = 4
			#define NFOLD 4
		#elif (LOACC == 5)
			#warning LOACC = 5
			#define NFOLD 5
		#else
			#error If LOACC defined for build of radix1008_ditN_cy_dif1.c, must be given value 0,1,2,3,4 or 5!
		#endif

		#ifdef USE_AVX512
		// For NFOLD > 3,  RADIX not divisible by 2^(3+NFOLD), so use a more-general inner-loop scheme which can handle that:
		  #if NFOLD == 0
			const int nexec[] = {126};
		  #elif NFOLD == 1
			const int nexec[] = {63,63};
		  #elif NFOLD == 2
			const int nexec[] = {32,31,32,31};
		  #elif NFOLD == 3
			const int nexec[] = {16,16,16,15,16,16,16,15};
		  #elif NFOLD == 4
			const int nexec[] = {8,8,8,8,8,8,8,7,8,8,8,8,8,8,8,7};
		  #elif NFOLD == 5
			const int nexec[] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3};
		  #else
			#error NFOLD may only range from 0-5!
		  #endif
		#elif defined(USE_AVX)
		// For NFOLD > 3,  RADIX not divisible by 2^(3+NFOLD), so use a more-general inner-loop scheme which can handle that:
		  #if NFOLD == 0
			const int nexec[] = {252};
		  #elif NFOLD == 1
			const int nexec[] = {126,126};
		  #elif NFOLD == 2
			const int nexec[] = {63,63,63,63};
		  #elif NFOLD == 3
			const int nexec[] = {32,31,32,31,32,31,32,31};
		  #elif NFOLD == 4
			const int nexec[] = {16,16,16,15,16,16,16,15,16,16,16,15,16,16,16,15};
		  #elif NFOLD == 5
			const int nexec[] = {8,8,8,8,8,8,8,7,8,8,8,8,8,8,8,7,8,8,8,8,8,8,8,7,8,8,8,8,8,8,8,7};
		  #else
			#error NFOLD may only range from 0-5!
		  #endif
		#endif

		tm0 = s1p00; tm1 = cy_r; // tm2 = cy_i;	*** replace with literal-byte-offset in macro call to save a reg
		ic_idx = 0; jc_idx = 1; kc_idx = 2; lc_idx = 3;
	  #ifdef USE_AVX512
		mc_idx = 4; nc_idx = 5; oc_idx = 6; pc_idx = 7;
	  #endif

		uint32 naccum = 0;	// Stores sum of [0-ntmp]th elements of nexec[]
		for(ntmp = 0; ntmp < (1 << NFOLD); ++ntmp)
		{
			// E.g.: NFOLD = 1 (==> 2^NFOLD = 2-subchains) means L takes its value
			// from (j) at start of 1st inner-loop exec, and from (j + n/2) at start of 2nd:
		//	l = (j + ntmp*(n>>NFOLD)) >> 1;	*** Only works if RADIX divisible by 2^(lg(RE_IM_STRIDE)+NFOLD)
			l = (j + naccum*NDIVR*RE_IM_STRIDE) >> 1;	naccum += nexec[ntmp];

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

			for(l = 0; l < nexec[ntmp]; l++) {
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
													/* (cy_i_cy_r) --vvvvvv  vvvvvvvvvvvvvvvvv--[1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
				SSE2_fermat_carry_norm_errcheck_X8_loacc(tm0,tmp,tm1,0x1f80, 0x3c0,0x780,0xb40, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, tm2,p1,p2,p3,p4, addr);
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

			for(l = 0; l < nexec[ntmp]; l++) {
				k1 = icycle[ic_idx];
				k2 = icycle[jc_idx];	k5 = jcycle[ic_idx];
				k3 = icycle[kc_idx];	k6 = kcycle[ic_idx];
				k4 = icycle[lc_idx];	k7 = lcycle[ic_idx];
				// Each AVX carry macro call also processes 4 prefetches of main-array data
				tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
													/* (cy_i_cy_r) --vvvvvv  vvvvvvvvvvvvvvvvv--[1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
				SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x1f80, 0x7e0,0xfc0,0x17a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7, tm2,p1,p2,p3, addr);
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

/*...The radix-1008 DIF pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (1008 64-bit complex, i.e 2016 64-bit reals) and do 16 radix-63 transforms...
		tmp = r00;
		for(kk = 0; kk < 16; ++kk)
		{
			iptr = dif_p10_cperms+(kk<<2);
			for(l = 0; l < ODD_RADIX; l++)
			{
				ioff[l] = iptr[l]<<5;	// SIMD: Local-mem analog of phi[l] ==> l*16, cast-to-vec-dbl needs another 2x
			}
			SSE2_RADIX_63_DIF( FALSE, thr_id,
				s1p00 + (kk<<1), ioff,
				tmp, toff
			);
			tmp += (ODD_RADIX<<1);
		}
	//...and now do 63 radix-16 transforms:
		tmp = r00;
		for(l = 0; l < ODD_RADIX; ++l) {
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
			addr = &a[j1] + phi[dif_phi[l]];
			SSE2_RADIX16_DIF_0TWIDDLE(
				tmp,OFF1,OFF2,OFF3,OFF4,
				isrt2,two,
				addr,po_ptr
			);
			tmp += 0x2;
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

	#else	/* !USE_SSE2 */

	//...gather the needed data (1008 64-bit complex, i.e 2016 64-bit reals) and do 16 radix-63 transforms...
		tptr = t;
		for(kk = 0; kk < 16; ++kk)
		{
			iptr = dif_p10_cperms+(kk<<2);
			for(l = 0; l < ODD_RADIX; l++)
			{
				ioff[l] = phi[iptr[l]];
			}
			RADIX_63_DIF(
				a+j1+plo[kk], ioff, RE_IM_STRIDE,
				(double *)tptr, toff, 1
			);
			tptr += ODD_RADIX;
		}
	//...and now do 63 radix-16 transforms:
		tptr = t;
		for(l = 0; l < ODD_RADIX; ++l) {
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
			jp = phi[dif_phi[l]];
			jt = j1 + jp; jp += j2;
			RADIX_16_DIF(
				tptr->re,tptr->im,(tptr+0x03f)->re,(tptr+0x03f)->im,(tptr+0x07e)->re,(tptr+0x07e)->im,(tptr+0x0bd)->re,(tptr+0x0bd)->im,(tptr+0x0fc)->re,(tptr+0x0fc)->im,(tptr+0x13b)->re,(tptr+0x13b)->im,(tptr+0x17a)->re,(tptr+0x17a)->im,(tptr+0x1b9)->re,(tptr+0x1b9)->im,(tptr+0x1f8)->re,(tptr+0x1f8)->im,(tptr+0x237)->re,(tptr+0x237)->im,(tptr+0x276)->re,(tptr+0x276)->im,(tptr+0x2b5)->re,(tptr+0x2b5)->im,(tptr+0x2f4)->re,(tptr+0x2f4)->im,(tptr+0x333)->re,(tptr+0x333)->im,(tptr+0x372)->re,(tptr+0x372)->im,(tptr+0x3b1)->re,(tptr+0x3b1)->im,
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
}	/* end for(k=1; k <= khi; k++) */

