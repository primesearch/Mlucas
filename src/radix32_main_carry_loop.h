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
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-32 DIT pass is here:	*/

	#ifdef USE_SSE2

		add0 = &a[j1    ];	itmp = (int *)arr_offsets;
		SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, r00, isrt2);

	#else	// USE_SSE2 = False:

		RADIX_32_DIT(\
			a+j1,arr_offsets,RE_IM_STRIDE,\
			a+j1,arr_offsets,RE_IM_STRIDE	/* In-place DFT here, so no need for separate input/output index offsets */\
		);

	#endif	// USE_SSE2 ?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to [radix] separate blocks of the A-array, we need [radix] separate carries.	*/

/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

	if(MODULUS_TYPE == MODULUS_TYPE_GENFFTMUL)
	{
		// TODO: 'addi' use is commented out since 'genfftmul_carry_norm_pow2_errcheck' doesn't use its cy parameter in its current implementation
		l = 0; addr = (double *)cy_r; //addi = (double *)cy_i;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
							// Indices in rightmost col are debug-usage only: vvv
			genfftmul_carry_norm_pow2_errcheck(a[jt    ],a[jp    ],*addr,*addi,l); ++l; ++addr; //++addi;
			genfftmul_carry_norm_pow2_errcheck(a[jt+p01],a[jp+p01],*addr,*addi,l); ++l; ++addr; //++addi;
			genfftmul_carry_norm_pow2_errcheck(a[jt+p02],a[jp+p02],*addr,*addi,l); ++l; ++addr; //++addi;
			genfftmul_carry_norm_pow2_errcheck(a[jt+p03],a[jp+p03],*addr,*addi,l); ++l; ++addr; //++addi;
		}
	// continue;	// For gen-FFT case, only expect one mul at a time, rather than primality-test-style
					// one-mul-after-another, thus skip the ensuing first-passof-fFFT stuff
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		// Check if current index-interval contains the target index for rotated-residue carry injection.
		// In data-init we set target_idx = -1 on wraparound-carry mini-pass, so if() only taken on full pass:
		if(target_idx == j) {
		#ifdef USE_SSE2
			addr = (double *)r00 + target_set;
			*addr += target_cy*(n>>1);	// target_cy = [-2 << within-word-shift]*[DWT weight]*n/2, i.e. includes fwd DWT weight and n/2 factor
		#else
			// target_set in [0,2*RADIX); tidx_mod_stride [even|odd] means shifted-carry goes into [Re|Im] part of the complex FFT datum:
			l = target_set&1;	target_set >>= 1;
			a[j1+poff[target_set>>2]+p0123[target_set&3]+l] += target_cy*(n>>1);
		#endif
			target_idx = -1;
		}

	#ifdef USE_AVX
		// For AVX512-and-beyond we support only the fast Mers-carry macros.
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

	 #ifdef USE_AVX512
	  if(1) {	// No HIACC mode for AVX-512
	 #else
	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining
	 #endif
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
		i = (!j);
		addr = &prp_mult;

	  #ifdef USE_AVX512

	   #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
		AVX_cmplx_carry_fast_pow2_wtsinit_X16(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_nm1)
	   #else
		const uint32 nloop = RADIX>>3;
		AVX_cmplx_carry_fast_pow2_wtsinit_X8 (add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_nm1)
	   #endif
		tmp = r00; tm0 = cy_r; itmp = bjmodn;
		for(l = 0; l < nloop; l++) {
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[2*l]);	// poff[] = p0,4,8,...
		  #ifdef CARRY_16_WAY
			AVX_cmplx_carry_fast_pow2_errcheck_X16(tmp, tm0,     itmp, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, tm2,p01,p02,p03,p04, addr);
			tmp += 32; tm0 += 2;        itmp += 16; i = 0;
		  #else
			AVX_cmplx_carry_fast_pow2_errcheck_X8 (tmp, tm0,     itmp, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, tm2,p01,p02,p03,p04, addr);
			tmp += 16; tm0 += 1;         itmp += 8; i = 0;
		  #endif
		}

	  #elif (defined(USE_AVX2) && defined(GCC_5PLUS))

		// gcc 4.x !support the needed AVX2 8-way int instructions (while still being fine for the floating-FMA used for the FFT),
		// so require an added compile-time define to enable this version of 8-way. Based on my [Has,Broad]well timings this is no
		// faster than the version using half-width 128-bit arithmetic for the integer math attendant to Mersenne-mod IBDWT, so
		// leave these calls around just in this one radix-16 carry routine which I used to do my comparative timing tests:
		AVX_cmplx_carry_fast_pow2_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_nm1)
		tmp = r00; tm0 = cy_r; tm1 = cy_r+1; itmp = bjmodn;
		for(l = 0; l < RADIX>>3; l++) {
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[2*l]);	// poff[] = p0,4,8,...
			AVX2_cmplx_carry_fast_pow2_errcheck_X8(tmp, tm0,tm1, itmp, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, tm2,p01,p02,p03,p04, addr);
			tmp += 16; tm0+=2; tm1+=2; itmp += 8; i = 0;
		}

	  #else	// Version targeting AVX-only and AVX2-under-gcc 4.x, which has FMA3 support, but not AVX2 8-way-int support

		AVX_cmplx_carry_fast_pow2_wtsinit_X8(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_nm1)
		tmp = r00; tm0 = cy_r; tm1 = cy_r+1; itmp = bjmodn; itm2 = itmp+4;
		for(l = 0; l < RADIX>>3; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[2*l]);	// poff[] = p0,4,8,...
			AVX_cmplx_carry_fast_pow2_errcheck_X8(tmp, tm0,tm1, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw,tm2,p01,p02,p03,p04, addr);
			tmp += 16; tm0+=2; tm1+=2; itmp += 8; itm2 += 8; i = 0;
		}

	  #endif	// ifdef GCC_5PLUS ?

	  } else {	// HiACC:

	/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		i = (!j);
		addr = &prp_mult;
		tm1 = r00; tmp = cy_r; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[l]);	// poff[] = p0,4,8,...
			AVX_cmplx_carry_norm_pow2_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, tm2,p01,p02,p03, addr);
			tm1 += 8; tmp++; itmp += 4; i = 0;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  }	// USE_AVX: (8-way or 4-way LOACC) or (4-way HIACC) ?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	 #ifdef USE_ARM_V8_SIMD
	  if(1) {	// No HIACC mode for ARMv8
	 #else
	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining
	 #endif
		/*** wt_re,wi_re,wt_im,wi_im inits. Cf. radix16_main_carry_loop.h for scalar-macro prototyping of this: ***/
		uint32 k0,k1,k2,k3;
		l = j & (nwt-1);	nwtml = nwt-l;
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwtml  ];
		sinwtm1 = si[nwtml-1];
		wtl     = wt0[    l  ];
		wtn     = wt0[nwtml  ]*scale;
		wtlp1   = wt0[    l+1];
		wtnm1   = wt0[nwtml-1]*scale;

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

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
		// *But*: since the init macro does an on-the-fly version of this between j,j+2 portions, external code co2=co3 must come *after* both ctmp-data octets are inited.
		SSE2_cmplx_carry_fast_pow2_wtsinit(add1,add2,add3, bjmodn, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, k0,k1,k2,k3, sse_bw,sse_nm1)

		i = (!j);
		addr = &prp_mult;
		tm1 = r00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			SSE2_cmplx_carry_fast_pow2_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p01,p02,p03, addr);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
		}

	  } else {	// HiACC:

		l= j & (nwt-1);	nwtml = nwt-l;	/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwtml  ];			/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwtml-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwtml  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwtml-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = ctmp->im = wtl;		++ctmp;
		ctmp->re = ctmp->im = wtn;		++ctmp;
		ctmp->re = ctmp->im = wtlp1;	++ctmp;
		ctmp->re = ctmp->im = wtnm1;

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		i = (!j);
		addr = &prp_mult;
		tm1 = r00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_pow2_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p01, addr);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
		}

		l= (j+2) & (nwt-1);	nwtml = nwt-l;	/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwtml  ];			/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwtml-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwtml  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwtml-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = ctmp->im = wtl;		++ctmp;
		ctmp->re = ctmp->im = wtn;		++ctmp;
		ctmp->re = ctmp->im = wtlp1;	++ctmp;
		ctmp->re = ctmp->im = wtnm1;

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = r00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_pow2_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p02,p03, addr);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#else	// Scalar-double mode:

		l= j & (nwt-1);	nwtml = nwt-l;	/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwtml  ];			/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwtml-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwtml  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwtml-1]*scale;	/* ...and here.	*/

	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; addr = cy_r; itmp = bjmodn;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			// Re-init weights every 4th macro invocation to keep errors under control:
			cmplx_carry_norm_pow2_errcheck0(a[jt    ],a[jp    ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p01],a[jp+p01],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p02],a[jp+p02],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p03],a[jp+p03],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		}

	  } else {	// HiACC:

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; addr = cy_r; itmp = bjmodn;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			cmplx_carry_norm_pow2_errcheck0(a[jt    ],a[jp    ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p01],a[jp+p01],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p02],a[jp+p02],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p03],a[jp+p03],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		}

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?
	}
	else	/* Fermat-mod carry in SIMD mode */
	{
	#ifdef USE_SSE2
		add0 = &prp_mult;
	#endif
	#ifdef USE_AVX512
		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.
		tmp = half_arr+2;	// In Fermat-mod mode, use first 2 vector slots above half_arr for base and 1/base
		VEC_DBL_INIT(tmp, scale);
		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		tmp = base_negacyclic_root;	tm2 = tmp+1;

		// Hi-accuracy version needs 2 copies of each base root, one for each invocation of the SSE2_fermat_carry_norm_pow2 carry macri:
		l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 16) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 16) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 16) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 16) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 16) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 16) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 16) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 16) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}

		// AVX-custom 8-way carry macro - each contains 8 of the [RADIX] stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6,j+8,j+10,j+12,j+14:
		tm0 = r00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x1000;
		for(i = 0; i < RADIX>>3; i++) {	// RADIX/8 loop passes
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			addr = a + j1 + pfetch_dist + poff[2*i];	// Can't use tm2 as prefetch base-ptr for pow2 Fermat-mod carry-macro since that's used for cy_i
			SSE2_fermat_carry_norm_pow2_errcheck_X8(tm0,tmp,l,tm1,tm2,half_arr,sign_mask, addr,p01,p02,p03,p04, add0);
			tm0 += 16; tm1++; tm2++; tmp += 16; l -= 0x380;
		}

	#elif defined(USE_AVX)

		// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.
		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

		tmp = half_arr+2;
		VEC_DBL_INIT(tmp, scale);
		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		tmp = base_negacyclic_root;	tm2 = tmp+1;

		// Hi-accuracy version needs 8 copies of each base root:
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

		// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:

	// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl) = 2^12
	// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
	// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
	// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
		tm0 = r00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x800;
		for(i = 0; i < RADIX>>2; i++) {	// RADIX/4 loop passes
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			addr = a + j1 + pfetch_dist + poff[i];	// Can't use tm2 as prefetch base-ptr for pow2 Fermat-mod carry-macro since that's used for cy_i
			SSE2_fermat_carry_norm_pow2_errcheck_X4(tm0,tmp,l,tm1,tm2,half_arr,sign_mask, addr,p01,p02,p03, add0);
			tm0 += 8; tm1++; tm2++; tmp += 8; l -= 0xc0;
		}

	#elif defined(USE_SSE2)

		// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.
		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

		tmp = half_arr+2;
		VEC_DBL_INIT(tmp, scale);
		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

		tm1 = r00; tmp = cy_r;	/* <*** Again rely on contiguity of cy_r,i here ***/
		for(l = 0; l < RADIX>>1; l++) {	// RADIX/2 loop passes
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[l>>1]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
			tm2 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_fermat_carry_norm_pow2_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, tm2,p01, add0);
			tm1 += 4; tmp += 2;
		}

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2);
		ntmp = 0; addr = cy_r; addi = cy_i;
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];	// poff[] = p04,08,...
			fermat_carry_norm_pow2_errcheck(a[jt    ],a[jp    ],*addr,*addi,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi;
			fermat_carry_norm_pow2_errcheck(a[jt+p01],a[jp+p01],*addr,*addi,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi;
			fermat_carry_norm_pow2_errcheck(a[jt+p02],a[jp+p02],*addr,*addi,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi;
			fermat_carry_norm_pow2_errcheck(a[jt+p03],a[jp+p03],*addr,*addi,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi;
		}

	#endif	/* #ifdef USE_SSE2 */

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-32 DIF pass is here:	*/

	#ifdef USE_SSE2

		add0 = &a[j1    ];	itmp = (int *)arr_offsets;
		SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, r00, isrt2);

	#else	/* #ifdef USE_SSE2 */

		RADIX_32_DIF(\
			a+j1,arr_offsets,RE_IM_STRIDE,\
			a+j1,arr_offsets,RE_IM_STRIDE	/* In-place DFT here, so no need for separate input/output index offsets */\
		);

	#endif	/* #ifdef USE_SSE2 */

	}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(int k=1; k <= khi; k++) */
