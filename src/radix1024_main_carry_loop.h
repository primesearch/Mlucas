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
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

/*...The radix-1024 DIT pass is here:	*/

	#ifdef USE_SSE2

	// Gather the needed data and do 16 twiddleless length-64 subtransforms, with p-offsets in-order:

		for(l = 0; l < 16; l++) {
			jt = j1 + dit_poffs[l];	// poffs[] = p40,p80,...,p3c0
			SSE2_RADIX_64_DIT( FALSE, thr_id,
				// Inputs: Base address plus index offsets:
				a+jt,dit_i_offsets,
				// Intermediates-storage pointer:
				vd00,
				// Output pointer: Base ptr of 64-vec_cmplx [128 vec_dbl, whence the (l<<7)] local-mem block:
		//		r00 + (l<<7),o_offsets
				(vec_dbl *)(r00 + (l<<7)), o_offsets
			);
		}

	/*...and now do 64 radix-16 subtransforms, including the internal twiddle factors - we use the same positive-power
	roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
	in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f],
	with each of the foregoing 16 indices being head of a (i0,i0+p20,i0+p10,i0+p30) quartet:
	*/
	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4;
		OFF1 = 0x0800;
		OFF2 = 0x1000;
		OFF3 = 0x1800;
		OFF4 = 0x2000;
	  #elif defined(USE_AVX512)
		#define OFF1	4*0x0800
		#define OFF2	4*0x1000
		#define OFF3	4*0x1800
		#define OFF4	4*0x2000
	  #elif defined(USE_AVX)
		#define OFF1	2*0x0800
		#define OFF2	2*0x1000
		#define OFF3	2*0x1800
		#define OFF4	2*0x2000
	  #else
		#define OFF1	0x0800	// Base stride = 64 vec_cpmplx = 128 vec_dbl ==> 0x80<<4 = 0x800 bytes in SSE2 mode
		#define OFF2	0x1000
		#define OFF3	0x1800
		#define OFF4	0x2000
	  #endif

	  #if defined(USE_AVX2) && !defined(USE_IMCI512)

		// Due to tangent-twiddles scheme and resulting singularity of tangent(arg(I)) = 1/0,
		// only last 62 of the 63 with-twiddles DFTs allow use of FMA-based macros under Intel AVX2/FMA3:
		for(l = 0; l < 2; l++) {
			jt = reverse(l,6)<<1;
			tm1 = r00 + jt; tm2 = s1p00 + jt; tmp = twid00 + (jt<<4)-jt;	// Twid-offsets are multiples of 30 vec_dbl
			SSE2_RADIX16_DIT_TWIDDLE_OOP(
				tm1,OFF1,OFF2,OFF3,OFF4, tm2,OFF1,OFF2,OFF3,OFF4, isrt2, tmp
			);
		}
		for(l = 2; l < 64; l++) {
			jt = reverse(l,6)<<1;
			tm1 = r00 + jt; tm2 = s1p00 + jt; tmp = twid00 + (jt<<4)-jt;	// Twid-offsets are multiples of 30 vec_dbl
			SSE2_RADIX16_DIT_FMA_OOP(
				tm1,OFF1,OFF2,OFF3,OFF4, tm2,OFF1,OFF2,OFF3,OFF4, tmp
			);
		}

	  #else	// Non-FMA version:

		// Block 0: has all-unity twiddles, but not worth doing separately here
		// in the "Why add extra code to save (a few %)/64?" sense:
		for(l = 0; l < 64; l++) {
			jt = reverse(l,6)<<1;
			tm1 = r00 + jt; tm2 = s1p00 + jt; tmp = twid00 + (jt<<4)-jt;	// Twid-offsets are multiples of 30 vec_dbl
			SSE2_RADIX16_DIT_TWIDDLE_OOP(
				tm1,OFF1,OFF2,OFF3,OFF4, tm2,OFF1,OFF2,OFF3,OFF4, isrt2, tmp
			);
		}

	  #endif	// FMA/AVX2 ?

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1	// DIF and DIT share same radix-16 DFT strides in this case, so could def
		#undef OFF2	// once at top and undef at bottom, but keep it clean and local-def.
		#undef OFF3
		#undef OFF4
	  #endif

	#else

	// Gather the needed data and do 16 twiddleless length-64 subtransforms, with p-offsets in-order:
		for(l = 0, jp = 0; l < 16; l++, jp += 64) {
			jt = j1 + dit_poffs[l];	// poffs[] = p40,p80,...,p3c0
			RADIX_64_DIT((a+jt),dit_i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		}
	//...and now do 64 radix-16 subtransforms, including the internal twiddle factors
		// Block 0: has all-unity twiddles
		tptr = t;
		jt = j1;	jp = j2;
		ju = jt+p200;	jv = jp+p200;
		// Twiddleless DIF bit-reverses its outputs, so a_p* terms appear in BR-order [swap index pairs 1/8,2/4,3/c,5/a,7/e.b/d]:
		RADIX_16_DIT(
			tptr->re,tptr->im,(tptr+0x040)->re,(tptr+0x040)->im,(tptr+0x080)->re,(tptr+0x080)->im,(tptr+0x0c0)->re,(tptr+0x0c0)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x200)->re,(tptr+0x200)->im,(tptr+0x240)->re,(tptr+0x240)->im,(tptr+0x280)->re,(tptr+0x280)->im,(tptr+0x2c0)->re,(tptr+0x2c0)->im,(tptr+0x300)->re,(tptr+0x300)->im,(tptr+0x340)->re,(tptr+0x340)->im,(tptr+0x380)->re,(tptr+0x380)->im,(tptr+0x3c0)->re,(tptr+0x3c0)->im,
			a[jt     ],a[jp     ],a[jt+ p40],a[jp+ p40],a[jt+ p80],a[jp+ p80],a[jt+ pc0],a[jp+ pc0],a[jt+p100],a[jp+p100],a[jt+p140],a[jp+p140],a[jt+p180],a[jp+p180],a[jt+p1c0],a[jp+p1c0],a[ju     ],a[jv     ],a[ju+ p40],a[jv+ p40],a[ju+ p80],a[jv+ p80],a[ju+ pc0],a[jv+ pc0],a[ju+p100],a[jv+p100],a[ju+p140],a[jv+p140],a[ju+p180],a[jv+p180],a[ju+p1c0],a[jv+p1c0],
			c16,s16
		);

		// Remaining 63 sets of macro calls done in loop:
		for(ntmp = 1; ntmp < 64; ntmp++) {
			tptr = t + reverse(ntmp,6);
			jt = j1 + dit_po_br[ntmp]; jp = j2 + dit_po_br[ntmp];	// po_br[] = p[084c2a6e195d3b7f]
			ju = jt+p200;	jv = jp+p200;
			const double *addr = DFT1024_TWIDDLES[ntmp], *addi = addr+1;	// Pointer to required row of 2-D twiddles array
			RADIX_16_DIT_TWIDDLE_OOP(
				tptr->re,tptr->im,(tptr+0x040)->re,(tptr+0x040)->im,(tptr+0x080)->re,(tptr+0x080)->im,(tptr+0x0c0)->re,(tptr+0x0c0)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x200)->re,(tptr+0x200)->im,(tptr+0x240)->re,(tptr+0x240)->im,(tptr+0x280)->re,(tptr+0x280)->im,(tptr+0x2c0)->re,(tptr+0x2c0)->im,(tptr+0x300)->re,(tptr+0x300)->im,(tptr+0x340)->re,(tptr+0x340)->im,(tptr+0x380)->re,(tptr+0x380)->im,(tptr+0x3c0)->re,(tptr+0x3c0)->im,
				a[jt     ],a[jp     ],a[jt+ p40],a[jp+ p40],a[jt+ p80],a[jp+ p80],a[jt+ pc0],a[jp+ pc0],a[jt+p100],a[jp+p100],a[jt+p140],a[jp+p140],a[jt+p180],a[jp+p180],a[jt+p1c0],a[jp+p1c0],a[ju     ],a[jv     ],a[ju+ p40],a[jv+ p40],a[ju+ p80],a[jv+ p80],a[ju+ pc0],a[jv+ pc0],a[ju+p100],a[jv+p100],a[ju+p140],a[jv+p140],a[ju+p180],a[jv+p180],a[ju+p1c0],a[jv+p1c0],
				*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
				c16,s16
			);
		}

	#endif	// USE_SSE2?

/*...Now do the carries. Since the outputs would
normally be getting dispatched to [radix] separate blocks of the A-array, we need [radix] separate carries.	*/

/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		// Check if current index-interval contains the target index for rotated-residue carry injection.
		// In data-init we set target_idx = -1 on wraparound-carry mini-pass, so if() only taken on full pass:
		if(target_idx == j) {
		#ifdef USE_SSE2
			double *addr_ = (double *)s1p00 + target_set;
			*addr_ += target_cy*(n>>1);	// target_cy = [-2 << within-word-shift]*[DWT weight]*n/2, i.e. includes fwd DWT weight and n/2 factor
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
		uint32 ii,loop, co2save = co2;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every 8th pass or better:
		// incr must divide nloop [RADIX/8 = 128 or RADIX/16 = 64, depending on whether we use 8-or-16-way carry macros]!
	  #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
	  #else
		const uint32 nloop = RADIX>>3;
	  #endif
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy_r; tm2 = cy_r+1; itmp = bjmodn;
	  #ifndef USE_AVX512
		itm2 = bjmodn+4;	// itm2 not used in AVX-512 mode
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
			AVX_cmplx_carry_fast_pow2_wtsinit_X16(add1,add2,add3, itmp, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_nm1)
		  #else
			AVX_cmplx_carry_fast_pow2_wtsinit_X8 (add1,add2,add3, itmp, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_nm1)
		  #endif

			for(l = loop; l < loop+incr; l++) {
				// Each AVX carry macro call also processes 8 prefetches of main-array data
				add0 = a + j1 + pfetch_dist + poff[l+l];
			#ifdef USE_AVX512
			  #ifdef CARRY_16_WAY
				AVX_cmplx_carry_fast_pow2_errcheck_X16(tmp, tm1,     itmp, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 32; tm1 += 2;           itmp += 16; i = 0;
			  #else
				AVX_cmplx_carry_fast_pow2_errcheck_X8 (tmp, tm1,     itmp, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 16; tm1 += 1;           itmp += 8; i = 0;
			  #endif
			#else
				AVX_cmplx_carry_fast_pow2_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr);
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
			add0 = a + j1 + pfetch_dist + poff[l];
			AVX_cmplx_carry_norm_pow2_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);
			tm1 += 8; tmp++; itmp += 4; i = 0;
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
		uint32 k0,k1,k2,k3, ii,loop,nwtml, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every 8th pass:
		for(loop = 0; loop < RADIX>>2; loop += incr)
		{
			ii = loop << 2;	// Reflects 4 independent carry chains being done in eah SSE2_cmplx_carry_fast_pow2_errcheck call
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
			SSE2_cmplx_carry_fast_pow2_wtsinit(add1,add2,add3, add0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, k0,k1,k2,k3, sse_bw,sse_nm1)

			for(l = loop; l < loop+incr; l++) {
				// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data:
				add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
				SSE2_cmplx_carry_fast_pow2_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);
				tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
			}
		}

	  } else {	// HiACC:

		l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = ctmp->im = wtl;		++ctmp;
		ctmp->re = ctmp->im = wtn;		++ctmp;
		ctmp->re = ctmp->im = wtlp1;	++ctmp;
		ctmp->re = ctmp->im = wtnm1;

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		addr = &prp_mult;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_pow2_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1, addr);
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
		ctmp->re = ctmp->im = wtl;		++ctmp;
		ctmp->re = ctmp->im = wtn;		++ctmp;
		ctmp->re = ctmp->im = wtlp1;	++ctmp;
		ctmp->re = ctmp->im = wtnm1;

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_pow2_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p2,p3, addr);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#else	// Scalar-double mode:

		l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; itmp = bjmodn;
		double *addr_ = cy_r;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			// Re-init weights every 4th macro invocation to keep errors under control:
			cmplx_carry_norm_pow2_errcheck0(a[jt   ],a[jp   ],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p1],a[jp+p1],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p2],a[jp+p2],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p3],a[jp+p3],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
		}

	  } else {	// HiACC:

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; itmp = bjmodn;
		double *addr_ = cy_r;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			cmplx_carry_norm_pow2_errcheck0(a[jt   ],a[jp   ],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p1],a[jp+p1],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p2],a[jp+p2],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p3],a[jp+p3],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
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
		tm0 = s1p00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x20000;
		for(i = 0; i < RADIX>>3; i++) {	// RADIX/8 loop passes
			// Each AVX carry macro call also processes 8 prefetches of main-array data
			addr = a + j1 + pfetch_dist + poff[2*i];	// Can't use tm2 as prefetch base-ptr for pow2 Fermat-mod carry-macro since that's used for cy_i
			SSE2_fermat_carry_norm_pow2_errcheck_X8(tm0,tmp,l,tm1,tm2,half_arr,sign_mask, addr,p1,p2,p3,p4, add0);
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

	// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl) = 2^16
	// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
	// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
	// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
		tm0 = s1p00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x10000;
		for(i = 0; i < RADIX>>2; i++) {	// RADIX/4 loop passes
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			addr = a + j1 + pfetch_dist + poff[i];	// Can't use tm2 as prefetch base-ptr for pow2 Fermat-mod carry-macro since that's used for cy_i
			SSE2_fermat_carry_norm_pow2_errcheck_X4(tm0,tmp,l,tm1,tm2,half_arr,sign_mask, addr,p1,p2,p3, add0);
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

		tm1 = s1p00; tmp = cy_r;	/* <*** Again rely on contiguity of cy_r,i here ***/
		for(l = 0; l < RADIX>>1; l++) {	// RADIX/2 loop passes
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			addr = a + j1 + pfetch_dist + poff[l>>1];	// poff[] = p0,4,8,...
			addr += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_fermat_carry_norm_pow2_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p2, add0);
			tm1 += 4; tmp += 2;
		}

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		ntmp = 0;
		double *addr_ = cy_r, *addi_ = cy_i;
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];
			fermat_carry_norm_pow2_errcheck(a[jt   ],a[jp   ],*addr_,*addi_,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr_; ++addi_;
			fermat_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],*addr_,*addi_,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr_; ++addi_;
			fermat_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],*addr_,*addi_,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr_; ++addi_;
			fermat_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],*addr_,*addi_,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr_; ++addi_;
		}

	#endif	/* #ifdef USE_SSE2 */

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-1024 DIF pass is here:	*/

	#ifdef USE_SSE2

	// Gather the needed data and do 16 twiddleless length-64 subtransforms, with p-offsets in br16 order: 084c2a6e195d3b7f:

		for(l = 0; l < 16; l++) {
			jt = reverse(l,4);	// po_br[] = p[084c2a6e195d3b7f], here sans 'p' because SIMD takes ins from a contig-local memblock
			SSE2_RADIX_64_DIF( FALSE, thr_id,
				4,	// set = trailz(N) - trailz(64)
				// Input pointer; no offsets array in pow2-radix case:
				(double *)(s1p00 + (jt<<1)), 0x0,
				// Intermediates-storage pointer:
				vd00,
				// Outputs: Base address plus index offsets:
				(double *)(r00 + (l<<7)), dif_o_offsets
			);
		}

	//...and now do 64 radix-16 subtransforms, including the internal twiddles:

	  #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x0800;
		OFF2 = 0x1000;
		OFF3 = 0x1800;
		OFF4 = 0x2000;
	  #elif defined(USE_AVX512)
		#define OFF1	4*0x0800
		#define OFF2	4*0x1000
		#define OFF3	4*0x1800
		#define OFF4	4*0x2000
	  #elif defined(USE_AVX)
		#define OFF1	2*0x0800
		#define OFF2	2*0x1000
		#define OFF3	2*0x1800
		#define OFF4	2*0x2000
	  #else
		#define OFF1	0x0800	// Base stride = 64 vec_cpmplx = 128 vec_dbl ==> 0x80<<4 = 0x800 bytes in SSE2 mode
		#define OFF2	0x1000
		#define OFF3	0x1800
		#define OFF4	0x2000
	  #endif
		// SSE2_RADIX16_DIF_TWIDDLE_OOP needs these:
		const int off_arr[16] = {0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf};
		const int*off_ptr = &(off_arr[0]);

	  #if defined(USE_AVX2) && !defined(USE_IMCI512)

		// Due to tangent-twiddles scheme and resulting singularity of tangent(arg(I)) = 1/0,
		// only last 62 of the 63 with-twiddles DFTs allow use of FMA-based macros under Intel AVX2/FMA3:
		tm1 = r00;
		for(l = 0; l < 2; l++) {
			ntmp = reverse(l,6)<<1;
			tmp = twid00 + (ntmp<<4)-ntmp;	// Twid-offsets are multiples of 30 vec_dbl
			add0 = &a[j1] + dif_i_offsets[l];	// poffs[] = p10,p20,...,p3f0
			SSE2_RADIX16_DIF_TWIDDLE_OOP(
				tm1,OFF1,OFF4,
				add0,off_ptr,
				isrt2, tmp
			);	tm1 += 2;
		}
		for(l = 2; l < 64; l++) {
			ntmp = reverse(l,6)<<1;
			tmp = twid00 + (ntmp<<4)-ntmp;	// Twid-offsets are multiples of 30 vec_dbl
			add0 = &a[j1] + dif_i_offsets[l];	// poffs[] = p10,p20,...,p3f0
				add1 = add0+p1; add2 = add0+p2; add3 = add0+p3; add4 = add0+p4; add5 = add0+p5; add6 = add0+p6; add7 = add0+p7;
			add8 = add0+p8; add9 = add1+p8; adda = add2+p8; addb = add3+p8; addc = add4+p8; addd = add5+p8; adde = add6+p8; addf = add7+p8;
			SSE2_RADIX16_DIF_FMA_OOP(
				tm1,OFF1,OFF2,OFF3,OFF4,
				add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf,
				tmp
			);	tm1 += 2;
		}

	  #else	// Non-FMA version:

		tm1 = r00;
		for(l = 0; l < 64; l++) {
			ntmp = reverse(l,6)<<1;
			tmp = twid00 + (ntmp<<4)-ntmp;	// Twid-offsets are multiples of 30 vec_dbl
			add0 = &a[j1] + dif_i_offsets[l];	// poffs[] = p10,p20,...,p3f0
			SSE2_RADIX16_DIF_TWIDDLE_OOP(
				tm1,OFF1,OFF4,
				add0,off_ptr,
				isrt2, tmp
			);	tm1 += 2;
		}

	  #endif	// FMA/AVX2 ?

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

	#else	/* !USE_SSE2 */

	// Gather the needed data and do 16 twiddleless length-64 subtransforms, with p-offsets in br16 order: 084c2a6e195d3b7f:
		for(l = 0, jp = 0; l < 16; l++, jp += 64) {
			jt = j1 + dif_po_br[l];	// po_br[] = p[084c2a6e195d3b7f]
			//	NOTE that RADIX_64_DIF outputs are IN-ORDER rather than BR:
			RADIX_64_DIF((a+jt),dif_i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		}
	//...and now do 64 radix-16 subtransforms, including the internal twiddles:
		// Block 0: has all-unity twiddles
		tptr = t;
		jt = j1;	jp = j2;
	/*
		// Twiddleless DIF bit-reverses its outputs, so a_p* terms appear in BR-order [swap index pairs 1/8,2/4,3/c,5/a,7/e.b/d]:
		RADIX_16_DIF(
			tptr->re,tptr->im,(tptr+0x200)->re,(tptr+0x200)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x300)->re,(tptr+0x300)->im,(tptr+0x080)->re,(tptr+0x080)->im,(tptr+0x280)->re,(tptr+0x280)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x380)->re,(tptr+0x380)->im,(tptr+0x040)->re,(tptr+0x040)->im,(tptr+0x240)->re,(tptr+0x240)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x340)->re,(tptr+0x340)->im,(tptr+0x0c0)->re,(tptr+0x0c0)->im,(tptr+0x2c0)->re,(tptr+0x2c0)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x3c0)->re,(tptr+0x3c0)->im,
			a[jt],a[jp],a[jt+p1],a[jp+p1],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],
			c16,s16
		);	tptr++;
		// Remaining 63 sets of macro calls done in loop:
		for(l = 1; l < 64; l++) {
	*/
		for(l = 0; l < 64; l++) {
			jt = j1 + dif_i_offsets[l]; jp = j2 + dif_i_offsets[l];	// poffs[] = p10,p20,...,p3f0
			const double *addr = DFT1024_TWIDDLES[l], *addi = addr+1;	// Pointer to required row of 2-D twiddles array
			RADIX_16_DIF_TWIDDLE_OOP(
				tptr->re,tptr->im,(tptr+0x200)->re,(tptr+0x200)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x300)->re,(tptr+0x300)->im,(tptr+0x080)->re,(tptr+0x080)->im,(tptr+0x280)->re,(tptr+0x280)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x380)->re,(tptr+0x380)->im,(tptr+0x040)->re,(tptr+0x040)->im,(tptr+0x240)->re,(tptr+0x240)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x340)->re,(tptr+0x340)->im,(tptr+0x0c0)->re,(tptr+0x0c0)->im,(tptr+0x2c0)->re,(tptr+0x2c0)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x3c0)->re,(tptr+0x3c0)->im,
				a[jt],a[jp],a[jt+p1],a[jp+p1],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],
				*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
				c16,s16
			);	tptr++;
		}

	#endif	// USE_SSE2?

	}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(int k=1; k <= khi; k++) */

