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

/*...The radix-128 DIT pass is here:	*/

#ifdef USE_SSE2

/* Gather the needed data (128 64-bit complex, i.e. 256 64-bit reals) and do 8 twiddleless length-16 subtransforms: */

	#if 1

		int dit_offsets[RADIX];
		dit_offsets[0x0] =   0; dit_offsets[0x1] = p01; dit_offsets[0x2] = p02; dit_offsets[0x3] = p03; dit_offsets[0x4] = p04; dit_offsets[0x5] = p05; dit_offsets[0x6] = p06; dit_offsets[0x7] = p07;
		dit_offsets[0x8] = p08; dit_offsets[0x9] = p09; dit_offsets[0xa] = p0a; dit_offsets[0xb] = p0b; dit_offsets[0xc] = p0c; dit_offsets[0xd] = p0d; dit_offsets[0xe] = p0e; dit_offsets[0xf] = p0f;
		for(l = 0; l < 16; l++) {
			dit_offsets[l+0x10] = dit_offsets[l] + p10;
			dit_offsets[l+0x20] = dit_offsets[l] + p20;
			dit_offsets[l+0x30] = dit_offsets[l] + p30;
			dit_offsets[l+0x40] = dit_offsets[l] + p40;
			dit_offsets[l+0x50] = dit_offsets[l] + p50;
			dit_offsets[l+0x60] = dit_offsets[l] + p60;
			dit_offsets[l+0x70] = dit_offsets[l] + p70;
		}
		SSE2_RADIX128_DIT(	// Inputs in s1p00[ 00- 7f]
			(a+j1     ),dit_offsets,
			r00, isrt2,two,twid0,
			s1p00
		);

	#else
	  #if !COMPACT_OBJ
		#error As of 2021, COMPACT_OBJ must be def'd nonzero in this routine's .c file!
	  #endif
	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7;
		OFF1 = 0x20;
		OFF2 = 0x40;
		OFF3 = 0x60;
		OFF4 = 0x80;
	  #elif defined(USE_AVX512)
		#define OFF1	4*0x20
		#define OFF2	4*0x40
		#define OFF3	4*0x60
		#define OFF4	4*0x80
	  #elif defined(USE_AVX)
		#define OFF1	2*0x20
		#define OFF2	2*0x40
		#define OFF3	2*0x60
		#define OFF4	2*0x80
	  #else
		#define OFF1	0x20
		#define OFF2	0x40
		#define OFF3	0x60
		#define OFF4	0x80
	  #endif

		// Start of isrt2/cc16/ss16 triplet needed for radix-16 SSE2 DFT macros [use instead of isrt2 here]
	  #if COMPACT_OBJ
		tmp = isrt2;	// In COMPACT_OBJ mode, use that [isrt2] ptr heads up the needed isrt2/cc16/ss16 pointer-triplet
						// (actually a quartet in AVX2 mode, since there need sqrt2 preceding isrt2)
	  #else
		tmp = sqrt2+1;	// In non-compact-obj-code mode, access the same (unnamed) pointer this way
	  #endif
		tm1 = r00;
		for(l = 0; l < 8; l++) {
			addr = &a[j1] + poff[l<<2];	// poff[4*l] = p10,p20,...,pf0
			SSE2_RADIX16_DIT_0TWIDDLE(
				addr,po_lin,
				tmp,two,
				tm1,OFF1,OFF2,OFF3,OFF4
			);
			tm1 += 32;
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

	/*...and now do 16 radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */

	/* Block 0: */
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
			r00,r10,r20,r30,r40,r50,r60,r70,
			s1p00,s1p70,s1p60,s1p50,s1p40,s1p30,s1p20,s1p10, isrt2,two
		);
	  #else
		SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
			r00,r10,r20,r30,r40,r50,r60,r70,
			s1p00,s1p70,s1p60,s1p50,s1p40,s1p30,s1p20,s1p10, isrt2
		);
	  #endif

	/* Note: Each of the 7 sincos-pairs sent to SSE2_RADIX8_DIT_TWIDDLE_OOP is actually a vector-data triplet
	whose first (unnamed) member contains the isrt2 needed for radix-8.
	This is a workaround of GCC's 30-arg limit for inline ASM macros, which proves a royal pain here.
	*/
	/*** v20: Reduced-#args code mods to accommodate Clang builds on Apple M1 moot the 30-arg issue, but
	leave above convention in place for now. Don't need the OFF3-arg in SSE2_RADIX8_DIT_TWIDDLE_OOP(): ***/
	   #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x200;
		OFF2 = 0x400;
		OFF4 = 0x800;
	  #elif defined(USE_AVX512)
		#define OFF1	0x0800
		#define OFF2	0x1000
		#define OFF4	0x2000
	  #elif defined(USE_AVX)
		#define OFF1	0x0400
		#define OFF2	0x0800
		#define OFF4	0x1000
	  #else
		#define OFF1	0x200
		#define OFF2	0x400
		#define OFF4	0x800
	   #endif
		uint64 dit_o_ptr_stride = (uint64)OFF1;
		// Remaining 15 sets of macro calls done in loop:
		for(l = 1; l < 16; l++) {
			k1 = reverse(l,4);
			// Twid-offsets are multiples of 21 vec_dbl; +1 to point to cos [middle] term of each [isrt2,c,s] triplet
			tm1 = r00 + (k1<<1); tm2 = s1p00 + (k1<<1); tmp = twid0 + (k1<<4)+(k1<<2)+k1;
			// 7 [c,s] pointer-pairs, +=3 between pairs; c0 = tmp+1 to point to cos [middle] term of each [isrt2,c,s] triplet:
			vec_dbl *w[14], **twid_ptrs = &(w[0]);	// Array of 14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DIT_TWIDDLE_OOP
			w[0x0] = tmp+ 1; w[0x1] = tmp+ 2;
			w[0x2] = tmp+ 4; w[0x3] = tmp+ 5;
			w[0x4] = tmp+ 7; w[0x5] = tmp+ 8;
			w[0x6] = tmp+10; w[0x7] = tmp+11;
			w[0x8] = tmp+13; w[0x9] = tmp+14;
			w[0xa] = tmp+16; w[0xb] = tmp+17;
			w[0xc] = tmp+19; w[0xd] = tmp+20;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				tm1,OFF1,
				tm2,dit_o_ptr_stride,
				twid_ptrs, two
			);
		}
	#endif	// 1/0?

#else	// USE_SSE2 = False, or non-64-bit-GCC:

/* Gather the needed data (128 64-bit complex, i.e. 256 64-bit reals) and do 8 twiddleless length-16 subtransforms: */

	tptr = t;
	for(l = 0; l < 8; l++) {
		k2 = poff[l<<2];	// poff[4*l] = p00,10,...,60,70
		jt = j1 + k2; jp = j2 + k2;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im,(tptr+9)->re,(tptr+9)->im,(tptr+10)->re,(tptr+10)->im,(tptr+11)->re,(tptr+11)->im,(tptr+12)->re,(tptr+12)->im,(tptr+13)->re,(tptr+13)->im,(tptr+14)->re,(tptr+14)->im,(tptr+15)->re,(tptr+15)->im,
			c16,s16
		);
		tptr += 16;
	}

/*...and now do 16 radix-8 subtransforms, including the internal twiddle factors - we use the same positive-power
roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f],
as are the index offsets of each sets of complex outputs in the A-array: [jt,jp] + p10*[0,4,2,6,1,5,3,7]:
*/
	/* Block 0: */
	tptr = t;	jt = j1;	jp = j2;
	/* 0-index block has all-unity twiddles: Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in-order here: */
	RADIX_08_DIT_OOP(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,
		a[jt],a[jp],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70]
	);

	// Remaining 15 sets of macro calls done in loop:
	for(l = 1; l < 16; l++) {
		tptr = t + reverse(l,4);
		k2 = po_br[l];	// po_br[] = p[084c2a6e195d3b7f]
		jt = j1 + k2; jp = j2 + k2;
		const double *addr = DFT128_TWIDDLES[l], *addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_08_DIT_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c)
		);
	}

#endif	// USE_SSE2 ?

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
		// incr must divide nloop [RADIX/8 = 16 or RADIX/16 = 8, depending on whether we use 8-or-16-way carry macros]!
	  #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
	  #else
		const uint32 nloop = RADIX>>3;
	  #endif
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy_r; itmp = bjmodn;
	  #ifndef USE_AVX512
		tm2 = cy_r+1; itm2 = bjmodn+4; 	// tm2,itm2 not used in AVX-512 mode
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
				AVX_cmplx_carry_fast_pow2_errcheck_X16(tmp, tm1,     itmp, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p01,p02,p03,p04, addr);
				tmp += 32; tm1 += 2;           itmp += 16; i = 0;
			  #else
				AVX_cmplx_carry_fast_pow2_errcheck_X8 (tmp, tm1,     itmp, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p01,p02,p03,p04, addr);
				tmp += 16; tm1 += 1;           itmp += 8; i = 0;
			  #endif
			#else
				AVX_cmplx_carry_fast_pow2_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p01,p02,p03,p04, addr);
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
			AVX_cmplx_carry_norm_pow2_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p01,p02,p03, addr);
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
			SSE2_cmplx_carry_fast_pow2_wtsinit(add1,add2,add3, add0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, k0,k1,k2,k3, sse_bw,sse_nm1)

			for(l = loop; l < loop+incr; l++) {
				// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data:
				add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
				SSE2_cmplx_carry_fast_pow2_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p01,p02,p03, addr);
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
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_pow2_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p01, addr);
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
			add0 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_pow2_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p02,p03, addr);
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
			cmplx_carry_norm_pow2_errcheck0(a[jt    ],a[jp    ],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p01],a[jp+p01],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p02],a[jp+p02],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_pow2_errcheck (a[jt+p03],a[jp+p03],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
		}

	  } else {	// HiACC:

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; itmp = bjmodn;
		double *addr_ = cy_r;
		for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...
			cmplx_carry_norm_pow2_errcheck0(a[jt    ],a[jp    ],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p01],a[jp+p01],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p02],a[jp+p02],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_pow2_errcheck (a[jt+p03],a[jp+p03],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
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
		tm0 = s1p00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x4000;
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
		tm0 = s1p00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x2000;
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

		tm1 = s1p00; tmp = cy_r;	/* <*** Again rely on contiguity of cy_r,i here ***/
		for(l = 0; l < RADIX>>1; l++) {	// RADIX/2 loop passes
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[l>>1]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
			tm2 += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_fermat_carry_norm_pow2_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, tm2,p01, add0);
			tm1 += 4; tmp += 2;
		}

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2);
		ntmp = 0;
		double *addr_ = cy_r, *addi_ = cy_i;
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];	// poff[] = p04,08,...
			fermat_carry_norm_pow2_errcheck(a[jt    ],a[jp    ],*addr_,*addi_,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr_; ++addi_;
			fermat_carry_norm_pow2_errcheck(a[jt+p01],a[jp+p01],*addr_,*addi_,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr_; ++addi_;
			fermat_carry_norm_pow2_errcheck(a[jt+p02],a[jp+p02],*addr_,*addi_,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr_; ++addi_;
			fermat_carry_norm_pow2_errcheck(a[jt+p03],a[jp+p03],*addr_,*addi_,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr_; ++addi_;
		}

	#endif	/* #ifdef USE_SSE2 */

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-128 DIF pass is here:	*/

#ifdef USE_SSE2

	#if 1

		int dif_offsets[RADIX];
		dif_offsets[0x0] =   0; dif_offsets[0x1] = p01; dif_offsets[0x2] = p02; dif_offsets[0x3] = p03; dif_offsets[0x4] = p04; dif_offsets[0x5] = p05; dif_offsets[0x6] = p06; dif_offsets[0x7] = p07;
		dif_offsets[0x8] = p08; dif_offsets[0x9] = p09; dif_offsets[0xa] = p0a; dif_offsets[0xb] = p0b; dif_offsets[0xc] = p0c; dif_offsets[0xd] = p0d; dif_offsets[0xe] = p0e; dif_offsets[0xf] = p0f;
		for(l = 0; l < 16; l++) {
			dif_offsets[l+0x10] = dif_offsets[l] + p10;
			dif_offsets[l+0x20] = dif_offsets[l] + p20;
			dif_offsets[l+0x30] = dif_offsets[l] + p30;
			dif_offsets[l+0x40] = dif_offsets[l] + p40;
			dif_offsets[l+0x50] = dif_offsets[l] + p50;
			dif_offsets[l+0x60] = dif_offsets[l] + p60;
			dif_offsets[l+0x70] = dif_offsets[l] + p70;
		}
		SSE2_RADIX128_DIF(	// Inputs in s1p00[ 00- 7f]
			s1p00,
			r00, isrt2,two,twid0,
			(a+j1     ),dif_offsets
		);

	#else

	// Gather the needed data and do 8 twiddleless length-16 subtransforms, with p-offsets in br8 order: 04261537:
	// NOTE that unlike the RADIX_08_DIF_OOP() macro used for pass 1 of the radix-64 DFT, RADIX_16_DIF outputs are IN-ORDER rather than BR:

	  #ifndef USE_ARM_V8_SIMD	// These were last def'd prior to Pass 2 of above carry-preceding 128-DIT
		#undef OFF1
		#undef OFF2
		#undef OFF4
	  #endif
	  #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x100;
		OFF2 = 0x200;
		OFF3 = 0x300;
		OFF4 = 0x400;
	  #elif defined(USE_AVX512)
		#define OFF1	4*0x100
		#define OFF2	4*0x200
		#define OFF3	4*0x300
		#define OFF4	4*0x400
	  #elif defined(USE_AVX)
		#define OFF1	0x200
		#define OFF2	0x400
		#define OFF3	0x600
		#define OFF4	0x800
	  #else
		#define OFF1	0x100
		#define OFF2	0x200
		#define OFF3	0x300
		#define OFF4	0x400
	  #endif

	// NOTE that unlike the RADIX_08_DIF_OOP() macro used for pass 1 of the radix-64 DFT, RADIX_16_DIF outputs are IN-ORDER rather than BR:
		// Start of isrt2/cc16/ss16 triplet needed for radix-16 SSE2 DFT macros [use instead of isrt2 here]
	  #if COMPACT_OBJ
		tmp = isrt2;	// In COMPACT_OBJ mode, use that [isrt2] ptr heads up the needed isrt2/cc16/ss16 pointer-triplet
						// (actually a quartet in AVX2 mode, since there need sqrt2 preceding isrt2)
	  #else
		tmp = sqrt2+1;	// In non-compact-obj-code mode, access the same (unnamed) pointer this way
	  #endif
		tm1 = r00;
		for(l = 0; l < 8; l++) {
			k1 = reverse(l,3)<<1;
			tm2 = s1p00 + k1;
			SSE2_RADIX16_DIF_0TWIDDLE_B(tm2,OFF1,OFF2,OFF3,OFF4, tmp,two, tm1);
			tm1 += 32;
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

	/*...and now do 16 radix-8 subtransforms, including the internal twiddle factors: */
	/* Note: Each of the 7 sincos-pairs sent to SSE2_RADIX8_DIF_TWIDDLE_OOP is actually a vector-data triplet
	whose first (unnamed) member contains the isrt2 needed for radix-8.
	This is a workaround of GCC's 30-arg limit for inline ASM macros, which proves a royal pain here.
	*/
	   #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x200;
		OFF2 = 0x400;
		OFF3 = 0x600;
		OFF4 = 0x800;
		OFF5 = 0xa00;
		OFF6 = 0xc00;
		OFF7 = 0xe00;
	  #elif defined(USE_AVX512)
		#define OFF1	0x0800
		#define OFF2	0x1000
		#define OFF3	0x1800
		#define OFF4	0x2000
		#define OFF5	0x2800
		#define OFF6	0x3000
		#define OFF7	0x3800
	  #elif defined(USE_AVX)
		#define OFF1	0x0400
		#define OFF2	0x0800
		#define OFF3	0x0c00
		#define OFF4	0x1000
		#define OFF5	0x1400
		#define OFF6	0x1800
		#define OFF7	0x1c00
	  #else
		#define OFF1	0x200
		#define OFF2	0x400
		#define OFF3	0x600
		#define OFF4	0x800
		#define OFF5	0xa00
		#define OFF6	0xc00
		#define OFF7	0xe00
	   #endif
		/* Block 0: */
		add0 = &a[j1]      ; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
		/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
		so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
		SSE2_RADIX8_DIF_0TWIDDLE(
			r00, OFF4,OFF2,OFF6,OFF1,OFF5,OFF3,OFF7,	//	r00,r40,r20,r60,r10,r50,r30,r70
		  #ifdef USE_AVX2
			add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two
		  #else
			add0,add1,add2,add3,add4,add5,add6,add7, isrt2
		  #endif
		);

		uint32 *ui32_ptr = &(po_lin[0]);	// Kludge for Arm-SIMD builds - with just po_lin as macro arg, gcc was feeding the low 8 bytes of po_lin
		// po_lin[1,0]as arg to the macro, rather than po_lin-as-pointer. Can't feed the latter directly as arg to macro, because that causes gcc to
		// emit "error: memory input 5 is not directly addressable" error.

		// Remaining 15 sets of macro calls done in loop:
		for(l = 1; l < 16; l++) {
			k1 = reverse(l,4);
			// I-ptrs and sincos-ptrs; Twid-offsets are multiples of 21 vec_dbl; +1 to point to cos [middle] term of each [isrt2,c,s] triplet
			tm1 = r00 + (l<<1); tmp = twid0 + (k1<<4)+(k1<<2)+k1;
			vec_dbl *w[14], **twid_ptrs = &(w[0]);	// Array of 14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DIT_TWIDDLE_OOP
			w[0x0] = tmp+ 1; w[0x1] = tmp+ 2;
			w[0x2] = tmp+ 4; w[0x3] = tmp+ 5;
			w[0x4] = tmp+ 7; w[0x5] = tmp+ 8;
			w[0x6] = tmp+10; w[0x7] = tmp+11;
			w[0x8] = tmp+13; w[0x9] = tmp+14;
			w[0xa] = tmp+16; w[0xb] = tmp+17;
			w[0xc] = tmp+19; w[0xd] = tmp+20;
			// O-ptrs: a[j1] offset here = p08,p10,p18,...
			add0 = &a[j1] + poff[l<<1];	// Dec 2020: Apple M1/Clang patch cuts SSE2_RADIX8_DIF_TWIDDLE_OOP #args from 30 to 21
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				tm1,OFF1,
				add0,ui32_ptr,
				twid_ptrs, two
			);
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
		#undef OFF5
		#undef OFF6
		#undef OFF7
	  #endif

	#endif	// 1/0?

#else	// USE_SSE2 = False:

// Gather the needed data and do 8 twiddleless length-16 subtransforms, with p-offsets in br8 order: 04261537:
// NOTE that unlike the RADIX_08_DIF_OOP() macro used for pass 1 of the radix-64 DFT, RADIX_16_DIF outputs are IN-ORDER rather than BR:

	tptr = t;
	for(l = 0; l < 8; l++) {
		k2 = po_br[l+l];	// po_br[2*l] = p[04261537]
		jt = j1 + k2; jp = j2 + k2;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im,(tptr+9)->re,(tptr+9)->im,(tptr+10)->re,(tptr+10)->im,(tptr+11)->re,(tptr+11)->im,(tptr+12)->re,(tptr+12)->im,(tptr+13)->re,(tptr+13)->im,(tptr+14)->re,(tptr+14)->im,(tptr+15)->re,(tptr+15)->im,
			c16,s16
		);
		tptr += 16;
	}

/*...and now do 16 radix-8 subtransforms, including the internal twiddle factors: */
	// Block 0: t*0ri, has all-unity twiddles
	tptr = t;	jt = j1;	jp = j2;
	// Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in BR-order here [swap index pairs 1/4 and 3/6]:
	RADIX_08_DIF_OOP(
		tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x70)->re,(tptr+0x70)->im,
		a[jt],a[jp],a[jt+p04],a[jp+p04],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07]
	);	tptr++;

	// Remaining 15 sets of macro calls done in loop:
	for(l = 1; l < 16; l++) {
		k2 = poff[l+l];	// poffs[2*l] = p00,08,10,...,70,78
		jt = j1 + k2; jp = j2 + k2;
		const double *addr = DFT128_TWIDDLES[l], *addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_08_DIF_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c)
		);	tptr++;
	}

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
