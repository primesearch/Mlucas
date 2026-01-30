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
	for(j = jstart; j < jhi; j += stride)
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	/*...The radix-320 DIT pass is here:	*/
	#ifdef USE_SSE2

	/*...gather the needed data (160 64-bit complex, i.e. 320 64-bit reals) and do 5 radix-32 transforms...*/
	/*	tmp = r00;
		for(l = 0; l < ODD_RADIX; l++) {
			add0 = &a[j1+dit_phi[l]];	itmp = (int *)dit_offsets+(l<<5);
			SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	tmp += 64;
		}*/
	//...gather the needed data (320 64-bit complex) and do 5 radix-64 transforms:
		// Use s1p00-... for scratch for these 5 DFTs:
	  #if 0	//*** Parametrized-loop version slower, likely bc (unlike 32-DFT and smaller) the 64-DFT is a function call rather an inlined macro
		tmp = r00;
		for(l = 0; l < ODD_RADIX; l++) {
			jt = j1+poff[l<<4]; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,dit_i_offsets+(l<<6), s1p00, tmp,t_offsets);	tmp += 0x80;
		}
	  #else
		jt = j1    ; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,(dit_i_offsets+0x000), s1p00, r00      ,t_offsets);
		jt = j1+p40; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,(dit_i_offsets+0x040), s1p00, r00+0x080,t_offsets);
		jt = j1+p80; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,(dit_i_offsets+0x080), s1p00, r00+0x100,t_offsets);
		jt = j1+pc0; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,(dit_i_offsets+0x0c0), s1p00, r00+0x180,t_offsets);
		jt = j1+pg0; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,(dit_i_offsets+0x100), s1p00, r00+0x200,t_offsets);
	  #endif
	//...and now do 64 radix-5 transforms:

	   #if (OS_BITS == 64) && defined(USE_AVX2) && !defined(USE_AVX512)	//*** Nov 2018: Timing tests on Core2[sse2] and KNL[avx512] both show using undoubled radix-5 macro is faster ***

		tmp = r00;
		for(l = 0; l < 64; l += 2) {
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;	va1 = tmp + 0x80;	va2 = tmp + 0x100;	va3 = tmp + 0x180;	va4 = tmp + 0x200;
			// Output pointers are into s1p** memblock:
			int lm64 = 64-l, lne0 = (l != 0);
			int k = (lm64>>4) & (-lne0); k += (k<<3);// Which-5-perm-index, *= 9
			int ic = ((l + lne0) % 5) + k;			// Leftward cshift count plus 9-elts offset per increment in k
			k = plo[lm64 & 15];						// Re-use k for p-offset index, plo[(64-i) % 16]
			k0 = dit_p20_cperms[ic]; k1 = dit_p20_cperms[ic+1]; k2 = dit_p20_cperms[ic+2]; k3 = dit_p20_cperms[ic+3]; k4 = dit_p20_cperms[ic+4];
			tm2 = s1p00 + k;	vb0 = tm2 + k0;	vb1 = tm2 + k1;	vb2 = tm2 + k2;	vb3 = tm2 + k3;	vb4 = tm2 + k4;
		// Now process 2nd block of indices before calling doubled-data version of radix-5 DFT macro:
			tmp += 2;
			wa0 = tmp;	wa1 = tmp + 0x80;	wa2 = tmp + 0x100;	wa3 = tmp + 0x180;	wa4 = tmp + 0x200;
			lm64--; lne0 = 1;	// Replace l --> l+1, i.e. l cannot == 0 in 2nd block
			k = (lm64>>4) & (-lne0); k += (k<<3);	// Which-5-perm-index, *= 9
			ic = (((l+1) + lne0) % 5) + k;			// Leftward cshift count plus 9-elts offset per increment in k
			k = plo[lm64 & 15];						// Re-use k for p-offset index, plo[(64-i) % 16]
			k0 = dit_p20_cperms[ic]; k1 = dit_p20_cperms[ic+1]; k2 = dit_p20_cperms[ic+2]; k3 = dit_p20_cperms[ic+3]; k4 = dit_p20_cperms[ic+4];
			tm2 = s1p00 + k;	wb0 = tm2 + k0;	wb1 = tm2 + k1;	wb2 = tm2 + k2;	wb3 = tm2 + k3;	wb4 = tm2 + k4;

			SSE2_RADIX_05_DFT_0TWIDDLE_X2( ycc1,two,
				va0,va1,va2,va3,va4, vb0,vb1,vb2,vb3,vb4,
				wa0,wa1,wa2,wa3,wa4, wb0,wb1,wb2,wb3,wb4
			); tmp += 2;
		}

	   #else

		tmp = r00;
		for(l = 0; l < 64; l++) {
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;	va1 = tmp + 0x080;	va2 = tmp + 0x100;	va3 = tmp + 0x180;	va4 = tmp + 0x200;
			// Output pointers are into s1p** memblock:
			int lm64 = 64-l, lne0 = (l != 0);
			int k = (lm64>>4) & (-lne0); k += (k<<3);// Which-5-perm-index, *= 9
		#if 1	// On-the-fly mod faster in my Core2 tests:
			int ic = ((l + lne0) % 5) + k;			// Leftward cshift count plus 9-elts offset per increment in k
		#else	// (l + lne0) ranges from 0-64 (skipping 2), try a predef of the resulting (mod 5) values:
			const uint8 mod5[64] = {0,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4};
			int ic = mod5[l] + k;
		#endif
			k = plo[lm64 & 15];						// Re-use k for p-offset index, plo[(64-i) % 16]
			k0 = dit_p20_cperms[ic]; k1 = dit_p20_cperms[ic+1]; k2 = dit_p20_cperms[ic+2]; k3 = dit_p20_cperms[ic+3]; k4 = dit_p20_cperms[ic+4];
			vb0 = s1p00 + (k0+k);
			vb1 = s1p00 + (k1+k);
			vb2 = s1p00 + (k2+k);
			vb3 = s1p00 + (k3+k);
			vb4 = s1p00 + (k4+k);
			SSE2_RADIX_05_DFT_0TWIDDLE(
				va0,va1,va2,va3,va4,
				ycc1,
				vb0,vb1,vb2,vb3,vb4
			);	tmp += 2;
		}

	  #endif

	#else	// USE_SSE2 = False:

	//...gather the needed data (320 64-bit complex) and do 5 radix-64 transforms:
		jt = j1    ;	RADIX_64_DIT(a+jt,(dit_i_offsets+0x000),RE_IM_STRIDE, (double *)(t+0x000),t_offsets,1);	// Outputs in t[00-63]
		jt = j1+p40;	RADIX_64_DIT(a+jt,(dit_i_offsets+0x040),RE_IM_STRIDE, (double *)(t+0x040),t_offsets,1);	// Outputs in t[64-127]
		jt = j1+p80;	RADIX_64_DIT(a+jt,(dit_i_offsets+0x080),RE_IM_STRIDE, (double *)(t+0x080),t_offsets,1);	// Outputs in t[64-127]
		jt = j1+pc0;	RADIX_64_DIT(a+jt,(dit_i_offsets+0x0c0),RE_IM_STRIDE, (double *)(t+0x0c0),t_offsets,1);	// Outputs in t[64-127]
		jt = j1+pg0;	RADIX_64_DIT(a+jt,(dit_i_offsets+0x100),RE_IM_STRIDE, (double *)(t+0x100),t_offsets,1);	// Outputs in t[64-127]
	//...and now do 64 radix-5 transforms, with the columns of t*[r,i] output pairs in the above 5x radix-64 set now acting as input rows:
		tptr = t;
		for(l = 0; l < 64; l++) {
		#if 0	// Initial try, used to compute needed output-permutation:
			int k = plo[l&15] + (phi[l>>4]);	//*** Must use *local*-defined k to not collide with the above outer-loop index ***
				k0 =   0;	k1 = p40;	k2 = p80;	k3 = pc0;	k4 = pg0;
		#else
			int lm64 = 64-l, lne0 = (l != 0);
			int k = (lm64>>4) & (-lne0); k += (k<<3);// Which-5-perm-index, *= 9
			int ic = ((l + lne0) % 5) + k;			// Leftward cshift count plus 9-elts offset per increment in k
			k = plo[lm64 & 15];						// Re-use k for p-offset index, plo[(64-i) % 16]
			k0 = dit_p20_cperms[ic]; k1 = dit_p20_cperms[ic+1]; k2 = dit_p20_cperms[ic+2]; k3 = dit_p20_cperms[ic+3]; k4 = dit_p20_cperms[ic+4];
		#endif
			jt = j1+k; jp = j2+k;
			RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,
				tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x100)->re,(tptr+0x100)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],
				rt,it);
			tptr++;
		}

	#endif	/* USE_SSE2 */

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

	  if(incr) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of incr
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		uint32 ii,loop, co2save = co2;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every 8th pass or better:
		// incr must divide nloop [RADIX/8 = 40 or RADIX/16 = 20, depending on whether we use 8-or-16-way carry macros]!
	  #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
	  #else
		const uint32 nloop = RADIX>>3;
	  #endif
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy; itmp = bjmodn;
	  #ifndef USE_AVX512
		tm2 = cy+1; itm2 = bjmodn+4;	// tm2,itm2 not used in AVX-512 mode
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
		tm1 = s1p00; tmp = cy; itmp = bjmodn;
		addr = &prp_mult;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];
			AVX_cmplx_carry_norm_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr);
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
		uint32 ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
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

/*...The radix-320 DIF pass is here:	*/

	#ifdef USE_SSE2
//if(!j && full_pass)printf("s1p[0-3] = %20.15f,%20.15f,%20.15f,%20.15f\n",s1p00->d0,(s1p00+1)->d0,s1p00->d1,(s1p00+1)->d1);
	/*...gather the needed data (320 64-bit complex, i.e. 640 64-bit reals) and do 64 radix-5 transforms...*/

	   #if (OS_BITS == 64) && defined(USE_AVX2) && !defined(USE_AVX512)

		tmp = r00;
		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head  of each trio
		k0 = 0;	k1 = 0x100<<1;	k2 = 0xc0<<1;	k3 = 0x80<<1;	k4 = 0x40<<1;
		for(l = 0; l < 64; l += 2) {
			// Input pointers are into s1p** memblock:
			int k = ilo<<1;
			tm2 = s1p00 + k;	vb0 = tm2 + k0;	vb1 = tm2 + k1;	vb2 = tm2 + k2;	vb3 = tm2 + k3;	vb4 = tm2 + k4;
			// Output-ptrs [va/vb swap roles here vs DIT] are regular-stride offsets of r00:
			va0 = tmp;	va1 = tmp + 0x80;	va2 = tmp + 0x100;	va3 = tmp + 0x180;	va4 = tmp + 0x200;
			ilo -= 5;
			if(ilo < 0) {
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 20);
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 20);
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 20);
				k3 = k2 - 4;	hi_neg = k3 < 0;	k3 += ((-hi_neg) & 20);
				k4 = k3 - 4;	hi_neg = k4 < 0;	k4 += ((-hi_neg) & 20);
				k0 = ihi<<5;	k1 = k1<<5;	k2 = k2<<5;	k3 = k3<<5;	k4 = k4<<5;
			}
		// Now process 2nd block of indices before calling doubled-data version of radix-5 DFT macro:
			tmp += 2;
			k = ilo<<1;
			tm2 = s1p00 + k;	wb0 = tm2 + k0;	wb1 = tm2 + k1;	wb2 = tm2 + k2;	wb3 = tm2 + k3;	wb4 = tm2 + k4;
			wa0 = tmp;	wa1 = tmp + 0x80;	wa2 = tmp + 0x100;	wa3 = tmp + 0x180;	wa4 = tmp + 0x200;
			ilo -= 5;
			if(ilo < 0) {
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 20);
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 20);
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 20);
				k3 = k2 - 4;	hi_neg = k3 < 0;	k3 += ((-hi_neg) & 20);
				k4 = k3 - 4;	hi_neg = k4 < 0;	k4 += ((-hi_neg) & 20);
				k0 = ihi<<5;	k1 = k1<<5;	k2 = k2<<5;	k3 = k3<<5;	k4 = k4<<5;
			}
			SSE2_RADIX_05_DFT_0TWIDDLE_X2( ycc1,two,
				vb0,vb1,vb2,vb3,vb4, va0,va1,va2,va3,va4,
				wb0,wb1,wb2,wb3,wb4, wa0,wa1,wa2,wa3,wa4
			); tmp += 2;
		}

	   #else

		tmp = r00;
		// SIMD: Scalar-code's plo[ilo] --> ilo<<1; phi[ihi] --> ihi<<5
		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head  of each trio
		k0 = 0;	k1 = 0x100<<1;	k2 = 0xc0<<1;	k3 = 0x80<<1;	k4 = 0x40<<1;
		for(l = 0; l < 64; l++) {
			// Input pointers are into s1p** memblock:
			int k = ilo<<1;
			tm2 = s1p00 + k;	vb0 = tm2 + k0;	vb1 = tm2 + k1;	vb2 = tm2 + k2;	vb3 = tm2 + k3;	vb4 = tm2 + k4;
			// Output-ptrs [va/vb swap roles here vs DIT] are regular-stride offsets of r00:
			va0 = tmp;	va1 = tmp + 0x80;	va2 = tmp + 0x100;	va3 = tmp + 0x180;	va4 = tmp + 0x200;
			SSE2_RADIX_05_DFT_0TWIDDLE(
				vb0,vb1,vb2,vb3,vb4,
				ycc1,
				va0,va1,va2,va3,va4
			);	tmp += 2;
			ilo -= 5;
			// If low nibble underflows need to add 16 to it, decrement high nibble ihi, and recompute k1-4.
			// Use branch here because dependent-op count > 10:
			if(ilo < 0) {
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 20);	// 1st high-part offset
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 2nd offset
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 3rd offset
				k3 = k2 - 4;	hi_neg = k3 < 0;	k3 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 4th offset
				k4 = k3 - 4;	hi_neg = k4 < 0;	k4 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 5th offset
				// Use distinct ihi and k0 = phi[ihi] for 1st high-part offset because basic-index ihi must survive from one loop pass to next:
				k0 = ihi<<5;	k1 = k1<<5;	k2 = k2<<5;	k3 = k3<<5;	k4 = k4<<5;
			}
		}

	  #endif

		// Use s1p00-3f for scratch for these 5 DFTs ... since transform length N = odd*64, pow2-shift arg = trailz(N) - trailz(64) = 0:
	  #if 0	//*** See note for 64-DIT at top of file; similar looped-version-slower seen here
		tmp = r00;
		for(l = 0; l < ODD_RADIX; l++) {
			jt = j1+poff[(80-(l<<4)) & -(l!=0)]; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)tmp,t_offsets, s1p00, a+jt,dif_o_offsets+(l<<6));	tmp += 0x80;
		}
	  #else
		jt = j1    ; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)(r00      ),t_offsets, s1p00, a+jt,dif_o_offsets+0x000);
		jt = j1+pg0; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)(r00+0x080),t_offsets, s1p00, a+jt,dif_o_offsets+0x040);
		jt = j1+pc0; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)(r00+0x100),t_offsets, s1p00, a+jt,dif_o_offsets+0x080);
		jt = j1+p80; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)(r00+0x180),t_offsets, s1p00, a+jt,dif_o_offsets+0x0c0);
		jt = j1+p40; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)(r00+0x200),t_offsets, s1p00, a+jt,dif_o_offsets+0x100);
	  #endif

	#else	// USE_SSE2 = False:

	//...gather the needed data (320 64-bit complex) and do 64 radix-5 transforms - We want unit-strides in the radix64-DFT macros, so use large output strides here:
		tptr = t;
		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head  of each trio
		k0 =   0;	k1 = pg0;	k2 = pc0;	k3 = p80;	k4 = p40;
		for(l = 0; l < 64; l++) {
			jp = plo[ilo];	// Same for each elt of the current-offset trio
			jt = j1 + jp; jp += j2;
			RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],
				tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x100)->re,(tptr+0x100)->im,
				rt,it);
			tptr++;
			ilo -= 5;
			// If low nibble underflows need to add 16 to it, decrement high nibble ihi, and recompute k1-4.
			// Use branch here because dependent-op count > 10:
			if(ilo < 0) {
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 20);	// 1st high-part offset
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 2nd offset
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 3rd offset
				k3 = k2 - 4;	hi_neg = k3 < 0;	k3 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 4th offset
				k4 = k3 - 4;	hi_neg = k4 < 0;	k4 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 5th offset
				// Use distinct ihi and k0 = phi[ihi] for 1st high-part offset because basic-index ihi must survive from one loop pass to next:
				k0 = phi[ihi];	k1 = phi[k1];	k2 = phi[k2];	k3 = phi[k3];	k4 = phi[k4];
			}
		}
	//...and now do 5 radix-64 transforms:
		//	NOTE that RADIX_64_DIF outputs are IN-ORDER rather than BR:
		jt = j1    ;	RADIX_64_DIF((double *)(t+0x000),t_offsets,1, (a+jt),(dif_o_offsets+0x000),RE_IM_STRIDE);	// Inputs in t[00-63]
		jt = j1+pg0;	RADIX_64_DIF((double *)(t+0x040),t_offsets,1, (a+jt),(dif_o_offsets+0x040),RE_IM_STRIDE);	// Inputs in t[64-127]
		jt = j1+pc0;	RADIX_64_DIF((double *)(t+0x080),t_offsets,1, (a+jt),(dif_o_offsets+0x080),RE_IM_STRIDE);	// Inputs in t[128-191]
		jt = j1+p80;	RADIX_64_DIF((double *)(t+0x0c0),t_offsets,1, (a+jt),(dif_o_offsets+0x0c0),RE_IM_STRIDE);	// Inputs in t[192-255]
		jt = j1+p40;	RADIX_64_DIF((double *)(t+0x100),t_offsets,1, (a+jt),(dif_o_offsets+0x100),RE_IM_STRIDE);	// Inputs in t[256-319]

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
