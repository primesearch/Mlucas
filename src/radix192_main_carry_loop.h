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

	/*...The radix-192 DIT pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (192 64-bit complex) and do 3 radix-64 transforms:
		// Use s1p00-3f for scratch for these 3 DFTs:
		jt = j1    ; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,(dit_i_offsets+0x00), s1p00, r00,t_offsets);
		jt = j1+p80; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,(dit_i_offsets+0x40), s1p00, r40,t_offsets);
		jt = j1+p40; SSE2_RADIX_64_DIT( FALSE, thr_id, a+jt,(dit_i_offsets+0x80), s1p00, r80,t_offsets);

	//...and now do 64 radix-3 transforms:
	  #if OS_BITS == 64

		tmp =   r00;	tm1 =   r00+0x80;	tm2 =   r00+0x100;	// Need to double the offsets from the p-stripped ones
		va0 = s1p00;	va1 = s1p00+0x80;	va2 = s1p00+0x100;	// since vec_dbl pointers but complex data
		ilo = 0;	// int tracking the low nibble of the p* offset at head of each trio
		idx = 2; jdx = 0; iptr = &p10_3perm[jdx][0];	// jdx is [A-D] index, idx is 0-2, so this formally corr. to A2...
		// ... but set the k0-2 offset values corr. to A0 [rather than A2] on loop entry to take care of leading A0:
		k0 =   0;	k1 = 0x80;	k2 = 0x100;
		for(l = 0; l < 64; l += 2) {
		// *** 1st set of DFT-3 data: ***
			jp = ilo<<1;
			vb0 = va0+jp;	vb1 = va1+jp;	vb2 = va2+jp;
			idx = (idx + 1 + (idx == 2)) & 3;
			ilo--;
			if(ilo < 0) {	// This branch only taken when L == 0 (mod 16), i.e. every 16th loop exec
				ilo += 16;
				jdx = (jdx + 3)&3;
				iptr = &p10_3perm[jdx][0];
			}
		// *** 2nd set of 2 DFT-3 data: ***
			k0 = *(iptr+idx);	k1 = *(iptr+idx+1);	k2 = *(iptr+idx+2);
			vc0 = tmp+2; vc1 = tm1+2; vc2 = tm2+2;						// 2nd set of iptrs
			jp = ilo<<1;
			vd0 = s1p00+k0+jp;	vd1 = s1p00+k1+jp;	vd2 = s1p00+k2+jp;	// 2nd set of optrs
			idx = (idx + 1 + (idx == 2)) & 3;
			ilo--;
			k0 = *(iptr+idx);	k1 = *(iptr+idx+1);	k2 = *(iptr+idx+2);
			SSE2_RADIX_03_DFT_X2(cc0, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
		// Set up for next loop pass:
			tmp+=4; tm1+=4; tm2+=4;
			va0 = s1p00+k0;	va1 = s1p00+k1;	va2 = s1p00+k2;
		}

	  #else

		tmp =   r00;	tm1 =   r00+0x80;	tm2 =   r00+0x100;	// Need to double the offsets from the p-stripped ones
		va0 = s1p00;	va1 = s1p00+0x80;	va2 = s1p00+0x100;	// since vec_dbl pointers but complex data
		ilo = 0;	// int tracking the low nibble of the p* offset at head of each trio
		idx = 2; jdx = 0; iptr = &p10_3perm[jdx][0];	// jdx is [A-D] index, idx is 0-2, so this formally corr. to A2...
		// ... but set the k0-2 offset values corr. to A0 [rather than A2] on loop entry to take care of leading A0:
		k0 =   0;	k1 = 0x80;	k2 = 0x100;
		/*
		Here is the index patterning for the loop:
			l =  0: idx,ilo =  2, 0		ilo < 0 branch! reset ilo = 15
			l =  1: idx,ilo =  0,15
			l =  2: idx,ilo =  1,14
			l =  3: idx,ilo =  2,13
			l =  4: idx,ilo =  0,12
			l =  5: idx,ilo =  1,11
			l =  6: idx,ilo =  2,10
			l =  7: idx,ilo =  0, 9
			l =  8: idx,ilo =  1, 8
			l =  9: idx,ilo =  2, 7
			l = 10: idx,ilo =  0, 6
			l = 11: idx,ilo =  1, 5
			l = 12: idx,ilo =  2, 4
			l = 13: idx,ilo =  0, 3
			l = 14: idx,ilo =  1, 2
			l = 15: idx,ilo =  2, 1
			l = 16: idx,ilo =  0, 0		ilo < 0 branch! reset ilo = 15
			l = 17: idx,ilo =  1,15
			l = 18: idx,ilo =  2,14
			l = 19: idx,ilo =  0,13
			l = 20: idx,ilo =  1,12
			l = 21: idx,ilo =  2,11
			l = 22: idx,ilo =  0,10
			l = 23: idx,ilo =  1, 9
			l = 24: idx,ilo =  2, 8
			l = 25: idx,ilo =  0, 7
			l = 26: idx,ilo =  1, 6
			l = 27: idx,ilo =  2, 5
			l = 28: idx,ilo =  0, 4
			l = 29: idx,ilo =  1, 3
			l = 30: idx,ilo =  2, 2
			l = 31: idx,ilo =  0, 1
			l = 32: idx,ilo =  1, 0		ilo < 0 branch! reset ilo = 15
			l = 33: idx,ilo =  2,15
			l = 34: idx,ilo =  0,14
			l = 35: idx,ilo =  1,13
			l = 36: idx,ilo =  2,12
			l = 37: idx,ilo =  0,11
			l = 38: idx,ilo =  1,10
			l = 39: idx,ilo =  2, 9
			l = 40: idx,ilo =  0, 8
			l = 41: idx,ilo =  1, 7
			l = 42: idx,ilo =  2, 6
			l = 43: idx,ilo =  0, 5
			l = 44: idx,ilo =  1, 4
			l = 45: idx,ilo =  2, 3
			l = 46: idx,ilo =  0, 2
			l = 47: idx,ilo =  1, 1
			l = 48: idx,ilo =  2, 0		ilo < 0 branch! reset ilo = 15
			l = 49: idx,ilo =  0,15
			l = 50: idx,ilo =  1,14
			l = 51: idx,ilo =  2,13
			l = 52: idx,ilo =  0,12
			l = 53: idx,ilo =  1,11
			l = 54: idx,ilo =  2,10
			l = 55: idx,ilo =  0, 9
			l = 56: idx,ilo =  1, 8
			l = 57: idx,ilo =  2, 7
			l = 58: idx,ilo =  0, 6
			l = 59: idx,ilo =  1, 5
			l = 60: idx,ilo =  2, 4
			l = 61: idx,ilo =  0, 3
			l = 62: idx,ilo =  1, 2
			l = 63: idx,ilo =  2, 1
		*/
		for(l = 0; l < 64; l++) {
			jp = ilo<<1;	// Same for each elt of the current-offset trio
			vb0 = va0+jp;	vb1 = va1+jp;	vb2 = va2+jp;
			SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc0, vb0,vb1,vb2);
			tmp+=2; tm1+=2; tm2+=2;
			// Incr [A-D]-perm cshift counter (mod 3):
			idx++;
			idx += (idx == 3);	// When idx hits 3, bump up to 4...
			idx &= 3;	// ...and reduce (mod 4), which is cheaper than (mod 3)
			// If low nibble underflows, += 16 and decrement [A-D] index (mod 4), which ==> jdx += 3 (mod 4):
			ilo--;
			if(ilo < 0) {	// This branch only taken every 16th loop exec
				ilo += 16;
				jdx = (jdx + 3)&3;
				iptr = &p10_3perm[jdx][0];
			}
			k0 = *(iptr+idx);	k1 = *(iptr+idx+1);	k2 = *(iptr+idx+2);
			va0 = s1p00+k0;	va1 = s1p00+k1;	va2 = s1p00+k2;
		}

	  #endif

	#else	/* !USE_SSE2 */

	//...gather the needed data (192 64-bit complex) and do 3 radix-64 transforms:
		jt = j1    ;	RADIX_64_DIT(a+jt,(dit_i_offsets+0x00),RE_IM_STRIDE, (double *)(t+0x00),t_offsets,1);	// Outputs in t[00-63]
		jt = j1+p80;	RADIX_64_DIT(a+jt,(dit_i_offsets+0x40),RE_IM_STRIDE, (double *)(t+0x40),t_offsets,1);	// Outputs in t[64-127]
		jt = j1+p40;	RADIX_64_DIT(a+jt,(dit_i_offsets+0x80),RE_IM_STRIDE, (double *)(t+0x80),t_offsets,1);	// Outputs in t[64-127]

	//...and now do 64 radix-3 transforms:
		tptr = t;
		ilo = 0;	// int tracking the low nibble of the p* offset at head of each trio
		idx = 2; jdx = 0; iptr = &p10_3perm[jdx][0];	// jdx is [A-D] index, idx is 0-2, so this formally corr. to A2...
		// ... but set the k0-2 offset values corr. to A0 [rather than A2] on loop entry to take care of leading A0:
		k0 =   0;	k1 = p40;	k2 = p80;
		for(l = 0; l < 64; l++) {
			jp = plo[ilo];	// Same for each elt of the current-offset trio
			jt = j1 + jp; jp += j2;
			RADIX_03_DFT(s,c3m1,
				tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x80)->re,(tptr+0x80)->im,
				t00,t01,t02,t03,t04,t05,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]
			);	tptr++;
			// Incr [A-D]-perm cshift counter (mod 3):
			idx++;
			idx += (idx == 3);	// When idx hits 3, bump up to 4...
			idx &= 3;	// ...and reduce (mod 4), which is cheaper than (mod 3)
			// If low nibble underflows, += 16 and decrement [A-D] index (mod 4), which ==> jdx += 3 (mod 4):
			ilo--;
			if(ilo < 0) {	// This branch only taken every 16th loop exec
				ilo += 16;
				jdx = (jdx + 3)&3;
				iptr = &p10_3perm[jdx][0];
			}
			k0 = *(iptr+idx);	k1 = *(iptr+idx+1);	k2 = *(iptr+idx+2);
		}

	#endif

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 48 separate blocks of the A-array, we need 48 separate carries.	*/

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

	  if(incr) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of the low
						// inc_arr element to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		uint32 ii,loop, co2save = co2;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every 8th pass or better:
		// incr must divide nloop [RADIX/8 = 24 or RADIX/16 = 12, depending on whether we use 8-or-16-way carry macros]!
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
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[l]);
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
		// DIF/DIT loops use k1,k2 and there those mustbbe signed, so use i0-3 here:
		uint32 i0,i1,i2,i3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;

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

	/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];
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

	#endif	// USE_AVX?

	/*...The radix-192 DIF pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (192 64-bit complex) and do 16 radix-3 transforms:
	  #if OS_BITS == 64

		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head of each trio
		va0 = s1p00;	va1 = s1p00+0x100;	va2 = s1p00+0x080;
		tmp =   r00;	tm1 =   r00+0x080;	tm2 =   r00+0x100;
		for(l = 0; l < 64; l += 2) {
		// *** 1st set of DFT-3 data: ***
			jp = ilo<<1;	// Same for each elt of the current-offset trio
			vb0 = va0+jp;	vb1 = va1+jp;	vb2 = va2+jp;
			ilo -= 3;
			// If low nibble underflows need to add 16 to it, decrement high nibble ihi, and recompute k1,k2.
			// Use branch here because dependent-op count > 10:
			if(ilo < 0) {	// This branch taken every 5th or 6th loop exec, but never twice in a given loop pass
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 12);	// 1st high-part offset
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 2nd offset
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 3rd offset
				k0 = ihi<<5;	k1 = k1<<5;	k2 = k2<<5;	// Use distinct ihi and k0 = ihi<<5 for 1st high-part offset because basic-index ihi must survive from one loop pass to next
				va0 = s1p00+k0;	va1 = s1p00+k1;	va2 = s1p00+k2;
			// *** 2nd set of DFT-3 data: ***
				jp = ilo<<1;	// Same for each elt of the current-offset trio
				vc0 = va0+jp;	vc1 = va1+jp;	vc2 = va2+jp;				// 2nd set of iptrs
				vd0 = tmp+2; vd1 = tm1+2; vd2 = tm2+2;						// 2nd set of optrs
				ilo -= 3;
				SSE2_RADIX_03_DFT_X2(cc0, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			}
			else
			{
			// *** 2nd set of DFT-3 data: ***
				jp = ilo<<1;	// Same for each elt of the current-offset trio
				vc0 = va0+jp;	vc1 = va1+jp;	vc2 = va2+jp;				// 2nd set of iptrs
				vd0 = tmp+2; vd1 = tm1+2; vd2 = tm2+2;						// 2nd set of optrs
				ilo -= 3;
				if(ilo < 0) {	// This branch taken every 5th or 6th loop exec, thus need in both DFT-3 sections
					ilo += 16;
					ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 12);	// 1st high-part offset
					k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 2nd offset
					k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 3rd offset
					k0 = ihi<<5;	k1 = k1<<5;	k2 = k2<<5;	// Use distinct ihi and k0 = ihi<<5 for 1st high-part offset because basic-index ihi must survive from one loop pass to next
					va0 = s1p00+k0;	va1 = s1p00+k1;	va2 = s1p00+k2;
				}
				SSE2_RADIX_03_DFT_X2(cc0, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			}
			tmp+=4; tm1+=4; tm2+=4;
		}

	  #else

		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head of each trio
		va0 = s1p00;	va1 = s1p00+0x100;	va2 = s1p00+0x080;
		tmp =   r00;	tm1 =   r00+0x080;	tm2 =   r00+0x100;
		/*
		Here is the index patterning for the loop:
			l =  0: ihi,ilo =  0, 0
			ilo < 0 branch! reset ilo = -3
			l =  1: ihi,ilo = 11,13
			l =  2: ihi,ilo = 11,10
			l =  3: ihi,ilo = 11, 7
			l =  4: ihi,ilo = 11, 4
			l =  5: ihi,ilo = 11, 1
			ilo < 0 branch! reset ilo = -2
			l =  6: ihi,ilo = 10,14
			l =  7: ihi,ilo = 10,11
			l =  8: ihi,ilo = 10, 8
			l =  9: ihi,ilo = 10, 5
			l = 10: ihi,ilo = 10, 2
			ilo < 0 branch! reset ilo = -1
			l = 11: ihi,ilo =  9,15
			l = 12: ihi,ilo =  9,12
			l = 13: ihi,ilo =  9, 9
			l = 14: ihi,ilo =  9, 6
			l = 15: ihi,ilo =  9, 3
			l = 16: ihi,ilo =  9, 0
			ilo < 0 branch! reset ilo = -3
			l = 17: ihi,ilo =  8,13
			l = 18: ihi,ilo =  8,10
			l = 19: ihi,ilo =  8, 7
			l = 20: ihi,ilo =  8, 4
			l = 21: ihi,ilo =  8, 1
			ilo < 0 branch! reset ilo = -2
			l = 22: ihi,ilo =  7,14
			l = 23: ihi,ilo =  7,11
			l = 24: ihi,ilo =  7, 8
			l = 25: ihi,ilo =  7, 5
			l = 26: ihi,ilo =  7, 2
			ilo < 0 branch! reset ilo = -1
			l = 27: ihi,ilo =  6,15
			l = 28: ihi,ilo =  6,12
			l = 29: ihi,ilo =  6, 9
			l = 30: ihi,ilo =  6, 6
			l = 31: ihi,ilo =  6, 3
			l = 32: ihi,ilo =  6, 0
			ilo < 0 branch! reset ilo = -3
			l = 33: ihi,ilo =  5,13
			l = 34: ihi,ilo =  5,10
			l = 35: ihi,ilo =  5, 7
			l = 36: ihi,ilo =  5, 4
			l = 37: ihi,ilo =  5, 1
			ilo < 0 branch! reset ilo = -2
			l = 38: ihi,ilo =  4,14
			l = 39: ihi,ilo =  4,11
			l = 40: ihi,ilo =  4, 8
			l = 41: ihi,ilo =  4, 5
			l = 42: ihi,ilo =  4, 2
			ilo < 0 branch! reset ilo = -1
			l = 43: ihi,ilo =  3,15
			l = 44: ihi,ilo =  3,12
			l = 45: ihi,ilo =  3, 9
			l = 46: ihi,ilo =  3, 6
			l = 47: ihi,ilo =  3, 3
			l = 48: ihi,ilo =  3, 0
			ilo < 0 branch! reset ilo = -3
			l = 49: ihi,ilo =  2,13
			l = 50: ihi,ilo =  2,10
			l = 51: ihi,ilo =  2, 7
			l = 52: ihi,ilo =  2, 4
			l = 53: ihi,ilo =  2, 1
			ilo < 0 branch! reset ilo = -2
			l = 54: ihi,ilo =  1,14
			l = 55: ihi,ilo =  1,11
			l = 56: ihi,ilo =  1, 8
			l = 57: ihi,ilo =  1, 5
			l = 58: ihi,ilo =  1, 2
			ilo < 0 branch! reset ilo = -1
			l = 59: ihi,ilo =  0,15
			l = 60: ihi,ilo =  0,12
			l = 61: ihi,ilo =  0, 9
			l = 62: ihi,ilo =  0, 6
			l = 63: ihi,ilo =  0, 3
		*/
		for(l = 0; l < 64; l++) {
			jp = ilo<<1;	// Same for each elt of the current-offset trio
			vb0 = va0+jp;	vb1 = va1+jp;	vb2 = va2+jp;
			SSE2_RADIX_03_DFT(vb0,vb1,vb2, cc0, tmp,tm1,tm2);
			tmp+=2; tm1+=2; tm2+=2;
			ilo -= 3;
			// If low nibble underflows need to add 16 to it, decrement high nibble ihi, and recompute k1,k2.
			// Use branch here because dependent-op count > 10:
			if(ilo < 0) {	// This branch only taken every 16th loop exec
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 12);	// 1st high-part offset
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 2nd offset
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 3rd offset
				k0 = ihi<<5;	k1 = k1<<5;	k2 = k2<<5;	// Use distinct ihi and k0 = ihi<<5 for 1st high-part offset because basic-index ihi must survive from one loop pass to next
				va0 = s1p00+k0;	va1 = s1p00+k1;	va2 = s1p00+k2;
			}
		}

	  #endif	// 32/64-bit?

	//...and now do 3 radix-64 transforms:

		// Use s1p00-3f for scratch for these 3 DFTs ... since transform length N = odd*64,
		// the leading pow2-shift arg = trailz(N) - trailz(64) = 0:
		jt = j1    ; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)r00,t_offsets, s1p00, a+jt,(dif_o_offsets+0x00));
		jt = j1+p80; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)r40,t_offsets, s1p00, a+jt,(dif_o_offsets+0x40));
		jt = j1+p40; SSE2_RADIX_64_DIF( FALSE, thr_id, 0, (double*)r80,t_offsets, s1p00, a+jt,(dif_o_offsets+0x80));

	#else	/* !USE_SSE2 */

	//...gather the needed data (192 64-bit complex) and do 16 radix-3 transforms:
		tptr = t;
		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head of each trio
		k0 =   0;	k1 = p80;	k2 = p40;
		for(l = 0; l < 64; l++) {
			jp = plo[ilo];	// Same for each elt of the current-offset trio
			jt = j1 + jp; jp += j2;
			RADIX_03_DFT(s,c3m1,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],
				t00,t01,t02,t03,t04,t05,
				tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x80)->re,(tptr+0x80)->im
			);	tptr++;
			ilo -= 3;
			// If low nibble underflows need to add 16 to it, decrement high nibble ihi, and recompute k1,k2.
			// Use branch here because dependent-op count > 10:
			if(ilo < 0) {	// This branch only taken every 16th loop exec
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 12);	// 1st high-part offset
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 2nd offset
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 3rd offset
				k0 = phi[ihi];	k1 = phi[k1];	k2 = phi[k2];	// Use distinct ihi and k0 = phi[ihi] for 1st high-part offset because basic-index ihi must survive from one loop pass to next
			}
		}
	//...and now do 3 radix-64 transforms:
		jt = j1    ;	RADIX_64_DIF((double *)(t+0x00),t_offsets,1, (a+jt),(dif_o_offsets+0x00),RE_IM_STRIDE);	// Inouts in t[00-63]
		jt = j1+p80;	RADIX_64_DIF((double *)(t+0x40),t_offsets,1, (a+jt),(dif_o_offsets+0x40),RE_IM_STRIDE);	// Inouts in t[64-127]
		jt = j1+p40;	RADIX_64_DIF((double *)(t+0x80),t_offsets,1, (a+jt),(dif_o_offsets+0x80),RE_IM_STRIDE);	// Inouts in t[128-191]

	#endif	/* !USE_SSE2 */
	}

	jstart += nwt;
	jhi    += nwt;

	col += RADIX;
	co3 -= RADIX;
}	/* end for(int k=1; k <= khi; k++) */
