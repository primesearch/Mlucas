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
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-992 DIT pass is here:	*/
	#ifdef USE_SSE2

		add0 = &prp_mult;
		#error No SIMD support in this experimntal code as yet!

	#else	/* !USE_SSE2 */

	/*...gather the needed data (992 64-bit complex, i.e. 1984 64-bit reals) and do 31 radix-32 transforms...*/
		for(int kk = 0; kk < 32; ++kk) {
			jj[kk] = ((kk<<5)-kk)<<1;	/* (kk*62) = (kk*31)<<1 */
		}
		iptr = o;	// Name offs-array here 'o' rather than 'q' so can use same length-32 array for DIF and DIT index-offsets
		for(l = 0; l < 31; ++l) {
			idx = (31-l) & (-(l>0));	// Base-offset Index into phi[] for the current row = 00,3c,3a,...,4,2
										// (Since our phi array here is in mults of 32, phi = 0x3c0 ==> idx = 3c/2 = 1e)
			is_odd = (l&1) + (l==0); is_even = (~is_odd)&1;	// 0-row phi-pair relative offsets behave like odd-index row
			// Compute index in the perm16 array needed for current row:
			pidx = 1+((l+1)>>1);	//<*** odd  row (and 0): 1st []-idx = (row+1)/2, 2nd []-idx = 1+(row+1)/2
									//<*** even row: []-idx = 1+(row+1)/2 for both []
			// First 16 offsets:
			jp = phi[idx] + plo[is_even<<4];  // = phi[idx]+p10 for even-idx rows, = phi[idx] for odd]
			i64 = dit_perm16[(pidx-is_odd)&0xf];
			// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];		o[0x0] = jp + k0;
			k1 = plo[(i64 >> 56)&0xf];		o[0x1] = jp + k1;
			k2 = plo[(i64 >> 52)&0xf];		o[0x2] = jp + k2;
			k3 = plo[(i64 >> 48)&0xf];		o[0x3] = jp + k3;
			k4 = plo[(i64 >> 44)&0xf];		o[0x4] = jp + k4;
			k5 = plo[(i64 >> 40)&0xf];		o[0x5] = jp + k5;
			k6 = plo[(i64 >> 36)&0xf];		o[0x6] = jp + k6;
			k7 = plo[(i64 >> 32)&0xf];		o[0x7] = jp + k7;
			k8 = plo[(i64 >> 28)&0xf];		o[0x8] = jp + k8;
			k9 = plo[(i64 >> 24)&0xf];		o[0x9] = jp + k9;
			ka = plo[(i64 >> 20)&0xf];		o[0xa] = jp + ka;
			kb = plo[(i64 >> 16)&0xf];		o[0xb] = jp + kb;
			kc = plo[(i64 >> 12)&0xf];		o[0xc] = jp + kc;
			kd = plo[(i64 >>  8)&0xf];		o[0xd] = jp + kd;
			ke = plo[(i64 >>  4)&0xf];		o[0xe] = jp + ke;
			kf = plo[(i64      )&0xf];		o[0xf] = jp + kf;
			// Second 16 offsets:
			jp = phi[idx] + plo[is_odd<<4];  // = phi[idx] for even-idx rows, = phi[idx]+p10 for odd]
			i64 = dit_perm16[pidx&0xf];	// 2nd []-idx = (row+1)/2 for both even and odd rows
			// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];		o[0x10] = jp + k0;
			k1 = plo[(i64 >> 56)&0xf];		o[0x11] = jp + k1;
			k2 = plo[(i64 >> 52)&0xf];		o[0x12] = jp + k2;
			k3 = plo[(i64 >> 48)&0xf];		o[0x13] = jp + k3;
			k4 = plo[(i64 >> 44)&0xf];		o[0x14] = jp + k4;
			k5 = plo[(i64 >> 40)&0xf];		o[0x15] = jp + k5;
			k6 = plo[(i64 >> 36)&0xf];		o[0x16] = jp + k6;
			k7 = plo[(i64 >> 32)&0xf];		o[0x17] = jp + k7;
			k8 = plo[(i64 >> 28)&0xf];		o[0x18] = jp + k8;
			k9 = plo[(i64 >> 24)&0xf];		o[0x19] = jp + k9;
			ka = plo[(i64 >> 20)&0xf];		o[0x1a] = jp + ka;
			kb = plo[(i64 >> 16)&0xf];		o[0x1b] = jp + kb;
			kc = plo[(i64 >> 12)&0xf];		o[0x1c] = jp + kc;
			kd = plo[(i64 >>  8)&0xf];		o[0x1d] = jp + kd;
			ke = plo[(i64 >>  4)&0xf];		o[0x1e] = jp + ke;
			kf = plo[(i64      )&0xf];		o[0x1f] = jp + kf;
			RADIX_32_DIT(
				a+j1,iptr,RE_IM_STRIDE,
				(double *)t,jj,1
			);
			for(int kk = 0; kk < 32; ++kk)
			{
				jj[kk] += 2;
			}
		}
		//...and now do 32 radix-31 transforms, with output permutation as described in radix992_dit_pass1:
		tptr = t;
		for(int kk = 0; kk < 32; ++kk)
		{
			mask = (-(kk>0));
			lshift = (kk-1) & mask;
			iptr = phi + lshift;
			RADIX_31_DIT(
				(double *)tptr,
				a+j1+plo[(32-kk) & mask], iptr
			);
			tptr += 31;
		}

	#endif	// USE_SSE2?

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 992 separate blocks of the A-array, we need 992 separate carries.	*/

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
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			AVX_cmplx_carry_norm_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			tm1 += 8; tmp += 1; itmp += 4; i = 0;
		}

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1). */

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)
	  // SSE2 build not supported (in the sense of yielding working code) for this function, but need it to compile
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
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
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
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#else	/* Scalar-double mode: */

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
		cmplx_carry_norm_errcheck0(a[j1  ],a[j2   ],*addr,*itmp,0,prp_mult); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p1],a[j2+p1],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p2],a[j2+p2],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p3],a[j2+p3],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		// Remaining 16-tets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];
			cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],*addr,*itmp,l,prp_mult); ++l; ++addr; ++itmp;
		}
		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	/* #ifdef USE_SSE2 */

	}
	else	/* MODULUS_TYPE_FERMAT */
	{

	#ifdef USE_AVX

		/* For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1. */

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

		tmp = base_negacyclic_root;	tm2 = tmp+1;

	  #if HIACC
		/* Hi-accuracy version needs 8 copies of each base root: */
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

	  #else	/* HIACC = false: */

		/* Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
		fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6: */
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

		/* The above need some inits to prepare for the AVX version of the Fermat-mod carry macro: */
		SSE2_fermat_carry_init_loacc(base_negacyclic_root);

	  #endif

		/* AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
		(processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6: */
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
		while(tm0 < x00) {	// Can't use l for loop index here since need it for byte offset in carry macro call
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x780, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic], tm2,p1,p2,p3, add0);
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
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tm1-cy_r)]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x780, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic], tm2,p1,p2,p3, add0);
			tm0 += 8; tm1++;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}


	  #endif	/* HIACC? */

	#elif defined(USE_SSE2)

		/* For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c. */

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

	  #ifdef MULTITHREAD	// Run out of registers here in serial-build mode, so use (threaded or not?) to toggle carry-macro version selection here:

		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],jcycle[ic],icycle[ic+1],jcycle[ic+1] of the non-looped version with icycle[ic],jcycle[ic],icycle[jc],jcycle[jc]:
		ic = 0; jc = 1;
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while((int)(tmp-cy_r) < RADIX) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tmp-cy_r)>>2]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
			tm2 += (-((int)((tmp-cy_r)>>1)&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_fermat_carry_norm_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,icycle[ic],jcycle[ic],icycle[jc],jcycle[jc], tm2,p1, add0);
			tm1 += 4; tmp += 2;
			MOD_ADD32(ic, 2, ODD_RADIX, ic);
			MOD_ADD32(jc, 2, ODD_RADIX, jc);
		}

	  #else // Mar 2014: Worked around the out-of-regs compiler issues with the _X2 version of this macro (the
			// code in carry_gcc64.h has details), but keep non-X2 version in case hit out-of-regs again at some point

		ic = 0;	// ic = idx into [i|j]cycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while((int)(tmp-cy_r) < RADIX) {
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[(int)(tmp-cy_r)>>2]);	// poff[] = p0,4,8,...; (tm1-cy_r) acts as a linear loop index running from 0,...,RADIX-1 here.
			tm2 += p1*((int)(tmp-cy_r)&0x3);	// Added offset cycles among p0,1,2,3
			SSE2_fermat_carry_norm_errcheck(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,icycle[ic],jcycle[ic], tm2, add0);
			tm1 += 2; tmp++;
			MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}

	  #endif

	#else	/* Scalar-double mode: */

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		ntmp = 0; addr = cy_r; addi = cy_i; ic = 0;	// ic = idx into icycle mini-array, gets incremented (mod ODD_RADIX) between macro calls
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];
			fermat_carry_norm_errcheckB(a[jt   ],a[jp   ],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p1],a[jp+p1],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p2],a[jp+p2],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p3],a[jp+p3],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}

		for(n = 0; n < ODD_RADIX; n++) {
			icycle[n] += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
			icycle[n] += ( (-(int)((uint32)icycle[n] >> 31)) & nwt);
		}

	#endif	/* #ifdef USE_SSE2 */

	/* Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only: */
	#ifdef USE_SSE2

		for(n = 0; n < ODD_RADIX; n++) {
			icycle[n] += wts_idx_inc2;		icycle[n] += ( (-(icycle[n] < 0)) & nwt16);
			jcycle[n] += wts_idx_inc2;		jcycle[n] += ( (-(jcycle[n] < 0)) & nwt16);
		#ifdef USE_AVX
			kcycle[n] += wts_idx_inc2;		kcycle[n] += ( (-(kcycle[n] < 0)) & nwt16);
			lcycle[n] += wts_idx_inc2;		lcycle[n] += ( (-(lcycle[n] < 0)) & nwt16);
		#endif
		}
	#endif

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-992 DIF pass is here:	*/

	#ifdef USE_SSE2
		// no sse2!
	#else	/* !USE_SSE2 */

		/*...gather the needed data (992 64-bit complex, i.e 1984 64-bit reals) and do 32 radix-31 transforms...*/
		tptr = t;
		for(l = 0; l < 32; ++l)
		{
			iptr = dif_p20_cperms+l;
			RADIX_31_DIF(
				a+j1+plo[l], iptr,
				(double *)tptr
			);
			tptr += 31;
		}
	//...and now do 31 radix-32 transforms:
		for(int kk = 0; kk < 32; ++kk) {
			jj[kk] = ((kk<<5)-kk)<<1;	// DFT macro takes *real*-double inputs, thus compute doubled offsets k*62
		}
		iptr = o;
		for(l = 0; l < 31; ++l) {
			idx = (31-l) & (-(l>0));	// Base-offset Index into phi[] for the current row = 00,3c,3a,...,4,2
										// (Since our phi array here is in mults of 32, phi = 0x3c0 ==> idx = 3c/2 = 1e)
			is_odd = (l&1) + (l==0); is_even = (~is_odd)&1;	// 0-row phi-pair relative offsets behave like odd-index row
			// Compute index in the perm16 array needed for current row:
			pidx = (l+1)>>1;	//<*** odd  row (and 0): []-idx = (row+1)/2 for both []
								//<*** even row: 1st []-idx = 1+(row+1)/2, 2nd []-idx = (row+1)/2
			// First 16 offsets:
			jp = phi[idx] + plo[is_even<<4];  // = phi[idx]+p10 for even-idx rows, = phi[idx] for odd]
			i64 = dif_perm16[(pidx+is_even)&0xf];
			// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];		o[0x0] = jp + k0;
			k1 = plo[(i64 >> 56)&0xf];		o[0x1] = jp + k1;
			k2 = plo[(i64 >> 52)&0xf];		o[0x2] = jp + k2;
			k3 = plo[(i64 >> 48)&0xf];		o[0x3] = jp + k3;
			k4 = plo[(i64 >> 44)&0xf];		o[0x4] = jp + k4;
			k5 = plo[(i64 >> 40)&0xf];		o[0x5] = jp + k5;
			k6 = plo[(i64 >> 36)&0xf];		o[0x6] = jp + k6;
			k7 = plo[(i64 >> 32)&0xf];		o[0x7] = jp + k7;
			k8 = plo[(i64 >> 28)&0xf];		o[0x8] = jp + k8;
			k9 = plo[(i64 >> 24)&0xf];		o[0x9] = jp + k9;
			ka = plo[(i64 >> 20)&0xf];		o[0xa] = jp + ka;
			kb = plo[(i64 >> 16)&0xf];		o[0xb] = jp + kb;
			kc = plo[(i64 >> 12)&0xf];		o[0xc] = jp + kc;
			kd = plo[(i64 >>  8)&0xf];		o[0xd] = jp + kd;
			ke = plo[(i64 >>  4)&0xf];		o[0xe] = jp + ke;
			kf = plo[(i64      )&0xf];		o[0xf] = jp + kf;
			// Second 16 offsets:
			jp = phi[idx] + plo[is_odd<<4];  // = phi[idx] for even-idx rows, = phi[idx]+p10 for odd]
			i64 = dif_perm16[pidx];	// 2nd []-idx = (row+1)/2 for both even and odd rows
			// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];		o[0x10] = jp + k0;
			k1 = plo[(i64 >> 56)&0xf];		o[0x11] = jp + k1;
			k2 = plo[(i64 >> 52)&0xf];		o[0x12] = jp + k2;
			k3 = plo[(i64 >> 48)&0xf];		o[0x13] = jp + k3;
			k4 = plo[(i64 >> 44)&0xf];		o[0x14] = jp + k4;
			k5 = plo[(i64 >> 40)&0xf];		o[0x15] = jp + k5;
			k6 = plo[(i64 >> 36)&0xf];		o[0x16] = jp + k6;
			k7 = plo[(i64 >> 32)&0xf];		o[0x17] = jp + k7;
			k8 = plo[(i64 >> 28)&0xf];		o[0x18] = jp + k8;
			k9 = plo[(i64 >> 24)&0xf];		o[0x19] = jp + k9;
			ka = plo[(i64 >> 20)&0xf];		o[0x1a] = jp + ka;
			kb = plo[(i64 >> 16)&0xf];		o[0x1b] = jp + kb;
			kc = plo[(i64 >> 12)&0xf];		o[0x1c] = jp + kc;
			kd = plo[(i64 >>  8)&0xf];		o[0x1d] = jp + kd;
			ke = plo[(i64 >>  4)&0xf];		o[0x1e] = jp + ke;
			kf = plo[(i64      )&0xf];		o[0x1f] = jp + kf;
			RADIX_32_DIF(
				(double *)t,jj,   1,
				a+j1,o,RE_IM_STRIDE
			);
			for(int kk = 0; kk < 32; ++kk) {
				jj[kk] += 2;
			}
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

