/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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
	for(j = jstart; j < jhi; j += stride)
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-4032 DIT pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (4032 64-bit complex, i.e. 8064 64-bit reals) and do 63 radix-64 transforms...
		tmp = r00;
		for(kk = 0; kk < ODD_RADIX; ++kk) {
			iptr = dif64_oidx_lo + (kk<<6);
			for(l = 0; l < 64; l += 4)
			{	// Since default opt-moe is now -Os, manually unroll this loop 4x:
				io_offsets[l  ] = plo[*(iptr+dit64_iidx_lo[l  ])];
				io_offsets[l+1] = plo[*(iptr+dit64_iidx_lo[l+1])];
				io_offsets[l+2] = plo[*(iptr+dit64_iidx_lo[l+2])];
				io_offsets[l+3] = plo[*(iptr+dit64_iidx_lo[l+3])];
			}
			jt = j1 + phi[dft_phi[kk]];
			SSE2_RADIX_64_DIT( FALSE, thr_id,
				a+jt,io_offsets,
				s1p00,	// local scratch storage
				tmp,t_offsets
			); tmp += 2;
		}
	//...and now do 64 radix-63 transforms. The required output permutation is as follows:
		tmp = r00;
		for(kk = 0; kk < 64; ++kk)
		{
			jt = (64-kk) & (-(kk > 0));
			iptr = dft_p40_cperms + (jt & (-(kk != 1)));
			for(l = 0; l < ODD_RADIX; l++)
			{
				io_offsets[l] = (int)iptr[l]<<7;	// SIMD: Local-mem analog of phi[l] ==> l*64, cast-to-vec-dbl needs another 2x
			}
			SSE2_RADIX_63_DIT( FALSE, thr_id,
				tmp, toff,
				s1p00 + (jt<<1), io_offsets
			);
			tmp += (ODD_RADIX<<1);
		}

	#else	/* !USE_SSE2 */

	//...gather the needed data (4032 64-bit complex, i.e. 8032 64-bit reals) and do 63 radix-64 transforms...
		tptr = t;
		for(kk = 0; kk < ODD_RADIX; ++kk) {
			iptr = dif64_oidx_lo + (kk<<6);
			for(l = 0; l < 64; l += 4)
			{	// Since default opt-moe is now -Os, manually unroll this loop 4x:
				io_offsets[l  ] = plo[*(iptr+dit64_iidx_lo[l  ])];
				io_offsets[l+1] = plo[*(iptr+dit64_iidx_lo[l+1])];
				io_offsets[l+2] = plo[*(iptr+dit64_iidx_lo[l+2])];
				io_offsets[l+3] = plo[*(iptr+dit64_iidx_lo[l+3])];
			}
			jt = j1 + phi[dft_phi[kk]];
			RADIX_64_DIT(a+jt,io_offsets,RE_IM_STRIDE, (double *)tptr,t_offsets,1);
			tptr++;
		}
	//...and now do 64 radix-63 transforms. The required output permutation is as follows:
		tptr = t;
		for(kk = 0; kk < 64; ++kk)
		{
			jt = (64-kk) & (-(kk > 0));
			iptr = dft_p40_cperms + (jt & (-(kk != 1)));
			for(l = 0; l < ODD_RADIX; l++)
			{
				io_offsets[l] = phi[iptr[l]];
			}
			RADIX_63_DIT(
				(double *)tptr, toff, 1,
				a+j1+plo[jt], io_offsets, RE_IM_STRIDE
			);
			tptr += ODD_RADIX;
		}

	#endif	// USE_SSE2?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need 28 separate carries.	*/

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
		AVX_cmplx_carry_norm_errcheck0_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		tm1 += 8; tmp += 1; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			AVX_cmplx_carry_norm_errcheck1_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			tm1 += 8; tmp += 1; itmp += 4;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

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
		SSE2_cmplx_carry_norm_errcheck0_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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

		/*...set0 is slightly different from others; divide work into blocks of RADIX/4 macro calls, 1st set of which gets pulled out of loop: */		
		l = 0; addr = cy_r; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1   ],a[j2   ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p1],a[j2+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p2],a[j2+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p3],a[j2+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];
			cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

	}
	else	/* MODULUS_TYPE_FERMAT */
	{

	#ifdef USE_AVX
		int k3,k4,k5,k6,k7;
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

	  #else	// HIACC = false:

		// Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
		// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
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

		// The above need some inits to prepare for the AVX version of the Fermat-mod carry macro:
		SSE2_fermat_carry_init_loacc(base_negacyclic_root);

	  #endif

		// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indices j,j+2,j+4,j+6:
	  #if HIACC

		/* The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
		which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls: */

		tm0 = s1p00; tmp = base_negacyclic_root; l = 0x3f000;
		tm1 = cy_r; // tm2 = cy_i;	*** replace with literal-byte-offset in macro call to save a reg
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic = 0; jc = 1; kc = 2; lc = 3;
		while(tm0 < cy_r)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x7e00, 0x1f80,0x3f00,0x5e80, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
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
		for(l = 0; l < RADIX>>2; l++)	// RADIX/4 loop passes
		{
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x7e00, 0x1f80,0x3f00,0x5e80, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
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

	  #if (OS_BITS == 64)

		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],jcycle[ic],icycle[ic+1],jcycle[ic+1] of the non-looped version with icycle[ic],jcycle[ic],icycle[jc],jcycle[jc]:
		ic = 0; jc = 1;
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while(tm1 < cy_r) {
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];
			k2 = jcycle[ic];
			int k3 = icycle[jc];
			int k4 = jcycle[jc];
			SSE2_fermat_carry_norm_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2,k3,k4);
			tm1 += 4; tmp += 2;
			MOD_ADD32(ic, 2, ODD_RADIX, ic);
			MOD_ADD32(jc, 2, ODD_RADIX, jc);
		}

	  #else // Mar 2014: Worked around the out-of-regs compiler issues with the _X2 version of this macro (the
			// code in carry_gcc64.h has details), but keep non-X2 version in case hit out-of-regs again at some point

		ic = 0;	// ic = idx into [i|j]cycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX << 4;	// 32-bit version needs preshifted << 4 input value
		while(tm1 < cy_r) {
			//Sep 2014: Even with reduced-register version of the 32-bit Fermat-mod carry macro,
			// GCC runs out of registers on this one, without some playing-around-with-alternate code-sequences ...
			// Pulling the array-refs out of the carry-macro call like so solves the problem:
			k1 = icycle[ic];
			k2 = jcycle[ic];
			SSE2_fermat_carry_norm_errcheck(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2);
			tm1 += 2; tmp++;
			MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}

	  #endif

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		ntmp = 0; addr = cy_r; addi = cy_i; ic = 0;	// ic = idx into icycle mini-array, gets incremented (mod ODD_RADIX) between macro calls
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];
			fermat_carry_norm_errcheckB(a[jt   ],a[jp   ],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p1],a[jp+p1],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p2],a[jp+p2],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p3],a[jp+p3],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
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
		}
	#endif

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-4032 DIF pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (4032 64-bit complex, i.e 8064 64-bit reals) and do 64 radix-63 transforms...
		tmp = r00;
		for(kk = 0; kk < 64; ++kk)
		{
			iptr = dft_p40_cperms + kk;
			for(l = 0; l < ODD_RADIX; l += 4)
			{	// Since default opt-moe is now -Os, manually unroll this loop 4x:
				io_offsets[l  ] = (int)iptr[l  ]<<7;	// SIMD: Local-mem analog of phi[l] ==> l*64, cast-to-vec-dbl needs another 2x
				io_offsets[l+1] = (int)iptr[l+1]<<7;
				io_offsets[l+2] = (int)iptr[l+2]<<7;
				io_offsets[l+3] = (int)iptr[l+3]<<7;	// Don't care about 1-extra-init here, since we dimensioned io_offsets to 64
			}
			SSE2_RADIX_63_DIF( FALSE, thr_id,
				s1p00 + (kk<<1), io_offsets,
				tmp, toff
			);
			tmp += (ODD_RADIX<<1);
		}
	//...and now do 63 radix-64 transforms:
		tmp = r00;
		for(kk = 0; kk < ODD_RADIX; ++kk) {
			iptr = dif64_oidx_lo + (kk<<6);
			for(l = 0; l < 64; l += 4)
			{	// Since default opt-moe is now -Os, manually unroll this loop 4x:
				io_offsets[l  ] = plo[*iptr++];
				io_offsets[l+1] = plo[*iptr++];
				io_offsets[l+2] = plo[*iptr++];
				io_offsets[l+3] = plo[*iptr++];
			}
			jt = j1 + phi[dft_phi[kk]];
			// Use s1p00-3f for scratch for this 3 DFT ... since transform length N = odd*64,
			// the leading pow2-shift arg = trailz(N) - trailz(64) = 0:
			SSE2_RADIX_64_DIF( FALSE, thr_id,
				0,
				tmp,t_offsets,
				s1p00,	// tmp-storage
				a+jt,io_offsets
			); tmp += 2;
		}

	#else	/* !USE_SSE2 */

	//...gather the needed data (4032 64-bit complex, i.e 8064 64-bit reals) and do 64 radix-63 transforms...
		tptr = t;
		for(kk = 0; kk < 64; ++kk)
		{
			iptr = dft_p40_cperms + kk;
			for(l = 0; l < ODD_RADIX; l += 4)
			{	// Since default opt-moe is now -Os, manually unroll this loop 4x:
				io_offsets[l  ] = phi[iptr[l  ]];
				io_offsets[l+1] = phi[iptr[l+1]];
				io_offsets[l+2] = phi[iptr[l+2]];
				io_offsets[l+3] = phi[iptr[l+3]];	// Don't care about 1-extra-init here, since we dimensioned io_offsets to 64
			}
			RADIX_63_DIF(
				a+j1+plo[kk], io_offsets, RE_IM_STRIDE,
				(double *)tptr, toff, 1
			);
			tptr += ODD_RADIX;
		}
	//...and now do 63 radix-64 transforms:
		tptr = t;
		for(kk = 0; kk < ODD_RADIX; ++kk) {
			iptr = dif64_oidx_lo + (kk<<6);
			for(l = 0; l < 64; l += 4)
			{	// Since default opt-moe is now -Os, manually unroll this loop 4x:
				io_offsets[l  ] = plo[*iptr++];
				io_offsets[l+1] = plo[*iptr++];
				io_offsets[l+2] = plo[*iptr++];
				io_offsets[l+3] = plo[*iptr++];
			}
			jt = j1 + phi[dft_phi[kk]];
			RADIX_64_DIF((double *)tptr,t_offsets,1, (a+jt),io_offsets,RE_IM_STRIDE);
			tptr++;
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

