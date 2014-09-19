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

	/*...The radix-64 DIT pass is here:	*/

#ifdef USE_SSE2

  #if USE_SCALAR_DFT_MACRO

	SSE2_RADIX_64_DIT( FALSE, thr_id,
		(a+j1),dft_offsets,
		// Intermediates base pointer:
		r00,
		// Output pointer: Base ptr of 8 pairs of vec_dbl local-mem:
		s1p00, c_offsets
	);

  #else

	/* NOTE: THE SIMD MACRO HERE RETURNS OUTPUTS 1-7 INDEX-REVERSED W.R.TO THE SCALAR ANALOG -- THIS REVERSAL IS
	REFLECTED IN THE INPUT-ARGS TO THE ENSUING DFTS-WITH-TWIDDLES, I.E. INSTEAD OF THOSE BEING ORDERED 0,4,2,6,1,5,3,7,
	WE HAVE INPUT ORDERING 0,4,6,2,7,3,5,1: */
	tmp = r00;
	for(l = 0; l < 8; l++) {
		add0 = &a[j1] + poff[l+l];	// poff[2*l] = p00,08,...,38
		add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
		SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, tmp, isrt2)
		tmp += 16;
	}

	/******************* AVX debug stuff: *******************/
#if 0
if(full_pass) {
	tmp = r00;	tm2=tmp+1;
	printf("DIT mid-Outputs for j1 = %d:\n",j1);
	for(jp = 0; jp < 2*RADIX; jp += 16) {
		printf("%3x = %20.10e %20.10e\n",(jp>>1)+0,(tmp+jp+0x0)->d0,(tm2+jp+0x0)->d0);
		printf("%3x = %20.10e %20.10e\n",(jp>>1)+1,(tmp+jp+0xE)->d0,(tm2+jp+0xE)->d0);
		printf("%3x = %20.10e %20.10e\n",(jp>>1)+2,(tmp+jp+0xC)->d0,(tm2+jp+0xC)->d0);
		printf("%3x = %20.10e %20.10e\n",(jp>>1)+3,(tmp+jp+0xA)->d0,(tm2+jp+0xA)->d0);
		printf("%3x = %20.10e %20.10e\n",(jp>>1)+4,(tmp+jp+0x8)->d0,(tm2+jp+0x8)->d0);
		printf("%3x = %20.10e %20.10e\n",(jp>>1)+5,(tmp+jp+0x6)->d0,(tm2+jp+0x6)->d0);
		printf("%3x = %20.10e %20.10e\n",(jp>>1)+6,(tmp+jp+0x4)->d0,(tm2+jp+0x4)->d0);
		printf("%3x = %20.10e %20.10e\n",(jp>>1)+7,(tmp+jp+0x2)->d0,(tm2+jp+0x2)->d0);
	}
}
exit(0);
#endif

/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */
// Note: 1st of the 15 sincos args in each call to SSE2_RADIX8_DIT_TWIDDLE_OOP is the basic isrt2 needed for
// radix-8. This is a workaround of GCC's 30-arg limit for inline ASM macros, which proves a royal pain here.
//...and another kludge for the 30-arg limit: put a copy of (vec_dbl)2.0 into the first of each s1p*r outputs:
	VEC_DBL_INIT(s1p01,2.0);
	VEC_DBL_INIT(s1p02,2.0);
	VEC_DBL_INIT(s1p03,2.0);
	VEC_DBL_INIT(s1p04,2.0);
	VEC_DBL_INIT(s1p05,2.0);
	VEC_DBL_INIT(s1p06,2.0);
	VEC_DBL_INIT(s1p07,2.0);

	/* NOTE: THE ABOVE 0-TWIDDLES SIMD MACRO HERE RETURNS OUTPUTS INDEX-REVERSED W.R.TO THE SCALAR ANALOG -- THIS REVERSAL IS
	REFLECTED IN THE INPUT-ARGS TO THE ENSUING DFTS-WITH-TWIDDLES, I.E. INSTEAD OF THOSE BEING ORDERED 0,4,2,6,1,5,3,7,
	WE HAVE INPUT ORDERING 0,4,6,2,7,3,5,1:
	*/
	// Block 0: jt = j1; jp = j2, All unity twiddles:
	SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
		r00,r10,r20,r30,r40,r50,r60,r70,
		s1p00,s1p38,s1p30,s1p28,s1p20,s1p18,s1p10,s1p08, isrt2
	);
	// Block 4: jt = j1 + p04;	jp = j2 + p04;
	SSE2_RADIX8_DIT_TWIDDLE_OOP(
		r08,r18,r28,r38,r48,r58,r68,r78,
		s1p04,s1p24,s1p14,s1p34,s1p0c,s1p2c,s1p1c,s1p3c,
		/* c1,s1	c2,s2		c3,s3		c4,s4		c5,s5	c6,s6	c7,s7 (in terms of roots-naming in this macro) */
		ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
	);
#if 0
tmp = r08;	tm2=tmp+1;
printf("DIT CMUL Outputs::\n");
for(jp = 0; jp < 2*RADIX; jp += 16) {
	if(jp&16)printf("%1x = %20.10e %20.10e\n",jp>>4,(tmp+jp)->d0,(tm2+jp)->d0);
}
#endif
	// Block 2: jt = j1 + p02;	jp = j2 + p02;
	SSE2_RADIX8_DIT_TWIDDLE_OOP(
		r0C,r1C,r2C,r3C,r4C,r5C,r6C,r7C,	// 2/6 swap ==> r*4 / r*C swap
		s1p02,s1p22,s1p12,s1p32,s1p0a,s1p2a,s1p1a,s1p3a,
		isrt2,isrt2, cc4,ss4, ss4,cc4, cc2,ss2, ss6,cc6, cc6,ss6, ss2,cc2
	);
#if 0
tmp = r0C;	tm2=tmp+1;
printf("DIT CMUL Outputs:\n");
for(jp = 0; jp < 2*RADIX; jp += 16) {
	if(jp&16)printf("%1x = %20.10e %20.10e\n",jp>>4,(tmp+jp)->d0,(tm2+jp)->d0);
}
#endif
	// Block 6: jt = j1 + p06;	jp = j2 + p06;
	SSE2_RADIX8_DIT_TWIDDLE_OOP(
		r04,r14,r24,r34,r44,r54,r64,r74,	// 2/6 swap ==> r*4 / r*C swap
		s1p06,s1p26,s1p16,s1p36,s1p0e,s1p2e,s1p1e,s1p3e,
		nisrt2,isrt2, ss4,cc4, ncc4,nss4, cc6,ss6, ncc2,ss2, nss2,cc2, nss6,ncc6
	);
#if 0
tmp = r04;	tm2=tmp+1;
printf("DIT CMUL Outputs:\n");
for(jp = 0; jp < 2*RADIX; jp += 16) {
	if(jp&16)printf("%1x = %20.10e %20.10e\n",jp>>4,(tmp+jp)->d0,(tm2+jp)->d0);
}
#endif
	// Block 1: jt = j1 + p01;	jp = j2 + p01;
	SSE2_RADIX8_DIT_TWIDDLE_OOP(
		r0E,r1E,r2E,r3E,r4E,r5E,r6E,r7E,	// 1/7 swap ==> r*2 / r*E swap
		s1p01,s1p21,s1p11,s1p31,s1p09,s1p29,s1p19,s1p39,
		cc4,ss4, cc2,ss2, cc6,ss6, cc1,ss1, cc5,ss5, cc3,ss3, cc7,ss7
	);
#if 0
tmp = r0E;	tm2=tmp+1;
printf("DIT CMUL Outputs:\n");
for(jp = 0; jp < 2*RADIX; jp += 16) {
	if(jp&16)printf("%1x = %20.10e %20.10e\n",jp>>4,(tmp+jp)->d0,(tm2+jp)->d0);
}
#endif
	// Block 5: jt = j1 + p05;	jp = j2 + p05;
	SSE2_RADIX8_DIT_TWIDDLE_OOP(
		r06,r16,r26,r36,r46,r56,r66,r76,
		s1p05,s1p25,s1p15,s1p35,s1p0d,s1p2d,s1p1d,s1p3d,
		nss4,cc4, ss6,cc6, ncc2,ss2, cc5,ss5, ncc7,ss7, ss1,cc1, ncc3,nss3
	);
#if 0
tmp = r06;	tm2=tmp+1;
printf("DIT CMUL Outputs:\n");
for(jp = 0; jp < 2*RADIX; jp += 16) {
	if(jp&16)printf("%1x = %20.10e %20.10e\n",jp>>4,(tmp+jp)->d0,(tm2+jp)->d0);
}
#endif
	// Block 3: jt = j1 + p03;	jp = j2 + p03;
	SSE2_RADIX8_DIT_TWIDDLE_OOP(
		r0A,r1A,r2A,r3A,r4A,r5A,r6A,r7A,	// 3/5 swap ==> r*6 / r*A swap
		s1p03,s1p23,s1p13,s1p33,s1p0b,s1p2b,s1p1b,s1p3b,
		ss4,cc4, cc6,ss6, nss2,cc2, cc3,ss3, ss1,cc1, ss7,cc7, nss5,cc5
	);
#if 0
tmp = r0A;	tm2=tmp+1;
printf("DIT CMUL Outputs:\n");
for(jp = 0; jp < 2*RADIX; jp += 16) {
	if(jp&16)printf("%1x = %20.10e %20.10e\n",jp>>4,(tmp+jp)->d0,(tm2+jp)->d0);
}
#endif
	// Block 7: jt = j1 + p07;	jp = j2 + p07;
	SSE2_RADIX8_DIT_TWIDDLE_OOP(
		r02,r12,r22,r32,r42,r52,r62,r72,	// 1/7 swap ==> r*2 / r*E swap
		s1p07,s1p27,s1p17,s1p37,s1p0f,s1p2f,s1p1f,s1p3f,
		ncc4,ss4, ss2,cc2, nss6,ncc6, cc7,ss7, ncc3,nss3, nss5,cc5, ss1,ncc1
	);
#if 0
tmp = r02;	tm2=tmp+1;
printf("DIT CMUL Outputs:\n");
for(jp = 0; jp < 2*RADIX; jp += 16) {
	if(jp&16)printf("%1x = %20.10e %20.10e\n",jp>>4,(tmp+jp)->d0,(tm2+jp)->d0);
}
exit(0);
#endif
	/******************* AVX debug stuff: *******************/
#if 0
//Outputs s1p*[a-e] bad ... corr. to 2nd-half outs of Blocks 2,4,6,3,5 - why not Blocks 1,7 ???
//More:
//- 0a/2a and 1a/3a swapped
//- 0e/2e and 1e/3e are duplicates (but no match to ref-outs)
if(full_pass) {
	tmp = s1p00;	tm2=tmp+1;
	printf("DIT Outputs for j1 = %d:\n",j1);
	for(jp = 0; jp < 2*RADIX; jp += 2) {
		printf("%3x = %20.10e %20.10e\n",jp>>1,(tmp+jp)->d0,(tm2+jp)->d0);
	}
}
exit(0);
#endif


  #endif

#else	// USE_SSE2 = False, or non-64-bit-GCC:

/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */

	// Gather the needed data and do 8 twiddleless length-8 subtransforms, with p-offsets in-order:
	tptr = t;
	for(l = 0; l < 8; l++) {
		k2 = poff[l+l];	// poffs[2*l] = p00,08,10,...,30,38
		jt = j1 + k2; jp = j2 + k2;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im
		);
		tptr += 8;
	}
	/******************* AVX debug stuff: *******************/
#if 0
tptr = t;
printf("DIT mid-Outputs for j1 = %d:\n",j1);
for(jp = 0; jp < RADIX; jp++) {
	printf("%3x = %20.10e %20.10e\n",jp,(tptr+jp)->re,(tptr+jp)->im);
}
exit(0);
#endif

/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */

	/* Block 0: */ tptr = t;	jt = j1;	jp = j2;
	/* 0-index block has all-unity twiddles: Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in-order here: */
	RADIX_08_DIT_OOP(
		tptr->re,tptr->im,(tptr+0x08)->re,(tptr+0x08)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x18)->re,(tptr+0x18)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x28)->re,(tptr+0x28)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x38)->re,(tptr+0x38)->im,
		a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38]
	);

	// Remaining 7 sets of macro calls done in loop:
	for(l = 1; l < 8; l++) {
		tptr = t + reverse(l,8);
		k2 = po_br[l];	// po_br[] = p[04261537]
		jt = j1 + k2; jp = j2 + k2;
		addr = DFT64_TWIDDLES[l-1]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_08_DIT_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x08)->re,(tptr+0x08)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x18)->re,(tptr+0x18)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x28)->re,(tptr+0x28)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x38)->re,(tptr+0x38)->im,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c)
		);
//if(l==2)exit(0);
	}
	/******************* AVX debug stuff: *******************/
#if 0
printf("DIT Outputs for j1 = %d:\n",j1);
for(l = 0; l < RADIX>>2; l++) {
	jt = j1 + poff[l];	jp = jt + RE_IM_STRIDE;
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x0,a[jt    ],a[jp    ]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x1,a[jt+p01],a[jp+p01]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x2,a[jt+p02],a[jp+p02]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x3,a[jt+p03],a[jp+p03]);
}
exit(0);
#endif

#endif	// USE_SSE2 ?

/*...Now do the carries. Since the outputs would
normally be getting dispatched to [radix] separate blocks of the A-array, we need [radix] separate carries.	*/

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
		AVX_cmplx_carry_norm_pow2_errcheck0_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		tm1 += 8; tmp++; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			AVX_cmplx_carry_norm_pow2_errcheck1_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			tm1 += 8; tmp++; itmp += 4;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

		l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = ctmp->im = wtl;		++ctmp;
		ctmp->re = ctmp->im = wtn;		++ctmp;
		ctmp->re = ctmp->im = wtlp1;	++ctmp;
		ctmp->re = ctmp->im = wtnm1;

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		SSE2_cmplx_carry_norm_pow2_errcheck0_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_pow2_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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
			SSE2_cmplx_carry_norm_pow2_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

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

		/*...set0 is slightly different from others; divide work into 16 blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		
		l = 0; addr = cy_r; itmp = bjmodn;
	   cmplx_carry_norm_pow2_errcheck0(a[j1    ],a[j2    ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_pow2_errcheck(a[j1+p01],a[j2+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_pow2_errcheck(a[j1+p02],a[j2+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_pow2_errcheck(a[j1+p03],a[j2+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,08,...,60
			cmplx_carry_norm_pow2_errcheck(a[jt    ],a[jp    ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_pow2_errcheck(a[jt+p01],a[jp+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_pow2_errcheck(a[jt+p02],a[jp+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_pow2_errcheck(a[jt+p03],a[jp+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?
	}
	else	/* Fermat-mod carry in SIMD mode */
	{
	#ifdef USE_AVX

		// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.
		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

		tmp = half_arr+2;
		VEC_DBL_INIT(tmp, scale);
		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

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
		tm0 = s1p00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x1000;
		for(i = 0; i < RADIX>>2; i++) {	// RADIX/4 loop passes
			SSE2_fermat_carry_norm_pow2_errcheck_X4(tm0,tmp,l,tm1,tm2,half_arr,sign_mask);
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
	  #if (OS_BITS == 32)
		for(l = 0; l < RADIX; l++) {	// RADIX loop passes
			SSE2_fermat_carry_norm_pow2_errcheck   (tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
			tm1 += 2; tmp++;
		}
	  #else	// 64-bit SSE2
		for(l = 0; l < RADIX>>1; l++) {	// RADIX/2 loop passes
			SSE2_fermat_carry_norm_pow2_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
			tm1 += 4; tmp += 2;
		}
	  #endif

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		ntmp = 0; addr = cy_r; addi = cy_i;
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];	// poff[] = p04,08,...,60
			fermat_carry_norm_pow2_errcheck(a[jt    ],a[jp    ],*addr,*addi,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi;
			fermat_carry_norm_pow2_errcheck(a[jt+p01],a[jp+p01],*addr,*addi,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi;
			fermat_carry_norm_pow2_errcheck(a[jt+p02],a[jp+p02],*addr,*addi,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi;
			fermat_carry_norm_pow2_errcheck(a[jt+p03],a[jp+p03],*addr,*addi,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi;
		}

	#endif	/* #ifdef USE_SSE2 */

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-64 DIF pass is here:	*/

#ifdef USE_SSE2

/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors - we use the same positive-power
roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f].
*/
/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */

  #if USE_SCALAR_DFT_MACRO

	SSE2_RADIX_64_DIF( FALSE, thr_id,
		0,
		// Input pointer: Base ptr of 8 local-mem:
		s1p00, c_offsets,	// If want to use this, need a valid i_offsets array/pointer here!
		// Intermediates base pointer:
		r00,
		// Outputs: Base address plus index offset lo/hi-half arrays:
		(a+j1),dft_offsets
	);

  #else

   #ifdef USE_AVX
	#define OFF1	0x200
	#define OFF2	0x400
	#define OFF3	0x600
	#define OFF4	0x800
	#define OFF5	0xa00
	#define OFF6	0xc00
	#define OFF7	0xe00
   #else
	#define OFF1	0x100
	#define OFF2	0x200
	#define OFF3	0x300
	#define OFF4	0x400
	#define OFF5	0x500
	#define OFF6	0x600
	#define OFF7	0x700
   #endif

#if 0
if(full_pass) {
	tmp = s1p00;	tm2=tmp+1;
	for(l = 0; l < RADIX>>2; l++) {
		jt = j1 + poff[l];	jp = jt + RE_IM_STRIDE;
		(tmp+0x0)->d0 = a[jt    ]; (tm2+0x0)->d0 = a[jp    ];
		(tmp+0x2)->d0 = a[jt+p01]; (tm2+0x2)->d0 = a[jp+p01];
		(tmp+0x4)->d0 = a[jt+p02]; (tm2+0x4)->d0 = a[jp+p02];
		(tmp+0x6)->d0 = a[jt+p03]; (tm2+0x6)->d0 = a[jp+p03];
		tmp += 8;	tm2 += 8;
	}
}
#endif
#if 0
if(full_pass) {
	tmp = s1p00;	tm2=tmp+1;
	printf("DIF Inputs for j1 = %d:\n",j1);
	for(jp = 0; jp < 2*RADIX; jp += 2) {
		printf("%3x = %20.10e %20.10e\n",jp>>1,(tmp+jp)->d0,(tm2+jp)->d0);
	}
}
exit(0);
#endif

   // Unlike for DIT, reducing obj-code size here by wrapping SSE2_RADIX8_DIF_0TWIDDLE in a loop hurt performance on SSE2:
   #if 1
	// Outout (s-pointer) offsets are normal BRed - note r's in this routine are separated as 1*index-stride, not 2* as for others:
	vec_dbl *out0, *out1, *out2, *out3, *out4, *out5, *out6, *out7;
	out0 = r00; out1 = out0+8; out2 = out0+4; out3 = out0+12; out4 = out0+2; out5 = out0+10; out6 = out0+6; out7 = out0+14;
	for(l = 0; l < 8; l++) {
		k1 = reverse(l,8)<<1;
		tmp = s1p00 + k1;
		SSE2_RADIX8_DIF_0TWIDDLE(
			tmp,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
			out0,out1,out2,out3,out4,out5,out6,out7, isrt2
		);
		out0 += 16; out1 += 16; out2 += 16; out3 += 16; out4 += 16; out5 += 16; out6 += 16; out7 += 16;
	}
	/******************* AVX debug stuff: *******************/
#if 0
if(full_pass) {
	tmp = r00;	tm2=tmp+1;
	printf("DIF mid-Outputs for j1 = %d:\n",j1);
	for(jp = 0; jp < 2*RADIX; jp += 2) {
		printf("%3x = %20.10e %20.10e\n",jp>>1,(tmp+jp)->d0,(tm2+jp)->d0);
	}
}
exit(0);
#endif

   #else

	//...Block 0: jt = j1;	jp = j2;
// Relative to RADIX_08_DIF_OOP, this SIMD macro produces outputs in BR order [04261537], so swap r-pointer pairs 2/8,6/C to handle that:
	SSE2_RADIX8_DIF_0TWIDDLE(
		s1p00,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
		r00,r08,r04,r0C,r02,r0A,r06,r0E, isrt2
	);
	//...Block 1: jt = j1 + p04;	jp = j2 + p04;
	SSE2_RADIX8_DIF_0TWIDDLE(
		s1p04,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
		r10,r18,r14,r1C,r12,r1A,r16,r1E, isrt2
	);
	//...Block 2: jt = j1 + p02;	jp = j2 + p02;
	SSE2_RADIX8_DIF_0TWIDDLE(
		s1p02,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
		r20,r28,r24,r2C,r22,r2A,r26,r2E, isrt2
	);
	//...Block 3: jt = j1 + p06;	jp = j2 + p06;
	SSE2_RADIX8_DIF_0TWIDDLE(
		s1p06,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
		r30,r38,r34,r3C,r32,r3A,r36,r3E, isrt2
	);
	//...Block 4: jt = j1 + p01;	jp = j2 + p01;
	SSE2_RADIX8_DIF_0TWIDDLE(
		s1p01,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
		r40,r48,r44,r4C,r42,r4A,r46,r4E, isrt2
	);
	//...Block 5: jt = j1 + p05;	jp = j2 + p05;
	SSE2_RADIX8_DIF_0TWIDDLE(
		s1p05,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
		r50,r58,r54,r5C,r52,r5A,r56,r5E, isrt2
	);
	//...Block 6: jt = j1 + p03;	jp = j2 + p03;
	SSE2_RADIX8_DIF_0TWIDDLE(
		s1p03,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
		r60,r68,r64,r6C,r62,r6A,r66,r6E, isrt2
	);
	//...Block 7: jt = j1 + p07;	jp = j2 + p07;
	SSE2_RADIX8_DIF_0TWIDDLE(
		s1p07,OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,
		r70,r78,r74,r7C,r72,r7A,r76,r7E, isrt2
	);

    #endif	// 0?

/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

	/* Block 0: */
	add0 = &a[j1]      ; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
	so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
	SSE2_RADIX8_DIF_0TWIDDLE(
		r00, OFF4,OFF2,OFF6,OFF1,OFF5,OFF3,OFF7,
		add0,add1,add2,add3,add4,add5,add6,add7, isrt2
	);
	/* Block 4: */
	add0 = &a[j1] + p08; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	SSE2_RADIX8_DIF_TWIDDLE_OOP(
		r08,r48,r28,r68,r18,r58,r38,r78,
		add0,add1,add2,add3,add4,add5,add6,add7,
		ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
	);
	/* Block 2: */
	add0 = &a[j1] + p10; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	SSE2_RADIX8_DIF_TWIDDLE_OOP(
		r04,r44,r24,r64,r14,r54,r34,r74,
		add0,add1,add2,add3,add4,add5,add6,add7,
		isrt2,isrt2, cc4,ss4, ss4,cc4, cc2,ss2, ss6,cc6, cc6,ss6, ss2,cc2
	);
	/* Block 6: */
	add0 = &a[j1] + p18; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	SSE2_RADIX8_DIF_TWIDDLE_OOP(
		r0C,r4C,r2C,r6C,r1C,r5C,r3C,r7C,
		add0,add1,add2,add3,add4,add5,add6,add7,
		nisrt2,isrt2, ss4,cc4, ncc4,nss4, cc6,ss6, ncc2,ss2 ,nss2,cc2, nss6,ncc6
	);
	/* Block 1: */
	add0 = &a[j1] + p20; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	SSE2_RADIX8_DIF_TWIDDLE_OOP(
		r02,r42,r22,r62,r12,r52,r32,r72,
		add0,add1,add2,add3,add4,add5,add6,add7,
		cc4,ss4, cc2,ss2, cc6,ss6, cc1,ss1, cc5,ss5, cc3,ss3, cc7,ss7
	);
	/* Block 5: */
	add0 = &a[j1] + p28; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	SSE2_RADIX8_DIF_TWIDDLE_OOP(
		r0A,r4A,r2A,r6A,r1A,r5A,r3A,r7A,
		add0,add1,add2,add3,add4,add5,add6,add7,
		nss4,cc4, ss6,cc6, ncc2,ss2, cc5,ss5, ncc7,ss7, ss1,cc1, ncc3,nss3
	);
	/* Block 3: */
	add0 = &a[j1] + p30; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	SSE2_RADIX8_DIF_TWIDDLE_OOP(
		r06,r46,r26,r66,r16,r56,r36,r76,
		add0,add1,add2,add3,add4,add5,add6,add7,
		ss4,cc4, cc6,ss6, nss2,cc2, cc3,ss3, ss1,cc1, ss7,cc7, nss5,cc5
	);
	/* Block 7: */
	add0 = &a[j1] + p38; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
	SSE2_RADIX8_DIF_TWIDDLE_OOP(
		r0E,r4E,r2E,r6E,r1E,r5E,r3E,r7E,
		add0,add1,add2,add3,add4,add5,add6,add7,
		ncc4,ss4, ss2,cc2, nss6,ncc6, cc7,ss7, ncc3,nss3, nss5,cc5, ss1,ncc1
	);

	#undef OFF1
	#undef OFF2
	#undef OFF3
	#undef OFF4
	#undef OFF5
	#undef OFF6
	#undef OFF7

  #endif	// USE_SCALAR_DFT_MACRO?

#else	// USE_SSE2 = False:

/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */

	/******************* AVX debug stuff: *******************/
#if 0
printf("DIF Inputs for j1 = %d:\n",j1);
for(l = 0; l < RADIX>>2; l++) {
	jt = j1 + poff[l];	jp = jt + RE_IM_STRIDE;
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x0,a[jt    ],a[jp    ]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x1,a[jt+p01],a[jp+p01]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x2,a[jt+p02],a[jp+p02]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x3,a[jt+p03],a[jp+p03]);
}
exit(0);
#endif

	tptr = t;
	for(l = 0; l < 8; l++) {
		k2 = po_br[l];	// po_br[] = p[04261537]
		jt = j1 + k2; jp = j2 + k2;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im
		);
		tptr += 8;
	}

	/******************* AVX debug stuff: *******************/
#if 0
tptr = t;
printf("DIF mid-Outputs for j1 = %d:\n",j1);
for(jp = 0; jp < RADIX; jp++) {
	printf("%3x = %20.10e %20.10e\n",jp,(tptr+jp)->re,(tptr+jp)->im);
}
exit(0);
#endif

/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

	/* Block 0: */ tptr = t;	jt = j1;	jp = j2;
	/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
	so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
	RADIX_08_DIF_OOP(
		tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x08)->re,(tptr+0x08)->im,(tptr+0x28)->re,(tptr+0x28)->im,(tptr+0x18)->re,(tptr+0x18)->im,(tptr+0x38)->re,(tptr+0x38)->im,
		a[jt],a[jp],a[jt+p04],a[jp+p04],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07]
	);

	// Remaining 7 sets of macro calls done in loop:
	for(l = 1; l < 8; l++) {
		tptr = t + reverse(l,8);
		k2 = poff[l+l];	// poffs[2*l] = p00,08,10,...,30,38
		jt = j1 + k2; jp = j2 + k2;
		addr = DFT64_TWIDDLES[l-1]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_08_DIF_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x08)->re,(tptr+0x08)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x18)->re,(tptr+0x18)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x28)->re,(tptr+0x28)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x38)->re,(tptr+0x38)->im,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c)
		);
	}

#endif	/* #ifdef USE_SSE2 */

	/******************* AVX debug stuff: *******************/
#if 0
printf("DIF Outputs for j1 = %d:\n",j1);
for(l = 0; l < RADIX>>2; l++) {
	jt = j1 + poff[l];	jp = jt + RE_IM_STRIDE;
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x0,a[jt    ],a[jp    ]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x1,a[jt+p01],a[jp+p01]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x2,a[jt+p02],a[jp+p02]);
	printf("%3x = %20.10e %20.10e\n",(l<<2)+0x3,a[jt+p03],a[jp+p03]);
}
exit(0);
#endif

	}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(k=1; k <= khi; k++) */

