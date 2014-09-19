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

		/*
		!...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do a radix-20 DIT transform...
		*/
#ifdef USE_SSE2

		add0 = &a[j1    ];
	//	add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p03;	add3 = add0+p02;

	#if defined(COMPILER_TYPE_MSVC)

		/* Outputs in SSE2 modes are temps 2*5*16 = 10*16 = 0x0a0 bytes apart: */
		__asm	mov	eax, add0	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C */
		__asm	mov	edx, add0
		__asm	mov	esi, p01	/* esi will store power-of-2 multiples of p01 throughout */
		__asm	shl	esi, 3		/* Pointer offset for floating doubles */
		__asm	add edx, esi
		__asm	mov ebx, edx	/* add1 = add0+p01 */
		__asm	add edx, esi
		__asm	mov ecx, edx	/* add3 = add0+p02 */
		__asm	add edx, esi	/* add2 = add0+p03 */
		SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x0a0, 0x140, r00)

	/*	add0,1,2,3 = &a[j1+p04]+p3,2,1,0 */
		__asm	mov	esi, p04
		__asm	shl	esi, 3
		__asm	add	eax, esi	/* &a[j1+p04] */
		__asm	add ebx, esi
		__asm	add ecx, esi
		__asm	add edx, esi
		SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x0a0, 0x140, r02)

	/*	add0,1,2,3 = &a[j1+p08]+p1,0,2,3 */
		__asm	mov	esi, p04
		__asm	shl	esi, 3
		__asm	add	eax, esi	/* &a[j1+p08] */
		__asm	add ebx, esi
		__asm	add ecx, esi
		__asm	add edx, esi
		SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x0a0, 0x140, r04)

	/*	add0,1,2,3 = &a[j1+p12]+p2,3,0,1 */
		__asm	mov	esi, p04
		__asm	shl	esi, 3
		__asm	add	eax, esi	/* &a[j1+p12] */
		__asm	add ebx, esi
		__asm	add ecx, esi
		__asm	add edx, esi
		SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x0a0, 0x140, r06)

	/*	add0,1,2,3 = &a[j1+p16]+p0,1,3,2 */
		__asm	mov	esi, p04
		__asm	shl	esi, 3
		__asm	add	eax, esi	/* &a[j1+p16] */
		__asm	add ebx, esi
		__asm	add ecx, esi
		__asm	add edx, esi
		SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x0a0, 0x140, r08)

		/* Radix-5 DFT uses adjacent temps, i.e. stride = 2*16 bytes: */
		SSE2_RADIX_05_DFT_0TWIDDLE(r00,r02,r04,r06,r08,cc1,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r10,r12,r14,r16,r18,cc1,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r20,r22,r24,r26,r28,cc1,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r)
		SSE2_RADIX_05_DFT_0TWIDDLE(r30,r32,r34,r36,r38,cc1,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r)

	#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		SSE2_RADIX20_DIT_NOTWIDDLE(add0,p01,p04,r00,r10,r20,r30,cc1,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r);

	#endif

#else
						 /*          inputs           */ /*                                      outputs                                      */
		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00,t01,t10,t11,t20,t21,t30,t31,rt,it);
		RADIX_04_DIT(a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],t02,t03,t12,t13,t22,t23,t32,t33,rt,it);
		RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],t04,t05,t14,t15,t24,t25,t34,t35,rt,it);
		RADIX_04_DIT(a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],t06,t07,t16,t17,t26,t27,t36,t37,rt,it);
		RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],t08,t09,t18,t19,t28,t29,t38,t39,rt,it);

	/*...and now do 4 radix-5 transforms...*/
						 /*                                                inputs                                                   */ /*                 outputs                   */
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,rt,it);
		RADIX_05_DFT(uc1,uc2,us1,us2,us3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,rt,it);

#endif

/*...Now do the carries. Since the outputs would
normally be getting dispatched to 20 separate blocks of the A-array, we need 20 separate carries.	*/

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

		AVX_cmplx_carry_norm_errcheck0_X4(s1p00r,add1,add2,add3,cy00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p04r,add1,add2,add3,cy04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p08r,add1,add2,add3,cy08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p12r,add1,add2,add3,cy12,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		AVX_cmplx_carry_norm_errcheck1_X4(s1p16r,add1,add2,add3,cy16,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

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

	  /************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

	  #if defined(COMPILER_TYPE_MSVC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
	   #else
		SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
	   #endif

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
	   #else
		SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
	   #endif

	  #endif

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

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

	  #if defined(COMPILER_TYPE_MSVC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16);
	   #else
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16);
	   #endif

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
	   #else
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
	   #endif

	  #endif

		i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

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

		/*...set0 is slightly different from others:	*/
		cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00   );
		cmplx_carry_norm_errcheck (a1p01r,a1p01i,cy01,bjmodn01,1 );
		cmplx_carry_norm_errcheck (a1p02r,a1p02i,cy02,bjmodn02,2 );
		cmplx_carry_norm_errcheck (a1p03r,a1p03i,cy03,bjmodn03,3 );
		cmplx_carry_norm_errcheck (a1p04r,a1p04i,cy04,bjmodn04,4 );
		cmplx_carry_norm_errcheck (a1p05r,a1p05i,cy05,bjmodn05,5 );
		cmplx_carry_norm_errcheck (a1p06r,a1p06i,cy06,bjmodn06,6 );
		cmplx_carry_norm_errcheck (a1p07r,a1p07i,cy07,bjmodn07,7 );
		cmplx_carry_norm_errcheck (a1p08r,a1p08i,cy08,bjmodn08,8 );
		cmplx_carry_norm_errcheck (a1p09r,a1p09i,cy09,bjmodn09,9 );
		cmplx_carry_norm_errcheck (a1p10r,a1p10i,cy10,bjmodn10,10);
		cmplx_carry_norm_errcheck (a1p11r,a1p11i,cy11,bjmodn11,11);
		cmplx_carry_norm_errcheck (a1p12r,a1p12i,cy12,bjmodn12,12);
		cmplx_carry_norm_errcheck (a1p13r,a1p13i,cy13,bjmodn13,13);
		cmplx_carry_norm_errcheck (a1p14r,a1p14i,cy14,bjmodn14,14);
		cmplx_carry_norm_errcheck (a1p15r,a1p15i,cy15,bjmodn15,15);
		cmplx_carry_norm_errcheck (a1p16r,a1p16i,cy16,bjmodn16,16);
		cmplx_carry_norm_errcheck (a1p17r,a1p17i,cy17,bjmodn17,17);
		cmplx_carry_norm_errcheck (a1p18r,a1p18i,cy18,bjmodn18,18);
		cmplx_carry_norm_errcheck (a1p19r,a1p19i,cy19,bjmodn19,19);

		i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
		co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
			 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

/*...The radix-20 DIF pass is here:	*/

	#ifdef USE_SSE2

	  #if defined(COMPILER_TYPE_MSVC)

		/* Index patterns of s1p-terms here are the same as in DIT DFT: */
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,cc1,r00,r02,r04,r08,r06)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,cc1,r10,r12,r14,r18,r16)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,cc1,r20,r22,r24,r28,r26)
		SSE2_RADIX_05_DFT_0TWIDDLE(s1p05r,s1p01r,s1p17r,s1p13r,s1p09r,cc1,r30,r32,r34,r38,r36)

	//	add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p03;	add3 = add0+p02;
		__asm	mov	eax, add0	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */
		__asm	mov	edx, add0
		__asm	mov	esi, p01	/* esi will store power-of-2 multiples of p01 throughout */
		__asm	shl	esi, 3		/* Pointer offset for floating doubles */
		__asm	add edx, esi
		__asm	mov ebx, edx	/* add0+p01 */
		__asm	add edx, esi
		__asm	mov ecx, edx	/* add0+p02 */
		__asm	add edx, esi	/* add0+p03 */
		SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x0a0, 0x140, eax,ebx,edx,ecx)

	//	add0,1,2,3 = &a[j1+p16]+p0,1,3,2
		__asm	mov	esi, p16
		__asm	shl	esi, 3
		__asm	add	eax, esi// &a[j1+p16]
		__asm	add ebx, esi
		__asm	add ecx, esi
		__asm	add edx, esi
		SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, 0x0a0, 0x140, eax,ebx,edx,ecx)

	//	add0,1,2,3 = &a[j1+p12]+p2,3,0,1
		__asm	mov	esi, p04
		__asm	shl	esi, 3
		__asm	sub	eax, esi// &a[j1+p12]
		__asm	sub ebx, esi
		__asm	sub ecx, esi
		__asm	sub edx, esi
		SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, 0x0a0, 0x140, ecx,edx,eax,ebx)

	//	add0,1,2,3 = &a[j1+p04]+p3,2,1,0
		__asm	mov	esi, p08
		__asm	shl	esi, 3
		__asm	sub	eax, esi// &a[j1+p04]
		__asm	sub ebx, esi
		__asm	sub ecx, esi
		__asm	sub edx, esi
		SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, 0x0a0, 0x140, edx,ecx,ebx,eax)

	//	add0,1,2,3 = &a[j1+p08]+p1,0,2,3
		__asm	mov	esi, p04
		__asm	shl	esi, 3
		__asm	add	eax, esi// &a[j1+p08]
		__asm	add ebx, esi
		__asm	add ecx, esi
		__asm	add edx, esi
		SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, 0x0a0, 0x140, ebx,eax,ecx,edx)

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		SSE2_RADIX20_DIF_NOTWIDDLE(add0,p01,p04,p08,p16,s1p00r,s1p16r,s1p12r,s1p08r,s1p04r,s1p15r,s1p11r,s1p07r,s1p03r,s1p19r,s1p10r,s1p06r,s1p02r,s1p18r,s1p14r,s1p05r,s1p01r,s1p17r,s1p13r,s1p09r,cc1,r00,r10,r20,r30)

	  #endif

	#else	// Scalar-double mode:

	/*...The radix-20 DIF pass is here:	*/
	  #if PFETCH
		addr = &a[j1];
		prefetch_p_doubles(addr);
	  #endif

	/*...gather the needed data (20 64-bit complex, i.e. 40 64-bit reals) and do 4 radix-5 transforms...*/
						 /*                                                inputs                                                   */ /*                 outputs                   */
	  #if PFETCH																																/*[--y3-] [--y4-] <<<<< swap last 2 outputs to undo swap of these in macro */
		RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it,addr,addp,p01,p02);
		RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it,addr,addp,p03,p04);
		RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it,addr,addp,p05,p06);
		RADIX_05_DFT_PFETCH(uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it,addr,addp,p07,p08);
	  #else
		RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p00r,a1p00i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t08,t09,t06,t07,rt,it);
		RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p15r,a1p15i,a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p19r,a1p19i,t10,t11,t12,t13,t14,t15,t18,t19,t16,t17,rt,it);
		RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p18r,a1p18i,a1p14r,a1p14i,t20,t21,t22,t23,t24,t25,t28,t29,t26,t27,rt,it);
		RADIX_05_DFT       (uc1,uc2,us1,us2,us3,a1p05r,a1p05i,a1p01r,a1p01i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,t30,t31,t32,t33,t34,t35,t38,t39,t36,t37,rt,it);
	  #endif

	/*...and now do 5 radix-4 transforms...*/
						 /*          inputs           */ /*                                      outputs                                      */
	  #if PFETCH
		addp = addr+p09;
		prefetch_p_doubles(addp);

		RADIX_04_DIF_PFETCH(t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it,addr,addp,p10,p11);
		RADIX_04_DIF_PFETCH(t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it,addr,addp,p12,p13);
		RADIX_04_DIF_PFETCH(t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it,addr,addp,p14,p15);
		RADIX_04_DIF_PFETCH(t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it,addr,addp,p16,p17);
		RADIX_04_DIF_PFETCH(t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it,addr,addp,p18,p19);
	  #else
		RADIX_04_DIF       (t00,t01,t10,t11,t20,t21,t30,t31,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],rt,it);
		RADIX_04_DIF       (t02,t03,t12,t13,t22,t23,t32,t33,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],rt,it);
		RADIX_04_DIF       (t04,t05,t14,t15,t24,t25,t34,t35,a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
		RADIX_04_DIF       (t06,t07,t16,t17,t26,t27,t36,t37,a[j1+p07],a[j2+p07],a[j1+p06],a[j2+p06],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
		RADIX_04_DIF       (t08,t09,t18,t19,t28,t29,t38,t39,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],rt,it);
	  #endif

	#endif	/* USE_SSE2 */
	}

	jstart += nwt;
	jhi    += nwt;
	col += RADIX;
	co3 -= RADIX;

}	/* end for(k=1; k <= khi; k++) */

