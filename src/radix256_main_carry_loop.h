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

/*...The radix-256 DIT pass is here:	*/

#ifdef USE_SSE2

	/* Gather the needed data (256 64-bit complex, i.e. 512 64-bit reals) and do 8 twiddleless length-16 subtransforms: */

 #if USE_SCALAR_DFT_MACRO
  #ifdef USE_ARM_V8_SIMD
	#error Unsupported in ARM v8 SIMD build mode!
  #elif defined(USE_AVX512)
	#error Unsupported in AVX-512 build mode!
  #endif
	SSE2_RADIX256_DIT(
		// Inputs: Base address plus 30 index offsets:
		(a+j1),
		dft_offsets_lo,(uint32)0,dft_offsets_hi,
		// Intermediates base pointer:
		r00,
		// Pointers to base-roots data and first of 16 twiddle vectors:
		isrt2,two, twid0,
		// Output pointer: Base ptr of 16 local-mem:
		s1p00
	);

 #else

  #ifdef USE_ARM_V8_SIMD
	uint32 OFF1,OFF2,OFF3,OFF4;
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
	uint32 *ui32_ptr = &(po_lin[0]);	// Kludge - with just po_lin as macro arg, gcc was feeding the low 8 bytes of po_lin
	// po_lin[1,0]as arg to the macro, rather than po_lin-as-pointer. Can't feed the latter directly as arg to macro, because that causes gcc to
	// emit "error: memory input 5 is not directly addressable" error.
	tmp = r00;
	for(l = 0; l < 16; l++) {
		addr = &a[j1] + poff[l<<2];	// poff[4*i] = p10,p20,...,pf0
		SSE2_RADIX16_DIT_0TWIDDLE(
			addr,ui32_ptr,
			isrt2,two,
			tmp,OFF1,OFF2,OFF3,OFF4
		); tmp += 32;
	}

  #ifndef USE_ARM_V8_SIMD
	#undef OFF1
	#undef OFF2
	#undef OFF3
	#undef OFF4
  #endif

/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors - we use the same positive-power
roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f].
*/
  #ifdef USE_ARM_V8_SIMD
	OFF1 = 0x200;
	OFF2 = 0x400;
	OFF3 = 0x600;
	OFF4 = 0x800;
  #elif defined(USE_AVX512)
	#define OFF1	4*0x200
	#define OFF2	4*0x400
	#define OFF3	4*0x600
	#define OFF4	4*0x800
  #elif defined(USE_AVX)
	#define OFF1	2*0x200
	#define OFF2	2*0x400
	#define OFF3	2*0x600
	#define OFF4	2*0x800
  #else
	#define OFF1	0x200
	#define OFF2	0x400
	#define OFF3	0x600
	#define OFF4	0x800
  #endif
// Block 0: All unity twiddles:
	// Jan 2020: For reduced-#args version of macro, need base-pointer r00 and array of 16
	// double-array-equivalent index offsets corr. to the address strides between r00 and r10,r20...,
	// which are 32 vec_dbl = 32*RE_IM_STRIDE doubles apart:
	const int dd = RE_IM_STRIDE<<5, offs[16] = {0,dd,2*dd,3*dd,4*dd,5*dd,6*dd,7*dd,8*dd,9*dd,10*dd,11*dd,12*dd,13*dd,14*dd,15*dd};
	const int *off_ptr = &(offs[0]);
	SSE2_RADIX16_DIT_0TWIDDLE(
		r00,off_ptr,
		isrt2,two,
		s1p00,OFF1,OFF2,OFF3,OFF4
	);

  #if defined(USE_AVX2) && !defined(USE_IMCI512)

/*** Only last 14 of the 15 with-twiddles DFTs allow use of FMA-based macros under Intel AVX2/FMA3: ***/
// Block 8: BR twiddles = {  I.{},  C^ 8,-~C^ 8,  C^ 4,*~C^ 4, *C^ 4,-~C^ 4,  C^ 2,*~C^ 2, *C^ 6,-~C^ 6,  C^ 6,*~C^ 6, *C^ 2,-~C^ 2}
	//jt = j1 + p08;	jp = j2 + p08;
	SSE2_RADIX16_DIT_TWIDDLE_OOP(
		r08,OFF1,OFF2,OFF3,OFF4, s1p08,OFF1,OFF2,OFF3,OFF4, isrt2, twid8
	);

	// Remaining 14 sets of macro calls done in loop:
	for(l = 2; l < 16; l++) {
		k1 = reverse(l,4)<<1;
		tm1 = r00 + k1; tm2 = s1p00 + k1; tmp = twid0 + (k1<<4)-k1;	// Twid-offsets are multiples of 30 vec_dbl
		SSE2_RADIX16_DIT_FMA_OOP(
			tm1,OFF1,OFF2,OFF3,OFF4, tm2,OFF1,OFF2,OFF3,OFF4, tmp
		);
	}

  #else	// Non-FMA version:

	// Remaining 15 sets of macro calls done in loop:
	for(l = 1; l < 16; l++) {
		k1 = reverse(l,4)<<1;
		tm1 = r00 + k1; tm2 = s1p00 + k1; tmp = twid0 + (k1<<4)-k1;	// Twid-offsets are multiples of 30 vec_dbl
		SSE2_RADIX16_DIT_TWIDDLE_OOP(
			tm1,OFF1,OFF2,OFF3,OFF4, tm2,OFF1,OFF2,OFF3,OFF4, isrt2, tmp
		);
	}

  #endif	// FMA/AVX2 ?

  #ifndef USE_ARM_V8_SIMD
	#undef OFF1
	#undef OFF2
	#undef OFF3
	#undef OFF4
  #endif

 #endif	// USE_SCALAR_DFT_MACRO ?

#else	// USE_SSE2 = False, or non-64-bit-GCC:

 #if USE_SCALAR_DFT_MACRO

	RADIX_256_DIT(
		(a+j1),RE_IM_STRIDE,dft_offsets_lo,(uint32)0,dft_offsets_hi,
		(a+j1),RE_IM_STRIDE,dft_offsets_lo,dft_offsets_hi
	);

 #else

	// Gather the needed data and do 16 twiddleless length-16 subtransforms, with p-offsets in-order:
	tptr = t;
	for(ntmp = 0; ntmp < 16; ntmp++) {
		jt = j1 + poff[ntmp<<2]; jp = j2 + poff[ntmp<<2];	// poff[4*ntmp] = p10,p20,...,pf0
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
			c16,s16	// These = (c16,s16) typically def'd for use in the radix-16 DFT
		);	tptr += 16;
	}

/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors - we use the same positive-power
roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f].
*/
// Block 0: All unity twiddles:
	tptr = t; jt = j1;	jp = j2;
	RADIX_16_DIT(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c16,s16
	);

/*** For the last 14 of the 15 with-twiddles DFTs to test with-FMA option in prep. for porting to Intel AVX2/FMA3: ***/
  #ifdef USE_FMA	// Must define - or not - @compile time

/* Can't use FMA-based DFT for first set of non-unity twiddles - the leading [0+i] twiddle leads to div-by-0 in our tangent-based FMA-DFT scheme:
// Block 8: BR twiddles = {  I.{},  C^ 8,-~C^ 8,  C^ 4,*~C^ 4, *C^ 4,-~C^ 4,  C^ 2,*~C^ 2, *C^ 6,-~C^ 6,  C^ 6,*~C^ 6, *C^ 2,-~C^ 2}
	tptr = t + 8; jt = j1 + p08;	jp = j2 + p08;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c1_1,s1_1,c1_2,s1_2,c1_3,s1_3,c1_4,s1_4,c1_5,s1_5,c1_6,s1_6,c1_7,s1_7,c1_8,s1_8,c1_9,s1_9,c1_10,s1_10,c1_11,s1_11,c1_12,s1_12,c1_13,s1_13,c1_14,s1_14,c1_15,s1_15,c1_1_c,tan,c1_1i2,c1_2i2
	);
*/
// Block 8: BR twiddles = {  I.{},  C^ 8,-~C^ 8,  C^ 4,*~C^ 4, *C^ 4,-~C^ 4,  C^ 2,*~C^ 2, *C^ 6,-~C^ 6,  C^ 6,*~C^ 6, *C^ 2,-~C^ 2}
	tptr = t + 8; jt = j1 + p08;	jp = j2 + p08;
	RADIX_16_DIT_TWIDDLE_OOP(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		0,1, ISRT2,ISRT2, -ISRT2,ISRT2, c16,s16, -s16,c16, s16,c16, -c16,s16, c32_1,s32_1, -s32_1,c32_1, s32_3,c32_3, -c32_3,s32_3, c32_3,s32_3, -s32_3,c32_3, s32_1,c32_1, -c32_1,s32_1,
		c16,s16
	);
// Block 4: BR twiddles = {  C^ 8,  C^ 4, *C^ 4,  C^ 2, *C^ 6,  C^ 6, *C^ 2,  C^ 1, *C^ 7,  C^ 5, *C^ 3,  C^ 3, *C^ 5,  C^ 7, *C^ 1}
	tptr = t + 4; jt = j1 + p04;	jp = j2 + p04;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c2_1,s2_1,c2_2,s2_2,c2_3,s2_3,c2_4,s2_4,c2_5,s2_5,c2_6,s2_6,c2_7,s2_7,c2_8,s2_8,c2_9,s2_9,c2_10,s2_10,c2_11,s2_11,c2_12,s2_12,c2_13,s2_13,c2_14,s2_14,c2_15,s2_15,c2_1_c,tan,c2_1i2,c2_2i2
	);
// Block c: BR twiddles = {-~C^ 8, *C^ 4, -C^ 4,  C^ 6,-~C^ 2,*~C^ 2,-*C^ 6,  C^ 3,-~C^ 5, *C^ 1, -C^ 7, *C^ 7, -C^ 1,*~C^ 5,-*C^ 3}
	tptr = t + 12; jt = j1 + p0c;	jp = j2 + p0c;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c3_1,s3_1,c3_2,s3_2,c3_3,s3_3,c3_4,s3_4,c3_5,s3_5,c3_6,s3_6,c3_7,s3_7,c3_8,s3_8,c3_9,s3_9,c3_10,s3_10,c3_11,s3_11,c3_12,s3_12,c3_13,s3_13,c3_14,s3_14,c3_15,s3_15,c3_1_c,tan,c3_1i2,c3_2i2
	);
// Block 2: BR twiddles = {  C^ 4,  C^ 2,  C^ 6,  C^ 1,  C^ 5,  C^ 3,  C^ 7,  D^ 1,  D^ 9,  D^ 5,  D^ d,  D^ 3,  D^ b,  D^ 7,  D^ f}
	tptr = t + 2; jt = j1 + p02;	jp = j2 + p02;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c4_1,s4_1,c4_2,s4_2,c4_3,s4_3,c4_4,s4_4,c4_5,s4_5,c4_6,s4_6,c4_7,s4_7,c4_8,s4_8,c4_9,s4_9,c4_10,s4_10,c4_11,s4_11,c4_12,s4_12,c4_13,s4_13,c4_14,s4_14,c4_15,s4_15,c4_1_c,tan,c4_1i2,c4_2i2
	);
// Block a: BR twiddles = {*~C^ 4, *C^ 6,-~C^ 2,  C^ 5,-~C^ 7, *C^ 1, -C^ 3,  D^ 5,*~D^ d, *D^ 7, -D^ 1,  D^ f,-~D^ 9,*~D^ 3, -D^ b}
	tptr = t + 10; jt = j1 + p0a;	jp = j2 + p0a;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c5_1,s5_1,c5_2,s5_2,c5_3,s5_3,c5_4,s5_4,c5_5,s5_5,c5_6,s5_6,c5_7,s5_7,c5_8,s5_8,c5_9,s5_9,c5_10,s5_10,c5_11,s5_11,c5_12,s5_12,c5_13,s5_13,c5_14,s5_14,c5_15,s5_15,c5_1_c,tan,c5_1i2,c5_2i2
	);
// Block 6: BR twiddles = { *C^ 4,  C^ 6,*~C^ 2,  C^ 3, *C^ 1, *C^ 7,*~C^ 5,  D^ 3, *D^ 5,  D^ f,*~D^ 7,  D^ 9,*~D^ 1, *D^ b,*~D^ d}
	tptr = t + 6; jt = j1 + p06;	jp = j2 + p06;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c6_1,s6_1,c6_2,s6_2,c6_3,s6_3,c6_4,s6_4,c6_5,s6_5,c6_6,s6_6,c6_7,s6_7,c6_8,s6_8,c6_9,s6_9,c6_10,s6_10,c6_11,s6_11,c6_12,s6_12,c6_13,s6_13,c6_14,s6_14,c6_15,s6_15,c6_1_c,tan,c6_1i2,c6_2i2
	);
// Block e: BR twiddles = {-~C^ 4, *C^ 2,-*C^ 6,  C^ 7, -C^ 3,*~C^ 5,~*C^ 1,  D^ 7,-~D^ 1,*~D^ 3,-*D^ 5, *D^ b, -D^ d,-~D^ f,~*D^ 9}
	tptr = t + 14; jt = j1 + p0e;	jp = j2 + p0e;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c7_1,s7_1,c7_2,s7_2,c7_3,s7_3,c7_4,s7_4,c7_5,s7_5,c7_6,s7_6,c7_7,s7_7,c7_8,s7_8,c7_9,s7_9,c7_10,s7_10,c7_11,s7_11,c7_12,s7_12,c7_13,s7_13,c7_14,s7_14,c7_15,s7_15,c7_1_c,tan,c7_1i2,c7_2i2
	);
// Block 1: BR twiddles = {  C^ 2,  C^ 1,  C^ 3,  D^ 1,  D^ 5,  D^ 3,  D^ 7,  E^ 1,  E^ 9,  E^ 5,  E^ d,  E^ 3,  E^ b,  E^ 7,  E^ f}
	tptr = t + 1; jt = j1 + p01;	jp = j2 + p01;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c8_1,s8_1,c8_2,s8_2,c8_3,s8_3,c8_4,s8_4,c8_5,s8_5,c8_6,s8_6,c8_7,s8_7,c8_8,s8_8,c8_9,s8_9,c8_10,s8_10,c8_11,s8_11,c8_12,s8_12,c8_13,s8_13,c8_14,s8_14,c8_15,s8_15,c8_1_c,tan,c8_1i2,c8_2i2
	);
// Block 9: BR twiddles = {*~C^ 2, *C^ 7,-~C^ 5,  D^ 9,*~D^ d, *D^ 5,-~D^ 1,  E^ 9,*~E^ h, *E^ j,-~E^ b,  E^ r,-~E^ t, *E^ 1, -E^ 7}
	tptr = t + 9; jt = j1 + p09;	jp = j2 + p09;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		c9_1,s9_1,c9_2,s9_2,c9_3,s9_3,c9_4,s9_4,c9_5,s9_5,c9_6,s9_6,c9_7,s9_7,c9_8,s9_8,c9_9,s9_9,c9_10,s9_10,c9_11,s9_11,c9_12,s9_12,c9_13,s9_13,c9_14,s9_14,c9_15,s9_15,c9_1_c,tan,c9_1i2,c9_2i2
	);
// Block 5: BR twiddles = { *C^ 6,  C^ 5, *C^ 1,  D^ 5, *D^ 7,  D^ f,*~D^ 3,  E^ 5, *E^ j,  E^ p,*~E^ 1,  E^ f, *E^ 9, *E^ t,*~E^ b}
	tptr = t + 5; jt = j1 + p05;	jp = j2 + p05;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		ca_1,sa_1,ca_2,sa_2,ca_3,sa_3,ca_4,sa_4,ca_5,sa_5,ca_6,sa_6,ca_7,sa_7,ca_8,sa_8,ca_9,sa_9,ca_10,sa_10,ca_11,sa_11,ca_12,sa_12,ca_13,sa_13,ca_14,sa_14,ca_15,sa_15,ca_1_c,tan,ca_1i2,ca_2i2
	);
// Block d: BR twiddles = {-~C^ 6, *C^ 3, -C^ 7,  D^ d, -D^ 1,*~D^ 7,-*D^ 5,  E^ d,-~E^ b,*~E^ 1,-*E^ n, *E^ p, -E^ f,*~E^ r,~*E^ 3}
	tptr = t + 13; jt = j1 + p0d;	jp = j2 + p0d;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		cb_1,sb_1,cb_2,sb_2,cb_3,sb_3,cb_4,sb_4,cb_5,sb_5,cb_6,sb_6,cb_7,sb_7,cb_8,sb_8,cb_9,sb_9,cb_10,sb_10,cb_11,sb_11,cb_12,sb_12,cb_13,sb_13,cb_14,sb_14,cb_15,sb_15,cb_1_c,tan,cb_1i2,cb_2i2
	);
// Block 3: BR twiddles = {  C^ 6,  C^ 3, *C^ 7,  D^ 3,  D^ f,  D^ 9, *D^ b,  E^ 3,  E^ r,  E^ f, *E^ p,  E^ 9, *E^ v,  E^ l, *E^ j}
	tptr = t + 3; jt = j1 + p03;	jp = j2 + p03;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		cc_1,sc_1,cc_2,sc_2,cc_3,sc_3,cc_4,sc_4,cc_5,sc_5,cc_6,sc_6,cc_7,sc_7,cc_8,sc_8,cc_9,sc_9,cc_10,sc_10,cc_11,sc_11,cc_12,sc_12,cc_13,sc_13,cc_14,sc_14,cc_15,sc_15,cc_1_c,tan,cc_1i2,cc_2i2
	);
// Block b: BR twiddles = {*~C^ 6, *C^ 5, -C^ 1,  D^ b,-~D^ 9,*~D^ 1, -D^ d,  E^ b,-~E^ t, *E^ 9, -E^ f, *E^ v,-~E^ 7,*~E^ d,-*E^ r}
	tptr = t + 11; jt = j1 + p0b;	jp = j2 + p0b;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		cd_1,sd_1,cd_2,sd_2,cd_3,sd_3,cd_4,sd_4,cd_5,sd_5,cd_6,sd_6,cd_7,sd_7,cd_8,sd_8,cd_9,sd_9,cd_10,sd_10,cd_11,sd_11,cd_12,sd_12,cd_13,sd_13,cd_14,sd_14,cd_15,sd_15,cd_1_c,tan,cd_1i2,cd_2i2
	);
// Block 7: BR twiddles = { *C^ 2,  C^ 7,*~C^ 5,  D^ 7,*~D^ 3, *D^ b,-~D^ f,  E^ 7, *E^ 1, *E^ t,*~E^ r,  E^ l,*~E^ d, *E^ f,-~E^ n}
	tptr = t + 7; jt = j1 + p07;	jp = j2 + p07;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		ce_1,se_1,ce_2,se_2,ce_3,se_3,ce_4,se_4,ce_5,se_5,ce_6,se_6,ce_7,se_7,ce_8,se_8,ce_9,se_9,ce_10,se_10,ce_11,se_11,ce_12,se_12,ce_13,se_13,ce_14,se_14,ce_15,se_15,ce_1_c,tan,ce_1i2,ce_2i2
	);
// Block f: BR twiddles = {-~C^ 2, *C^ 1,-*C^ 3,  D^ f, -D^ b,*~D^ d,~*D^ 9,  E^ f, -E^ 7,*~E^ b,~*E^ 3, *E^ j,-*E^ r,-~E^ n, ~E^ v}
	tptr = t + 15; jt = j1 + p0f;	jp = j2 + p0f;
	RADIX_16_DIT_FMA_PRETWIDDLE(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
		cf_1,sf_1,cf_2,sf_2,cf_3,sf_3,cf_4,sf_4,cf_5,sf_5,cf_6,sf_6,cf_7,sf_7,cf_8,sf_8,cf_9,sf_9,cf_10,sf_10,cf_11,sf_11,cf_12,sf_12,cf_13,sf_13,cf_14,sf_14,cf_15,sf_15,cf_1_c,tan,cf_1i2,cf_2i2
	);

  #else

	// Remaining 15 sets of macro calls done in loop:
	for(ntmp = 1; ntmp < 16; ntmp++) {
		tptr = t + reverse(ntmp,4);
		jt = j1 + po_br[ntmp]; jp = j2 + po_br[ntmp];	// po_br[] = p[084c2a6e195d3b7f]
		const double *addr = DFT256_TWIDDLES[ntmp], *addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_16_DIT_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
			a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
			c16,s16
		);
	}

  #endif // USE_FMA ?

 #endif	// USE_SCALAR_DFT_MACRO ?

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
		// incr must divide nloop [RADIX/8 = 32 or RADIX/16 = 64, depending on whether we use 8-or-16-way carry macros]!
	  #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
	  #else
		const uint32 nloop = RADIX>>3;
	  #endif
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tmp = s1p00; tm1 = cy_r; itmp = bjmodn;
	  #ifndef USE_AVX512
		tm2 = cy_r+1; itm2 = bjmodn+4;	// tm2,itm2 not used in AVX-512 mode
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
		i = (!j);
		addr = &prp_mult;
		tm1 = s1p00; tmp = cy_r; itmp = bjmodn;
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
		tm0 = s1p00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x8000;
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
		tm0 = s1p00; tmp = base_negacyclic_root; tm1 = cy_r; tm2 = cy_i; l = 0x4000;
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
			addr = a + j1 + pfetch_dist + poff[l>>1];	// poff[] = p0,4,8,...
			addr += (-(l&0x1)) & p02;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_fermat_carry_norm_pow2_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p01, add0);
			tm1 += 4; tmp += 2;
		}

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
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

/*...The radix-256 DIF pass is here:	*/

#ifdef USE_SSE2

	// Gather the needed data and do 8 twiddleless length-16 subtransforms, with p-offsets in br8 order: 04261537:

 #if USE_SCALAR_DFT_MACRO

	SSE2_RADIX256_DIF(
		// Input pointer: Base ptr of 16 local-mem:
		s1p00,
		// Intermediates base pointer:
		r00,
		// Pointers to base-roots data and first of 16 twiddle vectors:
		isrt2,two, twid0,
		// Outputs: Base address plus index offset lo/hi-half arrays:
		(a+j1),
		dft_offsets_lo,(uint32)0,dft_offsets_hi
	);

 #else

  #ifdef USE_ARM_V8_SIMD
	OFF1 = 0x200;
	OFF2 = 0x400;
	OFF3 = 0x600;
	OFF4 = 0x800;
  #elif defined(USE_AVX512)
	#define OFF1	4*0x200
	#define OFF2	4*0x400
	#define OFF3	4*0x600
	#define OFF4	4*0x800
  #elif defined(USE_AVX)
	#define OFF1	2*0x200
	#define OFF2	2*0x400
	#define OFF3	2*0x600
	#define OFF4	2*0x800
  #else
	#define OFF1	0x200
	#define OFF2	0x400
	#define OFF3	0x600
	#define OFF4	0x800
  #endif
// NOTE that unlike the RADIX_08_DIF_OOP() macro used for pass 1 of the radix-64 DFT, RADIX_16_DIF outputs are IN-ORDER rather than BR:
	tmp = r00;
	for(l = 0; l < 16; l++) {
		k1 = reverse(l,4)<<1;
		tm2 = s1p00 + k1;
		SSE2_RADIX16_DIF_0TWIDDLE_B(tm2,OFF1,OFF2,OFF3,OFF4, isrt2,two, tmp);
		tmp += 32;
	}

  #ifndef USE_ARM_V8_SIMD
	#undef OFF1
	#undef OFF2
	#undef OFF3
	#undef OFF4
  #endif

/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors: */

  #ifdef USE_ARM_V8_SIMD
	OFF1 = 0x200;
	OFF2 = 0x400;
	OFF3 = 0x600;
	OFF4 = 0x800;
  #elif defined(USE_AVX512)
	#define OFF1	4*0x200
	#define OFF2	4*0x400
	#define OFF3	4*0x600
	#define OFF4	4*0x800
  #elif defined(USE_AVX)
	#define OFF1	2*0x200
	#define OFF2	2*0x400
	#define OFF3	2*0x600
	#define OFF4	2*0x800
  #else
	#define OFF1	0x200
	#define OFF2	0x400
	#define OFF3	0x600
	#define OFF4	0x800
  #endif
	const int off_arr[16] = {0,p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f};
	off_ptr = &(off_arr[0]);	// off_ptr declared at top in DIT section
// Block 0: has all-unity twiddles
	add0 = &a[j1];
	SSE2_RADIX16_DIF_TWIDDLE_OOP(
		r00,OFF1,OFF4,
		add0,off_ptr,
		isrt2,twid0
	);

  #if defined(USE_AVX2) && !defined(USE_IMCI512)

/*** Only last 14 of the 15 with-twiddles DFTs allow use of FMA-based macros under Intel AVX2/FMA3: ***/
// Block 8: BR twiddles = {  I.{},  C^ 8,-~C^ 8,  C^ 4,*~C^ 4, *C^ 4,-~C^ 4,  C^ 2,*~C^ 2, *C^ 6,-~C^ 6,  C^ 6,*~C^ 6, *C^ 2,-~C^ 2}
	add0 = &a[j1] + p10;
	SSE2_RADIX16_DIF_TWIDDLE_OOP(
		r01,OFF1,OFF4,
		add0,off_ptr,
		isrt2,twid8
	);

	// Remaining 14 sets of macro calls done in loop:
	tm1 = r02;
	for(l = 2; l < 16; l++) {
		k1 = reverse(l,4)<<1;
		tmp = twid0 + (k1<<4)-k1;	// Twid-offsets are multiples of 30 vec_dbl
		add0 = &a[j1] + poff[l<<2];	// poff[4*i] = p10,p20,...,pf0
			add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
		add8 = add0+p08; add9 = add1+p08; adda = add2+p08; addb = add3+p08; addc = add4+p08; addd = add5+p08; adde = add6+p08; addf = add7+p08;
		SSE2_RADIX16_DIF_FMA_OOP(
			tm1,OFF1,OFF2,OFF3,OFF4, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, tmp
		);	tm1 += 2;
	}

  #else	// Non-FMA version:

	// Remaining 15 sets of macro calls done in loop:
	tm1 = r01;
	for(l = 1; l < 16; l++) {
		k1 = reverse(l,4)<<1;
		tmp = twid0 + (k1<<4)-k1;	// Twid-offsets are multiples of 30 vec_dbl
		add0 = &a[j1] + poff[l<<2];	// poff[4*i] = p10,p20,...,pf0
		SSE2_RADIX16_DIF_TWIDDLE_OOP(
			tm1,OFF1,OFF4,
			add0,off_ptr,
			isrt2,tmp
		);	tm1 += 2;
	}

  #endif	// FMA/AVX2 ?

  #ifndef USE_ARM_V8_SIMD
	#undef OFF1
	#undef OFF2
	#undef OFF3
	#undef OFF4
  #endif

 #endif	// USE_SCALAR_DFT_MACRO ?

#else	// USE_SSE2 = False:

 #if USE_SCALAR_DFT_MACRO

	RADIX_256_DIF(
		(a+j1),RE_IM_STRIDE,dft_offsets_lo,dft_offsets_hi,
		(a+j1),RE_IM_STRIDE,dft_offsets_lo,(uint32)0,dft_offsets_hi
	);

 #else

	// Gather the needed data and do 16 twiddleless length-16 subtransforms, with p-offsets in br8 order: 084c2a6e195d3b7f:
	// NOTE that RADIX_16_DIF outputs are IN-ORDER rather than BR:
	tptr = t;
	for(ntmp = 0; ntmp < 16; ntmp++) {
		jt = j1 + po_br[ntmp]; jp = j2 + po_br[ntmp];	// po_br[] = p[084c2a6e195d3b7f]
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70],a[jt+p80],a[jp+p80],a[jt+p90],a[jp+p90],a[jt+pa0],a[jp+pa0],a[jt+pb0],a[jp+pb0],a[jt+pc0],a[jp+pc0],a[jt+pd0],a[jp+pd0],a[jt+pe0],a[jp+pe0],a[jt+pf0],a[jp+pf0],
			tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
			c16,s16	// These = (c16,s16) typically def'd for use in the radix-16 DFT
		);	tptr += 16;
	}

/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors: */
/*
Note: As used here, Inputs to RADIX_16_DIF_TWIDDLE_OOP maxro ordered as
	__A0-f = t0@,8@,4@,c@,2@,a@,6@,e@,1@,9@,5@,d@,3@,b@,7@,f@
So internal DFT-4 blocks use data in-order!
*/
	// Block 0: has all-unity twiddles
	tptr = t;
	jt = j1;	jp = j2;
	// Twiddleless DIF bit-reverses its outputs, so a_p* terms appear in BR-order [swap index pairs 1/8,2/4,3/c,5/a,7/e.b/d]:
	RADIX_16_DIF(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c16,s16
	);	tptr++;

/*** For the last 14 of the 15 with-twiddles DFTs to test with-FMA option in prep. for porting to Intel AVX2/FMA3: ***/
  #ifdef USE_FMA	// Must define - or not - @compile time

/* Can't use FMA-based DFT for first set of non-unity twiddles - the leading [0+i] twiddle leads to div-by-0 in our tangent-based FMA-DFT scheme:
// Block 8: BR twiddles = {  I.{},  C^ 8,-~C^ 8,  C^ 4,*~C^ 4, *C^ 4,-~C^ 4,  C^ 2,*~C^ 2, *C^ 6,-~C^ 6,  C^ 6,*~C^ 6, *C^ 2,-~C^ 2}
	jt = j1 + p10;	jp = j2 + p10;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c1_1,s1_1,c1_2,s1_2,c1_3,s1_3,c1_4,s1_4,c1_5,s1_5,c1_6,s1_6,c1_7,s1_7,c1_8,s1_8,c1_9,s1_9,c1_10,s1_10,c1_11,s1_11,c1_12,s1_12,c1_13,s1_13,c1_14,s1_14,c1_15,s1_15,c1_1_c,tan,c1_1i2,c1_2i2
	);	tptr++;
*/
// Block 8: BR twiddles = {  I.{},  C^ 8,-~C^ 8,  C^ 4,*~C^ 4, *C^ 4,-~C^ 4,  C^ 2,*~C^ 2, *C^ 6,-~C^ 6,  C^ 6,*~C^ 6, *C^ 2,-~C^ 2}
	jt = j1 + p10;	jp = j2 + p10;
	RADIX_16_DIF_TWIDDLE_OOP(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		0,1, ISRT2,ISRT2, -ISRT2,ISRT2, c16,s16, -s16,c16, s16,c16, -c16,s16, c32_1,s32_1, -s32_1,c32_1, s32_3,c32_3, -c32_3,s32_3, c32_3,s32_3, -s32_3,c32_3, s32_1,c32_1, -c32_1,s32_1,
		c16,s16
	);	tptr++;
// Block 4: BR twiddles = {  C^ 8,  C^ 4, *C^ 4,  C^ 2, *C^ 6,  C^ 6, *C^ 2,  C^ 1, *C^ 7,  C^ 5, *C^ 3,  C^ 3, *C^ 5,  C^ 7, *C^ 1}
	jt = j1 + p20;	jp = j2 + p20;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c2_1,s2_1,c2_2,s2_2,c2_3,s2_3,c2_4,s2_4,c2_5,s2_5,c2_6,s2_6,c2_7,s2_7,c2_8,s2_8,c2_9,s2_9,c2_10,s2_10,c2_11,s2_11,c2_12,s2_12,c2_13,s2_13,c2_14,s2_14,c2_15,s2_15,c2_1_c,tan,c2_1i2,c2_2i2
	);	tptr++;
// Block c: BR twiddles = {-~C^ 8, *C^ 4, -C^ 4,  C^ 6,-~C^ 2,*~C^ 2,-*C^ 6,  C^ 3,-~C^ 5, *C^ 1, -C^ 7, *C^ 7, -C^ 1,*~C^ 5,-*C^ 3}
	jt = j1 + p30;	jp = j2 + p30;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c3_1,s3_1,c3_2,s3_2,c3_3,s3_3,c3_4,s3_4,c3_5,s3_5,c3_6,s3_6,c3_7,s3_7,c3_8,s3_8,c3_9,s3_9,c3_10,s3_10,c3_11,s3_11,c3_12,s3_12,c3_13,s3_13,c3_14,s3_14,c3_15,s3_15,c3_1_c,tan,c3_1i2,c3_2i2
	);	tptr++;
// Block 2: BR twiddles = {  C^ 4,  C^ 2,  C^ 6,  C^ 1,  C^ 5,  C^ 3,  C^ 7,  D^ 1,  D^ 9,  D^ 5,  D^ d,  D^ 3,  D^ b,  D^ 7,  D^ f}
	jt = j1 + p40;	jp = j2 + p40;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c4_1,s4_1,c4_2,s4_2,c4_3,s4_3,c4_4,s4_4,c4_5,s4_5,c4_6,s4_6,c4_7,s4_7,c4_8,s4_8,c4_9,s4_9,c4_10,s4_10,c4_11,s4_11,c4_12,s4_12,c4_13,s4_13,c4_14,s4_14,c4_15,s4_15,c4_1_c,tan,c4_1i2,c4_2i2
	);	tptr++;
// Block a: BR twiddles = {*~C^ 4, *C^ 6,-~C^ 2,  C^ 5,-~C^ 7, *C^ 1, -C^ 3,  D^ 5,*~D^ d, *D^ 7, -D^ 1,  D^ f,-~D^ 9,*~D^ 3, -D^ b}
	jt = j1 + p50;	jp = j2 + p50;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c5_1,s5_1,c5_2,s5_2,c5_3,s5_3,c5_4,s5_4,c5_5,s5_5,c5_6,s5_6,c5_7,s5_7,c5_8,s5_8,c5_9,s5_9,c5_10,s5_10,c5_11,s5_11,c5_12,s5_12,c5_13,s5_13,c5_14,s5_14,c5_15,s5_15,c5_1_c,tan,c5_1i2,c5_2i2
	);	tptr++;
// Block 6: BR twiddles = { *C^ 4,  C^ 6,*~C^ 2,  C^ 3, *C^ 1, *C^ 7,*~C^ 5,  D^ 3, *D^ 5,  D^ f,*~D^ 7,  D^ 9,*~D^ 1, *D^ b,*~D^ d}
	jt = j1 + p60;	jp = j2 + p60;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c6_1,s6_1,c6_2,s6_2,c6_3,s6_3,c6_4,s6_4,c6_5,s6_5,c6_6,s6_6,c6_7,s6_7,c6_8,s6_8,c6_9,s6_9,c6_10,s6_10,c6_11,s6_11,c6_12,s6_12,c6_13,s6_13,c6_14,s6_14,c6_15,s6_15,c6_1_c,tan,c6_1i2,c6_2i2
	);	tptr++;
// Block e: BR twiddles = {-~C^ 4, *C^ 2,-*C^ 6,  C^ 7, -C^ 3,*~C^ 5,~*C^ 1,  D^ 7,-~D^ 1,*~D^ 3,-*D^ 5, *D^ b, -D^ d,-~D^ f,~*D^ 9}
	jt = j1 + p70;	jp = j2 + p70;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c7_1,s7_1,c7_2,s7_2,c7_3,s7_3,c7_4,s7_4,c7_5,s7_5,c7_6,s7_6,c7_7,s7_7,c7_8,s7_8,c7_9,s7_9,c7_10,s7_10,c7_11,s7_11,c7_12,s7_12,c7_13,s7_13,c7_14,s7_14,c7_15,s7_15,c7_1_c,tan,c7_1i2,c7_2i2
	);	tptr++;
// Block 1: BR twiddles = {  C^ 2,  C^ 1,  C^ 3,  D^ 1,  D^ 5,  D^ 3,  D^ 7,  E^ 1,  E^ 9,  E^ 5,  E^ d,  E^ 3,  E^ b,  E^ 7,  E^ f}
	jt = j1 + p80;	jp = j2 + p80;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c8_1,s8_1,c8_2,s8_2,c8_3,s8_3,c8_4,s8_4,c8_5,s8_5,c8_6,s8_6,c8_7,s8_7,c8_8,s8_8,c8_9,s8_9,c8_10,s8_10,c8_11,s8_11,c8_12,s8_12,c8_13,s8_13,c8_14,s8_14,c8_15,s8_15,c8_1_c,tan,c8_1i2,c8_2i2
	);	tptr++;
// Block 9: BR twiddles = {*~C^ 2, *C^ 7,-~C^ 5,  D^ 9,*~D^ d, *D^ 5,-~D^ 1,  E^ 9,*~E^ h, *E^ j,-~E^ b,  E^ r,-~E^ t, *E^ 1, -E^ 7}
	jt = j1 + p90;	jp = j2 + p90;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		c9_1,s9_1,c9_2,s9_2,c9_3,s9_3,c9_4,s9_4,c9_5,s9_5,c9_6,s9_6,c9_7,s9_7,c9_8,s9_8,c9_9,s9_9,c9_10,s9_10,c9_11,s9_11,c9_12,s9_12,c9_13,s9_13,c9_14,s9_14,c9_15,s9_15,c9_1_c,tan,c9_1i2,c9_2i2
	);	tptr++;
// Block 5: BR twiddles = { *C^ 6,  C^ 5, *C^ 1,  D^ 5, *D^ 7,  D^ f,*~D^ 3,  E^ 5, *E^ j,  E^ p,*~E^ 1,  E^ f, *E^ 9, *E^ t,*~E^ b}
	jt = j1 + pa0;	jp = j2 + pa0;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		ca_1,sa_1,ca_2,sa_2,ca_3,sa_3,ca_4,sa_4,ca_5,sa_5,ca_6,sa_6,ca_7,sa_7,ca_8,sa_8,ca_9,sa_9,ca_10,sa_10,ca_11,sa_11,ca_12,sa_12,ca_13,sa_13,ca_14,sa_14,ca_15,sa_15,ca_1_c,tan,ca_1i2,ca_2i2
	);	tptr++;
// Block d: BR twiddles = {-~C^ 6, *C^ 3, -C^ 7,  D^ d, -D^ 1,*~D^ 7,-*D^ 5,  E^ d,-~E^ b,*~E^ 1,-*E^ n, *E^ p, -E^ f,*~E^ r,~*E^ 3}
	jt = j1 + pb0;	jp = j2 + pb0;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		cb_1,sb_1,cb_2,sb_2,cb_3,sb_3,cb_4,sb_4,cb_5,sb_5,cb_6,sb_6,cb_7,sb_7,cb_8,sb_8,cb_9,sb_9,cb_10,sb_10,cb_11,sb_11,cb_12,sb_12,cb_13,sb_13,cb_14,sb_14,cb_15,sb_15,cb_1_c,tan,cb_1i2,cb_2i2
	);	tptr++;
// Block 3: BR twiddles = {  C^ 6,  C^ 3, *C^ 7,  D^ 3,  D^ f,  D^ 9, *D^ b,  E^ 3,  E^ r,  E^ f, *E^ p,  E^ 9, *E^ v,  E^ l, *E^ j}
	jt = j1 + pc0;	jp = j2 + pc0;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		cc_1,sc_1,cc_2,sc_2,cc_3,sc_3,cc_4,sc_4,cc_5,sc_5,cc_6,sc_6,cc_7,sc_7,cc_8,sc_8,cc_9,sc_9,cc_10,sc_10,cc_11,sc_11,cc_12,sc_12,cc_13,sc_13,cc_14,sc_14,cc_15,sc_15,cc_1_c,tan,cc_1i2,cc_2i2
	);	tptr++;
// Block b: BR twiddles = {*~C^ 6, *C^ 5, -C^ 1,  D^ b,-~D^ 9,*~D^ 1, -D^ d,  E^ b,-~E^ t, *E^ 9, -E^ f, *E^ v,-~E^ 7,*~E^ d,-*E^ r}
	jt = j1 + pd0;	jp = j2 + pd0;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		cd_1,sd_1,cd_2,sd_2,cd_3,sd_3,cd_4,sd_4,cd_5,sd_5,cd_6,sd_6,cd_7,sd_7,cd_8,sd_8,cd_9,sd_9,cd_10,sd_10,cd_11,sd_11,cd_12,sd_12,cd_13,sd_13,cd_14,sd_14,cd_15,sd_15,cd_1_c,tan,cd_1i2,cd_2i2
	);	tptr++;
// Block 7: BR twiddles = { *C^ 2,  C^ 7,*~C^ 5,  D^ 7,*~D^ 3, *D^ b,-~D^ f,  E^ 7, *E^ 1, *E^ t,*~E^ r,  E^ l,*~E^ d, *E^ f,-~E^ n}
	jt = j1 + pe0;	jp = j2 + pe0;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		ce_1,se_1,ce_2,se_2,ce_3,se_3,ce_4,se_4,ce_5,se_5,ce_6,se_6,ce_7,se_7,ce_8,se_8,ce_9,se_9,ce_10,se_10,ce_11,se_11,ce_12,se_12,ce_13,se_13,ce_14,se_14,ce_15,se_15,ce_1_c,tan,ce_1i2,ce_2i2
	);	tptr++;
// Block f: BR twiddles = {-~C^ 2, *C^ 1,-*C^ 3,  D^ f, -D^ b,*~D^ d,~*D^ 9,  E^ f, -E^ 7,*~E^ b,~*E^ 3, *E^ j,-*E^ r,-~E^ n, ~E^ v}
	jt = j1 + pf0;	jp = j2 + pf0;
	RADIX_16_DIF_FMA(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
		cf_1,sf_1,cf_2,sf_2,cf_3,sf_3,cf_4,sf_4,cf_5,sf_5,cf_6,sf_6,cf_7,sf_7,cf_8,sf_8,cf_9,sf_9,cf_10,sf_10,cf_11,sf_11,cf_12,sf_12,cf_13,sf_13,cf_14,sf_14,cf_15,sf_15,cf_1_c,tan,cf_1i2,cf_2i2
	);	tptr++;

  #else		// USE_FMA = False

	// Remaining 15 sets of macro calls done in loop:
	for(ntmp = 1; ntmp < 16; ntmp++) {
		jt = j1 + poff[ntmp<<2]; jp = j2 + poff[ntmp<<2];	// poff[4*i] = p10,p20,...,pf0
		const double *addr = DFT256_TWIDDLES[ntmp], *addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_16_DIF_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
			c16,s16
		);	tptr++;
	}

  #endif	// USE_FMA ?

 #endif	// USE_SCALAR_DFT_MACRO ?

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
