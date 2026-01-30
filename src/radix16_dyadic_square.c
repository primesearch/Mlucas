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

#include "Mlucas.h"

#define EPS 1e-10

#ifndef PFETCH_DIST
  #ifdef USE_AVX
	#define PFETCH_DIST	4096	// This seems to work best on my Haswell
  #else
	#define PFETCH_DIST	1024
  #endif
#endif

#ifdef USE_SSE2

	#include "sse2_macro_gcc64.h"
	#include "radix16_utils_asm.h"

	  // Using != 0 *value* of flag to control preprocessor inlining allows us to override default setting at compile time:
	  #ifndef FULLY_FUSED
		#define FULLY_FUSED		0	// For 64-bit GCC-mode builds, this invokes a fused DIF/Square/DIT macro.
	  #endif
	  #if FULLY_FUSED
		#warning Using FULLY_FUSED version of radix-16 DIF/dyadic-square/DIT algorithm.
		#error FULLY_FUSED path does not support Gerbicz-check!
		#include "radix16_dyadic_square_gcc64.h"
	  #else
		#include "radix16_wrapper_square_gcc64.h"
	  #endif

#endif

#ifdef MULTITHREAD
	#ifndef USE_PTHREAD
		#error Multithreading currently only supported using Pthreads!
	#endif
#endif

/***************/

/*
!   NOTE: In the following commentary, N refers to the COMPLEX vector length (N2 in the code),
!   which is half the real vector length.
!
!...Acronyms:   DIT = Decimation In Time
!               DIF = Decimation In Frequency
!               FFT = Fast Fourier Transform, i.e. a discrete FT over the complex numbers
!
!...Final pass of the DIF forward transform, pointwise squaring and initial pass of
!   the inverse transform, all performed in one fell swoop (or better, one swell loop :)

The scratch array (2nd input argument) is only needed for data table initializations, i.e. if init_sse2 = TRUE.
*/
/* v20: To support Gerbicz check (and later, p-1) need 2-input FFT-mul - added fwd_fft_only flag, which takes 3 distinct possible values:
	  0 - Do normal FFT-autosquare(a) for (ilo-ihi) successive iterations;
	  1 - Do forward FFT(a) only and store result in a[], skipping dyadic-square and inv-FFT steps. In this case ilo = ihi;
	> 1 - The fwd_fft_only arg contains a pointer to an array b[] in previously-forward-FFTed form, i.e. FFT(b). In this
		  case we do the final-radix pass of FFT(a), then dyadic-multiply FFT(a) * FFT(b) (storing the result in a[]) in place
		  of the usual pair_square step, then do the initial-radix pass of the iFFT of the result.
*/
void radix16_dyadic_square(
	double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[],
	int ii, int nradices_prim, int radix_prim[], int incr, int init_sse2, int thr_id, uint64 fwd_fft_only, double c_arr[]
)
{
	const char func[] = "radix16_dyadic_square";
#ifdef USE_AVX
	const int pfetch_dist = PFETCH_DIST;
#endif
	const int stride = (int)RE_IM_STRIDE << 5, stridh = (stride>>1);	// main-array loop stride = 32*RE_IM_STRIDE
	static int max_threads = 0;
	static int nsave = 0;
	static int rad0save = 0, ndivrad0 = 0;
/*	static int *index = 0x0;	OBSOLETE: full-N2/16-length Bit-reversal index array. */
	static int *index0 = 0x0, *index1 = 0x0, *index_ptmp0 = 0x0, *index_ptmp1 = 0x0;
		   int index0_idx=-1, index1_idx=-1;
	static int index0_mod=-1, index1_mod=-1;
	int nradices_prim_radix0;
	int i,j,j1,iroot/* ,k */;
#ifndef USE_AVX512
	int l,k1,k2;
#endif
#ifdef USE_SSE2
	int nbytes;
#endif
#ifndef USE_SSE2
	int j2;
#endif
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)] */
#ifndef USE_SSE2
	double re0,im0,re1,im1;
#endif
	double rt,it;

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error X86 SIMD code not supported for this compiler!
  #endif
	static uint32 *sm_arr = 0x0,*sm_ptr;	// Base-ptr to arrays of k1,k2-index-vectors used in SIMD roots-computation.
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0,*add1;	/* Addresses into array sections */
  #ifndef USE_ARM_V8_SIMD
	double *bdd0;		/* Addresses into array sections */
  #endif
  #ifndef USE_ARM_V8_SIMD
	const double *cdd0;
  #endif
  #ifdef USE_AVX
	double *add2,*add3;
  #endif
  #ifdef USE_AVX512
	double *add4,*add5,*add6,*add7;
  #endif
	vec_dbl *c_tmp,*s_tmp;
  #ifdef USE_IMCI512
	vec_dbl *tmp,*tm1;
  #endif
  #ifdef MULTITHREAD
	// Base addresses for discrete per-thread local stores ... 'r' for double-float data, 'i' for int:
	static vec_dbl *__r0;
	static uint32  *__i0;
	// In || mode, only above base-pointers (shared by all threads) are static:
	uint32  *k1_arr,*k2_arr;
	vec_dbl *cc0,*ss0,*isrt2, *two,*one, *r1,*r9,*r17,*r25,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15;
  #else
	static uint32  *k1_arr,*k2_arr;
	static vec_dbl *cc0,*ss0,*isrt2, *two,*one, *r1,*r9,*r17,*r25,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15;
  #endif

#else

  #if PFETCH
	double *add0;
  #endif
	double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
	,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1pAr,aj1pBr,aj1pCr,aj1pDr,aj1pEr,aj1pFr
	,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1pAi,aj1pBi,aj1pCi,aj1pDi,aj1pEi,aj1pFi
	,bj1p0r,bj1p1r,bj1p2r,bj1p3r,bj1p4r,bj1p5r,bj1p6r,bj1p7r,bj1p8r,bj1p9r,bj1pAr,bj1pBr,bj1pCr,bj1pDr,bj1pEr,bj1pFr
	,bj1p0i,bj1p1i,bj1p2i,bj1p3i,bj1p4i,bj1p5i,bj1p6i,bj1p7i,bj1p8i,bj1p9i,bj1pAi,bj1pBi,bj1pCi,bj1pDi,bj1pEi,bj1pFi;

#endif

	// For case of 2-input modmul, pointer to 2nd (presumed already fwd-FFTed) input array is supplied via fwd_fft_only argument:
	/* v20: Bits 2:3 of fwd_fft = 3 means "dyadic-multiply of 2 inputs a and b, with both already forward-FFTed:
		(double *)a = FFT(a), (double *)(fwd_fft_only - mode_flag) = FFT(b).
	In this case we require bits 0:1 == 0, and fwd_fft = & ~0xC yields pointer to FFT(b), and we skip over fwd-FFT directly to
	the dyadic-multiply FFT(a) * FFT(b) step, then iFFT the product, storing the result in a[].
	*/
#if !defined(USE_SSE2) || (!FULLY_FUSED && !defined(USE_ARM_V8_SIMD))
	double *b = 0x0;
#endif
	if(fwd_fft_only >> 2) {		// Do the following 3 steps in both cases - if bits 2:3 == 0 the ANDs are no-ops...
		// The submul-auxiliary array c_arr[], if present, in already in proper double[] form, but b[] needs low-bits-cleared and casting
	  #if !defined(USE_SSE2) || (!FULLY_FUSED && !defined(USE_ARM_V8_SIMD))
		b = (double *)(fwd_fft_only & ~0xCull);
	  #endif
		// BUT, if bits 2:3 == 0, must avoid zeroing fwd_fft_only since "do 2-input dyadic-mul following fwd-FFT" relies on that != 0:
		if(fwd_fft_only & 0xC) {
			ASSERT((fwd_fft_only & 0xF) == 0xC,"Illegal value for bits 2:3 of fwd_fft_only!");	// Otherwise bits 2:3 should've been zeroed prior to entry
			fwd_fft_only = 3ull;
		}
	}

/*...If a new runlength or first-pass radix, it is assumed this function has been first-called with init_sse2 = true to
     initialize static data and lcoal-storage prior to actual use in computing a transform-based result.
*/

	/* In order to eschew complex thread-block-and-sync logic related to the local-store-init step in multithread mode,
	switch to a special-init-mode-call paradigm, in which the function is inited once (for as many threads as needed)
	prior to being executed:
	*/
/**************************************************************************************************************************************/
/*** To-Do: Need to add code to allow for re-init when any of the FFT-related params or #threads changes during course of execution ***/
/**************************************************************************************************************************************/
	if((rad0save != radix0) || (nsave != n))
	{
		ASSERT(thr_id == -1, "Init-mode call must be outside of any multithreading!");
		nsave = n;
		ASSERT(N2 == n/2, "N2 bad!");
		rad0save = radix0;
		ndivrad0 = n/radix0;
		for(j = 0; j < ndivrad0; j += stride)
		{
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
			if( (j1+stridh) != (j+stridh) + ( ((j+stridh) >> DAT_BITS) << PAD_BITS ) ) {
				printf("j, j1, stride/2 = %d,%d,%d, jpad = %d\n",j,j1, stridh, (j+stridh) + (((j+stridh) >> DAT_BITS) << PAD_BITS) );
				ASSERT(0 , "add1 calculation violates padded index rules!");
			}
		}
		if(index_ptmp0) {
			free((void *)index_ptmp0);	index_ptmp0=0x0;	index0=0x0;
			free((void *)index_ptmp1);	index_ptmp1=0x0;	index1=0x0;
		}
		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Allocate and initialize an index array containing N/16 indices...

		index_ptmp = ALLOC_INT(N2/16);
		index = ALIGN_INT(index_ptmp);
		if(!index){ sprintf(cbuf,"ERROR: unable to allocate array ITMP in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		for(i=0; i < N2/16; i++)
		{
			index[i]=i;
		}
		*/
		index0_mod =        radix0;
		index1_mod = (n>>5)/radix0;	/* complex length requires an additional divide by 2 */

		index_ptmp0 = ALLOC_INT(index_ptmp0, index0_mod);
		if(!index_ptmp0){ sprintf(cbuf,"ERROR: unable to allocate array INDEX_PTMP0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		index0 = ALIGN_INT(index_ptmp0);

		index_ptmp1 = ALLOC_INT(index_ptmp1, index1_mod);
		if(!index_ptmp1){ sprintf(cbuf,"ERROR: unable to allocate array INDEX_PTMP1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		index1 = ALIGN_INT(index_ptmp1);

		for(i=0; i < index0_mod; i++){index0[i]=       i;}
		for(i=0; i < index1_mod; i++){index1[i]=radix0*i;}

		/*...then bit-reverse INDEX with respect to N/16. Radices are sent to BR routine
		in opposite order from which they are processed in forward FFT.

		bit_reverse_int(index, N2/16, nradices_prim-4, &radix_prim[nradices_prim-5], -1,(int *)arr_scratch);
		*/

		i = 1;
		for(nradices_prim_radix0 = 1; nradices_prim_radix0 < nradices_prim; nradices_prim_radix0++)
		{
			i *= radix_prim[nradices_prim_radix0-1];
			if(i == radix0)
				break;
		}
		if(nradices_prim_radix0 >= nradices_prim) { sprintf(cbuf,"ERROR: nradices_prim_radix0 must be < nradices_prim in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

		bit_reverse_int(index0, index0_mod,                 nradices_prim_radix0, &radix_prim[nradices_prim_radix0-1], -1,(int *)arr_scratch);
		bit_reverse_int(index1, index1_mod, nradices_prim-4-nradices_prim_radix0, &radix_prim[nradices_prim       -5], -1,(int *)arr_scratch);

		if(!init_sse2)	// SSE2 stuff pvsly inited, but new leading radix (e.g. self-test mode)
			return;
	}

	/* Here this variable is somewhat misnamed because it is used to init both non-SIMD and SIMD-specific data */
	if(init_sse2)	// Just check nonzero here, to allow the *value* of init_sse2 to store #threads
	{
		if(init_sse2 <= max_threads)	// current alloc sufficient
			return;

		ASSERT(thr_id == -1, "Init-mode call must be outside of any multithreading!");
		max_threads = init_sse2;
	#ifndef COMPILER_TYPE_GCC
		ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif

	#ifdef USE_SSE2
	//	fprintf(stderr, "%s: pfetch_dist = %d\n",func,pfetch_dist);
		if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)sm_arr);	sm_arr=0x0;
			free((void *)sc_arr);	sc_arr=0x0;
		}
		// Index vectors used in SIMD roots-computation.
		sm_arr = ALLOC_UINT(sm_arr, max_threads*10*RE_IM_STRIDE + 16);	if(!sm_arr){ sprintf(cbuf, "ERROR: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sm_ptr = ALIGN_UINT(sm_arr);
		ASSERT(((uintptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");
		// Twiddles-array:
		sc_arr = ALLOC_VEC_DBL(sc_arr, 72*max_threads + 100);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((uintptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	last 30 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
	  #ifdef MULTITHREAD
		__i0  = sm_ptr;
		__r0  = sc_ptr;
		isrt2 = sc_ptr + 0x20;
		cc0   = sc_ptr + 0x21;
		ss0   = sc_ptr + 0x22;
		one   = sc_ptr + 0x43;	// Mar 2017: Somehow order of one/two-ptrs switched here vs. radix_16_wrapper_square ... swap to match.
		two   = sc_ptr + 0x44;
		for(i = 0; i < max_threads; ++i) {
			/* These remain fixed within each per-thread local store: */
			VEC_DBL_INIT(isrt2, ISRT2);
			VEC_DBL_INIT(one  , 1.0);
			VEC_DBL_INIT(two  , 2.0);
			VEC_DBL_INIT(cc0  , c);
			VEC_DBL_INIT(ss0  , s);
			isrt2 += 72;	/* Move on to next thread's local store */
			cc0   += 72;
			ss0   += 72;
			one   += 72;
			two   += 72;
		}
	  #elif defined(COMPILER_TYPE_GCC)
		k1_arr = sm_ptr;
		k2_arr = sm_ptr + 5*RE_IM_STRIDE;
		r1    = sc_ptr;
		r9    = sc_ptr + 0x08;
		r17   = sc_ptr + 0x10;
		r25   = sc_ptr + 0x18;
		isrt2 = sc_ptr + 0x20;
		cc0   = sc_ptr + 0x21;
		ss0   = sc_ptr + 0x22;
		c8    = sc_ptr + 0x25;	/* Need some added root-addresses here for expedited roots computation via SSE2_CMUL_EXPO macros */
		c4    = sc_ptr + 0x27;
		c12   = sc_ptr + 0x29;
		c2    = sc_ptr + 0x2b;
		c10   = sc_ptr + 0x2d;
		c6    = sc_ptr + 0x2f;
		c14   = sc_ptr + 0x31;
		c1    = sc_ptr + 0x33;
		c9    = sc_ptr + 0x35;
		c5    = sc_ptr + 0x37;
		c13   = sc_ptr + 0x39;
		c3    = sc_ptr + 0x3b;
		c11   = sc_ptr + 0x3d;
		c7    = sc_ptr + 0x3f;
		c15   = sc_ptr + 0x41;
		one   = sc_ptr + 0x43;
		two   = sc_ptr + 0x44;
		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);
		VEC_DBL_INIT(one  , 1.0);
		VEC_DBL_INIT(two  , 2.0);
		VEC_DBL_INIT(cc0  , c);
		VEC_DBL_INIT(ss0  , s);
	  #endif
	#endif	// USE_SSE2
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
#ifdef MULTITHREAD
	ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
  #ifdef USE_SSE2
	k1_arr =   __i0 + thr_id*10*RE_IM_STRIDE;
	k2_arr = k1_arr + 5*RE_IM_STRIDE;
	r1    = __r0 + thr_id*72;
	r9    = r1 + 0x08;
	r17   = r1 + 0x10;
	r25   = r1 + 0x18;
	isrt2 = r1 + 0x20;
	cc0   = r1 + 0x21;
	c8    = r1 + 0x25;	/* Need some added root-addresses here for expedited roots computation via SSE2_CMUL_EXPO macros */
	c4    = r1 + 0x27;
	c12   = r1 + 0x29;
	c2    = r1 + 0x2b;
	c10   = r1 + 0x2d;
	c6    = r1 + 0x2f;
	c14   = r1 + 0x31;
	c1    = r1 + 0x33;
	c9    = r1 + 0x35;
	c5    = r1 + 0x37;
	c13   = r1 + 0x39;
	c3    = r1 + 0x3b;
	c11   = r1 + 0x3d;
	c7    = r1 + 0x3f;
	c15   = r1 + 0x41;
  #endif
#endif

	/*...If a new runlength, should not get to this point: */
	ASSERT(n == nsave,"n != nsave");
	ASSERT(incr == 32,"incr != 32");
	ASSERT(ndivrad0 == n/radix0,"bad value for ndivrad0!");
	/*
	k = ii*(ndivrad0 >> 5);
	*/
	index0_idx = ii;
	index1_idx = 0;
//printf("\nJ = [%d-%d]:\n",0,ndivrad0);
	for(j = 0; j < ndivrad0; j += stride)
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		j2 = j1+RE_IM_STRIDE;
	#endif

	#ifndef USE_SSE2	// Scalar-double mode:

		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c1 =rt;	s1 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += (iroot << 1);	/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c2 =rt;	s2 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += (iroot << 2);	/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c4 =rt;	s4 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += (iroot << 2) + iroot;	/* 13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c8 =rt;	s8 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c13 =rt;	s13 =it;

		/* c3,5 */
		t1=c1 *c4 ;	t2=c1 *s4 ;	rt=s1 *c4 ;	it=s1 *s4;
		c3 =t1 +it;	s3 =t2 -rt;	c5 =t1 -it;	s5 =t2 +rt;

		/* c6,7,9,10 */
		t1=c1 *c8 ;	t2=c1 *s8 ;	rt=s1 *c8 ;	it=s1 *s8;
		c7 =t1 +it;	s7 =t2 -rt;	c9 =t1 -it;	s9 =t2 +rt;

		t1=c2 *c8 ;	t2=c2 *s8 ;	rt=s2 *c8 ;	it=s2 *s8;
		c6 =t1 +it;	s6 =t2 -rt;	c10=t1 -it;	s10=t2 +rt;

		/* c11,12,14,15 */
		t1=c1 *c13;	t2=c1 *s13;	rt=s1 *c13;	it=s1 *s13;
		c12=t1 +it;	s12=t2 -rt;	c14=t1 -it;	s14=t2 +rt;

		t1=c2 *c13;	t2=c2 *s13;	rt=s2 *c13;	it=s2 *s13;
		c11=t1 +it;	s11=t2 -rt;	c15=t1 -it;	s15=t2 +rt;

	#elif 0//def USE_AVX2	// Special AVX2 code: do all 4 elts of sincos quartet using gather-load and SIMD arithmetic
	  #error To-do!
		struct uint32x4 l_4,iroot_4;

	/*=======================
	AVX2: Try doing the roots using VGATHERDPD/VGATHERQPD gather-load and inline ASM for the CMUL:

	VGATHERDPD/VGATHERQPD Description:

	The instruction conditionally loads up to 2 or 4 double-precision floating-point values from memory addresses
	specified by the memory operand (the second operand) and using qword indices. The memory operand uses the VSIB
	form of the SIB byte to specify a general purpose register operand as the common base, a vector register for
	an array of indices relative to the base and a constant scale factor.

	The mask operand (the third operand) specifies the conditional load operation from each memory address and the
	corresponding update of each data element of the destination operand (the first operand). Conditionality is
	specified by the most significant bit of each data element of the mask register. If an element’s mask bit is
	not set, the corresponding element of the destination register is left unchanged. The width of data element in
	the destination register and mask register are identical. The entire mask register will be set to zero by this
	instruction unless the instruction causes an exception.

	Using dword indices in the lower half of the mask register, the instruction conditionally loads up to 2 or 4
	double-precision floating-point values from the VSIB addressing memory operand, and updates the destination register.

	VGATHERDPD/VGATHERQPD Operation:

	DEST <-- SRC1;
	BASE_ADDR: base register encoded in VSIB addressing;
	VINDEX: the vector index register encoded by VSIB addressing;
	SCALE: scale factor encoded by SIB:[7:6];
	DISP: optional 1, 2, 4 byte displacement;
	MASK <-- SRC3;

	EWM:
	o We have 2 32-bit index quartets (i.e. each quartet = 128 bits) based on quads of k1 and k2 values
	o For the k1,2-quartets the load base addresses = rt0,1, resp.

	VGATHERDPD (VEX.256 version):

		FOR j<--0 to 3 i<--j * 64;
			IF MASK[63+i] THEN MASK[i +63:i]<--0xFFFFFFFF_FFFFFFFF; // extend from most significant bit
			ELSE MASK[i +63:i]<--0;
		FI; ENDFOR
		FOR j<--0 to 3 k<--j * 32;
			i<--j * 64; DATA_ADDR<--BASE_ADDR + (SignExtend(VINDEX1[k+31:k])*SCALE + DISP; IF MASK[63+i] THEN
			DEST[i +63:i]<--FETCH_64BITS(DATA_ADDR); // a fault exits the loop
		FI;
		MASK[i +63:i]<--0; ENDFOR

	And basics of VSIB:

	4.2	VECTOR SIB (VSIB) MEMORY ADDRESSING
	In AVX2, an SIB byte that follows the ModR/M byte can support VSIB memory addressing to an array of linear
	addresses. VSIB addressing is only supported in a subset of AVX2 instructions. VSIB memory addressing requires
	32-bit or 64-bit effec- tive address. In 32-bit mode, VSIB addressing is not supported when address size attribute
	is overridden to 16 bits. In 16-bit protected mode, VSIB memory addressing is permitted if address size attribute
	is overridden to 32 bits. Additionally, VSIB memory addressing is supported only with VEX prefix.

	In VSIB memory addressing, the SIB byte consists of:
	•	The scale field (bit 7:6) specifies the scale factor.
	•	The index field (bits 5:3) specifies the register number of the vector index register,
		each element in the vector register specifies an index.
	•	The base field (bits 2:0) specifies the register number of the base register.

	=======================
	*/
		//=== Need 4-way vector-SIMD version of this ===
		iroot_4->d0 = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		iroot_4->d1 = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		iroot_4->d2 = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		iroot_4->d3 = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
	//=== start asm here ===
	__asm__ volatile (\
		"movq	%[__iroot_4],%%rax		\n\t"\
		"vmovaps	 (%%rax),%%xmm0		\n\t"/* iroot quartet is 4x32-bit int */\
		"vmovaps	 %%xmm0,%%xmm1		\n\t"/* l = iroot */\
		:					// outputs: none
		: [__iroot_4] "m" (iroot_4)	// All inputs from memory addresses here
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	// Clobbered registers
	);
//=====================
//*** need full inline-asm for these 5 macro calls ***
		// Now generate the remaining roots from the smaller "seed set":
		SSE2_CMUL_EXPO(c1,c4 ,c3 ,c5 );
		SSE2_CMUL_EXPO(c1,c8 ,c7 ,c9 );
		SSE2_CMUL_EXPO(c2,c8 ,c6 ,c10);
		SSE2_CMUL_EXPO(c1,c13,c12,c14);
		SSE2_CMUL_EXPO(c2,c13,c11,c15);

	#else	// SIMD:

		/* Due to roots-locality considerations, roots (c,s)[0-15] are offset w.r.to the thread-local ptr pair as
		(cc0,ss0) + 0x[2,12,a,1a,6,16,e,1e,4,14,c,1c,8,18,10,20]. Here, due to the need to compute a new set of roots
		for each set of inputs, we use a streamlined sequence which computes only the [0,1,2,4,8,13]th roots with
		maximal accuracy (i.e. using 2-table-multiply), then generates the remaining ones from those. Thus the needed
		pointer offsets below are (cc0,ss0) + 0x[2,12,a,6,4,18]:
		*/
		c_tmp = cc0 + 0x02; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);

		// In SSE2 mode need 2 sets of sincos - since we may use different indexing patterns for storing k1,2_arr data
		// in AVX+ mode, don't share code for first 2 sets with those wider SIMD cases:
	  #if defined(USE_SSE2) && !defined(USE_AVX)

		// 1st set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[0] = k1<<4;	k2_arr[0] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[2] = k1<<4;	k2_arr[2] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[4] = k1<<4;	k2_arr[4] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[6] = k1<<4;	k2_arr[6] = k2<<4;
		l += (iroot << 2) + iroot;	/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[8] = k1<<4;	k2_arr[8] = k2<<4;

		// 2nd set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[1] = k1<<4;	k2_arr[1] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[3] = k1<<4;	k2_arr[3] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[5] = k1<<4;	k2_arr[5] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[7] = k1<<4;	k2_arr[7] = k2<<4;
		l += (iroot << 2) + iroot;	/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[9] = k1<<4;	k2_arr[9] = k2<<4;

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX16_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #elif defined(USE_AVX) && !defined(USE_AVX512)	// AVX/AVX2:

	  // In AVX mode need 4 sets of sincos:
		// 1st set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 0] = k1<<4;	k2_arr[ 0] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 4] = k1<<4;	k2_arr[ 4] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 8] = k1<<4;	k2_arr[ 8] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[12] = k1<<4;	k2_arr[12] = k2<<4;
		l += (iroot << 2) + iroot;	/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[16] = k1<<4;	k2_arr[16] = k2<<4;

		// 2nd set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 1] = k1<<4;	k2_arr[ 1] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 5] = k1<<4;	k2_arr[ 5] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 9] = k1<<4;	k2_arr[ 9] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[13] = k1<<4;	k2_arr[13] = k2<<4;
		l += (iroot << 2) + iroot;	/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[17] = k1<<4;	k2_arr[17] = k2<<4;

		// 3rd set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 2] = k1<<4;	k2_arr[ 2] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 6] = k1<<4;	k2_arr[ 6] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[10] = k1<<4;	k2_arr[10] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[14] = k1<<4;	k2_arr[14] = k2<<4;
		l += (iroot << 2) + iroot;	/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[18] = k1<<4;	k2_arr[18] = k2<<4;

		// 4th set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 3] = k1<<4;	k2_arr[ 3] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 7] = k1<<4;	k2_arr[ 7] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[11] = k1<<4;	k2_arr[11] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[15] = k1<<4;	k2_arr[15] = k2<<4;
		l += (iroot << 2) + iroot;	/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[19] = k1<<4;	k2_arr[19] = k2<<4;

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX16_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #elif defined(USE_AVX512)	// AVX512: need 8 sets of sincos gathered into 512-bit vectors-of-doubles
		/*
		Three options for the gather-four-128-bit-complex-double-data loads of our zmm-regs here:

		[1] In AVX512 the needed analog of VINSERTF128 is VINSERTF32x4, which loads 128-bits from the 2nd source operand
		(xmm/mem) into the desired 128-bit quartile of the destination zmm-reg, filling the other 3 quartiles from the
		1st source-operand (middle-operand) zmm-reg. Can do one load-low-128-bits and 3 VINSERTF32x4s for each zmm-reg.

		[2] Or, we can do 2 of the load-low-128-bit/VINSERTF128 pairs we use in the AVX case, each to fill a ymm-reg,
		then a VINSERTF64x4 to merge the 2 ymm-regs into a zmm-reg.

		[3] Or, we can use an AVX-512 gather-load (operands in reverse of GCC inline-asm order here):

			VGATHERDPD zmm1{k1},vm32y	Using signed qword indices, gather float64 data into zmm1 using OPMASK register
			k1 as completion mask.

		Description
			A set of 8 double-precision floating-point memory locations pointed by base address
			BASE_ADDR and index vector V_INDEX with scale SCALE are gathered. The result is
			written into a vector register. The elements are specified via the VSIB (i.e.,
			the index register vm32y is a vector register, holding packed indices). Elements
			will only be loaded if their corresponding mask bit in the OPMASK register is one.
			If an element’s mask bit is not set, the corresponding element of the destination
			register is left unchanged. The entire mask register will be set to zero by this
			instruction unless it triggers an exception.

		AT&T/GCC enhanced-inline-asm syntax example - note the required '%' escapes preceding
		each curly brace in the opmask-reg specifier:

			vgatherdpd 0x8(%%rax,%%ymm0),%%zmm1%{%%k1%}

		Our base-address here will be &rt0,1 stored in a GPR; the vector-offset data will be elements of the k1,2_arr arrays.
		How to efficiently load the latter - which are themselves vector-strided - into a zmm?
		IDEA: rearrange the k1,2_arr data so such (currently) uniform-vector-strided elements instead occupy adjacent slots,
		then simply load the needed block of elements from k1,2_arr - each half the width of the double it will generate the
		address of upon addition to the common base-address - into a ymm.
		*/
	  #ifdef USE_IMCI512	// Vectorized version below a no-go on 1st-gen Xeon Phi

		#warning AVX-512: Using slower non-ASM-macro version of sincos-indexing in radix16_dyadic_square.c.
		/*
		IMCI512: SSE2_RADIX16_CALC_TWIDDLES_LOACC reads k[0|1]_arr[0-39] in 8-index chunks,
			does gather-loads associated rt[0|1] elts, computes twiddles, writes to cc0+[0x4-...].
			Pull the gather-loads out here and do in C, just do the 512-bit vector CMULs in asm-macro:
		*/
		// 1st set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp = cc0+0x4; tm1 = cc0+0x5;	c_tmp = cc0+0x6; s_tmp = cc0+0x7;
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 2nd set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp = cc0+0x4; tm1 = cc0+0x5;	c_tmp = cc0+0x6; s_tmp = cc0+0x7;
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 3rd set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp = cc0+0x4; tm1 = cc0+0x5;	c_tmp = cc0+0x6; s_tmp = cc0+0x7;
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 4th set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp = cc0+0x4; tm1 = cc0+0x5;	c_tmp = cc0+0x6; s_tmp = cc0+0x7;
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 5th set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp = cc0+0x4; tm1 = cc0+0x5;	c_tmp = cc0+0x6; s_tmp = cc0+0x7;
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 6th set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp = cc0+0x4; tm1 = cc0+0x5;	c_tmp = cc0+0x6; s_tmp = cc0+0x7;
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 7th set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp = cc0+0x4; tm1 = cc0+0x5;	c_tmp = cc0+0x6; s_tmp = cc0+0x7;
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 8th set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp = cc0+0x4; tm1 = cc0+0x5;	c_tmp = cc0+0x6; s_tmp = cc0+0x7;
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

	  #else
		#warning Using AVX-512 code in radix16_dyadic_square.c

	   #if 0	// Set = 1 to enable timing-loop
		/*
			16 Jan 2017: 10^7-execs timings on single 1.3GHz core of KNL:
			Description																T(sec)
			1. Initial version of macro												0.740 [96 cycles]
			2. Use concat-regs to cut final indexing stage from 5 steps to 3		0.640 [83 cycles]
		*/
		/* Time-related stuff: *MUST* put the timing-loop clock() call inside this #if, since in production
		(non-timing) build mode it severely constrains the achievable parallelism: */
		clock_t clock1 = clock();
		double tdiff;
		for(l = 0; l < 10000000; l++) {	// Timing loop for dev-phase of code below
		/******* NOTE! Must comment out 'write updated values of index[0,1]_idx' in asm below when doing timings ******/
	   #endif
		// Vectorized form of C-code index computations below:
		__asm__ volatile (\
			"movq	%[__index0]    ,%%rax 	\n\t"\
			"movq	%[__index1]    ,%%rbx 	\n\t"\
			"movl	%[__index0_idx],%%ecx 	\n\t"\
			"movl	%[__index1_idx],%%edx 	\n\t"\
			"movl	%[__index1_mod],%%esi 	\n\t"\
			"movl	$1,%%edi				\n\t"\
			"vpbroadcastd	%%ecx,%%zmm0	\n\t"/* index0_idx x 16 ... only need x 8, but AVX-512F only supports full 512-bit vector width */\
			"vpbroadcastd	%%edx,%%zmm1	\n\t"/* index1_idx x 16 */\
			"vpbroadcastd	%%edi,%%zmm2	\n\t"/* 1 x 16 */\
			"vpbroadcastd	%%esi,%%zmm3	\n\t"/* index1_mod x 16 */\
		"movq	$0x0706050403020100,%%rsi	\n\t"/* 64-bit register w/byte offsets 0-7, bytes ordered left-to-right in decreasing significance */\
			"vmovq		%%rsi,%%xmm8 		\n\t"/* Copy byte pattern to low qword (64 bits) of ymm8 [NB: AVX-512 only supports MOVQ to/from 128-bit vector regs] */\
			"vpmovzxbd	%%xmm8,%%ymm8		\n\t"/* vector-index offsets: ymm8 = [0,1,2,3,4,5,6,7] in 32-bit form in low 8 dwords */\
			"vpaddd	%%ymm8,%%ymm1,%%ymm1	\n\t"/* index1_idx + [0-7] */\
			/* VPCMPD: Effect (zmm1 <= zmm3 ?) here via (zmm1 [not less than] zmm3 ?) ==> imm8 = 5: */\
			"vpcmpd	$5,%%zmm3,%%zmm1,%%k1	\n\t"/* (++index1_idx[zmm1] >= index1_mod[zmm3]) ? ... 8-fold version has index1_idx + [0-7] on LHS of the >= */\
			/* Writemask-version of padd/psub only supported at full register width in AVX-512F: */\
		"vpaddd	%%zmm2,%%zmm0,%%zmm0%{%%k1%}\n\t"/* if(++index1_idx >= index1_mod) ++index0_idx; */\
		"vpsubd	%%zmm3,%%zmm1,%%zmm1%{%%k1%}\n\t"/* if(++index1_idx >= index1_mod) index1_idx -= index1_mod; */\
			"movl	$-1,%%esi	\n\t"/* Init opmasks k1,2 (Only need the low byte of esi) */\
		"kmovw	%%esi,%%k1	\n\t	kmovw	%%esi,%%k2	\n\t"/* Masks k1,2 = 0x11111111 */\
			"vpgatherdd (%%rax,%%zmm0,4),%%zmm2%{%%k1%}	\n\t"/* index0[index0_idx] */\
			"vpgatherdd (%%rbx,%%zmm1,4),%%zmm3%{%%k2%}	\n\t"/* index1[index1_idx] */\
		"movl $13,%%edi \n\t vpbroadcastd %%edi,%%zmm5	\n\t"/* (uint32)13 x 16 ... only need x 8, but AVX-512F only supports full 512-bit vector width */\
			/* Move high dwords of ymm0,1 into low dword slots...due to inanity of Intel ISA, this needs 2 steps, first
			moving the high 128 bits into lo128, then the high 32 bits of that 128-bit chunk into the low dword slot: */\
			/* vextracti32x4 - AVX-512F only supports full 512-bit vector width (zmm -> xmm) version: */\
			"vextracti32x4 $1,%%zmm0,%%xmm0	\n\t	vextracti32x4 $1,%%zmm1,%%xmm1	\n\t"\
			"vpsrldq	$12,%%xmm0,%%xmm0	\n\t	vpsrldq	$12,%%xmm1,%%xmm1	\n\t"\
		/* ... and write updated values of index[0,1]_idx back to local vars: */\
		"vmovd	%%xmm0,%[__index0_idx] 	\n\t	vmovd	%%xmm1,%[__index1_idx] 	\n\t"\
			"vpaddd	%%ymm2,%%ymm3,%%ymm0	\n\t"/* l = iroot = index0[index0_idx] + index1[index1_idx] */\
		/*
		IDEA: Use full 512-bit vector width by placing the 8 1*iroot data in low half of zmm, then with 8 2*iroot data in ymm,
		copy the latter into upper half of zmm via VINSERTF32X8 $1,ymm,zmm,zmm. [Latency 3-6 cycles] **** better option available?
		Then:
		- mask/shift and write results into k1,2_arr[0-15];
		- compute zmm<<2 to get [4*iroot,8*iroot] in zmm, mask/shift yields k1,2_arr[16-31];
		- compute 13*ymm, mask/shift yields k1,2_arr[32-39].
		*/\
			"vpmulld %%ymm0,%%ymm5,%%ymm5	\n\t"/* ymm5 = 13*iroot */\
			"vpaddd	%%ymm0,%%ymm0,%%ymm4	\n\t"/* 2*iroot... */\
		"vinserti64x4 $1,%%ymm4,%%zmm0,%%zmm4	\n\t"/* Concatenate the [iroot,2*iroot] uint32-octets into a single zmm;
													use ...64x4 form here since ...32x8 needs AVX-512DQ enhancements */\
			"movl	%[__NRTM1]   ,%%esi 	\n\t"\
			"movl	%[__NRT_BITS],%%edi 	\n\t"\
			"movq	%[__k1_arr]  ,%%rax 	\n\t"\
			"movq	%[__k2_arr]  ,%%rbx 	\n\t"/* NB: k2_arr only 32-byte aligned, hence use vmovups for 64-byte (rbx) stores below */\
			"vpbroadcastd	%%esi,%%zmm2	\n\t"/* NRTM1    x 16 */\
			"vpbroadcastd	%%edi,%%zmm3	\n\t"/* NRT_BITS x 16 */\
			/* [1,2]*iroot: */\
			"vpandd	%%zmm2,%%zmm4,%%zmm0	\n\t"/* k1=(l & NRTM1) */\
			"vpsrlvd %%zmm3,%%zmm4,%%zmm1	\n\t"/* k2=(l >> NRT_BITS) */\
			"vpslld	$4,%%zmm0,%%zmm0		\n\t	vpslld	$4,%%zmm1,%%zmm1	\n\t"/* k1,2<<4 */\
			"vmovaps	%%zmm0,0x00(%%rax)	\n\t	vmovups	%%zmm1,0x00(%%rbx)	\n\t"/* Write k1,2_arr[ 0-15] = k1,2<<4 */\
			/* [4,8]*iroot: */\
			"vpslld	$2,%%zmm4,%%zmm4		\n\t"/* l <<= 2; result = [4,8]*iroot */\
			"vpandd	%%zmm2,%%zmm4,%%zmm0	\n\t"\
			"vpsrlvd %%zmm3,%%zmm4,%%zmm1	\n\t"\
			"vpslld	$4,%%zmm0,%%zmm0		\n\t	vpslld	$4,%%zmm1,%%zmm1	\n\t"\
			"vmovaps	%%zmm0,0x40(%%rax)	\n\t	vmovups	%%zmm1,0x40(%%rbx)	\n\t"/* Write k1,2_arr[16-31] */\
			/* 13*iroot, using half-register-width operands: */\
			"vpand	%%ymm2,%%ymm5,%%ymm0	\n\t"\
			"vpsrlvd %%ymm3,%%ymm5,%%ymm1	\n\t"\
			"vpslld	$4,%%ymm0,%%ymm0		\n\t	vpslld	$4,%%ymm1,%%ymm1	\n\t"\
			"vmovaps	%%ymm0,0x80(%%rax)	\n\t	vmovaps	%%ymm1,0x80(%%rbx)	\n\t"/* Write k1,2_arr[32-39] */\
		:	: [__iroot] "m" (iroot)	/* No outputs; All inputs from memory addresses here */\
			 ,[__NRT_BITS] "m" (NRT_BITS)	\
			 ,[__NRTM1] "m" (NRTM1)	\
			 ,[__index0] "m" (index0)	\
			 ,[__index1] "m" (index1)	\
			 ,[__index0_idx] "m" (index0_idx)	\
			 ,[__index1_idx] "m" (index1_idx)	\
			 ,[__index1_mod] "m" (index1_mod)	\
			 ,[__k1_arr] "m" (k1_arr)	\
			 ,[__k2_arr] "m" (k2_arr)	\
			: "cc","memory","cl","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	/* Clobbered registers */\
		);
	   #if 0	// Set = 1 to enable timing-loop
		}	tdiff = (double)(clock() - clock1);	printf("Time for %u 8x8 execs of AVX-512 twiddle-index macro =%s\n",l,get_time_str(tdiff));	exit(0);
	   #endif

		// Need one more (scalar-mode) update of these in preparation for the next 8-fold chunk:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }

	  #endif	// (IMCI512 or AVX512?) toggle

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX16_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #endif	// SIMD mode?

	#endif	// SIMD ?

/*   Radix-16 DIF FFT: gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms... */

	/* In AVX512 mode, the input data are arranged in memory like so, where we view things in 64-byte chunks:

		&a[j1]	a0-7.re	...+32]	b0-7.re	...+64]	c0-7.re	...+96]	d0-7.re	...+128	e0-7.re	...+160	f0-7.re	...+192	g0-7.re	...+224	h0-7.re
		+0x040	a0-7.im	 +0x040	b0-7.im	 +0x040	c0-7.im	 +0x040	d0-7.im	 +0x040	e0-7.im	 +0x040	f0-7.im	 +0x040	g0-7.im	 +0x040	h0-7.im
		+0x080	a8-f.re	 +0x080	b8-f.re	 +0x080	c8-f.re	 +0x080	d8-f.re	 +0x080	e8-f.re	 +0x080	f8-f.re	 +0x080	g8-f.re	 +0x080	h8-f.re
		+0x0c0	a8-f.im	 +0x0c0	b8-f.im	 +0x0c0	c8-f.im	 +0x0c0	d8-f.im	 +0x0c0	e8-f.im	 +0x0c0	f8-f.im	 +0x0c0	g8-f.im	 +0x0c0	h8-f.im
	for a total of 16 vector-complex data (= 256 float64), 2 in each a[j1+...] column.
	Following the same outline as the AVX case described in more detail below this leads to an interleaving consisting
	of four 8x8 matrix-of-float64 transpositions, one pair for the 0-7.re,im data, another pair for the 8-f.re,im data.
	*/

	/* In AVX mode, the input data are arranged in memory like so, where we view things in 32-byte chunks:

		&a[j1]	a0.re,a1.re,a2.re,a3.re	&a[j1+32]	b0.re,b1.re,b2.re,b3.re	 &a[j1+64]	c0.re,c1.re,c2.re,c3.re	 &a[j1+96]	d0.re,d1.re,d2.re,d3.re
		+0x020	a0.im,a1.im,a2.im,a3.im		+0x020	b0.im,b1.im,b2.im,b3.im		+0x020	c0.im,c1.im,c2.im,c3.im		+0x020	d0.im,d1.im,d2.im,d3.im
		+0x040	a4.re,a5.re,a6.re,a7.re		+0x040	b4.re,b5.re,b6.re,b7.re		+0x040	c4.re,c5.re,c6.re,c7.re		+0x040	d4.re,d5.re,d6.re,d7.re
		+0x060	a4.im,a5.im,a6.im,a7.im		+0x060	b4.im,b5.im,b6.im,b7.im		+0x060	c4.im,c5.im,c6.im,c7.im		+0x060	d4.im,d5.im,d6.im,d7.im
		+0x080	a8.re,a9.re,aA.re,aB.re		+0x080	b8.re,b9.re,bA.re,bB.re		+0x080	c8.re,c9.re,cA.re,cB.re		+0x080	d8.re,d9.re,dA.re,dB.re
		+0x0a0	a8.im,a9.im,aA.im,aB.im		+0x0a0	b8.im,b9.im,bA.im,bB.im		+0x0a0	c8.im,c9.im,cA.im,cB.im		+0x0a0	d8.im,d9.im,dA.im,dB.im
		+0x0c0	aC.re,aD.re,aE.re,aF.re		+0x0c0	bC.re,bD.re,bE.re,bF.re		+0x0c0	cC.re,cD.re,cE.re,cF.re		+0x0c0	dC.re,dD.re,dE.re,dF.re
		+0x0e0	aC.im,aD.im,aE.im,aF.im		+0x0e0	bC.im,bD.im,bE.im,bF.im		+0x0e0	cC.im,cD.im,cE.im,cF.im		+0x0e0	dC.im,dD.im,dE.im,dF.im
	for a total of 16 vector-complex data (= 128 float64), 4 in each a[j1+...] column.

	Thus, we have 16 input complex data quartets, whose subelements need to be properly interleaved prior to the radix-16 complex DFT.

	We can break this into a bunch of smaller interleave steps, in each of which have 4 register quartets-of-doubles
	a0-3, b0-3, c0-3, d0-3, which we need to interleave like so (= a 4x4 matrix transposition, if one considers each
	input register to hold a row of the matrix):

		[a0,a1,a2,a3]         |\     [a0,b0,c0,d0]
		[b0,b1,b2,b3]    -----  \    [a1,b1,c1,d1]
		[c0,c1,c2,c3]    -----  /    [a2,b2,c2,d2]
		[d0,d1,d2,d3]         |/     [a3,b3,c3,d3]

	George says he uses this sequence for this transposition operation:

		vmovapd	ymm1, [srcreg]				;; R1											a0,a1,a2,a3
		vmovapd	ymm7, [srcreg+d1]			;; R2											b0,b1,b2,b3
		vshufpd	ymm0, ymm1, ymm7, 15		;; Shuffle R1 and R2 to create R1/R2 hi			ymm0:a1,b1,a3,b3
		vshufpd	ymm1, ymm1, ymm7, 0			;; Shuffle R1 and R2 to create R1/R2 low		ymm1:a0,b0,a2,b2

		vmovapd	ymm2, [srcreg+d2]			;; R3											c0,c1,c2,c3
		vmovapd	ymm7, [srcreg+d2+d1]		;; R4											d0,d1,d2,d3
		vshufpd	ymm3, ymm2, ymm7, 15		;; Shuffle R3 and R4 to create R3/R4 hi			ymm3:c1,d1,c3,d3
		vshufpd	ymm2, ymm2, ymm7, 0			;; Shuffle R3 and R4 to create R3/R4 low		ymm2:c0,d0,c2,d2

		vperm2f128 ymm4, ymm0, ymm3, 32		;; Shuffle R1/R2 hi and R3/R4 hi (new R2)		ymm4:[a1,b1][c1,d1]	imm0:1=0,4:5=2 -> [in1.lo,in2.lo] = [ymm0.lo,ymm3.lo]
		vperm2f128 ymm0, ymm0, ymm3, 49		;; Shuffle R1/R2 hi and R3/R4 hi (new R4)		ymm0:[a3,b3][c3,d3]	imm0:1=1,4:5=3 -> [in1.hi,in2.hi] = [ymm0.hi,ymm3.hi]

		vperm2f128 ymm3, ymm1, ymm2, 32		;; Shuffle R1/R2 low and R3/R4 low (new R1)		ymm3:[a0,b0][c0,d0]	imm0:1=0,4:5=2 -> [in1.lo,in2.lo] = [ymm1.lo,ymm2.lo]
		vperm2f128 ymm1, ymm1, ymm2, 49		;; Shuffle R1/R2 low and R3/R4 low (new R3)		ymm1:[a2,b2][c2,d2]	imm0:1=1,4:5=3 -> [in1.hi,in2.hi] = [ymm1.hi,ymm2.hi]

	I want output-reg index to match subscripts of final a-d set stored in the register, so swap ymm indices 0134 -> 3201:

		// Data exit in ymmA-D; ymmX,Y are any 2 other registers
		vmovapd	ymmC, [srcreg]				// a0,a1,a2,a3
		vmovapd	ymmX, [srcreg+d1]			// b0,b1,b2,b3
		vshufpd	ymmD, ymmC, ymmX, 15		// ymmD:a1,b1,a3,b3
		vshufpd	ymmC, ymmC, ymmX, 0			// ymmC:a0,b0,a2,b2

		vmovapd	ymmY, [srcreg+d2]			// c0,c1,c2,c3
		vmovapd	ymmX, [srcreg+d2+d1]		// d0,d1,d2,d3
		vshufpd	ymmA, ymmY, ymmX, 15		// ymmA:c1,d1,c3,d3
		vshufpd	ymmY, ymmY, ymmX, 0			// ymmY:c0,d0,c2,d2

		vperm2f128 ymmB, ymmD, ymmA, 32		// ymmB:[a1,b1][c1,d1]
		vperm2f128 ymmD, ymmD, ymmA, 49		// ymmD:[a3,b3][c3,d3]
		vperm2f128 ymmA, ymmC, ymmY, 32		// ymmA:[a0,b0][c0,d0]
		vperm2f128 ymmC, ymmC, ymmY, 49		// ymmC:[a2,b2][c2,d2]

	If we were writing the shuffled data quartets back to main-array, we would do it like so:

		&a[j1]	a0.re,a1.re,a2.re,a3.re	&a[j1+32]	b0.re,b1.re,b2.re,b3.re	 &a[j1+64]	c0.re,c1.re,c2.re,c3.re	 &a[j1+96]	d0.re,d1.re,d2.re,d3.re
		+0x020	a0.im,a1.im,a2.im,a3.im		+0x020	b0.im,b1.im,b2.im,b3.im		+0x020	c0.im,c1.im,c2.im,c3.im		+0x020	d0.im,d1.im,d2.im,d3.im

		==>		a0.re,b0.re,c0.re,d0.re				a1.re,b1.re,c1.re,d1.re				a2.re,b2.re,c2.re,d2.re				a3.re,b3.re,c3.re,d3.re
				a0.im,b0.im,c0.im,d0.im				a1.im,b1.im,c1.im,d1.im				a2.im,b2.im,c2.im,d2.im				a3.im,b3.im,c3.im,d3.im
				[radix-4 DFT Block #1]				[radix-4 DFT Block #3]				[radix-4 DFT Block #2]				[radix-4 DFT Block #4]
	----------
		+0x040	a4.re,a5.re,a6.re,a7.re		+0x040	b4.re,b5.re,b6.re,b7.re		+0x040	c4.re,c5.re,c6.re,c7.re		+0x040	d4.re,d5.re,d6.re,d7.re
		+0x060	a4.im,a5.im,a6.im,a7.im		+0x060	b4.im,b5.im,b6.im,b7.im		+0x060	c4.im,c5.im,c6.im,c7.im		+0x060	d4.im,d5.im,d6.im,d7.im

		==>		a4.re,b4.re,c4.re,d4.re				a5.re,b5.re,c5.re,d5.re				a6.re,b6.re,c6.re,d6.re				a7.re,b7.re,c7.re,d7.re
				a4.im,b4.im,c4.im,d4.im				a5.im,b5.im,c5.im,d5.im				a6.im,b6.im,c6.im,d6.im				a7.im,b7.im,c7.im,d7.im
				[radix-4 DFT Block #1]				[radix-4 DFT Block #3]				[radix-4 DFT Block #2]				[radix-4 DFT Block #4]
	----------
		+0x080	a8.re,a9.re,aA.re,aB.re		+0x080	b8.re,b9.re,bA.re,bB.re		+0x080	c8.re,c9.re,cA.re,cB.re		+0x080	d8.re,d9.re,dA.re,dB.re
		+0x0a0	a8.im,a9.im,aA.im,aB.im		+0x0a0	b8.im,b9.im,bA.im,bB.im		+0x0a0	c8.im,c9.im,cA.im,cB.im		+0x0a0	d8.im,d9.im,dA.im,dB.im

		==>		a8.re,b8.re,c8.re,d8.re				a9.re,b9.re,c9.re,d9.re				aA.re,bA.re,cA.re,dA.re				aB.re,bB.re,cB.re,dB.re
				a8.im,b8.im,c8.im,d8.im				a9.im,b9.im,c9.im,d9.im				aA.im,bA.im,cA.im,dA.im				aB.im,bB.im,cB.im,dB.im
				[radix-4 DFT Block #1]				[radix-4 DFT Block #3]				[radix-4 DFT Block #2]				[radix-4 DFT Block #4]
	----------
		+0x0c0	aC.re,aD.re,aE.re,aF.re		+0x0c0	bC.re,bD.re,bE.re,bF.re		+0x0c0	cC.re,cD.re,cE.re,cF.re		+0x0c0	dC.re,dD.re,dE.re,dF.re
		+0x0e0	aC.im,aD.im,aE.im,aF.im		+0x0e0	bC.im,bD.im,bE.im,bF.im		+0x0e0	cC.im,cD.im,cE.im,cF.im		+0x0e0	dC.im,dD.im,dE.im,dF.im

		==>		aC.re,bC.re,cC.re,dC.re				aD.re,bD.re,cD.re,dD.re				aE.re,bE.re,cE.re,dE.re				aF.re,bF.re,cF.re,dF.re
				aC.im,bC.im,cC.im,dC.im				aD.im,bD.im,cD.im,dD.im				aE.im,bE.im,cE.im,dE.im				aF.im,bF.im,cF.im,dF.im
				[radix-4 DFT Block #1]				[radix-4 DFT Block #3]				[radix-4 DFT Block #2]				[radix-4 DFT Block #4]

	But we prefer to write the outputs into the same local-mem storage used by the ensuing radix-4/16 DFTs.
	Those 4 radix-4 DFTs use the data underscored above as [radix-4 DFT block X], thus all data belonging to the same
	block get written to a contiguous block of local memory: Blocks 1-4 go into addresses starting at r1+0,0x80,0x100,0x180.
	*/

	/* In SSE2 mode, the input data are arranged in memory like so, where we view things in 16-byte chunks:

		&a[j1]	a0.re,a1.re  &a[j1+32]	b0.re,b1.re
		+0x010	a0.im,a1.im		+0x010	b0.im,b1.im
		+0x020	a2.re,a3.re		+0x020	b2.re,b3.re
		+0x030	a2.im,a3.im		+0x030	b2.im,b3.im
		+0x040	a4.re,a5.re		+0x040	b4.re,b5.re
		+0x050	a4.im,a5.im		+0x050	b4.im,b5.im
		+0x060	a6.re,a7.re		+0x060	b6.re,b7.re
		+0x070	a6.im,a7.im		+0x070	b6.im,b7.im
		+0x080	a8.re,a9.re		+0x080	b8.re,b9.re
		+0x090	a8.im,a9.im		+0x090	b8.im,b9.im
		+0x0a0	aA.re,aB.re		+0x0a0	bA.re,bB.re
		+0x0b0	aA.im,aB.im		+0x0b0	bA.im,bB.im
		+0x0c0	aC.re,aD.re		+0x0c0	bC.re,bD.re
		+0x0d0	aC.im,aD.im		+0x0d0	bC.im,bD.im
		+0x0e0	aE.re,aF.re		+0x0e0	bE.re,bF.re
		+0x0f0	aE.im,aF.im		+0x0f0	bE.im,bF.im, for a total of 16 vector-complex data (= 64 float64),
													8 vector-complex (= 32 float64) in each a[j1+...] column.

	We need to interleave these pairwise so as to swap the high double of each A-pair with the low double of the corresponding B-pair, e.g

		low		high	low		high
		[a0.re,a1.re]	[b0.re,b1.re]
		   |      \       /      |
		   |        \   /        |
		   |          x          |
		   |        /   \        |
		   V      /       \      V
		[a0.re,b0.re]	[a1.re,b1.re]
	=	[ A.lo, B.lo]	[ A.hi, B.hi]
			(1)				(2)

	The instructions needed for this permutation [assuming A-term in memA, B-term in memB] is

		(0)		movaps		xmm0,memA
		(1)		movaps		xmm1,memA
		(1)		unpcklpd	xmm0,xmmB		[xmm0.lo,xmm0.hi] <- [A.lo,B.lo]
		(2)		unpckhpd	xmm1,xmmB		[xmm1.lo,xmm1.hi] <- [A.hi,B.hi] ,

	Alternatively, we can use the SHUFPD instruction.
	SHUFPD xmm,mem,imm produces [xmm.lo,xmm.hi] <- [xmm.(lo or hi),xmm.(lo or hi), depending on the value of imm:

		IMM:
		0		shufpd	xmm,mem		[xmm.lo,xmm.hi] <- [xmm.lo,mem.lo]
		1		shufpd	xmm,mem		[xmm.lo,xmm.hi] <- [xmm.hi,mem.lo]
		2		shufpd	xmm,mem		[xmm.lo,xmm.hi] <- [xmm.lo,mem.hi]
		3		shufpd	xmm,mem		[xmm.lo,xmm.hi] <- [xmm.hi,mem.hi]

	So what we need is

		(0)		movaps		xmm0,memA
		(1)		movaps		xmm1,memA
		(2)		shufpd		xmm0,memB,0		[xmm0.lo,xmm0.hi] <- [A.lo,B.lo]
		(3)		shufpd		xmm1,memB,3		[xmm1.lo,xmm1.hi] <- [A.hi,B.hi]

	It's not clear whether there is a preference for one or the other instruction sequence based on resulting performance.

	It is useful to visualize the interleaving which occurs for the various radix-4 sub-DFTs (labeled [1]-[4]),
	in terms of inputs from the main data array and outputs into the local storage. We do this for 2 consecutive
	set of data (i.e. 2 loop passes), to show the same data which will be getting processed in 1 loop pass in AVX mode:

						Pass 1 INPUTS:														Pass 2 INPUTS:
			add0:					add1:											add2:					add3:
	[1re]	a + 0	a0.re,a1.re		a + 32	b0.re,b1.re	[1*,3]				[1re]	a + 64	c0.re,c1.re		a + 96	d0.re,d1.re	[1*,3]
	[1im]	+0x010	a0.im,a1.im		+0x010	b0.im,b1.im	[1*,3]				[1im]	+0x010	c0.im,c1.im		+0x010	d0.im,d1.im	[1*,3]
	[2re]	+0x020	a2.re,a3.re		+0x020	b2.re,b3.re		[2*,4]			[2re]	+0x020	c2.re,c3.re		+0x020	d2.re,d3.re		[2*,4]
	[2im]	+0x030	a2.im,a3.im		+0x030	b2.im,b3.im		[2*,4]			[2im]	+0x030	c2.im,c3.im		+0x030	d2.im,d3.im		[2*,4]
	[1re]	+0x040	a4.re,a5.re		+0x040	b4.re,b5.re	[1*,3]				[1re]	+0x040	c4.re,c5.re		+0x040	d4.re,d5.re	[1*,3]
	[1im]	+0x050	a4.im,a5.im		+0x050	b4.im,b5.im	[1*,3]				[1im]	+0x050	c4.im,c5.im		+0x050	d4.im,d5.im	[1*,3]
	[2re]	+0x060	a6.re,a7.re		+0x060	b6.re,b7.re		[2*,4]			[2re]	+0x060	c6.re,c7.re		+0x060	d6.re,d7.re		[2*,4]
	[2im]	+0x070	a6.im,a7.im		+0x070	b6.im,b7.im		[2*,4]			[2im]	+0x070	c6.im,c7.im		+0x070	d6.im,d7.im		[2*,4]
	[1re]	+0x080	a8.re,a9.re		+0x080	b8.re,b9.re	[1*,3]				[1re]	+0x080	c8.re,c9.re		+0x080	d8.re,d9.re	[1*,3]
	[1im]	+0x090	a8.im,a9.im		+0x090	b8.im,b9.im	[1*,3]				[1im]	+0x090	c8.im,c9.im		+0x090	d8.im,d9.im	[1*,3]
	[2re]	+0x0a0	aA.re,aB.re		+0x0a0	bA.re,bB.re		[2*,4]			[2re]	+0x0a0	cA.re,cB.re		+0x0a0	dA.re,dB.re		[2*,4]
	[2im]	+0x0b0	aA.im,aB.im		+0x0b0	bA.im,bB.im		[2*,4]			[2im]	+0x0b0	cA.im,cB.im		+0x0b0	dA.im,dB.im		[2*,4]
	[1re]	+0x0c0	aC.re,aD.re		+0x0c0	bC.re,bD.re	[1*,3]				[1re]	+0x0c0	cC.re,cD.re		+0x0c0	dC.re,dD.re	[1*,3]
	[1im]	+0x0d0	aC.im,aD.im		+0x0d0	bC.im,bD.im	[1*,3]				[1im]	+0x0d0	cC.im,cD.im		+0x0d0	dC.im,dD.im	[1*,3]
	[2re]	+0x0e0	aE.re,aF.re		+0x0e0	bE.re,bF.re		[2*,4]			[2re]	+0x0e0	cE.re,cF.re		+0x0e0	dE.re,dF.re		[2*,4]
	[2im]	+0x0f0	aE.im,aF.im		+0x0f0	bE.im,bF.im		[2*,4]			[2im]	+0x0f0	cE.im,cF.im		+0x0f0	dE.im,dF.im		[2*,4]

						Pass 1 OUTPUTS:															Pass 2 OUTPUTS:
	[1re]	r0+  0	a0.re,b0.re		r0+x100	a1.re,b1.re	[1*,3]				[1re]	r0+x200	c0.re,d0.re		r0+x300	c1.re,d1.re	[1*,3]
	[1im]	+0x010	a0.im,b0.im		+0x010	a1.im,b1.im	[1*,3]				[1im]	+0x010	c0.im,d0.im		+0x010	c1.im,d1.im	[1*,3]
	[2re]	+0x020	a2.re,b2.re		+0x020	a3.re,b3.re		[2*,4]			[2re]	+0x020	c2.re,d2.re		+0x020	c3.re,d3.re		[2*,4]
	[2im]	+0x030	a2.im,b2.im		+0x030	a3.im,b3.im		[2*,4]			[2im]	+0x030	c2.im,d2.im		+0x030	c3.im,d3.im		[2*,4]
	[1re]	+0x040	a4.re,b4.re		+0x040	a5.re,b5.re	[1*,3]				[1re]	+0x040	c4.re,d4.re		+0x040	c5.re,d5.re	[1*,3]
	[1im]	+0x050	a4.im,b4.im		+0x050	a5.im,b5.im	[1*,3]				[1im]	+0x050	c4.im,d4.im		+0x050	c5.im,d5.im	[1*,3]
	[2re]	+0x060	a6.re,b6.re		+0x060	a7.re,b7.re		[2*,4]			[2re]	+0x060	c6.re,d6.re		+0x060	c7.re,d7.re		[2*,4]
	[2im]	+0x070	a6.im,b6.im		+0x070	a7.im,b7.im		[2*,4]			[2im]	+0x070	c6.im,d6.im		+0x070	c7.im,d7.im		[2*,4]
	[1re]	+0x080	a8.re,b8.re		+0x080	a9.re,b9.re	[1*,3]				[1re]	+0x080	c8.re,d8.re		+0x080	c9.re,d9.re	[1*,3]
	[1im]	+0x090	a8.im,b8.im		+0x090	a9.im,b9.im	[1*,3]				[1im]	+0x090	c8.im,d8.im		+0x090	c9.im,d9.im	[1*,3]
	[2re]	+0x0a0	aA.re,bA.re		+0x0a0	aB.re,bB.re		[2*,4]			[2re]	+0x0a0	cA.re,dA.re		+0x0a0	cB.re,dB.re		[2*,4]
	[2im]	+0x0b0	aA.im,bA.im		+0x0b0	aB.im,bB.im		[2*,4]			[2im]	+0x0b0	cA.im,dA.im		+0x0b0	cB.im,dB.im		[2*,4]
	[1re]	+0x0c0	aC.re,bC.re		+0x0c0	aD.re,bD.re	[1*,3]				[1re]	+0x0c0	cC.re,dC.re		+0x0c0	cD.re,dD.re	[1*,3]
	[1im]	+0x0d0	aC.im,bC.im		+0x0d0	aD.im,bD.im	[1*,3]				[1im]	+0x0d0	cC.im,dC.im		+0x0d0	cD.im,dD.im	[1*,3]
	[2re]	+0x0e0	aE.re,bE.re		+0x0e0	aF.re,bF.re		[2*,4]			[2re]	+0x0e0	cE.re,dE.re		+0x0e0	cF.re,dF.re		[2*,4]
	[2im]	+0x0f0	aE.im,bE.im		+0x0f0	aF.im,bF.im		[2*,4]			[2im]	+0x0f0	cE.im,dE.im		+0x0f0	cF.im,dF.im		[2*,4]

	Here, e.g. [1*,3] means "interleaved with Block 1 inputs, high parts stored in local mem and used for Block 3 radix-4 DFT".

Now compare with AVX mode and its data quartets:
																		INPUTS:
			add0:								add0+0x100:							add1:								add1+0x100:
	[1-4re]	a +  0	a0.re,a1.re,a2.re,a3.re		a + 32	b0.re,b1.re,b2.re,b3.re		a + 64	c0.re,c1.re,c2.re,c3.re		a + 96	d0.re,d1.re,d2.re,d3.re
	[1-4im]	+0x020	a0.im,a1.im,a2.im,a3.im		+0x020	b0.im,b1.im,b2.im,b3.im		+0x020	c0.im,c1.im,c2.im,c3.im		+0x020	d0.im,d1.im,d2.im,d3.im
	[1-4re]	+0x040	a4.re,a5.re,a6.re,a7.re		+0x040	b4.re,b5.re,b6.re,b7.re		+0x040	c4.re,c5.re,c6.re,c7.re		+0x040	d4.re,d5.re,d6.re,d7.re
	[1-4im]	+0x060	a4.im,a5.im,a6.im,a7.im		+0x060	b4.im,b5.im,b6.im,b7.im		+0x060	c4.im,c5.im,c6.im,c7.im		+0x060	d4.im,d5.im,d6.im,d7.im
	[1-4re]	+0x080	a8.re,a9.re,aA.re,aB.re		+0x080	b8.re,b9.re,bA.re,bB.re		+0x080	c8.re,c9.re,cA.re,cB.re		+0x080	d8.re,d9.re,dA.re,dB.re
	[1-4im]	+0x0a0	a8.im,a9.im,aA.im,aB.im		+0x0a0	b8.im,b9.im,bA.im,bB.im		+0x0a0	c8.im,c9.im,cA.im,cB.im		+0x0a0	d8.im,d9.im,dA.im,dB.im
	[1-4re]	+0x0c0	aC.re,aD.re,aE.re,aF.re		+0x0c0	bC.re,bD.re,bE.re,bF.re		+0x0c0	cC.re,cD.re,cE.re,cF.re		+0x0c0	dC.re,dD.re,dE.re,dF.re
	[1-4im]	+0x0e0	aC.im,aD.im,aE.im,aF.im		+0x0e0	bC.im,bD.im,bE.im,bF.im		+0x0e0	cC.im,cD.im,cE.im,cF.im		+0x0e0	dC.im,dD.im,dE.im,dF.im

Notice how the data used in all 4 sub-DFTs are all mixed together on input.
Taking the 4 quartets in the first row of the above and 4-way-interleaving them yields

	a0.re,b0.re,c0.re,d0.re	[1re]
	a1.re,b1.re,c1.re,d1.re	[3re]
	a2.re,b2.re,c2.re,d2.re	[2re]
	a3.re,b3.re,c3.re,d3.re	[4re]

The 0 and 2-set of which are used in the Block 1 and 2 4-DFTs, and the 1 and 3-sets in the Block 3 and 4 4-DFTs.
So this address pattern is rather different than in the SSE2 2-way case
Thus we should simply combine the separate 2-way interleaving steps in the side-by-side instruction streams of the 64-bit
inline asm code into 4-way interleaves in the AVX version, with outputs dispatched to the same registers and local-store
locations, just with local-store address offsets doubled due to the double data width, obviously.
*/

#ifdef USE_SSE2	// SSE2 and AVX share same code flow, just with different versions of the ASM compute macros

		add0 = &a[j1];
		add1 = &a[j1+stridh];	// stridh = 128/64/32/16 for AVX512/AVX/SSE2/scalar-double, respectively
	#ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX mode
	  					// because (add[1,3,5,7]-add[0,2,4,6]) have opposite signs for Fermat and Mersenne-mod:
		add1 = add0 +  32;
		add2 = add0 +  64;
		add3 = add0 +  96;
		add4 = add0 + 128;
		add5 = add0 + 160;
		add6 = add0 + 192;
		add7 = add0 + 224;
	#elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  					// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:
		add1 = add0 + 32;
		add2 = add0 + 64;
		add3 = add0 + 96;
	#endif

	#if FULLY_FUSED

		SSE2_RADIX16_DIF_DYADIC_DIT(add0,add1,r1,isrt2,pfetch_dist);

	#else

	  if(fwd_fft_only != 3)	// v20: add support for both-inputs-already-fwd-FFTed case - fwd_fft_only == 3 means skip fwd-FFT
	  {
	  #ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX mode
	  					// because (add[1,3,5,7]-add[0,2,4,6]) have opposite signs for Fermat and Mersenne-mod:
		// process 8 main-array blocks of 4 vec_dbl [= 4 x 8 = 32 doubles each] in AVX512 mode, total = 256 float64
		SSE2_RADIX16_WRAPPER_DIF(add0,add1,add2,add3,add4,add5,add6,add7
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,pfetch_dist)
	  #elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  					// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:
		// process 4 main-array blocks of 8 vec_dbl [= 8 x 4 = 32 doubles each] in AVX mode, total = 128 float64
		SSE2_RADIX16_WRAPPER_DIF(add0,add1,add2,add3
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,pfetch_dist)
	  #else	// SSE2:
		SSE2_RADIX16_WRAPPER_DIF(add0,add1
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)
	  #endif
	  }

	  // v19: If fwd_fft_only = 1, write fwd-FFT result back to input array, skipping dyadic-square and inv-FFT steps:
	  if(fwd_fft_only == 1)
	  {
		nbytes = ((intptr_t)r17 - (intptr_t)r1)<<1;
		memcpy(add0, r1, nbytes);	// add0 = a + j1pad;
		goto loop;	// Skip dyadic-mul and iFFT

	  } else if (fwd_fft_only) {	// 2-input modmul: fwdFFT data dyadic-mul'ed with precomputed 2nd-vector stored in fwdFFTed form at address stored in (uint64)-cast form in fwd_fft_only

	  if(fwd_fft_only == 3) {	// v20: Both inputs enter fwd-FFTed, must copy from main arrays a[],b[] to local-vars
		nbytes = ((intptr_t)r17 - (intptr_t)r1)<<1;
		memcpy(r1, add0, nbytes);	// add0 = a + j1pad;
	  }
	  #ifndef USE_ARM_V8_SIMD
		bdd0 = b + j1;	// No reverse-running addresses as in Mers-mod/real-vector-FFT case, just need single B-array base address
		cdd0 = c_arr + j1;
	  #endif

		if(c_arr) {	// c_arr != 0x0: a * (b - c)

		#ifdef USE_ARM_V8_SIMD

		// No Fermat-mod support on ARMv8, just supply a stub macro:
		  __asm__ volatile (\
			"ldr x0,%[__r1]	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			: "cc","memory","x0"	/* Clobbered registers */\
		  );

		#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

		  __asm__ volatile (\
			"movq	%[__r1],%%rax		\n\t"\
			"movq	%[__bdd0],%%rbx		\n\t"\
			"movq	%[__cdd0],%%rcx		\n\t"\
			/* x0.y0: */							/* x8.y8: */\
			"vmovaps		 (%%rbx),%%zmm2	\n\t	vmovaps	0x400(%%rbx),%%zmm7	\n\t"/* y.re */\
			"vmovaps	0x040(%%rbx),%%zmm3	\n\t	vmovaps	0x440(%%rbx),%%zmm8	\n\t"/* y.im */\
		"vsubpd		 (%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x400(%%rcx),%%zmm7,%%zmm7	\n\t"/* y.re -= z.re */\
		"vsubpd	0x040(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x440(%%rcx),%%zmm8,%%zmm8	\n\t"/* y.im -= z.im */\
			"vmovaps		 (%%rax),%%zmm0	\n\t	vmovaps	0x400(%%rax),%%zmm5	\n\t"/* x.re */\
			"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x440(%%rax),%%zmm6	\n\t"/* x.im */\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"/* copy x.re */\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"/* x.re *= y.re */\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"/* y.re *= x.im */\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"/* x.im *= y.im */\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"/* y.im *= x.re[copy] */\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm5,0x400(%%rax)	\n\t"/* write z.re */\
			"vmovaps	%%zmm2,0x040(%%rax)	\n\t	vmovaps	%%zmm7,0x440(%%rax)	\n\t"/* write z.im */\
			/* x1.y1: */							/* x9.y9: */\
			"vmovaps	0x080(%%rbx),%%zmm2	\n\t	vmovaps	0x480(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x0c0(%%rbx),%%zmm3	\n\t	vmovaps	0x4c0(%%rbx),%%zmm8	\n\t"\
		"vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x480(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x4c0(%%rcx),%%zmm8,%%zmm8	\n\t"\
			"vmovaps	0x080(%%rax),%%zmm0	\n\t	vmovaps	0x480(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x0c0(%%rax),%%zmm1	\n\t	vmovaps	0x4c0(%%rax),%%zmm6	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x080(%%rax)	\n\t	vmovaps	%%zmm5,0x480(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x0c0(%%rax)	\n\t	vmovaps	%%zmm7,0x4c0(%%rax)	\n\t"\
			/* x2.y2: */							/* xA.yA: */\
			"vmovaps	0x100(%%rbx),%%zmm2	\n\t	vmovaps	0x500(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x140(%%rbx),%%zmm3	\n\t	vmovaps	0x540(%%rbx),%%zmm8	\n\t"\
		"vsubpd	0x100(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x500(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	0x140(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x540(%%rcx),%%zmm8,%%zmm8	\n\t"\
			"vmovaps	0x100(%%rax),%%zmm0	\n\t	vmovaps	0x500(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x140(%%rax),%%zmm1	\n\t	vmovaps	0x540(%%rax),%%zmm6	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x100(%%rax)	\n\t	vmovaps	%%zmm5,0x500(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x140(%%rax)	\n\t	vmovaps	%%zmm7,0x540(%%rax)	\n\t"\
			/* x3.y3: */							/* xB.yB: */\
			"vmovaps	0x180(%%rbx),%%zmm2	\n\t	vmovaps	0x580(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x1c0(%%rbx),%%zmm3	\n\t	vmovaps	0x5c0(%%rbx),%%zmm8	\n\t"\
		"vsubpd	0x180(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x580(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	0x1c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x5c0(%%rcx),%%zmm8,%%zmm8	\n\t"\
			"vmovaps	0x180(%%rax),%%zmm0	\n\t	vmovaps	0x580(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x1c0(%%rax),%%zmm1	\n\t	vmovaps	0x5c0(%%rax),%%zmm6	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x180(%%rax)	\n\t	vmovaps	%%zmm5,0x580(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x1c0(%%rax)	\n\t	vmovaps	%%zmm7,0x5c0(%%rax)	\n\t"\
			/* x4.y4: */							/* xC.yC: */\
			"vmovaps	0x200(%%rbx),%%zmm2	\n\t	vmovaps	0x600(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x240(%%rbx),%%zmm3	\n\t	vmovaps	0x640(%%rbx),%%zmm8	\n\t"\
		"vsubpd	0x200(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x600(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	0x240(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x640(%%rcx),%%zmm8,%%zmm8	\n\t"\
			"vmovaps	0x200(%%rax),%%zmm0	\n\t	vmovaps	0x600(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x240(%%rax),%%zmm1	\n\t	vmovaps	0x640(%%rax),%%zmm6	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x200(%%rax)	\n\t	vmovaps	%%zmm5,0x600(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x240(%%rax)	\n\t	vmovaps	%%zmm7,0x640(%%rax)	\n\t"\
			/* x5.y5: */							/* xD.yD: */\
			"vmovaps	0x280(%%rbx),%%zmm2	\n\t	vmovaps	0x680(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x2c0(%%rbx),%%zmm3	\n\t	vmovaps	0x6c0(%%rbx),%%zmm8	\n\t"\
		"vsubpd	0x280(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x680(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	0x2c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x6c0(%%rcx),%%zmm8,%%zmm8	\n\t"\
			"vmovaps	0x280(%%rax),%%zmm0	\n\t	vmovaps	0x680(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x2c0(%%rax),%%zmm1	\n\t	vmovaps	0x6c0(%%rax),%%zmm6	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x280(%%rax)	\n\t	vmovaps	%%zmm5,0x680(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x2c0(%%rax)	\n\t	vmovaps	%%zmm7,0x6c0(%%rax)	\n\t"\
			/* x6.y6: */							/* xE.yE: */\
			"vmovaps	0x300(%%rbx),%%zmm2	\n\t	vmovaps	0x700(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x340(%%rbx),%%zmm3	\n\t	vmovaps	0x740(%%rbx),%%zmm8	\n\t"\
		"vsubpd	0x300(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x700(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	0x340(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x740(%%rcx),%%zmm8,%%zmm8	\n\t"\
			"vmovaps	0x300(%%rax),%%zmm0	\n\t	vmovaps	0x700(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x340(%%rax),%%zmm1	\n\t	vmovaps	0x740(%%rax),%%zmm6	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x300(%%rax)	\n\t	vmovaps	%%zmm5,0x700(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x340(%%rax)	\n\t	vmovaps	%%zmm7,0x740(%%rax)	\n\t"\
			/* x7.y7: */							/* xF.yF: */\
			"vmovaps	0x380(%%rbx),%%zmm2	\n\t	vmovaps	0x780(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x3c0(%%rbx),%%zmm3	\n\t	vmovaps	0x7c0(%%rbx),%%zmm8	\n\t"\
		"vsubpd	0x380(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x780(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vsubpd	0x3c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x7c0(%%rcx),%%zmm8,%%zmm8	\n\t"\
			"vmovaps	0x380(%%rax),%%zmm0	\n\t	vmovaps	0x780(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x3c0(%%rax),%%zmm1	\n\t	vmovaps	0x7c0(%%rax),%%zmm6	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x380(%%rax)	\n\t	vmovaps	%%zmm5,0x780(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x3c0(%%rax)	\n\t	vmovaps	%%zmm7,0x7c0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			 ,[__cdd0] "m" (cdd0)	// C-array base-address
			: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	// Clobbered registers
		  );

		#elif defined(USE_AVX)

		  __asm__ volatile (\
			"movq	%[__r1],%%rax		\n\t"\
			"movq	%[__bdd0],%%rbx		\n\t"\
			"movq	%[__cdd0],%%rcx		\n\t"\
			/* x0.y0: */							/* x8.y8: */\
			"vmovaps		 (%%rbx),%%ymm2	\n\t	vmovaps	0x200(%%rbx),%%ymm7	\n\t"/* y.re */\
			"vmovaps	0x020(%%rbx),%%ymm3	\n\t	vmovaps	0x220(%%rbx),%%ymm8	\n\t"/* y.im */\
		"vsubpd	     (%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x200(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x220(%%rcx),%%ymm8,%%ymm8	\n\t"\
			"vmovaps		 (%%rax),%%ymm0	\n\t	vmovaps	0x200(%%rax),%%ymm5	\n\t"/* x.re */\
			"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x220(%%rax),%%ymm6	\n\t"/* x.im */\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"/* copy x.re */\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"/* x.re *= y.re */\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"/* y.re *= x.im */\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"/* x.im *= y.im */\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"/* y.im *= x.re[copy] */\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm5,0x200(%%rax)	\n\t"/* write z.re */\
			"vmovaps	%%ymm2,0x020(%%rax)	\n\t	vmovaps	%%ymm7,0x220(%%rax)	\n\t"/* write z.im */\
			/* x1.y1: */							/* x9.y9: */\
			"vmovaps	0x040(%%rbx),%%ymm2	\n\t	vmovaps	0x240(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x060(%%rbx),%%ymm3	\n\t	vmovaps	0x260(%%rbx),%%ymm8	\n\t"\
		"vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x240(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x260(%%rcx),%%ymm8,%%ymm8	\n\t"\
			"vmovaps	0x040(%%rax),%%ymm0	\n\t	vmovaps	0x240(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x060(%%rax),%%ymm1	\n\t	vmovaps	0x260(%%rax),%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x040(%%rax)	\n\t	vmovaps	%%ymm5,0x240(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x060(%%rax)	\n\t	vmovaps	%%ymm7,0x260(%%rax)	\n\t"\
			/* x2.y2: */							/* xA.yA: */\
			"vmovaps	0x080(%%rbx),%%ymm2	\n\t	vmovaps	0x280(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x0a0(%%rbx),%%ymm3	\n\t	vmovaps	0x2a0(%%rbx),%%ymm8	\n\t"\
		"vsubpd	0x080(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x280(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	0x0a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x2a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
			"vmovaps	0x080(%%rax),%%ymm0	\n\t	vmovaps	0x280(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x0a0(%%rax),%%ymm1	\n\t	vmovaps	0x2a0(%%rax),%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x080(%%rax)	\n\t	vmovaps	%%ymm5,0x280(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x0a0(%%rax)	\n\t	vmovaps	%%ymm7,0x2a0(%%rax)	\n\t"\
			/* x3.y3: */							/* xB.yB: */\
			"vmovaps	0x0c0(%%rbx),%%ymm2	\n\t	vmovaps	0x2c0(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x0e0(%%rbx),%%ymm3	\n\t	vmovaps	0x2e0(%%rbx),%%ymm8	\n\t"\
		"vsubpd	0x0c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x2c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	0x0e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x2e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
			"vmovaps	0x0c0(%%rax),%%ymm0	\n\t	vmovaps	0x2c0(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x0e0(%%rax),%%ymm1	\n\t	vmovaps	0x2e0(%%rax),%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x0c0(%%rax)	\n\t	vmovaps	%%ymm5,0x2c0(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x0e0(%%rax)	\n\t	vmovaps	%%ymm7,0x2e0(%%rax)	\n\t"\
			/* x4.y4: */							/* xC.yC: */\
			"vmovaps	0x100(%%rbx),%%ymm2	\n\t	vmovaps	0x300(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x120(%%rbx),%%ymm3	\n\t	vmovaps	0x320(%%rbx),%%ymm8	\n\t"\
		"vsubpd	0x100(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x300(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	0x120(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x320(%%rcx),%%ymm8,%%ymm8	\n\t"\
			"vmovaps	0x100(%%rax),%%ymm0	\n\t	vmovaps	0x300(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x120(%%rax),%%ymm1	\n\t	vmovaps	0x320(%%rax),%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x100(%%rax)	\n\t	vmovaps	%%ymm5,0x300(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x120(%%rax)	\n\t	vmovaps	%%ymm7,0x320(%%rax)	\n\t"\
			/* x5.y5: */							/* xD.yD: */\
			"vmovaps	0x140(%%rbx),%%ymm2	\n\t	vmovaps	0x340(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x160(%%rbx),%%ymm3	\n\t	vmovaps	0x360(%%rbx),%%ymm8	\n\t"\
		"vsubpd	0x140(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x340(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	0x160(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x360(%%rcx),%%ymm8,%%ymm8	\n\t"\
			"vmovaps	0x140(%%rax),%%ymm0	\n\t	vmovaps	0x340(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x160(%%rax),%%ymm1	\n\t	vmovaps	0x360(%%rax),%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x140(%%rax)	\n\t	vmovaps	%%ymm5,0x340(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x160(%%rax)	\n\t	vmovaps	%%ymm7,0x360(%%rax)	\n\t"\
			/* x6.y6: */							/* xE.yE: */\
			"vmovaps	0x180(%%rbx),%%ymm2	\n\t	vmovaps	0x380(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x1a0(%%rbx),%%ymm3	\n\t	vmovaps	0x3a0(%%rbx),%%ymm8	\n\t"\
		"vsubpd	0x180(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x380(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	0x1a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x3a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
			"vmovaps	0x180(%%rax),%%ymm0	\n\t	vmovaps	0x380(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x1a0(%%rax),%%ymm1	\n\t	vmovaps	0x3a0(%%rax),%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x180(%%rax)	\n\t	vmovaps	%%ymm5,0x380(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x1a0(%%rax)	\n\t	vmovaps	%%ymm7,0x3a0(%%rax)	\n\t"\
			/* x7.y7: */							/* xF.yF: */\
			"vmovaps	0x1c0(%%rbx),%%ymm2	\n\t	vmovaps	0x3c0(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x1e0(%%rbx),%%ymm3	\n\t	vmovaps	0x3e0(%%rbx),%%ymm8	\n\t"\
		"vsubpd	0x1c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x3c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
		"vsubpd	0x1e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x3e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
			"vmovaps	0x1c0(%%rax),%%ymm0	\n\t	vmovaps	0x3c0(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x1e0(%%rax),%%ymm1	\n\t	vmovaps	0x3e0(%%rax),%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x1c0(%%rax)	\n\t	vmovaps	%%ymm5,0x3c0(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x1e0(%%rax)	\n\t	vmovaps	%%ymm7,0x3e0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			 ,[__cdd0] "m" (cdd0)	// C-array base-address
			: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	// Clobbered registers
		  );

		#else	// 64-bit SSE2:

		// Inputs are X = [x.re,x.im], Y = [y.re,y.im]
		// Output is Z = [x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re], overwriting X
		  __asm__ volatile (\
			"movq	%[__r1],%%rax		\n\t"\
			"movq	%[__bdd0],%%rbx		\n\t"\
			"movq	%[__cdd0],%%rcx		\n\t"\
			/* x0.y0: */						/* x8.y8: */\
			"movaps		 (%%rbx),%%xmm2	\n\t	movaps	0x100(%%rbx),%%xmm7	\n\t"/* y.re */\
			"movaps	0x010(%%rbx),%%xmm3	\n\t	movaps	0x110(%%rbx),%%xmm8	\n\t"/* y.im */\
			"subpd	     (%%rcx),%%xmm2	\n\t	subpd	0x100(%%rcx),%%xmm7	\n\t"/* y.re -= z.re */\
			"subpd	0x010(%%rcx),%%xmm3	\n\t	subpd	0x110(%%rcx),%%xmm8	\n\t"/* y.im -= z.im */\
			"movaps		 (%%rax),%%xmm0	\n\t	movaps	0x100(%%rax),%%xmm5	\n\t"/* x.re */\
			"movaps	0x010(%%rax),%%xmm1	\n\t	movaps	0x110(%%rax),%%xmm6	\n\t"/* x.im */\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"/* copy x.re */\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"/* x.re *= y.re */\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"/* y.re *= x.im */\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"/* x.im *= y.im */\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"/* y.im *= x.re[copy] */\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm5,0x100(%%rax)	\n\t"/* write z.re */\
			"movaps	%%xmm2,0x010(%%rax)	\n\t	movaps	%%xmm7,0x110(%%rax)	\n\t"/* write z.im */\
			/* x1.y1: */						/* x9.y9: */\
			"movaps	0x020(%%rbx),%%xmm2	\n\t	movaps	0x120(%%rbx),%%xmm7	\n\t"\
			"movaps	0x030(%%rbx),%%xmm3	\n\t	movaps	0x130(%%rbx),%%xmm8	\n\t"\
			"subpd	0x020(%%rcx),%%xmm2	\n\t	subpd	0x120(%%rcx),%%xmm7	\n\t"\
			"subpd	0x030(%%rcx),%%xmm3	\n\t	subpd	0x130(%%rcx),%%xmm8	\n\t"\
			"movaps	0x020(%%rax),%%xmm0	\n\t	movaps	0x120(%%rax),%%xmm5	\n\t"\
			"movaps	0x030(%%rax),%%xmm1	\n\t	movaps	0x130(%%rax),%%xmm6	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x020(%%rax)	\n\t	movaps	%%xmm5,0x120(%%rax)	\n\t"\
			"movaps	%%xmm2,0x030(%%rax)	\n\t	movaps	%%xmm7,0x130(%%rax)	\n\t"\
			/* x2.y2: */						/* xA.yA: */\
			"movaps	0x040(%%rbx),%%xmm2	\n\t	movaps	0x140(%%rbx),%%xmm7	\n\t"\
			"movaps	0x050(%%rbx),%%xmm3	\n\t	movaps	0x150(%%rbx),%%xmm8	\n\t"\
			"subpd	0x040(%%rcx),%%xmm2	\n\t	subpd	0x140(%%rcx),%%xmm7	\n\t"\
			"subpd	0x050(%%rcx),%%xmm3	\n\t	subpd	0x150(%%rcx),%%xmm8	\n\t"\
			"movaps	0x040(%%rax),%%xmm0	\n\t	movaps	0x140(%%rax),%%xmm5	\n\t"\
			"movaps	0x050(%%rax),%%xmm1	\n\t	movaps	0x150(%%rax),%%xmm6	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x040(%%rax)	\n\t	movaps	%%xmm5,0x140(%%rax)	\n\t"\
			"movaps	%%xmm2,0x050(%%rax)	\n\t	movaps	%%xmm7,0x150(%%rax)	\n\t"\
			/* x3.y3: */						/* xB.yB: */\
			"movaps	0x060(%%rbx),%%xmm2	\n\t	movaps	0x160(%%rbx),%%xmm7	\n\t"\
			"movaps	0x070(%%rbx),%%xmm3	\n\t	movaps	0x170(%%rbx),%%xmm8	\n\t"\
			"subpd	0x060(%%rcx),%%xmm2	\n\t	subpd	0x160(%%rcx),%%xmm7	\n\t"\
			"subpd	0x070(%%rcx),%%xmm3	\n\t	subpd	0x170(%%rcx),%%xmm8	\n\t"\
			"movaps	0x060(%%rax),%%xmm0	\n\t	movaps	0x160(%%rax),%%xmm5	\n\t"\
			"movaps	0x070(%%rax),%%xmm1	\n\t	movaps	0x170(%%rax),%%xmm6	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x060(%%rax)	\n\t	movaps	%%xmm5,0x160(%%rax)	\n\t"\
			"movaps	%%xmm2,0x070(%%rax)	\n\t	movaps	%%xmm7,0x170(%%rax)	\n\t"\
			/* x4.y4: */						/* xC.yC: */\
			"movaps	0x080(%%rbx),%%xmm2	\n\t	movaps	0x180(%%rbx),%%xmm7	\n\t"\
			"movaps	0x090(%%rbx),%%xmm3	\n\t	movaps	0x190(%%rbx),%%xmm8	\n\t"\
			"subpd	0x080(%%rcx),%%xmm2	\n\t	subpd	0x180(%%rcx),%%xmm7	\n\t"\
			"subpd	0x090(%%rcx),%%xmm3	\n\t	subpd	0x190(%%rcx),%%xmm8	\n\t"\
			"movaps	0x080(%%rax),%%xmm0	\n\t	movaps	0x180(%%rax),%%xmm5	\n\t"\
			"movaps	0x090(%%rax),%%xmm1	\n\t	movaps	0x190(%%rax),%%xmm6	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x080(%%rax)	\n\t	movaps	%%xmm5,0x180(%%rax)	\n\t"\
			"movaps	%%xmm2,0x090(%%rax)	\n\t	movaps	%%xmm7,0x190(%%rax)	\n\t"\
			/* x5.y5: */						/* xD.yD: */\
			"movaps	0x0a0(%%rbx),%%xmm2	\n\t	movaps	0x1a0(%%rbx),%%xmm7	\n\t"\
			"movaps	0x0b0(%%rbx),%%xmm3	\n\t	movaps	0x1b0(%%rbx),%%xmm8	\n\t"\
			"subpd	0x0a0(%%rcx),%%xmm2	\n\t	subpd	0x1a0(%%rcx),%%xmm7	\n\t"\
			"subpd	0x0b0(%%rcx),%%xmm3	\n\t	subpd	0x1b0(%%rcx),%%xmm8	\n\t"\
			"movaps	0x0a0(%%rax),%%xmm0	\n\t	movaps	0x1a0(%%rax),%%xmm5	\n\t"\
			"movaps	0x0b0(%%rax),%%xmm1	\n\t	movaps	0x1b0(%%rax),%%xmm6	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x0a0(%%rax)	\n\t	movaps	%%xmm5,0x1a0(%%rax)	\n\t"\
			"movaps	%%xmm2,0x0b0(%%rax)	\n\t	movaps	%%xmm7,0x1b0(%%rax)	\n\t"\
			/* x6.y6: */						/* xE.yE: */\
			"movaps	0x0c0(%%rbx),%%xmm2	\n\t	movaps	0x1c0(%%rbx),%%xmm7	\n\t"\
			"movaps	0x0d0(%%rbx),%%xmm3	\n\t	movaps	0x1d0(%%rbx),%%xmm8	\n\t"\
			"subpd	0x0c0(%%rcx),%%xmm2	\n\t	subpd	0x1c0(%%rcx),%%xmm7	\n\t"\
			"subpd	0x0d0(%%rcx),%%xmm3	\n\t	subpd	0x1d0(%%rcx),%%xmm8	\n\t"\
			"movaps	0x0c0(%%rax),%%xmm0	\n\t	movaps	0x1c0(%%rax),%%xmm5	\n\t"\
			"movaps	0x0d0(%%rax),%%xmm1	\n\t	movaps	0x1d0(%%rax),%%xmm6	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x0c0(%%rax)	\n\t	movaps	%%xmm5,0x1c0(%%rax)	\n\t"\
			"movaps	%%xmm2,0x0d0(%%rax)	\n\t	movaps	%%xmm7,0x1d0(%%rax)	\n\t"\
			/* x7.y7: */						/* xF.yF: */\
			"movaps	0x0e0(%%rbx),%%xmm2	\n\t	movaps	0x1e0(%%rbx),%%xmm7	\n\t"\
			"movaps	0x0f0(%%rbx),%%xmm3	\n\t	movaps	0x1f0(%%rbx),%%xmm8	\n\t"\
			"subpd	0x0e0(%%rcx),%%xmm2	\n\t	subpd	0x1e0(%%rcx),%%xmm7	\n\t"\
			"subpd	0x0f0(%%rcx),%%xmm3	\n\t	subpd	0x1f0(%%rcx),%%xmm8	\n\t"\
			"movaps	0x0e0(%%rax),%%xmm0	\n\t	movaps	0x1e0(%%rax),%%xmm5	\n\t"\
			"movaps	0x0f0(%%rax),%%xmm1	\n\t	movaps	0x1f0(%%rax),%%xmm6	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x0e0(%%rax)	\n\t	movaps	%%xmm5,0x1e0(%%rax)	\n\t"\
			"movaps	%%xmm2,0x0f0(%%rax)	\n\t	movaps	%%xmm7,0x1f0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			 ,[__cdd0] "m" (cdd0)	// C-array base-address
			: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	// Clobbered registers
		  );

		#endif	// AVX or SSE2?

		} else {	// c_arr = 0x0: a * b

		#ifdef USE_ARM_V8_SIMD

		// No Fermat-mod support on ARMv8, just supply a stub macro:
		  __asm__ volatile (\
			"ldr x0,%[__r1]	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			: "cc","memory","x0"	/* Clobbered registers */\
		  );

		#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

		  __asm__ volatile (\
			"movq	%[__r1],%%rax		\n\t"\
			"movq	%[__bdd0],%%rbx		\n\t"\
			/* x0.y0: */							/* x8.y8: */\
			"vmovaps		 (%%rax),%%zmm0	\n\t	vmovaps	0x400(%%rax),%%zmm5	\n\t"/* x.re */\
			"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x440(%%rax),%%zmm6	\n\t"/* x.im */\
			"vmovaps		 (%%rbx),%%zmm2	\n\t	vmovaps	0x400(%%rbx),%%zmm7	\n\t"/* y.re */\
			"vmovaps	0x040(%%rbx),%%zmm3	\n\t	vmovaps	0x440(%%rbx),%%zmm8	\n\t"/* y.im */\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"/* copy x.re */\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"/* x.re *= y.re */\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"/* y.re *= x.im */\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"/* x.im *= y.im */\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"/* y.im *= x.re[copy] */\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm5,0x400(%%rax)	\n\t"/* write z.re */\
			"vmovaps	%%zmm2,0x040(%%rax)	\n\t	vmovaps	%%zmm7,0x440(%%rax)	\n\t"/* write z.im */\
			/* x1.y1: */							/* x9.y9: */\
			"vmovaps	0x080(%%rax),%%zmm0	\n\t	vmovaps	0x480(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x0c0(%%rax),%%zmm1	\n\t	vmovaps	0x4c0(%%rax),%%zmm6	\n\t"\
			"vmovaps	0x080(%%rbx),%%zmm2	\n\t	vmovaps	0x480(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x0c0(%%rbx),%%zmm3	\n\t	vmovaps	0x4c0(%%rbx),%%zmm8	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x080(%%rax)	\n\t	vmovaps	%%zmm5,0x480(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x0c0(%%rax)	\n\t	vmovaps	%%zmm7,0x4c0(%%rax)	\n\t"\
			/* x2.y2: */							/* xA.yA: */\
			"vmovaps	0x100(%%rax),%%zmm0	\n\t	vmovaps	0x500(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x140(%%rax),%%zmm1	\n\t	vmovaps	0x540(%%rax),%%zmm6	\n\t"\
			"vmovaps	0x100(%%rbx),%%zmm2	\n\t	vmovaps	0x500(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x140(%%rbx),%%zmm3	\n\t	vmovaps	0x540(%%rbx),%%zmm8	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x100(%%rax)	\n\t	vmovaps	%%zmm5,0x500(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x140(%%rax)	\n\t	vmovaps	%%zmm7,0x540(%%rax)	\n\t"\
			/* x3.y3: */							/* xB.yB: */\
			"vmovaps	0x180(%%rax),%%zmm0	\n\t	vmovaps	0x580(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x1c0(%%rax),%%zmm1	\n\t	vmovaps	0x5c0(%%rax),%%zmm6	\n\t"\
			"vmovaps	0x180(%%rbx),%%zmm2	\n\t	vmovaps	0x580(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x1c0(%%rbx),%%zmm3	\n\t	vmovaps	0x5c0(%%rbx),%%zmm8	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x180(%%rax)	\n\t	vmovaps	%%zmm5,0x580(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x1c0(%%rax)	\n\t	vmovaps	%%zmm7,0x5c0(%%rax)	\n\t"\
			/* x4.y4: */							/* xC.yC: */\
			"vmovaps	0x200(%%rax),%%zmm0	\n\t	vmovaps	0x600(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x240(%%rax),%%zmm1	\n\t	vmovaps	0x640(%%rax),%%zmm6	\n\t"\
			"vmovaps	0x200(%%rbx),%%zmm2	\n\t	vmovaps	0x600(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x240(%%rbx),%%zmm3	\n\t	vmovaps	0x640(%%rbx),%%zmm8	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x200(%%rax)	\n\t	vmovaps	%%zmm5,0x600(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x240(%%rax)	\n\t	vmovaps	%%zmm7,0x640(%%rax)	\n\t"\
			/* x5.y5: */							/* xD.yD: */\
			"vmovaps	0x280(%%rax),%%zmm0	\n\t	vmovaps	0x680(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x2c0(%%rax),%%zmm1	\n\t	vmovaps	0x6c0(%%rax),%%zmm6	\n\t"\
			"vmovaps	0x280(%%rbx),%%zmm2	\n\t	vmovaps	0x680(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x2c0(%%rbx),%%zmm3	\n\t	vmovaps	0x6c0(%%rbx),%%zmm8	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x280(%%rax)	\n\t	vmovaps	%%zmm5,0x680(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x2c0(%%rax)	\n\t	vmovaps	%%zmm7,0x6c0(%%rax)	\n\t"\
			/* x6.y6: */							/* xE.yE: */\
			"vmovaps	0x300(%%rax),%%zmm0	\n\t	vmovaps	0x700(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x340(%%rax),%%zmm1	\n\t	vmovaps	0x740(%%rax),%%zmm6	\n\t"\
			"vmovaps	0x300(%%rbx),%%zmm2	\n\t	vmovaps	0x700(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x340(%%rbx),%%zmm3	\n\t	vmovaps	0x740(%%rbx),%%zmm8	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x300(%%rax)	\n\t	vmovaps	%%zmm5,0x700(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x340(%%rax)	\n\t	vmovaps	%%zmm7,0x740(%%rax)	\n\t"\
			/* x7.y7: */							/* xF.yF: */\
			"vmovaps	0x380(%%rax),%%zmm0	\n\t	vmovaps	0x780(%%rax),%%zmm5	\n\t"\
			"vmovaps	0x3c0(%%rax),%%zmm1	\n\t	vmovaps	0x7c0(%%rax),%%zmm6	\n\t"\
			"vmovaps	0x380(%%rbx),%%zmm2	\n\t	vmovaps	0x780(%%rbx),%%zmm7	\n\t"\
			"vmovaps	0x3c0(%%rbx),%%zmm3	\n\t	vmovaps	0x7c0(%%rbx),%%zmm8	\n\t"\
			"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
			"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
			"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
			"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
			"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
			"vmovaps	%%zmm0,0x380(%%rax)	\n\t	vmovaps	%%zmm5,0x780(%%rax)	\n\t"\
			"vmovaps	%%zmm2,0x3c0(%%rax)	\n\t	vmovaps	%%zmm7,0x7c0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	// Clobbered registers
		  );

		#elif defined(USE_AVX)

		  __asm__ volatile (\
			"movq	%[__r1],%%rax		\n\t"\
			"movq	%[__bdd0],%%rbx		\n\t"\
			/* x0.y0: */							/* x8.y8: */\
			"vmovaps		 (%%rax),%%ymm0	\n\t	vmovaps	0x200(%%rax),%%ymm5	\n\t"/* x.re */\
			"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x220(%%rax),%%ymm6	\n\t"/* x.im */\
			"vmovaps		 (%%rbx),%%ymm2	\n\t	vmovaps	0x200(%%rbx),%%ymm7	\n\t"/* y.re */\
			"vmovaps	0x020(%%rbx),%%ymm3	\n\t	vmovaps	0x220(%%rbx),%%ymm8	\n\t"/* y.im */\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"/* copy x.re */\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"/* x.re *= y.re */\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"/* y.re *= x.im */\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"/* x.im *= y.im */\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"/* y.im *= x.re[copy] */\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm5,0x200(%%rax)	\n\t"/* write z.re */\
			"vmovaps	%%ymm2,0x020(%%rax)	\n\t	vmovaps	%%ymm7,0x220(%%rax)	\n\t"/* write z.im */\
			/* x1.y1: */							/* x9.y9: */\
			"vmovaps	0x040(%%rax),%%ymm0	\n\t	vmovaps	0x240(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x060(%%rax),%%ymm1	\n\t	vmovaps	0x260(%%rax),%%ymm6	\n\t"\
			"vmovaps	0x040(%%rbx),%%ymm2	\n\t	vmovaps	0x240(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x060(%%rbx),%%ymm3	\n\t	vmovaps	0x260(%%rbx),%%ymm8	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x040(%%rax)	\n\t	vmovaps	%%ymm5,0x240(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x060(%%rax)	\n\t	vmovaps	%%ymm7,0x260(%%rax)	\n\t"\
			/* x2.y2: */							/* xA.yA: */\
			"vmovaps	0x080(%%rax),%%ymm0	\n\t	vmovaps	0x280(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x0a0(%%rax),%%ymm1	\n\t	vmovaps	0x2a0(%%rax),%%ymm6	\n\t"\
			"vmovaps	0x080(%%rbx),%%ymm2	\n\t	vmovaps	0x280(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x0a0(%%rbx),%%ymm3	\n\t	vmovaps	0x2a0(%%rbx),%%ymm8	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x080(%%rax)	\n\t	vmovaps	%%ymm5,0x280(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x0a0(%%rax)	\n\t	vmovaps	%%ymm7,0x2a0(%%rax)	\n\t"\
			/* x3.y3: */							/* xB.yB: */\
			"vmovaps	0x0c0(%%rax),%%ymm0	\n\t	vmovaps	0x2c0(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x0e0(%%rax),%%ymm1	\n\t	vmovaps	0x2e0(%%rax),%%ymm6	\n\t"\
			"vmovaps	0x0c0(%%rbx),%%ymm2	\n\t	vmovaps	0x2c0(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x0e0(%%rbx),%%ymm3	\n\t	vmovaps	0x2e0(%%rbx),%%ymm8	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x0c0(%%rax)	\n\t	vmovaps	%%ymm5,0x2c0(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x0e0(%%rax)	\n\t	vmovaps	%%ymm7,0x2e0(%%rax)	\n\t"\
			/* x4.y4: */							/* xC.yC: */\
			"vmovaps	0x100(%%rax),%%ymm0	\n\t	vmovaps	0x300(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x120(%%rax),%%ymm1	\n\t	vmovaps	0x320(%%rax),%%ymm6	\n\t"\
			"vmovaps	0x100(%%rbx),%%ymm2	\n\t	vmovaps	0x300(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x120(%%rbx),%%ymm3	\n\t	vmovaps	0x320(%%rbx),%%ymm8	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x100(%%rax)	\n\t	vmovaps	%%ymm5,0x300(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x120(%%rax)	\n\t	vmovaps	%%ymm7,0x320(%%rax)	\n\t"\
			/* x5.y5: */							/* xD.yD: */\
			"vmovaps	0x140(%%rax),%%ymm0	\n\t	vmovaps	0x340(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x160(%%rax),%%ymm1	\n\t	vmovaps	0x360(%%rax),%%ymm6	\n\t"\
			"vmovaps	0x140(%%rbx),%%ymm2	\n\t	vmovaps	0x340(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x160(%%rbx),%%ymm3	\n\t	vmovaps	0x360(%%rbx),%%ymm8	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x140(%%rax)	\n\t	vmovaps	%%ymm5,0x340(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x160(%%rax)	\n\t	vmovaps	%%ymm7,0x360(%%rax)	\n\t"\
			/* x6.y6: */							/* xE.yE: */\
			"vmovaps	0x180(%%rax),%%ymm0	\n\t	vmovaps	0x380(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x1a0(%%rax),%%ymm1	\n\t	vmovaps	0x3a0(%%rax),%%ymm6	\n\t"\
			"vmovaps	0x180(%%rbx),%%ymm2	\n\t	vmovaps	0x380(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x1a0(%%rbx),%%ymm3	\n\t	vmovaps	0x3a0(%%rbx),%%ymm8	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x180(%%rax)	\n\t	vmovaps	%%ymm5,0x380(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x1a0(%%rax)	\n\t	vmovaps	%%ymm7,0x3a0(%%rax)	\n\t"\
			/* x7.y7: */							/* xF.yF: */\
			"vmovaps	0x1c0(%%rax),%%ymm0	\n\t	vmovaps	0x3c0(%%rax),%%ymm5	\n\t"\
			"vmovaps	0x1e0(%%rax),%%ymm1	\n\t	vmovaps	0x3e0(%%rax),%%ymm6	\n\t"\
			"vmovaps	0x1c0(%%rbx),%%ymm2	\n\t	vmovaps	0x3c0(%%rbx),%%ymm7	\n\t"\
			"vmovaps	0x1e0(%%rbx),%%ymm3	\n\t	vmovaps	0x3e0(%%rbx),%%ymm8	\n\t"\
			"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
			"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
			"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
			"vmovaps	%%ymm0,0x1c0(%%rax)	\n\t	vmovaps	%%ymm5,0x3c0(%%rax)	\n\t"\
			"vmovaps	%%ymm2,0x1e0(%%rax)	\n\t	vmovaps	%%ymm7,0x3e0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	// Clobbered registers
		  );

		#else	// 64-bit SSE2:

		// Inputs are X = [x.re,x.im], Y = [y.re,y.im]
		// Output is Z = [x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re], overwriting X
		  __asm__ volatile (\
			"movq	%[__r1],%%rax		\n\t"\
			"movq	%[__bdd0],%%rbx		\n\t"\
			/* x0.y0: */						/* x8.y8: */\
			"movaps		 (%%rax),%%xmm0	\n\t	movaps	0x100(%%rax),%%xmm5	\n\t"/* x.re */\
			"movaps	0x010(%%rax),%%xmm1	\n\t	movaps	0x110(%%rax),%%xmm6	\n\t"/* x.im */\
			"movaps		 (%%rbx),%%xmm2	\n\t	movaps	0x100(%%rbx),%%xmm7	\n\t"/* y.re */\
			"movaps	0x010(%%rbx),%%xmm3	\n\t	movaps	0x110(%%rbx),%%xmm8	\n\t"/* y.im */\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"/* copy x.re */\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"/* x.re *= y.re */\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"/* y.re *= x.im */\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"/* x.im *= y.im */\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"/* y.im *= x.re[copy] */\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm5,0x100(%%rax)	\n\t"/* write z.re */\
			"movaps	%%xmm2,0x010(%%rax)	\n\t	movaps	%%xmm7,0x110(%%rax)	\n\t"/* write z.im */\
			/* x1.y1: */						/* x9.y9: */\
			"movaps	0x020(%%rax),%%xmm0	\n\t	movaps	0x120(%%rax),%%xmm5	\n\t"\
			"movaps	0x030(%%rax),%%xmm1	\n\t	movaps	0x130(%%rax),%%xmm6	\n\t"\
			"movaps	0x020(%%rbx),%%xmm2	\n\t	movaps	0x120(%%rbx),%%xmm7	\n\t"\
			"movaps	0x030(%%rbx),%%xmm3	\n\t	movaps	0x130(%%rbx),%%xmm8	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x020(%%rax)	\n\t	movaps	%%xmm5,0x120(%%rax)	\n\t"\
			"movaps	%%xmm2,0x030(%%rax)	\n\t	movaps	%%xmm7,0x130(%%rax)	\n\t"\
			/* x2.y2: */						/* xA.yA: */\
			"movaps	0x040(%%rax),%%xmm0	\n\t	movaps	0x140(%%rax),%%xmm5	\n\t"\
			"movaps	0x050(%%rax),%%xmm1	\n\t	movaps	0x150(%%rax),%%xmm6	\n\t"\
			"movaps	0x040(%%rbx),%%xmm2	\n\t	movaps	0x140(%%rbx),%%xmm7	\n\t"\
			"movaps	0x050(%%rbx),%%xmm3	\n\t	movaps	0x150(%%rbx),%%xmm8	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x040(%%rax)	\n\t	movaps	%%xmm5,0x140(%%rax)	\n\t"\
			"movaps	%%xmm2,0x050(%%rax)	\n\t	movaps	%%xmm7,0x150(%%rax)	\n\t"\
			/* x3.y3: */						/* xB.yB: */\
			"movaps	0x060(%%rax),%%xmm0	\n\t	movaps	0x160(%%rax),%%xmm5	\n\t"\
			"movaps	0x070(%%rax),%%xmm1	\n\t	movaps	0x170(%%rax),%%xmm6	\n\t"\
			"movaps	0x060(%%rbx),%%xmm2	\n\t	movaps	0x160(%%rbx),%%xmm7	\n\t"\
			"movaps	0x070(%%rbx),%%xmm3	\n\t	movaps	0x170(%%rbx),%%xmm8	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x060(%%rax)	\n\t	movaps	%%xmm5,0x160(%%rax)	\n\t"\
			"movaps	%%xmm2,0x070(%%rax)	\n\t	movaps	%%xmm7,0x170(%%rax)	\n\t"\
			/* x4.y4: */						/* xC.yC: */\
			"movaps	0x080(%%rax),%%xmm0	\n\t	movaps	0x180(%%rax),%%xmm5	\n\t"\
			"movaps	0x090(%%rax),%%xmm1	\n\t	movaps	0x190(%%rax),%%xmm6	\n\t"\
			"movaps	0x080(%%rbx),%%xmm2	\n\t	movaps	0x180(%%rbx),%%xmm7	\n\t"\
			"movaps	0x090(%%rbx),%%xmm3	\n\t	movaps	0x190(%%rbx),%%xmm8	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x080(%%rax)	\n\t	movaps	%%xmm5,0x180(%%rax)	\n\t"\
			"movaps	%%xmm2,0x090(%%rax)	\n\t	movaps	%%xmm7,0x190(%%rax)	\n\t"\
			/* x5.y5: */						/* xD.yD: */\
			"movaps	0x0a0(%%rax),%%xmm0	\n\t	movaps	0x1a0(%%rax),%%xmm5	\n\t"\
			"movaps	0x0b0(%%rax),%%xmm1	\n\t	movaps	0x1b0(%%rax),%%xmm6	\n\t"\
			"movaps	0x0a0(%%rbx),%%xmm2	\n\t	movaps	0x1a0(%%rbx),%%xmm7	\n\t"\
			"movaps	0x0b0(%%rbx),%%xmm3	\n\t	movaps	0x1b0(%%rbx),%%xmm8	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x0a0(%%rax)	\n\t	movaps	%%xmm5,0x1a0(%%rax)	\n\t"\
			"movaps	%%xmm2,0x0b0(%%rax)	\n\t	movaps	%%xmm7,0x1b0(%%rax)	\n\t"\
			/* x6.y6: */						/* xE.yE: */\
			"movaps	0x0c0(%%rax),%%xmm0	\n\t	movaps	0x1c0(%%rax),%%xmm5	\n\t"\
			"movaps	0x0d0(%%rax),%%xmm1	\n\t	movaps	0x1d0(%%rax),%%xmm6	\n\t"\
			"movaps	0x0c0(%%rbx),%%xmm2	\n\t	movaps	0x1c0(%%rbx),%%xmm7	\n\t"\
			"movaps	0x0d0(%%rbx),%%xmm3	\n\t	movaps	0x1d0(%%rbx),%%xmm8	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x0c0(%%rax)	\n\t	movaps	%%xmm5,0x1c0(%%rax)	\n\t"\
			"movaps	%%xmm2,0x0d0(%%rax)	\n\t	movaps	%%xmm7,0x1d0(%%rax)	\n\t"\
			/* x7.y7: */						/* xF.yF: */\
			"movaps	0x0e0(%%rax),%%xmm0	\n\t	movaps	0x1e0(%%rax),%%xmm5	\n\t"\
			"movaps	0x0f0(%%rax),%%xmm1	\n\t	movaps	0x1f0(%%rax),%%xmm6	\n\t"\
			"movaps	0x0e0(%%rbx),%%xmm2	\n\t	movaps	0x1e0(%%rbx),%%xmm7	\n\t"\
			"movaps	0x0f0(%%rbx),%%xmm3	\n\t	movaps	0x1f0(%%rbx),%%xmm8	\n\t"\
			"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
			"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
			"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
			"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
			"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
			"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
			"movaps	%%xmm0,0x0e0(%%rax)	\n\t	movaps	%%xmm5,0x1e0(%%rax)	\n\t"\
			"movaps	%%xmm2,0x0f0(%%rax)	\n\t	movaps	%%xmm7,0x1f0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	// Clobbered registers
		  );

		#endif	// AVX or SSE2?

		}	// endif(c_arr)

	} else {	// fwd_fft_only = 0: Normal execution, dyadic-squaring followed by iFFT:

	  #ifdef USE_ARM_V8_SIMD
		// No Fermat-mod support on ARMv8, just supply a stub macro:
		__asm__ volatile (\
			"ldr x0,%[__r1]	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			: "cc","memory","x0"	/* Clobbered registers */\
		);
	  #elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro
		__asm__ volatile (\
			"movq	%[__r1],%%rax			\n\t"\
			/* z0^2: */								/* z8^2: */							/* z4^2: */								/* zC^2: */\
			"vmovaps		 (%%rax),%%zmm3	\n\t	vmovaps	0x400(%%rax),%%zmm7	\n\t	vmovaps	0x200(%%rax),%%zmm11	\n\t	vmovaps	0x600(%%rax),%%zmm15	\n\t"\
			"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x440(%%rax),%%zmm5	\n\t	vmovaps	0x240(%%rax),%%zmm9		\n\t	vmovaps	0x640(%%rax),%%zmm13	\n\t"\
			"vsubpd	%%zmm1,%%zmm3,%%zmm2	\n\t	vsubpd	%%zmm5,%%zmm7,%%zmm6\n\t	vsubpd	%%zmm9,%%zmm11,%%zmm10	\n\t	vsubpd	%%zmm13,%%zmm15,%%zmm14	\n\t"\
			"vaddpd	%%zmm1,%%zmm3,%%zmm0	\n\t	vaddpd	%%zmm5,%%zmm7,%%zmm4\n\t	vaddpd	%%zmm9,%%zmm11,%%zmm8	\n\t	vaddpd	%%zmm13,%%zmm15,%%zmm12	\n\t"\
			"vaddpd	%%zmm1,%%zmm1,%%zmm1	\n\t	vaddpd	%%zmm5,%%zmm5,%%zmm5\n\t	vaddpd	%%zmm9,%%zmm9,%%zmm9	\n\t	vaddpd	%%zmm13,%%zmm13,%%zmm13	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t	vmulpd	%%zmm10,%%zmm8,%%zmm8	\n\t	vmulpd	%%zmm14,%%zmm12,%%zmm12	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t	vmulpd	%%zmm11,%%zmm9,%%zmm9	\n\t	vmulpd	%%zmm15,%%zmm13,%%zmm13	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm4,0x400(%%rax)	\n\t	vmovaps	%%zmm8,0x200(%%rax)		\n\t	vmovaps	%%zmm12,0x600(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm5,0x440(%%rax)	\n\t	vmovaps	%%zmm9,0x240(%%rax)		\n\t	vmovaps	%%zmm13,0x640(%%rax)	\n\t"\
			/* z1^2: */								/* z9^2: */							/* z5^2: */								/* zD^2: */\
			"vmovaps	0x080(%%rax),%%zmm3	\n\t	vmovaps	0x480(%%rax),%%zmm7	\n\t	vmovaps	0x280(%%rax),%%zmm11	\n\t	vmovaps	0x680(%%rax),%%zmm15	\n\t"\
			"vmovaps	0x0c0(%%rax),%%zmm1	\n\t	vmovaps	0x4c0(%%rax),%%zmm5	\n\t	vmovaps	0x2c0(%%rax),%%zmm9		\n\t	vmovaps	0x6c0(%%rax),%%zmm13	\n\t"\
			"vsubpd	%%zmm1,%%zmm3,%%zmm2	\n\t	vsubpd	%%zmm5,%%zmm7,%%zmm6\n\t	vsubpd	%%zmm9,%%zmm11,%%zmm10	\n\t	vsubpd	%%zmm13,%%zmm15,%%zmm14	\n\t"\
			"vaddpd	%%zmm1,%%zmm3,%%zmm0	\n\t	vaddpd	%%zmm5,%%zmm7,%%zmm4\n\t	vaddpd	%%zmm9,%%zmm11,%%zmm8	\n\t	vaddpd	%%zmm13,%%zmm15,%%zmm12	\n\t"\
			"vaddpd	%%zmm1,%%zmm1,%%zmm1	\n\t	vaddpd	%%zmm5,%%zmm5,%%zmm5\n\t	vaddpd	%%zmm9,%%zmm9,%%zmm9	\n\t	vaddpd	%%zmm13,%%zmm13,%%zmm13	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t	vmulpd	%%zmm10,%%zmm8,%%zmm8	\n\t	vmulpd	%%zmm14,%%zmm12,%%zmm12	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t	vmulpd	%%zmm11,%%zmm9,%%zmm9	\n\t	vmulpd	%%zmm15,%%zmm13,%%zmm13	\n\t"\
			"vmovaps	%%zmm0,0x080(%%rax)	\n\t	vmovaps	%%zmm4,0x480(%%rax)	\n\t	vmovaps	%%zmm8,0x280(%%rax)		\n\t	vmovaps	%%zmm12,0x680(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x0c0(%%rax)	\n\t	vmovaps	%%zmm5,0x4c0(%%rax)	\n\t	vmovaps	%%zmm9,0x2c0(%%rax)		\n\t	vmovaps	%%zmm13,0x6c0(%%rax)	\n\t"\
			/* z2^2: */								/* zA^2: */							/* z6^2: */								/* zE^2: */\
			"vmovaps	0x100(%%rax),%%zmm3	\n\t	vmovaps	0x500(%%rax),%%zmm7	\n\t	vmovaps	0x300(%%rax),%%zmm11	\n\t	vmovaps	0x700(%%rax),%%zmm15	\n\t"\
			"vmovaps	0x140(%%rax),%%zmm1	\n\t	vmovaps	0x540(%%rax),%%zmm5	\n\t	vmovaps	0x340(%%rax),%%zmm9		\n\t	vmovaps	0x740(%%rax),%%zmm13	\n\t"\
			"vsubpd	%%zmm1,%%zmm3,%%zmm2	\n\t	vsubpd	%%zmm5,%%zmm7,%%zmm6\n\t	vsubpd	%%zmm9,%%zmm11,%%zmm10	\n\t	vsubpd	%%zmm13,%%zmm15,%%zmm14	\n\t"\
			"vaddpd	%%zmm1,%%zmm3,%%zmm0	\n\t	vaddpd	%%zmm5,%%zmm7,%%zmm4\n\t	vaddpd	%%zmm9,%%zmm11,%%zmm8	\n\t	vaddpd	%%zmm13,%%zmm15,%%zmm12	\n\t"\
			"vaddpd	%%zmm1,%%zmm1,%%zmm1	\n\t	vaddpd	%%zmm5,%%zmm5,%%zmm5\n\t	vaddpd	%%zmm9,%%zmm9,%%zmm9	\n\t	vaddpd	%%zmm13,%%zmm13,%%zmm13	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t	vmulpd	%%zmm10,%%zmm8,%%zmm8	\n\t	vmulpd	%%zmm14,%%zmm12,%%zmm12	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t	vmulpd	%%zmm11,%%zmm9,%%zmm9	\n\t	vmulpd	%%zmm15,%%zmm13,%%zmm13	\n\t"\
			"vmovaps	%%zmm0,0x100(%%rax)	\n\t	vmovaps	%%zmm4,0x500(%%rax)	\n\t	vmovaps	%%zmm8,0x300(%%rax)		\n\t	vmovaps	%%zmm12,0x700(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x140(%%rax)	\n\t	vmovaps	%%zmm5,0x540(%%rax)	\n\t	vmovaps	%%zmm9,0x340(%%rax)		\n\t	vmovaps	%%zmm13,0x740(%%rax)	\n\t"\
			/* z3^2: */								/* zB^2: */							/* z7^2: */								/* zF^2: */\
			"vmovaps	0x180(%%rax),%%zmm3	\n\t	vmovaps	0x580(%%rax),%%zmm7	\n\t	vmovaps	0x380(%%rax),%%zmm11	\n\t	vmovaps	0x780(%%rax),%%zmm15	\n\t"\
			"vmovaps	0x1c0(%%rax),%%zmm1	\n\t	vmovaps	0x5c0(%%rax),%%zmm5	\n\t	vmovaps	0x3c0(%%rax),%%zmm9		\n\t	vmovaps	0x7c0(%%rax),%%zmm13	\n\t"\
			"vsubpd	%%zmm1,%%zmm3,%%zmm2	\n\t	vsubpd	%%zmm5,%%zmm7,%%zmm6\n\t	vsubpd	%%zmm9,%%zmm11,%%zmm10	\n\t	vsubpd	%%zmm13,%%zmm15,%%zmm14	\n\t"\
			"vaddpd	%%zmm1,%%zmm3,%%zmm0	\n\t	vaddpd	%%zmm5,%%zmm7,%%zmm4\n\t	vaddpd	%%zmm9,%%zmm11,%%zmm8	\n\t	vaddpd	%%zmm13,%%zmm15,%%zmm12	\n\t"\
			"vaddpd	%%zmm1,%%zmm1,%%zmm1	\n\t	vaddpd	%%zmm5,%%zmm5,%%zmm5\n\t	vaddpd	%%zmm9,%%zmm9,%%zmm9	\n\t	vaddpd	%%zmm13,%%zmm13,%%zmm13	\n\t"\
			"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm6,%%zmm4,%%zmm4\n\t	vmulpd	%%zmm10,%%zmm8,%%zmm8	\n\t	vmulpd	%%zmm14,%%zmm12,%%zmm12	\n\t"\
			"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5\n\t	vmulpd	%%zmm11,%%zmm9,%%zmm9	\n\t	vmulpd	%%zmm15,%%zmm13,%%zmm13	\n\t"\
			"vmovaps	%%zmm0,0x180(%%rax)	\n\t	vmovaps	%%zmm4,0x580(%%rax)	\n\t	vmovaps	%%zmm8,0x380(%%rax)		\n\t	vmovaps	%%zmm12,0x780(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x1c0(%%rax)	\n\t	vmovaps	%%zmm5,0x5c0(%%rax)	\n\t	vmovaps	%%zmm9,0x3c0(%%rax)		\n\t	vmovaps	%%zmm13,0x7c0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	// Clobbered registers
		);
	  #elif defined(USE_AVX)
		__asm__ volatile (\
			"movq	%[__r1],%%rax			\n\t"\
			"/* z0^2: */					\n\t		/* z8^2: */					\n\t"\
			"vmovaps		 (%%rax),%%ymm0	\n\t		vmovaps	0x200(%%rax),%%ymm4	\n\t"\
			"vmovaps	0x020(%%rax),%%ymm1	\n\t		vmovaps	0x220(%%rax),%%ymm5	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t		vmovaps		  %%ymm4,%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm3	\n\t		vmovaps		  %%ymm4,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
			"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t		vsubpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t		vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t		vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t		vmovaps	%%ymm4,0x200(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t		vmovaps	%%ymm5,0x220(%%rax)	\n\t"\
			"/* z1^2: */					\n\t		/* z9^2: */					\n\t"\
			"vmovaps	0x040(%%rax),%%ymm0	\n\t		vmovaps	0x240(%%rax),%%ymm4	\n\t"\
			"vmovaps	0x060(%%rax),%%ymm1	\n\t		vmovaps	0x260(%%rax),%%ymm5	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t		vmovaps		  %%ymm4,%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm3	\n\t		vmovaps		  %%ymm4,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
			"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t		vsubpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t		vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t		vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmovaps	%%ymm0,0x040(%%rax)	\n\t		vmovaps	%%ymm4,0x240(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x060(%%rax)	\n\t		vmovaps	%%ymm5,0x260(%%rax)	\n\t"\
			"/* z2^2: */					\n\t		/* zA^2: */					\n\t"\
			"vmovaps	0x080(%%rax),%%ymm0	\n\t		vmovaps	0x280(%%rax),%%ymm4	\n\t"\
			"vmovaps	0x0a0(%%rax),%%ymm1	\n\t		vmovaps	0x2a0(%%rax),%%ymm5	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t		vmovaps		  %%ymm4,%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm3	\n\t		vmovaps		  %%ymm4,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
			"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t		vsubpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t		vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t		vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmovaps	%%ymm0,0x080(%%rax)	\n\t		vmovaps	%%ymm4,0x280(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rax)	\n\t		vmovaps	%%ymm5,0x2a0(%%rax)	\n\t"\
			"/* z3^2: */					\n\t		/* zB^2: */					\n\t"\
			"vmovaps	0x0c0(%%rax),%%ymm0	\n\t		vmovaps	0x2c0(%%rax),%%ymm4	\n\t"\
			"vmovaps	0x0e0(%%rax),%%ymm1	\n\t		vmovaps	0x2e0(%%rax),%%ymm5	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t		vmovaps		  %%ymm4,%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm3	\n\t		vmovaps		  %%ymm4,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
			"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t		vsubpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t		vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t		vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmovaps	%%ymm0,0x0c0(%%rax)	\n\t		vmovaps	%%ymm4,0x2c0(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x0e0(%%rax)	\n\t		vmovaps	%%ymm5,0x2e0(%%rax)	\n\t"\
			"/* z4^2: */					\n\t		/* zC^2: */					\n\t"\
			"vmovaps	0x100(%%rax),%%ymm0	\n\t		vmovaps	0x300(%%rax),%%ymm4	\n\t"\
			"vmovaps	0x120(%%rax),%%ymm1	\n\t		vmovaps	0x320(%%rax),%%ymm5	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t		vmovaps		  %%ymm4,%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm3	\n\t		vmovaps		  %%ymm4,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
			"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t		vsubpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t		vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t		vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmovaps	%%ymm0,0x100(%%rax)	\n\t		vmovaps	%%ymm4,0x300(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x120(%%rax)	\n\t		vmovaps	%%ymm5,0x320(%%rax)	\n\t"\
			"/* z5^2: */					\n\t		/* zD^2: */					\n\t"\
			"vmovaps	0x140(%%rax),%%ymm0	\n\t		vmovaps	0x340(%%rax),%%ymm4	\n\t"\
			"vmovaps	0x160(%%rax),%%ymm1	\n\t		vmovaps	0x360(%%rax),%%ymm5	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t		vmovaps		  %%ymm4,%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm3	\n\t		vmovaps		  %%ymm4,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
			"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t		vsubpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t		vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t		vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmovaps	%%ymm0,0x140(%%rax)	\n\t		vmovaps	%%ymm4,0x340(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x160(%%rax)	\n\t		vmovaps	%%ymm5,0x360(%%rax)	\n\t"\
			"/* z6^2: */					\n\t		/* zE^2: */					\n\t"\
			"vmovaps	0x180(%%rax),%%ymm0	\n\t		vmovaps	0x380(%%rax),%%ymm4	\n\t"\
			"vmovaps	0x1a0(%%rax),%%ymm1	\n\t		vmovaps	0x3a0(%%rax),%%ymm5	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t		vmovaps		  %%ymm4,%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm3	\n\t		vmovaps		  %%ymm4,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
			"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t		vsubpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t		vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t		vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmovaps	%%ymm0,0x180(%%rax)	\n\t		vmovaps	%%ymm4,0x380(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x1a0(%%rax)	\n\t		vmovaps	%%ymm5,0x3a0(%%rax)	\n\t"\
			"/* z7^2: */					\n\t		/* zF^2: */					\n\t"\
			"vmovaps	0x1c0(%%rax),%%ymm0	\n\t		vmovaps	0x3c0(%%rax),%%ymm4	\n\t"\
			"vmovaps	0x1e0(%%rax),%%ymm1	\n\t		vmovaps	0x3e0(%%rax),%%ymm5	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t		vmovaps		  %%ymm4,%%ymm6	\n\t"\
			"vmovaps		  %%ymm0,%%ymm3	\n\t		vmovaps		  %%ymm4,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm0,%%ymm0	\n\t		vaddpd	%%ymm5,%%ymm4,%%ymm4	\n\t"\
			"vsubpd	%%ymm1,%%ymm2,%%ymm2	\n\t		vsubpd	%%ymm5,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t		vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t		vmulpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t		vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"vmovaps	%%ymm0,0x1c0(%%rax)	\n\t		vmovaps	%%ymm4,0x3c0(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x1e0(%%rax)	\n\t		vmovaps	%%ymm5,0x3e0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	// Clobbered registers
		);
	  #elif OS_BITS == 64
		__asm__ volatile (\
			"movq	%[__r1],%%rax		\n\t"\
			"/* z0^2: */				\n\t		/* z8^2: */				\n\t"\
			"movaps		 (%%rax),%%xmm0	\n\t		movaps	0x100(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t		movaps	0x110(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t		movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t		addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t		subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t		addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t		mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd		 (%%rax),%%xmm1	\n\t		mulpd	0x100(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t		movaps	%%xmm3,0x100(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t		movaps	%%xmm4,0x110(%%rax)	\n\t"\
			"/* z1^2: */				\n\t		/* z9^2: */				\n\t"\
			"movaps	0x020(%%rax),%%xmm0	\n\t		movaps	0x120(%%rax),%%xmm3	\n\t"\
			"movaps	0x030(%%rax),%%xmm1	\n\t		movaps	0x130(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t		movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t		addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t		subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t		addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t		mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x020(%%rax),%%xmm1	\n\t		mulpd	0x120(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x020(%%rax)	\n\t		movaps	%%xmm3,0x120(%%rax)	\n\t"\
			"movaps	%%xmm1,0x030(%%rax)	\n\t		movaps	%%xmm4,0x130(%%rax)	\n\t"\
			"/* z2^2: */				\n\t		/* zA^2: */				\n\t"\
			"movaps	0x040(%%rax),%%xmm0	\n\t		movaps	0x140(%%rax),%%xmm3	\n\t"\
			"movaps	0x050(%%rax),%%xmm1	\n\t		movaps	0x150(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t		movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t		addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t		subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t		addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t		mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x040(%%rax),%%xmm1	\n\t		mulpd	0x140(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x040(%%rax)	\n\t		movaps	%%xmm3,0x140(%%rax)	\n\t"\
			"movaps	%%xmm1,0x050(%%rax)	\n\t		movaps	%%xmm4,0x150(%%rax)	\n\t"\
			"/* z3^2: */				\n\t		/* zB^2: */				\n\t"\
			"movaps	0x060(%%rax),%%xmm0	\n\t		movaps	0x160(%%rax),%%xmm3	\n\t"\
			"movaps	0x070(%%rax),%%xmm1	\n\t		movaps	0x170(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t		movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t		addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t		subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t		addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t		mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x060(%%rax),%%xmm1	\n\t		mulpd	0x160(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x060(%%rax)	\n\t		movaps	%%xmm3,0x160(%%rax)	\n\t"\
			"movaps	%%xmm1,0x070(%%rax)	\n\t		movaps	%%xmm4,0x170(%%rax)	\n\t"\
			"/* z4^2: */				\n\t		/* zC^2: */				\n\t"\
			"movaps	0x080(%%rax),%%xmm0	\n\t		movaps	0x180(%%rax),%%xmm3	\n\t"\
			"movaps	0x090(%%rax),%%xmm1	\n\t		movaps	0x190(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t		movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t		addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t		subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t		addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t		mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x080(%%rax),%%xmm1	\n\t		mulpd	0x180(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x080(%%rax)	\n\t		movaps	%%xmm3,0x180(%%rax)	\n\t"\
			"movaps	%%xmm1,0x090(%%rax)	\n\t		movaps	%%xmm4,0x190(%%rax)	\n\t"\
			"/* z5^2: */				\n\t		/* zD^2: */				\n\t"\
			"movaps	0x0a0(%%rax),%%xmm0	\n\t		movaps	0x1a0(%%rax),%%xmm3	\n\t"\
			"movaps	0x0b0(%%rax),%%xmm1	\n\t		movaps	0x1b0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t		movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t		addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t		subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t		addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t		mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0a0(%%rax),%%xmm1	\n\t		mulpd	0x1a0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0a0(%%rax)	\n\t		movaps	%%xmm3,0x1a0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0b0(%%rax)	\n\t		movaps	%%xmm4,0x1b0(%%rax)	\n\t"\
			"/* z6^2: */				\n\t		/* zE^2: */				\n\t"\
			"movaps	0x0c0(%%rax),%%xmm0	\n\t		movaps	0x1c0(%%rax),%%xmm3	\n\t"\
			"movaps	0x0d0(%%rax),%%xmm1	\n\t		movaps	0x1d0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t		movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t		addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t		subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t		addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t		mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0c0(%%rax),%%xmm1	\n\t		mulpd	0x1c0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0c0(%%rax)	\n\t		movaps	%%xmm3,0x1c0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0d0(%%rax)	\n\t		movaps	%%xmm4,0x1d0(%%rax)	\n\t"\
			"/* z7^2: */				\n\t		/* zF^2: */				\n\t"\
			"movaps	0x0e0(%%rax),%%xmm0	\n\t		movaps	0x1e0(%%rax),%%xmm3	\n\t"\
			"movaps	0x0f0(%%rax),%%xmm1	\n\t		movaps	0x1f0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t		movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t		addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t		subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t		addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t		mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0e0(%%rax),%%xmm1	\n\t		mulpd	0x1e0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0e0(%%rax)	\n\t		movaps	%%xmm3,0x1e0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0f0(%%rax)	\n\t		movaps	%%xmm4,0x1f0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	// Clobbered registers
		);
	  #else	// OS_BITS = 32
		#error 32-bit OSes no longer supported for SIMD builds!
	  #endif	// avx512 / avx / sse2_64 / sse2_32 ?

	  }	// endif(fwd_fft_only == 1)

	  /*...And do an inverse DIT radix-16 pass on the squared-data blocks. */

	  #ifdef USE_AVX512
		SSE2_RADIX16_WRAPPER_DIT(add0,add1,add2,add3,add4,add5,add6,add7
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,pfetch_dist)
	  #elif defined(USE_AVX)
		SSE2_RADIX16_WRAPPER_DIT(add0,add1,add2,add3
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,pfetch_dist)
	  #else	// SSE2:
		SSE2_RADIX16_WRAPPER_DIT(add0,add1
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)
	  #endif // AVX or SSE2?

	#endif // FULLY_FUSED?

#else	/* if(!USE_SSE2) */

	if(fwd_fft_only == 3) {
		goto skip_fwd_fft;	// v20: jump-to-point for both-inputs-already-fwd-FFTed case
	}

	/*...Block 1: */
		t1 =a[j1   ];					t2 =a[j2   ];
		rt =a[j1+16]*c8 -a[j2+16]*s8 ;	it =a[j2+16]*c8 +a[j1+16]*s8;
		t3 =t1 -rt;	t1 =t1 +rt;
		t4 =t2 -it;	t2 =t2 +it;

		t5 =a[j1+8 ]*c4 -a[j2+8 ]*s4;	t6 =a[j2+8 ]*c4 +a[j1+8 ]*s4;
		rt =a[j1+24]*c12-a[j2+24]*s12;	it =a[j2+24]*c12+a[j1+24]*s12;
		t7 =t5 -rt;	t5 =t5 +rt;
		t8 =t6 -it;	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;		t1 =t1 +rt;
		it =t6;	t6 =t2 -it;		t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;		t3 =t3 -t8;
				t8 =t4 -rt;		t4 =t4 +rt;

	/*...Block 2: */
		t9 =a[j1+4 ]*c2 -a[j2+4 ]*s2 ;	t10=a[j2+4 ]*c2 +a[j1+4 ]*s2 ;
		rt =a[j1+20]*c10-a[j2+20]*s10;	it =a[j2+20]*c10+a[j1+20]*s10;
		t11=t9 -rt;	t9 =t9 +rt;
		t12=t10-it;	t10=t10+it;

		t13=a[j1+12]*c6 -a[j2+12]*s6 ;	t14=a[j2+12]*c6 +a[j1+12]*s6 ;
		rt =a[j1+28]*c14-a[j2+28]*s14;	it =a[j2+28]*c14+a[j1+28]*s14;
		t15=t13-rt;	t13=t13+rt;
		t16=t14-it;	t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11+t16;t11=t11-t16;
					t16=t12-rt;	t12=t12+rt;

	/*...Block 3: */
		t17=a[j1+2 ]*c1 -a[j2+2 ]*s1 ;	t18=a[j2+2 ]*c1 +a[j1+2 ]*s1 ;
		rt =a[j1+18]*c9 -a[j2+18]*s9 ;	it =a[j2+18]*c9 +a[j1+18]*s9 ;
		t19=t17-rt;	t17=t17+rt;
		t20=t18-it;	t18=t18+it;

		t21=a[j1+10]*c5 -a[j2+10]*s5 ;	t22=a[j2+10]*c5 +a[j1+10]*s5 ;
		rt =a[j1+26]*c13-a[j2+26]*s13;	it =a[j2+26]*c13+a[j1+26]*s13;
		t23=t21-rt;	t21=t21+rt;
		t24=t22-it;	t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19+t24;t19=t19-t24;
					t24=t20-rt;	t20=t20+rt;

	/*...Block 4: */
		t25=a[j1+6 ]*c3 -a[j2+6 ]*s3 ;	t26=a[j2+6 ]*c3 +a[j1+6 ]*s3 ;
		rt =a[j1+22]*c11-a[j2+22]*s11;	it =a[j2+22]*c11+a[j1+22]*s11;
		t27=t25-rt;	t25=t25+rt;
		t28=t26-it;	t26=t26+it;

		t29=a[j1+14]*c7 -a[j2+14]*s7 ;	t30=a[j2+14]*c7 +a[j1+14]*s7 ;
		rt =a[j1+30]*c15-a[j2+30]*s15;	it =a[j2+30]*c15+a[j1+30]*s15;
		t31=t29-rt;	t29=t29+rt;
		t32=t30-it;	t30=t30+it;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27+t32;t27=t27-t32;
					t32=t28-rt;	t28=t28+rt;

	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(i* 1*twopi/16) =       ( c, s), exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 3*twopi/16) =       ( s, c) (for inputs to transform block 2)
	!	1, exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 4*twopi/16) =       ( 0, 1), exp(i* 6*twopi/16) = ISRT2*(-1, 1) (for inputs to transform block 3)
	!	1, exp(i* 3*twopi/16) =       ( s, c), exp(i* 6*twopi/16) = ISRT2*(-1, 1), exp(i* 9*twopi/16) =       (-c,-s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1  +p0:15:1] are replaced by t0:30:2,
	!						 a[j2+p0:15:1] are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
		rt =t9;		t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		aj1p0r =t1 +t17;	aj1p0i =t2 +t18;
		aj1p1r =t1 -t17;	aj1p1i =t2 -t18;

		aj1p2r =t9 -t26;	aj1p2i =t10+t25;
		aj1p3r =t9 +t26;	aj1p3i =t10-t25;
	/*...Block 3: t5,13,21,29 */
		rt =t13;	t13=t5 +t14;t5 =t5 -t14;
					t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;	t22=(t21+t22)*ISRT2;	t21=rt;
		rt =(t30+t29)*ISRT2;	it =(t30-t29)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		aj1p4r =t5 +t21;	aj1p4i =t6 +t22;
		aj1p5r =t5 -t21;	aj1p5i =t6 -t22;

		aj1p6r =t13-t30;	aj1p6i =t14+t29;
		aj1p7r =t13+t30;	aj1p7i =t14-t29;

	/*...Block 2: t3,11,19,27 */
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;
		rt =t27*s - t28*c;	it =t28*s + t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		aj1p8r =t3 +t19;	aj1p8i =t4 +t20;
		aj1p9r =t3 -t19;	aj1p9i =t4 -t20;

		aj1pAr=t11-t28;	aj1pAi=t12+t27;
		aj1pBr=t11+t28;	aj1pBi=t12-t27;

	/*...Block 4: t7,15,23,31 */
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;
		rt =t31*c - t32*s;	it =t32*c + t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		aj1pCr=t7 +t23;	aj1pCi=t8 +t24;
		aj1pDr=t7 -t23;	aj1pDi=t8 -t24;

		aj1pEr=t15-t32;	aj1pEi=t16+t31;
		aj1pFr=t15+t32;	aj1pFi=t16-t31;

skip_fwd_fft:	// v20: jump-to-point for both-inputs-already-fwd-FFTed case

	// v20: If fwd_fft_only = 1, write fwd-FFT result back to input array, skipping dyadic-square and inv-FFT steps:
	if(fwd_fft_only == 1)
	{
	/*...Block 1: t1,9,17,25 */
		a[j1   ] = aj1p0r;	a[j2   ] = aj1p0i;
		a[j1+16] = aj1p1r;	a[j2+16] = aj1p1i;
		a[j1+8 ] = aj1p2r;	a[j2+8 ] = aj1p2i;
		a[j1+24] = aj1p3r;	a[j2+24] = aj1p3i;
	/*...Block 3: t5,13,21,29 */
		a[j1+4 ] = aj1p4r;	a[j2+4 ] = aj1p4i;
		a[j1+20] = aj1p5r;	a[j2+20] = aj1p5i;
		a[j1+12] = aj1p6r;	a[j2+12] = aj1p6i;
		a[j1+28] = aj1p7r;	a[j2+28] = aj1p7i;
	/*...Block 2: t3,11,19,27 */
		a[j1+2 ] = aj1p8r;	a[j2+2 ] = aj1p8i;
		a[j1+18] = aj1p9r;	a[j2+18] = aj1p9i;
		a[j1+10] = aj1pAr;	a[j2+10] = aj1pAi;
		a[j1+26] = aj1pBr;	a[j2+26] = aj1pBi;
	/*...Block 4: t7,15,23,31 */
		a[j1+6 ] = aj1pCr;	a[j2+6 ] = aj1pCi;
		a[j1+22] = aj1pDr;	a[j2+22] = aj1pDi;
		a[j1+14] = aj1pEr;	a[j2+14] = aj1pEi;
		a[j1+30] = aj1pFr;	a[j2+30] = aj1pFi;
		goto loop;	// Skip dyadic-mul and iFFT

	} else if (fwd_fft_only) {	// 2-input modmul: fwdFFT data dyadic-mul'ed with precomputed 2nd-vector stored in fwd-FFTed form in b[]:
	 	if(fwd_fft_only == 3) {	// v20: Both inputs enter fwd-FFTed, must copy from main arrays a[],b[] to local-vars
		/*...Block 1: t1,9,17,25 */
			aj1p0r = a[j1   ];	aj1p0i = a[j2   ];
			aj1p1r = a[j1+16];	aj1p1i = a[j2+16];
			aj1p2r = a[j1+8 ];	aj1p2i = a[j2+8 ];
			aj1p3r = a[j1+24];	aj1p3i = a[j2+24];
		/*...Block 3: t5,13,21,29 */
			aj1p4r = a[j1+4 ];	aj1p4i = a[j2+4 ];
			aj1p5r = a[j1+20];	aj1p5i = a[j2+20];
			aj1p6r = a[j1+12];	aj1p6i = a[j2+12];
			aj1p7r = a[j1+28];	aj1p7i = a[j2+28];
		/*...Block 2: t3,11,19,27 */
			aj1p8r = a[j1+2 ];	aj1p8i = a[j2+2 ];
			aj1p9r = a[j1+18];	aj1p9i = a[j2+18];
			aj1pAr = a[j1+10];	aj1pAi = a[j2+10];
			aj1pBr = a[j1+26];	aj1pBi = a[j2+26];
		/*...Block 4: t7,15,23,31 */
			aj1pCr = a[j1+6 ];	aj1pCi = a[j2+6 ];
			aj1pDr = a[j1+22];	aj1pDi = a[j2+22];
			aj1pEr = a[j1+14];	aj1pEi = a[j2+14];
			aj1pFr = a[j1+30];	aj1pFi = a[j2+30];
		}

		if(c_arr) {	// c_arr != 0x0: a * (b - c)

		/*...Block 1: t1,9,17,25 */
			bj1p0r = b[j1   ] - c_arr[j1   ];	bj1p0i = b[j2   ] - c_arr[j2   ];
			bj1p1r = b[j1+16] - c_arr[j1+16];	bj1p1i = b[j2+16] - c_arr[j2+16];
			bj1p2r = b[j1+8 ] - c_arr[j1+8 ];	bj1p2i = b[j2+8 ] - c_arr[j2+8 ];
			bj1p3r = b[j1+24] - c_arr[j1+24];	bj1p3i = b[j2+24] - c_arr[j2+24];
		/*...Block 3: t5,13,21,29 */
			bj1p4r = b[j1+4 ] - c_arr[j1+4 ];	bj1p4i = b[j2+4 ] - c_arr[j2+4 ];
			bj1p5r = b[j1+20] - c_arr[j1+20];	bj1p5i = b[j2+20] - c_arr[j2+20];
			bj1p6r = b[j1+12] - c_arr[j1+12];	bj1p6i = b[j2+12] - c_arr[j2+12];
			bj1p7r = b[j1+28] - c_arr[j1+28];	bj1p7i = b[j2+28] - c_arr[j2+28];
		/*...Block 2: t3,11,19,27 */
			bj1p8r = b[j1+2 ] - c_arr[j1+2 ];	bj1p8i = b[j2+2 ] - c_arr[j2+2 ];
			bj1p9r = b[j1+18] - c_arr[j1+18];	bj1p9i = b[j2+18] - c_arr[j2+18];
			bj1pAr = b[j1+10] - c_arr[j1+10];	bj1pAi = b[j2+10] - c_arr[j2+10];
			bj1pBr = b[j1+26] - c_arr[j1+26];	bj1pBi = b[j2+26] - c_arr[j2+26];
		/*...Block 4: t7,15,23,31 */
			bj1pCr = b[j1+6 ] - c_arr[j1+6 ];	bj1pCi = b[j2+6 ] - c_arr[j2+6 ];
			bj1pDr = b[j1+22] - c_arr[j1+22];	bj1pDi = b[j2+22] - c_arr[j2+22];
			bj1pEr = b[j1+14] - c_arr[j1+14];	bj1pEi = b[j2+14] - c_arr[j2+14];
			bj1pFr = b[j1+30] - c_arr[j1+30];	bj1pFi = b[j2+30] - c_arr[j2+30];

		} else {	// c_arr = 0x0: a * b

		/*...Block 1: t1,9,17,25 */
			bj1p0r = b[j1   ];	bj1p0i = b[j2   ];
			bj1p1r = b[j1+16];	bj1p1i = b[j2+16];
			bj1p2r = b[j1+8 ];	bj1p2i = b[j2+8 ];
			bj1p3r = b[j1+24];	bj1p3i = b[j2+24];
		/*...Block 3: t5,13,21,29 */
			bj1p4r = b[j1+4 ];	bj1p4i = b[j2+4 ];
			bj1p5r = b[j1+20];	bj1p5i = b[j2+20];
			bj1p6r = b[j1+12];	bj1p6i = b[j2+12];
			bj1p7r = b[j1+28];	bj1p7i = b[j2+28];
		/*...Block 2: t3,11,19,27 */
			bj1p8r = b[j1+2 ];	bj1p8i = b[j2+2 ];
			bj1p9r = b[j1+18];	bj1p9i = b[j2+18];
			bj1pAr = b[j1+10];	bj1pAi = b[j2+10];
			bj1pBr = b[j1+26];	bj1pBi = b[j2+26];
		/*...Block 4: t7,15,23,31 */
			bj1pCr = b[j1+6 ];	bj1pCi = b[j2+6 ];
			bj1pDr = b[j1+22];	bj1pDi = b[j2+22];
			bj1pEr = b[j1+14];	bj1pEi = b[j2+14];
			bj1pFr = b[j1+30];	bj1pFi = b[j2+30];

		}	// endif(c_arr)

	/*...Dyadic mul of the forward FFT outputs: (ar + I.ai)*(br + I.bi) = (ar.br - ai.bi) + I.(ar.bi + ai.br) */
		rt = aj1p0r;	aj1p0r = (aj1p0r*bj1p0r - aj1p0i*bj1p0i);	aj1p0i = (rt*bj1p0i + aj1p0i*bj1p0r);
		rt = aj1p1r;	aj1p1r = (aj1p1r*bj1p1r - aj1p1i*bj1p1i);	aj1p1i = (rt*bj1p1i + aj1p1i*bj1p1r);
		rt = aj1p2r;	aj1p2r = (aj1p2r*bj1p2r - aj1p2i*bj1p2i);	aj1p2i = (rt*bj1p2i + aj1p2i*bj1p2r);
		rt = aj1p3r;	aj1p3r = (aj1p3r*bj1p3r - aj1p3i*bj1p3i);	aj1p3i = (rt*bj1p3i + aj1p3i*bj1p3r);
		rt = aj1p4r;	aj1p4r = (aj1p4r*bj1p4r - aj1p4i*bj1p4i);	aj1p4i = (rt*bj1p4i + aj1p4i*bj1p4r);
		rt = aj1p5r;	aj1p5r = (aj1p5r*bj1p5r - aj1p5i*bj1p5i);	aj1p5i = (rt*bj1p5i + aj1p5i*bj1p5r);
		rt = aj1p6r;	aj1p6r = (aj1p6r*bj1p6r - aj1p6i*bj1p6i);	aj1p6i = (rt*bj1p6i + aj1p6i*bj1p6r);
		rt = aj1p7r;	aj1p7r = (aj1p7r*bj1p7r - aj1p7i*bj1p7i);	aj1p7i = (rt*bj1p7i + aj1p7i*bj1p7r);
		rt = aj1p8r;	aj1p8r = (aj1p8r*bj1p8r - aj1p8i*bj1p8i);	aj1p8i = (rt*bj1p8i + aj1p8i*bj1p8r);
		rt = aj1p9r;	aj1p9r = (aj1p9r*bj1p9r - aj1p9i*bj1p9i);	aj1p9i = (rt*bj1p9i + aj1p9i*bj1p9r);
		rt = aj1pAr;	aj1pAr = (aj1pAr*bj1pAr - aj1pAi*bj1pAi);	aj1pAi = (rt*bj1pAi + aj1pAi*bj1pAr);
		rt = aj1pBr;	aj1pBr = (aj1pBr*bj1pBr - aj1pBi*bj1pBi);	aj1pBi = (rt*bj1pBi + aj1pBi*bj1pBr);
		rt = aj1pCr;	aj1pCr = (aj1pCr*bj1pCr - aj1pCi*bj1pCi);	aj1pCi = (rt*bj1pCi + aj1pCi*bj1pCr);
		rt = aj1pDr;	aj1pDr = (aj1pDr*bj1pDr - aj1pDi*bj1pDi);	aj1pDi = (rt*bj1pDi + aj1pDi*bj1pDr);
		rt = aj1pEr;	aj1pEr = (aj1pEr*bj1pEr - aj1pEi*bj1pEi);	aj1pEi = (rt*bj1pEi + aj1pEi*bj1pEr);
		rt = aj1pFr;	aj1pFr = (aj1pFr*bj1pFr - aj1pFi*bj1pFi);	aj1pFi = (rt*bj1pFi + aj1pFi*bj1pFr);

	} else {	// fwd_fft_only = 0: Normal execution, dyadic-squaring followed by iFFT:

	/*...Dyadic square of the forward FFT outputs: */
		rt = aj1p0r*aj1p0i;	aj1p0r = (aj1p0r + aj1p0i)*(aj1p0r - aj1p0i);	aj1p0i = (rt + rt);
		rt = aj1p1r*aj1p1i;	aj1p1r = (aj1p1r + aj1p1i)*(aj1p1r - aj1p1i);	aj1p1i = (rt + rt);
		rt = aj1p2r*aj1p2i;	aj1p2r = (aj1p2r + aj1p2i)*(aj1p2r - aj1p2i);	aj1p2i = (rt + rt);
		rt = aj1p3r*aj1p3i;	aj1p3r = (aj1p3r + aj1p3i)*(aj1p3r - aj1p3i);	aj1p3i = (rt + rt);
		rt = aj1p4r*aj1p4i;	aj1p4r = (aj1p4r + aj1p4i)*(aj1p4r - aj1p4i);	aj1p4i = (rt + rt);
		rt = aj1p5r*aj1p5i;	aj1p5r = (aj1p5r + aj1p5i)*(aj1p5r - aj1p5i);	aj1p5i = (rt + rt);
		rt = aj1p6r*aj1p6i;	aj1p6r = (aj1p6r + aj1p6i)*(aj1p6r - aj1p6i);	aj1p6i = (rt + rt);
		rt = aj1p7r*aj1p7i;	aj1p7r = (aj1p7r + aj1p7i)*(aj1p7r - aj1p7i);	aj1p7i = (rt + rt);
		rt = aj1p8r*aj1p8i;	aj1p8r = (aj1p8r + aj1p8i)*(aj1p8r - aj1p8i);	aj1p8i = (rt + rt);
		rt = aj1p9r*aj1p9i;	aj1p9r = (aj1p9r + aj1p9i)*(aj1p9r - aj1p9i);	aj1p9i = (rt + rt);
		rt = aj1pAr*aj1pAi;	aj1pAr = (aj1pAr + aj1pAi)*(aj1pAr - aj1pAi);	aj1pAi = (rt + rt);
		rt = aj1pBr*aj1pBi;	aj1pBr = (aj1pBr + aj1pBi)*(aj1pBr - aj1pBi);	aj1pBi = (rt + rt);
		rt = aj1pCr*aj1pCi;	aj1pCr = (aj1pCr + aj1pCi)*(aj1pCr - aj1pCi);	aj1pCi = (rt + rt);
		rt = aj1pDr*aj1pDi;	aj1pDr = (aj1pDr + aj1pDi)*(aj1pDr - aj1pDi);	aj1pDi = (rt + rt);
		rt = aj1pEr*aj1pEi;	aj1pEr = (aj1pEr + aj1pEi)*(aj1pEr - aj1pEi);	aj1pEi = (rt + rt);
		rt = aj1pFr*aj1pFi;	aj1pFr = (aj1pFr + aj1pFi)*(aj1pFr - aj1pFi);	aj1pFi = (rt + rt);

	}	// endif(fwd_fft_only == 1)

	/*********************************************************************/
	/*...And do an inverse DIT radix-16 pass on the squared-data       : */
	/*********************************************************************/
	#if PFETCH
	add0 = &a[j1+32];
	#endif
	/*   gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 iDIT transforms... */

	/*...Block 1: */
	#if PFETCH
	prefetch_p_doubles(add0);
	add0 += 4;
	#endif
		t3 =aj1p0r -aj1p1r ;	t4 =aj1p0i -aj1p1i ;
		t1 =aj1p0r +aj1p1r ;	t2 =aj1p0i +aj1p1i ;

		t7 =aj1p2r -aj1p3r ;	t8 =aj1p2i -aj1p3i ;
		t5 =aj1p2r +aj1p3r ;	t6 =aj1p2i +aj1p3i ;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
		t8 =t4 +rt;	t4 =t4 -rt;

	/*...Block 2: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(add0);
		#endif
	add0 += 4;
	#endif
		t11=aj1p4r -aj1p5r ;	t12=aj1p4i -aj1p5i ;
		t9 =aj1p4r +aj1p5r ;	t10=aj1p4i +aj1p5i ;

		t15=aj1p6r -aj1p7r ;	t16=aj1p6i -aj1p7i ;
		t13=aj1p6r +aj1p7r ;	t14=aj1p6i +aj1p7i ;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11-t16;	t11=t11+t16;
		t16=t12+rt;	t12=t12-rt;

	/*...Block 3: */
	#if PFETCH
	prefetch_p_doubles(add0);
	add0 += 4;
	#endif
		t19=aj1p8r -aj1p9r ;	t20=aj1p8i -aj1p9i ;
		t17=aj1p8r +aj1p9r ;	t18=aj1p8i +aj1p9i ;

		t23=aj1pAr-aj1pBr;	t24=aj1pAi-aj1pBi;
		t21=aj1pAr+aj1pBr;	t22=aj1pAi+aj1pBi;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19-t24;	t19=t19+t24;
		t24=t20+rt;	t20=t20-rt;

	/*...Block 4: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(add0);
		#endif
	add0 += 4;
	#endif
		t27=aj1pCr-aj1pDr;	t28=aj1pCi-aj1pDi;
		t25=aj1pCr+aj1pDr;	t26=aj1pCi+aj1pDi;

		t31=aj1pEr-aj1pFr;	t32=aj1pEi-aj1pFi;
		t29=aj1pEr+aj1pFr;	t30=aj1pEi+aj1pFi;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27-t32;	t27=t27+t32;
		t32=t28+rt;	t28=t28-rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
	!	1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
	!	1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1  +p0:15:1] are replaced by t0:30:2,
	!						 a[j2+p0:15:1] are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
	#if PFETCH
	prefetch_p_doubles(add0);
	add0 += 4;
	#endif
		rt =t9;		t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[j1   ]=t1+t17;				a[j2   ]=t2+t18;
		t1	=t1-t17;					t2	=t2-t18;
		a[j1+16]=t1 *c8 +t2 *s8 ;	a[j2+16]=t2 *c8 -t1 *s8;

		rt	=t9 +t26;					it	=t10-t25;
		t9	=t9 -t26;					t10	=t10+t25;
		a[j1+8 ]=rt *c4 +it *s4 ;	a[j2+8 ]=it *c4 -rt *s4;
		a[j1+24]=t9 *c12+t10*s12;	a[j2+24]=t10*c12-t9 *s12;

	/*...Block 3: t5,13,21,29 */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(add0);
		#endif
	add0 += 4;
	#endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;
					t14=t6 +rt;		t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;	t22=(t22-t21)*ISRT2;	t21=rt;
		rt =(t29-t30)*ISRT2;	it =(t29+t30)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		rt	=t5 +t21;					it	=t6 +t22;
		t5	=t5 -t21;					t6	=t6 -t22;
		a[j1+4 ]=rt *c2 +it *s2 ;	a[j2+4 ]=it *c2 -rt *s2 ;
		a[j1+20]=t5 *c10+t6 *s10;	a[j2+20]=t6 *c10-t5 *s10;

		rt	=t13+t30;					it	=t14-t29;
		t13	=t13-t30;					t14	=t14+t29;
		a[j1+12]=rt *c6 +it *s6 ;	a[j2+12]=it *c6 -rt *s6 ;
		a[j1+28]=t13*c14+t14*s14;	a[j2+28]=t14*c14-t13*s14;

	/*...Block 2: t3,11,19,27 */
	#if PFETCH
	prefetch_p_doubles(add0);
	add0 += 4;
	#endif
		rt =(t12+t11)*ISRT2;	it =(t12-t11)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;
		rt =t27*s + t28*c;	it =t28*s - t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		rt	=t3 +t19;					it	=t4 +t20;
		t3	=t3 -t19;					t4	=t4 -t20;
		a[j1+2 ]=rt *c1 +it *s1 ;	a[j2+2 ]=it *c1 -rt *s1 ;
		a[j1+18]=t3 *c9 +t4 *s9 ;	a[j2+18]=t4 *c9 -t3 *s9 ;

		rt	=t11+t28;					it	=t12-t27;
		t11	=t11-t28;					t12	=t12+t27;
		a[j1+10]=rt *c5 +it *s5 ;	a[j2+10]=it *c5 -rt *s5 ;
		a[j1+26]=t11*c13+t12*s13;	a[j2+26]=t12*c13-t11*s13;

	/*...Block 4: t7,15,23,31 */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(add0);
		#endif
	add0 += 4;
	#endif
		rt =(t15-t16)*ISRT2;	it =(t15+t16)*ISRT2;
		rt =(t15-t16)*ISRT2;	it =(t15+t16)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;
		rt =t31*c + t32*s;	it =t32*c - t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		rt	=t7 +t23;					it	=t8 +t24;
		t7	=t7 -t23;					t8	=t8 -t24;
		a[j1+6 ]=rt *c3 +it *s3 ;	a[j2+6 ]=it *c3 -rt *s3 ;
		a[j1+22]=t7 *c11+t8 *s11;	a[j2+22]=t8 *c11-t7 *s11;

		rt	=t15+t32;					it	=t16-t31;
		t15	=t15-t32;					t16	=t16+t31;
		a[j1+14]=rt *c7 +it *s7 ;	a[j2+14]=it *c7 -rt *s7 ;
		a[j1+30]=t15*c15+t16*s15;	a[j2+30]=t16*c15-t15*s15;
#endif	/* USE_SSE2 */
	loop:
		continue;	// Need a statement following the jump-target label 'loop' here
	}	/* endfor() */

}

#undef PFETCH_DIST
