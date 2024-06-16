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
#include "pair_square.h"

#ifdef USE_SSE2

	#include "sse2_macro_gcc64.h"
	#include "radix32_utils_asm.h"

	#ifndef COMPILER_TYPE_GCC
		#error X86 SIMD build requires GCC-compatible compiler!
	#endif	/* GCC-style inline ASM: */

	#include "radix32_wrapper_square_gcc64.h"


#endif	/* USE_SSE2 */

/***************/

/* v19: To support Gerbicz check (and later, p-1) need 2-input FFT-mul - added fwd_fft_only flag, which takes 3 distinct possible values:
	  0 - Do normal FFT-autosquare(a) for (ilo-ihi) successive iterations;
	  1 - Do forward FFT(a) only and store result in a[], skipping dyadic-square and inv-FFT steps. In this case ilo = ihi;
	> 1 - The fwd_fft_only arg contains a pointer to an array b[] in previously-forward-FFTed form, i.e. FFT(b). In this
		  case we do the final-radix fwd-FFT pass of FFT(a), then dyadic-multiply FFT(a) * FFT(b) (storing the result in a[]) in place
		  of the usual pair_square step, then do the initial-radix pass of the iFFT of the result.
*/
void radix32_wrapper_square(
	double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[],
	int nradices_prim, int radix_prim[], int ws_i, int ws_j1, int ws_j2, int ws_j2_start, int ws_k, int ws_m, int ws_blocklen, int ws_blocklen_sum,
	int init_sse2, int thr_id, uint64 fwd_fft_only, double c_arr[]
)
{
	const char func[] = "radix32_wrapper_square";
/*
!   NOTE: In the following commentary, N refers to the COMPLEX vector length (N2 in the code),
!   which is half the real vector length.
!
!...Acronyms:   DIT = Decimation In Time
!               DIF = Decimation In Frequency
!               FFT = Fast Fourier Transform, i.e. a discrete FT over the complex numbers
!
!...Complex--->Real wrapper for data returned by forward transform, pointwise squaring and inverse wrapper
!   prior to inverse transform, all performed in one fell swoop (or better, one swell loop :)
!
!   On entry, the A-array contains floating data resulting from a forward DIF FFT, i.e. in bit-reversed order.
!
!   Quite a few operations are eliminated via use of algebraic identities in the derivation of the squaring
!   formula: the standard version of the algorithm (see, e.g., Numerical Recipes)
!   takes (ordered) forward FFT outputs H[0,...,N-1] and uses the following formula
!   to generate the N nonredundant terms F[0,...,N-1] of the DFT of the (true) real input signal:
!
!   F[j] = {( H[j]+H~[N-j] ) - I*exp(+pi*I*j/N)*( H[j]+H~[N-j] )}/2,
!
!       j=0, ... , N-1.        (Complex forward wrapper step)
!
!   Then one does a pointwise squaring to obtain terms G[0,...,N-1],
!   and does an inverse complex wrapper prior to entering the inverse FFT:
!
!   I[j] = {( G[j]+G~[N-j] ) + I*exp(-pi*I*j/N)*( G[j]-G~[N-j] )}/2,
!
!       j=0, ... , N-1.
!
!   By combining the 3 steps into one and doing some algebra, one
!   can get directly from the H's to the I's via
!
!   I[j] = H[j]^2 + {1 + exp(2*pi*I*j/N)}*{H[j]-H~[N-j]}^2/4,
!
!       j=0, ... , N-1,
!
!   and of course there are lots of symmetries in the calculation of
!   the data pairs I[j] and I[N-j] which can be (and are) exploited.
!
!   When the data in question are bit-reversed rather than ordered, it is
!   a nontrivial problem as to how to effect the above in a cache-friendly manner.
!   The problem reduces to figuring out what happens to (j, N-j) index pairs
!   under bit-reversal reordering. Happily this problem has an elegant solution,
!   which amounts to bit-reversing the array of roots of unity (i.e. the exp(2*pi*I*j/N)
!   in the last expression above, since this issue only affects the floating data)
!   and processing the bit-reversed data in the A-array (the floating-point H[j]'s)
!   in a sequence of smaller (but contiguous) sub-blocks - in each of the sub-blocks,
!   one array index increases monotonically, and a second (pointing to the bit-reversed
!   H[N-j]'s) decreases monotonically. Thus, things are very much like in the ordered-
!   data (j, N-j) index scheme, except that we have O(log2(N)) blocks, and the (k)th block
!   has index pairs (j, N'(k)-j). Further details on this may be found below, under the
!   heading "SOLVING THE CACHE FLOW PROBLEM...".
*/
#ifdef USE_SSE2
	const int stride = (int)RE_IM_STRIDE << 5;	// main-array loop stride = 64 for sse2, 128 for avx, 256 for AVX-512
#else
	const int stride = 64;	// In this particular routine, scalar mode has same stride as SSE2
#endif
	static int max_threads = 0;
	static int nsave = 0;
	static int *index = 0x0, *index_ptmp = 0x0;	/* N2/32-length Bit-reversal index array. */
	int *itmp = 0x0;
#if PFETCH
	double *addr;
#endif
	int rdum,idum,j1pad,j2pad,kp,l,iroot,k1,k2;
	int i,j1,j2,j2_start,k,m,blocklen,blocklen_sum,nbytes;
	const double c     = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
				,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784		/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
				,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;	/* exp(3*i*twopi/32)	*/
	double rt,it,re = 0.0, im= 0.0,dmult;
	double re0,im0,re1,im1;
	// sincos for 1st set of scalar-mode inputs
	double cA01,cA02,cA03,cA04,cA05,cA06,cA07,cA08,cA09,cA0A,cA0B,cA0C,cA0D,cA0E,cA0F,cA10,cA11,cA12,cA13,cA14,cA15,cA16,cA17,cA18,cA19,cA1A,cA1B,cA1C,cA1D,cA1E,cA1F
		,sA01,sA02,sA03,sA04,sA05,sA06,sA07,sA08,sA09,sA0A,sA0B,sA0C,sA0D,sA0E,sA0F,sA10,sA11,sA12,sA13,sA14,sA15,sA16,sA17,sA18,sA19,sA1A,sA1B,sA1C,sA1D,sA1E,sA1F;
	// sincos for 2nd set of scalar-mode inputs
	double cB01,cB02,cB03,cB04,cB05,cB06,cB07,cB08,cB09,cB0A,cB0B,cB0C,cB0D,cB0E,cB0F,cB10,cB11,cB12,cB13,cB14,cB15,cB16,cB17,cB18,cB19,cB1A,cB1B,cB1C,cB1D,cB1E,cB1F
		,sB01,sB02,sB03,sB04,sB05,sB06,sB07,sB08,sB09,sB0A,sB0B,sB0C,sB0D,sB0E,sB0F,sB10,sB11,sB12,sB13,sB14,sB15,sB16,sB17,sB18,sB19,sB1A,sB1B,sB1C,sB1D,sB1E,sB1F;
	double t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
	,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
	,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
	,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p0Ar,a1p0Br,a1p0Cr,a1p0Dr,a1p0Er,a1p0Fr
	,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p1Ar,a1p1Br,a1p1Cr,a1p1Dr,a1p1Er,a1p1Fr
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p0Ai,a1p0Bi,a1p0Ci,a1p0Di,a1p0Ei,a1p0Fi
	,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p1Ai,a1p1Bi,a1p1Ci,a1p1Di,a1p1Ei,a1p1Fi
	,a2p00r,a2p01r,a2p02r,a2p03r,a2p04r,a2p05r,a2p06r,a2p07r,a2p08r,a2p09r,a2p0Ar,a2p0Br,a2p0Cr,a2p0Dr,a2p0Er,a2p0Fr
	,a2p10r,a2p11r,a2p12r,a2p13r,a2p14r,a2p15r,a2p16r,a2p17r,a2p18r,a2p19r,a2p1Ar,a2p1Br,a2p1Cr,a2p1Dr,a2p1Er,a2p1Fr
	,a2p00i,a2p01i,a2p02i,a2p03i,a2p04i,a2p05i,a2p06i,a2p07i,a2p08i,a2p09i,a2p0Ai,a2p0Bi,a2p0Ci,a2p0Di,a2p0Ei,a2p0Fi
	,a2p10i,a2p11i,a2p12i,a2p13i,a2p14i,a2p15i,a2p16i,a2p17i,a2p18i,a2p19i,a2p1Ai,a2p1Bi,a2p1Ci,a2p1Di,a2p1Ei,a2p1Fi
	,b1p00r,b1p01r,b1p02r,b1p03r,b1p04r,b1p05r,b1p06r,b1p07r,b1p08r,b1p09r,b1p0Ar,b1p0Br,b1p0Cr,b1p0Dr,b1p0Er,b1p0Fr
	,b1p10r,b1p11r,b1p12r,b1p13r,b1p14r,b1p15r,b1p16r,b1p17r,b1p18r,b1p19r,b1p1Ar,b1p1Br,b1p1Cr,b1p1Dr,b1p1Er,b1p1Fr
	,b1p00i,b1p01i,b1p02i,b1p03i,b1p04i,b1p05i,b1p06i,b1p07i,b1p08i,b1p09i,b1p0Ai,b1p0Bi,b1p0Ci,b1p0Di,b1p0Ei,b1p0Fi
	,b1p10i,b1p11i,b1p12i,b1p13i,b1p14i,b1p15i,b1p16i,b1p17i,b1p18i,b1p19i,b1p1Ai,b1p1Bi,b1p1Ci,b1p1Di,b1p1Ei,b1p1Fi
	,b2p00r,b2p01r,b2p02r,b2p03r,b2p04r,b2p05r,b2p06r,b2p07r,b2p08r,b2p09r,b2p0Ar,b2p0Br,b2p0Cr,b2p0Dr,b2p0Er,b2p0Fr
	,b2p10r,b2p11r,b2p12r,b2p13r,b2p14r,b2p15r,b2p16r,b2p17r,b2p18r,b2p19r,b2p1Ar,b2p1Br,b2p1Cr,b2p1Dr,b2p1Er,b2p1Fr
	,b2p00i,b2p01i,b2p02i,b2p03i,b2p04i,b2p05i,b2p06i,b2p07i,b2p08i,b2p09i,b2p0Ai,b2p0Bi,b2p0Ci,b2p0Di,b2p0Ei,b2p0Fi
	,b2p10i,b2p11i,b2p12i,b2p13i,b2p14i,b2p15i,b2p16i,b2p17i,b2p18i,b2p19i,b2p1Ai,b2p1Bi,b2p1Ci,b2p1Di,b2p1Ei,b2p1Fi;

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error X86 SIMD code not supported for this compiler!
  #endif
		// The AVX512 compute-sincos-mults code needs 2 elements per complex-double-load, so use 14*RE_IM_STRIDE
		// to alloc storage here for all cases, even though that leaves upper array halves unused for sub-AVX512.
	static uint32 *sm_arr = 0x0,*sm_ptr;	// Base-ptr to arrays of k1,k2-index-vectors used in SIMD roots-computation.
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0,*add1, *bdd0,*bdd1, *cdd0,*cdd1;	/* Addresses into array sections */
  #ifdef USE_AVX
	double *add2,*add3, *bdd2,*bdd3, *cdd2,*cdd3;
  #endif
  #ifdef USE_AVX512
	double *add4,*add5,*add6,*add7, *bdd4,*bdd5,*bdd6,*bdd7, *cdd4,*cdd5,*cdd6,*cdd7;
  #endif
	vec_dbl *tmp,*tm1, *c_tmp,*s_tmp, *bpt0,*bpt1,*bpt2,*bpt3;
  #ifdef USE_AVX2
	vec_dbl *bpt4,*bpt5,*bpt6,*bpt7;
  #endif

  #ifdef MULTITHREAD
	// Base addresses for discrete per-thread local stores ... 'r' for double-float data, 'i' for int:
	static vec_dbl *__r0;
	static uint32  *__i0;
	// In || mode, only above base-pointer (shared by all threads) is static:
	uint32  *k1_arr, *k2_arr;
	vec_dbl *isrt2,*sqrt2,*one,*two, *cc0,*ss0,*cc1,*ss1,*cc3,*ss3
		,*forth,*tmp0,*tmp1,*tmp2,*tmp3,*tmp4,*tmp5,*tmp6,*tmp7
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*r00,*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E
		,*r10,*r12,*r14,*r16,*r18,*r1A,*r1C,*r1E
		,*r20,*r22,*r24,*r26,*r28,*r2A,*r2C,*r2E
		,*r30,*r32,*r34,*r36,*r38,*r3A,*r3C,*r3E, *bminusc_ptr;
  #else		// Same list of ptrs as above (except for base-ptr), but now make them static:
	static uint32  *k1_arr, *k2_arr;
	static vec_dbl *isrt2,*sqrt2,*one,*two, *cc0,*ss0,*cc1,*ss1,*cc3,*ss3
		,*forth,*tmp0,*tmp1,*tmp2,*tmp3,*tmp4,*tmp5,*tmp6,*tmp7
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*r00,*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E
		,*r10,*r12,*r14,*r16,*r18,*r1A,*r1C,*r1E
		,*r20,*r22,*r24,*r26,*r28,*r2A,*r2C,*r2E
		,*r30,*r32,*r34,*r36,*r38,*r3A,*r3C,*r3E, *bminusc_ptr;
  #endif

#endif

/*...If a new runlength or first-pass radix, it is assumed this function has been first-called with init_sse2 = true to
     initialize static data and lcoal-storage prior to actual use in computing a transform-based result.
*/
	// For case of 2-input modmul, pointer to 2nd (presumed already fwd-FFTed) input array is supplied via fwd_fft_only argument:
	/* v20: Bits 2:3 of fwd_fft = 3 means "dyadic-multiply of 2 inputs a and b, with both already forward-FFTed:
		(double *)a = FFT(a), (double *)(fwd_fft_only - mode_flag) = FFT(b).
	In this case we require bits 0:1 == 0, and fwd_fft = & ~0xC yields pointer to FFT(b), and we skip over fwd-FFT directly to
	the dyadic-multiply FFT(a) * FFT(b) step, then iFFT the product, storing the result in a[].
	*/
	double *b = 0x0;
	if(fwd_fft_only >> 2) {		// Do the following 3 steps in both cases - if bits 2:3 == 0 the ANDs are no-ops...
		// The submul-auxiliary array c_arr[], if present, in already in proper double[] form, but b[] needs low-bits-cleared and casting
		b = (double *)(fwd_fft_only & ~0xCull);
		// BUT, if bits 2:3 == 0, must avoid zeroing fwd_fft_only since "do 2-input dyadic-mul following fwd-FFT" relies on that != 0:
		if(fwd_fft_only & 0xC) {
			ASSERT((fwd_fft_only & 0xF) == 0xC,"Illegal value for bits 2:3 of fwd_fft_only!");	// Otherwise bits 2:3 should've been zeroed prior to entry
			fwd_fft_only = 3ull;
		}
	}

	/* In order to eschew complex thread-block-and-sync logic related to the local-store-init step in multithread mode,
	switch to a special-init-mode-call paradigm, in which the function is inited once (for as many threads as needed)
	prior to being executed:
	*/
/**************************************************************************************************************************************/
/*** To-Do: Need to add code to allow for re-init when any of the FFT-related params or #threads changes during course of execution ***/
/**************************************************************************************************************************************/

	/* Here this variable is somewhat misnamed because it is used to init both non-SIMD and SIMD-specific data */
	if(init_sse2)	// Just check nonzero here, to allow the *value* of init_sse2 to store #threads
	{
		nsave = n;
		if(init_sse2 > max_threads)	// current SIMD local-alloc insufficient
		{
			ASSERT(thr_id == -1, "Init-mode call must be outside of any multithreading!");
			max_threads = init_sse2;
		#ifndef COMPILER_TYPE_GCC
			ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
		#endif
		//	printf("%Ns: max_threads = %d, NTHREADS = %d\n",func, max_threads, NTHREADS);

		#ifdef USE_SSE2
		if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)sm_arr);	sm_arr=0x0;
			free((void *)sc_arr);	sc_arr=0x0;
		}
		// Index vectors used in SIMD roots-computation.
		// The AVX512 compute-sincos-mults code needs 2 elements per complex-double-load, so use 14*RE_IM_STRIDE per array
		// to alloc storage here for all cases, even though that leaves upper array halves unused for sub-AVX512.
		sm_arr = ALLOC_INT(sm_arr, max_threads*28*RE_IM_STRIDE + 16);	if(!sm_arr){ sprintf(cbuf, "ERROR: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sm_ptr = ALIGN_INT(sm_arr);
		ASSERT(((uintptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");
		// Twiddles-array: Need 0x92 slots for data, plus need to leave room to pad-align.
		// v20: To support inline a*(b-c) for p-1 stage 2, need 2*RADIX = 64 added vec_dbl, thus 0x98 ==> 0xd8:
		sc_arr = ALLOC_VEC_DBL(sc_arr, 0xd8*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

		/* Use low 64 vec_dbl slots of sc_arr for temporaries, next 8 for scratch, next 7 for the nontrivial complex 16th roots,
		next 62 for the doubled sincos twiddles, next 4 for [1.0,2.0,0.25,sqrt2] and at least 3 more to allow for 64-byte alignment of the array.

		*** NOTE ***: Offsets below must match those in radix32_dyadic_square,
					since routines share DFT macros which use many literal byte offsets to reduce argument count.
		*/
		  #ifdef MULTITHREAD
			__i0  = sm_ptr;
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x48;
			cc0	  = sc_ptr + 0x49;
			ss0	  = sc_ptr + 0x4a;
			cc1	  = sc_ptr + 0x4b;
			ss1	  = sc_ptr + 0x4c;
			cc3	  = sc_ptr + 0x4d;
			ss3	  = sc_ptr + 0x4e;
			two   = sc_ptr + 0x8f;	// PAIR_SQUARE_4_AVX2 not used for Fermat-mod (i.e. no need for *forth),
			forth = sc_ptr + 0x90;	// but note SSE2_RADIX32_WRAPPER_DIF/T asm macros are shared by Mersenne & Fermat,
			sqrt2 = sc_ptr + 0x91;	// and assume two = isrt2 + 0x47, sqrt2 = two + 2, one = two + 3, i.e. need to
			one	  = sc_ptr + 0x92;	// "leave slot for forth just above two".
			bminusc_ptr = sc_ptr + 0x98;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
				VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );	VEC_DBL_INIT(forth, 0.25 );
				VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
				VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
				VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
				/* Move on to next thread's local store */
				isrt2 += 0xd8;
				cc0   += 0xd8;
				ss0   += 0xd8;
				cc1   += 0xd8;
				ss1   += 0xd8;
				cc3   += 0xd8;
				ss3   += 0xd8;
				one   += 0xd8;
				two   += 0xd8;
				forth += 0xd8;
				sqrt2 += 0xd8;
				bminusc_ptr = sc_ptr + 0xd8;
			}
		  #else
			k1_arr = sm_ptr;
			k2_arr = sm_ptr + 14*RE_IM_STRIDE;
			r00		= sc_ptr + 0x00;	isrt2	= sc_ptr + 0x48;
			r02		= sc_ptr + 0x02;	cc0		= isrt2 + 0x01;	ss0   = cc0 + 0x01;
			r04		= sc_ptr + 0x04;	cc1		= cc0 + 0x02;	ss1   = cc0 + 0x03;
			r06		= sc_ptr + 0x06;	cc3		= cc0 + 0x04;	ss3   = cc0 + 0x05;
			r08		= sc_ptr + 0x08;	c00		= cc0 + 0x06;
			r0A		= sc_ptr + 0x0a;	c10		= cc0 + 0x08;
			r0C		= sc_ptr + 0x0c;	c08		= cc0 + 0x0a;
			r0E		= sc_ptr + 0x0e;	c18		= cc0 + 0x0c;
			r10		= sc_ptr + 0x10;	c04		= cc0 + 0x0e;
			r12		= sc_ptr + 0x12;	c14		= cc0 + 0x10;
			r14		= sc_ptr + 0x14;	c0C		= cc0 + 0x12;
			r16		= sc_ptr + 0x16;	c1C		= cc0 + 0x14;
			r18		= sc_ptr + 0x18;	c02		= cc0 + 0x16;
			r1A		= sc_ptr + 0x1a;	c12		= cc0 + 0x18;
			r1C		= sc_ptr + 0x1c;	c0A		= cc0 + 0x1a;
			r1E		= sc_ptr + 0x1e;	c1A		= cc0 + 0x1c;
			r20		= sc_ptr + 0x20;	c06		= cc0 + 0x1e;
			r22		= sc_ptr + 0x22;	c16		= cc0 + 0x20;
			r24		= sc_ptr + 0x24;	c0E		= cc0 + 0x22;
			r26		= sc_ptr + 0x26;	c1E		= cc0 + 0x24;
			r28		= sc_ptr + 0x28;	c01		= cc0 + 0x26;
			r2A		= sc_ptr + 0x2a;	c11		= cc0 + 0x28;
			r2C		= sc_ptr + 0x2c;	c09		= cc0 + 0x2a;
			r2E		= sc_ptr + 0x2e;	c19		= cc0 + 0x2c;
			r30		= sc_ptr + 0x30;	c05		= cc0 + 0x2e;
			r32		= sc_ptr + 0x32;	c15		= cc0 + 0x30;
			r34		= sc_ptr + 0x34;	c0D		= cc0 + 0x32;
			r36		= sc_ptr + 0x36;	c1D		= cc0 + 0x34;
			r38		= sc_ptr + 0x38;	c03		= cc0 + 0x36;
			r3A		= sc_ptr + 0x3a;	c13		= cc0 + 0x38;
			r3C		= sc_ptr + 0x3c;	c0B		= cc0 + 0x3a;
			r3E		= sc_ptr + 0x3e;	c1B		= cc0 + 0x3c;
			tmp0	= sc_ptr + 0x40;	c07		= cc0 + 0x3e;
			tmp1	= sc_ptr + 0x41;	c17		= cc0 + 0x40;
			tmp2	= sc_ptr + 0x42;	c0F		= cc0 + 0x42;
			tmp3	= sc_ptr + 0x43;	c1F		= cc0 + 0x44;
			tmp4	= sc_ptr + 0x44;	two  	= cc0 + 0x46;
			tmp5	= sc_ptr + 0x45;	forth	= cc0 + 0x47;
			tmp6	= sc_ptr + 0x46;	sqrt2	= cc0 + 0x48;
			tmp7	= sc_ptr + 0x47;	one	 	= cc0 + 0x49;
			bminusc_ptr = sc_ptr + 0x98;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
			VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );	VEC_DBL_INIT(forth, 0.25 );
			VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
			VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
			VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
		  #endif
		#endif	// USE_SSE2
		}	// SIMD alloc block

		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Initialize a scratch array containing N2/16 indices - again use big
		!   as-yet-unused A-array for this, but making sure the section of A used
		!   for the itmp space and that sent to the bit_reverse_int for scratch space
		!   don't overlap:
		*/
		ASSERT(N2 == n/2, "N2 bad!");
		itmp = (int *)&arr_scratch[N2/32];	/* Conservatively assume an int might be as long as 8 bytes here */
		for(i=0; i < N2/32; i++)
		{
			itmp[i]=i;
		}

		/*...then bit-reverse INDEX with respect to N/16. Radices are sent to BR routine
		in opposite order from which they are processed in forward FFT.
		*/
		bit_reverse_int(itmp, N2/32, nradices_prim-5, &radix_prim[nradices_prim-6], -1, (int*)arr_scratch);
	/*
	!...trig array for wrapper/squaring part is here. We will make use of the fact that the sincos
	!   data will be in bit-reversed order, together with the symmetries of the sin and cosine
	!   functions about pi/4 and pi/2, to cut our table length by a factor of 16 here. Here's how:
	!   the ordered sincos data (prior to bit-reversal of the data and reduction of table length) are a set of
	!   complex exponentials exp(I*j*2*pi/N), where j=0,...,N/2-1 and N is the complex FFT length.
	!   That is, the argument (complex polar angle) of the exponential is in [0,pi). Within the
	!   fused radix32/wrapper-square/radix32inv pass routine, we process 16 such complex exponentials
	!   at a time, together with 32 complex FFT array data. Since these 16 sincos data are bit-reversal-
	!   reordered, we can treat them as eight pairs, each of the form (angle, angle+pi/2):
	!
	!   x0+(0,pi/2), x1+(0,pi/2), x2+(0,pi/2), x3+(0,pi/2), x4+(0,pi/2), x5+(0,pi/2), x6+(0,pi/2), x7+(0,pi/2).
	!
	!   Given exp(I*angle) = (re,im), we use the complex symmetry  exp(I*(angle+pi/2)) = (-im,re)
	!   to get the second exponential of each pair. Since, for given x0, the other sincos arguments have form
	!
	!   x1=x0+pi/4, x2=x0+pi/8, x3=x0+3*pi/8, x4=x0+pi/16, x5=x0+5*pi/16, x=x0+3*pi/16, x=x0+7*pi/16
	!
	!   (i.e. x_{0,1,2,3,4,5,6,7} = x0 + {0,4,2,6,1,5,3,7)}*pi/16), we further let
	!
	!      exp(I*  pi/8 ) = (cos(  pi/4 ) , sin(  pi/4 )) = (1/sqrt(2), 1/sqrt(2)):= (ISRT2, ISRT2);
	!
	!      exp(I*  pi/8 ) = (cos(  pi/8 ) , sin(  pi/8 )):= (c    ,s    );
	!      exp(I*3*pi/8 ) = (cos(3*pi/8 ) , sin(3*pi/8 )) = (s    ,c    );
	!
	!      exp(I*  pi/16) = (cos(  pi/16) , sin(  pi/16)):= (c32_1,s32_1);
	!      exp(I*5*pi/16) = (cos(5*pi/16) , sin(5*pi/16)):= (s32_3,c32_3);
	!      exp(I*3*pi/16) = (cos(3*pi/16) , sin(3*pi/16)):= (c32_3,s32_3);
	!      exp(I*7*pi/16) = (cos(7*pi/16) , sin(7*pi/16)):= (s32_1,c32_1) .
	!
	!   Given only exp(I*x0) = (re,im), this allows us to get the remaining seven complex exponentials via
	!
	!      exp(I*(x0       )) = ( re, im)
	!      exp(I*(x0+  pi/2)) = I*exp(I*x0) = (-im, re)
	!
	!      exp(I*(x1       )) = (re-im, re+im)*ISRT2
	!      exp(I*(x1+  pi/2)) = I*exp(I*x1)
	!
	!      exp(I*(x2       )) = ( re, im)*(c    ,s    ) = (re*c    -im*s    , re*s    +im*c    )
	!      exp(I*(x2+  pi/2)) = I*exp(I*x2)
	!
	!      exp(I*(x3       )) = ( re, im)*(s    ,c    ) = (re*s    -im*c    , re*c    +im*s    )
	!      exp(I*(x3+  pi/2)) = I*exp(I*x3),
	!
	!      exp(I*(x4       )) = ( re, im)*(c32_1,s32_1) = (re*c32_1-im*s32_1, re*s32_1+im*c32_1)
	!      exp(I*(x4+  pi/2)) = I*exp(I*x4)
	!
	!      exp(I*(x5       )) = ( re, im)*(s32_3,c32_3) = (re*s32_3-im*c32_3, re*c32_3+im*s32_3)
	!      exp(I*(x5+  pi/2)) = I*exp(I*x5),
	!
	!      exp(I*(x6       )) = ( re, im)*(c32_3,s32_3) = (re*c32_3-im*s32_3, re*s32_3+im*c32_3)
	!      exp(I*(x6+  pi/2)) = I*exp(I*x6)
	!
	!      exp(I*(x7       )) = ( re, im)*(s32_1,c32_1) = (re*s32_1-im*c32_1, re*c32_1+im*s32_1)
	!      exp(I*(x7+  pi/2)) = I*exp(I*x7),
	!
	!   where (a,b) = a+I*b and (a,b)*(c,d) denotes a complex product, i.e. using both explicit-I* and complex-
	!   as-real-pair notation, (a,b)*(c,d) = (a+I*b)*(c+I*d) = (a*c - b*d) + I*(a*d + b*c) = (a*c - b*d, a*d + b*c).
	!
	!   We already have the constants ISRT2, c, s, c32_1, s32_1, c32_3, s32_3 in registers for the radix-32 transform.
	!   By saving the 12 scalar products re*c, im*s, re*s, im*c, re*c32_1, im*s32_1, re*s32_1, im*c32_1, re*c32_3, im*s32_3, re*s32_3, im*c32_3,
	!   we can get each of exp(I*x2|x3), exp(I*x4|x5), exp(I*x6|x7) using just 4 FMUL and 4 FADD, so our total
	!   cost for the 16 complex sincos data is 2 loads from memory, 14 FMUL, 14 FADD and 8 negations.
	!
	!   We further use the fact that for every two blocks of 32 FFT data processed together, the wrapper sincos
	!   data for the second (upper) block are just the reflection about pi/2 of those for the first block.
	!   This cuts the table length in half again.
	*/

	/*
	!...Allocate and initialize an index array containing N/32 indices to store the sequentially-rearranged
	!   FFT sincos data indices. (Only need N/32 of these since we need one base sincos index for every 32 complex data).
	!   We don't need a separate sincos array for the real/complex wrapper phase, since this uses the same sincos datum
	!   as is used for the the first of each of the two blocks of 32 complex FFT data.
	*/
		if(index_ptmp != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)index_ptmp);	index_ptmp=0x0;
		}
		index_ptmp = ALLOC_INT(index_ptmp, N2/32);
		ASSERT(index_ptmp != 0,"ERROR: unable to allocate array INDEX!");
		index      = ALIGN_INT(index_ptmp);
	/*
	!...Now rearrange FFT sincos indices using the main loop structure as a template.
	!   The first length-2 block is a little different, since there we process the 0-31 and 32-63 array
	!   elements separately.
	*/
		index[0]=itmp [0];
		index[1]=itmp [1];

		k1 =32;	/* init uncompressed wrapper sincos array index */
		k  =2;	/* init   compressed wrapper sincos array index */

		blocklen=32;
		blocklen_sum=32;
		j2_start=192;

		for(i = nradices_prim-7; i >= 0; i--)	/* Radices get processed in reverse order here as in forward FFT. */
		{
		  kp = k1 + blocklen;

		  for(m = 0; m < (blocklen-1)>>1; m +=16)	/* Since now process TWO 32-element sets per loop execution, only execute the loop half as many times as before. */
		  {
	/*...grab the next 2 FFT sincos indices: one for the lower (j1, in the actual loop) data block, one for the upper (j2) data block... */

			index[k  ]=itmp[(k1   )>>4];
			index[k+1]=itmp[(kp-16)>>4];

			k1 = k1+16;
			kp = kp-16;

			k  = k + 2;
		  }

		  k1 = k1 + (blocklen >> 1);

		  if(j2_start == n-64)break;

		  blocklen_sum = blocklen_sum + blocklen;
		  ASSERT(i != 0,"ERROR 10!");
		  blocklen = (radix_prim[i-1]-1)*blocklen_sum;

		  j2_start = j2_start+(blocklen<<2);
		}

		j1 = 0;

		/* Restore zeros here, to prevent any barfing due to interpretation of the above integer values as floats,
		in case at some future point we find it useful to be able to re-use a part of the main a-array for scratch: */
		for(i=0; i < N2/32; i++)
		{
			itmp[i]=0;
		}

		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
#ifdef MULTITHREAD
	ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
  #ifdef USE_SSE2
	k1_arr = __i0 + thr_id*28*RE_IM_STRIDE;
	k2_arr = k1_arr +      14*RE_IM_STRIDE;
	r00 = __r0 + thr_id*0xd8;isrt2	= r00 + 0x48;
	r02		= r00 + 0x02;	cc0		= r00 + 0x49;
	r04		= r00 + 0x04;	cc1		= r00 + 0x4b;
	r06		= r00 + 0x06;	cc3		= r00 + 0x4d;
	r08		= r00 + 0x08;	c00		= r00 + 0x4f;
	r0A		= r00 + 0x0a;	c10		= r00 + 0x51;
	r0C		= r00 + 0x0c;	c08		= r00 + 0x53;
	r0E		= r00 + 0x0e;	c18		= r00 + 0x55;
	r10		= r00 + 0x10;	c04		= r00 + 0x57;
	r12		= r00 + 0x12;	c14		= r00 + 0x59;
	r14		= r00 + 0x14;	c0C		= r00 + 0x5b;
	r16		= r00 + 0x16;	c1C		= r00 + 0x5d;
	r18		= r00 + 0x18;	c02		= r00 + 0x5f;
	r1A		= r00 + 0x1a;	c12		= r00 + 0x61;
	r1C		= r00 + 0x1c;	c0A		= r00 + 0x63;
	r1E		= r00 + 0x1e;	c1A		= r00 + 0x65;
	r20		= r00 + 0x20;	c06		= r00 + 0x67;
	r22		= r00 + 0x22;	c16		= r00 + 0x69;
	r24		= r00 + 0x24;	c0E		= r00 + 0x6b;
	r26		= r00 + 0x26;	c1E		= r00 + 0x6d;
	r28		= r00 + 0x28;	c01		= r00 + 0x6f;
	r2A		= r00 + 0x2a;	c11		= r00 + 0x71;
	r2C		= r00 + 0x2c;	c09		= r00 + 0x73;
	r2E		= r00 + 0x2e;	c19		= r00 + 0x75;
	r30		= r00 + 0x30;	c05		= r00 + 0x77;
	r32		= r00 + 0x32;	c15		= r00 + 0x79;
	r34		= r00 + 0x34;	c0D		= r00 + 0x7b;
	r36		= r00 + 0x36;	c1D		= r00 + 0x7d;
	r38		= r00 + 0x38;	c03		= r00 + 0x7f;
	r3A		= r00 + 0x3a;	c13		= r00 + 0x81;
	r3C		= r00 + 0x3c;	c0B		= r00 + 0x83;
	r3E		= r00 + 0x3e;	c1B		= r00 + 0x85;
	tmp0	= r00 + 0x40;	c07		= r00 + 0x87;
	tmp1	= r00 + 0x41;	c17		= r00 + 0x89;
	tmp2	= r00 + 0x42;	c0F		= r00 + 0x8b;
	tmp3	= r00 + 0x43;	c1F		= r00 + 0x8d;
	tmp4	= r00 + 0x44;	two  	= r00 + 0x8f;
	tmp5	= r00 + 0x45;	forth	= r00 + 0x90;
	tmp6	= r00 + 0x46;	sqrt2	= r00 + 0x91;
	tmp7	= r00 + 0x47;	one	 	= r00 + 0x92;
	bminusc_ptr = r00 + 0x98;
  #endif
#endif
	/*...If a new runlength, should not get to this point: */
	ASSERT(n == nsave,"n != nsave");

/*
!   SOLVING THE CACHE FLOW PROBLEM FOR BIT-REVERSED ARRAY DATA:
!
!   In a previous version of this routine, the data entering the wrapper/square
!   step were ordered, and (j, N-j) -indexed data pairs were processed
!   together with data from a simple ordered sincos array W0, which
!   stored W0(j) = [1+exp(2*I*pi*j/N)], j=0,...,N/2-1, where n is the
!   complex vector length. Life was easy in the wrapper/square step, but
!   the explicit bit-reversals needed to arrange data for the DIT forward
!   transform (and after the wrapper/square for the DIT inverse transform)
!   had a significant adverse effect on performance of the code overall.
!
!   In this version, we do a DIF forward transform, so data enter the wrapper/square
!   step in bit-reversed order. Lacking a simpler pattern, I at first processed
!   them in the same order as in the old version, i.e. simply fetched
!   data in slots (index(j),index(N-j)) and processed as before. However,
!   if one looks at the actual access patterns that result from this non-strategy,
!   things are ugly: whereas j and N-j have nice predictable behavior,
!   index(j) and index(N-j) bounce all over the place. In fact, it is hard
!   to devise a worse stategy than this, from a caching perspective.
!
!   So I thought, why not loop over index(j), rather than j?
!   This requires us to bit-reverse the data in the W0 (sincos) array
!   but that imposes no penalty since these data are precomputed anyway.
!   Thus, at least the index of the first data element of each pair
!   behaves in a predictable fashion.
!
!   Well, when one does things this way, it turns out the second index also
!   behaves in a predictable way, although not quite as simply as for ordered
!   data. We note that under a bit reversal, the index range j=0,...,N/2-1
!   maps to the evens and j=N/2,...,N-1 maps to the odds, respectively.
!   This tells us we should loop over index(j) and increment by two.
!   When we do this, a surprisingly nice pattern emerges for index(N-j).
!   Even more surprisingly, this pattern also persists for N an odd prime
!   times a power of 2, with a tiny amount of added bookeeping required.
!
!   What happens is that as the first (even) index advances monotonically
!   from 0 to N-2, the second (odd) index runs through O(log2(N)) contiguous
!   odd-index subranges, but runs through each one in REVERSE order, i.e.
!   in decrements of two. The precise number of these odd-index subranges
!   (which overlap a corresponding even-index subrange to form what I call
!   a BLOCK) is equal to log2(power-of-2 part of N) + 1, where
!
!	N = P*2^M, where P is a small prime first radix > 1 (possibly 2.)
!
!   Here are four examples, for N = 16, 24, 20 and 28, assuming a 1-D, unpadded
!   zero-offset data array. I've broken chunks of data into blocks that reveal the pattern.
!   Note the correlation between the bit pattern of N-1 (written vertically, with least
!   significant bit at the top) and the pattern of block lengths
!  (block length = half the total number of complex data processed):
!
!       N=16 (FIRST RADIX = 2):                 N=24 (FIRST RADIX = 3):                 N=20 (FIRST RADIX = 5):                 N=28 (FIRST RADIX = 7):
!
!       N - 1 = 15 = 1111_2:                    N - 1 = 23 = 10111_2:                   N - 1 = 19 = 10011_2:                   N - 1 = 27 = 11011_2:
!       j1 =    j2 =            Current bits    j1 =    j2 =            Current bits    j1 =    j2 =            Current bits    j1 =    j2 =            Current bits
!       index(j)index(N-j)      of (N - 1)_2:   index(j)index(N-j)      of (N - 1)_2:   index(j)index(N-j)      of (N - 1)_2:   index(j)index(N-j)      of (N - 1)_2:
!       ------  ------          -------------   ------  ------          -------------   ------  ------          -------------   ------  ------          -------------
!       0,1 (j=0,N/2) done separately   1       0,1 (j=0,N/2) done separately   1       0,1 (j=0,N/2) done separately   1       0,1 (j=0,N/2) done separately   1
!       ------  ------                          ------  ------                          ------  ------                          ------  ------
!       2       3       length-2 block  1       2       3       length-2 block  1       2       3       length-2 block  1       2       3       length-2 block  1
!       ------  ------                          ------  ------	                        ------  ------                          ------  ------
!       4       7       length-4 block  1       4       7       length-4 block  1       4       19                              4       27
!       6       5                               6       5                               6       17                              6       25
!       ------  ------                          ------  ------                          8       15                              8       23
!       8       15                              8       23                              10      13      length-16 block 100     10      21
!       10      13      length-8 block  1       10      21                              12      11                              12      19
!       12      11                              12      19                              14      9                               14      17      length-24 block 110
!       14      9                               14      17      length-16 block 10      16      7                               16      15
!       ------  ------                          16      15                              18      5                               18      13
!                                               18      13                              ------  ------                          20      11
!                                               20      11                                                                      22      9
!                                               22      9                                                                       24      7
!                                               ------  ------                                                                  26      5
!   So the formula for calculating block sizes is as follows: At the end of each block, we
!   multiply the accumulated blocklength (which is just double the current blocklength) by
!   the next value of (current bits), which is = (radix_next - 1), to get the length
!   of the next block. To get the value of j2_start for the next block,
!   we add the length of the next block to j2_start for the current block.
!
!...Exploiting the pattern makes for a much cache-friendlier code, which follows. Note
!   that if we are processing ordered modular data alongside bit-reversed floating data,
!   we actually have two index pairs to keep track of, e.g. for N=16:
!
!   Floating index 1 (j1):  2, 4, 6, 8,10,12,14
!   Floating index 2 (j2):  3, 7, 5,15,13,11, 9	(0 and 1 done separately)
!
!   Modular  index 1 (j3):  1, 2, 3, 4, 5, 6, 7
!   Modular  index 2 (j4): 15,14,13,12,11,10, 9	(0 and N/2 done separately)
!
!   so the floating index 1 is just double its modular counterpart. The index 2's must be
!   stored separately.
!
!
!   FUSING THE WRAPPER-SQUARE STEP WITH THE RADIX-16 PASSES PRECEDING AND FOLLOWING IT:
!
!   OK, so you've made it this far, and justifiably are hoping to see some actual, honest-
!   to-goodness code. Alas, things must get just a bit more complicated before that happy
!   moment arrives. That's because in order to squeeze maximum performance out of our code,
!   we should seize any opportunities to minimize unnecessary data movement. One such
!   presents itself here. We have a radix-16 DIF pass immediately preceding and a radix-16
!   DIT pass immediately following the wrapper/square, so we do all three of the above steps
!   in a single in-place sweep, thus reducing the number of passes through the big data
!   arrays from three to one. We do the first 16 complex data (and their images under the (j,N'-j)
!   correlation) separately, after which blocklengths are always a multiple of 16, allowing us to simplify the indexing.
!
!   The general loop structure looks as follows:
!
!   loop:
!   	do m = 1,blocklen
!
!   	  call pair_square(a(j1),a(j1+1),a(j2),a(j2+1),w(1,j),w(2,j))
!
!   	  j  = j +1
!   	  j1 = j1+4
!   	  j2 = j2-4
!   	enddo
!   	...
!   	blocklen = 2*blocklen
!   	j2_start = j2_start+ishft(blocklen,2)
!   	j2=j2_start
!   	cycle
!   endloop
!
!   So the first 16 complex elements and their images are processed as follows:
!
!   exe0: process a(0,1) and a(2,3) separately;
!
!   exe1:
!   	j=1
!   	blocklen = 1
!   	j1=4
!   	j2_start=6
!
!   	  call pair_square(a( 4),a( 5),a( 6),a( 7),w(:,1))
!
!   	  j  = 2
!   	  j1 = 8
!   	  j2 = 2
!   	enddo
!
!   	blocklen = 2
!   	j2_start = 14
!   end1
!
!   exe2:
!   	j=2
!   	blocklen = 2
!   	j1=8
!   	j2_start=14
!
!   	  call pair_square(a( 8),a( 9),a(14),a(15),w(:,2))
!   	  call pair_square(a(12),a(13),a(10),a(11),w(:,3))
!
!   	  j  = 4
!   	  j1 = 16
!   	  j2 = 6
!   	enddo
!
!   	blocklen = 4
!   	j2_start = 30
!   end2
!
!   exe3:
!   	j=4
!   	blocklen = 4
!   	j1=16
!   	j2_start=30
!
!   	  call pair_square(a(16),a(17),a(30),a(31),w(:,4))
!   	  call pair_square(a(20),a(21),a(26),a(27),w(:,5))
!   	  call pair_square(a(24),a(25),a(22),a(23),w(:,6))
!   	  call pair_square(a(28),a(29),a(18),a(19),w(:,7))
!
!   	  j  = 8
!   	  j1 = 32
!   	  j2 = 14
!   	enddo
!
!   	blocklen = 8
!   	j2_start = 62
!   end3
!
!   exe4:
!   	j=8
!   	blocklen = 8
!   	j1=32
!   	j2_start=62
!
!   	  call pair_square(a(32),a(33),a(62),a(63),w(:,8))
!   	  call pair_square(a(36),a(37),a(58),a(59),w(:,9))
!   	  call pair_square(a(40),a(41),a(54),a(55),w(:,10))
!   	  call pair_square(a(44),a(45),a(50),a(51),w(:,11))
!   	  call pair_square(a(48),a(49),a(46),a(47),w(:,12))
!   	  call pair_square(a(52),a(53),a(42),a(43),w(:,13))
!   	  call pair_square(a(56),a(57),a(38),a(39),w(:,14))
!   	  call pair_square(a(60),a(61),a(34),a(35),w(:,15))
!
!   	  j  = 16
!   	  j1 = 64
!   	  j2 = 30
!   	enddo
!
!   	blocklen = 16
!   	j2_start = 126
!   end4
!
!   exe5:
!   	j=16
!   	blocklen = 16
!   	j1=64
!   	j2_start=126
!
!   	  call pair_square(a( 64),a( 65),a(126),a(127),w(:,16))
!   	  call pair_square(a( 68),a( 69),a(122),a(123),w(:,17))
!   	  call pair_square(a( 72),a( 73),a(118),a(119),w(:,18))
!   	  call pair_square(a( 76),a( 77),a(114),a(115),w(:,19))
!   	  call pair_square(a( 80),a( 81),a(110),a(111),w(:,20))
!   	  call pair_square(a( 84),a( 85),a(106),a(107),w(:,21))
!   	  call pair_square(a( 88),a( 89),a(102),a(103),w(:,22))
!   	  call pair_square(a( 92),a( 93),a( 98),a( 99),w(:,23))
!   	  call pair_square(a( 96),a( 97),a( 94),a( 95),w(:,24))
!   	  call pair_square(a(100),a(101),a( 90),a( 91),w(:,25))
!   	  call pair_square(a(104),a(105),a( 86),a( 87),w(:,26))
!   	  call pair_square(a(108),a(109),a( 82),a( 83),w(:,27))
!   	  call pair_square(a(112),a(113),a( 78),a( 79),w(:,28))
!   	  call pair_square(a(116),a(117),a( 74),a( 75),w(:,29))
!   	  call pair_square(a(120),a(121),a( 70),a( 71),w(:,30))
!   	  call pair_square(a(124),a(125),a( 66),a( 67),w(:,31))
!
!   	  j  = 32
!   	  j1 = 128
!   	  j2 = 62
!   	enddo
!
!   	blocklen = 32
!   	j2_start = 254
!   end5
!
!   From here onward the blocklength will always be an integer multiple of 32, i.e. we can process each block using pairs of nonoverlapping
!   blocks of 32 complex data each, which is compatible to fusion with radix-32 pass routines. For example, the next (fifth) loop execution
!   looks like:
!
!   exe6:
!   	j=32
!   	blocklen = 32
!   	j1=128
!   	j2_start=254
!
!   	  call pair_square(a(128),a(129),a(254),a(255),w(:,32))
!   	  call pair_square(a(132),a(133),a(250),a(251),w(:,33))
!   	  call pair_square(a(136),a(137),a(246),a(247),w(:,34))
!   	  call pair_square(a(140),a(141),a(242),a(243),w(:,35))
!   	  call pair_square(a(144),a(145),a(238),a(239),w(:,36))
!   	  call pair_square(a(148),a(149),a(234),a(235),w(:,37))
!   	  call pair_square(a(152),a(153),a(230),a(231),w(:,38))
!   	  call pair_square(a(156),a(157),a(226),a(227),w(:,39))
!   	  call pair_square(a(160),a(161),a(222),a(223),w(:,40))
!   	  call pair_square(a(164),a(165),a(218),a(219),w(:,41))
!   	  call pair_square(a(168),a(169),a(214),a(215),w(:,42))
!   	  call pair_square(a(172),a(173),a(210),a(211),w(:,43))
!   	  call pair_square(a(176),a(177),a(206),a(207),w(:,44))
!   	  call pair_square(a(180),a(181),a(202),a(203),w(:,45))
!   	  call pair_square(a(184),a(185),a(198),a(199),w(:,46))
!   	  call pair_square(a(188),a(189),a(194),a(195),w(:,47))
!
!   	  call pair_square(a(192),a(193),a(190),a(191),w(:,48))
!   	  call pair_square(a(196),a(197),a(186),a(187),w(:,49))
!   	  call pair_square(a(200),a(201),a(182),a(183),w(:,50))
!   	  call pair_square(a(204),a(205),a(178),a(179),w(:,51))
!   	  call pair_square(a(208),a(209),a(174),a(175),w(:,52))
!   	  call pair_square(a(212),a(213),a(170),a(171),w(:,53))
!   	  call pair_square(a(216),a(217),a(166),a(167),w(:,54))
!   	  call pair_square(a(220),a(221),a(162),a(163),w(:,55))
!   	  call pair_square(a(224),a(225),a(158),a(159),w(:,56))
!   	  call pair_square(a(228),a(229),a(154),a(155),w(:,57))
!   	  call pair_square(a(232),a(233),a(150),a(151),w(:,58))
!   	  call pair_square(a(236),a(237),a(146),a(147),w(:,59))
!   	  call pair_square(a(240),a(241),a(142),a(143),w(:,60))
!   	  call pair_square(a(244),a(245),a(138),a(139),w(:,61))
!   	  call pair_square(a(248),a(249),a(134),a(135),w(:,62))
!   	  call pair_square(a(252),a(253),a(130),a(131),w(:,63))
!
!   	...and these 64 complex data can be processed as follows:
!	1) do radix-32 DIF (using sincos data  0-31) to get a(128:191);
!	2) do radix-32 DIF (using sincos data 32-63) to get a(192:255);
!	3) combine the 64 resulting complex array data as shown above;
!	4) do radix-32 DIT on a(128:191);
!	5) do radix-32 DIT on a(192:255).
!   	The radix-32 passes share many register data, providing added savings.
!   	  j  = 64
!   	  j1 = 256
!   	  j2 = 126
!   	enddo
!
!   	blocklen = 64
!   	j2_start = 510
!   end6
*/
/* Init the loop-control variables: */

	i            = ws_i           ;
	j1           = ws_j1          ;
	j2           = ws_j2          ;
	j2_start     = ws_j2_start    ;
	k            = ws_k           ;
	m            = ws_m           ;
	blocklen     = ws_blocklen    ;
	blocklen_sum = ws_blocklen_sum;

//	fprintf(stderr,"%s: stride = %d\n",func,stride);
//	fprintf(stderr,"%s: On entry: j1,j2 = %u, %u, nradices_prim = %u, blocklen = %u\n",func,j1,j2,nradices_prim,blocklen);

	/* If j1 == 0 we need to init the loop counters; otherwise, just jump
	   right in and pick up where we left off on the previous pair of blocks:
	*/
	if(j1 > 0) {
	//	fprintf(stderr,"%s: Jumping into loop!\n",func);
		goto jump_in;
	}

/*
!...All but the first two radix-16 blocks are done on Mr. Henry Ford's assembly line. From there onward the blocklength
!   will always be an integer multiple of 32, i.e. we can process each block using pairs of nonoverlapping blocks of 32
!   complex data each, which is compatible to fusion with radix-32 pass routines.
*/

for(i = nradices_prim-6; i >= 0; i-- )	/* Main loop: lower bound = nradices_prim-radix_now. */
{						/* Remember, radices get processed in reverse order here as in forward FFT. */

#ifdef USE_AVX512
	for(m = 0; m < (blocklen-1)>>1; m += 64) /* In AVX-512, process eight 32-complex-double datasets per loop execution, thus only execute the loop half as many times as for AVX case. */
#elif defined(USE_AVX)
	for(m = 0; m < (blocklen-1)>>1; m += 32) /* In AVX mode, process four 32-complex-double datasets per loop execution, thus only execute the loop half as many times as for scalar/SSE2 case. */
#else
	for(m = 0; m < (blocklen-1)>>1; m += 16) /* Both scalar and SSE2 modes process two 32-complex-double datasets per loop, only execute the loop half as many times as before. */
#endif
	{
		// This tells us when we've reached the end of the current data block:
		// Apr 2014: Must store intermediate product j1*radix0 in a 64-bit int to prevent overflow!
		if(j1 && ((uint64)j1*radix0)%n == 0)
		{
		//	fprintf(stderr,"(%s: j1 && j1*radix0 == 0 (mod n)) check hit: returning\n",func);
			return;
		}

jump_in:	/* Entry point for all blocks but the first. */

	  j1pad = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* floating padded-array 1st element index is here */
	  j2pad = j2 + ( (j2 >> DAT_BITS) << PAD_BITS );	/* floating padded-array 2nd element index is here */

	#ifndef USE_SSE2	// Scalar-double mode:

	/*************************************************************/
	/*                  1st set of sincos:                       */
	/*************************************************************/
		iroot = index[k++];
		l = iroot;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cA01 =rt;	sA01 =it;

	re= rt;	im= it;	/* Save for the wrapper Step... */

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cA02 =rt;	sA02 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += 4*iroot;			/* 7*iroot	*/
		iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cA03 =rt;	sA03 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cA07 =rt;	sA07 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cA0E =rt;	sA0E =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cA15 =rt;	sA15 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cA1C =rt;	sA1C =it;

		/*(c,s)4-10 (decimal indices), 4-0A (hexadecimal indices):	*/
		t00=cA01*cA07;	t01=cA01*sA07;	t02=sA01*cA07;	t03=sA01*sA07;
		cA06=t00+t03;	sA06=t01-t02;	cA08=t00-t03;	sA08=t01+t02;

		t00=cA02*cA07;	t01=cA02*sA07;	t02=sA02*cA07;	t03=sA02*sA07;
		cA05=t00+t03;	sA05=t01-t02;	cA09=t00-t03;	sA09=t01+t02;

		t00=cA03*cA07;	t01=cA03*sA07;	t02=sA03*cA07;	t03=sA03*sA07;
		cA04=t00+t03;	sA04=t01-t02;	cA0A=t00-t03;	sA0A=t01+t02;

		/*(c,s)11-17 (decimal indices), 0B-11 (hexadecimal indices):;	*/
		t00=cA01*cA0E;	t01=cA01*sA0E;	t02=sA01*cA0E;	t03=sA01*sA0E;
		cA0D=t00+t03;	sA0D=t01-t02;	cA0F=t00-t03;	sA0F=t01+t02;

		t00=cA02*cA0E;	t01=cA02*sA0E;	t02=sA02*cA0E;	t03=sA02*sA0E;
		cA0C=t00+t03;	sA0C=t01-t02;	cA10=t00-t03;	sA10=t01+t02;

		t00=cA03*cA0E;	t01=cA03*sA0E;	t02=sA03*cA0E;	t03=sA03*sA0E;
		cA0B=t00+t03;	sA0B=t01-t02;	cA11=t00-t03;	sA11=t01+t02;

		/*(c,s)18-24 (decimal indices), 12-18 (hexadecimal indices):	*/
		t00=cA01*cA15;	t01=cA01*sA15;	t02=sA01*cA15;	t03=sA01*sA15;
		cA14=t00+t03;	sA14=t01-t02;	cA16=t00-t03;	sA16=t01+t02;

		t00=cA02*cA15;	t01=cA02*sA15;	t02=sA02*cA15;	t03=sA02*sA15;
		cA13=t00+t03;	sA13=t01-t02;	cA17=t00-t03;	sA17=t01+t02;

		t00=cA03*cA15;	t01=cA03*sA15;	t02=sA03*cA15;	t03=sA03*sA15;
		cA12=t00+t03;	sA12=t01-t02;	cA18=t00-t03;	sA18=t01+t02;

		/*(c,s)25-31 (decimal indices), 19-1F (hexadecimal indices):	*/
		t00=cA01*cA1C;	t01=cA01*sA1C;	t02=sA01*cA1C;	t03=sA01*sA1C;
		cA1B=t00+t03;	sA1B=t01-t02;	cA1D=t00-t03;	sA1D=t01+t02;

		t00=cA02*cA1C;	t01=cA02*sA1C;	t02=sA02*cA1C;	t03=sA02*sA1C;
		cA1A=t00+t03;	sA1A=t01-t02;	cA1E=t00-t03;	sA1E=t01+t02;

		t00=cA03*cA1C;	t01=cA03*sA1C;	t02=sA03*cA1C;	t03=sA03*sA1C;
		cA19=t00+t03;	sA19=t01-t02;	cA1F=t00-t03;	sA1F=t01+t02;

	/*************************************************************/
	/*                  2nd set of sincos:                       */
	/*************************************************************/
		iroot = index[k++];
		l = iroot;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;			/* 2*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cB01 =rt;	sB01 =it;

	if(j1 == 0){ re= rt;	im= it; }   /* The j1 = 0 case is special... */

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cB02 =rt;	sB02 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += 4*iroot;			/* 7*iroot	*/
		iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cB03 =rt;	sB03 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cB07 =rt;	sB07 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cB0E =rt;	sB0E =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cB15 =rt;	sB15 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		cB1C =rt;	sB1C =it;

		/*(c,s)4-10 (decimal indices), 4-0A (hexadecimal indices):	*/
		t00=cB01*cB07;	t01=cB01*sB07;	t02=sB01*cB07;	t03=sB01*sB07;
		cB06=t00+t03;	sB06=t01-t02;	cB08=t00-t03;	sB08=t01+t02;

		t00=cB02*cB07;	t01=cB02*sB07;	t02=sB02*cB07;	t03=sB02*sB07;
		cB05=t00+t03;	sB05=t01-t02;	cB09=t00-t03;	sB09=t01+t02;

		t00=cB03*cB07;	t01=cB03*sB07;	t02=sB03*cB07;	t03=sB03*sB07;
		cB04=t00+t03;	sB04=t01-t02;	cB0A=t00-t03;	sB0A=t01+t02;

		/*(c,s)11-17 (decimal indices), 0B-11 (hexadecimal indices):;	*/
		t00=cB01*cB0E;	t01=cB01*sB0E;	t02=sB01*cB0E;	t03=sB01*sB0E;
		cB0D=t00+t03;	sB0D=t01-t02;	cB0F=t00-t03;	sB0F=t01+t02;

		t00=cB02*cB0E;	t01=cB02*sB0E;	t02=sB02*cB0E;	t03=sB02*sB0E;
		cB0C=t00+t03;	sB0C=t01-t02;	cB10=t00-t03;	sB10=t01+t02;

		t00=cB03*cB0E;	t01=cB03*sB0E;	t02=sB03*cB0E;	t03=sB03*sB0E;
		cB0B=t00+t03;	sB0B=t01-t02;	cB11=t00-t03;	sB11=t01+t02;

		/*(c,s)18-24 (decimal indices), 12-18 (hexadecimal indices):	*/
		t00=cB01*cB15;	t01=cB01*sB15;	t02=sB01*cB15;	t03=sB01*sB15;
		cB14=t00+t03;	sB14=t01-t02;	cB16=t00-t03;	sB16=t01+t02;

		t00=cB02*cB15;	t01=cB02*sB15;	t02=sB02*cB15;	t03=sB02*sB15;
		cB13=t00+t03;	sB13=t01-t02;	cB17=t00-t03;	sB17=t01+t02;

		t00=cB03*cB15;	t01=cB03*sB15;	t02=sB03*cB15;	t03=sB03*sB15;
		cB12=t00+t03;	sB12=t01-t02;	cB18=t00-t03;	sB18=t01+t02;

		/*(c,s)25-31 (decimal indices), 19-1F (hexadecimal indices):	*/
		t00=cB01*cB1C;	t01=cB01*sB1C;	t02=sB01*cB1C;	t03=sB01*sB1C;
		cB1B=t00+t03;	sB1B=t01-t02;	cB1D=t00-t03;	sB1D=t01+t02;

		t00=cB02*cB1C;	t01=cB02*sB1C;	t02=sB02*cB1C;	t03=sB02*sB1C;
		cB1A=t00+t03;	sB1A=t01-t02;	cB1E=t00-t03;	sB1E=t01+t02;

		t00=cB03*cB1C;	t01=cB03*sB1C;	t02=sB03*cB1C;	t03=sB03*sB1C;
		cB19=t00+t03;	sB19=t01-t02;	cB1F=t00-t03;	sB1F=t01+t02;

	#else	// SIMD:

		/* Due to roots-locality considerations, roots (c,s)[0-31] are offset w.r.to the thread-local ptr pair as
					cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
		(cc0,ss0) + 0x[06,26,16,36,0e,2e,1e,3e,0a,2a,1a,3a,12,32,22,42,08,28,18,38,10,30,20,40,0c,2c,1c,3c,14,34,24,44].
		Here, due to the need to compute a new set of roots for each set of inputs, we use a streamlined sequence which
		computes only the [ 0, 1, 2, 3, 7,14,21,28]th roots with maximal accuracy (i.e. using 2-table-multiply),
		then generates the remaining ones from those. Thus the needed pointer offsets below are
			(cc0,ss0) + 0x[06,26,16,36,3e,22,30,14]:
		*/
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);

	  #if defined(USE_SSE2) && !defined(USE_AVX)
	// In SSE2 mode need 2 sets of sincos - since we may use different indexing patterns for storing k1,2_arr data
	// in AVX+ mode, don't share code for first 2 sets with those wider SIMD cases:

		// 1st set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 0] = k1<<4;	k2_arr[ 0] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 2] = k1<<4;	k2_arr[ 2] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 4] = k1<<4;	k2_arr[ 4] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 6] = k1<<4;	k2_arr[ 6] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 8] = k1<<4;	k2_arr[ 8] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[10] = k1<<4;	k2_arr[10] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[12] = k1<<4;	k2_arr[12] = k2<<4;

		// 2nd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 1] = k1<<4;	k2_arr[ 1] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 3] = k1<<4;	k2_arr[ 3] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 5] = k1<<4;	k2_arr[ 5] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 7] = k1<<4;	k2_arr[ 7] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 9] = k1<<4;	k2_arr[ 9] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[11] = k1<<4;	k2_arr[11] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[13] = k1<<4;	k2_arr[13] = k2<<4;

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX32_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #elif !defined(USE_AVX512)	// AVX/AVX2:

	  // In AVX mode need 4 sets of sincos:
		// Need to explicitly 0 these for the first-blocks case in AVX mode to make sure the
		// unused-in-those-cases 3rd/4th-set indices stay nice and non-segfault-y:
		if(j1 <= 320) {
			memset(k1_arr, 0, 14*RE_IM_STRIDE*sizeof(uint32));
			memset(k2_arr, 0, 14*RE_IM_STRIDE*sizeof(uint32));
		}
  #if 0
	printf("%s: Computing Twiddles with k*_arr-index K = %d, Iroot = %d\n",func,k,index[k]);
  #endif
		// 1st set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 0] = k1<<4;	k2_arr[ 0] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 4] = k1<<4;	k2_arr[ 4] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 8] = k1<<4;	k2_arr[ 8] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[12] = k1<<4;	k2_arr[12] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[16] = k1<<4;	k2_arr[16] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[20] = k1<<4;	k2_arr[20] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[24] = k1<<4;	k2_arr[24] = k2<<4;

		// 2nd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 1] = k1<<4;	k2_arr[ 1] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 5] = k1<<4;	k2_arr[ 5] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 9] = k1<<4;	k2_arr[ 9] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[13] = k1<<4;	k2_arr[13] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[17] = k1<<4;	k2_arr[17] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[21] = k1<<4;	k2_arr[21] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[25] = k1<<4;	k2_arr[25] = k2<<4;

	  if(j1 > 128)	// Sincos data for the 2 starting scalar-mode data blocks get done in SSE2 mode, i.e. only using d0,d1
	  {
		// 3rd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 2] = k1<<4;	k2_arr[ 2] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 6] = k1<<4;	k2_arr[ 6] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[10] = k1<<4;	k2_arr[10] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[14] = k1<<4;	k2_arr[14] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[18] = k1<<4;	k2_arr[18] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[22] = k1<<4;	k2_arr[22] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[26] = k1<<4;	k2_arr[26] = k2<<4;

		// 4th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 3] = k1<<4;	k2_arr[ 3] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 7] = k1<<4;	k2_arr[ 7] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[11] = k1<<4;	k2_arr[11] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[15] = k1<<4;	k2_arr[15] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[19] = k1<<4;	k2_arr[19] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[23] = k1<<4;	k2_arr[23] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[27] = k1<<4;	k2_arr[27] = k2<<4;
	  }	// endif(j1 > 128)

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX32_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #elif defined(USE_AVX512)	// AVX512:

	  // In AVX-512 mode need 8 sets of sincos:
		// Need to explicitly 0 these for the first-few-blocks-done-via-scalar-code case in AVX mode to make sure the
		// unused-in-those-cases 3rd/4th-set indices stay nice and non-segfault-y:
		if(j1 <= 320) {
			memset(k1_arr, 0, 14*RE_IM_STRIDE*sizeof(uint32));
			memset(k2_arr, 0, 14*RE_IM_STRIDE*sizeof(uint32));
		}

	   #ifdef USE_IMCI512	// Vectorized version below a no-go on 1st-gen Xeon Phi

		#warning AVX-512: Using slower non-ASM-macro version of sincos-indexing in radix32_dyadic_square.c.
		/*
		IMCI512: SSE2_RADIX16_CALC_TWIDDLES_LOACC reads k[0|1]_arr[0-55] in 8-index chunks,
			does gather-loads associated rt[0|1] elts, computes twiddles, writes to cc0+[0x8-...].
			Pull the gather-loads out here and do in C, just do the 512-bit vector CMULs in asm-macro:
		*/
		// 1st set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		tmp = cc0+0x8; tm1 = cc0+0x9;	c_tmp = cc0+0xa; s_tmp = cc0+0xb;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d0 = rt0[k1].re; tm1->d0 = rt0[k1].im;	c_tmp->d0 = rt1[k2].re; s_tmp->d0 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 2nd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		tmp = cc0+0x8; tm1 = cc0+0x9;	c_tmp = cc0+0xa; s_tmp = cc0+0xb;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d1 = rt0[k1].re; tm1->d1 = rt0[k1].im;	c_tmp->d1 = rt1[k2].re; s_tmp->d1 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

	  if(j1 > 128) {	// Sincos data for the initial done-in-scalar-mode data blocks get done in SSE2 mode, i.e. only using d0,d1
		// 3rd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		tmp = cc0+0x8; tm1 = cc0+0x9;	c_tmp = cc0+0xa; s_tmp = cc0+0xb;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d2 = rt0[k1].re; tm1->d2 = rt0[k1].im;	c_tmp->d2 = rt1[k2].re; s_tmp->d2 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 4th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		tmp = cc0+0x8; tm1 = cc0+0x9;	c_tmp = cc0+0xa; s_tmp = cc0+0xb;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d3 = rt0[k1].re; tm1->d3 = rt0[k1].im;	c_tmp->d3 = rt1[k2].re; s_tmp->d3 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
	  }	// endif(j1 > 128)
	  if(j1 > 320) {	// Must skip these sincos-inits until shift to SIMD-mode DFT, to keep array index k from being incremented
		// 5th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		tmp = cc0+0x8; tm1 = cc0+0x9;	c_tmp = cc0+0xa; s_tmp = cc0+0xb;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d4 = rt0[k1].re; tm1->d4 = rt0[k1].im;	c_tmp->d4 = rt1[k2].re; s_tmp->d4 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 6th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		tmp = cc0+0x8; tm1 = cc0+0x9;	c_tmp = cc0+0xa; s_tmp = cc0+0xb;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d5 = rt0[k1].re; tm1->d5 = rt0[k1].im;	c_tmp->d5 = rt1[k2].re; s_tmp->d5 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 7th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		tmp = cc0+0x8; tm1 = cc0+0x9;	c_tmp = cc0+0xa; s_tmp = cc0+0xb;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d6 = rt0[k1].re; tm1->d6 = rt0[k1].im;	c_tmp->d6 = rt1[k2].re; s_tmp->d6 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;

		// 8th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		tmp = cc0+0x8; tm1 = cc0+0x9;	c_tmp = cc0+0xa; s_tmp = cc0+0xb;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		tmp->d7 = rt0[k1].re; tm1->d7 = rt0[k1].im;	c_tmp->d7 = rt1[k2].re; s_tmp->d7 = rt1[k2].im;	tmp += 4; tm1 += 4; c_tmp += 4; s_tmp += 4;
	  }	// endif(j1 > 320)

	   #else

		#warning Using AVX-512 code in radix32_dyadic_square.c
		// 1st set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 0] = k1<<4;	k2_arr[ 0] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 8] = k1<<4;	k2_arr[ 8] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[16] = k1<<4;	k2_arr[16] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[24] = k1<<4;	k2_arr[24] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[32] = k1<<4;	k2_arr[32] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[40] = k1<<4;	k2_arr[40] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[48] = k1<<4;	k2_arr[48] = k2<<4;

		// 2nd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 1] = k1<<4;	k2_arr[ 1] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 9] = k1<<4;	k2_arr[ 9] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[17] = k1<<4;	k2_arr[17] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[25] = k1<<4;	k2_arr[25] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[33] = k1<<4;	k2_arr[33] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[41] = k1<<4;	k2_arr[41] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[49] = k1<<4;	k2_arr[49] = k2<<4;

	  if(j1 > 128) {	// Sincos data for the initial done-in-scalar-mode data blocks get done in SSE2 mode, i.e. only using d0,d1
		// 3rd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 2] = k1<<4;	k2_arr[ 2] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[10] = k1<<4;	k2_arr[10] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[18] = k1<<4;	k2_arr[18] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[26] = k1<<4;	k2_arr[26] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[34] = k1<<4;	k2_arr[34] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[42] = k1<<4;	k2_arr[42] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[50] = k1<<4;	k2_arr[50] = k2<<4;

		// 4th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 3] = k1<<4;	k2_arr[ 3] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[11] = k1<<4;	k2_arr[11] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[19] = k1<<4;	k2_arr[19] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[27] = k1<<4;	k2_arr[27] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[35] = k1<<4;	k2_arr[35] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[43] = k1<<4;	k2_arr[43] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[51] = k1<<4;	k2_arr[51] = k2<<4;
	  }	// endif(j1 > 128)
	  if(j1 > 320) {	// Must skip these sincos-inits until shift to SIMD-mode DFT, to keep array index k from being incremented
		// 5th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 4] = k1<<4;	k2_arr[ 4] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[12] = k1<<4;	k2_arr[12] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[20] = k1<<4;	k2_arr[20] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[28] = k1<<4;	k2_arr[28] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[36] = k1<<4;	k2_arr[36] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[44] = k1<<4;	k2_arr[44] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[52] = k1<<4;	k2_arr[52] = k2<<4;

		// 6th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 5] = k1<<4;	k2_arr[ 5] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[13] = k1<<4;	k2_arr[13] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[21] = k1<<4;	k2_arr[21] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[29] = k1<<4;	k2_arr[29] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[37] = k1<<4;	k2_arr[37] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[45] = k1<<4;	k2_arr[45] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[53] = k1<<4;	k2_arr[53] = k2<<4;

		// 7th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 6] = k1<<4;	k2_arr[ 6] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[14] = k1<<4;	k2_arr[14] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[22] = k1<<4;	k2_arr[22] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[30] = k1<<4;	k2_arr[30] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[38] = k1<<4;	k2_arr[38] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[46] = k1<<4;	k2_arr[46] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[54] = k1<<4;	k2_arr[54] = k2<<4;

		// 8th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 7] = k1<<4;	k2_arr[ 7] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[15] = k1<<4;	k2_arr[15] = k2<<4;
		l += iroot;			/* 3*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[23] = k1<<4;	k2_arr[23] = k2<<4;
		l += (iroot << 2);	/* 7*iroot */	iroot = l;
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[31] = k1<<4;	k2_arr[31] = k2<<4;
		l += iroot;			/* 14*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[39] = k1<<4;	k2_arr[39] = k2<<4;
		l += iroot;			/* 21*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[47] = k1<<4;	k2_arr[47] = k2<<4;
		l += iroot;			/* 28*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[55] = k1<<4;	k2_arr[55] = k2<<4;
	  }	// endif(j1 > 320)

	   #endif	// (IMCI512 or AVX512?) toggle

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX32_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #endif	// SIMD mode?

	#endif	// SIMD ?

	#ifdef USE_SSE2	// Both SSE2 and AVX share this:

	  #ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX512 mode

		// process 8 main-array blocks of [8 vec_dbl = 8 x 8 = 64 doubles] each in AVX512 mode, total = 64 vec_dbl = 32 vec_cmplx
		add0 = a + j1pad;
		add2 = add0 + 64;	// add2 = add0 + [64 doubles, equiv to 8 AVX-512 registers]
		add4 = add2 + 64;
		add6 = add4 + 64;
		add1 = a + j2pad;	// Last 4 offsets run in descending order for Mers-mod
		add3 = add1 - 64;
		add5 = add3 - 64;
		add7 = add5 - 64;

	  #elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  						// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:

		// process 4 main-array blocks of [16 vec_dbl = 16 x 4 = 64 doubles] each in AVX mode, total = 64 vec_dbl = 32 vec_cmplx
		add0 = a + j1pad;
		add1 = a + j2pad;
		add2 = add0 + 64;	// add2 = add0 + [64 doubles, equiv to 16 AVX registers]
		add3 = add1 - 64;	// Last 2 offsets run in descending order for Mers-mod

	  #else	// SSE2:

		// process 2 main-array blocks of [32 vec_dbl = 32 x 2 = 64 doubles] each in SSE2 mode, total = 64 vec_dbl = 32 vec_cmplx
		add0 = a + j1pad;
		add1 = a + j2pad;

	  #endif

	jump_new:
	  #ifdef USE_AVX
		// SSE2/AVX/AVX-512 need SIMD data processing to begin at j1 = 64,128,256, respectively. Here, 160
		// is the scalar-DFT-mode j1-value which is followed by j1 = 256 after block-index updating:
		if(j1 <= 320) {
		  if((j1 & 127) == 0) {	// Differentiate between j1 == 0 and 64 (mod 128) cases:
	  #else	// SSE2 begins SIMD-mode data processing at j1 = 256, which follows j1 = 128 after block-index-updating.
	  		// (Could also use same value for AVX, but using same breakover for AVX,AVX-512 makes for easier debug of the latter.
		if(j1 <= 128) {
	  #endif
			cA01 = c01->d0;	sA01 = (c01+1)->d0;		cB01 = c01->d1;	sB01 = (c01+1)->d1;
			cA02 = c02->d0;	sA02 = (c02+1)->d0;		cB02 = c02->d1;	sB02 = (c02+1)->d1;
			cA03 = c03->d0;	sA03 = (c03+1)->d0;		cB03 = c03->d1;	sB03 = (c03+1)->d1;
			cA04 = c04->d0;	sA04 = (c04+1)->d0;		cB04 = c04->d1;	sB04 = (c04+1)->d1;
			cA05 = c05->d0;	sA05 = (c05+1)->d0;		cB05 = c05->d1;	sB05 = (c05+1)->d1;
			cA06 = c06->d0;	sA06 = (c06+1)->d0;		cB06 = c06->d1;	sB06 = (c06+1)->d1;
			cA07 = c07->d0;	sA07 = (c07+1)->d0;		cB07 = c07->d1;	sB07 = (c07+1)->d1;
			cA08 = c08->d0;	sA08 = (c08+1)->d0;		cB08 = c08->d1;	sB08 = (c08+1)->d1;
			cA09 = c09->d0;	sA09 = (c09+1)->d0;		cB09 = c09->d1;	sB09 = (c09+1)->d1;
			cA0A = c0A->d0;	sA0A = (c0A+1)->d0;		cB0A = c0A->d1;	sB0A = (c0A+1)->d1;
			cA0B = c0B->d0;	sA0B = (c0B+1)->d0;		cB0B = c0B->d1;	sB0B = (c0B+1)->d1;
			cA0C = c0C->d0;	sA0C = (c0C+1)->d0;		cB0C = c0C->d1;	sB0C = (c0C+1)->d1;
			cA0D = c0D->d0;	sA0D = (c0D+1)->d0;		cB0D = c0D->d1;	sB0D = (c0D+1)->d1;
			cA0E = c0E->d0;	sA0E = (c0E+1)->d0;		cB0E = c0E->d1;	sB0E = (c0E+1)->d1;
			cA0F = c0F->d0;	sA0F = (c0F+1)->d0;		cB0F = c0F->d1;	sB0F = (c0F+1)->d1;
			cA10 = c10->d0;	sA10 = (c10+1)->d0;		cB10 = c10->d1;	sB10 = (c10+1)->d1;
			cA11 = c11->d0;	sA11 = (c11+1)->d0;		cB11 = c11->d1;	sB11 = (c11+1)->d1;
			cA12 = c12->d0;	sA12 = (c12+1)->d0;		cB12 = c12->d1;	sB12 = (c12+1)->d1;
			cA13 = c13->d0;	sA13 = (c13+1)->d0;		cB13 = c13->d1;	sB13 = (c13+1)->d1;
			cA14 = c14->d0;	sA14 = (c14+1)->d0;		cB14 = c14->d1;	sB14 = (c14+1)->d1;
			cA15 = c15->d0;	sA15 = (c15+1)->d0;		cB15 = c15->d1;	sB15 = (c15+1)->d1;
			cA16 = c16->d0;	sA16 = (c16+1)->d0;		cB16 = c16->d1;	sB16 = (c16+1)->d1;
			cA17 = c17->d0;	sA17 = (c17+1)->d0;		cB17 = c17->d1;	sB17 = (c17+1)->d1;
			cA18 = c18->d0;	sA18 = (c18+1)->d0;		cB18 = c18->d1;	sB18 = (c18+1)->d1;
			cA19 = c19->d0;	sA19 = (c19+1)->d0;		cB19 = c19->d1;	sB19 = (c19+1)->d1;
			cA1A = c1A->d0;	sA1A = (c1A+1)->d0;		cB1A = c1A->d1;	sB1A = (c1A+1)->d1;
			cA1B = c1B->d0;	sA1B = (c1B+1)->d0;		cB1B = c1B->d1;	sB1B = (c1B+1)->d1;
			cA1C = c1C->d0;	sA1C = (c1C+1)->d0;		cB1C = c1C->d1;	sB1C = (c1C+1)->d1;
			cA1D = c1D->d0;	sA1D = (c1D+1)->d0;		cB1D = c1D->d1;	sB1D = (c1D+1)->d1;
			cA1E = c1E->d0;	sA1E = (c1E+1)->d0;		cB1E = c1E->d1;	sB1E = (c1E+1)->d1;
			cA1F = c1F->d0;	sA1F = (c1F+1)->d0;		cB1F = c1F->d1;	sB1F = (c1F+1)->d1;
	  #ifdef USE_AVX
		  } else {	// j1 == 64 (mod 128):
			cA01 = c01->d2;	sA01 = (c01+1)->d2;		cB01 = c01->d3;	sB01 = (c01+1)->d3;
			cA02 = c02->d2;	sA02 = (c02+1)->d2;		cB02 = c02->d3;	sB02 = (c02+1)->d3;
			cA03 = c03->d2;	sA03 = (c03+1)->d2;		cB03 = c03->d3;	sB03 = (c03+1)->d3;
			cA04 = c04->d2;	sA04 = (c04+1)->d2;		cB04 = c04->d3;	sB04 = (c04+1)->d3;
			cA05 = c05->d2;	sA05 = (c05+1)->d2;		cB05 = c05->d3;	sB05 = (c05+1)->d3;
			cA06 = c06->d2;	sA06 = (c06+1)->d2;		cB06 = c06->d3;	sB06 = (c06+1)->d3;
			cA07 = c07->d2;	sA07 = (c07+1)->d2;		cB07 = c07->d3;	sB07 = (c07+1)->d3;
			cA08 = c08->d2;	sA08 = (c08+1)->d2;		cB08 = c08->d3;	sB08 = (c08+1)->d3;
			cA09 = c09->d2;	sA09 = (c09+1)->d2;		cB09 = c09->d3;	sB09 = (c09+1)->d3;
			cA0A = c0A->d2;	sA0A = (c0A+1)->d2;		cB0A = c0A->d3;	sB0A = (c0A+1)->d3;
			cA0B = c0B->d2;	sA0B = (c0B+1)->d2;		cB0B = c0B->d3;	sB0B = (c0B+1)->d3;
			cA0C = c0C->d2;	sA0C = (c0C+1)->d2;		cB0C = c0C->d3;	sB0C = (c0C+1)->d3;
			cA0D = c0D->d2;	sA0D = (c0D+1)->d2;		cB0D = c0D->d3;	sB0D = (c0D+1)->d3;
			cA0E = c0E->d2;	sA0E = (c0E+1)->d2;		cB0E = c0E->d3;	sB0E = (c0E+1)->d3;
			cA0F = c0F->d2;	sA0F = (c0F+1)->d2;		cB0F = c0F->d3;	sB0F = (c0F+1)->d3;
			cA10 = c10->d2;	sA10 = (c10+1)->d2;		cB10 = c10->d3;	sB10 = (c10+1)->d3;
			cA11 = c11->d2;	sA11 = (c11+1)->d2;		cB11 = c11->d3;	sB11 = (c11+1)->d3;
			cA12 = c12->d2;	sA12 = (c12+1)->d2;		cB12 = c12->d3;	sB12 = (c12+1)->d3;
			cA13 = c13->d2;	sA13 = (c13+1)->d2;		cB13 = c13->d3;	sB13 = (c13+1)->d3;
			cA14 = c14->d2;	sA14 = (c14+1)->d2;		cB14 = c14->d3;	sB14 = (c14+1)->d3;
			cA15 = c15->d2;	sA15 = (c15+1)->d2;		cB15 = c15->d3;	sB15 = (c15+1)->d3;
			cA16 = c16->d2;	sA16 = (c16+1)->d2;		cB16 = c16->d3;	sB16 = (c16+1)->d3;
			cA17 = c17->d2;	sA17 = (c17+1)->d2;		cB17 = c17->d3;	sB17 = (c17+1)->d3;
			cA18 = c18->d2;	sA18 = (c18+1)->d2;		cB18 = c18->d3;	sB18 = (c18+1)->d3;
			cA19 = c19->d2;	sA19 = (c19+1)->d2;		cB19 = c19->d3;	sB19 = (c19+1)->d3;
			cA1A = c1A->d2;	sA1A = (c1A+1)->d2;		cB1A = c1A->d3;	sB1A = (c1A+1)->d3;
			cA1B = c1B->d2;	sA1B = (c1B+1)->d2;		cB1B = c1B->d3;	sB1B = (c1B+1)->d3;
			cA1C = c1C->d2;	sA1C = (c1C+1)->d2;		cB1C = c1C->d3;	sB1C = (c1C+1)->d3;
			cA1D = c1D->d2;	sA1D = (c1D+1)->d2;		cB1D = c1D->d3;	sB1D = (c1D+1)->d3;
			cA1E = c1E->d2;	sA1E = (c1E+1)->d2;		cB1E = c1E->d3;	sB1E = (c1E+1)->d3;
			cA1F = c1F->d2;	sA1F = (c1F+1)->d2;		cB1F = c1F->d3;	sB1F = (c1F+1)->d3;
		  }
	  #endif
		if(j1 == 0){ re = cB01; im = sB01; } else { re = cA01; im = sA01; }   /* Save for the wrapper Step ... The j1 = 0 case is special. */

		//************ NOTE: The above if(j1 <= 128 [or 320]) clause is still open! In SIMD mode we use it to wrap
		// all the scalar-mode FFT/dyad-mul/iFFT code, which allows scalar-double and SIMD builds to share code.
		//**********************************************************************************************************

	#endif

	if(fwd_fft_only == 3)
		goto skip_fwd_fft;	// v20: jump-to-point for both-inputs-already-fwd-FFTed case

	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/
	/*
	Data layout comparison:
	A-index:0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31
Scalar-dbl:	r0	i0	r1	i1	r2	i2	r3	i3	r4	i4	r5	i5	r6	i6	r7	i7	r8	i8	r9	i9	r10	i10	r11	i11	r12	i12	r13	i13	r14	i14	r15	i15
	SSE:	r0	r1	i0	i1	r2	r3	i2	i3	r4	r5	i4	i5	r6	r7	i6	i7	r8	r9	i8	i9	r10	r11	i10	i11	r12	r13	i12	i13	r14	r15	i14	i15
	AVX:	r0	r1	r2	r3	i0	i1	i2	i3	r4	r5	r6	r7	i4	i5	i6	i7	r8	r9	r10	r11	i8	i9	i10	i11	r12	r13	r14	r15	i12	i13	i14	i15
	AVX-512:r0	r1	r2	r3	r4	r5	r6	r7	i0	i1	i2	i3	i4	i5	i6	i7	r8	r9	r10	r11	r12	r13	r14	r15	i8	i9	i10	i11	i12	i13	i14	i15
							[And same pattern for the next 16 complex inputs, just with the r,i-indices += 16,
							i.e. in computing idx-adjustments for rJ -> rI jumps just use I,J (mod 16) values.]

	Within the opening pass of four radix-8 sub-DFTs we combine complex input octets
		[0,16,8,24,4,20,12,28], [2,18,10,26,6,22,14,30], [1,17,9,25,5,21,13,29] and [3,19,11,27,7,23,15,31].
	Must fiddle the linear-idx jump w.r.to the index strides of the same inputs in the non-SIMD layout ... the #ifs below accomplish that:
	*/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;

	/*...Block 1:	*/
		t00=a[rdum   ];							t01=a[idum   ];							// z0
		rt =a[rdum+32]*cA10 - a[idum+32]*sA10;	it =a[idum+32]*cA10 + a[rdum+32]*sA10;	// z16
		t02=t00-rt;		t03=t01-it;
		t00=t00+rt;		t01=t01+it;
		t04=a[rdum+16]*cA08 - a[idum+16]*sA08;	t05=a[idum+16]*cA08 + a[rdum+16]*sA08;	// z8
		rt =a[rdum+48]*cA18 - a[idum+48]*sA18;	it =a[idum+48]*cA18 + a[rdum+48]*sA18;	// z24
		t06=t04-rt;		t07=t05-it;
		t04=t04+rt;		t05=t05+it;
		rt =t04;		it =t05;
		t04=t00-rt;		t05=t01-it;
		t00=t00+rt;		t01=t01+it;
		rt =t06;		it =t07;
		t06=t02+it;		t07=t03-rt;
		t02=t02-it;		t03=t03+rt;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;	// r24 -> r4 is -= 40 in scalar-double/SSE2/AVX modes, -= 44 in AVX-512 mode
	#endif
		t08=a[rdum+8 ]*cA04 - a[idum+8 ]*sA04;	t09=a[idum+8 ]*cA04 + a[rdum+8 ]*sA04;	// z4
		rt =a[rdum+40]*cA14 - a[idum+40]*sA14;	it =a[idum+40]*cA14 + a[rdum+40]*sA14;	// z20
		t0A=t08-rt;		t0B=t09-it;
		t08=t08+rt;		t09=t09+it;
		t0C=a[rdum+24]*cA0C - a[idum+24]*sA0C;	t0D=a[idum+24]*cA0C + a[rdum+24]*sA0C;	// z12
		rt =a[rdum+56]*cA1C - a[idum+56]*sA1C;	it =a[idum+56]*cA1C + a[rdum+56]*sA1C;	// z28
		t0E=t0C-rt;		t0F=t0D-it;
		t0C=t0C+rt;		t0D=t0D+it;
		rt =t0C;		it =t0D;
		t0C=t08-rt;		t0D=t09-it;
		t08=t08+rt;		t09=t09+it;
		rt =t0E;		it =t0F;
		t0E=t0A+it;		t0F=t0B-rt;
		t0A=t0A-it;		t0B=t0B+rt;
		rt =t08;		it =t09;
		t08=t00-rt;		t09=t01-it;
		t00=t00+rt;		t01=t01+it;
		rt =t0C;		it =t0D;
		t0C=t04+it;		t0D=t05-rt;
		t04=t04-it;		t05=t05+rt;
		rt =(t0A-t0B)*ISRT2;it =(t0A+t0B)*ISRT2;
		t0A=t02-rt;		t0B=t03-it;
		t02=t02+rt;		t03=t03+it;
		rt =(t0E+t0F)*ISRT2;it =(t0F-t0E)*ISRT2;
		t0E=t06+rt;		t0F=t07+it;
		t06=t06-rt;		t07=t07-it;
	/*...Block 2:	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;	// r28 -> r2 jump in AVX-512 is 4 fewer leftward than in AVX, hence idx-adjustment here is AVX value += 4
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;	// r28 -> r2 jump same in scalar-double and SSE2 modes, AVX is 2 more places leftward
	#endif
		t10=a[rdum+4 ]*cA02 - a[idum+4 ]*sA02;	t11=a[idum+4 ]*cA02 + a[rdum+4 ]*sA02;	// z2
		rt =a[rdum+36]*cA12 - a[idum+36]*sA12;	it =a[idum+36]*cA12 + a[rdum+36]*sA12;	// z18
		t12=t10-rt;		t13=t11-it;
		t10=t10+rt;		t11=t11+it;
		t14=a[rdum+20]*cA0A - a[idum+20]*sA0A;	t15=a[idum+20]*cA0A + a[rdum+20]*sA0A;	// z10
		rt =a[rdum+52]*cA1A - a[idum+52]*sA1A;	it =a[idum+52]*cA1A + a[rdum+52]*sA1A;	// z26
		t16=t14-rt;		t17=t15-it;
		t14=t14+rt;		t15=t15+it;
		rt =t14;		it =t15;
		t14=t10-rt;		t15=t11-it;
		t10=t10+rt;		t11=t11+it;
		rt =t16;		it =t17;
		t16=t12+it;		t17=t13-rt;
		t12=t12-it;		t13=t13+rt;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t18=a[rdum+12]*cA06 - a[idum+12]*sA06;	t19=a[idum+12]*cA06 + a[rdum+12]*sA06;	// z6
		rt =a[rdum+44]*cA16 - a[idum+44]*sA16;	it =a[idum+44]*cA16 + a[rdum+44]*sA16;	// z22
		t1A=t18-rt;		t1B=t19-it;
		t18=t18+rt;		t19=t19+it;
		t1C=a[rdum+28]*cA0E - a[idum+28]*sA0E;	t1D=a[idum+28]*cA0E + a[rdum+28]*sA0E;	// z14
		rt =a[rdum+60]*cA1E - a[idum+60]*sA1E;	it =a[idum+60]*cA1E + a[rdum+60]*sA1E;	// z30
		t1E=t1C-rt;		t1F=t1D-it;
		t1C=t1C+rt;		t1D=t1D+it;
		rt =t1C;		it =t1D;
		t1C=t18-rt;		t1D=t19-it;
		t18=t18+rt;		t19=t19+it;
		rt =t1E;		it =t1F;
		t1E=t1A+it;		t1F=t1B-rt;
		t1A=t1A-it;		t1B=t1B+rt;
		rt =t18;		it =t19;
		t18=t10-rt;		t19=t11-it;
		t10=t10+rt;		t11=t11+it;
		rt =t1C;		it =t1D;
		t1C=t14+it;		t1D=t15-rt;
		t14=t14-it;		t15=t15+rt;
		rt =(t1A-t1B)*ISRT2;it =(t1A+t1B)*ISRT2;
		t1A=t12-rt;		t1B=t13-it;
		t12=t12+rt;		t13=t13+it;
		rt =(t1E+t1F)*ISRT2;it =(t1F-t1E)*ISRT2;
		t1E=t16+rt;		t1F=t17+it;
		t16=t16-rt;		t17=t17-it;
	/*...Block 3:	*/
	#ifdef USE_AVX512
		rdum += 5;	idum += 5;	// r30 -> r1 jump in AVX-512 is 4 fewer leftward than in AVX, hence idx-adjustment here is AVX value += 4
	#elif defined(USE_AVX)
		++rdum;	++idum;			// r30 -> r1 jump same in AVX mode is 2 fewer places leftward than SSE2
	#elif defined(USE_SSE2)
		--rdum;	--idum;			// r30 -> r1 jump in SSE2 mode is 1 more places leftward than in scalar-double mode
	#endif
		t20=a[rdum+2 ]*cA01 - a[idum+2 ]*sA01;	t21=a[idum+2 ]*cA01 + a[rdum+2 ]*sA01;	// z1
		rt =a[rdum+34]*cA11 - a[idum+34]*sA11;	it =a[idum+34]*cA11 + a[rdum+34]*sA11;	// z17
		t22=t20-rt;		t23=t21-it;
		t20=t20+rt;		t21=t21+it;
		t24=a[rdum+18]*cA09 - a[idum+18]*sA09;	t25=a[idum+18]*cA09 + a[rdum+18]*sA09;	// z9
		rt =a[rdum+50]*cA19 - a[idum+50]*sA19;	it =a[idum+50]*cA19 + a[rdum+50]*sA19;	// z25
		t26=t24-rt;		t27=t25-it;
		t24=t24+rt;		t25=t25+it;
		rt =t24;		it =t25;
		t24=t20-rt;		t25=t21-it;
		t20=t20+rt;		t21=t21+it;
		rt =t26;		it =t27;
		t26=t22+it;		t27=t23-rt;
		t22=t22-it;		t23=t23+rt;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t28=a[rdum+10]*cA05 - a[idum+10]*sA05;	t29=a[idum+10]*cA05 + a[rdum+10]*sA05;	// z5
		rt =a[rdum+42]*cA15 - a[idum+42]*sA15;	it =a[idum+42]*cA15 + a[rdum+42]*sA15;	// z21
		t2A=t28-rt;		t2B=t29-it;
		t28=t28+rt;		t29=t29+it;
		t2C=a[rdum+26]*cA0D - a[idum+26]*sA0D;	t2D=a[idum+26]*cA0D + a[rdum+26]*sA0D;	// z13
		rt =a[rdum+58]*cA1D - a[idum+58]*sA1D;	it =a[idum+58]*cA1D + a[rdum+58]*sA1D;	// z29
		t2E=t2C-rt;		t2F=t2D-it;
		t2C=t2C+rt;		t2D=t2D+it;
		rt =t2C;		it =t2D;
		t2C=t28-rt;		t2D=t29-it;
		t28=t28+rt;		t29=t29+it;
		rt =t2E;		it =t2F;
		t2E=t2A+it;		t2F=t2B-rt;
		t2A=t2A-it;		t2B=t2B+rt;
		rt =t28;		it =t29;
		t28=t20-rt;		t29=t21-it;
		t20=t20+rt;		t21=t21+it;
		rt =t2C;		it =t2D;
		t2C=t24+it;		t2D=t25-rt;
		t24=t24-it;		t25=t25+rt;
		rt =(t2A-t2B)*ISRT2;it =(t2A+t2B)*ISRT2;
		t2A=t22-rt;		t2B=t23-it;
		t22=t22+rt;		t23=t23+it;
		rt =(t2E+t2F)*ISRT2;it =(t2F-t2E)*ISRT2;
		t2E=t26+rt;		t2F=t27+it;
		t26=t26-rt;		t27=t27-it;
	/*...Block 4:	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;	// r29 -> r3 jump in AVX-512 is 4 fewer leftward than in AVX, hence idx-adjustment here is AVX value += 4
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;	// r29 -> r3 jump same in scalar-double and SSE2 modes, AVX is 2 more places leftward
	#endif
		t30=a[rdum+6 ]*cA03 - a[idum+6 ]*sA03;	t31=a[idum+6 ]*cA03 + a[rdum+6 ]*sA03;	// z3
		rt =a[rdum+38]*cA13 - a[idum+38]*sA13;	it =a[idum+38]*cA13 + a[rdum+38]*sA13;	// z19
		t32=t30-rt;		t33=t31-it;
		t30=t30+rt;		t31=t31+it;
		t34=a[rdum+22]*cA0B - a[idum+22]*sA0B;	t35=a[idum+22]*cA0B + a[rdum+22]*sA0B;	// z11
		rt =a[rdum+54]*cA1B - a[idum+54]*sA1B;	it =a[idum+54]*cA1B + a[rdum+54]*sA1B;	// z27
		t36=t34-rt;		t37=t35-it;
		t34=t34+rt;		t35=t35+it;
		rt =t34;		it =t35;
		t34=t30-rt;		t35=t31-it;
		t30=t30+rt;		t31=t31+it;
		rt =t36;		it =t37;
		t36=t32+it;		t37=t33-rt;
		t32=t32-it;		t33=t33+rt;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t38=a[rdum+14]*cA07 - a[idum+14]*sA07;	t39=a[idum+14]*cA07 + a[rdum+14]*sA07;	// z7
		rt =a[rdum+46]*cA17 - a[idum+46]*sA17;	it =a[idum+46]*cA17 + a[rdum+46]*sA17;	// z23
		t3A=t38-rt;		t3B=t39-it;
		t38=t38+rt;		t39=t39+it;
		t3C=a[rdum+30]*cA0F - a[idum+30]*sA0F;	t3D=a[idum+30]*cA0F + a[rdum+30]*sA0F;	// z15
		rt =a[rdum+62]*cA1F - a[idum+62]*sA1F;	it =a[idum+62]*cA1F + a[rdum+62]*sA1F;	// z31
		t3E=t3C-rt;		t3F=t3D-it;
		t3C=t3C+rt;		t3D=t3D+it;
		rt =t3C;		it =t3D;
		t3C=t38-rt;		t3D=t39-it;
		t38=t38+rt;		t39=t39+it;
		rt =t3E;		it =t3F;
		t3E=t3A+it;		t3F=t3B-rt;
		t3A=t3A-it;		t3B=t3B+rt;
		rt =t38;		it =t39;
		t38=t30-rt;		t39=t31-it;
		t30=t30+rt;		t31=t31+it;
		rt =t3C;		it =t3D;
		t3C=t34+it;		t3D=t35-rt;
		t34=t34-it;		t35=t35+rt;
		rt =(t3A-t3B)*ISRT2;it =(t3A+t3B)*ISRT2;
		t3A=t32-rt;		t3B=t33-it;
		t32=t32+rt;		t33=t33+it;
		rt =(t3E+t3F)*ISRT2;it =(t3F-t3E)*ISRT2;
		t3E=t36+rt;		t3F=t37+it;
		t36=t36-rt;		t37=t37-it;
	//...and now do eight radix-4 transforms, including the internal twiddle factors:
	/*...Block 1: t00,t10,t20,t30	*/
		rt =t10;	t10=t00-rt;	t00=t00+rt;
		it =t11;	t11=t01-it;	t01=t01+it;
		rt =t30;	t30=t20-rt;	t20=t20+rt;
		it =t31;	t31=t21-it;	t21=t21+it;
		a1p00r = t00+t20;		a1p00i = t01+t21;
		a1p01r = t00-t20;		a1p01i = t01-t21;
		a1p02r = t10-t31;		a1p02i = t11+t30;
		a1p03r = t10+t31;		a1p03i = t11-t30;
	/*...Block 5: t08,t18,t28,t38	*/
		rt =t18;	t18=t08+t19;	t08=t08-t19;
			t19=t09-rt;	t09=t09+rt;
		rt =(t28-t29)*ISRT2;	t29=(t28+t29)*ISRT2;		t28=rt;
		rt =(t39+t38)*ISRT2;	it =(t39-t38)*ISRT2;
		t38=t28+rt;			t28=t28-rt;
		t39=t29+it;			t29=t29-it;
		a1p04r = t08+t28;		a1p04i = t09+t29;
		a1p05r = t08-t28;		a1p05i = t09-t29;
		a1p06r = t18-t39;		a1p06i = t19+t38;
		a1p07r = t18+t39;		a1p07i = t19-t38;
	/*...Block 3: t04,t14,t24,t34	*/
		rt =(t14-t15)*ISRT2;	it =(t14+t15)*ISRT2;
		t14=t04-rt;			t04=t04+rt;
		t15=t05-it;			t05=t05+it;
		rt =t24*c - t25*s;		t25=t25*c + t24*s;		t24=rt;
		rt =t34*s - t35*c;		it =t35*s + t34*c;
		t34=t24-rt;			t24=t24+rt;
		t35=t25-it;			t25=t25+it;
		a1p08r = t04+t24;		a1p08i = t05+t25;
		a1p09r = t04-t24;		a1p09i = t05-t25;
		a1p0Ar = t14-t35;		a1p0Ai = t15+t34;
		a1p0Br = t14+t35;		a1p0Bi = t15-t34;
	/*...Block 7: t0C,t1C,t2C,t3C	*/
		rt =(t1D+t1C)*ISRT2;	it =(t1D-t1C)*ISRT2;
		t1C=t0C+rt;			t0C=t0C-rt;
		t1D=t0D+it;			t0D=t0D-it;
		rt =t2C*s - t2D*c;		t2D=t2D*s + t2C*c;		t2C=rt;
		rt =t3C*c - t3D*s;		it =t3D*c + t3C*s;
		t3C=t2C+rt;			t2C=t2C-rt;
		t3D=t2D+it;			t2D=t2D-it;
		a1p0Cr = t0C+t2C;		a1p0Ci = t0D+t2D;
		a1p0Dr = t0C-t2C;		a1p0Di = t0D-t2D;
		a1p0Er = t1C-t3D;		a1p0Ei = t1D+t3C;
		a1p0Fr = t1C+t3D;		a1p0Fi = t1D-t3C;
	/*...Block 2: t02,t12,t22,t32	*/
		rt =t12*c - t13*s;		it =t13*c + t12*s;
		t12=t02-rt;			t02=t02+rt;
		t13=t03-it;			t03=t03+it;
		rt =t22*c32_1 - t23*s32_1;	t23=t23*c32_1 + t22*s32_1;	t22=rt;
		rt =t32*c32_3 - t33*s32_3;	it =t33*c32_3 + t32*s32_3;
		t32=t22-rt;			t22=t22+rt;
		t33=t23-it;			t23=t23+it;
		a1p10r = t02+t22;		a1p10i = t03+t23;
		a1p11r = t02-t22;		a1p11i = t03-t23;
		a1p12r = t12-t33;		a1p12i = t13+t32;
		a1p13r = t12+t33;		a1p13i = t13-t32;
	/*...Block 6: t0A,t1A,t2A,t3A	*/
		rt =t1A*s + t1B*c;		it =t1B*s - t1A*c;
		t1A=t0A+rt;			t0A =t0A-rt;
		t1B=t0B+it;			t0B =t0B-it;
		rt =t2A*s32_3 - t2B*c32_3;	t2B=t2B*s32_3 + t2A*c32_3;	t2A=rt;
		rt =t3A*c32_1 + t3B*s32_1;	it =t3B*c32_1 - t3A*s32_1;
		t3A=t2A+rt;			t2A=t2A-rt;
		t3B=t2B+it;			t2B=t2B-it;
		a1p14r = t0A+t2A;		a1p14i = t0B+t2B;
		a1p15r = t0A-t2A;		a1p15i = t0B-t2B;
		a1p16r = t1A-t3B;		a1p16i = t1B+t3A;
		a1p17r = t1A+t3B;		a1p17i = t1B-t3A;
	/*...Block 4: t06,t16,t26,t36	*/
		rt =t16*s - t17*c;		it =t17*s + t16*c;
		t16=t06-rt;			t06 =t06+rt;
		t17=t07-it;			t07 =t07+it;
		rt =t26*c32_3 - t27*s32_3;	t27=t27*c32_3 + t26*s32_3;	t26=rt;
		rt =t36*s32_1 + t37*c32_1;	it =t37*s32_1 - t36*c32_1;
		t36=t26+rt;			t26=t26-rt;
		t37=t27+it;			t27=t27-it;
		a1p18r = t06+t26;		a1p18i = t07+t27;
		a1p19r = t06-t26;		a1p19i = t07-t27;
		a1p1Ar = t16-t37;		a1p1Ai = t17+t36;
		a1p1Br = t16+t37;		a1p1Bi = t17-t36;
	/*...Block 8: t0E,t1E,t2E,t3E	*/
		rt =t1E*c + t1F*s;		it =t1F*c - t1E*s;
		t1E=t0E+rt;			t0E =t0E-rt;
		t1F=t0F+it;			t0F =t0F-it;
		rt =t2E*s32_1 - t2F*c32_1;	t2F=t2F*s32_1 + t2E*c32_1;	t2E=rt;
		rt =t3E*s32_3 - t3F*c32_3;	it =t3F*s32_3 + t3E*c32_3;
		t3E=t2E+rt;			t2E=t2E-rt;
		t3F=t2F+it;			t2F=t2F-it;
		a1p1Cr = t0E+t2E;		a1p1Ci = t0F+t2F;
		a1p1Dr = t0E-t2E;		a1p1Di = t0F-t2F;
		a1p1Er = t1E-t3F;		a1p1Ei = t1F+t3E;
		a1p1Fr = t1E+t3F;		a1p1Fi = t1F-t3E;
	/*************************************************************/
	/*                  2nd set of inputs:                       */
	/*************************************************************/
		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;
	/*...Block 1:	*/
		t00=a[rdum   ];				t01=a[idum   ];
		rt =a[rdum+32]*cB10 - a[idum+32]*sB10;	it =a[idum+32]*cB10 + a[rdum+32]*sB10;
		t02=t00-rt;		t03=t01-it;
		t00=t00+rt;		t01=t01+it;
		t04=a[rdum+16]*cB08 - a[idum+16]*sB08;	t05=a[idum+16]*cB08 + a[rdum+16]*sB08;
		rt =a[rdum+48]*cB18 - a[idum+48]*sB18;	it =a[idum+48]*cB18 + a[rdum+48]*sB18;
		t06=t04-rt;		t07=t05-it;
		t04=t04+rt;		t05=t05+it;
		rt =t04;		it =t05;
		t04=t00-rt;		t05=t01-it;
		t00=t00+rt;		t01=t01+it;
		rt =t06;		it =t07;
		t06=t02+it;		t07=t03-rt;
		t02=t02-it;		t03=t03+rt;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;	// r24 -> r4 is -= 40 in scalar-double/SSE2/AVX modes, -= 44 in AVX-512 mode
	#endif
		t08=a[rdum+8 ]*cB04 - a[idum+8 ]*sB04;	t09=a[idum+8 ]*cB04 + a[rdum+8 ]*sB04;
		rt =a[rdum+40]*cB14 - a[idum+40]*sB14;	it =a[idum+40]*cB14 + a[rdum+40]*sB14;
		t0A=t08-rt;		t0B=t09-it;
		t08=t08+rt;		t09=t09+it;
		t0C=a[rdum+24]*cB0C - a[idum+24]*sB0C;	t0D=a[idum+24]*cB0C + a[rdum+24]*sB0C;
		rt =a[rdum+56]*cB1C - a[idum+56]*sB1C;	it =a[idum+56]*cB1C + a[rdum+56]*sB1C;
		t0E=t0C-rt;		t0F=t0D-it;
		t0C=t0C+rt;		t0D=t0D+it;
		rt =t0C;		it =t0D;
		t0C=t08-rt;		t0D=t09-it;
		t08=t08+rt;		t09=t09+it;
		rt =t0E;		it =t0F;
		t0E=t0A+it;		t0F=t0B-rt;
		t0A=t0A-it;		t0B=t0B+rt;
		rt =t08;		it =t09;
		t08=t00-rt;		t09=t01-it;
		t00=t00+rt;		t01=t01+it;
		rt =t0C;		it =t0D;
		t0C=t04+it;		t0D=t05-rt;
		t04=t04-it;		t05=t05+rt;
		rt =(t0A-t0B)*ISRT2;it =(t0A+t0B)*ISRT2;
		t0A=t02-rt;		t0B=t03-it;
		t02=t02+rt;		t03=t03+it;
		rt =(t0E+t0F)*ISRT2;it =(t0F-t0E)*ISRT2;
		t0E=t06+rt;		t0F=t07+it;
		t06=t06-rt;		t07=t07-it;
	/*...Block 2:	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;	// r28 -> r2 jump in AVX-512 is 4 fewer leftward than in AVX, hence idx-adjustment here is AVX value += 4
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;	// r28 -> r2 jump same in scalar-double and SSE2 modes, AVX is 2 more places leftward
	#endif
		t10=a[rdum+4 ]*cB02 - a[idum+4 ]*sB02;	t11=a[idum+4 ]*cB02 + a[rdum+4 ]*sB02;
		rt =a[rdum+36]*cB12 - a[idum+36]*sB12;	it =a[idum+36]*cB12 + a[rdum+36]*sB12;
		t12=t10-rt;		t13=t11-it;
		t10=t10+rt;		t11=t11+it;
		t14=a[rdum+20]*cB0A - a[idum+20]*sB0A;	t15=a[idum+20]*cB0A + a[rdum+20]*sB0A;
		rt =a[rdum+52]*cB1A - a[idum+52]*sB1A;	it =a[idum+52]*cB1A + a[rdum+52]*sB1A;
		t16=t14-rt;		t17=t15-it;
		t14=t14+rt;		t15=t15+it;
		rt =t14;		it =t15;
		t14=t10-rt;		t15=t11-it;
		t10=t10+rt;		t11=t11+it;
		rt =t16;		it =t17;
		t16=t12+it;		t17=t13-rt;
		t12=t12-it;		t13=t13+rt;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;	// r24 -> r4 jump same in scalar-double/SSE2/AVX modes, 4 more places leftward in AVX-512 mode
	#endif
		t18=a[rdum+12]*cB06 - a[idum+12]*sB06;	t19=a[idum+12]*cB06 + a[rdum+12]*sB06;
		rt =a[rdum+44]*cB16 - a[idum+44]*sB16;	it =a[idum+44]*cB16 + a[rdum+44]*sB16;
		t1A=t18-rt;		t1B=t19-it;
		t18=t18+rt;		t19=t19+it;
		t1C=a[rdum+28]*cB0E - a[idum+28]*sB0E;	t1D=a[idum+28]*cB0E + a[rdum+28]*sB0E;
		rt =a[rdum+60]*cB1E - a[idum+60]*sB1E;	it =a[idum+60]*cB1E + a[rdum+60]*sB1E;
		t1E=t1C-rt;		t1F=t1D-it;
		t1C=t1C+rt;		t1D=t1D+it;
		rt =t1C;		it =t1D;
		t1C=t18-rt;		t1D=t19-it;
		t18=t18+rt;		t19=t19+it;
		rt =t1E;		it =t1F;
		t1E=t1A+it;		t1F=t1B-rt;
		t1A=t1A-it;		t1B=t1B+rt;
		rt =t18;		it =t19;
		t18=t10-rt;		t19=t11-it;
		t10=t10+rt;		t11=t11+it;
		rt =t1C;		it =t1D;
		t1C=t14+it;		t1D=t15-rt;
		t14=t14-it;		t15=t15+rt;
		rt =(t1A-t1B)*ISRT2;it =(t1A+t1B)*ISRT2;
		t1A=t12-rt;		t1B=t13-it;
		t12=t12+rt;		t13=t13+it;
		rt =(t1E+t1F)*ISRT2;it =(t1F-t1E)*ISRT2;
		t1E=t16+rt;		t1F=t17+it;
		t16=t16-rt;		t17=t17-it;
	/*...Block 3:	*/
	#ifdef USE_AVX512
		rdum += 5;	idum += 5;	// r30 -> r1 jump in AVX-512 is 4 fewer leftward than in AVX, hence idx-adjustment here is AVX value += 4
	#elif defined(USE_AVX)
		++rdum;	++idum;			// r30 -> r1 jump same in AVX mode is 2 fewer places leftward than SSE2
	#elif defined(USE_SSE2)
		--rdum;	--idum;			// r30 -> r1 jump in SSE2 mode is 1 more places leftward than in scalar-double mode
	#endif
		t20=a[rdum+2 ]*cB01 - a[idum+2 ]*sB01;	t21=a[idum+2 ]*cB01 + a[rdum+2 ]*sB01;
		rt =a[rdum+34]*cB11 - a[idum+34]*sB11;	it =a[idum+34]*cB11 + a[rdum+34]*sB11;
		t22=t20-rt;		t23=t21-it;
		t20=t20+rt;		t21=t21+it;
		t24=a[rdum+18]*cB09 - a[idum+18]*sB09;	t25=a[idum+18]*cB09 + a[rdum+18]*sB09;
		rt =a[rdum+50]*cB19 - a[idum+50]*sB19;	it =a[idum+50]*cB19 + a[rdum+50]*sB19;
		t26=t24-rt;		t27=t25-it;
		t24=t24+rt;		t25=t25+it;
		rt =t24;		it =t25;
		t24=t20-rt;		t25=t21-it;
		t20=t20+rt;		t21=t21+it;
		rt =t26;		it =t27;
		t26=t22+it;		t27=t23-rt;
		t22=t22-it;		t23=t23+rt;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t28=a[rdum+10]*cB05 - a[idum+10]*sB05;	t29=a[idum+10]*cB05 + a[rdum+10]*sB05;
		rt =a[rdum+42]*cB15 - a[idum+42]*sB15;	it =a[idum+42]*cB15 + a[rdum+42]*sB15;
		t2A=t28-rt;		t2B=t29-it;
		t28=t28+rt;		t29=t29+it;
		t2C=a[rdum+26]*cB0D - a[idum+26]*sB0D;	t2D=a[idum+26]*cB0D + a[rdum+26]*sB0D;
		rt =a[rdum+58]*cB1D - a[idum+58]*sB1D;	it =a[idum+58]*cB1D + a[rdum+58]*sB1D;
		t2E=t2C-rt;		t2F=t2D-it;
		t2C=t2C+rt;		t2D=t2D+it;
		rt =t2C;		it =t2D;
		t2C=t28-rt;		t2D=t29-it;
		t28=t28+rt;		t29=t29+it;
		rt =t2E;		it =t2F;
		t2E=t2A+it;		t2F=t2B-rt;
		t2A=t2A-it;		t2B=t2B+rt;
		rt =t28;		it =t29;
		t28=t20-rt;		t29=t21-it;
		t20=t20+rt;		t21=t21+it;
		rt =t2C;		it =t2D;
		t2C=t24+it;		t2D=t25-rt;
		t24=t24-it;		t25=t25+rt;
		rt =(t2A-t2B)*ISRT2;it =(t2A+t2B)*ISRT2;
		t2A=t22-rt;		t2B=t23-it;
		t22=t22+rt;		t23=t23+it;
		rt =(t2E+t2F)*ISRT2;it =(t2F-t2E)*ISRT2;
		t2E=t26+rt;		t2F=t27+it;
		t26=t26-rt;		t27=t27-it;
	/*...Block 4:	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;	// r29 -> r3 jump in AVX-512 is 4 fewer leftward than in AVX, hence idx-adjustment here is AVX value += 4
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;	// r29 -> r3 jump same in scalar-double and SSE2 modes, AVX is 2 more places leftward
	#endif
		t30=a[rdum+6 ]*cB03 - a[idum+6 ]*sB03;	t31=a[idum+6 ]*cB03 + a[rdum+6 ]*sB03;
		rt =a[rdum+38]*cB13 - a[idum+38]*sB13;	it =a[idum+38]*cB13 + a[rdum+38]*sB13;
		t32=t30-rt;		t33=t31-it;
		t30=t30+rt;		t31=t31+it;
		t34=a[rdum+22]*cB0B - a[idum+22]*sB0B;	t35=a[idum+22]*cB0B + a[rdum+22]*sB0B;
		rt =a[rdum+54]*cB1B - a[idum+54]*sB1B;	it =a[idum+54]*cB1B + a[rdum+54]*sB1B;
		t36=t34-rt;		t37=t35-it;
		t34=t34+rt;		t35=t35+it;
		rt =t34;		it =t35;
		t34=t30-rt;		t35=t31-it;
		t30=t30+rt;		t31=t31+it;
		rt =t36;		it =t37;
		t36=t32+it;		t37=t33-rt;
		t32=t32-it;		t33=t33+rt;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t38=a[rdum+14]*cB07 - a[idum+14]*sB07;	t39=a[idum+14]*cB07 + a[rdum+14]*sB07;
		rt =a[rdum+46]*cB17 - a[idum+46]*sB17;	it =a[idum+46]*cB17 + a[rdum+46]*sB17;
		t3A=t38-rt;		t3B=t39-it;
		t38=t38+rt;		t39=t39+it;
		t3C=a[rdum+30]*cB0F - a[idum+30]*sB0F;	t3D=a[idum+30]*cB0F + a[rdum+30]*sB0F;
		rt =a[rdum+62]*cB1F - a[idum+62]*sB1F;	it =a[idum+62]*cB1F + a[rdum+62]*sB1F;
		t3E=t3C-rt;		t3F=t3D-it;
		t3C=t3C+rt;		t3D=t3D+it;
		rt =t3C;		it =t3D;
		t3C=t38-rt;		t3D=t39-it;
		t38=t38+rt;		t39=t39+it;
		rt =t3E;		it =t3F;
		t3E=t3A+it;		t3F=t3B-rt;
		t3A=t3A-it;		t3B=t3B+rt;
		rt =t38;		it =t39;
		t38=t30-rt;		t39=t31-it;
		t30=t30+rt;		t31=t31+it;
		rt =t3C;		it =t3D;
		t3C=t34+it;		t3D=t35-rt;
		t34=t34-it;		t35=t35+rt;
		rt =(t3A-t3B)*ISRT2;it =(t3A+t3B)*ISRT2;
		t3A=t32-rt;		t3B=t33-it;
		t32=t32+rt;		t33=t33+it;
		rt =(t3E+t3F)*ISRT2;it =(t3F-t3E)*ISRT2;
		t3E=t36+rt;		t3F=t37+it;
		t36=t36-rt;		t37=t37-it;
	//...and now do eight radix-4 transforms, including the internal twiddle factors:
	/*...Block 1: t00,t10,t20,t30	*/
		rt =t10;	t10=t00-rt;	t00=t00+rt;
		it =t11;	t11=t01-it;	t01=t01+it;
		rt =t30;	t30=t20-rt;	t20=t20+rt;
		it =t31;	t31=t21-it;	t21=t21+it;
		a2p00r = t00+t20;		a2p00i = t01+t21;
		a2p01r = t00-t20;		a2p01i = t01-t21;
		a2p02r = t10-t31;		a2p02i = t11+t30;
		a2p03r = t10+t31;		a2p03i = t11-t30;
	/*...Block 5: t08,t18,t28,t38	*/
		rt =t18;	t18=t08+t19;	t08=t08-t19;
			t19=t09-rt;	t09=t09+rt;
		rt =(t28-t29)*ISRT2;	t29=(t28+t29)*ISRT2;		t28=rt;
		rt =(t39+t38)*ISRT2;	it =(t39-t38)*ISRT2;
		t38=t28+rt;			t28=t28-rt;
		t39=t29+it;			t29=t29-it;
		a2p04r = t08+t28;		a2p04i = t09+t29;
		a2p05r = t08-t28;		a2p05i = t09-t29;
		a2p06r = t18-t39;		a2p06i = t19+t38;
		a2p07r = t18+t39;		a2p07i = t19-t38;
	/*...Block 3: t04,t14,t24,t34	*/
		rt =(t14-t15)*ISRT2;	it =(t14+t15)*ISRT2;
		t14=t04-rt;			t04=t04+rt;
		t15=t05-it;			t05=t05+it;
		rt =t24*c - t25*s;		t25=t25*c + t24*s;		t24=rt;
		rt =t34*s - t35*c;		it =t35*s + t34*c;
		t34=t24-rt;			t24=t24+rt;
		t35=t25-it;			t25=t25+it;
		a2p08r = t04+t24;		a2p08i = t05+t25;
		a2p09r = t04-t24;		a2p09i = t05-t25;
		a2p0Ar = t14-t35;		a2p0Ai = t15+t34;
		a2p0Br = t14+t35;		a2p0Bi = t15-t34;
	/*...Block 7: t0C,t1C,t2C,t3C	*/
		rt =(t1D+t1C)*ISRT2;	it =(t1D-t1C)*ISRT2;
		t1C=t0C+rt;			t0C=t0C-rt;
		t1D=t0D+it;			t0D=t0D-it;
		rt =t2C*s - t2D*c;		t2D=t2D*s + t2C*c;		t2C=rt;
		rt =t3C*c - t3D*s;		it =t3D*c + t3C*s;
		t3C=t2C+rt;			t2C=t2C-rt;
		t3D=t2D+it;			t2D=t2D-it;
		a2p0Cr = t0C+t2C;		a2p0Ci = t0D+t2D;
		a2p0Dr = t0C-t2C;		a2p0Di = t0D-t2D;
		a2p0Er = t1C-t3D;		a2p0Ei = t1D+t3C;
		a2p0Fr = t1C+t3D;		a2p0Fi = t1D-t3C;
	/*...Block 2: t02,t12,t22,t32	*/
		rt =t12*c - t13*s;		it =t13*c + t12*s;
		t12=t02-rt;			t02=t02+rt;
		t13=t03-it;			t03=t03+it;
		rt =t22*c32_1 - t23*s32_1;	t23=t23*c32_1 + t22*s32_1;	t22=rt;
		rt =t32*c32_3 - t33*s32_3;	it =t33*c32_3 + t32*s32_3;
		t32=t22-rt;			t22=t22+rt;
		t33=t23-it;			t23=t23+it;
		a2p10r = t02+t22;		a2p10i = t03+t23;
		a2p11r = t02-t22;		a2p11i = t03-t23;
		a2p12r = t12-t33;		a2p12i = t13+t32;
		a2p13r = t12+t33;		a2p13i = t13-t32;
	/*...Block 6: t0A,t1A,t2A,t3A	*/
		rt =t1A*s + t1B*c;		it =t1B*s - t1A*c;
		t1A=t0A+rt;			t0A =t0A-rt;
		t1B=t0B+it;			t0B =t0B-it;
		rt =t2A*s32_3 - t2B*c32_3;	t2B=t2B*s32_3 + t2A*c32_3;	t2A=rt;
		rt =t3A*c32_1 + t3B*s32_1;	it =t3B*c32_1 - t3A*s32_1;
		t3A=t2A+rt;			t2A=t2A-rt;
		t3B=t2B+it;			t2B=t2B-it;
		a2p14r = t0A+t2A;		a2p14i = t0B+t2B;
		a2p15r = t0A-t2A;		a2p15i = t0B-t2B;
		a2p16r = t1A-t3B;		a2p16i = t1B+t3A;
		a2p17r = t1A+t3B;		a2p17i = t1B-t3A;
	/*...Block 4: t06,t16,t26,t36	*/
		rt =t16*s - t17*c;		it =t17*s + t16*c;
		t16=t06-rt;			t06 =t06+rt;
		t17=t07-it;			t07 =t07+it;
		rt =t26*c32_3 - t27*s32_3;	t27=t27*c32_3 + t26*s32_3;	t26=rt;
		rt =t36*s32_1 + t37*c32_1;	it =t37*s32_1 - t36*c32_1;
		t36=t26+rt;			t26=t26-rt;
		t37=t27+it;			t27=t27-it;
		a2p18r = t06+t26;		a2p18i = t07+t27;
		a2p19r = t06-t26;		a2p19i = t07-t27;
		a2p1Ar = t16-t37;		a2p1Ai = t17+t36;
		a2p1Br = t16+t37;		a2p1Bi = t17-t36;
	/*...Block 8: t0E,t1E,t2E,t3E	*/
		rt =t1E*c + t1F*s;		it =t1F*c - t1E*s;
		t1E=t0E+rt;			t0E =t0E-rt;
		t1F=t0F+it;			t0F =t0F-it;
		rt =t2E*s32_1 - t2F*c32_1;	t2F=t2F*s32_1 + t2E*c32_1;	t2E=rt;
		rt =t3E*s32_3 - t3F*c32_3;	it =t3F*s32_3 + t3E*c32_3;
		t3E=t2E+rt;			t2E=t2E-rt;
		t3F=t2F+it;			t2F=t2F-it;
		a2p1Cr = t0E+t2E;		a2p1Ci = t0F+t2F;
		a2p1Dr = t0E-t2E;		a2p1Di = t0F-t2F;
		a2p1Er = t1E-t3F;		a2p1Ei = t1F+t3E;
		a2p1Fr = t1E+t3F;		a2p1Fi = t1F-t3E;

skip_fwd_fft:	// v20: jump-to-point for both-inputs-already-fwd-FFTed case

	// v19: If fwd_fft_only = 1, write fwd-FFT result back to input array, skipping dyadic-square and inv-FFT steps:
	if(fwd_fft_only == 1)
	{
	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;
	/*...Block 1: t00,t10,t20,t30	*/
		a[rdum   ] = a1p00r;	a[idum   ] = a1p00i;
		a[rdum+32] = a1p01r;	a[idum+32] = a1p01i;
		a[rdum+16] = a1p02r;	a[idum+16] = a1p02i;
		a[rdum+48] = a1p03r;	a[idum+48] = a1p03i;
	/*...Block 5: t08,t18,t28,t38	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		a[rdum+8 ] = a1p04r;	a[idum+8 ] = a1p04i;
		a[rdum+40] = a1p05r;	a[idum+40] = a1p05i;
		a[rdum+24] = a1p06r;	a[idum+24] = a1p06i;
		a[rdum+56] = a1p07r;	a[idum+56] = a1p07i;
	/*...Block 3: t04,t14,t24,t34	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
		a[rdum+4 ] = a1p08r;	a[idum+4 ] = a1p08i;
		a[rdum+36] = a1p09r;	a[idum+36] = a1p09i;
		a[rdum+20] = a1p0Ar;	a[idum+20] = a1p0Ai;
		a[rdum+52] = a1p0Br;	a[idum+52] = a1p0Bi;
	/*...Block 7: t0C,t1C,t2C,t3C	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		a[rdum+12] = a1p0Cr;	a[idum+12] = a1p0Ci;
		a[rdum+44] = a1p0Dr;	a[idum+44] = a1p0Di;
		a[rdum+28] = a1p0Er;	a[idum+28] = a1p0Ei;
		a[rdum+60] = a1p0Fr;	a[idum+60] = a1p0Fi;
	/*...Block 2: t02,t12,t22,t32	*/
	#ifdef USE_AVX512
		rdum += 5;	idum += 5;
	#elif defined(USE_AVX)
		++rdum;	++idum;
	#elif defined(USE_SSE2)
		--rdum;	--idum;
	#endif
		a[rdum+2 ] = a1p10r;	a[idum+2 ] = a1p10i;
		a[rdum+34] = a1p11r;	a[idum+34] = a1p11i;
		a[rdum+18] = a1p12r;	a[idum+18] = a1p12i;
		a[rdum+50] = a1p13r;	a[idum+50] = a1p13i;
	/*...Block 6: t0A,t1A,t2A,t3A	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		a[rdum+10] = a1p14r;	a[idum+10] = a1p14i;
		a[rdum+42] = a1p15r;	a[idum+42] = a1p15i;
		a[rdum+26] = a1p16r;	a[idum+26] = a1p16i;
		a[rdum+58] = a1p17r;	a[idum+58] = a1p17i;
	/*...Block 4: t06,t16,t26,t36	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
		a[rdum+6 ] = a1p18r;	a[idum+6 ] = a1p18i;
		a[rdum+38] = a1p19r;	a[idum+38] = a1p19i;
		a[rdum+22] = a1p1Ar;	a[idum+22] = a1p1Ai;
		a[rdum+54] = a1p1Br;	a[idum+54] = a1p1Bi;
	/*...Block 8: t0E,t1E,t2E,t3E	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		a[rdum+14] = a1p1Cr;	a[idum+14] = a1p1Ci;
		a[rdum+46] = a1p1Dr;	a[idum+46] = a1p1Di;
		a[rdum+30] = a1p1Er;	a[idum+30] = a1p1Ei;
		a[rdum+62] = a1p1Fr;	a[idum+62] = a1p1Fi;
	/*************************************************************/
	/*                  2nd set of inputs:                       */
	/*************************************************************/
		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;
	/*...Block 1: t00,t10,t20,t30	*/
		a[rdum   ] = a2p00r;	a[idum   ] = a2p00i;
		a[rdum+32] = a2p01r;	a[idum+32] = a2p01i;
		a[rdum+16] = a2p02r;	a[idum+16] = a2p02i;
		a[rdum+48] = a2p03r;	a[idum+48] = a2p03i;
	/*...Block 5: t08,t18,t28,t38	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		a[rdum+8 ] = a2p04r;	a[idum+8 ] = a2p04i;
		a[rdum+40] = a2p05r;	a[idum+40] = a2p05i;
		a[rdum+24] = a2p06r;	a[idum+24] = a2p06i;
		a[rdum+56] = a2p07r;	a[idum+56] = a2p07i;
	/*...Block 3: t04,t14,t24,t34	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
		a[rdum+4 ] = a2p08r;	a[idum+4 ] = a2p08i;
		a[rdum+36] = a2p09r;	a[idum+36] = a2p09i;
		a[rdum+20] = a2p0Ar;	a[idum+20] = a2p0Ai;
		a[rdum+52] = a2p0Br;	a[idum+52] = a2p0Bi;
	/*...Block 7: t0C,t1C,t2C,t3C	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		a[rdum+12] = a2p0Cr;	a[idum+12] = a2p0Ci;
		a[rdum+44] = a2p0Dr;	a[idum+44] = a2p0Di;
		a[rdum+28] = a2p0Er;	a[idum+28] = a2p0Ei;
		a[rdum+60] = a2p0Fr;	a[idum+60] = a2p0Fi;
	/*...Block 2: t02,t12,t22,t32	*/
	#ifdef USE_AVX512
		rdum += 5;	idum += 5;
	#elif defined(USE_AVX)
		++rdum;	++idum;
	#elif defined(USE_SSE2)
		--rdum;	--idum;
	#endif
		a[rdum+2 ] = a2p10r;	a[idum+2 ] = a2p10i;
		a[rdum+34] = a2p11r;	a[idum+34] = a2p11i;
		a[rdum+18] = a2p12r;	a[idum+18] = a2p12i;
		a[rdum+50] = a2p13r;	a[idum+50] = a2p13i;
	/*...Block 6: t0A,t1A,t2A,t3A	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		a[rdum+10] = a2p14r;	a[idum+10] = a2p14i;
		a[rdum+42] = a2p15r;	a[idum+42] = a2p15i;
		a[rdum+26] = a2p16r;	a[idum+26] = a2p16i;
		a[rdum+58] = a2p17r;	a[idum+58] = a2p17i;
	/*...Block 4: t06,t16,t26,t36	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
		a[rdum+6 ] = a2p18r;	a[idum+6 ] = a2p18i;
		a[rdum+38] = a2p19r;	a[idum+38] = a2p19i;
		a[rdum+22] = a2p1Ar;	a[idum+22] = a2p1Ai;
		a[rdum+54] = a2p1Br;	a[idum+54] = a2p1Bi;
	/*...Block 8: t0E,t1E,t2E,t3E	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		a[rdum+14] = a2p1Cr;	a[idum+14] = a2p1Ci;
		a[rdum+46] = a2p1Dr;	a[idum+46] = a2p1Di;
		a[rdum+30] = a2p1Er;	a[idum+30] = a2p1Ei;
		a[rdum+62] = a2p1Fr;	a[idum+62] = a2p1Fi;
		goto loop;	// Skip dyadic-mul and iFFT

	} else if (fwd_fft_only) {	// 2-input modmul: fwdFFT data dyadic-mul'ed with precomputed 2nd-vector stored in fwd-FFTed form in b[]:

	  if(fwd_fft_only == 3)	// v20: Both inputs enter fwd-FFTed, must copy from main arrays a[],b[] to local-vars
	  {
	  /*************************************************************/
	  /*                  1st set of inputs:                       */
	  /*************************************************************/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;
	  /*...Block 1: t00,t10,t20,t30	*/
		a1p00r = a[rdum   ];	a1p00i = a[idum   ];
		a1p01r = a[rdum+32];	a1p01i = a[idum+32];
		a1p02r = a[rdum+16];	a1p02i = a[idum+16];
		a1p03r = a[rdum+48];	a1p03i = a[idum+48];
	  /*...Block 5: t08,t18,t28,t38	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		a1p04r = a[rdum+8 ];	a1p04i = a[idum+8 ];
		a1p05r = a[rdum+40];	a1p05i = a[idum+40];
		a1p06r = a[rdum+24];	a1p06i = a[idum+24];
		a1p07r = a[rdum+56];	a1p07i = a[idum+56];
	  /*...Block 3: t04,t14,t24,t34	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		a1p08r = a[rdum+4 ];	a1p08i = a[idum+4 ];
		a1p09r = a[rdum+36];	a1p09i = a[idum+36];
		a1p0Ar = a[rdum+20];	a1p0Ai = a[idum+20];
		a1p0Br = a[rdum+52];	a1p0Bi = a[idum+52];
	  /*...Block 7: t0C,t1C,t2C,t3C	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		a1p0Cr = a[rdum+12];	a1p0Ci = a[idum+12];
		a1p0Dr = a[rdum+44];	a1p0Di = a[idum+44];
		a1p0Er = a[rdum+28];	a1p0Ei = a[idum+28];
		a1p0Fr = a[rdum+60];	a1p0Fi = a[idum+60];
	  /*...Block 2: t02,t12,t22,t32	*/
	  #ifdef USE_AVX512
		rdum += 5;	idum += 5;
	  #elif defined(USE_AVX)
		++rdum;	++idum;
	  #elif defined(USE_SSE2)
		--rdum;	--idum;
	  #endif
		a1p10r = a[rdum+2 ];	a1p10i = a[idum+2 ];
		a1p11r = a[rdum+34];	a1p11i = a[idum+34];
		a1p12r = a[rdum+18];	a1p12i = a[idum+18];
		a1p13r = a[rdum+50];	a1p13i = a[idum+50];
	  /*...Block 6: t0A,t1A,t2A,t3A	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		a1p14r = a[rdum+10];	a1p14i = a[idum+10];
		a1p15r = a[rdum+42];	a1p15i = a[idum+42];
		a1p16r = a[rdum+26];	a1p16i = a[idum+26];
		a1p17r = a[rdum+58];	a1p17i = a[idum+58];
	  /*...Block 4: t06,t16,t26,t36	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		a1p18r = a[rdum+6 ];	a1p18i = a[idum+6 ];
		a1p19r = a[rdum+38];	a1p19i = a[idum+38];
		a1p1Ar = a[rdum+22];	a1p1Ai = a[idum+22];
		a1p1Br = a[rdum+54];	a1p1Bi = a[idum+54];
	  /*...Block 8: t0E,t1E,t2E,t3E	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		a1p1Cr = a[rdum+14];	a1p1Ci = a[idum+14];
		a1p1Dr = a[rdum+46];	a1p1Di = a[idum+46];
		a1p1Er = a[rdum+30];	a1p1Ei = a[idum+30];
		a1p1Fr = a[rdum+62];	a1p1Fi = a[idum+62];
	  /*************************************************************/
	  /*                  2nd set of inputs:                       */
	  /*************************************************************/
		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;
	  /*...Block 1: t00,t10,t20,t30	*/
		a2p00r = a[rdum   ];	a2p00i = a[idum   ];
		a2p01r = a[rdum+32];	a2p01i = a[idum+32];
		a2p02r = a[rdum+16];	a2p02i = a[idum+16];
		a2p03r = a[rdum+48];	a2p03i = a[idum+48];
	  /*...Block 5: t08,t18,t28,t38	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		a2p04r = a[rdum+8 ];	a2p04i = a[idum+8 ];
		a2p05r = a[rdum+40];	a2p05i = a[idum+40];
		a2p06r = a[rdum+24];	a2p06i = a[idum+24];
		a2p07r = a[rdum+56];	a2p07i = a[idum+56];
	  /*...Block 3: t04,t14,t24,t34	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		a2p08r = a[rdum+4 ];	a2p08i = a[idum+4 ];
		a2p09r = a[rdum+36];	a2p09i = a[idum+36];
		a2p0Ar = a[rdum+20];	a2p0Ai = a[idum+20];
		a2p0Br = a[rdum+52];	a2p0Bi = a[idum+52];
	  /*...Block 7: t0C,t1C,t2C,t3C	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		a2p0Cr = a[rdum+12];	a2p0Ci = a[idum+12];
		a2p0Dr = a[rdum+44];	a2p0Di = a[idum+44];
		a2p0Er = a[rdum+28];	a2p0Ei = a[idum+28];
		a2p0Fr = a[rdum+60];	a2p0Fi = a[idum+60];
	  /*...Block 2: t02,t12,t22,t32	*/
	  #ifdef USE_AVX512
		rdum += 5;	idum += 5;
	  #elif defined(USE_AVX)
		++rdum;	++idum;
	  #elif defined(USE_SSE2)
		--rdum;	--idum;
	  #endif
		a2p10r = a[rdum+2 ];	a2p10i = a[idum+2 ];
		a2p11r = a[rdum+34];	a2p11i = a[idum+34];
		a2p12r = a[rdum+18];	a2p12i = a[idum+18];
		a2p13r = a[rdum+50];	a2p13i = a[idum+50];
	  /*...Block 6: t0A,t1A,t2A,t3A	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		a2p14r = a[rdum+10];	a2p14i = a[idum+10];
		a2p15r = a[rdum+42];	a2p15i = a[idum+42];
		a2p16r = a[rdum+26];	a2p16i = a[idum+26];
		a2p17r = a[rdum+58];	a2p17i = a[idum+58];
	  /*...Block 4: t06,t16,t26,t36	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		a2p18r = a[rdum+6 ];	a2p18i = a[idum+6 ];
		a2p19r = a[rdum+38];	a2p19i = a[idum+38];
		a2p1Ar = a[rdum+22];	a2p1Ai = a[idum+22];
		a2p1Br = a[rdum+54];	a2p1Bi = a[idum+54];
	  /*...Block 8: t0E,t1E,t2E,t3E	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		a2p1Cr = a[rdum+14];	a2p1Ci = a[idum+14];
		a2p1Dr = a[rdum+46];	a2p1Di = a[idum+46];
		a2p1Er = a[rdum+30];	a2p1Ei = a[idum+30];
		a2p1Fr = a[rdum+62];	a2p1Fi = a[idum+62];
	  }

	  if(c_arr) {	// c_arr != 0x0: a * (b - c)

	  /*************************************************************/
	  /*                  1st set of inputs:                       */
	  /*************************************************************/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;
	  /*...Block 1: t00,t10,t20,t30	*/
		b1p00r = b[rdum   ] - c_arr[rdum   ];	b1p00i = b[idum   ] - c_arr[idum   ];
		b1p01r = b[rdum+32] - c_arr[rdum+32];	b1p01i = b[idum+32] - c_arr[idum+32];
		b1p02r = b[rdum+16] - c_arr[rdum+16];	b1p02i = b[idum+16] - c_arr[idum+16];
		b1p03r = b[rdum+48] - c_arr[rdum+48];	b1p03i = b[idum+48] - c_arr[idum+48];
	  /*...Block 5: t08,t18,t28,t38	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b1p04r = b[rdum+8 ] - c_arr[rdum+8 ];	b1p04i = b[idum+8 ] - c_arr[idum+8 ];
		b1p05r = b[rdum+40] - c_arr[rdum+40];	b1p05i = b[idum+40] - c_arr[idum+40];
		b1p06r = b[rdum+24] - c_arr[rdum+24];	b1p06i = b[idum+24] - c_arr[idum+24];
		b1p07r = b[rdum+56] - c_arr[rdum+56];	b1p07i = b[idum+56] - c_arr[idum+56];
	  /*...Block 3: t04,t14,t24,t34	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		b1p08r = b[rdum+4 ] - c_arr[rdum+4 ];	b1p08i = b[idum+4 ] - c_arr[idum+4 ];
		b1p09r = b[rdum+36] - c_arr[rdum+36];	b1p09i = b[idum+36] - c_arr[idum+36];
		b1p0Ar = b[rdum+20] - c_arr[rdum+20];	b1p0Ai = b[idum+20] - c_arr[idum+20];
		b1p0Br = b[rdum+52] - c_arr[rdum+52];	b1p0Bi = b[idum+52] - c_arr[idum+52];
	  /*...Block 7: t0C,t1C,t2C,t3C	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b1p0Cr = b[rdum+12] - c_arr[rdum+12];	b1p0Ci = b[idum+12] - c_arr[idum+12];
		b1p0Dr = b[rdum+44] - c_arr[rdum+44];	b1p0Di = b[idum+44] - c_arr[idum+44];
		b1p0Er = b[rdum+28] - c_arr[rdum+28];	b1p0Ei = b[idum+28] - c_arr[idum+28];
		b1p0Fr = b[rdum+60] - c_arr[rdum+60];	b1p0Fi = b[idum+60] - c_arr[idum+60];
	  /*...Block 2: t02,t12,t22,t32	*/
	  #ifdef USE_AVX512
		rdum += 5;	idum += 5;
	  #elif defined(USE_AVX)
		++rdum;	++idum;
	  #elif defined(USE_SSE2)
		--rdum;	--idum;
	  #endif
		b1p10r = b[rdum+2 ] - c_arr[rdum+2 ];	b1p10i = b[idum+2 ] - c_arr[idum+2 ];
		b1p11r = b[rdum+34] - c_arr[rdum+34];	b1p11i = b[idum+34] - c_arr[idum+34];
		b1p12r = b[rdum+18] - c_arr[rdum+18];	b1p12i = b[idum+18] - c_arr[idum+18];
		b1p13r = b[rdum+50] - c_arr[rdum+50];	b1p13i = b[idum+50] - c_arr[idum+50];
	  /*...Block 6: t0A,t1A,t2A,t3A	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b1p14r = b[rdum+10] - c_arr[rdum+10];	b1p14i = b[idum+10] - c_arr[idum+10];
		b1p15r = b[rdum+42] - c_arr[rdum+42];	b1p15i = b[idum+42] - c_arr[idum+42];
		b1p16r = b[rdum+26] - c_arr[rdum+26];	b1p16i = b[idum+26] - c_arr[idum+26];
		b1p17r = b[rdum+58] - c_arr[rdum+58];	b1p17i = b[idum+58] - c_arr[idum+58];
	  /*...Block 4: t06,t16,t26,t36	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		b1p18r = b[rdum+6 ] - c_arr[rdum+6 ];	b1p18i = b[idum+6 ] - c_arr[idum+6 ];
		b1p19r = b[rdum+38] - c_arr[rdum+38];	b1p19i = b[idum+38] - c_arr[idum+38];
		b1p1Ar = b[rdum+22] - c_arr[rdum+22];	b1p1Ai = b[idum+22] - c_arr[idum+22];
		b1p1Br = b[rdum+54] - c_arr[rdum+54];	b1p1Bi = b[idum+54] - c_arr[idum+54];
	  /*...Block 8: t0E,t1E,t2E,t3E	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b1p1Cr = b[rdum+14] - c_arr[rdum+14];	b1p1Ci = b[idum+14] - c_arr[idum+14];
		b1p1Dr = b[rdum+46] - c_arr[rdum+46];	b1p1Di = b[idum+46] - c_arr[idum+46];
		b1p1Er = b[rdum+30] - c_arr[rdum+30];	b1p1Ei = b[idum+30] - c_arr[idum+30];
		b1p1Fr = b[rdum+62] - c_arr[rdum+62];	b1p1Fi = b[idum+62] - c_arr[idum+62];
	  /*************************************************************/
	  /*                  2nd set of inputs:                       */
	  /*************************************************************/
		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;
	  /*...Block 1: t00,t10,t20,t30	*/
		b2p00r = b[rdum   ] - c_arr[rdum   ];	b2p00i = b[idum   ] - c_arr[idum   ];
		b2p01r = b[rdum+32] - c_arr[rdum+32];	b2p01i = b[idum+32] - c_arr[idum+32];
		b2p02r = b[rdum+16] - c_arr[rdum+16];	b2p02i = b[idum+16] - c_arr[idum+16];
		b2p03r = b[rdum+48] - c_arr[rdum+48];	b2p03i = b[idum+48] - c_arr[idum+48];
	  /*...Block 5: t08,t18,t28,t38	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b2p04r = b[rdum+8 ] - c_arr[rdum+8 ];	b2p04i = b[idum+8 ] - c_arr[idum+8 ];
		b2p05r = b[rdum+40] - c_arr[rdum+40];	b2p05i = b[idum+40] - c_arr[idum+40];
		b2p06r = b[rdum+24] - c_arr[rdum+24];	b2p06i = b[idum+24] - c_arr[idum+24];
		b2p07r = b[rdum+56] - c_arr[rdum+56];	b2p07i = b[idum+56] - c_arr[idum+56];
	  /*...Block 3: t04,t14,t24,t34	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		b2p08r = b[rdum+4 ] - c_arr[rdum+4 ];	b2p08i = b[idum+4 ] - c_arr[idum+4 ];
		b2p09r = b[rdum+36] - c_arr[rdum+36];	b2p09i = b[idum+36] - c_arr[idum+36];
		b2p0Ar = b[rdum+20] - c_arr[rdum+20];	b2p0Ai = b[idum+20] - c_arr[idum+20];
		b2p0Br = b[rdum+52] - c_arr[rdum+52];	b2p0Bi = b[idum+52] - c_arr[idum+52];
	  /*...Block 7: t0C,t1C,t2C,t3C	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b2p0Cr = b[rdum+12] - c_arr[rdum+12];	b2p0Ci = b[idum+12] - c_arr[idum+12];
		b2p0Dr = b[rdum+44] - c_arr[rdum+44];	b2p0Di = b[idum+44] - c_arr[idum+44];
		b2p0Er = b[rdum+28] - c_arr[rdum+28];	b2p0Ei = b[idum+28] - c_arr[idum+28];
		b2p0Fr = b[rdum+60] - c_arr[rdum+60];	b2p0Fi = b[idum+60] - c_arr[idum+60];
	  /*...Block 2: t02,t12,t22,t32	*/
	  #ifdef USE_AVX512
		rdum += 5;	idum += 5;
	  #elif defined(USE_AVX)
		++rdum;	++idum;
	  #elif defined(USE_SSE2)
		--rdum;	--idum;
	  #endif
		b2p10r = b[rdum+2 ] - c_arr[rdum+2 ];	b2p10i = b[idum+2 ] - c_arr[idum+2 ];
		b2p11r = b[rdum+34] - c_arr[rdum+34];	b2p11i = b[idum+34] - c_arr[idum+34];
		b2p12r = b[rdum+18] - c_arr[rdum+18];	b2p12i = b[idum+18] - c_arr[idum+18];
		b2p13r = b[rdum+50] - c_arr[rdum+50];	b2p13i = b[idum+50] - c_arr[idum+50];
	  /*...Block 6: t0A,t1A,t2A,t3A	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b2p14r = b[rdum+10] - c_arr[rdum+10];	b2p14i = b[idum+10] - c_arr[idum+10];
		b2p15r = b[rdum+42] - c_arr[rdum+42];	b2p15i = b[idum+42] - c_arr[idum+42];
		b2p16r = b[rdum+26] - c_arr[rdum+26];	b2p16i = b[idum+26] - c_arr[idum+26];
		b2p17r = b[rdum+58] - c_arr[rdum+58];	b2p17i = b[idum+58] - c_arr[idum+58];
	  /*...Block 4: t06,t16,t26,t36	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		b2p18r = b[rdum+6 ] - c_arr[rdum+6 ];	b2p18i = b[idum+6 ] - c_arr[idum+6 ];
		b2p19r = b[rdum+38] - c_arr[rdum+38];	b2p19i = b[idum+38] - c_arr[idum+38];
		b2p1Ar = b[rdum+22] - c_arr[rdum+22];	b2p1Ai = b[idum+22] - c_arr[idum+22];
		b2p1Br = b[rdum+54] - c_arr[rdum+54];	b2p1Bi = b[idum+54] - c_arr[idum+54];
	  /*...Block 8: t0E,t1E,t2E,t3E	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b2p1Cr = b[rdum+14] - c_arr[rdum+14];	b2p1Ci = b[idum+14] - c_arr[idum+14];
		b2p1Dr = b[rdum+46] - c_arr[rdum+46];	b2p1Di = b[idum+46] - c_arr[idum+46];
		b2p1Er = b[rdum+30] - c_arr[rdum+30];	b2p1Ei = b[idum+30] - c_arr[idum+30];
		b2p1Fr = b[rdum+62] - c_arr[rdum+62];	b2p1Fi = b[idum+62] - c_arr[idum+62];

	  } else {	// c_arr = 0x0: a * b

	  /*************************************************************/
	  /*                  1st set of inputs:                       */
	  /*************************************************************/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;
	  /*...Block 1: t00,t10,t20,t30	*/
		b1p00r = b[rdum   ];	b1p00i = b[idum   ];
		b1p01r = b[rdum+32];	b1p01i = b[idum+32];
		b1p02r = b[rdum+16];	b1p02i = b[idum+16];
		b1p03r = b[rdum+48];	b1p03i = b[idum+48];
	  /*...Block 5: t08,t18,t28,t38	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b1p04r = b[rdum+8 ];	b1p04i = b[idum+8 ];
		b1p05r = b[rdum+40];	b1p05i = b[idum+40];
		b1p06r = b[rdum+24];	b1p06i = b[idum+24];
		b1p07r = b[rdum+56];	b1p07i = b[idum+56];
	  /*...Block 3: t04,t14,t24,t34	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		b1p08r = b[rdum+4 ];	b1p08i = b[idum+4 ];
		b1p09r = b[rdum+36];	b1p09i = b[idum+36];
		b1p0Ar = b[rdum+20];	b1p0Ai = b[idum+20];
		b1p0Br = b[rdum+52];	b1p0Bi = b[idum+52];
	  /*...Block 7: t0C,t1C,t2C,t3C	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b1p0Cr = b[rdum+12];	b1p0Ci = b[idum+12];
		b1p0Dr = b[rdum+44];	b1p0Di = b[idum+44];
		b1p0Er = b[rdum+28];	b1p0Ei = b[idum+28];
		b1p0Fr = b[rdum+60];	b1p0Fi = b[idum+60];
	  /*...Block 2: t02,t12,t22,t32	*/
	  #ifdef USE_AVX512
		rdum += 5;	idum += 5;
	  #elif defined(USE_AVX)
		++rdum;	++idum;
	  #elif defined(USE_SSE2)
		--rdum;	--idum;
	  #endif
		b1p10r = b[rdum+2 ];	b1p10i = b[idum+2 ];
		b1p11r = b[rdum+34];	b1p11i = b[idum+34];
		b1p12r = b[rdum+18];	b1p12i = b[idum+18];
		b1p13r = b[rdum+50];	b1p13i = b[idum+50];
	  /*...Block 6: t0A,t1A,t2A,t3A	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b1p14r = b[rdum+10];	b1p14i = b[idum+10];
		b1p15r = b[rdum+42];	b1p15i = b[idum+42];
		b1p16r = b[rdum+26];	b1p16i = b[idum+26];
		b1p17r = b[rdum+58];	b1p17i = b[idum+58];
	  /*...Block 4: t06,t16,t26,t36	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		b1p18r = b[rdum+6 ];	b1p18i = b[idum+6 ];
		b1p19r = b[rdum+38];	b1p19i = b[idum+38];
		b1p1Ar = b[rdum+22];	b1p1Ai = b[idum+22];
		b1p1Br = b[rdum+54];	b1p1Bi = b[idum+54];
	  /*...Block 8: t0E,t1E,t2E,t3E	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b1p1Cr = b[rdum+14];	b1p1Ci = b[idum+14];
		b1p1Dr = b[rdum+46];	b1p1Di = b[idum+46];
		b1p1Er = b[rdum+30];	b1p1Ei = b[idum+30];
		b1p1Fr = b[rdum+62];	b1p1Fi = b[idum+62];
	  /*************************************************************/
	  /*                  2nd set of inputs:                       */
	  /*************************************************************/
		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;
	  /*...Block 1: t00,t10,t20,t30	*/
		b2p00r = b[rdum   ];	b2p00i = b[idum   ];
		b2p01r = b[rdum+32];	b2p01i = b[idum+32];
		b2p02r = b[rdum+16];	b2p02i = b[idum+16];
		b2p03r = b[rdum+48];	b2p03i = b[idum+48];
	  /*...Block 5: t08,t18,t28,t38	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b2p04r = b[rdum+8 ];	b2p04i = b[idum+8 ];
		b2p05r = b[rdum+40];	b2p05i = b[idum+40];
		b2p06r = b[rdum+24];	b2p06i = b[idum+24];
		b2p07r = b[rdum+56];	b2p07i = b[idum+56];
	  /*...Block 3: t04,t14,t24,t34	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		b2p08r = b[rdum+4 ];	b2p08i = b[idum+4 ];
		b2p09r = b[rdum+36];	b2p09i = b[idum+36];
		b2p0Ar = b[rdum+20];	b2p0Ai = b[idum+20];
		b2p0Br = b[rdum+52];	b2p0Bi = b[idum+52];
	  /*...Block 7: t0C,t1C,t2C,t3C	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b2p0Cr = b[rdum+12];	b2p0Ci = b[idum+12];
		b2p0Dr = b[rdum+44];	b2p0Di = b[idum+44];
		b2p0Er = b[rdum+28];	b2p0Ei = b[idum+28];
		b2p0Fr = b[rdum+60];	b2p0Fi = b[idum+60];
	  /*...Block 2: t02,t12,t22,t32	*/
	  #ifdef USE_AVX512
		rdum += 5;	idum += 5;
	  #elif defined(USE_AVX)
		++rdum;	++idum;
	  #elif defined(USE_SSE2)
		--rdum;	--idum;
	  #endif
		b2p10r = b[rdum+2 ];	b2p10i = b[idum+2 ];
		b2p11r = b[rdum+34];	b2p11i = b[idum+34];
		b2p12r = b[rdum+18];	b2p12i = b[idum+18];
		b2p13r = b[rdum+50];	b2p13i = b[idum+50];
	  /*...Block 6: t0A,t1A,t2A,t3A	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b2p14r = b[rdum+10];	b2p14i = b[idum+10];
		b2p15r = b[rdum+42];	b2p15i = b[idum+42];
		b2p16r = b[rdum+26];	b2p16i = b[idum+26];
		b2p17r = b[rdum+58];	b2p17i = b[idum+58];
	  /*...Block 4: t06,t16,t26,t36	*/
	  #ifdef USE_AVX512
		rdum += 2;	idum += 2;
	  #elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	  #endif
		b2p18r = b[rdum+6 ];	b2p18i = b[idum+6 ];
		b2p19r = b[rdum+38];	b2p19i = b[idum+38];
		b2p1Ar = b[rdum+22];	b2p1Ai = b[idum+22];
		b2p1Br = b[rdum+54];	b2p1Bi = b[idum+54];
	  /*...Block 8: t0E,t1E,t2E,t3E	*/
	  #ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	  #endif
		b2p1Cr = b[rdum+14];	b2p1Ci = b[idum+14];
		b2p1Dr = b[rdum+46];	b2p1Di = b[idum+46];
		b2p1Er = b[rdum+30];	b2p1Ei = b[idum+30];
		b2p1Fr = b[rdum+62];	b2p1Fi = b[idum+62];

	  }	// endif(c_arr)

	  if(j1 == 0)	/* NB: mustn't use re, im as temps in exe1-3 section, since those contain saved sincos data for exe4 block. */
	  {				// Jul 2015: Now moot since to avoid name-clashes in FGT61-mode, have renamed the doubles re,im -> RT,IT
		/*...j = 0 (for which M(0) = Re{H(0)}*Re{I(0)} + i*Im{H(0)}*Im{I(0)}) is done separately...*/
		rt = a1p00r;
		a1p00r = (rt + a1p00i)*(b1p00r + b1p00i);
		a1p00i = (rt - a1p00i)*(b1p00r - b1p00i);
		rt = a1p00r;
		a1p00r = 0.5*(rt + a1p00i);
		a1p00i = 0.5*(rt - a1p00i);
		/*
		!...as is j = N/2 (for which M(j) = H(j)*I(j)). Note that under bit-reversal the N/2 element gets mapped into
		!   the second complex data slot, i.e. is adjacent to the starting element.
		*/
		rt = a1p01r;
		a1p01r = rt*b1p01r - a1p01i*b1p01i;
		a1p01i = rt*b1p01i + a1p01i*b1p01r;
		// exe1:
		pair_mul(&a1p02r,&a1p02i,&a1p03r,&a1p03i, b1p02r,b1p02i,b1p03r,b1p03i,0.0,1.0);		rt=ISRT2;	it=ISRT2;
		// exe2:
		pair_mul(&a1p04r,&a1p04i,&a1p07r,&a1p07i, b1p04r,b1p04i,b1p07r,b1p07i, rt, it);		t03=c; t04=0; t05=s; t06=0;
		pair_mul(&a1p06r,&a1p06i,&a1p05r,&a1p05i, b1p06r,b1p06i,b1p05r,b1p05i,-it, rt);		rt=t03-t04;		it=t05+t06;
		// exe3:
		pair_mul(&a1p08r,&a1p08i,&a1p0Fr,&a1p0Fi, b1p08r,b1p08i,b1p0Fr,b1p0Fi, rt, it);
		pair_mul(&a1p0Ar,&a1p0Ai,&a1p0Dr,&a1p0Di, b1p0Ar,b1p0Ai,b1p0Dr,b1p0Di,-it, rt);		rt=t05-t06;		it=t03+t04;
		pair_mul(&a1p0Cr,&a1p0Ci,&a1p0Br,&a1p0Bi, b1p0Cr,b1p0Ci,b1p0Br,b1p0Bi, rt, it);		t03=c32_1; t04=0; t05=s32_1; t06=0;
		pair_mul(&a1p0Er,&a1p0Ei,&a1p09r,&a1p09i, b1p0Er,b1p0Ei,b1p09r,b1p09i,-it, rt);		rt=t03-t04;		it=t05+t06;
		// exe4:
		pair_mul(&a1p10r,&a1p10i,&a1p1Fr,&a1p1Fi, b1p10r,b1p10i,b1p1Fr,b1p1Fi, rt, it);		t07=c32_3; t08=0; t09=s32_3; t0A=0;
		pair_mul(&a1p12r,&a1p12i,&a1p1Dr,&a1p1Di, b1p12r,b1p12i,b1p1Dr,b1p1Di,-it, rt);		rt=t09-t0A;		it=t07+t08;
		pair_mul(&a1p14r,&a1p14i,&a1p1Br,&a1p1Bi, b1p14r,b1p14i,b1p1Br,b1p1Bi, rt, it);
		pair_mul(&a1p16r,&a1p16i,&a1p19r,&a1p19i, b1p16r,b1p16i,b1p19r,b1p19i,-it, rt);		rt=t07-t08;		it=t09+t0A;
		pair_mul(&a1p18r,&a1p18i,&a1p17r,&a1p17i, b1p18r,b1p18i,b1p17r,b1p17i, rt, it);
		pair_mul(&a1p1Ar,&a1p1Ai,&a1p15r,&a1p15i, b1p1Ar,b1p1Ai,b1p15r,b1p15i,-it, rt);		rt=t05-t06;		it=t03+t04;
		pair_mul(&a1p1Cr,&a1p1Ci,&a1p13r,&a1p13i, b1p1Cr,b1p1Ci,b1p13r,b1p13i, rt, it);
		pair_mul(&a1p1Er,&a1p1Ei,&a1p11r,&a1p11i, b1p1Er,b1p1Ei,b1p11r,b1p11i,-it, rt);
		// exe5:
		pair_mul(&a2p00r,&a2p00i,&a2p1Fr,&a2p1Fi, b2p00r,b2p00i,b2p1Fr,b2p1Fi, re, im);
		pair_mul(&a2p02r,&a2p02i,&a2p1Dr,&a2p1Di, b2p02r,b2p02i,b2p1Dr,b2p1Di,-im, re);		rt=(re-im)*ISRT2;	it=(re+im)*ISRT2;
		pair_mul(&a2p04r,&a2p04i,&a2p1Br,&a2p1Bi, b2p04r,b2p04i,b2p1Br,b2p1Bi, rt, it);		t03=re*c; t04=im*s; t05=re*s; t06=im*c;
		pair_mul(&a2p06r,&a2p06i,&a2p19r,&a2p19i, b2p06r,b2p06i,b2p19r,b2p19i,-it, rt);		rt=t03-t04;		it=t05+t06;
		pair_mul(&a2p08r,&a2p08i,&a2p17r,&a2p17i, b2p08r,b2p08i,b2p17r,b2p17i, rt, it);
		pair_mul(&a2p0Ar,&a2p0Ai,&a2p15r,&a2p15i, b2p0Ar,b2p0Ai,b2p15r,b2p15i,-it, rt);		rt=t05-t06;		it=t03+t04;
		pair_mul(&a2p0Cr,&a2p0Ci,&a2p13r,&a2p13i, b2p0Cr,b2p0Ci,b2p13r,b2p13i, rt, it);		t03=re*c32_1; t04=im*s32_1; t05=re*s32_1; t06=im*c32_1;
		pair_mul(&a2p0Er,&a2p0Ei,&a2p11r,&a2p11i, b2p0Er,b2p0Ei,b2p11r,b2p11i,-it, rt);		rt=t03-t04;		it=t05+t06;
		pair_mul(&a2p10r,&a2p10i,&a2p0Fr,&a2p0Fi, b2p10r,b2p10i,b2p0Fr,b2p0Fi, rt, it);		t07=re*c32_3; t08=im*s32_3; t09=re*s32_3; t0A=im*c32_3;
		pair_mul(&a2p12r,&a2p12i,&a2p0Dr,&a2p0Di, b2p12r,b2p12i,b2p0Dr,b2p0Di,-it, rt);		rt=t09-t0A;		it=t07+t08;
		pair_mul(&a2p14r,&a2p14i,&a2p0Br,&a2p0Bi, b2p14r,b2p14i,b2p0Br,b2p0Bi, rt, it);
		pair_mul(&a2p16r,&a2p16i,&a2p09r,&a2p09i, b2p16r,b2p16i,b2p09r,b2p09i,-it, rt);		rt=t07-t08;		it=t09+t0A;
		pair_mul(&a2p18r,&a2p18i,&a2p07r,&a2p07i, b2p18r,b2p18i,b2p07r,b2p07i, rt, it);
		pair_mul(&a2p1Ar,&a2p1Ai,&a2p05r,&a2p05i, b2p1Ar,b2p1Ai,b2p05r,b2p05i,-it, rt);		rt=t05-t06;		it=t03+t04;
		pair_mul(&a2p1Cr,&a2p1Ci,&a2p03r,&a2p03i, b2p1Cr,b2p1Ci,b2p03r,b2p03i, rt, it);
		pair_mul(&a2p1Er,&a2p1Ei,&a2p01r,&a2p01i, b2p1Er,b2p1Ei,b2p01r,b2p01i,-it, rt);
	  }
	  else
	  {
		t01=(re-im)*ISRT2;	t02=(re+im)*ISRT2;
		re0=re*c;	im0=im*s;	re1=re*s;	im1=im*c;
		t03=re0-im0;	t04=re1+im1;
		t05=re1-im1;	t06=re0+im0;
		re0=re*c32_1;	im0=im*s32_1;	re1=re*s32_1;	im1=im*c32_1;
		t07=re0-im0;	t08=re1+im1;
		t09=re1-im1;	t0A=re0+im0;
		re0=re*c32_3;	im0=im*s32_3;	re1=re*s32_3;	im1=im*c32_3;
		t0B=re1-im1;	t0C=re0+im0;
		t0D=re0-im0;	t0E=re1+im1;
		/*
		if(j1 == 256) {
			printf("j1,j2 = [%u],[%u]: Inputs to PAIR_MUL_4:\n",j1pad,j2pad);
			printf("re,im,t1-6: %18.10e,%18.10e, %18.10e,%18.10e, %18.10e,%18.10e, %18.10e,%18.10e\n",re ,im ,t00,t01,t02,t03,t04,t05,t06);
			printf("      t7-E: %18.10e,%18.10e, %18.10e,%18.10e, %18.10e,%18.10e, %18.10e,%18.10e\n",t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F);
		// BR pattern: [0,16,8,24] + [0,4,2,6,1,5,3,7]
			printf("0 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p00r,a2p00r,a1p00i,a2p00i,b1p00r,b2p00r,b1p00i,b2p00i);
			printf("1 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p10r,a2p10r,a1p10i,a2p10i,b1p10r,b2p10r,b1p10i,b2p10i);
			printf("2 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p08r,a2p08r,a1p08i,a2p08i,b1p08r,b2p08r,b1p08i,b2p08i);
			printf("3 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p18r,a2p18r,a1p18i,a2p18i,b1p18r,b2p18r,b1p18i,b2p18i);
			printf("4 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p04r,a2p04r,a1p04i,a2p04i,b1p04r,b2p04r,b1p04i,b2p04i);
			printf("5 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p14r,a2p14r,a1p14i,a2p14i,b1p14r,b2p14r,b1p14i,b2p14i);
			printf("6 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p0Cr,a2p0Cr,a1p0Ci,a2p0Ci,b1p0Cr,b2p0Cr,b1p0Ci,b2p0Ci);
			printf("7 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p1Cr,a2p1Cr,a1p1Ci,a2p1Ci,b1p1Cr,b2p1Cr,b1p1Ci,b2p1Ci);
			printf("8 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p02r,a2p02r,a1p02i,a2p02i,b1p02r,b2p02r,b1p02i,b2p02i);
			printf("9 : %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p12r,a2p12r,a1p12i,a2p12i,b1p12r,b2p12r,b1p12i,b2p12i);
			printf("10: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p0Ar,a2p0Ar,a1p0Ai,a2p0Ai,b1p0Ar,b2p0Ar,b1p0Ai,b2p0Ai);
			printf("11: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p1Ar,a2p1Ar,a1p1Ai,a2p1Ai,b1p1Ar,b2p1Ar,b1p1Ai,b2p1Ai);
			printf("12: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p06r,a2p06r,a1p06i,a2p06i,b1p06r,b2p06r,b1p06i,b2p06i);
			printf("13: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p16r,a2p16r,a1p16i,a2p16i,b1p16r,b2p16r,b1p16i,b2p16i);
			printf("14: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p0Er,a2p0Er,a1p0Ei,a2p0Ei,b1p0Er,b2p0Er,b1p0Ei,b2p0Ei);
			printf("15: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p1Er,a2p1Er,a1p1Ei,a2p1Ei,b1p1Er,b2p1Er,b1p1Ei,b2p1Ei);
			printf("16: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p01r,a2p01r,a1p01i,a2p01i,b1p01r,b2p01r,b1p01i,b2p01i);
			printf("17: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p11r,a2p11r,a1p11i,a2p11i,b1p11r,b2p11r,b1p11i,b2p11i);
			printf("18: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p09r,a2p09r,a1p09i,a2p09i,b1p09r,b2p09r,b1p09i,b2p09i);
			printf("19: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p19r,a2p19r,a1p19i,a2p19i,b1p19r,b2p19r,b1p19i,b2p19i);
			printf("20: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p05r,a2p05r,a1p05i,a2p05i,b1p05r,b2p05r,b1p05i,b2p05i);
			printf("21: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p15r,a2p15r,a1p15i,a2p15i,b1p15r,b2p15r,b1p15i,b2p15i);
			printf("22: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p0Dr,a2p0Dr,a1p0Di,a2p0Di,b1p0Dr,b2p0Dr,b1p0Di,b2p0Di);
			printf("23: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p1Dr,a2p1Dr,a1p1Di,a2p1Di,b1p1Dr,b2p1Dr,b1p1Di,b2p1Di);
			printf("24: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p03r,a2p03r,a1p03i,a2p03i,b1p03r,b2p03r,b1p03i,b2p03i);
			printf("25: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p13r,a2p13r,a1p13i,a2p13i,b1p13r,b2p13r,b1p13i,b2p13i);
			printf("26: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p0Br,a2p0Br,a1p0Bi,a2p0Bi,b1p0Br,b2p0Br,b1p0Bi,b2p0Bi);
			printf("27: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p1Br,a2p1Br,a1p1Bi,a2p1Bi,b1p1Br,b2p1Br,b1p1Bi,b2p1Bi);
			printf("28: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p07r,a2p07r,a1p07i,a2p07i,b1p07r,b2p07r,b1p07i,b2p07i);
			printf("29: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p17r,a2p17r,a1p17i,a2p17i,b1p17r,b2p17r,b1p17i,b2p17i);
			printf("30: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p0Fr,a2p0Fr,a1p0Fi,a2p0Fi,b1p0Fr,b2p0Fr,b1p0Fi,b2p0Fi);
			printf("31: %18.10e,%18.10e,%18.10e,%18.10e | %18.10e,%18.10e,%18.10e,%18.10e\n",a1p1Fr,a2p1Fr,a1p1Fi,a2p1Fi,b1p1Fr,b2p1Fr,b1p1Fi,b2p1Fi);
		//	exit(0);
		}
		*/
		PAIR_MUL_4( a1p00r,a1p00i,a2p1Fr,a2p1Fi, b1p00r,b1p00i,b2p1Fr,b2p1Fi,
					a2p00r,a2p00i,a1p1Fr,a1p1Fi, b2p00r,b2p00i,b1p1Fr,b1p1Fi,
					a1p02r,a1p02i,a2p1Dr,a2p1Di, b1p02r,b1p02i,b2p1Dr,b2p1Di,
					a2p02r,a2p02i,a1p1Dr,a1p1Di, b2p02r,b2p02i,b1p1Dr,b1p1Di,  re, im, t09,t0A);
		PAIR_MUL_4( a1p04r,a1p04i,a2p1Br,a2p1Bi, b1p04r,b1p04i,b2p1Br,b2p1Bi,
					a2p04r,a2p04i,a1p1Br,a1p1Bi, b2p04r,b2p04i,b1p1Br,b1p1Bi,
					a1p06r,a1p06i,a2p19r,a2p19i, b1p06r,b1p06i,b2p19r,b2p19i,
					a2p06r,a2p06i,a1p19r,a1p19i, b2p06r,b2p06i,b1p19r,b1p19i, t01,t02, t0D,t0E);
		PAIR_MUL_4( a1p08r,a1p08i,a2p17r,a2p17i, b1p08r,b1p08i,b2p17r,b2p17i,
					a2p08r,a2p08i,a1p17r,a1p17i, b2p08r,b2p08i,b1p17r,b1p17i,
					a1p0Ar,a1p0Ai,a2p15r,a2p15i, b1p0Ar,b1p0Ai,b2p15r,b2p15i,
					a2p0Ar,a2p0Ai,a1p15r,a1p15i, b2p0Ar,b2p0Ai,b1p15r,b1p15i, t03,t04, t0B,t0C);
		PAIR_MUL_4( a1p0Cr,a1p0Ci,a2p13r,a2p13i, b1p0Cr,b1p0Ci,b2p13r,b2p13i,
					a2p0Cr,a2p0Ci,a1p13r,a1p13i, b2p0Cr,b2p0Ci,b1p13r,b1p13i,
					a1p0Er,a1p0Ei,a2p11r,a2p11i, b1p0Er,b1p0Ei,b2p11r,b2p11i,
					a2p0Er,a2p0Ei,a1p11r,a1p11i, b2p0Er,b2p0Ei,b1p11r,b1p11i, t05,t06, t07,t08);
		PAIR_MUL_4( a1p10r,a1p10i,a2p0Fr,a2p0Fi, b1p10r,b1p10i,b2p0Fr,b2p0Fi,
					a2p10r,a2p10i,a1p0Fr,a1p0Fi, b2p10r,b2p10i,b1p0Fr,b1p0Fi,
					a1p12r,a1p12i,a2p0Dr,a2p0Di, b1p12r,b1p12i,b2p0Dr,b2p0Di,
					a2p12r,a2p12i,a1p0Dr,a1p0Di, b2p12r,b2p12i,b1p0Dr,b1p0Di, t07,t08, t05,t06);
		PAIR_MUL_4( a1p14r,a1p14i,a2p0Br,a2p0Bi, b1p14r,b1p14i,b2p0Br,b2p0Bi,
					a2p14r,a2p14i,a1p0Br,a1p0Bi, b2p14r,b2p14i,b1p0Br,b1p0Bi,
					a1p16r,a1p16i,a2p09r,a2p09i, b1p16r,b1p16i,b2p09r,b2p09i,
					a2p16r,a2p16i,a1p09r,a1p09i, b2p16r,b2p16i,b1p09r,b1p09i, t0B,t0C, t03,t04);
		PAIR_MUL_4( a1p18r,a1p18i,a2p07r,a2p07i, b1p18r,b1p18i,b2p07r,b2p07i,
					a2p18r,a2p18i,a1p07r,a1p07i, b2p18r,b2p18i,b1p07r,b1p07i,
					a1p1Ar,a1p1Ai,a2p05r,a2p05i, b1p1Ar,b1p1Ai,b2p05r,b2p05i,
					a2p1Ar,a2p1Ai,a1p05r,a1p05i, b2p1Ar,b2p1Ai,b1p05r,b1p05i, t0D,t0E, t01,t02);
		PAIR_MUL_4( a1p1Cr,a1p1Ci,a2p03r,a2p03i, b1p1Cr,b1p1Ci,b2p03r,b2p03i,
					a2p1Cr,a2p1Ci,a1p03r,a1p03i, b2p1Cr,b2p1Ci,b1p03r,b1p03i,
					a1p1Er,a1p1Ei,a2p01r,a2p01i, b1p1Er,b1p1Ei,b2p01r,b2p01i,
					a2p1Er,a2p1Ei,a1p01r,a1p01i, b2p1Er,b2p1Ei,b1p01r,b1p01i, t09,t0A,  re, im);
	  }

	} else {	// fwd_fft_only = 0: Normal execution, dyadic-squaring followed by iFFT:

	/*...send the pairs of complex elements which are to be combined and sincos temporaries needed for the squaring to a
     small subroutine. The j1 = 0 case is again exceptional.	*/

	  if(j1==0)    	/* NB: mustn't use re, im as temps in exe1-3 section, since those contain saved sincos data for exe4 block.	*/
	  {
	/*...j=0 (for which I(0)=Re{H(0)}^2+i*Im{H(0)}^2) is done separately...	*/
		rt=a1p00r;
		a1p00r=(rt+a1p00i)*(rt+a1p00i);
		a1p00i=(rt-a1p00i)*(rt-a1p00i);
		rt=a1p00r;
		a1p00r=0.5*(rt+a1p00i);
		a1p00i=0.5*(rt-a1p00i);
	/*...as is j=N/2 (for which I(j)=H(j)^2). Note that under bit-reversal the N/2 element gets mapped into
     the second complex data slot, i.e. is adjacent to the starting element.	*/
		rt =a1p01r*a1p01i;
		a1p01r=(a1p01r+a1p01i)*(a1p01r-a1p01i);
		a1p01i=rt+rt;
		pair_square(&a1p02r,&a1p02i,&a1p03r,&a1p03i,0.0,1.0);		/* exe1	*/
		rt=ISRT2;	it=ISRT2;
		pair_square(&a1p04r,&a1p04i,&a1p07r,&a1p07i, rt, it);		/* exe2	*/
		pair_square(&a1p06r,&a1p06i,&a1p05r,&a1p05i,-it, rt);		/* exe2	*/
		t03=c; t04=0; t05=s; t06=0;
		rt=t03-t04;		it=t05+t06;
		pair_square(&a1p08r,&a1p08i,&a1p0Fr,&a1p0Fi, rt, it);		/* exe3	*/
		pair_square(&a1p0Ar,&a1p0Ai,&a1p0Dr,&a1p0Di,-it, rt);		/* exe3	*/
		rt=t05-t06;		it=t03+t04;
		pair_square(&a1p0Cr,&a1p0Ci,&a1p0Br,&a1p0Bi, rt, it);		/* exe3	*/
		pair_square(&a1p0Er,&a1p0Ei,&a1p09r,&a1p09i,-it, rt);		/* exe3	*/
		t03=c32_1; t04=0; t05=s32_1; t06=0;
		rt=t03-t04;		it=t05+t06;
		pair_square(&a1p10r,&a1p10i,&a1p1Fr,&a1p1Fi, rt, it);		/* exe4	*/
		pair_square(&a1p12r,&a1p12i,&a1p1Dr,&a1p1Di,-it, rt);		/* exe4	*/
		t07=c32_3; t08=0; t09=s32_3; t0A=0;
		rt=t09-t0A;		it=t07+t08;
		pair_square(&a1p14r,&a1p14i,&a1p1Br,&a1p1Bi, rt, it);		/* exe4	*/
		pair_square(&a1p16r,&a1p16i,&a1p19r,&a1p19i,-it, rt);		/* exe4	*/
		rt=t07-t08;		it=t09+t0A;
		pair_square(&a1p18r,&a1p18i,&a1p17r,&a1p17i, rt, it);		/* exe4	*/
		pair_square(&a1p1Ar,&a1p1Ai,&a1p15r,&a1p15i,-it, rt);		/* exe4	*/
		rt=t05-t06;		it=t03+t04;
		pair_square(&a1p1Cr,&a1p1Ci,&a1p13r,&a1p13i, rt, it);		/* exe4	*/
		pair_square(&a1p1Er,&a1p1Ei,&a1p11r,&a1p11i,-it, rt);		/* exe4	*/
	/*   exe5:	*/
		pair_square(&a2p00r,&a2p00i,&a2p1Fr,&a2p1Fi, re, im);
		pair_square(&a2p02r,&a2p02i,&a2p1Dr,&a2p1Di,-im, re);
		rt=(re-im)*ISRT2;	it=(re+im)*ISRT2;
		pair_square(&a2p04r,&a2p04i,&a2p1Br,&a2p1Bi, rt, it);
		pair_square(&a2p06r,&a2p06i,&a2p19r,&a2p19i,-it, rt);
		t03=re*c; t04=im*s; t05=re*s; t06=im*c;
		rt=t03-t04;		it=t05+t06;
		pair_square(&a2p08r,&a2p08i,&a2p17r,&a2p17i, rt, it);
		pair_square(&a2p0Ar,&a2p0Ai,&a2p15r,&a2p15i,-it, rt);
		rt=t05-t06;		it=t03+t04;
		pair_square(&a2p0Cr,&a2p0Ci,&a2p13r,&a2p13i, rt, it);
		pair_square(&a2p0Er,&a2p0Ei,&a2p11r,&a2p11i,-it, rt);
		t03=re*c32_1; t04=im*s32_1; t05=re*s32_1; t06=im*c32_1;
		rt=t03-t04;		it=t05+t06;
		pair_square(&a2p10r,&a2p10i,&a2p0Fr,&a2p0Fi, rt, it);
		pair_square(&a2p12r,&a2p12i,&a2p0Dr,&a2p0Di,-it, rt);
		t07=re*c32_3; t08=im*s32_3; t09=re*s32_3; t0A=im*c32_3;
		rt=t09-t0A;		it=t07+t08;
		pair_square(&a2p14r,&a2p14i,&a2p0Br,&a2p0Bi, rt, it);
		pair_square(&a2p16r,&a2p16i,&a2p09r,&a2p09i,-it, rt);
		rt=t07-t08;		it=t09+t0A;
		pair_square(&a2p18r,&a2p18i,&a2p07r,&a2p07i, rt, it);
		pair_square(&a2p1Ar,&a2p1Ai,&a2p05r,&a2p05i,-it, rt);
		rt=t05-t06;		it=t03+t04;
		pair_square(&a2p1Cr,&a2p1Ci,&a2p03r,&a2p03i, rt, it);
		pair_square(&a2p1Er,&a2p1Ei,&a2p01r,&a2p01i,-it, rt);
	  }
	  else
	  {
		t01=(re-im)*ISRT2;	t02=(re+im)*ISRT2;
		re0=re*c;	im0=im*s;	re1=re*s;	im1=im*c;
		t03=re0-im0;	t04=re1+im1;
		t05=re1-im1;	t06=re0+im0;
		re0=re*c32_1;	im0=im*s32_1;	re1=re*s32_1;	im1=im*c32_1;
		t07=re0-im0;	t08=re1+im1;
		t09=re1-im1;	t0A=re0+im0;
		re0=re*c32_3;	im0=im*s32_3;	re1=re*s32_3;	im1=im*c32_3;
		t0B=re1-im1;	t0C=re0+im0;
		t0D=re0-im0;	t0E=re1+im1;
		PAIR_SQUARE_4(a1p00r,a1p00i,a2p1Fr,a2p1Fi,a1p02r,a1p02i,a2p1Dr,a2p1Di, re, im
					, a2p02r,a2p02i,a1p1Dr,a1p1Di,a2p00r,a2p00i,a1p1Fr,a1p1Fi,t09,t0A);
		PAIR_SQUARE_4(a1p04r,a1p04i,a2p1Br,a2p1Bi,a1p06r,a1p06i,a2p19r,a2p19i,t01,t02
					, a2p06r,a2p06i,a1p19r,a1p19i,a2p04r,a2p04i,a1p1Br,a1p1Bi,t0D,t0E);
		PAIR_SQUARE_4(a1p08r,a1p08i,a2p17r,a2p17i,a1p0Ar,a1p0Ai,a2p15r,a2p15i,t03,t04
					, a2p0Ar,a2p0Ai,a1p15r,a1p15i,a2p08r,a2p08i,a1p17r,a1p17i,t0B,t0C);
		PAIR_SQUARE_4(a1p0Cr,a1p0Ci,a2p13r,a2p13i,a1p0Er,a1p0Ei,a2p11r,a2p11i,t05,t06
					, a2p0Er,a2p0Ei,a1p11r,a1p11i,a2p0Cr,a2p0Ci,a1p13r,a1p13i,t07,t08);
		PAIR_SQUARE_4(a1p10r,a1p10i,a2p0Fr,a2p0Fi,a1p12r,a1p12i,a2p0Dr,a2p0Di,t07,t08
					, a2p12r,a2p12i,a1p0Dr,a1p0Di,a2p10r,a2p10i,a1p0Fr,a1p0Fi,t05,t06);
		PAIR_SQUARE_4(a1p14r,a1p14i,a2p0Br,a2p0Bi,a1p16r,a1p16i,a2p09r,a2p09i,t0B,t0C
					, a2p16r,a2p16i,a1p09r,a1p09i,a2p14r,a2p14i,a1p0Br,a1p0Bi,t03,t04);
		PAIR_SQUARE_4(a1p18r,a1p18i,a2p07r,a2p07i,a1p1Ar,a1p1Ai,a2p05r,a2p05i,t0D,t0E
					, a2p1Ar,a2p1Ai,a1p05r,a1p05i,a2p18r,a2p18i,a1p07r,a1p07i,t01,t02);
		PAIR_SQUARE_4(a1p1Cr,a1p1Ci,a2p03r,a2p03i,a1p1Er,a1p1Ei,a2p01r,a2p01i,t09,t0A
					, a2p1Er,a2p1Ei,a1p01r,a1p01i,a2p1Cr,a2p1Ci,a1p03r,a1p03i, re, im);
	  }	/* endif(j1==0)	*/

	}	// endif(fwd_fft_only == 1)

	/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/

	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;
	#if PFETCH
	addr = &a[rdum+64];
	#endif
		/*   gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 IDIT transforms... */

		/*...Block 1: */
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		t02=a1p00r-a1p01r;	t03=a1p00i-a1p01i;
		t00=a1p00r+a1p01r;	t01=a1p00i+a1p01i;
		t06=a1p02r-a1p03r;	t07=a1p02i-a1p03i;
		t04=a1p02r+a1p03r;	t05=a1p02i+a1p03i;
		rt =t04;		it =t05;
		t04=t00-rt;		t05=t01-it;
		t00=t00+rt;		t01=t01+it;
		rt =t06;		it =t07;
		t06=t02-it;		t07=t03+rt;
		t02=t02+it;		t03=t03-rt;
		t0A=a1p04r-a1p05r;	t0B=a1p04i-a1p05i;
		t08=a1p04r+a1p05r;	t09=a1p04i+a1p05i;
		t0E=a1p06r-a1p07r;	t0F=a1p06i-a1p07i;
		t0C=a1p06r+a1p07r;	t0D=a1p06i+a1p07i;
		rt =t0C;		it =t0D;
		t0C=t08-rt;		t0D=t09-it;
		t08=t08+rt;		t09=t09+it;
		rt =t0E;		it =t0F;
		t0E=t0A-it;		t0F=t0B+rt;
		t0A=t0A+it;		t0B=t0B-rt;
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t08;		it =t09;
		t08=t00-rt;		t09=t01-it;
		t00=t00+rt;		t01=t01+it;
		rt =t0C;		it =t0D;
		t0C=t04-it;		t0D=t05+rt;
		t04=t04+it;		t05=t05-rt;
		rt =(t0A+t0B)*ISRT2;it =(t0A-t0B)*ISRT2;
		t0A=t02-rt;		t0B=t03+it;
		t02=t02+rt;		t03=t03-it;
		rt =(t0E-t0F)*ISRT2;it =(t0F+t0E)*ISRT2;
		t0E=t06+rt;		t0F=t07+it;
		t06=t06-rt;		t07=t07-it;
	/*...Block 2:	*/
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		t12=a1p08r-a1p09r;	t13=a1p08i-a1p09i;
		t10=a1p08r+a1p09r;	t11=a1p08i+a1p09i;
		t16=a1p0Ar-a1p0Br;	t17=a1p0Ai-a1p0Bi;
		t14=a1p0Ar+a1p0Br;	t15=a1p0Ai+a1p0Bi;
		rt =t14;		it =t15;
		t14=t10-rt;		t15=t11-it;
		t10=t10+rt;		t11=t11+it;
		rt =t16;		it =t17;
		t16=t12-it;		t17=t13+rt;
		t12=t12+it;		t13=t13-rt;
		t1A=a1p0Cr-a1p0Dr;	t1B=a1p0Ci-a1p0Di;
		t18=a1p0Cr+a1p0Dr;	t19=a1p0Ci+a1p0Di;
		t1E=a1p0Er-a1p0Fr;	t1F=a1p0Ei-a1p0Fi;
		t1C=a1p0Er+a1p0Fr;	t1D=a1p0Ei+a1p0Fi;
		rt =t1C;		it =t1D;
		t1C=t18-rt;		t1D=t19-it;
		t18=t18+rt;		t19=t19+it;
		rt =t1E;		it =t1F;
		t1E=t1A-it;		t1F=t1B+rt;
		t1A=t1A+it;		t1B=t1B-rt;
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t18;		it =t19;
		t18=t10-rt;		t19=t11-it;
		t10=t10+rt;		t11=t11+it;
		rt =t1C;		it =t1D;
		t1C=t14-it;		t1D=t15+rt;
		t14=t14+it;		t15=t15-rt;
		rt =(t1A+t1B)*ISRT2;it =(t1A-t1B)*ISRT2;
		t1A=t12-rt;		t1B=t13+it;
		t12=t12+rt;		t13=t13-it;
		rt =(t1E-t1F)*ISRT2;it =(t1F+t1E)*ISRT2;
		t1E=t16+rt;		t1F=t17+it;
		t16=t16-rt;		t17=t17-it;
	/*...Block 3:	*/
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		t22=a1p10r-a1p11r;	t23=a1p10i-a1p11i;
		t20=a1p10r+a1p11r;	t21=a1p10i+a1p11i;
		t26=a1p12r-a1p13r;	t27=a1p12i-a1p13i;
		t24=a1p12r+a1p13r;	t25=a1p12i+a1p13i;
		rt =t24;		it =t25;
		t24=t20-rt;		t25=t21-it;
		t20=t20+rt;		t21=t21+it;
		rt =t26;		it =t27;
		t26=t22-it;		t27=t23+rt;
		t22=t22+it;		t23=t23-rt;
		t2A=a1p14r-a1p15r;	t2B=a1p14i-a1p15i;
		t28=a1p14r+a1p15r;	t29=a1p14i+a1p15i;
		t2E=a1p16r-a1p17r;	t2F=a1p16i-a1p17i;
		t2C=a1p16r+a1p17r;	t2D=a1p16i+a1p17i;
		rt =t2C;		it =t2D;
		t2C=t28-rt;		t2D=t29-it;
		t28=t28+rt;		t29=t29+it;
		rt =t2E;		it =t2F;
		t2E=t2A-it;		t2F=t2B+rt;
		t2A=t2A+it;		t2B=t2B-rt;
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t28;		it =t29;
		t28=t20-rt;		t29=t21-it;
		t20=t20+rt;		t21=t21+it;
		rt =t2C;		it =t2D;
		t2C=t24-it;		t2D=t25+rt;
		t24=t24+it;		t25=t25-rt;
		rt =(t2A+t2B)*ISRT2;it =(t2A-t2B)*ISRT2;
		t2A=t22-rt;		t2B=t23+it;
		t22=t22+rt;		t23=t23-it;
		rt =(t2E-t2F)*ISRT2;it =(t2F+t2E)*ISRT2;
		t2E=t26+rt;		t2F=t27+it;
		t26=t26-rt;		t27=t27-it;
	/*...Block 4:	*/
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		t32=a1p18r-a1p19r;	t33=a1p18i-a1p19i;
		t30=a1p18r+a1p19r;	t31=a1p18i+a1p19i;
		t36=a1p1Ar-a1p1Br;	t37=a1p1Ai-a1p1Bi;
		t34=a1p1Ar+a1p1Br;	t35=a1p1Ai+a1p1Bi;
		rt =t34;		it =t35;
		t34=t30-rt;		t35=t31-it;
		t30=t30+rt;		t31=t31+it;
		rt =t36;		it =t37;
		t36=t32-it;		t37=t33+rt;
		t32=t32+it;		t33=t33-rt;
		t3A=a1p1Cr-a1p1Dr;	t3B=a1p1Ci-a1p1Di;
		t38=a1p1Cr+a1p1Dr;	t39=a1p1Ci+a1p1Di;
		t3E=a1p1Er-a1p1Fr;	t3F=a1p1Ei-a1p1Fi;
		t3C=a1p1Er+a1p1Fr;	t3D=a1p1Ei+a1p1Fi;
		rt =t3C;		it =t3D;
		t3C=t38-rt;		t3D=t39-it;
		t38=t38+rt;		t39=t39+it;
		rt =t3E;		it =t3F;
		t3E=t3A-it;		t3F=t3B+rt;
		t3A=t3A+it;		t3B=t3B-rt;
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t38;		it =t39;
		t38=t30-rt;		t39=t31-it;
		t30=t30+rt;		t31=t31+it;
		rt =t3C;		it =t3D;
		t3C=t34-it;		t3D=t35+rt;
		t34=t34+it;		t35=t35-rt;
		rt =(t3A+t3B)*ISRT2;it =(t3A-t3B)*ISRT2;
		t3A=t32-rt;		t3B=t33+it;
		t32=t32+rt;		t33=t33-it;
		rt =(t3E-t3F)*ISRT2;it =(t3F+t3E)*ISRT2;
		t3E=t36+rt;		t3F=t37+it;
		t36=t36-rt;		t37=t37-it;
	//...and now do eight radix-4 transforms, including the internal twiddle factors:
	/*...Block 1: t00,t10,t20,t30	*/
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		rt =t10;	t10=t00-rt;	t00=t00+rt;
		it =t11;	t11=t01-it;	t01=t01+it;
		rt =t30;	t30=t20-rt;	t20=t20+rt;
		it =t31;	t31=t21-it;	t21=t21+it;
		a[rdum   ]=t00+t20;			a[idum   ]=t01+t21;
		t00       =t00-t20;			t01       =t01-t21;
		a[rdum+32]=t00*cA10+t01*sA10;	a[idum+32]=t01*cA10-t00*sA10;
		rt        =t10+t31;			it        =t11-t30;
		t10       =t10-t31;			t11       =t11+t30;
		a[rdum+16]=rt *cA08+it *sA08;	a[idum+16]=it *cA08-rt *sA08;
		a[rdum+48]=t10*cA18+t11*sA18;	a[idum+48]=t11*cA18-t10*sA18;
	/*...Block 5: t08,t18,t28,t38	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t18;	t18=t08-t19;	t08=t08+t19;
					t19=t09+rt;		t09=t09-rt;
		rt =(t29+t28)*ISRT2;		t29=(t29-t28)*ISRT2;		t28=rt;
		rt =(t38-t39)*ISRT2;		it =(t38+t39)*ISRT2;
		t38=t28+rt;					t28=t28-rt;
		t39=t29+it;					t29=t29-it;
		rt        =t08+t28;			it        =t09+t29;
		t08       =t08-t28;			t09       =t09-t29;
		a[rdum+8 ]=rt *cA04+it *sA04;	a[idum+8 ]=it *cA04-rt *sA04;
		a[rdum+40]=t08*cA14+t09*sA14;	a[idum+40]=t09*cA14-t08*sA14;
		rt        =t18+t39;			it        =t19-t38;
		t18       =t18-t39;			t19       =t19+t38;
		a[rdum+24]=rt *cA0C+it *sA0C;	a[idum+24]=it *cA0C-rt *sA0C;
		a[rdum+56]=t18*cA1C+t19*sA1C;	a[idum+56]=t19*cA1C-t18*sA1C;
	/*...Block 3: t04,t14,t24,t34	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		rt =(t15+t14)*ISRT2;		it =(t15-t14)*ISRT2;
		t14=t04-rt;					t04=t04+rt;
		t15=t05-it;					t05=t05+it;
		rt =t24*c + t25*s;			t25=t25*c - t24*s;		t24=rt;
		rt =t34*s + t35*c;			it =t35*s - t34*c;
		t34=t24-rt;					t24=t24+rt;
		t35=t25-it;					t25=t25+it;
		rt        =t04+t24;			it        =t05+t25;
		t04       =t04-t24;			t05       =t05-t25;
		a[rdum+4 ]=rt *cA02+it *sA02;	a[idum+4 ]=it *cA02-rt *sA02;
		a[rdum+36]=t04*cA12+t05*sA12;	a[idum+36]=t05*cA12-t04*sA12;
		rt        =t14+t35;			it        =t15-t34;
		t14       =t14-t35;			t15       =t15+t34;
		a[rdum+20]=rt *cA0A+it *sA0A;	a[idum+20]=it *cA0A-rt *sA0A;
		a[rdum+52]=t14*cA1A+t15*sA1A;	a[idum+52]=t15*cA1A-t14*sA1A;
	/*...Block 7: t0C,t1C,t2C,t3C	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =(t1C-t1D)*ISRT2;		it =(t1C+t1D)*ISRT2;
		t1C=t0C+rt;					t0C=t0C-rt;
		t1D=t0D+it;					t0D=t0D-it;
		rt =t2C*s + t2D*c;			t2D=t2D*s - t2C*c;		t2C=rt;
		rt =t3C*c + t3D*s;			it =t3D*c - t3C*s;
		t3C=t2C+rt;					t2C=t2C-rt;
		t3D=t2D+it;					t2D=t2D-it;
		rt        =t0C+t2C;			it        =t0D+t2D;
		t0C       =t0C-t2C;			t0D       =t0D-t2D;
		a[rdum+12]=rt *cA06+it *sA06;	a[idum+12]=it *cA06-rt *sA06;
		a[rdum+44]=t0C*cA16+t0D*sA16;	a[idum+44]=t0D*cA16-t0C*sA16;
		rt        =t1C+t3D;			it        =t1D-t3C;
		t1C       =t1C-t3D;			t1D       =t1D+t3C;
		a[rdum+28]=rt *cA0E+it *sA0E;	a[idum+28]=it *cA0E-rt *sA0E;
		a[rdum+60]=t1C*cA1E+t1D*sA1E;	a[idum+60]=t1D*cA1E-t1C*sA1E;
	/*...Block 2: t02,t12,t22,t32	*/
	#ifdef USE_AVX512
		rdum += 5;	idum += 5;
	#elif defined(USE_AVX)
		++rdum;	++idum;
	#elif defined(USE_SSE2)
		--rdum;	--idum;
	#endif
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		rt =t12*c + t13*s;			it =t13*c - t12*s;
		t12=t02-rt;					t02=t02+rt;
		t13=t03-it;					t03=t03+it;
		rt =t22*c32_1 + t23*s32_1;	t23=t23*c32_1 - t22*s32_1;	t22=rt;
		rt =t32*c32_3 + t33*s32_3;	it =t33*c32_3 - t32*s32_3;
		t32=t22-rt;					t22=t22+rt;
		t33=t23-it;					t23=t23+it;
		rt        =t02+t22;			it        =t03+t23;
		t02       =t02-t22;			t03       =t03-t23;
		a[rdum+2 ]=rt *cA01+it *sA01;	a[idum+2 ]=it *cA01-rt *sA01;
		a[rdum+34]=t02*cA11+t03*sA11;	a[idum+34]=t03*cA11-t02*sA11;
		rt        =t12+t33;			it        =t13-t32;
		t12       =t12-t33;			t13       =t13+t32;
		a[rdum+18]=rt *cA09+it *sA09;	a[idum+18]=it *cA09-rt *sA09;
		a[rdum+50]=t12*cA19+t13*sA19;	a[idum+50]=t13*cA19-t12*sA19;
	/*...Block 6: t0A,t1A,t2A,t3A	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t1A*s - t1B*c;			it =t1B*s + t1A*c;
		t1A=t0A+rt;					t0A =t0A-rt;
		t1B=t0B+it;					t0B =t0B-it;
		rt =t2A*s32_3 + t2B*c32_3;	t2B=t2B*s32_3 - t2A*c32_3;	t2A=rt;
		rt =t3A*c32_1 - t3B*s32_1;	it =t3B*c32_1 + t3A*s32_1;
		t3A=t2A+rt;					t2A=t2A-rt;
		t3B=t2B+it;					t2B=t2B-it;
		rt        =t0A+t2A;			it        =t0B+t2B;
		t0A       =t0A-t2A;			t0B       =t0B-t2B;
		a[rdum+10]=rt *cA05+it *sA05;	a[idum+10]=it *cA05-rt *sA05;
		a[rdum+42]=t0A*cA15+t0B*sA15;	a[idum+42]=t0B*cA15-t0A*sA15;
		rt        =t1A+t3B;			it        =t1B-t3A;
		t1A       =t1A-t3B;			t1B       =t1B+t3A;
		a[rdum+26]=rt *cA0D+it *sA0D;	a[idum+26]=it *cA0D-rt *sA0D;
		a[rdum+58]=t1A*cA1D+t1B*sA1D;	a[idum+58]=t1B*cA1D-t1A*sA1D;
	/*...Block 4: t06,t16,t26,t36	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		rt =t16*s + t17*c;			it =t17*s - t16*c;
		t16=t06-rt;					t06 =t06+rt;
		t17=t07-it;					t07 =t07+it;
		rt =t26*c32_3 + t27*s32_3;	t27=t27*c32_3 - t26*s32_3;	t26=rt;
		rt =t36*s32_1 - t37*c32_1;	it =t37*s32_1 + t36*c32_1;
		t36=t26+rt;					t26=t26-rt;
		t37=t27+it;					t27=t27-it;
		rt        =t06+t26;			it        =t07+t27;
		t06       =t06-t26;			t07       =t07-t27;
		a[rdum+6 ]=rt *cA03+it *sA03;	a[idum+6 ]=it *cA03-rt *sA03;
		a[rdum+38]=t06*cA13+t07*sA13;	a[idum+38]=t07*cA13-t06*sA13;
		rt        =t16+t37;			it        =t17-t36;
		t16       =t16-t37;			t17       =t17+t36;
		a[rdum+22]=rt *cA0B+it *sA0B;	a[idum+22]=it *cA0B-rt *sA0B;
		a[rdum+54]=t16*cA1B+t17*sA1B;	a[idum+54]=t17*cA1B-t16*sA1B;
	/*...Block 8: t0E,t1E,t2E,t3E	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t1E*c - t1F*s;			it =t1F*c + t1E*s;
		t1E=t0E+rt;					t0E =t0E-rt;
		t1F=t0F+it;					t0F =t0F-it;
		rt =t2E*s32_1 + t2F*c32_1;	t2F=t2F*s32_1 - t2E*c32_1;	t2E=rt;
		rt =t3E*s32_3 + t3F*c32_3;	it =t3F*s32_3 - t3E*c32_3;
		t3E=t2E+rt;					t2E=t2E-rt;
		t3F=t2F+it;					t2F=t2F-it;
		rt        =t0E+t2E;			it        =t0F+t2F;
		t0E       =t0E-t2E;			t0F       =t0F-t2F;
		a[rdum+14]=rt *cA07+it *sA07;	a[idum+14]=it *cA07-rt *sA07;
		a[rdum+46]=t0E*cA17+t0F*sA17;	a[idum+46]=t0F*cA17-t0E*sA17;
		rt        =t1E+t3F;			it        =t1F-t3E;
		t1E       =t1E-t3F;			t1F       =t1F+t3E;
		a[rdum+30]=rt *cA0F+it *sA0F;	a[idum+30]=it *cA0F-rt *sA0F;
		a[rdum+62]=t1E*cA1F+t1F*sA1F;	a[idum+62]=t1F*cA1F-t1E*sA1F;
	/*************************************************************/
	/*                  2nd set of inputs:                       */
	/*************************************************************/
		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;
	#if PFETCH
	addr = &a[rdum-64];
	#endif
		/*...Block 1: */
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
	/*...Block 1:	*/
		t02=a2p00r-a2p01r;	t03=a2p00i-a2p01i;
		t00=a2p00r+a2p01r;	t01=a2p00i+a2p01i;
		t06=a2p02r-a2p03r;	t07=a2p02i-a2p03i;
		t04=a2p02r+a2p03r;	t05=a2p02i+a2p03i;
		rt =t04;		it =t05;
		t04=t00-rt;		t05=t01-it;
		t00=t00+rt;		t01=t01+it;
		rt =t06;		it =t07;
		t06=t02-it;		t07=t03+rt;
		t02=t02+it;		t03=t03-rt;
		t0A=a2p04r-a2p05r;	t0B=a2p04i-a2p05i;
		t08=a2p04r+a2p05r;	t09=a2p04i+a2p05i;
		t0E=a2p06r-a2p07r;	t0F=a2p06i-a2p07i;
		t0C=a2p06r+a2p07r;	t0D=a2p06i+a2p07i;
		rt =t0C;		it =t0D;
		t0C=t08-rt;		t0D=t09-it;
		t08=t08+rt;		t09=t09+it;
		rt =t0E;		it =t0F;
		t0E=t0A-it;		t0F=t0B+rt;
		t0A=t0A+it;		t0B=t0B-rt;
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t08;		it =t09;
		t08=t00-rt;		t09=t01-it;
		t00=t00+rt;		t01=t01+it;
		rt =t0C;		it =t0D;
		t0C=t04-it;		t0D=t05+rt;
		t04=t04+it;		t05=t05-rt;
		rt =(t0A+t0B)*ISRT2;it =(t0A-t0B)*ISRT2;
		t0A=t02-rt;		t0B=t03+it;
		t02=t02+rt;		t03=t03-it;
		rt =(t0E-t0F)*ISRT2;it =(t0F+t0E)*ISRT2;
		t0E=t06+rt;		t0F=t07+it;
		t06=t06-rt;		t07=t07-it;
	/*...Block 2:	*/
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		t12=a2p08r-a2p09r;	t13=a2p08i-a2p09i;
		t10=a2p08r+a2p09r;	t11=a2p08i+a2p09i;
		t16=a2p0Ar-a2p0Br;	t17=a2p0Ai-a2p0Bi;
		t14=a2p0Ar+a2p0Br;	t15=a2p0Ai+a2p0Bi;
		rt =t14;		it =t15;
		t14=t10-rt;		t15=t11-it;
		t10=t10+rt;		t11=t11+it;
		rt =t16;		it =t17;
		t16=t12-it;		t17=t13+rt;
		t12=t12+it;		t13=t13-rt;
		t1A=a2p0Cr-a2p0Dr;	t1B=a2p0Ci-a2p0Di;
		t18=a2p0Cr+a2p0Dr;	t19=a2p0Ci+a2p0Di;
		t1E=a2p0Er-a2p0Fr;	t1F=a2p0Ei-a2p0Fi;
		t1C=a2p0Er+a2p0Fr;	t1D=a2p0Ei+a2p0Fi;
		rt =t1C;		it =t1D;
		t1C=t18-rt;		t1D=t19-it;
		t18=t18+rt;		t19=t19+it;
		rt =t1E;		it =t1F;
		t1E=t1A-it;		t1F=t1B+rt;
		t1A=t1A+it;		t1B=t1B-rt;
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t18;		it =t19;
		t18=t10-rt;		t19=t11-it;
		t10=t10+rt;		t11=t11+it;
		rt =t1C;		it =t1D;
		t1C=t14-it;		t1D=t15+rt;
		t14=t14+it;		t15=t15-rt;
		rt =(t1A+t1B)*ISRT2;it =(t1A-t1B)*ISRT2;
		t1A=t12-rt;		t1B=t13+it;
		t12=t12+rt;		t13=t13-it;
		rt =(t1E-t1F)*ISRT2;it =(t1F+t1E)*ISRT2;
		t1E=t16+rt;		t1F=t17+it;
		t16=t16-rt;		t17=t17-it;
	/*...Block 3:	*/
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		t22=a2p10r-a2p11r;	t23=a2p10i-a2p11i;
		t20=a2p10r+a2p11r;	t21=a2p10i+a2p11i;
		t26=a2p12r-a2p13r;	t27=a2p12i-a2p13i;
		t24=a2p12r+a2p13r;	t25=a2p12i+a2p13i;
		rt =t24;		it =t25;
		t24=t20-rt;		t25=t21-it;
		t20=t20+rt;		t21=t21+it;
		rt =t26;		it =t27;
		t26=t22-it;		t27=t23+rt;
		t22=t22+it;		t23=t23-rt;
		t2A=a2p14r-a2p15r;	t2B=a2p14i-a2p15i;
		t28=a2p14r+a2p15r;	t29=a2p14i+a2p15i;
		t2E=a2p16r-a2p17r;	t2F=a2p16i-a2p17i;
		t2C=a2p16r+a2p17r;	t2D=a2p16i+a2p17i;
		rt =t2C;		it =t2D;
		t2C=t28-rt;		t2D=t29-it;
		t28=t28+rt;		t29=t29+it;
		rt =t2E;		it =t2F;
		t2E=t2A-it;		t2F=t2B+rt;
		t2A=t2A+it;		t2B=t2B-rt;
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t28;		it =t29;
		t28=t20-rt;		t29=t21-it;
		t20=t20+rt;		t21=t21+it;
		rt =t2C;		it =t2D;
		t2C=t24-it;		t2D=t25+rt;
		t24=t24+it;		t25=t25-rt;
		rt =(t2A+t2B)*ISRT2;it =(t2A-t2B)*ISRT2;
		t2A=t22-rt;		t2B=t23+it;
		t22=t22+rt;		t23=t23-it;
		rt =(t2E-t2F)*ISRT2;it =(t2F+t2E)*ISRT2;
		t2E=t26+rt;		t2F=t27+it;
		t26=t26-rt;		t27=t27-it;
	/*...Block 4:	*/
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		t32=a2p18r-a2p19r;	t33=a2p18i-a2p19i;
		t30=a2p18r+a2p19r;	t31=a2p18i+a2p19i;
		t36=a2p1Ar-a2p1Br;	t37=a2p1Ai-a2p1Bi;
		t34=a2p1Ar+a2p1Br;	t35=a2p1Ai+a2p1Bi;
		rt =t34;		it =t35;
		t34=t30-rt;		t35=t31-it;
		t30=t30+rt;		t31=t31+it;
		rt =t36;		it =t37;
		t36=t32-it;		t37=t33+rt;
		t32=t32+it;		t33=t33-rt;
		t3A=a2p1Cr-a2p1Dr;	t3B=a2p1Ci-a2p1Di;
		t38=a2p1Cr+a2p1Dr;	t39=a2p1Ci+a2p1Di;
		t3E=a2p1Er-a2p1Fr;	t3F=a2p1Ei-a2p1Fi;
		t3C=a2p1Er+a2p1Fr;	t3D=a2p1Ei+a2p1Fi;
		rt =t3C;		it =t3D;
		t3C=t38-rt;		t3D=t39-it;
		t38=t38+rt;		t39=t39+it;
		rt =t3E;		it =t3F;
		t3E=t3A-it;		t3F=t3B+rt;
		t3A=t3A+it;		t3B=t3B-rt;
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t38;		it =t39;
		t38=t30-rt;		t39=t31-it;
		t30=t30+rt;		t31=t31+it;
		rt =t3C;		it =t3D;
		t3C=t34-it;		t3D=t35+rt;
		t34=t34+it;		t35=t35-rt;
		rt =(t3A+t3B)*ISRT2;it =(t3A-t3B)*ISRT2;
		t3A=t32-rt;		t3B=t33+it;
		t32=t32+rt;		t33=t33-it;
		rt =(t3E-t3F)*ISRT2;it =(t3F+t3E)*ISRT2;
		t3E=t36+rt;		t3F=t37+it;
		t36=t36-rt;		t37=t37-it;
	//...and now do eight radix-4 transforms, including the internal twiddle factors:
	/*...Block 1: t00,t10,t20,t30	*/
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		rt =t10;	t10=t00-rt;	t00=t00+rt;
		it =t11;	t11=t01-it;	t01=t01+it;
		rt =t30;	t30=t20-rt;	t20=t20+rt;
		it =t31;	t31=t21-it;	t21=t21+it;
		a[rdum   ]=t00+t20;			a[idum   ]=t01+t21;
		t00       =t00-t20;			t01       =t01-t21;
		a[rdum+32]=t00*cB10+t01*sB10;	a[idum+32]=t01*cB10-t00*sB10;
		rt        =t10+t31;			it        =t11-t30;
		t10       =t10-t31;			t11       =t11+t30;
		a[rdum+16]=rt *cB08+it *sB08;	a[idum+16]=it *cB08-rt *sB08;
		a[rdum+48]=t10*cB18+t11*sB18;	a[idum+48]=t11*cB18-t10*sB18;
	/*...Block 5: t08,t18,t28,t38	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t18;	t18=t08-t19;	t08=t08+t19;
					t19=t09+rt;		t09=t09-rt;
		rt =(t29+t28)*ISRT2;		t29=(t29-t28)*ISRT2;		t28=rt;
		rt =(t38-t39)*ISRT2;		it =(t38+t39)*ISRT2;
		t38=t28+rt;					t28=t28-rt;
		t39=t29+it;					t29=t29-it;
		rt        =t08+t28;			it        =t09+t29;
		t08       =t08-t28;			t09       =t09-t29;
		a[rdum+8 ]=rt *cB04+it *sB04;	a[idum+8 ]=it *cB04-rt *sB04;
		a[rdum+40]=t08*cB14+t09*sB14;	a[idum+40]=t09*cB14-t08*sB14;
		rt        =t18+t39;			it        =t19-t38;
		t18       =t18-t39;			t19       =t19+t38;
		a[rdum+24]=rt *cB0C+it *sB0C;	a[idum+24]=it *cB0C-rt *sB0C;
		a[rdum+56]=t18*cB1C+t19*sB1C;	a[idum+56]=t19*cB1C-t18*sB1C;
	/*...Block 3: t04,t14,t24,t34	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		rt =(t15+t14)*ISRT2;		it =(t15-t14)*ISRT2;
		t14=t04-rt;					t04=t04+rt;
		t15=t05-it;					t05=t05+it;
		rt =t24*c + t25*s;			t25=t25*c - t24*s;		t24=rt;
		rt =t34*s + t35*c;			it =t35*s - t34*c;
		t34=t24-rt;					t24=t24+rt;
		t35=t25-it;					t25=t25+it;
		rt        =t04+t24;			it        =t05+t25;
		t04       =t04-t24;			t05       =t05-t25;
		a[rdum+4 ]=rt *cB02+it *sB02;	a[idum+4 ]=it *cB02-rt *sB02;
		a[rdum+36]=t04*cB12+t05*sB12;	a[idum+36]=t05*cB12-t04*sB12;
		rt        =t14+t35;			it        =t15-t34;
		t14       =t14-t35;			t15       =t15+t34;
		a[rdum+20]=rt *cB0A+it *sB0A;	a[idum+20]=it *cB0A-rt *sB0A;
		a[rdum+52]=t14*cB1A+t15*sB1A;	a[idum+52]=t15*cB1A-t14*sB1A;
	/*...Block 7: t0C,t1C,t2C,t3C	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =(t1C-t1D)*ISRT2;		it =(t1C+t1D)*ISRT2;
		t1C=t0C+rt;					t0C=t0C-rt;
		t1D=t0D+it;					t0D=t0D-it;
		rt =t2C*s + t2D*c;			t2D=t2D*s - t2C*c;		t2C=rt;
		rt =t3C*c + t3D*s;			it =t3D*c - t3C*s;
		t3C=t2C+rt;					t2C=t2C-rt;
		t3D=t2D+it;					t2D=t2D-it;
		rt        =t0C+t2C;			it        =t0D+t2D;
		t0C       =t0C-t2C;			t0D       =t0D-t2D;
		a[rdum+12]=rt *cB06+it *sB06;	a[idum+12]=it *cB06-rt *sB06;
		a[rdum+44]=t0C*cB16+t0D*sB16;	a[idum+44]=t0D*cB16-t0C*sB16;
		rt        =t1C+t3D;			it        =t1D-t3C;
		t1C       =t1C-t3D;			t1D       =t1D+t3C;
		a[rdum+28]=rt *cB0E+it *sB0E;	a[idum+28]=it *cB0E-rt *sB0E;
		a[rdum+60]=t1C*cB1E+t1D*sB1E;	a[idum+60]=t1D*cB1E-t1C*sB1E;
	/*...Block 2: t02,t12,t22,t32	*/
	#ifdef USE_AVX512
		rdum += 5;	idum += 5;
	#elif defined(USE_AVX)
		++rdum;	++idum;
	#elif defined(USE_SSE2)
		--rdum;	--idum;
	#endif
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		rt =t12*c + t13*s;			it =t13*c - t12*s;
		t12=t02-rt;					t02=t02+rt;
		t13=t03-it;					t03=t03+it;
		rt =t22*c32_1 + t23*s32_1;	t23=t23*c32_1 - t22*s32_1;	t22=rt;
		rt =t32*c32_3 + t33*s32_3;	it =t33*c32_3 - t32*s32_3;
		t32=t22-rt;					t22=t22+rt;
		t33=t23-it;					t23=t23+it;
		rt        =t02+t22;			it        =t03+t23;
		t02       =t02-t22;			t03       =t03-t23;
		a[rdum+2 ]=rt *cB01+it *sB01;	a[idum+2 ]=it *cB01-rt *sB01;
		a[rdum+34]=t02*cB11+t03*sB11;	a[idum+34]=t03*cB11-t02*sB11;
		rt        =t12+t33;			it        =t13-t32;
		t12       =t12-t33;			t13       =t13+t32;
		a[rdum+18]=rt *cB09+it *sB09;	a[idum+18]=it *cB09-rt *sB09;
		a[rdum+50]=t12*cB19+t13*sB19;	a[idum+50]=t13*cB19-t12*sB19;
	/*...Block 6: t0A,t1A,t2A,t3A	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t1A*s - t1B*c;			it =t1B*s + t1A*c;
		t1A=t0A+rt;					t0A =t0A-rt;
		t1B=t0B+it;					t0B =t0B-it;
		rt =t2A*s32_3 + t2B*c32_3;	t2B=t2B*s32_3 - t2A*c32_3;	t2A=rt;
		rt =t3A*c32_1 - t3B*s32_1;	it =t3B*c32_1 + t3A*s32_1;
		t3A=t2A+rt;					t2A=t2A-rt;
		t3B=t2B+it;					t2B=t2B-it;
		rt        =t0A+t2A;			it        =t0B+t2B;
		t0A       =t0A-t2A;			t0B       =t0B-t2B;
		a[rdum+10]=rt *cB05+it *sB05;	a[idum+10]=it *cB05-rt *sB05;
		a[rdum+42]=t0A*cB15+t0B*sB15;	a[idum+42]=t0B*cB15-t0A*sB15;
		rt        =t1A+t3B;			it        =t1B-t3A;
		t1A       =t1A-t3B;			t1B       =t1B+t3A;
		a[rdum+26]=rt *cB0D+it *sB0D;	a[idum+26]=it *cB0D-rt *sB0D;
		a[rdum+58]=t1A*cB1D+t1B*sB1D;	a[idum+58]=t1B*cB1D-t1A*sB1D;
	/*...Block 4: t06,t16,t26,t36	*/
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
	#if PFETCH
	prefetch_p_doubles(addr);
	addr += 4;
	#endif
		rt =t16*s + t17*c;			it =t17*s - t16*c;
		t16=t06-rt;					t06 =t06+rt;
		t17=t07-it;					t07 =t07+it;
		rt =t26*c32_3 + t27*s32_3;	t27=t27*c32_3 - t26*s32_3;	t26=rt;
		rt =t36*s32_1 - t37*c32_1;	it =t37*s32_1 + t36*c32_1;
		t36=t26+rt;					t26=t26-rt;
		t37=t27+it;					t27=t27-it;
		rt        =t06+t26;			it        =t07+t27;
		t06       =t06-t26;			t07       =t07-t27;
		a[rdum+6 ]=rt *cB03+it *sB03;	a[idum+6 ]=it *cB03-rt *sB03;
		a[rdum+38]=t06*cB13+t07*sB13;	a[idum+38]=t07*cB13-t06*sB13;
		rt        =t16+t37;			it        =t17-t36;
		t16       =t16-t37;			t17       =t17+t36;
		a[rdum+22]=rt *cB0B+it *sB0B;	a[idum+22]=it *cB0B-rt *sB0B;
		a[rdum+54]=t16*cB1B+t17*sB1B;	a[idum+54]=t17*cB1B-t16*sB1B;
	/*...Block 8: t0E,t1E,t2E,t3E	*/
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(addr);
		#endif
	addr += 4;
	#endif
		rt =t1E*c - t1F*s;			it =t1F*c + t1E*s;
		t1E=t0E+rt;					t0E =t0E-rt;
		t1F=t0F+it;					t0F =t0F-it;
		rt =t2E*s32_1 + t2F*c32_1;	t2F=t2F*s32_1 - t2E*c32_1;	t2E=rt;
		rt =t3E*s32_3 + t3F*c32_3;	it =t3F*s32_3 - t3E*c32_3;
		t3E=t2E+rt;					t2E=t2E-rt;
		t3F=t2F+it;					t2F=t2F-it;
		rt        =t0E+t2E;			it        =t0F+t2F;
		t0E       =t0E-t2E;			t0F       =t0F-t2F;
		a[rdum+14]=rt *cB07+it *sB07;	a[idum+14]=it *cB07-rt *sB07;
		a[rdum+46]=t0E*cB17+t0F*sB17;	a[idum+46]=t0F*cB17-t0E*sB17;
		rt        =t1E+t3F;			it        =t1F-t3E;
		t1E       =t1E-t3F;			t1F       =t1F+t3E;
		a[rdum+30]=rt *cB0F+it *sB0F;	a[idum+30]=it *cB0F-rt *sB0F;
		a[rdum+62]=t1E*cB1F+t1F*sB1F;	a[idum+62]=t1F*cB1F-t1E*sB1F;

#ifdef USE_SSE2

	} else {	// (j1 == 0) = false:

		re = c01->d0;	im = (c01+1)->d0;	// Save copies of these 2 for wrapper step

	/* [cf. the radix-16 analog if this file for notes on the AVX data-alyout]
	In sse2 mode, the input data are laid out in memory like so, where we view things in 16-byte chunks:

		&a[j1]	a0.re,a1.re		&a[j2]	b0.re,b1.re
		+0x010	a0.im,a1.im		+0x010	b0.im,b1.im
		+0x020	a2.re,a3.re		+0x020	b2.re,b3.re
		+0x030	a2.im,a3.im		+0x030	b2.im,b3.im
		+0x040	a4.re,a5.re		+0x040	b4.re,b5.re
		+0x050	a4.im,a5.im		+0x050	b4.im,b5.im
		...

	We need to interleave these pairwise so as to swap the high word of each A-pair with the low word of the corresponding B-pair, e.g

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
	*/
	if(fwd_fft_only != 3)	// v20: add support for both-inputs-already-fwd-FFTed case
	{
	  #ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX mode
	  					// because (add[1,3,5,7]-add[0,2,4,6]) have opposite signs for Fermat and Mersenne-mod:
		// process 8 main-array blocks of 8 vec_dbl [= 8 x 8 = 64 doubles each] in AVX mode, total = 512 float64
		SSE2_RADIX32_WRAPPER_DIF(add0,add1,add2,add3,add4,add5,add6,add7
													,r00,r10,r20,r30,isrt2,cc0,c00,c01,c02,c03,c05,c07)
	  #elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  					// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:
		// process 4 main-array blocks of 16 vec_dbl [= 16 x 4 = 64 doubles each] in AVX mode, total = 256 float64
		SSE2_RADIX32_WRAPPER_DIF(add0,add1,add2,add3,r00,r10,r20,r30,isrt2,cc0,c00,c01,c02,c03,c05,c07)
	  #else	// SSE2:
		SSE2_RADIX32_WRAPPER_DIF(add0,add1,          r00,r10,r20,r30,isrt2,cc0,c00,c01,c02,c03,c05,c07)
	  #endif
	}

	// v19: If fwd_fft_only = 1, write fwd-FFT result back to input array, skipping dyadic-square and inv-FFT steps:
	if(fwd_fft_only == 1)
	{
	  #ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX512 mode
		// process 8 main-array blocks of [8 vec_dbl = 8 x 8 = 64 doubles] each in AVX512 mode, total = 64 vec_dbl = 32 vec_cmplx
		nbytes = (intptr_t)r08 - (intptr_t)r00;
		memcpy(add0, r00, nbytes);	// add0 = a + j1pad;
		memcpy(add2, r10, nbytes);	// add2 = add0 + 64;	// add2 = add0 + [64 doubles, equiv to 8 AVX-512 registers]
		memcpy(add4, r20, nbytes);	// add4 = add2 + 64;
		memcpy(add6, r30, nbytes);	// add6 = add4 + 64;
		memcpy(add1, r08, nbytes);	// add1 = a + j2pad;
		memcpy(add3, r18, nbytes);	// add3 = add1 - 64;	// Last 4 offsets run in descending order for Mers-mod
		memcpy(add5, r28, nbytes);	// add5 = add3 - 64;
		memcpy(add7, r38, nbytes);	// add7 = add5 - 64;
	  #elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  						// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:
		// process 4 main-array blocks of [16 vec_dbl = 16 x 4 = 64 doubles] each in AVX mode, total = 64 vec_dbl = 32 vec_cmplx
		nbytes = (intptr_t)r10 - (intptr_t)r00;
		memcpy(add0, r00, nbytes);	// add0 = a + j1pad;
		memcpy(add1, r10, nbytes);	// add1 = a + j2pad;
		memcpy(add2, r20, nbytes);	// add2 = add0 + 64;	// add2 = add0 + [64 doubles, equiv to 8 AVX registers]
		memcpy(add3, r30, nbytes);	// add3 = add1 - 64;	// Last 2 offsets run in descending order for Mers-mod
	  #else	// SSE2:
		// process 2 main-array blocks of [32 vec_dbl = 32 x 2 = 64 doubles] each in SSE2 mode, total = 64 vec_dbl = 32 vec_cmplx
		nbytes = (intptr_t)r20 - (intptr_t)r00;
		memcpy(add0, r00, nbytes);	// add0 = a + j1pad;
		memcpy(add1, r20, nbytes);	// add1 = a + j2pad;
	  #endif
		goto loop;	// Skip dyadic-mul and iFFT

	} else {
	  // Structure remaining 2 options as nested if() subclause to allow sharing of several hundred lines of sincos-multiplier-related code:
	  if (fwd_fft_only) {	// 2-input modmul: fwdFFT data dyadic-mul'ed with precomputed 2nd-vector stored in fwdFFTed form at address stored in (uint64)-cast form in fwd_fft_only

	   if(fwd_fft_only == 3)	// v20: Both inputs enter fwd-FFTed, must copy from main arrays a[],b[] to local-vars
	   {
	   #ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX512 mode
		// process 8 main-array blocks of [8 vec_dbl = 8 x 8 = 64 doubles] each in AVX512 mode, total = 64 vec_dbl = 32 vec_cmplx
		nbytes = (intptr_t)r08 - (intptr_t)r00;
		memcpy(r00, add0, nbytes);	// add0 = a + j1pad;
		memcpy(r10, add2, nbytes);	// add2 = add0 + 64;	// add2 = add0 + [64 doubles, equiv to 8 AVX-512 registers]
		memcpy(r20, add4, nbytes);	// add4 = add2 + 64;
		memcpy(r30, add6, nbytes);	// add6 = add4 + 64;
		memcpy(r08, add1, nbytes);	// add1 = a + j2pad;
		memcpy(r18, add3, nbytes);	// add3 = add1 - 64;	// Last 4 offsets run in descending order for Mers-mod
		memcpy(r28, add5, nbytes);	// add5 = add3 - 64;
		memcpy(r38, add7, nbytes);	// add7 = add5 - 64;
	   #elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  						// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:
		// process 4 main-array blocks of [16 vec_dbl = 16 x 4 = 64 doubles] each in AVX mode, total = 64 vec_dbl = 32 vec_cmplx
		nbytes = (intptr_t)r10 - (intptr_t)r00;
		memcpy(r00, add0, nbytes);	// add0 = a + j1pad;
		memcpy(r10, add1, nbytes);	// add1 = a + j2pad;
		memcpy(r20, add2, nbytes);	// add2 = add0 + 64;	// add2 = add0 + [64 doubles, equiv to 8 AVX registers]
		memcpy(r30, add3, nbytes);	// add3 = add1 - 64;	// Last 2 offsets run in descending order for Mers-mod
	   #else	// SSE2:
		// process 2 main-array blocks of [32 vec_dbl = 32 x 2 = 64 doubles] each in SSE2 mode, total = 64 vec_dbl = 32 vec_cmplx
		nbytes = (intptr_t)r20 - (intptr_t)r00;
		memcpy(r00, add0, nbytes);	// add0 = a + j1pad;
		memcpy(r20, add1, nbytes);	// add1 = a + j2pad;
	   #endif
	   }

	  #ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX512 mode
		// process 8 main-array blocks of [8 vec_dbl = 8 x 8 = 64 doubles] each in AVX512 mode, total = 64 vec_dbl = 32 vec_cmplx
		bdd0 = b + j1pad;
		bdd2 = bdd0 + 64;	// bdd2 = bdd0 + [64 doubles, equiv to 4 AVX-512 registers]
		bdd4 = bdd2 + 64;
		bdd6 = bdd4 + 64;
		bdd1 = b + j2pad;
		bdd3 = bdd1 - 64;	// Last 4 offsets run in descending order for Mers-mod
		bdd5 = bdd3 - 64;
		bdd7 = bdd5 - 64;
	  #elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  						// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:
		// process 4 main-array blocks of [16 vec_dbl = 16 x 4 = 64 doubles] each in AVX mode, total = 64 vec_dbl = 32 vec_cmplx
		bdd0 = b + j1pad;
		bdd1 = b + j2pad;
		bdd2 = bdd0 + 64;	// bdd2 = bdd0 + [64 doubles, equiv to 8 AVX registers]
		bdd3 = bdd1 - 64;	// Last 2 offsets run in descending order for Mers-mod
	  #else	// SSE2:
		// process 2 main-array blocks of [32 vec_dbl = 32 x 2 = 64 doubles] each in SSE2 mode, total = 64 vec_dbl = 32 vec_cmplx
		bdd0 = b + j1pad;
		bdd1 = b + j2pad;
	  #endif

		if(c_arr) {	// c_arr != 0x0: a * (b - c)

		#ifdef USE_AVX512
			cdd0 = c_arr + j1pad;
			cdd2 = cdd0 + 64;
			cdd4 = cdd2 + 64;
			cdd6 = cdd4 + 64;
			cdd1 = c_arr + j2pad;
			cdd3 = cdd1 - 64;
			cdd5 = cdd3 - 64;
			cdd7 = cdd5 - 64;
		#elif defined(USE_AVX)
			cdd0 = c_arr + j1pad;
			cdd1 = c_arr + j2pad;
			cdd2 = cdd0 + 64;
			cdd3 = cdd1 - 64;
		#else	// SSE2:
			cdd0 = c_arr + j1pad;
			cdd1 = c_arr + j2pad;
		#endif
		/********** In the inline-(B-C) macros, r1 and bminusc_ptr are the I/O base-addresses;
				bdd0,cdd0 the B and C-array-section addresses: **********/
		#ifdef USE_ARM_V8_SIMD
		  __asm__ volatile (\
			"ldr	x0,%[__bminusc_ptr]	\n\t"\
			"ldr	x1,%[__bdd0]		\n\t	ldr	x2,%[__cdd0]	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t"\
			"ldr	x1,%[__bdd1]		\n\t	ldr	x2,%[__cdd1]\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			"add x0,x0,#0x40 \n\t add x1,x1,#0x40 \n\t add x2,x2,#0x40	\n\t"\
			"ldp	q0,q1,[x1]			\n\t	ldp	q2,q3,[x1,#0x20]	\n\t"\
			"ldp	q4,q5,[x2]			\n\t	ldp	q6,q7,[x2,#0x20]	\n\t"\
			"fsub	v0.2d,v0.2d,v4.2d	\n\t	fsub	v1.2d,v1.2d,v5.2d	\n\t"\
			"fsub	v2.2d,v2.2d,v6.2d	\n\t	fsub	v3.2d,v3.2d,v7.2d	\n\t"\
			"stp	q0,q1,[x0]			\n\t	stp	q2,q3,[x0,#0x20]	\n\t"\
			:					// outputs: none
			: [__bminusc_ptr] "m" (bminusc_ptr)	// Local temp-storage array to hold B-C outputs
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			 ,[__bdd1] "m" (bdd1)
			 ,[__cdd0] "m" (cdd0)
			 ,[__cdd1] "m" (cdd1)
			: "cc","memory","x0","x1","x2", "v0","v1","v2","v3","v4","v5","v6","v7"	/* Clobbered registers */\
		  );
		  // Re-point bdd* to bminusc_ptr-offsets - bdd are double*, thus each bminusc_ptr-incr = 2 doubles,
		  // total local-table size = 2*radix = 32 vec_dbl = 64 double. Within this contig-table bdd* indices are in-order:
		  bdd0 = (double *)bminusc_ptr;
		  bdd1 = bdd0 + 64;

		#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

		  __asm__ volatile (\
		"movq	%[__bminusc_ptr],%%rax	\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t	movq	%[__cdd0],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax	\n\t	addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"movq	%[__bdd1],%%rbx		\n\t	movq	%[__cdd1],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax	\n\t	addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"movq	%[__bdd2],%%rbx		\n\t	movq	%[__cdd2],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax	\n\t	addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"movq	%[__bdd3],%%rbx		\n\t	movq	%[__cdd3],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax	\n\t	addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"movq	%[__bdd4],%%rbx		\n\t	movq	%[__cdd4],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax	\n\t	addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"movq	%[__bdd5],%%rbx		\n\t	movq	%[__cdd5],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax	\n\t	addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"movq	%[__bdd6],%%rbx		\n\t	movq	%[__cdd6],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax	\n\t	addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax		\n\t"\
		"movq	%[__bdd7],%%rbx		\n\t	movq	%[__cdd7],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
		"addq	$0x100,%%rax	\n\t	addq	$0x100,%%rbx	\n\t	addq	$0x100,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%zmm0	\n\t	vmovaps	0x080(%%rbx),%%zmm2	\n\t"\
			"vmovaps	0x040(%%rbx),%%zmm1	\n\t	vmovaps	0x0c0(%%rbx),%%zmm3	\n\t"\
		"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t"\
		"vsubpd	0x040(%%rcx),%%zmm1,%%zmm1	\n\t	vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t"\
			"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm2,0x080(%%rax)	\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm3,0x0c0(%%rax)	\n\t"\
			:					// outputs: none
			: [__bminusc_ptr] "m" (bminusc_ptr)	// Local temp-storage array to hold B-C outputs
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			 ,[__bdd1] "m" (bdd1)
			 ,[__bdd2] "m" (bdd2)
			 ,[__bdd3] "m" (bdd3)
			 ,[__bdd4] "m" (bdd4)
			 ,[__bdd5] "m" (bdd5)
			 ,[__bdd6] "m" (bdd6)
			 ,[__bdd7] "m" (bdd7)
			 ,[__cdd0] "m" (cdd0)
			 ,[__cdd1] "m" (cdd1)
			 ,[__cdd2] "m" (cdd2)
			 ,[__cdd3] "m" (cdd3)
			 ,[__cdd4] "m" (cdd4)
			 ,[__cdd5] "m" (cdd5)
			 ,[__cdd6] "m" (cdd6)
			 ,[__cdd7] "m" (cdd7)
			: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3"	// Clobbered registers
		  );
		  // Re-point bdd* to bminusc_ptr-offsets - bdd are double*, thus each bminusc_ptr-incr = 8 doubles,
		  // total local-table size = 2*radix = 32 vec_dbl = 256 double. Within this contig-table bdd* indices are in-order:
		  bdd0 = (double *)bminusc_ptr;
		  bdd1 = bdd0 + 64;
		  bdd2 = bdd1 + 64;
		  bdd3 = bdd2 + 64;
		  bdd4 = bdd3 + 64;
		  bdd5 = bdd4 + 64;
		  bdd6 = bdd5 + 64;
		  bdd7 = bdd6 + 64;

		#elif defined(USE_AVX)

		  __asm__ volatile (\
		"movq	%[__bminusc_ptr],%%rax	\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t	movq	%[__cdd0],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x080,%%rax		\n\t"\
		"movq	%[__bdd1],%%rbx		\n\t	movq	%[__cdd1],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x080,%%rax		\n\t"\
		"movq	%[__bdd2],%%rbx		\n\t	movq	%[__cdd2],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x080,%%rax		\n\t"\
		"movq	%[__bdd3],%%rbx		\n\t	movq	%[__cdd3],%%rcx		\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
		"addq	$0x80,%%rax	\n\t	addq	$0x80,%%rbx	\n\t	addq	$0x80,%%rcx	\n\t"\
			"vmovaps		 (%%rbx),%%ymm0	\n\t	vmovaps	0x040(%%rbx),%%ymm2	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t	vmovaps	0x060(%%rbx),%%ymm3	\n\t"\
		"vsubpd		 (%%rcx),%%ymm0,%%ymm0	\n\t	vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t"\
		"vsubpd	0x020(%%rcx),%%ymm1,%%ymm1	\n\t	vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm2,0x040(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm3,0x060(%%rax)	\n\t"\
			:					// outputs: none
			: [__bminusc_ptr] "m" (bminusc_ptr)	// Local temp-storage array to hold B-C outputs
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			 ,[__bdd1] "m" (bdd1)
			 ,[__bdd2] "m" (bdd2)
			 ,[__bdd3] "m" (bdd3)
			 ,[__cdd0] "m" (cdd0)
			 ,[__cdd1] "m" (cdd1)
			 ,[__cdd2] "m" (cdd2)
			 ,[__cdd3] "m" (cdd3)
			: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3"	// Clobbered registers
		  );
		  // Re-point bdd* to bminusc_ptr-offsets - bdd are double*, thus each bminusc_ptr-incr = 4 doubles,
		  // total local-table size = 2*radix = 32 vec_dbl = 128 double. Within this contig-table bdd* indices are in-order:
		  bdd0 = (double *)bminusc_ptr;
		  bdd1 = bdd0 + 64;
		  bdd2 = bdd1 + 64;
		  bdd3 = bdd2 + 64;

		#else	// 64-bit SSE2:

		  __asm__ volatile (\
		"movq	%[__bminusc_ptr],%%rax	\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t	movq	%[__cdd0],%%rcx		\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x040,%%rax		\n\t"\
		"movq	%[__bdd1],%%rbx		\n\t	movq	%[__cdd1],%%rcx		\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
		"addq	$0x40,%%rax	\n\t	addq	$0x40,%%rbx	\n\t	addq	$0x40,%%rcx	\n\t"\
			"movaps		 (%%rbx),%%xmm0	\n\t	movaps	0x020(%%rbx),%%xmm2	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t	movaps	0x030(%%rbx),%%xmm3	\n\t"\
			"subpd		 (%%rcx),%%xmm0	\n\t	subpd	0x020(%%rcx),%%xmm2	\n\t"\
			"subpd	0x010(%%rcx),%%xmm1	\n\t	subpd	0x030(%%rcx),%%xmm3	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm2,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm3,0x030(%%rax)	\n\t"\
			:					// outputs: none
			: [__bminusc_ptr] "m" (bminusc_ptr)	// Local temp-storage array to hold B-C outputs
			 ,[__bdd0] "m" (bdd0)	// B-array base-address
			 ,[__bdd1] "m" (bdd1)
			 ,[__cdd0] "m" (cdd0)
			 ,[__cdd1] "m" (cdd1)
			: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3"	// Clobbered registers
		  );
		  // Re-point bdd* to bminusc_ptr-offsets - bdd are double*, thus each bminusc_ptr-incr = 2 doubles,
		  // total local-table size = 2*radix = 32 vec_dbl = 64 double. Within this contig-table bdd* indices are in-order:
		  bdd0 = (double *)bminusc_ptr;
		  bdd1 = bdd0 + 64;

		#endif	// AVX or SSE2?

		} else {	// c_arr = 0x0: a * b

		}	// endif(c_arr)

		/* In SSE2 mode, data laid out in memory as described in "Normal execution" section below, with b-inputs similar to
		a-inputs, just with r00-31 and r32-63 replaced by offsets relative to b-array addresses bdd0 and bdd1, respectively:
							A-data:									B-data[SSE2]:				B-data[AVX]:			B-data[AVX-512]:
			r00:a1p00r,a2p00r		r01:a1p00i,a2p00i			bdd0+0x000,bdd0+0x010		bdd0+0x000,bdd0+0x020		bdd0+0x000,bdd0+0x040
			r02:a1p10r,a2p10r		r03:a1p10i,a2p10i			bdd0+0x020,bdd0+0x030		bdd0+0x040,bdd0+0x060		bdd0+0x080,bdd0+0x0c0
			r04:a1p08r,a2p08r		r05:a1p08i,a2p08i			bdd0+0x040,bdd0+0x050		bdd0+0x080,bdd0+0x0a0		bdd0+0x100,bdd1+0x140
			r06:a1p18r,a2p18r		r07:a1p18i,a2p18i			bdd0+0x060,bdd0+0x070		bdd0+0x0c0,bdd0+0x0e0		bdd0+0x180,bdd1+0x1c0
			r08:a1p04r,a2p04r		r09:a1p04i,a2p04i			bdd0+0x080,bdd0+0x090		bdd0+0x100,bdd0+0x120		bdd1+0x000,bdd2+0x040
			r0A:a1p14r,a2p14r		r0B:a1p14i,a2p14i			bdd0+0x0a0,bdd0+0x0b0		bdd0+0x140,bdd0+0x160		bdd1+0x080,bdd2+0x0c0
			r0C:a1p0Cr,a2p0Cr		r0D:a1p0Ci,a2p0Ci			bdd0+0x0c0,bdd0+0x0d0		bdd0+0x180,bdd0+0x1a0		bdd1+0x100,bdd3+0x140
			r0E:a1p1Cr,a2p1Cr		r0F:a1p1Ci,a2p1Ci			bdd0+0x0e0,bdd0+0x0f0		bdd0+0x1c0,bdd0+0x1e0		bdd1+0x180,bdd3+0x1c0
			r10:a1p02r,a2p02r		r11:a1p02i,a2p02i			bdd0+0x100,bdd0+0x110		bdd1+0x000,bdd1+0x020		bdd2+0x000,bdd4+0x040
			r12:a1p12r,a2p12r		r13:a1p12i,a2p12i			bdd0+0x120,bdd0+0x130		bdd1+0x040,bdd1+0x060		bdd2+0x080,bdd4+0x0c0
			r14:a1p0Ar,a2p0Ar		r15:a1p0Ai,a2p0Ai			bdd0+0x140,bdd0+0x150		bdd1+0x080,bdd1+0x0a0		bdd2+0x100,bdd5+0x140
			r16:a1p1Ar,a2p1Ar		r17:a1p1Ai,a2p1Ai			bdd0+0x160,bdd0+0x170		bdd1+0x0c0,bdd1+0x0e0		bdd2+0x180,bdd5+0x1c0
			r18:a1p06r,a2p06r		r19:a1p06i,a2p06i			bdd0+0x180,bdd0+0x190		bdd1+0x100,bdd1+0x120		bdd3+0x000,bdd6+0x040
			r1A:a1p16r,a2p16r		r1B:a1p16i,a2p16i			bdd0+0x1a0,bdd0+0x1b0		bdd1+0x140,bdd1+0x160		bdd3+0x080,bdd6+0x0c0
			r1C:a1p0Er,a2p0Er		r1D:a1p0Ei,a2p0Ei			bdd0+0x1c0,bdd0+0x1d0		bdd1+0x180,bdd1+0x1a0		bdd3+0x100,bdd7+0x140
			r1E:a1p1Er,a2p1Er		r1F:a1p1Ei,a2p1Ei .			bdd0+0x1e0,bdd0+0x1f0		bdd1+0x1c0,bdd1+0x1e0		bdd3+0x180,bdd7+0x1c0
			r20:a1p01r,a2p01r		r21:a1p01i,a2p01i			bdd1+0x000,bdd1+0x010		bdd2+0x000,bdd2+0x020		bdd4+0x000,bdd0+0x040
			r22:a1p11r,a2p11r		r23:a1p11i,a2p11i			bdd1+0x020,bdd1+0x030		bdd2+0x040,bdd2+0x060		bdd4+0x080,bdd0+0x0c0
			r24:a1p09r,a2p09r		r25:a1p09i,a2p09i			bdd1+0x040,bdd1+0x050		bdd2+0x080,bdd2+0x0a0		bdd4+0x100,bdd1+0x140
			r26:a1p19r,a2p19r		r27:a1p19i,a2p19i			bdd1+0x060,bdd1+0x070		bdd2+0x0c0,bdd2+0x0e0		bdd4+0x180,bdd1+0x1c0
			r28:a1p05r,a2p05r		r29:a1p05i,a2p05i			bdd1+0x080,bdd1+0x090		bdd2+0x100,bdd2+0x120		bdd5+0x000,bdd2+0x040
			r2A:a1p15r,a2p15r		r2B:a1p15i,a2p15i			bdd1+0x0a0,bdd1+0x0b0		bdd2+0x140,bdd2+0x160		bdd5+0x080,bdd2+0x0c0
			r2C:a1p0Dr,a2p0Dr		r2D:a1p0Di,a2p0Di			bdd1+0x0c0,bdd1+0x0d0		bdd2+0x180,bdd2+0x1a0		bdd5+0x100,bdd3+0x140
			r2E:a1p1Dr,a2p1Dr		r2F:a1p1Di,a2p1Di			bdd1+0x0e0,bdd1+0x0f0		bdd2+0x1c0,bdd2+0x1e0		bdd5+0x180,bdd3+0x1c0
			r30:a1p03r,a2p03r		r31:a1p03i,a2p03i			bdd1+0x100,bdd1+0x110		bdd3+0x000,bdd3+0x020		bdd6+0x000,bdd4+0x040
			r32:a1p13r,a2p13r		r33:a1p13i,a2p13i			bdd1+0x120,bdd1+0x130		bdd3+0x040,bdd3+0x060		bdd6+0x080,bdd4+0x0c0
			r34:a1p0Br,a2p0Br		r35:a1p0Bi,a2p0Bi			bdd1+0x140,bdd1+0x150		bdd3+0x080,bdd3+0x0a0		bdd6+0x100,bdd5+0x140
			r36:a1p1Br,a2p1Br		r37:a1p1Bi,a2p1Bi			bdd1+0x160,bdd1+0x170		bdd3+0x0c0,bdd3+0x0e0		bdd6+0x180,bdd5+0x1c0
			r38:a1p07r,a2p07r		r39:a1p07i,a2p07i			bdd1+0x180,bdd1+0x190		bdd3+0x100,bdd3+0x120		bdd7+0x000,bdd6+0x040
			r3A:a1p17r,a2p17r		r3B:a1p17i,a2p17i			bdd1+0x1a0,bdd1+0x1b0		bdd3+0x140,bdd3+0x160		bdd7+0x080,bdd6+0x0c0
			r3C:a1p0Fr,a2p0Fr		r3D:a1p0Fi,a2p0Fi			bdd1+0x1c0,bdd1+0x1d0		bdd3+0x180,bdd3+0x1a0		bdd7+0x100,bdd7+0x140
			r3E:a1p1Fr,a2p1Fr		r3F:a1p1Fi,a2p1Fi .			bdd1+0x1e0,bdd1+0x1f0		bdd3+0x1c0,bdd3+0x1e0		bdd7+0x180,bdd7+0x1c0
		The modified call sequence below takes advantage of that, by processing data which are in 8 XMM registers in-place.
		*/
		// Premultiply RT,IT by 0.25 here to avoid awkward mul-by-0.25s in PAIR_MUL_4 macros:
		dmult = 0.25;
	  } else {
		dmult = 1.0;
	  }	// endif(fwd_fft_only != 0 or == 0)

		re *= dmult;	im *= dmult;
		t01=(re-im)*ISRT2;	t02=(re+im)*ISRT2;

		re0=re*c;	im0=im*s;	re1=re*s;	im1=im*c;
		t03=re0-im0;	t04=re1+im1;
		t05=re1-im1;	t06=re0+im0;

		re0=re*c32_1;	im0=im*s32_1;	re1=re*s32_1;	im1=im*c32_1;
		t07=re0-im0;	t08=re1+im1;
		t09=re1-im1;	t0A=re0+im0;

		re0=re*c32_3;	im0=im*s32_3;	re1=re*s32_3;	im1=im*c32_3;
		t0B=re1-im1;	t0C=re0+im0;
		t0D=re0-im0;	t0E=re1+im1;

		tmp0->d0 = re;		tmp1->d0 = im;
		tmp0->d1 =t0A;		tmp1->d1 =t09;

		tmp2->d0 =t01;		tmp3->d0 =t02;
		tmp2->d1 =t0E;		tmp3->d1 =t0D;

		tmp4->d0 =t03;		tmp5->d0 =t04;
		tmp4->d1 =t0C;		tmp5->d1 =t0B;

		tmp6->d0 =t05;		tmp7->d0 =t06;
		tmp6->d1 =t08;		tmp7->d1 =t07;

	  #ifdef USE_AVX
		re = dmult*c01->d2;	im = dmult*(c01+1)->d2;	// Notice we only use even-indexed ->d* values to set RT,IT in SIMD mode
		t01=(re-im)*ISRT2;	t02=(re+im)*ISRT2;

		re0=re*c;	im0=im*s;	re1=re*s;	im1=im*c;
		t03=re0-im0;	t04=re1+im1;
		t05=re1-im1;	t06=re0+im0;

		re0=re*c32_1;	im0=im*s32_1;	re1=re*s32_1;	im1=im*c32_1;
		t07=re0-im0;	t08=re1+im1;
		t09=re1-im1;	t0A=re0+im0;

		re0=re*c32_3;	im0=im*s32_3;	re1=re*s32_3;	im1=im*c32_3;
		t0B=re1-im1;	t0C=re0+im0;
		t0D=re0-im0;	t0E=re1+im1;

		tmp0->d2 = re;		tmp1->d2 = im;
		tmp0->d3 =t0A;		tmp1->d3 =t09;

		tmp2->d2 =t01;		tmp3->d2 =t02;
		tmp2->d3 =t0E;		tmp3->d3 =t0D;

		tmp4->d2 =t03;		tmp5->d2 =t04;
		tmp4->d3 =t0C;		tmp5->d3 =t0B;

		tmp6->d2 =t05;		tmp7->d2 =t06;
		tmp6->d3 =t08;		tmp7->d3 =t07;
	  #endif

	  #ifdef USE_AVX512
		re = dmult*c01->d4;	im = dmult*(c01+1)->d4;	// AVX-512 uses d4...
		t01=(re-im)*ISRT2;	t02=(re+im)*ISRT2;

		re0=re*c;	im0=im*s;	re1=re*s;	im1=im*c;
		t03=re0-im0;	t04=re1+im1;
		t05=re1-im1;	t06=re0+im0;

		re0=re*c32_1;	im0=im*s32_1;	re1=re*s32_1;	im1=im*c32_1;
		t07=re0-im0;	t08=re1+im1;
		t09=re1-im1;	t0A=re0+im0;

		re0=re*c32_3;	im0=im*s32_3;	re1=re*s32_3;	im1=im*c32_3;
		t0B=re1-im1;	t0C=re0+im0;
		t0D=re0-im0;	t0E=re1+im1;

		tmp0->d4 = re;		tmp1->d4 = im;
		tmp0->d5 =t0A;		tmp1->d5 =t09;

		tmp2->d4 =t01;		tmp3->d4 =t02;
		tmp2->d5 =t0E;		tmp3->d5 =t0D;

		tmp4->d4 =t03;		tmp5->d4 =t04;
		tmp4->d5 =t0C;		tmp5->d5 =t0B;

		tmp6->d4 =t05;		tmp7->d4 =t06;
		tmp6->d5 =t08;		tmp7->d5 =t07;

		re = dmult*c01->d6;	im = dmult*(c01+1)->d6;	// ...and d6.
		t01=(re-im)*ISRT2;	t02=(re+im)*ISRT2;

		re0=re*c;	im0=im*s;	re1=re*s;	im1=im*c;
		t03=re0-im0;	t04=re1+im1;
		t05=re1-im1;	t06=re0+im0;

		re0=re*c32_1;	im0=im*s32_1;	re1=re*s32_1;	im1=im*c32_1;
		t07=re0-im0;	t08=re1+im1;
		t09=re1-im1;	t0A=re0+im0;

		re0=re*c32_3;	im0=im*s32_3;	re1=re*s32_3;	im1=im*c32_3;
		t0B=re1-im1;	t0C=re0+im0;
		t0D=re0-im0;	t0E=re1+im1;

		tmp0->d6 = re;		tmp1->d6 = im;
		tmp0->d7 =t0A;		tmp1->d7 =t09;

		tmp2->d6 =t01;		tmp3->d6 =t02;
		tmp2->d7 =t0E;		tmp3->d7 =t0D;

		tmp4->d6 =t03;		tmp5->d6 =t04;
		tmp4->d7 =t0C;		tmp5->d7 =t0B;

		tmp6->d6 =t05;		tmp7->d6 =t06;
		tmp6->d7 =t08;		tmp7->d7 =t07;
	  #endif

	  if (fwd_fft_only) {	// 2-input modmul: fwdFFT data dyadic-mul'ed with precomputed 2nd-vector stored in fwdFFTed form at address stored in (uint64)-cast form in fwd_fft_only

	  #ifdef USE_AVX512		// Register blocks beginning with r00,08,10,18,20,28,30,38 map to memlocs bdd0-7, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd0+0x0; bpt1 = (vec_dbl*)bdd2+0x0; bpt2 = (vec_dbl*)bdd5+0x6; bpt3 = (vec_dbl*)bdd7+0x6;	// 00,10,2E,3E
	  #elif defined(USE_AVX)// Register blocks beginning with r00,10,20,30 map to memlocs bdd0-3, no offsets >= 16:
		bpt0 = (vec_dbl*)bdd0+0x0; bpt1 = (vec_dbl*)bdd1+0x0; bpt2 = (vec_dbl*)bdd2+0xe; bpt3 = (vec_dbl*)bdd3+0xe;	// 00,10,2E,3E
	  #else			// SSE2: Register blocks beginning with r00,r20 map to memlocs bdd0-1, no offsets >= 32:
		bpt0 = (vec_dbl*)bdd0+0x0; bpt1 = (vec_dbl*)bdd0+0x10; bpt2 = (vec_dbl*)bdd1+0xe; bpt3 = (vec_dbl*)bdd1+0x1e;	// 00,10,2E,3E
	  #endif
		PAIR_MUL_4_SSE2(r00,r10,r2E,r3E, bpt0,bpt1,bpt2,bpt3, tmp0,tmp1,forth);

	  #ifdef USE_AVX512		// Register blocks beginning with r00,08,10,18,20,28,30,38 map to memlocs bdd0-7, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd1+0x0; bpt1 = (vec_dbl*)bdd3+0x0; bpt2 = (vec_dbl*)bdd4+0x6; bpt3 = (vec_dbl*)bdd6+0x6;	// 08,18,26,36
	  #elif defined(USE_AVX)// Register blocks beginning with r00,10,20,30 map to memlocs bdd0-3, no offsets >= 16:
		bpt0 = (vec_dbl*)bdd0+0x8; bpt1 = (vec_dbl*)bdd1+0x8; bpt2 = (vec_dbl*)bdd2+0x6; bpt3 = (vec_dbl*)bdd3+0x6;	// 08,18,26,36
	  #else			// SSE2: Register blocks beginning with r00,r20 map to memlocs bdd0-1, no offsets >= 32:
		bpt0 = (vec_dbl*)bdd0+0x8; bpt1 = (vec_dbl*)bdd0+0x18; bpt2 = (vec_dbl*)bdd1+0x6; bpt3 = (vec_dbl*)bdd1+0x16;	// 08,18,26,36
	  #endif
		PAIR_MUL_4_SSE2(r08,r18,r26,r36, bpt0,bpt1,bpt2,bpt3, tmp2,tmp3,forth);

	  #ifdef USE_AVX512		// Register blocks beginning with r00,08,10,18,20,28,30,38 map to memlocs bdd0-7, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd0+0x4; bpt1 = (vec_dbl*)bdd2+0x4; bpt2 = (vec_dbl*)bdd5+0x2; bpt3 = (vec_dbl*)bdd7+0x2;	// 04,14,2A,3A
	  #elif defined(USE_AVX)// Register blocks beginning with r1,9,17,25 map to memlocs bdd0-3, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd0+0x4; bpt1 = (vec_dbl*)bdd1+0x4; bpt2 = (vec_dbl*)bdd2+0xa; bpt3 = (vec_dbl*)bdd3+0xa;	// 04,14,2A,3A
	  #else			// SSE2: Register blocks beginning with r1,17 map to memlocs bdd0-1, no offsets >= 16:
		bpt0 = (vec_dbl*)bdd0+0x4; bpt1 = (vec_dbl*)bdd0+0x14; bpt2 = (vec_dbl*)bdd1+0xa; bpt3 = (vec_dbl*)bdd1+0x1a;	// 04,14,2A,3A
	  #endif
		PAIR_MUL_4_SSE2(r04,r14,r2A,r3A, bpt0,bpt1,bpt2,bpt3, tmp4,tmp5,forth);

	  #ifdef USE_AVX512		// Register blocks beginning with r00,08,10,18,20,28,30,38 map to memlocs bdd0-7, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd1+0x4; bpt1 = (vec_dbl*)bdd3+0x4; bpt2 = (vec_dbl*)bdd4+0x2; bpt3 = (vec_dbl*)bdd6+0x2;	// 0C,1C,22,32
	  #elif defined(USE_AVX)// Register blocks beginning with r1,9,17,25 map to memlocs bdd0-3, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd0+0xc; bpt1 = (vec_dbl*)bdd1+0xc; bpt2 = (vec_dbl*)bdd2+0x2; bpt3 = (vec_dbl*)bdd3+0x2;	// 0C,1C,22,32
	  #else			// SSE2: Register blocks beginning with r1,17 map to memlocs bdd0-1, no offsets >= 16:
		bpt0 = (vec_dbl*)bdd0+0xc; bpt1 = (vec_dbl*)bdd0+0x1c; bpt2 = (vec_dbl*)bdd1+0x2; bpt3 = (vec_dbl*)bdd1+0x12;	// 0C,1C,22,32
	  #endif
		PAIR_MUL_4_SSE2(r0C,r1C,r22,r32, bpt0,bpt1,bpt2,bpt3, tmp6,tmp7,forth);

	  } else {	// fwd_fft_only = 0: Normal execution, dyadic-squaring followed by iFFT:

		/*
		!...send the pairs of complex elements which are to be combined and sincos temporaries needed for the squaring to a
		!   small subroutine. The j1 = 0 case is again exceptional. [For this reason we don't supply SSE2 code it - not worth the work].
		*/
	  #ifdef USE_AVX2
		PAIR_SQUARE_4_AVX2(r00,r10,r2E,r3E, tmp0,tmp1,
						   r08,r18,r26,r36, tmp2,tmp3,forth);
		PAIR_SQUARE_4_AVX2(r04,r14,r2A,r3A, tmp4,tmp5,
						   r0C,r1C,r22,r32, tmp6,tmp7,forth);
	  #else
		PAIR_SQUARE_4_SSE2(r00,r10,r2E,r3E, tmp0,tmp1,forth);
		PAIR_SQUARE_4_SSE2(r08,r18,r26,r36, tmp2,tmp3,forth);
		PAIR_SQUARE_4_SSE2(r04,r14,r2A,r3A, tmp4,tmp5,forth);
		PAIR_SQUARE_4_SSE2(r0C,r1C,r22,r32, tmp6,tmp7,forth);
	  #endif

	  }	// endif(fwd_fft_only != 0 or == 0)

	  // Permute the sincos-multiplier data:
	  #ifdef USE_ARM_V8_SIMD

		__asm__ volatile (\
			"ldr	x0,%[tmp0]			\n\t"\
			"ldp	q0,q1,[x0      ]	\n\t"\
			"ldp	q2,q3,[x0,#0x20]	\n\t"\
			"ldp	q4,q5,[x0,#0x40]	\n\t"\
			"ldp	q6,q7,[x0,#0x60]	\n\t"\
			"ext v0.16b,v0.16b,v0.16b,#8	\n\t"\
			"ext v1.16b,v1.16b,v1.16b,#8	\n\t"\
			"ext v2.16b,v2.16b,v2.16b,#8	\n\t"\
			"ext v3.16b,v3.16b,v3.16b,#8	\n\t"\
			"ext v4.16b,v4.16b,v4.16b,#8	\n\t"\
			"ext v5.16b,v5.16b,v5.16b,#8	\n\t"\
			"ext v6.16b,v6.16b,v6.16b,#8	\n\t"\
			"ext v7.16b,v7.16b,v7.16b,#8	\n\t"\
			"stp	q0,q1,[x0      ]	\n\t"\
			"stp	q2,q3,[x0,#0x20]	\n\t"\
			"stp	q4,q5,[x0,#0x40]	\n\t"\
			"stp	q6,q7,[x0,#0x60]	\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			: "cc","memory","x0","x1","v0","v1","v2","v3","v4","v5","v6","v7"	/* Clobbered registers */\
		);

	  #elif defined(USE_IMCI512)
		/* 1st-gen Xeon Phi (k1om) has no VSHUFPD:
		VSHUFPD with imm8 = 0x55 = 01010101_2 and src2 = src1 = dest simply swaps even/odd doubles of each each 128-bit
		[lo,hi]-double-pair. On k1om, effect this via VMOVAPD with register-to-self-copy using swap-inner-pair {cdab} swizzle:
		*/
		__asm__ volatile (\
			"movq	%[tmp0],%%rax\n\t"\
			"vmovaps	     (%%rax),%%zmm0\n\t"\
			"vmovaps	0x040(%%rax),%%zmm1\n\t"\
			"vmovaps	0x080(%%rax),%%zmm2\n\t"\
			"vmovaps	0x0c0(%%rax),%%zmm3\n\t"\
			"vmovaps	0x100(%%rax),%%zmm4\n\t"\
			"vmovaps	0x140(%%rax),%%zmm5\n\t"\
			"vmovaps	0x180(%%rax),%%zmm6\n\t"\
			"vmovaps	0x1c0(%%rax),%%zmm7\n\t"\
			"vmovapd	%%zmm0%{cdab%},%%zmm0\n\t"\
			"vmovapd	%%zmm1%{cdab%},%%zmm1\n\t"\
			"vmovapd	%%zmm2%{cdab%},%%zmm2\n\t"\
			"vmovapd	%%zmm3%{cdab%},%%zmm3\n\t"\
			"vmovapd	%%zmm4%{cdab%},%%zmm4\n\t"\
			"vmovapd	%%zmm5%{cdab%},%%zmm5\n\t"\
			"vmovapd	%%zmm6%{cdab%},%%zmm6\n\t"\
			"vmovapd	%%zmm7%{cdab%},%%zmm7\n\t"\
			"vmovaps	%%zmm0,     (%%rax)\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)\n\t"\
			"vmovaps	%%zmm2,0x080(%%rax)\n\t"\
			"vmovaps	%%zmm3,0x0c0(%%rax)\n\t"\
			"vmovaps	%%zmm4,0x100(%%rax)\n\t"\
			"vmovaps	%%zmm5,0x140(%%rax)\n\t"\
			"vmovaps	%%zmm6,0x180(%%rax)\n\t"\
			"vmovaps	%%zmm7,0x1c0(%%rax)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);

	  #elif defined(USE_AVX512)
		// AVX-512 version has shufpd immediate = 0x55 = 01010101_2, which is the fourfold analog of the SSE2 imm8 = 1 = 01_2:
		__asm__ volatile (\
			"movq	%[tmp0],%%rax\n\t"\
			"vmovaps	     (%%rax),%%zmm0\n\t"\
			"vmovaps	0x040(%%rax),%%zmm1\n\t"\
			"vmovaps	0x080(%%rax),%%zmm2\n\t"\
			"vmovaps	0x0c0(%%rax),%%zmm3\n\t"\
			"vmovaps	0x100(%%rax),%%zmm4\n\t"\
			"vmovaps	0x140(%%rax),%%zmm5\n\t"\
			"vmovaps	0x180(%%rax),%%zmm6\n\t"\
			"vmovaps	0x1c0(%%rax),%%zmm7\n\t"\
			"vshufpd	$0x55,%%zmm0,%%zmm0,%%zmm0\n\t"\
			"vshufpd	$0x55,%%zmm1,%%zmm1,%%zmm1\n\t"\
			"vshufpd	$0x55,%%zmm2,%%zmm2,%%zmm2\n\t"\
			"vshufpd	$0x55,%%zmm3,%%zmm3,%%zmm3\n\t"\
			"vshufpd	$0x55,%%zmm4,%%zmm4,%%zmm4\n\t"\
			"vshufpd	$0x55,%%zmm5,%%zmm5,%%zmm5\n\t"\
			"vshufpd	$0x55,%%zmm6,%%zmm6,%%zmm6\n\t"\
			"vshufpd	$0x55,%%zmm7,%%zmm7,%%zmm7\n\t"\
			"vmovaps	%%zmm0,     (%%rax)\n\t"\
			"vmovaps	%%zmm1,0x040(%%rax)\n\t"\
			"vmovaps	%%zmm2,0x080(%%rax)\n\t"\
			"vmovaps	%%zmm3,0x0c0(%%rax)\n\t"\
			"vmovaps	%%zmm4,0x100(%%rax)\n\t"\
			"vmovaps	%%zmm5,0x140(%%rax)\n\t"\
			"vmovaps	%%zmm6,0x180(%%rax)\n\t"\
			"vmovaps	%%zmm7,0x1c0(%%rax)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);

	  #elif defined(USE_AVX)	// 64-bit AVX build

		// AVX version has shufpd immediate = 5 = 0101_2, which is the doubled analog of the SSE2 imm8 = 1 = 01_2:
		__asm__ volatile (\
			"movq	%[tmp0],%%rax\n\t"\
			"vmovaps	    (%%rax),%%ymm0\n\t"\
			"vmovaps	0x20(%%rax),%%ymm1\n\t"\
			"vmovaps	0x40(%%rax),%%ymm2\n\t"\
			"vmovaps	0x60(%%rax),%%ymm3\n\t"\
			"vmovaps	0x80(%%rax),%%ymm4\n\t"\
			"vmovaps	0xa0(%%rax),%%ymm5\n\t"\
			"vmovaps	0xc0(%%rax),%%ymm6\n\t"\
			"vmovaps	0xe0(%%rax),%%ymm7\n\t"\
			"vshufpd	$5,%%ymm0,%%ymm0,%%ymm0\n\t"\
			"vshufpd	$5,%%ymm1,%%ymm1,%%ymm1\n\t"\
			"vshufpd	$5,%%ymm2,%%ymm2,%%ymm2\n\t"\
			"vshufpd	$5,%%ymm3,%%ymm3,%%ymm3\n\t"\
			"vshufpd	$5,%%ymm4,%%ymm4,%%ymm4\n\t"\
			"vshufpd	$5,%%ymm5,%%ymm5,%%ymm5\n\t"\
			"vshufpd	$5,%%ymm6,%%ymm6,%%ymm6\n\t"\
			"vshufpd	$5,%%ymm7,%%ymm7,%%ymm7\n\t"\
			"vmovaps	%%ymm0,    (%%rax)\n\t"\
			"vmovaps	%%ymm1,0x20(%%rax)\n\t"\
			"vmovaps	%%ymm2,0x40(%%rax)\n\t"\
			"vmovaps	%%ymm3,0x60(%%rax)\n\t"\
			"vmovaps	%%ymm4,0x80(%%rax)\n\t"\
			"vmovaps	%%ymm5,0xa0(%%rax)\n\t"\
			"vmovaps	%%ymm6,0xc0(%%rax)\n\t"\
			"vmovaps	%%ymm7,0xe0(%%rax)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);

	  #elif OS_BITS == 64		// 64-bit SSE2 build

		__asm__ volatile (\
			"movq	%[tmp0],%%rax\n\t"\
			"movaps	    (%%rax),%%xmm0\n\t"\
			"movaps	0x10(%%rax),%%xmm1\n\t"\
			"movaps	0x20(%%rax),%%xmm2\n\t"\
			"movaps	0x30(%%rax),%%xmm3\n\t"\
			"movaps	0x40(%%rax),%%xmm4\n\t"\
			"movaps	0x50(%%rax),%%xmm5\n\t"\
			"movaps	0x60(%%rax),%%xmm6\n\t"\
			"movaps	0x70(%%rax),%%xmm7\n\t"\
			"shufpd	$1	,%%xmm0	,%%xmm0\n\t"\
			"shufpd	$1	,%%xmm1	,%%xmm1\n\t"\
			"shufpd	$1	,%%xmm2	,%%xmm2\n\t"\
			"shufpd	$1	,%%xmm3	,%%xmm3\n\t"\
			"shufpd	$1	,%%xmm4	,%%xmm4\n\t"\
			"shufpd	$1	,%%xmm5	,%%xmm5\n\t"\
			"shufpd	$1	,%%xmm6	,%%xmm6\n\t"\
			"shufpd	$1	,%%xmm7	,%%xmm7\n\t"\
			"movaps	%%xmm0,    (%%rax)\n\t"\
			"movaps	%%xmm1,0x10(%%rax)\n\t"\
			"movaps	%%xmm2,0x20(%%rax)\n\t"\
			"movaps	%%xmm3,0x30(%%rax)\n\t"\
			"movaps	%%xmm4,0x40(%%rax)\n\t"\
			"movaps	%%xmm5,0x50(%%rax)\n\t"\
			"movaps	%%xmm6,0x60(%%rax)\n\t"\
			"movaps	%%xmm7,0x70(%%rax)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"				// Clobbered registers
		);

	  #else					// 32-bit SSE2 build
		#error 32-bit SSE2 build no longer supported as of Mlucas v19!
	  #endif

	  if (fwd_fft_only) {	// 2-input modmul: fwdFFT data dyadic-mul'ed with precomputed 2nd-vector stored in fwdFFTed form at address stored in (uint64)-cast form in fwd_fft_only

	  #ifdef USE_AVX512		// Register blocks beginning with r00,08,10,18,20,28,30,38 map to memlocs bdd0-7, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd0+0x2; bpt1 = (vec_dbl*)bdd2+0x2; bpt2 = (vec_dbl*)bdd5+0x4; bpt3 = (vec_dbl*)bdd7+0x4;	// 02,12,2C,3C
	  #elif defined(USE_AVX)// Register blocks beginning with r00,10,20,30 map to memlocs bdd0-3, no offsets >= 16:
		bpt0 = (vec_dbl*)bdd0+0x2; bpt1 = (vec_dbl*)bdd1+0x2; bpt2 = (vec_dbl*)bdd2+0xc; bpt3 = (vec_dbl*)bdd3+0xc;	// 02,12,2C,3C
	  #else			// SSE2: Register blocks beginning with r00,r20 map to memlocs bdd0-1, no offsets >= 32:
		bpt0 = (vec_dbl*)bdd0+0x2; bpt1 = (vec_dbl*)bdd0+0x12; bpt2 = (vec_dbl*)bdd1+0xc; bpt3 = (vec_dbl*)bdd1+0x1c;	// 02,12,2C,3C
	  #endif
		PAIR_MUL_4_SSE2(r02,r12,r2C,r3C, bpt0,bpt1,bpt2,bpt3, tmp7,tmp6,forth);

	  #ifdef USE_AVX512		// Register blocks beginning with r00,08,10,18,20,28,30,38 map to memlocs bdd0-7, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd1+0x2; bpt1 = (vec_dbl*)bdd3+0x2; bpt2 = (vec_dbl*)bdd4+0x4; bpt3 = (vec_dbl*)bdd6+0x4;	// 0A,1A,24,34
	  #elif defined(USE_AVX)// Register blocks beginning with r00,10,20,30 map to memlocs bdd0-3, no offsets >= 16:
		bpt0 = (vec_dbl*)bdd0+0xa; bpt1 = (vec_dbl*)bdd1+0xa; bpt2 = (vec_dbl*)bdd2+0x4; bpt3 = (vec_dbl*)bdd3+0x4;	// 0A,1A,24,34
	  #else			// SSE2: Register blocks beginning with r00,r20 map to memlocs bdd0-1, no offsets >= 32:
		bpt0 = (vec_dbl*)bdd0+0xa; bpt1 = (vec_dbl*)bdd0+0x1a; bpt2 = (vec_dbl*)bdd1+0x4; bpt3 = (vec_dbl*)bdd1+0x14;	// 0A,1A,24,34
	  #endif
		PAIR_MUL_4_SSE2(r0A,r1A,r24,r34, bpt0,bpt1,bpt2,bpt3, tmp5,tmp4,forth);

	  #ifdef USE_AVX512		// Register blocks beginning with r00,08,10,18,20,28,30,38 map to memlocs bdd0-7, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd0+0x6; bpt1 = (vec_dbl*)bdd2+0x6; bpt2 = (vec_dbl*)bdd5+0x0; bpt3 = (vec_dbl*)bdd7+0x0;	// 06,16,28,38
	  #elif defined(USE_AVX)// Register blocks beginning with r00,10,20,30 map to memlocs bdd0-3, no offsets >= 16:
		bpt0 = (vec_dbl*)bdd0+0x6; bpt1 = (vec_dbl*)bdd1+0x6; bpt2 = (vec_dbl*)bdd2+0x8; bpt3 = (vec_dbl*)bdd3+0x8;	// 06,16,28,38
	  #else			// SSE2: Register blocks beginning with r00,r20 map to memlocs bdd0-1, no offsets >= 32:
		bpt0 = (vec_dbl*)bdd0+0x6; bpt1 = (vec_dbl*)bdd0+0x16; bpt2 = (vec_dbl*)bdd1+0x8; bpt3 = (vec_dbl*)bdd1+0x18;	// 06,16,28,38
	  #endif
		PAIR_MUL_4_SSE2(r06,r16,r28,r38, bpt0,bpt1,bpt2,bpt3, tmp3,tmp2,forth);

	  #ifdef USE_AVX512		// Register blocks beginning with r00,08,10,18,20,28,30,38 map to memlocs bdd0-7, no offsets >= 8:
		bpt0 = (vec_dbl*)bdd1+0x6; bpt1 = (vec_dbl*)bdd3+0x6; bpt2 = (vec_dbl*)bdd4+0x0; bpt3 = (vec_dbl*)bdd6+0x0;	// 0E,1E,20,30
	  #elif defined(USE_AVX)// Register blocks beginning with r00,10,20,30 map to memlocs bdd0-3, no offsets >= 16:
		bpt0 = (vec_dbl*)bdd0+0xe; bpt1 = (vec_dbl*)bdd1+0xe; bpt2 = (vec_dbl*)bdd2+0x0; bpt3 = (vec_dbl*)bdd3+0x0;	// 0E,1E,20,30
	  #else			// SSE2: Register blocks beginning with r00,r20 map to memlocs bdd0-1, no offsets >= 32:
		bpt0 = (vec_dbl*)bdd0+0xe; bpt1 = (vec_dbl*)bdd0+0x1e; bpt2 = (vec_dbl*)bdd1+0x0; bpt3 = (vec_dbl*)bdd1+0x10;	// 0E,1E,20,30
	  #endif
		PAIR_MUL_4_SSE2(r0E,r1E,r20,r30, bpt0,bpt1,bpt2,bpt3, tmp1,tmp0,forth);

	  } else {	// fwd_fft_only = 0: Normal execution, dyadic-squaring followed by iFFT:

	  #ifdef USE_AVX2
		PAIR_SQUARE_4_AVX2(r02,r12,r2C,r3C, tmp7,tmp6,
						   r0A,r1A,r24,r34, tmp5,tmp4,forth);
		PAIR_SQUARE_4_AVX2(r06,r16,r28,r38, tmp3,tmp2,
						   r0E,r1E,r20,r30, tmp1,tmp0,forth);
	  #else
		PAIR_SQUARE_4_SSE2(r02,r12,r2C,r3C, tmp7,tmp6,forth);
		PAIR_SQUARE_4_SSE2(r0A,r1A,r24,r34, tmp5,tmp4,forth);
		PAIR_SQUARE_4_SSE2(r06,r16,r28,r38, tmp3,tmp2,forth);
		PAIR_SQUARE_4_SSE2(r0E,r1E,r20,r30, tmp1,tmp0,forth);
	  #endif

	  }	// endif(fwd_fft_only != 0 or == 0)

	}	// endif(fwd_fft_only == 1)

	/************************************************************************/
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/
	/************************************************************************/

	#ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX mode
	  					// because (add[1,3,5,7]-add[0,2,4,6]) have opposite signs for Fermat and Mersenne-mod:

		SSE2_RADIX32_WRAPPER_DIT(add0,add1,add2,add3,add4,add5,add6,add7
						,isrt2,r00,r08,r10,r20,r28,r30,c01,c02,c04,c06,c0A,c10,c12,c1A)

	#elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode

		SSE2_RADIX32_WRAPPER_DIT(add0,add1,add2,add3
						,isrt2,r00,r08,r10,r20,r28,r30,c01,c02,c04,c06,c08,c0A,c0C,c0E,c10,c12,c14,c16,c18,c1A,c1C,c1E)

	#else	// SSE2:

		SSE2_RADIX32_WRAPPER_DIT(add0,add1
		,isrt2,r00,r08,r10                ,c01,c02    ,c04    ,c06    ,c08,c0A,c0C,c0E,c10,c12,c14,c16,c18,c1A,c1C,c1E)

	#endif
	}	// endif(j1 == 0)

#endif	/* USE_SSE2 */

/*...Update the data (j1 and j2) array indices. */
loop:
	#ifdef USE_AVX
	  if(j1 <= 320) {
	#else
	  if(j1 <= 128) {
	#endif
		// Use scalar code (with index offsets properly fiddled) for j1 == 0 case in SIMD mode
		j1 = j1+64;
		j2 = j2-64;
		// Allow more than the default number of data blocks to be handled using the scalar-double DFT code -
		// this conditional will only get taken in the AVX-512 case with its (j1 <= 160) - modified conditionals:
	  #ifdef USE_AVX
		if(j1 < j2) {
		//	printf("j1,j2 = %d,%d; (j1 < j2), skip loop-control and goto jump_new.\n",j1,j2);
			j1pad = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* floating padded-array 1st element index is here */
			j2pad = j2 + ( (j2 >> DAT_BITS) << PAD_BITS );	/* floating padded-array 2nd element index is here */
			goto jump_new;
		}
	  #endif
	  } else {
		// Scalar and SSE2 mode both use same increment of 32; avx uses 64:
		j1 = j1+stride;
		j2 = j2-stride;
	  }
	}	/* endfor(m-loop) */

/*
!...Since the foregoing loop only gets executed half as many times as in the simple version, to properly position
!   ourselves in the data array for the start of the next block, need to bump up j1 by as much as would occur in a
!   second execution of the above loop. The exception is the first loop execution, where j1 needs to be doubled (64 x 2).
*/
update_blocklen:

	j1 = j1+(blocklen << 1);
	if(j2_start == n-64)
	{
		j1 = 0;
		return;
	}

/*...Reset half-complex-blocklength for next pass. If K >> 1 has a zero trailing bit, we multiply the blocklength by K >> 1 in preparation for the final block.	*/

	blocklen_sum = blocklen_sum + blocklen;
	if(i > 0) blocklen = (radix_prim[i-1]-1)*blocklen_sum;	// i = 0 update-value of blobklen not used, but wrap in if(i > 0) to get rid of UMR

/*...Next j2_start is previous one plus the (real) length of the current block = 4*(half-complex-blocklength) */

	j2_start = j2_start+(blocklen<<2);
	j2=j2_start;				/* Reset j2 for start of the next block. */

}	/* End of Main loop */

}

