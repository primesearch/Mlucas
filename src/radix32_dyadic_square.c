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

#ifdef USE_SSE2

	#include "sse2_macro_gcc64.h"
	#include "radix32_utils_asm.h"
	#include "radix32_dyadic_square_gcc64.h"

	#ifndef COMPILER_TYPE_GCC
		#error X86 SIMD build requires GCC-compatible compiler!
	#endif

	#include "radix32_wrapper_square_gcc64.h"

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
void radix32_dyadic_square(
	double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[],
	int ii, int nradices_prim, int radix_prim[], int incr, int init_sse2, int thr_id, uint64 fwd_fft_only, double c_arr[]
)
{
	const char func[] = "radix32_dyadic_square";
	const int stride = (int)RE_IM_STRIDE << 6, stridh = (stride>>1);	// main-array loop stride = 32*RE_IM_STRIDE
	static int max_threads = 0;
	static int nsave = 0;
	static int rad0save = 0, ndivrad0 = 0, ndivrad0m1 = 0;
/*	static int *index = 0x0;	OBSOLETE: full-N2/16-length Bit-reversal index array. */
	static int *index0 = 0x0, *index1 = 0x0, *index_ptmp0 = 0x0, *index_ptmp1 = 0x0;
	// Nov 2017: Making these static to support synthesized final-pass radices is not an option as it breaks multiheading,
	// so instead alloc an extra slot in sm_arr and stash each thread's entry/exit values of these there:
#ifdef MULTITHREAD
//	int *idx0_ptr, *idx1_ptr;	// Nov 2017: Experimental code, no gain, disable
#endif
	       int index0_idx= 0, index1_idx= 0;
	static int index0_mod=-1, index1_mod=-1;
	int nradices_prim_radix0;
	int i,j,j1,j2,l,iroot,k1,k2,nbytes;
	const double c = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/

	double re0,im0,re1,im1,rt,it;

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error X86 SIMD code not supported for this compiler!
  #endif
	static uint32 *sm_arr = 0x0,*sm_ptr;	// Base-ptr to arrays of k1,k2-index-vectors used in SIMD roots-computation.
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0,*add1, *bdd0;	/* Addresses into array sections */
	const double *cdd0;
  #ifdef USE_AVX
	double *add2,*add3;
  #endif
  #ifdef USE_AVX512
	double *add4,*add5,*add6,*add7;
  #endif
	vec_dbl *c_tmp,*s_tmp;
	vec_dbl *tmp,*tm1;

  #ifdef MULTITHREAD
	// Base addresses for discrete per-thread local stores ... 'r' for double-float data, 'i' for int:
	static vec_dbl *__r0;
	static uint32  *__i0;
	// In || mode, only above base-pointers (shared by all threads) are static:
	uint32  *k1_arr, *k2_arr;
	vec_dbl *isrt2,*sqrt2, *one,*two, *cc0,*ss0,*cc1,*ss1,*cc3,*ss3
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*r00,*r08,*r10,*r18,*r20,*r28,*r30,*r38;
  #else
	static uint32  *k1_arr, *k2_arr;
	static vec_dbl *isrt2,*sqrt2, *one,*two, *cc0,*ss0,*cc1,*ss1,*cc3,*ss3
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0A,*r0B,*r0C,*r0D,*r0E,*r0F
		,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1A,*r1B,*r1C,*r1D,*r1E,*r1F
		,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2A,*r2B,*r2C,*r2D,*r2E,*r2F
		,*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3A,*r3B,*r3C,*r3D,*r3E,*r3F;
  #endif

#else

  #if PFETCH
	double *add0;
  #endif
	double   c01,c02,c03,c04,c05,c06,c07,c08,c09,c0A,c0B,c0C,c0D,c0E,c0F
	,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c1A,c1B,c1C,c1D,c1E,c1F
			,s01,s02,s03,s04,s05,s06,s07,s08,s09,s0A,s0B,s0C,s0D,s0E,s0F
	,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s1A,s1B,s1C,s1D,s1E,s1F
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
	,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
	,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
	,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p0Ar,a1p0Br,a1p0Cr,a1p0Dr,a1p0Er,a1p0Fr
	,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p1Ar,a1p1Br,a1p1Cr,a1p1Dr,a1p1Er,a1p1Fr
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p0Ai,a1p0Bi,a1p0Ci,a1p0Di,a1p0Ei,a1p0Fi
	,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p1Ai,a1p1Bi,a1p1Ci,a1p1Di,a1p1Ei,a1p1Fi
	,b1p00r,b1p01r,b1p02r,b1p03r,b1p04r,b1p05r,b1p06r,b1p07r,b1p08r,b1p09r,b1p0Ar,b1p0Br,b1p0Cr,b1p0Dr,b1p0Er,b1p0Fr
	,b1p10r,b1p11r,b1p12r,b1p13r,b1p14r,b1p15r,b1p16r,b1p17r,b1p18r,b1p19r,b1p1Ar,b1p1Br,b1p1Cr,b1p1Dr,b1p1Er,b1p1Fr
	,b1p00i,b1p01i,b1p02i,b1p03i,b1p04i,b1p05i,b1p06i,b1p07i,b1p08i,b1p09i,b1p0Ai,b1p0Bi,b1p0Ci,b1p0Di,b1p0Ei,b1p0Fi
	,b1p10i,b1p11i,b1p12i,b1p13i,b1p14i,b1p15i,b1p16i,b1p17i,b1p18i,b1p19i,b1p1Ai,b1p1Bi,b1p1Ci,b1p1Di,b1p1Ei,b1p1Fi;

#endif

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
		ndivrad0 = n/radix0;	ndivrad0m1 = ndivrad0-1;	// ndivrad0 always a power of 2, so can do a fast-mod via & (ndivrad0-1)
		for(j = 0; j < ndivrad0; j += stride)
		{
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
			if( (j1+stridh) != (j+stridh) + ( ((j+stridh) >> DAT_BITS) << PAD_BITS ) ) {
				printf("j, j1, stride/2 = %d,%d,%d, jpad = %d\n",j,j1, stridh, (j+stridh) + (((j+stridh) >> DAT_BITS) << PAD_BITS) );
				ASSERT(0 , "add1 calculation violates padded index rules!");
			}
		}
		// Nov 2017: For the non-synthetic final-pass radices (16 and 32) the default contiguous-data chunksize
		// is n/radix0 doubles. For synthesized radices we need a chunksize which corresponds to treating the
		// product of all the upstream (non-synthesized) radices as 'radix0':
		if(RADIX_VEC[NRADICES-1] > 32) {
			for(j = 1; j < NRADICES-1; j++) {
				ndivrad0 /= RADIX_VEC[j];
			}
		}
		if(index_ptmp0) {
			free((void *)index_ptmp0);	index_ptmp0=0x0;	index0=0x0;
			free((void *)index_ptmp1);	index_ptmp1=0x0;	index1=0x0;
		}
		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Allocate and initialize an index array containing N/32 indices...

		index_ptmp = ALLOC_INT(N2/32);
		if(!index_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array INDEX in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		index = ALIGN_INT(index_ptmp);
		if(!index){ sprintf(cbuf,"ERROR: unable to allocate array ITMP in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		for(i=0; i < N2/32; i++)
		{
			index[i]=i;
		}
		*/
		index0_mod =        radix0;
		index1_mod = (n>>6)/radix0;	/* complex length requires an additional divide by 2 */

		index_ptmp0 = ALLOC_INT(index_ptmp0, index0_mod);
		if(!index_ptmp0){ sprintf(cbuf,"ERROR: unable to allocate array INDEX_PTMP0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		index0 = ALIGN_INT(index_ptmp0);

		index_ptmp1 = ALLOC_INT(index_ptmp1, index1_mod);
		if(!index_ptmp1){ sprintf(cbuf,"ERROR: unable to allocate array INDEX_PTMP1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		index1 = ALIGN_INT(index_ptmp1);

		for(i=0; i < index0_mod; i++){index0[i]=       i;}
		for(i=0; i < index1_mod; i++){index1[i]=radix0*i;}

		/*...then bit-reverse INDEX with respect to N/32. Radices are sent to BR routine
		in opposite order from which they are processed in forward FFT.
		bit_reverse_int(index, N2/32, nradices_prim-5, &radix_prim[nradices_prim-6], -1, (int *)arr_scratch);
		*/

		i = 1;
		for(nradices_prim_radix0 = 1; nradices_prim_radix0 < nradices_prim; nradices_prim_radix0++)
		{
			i *= radix_prim[nradices_prim_radix0-1];
			if(i == radix0)
				break;
		}
		if(nradices_prim_radix0 >= nradices_prim) { sprintf(cbuf,"ERROR: nradices_prim_radix0 must be < nradices_prim in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

		bit_reverse_int(index0, index0_mod,                 nradices_prim_radix0, &radix_prim[nradices_prim_radix0-1], -1, (int *)arr_scratch);
		bit_reverse_int(index1, index1_mod, nradices_prim-5-nradices_prim_radix0, &radix_prim[nradices_prim       -6], -1, (int *)arr_scratch);
	/* fprintf(stderr, "index0[dim = %d] =",index0_mod);	for(i=0; i < index0_mod; i++){fprintf(stderr, " %d",index0[i]);}	fprintf(stderr, "\n"); */
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
		if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)sm_arr);	sm_arr=0x0;
			free((void *)sc_arr);	sc_arr=0x0;
		}
		// Index vectors used in SIMD roots-computation.
		// Nov 2017: Add pair of int-slots per thread here ----vv, to support synthesized final-pass radices >= 256.
	//	sm_arr = ALLOC_INT(sm_arr, max_threads*(14*RE_IM_STRIDE+2) + 16);	if(!sm_arr){ sprintf(cbuf, "ERROR: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sm_arr = ALLOC_INT(sm_arr, max_threads* 14*RE_IM_STRIDE    + 16);	if(!sm_arr){ sprintf(cbuf, "ERROR: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sm_ptr = ALIGN_INT(sm_arr);
		ASSERT(((uintptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");
		// Twiddles-array:
		sc_arr = ALLOC_VEC_DBL(sc_arr, 0x94*max_threads + 100);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 64 vec_dbl slots of sc_arr for temporaries, next 8 for scratch, next 7 for the nontrivial complex 16th roots,
	next 62 for the doubled sincos twiddles, next 4 for [1.0,2.0,{0.25, unused in fermat-mod mode},sqrt2] and at least 3 more to allow for 64-byte alignment of the array.

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
		two   = sc_ptr + 0x8f;	// PAIR_SQUARE_4_AVX2 not used for Fermat-mod (i.e. no need for *forth), but note the SSE2_RADIX32_WRAPPER_DIF/T
	//	forth = sc_ptr + 0x90;	// asm macros are shared by Mersenne & Fermat, and assume two = isrt2 + 0x47, sqrt2 = two + 2, one = two - 1,
		sqrt2 = sc_ptr + 0x91;	// i.e. need to "leave slot for forth just above two".
		one	  = sc_ptr + 0x92;
		for(i = 0; i < max_threads; ++i) {
			/* These remain fixed within each per-thread local store: */
			VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
			VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );	//VEC_DBL_INIT(forth, 0.25 );
			VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
			VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
			VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
			/* Move on to next thread's local store */
			isrt2 += 0x94;
			cc0   += 0x94;
			ss0   += 0x94;
			cc1   += 0x94;
			ss1   += 0x94;
			cc3   += 0x94;
			ss3   += 0x94;
			one   += 0x94;
			two   += 0x94;
		//	forth += 0x94;
			sqrt2 += 0x94;
		}
	  #else
		k1_arr = sm_ptr;
		k2_arr = sm_ptr + 7*RE_IM_STRIDE;
		r00	 = sc_ptr + 0x00;	isrt2 = sc_ptr + 0x48;
		r02	 = sc_ptr + 0x02;	cc0   = sc_ptr + 0x49;	ss0 = sc_ptr + 0x4a;
		r04	 = sc_ptr + 0x04;	cc1   = sc_ptr + 0x4b;	ss1 = sc_ptr + 0x4c;
		r06	 = sc_ptr + 0x06;	cc3   = sc_ptr + 0x4d;	ss3 = sc_ptr + 0x4e;
		r08	 = sc_ptr + 0x08;	c00   = sc_ptr + 0x4f;
		r0A	 = sc_ptr + 0x0a;	c10   = sc_ptr + 0x51;
		r0C	 = sc_ptr + 0x0c;	c08   = sc_ptr + 0x53;
		r0E	 = sc_ptr + 0x0e;	c18   = sc_ptr + 0x55;
		r10	 = sc_ptr + 0x10;	c04   = sc_ptr + 0x57;
		r12	 = sc_ptr + 0x12;	c14   = sc_ptr + 0x59;
		r14	 = sc_ptr + 0x14;	c0C   = sc_ptr + 0x5b;
		r16	 = sc_ptr + 0x16;	c1C   = sc_ptr + 0x5d;
		r18	 = sc_ptr + 0x18;	c02   = sc_ptr + 0x5f;
		r1A	 = sc_ptr + 0x1a;	c12   = sc_ptr + 0x61;
		r1C	 = sc_ptr + 0x1c;	c0A   = sc_ptr + 0x63;
		r1E	 = sc_ptr + 0x1e;	c1A   = sc_ptr + 0x65;
		r20	 = sc_ptr + 0x20;	c06   = sc_ptr + 0x67;
		r22	 = sc_ptr + 0x22;	c16   = sc_ptr + 0x69;
		r24	 = sc_ptr + 0x24;	c0E   = sc_ptr + 0x6b;
		r26	 = sc_ptr + 0x26;	c1E   = sc_ptr + 0x6d;
		r28	 = sc_ptr + 0x28;	c01   = sc_ptr + 0x6f;
		r2A	 = sc_ptr + 0x2a;	c11   = sc_ptr + 0x71;
		r2C	 = sc_ptr + 0x2c;	c09   = sc_ptr + 0x73;
		r2E	 = sc_ptr + 0x2e;	c19   = sc_ptr + 0x75;
		r30	 = sc_ptr + 0x30;	c05   = sc_ptr + 0x77;
		r32	 = sc_ptr + 0x32;	c15   = sc_ptr + 0x79;
		r34	 = sc_ptr + 0x34;	c0D   = sc_ptr + 0x7b;
		r36	 = sc_ptr + 0x36;	c1D   = sc_ptr + 0x7d;
		r38	 = sc_ptr + 0x38;	c03   = sc_ptr + 0x7f;
		r3A	 = sc_ptr + 0x3a;	c13   = sc_ptr + 0x81;
		r3C	 = sc_ptr + 0x3c;	c0B   = sc_ptr + 0x83;
		r3E	 = sc_ptr + 0x3e;	c1B   = sc_ptr + 0x85;
								c07   = sc_ptr + 0x87;
								c17   = sc_ptr + 0x89;
								c0F   = sc_ptr + 0x8b;
								c1F   = sc_ptr + 0x8d;
								two   = sc_ptr + 0x8f;
							//	forth = sc_ptr + 0x90;
								sqrt2 = sc_ptr + 0x91;
								one	  = sc_ptr + 0x92;
		/* These remain fixed within each per-thread local store: */
		VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
		VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );	//VEC_DBL_INIT(forth, 0.25 );
		VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
		VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
		VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
	  #endif
	#endif	// USE_SSE2
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
#ifdef MULTITHREAD
	ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
  #ifdef USE_SSE2
	k1_arr =   __i0 + thr_id*14*RE_IM_STRIDE;
	k2_arr = k1_arr + 7*RE_IM_STRIDE;
	// Nov 2017: If synthesized final-pass radix scheme showed any gain would want to move these so as to
	// also support them in non-SIMD mode, but no speedup observed, so leave in SIMD-only code:
//	idx0_ptr = k2_arr + 7*RE_IM_STRIDE;	idx1_ptr = idx0_ptr+1;	index0_idx = *idx0_ptr;	index1_idx = *idx1_ptr;
	r00 = __r0 + thr_id*148;	isrt2 = r00 + 0x48;
/*	r02  = r00 + 0x02;*/cc0 = isrt2 + 0x01;
/*	r04  = r00 + 0x04;*/cc1   = cc0 + 0x02;
/*	r06  = r00 + 0x06;*/cc3   = cc0 + 0x04;
	r08  = r00 + 0x08;	c00   = cc0 + 0x06;
/*	r0A  = r00 + 0x0a;*/c10   = cc0 + 0x08;
/*	r0C  = r00 + 0x0c;*/c08   = cc0 + 0x0a;
/*	r0E  = r00 + 0x0e;*/c18   = cc0 + 0x0c;
	r10  = r00 + 0x10;	c04   = cc0 + 0x0e;
/*	r12  = r00 + 0x12;*/c14   = cc0 + 0x10;
/*	r14  = r00 + 0x14;*/c0C   = cc0 + 0x12;
/*	r16  = r00 + 0x16;*/c1C   = cc0 + 0x14;
	r18  = r00 + 0x18;	c02   = cc0 + 0x16;
/*	r1A  = r00 + 0x1a;*/c12   = cc0 + 0x18;
/*	r1C  = r00 + 0x1c;*/c0A   = cc0 + 0x1a;
/*	r1E  = r00 + 0x1e;*/c1A   = cc0 + 0x1c;
	r20  = r00 + 0x20;	c06   = cc0 + 0x1e;
/*	r22  = r00 + 0x22;*/c16   = cc0 + 0x20;
/*	r24  = r00 + 0x24;*/c0E   = cc0 + 0x22;
/*	r26  = r00 + 0x26;*/c1E   = cc0 + 0x24;
	r28  = r00 + 0x28;	c01   = cc0 + 0x26;
/*	r2A  = r00 + 0x2a;*/c11   = cc0 + 0x28;
/*	r2C  = r00 + 0x2c;*/c09   = cc0 + 0x2a;
/*	r2E  = r00 + 0x2e;*/c19   = cc0 + 0x2c;
	r30  = r00 + 0x30;	c05   = cc0 + 0x2e;
/*	r32  = r00 + 0x32;*/c15   = cc0 + 0x30;
/*	r34  = r00 + 0x34;*/c0D   = cc0 + 0x32;
/*	r36  = r00 + 0x36;*/c1D   = cc0 + 0x34;
	r38  = r00 + 0x38;	c03   = cc0 + 0x36;
/*	r3A  = r00 + 0x3a;*/c13   = cc0 + 0x38;
/*	r3C  = r00 + 0x3c;*/c0B   = cc0 + 0x3a;
/*	r3E  = r00 + 0x3e;*/c1B   = cc0 + 0x3c;
						c07   = cc0 + 0x3e;
						c17   = cc0 + 0x40;
						c0F   = cc0 + 0x42;
						c1F   = cc0 + 0x44;
						two   = cc0 + 0x46;
					//	forth = cc0 + 0x47;
						sqrt2 = cc0 + 0x48;
						one	  = cc0 + 0x49;
  #endif
#endif

	/*...If a new runlength, should not get to this point: */
	ASSERT(n == nsave,"n != nsave");
	ASSERT(incr == 64,"incr == 64");
//	ASSERT(ndivrad0 == n/radix0,"bad value for ndivrad0!");	Synthesized final-pass radices break this
	/*
	k = ii*(ndivrad0 >> 6);
	*/
	// Nov 2017: To support synthesized final-pass radices, only reset iroot = 0 if our fiddled value of ndivrad0 is a multiple of n/radix0,
	// or index0_idx > (radix0-1). The latter is to ensure proper reinit at beginning of each iteration.
	//***Update: No gain, and would need to make threadsafe, so revert to old non-static vars unconditional re-init of these here:
//	if(index0_idx >= radix0 || (ndivrad0 & ndivrad0m1) == 0) {	//if(thr_id < 2 || thr_id == (radix0-1))printf("%s: Thread %u, Reset index0,1_idx\n",func,thr_id);
		index0_idx = ii;
		index1_idx = 0;
//	}	//if(thr_id < 2 || thr_id == (radix0-1))printf("%s: Thread %u, ndivrad0 = %u, stride = %u, index0,1_idx = %u,%u\n",func,thr_id,ndivrad0,stride,index0_idx,index1_idx);
	for(j = 0; j < ndivrad0; j += stride)
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
		j2 = j1+RE_IM_STRIDE;

	#ifndef USE_SSE2	// Scalar-double mode:

		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
		l = iroot;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c01 =rt;	s01 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c02 =rt;	s02 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += 4*iroot;			/* 7*iroot	*/
		iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c03 =rt;	s03 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c07 =rt;	s07 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c0E =rt;	s0E =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c15 =rt;	s15 =it;

		k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c1C =rt;	s1C =it;

		/*(c,s)4-10 (decimal indices), 4-0A (hexadecimal indices):	*/
		t00=c01*c07; t01=c01*s07; t02=s01*c07; t03=s01*s07;
		c06=t00+t03; s06=t01-t02; c08=t00-t03; s08=t01+t02;

		t00=c02*c07; t01=c02*s07; t02=s02*c07; t03=s02*s07;
		c05=t00+t03; s05=t01-t02; c09=t00-t03; s09=t01+t02;

		t00=c03*c07; t01=c03*s07; t02=s03*c07; t03=s03*s07;
		c04=t00+t03; s04=t01-t02; c0A=t00-t03; s0A=t01+t02;

		/*(c,s)11-17 (decimal indices), 0B-11 (hexadecimal indices):;	*/
		t00=c01*c0E; t01=c01*s0E; t02=s01*c0E; t03=s01*s0E;
		c0D=t00+t03; s0D=t01-t02; c0F=t00-t03; s0F=t01+t02;

		t00=c02*c0E; t01=c02*s0E; t02=s02*c0E; t03=s02*s0E;
		c0C=t00+t03; s0C=t01-t02; c10=t00-t03; s10=t01+t02;

		t00=c03*c0E; t01=c03*s0E; t02=s03*c0E; t03=s03*s0E;
		c0B=t00+t03; s0B=t01-t02; c11=t00-t03; s11=t01+t02;

		/*(c,s)18-24 (decimal indices), 12-18 (hexadecimal indices):	*/
		t00=c01*c15; t01=c01*s15; t02=s01*c15; t03=s01*s15;
		c14=t00+t03; s14=t01-t02; c16=t00-t03; s16=t01+t02;

		t00=c02*c15; t01=c02*s15; t02=s02*c15; t03=s02*s15;
		c13=t00+t03; s13=t01-t02; c17=t00-t03; s17=t01+t02;

		t00=c03*c15; t01=c03*s15; t02=s03*c15; t03=s03*s15;
		c12=t00+t03; s12=t01-t02; c18=t00-t03; s18=t01+t02;

		/*(c,s)25-31 (decimal indices), 19-1F (hexadecimal indices):	*/
		t00=c01*c1C; t01=c01*s1C; t02=s01*c1C; t03=s01*s1C;
		c1B=t00+t03; s1B=t01-t02; c1D=t00-t03; s1D=t01+t02;

		t00=c02*c1C; t01=c02*s1C; t02=s02*c1C; t03=s02*s1C;
		c1A=t00+t03; s1A=t01-t02; c1E=t00-t03; s1E=t01+t02;

		t00=c03*c1C; t01=c03*s1C; t02=s03*c1C; t03=s03*s1C;
		c19=t00+t03; s19=t01-t02; c1F=t00-t03; s1F=t01+t02;

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

		// In SSE2 mode need 2 sets of sincos - since we may use different indexing patterns for storing k1,2_arr data
		// in AVX+ mode, don't share code for first 2 sets with those wider SIMD cases:
	  #if defined(USE_SSE2) && !defined(USE_AVX)

		// 1st set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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

		// 3rd set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX32_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #elif defined(USE_AVX512)	// AVX512: need 8 sets of sincos gathered into 512-bit vectors-of-doubles

	   #ifdef USE_IMCI512	// Vectorized version below a no-go on 1st-gen Xeon Phi

		#warning AVX-512: Using slower non-ASM-macro version of sincos-indexing in radix32_dyadic_square.c.
		/*
		IMCI512: SSE2_RADIX16_CALC_TWIDDLES_LOACC reads k[0|1]_arr[0-55] in 8-index chunks,
			does gather-loads associated rt[0|1] elts, computes twiddles, writes to cc0+[0x8-...].
			Pull the gather-loads out here and do in C, just do the 512-bit vector CMULs in asm-macro:
		*/
		// 1st set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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

		// 3rd set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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

		// 5th set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }
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

	   #else

		#warning Using AVX-512 code in radix32_dyadic_square.c
		/*
		Cf. the radix-16 version of this routine for notes re. the following vector-index-computation code:
		*/
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
		- separately compute 3*ymm to obtain k1,2_arr[16-23];
		- multiply zmm by 7 to get [7,14]*iroot, mask/shift yields k1,2_arr[24-39];
		- put 3x8, 2x8 in a pair of ymm-regs, vinsert-merge, *= zmm to get [21,28]*iroot, mask/shift yields k1,2_arr[40-55].
		*/\
			"vpaddd	%%ymm0,%%ymm0,%%ymm4	\n\t"/* 2*iroot... */\
			"vpaddd	%%ymm0,%%ymm4,%%ymm5	\n\t"/* ymm5 = 3*iroot */\
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
		"movl $7,%%edi \n\t vpbroadcastd %%edi,%%zmm6	\n\t"/* (uint32)7 x 16 */\
		"vpmulld %%zmm6,%%zmm4,%%zmm4	\n\t"/* zmm4 = [7,14]*iroot; use full 16-element multiply, obviously */\
			"vpslld	$4,%%zmm0,%%zmm0		\n\t	vpslld	$4,%%zmm1,%%zmm1	\n\t"/* k1,2<<4 */\
			"vmovaps	%%zmm0,0x00(%%rax)	\n\t	vmovups	%%zmm1,0x00(%%rbx)	\n\t"/* Write k1,2_arr[ 0-15] = k1,2<<4 */\
			/* 3*iroot, using half-register-width operands: */\
			"vpand	%%ymm2,%%ymm5,%%ymm0	\n\t"\
			"vpsrlvd %%ymm3,%%ymm5,%%ymm1	\n\t"\
		"movl $3,%%edi \n\t vpbroadcastd %%edi,%%zmm5	\n\t"/* (uint32)3 x 16, 2x8 will overwrite high half */\
		"movl $2,%%edi \n\t vpbroadcastd %%edi,%%zmm6	\n\t"/* (uint32)2 x 16, only use low half */\
		"vinserti64x4 $1,%%ymm6,%%zmm5,%%zmm5	\n\t"/* Concatenate the [3x8,2x8] uint32-octets into a single zmm */\
			"vpslld	$4,%%ymm0,%%ymm0		\n\t	vpslld	$4,%%ymm1,%%ymm1	\n\t"\
			"vmovaps	%%ymm0,0x40(%%rax)	\n\t	vmovaps	%%ymm1,0x40(%%rbx)	\n\t"/* Write k1,2_arr[16-23] */\
			/* [7,14]*iroot - note all writes from hereon must be unaligned, since offsets no longer multiples of 0x40: */\
			"vpandd	%%zmm2,%%zmm4,%%zmm0	\n\t"\
			"vpsrlvd %%zmm3,%%zmm4,%%zmm1	\n\t"\
		"vpmulld %%zmm5,%%zmm4,%%zmm4	\n\t"/* zmm4 = [21,28]*iroot */\
			"vpslld	$4,%%zmm0,%%zmm0		\n\t	vpslld	$4,%%zmm1,%%zmm1	\n\t"\
			"vmovups	%%zmm0,0x60(%%rax)	\n\t	vmovups	%%zmm1,0x60(%%rbx)	\n\t"/* Write k1,2_arr[24-39] */\
			/* [21,28]*iroot: */\
			"vpandd	%%zmm2,%%zmm4,%%zmm0	\n\t"\
			"vpsrlvd %%zmm3,%%zmm4,%%zmm1	\n\t"\
			"vpslld	$4,%%zmm0,%%zmm0		\n\t	vpslld	$4,%%zmm1,%%zmm1	\n\t"\
			"vmovups	%%zmm0,0xa0(%%rax)	\n\t	vmovups	%%zmm1,0xa0(%%rbx)	\n\t"/* Write k1,2_arr[40-55] */\
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
			: "cc","memory","cl","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6"	/* Clobbered registers */\
		);

		// Need one more (scalar-mode) update of these in preparation for the next 8-fold chunk:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod) { index1_idx -= index1_mod;	++index0_idx; }

	   #endif	// (IMCI512 or AVX512?) toggle

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX32_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);
/*
printf("AVX512: Outputs of twiddle-comp:\n");
for(l = 0, tmp = cc0; l < 35; l++, tmp += 2) {	// 3 [cc,ss]-terms plus 32 complex twiddles
printf("c[%2d] = %18.15f,%18.15f,%18.15f,%18.15f,%18.15f,%18.15f,%18.15f,%18.15f\ns[%2d] = %18.15f,%18.15f,%18.15f,%18.15f,%18.15f,%18.15f,%18.15f,%18.15f\n",l-3,tmp->d0,tmp->d1,tmp->d2,tmp->d3,tmp->d4,tmp->d5,tmp->d6,tmp->d7,l-3,(tmp+1)->d0,(tmp+1)->d1,(tmp+1)->d2,(tmp+1)->d3,(tmp+1)->d4,(tmp+1)->d5,(tmp+1)->d6,(tmp+1)->d7);
} exit(0);
*/
	  #endif	// SIMD mode?

	#endif	// SIMD ?

	/* gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...We process the sincos data in bit-reversed order. */

	// See the radix-16 analog of this routine for details about the SSE2 and AVX data layouts

	#if defined(USE_SSE2)	// SSE2 and AVX share same code flow, just with different versions of the ASM compute macros

		add0 = &a[j1];
		add1 = &a[j1+stridh];
	//	printf("stride = %d, add0,1 = %" PRIX64 ", %" PRIX64 ", diff = %" PRIX64 "\n",stride,(int64)add0,(int64)add1, (int64)add1-(int64)add0);	exit(0);
	//	ASSERT((j1+stride) == (j+stride) + ( ((j+stride) >> DAT_BITS) << PAD_BITS ) , "add1 calculation violates padded index rules!");
	  #ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX mode
	  					// because (add[1,3,5,7]-add[0,2,4,6]) have opposite signs for Fermat and Mersenne-mod:
		add1 = add0 +  64;
		add2 = add0 + 128;
		add3 = add0 + 192;
		add4 = add0 + 256;
		add5 = add0 + 320;
		add6 = add0 + 384;
		add7 = add0 + 448;
	  #elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  					// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:
		add1 = add0 +  64;
		add2 = add0 + 128;
		add3 = add0 + 192;
	#endif

	if(fwd_fft_only != 3)	// v20: add support for both-inputs-already-fwd-FFTed case - fwd_fft_only == 3 means skip fwd-FFT
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
		nbytes = ((intptr_t)r20 - (intptr_t)r00)<<1;
		memcpy(add0, r00, nbytes);	// add0 = a + j1pad;
		goto loop;	// Skip dyadic-mul and iFFT

	} else if (fwd_fft_only) {	// 2-input modmul: fwdFFT data dyadic-mul'ed with precomputed 2nd-vector stored in fwdFFTed form at address stored in (uint64)-cast form in fwd_fft_only

	  if(fwd_fft_only == 3) {	// v20: Both inputs enter fwd-FFTed, must copy from main arrays a[],b[] to local-vars
		nbytes = ((intptr_t)r20 - (intptr_t)r00)<<1;
		memcpy(r00, add0, nbytes);	// add0 = a + j1pad;
	  }
		bdd0 = b + j1;	// No reverse-running addresses as in Mers-mod/real-vector-FFT case, just need single B-array base address

		// Use high bit of c_arr as a flag - If set, it means "do a mulsub a*b - c" on the data at the complemented address:
		if((int64)c_arr < 0) {	// (int64)c_arr < 0x0: a*b - c

			cdd0 = (double*)(~(uint64)c_arr) + j1;
			SIMD_MULSUB(r00,bdd0,cdd0)

		} else if(c_arr) {	// c_arr != 0x0: a * (b - c)

			cdd0 = c_arr + j1;
			SIMD_SUBMUL(r00,bdd0,cdd0)

		} else {	// c_arr = 0x0: a * b

			SIMD_MUL(r00,bdd0)

		}	// endif(c_arr)

	} else {	// fwd_fft_only = 0: Normal execution, dyadic-squaring followed by iFFT:

		SIMD_SQR(r00)

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

#else	// USE_SSE2 = False, scalar-double mode:

	if(fwd_fft_only == 3)
		goto skip_fwd_fft;	// v20: jump-to-point for both-inputs-already-fwd-FFTed case

		/*...Block 1:	*/
		t00=a[j1   ];						t01=a[j2     ];
		rt =a[j1+32]*c10 - a[j2+32]*s10;	it =a[j2+32]*c10 + a[j1+32]*s10;
		t02=t00-rt;		t03=t01-it;
		t00=t00+rt;		t01=t01+it;

		t04=a[j1+16]*c08 - a[j2+16]*s08;	t05=a[j2+16]*c08 + a[j1+16]*s08;
		rt =a[j1+48]*c18 - a[j2+48]*s18;	it =a[j2+48]*c18 + a[j1+48]*s18;
		t06=t04-rt;		t07=t05-it;
		t04=t04+rt;		t05=t05+it;

		rt =t04;		it =t05;
		t04=t00-rt;		t05=t01-it;
		t00=t00+rt;		t01=t01+it;

		rt =t06;		it =t07;
		t06=t02+it;		t07=t03-rt;
		t02=t02-it;		t03=t03+rt;

		t08=a[j1+8 ]*c04 - a[j2+8 ]*s04;	t09=a[j2+8 ]*c04 + a[j1+8 ]*s04;
		rt =a[j1+40]*c14 - a[j2+40]*s14;	it =a[j2+40]*c14 + a[j1+40]*s14;
		t0A=t08-rt;		t0B=t09-it;
		t08=t08+rt;		t09=t09+it;

		t0C=a[j1+24]*c0C - a[j2+24]*s0C;	t0D=a[j2+24]*c0C + a[j1+24]*s0C;
		rt =a[j1+56]*c1C - a[j2+56]*s1C;	it =a[j2+56]*c1C + a[j1+56]*s1C;
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
		t10=a[j1+4 ]*c02 - a[j2+4 ]*s02;	t11=a[j2+4 ]*c02 + a[j1+4 ]*s02;
		rt =a[j1+36]*c12 - a[j2+36]*s12;	it =a[j2+36]*c12 + a[j1+36]*s12;
		t12=t10-rt;		t13=t11-it;
		t10=t10+rt;		t11=t11+it;

		t14=a[j1+20]*c0A - a[j2+20]*s0A;	t15=a[j2+20]*c0A + a[j1+20]*s0A;
		rt =a[j1+52]*c1A - a[j2+52]*s1A;	it =a[j2+52]*c1A + a[j1+52]*s1A;
		t16=t14-rt;		t17=t15-it;
		t14=t14+rt;		t15=t15+it;

		rt =t14;		it =t15;
		t14=t10-rt;		t15=t11-it;
		t10=t10+rt;		t11=t11+it;

		rt =t16;		it =t17;
		t16=t12+it;		t17=t13-rt;
		t12=t12-it;		t13=t13+rt;

		t18=a[j1+12]*c06 - a[j2+12]*s06;	t19=a[j2+12]*c06 + a[j1+12]*s06;
		rt =a[j1+44]*c16 - a[j2+44]*s16;	it =a[j2+44]*c16 + a[j1+44]*s16;
		t1A=t18-rt;		t1B=t19-it;
		t18=t18+rt;		t19=t19+it;

		t1C=a[j1+28]*c0E - a[j2+28]*s0E;	t1D=a[j2+28]*c0E + a[j1+28]*s0E;
		rt =a[j1+60]*c1E - a[j2+60]*s1E;	it =a[j2+60]*c1E + a[j1+60]*s1E;
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
		t20=a[j1+2 ]*c01 - a[j2+2 ]*s01;	t21=a[j2+2 ]*c01 + a[j1+2 ]*s01;
		rt =a[j1+34]*c11 - a[j2+34]*s11;	it =a[j2+34]*c11 + a[j1+34]*s11;
		t22=t20-rt;		t23=t21-it;
		t20=t20+rt;		t21=t21+it;

		t24=a[j1+18]*c09 - a[j2+18]*s09;	t25=a[j2+18]*c09 + a[j1+18]*s09;
		rt =a[j1+50]*c19 - a[j2+50]*s19;	it =a[j2+50]*c19 + a[j1+50]*s19;
		t26=t24-rt;		t27=t25-it;
		t24=t24+rt;		t25=t25+it;

		rt =t24;		it =t25;
		t24=t20-rt;		t25=t21-it;
		t20=t20+rt;		t21=t21+it;

		rt =t26;		it =t27;
		t26=t22+it;		t27=t23-rt;
		t22=t22-it;		t23=t23+rt;

		t28=a[j1+10]*c05 - a[j2+10]*s05;	t29=a[j2+10]*c05 + a[j1+10]*s05;
		rt =a[j1+42]*c15 - a[j2+42]*s15;	it =a[j2+42]*c15 + a[j1+42]*s15;
		t2A=t28-rt;		t2B=t29-it;
		t28=t28+rt;		t29=t29+it;

		t2C=a[j1+26]*c0D - a[j2+26]*s0D;	t2D=a[j2+26]*c0D + a[j1+26]*s0D;
		rt =a[j1+58]*c1D - a[j2+58]*s1D;	it =a[j2+58]*c1D + a[j1+58]*s1D;
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
		t30=a[j1+6 ]*c03 - a[j2+6 ]*s03;	t31=a[j2+6 ]*c03 + a[j1+6 ]*s03;
		rt =a[j1+38]*c13 - a[j2+38]*s13;	it =a[j2+38]*c13 + a[j1+38]*s13;
		t32=t30-rt;		t33=t31-it;
		t30=t30+rt;		t31=t31+it;

		t34=a[j1+22]*c0B - a[j2+22]*s0B;	t35=a[j2+22]*c0B + a[j1+22]*s0B;
		rt =a[j1+54]*c1B - a[j2+54]*s1B;	it =a[j2+54]*c1B + a[j1+54]*s1B;
		t36=t34-rt;		t37=t35-it;
		t34=t34+rt;		t35=t35+it;

		rt =t34;		it =t35;
		t34=t30-rt;		t35=t31-it;
		t30=t30+rt;		t31=t31+it;

		rt =t36;		it =t37;
		t36=t32+it;		t37=t33-rt;
		t32=t32-it;		t33=t33+rt;

		t38=a[j1+14]*c07 - a[j2+14]*s07;	t39=a[j2+14]*c07 + a[j1+14]*s07;
		rt =a[j1+46]*c17 - a[j2+46]*s17;	it =a[j2+46]*c17 + a[j1+46]*s17;
		t3A=t38-rt;		t3B=t39-it;
		t38=t38+rt;		t39=t39+it;

		t3C=a[j1+30]*c0F - a[j2+30]*s0F;	t3D=a[j2+30]*c0F + a[j1+30]*s0F;
		rt =a[j1+62]*c1F - a[j2+62]*s1F;	it =a[j2+62]*c1F + a[j1+62]*s1F;
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

		/*...and now do eight radix-4 transforms, including the internal twiddle factors:
		1, exp(i* 1*twopi/32) =       ( c32_1, s32_1), exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 3*twopi/32) =       ( c32_3, s32_3) (for inputs to transform block 2),
		1, exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 3*twopi/32) =       ( s    , c    ) (for inputs to transform block 3),
		1, exp(i* 3*twopi/32) =       ( c32_3, s32_3), exp(i* 6*twopi/32) =       ( s    , c    ), exp(i* 9*twopi/32) =       (-s32_1, c32_1) (for inputs to transform block 4),
		1, exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 8*twopi/32) =       ( 0    , 1    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ) (for inputs to transform block 5),
		1, exp(i* 5*twopi/32) =       ( s32_3, c32_3), exp(i*10*twopi/32) =       (-s    , c    ), exp(i*15*twopi/32) =       (-c32_1, s32_1) (for inputs to transform block 6),
		1, exp(i* 6*twopi/32) =       ( s    , c    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ), exp(i*18*twopi/32) =       (-c    ,-s    ) (for inputs to transform block 7),
		1, exp(i* 7*twopi/32) =       ( s32_1, c32_1), exp(i*14*twopi/32) =       (-c    , s    ), exp(i*21*twopi/32) =       (-s32_3,-c32_3) (for inputs to transform block 8),
		and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.
		Within each block we process the 8 needed twiddles in bit-reversed order.	*/

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

skip_fwd_fft:	// v20: jump-to-point for both-inputs-already-fwd-FFTed case

	// v20: If fwd_fft_only = 1, write fwd-FFT result back to input array, skipping dyadic-square and inv-FFT steps:
	if(fwd_fft_only == 1)
	{
	// Block 1:
		a[j1   ] = a1p00r;	a[j2   ] = a1p00i;
		a[j1+32] = a1p01r;	a[j2+32] = a1p01i;
		a[j1+16] = a1p02r;	a[j2+16] = a1p02i;
		a[j1+48] = a1p03r;	a[j2+48] = a1p03i;
		a[j1+ 8] = a1p04r;	a[j2+ 8] = a1p04i;
		a[j1+40] = a1p05r;	a[j2+40] = a1p05i;
		a[j1+24] = a1p06r;	a[j2+24] = a1p06i;
		a[j1+56] = a1p07r;	a[j2+56] = a1p07i;
	// Block 2:
		a[j1+ 4] = a1p08r;	a[j2+ 4] = a1p08i;
		a[j1+36] = a1p09r;	a[j2+36] = a1p09i;
		a[j1+20] = a1p0Ar;	a[j2+20] = a1p0Ai;
		a[j1+52] = a1p0Br;	a[j2+52] = a1p0Bi;
		a[j1+12] = a1p0Cr;	a[j2+12] = a1p0Ci;
		a[j1+44] = a1p0Dr;	a[j2+44] = a1p0Di;
		a[j1+28] = a1p0Er;	a[j2+28] = a1p0Ei;
		a[j1+60] = a1p0Fr;	a[j2+60] = a1p0Fi;
	// Block 3:
		a[j1+ 2] = a1p10r;	a[j2+ 2] = a1p10i;
		a[j1+34] = a1p11r;	a[j2+34] = a1p11i;
		a[j1+18] = a1p12r;	a[j2+18] = a1p12i;
		a[j1+50] = a1p13r;	a[j2+50] = a1p13i;
		a[j1+10] = a1p14r;	a[j2+10] = a1p14i;
		a[j1+42] = a1p15r;	a[j2+42] = a1p15i;
		a[j1+26] = a1p16r;	a[j2+26] = a1p16i;
		a[j1+58] = a1p17r;	a[j2+58] = a1p17i;
	// Block 4:
		a[j1+ 6] = a1p18r;	a[j2+ 6] = a1p18i;
		a[j1+38] = a1p19r;	a[j2+38] = a1p19i;
		a[j1+22] = a1p1Ar;	a[j2+22] = a1p1Ai;
		a[j1+54] = a1p1Br;	a[j2+54] = a1p1Bi;
		a[j1+14] = a1p1Cr;	a[j2+14] = a1p1Ci;
		a[j1+46] = a1p1Dr;	a[j2+46] = a1p1Di;
		a[j1+30] = a1p1Er;	a[j2+30] = a1p1Ei;
		a[j1+62] = a1p1Fr;	a[j2+62] = a1p1Fi;
		goto loop;	// Skip dyadic-mul and iFFT

	} else if (fwd_fft_only) {	// 2-input modmul: fwdFFT data dyadic-mul'ed with precomputed 2nd-vector stored in fwd-FFTed form in b[]:
	  if(fwd_fft_only == 3) {	// v20: Both inputs enter fwd-FFTed, must copy from main arrays a[],b[] to local-vars
	  // Block 1:
		a1p00r = a[j1   ];	a1p00i = a[j2   ];
		a1p01r = a[j1+32];	a1p01i = a[j2+32];
		a1p02r = a[j1+16];	a1p02i = a[j2+16];
		a1p03r = a[j1+48];	a1p03i = a[j2+48];
		a1p04r = a[j1+ 8];	a1p04i = a[j2+ 8];
		a1p05r = a[j1+40];	a1p05i = a[j2+40];
		a1p06r = a[j1+24];	a1p06i = a[j2+24];
		a1p07r = a[j1+56];	a1p07i = a[j2+56];
	  // Block 2:
		a1p08r = a[j1+ 4];	a1p08i = a[j2+ 4];
		a1p09r = a[j1+36];	a1p09i = a[j2+36];
		a1p0Ar = a[j1+20];	a1p0Ai = a[j2+20];
		a1p0Br = a[j1+52];	a1p0Bi = a[j2+52];
		a1p0Cr = a[j1+12];	a1p0Ci = a[j2+12];
		a1p0Dr = a[j1+44];	a1p0Di = a[j2+44];
		a1p0Er = a[j1+28];	a1p0Ei = a[j2+28];
		a1p0Fr = a[j1+60];	a1p0Fi = a[j2+60];
	  // Block 3:
		a1p10r = a[j1+ 2];	a1p10i = a[j2+ 2];
		a1p11r = a[j1+34];	a1p11i = a[j2+34];
		a1p12r = a[j1+18];	a1p12i = a[j2+18];
		a1p13r = a[j1+50];	a1p13i = a[j2+50];
		a1p14r = a[j1+10];	a1p14i = a[j2+10];
		a1p15r = a[j1+42];	a1p15i = a[j2+42];
		a1p16r = a[j1+26];	a1p16i = a[j2+26];
		a1p17r = a[j1+58];	a1p17i = a[j2+58];
	  // Block 4:
		a1p18r = a[j1+ 6];	a1p18i = a[j2+ 6];
		a1p19r = a[j1+38];	a1p19i = a[j2+38];
		a1p1Ar = a[j1+22];	a1p1Ai = a[j2+22];
		a1p1Br = a[j1+54];	a1p1Bi = a[j2+54];
		a1p1Cr = a[j1+14];	a1p1Ci = a[j2+14];
		a1p1Dr = a[j1+46];	a1p1Di = a[j2+46];
		a1p1Er = a[j1+30];	a1p1Ei = a[j2+30];
		a1p1Fr = a[j1+62];	a1p1Fi = a[j2+62];
	  }

		if(c_arr) {	// c_arr != 0x0: a * (b - c)

		  // Block 1:
			b1p00r = b[j1   ] - c_arr[j1   ];	b1p00i = b[j2   ] - c_arr[j2   ];
			b1p01r = b[j1+32] - c_arr[j1+32];	b1p01i = b[j2+32] - c_arr[j2+32];
			b1p02r = b[j1+16] - c_arr[j1+16];	b1p02i = b[j2+16] - c_arr[j2+16];
			b1p03r = b[j1+48] - c_arr[j1+48];	b1p03i = b[j2+48] - c_arr[j2+48];
			b1p04r = b[j1+ 8] - c_arr[j1+ 8];	b1p04i = b[j2+ 8] - c_arr[j2+ 8];
			b1p05r = b[j1+40] - c_arr[j1+40];	b1p05i = b[j2+40] - c_arr[j2+40];
			b1p06r = b[j1+24] - c_arr[j1+24];	b1p06i = b[j2+24] - c_arr[j2+24];
			b1p07r = b[j1+56] - c_arr[j1+56];	b1p07i = b[j2+56] - c_arr[j2+56];
		  // Block 2:
			b1p08r = b[j1+ 4] - c_arr[j1+ 4];	b1p08i = b[j2+ 4] - c_arr[j2+ 4];
			b1p09r = b[j1+36] - c_arr[j1+36];	b1p09i = b[j2+36] - c_arr[j2+36];
			b1p0Ar = b[j1+20] - c_arr[j1+20];	b1p0Ai = b[j2+20] - c_arr[j2+20];
			b1p0Br = b[j1+52] - c_arr[j1+52];	b1p0Bi = b[j2+52] - c_arr[j2+52];
			b1p0Cr = b[j1+12] - c_arr[j1+12];	b1p0Ci = b[j2+12] - c_arr[j2+12];
			b1p0Dr = b[j1+44] - c_arr[j1+44];	b1p0Di = b[j2+44] - c_arr[j2+44];
			b1p0Er = b[j1+28] - c_arr[j1+28];	b1p0Ei = b[j2+28] - c_arr[j2+28];
			b1p0Fr = b[j1+60] - c_arr[j1+60];	b1p0Fi = b[j2+60] - c_arr[j2+60];
		  // Block 3:
			b1p10r = b[j1+ 2] - c_arr[j1+ 2];	b1p10i = b[j2+ 2] - c_arr[j2+ 2];
			b1p11r = b[j1+34] - c_arr[j1+34];	b1p11i = b[j2+34] - c_arr[j2+34];
			b1p12r = b[j1+18] - c_arr[j1+18];	b1p12i = b[j2+18] - c_arr[j2+18];
			b1p13r = b[j1+50] - c_arr[j1+50];	b1p13i = b[j2+50] - c_arr[j2+50];
			b1p14r = b[j1+10] - c_arr[j1+10];	b1p14i = b[j2+10] - c_arr[j2+10];
			b1p15r = b[j1+42] - c_arr[j1+42];	b1p15i = b[j2+42] - c_arr[j2+42];
			b1p16r = b[j1+26] - c_arr[j1+26];	b1p16i = b[j2+26] - c_arr[j2+26];
			b1p17r = b[j1+58] - c_arr[j1+58];	b1p17i = b[j2+58] - c_arr[j2+58];
		  // Block 4:
			b1p18r = b[j1+ 6] - c_arr[j1+ 6];	b1p18i = b[j2+ 6] - c_arr[j2+ 6];
			b1p19r = b[j1+38] - c_arr[j1+38];	b1p19i = b[j2+38] - c_arr[j2+38];
			b1p1Ar = b[j1+22] - c_arr[j1+22];	b1p1Ai = b[j2+22] - c_arr[j2+22];
			b1p1Br = b[j1+54] - c_arr[j1+54];	b1p1Bi = b[j2+54] - c_arr[j2+54];
			b1p1Cr = b[j1+14] - c_arr[j1+14];	b1p1Ci = b[j2+14] - c_arr[j2+14];
			b1p1Dr = b[j1+46] - c_arr[j1+46];	b1p1Di = b[j2+46] - c_arr[j2+46];
			b1p1Er = b[j1+30] - c_arr[j1+30];	b1p1Ei = b[j2+30] - c_arr[j2+30];
			b1p1Fr = b[j1+62] - c_arr[j1+62];	b1p1Fi = b[j2+62] - c_arr[j2+62];

		} else {	// c_arr = 0x0: a * b

		  // Block 1:
			b1p00r = b[j1   ];	b1p00i = b[j2   ];
			b1p01r = b[j1+32];	b1p01i = b[j2+32];
			b1p02r = b[j1+16];	b1p02i = b[j2+16];
			b1p03r = b[j1+48];	b1p03i = b[j2+48];
			b1p04r = b[j1+ 8];	b1p04i = b[j2+ 8];
			b1p05r = b[j1+40];	b1p05i = b[j2+40];
			b1p06r = b[j1+24];	b1p06i = b[j2+24];
			b1p07r = b[j1+56];	b1p07i = b[j2+56];
		  // Block 2:
			b1p08r = b[j1+ 4];	b1p08i = b[j2+ 4];
			b1p09r = b[j1+36];	b1p09i = b[j2+36];
			b1p0Ar = b[j1+20];	b1p0Ai = b[j2+20];
			b1p0Br = b[j1+52];	b1p0Bi = b[j2+52];
			b1p0Cr = b[j1+12];	b1p0Ci = b[j2+12];
			b1p0Dr = b[j1+44];	b1p0Di = b[j2+44];
			b1p0Er = b[j1+28];	b1p0Ei = b[j2+28];
			b1p0Fr = b[j1+60];	b1p0Fi = b[j2+60];
		  // Block 3:
			b1p10r = b[j1+ 2];	b1p10i = b[j2+ 2];
			b1p11r = b[j1+34];	b1p11i = b[j2+34];
			b1p12r = b[j1+18];	b1p12i = b[j2+18];
			b1p13r = b[j1+50];	b1p13i = b[j2+50];
			b1p14r = b[j1+10];	b1p14i = b[j2+10];
			b1p15r = b[j1+42];	b1p15i = b[j2+42];
			b1p16r = b[j1+26];	b1p16i = b[j2+26];
			b1p17r = b[j1+58];	b1p17i = b[j2+58];
		  // Block 4:
			b1p18r = b[j1+ 6];	b1p18i = b[j2+ 6];
			b1p19r = b[j1+38];	b1p19i = b[j2+38];
			b1p1Ar = b[j1+22];	b1p1Ai = b[j2+22];
			b1p1Br = b[j1+54];	b1p1Bi = b[j2+54];
			b1p1Cr = b[j1+14];	b1p1Ci = b[j2+14];
			b1p1Dr = b[j1+46];	b1p1Di = b[j2+46];
			b1p1Er = b[j1+30];	b1p1Ei = b[j2+30];
			b1p1Fr = b[j1+62];	b1p1Fi = b[j2+62];

		}	// endif(c_arr)

	/*...Dyadic mul of the forward FFT outputs: (ar + I.ai)*(br + I.bi) = (ar.br - ai.bi) + I.(ar.bi + ai.br) */
		rt = a1p00r;	a1p00r = (a1p00r*b1p00r - a1p00i*b1p00i);	a1p00i = (rt*b1p00i + a1p00i*b1p00r);
		rt = a1p01r;	a1p01r = (a1p01r*b1p01r - a1p01i*b1p01i);	a1p01i = (rt*b1p01i + a1p01i*b1p01r);
		rt = a1p02r;	a1p02r = (a1p02r*b1p02r - a1p02i*b1p02i);	a1p02i = (rt*b1p02i + a1p02i*b1p02r);
		rt = a1p03r;	a1p03r = (a1p03r*b1p03r - a1p03i*b1p03i);	a1p03i = (rt*b1p03i + a1p03i*b1p03r);
		rt = a1p04r;	a1p04r = (a1p04r*b1p04r - a1p04i*b1p04i);	a1p04i = (rt*b1p04i + a1p04i*b1p04r);
		rt = a1p05r;	a1p05r = (a1p05r*b1p05r - a1p05i*b1p05i);	a1p05i = (rt*b1p05i + a1p05i*b1p05r);
		rt = a1p06r;	a1p06r = (a1p06r*b1p06r - a1p06i*b1p06i);	a1p06i = (rt*b1p06i + a1p06i*b1p06r);
		rt = a1p07r;	a1p07r = (a1p07r*b1p07r - a1p07i*b1p07i);	a1p07i = (rt*b1p07i + a1p07i*b1p07r);
		rt = a1p08r;	a1p08r = (a1p08r*b1p08r - a1p08i*b1p08i);	a1p08i = (rt*b1p08i + a1p08i*b1p08r);
		rt = a1p09r;	a1p09r = (a1p09r*b1p09r - a1p09i*b1p09i);	a1p09i = (rt*b1p09i + a1p09i*b1p09r);
		rt = a1p0Ar;	a1p0Ar = (a1p0Ar*b1p0Ar - a1p0Ai*b1p0Ai);	a1p0Ai = (rt*b1p0Ai + a1p0Ai*b1p0Ar);
		rt = a1p0Br;	a1p0Br = (a1p0Br*b1p0Br - a1p0Bi*b1p0Bi);	a1p0Bi = (rt*b1p0Bi + a1p0Bi*b1p0Br);
		rt = a1p0Cr;	a1p0Cr = (a1p0Cr*b1p0Cr - a1p0Ci*b1p0Ci);	a1p0Ci = (rt*b1p0Ci + a1p0Ci*b1p0Cr);
		rt = a1p0Dr;	a1p0Dr = (a1p0Dr*b1p0Dr - a1p0Di*b1p0Di);	a1p0Di = (rt*b1p0Di + a1p0Di*b1p0Dr);
		rt = a1p0Er;	a1p0Er = (a1p0Er*b1p0Er - a1p0Ei*b1p0Ei);	a1p0Ei = (rt*b1p0Ei + a1p0Ei*b1p0Er);
		rt = a1p0Fr;	a1p0Fr = (a1p0Fr*b1p0Fr - a1p0Fi*b1p0Fi);	a1p0Fi = (rt*b1p0Fi + a1p0Fi*b1p0Fr);
		rt = a1p10r;	a1p10r = (a1p10r*b1p10r - a1p10i*b1p10i);	a1p10i = (rt*b1p10i + a1p10i*b1p10r);
		rt = a1p11r;	a1p11r = (a1p11r*b1p11r - a1p11i*b1p11i);	a1p11i = (rt*b1p11i + a1p11i*b1p11r);
		rt = a1p12r;	a1p12r = (a1p12r*b1p12r - a1p12i*b1p12i);	a1p12i = (rt*b1p12i + a1p12i*b1p12r);
		rt = a1p13r;	a1p13r = (a1p13r*b1p13r - a1p13i*b1p13i);	a1p13i = (rt*b1p13i + a1p13i*b1p13r);
		rt = a1p14r;	a1p14r = (a1p14r*b1p14r - a1p14i*b1p14i);	a1p14i = (rt*b1p14i + a1p14i*b1p14r);
		rt = a1p15r;	a1p15r = (a1p15r*b1p15r - a1p15i*b1p15i);	a1p15i = (rt*b1p15i + a1p15i*b1p15r);
		rt = a1p16r;	a1p16r = (a1p16r*b1p16r - a1p16i*b1p16i);	a1p16i = (rt*b1p16i + a1p16i*b1p16r);
		rt = a1p17r;	a1p17r = (a1p17r*b1p17r - a1p17i*b1p17i);	a1p17i = (rt*b1p17i + a1p17i*b1p17r);
		rt = a1p18r;	a1p18r = (a1p18r*b1p18r - a1p18i*b1p18i);	a1p18i = (rt*b1p18i + a1p18i*b1p18r);
		rt = a1p19r;	a1p19r = (a1p19r*b1p19r - a1p19i*b1p19i);	a1p19i = (rt*b1p19i + a1p19i*b1p19r);
		rt = a1p1Ar;	a1p1Ar = (a1p1Ar*b1p1Ar - a1p1Ai*b1p1Ai);	a1p1Ai = (rt*b1p1Ai + a1p1Ai*b1p1Ar);
		rt = a1p1Br;	a1p1Br = (a1p1Br*b1p1Br - a1p1Bi*b1p1Bi);	a1p1Bi = (rt*b1p1Bi + a1p1Bi*b1p1Br);
		rt = a1p1Cr;	a1p1Cr = (a1p1Cr*b1p1Cr - a1p1Ci*b1p1Ci);	a1p1Ci = (rt*b1p1Ci + a1p1Ci*b1p1Cr);
		rt = a1p1Dr;	a1p1Dr = (a1p1Dr*b1p1Dr - a1p1Di*b1p1Di);	a1p1Di = (rt*b1p1Di + a1p1Di*b1p1Dr);
		rt = a1p1Er;	a1p1Er = (a1p1Er*b1p1Er - a1p1Ei*b1p1Ei);	a1p1Ei = (rt*b1p1Ei + a1p1Ei*b1p1Er);
		rt = a1p1Fr;	a1p1Fr = (a1p1Fr*b1p1Fr - a1p1Fi*b1p1Fi);	a1p1Fi = (rt*b1p1Fi + a1p1Fi*b1p1Fr);

	} else {	// fwd_fft_only = 0: Normal execution, dyadic-squaring followed by iFFT:

	/*...Dyadic square of the forward FFT outputs: */
		rt = a1p00r*a1p00i;	a1p00r = (a1p00r + a1p00i)*(a1p00r - a1p00i);	a1p00i = (rt + rt);
		rt = a1p01r*a1p01i;	a1p01r = (a1p01r + a1p01i)*(a1p01r - a1p01i);	a1p01i = (rt + rt);
		rt = a1p02r*a1p02i;	a1p02r = (a1p02r + a1p02i)*(a1p02r - a1p02i);	a1p02i = (rt + rt);
		rt = a1p03r*a1p03i;	a1p03r = (a1p03r + a1p03i)*(a1p03r - a1p03i);	a1p03i = (rt + rt);
		rt = a1p04r*a1p04i;	a1p04r = (a1p04r + a1p04i)*(a1p04r - a1p04i);	a1p04i = (rt + rt);
		rt = a1p05r*a1p05i;	a1p05r = (a1p05r + a1p05i)*(a1p05r - a1p05i);	a1p05i = (rt + rt);
		rt = a1p06r*a1p06i;	a1p06r = (a1p06r + a1p06i)*(a1p06r - a1p06i);	a1p06i = (rt + rt);
		rt = a1p07r*a1p07i;	a1p07r = (a1p07r + a1p07i)*(a1p07r - a1p07i);	a1p07i = (rt + rt);
		rt = a1p08r*a1p08i;	a1p08r = (a1p08r + a1p08i)*(a1p08r - a1p08i);	a1p08i = (rt + rt);
		rt = a1p09r*a1p09i;	a1p09r = (a1p09r + a1p09i)*(a1p09r - a1p09i);	a1p09i = (rt + rt);
		rt = a1p0Ar*a1p0Ai;	a1p0Ar = (a1p0Ar + a1p0Ai)*(a1p0Ar - a1p0Ai);	a1p0Ai = (rt + rt);
		rt = a1p0Br*a1p0Bi;	a1p0Br = (a1p0Br + a1p0Bi)*(a1p0Br - a1p0Bi);	a1p0Bi = (rt + rt);
		rt = a1p0Cr*a1p0Ci;	a1p0Cr = (a1p0Cr + a1p0Ci)*(a1p0Cr - a1p0Ci);	a1p0Ci = (rt + rt);
		rt = a1p0Dr*a1p0Di;	a1p0Dr = (a1p0Dr + a1p0Di)*(a1p0Dr - a1p0Di);	a1p0Di = (rt + rt);
		rt = a1p0Er*a1p0Ei;	a1p0Er = (a1p0Er + a1p0Ei)*(a1p0Er - a1p0Ei);	a1p0Ei = (rt + rt);
		rt = a1p0Fr*a1p0Fi;	a1p0Fr = (a1p0Fr + a1p0Fi)*(a1p0Fr - a1p0Fi);	a1p0Fi = (rt + rt);
		rt = a1p10r*a1p10i;	a1p10r = (a1p10r + a1p10i)*(a1p10r - a1p10i);	a1p10i = (rt + rt);
		rt = a1p11r*a1p11i;	a1p11r = (a1p11r + a1p11i)*(a1p11r - a1p11i);	a1p11i = (rt + rt);
		rt = a1p12r*a1p12i;	a1p12r = (a1p12r + a1p12i)*(a1p12r - a1p12i);	a1p12i = (rt + rt);
		rt = a1p13r*a1p13i;	a1p13r = (a1p13r + a1p13i)*(a1p13r - a1p13i);	a1p13i = (rt + rt);
		rt = a1p14r*a1p14i;	a1p14r = (a1p14r + a1p14i)*(a1p14r - a1p14i);	a1p14i = (rt + rt);
		rt = a1p15r*a1p15i;	a1p15r = (a1p15r + a1p15i)*(a1p15r - a1p15i);	a1p15i = (rt + rt);
		rt = a1p16r*a1p16i;	a1p16r = (a1p16r + a1p16i)*(a1p16r - a1p16i);	a1p16i = (rt + rt);
		rt = a1p17r*a1p17i;	a1p17r = (a1p17r + a1p17i)*(a1p17r - a1p17i);	a1p17i = (rt + rt);
		rt = a1p18r*a1p18i;	a1p18r = (a1p18r + a1p18i)*(a1p18r - a1p18i);	a1p18i = (rt + rt);
		rt = a1p19r*a1p19i;	a1p19r = (a1p19r + a1p19i)*(a1p19r - a1p19i);	a1p19i = (rt + rt);
		rt = a1p1Ar*a1p1Ai;	a1p1Ar = (a1p1Ar + a1p1Ai)*(a1p1Ar - a1p1Ai);	a1p1Ai = (rt + rt);
		rt = a1p1Br*a1p1Bi;	a1p1Br = (a1p1Br + a1p1Bi)*(a1p1Br - a1p1Bi);	a1p1Bi = (rt + rt);
		rt = a1p1Cr*a1p1Ci;	a1p1Cr = (a1p1Cr + a1p1Ci)*(a1p1Cr - a1p1Ci);	a1p1Ci = (rt + rt);
		rt = a1p1Dr*a1p1Di;	a1p1Dr = (a1p1Dr + a1p1Di)*(a1p1Dr - a1p1Di);	a1p1Di = (rt + rt);
		rt = a1p1Er*a1p1Ei;	a1p1Er = (a1p1Er + a1p1Ei)*(a1p1Er - a1p1Ei);	a1p1Ei = (rt + rt);
		rt = a1p1Fr*a1p1Fi;	a1p1Fr = (a1p1Fr + a1p1Fi)*(a1p1Fr - a1p1Fi);	a1p1Fi = (rt + rt);

	}	// endif(fwd_fft_only == 1)

	/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/

		#if PFETCH
		add0 = &a[j1+64];
		#endif
	/*   gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 IDIT transforms... */

		/*...Block 1: */
		#if PFETCH
		prefetch_p_doubles(add0);
		add0 += 4;
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
		prefetch_p_doubles(add0);
		#endif
		add0 += 4;
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
		prefetch_p_doubles(add0);
		add0 += 4;
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
		prefetch_p_doubles(add0);
		#endif
		add0 += 4;
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
		prefetch_p_doubles(add0);
		add0 += 4;
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
		prefetch_p_doubles(add0);
		#endif
		add0 += 4;
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
		prefetch_p_doubles(add0);
		add0 += 4;
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
		prefetch_p_doubles(add0);
		#endif
		add0 += 4;
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

		/*...and now do eight radix-4 transforms, including the internal twiddle factors:
		1, exp(-i* 1*twopi/32) =       ( c32_1,-s32_1), exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 3*twopi/32) =       ( c32_3,-s32_3) (for inputs to transform block 2),
		1, exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 3*twopi/32) =       ( s    ,-c    ) (for inputs to transform block 3),
		1, exp(-i* 3*twopi/32) =       ( c32_3,-s32_3), exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i* 9*twopi/32) =       (-s32_1,-c32_1) (for inputs to transform block 4),
		1, exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 8*twopi/32) =       ( 0    ,-1    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ) (for inputs to transform block 5),
		1, exp(-i* 5*twopi/32) =       ( s32_3,-c32_3), exp(-i*10*twopi/32) =       (-s    ,-c    ), exp(-i*15*twopi/32) =       (-c32_1,-s32_1) (for inputs to transform block 6),
		1, exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ), exp(-i*18*twopi/32) =       (-c    , s    ) (for inputs to transform block 7),
		1, exp(-i* 7*twopi/32) =       ( s32_1,-c32_1), exp(-i*14*twopi/32) =       (-c    ,-s    ), exp(-i*21*twopi/32) =       (-s32_3, c32_3) (for inputs to transform block 8),
		and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.	*/

		/*...Block 1: t00,t10,t20,t30	*/
		#if PFETCH
		prefetch_p_doubles(add0);
		add0 += 4;
		#endif
		rt =t10;	t10=t00-rt;	t00=t00+rt;
		it =t11;	t11=t01-it;	t01=t01+it;

		rt =t30;	t30=t20-rt;	t20=t20+rt;
		it =t31;	t31=t21-it;	t21=t21+it;

		a[j1   ]=t00+t20;			a[j2   ]=t01+t21;
		t00       =t00-t20;				t01       =t01-t21;
		a[j1+32]=t00*c10+t01*s10;	a[j2+32]=t01*c10-t00*s10;

		rt        =t10+t31;				it        =t11-t30;
		t10       =t10-t31;				t11       =t11+t30;
		a[j1+16]=rt *c08+it *s08;	a[j2+16]=it *c08-rt *s08;
		a[j1+48]=t10*c18+t11*s18;	a[j2+48]=t11*c18-t10*s18;

		/*...Block 5: t08,t18,t28,t38	*/
		#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(add0);
		#endif
		add0 += 4;
		#endif
		rt =t18;	t18=t08-t19;		t08=t08+t19;
					t19=t09+rt;			t09=t09-rt;

		rt =(t29+t28)*ISRT2;			t29=(t29-t28)*ISRT2;		t28=rt;
		rt =(t38-t39)*ISRT2;			it =(t38+t39)*ISRT2;
		t38=t28+rt;						t28=t28-rt;
		t39=t29+it;						t29=t29-it;

		rt        =t08+t28;				it        =t09+t29;
		t08       =t08-t28;				t09       =t09-t29;
		a[j1+8 ]=rt *c04+it *s04;	a[j2+8 ]=it *c04-rt *s04;
		a[j1+40]=t08*c14+t09*s14;	a[j2+40]=t09*c14-t08*s14;

		rt        =t18+t39;				it        =t19-t38;
		t18       =t18-t39;				t19       =t19+t38;
		a[j1+24]=rt *c0C+it *s0C;	a[j2+24]=it *c0C-rt *s0C;
		a[j1+56]=t18*c1C+t19*s1C;	a[j2+56]=t19*c1C-t18*s1C;

		/*...Block 3: t04,t14,t24,t34	*/
		#if PFETCH
		prefetch_p_doubles(add0);
		add0 += 4;
		#endif
		rt =(t15+t14)*ISRT2;			it =(t15-t14)*ISRT2;
		t14=t04-rt;						t04=t04+rt;
		t15=t05-it;						t05=t05+it;

		rt =t24*c + t25*s;				t25=t25*c - t24*s;		t24=rt;
		rt =t34*s + t35*c;				it =t35*s - t34*c;
		t34=t24-rt;						t24=t24+rt;
		t35=t25-it;						t25=t25+it;

		rt        =t04+t24;				it        =t05+t25;
		t04       =t04-t24;				t05       =t05-t25;
		a[j1+4 ]=rt *c02+it *s02;	a[j2+4 ]=it *c02-rt *s02;
		a[j1+36]=t04*c12+t05*s12;	a[j2+36]=t05*c12-t04*s12;

		rt        =t14+t35;				it        =t15-t34;
		t14       =t14-t35;				t15       =t15+t34;
		a[j1+20]=rt *c0A+it *s0A;	a[j2+20]=it *c0A-rt *s0A;
		a[j1+52]=t14*c1A+t15*s1A;	a[j2+52]=t15*c1A-t14*s1A;

		/*...Block 7: t0C,t1C,t2C,t3C	*/
		#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(add0);
		#endif
		add0 += 4;
		#endif
		rt =(t1C-t1D)*ISRT2;			it =(t1C+t1D)*ISRT2;
		t1C=t0C+rt;						t0C=t0C-rt;
		t1D=t0D+it;						t0D=t0D-it;

		rt =t2C*s + t2D*c;				t2D=t2D*s - t2C*c;		t2C=rt;
		rt =t3C*c + t3D*s;				it =t3D*c - t3C*s;
		t3C=t2C+rt;						t2C=t2C-rt;
		t3D=t2D+it;						t2D=t2D-it;

		rt        =t0C+t2C;				it        =t0D+t2D;
		t0C       =t0C-t2C;				t0D       =t0D-t2D;
		a[j1+12]=rt *c06+it *s06;	a[j2+12]=it *c06-rt *s06;
		a[j1+44]=t0C*c16+t0D*s16;	a[j2+44]=t0D*c16-t0C*s16;

		rt        =t1C+t3D;				it        =t1D-t3C;
		t1C       =t1C-t3D;				t1D       =t1D+t3C;
		a[j1+28]=rt *c0E+it *s0E;	a[j2+28]=it *c0E-rt *s0E;
		a[j1+60]=t1C*c1E+t1D*s1E;	a[j2+60]=t1D*c1E-t1C*s1E;

		/*...Block 2: t02,t12,t22,t32	*/
		#if PFETCH
		prefetch_p_doubles(add0);
		add0 += 4;
		#endif
		rt =t12*c + t13*s;				it =t13*c - t12*s;
		t12=t02-rt;						t02=t02+rt;
		t13=t03-it;						t03=t03+it;

		rt =t22*c32_1 + t23*s32_1;		t23=t23*c32_1 - t22*s32_1;	t22=rt;
		rt =t32*c32_3 + t33*s32_3;		it =t33*c32_3 - t32*s32_3;
		t32=t22-rt;						t22=t22+rt;
		t33=t23-it;						t23=t23+it;

		rt        =t02+t22;				it        =t03+t23;
		t02       =t02-t22;				t03       =t03-t23;
		a[j1+2 ]=rt *c01+it *s01;	a[j2+2 ]=it *c01-rt *s01;
		a[j1+34]=t02*c11+t03*s11;	a[j2+34]=t03*c11-t02*s11;

		rt        =t12+t33;				it        =t13-t32;
		t12       =t12-t33;				t13       =t13+t32;
		a[j1+18]=rt *c09+it *s09;	a[j2+18]=it *c09-rt *s09;
		a[j1+50]=t12*c19+t13*s19;	a[j2+50]=t13*c19-t12*s19;

		/*...Block 6: t0A,t1A,t2A,t3A	*/
		#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(add0);
		#endif
		add0 += 4;
		#endif
		rt =t1A*s - t1B*c;				it =t1B*s + t1A*c;
		t1A=t0A+rt;						t0A =t0A-rt;
		t1B=t0B+it;						t0B =t0B-it;

		rt =t2A*s32_3 + t2B*c32_3;		t2B=t2B*s32_3 - t2A*c32_3;	t2A=rt;
		rt =t3A*c32_1 - t3B*s32_1;		it =t3B*c32_1 + t3A*s32_1;
		t3A=t2A+rt;						t2A=t2A-rt;
		t3B=t2B+it;						t2B=t2B-it;

		rt        =t0A+t2A;				it        =t0B+t2B;
		t0A       =t0A-t2A;				t0B       =t0B-t2B;
		a[j1+10]=rt *c05+it *s05;	a[j2+10]=it *c05-rt *s05;
		a[j1+42]=t0A*c15+t0B*s15;	a[j2+42]=t0B*c15-t0A*s15;

		rt        =t1A+t3B;				it        =t1B-t3A;
		t1A       =t1A-t3B;				t1B       =t1B+t3A;
		a[j1+26]=rt *c0D+it *s0D;	a[j2+26]=it *c0D-rt *s0D;
		a[j1+58]=t1A*c1D+t1B*s1D;	a[j2+58]=t1B*c1D-t1A*s1D;

		/*...Block 4: t06,t16,t26,t36	*/
		#if PFETCH
		prefetch_p_doubles(add0);
		add0 += 4;
		#endif
		rt =t16*s + t17*c;				it =t17*s - t16*c;
		t16=t06-rt;						t06 =t06+rt;
		t17=t07-it;						t07 =t07+it;

		rt =t26*c32_3 + t27*s32_3;		t27=t27*c32_3 - t26*s32_3;	t26=rt;
		rt =t36*s32_1 - t37*c32_1;		it =t37*s32_1 + t36*c32_1;
		t36=t26+rt;						t26=t26-rt;
		t37=t27+it;						t27=t27-it;

		rt        =t06+t26;				it        =t07+t27;
		t06       =t06-t26;				t07       =t07-t27;
		a[j1+6 ]=rt *c03+it *s03;	a[j2+6 ]=it *c03-rt *s03;
		a[j1+38]=t06*c13+t07*s13;	a[j2+38]=t07*c13-t06*s13;

		rt        =t16+t37;				it        =t17-t36;
		t16       =t16-t37;				t17       =t17+t36;
		a[j1+22]=rt *c0B+it *s0B;	a[j2+22]=it *c0B-rt *s0B;
		a[j1+54]=t16*c1B+t17*s1B;	a[j2+54]=t17*c1B-t16*s1B;

		/*...Block 8: t0E,t1E,t2E,t3E	*/
		#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
		prefetch_p_doubles(add0);
		#endif
		add0 += 4;
		#endif
		rt =t1E*c - t1F*s;				it =t1F*c + t1E*s;
		t1E=t0E+rt;						t0E =t0E-rt;
		t1F=t0F+it;						t0F =t0F-it;

		rt =t2E*s32_1 + t2F*c32_1;		t2F=t2F*s32_1 - t2E*c32_1;	t2E=rt;
		rt =t3E*s32_3 + t3F*c32_3;		it =t3F*s32_3 - t3E*c32_3;
		t3E=t2E+rt;						t2E=t2E-rt;
		t3F=t2F+it;						t2F=t2F-it;

		rt        =t0E+t2E;				it        =t0F+t2F;
		t0E       =t0E-t2E;				t0F       =t0F-t2F;
		a[j1+14]=rt *c07+it *s07;	a[j2+14]=it *c07-rt *s07;
		a[j1+46]=t0E*c17+t0F*s17;	a[j2+46]=t0F*c17-t0E*s17;

		rt        =t1E+t3F;				it        =t1F-t3E;
		t1E       =t1E-t3F;				t1F       =t1F+t3E;
		a[j1+30]=rt *c0F+it *s0F;	a[j2+30]=it *c0F-rt *s0F;
		a[j1+62]=t1E*c1F+t1F*s1F;	a[j2+62]=t1F*c1F-t1E*s1F;

#endif	/* USE_SSE2 */
	loop:
		continue;	// Need a statement following the jump-target label 'loop' here
	}	/* endfor() */
//if(thr_id < 2 || thr_id == (radix0-1))printf("%s: on-exit index0,1_idx = %u,%u\n",func,index0_idx,index1_idx);
#ifdef MULTITHREAD
	// Write exit values of these 2 ints back to thread-associated storage:
//	*idx0_ptr = index0_idx; *idx1_ptr = index1_idx;
#endif
}

