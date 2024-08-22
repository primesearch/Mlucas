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

// In SSE2 mode this variant is unambiguously faster; for ARMv8 (which also sets USE_SSE2) it is the only implemented 16-DFT.
// Tangent-DIF is default version in avx512 mode ... compile with -DREFACTOR_4DFT_3TWIDDLE to force non-tangent:
#if defined(USE_SSE2) && (defined(USE_IMCI512) || !defined(USE_AVX512))	// Use non-tangent for 1st-gen Xeon Phi
	#define REFACTOR_4DFT_3TWIDDLE
#endif
#ifdef REFACTOR_4DFT_3TWIDDLE
	#include "sse2_macro_gcc64.h"
#endif

#define RADIX 16

/* Use for testing higher-accuracy version of the twiddles computation */
#define HIACC 1
#define EPS	1e-10

#ifdef USE_AVX512
	#define DFT_V1
#elif defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE)
	#if !defined(DFT_V1) && !defined(DFT_V2)
		#define DFT_V1	// Compile-time toggle between 2 versions of the FMA-based DFT macros - Must
		// define one of these for AVX+FMA builds, default is V1 (slightly faster on my Haswell 4670)
	#endif
#endif

#ifndef PFETCH_DIST
  #ifdef USE_AVX
	#define PFETCH_DIST	4096	// This seems to work best on my Haswell
  #else
	#define PFETCH_DIST	32
  #endif
#endif

#ifdef USE_SSE2

	#ifdef REFACTOR_4DFT_3TWIDDLE
		#warning REFACTOR_4DFT_3TWIDDLE
	#elif defined(DFT_V2)
		#warning DFT_V2
	#else
		#warning DFT_V1
	#endif

  #ifndef COMPILER_TYPE_GCC
	#error SSE2 code not supported for this compiler!
  #endif

  #if(HIACC != 1)
	#error SIMD Mode requires HIACC flag to be set!
  #endif

	#include "radix16_dif_dit_pass_asm.h"

#endif

/***************/

/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform a single radix-16 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
*/
#ifdef USE_FGT61
void radix16_dif_pass	(double a[], uint64 b[], int n, struct complex rt0[], struct complex rt1[], uint128 mt0[], uint128 mt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
#else
void radix16_dif_pass	(double a[],             int n, struct complex rt0[], struct complex rt1[],                               int index[], int nloops, int incr, int init_sse2, int thr_id)
#endif
{
#ifdef USE_SSE2
	const int pfetch_dist = PFETCH_DIST;
	int pfetch_addr;	// Since had pre-existing L1-targeting pfetch here, add numerical suffix to differentiate cache levels being targeted by the various types of prefetching
	static int max_threads = 0;
#endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	// lg(stride):
	const int l2_stride = L2_SZ_VD-2;	// 16 doubles at a time
#endif
#ifdef USE_FGT61
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull, q2=q+q, q3=q2+q, q4=q2+q2, q5=q4+q;	// q = 2^61 - 1, and needed small multiples
	// primitive 16th root of unity, scaled by *8:
	const uint64 cm = 1693317751237720973ull<<3, sm = 2283815672160731785ull<<3;
#endif
	const double c = 0.9238795325112867561281831	/* exp[i*(twopi/16)]*/
#if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
		,s = 0.3826834323650897717284599			/* exp[i*(twopi/16)]*/
#endif
		;
#if defined(USE_SCALAR_DFT_MACRO) || defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE)	// FMA-based DFT needs the tangent
	const double tan = 0.41421356237309504879;
#endif
	int i,j,j1,jlo,jhi,m,iroot_prim,iroot,k1,k2;
#ifndef USE_SSE2
	int j2;
#endif
	int p1,p2,p3,p4,p8,p12;
#if !defined(USE_SSE2) || defined(USE_AVX512)
	double rt,it;
#endif
#if defined(MULTITHREAD) && defined(USE_SSE2)
	double dtmp;
#endif
	double re0,im0,re1,im1;
#if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
	// These needed both for scalar mode and for certain SIMD-mode inits:
	double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;
#endif

#ifdef USE_SSE2
  #if 0
	#error Experimental code: needs debug!
	// Experimental: Coeffs of 6-term [even|odd] Chebyshev approximation to [cos(x)|sin(x)], x in [0,Pi/4], IEEE64-doubles stored as uint64 bitfields:
	const uint64 cheb_c[6] = {
					0x3feffffffffffe0bull,	// d[ 0]
					0xbfdffffffffe3763ull,	// d[ 2]
					0x3fa5555554476425ull,	// d[ 4]
					0xbf56c16b2d3e61f5ull,	// d[ 6]
					0x3efa00e9685da9c4ull,	// d[ 8]
					0xbe923c5c15441a00ull};	// d[10]
	const uint64 cheb_s[6] = {
					0x3fefffffffffffd9ull,	// d[ 1]
					0xbfc5555555550efeull,	// d[ 3]
					0x3f81111110bde5e3ull,	// d[ 5]
					0xbf2a019f8a207d3full,	// d[ 7]
					0x3ec71d7317b8ee33ull,	// d[ 9]
					0xbe5a9507f3711e2dull};	// d[11]
  #endif
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0;	/* Addresses into array sections */

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
  #endif
  #ifndef MULTITHREAD
	static
  #endif
	  vec_dbl *cc0,*isrt2,*two,*pi4, *r1
	  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
		,*ss0
	  #endif
	  #if !defined(MULTITHREAD) || (!(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE)))
		,*c0thru15
	  #endif
		;
	// This variable might have to be reenabled if commented out Chebyshev code below is ever enabled.
	// In this case the variable must be 'static' if !MULTITHREAD.
  #if defined(MULTITHREAD) && !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
	uint64 *sign_mask;
  #endif

#else

  #ifndef USE_SCALAR_DFT_MACRO
	int jp,jt;
	double *addr, *addp;
	int prefetch_offset;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
  #endif
  #ifdef USE_FGT61
	uint64 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15
		  ,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32;
  #endif

#endif

#ifdef USE_SSE2

	/* In order to eschew complex thread-block-and-sync logic related to the local-store-init step in multithread mode,
	switch to a special-init-mode-call paradigm, in which the function is inited once (for as many threads as needed)
	prior to being executed:
	*/
	if(init_sse2)	// Just check (init_sse2 != 0) here, to allow the *value* of init_sse2 to store #threads
	{
	//	fprintf(stderr, "radix16_dif_dit_pass pfetch_dist = %d\n", pfetch_dist);
		max_threads = init_sse2;
	#ifndef COMPILER_TYPE_GCC
		ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(thr_id == -1, "Init-mode call must be outside of any multithreading!");
		if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)sc_arr);	sc_arr=0x0;
		}
		// v19 alloc'ed 72* ... v20 needs [1+1+4+8] = 18 more slots in SSE2 mode, [1+1+2+4] = 8 more in AVX/AVX2 mode, [1+1+1+2] = 5 more in AVX-512 mode,
		// just use 20 more slots in all cases for simplicity's sake. Further add 12 slots for doubled-into-vectors 6-term Chebyshev expansions of cos, sin:
		sc_arr = ALLOC_VEC_DBL(sc_arr, 104*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	last 30 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	v20: Add 16 slots for auxiliary data needed for vectorized twiddle computation.
	*/
		#ifdef MULTITHREAD
	//	if(max_threads > 1) {
			__r0     = sc_ptr;
			isrt2    = sc_ptr + 0x20;
			cc0      = sc_ptr + 0x21;
		  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
			ss0      = sc_ptr + 0x22;
		  #endif
			two      = sc_ptr + 0x43;
			/* v20:
				Vector	#SIMD slots needed:
				Const:		SSE2	AVX	AVX-512	parametrized:
				sign_mask	4		??		??
				2*pi/n		1		1		1	1
				0-15flt		4		2		1	1 << (6-L2_SZ_VD)
				0-15dbl		8		4		2	1 << (7-L2_SZ_VD)
			*/
		  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
			sign_mask = (uint64 *)(sc_ptr + 0x48);	// Start v20-added data after v19-alloc block; leave a few slots below the v20 stuff for spills, etc.
		  #endif
			pi4      = sc_ptr + 0x4c;
		  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
			c0thru15 = sc_ptr + 0x4d;	// First 1 << (6-L2_SZ_VD) slots hold (float)0-15, next 1 << (7-L2_SZ_VD) slots hold (double)0-15
		  #endif
			/* Chebyshev-expansion coeffs for cos & sin start at sc_ptr + 0x58 */
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);
			  #if defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE)
				VEC_DBL_INIT(two  , 1.0);	// Yeah, I *know* it's a misnomer :)
				// cc0,ss0 inited below for AVX2
			  #else
				VEC_DBL_INIT(two  , 2.0);
				VEC_DBL_INIT(cc0  , c);
				VEC_DBL_INIT(ss0  , s);
				// Set up the SIMD-tupled constants used by the vectorized twiddles-on-the-fly computation:
																								// lo64,hi64:
				*(sign_mask+0) = 0x7FFFFFFFFFFFFFFFull; *(sign_mask+1) = 0x7FFFFFFFFFFFFFFFull;	// +,+
				*(sign_mask+2) = 0xFFFFFFFFFFFFFFFFull; *(sign_mask+3) = 0x7FFFFFFFFFFFFFFFull;	// -,+
				*(sign_mask+4) = 0x7FFFFFFFFFFFFFFFull; *(sign_mask+5) = 0xFFFFFFFFFFFFFFFFull;	// +,-
				*(sign_mask+6) = 0xFFFFFFFFFFFFFFFFull; *(sign_mask+7) = 0xFFFFFFFFFFFFFFFFull;	// -,-
				// 128-bit memloc just below sign_mask reserved for spill of bytewise io[]-vec; next-lower 2 slots hold byte-mask vectors:
				*(sign_mask-4) = 0x1111111111111111ull; *(sign_mask-3) = 0ull;	// 0x1111111111111111 in lo64, 0 in hi64
				*(sign_mask-6) = 0x0000000f0000000full; *(sign_mask-5) = 0ull;	// 0x0000000f0000000f in lo64, 0 in hi64
				VEC_DBL_INIT(pi4   , (double) 0.78539816339744830961  );	// pi/4
			//	VEC_DBL_INIT(twopin, (double)12.56637061435917295376/n);	*** This must wait until regular thread-exec (i.e. non-init-mode) because in init mode, n is not supplied ... just combine with (double)0-15 per-thread init below
				// (float)0-15, (double)0-15:								GDB needs these, else e.g.'p *(vec_flt *)(sc_ptr+0x4a)' gives 'No symbol "vec_flt" in current context':
				float *flt_ptr = (float *)c0thru15;							//	vec_flt *vec_flt_ptr = (vec_flt *)flt_ptr;
			//	double*dbl_ptr = (double*)(c0thru15 + (1 << (6-L2_SZ_VD)));	//	vec_dbl *vec_dbl_ptr = (vec_dbl *)dbl_ptr;
				for(j = 0; j < 16; ++j) {
					*flt_ptr++ = (float)j;
			//		*dbl_ptr++ = (double)j;	<*** Just init these slots to [0-15]*twopin at same runtime-spot below - in fact instead of - where we init twopin
				}
			   #if 0
				#error Experimental code: needs debug!
				// Chebyshev-expansion coeffs for cos & sin start at (sign_mask+18) = (sc_ptr + 0x5a)...
				// use the uint64 sign_mask pointer to init const-vec_dbl as paired-uint64 bitfields:
				for(j =  0; j < 12; j+=2) { *(sign_mask+j+18) = *(sign_mask+j+19) = cheb_c[j>>1]; }
				for(j = 12; j < 24; j+=2) { *(sign_mask+j+18) = *(sign_mask+j+19) = cheb_s[j>>1]; }
			   #endif
			  #endif
				isrt2 += 104;	/* Move on to next thread's local store */
				cc0   += 104;
			  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
				ss0   += 104;
			  #endif
				two   += 104;
				pi4       += 104;
			  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
				c0thru15  += 104;
				sign_mask += 104;
			  #endif
			}

		#elif defined(COMPILER_TYPE_GCC)
			r1  = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;
										cc0 = sc_ptr + 0x21;
		  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
										ss0 = sc_ptr + 0x22;
		  #endif
										two = sc_ptr + 0x43;
		  #if 0 // for experimental code below
			sign_mask = (uint64 *)(sc_ptr + 0x48);	// Start v20-added data after v19-alloc block; leave a few slots below the v20 stuff for spills, etc.
		  #endif
			pi4      = sc_ptr + 0x4c;
			c0thru15 = sc_ptr + 0x4d;	// First 1 << (6-L2_SZ_VD) slots hold (float)0-15, next 1 << (7-L2_SZ_VD) slots hold (double)0-15
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
		  #if defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE)
			VEC_DBL_INIT(two  , 1.0);
			// cc0,ss0 inited below for AVX2
		  #else
			VEC_DBL_INIT(two  , 2.0);
			VEC_DBL_INIT(cc0  , c);
			VEC_DBL_INIT(ss0  , s);
		  #endif
			// Set up the SIMD-tupled constants used by the vectorized twiddles-on-the-fly computation:
			VEC_DBL_INIT(pi4   , (double) 0.78539816339744830961  );	// pi/4
			float *flt_ptr = (float *)c0thru15;							//	vec_flt *vec_flt_ptr = (vec_flt *)flt_ptr;
		//	double*dbl_ptr = (double*)(c0thru15 + (1 << (6-L2_SZ_VD)));	//	vec_dbl *vec_dbl_ptr = (vec_dbl *)dbl_ptr;
			for(j = 0; j < 16; ++j) {
				*flt_ptr++ = (float)j;
		//		*dbl_ptr++ = (double)j;	<*** Just init these slots to [0-15]*twopin at same runtime-spot below - in fact instead of - where we init twopin
			}
		   #if 0
			#error Experimental code: needs debug!
			// Chebyshev-expansion coeffs for cos & sin start at (sign_mask+18) = (sc_ptr + 0x5a)...
			// use the uint64 sign_mask pointer to init const-vec_dbl as paired-uint64 bitfields:
			for(j =  0; j < 12; j+=2) { *(sign_mask+j+18) = *(sign_mask+j+19) = cheb_c[j>>1]; }
			for(j = 12; j < 24; j+=2) { *(sign_mask+j+18) = *(sign_mask+j+19) = cheb_s[j>>1]; }
		   #endif
		#else
		  #error Non-GCC-compatible compilers not supported for SIMD builds!
		#endif
	//	}
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		dtmp = (double)12.56637061435917295376/n;	// twopin = 2*pi/[complex FFT length] = 2*pi/(n/2) = 4*pi/n
		r1 = __r0 + thr_id*104;
		isrt2 = r1 + 0x20;
		cc0   = r1 + 0x21;
		two   = r1 + 0x43;
		// Start v20-added data after v19-alloc block; unnamed sign_mask ptr is at += 0x48 and uses 4 slots:
		pi4   = r1 + 0x4c;
		// [0-15]*twopin:
		double*dbl_ptr = (double*)(pi4 + 1 + (1 << (6-L2_SZ_VD)));
		re0 = 4*dtmp; *(dbl_ptr+3) = 3*dtmp; *(dbl_ptr+2) = dtmp+dtmp; *(dbl_ptr+1) = dtmp; *dbl_ptr = 0.0;
		for(j = 0; j < 12; j += 4) {
			dbl_ptr += 4;
			*(dbl_ptr+0) = *(dbl_ptr-4) + re0;
			*(dbl_ptr+1) = *(dbl_ptr-3) + re0;
			*(dbl_ptr+2) = *(dbl_ptr-2) + re0;
			*(dbl_ptr+3) = *(dbl_ptr-1) + re0;
		}
	#endif

#endif

/*
!   Here's how our padded-array indexing works: given an unpadded-array index J,
!   we convert that to a padded-array index J1 using the formula
!
!	J1 = floor[J/(NBLOCK)] x [NBLOCK+PAD] + mod[J, NBLOCK],
!
!   where NBLOCK is the number of 8-byte array data in each contiguous-data block = 2^DAT_BITS,
!   and PAD is the number of 8-byte padding data between data block = 2^PAD_BITS.
!   Since N is an integer, we use the fact that mod[J, NBLOCK] = J - [NBLOCK] x floor[J/(NBLOCK)] to derive
!
!	J1 = J + PAD x floor[J/(NBLOCK)],
!
!   which expresses the fact that J1-J is PAD times the number of padding blocks we skip over.
!
!   Since NBLOCK is a power of 2, i.e. DAT_BITS = log_2(NBLOCK) is an integer, then floor[J/(NBLOCK)] = J >> DAT_BITS,
!   and J1 can be calculated using a shift, an integer multiply, and an add:
!
!	J1 = J + PAD x (J >> DAT_BITS).
!
!   If both NBLOCK and PAD are powers of 2 (as here), this can be done using two shifts and an add:
!
!	J1 = J + (J >> DAT_BITS) << PAD_BITS.
!
!   If PAD is 4 or 8, this can be done using a shift and a scaled add on Alpha.
!
!   Note that we only need to do this kind of calculation for the first elements of each block of 16 complex data.
!
!...The radix-16 pass is here. The data are processed in NLOOPS blocks of INCR real elements each
!   stride between elements processed together.
*/
	p1 = incr >> 4;
	p2 = p1 +p1;
	p3 = p2 +p1;
	p4 = p3 +p1;
	p8 = p4 +p4;
	p12= p8 +p4;

	p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
	p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
	/* Make sure to at least one-time pre-test any index arithmetic assumptions used to save cycles in the loop
	body (both C and ASM). Since such checks may be runlength-dependent, need to be cheap enough to leave on
	all the time, as here where we do them just once prior to entering the processing loop. Since DIF and DIT
	encounter the same sets or index strides (albeit in opposite order), can split such tests between them:
	*** 2014: Failure of this assertion led me to find dependence on it in my new AVX2/FMA-based DIT macro ***
	*** [But fix obviates said dependence, so no longer appropriate to enforce it.] ***
		ASSERT(p2  == p1+p1, "radix16_dif_pass: p2  != p1+p1!");
	*/
	iroot_prim=(incr >> 5);		/* (incr/2)/radix_now */
	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16).	 */
	{				/* NLOOPS satisfies the identity NLOOPS * INCR = P*N, which says that we process the entire */
					/* array each time this subroutine is executed (since P*N = vector length, sans padding.)   */

	/*	here are the needed sincos data - these are processed below in bit-reversed order. */
	  iroot = index[m]*iroot_prim;	// We are doing a length-N/2 complex FFT, and a
	  i = iroot;					// given value of iroot corresponds to angle theta = iroot*(2*Pi/(N/2))

/* In the DIF pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks.

v20:
	For Mlucas real-FFT of n doubles => cFFT of length n/2, which equals the number of equal-sized slices the interval [0,2*Pi)
	is subdivided into for computing twiddles. Thus, given an index [i] used for lookup into precomputed roots tables via
			k1=(iroot & NRTM1);	k2=(iroot >> NRT_BITS);
	the quadrant in which the corresponding twiddle lies is given by iq = i / (n/8).

	Exploitable symmetries which can be used to cut size of theta, based on the quadrant. Assume theta in [0, 2*pi).
	Let iq = int(theta/(pi/2) = 0-3 be the quadrant index, and gamma := theta - iq*pi/2 in [0, pi/2). Then:
									Precompute 2-arrays scos[0,1] = [cos(gamma),sin(gamma)] and sign[0,1] = [+1,-1], then Re,Im terms
	i	cos(theta)	sin(theta)		of output twiddle are [Re,Im] = [sign[s0]*scos[j],sign[s1]*scos[j^1]], where signs and index j
	-	-----------	-----------		are computed as:
	0	+cos(gamma)	+sin(gamma)			s0 = ((unsigned int)iq - 1) < 2
	1	-sin(gamma)	+cos(gamma)			s1 = iq > 1
	2	-cos(gamma)	-sin(gamma)			j = IS_ODD(iq)
	3	+sin(gamma)	-cos(gamma)

	We can even go one better by using the fact that all roots can be mapped to data in the first octant.
	E.g. for gamma in [pi/4,pi/2), let delta = (pi/2-gamma), thus delta in [0,pi/4), use cos,sin(gamma) = sin,cos(delta).
	Thus we have the following refinement based on octants - again precompute 2-arrays scos[0,1] = [cos(gamma),sin(gamma)]
	and sign[0,1] = [+1,-1].  Let io = int(theta/(pi/4) = 0-7 be the octant index, and

		gamma = theta - io*pi/4 if io even,
		gamma = ((io+1)*pi/4 - theta if io odd.
	We can handle both of these cases in branchless fashion by defining iodd = IS_ODD(io) and computing
		gamma = sign[iodd]*(theta - (io+iodd)*pi/4) .

	io	cos(theta)	sin(theta)		Then output twiddle = [Re,Im] = [sign[s0]*scos[j],sign[s1]*scos[j^1]], where signs and index j
	--	-----------	-----------		are computed as:
	0	+cos(gamma)	+sin(gamma)			s0 = ((unsigned int)io - 2) < 4
	1	+sin(gamma)	+cos(gamma)			s1 = io > 3
	2	-sin(gamma)	+cos(gamma)			j = IS_ODD((io+1)/2)
	3	-cos(gamma)	+sin(gamma)
	4	-cos(gamma)	-sin(gamma)
	5	-sin(gamma)	-cos(gamma)
	6	+sin(gamma)	-cos(gamma)
	7	+cos(gamma)	-sin(gamma)

For our SIMD implementation, however, the s0 computation is problemeatic because SSE2 only supports packed-compared for signed ints.
Alternative to s0 = ((unsigned int)io - 2) < 4 which works even if 16 io-values concatenated as nibbles of a uint64, call it io_cat:

	io		io+2	(io+2)>>2
	000		0010	00
	001		0011	00
	010*	0100	01	<*** s0 = low bit of these ***
	011*	0101	01
	100*	0110	01
	101*	0111	01
	110		1000	10
	111		1001	10
	uint64 io_cat = ([16 io-values concatenated as nibbles] + 0x2222222222222222ull) >> 2;

BUT, how to efficiently concatenate original 4-packed-dword vectors into this form?
Better to leave as 128-bit vectors and focus on concatenating our initial 4 vectors-of-dwords into a single vector-of-bytes;
notation below is low-to-high-[byte|word] within xmm-regs; '|' denotes dword boundaries:

	xmm0 = [w0|w1|w2|w3]
	xmm1 = [w4|w5|w6|w7]
	xmm2 = [w8|w9|wa|wb]
	xmm3 = [wc|wd|we|wf]
4x4 matrix-transpose the above:
	0,1,2,3		0,4,8,c
	4,5,6,7	==>	1,5,9,d
	8,9,a,b		2,6,a,e
	c,d,e,f		3,7,b,f

	movaps	xmm0,xmm4	// copy of xmm0
	punpckldq xmm1,xmm0	// xmm0 = [w0|w4|w1|w5]
	punpckhdq xmm1,xmm4	// xmm4 = [w2|w6|w3|w7]
	movaps	xmm2,xmm5	// copy of xmm2
	punpckldq xmm3,xmm2	// xmm2 = [w8|wc|w9|wd]
	punpckhdq xmm3,xmm5	// xmm5 = [wa|we|wb|wf]

	movaps	xmm0,xmm1	// copy of xmm0
	punpcklqdq xmm2,xmm0	// xmm0 = [w0|w4|w8|wc]
	punpcklqdq xmm2,xmm1	// xmm1 = [w1|w5|w9|wd]
	movaps	xmm4,xmm3	// copy of xmm4
	punpcklqdq xmm5,xmm4	// xmm4 = [w2|w6|wa|we]
	punpcklqdq xmm5,xmm3	// xmm3 = [w3|w7|wb|wf]
	movaps	xmm4,xmm2	// reg-copy xmm4 -> xmm2 to 'neatify' indexing
	// Now concatenate bytes:
	pslldq	$1,xmm1	// xmm1 = [ 0,w1, 0, 0| 0,w5, 0, 0| 0,w9, 0, 0| 0,wd, 0, 0]
	pslldq	$2,xmm2	// xmm2 = [ 0, 0,w2, 0| 0, 0,w6, 0| 0, 0,wa, 0| 0, 0,we, 0]
	pslldq	$3,xmm3	// xmm3 = [ 0, 0, 0,w3| 0, 0, 0,w7| 0, 0, 0,wb| 0, 0, 0,wf]
	paddd xmm1,xmm0
	paddd xmm3,xmm2
	paddd xmm2,xmm0	// xmm0 = [w0,w1,w2,w3,w4,w5,w6,w7,w8,w9,wa,wb,wc,wd,we,wf]

	But note, ultimately we need to extract pairs of adjacent-indexed bytes from any such packed-vetor-of-bytes
	result and up-convert into packed-doubles. For the above bytes-in-order layout, can use PEXTRW [3 cycle,1/cycle]
	to extract a target bytepair (PEXTRB needs SSE4.1) int either a 64-bit GPR or memloc. However, copying said bytes
	into 4-byte-apart slots of an xmm is a headache. But in fact there is no need fo any of this awkwardness - just
	skip all the matrix-transpose crap and bytewise-shift the latter 3 of
		xmm0 = [w0|w1|w2|w3]
		xmm1 = [w4|w5|w6|w7]
		xmm2 = [w8|w9|wa|wb]
		xmm3 = [wc|wd|we|wf]
	to obtain 16-byte vectors
		[0,-,-,-,1,-,-,-,2,-,-,-,3,-,-,-]
		[-,4,-,-,-,5,-,-,-,6,-,-,-,7,-,-]
		[-,-,8,-,-,-,9,-,-,-,a,-,-,-,b,-]
		[-,-,-,c,-,-,-,d,-,-,-,e,-,-,-,f]
	which we add together to get byte-vec [0,4,8,c,1,5,9,d,2,6,a,e,3,7,b,f]. In order to use CVTDQ2PD to convert
	selected byte-pairs like 4,5 to packed-doubles, just bytewise-shift as needed and AND with byte-vec
		[*,0,0,0,*,0,0,0,0,0,0,0,0,0,0,0}, where * denotes a byte == 11111111.
*/
#if HIACC
  #if 0
	#error Experimental code, not ready for release!
	// v20: SIMD -vectorized twiddles-on-the-fly computation:
	const double sign[2] = {+1.0,-1.0};
	float ff[16];
	double theta[16], gamma[16], scos[2];
	struct complex twiddle[16];
	int ii[16],io[16],iodd[16],iq[16],is0[16],is1[16],jj[16];
	#if 0	// Quadrants-based scheme:
		const double pi2 = (double)1.57079632679489661922, pi_mults[4] = {0,pi2,2*pi2,3*pi2}, twopin = 8*pi2/n, sign[2] = {+1.0,-1.0};
		int i15 = i*15, iq = i15/(n>>3),is0,is1;	// iq = 'index of quadrant'
		double theta = i15*twopin, gamma = theta - pi_mults[iq], scos[2] = {cos(gamma),sin(gamma)};
		is0 = ((unsigned int)iq - 1) < 2; is1 = (iq > 1);
		j = IS_ODD(iq);
	#else	// Octants-based scheme:
		//uint32 idiv = (n>>4), imult = (uint32)((0x0000000100000000ull+idiv-1)/(n>>4));	// Divisor idiv = (n>>4); imult = (2^32+idiv-1)/idiv
		float fndiv16 = (float)1.0/(float)(n>>4), *fndiv16_ptr = &fndiv16;
		// Loop to test various fast alternatives to j/(n>>4) for every j < n/2:
		for(j = 0; j < (n>>1); j++) {
			// This fails for e.g. j = 393205 and (n>>4) = 393216:
			//	ASSERT(__MULH32(j,imult) == j/(n>>4), "umulh32(j,imult) != j/(n>>4)");
			ASSERT((int)((float)j*fndiv16) == j/(n>>4), "(float)j*fndiv16 != j/(n>>4)");
		}
	i = 163397;
		const double pi4_dbl = (double)0.78539816339744830961, twopin_dbl = 16*pi4_dbl/n;
		double pi_mults[9]; for(j = 0; j < 9; j++) { pi_mults[j] = j*pi4_dbl; }
		for(j = 0; j < 16; j++) {
			ii[j] = i*j; io[j] = ii[j]/(n>>4); iodd[j] = io[j]&1;	// io = 'index of octant'
			ff[j] = (float)ii[j]*fndiv16;
			theta[j] = ii[j]*twopin_dbl; gamma[j] = sign[iodd[j]]*(theta[j] - pi_mults[io[j]+iodd[j]]);
			scos[0] = cos(gamma[j]); scos[1] = sin(gamma[j]);
			is0[j] = ((unsigned int)io[j] - 2) < 4; is1[j] = (io[j] > 3);
			jj[j] = IS_ODD((io[j]+1)>>1);
			ASSERT((int)ff[j] == io[j], "ff != io error!");
			twiddle[j].re = sign[is0[j]]*scos[jj[j]]; twiddle[j].im = sign[is1[j]]*scos[jj[j]^1];
		}
	#endif
	#if 1
	/*
		PMULLD is slow! Latency = 6, and a PMULLD can start only every other cycle.
		Need a faster way to compute [0-15]*iroot ... try some things involving fast permute/shift/adds.
		Op		Form	Latency	Thruput[how many can start per cycle]
		-----	----	----	----
		movd	r32,xmm	2		3	<*** cheaper to do 2-3 redundant movd's instead of movd followed by xmm-copy
		movddup	xmm,xmm	1		1	<*** formally for copying double in low half of src to both halves of dest, likely data-type-agnostic
		pshufd	2.xmm,i	1		2
		paddd	xmm,xmm	1		2	<*** prefer add to << 1 for doubling ***
		pslld	xmm,i	1		1
	We then need to do a batch of CVTDQ2PS vector-int32-to-float conversions, which need 5 cycles and can execute 1/cycle...that leads to
	this still-much-too-slow instruction sequence:
		movl	[__i],eax
		movq	[__sm_ptr],rdx
		leal	(eax,2),ebx			leal	(eax,4),ecx		// 2i,4i
		movd rax,xmm0  movd rbx,xmm1  movd rcx,xmm2 // [1,2,4]*i in low 32 bits of xmm0,1,2 ... movd latency 2 but 3/cycle
		pshufd	17,xmm0,xmm0	// [1,0,1,0]*i [MSW of vec @left], order-byte bit pattern [00 gives 1*i, 01 gives 0*i] = 00.01.00.01 = 17
		pshufd	 5,xmm1,xmm1	// [2,2,0,0]*i, order-byte bit pattern [00 gives 2*i, 01 gives 0*i] = 00.00.01.01 = 5
		pshufd	 0,xmm2,xmm2	// [4,4,4,4]*i, order-byte bit pattern [00 gives 4*i, 01 gives 0*i] = 00.00.00.00 = 0
		paddd		xmm1,xmm0	// [3,2,1,0]*i = [1,0,1,0]*i + [2,2,0,0]*i; result overwrites [1,0,1,0]*i
		movaps	xmm2,xmm1		// copy of [4,4,4,4]*i overwrites [2,2,0,0]*i
		movaps	xmm2,xmm3		// another copy of [4,4,4,4]*i
		paddd		xmm2,xmm2	// [8,8,8,8]*i = [4,4,4,4]*i + [4,4,4,4]*i
		paddd		xmm0,xmm1	// [7,6,5,4]*i = [3,2,1,0]*i + [4,4,4,4]*i; result overwrites copy of [4,4,4,4]*i
		paddd		xmm3,xmm3	// [8,8,8,8]*i = [4,4,4,4]*i + [4,4,4,4]*i
		paddd		xmm0,xmm2	// [b,a,9,8]*i = [3,2,1,0]*i + [8,8,8,8]*i
		paddd		xmm1,xmm3	// [f,e,d,c]*i = [7,6,5,4]*i + [8,8,8,8]*i
		// Vector analog of the above octants-based scheme, with short-length polynomial approximation to cos(gamma),sin(gamma):
		movss	[__fndiv16],xmm4
		pshufd	0,xmm4,xmm4	// Broadcast low 32 bits of xmm4 to all 4 slots of xmm4 ***Want PSHUFD, *not* SHUFPD, even though float data ***
		// There is no vector-integer DIV, so do the vector io = [ 0-15]*iroot/(n>>4) via SP-floats mul-by-reciprocal:
		cvtdq2ps	xmm0,xmm0	[similar for xmm1-3]
		mulps	xmm4,xmm0			"
		cvtps2dq	xmm0,xmm0		"

	Wait! The fundamental problem is that vector-int MUL is horribly slow on older x86_64 releases, to the extent it was supported at all.
	But vector-float|double MUL is fast - we can use vector-floats instead. Now e.g. 15*iroot may be larger that a float 24-bit mantissa
	can hold, but iroot will always have a large power-of-2 component, so that incurs no precision loss in the present context.
	Similarly with the fndiv16 = 1/(n>>4) term we multiply the [0-15]*iroot vector with ... the only potential precision loss there
	is any associated with the odd component of n, and since the result, which will get (int)-cast into the octant index, is < 8, it
	seems unlikely that said index would fail to match the exactly-computed one, except possibly in very rare (Q: how rare?) cases.
	So, we init a static vec-float = [15.0,14.0,...,1.0,0]
	*/
	__asm__ volatile (\
	"movl	%[__i],%%ebx		\n\t"\
	"movq	%[__fndiv16],%%rcx	\n\t"\
	"movq	%[__sc_ptr],%%rax	\n\t"\
	"movd	%%rbx,%%xmm4		\n\t"\
	"movss	(%%rcx),%%xmm5		\n\t"/* [0,0,0,1]*fndiv16 */\
		"pshufd	$0,%%xmm4,%%xmm4	\n\t"/* [1,1,1,1]*i */\
		"pshufd	$0,%%xmm5,%%xmm5	\n\t"/* Broadcast fndiv16 to all 4 (float) slots of xmm5 ***Want PSHUFD, *not* SHUFPD, even though float data ***/\
		"cvtdq2pd	%%xmm4,%%xmm7	\n\t"/* Bottom 2 copies of i get converted to doubles */\
		"cvtdq2ps	%%xmm4,%%xmm4	\n\t"/* xmm5 already has float-format data, no type conversion needed */\
	/* Here is the local-store layout for the auxiliary data needed by this macro:
		sign_mask = sc_ptr + 0x48;	// Start v20-added data after v19-alloc block
		pi4       = sc_ptr + 0x4c;
		c0thru15  = sc_ptr + 0x4d;
		// Starting at c0thru15,
		// first 4 = [1 << (6-L2_SZ_VD)] slots = sc_ptr + 0x[4a-4d] hold (float)0-15
		// next  8 = [1 << (7-L2_SZ_VD)] slots = sc_ptr + 0x[4e-56] hold (double)0-15
	*/\
	"movaps	0x4d0(%%rax),%%xmm0	\n\t"/* (float)[3,2,1,0] */\
	"movaps	0x4e0(%%rax),%%xmm1	\n\t"/* (float)[7,6,5,4] */\
	"movaps	0x4f0(%%rax),%%xmm2	\n\t"/* (float)[b,a,9,8] */\
	"movaps	0x500(%%rax),%%xmm3	\n\t"/* (float)[f,e,d,c] */\
		"mulps	%%xmm5,%%xmm4	\n\t"/* i*fndiv16, truncate-to-int gives io = i/(n>>4) */\
	"movaps	0x510(%%rax),%%xmm8	\n\t"/* (double)[1,0]*twopin */\
	"movaps	0x520(%%rax),%%xmm9	\n\t"/* (double)[3,2]*twopin */\
	"movaps	0x530(%%rax),%%xmm10\n\t"/* (double)[5,4]*twopin */\
	"movaps	0x540(%%rax),%%xmm11\n\t"/* (double)[7,6]*twopin */\
		"mulps	%%xmm4,%%xmm0	\n\t"/* (float)[3,2,1,0]*i*fndiv16 */\
		"mulps	%%xmm4,%%xmm1	\n\t"/* (float)[7,6,5,4]*i*fndiv16 */\
		"mulps	%%xmm4,%%xmm2	\n\t"/* (float)[b,a,9,8]*i*fndiv16 */\
		"mulps	%%xmm4,%%xmm3	\n\t"/* (float)[f,e,d,c]*i*fndiv16 */\
		"mulpd	%%xmm7,%%xmm8	\n\t"/* (double)[1,0]*twopin*i */\
		"mulpd	%%xmm7,%%xmm9	\n\t"/* (double)[3,2]*twopin*i */\
		"mulpd	%%xmm7,%%xmm10	\n\t"/* (double)[5,4]*twopin*i */\
		"mulpd	%%xmm7,%%xmm11	\n\t"/* (double)[7,6]*twopin*i */\
		/* (int)((float)i*fndiv16) gives io = i/(n>>4), the vector octant-index: */\
		"cvttps2dq	%%xmm0,%%xmm0	\n\t"\
		"cvttps2dq	%%xmm1,%%xmm1	\n\t"\
		"cvttps2dq	%%xmm2,%%xmm2	\n\t"\
		"cvttps2dq	%%xmm3,%%xmm3	\n\t"\
/*** Harmonize xmm-indexing between above and below ***/\
		/* Convert the 16 dword io-values to concatenated 16-byte vec: */\
		/* Viewed as 16-byte-vecs, xmm0=[-,-,-,3,-,-,-,2,-,-,-,1,-,-,-,0] */\
		"pslldq	$1,%%xmm1		\n\t"/* [-,-,7,-,-,-,6,-,-,-,5,-,-,-,4,-] */\
		"pslldq	$2,%%xmm2		\n\t"/* [-,b,-,-,-,a,-,-,-,9,-,-,-,8,-,-] */\
		"pslldq	$3,%%xmm3		\n\t"/* [f,-,-,-,e,-,-,-,d,-,-,-,c,-,-,-] */\
		"paddd %%xmm1,%%xmm0	\n\t"\
		"paddd %%xmm3,%%xmm2	\n\t"\
		"paddd %%xmm2,%%xmm0	\n\t"/* xmm0 = [f,b,7,3,e,a,6,2,d,9,5,1,c,8,4,0] */\
	/***************************************************************************************************
	would really like io-values packed 4-bits each into a int64 for later manipulations ... can we do this
	and still preserve the 4-bytes-apart-ness of sequential index pairs for the below cast-to-double?
	Yes we can: Take above vec-of-bytes and fold-with-shift upper half into lower:
		[f,b,7,3,e,a,6,2]<<4 + [d,9,5,1,c,8,4,0]
	then save copy of resulting 64-bit low-half to mem and use working copy for cast-to-double.
	***************************************************************************************************/\
		"movaps	%%xmm0,%%xmm1	\n\t"/* copy in xmm1 */
		"psrldq	$8,%%xmm1		\n\t"/* [-,-,-,-,-,-,-,-,f,b,7,3,e,a,6,2] */
		"pslld	$4,%%xmm1		\n\t"/* Upper halves of all bytes = 0, so can use any desired granularity for the <<= 4 */\
		"paddd %%xmm1,%%xmm0	\n\t"/* xmm0: lo64 = packed vector of nibbles [f,d,b,9,7,5,3,1,e,c,a,8,6,4,2,0] *** hi64 still = [f,b,7,3,e,a,6,2] ***/\
	/* Due to ambiguity re. MOVQ in GCC inline-asm can't copy to a GPR, so save copy of
	64-bit packed-nibbles io[]-vec to memloc just below sign_mask: */\
	"movsd	%%xmm0,0x478(%%rax)	\n\t"\
	"movaps	0x460(%%rax),%%xmm1	\n\t"/* next-lower slot holds 0x1111111111111111 in lo64, 0 in hi64: */\
	"movaps	0x550(%%rax),%%xmm12	\n\t"/* (double)[9,8]*twopin */\
	"movaps	0x560(%%rax),%%xmm13	\n\t"/* (double)[b,a]*twopin */\
	"movaps	0x570(%%rax),%%xmm14	\n\t"/* (double)[d,c]*twopin */\
	"movaps	0x580(%%rax),%%xmm15	\n\t"/* (double)[f,e]*twopin */\
		"mulpd	%%xmm7,%%xmm12	\n\t"/* (double)[9,8]*twopin*i */\
		"mulpd	%%xmm7,%%xmm13	\n\t"/* (double)[b,a]*twopin*i */\
		"mulpd	%%xmm7,%%xmm14	\n\t"/* (double)[d,c]*twopin*i */\
		"mulpd	%%xmm7,%%xmm15	\n\t"/* (double)[f,e]*twopin*i */\
	"movaps	0x450(%%rax),%%xmm7	\n\t"/* next-lower slot holds 0x0000000f0000000f in lo64, 0 in hi64: */\
		/* AND io (in lo64 of xmm0) with a low-half vector of 1-nibbled to effect iodd = io & 1: */\
		"pand	%%xmm0,%%xmm1	\n\t"/* lo64 = iodd in vec-of-nibbles */\
		"paddb	%%xmm1,%%xmm0	\n\t"/* io += iodd: Could actually use any of PADDB/PADDW/PADDD/PADDQ here since no chance of overflow */\
		/* write results: */\
/**** Need to duplicate each of the 16 scalar-double cos outputs before writing twinned copies to memory ****/\
		:					/* outputs: none */\
		: [__sc_ptr] "m" (sc_ptr),	/* All inputs from memory addresses here */\
		  [__i] "m" (i),\
		  [__n] "m" (n),\
		  [__fndiv16] "m" (fndiv16_ptr)\
		: "cc","memory","rax","rbx","rcx","r8","r9","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);
	#endif
  #endif	// #endif(0) - experimental code!

	#ifdef REFACTOR_4DFT_3TWIDDLE

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c1=re0*re1-im0*im1;	s1=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c2=re0*re1-im0*im1;	s2=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c3=re0*re1-im0*im1;	s3=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c4=re0*re1-im0*im1;	s4=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c5=re0*re1-im0*im1;	s5=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c6=re0*re1-im0*im1;	s6=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c7=re0*re1-im0*im1;	s7=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c8=re0*re1-im0*im1;	s8=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c9=re0*re1-im0*im1;	s9=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c10=re0*re1-im0*im1;	s10=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c11=re0*re1-im0*im1;	s11=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c12=re0*re1-im0*im1;	s12=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c13=re0*re1-im0*im1;	s13=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c14=re0*re1-im0*im1;	s14=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c15=re0*re1-im0*im1;	s15=re0*im1+im0*re1;

		// Cf. util.c:test_radix16_dft() for details on these combos:
		c5 =  (c1-s1)*ISRT2;	s5 = (c1+s1)*ISRT2;
		c6 = -s2;				s6 = c2;
		c7 = -(c3+s3)*ISRT2;	s7 = (c3-s3)*ISRT2;

		c9 = c1*c-s1*s;			s9 = c1*s+s1*c;
		c10= (c2-s2)*ISRT2;		s10= (c2+s2)*ISRT2;
		c11= c3*s-s3*c;			s11= c3*c+s3*s;

		c13= c1*s-s1*c;			s13= c1*c+s1*s;
		c14= -(c2+s2)*ISRT2;	s14= (c2-s2)*ISRT2;
		c15= -c3*c+s3*s;		s15= -(c3*s+s3*c);

	  #ifdef USE_SSE2
		/* Sincos data stored in terms of the following 5 contiguous-data triplets:
			c4,s4, c8,s8, cC,sC
			c1,s1, c2,s2, c3,s3
			c5,s5, c6,s6, c7,s7
			c9,s9, cA,sA, cB,sB
			cD,sD, cE,sE, cF,sF .
		Note that due to my layout of the SSE2_RADIX_04_DIF_3TWIDDLE_X2-macro arglist,
		we need to swap the order of the first 2 sincos-pairs of each triplet:
		*/
		vec_dbl *c_tmp = cc0, *s_tmp = c_tmp+1;	/* c0,s0 */
		VEC_DBL_INIT(c_tmp, c8 );	VEC_DBL_INIT(s_tmp, s8 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c4 );	VEC_DBL_INIT(s_tmp, s4 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c12);	VEC_DBL_INIT(s_tmp, s12);	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c2 );	VEC_DBL_INIT(s_tmp, s2 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c1 );	VEC_DBL_INIT(s_tmp, s1 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c3 );	VEC_DBL_INIT(s_tmp, s3 );	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c6 );	VEC_DBL_INIT(s_tmp, s6 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c5 );	VEC_DBL_INIT(s_tmp, s5 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c7 );	VEC_DBL_INIT(s_tmp, s7 );	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c10);	VEC_DBL_INIT(s_tmp, s10);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c9 );	VEC_DBL_INIT(s_tmp, s9 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c11);	VEC_DBL_INIT(s_tmp, s11);	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c14);	VEC_DBL_INIT(s_tmp, s14);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c13);	VEC_DBL_INIT(s_tmp, s13);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c15);	VEC_DBL_INIT(s_tmp, s15);	c_tmp+=2; s_tmp+=2;
	  #endif

	#elif defined(USE_AVX2)
		// In AVX2/FMA mode, since we need to replace most of the raw sincos data with derived ones,
		// simply place one copy of each computed double in a double-sized slot of the local memory.
		// We will be using AVX2/FMA-based Newtonian iterative inversion on the 16 doubles whose
		// multiplicative inverse is needed (the real part of the basic root of unity c and of the
		// 15 complex twiddles, c1-15), so store those in packed form in 4 AVX-register-sized
		// contiguous memory locations, and the others in a separate chunk of memory. After the
		// vector-iterative inversion we'll need to combine the 2 sets of data and place (in suitable
		// vector-register-sized broadcast form) into their final SIMD-suitable memory slots.

		add0 = (double *)cc0;		// add0 points to 16 cos-data-to-be-inverted; Need a double-ptr on lhs here
		double *add1 = add0 + 16;	// add1 points to block of memory temporarily used to store the corresponding sine data
		double *add2 = add0 + 32;	// add2 points to block of memory temporarily used to store the 11 [0-padded to 12]
									//	cosine data which need to be divided by other cosines (i.e. multiplied by inverses)
		/* The add2-addressed cosine ratios are arranged in 3 YMM-register/memory-sized slots like so;
		  once we have filled 4 YYMs with inverses 1/[c3,c1-15] and used those to get the 16 tangents (1st set = 1/c3
		  and discarded) we will do as described in the right column to set up for the cosine-ratios computation:

			double __c31 = __c3/__c1;
			double __c51 = __c5/__c1;
			double __c62 = __c6/__c2;
			[0 pad]						shuffle YMM with 1/[c3,c1,c2,c3] to get 1/[c1,c1,c2,c3], then *= [c3,c5,c6,0]

			double __c73 = __c7/__c3;
			double __c91 = __c9/__c1;
			double __cA2 = __cA/__c2;
			double __cB3 = __cB/__c3;	initialize YMM with 1/[c3,c1,c2,c3], then *= [c7,c9,cA,cB]

			double __cC4 = __cC/__c4;
			double __cD5 = __cD/__c5;
			double __cE6 = __cE/__c6;
			double __cF7 = __cF/__c7;	Multiply YMM with 1/[c4-7] *= [cC-F]
		*/
		*add0++ = 0.0;	// Since tan0 defined as const, use this pair of double slots to hold 1/c3 (via c3,1 on input, then invert c3 and multiply
		*add1++ = 1.0;	// them together), which extra 1/c3 copy saves some really awkward permuting, at least in terms of the idiotic x86 ISA.

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c1, for inversion
		*add1++ = it;	// s1  slot will hold __r1 = s1 /c1

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c2, for inversion
		*add1++ = it;	// s2  slot will hold __r2 = s2 /c2

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*(add0-3) = rt;	// c3, for inversion ...
		*add0++   = rt;	// place extra copy in 0-slot as described above - put on separate line to avoid ambiguity of *(add0-3) = *add0++ = rt
		*add1++ = it;	// s3  slot will hold __r3 = s3 /c3
		*add2++ = rt;	// c3, will get multiplied by 1/c1 to yield __c31
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c4, for inversion
		*add1++ = it;	// s4  slot will hold __r4 = s4 /c4

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c5, for inversion
		*add1++ = it;	// s5  slot will hold __r5 = s5 /c5
		*add2++ = rt;	// c5, will get multiplied by 1/c1 to yield __c51
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c6, for inversion
		*add1++ = it;	// s6  slot will hold __r6 = s6 /c6
		*add2++ = rt;	// c6, will get multiplied by 1/c2 to yield __c62
		*add2++ = 0.0;	// 0-pad will get multiplied by 1/c3 term, remains 0-pad.
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c7, for inversion
		*add1++ = it;	// s7  slot will hold __r7 = s7 /c7
		*add2++ = rt;	// c7, will get multiplied by 1/c3 to yield __c73
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c8, for inversion
		*add1++ = it;	// s8  slot will hold __r8 = s8 /c8

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c9, for inversion
		*add1++ = it;	// s9  slot will hold __r9 = s9 /c9
		*add2++ = rt;	// c9, will get multiplied by 1/c1 to yield __c91
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c10, for inversion
		*add1++ = it;	// s10 slot will hold __rA = s10/c10
		*add2++ = rt;	// cA, will get multiplied by 1/c2 to yield __cA2
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c11, for inversion
		*add1++ = it;	// s11 slot will hold __rB = s11/c11
		*add2++ = rt;	// cB, will get multiplied by 1/c3 to yield __cB3
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c12, for inversion
		*add1++ = it;	// s12 slot will hold __rC = s12/c12
		*add2++ = rt;	// cC, will get multiplied by 1/c4 to yield __cC4
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c13, for inversion
		*add1++ = it;	// s13 slot will hold __rD = s13/c13
		*add2++ = rt;	// cD, will get multiplied by 1/c5 to yield __cD5
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c14, for inversion
		*add1++ = it;	// s14 slot will hold __rE = s14/c14
		*add2++ = rt;	// cE, will get multiplied by 1/c6 to yield __cE6

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c15, for inversion
		*add1++ = it;	// s15 slot will hold __rF = s15/c15
		*add2++ = rt;	// cF, will get multiplied by 1/c7 to yield __cF7

		// This places us at add0 == c8 and add1 = c12.
		ASSERT(add0 == (double *)cc0+16 && add1 == (double *)cc0+32 && add2 == (double *)cc0+44, "add0,1,2 checksum failed in AVX2 sincos inits!");
		/*
		At this point, the 11 ymm-sized [32-byte] chunks starting at &cc0 contain the following scalar-double data:

		0:	c3,c1-3
		1:	c4-7
		2:	c8-11
		3:	c12-c15
		4:	1.0,s1-3
		5:	s4-7
		6:	s8-11
		7:	s12-s15
		8:	c3,5,6,[0-pad]
		9:	c7,9-B
		A:	cC-F
		*/

		// Now send the cosine terms to the inversion routine, which also does the combine-and-populate-SIMD-slots step.
	  #ifdef DFT_V2	// Toggle between 2 versions

		RADIX16_COMPUTE_FMA_SINCOS_DIF_2(cc0,two);	//*** Note 'two' contains 1.0x4 in this mode! ***
		add0 = (double *)cc0;

		/* 44 scalar-double data starting at add0 = cc0 laid out as below - * denote the 12 data unused by the V2 DIF DFT:

		a[0x00,0x01,0x02,0x03]: add0 + 0x[  0,  8, 10, 18]: [c3*,c1 ,c2 ,c3*] Cosines:
		a[0x04,0x05,0x06,0x07]: add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5*,c6*,c7*]
		a[0x08,0x09,0x0a,0x0b]: add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9*,cA*,cB*]
		a[0x0c,0x0d,0x0e,0x0f]: add0 + 0x[ 60, 68, 70, 78]: [cC*,cD*,cE*,cF*]
		a[0x10,0x11,0x12,0x13]: add0 + 0x[ 80, 88, 90, 98]: [-- ,r1 ,r2 ,r3 ] Tangents:
		a[0x14,0x15,0x16,0x17]: add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ]
		a[0x18,0x19,0x1a,0x1b]: add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ]
		a[0x1c,0x1d,0x1e,0x1f]: add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ]
		a[0x20,0x21,0x22,0x23]: add0 + 0x[100,108,110,118]: [c31,c51,c62,  0] Cosine ratios:
		a[0x24,0x25,0x26,0x27]: add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3]
		a[0x28,0x29,0x2a,0x2b]: add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7]

		V2 DIF needs 35 doubles - subtracting the 2 unused/0 slots and the 12 unused cosine-data
		in the above leaves just 30. We get to 35 by applying the following data-munges, expressed
		in terms of the hex double-array indices of the targeted datum:
		*/
		add0[0x00] = 1.0;				// 1. replace 0x00 with 1.0
		add0[0x03] = add0[0x01]*c;		// 2. replace 0x03 with c1*c
		add0[0x05] = add0[0x01]*ISRT2;	// 3. replace 0x05,0x06 with [c1,c2]*ISRT2
		add0[0x06] = add0[0x02]*ISRT2;
		add0[0x10] = tan;				// 4. replace 0x10 with tan = s/c
		/*
		Yielding the following layout-of-scalar-doubles, with dash-marked slots unused in V2 DIF
		and appended 'I' meaning *= ISRT2; we can see the expected 9 unused slots:

		a[0x00,0x01,0x02,0x03]: add0 + 0x[  0,  8, 10, 18]: [1.0,c1 ,c2 ,c1*c]
		a[0x04,0x05,0x06,0x07]: add0 + 0x[ 20, 28, 30, 38]: [c4 ,c1I,c2I,---]
		a[0x08,0x09,0x0a,0x0b]: add0 + 0x[ 40, 48, 50, 58]: [c8 ,---,---,---]
		a[0x0c,0x0d,0x0e,0x0f]: add0 + 0x[ 60, 68, 70, 78]: [---,---,---,---]
		a[0x10,0x11,0x12,0x13]: add0 + 0x[ 80, 88, 90, 98]: [s/c,r1 ,r2 ,r3 ] Tangents:
		a[0x14,0x15,0x16,0x17]: add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ]
		a[0x18,0x19,0x1a,0x1b]: add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ]
		a[0x1c,0x1d,0x1e,0x1f]: add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ]
		a[0x20,0x21,0x22,0x23]: add0 + 0x[100,108,110,118]: [c31,c51,c62,---] Cosine ratios:
		a[0x24,0x25,0x26,0x27]: add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3]
		a[0x28,0x29,0x2a,0x2b]: add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7]
		*/
	  #else	// Make DFT_V1 the default here:

		const double *cd_ptr0 = &c, *cd_ptr1 = &tan;	// GCC/Clang don't allow cd_address-taking inlined in arglist of macros, so do it here
		RADIX16_COMPUTE_FMA_SINCOS_DIF(cc0,two,cd_ptr0,cd_ptr1);

	  #endif

	#else	// AVX2 = False:

	  #ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-15] are offset w.r.to the thread-local ptr pair as
		(cc0,ss0) + 0x[2,12,a,1a, 6,16,e,1e, 4,14,c,1c, 8,18,10,20], i.e. BR-ordered [0,8,4,C, 2,A,6,E, 1,9,5,D, 3,B,7,F];:
		*/
		c_tmp = cc0 + 0x02; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_modq(mt0[k1].d0,mt0[k1].d1, mt1[k2].d0,mt1[k2].d1, &rm,&im);
		// Premultiply results by 8 so ensuing twiddle-MULs can use the streamlined CMUL_MODQ8 routine:
		a1 = qreduce_full(rm)<<3;	b1 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x12; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c1 =rt;		s1 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		csqr_modq(a1,b1, &rm,&im);	// Since MODMULs are roundoff-error-free, cheaper to use a LOACC-style sequence
						// of successive squarings to get E^2,4,8, a CMUL to get E^13, then CMUL_CONJ_MODQs to get rest.
		a2 = qreduce_full(rm)<<3;	b2 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0xa; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c2 =rt;		s2 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x1a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c3 =rt;		s3 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		csqr_modq(a2,b2, &rm,&im);
		a4 = qreduce_full(rm)<<3;	b4 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x6; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c4 =rt;		s4 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a1,b1, a4,b4, &a3,&b3, &a5,&b5);
		a3 = qreduce_full(a3)<<3;	b3 = qreduce_full(b3)<<3;
		a5 = qreduce_full(a5)<<3;	b5 = qreduce_full(b5)<<3;
	  #endif
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c5 =rt;		s5 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0xe; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c6 =rt;		s6 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x1e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c7 =rt;		s7 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		csqr_modq(a4,b4, &rm,&im);
		a8 = qreduce_full(rm)<<3;	b8 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x4; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c8 =rt;		s8 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a1,b1, a8,b8, &a7,&b7, &a9,&b9);
		a7 = qreduce_full(a7)<<3;	b7 = qreduce_full(b7)<<3;
		a9 = qreduce_full(a9)<<3;	b9 = qreduce_full(b9)<<3;
	  #endif
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c9 =rt;		s9 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a2,b2, a8,b8, &a6,&b6, &a10,&b10);
		a6  = qreduce_full(a6 )<<3;	b6  = qreduce_full(b6 )<<3;
		a10 = qreduce_full(a10)<<3;	b10 = qreduce_full(b10)<<3;
	  #endif
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0xc; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c10=rt;		s10=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x1c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c11=rt;		s11=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x8; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c12=rt;		s12=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_modq(mt0[k1].d0,mt0[k1].d1, mt1[k2].d0,mt1[k2].d1, &rm,&im);
		a13 = qreduce_full(rm)<<3;	b13 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x18; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c13=rt;		s13=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a1,b1, a13,b13, &a12,&b12, &a14,&b14);
		a12 = qreduce_full(a12)<<3;	b12 = qreduce_full(b12)<<3;
		a14 = qreduce_full(a14)<<3;	b14 = qreduce_full(b14)<<3;
	  #endif
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x10; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c14=rt;		s14=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a2,b2, a13,b13, &a11,&b11, &a15,&b15);
		a11 = qreduce_full(a11)<<3;	b11 = qreduce_full(b11)<<3;
		a15 = qreduce_full(a15)<<3;	b15 = qreduce_full(b15)<<3;
	  #endif
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x20; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c15=rt;		s15=it;
	  #endif

	#endif	// USE_AVX2?

#else
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c1=t1*t3-t2*t4;	s1=t1*t4+t2*t3;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += i;				/* 4*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c2=t1*t3-t2*t4;	s2=t1*t4+t2*t3;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += i;				/* 8*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c4=t1*t3-t2*t4;	s4=t1*t4+t2*t3;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += (iroot << 2)+iroot;		/* 13*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c8=t1*t3-t2*t4;	s8=t1*t4+t2*t3;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c13=t1*t3-t2*t4;	s13=t1*t4+t2*t3;

		/* c3,5 */
		t1=c1*c4; t2=c1*s4; t3=s1*c4; t4=s1*s4;
		c3=t1+t4; s3=t2-t3; c5=t1-t4; s5=t2+t3;

		/* c6,7,9,10 */
		t1=c1*c8; t2=c1*s8; t3=s1*c8; t4=s1*s8;
		c7=t1+t4; s7=t2-t3; c9=t1-t4; s9=t2+t3;

		t1=c2*c8; t2=c2*s8; t3=s2*c8; t4=s2*s8;
		c6=t1+t4; s6=t2-t3; c10=t1-t4; s10=t2+t3;

		/* c11,12,14,15 */
		t1=c1*c13; t2=c1*s13; t3=s1*c13; t4=s1*s13;
		c12=t1+t4; s12=t2-t3; c14=t1-t4; s14=t2+t3;

		t1=c2*c13; t2=c2*s13; t3=s2*c13; t4=s2*s13;
		c11=t1+t4; s11=t2-t3; c15=t1-t4; s15=t2+t3;
#endif

#ifndef USE_SSE2
  #ifdef USE_SCALAR_DFT_MACRO	// Must define - or not - @compile time
	// FMA-based version replaces sine terms with tangents:
	s1  /= c1;
	s2  /= c2;
	s3  /= c3;
	s4  /= c4;
	s5  /= c5;
	s6  /= c6;
	s7  /= c7;
	s8  /= c8;
	s9  /= c9;
	s10 /= c10;
	s11 /= c11;
	s12 /= c12;
	s13 /= c13;
	s14 /= c14;
	s15 /= c15;
	/* Cosine terms defined like so: */
	double c1_c = c1*c;	/* In FMA code, this will take the place of c */
	double c1i2 = c1*ISRT2;
	double c2i2 = c2*ISRT2;
	c9  /= c1;
	c10 /= c2;
	c11 /= c3;
	c12 /= c4;
	c13 /= c5;
	c14 /= c6;
	c15 /= c7;
	// Delay these, since need them unmodified as divisors in the above defs:
	c7  /= c3;
	c3  /= c1;
	c5  /= c1;
	c6  /= c2;
  #endif
#endif

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 4);

	/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
	Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
	but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
	(But only ever need to explicitly do this in debug mode).
	*/
	  for(j = jlo; j < jhi; j += stride)
	  {
		j1 = j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	#ifdef USE_SSE2

		pfetch_addr = p1*((j >> l2_stride) & 0x3);	// cycle prefetch-offset-address among p0,1,2,3
			// These get added to the base addresses p0,4,8,12 in the DFT macros, thus every 4 loop
			// executions we have covered prefetches from [current address] + [pfetch distance] + p0,1,2,...15 .

		add0 = &a[j1];

	  #ifdef USE_AVX2

	   #ifdef REFACTOR_4DFT_3TWIDDLE
		// Since pvsly had no non-suffixed 'SSE2_RADIX16_DIF_TWIDDLE' macro in AVX2 mode, no need for _V2 suffix here:
		SSE2_RADIX16_DIF_TWIDDLE(add0,p1,p2,p3,p4,p8,p12,r1,two,cc0,pfetch_addr,pfetch_dist);

	   #elif defined(DFT_V2)	// Toggle between 2 versions

		SSE2_RADIX16_DIF_TWIDDLE_2(add0,p1,p2,p3,p4,p8,p12,r1,  cc0,pfetch_addr,pfetch_dist);

	   #else	// Make DFT_V1 the default:

		SSE2_RADIX16_DIF_TWIDDLE_1(add0,p1,p2,p3,p4,p8,p12,r1,isrt2,pfetch_addr,pfetch_dist);

	   #endif

	  // Play with prefetch in 64-bit AVX/SSE2 versions - once find reliable pfetch scheme, propagate to 32-bit SSE2 and AVX2 macros:
	  #else

	   #ifdef REFACTOR_4DFT_3TWIDDLE

		SSE2_RADIX16_DIF_TWIDDLE_V2(add0,p1,p2,p3,p4,p8,p12,r1,two,cc0,pfetch_addr,pfetch_dist);

	   #else

		SSE2_RADIX16_DIF_TWIDDLE(add0,p1,p2,p3,p4,p8,p12,r1,isrt2,pfetch_addr,pfetch_dist);

	   #endif

	  #endif

	#else	// USE_SSE2 = False:

	  #ifdef USE_SCALAR_DFT_MACRO	// Must define - or not - @compile time

		// Test FMA-based DIF macro:
		RADIX_16_DIF_FMA(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p4+p1 ],a[j2+p4+p1 ],a[j1+p4+p2 ],a[j2+p4+p2 ],a[j1+p4+p3 ],a[j2+p4+p3 ],a[j1+p8 ],a[j2+p8 ],a[j1+p8+p1 ],a[j2+p8+p1 ],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+p12],a[j2+p12],a[j1+p12+p1],a[j2+p12+p1],a[j1+p12+p2],a[j2+p12+p2],a[j1+p12+p3],a[j2+p12+p3]
						,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p4+p1 ],a[j2+p4+p1 ],a[j1+p4+p2 ],a[j2+p4+p2 ],a[j1+p4+p3 ],a[j2+p4+p3 ],a[j1+p8 ],a[j2+p8 ],a[j1+p8+p1 ],a[j2+p8+p1 ],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+p12],a[j2+p12],a[j1+p12+p1],a[j2+p12+p1],a[j1+p12+p2],a[j2+p12+p2],a[j1+p12+p3],a[j2+p12+p3]
						,c1,s1,c2,s2,c3,s3,c4,s4,c5,s5,c6,s6,c7,s7,c8,s8,c9,s9,c10,s10,c11,s11,c12,s12,c13,s13,c14,s14,c15,s15
						,c1_c,tan,c1i2,c2i2)

	  #else		// USE_SCALAR_DFT_MACRO = False

	#ifdef USE_FGT61

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/
	  #ifdef PFETCH_AGGRESSIVE
		addp = &a[j1];
	  #elif PFETCH
		prefetch_offset = ((j >> 1) & 3)*p4 + 4;	/* Cycle among p0, p4, p8 and p12. */
	  #endif

	/*...Block 1: */
		jt = j1;		jp = j2;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t1 =a[jt	];						t2 =a[jp	];						m1 =b[jt	];				m2 =b[jp	];			//   0, b (treat as 0,4q for simplicity's sake)
		rt =a[jt+p8 ]*c8 -a[jp+p8 ]*s8;		it =a[jp+p8 ]*c8 +a[jt+p8 ]*s8;		cmul_modq8(b[jt+p8 ],b[jp+p8 ], a8,b8, &rm,&im);	//   0,4q
		t3 =t1 -rt;							t1 =t1 +rt;							m3 =qreduce(m1 -rm+q4);		m4 =qreduce(m2 -im+q4);	// -4q,4q	=> 0,b
		t4 =t2 -it;							t2 =t2 +it;							m1 =qreduce(m1 +rm   );		m2 =qreduce(m2 +im   );	//   0,8q	=> 0,b (rest similar)
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t5 =a[jt+p4 ]*c4 -a[jp+p4 ]*s4 ;	t6 =a[jp+p4 ]*c4 +a[jt+p4 ]*s4;		cmul_modq8(b[jt+p4 ],b[jp+p4 ],a4 ,b4 ,&m5 ,&m6 );
		rt =a[jt+p12]*c12-a[jp+p12]*s12;	it =a[jp+p12]*c12+a[jt+p12]*s12;	cmul_modq8(b[jt+p12],b[jp+p12],a12,b12,&rm ,&im );
		t7 =t5 -rt;							t5 =t5 +rt;							m7 =qreduce(m5 -rm+q4);		m8 =qreduce(m6 -im+q4);
		t8 =t6 -it;							t6 =t6 +it;							m5 =qreduce(m5 +rm   );		m6 =qreduce(m6 +im   );
																				rm =m5;						im =m6		;			// 0,b
		rt =t5;	t5 =t1 -rt;					t1 =t1 +rt;							m5 =m1 -rm;					m6 =m2 -im	;			// -:-b,b
		it =t6;	t6 =t2 -it;					t2 =t2 +it;							m1 =m1 +rm;					m2 =m2 +im	;			// +:0,2b (rest similar)
																				rm =m7;						im =m8		;
		rt =t7;	t7 =t3 +t8;					t3 =t3 -t8;							m7 =m3 +im;					m8 =m4 -rm	;			// -:-b,b
				t8 =t4 -rt;					t4 =t4 +rt;							m3 =m3 -im;					m4 =m4 +rm	;			// +:0,2b

	/*...Block 2: */
		jt = j1 + p2;		jp = j2 + p2;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t9 =a[jt    ]*c2 -a[jp    ]*s2 ;	t10=a[jp    ]*c2 +a[jt    ]*s2;		cmul_modq8(b[jt    ],b[jp    ],a2 ,b2 ,&m9 ,&m10);
		rt =a[jt+p8 ]*c10-a[jp+p8 ]*s10;	it =a[jp+p8 ]*c10+a[jt+p8 ]*s10;	cmul_modq8(b[jt+p8 ],b[jp+p8 ],a10,b10, &rm, &im);
		t11=t9 -rt;							t9 =t9 +rt;							m11=qreduce(m9 -rm+q4);		m12=qreduce(m10-im+q4);
		t12=t10-it;							t10=t10+it;							m9 =qreduce(m9 +rm   );		m10=qreduce(m10+im   );
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t13=a[jt+p4 ]*c6 -a[jp+p4 ]*s6 ;	t14=a[jp+p4 ]*c6 +a[jt+p4 ]*s6;		cmul_modq8(b[jt+p4 ],b[jp+p4 ],a6 ,b6 ,&m13,&m14);
		rt =a[jt+p12]*c14-a[jp+p12]*s14;	it =a[jp+p12]*c14+a[jt+p12]*s14;	cmul_modq8(b[jt+p12],b[jp+p12],a14,b14, &rm, &im);
		t15=t13-rt;							t13=t13+rt;							m15=qreduce(m13-rm+q4);		m16=qreduce(m14-im+q4);
		t16=t14-it;							t14=t14+it;							m13=qreduce(m13+rm   );		m14=qreduce(m14+im   );
																				rm =m13;					im =m14		;
		rt =t13;	t13=t9 -rt;				t9 =t9 +rt;							m13=m9 -rm;					m14=m10-im	;			// -:-b,b
		it =t14;	t14=t10-it;				t10=t10+it;							m9 =m9 +rm;					m10=m10+im	;			// +:0,2b
																				rm =m15;					im =m16		;
		rt =t15;	t15=t11+t16;			t11=t11-t16;						m15=m11+im;					m16=m12-rm	;			// -:-b,b
					t16=t12-rt;				t12=t12+rt;							m11=m11-im;					m12=m12+rm	;			// +:0,2b

	/*...Block 3: */
		jt = j1 + p1;		jp = j2 + p1;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p4;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t17=a[jt    ]*c1 -a[jp    ]*s1 ;	t18=a[jp    ]*c1 +a[jt    ]*s1;		cmul_modq8(b[jt    ],b[jp    ],a1 ,b1 ,&m17,&m18);
		rt =a[jt+p8 ]*c9 -a[jp+p8 ]*s9 ;	it =a[jp+p8 ]*c9 +a[jt+p8 ]*s9;		cmul_modq8(b[jt+p8 ],b[jp+p8 ],a9 ,b9 , &rm, &im);
		t19=t17-rt;							t17=t17+rt;							m19=qreduce(m17-rm+q4);		m20=qreduce(m18-im+q4);
		t20=t18-it;							t18=t18+it;							m17=qreduce(m17+rm   );		m18=qreduce(m18+im   );
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t21=a[jt+p4 ]*c5 -a[jp+p4 ]*s5 ;	t22=a[jp+p4 ]*c5 +a[jt+p4 ]*s5;		cmul_modq8(b[jt+p4 ],b[jp+p4 ],a5 ,b5 ,&m21,&m22);
		rt =a[jt+p12]*c13-a[jp+p12]*s13;	it =a[jp+p12]*c13+a[jt+p12]*s13;	cmul_modq8(b[jt+p12],b[jp+p12],a13,b13, &rm, &im);
		t23=t21-rt;							t21=t21+rt;							m23=qreduce(m21-rm+q4);		m24=qreduce(m22-im+q4);
		t24=t22-it;							t22=t22+it;							m21=qreduce(m21+rm   );		m22=qreduce(m22+im   );
																				rm =m21;					im =m22		;
		rt =t21;	t21=t17-rt;				t17=t17+rt;							m21=m17-rm;					m22=m18-im	;			// -:-b,b
		it =t22;	t22=t18-it;				t18=t18+it;							m17=m17+rm;					m18=m18+im	;			// +:0,2b
																				rm =m23;					im =m24		;
		rt =t23;	t23=t19+t24;			t19=t19-t24;						m23=qreduce(m19+im   );		m24=qreduce(m20-rm+q2);	// -:-b,b
					t24=t20-rt;				t20=t20+rt;							m19=qreduce(m19-im+q2);		m20=qreduce(m20+rm   );	// +:0,2b
																				// m19,20,23,24 are needed for CMUL, so reduce.
	/*...Block 4: */
		jt = j1 + p3;		jp = j2 + p3;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t25=a[jt    ]*c3 -a[jp    ]*s3 ;	t26=a[jp    ]*c3 +a[jt    ]*s3;		cmul_modq8(b[jt    ],b[jp    ],a3 ,b3 ,&m25,&m26);
		rt =a[jt+p8 ]*c11-a[jp+p8 ]*s11;	it =a[jp+p8 ]*c11+a[jt+p8 ]*s11;	cmul_modq8(b[jt+p8 ],b[jp+p8 ],a11,b11, &rm, &im);
		t27=t25-rt;							t25=t25+rt;							m27=qreduce(m25-rm+q4);		m28=qreduce(m26-im+q4);
		t28=t26-it;							t26=t26+it;							m25=qreduce(m25+rm   );		m26=qreduce(m26+im   );
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t29=a[jt+p4 ]*c7 -a[jp+p4 ]*s7 ;	t30=a[jp+p4 ]*c7 +a[jt+p4 ]*s7;		cmul_modq8(b[jt+p4 ],b[jp+p4 ],a7 ,b7 ,&m29,&m30);
		rt =a[jt+p12]*c15-a[jp+p12]*s15;	it =a[jp+p12]*c15+a[jt+p12]*s15;	cmul_modq8(b[jt+p12],b[jp+p12],a15,b15, &rm, &im);
		t31=t29-rt;							t29=t29+rt;							m31=qreduce(m29-rm+q4);		m32=qreduce(m30-im+q4);
		t32=t30-it;							t30=t30+it;							m29=qreduce(m29+rm   );		m30=qreduce(m30+im   );
																				rm =m29;					im =m30		;
		rt =t29;	t29=t25-rt;				t25=t25+rt;							m29=m25-rm;					m30=m26-im	;			// -:-b,b
		it =t30;	t30=t26-it;				t26=t26+it;							m25=m25+rm;					m26=m26+im	;			// +:0,2b
																				rm =m31;					im =m32		;
		rt =t31;	t31=t27+t32;			t27=t27-t32;						m31=qreduce(m27+im   );		m32=qreduce(m28-rm+q2);	// -:-b,b
					t32=t28-rt;				t28=t28+rt;							m27=qreduce(m27-im+q2);		m28=qreduce(m28+rm   );	// +:0,2b
																				// m27,28,31,32 are needed for CMUL, so reduce.
	/*
    !...and now do four more radix-4 transforms, including the internal twiddle factors:
    !    1, exp(i* 1*twopi/16) =        ( c, s), exp(i* 2*twopi/16) = isqrt2*( 1, 1), exp(i* 3*twopi/16) =        ( s, c) (for inputs to transform block 2)
    !    1, exp(i* 2*twopi/16) = isqrt2*( 1, 1), exp(i* 4*twopi/16) =        ( 0, 1), exp(i* 6*twopi/16) = isqrt2*(-1, 1) (for inputs to transform block 3)
    !    1, exp(i* 3*twopi/16) =        ( s, c), exp(i* 6*twopi/16) = isqrt2*(-1, 1), exp(i* 9*twopi/16) =        (-c,-s) (for inputs to transform block 4).
    !  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
    !   I.e. do similar as above, except inputs a[j1,2+p0:15:1) are replaced by t0:30:2,
    !                                           a[j1,2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
													/*===============
													Input bounds on modular terms:
														-:-b,b: m3,5,6,8,11,13,14,16,21,22,24,29,30,32
														+:0,2b: m1,2,4,7, 9,10,12,15,17,18,23,25,26,31
													These are all reduced (in 0,b) because they are inputs to CMUL:
																m19,20,23,24,27,28,31,32
													===============*/
	/*...Block 1: t1,9,17,25 */
		jt = j1;		jp = j2;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p8 ;
		prefetch_p_doubles(addr);
	  #endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;				rm =m9;	m9 =m1 -rm;	m1 =m1 +rm;	// m 1, 2 in   0,4b
		it =t10;t10=t2 -it;	t2 =t2 +it;				im =m10;m10=m2 -im;	m2 =m2 +im;	// m 9,10 in -2b,2b

		rt =t25;t25=t17-rt;	t17=t17+rt;				rm =m25;m25=m17-rm;	m17=m17+rm;	// m17,18 in   0,4b
		it =t26;t26=t18-it;	t18=t18+it;				im =m26;m26=m18-im;	m18=m18+im;	// m25,26 in -2b,2b
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		/* Debug: check for overflow of + terms: */	ASSERT(m1+m17 >= m1 && m2+m18 >= m2,"Overflow of [0,8b] term!");
		a[jt    ]= t1+t17;	a[jp    ]= t2+t18;		b[jt    ]=qreduce( m1+m17   );	b[jp    ]=qreduce( m2+m18   );	// + terms in   0,8b
		a[jt+p1 ]= t1-t17;	a[jp+p1 ]= t2-t18;		b[jt+p1 ]=qreduce( m1-m17+q4);	b[jp+p1 ]=qreduce( m2-m18+q4);	// - terms in -4b,4b
		// mpy by E^4=i is inlined here:
		a[jt+p2 ]=t9 -t26;	a[jp+p2 ]=t10+t25;		b[jt+p2 ]=qreduce(m9 -m26+q4);	b[jp+p2 ]=qreduce(m10+m25+q4);	// + terms in -4b,4b
		a[jt+p3 ]=t9 +t26;	a[jp+p3 ]=t10-t25;		b[jt+p3 ]=qreduce(m9 +m26+q4);	b[jp+p3 ]=qreduce(m10-m25+q4);	// - terms in -4b,4b

	/*...Block 3: t5,13,21,29 */
		jt = j1 + p4;		jp = j2 + p4;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =t13;t13=t5 +t14;t5 =t5 -t14;					rm =m13;m13=m5 +m14;m5 =m5 -m14;	// m5,6,13,14
				t14=t6 -rt;	t6 =t6 +rt;								m14=m6 -rm;	m6 =m6 +rm;		// all in -2b,2b
	// twiddle mpy by E^2:								// All 4 +- sums in -2b,2b; add q4 to ensure MUL inputs >= 0:
		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;			rm = mul_i2(m21-m22+q4);m22= mul_i2(m21+m22+q4);	// All 4 outs
t21=rt;	rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;	m21=rm;	rm = mul_i2(m30+m29+q4);im = mul_i2(m30-m29+q4);	// in [0,b30]
		t29=t21+rt;			t21=t21-rt;						m29=m21+rm;			m21=m21-rm;		// m21,22 in [-b30,b30]
		t30=t22+it;			t22=t22-it;						m30=m22+im;			m22=m22-im;		// m29,30 in [0,2*b30]
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		a[jt    ]= t5+t21;	a[jp    ]= t6+t22;		b[jt    ]=qreduce( m5+m21+q4);	b[jp    ]=qreduce( m6+m22+q4);	// 5 +-21 in [-2b,2b] + [-b30,b30] = [-2b-b30,2b+b30]
		a[jt+p1 ]= t5-t21;	a[jp+p1 ]= t6-t22;		b[jt+p1 ]=qreduce( m5-m21+q4);	b[jp+p1 ]=qreduce( m6-m22+q4);	// 6 +-22 same
		// mpy by E^4=i is inlined here:
		a[jt+p2 ]=t13-t30;	a[jp+p2 ]=t14+t29;		b[jt+p2 ]=qreduce(m13-m30+q5);	b[jp+p2 ]=qreduce(m14+m29+q3);	// + in [-2b,2b] + [0,2*b30] = [-2b,2*(b+b30)]
		a[jt+p3 ]=t13+t30;	a[jp+p3 ]=t14-t29;		b[jt+p3 ]=qreduce(m13+m30+q3);	b[jp+p3 ]=qreduce(m14-m29+q5);	// - in [-2b,2b] - [0,2*b30] = [-2*(b+b30),2b]

	/*...Block 2: t3,11,19,27 */
		jt = j1 + p8;		jp = j2 + p8;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p12;
		prefetch_p_doubles(addr);
	  #endif
	// twiddle mpy by E^2							// m11-m12 in -3b,b; m11+m12 in -b,3b; rm,im in 0,b30:
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;	rm =mul_i2(m11-m12+q4);	im =mul_i2(m11+m12+q4);
		t11=t3 -rt;			t3 =t3 +rt;				m11=m3 -rm;				m12=m4 -im;	// m11 in -b-b30,b; m12 in -b30,2b
		t12=t4 -it;			t4 =t4 +it;				m3 =m3 +rm;				m4 =m4 +im;	// m3  in -b,b+b30; m4 in 0,2b+b30

		rt =t19*c - t20*s;	t20=t20*c + t19*s;		cmul_modq8(m19,m20, cm,sm, &m19,&m20);	// 0,4q
t19=rt;	rt =t27*s - t28*c;	it =t28*s + t27*c;		cmul_modq8(m27,m28, sm,cm,  &rm, &im);	// 0,4q
		t27=t19-rt;			t19=t19+rt;				m27=qreduce(m19-rm+q4);	m19=qreduce(m19+rm);	// +:   0,8q ==> 0,b
		t28=t20-it;			t20=t20+it;				m28=qreduce(m20-im+q4);	m20=qreduce(m20+im);	// -: -4q,4q ==> 0,b
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;		prefetch_p_doubles(addr);
	  #endif
		a[jt    ]= t3+t19;	a[jp    ]= t4+t20;		b[jt    ]=qreduce(m3+m19+q2);	b[jp    ]=qreduce(m4+m20   );	// m3+m19 in -b,2b+b30; m4+m20 in  0,3b+b30
		a[jt+p1 ]= t3-t19;	a[jp+p1 ]= t4-t20;		b[jt+p1 ]=qreduce(m3-m19+q4);	b[jp+p1 ]=qreduce(m4-m20+q2);	// m3-m19 in -2b,b+b30; m4-m20 in -b,2b+b30
		// mpy by E^4=i is inlined here:
		a[jt+p2 ]=t11-t28;	a[jp+p2 ]=t12+t27;		b[jt+p2 ]=qreduce(m11-m28+q4);	b[jp+p2 ]=qreduce(m12+m27+q2);	// m11-28 in -2b-b30,b; m12+27 in   -b30,3b
		a[jt+p3 ]=t11+t28;	a[jp+p3 ]=t12-t27;		b[jt+p3 ]=qreduce(m11+m28+q4);	b[jp+p3 ]=qreduce(m12-m27+q4);	// m11+28 in -b-b30,2b; m12-27 in -b-b30,2b

	/*...Block 4: t7,15,23,31 */
		jt = j1 + p12;		jp = j2 + p12;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;		prefetch_p_doubles(addr);
	  #endif
		/* twiddle mpy by -E^6 is here... */		// m16+m15 in -b,3b; m16-m15 in -3b,b; rm,im in 0,b30
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;	rm =mul_i2(m16+m15+q4);	im =mul_i2(m16-m15+q4);
		t15=t7 +rt;			t7 =t7 -rt;				m15=m7 +rm;				m16=m8 +im;	// m15 in 0,2b+b30; m16 in -b,b+b30
		t16=t8 +it;			t8 =t8 -it;				m7 =m7 -rm;				m8 =m8 -im;	// m7  in -b30,2b ; m8  in -b-b30,b

		rt =t23*s - t24*c;	t24=t24*s + t23*c;		cmul_modq8(m23,m24, sm,cm, &m23,&m24);	// 0,4q
t23=rt;	rt =t31*c - t32*s;	it =t32*c + t31*s;		cmul_modq8(m31,m32, cm,sm,  &rm, &im);	// 0,4q
		t31=t23+rt;			t23=t23-rt;				m31=qreduce(m23+rm);	m23=qreduce(m23-rm+q4);	// +:   0,8q ==> 0,b
		t32=t24+it;			t24=t24-it;				m32=qreduce(m24+im);	m24=qreduce(m24-im+q4);	// -: -4q,4q ==> 0,b
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		a[jt    ]= t7+t23;	a[jp    ]= t8+t24;		b[jt    ]=qreduce( m7+m23+q4);	b[jp    ]=qreduce( m8+m24+q4);	// 7 +-23 in [-b30, 2b] +- [0,b], handle both with +q4
		a[jt+p1 ]= t7-t23;	a[jp+p1 ]= t8-t24;		b[jt+p1 ]=qreduce( m7-m23+q4);	b[jp+p1 ]=qreduce( m8-m24+q4);	// 8 +-24 in [-b-b30,b] +- [0,b], handle both with +q4
		// mpy by E^4=i is inlined here:
		a[jt+p2 ]=t15-t32;	a[jp+p2 ]=t16+t31;		b[jt+p2 ]=qreduce(m15-m32+q2);	b[jp+p2 ]=qreduce(m16+m31+q2);	// 15+-32 in [0,2b+b30] +- [0,b], handle both with +q2
		a[jt+p3 ]=t15+t32;	a[jp+p3 ]=t16-t31;		b[jt+p3 ]=qreduce(m15+m32+q4);	b[jp+p3 ]=qreduce(m16-m31+q4);	// 16+-31 in [-b,b+b30] +- [0,b], handle both with +q4
							// Cost (aside from basic add/sub): 19 CMUL_MODQ8,  8 MUL_I2, 80 QREDUCE
							// My original f90 version has same 19 CMUL_MODQ8, 12 MUL_I2, 72 QREDUCE, 8 fewer reductions;
							// Difference in MUL_I2 & QREDUCE counts due to differences in way we handle internal twiddles.

			/**********************************************/
	#else	// USE_FGT61 = False; Basic scalar-double mode:
			/**********************************************/

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/
	  #ifdef PFETCH_AGGRESSIVE
		addp = &a[j1];
	  #elif PFETCH
		prefetch_offset = ((j >> 1) & 3)*p4 + 4;	/* Cycle among p0, p4, p8 and p12. */
	  #endif

	/*...Block 1: */
		jt = j1;		jp = j2;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t1 =a[jt	];						t2 =a[jp	];
		rt =a[jt+p8 ]*c8 -a[jp+p8 ]*s8 ;	it =a[jp+p8 ]*c8 +a[jt+p8 ]*s8;
		t3 =t1 -rt;							t1 =t1 +rt;
		t4 =t2 -it;							t2 =t2 +it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t5 =a[jt+p4 ]*c4 -a[jp+p4 ]*s4 ;	t6 =a[jp+p4 ]*c4 +a[jt+p4 ]*s4;
		rt =a[jt+p12]*c12-a[jp+p12]*s12;	it =a[jp+p12]*c12+a[jt+p12]*s12;
		t7 =t5 -rt;							t5 =t5 +rt;
		t8 =t6 -it;							t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;					t1 =t1 +rt;
		it =t6;	t6 =t2 -it;					t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;					t3 =t3 -t8;
				t8 =t4 -rt;					t4 =t4 +rt;

	/*...Block 2: */
		jt = j1 + p2;		jp = j2 + p2;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t9 =a[jt    ]*c2 -a[jp    ]*s2 ;	t10=a[jp    ]*c2 +a[jt    ]*s2;
		rt =a[jt+p8 ]*c10-a[jp+p8 ]*s10;	it =a[jp+p8 ]*c10+a[jt+p8 ]*s10;
		t11=t9 -rt;							t9 =t9 +rt;
		t12=t10-it;							t10=t10+it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t13=a[jt+p4 ]*c6 -a[jp+p4 ]*s6 ;	t14=a[jp+p4 ]*c6 +a[jt+p4 ]*s6;
		rt =a[jt+p12]*c14-a[jp+p12]*s14;	it =a[jp+p12]*c14+a[jt+p12]*s14;
		t15=t13-rt;							t13=t13+rt;
		t16=t14-it;							t14=t14+it;

		rt =t13;	t13=t9 -rt;				t9 =t9 +rt;
		it =t14;	t14=t10-it;				t10=t10+it;

		rt =t15;	t15=t11+t16;			t11=t11-t16;
					t16=t12-rt;				t12=t12+rt;

	/*...Block 3: */
		jt = j1 + p1;		jp = j2 + p1;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p4;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t17=a[jt    ]*c1 -a[jp    ]*s1 ;	t18=a[jp    ]*c1 +a[jt    ]*s1;
		rt =a[jt+p8 ]*c9 -a[jp+p8 ]*s9 ;	it =a[jp+p8 ]*c9 +a[jt+p8 ]*s9;
		t19=t17-rt;							t17=t17+rt;
		t20=t18-it;							t18=t18+it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t21=a[jt+p4 ]*c5 -a[jp+p4 ]*s5 ;	t22=a[jp+p4 ]*c5 +a[jt+p4 ]*s5;
		rt =a[jt+p12]*c13-a[jp+p12]*s13;	it =a[jp+p12]*c13+a[jt+p12]*s13;
		t23=t21-rt;							t21=t21+rt;
		t24=t22-it;							t22=t22+it;

		rt =t21;	t21=t17-rt;				t17=t17+rt;
		it =t22;	t22=t18-it;				t18=t18+it;

		rt =t23;	t23=t19+t24;			t19=t19-t24;
					t24=t20-rt;				t20=t20+rt;

	/*...Block 4: */
		jt = j1 + p3;		jp = j2 + p3;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t25=a[jt    ]*c3 -a[jp    ]*s3 ;	t26=a[jp    ]*c3 +a[jt    ]*s3;
		rt =a[jt+p8 ]*c11-a[jp+p8 ]*s11;	it =a[jp+p8 ]*c11+a[jt+p8 ]*s11;
		t27=t25-rt;							t25=t25+rt;
		t28=t26-it;							t26=t26+it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t29=a[jt+p4 ]*c7 -a[jp+p4 ]*s7 ;	t30=a[jp+p4 ]*c7 +a[jt+p4 ]*s7;
		rt =a[jt+p12]*c15-a[jp+p12]*s15;	it =a[jp+p12]*c15+a[jt+p12]*s15;
		t31=t29-rt;							t29=t29+rt;
		t32=t30-it;							t30=t30+it;

		rt =t29;	t29=t25-rt;				t25=t25+rt;
		it =t30;	t30=t26-it;				t26=t26+it;

		rt =t31;	t31=t27+t32;			t27=t27-t32;
					t32=t28-rt;				t28=t28+rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(i* 1*twopi/16) =	       ( c, s), exp(i* 2*twopi/16) = isqrt2*( 1, 1), exp(i* 3*twopi/16) =        ( s, c) (for inputs to transform block 2)
	!	1, exp(i* 2*twopi/16) = isqrt2*( 1, 1), exp(i* 4*twopi/16) =        ( 0, 1), exp(i* 6*twopi/16) = isqrt2*(-1, 1) (for inputs to transform block 3)
	!	1, exp(i* 3*twopi/16) =        ( s, c), exp(i* 6*twopi/16) = isqrt2*(-1, 1), exp(i* 9*twopi/16) =        (-c,-s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1,2+p0:15:1) are replaced by t0:30:2,
	!										   a[j1,2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
		jt = j1;		jp = j2;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p8 ;
		prefetch_p_doubles(addr);
	  #endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;t10=t2 -it;	t2 =t2 +it;

		rt =t25;t25=t17-rt;	t17=t17+rt;
		it =t26;t26=t18-it;	t18=t18+it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		a[jt    ]=t1+t17;	a[jp    ]=t2+t18;
		a[jt+p1 ]=t1-t17;	a[jp+p1 ]=t2-t18;

		a[jt+p2 ]=t9 -t26;	a[jp+p2 ]=t10+t25;	/* mpy by E^4=i is inlined here... */
		a[jt+p3 ]=t9 +t26;	a[jp+p3 ]=t10-t25;

	/*...Block 3: t5,13,21,29 */
		jt = j1 + p4;		jp = j2 + p4;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =t13;t13=t5 +t14;t5 =t5 -t14;		/* twiddle mpy by E^4 = I */
				t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30 */
		rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;		/* twiddle mpy by -E^6 is here... */
		t29=t21+rt;			t21=t21-rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */
		t30=t22+it;			t22=t22-it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		a[jt    ]=t5+t21;	a[jp    ]=t6+t22;
		a[jt+p1 ]=t5-t21;	a[jp+p1 ]=t6-t22;

		a[jt+p2 ]=t13-t30;	a[jp+p2 ]=t14+t29;	/* mpy by E^4=i is inlined here... */
		a[jt+p3 ]=t13+t30;	a[jp+p3 ]=t14-t29;

	/*...Block 2: t3,11,19,27 */
		jt = j1 + p8;		jp = j2 + p8;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p12;
		prefetch_p_doubles(addr);
	  #endif
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;		/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12 */
		t11=t3 -rt;			t3 =t3 +rt;
		t12=t4 -it;			t4 =t4 +it;

		rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;	/* twiddle mpy by E^1 */
		rt =t27*s - t28*c;	it =t28*s + t27*c;		/* twiddle mpy by E^3 */
		t27=t19-rt;			t19=t19+rt;
		t28=t20-it;			t20=t20+it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		a[jt    ]=t3+t19;	a[jp    ]=t4+t20;
		a[jt+p1 ]=t3-t19;	a[jp+p1 ]=t4-t20;

		a[jt+p2 ]=t11-t28;	a[jp+p2 ]=t12+t27;	/* mpy by E^4=i is inlined here... */
		a[jt+p3 ]=t11+t28;	a[jp+p3 ]=t12-t27;

	/*...Block 4: t7,15,23,31 */
		jt = j1 + p12;		jp = j2 + p12;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;		/* twiddle mpy by -E^6 is here... */
		t15=t7 +rt;			t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */
		t16=t8 +it;			t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;	/* twiddle mpy by E^3 */
		rt =t31*c - t32*s;	it =t32*c + t31*s;		/* twiddle mpy by E^1 = -E^9... */
		t31=t23+rt;			t23=t23-rt;			/* ...and get E^9 by flipping signs here. */
		t32=t24+it;			t24=t24-it;			/* Note: t23+rt = t23*(s+1) */
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		a[jt    ]=t7+t23;	a[jp    ]=t8+t24;
		a[jt+p1 ]=t7-t23;	a[jp+p1 ]=t8-t24;

		a[jt+p2 ]=t15-t32;	a[jp+p2 ]=t16+t31;	/* mpy by E^4=i is inlined here... */
		a[jt+p3 ]=t15+t32;	a[jp+p3 ]=t16-t31;

	#endif	// USE_FGT61 ?

  #endif	// USE_SCALAR_DFT_MACRO ?

	#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

/***************/

/*
!...Acronym: DIT = Decimation In Time
!...Post-twiddles implementation of radix16_dit_pass.
*/
#ifdef USE_FGT61
void radix16_dit_pass	(double a[], uint64 b[], int n, struct complex rt0[], struct complex rt1[], uint128 mt0[], uint128 mt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
#else
void radix16_dit_pass	(double a[],             int n, struct complex rt0[], struct complex rt1[],                               int index[], int nloops, int incr, int init_sse2, int thr_id)
#endif
{
#ifdef USE_SSE2
	const int pfetch_dist = PFETCH_DIST;
	int pfetch_addr;
	static int max_threads = 0;
#endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	// lg(stride):
	const int l2_stride = L2_SZ_VD-2;	// 16 doubles at a time
#endif
#ifdef USE_FGT61
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull, q2=q+q, q3=q2+q, q4=q2+q2, q5=q4+q, q8=q4+q4;	// q = 2^61 - 1, and needed small multiples
	// primitive 16th root of unity, scaled by *8:
	const uint64 cm = 1693317751237720973ull<<3, sm = 2283815672160731785ull<<3;
#endif
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]*/
#if (!defined(USE_AVX2) && defined(USE_SCALAR_DFT_MACRO)) || (defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE) && defined(DFT_V2))	// FMA-based DFT needs the tangent
	const double tan = 0.41421356237309504879;
#endif
	int i,j,j1,jlo,jhi,m,iroot_prim,iroot,k1,k2;
#ifndef USE_SSE2
	int j2;
#endif
	int p1,p2,p3,p4,p8,p12;
#if !defined(USE_SSE2) || defined(USE_AVX512)
	double rt,it;
#endif
	double re0,im0,re1,im1;
#if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
	// These needed both for scalar mode and for certain SIMD-mode inits:
	double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;
#endif

#ifdef USE_SSE2

	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0;	/* Addresses into array sections */

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
  #endif
  #ifndef MULTITHREAD
	static
  #endif
	  vec_dbl *cc0, *isrt2, *two, *r1
	  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
		,*ss0
	  #endif
		;

#else

  #ifndef USE_SCALAR_DFT_MACRO
	int jp,jt;
	double *addr, *addp;
	int prefetch_offset;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
  #endif
  #ifdef USE_FGT61
	uint64 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15
		  ,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32;
  #endif

#endif

#ifndef COMPILER_TYPE_GCC
	ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
#endif

#ifdef USE_SSE2

	/* In order to eschew complex thread-block-and-sync logic related to the local-store-init step in multithread mode,
	switch to a special-init-mode-call paradigm, in which the function is inited once (for as many threads as needed)
	prior to being executed:
	*/
	if(init_sse2)	// Just check (init_sse2 != 0) here, to allow the *value* of init_sse2 to store #threads
	{
		max_threads = init_sse2;
	#ifndef COMPILER_TYPE_GCC
		ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(thr_id == -1, "Init-mode call must be outside of any multithreading!");
		if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)sc_arr);	sc_arr=0x0;
		}
		sc_arr = ALLOC_VEC_DBL(sc_arr, 72*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	last 30 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
		#ifdef MULTITHREAD
	//	if(max_threads > 1) {
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x20;
			cc0   = sc_ptr + 0x21;
		  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
			ss0   = sc_ptr + 0x22;
		  #endif
			two   = sc_ptr + 0x43;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);
			  #if defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE)
				VEC_DBL_INIT(two  , 1.0);	// Yeah, I *know" "it's a misnomer" :)
				// cc0,ss0 inited below for AVX2
			  #else
				VEC_DBL_INIT(two  , 2.0);
				VEC_DBL_INIT(cc0  , c);
				VEC_DBL_INIT(ss0  , s);
			  #endif
				isrt2 += 72;	/* Move on to next thread's local store */
				cc0   += 72;
			  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
				ss0   += 72;
			  #endif
				two   += 72;
			}
		#elif defined(COMPILER_TYPE_GCC)
			r1  = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;
										cc0 = sc_ptr + 0x21;
		  #if !(defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE))
										ss0 = sc_ptr + 0x22;
		  #endif
										two = sc_ptr + 0x43;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
		  #if defined(USE_AVX2) && !defined(REFACTOR_4DFT_3TWIDDLE)
			VEC_DBL_INIT(two  , 1.0);	// Yeah, I *know" "it's a misnomer" :)
			// cc0,ss0 inited below for AVX2
		  #else
			VEC_DBL_INIT(two  , 2.0);
			VEC_DBL_INIT(cc0  , c);
			VEC_DBL_INIT(ss0  , s);
		  #endif
		#else
		  #error Non-GCC-compatible compilers not supported for SIMD builds!
		#endif
	//	}
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		r1 = __r0 + thr_id*72;
		isrt2 = r1 + 0x20;
		cc0   = r1 + 0x21;
		two   = r1 + 0x43;
	#endif

#endif

	p1 = incr >> 4;
	p2 = p1 +p1;
	p3 = p2 +p1;
	p4 = p3 +p1;
	p8 = p4 +p4;
	p12= p8 +p4;
	p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
	p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
	// Make sure to at least one-time pre-test any index arithmetic assumptions used to save cycles in the loop
	// body (both C and ASM). Since such checks may be runlength-dependent, need to be cheap enough to leave on
	// all the time, as here where we do them just once prior to entering the processing loop. Since DIF and DIT
	// encounter the same sets or index strides (albeit in opposite order), can split such tests between them:
	ASSERT(p4  == p2+p2, "radix16_dit_pass: p4  != p2+p2!");
	ASSERT(p8  == p4+p4, "radix16_dit_pass: p8  != p4+p4!");
	ASSERT(p12 == p4+p8, "radix16_dit_pass: p12 != p4+p8!");

	iroot_prim=(incr >> 5);		/* (incr/2)/radix_now */

	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16).	 */
	{				/* NLOOPS satisfies the identity NLOOPS * INCR = P*N, which says that we process the entire */
					/* array each time this subroutine is executed (since P*N = vector length, sans padding.)   */

/*	here are the needed sincos data - these are processed below in bit-reversed order. */
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* In the post-twiddle-mul version of the DIT pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks.
*/
#if HIACC

	#ifdef REFACTOR_4DFT_3TWIDDLE

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c1=re0*re1-im0*im1;	s1=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c2=re0*re1-im0*im1;	s2=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c3=re0*re1-im0*im1;	s3=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c4=re0*re1-im0*im1;	s4=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c5=re0*re1-im0*im1;	s5=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c6=re0*re1-im0*im1;	s6=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c7=re0*re1-im0*im1;	s7=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c8=re0*re1-im0*im1;	s8=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c9=re0*re1-im0*im1;	s9=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c10=re0*re1-im0*im1;	s10=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c11=re0*re1-im0*im1;	s11=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c12=re0*re1-im0*im1;	s12=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c13=re0*re1-im0*im1;	s13=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c14=re0*re1-im0*im1;	s14=re0*im1+im0*re1;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		c15=re0*re1-im0*im1;	s15=re0*im1+im0*re1;

		// Cf. util.c:test_radix16_dft() for details on these combos:
		c5 =  (c1-s1)*ISRT2;	s5 = (c1+s1)*ISRT2;
		c6 = -s2;				s6 = c2;
		c7 = -(c3+s3)*ISRT2;	s7 = (c3-s3)*ISRT2;

		c9 = c1*c-s1*s;			s9 = c1*s+s1*c;
		c10= (c2-s2)*ISRT2;		s10= (c2+s2)*ISRT2;
		c11= c3*s-s3*c;			s11= c3*c+s3*s;

		c13= c1*s-s1*c;			s13= c1*c+s1*s;
		c14= -(c2+s2)*ISRT2;	s14= (c2-s2)*ISRT2;
		c15= -c3*c+s3*s;		s15= -(c3*s+s3*c);

	  #ifdef USE_SSE2
		/* Sincos data stored in terms of the following 5 contiguous-data triplets:
			c4,s4, c8,s8, cC,sC
			c1,s1, c2,s2, c3,s3
			c5,s5, c6,s6, c7,s7
			c9,s9, cA,sA, cB,sB
			cD,sD, cE,sE, cF,sF .
		Note that due to my layout of the SSE2_RADIX_04_DIF_3TWIDDLE_X2-macro arglist,
		we need to swap the order of the first 2 sincos-pairs of each triplet:
		*/
		vec_dbl *c_tmp = cc0, *s_tmp = c_tmp+1;	/* c0,s0 */
		VEC_DBL_INIT(c_tmp, c8 );	VEC_DBL_INIT(s_tmp, s8 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c4 );	VEC_DBL_INIT(s_tmp, s4 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c12);	VEC_DBL_INIT(s_tmp, s12);	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c2 );	VEC_DBL_INIT(s_tmp, s2 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c1 );	VEC_DBL_INIT(s_tmp, s1 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c3 );	VEC_DBL_INIT(s_tmp, s3 );	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c6 );	VEC_DBL_INIT(s_tmp, s6 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c5 );	VEC_DBL_INIT(s_tmp, s5 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c7 );	VEC_DBL_INIT(s_tmp, s7 );	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c10);	VEC_DBL_INIT(s_tmp, s10);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c9 );	VEC_DBL_INIT(s_tmp, s9 );	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c11);	VEC_DBL_INIT(s_tmp, s11);	c_tmp+=2; s_tmp+=2;

		VEC_DBL_INIT(c_tmp, c14);	VEC_DBL_INIT(s_tmp, s14);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c13);	VEC_DBL_INIT(s_tmp, s13);	c_tmp+=2; s_tmp+=2;
		VEC_DBL_INIT(c_tmp, c15);	VEC_DBL_INIT(s_tmp, s15);	c_tmp+=2; s_tmp+=2;
	  #endif

	#elif defined(USE_AVX2)
		// In AVX2/FMA mode, since we need to replace most of the raw sincos data with derived ones,
		// simply place one copy of each computed double in a double-sized slot of the local memory.
		// We will be using AVX2/FMA-based Newtonian iterative inversion on the 16 doubles whose
		// multiplicative inverse is needed (the real part of the basic root of unity c and of the
		// 15 complex twiddles, c1-15), so store those in packed form in 4 AVX-register-sized
		// contiguous memory locations, and the others in a separate chunk of memory. After the
		// vector-iterative inversion we'll need to combine the 2 sets of data and place (in quadruplicate)
		// into their final SIMD-suitable memory slots.

		add0 = (double *)cc0;			// add0 points to 16 cos-data-to-be-inverted; Need a double-ptr on lhs here
		double *add1 = add0 + 16;		// add1 points to block of memory temporarily used to store the corresponding sine data
		*add0++ = c;	// Since tan0 defined as const, we can init these directly, but init with c0,s0 anyway
		*add1++ = s;	// and use result as a check onthe accuracy of the FMA-based Newton iterative inversion.

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c1, for inversion
		*add1++ = it;	// s1  slot will hold __r1 = s1 /c1

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c2, for inversion
		*add1++ = it;	// s2  slot will hold __r2 = s2 /c2

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c3, for inversion
		*add1++ = it;	// s3  slot will hold __r3 = s3 /c3

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c4, for inversion
		*add1++ = it;	// s4  slot will hold __r4 = s4 /c4

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c5, for inversion
		*add1++ = it;	// s5  slot will hold __r5 = s5 /c5

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c6, for inversion
		*add1++ = it;	// s6  slot will hold __r6 = s6 /c6

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c7, for inversion
		*add1++ = it;	// s7  slot will hold __r7 = s7 /c7

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c8, for inversion
		*add1++ = it;	// s8  slot will hold __r8 = s8 /c8

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c9, for inversion
		*add1++ = it;	// s9  slot will hold __r9 = s9 /c9

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c10, for inversion
		*add1++ = it;	// s10 slot will hold __rA = s10/c10

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c11, for inversion
		*add1++ = it;	// s11 slot will hold __rB = s11/c11

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c12, for inversion
		*add1++ = it;	// s12 slot will hold __rC = s12/c12

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c13, for inversion
		*add1++ = it;	// s13 slot will hold __rD = s13/c13

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c14, for inversion
		*add1++ = it;	// s14 slot will hold __rE = s14/c14

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c15, for inversion
		*add1++ = it;	// s15 slot will hold __rF = s15/c15

		// This places us at add0 == c8 and add1 = c12.
		ASSERT(add0 == (double *)cc0+16 && add1 == (double *)cc0+32, "add0,1 checksum failed in AVX2 DIT sincos inits!");
		/*
		At this point, the 8 ymm-sized [32-byte] chunks starting at &cc0 contain the following scalar-double data:

		0:	c,c1-3		4:	s,s1-3
		1:	c4-7		5:	s4-7
		2:	c8-B		6:	s8-B
		3:	cC-F		7:	sC-F
		*/

		// Now send the cosine terms to the inversion routine, which also does the combine-and-populate-SIMD-slots step.
	   #ifdef DFT_V1	// Toggle between 2 versions

		RADIX16_COMPUTE_FMA_SINCOS_DIT(cc0,two);

	   #elif defined(DFT_V2)

		RADIX16_COMPUTE_FMA_SINCOS_DIT_2(cc0);
		add0 = (double *)cc0;
		add0[0x00] = c;
		add0[0x10] = tan;
		add0[0x20] = 1.0;
	//	ASSERT(*(add0-1) == ISRT2, "Scalar ISRT2 bad!");
		c_tmp = cc0 + 0x22;	// 1.0 x 4
	//	ASSERT(c_tmp->d0 == 1.0 && c_tmp->d0 == c_tmp->d1 && c_tmp->d0 == c_tmp->d2 && c_tmp->d0 == c_tmp->d3, "1.0 x 4 mismatch!");

		/* Scalar data starting at add0 = cc0 now laid out as below:

		a[-0x4,-0x3,-0x2,-0x1]: add0 - 0x[ 20, 18, 10, 08]: [   ISRT2 x 4   ]
		a[0x00,0x01,0x02,0x03]: add0 + 0x[  0,  8, 10, 18]: [__c,c1 ,c2 ,c3 ] Cosines:
		a[0x04,0x05,0x06,0x07]: add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ]
		a[0x08,0x09,0x0a,0x0b]: add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ]
		a[0x0c,0x0d,0x0e,0x0f]: add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ]
		a[0x10,0x11,0x12,0x13]: add0 + 0x[ 80, 88, 90, 98]: [_sc,r1 ,r2 ,r3 ] Tangents:
		a[0x14,0x15,0x16,0x17]: add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ]
		a[0x18,0x19,0x1a,0x1b]: add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ]
		a[0x1c,0x1d,0x1e,0x1f]: add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ]
		a[0x20]               : add0 + 0x[100]            : [1.0]
		*/
	  #endif

	#else	// AVX2 = False:

	  #ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-15] are offset w.r.to the thread-local ptr pair as
		(cc0,ss0) + 0x[2,12,a,1a, 4,14,c,1c, 6,16,e,1e, 8,18,10,20], i.e. ordered [0,4,8,C, 2,6,A,E, 1,5,9,D, 3,7,B,F];
		middle 2 quartet adress-offsets swapped [4,14,c,1c] <--> [6,16,e,1e] w.r.to DIF, which has effect of swapping
		the middle 2 elements of each quartet with respect to DIF's strict BR-ordering:
	*** May 2016: wtf was I thinking w.r.to this middle-two-roots-quartets-swap? Swapped address-offset quartets with
		idx == 4,6 (mod 8) in inits below and fiddled SSE2_RADIX16_DIT_TWIDDLE macro correspondingly, so now DIF,DIT
		share same layout: roots (c,s)[0-15] are offset w.r.to the thread-local ptr pair as
		(cc0,ss0) + 0x[2,12,a,1a, 6,16,e,1e, 4,14,c,1c, 8,18,10,20], i.e. BR-ordered [0,8,4,C, 2,A,6,E, 1,9,5,D, 3,B,7,F].
		*/
		c_tmp = cc0 + 0x02; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_modq(mt0[k1].d0,mt0[k1].d1, mt1[k2].d0,mt1[k2].d1, &rm,&im);
		// Premultiply results by 8 so ensuing twiddle-MULs can use the streamlined CMUL_MODQ8 routine:
		a1 = qreduce_full(rm)<<3;	b1 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x12; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c1 =rt;		s1 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		csqr_modq(a1,b1, &rm,&im);	// Since MODMULs are roundoff-error-free, cheaper to use a LOACC-style sequence
						// of successive squarings to get E^2,4,8, a CMUL to get E^13, then CMUL_CONJ_MODQs to get rest.
		a2 = qreduce_full(rm)<<3;	b2 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x0a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c2 =rt;		s2 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x1a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c3 =rt;		s3 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		csqr_modq(a2,b2, &rm,&im);
		a4 = qreduce_full(rm)<<3;	b4 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c4 =rt;		s4 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a1,b1, a4,b4, &a3,&b3, &a5,&b5);
		a3 = qreduce_full(a3)<<3;	b3 = qreduce_full(b3)<<3;
		a5 = qreduce_full(a5)<<3;	b5 = qreduce_full(b5)<<3;
	  #endif
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c5 =rt;		s5 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x0e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c6 =rt;		s6 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x1e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c7 =rt;		s7 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		csqr_modq(a4,b4, &rm,&im);
		a8 = qreduce_full(rm)<<3;	b8 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x04; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c8 =rt;		s8 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a1,b1, a8,b8, &a7,&b7, &a9,&b9);
		a7 = qreduce_full(a7)<<3;	b7 = qreduce_full(b7)<<3;
		a9 = qreduce_full(a9)<<3;	b9 = qreduce_full(b9)<<3;
	  #endif
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c9 =rt;		s9 =it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a2,b2, a8,b8, &a6,&b6, &a10,&b10);
		a6  = qreduce_full(a6 )<<3;	b6  = qreduce_full(b6 )<<3;
		a10 = qreduce_full(a10)<<3;	b10 = qreduce_full(b10)<<3;
	  #endif
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x0c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c10=rt;		s10=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x1c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c11=rt;		s11=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x08; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c12=rt;		s12=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_modq(mt0[k1].d0,mt0[k1].d1, mt1[k2].d0,mt1[k2].d1, &rm,&im);
		a13 = qreduce_full(rm)<<3;	b13 = qreduce_full(im)<<3;
	  #endif
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x18; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c13=rt;		s13=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a1,b1, a13,b13, &a12,&b12, &a14,&b14);
		a12 = qreduce_full(a12)<<3;	b12 = qreduce_full(b12)<<3;
		a14 = qreduce_full(a14)<<3;	b14 = qreduce_full(b14)<<3;
	  #endif
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x10; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c14=rt;		s14=it;
	  #endif

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
	  #ifdef USE_FGT61
		cmul_conj_modq(a2,b2, a13,b13, &a11,&b11, &a15,&b15);
		a11 = qreduce_full(a11)<<3;	b11 = qreduce_full(b11)<<3;
		a15 = qreduce_full(a15)<<3;	b15 = qreduce_full(b15)<<3;
	  #endif
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x20; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c15=rt;		s15=it;
	  #endif

	#endif	// USE_AVX2?

#else
		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c1=t1*t3-t2*t4;	s1=t1*t4+t2*t3;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += i;				/* 4*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c2=t1*t3-t2*t4;	s2=t1*t4+t2*t3;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += i;				/* 8*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c4=t1*t3-t2*t4;	s4=t1*t4+t2*t3;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		i += (iroot << 2)+iroot;		/* 13*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c8=t1*t3-t2*t4;	s8=t1*t4+t2*t3;

		k1=(i & NRTM1);	k2=(i >> NRT_BITS);
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c13=t1*t3-t2*t4;	s13=t1*t4+t2*t3;

		/* c3,5 */
		t1=c1*c4; t2=c1*s4; t3=s1*c4; t4=s1*s4;
		c3=t1+t4; s3=t2-t3; c5=t1-t4; s5=t2+t3;

		/* c6,7,9,10 */
		t1=c1*c8; t2=c1*s8; t3=s1*c8; t4=s1*s8;
		c7=t1+t4; s7=t2-t3; c9=t1-t4; s9=t2+t3;

		t1=c2*c8; t2=c2*s8; t3=s2*c8; t4=s2*s8;
		c6=t1+t4; s6=t2-t3; c10=t1-t4; s10=t2+t3;

		/* c11,12,14,15 */
		t1=c1*c13; t2=c1*s13; t3=s1*c13; t4=s1*s13;
		c12=t1+t4; s12=t2-t3; c14=t1-t4; s14=t2+t3;

		t1=c2*c13; t2=c2*s13; t3=s2*c13; t4=s2*s13;
		c11=t1+t4; s11=t2-t3; c15=t1-t4; s15=t2+t3;
#endif

#ifndef USE_SSE2
  #ifdef USE_SCALAR_DFT_MACRO	// Must define - or not - @compile time
	// FMA-based version replaces sine terms with tangents:
	s1  /= c1;
	s2  /= c2;
	s3  /= c3;
	s4  /= c4;
	s5  /= c5;
	s6  /= c6;
	s7  /= c7;
	s8  /= c8;
	s9  /= c9;
	s10 /= c10;
	s11 /= c11;
	s12 /= c12;
	s13 /= c13;
	s14 /= c14;
	s15 /= c15;
  #endif
#endif

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 4);

	/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
	Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
	but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
	(But only ever need to explicitly do this in debug mode).
	*/
	  for(j = jlo; j < jhi; j += stride)
	  {
		j1 = j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	#ifdef USE_SSE2

		pfetch_addr = p1*((j >> l2_stride) & 0x3);	// cycle prefetch-offset-address among p0,1,2,3
			// These get added to the base addresses p0,4,8,12 in thr DFT macros, thus every 4 loop
			// executions we have covered prefetches from [current address] + [pfetch distance] + p0,1,2,...15 .

		add0 = &a[j1];

	  #ifdef USE_AVX2

	   #ifdef REFACTOR_4DFT_3TWIDDLE
		// Since pvsly had no non-suffixed 'SSE2_RADIX16_DIF_TWIDDLE' macro in AVX2 mode, no need for _V2 suffix here:
		SSE2_RADIX16_DIT_TWIDDLE(add0,p1,p2,p3,p4,p8,p12,r1,two,cc0,pfetch_addr,pfetch_dist);

	   #elif defined(DFT_V1)	// Toggle between 2 versions

		SSE2_RADIX16_DIT_TWIDDLE_1(add0,p1,p2,p3,p4,p8,p12,r1,cc0,pfetch_addr,pfetch_dist);

	   #elif defined(DFT_V2)

		SSE2_RADIX16_DIT_TWIDDLE_2(add0,p1,p2,p3,p4,p8,p12,r1,cc0,pfetch_addr,pfetch_dist);

	   #else	// Allow non-FMA as the default option for comparative timing purposes

		SSE2_RADIX16_DIT_TWIDDLE_0(add0,p1,p2,p3,p4,p8,r1,isrt2,pfetch_addr,pfetch_dist);

	   #endif

	  // SSE2 code has different macro arglists for 32 and 64-bit modes:
	  #else

	   #ifdef REFACTOR_4DFT_3TWIDDLE

		SSE2_RADIX16_DIT_TWIDDLE_V2(add0,p1,p2,p3,p4,p8,p12,r1,two,cc0,pfetch_addr,pfetch_dist);

	   #else

		SSE2_RADIX16_DIT_TWIDDLE(add0,p1,p2,p3,p4,p8,r1,isrt2,pfetch_addr,pfetch_dist);

	   #endif

	  #endif

	#else	/* USE_SSE2 */

	  #ifdef USE_SCALAR_DFT_MACRO	// Must define - or not - @compile time

		// Test FMA-based DIT macro: 'sine' terms here are tangents!
		RADIX_16_DIT_FMA(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p4+p1 ],a[j2+p4+p1 ],a[j1+p4+p2 ],a[j2+p4+p2 ],a[j1+p4+p3 ],a[j2+p4+p3 ],a[j1+p8 ],a[j2+p8 ],a[j1+p8+p1 ],a[j2+p8+p1 ],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+p12],a[j2+p12],a[j1+p12+p1],a[j2+p12+p1],a[j1+p12+p2],a[j2+p12+p2],a[j1+p12+p3],a[j2+p12+p3]
						,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p4+p1 ],a[j2+p4+p1 ],a[j1+p4+p2 ],a[j2+p4+p2 ],a[j1+p4+p3 ],a[j2+p4+p3 ],a[j1+p8 ],a[j2+p8 ],a[j1+p8+p1 ],a[j2+p8+p1 ],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+p12],a[j2+p12],a[j1+p12+p1],a[j2+p12+p1],a[j1+p12+p2],a[j2+p12+p2],a[j1+p12+p3],a[j2+p12+p3]
						,c1,s1,c2,s2,c3,s3,c4,s4,c5,s5,c6,s6,c7,s7,c8,s8,c9,s9,c10,s10,c11,s11,c12,s12,c13,s13,c14,s14,c15,s15
						,c,tan)

	  #else		// USE_SCALAR_DFT_MACRO = False

	#ifdef USE_FGT61

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/
	  #ifdef PFETCH_AGGRESSIVE
		addp = &a[j1];
	  #elif PFETCH
		prefetch_offset = ((j >> 1) & 3)*p1 + 4;	/* Cycle among p0, p1, p2 and p3. */
	  #endif
	/*...Block 1: */
		jt = j1;		jp = j2;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t1 =a[jt ];		t2 =a[jp ];					m1 =b[jt ];		m2 =b[jp ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t3 =t1 -rt;		t1 =t1 +rt;					m3 =m1 -rm;		m1 =m1 +rm;		// 1,2 in 0,2b
		t4 =t2 -it;		t2 =t2 +it;					m4 =m2 -im;		m2 =m2 +im;		// 3,4 in -b,b
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t5 =a[jt+p2 ];	t6 =a[jp+p2 ];				m5 =b[jt+p2 ];	m6 =b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t7 =t5 -rt;		t5 =t5 +rt;					m7 =m5 -rm;		m5 =m5 +rm;		// 5,6 in 0,2b
		t8 =t6 -it;		t6 =t6 +it;					m8 =m6 -im;		m6 =m6 +im;		// 7,8 in -b,b

		rt =t5;	t5 =t1 -rt ;	t1 =t1 +rt;			rm =m5;	m5 =m1 -rm ;	m1 =m1 +rm;	// 1,2 in   0,4b
		it =t6;	t6 =t2 -it ;	t2 =t2 +it;			im =m6;	m6 =m2 -im ;	m2 =m2 +im;	// 5,6 in -2b,2b

		rt =t7;	t7 =t3 -t8 ;	t3 =t3 +t8;			rm =m7;	m7 =m3 -m8 ;	m3 =m3 +m8;	// 3,8 in -2b,2b
				t8 =t4 +rt ;	t4 =t4 -rt;					m8 =m4 +rm ;	m4 =m4 -rm;	// 4,7 same

	/*...Block 2: */
		jt = j1 + p4;		jp = j2 + p4;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t9 =a[jt    ];	t10=a[jp    ];				m9 =b[jt    ];	m10=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t11=t9 -rt;		t9 =t9 +rt;					m11=m9 -rm;		m9 =m9 +rm;		//  9,10 in 0,2b
		t12=t10-it;		t10=t10+it;					m12=m10-im;		m10=m10+im;		// 11,12 in -b,b
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t13=a[jt+p2 ];	t14=a[jp+p2 ];				m13=b[jt+p2 ];	m14=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t15=t13-rt;		t13=t13+rt;					m15=m13-rm;		m13=m13+rm;		// 13,14 in 0,2b
		t16=t14-it;		t14=t14+it;					m16=m14-im;		m14=m14+im;		// 15,16 in -b,b

		rt =t13;	t13=t9 -rt ;	t9 =t9 +rt;		rm =m13;	m13=m9 -rm ;	m9 =m9 +rm;
		it =t14;	t14=t10-it ;	t10=t10+it;		im =m14;	m14=m10-im ;	m10=m10+im;

		rt =t15;	t15=t11-t16;	t11=t11+t16;	rm =m15;	m15=m11-m16;	m11=m11+m16;
					t16=t12+rt ;	t12=t12-rt;					m16=m12+rm ;	m12=m12-rm;

	/*...Block 3: */
		jt = j1 + p8;		jp = j2 + p8;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p4 ;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t17=a[jt    ];	t18=a[jp    ];				m17=b[jt    ];	m18=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t19=t17-rt;		t17=t17+rt;					m19=m17-rm;		m17=m17+rm;		// 17,18 in 0,2b
		t20=t18-it;		t18=t18+it;					m20=m18-im;		m18=m18+im;		// 19,20 in -b,b
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t21=a[jt+p2 ];	t22=a[jp+p2 ];				m21=b[jt+p2 ];	m22=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t23=t21-rt;		t21=t21+rt;					m23=m21-rm;		m21=m21+rm;	// 21,22 in 0,2b
		t24=t22-it;		t22=t22+it;					m24=m22-im;		m22=m22+im;	// 23,24 in -b,b

		rt =t21;	t21=t17-rt ;	t17=t17+rt;		rm =m21;	m21=m17-rm ;	m17=m17+rm;	// 17,18 in   0,4b
		it =t22;	t22=t18-it ;	t18=t18+it;		im =m22;	m22=m18-im ;	m18=m18+im;	// 21,22 in -2b,2b
													// Prior to qreduce(...+q4), all 4 outs here in -2b,2b:
		rt =t23;	t23=t19-t24;	t19=t19+t24;	rm =m23;	m23=qreduce(m19-m24+q4);	m19=qreduce(m19+m24+q4);
					t24=t20+rt ;	t20=t20-rt;					m24=qreduce(m20+rm +q4);	m20=qreduce(m20-rm +q4);
													// m19,20,23,24 all needed for CMUL, so reduce.
	/*...Block 4: */
		jt = j1 + p12;		jp = j2 + p12;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t25=a[jt    ];	t26=a[jp    ];				m25=b[jt    ];	m26=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t27=t25-rt;		t25=t25+rt;					m27=m25-rm;		m25=m25+rm;	// 25,26 in 0,2b
		t28=t26-it;		t26=t26+it;					m28=m26-im;		m26=m26+im;	// 27,28 in -b,b
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		t29=a[jt+p2 ];	t30=a[jp+p2 ];				m29=b[jt+p2 ];	m30=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t31=t29-rt;		t29=t29+rt;					m31=m29-rm;		m29=m29+rm;	// 29,30 in 0,2b
		t32=t30-it;		t30=t30+it;					m32=m30-im;		m30=m30+im;	// 31,32 in -b,b

		rt =t29;	t29=t25-rt ;	t25=t25+rt;		rm =m29;	m29=m25-rm ;	m25=m25+rm;
		it =t30;	t30=t26-it ;	t26=t26+it;		im =m30;	m30=m26-im ;	m26=m26+im;
													// Prior to qreduce(...+q4), all 4 outs here in -2b,2b:
		rt =t31;	t31=t27-t32;	t27=t27+t32;	rm =m31;	m31=qreduce(m27-m32+q4);	m27=qreduce(m27+m32+q4);
					t32=t28+rt ;	t28=t28-rt;					m32=qreduce(m28+rm +q4);	m28=qreduce(m28-rm +q4);
													// m27,28,31,32 all needed for CMUL, so reduce.
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
	!	1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
	!	1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[jt  +p0:15:1) are replaced by t0:30:2,
	!										   a[jp+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
													/*===============
													Input bounds on modular terms:
														-2b,2b: m3-8, 11-16, 19-24, 27-32
														  0,4b: m1,2, 9 ,10, 17,18, 25,26
													These are all reduced (in 0,b) because they are inputs to CMUL:
																m19,20,23,24,27,28,31,32
													===============*/
	/*...Block 1: t1,9,17,25 */
		jt = j1;		jp = j2;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p8 ;
		prefetch_p_doubles(addr);
	  #endif							// We optimistically assume outputs bounded by 8b = 2^64 + 48 never > 2^64 in practice:
		rt =t9 ;	t9 =t1 -rt;	t1 =t1 +rt;							rm =m9 ;	m9 =qreduce(m1 -rm+q4);	m1 =qreduce(m1 +rm);	// 1, 2 in   0,8b
		it =t10;	t10=t2 -it;	t2 =t2 +it;							im =m10;	m10=qreduce(m2 -im+q4);	m2 =qreduce(m2 +im);	// 9,10 in -4b,4b

		rt =t25;	t25=t17-rt;	t17=t17+rt;							rm =m25;	m25=qreduce(m17-rm+q4);	m17=qreduce(m17+rm);
		it =t26;	t26=t18-it;	t18=t18+it;							im =m26;	m26=qreduce(m18-im+q4);	m18=qreduce(m18+im);

		a[jt    ]=t1+t17;			a[jp    ]=t2+t18;				b[jt    ] = qreduce(m1+m17   );	b[jp    ] = qreduce(m2+m18   );
		t1	     =t1-t17;			t2		 =t2-t18;				m1	      = qreduce(m1-m17+q2);	m2		  = qreduce(m2-m18+q2);
		a[jt+p8 ]=t1 *c8 +t2 *s8;	a[jp+p8 ]=t2 *c8 -t1 *s8;		cmul_modq8(m1,m2, a8,q8-b8, &rm,&im);
																	b[jt+p8 ] = qreduce(rm);	b[jp+p8 ] = qreduce(im);
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		/* mpy by E^-4 = -I is inlined here... */
		rt	   =t9 +t26;		it		=t10-t25;					rm  = qreduce(m9 +m26   );	im	 = qreduce(m10-m25+q2);
		t9	   =t9 -t26;		t10		=t10+t25;					m9  = qreduce(m9 -m26+q2);	m10	 = qreduce(m10+m25   );
		a[jt+p4 ]=rt *c4 +it *s4;	a[jp+p4 ]=it *c4 -rt *s4 ;		cmul_modq8(rm, im, a4 ,q8-b4 , &rm,&im);
		a[jt+p12]=t9 *c12+t10*s12;	a[jp+p12]=t10*c12-t9 *s12;		cmul_modq8(m9,m10, a12,q8-b12, &m9,&m10);
																	b[jt+p4 ] = qreduce(rm);	b[jp+p4 ] = qreduce( im);
																	b[jt+p12] = qreduce(m9);	b[jp+p12] = qreduce(m10);
	/*...Block 3: t5,13,21,29 */
		jt = j1 + p2;		jp = j2 + p2;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;					rm =m13;m13=qreduce(m5-m14+q4);	m5 =qreduce(m5 +m14+q4);	// all 4 outs in -4b,4b;
					t14=t6 +rt ;	t6 =t6 -rt ;							m14=qreduce(m6+rm +q4);	m6 =qreduce(m6 -rm +q4);	// reduce all 4 to 0,b.

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;			rm = mul_i2(m22+m21+q4);	m22= mul_i2(m22-m21+q4);	m21=rm;
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;					rm = mul_i2(m29-m30+q4);	im = mul_i2(m29+m30+q4);	// m21,22,rm,im in 0,b30
		t29=t21+rt;		t21=t21-rt;									m29 = m21+rm;		m21 = m21-rm;	// m21,22 in -b30,b30
		t30=t22+it;		t22=t22-it;									m30 = m22+im;		m22 = m22-im;	// m29,30 in 0,2*b30

		rt	   =t5 +t21;		it		 =t6 +t22;					rm = qreduce(m5 +m21+q2);	im =qreduce(m6 +m22+q2);	// rm,im in -b30,b+b30
		t5	   =t5 -t21;		t6		 =t6 -t22;					m5 = qreduce(m5 -m21+q2);	m6 =qreduce(m6 -m22+q2);	// m5,m6 in -b30,b+b30; reduce all to 0,b
		a[jt    ]=rt *c2 +it *s2 ;	a[jp    ]=it *c2 -rt *s2 ;		cmul_modq8(rm,im, a2 ,q8-b2 , &rm,&im);
		a[jt+p8 ]=t5 *c10+t6 *s10;	a[jp+p8 ]=t6 *c10-t5 *s10;		cmul_modq8(m5,m6, a10,q8-b10, &m5,&m6);
																	b[jt    ] = qreduce(rm);	b[jp    ] = qreduce(im);
																	b[jt+p8 ] = qreduce(m5);	b[jp+p8 ] = qreduce(m6);
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt	  =t13+t30;		it		=t14-t29;						rm  = qreduce(m13+m30+q4);	im  = qreduce(m14-m29+q4);	// rm,m14 in 0,b+2*b30
		t13	  =t13-t30;		t14		=t14+t29;						m13 = qreduce(m13-m30+q4);	m14 = qreduce(m14+m29+q4);	// im,m13 in -2*b30,b
		a[jt+p4 ]=rt *c6 +it *s6 ;	a[jp+p4 ]=it *c6 -rt *s6 ;		cmul_modq8( rm, im, a6 ,q8-b6 , &rm ,&im );
		a[jt+p12]=t13*c14+t14*s14;	a[jp+p12]=t14*c14-t13*s14;		cmul_modq8(m13,m14, a14,q8-b14, &m13,&m14);
																	b[jt+p4 ] = qreduce( rm);	b[jp+p4 ] = qreduce( im);
																	b[jt+p12] = qreduce(m13);	b[jp+p12] = qreduce(m14);
	/*...Block 2: t3,11,19,27 */
		jt = j1 + p1;		jp = j2 + p1;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;					rm = mul_i2(m12+m11+q4);im = mul_i2(m12-m11+q4);	// 0,b30
		t11 = t3 -rt;		t3 = t3 +rt;							m11 = m3 -rm;		m3 = m3 +rm;	//  3, 4 in -2b,2b+b30
		t12 = t4 -it;		t4 = t4 +it;							m12 = m4 -im;		m4 = m4 +im;	// 11,12 in -2b-b30,2b

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;				cmul_modq8(m19,m20, cm,q8-sm, &m19,&m20);
		rt =t27*s + t28*c;	it =t28*s - t27*c;						cmul_modq8(m27,m28, sm,q8-cm, &rm ,&im );
		t27 = t19-rt;		t19 = t19+rt;							m27 =qreduce(m19-rm+q4);	m19 =qreduce(m19+rm);	// +:   0,8q ==> 0,b
		t28 = t20-it;		t20 = t20+it;							m28 =qreduce(m20-im+q4);	m20 =qreduce(m20+im);	// -: -4q,4q ==> 0,b

		rt = t3 +t19;		it = t4 +t20;							rm = qreduce(m3 +m19+q3);		im = qreduce(m4 +m20+q3);	// rm,im in [-2b,2b+b30] + [0,b] = [-2b,3b+b30], add 3q and reduce
		t3 = t3 -t19;		t4 = t4 -t20;							m3 = qreduce(m3 -m19+q4);		m4 = qreduce(m4 -m20+q4);	// m3,m4 in [-2b,2b+b30] - [0,b] = [-3b,2b+b30], add 4q and reduce
		a[jt    ]=rt *c1 +it *s1;	a[jp    ]=it *c1 -rt *s1;		cmul_modq8(rm,im, a1,q8-b1, &rm,&im);
		a[jt+p8 ]=t3 *c9 +t4 *s9;	a[jp+p8 ]=t4 *c9 -t3 *s9;		cmul_modq8(m3,m4, a9,q8-b9, &m3,&m4);
																	b[jt    ] = qreduce(rm);	b[jp    ] = qreduce(im);
																	b[jt+p8 ] = qreduce(m3);	b[jp+p8 ] = qreduce(m4);
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
						// All 4 outs in [-2b-b30,2b] +- [-b,b] = [-3b-b30,3b]. Here we must +q5 prior to reducing since odds of
						// 3b+5q overflowing 64-bits much less than those of -3b-b30+4q not being enough to make result >= 0.
		rt  = t11+t28;		it  = t12-t27;							rm  = qreduce(m11+m28+q5);	im  = qreduce(m12-m27+q5);
		t11 = t11-t28;		t12 = t12+t27;							m11 = qreduce(m11-m28+q5);	m12 = qreduce(m12+m27+q5);
		a[jt+p4 ]=rt *c5 +it *s5;	a[jp+p4 ]=it *c5 -rt *s5;		cmul_modq8( rm, im, a5 ,q8-b5 ,  &rm, &im);
		a[jt+p12]=t11*c13+t12*s13;	a[jp+p12]=t12*c13-t11*s13;		cmul_modq8(m11,m12, a13,q8-b13, &m11,&m12);
																	b[jt+p4 ] = qreduce( rm);	b[jp+p4 ] = qreduce( im);
																	b[jt+p12] = qreduce(m11);	b[jp+p12] = qreduce(m12);
	/*...Block 4: t7,15,23,31 */
		jt = j1 + p3;		jp = j2 + p3;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;					rm = mul_i2(m15-m16+q4);im = mul_i2(m15+m16+q4);	// 0,b30
		t15 = t7 +rt;			t7 = t7 -rt;						m15 = m7 +rm;			m7 = m7 -rm;	//  7, 8 in -2b-b30,2b
		t16 = t8 +it;			t8 = t8 -it;						m16 = m8 +im;			m8 = m8 -im;	// 15,16 in -2b,2b+b30

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;				cmul_modq8(m23,m24, sm,q8-cm, &m23,&m24);
		rt =t31*c + t32*s;	it =t32*c - t31*s;						cmul_modq8(m31,m32, cm,q8-sm, &rm ,&im );
		t31 = t23+rt;			t23 = t23-rt;						m31 =qreduce(m23+rm);	m23 =qreduce(m23-rm+q4);// +:   0,8q ==> 0,b
		t32 = t24+it;			t24 = t24-it;						m32 =qreduce(m24+im);	m24 =qreduce(m24-im+q4);// -: -4q,4q ==> 0,b

		rt = t7 +t23;		it = t8 +t24;							rm = qreduce(m7 +m23+q4);	im = qreduce(m8 +m24+q4);	// +: [-2b-b30,2b] + [0,b] = [-2b-b30,3b], add 4q and reduce
		t7 = t7 -t23;		t8 = t8 -t24;							m7 = qreduce(m7 -m23+q5);	m8 = qreduce(m8 -m24+q5);	// -: [-2b-b30,2b] - [0,b] = [-3b-b30,2b], add 5q and reduce
		a[jt    ]=rt *c3 +it *s3;	a[jp    ]=it *c3 -rt *s3;		cmul_modq8(rm,im, a3 ,q8-b3 , &rm,&im);
		a[jt+p8 ]=t7 *c11+t8 *s11;	a[jp+p8 ]=t8 *c11-t7 *s11;		cmul_modq8(m7,m8, a11,q8-b11, &m7,&m8);
																	b[jt    ] = qreduce(rm);	b[jp    ] = qreduce(im);
																	b[jt+p8 ] = qreduce(m7);	b[jp+p8 ] = qreduce(m8);
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt  = t15+t32;		it  = t16-t31;							rm  = qreduce(m15+m32+q3);	im  = qreduce(m16-m31+q4);	// rm,m16 in [-2b,2b+b30] + [0,b] = [-2b,3b+b30], add 3q and reduce
		t15 = t15-t32;		t16 = t16+t31;							m15 = qreduce(m15-m32+q4);	m16 = qreduce(m16+m31+q3);	// im,m15 in [-2b,2b+b30] - [0,b] = [-3b,2b+b30], add 4q and reduce
		a[jt+p4 ]=rt *c7 +it *s7;	a[jp+p4 ]=it *c7 -rt *s7;		cmul_modq8( rm, im, a7 ,q8-b7 ,  &rm, &im);
		a[jt+p12]=t15*c15+t16*s15;	a[jp+p12]=t16*c15-t15*s15;		cmul_modq8(m15,m16, a15,q8-b15, &m15,&m16);
																	b[jt+p4 ] = qreduce( rm);	b[jp+p4 ] = qreduce( im);
																	b[jt+p12] = qreduce(m15);	b[jp+p12] = qreduce(m16);
							// Cost (aside from basic add/sub): 19 CMUL_MODQ8,  8 MUL_I2, 86 QREDUCE, 6 more reductions than DIF.

			/**********************************************/
	#else	// USE_FGT61 = False; Basic scalar-double mode:
			/**********************************************/

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/
	  #ifdef PFETCH_AGGRESSIVE
		addp = &a[j1];
	  #elif PFETCH
		prefetch_offset = ((j >> 1) & 3)*p1 + 4;	/* Cycle among p0, p1, p2 and p3. */
	  #endif

	/*...Block 1: */
		jt = j1;		jp = j2;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		t1 =a[jt ];	t2 =a[jp ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];
		t3 =t1 -rt;		t1 =t1 +rt;
		t4 =t2 -it;		t2 =t2 +it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =a[jt+p2 ];	t6 =a[jp+p2 ];	t5 =rt;
		rt =a[jt+p3 ];	it =a[jp+p3 ];
		t7 =t5 -rt;		t5 =t5 +rt;
		t8 =t6 -it;		t6 =t6 +it;

		rt =t5;	t5 =t1 -rt ;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it ;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8 ;	t3 =t3 +t8;
				t8 =t4 +rt ;	t4 =t4 -rt;

	/*...Block 2: */
		jt = j1 + p4;		jp = j2 + p4;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		rt =a[jt    ];	t10=a[jp    ];	t9 =rt;
		rt =a[jt+p1 ];	it =a[jp+p1 ];
		t11=t9 -rt;		t9 =t9 +rt;
		t12=t10-it;		t10=t10+it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =a[jt+p2 ];	t14=a[jp+p2 ];	t13=rt;
		rt =a[jt+p3 ];	it =a[jp+p3 ];
		t15=t13-rt;		t13=t13+rt;
		t16=t14-it;		t14=t14+it;

		rt =t13;	t13=t9 -rt ;	t9 =t9 +rt;
		it =t14;	t14=t10-it ;	t10=t10+it;

		rt =t15;	t15=t11-t16;	t11=t11+t16;
					t16=t12+rt ;	t12=t12-rt;

	/*...Block 3: */
		jt = j1 + p8;		jp = j2 + p8;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p4 ;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		rt =a[jt    ];	t18=a[jp    ];	t17=rt;
		rt =a[jt+p1 ];	it =a[jp+p1 ];
		t19=t17-rt;		t17=t17+rt;
		t20=t18-it;		t18=t18+it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =a[jt+p2 ];	t22=a[jp+p2 ];	t21=rt;
		rt =a[jt+p3 ];	it =a[jp+p3 ];
		t23=t21-rt;		t21=t21+rt;
		t24=t22-it;		t22=t22+it;

		rt =t21;	t21=t17-rt ;	t17=t17+rt;
		it =t22;	t22=t18-it ;	t18=t18+it;

		rt =t23;	t23=t19-t24;	t19=t19+t24;
					t24=t20+rt ;	t20=t20-rt;

	/*...Block 4: */
		jt = j1 + p12;		jp = j2 + p12;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	  #endif
		rt =a[jt    ];	t26=a[jp    ];	t25=rt;
		rt =a[jt+p1 ];	it =a[jp+p1 ];
		t27=t25-rt;		t25=t25+rt;
		t28=t26-it;		t26=t26+it;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =a[jt+p2 ];	t30=a[jp+p2 ];	t29=rt;
		rt =a[jt+p3 ];	it =a[jp+p3 ];
		t31=t29-rt;		t29=t29+rt;
		t32=t30-it;		t30=t30+it;

		rt =t29;	t29=t25-rt ;	t25=t25+rt;
		it =t30;	t30=t26-it ;	t26=t26+it;

		rt =t31;	t31=t27-t32;	t27=t27+t32;
					t32=t28+rt ;	t28=t28-rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(-i* 1*twopi/16) =	   ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =	   ( s,-c) (for inputs to transform block 2)
	!	1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =	   ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
	!	1, exp(-i* 3*twopi/16) =	   ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =	   (-c, s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[jt  +p0:15:1) are replaced by t0:30:2,
	!										   a[jp+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
		jt = j1;		jp = j2;
	  #ifdef PFETCH_AGGRESSIVE
		addr = addp + p8 ;
		prefetch_p_doubles(addr);
	  #endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[jt    ]=t1+t17;			a[jp    ]=t2+t18;
		t1	     =t1-t17;			t2		 =t2-t18;
		a[jt+p8 ]=t1 *c8 +t2 *s8;	a[jp+p8 ]=t2 *c8 -t1 *s8;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt	   =t9 +t26;		it		 =t10-t25;	/* mpy by E^-4 = -I is inlined here... */
		t9	   =t9 -t26;		t10		=t10+t25;
		a[jt+p4 ]=rt *c4 +it *s4;	a[jp+p4 ]=it *c4 -rt *s4;
		a[jt+p12]=t9 *c12+t10*s12;	a[jp+p12]=t10*c12-t9 *s12;

	/*...Block 3: t5,13,21,29 */
		jt = j1 + p2;		jp = j2 + p2;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;		/* twiddle mpy by E^4 =-I */
			t14=t6 +rt;	t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;	/* twiddle mpy by E^-2 */
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here... */
		t29=t21+rt;		t21=t21-rt;			/* ...and get E^-6 by flipping signs here. */
		t30=t22+it;		t22=t22-it;

		rt	   =t5 +t21;		it		 =t6 +t22;
		t5	   =t5 -t21;		t6		 =t6 -t22;
		a[jt    ]=rt *c2 +it *s2;	a[jp    ]=it *c2 -rt *s2;
		a[jt+p8 ]=t5 *c10+t6 *s10;	a[jp+p8 ]=t6 *c10-t5 *s10;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt	   =t13+t30;		it		 =t14-t29;	/* mpy by E^-4 = -I is inlined here... */
		t13	  =t13-t30;		t14		=t14+t29;
		a[jt+p4 ]=rt *c6 +it *s6;	a[jp+p4 ]=it *c6 -rt *s6;
		a[jt+p12]=t13*c14+t14*s14;	a[jp+p12]=t14*c14-t13*s14;

	/*...Block 2: t3,11,19,27 */
		jt = j1 + p1;		jp = j2 + p1;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;		/* twiddle mpy by E^-2 */
		t11=t3 -rt;		t3 =t3 +rt;
		t12=t4 -it;		t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;	/* twiddle mpy by E^-1 */
		rt =t27*s + t28*c;	it =t28*s - t27*c;		/* twiddle mpy by E^-3 */
		t27=t19-rt;		t19=t19+rt;
		t28=t20-it;		t20=t20+it;

		rt	   =t3 +t19;		it		 =t4 +t20;
		t3	   =t3 -t19;		t4		 =t4 -t20;
		a[jt    ]=rt *c1 +it *s1;	a[jp    ]=it *c1 -rt *s1;
		a[jt+p8 ]=t3 *c9 +t4 *s9;	a[jp+p8 ]=t4 *c9 -t3 *s9;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt	   =t11+t28;		it		 =t12-t27;	/* mpy by E^-4 = -I is inlined here... */
		t11	  =t11-t28;		t12		=t12+t27;
		a[jt+p4 ]=rt *c5 +it *s5;	a[jp+p4 ]=it *c5 -rt *s5;
		a[jt+p12]=t11*c13+t12*s13;	a[jp+p12]=t12*c13-t11*s13;

	/*...Block 4: t7,15,23,31 */
		jt = j1 + p3;		jp = j2 + p3;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here... */
		t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */
		t16=t8 +it;		t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;	/* twiddle mpy by E^-3 */
		rt =t31*c + t32*s;	it =t32*c - t31*s;		/* twiddle mpy by E^-1 = -E^-9... */
		t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here. */
		t32=t24+it;		t24=t24-it;

		rt	   =t7 +t23;		it		 =t8 +t24;
		t7	   =t7 -t23;		t8		 =t8 -t24;
		a[jt    ]=rt *c3 +it *s3;	a[jp    ]=it *c3 -rt *s3;
		a[jt+p8 ]=t7 *c11+t8 *s11;	a[jp+p8 ]=t8 *c11-t7 *s11;
	  #ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	  #endif
		rt	   =t15+t32;		it		 =t16-t31;	/* mpy by E^-4 = -I is inlined here... */
		t15	  =t15-t32;		t16		=t16+t31;
		a[jt+p4 ]=rt *c7 +it *s7;	a[jp+p4 ]=it *c7 -rt *s7;
		a[jt+p12]=t15*c15+t16*s15;	a[jp+p12]=t16*c15-t15*s15;

	#endif	// USE_FGT61 ?

  #endif	// USE_SCALAR_DFT_MACRO ?

#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

#undef RADIX
#undef PFETCH_DIST

