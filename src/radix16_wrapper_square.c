/*******************************************************************************
*                                                                              *
*   (C) 1997-2017 by Ernst W. Mayer.                                           *
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

#ifndef PFETCH_DIST
  #ifdef USE_AVX
	#define PFETCH_DIST	4096	// This seems to work best on my Haswell
  #else
	#define PFETCH_DIST	1024
  #endif
#endif

#ifdef USE_SSE2

	#include "sse2_macro.h"
	#include "radix16_utils_asm.h"

	#ifdef COMPILER_TYPE_GCC	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix16_wrapper_square_gcc32.h"

		#else

			#include "radix16_wrapper_square_gcc64.h"

		#endif

	#endif

#endif

/***************/

/*
Macro versions of these are in sse2_macro.h, since radix32_wrapper_square.c also needs to inline those:
SSE2 macros for this are in sse2_macro.h.
*/
void pair_square(double *x1, double *y1, double *x2, double *y2, double c, double s)
{
/*
!   Given complex scalars H[j] = (x1,y1) and H[N-j] = (x2,y2) along with complex exponential E = (c,s),
!   calculates I[j] = H[j]^2 + {1 + exp(2*pi*I*j/N)}*{H[j]-H~[N-j]}^2/4 and its complex conjugate I~,
!   returns the former in H[j] and the latter in H[N-j].
*/
	double rt0,rt1,rt2,rt3,it1,it2,it3;

	rt1 = *x1;	it1 = *y1;	rt2 = *x2;	it2 = *y2;

	// Calculate cross product terms:
	rt3 = rt1*rt2 + it1*it2; rt3 = rt3 + rt3;
	it3 = it1*rt2 - rt1*it2; it3 = it3 + it3;

	// Now calculate square terms and store back in the same temporaries:
	rt0 = (rt1 + it1)*(rt1 - it1); it1 = rt1*it1; it1 = it1 + it1; rt1 = rt0;
	rt0 = (rt2 + it2)*(rt2 - it2); it2 = rt2*it2; it2 = it2 + it2; rt2 = rt0;

	// Use that (H[j] - H~[N-j])^2 = H(j)^2 - 2*H(j)*H~(N-j) + H~(N-j)^2:
	rt3 = rt1 + rt2 - rt3;
	it3 = it1 - it2 - it3;
	rt0 = ((c + 1.0)*rt3 - s*it3)*0.25;
	it3 = (s*rt3 + (c + 1.0)*it3)*0.25;

	// And now complete and store the results:
	*x1 = (rt1 - rt0);
	*y1 = (it1 - it3);
	// N-j terms are as above, but with the replacements: rt1<-->rt2, it1<-->it2, it3|-->-it3:
	*x2 = (rt2 - rt0);
	*y2 = (it2 + it3);
}

/***************/

void radix16_wrapper_square(
	double a[],             int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[],
	int nradices_prim, int radix_prim[], int ws_i, int ws_j1, int ws_j2, int ws_j2_start, int ws_k, int ws_m, int ws_blocklen, int ws_blocklen_sum,
	int init_sse2, int thr_id
)
{
	const char func[] = "radix16_wrapper_square";
/*
!   NOTE: In the following commentary, N refers to the COMPLEX vector length (N2 in the code),
!   which is half the real vector length.
!
!...Acronyms:  DIT = Decimation In Time
!			   DIF = Decimation In Frequency
!			   FFT = Fast Fourier Transform, i.e. a discrete FT over the complex numbers
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
!	   j=0, ... , N-1.		(Complex forward wrapper step)
!
!   Then one does a pointwise squaring to obtain terms G[0,...,N-1],
!   and does an inverse complex wrapper prior to entering the inverse FFT:
!
!   I[j] = {( G[j]+G~[N-j] ) + I*exp(-pi*I*j/N)*( G[j]-G~[N-j] )}/2,
!
!	   j=0, ... , N-1.
!
!   By combining the 3 steps into one and doing some algebra, one
!   can get directly from the H's to the I's via
!
!   I[j] = H[j]^2 + {1 + exp(2*pi*I*j/N)}*{H[j]-H~[N-j]}^2/4,
!
!	   j=0, ... , N-1,
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

The scratch array (2nd input argument) is only needed for data table initializations, i.e. if first_entry = TRUE.
*/
	const int pfetch_dist = PFETCH_DIST;
#ifdef USE_SSE2
	const int stride = (int)RE_IM_STRIDE << 4;	// main-array loop stride = 32 for SSE2, 64 for AVX, 128 for AVX-512
#else
	const int stride = 32;	// In this particular routine, scalar mode has same stride as SSE2
#endif
	static int max_threads = 0;
	static int nsave = 0;
	static int *index = 0x0, *index_ptmp = 0x0;	/* N2/16-length Bit-reversal index array. */
	int *itmp = 0x0;
	int rdum,idum, j1pad,j2pad,kp,l,iroot,k1,k2;
	int i,j1,j2,j2_start,k,m,blocklen,blocklen_sum;
	/*int ndivrad0m1;*/
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)] */
	double rt,it, RT = 0.0,IT = 0.0;	// Caps used to highlight roots terms which need saving-for-squaring-step
	double re0,im0,re1,im1;
	double cA1,cA2,cA3,cA4,cA5,cA6,cA7,cA8,cA9,cA10,cA11,cA12,cA13,cA14,cA15,sA1,sA2,sA3,sA4,sA5,sA6,sA7,sA8,sA9,sA10,sA11,sA12,sA13,sA14,sA15;
	double cB1,cB2,cB3,cB4,cB5,cB6,cB7,cB8,cB9,cB10,cB11,cB12,cB13,cB14,cB15,sB1,sB2,sB3,sB4,sB5,sB6,sB7,sB8,sB9,sB10,sB11,sB12,sB13,sB14,sB15;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	double aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r
			,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i;
	double aj2p0r,aj2p1r,aj2p2r,aj2p3r,aj2p4r,aj2p5r,aj2p6r,aj2p7r,aj2p8r,aj2p9r,aj2p10r,aj2p11r,aj2p12r,aj2p13r,aj2p14r,aj2p15r
			,aj2p0i,aj2p1i,aj2p2i,aj2p3i,aj2p4i,aj2p5i,aj2p6i,aj2p7i,aj2p8i,aj2p9i,aj2p10i,aj2p11i,aj2p12i,aj2p13i,aj2p14i,aj2p15i;
#if PFETCH
	double *addr;
#endif

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error X86 SIMD code not supported for this compiler!
  #endif
	static uint32 *sm_arr = 0x0,*sm_ptr;	// Base-ptr to arrays of k1,k2-index-vectors used in SIMD roots-computation.
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0, *add1;	/* Addresses into array sections */
  #ifdef USE_AVX
	double *add2, *add3;
  #endif
  #ifdef USE_AVX512
	double *add4, *add5, *add6, *add7;
  #endif
	vec_dbl *c_tmp,*s_tmp;
	vec_dbl *tmp;

  #ifdef MULTITHREAD
	// Base addresses for discrete per-thread local stores ... 'r' for double-float data, 'i' for int:
	static vec_dbl *__r0;
	static uint32  *__i0;
	// In || mode, only above base-pointer (shared by all threads) is static:
	uint32  *k1_arr, *k2_arr;
	vec_dbl *cc0, *ss0, *isrt2,*one,*two,*forth, *tmp0, *tmp1, *tmp2, *tmp3
			,*r1,*r3,*r5,*r7,*r9,*r11,*r13,*r15,*r17,*r19,*r21,*r23,*r25,*r27,*r29,*r31
			,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15;
  #else
	// Same list of ptrs (sans base-address) as above, but now make them static:
	static uint32  *k1_arr, *k2_arr;
	static vec_dbl *cc0, *ss0, *isrt2,*one,*two,*forth, *tmp0, *tmp1, *tmp2, *tmp3
			,*r1,*r3,*r5,*r7,*r9,*r11,*r13,*r15,*r17,*r19,*r21,*r23,*r25,*r27,*r29,*r31
			,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15;
  #endif

#endif

/*...initialize things upon first entry */
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

	/* Here this variable is somewhat misnamed because it is used to init both non-SIMD and SIMD-specific data */
	if(init_sse2)	// Just check nonzero here, to allow the *value* of init_sse2 to store #threads
	{
		nsave = n;
		if(init_sse2 > max_threads)	// current SIMD local-alloc insufficient
		{
			ASSERT(HERE, thr_id == -1, "Init-mode call must be outside of any multithreading!");
			max_threads = init_sse2;
		#ifndef COMPILER_TYPE_GCC
			ASSERT(HERE, NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
		#endif
			ASSERT(HERE, max_threads >= NTHREADS, "Multithreading requires max_threads >= NTHREADS!");
	
		#ifdef USE_SSE2
			fprintf(stderr, "%s: pfetch_dist = %d\n",func,pfetch_dist);
			if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage
				free((void *)sm_arr);	sm_arr=0x0;
				free((void *)sc_arr);	sc_arr=0x0;
			}
			// Index vectors used in SIMD roots-computation.
			// The AVX512 compute-sincos-mults code needs 2 elements per complex-double-load, so use 10*RE_IM_STRIDE per array
			// to alloc storage here for all cases, even though that leaves upper array halves unused for sub-AVX512.
			sm_arr = ALLOC_INT(sm_arr, max_threads*20*RE_IM_STRIDE + 16);	if(!sm_arr){ sprintf(cbuf, "FATAL: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
			sm_ptr = ALIGN_INT(sm_arr);
			ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");
			// Twiddles-array: Need 0x47 slots for data, plus need to leave room to pad-align:
			sc_arr = ALLOC_VEC_DBL(sc_arr, 0x4c*max_threads);	ASSERT(HERE, sc_arr != 0,"FATAL: unable to allocate sc_arr!");
			sc_ptr = ALIGN_VEC_DBL(sc_arr);
			ASSERT(HERE, ((long)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		
			/* Use low 32 16-byte slots of sc_arr for temporaries, next 4 for const = 1/4 and nontrivial complex 16th roots,
			last 30 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
			*/
		  #ifdef MULTITHREAD
			__i0  = sm_ptr;
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x20;
			cc0   = sc_ptr + 0x21;
			ss0   = sc_ptr + 0x22;
			one   = sc_ptr + 0x43;
			two   = sc_ptr + 0x44;	// PAIR_SQUARE_4_AVX2 assumes two = forth-1
			forth = sc_ptr + 0x45;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);
				VEC_DBL_INIT(cc0  , c	);
				VEC_DBL_INIT(ss0  , s	);
				VEC_DBL_INIT(one  , 1.00);
				VEC_DBL_INIT(two  , 2.00);
				VEC_DBL_INIT(forth, 0.25);
				isrt2 += 0x4c;	/* Move on to next thread's local store */
				cc0   += 0x4c;
				ss0   += 0x4c;
				one   += 0x4c;
				two   += 0x4c;
				forth += 0x4c;
			}
		  #else
			k1_arr = sm_ptr;
			k2_arr = sm_ptr + 10*RE_IM_STRIDE;
			r1  = sc_ptr;			isrt2= r1  + 0x20;
			r3	= r1 + 0x02;		cc0	 = r1  + 0x21;	ss0	 = r1  + 0x22;
			r5	= r1 + 0x04;		c8	 = cc0 + 0x04;	// Insert pad-pair here to match offsets in radix16_dyadic_square
			r7	= r1 + 0x06;		c4	 = cc0 + 0x06;
			r9	= r1 + 0x08;		c12	 = cc0 + 0x08;
			r11	= r1 + 0x0a;		c2	 = cc0 + 0x0a;
			r13	= r1 + 0x0c;		c10	 = cc0 + 0x0c;
			r15	= r1 + 0x0e;		c6	 = cc0 + 0x0e;
			r17	= r1 + 0x10;		c14	 = cc0 + 0x10;
			r19	= r1 + 0x12;		c1	 = cc0 + 0x12;
			r21	= r1 + 0x14;		c9	 = cc0 + 0x14;
			r23	= r1 + 0x16;		c5	 = cc0 + 0x16;
			r25	= r1 + 0x18;		c13	 = cc0 + 0x18;
			r27	= r1 + 0x1a;		c3	 = cc0 + 0x1a;
			r29	= r1 + 0x1c;		c11	 = cc0 + 0x1c;
			r31	= r1 + 0x1e;		c7	 = cc0 + 0x1e;
									c15	 = cc0 + 0x20;
									one  = cc0 + 0x22;
									two  = cc0 + 0x23;	// PAIR_SQUARE_4_AVX2 assumes two = forth-1
									forth= cc0 + 0x24;
								tmp0 = cc0 + 0x02;	// stick these 2 tmp-ptrs into pad slots between ss0 and c8 above
								tmp1 = cc0 + 0x03;
									tmp2 = cc0 + 0x25;
									tmp3 = cc0 + 0x26;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
			VEC_DBL_INIT(cc0  , c	);
			VEC_DBL_INIT(ss0  , s	);
			VEC_DBL_INIT(one  , 1.00);
			VEC_DBL_INIT(two  , 2.00);
			VEC_DBL_INIT(forth, 0.25);
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
		itmp = (int *)&arr_scratch[N2/16];	/* Conservatively assume an int might be as long as 8 bytes here */
		for(i=0; i < N2/16; i++)
		{
			itmp[i]=i;
		}

		/*...then bit-reverse INDEX with respect to N/16. Radices are sent to BR routine
		in opposite order from which they are processed in forward FFT.
		*/
		bit_reverse_int(itmp, N2/16, nradices_prim-4, &radix_prim[nradices_prim-5], -1,(int *)arr_scratch);

	/*
	!...trig array for wrapper/squaring part is here. We will make use of the fact that the sincos
	!   data will be in bit-reversed order, together with the symmetries of the sineand cosine
	!   functions about pi/4 and pi/2, to cut our table length by a factor of 16 here. Here's how: the ordered sincos
	!   data (prior to bit-reversal of the data and reduction of table length) are a set of
	!   complex exponentials exp(I*j*2*pi/N), where j=0,...,N/2-1 and N is the complex FFT length.
	!   That is, the argument (complex polar angle) of the exponential is in [0,pi). Within the
	!   fused radix16/wrapper-square/radix16inv pass routine, we process 8 such complex exponentials
	!   at a time, together with 16 complex FFT array data. Since these 8 sincos data are bit-reversal-
	!   reordered, their angular arguments always have the form
	!
	!	x, x+pi/2, x+pi/4, x+3*pi/4, x+pi/8, x+5*pi/8, x+3*pi/8, x+7*pi/8,
	!
	!   which allows them to be written as four pairs, each of the form (angle, angle+pi/2):
	!
	!	(x0, x0+pi/2), (x1, x1+pi/2), (x2, x2+pi/2), (x3, x3+pi/2).
	!
	!   Given exp(I*angle) = (re,im), we use the complex symmetry  exp(I*(angle+pi/2)) = (-im,re)
	!   to get the second exponential of each pair. Since x1 = x0+pi/4, x2 = x0+pi/8 and x1 = x0+3*pi/8,
	!   we further let
	!		  ISRT2 :=	 1/sqrt(2) =    cos(pi/4) =    sin(pi/4),
	!		  (c,s) := exp(I*  pi/8) = (cos(  pi/8) , sin(  pi/8)),
	!	      and	  exp(I*3*pi/8) = (cos(3*pi/8) , sin(3*pi/8)) = (sin(pi/8) , cos(pi/8)) = (s,c)
	!
	!   Given only exp(I*x0) = (re,im), this allows us to get the remaining seven complex exponentials via
	!
	!		  exp(I*(x0       )) = ( re, im)
	!		  exp(I*(x0+  pi/2)) = I*exp(I*x0) = (-im, re)
	!
	!		  exp(I*(x1       )) = (re-im, re+im)*ISRT2
	!		  exp(I*(x1+  pi/2)) = I*exp(I*x1)
	!
	!		  exp(I*(x2       )) = ( re, im)*(c,s) = (re*c-im*s, re*s+im*c)
	!		  exp(I*(x2+  pi/2)) = I*exp(I*x2)
	!
	!		  exp(I*(x3       )) = ( re, im)*(s,c) = (re*s-im*c, re*c+im*s)
	!		  exp(I*(x3+  pi/2)) = I*exp(I*x3),
	!
	!   where (a,b) = a+I*b and (a,b)*(c,d) denotes a complex product, i.e. using both explicit-I* and complex-
	!   as-real-pair notation, (a,b)*(c,d) = (a+I*b)*(c+I*d) = (a*c - b*d) + I*(a*d+ b*c) = (a*c - b*d, a*d + b*c).
	!
	!   Since we already have the constants ISRT2, c and s in registers for the radix-16 transform, we need
	!   define no extra temporaries for the above. Further, by saving the four scalar products re*c, im*s, re*s, im*c,
	!   we can get both exp(I*x2) and exp(I*x3) using just four multiplies and fouradds, so our total cost for
	!   the eight complex sincos data is 2 loads from memory, 6 FMUL, 6 FADD and 4 negations, compared to
	!   16 loads from memory (and the accompanying register pressure) for the old version.
	!
	!   We further use the fact that for every two blocks of 16 FFT data processed together, the wrapper sincos
	!   data for the second (upper) block are just the reflection about pi/2 of those for the first block.
	!   This cuts the table length in half again.
	*/

	/*
	!...Allocate and initialize an index array containing N/16 indices to store the sequentially-rearranged
	!   FFT sincos data indices. (Only need N/16 of these since we need one base sincos index for every 16 complex data).
	!   We don't need a separate sincos array for the rea/complex wrapper phase, since this uses the same sincos datum
	!   as is used for the the first of each of the two blocks of 16 complex FFT data.
	*/
		if(index_ptmp != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)index_ptmp);	index_ptmp=0x0;
		}
		index_ptmp = ALLOC_INT(index_ptmp, N2/16);
		ASSERT(HERE, index_ptmp != 0,"FATAL: unable to allocate array INDEX!");
		index = ALIGN_INT(index_ptmp);
	/*
	!...Now rearrange FFT sincos indices using the main loop structure as a template.
	!   The first length-2 block is a little different, since there we process the 0-15 and 16-31 array
	!   elements separately.
	*/
		index[0]=itmp [0];
		index[1]=itmp [1];

		k1 =16;	/* init uncompressed wrapper sincos array index */
		k  =2;	/* init   compressed wrapper sincos array index */

		blocklen=16;
		blocklen_sum=16;
		j2_start=96;

		for(i = nradices_prim-6; i >= 0; i--)   /* Radices get processed in reverse order here as in forward FFT. */
		{
		  kp = k1 + blocklen;

		  for(m = 0; m < (blocklen-1)>>1; m += 8)	/* Since we now process TWO 16-element sets per loop execution, only execute the loop half as many times as before. */
		  {
	/*...grab the next 2 FFT sincos indices: one for the lower (j1, in the actual loop) data block, one for the upper (j2) data block... */

			index[k  ]=itmp[(k1  )>>3];
			index[k+1]=itmp[(kp-8)>>3];

			k1 = k1+ 8;
			kp = kp- 8;

			k  = k + 2;
		  }

		  k1 = k1 + (blocklen >> 1);

		  if(j2_start == n-32)break;

		  blocklen_sum = blocklen_sum + blocklen;
		  ASSERT(HERE, i != 0,"ERROR 10!");
		  blocklen = (radix_prim[i-1]-1)*blocklen_sum;

		  j2_start = j2_start+(blocklen<<2);
		}

		j1 = 0;

		/* Restore zeros here, to prevent any barfing due to interpretation of the above integer values as floats,
		in case at some future point we find it useful to be able to re-use a part of the main a-array for scratch: */
		for(i=0; i < N2/16; i++)
		{
			itmp[i]=0;
		}

		return;
	}	/* end of inits. */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
#ifdef MULTITHREAD
	ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
  #ifdef USE_SSE2
	k1_arr = __i0 + thr_id*20*RE_IM_STRIDE;
	k2_arr = k1_arr +      10*RE_IM_STRIDE;
	r1 = __r0 + thr_id*0x4c;cc0	= r1 + 0x21;
	r3	= r1 + 0x02;		ss0	= r1 + 0x22;
	r5	= r1 + 0x04;		c8	= r1 + 0x25;	// Insert a pad here to match offsets in radix16_dyadic_square
	r7	= r1 + 0x06;		c4	= r1 + 0x27;
	r9	= r1 + 0x08;		c12	= r1 + 0x29;
	r11	= r1 + 0x0a;		c2	= r1 + 0x2b;
	r13	= r1 + 0x0c;		c10	= r1 + 0x2d;
	r15	= r1 + 0x0e;		c6	= r1 + 0x2f;
	r17	= r1 + 0x10;		c14	= r1 + 0x31;
	r19	= r1 + 0x12;		c1	= r1 + 0x33;
	r21	= r1 + 0x14;		c9	= r1 + 0x35;
	r23	= r1 + 0x16;		c5	= r1 + 0x37;
	r25	= r1 + 0x18;		c13	= r1 + 0x39;
	r27	= r1 + 0x1a;		c3	= r1 + 0x3b;
	r29	= r1 + 0x1c;		c11	= r1 + 0x3d;
	r31	= r1 + 0x1e;		c7	= r1 + 0x3f;
	isrt2=r1 + 0x20;		c15	= r1 + 0x41;
							one   = r1 + 0x43;
							two   = r1 + 0x44;	// PAIR_SQUARE_4_AVX2 assumes two = forth-1
							forth = r1 + 0x45;
							tmp0 = r1 + 0x23;	/* stick these 2 tmp-ptrs into pad slots between ss0 and c8 above */
							tmp1 = r1 + 0x24;
							tmp2 = r1 + 0x46;
							tmp3 = r1 + 0x47;
  #endif
#endif
	/*...If a new runlength, should not get to this point: */
	ASSERT(HERE, n == nsave,"n != nsave");

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
!                                                                                                                               ------  ------
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
!   Floating index 2 (j2):  3, 7, 5,15,13,11, 9 (0 and 1 done separately)
!
!   Modular  index 1 (j3):  1, 2, 3, 4, 5, 6, 7
!   Modular  index 2 (j4): 15,14,13,12,11,10, 9 (0 and N/2 done separately)
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
!       do m = 1,blocklen
!
!         call pair_square(a(j1),a(j1+1),a(j2),a(j2+1),w(1,j),w(2,j))
!
!         j  = j +1
!         j1 = j1+4
!         j2 = j2-4
!       enddo
!       ...
!       blocklen = 2*blocklen
!       j2_start = j2_start+ishft(blocklen,2)
!       j2=j2_start
!       cycle
!   endloop
!
!   So the first 16 complex elements and their images are processed as follows:
!
!   exe0: process a(0,1) and a(2,3) separately;
!
!   exe1:
!       j=1
!       blocklen = 1
!       j1=4
!       j2_start=6
!
!         call pair_square(a( 4),a( 5),a( 6),a( 7),w(:,1))
!
!         j  = 2
!         j1 = 8
!         j2 = 2
!       enddo
!
!       blocklen = 2
!       j2_start = 14
!   end1
!
!   exe2:
!       j=2
!       blocklen = 2
!       j1=8
!       j2_start=14
!
!         call pair_square(a( 8),a( 9),a(14),a(15),w(:,2))
!         call pair_square(a(12),a(13),a(10),a(11),w(:,3))
!
!         j  = 4
!         j1 = 16
!         j2 = 6
!       enddo
!
!       blocklen = 4
!       j2_start = 30
!   end2
!
!   exe3:
!       j=4
!       blocklen = 4
!       j1=16
!       j2_start=30
!
!         call pair_square(a(16),a(17),a(30),a(31),w(:,4))
!         call pair_square(a(20),a(21),a(26),a(27),w(:,5))
!         call pair_square(a(24),a(25),a(22),a(23),w(:,6))
!         call pair_square(a(28),a(29),a(18),a(19),w(:,7))
!
!         j  = 8
!         j1 = 32
!         j2 = 14
!       enddo
!
!       blocklen = 8
!       j2_start = 62
!   end3
!
!   exe4:
!       j=8
!       blocklen = 8
!       j1=32
!       j2_start=62
!
!         call pair_square(a(32),a(33),a(62),a(63),w(:,8))
!         call pair_square(a(36),a(37),a(58),a(59),w(:,9))
!         call pair_square(a(40),a(41),a(54),a(55),w(:,10))
!         call pair_square(a(44),a(45),a(50),a(51),w(:,11))
!         call pair_square(a(48),a(49),a(46),a(47),w(:,12))
!         call pair_square(a(52),a(53),a(42),a(43),w(:,13))
!         call pair_square(a(56),a(57),a(38),a(39),w(:,14))
!         call pair_square(a(60),a(61),a(34),a(35),w(:,15))
!
!         j  = 16
!         j1 = 64
!         j2 = 30
!       enddo
!
!       blocklen = 16
!       j2_start = 126
!   end4
!
!   From here onward the blocklength will always be an integer multiple of 16, i.e. we can process each block using pairs of nonoverlapping
!   blocks of 16 complex data each, which is compatible to fusion with radix-16 pass routines. For example, the next (fifth) loop execution
!   looks like:
!
!   exe5:
!       j=16
!       blocklen = 16
!       j1=64
!       j2_start=126
!                                                                       ...and these 32 complex data can be processed as follows:
!         call pair_square(a( 64),a( 65),a(126),a(127),w(:,16))
!         call pair_square(a( 68),a( 69),a(122),a(123),w(:,17))         do radix-16 DIF pass (using DIF sincos  0-15) to get a( 64: 95) ;
!         call pair_square(a( 72),a( 73),a(118),a(119),w(:,18))         do radix-16 DIF pass (using DIF sincos 16-31) to get a( 96:127) ;
!         call pair_square(a( 76),a( 77),a(114),a(115),w(:,19))
!         call pair_square(a( 80),a( 81),a(110),a(111),w(:,20)) {need to see whether there's any nice pattern to the FFT sincos
!         call pair_square(a( 84),a( 85),a(106),a(107),w(:,21))  data here which would allow us to cut the storage of same}
!         call pair_square(a( 88),a( 89),a(102),a(103),w(:,22))
!         call pair_square(a( 92),a( 93),a( 98),a( 99),w(:,23))         combine the 32 resulting complex array data as shown at left    ;
!
!         call pair_square(a( 96),a( 97),a( 94),a( 95),w(:,24))
!         call pair_square(a(100),a(101),a( 90),a( 91),w(:,25))         do radix-16 DIT pass on a( 64: 95)      ;
!         call pair_square(a(104),a(105),a( 86),a( 87),w(:,26))         do radix-16 DIT pass on a( 96:127)      .
!         call pair_square(a(108),a(109),a( 82),a( 83),w(:,27))
!         call pair_square(a(112),a(113),a( 78),a( 79),w(:,28))
!         call pair_square(a(116),a(117),a( 74),a( 75),w(:,29))
!         call pair_square(a(120),a(121),a( 70),a( 71),w(:,30))
!         call pair_square(a(124),a(125),a( 66),a( 67),w(:,31))
!                                                                        The radix-16 passes share many register data, providing added savings.
!         j  = 32
!         j1 = 128
!         j2 = 62
!       enddo
!
!       blocklen = 32
!       j2_start = 254
!   end5
*/

	/*ndivrad0m1 = n/radix0 - 1;*/

/* Init the loop-control variables: */

	i            = ws_i           ;
	j1           = ws_j1          ;
	j2           = ws_j2          ;
	j2_start     = ws_j2_start    ;
	k            = ws_k           ;
	m            = ws_m           ;
	blocklen     = ws_blocklen    ;
	blocklen_sum = ws_blocklen_sum;

	/* If j1 == 0 we need to init the loop counters; otherwise, just jump
	   right in and pick up where we left off on the previous pair of blocks:
	*/
	if(j1 > 0) {
	//	fprintf(stderr,"j1,j2 = %d,%d: Jumping into loop!\n",j1,j2);
		goto jump_in;
	} else {
	#if 0
		rng_isaac_init(TRUE);
		double pow2_dmult = 1024.0*128.0;	// Restrict inputs to 18 bits, which in balanced-digit representation
											// means restricting multiplier of random-inputs-in-[-1,+1] below to 2^17
		int ipad, imax = MIN(1024,n);
		for(i = 0; i < imax; i += 16) {
			ipad = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
			// All the inits are w.r.to an un-SIMD-rearranged ...,re,im,re,im,... pattern:
		#ifdef USE_AVX512
			a[ipad+br16[ 0]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
			a[ipad+br16[ 1]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
			a[ipad+br16[ 2]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
			a[ipad+br16[ 3]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
			a[ipad+br16[ 4]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
			a[ipad+br16[ 5]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
			a[ipad+br16[ 6]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
			a[ipad+br16[ 7]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
			a[ipad+br16[ 8]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
			a[ipad+br16[ 9]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
			a[ipad+br16[10]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
			a[ipad+br16[11]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
			a[ipad+br16[12]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
			a[ipad+br16[13]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
			a[ipad+br16[14]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
			a[ipad+br16[15]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
		#elif defined(USE_AVX)
			a[ipad+br8[0]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
			a[ipad+br8[1]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
			a[ipad+br8[2]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
			a[ipad+br8[3]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
			a[ipad+br8[4]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
			a[ipad+br8[5]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
			a[ipad+br8[6]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
			a[ipad+br8[7]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
			a[ipad+br8[0]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
			a[ipad+br8[1]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
			a[ipad+br8[2]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
			a[ipad+br8[3]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
			a[ipad+br8[4]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
			a[ipad+br8[5]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
			a[ipad+br8[6]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
			a[ipad+br8[7]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
		#elif defined(USE_SSE2)
			a[ipad+br4[0]   ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
			a[ipad+br4[1]   ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
			a[ipad+br4[2]   ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
			a[ipad+br4[3]   ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
			a[ipad+br4[0]+4 ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
			a[ipad+br4[1]+4 ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
			a[ipad+br4[2]+4 ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
			a[ipad+br4[3]+4 ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
			a[ipad+br4[0]+8 ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
			a[ipad+br4[1]+8 ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
			a[ipad+br4[2]+8 ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
			a[ipad+br4[3]+8 ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
			a[ipad+br4[0]+12] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
			a[ipad+br4[1]+12] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
			a[ipad+br4[2]+12] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
			a[ipad+br4[3]+12] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
		#else
			#error Debug only enabled for SSE2 and above!
		#endif
		#ifdef USE_AVX512
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 0, a[ipad+br16[ 0]],i+ 0, a[ipad+ 0]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 1, a[ipad+br16[ 1]],i+ 1, a[ipad+ 1]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 2, a[ipad+br16[ 2]],i+ 2, a[ipad+ 2]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 3, a[ipad+br16[ 3]],i+ 3, a[ipad+ 3]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 4, a[ipad+br16[ 4]],i+ 4, a[ipad+ 4]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 5, a[ipad+br16[ 5]],i+ 5, a[ipad+ 5]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 6, a[ipad+br16[ 6]],i+ 6, a[ipad+ 6]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 7, a[ipad+br16[ 7]],i+ 7, a[ipad+ 7]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 8, a[ipad+br16[ 8]],i+ 8, a[ipad+ 8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+ 9, a[ipad+br16[ 9]],i+ 9, a[ipad+ 9]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+10, a[ipad+br16[10]],i+10, a[ipad+10]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+11, a[ipad+br16[11]],i+11, a[ipad+11]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+12, a[ipad+br16[12]],i+12, a[ipad+12]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+13, a[ipad+br16[13]],i+13, a[ipad+13]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+14, a[ipad+br16[14]],i+14, a[ipad+14]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+15, a[ipad+br16[15]],i+15, a[ipad+15]);
		#elif defined(USE_AVX)
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+0  ,a[ipad+br8[0]  ],i+0  ,a[ipad+0  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+1  ,a[ipad+br8[1]  ],i+1  ,a[ipad+1  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+2  ,a[ipad+br8[2]  ],i+2  ,a[ipad+2  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+3  ,a[ipad+br8[3]  ],i+3  ,a[ipad+3  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+4  ,a[ipad+br8[4]  ],i+4  ,a[ipad+4  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+5  ,a[ipad+br8[5]  ],i+5  ,a[ipad+5  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+6  ,a[ipad+br8[6]  ],i+6  ,a[ipad+6  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+7  ,a[ipad+br8[7]  ],i+7  ,a[ipad+7  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+0+8,a[ipad+br8[0]+8],i+0+8,a[ipad+0+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+1+8,a[ipad+br8[1]+8],i+1+8,a[ipad+1+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+2+8,a[ipad+br8[2]+8],i+2+8,a[ipad+2+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+3+8,a[ipad+br8[3]+8],i+3+8,a[ipad+3+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+4+8,a[ipad+br8[4]+8],i+4+8,a[ipad+4+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+5+8,a[ipad+br8[5]+8],i+5+8,a[ipad+5+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+6+8,a[ipad+br8[6]+8],i+6+8,a[ipad+6+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+7+8,a[ipad+br8[7]+8],i+7+8,a[ipad+7+8]);
		#elif defined(USE_SSE2)
		  if(0) {
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+0  ,a[ipad+br4[0]   ],i+0  ,a[ipad+0  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+1  ,a[ipad+br4[1]   ],i+1  ,a[ipad+1  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+2  ,a[ipad+br4[2]   ],i+2  ,a[ipad+2  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+3  ,a[ipad+br4[3]   ],i+3  ,a[ipad+3  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+4  ,a[ipad+br4[0]+4 ],i+4  ,a[ipad+4  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+5  ,a[ipad+br4[1]+4 ],i+5  ,a[ipad+5  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+6  ,a[ipad+br4[2]+4 ],i+6  ,a[ipad+6  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+7  ,a[ipad+br4[3]+4 ],i+7  ,a[ipad+7  ]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+0+8,a[ipad+br4[0]+8 ],i+0+8,a[ipad+0+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+1+8,a[ipad+br4[1]+8 ],i+1+8,a[ipad+1+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+2+8,a[ipad+br4[2]+8 ],i+2+8,a[ipad+2+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+3+8,a[ipad+br4[3]+8 ],i+3+8,a[ipad+3+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+4+8,a[ipad+br4[0]+12],i+4+8,a[ipad+4+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+5+8,a[ipad+br4[1]+12],i+5+8,a[ipad+5+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+6+8,a[ipad+br4[2]+12],i+6+8,a[ipad+6+8]);
			printf("A_in[%3d] = %20.10e; SIMD: A_in[%3d] = %20.10e\n",i+7+8,a[ipad+br4[3]+12],i+7+8,a[ipad+7+8]);
		  }
		#endif
		}
	#endif
	}

/*
!...All but the first two radix-16 blocks are done on Mr. Henry Ford's assembly line. From there onward the blocklength
!   will always be an integer multiple of 16, i.e. we can process each block using pairs of nonoverlapping blocks of 16
!   complex data each, which is compatible to fusion with radix-16 pass routines.
*/

for(i = nradices_prim-5; i >= 0; i-- )	/* Main loop: lower bound = nradices_prim-radix_now. */
{						/* Remember, radices get processed in reverse order here as in forward FFT. */
#ifdef USE_AVX512
	for(m = 0; m < (blocklen-1)>>1; m += 32) /* In AVX-512, process eight 16-complex-double datasets per loop execution, thus only execute the loop half as many times as for AVX case. */
#elif defined(USE_AVX)
	for(m = 0; m < (blocklen-1)>>1; m += 16) /* In AVX mode, process four 16-complex-double datasets per loop execution, thus only execute the loop half as many times as for scalar/SSE2 case. */
#else
	for(m = 0; m < (blocklen-1)>>1; m += 8) /* Both scalar and SSE2 modes process two 16-complex-double datasets per loop, only execute the loop half as many times as before. */
#endif
	{
		// This tells us when we've reached the end of the current data block:
		// Apr 2014: Must store intermediate product j1*radix0 in a 64-bit int to prevent overflow!
		if(j1 && ((uint64)j1*radix0)%n == 0)
		{
		//	fprintf(stderr,"(j1,j2 = %d,%d: j1 && j1*radix0 == 0 (mod n)) check hit: returning\n",j1,j2);
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

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA1 =rt;	sA1 =it;

	RT = rt;	IT = it;	// Save copies of these 2 for wrapper step

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 1);	/* 4*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA2 =rt;	sA2 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2);	/* 8*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA4 =rt;	sA4 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2) + iroot;	/* 13*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA8 =rt;	sA8 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA13=rt;	sA13=it;
		/* c3,5 */
		t1=cA1 *cA4 ;	t2=cA1 *sA4 ;	rt=sA1 *cA4 ;	it=sA1 *sA4;
		cA3 =t1 +it;	sA3 =t2 -rt;	cA5 =t1 -it;	sA5 =t2 +rt;

		/* c6,7,9,10 */
		t1=cA1 *cA8 ;	t2=cA1 *sA8 ;	rt=sA1 *cA8 ;	it=sA1 *sA8;
		cA7 =t1 +it;	sA7 =t2 -rt;	cA9 =t1 -it;	sA9 =t2 +rt;

		t1=cA2 *cA8 ;	t2=cA2 *sA8 ;	rt=sA2 *cA8 ;	it=sA2 *sA8;
		cA6 =t1 +it;	sA6 =t2 -rt;	cA10=t1 -it;	sA10=t2 +rt;

		/* c11,12,14,15 */
		t1=cA1 *cA13;	t2=cA1 *sA13;	rt=sA1 *cA13;	it=sA1 *sA13;
		cA12=t1 +it;	sA12=t2 -rt;	cA14=t1 -it;	sA14=t2 +rt;

		t1=cA2 *cA13;	t2=cA2 *sA13;	rt=sA2 *cA13;	it=sA2 *sA13;
		cA11=t1 +it;	sA11=t2 -rt;	cA15=t1 -it;	sA15=t2 +rt;

	/*************************************************************/
	/*                  2nd set of sincos:                       */
	/*************************************************************/
		iroot = index[k++];
		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB1 =rt;	sB1 =it;

	  if(j1 == 0) {   // The j1 = 0 case is special...need to overwrite above A-versions of saved trigs with B-data
		RT = rt;	IT = it;
	  }

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 1);	/* 4*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB2 =rt;	sB2 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2);	/* 8*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB4 =rt;	sB4 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2) + iroot;	/* 13*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB8 =rt;	sB8 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB13=rt;	sB13=it;
		/* c3,5 */
		t1=cB1 *cB4 ;	t2=cB1 *sB4 ;	rt=sB1 *cB4 ;	it=sB1 *sB4;
		cB3 =t1 +it;	sB3 =t2 -rt;	cB5 =t1 -it;	sB5 =t2 +rt;

		/* c6,7,9,10 */
		t1=cB1 *cB8 ;	t2=cB1 *sB8 ;	rt=sB1 *cB8 ;	it=sB1 *sB8;
		cB7 =t1 +it;	sB7 =t2 -rt;	cB9 =t1 -it;	sB9 =t2 +rt;

		t1=cB2 *cB8 ;	t2=cB2 *sB8 ;	rt=sB2 *cB8 ;	it=sB2 *sB8;
		cB6 =t1 +it;	sB6 =t2 -rt;	cB10=t1 -it;	sB10=t2 +rt;

		/* c11,12,14,15 */
		t1=cB1 *cB13;	t2=cB1 *sB13;	rt=sB1 *cB13;	it=sB1 *sB13;
		cB12=t1 +it;	sB12=t2 -rt;	cB14=t1 -it;	sB14=t2 +rt;

		t1=cB2 *cB13;	t2=cB2 *sB13;	rt=sB2 *cB13;	it=sB2 *sB13;
		cB11=t1 +it;	sB11=t2 -rt;	cB15=t1 -it;	sB15=t2 +rt;

	#else	// SIMD:

		/* Due to roots-locality considerations, roots (c,s)[1-15] - note no unused (c0,s0) pair here as in the
		Fermat-mod case -are offset w.r.to the thread-local ptr pair as

		             c 0  1  2  3  4  5  6  7  8  9  a  b  c  d  e  f
		(cc0,ss0) + 0x[-,10, 8,18, 4,14, c,1c, 2,12, a,1a, 6,16, e,1e].

		Here, due to the need to compute a new set of roots
		for each set of inputs, we use a streamlined sequence which computes only the [1,2,4,8,13]th roots with
		maximal accuracy (i.e. using 2-table-multiply), then generates the remaining ones from those. Thus the needed
		pointer offsets below are (cc0,ss0) + 0x[10, 8, 4, 2,16]:
		*/

	  #if defined(USE_SSE2) && !defined(USE_AVX)
	// In SSE2 mode need 2 sets of sincos - since we may use different indexing patterns for storing k1,2_arr data
	// in AVX+ mode, don't share code for first 2 sets with those wider SIMD cases:

		// 1st set:
		iroot = index[k++];
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
		iroot = index[k++];
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

	   #if OS_BITS == 32	// In 32-bit mode we only support SSE2 SIMD:

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX16_CALC_TWIDDLES_1_2_4_8_13(cc0,add0,add1,rt0,rt1);
		// Due to register paucity in 32-bit mode, use separate smaller asm macros to first compute the set of 5
		// 'anchor twiddles' using full 2-table complex multiplies, then the remaining 10, done in 5-pairs fashion:
		SSE2_CMUL_EXPO(c1,c4 ,c3 ,c5 );
		SSE2_CMUL_EXPO(c1,c8 ,c7 ,c9 );
		SSE2_CMUL_EXPO(c2,c8 ,c6 ,c10);
		SSE2_CMUL_EXPO(c1,c13,c12,c14);
		SSE2_CMUL_EXPO(c2,c13,c11,c15);

	   #else	// 64-bit SSE2:

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX16_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	   #endif	// 32-or-64-bit SSE2 ?

	  #elif !defined(USE_AVX512)	// AVX/AVX2:

	  // In AVX mode need 4 sets of sincos:
		// Need to explicitly 0 these for the first-few-blocks-done-via-scalar-code case in AVX mode to make sure the
		// unused-in-those-cases 3rd/4th-set indices stay nice and non-segfault-y:
		if(j1 <= 160) {
			memset(k1_arr, 0, 10*RE_IM_STRIDE*sizeof(uint32));
			memset(k2_arr, 0, 10*RE_IM_STRIDE*sizeof(uint32));
		}

		// 1st set:
		iroot = index[k++];
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
		iroot = index[k++];
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

	  if(j1 > 64)	// Sincos data for the initial done-in-scalar-mode data blocks get done in SSE2 mode, i.e. only using d0,d1
	  {
		// 3rd set:
		iroot = index[k++];
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
		iroot = index[k++];
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
	  }	// endif(j1 > 64)

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX16_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #elif defined(USE_AVX512)	// AVX512:

	  // In AVX-512 mode need 8 sets of sincos:
		// Need to explicitly 0 these for the first-few-blocks-done-via-scalar-code case in AVX mode to make sure the
		// unused-in-those-cases 3rd/4th-set indices stay nice and non-segfault-y:
		if(j1 <= 160) {
			memset(k1_arr, 0, 10*RE_IM_STRIDE*sizeof(uint32));
			memset(k2_arr, 0, 10*RE_IM_STRIDE*sizeof(uint32));
		}

	// *** To-do: Vectorize the index computations below! ***
		// 1st set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 0] = k1<<4;	k2_arr[ 0] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 8] = k1<<4;	k2_arr[ 8] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[16] = k1<<4;	k2_arr[16] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[24] = k1<<4;	k2_arr[24] = k2<<4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[32] = k1<<4;	k2_arr[32] = k2<<4;

		// 2nd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 1] = k1<<4;	k2_arr[ 1] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 9] = k1<<4;	k2_arr[ 9] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[17] = k1<<4;	k2_arr[17] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[25] = k1<<4;	k2_arr[25] = k2<<4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[33] = k1<<4;	k2_arr[33] = k2<<4;

	  if(j1 > 64) {	// Sincos data for the initial done-in-scalar-mode data blocks get done in SSE2 mode, i.e. only using d0,d1
		// 3rd set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 2] = k1<<4;	k2_arr[ 2] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[10] = k1<<4;	k2_arr[10] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[18] = k1<<4;	k2_arr[18] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[26] = k1<<4;	k2_arr[26] = k2<<4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[34] = k1<<4;	k2_arr[34] = k2<<4;

		// 4th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 3] = k1<<4;	k2_arr[ 3] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[11] = k1<<4;	k2_arr[11] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[19] = k1<<4;	k2_arr[19] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[27] = k1<<4;	k2_arr[27] = k2<<4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[35] = k1<<4;	k2_arr[35] = k2<<4;
	  }	// endif(j1 > 64)
	  if(j1 > 160) {	// Must skip these sincos-inits until shift to SIMD-mode DFT, to keep array index k from being incremented
		// 5th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 4] = k1<<4;	k2_arr[ 4] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[12] = k1<<4;	k2_arr[12] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[20] = k1<<4;	k2_arr[20] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[28] = k1<<4;	k2_arr[28] = k2<<4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[36] = k1<<4;	k2_arr[36] = k2<<4;

		// 6th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 5] = k1<<4;	k2_arr[ 5] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[13] = k1<<4;	k2_arr[13] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[21] = k1<<4;	k2_arr[21] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[29] = k1<<4;	k2_arr[29] = k2<<4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[37] = k1<<4;	k2_arr[37] = k2<<4;

		// 7th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 6] = k1<<4;	k2_arr[ 6] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[14] = k1<<4;	k2_arr[14] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[22] = k1<<4;	k2_arr[22] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[30] = k1<<4;	k2_arr[30] = k2<<4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[38] = k1<<4;	k2_arr[38] = k2<<4;

		// 8th set:
		iroot = index[k++];
		l = iroot;
		// Lshift all the vector k1,k2 indices to effect *= 16 needed for complex-double-strided access:
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[ 7] = k1<<4;	k2_arr[ 7] = k2<<4;
		l += iroot;			/* 2*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[15] = k1<<4;	k2_arr[15] = k2<<4;
		l += (iroot << 1);	/* 4*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[23] = k1<<4;	k2_arr[23] = k2<<4;
		l += (iroot << 2);	/* 8*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[31] = k1<<4;	k2_arr[31] = k2<<4;
		l += 5*iroot;		/* 13*iroot */
		k1=(l & NRTM1);	k2=(l >> NRT_BITS);	k1_arr[39] = k1<<4;	k2_arr[39] = k2<<4;
	  }	// endif(j1 > 160)

		// Stash head-of-array-ptrs in tmps to workaround GCC's "not directly addressable" macro arglist stupidity:
		add0 = (double *)k1_arr; add1 = (double *)k2_arr;	// Casts are only to get rid of compiler warnings
		SSE2_RADIX16_CALC_TWIDDLES_LOACC(cc0,add0,add1,rt0,rt1);

	  #endif	// SIMD mode?

	#endif	// SIMD ?

	#ifdef USE_SSE2	// Both SSE2 and AVX share this:

	  #ifdef USE_AVX512	// The generic pre-dyadic-square macro needs 8 main-array addresses in AVX512 mode

		// process 8 main-array blocks of [4 vec_dbl = 4 x 8 = 32 doubles] each in AVX512 mode, total = 32 vec_dbl = 16 vec_cmplx
		add0 = a + j1pad;
		add2 = add0 + 32;	// add2 = add0 + [32 doubles, equiv to 4 AVX-512 registers]
		add4 = add2 + 32;
		add6 = add4 + 32;
		add1 = a + j2pad;
		add3 = add1 - 32;	// Last 4 offsets run in descending order for Mers-mod
		add5 = add3 - 32;
		add7 = add5 - 32;

	  #elif defined(USE_AVX)	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  						// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:

		// process 4 main-array blocks of [8 vec_dbl = 8 x 4 = 32 doubles] each in AVX mode, total = 32 vec_dbl = 16 vec_cmplx
		add0 = a + j1pad;
		add1 = a + j2pad;
		add2 = add0 + 32;	// add2 = add0 + [32 doubles, equiv to 8 AVX registers]
		add3 = add1 - 32;	// Last 2 offsets run in descending order for Mers-mod

	  #else	// SSE2:

		add0 = a + j1pad;
		add1 = a + j2pad;

	  #endif

		// Use scalar code (with index offsets properly fiddled) for j1 == 0 case in SIMD mode:
	jump_new:
	  #ifdef USE_AVX
		// SSE2/AVX/AVX-512 need SIMD data processing to begin at j1 = 64,128,256, respectively. Here, 160
		// is the scalar-DFT-mode j1-value which is followed by j1 = 256 after block-index updating:
		if(j1 <= 160) {
		  if((j1 & 63) == 0) {	// Differentiate between j1 == 0 and 32 (mod 64) cases:
	  #else	// SSE2 begins SIMD-mode data processing at j1 = 128, which follows j1 = 64 after block-index-updating
	  		// (Could also use same value for AVX, but using same breakover for AVX,AVX-512 makes for easier debug of the latter.
		if(j1 <= 64) {
	  #endif
			cA1  = c1 ->d0;	sA1  = (c1 +1)->d0;		cB1  = c1 ->d1;	sB1  = (c1 +1)->d1;
			cA2  = c2 ->d0;	sA2  = (c2 +1)->d0;		cB2  = c2 ->d1;	sB2  = (c2 +1)->d1;
			cA3  = c3 ->d0;	sA3  = (c3 +1)->d0;		cB3  = c3 ->d1;	sB3  = (c3 +1)->d1;
			cA4  = c4 ->d0;	sA4  = (c4 +1)->d0;		cB4  = c4 ->d1;	sB4  = (c4 +1)->d1;
			cA5  = c5 ->d0;	sA5  = (c5 +1)->d0;		cB5  = c5 ->d1;	sB5  = (c5 +1)->d1;
			cA6  = c6 ->d0;	sA6  = (c6 +1)->d0;		cB6  = c6 ->d1;	sB6  = (c6 +1)->d1;
			cA7  = c7 ->d0;	sA7  = (c7 +1)->d0;		cB7  = c7 ->d1;	sB7  = (c7 +1)->d1;
			cA8  = c8 ->d0;	sA8  = (c8 +1)->d0;		cB8  = c8 ->d1;	sB8  = (c8 +1)->d1;
			cA9  = c9 ->d0;	sA9  = (c9 +1)->d0;		cB9  = c9 ->d1;	sB9  = (c9 +1)->d1;
			cA10 = c10->d0;	sA10 = (c10+1)->d0;		cB10 = c10->d1;	sB10 = (c10+1)->d1;
			cA11 = c11->d0;	sA11 = (c11+1)->d0;		cB11 = c11->d1;	sB11 = (c11+1)->d1;
			cA12 = c12->d0;	sA12 = (c12+1)->d0;		cB12 = c12->d1;	sB12 = (c12+1)->d1;
			cA13 = c13->d0;	sA13 = (c13+1)->d0;		cB13 = c13->d1;	sB13 = (c13+1)->d1;
			cA14 = c14->d0;	sA14 = (c14+1)->d0;		cB14 = c14->d1;	sB14 = (c14+1)->d1;
			cA15 = c15->d0;	sA15 = (c15+1)->d0;		cB15 = c15->d1;	sB15 = (c15+1)->d1;
	  #ifdef USE_AVX
		  } else {	// j1 == 32 (mod 64):
			cA1  = c1 ->d2;	sA1  = (c1 +1)->d2;		cB1  = c1 ->d3;	sB1  = (c1 +1)->d3;
			cA2  = c2 ->d2;	sA2  = (c2 +1)->d2;		cB2  = c2 ->d3;	sB2  = (c2 +1)->d3;
			cA3  = c3 ->d2;	sA3  = (c3 +1)->d2;		cB3  = c3 ->d3;	sB3  = (c3 +1)->d3;
			cA4  = c4 ->d2;	sA4  = (c4 +1)->d2;		cB4  = c4 ->d3;	sB4  = (c4 +1)->d3;
			cA5  = c5 ->d2;	sA5  = (c5 +1)->d2;		cB5  = c5 ->d3;	sB5  = (c5 +1)->d3;
			cA6  = c6 ->d2;	sA6  = (c6 +1)->d2;		cB6  = c6 ->d3;	sB6  = (c6 +1)->d3;
			cA7  = c7 ->d2;	sA7  = (c7 +1)->d2;		cB7  = c7 ->d3;	sB7  = (c7 +1)->d3;
			cA8  = c8 ->d2;	sA8  = (c8 +1)->d2;		cB8  = c8 ->d3;	sB8  = (c8 +1)->d3;
			cA9  = c9 ->d2;	sA9  = (c9 +1)->d2;		cB9  = c9 ->d3;	sB9  = (c9 +1)->d3;
			cA10 = c10->d2;	sA10 = (c10+1)->d2;		cB10 = c10->d3;	sB10 = (c10+1)->d3;
			cA11 = c11->d2;	sA11 = (c11+1)->d2;		cB11 = c11->d3;	sB11 = (c11+1)->d3;
			cA12 = c12->d2;	sA12 = (c12+1)->d2;		cB12 = c12->d3;	sB12 = (c12+1)->d3;
			cA13 = c13->d2;	sA13 = (c13+1)->d2;		cB13 = c13->d3;	sB13 = (c13+1)->d3;
			cA14 = c14->d2;	sA14 = (c14+1)->d2;		cB14 = c14->d3;	sB14 = (c14+1)->d3;
			cA15 = c15->d2;	sA15 = (c15+1)->d2;		cB15 = c15->d3;	sB15 = (c15+1)->d3;
		  }
	  #endif
		if(j1 == 0){ RT = cB1;	IT = sB1; } else { RT = cA1;	IT = sA1; }   /* Save for the wrapper Step ... The j1 = 0 case is special. */

		//************ NOTE: The above if(j1 <= 64 [or 160]) clause is still open! In SIMD mode we use it to wrap
		// all the scalar-mode FFT/dyad-mul/iFFT code, which allows scalar-double and SIMD builds to share code.
		//**********************************************************************************************************

	#endif	// USE_SSE2

	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/
	/*
	Data layout comparison:
	A-index:0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31
Scalar-dbl:	r0*	i0	r1	i1	r2	i2	r3	i3	r4*	i4	r5	i5	r6	i6	r7	i7	r8*	i8	r9	i9	r10	i10	r11	i11	r12*i12	r13	i13	r14	i14	r15	i15
	SSE:	r0*	r1	i0	i1	r2	r3	i2	i3	r4*	r5	i4	i5	r6	r7	i6	i7	r8	r9	i8	i9	r10	r11	i10	i11	r12	r13	i12	i13	r14	r15	i14	i15
	AVX:	r0*	r1	r2	r3	i0	i1	i2	i3	r4*	r5	r6	r7	i4	i5	i6	i7	r8	r9	r10	r11	i8	i9	i10	i11	r12	r13	r14	r15	i12	i13	i14	i15
	AVX-512:r0*	r1	r2	r3	r4*	r5	r6	r7	i0	i1	i2	i3	i4	i5	i6	i7	r8	r9	r10	r11	r12	r13	r14	r15	i8	i9	i10	i11	i12	i13	i14	i15

	Within the opening pass of four radix-4 sub-DFTs we combine complex inputs [0,8,4,12], [2,10,6,14], [1,9,5,13] and [3,11,7,15].
	Must fiddle the linear-idx jump w.r.to the +4,+2,+6 index strides of the non-SIMD layout, like so
	[ 0, 8, 4,12]: 4,12 need a-index -= 4 in avx-512 mode
	[ 2,10, 6,14]: 2,10 need a-index -= 2 in avx/avx-512 mode [i.e. avx-512 needs += 2 due to preceding -= 4]; 6,14 need idx -= 4 in avx-512
	[ 1, 9, 5,13]: 1, 9 need a-index -= 1 in sse2/avx/avx-512; 5,13 need idx -= 4 in avx-512
	[ 3,11, 7,15]: 3,11 need a-index -= 2 in avx/avx-512 mode [i.e. avx-512 needs += 2 due to preceding -= 4]; 7,15 need idx -= 4 in avx-512
	The #ifs below accomplish that:
	*/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;

	/*...Block 1: */
		t1 =a[rdum   ];							t2 =a[idum   ];						// z0
		rt =a[rdum+16]*cA8 -a[idum+16]*sA8 ;	it =a[idum+16]*cA8 +a[rdum+16]*sA8;	// z8
		t3 =t1 -rt;	t1 =t1 +rt;
		t4 =t2 -it;	t2 =t2 +it;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;	// r0 -> r4 is += 8 in scalar-double/SSE2/AVX modes, just += 4 in AVX-512 mode
	#endif
		t5 =a[rdum+8 ]*cA4 -a[idum+8 ]*sA4;		t6 =a[idum+8 ]*cA4 +a[rdum+8 ]*sA4;	// z4
		rt =a[rdum+24]*cA12-a[idum+24]*sA12;	it =a[idum+24]*cA12+a[rdum+24]*sA12;// z12
		t7 =t5 -rt;	t5 =t5 +rt;
		t8 =t6 -it;	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
				t8 =t4 -rt;	t4 =t4 +rt;

	/*...Block 2: */
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;	// r0 -> r2 just like in AVX mode, but AVX-512 needs added += 4 due to foregoing -= 4 used for z4,z12
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;	// r0 -> r2 is += 4 in scalar-double and SSE2 modes, just += 2 in AVX/AVX-512 mode
	#endif
		t9 =a[rdum+4 ]*cA2 -a[idum+4 ]*sA2 ;	t10=a[idum+4 ]*cA2 +a[rdum+4 ]*sA2 ;// z2
		rt =a[rdum+20]*cA10-a[idum+20]*sA10;	it =a[idum+20]*cA10+a[rdum+20]*sA10;// z10
		t11=t9 -rt;	t9 =t9 +rt;
		t12=t10-it;	t10=t10+it;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;	// r2 -> r6 is += 8 in scalar-double/SSE2/AVX modes, just += 4 in AVX-512 mode
	#endif
		t13=a[rdum+12]*cA6 -a[idum+12]*sA6 ;	t14=a[idum+12]*cA6 +a[rdum+12]*sA6 ;// z6
		rt =a[rdum+28]*cA14-a[idum+28]*sA14;	it =a[idum+28]*cA14+a[rdum+28]*sA14;// z14
		t15=t13-rt;	t13=t13+rt;
		t16=t14-it;	t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11+t16;t11=t11-t16;
					t16=t12-rt;	t12=t12+rt;

	/*...Block 3: */
	#ifdef USE_AVX512
		rdum += 5;	idum += 5;	// r2 -> r1 just like in AVX mode, but AVX-512 needs added += 4 due to foregoing -= 4 used for z6,14
	#elif defined(USE_AVX)
		++rdum;	++idum;			// r2 -> r1 is -= 2 in scalar-double mode, -= 3 in SSE2, -= 1 in AVX/AVX-512:
	#elif defined(USE_SSE2)
		--rdum;	--idum;
	#endif
		t17=a[rdum+2 ]*cA1 -a[idum+2 ]*sA1 ;	t18=a[idum+2 ]*cA1 +a[rdum+2 ]*sA1 ;// z1
		rt =a[rdum+18]*cA9 -a[idum+18]*sA9 ;	it =a[idum+18]*cA9 +a[rdum+18]*sA9 ;// z9
		t19=t17-rt;	t17=t17+rt;
		t20=t18-it;	t18=t18+it;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;	// r1 -> r5 is += 8 in scalar-double/SSE2/AVX modes, just += 4 in AVX-512 mode
	#endif
		t21=a[rdum+10]*cA5 -a[idum+10]*sA5 ;	t22=a[idum+10]*cA5 +a[rdum+10]*sA5 ;// z5
		rt =a[rdum+26]*cA13-a[idum+26]*sA13;	it =a[idum+26]*cA13+a[rdum+26]*sA13;// z13
		t23=t21-rt;	t21=t21+rt;
		t24=t22-it;	t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19+t24;t19=t19-t24;
					t24=t20-rt;	t20=t20+rt;

	/*...Block 4: */
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;	// r1 -> r3 just like in AVX mode, but AVX-512 needs added += 4 due to foregoing -= 4 used for z5,z13
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;	// r1 -> r3 is += 4 in scalar-double and SSE2, += 2 in AVX/AVX-512
	#endif
		t25=a[rdum+6 ]*cA3 -a[idum+6 ]*sA3 ;	t26=a[idum+6 ]*cA3 +a[rdum+6 ]*sA3 ;// z3
		rt =a[rdum+22]*cA11-a[idum+22]*sA11;	it =a[idum+22]*cA11+a[rdum+22]*sA11;// z11
		t27=t25-rt;	t25=t25+rt;
		t28=t26-it;	t26=t26+it;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;	// r3 -> r7 is += 8 in scalar-double/SSE2/AVX modes, just += 4 in AVX-512 mode
	#endif
		t29=a[rdum+14]*cA7 -a[idum+14]*sA7 ;	t30=a[idum+14]*cA7 +a[rdum+14]*sA7 ;// z7
		rt =a[rdum+30]*cA15-a[idum+30]*sA15;	it =a[idum+30]*cA15+a[rdum+30]*sA15;// z15
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
	!					   a[j1+1+p0:15:1] are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
		rt =t9;		t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		aj1p0r =t1+t17;	aj1p0i =t2+t18;
		aj1p1r =t1-t17;	aj1p1i =t2-t18;

		aj1p2r =t9 -t26;	aj1p2i =t10+t25;
		aj1p3r =t9 +t26;	aj1p3i =t10-t25;

	/*...Block 3: t5,13,21,29 */
		rt =t13;	t13=t5 +t14;t5 =t5 -t14;
					t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;
		rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		aj1p4r =t5+t21;	aj1p4i =t6+t22;
		aj1p5r =t5-t21;	aj1p5i =t6-t22;

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

		aj1p8r =t3+t19;	aj1p8i =t4+t20;
		aj1p9r =t3-t19;	aj1p9i =t4-t20;

		aj1p10r=t11-t28;	aj1p10i=t12+t27;
		aj1p11r=t11+t28;	aj1p11i=t12-t27;

	/*...Block 4: t7,15,23,31 */
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;
		rt =t31*c - t32*s;	it =t32*c + t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		aj1p12r=t7+t23;	aj1p12i=t8+t24;
		aj1p13r=t7-t23;	aj1p13i=t8-t24;

		aj1p14r=t15-t32;	aj1p14i=t16+t31;
		aj1p15r=t15+t32;	aj1p15i=t16-t31;

	/*************************************************************/
	/*                  2nd set of inputs:                       */
	/*************************************************************/

		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;

	/*...Block 1: */
		t1 =a[rdum   ];							t2 =a[idum   ];
		rt =a[rdum+16]*cB8 -a[idum+16]*sB8 ;	it =a[idum+16]*cB8 +a[rdum+16]*sB8;
		t3 =t1 -rt;	t1 =t1 +rt;
		t4 =t2 -it;	t2 =t2 +it;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t5 =a[rdum+8 ]*cB4 -a[idum+8 ]*sB4;		t6 =a[idum+8 ]*cB4 +a[rdum+8 ]*sB4;
		rt =a[rdum+24]*cB12-a[idum+24]*sB12;	it =a[idum+24]*cB12+a[rdum+24]*sB12;
		t7 =t5 -rt;	t5 =t5 +rt;
		t8 =t6 -it;	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
				t8 =t4 -rt;	t4 =t4 +rt;

	/*...Block 2: */
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
		t9 =a[rdum+4 ]*cB2 -a[idum+4 ]*sB2 ;	t10=a[idum+4 ]*cB2 +a[rdum+4 ]*sB2 ;
		rt =a[rdum+20]*cB10-a[idum+20]*sB10;	it =a[idum+20]*cB10+a[rdum+20]*sB10;
		t11=t9 -rt;	t9 =t9 +rt;
		t12=t10-it;	t10=t10+it;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t13=a[rdum+12]*cB6 -a[idum+12]*sB6 ;	t14=a[idum+12]*cB6 +a[rdum+12]*sB6 ;
		rt =a[rdum+28]*cB14-a[idum+28]*sB14;	it =a[idum+28]*cB14+a[rdum+28]*sB14;
		t15=t13-rt;	t13=t13+rt;
		t16=t14-it;	t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11+t16;t11=t11-t16;
					t16=t12-rt;	t12=t12+rt;

	/*...Block 3: */
	#ifdef USE_AVX512
		rdum += 5;	idum += 5;
	#elif defined(USE_AVX)
		++rdum;	++idum;
	#elif defined(USE_SSE2)
		--rdum;	--idum;
	#endif
		t17=a[rdum+2 ]*cB1 -a[idum+2 ]*sB1 ;	t18=a[idum+2 ]*cB1 +a[rdum+2 ]*sB1 ;
		rt =a[rdum+18]*cB9 -a[idum+18]*sB9 ;	it =a[idum+18]*cB9 +a[rdum+18]*sB9 ;
		t19=t17-rt;	t17=t17+rt;
		t20=t18-it;	t18=t18+it;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t21=a[rdum+10]*cB5 -a[idum+10]*sB5 ;	t22=a[idum+10]*cB5 +a[rdum+10]*sB5 ;
		rt =a[rdum+26]*cB13-a[idum+26]*sB13;	it =a[idum+26]*cB13+a[rdum+26]*sB13;
		t23=t21-rt;	t21=t21+rt;
		t24=t22-it;	t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19+t24;t19=t19-t24;
					t24=t20-rt;	t20=t20+rt;

	/*...Block 4: */
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
		t25=a[rdum+6 ]*cB3 -a[idum+6 ]*sB3 ;	t26=a[idum+6 ]*cB3 +a[rdum+6 ]*sB3 ;
		rt =a[rdum+22]*cB11-a[idum+22]*sB11;	it =a[idum+22]*cB11+a[rdum+22]*sB11;
		t27=t25-rt;	t25=t25+rt;
		t28=t26-it;	t26=t26+it;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		t29=a[rdum+14]*cB7 -a[idum+14]*sB7 ;	t30=a[idum+14]*cB7 +a[rdum+14]*sB7 ;
		rt =a[rdum+30]*cB15-a[idum+30]*sB15;	it =a[idum+30]*cB15+a[rdum+30]*sB15;
		t31=t29-rt;	t29=t29+rt;
		t32=t30-it;	t30=t30+it;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27+t32;t27=t27-t32;
					t32=t28-rt;	t28=t28+rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	*/
	/*...Block 1: t1,9,17,25 */
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		aj2p0r =t1+t17;	aj2p0i =t2+t18;
		aj2p1r =t1-t17;	aj2p1i =t2-t18;

		aj2p2r =t9 -t26;	aj2p2i =t10+t25;
		aj2p3r =t9 +t26;	aj2p3i =t10-t25;

	/*...Block 3: t5,13,21,29 */
		rt =t13;	t13=t5 +t14;	t5 =t5 -t14;
			t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;
		rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		aj2p4r =t5+t21;	aj2p4i =t6+t22;
		aj2p5r =t5-t21;	aj2p5i =t6-t22;

		aj2p6r =t13-t30;	aj2p6i =t14+t29;
		aj2p7r =t13+t30;	aj2p7i =t14-t29;

	/*...Block 2: t3,11,19,27 */
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;
		rt =t27*s - t28*c;	it =t28*s + t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		aj2p8r =t3+t19;	aj2p8i =t4+t20;
		aj2p9r =t3-t19;	aj2p9i =t4-t20;

		aj2p10r=t11-t28;	aj2p10i=t12+t27;
		aj2p11r=t11+t28;	aj2p11i=t12-t27;

	/*...Block 4: t7,15,23,31 */
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;
		rt =t31*c - t32*s;	it =t32*c + t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		aj2p12r=t7+t23;	aj2p12i=t8+t24;
		aj2p13r=t7-t23;	aj2p13i=t8-t24;

		aj2p14r=t15-t32;	aj2p14i=t16+t31;
		aj2p15r=t15+t32;	aj2p15i=t16-t31;

/*
!...send the pairs of complex elements which are to be combined and sincos temporaries needed for the squaring to a
!   small subroutine. The j1 = 0 case is again exceptional.
*/
	  if(j1 == 0)	/* NB: mustn't use re, im as temps in exe1-3 section, since those contain saved sincos data for exe4 block. */
	  {				// Jul 2015: Now moot since to avoid name-clashes in FGT61-mode, have renamed the doubles re,im -> RT,IT
		/*...j = 0 (for which I(0) = Re{H(0)}^2 + i*Im{H(0)}^2) is done separately... */
		rt = aj1p0r;
		aj1p0r = (rt + aj1p0i)*(rt + aj1p0i);
		aj1p0i = (rt - aj1p0i)*(rt - aj1p0i);
		rt = aj1p0r;
		aj1p0r = 0.5*(rt + aj1p0i);
		aj1p0i = 0.5*(rt - aj1p0i);
		/*
		!...as is j = N/2 (for which I(j) = H(j)^2). Note that under bit-reversal the N/2 element gets mapped into
		!   the second complex data slot, i.e. is adjacent to the starting element.
		*/
		rt = aj1p1r*aj1p1i;
		aj1p1r = (aj1p1r + aj1p1i)*(aj1p1r - aj1p1i);
		aj1p1i = rt + rt;

		pair_square(&aj1p2r ,&aj1p2i ,&aj1p3r ,&aj1p3i ,0.0,1.0);	/* exe1 */

		rt = it = ISRT2;
		pair_square(&aj1p4r ,&aj1p4i ,&aj1p7r ,&aj1p7i , rt, it);	/* exe2 */
		pair_square(&aj1p6r ,&aj1p6i ,&aj1p5r ,&aj1p5i , -it, rt);	/* exe2 */

		pair_square(&aj1p8r ,&aj1p8i ,&aj1p15r,&aj1p15i, c, s);		/* exe3 */
		pair_square(&aj1p10r,&aj1p10i,&aj1p13r,&aj1p13i, -s, c);	/* exe3 */

		pair_square(&aj1p12r,&aj1p12i,&aj1p11r,&aj1p11i, s, c);		/* exe3 */
		pair_square(&aj1p14r,&aj1p14i,&aj1p9r ,&aj1p9i , -c, s);	/* exe3 */

		/* exe4: */
		pair_square(&aj2p0r ,&aj2p0i ,&aj2p15r,&aj2p15i, RT, IT);
		pair_square(&aj2p2r ,&aj2p2i ,&aj2p13r,&aj2p13i, -IT, RT);

		rt = (RT - IT)*ISRT2;	it = (RT + IT)*ISRT2;
		pair_square(&aj2p4r ,&aj2p4i ,&aj2p11r,&aj2p11i, rt, it);
		pair_square(&aj2p6r ,&aj2p6i ,&aj2p9r ,&aj2p9i , -it, rt);

		t3 = RT*c;	t4 = IT*s;	t5 = RT*s;	t6 = IT*c;
		rt = t3 - t4;	it = t5 + t6;
		pair_square(&aj2p8r ,&aj2p8i ,&aj2p7r ,&aj2p7i , rt, it);
		pair_square(&aj2p10r,&aj2p10i,&aj2p5r ,&aj2p5i , -it, rt);

		rt = t5 - t6;	it = t3 + t4;
		pair_square(&aj2p12r,&aj2p12i,&aj2p3r ,&aj2p3i , rt, it);
		pair_square(&aj2p14r,&aj2p14i,&aj2p1r ,&aj2p1i , -it, rt);
	  }
	  else
	  {
		t1 = (RT - IT)*ISRT2;	t2 = (RT + IT)*ISRT2;

		re0 = RT*c;	im0 = IT*s;	re1 = RT*s;	im1 = IT*c;
		t3 = re0 - im0;	t4 = re1 + im1;
		t5 = re1 - im1;	t6 = re0 + im0;

		/*
		In SSE2 mode, the data are laid out in memory as

			t1:	aj1p0r,aj2p0r		t2: aj1p0i,aj2p0i
			t3:	aj1p1r,aj2p1r		t4: aj1p1i,aj2p1i
			t5:	aj1p2r,aj2p2r		t6: aj1p2i,aj2p2i
			t7:	aj1p3r,aj2p3r		t8: aj1p3i,aj2p3i
			t9:	aj1p4r,aj2p4r		t10:aj1p4i,aj2p4i
			t11:aj1p5r,aj2p5r		t12:aj1p5i,aj2p5i
			t13:aj1p6r,aj2p6r		t14:aj1p6i,aj2p6i
			t15:aj1p7r,aj2p7r		t16:aj1p7i,aj2p7i
			t17:aj1p8r,aj2p8r		t18:aj1p8i,aj2p8i
			t19:aj1p9r,aj2p9r		t20:aj1p9i,aj2p9i
			t21:aj1pAr,aj2pAr		t22:aj1pAi,aj2pAi
			t23:aj1pBr,aj2pBr		t24:aj1pBi,aj2pBi
			t25:aj1pCr,aj2pCr		t26:aj1pCi,aj2pCi
			t27:aj1pDr,aj2pDr		t28:aj1pDi,aj2pDi
			t29:aj1pEr,aj2pEr		t30:aj1pEi,aj2pEi
			t31:aj1pFr,aj2pFr		t32:aj1pFi,aj2pFi .

		The modified call sequence below takes advantage of that, by processing data which are in 8 XMM registers in-place.
		*/
		/* j1[ 0, 2,13,15] combined with j2[15,13, 2, 0]:
		PAIR_SQUARE2A(aj1p0r ,aj1p0i ,aj2p15r,aj2p15i,aj1p2r ,aj1p2i ,aj2p13r,aj2p13i, RT, IT);
		PAIR_SQUARE2B(aj2p2r ,aj2p2i ,aj1p13r,aj1p13i,aj2p0r ,aj2p0i ,aj1p15r,aj1p15i, t5, t6);
		*/
		PAIR_SQUARE_4(aj1p0r ,aj1p0i ,aj2p15r,aj2p15i,aj1p2r ,aj1p2i ,aj2p13r,aj2p13i, RT, IT
					, aj2p2r ,aj2p2i ,aj1p13r,aj1p13i,aj2p0r ,aj2p0i ,aj1p15r,aj1p15i, t5, t6);

		/* j1[ 4, 6, 9,11] combined with j2[11, 9, 6, 4]:
		PAIR_SQUARE2A(aj1p4r ,aj1p4i ,aj2p11r,aj2p11i,aj1p6r ,aj1p6i ,aj2p9r ,aj2p9i , t1, t2);
		PAIR_SQUARE2B(aj2p6r ,aj2p6i ,aj1p9r ,aj1p9i ,aj2p4r ,aj2p4i ,aj1p11r,aj1p11i, t3, t4);
		*/
		PAIR_SQUARE_4(aj1p4r ,aj1p4i ,aj2p11r,aj2p11i,aj1p6r ,aj1p6i ,aj2p9r ,aj2p9i , t1, t2
					, aj2p6r ,aj2p6i ,aj1p9r ,aj1p9i ,aj2p4r ,aj2p4i ,aj1p11r,aj1p11i, t3, t4);

		/* j1[ 8,10, 5, 7] combined with j2[ 7, 5,10, 8]:
		PAIR_SQUARE2A(aj1p8r ,aj1p8i ,aj2p7r ,aj2p7i ,aj1p10r,aj1p10i,aj2p5r ,aj2p5i , t3, t4);
		PAIR_SQUARE2B(aj2p10r,aj2p10i,aj1p5r ,aj1p5i ,aj2p8r ,aj2p8i ,aj1p7r ,aj1p7i , t1, t2);
		*/
		PAIR_SQUARE_4(aj1p8r ,aj1p8i ,aj2p7r ,aj2p7i ,aj1p10r,aj1p10i,aj2p5r ,aj2p5i , t3, t4
					, aj2p10r,aj2p10i,aj1p5r ,aj1p5i ,aj2p8r ,aj2p8i ,aj1p7r ,aj1p7i , t1, t2);

		/* j1[12,14, 1, 3] combined with j2[ 3, 1,14,12]:
		PAIR_SQUARE2A(aj1p12r,aj1p12i,aj2p3r ,aj2p3i ,aj1p14r,aj1p14i,aj2p1r ,aj2p1i , t5, t6);
		PAIR_SQUARE2B(aj2p14r,aj2p14i,aj1p1r ,aj1p1i ,aj2p12r,aj2p12i,aj1p3r ,aj1p3i , RT, IT);
		*/
		PAIR_SQUARE_4(aj1p12r,aj1p12i,aj2p3r ,aj2p3i ,aj1p14r,aj1p14i,aj2p1r ,aj2p1i , t5, t6
					, aj2p14r,aj2p14i,aj1p1r ,aj1p1i ,aj2p12r,aj2p12i,aj1p3r ,aj1p3i , RT, IT);
	  }

/*...And do an inverse DIT radix-16 pass on the squared-data blocks. */

	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;
	#if PFETCH
		addr = &a[rdum+32];
	#endif
	/*   gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 IDIT transforms... */

	/*...Block 1: */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		t3 =aj1p0r -aj1p1r;	t4 =aj1p0i -aj1p1i;
		t1 =aj1p0r +aj1p1r;	t2 =aj1p0i +aj1p1i;

		t7 =aj1p2r -aj1p3r;	t8 =aj1p2i -aj1p3i;
		t5 =aj1p2r +aj1p3r;	t6 =aj1p2i +aj1p3i;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
				t8 =t4 +rt;	t4 =t4 -rt;

	/*...Block 2: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		t11=aj1p4r -aj1p5r;	t12=aj1p4i-aj1p5i;
		t9 =aj1p4r +aj1p5r;	t10=aj1p4i+aj1p5i;

		t15=aj1p6r-aj1p7r;	t16=aj1p6i-aj1p7i;
		t13=aj1p6r+aj1p7r;	t14=aj1p6i+aj1p7i;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11-t16;t11=t11+t16;
					t16=t12+rt;	t12=t12-rt;

	/*...Block 3: */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		t19=aj1p8r-aj1p9r;	t20=aj1p8i-aj1p9i;
		t17=aj1p8r+aj1p9r;	t18=aj1p8i+aj1p9i;

		t23=aj1p10r-aj1p11r;	t24=aj1p10i-aj1p11i;
		t21=aj1p10r+aj1p11r;	t22=aj1p10i+aj1p11i;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19-t24;t19=t19+t24;
					t24=t20+rt;	t20=t20-rt;

	/*...Block 4: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		t27=aj1p12r-aj1p13r;	t28=aj1p12i-aj1p13i;
		t25=aj1p12r+aj1p13r;	t26=aj1p12i+aj1p13i;

		t31=aj1p14r-aj1p15r;	t32=aj1p14i-aj1p15i;
		t29=aj1p14r+aj1p15r;	t30=aj1p14i+aj1p15i;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27-t32;t27=t27+t32;
					t32=t28+rt;	t28=t28-rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
	!	1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
	!	1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1  +p0:15:1] are replaced by t0:30:2,
	!					   a[j1+1+p0:15:1] are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[rdum   ]=t1+t17;				a[idum   ]=t2+t18;
		t1	=t1-t17;					t2	=t2-t18;
		a[rdum+16]=t1 *cA8 +t2 *sA8 ;	a[idum+16]=t2 *cA8 -t1 *sA8;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		rt	=t9 +t26;					it	=t10-t25;
		t9	=t9 -t26;					t10	=t10+t25;
		a[rdum+8 ]=rt *cA4 +it *sA4 ;	a[idum+8 ]=it *cA4 -rt *sA4;
		a[rdum+24]=t9 *cA12+t10*sA12;	a[idum+24]=t10*cA12-t9 *sA12;

/*...Block 3: t5,13,21,29 */
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;
					t14=t6 +rt;		t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		rt	=t5 +t21;					it	=t6 +t22;
		t5	=t5 -t21;					t6	=t6 -t22;
		a[rdum+4 ]=rt *cA2 +it *sA2 ;	a[idum+4 ]=it *cA2 -rt *sA2 ;
		a[rdum+20]=t5 *cA10+t6 *sA10;	a[idum+20]=t6 *cA10-t5 *sA10;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		rt	=t13+t30;					it	=t14-t29;
		t13	=t13-t30;					t14	=t14+t29;
		a[rdum+12]=rt *cA6 +it *sA6 ;	a[idum+12]=it *cA6 -rt *sA6 ;
		a[rdum+28]=t13*cA14+t14*sA14;	a[idum+28]=t14*cA14-t13*sA14;

/*...Block 2: t3,11,19,27 */
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
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;
		rt =t27*s + t28*c;	it =t28*s - t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		rt	=t3 +t19;					it	=t4 +t20;
		t3	=t3 -t19;					t4	=t4 -t20;
		a[rdum+2 ]=rt *cA1 +it *sA1 ;	a[idum+2 ]=it *cA1 -rt *sA1 ;
		a[rdum+18]=t3 *cA9 +t4 *sA9 ;	a[idum+18]=t4 *cA9 -t3 *sA9 ;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		rt	=t11+t28;					it	=t12-t27;
		t11	=t11-t28;					t12	=t12+t27;
		a[rdum+10]=rt *cA5 +it *sA5 ;	a[idum+10]=it *cA5 -rt *sA5 ;
		a[rdum+26]=t11*cA13+t12*sA13;	a[idum+26]=t12*cA13-t11*sA13;

/*...Block 4: t7,15,23,31 */
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;
		rt =t31*c + t32*s;	it =t32*c - t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		rt	=t7 +t23;					it	=t8 +t24;
		t7	=t7 -t23;					t8	=t8 -t24;
		a[rdum+6 ]=rt *cA3 +it *sA3 ;	a[idum+6 ]=it *cA3 -rt *sA3 ;
		a[rdum+22]=t7 *cA11+t8 *sA11;	a[idum+22]=t8 *cA11-t7 *sA11;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		rt	=t15+t32;					it	=t16-t31;
		t15	=t15-t32;					t16	=t16+t31;
		a[rdum+14]=rt *cA7 +it *sA7 ;	a[idum+14]=it *cA7 -rt *sA7 ;
		a[rdum+30]=t15*cA15+t16*sA15;	a[idum+30]=t16*cA15-t15*sA15;

	/*************************************************************/
	/*                  2nd set of inputs:                       */
	/*************************************************************/

		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;
	#if PFETCH
		addr = &a[rdum-32];
	#endif
	/*...Block 1: */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		t3 =aj2p0r -aj2p1r;	t4 =aj2p0i -aj2p1i;
		t1 =aj2p0r +aj2p1r;	t2 =aj2p0i +aj2p1i;

		t7 =aj2p2r -aj2p3r;	t8 =aj2p2i -aj2p3i;
		t5 =aj2p2r +aj2p3r;	t6 =aj2p2i +aj2p3i;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
			t8 =t4 +rt;	t4 =t4 -rt;

	/*...Block 2: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		t11=aj2p4r -aj2p5r;	t12=aj2p4i-aj2p5i;
		t9 =aj2p4r +aj2p5r;	t10=aj2p4i+aj2p5i;

		t15=aj2p6r-aj2p7r;	t16=aj2p6i-aj2p7i;
		t13=aj2p6r+aj2p7r;	t14=aj2p6i+aj2p7i;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11-t16;	t11=t11+t16;
			t16=t12+rt;	t12=t12-rt;

	/*...Block 3: */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		t19=aj2p8r-aj2p9r;	t20=aj2p8i-aj2p9i;
		t17=aj2p8r+aj2p9r;	t18=aj2p8i+aj2p9i;

		t23=aj2p10r-aj2p11r;	t24=aj2p10i-aj2p11i;
		t21=aj2p10r+aj2p11r;	t22=aj2p10i+aj2p11i;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19-t24;	t19=t19+t24;
			t24=t20+rt;	t20=t20-rt;

	/*...Block 4: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		t27=aj2p12r-aj2p13r;	t28=aj2p12i-aj2p13i;
		t25=aj2p12r+aj2p13r;	t26=aj2p12i+aj2p13i;

		t31=aj2p14r-aj2p15r;	t32=aj2p14i-aj2p15i;
		t29=aj2p14r+aj2p15r;	t30=aj2p14i+aj2p15i;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27-t32;	t27=t27+t32;
			t32=t28+rt;	t28=t28-rt;

	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	*/
	/*...Block 1: t1,9,17,25 */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[rdum   ]=t1+t17;				a[idum   ]=t2+t18;
		t1	=t1-t17;					t2	=t2-t18;
		a[rdum+16]=t1 *cB8 +t2 *sB8 ;	a[idum+16]=t2 *cB8 -t1 *sB8;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		rt	=t9 +t26;					it	=t10-t25;
		t9	=t9 -t26;					t10	=t10+t25;
		a[rdum+8 ]=rt *cB4 +it *sB4;	a[idum+8 ]=it *cB4 -rt *sB4;
		a[rdum+24]=t9 *cB12+t10*sB12;	a[idum+24]=t10*cB12-t9 *sB12;

	/*...Block 3: t5,13,21,29 */
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;
			t14=t6 +rt;	t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		rt	=t5 +t21;					it	=t6 +t22;
		t5	=t5 -t21;					t6	=t6 -t22;
		a[rdum+4 ]=rt *cB2 +it *sB2 ;	a[idum+4 ]=it *cB2 -rt *sB2 ;
		a[rdum+20]=t5 *cB10+t6 *sB10;	a[idum+20]=t6 *cB10-t5 *sB10;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		rt	=t13+t30;					it	=t14-t29;
		t13	=t13-t30;					t14	=t14+t29;
		a[rdum+12]=rt *cB6 +it *sB6 ;	a[idum+12]=it *cB6 -rt *sB6 ;
		a[rdum+28]=t13*cB14+t14*sB14;	a[idum+28]=t14*cB14-t13*sB14;

	/*...Block 2: t3,11,19,27 */
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
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;
		rt =t27*s + t28*c;	it =t28*s - t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		rt	=t3 +t19;					it	=t4 +t20;
		t3	=t3 -t19;					t4	=t4 -t20;
		a[rdum+2 ]=rt *cB1 +it *sB1 ;	a[idum+2 ]=it *cB1 -rt *sB1 ;
		a[rdum+18]=t3 *cB9 +t4 *sB9 ;	a[idum+18]=t4 *cB9 -t3 *sB9 ;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		rt	=t11+t28;					it	=t12-t27;
		t11	=t11-t28;					t12	=t12+t27;
		a[rdum+10]=rt *cB5 +it *sB5 ;	a[idum+10]=it *cB5 -rt *sB5 ;
		a[rdum+26]=t11*cB13+t12*sB13;	a[idum+26]=t12*cB13-t11*sB13;

	/*...Block 4: t7,15,23,31 */
	#ifdef USE_AVX512
		rdum += 2;	idum += 2;
	#elif defined(USE_AVX)
		rdum -= 2;	idum -= 2;
	#endif
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;
		rt =t31*c + t32*s;	it =t32*c - t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		rt	=t7 +t23;					it	=t8 +t24;
		t7	=t7 -t23;					t8	=t8 -t24;
		a[rdum+6 ]=rt *cB3 +it *sB3 ;	a[idum+6 ]=it *cB3 -rt *sB3 ;
		a[rdum+22]=t7 *cB11+t8 *sB11;	a[idum+22]=t8 *cB11-t7 *sB11;
	#ifdef USE_AVX512
		rdum -= 4;	idum -= 4;
	#endif
		rt	=t15+t32;					it	=t16-t31;
		t15	=t15-t32;					t16	=t16+t31;
		a[rdum+14]=rt *cB7 +it *sB7 ;	a[idum+14]=it *cB7 -rt *sB7 ;
		a[rdum+30]=t15*cB15+t16*sB15;	a[idum+30]=t16*cB15-t15*sB15;

#ifdef USE_SSE2

	} else {	// (j1 == 0) = false:

		RT = c1->d0;	IT = (c1 +1)->d0;	// Save copies of these 2 for wrapper step

	/* In SSE2 mode, the input data are arranged in memory like so, where we view things in 16-byte chunks:

		&a[j1]	a0.re,a1.re		&a[j2]	b0.re,b1.re
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
		+0x0f0	aE.im,aF.im		+0x0f0	bE.im,bF.im

	In AVX mode, the input data are arranged in memory like so, where we view things in 32-byte chunks:

		&a[j1]	a0.re,a1.re,a2.re,a3.re		&a[j2]	b0.re,b1.re,b2.re,b3.re	 	&a[j3]	c0.re,c1.re,c2.re,c3.re	 	&a[j4]	d0.re,d1.re,d2.re,d3.re
		+0x020	a0.im,a1.im,a2.im,a3.im		+0x020	b0.im,b1.im,b2.im,b3.im		+0x020	c0.im,c1.im,c2.im,c3.im		+0x020	d0.im,d1.im,d2.im,d3.im
		+0x040	a4.re,a5.re,a6.re,a7.re		+0x040	b4.re,b5.re,b6.re,b7.re		+0x040	c4.re,c5.re,c6.re,c7.re		+0x040	d4.re,d5.re,d6.re,d7.re
		+0x060	a4.im,a5.im,a6.im,a7.im		+0x060	b4.im,b5.im,b6.im,b7.im		+0x060	c4.im,c5.im,c6.im,c7.im		+0x060	d4.im,d5.im,d6.im,d7.im
		+0x080	a8.re,a9.re,aA.re,aB.re		+0x080	b8.re,b9.re,bA.re,bB.re		+0x080	c8.re,c9.re,cA.re,cB.re		+0x080	d8.re,d9.re,dA.re,dB.re
		+0x0a0	a8.im,a9.im,aA.im,aB.im		+0x0a0	b8.im,b9.im,bA.im,bB.im		+0x0a0	c8.im,c9.im,cA.im,cB.im		+0x0a0	d8.im,d9.im,dA.im,dB.im
		+0x0c0	aC.re,aD.re,aE.re,aF.re		+0x0c0	bC.re,bD.re,bE.re,bF.re		+0x0c0	cC.re,cD.re,cE.re,cF.re		+0x0c0	dC.re,dD.re,dE.re,dF.re
		+0x0e0	aC.im,aD.im,aE.im,aF.im		+0x0e0	bC.im,bD.im,bE.im,bF.im		+0x0e0	cC.im,cD.im,cE.im,cF.im		+0x0e0	dC.im,dD.im,dE.im,dF.im
*/
/*
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
	#ifdef USE_AVX512		// process 8 main-array blocks of [4 vec_dbl = 4 x 8 = 32 doubles] each, total = 32 vec_dbl = 16 vec_cmplx

		SSE2_RADIX16_WRAPPER_DIF(add0,add1,add2,add3,add4,add5,add6,add7
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,pfetch_dist)

	#elif defined(USE_AVX)	// process 4 main-array blocks of [8 vec_dbl = 8 x 4 = 32 doubles] each, total = 32 vec_dbl = 16 vec_cmplx

		SSE2_RADIX16_WRAPPER_DIF(add0,add1,add2,add3
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,pfetch_dist)

	#else			// SSE2: process 2 main-array blocks of [16 vec_dbl = 16 x 2 = 32 doubles] each, total = 32 vec_dbl = 16 vec_cmplx

		SSE2_RADIX16_WRAPPER_DIF(add0,add1
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)

	#endif

/*
!...send the pairs of complex elements which are to be combined and sincos temporaries needed for the squaring to a
!   small subroutine. The j1 = 0 case is again exceptional. [For this reason we don't supply SSE2 code for it - not worth the work].
*/
		/*
		In SSE2 mode, the data are laid out in memory as

			t1 :aj1p0r,aj2p0r		t2 :aj1p0i,aj2p0i
			t3 :aj1p1r,aj2p1r		t4 :aj1p1i,aj2p1i
			t5 :aj1p2r,aj2p2r		t6 :aj1p2i,aj2p2i
			t7 :aj1p3r,aj2p3r		t8 :aj1p3i,aj2p3i
			t9 :aj1p4r,aj2p4r		t10:aj1p4i,aj2p4i
			t11:aj1p5r,aj2p5r		t12:aj1p5i,aj2p5i
			t13:aj1p6r,aj2p6r		t14:aj1p6i,aj2p6i
			t15:aj1p7r,aj2p7r		t16:aj1p7i,aj2p7i
			t17:aj1p8r,aj2p8r		t18:aj1p8i,aj2p8i
			t19:aj1p9r,aj2p9r		t20:aj1p9i,aj2p9i
			t21:aj1pAr,aj2pAr		t22:aj1pAi,aj2pAi
			t23:aj1pBr,aj2pBr		t24:aj1pBi,aj2pBi
			t25:aj1pCr,aj2pCr		t26:aj1pCi,aj2pCi
			t27:aj1pDr,aj2pDr		t28:aj1pDi,aj2pDi
			t29:aj1pEr,aj2pEr		t30:aj1pEi,aj2pEi
			t31:aj1pFr,aj2pFr		t32:aj1pFi,aj2pFi .

		The modified call sequence below takes advantage of that, by processing data which are in 8 XMM registers in-place.
		*/
		t1=(RT-IT)*ISRT2;	t2=(RT+IT)*ISRT2;

		re0=RT*c;	im0=IT*s;	re1=RT*s;	im1=IT*c;
		t3=re0-im0;	t4=re1+im1;
		t5=re1-im1;	t6=re0+im0;

		tmp0->d0 = RT;		tmp1->d0 = IT;
		tmp0->d1 = t6;		tmp1->d1 = t5;

		tmp2->d0 = t1;		tmp3->d0 = t2;
		tmp2->d1 = t4;		tmp3->d1 = t3;

	#ifdef USE_AVX
		RT = c1->d2;		IT = (c1+1)->d2;	// Notice we only use even-indexed ->d* values to set RT,IT in SIMD mode
		t1=(RT-IT)*ISRT2;	t2=(RT+IT)*ISRT2;

		re0=RT*c;	im0=IT*s;	re1=RT*s;	im1=IT*c;
		t3=re0-im0;	t4=re1+im1;
		t5=re1-im1;	t6=re0+im0;

		tmp0->d2 = RT;		tmp1->d2 = IT;
		tmp0->d3 = t6;		tmp1->d3 = t5;

		tmp2->d2 = t1;		tmp3->d2 = t2;
		tmp2->d3 = t4;		tmp3->d3 = t3;
	#endif

	#ifdef USE_AVX512
		RT = c1->d4;		IT = (c1+1)->d4;
		t1=(RT-IT)*ISRT2;	t2=(RT+IT)*ISRT2;

		re0=RT*c;	im0=IT*s;	re1=RT*s;	im1=IT*c;
		t3=re0-im0;	t4=re1+im1;
		t5=re1-im1;	t6=re0+im0;

		tmp0->d4 = RT;		tmp1->d4 = IT;
		tmp0->d5 = t6;		tmp1->d5 = t5;

		tmp2->d4 = t1;		tmp3->d4 = t2;
		tmp2->d5 = t4;		tmp3->d5 = t3;

		RT = c1->d6;		IT = (c1+1)->d6;
		t1=(RT-IT)*ISRT2;	t2=(RT+IT)*ISRT2;

		re0=RT*c;	im0=IT*s;	re1=RT*s;	im1=IT*c;
		t3=re0-im0;	t4=re1+im1;
		t5=re1-im1;	t6=re0+im0;

		tmp0->d6 = RT;		tmp1->d6 = IT;
		tmp0->d7 = t6;		tmp1->d7 = t5;

		tmp2->d6 = t1;		tmp3->d6 = t2;
		tmp2->d7 = t4;		tmp3->d7 = t3;
	#endif

	#ifdef USE_AVX2	// This also includes AVX-512:
		PAIR_SQUARE_4_AVX2(	r1, r9,r23,r31, tmp0,tmp1,
							r5,r13,r19,r27, tmp2,tmp3,forth);
	#else
		PAIR_SQUARE_4_SSE2( r1, r9,r23,r31, tmp0,tmp1,forth);
		PAIR_SQUARE_4_SSE2( r5,r13,r19,r27, tmp2,tmp3,forth);
	#endif

	// Permute the sincos-multiplier data:
	#ifdef USE_ARM_V8_SIMD

		/* Test Code:
			uint128 *itmp0;	itmp0->d1 = 0x000000a3000000a2ull;	itmp0->d0 = 0x000000a1000000a0ull;
			uint128 *itmp1;	itmp1->d1 = 0x000000b3000000b2ull;	itmp1->d0 = 0x000000b1000000b0ull;
		Here results of various shuffle-ops on these inputs:
			uzp1	v2.4s,v0.4s,v1.4s	// v2 = b2,b0,a2,a0
			uzp2	v3.4s,v0.4s,v1.4s	// v3 = b3,b1,a3,a1
			zip1	v4.4s,v0.4s,v1.4s	// v4 = b1,a1,b0,a0
			zip2	v5.4s,v0.4s,v1.4s	// v5 = b3,a3,b2,a2
			trn1	v6.4s,v0.4s,v1.4s	// v6 = b2,a2,b0,a0, i.e. interleave-even-indexed-subelements
			trn2	v7.4s,v0.4s,v1.4s	// v7 = b3,a3,b1,a1, i.e. interleave- odd-indexed-subelements
		Can use TRN1,2 to effect swap(lo64,hi64) on pairs-of-vec-regs; on input v0 = [a1,a0], v1 = [b1,b0]:
			trn1	v4.2d,v0.2d,v1.2d	// v4 = [b0,a0]
			trn2	v5.2d,v0.2d,v1.2d	// v5 = [b1,a1]
			trn1	v0.2d,v5.2d,v4.2d	// v4 = [a0,a1]
			trn2	v1.2d,v5.2d,v4.2d	// v5 = [b0,b1]
		*/
	  #if 0
		// The needed swap-64-bit-halves-of-a-vector-register here was easy in Aarch32, e.g. 'vswp d0,d1'
		// to swap hi & lo halves of q0|v0. That is out the window in Aarch64 due the register model no
		// longer giving us a short-register alias to the upper half of a vreg. So cast about for alternatives -
		// V1 works in vreg-pair fashion, but needs 2 shuffle-ops per vreg:
		__asm__ volatile (\
			"ldr	x0,%[tmp0]	\n\t"\
			"ldr	x1,%[tmp2]	\n\t"\
			"ldp	q0,q1,[x0]	\n\t"/* tmp0,1 occupy adjacent vec_dbl slots... */\
			"ldp	q2,q3,[x1]	\n\t"/* ...as do tmp2,3. (Though those *not* adj. to tmp0,1). */\
			"trn1	v4.2d,v0.2d,v1.2d	\n\t"\
			"trn2	v5.2d,v0.2d,v1.2d	\n\t"\
			"trn1	v0.2d,v5.2d,v4.2d	\n\t"\
			"trn2	v1.2d,v5.2d,v4.2d	\n\t"\
			"trn1	v4.2d,v2.2d,v3.2d	\n\t"\
			"trn2	v5.2d,v2.2d,v3.2d	\n\t"\
			"trn1	v2.2d,v5.2d,v4.2d	\n\t"\
			"trn2	v3.2d,v5.2d,v4.2d	\n\t"\
			"stp	q0,q1,[x0]	\n\t"\
			"stp	q2,q3,[x1]	\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			 ,[tmp2] "m" (tmp2)
			: "cc","memory","x0","x1","v0","v1","v2","v3","v4","v5"	/* Clobbered registers */\
		);
	  #elif 1
		// V2: Pinged ARM's Tom Womack re. possible single-instruction-per-swap alternatives,
		// His reply: "I have asked some of the gurus at work. The answer: make a temp vector
		// by concatenating Vsrc and Vsrc, then extract bytes (8+15)..8 from it into Vdest".
		// Repeater-loop-enclosed timings of the 2 variants shows this EXT-based one needing ~20% fewer cycles:
		__asm__ volatile (\
			"ldr	x0,%[tmp0]	\n\t"\
			"ldr	x1,%[tmp2]	\n\t"\
			"ldp	q0,q1,[x0]	\n\t"/* tmp0,1 occupy adjacent vec_dbl slots... */\
			"ldp	q2,q3,[x1]	\n\t"/* ...as do tmp2,3. (Though those *not* adj. to tmp0,1). */\
			"ext v0.16b,v0.16b,v0.16b,#8	\n\t"\
			"ext v1.16b,v1.16b,v1.16b,#8	\n\t"\
			"ext v2.16b,v2.16b,v2.16b,#8	\n\t"\
			"ext v3.16b,v3.16b,v3.16b,#8	\n\t"\
			"stp	q0,q1,[x0]	\n\t"\
			"stp	q2,q3,[x1]	\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			 ,[tmp2] "m" (tmp2)
			: "cc","memory","x0","x1","v0","v1","v2","v3"	/* Clobbered registers */\
		);
	  #else
		// In C code, just swap the d0/d1 element selectors from our above previous inite of tmp0-3:
		tmp0->d1 = RT;		tmp1->d1 = IT;
		tmp0->d0 = t6;		tmp1->d0 = t5;

		tmp2->d1 = t1;		tmp3->d1 = t2;
		tmp2->d0 = t4;		tmp3->d0 = t3;
	  #endif

	#elif defined(USE_AVX512)
		// AVX-512 version has shufpd immediate = 0x55 = 01010101_2, which is the fourfold analog of the SSE2 imm8 = 1 = 01_2:
		__asm__ volatile (\
			"movq	%[tmp0],%%rax\n\t"\
			"movq	%[tmp1],%%rbx\n\t"\
			"movq	%[tmp2],%%rcx\n\t"\
			"movq	%[tmp3],%%rdx\n\t"\
			"vmovaps	(%%rax),%%zmm0\n\t"\
			"vmovaps	(%%rbx),%%zmm1\n\t"\
			"vmovaps	(%%rcx),%%zmm2\n\t"\
			"vmovaps	(%%rdx),%%zmm3\n\t"\
			"vshufpd	$0x55,%%zmm0,%%zmm0,%%zmm0\n\t"\
			"vshufpd	$0x55,%%zmm1,%%zmm1,%%zmm1\n\t"\
			"vshufpd	$0x55,%%zmm2,%%zmm2,%%zmm2\n\t"\
			"vshufpd	$0x55,%%zmm3,%%zmm3,%%zmm3\n\t"\
			"vmovaps	%%zmm0,(%%rax)\n\t"\
			"vmovaps	%%zmm1,(%%rbx)\n\t"\
			"vmovaps	%%zmm2,(%%rcx)\n\t"\
			"vmovaps	%%zmm3,(%%rdx)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			 ,[tmp1] "m" (tmp1)
			 ,[tmp2] "m" (tmp2)
			 ,[tmp3] "m" (tmp3)
			: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3"		/* Clobbered registers */\
		);

	#elif defined(USE_AVX)
		// AVX version has shufpd immediate = 5 = 0101_2, which is the doubled analog of the SSE2 imm8 = 1 = 01_2:
		__asm__ volatile (\
			"movq	%[tmp0],%%rax\n\t"\
			"movq	%[tmp1],%%rbx\n\t"\
			"movq	%[tmp2],%%rcx\n\t"\
			"movq	%[tmp3],%%rdx\n\t"\
			"vmovaps	(%%rax),%%ymm0\n\t"\
			"vmovaps	(%%rbx),%%ymm1\n\t"\
			"vmovaps	(%%rcx),%%ymm2\n\t"\
			"vmovaps	(%%rdx),%%ymm3\n\t"\
			"vshufpd	$5,%%ymm0,%%ymm0,%%ymm0\n\t"\
			"vshufpd	$5,%%ymm1,%%ymm1,%%ymm1\n\t"\
			"vshufpd	$5,%%ymm2,%%ymm2,%%ymm2\n\t"\
			"vshufpd	$5,%%ymm3,%%ymm3,%%ymm3\n\t"\
			"vmovaps	%%ymm0,(%%rax)\n\t"\
			"vmovaps	%%ymm1,(%%rbx)\n\t"\
			"vmovaps	%%ymm2,(%%rcx)\n\t"\
			"vmovaps	%%ymm3,(%%rdx)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			 ,[tmp1] "m" (tmp1)
			 ,[tmp2] "m" (tmp2)
			 ,[tmp3] "m" (tmp3)
			: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3"		/* Clobbered registers */\
		);

	#elif(OS_BITS == 64)	// 64-bit SSE2 build

		__asm__ volatile (\
			"movq	%[tmp0],%%rax\n\t"\
			"movq	%[tmp1],%%rbx\n\t"\
			"movq	%[tmp2],%%rcx\n\t"\
			"movq	%[tmp3],%%rdx\n\t"\
			"movaps	(%%rax),%%xmm0\n\t"\
			"movaps	(%%rbx),%%xmm1\n\t"\
			"movaps	(%%rcx),%%xmm2\n\t"\
			"movaps	(%%rdx),%%xmm3\n\t"\
			"shufpd	$1	,%%xmm0	,%%xmm0\n\t"\
			"shufpd	$1	,%%xmm1	,%%xmm1\n\t"\
			"shufpd	$1	,%%xmm2	,%%xmm2\n\t"\
			"shufpd	$1	,%%xmm3	,%%xmm3\n\t"\
			"movaps	%%xmm0,(%%rax)\n\t"\
			"movaps	%%xmm1,(%%rbx)\n\t"\
			"movaps	%%xmm2,(%%rcx)\n\t"\
			"movaps	%%xmm3,(%%rdx)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			 ,[tmp1] "m" (tmp1)
			 ,[tmp2] "m" (tmp2)
			 ,[tmp3] "m" (tmp3)
			: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3"		/* Clobbered registers */\
		);

	#else					// 32-bit SSE2 build

		__asm__ volatile (\
			"movl	%[tmp0],%%eax\n\t"\
			"movl	%[tmp1],%%esi\n\t"\
			"movl	%[tmp2],%%ecx\n\t"\
			"movl	%[tmp3],%%edx\n\t"\
			"movaps	(%%eax),%%xmm0\n\t"\
			"movaps	(%%esi),%%xmm1\n\t"\
			"movaps	(%%ecx),%%xmm2\n\t"\
			"movaps	(%%edx),%%xmm3\n\t"\
			"shufpd	$1	,%%xmm0	,%%xmm0\n\t"\
			"shufpd	$1	,%%xmm1	,%%xmm1\n\t"\
			"shufpd	$1	,%%xmm2	,%%xmm2\n\t"\
			"shufpd	$1	,%%xmm3	,%%xmm3\n\t"\
			"movaps	%%xmm0,(%%eax)\n\t"\
			"movaps	%%xmm1,(%%esi)\n\t"\
			"movaps	%%xmm2,(%%ecx)\n\t"\
			"movaps	%%xmm3,(%%edx)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			 ,[tmp1] "m" (tmp1)
			 ,[tmp2] "m" (tmp2)
			 ,[tmp3] "m" (tmp3)
			: "cc","memory","eax","esi","ecx","edx","xmm0","xmm1","xmm2","xmm3"		// Clobbered registers
		);

	#endif

	#ifdef USE_AVX2
		PAIR_SQUARE_4_AVX2(	r3,r11,r21,r29, tmp3,tmp2,
							r7,r15,r17,r25, tmp1,tmp0,forth);
	#else
		PAIR_SQUARE_4_SSE2( r3,r11,r21,r29, tmp3,tmp2,forth);
		PAIR_SQUARE_4_SSE2( r7,r15,r17,r25, tmp1,tmp0,forth);
	#endif
	/******** NB: The cost of each PAIR_SQUARE_4 call costs ~1/4 the cost of a single radix-16 DIF or DIT pass,
			  so this entire sequence costs ~= 1 radix-16 pass, thus the entire function costs ~3 radix-16 passes.
	*********/

	/*********************************************************************/
	/*...And do an inverse DIT radix-16 pass on the squared-data blocks: */
	/*********************************************************************/
	#ifdef USE_AVX512

		SSE2_RADIX16_WRAPPER_DIT(add0,add1,add2,add3,add4,add5,add6,add7
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,pfetch_dist)

	#elif defined(USE_AVX)

		SSE2_RADIX16_WRAPPER_DIT(add0,add1,add2,add3
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,pfetch_dist)

	#else	// SSE2:

		SSE2_RADIX16_WRAPPER_DIT(add0,add1
								,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)
	  #if 0	// Debug:
		if(j1 < 1024) {
			int idbg;
			printf("j1 = %u: SSE2 DIT outputs:\n",j1);
			tmp = (vec_dbl*)add0;
			for(idbg = 0; idbg < 16; idbg++) { printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1); tmp+=2; }
			printf("j2 = %u: SSE2 DIT outputs:\n",j2);
			tmp = (vec_dbl*)add1;
			for(idbg = 0; idbg < 16; idbg++) { printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1); tmp+=2; }
			exit(0);
		}
	  #endif
	#endif // AVX or SSE2?

	}	// endif(j1 == 0)

#endif	/* USE_SSE2 */

/*...Update the data (j1 and j2) array indices. */
loop:
	#ifdef USE_AVX
	  if(j1 <= 160) {
	#else
	  if(j1 <= 64) {
	#endif
		// Use scalar code (with index offsets properly fiddled) for j1 == 0 case in SIMD mode
		j1 = j1+32;
		j2 = j2-32;
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
!   second execution of the above loop. The exception is the first loop execution, where j1 needs to be doubled (32 x 2).
*/

update_blocklen:

	j1 = j1+(blocklen << 1);
//	fprintf(stderr,"(j2_start == %d\n",j2_start);
	if(j2_start == n-32)//	if(j2_start == n-stride)
	{
		j1 = 0;
//	fprintf(stderr,"(j2_start == n-32) check hit: returning\n");	exit(0);
		return;
	}

	/*
	!...Since the foregoing loop only gets executed half as many times as in the simple version, to properly position ourselves in the data
	!   array for the start of the next block, need to bump up j1 by as much as would occur in a second execution of the above loop.
	!
	!   Examples:
	!
	!* N=768=3*2^8                                  complex data index ranges
	!*      =1100000000_2.          blocklength:    evens:  j1      odds:   j2        j2_next = j2_start + next blocklength
	!* k=N-1=1011111111, lobit(k>>1) = 1     2              2               3         7
	!*      k=101111111, lobit(k>>1) = 1,    4              4-6             7-5       15
	!*      k= 10111111, lobit(k>>1) = 1,    8              8-14            15-9      31
	!*      k=  1011111, lobit(k>>1) = 1,   16              16-30           31-17     63
	!*      k=   101111, lobit(k>>1) = 1,   32              32-62           63-33     127
	!*      k=    10111, lobit(k>>1) = 1,   64              64-126          127-65    255
	!*      k=     1011, lobit(k>>1) = 1,  128              128-254         255-129   767
	!*      k=      101, lobit(k>>1) = 0,  512 = 128*(k-1)  256-766         767-257   j2_start = N-1; STOP    <---NOTE: sum of blocklengths up to this point + 2 = N.
	!*                                (The added 2 is for elements 0 and 1, processed separately.)
	!
	!* N=1280=5*2^8
	!*      =10100000000_2: lobit(k>>1)             Sum:    j1              j2        j2 at start of next block
	!* k=N-1=10011111111      1              2         4    2               3         7
	!*        1001111111      1              4         8    4-6             7-5       15
	!*         100111111      1              8        16    8-14            15-9      31
	!*          10011111      1             16        32    16-30           31-17     63
	!*           1001111      1             32        64    32-62           63-33     127
	!*            100111      1             64       128    64-126          127-65    255
	!*             10011      1            128       256    128-254         255-129   1279
	!*              1001      0 128*(k-1)=1024      1280    256-1278        1279-257  j2_start = N-1; STOP
	!
	!* N=1792=7*2^8
	!*      =11100000000_2: lobit(k>>1)             Sum:    j1              j2        j2 at start of next block
	!* k=N-1=11011111111      1              2         4    2               3         7
	!*        1101111111      1              4         8    4-6             7-5       15
	!*         110111111      1              8        16    8-14            15-9      31
	!*          11011111      1             16        32    16-30           31-17     63
	!*           1101111      1             32        64    32-62           63-33     127
	!*            110111      1             64       128    64-126          127-65    255
	!*             11011      1            128       256    128-254         255-129   1791
	!*              1101      0 128*(k-1)=1536      1792    256-1790        1791-257  j2_start = N-1; STOP
	*/

	/*...Reset half-complex-blocklength for next pass. If K >> 1 has a zero trailing bit, we multiply the blocklength by K >> 1 in preparation for the final block.
	*/

	blocklen_sum = blocklen_sum + blocklen;
	if(i > 0) blocklen = (radix_prim[i-1]-1)*blocklen_sum;	// i = 0 update-value of blobklen not used, but wrap in if(i > 0) to get rid of UMR

	/*...Next j2_start is previous one plus the (real) length of the current block = 4*(half-complex-blocklength) */

	j2_start = j2_start+(blocklen<<2);
	j2=j2_start;			    /* Reset j2 for start of the next block. */

//	fprintf(stderr,"after update_blocklen: j1,j2 = %u, %u\n",j1,j2);

/*printf("newblock: blocklen = %8d blocklen_sum = %8d j2 = %8d\n",blocklen,blocklen_sum,j2);*/
}	 /* End of Main (i) loop */

}

#undef PFETCH_DIST
