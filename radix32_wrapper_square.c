/*******************************************************************************
*                                                                              *
*   (C) 1997-2013 by Ernst W. Mayer.                                           *
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

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix32_wrapper_square_gcc32.h"

		#else

			#include "radix32_wrapper_square_gcc64.h"

		#endif

	#endif

#endif	/* USE_SSE2 */

/***************/

void radix32_wrapper_square(double a[], int arr_scratch[],int n, int radix0, struct complex rt0[], struct complex rt1[], int nradices_prim, int radix_prim[], int ws_i, int ws_j1, int ws_j2, int ws_j2_start, int ws_k, int ws_m, int ws_blocklen, int ws_blocklen_sum, int init_sse2, int thr_id)
{
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
	const int stride = (int)RE_IM_STRIDE << 5;	// main-array loop stride = 64 for sse2, 128 for avx
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
	int i,j1,j2,j2_start,k,m,blocklen,blocklen_sum;
	const double c     = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
				,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784		/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
				,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;	/* exp(3*i*twopi/32)	*/
	double rt,it,re = 0.0, im= 0.0;
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
	,a2p10i,a2p11i,a2p12i,a2p13i,a2p14i,a2p15i,a2p16i,a2p17i,a2p18i,a2p19i,a2p1Ai,a2p1Bi,a2p1Ci,a2p1Di,a2p1Ei,a2p1Fi;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0, *add1;	/* Addresses into array sections */
  #ifdef USE_AVX
	double *add2, *add3;
	const int stridh = (stride>>1);
  #endif
	vec_dbl *c_tmp,*s_tmp;
	vec_dbl *tmp;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
	// In || mode, only above base-pointer (shared by all threads) is static:
	vec_dbl *isrt2,*sqrt2,*one,*two, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3
		, *forth, *tmp0, *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6, *tmp7
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*r00,*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E
		,*r10,*r12,*r14,*r16,*r18,*r1A,*r1C,*r1E
		,*r20,*r22,*r24,*r26,*r28,*r2A,*r2C,*r2E
		,*r30,*r32,*r34,*r36,*r38,*r3A,*r3C,*r3E;
  #else		// Same list of ptrs as above (except for base-ptr), but now make them static:
	static vec_dbl *isrt2,*sqrt2,*one,*two, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3
		, *forth, *tmp0, *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6, *tmp7
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*r00,*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E
		,*r10,*r12,*r14,*r16,*r18,*r1A,*r1C,*r1E
		,*r20,*r22,*r24,*r26,*r28,*r2A,*r2C,*r2E
		,*r30,*r32,*r34,*r36,*r38,*r3A,*r3C,*r3E;
  #endif

#endif

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
		max_threads = init_sse2;
	#ifndef COMPILER_TYPE_GCC
		ASSERT(HERE, NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
	//	printf("max_threads = %d, NTHREADS = %d\n",max_threads, NTHREADS);
		ASSERT(HERE, max_threads >= NTHREADS, "Multithreading requires max_threads >= NTHREADS!");

		nsave = n;

	#ifdef USE_SSE2

		ASSERT(HERE, thr_id == -1, "Init-mode call must be outside of any multithreading!");
		if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)sc_arr);	sc_arr=0x0;
		}
		sc_arr = ALLOC_VEC_DBL(sc_arr, 148*max_threads);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 64 vec_dbl slots of sc_arr for temporaries, next 8 for scratch, next 7 for the nontrivial complex 16th roots,
	next 62 for the doubled sincos twiddles, next 4 for [1.0,2.0,0.25,sqrt2] and at least 3 more to allow for 64-byte alignment of the array.

	*** NOTE ***: Offsets below must match those in radix32_dyadic_square,
				since routines share DFT macros which use many literal byte offsets to reduce argument count.
	*/
		#ifdef MULTITHREAD
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x48;
			cc0	  = sc_ptr + 0x49;
			ss0	  = sc_ptr + 0x4a;
			cc1	  = sc_ptr + 0x4b;
			ss1	  = sc_ptr + 0x4c;
			cc3	  = sc_ptr + 0x4d;
			ss3	  = sc_ptr + 0x4e;
			two   = sc_ptr + 0x8f;	// PAIR_SQUARE_4_AVX2 not used for Fermat-mod (i.e. no need for *forth), but note the SSE2_RADIX32_WRAPPER_DIF/T
			forth = sc_ptr + 0x90;	// asm macros are shared by Mersenne & Fermat, and assume two = isrt2 + 0x47, sqrt2 = two + 2, one = two + 3,
			sqrt2 = sc_ptr + 0x91;	// i.e. need to "leave slot for forth just above two".
			one	  = sc_ptr + 0x92;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
				VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );	VEC_DBL_INIT(forth, 0.25 );
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
				forth += 0x94;
				sqrt2 += 0x94;
			}
		#else
			r00		= sc_ptr + 0x00;	isrt2	= sc_ptr + 0x48;
			r02		= sc_ptr + 0x02;	cc0		= sc_ptr + 0x49;	ss0   = sc_ptr + 0x4a;
			r04		= sc_ptr + 0x04;	cc1		= sc_ptr + 0x4b;	ss1   = sc_ptr + 0x4c;
			r06		= sc_ptr + 0x06;	cc3		= sc_ptr + 0x4d;	ss3   = sc_ptr + 0x4e;
			r08		= sc_ptr + 0x08;	c00		= sc_ptr + 0x4f;
			r0A		= sc_ptr + 0x0a;	c10		= sc_ptr + 0x51;
			r0C		= sc_ptr + 0x0c;	c08		= sc_ptr + 0x53;
			r0E		= sc_ptr + 0x0e;	c18		= sc_ptr + 0x55;
			r10		= sc_ptr + 0x10;	c04		= sc_ptr + 0x57;
			r12		= sc_ptr + 0x12;	c14		= sc_ptr + 0x59;
			r14		= sc_ptr + 0x14;	c0C		= sc_ptr + 0x5b;
			r16		= sc_ptr + 0x16;	c1C		= sc_ptr + 0x5d;
			r18		= sc_ptr + 0x18;	c02		= sc_ptr + 0x5f;
			r1A		= sc_ptr + 0x1a;	c12		= sc_ptr + 0x61;
			r1C		= sc_ptr + 0x1c;	c0A		= sc_ptr + 0x63;
			r1E		= sc_ptr + 0x1e;	c1A		= sc_ptr + 0x65;
			r20		= sc_ptr + 0x20;	c06		= sc_ptr + 0x67;
			r22		= sc_ptr + 0x22;	c16		= sc_ptr + 0x69;
			r24		= sc_ptr + 0x24;	c0E		= sc_ptr + 0x6b;
			r26		= sc_ptr + 0x26;	c1E		= sc_ptr + 0x6d;
			r28		= sc_ptr + 0x28;	c01		= sc_ptr + 0x6f;
			r2A		= sc_ptr + 0x2a;	c11		= sc_ptr + 0x71;
			r2C		= sc_ptr + 0x2c;	c09		= sc_ptr + 0x73;
			r2E		= sc_ptr + 0x2e;	c19		= sc_ptr + 0x75;
			r30		= sc_ptr + 0x30;	c05		= sc_ptr + 0x77;
			r32		= sc_ptr + 0x32;	c15		= sc_ptr + 0x79;
			r34		= sc_ptr + 0x34;	c0D		= sc_ptr + 0x7b;
			r36		= sc_ptr + 0x36;	c1D		= sc_ptr + 0x7d;
			r38		= sc_ptr + 0x38;	c03		= sc_ptr + 0x7f;
			r3A		= sc_ptr + 0x3a;	c13		= sc_ptr + 0x81;
			r3C		= sc_ptr + 0x3c;	c0B		= sc_ptr + 0x83;
			r3E		= sc_ptr + 0x3e;	c1B		= sc_ptr + 0x85;
			tmp0	= sc_ptr + 0x40;	c07		= sc_ptr + 0x87;
			tmp1	= sc_ptr + 0x41;	c17		= sc_ptr + 0x89;
			tmp2	= sc_ptr + 0x42;	c0F		= sc_ptr + 0x8b;
			tmp3	= sc_ptr + 0x43;	c1F		= sc_ptr + 0x8d;
			tmp4	= sc_ptr + 0x44;	two  	= sc_ptr + 0x8f;
			tmp5	= sc_ptr + 0x45;	forth	= sc_ptr + 0x90;
			tmp6	= sc_ptr + 0x46;	sqrt2	= sc_ptr + 0x91;
			tmp7	= sc_ptr + 0x47;	one	 	= sc_ptr + 0x92;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
			VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );	VEC_DBL_INIT(forth, 0.25 );
			VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
			VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
			VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
		#endif

	#endif	// USE_SSE2

		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Initialize a scratch array containing N2/16 indices - again use big
		!   as-yet-unused A-array for this, but making sure the section of A used
		!   for the itmp space and that sent to the bit_reverse_int for scratch space
		!   don't overlap:
		*/
		ASSERT(HERE, N2 == n/2, "N2 bad!");
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
		ASSERT(HERE, index_ptmp != 0,"FATAL: unable to allocate array INDEX!");
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
		  ASSERT(HERE, i != 0,"ERROR 10!");
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
	ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
  #ifdef USE_SSE2
	r00 = __r0 + thr_id*148;isrt2	= r00 + 0x48;
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
!	N=16 (FIRST RADIX = 2):			N=24 (FIRST RADIX = 3):			N=20 (FIRST RADIX = 5):			N=28 (FIRST RADIX = 7):
!
!	N - 1 = 15 = 1111_2:			N - 1 = 23 = 10111_2:			N - 1 = 19 = 10011_2:			N - 1 = 27 = 11011_2:
!	j1 =	j2 =		Current bits	j1 =	j2 =		Current bits	j1 =	j2 =		Current bits	j1 =	j2 =		Current bits
!	index(j)index(N-j)	of (N - 1)_2:	index(j)index(N-j)	of (N - 1)_2:	index(j)index(N-j)	of (N - 1)_2:	index(j)index(N-j)	of (N - 1)_2:
!	------	------		-------------	------	------		-------------	------	------		-------------	------	------		-------------
!	0,1 (j=0,N/2) done separately	1	0,1 (j=0,N/2) done separately	1	0,1 (j=0,N/2) done separately	1	0,1 (j=0,N/2) done separately	1
!	------	------				------	------				------	------				------	------
!	2	3	length-2 block	1	2	3	length-2 block	1	2	3	length-2 block	1	2	3	length-2 block	1
!	------	------				------	------				------	------				------	------
!	4	7	length-4 block	1	4	7	length-4 block	1	4	19				4	27
!	6	5				6	5				6	17				6	25
!	------	------				------	------				8	15				8	23
!	8	15				8	23				10	13	length-16 block	100	10	21
!	10	13	length-8 block	1	10	21				12	11				12	19
!	12	11				12	19				14	9				14	17	length-24 block	110
!	14	9				14	17	length-16 block	10	16	7				16	15
!	------	------				16	15				18	5				18	13
!						18	13				------	------				20	11
!						20	11									22	9
!						22	9									24	7
!						------	------									26	5
!																------	------
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

//	fprintf(stderr,"stride = %d\n",stride);
//	fprintf(stderr,"On entry: j1,j2 = %u, %u, nradices_prim = %u, blocklen = %u\n",j1,j2,nradices_prim,blocklen);

	/* If j1 == 0 we need to init the loop counters; otherwise, just jump
	   right in and pick up where we left off on the previous pair of blocks:
	*/
	if(j1 > 0) {
	//	fprintf(stderr,"Jumping into loop!\n");
		goto jump_in;
	}

/*
!...All but the first two radix-16 blocks are done on Mr. Henry Ford's assembly line. From there onward the blocklength
!   will always be an integer multiple of 32, i.e. we can process each block using pairs of nonoverlapping blocks of 32
!   complex data each, which is compatible to fusion with radix-32 pass routines.
*/

for(i = nradices_prim-6; i >= 0; i-- )	/* Main loop: lower bound = nradices_prim-radix_now. */
{						/* Remember, radices get processed in reverse order here as in forward FFT. */

#ifdef USE_AVX
	for(m = 0; m < (blocklen-1)>>1; m += 32) /* In AVX mode, process two 32-element sets per loop execution, thus only execute the loop half as many times as for scalar/SSE2 case. */
#else
	for(m = 0; m < (blocklen-1)>>1; m += 16) /* Since we now process TWO 16-element sets per loop execution, only execute the loop half as many times as before. */
#endif
	{
		// This tells us when we've reached the end of the current data block:
		// Apr 2014: Must store intermediate product j1*radix0 in a 64-bit int to prevent overflow!
		if(j1 && ((uint64)j1*radix0)%n == 0)
		{
		//	fprintf(stderr,"(j1 && j1*radix0 == 0 (mod n)) check hit: returning\n");
			return;
		}

jump_in:	/* Entry point for all blocks but the first. */

	  j1pad = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* floating padded-array 1st element index is here */
	  j2pad = j2 + ( (j2 >> DAT_BITS) << PAD_BITS );	/* floating padded-array 2nd element index is here */

	/*************************************************************/
	/*                  1st set of sincos:                       */
	/*************************************************************/
		iroot = index[k++];
		l = iroot;

	#ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-31] are offset w.r.to the thread-local ptr pair as

		             c 00 10 08 18 04 14 0C 1C 02 12 0A 1A 06 16 0E 1E 01 11 09 19 05 15 0D 1D 03 13 0B 1B 07 17 0F 1F
		(cc0,ss0) + 0x[06,08,0a,0c,0e,10,12,14,16,18,1a,1c,1e,20,22,24,26,28,2a,2c,2e,30,32,34,36,38,3a,3c,3e,40,42,44]

		Here, due to the need to compute a new set of roots
		for each set of inputs, we use a streamlined sequence which computes only the [0,1,2,3,7,E,15,1C]th roots with
		maximal accuracy (i.e. using 2-table-multiply), then generates the remaining ones from those. Thus the needed
		pointer offsets below are (cc0,ss0) + 0x[06,26,16,36,3e,22,30,14]:
		*/
	#endif

	#ifdef USE_SSE2
		c_tmp = cc0+0x06; s_tmp = c_tmp+1;	// c0,s0
		VEC_DBL_INIT(c_tmp, 1.0);
		VEC_DBL_INIT(s_tmp, 0.0);
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;	// [c,s]1
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		cA01 =rt;	sA01 =it;
	#endif

	re= rt;	im= it;	/* Save for the wrapper Step... */

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;	// [c,s]2
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		cA02 =rt;	sA02 =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += 4*iroot;			/* 7*iroot	*/
	    iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;	// [c,s]3
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		cA03 =rt;	sA03 =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;	// [c,s]7
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		cA07 =rt;	sA07 =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;	// [c,s]E
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		cA0E =rt;	sA0E =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;	// [c,s]15
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		cA15 =rt;	sA15 =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;	// [c,s]1C
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		cA1C =rt;	sA1C =it;

		/*(c,s)4-10 (decimal indices), 4-0A (hexadecimal indices):	*/
	    t00=cA01*cA07; t01=cA01*sA07; t02=sA01*cA07; t03=sA01*sA07;
	    cA06=t00+t03; sA06=t01-t02; cA08=t00-t03; sA08=t01+t02;

	    t00=cA02*cA07; t01=cA02*sA07; t02=sA02*cA07; t03=sA02*sA07;
	    cA05=t00+t03; sA05=t01-t02; cA09=t00-t03; sA09=t01+t02;

	    t00=cA03*cA07; t01=cA03*sA07; t02=sA03*cA07; t03=sA03*sA07;
	    cA04=t00+t03; sA04=t01-t02; cA0A=t00-t03; sA0A=t01+t02;

		/*(c,s)11-17 (decimal indices), 0B-11 (hexadecimal indices):;	*/
	    t00=cA01*cA0E; t01=cA01*sA0E; t02=sA01*cA0E; t03=sA01*sA0E;
	    cA0D=t00+t03; sA0D=t01-t02; cA0F=t00-t03; sA0F=t01+t02;

	    t00=cA02*cA0E; t01=cA02*sA0E; t02=sA02*cA0E; t03=sA02*sA0E;
	    cA0C=t00+t03; sA0C=t01-t02; cA10=t00-t03; sA10=t01+t02;

	    t00=cA03*cA0E; t01=cA03*sA0E; t02=sA03*cA0E; t03=sA03*sA0E;
	    cA0B=t00+t03; sA0B=t01-t02; cA11=t00-t03; sA11=t01+t02;

		/*(c,s)18-24 (decimal indices), 12-18 (hexadecimal indices):	*/
	    t00=cA01*cA15; t01=cA01*sA15; t02=sA01*cA15; t03=sA01*sA15;
	    cA14=t00+t03; sA14=t01-t02; cA16=t00-t03; sA16=t01+t02;

	    t00=cA02*cA15; t01=cA02*sA15; t02=sA02*cA15; t03=sA02*sA15;
	    cA13=t00+t03; sA13=t01-t02; cA17=t00-t03; sA17=t01+t02;

	    t00=cA03*cA15; t01=cA03*sA15; t02=sA03*cA15; t03=sA03*sA15;
	    cA12=t00+t03; sA12=t01-t02; cA18=t00-t03; sA18=t01+t02;

		/*(c,s)25-31 (decimal indices), 19-1F (hexadecimal indices):	*/
	    t00=cA01*cA1C; t01=cA01*sA1C; t02=sA01*cA1C; t03=sA01*sA1C;
	    cA1B=t00+t03; sA1B=t01-t02; cA1D=t00-t03; sA1D=t01+t02;

	    t00=cA02*cA1C; t01=cA02*sA1C; t02=sA02*cA1C; t03=sA02*sA1C;
	    cA1A=t00+t03; sA1A=t01-t02; cA1E=t00-t03; sA1E=t01+t02;

	    t00=cA03*cA1C; t01=cA03*sA1C; t02=sA03*cA1C; t03=sA03*sA1C;
	    cA19=t00+t03; sA19=t01-t02; cA1F=t00-t03; sA1F=t01+t02;
	#endif

	/*************************************************************/
	/*                  2nd set of sincos:                       */
	/*************************************************************/
	    iroot = index[k++];
	    l = iroot;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;			/* 2*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;	// [c,s]1
		c_tmp->d1=rt;	s_tmp->d1=it;
	#else
		cB01 =rt;	sB01 =it;
	#endif

	if(j1 == 0){ re= rt;	im= it; }   /* The j1 = 0 case is special... */

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;	// [c,s]2
		c_tmp->d1=rt;	s_tmp->d1=it;
	#else
		cB02 =rt;	sB02 =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += 4*iroot;			/* 7*iroot	*/
	    iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;	// [c,s]3
		c_tmp->d1=rt;	s_tmp->d1=it;
	#else
		cB03 =rt;	sB03 =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;	// [c,s]7
		c_tmp->d1=rt;	s_tmp->d1=it;
	#else
		cB07 =rt;	sB07 =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;	// [c,s]E
		c_tmp->d1=rt;	s_tmp->d1=it;
	#else
		cB0E =rt;	sB0E =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;	// [c,s]15
		c_tmp->d1=rt;	s_tmp->d1=it;
	#else
		cB15 =rt;	sB15 =it;
	#endif

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;	// [c,s]1C
		c_tmp->d1=rt;	s_tmp->d1=it;
	#else
		cB1C =rt;	sB1C =it;

		/*(c,s)4-10 (decimal indices), 4-0A (hexadecimal indices):	*/
	    t00=cB01*cB07; t01=cB01*sB07; t02=sB01*cB07; t03=sB01*sB07;
	    cB06=t00+t03; sB06=t01-t02; cB08=t00-t03; sB08=t01+t02;

	    t00=cB02*cB07; t01=cB02*sB07; t02=sB02*cB07; t03=sB02*sB07;
	    cB05=t00+t03; sB05=t01-t02; cB09=t00-t03; sB09=t01+t02;

	    t00=cB03*cB07; t01=cB03*sB07; t02=sB03*cB07; t03=sB03*sB07;
	    cB04=t00+t03; sB04=t01-t02; cB0A=t00-t03; sB0A=t01+t02;

		/*(c,s)11-17 (decimal indices), 0B-11 (hexadecimal indices):;	*/
	    t00=cB01*cB0E; t01=cB01*sB0E; t02=sB01*cB0E; t03=sB01*sB0E;
	    cB0D=t00+t03; sB0D=t01-t02; cB0F=t00-t03; sB0F=t01+t02;

	    t00=cB02*cB0E; t01=cB02*sB0E; t02=sB02*cB0E; t03=sB02*sB0E;
	    cB0C=t00+t03; sB0C=t01-t02; cB10=t00-t03; sB10=t01+t02;

	    t00=cB03*cB0E; t01=cB03*sB0E; t02=sB03*cB0E; t03=sB03*sB0E;
	    cB0B=t00+t03; sB0B=t01-t02; cB11=t00-t03; sB11=t01+t02;

		/*(c,s)18-24 (decimal indices), 12-18 (hexadecimal indices):	*/
	    t00=cB01*cB15; t01=cB01*sB15; t02=sB01*cB15; t03=sB01*sB15;
	    cB14=t00+t03; sB14=t01-t02; cB16=t00-t03; sB16=t01+t02;

	    t00=cB02*cB15; t01=cB02*sB15; t02=sB02*cB15; t03=sB02*sB15;
	    cB13=t00+t03; sB13=t01-t02; cB17=t00-t03; sB17=t01+t02;

	    t00=cB03*cB15; t01=cB03*sB15; t02=sB03*cB15; t03=sB03*sB15;
	    cB12=t00+t03; sB12=t01-t02; cB18=t00-t03; sB18=t01+t02;

		/*(c,s)25-31 (decimal indices), 19-1F (hexadecimal indices):	*/
	    t00=cB01*cB1C; t01=cB01*sB1C; t02=sB01*cB1C; t03=sB01*sB1C;
	    cB1B=t00+t03; sB1B=t01-t02; cB1D=t00-t03; sB1D=t01+t02;

	    t00=cB02*cB1C; t01=cB02*sB1C; t02=sB02*cB1C; t03=sB02*sB1C;
	    cB1A=t00+t03; sB1A=t01-t02; cB1E=t00-t03; sB1E=t01+t02;

	    t00=cB03*cB1C; t01=cB03*sB1C; t02=sB03*cB1C; t03=sB03*sB1C;
	    cB19=t00+t03; sB19=t01-t02; cB1F=t00-t03; sB1F=t01+t02;
	#endif

	#ifdef USE_AVX

	  if(j1 > 128)	// Sincos data for the 2 starting scalar-mode data blocks get done in SSE2 mode, i.e. only using d0,d1
	  {
	/*************************************************************/
	/*                  3rd set of sincos:                       */
	/*************************************************************/
		iroot = index[k++];
		l = iroot;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;			/* 2*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;	// [c,s]1
		c_tmp->d2=rt;	s_tmp->d2=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;	// [c,s]2
		c_tmp->d2=rt;	s_tmp->d2=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += 4*iroot;			/* 7*iroot	*/
	    iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;	// [c,s]3
		c_tmp->d2=rt;	s_tmp->d2=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;	// [c,s]7
		c_tmp->d2=rt;	s_tmp->d2=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;	// [c,s]E
		c_tmp->d2=rt;	s_tmp->d2=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;	// [c,s]15
		c_tmp->d2=rt;	s_tmp->d2=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;	// [c,s]1C
		c_tmp->d2=rt;	s_tmp->d2=it;

	/*************************************************************/
	/*                  4th set of sincos:                       */
	/*************************************************************/
		iroot = index[k++];
		l = iroot;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;			/* 2*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;	// [c,s]1
		c_tmp->d3=rt;	s_tmp->d3=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;	// [c,s]2
		c_tmp->d3=rt;	s_tmp->d3=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += 4*iroot;			/* 7*iroot	*/
	    iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;	// [c,s]3
		c_tmp->d3=rt;	s_tmp->d3=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;	// [c,s]7
		c_tmp->d3=rt;	s_tmp->d3=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;	// [c,s]E
		c_tmp->d3=rt;	s_tmp->d3=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
	    l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;	// [c,s]15
		c_tmp->d3=rt;	s_tmp->d3=it;

	    k1=(l & NRTM1);
	    k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;	// [c,s]1C
		c_tmp->d3=rt;	s_tmp->d3=it;

	  }	// endif(j1 > 128)

	#endif	// USE_AVX ?

	#ifdef USE_SSE2	// Both SSE2 and AVX share this:

		SSE2_CMUL_EXPO(c01,c07,c06,c08)
		SSE2_CMUL_EXPO(c02,c07,c05,c09)
		SSE2_CMUL_EXPO(c03,c07,c04,c0A)

		SSE2_CMUL_EXPO(c01,c0E,c0D,c0F)
		SSE2_CMUL_EXPO(c02,c0E,c0C,c10)
		SSE2_CMUL_EXPO(c03,c0E,c0B,c11)

		SSE2_CMUL_EXPO(c01,c15,c14,c16)
		SSE2_CMUL_EXPO(c02,c15,c13,c17)
		SSE2_CMUL_EXPO(c03,c15,c12,c18)

		SSE2_CMUL_EXPO(c01,c1C,c1B,c1D)
		SSE2_CMUL_EXPO(c02,c1C,c1A,c1E)
		SSE2_CMUL_EXPO(c03,c1C,c19,c1F)

	  #ifdef USE_AVX	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  					// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:

		// process 4 main-array blocks of 8 vec_dbl = 8 x 4 = 32 doubles each in AVX mode
		add0 = a + j1pad;
		add1 = a + j2pad;
		add2 = add0 + stridh;	// add2 = add0 + [32 doubles, equiv to 8 AVX registers]
		add3 = add1 - stridh;	// Last 2 offsets run in descending order for Mers-mod

	  #else	// SSE2:

		add0 = a + j1pad;
		add1 = a + j2pad;

	  #endif

	if(j1 <= 128)	// Use scalar code (with index offsets properly fiddled) for j1 == 0 case in SIMD mode
	{
		cA01 = c01->d0;	sA01 = (c01+1)->d0;
		cA02 = c02->d0;	sA02 = (c02+1)->d0;
		cA03 = c03->d0;	sA03 = (c03+1)->d0;
		cA04 = c04->d0;	sA04 = (c04+1)->d0;
		cA05 = c05->d0;	sA05 = (c05+1)->d0;
		cA06 = c06->d0;	sA06 = (c06+1)->d0;
		cA07 = c07->d0;	sA07 = (c07+1)->d0;
		cA08 = c08->d0;	sA08 = (c08+1)->d0;
		cA09 = c09->d0;	sA09 = (c09+1)->d0;
		cA0A = c0A->d0;	sA0A = (c0A+1)->d0;
		cA0B = c0B->d0;	sA0B = (c0B+1)->d0;
		cA0C = c0C->d0;	sA0C = (c0C+1)->d0;
		cA0D = c0D->d0;	sA0D = (c0D+1)->d0;
		cA0E = c0E->d0;	sA0E = (c0E+1)->d0;
		cA0F = c0F->d0;	sA0F = (c0F+1)->d0;
		cA10 = c10->d0;	sA10 = (c10+1)->d0;
		cA11 = c11->d0;	sA11 = (c11+1)->d0;
		cA12 = c12->d0;	sA12 = (c12+1)->d0;
		cA13 = c13->d0;	sA13 = (c13+1)->d0;
		cA14 = c14->d0;	sA14 = (c14+1)->d0;
		cA15 = c15->d0;	sA15 = (c15+1)->d0;
		cA16 = c16->d0;	sA16 = (c16+1)->d0;
		cA17 = c17->d0;	sA17 = (c17+1)->d0;
		cA18 = c18->d0;	sA18 = (c18+1)->d0;
		cA19 = c19->d0;	sA19 = (c19+1)->d0;
		cA1A = c1A->d0;	sA1A = (c1A+1)->d0;
		cA1B = c1B->d0;	sA1B = (c1B+1)->d0;
		cA1C = c1C->d0;	sA1C = (c1C+1)->d0;
		cA1D = c1D->d0;	sA1D = (c1D+1)->d0;
		cA1E = c1E->d0;	sA1E = (c1E+1)->d0;
		cA1F = c1F->d0;	sA1F = (c1F+1)->d0;
	re= cA01;	im= sA01;	/* Save for the wrapper Step... */

		cB01 = c01->d1;	sB01 = (c01+1)->d1;
		cB02 = c02->d1;	sB02 = (c02+1)->d1;
		cB03 = c03->d1;	sB03 = (c03+1)->d1;
		cB04 = c04->d1;	sB04 = (c04+1)->d1;
		cB05 = c05->d1;	sB05 = (c05+1)->d1;
		cB06 = c06->d1;	sB06 = (c06+1)->d1;
		cB07 = c07->d1;	sB07 = (c07+1)->d1;
		cB08 = c08->d1;	sB08 = (c08+1)->d1;
		cB09 = c09->d1;	sB09 = (c09+1)->d1;
		cB0A = c0A->d1;	sB0A = (c0A+1)->d1;
		cB0B = c0B->d1;	sB0B = (c0B+1)->d1;
		cB0C = c0C->d1;	sB0C = (c0C+1)->d1;
		cB0D = c0D->d1;	sB0D = (c0D+1)->d1;
		cB0E = c0E->d1;	sB0E = (c0E+1)->d1;
		cB0F = c0F->d1;	sB0F = (c0F+1)->d1;
		cB10 = c10->d1;	sB10 = (c10+1)->d1;
		cB11 = c11->d1;	sB11 = (c11+1)->d1;
		cB12 = c12->d1;	sB12 = (c12+1)->d1;
		cB13 = c13->d1;	sB13 = (c13+1)->d1;
		cB14 = c14->d1;	sB14 = (c14+1)->d1;
		cB15 = c15->d1;	sB15 = (c15+1)->d1;
		cB16 = c16->d1;	sB16 = (c16+1)->d1;
		cB17 = c17->d1;	sB17 = (c17+1)->d1;
		cB18 = c18->d1;	sB18 = (c18+1)->d1;
		cB19 = c19->d1;	sB19 = (c19+1)->d1;
		cB1A = c1A->d1;	sB1A = (c1A+1)->d1;
		cB1B = c1B->d1;	sB1B = (c1B+1)->d1;
		cB1C = c1C->d1;	sB1C = (c1C+1)->d1;
		cB1D = c1D->d1;	sB1D = (c1D+1)->d1;
		cB1E = c1E->d1;	sB1E = (c1E+1)->d1;
		cB1F = c1F->d1;	sB1F = (c1F+1)->d1;
	if(j1 == 0){ re= cB01;	im= sB01; }   /* The j1 = 0 case is special... */

	#endif

	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/

		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;

/*...Block 1:	*/
	    t00=a[rdum   ];							t01=a[idum   ];
	    rt =a[rdum+32]*cA10 - a[idum+32]*sA10;	it =a[idum+32]*cA10 + a[rdum+32]*sA10;
	    t02=t00-rt;		t03=t01-it;
	    t00=t00+rt;		t01=t01+it;

	    t04=a[rdum+16]*cA08 - a[idum+16]*sA08;	t05=a[idum+16]*cA08 + a[rdum+16]*sA08;
	    rt =a[rdum+48]*cA18 - a[idum+48]*sA18;	it =a[idum+48]*cA18 + a[rdum+48]*sA18;
	    t06=t04-rt;		t07=t05-it;
	    t04=t04+rt;		t05=t05+it;

	    rt =t04;		it =t05;
	    t04=t00-rt;		t05=t01-it;
	    t00=t00+rt;		t01=t01+it;

	    rt =t06;		it =t07;
	    t06=t02+it;		t07=t03-rt;
	    t02=t02-it;		t03=t03+rt;

	    t08=a[rdum+8 ]*cA04 - a[idum+8 ]*sA04;	t09=a[idum+8 ]*cA04 + a[rdum+8 ]*sA04;
	    rt =a[rdum+40]*cA14 - a[idum+40]*sA14;	it =a[idum+40]*cA14 + a[rdum+40]*sA14;
	    t0A=t08-rt;		t0B=t09-it;
	    t08=t08+rt;		t09=t09+it;

	    t0C=a[rdum+24]*cA0C - a[idum+24]*sA0C;	t0D=a[idum+24]*cA0C + a[rdum+24]*sA0C;
	    rt =a[rdum+56]*cA1C - a[idum+56]*sA1C;	it =a[idum+56]*cA1C + a[rdum+56]*sA1C;
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
	#ifdef USE_AVX
		rdum -= 2;
		idum -= 2;
	#endif
	    t10=a[rdum+4 ]*cA02 - a[idum+4 ]*sA02;	t11=a[idum+4 ]*cA02 + a[rdum+4 ]*sA02;
	    rt =a[rdum+36]*cA12 - a[idum+36]*sA12;	it =a[idum+36]*cA12 + a[rdum+36]*sA12;
	    t12=t10-rt;		t13=t11-it;
	    t10=t10+rt;		t11=t11+it;

	    t14=a[rdum+20]*cA0A - a[idum+20]*sA0A;	t15=a[idum+20]*cA0A + a[rdum+20]*sA0A;
	    rt =a[rdum+52]*cA1A - a[idum+52]*sA1A;	it =a[idum+52]*cA1A + a[rdum+52]*sA1A;
	    t16=t14-rt;		t17=t15-it;
	    t14=t14+rt;		t15=t15+it;

	    rt =t14;		it =t15;
	    t14=t10-rt;		t15=t11-it;
	    t10=t10+rt;		t11=t11+it;

	    rt =t16;		it =t17;
	    t16=t12+it;		t17=t13-rt;
	    t12=t12-it;		t13=t13+rt;

	    t18=a[rdum+12]*cA06 - a[idum+12]*sA06;	t19=a[idum+12]*cA06 + a[rdum+12]*sA06;
	    rt =a[rdum+44]*cA16 - a[idum+44]*sA16;	it =a[idum+44]*cA16 + a[rdum+44]*sA16;
	    t1A=t18-rt;		t1B=t19-it;
	    t18=t18+rt;		t19=t19+it;

	    t1C=a[rdum+28]*cA0E - a[idum+28]*sA0E;	t1D=a[idum+28]*cA0E + a[rdum+28]*sA0E;
	    rt =a[rdum+60]*cA1E - a[idum+60]*sA1E;	it =a[idum+60]*cA1E + a[rdum+60]*sA1E;
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
	#ifdef USE_AVX
		++rdum;
		++idum;
	#elif defined(USE_SSE2)
		--rdum;
		--idum;
	#endif
	    t20=a[rdum+2 ]*cA01 - a[idum+2 ]*sA01;	t21=a[idum+2 ]*cA01 + a[rdum+2 ]*sA01;
	    rt =a[rdum+34]*cA11 - a[idum+34]*sA11;	it =a[idum+34]*cA11 + a[rdum+34]*sA11;
	    t22=t20-rt;		t23=t21-it;
	    t20=t20+rt;		t21=t21+it;

	    t24=a[rdum+18]*cA09 - a[idum+18]*sA09;	t25=a[idum+18]*cA09 + a[rdum+18]*sA09;
	    rt =a[rdum+50]*cA19 - a[idum+50]*sA19;	it =a[idum+50]*cA19 + a[rdum+50]*sA19;
	    t26=t24-rt;		t27=t25-it;
	    t24=t24+rt;		t25=t25+it;

	    rt =t24;		it =t25;
	    t24=t20-rt;		t25=t21-it;
	    t20=t20+rt;		t21=t21+it;

	    rt =t26;		it =t27;
	    t26=t22+it;		t27=t23-rt;
	    t22=t22-it;		t23=t23+rt;

	    t28=a[rdum+10]*cA05 - a[idum+10]*sA05;	t29=a[idum+10]*cA05 + a[rdum+10]*sA05;
	    rt =a[rdum+42]*cA15 - a[idum+42]*sA15;	it =a[idum+42]*cA15 + a[rdum+42]*sA15;
	    t2A=t28-rt;		t2B=t29-it;
	    t28=t28+rt;		t29=t29+it;

	    t2C=a[rdum+26]*cA0D - a[idum+26]*sA0D;	t2D=a[idum+26]*cA0D + a[rdum+26]*sA0D;
	    rt =a[rdum+58]*cA1D - a[idum+58]*sA1D;	it =a[idum+58]*cA1D + a[rdum+58]*sA1D;
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
	#ifdef USE_AVX
		rdum -= 2;
		idum -= 2;
	#endif
	    t30=a[rdum+6 ]*cA03 - a[idum+6 ]*sA03;	t31=a[idum+6 ]*cA03 + a[rdum+6 ]*sA03;
	    rt =a[rdum+38]*cA13 - a[idum+38]*sA13;	it =a[idum+38]*cA13 + a[rdum+38]*sA13;
	    t32=t30-rt;		t33=t31-it;
	    t30=t30+rt;		t31=t31+it;

	    t34=a[rdum+22]*cA0B - a[idum+22]*sA0B;	t35=a[idum+22]*cA0B + a[rdum+22]*sA0B;
	    rt =a[rdum+54]*cA1B - a[idum+54]*sA1B;	it =a[idum+54]*cA1B + a[rdum+54]*sA1B;
	    t36=t34-rt;		t37=t35-it;
	    t34=t34+rt;		t35=t35+it;

	    rt =t34;		it =t35;
	    t34=t30-rt;		t35=t31-it;
	    t30=t30+rt;		t31=t31+it;

	    rt =t36;		it =t37;
	    t36=t32+it;		t37=t33-rt;
	    t32=t32-it;		t33=t33+rt;

	    t38=a[rdum+14]*cA07 - a[idum+14]*sA07;	t39=a[idum+14]*cA07 + a[rdum+14]*sA07;
	    rt =a[rdum+46]*cA17 - a[idum+46]*sA17;	it =a[idum+46]*cA17 + a[rdum+46]*sA17;
	    t3A=t38-rt;		t3B=t39-it;
	    t38=t38+rt;		t39=t39+it;

	    t3C=a[rdum+30]*cA0F - a[idum+30]*sA0F;	t3D=a[idum+30]*cA0F + a[rdum+30]*sA0F;
	    rt =a[rdum+62]*cA1F - a[idum+62]*sA1F;	it =a[idum+62]*cA1F + a[rdum+62]*sA1F;
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
	#ifdef USE_AVX
		rdum -= 2;
		idum -= 2;
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
	#ifdef USE_AVX
		++rdum;
		++idum;
	#elif defined(USE_SSE2)
		--rdum;
		--idum;
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
	#ifdef USE_AVX
		rdum -= 2;
		idum -= 2;
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
	#ifdef USE_AVX
		rdum -= 2;
		idum -= 2;
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
	#ifdef USE_AVX
		++rdum;
		++idum;
	#elif defined(USE_SSE2)
		--rdum;
		--idum;
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
	#ifdef USE_AVX
		rdum -= 2;
		idum -= 2;
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
	#ifdef USE_AVX
		rdum -= 2;
		idum -= 2;
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
	#ifdef USE_AVX
		++rdum;
		++idum;
	#elif defined(USE_SSE2)
		--rdum;
		--idum;
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
	#ifdef USE_AVX
		rdum -= 2;
		idum -= 2;
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

	/* [cf. the radix-16 analog if this file for notes on the AVX data-alyout]
	In sse2 mode, the input data are arranged in memory like so, where we view things in 16-byte chunks:

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

	#if defined(COMPILER_TYPE_MSVC)

	/************************************************************************/
	/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */
	/************************************************************************/

	  #if(0)	/* The "else" block below is the fully expaanded and argcount-reduced version of this macro sequence, used to generate the GCC inline ASM */

	/*...Block 0:	*/
		__asm	mov	eax,add0
		__asm	mov	ebx,add1
		SSE2_RADIX4_DIF_4WRAPPER         (c00,c08,c10,c18,r00)
		SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00)

	/*...Block 1:	*/
		__asm	add	eax,0x20	/* add0 += 4 */
		__asm	add	ebx,0x20	/* add1 += 4 */
		SSE2_RADIX4_DIF_4WRAPPER         (c02,c0A,c12,c1A,r10)
		SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10)
		__asm	sub	eax,0x20	/* add0 -= 4 */
		__asm	sub	ebx,0x20	/* add1 -= 4 */

	/***************************************************************************************************************************
	Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks
	[operating on the odd-indexed elements from the unpck*pd commands which were stored to temporaries can use a common macro:
	***************************************************************************************************************************/
	/*...Block 2:	*/
		SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01)
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)

	/*...Block 3:	*/
		SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30)

	  #else

	/*...Block 0:	*/
		__asm	mov	eax,add0
		__asm	mov	ebx,add1
	/*****	SSE2_RADIX4_DIF_4WRAPPER         (c00,c08,c10,c18,r00) *****/
		__asm	mov	ecx,r00
		__asm	mov	edx,c00
		/* Do the p00,p10 combo: */
		/* For interleaved [j1,j2] version, replace [...] by the following: */
		/* Real parts: */
		__asm	movaps		xmm6,[eax     ]	/* a[j1+p0 ], this is the scratch xmm register */
		__asm	movaps		xmm0,[eax     ]	/* a[j1+p0 ], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx     ]	/* a[j2+p0 ] gets read twice */
		__asm	unpcklpd	xmm0,[ebx     ]	/* a[jt+p0 ] */
		__asm	movaps	[ecx+0x200],xmm6	/* Store hi real in t00+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x10]
		__asm	movaps		xmm1,[eax+0x10]
		__asm	unpckhpd	xmm7,[ebx+0x10]
		__asm	unpcklpd	xmm1,[ebx+0x10]	/* a[jp+p0 ] */
		__asm	movaps	[ecx+0x210],xmm7	/* Store hi imag in t01+32 */

		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy a[jt+p0] */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy a[jp+p0] */

		/* From here on, things are identical to the code in radix32_dif_pass: */
		__asm	mulpd	xmm0,[edx     ]	/* a[jt+p0]*c0 */
		__asm	mulpd	xmm1,[edx     ]	/* a[jp+p0]*c0 */
		__asm	mulpd	xmm2,[edx+0x10]	/* a[jt+p0]*s0 */
		__asm	mulpd	xmm3,[edx+0x10]	/* a[jp+p0]*s0 */
		__asm	addpd	xmm1,xmm2	/* xmm1 <- t01*/
		__asm	subpd	xmm0,xmm3	/* xmm0 <- t00*/
		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t00*/
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t01*/

		__asm	add	edx,0x020	/* c10 */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x100]	/* a[j1+p10], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x100]	/* a[j1+p10], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x100]	/* a[j2+p10] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x100]	/* a[jt+p10] */
		__asm	movaps	[ecx+0x220],xmm6	/* Store hi real in t11+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x110]
		__asm	movaps		xmm5,[eax+0x110]
		__asm	unpckhpd	xmm7,[ebx+0x110]
		__asm	unpcklpd	xmm5,[ebx+0x110]	/* a[jp+p10] */
		__asm	movaps	[ecx+0x230],xmm7	/* Store hi imag in t12+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p10] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p10] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p10]*c10 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p10]*c10 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p10]*s10 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p10]*s10 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */

		__asm	addpd	xmm0,xmm4	/* ~t1 <- t1 +rt */
		__asm	addpd	xmm1,xmm5	/* ~t2 <- t2 +it */
		__asm	subpd	xmm2,xmm4	/* ~t3 <- t1 -rt */
		__asm	subpd	xmm3,xmm5	/* ~t4 <- t2 -it	xmm4,5 free */

		/* Do the p08,p18 [hexadecimal here] combo - do p18 first so register assignments come out in same relative order as for p0,8 */
		__asm	add	edx,0x040	/* c18 */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x180]	/* a[j1+p18], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x180]	/* a[j1+p18], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x180]	/* a[j2+p18] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x180]	/* a[jt+p18] */
		__asm	movaps	[ecx+0x260],xmm6	/* Store hi real in t15+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x190]
		__asm	movaps		xmm5,[eax+0x190]
		__asm	unpckhpd	xmm7,[ebx+0x190]
		__asm	unpcklpd	xmm5,[ebx+0x190]	/* a[jp+p18] */
		__asm	movaps	[ecx+0x270],xmm7	/* Store hi imag in t16+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p18] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p18] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p18]*c18 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p18]*c18 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p18]*s18 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p18]*s18 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */
		__asm	movaps	[ecx+0x010],xmm5	/* Store it */
		__asm	movaps	[ecx      ],xmm4	/* Store rt */

		__asm	sub	edx,0x020	/* c08 */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x080]	/* a[j1+p08], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x080]	/* a[j1+p08], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x080]	/* a[j2+p08] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x080]	/* a[jt+p08] */
		__asm	movaps	[ecx+0x240],xmm6	/* Store hi real in t13+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x090]
		__asm	movaps		xmm5,[eax+0x090]
		__asm	unpckhpd	xmm7,[ebx+0x090]
		__asm	unpcklpd	xmm5,[ebx+0x090]	/* a[jp+p08] */
		__asm	movaps	[ecx+0x250],xmm7	/* Store hi imag in t14+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p08] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p08] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p08]*c08 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p08]*c08 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p08]*s08 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p08]*s08 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */

		__asm	subpd	xmm4,[ecx      ]	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm5,[ecx+0x010]	/* ~t8 <- t6 -it */
		__asm	addpd	xmm6,[ecx      ]	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm7,[ecx+0x010]	/* ~t6 <- t6 +it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;
		~t6 =t2 -t6;		~t2 =t2 +t6;
		*/
		__asm	subpd	xmm0,xmm6	/*~t5 */
		__asm	subpd	xmm1,xmm7	/*~t6 */
		__asm	movaps	[ecx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */
		__asm	movaps	[ecx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */
		__asm	addpd	xmm6,xmm6	/* 2*t5 */
		__asm	addpd	xmm7,xmm7	/* 2*t6 */
		__asm	addpd	xmm6,xmm0	/*~t1 */
		__asm	addpd	xmm7,xmm1	/*~t2 */
		__asm	movaps	[ecx      ],xmm6	/* a[jt    ] <- ~t1 */
		__asm	movaps	[ecx+0x010],xmm7	/* a[jp    ] <- ~t2 */

		/*
		~t7 =t3 +t8;		~t3 =t3 -t8;
		~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	subpd	xmm2,xmm5	/*~t3 */
		__asm	subpd	xmm3,xmm4	/*~t8 */
		__asm	movaps	[ecx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[ecx+0x070],xmm3	/* a[jp+p12] <- ~t8 */
		__asm	addpd	xmm5,xmm5	/* 2*t8 */
		__asm	addpd	xmm4,xmm4	/* 2*t7 */
		__asm	addpd	xmm5,xmm2	/*~t7 */
		__asm	addpd	xmm4,xmm3	/*~t4 */
		__asm	movaps	[ecx+0x060],xmm5	/* a[jt+p12] <- ~t7 */
		__asm	movaps	[ecx+0x030],xmm4	/* a[jp+p4 ] <- ~t4 */

	/*****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/
		__asm	add	ecx,0x080	/* r08 */
		__asm	add	edx,0x040	/* c04 */

		/* Do the p00,p10 combo: */
		/* For interleaved [j1,j2] version, replace [...] by the following:
		*/
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x40]	/* a[j1+p4 ], this is the scratch xmm register */
		__asm	movaps		xmm0,[eax+0x40]	/* a[j1+p4 ], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x40]	/* a[j2+p4 ] gets read twice */
		__asm	unpcklpd	xmm0,[ebx+0x40]	/* a[jt+p4 ] */
		__asm	movaps	[ecx+0x200],xmm6	/* Store hi real in t00+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x50]
		__asm	movaps		xmm1,[eax+0x50]
		__asm	unpckhpd	xmm7,[ebx+0x50]
		__asm	unpcklpd	xmm1,[ebx+0x50]	/* a[jp+p4 ] */
		__asm	movaps	[ecx+0x210],xmm7	/* Store hi imag in t01+32 */

		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy a[jt+p4] */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy a[jp+p4] */

		/* From here on, things are identical to the code in radix32_dif_pass: */
		__asm	mulpd	xmm0,[edx     ]	/* a[jt+p4]*c4 */
		__asm	mulpd	xmm1,[edx     ]	/* a[jp+p4]*c4 */
		__asm	mulpd	xmm2,[edx+0x10]	/* a[jt+p4]*s4 */
		__asm	mulpd	xmm3,[edx+0x10]	/* a[jp+p4]*s4 */
		__asm	addpd	xmm1,xmm2	/* xmm1 <- t01*/
		__asm	subpd	xmm0,xmm3	/* xmm0 <- t00*/
		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t00*/
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t01*/

		__asm	add	edx,0x020	/* c14 */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x140]	/* a[j1+p14], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x140]	/* a[j1+p14], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x140]	/* a[j2+p14] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x140]	/* a[jt+p14] */
		__asm	movaps	[ecx+0x220],xmm6	/* Store hi real in t11+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x150]
		__asm	movaps		xmm5,[eax+0x150]
		__asm	unpckhpd	xmm7,[ebx+0x150]
		__asm	unpcklpd	xmm5,[ebx+0x150]	/* a[jp+p14] */
		__asm	movaps	[ecx+0x230],xmm7	/* Store hi imag in t12+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p14] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p14] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p14]*c14 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p14]*c14 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p14]*s14 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p14]*s14 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */

		__asm	addpd	xmm0,xmm4	/* ~t13<- t13+rt */
		__asm	addpd	xmm1,xmm5	/* ~t14<- t14+it */
		__asm	subpd	xmm2,xmm4	/* ~t15<- t13-rt */
		__asm	subpd	xmm3,xmm5	/* ~t16<- t14-it	xmm4,5 free */

		/* Do the p08,p18 [hexadecimal here] combo - do p18 first so register assignments come out in same relative order as for p0,8 */
		__asm	add	edx,0x040	/* c1C */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x1c0]	/* a[j1+p1C], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x1c0]	/* a[j1+p1C], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x1c0]	/* a[j2+p1C] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x1c0]	/* a[jt+p1C] */
		__asm	movaps	[ecx+0x260],xmm6	/* Store hi real in t15+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x1d0]
		__asm	movaps		xmm5,[eax+0x1d0]
		__asm	unpckhpd	xmm7,[ebx+0x1d0]
		__asm	unpcklpd	xmm5,[ebx+0x1d0]	/* a[jp+p1C] */
		__asm	movaps	[ecx+0x270],xmm7	/* Store hi imag in t16+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p1C] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p1C] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p1C]*c1C */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p1C]*c1C */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p1C]*s1C */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p1C]*s1C */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */
		__asm	movaps	[ecx+0x010],xmm5	/* Store it */
		__asm	movaps	[ecx      ],xmm4	/* Store rt */

		__asm	sub	edx,0x020	/* c0C */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x0c0]	/* a[j1+p0C], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x0c0]	/* a[j1+p0C], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x0c0]	/* a[j2+p0C] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x0c0]	/* a[jt+p0C] */
		__asm	movaps	[ecx+0x240],xmm6	/* Store hi real in t13+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x0d0]
		__asm	movaps		xmm5,[eax+0x0d0]
		__asm	unpckhpd	xmm7,[ebx+0x0d0]
		__asm	unpcklpd	xmm5,[ebx+0x0d0]	/* a[jp+p0C] */
		__asm	movaps	[ecx+0x250],xmm7	/* Store hi imag in t14+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p0C] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p0C] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p0C]*c0C*/
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p0C]*c0C*/
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p0C]*s0C*/
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p0C]*s0C*/
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t14 */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t13 */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t14 */
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t13 */

		__asm	subpd	xmm4,[ecx      ]	/* ~t15<- t13-rt */
		__asm	subpd	xmm5,[ecx+0x010]	/* ~t16<- t14-it */
		__asm	addpd	xmm6,[ecx      ]	/* ~t13<- t13+rt */
		__asm	addpd	xmm7,[ecx+0x010]	/* ~t14<- t14+it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;
		~t6 =t2 -t6;		~t2 =t2 +t6;
		*/
		__asm	subpd	xmm0,xmm6	/*~t5 */
		__asm	subpd	xmm1,xmm7	/*~t6 */
		__asm	movaps	[ecx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */
		__asm	movaps	[ecx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */
		__asm	addpd	xmm6,xmm6	/* 2*t5 */
		__asm	addpd	xmm7,xmm7	/* 2*t6 */
		__asm	addpd	xmm6,xmm0	/*~t1 */
		__asm	addpd	xmm7,xmm1	/*~t2 */
		__asm	movaps	[ecx      ],xmm6	/* a[jt    ] <- ~t1 */
		__asm	movaps	[ecx+0x010],xmm7	/* a[jp    ] <- ~t2 */

		/*
		~t7 =t3 +t8;		~t3 =t3 -t8;
		~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	mov	edx, isrt2

		__asm	subpd	xmm2,xmm5	/*~t3 */
		__asm	subpd	xmm3,xmm4	/*~t8 */
		__asm	addpd	xmm5,xmm5	/* 2*t8 */
		__asm	addpd	xmm4,xmm4	/* 2*t7 */
		__asm	addpd	xmm5,xmm2	/*~t7 */
		__asm	addpd	xmm4,xmm3	/*~t4 */
		/*
		t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;
		t7 =(t7-t8)*ISRT2;	t8 =(t7+t8)*ISRT2;
		*/
		__asm	movaps	xmm6,xmm2	/* cpy t3 */
		__asm	movaps	xmm7,xmm5	/* cpy t7 */
		__asm	subpd	xmm2,xmm4	/* 3-4*/
		__asm	subpd	xmm5,xmm3	/* 7-8*/
		__asm	addpd	xmm6,xmm4	/* 3+4*/
		__asm	addpd	xmm7,xmm3	/* 7+8*/
		__asm	mulpd	xmm2,[edx]	/* (3-4)*ISRT2 */
		__asm	mulpd	xmm5,[edx]	/* (7-8)*ISRT2 */
		__asm	mulpd	xmm6,[edx]	/* (3+4)*ISRT2 */
		__asm	mulpd	xmm7,[edx]	/* (7+8)*ISRT2 */
		__asm	movaps	[ecx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[ecx+0x060],xmm5	/* a[jp+p12] <- ~t7 */
		__asm	movaps	[ecx+0x030],xmm6	/* a[jp+p4 ] <- ~t4 */
		__asm	movaps	[ecx+0x070],xmm7	/* a[jt+p12] <- ~t8 */

		/* Combine the 2 radix-4 subtransforms: */
	/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00)	*****/
		__asm	mov	esi, r00

		__asm	movaps	xmm0,[esi     ]	/* t0 */			__asm	movaps	xmm4,[esi+0x40]	/* t4 */
		__asm	movaps	xmm1,[esi+0x10]	/* t1 */			__asm	movaps	xmm5,[esi+0x50]	/* t5 */
		__asm	movaps	xmm2,[esi+0x80]	/* cpy t8 */		__asm	movaps	xmm7,[esi+0xd0]	/* td */
		__asm	movaps	xmm3,[esi+0x90]	/* cpy t9 */		__asm	movaps	xmm6,[esi+0xc0]	/* tc */
		__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* t4 = 4-d */
		__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* td = 5-c */
		__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */
		__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */
		__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* tc = 4+d */
		__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* t5 = 5+c */

		__asm	movaps	[esi+0x80],xmm0	/* t8 */			__asm	movaps	[esi+0x40],xmm4	/* t4 */
		__asm	movaps	[esi+0x90],xmm1	/* t9 */			__asm	movaps	[esi+0xd0],xmm5	/* td */
		__asm	movaps	[esi     ],xmm2	/* t0 */			__asm	movaps	[esi+0xc0],xmm7	/* tc */
		__asm	movaps	[esi+0x10],xmm3	/* t1 */			__asm	movaps	[esi+0x50],xmm6	/* t5 */

		__asm	movaps	xmm0,[esi+0x20]	/* t2 */			__asm	movaps	xmm4,[esi+0x60]	/* t6 */
		__asm	movaps	xmm1,[esi+0x30]	/* t3 */			__asm	movaps	xmm5,[esi+0x70]	/* t7 */
		__asm	movaps	xmm2,[esi+0xa0]	/* cpy ta */		__asm	movaps	xmm7,[esi+0xf0]	/* tf */
		__asm	movaps	xmm3,[esi+0xb0]	/* cpy tb */		__asm	movaps	xmm6,[esi+0xe0]	/* te */
		__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* t6 = 6-f */
		__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* tf = 7-e */
		__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */
		__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */
		__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* te = 6+f */
		__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* t7 = 7+e */

		__asm	movaps	[esi+0xa0],xmm0	/* ta */			__asm	movaps	[esi+0x60],xmm4	/* t6 */
		__asm	movaps	[esi+0xb0],xmm1	/* tb */			__asm	movaps	[esi+0xf0],xmm5	/* tf */
		__asm	movaps	[esi+0x20],xmm2	/* t2 */			__asm	movaps	[esi+0xe0],xmm7	/* te */
		__asm	movaps	[esi+0x30],xmm3	/* t3 */			__asm	movaps	[esi+0x70],xmm6	/* t7 */

	/*...Block 2:	*/
		__asm	add	eax,0x20	/* add0 += 4 */
		__asm	add	ebx,0x20	/* add1 += 4 */
	/*****	SSE2_RADIX4_DIF_4WRAPPER         (c02,c0A,c12,c1A,r10)	*****/
		__asm	mov	ecx,r10
		__asm	mov	edx,c02
		/* Do the p00,p10 combo: */
		/* For interleaved [j1,j2] version, replace [...] by the following: */
		/* Real parts: */
		__asm	movaps		xmm6,[eax     ]	/* a[j1+p0 ], this is the scratch xmm register */
		__asm	movaps		xmm0,[eax     ]	/* a[j1+p0 ], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx     ]	/* a[j2+p0 ] gets read twice */
		__asm	unpcklpd	xmm0,[ebx     ]	/* a[jt+p0 ] */
		__asm	movaps	[ecx+0x200],xmm6	/* Store hi real in t00+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x10]
		__asm	movaps		xmm1,[eax+0x10]
		__asm	unpckhpd	xmm7,[ebx+0x10]
		__asm	unpcklpd	xmm1,[ebx+0x10]	/* a[jp+p0 ] */
		__asm	movaps	[ecx+0x210],xmm7	/* Store hi imag in t01+32 */

		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy a[jt+p0] */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy a[jp+p0] */

		/* From here on, things are identical to the code in radix32_dif_pass: */
		__asm	mulpd	xmm0,[edx     ]	/* a[jt+p0]*c0 */
		__asm	mulpd	xmm1,[edx     ]	/* a[jp+p0]*c0 */
		__asm	mulpd	xmm2,[edx+0x10]	/* a[jt+p0]*s0 */
		__asm	mulpd	xmm3,[edx+0x10]	/* a[jp+p0]*s0 */
		__asm	addpd	xmm1,xmm2	/* xmm1 <- t01*/
		__asm	subpd	xmm0,xmm3	/* xmm0 <- t00*/
		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t00*/
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t01*/

		__asm	add	edx,0x020	/* c10 */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x100]	/* a[j1+p10], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x100]	/* a[j1+p10], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x100]	/* a[j2+p10] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x100]	/* a[jt+p10] */
		__asm	movaps	[ecx+0x220],xmm6	/* Store hi real in t11+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x110]
		__asm	movaps		xmm5,[eax+0x110]
		__asm	unpckhpd	xmm7,[ebx+0x110]
		__asm	unpcklpd	xmm5,[ebx+0x110]	/* a[jp+p10] */
		__asm	movaps	[ecx+0x230],xmm7	/* Store hi imag in t12+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p10] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p10] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p10]*c10 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p10]*c10 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p10]*s10 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p10]*s10 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */

		__asm	addpd	xmm0,xmm4	/* ~t1 <- t1 +rt */
		__asm	addpd	xmm1,xmm5	/* ~t2 <- t2 +it */
		__asm	subpd	xmm2,xmm4	/* ~t3 <- t1 -rt */
		__asm	subpd	xmm3,xmm5	/* ~t4 <- t2 -it	xmm4,5 free */

		/* Do the p08,p18 [hexadecimal here] combo - do p18 first so register assignments come out in same relative order as for p0,8 */
		__asm	add	edx,0x040	/* c18 */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x180]	/* a[j1+p18], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x180]	/* a[j1+p18], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x180]	/* a[j2+p18] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x180]	/* a[jt+p18] */
		__asm	movaps	[ecx+0x260],xmm6	/* Store hi real in t15+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x190]
		__asm	movaps		xmm5,[eax+0x190]
		__asm	unpckhpd	xmm7,[ebx+0x190]
		__asm	unpcklpd	xmm5,[ebx+0x190]	/* a[jp+p18] */
		__asm	movaps	[ecx+0x270],xmm7	/* Store hi imag in t16+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p18] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p18] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p18]*c18 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p18]*c18 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p18]*s18 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p18]*s18 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */
		__asm	movaps	[ecx+0x010],xmm5	/* Store it */
		__asm	movaps	[ecx      ],xmm4	/* Store rt */

		__asm	sub	edx,0x020	/* c08 */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x080]	/* a[j1+p08], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x080]	/* a[j1+p08], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x080]	/* a[j2+p08] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x080]	/* a[jt+p08] */
		__asm	movaps	[ecx+0x240],xmm6	/* Store hi real in t13+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x090]
		__asm	movaps		xmm5,[eax+0x090]
		__asm	unpckhpd	xmm7,[ebx+0x090]
		__asm	unpcklpd	xmm5,[ebx+0x090]	/* a[jp+p08] */
		__asm	movaps	[ecx+0x250],xmm7	/* Store hi imag in t14+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p08] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p08] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p08]*c08 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p08]*c08 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p08]*s08 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p08]*s08 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */

		__asm	subpd	xmm4,[ecx      ]	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm5,[ecx+0x010]	/* ~t8 <- t6 -it */
		__asm	addpd	xmm6,[ecx      ]	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm7,[ecx+0x010]	/* ~t6 <- t6 +it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;
		~t6 =t2 -t6;		~t2 =t2 +t6;
		*/
		__asm	subpd	xmm0,xmm6	/*~t5 */
		__asm	subpd	xmm1,xmm7	/*~t6 */
		__asm	movaps	[ecx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */
		__asm	movaps	[ecx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */
		__asm	addpd	xmm6,xmm6	/* 2*t5 */
		__asm	addpd	xmm7,xmm7	/* 2*t6 */
		__asm	addpd	xmm6,xmm0	/*~t1 */
		__asm	addpd	xmm7,xmm1	/*~t2 */
		__asm	movaps	[ecx      ],xmm6	/* a[jt    ] <- ~t1 */
		__asm	movaps	[ecx+0x010],xmm7	/* a[jp    ] <- ~t2 */

		/*
		~t7 =t3 +t8;		~t3 =t3 -t8;
		~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	subpd	xmm2,xmm5	/*~t3 */
		__asm	subpd	xmm3,xmm4	/*~t8 */
		__asm	movaps	[ecx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[ecx+0x070],xmm3	/* a[jp+p12] <- ~t8 */
		__asm	addpd	xmm5,xmm5	/* 2*t8 */
		__asm	addpd	xmm4,xmm4	/* 2*t7 */
		__asm	addpd	xmm5,xmm2	/*~t7 */
		__asm	addpd	xmm4,xmm3	/*~t4 */
		__asm	movaps	[ecx+0x060],xmm5	/* a[jt+p12] <- ~t7 */
		__asm	movaps	[ecx+0x030],xmm4	/* a[jp+p4 ] <- ~t4 */

	/*****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/
		__asm	add	ecx,0x080	/* r18 */
		__asm	add	edx,0x040	/* c06 */

		/* Do the p00,p10 combo: */
		/* For interleaved [j1,j2] version, replace [...] by the following:
		*/
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x40]	/* a[j1+p4 ], this is the scratch xmm register */
		__asm	movaps		xmm0,[eax+0x40]	/* a[j1+p4 ], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x40]	/* a[j2+p4 ] gets read twice */
		__asm	unpcklpd	xmm0,[ebx+0x40]	/* a[jt+p4 ] */
		__asm	movaps	[ecx+0x200],xmm6	/* Store hi real in t00+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x50]
		__asm	movaps		xmm1,[eax+0x50]
		__asm	unpckhpd	xmm7,[ebx+0x50]
		__asm	unpcklpd	xmm1,[ebx+0x50]	/* a[jp+p4 ] */
		__asm	movaps	[ecx+0x210],xmm7	/* Store hi imag in t01+32 */

		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy a[jt+p4] */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy a[jp+p4] */

		/* From here on, things are identical to the code in radix32_dif_pass: */
		__asm	mulpd	xmm0,[edx     ]	/* a[jt+p4]*c4 */
		__asm	mulpd	xmm1,[edx     ]	/* a[jp+p4]*c4 */
		__asm	mulpd	xmm2,[edx+0x10]	/* a[jt+p4]*s4 */
		__asm	mulpd	xmm3,[edx+0x10]	/* a[jp+p4]*s4 */
		__asm	addpd	xmm1,xmm2	/* xmm1 <- t01*/
		__asm	subpd	xmm0,xmm3	/* xmm0 <- t00*/
		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t00*/
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t01*/

		__asm	add	edx,0x020	/* c14 */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x140]	/* a[j1+p14], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x140]	/* a[j1+p14], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x140]	/* a[j2+p14] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x140]	/* a[jt+p14] */
		__asm	movaps	[ecx+0x220],xmm6	/* Store hi real in t11+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x150]
		__asm	movaps		xmm5,[eax+0x150]
		__asm	unpckhpd	xmm7,[ebx+0x150]
		__asm	unpcklpd	xmm5,[ebx+0x150]	/* a[jp+p14] */
		__asm	movaps	[ecx+0x230],xmm7	/* Store hi imag in t12+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p14] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p14] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p14]*c14 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p14]*c14 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p14]*s14 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p14]*s14 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */

		__asm	addpd	xmm0,xmm4	/* ~t13<- t13+rt */
		__asm	addpd	xmm1,xmm5	/* ~t14<- t14+it */
		__asm	subpd	xmm2,xmm4	/* ~t15<- t13-rt */
		__asm	subpd	xmm3,xmm5	/* ~t16<- t14-it	xmm4,5 free */

		/* Do the p08,p18 [hexadecimal here] combo - do p18 first so register assignments come out in same relative order as for p0,8 */
		__asm	add	edx,0x040	/* c1C */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x1c0]	/* a[j1+p1C], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x1c0]	/* a[j1+p1C], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x1c0]	/* a[j2+p1C] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x1c0]	/* a[jt+p1C] */
		__asm	movaps	[ecx+0x260],xmm6	/* Store hi real in t15+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x1d0]
		__asm	movaps		xmm5,[eax+0x1d0]
		__asm	unpckhpd	xmm7,[ebx+0x1d0]
		__asm	unpcklpd	xmm5,[ebx+0x1d0]	/* a[jp+p1C] */
		__asm	movaps	[ecx+0x270],xmm7	/* Store hi imag in t16+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p1C] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p1C] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p1C]*c1C */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p1C]*c1C */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p1C]*s1C */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p1C]*s1C */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */
		__asm	movaps	[ecx+0x010],xmm5	/* Store it */
		__asm	movaps	[ecx      ],xmm4	/* Store rt */

		__asm	sub	edx,0x020	/* c0C */
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x0c0]	/* a[j1+p0C], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x0c0]	/* a[j1+p0C], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x0c0]	/* a[j2+p0C] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x0c0]	/* a[jt+p0C] */
		__asm	movaps	[ecx+0x240],xmm6	/* Store hi real in t13+32 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x0d0]
		__asm	movaps		xmm5,[eax+0x0d0]
		__asm	unpckhpd	xmm7,[ebx+0x0d0]
		__asm	unpcklpd	xmm5,[ebx+0x0d0]	/* a[jp+p0C] */
		__asm	movaps	[ecx+0x250],xmm7	/* Store hi imag in t14+32 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p0C] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p0C] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p0C]*c0C*/
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p0C]*c0C*/
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p0C]*s0C*/
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p0C]*s0C*/
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t14 */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t13 */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t14 */
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t13 */

		__asm	subpd	xmm4,[ecx      ]	/* ~t15<- t13-rt */
		__asm	subpd	xmm5,[ecx+0x010]	/* ~t16<- t14-it */
		__asm	addpd	xmm6,[ecx      ]	/* ~t13<- t13+rt */
		__asm	addpd	xmm7,[ecx+0x010]	/* ~t14<- t14+it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;
		~t6 =t2 -t6;		~t2 =t2 +t6;
		*/
		__asm	subpd	xmm0,xmm6	/*~t5 */
		__asm	subpd	xmm1,xmm7	/*~t6 */
		__asm	movaps	[ecx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */
		__asm	movaps	[ecx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */
		__asm	addpd	xmm6,xmm6	/* 2*t5 */
		__asm	addpd	xmm7,xmm7	/* 2*t6 */
		__asm	addpd	xmm6,xmm0	/*~t1 */
		__asm	addpd	xmm7,xmm1	/*~t2 */
		__asm	movaps	[ecx      ],xmm6	/* a[jt    ] <- ~t1 */
		__asm	movaps	[ecx+0x010],xmm7	/* a[jp    ] <- ~t2 */

		/*
		~t7 =t3 +t8;		~t3 =t3 -t8;
		~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	mov	edx, isrt2

		__asm	subpd	xmm2,xmm5	/*~t3 */
		__asm	subpd	xmm3,xmm4	/*~t8 */
		__asm	addpd	xmm5,xmm5	/* 2*t8 */
		__asm	addpd	xmm4,xmm4	/* 2*t7 */
		__asm	addpd	xmm5,xmm2	/*~t7 */
		__asm	addpd	xmm4,xmm3	/*~t4 */
		/*
		t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;
		t7 =(t7-t8)*ISRT2;	t8 =(t7+t8)*ISRT2;
		*/
		__asm	movaps	xmm6,xmm2	/* cpy t3 */
		__asm	movaps	xmm7,xmm5	/* cpy t7 */
		__asm	subpd	xmm2,xmm4	/* 3-4*/
		__asm	subpd	xmm5,xmm3	/* 7-8*/
		__asm	addpd	xmm6,xmm4	/* 3+4*/
		__asm	addpd	xmm7,xmm3	/* 7+8*/
		__asm	mulpd	xmm2,[edx]	/* (3-4)*ISRT2 */
		__asm	mulpd	xmm5,[edx]	/* (7-8)*ISRT2 */
		__asm	mulpd	xmm6,[edx]	/* (3+4)*ISRT2 */
		__asm	mulpd	xmm7,[edx]	/* (7+8)*ISRT2 */
		__asm	movaps	[ecx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[ecx+0x060],xmm5	/* a[jp+p12] <- ~t7 */
		__asm	movaps	[ecx+0x030],xmm6	/* a[jp+p4 ] <- ~t4 */
		__asm	movaps	[ecx+0x070],xmm7	/* a[jt+p12] <- ~t8 */

		/* Combine the 2 radix-4 subtransforms: */
	/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10)	*****/
		__asm	mov	esi, r10

		__asm	movaps	xmm0,[esi     ]	/* t0 */			__asm	movaps	xmm4,[esi+0x40]	/* t4 */
		__asm	movaps	xmm1,[esi+0x10]	/* t1 */			__asm	movaps	xmm5,[esi+0x50]	/* t5 */
		__asm	movaps	xmm2,[esi+0x80]	/* cpy t8 */		__asm	movaps	xmm7,[esi+0xd0]	/* td */
		__asm	movaps	xmm3,[esi+0x90]	/* cpy t9 */		__asm	movaps	xmm6,[esi+0xc0]	/* tc */
		__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* t4 = 4-d */
		__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* td = 5-c */
		__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */
		__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */
		__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* tc = 4+d */
		__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* t5 = 5+c */

		__asm	movaps	[esi+0x80],xmm0	/* t8 */			__asm	movaps	[esi+0x40],xmm4	/* t4 */
		__asm	movaps	[esi+0x90],xmm1	/* t9 */			__asm	movaps	[esi+0xd0],xmm5	/* td */
		__asm	movaps	[esi     ],xmm2	/* t0 */			__asm	movaps	[esi+0xc0],xmm7	/* tc */
		__asm	movaps	[esi+0x10],xmm3	/* t1 */			__asm	movaps	[esi+0x50],xmm6	/* t5 */

		__asm	movaps	xmm0,[esi+0x20]	/* t2 */			__asm	movaps	xmm4,[esi+0x60]	/* t6 */
		__asm	movaps	xmm1,[esi+0x30]	/* t3 */			__asm	movaps	xmm5,[esi+0x70]	/* t7 */
		__asm	movaps	xmm2,[esi+0xa0]	/* cpy ta */		__asm	movaps	xmm7,[esi+0xf0]	/* tf */
		__asm	movaps	xmm3,[esi+0xb0]	/* cpy tb */		__asm	movaps	xmm6,[esi+0xe0]	/* te */
		__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* t6 = 6-f */
		__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* tf = 7-e */
		__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */
		__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */
		__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* te = 6+f */
		__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* t7 = 7+e */

		__asm	movaps	[esi+0xa0],xmm0	/* ta */			__asm	movaps	[esi+0x60],xmm4	/* t6 */
		__asm	movaps	[esi+0xb0],xmm1	/* tb */			__asm	movaps	[esi+0xf0],xmm5	/* tf */
		__asm	movaps	[esi+0x20],xmm2	/* t2 */			__asm	movaps	[esi+0xe0],xmm7	/* te */
		__asm	movaps	[esi+0x30],xmm3	/* t3 */			__asm	movaps	[esi+0x70],xmm6	/* t7 */

	/***************************************************************************************************************************
	Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks
	[operating on the odd-indexed elements from the unpck*pd commands which were stored to temporaries can use a common macro:
	***************************************************************************************************************************/
	/*...Block 3:	*/
	/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */
		__asm	add	ecx,0x080	/* r20 */
		/* Do the p0,p8 combo: */
		__asm	mov	ebx, c01
		__asm	mov	eax, ecx	/* r20 */
		__asm	add	ecx, 0x020	/* r22 */

		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */
		__asm	movaps	xmm6,[ebx     ]	/* c0 */
		__asm	movaps	xmm7,[ebx+0x10]	/* s0 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */

		__asm	mulpd	xmm0,xmm6		/* a[jt   ]*c0 */
		__asm	mulpd	xmm1,xmm6		/* a[jp   ]*c0 */
		__asm	mulpd	xmm2,xmm7		/* a[jt   ]*s0 */
		__asm	mulpd	xmm3,xmm7		/* a[jp   ]*s0	xmm6,7 free */
																	__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */
		__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */
																	__asm	mulpd	xmm4,[ebx+0x20] /* a[jt+p8 ]*c8 */
		__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[ebx+0x20]	/* a[jp+p8 ]*c8 */
																	__asm	mulpd	xmm6,[ebx+0x30]	/* a[jt+p8 ]*s8 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[ebx+0x30]	/* a[jp+p8 ]*s8 */
																	__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */

																	__asm	add	ecx, 0x040	/* r26 */
																	__asm	add	ebx, 0x60
																	__asm	movaps	xmm6,[ecx     ]	/* a[jt+p12] */
																	__asm	movaps	xmm7,[ecx+0x10]	/* a[jp+p12] */

		__asm	addpd	xmm0,xmm4		/* ~t1 <- t1 +rt */
		__asm	addpd	xmm1,xmm5		/* ~t2 <- t2 +it */
		__asm	subpd	xmm2,xmm4		/* ~t3 <- t1 -rt */
		__asm	subpd	xmm3,xmm5		/* ~t4 <- t2 -it	xmm4,5 free */

		/* Do the p4,12 combo: */
		__asm	movaps	xmm4,xmm6		/* xmm4 <- cpy a[jt+p12] */
		__asm	movaps	xmm5,xmm7		/* xmm5 <- cpy a[jp+p12] */

		__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p12]*c12 */
		__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p12]*c12 */
		__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p12]*s12 */
		__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p12]*s12 */
		__asm	mov	edx, eax	/* r20 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */
		__asm	movaps	[edx+0x010],xmm5	/* store it */
		__asm	movaps	[edx      ],xmm4	/* store rt */

		__asm	add	eax, 0x040	/* r24 */
		__asm	sub	ebx, 0x20
		__asm	movaps	xmm4,[eax     ]	/* a[jt+p4] */
		__asm	movaps	xmm5,[eax+0x10]	/* a[jp+p4] */
		__asm	movaps	xmm6,xmm4		/* xmm4 <- cpy a[jt+p4] */
		__asm	movaps	xmm7,xmm5		/* xmm5 <- cpy a[jp+p4] */

		__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p4]*c4 */
		__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p4]*c4 */
		__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p4]*s4 */
		__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p4]*s4 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */

		__asm	subpd	xmm4,[edx      ]	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm5,[edx+0x010]	/* ~t8 <- t6 -it */
		__asm	addpd	xmm6,[edx      ]	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm7,[edx+0x010]	/* ~t6 <- t6 +it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;
		~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	subpd	xmm0,xmm6	/*~t5 */									__asm	subpd	xmm2,xmm5	/*~t3 */
		__asm	subpd	xmm1,xmm7	/*~t6 */									__asm	subpd	xmm3,xmm4	/*~t8 */
		__asm	movaps	[edx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	movaps	[edx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[edx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	movaps	[edx+0x070],xmm3	/* a[jp+p12] <- ~t8 */
		__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */
		__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */
		__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */
		__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */
		__asm	movaps	[edx      ],xmm6	/* a[jt    ] <- ~t1 */				__asm	movaps	[edx+0x060],xmm5	/* a[jt+p12] <- ~t7 */
		__asm	movaps	[edx+0x010],xmm7	/* a[jp    ] <- ~t2 */				__asm	movaps	[edx+0x030],xmm4	/* a[jp+p4 ] <- ~t4 */

	/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/
		/* Do the p0,p8 combo: */
		__asm	mov	ebx, c05
		__asm	add	eax, 0x040	/* r28 */
		__asm	add	ecx, 0x040	/* r2A */
		__asm	mov	edx, eax	/* r28 */

		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */
		__asm	movaps	xmm6,[ebx     ]	/* c0 */
		__asm	movaps	xmm7,[ebx+0x10]	/* s0 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */

		__asm	mulpd	xmm0,xmm6		/* a[jt   ]*c0 */
		__asm	mulpd	xmm1,xmm6		/* a[jp   ]*c0 */
		__asm	mulpd	xmm2,xmm7		/* a[jt   ]*s0 */
		__asm	mulpd	xmm3,xmm7		/* a[jp   ]*s0	xmm6,7 free */
																	__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */
		__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */
																	__asm	mulpd	xmm4,[ebx+0x20] /* a[jt+p8 ]*c8 */
		__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[ebx+0x20]	/* a[jp+p8 ]*c8 */
																	__asm	mulpd	xmm6,[ebx+0x30]	/* a[jt+p8 ]*s8 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[ebx+0x30]	/* a[jp+p8 ]*s8 */
																	__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */

		__asm	addpd	xmm0,xmm4	/* ~t1 <- t1 +rt */
		__asm	addpd	xmm1,xmm5	/* ~t2 <- t2 +it */
		__asm	subpd	xmm2,xmm4	/* ~t3 <- t1 -rt */
		__asm	subpd	xmm3,xmm5	/* ~t4 <- t2 -it	xmm4,5 free */

		/* Do the p4,12 combo: */
		__asm	add	ebx, 0x60
		__asm	add	eax, 0x040	/* r2C */
		__asm	add	ecx, 0x040	/* r2E */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p12] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p12] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p12] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p12] */

		__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p12]*c12 */
		__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p12]*c12 */
		__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p12]*s12 */
		__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p12]*s12 */

		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */
		__asm	movaps	[edx+0x010],xmm5	/* tmp store it */
		__asm	movaps	[edx      ],xmm4	/* tmp store rt */

		__asm	sub	ebx, 0x20
		__asm	movaps	xmm4,[eax     ]	/* a[jt+p4] */
		__asm	movaps	xmm5,[eax+0x10]	/* a[jp+p4] */
		__asm	movaps	xmm6,[eax     ]	/* xmm2 <- cpy a[jt+p4] */
		__asm	movaps	xmm7,[eax+0x10]	/* xmm3 <- cpy a[jp+p4] */

		__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p4]*c4 */
		__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p4]*c4 */
		__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p4]*s4 */
		__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p4]*s4 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */

		__asm	subpd	xmm4,[edx      ]	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm5,[edx+0x010]	/* ~t8 <- t6 -it */
		__asm	addpd	xmm6,[edx      ]	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm7,[edx+0x010]	/* ~t6 <- t6 +it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;
		~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	subpd	xmm0,xmm6	/*~t5 */									__asm	mov	ecx, isrt2
		__asm	subpd	xmm1,xmm7	/*~t6 */
		__asm	movaps	[edx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	subpd	xmm2,xmm5	/*~t3 */
		__asm	movaps	[edx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	subpd	xmm3,xmm4	/*~t8 */
		__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */
		__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */
		__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */
		__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */
		__asm	movaps	[edx      ],xmm6	/* a[jt    ] <- ~t1 */
		__asm	movaps	[edx+0x010],xmm7	/* a[jp    ] <- ~t2 */
		/*
		t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;
		t7 =(t7-t8)*ISRT2;	t8 =(t7+t8)*ISRT2;
		*/
		__asm	movaps	xmm0,[ecx]	/* ISRT2 */
		__asm	movaps	xmm6,xmm2	/* cpy t3 */
		__asm	movaps	xmm7,xmm5	/* cpy t7 */
		__asm	subpd	xmm2,xmm4	/* 3-4*/
		__asm	subpd	xmm5,xmm3	/* 7-8*/
		__asm	addpd	xmm6,xmm4	/* 3+4*/
		__asm	addpd	xmm7,xmm3	/* 7+8*/
		__asm	mulpd	xmm2,xmm0	/* (3-4)*ISRT2 */
		__asm	mulpd	xmm5,xmm0	/* (7-8)*ISRT2 */
		__asm	mulpd	xmm6,xmm0	/* (3+4)*ISRT2 */
		__asm	mulpd	xmm7,xmm0	/* (7+8)*ISRT2 */
		__asm	movaps	[edx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[edx+0x060],xmm5	/* a[jp+p12] <- ~t7 */
		__asm	movaps	[edx+0x030],xmm6	/* a[jp+p4 ] <- ~t4 */
		__asm	movaps	[edx+0x070],xmm7	/* a[jt+p12] <- ~t8 */

		/* Combine the 2 radix-4 subtransforms: */
	/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)	*****/
		__asm	mov	esi, r20

		__asm	movaps	xmm0,[esi     ]	/* t0 */			__asm	movaps	xmm4,[esi+0x40]	/* t4 */
		__asm	movaps	xmm1,[esi+0x10]	/* t1 */			__asm	movaps	xmm5,[esi+0x50]	/* t5 */
		__asm	movaps	xmm2,[esi+0x80]	/* cpy t8 */		__asm	movaps	xmm7,[esi+0xd0]	/* td */
		__asm	movaps	xmm3,[esi+0x90]	/* cpy t9 */		__asm	movaps	xmm6,[esi+0xc0]	/* tc */
		__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* t4 = 4-d */
		__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* td = 5-c */
		__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */
		__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */
		__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* tc = 4+d */
		__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* t5 = 5+c */

		__asm	movaps	[esi+0x80],xmm0	/* t8 */			__asm	movaps	[esi+0x40],xmm4	/* t4 */
		__asm	movaps	[esi+0x90],xmm1	/* t9 */			__asm	movaps	[esi+0xd0],xmm5	/* td */
		__asm	movaps	[esi     ],xmm2	/* t0 */			__asm	movaps	[esi+0xc0],xmm7	/* tc */
		__asm	movaps	[esi+0x10],xmm3	/* t1 */			__asm	movaps	[esi+0x50],xmm6	/* t5 */

		__asm	movaps	xmm0,[esi+0x20]	/* t2 */			__asm	movaps	xmm4,[esi+0x60]	/* t6 */
		__asm	movaps	xmm1,[esi+0x30]	/* t3 */			__asm	movaps	xmm5,[esi+0x70]	/* t7 */
		__asm	movaps	xmm2,[esi+0xa0]	/* cpy ta */		__asm	movaps	xmm7,[esi+0xf0]	/* tf */
		__asm	movaps	xmm3,[esi+0xb0]	/* cpy tb */		__asm	movaps	xmm6,[esi+0xe0]	/* te */
		__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* t6 = 6-f */
		__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* tf = 7-e */
		__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */
		__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */
		__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* te = 6+f */
		__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* t7 = 7+e */

		__asm	movaps	[esi+0xa0],xmm0	/* ta */			__asm	movaps	[esi+0x60],xmm4	/* t6 */
		__asm	movaps	[esi+0xb0],xmm1	/* tb */			__asm	movaps	[esi+0xf0],xmm5	/* tf */
		__asm	movaps	[esi+0x20],xmm2	/* t2 */			__asm	movaps	[esi+0xe0],xmm7	/* te */
		__asm	movaps	[esi+0x30],xmm3	/* t3 */			__asm	movaps	[esi+0x70],xmm6	/* t7 */

	/*...Block 4:	*/
	/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/
		/* Do the p0,p8 combo: */
		__asm	mov	ebx, c03
		__asm	mov	eax, r30	/* r30 */
		__asm	mov	ecx, eax
		__asm	add	ecx, 0x020	/* r32 */

		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */
		__asm	movaps	xmm6,[ebx     ]	/* c0 */
		__asm	movaps	xmm7,[ebx+0x10]	/* s0 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */

		__asm	mulpd	xmm0,xmm6		/* a[jt   ]*c0 */
		__asm	mulpd	xmm1,xmm6		/* a[jp   ]*c0 */
		__asm	mulpd	xmm2,xmm7		/* a[jt   ]*s0 */
		__asm	mulpd	xmm3,xmm7		/* a[jp   ]*s0	xmm6,7 free */
																	__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */
		__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */
																	__asm	mulpd	xmm4,[ebx+0x20] /* a[jt+p8 ]*c8 */
		__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[ebx+0x20]	/* a[jp+p8 ]*c8 */
																	__asm	mulpd	xmm6,[ebx+0x30]	/* a[jt+p8 ]*s8 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[ebx+0x30]	/* a[jp+p8 ]*s8 */
																	__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */

																	__asm	add	ecx, 0x040	/* r36 */
																	__asm	add	ebx, 0x60
																	__asm	movaps	xmm6,[ecx     ]	/* a[jt+p12] */
																	__asm	movaps	xmm7,[ecx+0x10]	/* a[jp+p12] */

		__asm	addpd	xmm0,xmm4		/* ~t1 <- t1 +rt */
		__asm	addpd	xmm1,xmm5		/* ~t2 <- t2 +it */
		__asm	subpd	xmm2,xmm4		/* ~t3 <- t1 -rt */
		__asm	subpd	xmm3,xmm5		/* ~t4 <- t2 -it	xmm4,5 free */

		/* Do the p4,12 combo: */
		__asm	movaps	xmm4,xmm6		/* xmm4 <- cpy a[jt+p12] */
		__asm	movaps	xmm5,xmm7		/* xmm5 <- cpy a[jp+p12] */

		__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p12]*c12 */
		__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p12]*c12 */
		__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p12]*s12 */
		__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p12]*s12 */
		__asm	mov	edx, eax	/* r30 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */
		__asm	movaps	[edx+0x010],xmm5	/* store it */
		__asm	movaps	[edx      ],xmm4	/* store rt */

		__asm	add	eax, 0x040	/* r34 */
		__asm	sub	ebx, 0x20
		__asm	movaps	xmm4,[eax     ]	/* a[jt+p4] */
		__asm	movaps	xmm5,[eax+0x10]	/* a[jp+p4] */
		__asm	movaps	xmm6,xmm4		/* xmm4 <- cpy a[jt+p4] */
		__asm	movaps	xmm7,xmm5		/* xmm5 <- cpy a[jp+p4] */

		__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p4]*c4 */
		__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p4]*c4 */
		__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p4]*s4 */
		__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p4]*s4 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */

		__asm	subpd	xmm4,[edx      ]	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm5,[edx+0x010]	/* ~t8 <- t6 -it */
		__asm	addpd	xmm6,[edx      ]	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm7,[edx+0x010]	/* ~t6 <- t6 +it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;
		~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	subpd	xmm0,xmm6	/*~t5 */									__asm	subpd	xmm2,xmm5	/*~t3 */
		__asm	subpd	xmm1,xmm7	/*~t6 */									__asm	subpd	xmm3,xmm4	/*~t8 */
		__asm	movaps	[edx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	movaps	[edx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[edx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	movaps	[edx+0x070],xmm3	/* a[jp+p12] <- ~t8 */
		__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */
		__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */
		__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */
		__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */
		__asm	movaps	[edx      ],xmm6	/* a[jt    ] <- ~t1 */				__asm	movaps	[edx+0x060],xmm5	/* a[jt+p12] <- ~t7 */
		__asm	movaps	[edx+0x010],xmm7	/* a[jp    ] <- ~t2 */				__asm	movaps	[edx+0x030],xmm4	/* a[jp+p4 ] <- ~t4 */

	/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/
		/* Do the p0,p8 combo: */
		__asm	mov	ebx, c07
		__asm	add	eax, 0x040	/* r38 */
		__asm	add	ecx, 0x040	/* r3A */
		__asm	mov	edx, eax	/* r38 */

		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */
		__asm	movaps	xmm6,[ebx     ]	/* c0 */
		__asm	movaps	xmm7,[ebx+0x10]	/* s0 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */

		__asm	mulpd	xmm0,xmm6		/* a[jt   ]*c0 */
		__asm	mulpd	xmm1,xmm6		/* a[jp   ]*c0 */
		__asm	mulpd	xmm2,xmm7		/* a[jt   ]*s0 */
		__asm	mulpd	xmm3,xmm7		/* a[jp   ]*s0	xmm6,7 free */
																	__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */
		__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */
																	__asm	mulpd	xmm4,[ebx+0x20] /* a[jt+p8 ]*c8 */
		__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[ebx+0x20]	/* a[jp+p8 ]*c8 */
																	__asm	mulpd	xmm6,[ebx+0x30]	/* a[jt+p8 ]*s8 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[ebx+0x30]	/* a[jp+p8 ]*s8 */
																	__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */

		__asm	addpd	xmm0,xmm4	/* ~t1 <- t1 +rt */
		__asm	addpd	xmm1,xmm5	/* ~t2 <- t2 +it */
		__asm	subpd	xmm2,xmm4	/* ~t3 <- t1 -rt */
		__asm	subpd	xmm3,xmm5	/* ~t4 <- t2 -it	xmm4,5 free */

		/* Do the p4,12 combo: */
		__asm	add	ebx, 0x60
		__asm	add	eax, 0x040	/* r3C */
		__asm	add	ecx, 0x040	/* r3E */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p12] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p12] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p12] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p12] */

		__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p12]*c12 */
		__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p12]*c12 */
		__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p12]*s12 */
		__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p12]*s12 */

		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */
		__asm	movaps	[edx+0x010],xmm5	/* tmp store it */
		__asm	movaps	[edx      ],xmm4	/* tmp store rt */

		__asm	sub	ebx, 0x20
		__asm	movaps	xmm4,[eax     ]	/* a[jt+p4] */
		__asm	movaps	xmm5,[eax+0x10]	/* a[jp+p4] */
		__asm	movaps	xmm6,[eax     ]	/* xmm2 <- cpy a[jt+p4] */
		__asm	movaps	xmm7,[eax+0x10]	/* xmm3 <- cpy a[jp+p4] */

		__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p4]*c4 */
		__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p4]*c4 */
		__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p4]*s4 */
		__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p4]*s4 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */

		__asm	subpd	xmm4,[edx      ]	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm5,[edx+0x010]	/* ~t8 <- t6 -it */
		__asm	addpd	xmm6,[edx      ]	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm7,[edx+0x010]	/* ~t6 <- t6 +it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;
		~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	subpd	xmm0,xmm6	/*~t5 */									__asm	mov	ecx, isrt2
		__asm	subpd	xmm1,xmm7	/*~t6 */
		__asm	movaps	[edx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	subpd	xmm2,xmm5	/*~t3 */
		__asm	movaps	[edx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	subpd	xmm3,xmm4	/*~t8 */
		__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */
		__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */
		__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */
		__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */
		__asm	movaps	[edx      ],xmm6	/* a[jt    ] <- ~t1 */
		__asm	movaps	[edx+0x010],xmm7	/* a[jp    ] <- ~t2 */
		/*
		t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;
		t7 =(t7-t8)*ISRT2;	t8 =(t7+t8)*ISRT2;
		*/
		__asm	movaps	xmm0,[ecx]	/* ISRT2 */
		__asm	movaps	xmm6,xmm2	/* cpy t3 */
		__asm	movaps	xmm7,xmm5	/* cpy t7 */
		__asm	subpd	xmm2,xmm4	/* 3-4*/
		__asm	subpd	xmm5,xmm3	/* 7-8*/
		__asm	addpd	xmm6,xmm4	/* 3+4*/
		__asm	addpd	xmm7,xmm3	/* 7+8*/
		__asm	mulpd	xmm2,xmm0	/* (3-4)*ISRT2 */
		__asm	mulpd	xmm5,xmm0	/* (7-8)*ISRT2 */
		__asm	mulpd	xmm6,xmm0	/* (3+4)*ISRT2 */
		__asm	mulpd	xmm7,xmm0	/* (7+8)*ISRT2 */
		__asm	movaps	[edx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[edx+0x060],xmm5	/* a[jp+p12] <- ~t7 */
		__asm	movaps	[edx+0x030],xmm6	/* a[jp+p4 ] <- ~t4 */
		__asm	movaps	[edx+0x070],xmm7	/* a[jt+p12] <- ~t8 */

		/* Combine the 2 radix-4 subtransforms: */
	/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30)	*****/
		__asm	mov	esi, r30

		__asm	movaps	xmm0,[esi     ]	/* t0 */			__asm	movaps	xmm4,[esi+0x40]	/* t4 */
		__asm	movaps	xmm1,[esi+0x10]	/* t1 */			__asm	movaps	xmm5,[esi+0x50]	/* t5 */
		__asm	movaps	xmm2,[esi+0x80]	/* cpy t8 */		__asm	movaps	xmm7,[esi+0xd0]	/* td */
		__asm	movaps	xmm3,[esi+0x90]	/* cpy t9 */		__asm	movaps	xmm6,[esi+0xc0]	/* tc */
		__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* t4 = 4-d */
		__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* td = 5-c */
		__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */
		__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */
		__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* tc = 4+d */
		__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* t5 = 5+c */

		__asm	movaps	[esi+0x80],xmm0	/* t8 */			__asm	movaps	[esi+0x40],xmm4	/* t4 */
		__asm	movaps	[esi+0x90],xmm1	/* t9 */			__asm	movaps	[esi+0xd0],xmm5	/* td */
		__asm	movaps	[esi     ],xmm2	/* t0 */			__asm	movaps	[esi+0xc0],xmm7	/* tc */
		__asm	movaps	[esi+0x10],xmm3	/* t1 */			__asm	movaps	[esi+0x50],xmm6	/* t5 */

		__asm	movaps	xmm0,[esi+0x20]	/* t2 */			__asm	movaps	xmm4,[esi+0x60]	/* t6 */
		__asm	movaps	xmm1,[esi+0x30]	/* t3 */			__asm	movaps	xmm5,[esi+0x70]	/* t7 */
		__asm	movaps	xmm2,[esi+0xa0]	/* cpy ta */		__asm	movaps	xmm7,[esi+0xf0]	/* tf */
		__asm	movaps	xmm3,[esi+0xb0]	/* cpy tb */		__asm	movaps	xmm6,[esi+0xe0]	/* te */
		__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* t6 = 6-f */
		__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* tf = 7-e */
		__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */
		__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */
		__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* te = 6+f */
		__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* t7 = 7+e */

		__asm	movaps	[esi+0xa0],xmm0	/* ta */			__asm	movaps	[esi+0x60],xmm4	/* t6 */
		__asm	movaps	[esi+0xb0],xmm1	/* tb */			__asm	movaps	[esi+0xf0],xmm5	/* tf */
		__asm	movaps	[esi+0x20],xmm2	/* t2 */			__asm	movaps	[esi+0xe0],xmm7	/* te */
		__asm	movaps	[esi+0x30],xmm3	/* t3 */			__asm	movaps	[esi+0x70],xmm6	/* t7 */

	  #endif	/* if(1) */

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
		__asm	mov	eax, r00
		__asm	mov	ebx, r10
		__asm	mov	ecx, r20
		__asm	mov	edx, r30

		__asm	movaps	xmm0,[eax      ]	/* t00 */
		__asm	movaps	xmm1,[eax+0x010]	/* t01 */
		__asm	movaps	xmm2,[ebx      ]	/* t10 */
		__asm	movaps	xmm3,[ebx+0x010]	/* t11 */

		__asm	subpd	xmm0,[ebx      ]	/* t10=t00-rt */
		__asm	subpd	xmm1,[ebx+0x010]	/* t11=t01-it */
		__asm	addpd	xmm2,[eax      ]	/* t00=t00+rt */
		__asm	addpd	xmm3,[eax+0x010]	/* t01=t01+it */

		__asm	movaps	xmm4,[ecx      ]	/* t20 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t21 */
		__asm	movaps	xmm6,[edx      ]	/* t30 */
		__asm	movaps	xmm7,[edx+0x010]	/* t31 */

		__asm	subpd	xmm4,[edx      ]	/* t30=t20-rt */
		__asm	subpd	xmm5,[edx+0x010]	/* t31=t21-it */
		__asm	addpd	xmm6,[ecx      ]	/* t20=t20+rt */
		__asm	addpd	xmm7,[ecx+0x010]	/* t21=t21+it */

		__asm	subpd	xmm2,xmm6		/* t00 <- t00-t20 */
		__asm	subpd	xmm3,xmm7		/* t01 <- t01-t21 */
		__asm	addpd	xmm6,xmm6		/*          2*t20 */
		__asm	addpd	xmm7,xmm7		/*          2*t21 */
		__asm	movaps	[ecx     ],xmm2	/* a1p01r, store in r20 */
		__asm	movaps	[ecx+0x10],xmm3	/* a1p01i, store in r21 */
		__asm	addpd	xmm6,xmm2		/* t20 <- t00+t20 */
		__asm	addpd	xmm7,xmm3		/* t21 <- t01+t21 */
		__asm	movaps	[eax     ],xmm6	/* a1p00r, store in r00 */
		__asm	movaps	[eax+0x10],xmm7	/* a1p00i, store in r01 */

		__asm	subpd	xmm0,xmm5		/* t10 <- t10-t31 */
		__asm	subpd	xmm1,xmm4		/* t11 <- t11-t30 */
		__asm	addpd	xmm5,xmm5		/*          2*t31 */
		__asm	addpd	xmm4,xmm4		/*          2*t30 */
		__asm	movaps	[ebx     ],xmm0	/* a1p02r, store in r10 */
		__asm	movaps	[edx+0x10],xmm1	/* a1p03i, store in r31 */
		__asm	addpd	xmm5,xmm0		/* t31 <- t10+t31 */
		__asm	addpd	xmm4,xmm1		/* t30 <- t11+t30 */
		__asm	movaps	[edx     ],xmm5	/* a1p03r, store in r30 */
		__asm	movaps	[ebx+0x10],xmm4	/* a1p03i, store in r11 */

	/*...Block 5: t08,t18,t28,t38	*/
		__asm	add	eax, 0x080	/* r08 */
		__asm	add	ebx, 0x080	/* r18 */
		__asm	add	ecx, 0x080	/* r28 */
		__asm	add	edx, 0x080	/* r38 */
		__asm	mov	esi, isrt2

		__asm	movaps	xmm0,[eax      ]	/* t08 */
		__asm	movaps	xmm1,[eax+0x010]	/* t09 */
		__asm	movaps	xmm2,[ebx      ]	/* t18 */
		__asm	movaps	xmm3,[ebx+0x010]	/* t19 */

		__asm	subpd	xmm0,[ebx+0x010]	/* t08=t08-t19*/
		__asm	subpd	xmm1,[ebx      ]	/* t19=t09-t18*/
		__asm	addpd	xmm2,[eax+0x010]	/* t09=t18+t09*/
		__asm	addpd	xmm3,[eax      ]	/* t18=t19+t08*/

		__asm	movaps	xmm4,[ecx      ]	/* t28 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t29 */
		__asm	movaps	xmm6,[edx      ]	/* t38 */
		__asm	movaps	xmm7,[edx+0x010]	/* t39 */

		__asm	subpd	xmm4,[ecx+0x010]	/* t28-t29 */
		__asm	addpd	xmm5,[ecx      ]	/* t29+t28 */
		__asm	mulpd	xmm4,[esi]	/* t28 = (t28-t29)*ISRT2 */
		__asm	mulpd	xmm5,[esi]	/* t29 = (t29+t28)*ISRT2 */

		__asm	addpd	xmm6,[edx+0x010]	/* t38+t39 */
		__asm	subpd	xmm7,[edx      ]	/* t39-t38 */
		__asm	mulpd	xmm6,[esi]	/*  rt = (t38+t39)*ISRT2 */
		__asm	mulpd	xmm7,[esi]	/*  it = (t39-t38)*ISRT2 */

		__asm	subpd	xmm4,xmm6		/* t28=t28-rt */
		__asm	subpd	xmm5,xmm7		/* t29=t29-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/* t38=t28+rt */
		__asm	addpd	xmm7,xmm5		/* t39=t29+it */

		__asm	subpd	xmm0,xmm4		/* t08-t28 */
		__asm	subpd	xmm2,xmm5		/* t09-t29 */
		__asm	addpd	xmm4,xmm4		/*   2*t28 */
		__asm	addpd	xmm5,xmm5		/*   2*t29 */

		__asm	movaps	[ecx     ],xmm0	/* a1p05r, store in r28 */
		__asm	movaps	[ecx+0x10],xmm2	/* a1p05i, store in r29 */
		__asm	addpd	xmm4,xmm0		/* t08+t28 */
		__asm	addpd	xmm5,xmm2		/* t09+t29 */
		__asm	movaps	[eax     ],xmm4	/* a1p04r, store in r08 */
		__asm	movaps	[eax+0x10],xmm5	/* a1p04i, store in r09 */

		__asm	subpd	xmm3,xmm7		/* t18-t39 */
		__asm	subpd	xmm1,xmm6		/* t19-t38 */
		__asm	addpd	xmm7,xmm7		/*   2*t39 */
		__asm	addpd	xmm6,xmm6		/*   2*t38 */
		__asm	movaps	[ebx     ],xmm3	/* a1p06r, store in r18 */
		__asm	movaps	[edx+0x10],xmm1	/* a1p07i, store in r39 */
		__asm	addpd	xmm7,xmm3		/* t18+t39 */
		__asm	addpd	xmm6,xmm1		/* t19+t38 */
		__asm	movaps	[edx     ],xmm7	/* a1p07r, store in r38 */
		__asm	movaps	[ebx+0x10],xmm6	/* a1p06i, store in r19 */

	/*...Block 3: t04,t14,t24,t34	*/
		__asm	sub	eax, 0x040	/* r04 */
		__asm	sub	ebx, 0x040	/* r14 */
		__asm	sub	ecx, 0x040	/* r24 */
		__asm	sub	edx, 0x040	/* r34 */
		__asm	mov	edi, cc0

		__asm	movaps	xmm4,[ecx      ]	/* t24 */		__asm	movaps	xmm6,[edx      ]	/* t34 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t25 */		__asm	movaps	xmm7,[edx+0x010]	/* t35 */
		__asm	movaps	xmm0,[ecx      ]	/* copy t24 */	__asm	movaps	xmm2,[edx      ]	/* copy t34 */
		__asm	movaps	xmm1,[ecx+0x010]	/* copy t25 */	__asm	movaps	xmm3,[edx+0x010]	/* copy t35 */

		__asm	mulpd	xmm4,[edi     ]	/* t24*c */			__asm	mulpd	xmm6,[edi+0x10]	/* t34*s */
		__asm	mulpd	xmm1,[edi+0x10]	/* t25*s */			__asm	mulpd	xmm3,[edi     ]	/* t35*c */
		__asm	mulpd	xmm5,[edi     ]	/* t25*c */			__asm	mulpd	xmm7,[edi+0x10]	/* t35*s */
		__asm	mulpd	xmm0,[edi+0x10]	/* t24*s */			__asm	mulpd	xmm2,[edi     ]	/* t34*c */
		__asm	subpd	xmm4,xmm1	/* ~t24 */				__asm	subpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t25 */				__asm	addpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t34=t24-rt */
		__asm	subpd	xmm5,xmm7		/*~t35=t25-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t24=t24+rt */
		__asm	addpd	xmm7,xmm5		/*~t25=t25+it */

		__asm	movaps	xmm2,[ebx      ]	/* t14 */
		__asm	movaps	xmm3,[ebx+0x010]	/* t15 */
		__asm	subpd	xmm2,[ebx+0x010]	/* t14-t15 */
		__asm	addpd	xmm3,[ebx      ]	/* t15+t14 */
		__asm	mulpd	xmm2,[esi]	/* rt = (t14-t15)*ISRT2 */
		__asm	mulpd	xmm3,[esi]	/* it = (t15+t14)*ISRT2 */

		__asm	movaps	xmm0,[eax      ]	/* t04 */
		__asm	movaps	xmm1,[eax+0x010]	/* t05 */

		__asm	subpd	xmm0,xmm2			/*~t14=t04-rt */
		__asm	subpd	xmm1,xmm3			/*~t15=t05-it */
		__asm	addpd	xmm2,[eax      ]	/*~t04=rt +t04*/
		__asm	addpd	xmm3,[eax+0x010]	/*~t05=it +t05*/

		__asm	subpd	xmm2,xmm6		/* t04-t24 */
		__asm	subpd	xmm3,xmm7		/* t05-t25 */
		__asm	addpd	xmm6,xmm6		/*   2*t24 */
		__asm	addpd	xmm7,xmm7		/*   2*t25 */
		__asm	movaps	[ecx     ],xmm2	/* a1p09r, store in r24 */
		__asm	movaps	[ecx+0x10],xmm3	/* a1p09i, store in r25 */
		__asm	addpd	xmm6,xmm2		/* t04+t24 */
		__asm	addpd	xmm7,xmm3		/* t05+t25 */
		__asm	movaps	[eax     ],xmm6	/* a1p08r, store in r04 */
		__asm	movaps	[eax+0x10],xmm7	/* a1p08i, store in r05 */

		__asm	subpd	xmm0,xmm5		/* t14-t35 */
		__asm	subpd	xmm1,xmm4		/* t15-t34 */
		__asm	addpd	xmm5,xmm5		/*          2*t35 */
		__asm	addpd	xmm4,xmm4		/*          2*t34 */
		__asm	movaps	[ebx     ],xmm0	/* a1p0Ar, store in r14 */
		__asm	movaps	[edx+0x10],xmm1	/* a1p0Bi, store in r35 */
		__asm	addpd	xmm5,xmm0		/* t14+t35 */
		__asm	addpd	xmm4,xmm1		/* t15+t34 */
		__asm	movaps	[edx     ],xmm5	/* a1p0Br, store in r34 */
		__asm	movaps	[ebx+0x10],xmm4	/* a1p0Ai, store in r15 */

	/*...Block 7: t0C,t1C,t2C,t3C	*/
		__asm	add	eax, 0x080	/* r0C */
		__asm	add	ebx, 0x080	/* r1C */
		__asm	add	ecx, 0x080	/* r2C */
		__asm	add	edx, 0x080	/* r3C */

		__asm	movaps	xmm4,[ecx      ]	/* t2C */		__asm	movaps	xmm6,[edx      ]	/* t3C */
		__asm	movaps	xmm5,[ecx+0x010]	/* t2D */		__asm	movaps	xmm7,[edx+0x010]	/* t3D */
		__asm	movaps	xmm0,[ecx      ]	/* copy t2C */	__asm	movaps	xmm2,[edx      ]	/* copy t3C */
		__asm	movaps	xmm1,[ecx+0x010]	/* copy t2D */	__asm	movaps	xmm3,[edx+0x010]	/* copy t3D */

		__asm	mulpd	xmm4,[edi+0x10]	/* t2C*s */			__asm	mulpd	xmm6,[edi     ]	/* t3C*c */
		__asm	mulpd	xmm1,[edi     ]	/* t2D*c */			__asm	mulpd	xmm3,[edi+0x10]	/* t3D*s */
		__asm	mulpd	xmm5,[edi+0x10]	/* t2D*s */			__asm	mulpd	xmm7,[edi     ]	/* t3D*c */
		__asm	mulpd	xmm0,[edi     ]	/* t2C*c */			__asm	mulpd	xmm2,[edi+0x10]	/* t3C*s */
		__asm	subpd	xmm4,xmm1	/* ~t24 */				__asm	subpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t25 */				__asm	addpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t2C=t2C-rt */
		__asm	subpd	xmm5,xmm7		/*~t2D=t2D-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t3C=t2C+rt */
		__asm	addpd	xmm7,xmm5		/*~t3D=t2D+it */

		__asm	movaps	xmm2,[ebx      ]	/* t1C */
		__asm	movaps	xmm3,[ebx+0x010]	/* t1D */
		__asm	addpd	xmm2,[ebx+0x010]	/* t1C+t1D */
		__asm	subpd	xmm3,[ebx      ]	/* t1D-t1C */
		__asm	mulpd	xmm2,[esi]	/* rt = (t1C+t1D)*ISRT2 */
		__asm	mulpd	xmm3,[esi]	/* it = (t1D-t1C)*ISRT2 */

		__asm	movaps	xmm0,[eax      ]	/* t0C */
		__asm	movaps	xmm1,[eax+0x010]	/* t0D */

		__asm	subpd	xmm0,xmm2			/*~t0C=t0C-rt */
		__asm	subpd	xmm1,xmm3			/*~t0D=t0D-it */
		__asm	addpd	xmm2,[eax      ]	/*~t1C=rt +t0C*/
		__asm	addpd	xmm3,[eax+0x010]	/*~t1D=it +t0D*/

		__asm	subpd	xmm0,xmm4		/* t0C-t2C */
		__asm	subpd	xmm1,xmm5		/* t0D-t2D */
		__asm	addpd	xmm4,xmm4		/*   2*t2C */
		__asm	addpd	xmm5,xmm5		/*   2*t2D */
		__asm	movaps	[ecx     ],xmm0	/* a1p0Dr, store in r2C */
		__asm	movaps	[ecx+0x10],xmm1	/* a1p0Di, store in r2D */
		__asm	addpd	xmm4,xmm0		/* t0C+t2C */
		__asm	addpd	xmm5,xmm1		/* t0D+t2D */
		__asm	movaps	[eax     ],xmm4	/* a1p0Cr, store in r0C */
		__asm	movaps	[eax+0x10],xmm5	/* a1p0Ci, store in r0D */

		__asm	subpd	xmm2,xmm7		/* t1C-t3D */
		__asm	subpd	xmm3,xmm6		/* t1D-t3C */
		__asm	addpd	xmm7,xmm7		/*   2*t3D */
		__asm	addpd	xmm6,xmm6		/*   2*t3C */
		__asm	movaps	[ebx     ],xmm2	/* a1p0Er, store in r1C */
		__asm	movaps	[edx+0x10],xmm3	/* a1p0Fi, store in r3D */
		__asm	addpd	xmm7,xmm2		/* t1C+t3D */
		__asm	addpd	xmm6,xmm3		/* t1D+t3C */
		__asm	movaps	[edx     ],xmm7	/* a1p0Fr, store in r3C */
		__asm	movaps	[ebx+0x10],xmm6	/* a1p0Ei, store in r1D */

	/*...Block 2: t02,t12,t22,t32	*/
		__asm	sub	eax, 0x0a0	/* r02 */
		__asm	sub	ebx, 0x0a0	/* r12 */
		__asm	sub	ecx, 0x0a0	/* r22 */
		__asm	sub	edx, 0x0a0	/* r32 */
		__asm	add esi, 0x30	/* cc1 */

		__asm	movaps	xmm4,[ecx      ]	/* t22 */		__asm	movaps	xmm6,[edx      ]	/* t32 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t23 */		__asm	movaps	xmm7,[edx+0x010]	/* t33 */
		__asm	movaps	xmm0,[ecx      ]	/* copy t22 */	__asm	movaps	xmm2,[edx      ]	/* copy t32 */
		__asm	movaps	xmm1,[ecx+0x010]	/* copy t23 */	__asm	movaps	xmm3,[edx+0x010]	/* copy t33 */

		__asm	mulpd	xmm4,[esi     ]	/* t22*c32_1 */		__asm	mulpd	xmm6,[esi+0x20]	/* t32*c32_3 */
		__asm	mulpd	xmm1,[esi+0x10]	/* t23*s32_1 */		__asm	mulpd	xmm3,[esi+0x30]	/* t33*s32_3 */
		__asm	mulpd	xmm5,[esi     ]	/* t23*c32_1 */		__asm	mulpd	xmm7,[esi+0x20]	/* t33*c32_3 */
		__asm	mulpd	xmm0,[esi+0x10]	/* t22*s32_1 */		__asm	mulpd	xmm2,[esi+0x30]	/* t32*s32_3 */
		__asm	subpd	xmm4,xmm1	/* ~t22 */				__asm	subpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t23 */				__asm	addpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t32=t22-rt */
		__asm	subpd	xmm5,xmm7		/*~t33=t23-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t22=t22+rt */
		__asm	addpd	xmm7,xmm5		/*~t23=t23+it */

		__asm	movaps	xmm2,[ebx      ]	/* t12 */
		__asm	movaps	xmm0,[ebx+0x010]	/* t13 */
		__asm	movaps	xmm1,[ebx      ]	/* copy t12 */
		__asm	movaps	xmm3,[ebx+0x010]	/* copy t13 */

		__asm	mulpd	xmm2,[edi     ]	/* t12*c */
		__asm	mulpd	xmm0,[edi+0x10]	/* t13*s */
		__asm	mulpd	xmm3,[edi     ]	/* t13*c */
		__asm	mulpd	xmm1,[edi+0x10]	/* t12*s */
		__asm	subpd	xmm2,xmm0	/* rt */
		__asm	addpd	xmm3,xmm1	/* it */

		__asm	movaps	xmm0,[eax      ]	/* t02 */
		__asm	movaps	xmm1,[eax+0x010]	/* t03 */
		__asm	subpd	xmm0,xmm2		/*~t12=t02-rt */
		__asm	subpd	xmm1,xmm3		/*~t13=t03-it */
		__asm	addpd	xmm2,[eax      ]/*~t02=rt+t02 */
		__asm	addpd	xmm3,[eax+0x010]/*~t03=it+t03 */

		__asm	subpd	xmm2,xmm6		/* t02-t22 */
		__asm	subpd	xmm3,xmm7		/* t03-t23 */
		__asm	addpd	xmm6,xmm6		/*   2*t22 */
		__asm	addpd	xmm7,xmm7		/*   2*t23 */
		__asm	movaps	[ecx     ],xmm2	/* a1p11r, store in r22 */
		__asm	movaps	[ecx+0x10],xmm3	/* a1p11i, store in r23 */
		__asm	addpd	xmm6,xmm2		/* t02+t22 */
		__asm	addpd	xmm7,xmm3		/* t03+t23 */
		__asm	movaps	[eax     ],xmm6	/* a1p10r, store in r02 */
		__asm	movaps	[eax+0x10],xmm7	/* a1p10i, store in r03 */

		__asm	subpd	xmm0,xmm5		/* t12-t33 */
		__asm	subpd	xmm1,xmm4		/* t13-t32 */
		__asm	addpd	xmm5,xmm5		/*   2*t33 */
		__asm	addpd	xmm4,xmm4		/*   2*t32 */
		__asm	movaps	[ebx     ],xmm0	/* a1p12r, store in r12 */
		__asm	movaps	[edx+0x10],xmm1	/* a1p13i, store in r33 */
		__asm	addpd	xmm5,xmm0		/* t12+t33 */
		__asm	addpd	xmm4,xmm1		/* t13+t32 */
		__asm	movaps	[edx     ],xmm5	/* a1p13r, store in r32 */
		__asm	movaps	[ebx+0x10],xmm4	/* a1p12i, store in r13 */

	/*...Block 6: t0A,t1A,t2A,t3A	*/
		__asm	add	eax, 0x080	/* r0A */
		__asm	add	ebx, 0x080	/* r1A */
		__asm	add	ecx, 0x080	/* r2A */
		__asm	add	edx, 0x080	/* r3A */

		__asm	movaps	xmm4,[ecx      ]	/* t2A */		__asm	movaps	xmm6,[edx      ]	/* t3A */
		__asm	movaps	xmm5,[ecx+0x010]	/* t2B */		__asm	movaps	xmm7,[edx+0x010]	/* t3B */
		__asm	movaps	xmm0,[ecx      ]	/* copy t2A */	__asm	movaps	xmm2,[edx      ]	/* copy t3A */
		__asm	movaps	xmm1,[ecx+0x010]	/* copy t2B */	__asm	movaps	xmm3,[edx+0x010]	/* copy t3B */

		__asm	mulpd	xmm4,[esi+0x30]	/* t2A*s32_3 */		__asm	mulpd	xmm6,[esi     ]	/* t3A*c32_1 */
		__asm	mulpd	xmm1,[esi+0x20]	/* t2B*c32_3 */		__asm	mulpd	xmm3,[esi+0x10]	/* t3B*s32_1 */
		__asm	mulpd	xmm5,[esi+0x30]	/* t2B*s32_3 */		__asm	mulpd	xmm7,[esi     ]	/* t3B*c32_1 */
		__asm	mulpd	xmm0,[esi+0x20]	/* t2A*c32_3 */		__asm	mulpd	xmm2,[esi+0x10]	/* t3A*s32_1 */
		__asm	subpd	xmm4,xmm1	/* ~t2A */				__asm	addpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t2B */				__asm	subpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t2A=t2A-rt */
		__asm	subpd	xmm5,xmm7		/*~t2B=t2B-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t3A=t2A+rt */
		__asm	addpd	xmm7,xmm5		/*~t3B=t2B+it */

		__asm	movaps	xmm2,[ebx      ]	/* t1A */
		__asm	movaps	xmm0,[ebx+0x010]	/* t1B */
		__asm	movaps	xmm1,[ebx      ]	/* copy t1A */
		__asm	movaps	xmm3,[ebx+0x010]	/* copy t1B */

		__asm	mulpd	xmm2,[edi+0x10]	/* t1A*s */
		__asm	mulpd	xmm0,[edi     ]	/* t1B*c */
		__asm	mulpd	xmm3,[edi+0x10]	/* t1B*s */
		__asm	mulpd	xmm1,[edi     ]	/* t1A*c */
		__asm	addpd	xmm2,xmm0	/* rt */
		__asm	subpd	xmm3,xmm1	/* it */

		__asm	movaps	xmm0,[eax      ]	/* t0A */
		__asm	movaps	xmm1,[eax+0x010]	/* t0B */
		__asm	subpd	xmm0,xmm2		/*~t0A=t0A-rt */
		__asm	subpd	xmm1,xmm3		/*~t0B=t0B-it */
		__asm	addpd	xmm2,[eax      ]/*~t1A=rt+t0A */
		__asm	addpd	xmm3,[eax+0x010]/*~t1B=it+t0B */

		__asm	subpd	xmm0,xmm4		/* t0A-t2A */
		__asm	subpd	xmm1,xmm5		/* t0B-t2B */
		__asm	addpd	xmm4,xmm4		/*   2*t2A */
		__asm	addpd	xmm5,xmm5		/*   2*t2B */
		__asm	movaps	[ecx     ],xmm0	/* a1p15r, store in r2A */
		__asm	movaps	[ecx+0x10],xmm1	/* a1p15i, store in r2B */
		__asm	addpd	xmm4,xmm0		/* t0A+t2A */
		__asm	addpd	xmm5,xmm1		/* t0B+t2B */
		__asm	movaps	[eax     ],xmm4	/* a1p14r, store in r0A */
		__asm	movaps	[eax+0x10],xmm5	/* a1p14i, store in r0B */

		__asm	subpd	xmm2,xmm7		/* t1A-t3B */
		__asm	subpd	xmm3,xmm6		/* t1B-t3A */
		__asm	addpd	xmm7,xmm7		/*   2*t3B */
		__asm	addpd	xmm6,xmm6		/*   2*t3A */
		__asm	movaps	[ebx     ],xmm2	/* a1p16r, store in r1A */
		__asm	movaps	[edx+0x10],xmm3	/* a1p17i, store in r3B */
		__asm	addpd	xmm7,xmm2		/* t1A+t3B */
		__asm	addpd	xmm6,xmm3		/* t1B+t3A */
		__asm	movaps	[edx     ],xmm7	/* a1p17r, store in r3A */
		__asm	movaps	[ebx+0x10],xmm6	/* a1p16i, store in r1B */

	/*...Block 4: t06,t16,t26,t36	*/
		__asm	sub	eax, 0x040	/* r06 */
		__asm	sub	ebx, 0x040	/* r16 */
		__asm	sub	ecx, 0x040	/* r26 */
		__asm	sub	edx, 0x040	/* r36 */

		__asm	movaps	xmm4,[ecx      ]	/* t26 */		__asm	movaps	xmm6,[edx      ]	/* t36 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t27 */		__asm	movaps	xmm7,[edx+0x010]	/* t37 */
		__asm	movaps	xmm0,[ecx      ]	/* copy t26 */	__asm	movaps	xmm2,[edx      ]	/* copy t36 */
		__asm	movaps	xmm1,[ecx+0x010]	/* copy t27 */	__asm	movaps	xmm3,[edx+0x010]	/* copy t37 */

		__asm	mulpd	xmm4,[esi+0x20]	/* t26*s32_3 */		__asm	mulpd	xmm6,[esi+0x10]	/* t36*s32_1 */
		__asm	mulpd	xmm1,[esi+0x30]	/* t27*s32_3 */		__asm	mulpd	xmm3,[esi     ]	/* t37*c32_1 */
		__asm	mulpd	xmm5,[esi+0x20]	/* t27*c32_3 */		__asm	mulpd	xmm7,[esi+0x10]	/* t37*s32_1 */
		__asm	mulpd	xmm0,[esi+0x30]	/* t26*s32_3 */		__asm	mulpd	xmm2,[esi     ]	/* t36*c32_1 */
		__asm	subpd	xmm4,xmm1	/* ~t26 */				__asm	addpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t27 */				__asm	subpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t26=t26-rt */
		__asm	subpd	xmm5,xmm7		/*~t27=t27-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t36=t26+rt */
		__asm	addpd	xmm7,xmm5		/*~t37=t27+it */

		__asm	movaps	xmm2,[ebx      ]	/* t16 */
		__asm	movaps	xmm0,[ebx+0x010]	/* t17 */
		__asm	movaps	xmm1,[ebx      ]	/* copy t16 */
		__asm	movaps	xmm3,[ebx+0x010]	/* copy t17 */

		__asm	mulpd	xmm2,[edi+0x10]	/* t16*s */
		__asm	mulpd	xmm0,[edi     ]	/* t17*c */
		__asm	mulpd	xmm3,[edi+0x10]	/* t17*s */
		__asm	mulpd	xmm1,[edi     ]	/* t16*c */
		__asm	subpd	xmm2,xmm0	/* rt */
		__asm	addpd	xmm3,xmm1	/* it */

		__asm	movaps	xmm0,[eax      ]	/* t06 */
		__asm	movaps	xmm1,[eax+0x010]	/* t07 */
		__asm	subpd	xmm0,xmm2		/*~t16=t06-rt */
		__asm	subpd	xmm1,xmm3		/*~t17=t07-it */
		__asm	addpd	xmm2,[eax      ]/*~t06=rt+t06 */
		__asm	addpd	xmm3,[eax+0x010]/*~t07=it+t07 */

		__asm	subpd	xmm2,xmm4		/* t06-t26 */
		__asm	subpd	xmm3,xmm5		/* t07-t27 */
		__asm	addpd	xmm4,xmm4		/*   2*t26 */
		__asm	addpd	xmm5,xmm5		/*   2*t27 */
		__asm	movaps	[ecx     ],xmm2	/* a1p19r, store in r26 */
		__asm	movaps	[ecx+0x10],xmm3	/* a1p19i, store in r27 */
		__asm	addpd	xmm4,xmm2		/* t06+t26 */
		__asm	addpd	xmm5,xmm3		/* t07+t27 */
		__asm	movaps	[eax     ],xmm4	/* a1p18r, store in r06 */
		__asm	movaps	[eax+0x10],xmm5	/* a1p18i, store in r07 */

		__asm	subpd	xmm0,xmm7		/* t16-t37 */
		__asm	subpd	xmm1,xmm6		/* t17-t36 */
		__asm	addpd	xmm7,xmm7		/*   2*t37 */
		__asm	addpd	xmm6,xmm6		/*   2*t36 */
		__asm	movaps	[ebx     ],xmm0	/* a1p1Ar, store in r16 */
		__asm	movaps	[edx+0x10],xmm1	/* a1p1Bi, store in r37 */
		__asm	addpd	xmm7,xmm0		/* t16+t37 */
		__asm	addpd	xmm6,xmm1		/* t17+t36 */
		__asm	movaps	[edx     ],xmm7	/* a1p1Br, store in r36 */
		__asm	movaps	[ebx+0x10],xmm6	/* a1p1Ai, store in r17 */

	/*...Block 8: t0E,t1E,t2E,t3E	*/
		__asm	add	eax, 0x080	/* r0E */
		__asm	add	ebx, 0x080	/* r1E */
		__asm	add	ecx, 0x080	/* r2E */
		__asm	add	edx, 0x080	/* r3E */

		__asm	movaps	xmm4,[ecx      ]	/* t2E */		__asm	movaps	xmm6,[edx      ]	/* t3E */
		__asm	movaps	xmm5,[ecx+0x010]	/* t2F */		__asm	movaps	xmm7,[edx+0x010]	/* t3F */
		__asm	movaps	xmm0,[ecx      ]	/* copy t2E */	__asm	movaps	xmm2,[edx      ]	/* copy t3E */
		__asm	movaps	xmm1,[ecx+0x010]	/* copy t2F */	__asm	movaps	xmm3,[edx+0x010]	/* copy t3F */

		__asm	mulpd	xmm4,[esi+0x10]	/* t2E*s32_1 */		__asm	mulpd	xmm6,[esi+0x30]	/* t3E*c32_3 */
		__asm	mulpd	xmm1,[esi     ]	/* t2F*c32_1 */		__asm	mulpd	xmm3,[esi+0x20]	/* t3F*s32_3 */
		__asm	mulpd	xmm5,[esi+0x10]	/* t2F*s32_1 */		__asm	mulpd	xmm7,[esi+0x30]	/* t3F*c32_3 */
		__asm	mulpd	xmm0,[esi     ]	/* t2E*c32_1 */		__asm	mulpd	xmm2,[esi+0x20]	/* t3E*s32_3 */
		__asm	subpd	xmm4,xmm1	/* ~t2E */				__asm	subpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t2F */				__asm	addpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t2E=t2E-rt */
		__asm	subpd	xmm5,xmm7		/*~t2F=t2F-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t3E=t2E+rt */
		__asm	addpd	xmm7,xmm5		/*~t3F=t2F+it */

		__asm	movaps	xmm2,[ebx      ]	/* t1E */
		__asm	movaps	xmm0,[ebx+0x010]	/* t1F */
		__asm	movaps	xmm1,[ebx      ]	/* copy t1E */
		__asm	movaps	xmm3,[ebx+0x010]	/* copy t1F */

		__asm	mulpd	xmm2,[edi     ]	/* t1E*c */
		__asm	mulpd	xmm0,[edi+0x10]	/* t1F*s */
		__asm	mulpd	xmm3,[edi     ]	/* t1F*c */
		__asm	mulpd	xmm1,[edi+0x10]	/* t1E*s */
		__asm	addpd	xmm2,xmm0	/* rt */
		__asm	subpd	xmm3,xmm1	/* it */

		__asm	movaps	xmm0,[eax      ]	/* t0E */
		__asm	movaps	xmm1,[eax+0x010]	/* t0F */
		__asm	subpd	xmm0,xmm2		/*~t0E=t0E-rt */
		__asm	subpd	xmm1,xmm3		/*~t0F=t0F-it */
		__asm	addpd	xmm2,[eax      ]/*~t1E=rt+t0E */
		__asm	addpd	xmm3,[eax+0x010]/*~t1F=it+t0F */

		__asm	subpd	xmm0,xmm4		/* t0E-t2E */
		__asm	subpd	xmm1,xmm5		/* t0F-t2F */
		__asm	addpd	xmm4,xmm4		/*   2*t2E */
		__asm	addpd	xmm5,xmm5		/*   2*t2F */
		__asm	movaps	[ecx     ],xmm0	/* a1p1Dr, store in r2E */
		__asm	movaps	[ecx+0x10],xmm1	/* a1p1Di, store in r2F */
		__asm	addpd	xmm4,xmm0		/* t0E+t2E */
		__asm	addpd	xmm5,xmm1		/* t0F+t2F */
		__asm	movaps	[eax     ],xmm4	/* a1p1Cr, store in r0E */
		__asm	movaps	[eax+0x10],xmm5	/* a1p1Ci, store in r0F */

		__asm	subpd	xmm2,xmm7		/* t1E-t3F */
		__asm	subpd	xmm3,xmm6		/* t1F-t3E */
		__asm	addpd	xmm7,xmm7		/*   2*t3F */
		__asm	addpd	xmm6,xmm6		/*   2*t3E */
		__asm	movaps	[ebx     ],xmm2	/* a1pE2r, store in r1E */
		__asm	movaps	[edx+0x10],xmm3	/* a1pF2i, store in r3F */
		__asm	addpd	xmm7,xmm2		/* t1E+t3F */
		__asm	addpd	xmm6,xmm3		/* t1F+t3E */
		__asm	movaps	[edx     ],xmm7	/* a1pF3r, store in r3E */
		__asm	movaps	[ebx+0x10],xmm6	/* a1pE3i, store in r1F */

	#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

	  #ifdef USE_AVX	// process 4 main-array blocks of 8 vec_dbl = 8 x 4 = 32 doubles each in AVX mode:

		SSE2_RADIX32_WRAPPER_DIF(add0,add1,add2,add3,r00,r10,r20,r30,isrt2,cc0,c00,c01,c02,c03,c05,c07)

	  #else	// SSE2:

		SSE2_RADIX32_WRAPPER_DIF(add0,add1          ,r00,r10,r20,r30,isrt2,cc0,c00,c01,c02,c03,c05,c07)

	  #endif

	#endif

/*
!...send the pairs of complex elements which are to be combined and sincos temporaries needed for the squaring to a
!   small subroutine. The j1 = 0 case is again exceptional. [For this reason we don't supply SSE2 code it - not worth the work].
*/
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

		re = c01->d2;		im = (c01+1)->d2;
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

	#if defined(COMPILER_TYPE_MSVC)

		__asm	mov eax, tmp0
		__asm	movaps	xmm0,[eax     ]
		__asm	movaps	xmm1,[eax+0x10]
		__asm	movaps	xmm2,[eax+0x20]
		__asm	movaps	xmm3,[eax+0x30]
		__asm	movaps	xmm4,[eax+0x40]
		__asm	movaps	xmm5,[eax+0x50]
		__asm	movaps	xmm6,[eax+0x60]
		__asm	movaps	xmm7,[eax+0x70]
		__asm	shufpd	xmm0,xmm0,1
		__asm	shufpd	xmm1,xmm1,1
		__asm	shufpd	xmm2,xmm2,1
		__asm	shufpd	xmm3,xmm3,1
		__asm	shufpd	xmm4,xmm4,1
		__asm	shufpd	xmm5,xmm5,1
		__asm	shufpd	xmm6,xmm6,1
		__asm	shufpd	xmm7,xmm7,1
		__asm	movaps	[eax     ],xmm0
		__asm	movaps	[eax+0x10],xmm1
		__asm	movaps	[eax+0x20],xmm2
		__asm	movaps	[eax+0x30],xmm3
		__asm	movaps	[eax+0x40],xmm4
		__asm	movaps	[eax+0x50],xmm5
		__asm	movaps	[eax+0x60],xmm6
		__asm	movaps	[eax+0x70],xmm7

	#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		#if OS_BITS == 32

		__asm__ volatile (\
			"movl	%[tmp0],%%eax\n\t"\
			"movaps	    (%%eax),%%xmm0\n\t"\
			"movaps	0x10(%%eax),%%xmm1\n\t"\
			"movaps	0x20(%%eax),%%xmm2\n\t"\
			"movaps	0x30(%%eax),%%xmm3\n\t"\
			"movaps	0x40(%%eax),%%xmm4\n\t"\
			"movaps	0x50(%%eax),%%xmm5\n\t"\
			"movaps	0x60(%%eax),%%xmm6\n\t"\
			"movaps	0x70(%%eax),%%xmm7\n\t"\
			"shufpd	$1	,%%xmm0	,%%xmm0\n\t"\
			"shufpd	$1	,%%xmm1	,%%xmm1\n\t"\
			"shufpd	$1	,%%xmm2	,%%xmm2\n\t"\
			"shufpd	$1	,%%xmm3	,%%xmm3\n\t"\
			"shufpd	$1	,%%xmm4	,%%xmm4\n\t"\
			"shufpd	$1	,%%xmm5	,%%xmm5\n\t"\
			"shufpd	$1	,%%xmm6	,%%xmm6\n\t"\
			"shufpd	$1	,%%xmm7	,%%xmm7\n\t"\
			"movaps	%%xmm0,    (%%eax)\n\t"\
			"movaps	%%xmm1,0x10(%%eax)\n\t"\
			"movaps	%%xmm2,0x20(%%eax)\n\t"\
			"movaps	%%xmm3,0x30(%%eax)\n\t"\
			"movaps	%%xmm4,0x40(%%eax)\n\t"\
			"movaps	%%xmm5,0x50(%%eax)\n\t"\
			"movaps	%%xmm6,0x60(%%eax)\n\t"\
			"movaps	%%xmm7,0x70(%%eax)\n\t"\
			:					// outputs: none
			: [tmp0] "m" (tmp0)	// All inputs from memory addresses here
			: "cc","memory","eax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	// Clobbered registers
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

		#else					// 64-bit SSE2 build

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

		#endif

	#endif

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

	/************************************************************************/
	/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/
	/************************************************************************/

	/* Totals: 84 load/store [66 movaps], 68 add/subpd,  4 mulpd, 48 address-compute per radix-8 block here: */\

	#if defined(COMPILER_TYPE_MSVC)

	  #if(0)	/* The "else" block below is the fully expaanded and argcount-reduced version of this macro sequence, used to generate the GCC inline ASM */

	/*...Block 1: */
		/* DIT radix-8 subconvolution, sans twiddles. */
		/* eax,ebx,ecx,edx = r00,r20,r10,r30: */
		__asm	mov eax, r00
		__asm	mov ebx, eax
		__asm	mov ecx, eax
		__asm	mov edx, eax
		__asm	add ebx, 0x200
		__asm	add ecx, 0x100
		__asm	add edx, 0x300
		SSE2_RADIX4_DIT_IN_PLACE         ()
		/* eax,ebx,ecx,edx = r08,r28,r18,r38: */
		__asm	add eax, 0x080
		__asm	add ebx, 0x080
		__asm	add ecx, 0x080
		__asm	add edx, 0x080
		SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()
		__asm	mov esi, eax
		__asm	sub esi, 0x040	/* r04 */
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38)
												/* Totals: 84 load/store [66 movaps], 68 add/subpd,  4 mulpd, 48 address-compute */
	/*...Block 2:	*/
		/* DIT radix-8 subconvolution, sans twiddles. */
		/* eax,ebx,ecx,edx = r04,r24,r14,r34: */
		__asm	mov eax, esi
		__asm	mov ebx, esi
		__asm	mov ecx, esi
		__asm	mov edx, esi
		__asm	add ebx, 0x200
		__asm	add ecx, 0x100
		__asm	add edx, 0x300
		SSE2_RADIX4_DIT_IN_PLACE         ()
		/* eax,ebx,ecx,edx = r0C,r2C,r1C,r3C: */
		__asm	add eax, 0x080
		__asm	add ebx, 0x080
		__asm	add ecx, 0x080
		__asm	add edx, 0x080
		SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()
		__asm	mov esi, eax
		__asm	sub esi, 0x0a0	/* r02 */
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C)

	/*...Block 3:	*/
		/* DIT radix-8 subconvolution, sans twiddles. */
		/* eax,ebx,ecx,edx = r02,r22,r12,r32: */
		__asm	mov eax, esi
		__asm	mov ebx, esi
		__asm	mov ecx, esi
		__asm	mov edx, esi
		__asm	add ebx, 0x200
		__asm	add ecx, 0x100
		__asm	add edx, 0x300
		SSE2_RADIX4_DIT_IN_PLACE         ()
		/* eax,ebx,ecx,edx = r0A,r2A,r1A,r3A: */
		__asm	add eax, 0x080
		__asm	add ebx, 0x080
		__asm	add ecx, 0x080
		__asm	add edx, 0x080
		SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()
		__asm	mov esi, eax
		__asm	sub esi, 0x040	/* r04 */
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A)

	/*...Block 4:	*/
		/* DIT radix-8 subconvolution, sans twiddles. */
		/* eax,ebx,ecx,edx = r06,r26,r16,r36: */
		__asm	mov eax, esi
		__asm	mov ebx, esi
		__asm	mov ecx, esi
		__asm	mov edx, esi
		__asm	add ebx, 0x200
		__asm	add ecx, 0x100
		__asm	add edx, 0x300
		SSE2_RADIX4_DIT_IN_PLACE         ()
		/* eax,ebx,ecx,edx = r0E,r2E,r1E,r3E: */
		__asm	add eax, 0x080
		__asm	add ebx, 0x080
		__asm	add ecx, 0x080
		__asm	add edx, 0x080
		SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E)

	  #else

	/*...Block 1: */
		/* DIT radix-8 subconvolution, sans twiddles. */
		/* eax,ebx,ecx,edx = r00,r20,r10,r30: */
		__asm	mov	esi, isrt2
		__asm	mov eax, r00
		__asm	mov ebx, eax
		__asm	mov ecx, eax
		__asm	mov edx, eax
		__asm	add ebx, 0x200
		__asm	add ecx, 0x100
		__asm	add edx, 0x300
	/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/
		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */

		__asm	addpd	xmm0,[ebx     ]	/* t1 */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 */
		__asm	subpd	xmm2,[ebx     ]	/* t3 */
		__asm	subpd	xmm3,[ebx+0x10]	/* t4 */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */

		__asm	addpd	xmm4,[edx     ]	/* t5 */
		__asm	addpd	xmm5,[edx+0x10]	/* t6 */
		__asm	subpd	xmm6,[edx     ]	/* t7 */
		__asm	subpd	xmm7,[edx+0x10]	/* t8 */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */
		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */
		__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */
		__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */
		__asm	addpd	xmm4,xmm4	/*          2*t5 */
		__asm	addpd	xmm5,xmm5	/*          2*t6 */
		__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */
		__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */
		__asm	movaps	[eax      ],xmm4	/* <- ~t1 */
		__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */

		__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */
		__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */
		__asm	movaps	[edx      ],xmm2	/* <- ~t7 */
		__asm	movaps	[ecx+0x010],xmm3	/* <- ~t4 */
		__asm	addpd	xmm7,xmm7	/*          2*t8 */
		__asm	addpd	xmm6,xmm6	/*          2*t7 */
		__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */
		__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */
		__asm	movaps	[ecx      ],xmm7	/* <- ~t3 */
		__asm	movaps	[edx+0x010],xmm6	/* <- ~t8 */

		/* eax,ebx,ecx,edx = r08,r28,r18,r38: */
		__asm	add eax, 0x080
		__asm	add ebx, 0x080
		__asm	add ecx, 0x080
		__asm	add edx, 0x080
		/*****	SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()	*****/
		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */

		__asm	addpd	xmm0,[ebx     ]	/* t1 */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 */
		__asm	subpd	xmm2,[ebx     ]	/* t3 */
		__asm	subpd	xmm3,[ebx+0x10]	/* t4 */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */

		__asm	addpd	xmm4,[edx     ]	/* t5 */
		__asm	addpd	xmm5,[edx+0x10]	/* t6 */
		__asm	subpd	xmm6,[edx     ]	/* t7 */
		__asm	subpd	xmm7,[edx+0x10]	/* t8 */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */

		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */
		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */
		__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */
		__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */
		__asm	addpd	xmm4,xmm4	/*          2*t5 */
		__asm	addpd	xmm5,xmm5	/*          2*t6 */
		__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */
		__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */
		__asm	movaps	[eax      ],xmm4	/* <- ~t1 */
		__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */

		__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */
		__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */
		__asm	addpd	xmm7,xmm7	/*          2*t8 */
		__asm	addpd	xmm6,xmm6	/*          2*t7 */
		__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */
		__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */

		/*
		t3 =(t3+t4)*ISRT2;t4 =(t4-t3)*ISRT2;
		t7 =(t7+t8)*ISRT2;t8 =(t8-t7)*ISRT2;
		*/
		__asm	movaps	xmm0,xmm3	/* cpy t4 */
		__asm	movaps	xmm1,xmm6	/* cpy t8 */
		__asm	subpd	xmm3,xmm7	/* 4-3*/
		__asm	subpd	xmm6,xmm2	/* 8-7*/
		__asm	addpd	xmm0,xmm7	/* 4+3*/
		__asm	addpd	xmm1,xmm2	/* 8+7*/
		__asm	mulpd	xmm3,[esi]	/* (4-3)*ISRT2 */
		__asm	mulpd	xmm6,[esi]	/* (8-7)*ISRT2 */
		__asm	mulpd	xmm0,[esi]	/* (4+3)*ISRT2 */
		__asm	mulpd	xmm1,[esi]	/* (8+7)*ISRT2 */
		__asm	movaps	[ecx+0x010],xmm3	/* a[jp+p4 ] <- ~t4 */
		__asm	movaps	[edx+0x010],xmm6	/* a[jp+p12] <- ~t8 */
		__asm	movaps	[ecx      ],xmm0	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[edx      ],xmm1	/* a[jt+p12] <- ~t7 */

		/* Combine the 2 radix-4 subtransforms: */
	/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38)	*****/
		/* eax,ebx,ecx,edx contain r08,r28,r18,r38, resp: */
		__asm	movaps	xmm0,[eax-0x080]	/* t0 */		__asm	movaps	xmm4,[ebx-0x080]	/* t4 */
		__asm	movaps	xmm1,[eax-0x070]	/* t1 */		__asm	movaps	xmm5,[ebx-0x070]	/* t5 */
		__asm	movaps	xmm2,[eax      ]	/* cpy t8 */	__asm	movaps	xmm7,[ebx+0x010]	/* td */
		__asm	movaps	xmm3,[eax+0x010]	/* cpy t9 */	__asm	movaps	xmm6,[ebx      ]	/* tc */
		__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* tc = 4-d */
		__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* t5 = 5-c */
		__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */
		__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */
		__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* t4 = 4+d */
		__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* td = 5+c */

		__asm	movaps	[eax      ],xmm0	/* t8 */		__asm	movaps	[ebx      ],xmm4	/* t4 */
		__asm	movaps	[eax+0x010],xmm1	/* t9 */		__asm	movaps	[ebx-0x070],xmm5	/* td */
		__asm	movaps	[eax-0x080],xmm2	/* t0 */		__asm	movaps	[ebx-0x080],xmm7	/* tc */
		__asm	movaps	[eax-0x070],xmm3	/* t1 */		__asm	movaps	[ebx+0x010],xmm6	/* t5 */

		__asm	movaps	xmm0,[ecx-0x080]	/* t2 */		__asm	movaps	xmm4,[edx-0x080]	/* t6 */
		__asm	movaps	xmm1,[ecx-0x070]	/* t3 */		__asm	movaps	xmm5,[edx-0x070]	/* t7 */
		__asm	movaps	xmm2,[ecx      ]	/* cpy ta */	__asm	movaps	xmm7,[edx+0x010]	/* tf */
		__asm	movaps	xmm3,[ecx+0x010]	/* cpy tb */	__asm	movaps	xmm6,[edx      ]	/* te */
		__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* te = 6-f */
		__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* t7 = 7-e */
		__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */
		__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */
		__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* t6 = 6+f */
		__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* tf = 7+e */

		__asm	movaps	[ecx      ],xmm0	/* ta */		__asm	movaps	[edx      ],xmm4	/* t6 */
		__asm	movaps	[ecx+0x010],xmm1	/* tb */		__asm	movaps	[edx-0x070],xmm5	/* tf */
		__asm	movaps	[ecx-0x080],xmm2	/* t2 */		__asm	movaps	[edx-0x080],xmm7	/* te */
		__asm	movaps	[ecx-0x070],xmm3	/* t3 */		__asm	movaps	[edx+0x010],xmm6	/* t7 */

	/*...Block 2:	*/
		/* DIT radix-8 subconvolution, sans twiddles. */
		/* eax,ebx,ecx,edx = r04,r24,r14,r34: */
		__asm	sub eax, 0x040
		__asm	sub ebx, 0x040
		__asm	sub ecx, 0x040
		__asm	sub edx, 0x040
		/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/
		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */

		__asm	addpd	xmm0,[ebx     ]	/* t1 */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 */
		__asm	subpd	xmm2,[ebx     ]	/* t3 */
		__asm	subpd	xmm3,[ebx+0x10]	/* t4 */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */

		__asm	addpd	xmm4,[edx     ]	/* t5 */
		__asm	addpd	xmm5,[edx+0x10]	/* t6 */
		__asm	subpd	xmm6,[edx     ]	/* t7 */
		__asm	subpd	xmm7,[edx+0x10]	/* t8 */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */
		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */
		__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */
		__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */
		__asm	addpd	xmm4,xmm4	/*          2*t5 */
		__asm	addpd	xmm5,xmm5	/*          2*t6 */
		__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */
		__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */
		__asm	movaps	[eax      ],xmm4	/* <- ~t1 */
		__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */

		__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */
		__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */
		__asm	movaps	[edx      ],xmm2	/* <- ~t7 */
		__asm	movaps	[ecx+0x010],xmm3	/* <- ~t4 */
		__asm	addpd	xmm7,xmm7	/*          2*t8 */
		__asm	addpd	xmm6,xmm6	/*          2*t7 */
		__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */
		__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */
		__asm	movaps	[ecx      ],xmm7	/* <- ~t3 */
		__asm	movaps	[edx+0x010],xmm6	/* <- ~t8 */

		/* eax,ebx,ecx,edx = r0C,r2C,r1C,r3C: */
		__asm	add eax, 0x080
		__asm	add ebx, 0x080
		__asm	add ecx, 0x080
		__asm	add edx, 0x080
		/*****	SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()	*****/
		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */

		__asm	addpd	xmm0,[ebx     ]	/* t1 */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 */
		__asm	subpd	xmm2,[ebx     ]	/* t3 */
		__asm	subpd	xmm3,[ebx+0x10]	/* t4 */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */

		__asm	addpd	xmm4,[edx     ]	/* t5 */
		__asm	addpd	xmm5,[edx+0x10]	/* t6 */
		__asm	subpd	xmm6,[edx     ]	/* t7 */
		__asm	subpd	xmm7,[edx+0x10]	/* t8 */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */

		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */
		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */
		__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */
		__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */
		__asm	addpd	xmm4,xmm4	/*          2*t5 */
		__asm	addpd	xmm5,xmm5	/*          2*t6 */
		__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */
		__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */
		__asm	movaps	[eax      ],xmm4	/* <- ~t1 */
		__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */

		__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */
		__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */
		__asm	addpd	xmm7,xmm7	/*          2*t8 */
		__asm	addpd	xmm6,xmm6	/*          2*t7 */
		__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */
		__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */

		/*
		t3 =(t3+t4)*ISRT2;t4 =(t4-t3)*ISRT2;
		t7 =(t7+t8)*ISRT2;t8 =(t8-t7)*ISRT2;
		*/
		__asm	movaps	xmm0,xmm3	/* cpy t4 */
		__asm	movaps	xmm1,xmm6	/* cpy t8 */
		__asm	subpd	xmm3,xmm7	/* 4-3*/
		__asm	subpd	xmm6,xmm2	/* 8-7*/
		__asm	addpd	xmm0,xmm7	/* 4+3*/
		__asm	addpd	xmm1,xmm2	/* 8+7*/
		__asm	mulpd	xmm3,[esi]	/* (4-3)*ISRT2 */
		__asm	mulpd	xmm6,[esi]	/* (8-7)*ISRT2 */
		__asm	mulpd	xmm0,[esi]	/* (4+3)*ISRT2 */
		__asm	mulpd	xmm1,[esi]	/* (8+7)*ISRT2 */
		__asm	movaps	[ecx+0x010],xmm3	/* a[jp+p4 ] <- ~t4 */
		__asm	movaps	[edx+0x010],xmm6	/* a[jp+p12] <- ~t8 */
		__asm	movaps	[ecx      ],xmm0	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[edx      ],xmm1	/* a[jt+p12] <- ~t7 */

		/* Combine the 2 radix-4 subtransforms: */
	/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C)	*****/
		/* eax,ebx,ecx,edx contain r0C,r2C,r1C,r3C, resp: */
		__asm	movaps	xmm0,[eax-0x080]	/* t0 */		__asm	movaps	xmm4,[ebx-0x080]	/* t4 */
		__asm	movaps	xmm1,[eax-0x070]	/* t1 */		__asm	movaps	xmm5,[ebx-0x070]	/* t5 */
		__asm	movaps	xmm2,[eax      ]	/* cpy t8 */	__asm	movaps	xmm7,[ebx+0x010]	/* td */
		__asm	movaps	xmm3,[eax+0x010]	/* cpy t9 */	__asm	movaps	xmm6,[ebx      ]	/* tc */
		__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* tc = 4-d */
		__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* t5 = 5-c */
		__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */
		__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */
		__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* t4 = 4+d */
		__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* td = 5+c */

		__asm	movaps	[eax      ],xmm0	/* t8 */		__asm	movaps	[ebx      ],xmm4	/* t4 */
		__asm	movaps	[eax+0x010],xmm1	/* t9 */		__asm	movaps	[ebx-0x070],xmm5	/* td */
		__asm	movaps	[eax-0x080],xmm2	/* t0 */		__asm	movaps	[ebx-0x080],xmm7	/* tc */
		__asm	movaps	[eax-0x070],xmm3	/* t1 */		__asm	movaps	[ebx+0x010],xmm6	/* t5 */

		__asm	movaps	xmm0,[ecx-0x080]	/* t2 */		__asm	movaps	xmm4,[edx-0x080]	/* t6 */
		__asm	movaps	xmm1,[ecx-0x070]	/* t3 */		__asm	movaps	xmm5,[edx-0x070]	/* t7 */
		__asm	movaps	xmm2,[ecx      ]	/* cpy ta */	__asm	movaps	xmm7,[edx+0x010]	/* tf */
		__asm	movaps	xmm3,[ecx+0x010]	/* cpy tb */	__asm	movaps	xmm6,[edx      ]	/* te */
		__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* te = 6-f */
		__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* t7 = 7-e */
		__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */
		__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */
		__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* t6 = 6+f */
		__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* tf = 7+e */

		__asm	movaps	[ecx      ],xmm0	/* ta */		__asm	movaps	[edx      ],xmm4	/* t6 */
		__asm	movaps	[ecx+0x010],xmm1	/* tb */		__asm	movaps	[edx-0x070],xmm5	/* tf */
		__asm	movaps	[ecx-0x080],xmm2	/* t2 */		__asm	movaps	[edx-0x080],xmm7	/* te */
		__asm	movaps	[ecx-0x070],xmm3	/* t3 */		__asm	movaps	[edx+0x010],xmm6	/* t7 */

	/*...Block 3:	*/
		/* DIT radix-8 subconvolution, sans twiddles. */
		/* eax,ebx,ecx,edx = r02,r22,r12,r32: */
		__asm	sub eax, 0x0a0
		__asm	sub ebx, 0x0a0
		__asm	sub ecx, 0x0a0
		__asm	sub edx, 0x0a0
		/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/
		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */

		__asm	addpd	xmm0,[ebx     ]	/* t1 */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 */
		__asm	subpd	xmm2,[ebx     ]	/* t3 */
		__asm	subpd	xmm3,[ebx+0x10]	/* t4 */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */

		__asm	addpd	xmm4,[edx     ]	/* t5 */
		__asm	addpd	xmm5,[edx+0x10]	/* t6 */
		__asm	subpd	xmm6,[edx     ]	/* t7 */
		__asm	subpd	xmm7,[edx+0x10]	/* t8 */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */
		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */
		__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */
		__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */
		__asm	addpd	xmm4,xmm4	/*          2*t5 */
		__asm	addpd	xmm5,xmm5	/*          2*t6 */
		__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */
		__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */
		__asm	movaps	[eax      ],xmm4	/* <- ~t1 */
		__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */

		__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */
		__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */
		__asm	movaps	[edx      ],xmm2	/* <- ~t7 */
		__asm	movaps	[ecx+0x010],xmm3	/* <- ~t4 */
		__asm	addpd	xmm7,xmm7	/*          2*t8 */
		__asm	addpd	xmm6,xmm6	/*          2*t7 */
		__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */
		__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */
		__asm	movaps	[ecx      ],xmm7	/* <- ~t3 */
		__asm	movaps	[edx+0x010],xmm6	/* <- ~t8 */

		/* eax,ebx,ecx,edx = r0A,r2A,r1A,r3A: */
		__asm	add eax, 0x080
		__asm	add ebx, 0x080
		__asm	add ecx, 0x080
		__asm	add edx, 0x080
		/*****	SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()	*****/
		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */

		__asm	addpd	xmm0,[ebx     ]	/* t1 */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 */
		__asm	subpd	xmm2,[ebx     ]	/* t3 */
		__asm	subpd	xmm3,[ebx+0x10]	/* t4 */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */

		__asm	addpd	xmm4,[edx     ]	/* t5 */
		__asm	addpd	xmm5,[edx+0x10]	/* t6 */
		__asm	subpd	xmm6,[edx     ]	/* t7 */
		__asm	subpd	xmm7,[edx+0x10]	/* t8 */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */

		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */
		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */
		__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */
		__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */
		__asm	addpd	xmm4,xmm4	/*          2*t5 */
		__asm	addpd	xmm5,xmm5	/*          2*t6 */
		__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */
		__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */
		__asm	movaps	[eax      ],xmm4	/* <- ~t1 */
		__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */

		__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */
		__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */
		__asm	addpd	xmm7,xmm7	/*          2*t8 */
		__asm	addpd	xmm6,xmm6	/*          2*t7 */
		__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */
		__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */

		/*
		t3 =(t3+t4)*ISRT2;t4 =(t4-t3)*ISRT2;
		t7 =(t7+t8)*ISRT2;t8 =(t8-t7)*ISRT2;
		*/
		__asm	movaps	xmm0,xmm3	/* cpy t4 */
		__asm	movaps	xmm1,xmm6	/* cpy t8 */
		__asm	subpd	xmm3,xmm7	/* 4-3*/
		__asm	subpd	xmm6,xmm2	/* 8-7*/
		__asm	addpd	xmm0,xmm7	/* 4+3*/
		__asm	addpd	xmm1,xmm2	/* 8+7*/
		__asm	mulpd	xmm3,[esi]	/* (4-3)*ISRT2 */
		__asm	mulpd	xmm6,[esi]	/* (8-7)*ISRT2 */
		__asm	mulpd	xmm0,[esi]	/* (4+3)*ISRT2 */
		__asm	mulpd	xmm1,[esi]	/* (8+7)*ISRT2 */
		__asm	movaps	[ecx+0x010],xmm3	/* a[jp+p4 ] <- ~t4 */
		__asm	movaps	[edx+0x010],xmm6	/* a[jp+p12] <- ~t8 */
		__asm	movaps	[ecx      ],xmm0	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[edx      ],xmm1	/* a[jt+p12] <- ~t7 */

		/* Combine the 2 radix-4 subtransforms: */
	/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A)	*****/
		/* eax,ebx,ecx,edx contain r0A,r2A,r1A,r3A, resp: */
		__asm	movaps	xmm0,[eax-0x080]	/* t0 */		__asm	movaps	xmm4,[ebx-0x080]	/* t4 */
		__asm	movaps	xmm1,[eax-0x070]	/* t1 */		__asm	movaps	xmm5,[ebx-0x070]	/* t5 */
		__asm	movaps	xmm2,[eax      ]	/* cpy t8 */	__asm	movaps	xmm7,[ebx+0x010]	/* td */
		__asm	movaps	xmm3,[eax+0x010]	/* cpy t9 */	__asm	movaps	xmm6,[ebx      ]	/* tc */
		__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* tc = 4-d */
		__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* t5 = 5-c */
		__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */
		__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */
		__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* t4 = 4+d */
		__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* td = 5+c */

		__asm	movaps	[eax      ],xmm0	/* t8 */		__asm	movaps	[ebx      ],xmm4	/* t4 */
		__asm	movaps	[eax+0x010],xmm1	/* t9 */		__asm	movaps	[ebx-0x070],xmm5	/* td */
		__asm	movaps	[eax-0x080],xmm2	/* t0 */		__asm	movaps	[ebx-0x080],xmm7	/* tc */
		__asm	movaps	[eax-0x070],xmm3	/* t1 */		__asm	movaps	[ebx+0x010],xmm6	/* t5 */

		__asm	movaps	xmm0,[ecx-0x080]	/* t2 */		__asm	movaps	xmm4,[edx-0x080]	/* t6 */
		__asm	movaps	xmm1,[ecx-0x070]	/* t3 */		__asm	movaps	xmm5,[edx-0x070]	/* t7 */
		__asm	movaps	xmm2,[ecx      ]	/* cpy ta */	__asm	movaps	xmm7,[edx+0x010]	/* tf */
		__asm	movaps	xmm3,[ecx+0x010]	/* cpy tb */	__asm	movaps	xmm6,[edx      ]	/* te */
		__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* te = 6-f */
		__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* t7 = 7-e */
		__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */
		__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */
		__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* t6 = 6+f */
		__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* tf = 7+e */

		__asm	movaps	[ecx      ],xmm0	/* ta */		__asm	movaps	[edx      ],xmm4	/* t6 */
		__asm	movaps	[ecx+0x010],xmm1	/* tb */		__asm	movaps	[edx-0x070],xmm5	/* tf */
		__asm	movaps	[ecx-0x080],xmm2	/* t2 */		__asm	movaps	[edx-0x080],xmm7	/* te */
		__asm	movaps	[ecx-0x070],xmm3	/* t3 */		__asm	movaps	[edx+0x010],xmm6	/* t7 */

	/*...Block 4:	*/
		/* DIT radix-8 subconvolution, sans twiddles. */
		/* eax,ebx,ecx,edx = r06,r26,r16,r36: */
		__asm	sub eax, 0x040
		__asm	sub ebx, 0x040
		__asm	sub ecx, 0x040
		__asm	sub edx, 0x040
		/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/
		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */

		__asm	addpd	xmm0,[ebx     ]	/* t1 */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 */
		__asm	subpd	xmm2,[ebx     ]	/* t3 */
		__asm	subpd	xmm3,[ebx+0x10]	/* t4 */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */

		__asm	addpd	xmm4,[edx     ]	/* t5 */
		__asm	addpd	xmm5,[edx+0x10]	/* t6 */
		__asm	subpd	xmm6,[edx     ]	/* t7 */
		__asm	subpd	xmm7,[edx+0x10]	/* t8 */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */
		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */
		__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */
		__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */
		__asm	addpd	xmm4,xmm4	/*          2*t5 */
		__asm	addpd	xmm5,xmm5	/*          2*t6 */
		__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */
		__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */
		__asm	movaps	[eax      ],xmm4	/* <- ~t1 */
		__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */

		__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */
		__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */
		__asm	movaps	[edx      ],xmm2	/* <- ~t7 */
		__asm	movaps	[ecx+0x010],xmm3	/* <- ~t4 */
		__asm	addpd	xmm7,xmm7	/*          2*t8 */
		__asm	addpd	xmm6,xmm6	/*          2*t7 */
		__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */
		__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */
		__asm	movaps	[ecx      ],xmm7	/* <- ~t3 */
		__asm	movaps	[edx+0x010],xmm6	/* <- ~t8 */

		/* eax,ebx,ecx,edx = r0E,r2E,r1E,r3E: */
		__asm	add eax, 0x080
		__asm	add ebx, 0x080
		__asm	add ecx, 0x080
		__asm	add edx, 0x080
		/*****	SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()	*****/
		__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */

		__asm	addpd	xmm0,[ebx     ]	/* t1 */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 */
		__asm	subpd	xmm2,[ebx     ]	/* t3 */
		__asm	subpd	xmm3,[ebx+0x10]	/* t4 */

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */

		__asm	addpd	xmm4,[edx     ]	/* t5 */
		__asm	addpd	xmm5,[edx+0x10]	/* t6 */
		__asm	subpd	xmm6,[edx     ]	/* t7 */
		__asm	subpd	xmm7,[edx+0x10]	/* t8 */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */

		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */
		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */
		__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */
		__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */
		__asm	addpd	xmm4,xmm4	/*          2*t5 */
		__asm	addpd	xmm5,xmm5	/*          2*t6 */
		__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */
		__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */
		__asm	movaps	[eax      ],xmm4	/* <- ~t1 */
		__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */

		__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */
		__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */
		__asm	addpd	xmm7,xmm7	/*          2*t8 */
		__asm	addpd	xmm6,xmm6	/*          2*t7 */
		__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */
		__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */

		/*
		t3 =(t3+t4)*ISRT2;t4 =(t4-t3)*ISRT2;
		t7 =(t7+t8)*ISRT2;t8 =(t8-t7)*ISRT2;
		*/
		__asm	movaps	xmm0,xmm3	/* cpy t4 */
		__asm	movaps	xmm1,xmm6	/* cpy t8 */
		__asm	subpd	xmm3,xmm7	/* 4-3*/
		__asm	subpd	xmm6,xmm2	/* 8-7*/
		__asm	addpd	xmm0,xmm7	/* 4+3*/
		__asm	addpd	xmm1,xmm2	/* 8+7*/
		__asm	mulpd	xmm3,[esi]	/* (4-3)*ISRT2 */
		__asm	mulpd	xmm6,[esi]	/* (8-7)*ISRT2 */
		__asm	mulpd	xmm0,[esi]	/* (4+3)*ISRT2 */
		__asm	mulpd	xmm1,[esi]	/* (8+7)*ISRT2 */
		__asm	movaps	[ecx+0x010],xmm3	/* a[jp+p4 ] <- ~t4 */
		__asm	movaps	[edx+0x010],xmm6	/* a[jp+p12] <- ~t8 */
		__asm	movaps	[ecx      ],xmm0	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[edx      ],xmm1	/* a[jt+p12] <- ~t7 */

		/* Combine the 2 radix-4 subtransforms: */
	/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E)	*****/
		/* eax,ebx,ecx,edx contain r0E,r2E,r1E,r3E, resp: */
		__asm	movaps	xmm0,[eax-0x080]	/* t0 */		__asm	movaps	xmm4,[ebx-0x080]	/* t4 */
		__asm	movaps	xmm1,[eax-0x070]	/* t1 */		__asm	movaps	xmm5,[ebx-0x070]	/* t5 */
		__asm	movaps	xmm2,[eax      ]	/* cpy t8 */	__asm	movaps	xmm7,[ebx+0x010]	/* td */
		__asm	movaps	xmm3,[eax+0x010]	/* cpy t9 */	__asm	movaps	xmm6,[ebx      ]	/* tc */
		__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* tc = 4-d */
		__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* t5 = 5-c */
		__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */
		__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */
		__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* t4 = 4+d */
		__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* td = 5+c */

		__asm	movaps	[eax      ],xmm0	/* t8 */		__asm	movaps	[ebx      ],xmm4	/* t4 */
		__asm	movaps	[eax+0x010],xmm1	/* t9 */		__asm	movaps	[ebx-0x070],xmm5	/* td */
		__asm	movaps	[eax-0x080],xmm2	/* t0 */		__asm	movaps	[ebx-0x080],xmm7	/* tc */
		__asm	movaps	[eax-0x070],xmm3	/* t1 */		__asm	movaps	[ebx+0x010],xmm6	/* t5 */

		__asm	movaps	xmm0,[ecx-0x080]	/* t2 */		__asm	movaps	xmm4,[edx-0x080]	/* t6 */
		__asm	movaps	xmm1,[ecx-0x070]	/* t3 */		__asm	movaps	xmm5,[edx-0x070]	/* t7 */
		__asm	movaps	xmm2,[ecx      ]	/* cpy ta */	__asm	movaps	xmm7,[edx+0x010]	/* tf */
		__asm	movaps	xmm3,[ecx+0x010]	/* cpy tb */	__asm	movaps	xmm6,[edx      ]	/* te */
		__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* te = 6-f */
		__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* t7 = 7-e */
		__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */
		__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */
		__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* t6 = 6+f */
		__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* tf = 7+e */

		__asm	movaps	[ecx      ],xmm0	/* ta */		__asm	movaps	[edx      ],xmm4	/* t6 */
		__asm	movaps	[ecx+0x010],xmm1	/* tb */		__asm	movaps	[edx-0x070],xmm5	/* tf */
		__asm	movaps	[ecx-0x080],xmm2	/* t2 */		__asm	movaps	[edx-0x080],xmm7	/* te */
		__asm	movaps	[ecx-0x070],xmm3	/* t3 */		__asm	movaps	[edx+0x010],xmm6	/* t7 */

	  #endif	/* if(1) */

	/****************************************************************************************************
	!......Now do eight more radix-4 transforms, including the internal and external twiddle factors.   !
	!   Write even-index 16-byte output pairs to a[j1], odd-index to a[j2], unpack same as on inputs.   !
	!   We do the last 4 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.     !
	****************************************************************************************************/

	/* Main-array addresses still in add0,1, no need to re-init:
		add0 = &a[j1pad];
		add1 = &a[j2pad];
	*/
	/************************************************/
	/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16:	*/
	/************************************************/
		__asm	mov	esi, isrt2
		__asm	mov	eax, r10
		__asm	mov	ebx, add1	/* restore main-array index [now into j2 block], but keep r10 in eax for temporary storage */
		__asm	mov	ecx, esi
		__asm	mov	edx, esi
		__asm	mov	edi, c01
		__asm	add	ecx, 0x10	/* cc0 */
		__asm	add	edx, 0x30	/* cc1 */

		__asm	movaps	xmm4,[eax+0x020]	/* t22 */				__asm	movaps	xmm0,[eax+0x060]	/* t32 */
		__asm	movaps	xmm5,[eax+0x030]	/* t23 */				__asm	movaps	xmm1,[eax+0x070]	/* t33 */
		__asm	movaps	xmm6,[eax+0x020]	/* xmm2 <- cpy t22 */	__asm	movaps	xmm2,[eax+0x060]	/* xmm6 <- cpy t32 */
		__asm	movaps	xmm7,[eax+0x030]	/* xmm3 <- cpy t23 */	__asm	movaps	xmm3,[eax+0x070]	/* xmm7 <- cpy t33 */

		__asm	mulpd	xmm4,[edx     ]	/* t22*c32_1 */				__asm	mulpd	xmm0,[edx+0x20]	/* t32*c32_3 */
		__asm	mulpd	xmm5,[edx     ]	/* t23*c32_1 */				__asm	mulpd	xmm1,[edx+0x20]	/* t33*c32_3 */
		__asm	mulpd	xmm6,[edx+0x10]	/* t22*s32_1 */				__asm	mulpd	xmm2,[edx+0x30]	/* t32*s32_3 */
		__asm	mulpd	xmm7,[edx+0x10]	/* t23*s32_1 */				__asm	mulpd	xmm3,[edx+0x30]	/* t33*s32_3 */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t23 */				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t22 */				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t23 */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t22 */

		__asm	subpd	xmm4,xmm0	/* ~t32 <- t22-rt */
		__asm	subpd	xmm5,xmm1	/* ~t33 <- t23-it */
		__asm	addpd	xmm6,xmm0	/* ~t22 <- t22+rt */
		__asm	addpd	xmm7,xmm1	/* ~t23 <- t23+it */

		__asm	movaps	xmm2,[eax+0x040]	/* t12 */
		__asm	movaps	xmm3,[eax+0x050]	/* t13 */
		__asm	movaps	xmm0,[eax+0x040]	/* cpy t12 */
		__asm	movaps	xmm1,[eax+0x050]	/* cpy t13 */

		__asm	mulpd	xmm2,[ecx     ]	/* t12*c */
		__asm	mulpd	xmm1,[ecx+0x10]	/* t13*s */
		__asm	mulpd	xmm3,[ecx     ]	/* t13*c */
		__asm	mulpd	xmm0,[ecx+0x10]	/* t12*s */
		__asm	addpd	xmm2,xmm1	/* xmm2 <- rt */
		__asm	subpd	xmm3,xmm0	/* xmm3 <- it */

		__asm	movaps	xmm0,[eax      ]	/* t02 */
		__asm	movaps	xmm1,[eax+0x010]	/* t03 */

		__asm	subpd	xmm0,xmm2	/*~t12 <- t02- rt */
		__asm	subpd	xmm1,xmm3	/*~t13 <- t03- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t02 <- t02+ rt */
		__asm	addpd	xmm3,xmm1	/*~t03 <- t03+ it */
	/*
		rt       =t02+t22;			it       =t03+t23;
		t22      =t02-t22;			t23      =t03-t23;
		a[jt    ]=rt *c01+it *s01	a[jp    ]=it *c01-rt *s01
		a[jt+p8 ]=t22*c11+t23*s11;	a[jp+p8 ]=t22*c11-t23*s11;
	*/
		__asm	subpd	xmm2,xmm6		/*~t22 <- t02-t22 */
		__asm	subpd	xmm3,xmm7		/*~t23 <- t03-t23 */
		__asm	addpd	xmm6,xmm6		/*          2*t22 */
		__asm	addpd	xmm7,xmm7		/*          2*t23 */
		__asm	addpd	xmm6,xmm2		/* rt  <- t02+t22 */
		__asm	addpd	xmm7,xmm3		/* it  <- t03+t23 */
		__asm	movaps	[eax+0x020],xmm2	/* store ~t22 */
		__asm	movaps	[eax+0x030],xmm3	/* store ~t23 */
		__asm	movaps	xmm2,xmm6		/* rt copy */
		__asm	movaps	xmm3,xmm7		/* it copy */
		__asm	mulpd	xmm6,[edi     ]	/* rt *c01 */
		__asm	mulpd	xmm7,[edi     ]	/* it *c01 */
		__asm	mulpd	xmm2,[edi+0x10]	/* rt *s01 */
		__asm	mulpd	xmm3,[edi+0x10]	/* it *s01 */
		__asm	subpd	xmm7,xmm2	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm3	/* xmm6 <- re */

		__asm	movaps	[ebx+0x010],xmm7	/* a[jp    ] */
		__asm	movaps	[ebx      ],xmm6	/* a[jt    ] */

		__asm	movaps	xmm6,[eax+0x020]	/* load ~t22 */
		__asm	movaps	xmm7,[eax+0x030]	/* load ~t23 */
		__asm	movaps	xmm2,xmm6		/* re copy */
		__asm	movaps	xmm3,xmm7		/* im copy */
		__asm	mulpd	xmm6,[edi+0x20]	/* re *c11 */
		__asm	mulpd	xmm7,[edi+0x20]	/* im *c11 */
		__asm	mulpd	xmm2,[edi+0x30]	/* re *s11 */
		__asm	mulpd	xmm3,[edi+0x30]	/* im *s11 */
		__asm	subpd	xmm7,xmm2	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm3	/* xmm6 <- re */

		__asm	movaps	[ebx+0x110],xmm7	/* a[jp+p8 ] */
		__asm	movaps	[ebx+0x100],xmm6	/* a[jt+p8 ] */

	/* mpy by E^-4 = -I is inlined here...
		rt       =t12+t33;			it         =t13-t32;
		t33      =t12-t33;			t32        =t13+t32;
		a[jt+p08]=rt *c09+it *s09;	a[jp+p08]=it *c09-rt *s09;
		a[jt+p18]=t33*c19+t32*s19;	a[jp+p18]=t32*c19-t33*s19;
	*/
		__asm	add	edi, 0x040	/* c09 */
		__asm	subpd	xmm0,xmm5		/*~t33 <- t12-t33 */
		__asm	subpd	xmm1,xmm4		/* it  <- t13-t32 */
		__asm	addpd	xmm5,xmm5		/*          2*t33 */
		__asm	addpd	xmm4,xmm4		/*          2*t32 */
		__asm	addpd	xmm5,xmm0		/* rt  <- t12+t33 */
		__asm	addpd	xmm4,xmm1		/*~t32 <- t13+t32 */
		__asm	movaps	xmm2,xmm5		/* rt copy */
		__asm	movaps	xmm3,xmm1		/* it copy */
		__asm	mulpd	xmm5,[edi     ]	/* rt*c09 */
		__asm	mulpd	xmm1,[edi     ]	/* it*c09 */
		__asm	mulpd	xmm2,[edi+0x10]	/* rt*s09 */
		__asm	mulpd	xmm3,[edi+0x10]	/* it*s09 */
		__asm	subpd	xmm1,xmm2	/* xmm1 <- im */
		__asm	addpd	xmm5,xmm3	/* xmm5 <- re */

		__asm	movaps	[ebx+0x090],xmm1	/* a[jp+p4 ] */
		__asm	movaps	[ebx+0x080],xmm5	/* a[jt+p4 ] */

		__asm	movaps	xmm2,xmm0		/*t33 copy */
		__asm	movaps	xmm3,xmm4		/*t32 copy */
		__asm	mulpd	xmm0,[edi+0x20]	/*t33*c19 */
		__asm	mulpd	xmm4,[edi+0x20]	/*t32*c19 */
		__asm	mulpd	xmm2,[edi+0x30]	/*t33*s19 */
		__asm	mulpd	xmm3,[edi+0x30]	/*t32*s19 */
		__asm	subpd	xmm4,xmm2	/* xmm4 <- im */
		__asm	addpd	xmm0,xmm3	/* xmm0 <- re */

		__asm	movaps	[ebx+0x190],xmm4	/* a[jp+p12] */
		__asm	movaps	[ebx+0x180],xmm0	/* a[jt+p12] */

	/************************************************/
	/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:	*/
	/************************************************/
		__asm	add	eax, 0x080	/* r18 */

		__asm	movaps	xmm4,[eax+0x020]	/* t2A */				__asm	movaps	xmm0,[eax+0x060]	/* t3A */
		__asm	movaps	xmm5,[eax+0x030]	/* t2B */				__asm	movaps	xmm1,[eax+0x070]	/* t3B */
		__asm	movaps	xmm6,[eax+0x020]	/* xmm2 <- cpy t2A */	__asm	movaps	xmm2,[eax+0x060]	/* xmm6 <- cpy t3A */
		__asm	movaps	xmm7,[eax+0x030]	/* xmm3 <- cpy t2B */	__asm	movaps	xmm3,[eax+0x070]	/* xmm7 <- cpy t3B */

		__asm	mulpd	xmm4,[edx+0x30]	/* t2A*s32_3 */				__asm	mulpd	xmm0,[edx     ]	/* t3A*c32_1 */
		__asm	mulpd	xmm5,[edx+0x30]	/* t2B*s32_3 */				__asm	mulpd	xmm1,[edx     ]	/* t3B*c32_1 */
		__asm	mulpd	xmm6,[edx+0x20]	/* t2A*c32_3 */				__asm	mulpd	xmm2,[edx+0x10]	/* t3A*s32_1 */
		__asm	mulpd	xmm7,[edx+0x20]	/* t2B*c32_3 */				__asm	mulpd	xmm3,[edx+0x10]	/* t3B*s32_1 */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t2B */				__asm	addpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t2A */				__asm	subpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t2B */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t2A */

		__asm	addpd	xmm4,xmm0	/* ~t3A <- t2A+rt */
		__asm	addpd	xmm5,xmm1	/* ~t3B <- t2B+it */
		__asm	subpd	xmm6,xmm0	/* ~t2A <- t2A-rt */
		__asm	subpd	xmm7,xmm1	/* ~t2B <- t2B-it */

		__asm	movaps	xmm2,[eax+0x040]	/* t1A */
		__asm	movaps	xmm3,[eax+0x050]	/* t1B */
		__asm	movaps	xmm0,[eax+0x040]	/* cpy t1A */
		__asm	movaps	xmm1,[eax+0x050]	/* cpy t1B */

		__asm	mulpd	xmm2,[ecx+0x10]	/* t1A*s */
		__asm	mulpd	xmm1,[ecx     ]	/* t1B*c */
		__asm	mulpd	xmm3,[ecx+0x10]	/* t1B*s */
		__asm	mulpd	xmm0,[ecx     ]	/* t1A*c */
		__asm	subpd	xmm2,xmm1	/* xmmA <- rt */
		__asm	addpd	xmm3,xmm0	/* xmmB <- it */

		__asm	movaps	xmm0,[eax      ]	/* t0A */
		__asm	movaps	xmm1,[eax+0x010]	/* t0B */

		__asm	subpd	xmm0,xmm2	/*~t0A <- t0A- rt */
		__asm	subpd	xmm1,xmm3	/*~t0B <- t0B- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t1A <- t0A+ rt */
		__asm	addpd	xmm3,xmm1	/*~t1B <- t0B+ it */

	/*
		rt       =t0A+t2A;			it         =t0B+t2B;
		t2A      =t0A-t2A;			t2B        =t0B-t2B;
		a[jt    ]=rt *c05+it *s05;	a[jp    ]=it *c05-rt *s05;
		a[jt+p10]=t2A*c15+t2B*s15;	a[jp+p10]=t2B*c15-t2A*s15;
	*/
		__asm	add	edi, 0x040	/* c05 */
		__asm	subpd	xmm0,xmm6		/*~t2A <- t0A-t2A */
		__asm	subpd	xmm1,xmm7		/*~t2B <- t0B-t2B */
		__asm	addpd	xmm6,xmm6		/*          2*t2A */
		__asm	addpd	xmm7,xmm7		/*          2*t2B */
		__asm	addpd	xmm6,xmm0		/* rt  <- t0A+t2A */
		__asm	addpd	xmm7,xmm1		/* it  <- t0B+t2B */
		__asm	movaps	[eax+0x020],xmm0	/* store ~t2A */
		__asm	movaps	[eax+0x030],xmm1	/* store ~t2B */
		__asm	movaps	xmm0,xmm6		/* rt copy */
		__asm	movaps	xmm1,xmm7		/* it copy */
		__asm	mulpd	xmm6,[edi     ]	/* rt *c05 */
		__asm	mulpd	xmm7,[edi     ]	/* it *c05 */
		__asm	mulpd	xmm0,[edi+0x10]	/* rt *s05 */
		__asm	mulpd	xmm1,[edi+0x10]	/* it *s05 */
		__asm	subpd	xmm7,xmm0	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm1	/* xmm6 <- re */

		__asm	movaps	[ebx+0x050],xmm7	/* a[jp    ] */
		__asm	movaps	[ebx+0x040],xmm6	/* a[jt    ] */

		__asm	movaps	xmm6,[eax+0x020]	/* load ~t2A */
		__asm	movaps	xmm7,[eax+0x030]	/* load ~t2B */
		__asm	movaps	xmm0,xmm6		/* re copy */
		__asm	movaps	xmm1,xmm7		/* im copy */
		__asm	mulpd	xmm6,[edi+0x20]	/* re *c15 */
		__asm	mulpd	xmm7,[edi+0x20]	/* im *c15 */
		__asm	mulpd	xmm0,[edi+0x30]	/* re *s15 */
		__asm	mulpd	xmm1,[edi+0x30]	/* im *s15 */
		__asm	subpd	xmm7,xmm0	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm1	/* xmm6 <- re */

		__asm	movaps	[ebx+0x150],xmm7	/* a[jp+p8 ] */
		__asm	movaps	[ebx+0x140],xmm6	/* a[jt+p8 ] */
	/*
		rt       =t1A+t3B;		it         =t1B-t3A;
		t3B      =t1A-t3B;		t3A        =t1B+t3A;
		a[jt+p08]=rt *c0D+it *s0D;	a[jp+p08]=it *c0D-rt *s0D;
		a[jt+p18]=t3B*c1D+t3A*s1D;	a[jp+p18]=t3A*c1D-t3B*s1D;
	*/
		__asm	add	edi, 0x040	/* c0D */
		__asm	subpd	xmm2,xmm5		/*~t3B <- t1A-t3B */
		__asm	subpd	xmm3,xmm4		/* it  <- t1B-t3A */
		__asm	addpd	xmm5,xmm5		/*          2*t3B */
		__asm	addpd	xmm4,xmm4		/*          2*t3A */
		__asm	addpd	xmm5,xmm2		/* rt  <- t1A+t3B */
		__asm	addpd	xmm4,xmm3		/*~t3A <- t1B+t3A */
		__asm	movaps	xmm0,xmm5		/* rt copy */
		__asm	movaps	xmm1,xmm3		/* it copy */
		__asm	mulpd	xmm5,[edi     ]	/* rt*c0D */
		__asm	mulpd	xmm3,[edi     ]	/* it*c0D */
		__asm	mulpd	xmm0,[edi+0x10]	/* rt*s0D */
		__asm	mulpd	xmm1,[edi+0x10]	/* it*s0D */
		__asm	subpd	xmm3,xmm0	/* xmm3 <- im */
		__asm	addpd	xmm5,xmm1	/* xmm5 <- re */

		__asm	movaps	[ebx+0x0d0],xmm3	/* a[jp+p4 ] */
		__asm	movaps	[ebx+0x0c0],xmm5	/* a[jt+p4 ] */

		__asm	movaps	xmm0,xmm2		/*t3B copy */
		__asm	movaps	xmm1,xmm4		/*t3A copy */
		__asm	mulpd	xmm2,[edi+0x20]	/*t3B*c1D */
		__asm	mulpd	xmm4,[edi+0x20]	/*t3A*c1D */
		__asm	mulpd	xmm0,[edi+0x30]	/*t3B*s1D */
		__asm	mulpd	xmm1,[edi+0x30]	/*t3A*s1D */
		__asm	subpd	xmm4,xmm0	/* xmm4 <- im */
		__asm	addpd	xmm2,xmm1	/* xmm2 <- re */

		__asm	movaps	[ebx+0x1d0],xmm4	/* a[jp+p12] */
		__asm	movaps	[ebx+0x1c0],xmm2	/* a[jt+p12] */

	/************************************************/
	/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:	*/
	/************************************************/
		__asm	add	eax, 0x180	/* r30 */

		__asm	movaps	xmm4,[eax+0x020]	/* t26 */				__asm	movaps	xmm0,[eax+0x060]	/* t36 */
		__asm	movaps	xmm5,[eax+0x030]	/* t27 */				__asm	movaps	xmm1,[eax+0x070]	/* t37 */
		__asm	movaps	xmm6,[eax+0x020]	/* xmm2 <- cpy t26 */	__asm	movaps	xmm2,[eax+0x060]	/* xmm6 <- cpy t36 */
		__asm	movaps	xmm7,[eax+0x030]	/* xmm3 <- cpy t27 */	__asm	movaps	xmm3,[eax+0x070]	/* xmm7 <- cpy t37 */

		__asm	mulpd	xmm4,[edx+0x20]	/* t26*c32_3 */				__asm	mulpd	xmm0,[edx+0x10]	/* t36*s32_1 */
		__asm	mulpd	xmm5,[edx+0x20]	/* t27*c32_3 */				__asm	mulpd	xmm1,[edx+0x10]	/* t37*s32_1 */
		__asm	mulpd	xmm6,[edx+0x30]	/* t26*s32_3 */				__asm	mulpd	xmm2,[edx     ]	/* t36*c32_1 */
		__asm	mulpd	xmm7,[edx+0x30]	/* t27*s32_3 */				__asm	mulpd	xmm3,[edx     ]	/* t37*c32_1 */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t27 */				__asm	addpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t26 */				__asm	subpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t27 */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t26 */

		__asm	addpd	xmm4,xmm0	/* ~t36 <- t26+rt */
		__asm	addpd	xmm5,xmm1	/* ~t37 <- t27+it */
		__asm	subpd	xmm6,xmm0	/* ~t26 <- t26-rt */
		__asm	subpd	xmm7,xmm1	/* ~t27 <- t27-it */

		__asm	movaps	xmm2,[eax+0x040]	/* t16 */
		__asm	movaps	xmm3,[eax+0x050]	/* t17 */
		__asm	movaps	xmm0,[eax+0x040]	/* cpy t16 */
		__asm	movaps	xmm1,[eax+0x050]	/* cpy t17 */

		__asm	mulpd	xmm2,[ecx+0x10]	/* t16*s */
		__asm	mulpd	xmm1,[ecx     ]	/* t17*c */
		__asm	mulpd	xmm3,[ecx+0x10]	/* t17*s */
		__asm	mulpd	xmm0,[ecx     ]	/* t16*c */
		__asm	addpd	xmm2,xmm1	/* xmm2 <- rt */
		__asm	subpd	xmm3,xmm0	/* xmm3 <- it */

		__asm	movaps	xmm0,[eax      ]	/* t06 */
		__asm	movaps	xmm1,[eax+0x010]	/* t07 */

		__asm	subpd	xmm0,xmm2	/*~t16 <- t06- rt */
		__asm	subpd	xmm1,xmm3	/*~t17 <- t07- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t06 <- t06+ rt */
		__asm	addpd	xmm3,xmm1	/*~t07 <- t07+ it */

	/*
		rt       =t06+t26;		it         =t07+t27;
		t26      =t06-t26;		t27        =t07-t27;
		a[jt    ]=rt *c03+it *s03;	a[jp    ]=it *c03-rt *s03;
		a[jt+p10]=t26*c13+t27*s13;	a[jp+p10]=t27*c13-t26*s13;
	*/
		__asm	add	edi, 0x040	/* c03 */
		__asm	subpd	xmm2,xmm6		/*~t26 <- t06-t26 */
		__asm	subpd	xmm3,xmm7		/*~t27 <- t07-t27 */
		__asm	addpd	xmm6,xmm6		/*          2*t26 */
		__asm	addpd	xmm7,xmm7		/*          2*t27 */
		__asm	addpd	xmm6,xmm2		/* rt  <- t06+t26 */
		__asm	addpd	xmm7,xmm3		/* it  <- t07+t27 */
		__asm	movaps	[eax+0x020],xmm2	/* store ~t26 */
		__asm	movaps	[eax+0x030],xmm3	/* store ~t27 */
		__asm	movaps	xmm2,xmm6		/* rt copy */
		__asm	movaps	xmm3,xmm7		/* it copy */
		__asm	mulpd	xmm6,[edi     ]	/* rt *c03 */
		__asm	mulpd	xmm7,[edi     ]	/* it *c03 */
		__asm	mulpd	xmm2,[edi+0x10]	/* rt *s03 */
		__asm	mulpd	xmm3,[edi+0x10]	/* it *s03 */
		__asm	subpd	xmm7,xmm2	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm3	/* xmm6 <- re */

		__asm	movaps	[ebx+0x030],xmm7	/* a[jp    ] */
		__asm	movaps	[ebx+0x020],xmm6	/* a[jt    ] */

		__asm	movaps	xmm6,[eax+0x020]	/* load ~t26 */
		__asm	movaps	xmm7,[eax+0x030]	/* load ~t27 */
		__asm	movaps	xmm2,xmm6		/* re copy */
		__asm	movaps	xmm3,xmm7		/* im copy */
		__asm	mulpd	xmm6,[edi+0x20]	/* re *c13 */
		__asm	mulpd	xmm7,[edi+0x20]	/* im *c13 */
		__asm	mulpd	xmm2,[edi+0x30]	/* re *s13 */
		__asm	mulpd	xmm3,[edi+0x30]	/* im *s13 */
		__asm	subpd	xmm7,xmm2	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm3	/* xmm6 <- re */

		__asm	movaps	[ebx+0x130],xmm7	/* a[jp+p8 ] */
		__asm	movaps	[ebx+0x120],xmm6	/* a[jt+p8 ] */

	/* mpy by e^-4 = -I is inlined here...
		rt       =t16+t37;		it         =t17-t36;
		t37      =t16-t37;		t36        =t17+t36;
		a[jt+p08]=rt *c0B+it *s0B;	a[jp+p08]=it *c0B-rt *s0B;
		a[jt+p18]=t37*c1B+t36*s1B;	a[jp+p18]=t36*c1B-t37*s1B;
	*/
		__asm	add	edi, 0x040	/* c0B */
		__asm	subpd	xmm0,xmm5		/*~t37 <- t16-t37 */
		__asm	subpd	xmm1,xmm4		/* it  <- t17-t36 */
		__asm	addpd	xmm5,xmm5		/*          2*t37 */
		__asm	addpd	xmm4,xmm4		/*          2*t36 */
		__asm	addpd	xmm5,xmm0		/* rt  <- t16+t37 */
		__asm	addpd	xmm4,xmm1		/*~t36 <- t17+t36 */
		__asm	movaps	xmm2,xmm5		/* rt copy */
		__asm	movaps	xmm3,xmm1		/* it copy */
		__asm	mulpd	xmm5,[edi     ]	/* rt*c0B */
		__asm	mulpd	xmm1,[edi     ]	/* it*c0B */
		__asm	mulpd	xmm2,[edi+0x10]	/* rt*s0B */
		__asm	mulpd	xmm3,[edi+0x10]	/* it*s0B */
		__asm	subpd	xmm1,xmm2	/* xmm1 <- im */
		__asm	addpd	xmm5,xmm3	/* xmm5 <- re */

		__asm	movaps	[ebx+0x0b0],xmm1	/* a[jp+p4 ] */
		__asm	movaps	[ebx+0x0a0],xmm5	/* a[jt+p4 ] */

		__asm	movaps	xmm2,xmm0		/*t37 copy */
		__asm	movaps	xmm3,xmm4		/*t36 copy */
		__asm	mulpd	xmm0,[edi+0x20]	/*t37*c1B */
		__asm	mulpd	xmm4,[edi+0x20]	/*t36*c1B */
		__asm	mulpd	xmm2,[edi+0x30]	/*t37*s1B */
		__asm	mulpd	xmm3,[edi+0x30]	/*t36*s1B */
		__asm	subpd	xmm4,xmm2	/* xmm4 <- im */
		__asm	addpd	xmm0,xmm3	/* xmm0 <- re */

		__asm	movaps	[ebx+0x1b0],xmm4	/* a[jp+p12] */
		__asm	movaps	[ebx+0x1a0],xmm0	/* a[jt+p12] */

	/************************************************/
	/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:	*/
	/************************************************/
		__asm	add	eax, 0x080	/* r38 */

		__asm	movaps	xmm4,[eax+0x020]	/* t2E */				__asm	movaps	xmm0,[eax+0x060]	/* t3E */
		__asm	movaps	xmm5,[eax+0x030]	/* t2F */				__asm	movaps	xmm1,[eax+0x070]	/* t3F */
		__asm	movaps	xmm6,[eax+0x020]	/* xmm2 <- cpy t2E */	__asm	movaps	xmm2,[eax+0x060]	/* xmm6 <- cpy t3E */
		__asm	movaps	xmm7,[eax+0x030]	/* xmm3 <- cpy t2F */	__asm	movaps	xmm3,[eax+0x070]	/* xmm7 <- cpy t3F */

		__asm	mulpd	xmm4,[edx+0x10]	/* t2E*s32_1 */				__asm	mulpd	xmm0,[edx+0x30]	/* t3E*s32_3 */
		__asm	mulpd	xmm5,[edx+0x10]	/* t2F*s32_1 */				__asm	mulpd	xmm1,[edx+0x30]	/* t3F*s32_3 */
		__asm	mulpd	xmm6,[edx     ]	/* t2E*c32_1 */				__asm	mulpd	xmm2,[edx+0x20]	/* t3E*c32_3 */
		__asm	mulpd	xmm7,[edx     ]	/* t2F*c32_1 */				__asm	mulpd	xmm3,[edx+0x20]	/* t3F*c32_3 */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t2F */				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t2E */				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t2F */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t2E */

		__asm	addpd	xmm4,xmm0	/* ~t3E <- t2E+rt */
		__asm	addpd	xmm5,xmm1	/* ~t3F <- t2F+it */
		__asm	subpd	xmm6,xmm0	/* ~t2E <- t2E-rt */
		__asm	subpd	xmm7,xmm1	/* ~t2F <- t2F-it */

		__asm	movaps	xmm2,[eax+0x040]	/* t1E */
		__asm	movaps	xmm3,[eax+0x050]	/* t1F */
		__asm	movaps	xmm0,[eax+0x040]	/* cpy t1E */
		__asm	movaps	xmm1,[eax+0x050]	/* cpy t1F */

		__asm	mulpd	xmm2,[ecx     ]	/* t1E*c */
		__asm	mulpd	xmm1,[ecx+0x10]	/* t1F*s */
		__asm	mulpd	xmm3,[ecx     ]	/* t1F*c */
		__asm	mulpd	xmm0,[ecx+0x10]	/* t1E*s */
		__asm	subpd	xmm2,xmm1	/* xmmE <- rt */
		__asm	addpd	xmm3,xmm0	/* xmmF <- it */

		__asm	movaps	xmm0,[eax      ]	/* t0E */
		__asm	movaps	xmm1,[eax+0x010]	/* t0F */

		__asm	subpd	xmm0,xmm2	/*~t0E <- t0E- rt */
		__asm	subpd	xmm1,xmm3	/*~t0F <- t0F- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t1E <- t0E+ rt */
		__asm	addpd	xmm3,xmm1	/*~t1F <- t0F+ it */
	/*
		rt       =t0E+t2E;		it         =t0F+t2F;
		t0E      =t0E-t2E;		t0F        =t0F-t2F;
		a[jt    ]=rt *c07+it *s07;	a[jp    ]=it *c07-rt *s07;
		a[jt+p10]=t0E*c17+t0F*s17;	a[jp+p10]=t0F*c17-t0E*s17;
	*/
		__asm	add	edi, 0x040	/* c07 */
		__asm	subpd	xmm0,xmm6		/*~t2E <- t0E-t2E */
		__asm	subpd	xmm1,xmm7		/*~t2F <- t0F-t2F */
		__asm	addpd	xmm6,xmm6		/*          2*t2E */
		__asm	addpd	xmm7,xmm7		/*          2*t2F */
		__asm	addpd	xmm6,xmm0		/* rt  <- t0E+t2E */
		__asm	addpd	xmm7,xmm1		/* it  <- t0F+t2F */
		__asm	movaps	[eax+0x020],xmm0	/* store ~t2E */
		__asm	movaps	[eax+0x030],xmm1	/* store ~t2F */
		__asm	movaps	xmm0,xmm6		/* rt copy */
		__asm	movaps	xmm1,xmm7		/* it copy */
		__asm	mulpd	xmm6,[edi     ]	/* rt *c07 */
		__asm	mulpd	xmm7,[edi     ]	/* it *c07 */
		__asm	mulpd	xmm0,[edi+0x10]	/* rt *s07 */
		__asm	mulpd	xmm1,[edi+0x10]	/* it *s07 */
		__asm	subpd	xmm7,xmm0	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm1	/* xmm6 <- re */

		__asm	movaps	[ebx+0x070],xmm7	/* a[jp    ] */
		__asm	movaps	[ebx+0x060],xmm6	/* a[jt    ] */

		__asm	movaps	xmm6,[eax+0x020]	/* load ~t2E */
		__asm	movaps	xmm7,[eax+0x030]	/* load ~t2F */
		__asm	movaps	xmm0,xmm6		/* re copy */
		__asm	movaps	xmm1,xmm7		/* im copy */
		__asm	mulpd	xmm6,[edi+0x20]	/* re *c17 */
		__asm	mulpd	xmm7,[edi+0x20]	/* im *c17 */
		__asm	mulpd	xmm0,[edi+0x30]	/* re *s17 */
		__asm	mulpd	xmm1,[edi+0x30]	/* im *s17 */
		__asm	subpd	xmm7,xmm0	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm1	/* xmm6 <- re */

		__asm	movaps	[ebx+0x170],xmm7	/* a[jp+p8 ] */
		__asm	movaps	[ebx+0x160],xmm6	/* a[jt+p8 ] */
	/*
		rt       =t1E+t3F;		it         =t1F-t3E;
		t3F      =t1E-t3F;		t3E        =t1F+t3E;
		a[jt+p08]=rt *c0F+it *s0F;	a[jp+p08]=it *c0F-rt *s0F;
		a[jt+p18]=t3F*c1F+t3E*s1F;	a[jp+p18]=t3E*c1F-t3F*s1F;
	*/
		__asm	add	edi, 0x040	/* c0F */
		__asm	subpd	xmm2,xmm5		/*~t3F <- t1E-t3F */
		__asm	subpd	xmm3,xmm4		/* it  <- t1F-t3E */
		__asm	addpd	xmm5,xmm5		/*          2*t3F */
		__asm	addpd	xmm4,xmm4		/*          2*t3E */
		__asm	addpd	xmm5,xmm2		/* rt  <- t1E+t3F */
		__asm	addpd	xmm4,xmm3		/*~t3E <- t1F+t3E */
		__asm	movaps	xmm0,xmm5		/* rt copy */
		__asm	movaps	xmm1,xmm3		/* it copy */
		__asm	mulpd	xmm5,[edi     ]	/* rt*c0F */
		__asm	mulpd	xmm3,[edi     ]	/* it*c0F */
		__asm	mulpd	xmm0,[edi+0x10]	/* rt*s0F */
		__asm	mulpd	xmm1,[edi+0x10]	/* it*s0F */
		__asm	subpd	xmm3,xmm0	/* xmm3 <- im */
		__asm	addpd	xmm5,xmm1	/* xmm5 <- re */

		__asm	movaps	[ebx+0x0f0],xmm3	/* a[jp+p4 ] */
		__asm	movaps	[ebx+0x0e0],xmm5	/* a[jt+p4 ] */

		__asm	movaps	xmm0,xmm2		/*t3F copy */
		__asm	movaps	xmm1,xmm4		/*t3E copy */
		__asm	mulpd	xmm2,[edi+0x20]	/*t3F*c1F */
		__asm	mulpd	xmm4,[edi+0x20]	/*t3E*c1F */
		__asm	mulpd	xmm0,[edi+0x30]	/*t3F*s1F */
		__asm	mulpd	xmm1,[edi+0x30]	/*t3E*s1F */
		__asm	subpd	xmm4,xmm0	/* xmm4 <- im */
		__asm	addpd	xmm2,xmm1	/* xmm2 <- re */

		__asm	movaps	[ebx+0x1f0],xmm4	/* a[jp+p12] */
		__asm	movaps	[ebx+0x1e0],xmm2	/* a[jt+p12] */

	/************************************************/
	/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:	*/
	/************************************************/
		__asm	mov	edx, r00

		__asm	movaps	xmm0,[edx      ]	/* t00 */
		__asm	movaps	xmm1,[edx+0x010]	/* t01 */
		__asm	movaps	xmm2,[edx+0x040]	/* t10 */
		__asm	movaps	xmm3,[edx+0x050]	/* t11 */

		__asm	subpd	xmm0,[edx+0x040]	/*~t10=t00-t10*/
		__asm	subpd	xmm1,[edx+0x050]	/*~t11=t01-t11*/
		__asm	addpd	xmm2,[edx      ]	/*~t00=t10+t00*/
		__asm	addpd	xmm3,[edx+0x010]	/*~t01=t11+t01*/

		__asm	movaps	xmm4,[edx+0x020]	/* t20 */
		__asm	movaps	xmm5,[edx+0x030]	/* t21 */
		__asm	movaps	xmm6,[edx+0x060]	/* t30 */
		__asm	movaps	xmm7,[edx+0x070]	/* t31 */

		__asm	subpd	xmm4,[edx+0x060]	/*~t30=t20-t30*/
		__asm	subpd	xmm5,[edx+0x070]	/*~t31=t21-t31*/
		__asm	addpd	xmm6,[edx+0x020]	/*~t20=t30+t20*/
		__asm	addpd	xmm7,[edx+0x030]	/*~t21=t31+t21*/

		__asm	mov	eax, add0	/* restore main-array index */
		__asm	mov	ecx, c10
	/*
		a[jt    ]=t00+t20;			a[jp    ]=t01+t21;
		t20      =t00-t20;			t21      =t01-t21;
		a[jt+p10]=t20*c10+t21*s10;	a[jp+p10]=t21*c10-t20*s10;
	*/
		__asm	addpd	xmm2,xmm6		/* t1 +t17 */
		__asm	addpd	xmm3,xmm7		/* t2 +t18 */
		__asm	movaps	[edx      ],xmm2	/* tmp store non-interleaved a[jt+p0 ] */
		__asm	movaps	[edx+0x010],xmm3	/* tmp store non-interleaved a[jp+p0 ] */
		__asm	addpd	xmm6,xmm6		/*   2*t17 */
		__asm	addpd	xmm7,xmm7		/*   2*t18 */
		__asm	subpd	xmm2,xmm6		/*~t17 <- t1 -t17 */
		__asm	subpd	xmm3,xmm7		/*~t18 <- t2 -t18 */
		__asm	movaps	xmm6,xmm2		/*~t17 copy */
		__asm	movaps	xmm7,xmm3		/*~t18 copy */
		__asm	mulpd	xmm2,[ecx     ]	/*~t17*c10*/
		__asm	mulpd	xmm3,[ecx     ]	/*~t18*c10*/
		__asm	mulpd	xmm6,[ecx+0x10]	/*~t17*s10*/
		__asm	mulpd	xmm7,[ecx+0x10]	/*~t18*s10*/
		__asm	subpd	xmm3,xmm6	/* xmm3 <- it */
		__asm	addpd	xmm2,xmm7	/* xmm2 <- rt 	xmm6,7 free */

		__asm	mov	ebx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm7,xmm3	/* cpy a[jp+p10] */
		__asm	movaps		xmm6,xmm2	/* cpy a[jt+p10] */
		__asm	unpckhpd	xmm7,[ebx+0x110]
		__asm	unpcklpd	xmm3,[ebx+0x110]
		__asm	movaps	[ebx+0x110],xmm7	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm6,[ebx+0x100]
		__asm	unpcklpd	xmm2,[ebx+0x100]
		__asm	movaps	[ebx+0x100],xmm6	/* Store hi real in aj2 */

		/* lo halves go into aj1: */
		__asm	movaps	[eax+0x110],xmm3	/* a[jp+p10] */
		__asm	movaps	[eax+0x100],xmm2	/* a[jt+p10] */

		__asm	movaps	xmm3,[edx+0x010]	/* reload a[jp+p0 ] */
		__asm	movaps	xmm2,[edx      ]	/* reload a[jt+p0 ] */
		__asm	movaps		xmm7,xmm3	/* cpy a[jp    ] */
		__asm	movaps		xmm6,xmm2	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm7,[ebx+0x010]
		__asm	unpcklpd	xmm3,[ebx+0x010]
		__asm	movaps	[ebx+0x10],xmm7	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm6,[ebx      ]
		__asm	unpcklpd	xmm2,[ebx      ]
		__asm	movaps	[ebx     ],xmm6	/* Store hi real in aj2 */

		/* lo halves go into aj1: */
		__asm	movaps	[eax+0x010],xmm3	/* a[jp+p0 ] */
		__asm	movaps	[eax      ],xmm2	/* a[jt+p0 ] */

	/* mpy by E^-4 = -I is inlined here...
		rt       =t10+t31;			it         =t11-t30;
		t31      =t10-t31;			t30        =t11+t30;
		a[jt+p08]=rt *c08+it *s08;	a[jp+p08]=it *c08-rt *s08;
		a[jt+p18]=t31*c18+t30*s18;	a[jp+p18]=t30*c18-t31*s18;
	*/
		__asm	mov	ecx, c08
		__asm	addpd	xmm0,xmm5		/* rt  <- t9 +t26 */
		__asm	subpd	xmm1,xmm4		/* it  <- t10-t25 */
		__asm	movaps	xmm2,xmm0		/* rt  copy */
		__asm	movaps	xmm3,xmm1		/* it  copy */
		__asm	addpd	xmm5,xmm5		/*          2*t26 */
		__asm	addpd	xmm4,xmm4		/*          2*t25 */
		__asm	movaps	xmm6,xmm0		/* rt  copy */
		__asm	movaps	xmm7,xmm1		/* it  copy */
		__asm	mulpd	xmm2,[ecx     ]	/* rt *c4 */
		__asm	mulpd	xmm3,[ecx     ]	/* it *c4 */
		__asm	mulpd	xmm6,[ecx+0x10]	/* rt *s4 */
		__asm	mulpd	xmm7,[ecx+0x10]	/* it *s4 */
		__asm	subpd	xmm3,xmm6	/* xmm3 <- im */
		__asm	addpd	xmm2,xmm7	/* xmm2 <- re */

		__asm	movaps		xmm7,xmm3	/* cpy a[jp+p08] */
		__asm	movaps		xmm6,xmm2	/* cpy a[jt+p08] */
		__asm	unpckhpd	xmm7,[ebx+0x090]
		__asm	unpcklpd	xmm3,[ebx+0x090]
		__asm	movaps	[ebx+0x090],xmm7	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm6,[ebx+0x080]
		__asm	unpcklpd	xmm2,[ebx+0x080]
		__asm	movaps	[ebx+0x080],xmm6	/* Store hi real in aj2 */
		/* lo halves go into aj1: */
		__asm	movaps	[eax+0x090],xmm3	/* a[jp+p08] */
		__asm	movaps	[eax+0x080],xmm2	/* a[jt+p08] */

		__asm	mov	ecx, c18
		__asm	subpd	xmm0,xmm5		/*~t26 <- t9 -t26 */
		__asm	addpd	xmm1,xmm4		/*~t25 <- t10+t25 */
		__asm	movaps	xmm6,xmm0		/*~t26 copy */
		__asm	movaps	xmm7,xmm1		/*~t25 copy */
		__asm	mulpd	xmm0,[ecx     ]	/*~t26*c18*/
		__asm	mulpd	xmm1,[ecx     ]	/*~t25*c18*/
		__asm	mulpd	xmm6,[ecx+0x10]	/*~t26*s18*/
		__asm	mulpd	xmm7,[ecx+0x10]	/*~t25*s18*/
		__asm	subpd	xmm1,xmm6	/* xmm3 <- im */
		__asm	addpd	xmm0,xmm7	/* xmm2 <- re */

		__asm	movaps		xmm7,xmm1	/* cpy a[jp+p08] */
		__asm	movaps		xmm6,xmm0	/* cpy a[jt+p08] */
		__asm	unpckhpd	xmm7,[ebx+0x190]
		__asm	unpcklpd	xmm1,[ebx+0x190]
		__asm	movaps	[ebx+0x190],xmm7	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm6,[ebx+0x180]
		__asm	unpcklpd	xmm0,[ebx+0x180]
		__asm	movaps	[ebx+0x180],xmm6	/* Store hi real in aj2 */
		/* lo halves go into aj1: */
		__asm	movaps	[eax+0x190],xmm1	/* a[jp+p18] */
		__asm	movaps	[eax+0x180],xmm0	/* a[jt+p18] */

	/************************************************/
	/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E	*/
	/************************************************/
		__asm	mov	edx, r08
		__asm	movaps	xmm2,[esi]	/* isrt2 */

		__asm	movaps	xmm4,[edx+0x020]	/* t28 */
		__asm	movaps	xmm5,[edx+0x030]	/* t29 */
		__asm	movaps	xmm0,[edx+0x060]	/* t38 */
		__asm	movaps	xmm1,[edx+0x070]	/* t39 */

		__asm	addpd	xmm4,[edx+0x030]	/*~t28=t28+t29*/
		__asm	subpd	xmm5,[edx+0x020]	/*~t29=t29-t28*/
		__asm	subpd	xmm0,[edx+0x070]	/* rt =t38-t39*/
		__asm	addpd	xmm1,[edx+0x060]	/* it =t39+t38*/
		__asm	mulpd	xmm4,xmm2
		__asm	mulpd	xmm5,xmm2
		__asm	mulpd	xmm0,xmm2
		__asm	mulpd	xmm1,xmm2
		__asm	movaps	xmm6,xmm4			/* t28 copy */
		__asm	movaps	xmm7,xmm5			/* t29 copy */

		__asm	subpd	xmm4,xmm0			/*~t28=t28-rt */
		__asm	subpd	xmm5,xmm1			/*~t29=t29-it */
		__asm	addpd	xmm6,xmm0			/*~t38=t28+rt */
		__asm	addpd	xmm7,xmm1			/*~t39=t29+it */

		__asm	movaps	xmm0,[edx      ]	/* t08 */
		__asm	movaps	xmm1,[edx+0x010]	/* t09 */
		__asm	movaps	xmm2,[edx+0x040]	/* t18 */
		__asm	movaps	xmm3,[edx+0x050]	/* t19 */

		__asm	subpd	xmm0,[edx+0x050]	/*~t18=t08-t19*/
		__asm	subpd	xmm1,[edx+0x040]	/*~t09=t09-t18*/
		__asm	addpd	xmm3,[edx      ]	/*~t08=t19+t08*/
		__asm	addpd	xmm2,[edx+0x010]	/*~t19=t18+t09*/

		__asm	mov	eax, add0	/* restore main-array index, but keep r08 in edx for temporary storage */
	/*
		rt       =t08+t28;		it         =t09+t29;
		t28      =t08-t28;		t29        =t09-t29;
		a[jt    ]=rt *c04+it *s04;	a[jp    ]=it *c04-rt *s04;
		a[jt+p10]=t28*c14+t29*s14;	a[jp+p10]=t29*c14-t28*s14;
	*/
		__asm	mov	ecx, c04
		__asm	subpd	xmm3,xmm4		/*~t28 <- t08-t28 */
		__asm	subpd	xmm1,xmm5		/*~t29 <- t09-t29 */
		__asm	addpd	xmm4,xmm4		/*          2*t28 */
		__asm	addpd	xmm5,xmm5		/*          2*t29 */
		__asm	addpd	xmm4,xmm3		/* rt  <- t08+t28 */
		__asm	addpd	xmm5,xmm1		/* it  <- t09+t29 */
		__asm	movaps	[edx      ],xmm3	/* tmp store ~t28 */
		__asm	movaps	[edx+0x010],xmm1	/* tmp store ~t29 */
		__asm	movaps	xmm3,xmm4		/* rt copy */
		__asm	movaps	xmm1,xmm5		/* it copy */
		__asm	mulpd	xmm4,[ecx     ]	/* rt*c04 */
		__asm	mulpd	xmm5,[ecx     ]	/* it*c04 */
		__asm	mulpd	xmm3,[ecx+0x10]	/* rt*s04 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* it*s04 */
		__asm	subpd	xmm5,xmm3	/* xmm5 <- im */
		__asm	addpd	xmm4,xmm1	/* xmm4 <- re	xmm1,3 free */

		__asm	mov	ebx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm3,xmm5	/* cpy a[jp+p10] */
		__asm	movaps		xmm1,xmm4	/* cpy a[jt+p10] */
		__asm	unpckhpd	xmm3,[ebx+0x050]
		__asm	unpcklpd	xmm5,[ebx+0x050]
		__asm	movaps	[ebx+0x050],xmm3	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm1,[ebx+0x040]
		__asm	unpcklpd	xmm4,[ebx+0x040]
		__asm	movaps	[ebx+0x040],xmm1	/* Store hi real in aj2 */

		/* lo halves go into aj1: */
		__asm	movaps	[eax+0x050],xmm5	/* a[jp    ] */
		__asm	movaps	[eax+0x040],xmm4	/* a[jt    ] */

		__asm	mov	ecx, c14
		__asm	movaps	xmm4,[edx      ]	/* reload ~t28 */
		__asm	movaps	xmm5,[edx+0x010]	/* reload ~t29 */
		__asm	movaps	xmm3,xmm4		/* re copy */
		__asm	movaps	xmm1,xmm5		/* im copy */
		__asm	mulpd	xmm4,[ecx     ]	/* re*c14 */
		__asm	mulpd	xmm5,[ecx     ]	/* im*c14 */
		__asm	mulpd	xmm3,[ecx+0x10]	/* re*s14 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* im*s14 */
		__asm	subpd	xmm5,xmm3	/* xmm5 <- im */
		__asm	addpd	xmm4,xmm1	/* xmm4 <- re 	xmm1,3 free */

		__asm	movaps		xmm3,xmm5	/* cpy a[jp+p8 ] */
		__asm	movaps		xmm1,xmm4	/* cpy a[jt+p8 ] */
		__asm	unpckhpd	xmm3,[ebx+0x150]
		__asm	unpcklpd	xmm5,[ebx+0x150]
		__asm	movaps	[ebx+0x150],xmm3	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm1,[ebx+0x140]
		__asm	unpcklpd	xmm4,[ebx+0x140]
		__asm	movaps	[ebx+0x140],xmm1	/* Store hi real in aj2 */

		/* lo halves go into aj1: */
		__asm	movaps	[eax+0x150],xmm5	/* a[jp+p8 ] */
		__asm	movaps	[eax+0x140],xmm4	/* a[jt+p8 ] */

	/* mpy by E^-4 = -I is inlined here...
		rt       =t18+t39;		it         =t19-t38;
		t39      =t18-t39;		t38        =t19+t38;
		a[jt+p08]=rt *c0C+it *s0C;	a[jp+p08]=it *c0C-rt *s0C;
		a[jt+p18]=t39*c1C+t38*s1C;	a[jp+p18]=t38*c1C-t39*s1C;
	*/
		__asm	mov	ecx, c0C
		__asm	subpd	xmm0,xmm7		/*~t39 <- t18-t39 */
		__asm	subpd	xmm2,xmm6		/* it  <- t19-t38 */
		__asm	addpd	xmm7,xmm7		/*          2*t39 */
		__asm	addpd	xmm6,xmm6		/*          2*t38 */
		__asm	addpd	xmm7,xmm0		/* rt  <- t18+t39 */
		__asm	addpd	xmm6,xmm2		/*~t38 <- t19+t38 */
		__asm	movaps	xmm4,xmm7		/* rt copy */
		__asm	movaps	xmm5,xmm2		/* it copy */
		__asm	mulpd	xmm7,[ecx     ]	/* rt*c0C */
		__asm	mulpd	xmm2,[ecx     ]	/* it*c0C */
		__asm	mulpd	xmm4,[ecx+0x10]	/* rt*s0C */
		__asm	mulpd	xmm5,[ecx+0x10]	/* it*s0C */
		__asm	subpd	xmm2,xmm4	/* xmm2 <- im */
		__asm	addpd	xmm7,xmm5	/* xmm7 <- re	xmm4,5 free */

		__asm	movaps		xmm5,xmm2	/* cpy a[jp    ] */
		__asm	movaps		xmm4,xmm7	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm5,[ebx+0x0d0]
		__asm	unpcklpd	xmm2,[ebx+0x0d0]
		__asm	movaps	[ebx+0x0d0],xmm5	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm4,[ebx+0x0c0]
		__asm	unpcklpd	xmm7,[ebx+0x0c0]
		__asm	movaps	[ebx+0x0c0],xmm4	/* Store hi real in aj2 */

		/* lo halves go into aj1: */
		__asm	movaps	[eax+0x0d0],xmm2	/* a[jp+p4 ] */
		__asm	movaps	[eax+0x0c0],xmm7	/* a[jt+p4 ] */

		__asm	mov	ecx, c1C
		__asm	movaps	xmm4,xmm0		/* t39 copy */
		__asm	movaps	xmm5,xmm6		/* t38 copy */
		__asm	mulpd	xmm0,[ecx     ]	/* t39*c1C */
		__asm	mulpd	xmm6,[ecx     ]	/* t38*c1C */
		__asm	mulpd	xmm4,[ecx+0x10]	/* t39*s1C */
		__asm	mulpd	xmm5,[ecx+0x10]	/* t38*s1C */
		__asm	subpd	xmm6,xmm4	/* xmm6 <- im */
		__asm	addpd	xmm0,xmm5	/* xmm7 <- re	xmm4,5 free */

		__asm	movaps		xmm5,xmm6	/* cpy a[jp    ] */
		__asm	movaps		xmm4,xmm0	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm5,[ebx+0x1d0]
		__asm	unpcklpd	xmm6,[ebx+0x1d0]
		__asm	movaps	[ebx+0x1d0],xmm5	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm4,[ebx+0x1c0]
		__asm	unpcklpd	xmm0,[ebx+0x1c0]
		__asm	movaps	[ebx+0x1c0],xmm4	/* Store hi real in aj2 */

		/* lo halves go into aj1: */
		__asm	movaps	[eax+0x1d0],xmm6	/* a[jp+p12] */
		__asm	movaps	[eax+0x1c0],xmm0	/* a[jt+p12] */

	/************************************************/
	/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26	*/
	/************************************************/
		__asm	mov	edx, r20
		__asm	mov	ecx, esi
		__asm	add	ecx, 0x10	/* cc0 */

		__asm	movaps	xmm4,[edx+0x020]	/* t24 */				__asm	movaps	xmm0,[edx+0x060]	/* t34 */
		__asm	movaps	xmm5,[edx+0x030]	/* t25 */				__asm	movaps	xmm1,[edx+0x070]	/* t35 */
		__asm	movaps	xmm6,[edx+0x020]	/* xmm2 <- cpy t24 */	__asm	movaps	xmm2,[edx+0x060]	/* xmm6 <- cpy t34 */
		__asm	movaps	xmm7,[edx+0x030]	/* xmm3 <- cpy t25 */	__asm	movaps	xmm3,[edx+0x070]	/* xmm7 <- cpy t35 */

		__asm	mulpd	xmm4,[ecx     ]	/* t24*c */					__asm	mulpd	xmm0,[ecx+0x10]	/* t34*s */
		__asm	mulpd	xmm5,[ecx     ]	/* t25*c */					__asm	mulpd	xmm1,[ecx+0x10]	/* t35*s */
		__asm	mulpd	xmm6,[ecx+0x10]	/* t24*s */					__asm	mulpd	xmm2,[ecx     ]	/* t34*c */
		__asm	mulpd	xmm7,[ecx+0x10]	/* t25*s */					__asm	mulpd	xmm3,[ecx     ]	/* t35*c */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t25*/				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t24*/				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t25*/
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t24*/

		__asm	addpd	xmm4,xmm0	/* ~t24 <- t24+rt */
		__asm	addpd	xmm5,xmm1	/* ~t25 <- t25+it */
		__asm	subpd	xmm6,xmm0	/* ~t34 <- t24-rt */
		__asm	subpd	xmm7,xmm1	/* ~t35 <- t25-it */

		__asm	movaps	xmm2,[edx+0x040]	/* t14 */
		__asm	movaps	xmm3,[edx+0x050]	/* t15 */
		__asm	movaps	xmm0,[edx      ]	/* t04 */
		__asm	movaps	xmm1,[edx+0x010]	/* t05 */
		__asm	addpd	xmm2,[edx+0x050]	/*~t14=t14+t15*/
		__asm	subpd	xmm3,[edx+0x040]	/*~t15=t15-t14*/
		__asm	mulpd	xmm2,[esi]	/* rt */
		__asm	mulpd	xmm3,[esi]	/* it */

		__asm	subpd	xmm0,xmm2	/*~t14 <- t04- rt */
		__asm	subpd	xmm1,xmm3	/*~t15 <- t05- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t04 <- t04+ rt */
		__asm	addpd	xmm3,xmm1	/*~t05 <- t05+ it */

		__asm	mov	eax, add0	/* restore main-array index0, but keep r02 in edx for temporary storage */
	/*
		rt       =t04+t24;		it         =t05+t25;
		t24      =t04-t24;		t25        =t05-t25;
		a[jt    ]=rt *c02+it *s02;	a[jp    ]=it *c02-rt *s02;
		a[jt+p10]=t24*c12+t25*s12;	a[jp+p10]=t25*c12-t24*s12;
	*/
		__asm	mov	ecx, c02
		__asm	subpd	xmm2,xmm4		/*~t24 <- t04-t24 */
		__asm	subpd	xmm3,xmm5		/*~t25 <- t05-t25 */
		__asm	addpd	xmm4,xmm4		/*          2*t24 */
		__asm	addpd	xmm5,xmm5		/*          2*t25 */
		__asm	addpd	xmm4,xmm2		/* rt  <- t04+t24 */
		__asm	addpd	xmm5,xmm3		/* it  <- t05+t25 */
		__asm	movaps	[edx      ],xmm2	/* tmp store ~t24 */
		__asm	movaps	[edx+0x010],xmm3	/* tmp store ~t25 */
		__asm	movaps	xmm2,xmm4		/* rt copy */
		__asm	movaps	xmm3,xmm5		/* it copy */
		__asm	mulpd	xmm4,[ecx     ]	/* rt *c02 */
		__asm	mulpd	xmm5,[ecx     ]	/* it *c02 */
		__asm	mulpd	xmm2,[ecx+0x10]	/* rt *s02 */
		__asm	mulpd	xmm3,[ecx+0x10]	/* it *s02 */
		__asm	subpd	xmm5,xmm2	/* xmm5 <- im */
		__asm	addpd	xmm4,xmm3	/* xmm4 <- re	xmm2,3 free */

		__asm	mov	ebx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm3,xmm5	/* cpy a[jp    ] */
		__asm	movaps		xmm2,xmm4	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm3,[ebx+0x030]
		__asm	unpcklpd	xmm5,[ebx+0x030]
		__asm	movaps	[ebx+0x030],xmm3	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm2,[ebx+0x020]
		__asm	unpcklpd	xmm4,[ebx+0x020]
		__asm	movaps	[ebx+0x020],xmm2	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x030],xmm5	/* a[jp    ] */
		__asm	movaps	[eax+0x020],xmm4	/* a[jt    ] */

		__asm	mov	ecx, c12
		__asm	movaps	xmm4,[edx      ]	/* load ~t24 */
		__asm	movaps	xmm5,[edx+0x010]	/* load ~t25 */
		__asm	movaps	xmm2,xmm4		/* re copy */
		__asm	movaps	xmm3,xmm5		/* im copy */
		__asm	mulpd	xmm4,[ecx     ]	/* re *c12 */
		__asm	mulpd	xmm5,[ecx     ]	/* im *c12 */
		__asm	mulpd	xmm2,[ecx+0x10]	/* re *s12 */
		__asm	mulpd	xmm3,[ecx+0x10]	/* im *s12 */
		__asm	subpd	xmm5,xmm2	/* xmm5 <- im */
		__asm	addpd	xmm4,xmm3	/* xmm4 <- re	xmm2,3 free */

		__asm	movaps		xmm3,xmm5	/* cpy a[jp    ] */
		__asm	movaps		xmm2,xmm4	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm3,[ebx+0x130]
		__asm	unpcklpd	xmm5,[ebx+0x130]
		__asm	movaps	[ebx+0x130],xmm3	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm2,[ebx+0x120]
		__asm	unpcklpd	xmm4,[ebx+0x120]
		__asm	movaps	[ebx+0x120],xmm2	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x130],xmm5	/* a[jp+p8 ] */
		__asm	movaps	[eax+0x120],xmm4	/* a[jt+p8 ] */

	/* mpy by E^-4 = -I is inlined here...
		rt       =t14+t35;		it         =t15-t34;
		t35      =t14-t35;		t34        =t15+t34;
		a[jt+p08]=rt *c0A+it *s0A;	a[jp+p08]=it *c0A-rt *s0A;
		a[jt+p18]=t35*c1A+t34*s1A;	a[jp+p18]=t34*c1A-t35*s1A;
	*/
		__asm	mov	ecx, c0A
		__asm	subpd	xmm0,xmm7		/*~t35 <- t14-t35 */
		__asm	subpd	xmm1,xmm6		/* it  <- t15-t34 */
		__asm	addpd	xmm7,xmm7		/*          2*t35 */
		__asm	addpd	xmm6,xmm6		/*          2*t34 */
		__asm	addpd	xmm7,xmm0		/* rt  <- t14+t35 */
		__asm	addpd	xmm6,xmm1		/*~t34 <- t15+t34 */
		__asm	movaps	xmm4,xmm7		/* rt copy */
		__asm	movaps	xmm5,xmm1		/* it copy */
		__asm	mulpd	xmm7,[ecx     ]	/* rt*c0A */
		__asm	mulpd	xmm1,[ecx     ]	/* it*c0A */
		__asm	mulpd	xmm4,[ecx+0x10]	/* rt*s0A */
		__asm	mulpd	xmm5,[ecx+0x10]	/* it*s0A */
		__asm	subpd	xmm1,xmm4	/* xmm1 <- im */
		__asm	addpd	xmm7,xmm5	/* xmm7 <- re	xmm4,5 free */

		__asm	movaps		xmm5,xmm1	/* cpy a[jp    ] */
		__asm	movaps		xmm4,xmm7	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm5,[ebx+0x0b0]
		__asm	unpcklpd	xmm1,[ebx+0x0b0]
		__asm	movaps	[ebx+0x0b0],xmm5	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm4,[ebx+0x0a0]
		__asm	unpcklpd	xmm7,[ebx+0x0a0]
		__asm	movaps	[ebx+0x0a0],xmm4	/* Store hi real in aj2 */
		__asm	movaps	[eax+0x0b0],xmm1	/* a[jp+p4 ] */
		__asm	movaps	[eax+0x0a0],xmm7	/* a[jt+p4 ] */

		__asm	mov	ecx, c1A
		__asm	movaps	xmm4,xmm0		/*t35 copy */
		__asm	movaps	xmm5,xmm6		/*t34 copy */
		__asm	mulpd	xmm0,[ecx     ]	/*t35*c1A */
		__asm	mulpd	xmm6,[ecx     ]	/*t34*c1A */
		__asm	mulpd	xmm4,[ecx+0x10]	/*t35*s1A */
		__asm	mulpd	xmm5,[ecx+0x10]	/*t34*s1A */
		__asm	subpd	xmm6,xmm4	/* xmm6 <- im */
		__asm	addpd	xmm0,xmm5	/* xmm7 <- re	xmm4,5 free */

		__asm	movaps		xmm5,xmm6	/* cpy a[jp    ] */
		__asm	movaps		xmm4,xmm0	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm5,[ebx+0x1b0]
		__asm	unpcklpd	xmm6,[ebx+0x1b0]
		__asm	movaps	[ebx+0x1b0],xmm5	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm4,[ebx+0x1a0]
		__asm	unpcklpd	xmm0,[ebx+0x1a0]
		__asm	movaps	[ebx+0x1a0],xmm4	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x1b0],xmm6	/* a[jp+p12] */
		__asm	movaps	[eax+0x1a0],xmm0	/* a[jt+p12] */

	/************************************************/
	/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E	*/
	/************************************************/
		__asm	mov	edx, r28
		__asm	mov	ecx, esi
		__asm	add	ecx, 0x10	/* cc0 */

		__asm	movaps	xmm4,[edx+0x020]	/* t2C */				__asm	movaps	xmm0,[edx+0x060]	/* t3C */
		__asm	movaps	xmm5,[edx+0x030]	/* t2D */				__asm	movaps	xmm1,[edx+0x070]	/* t3D */
		__asm	movaps	xmm6,[edx+0x020]	/* xmm2 <- cpy t2C */	__asm	movaps	xmm2,[edx+0x060]	/* xmm6 <- cpy t3C */
		__asm	movaps	xmm7,[edx+0x030]	/* xmm3 <- cpy t2D */	__asm	movaps	xmm3,[edx+0x070]	/* xmm7 <- cpy t3D */

		__asm	mulpd	xmm4,[ecx+0x10]	/* t2C*s */					__asm	mulpd	xmm0,[ecx     ]	/* t3C*c */
		__asm	mulpd	xmm5,[ecx+0x10]	/* t2D*s */					__asm	mulpd	xmm1,[ecx     ]	/* t3D*c */
		__asm	mulpd	xmm6,[ecx     ]	/* t2C*c */					__asm	mulpd	xmm2,[ecx+0x10]	/* t3C*s */
		__asm	mulpd	xmm7,[ecx     ]	/* t2D*c */					__asm	mulpd	xmm3,[ecx+0x10]	/* t3D*s */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t2D */				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t2C */				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t2D */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t2C */

		__asm	addpd	xmm4,xmm0	/* ~t3C <- t2C+rt */
		__asm	addpd	xmm5,xmm1	/* ~t3D <- t2D+it */
		__asm	subpd	xmm6,xmm0	/* ~t2C <- t2C-rt */
		__asm	subpd	xmm7,xmm1	/* ~t2D <- t2D-it */

		__asm	movaps	xmm2,[edx+0x040]	/* t1C */
		__asm	movaps	xmm3,[edx+0x050]	/* t1D */
		__asm	movaps	xmm0,[edx      ]	/* t0C */
		__asm	movaps	xmm1,[edx+0x010]	/* t0D */
		__asm	subpd	xmm2,[edx+0x050]	/*~t1C=t1C-t1D */
		__asm	addpd	xmm3,[edx+0x040]	/*~t1D=t1D+t1C */
		__asm	mulpd	xmm2,[esi]	/* rt */
		__asm	mulpd	xmm3,[esi]	/* it */

		__asm	subpd	xmm0,xmm2	/*~t0C <- t0C- rt */
		__asm	subpd	xmm1,xmm3	/*~t0D <- t0D- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t1C <- t0C+ rt */
		__asm	addpd	xmm3,xmm1	/*~t1D <- t0D+ it */

		__asm	mov	eax, add0	/* restore main-array index0, but keep r28 in edx for temporary storage */
	/*
		rt       =t0C+t2C;			it       =t0D+t2D;
		t2C      =t0C-t2C;			t2D      =t0D-t2D;
		a[jt    ]=rt *c06+it *s06	a[jp    ]=it *c06-rt *s06
		a[jt+p8 ]=t2C*c16+t2D*s16;	a[jp+p8 ]=t2C*c16-t2D*s16;
	*/
		__asm	mov	ecx, c06
		__asm	subpd	xmm0,xmm6		/*~t2C <- t0C-t2C */
		__asm	subpd	xmm1,xmm7		/*~t2D <- t0D-t2D */
		__asm	addpd	xmm6,xmm6		/*          2*t2C */
		__asm	addpd	xmm7,xmm7		/*          2*t2D */
		__asm	addpd	xmm6,xmm0		/* rt  <- t0C+t2C */
		__asm	addpd	xmm7,xmm1		/* it  <- t0D+t2D */
		__asm	movaps	[edx      ],xmm0	/* tmp store ~t2C */
		__asm	movaps	[edx+0x010],xmm1	/* tmp store ~t2D */
		__asm	movaps	xmm0,xmm6		/* rt copy */
		__asm	movaps	xmm1,xmm7		/* it copy */
		__asm	mulpd	xmm6,[ecx     ]	/* rt *c06 */
		__asm	mulpd	xmm7,[ecx     ]	/* it *c06 */
		__asm	mulpd	xmm0,[ecx+0x10]	/* rt *s06 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* it *s06 */
		__asm	subpd	xmm7,xmm0	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm1	/* xmm6 <- re	xmm0,1 free */

		__asm	mov	ebx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm1,xmm7	/* cpy a[jp    ] */
		__asm	movaps		xmm0,xmm6	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm1,[ebx+0x070]
		__asm	unpcklpd	xmm7,[ebx+0x070]
		__asm	movaps	[ebx+0x070],xmm1	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm0,[ebx+0x060]
		__asm	unpcklpd	xmm6,[ebx+0x060]
		__asm	movaps	[ebx+0x060],xmm0	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x070],xmm7	/* a[jp    ] */
		__asm	movaps	[eax+0x060],xmm6	/* a[jt    ] */

		__asm	mov	ecx, c16
		__asm	movaps	xmm6,[edx      ]	/* load ~t2C */
		__asm	movaps	xmm7,[edx+0x010]	/* load ~t2D */
		__asm	movaps	xmm0,xmm6		/* re copy */
		__asm	movaps	xmm1,xmm7		/* im copy */
		__asm	mulpd	xmm6,[ecx     ]	/* re *c16 */
		__asm	mulpd	xmm7,[ecx     ]	/* im *c16 */
		__asm	mulpd	xmm0,[ecx+0x10]	/* re *s16 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* im *s16 */
		__asm	subpd	xmm7,xmm0	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm1	/* xmm6 <- re	xmm0,1 free */

		__asm	movaps		xmm1,xmm7	/* cpy a[jp    ] */
		__asm	movaps		xmm0,xmm6	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm1,[ebx+0x170]
		__asm	unpcklpd	xmm7,[ebx+0x170]
		__asm	movaps	[ebx+0x170],xmm1	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm0,[ebx+0x160]
		__asm	unpcklpd	xmm6,[ebx+0x160]
		__asm	movaps	[ebx+0x160],xmm0	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x170],xmm7	/* a[jp+p8 ] */
		__asm	movaps	[eax+0x160],xmm6	/* a[jt+p8 ] */

	/* mpy by E^-4 = -I is inlined here...
		rt       =t1C+t3D;		it         =t1D-t3C;
		t3D      =t1C-t3D;		t3C        =t1D+t3C;
		a[jt+p08]=rt *c0E+it *s0E;	a[jp+p08]=it *c0E-rt *s0E;
		a[jt+p18]=t3D*c1E+t3C*s1E;	a[jp+p18]=t3C*c1E-t3D*s1E;
	*/
		__asm	mov	ecx, c0E
		__asm	subpd	xmm2,xmm5		/*~t3D <- t1C-t3D */
		__asm	subpd	xmm3,xmm4		/* it  <- t1D-t3C */
		__asm	addpd	xmm5,xmm5		/*          2*t3D */
		__asm	addpd	xmm4,xmm4		/*          2*t3C */
		__asm	addpd	xmm5,xmm2		/* rt  <- t1C+t3D */
		__asm	addpd	xmm4,xmm3		/*~t3C <- t1D+t3C */
		__asm	movaps	xmm0,xmm5		/* rt copy */
		__asm	movaps	xmm1,xmm3		/* it copy */
		__asm	mulpd	xmm5,[ecx     ]	/* rt*c0E */
		__asm	mulpd	xmm3,[ecx     ]	/* it*c0E */
		__asm	mulpd	xmm0,[ecx+0x10]	/* rt*s0E */
		__asm	mulpd	xmm1,[ecx+0x10]	/* it*s0E */
		__asm	subpd	xmm3,xmm0	/* xmm3 <- im */
		__asm	addpd	xmm5,xmm1	/* xmm5 <- re	xmm0,1 free */

		__asm	movaps		xmm1,xmm3	/* cpy a[jp    ] */
		__asm	movaps		xmm0,xmm5	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm1,[ebx+0x0f0]
		__asm	unpcklpd	xmm3,[ebx+0x0f0]
		__asm	movaps	[ebx+0x0f0],xmm1	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm0,[ebx+0x0e0]
		__asm	unpcklpd	xmm5,[ebx+0x0e0]
		__asm	movaps	[ebx+0x0e0],xmm0	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x0f0],xmm3	/* a[jp+p4 ] */
		__asm	movaps	[eax+0x0e0],xmm5	/* a[jt+p4 ] */

		__asm	mov	ecx, c1E
		__asm	movaps	xmm0,xmm2		/*t3D copy */
		__asm	movaps	xmm1,xmm4		/*t3C copy */
		__asm	mulpd	xmm2,[ecx     ]	/*t3D*c1E */
		__asm	mulpd	xmm4,[ecx     ]	/*t3C*c1E */
		__asm	mulpd	xmm0,[ecx+0x10]	/*t3D*s1E */
		__asm	mulpd	xmm1,[ecx+0x10]	/*t3C*s1E */
		__asm	subpd	xmm4,xmm0	/* xmm4 <- im */
		__asm	addpd	xmm2,xmm1	/* xmm2 <- re	xmm0,1 free */

		__asm	movaps		xmm1,xmm4	/* cpy a[jp    ] */
		__asm	movaps		xmm0,xmm2	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm1,[ebx+0x1f0]
		__asm	unpcklpd	xmm4,[ebx+0x1f0]
		__asm	movaps	[ebx+0x1f0],xmm1	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm0,[ebx+0x1e0]
		__asm	unpcklpd	xmm2,[ebx+0x1e0]
		__asm	movaps	[ebx+0x1e0],xmm0	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x1f0],xmm4	/* a[jp+p12] */
		__asm	movaps	[eax+0x1e0],xmm2	/* a[jt+p12] */

	#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

	  #ifdef USE_AVX

		SSE2_RADIX32_WRAPPER_DIT(add0,add1,add2,add3,isrt2,r00,r08,r10,r20,r28,r30,c01,c02,c04,c06,c08,c0A,c0C,c0E,c10,c12,c14,c16,c18,c1A,c1C,c1E)

	  #else	// SSE2:

	   #if OS_BITS == 64
		SSE2_RADIX32_WRAPPER_DIT(add0,add1          ,isrt2,r00,r08,r10,r20,r28,r30    ,c01,c02    ,c04    ,c06    ,c08,c0A,c0C,c0E,c10,c12,c14,c16,c18,c1A,c1C,c1E)
	   #else	// 32-bit SSE2:
		SSE2_RADIX32_WRAPPER_DIT(add0,add1          ,isrt2,r00,r08,r10,r20,r28,r30,c00,c01,c02,c03,c04,c05,c06,c07,c08,c0A,c0C,c0E,c10,c12,c14,c16,c18,c1A,c1C,c1E)
	   #endif

	  #endif

	#endif

	}	// endif(j1 == 0)

#endif	/* USE_SSE2 */

/*...Update the data (j1 and j2) array indices. */
loop:
	  if(j1 <= 128) {
		// Use scalar code (with index offsets properly fiddled) for j1 == 0 case in SIMD mode
		j1 = j1+64;
		j2 = j2-64;
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

//	fprintf(stderr,"(j2_start == %d\n",j2_start);

	if(j2_start == n-64)
	{
		j1 = 0;
//	fprintf(stderr,"(j2_start == n-64) check hit: returning\n");	exit(0);
		return;
	}

/*...Reset half-complex-blocklength for next pass. If K >> 1 has a zero trailing bit, we multiply the blocklength by K >> 1 in preparation for the final block.	*/

	blocklen_sum = blocklen_sum + blocklen;
	blocklen = (radix_prim[i-1]-1)*blocklen_sum;

/*...Next j2_start is previous one plus the (real) length of the current block = 4*(half-complex-blocklength) */

	j2_start = j2_start+(blocklen<<2);
	j2=j2_start;				/* Reset j2 for start of the next block. */

}	/* End of Main loop */

}

