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

#ifdef USE_SSE2

	#include "sse2_macro.h"

	#ifdef COMPILER_TYPE_MSVC

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix32_wrapper_square_gcc32.h"

		#else

			#include "radix32_wrapper_square_gcc64.h"

		#endif

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
void radix32_dyadic_square(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int incr, int init_sse2, int thr_id)
{
	const int stride = (int)RE_IM_STRIDE << 6, stridh = (stride>>1);	// main-array loop stride = 32*RE_IM_STRIDE
	static int max_threads = 0;
	static int nsave = 0;
	static int rad0save = 0, ndivrad0 = 0;
/*	static int *index = 0x0;	OBSOLETE: full-N2/16-length Bit-reversal index array. */
	static int *index0 = 0x0, *index1 = 0x0, *index_ptmp0 = 0x0, *index_ptmp1 = 0x0;
		   int index0_idx=-1, index1_idx=-1;
	static int index0_mod=-1, index1_mod=-1;
	int nradices_prim_radix0;
	int i,j,j1,j2,l,iroot,k1,k2;
	const double c = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/

	double re0,im0,re1,im1,rt,it;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0, *add1;	/* Addresses into array sections */
  #ifdef USE_AVX
	double *add2, *add3;
  #endif
	vec_dbl *c_tmp,*s_tmp;
	vec_dbl *tmp;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
	// In || mode, only above base-pointer (shared by all threads) is static:
	vec_dbl *isrt2,*sqrt2, *one,*two, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*r00,*r08,*r10,*r18,*r20,*r28,*r30,*r38;
  #else
	static vec_dbl *isrt2,*sqrt2, *one,*two, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3,
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
	,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p1Ai,a1p1Bi,a1p1Ci,a1p1Di,a1p1Ei,a1p1Fi;

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
	if((rad0save != radix0) || (nsave != n))
	{
		ASSERT(HERE, thr_id == -1, "Init-mode call must be outside of any multithreading!");
		nsave = n;
		ASSERT(HERE, N2 == n/2, "N2 bad!");
		rad0save = radix0;
		ndivrad0 = n/radix0;
		if(index_ptmp0) {
			free((void *)index_ptmp0);	index_ptmp0=0x0;	index0=0x0;
			free((void *)index_ptmp1);	index_ptmp1=0x0;	index1=0x0;
		}
		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Allocate and initialize an index array containing N/32 indices...

		index_ptmp = ALLOC_INT(N2/32);
		if(!index_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array INDEX in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		index = ALIGN_INT(index_ptmp);
		if(!index){ sprintf(cbuf,"FATAL: unable to allocate array ITMP in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }

		for(i=0; i < N2/32; i++)
		{
			index[i]=i;
		}
		*/
		index0_mod =        radix0;
		index1_mod = (n>>6)/radix0;	/* complex length requires an additional divide by 2 */

		index_ptmp0 = ALLOC_INT(index_ptmp0, index0_mod);
		if(!index_ptmp0){ sprintf(cbuf,"FATAL: unable to allocate array INDEX_PTMP0 in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		index0 = ALIGN_INT(index_ptmp0);

		index_ptmp1 = ALLOC_INT(index_ptmp1, index1_mod);
		if(!index_ptmp1){ sprintf(cbuf,"FATAL: unable to allocate array INDEX_PTMP1 in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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
		ASSERT(HERE, nradices_prim_radix0 < nradices_prim,"radix32_dyadic_square.c: nradices_prim_radix0 < nradices_prim");

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

		ASSERT(HERE, thr_id == -1, "Init-mode call must be outside of any multithreading!");
		max_threads = init_sse2;
	#ifndef COMPILER_TYPE_GCC
		ASSERT(HERE, NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(HERE, max_threads >= NTHREADS, "Multithreading requires max_threads >= NTHREADS!");

	#ifdef USE_SSE2
		if(sc_arr != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)sc_arr);	sc_arr=0x0;
		}
		sc_arr = ALLOC_VEC_DBL(sc_arr, 0x94*max_threads + 100);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 64 vec_dbl slots of sc_arr for temporaries, next 8 for scratch, next 7 for the nontrivial complex 16th roots,
	next 62 for the doubled sincos twiddles, next 3 for [1.0,2.0,sqrt2] and at least 3 more to allow for 64-byte alignment of the array.

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
	ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
  #ifdef USE_SSE2
	r00  = __r0 + thr_id*148;isrt2	= r00 + 0x48;
/*	r02  = r00 + 0x02;*/cc0   = r00 + 0x49;
/*	r04  = r00 + 0x04;*/cc1   = r00 + 0x4b;
/*	r06  = r00 + 0x06;*/cc3   = r00 + 0x4d;
	r08  = r00 + 0x08;	c00   = r00 + 0x4f;
/*	r0A  = r00 + 0x0a;*/c10   = r00 + 0x51;
/*	r0C  = r00 + 0x0c;*/c08   = r00 + 0x53;
/*	r0E  = r00 + 0x0e;*/c18   = r00 + 0x55;
	r10  = r00 + 0x10;	c04   = r00 + 0x57;
/*	r12  = r00 + 0x12;*/c14   = r00 + 0x59;
/*	r14  = r00 + 0x14;*/c0C   = r00 + 0x5b;
/*	r16  = r00 + 0x16;*/c1C   = r00 + 0x5d;
	r18  = r00 + 0x18;	c02   = r00 + 0x5f;
/*	r1A  = r00 + 0x1a;*/c12   = r00 + 0x61;
/*	r1C  = r00 + 0x1c;*/c0A   = r00 + 0x63;
/*	r1E  = r00 + 0x1e;*/c1A   = r00 + 0x65;
	r20  = r00 + 0x20;	c06   = r00 + 0x67;
/*	r22  = r00 + 0x22;*/c16   = r00 + 0x69;
/*	r24  = r00 + 0x24;*/c0E   = r00 + 0x6b;
/*	r26  = r00 + 0x26;*/c1E   = r00 + 0x6d;
	r28  = r00 + 0x28;	c01   = r00 + 0x6f;
/*	r2A  = r00 + 0x2a;*/c11   = r00 + 0x71;
/*	r2C  = r00 + 0x2c;*/c09   = r00 + 0x73;
/*	r2E  = r00 + 0x2e;*/c19   = r00 + 0x75;
	r30  = r00 + 0x30;	c05   = r00 + 0x77;
/*	r32  = r00 + 0x32;*/c15   = r00 + 0x79;
/*	r34  = r00 + 0x34;*/c0D   = r00 + 0x7b;
/*	r36  = r00 + 0x36;*/c1D   = r00 + 0x7d;
	r38  = r00 + 0x38;	c03   = r00 + 0x7f;
/*	r3A  = r00 + 0x3a;*/c13   = r00 + 0x81;
/*	r3C  = r00 + 0x3c;*/c0B   = r00 + 0x83;
/*	r3E  = r00 + 0x3e;*/c1B   = r00 + 0x85;
						c07   = r00 + 0x87;
						c17   = r00 + 0x89;
						c0F   = r00 + 0x8b;
						c1F   = r00 + 0x8d;
						two   = r00 + 0x8f;
					//	forth = r00 + 0x90;
						sqrt2 = r00 + 0x91;
						one	  = r00 + 0x92;
  #endif
#endif
	/*...If a new runlength, should not get to this point: */
	ASSERT(HERE, n == nsave,"n != nsave");
	ASSERT(HERE, incr == 64,"incr == 64");
	ASSERT(HERE, ndivrad0 == n/radix0,"bad value for ndivrad0!");
	/*
	k = ii*(ndivrad0 >> 6);
	*/
	index0_idx = ii;
	index1_idx = 0;

	for(j = 0; j < ndivrad0; j += stride)
	{
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
		j2 = j1+RE_IM_STRIDE;
	/*
		iroot = index[k];
		ASSERT(HERE, iroot == index0[index0_idx] + index1[index1_idx],"radix32_dyadic_square.c: iroot == index0[index0_idx] + index1[index1_idx]");
		k = k + 1;	// increment sincos array index
	*/

		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod)
		{
			index1_idx -= index1_mod;
		ASSERT(HERE, index0_idx < index0_mod,"radix32_dyadic_square.c: index0_idx < index0_mod");
			++index0_idx;
		}

		l = iroot;

	#ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-31] are offset w.r.to the thread-local ptr pair as
					cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
		(cc0,ss0) + 0x[06,26,16,36,0e,2e,1e,3e,0a,2a,1a,3a,12,32,22,42,08,28,18,38,10,30,20,40,0b,2b,1b,3b,14,34,24,44].
		Here, due to the need to compute a new set of roots for each set of inputs, we use a streamlined sequence which
		computes only the [ 0, 1, 2, 3, 7,14,21,28]th roots with maximal accuracy (i.e. using 2-table-multiply),
		then generates the remaining ones from those. Thus the needed pointer offsets below are
			(cc0,ss0) + 0x[06,26,16,36,3e,22,30,14]:
		*/
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;
		c_tmp->d0=rt;	s_tmp->d0=it;	/* For all but the trivial 0th root, re and im parts in SSE2 mode use a separate set of roots, im-parts done below */
	#else
		c01 =rt;	s01 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		c02 =rt;	s02 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += 4*iroot;			/* 7*iroot	*/
		iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		c03 =rt;	s03 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		c07 =rt;	s07 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		c0E =rt;	s0E =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		c15 =rt;	s15 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;
		c_tmp->d0=rt;	s_tmp->d0=it;
	#else
		c1C =rt;	s1C =it;
	#endif

	/* In SSE2 mode, also need next set of sincoa to put into "Imaginary" slots: */
	#ifdef USE_SSE2
		// 2nd set:
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod)
		{
			index1_idx -= index1_mod;
		ASSERT(HERE, index0_idx < index0_mod,"radix32_dyadic_square.c: index0_idx < index0_mod");
			++index0_idx;
		}

		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;			/* 2*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;	/* c1,s1 */
		c_tmp->d1=rt;	s_tmp->d1=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;	/* c2,s2 */
		c_tmp->d1=rt;	s_tmp->d1=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += 4*iroot;			/* 7*iroot	*/
		iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;	/* c3,s3 */
		c_tmp->d1=rt;	s_tmp->d1=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;	/* c7,s7 */
		c_tmp->d1=rt;	s_tmp->d1=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;	/* c14,s14 */
		c_tmp->d1=rt;	s_tmp->d1=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;	/* c21,s21 */
		c_tmp->d1=rt;	s_tmp->d1=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;	/* c28,s28 */
		c_tmp->d1=rt;	s_tmp->d1=it;

	  /* In AVX mode, also need next 2 sets of sincos: */
	  #ifdef USE_AVX
		// 3rd set (of 4):
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod)
		{
			index1_idx -= index1_mod;
		ASSERT(HERE, index0_idx < index0_mod,"radix32_dyadic_square.c: index0_idx < index0_mod");
			++index0_idx;
		}

		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;			/* 2*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;	/* c1,s1 */
		c_tmp->d2=rt;	s_tmp->d2=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;	/* c2,s2 */
		c_tmp->d2=rt;	s_tmp->d2=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += 4*iroot;			/* 7*iroot	*/
		iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;	/* c3,s3 */
		c_tmp->d2=rt;	s_tmp->d2=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;	/* c7,s7 */
		c_tmp->d2=rt;	s_tmp->d2=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;	/* c14,s14 */
		c_tmp->d2=rt;	s_tmp->d2=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;	/* c21,s21 */
		c_tmp->d2=rt;	s_tmp->d2=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;	/* c28,s28 */
		c_tmp->d2=rt;	s_tmp->d2=it;

		// 4th set (of 4):
		iroot = index0[index0_idx] + index1[index1_idx];
		if(++index1_idx >= index1_mod)
		{
			index1_idx -= index1_mod;
		ASSERT(HERE, index0_idx < index0_mod,"radix32_dyadic_square.c: index0_idx < index0_mod");
			++index0_idx;
		}

		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;			/* 2*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;	/* c1,s1 */
		c_tmp->d3=rt;	s_tmp->d3=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 3*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;	/* c2,s2 */
		c_tmp->d3=rt;	s_tmp->d3=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += 4*iroot;			/* 7*iroot	*/
		iroot = l;				/* 7*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;	/* c3,s3 */
		c_tmp->d3=rt;	s_tmp->d3=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 14*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;	/* c7,s7 */
		c_tmp->d3=rt;	s_tmp->d3=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 21*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;	/* c14,s14 */
		c_tmp->d3=rt;	s_tmp->d3=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;				/* 28*iroot	*/
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;	/* c21,s21 */
		c_tmp->d3=rt;	s_tmp->d3=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;	/* c28,s28 */
		c_tmp->d3=rt;	s_tmp->d3=it;

	  #endif	// AVX?

		// Now generate the remaining roots from the smaller "seed set":
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

	#else
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
	#endif

/*       gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...
         We process the sincos data in bit-reversed order.	*/

	// See the radix-16 analog of this routine for details about the SSE2 and AVX data layouts

#if defined(USE_SSE2)	// SSE2 and AVX share same code flow, just with different versions of the ASM compute macros

		add0 = &a[j1];
		add1 = &a[j1+stridh];
	//	printf("stride = %d, add0,1 = %llX, %llX, diff = %llX\n",stride,(int64)add0,(int64)add1, (int64)add1-(int64)add0);	exit(0);
	//	ASSERT(HERE, (j1+stride) == (j+stride) + ( ((j+stride) >> DAT_BITS) << PAD_BITS ) , "add1 calculation violates padded index rules!");

	#ifdef COMPILER_TYPE_MSVC

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

	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */

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

	#else	/* GCC-style inline ASM: */

	  #ifdef USE_AVX	// The generic pre-dyadic-square macro needs 4 main-array addresses in AVX mode
	  					// because (add1-add0) and (add3-add2) have opposite signs for fermat and mersenne-mod:

		// process 4 main-array blocks of 8 vec_dbl = 8 x 4 = 32 doubles each in AVX mode
		add1 = add0 +  64;
		add2 = add0 + 128;
		add3 = add0 + 192;
		SSE2_RADIX32_WRAPPER_DIF(add0,add1,add2,add3,r00,r10,r20,r30,isrt2,cc0,c00,c01,c02,c03,c05,c07)

	  #else	// SSE2:

		SSE2_RADIX32_WRAPPER_DIF(add0,add1,          r00,r10,r20,r30,isrt2,cc0,c00,c01,c02,c03,c05,c07)

	  #endif

	#endif

	#ifdef COMPILER_TYPE_MSVC

		__asm	mov	eax, r00

		__asm	movaps	xmm0,[eax      ]				__asm	movaps	xmm3,[eax+0x100]
		__asm	movaps	xmm1,[eax+0x010]				__asm	movaps	xmm4,[eax+0x110]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax      ]				__asm	mulpd	xmm4,[eax+0x100]
		__asm	movaps	[eax      ],xmm0				__asm	movaps	[eax+0x100],xmm3
		__asm	movaps	[eax+0x010],xmm1				__asm	movaps	[eax+0x110],xmm4

		__asm	movaps	xmm0,[eax+0x020]				__asm	movaps	xmm3,[eax+0x120]
		__asm	movaps	xmm1,[eax+0x030]				__asm	movaps	xmm4,[eax+0x130]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x020]				__asm	mulpd	xmm4,[eax+0x120]
		__asm	movaps	[eax+0x020],xmm0				__asm	movaps	[eax+0x120],xmm3
		__asm	movaps	[eax+0x030],xmm1				__asm	movaps	[eax+0x130],xmm4

		__asm	movaps	xmm0,[eax+0x040]				__asm	movaps	xmm3,[eax+0x140]
		__asm	movaps	xmm1,[eax+0x050]				__asm	movaps	xmm4,[eax+0x150]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x040]				__asm	mulpd	xmm4,[eax+0x140]
		__asm	movaps	[eax+0x040],xmm0				__asm	movaps	[eax+0x140],xmm3
		__asm	movaps	[eax+0x050],xmm1				__asm	movaps	[eax+0x150],xmm4

		__asm	movaps	xmm0,[eax+0x060]				__asm	movaps	xmm3,[eax+0x160]
		__asm	movaps	xmm1,[eax+0x070]				__asm	movaps	xmm4,[eax+0x170]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x060]				__asm	mulpd	xmm4,[eax+0x160]
		__asm	movaps	[eax+0x060],xmm0				__asm	movaps	[eax+0x160],xmm3
		__asm	movaps	[eax+0x070],xmm1				__asm	movaps	[eax+0x170],xmm4

		__asm	movaps	xmm0,[eax+0x080]				__asm	movaps	xmm3,[eax+0x180]
		__asm	movaps	xmm1,[eax+0x090]				__asm	movaps	xmm4,[eax+0x190]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x080]				__asm	mulpd	xmm4,[eax+0x180]
		__asm	movaps	[eax+0x080],xmm0				__asm	movaps	[eax+0x180],xmm3
		__asm	movaps	[eax+0x090],xmm1				__asm	movaps	[eax+0x190],xmm4

		__asm	movaps	xmm0,[eax+0x0a0]				__asm	movaps	xmm3,[eax+0x1a0]
		__asm	movaps	xmm1,[eax+0x0b0]				__asm	movaps	xmm4,[eax+0x1b0]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x0a0]				__asm	mulpd	xmm4,[eax+0x1a0]
		__asm	movaps	[eax+0x0a0],xmm0				__asm	movaps	[eax+0x1a0],xmm3
		__asm	movaps	[eax+0x0b0],xmm1				__asm	movaps	[eax+0x1b0],xmm4

		__asm	movaps	xmm0,[eax+0x0c0]				__asm	movaps	xmm3,[eax+0x1c0]
		__asm	movaps	xmm1,[eax+0x0d0]				__asm	movaps	xmm4,[eax+0x1d0]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x0c0]				__asm	mulpd	xmm4,[eax+0x1c0]
		__asm	movaps	[eax+0x0c0],xmm0				__asm	movaps	[eax+0x1c0],xmm3
		__asm	movaps	[eax+0x0d0],xmm1				__asm	movaps	[eax+0x1d0],xmm4

		__asm	movaps	xmm0,[eax+0x0e0]				__asm	movaps	xmm3,[eax+0x1e0]
		__asm	movaps	xmm1,[eax+0x0f0]				__asm	movaps	xmm4,[eax+0x1f0]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x0e0]				__asm	mulpd	xmm4,[eax+0x1e0]
		__asm	movaps	[eax+0x0e0],xmm0				__asm	movaps	[eax+0x1e0],xmm3
		__asm	movaps	[eax+0x0f0],xmm1				__asm	movaps	[eax+0x1f0],xmm4

		__asm	movaps	xmm0,[eax+0x200]				__asm	movaps	xmm3,[eax+0x300]
		__asm	movaps	xmm1,[eax+0x210]				__asm	movaps	xmm4,[eax+0x310]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x200]				__asm	mulpd	xmm4,[eax+0x300]
		__asm	movaps	[eax+0x200],xmm0				__asm	movaps	[eax+0x300],xmm3
		__asm	movaps	[eax+0x210],xmm1				__asm	movaps	[eax+0x310],xmm4

		__asm	movaps	xmm0,[eax+0x220]				__asm	movaps	xmm3,[eax+0x320]
		__asm	movaps	xmm1,[eax+0x230]				__asm	movaps	xmm4,[eax+0x330]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x220]				__asm	mulpd	xmm4,[eax+0x320]
		__asm	movaps	[eax+0x220],xmm0				__asm	movaps	[eax+0x320],xmm3
		__asm	movaps	[eax+0x230],xmm1				__asm	movaps	[eax+0x330],xmm4

		__asm	movaps	xmm0,[eax+0x240]				__asm	movaps	xmm3,[eax+0x340]
		__asm	movaps	xmm1,[eax+0x250]				__asm	movaps	xmm4,[eax+0x350]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x240]				__asm	mulpd	xmm4,[eax+0x340]
		__asm	movaps	[eax+0x240],xmm0				__asm	movaps	[eax+0x340],xmm3
		__asm	movaps	[eax+0x250],xmm1				__asm	movaps	[eax+0x350],xmm4

		__asm	movaps	xmm0,[eax+0x260]				__asm	movaps	xmm3,[eax+0x360]
		__asm	movaps	xmm1,[eax+0x270]				__asm	movaps	xmm4,[eax+0x370]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x260]				__asm	mulpd	xmm4,[eax+0x360]
		__asm	movaps	[eax+0x260],xmm0				__asm	movaps	[eax+0x360],xmm3
		__asm	movaps	[eax+0x270],xmm1				__asm	movaps	[eax+0x370],xmm4

		__asm	movaps	xmm0,[eax+0x280]				__asm	movaps	xmm3,[eax+0x380]
		__asm	movaps	xmm1,[eax+0x290]				__asm	movaps	xmm4,[eax+0x390]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x280]				__asm	mulpd	xmm4,[eax+0x380]
		__asm	movaps	[eax+0x280],xmm0				__asm	movaps	[eax+0x380],xmm3
		__asm	movaps	[eax+0x290],xmm1				__asm	movaps	[eax+0x390],xmm4

		__asm	movaps	xmm0,[eax+0x2a0]				__asm	movaps	xmm3,[eax+0x3a0]
		__asm	movaps	xmm1,[eax+0x2b0]				__asm	movaps	xmm4,[eax+0x3b0]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x2a0]				__asm	mulpd	xmm4,[eax+0x3a0]
		__asm	movaps	[eax+0x2a0],xmm0				__asm	movaps	[eax+0x3a0],xmm3
		__asm	movaps	[eax+0x2b0],xmm1				__asm	movaps	[eax+0x3b0],xmm4

		__asm	movaps	xmm0,[eax+0x2c0]				__asm	movaps	xmm3,[eax+0x3c0]
		__asm	movaps	xmm1,[eax+0x2d0]				__asm	movaps	xmm4,[eax+0x3d0]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x2c0]				__asm	mulpd	xmm4,[eax+0x3c0]
		__asm	movaps	[eax+0x2c0],xmm0				__asm	movaps	[eax+0x3c0],xmm3
		__asm	movaps	[eax+0x2d0],xmm1				__asm	movaps	[eax+0x3d0],xmm4

		__asm	movaps	xmm0,[eax+0x2e0]				__asm	movaps	xmm3,[eax+0x3e0]
		__asm	movaps	xmm1,[eax+0x2f0]				__asm	movaps	xmm4,[eax+0x3f0]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x2e0]				__asm	mulpd	xmm4,[eax+0x3e0]
		__asm	movaps	[eax+0x2e0],xmm0				__asm	movaps	[eax+0x3e0],xmm3
		__asm	movaps	[eax+0x2f0],xmm1				__asm	movaps	[eax+0x3f0],xmm4

	#else	/* GCC-style inline ASM: */

	  #if OS_BITS == 32

		__asm__ volatile (\
			"movl	%[__r00],%%eax		\n\t"\
			"/* z0^2: */				\n\t	/* z0^2: */				\n\t"\
			"movaps		 (%%eax),%%xmm0	\n\t	movaps	0x200(%%eax),%%xmm3	\n\t"\
			"movaps	0x010(%%eax),%%xmm1	\n\t	movaps	0x210(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd		 (%%eax),%%xmm1	\n\t	mulpd	0x200(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,     (%%eax)	\n\t	movaps	%%xmm3,0x200(%%eax)	\n\t"\
			"movaps	%%xmm1,0x010(%%eax)	\n\t	movaps	%%xmm4,0x210(%%eax)	\n\t"\
			"/* z1^2: */				\n\t	/* z1^2: */				\n\t"\
			"movaps	0x020(%%eax),%%xmm0	\n\t	movaps	0x220(%%eax),%%xmm3	\n\t"\
			"movaps	0x030(%%eax),%%xmm1	\n\t	movaps	0x230(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x020(%%eax),%%xmm1	\n\t	mulpd	0x220(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x020(%%eax)	\n\t	movaps	%%xmm3,0x220(%%eax)	\n\t"\
			"movaps	%%xmm1,0x030(%%eax)	\n\t	movaps	%%xmm4,0x230(%%eax)	\n\t"\
			"/* z2^2: */				\n\t	/* z2^2: */				\n\t"\
			"movaps	0x040(%%eax),%%xmm0	\n\t	movaps	0x240(%%eax),%%xmm3	\n\t"\
			"movaps	0x050(%%eax),%%xmm1	\n\t	movaps	0x250(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x040(%%eax),%%xmm1	\n\t	mulpd	0x240(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x040(%%eax)	\n\t	movaps	%%xmm3,0x240(%%eax)	\n\t"\
			"movaps	%%xmm1,0x050(%%eax)	\n\t	movaps	%%xmm4,0x250(%%eax)	\n\t"\
			"/* z3^2: */				\n\t	/* z3^2: */				\n\t"\
			"movaps	0x060(%%eax),%%xmm0	\n\t	movaps	0x260(%%eax),%%xmm3	\n\t"\
			"movaps	0x070(%%eax),%%xmm1	\n\t	movaps	0x270(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x060(%%eax),%%xmm1	\n\t	mulpd	0x260(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x060(%%eax)	\n\t	movaps	%%xmm3,0x260(%%eax)	\n\t"\
			"movaps	%%xmm1,0x070(%%eax)	\n\t	movaps	%%xmm4,0x270(%%eax)	\n\t"\
			"/* z4^2: */				\n\t	/* z4^2: */				\n\t"\
			"movaps	0x080(%%eax),%%xmm0	\n\t	movaps	0x280(%%eax),%%xmm3	\n\t"\
			"movaps	0x090(%%eax),%%xmm1	\n\t	movaps	0x290(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x080(%%eax),%%xmm1	\n\t	mulpd	0x280(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x080(%%eax)	\n\t	movaps	%%xmm3,0x280(%%eax)	\n\t"\
			"movaps	%%xmm1,0x090(%%eax)	\n\t	movaps	%%xmm4,0x290(%%eax)	\n\t"\
			"/* z5^2: */				\n\t	/* z5^2: */				\n\t"\
			"movaps	0x0a0(%%eax),%%xmm0	\n\t	movaps	0x2a0(%%eax),%%xmm3	\n\t"\
			"movaps	0x0b0(%%eax),%%xmm1	\n\t	movaps	0x2b0(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0a0(%%eax),%%xmm1	\n\t	mulpd	0x2a0(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0a0(%%eax)	\n\t	movaps	%%xmm3,0x2a0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x0b0(%%eax)	\n\t	movaps	%%xmm4,0x2b0(%%eax)	\n\t"\
			"/* z6^2: */				\n\t	/* z6^2: */				\n\t"\
			"movaps	0x0c0(%%eax),%%xmm0	\n\t	movaps	0x2c0(%%eax),%%xmm3	\n\t"\
			"movaps	0x0d0(%%eax),%%xmm1	\n\t	movaps	0x2d0(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0c0(%%eax),%%xmm1	\n\t	mulpd	0x2c0(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0c0(%%eax)	\n\t	movaps	%%xmm3,0x2c0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x0d0(%%eax)	\n\t	movaps	%%xmm4,0x2d0(%%eax)	\n\t"\
			"/* z7^2: */				\n\t	/* z7^2: */				\n\t"\
			"movaps	0x0e0(%%eax),%%xmm0	\n\t	movaps	0x2e0(%%eax),%%xmm3	\n\t"\
			"movaps	0x0f0(%%eax),%%xmm1	\n\t	movaps	0x2f0(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0e0(%%eax),%%xmm1	\n\t	mulpd	0x2e0(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0e0(%%eax)	\n\t	movaps	%%xmm3,0x2e0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x0f0(%%eax)	\n\t	movaps	%%xmm4,0x2f0(%%eax)	\n\t"\
			"/* z8^2: */				\n\t	/* z8^2: */				\n\t"\
			"movaps	0x100(%%eax),%%xmm0	\n\t	movaps	0x300(%%eax),%%xmm3	\n\t"\
			"movaps	0x110(%%eax),%%xmm1	\n\t	movaps	0x310(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x100(%%eax),%%xmm1	\n\t	mulpd	0x300(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x100(%%eax)	\n\t	movaps	%%xmm3,0x300(%%eax)	\n\t"\
			"movaps	%%xmm1,0x110(%%eax)	\n\t	movaps	%%xmm4,0x310(%%eax)	\n\t"\
			"/* z9^2: */				\n\t	/* z9^2: */				\n\t"\
			"movaps	0x120(%%eax),%%xmm0	\n\t	movaps	0x320(%%eax),%%xmm3	\n\t"\
			"movaps	0x130(%%eax),%%xmm1	\n\t	movaps	0x330(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x120(%%eax),%%xmm1	\n\t	mulpd	0x320(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x120(%%eax)	\n\t	movaps	%%xmm3,0x320(%%eax)	\n\t"\
			"movaps	%%xmm1,0x130(%%eax)	\n\t	movaps	%%xmm4,0x330(%%eax)	\n\t"\
			"/* zA^2: */				\n\t	/* zA^2: */				\n\t"\
			"movaps	0x140(%%eax),%%xmm0	\n\t	movaps	0x340(%%eax),%%xmm3	\n\t"\
			"movaps	0x150(%%eax),%%xmm1	\n\t	movaps	0x350(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x140(%%eax),%%xmm1	\n\t	mulpd	0x340(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x140(%%eax)	\n\t	movaps	%%xmm3,0x340(%%eax)	\n\t"\
			"movaps	%%xmm1,0x150(%%eax)	\n\t	movaps	%%xmm4,0x350(%%eax)	\n\t"\
			"/* zB^2: */				\n\t	/* zB^2: */				\n\t"\
			"movaps	0x160(%%eax),%%xmm0	\n\t	movaps	0x360(%%eax),%%xmm3	\n\t"\
			"movaps	0x170(%%eax),%%xmm1	\n\t	movaps	0x370(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x160(%%eax),%%xmm1	\n\t	mulpd	0x360(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x160(%%eax)	\n\t	movaps	%%xmm3,0x360(%%eax)	\n\t"\
			"movaps	%%xmm1,0x170(%%eax)	\n\t	movaps	%%xmm4,0x370(%%eax)	\n\t"\
			"/* zC^2: */				\n\t	/* zC^2: */				\n\t"\
			"movaps	0x180(%%eax),%%xmm0	\n\t	movaps	0x380(%%eax),%%xmm3	\n\t"\
			"movaps	0x190(%%eax),%%xmm1	\n\t	movaps	0x390(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x180(%%eax),%%xmm1	\n\t	mulpd	0x380(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x180(%%eax)	\n\t	movaps	%%xmm3,0x380(%%eax)	\n\t"\
			"movaps	%%xmm1,0x190(%%eax)	\n\t	movaps	%%xmm4,0x390(%%eax)	\n\t"\
			"/* zD^2: */				\n\t	/* zD^2: */				\n\t"\
			"movaps	0x1a0(%%eax),%%xmm0	\n\t	movaps	0x3a0(%%eax),%%xmm3	\n\t"\
			"movaps	0x1b0(%%eax),%%xmm1	\n\t	movaps	0x3b0(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x1a0(%%eax),%%xmm1	\n\t	mulpd	0x3a0(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x1a0(%%eax)	\n\t	movaps	%%xmm3,0x3a0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x1b0(%%eax)	\n\t	movaps	%%xmm4,0x3b0(%%eax)	\n\t"\
			"/* zE^2: */				\n\t	/* zE^2: */				\n\t"\
			"movaps	0x1c0(%%eax),%%xmm0	\n\t	movaps	0x3c0(%%eax),%%xmm3	\n\t"\
			"movaps	0x1d0(%%eax),%%xmm1	\n\t	movaps	0x3d0(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x1c0(%%eax),%%xmm1	\n\t	mulpd	0x3c0(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x1c0(%%eax)	\n\t	movaps	%%xmm3,0x3c0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x1d0(%%eax)	\n\t	movaps	%%xmm4,0x3d0(%%eax)	\n\t"\
			"/* zF^2: */				\n\t	/* zF^2: */				\n\t"\
			"movaps	0x1e0(%%eax),%%xmm0	\n\t	movaps	0x3e0(%%eax),%%xmm3	\n\t"\
			"movaps	0x1f0(%%eax),%%xmm1	\n\t	movaps	0x3f0(%%eax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x1e0(%%eax),%%xmm1	\n\t	mulpd	0x3e0(%%eax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x1e0(%%eax)	\n\t	movaps	%%xmm3,0x3e0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x1f0(%%eax)	\n\t	movaps	%%xmm4,0x3f0(%%eax)	\n\t"\
			:					// outputs: none
			: [__r00] "m" (r00)	// All inputs from memory addresses here
			: "cc","memory","eax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	// Clobbered registers
		);

	  #else

		#ifdef USE_AVX

		__asm__ volatile (\
			"movq	%[__r00],%%rax			\n\t"\
			"/* z0^2: */					\n\t	/* z0^2: */				\n\t"\
			"vmovaps		 (%%rax),%%ymm0	\n\t	vmovaps	0x400(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x420(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd		 (%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x400(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm3,0x400(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm4,0x420(%%rax)	\n\t"\
			"/* z1^2: */					\n\t	/* z1^2: */				\n\t"\
			"vmovaps	0x040(%%rax),%%ymm0	\n\t	vmovaps	0x440(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x060(%%rax),%%ymm1	\n\t	vmovaps	0x460(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x040(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x440(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x040(%%rax)	\n\t	vmovaps	%%ymm3,0x440(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x060(%%rax)	\n\t	vmovaps	%%ymm4,0x460(%%rax)	\n\t"\
			"/* z2^2: */					\n\t	/* z2^2: */				\n\t"\
			"vmovaps	0x080(%%rax),%%ymm0	\n\t	vmovaps	0x480(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x0a0(%%rax),%%ymm1	\n\t	vmovaps	0x4a0(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x080(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x480(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x080(%%rax)	\n\t	vmovaps	%%ymm3,0x480(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x0a0(%%rax)	\n\t	vmovaps	%%ymm4,0x4a0(%%rax)	\n\t"\
			"/* z3^2: */					\n\t	/* z3^2: */				\n\t"\
			"vmovaps	0x0c0(%%rax),%%ymm0	\n\t	vmovaps	0x4c0(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x0e0(%%rax),%%ymm1	\n\t	vmovaps	0x4e0(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x0c0(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x4c0(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x0c0(%%rax)	\n\t	vmovaps	%%ymm3,0x4c0(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x0e0(%%rax)	\n\t	vmovaps	%%ymm4,0x4e0(%%rax)	\n\t"\
			"/* z4^2: */					\n\t	/* z4^2: */				\n\t"\
			"vmovaps	0x100(%%rax),%%ymm0	\n\t	vmovaps	0x500(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x120(%%rax),%%ymm1	\n\t	vmovaps	0x520(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x100(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x500(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x100(%%rax)	\n\t	vmovaps	%%ymm3,0x500(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x120(%%rax)	\n\t	vmovaps	%%ymm4,0x520(%%rax)	\n\t"\
			"/* z5^2: */					\n\t	/* z5^2: */				\n\t"\
			"vmovaps	0x140(%%rax),%%ymm0	\n\t	vmovaps	0x540(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x160(%%rax),%%ymm1	\n\t	vmovaps	0x560(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x140(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x540(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x140(%%rax)	\n\t	vmovaps	%%ymm3,0x540(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x160(%%rax)	\n\t	vmovaps	%%ymm4,0x560(%%rax)	\n\t"\
			"/* z6^2: */					\n\t	/* z6^2: */				\n\t"\
			"vmovaps	0x180(%%rax),%%ymm0	\n\t	vmovaps	0x580(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x1a0(%%rax),%%ymm1	\n\t	vmovaps	0x5a0(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x180(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x580(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x180(%%rax)	\n\t	vmovaps	%%ymm3,0x580(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x1a0(%%rax)	\n\t	vmovaps	%%ymm4,0x5a0(%%rax)	\n\t"\
			"/* z7^2: */					\n\t	/* z7^2: */				\n\t"\
			"vmovaps	0x1c0(%%rax),%%ymm0	\n\t	vmovaps	0x5c0(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x1e0(%%rax),%%ymm1	\n\t	vmovaps	0x5e0(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x1c0(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x5c0(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x1c0(%%rax)	\n\t	vmovaps	%%ymm3,0x5c0(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x1e0(%%rax)	\n\t	vmovaps	%%ymm4,0x5e0(%%rax)	\n\t"\
			"/* z8^2: */					\n\t	/* z8^2: */				\n\t"\
			"vmovaps	0x200(%%rax),%%ymm0	\n\t	vmovaps	0x600(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x220(%%rax),%%ymm1	\n\t	vmovaps	0x620(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x200(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x600(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x200(%%rax)	\n\t	vmovaps	%%ymm3,0x600(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x220(%%rax)	\n\t	vmovaps	%%ymm4,0x620(%%rax)	\n\t"\
			"/* z9^2: */					\n\t	/* z9^2: */				\n\t"\
			"vmovaps	0x240(%%rax),%%ymm0	\n\t	vmovaps	0x640(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x260(%%rax),%%ymm1	\n\t	vmovaps	0x660(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x240(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x640(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x240(%%rax)	\n\t	vmovaps	%%ymm3,0x640(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x260(%%rax)	\n\t	vmovaps	%%ymm4,0x660(%%rax)	\n\t"\
			"/* zA^2: */					\n\t	/* zA^2: */				\n\t"\
			"vmovaps	0x280(%%rax),%%ymm0	\n\t	vmovaps	0x680(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x2a0(%%rax),%%ymm1	\n\t	vmovaps	0x6a0(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x280(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x680(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x280(%%rax)	\n\t	vmovaps	%%ymm3,0x680(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x2a0(%%rax)	\n\t	vmovaps	%%ymm4,0x6a0(%%rax)	\n\t"\
			"/* zB^2: */					\n\t	/* zB^2: */				\n\t"\
			"vmovaps	0x2c0(%%rax),%%ymm0	\n\t	vmovaps	0x6c0(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x2e0(%%rax),%%ymm1	\n\t	vmovaps	0x6e0(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x2c0(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x6c0(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x2c0(%%rax)	\n\t	vmovaps	%%ymm3,0x6c0(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x2e0(%%rax)	\n\t	vmovaps	%%ymm4,0x6e0(%%rax)	\n\t"\
			"/* zC^2: */					\n\t	/* zC^2: */				\n\t"\
			"vmovaps	0x300(%%rax),%%ymm0	\n\t	vmovaps	0x700(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x320(%%rax),%%ymm1	\n\t	vmovaps	0x720(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x300(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x700(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x300(%%rax)	\n\t	vmovaps	%%ymm3,0x700(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x320(%%rax)	\n\t	vmovaps	%%ymm4,0x720(%%rax)	\n\t"\
			"/* zD^2: */					\n\t	/* zD^2: */				\n\t"\
			"vmovaps	0x340(%%rax),%%ymm0	\n\t	vmovaps	0x740(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x360(%%rax),%%ymm1	\n\t	vmovaps	0x760(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x340(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x740(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x340(%%rax)	\n\t	vmovaps	%%ymm3,0x740(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x360(%%rax)	\n\t	vmovaps	%%ymm4,0x760(%%rax)	\n\t"\
			"/* zE^2: */					\n\t	/* zE^2: */				\n\t"\
			"vmovaps	0x380(%%rax),%%ymm0	\n\t	vmovaps	0x780(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x3a0(%%rax),%%ymm1	\n\t	vmovaps	0x7a0(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x380(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x780(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x380(%%rax)	\n\t	vmovaps	%%ymm3,0x780(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x3a0(%%rax)	\n\t	vmovaps	%%ymm4,0x7a0(%%rax)	\n\t"\
			"/* zF^2: */					\n\t	/* zF^2: */				\n\t"\
			"vmovaps	0x3c0(%%rax),%%ymm0	\n\t	vmovaps	0x7c0(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x3e0(%%rax),%%ymm1	\n\t	vmovaps	0x7e0(%%rax),%%ymm4	\n\t"\
			"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
			"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
			"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
			"vmulpd	0x3c0(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x7c0(%%rax),%%ymm4,%%ymm4		\n\t"\
			"vmovaps	%%ymm0,0x3c0(%%rax)	\n\t	vmovaps	%%ymm3,0x7c0(%%rax)	\n\t"\
			"vmovaps	%%ymm1,0x3e0(%%rax)	\n\t	vmovaps	%%ymm4,0x7e0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r00] "m" (r00)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	// Clobbered registers
		);

		#elif defined(USE_SSE2)

		__asm__ volatile (\
			"movq	%[__r00],%%rax		\n\t"\
			"/* z0^2: */				\n\t	/* z0^2: */				\n\t"\
			"movaps		 (%%rax),%%xmm0	\n\t	movaps	0x200(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t	movaps	0x210(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd		 (%%rax),%%xmm1	\n\t	mulpd	0x200(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm3,0x200(%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm4,0x210(%%rax)	\n\t"\
			"/* z1^2: */				\n\t	/* z1^2: */				\n\t"\
			"movaps	0x020(%%rax),%%xmm0	\n\t	movaps	0x220(%%rax),%%xmm3	\n\t"\
			"movaps	0x030(%%rax),%%xmm1	\n\t	movaps	0x230(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x020(%%rax),%%xmm1	\n\t	mulpd	0x220(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x020(%%rax)	\n\t	movaps	%%xmm3,0x220(%%rax)	\n\t"\
			"movaps	%%xmm1,0x030(%%rax)	\n\t	movaps	%%xmm4,0x230(%%rax)	\n\t"\
			"/* z2^2: */				\n\t	/* z2^2: */				\n\t"\
			"movaps	0x040(%%rax),%%xmm0	\n\t	movaps	0x240(%%rax),%%xmm3	\n\t"\
			"movaps	0x050(%%rax),%%xmm1	\n\t	movaps	0x250(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x040(%%rax),%%xmm1	\n\t	mulpd	0x240(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x040(%%rax)	\n\t	movaps	%%xmm3,0x240(%%rax)	\n\t"\
			"movaps	%%xmm1,0x050(%%rax)	\n\t	movaps	%%xmm4,0x250(%%rax)	\n\t"\
			"/* z3^2: */				\n\t	/* z3^2: */				\n\t"\
			"movaps	0x060(%%rax),%%xmm0	\n\t	movaps	0x260(%%rax),%%xmm3	\n\t"\
			"movaps	0x070(%%rax),%%xmm1	\n\t	movaps	0x270(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x060(%%rax),%%xmm1	\n\t	mulpd	0x260(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x060(%%rax)	\n\t	movaps	%%xmm3,0x260(%%rax)	\n\t"\
			"movaps	%%xmm1,0x070(%%rax)	\n\t	movaps	%%xmm4,0x270(%%rax)	\n\t"\
			"/* z4^2: */				\n\t	/* z4^2: */				\n\t"\
			"movaps	0x080(%%rax),%%xmm0	\n\t	movaps	0x280(%%rax),%%xmm3	\n\t"\
			"movaps	0x090(%%rax),%%xmm1	\n\t	movaps	0x290(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x080(%%rax),%%xmm1	\n\t	mulpd	0x280(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x080(%%rax)	\n\t	movaps	%%xmm3,0x280(%%rax)	\n\t"\
			"movaps	%%xmm1,0x090(%%rax)	\n\t	movaps	%%xmm4,0x290(%%rax)	\n\t"\
			"/* z5^2: */				\n\t	/* z5^2: */				\n\t"\
			"movaps	0x0a0(%%rax),%%xmm0	\n\t	movaps	0x2a0(%%rax),%%xmm3	\n\t"\
			"movaps	0x0b0(%%rax),%%xmm1	\n\t	movaps	0x2b0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0a0(%%rax),%%xmm1	\n\t	mulpd	0x2a0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0a0(%%rax)	\n\t	movaps	%%xmm3,0x2a0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0b0(%%rax)	\n\t	movaps	%%xmm4,0x2b0(%%rax)	\n\t"\
			"/* z6^2: */				\n\t	/* z6^2: */				\n\t"\
			"movaps	0x0c0(%%rax),%%xmm0	\n\t	movaps	0x2c0(%%rax),%%xmm3	\n\t"\
			"movaps	0x0d0(%%rax),%%xmm1	\n\t	movaps	0x2d0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0c0(%%rax),%%xmm1	\n\t	mulpd	0x2c0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0c0(%%rax)	\n\t	movaps	%%xmm3,0x2c0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0d0(%%rax)	\n\t	movaps	%%xmm4,0x2d0(%%rax)	\n\t"\
			"/* z7^2: */				\n\t	/* z7^2: */				\n\t"\
			"movaps	0x0e0(%%rax),%%xmm0	\n\t	movaps	0x2e0(%%rax),%%xmm3	\n\t"\
			"movaps	0x0f0(%%rax),%%xmm1	\n\t	movaps	0x2f0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x0e0(%%rax),%%xmm1	\n\t	mulpd	0x2e0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x0e0(%%rax)	\n\t	movaps	%%xmm3,0x2e0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0f0(%%rax)	\n\t	movaps	%%xmm4,0x2f0(%%rax)	\n\t"\
			"/* z8^2: */				\n\t	/* z8^2: */				\n\t"\
			"movaps	0x100(%%rax),%%xmm0	\n\t	movaps	0x300(%%rax),%%xmm3	\n\t"\
			"movaps	0x110(%%rax),%%xmm1	\n\t	movaps	0x310(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x100(%%rax),%%xmm1	\n\t	mulpd	0x300(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x100(%%rax)	\n\t	movaps	%%xmm3,0x300(%%rax)	\n\t"\
			"movaps	%%xmm1,0x110(%%rax)	\n\t	movaps	%%xmm4,0x310(%%rax)	\n\t"\
			"/* z9^2: */				\n\t	/* z9^2: */				\n\t"\
			"movaps	0x120(%%rax),%%xmm0	\n\t	movaps	0x320(%%rax),%%xmm3	\n\t"\
			"movaps	0x130(%%rax),%%xmm1	\n\t	movaps	0x330(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x120(%%rax),%%xmm1	\n\t	mulpd	0x320(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x120(%%rax)	\n\t	movaps	%%xmm3,0x320(%%rax)	\n\t"\
			"movaps	%%xmm1,0x130(%%rax)	\n\t	movaps	%%xmm4,0x330(%%rax)	\n\t"\
			"/* zA^2: */				\n\t	/* zA^2: */				\n\t"\
			"movaps	0x140(%%rax),%%xmm0	\n\t	movaps	0x340(%%rax),%%xmm3	\n\t"\
			"movaps	0x150(%%rax),%%xmm1	\n\t	movaps	0x350(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x140(%%rax),%%xmm1	\n\t	mulpd	0x340(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x140(%%rax)	\n\t	movaps	%%xmm3,0x340(%%rax)	\n\t"\
			"movaps	%%xmm1,0x150(%%rax)	\n\t	movaps	%%xmm4,0x350(%%rax)	\n\t"\
			"/* zB^2: */				\n\t	/* zB^2: */				\n\t"\
			"movaps	0x160(%%rax),%%xmm0	\n\t	movaps	0x360(%%rax),%%xmm3	\n\t"\
			"movaps	0x170(%%rax),%%xmm1	\n\t	movaps	0x370(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x160(%%rax),%%xmm1	\n\t	mulpd	0x360(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x160(%%rax)	\n\t	movaps	%%xmm3,0x360(%%rax)	\n\t"\
			"movaps	%%xmm1,0x170(%%rax)	\n\t	movaps	%%xmm4,0x370(%%rax)	\n\t"\
			"/* zC^2: */				\n\t	/* zC^2: */				\n\t"\
			"movaps	0x180(%%rax),%%xmm0	\n\t	movaps	0x380(%%rax),%%xmm3	\n\t"\
			"movaps	0x190(%%rax),%%xmm1	\n\t	movaps	0x390(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x180(%%rax),%%xmm1	\n\t	mulpd	0x380(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x180(%%rax)	\n\t	movaps	%%xmm3,0x380(%%rax)	\n\t"\
			"movaps	%%xmm1,0x190(%%rax)	\n\t	movaps	%%xmm4,0x390(%%rax)	\n\t"\
			"/* zD^2: */				\n\t	/* zD^2: */				\n\t"\
			"movaps	0x1a0(%%rax),%%xmm0	\n\t	movaps	0x3a0(%%rax),%%xmm3	\n\t"\
			"movaps	0x1b0(%%rax),%%xmm1	\n\t	movaps	0x3b0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x1a0(%%rax),%%xmm1	\n\t	mulpd	0x3a0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x1a0(%%rax)	\n\t	movaps	%%xmm3,0x3a0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x1b0(%%rax)	\n\t	movaps	%%xmm4,0x3b0(%%rax)	\n\t"\
			"/* zE^2: */				\n\t	/* zE^2: */				\n\t"\
			"movaps	0x1c0(%%rax),%%xmm0	\n\t	movaps	0x3c0(%%rax),%%xmm3	\n\t"\
			"movaps	0x1d0(%%rax),%%xmm1	\n\t	movaps	0x3d0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x1c0(%%rax),%%xmm1	\n\t	mulpd	0x3c0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x1c0(%%rax)	\n\t	movaps	%%xmm3,0x3c0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x1d0(%%rax)	\n\t	movaps	%%xmm4,0x3d0(%%rax)	\n\t"\
			"/* zF^2: */				\n\t	/* zF^2: */				\n\t"\
			"movaps	0x1e0(%%rax),%%xmm0	\n\t	movaps	0x3e0(%%rax),%%xmm3	\n\t"\
			"movaps	0x1f0(%%rax),%%xmm1	\n\t	movaps	0x3f0(%%rax),%%xmm4	\n\t"\
			"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
			"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
			"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
			"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
			"mulpd	0x1e0(%%rax),%%xmm1	\n\t	mulpd	0x3e0(%%rax),%%xmm4	\n\t"\
			"movaps	%%xmm0,0x1e0(%%rax)	\n\t	movaps	%%xmm3,0x3e0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x1f0(%%rax)	\n\t	movaps	%%xmm4,0x3f0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r00] "m" (r00)	// All inputs from memory addresses here
			: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	// Clobbered registers
		);

	   #endif	// AVX?

	  #endif	// 32 or 64-it GCC?

	#endif

	/*...And do an inverse DIT radix-16 pass on the squared-data blocks. */

	#if defined(COMPILER_TYPE_MSVC)

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

	#else	/* GCC-style inline ASM: */

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

#else	/* USE_SSE2 */

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

/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/

/* 1st set of inputs: */
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

	}	/* endfor() */

}

