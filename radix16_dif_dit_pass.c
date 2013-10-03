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

/* Use for testing higher-accuracy version of the twiddles computation */
#define HIACC 1
#define EPS	1e-10

#ifndef PFETCH_DIST
	#define PFETCH_DIST	8
#endif

#ifdef USE_SSE2

	#if(HIACC != 1)
		#error SIMD Mode requires HIACC flag to be set!
	#endif

	#ifdef COMPILER_TYPE_MSVC
		#include "sse2_macro.h"
	#endif

/*	Recipe for MSVC --> GCC inline ASM conversion:

Before you begin tranmslation:
	- Max. number of input variables GCC allows = 30 ... if you're using more than that,
	trying reducing the count e.g. by using var2 = var1 + memoffset in the ASM.
	DO THIS USING THE MSVC CODE, i.e. only *after* you've successfully reduced
	the inline ASM macro arg count should you proceed with syntax translation.
	That allows you to work through small chunks of inline ASM at a time, doing
	quick-build-and-debug to check the changes, i.e. greatly eases debug.

	0. Remove all but most-crucial comments to ease conversion, as follows:
		[blockmode] space all "keeper" comments to extreme left
		multistatement __asm lines --> one __asm pre line, realign __asm to left-justify, delete __asm\t
		For non-keeper comments: /* --> @@
		[regexp] @@ --> \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t@@
		(delete all @@... stuff)\
		\t\n --> \n (repeat until no more trailing tabs)
		[/regexp]
		Repeat /* --> @@, [regexp] @@ --> \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t@@ steps for keeper comments, to move out of harm's way.
	1. [...] --> (...)
	2. ALU ops [e.g. mov, add, shl] --> spaces to tabs, then append "l" [if e*x] or "q" [if r*x] to instruction name
	3. Numeric literals in above kinds of instructions: Prepend "$" [",0x" --> ",$0x"]
	4. Address offsets of form (...+0x100) --> 0x100(...), (...-0x100) --> -0x100(...)
	5. External variable names get wrapped in %[]
	6. Line up commas in vertically stacked columns, then reverse operand order columnwise [for both 2 and 3-operand instructions].
	7. Prepend "%%" to all register names
	8. Only e*x/r*x registers appear in clobber list, not special regs like mmx and xmm.

Additional Notes:
	- Need to strip off any leading white space from named vars inside [], e.g. for "movl %[  c4],%%ecx \n\t" get "undefined named operand '  c4'" error;
	- Offsets with explicit + sign, e.g. "+0x10(%%eax)", not allowed
*/

	#if(defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))

		#if OS_BITS == 32

			#include "radix16_dif_dit_pass_gcc32.h"

		#else

			#include "radix16_dif_dit_pass_gcc64.h"

		#endif

	#endif

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
void radix16_dif_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
{
	const int RADIX = 16;
	const int pfetch_dist = PFETCH_DIST;
	int pfetch_addr;
	static int max_threads = 0;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride (doubles) = 2*RE_IM_STRIDE
	// lg(stride):
  #ifdef USE_AVX
	const int l2_stride = 3;	// 8 doubles at a time
  #elif defined(USE_SSE2)
	const int l2_stride = 2;	// 4 doubles at a time
  #else
	const int l2_stride = 1;	// 2 doubles at a time in scalar [non-SIMD] mode
  #endif

	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]*/
#ifdef USE_AVX2	// FMA-based DFT in AVX2 mode needs the tangent
	static double tan = 0.41421356237309504879;
#endif
	int i,j,j1,j2,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	int p1,p2,p3,p4,p8,p12;
	double rt,it,dtmp;
	double re0,im0,re1,im1;
	uint64 tmp64;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0, *add1, *add2;	/* Addresses into array sections */
#ifdef COMPILER_TYPE_MSVC
	double *add3;
#endif
	vec_dbl *c_tmp,*s_tmp;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
	vec_dbl *cc0, *ss0, *isrt2, *two, *r1;	// In || mode, only above base-pointer (shared by all threads) is static:
  #elif defined(COMPILER_TYPE_GCC)
	static vec_dbl *cc0, *ss0, *isrt2, *two, *r1;
  #else
	static vec_dbl *cc0, *ss0, *isrt2, *two;
	static vec_dbl *c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15,*s0,*s1,*s2,*s3,*s4,*s5,*s6,*s7,*s8,*s9,*s10,*s11,*s12,*s13,*s14,*s15;
	static vec_dbl *r1,*r2,*r3,*r4,*r5,*r6,*r7,*r8,*r9,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r30,*r31,*r32;
  #endif

#else

	int jp,jt;
	double *addr, *addp;
	int prefetch_offset;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
	,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;

#endif

#ifdef USE_SSE2

	/* In order to eschew complex thread-block-and-sync logic related to the local-store-init step in multithread mode,
	switch to a special-init-mode-call paradigm, in which the function is inited once (for as many threads as needed)
	prior to being executed:
	*/
	if(init_sse2)	// Just check (init_sse2 != 0) here, to allow the *value* of init_sse2 to store #threads
	{
		fprintf(stderr, "radix16_dif_dit_pass pfetch_dist = %d\n", pfetch_dist);
		max_threads = init_sse2;
	#ifndef COMPILER_TYPE_GCC
		ASSERT(HERE, NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(HERE, max_threads >= NTHREADS, "Multithreading requires max_threads >= NTHREADS!");
		ASSERT(HERE, thr_id == -1, "Init-mode call must be outside of any multithreading!");
		sc_arr = ALLOC_VEC_DBL(sc_arr, 72*max_threads);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	last 30 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
		#ifdef MULTITHREAD
	//	if(max_threads > 1) {
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x20;
			cc0   = sc_ptr + 0x21;
			ss0   = sc_ptr + 0x22;
			two   = sc_ptr + 0x43;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);
			  #ifdef USE_AVX2
				VEC_DBL_INIT(two  , 1.0);	// Yeah, I *know" "it's a misnomer" :)
				// cc0,ss0 inited below for AVX2
			  #else
				VEC_DBL_INIT(two  , 2.0);
				VEC_DBL_INIT(cc0  , c);
				VEC_DBL_INIT(ss0  , s);
			  #endif
				isrt2 += 72;	/* Move on to next thread's local store */
				cc0   += 72;
				ss0   += 72;
				two   += 72;
			}

		#elif defined(COMPILER_TYPE_GCC)

			r1  = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;
										cc0 = sc_ptr + 0x21;
										ss0 = sc_ptr + 0x22;
										two = sc_ptr + 0x43;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
		  #ifdef USE_AVX2
			VEC_DBL_INIT(two  , 1.0);
			// cc0,ss0 inited below for AVX2
		  #else
			VEC_DBL_INIT(two  , 2.0);
			VEC_DBL_INIT(cc0  , c);
			VEC_DBL_INIT(ss0  , s);
		  #endif

		#else
	//	} else {
																// In AVX2 mode, the 34 sincos terms between isrt2 and two will be inited
																// as before but will end up containing derived quantitiess as shown below.
			r1  = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;	// The derived-term naming reflects that in the RADIX_16_DIF_FMA macro in dft_macro.h:
			r2  = sc_ptr + 0x01;		cc0 = sc_ptr + 0x21;	// __c1_c = c1*__c
			r3  = sc_ptr + 0x02;		ss0 = sc_ptr + 0x22;	// __sc = __s/__c	[note this const across all blocks]
			r4  = sc_ptr + 0x03;		c0  = sc_ptr + 0x23;	// __c1i2 = c1*ISRT2
			r5  = sc_ptr + 0x04;		s0  = sc_ptr + 0x24;	// __c2i2 = c2*ISRT2
			r6  = sc_ptr + 0x05;		c8  = sc_ptr + 0x25;	// __c8 [unchanged]
			r7  = sc_ptr + 0x06;		s8  = sc_ptr + 0x26;	// __r8 = s8 /c8 
			r8  = sc_ptr + 0x07;		c4  = sc_ptr + 0x27;	// __c4 [unchanged]
			r9  = sc_ptr + 0x08;		s4  = sc_ptr + 0x28;	// __r4 = s4 /c4 
			r10 = sc_ptr + 0x09;		c12 = sc_ptr + 0x29;	// __cC4 = __cC/__c4
			r11 = sc_ptr + 0x0a;		s12 = sc_ptr + 0x2a;	// __rC = s12/c12
			r12 = sc_ptr + 0x0b;		c2  = sc_ptr + 0x2b;	// __c2 [unchanged]
			r13 = sc_ptr + 0x0c;		s2  = sc_ptr + 0x2c;	// __r2 = s2 /c2 
			r14 = sc_ptr + 0x0d;		c10 = sc_ptr + 0x2d;	// __cA2 = __cA/__c2
			r15 = sc_ptr + 0x0e;		s10 = sc_ptr + 0x2e;	// __rA = s10/c10
			r16 = sc_ptr + 0x0f;		c6  = sc_ptr + 0x2f;	// __c62 = __c6/__c2
			r17 = sc_ptr + 0x10;		s6  = sc_ptr + 0x30;	// __r6 = s6 /c6 
			r18 = sc_ptr + 0x11;		c14 = sc_ptr + 0x31;	// __cE6 = __cE/__c6
			r19 = sc_ptr + 0x12;		s14 = sc_ptr + 0x32;	// __rE = s14/c14
			r20 = sc_ptr + 0x13;		c1  = sc_ptr + 0x33;	// __c1 [unchanged]
			r21 = sc_ptr + 0x14;		s1  = sc_ptr + 0x34;	// __r1 = s1 /c1 
			r22 = sc_ptr + 0x15;		c9  = sc_ptr + 0x35;	// __c91 = __c9/__c1
			r23 = sc_ptr + 0x16;		s9  = sc_ptr + 0x36;	// __r9 = s9 /c9 
			r24 = sc_ptr + 0x17;		c5  = sc_ptr + 0x37;	// __c51 = __c5/__c1
			r25 = sc_ptr + 0x18;		s5  = sc_ptr + 0x38;	// __r5 = s5 /c5 
			r26 = sc_ptr + 0x19;		c13 = sc_ptr + 0x39;	// __cD5 = __cD/__c5
			r27 = sc_ptr + 0x1a;		s13 = sc_ptr + 0x3a;	// __rD = s13/c13
			r28 = sc_ptr + 0x1b;		c3  = sc_ptr + 0x3b;	// __c31 = __c3/__c1
			r29 = sc_ptr + 0x1c;		s3  = sc_ptr + 0x3c;	// __r3 = s3 /c3 
			r30 = sc_ptr + 0x1d;		c11 = sc_ptr + 0x3d;	// __cB3 = __cB/__c3
			r31 = sc_ptr + 0x1e;		s11 = sc_ptr + 0x3e;	// __rB = s11/c11
			r32 = sc_ptr + 0x1f;		c7  = sc_ptr + 0x3f;	// __c73 = __c7/__c3
										s7  = sc_ptr + 0x40;	// __r7 = s7 /c7 
										c15 = sc_ptr + 0x41;	// __cF7 = __cF/__c7
										s15 = sc_ptr + 0x42;	// __rF = s15/c15
										two = sc_ptr + 0x43;	// Holds 1.0 in AVX2 modeThese will need to be computed afresh for each new set of roots of unity, obviously.
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
		  #ifdef USE_AVX2
			#error AVX2 build not supported for non-GCC-compliant compilers!
		  #else
			VEC_DBL_INIT(two  , 2.0);
			VEC_DBL_INIT(cc0  , c);
			VEC_DBL_INIT(ss0  , s);
		  #endif
		#endif
	//	}
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		r1 = __r0 + thr_id*72;
		isrt2 = r1 + 0x20;
		cc0   = r1 + 0x21;
		two   = r1 + 0x43;
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
	// Make sure to at least one-time pre-test any index arithmetic assumptions used to save cycles in the loop
	// body (both C and AS). Since such checks may be runlength-dependent, need to be cheap enough to leave on
	// all the time, as here where we do them just once prior to entering the processing loop. Since DIF and DIT
	// encounter the same sets or index strides (albeit in opposite order), can split such tests between them:
//	ASSERT(HERE, p2  == p1+p1, "radix16_dif_pass: p2  != p1+p1!");	<*** The failure of this assertion led me to find the dependence on it in my new AVX2/FMA-based DIT macro ***

	iroot_prim=(incr >> 5);		/* (incr/2)/radix_now */

	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16).	 */
	{				/* NLOOPS satisfies the identity NLOOPS * INCR = P*N, which says that we process the entire */
					/* array each time this subroutine is executed (since P*N = vector length, sans padding.)   */

/*	here are the needed sincos data - these are processed below in bit-reversed order. */
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* In the DIF pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks.
*/
#if HIACC

	// In AVX2/FMA mode, since we need to replace most of the raw sincos data with derived ones,
	// simply place one copy of each computed double in a double-sized slot of the local memory.
	// We will be using AVX2/FMA-based Newtonian iterative inversion on the 16 doubles whose
	// multiplicative inverse is needed (the real part of the basic root of unity c and of the
	// 15 complex twiddles, c1-15), so store those in packed form in 4 AVX-register-sized
	// contiguous memory locations, and the others in a separate chunk of memory. After the
	// vector-iterative inversion we'll need to combine the 2 sets of data and place (in quadruplicate)
	// into their final SIMD-suitable memory slots.
	#ifdef USE_AVX2

		add0 = (double *)cc0;	// add0 points to 16 cos-data-to-be-inverted; Need a double-ptr on lhs here
		add1 = add0 + 16;	// add1 points to block of memory temporarily used to store the corresponding sine data
		add2 = add0 + 32;	// add2 points to block of memory temporarily used to store the 11 [0-padded to 12]
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c1, for inversion
		*add1++ = it;	// s1  slot will hold __r1 = s1 /c1 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c2, for inversion
		*add1++ = it;	// s2  slot will hold __r2 = s2 /c2 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*(add0-3) = *add0++ = rt;	// c3, for inversion ... place extra copy in 0-slot as described above
		*add1++ = it;	// s3  slot will hold __r3 = s3 /c3 
		*add2++ = rt;	// c3, will get multiplied by 1/c1 to yield __c31
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c4, for inversion
		*add1++ = it;	// s4  slot will hold __r4 = s4 /c4 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c5, for inversion
		*add1++ = it;	// s5  slot will hold __r5 = s5 /c5 
		*add2++ = rt;	// c5, will get multiplied by 1/c1 to yield __c51
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c6, for inversion
		*add1++ = it;	// s6  slot will hold __r6 = s6 /c6 
		*add2++ = rt;	// c6, will get multiplied by 1/c2 to yield __c62
		*add2++ = 0.0;	// 0-pad will get multiplied by 1/c3 term, remains 0-pad.
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c7, for inversion
		*add1++ = it;	// s7  slot will hold __r7 = s7 /c7 
		*add2++ = rt;	// c7, will get multiplied by 1/c3 to yield __c73
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c8, for inversion
		*add1++ = it;	// s8  slot will hold __r8 = s8 /c8 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c9, for inversion
		*add1++ = it;	// s9  slot will hold __r9 = s9 /c9 
		*add2++ = rt;	// c9, will get multiplied by 1/c1 to yield __c91
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c10, for inversion
		*add1++ = it;	// s10 slot will hold __rA = s10/c10
		*add2++ = rt;	// cA, will get multiplied by 1/c2 to yield __cA2
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c11, for inversion
		*add1++ = it;	// s11 slot will hold __rB = s11/c11
		*add2++ = rt;	// cB, will get multiplied by 1/c3 to yield __cB3
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c12, for inversion
		*add1++ = it;	// s12 slot will hold __rC = s12/c12
		*add2++ = rt;	// cC, will get multiplied by 1/c4 to yield __cC4
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c13, for inversion
		*add1++ = it;	// s13 slot will hold __rD = s13/c13
		*add2++ = rt;	// cD, will get multiplied by 1/c5 to yield __cD5
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c14, for inversion
		*add1++ = it;	// s14 slot will hold __rE = s14/c14
		*add2++ = rt;	// cE, will get multiplied by 1/c6 to yield __cE6

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c15, for inversion
		*add1++ = it;	// s15 slot will hold __rF = s15/c15
		*add2++ = rt;	// cF, will get multiplied by 1/c7 to yield __cF7

		// This places us at add0 == c8 and add1 = c12.
		ASSERT(HERE, add0 == (double *)cc0+16 && add1 == (double *)cc0+32 && add2 == (double *)cc0+44, "add0,1,2 checksum failed in AVX2 sincos inits!");
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

		RADIX16_COMPUTE_FMA_SINCOS_DIF_2(cc0,two);
		add0 = (double *)cc0;

		/* Scalar data starting at add0 = cc0 laid out as below:

		a[0x00,0x01,0x02,0x03]: add0 + 0x[  0,  8, 10, 18]: [c3 ,c1 ,c2 ,c3 ] Cosines:
		a[0x04,0x05,0x06,0x07]: add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ]
		a[0x08,0x09,0x0a,0x0b]: add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ]
		a[0x0c,0x0d,0x0e,0x0f]: add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ]
		a[0x10,0x11,0x12,0x13]: add0 + 0x[ 80, 88, 90, 98]: [-- ,r1 ,r2 ,r3 ] Tangents:
		a[0x14,0x15,0x16,0x17]: add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ]
		a[0x18,0x19,0x1a,0x1b]: add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ]
		a[0x1c,0x1d,0x1e,0x1f]: add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ]
		a[0x20,0x21,0x22,0x23]: add0 + 0x[100,108,110,118]: [c31,c51,c62,  0] Cosine ratios:
		a[0x24,0x25,0x26,0x27]: add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3]
		a[0x28,0x29,0x2a,0x2b]: add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7]

		Ensuing C code massages the above into a scalar-data analog of the 4-copy layout.
		*/
		// put all overwrites-of-no-longer-needed data first, to minimize conflicts later.
		// Data which will not be used in the FMA-based radix-16 DIF are at indices 0x[0,3,5-7,9-15,16]:
		// Arrange the rest so RHS (read-elt) indices are ascending, then manually move
		// up those which overwrite indices appearing further down in the RHS:
		add0[0x00] = add0[0x01];	// c1, copy to *= __c
		add0[0x03] = add0[0x02];	// c2, copy to *= ISRT2
		add0[0x0a] = add0[0x02];	// copy c2 to final loc before overwriting
		add0[0x02] = add0[0x01];	// c1, copy to *= ISRT2
		add0[0x01] = tan;
		add0[0x06] = add0[0x04];	// c4 	
		add0[0x04] = add0[0x08];	// c8	
		add0[0x05] = add0[0x18];	// r8
		add0[0x07] = add0[0x14];	// r4
		add0[0x08] = add0[0x28];	// cC4	
		add0[0x09] = add0[0x1c];	// rC
		add0[0x0b] = add0[0x12];	// r2
		add0[0x0c] = add0[0x26];	// cA2	
		add0[0x0d] = add0[0x1a];	// rA
		add0[0x0e] = add0[0x22];	// c62	
		add0[0x0f] = add0[0x16];	// r6
		add0[0x10] = add0[0x2a];	// cE6			
		add0[0x16] = add0[0x21];	// c51	
		add0[0x21] = add0[0x1f];	// rF
		add0[0x1f] = add0[0x17];	// r7
		add0[0x17] = add0[0x15];	// r5
		add0[0x15] = add0[0x19];	// r9
		add0[0x19] = add0[0x1d];	// rD
		add0[0x1d] = add0[0x1b];	// rB
		add0[0x1b] = add0[0x13];	// r3
		add0[0x13] = add0[0x11];	// r1
		add0[0x11] = add0[0x1e];	// rE
		add0[0x12] = add0[0x00];	// c1 now in [0]
		add0[0x14] = add0[0x25];	// c91	
		add0[0x18] = add0[0x29];	// cD5	
		add0[0x1a] = add0[0x20];	// c31	
		add0[0x1c] = add0[0x27];	// cB3	
		add0[0x1e] = add0[0x24];	// c73	
		add0[0x20] = add0[0x2b];	// cF7	
		// Now mpy elts in slots 0,2,3 by __c, ISRT2, ISRT2, respectively:
		add0[0x00] *= c;
		add0[0x02] *= ISRT2;
		add0[0x03] *= ISRT2;
		// And stick a 1.0 at the end of the above block-of-doubles:
		add0[0x22] = 1.0;

		/* Yielding the following layout-of-scalar-doubles, with data above index 0x22 unused in the DIF DFT:

		a[0x00,0x01,0x02,0x03]: add0 + 0x[  0,  8, 10, 18]: [c1*c,s/c,c1*ISRT2,c2*ISRT2]
		a[0x04,0x05,0x06,0x07]: add0 + 0x[ 20, 28, 30, 38]: [c8 ,r8 ,c4 ,r4 ]
		a[0x08,0x09,0x0a,0x0b]: add0 + 0x[ 40, 48, 50, 58]: [cC4,rC ,c2 ,r2 ]
		a[0x0c,0x0d,0x0e,0x0f]: add0 + 0x[ 60, 68, 70, 78]: [cA2,rA ,c62,r6 ]
		a[0x10,0x11,0x12,0x13]: add0 + 0x[ 80, 88, 90, 98]: [cE6,rE ,c1 ,r1 ]
		a[0x14,0x15,0x16,0x17]: add0 + 0x[ a0, a8, b0, b8]: [c91,r9 ,c51,r5 ]
		a[0x18,0x19,0x1a,0x1b]: add0 + 0x[ c0, c8, d0, d8]: [cD5,rD ,c31,r3 ]
		a[0x1c,0x1d,0x1e,0x1f]: add0 + 0x[ e0, e8, f0, f8]: [cB3,rB ,c73,r7 ]
		a[0x20,0x21,0x22,0x23]: add0 + 0x[100,108,110,118]: [cF7,rF ,1.0,  0]
		a[0x24,0x25,0x26,0x27]: add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3]
		a[0x28,0x29,0x2a,0x2b]: add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7]
		*/
	  #else	// Make DFT_V1 the default here:

		add0 = &c; add1 = &tan;	// GCC/Clang don't allow address-taking inlined in arglist of macros, so do it here
		RADIX16_COMPUTE_FMA_SINCOS_DIF(cc0,two,add0,add1);

	  #endif

	#else	// AVX2 = False:

	  #ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-15] are offset w.r.to the thread-local ptr pair as
		(cc0,ss0) + 0x[2,12,a,1a,6,16,e,1e,4,14,c,1c,8,18,10,20]:
		*/
		c_tmp = cc0 + 0x02; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c1=t1*t3-t2*t4;	s1=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i;				/* 4*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c2=t1*t3-t2*t4;	s2=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i;				/* 8*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c4=t1*t3-t2*t4;	s4=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += (iroot << 2)+iroot;		/* 13*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c8=t1*t3-t2*t4;	s8=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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
		j2 = j1 + RE_IM_STRIDE;
		pfetch_addr = p1*((j >> l2_stride) & 0x3);	// cycle prefetch-offset-address among p0,1,2,3
			// These get added to the base addresses p0,4,8,12 in thr DFT macros, thus every 4 loop
			// executions we have covered prefetches from [current address] + [pfetch distance] + p0,1,2,...15 .

#ifdef USE_SSE2

	#ifdef COMPILER_TYPE_MSVC

	  #if 1	// if(1) - test out pure-asm version

		add0 = &a[j1];
		__asm	mov	eax, add0
		__asm	mov	ebx, p4		// Can't get these via simple load-one-and-shift-as-needed due to array padding scheme
		__asm	mov	ecx, p8
		__asm	mov	edx, p12
		__asm	shl	ebx,  3
		__asm	shl	ecx,  3
		__asm	shl	edx,  3
		__asm	add	ebx, eax
		__asm	add	ecx, eax
		__asm	add	edx, eax
		SSE2_RADIX4_DIF_4TWIDDLE_B(r1 ,c0)
		__asm	mov	edi, p2		// Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p1+p1 not guaranteed.
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p2];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_B(r9 ,c2)
		__asm	sub	eax, edi	// &a[j1];
		__asm	sub	ebx, edi
		__asm	sub	ecx, edi
		__asm	sub	edx, edi
		__asm	mov	edi, p1
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p1];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_B(r17,c1)
		__asm	sub	eax, edi	// &a[j1];
		__asm	sub	ebx, edi
		__asm	sub	ecx, edi
		__asm	sub	edx, edi
		__asm	mov	edi, p3
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p3];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_B(r25,c3)

	  #else

	/*...Block 1: */
		add0 = &a[j1];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;

		/* Do the p4,12 combo first. 	Cost: 24 MOVapd, 28 ADD/SUBpd, 12 MULpd */
		__asm	mov	eax, add1
		__asm	mov	ebx, c4
		__asm	mov	ecx, add3
		__asm	mov	edx, c12

		__asm	movaps	xmm0,[eax     ]	/* a[jt+p4] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p12] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp+p4] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p12] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt+p4] */	__asm	movaps	xmm6,[ecx     ]	/* xmm6 <- cpy a[jt+p12] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp+p4] */	__asm	movaps	xmm7,[ecx+0x10]	/* xmm7 <- cpy a[jp+p12] */

		__asm	mulpd	xmm0,[ebx     ]	/* a[jt+p4]*c4 */			__asm	mulpd	xmm4,[edx     ]	/* a[jt+p12]*c12 */
		__asm	mulpd	xmm1,[ebx     ]	/* a[jp+p4]*c4 */			__asm	mulpd	xmm5,[edx     ]	/* a[jp+p12]*c12 */
		__asm	mulpd	xmm2,[ebx+0x10]	/* a[jt+p4]*s4 */			__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p12]*s12 */
		__asm	mulpd	xmm3,[ebx+0x10]	/* a[jp+p4]*s4 */			__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p12]*s12 */
		__asm	addpd	xmm1,xmm2	/* xmm1 <- t6 */				__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm0,xmm3	/* xmm0 <- t5 */				__asm	subpd	xmm4,xmm7	/* xmm4 <- rt 	xmm6,7 free */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t6 */
		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t5 */

		__asm	addpd	xmm0,xmm4	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm1,xmm5	/* ~t6 <- t6 +it */
		__asm	subpd	xmm2,xmm4	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm3,xmm5	/* ~t8 <- t6 -it	xmm4,5 free */

		/* Now do the p0,8 combo: */
		__asm	mov	eax, add0
		__asm	mov	ecx, add2
		__asm	mov	edx, c8

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm6 <- cpy a[jt+p8] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm7 <- cpy a[jp+p8] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p8]*c8 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p8]*c8 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p8]*s8 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p8]*s8 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt 	xmm6,7 free - stick t1,t2 in those */

		__asm	movaps	xmm6,[eax     ]	/* t1 = a[jt   ] */
		__asm	movaps	xmm7,[eax+0x10]	/* t2 = a[jp   ] */

		__asm	subpd	xmm6,xmm4	/* ~t3 <- t1 -rt */
		__asm	subpd	xmm7,xmm5	/* ~t4 <- t2 -it */
		__asm	addpd	xmm4,xmm4	/*          2*rt */
		__asm	addpd	xmm5,xmm5	/*          2*it */
		__asm	addpd	xmm4,xmm6	/* ~t1 <- t1 +rt */
		__asm	addpd	xmm5,xmm7	/* ~t2 <- t2 +it	xmm4,5 free */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		__asm	mov	eax, r1
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;
		~t6 =t2 -t6;		~t2 =t2 +t6;
		*/
		__asm	subpd	xmm4,xmm0	/*~t5 =t1 -t5 */
		__asm	subpd	xmm5,xmm1	/*~t6 =t2 -t6 */
		__asm	movaps	[eax+0x040],xmm4	/* a[jt+p8 ] <- ~t5 */
		__asm	movaps	[eax+0x050],xmm5	/* a[jp+p8 ] <- ~t6 */
		__asm	addpd	xmm0,xmm0	/* 2*t5 */
		__asm	addpd	xmm1,xmm1	/* 2*t6 */
		__asm	addpd	xmm0,xmm4	/*~t1 =t1 +t5 */
		__asm	addpd	xmm1,xmm5	/*~t2 =t2 +t6 */
		__asm	movaps	[eax      ],xmm0	/* a[jt    ] <- ~t1 */
		__asm	movaps	[eax+0x010],xmm1	/* a[jp    ] <- ~t2 */

		/*
		~t7 =t3 +t8;		~t3 =t3 -t8;
		~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	subpd	xmm6,xmm3	/*~t3 =t3 -t8 */
		__asm	subpd	xmm7,xmm2	/*~t8 =t4 -t7 */
		__asm	movaps	[eax+0x020],xmm6	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[eax+0x070],xmm7	/* a[jp+p12] <- ~t8 */
		__asm	addpd	xmm3,xmm3	/* 2*t8 */
		__asm	addpd	xmm2,xmm2	/* 2*t7 */
		__asm	addpd	xmm3,xmm6	/*~t7 =t3 +t8 */
		__asm	addpd	xmm2,xmm7	/*~t4 =t4 +t7 */
		__asm	movaps	[eax+0x060],xmm3	/* a[jt+p12] <- ~t7 */
		__asm	movaps	[eax+0x030],xmm2	/* a[jp+p4 ] <- ~t4 */

	/*...Block 2: */
		add0 = &a[j1+p2];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;
		/* 	Cost: 30 MOVapd, 28 ADD/SUBpd, 16 MULpd */
		SSE2_RADIX4_DIF_4TWIDDLE(add0,add1,add2,add3,r9 ,c2)

	/*...Block 3: */
		add0 = &a[j1+p1];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;
		/* 	Cost: 30 MOVapd, 28 ADD/SUBpd, 16 MULpd */
		SSE2_RADIX4_DIF_4TWIDDLE(add0,add1,add2,add3,r17,c1 )

	/*...Block 4: */
		add0 = &a[j1+p3];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;
		/* 	Cost: 30 MOVapd, 28 ADD/SUBpd, 16 MULpd */
		SSE2_RADIX4_DIF_4TWIDDLE(add0,add1,add2,add3,r25,c3)

	  #endif	/* if(1) */

	/**************************************************************************************
	!...and now do four more radix-4 transforms, including the internal twiddle factors:  !
	**************************************************************************************/
	/*...Block 1: t1,9,17,25	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
	  #if 1
		__asm	mov eax, add0	// &a[j1]
		__asm	mov ebx, p1
		__asm	mov ecx, p2
		__asm	mov edx, p3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
	  #else
		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	  #endif
		__asm	mov	edi, r1
		__asm	movaps	xmm0,[edi      ]	/* t1  */			__asm	movaps	xmm4,[edi+0x100]	/* t17 */
		__asm	movaps	xmm1,[edi+0x010]	/* t2  */			__asm	movaps	xmm5,[edi+0x110]	/* t18 */
		__asm	movaps	xmm2,[edi+0x080]	/* t9  */			__asm	movaps	xmm6,[edi+0x180]	/* t25 */
		__asm	movaps	xmm3,[edi+0x090]	/* t10 */			__asm	movaps	xmm7,[edi+0x190]	/* t26 */

		__asm	subpd	xmm0,xmm2		/* ~t9 = t1-t9 */		__asm	subpd	xmm4,xmm6		/* ~t25=t17-t25*/
		__asm	subpd	xmm1,xmm3		/* ~t10= t2-t10*/		__asm	subpd	xmm5,xmm7		/* ~t26=t18-t26*/
		__asm	addpd	xmm2,xmm2		/*        2*t9 */		__asm	addpd	xmm6,xmm6		/*        2*t25*/
		__asm	addpd	xmm3,xmm3		/*        2*t10*/		__asm	addpd	xmm7,xmm7		/*        2*t26*/
		__asm	addpd	xmm2,xmm0		/* ~t1 = t1+t9 */		__asm	addpd	xmm6,xmm4		/* ~t17=t17+t25*/
		__asm	addpd	xmm3,xmm1		/* ~t2 = t2+t10*/		__asm	addpd	xmm7,xmm5		/* ~t18=t18+t26*/

		__asm	subpd	xmm2,xmm6		/* t1  <- t1 -t17 */
		__asm	subpd	xmm3,xmm7		/* t2  <- t2 -t18 */
		__asm	addpd	xmm6,xmm6		/*          2*t17 */
		__asm	addpd	xmm7,xmm7		/*          2*t18 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */
		__asm	addpd	xmm6,xmm2		/* t17 <- t1 +t17 */
		__asm	addpd	xmm7,xmm3		/* t18 <- t2 +t18 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */

		__asm	subpd	xmm0,xmm5		/* t9  <- t9 -t26 */
		__asm	subpd	xmm1,xmm4		/* t10 <- t10-t25 */
		__asm	addpd	xmm5,xmm5		/*          2*t26 */
		__asm	addpd	xmm4,xmm4		/*          2*t25 */
		__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm5,xmm0		/* t26 <- t9 +t26 */
		__asm	addpd	xmm4,xmm1		/* t25 <- t10+t25 */
		__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 3: t5,13,21,29	Cost: 16 MOVapd, 26 ADD/SUBpd,  4 MULpd */
	  #if 1
		__asm	mov	edi, p4
		__asm	shl	edi, 3
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	  #else
		add0 = &a[j1+p4];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	  #endif
		__asm	mov	edi, r5
		__asm	mov	esi, isrt2
		__asm	movaps	xmm3,[esi ]	/* isrt2 */
																__asm	movaps	xmm4,[edi+0x100]	/* t21 */
																__asm	movaps	xmm5,[edi+0x110]	/* t22 */
																__asm	movaps	xmm6,[edi+0x180]	/* t29 */
																__asm	movaps	xmm7,[edi+0x190]	/* t30 */
																__asm	mulpd	xmm4,xmm3	/* t21 *ISRT2 */
		__asm	movaps	xmm0,[edi      ]	/* t5  */			__asm	mulpd	xmm5,xmm3	/* t22 *ISRT2 */
		__asm	movaps	xmm1,[edi+0x010]	/* t6  */			__asm	mulpd	xmm6,xmm3	/* t29 *ISRT2 */
		__asm	movaps	xmm2,[edi+0x080]	/* t13 */			__asm	mulpd	xmm7,xmm3	/* t30 *ISRT2 */
		__asm	movaps	xmm3,[edi+0x090]	/* t14; this must execute after the last mul-by-ISRT2 above */

		__asm	subpd	xmm0,xmm3		/* ~t5 = t5 -t14*/		__asm	subpd	xmm4,xmm5		/* ~t21=t21-t22*/
		__asm	subpd	xmm1,xmm2		/* ~t14= t6 -t13*/		__asm	subpd	xmm7,xmm6		/*  it =t30-t29*/
		__asm	addpd	xmm3,xmm3		/*         2*t14*/		__asm	addpd	xmm5,xmm5		/*        2*t22*/
		__asm	addpd	xmm2,xmm2		/*         2*t13*/		__asm	addpd	xmm6,xmm6		/*        2*t29*/
		__asm	addpd	xmm3,xmm0		/* ~t13= t14+t5 */		__asm	addpd	xmm5,xmm4		/* ~t22=t22+t21*/
		__asm	addpd	xmm2,xmm1		/* ~t6 = t13+t6 */		__asm	addpd	xmm6,xmm7		/*  rt =t29+t30*/

		__asm	subpd	xmm4,xmm6		/* t21=t21-rt */
		__asm	subpd	xmm5,xmm7		/* t22=t22-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/* t29=t21+rt */
		__asm	addpd	xmm7,xmm5		/* t30=t22+it */

		__asm	subpd	xmm0,xmm4		/* t5 -t21 */			__asm	subpd	xmm3,xmm7		/* t13-t30 */
		__asm	subpd	xmm2,xmm5		/* t6 -t22 */			__asm	subpd	xmm1,xmm6		/* t14-t29 */
		__asm	addpd	xmm4,xmm4		/*   2*t21 */			__asm	addpd	xmm7,xmm7		/*   2*t30 */
		__asm	addpd	xmm5,xmm5		/*   2*t22 */			__asm	addpd	xmm6,xmm6		/*   2*t29 */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm3	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm2	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm0		/* t5 +t21 */			__asm	addpd	xmm7,xmm3		/* t13+t30 */
		__asm	addpd	xmm5,xmm2		/* t6 +t22 */			__asm	addpd	xmm6,xmm1		/* t14+t29 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

	/*...Block 2: t3,11,19,27	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
	  #if 1
		__asm	mov	edi, p4
		__asm	shl	edi, 3
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	  #else
		add0 = &a[j1+p8];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	  #endif
		__asm	mov	edi, r3
		__asm	mov	esi, cc0
		__asm	movaps	xmm4,[edi+0x100]	/* t19 */
		__asm	movaps	xmm5,[edi+0x110]	/* t20 */
		__asm	movaps	xmm3,[esi      ]	/* c */
		__asm	movaps	xmm2,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4		/* copy t19 */
		__asm	movaps	xmm7,xmm5		/* copy t20 */

		__asm	mulpd	xmm4,xmm3		/* t19*c */
		__asm	mulpd	xmm5,xmm3		/* t20*c */
		__asm	mulpd	xmm6,xmm2		/* t19*s */				__asm	movaps	xmm0,[edi+0x180]	/* t27 */
		__asm	mulpd	xmm7,xmm2		/* t20*s */				__asm	movaps	xmm1,[edi+0x190]	/* t28 */
		__asm	addpd	xmm5,xmm6	/* ~t20 */					__asm	movaps	xmm6,xmm0		/* copy t27 */
		__asm	subpd	xmm4,xmm7	/* ~t19 */					__asm	movaps	xmm7,xmm1		/* copy t28 */

																__asm	mulpd	xmm6,xmm2		/* t27*s */
																__asm	mulpd	xmm7,xmm2		/* t28*s */
																__asm	mulpd	xmm0,xmm3		/* t27*c */
																__asm	mulpd	xmm1,xmm3		/* t28*c */
																__asm	addpd	xmm7,xmm0	/* it */
																__asm	subpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t19 */
		__asm	movaps	xmm3,xmm5		/* copy t20 */
		__asm	subpd	xmm4,xmm6		/*~t27=t19-rt */
		__asm	subpd	xmm5,xmm7		/*~t28=t20-it */
		__asm	addpd	xmm6,xmm2		/*~t19=t19+rt */
		__asm	addpd	xmm7,xmm3		/*~t20=t20+it */

		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edi+0x080]	/* t11 */
		__asm	movaps	xmm3,[edi+0x090]	/* t12 */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm2	/* cpy t11 */
		__asm	subpd	xmm2,xmm3	/*~t11=t11-t12 */
		__asm	addpd	xmm3,xmm0	/*~t12=t12+t11 */
		__asm	mulpd	xmm2,xmm1	/* rt = (t11-t12)*ISRT2 */
		__asm	mulpd	xmm3,xmm1	/* it = (t12+t11)*ISRT2 */

		__asm	movaps	xmm0,[edi      ]	/* t3  */
		__asm	movaps	xmm1,[edi+0x010]	/* t4  */

		__asm	subpd	xmm0,xmm2			/*~t11=t3 -rt */
		__asm	subpd	xmm1,xmm3			/*~t12=t4 -it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t3 =rt +t3 */
		__asm	addpd	xmm3,xmm1			/*~t4 =it +t4 */

		__asm	subpd	xmm2,xmm6		/* t3 -t19 */			__asm	subpd	xmm0,xmm5		/* t11-t28 */
		__asm	subpd	xmm3,xmm7		/* t4 -t20 */			__asm	subpd	xmm1,xmm4		/* t12-t27 */
		__asm	addpd	xmm6,xmm6		/*   2*t19 */			__asm	addpd	xmm5,xmm5		/*          2*t28 */
		__asm	addpd	xmm7,xmm7		/*   2*t20 */			__asm	addpd	xmm4,xmm4		/*          2*t27 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm6,xmm2		/* t3 +t19 */			__asm	addpd	xmm5,xmm0		/* t11+t28 */
		__asm	addpd	xmm7,xmm3		/* t4 +t20 */			__asm	addpd	xmm4,xmm1		/* t12+t27 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 4: t7,15,23,31	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
	  #if 1
		__asm	mov	edi, p4
		__asm	shl	edi, 3
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	  #else
		add0 = &a[j1+p12];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	  #endif
		__asm	mov	edi, r7
		__asm	mov	esi, cc0
		__asm	movaps	xmm4,[edi+0x100]	/* t23 */
		__asm	movaps	xmm5,[edi+0x110]	/* t24 */
		__asm	movaps	xmm2,[esi      ]	/* c */
		__asm	movaps	xmm3,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4		/* copy t23 */
		__asm	movaps	xmm7,xmm5		/* copy t24 */

		__asm	mulpd	xmm4,xmm3		/* t23*s */
		__asm	mulpd	xmm5,xmm3		/* t24*s */
		__asm	mulpd	xmm6,xmm2		/* t23*c */				__asm	movaps	xmm0,[edi+0x180]	/* t31 */
		__asm	mulpd	xmm7,xmm2		/* t24*c */				__asm	movaps	xmm1,[edi+0x190]	/* t32 */
		__asm	addpd	xmm5,xmm6	/* ~t24 */					__asm	movaps	xmm6,xmm0		/* copy t31 */
		__asm	subpd	xmm4,xmm7	/* ~t23 */					__asm	movaps	xmm7,xmm1		/* copy t32 */

																__asm	mulpd	xmm6,xmm2		/* t31*c */
																__asm	mulpd	xmm7,xmm2		/* t32*c */
																__asm	mulpd	xmm0,xmm3		/* t31*s */
																__asm	mulpd	xmm1,xmm3		/* t32*s */
																__asm	addpd	xmm7,xmm0	/* it */
																__asm	subpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t23 */
		__asm	movaps	xmm3,xmm5		/* copy t24 */
		__asm	subpd	xmm4,xmm6		/*~t23=t23-rt */
		__asm	subpd	xmm5,xmm7		/*~t24=t24-it */
		__asm	addpd	xmm6,xmm2		/*~t31=t23+rt */
		__asm	addpd	xmm7,xmm3		/*~t32=t24+it */

		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edi+0x080]	/* t15 */
		__asm	movaps	xmm3,[edi+0x090]	/* t16 */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm2	/* cpy t15 */
		__asm	addpd	xmm2,xmm3	/*~t15=t15+t16 */
		__asm	subpd	xmm3,xmm0	/*~t16=t16-t15 */
		__asm	mulpd	xmm2,xmm1	/* rt = (t15+t16)*ISRT2 */
		__asm	mulpd	xmm3,xmm1	/* it = (t16-t15)*ISRT2 */

		__asm	movaps	xmm0,[edi      ]	/* t7  */
		__asm	movaps	xmm1,[edi+0x010]	/* t8  */

		__asm	subpd	xmm0,xmm2			/*~t7 =t7 -rt */
		__asm	subpd	xmm1,xmm3			/*~t8 =t8 -it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t15=rt +t7 */
		__asm	addpd	xmm3,xmm1			/*~t16=it +t8 */

		__asm	subpd	xmm0,xmm4		/* t7 -t23 */			__asm	subpd	xmm2,xmm7		/* t15-t32 */
		__asm	subpd	xmm1,xmm5		/* t8 -t24 */			__asm	subpd	xmm3,xmm6		/* t16-t31 */
		__asm	addpd	xmm4,xmm4		/*   2*t23 */			__asm	addpd	xmm7,xmm7		/*   2*t32 */
		__asm	addpd	xmm5,xmm5		/*   2*t24 */			__asm	addpd	xmm6,xmm6		/*   2*t31 */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm2	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm1	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm3	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm0		/* t7 +t23 */			__asm	addpd	xmm7,xmm2		/* t15+t32 */
		__asm	addpd	xmm5,xmm1		/* t8 +t24 */			__asm	addpd	xmm6,xmm3		/* t16+t31 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

		/***************************************************/
		/* Total Cost: 182 MOVapd, 214 ADD/SUBpd, 84 MULpd */
		/***************************************************/

	#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		add0 = &a[j1];

	  #ifdef USE_AVX2

	   #ifdef DFT_V2	// Toggle between 2 versions

		SSE2_RADIX16_DIF_TWIDDLE_2(add0,p1,p2,p3,p4,p8,p12,r1,  cc0,pfetch_addr,pfetch_dist);

	   #else

		SSE2_RADIX16_DIF_TWIDDLE  (add0,p1,p2,p3,p4,p8,p12,r1,isrt2,pfetch_addr,pfetch_dist);

	   #endif

	  // Play with prefetch in 64-bit AVX/SSE2 versions - once find reliable pfetch scheme, propagate to 32-bit SSE2 and AVX2 macros:
	  #elif OS_BITS == 64
		SSE2_RADIX16_DIF_TWIDDLE(add0,p1,p2,p3,p4,p8,p12,r1,isrt2,pfetch_addr,pfetch_dist);
	  #else	// 32-bit SSE2:
		SSE2_RADIX16_DIF_TWIDDLE(add0,p1,p2,p3,p4,p8,p12,r1,isrt2);
	  #endif

	#endif

#else	/* USE_SSE2 */

  #ifdef USE_SCALAR_DFT_MACRO	// Must define - or not - @compile time

		// Test FMA-based DIF macro:
		RADIX_16_DIF_FMA(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p4+p1 ],a[j2+p4+p1 ],a[j1+p4+p2 ],a[j2+p4+p2 ],a[j1+p4+p3 ],a[j2+p4+p3 ],a[j1+p8 ],a[j2+p8 ],a[j1+p8+p1 ],a[j2+p8+p1 ],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+p12],a[j2+p12],a[j1+p12+p1],a[j2+p12+p1],a[j1+p12+p2],a[j2+p12+p2],a[j1+p12+p3],a[j2+p12+p3]
						,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p4+p1 ],a[j2+p4+p1 ],a[j1+p4+p2 ],a[j2+p4+p2 ],a[j1+p4+p3 ],a[j2+p4+p3 ],a[j1+p8 ],a[j2+p8 ],a[j1+p8+p1 ],a[j2+p8+p1 ],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+p12],a[j2+p12],a[j1+p12+p1],a[j2+p12+p1],a[j1+p12+p2],a[j2+p12+p2],a[j1+p12+p3],a[j2+p12+p3]
						,c1,s1,c2,s2,c3,s3,c4,s4,c5,s5,c6,s6,c7,s7,c8,s8,c9,s9,c10,s10,c11,s11,c12,s12,c13,s13,c14,s14,c15,s15, c,s)

  #else		// USE_SCALAR_DFT_MACRO = False

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/
	#ifdef PFETCH_AGGRESSIVE
		addp = &a[j1];
	#elif PFETCH
		prefetch_offset = ((j >> 1) & 3)*p4 + 4;	/* Cycle among p0, p4, p8 and p12. */
	#endif
	/*...Block 1: */

		jt = j1;
		jp = j2;

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

		jt = j1 + p2;
		jp = j2 + p2;

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

		jt = j1 + p1;
		jp = j2 + p1;

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

		jt = j1 + p3;
		jp = j2 + p3;

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
	!	1, exp(i* 1*twopi/16) =		( c, s), exp(i* 2*twopi/16) = isqrt2*( 1, 1), exp(i* 3*twopi/16) =		( s, c) (for inputs to transform block 2)
	!	1, exp(i* 2*twopi/16) = isqrt2*( 1, 1), exp(i* 4*twopi/16) =		( 0, 1), exp(i* 6*twopi/16) = isqrt2*(-1, 1) (for inputs to transform block 3)
	!	1, exp(i* 3*twopi/16) =		( s, c), exp(i* 6*twopi/16) = isqrt2*(-1, 1), exp(i* 9*twopi/16) =		(-c,-s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1,2+p0:15:1) are replaced by t0:30:2,
	!										   a[j1,2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */

		jt = j1;
		jp = j2;

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

		jt = j1 + p4;
		jp = j2 + p4;

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

		jt = j1 + p8;
		jp = j2 + p8;

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

		jt = j1 + p12;
		jp = j2 + p12;

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
void radix16_dit_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
{
	const int RADIX = 16;
	const int pfetch_dist = PFETCH_DIST;
	int pfetch_addr;
	static int max_threads = 0;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	// lg(stride):
  #ifdef USE_AVX
	const int l2_stride = 3;	// 8 doubles at a time
  #elif defined(USE_SSE2)
	const int l2_stride = 2;	// 4 doubles at a time
  #else
	const int l2_stride = 1;	// 2 doubles at a time in scalar [non-SIMD] mode
  #endif

	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]*/
#ifdef USE_AVX2	// FMA-based DFT in AVX2 mode needs the tangent
	static double tan = 0.41421356237309504879;
#endif
	int i,j,j1,j2,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	int p1,p2,p3,p4,p8,p12;
	double rt,it,dtmp;
	double re0,im0,re1,im1;
	uint64 tmp64;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0, *add1;	/* Addresses into array sections */
#ifdef COMPILER_TYPE_MSVC
	double *add2, *add3;
#endif
	vec_dbl *c_tmp,*s_tmp;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
	vec_dbl *cc0, *ss0, *isrt2, *two, *r1;	// In || mode, only above base-pointer (shared by all threads) is static:
  #elif defined(COMPILER_TYPE_GCC)
	static vec_dbl *cc0, *ss0, *isrt2, *two, *r1;
  #else
	static vec_dbl *cc0, *ss0, *isrt2, *two;
	static vec_dbl *c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15,*s0,*s1,*s2,*s3,*s4,*s5,*s6,*s7,*s8,*s9,*s10,*s11,*s12,*s13,*s14,*s15;
	static vec_dbl *r1,*r2,*r3,*r4,*r5,*r6,*r7,*r8,*r9,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r30,*r31,*r32;
  #endif

#else

	int jp,jt;
	double *addr, *addp;
	int prefetch_offset;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
	,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;

#endif

#ifndef COMPILER_TYPE_GCC
	ASSERT(HERE, NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
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
		ASSERT(HERE, NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(HERE, max_threads >= NTHREADS, "Multithreading requires max_threads >= NTHREADS!");
		ASSERT(HERE, thr_id == -1, "Init-mode call must be outside of any multithreading!");
		sc_arr = ALLOC_VEC_DBL(sc_arr, 72*max_threads);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	last 30 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
		#ifdef MULTITHREAD
	//	if(max_threads > 1) {
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x20;
			cc0   = sc_ptr + 0x21;
			ss0   = sc_ptr + 0x22;
			two   = sc_ptr + 0x43;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);
			  #if defined(USE_AVX2) && (defined(DFT_V1) || defined(DFT_V2))
				VEC_DBL_INIT(two  , 1.0);	// Yeah, I *know" "it's a misnomer" :)
				// cc0,ss0 inited below for AVX2
			  #else
				VEC_DBL_INIT(two  , 2.0);
				VEC_DBL_INIT(cc0  , c);
				VEC_DBL_INIT(ss0  , s);
			  #endif
				isrt2 += 72;	/* Move on to next thread's local store */
				cc0   += 72;
				ss0   += 72;
				two   += 72;
			}
		#elif defined(COMPILER_TYPE_GCC)
			r1  = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;
										cc0 = sc_ptr + 0x21;
										ss0 = sc_ptr + 0x22;
										two = sc_ptr + 0x43;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
		  #ifdef USE_AVX2
			VEC_DBL_INIT(two  , 1.0);	// Yeah, I *know" "it's a misnomer" :)
			// cc0,ss0 inited below for AVX2
		  #else
			VEC_DBL_INIT(two  , 2.0);
			VEC_DBL_INIT(cc0  , c);
			VEC_DBL_INIT(ss0  , s);
		  #endif
		#else
	//	} else {					// DIT has roots-pair indices 4/8, 6/10, 5/9, 7/11 swapped w.r.to DIF.
									// In AVX2 mode, the sine slots hold tangents, e.g. s4 holds s4/c4:
			r1  = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;
			r2  = sc_ptr + 0x01;		cc0 = sc_ptr + 0x21;
			r3  = sc_ptr + 0x02;		ss0 = sc_ptr + 0x22;
			r4  = sc_ptr + 0x03;		c0  = sc_ptr + 0x23;
			r5  = sc_ptr + 0x04;		s0  = sc_ptr + 0x24;
			r6  = sc_ptr + 0x05;		c4  = sc_ptr + 0x25;
			r7  = sc_ptr + 0x06;		s4  = sc_ptr + 0x26;
			r8  = sc_ptr + 0x07;		c8  = sc_ptr + 0x27;
			r9  = sc_ptr + 0x08;		s8  = sc_ptr + 0x28;
			r10 = sc_ptr + 0x09;		c12 = sc_ptr + 0x29;
			r11 = sc_ptr + 0x0a;		s12 = sc_ptr + 0x2a;
			r12 = sc_ptr + 0x0b;		c2  = sc_ptr + 0x2b;
			r13 = sc_ptr + 0x0c;		s2  = sc_ptr + 0x2c;
			r14 = sc_ptr + 0x0d;		c6  = sc_ptr + 0x2d;
			r15 = sc_ptr + 0x0e;		s6  = sc_ptr + 0x2e;
			r16 = sc_ptr + 0x0f;		c10 = sc_ptr + 0x2f;
			r17 = sc_ptr + 0x10;		s10 = sc_ptr + 0x30;
			r18 = sc_ptr + 0x11;		c14 = sc_ptr + 0x31;
			r19 = sc_ptr + 0x12;		s14 = sc_ptr + 0x32;
			r20 = sc_ptr + 0x13;		c1  = sc_ptr + 0x33;
			r21 = sc_ptr + 0x14;		s1  = sc_ptr + 0x34;
			r22 = sc_ptr + 0x15;		c5  = sc_ptr + 0x35;
			r23 = sc_ptr + 0x16;		s5  = sc_ptr + 0x36;
			r24 = sc_ptr + 0x17;		c9  = sc_ptr + 0x37;
			r25 = sc_ptr + 0x18;		s9  = sc_ptr + 0x38;
			r26 = sc_ptr + 0x19;		c13 = sc_ptr + 0x39;
			r27 = sc_ptr + 0x1a;		s13 = sc_ptr + 0x3a;
			r28 = sc_ptr + 0x1b;		c3  = sc_ptr + 0x3b;
			r29 = sc_ptr + 0x1c;		s3  = sc_ptr + 0x3c;
			r30 = sc_ptr + 0x1d;		c7  = sc_ptr + 0x3d;
			r31 = sc_ptr + 0x1e;		s7  = sc_ptr + 0x3e;
			r32 = sc_ptr + 0x1f;		c11 = sc_ptr + 0x3f;
										s11 = sc_ptr + 0x40;
										c15 = sc_ptr + 0x41;
										s15 = sc_ptr + 0x42;
										two = sc_ptr + 0x43;	// Holds 1.0 in AVX2 mode
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
		  #ifdef USE_AVX2
			VEC_DBL_INIT(two  , 1.0);	// Yeah, I *know" "it's a misnomer" :)
			// cc0,ss0 inited below for AVX2
		  #else
			VEC_DBL_INIT(two  , 2.0);
			VEC_DBL_INIT(cc0  , c);
			VEC_DBL_INIT(ss0  , s);
		  #endif
		#endif
	//	}
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
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
	ASSERT(HERE, p4  == p2+p2, "radix16_dit_pass: p4  != p2+p2!");
	ASSERT(HERE, p8  == p4+p4, "radix16_dit_pass: p8  != p4+p4!");
	ASSERT(HERE, p12 == p4+p8, "radix16_dit_pass: p12 != p4+p8!");

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

	// In AVX2/FMA mode, since we need to replace most of the raw sincos data with derived ones,
	// simply place one copy of each computed double in a double-sized slot of the local memory.
	// We will be using AVX2/FMA-based Newtonian iterative inversion on the 16 doubles whose
	// multiplicative inverse is needed (the real part of the basic root of unity c and of the
	// 15 complex twiddles, c1-15), so store those in packed form in 4 AVX-register-sized
	// contiguous memory locations, and the others in a separate chunk of memory. After the
	// vector-iterative inversion we'll need to combine the 2 sets of data and place (in quadruplicate)
	// into their final SIMD-suitable memory slots.
	#if defined(USE_AVX2) && (defined(DFT_V1) || defined(DFT_V2))

		add0 = (double *)cc0;	// add0 points to 16 cos-data-to-be-inverted; Need a double-ptr on lhs here
		add1 = add0 + 16;		// add1 points to block of memory temporarily used to store the corresponding sine data
		*add0++ = c;	// Since tan0 defined as const, we can init these directly, but init with c0,s0 anyway
		*add1++ = s;	// and use result as a check onthe accuracy of the FMA-based Newton iterative inversion.

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c1, for inversion
		*add1++ = it;	// s1  slot will hold __r1 = s1 /c1 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c2, for inversion
		*add1++ = it;	// s2  slot will hold __r2 = s2 /c2 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c3, for inversion
		*add1++ = it;	// s3  slot will hold __r3 = s3 /c3 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c4, for inversion
		*add1++ = it;	// s4  slot will hold __r4 = s4 /c4 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c5, for inversion
		*add1++ = it;	// s5  slot will hold __r5 = s5 /c5 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c6, for inversion
		*add1++ = it;	// s6  slot will hold __r6 = s6 /c6 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c7, for inversion
		*add1++ = it;	// s7  slot will hold __r7 = s7 /c7 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c8, for inversion
		*add1++ = it;	// s8  slot will hold __r8 = s8 /c8 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c9, for inversion
		*add1++ = it;	// s9  slot will hold __r9 = s9 /c9 

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c10, for inversion
		*add1++ = it;	// s10 slot will hold __rA = s10/c10

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c11, for inversion
		*add1++ = it;	// s11 slot will hold __rB = s11/c11

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c12, for inversion
		*add1++ = it;	// s12 slot will hold __rC = s12/c12

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c13, for inversion
		*add1++ = it;	// s13 slot will hold __rD = s13/c13

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c14, for inversion
		*add1++ = it;	// s14 slot will hold __rE = s14/c14

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		*add0++ = rt;	// c15, for inversion
		*add1++ = it;	// s15 slot will hold __rF = s15/c15

		// This places us at add0 == c8 and add1 = c12.
		ASSERT(HERE, add0 == (double *)cc0+16 && add1 == (double *)cc0+32, "add0,1 checksum failed in AVX2 DIT sincos inits!");
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
		ASSERT(HERE, *(add0-1) == ISRT2, "Scalar ISRT2 bad!");
		c_tmp = cc0 + 0x22;	// 1.0 x 4
		ASSERT(HERE, c_tmp->d0 == 1.0 && c_tmp->d0 == c_tmp->d1 && c_tmp->d0 == c_tmp->d2 && c_tmp->d0 == c_tmp->d3, "1.0 x 4 mismatch!");

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
		(cc0,ss0) + 0x[2,12,a,1a,4,14,c,1c,6,16,e,1e,8,18,10,20]; middle 2 quartet-offsets swapped [4,14,c,1c] <--> [6,16,e,1e] w.r.to DIF:
		*/
		c_tmp = cc0 + 0x02; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x4; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c4 =rt;		s4 =it;
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c5 =rt;		s5 =it;
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0xc; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c6 =rt;		s6 =it;
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x1c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c7 =rt;		s7 =it;
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x6; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c8 =rt;		s8 =it;
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c9 =rt;		s9 =it;
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0xe; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c10=rt;		s10=it;
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	  #ifdef USE_SSE2
		c_tmp = cc0 + 0x1e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	  #else
		c11=rt;		s11=it;
	  #endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c1=t1*t3-t2*t4;	s1=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i;				/* 4*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c2=t1*t3-t2*t4;	s2=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i;				/* 8*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c4=t1*t3-t2*t4;	s4=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += (iroot << 2)+iroot;		/* 13*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c8=t1*t3-t2*t4;	s8=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
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
		j2 = j1 + RE_IM_STRIDE;
		pfetch_addr = p1*((j >> l2_stride) & 0x3);	// cycle prefetch-offset-address among p0,1,2,3
			// These get added to the base addresses p0,4,8,12 in thr DFT macros, thus every 4 loop
			// executions we have covered prefetches from [current address] + [pfetch distance] + p0,1,2,...15 .

	#ifdef USE_SSE2

	#ifdef COMPILER_TYPE_MSVC

	/*...Block 1: */
	#if 1
		add0 = &a[j1];
		__asm	mov eax, add0
		__asm	mov ebx, p1
		__asm	mov ecx, p2
		__asm	mov edx, p3
		__asm	mov	edi, p4		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		SSE2_RADIX4_DIT_0TWIDDLE_B(r1)
	#else
		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r1)
	#endif

	/*...Block 2: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r9)
	#else
		add0 = &a[j1+p4];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r9)
	#endif

	/*...Block 3: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r17)
	#else
		add0 = &a[j1+p8];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r17)
	#endif

	/*...Block 4: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r25)
	#else
		add0 = &a[j1+p12];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r25)
	#endif

	/****************************************************************************************************
	!...and now do four more radix-4 transforms, including the internal and external twiddle factors:   !
	****************************************************************************************************/

	/*...Block 1: Cost: 24 MOVapd, 26 ADD/SUBpd, 12 MULpd */
		__asm	mov eax, add0
		__asm	mov ebx, p4
		__asm	mov ecx, r1
		__asm	mov edx, r9
		__asm	mov edi, p8		/* edi will store copy of p8 throughout */
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	edi, 3
		__asm	add ebx, eax	/* add1 = add0+p4 */

		__asm	movaps	xmm2,[edx      ]	/* t10 */				__asm	movaps	xmm4,[edx+0x100]	/* t30 */
		__asm	movaps	xmm3,[edx+0x010]	/* t11 */				__asm	movaps	xmm5,[edx+0x110]	/* t31 */
		__asm	movaps	xmm0,[ecx      ]	/* t00 */				__asm	movaps	xmm6,[ecx+0x100]	/* t20 */
		__asm	movaps	xmm1,[ecx+0x010]	/* t01 */				__asm	movaps	xmm7,[ecx+0x110]	/* t21 */

		__asm	subpd	xmm0,xmm2			/*~t10=t00-t10*/		__asm	subpd	xmm6,xmm4			/*~t30=t20-t30*/
		__asm	subpd	xmm1,xmm3			/*~t11=t01-t11*/		__asm	subpd	xmm7,xmm5			/*~t31=t21-t31*/
		__asm	addpd	xmm2,xmm2			/*       2*t10*/		__asm	addpd	xmm4,xmm4			/*       2*t30*/
		__asm	addpd	xmm3,xmm3			/*       2*t11*/		__asm	addpd	xmm5,xmm5			/*       2*t31*/
		__asm	addpd	xmm2,xmm0			/*~t00=t00+t10*/		__asm	addpd	xmm4,xmm6			/*~t20=t20+t30*/
		__asm	addpd	xmm3,xmm1			/*~t01=t01+t11*/		__asm	addpd	xmm5,xmm7			/*~t21=t21+t31*/

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0)

	/*...Block 2: Cost: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */
		__asm	mov eax, add0
		__asm	mov esi, p1
		__asm	mov ebx, p4
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p1] */
		__asm	add ebx, eax	/* add1 = add0+p4 */
		__asm	mov ecx, r3
		__asm	mov edx, r11
		__asm	add	ecx, 0x100
		__asm	add	edx, 0x100
		__asm	mov	esi, cc0

		__asm	movaps	xmm4,[ecx      ]	/* t24 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t25 */
		__asm	movaps	xmm2,[esi      ]	/* c */
		__asm	movaps	xmm3,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t24 */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t25 */

		__asm	mulpd	xmm4,xmm2		/* t24*c */
		__asm	mulpd	xmm5,xmm2		/* t25*c */
		__asm	mulpd	xmm6,xmm3		/* t24*s */		__asm	movaps	xmm0,[edx      ]	/* t34 */
		__asm	mulpd	xmm7,xmm3		/* t25*s */		__asm	movaps	xmm1,[edx+0x010]	/* t35 */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t25 */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t34 */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t24 */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t35 */

														__asm	mulpd	xmm0,xmm3		/* t34*s */
														__asm	mulpd	xmm1,xmm3		/* t35*s */
														__asm	mulpd	xmm6,xmm2		/* t34*c */
														__asm	mulpd	xmm7,xmm2		/* t35*c */
														__asm	subpd	xmm1,xmm6	/* xmm5 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm4 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t25*/
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t24*/

		__asm	addpd	xmm4,xmm0	/* ~t24 <- t24+rt */
		__asm	addpd	xmm5,xmm1	/* ~t25 <- t25+it */
		__asm	subpd	xmm6,xmm0	/* ~t34 <- t24-rt */
		__asm	subpd	xmm7,xmm1	/* ~t35 <- t25-it */

		__asm	sub	ecx, 0x100
		__asm	sub	edx, 0x100
		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edx      ]	/* t14 */
		__asm	movaps	xmm3,[edx+0x010]	/* t15 */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm3	/* cpy t15 */
		__asm	subpd	xmm3,xmm2	/*~t15=t15-t14 */
		__asm	addpd	xmm2,xmm0	/*~t14=t14+t15 */
		__asm	mulpd	xmm2,xmm1	/* rt */
		__asm	mulpd	xmm3,xmm1	/* it */

		__asm	movaps	xmm0,[ecx      ]	/* t04 */
		__asm	movaps	xmm1,[ecx+0x010]	/* t05 */
		__asm	subpd	xmm0,xmm2	/*~t14 <- t04- rt */
		__asm	subpd	xmm1,xmm3	/*~t15 <- t05- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t04 <- t04+ rt */
		__asm	addpd	xmm3,xmm1	/*~t05 <- t05+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c1)

	/*...Block 3: Cost: 31 MOVapd, 32 ADD/SUBpd, 20 MULpd */
		__asm	mov eax, add0
		__asm	mov esi, p2
		__asm	mov ebx, p4
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p2] */
		__asm	add ebx, eax	/* add1 = add0+p4 */
		__asm	mov ecx, r5
		__asm	mov edx, r13
		__asm	add	ecx, 0x100
		__asm	add	edx, 0x100
		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[esi]	/* isrt2 */

		__asm	movaps	xmm4,[ecx      ]	/* t28 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t29 */
		__asm	movaps	xmm6,[edx      ]	/* t38 */
		__asm	movaps	xmm7,[edx+0x010]	/* t39 */
		__asm	sub	ecx, 0x100
		__asm	sub	edx, 0x100
		__asm	mulpd	xmm4,xmm2
		__asm	mulpd	xmm5,xmm2
		__asm	mulpd	xmm6,xmm2
		__asm	mulpd	xmm7,xmm2

		__asm	subpd	xmm5,xmm4			/*~t29=t29-t28*/		__asm	movaps	xmm0,[ecx      ]	/* t08 */
		__asm	subpd	xmm6,xmm7			/* rt =t38-t39*/		__asm	movaps	xmm2,[edx+0x010]	/* t19 */
		__asm	addpd	xmm4,xmm4			/*       2*t28*/		__asm	movaps	xmm3,[ecx+0x010]	/* t09 */
		__asm	addpd	xmm7,xmm7			/*       2*t39*/		__asm	movaps	xmm1,[edx      ]	/* t18 */
		__asm	addpd	xmm4,xmm5			/*~t28=t28+t29*/
		__asm	addpd	xmm7,xmm6			/* it =t39+t38*/

		__asm	subpd	xmm4,xmm6			/*~t28=t28-rt */		__asm	subpd	xmm0,xmm2			/*~t18=t08-t19*/
		__asm	subpd	xmm5,xmm7			/*~t29=t29-it */		__asm	subpd	xmm3,xmm1			/*~t09=t09-t18*/
		__asm	addpd	xmm6,xmm6			/*       2*rt */		__asm	addpd	xmm2,xmm2			/*       2*t08*/
		__asm	addpd	xmm7,xmm7			/*       2*it */		__asm	addpd	xmm1,xmm1			/*       2*t09*/
		__asm	addpd	xmm6,xmm4			/*~t38=t28+rt */		__asm	addpd	xmm2,xmm0			/*~t08=t19+t08*/
		__asm	addpd	xmm7,xmm5			/*~t39=t29+it */		__asm	addpd	xmm1,xmm3			/*~t19=t18+t09*/

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c2)

	/*...Block 4: Cost: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */
	#if 0
		add0 = &a[j1+p3];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;

		__asm	mov	eax, r7
		__asm	mov	edx, eax
		__asm	add	eax, 0x100
		__asm	add	edx, 0x180
		__asm	mov	ecx, cc0

		__asm	movaps	xmm4,[eax      ]	/* t2C */
		__asm	movaps	xmm5,[eax+0x010]	/* t2D */
		__asm	movaps	xmm3,[ecx      ]	/* c */
		__asm	movaps	xmm2,[ecx+0x010]	/* s */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t2C */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t2D */

		__asm	mulpd	xmm4,xmm2		/* t2C*s */
		__asm	mulpd	xmm5,xmm2		/* t2D*s */
		__asm	mulpd	xmm6,xmm3		/* t2C*c */		__asm	movaps	xmm0,[edx      ]	/* t3C */
		__asm	mulpd	xmm7,xmm3		/* t2D*c */		__asm	movaps	xmm1,[edx+0x010]	/* t3D */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t2D */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t3C */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t2C */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t3D */

														__asm	mulpd	xmm0,xmm3		/* t3C*c */
														__asm	mulpd	xmm1,xmm3		/* t3D*c */
														__asm	mulpd	xmm6,xmm2		/* t3C*s */
														__asm	mulpd	xmm7,xmm2		/* t3D*s */
														__asm	subpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t2D*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t2C*/

		__asm	addpd	xmm6,xmm0	/* ~t3C <- t2C+rt */
		__asm	addpd	xmm7,xmm1	/* ~t3D <- t2D+it */
		__asm	subpd	xmm4,xmm0	/* ~t2C <- t2C-rt */
		__asm	subpd	xmm5,xmm1	/* ~t2D <- t2D-it */

		__asm	sub	eax, 0x100
		__asm	sub	edx, 0x100
		__asm	mov	ecx, isrt2
		__asm	movaps	xmm0,[edx      ]	/* t1C */
		__asm	movaps	xmm1,[edx+0x010]	/* t1D */
		__asm	movaps	xmm3,[ecx]	/* isrt2 */
		__asm	movaps	xmm2,xmm0	/* cpy t1C */
		__asm	subpd	xmm0,xmm1	/*~t1C=t1C-t1D */
		__asm	addpd	xmm1,xmm2	/*~t1D=t1D+t1C */
		__asm	mulpd	xmm0,xmm3	/* it */
		__asm	mulpd	xmm1,xmm3	/* rt */

		__asm	movaps	xmm2,[eax      ]	/* t0C */
		__asm	movaps	xmm3,[eax+0x010]	/* t0D */
		__asm	subpd	xmm2,xmm0	/*~t0C <- t0C- rt */
		__asm	subpd	xmm3,xmm1	/*~t0D <- t0D- it */
		__asm	addpd	xmm0,xmm0	/*          2* rt */
		__asm	addpd	xmm1,xmm1	/*          2* it */
		__asm	addpd	xmm0,xmm2	/*~t1C <- t0C+ rt */
		__asm	addpd	xmm1,xmm3	/*~t1D <- t0D+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF(add0, add1, add2, add3, c3)
	#else
		__asm	mov eax, add0
		__asm	mov esi, p3
		__asm	mov ebx, p4
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p2] */
		__asm	add ebx, eax	/* add1 = add0+p4 */
		__asm	mov ecx, r7
		__asm	mov edx, r15
		__asm	add	ecx, 0x100
		__asm	add	edx, 0x100
		__asm	mov	esi, cc0

		__asm	movaps	xmm4,[ecx      ]	/* t2C */
		__asm	movaps	xmm5,[ecx+0x010]	/* t2D */
		__asm	movaps	xmm3,[esi      ]	/* c */
		__asm	movaps	xmm2,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t2C */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t2D */

		__asm	mulpd	xmm4,xmm2		/* t2C*s */
		__asm	mulpd	xmm5,xmm2		/* t2D*s */
		__asm	mulpd	xmm6,xmm3		/* t2C*c */		__asm	movaps	xmm0,[edx      ]	/* t3C */
		__asm	mulpd	xmm7,xmm3		/* t2D*c */		__asm	movaps	xmm1,[edx+0x010]	/* t3D */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t2D */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t3C */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t2C */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t3D */

														__asm	mulpd	xmm0,xmm3		/* t3C*c */
														__asm	mulpd	xmm1,xmm3		/* t3D*c */
														__asm	mulpd	xmm6,xmm2		/* t3C*s */
														__asm	mulpd	xmm7,xmm2		/* t3D*s */
														__asm	subpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t2D*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t2C*/

		__asm	addpd	xmm6,xmm0	/* ~t3C <- t2C+rt */
		__asm	addpd	xmm7,xmm1	/* ~t3D <- t2D+it */
		__asm	subpd	xmm4,xmm0	/* ~t2C <- t2C-rt */
		__asm	subpd	xmm5,xmm1	/* ~t2D <- t2D-it */

		__asm	sub	ecx, 0x100
		__asm	sub	edx, 0x100
		__asm	mov	esi, isrt2
		__asm	movaps	xmm0,[edx      ]	/* t1C */
		__asm	movaps	xmm1,[edx+0x010]	/* t1D */
		__asm	movaps	xmm3,[esi]	/* isrt2 */
		__asm	movaps	xmm2,xmm0	/* cpy t1C */
		__asm	subpd	xmm0,xmm1	/*~t1C=t1C-t1D */
		__asm	addpd	xmm1,xmm2	/*~t1D=t1D+t1C */
		__asm	mulpd	xmm0,xmm3	/* it */
		__asm	mulpd	xmm1,xmm3	/* rt */

		__asm	movaps	xmm2,[ecx      ]	/* t0C */
		__asm	movaps	xmm3,[ecx+0x010]	/* t0D */
		__asm	subpd	xmm2,xmm0	/*~t0C <- t0C- rt */
		__asm	subpd	xmm3,xmm1	/*~t0D <- t0D- it */
		__asm	addpd	xmm0,xmm0	/*          2* rt */
		__asm	addpd	xmm1,xmm1	/*          2* it */
		__asm	addpd	xmm0,xmm2	/*~t1C <- t0C+ rt */
		__asm	addpd	xmm1,xmm3	/*~t1D <- t0D+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c3)
	#endif

		/***************************************************/
		/* Total Cost: 187 MOVapd, 210 ADD/SUBpd, 84 MULpd */
		/***************************************************/

	#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		add0 = &a[j1];

	  #ifdef USE_AVX2

	   #ifdef DFT_V1	// Toggle between 2 versions

		SSE2_RADIX16_DIT_TWIDDLE_1(add0,p1,p2,p3,p4,p8,p12,r1,cc0,pfetch_addr,pfetch_dist);

	   #elif defined(DFT_V2)

		SSE2_RADIX16_DIT_TWIDDLE_2(add0,p1,p2,p3,p4,p8,p12,r1,cc0,pfetch_addr,pfetch_dist);

	   #else	// Allow non-FMA as the default option for comparative timing purposes

		SSE2_RADIX16_DIT_TWIDDLE(add0,p1,p2,p3,p4,p8,r1,isrt2,pfetch_addr,pfetch_dist);

	   #endif

	  // SSE2 code has different macro arglists for 32 and 64-bit modes:
	  #elif OS_BITS == 64
		SSE2_RADIX16_DIT_TWIDDLE(add0,p1,p2,p3,p4,p8,r1,isrt2);
	  #else
		SSE2_RADIX16_DIT_TWIDDLE(add0,p1,p2,p3,p4,p8,r1,r3,r5,r7,r9,r11,r13,r15,r17,r25,isrt2,cc0, c0, c1, c2, c3);
	  #endif

	#endif

#else	/* USE_SSE2 */

  #ifdef USE_SCALAR_DFT_MACRO	// Must define - or not - @compile time

		// Test FMA-based DIT macro:
		RADIX_16_DIT_FMA(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p4+p1 ],a[j2+p4+p1 ],a[j1+p4+p2 ],a[j2+p4+p2 ],a[j1+p4+p3 ],a[j2+p4+p3 ],a[j1+p8 ],a[j2+p8 ],a[j1+p8+p1 ],a[j2+p8+p1 ],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+p12],a[j2+p12],a[j1+p12+p1],a[j2+p12+p1],a[j1+p12+p2],a[j2+p12+p2],a[j1+p12+p3],a[j2+p12+p3]
						,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p4+p1 ],a[j2+p4+p1 ],a[j1+p4+p2 ],a[j2+p4+p2 ],a[j1+p4+p3 ],a[j2+p4+p3 ],a[j1+p8 ],a[j2+p8 ],a[j1+p8+p1 ],a[j2+p8+p1 ],a[j1+p8+p2],a[j2+p8+p2],a[j1+p8+p3],a[j2+p8+p3],a[j1+p12],a[j2+p12],a[j1+p12+p1],a[j2+p12+p1],a[j1+p12+p2],a[j2+p12+p2],a[j1+p12+p3],a[j2+p12+p3]
						,c1,s1,c2,s2,c3,s3,c4,s4,c5,s5,c6,s6,c7,s7,c8,s8,c9,s9,c10,s10,c11,s11,c12,s12,c13,s13,c14,s14,c15,s15, c,s)

  #else		// USE_SCALAR_DFT_MACRO = False

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/
	#ifdef PFETCH_AGGRESSIVE
		addp = &a[j1];
	#elif PFETCH
		prefetch_offset = ((j >> 1) & 3)*p1 + 4;	/* Cycle among p0, p1, p2 and p3. */
	#endif
	/*...Block 1: */

		jt = j1;
		jp = j2;

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

		jt = j1 + p4;
		jp = j2 + p4;

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

		jt = j1 + p8;
		jp = j2 + p8;

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

		jt = j1 + p12;
		jp = j2 + p12;

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

		jt = j1;
		jp = j2;

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

		jt = j1 + p2;
		jp = j2 + p2;

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

		jt = j1 + p1;
		jp = j2 + p1;

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

		jt = j1 + p3;
		jp = j2 + p3;

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

  #endif	// USE_SCALAR_DFT_MACRO ?

#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

