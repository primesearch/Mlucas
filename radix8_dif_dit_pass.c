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

#ifdef USE_SSE2

	#if(HIACC != 1)
		#error SIMD Mode requires HIACC flag to be set!
	#endif

	#ifdef COMPILER_TYPE_MSVC
		#include "sse2_macro.h"
	#endif

	#if defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix8_dif_dit_pass_gcc32.h"

		#else

			#include "radix8_dif_dit_pass_gcc64.h"

		#endif

	#endif

#endif

/***************/

/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform a single radix-8 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details.
*/
void radix8_dif_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
{
	static int max_threads = 0;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int i,j,j1,j2,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	int p1,p2,p3,p4,p5,p6,p7;
	double rt,it;
	double re0,im0,re1,im1;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	/* This allows us to align sincos temporaries on 16-byte boundaries (e.g. for SSE2) via pointers:
	In sse2 mode, alloc a complex[32] array (give it a few extra slots to allow ptr to initial element-to-be-used
	to be aligned on 16-byte bdry, and 2 more for the doubled 1/sqrt2 pair),
	and then store doubled pairs of sincos data - for each scalar sincos datum used
	in non-SSE2 mode, we have a double-pair consisting of 2 copies of the same datum, in SSE2 mode.
	This relies on the fact that in our FFT algorithm, the same set of sincos multipliers computed in the outer
	loop below is re-used for multiple (and an even number, in particular) inputs in the inner loop, i.e.
	in SSE2 mode we simply double the inner-loop stride, and process vector DFT inputs in pairs.
	*/
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */
	vec_dbl *c_tmp,*s_tmp;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
	vec_dbl *isrt2,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;	// In || mode, only above base-pointer (shared by all threads) is static:
  #elif defined(COMPILER_TYPE_GCC)
	static vec_dbl *isrt2,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
  #else
	static vec_dbl *r0,*r8;
	static vec_dbl *isrt2,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
  #endif

#else

	double *addr, *addp;	/* addresses for prefetch macros */
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16;
	double c1,c2,c3,c4,c5,c6,c7,s1,s2,s3,s4,s5,s6,s7;

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
		sc_arr = ALLOC_VEC_DBL(sc_arr, 36*max_threads);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 16 16-byte slots of sc_arr for temporaries, next 16 for the doubled sincos twiddles,
	next 1 for doubled 1/sqrt2, plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
		#ifdef MULTITHREAD
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x20;	// Need 16 local vec_dbl-sized slots for tmp-data, 16 for sincos twiddles, 1 for 1/sqrt2
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);
				isrt2 += 36;	/* Move on to next thread's local store */
			}
		#elif defined(COMPILER_TYPE_GCC)
			c0    = sc_ptr + 0x10;
			c4    = sc_ptr + 0x12;
			c2    = sc_ptr + 0x14;
			c6    = sc_ptr + 0x16;
			c1    = sc_ptr + 0x18;
			c5    = sc_ptr + 0x1a;
			c3    = sc_ptr + 0x1c;
			c7    = sc_ptr + 0x1e;
			isrt2 = sc_ptr + 0x20;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
		#else
			r0    = sc_ptr + 0x00;		c0    = sc_ptr + 0x10;
		//	r1    = sc_ptr + 0x01;		s0    = sc_ptr + 0x11;
		/*	r2    = sc_ptr + 0x02;	*/	c4    = sc_ptr + 0x12;
		//	r3    = sc_ptr + 0x03;		s4    = sc_ptr + 0x13;
		/*	r4    = sc_ptr + 0x04;	*/	c2    = sc_ptr + 0x14;
		//	r5    = sc_ptr + 0x05;		s2    = sc_ptr + 0x15;
		/*	r6    = sc_ptr + 0x06;	*/	c6    = sc_ptr + 0x16;
		//	r7    = sc_ptr + 0x07;		s6    = sc_ptr + 0x17;
			r8    = sc_ptr + 0x08;		c1    = sc_ptr + 0x18;
		//	r9    = sc_ptr + 0x09;		s1    = sc_ptr + 0x19;
		/*	ra    = sc_ptr + 0x0a;	*/	c5    = sc_ptr + 0x1a;
		//	rb    = sc_ptr + 0x0b;		s5    = sc_ptr + 0x1b;
		/*	rc    = sc_ptr + 0x0c;	*/	c3    = sc_ptr + 0x1c;
		//	rd    = sc_ptr + 0x0d;		s3    = sc_ptr + 0x1d;
		/*	re    = sc_ptr + 0x0e;	*/	c7    = sc_ptr + 0x1e;
		//	rf    = sc_ptr + 0x0f;		s7    = sc_ptr + 0x1f;
										isrt2 = sc_ptr + 0x20;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);
		#endif
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		c0 = __r0 + thr_id*36 + 0x10;
		c4    = c0 + 0x02;
		c2    = c0 + 0x04;
		c6    = c0 + 0x06;
		c1    = c0 + 0x08;
		c5    = c0 + 0x0a;
		c3    = c0 + 0x0c;
		c7    = c0 + 0x0e;
		isrt2 = c0 + 0x10;
	#endif

#endif

/*
!...The radix-8 pass is here. The data are processed in NLOOPS blocks of INCR real elements each,
!   where mloops = product of all preceding radices. The fundamental root of unity for this pass
!   is thus the primitive [(product of all preceding radices)*radix_now] =  [mloops*radix_now]th root,
!   and we get the proper table look-up index by multiplying the BR index by N/[mloops*radix_now] = (incr/2)/radix_now.
*/

	/*  stride between elements processed together. */
	p1 = incr>>3;
	p2 = p1 + p1;
	p3 = p2 + p1;
	p4 = p3 + p1;
	p5 = p4 + p1;
	p6 = p5 + p1;
	p7 = p6 + p1;

	p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
	p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );

	iroot_prim=(incr >> 4);	/* (incr/2)/radix_now */

/*
printf("radix8_dif_pass: nloops = %d  incr = %d  using %d th roots of unity.\n",nloops,incr,n/(2*iroot_prim));
*/
	for(m=0; m < nloops; m++)
	{
/*	here are the needed sincos data - these are processed below in bit-reversed order. */
	  iroot = index[m]*iroot_prim;
	  i = iroot;
/* We may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks. */
#if HIACC
	#ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-7] are offset w.r.to the thread-local ptr pair as
		(c0,s0) + 0x[0,8,4,c,2,a,6,e]:
		*/
		c_tmp = c0 + 0x0; s_tmp = c_tmp+1;	/* c0,s0 */
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
		c_tmp = c0 + 0x8; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0x4; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0xc; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0x2; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0xa; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0x6; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0xe; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c7 =rt;		s7 =it;
	#endif

#else

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c1=t1*t3-t2*t4;	s1=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i + iroot;
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c2=t1*t3-t2*t4;	s2=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;

		c5=t1*t3-t2*t4;	s5=t1*t4+t2*t3;

		t1=c1*c5; t2=c1*s5; t3=s1*c5; t4=s1*s5;
		c4=t1+t4; s4=t2-t3; c6=t1-t4; s6=t2+t3;

		t1=c2*c5; t2=c2*s5; t3=s2*c5; t4=s2*s5;
		c3=t1+t4; s3=t2-t3; c7=t1-t4; s7=t2+t3;
#endif

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 3);

	  for(j = jlo; j < jhi; j += stride)
	  {
		j1 = j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	#ifdef USE_SSE2

	  #if defined(COMPILER_TYPE_MSVC)	/* MSVC-style inline ASM: */

		/* Only need to store addresses of real element of each pair - im offset by 0x10
		(i.e. 2 doubles) from that, due to the [re0,re1,im0,im1] ordering we use in SSE2 mode:
		*/
		add0 = &a[j1];
	   #if 1
		SSE2_RADIX8_DIF_TWIDDLE(add0,p1,p2,p4,p6,r0,c0)
												/* 117 load/store [91 movaps], 82 add/subpd, 32 mulpd, 67 address-compute */

	   #elif 1
 		__asm	mov	eax, add0
		__asm	mov	ebx, p2		// Can't get these via simple load-one-and-shift-as-needed due to array padding scheme
		__asm	mov	ecx, p4
		__asm	mov	edx, p6
		__asm	shl	ebx,  3
		__asm	shl	ecx,  3
		__asm	shl	edx,  3
		__asm	add	ebx, eax
		__asm	add	ecx, eax
		__asm	add	edx, eax
		SSE2_RADIX4_DIF_4TWIDDLE_A         (r0,c4)	/* 36 load/store [30 movaps], 26 add/subpd, 14 mulpd, 21 address-compute */
		__asm	mov	edi, p1
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p1];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1)	/* 41 load/store [35 movaps], 32 add/subpd, 20 mulpd, 23 address-compute */
		/* Combine the 2 radix-4 subtransforms and write outputs back to main array: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0)		/* 32 load/store [32 movaps], 24 add/subpd,  0 mulpd, 23 address-compute */

												/* Totals:109 load/store [97 movaps], 82 add/subpd, 32 mulpd, 67 address-compute */
												/* Versus 140 load/store [100movaps], 82 add/subpd, 32 mulpd, 72 address-compute for the original implementation below. */
	   #else
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		add4 = add0+p4;
		add5 = add0+p5;
		add6 = add0+p6;
		add7 = add0+p7;

		/*
		t1 =a[j1   ];					t2 =a[j2   ];
		rt =a[j1+p4]*c4 -a[j2+p4]*s4;	it =a[j2+p4]*c4 +a[j1+p4]*s4;
		t3 =t1 -rt;				t4 =t2 -it;
		t1 =t1 +rt;				t2 =t2 +it;
		*/
		__asm	mov	eax, add0
		__asm	mov	ebx, add4
		__asm	mov	ecx,   c4	/* &c4 */

		__asm	movaps	xmm0,[eax     ]	/* t1 */
		__asm	movaps	xmm1,[eax+0x10]	/* t2 */
		__asm	movaps	xmm6,[eax     ]	/* cpy t1 */
		__asm	movaps	xmm7,[eax+0x10]	/* cpy t2 */

		__asm	movaps	xmm2,[ebx     ]	/* a[j1+p4] */
		__asm	movaps	xmm3,[ebx+0x10]	/* a[j2+p4] */
		__asm	movaps	xmm4,[ebx     ]	/* xmm4 <- cpy a[j1+p4] */
		__asm	movaps	xmm5,[ebx+0x10]	/* xmm5 <- cpy a[j2+p4] */
		__asm	mulpd	xmm2,[ecx     ]	/* a[j1+p4]*c4 */
		__asm	mulpd	xmm3,[ecx     ]	/* a[j2+p4]*c4 */
		__asm	mulpd	xmm4,[ecx+0x10]	/* a[j1+p4]*s4 */
		__asm	mulpd	xmm5,[ecx+0x10]	/* a[j2+p4]*s4 */
		__asm	subpd	xmm2,xmm5	/* xmm2 <- rt */
		__asm	addpd	xmm3,xmm4	/* xmm3 <- it */

		__asm	addpd	xmm0,xmm2	/* t1<-t1 +rt */
		__asm	addpd	xmm1,xmm3	/* t2<-t2 +it */
		__asm	subpd	xmm6,xmm2	/* t3 =t1 -rt */
		__asm	subpd	xmm7,xmm3	/* t4 =t2 -it */

		/* Store results back into original array slots: */
		__asm	movaps	[eax     ],xmm0	/* a[j1   ] <- t1 */
		__asm	movaps	[eax+0x10],xmm1	/* a[j2   ] <- t2 */
		__asm	movaps	[ebx     ],xmm6	/* a[j1+p4] <- t3 */
		__asm	movaps	[ebx+0x10],xmm7 /* a[j2+p4] <- t4 */

		/*
		t5 =a[j1+p2]*c2 -a[j2+p2]*s2;	t6 =a[j2+p2]*c2 +a[j1+p2]*s2;
		rt =a[j1+p6]*c6 -a[j2+p6]*s6;	it =a[j2+p6]*c6 +a[j1+p6]*s6;
		t7 =t5 -rt;				t8 =t6 -it;
		t5 =t5 +rt;				t6 =t6 +it;
		*/
		__asm	mov	eax, add2
		__asm	mov	ebx, add6
		__asm	mov	ecx,   c2
		__asm	mov	edx,   c6
		__asm	movaps	xmm0,[eax     ]	/* xmm0 <- a[j1+p2] */
		__asm	movaps	xmm2,[eax+0x10]	/* xmm2 <- a[j2+p2] */
		__asm	movaps	xmm1,[eax     ]	/* xmm1 <- cpy a[j1+p2] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[j2+p2] */
		__asm	mulpd	xmm0,[ecx     ]	/*~xmm0 <- a[j1+p2]*c2 */
		__asm	mulpd	xmm2,[ecx+0x10]	/*~xmm2 <- a[j2+p2]*s2 */
		__asm	mulpd	xmm1,[ecx+0x10]	/*~xmm1 <- a[j1+p2]*s2 */
		__asm	mulpd	xmm3,[ecx     ]	/*~xmm3 <- a[j2+p2]*c2 */
		__asm	subpd	xmm0,xmm2	/*     xmm0 <- t5  */
		__asm	addpd	xmm1,xmm3	/*     xmm1 <- t6  */

		__asm	movaps	xmm2,[ebx     ]	/* xmm2 <- a[j1+p6] */
		__asm	movaps	xmm3,[ebx+0x10]	/* xmm3 <- a[j2+p6] */
		__asm	movaps	xmm4,[ebx     ]	/* xmm4 <- cpy a[j1+p6] */
		__asm	movaps	xmm5,[ebx+0x10]	/* xmm5 <- cpy a[j2+p6] */
		__asm	mulpd	xmm2,[edx     ]	/*~xmm2 <- a[j1+p6]*c6 */
		__asm	mulpd	xmm3,[edx+0x10]	/*~xmm3 <- a[j2+p6]*s6 */
		__asm	mulpd	xmm4,[edx+0x10]	/*~xmm4 <- a[j1+p6]*s6 */
		__asm	mulpd	xmm5,[edx     ]	/*~xmm5 <- a[j2+p6]*c6 */
		__asm	subpd	xmm2,xmm3	/*     xmm2 <- rt */
		__asm	addpd	xmm4,xmm5	/*     xmm4 <- it */

		__asm	movaps	xmm3,xmm2	/*     xmm3 <- cpy rt */
		__asm	movaps	xmm5,xmm4	/*     xmm5 <- cpy it */
		__asm	addpd	xmm2,xmm0	/*    ~xmm2 <- t5 =t5 +rt */
		__asm	subpd	xmm0,xmm3	/*    ~xmm0 <- t7 =t5 -rt */
		__asm	addpd	xmm4,xmm1	/*    ~xmm4 <- t6 =t6 +it */
		__asm	subpd	xmm1,xmm5	/*    ~xmm1 <- t8 =t6 -it */

		/* Store results back into original array slots: */
		__asm	movaps	[eax     ],xmm2	/* a[j1+p2] <- t5 */
		__asm	movaps	[eax+0x10],xmm4	/* a[j2+p2] <- t6 */
		__asm	movaps	[ebx     ],xmm0	/* a[j1+p6] <- t7 */
		__asm	movaps	[ebx+0x10],xmm1	/* a[j2+p6] <- t8 */

		/*
		t9 =a[j1+p1]*c1 -a[j2+p1]*s1;	t10=a[j2+p1]*c1 +a[j1+p1]*s1;
		rt =a[j1+p5]*c5 -a[j2+p5]*s5;	it =a[j2+p5]*c5 +a[j1+p5]*s5;
		t11=t9 -rt;				t12=t10-it;
		t9 =t9 +rt;				t10=t10+it;
		*/
		__asm	mov	eax, add1
		__asm	mov	ebx, add5
		__asm	mov	ecx,   c1
		__asm	mov	edx,   c5
		__asm	movaps	xmm0,[eax     ]	/* xmm0 <- a[j1+p1] */
		__asm	movaps	xmm2,[eax+0x10]	/* xmm2 <- a[j2+p1] */
		__asm	movaps	xmm1,[eax     ]	/* xmm1 <- cpy a[j1+p1] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[j2+p1] */
		__asm	mulpd	xmm0,[ecx     ]	/*~xmm0 <- a[j1+p1]*c1 */
		__asm	mulpd	xmm2,[ecx+0x10]	/*~xmm2 <- a[j2+p1]*s1 */
		__asm	mulpd	xmm1,[ecx+0x10]	/*~xmm1 <- a[j1+p1]*s1 */
		__asm	mulpd	xmm3,[ecx     ]	/*~xmm3 <- a[j2+p1]*c1 */
		__asm	subpd	xmm0,xmm2	/*     xmm0 <- t5  */
		__asm	addpd	xmm1,xmm3	/*     xmm1 <- t6  */

		__asm	movaps	xmm2,[ebx     ]	/* xmm2 <- a[j1+p5] */
		__asm	movaps	xmm3,[ebx+0x10]	/* xmm3 <- a[j2+p5] */
		__asm	movaps	xmm4,[ebx     ]	/* xmm4 <- cpy a[j1+p5] */
		__asm	movaps	xmm5,[ebx+0x10]	/* xmm5 <- cpy a[j2+p5] */
		__asm	mulpd	xmm2,[edx     ]	/*    ~xmm2 <- a[j1+p5]*c5 */
		__asm	mulpd	xmm3,[edx+0x10]	/*    ~xmm3 <- a[j2+p5]*s5 */
		__asm	mulpd	xmm4,[edx+0x10]	/*    ~xmm4 <- a[j1+p5]*s5 */
		__asm	mulpd	xmm5,[edx     ]	/*    ~xmm5 <- a[j2+p5]*c5 */
		__asm	subpd	xmm2,xmm3	/*     xmm2 <- rt */
		__asm	addpd	xmm4,xmm5	/*     xmm4 <- it */

		__asm	movaps	xmm3,xmm2	/*     xmm3 <- cpy rt */
		__asm	movaps	xmm5,xmm4	/*     xmm5 <- cpy it */
		__asm	addpd	xmm2,xmm0	/*    ~xmm2 <- t9 =t9 +rt */
		__asm	subpd	xmm0,xmm3	/*    ~xmm0 <- t11=t9 -rt */
		__asm	addpd	xmm4,xmm1	/*    ~xmm4 <- t10=t10+it */
		__asm	subpd	xmm1,xmm5	/*    ~xmm1 <- t12=t10-it */

		/* Store results back into original array slots: */
		__asm	movaps	[eax     ],xmm2	/* a[j1+p2] <- t9  */
		__asm	movaps	[eax+0x10],xmm4	/* a[j2+p2] <- t10 */
		__asm	movaps	[ebx     ],xmm0	/* a[j1+p6] <- t11 */
		__asm	movaps	[ebx+0x10],xmm1	/* a[j2+p6] <- t12 */

		/*
		t13=a[j1+p3]*c3 -a[j2+p3]*s3;	t14=a[j2+p3]*c3 +a[j1+p3]*s3;
		rt =a[j1+p7]*c7 -a[j2+p7]*s7;	it =a[j2+p7]*c7 +a[j1+p7]*s7;
		t15=t13-rt;				t16=t14-it;
		t13=t13+rt;				t14=t14+it;
		*/
		__asm	mov	eax, add3
		__asm	mov	ebx, add7
		__asm	mov	ecx,   c3
		__asm	mov	edx,   c7
		__asm	movaps	xmm0,[eax     ]	/* xmm0 <- a[j1+p2] */
		__asm	movaps	xmm2,[eax+0x10]	/* xmm2 <- a[j2+p2] */
		__asm	movaps	xmm1,[eax     ]	/* xmm1 <- cpy a[j1+p2] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[j2+p2] */
		__asm	mulpd	xmm0,[ecx     ]	/*    ~xmm0 <- a[j1+p2]*c2 */
		__asm	mulpd	xmm2,[ecx+0x10]	/*    ~xmm2 <- a[j2+p2]*s2 */
		__asm	mulpd	xmm1,[ecx+0x10]	/*    ~xmm1 <- a[j1+p2]*s2 */
		__asm	mulpd	xmm3,[ecx     ]	/*    ~xmm3 <- a[j2+p2]*c2 */
		__asm	subpd	xmm0,xmm2	/*     xmm0 <- t5  */
		__asm	addpd	xmm1,xmm3	/*     xmm1 <- t6  */

		__asm	movaps	xmm2,[ebx     ]	/* xmm2 <- a[j1+p6] */
		__asm	movaps	xmm3,[ebx+0x10]	/* xmm3 <- a[j2+p6] */
		__asm	movaps	xmm4,[ebx     ]	/* xmm4 <- cpy a[j1+p6] */
		__asm	movaps	xmm5,[ebx+0x10]	/* xmm5 <- cpy a[j2+p6] */
		__asm	mulpd	xmm2,[edx     ]	/*    ~xmm2 <- a[j1+p6]*c6 */
		__asm	mulpd	xmm3,[edx+0x10]	/*    ~xmm3 <- a[j2+p6]*s6 */
		__asm	mulpd	xmm4,[edx+0x10]	/*    ~xmm4 <- a[j1+p6]*s6 */
		__asm	mulpd	xmm5,[edx     ]	/*    ~xmm5 <- a[j2+p6]*c6 */
		__asm	subpd	xmm2,xmm3	/*     xmm2 <- rt */
		__asm	addpd	xmm4,xmm5	/*     xmm4 <- it */

		__asm	movaps	xmm3,xmm2	/*     xmm3 <- cpy rt */
		__asm	movaps	xmm5,xmm4	/*     xmm5 <- cpy it */
		__asm	addpd	xmm2,xmm0	/*    ~xmm2 <- t13=t13+rt */
		__asm	subpd	xmm0,xmm3	/*    ~xmm0 <- t15=t13-rt */
		__asm	addpd	xmm4,xmm1	/*    ~xmm4 <- t14=t14+it */
		__asm	subpd	xmm1,xmm5	/*    ~xmm1 <- t16=t14-it */

		/* Store results back into original array slots: */
		__asm	movaps	[eax     ],xmm2	/* a[j1+p2] <- t13 */
		__asm	movaps	[eax+0x10],xmm4	/* a[j2+p2] <- t14 */
		__asm	movaps	[ebx     ],xmm0	/* a[j1+p6] <- t15 */
		__asm	movaps	[ebx+0x10],xmm1	/* a[j2+p6] <- t16 */

		/*************************************************/
		/* combine to get #1 of 2 length-4 transforms... */
		/*************************************************/

		/*
		rt =t5;					it =t6;
		t5 =t1 -rt;				t6 =t2 -it;
		t1 =t1 +rt;				t2 =t2 +it;

		rt =t13;				it =t14;
		t13=t9 -rt;				t14=t10-it;
		t9 =t9 +rt;				t10=t10+it;
		*/
		__asm	mov	eax, add0			/*&t1 */
		__asm	mov	ebx, add2			/*&t5 */
		__asm	movaps	xmm0,[eax]		/* t1 */
		__asm	movaps	xmm1,[eax+0x10]	/* t2 */
		__asm	movaps	xmm4,xmm0	/* copy of t1 */
		__asm	movaps	xmm5,xmm1	/* copy of t2 */
		/* */
		__asm	addpd	xmm0,[ebx]		/* t1 =t1 +rt */
		__asm	subpd	xmm4,[ebx]		/* t5 =t1 -rt */
		/* */
		__asm	addpd	xmm1,[ebx+0x10]	/* t2 =t2 +it */
		__asm	subpd	xmm5,[ebx+0x10]	/* t6 =t2 -it */
		/* Store copies of t1,2,5,6 into destination array slots: */
		__asm	movaps	[eax     ],xmm0	/* a[j1   ]=t1 */
		__asm	movaps	[eax+0x10],xmm1	/* a[j2   ]=t2 */
		__asm	movaps	[ebx     ],xmm4	/* a[j1+p2]=t5 */
		__asm	movaps	[ebx+0x10],xmm5	/* a[j2+p2]=t6 */

		__asm	mov	ecx, add1			/*&t9 */
		__asm	mov	edx, add3			/*&t13*/
		__asm	movaps	xmm2,[ecx     ]	/* t9 */
		__asm	movaps	xmm3,[ecx+0x10]	/* t10*/
		__asm	movaps	xmm6,xmm2	/* copy of t9 */
		__asm	movaps	xmm7,xmm3	/* copy of t10*/
		/* */
		__asm	addpd	xmm2,[edx]		/* t9 =t9 +rt */
		__asm	subpd	xmm6,[edx]		/* t13=t9 -rt */
		/* */
		__asm	addpd	xmm3,[edx+0x10]	/* t10=t10+it */
		__asm	subpd	xmm7,[edx+0x10]	/* t14=t10-it */

		/* now combine the two half-transforms */
		/*
		a[j1   ]=t1+t9;			a[j2   ]=t2+t10;
		a[j1+p1]=t1-t9;			a[j2+p1]=t2-t10;

		a[j1+p2]=t5-t14;		a[j2+p2]=t6+t13;
		a[j1+p3]=t5+t14;		a[j2+p3]=t6-t13;
		*/
		__asm	subpd	xmm0,xmm2	/* t1 <- t1-t9  */
		__asm	subpd	xmm1,xmm3	/* t2 <- t2-t10 */
		__asm	subpd	xmm4,xmm7	/* t5 <- t5-t14 */
		__asm	subpd	xmm5,xmm6	/* t6 <- t6-t13 */

		__asm	addpd	xmm2,[eax     ]	/* t9 <- t9 +t1 */
		__asm	addpd	xmm3,[eax+0x10]	/* t10<- t10+t2 */
		__asm	addpd	xmm7,[ebx     ]	/* t14<- t14+t5 */
		__asm	addpd	xmm6,[ebx+0x10]	/* t13<- t13+t6 */

		__asm	movaps	[eax     ],xmm2	/* a[j1   ]=t1+t9  */
		__asm	movaps	[eax+0x10],xmm3	/* a[j2   ]=t2+t10 */
		__asm	movaps	[ebx     ],xmm4	/* a[j1+p2]=t5-t14 */
		__asm	movaps	[ebx+0x10],xmm6	/* a[j2+p2]=t6+t13 */

		__asm	movaps	[ecx     ],xmm0	/* a[j1+p1]=t1-t9  */
		__asm	movaps	[ecx+0x10],xmm1	/* a[j2+p1]=t2-t10 */
		__asm	movaps	[edx     ],xmm7	/* a[j1+p3]=t5+t14 */
		__asm	movaps	[edx+0x10],xmm5	/* a[j2+p3]=t6-t13 */

		/*************************************************/
		/* combine to get #2 of 2 length-4 transforms... */
		/*************************************************/
		/*
		rt =t7;					it =t8;
		t7 =t3 +it;				t8 =t4 -rt;
		t3 =t3 -it;				t4 =t4 +rt;

		rt =t15;				it =t16;
		t15=t11+it;				t16=t12-rt;
		t11=t11-it;				t12=t12+rt;
		*/
		__asm	mov	eax, add4			/*&t3 */
		__asm	mov	ebx, add6			/*&t7 */
		__asm	movaps	xmm0,[eax     ]	/* t3 */
		__asm	movaps	xmm1,[eax+0x10]	/* t4 */
		__asm	movaps	xmm4,xmm0	/* copy of t3 */
		__asm	movaps	xmm5,xmm1	/* copy of t4 */
		/* */
		__asm	subpd	xmm0,[ebx+0x10]	/* t3 =t3 -it */
		__asm	addpd	xmm4,[ebx+0x10]	/* t7 =t3 +it */
		/* */
		__asm	addpd	xmm1,[ebx]		/* t4 =t4 +rt */
		__asm	subpd	xmm5,[ebx]		/* t8 =t4 -rt */
		/* Store copies of t3,4,7,8 into destination array slots: */
		__asm	movaps	[eax     ],xmm0	/* a[j1+p4]=t3 */
		__asm	movaps	[eax+0x10],xmm1	/* a[j2+p4]=t4 */
		__asm	movaps	[ebx     ],xmm4	/* a[j1+p6]=t7 */
		__asm	movaps	[ebx+0x10],xmm5	/* a[j2+p6]=t8 */

		__asm	mov	ecx, add5			/*&t11*/
		__asm	mov	edx, add7			/*&t15*/
		__asm	movaps	xmm2,[ecx     ]	/* t11*/
		__asm	movaps	xmm3,[ecx+0x10]	/* t12*/
		__asm	movaps	xmm6,xmm2	/* copy of t11*/
		__asm	movaps	xmm7,xmm3	/* copy of t12*/
		/* */
		__asm	subpd	xmm2,[edx+0x10]	/* t11=t11-it */
		__asm	addpd	xmm6,[edx+0x10]	/* t15=t11+it */
		/* */
		__asm	addpd	xmm3,[edx]		/* t12=t12+rt */
		__asm	subpd	xmm7,[edx]		/* t16=t12-tt */

		__asm	movaps	[edx     ],xmm4	/* a[j1+p7]=t7 (t7 ends up in xmm6) */
		__asm	movaps	[edx+0x10],xmm5	/* a[j2+p7]=t8 (t8 ends up in xmm3) */

		/* now combine the two half-transforms */
		/*
		rt =(t11-t12)*ISRT2;	it =(t11+t12)*ISRT2;
		a[j1+p4]=t3+rt;			a[j2+p4]=t4+it;
		a[j1+p5]=t3-rt;			a[j2+p5]=t4-it;

		rt =(t15+t16)*ISRT2;	it =(t16-t15)*ISRT2;
		a[j1+p6]=t7-rt;			a[j2+p6]=t8-it;
		a[j1+p7]=t7+rt;			a[j2+p7]=t8+it;
		*/
		/* Use xmm5 as tmp register here, then restore from [edx+0x10]: */
		__asm	mov	eax,isrt2	/* &ISRT2, overwrite add4 pointer value stored in register eax */

		__asm	movaps	xmm5,xmm2	/* copy of t11*/
		__asm	subpd	xmm2,xmm3	/* (t11-t12)  */
		__asm	addpd	xmm5,xmm3	/* (t11+t12)  */
		__asm	mulpd	xmm2,[eax]	/* rt0=(t11-t12)*ISRT2 */
		__asm	mulpd	xmm5,[eax]	/* it0=(t11+t12)*ISRT2 */
		__asm	movaps	xmm3,[edx+0x10]	/* xmm3 now free */

		/* Use xmm4 as tmp register here, then restore from [ecx]: */
		__asm	movaps	xmm4,xmm7	/* copy of t16*/
		__asm	addpd	xmm4,xmm6	/* (t16+t15)  */
		__asm	subpd	xmm7,xmm6	/* (t16-t15)  */
		__asm	mulpd	xmm4,[eax]	/* rt1=(t16+t15)*ISRT2 */
		__asm	mulpd	xmm7,[eax]	/* it1=(t16-t15)*ISRT2 */
		__asm	movaps	xmm6,[edx]	/* xmm6 now free */

		__asm	mov	eax, add4	/* Restore array pointer */

		__asm	subpd	xmm0,xmm2	/* t3 <- t3-rt0 */
		__asm	subpd	xmm1,xmm5	/* t4 <- t4-it0 */
		__asm	subpd	xmm6,xmm4	/* t7 <- t7-rt1 */
		__asm	subpd	xmm3,xmm7	/* t8 <- t8-it1 */

		__asm	addpd	xmm2,[eax     ]	/* rt0<- rt0+t3 */
		__asm	addpd	xmm5,[eax+0x10]	/* it0<- it0+t4 */
		__asm	addpd	xmm4,[ebx     ]	/* rt1<- rt1+t7 */
		__asm	addpd	xmm7,[ebx+0x10]	/* it1<- it1+t8 */

		__asm	movaps	[eax     ],xmm2	/* a[j1+p4]=t3+rt0 */
		__asm	movaps	[eax+0x10],xmm5	/* a[j2+p4]=t4+it0 */
		__asm	movaps	[ebx     ],xmm6	/* a[j1+p6]=t7-rt1 */
		__asm	movaps	[ebx+0x10],xmm3	/* a[j2+p6]=t8-it1 */

		__asm	movaps	[ecx     ],xmm0	/* a[j1+p5]=t3-rt0 */
		__asm	movaps	[ecx+0x10],xmm1	/* a[j2+p5]=t4-it0 */
		__asm	movaps	[edx     ],xmm4	/* a[j1+p7]=t7+rt1 */
		__asm	movaps	[edx+0x10],xmm7	/* a[j2+p7]=t8+it1 */
														/* Totals: 140 load/store, 82 add/subpd, 32 mulpd, 72 address-compute */
	   #endif	/* pure-ASM */

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		add4 = add0+p4;
		add5 = add0+p5;
		add6 = add0+p6;
		add7 = add0+p7;

		SSE2_RADIX8_DIF_TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,isrt2,c1,c2,c3,c4,c5,c6,c7);

	  #endif	/* MSVC or GCC */

	#else	// !USE_SSE2

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and combine to get the 4 length-2 transforms... */
	#if PFETCH
	addr = &a[j1];
	prefetch_p_doubles(addr);
	#endif
		t1 =a[j1   ];					t2 =a[j2   ];				/* r^0*x0 */
		rt =a[j1+p4]*c4 -a[j2+p4]*s4;	it =a[j2+p4]*c4 +a[j1+p4]*s4;		/* r^4*x4 */
		t3 =t1 -rt;				t4 =t2 -it;				/* y1 = x0-r^4*x4 */
		t1 =t1 +rt;				t2 =t2 +it;				/* y0 = x0+r^4*x4 */
	#if PFETCH
	addp = addr+p1;
	prefetch_p_doubles(addp);
	#endif
		t5 =a[j1+p2]*c2 -a[j2+p2]*s2;	t6 =a[j2+p2]*c2 +a[j1+p2]*s2;		/* r^2*x2 */
		rt =a[j1+p6]*c6 -a[j2+p6]*s6;	it =a[j2+p6]*c6 +a[j1+p6]*s6;		/* r^6*x6 */
		t7 =t5 -rt;				t8 =t6 -it;				/* r^2*y3 = r^2*(x2-r^4*x6) */
		t5 =t5 +rt;				t6 =t6 +it;				/* r^2*y2 = r^2*(x2+r^4*x6) */
	#if PFETCH
	addp = addr+p2;
	prefetch_p_doubles(addp);
	#endif
		t9 =a[j1+p1]*c1 -a[j2+p1]*s1;	t10=a[j2+p1]*c1 +a[j1+p1]*s1;		/* r^1*x1 */
		rt =a[j1+p5]*c5 -a[j2+p5]*s5;	it =a[j2+p5]*c5 +a[j1+p5]*s5;		/* r^5*x5 */
		t11=t9 -rt;				t12=t10-it;				/* r^1*y5 = r^1*(x1-r^4*x5) */
		t9 =t9 +rt;				t10=t10+it;				/* r^1*y4 = r^1*(x1+r^4*x5) */
	#if PFETCH
	addp = addr+p3;
	prefetch_p_doubles(addp);
	#endif
		t13=a[j1+p3]*c3 -a[j2+p3]*s3;	t14=a[j2+p3]*c3 +a[j1+p3]*s3;		/* r^3*x3 */
		rt =a[j1+p7]*c7 -a[j2+p7]*s7;	it =a[j2+p7]*c7 +a[j1+p7]*s7;		/* r^7*x7 */
		t15=t13-rt;				t16=t14-it;				/* r^3*y7 = r^3*(x3-r^4*x7) */
		t13=t13+rt;				t14=t14+it;				/* r^3*y6 = r^3*(x3+r^4*x7) */
	#if PFETCH
	addp = addr+p4;
	prefetch_p_doubles(addp);
	#endif
	/* combine to get the 2 length-4 transforms... */
		rt =t5;					it =t6;					/* r^2*y2 */
		t5 =t1 -rt;				t6 =t2 -it;				/* z1 = y0-r^2*y2 */
		t1 =t1 +rt;				t2 =t2 +it;				/* z0 = y0-r^2*y2 */

		rt =t7;					it =t8;					/* r^2*y3 */
		t7 =t3 +it;				t8 =t4 -rt;				/* z3 = y1-r^2*y3*I */
		t3 =t3 -it;				t4 =t4 +rt;				/* z2 = y1+r^2*y3*I */
	#if PFETCH
	addp = addr+p5;
	prefetch_p_doubles(addp);
	#endif
		rt =t13;				it =t14;				/* r^3*y6 */
		t13=t9 -rt;				t14=t10-it;				/* r^1*z5 = r^1*(y4-r^2*y6) */
		t9 =t9 +rt;				t10=t10+it;				/* r^1*z4 = r^1*(y4+r^2*y6) */

		rt =t15;				it =t16;				/* r^3*y7 */
		t15=t11+it;				t16=t12-rt;				/* r^1*z7 = r^1*(y5-r^4*y7*I) */
		t11=t11-it;				t12=t12+rt;				/* r^1*z6 = r^1*(y5+r^4*y7*I) */
	#if PFETCH
	addp = addr+p6;
	prefetch_p_doubles(addp);
	#endif
	/* now combine the two half-transforms */
		a[j1   ]=t1+t9;			a[j2   ]=t2+t10;
		a[j1+p1]=t1-t9;			a[j2+p1]=t2-t10;

		a[j1+p2]=t5-t14;		a[j2+p2]=t6+t13;
		a[j1+p3]=t5+t14;		a[j2+p3]=t6-t13;
	#if PFETCH
	addp = addr+p7;
	prefetch_p_doubles(addp);
	#endif
		rt =(t11-t12)*ISRT2;	it =(t11+t12)*ISRT2;
		a[j1+p4]=t3+rt;			a[j2+p4]=t4+it;
		a[j1+p5]=t3-rt;			a[j2+p5]=t4-it;

		rt =(t15+t16)*ISRT2;	it =(t16-t15)*ISRT2;
		a[j1+p6]=t7-rt;			a[j2+p6]=t8-it;
		a[j1+p7]=t7+rt;			a[j2+p7]=t8+it;
#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */
}

/***************/

/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform a single radix-8 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details.
*/
void radix8_dit_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
{
	static int max_threads = 0;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int i,j,j1,j2,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	int p1,p2,p3,p4,p5,p6,p7;
	double rt,it;
	double re0,im0,re1,im1;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	/* This allows us to align sincos temporaries on 16-byte boundaries (e.g. for SSE2) via pointers:
	In sse2 mode, alloc a complex[32] array (give it a few extra slots to allow ptr to initial element-to-be-used
	to be aligned on 16-byte bdry, and 2 more for the doubled 1/sqrt2 pair),
	and then store doubled pairs of sincos data - for each scalar sincos datum used
	in non-SSE2 mode, we have a double-pair consisting of 2 copies of the same datum, in SSE2 mode.
	This relies on the fact that in our FFT algorithm, the same set of sincos multipliers computed in the outer
	loop below is re-used for multiple (and an even number, in particular) inputs in the inner loop, i.e.
	in SSE2 mode we simply double the inner-loop stride, and process vector DFT inputs in pairs.
	*/
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */
	vec_dbl *c_tmp,*s_tmp;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
	vec_dbl *isrt2,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;	// In || mode, only above base-pointer (shared by all threads) is static:
  #elif defined(COMPILER_TYPE_GCC)
	static vec_dbl *isrt2,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
  #else
	static vec_dbl *r0,*r8;
	static vec_dbl *isrt2,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
  #endif

#else

	double *addr, *addp;	/* addresses for prefetch macros */
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16;
	double c1,c2,c3,c4,c5,c6,c7,s1,s2,s3,s4,s5,s6,s7;

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
		sc_arr = ALLOC_VEC_DBL(sc_arr, 36*max_threads);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 16 16-byte slots of sc_arr for temporaries, next 16 for the doubled sincos twiddles,
	next 1 for doubled 1/sqrt2, plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef MULTITHREAD
		__r0  = sc_ptr;
		isrt2 = sc_ptr + 0x20;
		for(i = 0; i < max_threads; ++i) {
			/* These remain fixed within each per-thread local store: */
			VEC_DBL_INIT(isrt2, ISRT2);
			isrt2 += 36;	/* Move on to next thread's local store */
		}
	#elif defined(COMPILER_TYPE_GCC)
		c0    = sc_ptr + 0x10;
		c4    = sc_ptr + 0x12;
		c2    = sc_ptr + 0x14;
		c6    = sc_ptr + 0x16;
		c1    = sc_ptr + 0x18;
		c5    = sc_ptr + 0x1a;
		c3    = sc_ptr + 0x1c;
		c7    = sc_ptr + 0x1e;
		isrt2 = sc_ptr + 0x20;
		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);
	#else
		r0    = sc_ptr + 0x00;		c0    = sc_ptr + 0x10;
	//	r1    = sc_ptr + 0x01;		s0    = sc_ptr + 0x11;
	/*	r2    = sc_ptr + 0x02;	*/	c4    = sc_ptr + 0x12;
	//	r3    = sc_ptr + 0x03;		s4    = sc_ptr + 0x13;
	/*	r4    = sc_ptr + 0x04;	*/	c2    = sc_ptr + 0x14;
	//	r5    = sc_ptr + 0x05;		s2    = sc_ptr + 0x15;
	/*	r6    = sc_ptr + 0x06;	*/	c6    = sc_ptr + 0x16;
	//	r7    = sc_ptr + 0x07;		s6    = sc_ptr + 0x17;
		r8    = sc_ptr + 0x08;		c1    = sc_ptr + 0x18;
	//	r9    = sc_ptr + 0x09;		s1    = sc_ptr + 0x19;
	/*	ra    = sc_ptr + 0x0a;	*/	c5    = sc_ptr + 0x1a;
	//	rb    = sc_ptr + 0x0b;		s5    = sc_ptr + 0x1b;
	/*	rc    = sc_ptr + 0x0c;	*/	c3    = sc_ptr + 0x1c;
	//	rd    = sc_ptr + 0x0d;		s3    = sc_ptr + 0x1d;
	/*	re    = sc_ptr + 0x0e;	*/	c7    = sc_ptr + 0x1e;
	//	rf    = sc_ptr + 0x0f;		s7    = sc_ptr + 0x1f;
									isrt2 = sc_ptr + 0x20;
		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);
	#endif
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		c0 = __r0 + thr_id*36 + 0x10;
		c4    = c0 + 0x02;
		c2    = c0 + 0x04;
		c6    = c0 + 0x06;
		c1    = c0 + 0x08;
		c5    = c0 + 0x0a;
		c3    = c0 + 0x0c;
		c7    = c0 + 0x0e;
		isrt2 = c0 + 0x10;
	#endif

#endif

/*
!...The radix-8 pass is here. The data are processed in NLOOPS blocks of INCR real elements each,
!   where mloops = product of all preceding radices. The fundamental root of unity for this pass
!   is thus the primitive [(product of all preceding radices)*radix_now] =  [mloops*radix_now]th root,
!   and we get the proper table look-up index by multiplying the BR index by N/[mloops*radix_now] = (incr/2)/radix_now.
*/

	/*  stride between elements processed together. */
	p1 = incr>>3;
	p2 = p1 + p1;
	p3 = p2 + p1;
	p4 = p3 + p1;
	p5 = p4 + p1;
	p6 = p5 + p1;
	p7 = p6 + p1;

	p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
	p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
	p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );

	iroot_prim=(incr >> 4);		/* (incr/2)/radix_now */

	for(m=0; m < nloops; m++)
	{
	/* Here are the needed sincos data - these are processed below in bit-reversed order. */
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* We may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks. */
#if HIACC
	#ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-7] are offset w.r.to the thread-local ptr pair as
		(c0,s0) + 0x[0,8,4,c,2,a,6,e]:
		*/
		c_tmp = c0 + 0x0; s_tmp = c_tmp+1;	/* c0,s0 */
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
		c_tmp = c0 + 0x8; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0x4; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0xc; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0x2; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0xa; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0x6; s_tmp = c_tmp+1;
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
		c_tmp = c0 + 0xe; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c7 =rt;		s7 =it;
	#endif

#else

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c1=t1*t3-t2*t4;	s1=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i + iroot;
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c2=t1*t3-t2*t4;	s2=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c5=t1*t3-t2*t4;	s5=t1*t4+t2*t3;

		t1=c1*c5; t2=c1*s5; t3=s1*c5; t4=s1*s5;
		c4=t1+t4; s4=t2-t3; c6=t1-t4; s6=t2+t3;

		t1=c2*c5; t2=c2*s5; t3=s2*c5; t4=s2*s5;
		c3=t1+t4; s3=t2-t3; c7=t1-t4; s7=t2+t3;
#endif

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 3);

	  for(j = jlo; j < jhi; j += stride)
	  {
		j1 = j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	#ifdef USE_SSE2

	  #if defined(COMPILER_TYPE_MSVC)	/* MSVC-style inline ASM: */

		/* Only need to store addresses of real element of each pair - im offset by 0x10
		(i.e. 2 doubles) from that, due to the [re0,re1,im0,im1] ordering we use in SSE2 mode:
		*/
		add0 = &a[j1];
	   #if 1
		SSE2_RADIX8_DIT_TWIDDLE(add0,p4,p5,p6,p7,r0,c4)
												/*  97 load/store, 82 add/subpd, 32 mulpd, 67 address-compute */
	   #elif 0
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		add4 = add0+p4;
		add5 = add0+p5;
		add6 = add0+p6;
		add7 = add0+p7;

		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r0)
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r8)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r0,r2,r4,r6,r8,ra,rc,re)

	   #elif 0	// if(1) - test out pure-asm version
		__asm	mov	eax, add0
		__asm	mov	ebx, p2		// Can't get these via simple load-one-and-shift-as-needed due to array padding scheme
		__asm	mov	ecx, p4
		__asm	mov	edx, p6
		__asm	shl	ebx,  3
		__asm	shl	ecx,  3
		__asm	shl	edx,  3
		__asm	add	ebx, eax
		__asm	add	ecx, eax
		__asm	add	edx, eax
		SSE2_RADIX4_DIF_4TWIDDLE_A         (r0,c4)	/* 36 load/store, 26 add/subpd, 14 mulpd, 21 address-compute */
		__asm	mov	edi, p1
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p1];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1)	/* 41 load/store, 32 add/subpd, 20 mulpd, 23 address-compute */
		/* Combine the 2 radix-4 subtransforms and write outputs back to main array: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0)		/* 32 load/store, 24 add/subpd,  0 mulpd, 23 address-compute */

												/* Totals:109 load/store, 82 add/subpd, 32 mulpd, 67 address-compute */
												/* Versus 140 load/store, 82 add/subpd, 32 mulpd, 72 address-compute for the original implementation below. */
	   #elif 1

		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		add4 = add0+p4;
		add5 = add0+p5;
		add6 = add0+p6;
		add7 = add0+p7;

		/*** 2nd of the 2 length-4 subtransforms gets done first, due to e.g. t1-+t9 combos in final step: ***/
			/*
			t9 =a[j1+p4];	t10=a[j2+p4];
			t11=a[j1+p5];	t12=a[j2+p5];
			~t11=t9 -t11;	~t12=t10-t12;
			~t9 =t9 +t11;	~t10=t10+t12;
			*/
			__asm	mov	eax, add4
			__asm	mov	ebx, add5
			__asm	movaps	xmm0,[eax]		/* xmm0 <- a[j1+p4] = t9 */
			__asm	movaps	xmm1,[eax+0x10]	/* xmm1 <- a[j2+p4] = t10*/
			__asm	movaps	xmm2,xmm0		/* xmm4 <- copy of t9    */
			__asm	movaps	xmm3,xmm1		/* xmm5 <- copy of t10   */
			__asm	addpd	xmm2,[ebx]		/* xmm2 <- t9  */
			__asm	addpd	xmm3,[ebx+0x10]	/* xmm3 <- t10 */
			__asm	subpd	xmm0,[ebx]		/* xmm0 <- t11 */
			__asm	subpd	xmm1,[ebx+0x10]	/* xmm1 <- t12 */
			/*
			t13=a[j1+p6];	t14=a[j2+p6];
			rt =a[j1+p7];	it =a[j2+p7];
			~t15=t13-t15	~t16=t14-t16
			~t13=t13+t15	~t14=t14+t16
			*/
			__asm	mov	ecx, add6
			__asm	mov	edx, add7
			__asm	movaps	xmm4,[ecx]		/* xmm4 <- a[j1+p6] = t13*/
			__asm	movaps	xmm5,[ecx+0x10]	/* xmm5 <- a[j2+p6] = t14*/
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t13   */
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t14   */
			__asm	addpd	xmm6,[edx]		/* xmm6 <- t13 */
			__asm	addpd	xmm7,[edx+0x10]	/* xmm7 <- t14 */
			__asm	subpd	xmm4,[edx]		/* xmm4 <- t15 */
			__asm	subpd	xmm5,[edx+0x10]	/* xmm5 <- t16 */
			/* Copy t13,14 into main-array slot add6 */
			__asm	movaps	[ecx     ],xmm6
			__asm	movaps	[ecx+0x10],xmm7
			/* Copy t15,16 into main-array slot add7 */
			__asm	movaps	[edx     ],xmm4
			__asm	movaps	[edx+0x10],xmm5

		/** GPRs: ***** SSE Regs: ***** Main array: ********\
		*	eax, add4	xmm0 <- t11		add0 <- unused		*
		*	ebx, add5	xmm1 <- t12		add1 <- unused 		*
		*	ecx, add6	xmm2 <- t9 		add2 <- unused		*
		*	edx, add7	xmm3 <- t10		add3 <- unused 		*
		*				xmm4 <- t15		add4 <- unused		*
		*				xmm5 <- t16		add5 <- unused		*
		*				xmm6 <- t13		add6 <- t13,14		*
		*				xmm7 <- t14		add7 <- t15,16		*
		\***************************************************/

		/*
		rt =t13;	t13=t9 -rt ;	t9 =t9 +rt ;	copies of t13 in add6     , xmm6
		it =t14;	t14=t10-it ;	t10=t10+it ;	copies of t14 in add6+0x10, xmm7

		rt =t15;	t15=t11-t16;	t11=t11+t16;	copies of t15 in add7     , xmm4
					t16=t12+rt ;	t12=t12-rt ;	copies of t16 in add7+0x10, xmm5
		*/
		/* Move outputs t11,12 into a[j1,j2+p5], first doing the addsub and mul by ISRT2: */
		/* Move outputs t15,16 into a[j1,j2+p7], first doing the addsub and mul by ISRT2: */

			__asm	addpd	xmm6,xmm2		/* xmm6 <- ~t9  */
			__asm	addpd	xmm7,xmm3		/* xmm7 <- ~t10 */
			__asm	subpd	xmm2,[ecx     ]	/* xmm2 <- ~t13 */
			__asm	subpd	xmm3,[ecx+0x10]	/* xmm3 <- ~t14 */
			/* Move t13,14 into a[j1,j2+p6] */
			__asm	movaps	[ecx     ],xmm2	/* add6r <- ~t13 */
			__asm	movaps	[ecx+0x10],xmm3	/* add6i <- ~t14 */

			__asm	movaps	xmm2,xmm4	/* xmm2 <- copy of t15 */
			__asm	movaps	xmm3,xmm5	/* xmm3 <- copy of t16 */
			__asm	addpd	xmm5,xmm0	/* xmm5 <-~t11 */
			__asm	subpd	xmm0,xmm3	/* xmm0 <-~t15 */
			__asm	addpd	xmm4,xmm1	/* xmm4 <-~t16 */
			__asm	subpd	xmm1,xmm2	/* xmm1 <-~t12 */

			__asm	movaps	xmm2,xmm5	/* xmm2 <- copy of~t11 */
			__asm	movaps	xmm3,xmm1	/* xmm3 <- copy of~t12 */
			__asm	addpd	xmm5,xmm1	/* xmm5 <-~(t11+t12), xmm1 FREE */
			__asm	mov	ecx,isrt2
			__asm	movaps	xmm1,[ecx]	/* xmm1 <- ISRT2 */
			__asm	subpd	xmm2,xmm3	/* xmm2 <-~(t11-t12), xmm3 FREE */
			__asm	mulpd	xmm5,xmm1	/* xmm5 <- (t11+t12)*ISRT2 */
			__asm	mulpd	xmm2,xmm1	/* xmm2 <- (t11-t12)*ISRT2 */
			__asm	movaps	xmm3,xmm0	/* xmm3 <- copy of~t15 */

			__asm	movaps	[ebx     ],xmm5	/* add5r<- (t11+t12)*ISRT2, xmm5 FREE */
			__asm	movaps	xmm5,xmm4		/* xmm5 <- copy of~t16 */
			__asm	addpd	xmm0,xmm4	/* xmm0 <-~(t15+t16) */
			__asm	movaps	[ebx+0x10],xmm2	/* add5i<- (t11-t12)*ISRT2 */
			__asm	subpd	xmm3,xmm5	/* xmm3 <-~(t15-t16) */
			__asm	mulpd	xmm0,xmm1	/* xmm0 <- (t15+t16)*ISRT2 */
			__asm	mulpd	xmm3,xmm1	/* xmm3 <- (t15-t16)*ISRT2 */
			__asm	movaps	[edx     ],xmm0	/* add7r<- (t15+t16)*ISRT2 */
			__asm	movaps	[edx+0x10],xmm3	/* add7i<- (t15-t16)*ISRT2 */

		/** GPRs: ***** SSE Regs: ***** Main array: ******************\
		*    eax, add4   xmm0 unused        add0 <- unused            *
		*    ebx, add5   xmm1 unused        add1 <- unused            *
		*    ecx, isrt2  xmm2 unused        add2 <- unused            *
		*    edx, add7   xmm3 unused        add3 <- unused            *
		*                xmm4 unused        add4 <- unused            *
		*                xmm5 unused        add5 <- (t11+-t12)*ISRT2  *
		*                xmm6 t9            add6 <-  t13,14           *
		*                xmm7 t10           add7 <- (t15+-t16)*ISRT2  *
		\*************************************************************/

		/**************** 1st of the 2 length-4 subtransforms... **************/
		/*
		t1 =a[j1   ];	t2 =a[j2   ];
		rt =a[j1+p1];	it =a[j2+p1];
		t3 =t1 -rt;		t4 =t2 -it;
		t1 =t1 +rt;		t2 =t2 +it;
		*/
		__asm	mov	eax, add0
		__asm	mov	ebx, add1

		__asm	movaps	xmm0,[eax]		/* xmm0 <- a[j1   ] = t1 */
		__asm	movaps	xmm1,[eax+0x10]	/* xmm1 <- a[j2   ] = t2 */
		__asm	movaps	xmm2,xmm0		/* xmm2 <- copy of t1    */
		__asm	movaps	xmm3,xmm1		/* xmm3 <- copy of t2    */
		__asm	addpd	xmm2,[ebx]		/*~xmm2 <- t1  */
		__asm	addpd	xmm3,[ebx+0x10]	/*~xmm3 <- t2  */
		__asm	subpd	xmm0,[ebx]		/*~xmm0 <- t3  */
		__asm	subpd	xmm1,[ebx+0x10]	/*~xmm1 <- t4  */

		/* Move   t9,  t10 into a[j1,j2   ] in anticipation of final outputs t1+-t9, t2+-t10 which will go there: */
		__asm	movaps	[eax     ],xmm6
		__asm	movaps	[eax+0x10],xmm7	/* add0 <-  t9,t10 */
		__asm	addpd	xmm6,xmm6		/* xmm6 <- 2*t9  */
		__asm	addpd	xmm7,xmm7		/* xmm7 <- 2*t10 */
		/* Move 2*t9,2*t10 into a[j1,j2+p4] in anticipation of final outputs t1+-t9, t2+-t10 which will go there: */
		__asm	mov	ecx, add4
		__asm	movaps	[ecx     ],xmm6	/* add4 <- 2*t9,2*t10  */
		__asm	movaps	[ecx+0x10],xmm7	/* xmm4-7 FREE */
		/*
		t5 =a[j1+p2];	t6 =a[j2+p2];
		rt =a[j1+p3];	it =a[j2+p3];
		t7 =t5 -rt;		t8 =t6 -it;
		t5 =t5 +rt;		t6 =t6 +it;
		*/
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
		__asm	movaps	xmm4,[ecx]		/* xmm4 <- a[j1+p2] = t5 */
		__asm	movaps	xmm5,[ecx+0x10]	/* xmm5 <- a[j2+p2] = t6 */
		__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t5    */
		__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t6    */
		__asm	addpd	xmm6,[edx]		/*~xmm6 <- t5  */
		__asm	addpd	xmm7,[edx+0x10]	/*~xmm7 <- t6  */
		__asm	subpd	xmm4,[edx]		/*~xmm4 <- t7  */
		__asm	subpd	xmm5,[edx+0x10]	/*~xmm5 <- t8  */
		/* Copy t5,6 into main-array slots a[j1,j2+p2] */
		__asm	movaps	[ecx     ],xmm6
		__asm	movaps	[ecx+0x10],xmm7	/* add2 <-  t5,t6 */

		/** GPRs: ***** SSE Regs: ***** Main array: ******************\
		*    eax, add0   xmm0 t3            add0 <-  t9,t10           *
		*    ebx, add1   xmm1 t4            add1 <- unused            *
		*    ecx, add2   xmm2 t1            add2 <-  t5,t6            *
		*    edx, add3   xmm3 t2            add3 <- unused            *
		*                xmm4 t7            add4 <- 2*t9,2*t10        *
		*                xmm5 t8            add5 <- (t11+-t12)*ISRT2  *
		*                xmm6 t5            add6 <-  t13,14           *
		*                xmm7 t6            add7 <- (t15+-t16)*ISRT2  *
		\*************************************************************/

		/*
		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;	copies of t5 in add2     , xmm6
		it =t6;	t6 =t2 -it;	t2 =t2 +it;	copies of t6 in add2+0x10, xmm7

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
				t8 =t4 +rt;	t4 =t4 -rt;
		*/
		__asm	addpd	xmm6,xmm2		/* xmm6 <- ~t1 */
		__asm	addpd	xmm7,xmm3		/* xmm7 <- ~t2 */
		__asm	subpd	xmm2,[ecx     ]	/* xmm2 <- ~t5 */
		__asm	subpd	xmm3,[ecx+0x10]	/* xmm3 <- ~t6 */

		/* Compute and dump first 2 outputs now, in order to free up 2 registers: */
		__asm	addpd	xmm6,[eax     ]	/* t1+t9 */
		__asm	addpd	xmm7,[eax+0x10]	/* t2+t10*/
		__asm	movaps	[eax     ],xmm6	/* a[j1   ], DONE. */
		__asm	movaps	[eax+0x10],xmm7	/* a[j2   ], DONE. */

		__asm	mov	ebx, add4
		__asm	subpd	xmm6,[ebx     ]	/* t1-t9  = [t1+t9 ] - 2*t9  */
		__asm	subpd	xmm7,[ebx+0x10]	/* t2-t10 = [t2+t10] - 2*t10 */
		__asm	movaps	[ebx     ],xmm6	/* Spill t1-t9  to a[j1+p4]. */
		__asm	movaps	[ebx+0x10],xmm7	/* Spill t2-t10 to a[j2+p4]. */
		__asm	movaps	[edx     ],xmm6	/* Copy: t1-t9  -> a[j1+p3]. */
		__asm	movaps	[edx+0x10],xmm7	/* Copy: t2-t10 -> a[j2+p3]. */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- copy of t7 */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- copy of t8 */
		__asm	addpd	xmm5,xmm0	/* xmm5 <- ~t3 */
		__asm	subpd	xmm0,xmm7	/* xmm0 <- ~t7 */
		__asm	addpd	xmm4,xmm1	/* xmm4 <- ~t8 */
		__asm	subpd	xmm1,xmm6	/* xmm1 <- ~t4 */

/** GPRs: ***** SSE Regs: ***** Main array: ******************\
*    eax, add0   xmm0 t7            add0 <- DONE              *
*    ebx, add4   xmm1 t4            add1 <- unused            *
*    ecx, ----   xmm2 t5            add2 <- unused            *
*    edx, add3   xmm3 t6            add3 <- t1-t9,t2-t10      *
*                xmm4 t8            add4 <- t1-t9,t2-t10      *
*                xmm5 t3            add5 <- (t11+-t12)*ISRT2  *
*                xmm6 unused        add6 <-  t13,14           *
*                xmm7 unused        add7 <- (t15+-t16)*ISRT2  *
\*************************************************************/

/* Now combine the two half-transforms & store outputs back into original array slots: */
		/*
		a[j1   ]=t1+t9;				a[j2   ]=t2+t10;	already done
		~t1     =t1-t9;				~t2     =t2-t10;	copies in add3[edx],add4[ebx]
		a[j1+p4]=~t1*c4 +~t2*s4;	a[j2+p4]=~t2*c4 -~t1*s4;
		*/
		__asm	mov	ecx,   c4	/* &c4 */
		__asm	movaps	xmm6,[ebx     ]	/* t1-t9  <- a[j1+p4]. */
		__asm	movaps	xmm7,[ebx+0x10]	/* t2-t10 <- a[j2+p4]. */
		__asm	mulpd	xmm6,[ecx+0x10]	/*~t1*s4 */
		__asm	mulpd	xmm7,[ecx+0x10]	/*~t2*s4 */
		__asm	movaps	[ebx     ],xmm6	/* a[j1+p4] <- ~t1*s4 */
		__asm	movaps	[ebx+0x10],xmm7	/* a[j2+p4] <- ~t2*s4 */

		__asm	movaps	xmm6,[edx     ]	/* xmm6 <- cpy ~t1 */
		__asm	movaps	xmm7,[edx+0x10]	/* xmm7 <- cpy ~t2 */
		__asm	mulpd	xmm6,[ecx     ]	/*~t1*c4 */
		__asm	mulpd	xmm7,[ecx     ]	/*~t2*c4 */
		__asm	addpd	xmm6,[ebx+0x10]	/* ~t1*c4 +~t2*s4 */
		__asm	subpd	xmm7,[ebx     ]	/* ~t2*c4 -~t1*s4 */

		__asm	movaps	[ebx     ],xmm6	/* a[j1+p4] <- ~t1*c4 +~t2*s4 */
		__asm	movaps	[ebx+0x10],xmm7	/* a[j2+p4] <- ~t2*c4 -~t1*s4 */

		/*
		rt=(t11+t12)*ISRT2;			it=(t11-t12)*ISRT2;	precomputed, in add5
		t11     =t3+rt;				t12       =t4-it;
		t3      =t3-rt;				t4        =t4+it;
		a[j1+p1]=t11*c1 +t12*s1;	a[j2+p1]=t12*c1 -t11*s1;
		a[j1+p5]=t3 *c5 +t4 *s5;	a[j2+p5]=t4 *c5 -t3 *s5;
		*/
		__asm	mov	eax, add1
		__asm	mov	ebx, add5
		__asm	mov	ecx,   c1
		__asm	mov	edx,   c5
		__asm	movaps	xmm6,xmm5		/* xmm6 <- copy of t3 */
		__asm	movaps	xmm7,xmm1		/* xmm7 <- copy of t4 */
		__asm	addpd	xmm5,[ebx     ]	/* ~t11 */
		__asm	subpd	xmm1,[ebx+0x10]	/* ~t12 */
		__asm	subpd	xmm6,[ebx     ]	/* ~t3  */
		__asm	addpd	xmm7,[ebx+0x10]	/* ~t4  */
		__asm	movaps	[ebx     ],xmm6	/* spill ~t3,t4 to add5 */
		__asm	movaps	[ebx+0x10],xmm7
		__asm	movaps	xmm6,xmm5	/* xmm6 <- cpy ~t11*/
		__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy ~t12*/

		__asm	mulpd	xmm5,[ecx     ]	/* t11*c1 */
		__asm	mulpd	xmm1,[ecx     ]	/* t12*c1 */
		__asm	mulpd	xmm6,[ecx+0x10]	/* t11*s1 */
		__asm	mulpd	xmm7,[ecx+0x10]	/* t12*s1 */
		__asm	subpd	xmm1,xmm6	/* t12*c1 - t11*s1 */
		__asm	addpd	xmm5,xmm7	/* t11*c1 + t12*s1 */
		__asm	movaps	[eax+0x10],xmm1	/* a[j2+p1] */
		__asm	movaps	[eax     ],xmm5	/* a[j1+p1] */

		__asm	movaps	xmm5,[ebx     ]	/* t3 */
		__asm	movaps	xmm1,[ebx+0x10]	/* t4 */
		__asm	movaps	xmm6,xmm5	/* xmm6 <- copy t3 */
		__asm	movaps	xmm7,xmm1	/* xmm7 <- copy t4 */

		__asm	mulpd	xmm5,[edx     ]	/* t3 *c5 */
		__asm	mulpd	xmm1,[edx     ]	/* t4 *c5 */
		__asm	mulpd	xmm6,[edx+0x10]	/* t3 *s5 */
		__asm	mulpd	xmm7,[edx+0x10]	/* t4 *s5 */
		__asm	subpd	xmm1,xmm6	/* t4*c5 - t3*s5 */
		__asm	addpd	xmm5,xmm7	/* t3*c5 + t4*s5 */
		__asm	movaps	[ebx+0x10],xmm1	/* a[j2+p5] */
		__asm	movaps	[ebx     ],xmm5	/* a[j1+p5] */	/*** xmm1,5,6,7 FREE ***/

		/*
		rt      =t5+t14;			it        =t6-t13;
		t5      =t5-t14;			t6        =t6+t13;
		a[j1+p2]=rt *c2 +it *s2;	a[j2+p2]=it *c2 -rt *s2;
		a[j1+p6]=t5 *c6 +t6 *s6;	a[j2+p6]=t6 *c6 -t5 *s6;
		*/
		__asm	mov	eax, add2
		__asm	mov	ebx, add6
		__asm	mov	ecx,   c2
		__asm	mov	edx,   c6
		__asm	movaps	xmm6,xmm2		/* xmm6 <- copy of t5 */
		__asm	movaps	xmm7,xmm3		/* xmm7 <- copy of t6 */
		__asm	addpd	xmm2,[ebx+0x10]	/*  rt */
		__asm	subpd	xmm3,[ebx     ]	/*  it */
		__asm	subpd	xmm6,[ebx+0x10]	/* ~t5  */
		__asm	addpd	xmm7,[ebx     ]	/* ~t6  */
		__asm	movaps	xmm1,xmm2	/* xmm1 <- cpy rt */
		__asm	movaps	xmm5,xmm3	/* xmm5 <- cpy it */

		__asm	mulpd	xmm2,[ecx     ]	/* rt*c2 */
		__asm	mulpd	xmm3,[ecx     ]	/* it*c2 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* rt*s2 */
		__asm	mulpd	xmm5,[ecx+0x10]	/* it*s2 */
		__asm	subpd	xmm3,xmm1	/* it*c2 - rt*s2 */
		__asm	addpd	xmm2,xmm5	/* rt*c2 + it*s2 */
		__asm	movaps	[eax+0x10],xmm3	/* a[j2+p2] */
		__asm	movaps	[eax     ],xmm2	/* a[j1+p2] */

		__asm	movaps	xmm1,xmm6	/* xmm1 <- cpy t5 */
		__asm	movaps	xmm5,xmm7	/* xmm5 <- cpy t6 */

		__asm	mulpd	xmm6,[edx     ]	/* t5*c6 */
		__asm	mulpd	xmm7,[edx     ]	/* t6*c6 */
		__asm	mulpd	xmm1,[edx+0x10]	/* t5*s6 */
		__asm	mulpd	xmm5,[edx+0x10]	/* t6*s6 */
		__asm	subpd	xmm7,xmm1	/* t6*c6 - t5*s6 */
		__asm	addpd	xmm6,xmm5	/* t5*c6 + t6*s6 */
		__asm	movaps	[ebx+0x10],xmm7	/* a[j2+p6] */
		__asm	movaps	[ebx     ],xmm6	/* a[j1+p6] */

		/*
		rt=(t15-t16)*ISRT2;			it=(t15+t16)*ISRT2;	precomputed, it,rt] in add7; NOTE reversed order!
		t15     =t7-rt;				t16       =t8-it;
		t7      =t7+rt;				t8        =t8+it;
		a[j1+p3]=t15*c3 +t16*s3;	a[j2+p3]=t16*c3 -t15*s3;
		a[j1+p7]=t7 *c7 +t8 *s7;	a[j2+p7]=t8 *c7 -t7 *s7;
		*/
		__asm	mov	eax, add3
		__asm	mov	ebx, add7
		__asm	mov	ecx,   c3
		__asm	mov	edx,   c7
		__asm	movaps	xmm6,xmm0		/* xmm6 <- copy of t7 */
		__asm	movaps	xmm7,xmm4		/* xmm7 <- copy of t8 */
		__asm	subpd	xmm0,[ebx+0x10]	/* ~t15 */
		__asm	subpd	xmm4,[ebx     ]	/* ~t16 */
		__asm	addpd	xmm6,[ebx+0x10]	/* ~t7  */
		__asm	addpd	xmm7,[ebx     ]	/* ~t8  */

		__asm	movaps	xmm1,xmm0	/* xmm6 <- cpy t15*/
		__asm	movaps	xmm5,xmm4	/* xmm7 <- cpy t16*/

		__asm	mulpd	xmm0,[ecx     ]	/* t15*c3 */
		__asm	mulpd	xmm4,[ecx     ]	/* t16*c3 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* t15*s3 */
		__asm	mulpd	xmm5,[ecx+0x10]	/* t16*s3 */
		__asm	subpd	xmm4,xmm1	/* t16*c3 -t15*s3 */
		__asm	addpd	xmm0,xmm5	/* t15*c3 +t16*s3 */
		__asm	movaps	[eax+0x10],xmm4	/* a[j2+p3] */
		__asm	movaps	[eax     ],xmm0	/* a[j1+p3] */

		__asm	movaps	xmm1,xmm6	/* xmm1 <- cpy t7 */
		__asm	movaps	xmm5,xmm7	/* xmm5 <- cpy t8 */

		__asm	mulpd	xmm6,[edx     ]	/* t7*c7 */
		__asm	mulpd	xmm7,[edx     ]	/* t8*c7 */
		__asm	mulpd	xmm1,[edx+0x10]	/* t7*s7 */
		__asm	mulpd	xmm5,[edx+0x10]	/* t8*s7 */
		__asm	subpd	xmm7,xmm1	/* t8*c7 - t7*s7 */
		__asm	addpd	xmm6,xmm5	/* t7*c7 + t8*s7 */
		__asm	movaps	[ebx+0x10],xmm7	/* a[j2+p7] */
		__asm	movaps	[ebx     ],xmm6	/* a[j1+p7] */
												/* Totals:121 load/store [64 movaps], 68 add/subpd, 32 mulpd, 60 address-compute */
														/* 42 fewer movaps than DIF! [I.e. it pays to do careful register mgmt] */
	   #endif	/* #if 1 */

	  #elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		add4 = add0+p4;
		add5 = add0+p5;
		add6 = add0+p6;
		add7 = add0+p7;

		SSE2_RADIX8_DIT_TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,isrt2,c1,c2,c3,c4,c5,c6,c7);

	  #endif	/* MSVC or GCC */

#else	// USE_SSE2 ?

	/*      gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and get the 4 length-2 transforms... */
	#if PFETCH
	addr = &a[j1];
	prefetch_p_doubles(addr);
	#endif
		t1 =a[j1   ];	t2 =a[j2   ];		/* r^0*x0 */
		rt =a[j1+p1];	it =a[j2+p1];		/* r^4*x4 */
		t3 =t1 -rt;		t4 =t2 -it;		/* r^0*y1 = r^0*x0 - r^4*x4 */
		t1 =t1 +rt;		t2 =t2 +it;		/* r^0*y0 = r^0*x0 + r^4*x4 */
	#if PFETCH
	addp = addr+p1;
	prefetch_p_doubles(addp);
	#endif
		t5 =a[j1+p2];	t6 =a[j2+p2];		/* r^2*x2 */
		rt =a[j1+p3];	it =a[j2+p3];		/* r^6*x6 */
		t7 =t5 -rt;		t8 =t6 -it;		/* r^2*y3 = r^2*x2 - r^6*x6 */
		t5 =t5 +rt;		t6 =t6 +it;		/* r^2*y2 = r^2*x2 + r^6*x6 */
	#if PFETCH
	addp = addr+p2;
	prefetch_p_doubles(addp);
	#endif
		t9 =a[j1+p4];	t10=a[j2+p4];		/* r^1*x1 */
		rt =a[j1+p5];	it =a[j2+p5];		/* r^5*x5 */
		t11=t9 -rt;		t12=t10-it;		/* r^1*y5 = r^1*x1 - r^5*x5 */
		t9 =t9 +rt;		t10=t10+it;		/* r^1*y4 = r^1*x1 + r^5*x5 */
	#if PFETCH
	addp = addr+p3;
	prefetch_p_doubles(addp);
	#endif
		t13=a[j1+p6];	t14=a[j2+p6];		/* r^3*x3 */
		rt =a[j1+p7];	it =a[j2+p7];		/* r^7*x7 */
		t15=t13-rt;		t16=t14-it;		/* r^3*y7 = r^3*x3 - r^7*x7 */
		t13=t13+rt;		t14=t14+it;		/* r^3*y6 = r^3*x3 + r^7*x7 */
	#if PFETCH
	addp = addr+p4;
	prefetch_p_doubles(addp);
	#endif
/*      combine to get the 2 length-4 transform... */
		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
			t8 =t4 +rt;	t4 =t4 -rt;
	#if PFETCH
	addp = addr+p5;
	prefetch_p_doubles(addp);
	#endif
		rt =t13;	t13=t9 -rt ;	t9 =t9 +rt ;
		it =t14;	t14=t10-it ;	t10=t10+it ;

		rt =t15;	t15=t11-t16;	t11=t11+t16;
					t16=t12+rt ;	t12=t12-rt ;
	#if PFETCH
	addp = addr+p6;
	prefetch_p_doubles(addp);
	#endif
/*      now combine the two half-transforms */
		a[j1   ]=t1+t9;				a[j2   ]=t2+t10;
		t1      =t1-t9;				t2      =t2-t10;
		a[j1+p4]=t1 *c4 +t2 *s4;	a[j2+p4]=t2 *c4 -t1 *s4;

		rt=(t11+t12)*ISRT2;			it=(t11-t12)*ISRT2;
		t11     =t3+rt;				t12       =t4-it;
		t3      =t3-rt;				t4        =t4+it;
		a[j1+p1]=t11*c1 +t12*s1;	a[j2+p1]=t12*c1 -t11*s1;
		a[j1+p5]=t3 *c5 +t4 *s5;	a[j2+p5]=t4 *c5 -t3 *s5;
	#if PFETCH
	addp = addr+p7;
	prefetch_p_doubles(addp);
	#endif
		rt      =t5+t14;			it        =t6-t13;
		t5      =t5-t14;			t6        =t6+t13;
		a[j1+p2]=rt *c2 +it *s2;	a[j2+p2]=it *c2 -rt *s2;
		a[j1+p6]=t5 *c6 +t6 *s6;	a[j2+p6]=t6 *c6 -t5 *s6;

		rt=(t15-t16)*ISRT2;			it=(t15+t16)*ISRT2;
		t15     =t7-rt;				t16       =t8-it;
		t7      =t7+rt;				t8        =t8+it;
		a[j1+p3]=t15*c3 +t16*s3;	a[j2+p3]=t16*c3 -t15*s3;
		a[j1+p7]=t7 *c7 +t8 *s7;	a[j2+p7]=t8 *c7 -t7 *s7;
#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

