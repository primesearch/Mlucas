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

/* Use for testing higher-accuracy version of the twiddles computation */
#define HIACC 1

#ifdef USE_SSE2

	#if(HIACC != 1)
		#error SIMD Mode requires HIACC flag to be set!
	#endif

	#include "radix8_dif_dit_pass_asm.h"

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
#ifdef USE_SSE2
	static int max_threads = 0;
#endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int i,j,j1,jlo,jhi,m,iroot_prim,iroot,k1,k2;
#ifndef USE_SSE2
	int j2;
#endif
	int p1,p2,p3,p4,p5,p6,p7;
#ifdef USE_SSE2
	double dtmp;
#endif
	double rt,it;
	double re0,im0,re1,im1;

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error SSE2 code not supported for this compiler!
  #endif

	// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
	// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
	const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;

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
	vec_dbl *sqrt2,*isrt2,*one,*two,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;	// In || mode, only above base-pointer (shared by all threads) is static:
  #elif defined(COMPILER_TYPE_GCC)
	static vec_dbl *sqrt2,*isrt2,*one,*two,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
  #else
	static vec_dbl *r0,*r8;
	static vec_dbl *sqrt2,*isrt2,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
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
		ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(thr_id == -1, "Init-mode call must be outside of any multithreading!");
		sc_arr = ALLOC_VEC_DBL(sc_arr, 36*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 16 16-byte slots of sc_arr for temporaries, next 16 for the doubled sincos twiddles,
	next 1 for doubled 1/sqrt2, plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef MULTITHREAD
		__r0  = sc_ptr;
		sqrt2 = sc_ptr + 0x20;	// Need 16 local vec_dbl-sized slots for tmp-data, 16 for sincos twiddles, 4 for [sqrt2,isrt2,1,0,2.0] quartet
		isrt2 = sc_ptr + 0x21;
		one   = sc_ptr + 0x22;
		two   = sc_ptr + 0x23;
		for(i = 0; i < max_threads; ++i) {
			/* These remain fixed within each per-thread local store: */
		  #if 1
			// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
			dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
			dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
		  #else
			VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
		  #endif
			VEC_DBL_INIT(one  , 1.0);
			VEC_DBL_INIT(two  , 2.0);
			isrt2 += 36;	/* Move on to next thread's local store */
			one   += 36;
			two   += 36;
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
		sqrt2 = sc_ptr + 0x20;	// Need 16 local vec_dbl-sized slots for tmp-data, 16 for sincos twiddles, 4 for [sqrt2,isrt2,1,0,2.0] quartet
		isrt2 = sc_ptr + 0x21;
		one   = sc_ptr + 0x22;
		two   = sc_ptr + 0x23;
		/* These remain fixed: */
	  #if 1
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
	  #endif
		VEC_DBL_INIT(one  , 1.0);
		VEC_DBL_INIT(two  , 2.0);
	#else
		#error Non-GCC-compatible compilers not supported for SIMD builds!
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
									isrt2 = sc_ptr + 0x21;
		/* These remain fixed: */
	  #if 1
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
	  #endif
	#endif
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
  #ifdef MULTITHREAD
	ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
	c0 = __r0 + thr_id*36 + 0x10;
	c4    = c0 + 0x02;
	c2    = c0 + 0x04;
	c6    = c0 + 0x06;
	c1    = c0 + 0x08;
	c5    = c0 + 0x0a;
	c3    = c0 + 0x0c;
	c7    = c0 + 0x0e;
	sqrt2 = c0 + 0x10;	// Need 16 local vec_dbl-sized slots for tmp-data, 16 for sincos twiddles, 4 for [sqrt2,isrt2,1,0,2.0] quartet
	isrt2 = c0 + 0x11;
	one   = c0 + 0x12;
	two   = c0 + 0x13;
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
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	#ifdef USE_SSE2

		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		add4 = add0+p4;
		add5 = add0+p5;
		add6 = add0+p6;
		add7 = add0+p7;

		SSE2_RADIX8_DIF_TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,isrt2,c1,c2,c3,c4,c5,c6,c7);

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
#ifdef USE_SSE2
	static int max_threads = 0;
#endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int i,j,j1,jlo,jhi,m,iroot_prim,iroot,k1,k2;
#ifndef USE_SSE2
	int j2;
#endif
	int p1,p2,p3,p4,p5,p6,p7;
#ifdef USE_SSE2
	double dtmp;
#endif
	double rt,it;
	double re0,im0,re1,im1;

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error SSE2 code not supported for this compiler!
  #endif

	// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
	// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
	const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;

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
	vec_dbl *sqrt2,*isrt2,*one,*two,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;	// In || mode, only above base-pointer (shared by all threads) is static:
  #elif defined(COMPILER_TYPE_GCC)
	static vec_dbl *sqrt2,*isrt2,*one,*two,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
  #else
	static vec_dbl *r0,*r8;
	static vec_dbl *sqrt2,*isrt2,*c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
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
		ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(thr_id == -1, "Init-mode call must be outside of any multithreading!");
		sc_arr = ALLOC_VEC_DBL(sc_arr, 36*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 16 16-byte slots of sc_arr for temporaries, next 16 for the doubled sincos twiddles,
	next 1 for doubled 1/sqrt2, plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef MULTITHREAD
		__r0  = sc_ptr;
		sqrt2 = sc_ptr + 0x20;	// Need 16 local vec_dbl-sized slots for tmp-data, 16 for sincos twiddles, 4 for [sqrt2,isrt2,1,0,2.0] quartet
		isrt2 = sc_ptr + 0x21;
		one   = sc_ptr + 0x22;
		two   = sc_ptr + 0x23;
		for(i = 0; i < max_threads; ++i) {
			/* These remain fixed within each per-thread local store: */
		  #if 1
			// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
			dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
			dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
		  #else
			VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
		  #endif
			VEC_DBL_INIT(one  , 1.0);
			VEC_DBL_INIT(two  , 2.0);
			isrt2 += 36;	/* Move on to next thread's local store */
			one   += 36;
			two   += 36;
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
		sqrt2 = sc_ptr + 0x20;	// Need 16 local vec_dbl-sized slots for tmp-data, 16 for sincos twiddles, 4 for [sqrt2,isrt2,1,0,2.0] quartet
		isrt2 = sc_ptr + 0x21;
		one   = sc_ptr + 0x22;
		two   = sc_ptr + 0x23;
		/* These remain fixed: */
	  #if 1
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
	  #endif
		VEC_DBL_INIT(one  , 1.0);
		VEC_DBL_INIT(two  , 2.0);
	#else
		#error Non-GCC-compatible compilers not supported for SIMD builds!
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
		sqrt2 = sc_ptr + 0x20;	// Need 16 local vec_dbl-sized slots for tmp-data, 16 for sincos twiddles, 4 for [sqrt2,isrt2,1,0,2.0] quartet
		isrt2 = sc_ptr + 0x21;
		one   = sc_ptr + 0x22;
		two   = sc_ptr + 0x23;
		/* These remain fixed: */
	  #if 1
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
	  #endif
	#endif
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
  #ifdef MULTITHREAD
	ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
	c0 = __r0 + thr_id*36 + 0x10;
	c4    = c0 + 0x02;
	c2    = c0 + 0x04;
	c6    = c0 + 0x06;
	c1    = c0 + 0x08;
	c5    = c0 + 0x0a;
	c3    = c0 + 0x0c;
	c7    = c0 + 0x0e;
	sqrt2 = c0 + 0x10;	// Need 16 local vec_dbl-sized slots for tmp-data, 16 for sincos twiddles, 4 for [sqrt2,isrt2,1,0,2.0] quartet
	isrt2 = c0 + 0x11;
	one   = c0 + 0x12;
	two   = c0 + 0x13;
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
/******************* AVX debug stuff: *******************/
#if 0
if(1) {
	int ipad;
	// Use RNG to populate data array:
	rng_isaac_init(TRUE);
	double dtmp = 1024.0*1024.0*1024.0*1024.0;
	j1 = jlo;
	j1 += ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
  #ifdef USE_SSE2
	ipad =  0;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[j1+ipad+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p1;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[j1+ipad+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p2;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[j1+ipad+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p3;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[j1+ipad+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p4;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[j1+ipad+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p5;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[j1+ipad+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p6;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[j1+ipad+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p7;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[j1+ipad+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
  #else
	ipad =  0;	for(i = 0; i < 4; i++) { a[j1+ipad+     i ] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p1;	for(i = 0; i < 4; i++) { a[j1+ipad+     i ] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p2;	for(i = 0; i < 4; i++) { a[j1+ipad+     i ] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p3;	for(i = 0; i < 4; i++) { a[j1+ipad+     i ] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p4;	for(i = 0; i < 4; i++) { a[j1+ipad+     i ] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p5;	for(i = 0; i < 4; i++) { a[j1+ipad+     i ] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p6;	for(i = 0; i < 4; i++) { a[j1+ipad+     i ] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	ipad = p7;	for(i = 0; i < 4; i++) { a[j1+ipad+     i ] = dtmp*rng_isaac_rand_double_norm_pm1(); }
  #endif
  #if 0
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 0, a[ipad+br16[ 0]],ipad+ 0, a[ipad+ 0]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 1, a[ipad+br16[ 1]],ipad+ 1, a[ipad+ 1]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 2, a[ipad+br16[ 2]],ipad+ 2, a[ipad+ 2]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 3, a[ipad+br16[ 3]],ipad+ 3, a[ipad+ 3]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 4, a[ipad+br16[ 4]],ipad+ 4, a[ipad+ 4]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 5, a[ipad+br16[ 5]],ipad+ 5, a[ipad+ 5]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 6, a[ipad+br16[ 6]],ipad+ 6, a[ipad+ 6]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 7, a[ipad+br16[ 7]],ipad+ 7, a[ipad+ 7]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 8, a[ipad+br16[ 8]],ipad+ 8, a[ipad+ 8]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+ 9, a[ipad+br16[ 9]],ipad+ 9, a[ipad+ 9]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+10, a[ipad+br16[10]],ipad+10, a[ipad+10]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+11, a[ipad+br16[11]],ipad+11, a[ipad+11]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+12, a[ipad+br16[12]],ipad+12, a[ipad+12]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+13, a[ipad+br16[13]],ipad+13, a[ipad+13]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+14, a[ipad+br16[14]],ipad+14, a[ipad+14]);
	printf("A_in[%2d] = %20.10e; SIMD: A_in[%2d] = %20.10e\n",ipad+15, a[ipad+br16[15]],ipad+15, a[ipad+15]);
  #endif
//exit(0);
}
#endif
/********************************************************/

	  for(j = jlo; j < jhi; j += stride)
	  {
		j1 = j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	#ifdef USE_SSE2

		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		add4 = add0+p4;
		add5 = add0+p5;
		add6 = add0+p6;
		add7 = add0+p7;

		SSE2_RADIX8_DIT_TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7,isrt2,c1,c2,c3,c4,c5,c6,c7);
	  #if 0
	  if(c1->d0 != 1.0) {	// Move debug to first case w/nontrivial twiddles
		int idbg = 0;	vec_dbl*tmp;
	//	printf("j1 = %u: SSE2 DIT Intermediates:\n",j1);
		printf("j1 = %u, c1 = %15.10f: SSE2 DIT Outputs:\n",j1,c1->d0);
		tmp = add0       ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg++;
		tmp = add0+p1    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg++;
		tmp = add0+p2    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg++;
		tmp = add0+p3    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg++;
		tmp = add0+p4    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg++;
		tmp = add0+p5    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg++;
		tmp = add0+p6    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg++;
		tmp = add0+p7    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		exit(0);
	  }
	  #endif

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
	  #if 0
		int idbg = 0;
		printf("j1 = %u: SSE2 DIT Intermediates:\n",j1);
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t1,t2);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t3,t4);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t5,t6);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t7,t8);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t9,t10);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t11,t12);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t13,t14);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t15,t16);
		exit(0);
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
	  #if 0
	  if(c1->d0 != 1.0) {	// Move debug to first case w/nontrivial twiddles
		int idbg = 0;
		printf("j1 = %u, c1 = %15.10f: Scalar-double DIT Outputs:\n",j1);
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[j1   ],a[j2   ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[j1+p1],a[j2+p1]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[j1+p2],a[j2+p2]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[j1+p3],a[j2+p3]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[j1+p4],a[j2+p4]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[j1+p5],a[j2+p5]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[j1+p6],a[j2+p6]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[j1+p7],a[j2+p7]);
		exit(0);
	  }
	  #endif

	#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

