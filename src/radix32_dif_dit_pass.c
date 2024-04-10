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

	#include "radix32_dif_dit_pass_asm.h"

#endif

/***************/

/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform a radix-32 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details.
*/
void radix32_dif_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
{
	static int max_threads = 0;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int i,j,j1,j2,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	int p01,p02,p03,p04,p08,p0C,p10,p14,p18,p1C;
	const double c = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/
	double rt,it,re0,im0,re1,im1;

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error SSE2 code not supported for this compiler!
  #endif

	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0;	/* Addresses into array sections */
	vec_dbl *c_tmp,*s_tmp;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;					// Base address for discrete per-thread local stores
	// In || mode, only above base-pointer (shared by all threads) is static:
	vec_dbl *isrt2,*sqrt2, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3, *one,*two, *r00;
  #else
	static vec_dbl *isrt2,*sqrt2, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3, *one,*two, *r00;
  #endif

#else

	int jp,jt;
	double *addr, *addp;
	int prefetch_offset;
	double	 c01,c02,c03,c04,c05,c06,c07,c08,c09,c0A,c0B,c0C,c0D,c0E,c0F
		,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c1A,c1B,c1C,c1D,c1E,c1F
		    ,s01,s02,s03,s04,s05,s06,s07,s08,s09,s0A,s0B,s0C,s0D,s0E,s0F
		,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s1A,s1B,s1C,s1D,s1E,s1F
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F;

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
		sc_arr = ALLOC_VEC_DBL(sc_arr, 0x90*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 64 16-byte slots of sc_arr for temporaries, next 7 for the nontrivial complex 32nd roots,
	last 64 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
		#ifdef MULTITHREAD
	//	if(max_threads > 1) {
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x40;
			cc0	  = sc_ptr + 0x41;
			ss0	  = sc_ptr + 0x42;
			cc1	  = sc_ptr + 0x43;
			ss1	  = sc_ptr + 0x44;
			cc3	  = sc_ptr + 0x45;
			ss3	  = sc_ptr + 0x46;
			one   = sc_ptr + 0x87;
			two   = sc_ptr + 0x88;
			sqrt2 = sc_ptr + 0x89;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
				VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );
				VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
				VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
				VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
				/* Move on to next thread's local store */
				isrt2 += 0x90;
				cc0   += 0x90;
				ss0   += 0x90;
				cc1   += 0x90;
				ss1   += 0x90;
				cc3   += 0x90;
				ss3   += 0x90;
				one   += 0x90;
				two   += 0x90;
				sqrt2 += 0x90;
			}
		#else
			r00   = sc_ptr;
			isrt2 = sc_ptr + 0x40;
			cc0	  = sc_ptr + 0x41;
			ss0	  = sc_ptr + 0x42;
			cc1	  = sc_ptr + 0x43;
			ss1	  = sc_ptr + 0x44;
			cc3	  = sc_ptr + 0x45;
			ss3	  = sc_ptr + 0x46;
			one   = sc_ptr + 0x87;
			two   = sc_ptr + 0x88;
			sqrt2 = sc_ptr + 0x89;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
			VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );
			VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
			VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
			VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
		#endif
	//	}
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		r00 = __r0 + thr_id*0x90;
		cc0	= r00 + 0x41;
	#endif

#endif

	p01 = incr >> 5;
	p02 = p01 +p01;
	p03 = p02 +p01;
	p04 = p03 +p01;
	p08 = p04 +p04;
	p0C = p08 +p04;
	p10 = p0C +p04;
	p14 = p10 +p04;
	p18 = p14 +p04;
	p1C = p18 +p04;

	p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
	p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
	p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
	p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
	p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
	p0C = p0C + ( (p0C >> DAT_BITS) << PAD_BITS );
	p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
	p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
	p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
	p1C = p1C + ( (p1C >> DAT_BITS) << PAD_BITS );
	ASSERT(p04+p04 == p08, "p04+p04 != p08");
	ASSERT(p04+p08 == p0C, "p04+p08 != p0C");

/*...The radix-32 pass is here.	*/

	iroot_prim=(incr >> 6);		/* (incr/2)/radix_now		*/

	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16).x	*/
	{				/* NLOOPS satisfies the identity NLOOPS * INCR = P*N, which says that we process the entire
					   array each time this subroutine is executed (since P*N = vector length, sans paddring.)	*/

/*	here are the needed sincos data - these are processed below in bit-reversed order.	*/
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* In the DIF pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks. */
#if HIACC
	#ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-31] are offset w.r.to the thread-local ptr pair as
					cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
		(cc0,ss0) + 0x[06,26,16,36|0e,2e,1e,3e|0a,2a,1a,3a|12,32,22,42|08,28,18,38|10,30,20,40|0c,2c,1c,3c|14,34,24,44].
		*/
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c01=rt;		s01=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c02=rt;		s02=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c03=rt;		s03=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c04=rt;		s04=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c05=rt;		s05=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c06=rt;		s06=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c07=rt;		s07=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c08=rt;		s08=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c09=rt;		s09=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0A=rt;		s0A=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0B=rt;		s0B=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x12; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0C=rt;		s0C=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x32; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0D=rt;		s0D=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0E=rt;		s0E=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x42; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0F=rt;		s0F=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x08; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c10=rt;		s10=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x28; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c11=rt;		s11=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x18; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c12=rt;		s12=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x38; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c13=rt;		s13=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
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
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c15=rt;		s15=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x20; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c16=rt;		s16=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x40; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c17=rt;		s17=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c18=rt;		s18=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c19=rt;		s19=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1A=rt;		s1A=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1B=rt;		s1B=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1C=rt;		s1C=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x34; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1D=rt;		s1D=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x24; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1E=rt;		s1E=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x44; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1F=rt;		s1F=it;
	#endif

#else
	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;			/* 2*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c01=t00*rt -t01*it;	s01=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 3*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c02=t00*rt -t01*it;	s02=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += 4*iroot;				/* 7*iroot	*/
	    iroot = i;				/* 7*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c03=t00*rt -t01*it;	s03=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 14*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c07=t00*rt -t01*it;	s07=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 21*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c0E=t00*rt -t01*it;	s0E=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 28*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c15=t00*rt -t01*it;	s15=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c1C=t00*rt -t01*it;	s1C=t00*it +t01*rt;

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

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 5);

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

#ifdef USE_SSE2

	/* Gather needed data (32 64-bit complex, i.e. 64 64-bit reals) and do first set of four length-8 transforms,
		processing sincos data in bit-reversed order.	*/

		add0 = &a[j1];
		SSE2_RADIX32_DIF_TWIDDLE(add0,p01,p02,p03,p04,p08,p0C,p10,p18,r00)

#else	// !USE_SEE2: Scalar (Non-SIMD) code:

	/*       gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...
			 We process the sincos data in bit-reversed order.	*/
	#ifdef PFETCH_AGGRESSIVE
		addr = &a[j1];
	#elif PFETCH
		prefetch_offset = ((j >> 1) & 7)*p04 + 4;	/* Cycle among p00, p04, p08, p0C, p10, p14, p18 and p1C. */
	#endif
	/*...Block 1: */

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t00=a[jt    ];									t01=a[jp    ];
			rt =a[jt+p10]*c10 - a[jp+p10]*s10;	it =a[jp+p10]*c10 + a[jt+p10]*s10;
			t02=t00-rt;		t03=t01-it;
			t00=t00+rt;		t01=t01+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t04=a[jt+p08]*c08 - a[jp+p08]*s08;	t05=a[jp+p08]*c08 + a[jt+p08]*s08;
			rt =a[jt+p18]*c18 - a[jp+p18]*s18;	it =a[jp+p18]*c18 + a[jt+p18]*s18;
			t06=t04-rt;		t07=t05-it;
			t04=t04+rt;		t05=t05+it;

			rt =t04;		it =t05;
			t04=t00-rt;		t05=t01-it;
			t00=t00+rt;		t01=t01+it;

			rt =t06;		it =t07;
			t06=t02+it;		t07=t03-rt;
			t02=t02-it;		t03=t03+rt;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t08=a[jt+p04]*c04 - a[jp+p04]*s04;	t09=a[jp+p04]*c04 + a[jt+p04]*s04;
			rt =a[jt+p14]*c14 - a[jp+p14]*s14;	it =a[jp+p14]*c14 + a[jt+p14]*s14;
			t0A=t08-rt;		t0B=t09-it;
			t08=t08+rt;		t09=t09+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t0C=a[jt+p0C]*c0C - a[jp+p0C]*s0C;	t0D=a[jp+p0C]*c0C + a[jt+p0C]*s0C;
			rt =a[jt+p1C]*c1C - a[jp+p1C]*s1C;	it =a[jp+p1C]*c1C + a[jt+p1C]*s1C;
			t0E=t0C-rt;		t0F=t0D-it;
			t0C=t0C+rt;		t0D=t0D+it;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p04;
		prefetch_p_doubles(addp);
	#endif
			rt =t0C;		it =t0D;
			t0C=t08-rt;		t0D=t09-it;
			t08=t08+rt;		t09=t09+it;

			rt =t0E;		it =t0F;
			t0E=t0A+it;		t0F=t0B-rt;
			t0A=t0A-it;		t0B=t0B+rt;

			rt =t08;		it =t09;
			t08=t00-rt;		t09=t01-it;
			t00=t00+rt;		t01=t01+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
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

		jt = j1 + p02;
		jp = j2 + p02;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t10=a[jt    ]*c02 - a[jp    ]*s02;	t11=a[jp    ]*c02 + a[jt    ]*s02;
			rt =a[jt+p10]*c12 - a[jp+p10]*s12;	it =a[jp+p10]*c12 + a[jt+p10]*s12;
			t12=t10-rt;		t13=t11-it;
			t10=t10+rt;		t11=t11+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t14=a[jt+p08]*c0A - a[jp+p08]*s0A;	t15=a[jp+p08]*c0A + a[jt+p08]*s0A;
			rt =a[jt+p18]*c1A - a[jp+p18]*s1A;	it =a[jp+p18]*c1A + a[jt+p18]*s1A;
			t16=t14-rt;		t17=t15-it;
			t14=t14+rt;		t15=t15+it;

			rt =t14;		it =t15;
			t14=t10-rt;		t15=t11-it;
			t10=t10+rt;		t11=t11+it;

			rt =t16;		it =t17;
			t16=t12+it;		t17=t13-rt;
			t12=t12-it;		t13=t13+rt;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p08;
		prefetch_p_doubles(addp);
	#endif
			t18=a[jt+p04]*c06 - a[jp+p04]*s06;	t19=a[jp+p04]*c06 + a[jt+p04]*s06;
			rt =a[jt+p14]*c16 - a[jp+p14]*s16;	it =a[jp+p14]*c16 + a[jt+p14]*s16;
			t1A=t18-rt;		t1B=t19-it;
			t18=t18+rt;		t19=t19+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t1C=a[jt+p0C]*c0E - a[jp+p0C]*s0E;	t1D=a[jp+p0C]*c0E + a[jt+p0C]*s0E;
			rt =a[jt+p1C]*c1E - a[jp+p1C]*s1E;	it =a[jp+p1C]*c1E + a[jt+p1C]*s1E;
			t1E=t1C-rt;		t1F=t1D-it;
			t1C=t1C+rt;		t1D=t1D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t1C;		it =t1D;
			t1C=t18-rt;		t1D=t19-it;
			t18=t18+rt;		t19=t19+it;

			rt =t1E;		it =t1F;
			t1E=t1A+it;		t1F=t1B-rt;
			t1A=t1A-it;		t1B=t1B+rt;

			rt =t18;		it =t19;
			t18=t10-rt;		t19=t11-it;
			t10=t10+rt;		t11=t11+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
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

		jt = j1 + p01;
		jp = j2 + p01;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p0C;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t20=a[jt    ]*c01 - a[jp    ]*s01;	t21=a[jp    ]*c01 + a[jt    ]*s01;
			rt =a[jt+p10]*c11 - a[jp+p10]*s11;	it =a[jp+p10]*c11 + a[jt+p10]*s11;
			t22=t20-rt;		t23=t21-it;
			t20=t20+rt;		t21=t21+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t24=a[jt+p08]*c09 - a[jp+p08]*s09;	t25=a[jp+p08]*c09 + a[jt+p08]*s09;
			rt =a[jt+p18]*c19 - a[jp+p18]*s19;	it =a[jp+p18]*c19 + a[jt+p18]*s19;
			t26=t24-rt;		t27=t25-it;
			t24=t24+rt;		t25=t25+it;

			rt =t24;		it =t25;
			t24=t20-rt;		t25=t21-it;
			t20=t20+rt;		t21=t21+it;

			rt =t26;		it =t27;
			t26=t22+it;		t27=t23-rt;
			t22=t22-it;		t23=t23+rt;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t28=a[jt+p04]*c05 - a[jp+p04]*s05;	t29=a[jp+p04]*c05 + a[jt+p04]*s05;
			rt =a[jt+p14]*c15 - a[jp+p14]*s15;	it =a[jp+p14]*c15 + a[jt+p14]*s15;
			t2A=t28-rt;		t2B=t29-it;
			t28=t28+rt;		t29=t29+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t2C=a[jt+p0C]*c0D - a[jp+p0C]*s0D;	t2D=a[jp+p0C]*c0D + a[jt+p0C]*s0D;
			rt =a[jt+p1C]*c1D - a[jp+p1C]*s1D;	it =a[jp+p1C]*c1D + a[jt+p1C]*s1D;
			t2E=t2C-rt;		t2F=t2D-it;
			t2C=t2C+rt;		t2D=t2D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p10;
		prefetch_p_doubles(addp);
	#endif

			rt =t2C;		it =t2D;
			t2C=t28-rt;		t2D=t29-it;
			t28=t28+rt;		t29=t29+it;

			rt =t2E;		it =t2F;
			t2E=t2A+it;		t2F=t2B-rt;
			t2A=t2A-it;		t2B=t2B+rt;

			rt =t28;		it =t29;
			t28=t20-rt;		t29=t21-it;
			t20=t20+rt;		t21=t21+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
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

		jt = j1 + p03;
		jp = j2 + p03;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t30=a[jt    ]*c03 - a[jp    ]*s03;	t31=a[jp    ]*c03 + a[jt    ]*s03;
			rt =a[jt+p10]*c13 - a[jp+p10]*s13;	it =a[jp+p10]*c13 + a[jt+p10]*s13;
			t32=t30-rt;		t33=t31-it;
			t30=t30+rt;		t31=t31+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t34=a[jt+p08]*c0B - a[jp+p08]*s0B;	t35=a[jp+p08]*c0B + a[jt+p08]*s0B;
			rt =a[jt+p18]*c1B - a[jp+p18]*s1B;	it =a[jp+p18]*c1B + a[jt+p18]*s1B;
			t36=t34-rt;		t37=t35-it;
			t34=t34+rt;		t35=t35+it;

			rt =t34;		it =t35;
			t34=t30-rt;		t35=t31-it;
			t30=t30+rt;		t31=t31+it;

			rt =t36;		it =t37;
			t36=t32+it;		t37=t33-rt;
			t32=t32-it;		t33=t33+rt;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p14;
		prefetch_p_doubles(addp);
	#endif
			t38=a[jt+p04]*c07 - a[jp+p04]*s07;	t39=a[jp+p04]*c07 + a[jt+p04]*s07;
			rt =a[jt+p14]*c17 - a[jp+p14]*s17;	it =a[jp+p14]*c17 + a[jt+p14]*s17;
			t3A=t38-rt;		t3B=t39-it;
			t38=t38+rt;		t39=t39+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t3C=a[jt+p0C]*c0F - a[jp+p0C]*s0F;	t3D=a[jp+p0C]*c0F + a[jt+p0C]*s0F;
			rt =a[jt+p1C]*c1F - a[jp+p1C]*s1F;	it =a[jp+p1C]*c1F + a[jt+p1C]*s1F;
			t3E=t3C-rt;		t3F=t3D-it;
			t3C=t3C+rt;		t3D=t3D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t3C;		it =t3D;
			t3C=t38-rt;		t3D=t39-it;
			t38=t38+rt;		t39=t39+it;

			rt =t3E;		it =t3F;
			t3E=t3A+it;		t3F=t3B-rt;
			t3A=t3A-it;		t3B=t3B+rt;

			rt =t38;		it =t39;
			t38=t30-rt;		t39=t31-it;
			t30=t30+rt;		t31=t31+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t3C;		it =t3D;
			t3C=t34+it;		t3D=t35-rt;
			t34=t34-it;		t35=t35+rt;

			rt =(t3A-t3B)*ISRT2;it =(t3A+t3B)*ISRT2;
			t3A=t32-rt;		t3B=t33-it;
			t32=t32+rt;		t33=t33+it;

			rt =(t3E+t3F)*ISRT2;it =(t3F-t3E)*ISRT2;
			t3E=t36+rt;		t3F=t37+it;
			t36=t36-rt;		t37=t37-it;
	  #if 0
	  if(c01 != 1.0) {	// Move debug to first case w/nontrivial twiddles
		int idbg = 0;
		printf("j1 = %u, c1 = %15.10f: Scalar-double DIF Intermediates:\n",j1,c01);
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t00,t01);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t02,t03);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t04,t05);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t06,t07);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t08,t09);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t0A,t0B);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t0C,t0D);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t0E,t0F);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t10,t11);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t12,t13);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t14,t15);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t16,t17);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t18,t19);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t1A,t1B);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t1C,t1D);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t1E,t1F);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t20,t21);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t22,t23);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t24,t25);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t26,t27);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t28,t29);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t2A,t2B);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t2C,t2D);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t2E,t2F);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t30,t31);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t32,t33);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t34,t35);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t36,t37);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t38,t39);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t3A,t3B);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t3C,t3D);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,t3E,t3F);	idbg++;
		exit(0);
	  }
	  #endif

	/*...and now do eight radix-4 transforms, including the internal twiddle factors:
		1, exp(i* 1*twopi/32) =       ( c32_1, s32_1), exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 3*twopi/32) =       ( c32_3, s32_3) (for inputs to transform block 2),
		1, exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 6*twopi/32) =       ( s    , c    ) (for inputs to transform block 3),
		1, exp(i* 3*twopi/32) =       ( c32_3, s32_3), exp(i* 6*twopi/32) =       ( s    , c    ), exp(i* 9*twopi/32) =       (-s32_1, c32_1) (for inputs to transform block 4),
		1, exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 8*twopi/32) =       ( 0    , 1    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ) (for inputs to transform block 5),
		1, exp(i* 5*twopi/32) =       ( s32_3, c32_3), exp(i*10*twopi/32) =       (-s    , c    ), exp(i*15*twopi/32) =       (-c32_1, s32_1) (for inputs to transform block 6),
		1, exp(i* 6*twopi/32) =       ( s    , c    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ), exp(i*18*twopi/32) =       (-c    ,-s    ) (for inputs to transform block 7),
		1, exp(i* 7*twopi/32) =       ( s32_1, c32_1), exp(i*14*twopi/32) =       (-c    , s    ), exp(i*21*twopi/32) =       (-s32_3,-c32_3) (for inputs to transform block 8),
		 and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.	*/

	/*...Block 1: t00,t10,t20,t30	*/

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p18;
		prefetch_p_doubles(addp);
	#endif
			rt =t10;	t10=t00-rt;	t00=t00+rt;
			it =t11;	t11=t01-it;	t01=t01+it;

			rt =t30;	t30=t20-rt;	t20=t20+rt;
			it =t31;	t31=t21-it;	t21=t21+it;

			a[jt    ]=t00+t20;		a[jp    ]=t01+t21;
			a[jt+p01]=t00-t20;		a[jp+p01]=t01-t21;

			a[jt+p02]=t10-t31;		a[jp+p02]=t11+t30;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t10+t31;		a[jp+p03]=t11-t30;

	/*...Block 5: t08,t18,t28,t38	*/

		jt = j1 + p04;
		jp = j2 + p04;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t18;	t18=t08+t19;	t08=t08-t19;		/* twiddle mpy by E^4 = I	*/
				t19=t09-rt;	t09=t09+rt;

			rt =(t28-t29)*ISRT2;	t29=(t28+t29)*ISRT2;		t28=rt;	/* twiddle mpy by E^8	*/
			rt =(t39+t38)*ISRT2;	it =(t39-t38)*ISRT2;			/* twiddle mpy by -E^12 is here...	*/
			t38=t28+rt;			t28=t28-rt;				/* ...and get E^12 by flipping signs here.	*/
			t39=t29+it;			t29=t29-it;

			a[jt    ]=t08+t28;		a[jp    ]=t09+t29;
			a[jt+p01]=t08-t28;		a[jp+p01]=t09-t29;

			a[jt+p02]=t18-t39;		a[jp+p02]=t19+t38;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t18+t39;		a[jp+p03]=t19-t38;

	/*...Block 3: t04,t14,t24,t34	*/

		jt = j1 + p08;
		jp = j2 + p08;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =(t14-t15)*ISRT2;	it =(t14+t15)*ISRT2;			/* twiddle mpy by E^4	*/
			t14=t04-rt;			t04=t04+rt;
			t15=t05-it;			t05=t05+it;

			rt =t24*c - t25*s;		t25=t25*c + t24*s;		t24=rt;	/* twiddle mpy by E^2	*/
			rt =t34*s - t35*c;		it =t35*s + t34*c;			/* twiddle mpy by E^6	*/
			t34=t24-rt;			t24=t24+rt;
			t35=t25-it;			t25=t25+it;

			a[jt    ]=t04+t24;		a[jp    ]=t05+t25;
			a[jt+p01]=t04-t24;		a[jp+p01]=t05-t25;

			a[jt+p02]=t14-t35;		a[jp+p02]=t15+t34;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t14+t35;		a[jp+p03]=t15-t34;

	/*...Block 7: t0C,t1C,t2C,t3C	*/

		jt = j1 + p0C;
		jp = j2 + p0C;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =(t1D+t1C)*ISRT2;	it =(t1D-t1C)*ISRT2;			/* twiddle mpy by -E^12 is here...	*/
			t1C=t0C+rt;			t0C=t0C-rt;				/* ...and get E^12 by flipping signs here.	*/
			t1D=t0D+it;			t0D=t0D-it;

			rt =t2C*s - t2D*c;		t2D=t2D*s + t2C*c;		t2C=rt;	/* twiddle mpy by E^6	*/
			rt =t3C*c - t3D*s;		it =t3D*c + t3C*s;			/* twiddle mpy by E^18 is here...	*/
			t3C=t2C+rt;			t2C=t2C-rt;				/* ...and get E^18 by flipping signs here.	*/
			t3D=t2D+it;			t2D=t2D-it;

			a[jt    ]=t0C+t2C;		a[jp    ]=t0D+t2D;
			a[jt+p01]=t0C-t2C;		a[jp+p01]=t0D-t2D;

			a[jt+p02]=t1C-t3D;		a[jp+p02]=t1D+t3C;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t1C+t3D;		a[jp+p03]=t1D-t3C;

	/*...Block 2: t02,t12,t22,t32	*/

		jt = j1 + p10;
		jp = j2 + p10;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p1C;
		prefetch_p_doubles(addp);
	#endif
			rt =t12*c - t13*s;		it =t13*c + t12*s;			/* twiddle mpy by E^2	*/
			t12=t02-rt;			t02=t02+rt;
			t13=t03-it;			t03=t03+it;

			rt =t22*c32_1 - t23*s32_1;	t23=t23*c32_1 + t22*s32_1;	t22=rt;	/* twiddle mpy by E^1	*/
			rt =t32*c32_3 - t33*s32_3;	it =t33*c32_3 + t32*s32_3;		/* twiddle mpy by E^3	*/
			t32=t22-rt;			t22=t22+rt;
			t33=t23-it;			t23=t23+it;

			a[jt    ]=t02+t22;		a[jp    ]=t03+t23;
			a[jt+p01]=t02-t22;		a[jp+p01]=t03-t23;

			a[jt+p02]=t12-t33;		a[jp+p02]=t13+t32;	/* mpy by E^4=i is inlined here...;	*/
			a[jt+p03]=t12+t33;		a[jp+p03]=t13-t32;

	/*...Block 6: t0A,t1A,t2A,t3A	*/

		jt = j1 + p14;
		jp = j2 + p14;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1A*s + t1B*c;		it =t1B*s - t1A*c;			/* twiddle mpy by -E^10 is here...	*/
			t1A=t0A+rt;			t0A =t0A-rt;				/* ...and get E^10 by flipping signs here.	*/
			t1B=t0B+it;			t0B =t0B-it;

			rt =t2A*s32_3 - t2B*c32_3;	t2B=t2B*s32_3 + t2A*c32_3;	t2A=rt;	/* twiddle mpy by E^5	*/
			rt =t3A*c32_1 + t3B*s32_1;	it =t3B*c32_1 - t3A*s32_1;	 	/* twiddle mpy by -E^15 is here...	*/
			t3A=t2A+rt;			t2A=t2A-rt;				/* ...and get E^15 by flipping signs here.	*/
			t3B=t2B+it;			t2B=t2B-it;

			a[jt    ]=t0A+t2A;		a[jp    ]=t0B+t2B;
			a[jt+p01]=t0A-t2A;		a[jp+p01]=t0B-t2B;

			a[jt+p02]=t1A-t3B;		a[jp+p02]=t1B+t3A;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t1A+t3B;		a[jp+p03]=t1B-t3A;

	/*...Block 4: t06,t16,t26,t36	*/

		jt = j1 + p18;
		jp = j2 + p18;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t16*s - t17*c;		it =t17*s + t16*c;			/* twiddle mpy by E^6	*/
			t16=t06-rt;			t06 =t06+rt;
			t17=t07-it;			t07 =t07+it;

			rt =t26*c32_3 - t27*s32_3;	t27=t27*c32_3 + t26*s32_3;	t26=rt;	/* twiddle mpy by E^3	*/
			rt =t36*s32_1 + t37*c32_1;	it =t37*s32_1 - t36*c32_1;		/* twiddle mpy by -E^9 is here...	*/
			t36=t26+rt;			t26=t26-rt;				/* ...and get E^9 by flipping signs here.	*/
			t37=t27+it;			t27=t27-it;

			a[jt    ]=t06+t26;		a[jp    ]=t07+t27;
			a[jt+p01]=t06-t26;		a[jp+p01]=t07-t27;

			a[jt+p02]=t16-t37;		a[jp+p02]=t17+t36;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t16+t37;		a[jp+p03]=t17-t36;

	/*...Block 8: t0E,t1E,t2E,t3E	*/

		jt = j1 + p1C;
		jp = j2 + p1C;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1E*c + t1F*s;		it =t1F*c - t1E*s;			/* twiddle mpy by -E^14 is here...	*/
			t1E=t0E+rt;			t0E =t0E-rt;				/* ...and get E^14 by flipping signs here.	*/
			t1F=t0F+it;			t0F =t0F-it;

			rt =t2E*s32_1 - t2F*c32_1;	t2F=t2F*s32_1 + t2E*c32_1;	t2E=rt;	/* twiddle mpy by E^7	*/
			rt =t3E*s32_3 - t3F*c32_3;	it =t3F*s32_3 + t3E*c32_3;		/* twiddle mpy by -E^21 is here...	*/
			t3E=t2E+rt;			t2E=t2E-rt;				/* ...and get E^21 by flipping signs here.	*/
			t3F=t2F+it;			t2F=t2F-it;

			a[jt    ]=t0E+t2E;		a[jp    ]=t0F+t2F;
			a[jt+p01]=t0E-t2E;		a[jp+p01]=t0F-t2F;

			a[jt+p02]=t1E-t3F;		a[jp+p02]=t1F+t3E;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t1E+t3F;		a[jp+p03]=t1F-t3E;
	  #if 0
	  if(c01 != 1.0) {	// Move debug to first case w/nontrivial twiddles
		int idbg = 0;
		printf("j1 = %u, c1 = %15.10f: Scalar-double DIF Outputs:\n",j1,c01);
		jt = j1      ;	jp = j2      ;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt    ],a[jp    ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p01],a[jp+p01]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p02],a[jp+p02]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p03],a[jp+p03]);	idbg++;
		jt = j1 + p04;	jp = j2 + p04;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt    ],a[jp    ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p01],a[jp+p01]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p02],a[jp+p02]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p03],a[jp+p03]);	idbg++;
		jt = j1 + p08;	jp = j2 + p08;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt    ],a[jp    ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p01],a[jp+p01]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p02],a[jp+p02]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p03],a[jp+p03]);	idbg++;
		jt = j1 + p0C;	jp = j2 + p0C;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt    ],a[jp    ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p01],a[jp+p01]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p02],a[jp+p02]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p03],a[jp+p03]);	idbg++;
		jt = j1 + p10;	jp = j2 + p10;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt    ],a[jp    ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p01],a[jp+p01]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p02],a[jp+p02]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p03],a[jp+p03]);	idbg++;
		jt = j1 + p14;	jp = j2 + p14;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt    ],a[jp    ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p01],a[jp+p01]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p02],a[jp+p02]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p03],a[jp+p03]);	idbg++;
		jt = j1 + p18;	jp = j2 + p18;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt    ],a[jp    ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p01],a[jp+p01]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p02],a[jp+p02]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p03],a[jp+p03]);	idbg++;
		jt = j1 + p1C;	jp = j2 + p1C;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt    ],a[jp    ]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p01],a[jp+p01]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p02],a[jp+p02]);	idbg++;
		printf("\t{re,im}[%2u] = %20.10e,%20.10e\n",idbg,a[jt+p03],a[jp+p03]);	idbg++;
		exit(0);
	  }
	  #endif

#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; ...) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

/***************/

/*
!...Post-twiddles implementation of radix32_dit_pass.
*/
void radix32_dit_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
{
	static int max_threads = 0;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int i,j,j1,j2,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	int p01,p02,p03,p04,p05,p06,p07,p08,p10,p18;
	const double c = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/
	double rt,it,re0,im0,re1,im1;

#ifdef USE_SSE2

  #ifndef COMPILER_TYPE_GCC
	#error SSE2 code not supported for this compiler!
  #endif

	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	double *add0;	/* Addresses into array sections */
	vec_dbl *c_tmp,*s_tmp;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	/* Base address for discrete per-thread local stores */
	// In || mode, only above base-pointer (shared by all threads) is static:
	vec_dbl *isrt2, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3, *one,*two,*sqrt2, *r00;
  #else
	static vec_dbl *isrt2, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3, *one,*two,*sqrt2, *r00;
  #endif

#else

	int jp,jt;
	double *addr, *addp;
	int prefetch_offset;
	double	 c01,c02,c03,c04,c05,c06,c07,c08,c09,c0A,c0B,c0C,c0D,c0E,c0F
		,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c1A,c1B,c1C,c1D,c1E,c1F
		    ,s01,s02,s03,s04,s05,s06,s07,s08,s09,s0A,s0B,s0C,s0D,s0E,s0F
		,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s1A,s1B,s1C,s1D,s1E,s1F
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F;

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
		sc_arr = ALLOC_VEC_DBL(sc_arr, 0x90*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 64 16-byte slots of sc_arr for temporaries, next 7 for the nontrivial complex 32nd roots,
	last 64 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
	#ifdef MULTITHREAD
		__r0  = sc_ptr;
		isrt2 = sc_ptr + 0x40;
		cc0	  = sc_ptr + 0x41;
		ss0	  = sc_ptr + 0x42;
		cc1	  = sc_ptr + 0x43;
		ss1	  = sc_ptr + 0x44;
		cc3	  = sc_ptr + 0x45;
		ss3	  = sc_ptr + 0x46;
		one   = sc_ptr + 0x87;
		two   = sc_ptr + 0x88;
		sqrt2 = sc_ptr + 0x89;
		for(i = 0; i < max_threads; ++i) {
			/* These remain fixed within each per-thread local store: */
			VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
			VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );
			VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
			VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
			VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
			/* Move on to next thread's local store */
			isrt2 += 0x90;
			cc0   += 0x90;
			ss0   += 0x90;
			cc1   += 0x90;
			ss1   += 0x90;
			cc3   += 0x90;
			ss3   += 0x90;
			one   += 0x90;
			two   += 0x90;
			sqrt2 += 0x90;
		}
		#else
			r00   = sc_ptr;
			isrt2 = sc_ptr + 0x40;
			cc0	  = sc_ptr + 0x41;
			ss0	  = sc_ptr + 0x42;
			cc1	  = sc_ptr + 0x43;
			ss1	  = sc_ptr + 0x44;
			cc3	  = sc_ptr + 0x45;
			ss3	  = sc_ptr + 0x46;
			one   = sc_ptr + 0x87;
			two   = sc_ptr + 0x88;
			sqrt2 = sc_ptr + 0x89;
			/* These remain fixed: */
			VEC_DBL_INIT(isrt2, ISRT2);		VEC_DBL_INIT(sqrt2, SQRT2);
			VEC_DBL_INIT(one  , 1.0  );		VEC_DBL_INIT(two, 2.0  );
			VEC_DBL_INIT(cc0  , c    );		VEC_DBL_INIT(ss0, s    );
			VEC_DBL_INIT(cc1  , c32_1);		VEC_DBL_INIT(ss1, s32_1);
			VEC_DBL_INIT(cc3  , c32_3);		VEC_DBL_INIT(ss3, s32_3);
	#endif
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
  #ifdef MULTITHREAD
	ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
	r00 = __r0 + thr_id*0x90;
	isrt2 = r00 + 0x40;
	cc0	= isrt2 + 1;
  #endif

#endif	// USE_SSE2 ?

	p01 = incr >> 5;
	p02 = p01 +p01;
	p03 = p02 +p01;
	p04 = p03 +p01;
	p05 = p04 +p01;
	p06 = p05 +p01;
	p07 = p06 +p01;
	p08 = p07 +p01;
	p10 = p08 +p08;
	p18 = p10 +p08;

	p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
	p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
	p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
	p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
	p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
	p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
	p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
	p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
	p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
	p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );

/*...The radix-32 pass is here.	*/

	iroot_prim=(incr >> 6);		/* (incr/2)/radix_now */

	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16). NLOOPS satisfies the */
	{				/* identity NLOOPS * INCR = P*N, which says that we process the entire array each time this subroutine is executed (since P*N = vector length, sans padding.) */
	/*	here are the needed sincos data - these are processed below in bit-reversed order.	*/
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* In the DIT pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks. */
#if HIACC
	#ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-31] are offset w.r.to the thread-local ptr pair as
					cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F] are in
		(cc0,ss0) + 0x[06,26,16,36|0e,2e,1e,3e|08,28,18,38|10,30,20,40|0a,2a,1a,3a|12,32,22,42|0c,2c,1c,3c|14,34,24,44].
	*** NOTE: This is the same pattern as for DIF, but with the middle 2 roots octets [08-40] and [0a-42] swapped ***
		Or, in terms of the "reverse directory" lookup:
		(cc0,ss0) + 0x[06,08,0a,0c,0e,10,12,14,16,18,1a,1c,1e,20,22,24,26,28,2a,2c,2e,30,32,34,36,38,3a,3c,3e,40,42,44] hold
					cc[00 08 10 18 04 0C 14 1C 02 0A 12 1A 06 0E 16 1E 01 09 11 19 05 0D 15 1D 03 0B 13 1B 07 0F 17 1F]
		*/
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;	/* c0,s0 */
		rt = 1.0; it = 0.0;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c01=rt;		s01=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c02=rt;		s02=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c03=rt;		s03=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c04=rt;		s04=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c05=rt;		s05=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c06=rt;		s06=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c07=rt;		s07=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x08; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c08=rt;		s08=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x28; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c09=rt;		s09=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x18; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0A=rt;		s0A=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x38; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0B=rt;		s0B=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x10; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0C=rt;		s0C=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0D=rt;		s0D=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x20; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0E=rt;		s0E=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x40; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c0F=rt;		s0F=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c10=rt;		s10=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c11=rt;		s11=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c12=rt;		s12=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3a; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c13=rt;		s13=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x12; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c14=rt;		s14=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x32; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c15=rt;		s15=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c16=rt;		s16=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x42; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c17=rt;		s17=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c18=rt;		s18=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c19=rt;		s19=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1A=rt;		s1A=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3c; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1B=rt;		s1B=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1C=rt;		s1C=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x34; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1D=rt;		s1D=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x24; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1E=rt;		s1E=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x44; s_tmp = c_tmp+1;
		VEC_DBL_INIT(c_tmp, rt);
		VEC_DBL_INIT(s_tmp, it);
	#else
		c1F=rt;		s1F=it;
	#endif

#else

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;			/* 2*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c01=t00*rt -t01*it;	s01=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 3*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c02=t00*rt -t01*it;	s02=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += 4*iroot;				/* 7*iroot	*/
	    iroot = i;				/* 7*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c03=t00*rt -t01*it;	s03=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 14*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c07=t00*rt -t01*it;	s07=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 21*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c0E=t00*rt -t01*it;	s0E=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 28*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c15=t00*rt -t01*it;	s15=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c1C=t00*rt -t01*it;	s1C=t00*it +t01*rt;

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

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 5);

/******************* AVX debug stuff: *******************/
#if 0
if(1) {
	int idbg;
	// Use RNG to populate data array:
	rng_isaac_init(TRUE);
	double dtmp = 1024.0*1024.0*1024.0*1024.0;
	int jt = jlo;
	jt += ( (jt >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
  #ifdef USE_AVX512
	idbg =   0;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br16[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
  #elif defined(USE_AVX)
	idbg =   0;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br8[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
  #elif defined(USE_SSE2)
	idbg =   0;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07;		for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p08;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p10;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p18;	for(i = 0; i < 2*RE_IM_STRIDE; i++) { a[jt+idbg+br4[ i]] = dtmp*rng_isaac_rand_double_norm_pm1(); }
  #else	// Force 4-slot inits [2*RE_IM_STRIDE replaced by 4] here to enable comparison vs SSE2 code:
	idbg =   0;		for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01;		for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02;		for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03;		for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04;		for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05;		for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06;		for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07;		for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p08;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p08;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p08;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p08;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p08;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p08;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p08;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p08;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p10;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p10;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p10;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p10;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p10;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p10;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p10;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p10;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg =   0+p18;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p01+p18;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p02+p18;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p03+p18;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p04+p18;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p05+p18;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p06+p18;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
	idbg = p07+p18;	for(i = 0; i < 4; i++) { a[jt+idbg+i] = dtmp*rng_isaac_rand_double_norm_pm1(); }
  #endif
}
#endif
/********************************************************/

	/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
	Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
	but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
	(But only ever need to explicitly do this in debug mode).
	*/
	  for(j = jlo; j < jhi; j += stride)
	  {
		j1 = j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

#ifdef USE_SSE2

	/* Gather needed data (32 64-bit complex, i.e. 64 64-bit reals) and do first set of four length-8 transforms,
		processing sincos data in bit-reversed order.	*/
		add0 = &a[j1];
		SSE2_RADIX32_DIT_TWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p10,r00,isrt2)

	  #if 0
	  if((cc0 + 0x26)->d0 != 1.0) {	// Move debug to first case w/nontrivial twiddles
		int idbg;	vec_dbl*tmp = r00;
		printf("j1 = %u, c1 = %15.10f: SSE2 32-DIT Intermediates:\n",j1,(cc0 + 0x26)->d0);
		for(idbg = 0; idbg < 32; idbg++, tmp += 2) {
			printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		}
		exit(0);
	  }
	  #elif 0
	  if((cc0 + 0x26)->d0 != 1.0) {	// Move debug to first case w/nontrivial twiddles
		int jt,jp,idbg = 0;	vec_dbl*tmp;
		printf("j1 = %u, c1 = %15.10f: SSE2 DIT Outputs:\n",j1,(cc0 + 0x26)->d0);
		jt = 0;		idbg = 0;
		tmp = add0+jt    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p08;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p10;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p18;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		jt = p04;	idbg = 4;
		tmp = add0+jt    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p08;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p10;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p18;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		jt = p01;	idbg = 1;
		tmp = add0+jt    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p08;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p10;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p18;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		jt = p05;	idbg = 5;
		tmp = add0+jt    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p08;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p10;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p18;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		jt = p02;	idbg = 2;
		tmp = add0+jt    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p08;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p10;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p18;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		jt = p06;	idbg = 6;
		tmp = add0+jt    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p08;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p10;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p18;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		jt = p03;	idbg = 3;
		tmp = add0+jt    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p08;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p10;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p18;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		jt = p07;	idbg = 7;
		tmp = add0+jt    ;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p08;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p10;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);	idbg+=8;
		tmp = add0+jt+p18;	printf("\t{re,im}[%2u] = %20.10e,%20.10e,%20.10e,%20.10e\n",idbg,tmp->d0,tmp->d1,(tmp+1)->d0,(tmp+1)->d1);
		exit(0);
	  }
	  #endif

#else	/* USE_SSE2 */

	/*       gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...
			 We process the sincos data in bit-reversed order.	*/
	#ifdef PFETCH_AGGRESSIVE
		addr = &a[j1];
	#elif PFETCH
		prefetch_offset = ((j >> 1) & 7)*p04 + 4;	/* Cycle among p00, p04, p08, p0C, p10, p14, p18 and p1C. */
	#endif

	/*...Block 1: */

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t00=a[jt    ];	t01=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t02=t00-rt;		t03=t01-it;
			t00=t00+rt;		t01=t01+it;

			t04=a[jt+p02];	t05=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t06=t04-rt;		t07=t05-it;
			t04=t04+rt;		t05=t05+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t04;		it =t05;
			t04=t00-rt;		t05=t01-it;
			t00=t00+rt;		t01=t01+it;

			rt =t06;		it =t07;
			t06=t02-it;		t07=t03+rt;
			t02=t02+it;		t03=t03-rt;

			t08=a[jt+p04];	t09=a[jp+p04];
			rt =a[jt+p05];	it =a[jp+p05];
			t0A=t08-rt;		t0B=t09-it;
			t08=t08+rt;		t09=t09+it;

			t0C=a[jt+p06];	t0D=a[jp+p06];
			rt =a[jt+p07];	it =a[jp+p07];
			t0E=t0C-rt;		t0F=t0D-it;
			t0C=t0C+rt;		t0D=t0D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t0C;		it =t0D;
			t0C=t08-rt;		t0D=t09-it;
			t08=t08+rt;		t09=t09+it;

			rt =t0E;		it =t0F;
			t0E=t0A-it;		t0F=t0B+rt;
			t0A=t0A+it;		t0B=t0B-rt;

			rt =t08;		it =t09;
			t08=t00-rt;		t09=t01-it;
			t00=t00+rt;		t01=t01+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t0C;		it =t0D;
			t0C=t04-it;		t0D=t05+rt;
			t04=t04+it;		t05=t05-rt;

			rt =(t0A+t0B)*ISRT2;it =(t0A-t0B)*ISRT2;
			t0A=t02-rt;		t0B=t03+it;
			t02=t02+rt;		t03=t03-it;

			rt =(t0E-t0F)*ISRT2;it =(t0F+t0E)*ISRT2;
			t0E=t06+rt;		t0F=t07+it;
			t06=t06-rt;		t07=t07-it;

	/*...Block 2:;	*/

		jt = j1 + p08;
		jp = j2 + p08;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t10=a[jt    ];	t11=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t12=t10-rt;		t13=t11-it;
			t10=t10+rt;		t11=t11+it;

			t14=a[jt+p02];	t15=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t16=t14-rt;		t17=t15-it;
			t14=t14+rt;		t15=t15+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t14;		it =t15;
			t14=t10-rt;		t15=t11-it;
			t10=t10+rt;		t11=t11+it;

			rt =t16;		it =t17;
			t16=t12-it;		t17=t13+rt;
			t12=t12+it;		t13=t13-rt;

			t18=a[jt+p04];	t19=a[jp+p04];
			rt =a[jt+p05];	it =a[jp+p05];
			t1A=t18-rt;		t1B=t19-it;
			t18=t18+rt;		t19=t19+it;

			t1C=a[jt+p06];	t1D=a[jp+p06];
			rt =a[jt+p07];	it =a[jp+p07];
			t1E=t1C-rt;		t1F=t1D-it;
			t1C=t1C+rt;		t1D=t1D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t1C;		it =t1D;
			t1C=t18-rt;		t1D=t19-it;
			t18=t18+rt;		t19=t19+it;

			rt =t1E;		it =t1F;
			t1E=t1A-it;		t1F=t1B+rt;
			t1A=t1A+it;		t1B=t1B-rt;

			rt =t18;		it =t19;
			t18=t10-rt;		t19=t11-it;
			t10=t10+rt;		t11=t11+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
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

		jt = j1 + p10;
		jp = j2 + p10;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t20=a[jt    ];	t21=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t22=t20-rt;		t23=t21-it;
			t20=t20+rt;		t21=t21+it;

			t24=a[jt+p02];	t25=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t26=t24-rt;		t27=t25-it;
			t24=t24+rt;		t25=t25+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t24;		it =t25;
			t24=t20-rt;		t25=t21-it;
			t20=t20+rt;		t21=t21+it;

			rt =t26;		it =t27;
			t26=t22-it;		t27=t23+rt;
			t22=t22+it;		t23=t23-rt;

			t28=a[jt+p04];	t29=a[jp+p04];
			rt =a[jt+p05];	it =a[jp+p05];
			t2A=t28-rt;		t2B=t29-it;
			t28=t28+rt;		t29=t29+it;

			t2C=a[jt+p06];	t2D=a[jp+p06];
			rt =a[jt+p07];	it =a[jp+p07];
			t2E=t2C-rt;		t2F=t2D-it;
			t2C=t2C+rt;		t2D=t2D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t2C;		it =t2D;
			t2C=t28-rt;		t2D=t29-it;
			t28=t28+rt;		t29=t29+it;

			rt =t2E;		it =t2F;
			t2E=t2A-it;		t2F=t2B+rt;
			t2A=t2A+it;		t2B=t2B-rt;

			rt =t28;		it =t29;
			t28=t20-rt;		t29=t21-it;
			t20=t20+rt;		t21=t21+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
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

		jt = j1 + p18;
		jp = j2 + p18;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t30=a[jt    ];	t31=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t32=t30-rt;		t33=t31-it;
			t30=t30+rt;		t31=t31+it;

			t34=a[jt+p02];	t35=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t36=t34-rt;		t37=t35-it;
			t34=t34+rt;		t35=t35+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t34;		it =t35;
			t34=t30-rt;		t35=t31-it;
			t30=t30+rt;		t31=t31+it;

			rt =t36;		it =t37;
			t36=t32-it;		t37=t33+rt;
			t32=t32+it;		t33=t33-rt;

			t38=a[jt+p04];	t39=a[jp+p04];
			rt =a[jt+p05];	it =a[jp+p05];
			t3A=t38-rt;		t3B=t39-it;
			t38=t38+rt;		t39=t39+it;

			t3C=a[jt+p06];	t3D=a[jp+p06];
			rt =a[jt+p07];	it =a[jp+p07];
			t3E=t3C-rt;		t3F=t3D-it;
			t3C=t3C+rt;		t3D=t3D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t3C;		it =t3D;
			t3C=t38-rt;		t3D=t39-it;
			t38=t38+rt;		t39=t39+it;

			rt =t3E;		it =t3F;
			t3E=t3A-it;		t3F=t3B+rt;
			t3A=t3A+it;		t3B=t3B-rt;

			rt =t38;		it =t39;
			t38=t30-rt;		t39=t31-it;
			t30=t30+rt;		t31=t31+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
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
		1, exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 6*twopi/32) =       ( s    ,-c    ) (for inputs to transform block 3),
		1, exp(-i* 3*twopi/32) =       ( c32_3,-s32_3), exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i* 9*twopi/32) =       (-s32_1,-c32_1) (for inputs to transform block 4),
		1, exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 8*twopi/32) =       ( 0    ,-1    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ) (for inputs to transform block 5),
		1, exp(-i* 5*twopi/32) =       ( s32_3,-c32_3), exp(-i*10*twopi/32) =       (-s    ,-c    ), exp(-i*15*twopi/32) =       (-c32_1,-s32_1) (for inputs to transform block 6),
		1, exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ), exp(-i*18*twopi/32) =       (-c    , s    ) (for inputs to transform block 7),
		1, exp(-i* 7*twopi/32) =       ( s32_1,-c32_1), exp(-i*14*twopi/32) =       (-c    ,-s    ), exp(-i*21*twopi/32) =       (-s32_3, c32_3) (for inputs to transform block 8),
		 and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.	*/

	/*...Block 1: t00,t10,t20,t30	*/

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p10;
		prefetch_p_doubles(addp);
	#endif
			rt =t10;	t10=t00-rt;	t00=t00+rt;
			it =t11;	t11=t01-it;	t01=t01+it;

			rt =t30;	t30=t20-rt;	t20=t20+rt;
			it =t31;	t31=t21-it;	t21=t21+it;

			a[jt    ]=t00+t20;		a[jp    ]=t01+t21;
			t00      =t00-t20;		t01        =t01-t21;
			a[jt+p10]=t00*c10+t01*s10;	a[jp+p10]=t01*c10-t00*s10;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t10+t31;		it         =t11-t30;	/* mpy by E^-4 = -I is inlined here...	*/
			t10      =t10-t31;		t11        =t11+t30;
			a[jt+p08]=rt *c08+it *s08;	a[jp+p08]=it *c08-rt *s08;
			a[jt+p18]=t10*c18+t11*s18;	a[jp+p18]=t11*c18-t10*s18;

	/*...Block 5: t08,t18,t28,t38	*/

		jt = j1 + p04;
		jp = j2 + p04;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t18;	t18=t08-t19;	t08=t08+t19;		/* twiddle mpy by E^8 =-I	*/
				t19=t09+rt;	t09=t09-rt;

			rt =(t29+t28)*ISRT2;	t29=(t29-t28)*ISRT2;		t28=rt;	/* twiddle mpy by E^-4	*/
			rt =(t38-t39)*ISRT2;	it =(t38+t39)*ISRT2;			/* twiddle mpy by E^4 = -E^-12 is here...	*/
			t38=t28+rt;			t28=t28-rt;				/* ...and get E^-12 by flipping signs here.	*/
			t39=t29+it;			t29=t29-it;

			rt       =t08+t28;		it         =t09+t29;
			t08      =t08-t28;		t09        =t09-t29;
			a[jt    ]=rt *c04+it *s04;	a[jp    ]=it *c04-rt *s04;
			a[jt+p10]=t08*c14+t09*s14;	a[jp+p10]=t09*c14-t08*s14;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t18+t39;		it         =t19-t38;	/* mpy by E^-4 = -I is inlined here...	*/
			t18      =t18-t39;		t19        =t19+t38;
			a[jt+p08]=rt *c0C+it *s0C;	a[jp+p08]=it *c0C-rt *s0C;
			a[jt+p18]=t18*c1C+t19*s1C;	a[jp+p18]=t19*c1C-t18*s1C;

	/*...Block 3: t04,t14,t24,t34	*/

		jt = j1 + p02;
		jp = j2 + p02;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =(t15+t14)*ISRT2;	it =(t15-t14)*ISRT2;			/* twiddle mpy by E^-4	*/
			t14=t04-rt;			t04=t04+rt;
			t15=t05-it;			t05=t05+it;

			rt =t24*c + t25*s;		t25=t25*c - t24*s;		t24=rt;	/* twiddle mpy by E^-2	*/
			rt =t34*s + t35*c;		it =t35*s - t34*c;			/* twiddle mpy by E^-6	*/
			t34=t24-rt;			t24=t24+rt;
			t35=t25-it;			t25=t25+it;

			rt       =t04+t24;		it         =t05+t25;
			t04      =t04-t24;		t05        =t05-t25;
			a[jt    ]=rt *c02+it *s02;	a[jp    ]=it *c02-rt *s02;
			a[jt+p10]=t04*c12+t05*s12;	a[jp+p10]=t05*c12-t04*s12;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t14+t35;		it         =t15-t34;	/* mpy by E^-4 = -I is inlined here...	*/
			t14      =t14-t35;		t15        =t15+t34;
			a[jt+p08]=rt *c0A+it *s0A;	a[jp+p08]=it *c0A-rt *s0A;
			a[jt+p18]=t14*c1A+t15*s1A;	a[jp+p18]=t15*c1A-t14*s1A;

	/*...Block 7: t0C,t1C,t2C,t3C	*/

		jt = j1 + p06;
		jp = j2 + p06;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =(t1C-t1D)*ISRT2;	it =(t1C+t1D)*ISRT2;			/* twiddle mpy by E^4 = -E^-12 is here...	*/
			t1C=t0C+rt;			t0C=t0C-rt;				/* ...and get E^-12 by flipping signs here.	*/
			t1D=t0D+it;			t0D=t0D-it;

			rt =t2C*s + t2D*c;		t2D=t2D*s - t2C*c;		t2C=rt;	/* twiddle mpy by E^-6	*/
			rt =t3C*c + t3D*s;		it =t3D*c - t3C*s;			/* twiddle mpy by E^-18 is here...	*/
			t3C=t2C+rt;			t2C=t2C-rt;				/* ...and get E^-18 by flipping signs here.	*/
			t3D=t2D+it;			t2D=t2D-it;

			rt       =t0C+t2C;		it         =t0D+t2D;
			t0C      =t0C-t2C;		t0D        =t0D-t2D;
			a[jt    ]=rt *c06+it *s06;	a[jp    ]=it *c06-rt *s06;
			a[jt+p10]=t0C*c16+t0D*s16;	a[jp+p10]=t0D*c16-t0C*s16;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t1C+t3D;		it         =t1D-t3C;	/* mpy by E^-4 = -I is inlined here...	*/
			t1C      =t1C-t3D;		t1D        =t1D+t3C;
			a[jt+p08]=rt *c0E+it *s0E;	a[jp+p08]=it *c0E-rt *s0E;
			a[jt+p18]=t1C*c1E+t1D*s1E;	a[jp+p18]=t1D*c1E-t1C*s1E;

	/*...Block 2: t02,t12,t22,t32	*/

		jt = j1 + p01;
		jp = j2 + p01;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p18;
		prefetch_p_doubles(addp);
	#endif
			rt =t12*c + t13*s;		it =t13*c - t12*s;			/* twiddle mpy by E^-2	*/
			t12=t02-rt;			t02=t02+rt;
			t13=t03-it;			t03=t03+it;

			rt =t22*c32_1 + t23*s32_1;	t23=t23*c32_1 - t22*s32_1;	t22=rt;	/* twiddle mpy by E^-1	*/
			rt =t32*c32_3 + t33*s32_3;	it =t33*c32_3 - t32*s32_3;		/* twiddle mpy by E^-3	*/
			t32=t22-rt;			t22=t22+rt;
			t33=t23-it;			t23=t23+it;

			rt       =t02+t22;		it         =t03+t23;
			t02      =t02-t22;		t03        =t03-t23;
			a[jt    ]=rt *c01+it *s01;	a[jp    ]=it *c01-rt *s01;
			a[jt+p10]=t02*c11+t03*s11;	a[jp+p10]=t03*c11-t02*s11;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t12+t33;		it         =t13-t32;	/* mpy by E^-4 = -I is inlined here...	*/
			t12      =t12-t33;		t13        =t13+t32;
			a[jt+p08]=rt *c09+it *s09;	a[jp+p08]=it *c09-rt *s09;
			a[jt+p18]=t12*c19+t13*s19;	a[jp+p18]=t13*c19-t12*s19;

	/*...Block 6: t0A,t1A,t2A,t3A	*/

		jt = j1 + p05;
		jp = j2 + p05;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1A*s - t1B*c;		it =t1B*s + t1A*c;			/* twiddle mpy by -E^-10 is here...	*/
			t1A=t0A+rt;			t0A =t0A-rt;				/* ...and get E^-10 by flipping signs here.	*/
			t1B=t0B+it;			t0B =t0B-it;

			rt =t2A*s32_3 + t2B*c32_3;	t2B=t2B*s32_3 - t2A*c32_3;	t2A=rt;	/* twiddle mpy by E^-5	*/
			rt =t3A*c32_1 - t3B*s32_1;	it =t3B*c32_1 + t3A*s32_1;		/* twiddle mpy by -E^-15 is here...	*/
			t3A=t2A+rt;			t2A=t2A-rt;				/* ...and get E^-15 by flipping signs here.	*/
			t3B=t2B+it;			t2B=t2B-it;

			rt       =t0A+t2A;		it         =t0B+t2B;
			t0A      =t0A-t2A;		t0B        =t0B-t2B;
			a[jt    ]=rt *c05+it *s05;	a[jp    ]=it *c05-rt *s05;
			a[jt+p10]=t0A*c15+t0B*s15;	a[jp+p10]=t0B*c15-t0A*s15;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t1A+t3B;		it         =t1B-t3A;	/* mpy by E^-4 = -I is inlined here...	*/
			t1A      =t1A-t3B;		t1B        =t1B+t3A;
			a[jt+p08]=rt *c0D+it *s0D;	a[jp+p08]=it *c0D-rt *s0D;
			a[jt+p18]=t1A*c1D+t1B*s1D;	a[jp+p18]=t1B*c1D-t1A*s1D;

	/*...Block 4: t06,t16,t26,t36	*/

		jt = j1 + p03;
		jp = j2 + p03;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t16*s + t17*c;		it =t17*s - t16*c;			/* twiddle mpy by E^-6	*/
			t16=t06-rt;			t06 =t06+rt;
			t17=t07-it;			t07 =t07+it;

			rt =t26*c32_3 + t27*s32_3;	t27=t27*c32_3 - t26*s32_3;	t26=rt;	/* twiddle mpy by E^-3	*/
			rt =t36*s32_1 - t37*c32_1;	it =t37*s32_1 + t36*c32_1;		/* twiddle mpy by -E^-9 is here...	*/
			t36=t26+rt;			t26=t26-rt;				/* ...and get E^-9 by flipping signs here.	*/
			t37=t27+it;			t27=t27-it;

			rt       =t06+t26;		it         =t07+t27;
			t06      =t06-t26;		t07        =t07-t27;
			a[jt    ]=rt *c03+it *s03;	a[jp    ]=it *c03-rt *s03;
			a[jt+p10]=t06*c13+t07*s13;	a[jp+p10]=t07*c13-t06*s13;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t16+t37;		it         =t17-t36;	/* mpy by E^-4 = -I is inlined here...	*/
			t16      =t16-t37;		t17        =t17+t36;
			a[jt+p08]=rt *c0B+it *s0B;	a[jp+p08]=it *c0B-rt *s0B;
			a[jt+p18]=t16*c1B+t17*s1B;	a[jp+p18]=t17*c1B-t16*s1B;

	/*...Block 8: t0E,t1E,t2E,t3E	*/

		jt = j1 + p07;
		jp = j2 + p07;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1E*c - t1F*s;		it =t1F*c + t1E*s;			/* twiddle mpy by -E^-14 is here...	*/
			t1E=t0E+rt;			t0E =t0E-rt;				/* ...and get E^-14 by flipping signs here.	*/
			t1F=t0F+it;			t0F =t0F-it;

			rt =t2E*s32_1 + t2F*c32_1;	t2F=t2F*s32_1 - t2E*c32_1;	t2E=rt;	/* twiddle mpy by E^-7	*/
			rt =t3E*s32_3 + t3F*c32_3;	it =t3F*s32_3 - t3E*c32_3;		/* twiddle mpy by -E^-21 is here...	*/
			t3E=t2E+rt;			t2E=t2E-rt;				/* ...and get E^-21 by flipping signs here.	*/
			t3F=t2F+it;			t2F=t2F-it;

			rt       =t0E+t2E;		it         =t0F+t2F;
			t0E      =t0E-t2E;		t0F        =t0F-t2F;
			a[jt    ]=rt *c07+it *s07;	a[jp    ]=it *c07-rt *s07;
			a[jt+p10]=t0E*c17+t0F*s17;	a[jp+p10]=t0F*c17-t0E*s17;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t1E+t3F;		it         =t1F-t3E;	/* mpy by E^-4 = -I is inlined here...	*/
			t1E      =t1E-t3F;		t1F        =t1F+t3E;
			a[jt+p08]=rt *c0F+it *s0F;	a[jp+p08]=it *c0F-rt *s0F;
			a[jt+p18]=t1E*c1F+t1F*s1F;	a[jp+p18]=t1F*c1F-t1E*s1F;

#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */
}

