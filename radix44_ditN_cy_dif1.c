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

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

#ifdef USE_SSE2

	const int radix44_creals_in_local_store = 244;

  #ifdef USE_PTHREAD

	#ifndef USE_SSE2
		#error Pthreading only available in SSE2 mode!
	#endif

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
		int tid;
		int ndivr;
	
		int khi;
		int i;
		int jstart;
		int jhi;
		int col;
		int co2;
		int co3;
		int sw;
		int nwt;

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *s1p00r;
		struct complex *half_arr;

		int bjmodn00;
		int bjmodn01;
		int bjmodn02;
		int bjmodn03;
		int bjmodn04;
		int bjmodn05;
		int bjmodn06;
		int bjmodn07;
		int bjmodn08;
		int bjmodn09;
		int bjmodn10;
		int bjmodn11;
		int bjmodn12;
		int bjmodn13;
		int bjmodn14;
		int bjmodn15;
		int bjmodn16;
		int bjmodn17;
		int bjmodn18;
		int bjmodn19;
		int bjmodn20;
		int bjmodn21;
		int bjmodn22;
		int bjmodn23;
		int bjmodn24;
		int bjmodn25;
		int bjmodn26;
		int bjmodn27;
		int bjmodn28;
		int bjmodn29;
		int bjmodn30;
		int bjmodn31;
		int bjmodn32;
		int bjmodn33;
		int bjmodn34;
		int bjmodn35;
		int bjmodn36;
		int bjmodn37;
		int bjmodn38;
		int bjmodn39;
		int bjmodn40;
		int bjmodn41;
		int bjmodn42;
		int bjmodn43;
		/* carries: */
		double cy00;
		double cy01;
		double cy02;
		double cy03;
		double cy04;
		double cy05;
		double cy06;
		double cy07;
		double cy08;
		double cy09;
		double cy10;
		double cy11;
		double cy12;
		double cy13;
		double cy14;
		double cy15;
		double cy16;
		double cy17;
		double cy18;
		double cy19;
		double cy20;
		double cy21;
		double cy22;
		double cy23;
		double cy24;
		double cy25;
		double cy26;
		double cy27;
		double cy28;
		double cy29;
		double cy30;
		double cy31;
		double cy32;
		double cy33;
		double cy34;
		double cy35;
		double cy36;
		double cy37;
		double cy38;
		double cy39;
		double cy40;
		double cy41;
		double cy42;
		double cy43;
	};

  #endif

//	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_MSVC)

		/*...Radix-11 DFT: Inputs in memory locations __I0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA],\
		where r0 is a memory address and the i's are literal [byte] offsets. Outputs similarly go into memory locations\
		__O0 + [__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA], assumed disjoint with inputs:\

                Total SSE opcounts: 118 MOVAPS, 182 ADD/SUBPD, 44 MULPD
		*/\
		#define SSE2_RADIX_11_DFT(__I0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA, __cc, __O0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA)\
		{\
		/*\
			t3  = __A1r +__Aar;			t4  = __A1i + __Aai;	// x1 + x10	; load a1-a4 into xmm1-4, a7-aa into xmm0,7,6,5 (in that order), \
			t5  = __A2r +__A9r;			t6  = __A2i + __A9i;	// x2 + x9	; compute t21,19,17,15 (output in xmm1-4, dump to memory), \
			t7  = __A3r +__A8r;			t8  = __A3i + __A8i;	// x3 + x8	; and t3,5,7,9 (output in xmm1-4), \
			t9  = __A4r +__A7r;			t10 = __A4i + __A7i;	// x4 + x7	; load a5,6 into xmm5,6, compute t11,13, \
			t11 = __A5r +__A6r;			t12 = __A5i + __A6i;	// x5 + x6	; dump t13 to memory. \
			t13 = __A5r -__A6r;			t14 = __A5i - __A6i;	// x5 - x6	\
			t15 = __A4r -__A7r;			t16 = __A4i - __A7i;	// x4 - x7	\
			t17 = __A3r -__A8r;			t18 = __A3i - __A8i;	// x3 - x8	\
			t19 = __A2r -__A9r;			t20 = __A2i - __A9i;	// x2 - x9	\
			t21 = __A1r -__Aar;			t22 = __A1i - __Aai;	// x1 - x10	\
			(t1)=__B0r = __A0r;			(t2)=__B0i = __A0i;		// x0, store in xmm0 \

		We calculate the cr and ci terms below first and it is natural to store these intermediates in the B[1-5] output slots
		and the t13,15,17,19,21 and t14,16,18,20,22 terms (needed for the computation of the sr and si terms) in the B[6-a]r,i output slots
		(resp.) to free up registers for the sr and si computations. But if we e.g. compute sr1-5 and combine with ci1-5 to get the Bi outputs,
		writing B[1-a]i (or more specifically B[6-a]i, we end up overwriting the t14-22 terms stored there which are still needed
		for the computation of the si's (and hence the Br-outputs.)

			B1r = cr1 - si1;			B1i = ci1 + sr1;	// X1 = C1 + I*S1	\
			B2r = cr2 - si2;			B2i = ci2 + sr2;	// X2 = C2 + I*S2	\
			B3r = cr3 - si3;			B3i = ci3 + sr3;	// X3 = C3 + I*S3	\
			B4r = cr4 - si4;			B4i = ci4 + sr4;	// X4 = C4 + I*S4	\
			B5r = cr5 - si5;			B5i = ci5 + sr5;	// X5 = C5 + I*S5	\
			B6r = cr5 + si5;			B6i = ci5 - sr5;	// X6 =	C5 - I*S5	\
			B7r = cr4 + si4;			B7i = ci4 - sr4;	// X7 =	C4 - I*S4	\
			B8r = cr3 + si3;			B8i = ci3 - sr3;	// X8 =	C3 - I*S3	\
			B9r = cr2 + si2;			B9i = ci2 - sr2;	// X9 =	C2 - I*S2	\
			Bar = cr1 + si1;			Bai = ci1 - sr1;	// X10=	C1 - I*S1	\

		SOLUTION: Instead, store the t13,15,17,19,21 and t14,16,18,20,22 terms in the B[6-a]i,r output slots, respectively.
		The resulting compute sequence then is:

			1. Compute cr1-5, store in B[1-5]r;
			2. Compute t13-21,store in B[6-a]i;
			3. Compute ci1-5, store in B[1-5]i;
			4. Compute t14-22,store in B[6-a]r;
			5. Compute sr1-5 (using the t13-21 terms stored in B[6-a]i, combine with the ci1-5 terms stored in B[1-5]i, write results to B[0-a]i;
			6. Compute si1-5 (using the t14-22 terms stored in B[6-a]r, combine with the cr1-5 terms stored in B[1-5]r, write results to B[0-a]r;
		*/\
		/********************************************/\
		/*       Here are the 5 cosine terms:       */\
		/********************************************/\
			__asm	mov	eax, __I0	\
			__asm	mov	ebx, __cc	\
			__asm	mov	ecx, __O0	\
			__asm	movaps	xmm1,[eax+__i1]	/* __A1r */\
			__asm	movaps	xmm5,[eax+__iA]	/* __Aar */\
			__asm	movaps	xmm2,[eax+__i2]	/* __A2r */\
			__asm	movaps	xmm6,[eax+__i9]	/* __A9r */\
			__asm	movaps	xmm3,[eax+__i3]	/* __A3r */\
			__asm	movaps	xmm7,[eax+__i8]	/* __A8r */\
			__asm	movaps	xmm4,[eax+__i4]	/* __A4r */\
			__asm	movaps	xmm0,[eax+__i7]	/* __A7r */\
			__asm	add	ecx, 0x10	/* t13-21 stored in B[6-a]i, not B[6-a]r. */\
			__asm	subpd	xmm1,xmm5	\
			__asm	subpd	xmm2,xmm6	\
			__asm	subpd	xmm3,xmm7	\
			__asm	subpd	xmm4,xmm0	\
			__asm	movaps	[ecx+__oA],xmm1	/* t21 */\
			__asm	movaps	[ecx+__o9],xmm2	/* t19 */\
			__asm	movaps	[ecx+__o8],xmm3	/* t17 */\
			__asm	movaps	[ecx+__o7],xmm4	/* t15 */\
			__asm	addpd	xmm5,xmm5	\
			__asm	addpd	xmm6,xmm6	\
			__asm	addpd	xmm7,xmm7	\
			__asm	addpd	xmm0,xmm0	\
			__asm	addpd	xmm1,xmm5	/* t3  */\
			__asm	addpd	xmm2,xmm6	/* t5  */\
			__asm	movaps	xmm5,[eax+__i5]	/* __A5r */\
			__asm	movaps	xmm6,[eax+__i6]	/* __A6r */\
			__asm	addpd	xmm3,xmm7	/* t7  */\
			__asm	addpd	xmm4,xmm0	/* t9  */\
			__asm	subpd	xmm5,xmm6	/* t13 */\
			__asm	movaps	[ecx+__o6],xmm5	\
			__asm	addpd	xmm6,xmm6	\
			__asm	movaps	xmm0,[eax     ]	/* t1/__B0r */\
			__asm	addpd	xmm5,xmm6	/* t11 */\
		/********************************************/\
		/*               Real Parts:                */\
		/********************************************/\
			__asm	sub	ecx, 0x10	/* Now switch ptr from B[6-a]i to B[6-a]r. */\
		/* b0/t1,t3,t5,t7,t9,t11 in xmm0-5, xmm6,7 free \
			c1 = t3-t5;					s1 = t4-t6;		// store in c1/t3, make 2 copies (call them u3,v3), mul a0*c1 and store \
			c2 = t11-t5;				s2 = t12-t6;	// store in c2/t11 \
			c4 = t7-t5;					s4 = t8-t6;		// store in c4/t7 \
			c5 = t9-t5;					s5 = t10-t6;	// store in c5/t9, then mul 5*t5 \
		*/\
			__asm	subpd	xmm1,xmm2	/* c1/t3  */\
			__asm	subpd	xmm5,xmm2	/* c2/t11 */\
			__asm	movaps	xmm6,xmm1	/* u3 = cpy c1 */\
			__asm	subpd	xmm3,xmm2	/* c4/t7  */\
			__asm	movaps	xmm7,xmm1	/* v3 = cpy c1 */\
			__asm	mulpd	xmm1,[ebx     ]	/* a0*c1 */\
			__asm	movaps	[ecx+__o1],xmm1	/* store a0*c1; xmm1 FREE */\
			__asm	subpd	xmm4,xmm2	/* c5/t9  */\
			__asm	mulpd	xmm2,[ebx-0x10]	/* 5*t5 */\
		/*	c3 = c1+c2;					s3 = s1+s2;		// t3+t11-2*t5; store in c3/u3 */\
			__asm	addpd	xmm6,xmm5	/* c3  */\
		/*	c7 = c1-c4;					s7 = s1-s4;		// t3-t7; store in v3, mul a6*c7 and store. */\
			__asm	subpd	xmm7,xmm3	/* c7  */\
			__asm	mulpd	xmm7,[ebx+0x60]	/* a6*c7 */\
			__asm	movaps	[ecx+__o3],xmm7	/* store a6*c7; xmm7 FREE */\
		/*	c6 = c4+c5;					s6 = s4+s5;		// copy c4, mul a3*c4 and store. */\
			__asm	movaps	xmm1,xmm3	/* copy c4 */\
			__asm	addpd	xmm3,xmm4	/* c6  */\
			__asm	mulpd	xmm1,[ebx+0x30]	/* a3*c4 */\
			__asm	movaps	[ecx+__o2],xmm1	/* store a3*c4; xmm1 FREE */\
		/*	c8 = c2-c5;					s8 = s2-s5;		// t11-t9; copy c2, mul a7*c8 and store. */\
			__asm	movaps	xmm7,xmm5	/* copy c2 */\
			__asm	subpd	xmm5,xmm4	/* c8  */\
			__asm	mulpd	xmm4,[ebx+0x40]	/* a4*c5 */\
			__asm	mulpd	xmm5,[ebx+0x70]	/* a7*c8 */\
			__asm	mulpd	xmm7,[ebx+0x10]	/* a1*c2 */\
		/*	c9 = c3-c6;					s9 = s3-s6;		// copy c3, mul a8*c9 and store. */\
			__asm	addpd	xmm2,xmm3	/* 5*t5+c6 */\
			__asm	movaps	xmm1,xmm6	/* copy c3 */\
			__asm	subpd	xmm6,xmm3	/* c9  */\
			__asm	mulpd	xmm6,[ebx+0x80]	/* c9 = a8*c9 */\
			__asm	mulpd	xmm3,[ebx+0x50]	/* c6 = a5*c6 */\
		/*	c10= c3+c6+5*t5;			s10= s3+s6+5*t6 // c10= t3+t5+t7+t9+t11;	 s10= t4+t6+t8+t10+t12; store in t6. */\
			__asm	addpd	xmm2,xmm1	/* c10 */\
			__asm	mulpd	xmm1,[ebx+0x20]	/* c3 = a2*c3 */\
		/*\
			__B0r += c10;				__B0i += s10;	// X0/t2, store result to memory \
			c3 = a2*c3;					s3 = a2*s3;	\
			c6 = a5*c6;					s6 = a5*s6;	\
			c9 = a8*c9;					s9 = a8*s9;	\
			c10= a9*c10+__B0r;			s10= a9*s10+__B0i;	\
		*/\
			__asm	addpd	xmm0,xmm2	/* __B0r */\
			__asm	mulpd	xmm2,[ebx+0x90]	/* a9*c10*/\
			__asm	movaps	[ecx     ],xmm0	/* store __B0r */\
			__asm	addpd	xmm2,xmm0	/* c10; xmm0 FREE */\
		\
		/* a[0,3,6]*c[1,4,7] in ecx+__o[1,2,3]; a[1,4,7]*c[2,5,8] in xmm[7,4,5]; c[3,6,9] in xmm[1,3,6]; c10 in xmm2; xmm0 free */\
			__asm	addpd	xmm7,xmm1		/* c1 = a1*c2+c3 */\
			__asm	addpd	xmm1,[ecx+__o1]	/* c2 = a0*c1+c3 */\
			__asm	addpd	xmm4,xmm3		/* c5 = a4*c5+c6 */\
			__asm	addpd	xmm3,[ecx+__o2]	/* c4 = a3*c4+c6 */\
			__asm	addpd	xmm5,xmm6		/* c8 = a7*c8+c9 */\
			__asm	addpd	xmm6,[ecx+__o3]	/* c7 = a6*c7+c9 */\
		/*\
			cr1 = c10+c1-c7;			ci1 = s10+s1-s7;\
			cr2 = c10-c1-c2-c4-c5;		ci2 = s10-s1-s2-s4-s5;\
			cr3 = c10+c4+c7;			ci3 = s10+s4+s7;\
			cr4 = c10+c5+c8;			ci4 = s10+s5+s8;\
			cr5 = c10+c2-c8;			ci5 = s10+s2-s8;\
		*/\
			__asm	movaps	xmm0,xmm2	/* copy c10*/\
			__asm	subpd	xmm2,xmm1	/* c10-c1 */\
			__asm	addpd	xmm1,xmm0	/* c10+c1 */\
			__asm	subpd	xmm2,xmm7	/* c10-c1-c2 */\
			__asm	addpd	xmm7,xmm0	/* c10+c2 */\
			__asm	subpd	xmm1,xmm6	/* cr1 = c10+c1-c7 */\
			__asm	subpd	xmm2,xmm3	/* c10-c1-c2-c4 */\
			__asm	movaps	[ecx+__o1],xmm1	/* store cr1 */\
			__asm	addpd	xmm3,xmm0	/* c10+c4 */\
			__asm	subpd	xmm7,xmm5	/* cr5 = c10+c2-c8 */\
			__asm	subpd	xmm2,xmm4	/* cr2 = c10-c1-c2-c4-c5 */\
			__asm	movaps	[ecx+__o5],xmm7	/* store cr5 */\
			__asm	movaps	[ecx+__o2],xmm2	/* store cr2 */\
			__asm	addpd	xmm4,xmm0	/* c10+c5 */\
			__asm	addpd	xmm3,xmm6	/* cr3 = c10+c4+c7 */\
			__asm	addpd	xmm4,xmm5	/* cr4 = c10+c5+c8 */\
			__asm	movaps	[ecx+__o3],xmm3	/* store cr3 */\
			__asm	movaps	[ecx+__o4],xmm4	/* store cr4 */\
		/********************************************/\
		/*          Imaginary Parts:                */\
		/********************************************/\
			__asm	add	eax, 0x10	\
			/* Keep ecx pointing to real outputs, since t14-22 get stored in B[6-a]r, not B[6-a]i. */\
			__asm	movaps	xmm1,[eax+__i1]	/* __A1r */\
			__asm	movaps	xmm5,[eax+__iA]	/* __Aar */\
			__asm	movaps	xmm2,[eax+__i2]	/* __A2r */\
			__asm	movaps	xmm6,[eax+__i9]	/* __A9r */\
			__asm	movaps	xmm3,[eax+__i3]	/* __A3r */\
			__asm	movaps	xmm7,[eax+__i8]	/* __A8r */\
			__asm	movaps	xmm4,[eax+__i4]	/* __A4r */\
			__asm	movaps	xmm0,[eax+__i7]	/* __A7r */\
			__asm	subpd	xmm1,xmm5	\
			__asm	subpd	xmm2,xmm6	\
			__asm	subpd	xmm3,xmm7	\
			__asm	subpd	xmm4,xmm0	\
			__asm	movaps	[ecx+__oA],xmm1	/* t21 */\
			__asm	movaps	[ecx+__o9],xmm2	/* t19 */\
			__asm	movaps	[ecx+__o8],xmm3	/* t17 */\
			__asm	movaps	[ecx+__o7],xmm4	/* t15 */\
			__asm	addpd	xmm5,xmm5	\
			__asm	addpd	xmm6,xmm6	\
			__asm	addpd	xmm7,xmm7	\
			__asm	addpd	xmm0,xmm0	\
			__asm	addpd	xmm1,xmm5	/* t3  */\
			__asm	addpd	xmm2,xmm6	/* t5  */\
			__asm	movaps	xmm5,[eax+__i5]	/* __A5r */\
			__asm	movaps	xmm6,[eax+__i6]	/* __A6r */\
			__asm	addpd	xmm3,xmm7	/* t7  */\
			__asm	addpd	xmm4,xmm0	/* t9  */\
			__asm	subpd	xmm5,xmm6	/* t13 */\
			__asm	movaps	[ecx+__o6],xmm5	\
			__asm	addpd	xmm6,xmm6	\
			__asm	movaps	xmm0,[eax     ]	/* t1/__B0r */\
			__asm	addpd	xmm5,xmm6	/* t11 */\
			__asm	add	ecx, 0x10	/* Now switch ptr from B[6-a]r to B[6-a]i. */\
		/* b0/t1,t3,t5,t7,t9,t11 in xmm0-5, xmm6,7 free */\
			__asm	subpd	xmm1,xmm2	/* c1/t3  */\
			__asm	subpd	xmm5,xmm2	/* c2/t11 */\
			__asm	movaps	xmm6,xmm1	/* u3 = cpy c1 */\
			__asm	subpd	xmm3,xmm2	/* c4/t7  */\
			__asm	movaps	xmm7,xmm1	/* v3 = cpy c1 */\
			__asm	mulpd	xmm1,[ebx     ]	/* a0*c1 */\
			__asm	movaps	[ecx+__o1],xmm1	/* store a0*c1; xmm1 FREE */\
			__asm	subpd	xmm4,xmm2	/* c5/t9  */\
			__asm	mulpd	xmm2,[ebx-0x10]	/* 5*t5 */\
			__asm	addpd	xmm6,xmm5	/* c3  */\
			__asm	subpd	xmm7,xmm3	/* c7  */\
			__asm	mulpd	xmm7,[ebx+0x60]	/* a6*c7 */\
			__asm	movaps	[ecx+__o3],xmm7	/* store a6*c7; xmm7 FREE */\
			__asm	movaps	xmm1,xmm3	/* copy c4 */\
			__asm	addpd	xmm3,xmm4	/* c6  */\
			__asm	mulpd	xmm1,[ebx+0x30]	/* a3*c4 */\
			__asm	movaps	[ecx+__o2],xmm1	/* store a3*c4; xmm1 FREE */\
			__asm	movaps	xmm7,xmm5	/* copy c2 */\
			__asm	subpd	xmm5,xmm4	/* c8  */\
			__asm	mulpd	xmm4,[ebx+0x40]	/* a4*c5 */\
			__asm	mulpd	xmm5,[ebx+0x70]	/* a7*c8 */\
			__asm	mulpd	xmm7,[ebx+0x10]	/* a1*c2 */\
			__asm	addpd	xmm2,xmm3	/* 5*t5+c6 */\
			__asm	movaps	xmm1,xmm6	/* copy c3 */\
			__asm	subpd	xmm6,xmm3	/* c9  */\
			__asm	mulpd	xmm6,[ebx+0x80]	/* c9 = a8*c9 */\
			__asm	mulpd	xmm3,[ebx+0x50]	/* c6 = a5*c6 */\
			__asm	addpd	xmm2,xmm1	/* c10 */\
			__asm	mulpd	xmm1,[ebx+0x20]	/* c3 = a2*c3 */\
			__asm	addpd	xmm0,xmm2	/* __B0r */\
			__asm	mulpd	xmm2,[ebx+0x90]	/* a9*c10*/\
			__asm	movaps	[ecx     ],xmm0	/* store __B0r */\
			__asm	addpd	xmm2,xmm0	/* c10; xmm0 FREE */\
		/* a[0,3,6]*c[1,4,7] in ecx+__o[1,2,3]; a[1,4,7]*c[2,5,8] in xmm[7,4,5]; c[3,6,9] in xmm[1,3,6]; c10 in xmm2; xmm0 free */\
			__asm	addpd	xmm7,xmm1		/* c1 = a1*c2+c3 */\
			__asm	addpd	xmm1,[ecx+__o1]	/* c2 = a0*c1+c3 */\
			__asm	addpd	xmm4,xmm3		/* c5 = a4*c5+c6 */\
			__asm	addpd	xmm3,[ecx+__o2]	/* c4 = a3*c4+c6 */\
			__asm	addpd	xmm5,xmm6		/* c8 = a7*c8+c9 */\
			__asm	addpd	xmm6,[ecx+__o3]	/* c7 = a6*c7+c9 */\
			__asm	movaps	xmm0,xmm2	/* copy c10*/\
			__asm	subpd	xmm2,xmm1	/* c10-c1 */\
			__asm	addpd	xmm1,xmm0	/* c10+c1 */\
			__asm	subpd	xmm2,xmm7	/* c10-c1-c2 */\
			__asm	addpd	xmm7,xmm0	/* c10+c2 */\
			__asm	subpd	xmm1,xmm6	/* cr1 = c10+c1-c7 */\
			__asm	subpd	xmm2,xmm3	/* c10-c1-c2-c4 */\
			__asm	movaps	[ecx+__o1],xmm1	/* store cr1 */\
			__asm	addpd	xmm3,xmm0	/* c10+c4 */\
			__asm	subpd	xmm7,xmm5	/* cr5 = c10+c2-c8 */\
			__asm	subpd	xmm2,xmm4	/* cr2 = c10-c1-c2-c4-c5 */\
			__asm	movaps	[ecx+__o5],xmm7	/* store cr5 */\
			__asm	movaps	[ecx+__o2],xmm2	/* store cr2 */\
			__asm	addpd	xmm4,xmm0	/* c10+c5 */\
			__asm	addpd	xmm3,xmm6	/* cr3 = c10+c4+c7 */\
			__asm	addpd	xmm4,xmm5	/* cr4 = c10+c5+c8 */\
			__asm	movaps	[ecx+__o3],xmm3	/* store cr3 */\
			__asm	movaps	[ecx+__o4],xmm4	/* store cr4 */\
	\
		/************************************************************************************************************************************/\
		/* Here are the 5 sine terms: Similar sequence to cosine terms, but with t21,19,17,15,13 replacing t3,5,7,9,11 and some sign flips: */\
		/************************************************************************************************************************************/\
			__asm	sub	eax, 0x10	\
			__asm	add	ebx, 0xa0	/* Increment sincos ptr to point to B-terms; 5-constant will be -0xa0 offset w.r.to ebx... */\
			/* Keep ecx pointing to imag outputs, since t13-21 are stored in B[6-a]i, not B[6-a]r. */\
		/********************************************/\
		/*               Real Parts:                */\
		/********************************************/\
			__asm	movaps	xmm1,[ecx+__oA]	/* t21 */\
			__asm	movaps	xmm2,[ecx+__o9]	/* t19 */\
			__asm	movaps	xmm3,[ecx+__o8]	/* t17 */\
			__asm	movaps	xmm4,[ecx+__o7]	/* t15 */\
			__asm	movaps	xmm5,[ecx+__o6]	/* t13 */\
			__asm	addpd	xmm1,xmm2	/* c1/t3  */\
			__asm	addpd	xmm5,xmm2	/* c2/t11 */\
			__asm	movaps	xmm6,xmm1	/* u3 = cpy c1 */\
			__asm	addpd	xmm3,xmm2	/* c4/t7  */\
			__asm	movaps	xmm7,xmm1	/* v3 = cpy c1 */\
			__asm	mulpd	xmm1,[ebx     ]	/* b0*c1 */\
			__asm	movaps	[ecx+__o6],xmm1	/* store b0*c1; xmm1 FREE */\
			__asm	addpd	xmm4,xmm2	/* c5/t9  */\
			__asm	mulpd	xmm2,[ebx-0xb0]	/* 5*t5 */\
		/*	c3 = c1+c2;					s3 = s1+s2;		 t3+t11-2*t5; store in c3/u3 */\
			__asm	addpd	xmm6,xmm5	/* c3  */\
		/*	c7 = c1-c4;					s7 = s1-s4;		 t3-t7; store in v3, mul b6*c7 and store. */\
			__asm	subpd	xmm7,xmm3	/* c7  */\
			__asm	mulpd	xmm7,[ebx+0x60]	/* b6*c7 */\
			__asm	movaps	[ecx+__o8],xmm7	/* store b6*c7; xmm7 FREE */\
		/*	c6 = c4+c5;					s6 = s4+s5;		 copy c4, mul b3*c4 and store. */\
			__asm	movaps	xmm1,xmm3	/* copy c4 */\
			__asm	addpd	xmm3,xmm4	/* c6  */\
			__asm	mulpd	xmm1,[ebx+0x30]	/* b3*c4 */\
			__asm	movaps	[ecx+__o7],xmm1	/* store b3*c4; xmm1 FREE */\
		/*	c8 = c2-c5;					s8 = s2-s5;		 t11-t9; copy c2, mul b7*c8 and store. */\
			__asm	movaps	xmm7,xmm5	/* copy c2 */\
			__asm	subpd	xmm5,xmm4	/* c8  */\
			__asm	mulpd	xmm4,[ebx+0x40]	/* b4*c5 */\
			__asm	mulpd	xmm5,[ebx+0x70]	/* b7*c8 */\
			__asm	mulpd	xmm7,[ebx+0x10]	/* b1*c2 */\
		/*	c9 = c3-c6;					s9 = s3-s6;		 copy c3, mul b8*c9 and store. */\
			__asm	subpd	xmm2,xmm3	/* 5*t5-c6 */\
			__asm	movaps	xmm1,xmm6	/* copy c3 */\
			__asm	subpd	xmm6,xmm3	/* c9  */\
			__asm	mulpd	xmm6,[ebx+0x80]	/* c9 = b8*c9 */\
			__asm	mulpd	xmm3,[ebx+0x50]	/* c6 = b5*c6 */\
		/*	c10= c3+c6-5*t5;			s10= s3+s6-5*t6 */\
			__asm	subpd	xmm2,xmm1	/* -c10 */\
			__asm	mulpd	xmm1,[ebx+0x20]	/* c3 = b2*c3 */\
		/*\
			c3 = b2*c3;					s3 = b2*s3;	\
			c6 = b5*c6;					s6 = b5*s6;	\
			c9 = b8*c9;					s9 = b8*s9;	\
			c10= (-b9)*(-c10);			s10= (-b9)*(-c10);	\
		*/\
			__asm	mulpd	xmm2,[ebx+0x90]	/* c10 = b9*c10*/\
		\
		/* a[0,3,6]*c[1,4,7] in ecx+__o[1,2,3]; a[1,4,7]*c[2,5,8] in xmm[7,4,5]; c[3,6,9] in xmm[1,3,6]; c10 in xmm2; xmm0 free */\
			__asm	addpd	xmm7,xmm1		/* c2 = b1*c2+c3 */\
			__asm	addpd	xmm1,[ecx+__o6]	/* c1 = b0*c1+c3 */\
			__asm	addpd	xmm4,xmm3		/* c5 = b4*c5+c6 */\
			__asm	addpd	xmm3,[ecx+__o7]	/* c4 = b3*c4+c6 */\
			__asm	addpd	xmm5,xmm6		/* c8 = b7*c8+c9 */\
			__asm	addpd	xmm6,[ecx+__o8]	/* c7 = b6*c7+c9 */\
		/*\
			sr1 = c10+c1-c7;			si1 = s10+s1-s7;\
			sr2 = c1+c2+c4+c5-c10;		si2 = s1+s2+s4+s5-s10;\
			sr3 = c10+c4+c7;			si3 = s10+s4+s7;\
			sr4 = c10+c5+c8;			si4 = s10+s5+s8;\
			sr5 = c10+c2-c8;			si5 = s10+s2-s8;\
		*/\
			__asm	xorpd	xmm0,xmm0	/* 0.0 */\
			__asm	subpd	xmm0,xmm2	/* copy (-c10) */\
			__asm	addpd	xmm0,xmm1	/* c1-c10 */\
			__asm	addpd	xmm1,xmm2	/* c10+c1 */\
			__asm	addpd	xmm0,xmm7	/* c1+c2-c10 */\
			__asm	addpd	xmm7,xmm2	/* c10+c2 */\
			__asm	subpd	xmm1,xmm6	/* sr1 = c10+c1-c7 */\
			__asm	addpd	xmm0,xmm3	/* c1+c2+c4-c10 */\
		/*	__asm	movaps	[ecx+__o6],xmm1	// store sr1 */\
			__asm	addpd	xmm3,xmm2	/* c10+c4 */\
			__asm	subpd	xmm7,xmm5	/* sr5 = c10+c2-c8 */\
			__asm	addpd	xmm0,xmm4	/* sr2 = c1+c2+c4+c5-c10 */\
		/*	__asm	movaps	[ecx+__oA],xmm7	// store sr5 */\
		/*	__asm	movaps	[ecx+__o7],xmm0	// store sr2 */\
			__asm	addpd	xmm4,xmm2	/* c10+c5 */\
			__asm	addpd	xmm3,xmm6	/* sr3 = c10+c4+c7 */\
			__asm	addpd	xmm4,xmm5	/* sr4 = c10+c5+c8 */\
		/*	__asm	movaps	[ecx+__o8],xmm3	// store sr3 */\
		/*	__asm	movaps	[ecx+__o9],xmm4	// store sr4 */\
		/*\
			sr1,2,3,4,5 in xmm1,0,3,4,7:

			B1i = ci1 + sr1;	// X1 = C1 + I*S1	\
			B2i = ci2 + sr2;	// X2 = C2 + I*S2	\
			B3i = ci3 + sr3;	// X3 = C3 + I*S3	\
			B4i = ci4 + sr4;	// X4 = C4 + I*S4	\
			B5i = ci5 + sr5;	// X5 = C5 + I*S5	\
			B6i = ci5 - sr5;	// X6 =	C5 - I*S5	\
			B7i = ci4 - sr4;	// X7 =	C4 - I*S4	\
			B8i = ci3 - sr3;	// X8 =	C3 - I*S3	\
			B9i = ci2 - sr2;	// X9 =	C2 - I*S2	\
			Bai = ci1 - sr1;	// X10=	C1 - I*S1	\
		*/\
		/* sr-terms get combined with ci's and output in Bi's, so no need to increment ecx: */\
		/* xmm2,5,6 FREE: */\
			__asm	movaps	xmm2,xmm1	/* cpy sr1 */\
			__asm	addpd	xmm1,[ecx+__o1]	/* B1i */\
			__asm	addpd	xmm2,xmm2	/* 2 x sr1 */\
			__asm	movaps	[ecx+__o1],xmm1	/* store B1i */\
			__asm	subpd	xmm1,xmm2		/* Bai */\
			__asm	movaps	[ecx+__oA],xmm1	/* store Bai */\
		\
			__asm	movaps	xmm5,xmm0	/* cpy sr2 */\
			__asm	addpd	xmm0,[ecx+__o2]	/* B2i */\
			__asm	addpd	xmm5,xmm5	/* 2 x sr2 */\
			__asm	movaps	[ecx+__o2],xmm0	/* store B2i */\
			__asm	subpd	xmm0,xmm5		/* B9i */\
			__asm	movaps	[ecx+__o9],xmm0	/* store B9i */\
		\
			__asm	movaps	xmm6,xmm3	/* cpy sr3 */\
			__asm	addpd	xmm3,[ecx+__o3]	/* B3i */\
			__asm	addpd	xmm6,xmm6	/* 2 x sr3 */\
			__asm	movaps	[ecx+__o3],xmm3	/* store B3i */\
			__asm	subpd	xmm3,xmm6		/* B8i */\
			__asm	movaps	[ecx+__o8],xmm3	/* store B8i */\
		\
			__asm	movaps	xmm2,xmm4	/* cpy sr4 */\
			__asm	addpd	xmm4,[ecx+__o4]	/* B4i */\
			__asm	addpd	xmm2,xmm2	/* 2 x sr4 */\
			__asm	movaps	[ecx+__o4],xmm4	/* store B4i */\
			__asm	subpd	xmm4,xmm2		/* B7i */\
			__asm	movaps	[ecx+__o7],xmm4	/* store B7i */\
		\
			__asm	movaps	xmm5,xmm7	/* cpy sr5 */\
			__asm	addpd	xmm7,[ecx+__o5]	/* B5i */\
			__asm	addpd	xmm5,xmm5	/* 2 x sr5 */\
			__asm	movaps	[ecx+__o5],xmm7	/* store B5i */\
			__asm	subpd	xmm7,xmm5		/* B6i */\
			__asm	movaps	[ecx+__o6],xmm7	/* store B6i */\
		\
		/********************************************/\
		/*          Imaginary Parts:                */\
		/********************************************/\
			/* Decrement ecx to point to real outputs, since t14-22 are stored in B[6-a]r, not B[6-a]i. */\
			__asm	sub	ecx, 0x10	\
			__asm	movaps	xmm1,[ecx+__oA]	/* t22 */\
			__asm	movaps	xmm2,[ecx+__o9]	/* t20 */\
			__asm	movaps	xmm3,[ecx+__o8]	/* t18 */\
			__asm	movaps	xmm4,[ecx+__o7]	/* t16 */\
			__asm	movaps	xmm5,[ecx+__o6]	/* t14 */\
			__asm	addpd	xmm1,xmm2	/* c1/t3  */\
			__asm	addpd	xmm5,xmm2	/* c2/t11 */\
			__asm	movaps	xmm6,xmm1	/* u3 = cpy c1 */\
			__asm	addpd	xmm3,xmm2	/* c4/t7  */\
			__asm	movaps	xmm7,xmm1	/* v3 = cpy c1 */\
			__asm	mulpd	xmm1,[ebx     ]	/* b0*c1 */\
			__asm	movaps	[ecx+__o6],xmm1	/* store b0*c1; xmm1 FREE */\
			__asm	addpd	xmm4,xmm2	/* c5/t9  */\
			__asm	mulpd	xmm2,[ebx-0xb0]	/* 5*t5 */\
		/*	c3 = c1+c2;					s3 = s1+s2;		 t3+t11-2*t5; store in c3/u3 */\
			__asm	addpd	xmm6,xmm5	/* c3  */\
		/*	c7 = c1-c4;					s7 = s1-s4;		 t3-t7; store in v3, mul b6*c7 and store. */\
			__asm	subpd	xmm7,xmm3	/* c7  */\
			__asm	mulpd	xmm7,[ebx+0x60]	/* b6*c7 */\
			__asm	movaps	[ecx+__o8],xmm7	/* store b6*c7; xmm7 FREE */\
		/*	c6 = c4+c5;					s6 = s4+s5;		 copy c4, mul b3*c4 and store. */\
			__asm	movaps	xmm1,xmm3	/* copy c4 */\
			__asm	addpd	xmm3,xmm4	/* c6  */\
			__asm	mulpd	xmm1,[ebx+0x30]	/* b3*c4 */\
			__asm	movaps	[ecx+__o7],xmm1	/* store b3*c4; xmm1 FREE */\
		/*	c8 = c2-c5;					s8 = s2-s5;		 t11-t9; copy c2, mul b7*c8 and store. */\
			__asm	movaps	xmm7,xmm5	/* copy c2 */\
			__asm	subpd	xmm5,xmm4	/* c8  */\
			__asm	mulpd	xmm4,[ebx+0x40]	/* b4*c5 */\
			__asm	mulpd	xmm5,[ebx+0x70]	/* b7*c8 */\
			__asm	mulpd	xmm7,[ebx+0x10]	/* b1*c2 */\
		/*	c9 = c3-c6;					s9 = s3-s6;		 copy c3, mul b8*c9 and store. */\
			__asm	subpd	xmm2,xmm3	/* 5*t5-c6 */\
			__asm	movaps	xmm1,xmm6	/* copy c3 */\
			__asm	subpd	xmm6,xmm3	/* c9  */\
			__asm	mulpd	xmm6,[ebx+0x80]	/* c9 = b8*c9 */\
			__asm	mulpd	xmm3,[ebx+0x50]	/* c6 = b5*c6 */\
		/*	c10= c3+c6-5*t5;			s10= s3+s6-5*t6 */\
			__asm	subpd	xmm2,xmm1	/* -c10 */\
			__asm	mulpd	xmm1,[ebx+0x20]	/* c3 = b2*c3 */\
		/*\
			c3 = b2*c3;					s3 = b2*s3;	\
			c6 = b5*c6;					s6 = b5*s6;	\
			c9 = b8*c9;					s9 = b8*s9;	\
			c10= (-b9)*(-c10);			s10= (-b9)*(-c10);	\
		*/\
			__asm	mulpd	xmm2,[ebx+0x90]	/* c10 = b9*c10*/\
		\
		/* a[0,3,6]*c[1,4,7] in ecx+__o[1,2,3]; a[1,4,7]*c[2,5,8] in xmm[7,4,5]; c[3,6,9] in xmm[1,3,6]; c10 in xmm2; xmm0 free */\
			__asm	addpd	xmm7,xmm1		/* c2 = b1*c2+c3 */\
			__asm	addpd	xmm1,[ecx+__o6]	/* c1 = b0*c1+c3 */\
			__asm	addpd	xmm4,xmm3		/* c5 = b4*c5+c6 */\
			__asm	addpd	xmm3,[ecx+__o7]	/* c4 = b3*c4+c6 */\
			__asm	addpd	xmm5,xmm6		/* c8 = b7*c8+c9 */\
			__asm	addpd	xmm6,[ecx+__o8]	/* c7 = b6*c7+c9 */\
		/*\
			sr1 = c10+c1-c7;			si1 = s10+s1-s7;\
			sr2 = c1+c2+c4+c5-c10;		si2 = s1+s2+s4+s5-s10;\
			sr3 = c10+c4+c7;			si3 = s10+s4+s7;\
			sr4 = c10+c5+c8;			si4 = s10+s5+s8;\
			sr5 = c10+c2-c8;			si5 = s10+s2-s8;\
		*/\
			__asm	xorpd	xmm0,xmm0	/* 0.0 */\
			__asm	subpd	xmm0,xmm2	/* copy (-c10) */\
			__asm	addpd	xmm0,xmm1	/* c1-c10 */\
			__asm	addpd	xmm1,xmm2	/* c10+c1 */\
			__asm	addpd	xmm0,xmm7	/* c1+c2-c10 */\
			__asm	addpd	xmm7,xmm2	/* c10+c2 */\
			__asm	subpd	xmm1,xmm6	/* sr1 = c10+c1-c7 */\
			__asm	addpd	xmm0,xmm3	/* c1+c2+c4-c10 */\
		/*	__asm	movaps	[ecx+__o6],xmm1	// store sr1 */\
			__asm	addpd	xmm3,xmm2	/* c10+c4 */\
			__asm	subpd	xmm7,xmm5	/* sr5 = c10+c2-c8 */\
			__asm	addpd	xmm0,xmm4	/* sr2 = c1+c2+c4+c5-c10 */\
		/*	__asm	movaps	[ecx+__oA],xmm7	// store sr5 */\
		/*	__asm	movaps	[ecx+__o7],xmm0	// store sr2 */\
			__asm	addpd	xmm4,xmm2	/* c10+c5 */\
			__asm	addpd	xmm3,xmm6	/* sr3 = c10+c4+c7 */\
			__asm	addpd	xmm4,xmm5	/* sr4 = c10+c5+c8 */\
		/*	__asm	movaps	[ecx+__o8],xmm3	// store sr3 */\
		/*	__asm	movaps	[ecx+__o9],xmm4	// store sr4 */\
		/*\
			si1,2,3,4,5 in xmm1,0,3,4,7:

			B1r = cr1 - si1;	// X1 = C1 + I*S1	\
			B2r = cr2 - si2;	// X2 = C2 + I*S2	\
			B3r = cr3 - si3;	// X3 = C3 + I*S3	\
			B4r = cr4 - si4;	// X4 = C4 + I*S4	\
			B5r = cr5 - si5;	// X5 = C5 + I*S5	\
			B6r = cr5 + si5;	// X6 =	C5 - I*S5	\
			B7r = cr4 + si4;	// X7 =	C4 - I*S4	\
			B8r = cr3 + si3;	// X8 =	C3 - I*S3	\
			B9r = cr2 + si2;	// X9 =	C2 - I*S2	\
			Bar = cr1 + si1;	// X10=	C1 - I*S1	\
		*/\
		/* si-terms get combined with cr's and output in Br's, so no need to decrement ecx: */\
		/* xmm2,5,6 FREE: */\
			__asm	movaps	xmm2,xmm1	/* cpy si1 */\
			__asm	addpd	xmm1,[ecx+__o1]	/* Bar */\
			__asm	addpd	xmm2,xmm2	/* 2 x si1 */\
			__asm	movaps	[ecx+__oA],xmm1	/* store Bar */\
			__asm	subpd	xmm1,xmm2		/* B1r */\
			__asm	movaps	[ecx+__o1],xmm1	/* store B1r */\
		\
			__asm	movaps	xmm5,xmm0	/* cpy si2 */\
			__asm	addpd	xmm0,[ecx+__o2]	/* B9r */\
			__asm	addpd	xmm5,xmm5	/* 2 x si2 */\
			__asm	movaps	[ecx+__o9],xmm0	/* store B9r */\
			__asm	subpd	xmm0,xmm5		/* B2r */\
			__asm	movaps	[ecx+__o2],xmm0	/* store B2r */\
		\
			__asm	movaps	xmm6,xmm3	/* cpy si3 */\
			__asm	addpd	xmm3,[ecx+__o3]	/* B8r */\
			__asm	addpd	xmm6,xmm6	/* 2 x si3 */\
			__asm	movaps	[ecx+__o8],xmm3	/* store B8r */\
			__asm	subpd	xmm3,xmm6		/* B3r */\
			__asm	movaps	[ecx+__o3],xmm3	/* store B3r */\
		\
			__asm	movaps	xmm2,xmm4	/* cpy si4 */\
			__asm	addpd	xmm4,[ecx+__o4]	/* B7r */\
			__asm	addpd	xmm2,xmm2	/* 2 x si4 */\
			__asm	movaps	[ecx+__o7],xmm4	/* store B7r */\
			__asm	subpd	xmm4,xmm2		/* B4r */\
			__asm	movaps	[ecx+__o4],xmm4	/* store B4r */\
		\
			__asm	movaps	xmm5,xmm7	/* cpy si5 */\
			__asm	addpd	xmm7,[ecx+__o5]	/* B6r */\
			__asm	addpd	xmm5,xmm5	/* 2 x si5 */\
			__asm	movaps	[ecx+__o6],xmm7	/* store B6r */\
			__asm	subpd	xmm7,xmm5		/* B5r */\
			__asm	movaps	[ecx+__o5],xmm7	/* store B5r */\
		}

	#else	/* GCC-style inline ASM: */

	  #if OS_BITS == 32

		#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
		{\
		__asm__ volatile (\
			"/********************************************/\n\t"\
			"/*       Here are the 5 cosine terms:       */\n\t"\
			"/********************************************/\n\t"\
			"movl	%[__I0],%%eax			\n\t"\
			"movl	%[__cc],%%ebx			\n\t"\
			"movl	%[__O0],%%ecx			\n\t"\
			"movaps	%c[__i1](%%eax),%%xmm1	\n\t"\
			"movaps	%c[__iA](%%eax),%%xmm5	\n\t"\
			"movaps	%c[__i2](%%eax),%%xmm2	\n\t"\
			"movaps	%c[__i9](%%eax),%%xmm6	\n\t"\
			"movaps	%c[__i3](%%eax),%%xmm3	\n\t"\
			"movaps	%c[__i8](%%eax),%%xmm7	\n\t"\
			"movaps	%c[__i4](%%eax),%%xmm4	\n\t"\
			"movaps	%c[__i7](%%eax),%%xmm0	\n\t"\
			"addl	$0x10,%%ecx				\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t"\
			"subpd	%%xmm0,%%xmm4			\n\t"\
			"movaps	%%xmm1,%c[__oA](%%ecx)	\n\t"\
			"movaps	%%xmm2,%c[__o9](%%ecx)	\n\t"\
			"movaps	%%xmm3,%c[__o8](%%ecx)	\n\t"\
			"movaps	%%xmm4,%c[__o7](%%ecx)	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t"\
			"movaps	%c[__i5](%%eax),%%xmm5	\n\t"\
			"movaps	%c[__i6](%%eax),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"subpd	%%xmm6,%%xmm5			\n\t"\
			"movaps	%%xmm5,%c[__o6](%%ecx)	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	(%%eax),%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"/********************************************/\n\t"\
			"/*               Real Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"subl	$0x10,%%ecx				\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"subpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%ebx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o1](%%ecx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0x10(%%ebx),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%ebx),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o3](%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%ebx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%ebx),%%xmm4		\n\t"\
			"mulpd	 0x70(%%ebx),%%xmm5		\n\t"\
			"mulpd	 0x10(%%ebx),%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%ebx),%%xmm6		\n\t"\
			"mulpd	 0x50(%%ebx),%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%ebx),%%xmm1		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	 0x90(%%ebx),%%xmm2		\n\t"\
			"movaps	%%xmm0,(%%ecx)			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%c[__o1](%%ecx),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	%c[__o3](%%ecx),%%xmm6	\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__o1](%%ecx)	\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"movaps	%%xmm7,%c[__o5](%%ecx)	\n\t"\
			"movaps	%%xmm2,%c[__o2](%%ecx)	\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%ecx)	\n\t"\
			"movaps	%%xmm4,%c[__o4](%%ecx)	\n\t"\
			"/********************************************/\n\t"\
			"/*          Imaginary Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"addl	$0x10,%%eax				\n\t"\
			"movaps	%c[__i1](%%eax),%%xmm1	\n\t"\
			"movaps	%c[__iA](%%eax),%%xmm5	\n\t"\
			"movaps	%c[__i2](%%eax),%%xmm2	\n\t"\
			"movaps	%c[__i9](%%eax),%%xmm6	\n\t"\
			"movaps	%c[__i3](%%eax),%%xmm3	\n\t"\
			"movaps	%c[__i8](%%eax),%%xmm7	\n\t"\
			"movaps	%c[__i4](%%eax),%%xmm4	\n\t"\
			"movaps	%c[__i7](%%eax),%%xmm0	\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t"\
			"subpd	%%xmm0,%%xmm4			\n\t"\
			"movaps	%%xmm1,%c[__oA](%%ecx)	\n\t"\
			"movaps	%%xmm2,%c[__o9](%%ecx)	\n\t"\
			"movaps	%%xmm3,%c[__o8](%%ecx)	\n\t"\
			"movaps	%%xmm4,%c[__o7](%%ecx)	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t"\
			"movaps	%c[__i5](%%eax),%%xmm5	\n\t"\
			"movaps	%c[__i6](%%eax),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"subpd	%%xmm6,%%xmm5			\n\t"\
			"movaps	%%xmm5,%c[__o6](%%ecx)	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	(%%eax),%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addl	$0x10,%%ecx				\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"subpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%ebx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o1](%%ecx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0x10(%%ebx),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%ebx),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o3](%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%ebx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%ebx),%%xmm4		\n\t"\
			"mulpd	 0x70(%%ebx),%%xmm5		\n\t"\
			"mulpd	 0x10(%%ebx),%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%ebx),%%xmm6		\n\t"\
			"mulpd	 0x50(%%ebx),%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%ebx),%%xmm1		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	 0x90(%%ebx),%%xmm2		\n\t"\
			"movaps	%%xmm0,(%%ecx)			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%c[__o1](%%ecx),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	%c[__o3](%%ecx),%%xmm6	\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__o1](%%ecx)	\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"movaps	%%xmm7,%c[__o5](%%ecx)	\n\t"\
			"movaps	%%xmm2,%c[__o2](%%ecx)	\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%ecx)	\n\t"\
			"movaps	%%xmm4,%c[__o4](%%ecx)	\n\t"\
			"/************************************************************************************************************************************/\n\t"\
			"/* Here are the 5 sine terms: Similar sequence to cosine terms, but with t21,19,17,15,13 replacing t3,5,7,9,11 and some sign flips: */\n\t"\
			"/************************************************************************************************************************************/\n\t"\
			"subl	$0x10,%%eax				\n\t"\
			"addl	$0xa0,%%ebx				\n\t"\
			"/********************************************/\n\t"\
			"/*               Real Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"movaps	%c[__oA](%%ecx),%%xmm1	\n\t"\
			"movaps	%c[__o9](%%ecx),%%xmm2	\n\t"\
			"movaps	%c[__o8](%%ecx),%%xmm3	\n\t"\
			"movaps	%c[__o7](%%ecx),%%xmm4	\n\t"\
			"movaps	%c[__o6](%%ecx),%%xmm5	\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%ebx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o6](%%ecx)	\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0xb0(%%ebx),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%ebx),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o8](%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%ebx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o7](%%ecx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%ebx),%%xmm4		\n\t"\
			"mulpd	 0x70(%%ebx),%%xmm5		\n\t"\
			"mulpd	 0x10(%%ebx),%%xmm7		\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%ebx),%%xmm6		\n\t"\
			"mulpd	 0x50(%%ebx),%%xmm3		\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%ebx),%%xmm1		\n\t"\
			"mulpd	 0x90(%%ebx),%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%c[__o6](%%ecx),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	%c[__o7](%%ecx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	%c[__o8](%%ecx),%%xmm6	\n\t"\
			"xorpd	%%xmm0,%%xmm0			\n\t"\
			"subpd	%%xmm2,%%xmm0			\n\t"\
			"addpd	%%xmm1,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm7,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"addpd	%%xmm3,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm1,%%xmm2			\n\t"\
			"addpd	%c[__o1](%%ecx),%%xmm1	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__o1](%%ecx)	\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"movaps	%%xmm1,%c[__oA](%%ecx)	\n\t"\
			"movaps	%%xmm0,%%xmm5			\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm0	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm0,%c[__o2](%%ecx)	\n\t"\
			"subpd	%%xmm5,%%xmm0			\n\t"\
			"movaps	%%xmm0,%c[__o9](%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm6			\n\t"\
			"addpd	%c[__o3](%%ecx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%ecx)	\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm3,%c[__o8](%%ecx)	\n\t"\
			"movaps	%%xmm4,%%xmm2			\n\t"\
			"addpd	%c[__o4](%%ecx),%%xmm4	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm4,%c[__o4](%%ecx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"movaps	%%xmm4,%c[__o7](%%ecx)	\n\t"\
			"movaps	%%xmm7,%%xmm5			\n\t"\
			"addpd	%c[__o5](%%ecx),%%xmm7	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm7,%c[__o5](%%ecx)	\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"movaps	%%xmm7,%c[__o6](%%ecx)	\n\t"\
			"/********************************************/\n\t"\
			"/*          Imaginary Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"subl	$0x10,%%ecx				\n\t"\
			"movaps	%c[__oA](%%ecx),%%xmm1	\n\t"\
			"movaps	%c[__o9](%%ecx),%%xmm2	\n\t"\
			"movaps	%c[__o8](%%ecx),%%xmm3	\n\t"\
			"movaps	%c[__o7](%%ecx),%%xmm4	\n\t"\
			"movaps	%c[__o6](%%ecx),%%xmm5	\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%ebx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o6](%%ecx)	\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0xb0(%%ebx),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%ebx),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o8](%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%ebx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o7](%%ecx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%ebx),%%xmm4		\n\t"\
			"mulpd	 0x70(%%ebx),%%xmm5		\n\t"\
			"mulpd	 0x10(%%ebx),%%xmm7		\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%ebx),%%xmm6		\n\t"\
			"mulpd	 0x50(%%ebx),%%xmm3		\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%ebx),%%xmm1		\n\t"\
			"mulpd	 0x90(%%ebx),%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%c[__o6](%%ecx),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	%c[__o7](%%ecx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	%c[__o8](%%ecx),%%xmm6	\n\t"\
			"xorpd	%%xmm0,%%xmm0			\n\t"\
			"subpd	%%xmm2,%%xmm0			\n\t"\
			"addpd	%%xmm1,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm7,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"addpd	%%xmm3,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm1,%%xmm2			\n\t"\
			"addpd	%c[__o1](%%ecx),%%xmm1	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__oA](%%ecx)	\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"movaps	%%xmm1,%c[__o1](%%ecx)	\n\t"\
			"movaps	%%xmm0,%%xmm5			\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm0	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm0,%c[__o9](%%ecx)	\n\t"\
			"subpd	%%xmm5,%%xmm0			\n\t"\
			"movaps	%%xmm0,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm6			\n\t"\
			"addpd	%c[__o3](%%ecx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	%%xmm3,%c[__o8](%%ecx)	\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%ecx)	\n\t"\
			"movaps	%%xmm4,%%xmm2			\n\t"\
			"addpd	%c[__o4](%%ecx),%%xmm4	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm4,%c[__o7](%%ecx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"movaps	%%xmm4,%c[__o4](%%ecx)	\n\t"\
			"movaps	%%xmm7,%%xmm5			\n\t"\
			"addpd	%c[__o5](%%ecx),%%xmm7	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm7,%c[__o6](%%ecx)	\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"movaps	%%xmm7,%c[__o5](%%ecx)	\n\t"\
			:					/* outputs: none */\
			: [__I0] "m" (XI0)	/* All inputs from memory addresses here */\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__i3] "e" (Xi3)\
			 ,[__i4] "e" (Xi4)\
			 ,[__i5] "e" (Xi5)\
			 ,[__i6] "e" (Xi6)\
			 ,[__i7] "e" (Xi7)\
			 ,[__i8] "e" (Xi8)\
			 ,[__i9] "e" (Xi9)\
			 ,[__iA] "e" (XiA)\
			 ,[__cc] "m" (Xcc)\
			 ,[__O0] "m" (XO0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			 ,[__o3] "e" (Xo3)\
			 ,[__o4] "e" (Xo4)\
			 ,[__o5] "e" (Xo5)\
			 ,[__o6] "e" (Xo6)\
			 ,[__o7] "e" (Xo7)\
			 ,[__o8] "e" (Xo8)\
			 ,[__o9] "e" (Xo9)\
			 ,[__oA] "e" (XoA)\
			: "eax","ebx","ecx"		/* Clobbered registers */\
			);\
		}

	  #else

		#define GCC_ASM_FULL_INLINE	0	// 0 to use small-macro form below, 1 to inline the fused macros as single big blob of asm (64-bit only)

		#if GCC_ASM_FULL_INLINE

			#include "radix44_ditN_cy_dif1_gcc64.h"

		#endif	// GCC_ASM_FULL_INLINE

		#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
		{\
		__asm__ volatile (\
			"/********************************************/\n\t"\
			"/*       Here are the 5 cosine terms:       */\n\t"\
			"/********************************************/\n\t"\
			"movq	%[__I0],%%rax			\n\t"\
			"movq	%[__cc],%%rbx			\n\t"\
			"movq	%[__O0],%%rcx			\n\t"\
			"movaps	%c[__i1](%%rax),%%xmm1	\n\t"\
			"movaps	%c[__iA](%%rax),%%xmm5	\n\t"\
			"movaps	%c[__i2](%%rax),%%xmm2	\n\t"\
			"movaps	%c[__i9](%%rax),%%xmm6	\n\t"\
			"movaps	%c[__i3](%%rax),%%xmm3	\n\t"\
			"movaps	%c[__i8](%%rax),%%xmm7	\n\t"\
			"movaps	%c[__i4](%%rax),%%xmm4	\n\t"\
			"movaps	%c[__i7](%%rax),%%xmm0	\n\t"\
			"addq	$0x10,%%rcx				\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t"\
			"subpd	%%xmm0,%%xmm4			\n\t"\
			"movaps	%%xmm1,%c[__oA](%%rcx)	\n\t"\
			"movaps	%%xmm2,%c[__o9](%%rcx)	\n\t"\
			"movaps	%%xmm3,%c[__o8](%%rcx)	\n\t"\
			"movaps	%%xmm4,%c[__o7](%%rcx)	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t"\
			"movaps	%c[__i5](%%rax),%%xmm5	\n\t"\
			"movaps	%c[__i6](%%rax),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"subpd	%%xmm6,%%xmm5			\n\t"\
			"movaps	%%xmm5,%c[__o6](%%rcx)	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	(%%rax),%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"/********************************************/\n\t"\
			"/*               Real Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"subq	$0x10,%%rcx				\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"subpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%rbx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o1](%%rcx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0x10(%%rbx),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%rbx),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o3](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%rbx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o2](%%rcx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%rbx),%%xmm4		\n\t"\
			"mulpd	 0x70(%%rbx),%%xmm5		\n\t"\
			"mulpd	 0x10(%%rbx),%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%rbx),%%xmm6		\n\t"\
			"mulpd	 0x50(%%rbx),%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%rbx),%%xmm1		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	 0x90(%%rbx),%%xmm2		\n\t"\
			"movaps	%%xmm0,(%%rcx)			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%c[__o1](%%rcx),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	%c[__o3](%%rcx),%%xmm6	\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__o1](%%rcx)	\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"movaps	%%xmm7,%c[__o5](%%rcx)	\n\t"\
			"movaps	%%xmm2,%c[__o2](%%rcx)	\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%rcx)	\n\t"\
			"movaps	%%xmm4,%c[__o4](%%rcx)	\n\t"\
			"/********************************************/\n\t"\
			"/*          Imaginary Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"addq	$0x10,%%rax				\n\t"\
			"movaps	%c[__i1](%%rax),%%xmm1	\n\t"\
			"movaps	%c[__iA](%%rax),%%xmm5	\n\t"\
			"movaps	%c[__i2](%%rax),%%xmm2	\n\t"\
			"movaps	%c[__i9](%%rax),%%xmm6	\n\t"\
			"movaps	%c[__i3](%%rax),%%xmm3	\n\t"\
			"movaps	%c[__i8](%%rax),%%xmm7	\n\t"\
			"movaps	%c[__i4](%%rax),%%xmm4	\n\t"\
			"movaps	%c[__i7](%%rax),%%xmm0	\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t"\
			"subpd	%%xmm0,%%xmm4			\n\t"\
			"movaps	%%xmm1,%c[__oA](%%rcx)	\n\t"\
			"movaps	%%xmm2,%c[__o9](%%rcx)	\n\t"\
			"movaps	%%xmm3,%c[__o8](%%rcx)	\n\t"\
			"movaps	%%xmm4,%c[__o7](%%rcx)	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t"\
			"movaps	%c[__i5](%%rax),%%xmm5	\n\t"\
			"movaps	%c[__i6](%%rax),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"subpd	%%xmm6,%%xmm5			\n\t"\
			"movaps	%%xmm5,%c[__o6](%%rcx)	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	(%%rax),%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addq	$0x10,%%rcx				\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"subpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%rbx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o1](%%rcx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0x10(%%rbx),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%rbx),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o3](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%rbx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o2](%%rcx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%rbx),%%xmm4		\n\t"\
			"mulpd	 0x70(%%rbx),%%xmm5		\n\t"\
			"mulpd	 0x10(%%rbx),%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%rbx),%%xmm6		\n\t"\
			"mulpd	 0x50(%%rbx),%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%rbx),%%xmm1		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	 0x90(%%rbx),%%xmm2		\n\t"\
			"movaps	%%xmm0,(%%rcx)			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%c[__o1](%%rcx),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	%c[__o3](%%rcx),%%xmm6	\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__o1](%%rcx)	\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"movaps	%%xmm7,%c[__o5](%%rcx)	\n\t"\
			"movaps	%%xmm2,%c[__o2](%%rcx)	\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%rcx)	\n\t"\
			"movaps	%%xmm4,%c[__o4](%%rcx)	\n\t"\
			"/************************************************************************************************************************************/\n\t"\
			"/* Here are the 5 sine terms: Similar sequence to cosine terms, but with t21,19,17,15,13 replacing t3,5,7,9,11 and some sign flips: */\n\t"\
			"/************************************************************************************************************************************/\n\t"\
			"subq	$0x10,%%rax				\n\t"\
			"addq	$0xa0,%%rbx				\n\t"\
			"/********************************************/\n\t"\
			"/*               Real Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"movaps	%c[__oA](%%rcx),%%xmm1	\n\t"\
			"movaps	%c[__o9](%%rcx),%%xmm2	\n\t"\
			"movaps	%c[__o8](%%rcx),%%xmm3	\n\t"\
			"movaps	%c[__o7](%%rcx),%%xmm4	\n\t"\
			"movaps	%c[__o6](%%rcx),%%xmm5	\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%rbx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o6](%%rcx)	\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0xb0(%%rbx),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%rbx),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o8](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%rbx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o7](%%rcx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%rbx),%%xmm4		\n\t"\
			"mulpd	 0x70(%%rbx),%%xmm5		\n\t"\
			"mulpd	 0x10(%%rbx),%%xmm7		\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%rbx),%%xmm6		\n\t"\
			"mulpd	 0x50(%%rbx),%%xmm3		\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%rbx),%%xmm1		\n\t"\
			"mulpd	 0x90(%%rbx),%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%c[__o6](%%rcx),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	%c[__o7](%%rcx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	%c[__o8](%%rcx),%%xmm6	\n\t"\
			"xorpd	%%xmm0,%%xmm0			\n\t"\
			"subpd	%%xmm2,%%xmm0			\n\t"\
			"addpd	%%xmm1,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm7,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"addpd	%%xmm3,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm1,%%xmm2			\n\t"\
			"addpd	%c[__o1](%%rcx),%%xmm1	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__o1](%%rcx)	\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"movaps	%%xmm1,%c[__oA](%%rcx)	\n\t"\
			"movaps	%%xmm0,%%xmm5			\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm0	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm0,%c[__o2](%%rcx)	\n\t"\
			"subpd	%%xmm5,%%xmm0			\n\t"\
			"movaps	%%xmm0,%c[__o9](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm6			\n\t"\
			"addpd	%c[__o3](%%rcx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%rcx)	\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm3,%c[__o8](%%rcx)	\n\t"\
			"movaps	%%xmm4,%%xmm2			\n\t"\
			"addpd	%c[__o4](%%rcx),%%xmm4	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm4,%c[__o4](%%rcx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"movaps	%%xmm4,%c[__o7](%%rcx)	\n\t"\
			"movaps	%%xmm7,%%xmm5			\n\t"\
			"addpd	%c[__o5](%%rcx),%%xmm7	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm7,%c[__o5](%%rcx)	\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"movaps	%%xmm7,%c[__o6](%%rcx)	\n\t"\
			"/********************************************/\n\t"\
			"/*          Imaginary Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"subq	$0x10,%%rcx				\n\t"\
			"movaps	%c[__oA](%%rcx),%%xmm1	\n\t"\
			"movaps	%c[__o9](%%rcx),%%xmm2	\n\t"\
			"movaps	%c[__o8](%%rcx),%%xmm3	\n\t"\
			"movaps	%c[__o7](%%rcx),%%xmm4	\n\t"\
			"movaps	%c[__o6](%%rcx),%%xmm5	\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%rbx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o6](%%rcx)	\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0xb0(%%rbx),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%rbx),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o8](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%rbx),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o7](%%rcx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%rbx),%%xmm4		\n\t"\
			"mulpd	 0x70(%%rbx),%%xmm5		\n\t"\
			"mulpd	 0x10(%%rbx),%%xmm7		\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%rbx),%%xmm6		\n\t"\
			"mulpd	 0x50(%%rbx),%%xmm3		\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%rbx),%%xmm1		\n\t"\
			"mulpd	 0x90(%%rbx),%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%c[__o6](%%rcx),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	%c[__o7](%%rcx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	%c[__o8](%%rcx),%%xmm6	\n\t"\
			"xorpd	%%xmm0,%%xmm0			\n\t"\
			"subpd	%%xmm2,%%xmm0			\n\t"\
			"addpd	%%xmm1,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm7,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"addpd	%%xmm3,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm1,%%xmm2			\n\t"\
			"addpd	%c[__o1](%%rcx),%%xmm1	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__oA](%%rcx)	\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"movaps	%%xmm1,%c[__o1](%%rcx)	\n\t"\
			"movaps	%%xmm0,%%xmm5			\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm0	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm0,%c[__o9](%%rcx)	\n\t"\
			"subpd	%%xmm5,%%xmm0			\n\t"\
			"movaps	%%xmm0,%c[__o2](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm6			\n\t"\
			"addpd	%c[__o3](%%rcx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	%%xmm3,%c[__o8](%%rcx)	\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%rcx)	\n\t"\
			"movaps	%%xmm4,%%xmm2			\n\t"\
			"addpd	%c[__o4](%%rcx),%%xmm4	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm4,%c[__o7](%%rcx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"movaps	%%xmm4,%c[__o4](%%rcx)	\n\t"\
			"movaps	%%xmm7,%%xmm5			\n\t"\
			"addpd	%c[__o5](%%rcx),%%xmm7	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm7,%c[__o6](%%rcx)	\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"movaps	%%xmm7,%c[__o5](%%rcx)	\n\t"\
			:					/* outputs: none */\
			: [__I0] "m" (XI0)	/* All inputs from memory addresses here */\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__i3] "e" (Xi3)\
			 ,[__i4] "e" (Xi4)\
			 ,[__i5] "e" (Xi5)\
			 ,[__i6] "e" (Xi6)\
			 ,[__i7] "e" (Xi7)\
			 ,[__i8] "e" (Xi8)\
			 ,[__i9] "e" (Xi9)\
			 ,[__iA] "e" (XiA)\
			 ,[__cc] "m" (Xcc)\
			 ,[__O0] "m" (XO0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			 ,[__o3] "e" (Xo3)\
			 ,[__o4] "e" (Xo4)\
			 ,[__o5] "e" (Xo5)\
			 ,[__o6] "e" (Xo6)\
			 ,[__o7] "e" (Xo7)\
			 ,[__o8] "e" (Xo8)\
			 ,[__o9] "e" (Xo9)\
			 ,[__oA] "e" (XoA)\
			: "rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
			);\
		}

	  #endif

	#endif

#endif

/**************/

int radix44_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-44 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-44 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const uint32 RADIX = 44;
	const double crnd = 3.0*0x4000000*0x2000000;
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40;
	/*...Fast length-5 cyclic convolution scheme needs the following: */
	static double 	a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1*/

			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
	static double radix_inv, n2inv;
	double scale, maxerr = 0.0;
	double a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p36r,a1p37r,a1p38r,a1p39r,a1p40r,a1p41r,a1p42r,a1p43r;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static int cslots_in_local_store;
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD

	#ifdef USE_PTHREAD
		static struct complex *__r0;	// Base address for discrete per-thread local stores
		static struct cy_thread_data_t *tdat = 0x0;
		// Threadpool-based dispatch stuff:
		static int main_work_units = 0, pool_work_units = 0;
		static struct threadpool *tpool = 0x0;
		static int task_is_blocking = TRUE;
		static thread_control_t thread_control = {0,0,0};
		// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
		static task_control_t   task_control = {NULL, (void*)cy44_process_chunk, NULL, 0x0};
	#endif

  #else
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
  #endif

	static struct complex *two,*five, *ua0,*ua1,*ua2,*ua3,*ua4,*ua5,*ua6,*ua7,*ua8,*ua9, *ub0,*ub1,*ub2,*ub3,*ub4,*ub5,*ub6,*ub7,*ub8,*ub9, *max_err, *sse2_rnd, *half_arr, *tmp,*tmp2;
	static struct complex
	 *t00r,*t01r,*t02r,*t03r,*t04r,*t05r,*t06r,*t07r,*t08r,*t09r,*t0ar,*t00i,*t01i,*t02i,*t03i,*t04i,*t05i,*t06i,*t07i,*t08i,*t09i,*t0ai
	,*t10r,*t11r,*t12r,*t13r,*t14r,*t15r,*t16r,*t17r,*t18r,*t19r,*t1ar,*t10i,*t11i,*t12i,*t13i,*t14i,*t15i,*t16i,*t17i,*t18i,*t19i,*t1ai
	,*t20r,*t21r,*t22r,*t23r,*t24r,*t25r,*t26r,*t27r,*t28r,*t29r,*t2ar,*t20i,*t21i,*t22i,*t23i,*t24i,*t25i,*t26i,*t27i,*t28i,*t29i,*t2ai
	,*t30r,*t31r,*t32r,*t33r,*t34r,*t35r,*t36r,*t37r,*t38r,*t39r,*t3ar,*t30i,*t31i,*t32i,*t33i,*t34i,*t35i,*t36i,*t37i,*t38i,*t39i,*t3ai
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p40r,*s1p41r,*s1p42r,*s1p43r
	,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i,*s1p20i,*s1p21i,*s1p22i,*s1p23i,*s1p24i,*s1p25i,*s1p26i,*s1p27i,*s1p28i,*s1p29i,*s1p30i,*s1p31i,*s1p32i,*s1p33i,*s1p34i,*s1p35i,*s1p36i,*s1p37i,*s1p38i,*s1p39i,*s1p40i,*s1p41i,*s1p42i,*s1p43i;
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39,*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43;
	static struct complex *cy00,*cy02,*cy04,*cy06,*cy08,*cy10,*cy12,*cy14,*cy16,*cy18,*cy20,*cy22,*cy24,*cy26,*cy28,*cy30,*cy32,*cy34,*cy36,*cy38,*cy40,*cy42;

#else

	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double rt,it,temp,frac;
	double wt,wtinv,wtA,wtB,wtC;
	int m,m2;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21
		,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,bjmodn40,bjmodn41,bjmodn42,bjmodn43;
	double a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,a1p36i,a1p37i,a1p38i,a1p39i,a1p40i,a1p41i,a1p42i,a1p43i;
	double cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23,cy24,cy25,cy26,cy27,cy28,cy29,cy30,cy31,cy32,cy33,cy34,cy35,cy36,cy37,cy38,cy39,cy40,cy41,cy42,cy43
		,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai
		,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai
		,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai
		,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0,*_bjmodn28 = 0x0,*_bjmodn29 = 0x0,*_bjmodn30 = 0x0,*_bjmodn31 = 0x0,*_bjmodn32 = 0x0,*_bjmodn33 = 0x0,*_bjmodn34 = 0x0,*_bjmodn35 = 0x0,*_bjmodn36 = 0x0,*_bjmodn37 = 0x0,*_bjmodn38 = 0x0,*_bjmodn39 = 0x0,*_bjmodn40 = 0x0,*_bjmodn41 = 0x0,*_bjmodn42 = 0x0,*_bjmodn43 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy00 = 0x0,*_cy01 = 0x0,*_cy02 = 0x0,*_cy03 = 0x0,*_cy04 = 0x0,*_cy05 = 0x0,*_cy06 = 0x0,*_cy07 = 0x0,*_cy08 = 0x0,*_cy09 = 0x0,*_cy10 = 0x0,*_cy11 = 0x0,*_cy12 = 0x0,*_cy13 = 0x0,*_cy14 = 0x0,*_cy15 = 0x0,*_cy16 = 0x0,*_cy17 = 0x0,*_cy18 = 0x0,*_cy19 = 0x0,*_cy20 = 0x0,*_cy21 = 0x0,*_cy22 = 0x0,*_cy23 = 0x0,*_cy24 = 0x0,*_cy25 = 0x0,*_cy26 = 0x0,*_cy27 = 0x0,*_cy28 = 0x0,*_cy29 = 0x0,*_cy30 = 0x0,*_cy31 = 0x0,*_cy32 = 0x0,*_cy33 = 0x0,*_cy34 = 0x0,*_cy35 = 0x0,*_cy36 = 0x0,*_cy37 = 0x0,*_cy38 = 0x0,*_cy39 = 0x0,*_cy40 = 0x0,*_cy41 = 0x0,*_cy42 = 0x0,*_cy43 = 0x0;

/*...change n44 and n_div_wt to non-static to work around a gcc compiler bug. */
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "radix24_ditN_cy_dif1: Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/RADIX in radix24_ditN_cy_dif1.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the smallest
		power of 2 that is >= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else
		{
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> (i-1);
		}

		if(CY_THREADS > MAX_THREADS)
		{
		//	CY_THREADS = MAX_THREADS;
			fprintf(stderr,"WARN: CY_THREADS = %d exceeds number of cores = %d\n", CY_THREADS, MAX_THREADS);
		}
		ASSERT(HERE, CY_THREADS >= NTHREADS,"CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, NDIVR    %CY_THREADS == 0,"NDIVR    %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"n_div_nwt%CY_THREADS != 0");
		}

	  #ifdef USE_PTHREAD

		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#ifdef OS_TYPE_MACOSX

			if(CY_THREADS > 1) {
				main_work_units = CY_THREADS/2;
				pool_work_units = CY_THREADS - main_work_units;
				ASSERT(HERE, 0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
				printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
			} else {
				main_work_units = 1;
				printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
			}

		#else

			pool_work_units = CY_THREADS;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

		#endif

		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);

	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix44_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_COMPLEX(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix44_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 88x2 16-byte slots of sc_arr for temporaries, next 21 for the constants needed by the radix-11 DFT,
	next RADIX/2 = 22 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
								tmp = sc_ptr + 0x58;	tmp2 = tmp + 0x58;
		s1p00r = sc_ptr + 0x00;	t00r = tmp + 0x00;	cy00  = tmp2+ 0x00;
		s1p00i = sc_ptr + 0x01;	t00i = tmp + 0x01;	cy02  = tmp2+ 0x01;
		s1p01r = sc_ptr + 0x02;	t01r = tmp + 0x02;	cy04  = tmp2+ 0x02;
		s1p01i = sc_ptr + 0x03;	t01i = tmp + 0x03;	cy06  = tmp2+ 0x03;
		s1p02r = sc_ptr + 0x04;	t02r = tmp + 0x04;	cy08  = tmp2+ 0x04;
		s1p02i = sc_ptr + 0x05;	t02i = tmp + 0x05;	cy10  = tmp2+ 0x05;
		s1p03r = sc_ptr + 0x06;	t03r = tmp + 0x06;	cy12  = tmp2+ 0x06;
		s1p03i = sc_ptr + 0x07;	t03i = tmp + 0x07;	cy14  = tmp2+ 0x07;
		s1p04r = sc_ptr + 0x08;	t04r = tmp + 0x08;	cy16  = tmp2+ 0x08;
		s1p04i = sc_ptr + 0x09;	t04i = tmp + 0x09;	cy18  = tmp2+ 0x09;
		s1p05r = sc_ptr + 0x0a;	t05r = tmp + 0x0a;	cy20  = tmp2+ 0x0a;
		s1p05i = sc_ptr + 0x0b;	t05i = tmp + 0x0b;	cy22  = tmp2+ 0x0b;
		s1p06r = sc_ptr + 0x0c;	t06r = tmp + 0x0c;	cy24  = tmp2+ 0x0c;
		s1p06i = sc_ptr + 0x0d;	t06i = tmp + 0x0d;	cy26  = tmp2+ 0x0d;
		s1p07r = sc_ptr + 0x0e;	t07r = tmp + 0x0e;	cy28  = tmp2+ 0x0e;
		s1p07i = sc_ptr + 0x0f;	t07i = tmp + 0x0f;	cy30  = tmp2+ 0x0f;
		s1p08r = sc_ptr + 0x10;	t08r = tmp + 0x10;	cy32  = tmp2+ 0x10;
		s1p08i = sc_ptr + 0x11;	t08i = tmp + 0x11;	cy34  = tmp2+ 0x11;
		s1p09r = sc_ptr + 0x12;	t09r = tmp + 0x12;	cy36  = tmp2+ 0x12;
		s1p09i = sc_ptr + 0x13;	t09i = tmp + 0x13;	cy38  = tmp2+ 0x13;
		s1p10r = sc_ptr + 0x14;	t0ar = tmp + 0x14;	cy40  = tmp2+ 0x14;
		s1p10i = sc_ptr + 0x15;	t0ai = tmp + 0x15;	cy42  = tmp2+ 0x15;
		s1p11r = sc_ptr + 0x16;	t10r = tmp + 0x16;	two     = tmp2+ 0x16;
		s1p11i = sc_ptr + 0x17;	t10i = tmp + 0x17;	five    = tmp2+ 0x17;
		s1p12r = sc_ptr + 0x18;	t11r = tmp + 0x18;	ua0     = tmp2+ 0x18;
		s1p12i = sc_ptr + 0x19;	t11i = tmp + 0x19;	ua1     = tmp2+ 0x19;
		s1p13r = sc_ptr + 0x1a;	t12r = tmp + 0x1a;	ua2     = tmp2+ 0x1a;
		s1p13i = sc_ptr + 0x1b;	t12i = tmp + 0x1b;	ua3     = tmp2+ 0x1b;
		s1p14r = sc_ptr + 0x1c;	t13r = tmp + 0x1c;	ua4     = tmp2+ 0x1c;
		s1p14i = sc_ptr + 0x1d;	t13i = tmp + 0x1d;	ua5     = tmp2+ 0x1d;
		s1p15r = sc_ptr + 0x1e;	t14r = tmp + 0x1e;	ua6     = tmp2+ 0x1e;
		s1p15i = sc_ptr + 0x1f;	t14i = tmp + 0x1f;	ua7     = tmp2+ 0x1f;
		s1p16r = sc_ptr + 0x20;	t15r = tmp + 0x20;	ua8     = tmp2+ 0x20;
		s1p16i = sc_ptr + 0x21;	t15i = tmp + 0x21;	ua9     = tmp2+ 0x21;
		s1p17r = sc_ptr + 0x22;	t16r = tmp + 0x22;	ub0     = tmp2+ 0x22;
		s1p17i = sc_ptr + 0x23;	t16i = tmp + 0x23;	ub1     = tmp2+ 0x23;
		s1p18r = sc_ptr + 0x24;	t17r = tmp + 0x24;	ub2     = tmp2+ 0x24;
		s1p18i = sc_ptr + 0x25;	t17i = tmp + 0x25;	ub3     = tmp2+ 0x25;
		s1p19r = sc_ptr + 0x26;	t18r = tmp + 0x26;	ub4     = tmp2+ 0x26;
		s1p19i = sc_ptr + 0x27;	t18i = tmp + 0x27;	ub5     = tmp2+ 0x27;
		s1p20r = sc_ptr + 0x28;	t19r = tmp + 0x28;	ub6     = tmp2+ 0x28;
		s1p20i = sc_ptr + 0x29;	t19i = tmp + 0x29;	ub7     = tmp2+ 0x29;
		s1p21r = sc_ptr + 0x2a;	t1ar = tmp + 0x2a;	ub8     = tmp2+ 0x2a;
		s1p21i = sc_ptr + 0x2b;	t1ai = tmp + 0x2b;	ub9     = tmp2+ 0x2b;
		s1p22r = sc_ptr + 0x2c;	t20r = tmp + 0x2c;	max_err = tmp2+ 0x2c;
		s1p22i = sc_ptr + 0x2d;	t20i = tmp + 0x2d;	sse2_rnd= tmp2+ 0x2d;
		s1p23r = sc_ptr + 0x2e;	t21r = tmp + 0x2e;	half_arr= tmp2+ 0x2e;	/* This table needs 20x16 bytes */
		s1p23i = sc_ptr + 0x2f;	t21i = tmp + 0x2f;
		s1p24r = sc_ptr + 0x30;	t22r = tmp + 0x30;
		s1p24i = sc_ptr + 0x31;	t22i = tmp + 0x31;
		s1p25r = sc_ptr + 0x32;	t23r = tmp + 0x32;
		s1p25i = sc_ptr + 0x33;	t23i = tmp + 0x33;
		s1p26r = sc_ptr + 0x34;	t24r = tmp + 0x34;
		s1p26i = sc_ptr + 0x35;	t24i = tmp + 0x35;
		s1p27r = sc_ptr + 0x36;	t25r = tmp + 0x36;
		s1p27i = sc_ptr + 0x37;	t25i = tmp + 0x37;
		s1p28r = sc_ptr + 0x38;	t26r = tmp + 0x38;
		s1p28i = sc_ptr + 0x39;	t26i = tmp + 0x39;
		s1p29r = sc_ptr + 0x3a;	t27r = tmp + 0x3a;
		s1p29i = sc_ptr + 0x3b;	t27i = tmp + 0x3b;
		s1p30r = sc_ptr + 0x3c;	t28r = tmp + 0x3c;
		s1p30i = sc_ptr + 0x3d;	t28i = tmp + 0x3d;
		s1p31r = sc_ptr + 0x3e;	t29r = tmp + 0x3e;
		s1p31i = sc_ptr + 0x3f;	t29i = tmp + 0x3f;
		s1p32r = sc_ptr + 0x40;	t2ar = tmp + 0x40;
		s1p32i = sc_ptr + 0x41;	t2ai = tmp + 0x41;
		s1p33r = sc_ptr + 0x42;	t30r = tmp + 0x42;
		s1p33i = sc_ptr + 0x43;	t30i = tmp + 0x43;
		s1p34r = sc_ptr + 0x44;	t31r = tmp + 0x44;
		s1p34i = sc_ptr + 0x45;	t31i = tmp + 0x45;
		s1p35r = sc_ptr + 0x46;	t32r = tmp + 0x46;
		s1p35i = sc_ptr + 0x47;	t32i = tmp + 0x47;
		s1p36r = sc_ptr + 0x48;	t33r = tmp + 0x48;
		s1p36i = sc_ptr + 0x49;	t33i = tmp + 0x49;
		s1p37r = sc_ptr + 0x4a;	t34r = tmp + 0x4a;
		s1p37i = sc_ptr + 0x4b;	t34i = tmp + 0x4b;
		s1p38r = sc_ptr + 0x4c;	t35r = tmp + 0x4c;
		s1p38i = sc_ptr + 0x4d;	t35i = tmp + 0x4d;
		s1p39r = sc_ptr + 0x4e;	t36r = tmp + 0x4e;
		s1p39i = sc_ptr + 0x4f;	t36i = tmp + 0x4f;
		s1p40r = sc_ptr + 0x50;	t37r = tmp + 0x50;
		s1p40i = sc_ptr + 0x51;	t37i = tmp + 0x51;
		s1p41r = sc_ptr + 0x52;	t38r = tmp + 0x52;
		s1p41i = sc_ptr + 0x53;	t38i = tmp + 0x53;
		s1p42r = sc_ptr + 0x54;	t39r = tmp + 0x54;
		s1p42i = sc_ptr + 0x55;	t39i = tmp + 0x55;
		s1p43r = sc_ptr + 0x56;	t3ar = tmp + 0x56;
		s1p43i = sc_ptr + 0x57;	t3ai = tmp + 0x57;
		ASSERT(HERE, (radix44_creals_in_local_store << 4) >= ((long)half_arr - (long)s1p00r) + (20 << 4), "radix44_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		two ->re = two ->im = 2.0;
		five->re = five->im = 5.0;
		ua0 ->re = ua0 ->im = a0;
		ua1 ->re = ua1 ->im = a1;
		ua2 ->re = ua2 ->im = a2;
		ua3 ->re = ua3 ->im = a3;
		ua4 ->re = ua4 ->im = a4;
		ua5 ->re = ua5 ->im = a5;
		ua6 ->re = ua6 ->im = a6;
		ua7 ->re = ua7 ->im = a7;
		ua8 ->re = ua8 ->im = a8;
		ua9 ->re = ua9 ->im = a9;
		ub0 ->re = ub0 ->im = b0;
		ub1 ->re = ub1 ->im = b1;
		ub2 ->re = ub2 ->im = b2;
		ub3 ->re = ub3 ->im = b3;
		ub4 ->re = ub4 ->im = b4;
		ub5 ->re = ub5 ->im = b5;
		ub6 ->re = ub6 ->im = b6;
		ub7 ->re = ub7 ->im = b7;
		ub8 ->re = ub8 ->im = b8;
		ub9 ->re = ub9 ->im =-b9;	/* Flip sign to simplify code re-use in radix-11 SSE2 macro */

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		sse2_rnd->re = sse2_rnd->im = crnd;

		/* SSE2 version of the one_half array - we have a 2-bit lookup, low bit is from the low word of the carry pair,
		high bit from the high, i.e. based on this lookup index [listed with LSB at right], we have:

			index	half_lo	half_hi
			00		1.0		1.0
			01		.50		1.0
			10		1.0		.50
			11		.50		.50

		The inverse-weights computation uses a similar table, but with all entries multiplied by .50:

			index2	half_lo	half_hi
			00		.50		.50
			01		.25		.50
			10		.50		.25
			11		.25		.25

		We do similarly for the base[] and baseinv[] table lookups - each of these get 4 further slots in half_arr.
		We also allocate a further 4 16-byte slots [uninitialized] for storage of the wtl,wtn,wtlp1,wtnm1 locals.
		*/
		tmp = half_arr;
		/* Forward-weight multipliers: */
		tmp->re = 1.0;	tmp->im = 1.0;	++tmp;
		tmp->re = .50;	tmp->im = 1.0;	++tmp;
		tmp->re = 1.0;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .50;	++tmp;
		/* Inverse-weight multipliers: */
		tmp->re = .50;	tmp->im = .50;	++tmp;
		tmp->re = .25;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .25;	++tmp;
		tmp->re = .25;	tmp->im = .25;	++tmp;
		/* Forward-base[] multipliers: */
		tmp->re = base   [0];	tmp->im = base   [0];	++tmp;
		tmp->re = base   [1];	tmp->im = base   [0];	++tmp;
		tmp->re = base   [0];	tmp->im = base   [1];	++tmp;
		tmp->re = base   [1];	tmp->im = base   [1];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->re = baseinv[0];	tmp->im = baseinv[0];	++tmp;
		tmp->re = baseinv[1];	tmp->im = baseinv[0];	++tmp;
		tmp->re = baseinv[0];	tmp->im = baseinv[1];	++tmp;
		tmp->re = baseinv[1];	tmp->im = baseinv[1];	++tmp;

		/* Floating-point sign mask used for FABS on packed doubles: */
		sign_mask = sm_ptr;
		*sign_mask++ = (uint64)0x7FFFFFFFFFFFFFFFull;
		*sign_mask-- = (uint64)0x7FFFFFFFFFFFFFFFull;

	#if 0
		// Set up the quadrupled-32-bit-int SSE constants used by the carry macros:
		sse_bw  = sm_ptr + 2;
		__asm	mov	eax, bw
		__asm	mov	ebx, sse_bw
		__asm	movd	xmm0,eax	/* Move actual *value* of reg eax into low 32 bits of xmm0 */
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0

		sse_sw  = sm_ptr + 4;
		__asm	lea	eax, sw
		__asm	mov	ebx, sse_sw
		__asm	movd	xmm0,[eax]	/* Variant 2: Move contents of address pointed to by reg eax into low 32 bits of xmm0 */
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0

		sse_n   = sm_ptr + 6;
		__asm	lea	eax, n
		__asm	mov	ebx, sse_n
		__asm	movd	xmm0,[eax]
		__asm	pshufd	xmm0,xmm0,0	// Broadcast low 32 bits of xmm0 to all 4 slots of xmm0
		__asm	movaps	[ebx],xmm0
	#else
		sse_bw  = sm_ptr + 2;
		tmp64 = (uint64)bw;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_bw++ = tmp64;
		*sse_bw-- = tmp64;

		sse_sw  = sm_ptr + 4;
		tmp64 = (uint64)sw;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_sw++ = tmp64;
		*sse_sw-- = tmp64;

		sse_n   = sm_ptr + 6;
		tmp64 = (uint64)n;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_n++ = tmp64;
		*sse_n-- = tmp64;
	#endif

#ifdef USE_PTHREAD
	/* Populate the elements of the thread-specific data structs which don't change after init: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
	// int data:
		tdat[ithread].tid = ithread;
		tdat[ithread].ndivr = NDIVR;

		tdat[ithread].sw  = sw;
		tdat[ithread].nwt = nwt;

	// pointer data:
		tdat[ithread].arrdat = a;			/* Main data array */
		tdat[ithread].wt0 = wt0;
		tdat[ithread].wt1 = wt1;
		tdat[ithread].si  = si;
		tdat[ithread].s1p00r = __r0 + ithread*cslots_in_local_store;
		tdat[ithread].half_arr = (long)tdat[ithread].s1p00r + ((long)half_arr - (long)s1p00r);
	}
#endif

		bjmodn00 = (uint32*)(sm_ptr + 8);
		bjmodn01 = bjmodn00 + 1;
		bjmodn02 = bjmodn01 + 1;
		bjmodn03 = bjmodn02 + 1;
		bjmodn04 = bjmodn03 + 1;
		bjmodn05 = bjmodn04 + 1;
		bjmodn06 = bjmodn05 + 1;
		bjmodn07 = bjmodn06 + 1;
		bjmodn08 = bjmodn07 + 1;
		bjmodn09 = bjmodn08 + 1;
		bjmodn10 = bjmodn09 + 1;
		bjmodn11 = bjmodn10 + 1;
		bjmodn12 = bjmodn11 + 1;
		bjmodn13 = bjmodn12 + 1;
		bjmodn14 = bjmodn13 + 1;
		bjmodn15 = bjmodn14 + 1;
		bjmodn16 = bjmodn15 + 1;
		bjmodn17 = bjmodn16 + 1;
		bjmodn18 = bjmodn17 + 1;
		bjmodn19 = bjmodn18 + 1;
		bjmodn20 = bjmodn19 + 1;
		bjmodn21 = bjmodn20 + 1;
		bjmodn22 = bjmodn21 + 1;
		bjmodn23 = bjmodn22 + 1;
		bjmodn24 = bjmodn23 + 1;
		bjmodn25 = bjmodn24 + 1;
		bjmodn26 = bjmodn25 + 1;
		bjmodn27 = bjmodn26 + 1;
		bjmodn28 = bjmodn27 + 1;
		bjmodn29 = bjmodn28 + 1;
		bjmodn30 = bjmodn29 + 1;
		bjmodn31 = bjmodn30 + 1;
		bjmodn32 = bjmodn31 + 1;
		bjmodn33 = bjmodn32 + 1;
		bjmodn34 = bjmodn33 + 1;
		bjmodn35 = bjmodn34 + 1;
		bjmodn36 = bjmodn35 + 1;
		bjmodn37 = bjmodn36 + 1;
		bjmodn38 = bjmodn37 + 1;
		bjmodn39 = bjmodn38 + 1;
		bjmodn40 = bjmodn39 + 1;
		bjmodn41 = bjmodn40 + 1;
		bjmodn42 = bjmodn41 + 1;
		bjmodn43 = bjmodn42 + 1;

	  #ifdef USE_PTHREAD
		tmp = __r0 + cslots_in_local_store;
		/* Init thread 1-CY_THREADS's local stores and pointers: */
		for(i = 1; i < CY_THREADS; ++i) {
			/* Only care about the constants for each thread here, but easier to just copy the entire thread0 local store: */
			memcpy(tmp, __r0, cslots_in_local_store<<4);	// bytewise copy treats complex and uint64 subdata the same
			tmp += cslots_in_local_store;
		}
	  #endif

	#endif	// USE_SSE2

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		/* Check the various add/sub combinations used in the fused radix-4 asm macros to make sure padding doesn't break them: */
		ASSERT(HERE, p01+p01 == p02, "p01+p01 != p02");
		ASSERT(HERE, p28-p16 == p12, "p28-p16 != p12");
		ASSERT(HERE, p12+p28 == p40, "p12+p28 != p40");
		ASSERT(HERE, p40-p16 == p24, "p40-p16 != p24");
		ASSERT(HERE, p24-p16 == p08, "p24-p16 != p08");
		ASSERT(HERE, p08+p28 == p36, "p08+p28 != p36");
		ASSERT(HERE, p36-p16 == p20, "p36-p16 != p20");
		ASSERT(HERE, p20-p16 == p04, "p20-p16 != p04");
		ASSERT(HERE, p04+p28 == p32, "p04+p28 != p32");
		ASSERT(HERE, p32-p16 == p16, "p32-p16 != p16");

		if(_cy00)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn00); _bjmodn00 = 0x0;
			free((void *)_bjmodn01); _bjmodn01 = 0x0;
			free((void *)_bjmodn02); _bjmodn02 = 0x0;
			free((void *)_bjmodn03); _bjmodn03 = 0x0;
			free((void *)_bjmodn04); _bjmodn04 = 0x0;
			free((void *)_bjmodn05); _bjmodn05 = 0x0;
			free((void *)_bjmodn06); _bjmodn06 = 0x0;
			free((void *)_bjmodn07); _bjmodn07 = 0x0;
			free((void *)_bjmodn08); _bjmodn08 = 0x0;
			free((void *)_bjmodn09); _bjmodn09 = 0x0;
			free((void *)_bjmodn10); _bjmodn10 = 0x0;
			free((void *)_bjmodn11); _bjmodn11 = 0x0;
			free((void *)_bjmodn12); _bjmodn12 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;
			free((void *)_bjmodn20); _bjmodn20 = 0x0;
			free((void *)_bjmodn21); _bjmodn21 = 0x0;
			free((void *)_bjmodn22); _bjmodn22 = 0x0;
			free((void *)_bjmodn23); _bjmodn23 = 0x0;
			free((void *)_bjmodn24); _bjmodn24 = 0x0;
			free((void *)_bjmodn25); _bjmodn25 = 0x0;
			free((void *)_bjmodn26); _bjmodn26 = 0x0;
			free((void *)_bjmodn27); _bjmodn27 = 0x0;
			free((void *)_bjmodn28); _bjmodn28 = 0x0;
			free((void *)_bjmodn29); _bjmodn29 = 0x0;
			free((void *)_bjmodn30); _bjmodn30 = 0x0;
			free((void *)_bjmodn31); _bjmodn31 = 0x0;
			free((void *)_bjmodn32); _bjmodn32 = 0x0;
			free((void *)_bjmodn33); _bjmodn33 = 0x0;
			free((void *)_bjmodn34); _bjmodn34 = 0x0;
			free((void *)_bjmodn35); _bjmodn35 = 0x0;
			free((void *)_bjmodn36); _bjmodn36 = 0x0;
			free((void *)_bjmodn37); _bjmodn37 = 0x0;
			free((void *)_bjmodn38); _bjmodn38 = 0x0;
			free((void *)_bjmodn39); _bjmodn39 = 0x0;
			free((void *)_bjmodn40); _bjmodn40 = 0x0;
			free((void *)_bjmodn41); _bjmodn41 = 0x0;
			free((void *)_bjmodn42); _bjmodn42 = 0x0;
			free((void *)_bjmodn43); _bjmodn43 = 0x0;

			free((void *)_cy00); _cy00 = 0x0;
			free((void *)_cy01); _cy01 = 0x0;
			free((void *)_cy02); _cy02 = 0x0;
			free((void *)_cy03); _cy03 = 0x0;
			free((void *)_cy04); _cy04 = 0x0;
			free((void *)_cy05); _cy05 = 0x0;
			free((void *)_cy06); _cy06 = 0x0;
			free((void *)_cy07); _cy07 = 0x0;
			free((void *)_cy08); _cy08 = 0x0;
			free((void *)_cy09); _cy09 = 0x0;
			free((void *)_cy10); _cy10 = 0x0;
			free((void *)_cy11); _cy11 = 0x0;
			free((void *)_cy12); _cy12 = 0x0;
			free((void *)_cy13); _cy13 = 0x0;
			free((void *)_cy14); _cy14 = 0x0;
			free((void *)_cy15); _cy15 = 0x0;
			free((void *)_cy16); _cy16 = 0x0;
			free((void *)_cy17); _cy17 = 0x0;
			free((void *)_cy18); _cy18 = 0x0;
			free((void *)_cy19); _cy19 = 0x0;
			free((void *)_cy20); _cy20 = 0x0;
			free((void *)_cy21); _cy21 = 0x0;
			free((void *)_cy22); _cy22 = 0x0;
			free((void *)_cy23); _cy23 = 0x0;
			free((void *)_cy24); _cy24 = 0x0;
			free((void *)_cy25); _cy25 = 0x0;
			free((void *)_cy26); _cy26 = 0x0;
			free((void *)_cy27); _cy27 = 0x0;
			free((void *)_cy28); _cy28 = 0x0;
			free((void *)_cy29); _cy29 = 0x0;
			free((void *)_cy30); _cy30 = 0x0;
			free((void *)_cy31); _cy31 = 0x0;
			free((void *)_cy32); _cy32 = 0x0;
			free((void *)_cy33); _cy33 = 0x0;
			free((void *)_cy34); _cy34 = 0x0;
			free((void *)_cy35); _cy35 = 0x0;
			free((void *)_cy36); _cy36 = 0x0;
			free((void *)_cy37); _cy37 = 0x0;
			free((void *)_cy38); _cy38 = 0x0;
			free((void *)_cy39); _cy39 = 0x0;
			free((void *)_cy40); _cy40 = 0x0;
			free((void *)_cy41); _cy41 = 0x0;
			free((void *)_cy42); _cy42 = 0x0;
			free((void *)_cy43); _cy43 = 0x0;

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		j = CY_THREADS*sizeof(int);
		_i       	= (int *)malloc(j);	ptr_prod += (uint32)(_i== 0x0);
		_bjmodn00	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn00== 0x0);
		_bjmodn01	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn01== 0x0);
		_bjmodn02	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn02== 0x0);
		_bjmodn03	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn03== 0x0);
		_bjmodn04	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn04== 0x0);
		_bjmodn05	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn05== 0x0);
		_bjmodn06	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn06== 0x0);
		_bjmodn07	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn07== 0x0);
		_bjmodn08	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn08== 0x0);
		_bjmodn09	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn09== 0x0);
		_bjmodn10	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn10== 0x0);
		_bjmodn11	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn11== 0x0);
		_bjmodn12	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn12== 0x0);
		_bjmodn13	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn13== 0x0);
		_bjmodn14	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn14== 0x0);
		_bjmodn15	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn15== 0x0);
		_bjmodn16	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn16== 0x0);
		_bjmodn17	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn17== 0x0);
		_bjmodn18	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn18== 0x0);
		_bjmodn19	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn19== 0x0);
		_bjmodn20	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_bjmodn24	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_bjmodn28	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn28== 0x0);
		_bjmodn29	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn29== 0x0);
		_bjmodn30	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn30== 0x0);
		_bjmodn31	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn31== 0x0);
		_bjmodn32	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn32== 0x0);
		_bjmodn33	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn33== 0x0);
		_bjmodn34	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn34== 0x0);
		_bjmodn35	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn35== 0x0);
		_bjmodn36	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn36== 0x0);
		_bjmodn37	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn37== 0x0);
		_bjmodn38	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn38== 0x0);
		_bjmodn39	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn39== 0x0);
		_bjmodn40	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn40== 0x0);
		_bjmodn41	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn41== 0x0);
		_bjmodn42	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn42== 0x0);
		_bjmodn43	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn43== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy00== 0x0);
		_cy01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy01== 0x0);
		_cy02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy02== 0x0);
		_cy03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy03== 0x0);
		_cy04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy04== 0x0);
		_cy05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy05== 0x0);
		_cy06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy06== 0x0);
		_cy07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy07== 0x0);
		_cy08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy08== 0x0);
		_cy09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy09== 0x0);
		_cy10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy10== 0x0);
		_cy11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy11== 0x0);
		_cy12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy12== 0x0);
		_cy13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy13== 0x0);
		_cy14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy14== 0x0);
		_cy15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy15== 0x0);
		_cy16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy16== 0x0);
		_cy17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy17== 0x0);
		_cy18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy18== 0x0);
		_cy19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy19== 0x0);
		_cy20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy20== 0x0);
		_cy21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy21== 0x0);
		_cy22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy22== 0x0);
		_cy23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy23== 0x0);
		_cy24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy24== 0x0);
		_cy25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy25== 0x0);
		_cy26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy26== 0x0);
		_cy27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy27== 0x0);
		_cy28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy28== 0x0);
		_cy29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy29== 0x0);
		_cy30	= (double *)malloc(j);	ptr_prod += (uint32)(_cy30== 0x0);
		_cy31	= (double *)malloc(j);	ptr_prod += (uint32)(_cy31== 0x0);
		_cy32	= (double *)malloc(j);	ptr_prod += (uint32)(_cy32== 0x0);
		_cy33	= (double *)malloc(j);	ptr_prod += (uint32)(_cy33== 0x0);
		_cy34	= (double *)malloc(j);	ptr_prod += (uint32)(_cy34== 0x0);
		_cy35	= (double *)malloc(j);	ptr_prod += (uint32)(_cy35== 0x0);
		_cy36	= (double *)malloc(j);	ptr_prod += (uint32)(_cy36== 0x0);
		_cy37	= (double *)malloc(j);	ptr_prod += (uint32)(_cy37== 0x0);
		_cy38	= (double *)malloc(j);	ptr_prod += (uint32)(_cy38== 0x0);
		_cy39	= (double *)malloc(j);	ptr_prod += (uint32)(_cy39== 0x0);
		_cy40	= (double *)malloc(j);	ptr_prod += (uint32)(_cy40== 0x0);
		_cy41	= (double *)malloc(j);	ptr_prod += (uint32)(_cy41== 0x0);
		_cy42	= (double *)malloc(j);	ptr_prod += (uint32)(_cy42== 0x0);
		_cy43	= (double *)malloc(j);	ptr_prod += (uint32)(_cy43== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix44_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/44-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix44_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		jhi = NDIVR/CY_THREADS;

		for(j=0; j < jhi; j++)
		{
			_bjmodnini[1] -= sw; _bjmodnini[1] = _bjmodnini[1] + ( (-(int)((uint32)_bjmodnini[1] >> 31)) & n);
		}

		if(CY_THREADS > 1)
		{
			for(ithread = 2; ithread <= CY_THREADS; ithread++)
			{
				_bjmodnini[ithread] = _bjmodnini[ithread-1] + _bjmodnini[1] - n; _bjmodnini[ithread] = _bjmodnini[ithread] + ( (-(int)((uint32)_bjmodnini[ithread] >> 31)) & n);
			}
		}
		/* Check upper element against scalar value, as precomputed in single-thread mode: */
		bjmodnini=0;
		for(j=0; j < jhi*CY_THREADS; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	}	/* endif(first_entry) */

/*...The radix-44 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy00[ithread] = 0;
		_cy01[ithread] = 0;
		_cy02[ithread] = 0;
		_cy03[ithread] = 0;
		_cy04[ithread] = 0;
		_cy05[ithread] = 0;
		_cy06[ithread] = 0;
		_cy07[ithread] = 0;
		_cy08[ithread] = 0;
		_cy09[ithread] = 0;
		_cy10[ithread] = 0;
		_cy11[ithread] = 0;
		_cy12[ithread] = 0;
		_cy13[ithread] = 0;
		_cy14[ithread] = 0;
		_cy15[ithread] = 0;
		_cy16[ithread] = 0;
		_cy17[ithread] = 0;
		_cy18[ithread] = 0;
		_cy19[ithread] = 0;
		_cy20[ithread] = 0;
		_cy21[ithread] = 0;
		_cy22[ithread] = 0;
		_cy23[ithread] = 0;
		_cy24[ithread] = 0;
		_cy25[ithread] = 0;
		_cy26[ithread] = 0;
		_cy27[ithread] = 0;
		_cy28[ithread] = 0;
		_cy29[ithread] = 0;
		_cy30[ithread] = 0;
		_cy31[ithread] = 0;
		_cy32[ithread] = 0;
		_cy33[ithread] = 0;
		_cy34[ithread] = 0;
		_cy35[ithread] = 0;
		_cy36[ithread] = 0;
		_cy37[ithread] = 0;
		_cy38[ithread] = 0;
		_cy39[ithread] = 0;
		_cy40[ithread] = 0;
		_cy41[ithread] = 0;
		_cy42[ithread] = 0;
		_cy43[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy00[      0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = 0.0;
	}

for(outer=0; outer <= 1; outer++)
{
	_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (_i[0] = 1).	*/

	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	/*
	Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
	then simply overwrite it with 1 prior to starting the k-loop.
	*/
	khi = n_div_nwt/CY_THREADS;
	j = _bjmodnini[CY_THREADS];
	// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn00[ithread] = _bjmodnini[ithread];
		MOD_ADD32(_bjmodn00[ithread], j, n, _bjmodn01[ithread]);
		MOD_ADD32(_bjmodn01[ithread], j, n, _bjmodn02[ithread]);
		MOD_ADD32(_bjmodn02[ithread], j, n, _bjmodn03[ithread]);
		MOD_ADD32(_bjmodn03[ithread], j, n, _bjmodn04[ithread]);
		MOD_ADD32(_bjmodn04[ithread], j, n, _bjmodn05[ithread]);
		MOD_ADD32(_bjmodn05[ithread], j, n, _bjmodn06[ithread]);
		MOD_ADD32(_bjmodn06[ithread], j, n, _bjmodn07[ithread]);
		MOD_ADD32(_bjmodn07[ithread], j, n, _bjmodn08[ithread]);
		MOD_ADD32(_bjmodn08[ithread], j, n, _bjmodn09[ithread]);
		MOD_ADD32(_bjmodn09[ithread], j, n, _bjmodn10[ithread]);
		MOD_ADD32(_bjmodn10[ithread], j, n, _bjmodn11[ithread]);
		MOD_ADD32(_bjmodn11[ithread], j, n, _bjmodn12[ithread]);
		MOD_ADD32(_bjmodn12[ithread], j, n, _bjmodn13[ithread]);
		MOD_ADD32(_bjmodn13[ithread], j, n, _bjmodn14[ithread]);
		MOD_ADD32(_bjmodn14[ithread], j, n, _bjmodn15[ithread]);
		MOD_ADD32(_bjmodn15[ithread], j, n, _bjmodn16[ithread]);
		MOD_ADD32(_bjmodn16[ithread], j, n, _bjmodn17[ithread]);
		MOD_ADD32(_bjmodn17[ithread], j, n, _bjmodn18[ithread]);
		MOD_ADD32(_bjmodn18[ithread], j, n, _bjmodn19[ithread]);
		MOD_ADD32(_bjmodn19[ithread], j, n, _bjmodn20[ithread]);
		MOD_ADD32(_bjmodn20[ithread], j, n, _bjmodn21[ithread]);
		MOD_ADD32(_bjmodn21[ithread], j, n, _bjmodn22[ithread]);
		MOD_ADD32(_bjmodn22[ithread], j, n, _bjmodn23[ithread]);
		MOD_ADD32(_bjmodn23[ithread], j, n, _bjmodn24[ithread]);
		MOD_ADD32(_bjmodn24[ithread], j, n, _bjmodn25[ithread]);
		MOD_ADD32(_bjmodn25[ithread], j, n, _bjmodn26[ithread]);
		MOD_ADD32(_bjmodn26[ithread], j, n, _bjmodn27[ithread]);
		MOD_ADD32(_bjmodn27[ithread], j, n, _bjmodn28[ithread]);
		MOD_ADD32(_bjmodn28[ithread], j, n, _bjmodn29[ithread]);
		MOD_ADD32(_bjmodn29[ithread], j, n, _bjmodn30[ithread]);
		MOD_ADD32(_bjmodn30[ithread], j, n, _bjmodn31[ithread]);
		MOD_ADD32(_bjmodn31[ithread], j, n, _bjmodn32[ithread]);
		MOD_ADD32(_bjmodn32[ithread], j, n, _bjmodn33[ithread]);
		MOD_ADD32(_bjmodn33[ithread], j, n, _bjmodn34[ithread]);
		MOD_ADD32(_bjmodn34[ithread], j, n, _bjmodn35[ithread]);
		MOD_ADD32(_bjmodn35[ithread], j, n, _bjmodn36[ithread]);
		MOD_ADD32(_bjmodn36[ithread], j, n, _bjmodn37[ithread]);
		MOD_ADD32(_bjmodn37[ithread], j, n, _bjmodn38[ithread]);
		MOD_ADD32(_bjmodn38[ithread], j, n, _bjmodn39[ithread]);
		MOD_ADD32(_bjmodn39[ithread], j, n, _bjmodn40[ithread]);
		MOD_ADD32(_bjmodn40[ithread], j, n, _bjmodn41[ithread]);
		MOD_ADD32(_bjmodn41[ithread], j, n, _bjmodn42[ithread]);
		MOD_ADD32(_bjmodn42[ithread], j, n, _bjmodn43[ithread]);

		_jstart[ithread] = ithread*NDIVR/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
	}

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	#if FFT_DEBUG
		fprintf(dbg_file, "radix44_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifdef USE_OMP
	#error radix44_ditN_cy_dif1: OpenMP private list needs updating!
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,j2,jt,jp,jstart,jhi,k,khi,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,rt,it,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23,cy24,cy25,cy26,cy27,cy28,cy29,cy30,cy31,cy32,cy33,cy34,cy35,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,cy_i28,cy_i29,cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35) default(shared) schedule(static)
*******************to-do: update private-vars list!*******************
#endif

#ifdef USE_PTHREAD
	/* Populate the thread-specific data structs - use the invariant terms as memchecks: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
	// int data:
		ASSERT(HERE, tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");

		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(HERE, tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].nwt == nwt, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = _maxerr[ithread];
		tdat[ithread].scale = scale;

	// pointer data:
		ASSERT(HERE, tdat[ithread].arrdat == a, "thread-local memcheck fail!");			/* Main data array */
		ASSERT(HERE, tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].s1p00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->re == crnd && (tmp-1)->im == crnd), "thread-local memcheck failed!");
		ASSERT(HERE, (tmp+10)->re * (tmp+14)->re == 1.0 && (tmp+10)->im * (tmp+14)->im == 1.0, "thread-local memcheck failed!");

		tdat[ithread].bjmodn00 = _bjmodn00[ithread];
		tdat[ithread].bjmodn01 = _bjmodn01[ithread];
		tdat[ithread].bjmodn02 = _bjmodn02[ithread];
		tdat[ithread].bjmodn03 = _bjmodn03[ithread];
		tdat[ithread].bjmodn04 = _bjmodn04[ithread];
		tdat[ithread].bjmodn05 = _bjmodn05[ithread];
		tdat[ithread].bjmodn06 = _bjmodn06[ithread];
		tdat[ithread].bjmodn07 = _bjmodn07[ithread];
		tdat[ithread].bjmodn08 = _bjmodn08[ithread];
		tdat[ithread].bjmodn09 = _bjmodn09[ithread];
		tdat[ithread].bjmodn10 = _bjmodn10[ithread];
		tdat[ithread].bjmodn11 = _bjmodn11[ithread];
		tdat[ithread].bjmodn12 = _bjmodn12[ithread];
		tdat[ithread].bjmodn13 = _bjmodn13[ithread];
		tdat[ithread].bjmodn14 = _bjmodn14[ithread];
		tdat[ithread].bjmodn15 = _bjmodn15[ithread];
		tdat[ithread].bjmodn16 = _bjmodn16[ithread];
		tdat[ithread].bjmodn17 = _bjmodn17[ithread];
		tdat[ithread].bjmodn18 = _bjmodn18[ithread];
		tdat[ithread].bjmodn19 = _bjmodn19[ithread];
		tdat[ithread].bjmodn20 = _bjmodn20[ithread];
		tdat[ithread].bjmodn21 = _bjmodn21[ithread];
		tdat[ithread].bjmodn22 = _bjmodn22[ithread];
		tdat[ithread].bjmodn23 = _bjmodn23[ithread];
		tdat[ithread].bjmodn24 = _bjmodn24[ithread];
		tdat[ithread].bjmodn25 = _bjmodn25[ithread];
		tdat[ithread].bjmodn26 = _bjmodn26[ithread];
		tdat[ithread].bjmodn27 = _bjmodn27[ithread];
		tdat[ithread].bjmodn28 = _bjmodn28[ithread];
		tdat[ithread].bjmodn29 = _bjmodn29[ithread];
		tdat[ithread].bjmodn30 = _bjmodn30[ithread];
		tdat[ithread].bjmodn31 = _bjmodn31[ithread];
		tdat[ithread].bjmodn32 = _bjmodn32[ithread];
		tdat[ithread].bjmodn33 = _bjmodn33[ithread];
		tdat[ithread].bjmodn34 = _bjmodn34[ithread];
		tdat[ithread].bjmodn35 = _bjmodn35[ithread];
		tdat[ithread].bjmodn36 = _bjmodn36[ithread];
		tdat[ithread].bjmodn37 = _bjmodn37[ithread];
		tdat[ithread].bjmodn38 = _bjmodn38[ithread];
		tdat[ithread].bjmodn39 = _bjmodn39[ithread];
		tdat[ithread].bjmodn40 = _bjmodn40[ithread];
		tdat[ithread].bjmodn41 = _bjmodn41[ithread];
		tdat[ithread].bjmodn42 = _bjmodn42[ithread];
		tdat[ithread].bjmodn43 = _bjmodn43[ithread];
		/* init carries	*/
		tdat[ithread].cy00 = _cy00[ithread];
		tdat[ithread].cy01 = _cy01[ithread];
		tdat[ithread].cy02 = _cy02[ithread];
		tdat[ithread].cy03 = _cy03[ithread];
		tdat[ithread].cy04 = _cy04[ithread];
		tdat[ithread].cy05 = _cy05[ithread];
		tdat[ithread].cy06 = _cy06[ithread];
		tdat[ithread].cy07 = _cy07[ithread];
		tdat[ithread].cy08 = _cy08[ithread];
		tdat[ithread].cy09 = _cy09[ithread];
		tdat[ithread].cy10 = _cy10[ithread];
		tdat[ithread].cy11 = _cy11[ithread];
		tdat[ithread].cy12 = _cy12[ithread];
		tdat[ithread].cy13 = _cy13[ithread];
		tdat[ithread].cy14 = _cy14[ithread];
		tdat[ithread].cy15 = _cy15[ithread];
		tdat[ithread].cy16 = _cy16[ithread];
		tdat[ithread].cy17 = _cy17[ithread];
		tdat[ithread].cy18 = _cy18[ithread];
		tdat[ithread].cy19 = _cy19[ithread];
		tdat[ithread].cy20 = _cy20[ithread];
		tdat[ithread].cy21 = _cy21[ithread];
		tdat[ithread].cy22 = _cy22[ithread];
		tdat[ithread].cy23 = _cy23[ithread];
		tdat[ithread].cy24 = _cy24[ithread];
		tdat[ithread].cy25 = _cy25[ithread];
		tdat[ithread].cy26 = _cy26[ithread];
		tdat[ithread].cy27 = _cy27[ithread];
		tdat[ithread].cy28 = _cy28[ithread];
		tdat[ithread].cy29 = _cy29[ithread];
		tdat[ithread].cy30 = _cy30[ithread];
		tdat[ithread].cy31 = _cy31[ithread];
		tdat[ithread].cy32 = _cy32[ithread];
		tdat[ithread].cy33 = _cy33[ithread];
		tdat[ithread].cy34 = _cy34[ithread];
		tdat[ithread].cy35 = _cy35[ithread];
		tdat[ithread].cy36 = _cy36[ithread];
		tdat[ithread].cy37 = _cy37[ithread];
		tdat[ithread].cy38 = _cy38[ithread];
		tdat[ithread].cy39 = _cy39[ithread];
		tdat[ithread].cy40 = _cy40[ithread];
		tdat[ithread].cy41 = _cy41[ithread];
		tdat[ithread].cy42 = _cy42[ithread];
		tdat[ithread].cy43 = _cy43[ithread];
	}
#endif

#ifdef USE_PTHREAD

	// If also using main thread to do work units, that task-dispatch occurs after all the threadpool-task launches:
	for(ithread = 0; ithread < pool_work_units; ithread++)
	{
		task_control.data = (void*)(&tdat[ithread]);
		threadpool_add_task(tpool, &task_control, task_is_blocking);

#else

    for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		/***** DEC/HP CC doesn't properly copy init value of maxerr = 0 into threads,
		so need to set once again explicitly for each: *****/
		maxerr = 0.0;
	#ifdef USE_SSE2
		max_err->re = 0.0;	max_err->im = 0.0;
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		col = _col[ithread];
		co2 = _co2[ithread];
		co3 = _co3[ithread];

	#ifdef USE_SSE2								/* init carries	*/
		*bjmodn00 = _bjmodn00[ithread];			cy00->re = _cy00[ithread];
		*bjmodn01 = _bjmodn01[ithread];			cy00->im = _cy01[ithread];
		*bjmodn02 = _bjmodn02[ithread];			cy02->re = _cy02[ithread];
		*bjmodn03 = _bjmodn03[ithread];			cy02->im = _cy03[ithread];
		*bjmodn04 = _bjmodn04[ithread];			cy04->re = _cy04[ithread];
		*bjmodn05 = _bjmodn05[ithread];			cy04->im = _cy05[ithread];
		*bjmodn06 = _bjmodn06[ithread];			cy06->re = _cy06[ithread];
		*bjmodn07 = _bjmodn07[ithread];			cy06->im = _cy07[ithread];
		*bjmodn08 = _bjmodn08[ithread];			cy08->re = _cy08[ithread];
		*bjmodn09 = _bjmodn09[ithread];			cy08->im = _cy09[ithread];
		*bjmodn10 = _bjmodn10[ithread];			cy10->re = _cy10[ithread];
		*bjmodn11 = _bjmodn11[ithread];			cy10->im = _cy11[ithread];
		*bjmodn12 = _bjmodn12[ithread];			cy12->re = _cy12[ithread];
		*bjmodn13 = _bjmodn13[ithread];			cy12->im = _cy13[ithread];
		*bjmodn14 = _bjmodn14[ithread];			cy14->re = _cy14[ithread];
		*bjmodn15 = _bjmodn15[ithread];			cy14->im = _cy15[ithread];
		*bjmodn16 = _bjmodn16[ithread];			cy16->re = _cy16[ithread];
		*bjmodn17 = _bjmodn17[ithread];			cy16->im = _cy17[ithread];
		*bjmodn18 = _bjmodn18[ithread];			cy18->re = _cy18[ithread];
		*bjmodn19 = _bjmodn19[ithread];			cy18->im = _cy19[ithread];
		*bjmodn20 = _bjmodn20[ithread];			cy20->re = _cy20[ithread];
		*bjmodn21 = _bjmodn21[ithread];			cy20->im = _cy21[ithread];
		*bjmodn22 = _bjmodn22[ithread];			cy22->re = _cy22[ithread];
		*bjmodn23 = _bjmodn23[ithread];			cy22->im = _cy23[ithread];
		*bjmodn24 = _bjmodn24[ithread];			cy24->re = _cy24[ithread];
		*bjmodn25 = _bjmodn25[ithread];			cy24->im = _cy25[ithread];
		*bjmodn26 = _bjmodn26[ithread];			cy26->re = _cy26[ithread];
		*bjmodn27 = _bjmodn27[ithread];			cy26->im = _cy27[ithread];
		*bjmodn28 = _bjmodn28[ithread];			cy28->re = _cy28[ithread];
		*bjmodn29 = _bjmodn29[ithread];			cy28->im = _cy29[ithread];
		*bjmodn30 = _bjmodn30[ithread];			cy30->re = _cy30[ithread];
		*bjmodn31 = _bjmodn31[ithread];			cy30->im = _cy31[ithread];
		*bjmodn32 = _bjmodn32[ithread];			cy32->re = _cy32[ithread];
		*bjmodn33 = _bjmodn33[ithread];			cy32->im = _cy33[ithread];
		*bjmodn34 = _bjmodn34[ithread];			cy34->re = _cy34[ithread];
		*bjmodn35 = _bjmodn35[ithread];			cy34->im = _cy35[ithread];
		*bjmodn36 = _bjmodn36[ithread];			cy36->re = _cy36[ithread];
		*bjmodn37 = _bjmodn37[ithread];			cy36->im = _cy37[ithread];
		*bjmodn38 = _bjmodn38[ithread];			cy38->re = _cy38[ithread];
		*bjmodn39 = _bjmodn39[ithread];			cy38->im = _cy39[ithread];
		*bjmodn40 = _bjmodn40[ithread];			cy40->re = _cy40[ithread];
		*bjmodn41 = _bjmodn41[ithread];			cy40->im = _cy41[ithread];
		*bjmodn42 = _bjmodn42[ithread];			cy42->re = _cy42[ithread];
		*bjmodn43 = _bjmodn43[ithread];			cy42->im = _cy43[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];			cy00 = _cy00[ithread];
		bjmodn01 = _bjmodn01[ithread];			cy01 = _cy01[ithread];
		bjmodn02 = _bjmodn02[ithread];			cy02 = _cy02[ithread];
		bjmodn03 = _bjmodn03[ithread];			cy03 = _cy03[ithread];
		bjmodn04 = _bjmodn04[ithread];			cy04 = _cy04[ithread];
		bjmodn05 = _bjmodn05[ithread];			cy05 = _cy05[ithread];
		bjmodn06 = _bjmodn06[ithread];			cy06 = _cy06[ithread];
		bjmodn07 = _bjmodn07[ithread];			cy07 = _cy07[ithread];
		bjmodn08 = _bjmodn08[ithread];			cy08 = _cy08[ithread];
		bjmodn09 = _bjmodn09[ithread];			cy09 = _cy09[ithread];
		bjmodn10 = _bjmodn10[ithread];			cy10 = _cy10[ithread];
		bjmodn11 = _bjmodn11[ithread];			cy11 = _cy11[ithread];
		bjmodn12 = _bjmodn12[ithread];			cy12 = _cy12[ithread];
		bjmodn13 = _bjmodn13[ithread];			cy13 = _cy13[ithread];
		bjmodn14 = _bjmodn14[ithread];			cy14 = _cy14[ithread];
		bjmodn15 = _bjmodn15[ithread];			cy15 = _cy15[ithread];
		bjmodn16 = _bjmodn16[ithread];			cy16 = _cy16[ithread];
		bjmodn17 = _bjmodn17[ithread];			cy17 = _cy17[ithread];
		bjmodn18 = _bjmodn18[ithread];			cy18 = _cy18[ithread];
		bjmodn19 = _bjmodn19[ithread];			cy19 = _cy19[ithread];
		bjmodn20 = _bjmodn20[ithread];			cy20 = _cy20[ithread];
		bjmodn21 = _bjmodn21[ithread];			cy21 = _cy21[ithread];
		bjmodn22 = _bjmodn22[ithread];			cy22 = _cy22[ithread];
		bjmodn23 = _bjmodn23[ithread];			cy23 = _cy23[ithread];
		bjmodn24 = _bjmodn24[ithread];			cy24 = _cy24[ithread];
		bjmodn25 = _bjmodn25[ithread];			cy25 = _cy25[ithread];
		bjmodn26 = _bjmodn26[ithread];			cy26 = _cy26[ithread];
		bjmodn27 = _bjmodn27[ithread];			cy27 = _cy27[ithread];
		bjmodn28 = _bjmodn28[ithread];			cy28 = _cy28[ithread];
		bjmodn29 = _bjmodn29[ithread];			cy29 = _cy29[ithread];
		bjmodn30 = _bjmodn30[ithread];			cy30 = _cy30[ithread];
		bjmodn31 = _bjmodn31[ithread];			cy31 = _cy31[ithread];
		bjmodn32 = _bjmodn32[ithread];			cy32 = _cy32[ithread];
		bjmodn33 = _bjmodn33[ithread];			cy33 = _cy33[ithread];
		bjmodn34 = _bjmodn34[ithread];			cy34 = _cy34[ithread];
		bjmodn35 = _bjmodn35[ithread];			cy35 = _cy35[ithread];
		bjmodn36 = _bjmodn36[ithread];			cy36 = _cy36[ithread];
		bjmodn37 = _bjmodn37[ithread];			cy37 = _cy37[ithread];
		bjmodn38 = _bjmodn38[ithread];			cy38 = _cy38[ithread];
		bjmodn39 = _bjmodn39[ithread];			cy39 = _cy39[ithread];
		bjmodn40 = _bjmodn40[ithread];			cy40 = _cy40[ithread];
		bjmodn41 = _bjmodn41[ithread];			cy41 = _cy41[ithread];
		bjmodn42 = _bjmodn42[ithread];			cy42 = _cy42[ithread];
		bjmodn43 = _bjmodn43[ithread];			cy43 = _cy43[ithread];
	#endif

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
		#ifdef USE_SSE2
			for(j = jstart; j < jhi; j += 4)
			{
			/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
			*/
				j1 = (j & mask01) + br4[j&3];
		#elif defined(USE_SSE2)	/* This allows us to use #if 0 above and disable sse2-based *computation*, while still using sse2-style data layout */
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
				j1 = (j & mask01) + br4[j&3];
		#else
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
				j1 =  j;
		#endif
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

		#ifdef USE_SSE2

		  #if GCC_ASM_FULL_INLINE	// GCC or SUNC implied

			add0 = &a[j1    ];
			SSE2_RADIX44_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,t00r, ua0, s1p00r);

		  #elif !defined(COMPILER_TYPE_MSVC) && !GCC_ASM_FULL_INLINE

			/* Outputs in SSE2 modes are temps 2*11*16 = 22*16 = 0x160 bytes apart: */
			add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t00r, 0x160)
			add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t01r, 0x160)
			add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t02r, 0x160)
			add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t03r, 0x160)
			add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t04r, 0x160)
			add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t05r, 0x160)
			add2 = &a[j1+p36];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t06r, 0x160)
			add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t07r, 0x160)
			add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t08r, 0x160)
			add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t09r, 0x160)
			add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t0ar, 0x160)

			/* Radix-11 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 12 (12*32 bytes = 0x180) or XX -= 32 (-32*32 bytes = -0x400) between successive outputs: */
																								/*   a1p00r,a1p12r,a1p24r,a1p36r,a1p04r,a1p16r,a1p28r,a1p40r,a1p08r,a1p20r,a1p32r */
			SSE2_RADIX_11_DFT(t00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p00r, 0x180, 0x300, 0x480, 0x080, 0x200, 0x380, 0x500, 0x100, 0x280, 0x400)
																								/*   a1p11r,a1p23r,a1p35r,a1p03r,a1p15r,a1p27r,a1p39r,a1p07r,a1p19r,a1p31r,a1p43r */
			SSE2_RADIX_11_DFT(t10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p11r, 0x180, 0x300,-0x100, 0x080, 0x200, 0x380,-0x080, 0x100, 0x280, 0x400)
																								/*   a1p22r,a1p34r,a1p02r,a1p14r,a1p26r,a1p38r,a1p06r,a1p18r,a1p30r,a1p42r,a1p10r */
			SSE2_RADIX_11_DFT(t20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p22r, 0x180,-0x280,-0x100, 0x080, 0x200,-0x200,-0x080, 0x100, 0x280,-0x180)
																								/*   a1p33r,a1p01r,a1p13r,a1p25r,a1p37r,a1p05r,a1p17r,a1p29r,a1p41r,a1p09r,a1p21r */
			SSE2_RADIX_11_DFT(t30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p33r,-0x400,-0x280,-0x100, 0x080,-0x380,-0x200,-0x080, 0x100,-0x300,-0x180)

		  #else

			add0 = &a[j1    ];
		//	add0,1,2,3 = &a[j1+p00]+p0,1,2,3
			__asm	mov	eax, add0	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_B/C */
			__asm	mov	edx, add0
			__asm	mov	esi, p01	/* esi will store power-of-2 multiples of p01 throughout */
			__asm	shl	esi, 3		/* Pointer offset for floating doubles */
			__asm	add edx, esi
			__asm	mov ebx, edx	/* add1 = add0+p01 */
			__asm	add edx, esi
			__asm	mov ecx, edx	/* add3 = add0+p02 */
			__asm	add edx, esi	/* add2 = add0+p03 */
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x160, 0x2c0, t00r)
		//	add3,2,1,0 = &a[j1+p28]+p0,1,2,3
			__asm	mov	esi, p28
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p28]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x160, 0x2c0, t01r)
		//	add3,2,1,0 = &a[j1+p12]+p0,1,2,3
			__asm	mov	esi, p16
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p12]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(edx,ecx,ebx,eax, 0x160, 0x2c0, t02r)
		//	add1,0,2,3 = &a[j1+p40]+p0,1,2,3
			__asm	mov	esi, p28
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p40]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x160, 0x2c0, t03r)
		//	add1,0,2,3 = &a[j1+p24]+p0,1,2,3
			__asm	mov	esi, p16
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p24]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x160, 0x2c0, t04r)
		//	add1,0,2,3 = &a[j1+p08]+p0,1,2,3
			__asm	mov	esi, p16
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p08]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ebx,eax,ecx,edx, 0x160, 0x2c0, t05r)
		//	add2,3,0,1 = &a[j1+p36]+p0,1,2,3
			__asm	mov	esi, p28
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p36]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x160, 0x2c0, t06r)
		//	add2,3,0,1 = &a[j1+p20]+p0,1,2,3
			__asm	mov	esi, p16
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p20]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x160, 0x2c0, t07r)
		//	add2,3,0,1 = &a[j1+p04]+p0,1,2,3
			__asm	mov	esi, p16
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p04]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(ecx,edx,eax,ebx, 0x160, 0x2c0, t08r)
		//	add0,1,3,2 = &a[j1+p32]+p0,1,2,3
			__asm	mov	esi, p28
			__asm	shl	esi, 3
			__asm	add	eax, esi// &a[j1+p32]
			__asm	add ebx, esi
			__asm	add ecx, esi
			__asm	add edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x160, 0x2c0, t09r)
		//	add0,1,3,2 = &a[j1+p16]+p0,1,2,3
			__asm	mov	esi, p16
			__asm	shl	esi, 3
			__asm	sub	eax, esi// &a[j1+p16]
			__asm	sub ebx, esi
			__asm	sub ecx, esi
			__asm	sub edx, esi
			SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x160, 0x2c0, t0ar)

			/* Radix-11 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 12 (24*16 bytes = 0x180) between successive outputs: */
																								/*   a1p00r,a1p12r,a1p24r,a1p36r,a1p04r,a1p16r,a1p28r,a1p40r,a1p08r,a1p20r,a1p32r */
			SSE2_RADIX_11_DFT(t00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p00r, 0x180, 0x300, 0x480, 0x080, 0x200, 0x380, 0x500, 0x100, 0x280, 0x400)
																								/*   a1p11r,a1p23r,a1p35r,a1p03r,a1p15r,a1p27r,a1p39r,a1p07r,a1p19r,a1p31r,a1p43r */
			SSE2_RADIX_11_DFT(t10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p11r, 0x180, 0x300,-0x100, 0x080, 0x200, 0x380,-0x080, 0x100, 0x280, 0x400)
																								/*   a1p22r,a1p34r,a1p02r,a1p14r,a1p26r,a1p38r,a1p06r,a1p18r,a1p30r,a1p42r,a1p10r */
			SSE2_RADIX_11_DFT(t20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p22r, 0x180,-0x280,-0x100, 0x080, 0x200,-0x200,-0x080, 0x100, 0x280,-0x180)
																								/*   a1p33r,a1p01r,a1p13r,a1p25r,a1p37r,a1p05r,a1p17r,a1p29r,a1p41r,a1p09r,a1p21r */
			SSE2_RADIX_11_DFT(t30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p33r,-0x400,-0x280,-0x100, 0x080,-0x380,-0x200,-0x080, 0x100,-0x300,-0x180)

		  #endif

		#else	/* !USE_SSE2 */

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 11 radix-4 transforms...*/
						/*                                   inputs                                   */ /*              outputs              */
			RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,rt,it);	jt = j1+p28; jp = j2+p28;
			RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,rt,it);	jt = j1+p12; jp = j2+p12;
			RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,rt,it);	jt = j1+p40; jp = j2+p40;
			RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,rt,it);	jt = j1+p24; jp = j2+p24;
			RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,rt,it);	jt = j1+p08; jp = j2+p08;
			RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,rt,it);	jt = j1+p36; jp = j2+p36;
			RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,rt,it);	jt = j1+p20; jp = j2+p20;
			RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,rt,it);	jt = j1+p04; jp = j2+p04;
			RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,rt,it);	jt = j1+p32; jp = j2+p32;
			RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,rt,it);	jt = j1+p16; jp = j2+p16;
			RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,rt,it);

			/*...and now do 4 radix-11 transforms...*/
						/*                                                  inputs                                                  */ /*                                                                                            outputs                                                                                                                    */
			jt = j1; jp = j2;

			RADIX_11_DFT(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,a1p00r,a1p00i,a1p12r,a1p12i,a1p24r,a1p24i,a1p36r,a1p36i,a1p04r,a1p04i,a1p16r,a1p16i,a1p28r,a1p28i,a1p40r,a1p40i,a1p08r,a1p08i,a1p20r,a1p20i,a1p32r,a1p32i);
			RADIX_11_DFT(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,a1p11r,a1p11i,a1p23r,a1p23i,a1p35r,a1p35i,a1p03r,a1p03i,a1p15r,a1p15i,a1p27r,a1p27i,a1p39r,a1p39i,a1p07r,a1p07i,a1p19r,a1p19i,a1p31r,a1p31i,a1p43r,a1p43i);
			RADIX_11_DFT(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,a1p22r,a1p22i,a1p34r,a1p34i,a1p02r,a1p02i,a1p14r,a1p14i,a1p26r,a1p26i,a1p38r,a1p38i,a1p06r,a1p06i,a1p18r,a1p18i,a1p30r,a1p30i,a1p42r,a1p42i,a1p10r,a1p10i);
			RADIX_11_DFT(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,a1p33r,a1p33i,a1p01r,a1p01i,a1p13r,a1p13i,a1p25r,a1p25i,a1p37r,a1p37i,a1p05r,a1p05i,a1p17r,a1p17i,a1p29r,a1p29i,a1p41r,a1p41i,a1p09r,a1p09i,a1p21r,a1p21i);

		#endif

	/*...Now do the carries. Since the outputs would
		normally be getting dispatched to 44 separate blocks of the A-array, we need 44 separate carries.	*/

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
		#ifdef USE_SSE2

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

		  #if defined(COMPILER_TYPE_MSVC)

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy20,cy22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy24,cy26,bjmodn24);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p28r,add1,add2,add3,cy28,cy30,bjmodn28);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p32r,add1,add2,add3,cy32,cy34,bjmodn32);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p36r,add1,add2,add3,cy36,cy38,bjmodn36);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p40r,add1,add2,add3,cy40,cy42,bjmodn40);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy20,cy22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy24,cy26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy28,cy30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy32,cy34,bjmodn32);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy36,cy38,bjmodn36);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p40r,add1,add2,add3,cy40,cy42,bjmodn40);
			  #endif

		  #else	/* GCC-style inline ASM: */

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p28r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p32r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p36r,add1,add2,add3,cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p40r,add1,add2,add3,cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p40r,add1,add2,add3,cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				/* Bizarre - when I disabled the diagnostic prints above and below, the resulting GCC build immediately gave
					fatal roundoff errors starting on iteration #5 - so insert the bogus [never taken] if() here as a workaround.
					Equally bizarre, inserting the bogus if() *before* the 4 carry-macro calls above gave the correct result as well,
					but ran fully 10% slower. Good old GCC...
				*/
				if(j < 0)
				{
					fprintf(stderr, "Iter %3d\n",iter);
				}

		  #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

		  #if defined(COMPILER_TYPE_MSVC)

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy20,cy22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy24,cy26,bjmodn24);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p28r,add1,add2,     cy28,cy30,bjmodn28);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p32r,add1,add2,     cy32,cy34,bjmodn32);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p36r,add1,add2,     cy36,cy38,bjmodn36);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p40r,add1,add2,     cy40,cy42,bjmodn40);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy20,cy22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy24,cy26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p28r,add1,add2,     cy28,cy30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p32r,add1,add2,     cy32,cy34,bjmodn32);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p36r,add1,add2,     cy36,cy38,bjmodn36);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p40r,add1,add2,     cy40,cy42,bjmodn40);
			  #endif

		  #else	/* GCC-style inline ASM: */

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p28r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p32r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p36r,add1,add2,     cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p40r,add1,add2,     cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p28r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p32r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p36r,add1,add2,     cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p40r,add1,add2,     cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

		  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

		#else	/* #ifdef USE_SSE2 */

				/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy01,bjmodn01,1 );
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy02,bjmodn02,2 );
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy03,bjmodn03,3 );
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy04,bjmodn04,4 );
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy05,bjmodn05,5 );
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy06,bjmodn06,6 );
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy07,bjmodn07,7 );
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy08,bjmodn08,8 );
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy09,bjmodn09,9 );
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy23,bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy24,bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy25,bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy26,bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy27,bjmodn27,27);
				cmplx_carry_norm_errcheck(a1p28r,a1p28i,cy28,bjmodn28,28);
				cmplx_carry_norm_errcheck(a1p29r,a1p29i,cy29,bjmodn29,29);
				cmplx_carry_norm_errcheck(a1p30r,a1p30i,cy30,bjmodn30,30);
				cmplx_carry_norm_errcheck(a1p31r,a1p31i,cy31,bjmodn31,31);
				cmplx_carry_norm_errcheck(a1p32r,a1p32i,cy32,bjmodn32,32);
				cmplx_carry_norm_errcheck(a1p33r,a1p33i,cy33,bjmodn33,33);
				cmplx_carry_norm_errcheck(a1p34r,a1p34i,cy34,bjmodn34,34);
				cmplx_carry_norm_errcheck(a1p35r,a1p35i,cy35,bjmodn35,35);
				cmplx_carry_norm_errcheck(a1p36r,a1p36i,cy36,bjmodn36,36);
				cmplx_carry_norm_errcheck(a1p37r,a1p37i,cy37,bjmodn37,37);
				cmplx_carry_norm_errcheck(a1p38r,a1p38i,cy38,bjmodn38,38);
				cmplx_carry_norm_errcheck(a1p39r,a1p39i,cy39,bjmodn39,39);
				cmplx_carry_norm_errcheck(a1p40r,a1p40i,cy40,bjmodn40,40);
				cmplx_carry_norm_errcheck(a1p41r,a1p41i,cy41,bjmodn41,41);
				cmplx_carry_norm_errcheck(a1p42r,a1p42i,cy42,bjmodn42,42);
				cmplx_carry_norm_errcheck(a1p43r,a1p43i,cy43,bjmodn43,43);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	/* #ifdef USE_SSE2 */

			/*...The radix-44 DIF pass is here:	*/

		#ifdef USE_SSE2

			/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-11 transforms...*/
			/* Radix-11 DFT inputs are (cyclic) with pXXr having XX += 40 (40*32 bytes = +0x500) or XX -= 4 (4*32 bytes = -0x080), outputs are adjacent 32-byte-separated temps: */
							/*a1p00r,a1p40r,a1p36r,a1p32r,a1p28r,a1p24r,a1p20r,a1p16r,a1p12r,a1p08r,a1p04r */
			SSE2_RADIX_11_DFT(s1p00r, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
							/*a1p33r,a1p29r,a1p25r,a1p21r,a1p17r,a1p13r,a1p09r,a1p05r,a1p01r,a1p41r,a1p37r */
			SSE2_RADIX_11_DFT(s1p33r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300,-0x380,-0x400, 0x100, 0x080, ua0, t10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
							/*a1p22r,a1p18r,a1p14r,a1p10r,a1p06r,a1p02r,a1p42r,a1p38r,a1p34r,a1p30r,a1p26r */
			SSE2_RADIX_11_DFT(s1p22r,-0x080,-0x100,-0x180,-0x200,-0x280, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
							/*a1p11r,a1p07r,a1p03r,a1p43r,a1p39r,a1p35r,a1p31r,a1p27r,a1p23r,a1p19r,a1p15r */
			SSE2_RADIX_11_DFT(s1p11r,-0x080,-0x100, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)

			/*...and now do 11 radix-4 transforms...*/
			/* Inputs in SSE2 modes are temps 2*11*16 = 22*16 = 0x160 bytes apart. Notice how the add* indices after the first row repeat with period 4: */
																															/*     outputs      */ /* inputs */
			add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t00r, 0x160)
			add1 = &a[j1+p40];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t01r, 0x160)
			add3 = &a[j1+p36];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t02r, 0x160)
			add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t03r, 0x160)
			add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t04r, 0x160)
			add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t05r, 0x160)
			add3 = &a[j1+p20];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t06r, 0x160)
			add0 = &a[j1+p16];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t07r, 0x160)
			add2 = &a[j1+p12];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t08r, 0x160)
			add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t09r, 0x160)
			add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t0ar, 0x160)

		#else	/* !USE_SSE2 */

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-9 transforms...*/
			RADIX_11_DFT(a1p00r,a1p00i,a1p40r,a1p40i,a1p36r,a1p36i,a1p32r,a1p32i,a1p28r,a1p28i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai);
			RADIX_11_DFT(a1p33r,a1p33i,a1p29r,a1p29i,a1p25r,a1p25i,a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p41r,a1p41i,a1p37r,a1p37i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai);
			RADIX_11_DFT(a1p22r,a1p22i,a1p18r,a1p18i,a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p42r,a1p42i,a1p38r,a1p38i,a1p34r,a1p34i,a1p30r,a1p30i,a1p26r,a1p26i,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai);
			RADIX_11_DFT(a1p11r,a1p11i,a1p07r,a1p07i,a1p03r,a1p03i,a1p43r,a1p43i,a1p39r,a1p39i,a1p35r,a1p35i,a1p31r,a1p31i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai);

			/*...and now do 11 radix-4 transforms...*/
						/*          inputs                    */ /*                                      outputs                              */
			RADIX_04_DIF(t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);	jt = j1+p40; jp = j2+p40;
			RADIX_04_DIF(t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p36; jp = j2+p36;
			RADIX_04_DIF(t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p32; jp = j2+p32;
			RADIX_04_DIF(t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p28; jp = j2+p28;
			RADIX_04_DIF(t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p24; jp = j2+p24;
			RADIX_04_DIF(t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p20; jp = j2+p20;
			RADIX_04_DIF(t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p16; jp = j2+p16;
			RADIX_04_DIF(t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p12; jp = j2+p12;
			RADIX_04_DIF(t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p08; jp = j2+p08;
			RADIX_04_DIF(t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p04; jp = j2+p04;
			RADIX_04_DIF(t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

		#endif

			}

			jstart += nwt;
			jhi    += nwt;

			col += RADIX;
			co3 -= RADIX;
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_SSE2
		_cy00[ithread] = cy00->re;
		_cy01[ithread] = cy00->im;
		_cy02[ithread] = cy02->re;
		_cy03[ithread] = cy02->im;
		_cy04[ithread] = cy04->re;
		_cy05[ithread] = cy04->im;
		_cy06[ithread] = cy06->re;
		_cy07[ithread] = cy06->im;
		_cy08[ithread] = cy08->re;
		_cy09[ithread] = cy08->im;
		_cy10[ithread] = cy10->re;
		_cy11[ithread] = cy10->im;
		_cy12[ithread] = cy12->re;
		_cy13[ithread] = cy12->im;
		_cy14[ithread] = cy14->re;
		_cy15[ithread] = cy14->im;
		_cy16[ithread] = cy16->re;
		_cy17[ithread] = cy16->im;
		_cy18[ithread] = cy18->re;
		_cy19[ithread] = cy18->im;
		_cy20[ithread] = cy20->re;
		_cy21[ithread] = cy20->im;
		_cy22[ithread] = cy22->re;
		_cy23[ithread] = cy22->im;
		_cy24[ithread] = cy24->re;
		_cy25[ithread] = cy24->im;
		_cy26[ithread] = cy26->re;
		_cy27[ithread] = cy26->im;
		_cy28[ithread] = cy28->re;
		_cy29[ithread] = cy28->im;
		_cy30[ithread] = cy30->re;
		_cy31[ithread] = cy30->im;
		_cy32[ithread] = cy32->re;
		_cy33[ithread] = cy32->im;
		_cy34[ithread] = cy34->re;
		_cy35[ithread] = cy34->im;
		_cy36[ithread] = cy36->re;
		_cy37[ithread] = cy36->im;
		_cy38[ithread] = cy38->re;
		_cy39[ithread] = cy38->im;
		_cy40[ithread] = cy40->re;
		_cy41[ithread] = cy40->im;
		_cy42[ithread] = cy42->re;
		_cy43[ithread] = cy42->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
	#else
		_cy00[ithread] = cy00;
		_cy01[ithread] = cy01;
		_cy02[ithread] = cy02;
		_cy03[ithread] = cy03;
		_cy04[ithread] = cy04;
		_cy05[ithread] = cy05;
		_cy06[ithread] = cy06;
		_cy07[ithread] = cy07;
		_cy08[ithread] = cy08;
		_cy09[ithread] = cy09;
		_cy10[ithread] = cy10;
		_cy11[ithread] = cy11;
		_cy12[ithread] = cy12;
		_cy13[ithread] = cy13;
		_cy14[ithread] = cy14;
		_cy15[ithread] = cy15;
		_cy16[ithread] = cy16;
		_cy17[ithread] = cy17;
		_cy18[ithread] = cy18;
		_cy19[ithread] = cy19;
		_cy20[ithread] = cy20;
		_cy21[ithread] = cy21;
		_cy22[ithread] = cy22;
		_cy23[ithread] = cy23;
		_cy24[ithread] = cy24;
		_cy25[ithread] = cy25;
		_cy26[ithread] = cy26;
		_cy27[ithread] = cy27;
		_cy28[ithread] = cy28;
		_cy29[ithread] = cy29;
		_cy30[ithread] = cy30;
		_cy31[ithread] = cy31;
		_cy32[ithread] = cy32;
		_cy33[ithread] = cy33;
		_cy34[ithread] = cy34;
		_cy35[ithread] = cy35;
		_cy36[ithread] = cy36;
		_cy37[ithread] = cy37;
		_cy38[ithread] = cy38;
		_cy39[ithread] = cy39;
		_cy40[ithread] = cy40;
		_cy41[ithread] = cy41;
		_cy42[ithread] = cy42;
		_cy43[ithread] = cy43;
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #ifdef OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(HERE, 0x0 == cy44_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix44_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}

		_cy00[ithread] = tdat[ithread].cy00;
		_cy01[ithread] = tdat[ithread].cy01;
		_cy02[ithread] = tdat[ithread].cy02;
		_cy03[ithread] = tdat[ithread].cy03;
		_cy04[ithread] = tdat[ithread].cy04;
		_cy05[ithread] = tdat[ithread].cy05;
		_cy06[ithread] = tdat[ithread].cy06;
		_cy07[ithread] = tdat[ithread].cy07;
		_cy08[ithread] = tdat[ithread].cy08;
		_cy09[ithread] = tdat[ithread].cy09;
		_cy10[ithread] = tdat[ithread].cy10;
		_cy11[ithread] = tdat[ithread].cy11;
		_cy12[ithread] = tdat[ithread].cy12;
		_cy13[ithread] = tdat[ithread].cy13;
		_cy14[ithread] = tdat[ithread].cy14;
		_cy15[ithread] = tdat[ithread].cy15;
		_cy16[ithread] = tdat[ithread].cy16;
		_cy17[ithread] = tdat[ithread].cy17;
		_cy18[ithread] = tdat[ithread].cy18;
		_cy19[ithread] = tdat[ithread].cy19;
		_cy20[ithread] = tdat[ithread].cy20;
		_cy21[ithread] = tdat[ithread].cy21;
		_cy22[ithread] = tdat[ithread].cy22;
		_cy23[ithread] = tdat[ithread].cy23;
		_cy24[ithread] = tdat[ithread].cy24;
		_cy25[ithread] = tdat[ithread].cy25;
		_cy26[ithread] = tdat[ithread].cy26;
		_cy27[ithread] = tdat[ithread].cy27;
		_cy28[ithread] = tdat[ithread].cy28;
		_cy29[ithread] = tdat[ithread].cy29;
		_cy30[ithread] = tdat[ithread].cy30;
		_cy31[ithread] = tdat[ithread].cy31;
		_cy32[ithread] = tdat[ithread].cy32;
		_cy33[ithread] = tdat[ithread].cy33;
		_cy34[ithread] = tdat[ithread].cy34;
		_cy35[ithread] = tdat[ithread].cy35;
		_cy36[ithread] = tdat[ithread].cy36;
		_cy37[ithread] = tdat[ithread].cy37;
		_cy38[ithread] = tdat[ithread].cy38;
		_cy39[ithread] = tdat[ithread].cy39;
		_cy40[ithread] = tdat[ithread].cy40;
		_cy41[ithread] = tdat[ithread].cy41;
		_cy42[ithread] = tdat[ithread].cy42;
		_cy43[ithread] = tdat[ithread].cy43;
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-44 forward DIF FFT of the first block of 44 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 44 outputs of (1);
	!   (3) Reweight and perform a radix-44 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 44 elements and repeat (1-4).
	*/
	a1p00r= _cy00[CY_THREADS - 1];	/* Note base-11 indexing of temps here */
	a1p01r= _cy01[CY_THREADS - 1];
	a1p02r= _cy02[CY_THREADS - 1];
	a1p03r= _cy03[CY_THREADS - 1];
	a1p04r= _cy04[CY_THREADS - 1];
	a1p05r= _cy05[CY_THREADS - 1];
	a1p06r= _cy06[CY_THREADS - 1];
	a1p07r= _cy07[CY_THREADS - 1];
	a1p08r= _cy08[CY_THREADS - 1];
	a1p09r= _cy09[CY_THREADS - 1];
	a1p10r= _cy10[CY_THREADS - 1];
	a1p11r= _cy11[CY_THREADS - 1];
	a1p12r= _cy12[CY_THREADS - 1];
	a1p13r= _cy13[CY_THREADS - 1];
	a1p14r= _cy14[CY_THREADS - 1];
	a1p15r= _cy15[CY_THREADS - 1];
	a1p16r= _cy16[CY_THREADS - 1];
	a1p17r= _cy17[CY_THREADS - 1];
	a1p18r= _cy18[CY_THREADS - 1];
	a1p19r= _cy19[CY_THREADS - 1];
	a1p20r= _cy20[CY_THREADS - 1];
	a1p21r= _cy21[CY_THREADS - 1];
	a1p22r= _cy22[CY_THREADS - 1];
	a1p23r= _cy23[CY_THREADS - 1];
	a1p24r= _cy24[CY_THREADS - 1];
	a1p25r= _cy25[CY_THREADS - 1];
	a1p26r= _cy26[CY_THREADS - 1];
	a1p27r= _cy27[CY_THREADS - 1];
	a1p28r= _cy28[CY_THREADS - 1];
	a1p29r= _cy29[CY_THREADS - 1];
	a1p30r= _cy30[CY_THREADS - 1];
	a1p31r= _cy31[CY_THREADS - 1];
	a1p32r= _cy32[CY_THREADS - 1];
	a1p33r= _cy33[CY_THREADS - 1];
	a1p34r= _cy34[CY_THREADS - 1];
	a1p35r= _cy35[CY_THREADS - 1];
	a1p36r= _cy36[CY_THREADS - 1];
	a1p37r= _cy37[CY_THREADS - 1];
	a1p38r= _cy38[CY_THREADS - 1];
	a1p39r= _cy39[CY_THREADS - 1];
	a1p40r= _cy40[CY_THREADS - 1];
	a1p41r= _cy41[CY_THREADS - 1];
	a1p42r= _cy42[CY_THREADS - 1];
	a1p43r= _cy43[CY_THREADS - 1];

	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		ASSERT(HERE, CY_THREADS > 1,"radix44_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
		_cy00[ithread] = _cy00[ithread-1];
		_cy01[ithread] = _cy01[ithread-1];
		_cy02[ithread] = _cy02[ithread-1];
		_cy03[ithread] = _cy03[ithread-1];
		_cy04[ithread] = _cy04[ithread-1];
		_cy05[ithread] = _cy05[ithread-1];
		_cy06[ithread] = _cy06[ithread-1];
		_cy07[ithread] = _cy07[ithread-1];
		_cy08[ithread] = _cy08[ithread-1];
		_cy09[ithread] = _cy09[ithread-1];
		_cy10[ithread] = _cy10[ithread-1];
		_cy11[ithread] = _cy11[ithread-1];
		_cy12[ithread] = _cy12[ithread-1];
		_cy13[ithread] = _cy13[ithread-1];
		_cy14[ithread] = _cy14[ithread-1];
		_cy15[ithread] = _cy15[ithread-1];
		_cy16[ithread] = _cy16[ithread-1];
		_cy17[ithread] = _cy17[ithread-1];
		_cy18[ithread] = _cy18[ithread-1];
		_cy19[ithread] = _cy19[ithread-1];
		_cy20[ithread] = _cy20[ithread-1];
		_cy21[ithread] = _cy21[ithread-1];
		_cy22[ithread] = _cy22[ithread-1];
		_cy23[ithread] = _cy23[ithread-1];
		_cy24[ithread] = _cy24[ithread-1];
		_cy25[ithread] = _cy25[ithread-1];
		_cy26[ithread] = _cy26[ithread-1];
		_cy27[ithread] = _cy27[ithread-1];
		_cy28[ithread] = _cy28[ithread-1];
		_cy29[ithread] = _cy29[ithread-1];
		_cy30[ithread] = _cy30[ithread-1];
		_cy31[ithread] = _cy31[ithread-1];
		_cy32[ithread] = _cy32[ithread-1];
		_cy33[ithread] = _cy33[ithread-1];
		_cy34[ithread] = _cy34[ithread-1];
		_cy35[ithread] = _cy35[ithread-1];
		_cy36[ithread] = _cy36[ithread-1];
		_cy37[ithread] = _cy37[ithread-1];
		_cy38[ithread] = _cy38[ithread-1];
		_cy39[ithread] = _cy39[ithread-1];
		_cy40[ithread] = _cy40[ithread-1];
		_cy41[ithread] = _cy41[ithread-1];
		_cy42[ithread] = _cy42[ithread-1];
		_cy43[ithread] = _cy43[ithread-1];
	}

	_cy00[0] =+a1p43r;	/* ...The wraparound carry is here: */
	_cy01[0] = a1p00r;
	_cy02[0] = a1p01r;
	_cy03[0] = a1p02r;
	_cy04[0] = a1p03r;
	_cy05[0] = a1p04r;
	_cy06[0] = a1p05r;
	_cy07[0] = a1p06r;
	_cy08[0] = a1p07r;
	_cy09[0] = a1p08r;
	_cy10[0] = a1p09r;
	_cy11[0] = a1p10r;
	_cy12[0] = a1p11r;
	_cy13[0] = a1p12r;
	_cy14[0] = a1p13r;
	_cy15[0] = a1p14r;
	_cy16[0] = a1p15r;
	_cy17[0] = a1p16r;
	_cy18[0] = a1p17r;
	_cy19[0] = a1p18r;
	_cy20[0] = a1p19r;
	_cy21[0] = a1p20r;
	_cy22[0] = a1p21r;
	_cy23[0] = a1p22r;
	_cy24[0] = a1p23r;
	_cy25[0] = a1p24r;
	_cy26[0] = a1p25r;
	_cy27[0] = a1p26r;
	_cy28[0] = a1p27r;
	_cy29[0] = a1p28r;
	_cy30[0] = a1p29r;
	_cy31[0] = a1p30r;
	_cy32[0] = a1p31r;
	_cy33[0] = a1p32r;
	_cy34[0] = a1p33r;
	_cy35[0] = a1p34r;
	_cy36[0] = a1p35r;
	_cy37[0] = a1p36r;
	_cy38[0] = a1p37r;
	_cy39[0] = a1p38r;
	_cy40[0] = a1p39r;
	_cy41[0] = a1p40r;
	_cy42[0] = a1p41r;
	_cy43[0] = a1p42r;

	full_pass = 0;
	scale = 1;
	j_jhi = 7;

    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			jt = j;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p04;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p08;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p12;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p16;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p20;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p24;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p28;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p32;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p36;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p40;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
		}
	}
}	/* endfor(outer) */

    a1p00r = 0;
    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		a1p00r += fabs(_cy00[0])+fabs(_cy01[0])+fabs(_cy02[0])+fabs(_cy03[0])+fabs(_cy04[0])+fabs(_cy05[0])+fabs(_cy06[0])+fabs(_cy07[0])+fabs(_cy08[0])+fabs(_cy09[0])+fabs(_cy10[0])+fabs(_cy11[0])+fabs(_cy12[0])+fabs(_cy13[0])+fabs(_cy14[0])+fabs(_cy15[0])+fabs(_cy16[0])+fabs(_cy17[0])+fabs(_cy18[0])+fabs(_cy19[0]);
		a1p00r += fabs(_cy20[0])+fabs(_cy21[0])+fabs(_cy22[0])+fabs(_cy23[0])+fabs(_cy24[0])+fabs(_cy25[0])+fabs(_cy26[0])+fabs(_cy27[0])+fabs(_cy28[0])+fabs(_cy29[0])+fabs(_cy30[0])+fabs(_cy31[0])+fabs(_cy32[0])+fabs(_cy33[0])+fabs(_cy34[0])+fabs(_cy35[0])+fabs(_cy36[0])+fabs(_cy37[0])+fabs(_cy38[0])+fabs(_cy39[0]);
		a1p00r += fabs(_cy40[0])+fabs(_cy41[0])+fabs(_cy42[0])+fabs(_cy43[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
    }

	if(a1p00r != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix44_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	return(0);
}

/***************/

void radix44_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-44 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,jt,jp,j1,j2;
	static int n44,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40;
	static double 	a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1*/

			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
	double rt,it,
		t00r,t01r,t02r,t03r,t04r,t05r,t06r,t07r,t08r,t09r,t0ar,t00i,t01i,t02i,t03i,t04i,t05i,t06i,t07i,t08i,t09i,t0ai,
		t10r,t11r,t12r,t13r,t14r,t15r,t16r,t17r,t18r,t19r,t1ar,t10i,t11i,t12i,t13i,t14i,t15i,t16i,t17i,t18i,t19i,t1ai,
		t20r,t21r,t22r,t23r,t24r,t25r,t26r,t27r,t28r,t29r,t2ar,t20i,t21i,t22i,t23i,t24i,t25i,t26i,t27i,t28i,t29i,t2ai,
		t30r,t31r,t32r,t33r,t34r,t35r,t36r,t37r,t38r,t39r,t3ar,t30i,t31i,t32i,t33i,t34i,t35i,t36i,t37i,t38i,t39i,t3ai;

	if(!first_entry && (n/44) != n44)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n44=n/44;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n44;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-44 pass is here.	*/

    for(j=0; j < n44; j += 2)
    {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-9 transforms...*/
	/*
	Twiddleless version arranges 4 sets of radix-11 DFT inputs as follows: 0 in upper left corner, decrement 4 horizontally and 11 vertically:

		RADIX_11_DFT(00,40,36,32,28,24,20,16,12,08,04)
		RADIX_11_DFT(33,29,25,21,17,13,09,05,01,41,37)
		RADIX_11_DFT(22,18,14,10,06,02,42,38,34,30,26)
		RADIX_11_DFT(11,07,03,43,39,35,31,27,23,19,15)

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-4 DFT outputs.
	*/
		RADIX_11_DFT(a[j1    ],a[j2    ],a[j1+p40],a[j2+p40],a[j1+p36],a[j2+p36],a[j1+p32],a[j2+p32],a[j1+p28],a[j2+p28],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04],t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai);	jt = j1+p01; jp = j2+p01;
		RADIX_11_DFT(a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai);	jt = j1+p02; jp = j2+p02;
		RADIX_11_DFT(a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai);	jt = j1+p03; jp = j2+p03;
		RADIX_11_DFT(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai);

		/*...and now do 11 radix-4 transforms...*/
					/*          inputs                    */ /*                                      outputs                              */
		RADIX_04_DIF(t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_04_DIF(t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p36; jp = j2+p36;
		RADIX_04_DIF(t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p32; jp = j2+p32;
		RADIX_04_DIF(t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_04_DIF(t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p24; jp = j2+p24;
		RADIX_04_DIF(t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_04_DIF(t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_04_DIF(t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p12; jp = j2+p12;
		RADIX_04_DIF(t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_04_DIF(t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_04_DIF(t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

		/* Totals: 4*radix11 + 11*radix04 = 4*(152 FADD, 44 FMUL)	+ 11*(16 FADD, 0 FMUL) = 784 FADD, 176 FMUL	*/
	}
}
/*
Unpermuted DIF outputs in left column, and the output index where each needs to go in []
in order to match the bit-reversed matrix-multiply DFT outputs in the right set of columns:

   0 [   0] 236.0000000000   189.0000000000  (   0)  236.0000000000   189.0000000000
   1 [   1] -32.0000000000   -27.0000000000  (  22)  -32.0000000000   -27.0000000000
   2 [   2] -10.0000000000    -5.0000000000  (  11)  -10.0000000000    -5.0000000000
   3 [   3] -18.0000000000   -25.0000000000  (  33)  -18.0000000000   -25.0000000000

   4 [  41]  11.1512215466    24.0242163196  (   1)   -3.5575051052    -7.5004478474
   5 [  40]  -3.6545658269    18.0551079573  (  23)    3.1041298533    -7.8963424647
   6 [  43]   7.4820065980   -12.0231995596  (  12)   14.9388279007   -12.3752085181
   7 [  42]  16.9375513781    23.8501660484  (  34)    8.8463140932   -16.0678583967

   8 [  38]  39.0936733850   -32.3447548001  (   2)   -9.4839628566     6.8049779883
   9 [  39]  -0.5216089174    15.2315827674  (  24)   29.2378241518    33.2630685362
  10 [  37] -48.7357676843    10.4812779061  (  13)   22.0526961791    -6.2994315330
  11 [  36]  -3.7958386977    14.7194188005  (  35)  -12.8767451770    -6.2427936444

  12 [  32] -10.8100095734    21.3984549366  (   3)  -19.1813630782   -13.9272285822
  13 [  33] -19.1224413133   -21.5236610012  (  25)  -10.4370837863   -26.2672878420
  14 [  34]   1.1660991398    28.0446190599  (  14)   28.5057795698   -23.4048812091
  15 [  35]  -1.1477307972    -3.9728027190  (  36)   13.2967408741    10.4272914314

  16 [  31]  -4.2498272043    10.6726404285  (   4)  -15.5870869025   -44.7914305956
  17 [  30]  -0.9626434619     0.5229653520  (  26)   -1.0423591098    40.9542963604
  18 [  28] -25.0400889167   -18.4141433297  (  15)    6.7261054808    -4.8184504585
  19 [  29] -28.7125623109   -32.1535503331  (  37)   17.0043942516    -9.8950027594

  20 [  25]  -4.1708965255     8.4490527718  (   5)  -24.4303401948   -18.0065547227
  21 [  24]  12.2156779504   -12.9675770627  (  27)    1.4536458478   -23.1558811910
  22 [  27] -14.0980220633   -24.0478577359  (  16)  -33.9004676525    34.2766694898
  23 [  26]  12.0864190834    -8.4745070923  (  38)    6.2198098725    -2.6049527558

  24 [  22] -33.9004676525    34.2766694898  (   6)   12.2156779504   -12.9675770627
  25 [  23]   6.2198098725    -2.6049527558  (  28)   -4.1708965255     8.4490527718
  26 [  21]   1.4536458478   -23.1558811910  (  17)   12.0864190834    -8.4745070923
  27 [  20] -24.4303401948   -18.0065547227  (  39)  -14.0980220633   -24.0478577359

  28 [  16] -15.5870869025   -44.7914305956  (   7)  -25.0400889167   -18.4141433297
  29 [  17]  -1.0423591098    40.9542963604  (  29)  -28.7125623109   -32.1535503331
  30 [  18]   6.7261054808    -4.8184504585  (  18)   -0.9626434619     0.5229653520
  31 [  19]  17.0043942516    -9.8950027594  (  40)   -4.2498272043    10.6726404285

  32 [  15]  13.2967408741    10.4272914314  (   8)  -10.8100095734    21.3984549366
  33 [  14]  28.5057795698   -23.4048812091  (  30)  -19.1224413133   -21.5236610012
  34 [  12] -19.1813630782   -13.9272285822  (  19)    1.1660991398    28.0446190599
  35 [  13] -10.4370837863   -26.2672878420  (  41)   -1.1477307972    -3.9728027190

  36 [   9]  29.2378241518    33.2630685362  (   9)   -3.7958386977    14.7194188005
  37 [   8]  -9.4839628566     6.8049779883  (  31)  -48.7357676843    10.4812779061
  38 [  11] -12.8767451770    -6.2427936444  (  20)   39.0936733850   -32.3447548001
  39 [  10]  22.0526961791    -6.2994315330  (  42)   -0.5216089174    15.2315827674

  40 [   6]  14.9388279007   -12.3752085181  (  10)   -3.6545658269    18.0551079573
  41 [   7]   8.8463140932   -16.0678583967  (  32)   11.1512215466    24.0242163196
  42 [   5]   3.1041298533    -7.8963424647  (  21)   16.9375513781    23.8501660484
  43 [   4]  -3.5575051052    -7.5004478474  (  43)    7.4820065980   -12.0231995596
*/
/***************/

void radix44_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-44 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix44_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int j,jt,jp,j1,j2;
	static int n44,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40;
	static double 	a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1*/

			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
	double rt,it,
		t00r,t01r,t02r,t03r,t04r,t05r,t06r,t07r,t08r,t09r,t0ar,t00i,t01i,t02i,t03i,t04i,t05i,t06i,t07i,t08i,t09i,t0ai,
		t10r,t11r,t12r,t13r,t14r,t15r,t16r,t17r,t18r,t19r,t1ar,t10i,t11i,t12i,t13i,t14i,t15i,t16i,t17i,t18i,t19i,t1ai,
		t20r,t21r,t22r,t23r,t24r,t25r,t26r,t27r,t28r,t29r,t2ar,t20i,t21i,t22i,t23i,t24i,t25i,t26i,t27i,t28i,t29i,t2ai,
		t30r,t31r,t32r,t33r,t34r,t35r,t36r,t37r,t38r,t39r,t3ar,t30i,t31i,t32i,t33i,t34i,t35i,t36i,t37i,t38i,t39i,t3ai;

	if(!first_entry && (n/44) != n44)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n44=n/44;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n44;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-44 pass is here.	*/

    for(j=0; j < n44; j += 2)
    {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43
			  -> 0,40,36,32,28,24,20,16,12, 8, 4|33,29,25,21,17,13, 9, 5, 1,41,37|22,18,14,10, 6, 2,42,38,34,30,26|11, 7, 3,43,39,35,31,27,23,19,15

		I.e. start out with first quartet of indices {0,11,22,33}, permute those according to
		(0,11,22,33}*43%44 = {0,33,22,11), then each is head of a length-11 list of indices with decrement 4.

		Remember, inputs to DIT are bit-reversed, so
		a[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43] contain
		x[ 0,22,11,33, 1,23,12,34, 2,24,13,35, 3,25,14,36, 4,26,15,37, 5,27,16,38, 6,28,17,39, 7,29,18,40, 8,30,19,41, 9,31,20,42,10,32,21,43], which get swapped to
	(Look at first nontrivial one, i.e. x[1 -> 40] ... in terms of a[] this translates to a[4 -> 31])
		x[ 0,22,33,11,40,18,29, 7,36,14,25, 3,32,10,21,43,28, 6,17,39,24, 2,13,35,20,42, 9,31,16,38, 5,27,12,34, 1,23, 8,30,41,19, 4,26,37,15], which means the a-indices get swapped as
		a[ 0, 1, 3, 2,31,30,29,28,15,14,13,12,41,40,42,43,25,24,26,27, 9, 8,10,11,38,39,36,37,22,23,20,21, 6, 7, 4, 5,32,33,35,34,16,17,19,18].
	*/
	/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 11 radix-4 transforms...*/
					/*                                   inputs                                   */ /*              outputs              */
		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,rt,it);	jt = j1+p12; jp = j2+p12;
		RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,rt,it);	jt = j1+p40; jp = j2+p40;
		RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,rt,it);	jt = j1+p24; jp = j2+p24;
		RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,rt,it);	jt = j1+p36; jp = j2+p36;
		RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,rt,it);	jt = j1+p04; jp = j2+p04;
		RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,rt,it);	jt = j1+p32; jp = j2+p32;
		RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,rt,it);	jt = j1+p16; jp = j2+p16;
		RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,rt,it);

		/*...and now do 4 radix-11 transforms...*/
					/*                                                  inputs                                                  */ /*                                                                                            outputs                                                                                                                    */
		jt = j1    ; jp = j2    ;	RADIX_11_DFT(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32]);
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40]);
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08]);
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20]);
	}
}
/*
Unpermuted DIF outputs in left column, and the output index where each needs to go in []
in order to match the matrix-multiply DFT outputs in the right set of columns:

   0 [00]   236.0000000000  -189.0000000000     236.0000000000  -189.0000000000
   1 [11]   -10.0000000000     5.0000000000      -3.5575051052     7.5004478474
   2 [22]   -32.0000000000    27.0000000000      -9.4839628566    -6.8049779883
   3 [33]   -18.0000000000    25.0000000000     -19.1813630782    13.9272285822

   4 [12]    14.9388279007    12.3752085181     -15.5870869025    44.7914305956
   5 [23]     3.1041298533     7.8963424647     -24.4303401948    18.0065547227
   6 [34]     8.8463140932    16.0678583967      12.2156779504    12.9675770627
   7 [01]    -3.5575051052     7.5004478474     -25.0400889167    18.4141433297

   8 [24]    29.2378241518   -33.2630685362     -10.8100095734   -21.3984549366
   9 [35]   -12.8767451770     6.2427936444      -3.7958386977   -14.7194188005
  10 [02]    -9.4839628566    -6.8049779883      -3.6545658269   -18.0551079573
  11 [13]    22.0526961791     6.2994315330     -10.0000000000     5.0000000000

  12 [36]    13.2967408741   -10.4272914314      14.9388279007    12.3752085181
  13 [03]   -19.1813630782    13.9272285822      22.0526961791     6.2994315330
  14 [14]    28.5057795698    23.4048812091      28.5057795698    23.4048812091
  15 [25]   -10.4370837863    26.2672878420       6.7261054808     4.8184504585

  16 [04]   -15.5870869025    44.7914305956     -33.9004676525   -34.2766694898
  17 [15]     6.7261054808     4.8184504585      12.0864190834     8.4745070923
  18 [26]    -1.0423591098   -40.9542963604      -0.9626434619    -0.5229653520
  19 [37]    17.0043942516     9.8950027594       1.1660991398   -28.0446190599

  20 [16]   -33.9004676525   -34.2766694898      39.0936733850    32.3447548001
  21 [27]     1.4536458478    23.1558811910      16.9375513781   -23.8501660484
  22 [38]     6.2198098725     2.6049527558     -32.0000000000    27.0000000000
  23 [05]   -24.4303401948    18.0065547227       3.1041298533     7.8963424647

  24 [28]    -4.1708965255    -8.4490527718      29.2378241518   -33.2630685362
  25 [39]   -14.0980220633    24.0478577359     -10.4370837863    26.2672878420
  26 [06]    12.2156779504    12.9675770627      -1.0423591098   -40.9542963604
  27 [17]    12.0864190834     8.4745070923       1.4536458478    23.1558811910

  28 [40]    -4.2498272043   -10.6726404285      -4.1708965255    -8.4490527718
  29 [07]   -25.0400889167    18.4141433297     -28.7125623109    32.1535503331
  30 [18]    -0.9626434619    -0.5229653520     -19.1224413133    21.5236610012
  31 [29]   -28.7125623109    32.1535503331     -48.7357676843   -10.4812779061

  32 [08]   -10.8100095734   -21.3984549366      11.1512215466   -24.0242163196
  33 [19]     1.1660991398   -28.0446190599     -18.0000000000    25.0000000000
  34 [30]   -19.1224413133    21.5236610012       8.8463140932    16.0678583967
  35 [41]    -1.1477307972     3.9728027190     -12.8767451770     6.2427936444

  36 [20]    39.0936733850    32.3447548001      13.2967408741   -10.4272914314
  37 [31]   -48.7357676843   -10.4812779061      17.0043942516     9.8950027594
  38 [42]    -0.5216089174   -15.2315827674       6.2198098725     2.6049527558
  39 [09]    -3.7958386977   -14.7194188005     -14.0980220633    24.0478577359

  40 [32]    11.1512215466   -24.0242163196      -4.2498272043   -10.6726404285
  41 [43]     7.4820065980    12.0231995596      -1.1477307972     3.9728027190
  42 [10]    -3.6545658269   -18.0551079573      -0.5216089174   -15.2315827674
  43 [21]    16.9375513781   -23.8501660484       7.4820065980    12.0231995596

Now because we output things such that the 0,1,2,3 index offsets are the *columns* of of 4 x 9 set of out-array-indices,
this means that the output permutation translates (in terms of of 4 radix-11 macro outputs above) translates into the following:

	00,12,24,36,04,16,28,40,08,20,32
	11,23,35,03,15,27,39,07,19,31,43
	22,34,02,14,26,38,06,18,30,42,10
	33,01,13,25,37,05,17,29,41,09,21
*/

/******************** Multithreaded function body: ***************************/

#ifdef USE_PTHREAD

	#ifndef USE_SSE2
		#error pthreaded carry code requires SSE2-enabled build!
	#endif
	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy44_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 44;
		const double crnd = 3.0*0x4000000*0x2000000;
		int j,j1,j2,k;
		int l,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40;
		double *add0, *add1, *add2, *add3;
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09
			,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19
			,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29
			,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39
			,*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43;
		struct complex *two,*five, *ua0,*ua1,*ua2,*ua3,*ua4,*ua5,*ua6,*ua7,*ua8,*ua9, *ub0,*ub1,*ub2,*ub3,*ub4,*ub5,*ub6,*ub7,*ub8,*ub9, *max_err, *sse2_rnd, *half_arr, *tmp,*tmp2;
		struct complex
		 *t00r,*t01r,*t02r,*t03r,*t04r,*t05r,*t06r,*t07r,*t08r,*t09r,*t0ar,*t00i,*t01i,*t02i,*t03i,*t04i,*t05i,*t06i,*t07i,*t08i,*t09i,*t0ai
		,*t10r,*t11r,*t12r,*t13r,*t14r,*t15r,*t16r,*t17r,*t18r,*t19r,*t1ar,*t10i,*t11i,*t12i,*t13i,*t14i,*t15i,*t16i,*t17i,*t18i,*t19i,*t1ai
		,*t20r,*t21r,*t22r,*t23r,*t24r,*t25r,*t26r,*t27r,*t28r,*t29r,*t2ar,*t20i,*t21i,*t22i,*t23i,*t24i,*t25i,*t26i,*t27i,*t28i,*t29i,*t2ai
		,*t30r,*t31r,*t32r,*t33r,*t34r,*t35r,*t36r,*t37r,*t38r,*t39r,*t3ar,*t30i,*t31i,*t32i,*t33i,*t34i,*t35i,*t36i,*t37i,*t38i,*t39i,*t3ai
		,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p40r,*s1p41r,*s1p42r,*s1p43r
		,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i,*s1p20i,*s1p21i,*s1p22i,*s1p23i,*s1p24i,*s1p25i,*s1p26i,*s1p27i,*s1p28i,*s1p29i,*s1p30i,*s1p31i,*s1p32i,*s1p33i,*s1p34i,*s1p35i,*s1p36i,*s1p37i,*s1p38i,*s1p39i,*s1p40i,*s1p41i,*s1p42i,*s1p43i;
		struct complex *cy00,*cy02,*cy04,*cy06,*cy08,*cy10,*cy12,*cy14,*cy16,*cy18,*cy20,*cy22,*cy24,*cy26,*cy28,*cy30,*cy32,*cy34,*cy36,*cy38,*cy40,*cy42;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int khi    = thread_arg->khi;
		int i      = thread_arg->i;	/* Pointer to the BASE and BASEINV arrays.	*/
		int jstart = thread_arg->jstart;
		int jhi    = thread_arg->jhi;
		int col = thread_arg->col;
		int co2 = thread_arg->co2;
		int co3 = thread_arg->co3;
		int sw  = thread_arg->sw;
		int nwt = thread_arg->nwt;

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );

		s1p00r	= thread_arg->s1p00r;
								tmp = s1p00r + 0x58;	tmp2 = tmp + 0x58;
		s1p00r = s1p00r + 0x00;	t00r = tmp + 0x00;	cy00  = tmp2+ 0x00;
		s1p00i = s1p00r + 0x01;	t00i = tmp + 0x01;	cy02  = tmp2+ 0x01;
		s1p01r = s1p00r + 0x02;	t01r = tmp + 0x02;	cy04  = tmp2+ 0x02;
		s1p01i = s1p00r + 0x03;	t01i = tmp + 0x03;	cy06  = tmp2+ 0x03;
		s1p02r = s1p00r + 0x04;	t02r = tmp + 0x04;	cy08  = tmp2+ 0x04;
		s1p02i = s1p00r + 0x05;	t02i = tmp + 0x05;	cy10  = tmp2+ 0x05;
		s1p03r = s1p00r + 0x06;	t03r = tmp + 0x06;	cy12  = tmp2+ 0x06;
		s1p03i = s1p00r + 0x07;	t03i = tmp + 0x07;	cy14  = tmp2+ 0x07;
		s1p04r = s1p00r + 0x08;	t04r = tmp + 0x08;	cy16  = tmp2+ 0x08;
		s1p04i = s1p00r + 0x09;	t04i = tmp + 0x09;	cy18  = tmp2+ 0x09;
		s1p05r = s1p00r + 0x0a;	t05r = tmp + 0x0a;	cy20  = tmp2+ 0x0a;
		s1p05i = s1p00r + 0x0b;	t05i = tmp + 0x0b;	cy22  = tmp2+ 0x0b;
		s1p06r = s1p00r + 0x0c;	t06r = tmp + 0x0c;	cy24  = tmp2+ 0x0c;
		s1p06i = s1p00r + 0x0d;	t06i = tmp + 0x0d;	cy26  = tmp2+ 0x0d;
		s1p07r = s1p00r + 0x0e;	t07r = tmp + 0x0e;	cy28  = tmp2+ 0x0e;
		s1p07i = s1p00r + 0x0f;	t07i = tmp + 0x0f;	cy30  = tmp2+ 0x0f;
		s1p08r = s1p00r + 0x10;	t08r = tmp + 0x10;	cy32  = tmp2+ 0x10;
		s1p08i = s1p00r + 0x11;	t08i = tmp + 0x11;	cy34  = tmp2+ 0x11;
		s1p09r = s1p00r + 0x12;	t09r = tmp + 0x12;	cy36  = tmp2+ 0x12;
		s1p09i = s1p00r + 0x13;	t09i = tmp + 0x13;	cy38  = tmp2+ 0x13;
		s1p10r = s1p00r + 0x14;	t0ar = tmp + 0x14;	cy40  = tmp2+ 0x14;
		s1p10i = s1p00r + 0x15;	t0ai = tmp + 0x15;	cy42  = tmp2+ 0x15;
		s1p11r = s1p00r + 0x16;	t10r = tmp + 0x16;	two     = tmp2+ 0x16;
		s1p11i = s1p00r + 0x17;	t10i = tmp + 0x17;	five    = tmp2+ 0x17;
		s1p12r = s1p00r + 0x18;	t11r = tmp + 0x18;	ua0     = tmp2+ 0x18;
		s1p12i = s1p00r + 0x19;	t11i = tmp + 0x19;	ua1     = tmp2+ 0x19;
		s1p13r = s1p00r + 0x1a;	t12r = tmp + 0x1a;	ua2     = tmp2+ 0x1a;
		s1p13i = s1p00r + 0x1b;	t12i = tmp + 0x1b;	ua3     = tmp2+ 0x1b;
		s1p14r = s1p00r + 0x1c;	t13r = tmp + 0x1c;	ua4     = tmp2+ 0x1c;
		s1p14i = s1p00r + 0x1d;	t13i = tmp + 0x1d;	ua5     = tmp2+ 0x1d;
		s1p15r = s1p00r + 0x1e;	t14r = tmp + 0x1e;	ua6     = tmp2+ 0x1e;
		s1p15i = s1p00r + 0x1f;	t14i = tmp + 0x1f;	ua7     = tmp2+ 0x1f;
		s1p16r = s1p00r + 0x20;	t15r = tmp + 0x20;	ua8     = tmp2+ 0x20;
		s1p16i = s1p00r + 0x21;	t15i = tmp + 0x21;	ua9     = tmp2+ 0x21;
		s1p17r = s1p00r + 0x22;	t16r = tmp + 0x22;	ub0     = tmp2+ 0x22;
		s1p17i = s1p00r + 0x23;	t16i = tmp + 0x23;	ub1     = tmp2+ 0x23;
		s1p18r = s1p00r + 0x24;	t17r = tmp + 0x24;	ub2     = tmp2+ 0x24;
		s1p18i = s1p00r + 0x25;	t17i = tmp + 0x25;	ub3     = tmp2+ 0x25;
		s1p19r = s1p00r + 0x26;	t18r = tmp + 0x26;	ub4     = tmp2+ 0x26;
		s1p19i = s1p00r + 0x27;	t18i = tmp + 0x27;	ub5     = tmp2+ 0x27;
		s1p20r = s1p00r + 0x28;	t19r = tmp + 0x28;	ub6     = tmp2+ 0x28;
		s1p20i = s1p00r + 0x29;	t19i = tmp + 0x29;	ub7     = tmp2+ 0x29;
		s1p21r = s1p00r + 0x2a;	t1ar = tmp + 0x2a;	ub8     = tmp2+ 0x2a;
		s1p21i = s1p00r + 0x2b;	t1ai = tmp + 0x2b;	ub9     = tmp2+ 0x2b;
		s1p22r = s1p00r + 0x2c;	t20r = tmp + 0x2c;	max_err = tmp2+ 0x2c;
		s1p22i = s1p00r + 0x2d;	t20i = tmp + 0x2d;	sse2_rnd= tmp2+ 0x2d;
		s1p23r = s1p00r + 0x2e;	t21r = tmp + 0x2e;	half_arr= tmp2+ 0x2e;	/* This table needs 20x16 bytes */
		s1p23i = s1p00r + 0x2f;	t21i = tmp + 0x2f;
		s1p24r = s1p00r + 0x30;	t22r = tmp + 0x30;
		s1p24i = s1p00r + 0x31;	t22i = tmp + 0x31;
		s1p25r = s1p00r + 0x32;	t23r = tmp + 0x32;
		s1p25i = s1p00r + 0x33;	t23i = tmp + 0x33;
		s1p26r = s1p00r + 0x34;	t24r = tmp + 0x34;
		s1p26i = s1p00r + 0x35;	t24i = tmp + 0x35;
		s1p27r = s1p00r + 0x36;	t25r = tmp + 0x36;
		s1p27i = s1p00r + 0x37;	t25i = tmp + 0x37;
		s1p28r = s1p00r + 0x38;	t26r = tmp + 0x38;
		s1p28i = s1p00r + 0x39;	t26i = tmp + 0x39;
		s1p29r = s1p00r + 0x3a;	t27r = tmp + 0x3a;
		s1p29i = s1p00r + 0x3b;	t27i = tmp + 0x3b;
		s1p30r = s1p00r + 0x3c;	t28r = tmp + 0x3c;
		s1p30i = s1p00r + 0x3d;	t28i = tmp + 0x3d;
		s1p31r = s1p00r + 0x3e;	t29r = tmp + 0x3e;
		s1p31i = s1p00r + 0x3f;	t29i = tmp + 0x3f;
		s1p32r = s1p00r + 0x40;	t2ar = tmp + 0x40;
		s1p32i = s1p00r + 0x41;	t2ai = tmp + 0x41;
		s1p33r = s1p00r + 0x42;	t30r = tmp + 0x42;
		s1p33i = s1p00r + 0x43;	t30i = tmp + 0x43;
		s1p34r = s1p00r + 0x44;	t31r = tmp + 0x44;
		s1p34i = s1p00r + 0x45;	t31i = tmp + 0x45;
		s1p35r = s1p00r + 0x46;	t32r = tmp + 0x46;
		s1p35i = s1p00r + 0x47;	t32i = tmp + 0x47;
		s1p36r = s1p00r + 0x48;	t33r = tmp + 0x48;
		s1p36i = s1p00r + 0x49;	t33i = tmp + 0x49;
		s1p37r = s1p00r + 0x4a;	t34r = tmp + 0x4a;
		s1p37i = s1p00r + 0x4b;	t34i = tmp + 0x4b;
		s1p38r = s1p00r + 0x4c;	t35r = tmp + 0x4c;
		s1p38i = s1p00r + 0x4d;	t35i = tmp + 0x4d;
		s1p39r = s1p00r + 0x4e;	t36r = tmp + 0x4e;
		s1p39i = s1p00r + 0x4f;	t36i = tmp + 0x4f;
		s1p40r = s1p00r + 0x50;	t37r = tmp + 0x50;
		s1p40i = s1p00r + 0x51;	t37i = tmp + 0x51;
		s1p41r = s1p00r + 0x52;	t38r = tmp + 0x52;
		s1p41i = s1p00r + 0x53;	t38i = tmp + 0x53;
		s1p42r = s1p00r + 0x54;	t39r = tmp + 0x54;
		s1p42i = s1p00r + 0x55;	t39i = tmp + 0x55;
		s1p43r = s1p00r + 0x56;	t3ar = tmp + 0x56;
		s1p43i = s1p00r + 0x57;	t3ai = tmp + 0x57;

		ASSERT(HERE, (s1p00r == thread_arg->s1p00r), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->re == crnd && sse2_rnd->im == crnd), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr+10)->re * (half_arr+14)->re == 1.0 && (half_arr+10)->im * (half_arr+14)->im == 1.0, "thread-local memcheck failed!");

		max_err->re = 0.0;	max_err->im = 0.0;

		sign_mask = (uint64*)(s1p00r + radix44_creals_in_local_store);
		sse_bw  = sign_mask + 2;
		sse_sw  = sign_mask + 4;
		sse_n   = sign_mask + 6;
		bjmodn00 = (int*)(sign_mask + 8);
		bjmodn01 = bjmodn00 +  1;
		bjmodn02 = bjmodn00 +  2;
		bjmodn03 = bjmodn00 +  3;
		bjmodn04 = bjmodn00 +  4;
		bjmodn05 = bjmodn00 +  5;
		bjmodn06 = bjmodn00 +  6;
		bjmodn07 = bjmodn00 +  7;
		bjmodn08 = bjmodn00 +  8;
		bjmodn09 = bjmodn00 +  9;
		bjmodn10 = bjmodn00 + 10;
		bjmodn11 = bjmodn00 + 11;
		bjmodn12 = bjmodn00 + 12;
		bjmodn13 = bjmodn00 + 13;
		bjmodn14 = bjmodn00 + 14;
		bjmodn15 = bjmodn00 + 15;
		bjmodn16 = bjmodn00 + 16;
		bjmodn17 = bjmodn00 + 17;
		bjmodn18 = bjmodn00 + 18;
		bjmodn19 = bjmodn00 + 19;
		bjmodn20 = bjmodn00 + 20;
		bjmodn21 = bjmodn00 + 21;
		bjmodn22 = bjmodn00 + 22;
		bjmodn23 = bjmodn00 + 23;
		bjmodn24 = bjmodn00 + 24;
		bjmodn25 = bjmodn00 + 25;
		bjmodn26 = bjmodn00 + 26;
		bjmodn27 = bjmodn00 + 27;
		bjmodn28 = bjmodn00 + 28;
		bjmodn29 = bjmodn00 + 29;
		bjmodn30 = bjmodn00 + 30;
		bjmodn31 = bjmodn00 + 31;
		bjmodn32 = bjmodn00 + 32;
		bjmodn33 = bjmodn00 + 33;
		bjmodn34 = bjmodn00 + 34;
		bjmodn35 = bjmodn00 + 35;
		bjmodn36 = bjmodn00 + 36;
		bjmodn37 = bjmodn00 + 37;
		bjmodn38 = bjmodn00 + 38;
		bjmodn39 = bjmodn00 + 39;
		bjmodn40 = bjmodn00 + 40;
		bjmodn41 = bjmodn00 + 41;
		bjmodn42 = bjmodn00 + 42;
		bjmodn43 = bjmodn00 + 43;

		/* Init DWT-indices: */	/* init carries	*/
		*bjmodn00 = thread_arg->bjmodn00;	cy00->re = thread_arg->cy00;
		*bjmodn01 = thread_arg->bjmodn01;	cy00->im = thread_arg->cy01;
		*bjmodn02 = thread_arg->bjmodn02;	cy02->re = thread_arg->cy02;
		*bjmodn03 = thread_arg->bjmodn03;	cy02->im = thread_arg->cy03;
		*bjmodn04 = thread_arg->bjmodn04;	cy04->re = thread_arg->cy04;
		*bjmodn05 = thread_arg->bjmodn05;	cy04->im = thread_arg->cy05;
		*bjmodn06 = thread_arg->bjmodn06;	cy06->re = thread_arg->cy06;
		*bjmodn07 = thread_arg->bjmodn07;	cy06->im = thread_arg->cy07;
		*bjmodn08 = thread_arg->bjmodn08;	cy08->re = thread_arg->cy08;
		*bjmodn09 = thread_arg->bjmodn09;	cy08->im = thread_arg->cy09;
		*bjmodn10 = thread_arg->bjmodn10;	cy10->re = thread_arg->cy10;
		*bjmodn11 = thread_arg->bjmodn11;	cy10->im = thread_arg->cy11;
		*bjmodn12 = thread_arg->bjmodn12;	cy12->re = thread_arg->cy12;
		*bjmodn13 = thread_arg->bjmodn13;	cy12->im = thread_arg->cy13;
		*bjmodn14 = thread_arg->bjmodn14;	cy14->re = thread_arg->cy14;
		*bjmodn15 = thread_arg->bjmodn15;	cy14->im = thread_arg->cy15;
		*bjmodn16 = thread_arg->bjmodn16;	cy16->re = thread_arg->cy16;
		*bjmodn17 = thread_arg->bjmodn17;	cy16->im = thread_arg->cy17;
		*bjmodn18 = thread_arg->bjmodn18;	cy18->re = thread_arg->cy18;
		*bjmodn19 = thread_arg->bjmodn19;	cy18->im = thread_arg->cy19;
		*bjmodn20 = thread_arg->bjmodn20;	cy20->re = thread_arg->cy20;
		*bjmodn21 = thread_arg->bjmodn21;	cy20->im = thread_arg->cy21;
		*bjmodn22 = thread_arg->bjmodn22;	cy22->re = thread_arg->cy22;
		*bjmodn23 = thread_arg->bjmodn23;	cy22->im = thread_arg->cy23;
		*bjmodn24 = thread_arg->bjmodn24;	cy24->re = thread_arg->cy24;
		*bjmodn25 = thread_arg->bjmodn25;	cy24->im = thread_arg->cy25;
		*bjmodn26 = thread_arg->bjmodn26;	cy26->re = thread_arg->cy26;
		*bjmodn27 = thread_arg->bjmodn27;	cy26->im = thread_arg->cy27;
		*bjmodn28 = thread_arg->bjmodn28;	cy28->re = thread_arg->cy28;
		*bjmodn29 = thread_arg->bjmodn29;	cy28->im = thread_arg->cy29;
		*bjmodn30 = thread_arg->bjmodn30;	cy30->re = thread_arg->cy30;
		*bjmodn31 = thread_arg->bjmodn31;	cy30->im = thread_arg->cy31;
		*bjmodn32 = thread_arg->bjmodn32;	cy32->re = thread_arg->cy32;
		*bjmodn33 = thread_arg->bjmodn33;	cy32->im = thread_arg->cy33;
		*bjmodn34 = thread_arg->bjmodn34;	cy34->re = thread_arg->cy34;
		*bjmodn35 = thread_arg->bjmodn35;	cy34->im = thread_arg->cy35;
		*bjmodn36 = thread_arg->bjmodn36;	cy36->re = thread_arg->cy36;
		*bjmodn37 = thread_arg->bjmodn37;	cy36->im = thread_arg->cy37;
		*bjmodn38 = thread_arg->bjmodn38;	cy38->re = thread_arg->cy38;
		*bjmodn39 = thread_arg->bjmodn39;	cy38->im = thread_arg->cy39;
		*bjmodn40 = thread_arg->bjmodn40;	cy40->re = thread_arg->cy40;
		*bjmodn41 = thread_arg->bjmodn41;	cy40->im = thread_arg->cy41;
		*bjmodn42 = thread_arg->bjmodn42;	cy42->re = thread_arg->cy42;
		*bjmodn43 = thread_arg->bjmodn43;	cy42->im = thread_arg->cy43;

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += 4)
			{
			/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
			*/
				j1 = (j & mask01) + br4[j&3];
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

			  #if GCC_ASM_FULL_INLINE	// GCC or SUNC implied
	
				add0 = &a[j1    ];
				SSE2_RADIX44_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,t00r, ua0, s1p00r);
	
			  #else
	
				/* Outputs in SSE2 modes are temps 2*11*16 = 22*16 = 0x160 bytes apart: */
				add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t00r, 0x160)
				add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t01r, 0x160)
				add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t02r, 0x160)
				add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t03r, 0x160)
				add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t04r, 0x160)
				add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t05r, 0x160)
				add2 = &a[j1+p36];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t06r, 0x160)
				add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t07r, 0x160)
				add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t08r, 0x160)
				add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t09r, 0x160)
				add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t0ar, 0x160)
	
				/* Radix-11 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 12 (12*32 bytes = 0x180) or XX -= 32 (-32*32 bytes = -0x400) between successive outputs: */
																									/*   a1p00r,a1p12r,a1p24r,a1p36r,a1p04r,a1p16r,a1p28r,a1p40r,a1p08r,a1p20r,a1p32r */
				SSE2_RADIX_11_DFT(t00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p00r, 0x180, 0x300, 0x480, 0x080, 0x200, 0x380, 0x500, 0x100, 0x280, 0x400)
																									/*   a1p11r,a1p23r,a1p35r,a1p03r,a1p15r,a1p27r,a1p39r,a1p07r,a1p19r,a1p31r,a1p43r */
				SSE2_RADIX_11_DFT(t10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p11r, 0x180, 0x300,-0x100, 0x080, 0x200, 0x380,-0x080, 0x100, 0x280, 0x400)
																									/*   a1p22r,a1p34r,a1p02r,a1p14r,a1p26r,a1p38r,a1p06r,a1p18r,a1p30r,a1p42r,a1p10r */
				SSE2_RADIX_11_DFT(t20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p22r, 0x180,-0x280,-0x100, 0x080, 0x200,-0x200,-0x080, 0x100, 0x280,-0x180)
																									/*   a1p33r,a1p01r,a1p13r,a1p25r,a1p37r,a1p05r,a1p17r,a1p29r,a1p41r,a1p09r,a1p21r */
				SSE2_RADIX_11_DFT(t30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140, ua0, s1p33r,-0x400,-0x280,-0x100, 0x080,-0x380,-0x200,-0x080, 0x100,-0x300,-0x180)
			#endif

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p28r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p32r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p36r,add1,add2,add3,cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p40r,add1,add2,add3,cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p40r,add1,add2,add3,cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				l= (j+2) & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p28r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p32r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p36r,add1,add2,     cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p40r,add1,add2,     cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p28r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p32r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p36r,add1,add2,     cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p40r,add1,add2,     cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

				/*...The radix-44 DIF pass is here:	*/
	
				/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-11 transforms...*/
				/* Radix-11 DFT inputs are (cyclic) with pXXr having XX += 40 (40*32 bytes = +0x500) or XX -= 4 (4*32 bytes = -0x080), outputs are adjacent 32-byte-separated temps: */
								/*a1p00r,a1p40r,a1p36r,a1p32r,a1p28r,a1p24r,a1p20r,a1p16r,a1p12r,a1p08r,a1p04r */
				SSE2_RADIX_11_DFT(s1p00r, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
								/*a1p33r,a1p29r,a1p25r,a1p21r,a1p17r,a1p13r,a1p09r,a1p05r,a1p01r,a1p41r,a1p37r */
				SSE2_RADIX_11_DFT(s1p33r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300,-0x380,-0x400, 0x100, 0x080, ua0, t10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
								/*a1p22r,a1p18r,a1p14r,a1p10r,a1p06r,a1p02r,a1p42r,a1p38r,a1p34r,a1p30r,a1p26r */
				SSE2_RADIX_11_DFT(s1p22r,-0x080,-0x100,-0x180,-0x200,-0x280, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
								/*a1p11r,a1p07r,a1p03r,a1p43r,a1p39r,a1p35r,a1p31r,a1p27r,a1p23r,a1p19r,a1p15r */
				SSE2_RADIX_11_DFT(s1p11r,-0x080,-0x100, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
	
				/*...and now do 11 radix-4 transforms...*/
				/* Inputs in SSE2 modes are temps 2*11*16 = 22*16 = 0x160 bytes apart. Notice how the add* indices after the first row repeat with period 4: */
																																/*     outputs      */ /* inputs */
				add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t00r, 0x160)
				add1 = &a[j1+p40];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t01r, 0x160)
				add3 = &a[j1+p36];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t02r, 0x160)
				add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t03r, 0x160)
				add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t04r, 0x160)
				add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t05r, 0x160)
				add3 = &a[j1+p20];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t06r, 0x160)
				add0 = &a[j1+p16];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t07r, 0x160)
				add2 = &a[j1+p12];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t08r, 0x160)
				add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t09r, 0x160)
				add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, t0ar, 0x160)
	
			}	/* end for(j=_jstart; j < _jhi; j += 2) */

			jstart += nwt;
			jhi    += nwt;

			col += RADIX;
			co3 -= RADIX;
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		thread_arg->cy00 = cy00->re;
		thread_arg->cy01 = cy00->im;
		thread_arg->cy02 = cy02->re;
		thread_arg->cy03 = cy02->im;
		thread_arg->cy04 = cy04->re;
		thread_arg->cy05 = cy04->im;
		thread_arg->cy06 = cy06->re;
		thread_arg->cy07 = cy06->im;
		thread_arg->cy08 = cy08->re;
		thread_arg->cy09 = cy08->im;
		thread_arg->cy10 = cy10->re;
		thread_arg->cy11 = cy10->im;
		thread_arg->cy12 = cy12->re;
		thread_arg->cy13 = cy12->im;
		thread_arg->cy14 = cy14->re;
		thread_arg->cy15 = cy14->im;
		thread_arg->cy16 = cy16->re;
		thread_arg->cy17 = cy16->im;
		thread_arg->cy18 = cy18->re;
		thread_arg->cy19 = cy18->im;
		thread_arg->cy20 = cy20->re;
		thread_arg->cy21 = cy20->im;
		thread_arg->cy22 = cy22->re;
		thread_arg->cy23 = cy22->im;
		thread_arg->cy24 = cy24->re;
		thread_arg->cy25 = cy24->im;
		thread_arg->cy26 = cy26->re;
		thread_arg->cy27 = cy26->im;
		thread_arg->cy28 = cy28->re;
		thread_arg->cy29 = cy28->im;
		thread_arg->cy30 = cy30->re;
		thread_arg->cy31 = cy30->im;
		thread_arg->cy32 = cy32->re;
		thread_arg->cy33 = cy32->im;
		thread_arg->cy34 = cy34->re;
		thread_arg->cy35 = cy34->im;
		thread_arg->cy36 = cy36->re;
		thread_arg->cy37 = cy36->im;
		thread_arg->cy38 = cy38->re;
		thread_arg->cy39 = cy38->im;
		thread_arg->cy40 = cy40->re;
		thread_arg->cy41 = cy40->im;
		thread_arg->cy42 = cy42->re;
		thread_arg->cy43 = cy42->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}

		return 0x0;
	}
#endif

