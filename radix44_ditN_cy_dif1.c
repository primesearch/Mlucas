/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_MSVC)

		/*...Radix-11 DFT: Inputs in memory locations __I0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA],\
		where r0 is a memory address and the i's are literal [byte] offsets. Outputs similarly go into memory locations\
		__O0 + [__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA], assumed disjoint with inputs:\
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

		#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
		{\
		__asm__ volatile (\
			"movl	%[__tmp]   ,%%eax	\n\t"\
			"movl	%[__stride],%%esi	\n\t"\
			"movl	%%eax,%%ebx			\n\t"\
			"addl	%%esi,%%ebx			/* add_in1  */\n\t"\
			"shll	$1,%%esi			/* stride*2 */\n\t"\
			"movaps	    (%%eax),%%xmm0	\n\t"\
			"movaps	    (%%ebx),%%xmm2	\n\t"\
			"movaps	0x10(%%eax),%%xmm1	\n\t"\
			"movaps	0x10(%%ebx),%%xmm3	\n\t"\
			"movaps	    (%%eax),%%xmm4	\n\t"\
			"movaps	    (%%ebx),%%xmm6	\n\t"\
			"movaps	0x10(%%eax),%%xmm5	\n\t"\
			"movaps	0x10(%%ebx),%%xmm7	\n\t"\
			"addl	%%esi,%%eax			/* add_in2  */\n\t"\
			"addl	%%esi,%%ebx			/* add_in3  */\n\t"\
			"addpd	    (%%eax),%%xmm0	\n\t"\
			"addpd	    (%%ebx),%%xmm2	\n\t"\
			"addpd	0x10(%%eax),%%xmm1	\n\t"\
			"addpd	0x10(%%ebx),%%xmm3	\n\t"\
			"subpd	    (%%eax),%%xmm4	\n\t"\
			"subpd	    (%%ebx),%%xmm6	\n\t"\
			"subpd	0x10(%%eax),%%xmm5	\n\t"\
			"subpd	0x10(%%ebx),%%xmm7	\n\t"\
			"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
			"movl	%[__add0],%%eax		\n\t"\
			"movl	%[__add1],%%ebx		\n\t"\
			"movl	%[__add2],%%ecx		\n\t"\
			"movl	%[__add3],%%edx		\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm4		\n\t"\
			"subpd	%%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm5		\n\t"\
			"movaps	%%xmm0,     (%%ebx)	\n\t"\
			"movaps	%%xmm4,     (%%ecx)	\n\t"\
			"movaps	%%xmm1,0x010(%%ebx)	\n\t"\
			"movaps	%%xmm5,0x010(%%edx)	\n\t"\
			"addpd	%%xmm2,%%xmm2		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm3		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm2		\n\t"\
			"addpd	%%xmm4,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm3		\n\t"\
			"addpd	%%xmm5,%%xmm6		\n\t"\
			"movaps	%%xmm2,     (%%eax)	\n\t"\
			"movaps	%%xmm7,     (%%edx)	\n\t"\
			"movaps	%%xmm3,0x010(%%eax)	\n\t"\
			"movaps	%%xmm6,0x010(%%ecx)	\n\t"\
			:					/* outputs: none */\
			: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
			 ,[__add1] "m" (Xadd1)\
			 ,[__add2] "m" (Xadd2)\
			 ,[__add3] "m" (Xadd3)\
			 ,[__tmp] "m" (Xtmp)\
			 ,[__stride] "e" (Xstride)\
			: "eax","ebx","ecx","edx","rsi"		/* Clobbered registers */\
		);\
		}

		#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
		{\
		__asm__ volatile (\
			"movl	%[__add0],%%eax		\n\t"\
			"movl	%[__add1],%%ebx		\n\t"\
			"movl	%[__add2],%%ecx		\n\t"\
			"movl	%[__add3],%%edx		\n\t"\
			"movaps	    (%%eax),%%xmm0	\n\t"\
			"movaps	    (%%ecx),%%xmm4	\n\t"\
			"movaps	0x10(%%eax),%%xmm1	\n\t"\
			"movaps	0x10(%%ecx),%%xmm5	\n\t"\
			"movaps	%%xmm0,%%xmm2			\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm1,%%xmm3			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"movl	%[__tmp]   ,%%eax	\n\t"\
			"movl	%[__stride],%%ecx	\n\t"\
			"addpd	    (%%ebx),%%xmm0	\n\t"\
			"addpd	    (%%edx),%%xmm4	\n\t"\
			"addpd	0x10(%%ebx),%%xmm1	\n\t"\
			"addpd	0x10(%%edx),%%xmm5	\n\t"\
			"subpd	    (%%ebx),%%xmm2	\n\t"\
			"subpd	    (%%edx),%%xmm6	\n\t"\
			"subpd	0x10(%%ebx),%%xmm3	\n\t"\
			"subpd	0x10(%%edx),%%xmm7	\n\t"\
			"movl	%%eax,%%ebx			\n\t"\
			"addl	%%ecx,%%ebx			\n\t"\
			"movl	%%ebx,%%edx			\n\t"\
			"addl	%%ecx,%%ecx			\n\t"\
			"addl	%%ecx,%%edx			\n\t"\
			"addl	%%eax,%%ecx			\n\t"\
			"/* Finish radix-4 butterfly and store results into temp-array slots: */\n\t"\
			"subpd	%%xmm4,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm0,     (%%ecx)	\n\t"\
			"movaps	%%xmm2,     (%%edx)	\n\t"\
			"movaps	%%xmm1,0x010(%%ecx)	\n\t"\
			"movaps	%%xmm3,0x010(%%ebx)	\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm6			\n\t"\
			"movaps	%%xmm4,     (%%eax)	\n\t"\
			"movaps	%%xmm7,     (%%ebx)	\n\t"\
			"movaps	%%xmm5,0x010(%%eax)	\n\t"\
			"movaps	%%xmm6,0x010(%%edx)	\n\t"\
			:					/* outputs: none */\
			: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
			 ,[__add1] "m" (Xadd1)\
			 ,[__add2] "m" (Xadd2)\
			 ,[__add3] "m" (Xadd3)\
			 ,[__tmp] "m" (Xtmp)\
			 ,[__stride] "e" (Xstride)\
			: "eax","ebx","ecx","edx","esi"		/* Clobbered registers */\
		);\
		}

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
			: "eax","ebx","rcx"		/* Clobbered registers */\
			);\
		}

	  #else

		#define GCC_ASM_FULL_INLINE	1	// 0 to use small-macro form below, 1 to inline the fused macros as single big blob of asm (64-bit only)

		#if GCC_ASM_FULL_INLINE

			#include "radix44_ditN_cy_dif1_gcc64.h"

		#endif	// GCC_ASM_FULL_INLINE

		#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
		{\
		__asm__ volatile (\
			"movq	%[__tmp]   ,%%rax	\n\t"\
			"movq	%[__stride],%%rsi	\n\t"\
			"movq	%%rax,%%r15			\n\t"\
			"addq	%%rsi,%%r15			/* add_in1  */\n\t"\
			"shlq	$1,%%rsi			/* stride*2 */\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t"\
			"movaps	    (%%r15),%%xmm2	\n\t"\
			"movaps	0x10(%%rax),%%xmm1	\n\t"\
			"movaps	0x10(%%r15),%%xmm3	\n\t"\
			"movaps	    (%%rax),%%xmm4	\n\t"\
			"movaps	    (%%r15),%%xmm6	\n\t"\
			"movaps	0x10(%%rax),%%xmm5	\n\t"\
			"movaps	0x10(%%r15),%%xmm7	\n\t"\
			"addq	%%rsi,%%rax			/* add_in2  */\n\t"\
			"addq	%%rsi,%%r15			/* add_in3  */\n\t"\
			"addpd	    (%%rax),%%xmm0	\n\t"\
			"addpd	    (%%r15),%%xmm2	\n\t"\
			"addpd	0x10(%%rax),%%xmm1	\n\t"\
			"addpd	0x10(%%r15),%%xmm3	\n\t"\
			"subpd	    (%%rax),%%xmm4	\n\t"\
			"subpd	    (%%r15),%%xmm6	\n\t"\
			"subpd	0x10(%%rax),%%xmm5	\n\t"\
			"subpd	0x10(%%r15),%%xmm7	\n\t"\
			"/* Finish radix-4 butterfly and store results into main-array slots: */\n\t"\
			"movq	%[__add0],%%rax		\n\t"\
			"movq	%[__add1],%%r15		\n\t"\
			"movq	%[__add2],%%rcx		\n\t"\
			"movq	%[__add3],%%rdx		\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm4		\n\t"\
			"subpd	%%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm5		\n\t"\
			"movaps	%%xmm0,     (%%r15)	\n\t"\
			"movaps	%%xmm4,     (%%rcx)	\n\t"\
			"movaps	%%xmm1,0x010(%%r15)	\n\t"\
			"movaps	%%xmm5,0x010(%%rdx)	\n\t"\
			"addpd	%%xmm2,%%xmm2		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm3		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm2		\n\t"\
			"addpd	%%xmm4,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm3		\n\t"\
			"addpd	%%xmm5,%%xmm6		\n\t"\
			"movaps	%%xmm2,     (%%rax)	\n\t"\
			"movaps	%%xmm7,     (%%rdx)	\n\t"\
			"movaps	%%xmm3,0x010(%%rax)	\n\t"\
			"movaps	%%xmm6,0x010(%%rcx)	\n\t"\
			:					/* outputs: none */\
			: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
			 ,[__add1] "m" (Xadd1)\
			 ,[__add2] "m" (Xadd2)\
			 ,[__add3] "m" (Xadd3)\
			 ,[__tmp] "m" (Xtmp)\
			 ,[__stride] "e" (Xstride)\
			: "rax","r15","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
		}

		#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(Xadd0, Xadd1, Xadd2, Xadd3, Xtmp, Xstride)\
		{\
		__asm__ volatile (\
			"movq	%[__add0],%%rax		\n\t"\
			"movq	%[__add1],%%r15		\n\t"\
			"movq	%[__add2],%%rcx		\n\t"\
			"movq	%[__add3],%%rdx		\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t"\
			"movaps	    (%%rcx),%%xmm4	\n\t"\
			"movaps	0x10(%%rax),%%xmm1	\n\t"\
			"movaps	0x10(%%rcx),%%xmm5	\n\t"\
			"movaps	%%xmm0,%%xmm2			\n\t"\
			"movaps	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm1,%%xmm3			\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"movq	%[__tmp]   ,%%rax	\n\t"\
			"movq	%[__stride],%%rcx	\n\t"\
			"addpd	    (%%r15),%%xmm0	\n\t"\
			"addpd	    (%%rdx),%%xmm4	\n\t"\
			"addpd	0x10(%%r15),%%xmm1	\n\t"\
			"addpd	0x10(%%rdx),%%xmm5	\n\t"\
			"subpd	    (%%r15),%%xmm2	\n\t"\
			"subpd	    (%%rdx),%%xmm6	\n\t"\
			"subpd	0x10(%%r15),%%xmm3	\n\t"\
			"subpd	0x10(%%rdx),%%xmm7	\n\t"\
			"movq	%%rax,%%r15			\n\t"\
			"addq	%%rcx,%%r15			\n\t"\
			"movq	%%r15,%%rdx			\n\t"\
			"addq	%%rcx,%%rcx			\n\t"\
			"addq	%%rcx,%%rdx			\n\t"\
			"addq	%%rax,%%rcx			\n\t"\
			"/* Finish radix-4 butterfly and store results into temp-array slots: */\n\t"\
			"subpd	%%xmm4,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm0,     (%%rcx)	\n\t"\
			"movaps	%%xmm2,     (%%rdx)	\n\t"\
			"movaps	%%xmm1,0x010(%%rcx)	\n\t"\
			"movaps	%%xmm3,0x010(%%r15)	\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm2,%%xmm7			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm6			\n\t"\
			"movaps	%%xmm4,     (%%rax)	\n\t"\
			"movaps	%%xmm7,     (%%r15)	\n\t"\
			"movaps	%%xmm5,0x010(%%rax)	\n\t"\
			"movaps	%%xmm6,0x010(%%rdx)	\n\t"\
			:					/* outputs: none */\
			: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
			 ,[__add1] "m" (Xadd1)\
			 ,[__add2] "m" (Xadd2)\
			 ,[__add3] "m" (Xadd3)\
			 ,[__tmp] "m" (Xtmp)\
			 ,[__stride] "e" (Xstride)\
			: "rax","r15","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
		}

		#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
		{\
		__asm__ volatile (\
			"/********************************************/\n\t"\
			"/*       Here are the 5 cosine terms:       */\n\t"\
			"/********************************************/\n\t"\
			"movq	%[__I0],%%rax			\n\t"\
			"movq	%[__cc],%%r15			\n\t"\
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
			"mulpd	     (%%r15),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o1](%%rcx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0x10(%%r15),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%r15),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o3](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%r15),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o2](%%rcx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%r15),%%xmm4		\n\t"\
			"mulpd	 0x70(%%r15),%%xmm5		\n\t"\
			"mulpd	 0x10(%%r15),%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%r15),%%xmm6		\n\t"\
			"mulpd	 0x50(%%r15),%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%r15),%%xmm1		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	 0x90(%%r15),%%xmm2		\n\t"\
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
			"mulpd	     (%%r15),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o1](%%rcx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0x10(%%r15),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%r15),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o3](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%r15),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o2](%%rcx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%r15),%%xmm4		\n\t"\
			"mulpd	 0x70(%%r15),%%xmm5		\n\t"\
			"mulpd	 0x10(%%r15),%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%r15),%%xmm6		\n\t"\
			"mulpd	 0x50(%%r15),%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%r15),%%xmm1		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	 0x90(%%r15),%%xmm2		\n\t"\
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
			"addq	$0xa0,%%r15				\n\t"\
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
			"mulpd	     (%%r15),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o6](%%rcx)	\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0xb0(%%r15),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%r15),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o8](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%r15),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o7](%%rcx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%r15),%%xmm4		\n\t"\
			"mulpd	 0x70(%%r15),%%xmm5		\n\t"\
			"mulpd	 0x10(%%r15),%%xmm7		\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%r15),%%xmm6		\n\t"\
			"mulpd	 0x50(%%r15),%%xmm3		\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%r15),%%xmm1		\n\t"\
			"mulpd	 0x90(%%r15),%%xmm2		\n\t"\
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
			"mulpd	     (%%r15),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o6](%%rcx)	\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0xb0(%%r15),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%r15),%%xmm7		\n\t"\
			"movaps	%%xmm7,%c[__o8](%%rcx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%r15),%%xmm1		\n\t"\
			"movaps	%%xmm1,%c[__o7](%%rcx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%r15),%%xmm4		\n\t"\
			"mulpd	 0x70(%%r15),%%xmm5		\n\t"\
			"mulpd	 0x10(%%r15),%%xmm7		\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%r15),%%xmm6		\n\t"\
			"mulpd	 0x50(%%r15),%%xmm3		\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%r15),%%xmm1		\n\t"\
			"mulpd	 0x90(%%r15),%%xmm2		\n\t"\
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
			: "rax","r15","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
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
	int n44,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,outer;
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
	double a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,a1p36i,a1p37i,a1p38i,a1p39i,a1p40i,a1p41i,a1p42i,a1p43i;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;

#ifdef USE_SSE2

	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static struct complex *two,*five, *ua0,*ua1,*ua2,*ua3,*ua4,*ua5,*ua6,*ua7,*ua8,*ua9, *ub0,*ub1,*ub2,*ub3,*ub4,*ub5,*ub6,*ub7,*ub8,*ub9, *max_err, *sse2_rnd, *half_arr, *tmp,*tmp2;
	static struct complex
	 *t00r,*t01r,*t02r,*t03r,*t04r,*t05r,*t06r,*t07r,*t08r,*t09r,*t0ar,*t00i,*t01i,*t02i,*t03i,*t04i,*t05i,*t06i,*t07i,*t08i,*t09i,*t0ai
	,*t10r,*t11r,*t12r,*t13r,*t14r,*t15r,*t16r,*t17r,*t18r,*t19r,*t1ar,*t10i,*t11i,*t12i,*t13i,*t14i,*t15i,*t16i,*t17i,*t18i,*t19i,*t1ai
	,*t20r,*t21r,*t22r,*t23r,*t24r,*t25r,*t26r,*t27r,*t28r,*t29r,*t2ar,*t20i,*t21i,*t22i,*t23i,*t24i,*t25i,*t26i,*t27i,*t28i,*t29i,*t2ai
	,*t30r,*t31r,*t32r,*t33r,*t34r,*t35r,*t36r,*t37r,*t38r,*t39r,*t3ar,*t30i,*t31i,*t32i,*t33i,*t34i,*t35i,*t36i,*t37i,*t38i,*t39i,*t3ai
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p40r,*s1p41r,*s1p42r,*s1p43r
	,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i,*s1p20i,*s1p21i,*s1p22i,*s1p23i,*s1p24i,*s1p25i,*s1p26i,*s1p27i,*s1p28i,*s1p29i,*s1p30i,*s1p31i,*s1p32i,*s1p33i,*s1p34i,*s1p35i,*s1p36i,*s1p37i,*s1p38i,*s1p39i,*s1p40i,*s1p41i,*s1p42i,*s1p43i;
	static uint64 *sm_arr = 0x0, *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39,*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43;
	static struct complex *cy_r00,*cy_r02,*cy_r04,*cy_r06,*cy_r08,*cy_r10,*cy_r12,*cy_r14,*cy_r16,*cy_r18,*cy_r20,*cy_r22,*cy_r24,*cy_r26,*cy_r28,*cy_r30,*cy_r32,*cy_r34,*cy_r36,*cy_r38,*cy_r40,*cy_r42;

#else

	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double rt,it,temp,frac;
	double wt,wtinv,wtA,wtB,wtC;
	int m,m2;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21
		,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,bjmodn40,bjmodn41,bjmodn42,bjmodn43;
	double cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_r36,cy_r37,cy_r38,cy_r39,cy_r40,cy_r41,cy_r42,cy_r43
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
	*_cy_r00 = 0x0,*_cy_r01 = 0x0,*_cy_r02 = 0x0,*_cy_r03 = 0x0,*_cy_r04 = 0x0,*_cy_r05 = 0x0,*_cy_r06 = 0x0,*_cy_r07 = 0x0,*_cy_r08 = 0x0,*_cy_r09 = 0x0,*_cy_r10 = 0x0,*_cy_r11 = 0x0,*_cy_r12 = 0x0,*_cy_r13 = 0x0,*_cy_r14 = 0x0,*_cy_r15 = 0x0,*_cy_r16 = 0x0,*_cy_r17 = 0x0,*_cy_r18 = 0x0,*_cy_r19 = 0x0,*_cy_r20 = 0x0,*_cy_r21 = 0x0,*_cy_r22 = 0x0,*_cy_r23 = 0x0,*_cy_r24 = 0x0,*_cy_r25 = 0x0,*_cy_r26 = 0x0,*_cy_r27 = 0x0,*_cy_r28 = 0x0,*_cy_r29 = 0x0,*_cy_r30 = 0x0,*_cy_r31 = 0x0,*_cy_r32 = 0x0,*_cy_r33 = 0x0,*_cy_r34 = 0x0,*_cy_r35 = 0x0,*_cy_r36 = 0x0,*_cy_r37 = 0x0,*_cy_r38 = 0x0,*_cy_r39 = 0x0,*_cy_r40 = 0x0,*_cy_r41 = 0x0,*_cy_r42 = 0x0,*_cy_r43 = 0x0;

/*...change n44 and n_div_wt to non-static to work around a gcc compiler bug. */
	n44   = n/44;
	n_div_nwt = n44 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n44)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/44 in radix44_ditN_cy_dif1.\n",iter);
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
			CY_THREADS = MAX_THREADS;

		ASSERT(HERE, CY_THREADS >= NTHREADS,"radix44_ditN_cy_dif1.c: CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"radix44_ditN_cy_dif1.c: CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, n44      %CY_THREADS == 0,"radix44_ditN_cy_dif1.c: n44      %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"radix44_ditN_cy_dif1.c: n_div_nwt%CY_THREADS != 0");
		}

	#ifdef MULTITHREAD
		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
	#endif

		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)44));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef USE_SSE2

		ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE), "SSE2 currently only supports Mersenne-mod!");
		sc_arr = ALLOC_COMPLEX(sc_arr, 248);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		/* Size here is [8 + radix/2 + 4] 8-byte elements */
		sm_arr = ALLOC_UINT64(sm_arr, 34);	if(!sm_arr){ sprintf(cbuf, "FATAL: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sm_ptr = ALIGN_UINT64(sm_arr);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 88x2 16-byte slots of sc_arr for temporaries, next 21 for the constants needed by the radix-11 DFT,
	next 22 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
									tmp    = sc_ptr + 0x58;	tmp2 = tmp + 0x58;
		t00r= sc_ptr + 0x00;		s1p00r = tmp + 0x00;	cy_r00  = tmp2+ 0x00;
		t00i= sc_ptr + 0x01;		s1p00i = tmp + 0x01;	cy_r02  = tmp2+ 0x01;
		t01r= sc_ptr + 0x02;		s1p01r = tmp + 0x02;	cy_r04  = tmp2+ 0x02;
		t01i= sc_ptr + 0x03;		s1p01i = tmp + 0x03;	cy_r06  = tmp2+ 0x03;
		t02r= sc_ptr + 0x04;		s1p02r = tmp + 0x04;	cy_r08  = tmp2+ 0x04;
		t02i= sc_ptr + 0x05;		s1p02i = tmp + 0x05;	cy_r10  = tmp2+ 0x05;
		t03r= sc_ptr + 0x06;		s1p03r = tmp + 0x06;	cy_r12  = tmp2+ 0x06;
		t03i= sc_ptr + 0x07;		s1p03i = tmp + 0x07;	cy_r14  = tmp2+ 0x07;
		t04r= sc_ptr + 0x08;		s1p04r = tmp + 0x08;	cy_r16  = tmp2+ 0x08;
		t04i= sc_ptr + 0x09;		s1p04i = tmp + 0x09;	cy_r18  = tmp2+ 0x09;
		t05r= sc_ptr + 0x0a;		s1p05r = tmp + 0x0a;	cy_r20  = tmp2+ 0x0a;
		t05i= sc_ptr + 0x0b;		s1p05i = tmp + 0x0b;	cy_r22  = tmp2+ 0x0b;
		t06r= sc_ptr + 0x0c;		s1p06r = tmp + 0x0c;	cy_r24  = tmp2+ 0x0c;
		t06i= sc_ptr + 0x0d;		s1p06i = tmp + 0x0d;	cy_r26  = tmp2+ 0x0d;
		t07r= sc_ptr + 0x0e;		s1p07r = tmp + 0x0e;	cy_r28  = tmp2+ 0x0e;
		t07i= sc_ptr + 0x0f;		s1p07i = tmp + 0x0f;	cy_r30  = tmp2+ 0x0f;
		t08r= sc_ptr + 0x10;		s1p08r = tmp + 0x10;	cy_r32  = tmp2+ 0x10;
		t08i= sc_ptr + 0x11;		s1p08i = tmp + 0x11;	cy_r34  = tmp2+ 0x11;
		t09r= sc_ptr + 0x12;		s1p09r = tmp + 0x12;	cy_r36  = tmp2+ 0x12;
		t09i= sc_ptr + 0x13;		s1p09i = tmp + 0x13;	cy_r38  = tmp2+ 0x13;
		t0ar= sc_ptr + 0x14;		s1p10r = tmp + 0x14;	cy_r40  = tmp2+ 0x14;
		t0ai= sc_ptr + 0x15;		s1p10i = tmp + 0x15;	cy_r42  = tmp2+ 0x15;
		t10r= sc_ptr + 0x16;		s1p11r = tmp + 0x16;	two     = tmp2+ 0x16;
		t10i= sc_ptr + 0x17;		s1p11i = tmp + 0x17;	five    = tmp2+ 0x17;
		t11r= sc_ptr + 0x18;		s1p12r = tmp + 0x18;	ua0     = tmp2+ 0x18;
		t11i= sc_ptr + 0x19;		s1p12i = tmp + 0x19;	ua1     = tmp2+ 0x19;
		t12r= sc_ptr + 0x1a;		s1p13r = tmp + 0x1a;	ua2     = tmp2+ 0x1a;
		t12i= sc_ptr + 0x1b;		s1p13i = tmp + 0x1b;	ua3     = tmp2+ 0x1b;
		t13r= sc_ptr + 0x1c;		s1p14r = tmp + 0x1c;	ua4     = tmp2+ 0x1c;
		t13i= sc_ptr + 0x1d;		s1p14i = tmp + 0x1d;	ua5     = tmp2+ 0x1d;
		t14r= sc_ptr + 0x1e;		s1p15r = tmp + 0x1e;	ua6     = tmp2+ 0x1e;
		t14i= sc_ptr + 0x1f;		s1p15i = tmp + 0x1f;	ua7     = tmp2+ 0x1f;
		t15r= sc_ptr + 0x20;		s1p16r = tmp + 0x20;	ua8     = tmp2+ 0x20;
		t15i= sc_ptr + 0x21;		s1p16i = tmp + 0x21;	ua9     = tmp2+ 0x21;
		t16r= sc_ptr + 0x22;		s1p17r = tmp + 0x22;	ub0     = tmp2+ 0x22;
		t16i= sc_ptr + 0x23;		s1p17i = tmp + 0x23;	ub1     = tmp2+ 0x23;
		t17r= sc_ptr + 0x24;		s1p18r = tmp + 0x24;	ub2     = tmp2+ 0x24;
		t17i= sc_ptr + 0x25;		s1p18i = tmp + 0x25;	ub3     = tmp2+ 0x25;
		t18r= sc_ptr + 0x26;		s1p19r = tmp + 0x26;	ub4     = tmp2+ 0x26;
		t18i= sc_ptr + 0x27;		s1p19i = tmp + 0x27;	ub5     = tmp2+ 0x27;
		t19r= sc_ptr + 0x28;		s1p20r = tmp + 0x28;	ub6     = tmp2+ 0x28;
		t19i= sc_ptr + 0x29;		s1p20i = tmp + 0x29;	ub7     = tmp2+ 0x29;
		t1ar= sc_ptr + 0x2a;		s1p21r = tmp + 0x2a;	ub8     = tmp2+ 0x2a;
		t1ai= sc_ptr + 0x2b;		s1p21i = tmp + 0x2b;	ub9     = tmp2+ 0x2b;
		t20r= sc_ptr + 0x2c;		s1p22r = tmp + 0x2c;	max_err = tmp2+ 0x2c;
		t20i= sc_ptr + 0x2d;		s1p22i = tmp + 0x2d;	sse2_rnd= tmp2+ 0x2d;
		t21r= sc_ptr + 0x2e;		s1p23r = tmp + 0x2e;	half_arr= tmp2+ 0x2e;	/* This table needs 20x16 bytes */
		t21i= sc_ptr + 0x2f;		s1p23i = tmp + 0x2f;
		t22r= sc_ptr + 0x30;		s1p24r = tmp + 0x30;
		t22i= sc_ptr + 0x31;		s1p24i = tmp + 0x31;
		t23r= sc_ptr + 0x32;		s1p25r = tmp + 0x32;
		t23i= sc_ptr + 0x33;		s1p25i = tmp + 0x33;
		t24r= sc_ptr + 0x34;		s1p26r = tmp + 0x34;
		t24i= sc_ptr + 0x35;		s1p26i = tmp + 0x35;
		t25r= sc_ptr + 0x36;		s1p27r = tmp + 0x36;
		t25i= sc_ptr + 0x37;		s1p27i = tmp + 0x37;
		t26r= sc_ptr + 0x38;		s1p28r = tmp + 0x38;
		t26i= sc_ptr + 0x39;		s1p28i = tmp + 0x39;
		t27r= sc_ptr + 0x3a;		s1p29r = tmp + 0x3a;
		t27i= sc_ptr + 0x3b;		s1p29i = tmp + 0x3b;
		t28r= sc_ptr + 0x3c;		s1p30r = tmp + 0x3c;
		t28i= sc_ptr + 0x3d;		s1p30i = tmp + 0x3d;
		t29r= sc_ptr + 0x3e;		s1p31r = tmp + 0x3e;
		t29i= sc_ptr + 0x3f;		s1p31i = tmp + 0x3f;
		t2ar= sc_ptr + 0x40;		s1p32r = tmp + 0x40;
		t2ai= sc_ptr + 0x41;		s1p32i = tmp + 0x41;
		t30r= sc_ptr + 0x42;		s1p33r = tmp + 0x42;
		t30i= sc_ptr + 0x43;		s1p33i = tmp + 0x43;
		t31r= sc_ptr + 0x44;		s1p34r = tmp + 0x44;
		t31i= sc_ptr + 0x45;		s1p34i = tmp + 0x45;
		t32r= sc_ptr + 0x46;		s1p35r = tmp + 0x46;
		t32i= sc_ptr + 0x47;		s1p35i = tmp + 0x47;
		t33r= sc_ptr + 0x48;		s1p36r = tmp + 0x48;
		t33i= sc_ptr + 0x49;		s1p36i = tmp + 0x49;
		t34r= sc_ptr + 0x4a;		s1p37r = tmp + 0x4a;
		t34i= sc_ptr + 0x4b;		s1p37i = tmp + 0x4b;
		t35r= sc_ptr + 0x4c;		s1p38r = tmp + 0x4c;
		t35i= sc_ptr + 0x4d;		s1p38i = tmp + 0x4d;
		t36r= sc_ptr + 0x4e;		s1p39r = tmp + 0x4e;
		t36i= sc_ptr + 0x4f;		s1p39i = tmp + 0x4f;
		t37r= sc_ptr + 0x50;		s1p40r = tmp + 0x50;
		t37i= sc_ptr + 0x51;		s1p40i = tmp + 0x51;
		t38r= sc_ptr + 0x52;		s1p41r = tmp + 0x52;
		t38i= sc_ptr + 0x53;		s1p41i = tmp + 0x53;
		t39r= sc_ptr + 0x54;		s1p42r = tmp + 0x54;
		t39i= sc_ptr + 0x55;		s1p42i = tmp + 0x55;
		t3ar= sc_ptr + 0x56;		s1p43r = tmp + 0x56;
		t3ai= sc_ptr + 0x57;		s1p43i = tmp + 0x57;

		/* These remain fixed: */
		two ->re =2.0;		two ->im =2.0;
		five->re =5.0;		five->im =5.0;
		ua0 ->re = a0;		ua0 ->im = a0;
		ua1 ->re = a1;		ua1 ->im = a1;
		ua2 ->re = a2;		ua2 ->im = a2;
		ua3 ->re = a3;		ua3 ->im = a3;
		ua4 ->re = a4;		ua4 ->im = a4;
		ua5 ->re = a5;		ua5 ->im = a5;
		ua6 ->re = a6;		ua6 ->im = a6;
		ua7 ->re = a7;		ua7 ->im = a7;
		ua8 ->re = a8;		ua8 ->im = a8;
		ua9 ->re = a9;		ua9 ->im = a9;
		ub0 ->re = b0;		ub0 ->im = b0;
		ub1 ->re = b1;		ub1 ->im = b1;
		ub2 ->re = b2;		ub2 ->im = b2;
		ub3 ->re = b3;		ub3 ->im = b3;
		ub4 ->re = b4;		ub4 ->im = b4;
		ub5 ->re = b5;		ub5 ->im = b5;
		ub6 ->re = b6;		ub6 ->im = b6;
		ub7 ->re = b7;		ub7 ->im = b7;
		ub8 ->re = b8;		ub8 ->im = b8;
		ub9 ->re =-b9;		ub9 ->im =-b9;	/* Flip sign to simplify code re-use in radix-11 SSE2 macro */

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		sse2_rnd->re = 3.0*0x4000000*0x2000000;
		sse2_rnd->im = 3.0*0x4000000*0x2000000;

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

	#endif

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

		if(_cy_r00)	/* If it's a new exponent of a range test, need to deallocate these. */
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

			free((void *)_cy_r00); _cy_r00 = 0x0;
			free((void *)_cy_r01); _cy_r01 = 0x0;
			free((void *)_cy_r02); _cy_r02 = 0x0;
			free((void *)_cy_r03); _cy_r03 = 0x0;
			free((void *)_cy_r04); _cy_r04 = 0x0;
			free((void *)_cy_r05); _cy_r05 = 0x0;
			free((void *)_cy_r06); _cy_r06 = 0x0;
			free((void *)_cy_r07); _cy_r07 = 0x0;
			free((void *)_cy_r08); _cy_r08 = 0x0;
			free((void *)_cy_r09); _cy_r09 = 0x0;
			free((void *)_cy_r10); _cy_r10 = 0x0;
			free((void *)_cy_r11); _cy_r11 = 0x0;
			free((void *)_cy_r12); _cy_r12 = 0x0;
			free((void *)_cy_r13); _cy_r13 = 0x0;
			free((void *)_cy_r14); _cy_r14 = 0x0;
			free((void *)_cy_r15); _cy_r15 = 0x0;
			free((void *)_cy_r16); _cy_r16 = 0x0;
			free((void *)_cy_r17); _cy_r17 = 0x0;
			free((void *)_cy_r18); _cy_r18 = 0x0;
			free((void *)_cy_r19); _cy_r19 = 0x0;
			free((void *)_cy_r20); _cy_r20 = 0x0;
			free((void *)_cy_r21); _cy_r21 = 0x0;
			free((void *)_cy_r22); _cy_r22 = 0x0;
			free((void *)_cy_r23); _cy_r23 = 0x0;
			free((void *)_cy_r24); _cy_r24 = 0x0;
			free((void *)_cy_r25); _cy_r25 = 0x0;
			free((void *)_cy_r26); _cy_r26 = 0x0;
			free((void *)_cy_r27); _cy_r27 = 0x0;
			free((void *)_cy_r28); _cy_r28 = 0x0;
			free((void *)_cy_r29); _cy_r29 = 0x0;
			free((void *)_cy_r30); _cy_r30 = 0x0;
			free((void *)_cy_r31); _cy_r31 = 0x0;
			free((void *)_cy_r32); _cy_r32 = 0x0;
			free((void *)_cy_r33); _cy_r33 = 0x0;
			free((void *)_cy_r34); _cy_r34 = 0x0;
			free((void *)_cy_r35); _cy_r35 = 0x0;
			free((void *)_cy_r36); _cy_r36 = 0x0;
			free((void *)_cy_r37); _cy_r37 = 0x0;
			free((void *)_cy_r38); _cy_r38 = 0x0;
			free((void *)_cy_r39); _cy_r39 = 0x0;
			free((void *)_cy_r40); _cy_r40 = 0x0;
			free((void *)_cy_r41); _cy_r41 = 0x0;
			free((void *)_cy_r42); _cy_r42 = 0x0;
			free((void *)_cy_r43); _cy_r43 = 0x0;

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		_i       	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_i== 0x0);
		_bjmodn00	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn00== 0x0);
		_bjmodn01	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn01== 0x0);
		_bjmodn02	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn02== 0x0);
		_bjmodn03	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn03== 0x0);
		_bjmodn04	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn04== 0x0);
		_bjmodn05	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn05== 0x0);
		_bjmodn06	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn06== 0x0);
		_bjmodn07	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn07== 0x0);
		_bjmodn08	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn08== 0x0);
		_bjmodn09	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn09== 0x0);
		_bjmodn10	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn10== 0x0);
		_bjmodn11	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn11== 0x0);
		_bjmodn12	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn12== 0x0);
		_bjmodn13	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn13== 0x0);
		_bjmodn14	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn14== 0x0);
		_bjmodn15	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn15== 0x0);
		_bjmodn16	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn16== 0x0);
		_bjmodn17	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn17== 0x0);
		_bjmodn18	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn18== 0x0);
		_bjmodn19	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn19== 0x0);
		_bjmodn20	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_bjmodn24	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_bjmodn28	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn28== 0x0);
		_bjmodn29	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn29== 0x0);
		_bjmodn30	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn30== 0x0);
		_bjmodn31	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn31== 0x0);
		_bjmodn32	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn32== 0x0);
		_bjmodn33	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn33== 0x0);
		_bjmodn34	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn34== 0x0);
		_bjmodn35	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn35== 0x0);
		_bjmodn36	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn36== 0x0);
		_bjmodn37	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn37== 0x0);
		_bjmodn38	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn38== 0x0);
		_bjmodn39	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn39== 0x0);
		_bjmodn40	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn40== 0x0);
		_bjmodn41	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn41== 0x0);
		_bjmodn42	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn42== 0x0);
		_bjmodn43	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn43== 0x0);
		_jstart  	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co3     == 0x0);

		_cy_r00	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r00== 0x0);
		_cy_r01	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r01== 0x0);
		_cy_r02	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r02== 0x0);
		_cy_r03	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r03== 0x0);
		_cy_r04	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r04== 0x0);
		_cy_r05	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r05== 0x0);
		_cy_r06	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r06== 0x0);
		_cy_r07	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r07== 0x0);
		_cy_r08	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r08== 0x0);
		_cy_r09	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r09== 0x0);
		_cy_r10	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r10== 0x0);
		_cy_r11	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r11== 0x0);
		_cy_r12	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r12== 0x0);
		_cy_r13	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r13== 0x0);
		_cy_r14	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r14== 0x0);
		_cy_r15	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r15== 0x0);
		_cy_r16	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r16== 0x0);
		_cy_r17	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r17== 0x0);
		_cy_r18	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r18== 0x0);
		_cy_r19	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r19== 0x0);
		_cy_r20	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r20== 0x0);
		_cy_r21	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r21== 0x0);
		_cy_r22	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r22== 0x0);
		_cy_r23	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r23== 0x0);
		_cy_r24	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r24== 0x0);
		_cy_r25	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r25== 0x0);
		_cy_r26	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r26== 0x0);
		_cy_r27	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r27== 0x0);
		_cy_r28	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r28== 0x0);
		_cy_r29	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r29== 0x0);
		_cy_r30	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r30== 0x0);
		_cy_r31	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r31== 0x0);
		_cy_r32	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r32== 0x0);
		_cy_r33	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r33== 0x0);
		_cy_r34	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r34== 0x0);
		_cy_r35	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r35== 0x0);
		_cy_r36	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r36== 0x0);
		_cy_r37	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r37== 0x0);
		_cy_r38	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r38== 0x0);
		_cy_r39	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r39== 0x0);
		_cy_r40	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r40== 0x0);
		_cy_r41	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r41== 0x0);
		_cy_r42	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r42== 0x0);
		_cy_r43	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r43== 0x0);

		_maxerr	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_maxerr== 0x0);

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
		jhi = n44/CY_THREADS;

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
		for(j=0; j < n44; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	}	/* endif(first_entry) */

/*...The radix-44 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r00[ithread] = 0;
		_cy_r01[ithread] = 0;
		_cy_r02[ithread] = 0;
		_cy_r03[ithread] = 0;
		_cy_r04[ithread] = 0;
		_cy_r05[ithread] = 0;
		_cy_r06[ithread] = 0;
		_cy_r07[ithread] = 0;
		_cy_r08[ithread] = 0;
		_cy_r09[ithread] = 0;
		_cy_r10[ithread] = 0;
		_cy_r11[ithread] = 0;
		_cy_r12[ithread] = 0;
		_cy_r13[ithread] = 0;
		_cy_r14[ithread] = 0;
		_cy_r15[ithread] = 0;
		_cy_r16[ithread] = 0;
		_cy_r17[ithread] = 0;
		_cy_r18[ithread] = 0;
		_cy_r19[ithread] = 0;
		_cy_r20[ithread] = 0;
		_cy_r21[ithread] = 0;
		_cy_r22[ithread] = 0;
		_cy_r23[ithread] = 0;
		_cy_r24[ithread] = 0;
		_cy_r25[ithread] = 0;
		_cy_r26[ithread] = 0;
		_cy_r27[ithread] = 0;
		_cy_r28[ithread] = 0;
		_cy_r29[ithread] = 0;
		_cy_r30[ithread] = 0;
		_cy_r31[ithread] = 0;
		_cy_r32[ithread] = 0;
		_cy_r33[ithread] = 0;
		_cy_r34[ithread] = 0;
		_cy_r35[ithread] = 0;
		_cy_r36[ithread] = 0;
		_cy_r37[ithread] = 0;
		_cy_r38[ithread] = 0;
		_cy_r39[ithread] = 0;
		_cy_r40[ithread] = 0;
		_cy_r41[ithread] = 0;
		_cy_r42[ithread] = 0;
		_cy_r43[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r00[      0] = -2;
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

	/*
	Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
	then simply overwrite it with 1 prior to starting the k-loop.
	*/
	khi = n_div_nwt/CY_THREADS;

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_jstart[ithread] = ithread*n44/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*44);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+44 - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-44;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
	}

	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn00[ithread] = _bjmodnini[ithread];
		_bjmodn01[ithread] = _bjmodn00[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn01[ithread] = _bjmodn01[ithread] + ( (-(int)((uint32)_bjmodn01[ithread] >> 31)) & n);
		_bjmodn02[ithread] = _bjmodn01[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn02[ithread] = _bjmodn02[ithread] + ( (-(int)((uint32)_bjmodn02[ithread] >> 31)) & n);
		_bjmodn03[ithread] = _bjmodn02[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn03[ithread] = _bjmodn03[ithread] + ( (-(int)((uint32)_bjmodn03[ithread] >> 31)) & n);
		_bjmodn04[ithread] = _bjmodn03[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn04[ithread] = _bjmodn04[ithread] + ( (-(int)((uint32)_bjmodn04[ithread] >> 31)) & n);
		_bjmodn05[ithread] = _bjmodn04[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn05[ithread] = _bjmodn05[ithread] + ( (-(int)((uint32)_bjmodn05[ithread] >> 31)) & n);
		_bjmodn06[ithread] = _bjmodn05[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn06[ithread] = _bjmodn06[ithread] + ( (-(int)((uint32)_bjmodn06[ithread] >> 31)) & n);
		_bjmodn07[ithread] = _bjmodn06[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn07[ithread] = _bjmodn07[ithread] + ( (-(int)((uint32)_bjmodn07[ithread] >> 31)) & n);
		_bjmodn08[ithread] = _bjmodn07[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn08[ithread] = _bjmodn08[ithread] + ( (-(int)((uint32)_bjmodn08[ithread] >> 31)) & n);
		_bjmodn09[ithread] = _bjmodn08[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn09[ithread] = _bjmodn09[ithread] + ( (-(int)((uint32)_bjmodn09[ithread] >> 31)) & n);
		_bjmodn10[ithread] = _bjmodn09[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn10[ithread] = _bjmodn10[ithread] + ( (-(int)((uint32)_bjmodn10[ithread] >> 31)) & n);
		_bjmodn11[ithread] = _bjmodn10[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn11[ithread] = _bjmodn11[ithread] + ( (-(int)((uint32)_bjmodn11[ithread] >> 31)) & n);
		_bjmodn12[ithread] = _bjmodn11[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn12[ithread] = _bjmodn12[ithread] + ( (-(int)((uint32)_bjmodn12[ithread] >> 31)) & n);
		_bjmodn13[ithread] = _bjmodn12[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn13[ithread] = _bjmodn13[ithread] + ( (-(int)((uint32)_bjmodn13[ithread] >> 31)) & n);
		_bjmodn14[ithread] = _bjmodn13[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn14[ithread] = _bjmodn14[ithread] + ( (-(int)((uint32)_bjmodn14[ithread] >> 31)) & n);
		_bjmodn15[ithread] = _bjmodn14[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn15[ithread] = _bjmodn15[ithread] + ( (-(int)((uint32)_bjmodn15[ithread] >> 31)) & n);
		_bjmodn16[ithread] = _bjmodn15[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn16[ithread] = _bjmodn16[ithread] + ( (-(int)((uint32)_bjmodn16[ithread] >> 31)) & n);
		_bjmodn17[ithread] = _bjmodn16[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn17[ithread] = _bjmodn17[ithread] + ( (-(int)((uint32)_bjmodn17[ithread] >> 31)) & n);
		_bjmodn18[ithread] = _bjmodn17[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn18[ithread] = _bjmodn18[ithread] + ( (-(int)((uint32)_bjmodn18[ithread] >> 31)) & n);
		_bjmodn19[ithread] = _bjmodn18[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn19[ithread] = _bjmodn19[ithread] + ( (-(int)((uint32)_bjmodn19[ithread] >> 31)) & n);
		_bjmodn20[ithread] = _bjmodn19[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn20[ithread] = _bjmodn20[ithread] + ( (-(int)((uint32)_bjmodn20[ithread] >> 31)) & n);
		_bjmodn21[ithread] = _bjmodn20[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn21[ithread] = _bjmodn21[ithread] + ( (-(int)((uint32)_bjmodn21[ithread] >> 31)) & n);
		_bjmodn22[ithread] = _bjmodn21[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn22[ithread] = _bjmodn22[ithread] + ( (-(int)((uint32)_bjmodn22[ithread] >> 31)) & n);
		_bjmodn23[ithread] = _bjmodn22[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn23[ithread] = _bjmodn23[ithread] + ( (-(int)((uint32)_bjmodn23[ithread] >> 31)) & n);
		_bjmodn24[ithread] = _bjmodn23[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn24[ithread] = _bjmodn24[ithread] + ( (-(int)((uint32)_bjmodn24[ithread] >> 31)) & n);
		_bjmodn25[ithread] = _bjmodn24[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn25[ithread] = _bjmodn25[ithread] + ( (-(int)((uint32)_bjmodn25[ithread] >> 31)) & n);
		_bjmodn26[ithread] = _bjmodn25[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn26[ithread] = _bjmodn26[ithread] + ( (-(int)((uint32)_bjmodn26[ithread] >> 31)) & n);
		_bjmodn27[ithread] = _bjmodn26[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn27[ithread] = _bjmodn27[ithread] + ( (-(int)((uint32)_bjmodn27[ithread] >> 31)) & n);
		_bjmodn28[ithread] = _bjmodn27[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn28[ithread] = _bjmodn28[ithread] + ( (-(int)((uint32)_bjmodn28[ithread] >> 31)) & n);
		_bjmodn29[ithread] = _bjmodn28[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn29[ithread] = _bjmodn29[ithread] + ( (-(int)((uint32)_bjmodn29[ithread] >> 31)) & n);
		_bjmodn30[ithread] = _bjmodn29[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn30[ithread] = _bjmodn30[ithread] + ( (-(int)((uint32)_bjmodn30[ithread] >> 31)) & n);
		_bjmodn31[ithread] = _bjmodn30[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn31[ithread] = _bjmodn31[ithread] + ( (-(int)((uint32)_bjmodn31[ithread] >> 31)) & n);
		_bjmodn32[ithread] = _bjmodn31[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn32[ithread] = _bjmodn32[ithread] + ( (-(int)((uint32)_bjmodn32[ithread] >> 31)) & n);
		_bjmodn33[ithread] = _bjmodn32[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn33[ithread] = _bjmodn33[ithread] + ( (-(int)((uint32)_bjmodn33[ithread] >> 31)) & n);
		_bjmodn34[ithread] = _bjmodn33[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn34[ithread] = _bjmodn34[ithread] + ( (-(int)((uint32)_bjmodn34[ithread] >> 31)) & n);
		_bjmodn35[ithread] = _bjmodn34[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn35[ithread] = _bjmodn35[ithread] + ( (-(int)((uint32)_bjmodn35[ithread] >> 31)) & n);
		_bjmodn36[ithread] = _bjmodn35[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn36[ithread] = _bjmodn36[ithread] + ( (-(int)((uint32)_bjmodn36[ithread] >> 31)) & n);
		_bjmodn37[ithread] = _bjmodn36[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn37[ithread] = _bjmodn37[ithread] + ( (-(int)((uint32)_bjmodn37[ithread] >> 31)) & n);
		_bjmodn38[ithread] = _bjmodn37[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn38[ithread] = _bjmodn38[ithread] + ( (-(int)((uint32)_bjmodn38[ithread] >> 31)) & n);
		_bjmodn39[ithread] = _bjmodn38[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn39[ithread] = _bjmodn39[ithread] + ( (-(int)((uint32)_bjmodn39[ithread] >> 31)) & n);
		_bjmodn40[ithread] = _bjmodn39[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn40[ithread] = _bjmodn40[ithread] + ( (-(int)((uint32)_bjmodn40[ithread] >> 31)) & n);
		_bjmodn41[ithread] = _bjmodn40[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn41[ithread] = _bjmodn41[ithread] + ( (-(int)((uint32)_bjmodn41[ithread] >> 31)) & n);
		_bjmodn42[ithread] = _bjmodn41[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn42[ithread] = _bjmodn42[ithread] + ( (-(int)((uint32)_bjmodn42[ithread] >> 31)) & n);
		_bjmodn43[ithread] = _bjmodn42[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn43[ithread] = _bjmodn43[ithread] + ( (-(int)((uint32)_bjmodn43[ithread] >> 31)) & n);
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
#ifdef MULTITHREAD
	#error radix44_ditN_cy_dif1: OpenMP private list needs updating!
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,k1,k2,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,rt,it,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,cy_i28,cy_i29,cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35) default(shared) schedule(static)
#endif

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
		*bjmodn00 = _bjmodn00[ithread];			cy_r00->re = _cy_r00[ithread];
		*bjmodn01 = _bjmodn01[ithread];			cy_r00->im = _cy_r01[ithread];
		*bjmodn02 = _bjmodn02[ithread];			cy_r02->re = _cy_r02[ithread];
		*bjmodn03 = _bjmodn03[ithread];			cy_r02->im = _cy_r03[ithread];
		*bjmodn04 = _bjmodn04[ithread];			cy_r04->re = _cy_r04[ithread];
		*bjmodn05 = _bjmodn05[ithread];			cy_r04->im = _cy_r05[ithread];
		*bjmodn06 = _bjmodn06[ithread];			cy_r06->re = _cy_r06[ithread];
		*bjmodn07 = _bjmodn07[ithread];			cy_r06->im = _cy_r07[ithread];
		*bjmodn08 = _bjmodn08[ithread];			cy_r08->re = _cy_r08[ithread];
		*bjmodn09 = _bjmodn09[ithread];			cy_r08->im = _cy_r09[ithread];
		*bjmodn10 = _bjmodn10[ithread];			cy_r10->re = _cy_r10[ithread];
		*bjmodn11 = _bjmodn11[ithread];			cy_r10->im = _cy_r11[ithread];
		*bjmodn12 = _bjmodn12[ithread];			cy_r12->re = _cy_r12[ithread];
		*bjmodn13 = _bjmodn13[ithread];			cy_r12->im = _cy_r13[ithread];
		*bjmodn14 = _bjmodn14[ithread];			cy_r14->re = _cy_r14[ithread];
		*bjmodn15 = _bjmodn15[ithread];			cy_r14->im = _cy_r15[ithread];
		*bjmodn16 = _bjmodn16[ithread];			cy_r16->re = _cy_r16[ithread];
		*bjmodn17 = _bjmodn17[ithread];			cy_r16->im = _cy_r17[ithread];
		*bjmodn18 = _bjmodn18[ithread];			cy_r18->re = _cy_r18[ithread];
		*bjmodn19 = _bjmodn19[ithread];			cy_r18->im = _cy_r19[ithread];
		*bjmodn20 = _bjmodn20[ithread];			cy_r20->re = _cy_r20[ithread];
		*bjmodn21 = _bjmodn21[ithread];			cy_r20->im = _cy_r21[ithread];
		*bjmodn22 = _bjmodn22[ithread];			cy_r22->re = _cy_r22[ithread];
		*bjmodn23 = _bjmodn23[ithread];			cy_r22->im = _cy_r23[ithread];
		*bjmodn24 = _bjmodn24[ithread];			cy_r24->re = _cy_r24[ithread];
		*bjmodn25 = _bjmodn25[ithread];			cy_r24->im = _cy_r25[ithread];
		*bjmodn26 = _bjmodn26[ithread];			cy_r26->re = _cy_r26[ithread];
		*bjmodn27 = _bjmodn27[ithread];			cy_r26->im = _cy_r27[ithread];
		*bjmodn28 = _bjmodn28[ithread];			cy_r28->re = _cy_r28[ithread];
		*bjmodn29 = _bjmodn29[ithread];			cy_r28->im = _cy_r29[ithread];
		*bjmodn30 = _bjmodn30[ithread];			cy_r30->re = _cy_r30[ithread];
		*bjmodn31 = _bjmodn31[ithread];			cy_r30->im = _cy_r31[ithread];
		*bjmodn32 = _bjmodn32[ithread];			cy_r32->re = _cy_r32[ithread];
		*bjmodn33 = _bjmodn33[ithread];			cy_r32->im = _cy_r33[ithread];
		*bjmodn34 = _bjmodn34[ithread];			cy_r34->re = _cy_r34[ithread];
		*bjmodn35 = _bjmodn35[ithread];			cy_r34->im = _cy_r35[ithread];
		*bjmodn36 = _bjmodn36[ithread];			cy_r36->re = _cy_r36[ithread];
		*bjmodn37 = _bjmodn37[ithread];			cy_r36->im = _cy_r37[ithread];
		*bjmodn38 = _bjmodn38[ithread];			cy_r38->re = _cy_r38[ithread];
		*bjmodn39 = _bjmodn39[ithread];			cy_r38->im = _cy_r39[ithread];
		*bjmodn40 = _bjmodn40[ithread];			cy_r40->re = _cy_r40[ithread];
		*bjmodn41 = _bjmodn41[ithread];			cy_r40->im = _cy_r41[ithread];
		*bjmodn42 = _bjmodn42[ithread];			cy_r42->re = _cy_r42[ithread];
		*bjmodn43 = _bjmodn43[ithread];			cy_r42->im = _cy_r43[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];			cy_r00 = _cy_r00[ithread];
		bjmodn01 = _bjmodn01[ithread];			cy_r01 = _cy_r01[ithread];
		bjmodn02 = _bjmodn02[ithread];			cy_r02 = _cy_r02[ithread];
		bjmodn03 = _bjmodn03[ithread];			cy_r03 = _cy_r03[ithread];
		bjmodn04 = _bjmodn04[ithread];			cy_r04 = _cy_r04[ithread];
		bjmodn05 = _bjmodn05[ithread];			cy_r05 = _cy_r05[ithread];
		bjmodn06 = _bjmodn06[ithread];			cy_r06 = _cy_r06[ithread];
		bjmodn07 = _bjmodn07[ithread];			cy_r07 = _cy_r07[ithread];
		bjmodn08 = _bjmodn08[ithread];			cy_r08 = _cy_r08[ithread];
		bjmodn09 = _bjmodn09[ithread];			cy_r09 = _cy_r09[ithread];
		bjmodn10 = _bjmodn10[ithread];			cy_r10 = _cy_r10[ithread];
		bjmodn11 = _bjmodn11[ithread];			cy_r11 = _cy_r11[ithread];
		bjmodn12 = _bjmodn12[ithread];			cy_r12 = _cy_r12[ithread];
		bjmodn13 = _bjmodn13[ithread];			cy_r13 = _cy_r13[ithread];
		bjmodn14 = _bjmodn14[ithread];			cy_r14 = _cy_r14[ithread];
		bjmodn15 = _bjmodn15[ithread];			cy_r15 = _cy_r15[ithread];
		bjmodn16 = _bjmodn16[ithread];			cy_r16 = _cy_r16[ithread];
		bjmodn17 = _bjmodn17[ithread];			cy_r17 = _cy_r17[ithread];
		bjmodn18 = _bjmodn18[ithread];			cy_r18 = _cy_r18[ithread];
		bjmodn19 = _bjmodn19[ithread];			cy_r19 = _cy_r19[ithread];
		bjmodn20 = _bjmodn20[ithread];			cy_r20 = _cy_r20[ithread];
		bjmodn21 = _bjmodn21[ithread];			cy_r21 = _cy_r21[ithread];
		bjmodn22 = _bjmodn22[ithread];			cy_r22 = _cy_r22[ithread];
		bjmodn23 = _bjmodn23[ithread];			cy_r23 = _cy_r23[ithread];
		bjmodn24 = _bjmodn24[ithread];			cy_r24 = _cy_r24[ithread];
		bjmodn25 = _bjmodn25[ithread];			cy_r25 = _cy_r25[ithread];
		bjmodn26 = _bjmodn26[ithread];			cy_r26 = _cy_r26[ithread];
		bjmodn27 = _bjmodn27[ithread];			cy_r27 = _cy_r27[ithread];
		bjmodn28 = _bjmodn28[ithread];			cy_r28 = _cy_r28[ithread];
		bjmodn29 = _bjmodn29[ithread];			cy_r29 = _cy_r29[ithread];
		bjmodn30 = _bjmodn30[ithread];			cy_r30 = _cy_r30[ithread];
		bjmodn31 = _bjmodn31[ithread];			cy_r31 = _cy_r31[ithread];
		bjmodn32 = _bjmodn32[ithread];			cy_r32 = _cy_r32[ithread];
		bjmodn33 = _bjmodn33[ithread];			cy_r33 = _cy_r33[ithread];
		bjmodn34 = _bjmodn34[ithread];			cy_r34 = _cy_r34[ithread];
		bjmodn35 = _bjmodn35[ithread];			cy_r35 = _cy_r35[ithread];
		bjmodn36 = _bjmodn36[ithread];			cy_r36 = _cy_r36[ithread];
		bjmodn37 = _bjmodn37[ithread];			cy_r37 = _cy_r37[ithread];
		bjmodn38 = _bjmodn38[ithread];			cy_r38 = _cy_r38[ithread];
		bjmodn39 = _bjmodn39[ithread];			cy_r39 = _cy_r39[ithread];
		bjmodn40 = _bjmodn40[ithread];			cy_r40 = _cy_r40[ithread];
		bjmodn41 = _bjmodn41[ithread];			cy_r41 = _cy_r41[ithread];
		bjmodn42 = _bjmodn42[ithread];			cy_r42 = _cy_r42[ithread];
		bjmodn43 = _bjmodn43[ithread];			cy_r43 = _cy_r43[ithread];
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

		#ifdef DEBUG_SSE2
			rng_isaac_init(TRUE);
			jt = j1;		jp = j2;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[00] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[01] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[02] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[03] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[04] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[05] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[06] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[07] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[08] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[09] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[10] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[11] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[12] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[13] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[14] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[15] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[16] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[17] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[18] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[29] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[20] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[21] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[22] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[23] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[24] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[25] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[26] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[27] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[28] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[29] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[30] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[31] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[32] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[33] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[34] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[35] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[36] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[37] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[38] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[39] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[40] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[41] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[42] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();		fprintf(stderr, "radix44: A_in[43] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "\n");
		#endif

		#ifdef USE_SSE2

		  #if GCC_ASM_FULL_INLINE	// GCC or SUNC implied

			add0 = &a[j1    ];
			SSE2_RADIX44_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,t00r, ua0, s1p00r);

		  #elif !defined(COMPILER_TYPE_MSVC) && !GCC_ASM_FULL_INLINE

			/* Outputs in SSE2 modes are temps 2*11*16 = 22*16 = 0x160 bytes apart: */
			add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t00r, 0x160)
			add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t08r, 0x160)
			add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t05r, 0x160)
			add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t02r, 0x160)
			add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t0ar, 0x160)
			add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t07r, 0x160)
			add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t04r, 0x160)
			add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t01r, 0x160)
			add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t09r, 0x160)
			add2 = &a[j1+p36];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t06r, 0x160)
			add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, t03r, 0x160)

			/* Radix-11 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 12 (24*16 bytes = 0x180) between successive outputs: */
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

		  #if 0//DEBUG_SSE2
			fprintf(stderr, "radix44_wrapper: tmp[00] = %24.5f, %24.5f\n",t00r->re,t00i->re);
			fprintf(stderr, "radix44_wrapper: tmp[10] = %24.5f, %24.5f\n",t10r->re,t10i->re);
			fprintf(stderr, "radix44_wrapper: tmp[20] = %24.5f, %24.5f\n",t20r->re,t20i->re);
			fprintf(stderr, "radix44_wrapper: tmp[30] = %24.5f, %24.5f\n",t30r->re,t30i->re);
			fprintf(stderr, "radix44_wrapper: tmp[01] = %24.5f, %24.5f\n",t01r->re,t01i->re);
			fprintf(stderr, "radix44_wrapper: tmp[11] = %24.5f, %24.5f\n",t11r->re,t11i->re);
			fprintf(stderr, "radix44_wrapper: tmp[21] = %24.5f, %24.5f\n",t21r->re,t21i->re);
			fprintf(stderr, "radix44_wrapper: tmp[31] = %24.5f, %24.5f\n",t31r->re,t31i->re);
			fprintf(stderr, "radix44_wrapper: tmp[02] = %24.5f, %24.5f\n",t02r->re,t02i->re);
			fprintf(stderr, "radix44_wrapper: tmp[12] = %24.5f, %24.5f\n",t12r->re,t12i->re);
			fprintf(stderr, "radix44_wrapper: tmp[22] = %24.5f, %24.5f\n",t22r->re,t22i->re);
			fprintf(stderr, "radix44_wrapper: tmp[32] = %24.5f, %24.5f\n",t32r->re,t32i->re);
			fprintf(stderr, "radix44_wrapper: tmp[03] = %24.5f, %24.5f\n",t03r->re,t03i->re);
			fprintf(stderr, "radix44_wrapper: tmp[13] = %24.5f, %24.5f\n",t13r->re,t13i->re);
			fprintf(stderr, "radix44_wrapper: tmp[23] = %24.5f, %24.5f\n",t23r->re,t23i->re);
			fprintf(stderr, "radix44_wrapper: tmp[33] = %24.5f, %24.5f\n",t33r->re,t33i->re);
			fprintf(stderr, "radix44_wrapper: tmp[04] = %24.5f, %24.5f\n",t04r->re,t04i->re);
			fprintf(stderr, "radix44_wrapper: tmp[14] = %24.5f, %24.5f\n",t14r->re,t14i->re);
			fprintf(stderr, "radix44_wrapper: tmp[24] = %24.5f, %24.5f\n",t24r->re,t24i->re);
			fprintf(stderr, "radix44_wrapper: tmp[34] = %24.5f, %24.5f\n",t34r->re,t34i->re);
			fprintf(stderr, "radix44_wrapper: tmp[05] = %24.5f, %24.5f\n",t05r->re,t05i->re);
			fprintf(stderr, "radix44_wrapper: tmp[15] = %24.5f, %24.5f\n",t15r->re,t15i->re);
			fprintf(stderr, "radix44_wrapper: tmp[25] = %24.5f, %24.5f\n",t25r->re,t25i->re);
			fprintf(stderr, "radix44_wrapper: tmp[35] = %24.5f, %24.5f\n",t35r->re,t35i->re);
			fprintf(stderr, "radix44_wrapper: tmp[06] = %24.5f, %24.5f\n",t06r->re,t06i->re);
			fprintf(stderr, "radix44_wrapper: tmp[16] = %24.5f, %24.5f\n",t16r->re,t16i->re);
			fprintf(stderr, "radix44_wrapper: tmp[26] = %24.5f, %24.5f\n",t26r->re,t26i->re);
			fprintf(stderr, "radix44_wrapper: tmp[36] = %24.5f, %24.5f\n",t36r->re,t36i->re);
			fprintf(stderr, "radix44_wrapper: tmp[07] = %24.5f, %24.5f\n",t07r->re,t07i->re);
			fprintf(stderr, "radix44_wrapper: tmp[17] = %24.5f, %24.5f\n",t17r->re,t17i->re);
			fprintf(stderr, "radix44_wrapper: tmp[27] = %24.5f, %24.5f\n",t27r->re,t27i->re);
			fprintf(stderr, "radix44_wrapper: tmp[37] = %24.5f, %24.5f\n",t37r->re,t37i->re);
			fprintf(stderr, "radix44_wrapper: tmp[08] = %24.5f, %24.5f\n",t08r->re,t08i->re);
			fprintf(stderr, "radix44_wrapper: tmp[18] = %24.5f, %24.5f\n",t18r->re,t18i->re);
			fprintf(stderr, "radix44_wrapper: tmp[28] = %24.5f, %24.5f\n",t28r->re,t28i->re);
			fprintf(stderr, "radix44_wrapper: tmp[38] = %24.5f, %24.5f\n",t38r->re,t38i->re);
			fprintf(stderr, "radix44_wrapper: tmp[09] = %24.5f, %24.5f\n",t09r->re,t09i->re);
			fprintf(stderr, "radix44_wrapper: tmp[19] = %24.5f, %24.5f\n",t19r->re,t19i->re);
			fprintf(stderr, "radix44_wrapper: tmp[29] = %24.5f, %24.5f\n",t29r->re,t29i->re);
			fprintf(stderr, "radix44_wrapper: tmp[39] = %24.5f, %24.5f\n",t39r->re,t39i->re);
			fprintf(stderr, "radix44_wrapper: tmp[0a] = %24.5f, %24.5f\n",t0ar->re,t0ai->re);
			fprintf(stderr, "radix44_wrapper: tmp[1a] = %24.5f, %24.5f\n",t1ar->re,t1ai->re);
			fprintf(stderr, "radix44_wrapper: tmp[2a] = %24.5f, %24.5f\n",t2ar->re,t2ai->re);
			fprintf(stderr, "radix44_wrapper: tmp[3a] = %24.5f, %24.5f\n",t3ar->re,t3ai->re);
			exit(0);
		  #endif

		  #ifdef DEBUG_SSE2
			fprintf(stderr, "a1p00 = %20.5f, %20.5f\n",s1p00r->re,s1p00i->re);
			fprintf(stderr, "a1p12 = %20.5f, %20.5f\n",s1p12r->re,s1p12i->re);
			fprintf(stderr, "a1p24 = %20.5f, %20.5f\n",s1p24r->re,s1p24i->re);
			fprintf(stderr, "a1p36 = %20.5f, %20.5f\n",s1p36r->re,s1p36i->re);
			fprintf(stderr, "a1p04 = %20.5f, %20.5f\n",s1p04r->re,s1p04i->re);
			fprintf(stderr, "a1p16 = %20.5f, %20.5f\n",s1p16r->re,s1p16i->re);
			fprintf(stderr, "a1p28 = %20.5f, %20.5f\n",s1p28r->re,s1p28i->re);
			fprintf(stderr, "a1p40 = %20.5f, %20.5f\n",s1p40r->re,s1p40i->re);
			fprintf(stderr, "a1p08 = %20.5f, %20.5f\n",s1p08r->re,s1p08i->re);
			fprintf(stderr, "a1p20 = %20.5f, %20.5f\n",s1p20r->re,s1p20i->re);
			fprintf(stderr, "a1p32 = %20.5f, %20.5f\n",s1p32r->re,s1p32i->re);	fprintf(stderr, "\n");
			fprintf(stderr, "a1p11 = %20.5f, %20.5f\n",s1p11r->re,s1p11i->re);
			fprintf(stderr, "a1p23 = %20.5f, %20.5f\n",s1p23r->re,s1p23i->re);
			fprintf(stderr, "a1p35 = %20.5f, %20.5f\n",s1p35r->re,s1p35i->re);
			fprintf(stderr, "a1p03 = %20.5f, %20.5f\n",s1p03r->re,s1p03i->re);
			fprintf(stderr, "a1p15 = %20.5f, %20.5f\n",s1p15r->re,s1p15i->re);
			fprintf(stderr, "a1p27 = %20.5f, %20.5f\n",s1p27r->re,s1p27i->re);
			fprintf(stderr, "a1p39 = %20.5f, %20.5f\n",s1p39r->re,s1p39i->re);
			fprintf(stderr, "a1p07 = %20.5f, %20.5f\n",s1p07r->re,s1p07i->re);
			fprintf(stderr, "a1p19 = %20.5f, %20.5f\n",s1p19r->re,s1p19i->re);
			fprintf(stderr, "a1p31 = %20.5f, %20.5f\n",s1p31r->re,s1p31i->re);
			fprintf(stderr, "a1p43 = %20.5f, %20.5f\n",s1p43r->re,s1p43i->re);	fprintf(stderr, "\n");
			fprintf(stderr, "a1p22 = %20.5f, %20.5f\n",s1p22r->re,s1p22i->re);
			fprintf(stderr, "a1p34 = %20.5f, %20.5f\n",s1p34r->re,s1p34i->re);
			fprintf(stderr, "a1p02 = %20.5f, %20.5f\n",s1p02r->re,s1p02i->re);
			fprintf(stderr, "a1p14 = %20.5f, %20.5f\n",s1p14r->re,s1p14i->re);
			fprintf(stderr, "a1p26 = %20.5f, %20.5f\n",s1p26r->re,s1p26i->re);
			fprintf(stderr, "a1p38 = %20.5f, %20.5f\n",s1p38r->re,s1p38i->re);
			fprintf(stderr, "a1p06 = %20.5f, %20.5f\n",s1p06r->re,s1p06i->re);
			fprintf(stderr, "a1p18 = %20.5f, %20.5f\n",s1p18r->re,s1p18i->re);
			fprintf(stderr, "a1p30 = %20.5f, %20.5f\n",s1p30r->re,s1p30i->re);
			fprintf(stderr, "a1p42 = %20.5f, %20.5f\n",s1p42r->re,s1p42i->re);
			fprintf(stderr, "a1p10 = %20.5f, %20.5f\n",s1p10r->re,s1p10i->re);	fprintf(stderr, "\n");
			fprintf(stderr, "a1p33 = %20.5f, %20.5f\n",s1p33r->re,s1p33i->re);
			fprintf(stderr, "a1p01 = %20.5f, %20.5f\n",s1p01r->re,s1p01i->re);
			fprintf(stderr, "a1p13 = %20.5f, %20.5f\n",s1p13r->re,s1p13i->re);
			fprintf(stderr, "a1p25 = %20.5f, %20.5f\n",s1p25r->re,s1p25i->re);
			fprintf(stderr, "a1p37 = %20.5f, %20.5f\n",s1p37r->re,s1p37i->re);
			fprintf(stderr, "a1p05 = %20.5f, %20.5f\n",s1p05r->re,s1p05i->re);
			fprintf(stderr, "a1p17 = %20.5f, %20.5f\n",s1p17r->re,s1p17i->re);
			fprintf(stderr, "a1p29 = %20.5f, %20.5f\n",s1p29r->re,s1p29i->re);
			fprintf(stderr, "a1p41 = %20.5f, %20.5f\n",s1p41r->re,s1p41i->re);
			fprintf(stderr, "a1p09 = %20.5f, %20.5f\n",s1p09r->re,s1p09i->re);
			fprintf(stderr, "a1p21 = %20.5f, %20.5f\n",s1p21r->re,s1p21i->re);	fprintf(stderr, "\n");
			exit(0);
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
#ifdef DEBUG_SSE2
	fprintf(stderr, "a1p00 = %20.5f, %20.5f\n",a1p00r,a1p00i);
	fprintf(stderr, "a1p12 = %20.5f, %20.5f\n",a1p12r,a1p12i);
	fprintf(stderr, "a1p24 = %20.5f, %20.5f\n",a1p24r,a1p24i);
	fprintf(stderr, "a1p36 = %20.5f, %20.5f\n",a1p36r,a1p36i);
	fprintf(stderr, "a1p04 = %20.5f, %20.5f\n",a1p04r,a1p04i);
	fprintf(stderr, "a1p16 = %20.5f, %20.5f\n",a1p16r,a1p16i);
	fprintf(stderr, "a1p28 = %20.5f, %20.5f\n",a1p28r,a1p28i);
	fprintf(stderr, "a1p40 = %20.5f, %20.5f\n",a1p40r,a1p40i);
	fprintf(stderr, "a1p08 = %20.5f, %20.5f\n",a1p08r,a1p08i);
	fprintf(stderr, "a1p20 = %20.5f, %20.5f\n",a1p20r,a1p20i);
	fprintf(stderr, "a1p32 = %20.5f, %20.5f\n",a1p32r,a1p32i);	fprintf(stderr, "\n");
	fprintf(stderr, "a1p11 = %20.5f, %20.5f\n",a1p11r,a1p11i);
	fprintf(stderr, "a1p23 = %20.5f, %20.5f\n",a1p23r,a1p23i);
	fprintf(stderr, "a1p35 = %20.5f, %20.5f\n",a1p35r,a1p35i);
	fprintf(stderr, "a1p03 = %20.5f, %20.5f\n",a1p03r,a1p03i);
	fprintf(stderr, "a1p15 = %20.5f, %20.5f\n",a1p15r,a1p15i);
	fprintf(stderr, "a1p27 = %20.5f, %20.5f\n",a1p27r,a1p27i);
	fprintf(stderr, "a1p39 = %20.5f, %20.5f\n",a1p39r,a1p39i);
	fprintf(stderr, "a1p07 = %20.5f, %20.5f\n",a1p07r,a1p07i);
	fprintf(stderr, "a1p19 = %20.5f, %20.5f\n",a1p19r,a1p19i);
	fprintf(stderr, "a1p31 = %20.5f, %20.5f\n",a1p31r,a1p31i);
	fprintf(stderr, "a1p43 = %20.5f, %20.5f\n",a1p43r,a1p43i);	fprintf(stderr, "\n");
	fprintf(stderr, "a1p22 = %20.5f, %20.5f\n",a1p22r,a1p22i);
	fprintf(stderr, "a1p34 = %20.5f, %20.5f\n",a1p34r,a1p34i);
	fprintf(stderr, "a1p02 = %20.5f, %20.5f\n",a1p02r,a1p02i);
	fprintf(stderr, "a1p14 = %20.5f, %20.5f\n",a1p14r,a1p14i);
	fprintf(stderr, "a1p26 = %20.5f, %20.5f\n",a1p26r,a1p26i);
	fprintf(stderr, "a1p38 = %20.5f, %20.5f\n",a1p38r,a1p38i);
	fprintf(stderr, "a1p06 = %20.5f, %20.5f\n",a1p06r,a1p06i);
	fprintf(stderr, "a1p18 = %20.5f, %20.5f\n",a1p18r,a1p18i);
	fprintf(stderr, "a1p30 = %20.5f, %20.5f\n",a1p30r,a1p30i);
	fprintf(stderr, "a1p42 = %20.5f, %20.5f\n",a1p42r,a1p42i);
	fprintf(stderr, "a1p10 = %20.5f, %20.5f\n",a1p10r,a1p10i);	fprintf(stderr, "\n");
	fprintf(stderr, "a1p33 = %20.5f, %20.5f\n",a1p33r,a1p33i);
	fprintf(stderr, "a1p01 = %20.5f, %20.5f\n",a1p01r,a1p01i);
	fprintf(stderr, "a1p13 = %20.5f, %20.5f\n",a1p13r,a1p13i);
	fprintf(stderr, "a1p25 = %20.5f, %20.5f\n",a1p25r,a1p25i);
	fprintf(stderr, "a1p37 = %20.5f, %20.5f\n",a1p37r,a1p37i);
	fprintf(stderr, "a1p05 = %20.5f, %20.5f\n",a1p05r,a1p05i);
	fprintf(stderr, "a1p17 = %20.5f, %20.5f\n",a1p17r,a1p17i);
	fprintf(stderr, "a1p29 = %20.5f, %20.5f\n",a1p29r,a1p29i);
	fprintf(stderr, "a1p41 = %20.5f, %20.5f\n",a1p41r,a1p41i);
	fprintf(stderr, "a1p09 = %20.5f, %20.5f\n",a1p09r,a1p09i);
	fprintf(stderr, "a1p21 = %20.5f, %20.5f\n",a1p21r,a1p21i);	fprintf(stderr, "\n");
#endif
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

				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy_r36,cy_r38,bjmodn36);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p40r,add1,add2,add3,cy_r40,cy_r42,bjmodn40);

		  #else	/* GCC-style inline ASM: */

				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p40r,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

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

				SSE2_cmplx_carry_norm_nocheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p36r,add1,add2,     cy_r36,cy_r38,bjmodn36);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p40r,add1,add2,     cy_r40,cy_r42,bjmodn40);

		  #else	/* GCC-style inline ASM: */

				SSE2_cmplx_carry_norm_nocheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p36r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p40r,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

		  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

		#else	/* #ifdef USE_SSE2 */

				/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);
				cmplx_carry_norm_errcheck(a1p28r,a1p28i,cy_r28,bjmodn28,28);
				cmplx_carry_norm_errcheck(a1p29r,a1p29i,cy_r29,bjmodn29,29);
				cmplx_carry_norm_errcheck(a1p30r,a1p30i,cy_r30,bjmodn30,30);
				cmplx_carry_norm_errcheck(a1p31r,a1p31i,cy_r31,bjmodn31,31);
				cmplx_carry_norm_errcheck(a1p32r,a1p32i,cy_r32,bjmodn32,32);
				cmplx_carry_norm_errcheck(a1p33r,a1p33i,cy_r33,bjmodn33,33);
				cmplx_carry_norm_errcheck(a1p34r,a1p34i,cy_r34,bjmodn34,34);
				cmplx_carry_norm_errcheck(a1p35r,a1p35i,cy_r35,bjmodn35,35);
				cmplx_carry_norm_errcheck(a1p36r,a1p36i,cy_r36,bjmodn36,36);
				cmplx_carry_norm_errcheck(a1p37r,a1p37i,cy_r37,bjmodn37,37);
				cmplx_carry_norm_errcheck(a1p38r,a1p38i,cy_r38,bjmodn38,38);
				cmplx_carry_norm_errcheck(a1p39r,a1p39i,cy_r39,bjmodn39,39);
				cmplx_carry_norm_errcheck(a1p40r,a1p40i,cy_r40,bjmodn40,40);
				cmplx_carry_norm_errcheck(a1p41r,a1p41i,cy_r41,bjmodn41,41);
				cmplx_carry_norm_errcheck(a1p42r,a1p42i,cy_r42,bjmodn42,42);
				cmplx_carry_norm_errcheck(a1p43r,a1p43i,cy_r43,bjmodn43,43);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	/* #ifdef USE_SSE2 */

			/*...The radix-44 DIF pass is here:	*/

		#ifdef DEBUG_SSE2
		  #ifdef USE_SSE2
			a1p00r = s1p00r->re;	a1p00i = s1p00i->re;
			a1p01r = s1p01r->re;	a1p01i = s1p01i->re;
			a1p02r = s1p02r->re;	a1p02i = s1p02i->re;
			a1p03r = s1p03r->re;	a1p03i = s1p03i->re;
			a1p04r = s1p04r->re;	a1p04i = s1p04i->re;
			a1p05r = s1p05r->re;	a1p05i = s1p05i->re;
			a1p06r = s1p06r->re;	a1p06i = s1p06i->re;
			a1p07r = s1p07r->re;	a1p07i = s1p07i->re;
			a1p08r = s1p08r->re;	a1p08i = s1p08i->re;
			a1p09r = s1p09r->re;	a1p09i = s1p09i->re;
			a1p10r = s1p10r->re;	a1p10i = s1p10i->re;
			a1p11r = s1p11r->re;	a1p11i = s1p11i->re;
			a1p12r = s1p12r->re;	a1p12i = s1p12i->re;
			a1p13r = s1p13r->re;	a1p13i = s1p13i->re;
			a1p14r = s1p14r->re;	a1p14i = s1p14i->re;
			a1p15r = s1p15r->re;	a1p15i = s1p15i->re;
			a1p16r = s1p16r->re;	a1p16i = s1p16i->re;
			a1p17r = s1p17r->re;	a1p17i = s1p17i->re;
			a1p18r = s1p18r->re;	a1p18i = s1p18i->re;
			a1p19r = s1p19r->re;	a1p19i = s1p19i->re;
			a1p20r = s1p20r->re;	a1p20i = s1p20i->re;
			a1p21r = s1p21r->re;	a1p21i = s1p21i->re;
			a1p22r = s1p22r->re;	a1p22i = s1p22i->re;
			a1p23r = s1p23r->re;	a1p23i = s1p23i->re;
			a1p24r = s1p24r->re;	a1p24i = s1p24i->re;
			a1p25r = s1p25r->re;	a1p25i = s1p25i->re;
			a1p26r = s1p26r->re;	a1p26i = s1p26i->re;
			a1p27r = s1p27r->re;	a1p27i = s1p27i->re;
			a1p28r = s1p28r->re;	a1p28i = s1p28i->re;
			a1p29r = s1p29r->re;	a1p29i = s1p29i->re;
			a1p30r = s1p30r->re;	a1p30i = s1p30i->re;
			a1p31r = s1p31r->re;	a1p31i = s1p31i->re;
			a1p32r = s1p32r->re;	a1p32i = s1p32i->re;
			a1p33r = s1p33r->re;	a1p33i = s1p33i->re;
			a1p34r = s1p34r->re;	a1p34i = s1p34i->re;
			a1p35r = s1p35r->re;	a1p35i = s1p35i->re;
			a1p36r = s1p36r->re;	a1p36i = s1p36i->re;
			a1p37r = s1p37r->re;	a1p37i = s1p37i->re;
			a1p38r = s1p38r->re;	a1p38i = s1p38i->re;
			a1p39r = s1p39r->re;	a1p39i = s1p39i->re;
			a1p40r = s1p40r->re;	a1p40i = s1p40i->re;
			a1p41r = s1p41r->re;	a1p41i = s1p41i->re;
			a1p42r = s1p42r->re;	a1p42i = s1p42i->re;
			a1p43r = s1p43r->re;	a1p43i = s1p43i->re;
		  #endif
			fprintf(stderr, "radix44_carry_out: s00= %20.5f, %20.5f\n",a1p00r,a1p00i);
			fprintf(stderr, "radix44_carry_out: s01= %20.5f, %20.5f\n",a1p01r,a1p01i);
			fprintf(stderr, "radix44_carry_out: s02= %20.5f, %20.5f\n",a1p02r,a1p02i);
			fprintf(stderr, "radix44_carry_out: s03= %20.5f, %20.5f\n",a1p03r,a1p03i);
			fprintf(stderr, "radix44_carry_out: s04= %20.5f, %20.5f\n",a1p04r,a1p04i);
			fprintf(stderr, "radix44_carry_out: s05= %20.5f, %20.5f\n",a1p05r,a1p05i);
			fprintf(stderr, "radix44_carry_out: s06= %20.5f, %20.5f\n",a1p06r,a1p06i);
			fprintf(stderr, "radix44_carry_out: s07= %20.5f, %20.5f\n",a1p07r,a1p07i);
			fprintf(stderr, "radix44_carry_out: s08= %20.5f, %20.5f\n",a1p08r,a1p08i);
			fprintf(stderr, "radix44_carry_out: s09= %20.5f, %20.5f\n",a1p09r,a1p09i);
			fprintf(stderr, "radix44_carry_out: s10= %20.5f, %20.5f\n",a1p10r,a1p10i);
			fprintf(stderr, "radix44_carry_out: s11= %20.5f, %20.5f\n",a1p11r,a1p11i);
			fprintf(stderr, "radix44_carry_out: s12= %20.5f, %20.5f\n",a1p12r,a1p12i);
			fprintf(stderr, "radix44_carry_out: s13= %20.5f, %20.5f\n",a1p13r,a1p13i);
			fprintf(stderr, "radix44_carry_out: s14= %20.5f, %20.5f\n",a1p14r,a1p14i);
			fprintf(stderr, "radix44_carry_out: s15= %20.5f, %20.5f\n",a1p15r,a1p15i);
			fprintf(stderr, "radix44_carry_out: s16= %20.5f, %20.5f\n",a1p16r,a1p16i);
			fprintf(stderr, "radix44_carry_out: s17= %20.5f, %20.5f\n",a1p17r,a1p17i);
			fprintf(stderr, "radix44_carry_out: s18= %20.5f, %20.5f\n",a1p18r,a1p18i);
			fprintf(stderr, "radix44_carry_out: s19= %20.5f, %20.5f\n",a1p19r,a1p19i);
			fprintf(stderr, "radix44_carry_out: s20= %20.5f, %20.5f\n",a1p20r,a1p20i);
			fprintf(stderr, "radix44_carry_out: s21= %20.5f, %20.5f\n",a1p21r,a1p21i);
			fprintf(stderr, "radix44_carry_out: s22= %20.5f, %20.5f\n",a1p22r,a1p22i);
			fprintf(stderr, "radix44_carry_out: s23= %20.5f, %20.5f\n",a1p23r,a1p23i);
			fprintf(stderr, "radix44_carry_out: s24= %20.5f, %20.5f\n",a1p24r,a1p24i);
			fprintf(stderr, "radix44_carry_out: s25= %20.5f, %20.5f\n",a1p25r,a1p25i);
			fprintf(stderr, "radix44_carry_out: s26= %20.5f, %20.5f\n",a1p26r,a1p26i);
			fprintf(stderr, "radix44_carry_out: s27= %20.5f, %20.5f\n",a1p27r,a1p27i);
			fprintf(stderr, "radix44_carry_out: s28= %20.5f, %20.5f\n",a1p28r,a1p28i);
			fprintf(stderr, "radix44_carry_out: s29= %20.5f, %20.5f\n",a1p29r,a1p29i);
			fprintf(stderr, "radix44_carry_out: s30= %20.5f, %20.5f\n",a1p30r,a1p30i);
			fprintf(stderr, "radix44_carry_out: s31= %20.5f, %20.5f\n",a1p31r,a1p31i);
			fprintf(stderr, "radix44_carry_out: s32= %20.5f, %20.5f\n",a1p32r,a1p32i);
			fprintf(stderr, "radix44_carry_out: s33= %20.5f, %20.5f\n",a1p33r,a1p33i);
			fprintf(stderr, "radix44_carry_out: s34= %20.5f, %20.5f\n",a1p34r,a1p34i);
			fprintf(stderr, "radix44_carry_out: s35= %20.5f, %20.5f\n",a1p35r,a1p35i);
			fprintf(stderr, "radix44_carry_out: s36= %20.5f, %20.5f\n",a1p36r,a1p36i);
			fprintf(stderr, "radix44_carry_out: s37= %20.5f, %20.5f\n",a1p37r,a1p37i);
			fprintf(stderr, "radix44_carry_out: s38= %20.5f, %20.5f\n",a1p38r,a1p38i);
			fprintf(stderr, "radix44_carry_out: s39= %20.5f, %20.5f\n",a1p39r,a1p39i);
			fprintf(stderr, "radix44_carry_out: s40= %20.5f, %20.5f\n",a1p40r,a1p40i);
			fprintf(stderr, "radix44_carry_out: s41= %20.5f, %20.5f\n",a1p41r,a1p41i);
			fprintf(stderr, "radix44_carry_out: s42= %20.5f, %20.5f\n",a1p42r,a1p42i);
			fprintf(stderr, "radix44_carry_out: s43= %20.5f, %20.5f\n",a1p43r,a1p43i);
		#endif

		#ifdef USE_SSE2

			/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-11 transforms...*/
			/* Radix-11 DFT inputs are (cyclic) with pXXr having XX -= 4 (4*16 bytes = 0x040), outputs are adjacent 32-byte-separated temps: */
							/*a1p00r,a1p40r,a1p36r,a1p32r,a1p28r,a1p24r,a1p20r,a1p16r,a1p12r,a1p08r,a1p04r */
			SSE2_RADIX_11_DFT(s1p00r, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
							/*a1p33r,a1p29r,a1p25r,a1p21r,a1p17r,a1p13r,a1p09r,a1p05r,a1p01r,a1p41r,a1p37r */
			SSE2_RADIX_11_DFT(s1p33r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300,-0x380,-0x400, 0x100, 0x080, ua0, t10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
							/*a1p22r,a1p18r,a1p14r,a1p10r,a1p06r,a1p02r,a1p42r,a1p38r,a1p34r,a1p30r,a1p26r */
			SSE2_RADIX_11_DFT(s1p22r,-0x080,-0x100,-0x180,-0x200,-0x280, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)
							/*a1p11r,a1p07r,a1p03r,a1p43r,a1p39r,a1p35r,a1p31r,a1p27r,a1p23r,a1p19r,a1p15r */
			SSE2_RADIX_11_DFT(s1p11r,-0x080,-0x100, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, ua0, t30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140)

			/*...and now do 11 radix-4 transforms...*/
			/* Inputs in SSE2 modes are temps 2*11*16 = 22*16 = 0x160 bytes apart: */
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

		#ifdef DEBUG_SSE2
			jt = j1;		jp = j2;
			fprintf(stderr, "radix44: A_out[00] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[01] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[02] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[03] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[04] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[05] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[06] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[07] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[08] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[09] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[10] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[11] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[12] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[13] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[14] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[15] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[16] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[17] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[18] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[29] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[20] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[21] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[22] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[23] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[24] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[25] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[26] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[27] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[28] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[29] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[30] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[31] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[32] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[33] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[34] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[35] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[36] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[37] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[38] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[39] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			fprintf(stderr, "radix44: A_out[40] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix44: A_out[41] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix44: A_out[42] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix44: A_out[43] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);	jt += p04;	jp += p04;
			exit(0);
		#endif
			}

			jstart += nwt;
			jhi    += nwt;

			col += 44;
			co3 -= 44;
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_SSE2
		_cy_r00[ithread] = cy_r00->re;
		_cy_r01[ithread] = cy_r00->im;
		_cy_r02[ithread] = cy_r02->re;
		_cy_r03[ithread] = cy_r02->im;
		_cy_r04[ithread] = cy_r04->re;
		_cy_r05[ithread] = cy_r04->im;
		_cy_r06[ithread] = cy_r06->re;
		_cy_r07[ithread] = cy_r06->im;
		_cy_r08[ithread] = cy_r08->re;
		_cy_r09[ithread] = cy_r08->im;
		_cy_r10[ithread] = cy_r10->re;
		_cy_r11[ithread] = cy_r10->im;
		_cy_r12[ithread] = cy_r12->re;
		_cy_r13[ithread] = cy_r12->im;
		_cy_r14[ithread] = cy_r14->re;
		_cy_r15[ithread] = cy_r14->im;
		_cy_r16[ithread] = cy_r16->re;
		_cy_r17[ithread] = cy_r16->im;
		_cy_r18[ithread] = cy_r18->re;
		_cy_r19[ithread] = cy_r18->im;
		_cy_r20[ithread] = cy_r20->re;
		_cy_r21[ithread] = cy_r20->im;
		_cy_r22[ithread] = cy_r22->re;
		_cy_r23[ithread] = cy_r22->im;
		_cy_r24[ithread] = cy_r24->re;
		_cy_r25[ithread] = cy_r24->im;
		_cy_r26[ithread] = cy_r26->re;
		_cy_r27[ithread] = cy_r26->im;
		_cy_r28[ithread] = cy_r28->re;
		_cy_r29[ithread] = cy_r28->im;
		_cy_r30[ithread] = cy_r30->re;
		_cy_r31[ithread] = cy_r30->im;
		_cy_r32[ithread] = cy_r32->re;
		_cy_r33[ithread] = cy_r32->im;
		_cy_r34[ithread] = cy_r34->re;
		_cy_r35[ithread] = cy_r34->im;
		_cy_r36[ithread] = cy_r36->re;
		_cy_r37[ithread] = cy_r36->im;
		_cy_r38[ithread] = cy_r38->re;
		_cy_r39[ithread] = cy_r38->im;
		_cy_r40[ithread] = cy_r40->re;
		_cy_r41[ithread] = cy_r40->im;
		_cy_r42[ithread] = cy_r42->re;
		_cy_r43[ithread] = cy_r42->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
	#else
		_cy_r00[ithread] = cy_r00;
		_cy_r01[ithread] = cy_r01;
		_cy_r02[ithread] = cy_r02;
		_cy_r03[ithread] = cy_r03;
		_cy_r04[ithread] = cy_r04;
		_cy_r05[ithread] = cy_r05;
		_cy_r06[ithread] = cy_r06;
		_cy_r07[ithread] = cy_r07;
		_cy_r08[ithread] = cy_r08;
		_cy_r09[ithread] = cy_r09;
		_cy_r10[ithread] = cy_r10;
		_cy_r11[ithread] = cy_r11;
		_cy_r12[ithread] = cy_r12;
		_cy_r13[ithread] = cy_r13;
		_cy_r14[ithread] = cy_r14;
		_cy_r15[ithread] = cy_r15;
		_cy_r16[ithread] = cy_r16;
		_cy_r17[ithread] = cy_r17;
		_cy_r18[ithread] = cy_r18;
		_cy_r19[ithread] = cy_r19;
		_cy_r20[ithread] = cy_r20;
		_cy_r21[ithread] = cy_r21;
		_cy_r22[ithread] = cy_r22;
		_cy_r23[ithread] = cy_r23;
		_cy_r24[ithread] = cy_r24;
		_cy_r25[ithread] = cy_r25;
		_cy_r26[ithread] = cy_r26;
		_cy_r27[ithread] = cy_r27;
		_cy_r28[ithread] = cy_r28;
		_cy_r29[ithread] = cy_r29;
		_cy_r30[ithread] = cy_r30;
		_cy_r31[ithread] = cy_r31;
		_cy_r32[ithread] = cy_r32;
		_cy_r33[ithread] = cy_r33;
		_cy_r34[ithread] = cy_r34;
		_cy_r35[ithread] = cy_r35;
		_cy_r36[ithread] = cy_r36;
		_cy_r37[ithread] = cy_r37;
		_cy_r38[ithread] = cy_r38;
		_cy_r39[ithread] = cy_r39;
		_cy_r40[ithread] = cy_r40;
		_cy_r41[ithread] = cy_r41;
		_cy_r42[ithread] = cy_r42;
		_cy_r43[ithread] = cy_r43;
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

	}	/******* END OF PARALLEL FOR-LOOP ********/

	if(!full_pass)break;

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-44 forward DIF FFT of the first block of 44 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 44 outputs of (1);
	!   (3) Reweight and perform a radix-44 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 44 elements and repeat (1-4).
	*/
	a1p00r= _cy_r00[CY_THREADS - 1];	/* Note base-11 indexing of temps here */
	a1p01r= _cy_r01[CY_THREADS - 1];
	a1p02r= _cy_r02[CY_THREADS - 1];
	a1p03r= _cy_r03[CY_THREADS - 1];
	a1p04r= _cy_r04[CY_THREADS - 1];
	a1p05r= _cy_r05[CY_THREADS - 1];
	a1p06r= _cy_r06[CY_THREADS - 1];
	a1p07r= _cy_r07[CY_THREADS - 1];
	a1p08r= _cy_r08[CY_THREADS - 1];
	a1p09r= _cy_r09[CY_THREADS - 1];
	a1p10r= _cy_r10[CY_THREADS - 1];
	a1p11r= _cy_r11[CY_THREADS - 1];
	a1p12r= _cy_r12[CY_THREADS - 1];
	a1p13r= _cy_r13[CY_THREADS - 1];
	a1p14r= _cy_r14[CY_THREADS - 1];
	a1p15r= _cy_r15[CY_THREADS - 1];
	a1p16r= _cy_r16[CY_THREADS - 1];
	a1p17r= _cy_r17[CY_THREADS - 1];
	a1p18r= _cy_r18[CY_THREADS - 1];
	a1p19r= _cy_r19[CY_THREADS - 1];
	a1p20r= _cy_r20[CY_THREADS - 1];
	a1p21r= _cy_r21[CY_THREADS - 1];
	a1p22r= _cy_r22[CY_THREADS - 1];
	a1p23r= _cy_r23[CY_THREADS - 1];
	a1p24r= _cy_r24[CY_THREADS - 1];
	a1p25r= _cy_r25[CY_THREADS - 1];
	a1p26r= _cy_r26[CY_THREADS - 1];
	a1p27r= _cy_r27[CY_THREADS - 1];
	a1p28r= _cy_r28[CY_THREADS - 1];
	a1p29r= _cy_r29[CY_THREADS - 1];
	a1p30r= _cy_r30[CY_THREADS - 1];
	a1p31r= _cy_r31[CY_THREADS - 1];
	a1p32r= _cy_r32[CY_THREADS - 1];
	a1p33r= _cy_r33[CY_THREADS - 1];
	a1p34r= _cy_r34[CY_THREADS - 1];
	a1p35r= _cy_r35[CY_THREADS - 1];
	a1p36r= _cy_r36[CY_THREADS - 1];
	a1p37r= _cy_r37[CY_THREADS - 1];
	a1p38r= _cy_r38[CY_THREADS - 1];
	a1p39r= _cy_r39[CY_THREADS - 1];
	a1p40r= _cy_r40[CY_THREADS - 1];
	a1p41r= _cy_r41[CY_THREADS - 1];
	a1p42r= _cy_r42[CY_THREADS - 1];
	a1p43r= _cy_r43[CY_THREADS - 1];

	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		ASSERT(HERE, CY_THREADS > 1,"radix44_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
		_cy_r00[ithread] = _cy_r00[ithread-1];
		_cy_r01[ithread] = _cy_r01[ithread-1];
		_cy_r02[ithread] = _cy_r02[ithread-1];
		_cy_r03[ithread] = _cy_r03[ithread-1];
		_cy_r04[ithread] = _cy_r04[ithread-1];
		_cy_r05[ithread] = _cy_r05[ithread-1];
		_cy_r06[ithread] = _cy_r06[ithread-1];
		_cy_r07[ithread] = _cy_r07[ithread-1];
		_cy_r08[ithread] = _cy_r08[ithread-1];
		_cy_r09[ithread] = _cy_r09[ithread-1];
		_cy_r10[ithread] = _cy_r10[ithread-1];
		_cy_r11[ithread] = _cy_r11[ithread-1];
		_cy_r12[ithread] = _cy_r12[ithread-1];
		_cy_r13[ithread] = _cy_r13[ithread-1];
		_cy_r14[ithread] = _cy_r14[ithread-1];
		_cy_r15[ithread] = _cy_r15[ithread-1];
		_cy_r16[ithread] = _cy_r16[ithread-1];
		_cy_r17[ithread] = _cy_r17[ithread-1];
		_cy_r18[ithread] = _cy_r18[ithread-1];
		_cy_r19[ithread] = _cy_r19[ithread-1];
		_cy_r20[ithread] = _cy_r20[ithread-1];
		_cy_r21[ithread] = _cy_r21[ithread-1];
		_cy_r22[ithread] = _cy_r22[ithread-1];
		_cy_r23[ithread] = _cy_r23[ithread-1];
		_cy_r24[ithread] = _cy_r24[ithread-1];
		_cy_r25[ithread] = _cy_r25[ithread-1];
		_cy_r26[ithread] = _cy_r26[ithread-1];
		_cy_r27[ithread] = _cy_r27[ithread-1];
		_cy_r28[ithread] = _cy_r28[ithread-1];
		_cy_r29[ithread] = _cy_r29[ithread-1];
		_cy_r30[ithread] = _cy_r30[ithread-1];
		_cy_r31[ithread] = _cy_r31[ithread-1];
		_cy_r32[ithread] = _cy_r32[ithread-1];
		_cy_r33[ithread] = _cy_r33[ithread-1];
		_cy_r34[ithread] = _cy_r34[ithread-1];
		_cy_r35[ithread] = _cy_r35[ithread-1];
		_cy_r36[ithread] = _cy_r36[ithread-1];
		_cy_r37[ithread] = _cy_r37[ithread-1];
		_cy_r38[ithread] = _cy_r38[ithread-1];
		_cy_r39[ithread] = _cy_r39[ithread-1];
		_cy_r40[ithread] = _cy_r40[ithread-1];
		_cy_r41[ithread] = _cy_r41[ithread-1];
		_cy_r42[ithread] = _cy_r42[ithread-1];
		_cy_r43[ithread] = _cy_r43[ithread-1];
	}

	_cy_r00[0] =+a1p43r;	/* ...The wraparound carry is here: */
	_cy_r01[0] = a1p00r;
	_cy_r02[0] = a1p01r;
	_cy_r03[0] = a1p02r;
	_cy_r04[0] = a1p03r;
	_cy_r05[0] = a1p04r;
	_cy_r06[0] = a1p05r;
	_cy_r07[0] = a1p06r;
	_cy_r08[0] = a1p07r;
	_cy_r09[0] = a1p08r;
	_cy_r10[0] = a1p09r;
	_cy_r11[0] = a1p10r;
	_cy_r12[0] = a1p11r;
	_cy_r13[0] = a1p12r;
	_cy_r14[0] = a1p13r;
	_cy_r15[0] = a1p14r;
	_cy_r16[0] = a1p15r;
	_cy_r17[0] = a1p16r;
	_cy_r18[0] = a1p17r;
	_cy_r19[0] = a1p18r;
	_cy_r20[0] = a1p19r;
	_cy_r21[0] = a1p20r;
	_cy_r22[0] = a1p21r;
	_cy_r23[0] = a1p22r;
	_cy_r24[0] = a1p23r;
	_cy_r25[0] = a1p24r;
	_cy_r26[0] = a1p25r;
	_cy_r27[0] = a1p26r;
	_cy_r28[0] = a1p27r;
	_cy_r29[0] = a1p28r;
	_cy_r30[0] = a1p29r;
	_cy_r31[0] = a1p30r;
	_cy_r32[0] = a1p31r;
	_cy_r33[0] = a1p32r;
	_cy_r34[0] = a1p33r;
	_cy_r35[0] = a1p34r;
	_cy_r36[0] = a1p35r;
	_cy_r37[0] = a1p36r;
	_cy_r38[0] = a1p37r;
	_cy_r39[0] = a1p38r;
	_cy_r40[0] = a1p39r;
	_cy_r41[0] = a1p40r;
	_cy_r42[0] = a1p41r;
	_cy_r43[0] = a1p42r;

	full_pass = 0;
	scale = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		j_jhi =15;
	}
	else
	{
		j_jhi = 7;
	}

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
		a1p00r += fabs(_cy_r00[0])+fabs(_cy_r01[0])+fabs(_cy_r02[0])+fabs(_cy_r03[0])+fabs(_cy_r04[0])+fabs(_cy_r05[0])+fabs(_cy_r06[0])+fabs(_cy_r07[0])+fabs(_cy_r08[0])+fabs(_cy_r09[0])+fabs(_cy_r10[0])+fabs(_cy_r11[0])+fabs(_cy_r12[0])+fabs(_cy_r13[0])+fabs(_cy_r14[0])+fabs(_cy_r15[0])+fabs(_cy_r16[0])+fabs(_cy_r17[0])+fabs(_cy_r18[0])+fabs(_cy_r19[0]);
		a1p00r += fabs(_cy_r20[0])+fabs(_cy_r21[0])+fabs(_cy_r22[0])+fabs(_cy_r23[0])+fabs(_cy_r24[0])+fabs(_cy_r25[0])+fabs(_cy_r26[0])+fabs(_cy_r27[0])+fabs(_cy_r28[0])+fabs(_cy_r29[0])+fabs(_cy_r30[0])+fabs(_cy_r31[0])+fabs(_cy_r32[0])+fabs(_cy_r33[0])+fabs(_cy_r34[0])+fabs(_cy_r35[0])+fabs(_cy_r36[0])+fabs(_cy_r37[0])+fabs(_cy_r38[0])+fabs(_cy_r39[0]);
		a1p00r += fabs(_cy_r40[0])+fabs(_cy_r41[0])+fabs(_cy_r42[0])+fabs(_cy_r43[0]);

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
		jt = j1; jp = j2;
		RADIX_11_DFT(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32]);	jt = j1+p03; jp = j2+p03;
		RADIX_11_DFT(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40]);	jt = j1+p02; jp = j2+p02;
		RADIX_11_DFT(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08]);	jt = j1+p01; jp = j2+p01;
		RADIX_11_DFT(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20]);
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

