/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef radix11_sse_macro_h_included
#define radix11_sse_macro_h_included

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_MSVC)

		/*...Radix-11 DFT: Inputs in memory locations __I0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA],\
		where r0 is a memory address and the i's are literal [byte] offsets. Outputs similarly go into memory locations\
		__O0 + [__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA], assumed disjoint with inputs:\

	            Total SSE opcounts: 118 MOVAPS, 182 ADD/SUBPD, 44 MULPD

		Mar 2014 UPDATE: Changed 64-bit GCC macros to use all memlocs, since byte offsets don't play nice with the kind
		of parametrized-pointer-offset approach we need for compact-object-code impl such as is used in radix 176.
		32-bit [both GCC and MSVC] no longer supported, so leave those macros as-was.
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

		// Simple 64-bit-ified version of the 32-bit ASM macro, using just xmm0-7
		#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
		{\
		__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
		/********************************************/\
		/*       Here are the 5 cosine terms:       */\
		/********************************************/\
			"movl	%[__cc],%%esi			\n\t"\
		"movl	%[__i1],%%eax	\n\t	movaps (%%eax),%%xmm1	\n\t"\
		"movl	%[__iA],%%ebx	\n\t	movaps (%%ebx),%%xmm5	\n\t"\
		"movl	%[__i2],%%ecx	\n\t	movaps (%%ecx),%%xmm2	\n\t"\
		"movl	%[__i9],%%edx	\n\t	movaps (%%edx),%%xmm6	\n\t"\
		"movl	%[__i3],%%eax	\n\t	movaps (%%eax),%%xmm3	\n\t"\
		"movl	%[__i8],%%ebx	\n\t	movaps (%%ebx),%%xmm7	\n\t"\
		"movl	%[__i4],%%ecx	\n\t	movaps (%%ecx),%%xmm4	\n\t"\
		"movl	%[__i7],%%edx	\n\t	movaps (%%edx),%%xmm0	\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t"\
			"subpd	%%xmm0,%%xmm4			\n\t"\
		"movl	%[__oA],%%eax	\n\t	movaps %%xmm1,0x10(%%eax)	\n\t"\
		"movl	%[__o9],%%ebx	\n\t	movaps %%xmm2,0x10(%%ebx)	\n\t"\
		"movl	%[__o8],%%ecx	\n\t	movaps %%xmm3,0x10(%%ecx)	\n\t"\
		"movl	%[__o7],%%edx	\n\t	movaps %%xmm4,0x10(%%edx)	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t"\
		"movl	%[__i5],%%eax	\n\t	movaps (%%eax),%%xmm5	\n\t"\
		"movl	%[__i6],%%ebx	\n\t	movaps (%%ebx),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"subpd	%%xmm6,%%xmm5			\n\t"\
		"movl	%[__o6],%%edi	\n\t	movaps %%xmm5,0x10(%%edi)	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movl	%[__I0],%%eax			\n\t"\
			"movaps	(%%eax),%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"/********************************************/\n\t"\
			"/*               Real Parts:                */\n\t"\
			"/********************************************/\n\t"\
		"movl	%[__o1],%%eax	\n\t"\
		"movl	%[__o2],%%ebx	\n\t"\
		"movl	%[__o3],%%ecx	\n\t"\
		"movl	%[__o4],%%edx	\n\t"\
		"movl	%[__o5],%%edi	\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"subpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%esi),%%xmm1		\n\t"\
			"movaps	%%xmm1,(%%eax)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0x10(%%esi),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%esi),%%xmm7		\n\t"\
			"movaps	%%xmm7,(%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%esi),%%xmm1		\n\t"\
			"movaps	%%xmm1,(%%ebx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%esi),%%xmm4		\n\t"\
			"mulpd	 0x70(%%esi),%%xmm5		\n\t"\
			"mulpd	 0x10(%%esi),%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%esi),%%xmm6		\n\t"\
			"mulpd	 0x50(%%esi),%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%esi),%%xmm1		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	 0x90(%%esi),%%xmm2		\n\t"\
			"movl	%[__O0],%%esi			\n\t"\
			"movaps	%%xmm0,(%%esi)			\n\t"\
			"movl	%[__cc],%%esi			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	(%%eax),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	(%%ebx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	(%%ecx),%%xmm6	\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm1,(%%eax)	\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"movaps	%%xmm7,(%%edi)	\n\t"\
			"movaps	%%xmm2,(%%ebx)	\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm3,(%%ecx)	\n\t"\
			"movaps	%%xmm4,(%%edx)	\n\t"\
			"/********************************************/\n\t"\
			"/*          Imaginary Parts:                */\n\t"\
			"/********************************************/\n\t"\
		"movl	%[__i1],%%eax	\n\t	movaps 0x10(%%eax),%%xmm1	\n\t"\
		"movl	%[__iA],%%ebx	\n\t	movaps 0x10(%%ebx),%%xmm5	\n\t"\
		"movl	%[__i2],%%ecx	\n\t	movaps 0x10(%%ecx),%%xmm2	\n\t"\
		"movl	%[__i9],%%edx	\n\t	movaps 0x10(%%edx),%%xmm6	\n\t"\
		"movl	%[__i3],%%eax	\n\t	movaps 0x10(%%eax),%%xmm3	\n\t"\
		"movl	%[__i8],%%ebx	\n\t	movaps 0x10(%%ebx),%%xmm7	\n\t"\
		"movl	%[__i4],%%ecx	\n\t	movaps 0x10(%%ecx),%%xmm4	\n\t"\
		"movl	%[__i7],%%edx	\n\t	movaps 0x10(%%edx),%%xmm0	\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm2			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t"\
			"subpd	%%xmm0,%%xmm4			\n\t"\
			/* Now this spill quartet is into Re slots of o-array: */\
		"movl	%[__oA],%%eax	\n\t	movaps %%xmm1,(%%eax)	\n\t"\
		"movl	%[__o9],%%ebx	\n\t	movaps %%xmm2,(%%ebx)	\n\t"\
		"movl	%[__o8],%%ecx	\n\t	movaps %%xmm3,(%%ecx)	\n\t"\
		"movl	%[__o7],%%edx	\n\t	movaps %%xmm4,(%%edx)	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t"\
		"movl	%[__i5],%%eax	\n\t	movaps 0x10(%%eax),%%xmm5	\n\t"\
		"movl	%[__i6],%%ebx	\n\t	movaps 0x10(%%ebx),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm3			\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"subpd	%%xmm6,%%xmm5			\n\t"\
		"movl	%[__o6],%%edi	\n\t	movaps %%xmm5,(%%edi)	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movl	%[__I0],%%eax			\n\t"\
			"movaps	0x10(%%eax),%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
		"movl	%[__o1],%%eax	\n\t"\
		"movl	%[__o2],%%ebx	\n\t"\
		"movl	%[__o3],%%ecx	\n\t"\
		"movl	%[__o4],%%edx	\n\t"\
		"movl	%[__o5],%%edi	\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"subpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%esi),%%xmm1		\n\t"\
			"movaps	%%xmm1,0x10(%%eax)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0x10(%%esi),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%esi),%%xmm7		\n\t"\
			"movaps	%%xmm7,0x10(%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%esi),%%xmm1		\n\t"\
			"movaps	%%xmm1,0x10(%%ebx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%esi),%%xmm4		\n\t"\
			"mulpd	 0x70(%%esi),%%xmm5		\n\t"\
			"mulpd	 0x10(%%esi),%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%esi),%%xmm6		\n\t"\
			"mulpd	 0x50(%%esi),%%xmm3		\n\t"\
			"addpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%esi),%%xmm1		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	 0x90(%%esi),%%xmm2		\n\t"\
			"movl	%[__O0],%%esi			\n\t"\
			"movaps	%%xmm0,0x10(%%esi)		\n\t"\
			"movl	%[__cc],%%esi			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	0x10(%%eax),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	0x10(%%ebx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	0x10(%%ecx),%%xmm6	\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"subpd	%%xmm7,%%xmm2			\n\t"\
			"addpd	%%xmm0,%%xmm7			\n\t"\
			"subpd	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm1,0x10(%%eax)	\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"movaps	%%xmm7,0x10(%%edi)	\n\t"\
			"movaps	%%xmm2,0x10(%%ebx)	\n\t"\
			"addpd	%%xmm0,%%xmm4			\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm4			\n\t"\
			"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
			"movaps	%%xmm4,0x10(%%edx)	\n\t"\
		/************************************************************************************************************************************/\
		/* Here are the 5 sine terms: Similar sequence to cosine terms, but with t21,19,17,15,13 replacing t3,5,7,9,11 and some sign flips: */\
		/************************************************************************************************************************************/\
			"addl	$0xa0,%%esi				\n\t"\
			"/********************************************/\n\t"\
			"/*               Real Parts:                */\n\t"\
			"/********************************************/\n\t"\
		"movl	%[__oA],%%eax	\n\t	movaps 0x10(%%eax),%%xmm1	\n\t"\
		"movl	%[__o9],%%ebx	\n\t	movaps 0x10(%%ebx),%%xmm2	\n\t"\
		"movl	%[__o8],%%ecx	\n\t	movaps 0x10(%%ecx),%%xmm3	\n\t"\
		"movl	%[__o7],%%edx	\n\t	movaps 0x10(%%edx),%%xmm4	\n\t"\
		"movl	%[__o6],%%edi	\n\t	movaps 0x10(%%edi),%%xmm5	\n\t"\
			/* r8-10 still free for o-addressing */\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%esi),%%xmm1		\n\t"\
			"movaps	%%xmm1,0x10(%%edi)	\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0xb0(%%esi),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%esi),%%xmm7		\n\t"\
			"movaps	%%xmm7,0x10(%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%esi),%%xmm1		\n\t"\
			"movaps	%%xmm1,0x10(%%edx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%esi),%%xmm4		\n\t"\
			"mulpd	 0x70(%%esi),%%xmm5		\n\t"\
			"mulpd	 0x10(%%esi),%%xmm7		\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%esi),%%xmm6		\n\t"\
			"mulpd	 0x50(%%esi),%%xmm3		\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%esi),%%xmm1		\n\t"\
			"mulpd	 0x90(%%esi),%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	0x10(%%edi),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	0x10(%%edx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	0x10(%%ecx),%%xmm6	\n\t"\
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
		"movl	%[__o1],%%esi	\n\t"\
			"addpd	0x10(%%esi),%%xmm1	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm1,0x10(%%esi)	\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"movaps	%%xmm1,0x10(%%eax)	\n\t"\
			"movaps	%%xmm0,%%xmm5			\n\t"\
		"movl	%[__o2],%%esi	\n\t"\
			"addpd	0x10(%%esi),%%xmm0	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm0,0x10(%%esi)	\n\t"\
			"subpd	%%xmm5,%%xmm0			\n\t"\
			"movaps	%%xmm0,0x10(%%ebx)	\n\t"\
			"movaps	%%xmm3,%%xmm6			\n\t"\
		"movl	%[__o3],%%esi	\n\t"\
			"addpd	0x10(%%esi),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	%%xmm3,0x10(%%esi)	\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm3,0x10(%%ecx)	\n\t"\
			"movaps	%%xmm4,%%xmm2			\n\t"\
		"movl	%[__o4],%%esi	\n\t"\
			"addpd	0x10(%%esi),%%xmm4	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm4,0x10(%%esi)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"movaps	%%xmm4,0x10(%%edx)	\n\t"\
			"movaps	%%xmm7,%%xmm5			\n\t"\
		"movl	%[__o5],%%esi	\n\t"\
			"addpd	0x10(%%esi),%%xmm7	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm7,0x10(%%esi)	\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"movaps	%%xmm7,0x10(%%edi)	\n\t"\
			"/********************************************/\n\t"\
			"/*          Imaginary Parts:                */\n\t"\
			"/********************************************/\n\t"\
			"movl	%[__cc],%%esi			\n\t"\
			"addl	$0xa0,%%esi				\n\t"\
		"movl	%[__oA],%%eax	\n\t	movaps (%%eax),%%xmm1	\n\t"\
		"movl	%[__o9],%%ebx	\n\t	movaps (%%ebx),%%xmm2	\n\t"\
		"movl	%[__o8],%%ecx	\n\t	movaps (%%ecx),%%xmm3	\n\t"\
		"movl	%[__o7],%%edx	\n\t	movaps (%%edx),%%xmm4	\n\t"\
		"movl	%[__o6],%%edi	\n\t	movaps (%%edi),%%xmm5	\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"addpd	%%xmm2,%%xmm5			\n\t"\
			"movaps	%%xmm1,%%xmm6			\n\t"\
			"addpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"mulpd	     (%%esi),%%xmm1		\n\t"\
			"movaps	%%xmm1,(%%edi)	\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	-0xb0(%%esi),%%xmm2		\n\t"\
			"addpd	%%xmm5,%%xmm6			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"mulpd	 0x60(%%esi),%%xmm7		\n\t"\
			"movaps	%%xmm7,(%%ecx)	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"mulpd	 0x30(%%esi),%%xmm1		\n\t"\
			"movaps	%%xmm1,(%%edx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"mulpd	 0x40(%%esi),%%xmm4		\n\t"\
			"mulpd	 0x70(%%esi),%%xmm5		\n\t"\
			"mulpd	 0x10(%%esi),%%xmm7		\n\t"\
			"subpd	%%xmm3,%%xmm2			\n\t"\
			"movaps	%%xmm6,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"mulpd	 0x80(%%esi),%%xmm6		\n\t"\
			"mulpd	 0x50(%%esi),%%xmm3		\n\t"\
			"subpd	%%xmm1,%%xmm2			\n\t"\
			"mulpd	 0x20(%%esi),%%xmm1		\n\t"\
			"mulpd	 0x90(%%esi),%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	(%%edi),%%xmm1	\n\t"\
			"addpd	%%xmm3,%%xmm4			\n\t"\
			"addpd	(%%edx),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm5			\n\t"\
			"addpd	(%%ecx),%%xmm6	\n\t"\
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
		"movl	%[__o1],%%esi	\n\t"\
			"addpd	(%%esi),%%xmm1	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm1,(%%eax)	\n\t"\
			"subpd	%%xmm2,%%xmm1			\n\t"\
			"movaps	%%xmm1,(%%esi)	\n\t"\
			"movaps	%%xmm0,%%xmm5			\n\t"\
		"movl	%[__o2],%%esi	\n\t"\
			"addpd	(%%esi),%%xmm0	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm0,(%%ebx)	\n\t"\
			"subpd	%%xmm5,%%xmm0			\n\t"\
			"movaps	%%xmm0,(%%esi)	\n\t"\
			"movaps	%%xmm3,%%xmm6			\n\t"\
		"movl	%[__o3],%%esi	\n\t"\
			"addpd	(%%esi),%%xmm3	\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"movaps	%%xmm3,(%%ecx)	\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"movaps	%%xmm3,(%%esi)	\n\t"\
			"movaps	%%xmm4,%%xmm2			\n\t"\
		"movl	%[__o4],%%esi	\n\t"\
			"addpd	(%%esi),%%xmm4	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm4,(%%edx)	\n\t"\
			"subpd	%%xmm2,%%xmm4			\n\t"\
			"movaps	%%xmm4,(%%esi)	\n\t"\
			"movaps	%%xmm7,%%xmm5			\n\t"\
		"movl	%[__o5],%%esi	\n\t"\
			"addpd	(%%esi),%%xmm7	\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"movaps	%%xmm7,(%%edi)	\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"movaps	%%xmm7,(%%esi)	\n\t"\
	"popl %%ebx	\n\t"\
			:					/* outputs: none */\
			: [__I0] "m" (XI0)	/* All inputs from memory addresses here */\
			 ,[__i1] "m" (Xi1)\
			 ,[__i2] "m" (Xi2)\
			 ,[__i3] "m" (Xi3)\
			 ,[__i4] "m" (Xi4)\
			 ,[__i5] "m" (Xi5)\
			 ,[__i6] "m" (Xi6)\
			 ,[__i7] "m" (Xi7)\
			 ,[__i8] "m" (Xi8)\
			 ,[__i9] "m" (Xi9)\
			 ,[__iA] "m" (XiA)\
			 ,[__cc] "m" (Xcc)\
			 ,[__O0] "m" (XO0)\
			 ,[__o1] "m" (Xo1)\
			 ,[__o2] "m" (Xo2)\
			 ,[__o3] "m" (Xo3)\
			 ,[__o4] "m" (Xo4)\
			 ,[__o5] "m" (Xo5)\
			 ,[__o6] "m" (Xo6)\
			 ,[__o7] "m" (Xo7)\
			 ,[__o8] "m" (Xo8)\
			 ,[__o9] "m" (Xo9)\
			 ,[__oA] "m" (XoA)\
			: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","edi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
			);\
		}

	  #else

		#ifdef USE_AVX
		//
		// NB: GCC builds gave "Error: invalid character '(' in mnemonic" for this macro ... I'd forgotten
		// to restore whites to several ...pd (... pairings after having removed it for purposes of col-
		// alignment during AVX-ification of what started as the SSE2 version of the macro:
		//
		// ******> A good way to track down such non-obvious syntax errors output by the assembler (typically w/o
		// helpful localization information such as line numbers) is to insert a single asm line with a different
		// kind of asm error - e.g. a vxorpd with just 2 register args in this case - and move 1 or 2 such around
		// in order to ever-more-closely bracket the first error of the type whose origin is being sought.
		//
			#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
			{\
			__asm__ volatile (\
				"/********************************************/\n\t"\
				"/*       Here are the 5 cosine terms:       */\n\t"\
				"/********************************************/\n\t"\
				"movq	%[__I0],%%rax			\n\t			movq	%[__O0],%%rcx		/* rax/rcx point to Re parts of I/Os */	\n\t"\
				"leaq	0x20(%%rax),%%rbx		\n\t			leaq	0x20(%%rcx),%%rdx	/* rbx/rdx point to Im parts of I/Os */	\n\t"\
				"movq	%[__cc],%%rsi			\n\t"\
				"/*****************************/\n\t			/******************************/\n\t"\
				"/*        Real Parts:         /\n\t			/*     Imaginary Parts:        /\n\t"\
				"/*****************************/\n\t			/******************************/\n\t"\
			"movq	%[__i1],%%r8 	\n\t	vmovaps (%%r8 ),%%ymm1	\n\t	vmovaps	0x20(%%r8 ),%%ymm9 	\n\t"\
			"movq	%[__iA],%%r9 	\n\t	vmovaps (%%r9 ),%%ymm5	\n\t	vmovaps	0x20(%%r9 ),%%ymm13	\n\t"\
			"movq	%[__i2],%%r10	\n\t	vmovaps (%%r10),%%ymm2	\n\t	vmovaps	0x20(%%r10),%%ymm10	\n\t"\
			"movq	%[__i9],%%r11	\n\t	vmovaps (%%r11),%%ymm6	\n\t	vmovaps	0x20(%%r11),%%ymm14	\n\t"\
			"movq	%[__i3],%%r12	\n\t	vmovaps (%%r12),%%ymm3	\n\t	vmovaps	0x20(%%r12),%%ymm11	\n\t"\
			"movq	%[__i8],%%r13	\n\t	vmovaps (%%r13),%%ymm7	\n\t	vmovaps	0x20(%%r13),%%ymm15	\n\t"\
			"movq	%[__i4],%%r14	\n\t	vmovaps (%%r14),%%ymm4	\n\t	vmovaps	0x20(%%r14),%%ymm12	\n\t"\
			"movq	%[__i7],%%r15	\n\t	vmovaps (%%r15),%%ymm0	\n\t	vmovaps	0x20(%%r15),%%ymm8 	\n\t"\
				"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vsubpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
				"vsubpd	%%ymm6,%%ymm2,%%ymm2	\n\t			vsubpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
				"vsubpd	%%ymm7,%%ymm3,%%ymm3	\n\t			vsubpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
				"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vsubpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
			"movq	%[__oA],%%r8 	\n\t	vmovaps %%ymm1,0x20(%%r8 )	\n\t	vmovaps %%ymm9 ,(%%r8 )	\n\t"\
			"movq	%[__o9],%%r9 	\n\t	vmovaps %%ymm2,0x20(%%r9 )	\n\t	vmovaps %%ymm10,(%%r9 )	\n\t"\
			"movq	%[__o8],%%r10	\n\t	vmovaps %%ymm3,0x20(%%r10)	\n\t	vmovaps %%ymm11,(%%r10)	\n\t"\
			"movq	%[__o7],%%r11	\n\t	vmovaps %%ymm4,0x20(%%r11)	\n\t	vmovaps %%ymm12,(%%r11)	\n\t"\
				"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"\
				"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
				"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t			vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
				"vaddpd	%%ymm0,%%ymm0,%%ymm0	\n\t			vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 			\n\t"\
				"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vaddpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
				"vaddpd	%%ymm6,%%ymm2,%%ymm2	\n\t			vaddpd	%%ymm14,%%ymm10,%%ymm10			\n\t"\
			"movq	%[__i5],%%r9 	\n\t	vmovaps (%%r9 ),%%ymm5	\n\t	vmovaps 0x20(%%r9 ),%%ymm13	\n\t"\
			"movq	%[__i6],%%r11	\n\t	vmovaps (%%r11),%%ymm6	\n\t	vmovaps 0x20(%%r11),%%ymm14	\n\t"\
				"vaddpd	%%ymm7,%%ymm3,%%ymm3	\n\t			vaddpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
				"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
				"vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t			vsubpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
			"movq	%[__o6],%%r8 	\n\t	vmovaps %%ymm5,0x20(%%r8 )	\n\t	vmovaps %%ymm13,(%%r8 )	\n\t"\
				"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
				"vmovaps	(%%rax),%%ymm0		\n\t			vmovaps	(%%rbx),%%ymm8 					\n\t"\
				"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t			vaddpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
			"movq	%[__o1],%%r11	\n\t"\
			"movq	%[__o2],%%r12	\n\t"\
			"movq	%[__o3],%%r13	\n\t"\
			"movq	%[__o4],%%r14	\n\t"\
			"movq	%[__o5],%%r15	\n\t"\
				"vsubpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
				"vsubpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
				"vmovaps	%%ymm1,%%ymm6				\n\t	vmovaps	%%ymm9 ,%%ymm14					\n\t"\
				"vsubpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm10,%%ymm11,%%ymm11			\n\t"\
				"vmovaps	%%ymm1,%%ymm7				\n\t	vmovaps	%%ymm9 ,%%ymm15					\n\t"\
				"vmulpd	      (%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	      (%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
				"vmovaps	%%ymm1,(%%r11)				\n\t	vmovaps	%%ymm9 ,0x20(%%r11)				\n\t"\
				"vsubpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
				"vmulpd	-0x020(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	-0x020(%%rsi),%%ymm10,%%ymm10	\n\t"\
				"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14			\n\t"\
				"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
				"vmulpd	 0x0c0(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"\
				"vmovaps	%%ymm7,(%%r13)				\n\t	vmovaps	%%ymm15,0x20(%%r13)				\n\t"\
				"vmovaps	%%ymm3,%%ymm1				\n\t	vmovaps	%%ymm11,%%ymm9 					\n\t"\
				"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11			\n\t"\
				"vmulpd	 0x060(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x060(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
				"vmovaps	%%ymm1,(%%r12)				\n\t	vmovaps	%%ymm9 ,0x20(%%r12)				\n\t"\
				"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15					\n\t"\
				"vsubpd	%%ymm4,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
				"vmulpd	 0x080(%%rsi),%%ymm4,%%ymm4		\n\t	vmulpd	 0x080(%%rsi),%%ymm12,%%ymm12	\n\t"\
				"vmulpd	 0x0e0(%%rsi),%%ymm5,%%ymm5		\n\t	vmulpd	 0x0e0(%%rsi),%%ymm13,%%ymm13	\n\t"\
				"vmulpd	 0x020(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x020(%%rsi),%%ymm15,%%ymm15	\n\t"\
				"vaddpd	%%ymm3,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm11,%%ymm10,%%ymm10			\n\t"\
				"vmovaps	%%ymm6,%%ymm1				\n\t	vmovaps	%%ymm14,%%ymm9 					\n\t"\
				"vsubpd	%%ymm3,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
				"vmulpd	 0x100(%%rsi),%%ymm6,%%ymm6		\n\t	vmulpd	 0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
				"vmulpd	 0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t	vmulpd	 0x0a0(%%rsi),%%ymm11,%%ymm11	\n\t"\
				"vaddpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"\
				"vmulpd	 0x040(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x040(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
				"vmulpd	 0x120(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	 0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
				"vmovaps	%%ymm0,(%%rcx)				\n\t	vmovaps	%%ymm8 ,(%%rdx)					\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
				"vaddpd	%%ymm1,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
				"vaddpd	(%%r11),%%ymm1,%%ymm1			\n\t	vaddpd	0x20(%%r11),%%ymm9 ,%%ymm9 		\n\t"\
				"vaddpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
				"vaddpd	(%%r12),%%ymm3,%%ymm3			\n\t	vaddpd	0x20(%%r12),%%ymm11,%%ymm11		\n\t"\
				"vaddpd	%%ymm6,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
				"vaddpd	(%%r13),%%ymm6,%%ymm6			\n\t	vaddpd	0x20(%%r13),%%ymm14,%%ymm14		\n\t"\
				"vmovaps	%%ymm2,%%ymm0				\n\t	vmovaps	%%ymm10,%%ymm8 					\n\t"\
				"vsubpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"\
				"vaddpd	%%ymm0,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 			\n\t"\
				"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
				"vaddpd	%%ymm0,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
				"vsubpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
				"vsubpd	%%ymm3,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10			\n\t"\
				"vmovaps	%%ymm1,(%%r11)				\n\t	vmovaps	%%ymm9 ,0x20(%%r11)				\n\t"\
				"vaddpd	%%ymm0,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm8 ,%%ymm11,%%ymm11			\n\t"\
				"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
				"vsubpd	%%ymm4,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
				"vmovaps	%%ymm7,(%%r15)				\n\t	vmovaps	%%ymm15,0x20(%%r15)				\n\t"\
				"vmovaps	%%ymm2,(%%r12)				\n\t	vmovaps	%%ymm10,0x20(%%r12)				\n\t"\
				"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
				"vaddpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
				"vaddpd	%%ymm5,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm13,%%ymm12,%%ymm12			\n\t"\
				"vmovaps	%%ymm3,(%%r13)				\n\t	vmovaps	%%ymm11,0x20(%%r13)				\n\t"\
				"vmovaps	%%ymm4,(%%r14)				\n\t	vmovaps	%%ymm12,0x20(%%r14)				\n\t"\
			"/***********************************************************************************************************/\n\t"\
			"/* Here are 5 sine terms: Similar to cosines, but t21,19,17,15,13 replace t3,5,7,9,11 and some sign flips: */\n\t"\
			"/***********************************************************************************************************/\n\t"\
				"addq	$0x140,%%rsi				\n\t"\
				"/*****************************/\n\t			/******************************/\n\t"\
				"/*        Real Parts:         /\n\t			/*     Imaginary Parts:        /\n\t"\
				"/*****************************/\n\t			/******************************/\n\t"\
			"movq	%[__oA],%%r11	\n\t	vmovaps (%%r11),%%ymm1	\n\t	vmovaps 0x20(%%r11),%%ymm9 	\n\t"\
			"movq	%[__o9],%%r12	\n\t	vmovaps (%%r12),%%ymm2	\n\t	vmovaps 0x20(%%r12),%%ymm10	\n\t"\
			"movq	%[__o8],%%r13	\n\t	vmovaps (%%r13),%%ymm3	\n\t	vmovaps 0x20(%%r13),%%ymm11	\n\t"\
			"movq	%[__o7],%%r14	\n\t	vmovaps (%%r14),%%ymm4	\n\t	vmovaps 0x20(%%r14),%%ymm12	\n\t"\
			"movq	%[__o6],%%r15	\n\t	vmovaps (%%r15),%%ymm5	\n\t	vmovaps 0x20(%%r15),%%ymm13	\n\t"\
				"vaddpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
				"vaddpd	%%ymm2,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
				"vmovaps	%%ymm1,%%ymm6				\n\t	vmovaps	%%ymm9 ,%%ymm14					\n\t"\
				"vaddpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11			\n\t"\
				"vmovaps	%%ymm1,%%ymm7				\n\t	vmovaps	%%ymm9 ,%%ymm15					\n\t"\
				"vmulpd	      (%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	      (%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
				"vmovaps	%%ymm1,(%%r15)				\n\t	vmovaps	%%ymm9 ,0x20(%%r15)				\n\t"\
				"vaddpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
				"vmulpd	-0x160(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	-0x160(%%rsi),%%ymm10,%%ymm10	\n\t"\
				"vaddpd	%%ymm5,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm13,%%ymm14,%%ymm14			\n\t"\
				"vsubpd	%%ymm3,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
				"vmulpd	 0x0c0(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"\
				"vmovaps	%%ymm7,(%%r13)				\n\t	vmovaps	%%ymm15,0x20(%%r13)				\n\t"\
				"vmovaps	%%ymm3,%%ymm1				\n\t	vmovaps	%%ymm11,%%ymm9 					\n\t"\
				"vaddpd	%%ymm4,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm12,%%ymm11,%%ymm11			\n\t"\
				"vmulpd	 0x060(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x060(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
				"vmovaps	%%ymm1,(%%r14)				\n\t	vmovaps	%%ymm9 ,0x20(%%r14)				\n\t"\
				"vmovaps	%%ymm5,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15					\n\t"\
				"vsubpd	%%ymm4,%%ymm5,%%ymm5			\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
				"vmulpd	 0x080(%%rsi),%%ymm4,%%ymm4		\n\t	vmulpd	 0x080(%%rsi),%%ymm12,%%ymm12	\n\t"\
				"vmulpd	 0x0e0(%%rsi),%%ymm5,%%ymm5		\n\t	vmulpd	 0x0e0(%%rsi),%%ymm13,%%ymm13	\n\t"\
				"vmulpd	 0x020(%%rsi),%%ymm7,%%ymm7		\n\t	vmulpd	 0x020(%%rsi),%%ymm15,%%ymm15	\n\t"\
				"vsubpd	%%ymm3,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm11,%%ymm10,%%ymm10			\n\t"\
				"vmovaps	%%ymm6,%%ymm1				\n\t	vmovaps	%%ymm14,%%ymm9 					\n\t"\
				"vsubpd	%%ymm3,%%ymm6,%%ymm6			\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
				"vmulpd	 0x100(%%rsi),%%ymm6,%%ymm6		\n\t	vmulpd	 0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
				"vmulpd	 0x0a0(%%rsi),%%ymm3,%%ymm3		\n\t	vmulpd	 0x0a0(%%rsi),%%ymm11,%%ymm11	\n\t"\
				"vsubpd	%%ymm1,%%ymm2,%%ymm2			\n\t	vsubpd	%%ymm9 ,%%ymm10,%%ymm10			\n\t"\
				"vmulpd	 0x040(%%rsi),%%ymm1,%%ymm1		\n\t	vmulpd	 0x040(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
				"vmulpd	 0x120(%%rsi),%%ymm2,%%ymm2		\n\t	vmulpd	 0x120(%%rsi),%%ymm10,%%ymm10	\n\t"\
				"vaddpd	%%ymm1,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
				"vaddpd	(%%r15),%%ymm1,%%ymm1			\n\t	vaddpd	0x20(%%r15),%%ymm9 ,%%ymm9 		\n\t"\
				"vaddpd	%%ymm3,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
				"vaddpd	(%%r14),%%ymm3,%%ymm3			\n\t	vaddpd	0x20(%%r14),%%ymm11,%%ymm11		\n\t"\
				"vaddpd	%%ymm6,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
				"vaddpd	(%%r13),%%ymm6,%%ymm6			\n\t	vaddpd	0x20(%%r13),%%ymm14,%%ymm14		\n\t"\
				"vxorpd	%%ymm0,%%ymm0,%%ymm0			\n\t	vxorpd	%%ymm8 ,%%ymm8 ,%%ymm8			\n\t"\
				"vsubpd	%%ymm2,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
				"vaddpd	%%ymm1,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm9 ,%%ymm8 ,%%ymm8 			\n\t"\
				"vaddpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
				"vaddpd	%%ymm7,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
				"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
				"vsubpd	%%ymm6,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
				"vaddpd	%%ymm3,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
				"vaddpd	%%ymm2,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm10,%%ymm11,%%ymm11			\n\t"\
				"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
				"vaddpd	%%ymm4,%%ymm0,%%ymm0			\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
				"vaddpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
				"vaddpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
				"vaddpd	%%ymm5,%%ymm4,%%ymm4			\n\t	vaddpd	%%ymm13,%%ymm12,%%ymm12			\n\t"\
				"vmovaps	%%ymm1,%%ymm2				\n\t	vmovaps	%%ymm9 ,%%ymm10					\n\t"\
			"movq	%[__o1],%%r8 					\n\t"\
				"vaddpd	(%%r8 ),%%ymm1,%%ymm1			\n\t	vaddpd	0x20(%%r8 ),%%ymm9 ,%%ymm9 		\n\t"\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10			\n\t"\
				"vmovaps	%%ymm1,(%%r11)				\n\t	vmovaps	%%ymm9 ,0x20(%%r8 )				\n\t"\
				"vsubpd	%%ymm2,%%ymm1,%%ymm1			\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
				"vmovaps	%%ymm1,(%%r8 )				\n\t	vmovaps	%%ymm9 ,0x20(%%r11)				\n\t"\
				"vmovaps	%%ymm0,%%ymm5				\n\t	vmovaps	%%ymm8 ,%%ymm13					\n\t"\
			"movq	%[__o2],%%r8 					\n\t"\
				"vaddpd	(%%r8 ),%%ymm0,%%ymm0			\n\t	vaddpd	0x20(%%r8 ),%%ymm8 ,%%ymm8 		\n\t"\
				"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"\
				"vmovaps	%%ymm0,(%%r12)				\n\t	vmovaps	%%ymm8 ,0x20(%%r8 )				\n\t"\
				"vsubpd	%%ymm5,%%ymm0,%%ymm0			\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 			\n\t"\
				"vmovaps	%%ymm0,(%%r8 )				\n\t	vmovaps	%%ymm8 ,0x20(%%r12)				\n\t"\
				"vmovaps	%%ymm3,%%ymm6				\n\t	vmovaps	%%ymm11,%%ymm14					\n\t"\
			"movq	%[__o3],%%r8 					\n\t"\
				"vaddpd	(%%r8 ),%%ymm3,%%ymm3			\n\t	vaddpd	0x20(%%r8 ),%%ymm11,%%ymm11		\n\t"\
				"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
				"vmovaps	%%ymm3,(%%r13)				\n\t	vmovaps	%%ymm11,0x20(%%r8 )				\n\t"\
				"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11			\n\t"\
				"vmovaps	%%ymm3,(%%r8 )				\n\t	vmovaps	%%ymm11,0x20(%%r13)				\n\t"\
				"vmovaps	%%ymm4,%%ymm2				\n\t	vmovaps	%%ymm12,%%ymm10					\n\t"\
			"movq	%[__o4],%%r8 					\n\t"\
				"vaddpd	(%%r8 ),%%ymm4,%%ymm4			\n\t	vaddpd	0x20(%%r8 ),%%ymm12,%%ymm12		\n\t"\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10			\n\t"\
				"vmovaps	%%ymm4,(%%r14)				\n\t	vmovaps	%%ymm12,0x20(%%r8 )				\n\t"\
				"vsubpd	%%ymm2,%%ymm4,%%ymm4			\n\t	vsubpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
				"vmovaps	%%ymm4,(%%r8 )				\n\t	vmovaps	%%ymm12,0x20(%%r14)				\n\t"\
				"vmovaps	%%ymm7,%%ymm5				\n\t	vmovaps	%%ymm15,%%ymm13					\n\t"\
			"movq	%[__o5],%%r8 					\n\t"\
				"vaddpd	(%%r8 ),%%ymm7,%%ymm7			\n\t	vaddpd	0x20(%%r8 ),%%ymm15,%%ymm15		\n\t"\
				"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"\
				"vmovaps	%%ymm7,(%%r15)				\n\t	vmovaps	%%ymm15,0x20(%%r8 )				\n\t"\
				"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
				"vmovaps	%%ymm7,(%%r8 )				\n\t	vmovaps	%%ymm15,0x20(%%r15)				\n\t"\
				:					/* outputs: none */\
				: [__I0] "m" (XI0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__i9] "m" (Xi9)\
				 ,[__iA] "m" (XiA)\
				 ,[__cc] "m" (Xcc)\
				 ,[__O0] "m" (XO0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				 ,[__o9] "m" (Xo9)\
				 ,[__oA] "m" (XoA)\
				: "cc","memory","rax","rbx","rcx","rdx","rsi","r8","r9","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
				);\
			}

		#else	// 64-bit SSE2:

		  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate 16-reg version of radix-11 DFT is faster

		  #if USE_64BIT_ASM_STYLE

			// 2-Instruction-stream-overlapped 64-bit-ified version of the 32-bit-style ASM macro, using all of xmm0-15
			#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
			{\
			__asm__ volatile (\
				"/********************************************/\n\t"\
				"/*       Here are the 5 cosine terms:       */\n\t"\
				"/********************************************/\n\t"\
				"movq	%[__I0],%%rax			\n\t			movq	%[__O0],%%rcx		/* rax/rcx point to Re parts of I/Os */	\n\t"\
				"leaq	0x10(%%rax),%%rbx		\n\t			leaq	0x10(%%rcx),%%rdx	/* rbx/rdx point to Im parts of I/Os */	\n\t"\
				"movq	%[__cc],%%rsi			\n\t"\
				"/*****************************/\n\t			/******************************/\n\t"\
				"/*        Real Parts:         /\n\t			/*     Imaginary Parts:        /\n\t"\
				"/*****************************/\n\t			/******************************/\n\t"\
			"movq	%[__i1],%%r8 	\n\t	movaps (%%r8 ),%%xmm1	\n\t	movaps	0x10(%%r8 ),%%xmm9 	\n\t"\
			"movq	%[__iA],%%r9 	\n\t	movaps (%%r9 ),%%xmm5	\n\t	movaps	0x10(%%r9 ),%%xmm13	\n\t"\
			"movq	%[__i2],%%r10	\n\t	movaps (%%r10),%%xmm2	\n\t	movaps	0x10(%%r10),%%xmm10	\n\t"\
			"movq	%[__i9],%%r11	\n\t	movaps (%%r11),%%xmm6	\n\t	movaps	0x10(%%r11),%%xmm14	\n\t"\
			"movq	%[__i3],%%r12	\n\t	movaps (%%r12),%%xmm3	\n\t	movaps	0x10(%%r12),%%xmm11	\n\t"\
			"movq	%[__i8],%%r13	\n\t	movaps (%%r13),%%xmm7	\n\t	movaps	0x10(%%r13),%%xmm15	\n\t"\
			"movq	%[__i4],%%r14	\n\t	movaps (%%r14),%%xmm4	\n\t	movaps	0x10(%%r14),%%xmm12	\n\t"\
			"movq	%[__i7],%%r15	\n\t	movaps (%%r15),%%xmm0	\n\t	movaps	0x10(%%r15),%%xmm8 	\n\t"\
				"subpd	%%xmm5,%%xmm1			\n\t			subpd	%%xmm13,%%xmm9 			\n\t"\
				"subpd	%%xmm6,%%xmm2			\n\t			subpd	%%xmm14,%%xmm10			\n\t"\
				"subpd	%%xmm7,%%xmm3			\n\t			subpd	%%xmm15,%%xmm11			\n\t"\
				"subpd	%%xmm0,%%xmm4			\n\t			subpd	%%xmm8 ,%%xmm12			\n\t"\
			"movq	%[__oA],%%r8 	\n\t	movaps %%xmm1,0x10(%%r8 )	\n\t	movaps %%xmm9 ,(%%r8 )	\n\t"\
			"movq	%[__o9],%%r9 	\n\t	movaps %%xmm2,0x10(%%r9 )	\n\t	movaps %%xmm10,(%%r9 )	\n\t"\
			"movq	%[__o8],%%r10	\n\t	movaps %%xmm3,0x10(%%r10)	\n\t	movaps %%xmm11,(%%r10)	\n\t"\
			"movq	%[__o7],%%r11	\n\t	movaps %%xmm4,0x10(%%r11)	\n\t	movaps %%xmm12,(%%r11)	\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t			addpd	%%xmm13,%%xmm13			\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t			addpd	%%xmm14,%%xmm14			\n\t"\
				"addpd	%%xmm7,%%xmm7			\n\t			addpd	%%xmm15,%%xmm15			\n\t"\
				"addpd	%%xmm0,%%xmm0			\n\t			addpd	%%xmm8 ,%%xmm8 			\n\t"\
				"addpd	%%xmm5,%%xmm1			\n\t			addpd	%%xmm13,%%xmm9 			\n\t"\
				"addpd	%%xmm6,%%xmm2			\n\t			addpd	%%xmm14,%%xmm10			\n\t"\
			"movq	%[__i5],%%r9 	\n\t	movaps (%%r9 ),%%xmm5	\n\t	movaps 0x10(%%r9 ),%%xmm13	\n\t"\
			"movq	%[__i6],%%r11	\n\t	movaps (%%r11),%%xmm6	\n\t	movaps 0x10(%%r11),%%xmm14	\n\t"\
				"addpd	%%xmm7,%%xmm3			\n\t			addpd	%%xmm15,%%xmm11			\n\t"\
				"addpd	%%xmm0,%%xmm4			\n\t			addpd	%%xmm8 ,%%xmm12			\n\t"\
				"subpd	%%xmm6,%%xmm5			\n\t			subpd	%%xmm14,%%xmm13			\n\t"\
			"movq	%[__o6],%%r8 	\n\t	movaps %%xmm5,0x10(%%r8 )	\n\t	movaps %%xmm13,(%%r8 )	\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t			addpd	%%xmm14,%%xmm14			\n\t"\
				"movaps	(%%rax),%%xmm0			\n\t			movaps	(%%rbx),%%xmm8 			\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t			addpd	%%xmm14,%%xmm13			\n\t"\
			"movq	%[__o1],%%r11	\n\t"\
			"movq	%[__o2],%%r12	\n\t"\
			"movq	%[__o3],%%r13	\n\t"\
			"movq	%[__o4],%%r14	\n\t"\
			"movq	%[__o5],%%r15	\n\t"\
				"subpd	%%xmm2,%%xmm1			\n\t			subpd	%%xmm10,%%xmm9 			\n\t"\
				"subpd	%%xmm2,%%xmm5			\n\t			subpd	%%xmm10,%%xmm13			\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t			movaps	%%xmm9 ,%%xmm14			\n\t"\
				"subpd	%%xmm2,%%xmm3			\n\t			subpd	%%xmm10,%%xmm11			\n\t"\
				"movaps	%%xmm1,%%xmm7			\n\t			movaps	%%xmm9 ,%%xmm15			\n\t"\
				"mulpd	     (%%rsi),%%xmm1		\n\t			mulpd	     (%%rsi),%%xmm9 	\n\t"\
				"movaps	%%xmm1,(%%r11)			\n\t			movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
				"subpd	%%xmm2,%%xmm4			\n\t			subpd	%%xmm10,%%xmm12			\n\t"\
				"mulpd	-0x10(%%rsi),%%xmm2		\n\t			mulpd	-0x10(%%rsi),%%xmm10	\n\t"\
				"addpd	%%xmm5,%%xmm6			\n\t			addpd	%%xmm13,%%xmm14			\n\t"\
				"subpd	%%xmm3,%%xmm7			\n\t			subpd	%%xmm11,%%xmm15			\n\t"\
				"mulpd	 0x60(%%rsi),%%xmm7		\n\t			mulpd	 0x60(%%rsi),%%xmm15	\n\t"\
				"movaps	%%xmm7,(%%r13)			\n\t			movaps	%%xmm15,0x10(%%r13)		\n\t"\
				"movaps	%%xmm3,%%xmm1			\n\t			movaps	%%xmm11,%%xmm9 			\n\t"\
				"addpd	%%xmm4,%%xmm3			\n\t			addpd	%%xmm12,%%xmm11			\n\t"\
				"mulpd	 0x30(%%rsi),%%xmm1		\n\t			mulpd	 0x30(%%rsi),%%xmm9 	\n\t"\
				"movaps	%%xmm1,(%%r12)			\n\t			movaps	%%xmm9 ,0x10(%%r12)		\n\t"\
				"movaps	%%xmm5,%%xmm7			\n\t			movaps	%%xmm13,%%xmm15			\n\t"\
				"subpd	%%xmm4,%%xmm5			\n\t			subpd	%%xmm12,%%xmm13			\n\t"\
				"mulpd	 0x40(%%rsi),%%xmm4		\n\t			mulpd	 0x40(%%rsi),%%xmm12	\n\t"\
				"mulpd	 0x70(%%rsi),%%xmm5		\n\t			mulpd	 0x70(%%rsi),%%xmm13	\n\t"\
				"mulpd	 0x10(%%rsi),%%xmm7		\n\t			mulpd	 0x10(%%rsi),%%xmm15	\n\t"\
				"addpd	%%xmm3,%%xmm2			\n\t			addpd	%%xmm11,%%xmm10			\n\t"\
				"movaps	%%xmm6,%%xmm1			\n\t			movaps	%%xmm14,%%xmm9 			\n\t"\
				"subpd	%%xmm3,%%xmm6			\n\t			subpd	%%xmm11,%%xmm14			\n\t"\
				"mulpd	 0x80(%%rsi),%%xmm6		\n\t			mulpd	 0x80(%%rsi),%%xmm14	\n\t"\
				"mulpd	 0x50(%%rsi),%%xmm3		\n\t			mulpd	 0x50(%%rsi),%%xmm11	\n\t"\
				"addpd	%%xmm1,%%xmm2			\n\t			addpd	%%xmm9 ,%%xmm10			\n\t"\
				"mulpd	 0x20(%%rsi),%%xmm1		\n\t			mulpd	 0x20(%%rsi),%%xmm9 	\n\t"\
				"addpd	%%xmm2,%%xmm0			\n\t			addpd	%%xmm10,%%xmm8 			\n\t"\
				"mulpd	 0x90(%%rsi),%%xmm2		\n\t			mulpd	 0x90(%%rsi),%%xmm10	\n\t"\
				"movaps	%%xmm0,(%%rcx)			\n\t			movaps	%%xmm8 ,(%%rdx)			\n\t"\
				"addpd	%%xmm0,%%xmm2			\n\t			addpd	%%xmm8 ,%%xmm10			\n\t"\
				"addpd	%%xmm1,%%xmm7			\n\t			addpd	%%xmm9 ,%%xmm15			\n\t"\
				"addpd	(%%r11),%%xmm1			\n\t			addpd	0x10(%%r11),%%xmm9 		\n\t"\
				"addpd	%%xmm3,%%xmm4			\n\t			addpd	%%xmm11,%%xmm12			\n\t"\
				"addpd	(%%r12),%%xmm3			\n\t			addpd	0x10(%%r12),%%xmm11		\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t			addpd	%%xmm14,%%xmm13			\n\t"\
				"addpd	(%%r13),%%xmm6			\n\t			addpd	0x10(%%r13),%%xmm14		\n\t"\
				"movaps	%%xmm2,%%xmm0			\n\t			movaps	%%xmm10,%%xmm8 			\n\t"\
				"subpd	%%xmm1,%%xmm2			\n\t			subpd	%%xmm9 ,%%xmm10			\n\t"\
				"addpd	%%xmm0,%%xmm1			\n\t			addpd	%%xmm8 ,%%xmm9 			\n\t"\
				"subpd	%%xmm7,%%xmm2			\n\t			subpd	%%xmm15,%%xmm10			\n\t"\
				"addpd	%%xmm0,%%xmm7			\n\t			addpd	%%xmm8 ,%%xmm15			\n\t"\
				"subpd	%%xmm6,%%xmm1			\n\t			subpd	%%xmm14,%%xmm9 			\n\t"\
				"subpd	%%xmm3,%%xmm2			\n\t			subpd	%%xmm11,%%xmm10			\n\t"\
				"movaps	%%xmm1,(%%r11)			\n\t			movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
				"addpd	%%xmm0,%%xmm3			\n\t			addpd	%%xmm8 ,%%xmm11			\n\t"\
				"subpd	%%xmm5,%%xmm7			\n\t			subpd	%%xmm13,%%xmm15			\n\t"\
				"subpd	%%xmm4,%%xmm2			\n\t			subpd	%%xmm12,%%xmm10			\n\t"\
				"movaps	%%xmm7,(%%r15)			\n\t			movaps	%%xmm15,0x10(%%r15)		\n\t"\
				"movaps	%%xmm2,(%%r12)			\n\t			movaps	%%xmm10,0x10(%%r12)		\n\t"\
				"addpd	%%xmm0,%%xmm4			\n\t			addpd	%%xmm8 ,%%xmm12			\n\t"\
				"addpd	%%xmm6,%%xmm3			\n\t			addpd	%%xmm14,%%xmm11			\n\t"\
				"addpd	%%xmm5,%%xmm4			\n\t			addpd	%%xmm13,%%xmm12			\n\t"\
				"movaps	%%xmm3,(%%r13)			\n\t			movaps	%%xmm11,0x10(%%r13)		\n\t"\
				"movaps	%%xmm4,(%%r14)			\n\t			movaps	%%xmm12,0x10(%%r14)		\n\t"\
			"/***********************************************************************************************************/\n\t"\
			"/* Here are 5 sine terms: Similar to cosines, but t21,19,17,15,13 replace t3,5,7,9,11 and some sign flips: */\n\t"\
			"/***********************************************************************************************************/\n\t"\
				"addq	$0xa0,%%rsi				\n\t"\
				"/*****************************/\n\t			/******************************/\n\t"\
				"/*        Real Parts:         /\n\t			/*     Imaginary Parts:        /\n\t"\
				"/*****************************/\n\t			/******************************/\n\t"\
			"movq	%[__oA],%%r11	\n\t	movaps (%%r11),%%xmm1	\n\t	movaps 0x10(%%r11),%%xmm9 	\n\t"\
			"movq	%[__o9],%%r12	\n\t	movaps (%%r12),%%xmm2	\n\t	movaps 0x10(%%r12),%%xmm10	\n\t"\
			"movq	%[__o8],%%r13	\n\t	movaps (%%r13),%%xmm3	\n\t	movaps 0x10(%%r13),%%xmm11	\n\t"\
			"movq	%[__o7],%%r14	\n\t	movaps (%%r14),%%xmm4	\n\t	movaps 0x10(%%r14),%%xmm12	\n\t"\
			"movq	%[__o6],%%r15	\n\t	movaps (%%r15),%%xmm5	\n\t	movaps 0x10(%%r15),%%xmm13	\n\t"\
				"addpd	%%xmm2,%%xmm1			\n\t			addpd	%%xmm10,%%xmm9 			\n\t"\
				"addpd	%%xmm2,%%xmm5			\n\t			addpd	%%xmm10,%%xmm13			\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t			movaps	%%xmm9 ,%%xmm14			\n\t"\
				"addpd	%%xmm2,%%xmm3			\n\t			addpd	%%xmm10,%%xmm11			\n\t"\
				"movaps	%%xmm1,%%xmm7			\n\t			movaps	%%xmm9 ,%%xmm15			\n\t"\
				"mulpd	     (%%rsi),%%xmm1		\n\t			mulpd	     (%%rsi),%%xmm9 	\n\t"\
				"movaps	%%xmm1,(%%r15)			\n\t			movaps	%%xmm9 ,0x10(%%r15)		\n\t"\
				"addpd	%%xmm2,%%xmm4			\n\t			addpd	%%xmm10,%%xmm12			\n\t"\
				"mulpd	-0xb0(%%rsi),%%xmm2		\n\t			mulpd	-0xb0(%%rsi),%%xmm10	\n\t"\
				"addpd	%%xmm5,%%xmm6			\n\t			addpd	%%xmm13,%%xmm14			\n\t"\
				"subpd	%%xmm3,%%xmm7			\n\t			subpd	%%xmm11,%%xmm15			\n\t"\
				"mulpd	 0x60(%%rsi),%%xmm7		\n\t			mulpd	 0x60(%%rsi),%%xmm15	\n\t"\
				"movaps	%%xmm7,(%%r13)			\n\t			movaps	%%xmm15,0x10(%%r13)		\n\t"\
				"movaps	%%xmm3,%%xmm1			\n\t			movaps	%%xmm11,%%xmm9 			\n\t"\
				"addpd	%%xmm4,%%xmm3			\n\t			addpd	%%xmm12,%%xmm11			\n\t"\
				"mulpd	 0x30(%%rsi),%%xmm1		\n\t			mulpd	 0x30(%%rsi),%%xmm9 	\n\t"\
				"movaps	%%xmm1,(%%r14)			\n\t			movaps	%%xmm9 ,0x10(%%r14)		\n\t"\
				"movaps	%%xmm5,%%xmm7			\n\t			movaps	%%xmm13,%%xmm15			\n\t"\
				"subpd	%%xmm4,%%xmm5			\n\t			subpd	%%xmm12,%%xmm13			\n\t"\
				"mulpd	 0x40(%%rsi),%%xmm4		\n\t			mulpd	 0x40(%%rsi),%%xmm12	\n\t"\
				"mulpd	 0x70(%%rsi),%%xmm5		\n\t			mulpd	 0x70(%%rsi),%%xmm13	\n\t"\
				"mulpd	 0x10(%%rsi),%%xmm7		\n\t			mulpd	 0x10(%%rsi),%%xmm15	\n\t"\
				"subpd	%%xmm3,%%xmm2			\n\t			subpd	%%xmm11,%%xmm10			\n\t"\
				"movaps	%%xmm6,%%xmm1			\n\t			movaps	%%xmm14,%%xmm9 			\n\t"\
				"subpd	%%xmm3,%%xmm6			\n\t			subpd	%%xmm11,%%xmm14			\n\t"\
				"mulpd	 0x80(%%rsi),%%xmm6		\n\t			mulpd	 0x80(%%rsi),%%xmm14	\n\t"\
				"mulpd	 0x50(%%rsi),%%xmm3		\n\t			mulpd	 0x50(%%rsi),%%xmm11	\n\t"\
				"subpd	%%xmm1,%%xmm2			\n\t			subpd	%%xmm9 ,%%xmm10			\n\t"\
				"mulpd	 0x20(%%rsi),%%xmm1		\n\t			mulpd	 0x20(%%rsi),%%xmm9 	\n\t"\
				"mulpd	 0x90(%%rsi),%%xmm2		\n\t			mulpd	 0x90(%%rsi),%%xmm10	\n\t"\
				"addpd	%%xmm1,%%xmm7			\n\t			addpd	%%xmm9 ,%%xmm15			\n\t"\
				"addpd	(%%r15),%%xmm1			\n\t			addpd	0x10(%%r15),%%xmm9 		\n\t"\
				"addpd	%%xmm3,%%xmm4			\n\t			addpd	%%xmm11,%%xmm12			\n\t"\
				"addpd	(%%r14),%%xmm3			\n\t			addpd	0x10(%%r14),%%xmm11		\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t			addpd	%%xmm14,%%xmm13			\n\t"\
				"addpd	(%%r13),%%xmm6			\n\t			addpd	0x10(%%r13),%%xmm14		\n\t"\
				"xorpd	%%xmm0,%%xmm0			\n\t			xorpd	%%xmm8 ,%%xmm8 			\n\t"\
				"subpd	%%xmm2,%%xmm0			\n\t			subpd	%%xmm10,%%xmm8 			\n\t"\
				"addpd	%%xmm1,%%xmm0			\n\t			addpd	%%xmm9 ,%%xmm8 			\n\t"\
				"addpd	%%xmm2,%%xmm1			\n\t			addpd	%%xmm10,%%xmm9 			\n\t"\
				"addpd	%%xmm7,%%xmm0			\n\t			addpd	%%xmm15,%%xmm8 			\n\t"\
				"addpd	%%xmm2,%%xmm7			\n\t			addpd	%%xmm10,%%xmm15			\n\t"\
				"subpd	%%xmm6,%%xmm1			\n\t			subpd	%%xmm14,%%xmm9 			\n\t"\
				"addpd	%%xmm3,%%xmm0			\n\t			addpd	%%xmm11,%%xmm8 			\n\t"\
				"addpd	%%xmm2,%%xmm3			\n\t			addpd	%%xmm10,%%xmm11			\n\t"\
				"subpd	%%xmm5,%%xmm7			\n\t			subpd	%%xmm13,%%xmm15			\n\t"\
				"addpd	%%xmm4,%%xmm0			\n\t			addpd	%%xmm12,%%xmm8 			\n\t"\
				"addpd	%%xmm2,%%xmm4			\n\t			addpd	%%xmm10,%%xmm12			\n\t"\
				"addpd	%%xmm6,%%xmm3			\n\t			addpd	%%xmm14,%%xmm11			\n\t"\
				"addpd	%%xmm5,%%xmm4			\n\t			addpd	%%xmm13,%%xmm12			\n\t"\
				"movaps	%%xmm1,%%xmm2			\n\t			movaps	%%xmm9 ,%%xmm10			\n\t"\
			"movq	%[__o1],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm1			\n\t			addpd	0x10(%%r8 ),%%xmm9 		\n\t"\
				"addpd	%%xmm2,%%xmm2			\n\t			addpd	%%xmm10,%%xmm10			\n\t"\
				"movaps	%%xmm1,(%%r11)			\n\t			movaps	%%xmm9 ,0x10(%%r8 )		\n\t"\
				"subpd	%%xmm2,%%xmm1			\n\t			subpd	%%xmm10,%%xmm9 			\n\t"\
				"movaps	%%xmm1,(%%r8 )			\n\t			movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
				"movaps	%%xmm0,%%xmm5			\n\t			movaps	%%xmm8 ,%%xmm13			\n\t"\
			"movq	%[__o2],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm0			\n\t			addpd	0x10(%%r8 ),%%xmm8 		\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t			addpd	%%xmm13,%%xmm13			\n\t"\
				"movaps	%%xmm0,(%%r12)			\n\t			movaps	%%xmm8 ,0x10(%%r8 )		\n\t"\
				"subpd	%%xmm5,%%xmm0			\n\t			subpd	%%xmm13,%%xmm8 			\n\t"\
				"movaps	%%xmm0,(%%r8 )			\n\t			movaps	%%xmm8 ,0x10(%%r12)		\n\t"\
				"movaps	%%xmm3,%%xmm6			\n\t			movaps	%%xmm11,%%xmm14			\n\t"\
			"movq	%[__o3],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm3			\n\t			addpd	0x10(%%r8 ),%%xmm11		\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t			addpd	%%xmm14,%%xmm14			\n\t"\
				"movaps	%%xmm3,(%%r13)			\n\t			movaps	%%xmm11,0x10(%%r8 )		\n\t"\
				"subpd	%%xmm6,%%xmm3			\n\t			subpd	%%xmm14,%%xmm11			\n\t"\
				"movaps	%%xmm3,(%%r8 )			\n\t			movaps	%%xmm11,0x10(%%r13)		\n\t"\
				"movaps	%%xmm4,%%xmm2			\n\t			movaps	%%xmm12,%%xmm10			\n\t"\
			"movq	%[__o4],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm4			\n\t			addpd	0x10(%%r8 ),%%xmm12		\n\t"\
				"addpd	%%xmm2,%%xmm2			\n\t			addpd	%%xmm10,%%xmm10			\n\t"\
				"movaps	%%xmm4,(%%r14)			\n\t			movaps	%%xmm12,0x10(%%r8 )		\n\t"\
				"subpd	%%xmm2,%%xmm4			\n\t			subpd	%%xmm10,%%xmm12			\n\t"\
				"movaps	%%xmm4,(%%r8 )			\n\t			movaps	%%xmm12,0x10(%%r14)		\n\t"\
				"movaps	%%xmm7,%%xmm5			\n\t			movaps	%%xmm15,%%xmm13			\n\t"\
			"movq	%[__o5],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm7			\n\t			addpd	0x10(%%r8 ),%%xmm15		\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t			addpd	%%xmm13,%%xmm13			\n\t"\
				"movaps	%%xmm7,(%%r15)			\n\t			movaps	%%xmm15,0x10(%%r8 )		\n\t"\
				"subpd	%%xmm5,%%xmm7			\n\t			subpd	%%xmm13,%%xmm15			\n\t"\
				"movaps	%%xmm7,(%%r8 )			\n\t			movaps	%%xmm15,0x10(%%r15)		\n\t"\
				:					/* outputs: none */\
				: [__I0] "m" (XI0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__i9] "m" (Xi9)\
				 ,[__iA] "m" (XiA)\
				 ,[__cc] "m" (Xcc)\
				 ,[__O0] "m" (XO0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				 ,[__o9] "m" (Xo9)\
				 ,[__oA] "m" (XoA)\
				: "cc","memory","rax","rbx","rcx","rdx","rsi","r8","r9","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
				);\
			}

		  #else // USE_64BIT_ASM_STYLE = false:

			// Simple 64-bit-ified version of the 32-bit ASM macro, using just xmm0-7
			#define SSE2_RADIX_11_DFT(XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA, Xcc, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA)\
			{\
			__asm__ volatile (\
				"/********************************************/\n\t"\
				"/*       Here are the 5 cosine terms:       */\n\t"\
				"/********************************************/\n\t"\
				"movq	%[__I0],%%rax			\n\t"\
				"movq	%[__cc],%%rbx			\n\t"\
				"movq	%[__O0],%%rcx			\n\t"\
			"movq	%[__i1],%%r8 	\n\t	movaps (%%r8 ),%%xmm1	\n\t"\
			"movq	%[__iA],%%r9 	\n\t	movaps (%%r9 ),%%xmm5	\n\t"\
			"movq	%[__i2],%%r10	\n\t	movaps (%%r10),%%xmm2	\n\t"\
			"movq	%[__i9],%%r11	\n\t	movaps (%%r11),%%xmm6	\n\t"\
			"movq	%[__i3],%%r12	\n\t	movaps (%%r12),%%xmm3	\n\t"\
			"movq	%[__i8],%%r13	\n\t	movaps (%%r13),%%xmm7	\n\t"\
			"movq	%[__i4],%%r14	\n\t	movaps (%%r14),%%xmm4	\n\t"\
			"movq	%[__i7],%%r15	\n\t	movaps (%%r15),%%xmm0	\n\t"\
				/*addq	$0x10,%%rcx	<*** This is so following temp-outs go into Im-slots of o-array; now add 0x10 to __o* addresses to make that happen */\
				"subpd	%%xmm5,%%xmm1			\n\t"\
				"subpd	%%xmm6,%%xmm2			\n\t"\
				"subpd	%%xmm7,%%xmm3			\n\t"\
				"subpd	%%xmm0,%%xmm4			\n\t"\
			"movq	%[__oA],%%r8 	\n\t	movaps %%xmm1,0x10(%%r8 )	\n\t"\
			"movq	%[__o9],%%r9 	\n\t	movaps %%xmm2,0x10(%%r9 )	\n\t"\
			"movq	%[__o8],%%r10	\n\t	movaps %%xmm3,0x10(%%r10)	\n\t"\
			"movq	%[__o7],%%r11	\n\t	movaps %%xmm4,0x10(%%r11)	\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t"\
				"addpd	%%xmm7,%%xmm7			\n\t"\
				"addpd	%%xmm0,%%xmm0			\n\t"\
				"addpd	%%xmm5,%%xmm1			\n\t"\
				"addpd	%%xmm6,%%xmm2			\n\t"\
			"movq	%[__i5],%%r9 	\n\t	movaps (%%r9 ),%%xmm5	\n\t"\
			"movq	%[__i6],%%r11	\n\t	movaps (%%r11),%%xmm6	\n\t"\
				"addpd	%%xmm7,%%xmm3			\n\t"\
				"addpd	%%xmm0,%%xmm4			\n\t"\
				"subpd	%%xmm6,%%xmm5			\n\t"\
			"movq	%[__o6],%%r8 	\n\t	movaps %%xmm5,0x10(%%r8 )	\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t"\
				"movaps	(%%rax),%%xmm0			\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t"\
				"/********************************************/\n\t"\
				"/*               Real Parts:                */\n\t"\
				"/********************************************/\n\t"\
			"movq	%[__o1],%%r11	\n\t"\
			"movq	%[__o2],%%r12	\n\t"\
			"movq	%[__o3],%%r13	\n\t"\
			"movq	%[__o4],%%r14	\n\t"\
			"movq	%[__o5],%%r15	\n\t"\
				/*subq	$0x10,%%rcx <*** Back to working with Re o-array slots */\
				"subpd	%%xmm2,%%xmm1			\n\t"\
				"subpd	%%xmm2,%%xmm5			\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t"\
				"subpd	%%xmm2,%%xmm3			\n\t"\
				"movaps	%%xmm1,%%xmm7			\n\t"\
				"mulpd	     (%%rbx),%%xmm1		\n\t"\
				"movaps	%%xmm1,(%%r11)	\n\t"\
				"subpd	%%xmm2,%%xmm4			\n\t"\
				"mulpd	-0x10(%%rbx),%%xmm2		\n\t"\
				"addpd	%%xmm5,%%xmm6			\n\t"\
				"subpd	%%xmm3,%%xmm7			\n\t"\
				"mulpd	 0x60(%%rbx),%%xmm7		\n\t"\
				"movaps	%%xmm7,(%%r13)	\n\t"\
				"movaps	%%xmm3,%%xmm1			\n\t"\
				"addpd	%%xmm4,%%xmm3			\n\t"\
				"mulpd	 0x30(%%rbx),%%xmm1		\n\t"\
				"movaps	%%xmm1,(%%r12)	\n\t"\
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
				"addpd	(%%r11),%%xmm1	\n\t"\
				"addpd	%%xmm3,%%xmm4			\n\t"\
				"addpd	(%%r12),%%xmm3	\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t"\
				"addpd	(%%r13),%%xmm6	\n\t"\
				"movaps	%%xmm2,%%xmm0			\n\t"\
				"subpd	%%xmm1,%%xmm2			\n\t"\
				"addpd	%%xmm0,%%xmm1			\n\t"\
				"subpd	%%xmm7,%%xmm2			\n\t"\
				"addpd	%%xmm0,%%xmm7			\n\t"\
				"subpd	%%xmm6,%%xmm1			\n\t"\
				"subpd	%%xmm3,%%xmm2			\n\t"\
				"movaps	%%xmm1,(%%r11)	\n\t"\
				"addpd	%%xmm0,%%xmm3			\n\t"\
				"subpd	%%xmm5,%%xmm7			\n\t"\
				"subpd	%%xmm4,%%xmm2			\n\t"\
				"movaps	%%xmm7,(%%r15)	\n\t"\
				"movaps	%%xmm2,(%%r12)	\n\t"\
				"addpd	%%xmm0,%%xmm4			\n\t"\
				"addpd	%%xmm6,%%xmm3			\n\t"\
				"addpd	%%xmm5,%%xmm4			\n\t"\
				"movaps	%%xmm3,(%%r13)	\n\t"\
				"movaps	%%xmm4,(%%r14)	\n\t"\
				"/********************************************/\n\t"\
				"/*          Imaginary Parts:                */\n\t"\
				"/********************************************/\n\t"\
				/*addq	$0x10,%%rax	<*** ensuing block works with Im i-data */\
			"movq	%[__i1],%%r8 	\n\t	movaps 0x10(%%r8 ),%%xmm1	\n\t"\
			"movq	%[__iA],%%r9 	\n\t	movaps 0x10(%%r9 ),%%xmm5	\n\t"\
			"movq	%[__i2],%%r10	\n\t	movaps 0x10(%%r10),%%xmm2	\n\t"\
			"movq	%[__i9],%%r11	\n\t	movaps 0x10(%%r11),%%xmm6	\n\t"\
			"movq	%[__i3],%%r12	\n\t	movaps 0x10(%%r12),%%xmm3	\n\t"\
			"movq	%[__i8],%%r13	\n\t	movaps 0x10(%%r13),%%xmm7	\n\t"\
			"movq	%[__i4],%%r14	\n\t	movaps 0x10(%%r14),%%xmm4	\n\t"\
			"movq	%[__i7],%%r15	\n\t	movaps 0x10(%%r15),%%xmm0	\n\t"\
				"subpd	%%xmm5,%%xmm1			\n\t"\
				"subpd	%%xmm6,%%xmm2			\n\t"\
				"subpd	%%xmm7,%%xmm3			\n\t"\
				"subpd	%%xmm0,%%xmm4			\n\t"\
				/* Now this spill quartet is into Re slots of o-array: */\
			"movq	%[__oA],%%r8 	\n\t	movaps %%xmm1,(%%r8 )	\n\t"\
			"movq	%[__o9],%%r9 	\n\t	movaps %%xmm2,(%%r9 )	\n\t"\
			"movq	%[__o8],%%r10	\n\t	movaps %%xmm3,(%%r10)	\n\t"\
			"movq	%[__o7],%%r11	\n\t	movaps %%xmm4,(%%r11)	\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t"\
				"addpd	%%xmm7,%%xmm7			\n\t"\
				"addpd	%%xmm0,%%xmm0			\n\t"\
				"addpd	%%xmm5,%%xmm1			\n\t"\
				"addpd	%%xmm6,%%xmm2			\n\t"\
			"movq	%[__i5],%%r9 	\n\t	movaps 0x10(%%r9 ),%%xmm5	\n\t"\
			"movq	%[__i6],%%r11	\n\t	movaps 0x10(%%r11),%%xmm6	\n\t"\
				"addpd	%%xmm7,%%xmm3			\n\t"\
				"addpd	%%xmm0,%%xmm4			\n\t"\
				"subpd	%%xmm6,%%xmm5			\n\t"\
			"movq	%[__o6],%%r8 	\n\t	movaps %%xmm5,(%%r8 )	\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t"\
				"movaps	0x10(%%rax),%%xmm0			\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t"\
				/*addq	$0x10,%%rcx	 <*** switch to Im slots of o-array: */\
			"movq	%[__o1],%%r11	\n\t"\
			"movq	%[__o2],%%r12	\n\t"\
			"movq	%[__o3],%%r13	\n\t"\
			"movq	%[__o4],%%r14	\n\t"\
			"movq	%[__o5],%%r15	\n\t"\
				"subpd	%%xmm2,%%xmm1			\n\t"\
				"subpd	%%xmm2,%%xmm5			\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t"\
				"subpd	%%xmm2,%%xmm3			\n\t"\
				"movaps	%%xmm1,%%xmm7			\n\t"\
				"mulpd	     (%%rbx),%%xmm1		\n\t"\
				"movaps	%%xmm1,0x10(%%r11)	\n\t"\
				"subpd	%%xmm2,%%xmm4			\n\t"\
				"mulpd	-0x10(%%rbx),%%xmm2		\n\t"\
				"addpd	%%xmm5,%%xmm6			\n\t"\
				"subpd	%%xmm3,%%xmm7			\n\t"\
				"mulpd	 0x60(%%rbx),%%xmm7		\n\t"\
				"movaps	%%xmm7,0x10(%%r13)	\n\t"\
				"movaps	%%xmm3,%%xmm1			\n\t"\
				"addpd	%%xmm4,%%xmm3			\n\t"\
				"mulpd	 0x30(%%rbx),%%xmm1		\n\t"\
				"movaps	%%xmm1,0x10(%%r12)	\n\t"\
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
				"movaps	%%xmm0,0x10(%%rcx)			\n\t"\
				"addpd	%%xmm0,%%xmm2			\n\t"\
				"addpd	%%xmm1,%%xmm7			\n\t"\
				"addpd	0x10(%%r11),%%xmm1	\n\t"\
				"addpd	%%xmm3,%%xmm4			\n\t"\
				"addpd	0x10(%%r12),%%xmm3	\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t"\
				"addpd	0x10(%%r13),%%xmm6	\n\t"\
				"movaps	%%xmm2,%%xmm0			\n\t"\
				"subpd	%%xmm1,%%xmm2			\n\t"\
				"addpd	%%xmm0,%%xmm1			\n\t"\
				"subpd	%%xmm7,%%xmm2			\n\t"\
				"addpd	%%xmm0,%%xmm7			\n\t"\
				"subpd	%%xmm6,%%xmm1			\n\t"\
				"subpd	%%xmm3,%%xmm2			\n\t"\
				"movaps	%%xmm1,0x10(%%r11)	\n\t"\
				"addpd	%%xmm0,%%xmm3			\n\t"\
				"subpd	%%xmm5,%%xmm7			\n\t"\
				"subpd	%%xmm4,%%xmm2			\n\t"\
				"movaps	%%xmm7,0x10(%%r15)	\n\t"\
				"movaps	%%xmm2,0x10(%%r12)	\n\t"\
				"addpd	%%xmm0,%%xmm4			\n\t"\
				"addpd	%%xmm6,%%xmm3			\n\t"\
				"addpd	%%xmm5,%%xmm4			\n\t"\
				"movaps	%%xmm3,0x10(%%r13)	\n\t"\
				"movaps	%%xmm4,0x10(%%r14)	\n\t"\
				"/************************************************************************************************************************************/\n\t"\
				"/* Here are the 5 sine terms: Similar sequence to cosine terms, but with t21,19,17,15,13 replacing t3,5,7,9,11 and some sign flips: */\n\t"\
				"/************************************************************************************************************************************/\n\t"\
				/*subq	$0x10,%%rax	<*** Re parts of i-array; on o-side still working with Im-data */\
				"addq	$0xa0,%%rbx				\n\t"\
				"/********************************************/\n\t"\
				"/*               Real Parts:                */\n\t"\
				"/********************************************/\n\t"\
			"movq	%[__oA],%%r11	\n\t	movaps 0x10(%%r11),%%xmm1	\n\t"\
			"movq	%[__o9],%%r12	\n\t	movaps 0x10(%%r12),%%xmm2	\n\t"\
			"movq	%[__o8],%%r13	\n\t	movaps 0x10(%%r13),%%xmm3	\n\t"\
			"movq	%[__o7],%%r14	\n\t	movaps 0x10(%%r14),%%xmm4	\n\t"\
			"movq	%[__o6],%%r15	\n\t	movaps 0x10(%%r15),%%xmm5	\n\t"\
				/* r8-10 still free for o-addressing */\
				"addpd	%%xmm2,%%xmm1			\n\t"\
				"addpd	%%xmm2,%%xmm5			\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t"\
				"addpd	%%xmm2,%%xmm3			\n\t"\
				"movaps	%%xmm1,%%xmm7			\n\t"\
				"mulpd	     (%%rbx),%%xmm1		\n\t"\
				"movaps	%%xmm1,0x10(%%r15)	\n\t"\
				"addpd	%%xmm2,%%xmm4			\n\t"\
				"mulpd	-0xb0(%%rbx),%%xmm2		\n\t"\
				"addpd	%%xmm5,%%xmm6			\n\t"\
				"subpd	%%xmm3,%%xmm7			\n\t"\
				"mulpd	 0x60(%%rbx),%%xmm7		\n\t"\
				"movaps	%%xmm7,0x10(%%r13)	\n\t"\
				"movaps	%%xmm3,%%xmm1			\n\t"\
				"addpd	%%xmm4,%%xmm3			\n\t"\
				"mulpd	 0x30(%%rbx),%%xmm1		\n\t"\
				"movaps	%%xmm1,0x10(%%r14)	\n\t"\
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
				"addpd	0x10(%%r15),%%xmm1	\n\t"\
				"addpd	%%xmm3,%%xmm4			\n\t"\
				"addpd	0x10(%%r14),%%xmm3	\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t"\
				"addpd	0x10(%%r13),%%xmm6	\n\t"\
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
			"movq	%[__o1],%%r8 	\n\t"\
				"addpd	0x10(%%r8 ),%%xmm1	\n\t"\
				"addpd	%%xmm2,%%xmm2			\n\t"\
				"movaps	%%xmm1,0x10(%%r8 )	\n\t"\
				"subpd	%%xmm2,%%xmm1			\n\t"\
				"movaps	%%xmm1,0x10(%%r11)	\n\t"\
				"movaps	%%xmm0,%%xmm5			\n\t"\
			"movq	%[__o2],%%r8 	\n\t"\
				"addpd	0x10(%%r8 ),%%xmm0	\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t"\
				"movaps	%%xmm0,0x10(%%r8 )	\n\t"\
				"subpd	%%xmm5,%%xmm0			\n\t"\
				"movaps	%%xmm0,0x10(%%r12)	\n\t"\
				"movaps	%%xmm3,%%xmm6			\n\t"\
			"movq	%[__o3],%%r8 	\n\t"\
				"addpd	0x10(%%r8 ),%%xmm3	\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t"\
				"movaps	%%xmm3,0x10(%%r8 )	\n\t"\
				"subpd	%%xmm6,%%xmm3			\n\t"\
				"movaps	%%xmm3,0x10(%%r13)	\n\t"\
				"movaps	%%xmm4,%%xmm2			\n\t"\
			"movq	%[__o4],%%r8 	\n\t"\
				"addpd	0x10(%%r8 ),%%xmm4	\n\t"\
				"addpd	%%xmm2,%%xmm2			\n\t"\
				"movaps	%%xmm4,0x10(%%r8 )	\n\t"\
				"subpd	%%xmm2,%%xmm4			\n\t"\
				"movaps	%%xmm4,0x10(%%r14)	\n\t"\
				"movaps	%%xmm7,%%xmm5			\n\t"\
			"movq	%[__o5],%%r8 	\n\t"\
				"addpd	0x10(%%r8 ),%%xmm7	\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t"\
				"movaps	%%xmm7,0x10(%%r8 )	\n\t"\
				"subpd	%%xmm5,%%xmm7			\n\t"\
				"movaps	%%xmm7,0x10(%%r15)	\n\t"\
				"/********************************************/\n\t"\
				"/*          Imaginary Parts:                */\n\t"\
				"/********************************************/\n\t"\
				/* subq	$0x10,%%rcx	<*** Temps needed for Im outputs stored in Re slots of o-array (I know it seems counterintuitive) */\
			"movq	%[__oA],%%r11	\n\t	movaps (%%r11),%%xmm1	\n\t"\
			"movq	%[__o9],%%r12	\n\t	movaps (%%r12),%%xmm2	\n\t"\
			"movq	%[__o8],%%r13	\n\t	movaps (%%r13),%%xmm3	\n\t"\
			"movq	%[__o7],%%r14	\n\t	movaps (%%r14),%%xmm4	\n\t"\
			"movq	%[__o6],%%r15	\n\t	movaps (%%r15),%%xmm5	\n\t"\
				"addpd	%%xmm2,%%xmm1			\n\t"\
				"addpd	%%xmm2,%%xmm5			\n\t"\
				"movaps	%%xmm1,%%xmm6			\n\t"\
				"addpd	%%xmm2,%%xmm3			\n\t"\
				"movaps	%%xmm1,%%xmm7			\n\t"\
				"mulpd	     (%%rbx),%%xmm1		\n\t"\
				"movaps	%%xmm1,(%%r15)	\n\t"\
				"addpd	%%xmm2,%%xmm4			\n\t"\
				"mulpd	-0xb0(%%rbx),%%xmm2		\n\t"\
				"addpd	%%xmm5,%%xmm6			\n\t"\
				"subpd	%%xmm3,%%xmm7			\n\t"\
				"mulpd	 0x60(%%rbx),%%xmm7		\n\t"\
				"movaps	%%xmm7,(%%r13)	\n\t"\
				"movaps	%%xmm3,%%xmm1			\n\t"\
				"addpd	%%xmm4,%%xmm3			\n\t"\
				"mulpd	 0x30(%%rbx),%%xmm1		\n\t"\
				"movaps	%%xmm1,(%%r14)	\n\t"\
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
				"addpd	(%%r15),%%xmm1	\n\t"\
				"addpd	%%xmm3,%%xmm4			\n\t"\
				"addpd	(%%r14),%%xmm3	\n\t"\
				"addpd	%%xmm6,%%xmm5			\n\t"\
				"addpd	(%%r13),%%xmm6	\n\t"\
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
			"movq	%[__o1],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm1	\n\t"\
				"addpd	%%xmm2,%%xmm2			\n\t"\
				"movaps	%%xmm1,(%%r11)	\n\t"\
				"subpd	%%xmm2,%%xmm1			\n\t"\
				"movaps	%%xmm1,(%%r8 )	\n\t"\
				"movaps	%%xmm0,%%xmm5			\n\t"\
			"movq	%[__o2],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm0	\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t"\
				"movaps	%%xmm0,(%%r12)	\n\t"\
				"subpd	%%xmm5,%%xmm0			\n\t"\
				"movaps	%%xmm0,(%%r8 )	\n\t"\
				"movaps	%%xmm3,%%xmm6			\n\t"\
			"movq	%[__o3],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm3	\n\t"\
				"addpd	%%xmm6,%%xmm6			\n\t"\
				"movaps	%%xmm3,(%%r13)	\n\t"\
				"subpd	%%xmm6,%%xmm3			\n\t"\
				"movaps	%%xmm3,(%%r8 )	\n\t"\
				"movaps	%%xmm4,%%xmm2			\n\t"\
			"movq	%[__o4],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm4	\n\t"\
				"addpd	%%xmm2,%%xmm2			\n\t"\
				"movaps	%%xmm4,(%%r14)	\n\t"\
				"subpd	%%xmm2,%%xmm4			\n\t"\
				"movaps	%%xmm4,(%%r8 )	\n\t"\
				"movaps	%%xmm7,%%xmm5			\n\t"\
			"movq	%[__o5],%%r8 	\n\t"\
				"addpd	(%%r8 ),%%xmm7	\n\t"\
				"addpd	%%xmm5,%%xmm5			\n\t"\
				"movaps	%%xmm7,(%%r15)	\n\t"\
				"subpd	%%xmm5,%%xmm7			\n\t"\
				"movaps	%%xmm7,(%%r8 )	\n\t"\
				:					/* outputs: none */\
				: [__I0] "m" (XI0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__i9] "m" (Xi9)\
				 ,[__iA] "m" (XiA)\
				 ,[__cc] "m" (Xcc)\
				 ,[__O0] "m" (XO0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				 ,[__o9] "m" (Xo9)\
				 ,[__oA] "m" (XoA)\
				: "cc","memory","rax","rbx","rcx","r8","r9","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
				);\
			}

		  #endif // USE_64BIT_ASM_STYLE

		#endif	// AVX and 64-bit SSE2

	  #endif	// IF(GCC), USE 32/64-BIT ASM STYLE

	#endif	// MSVC / GCC

#endif	/* radix11_sse_macro_h_included */

